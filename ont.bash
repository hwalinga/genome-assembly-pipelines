#!/bin/bash

###################
# PARSE ARGUMENTS #
###################

HELP=$(cat << 'EOF'
Usage:
$0 [--options] "FolderPath" OR/AND "FastqFiles"
--test:
    Show commands to run, but not execute them.
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
--nocov:
    Do not plot the coverage plots.
--help,-h
    Plot this help and exit.
EOF
)

echo "Parsing arguments"

NOMOVE=false
KEEP=false
BACTERIA=true
PHAGE=true
TEST=""
COVERAGE=true
PDF=true
ERRORLOG="error.log"

while [[ -n "$1" ]] && [[ "$1" =~ ^- ]]; do
    case "$1" in
        --test)
            TEST="--dry-run"
            ;;
        --keep)
            KEEP=true
            ;;
        --nomove)
            NOMOVE=true
            ;;
        --nocov)
            COVERAGE=false
            ;;
        -h|--help)
            >&2 echo "Printing help:"
            echo $HELP
            >&2 echo "Exiting"
            exit 1
            ;;
        --)
            break
            ;;
        -*)
            >&2 echo "Unrecognized option"
            >&2 echo "Exiting."
            exit 1
            ;;
    esac
    shift
done

FILES=()
DIRS=()
for i in "$@"; do
    if [[ -f "$i" ]]; then
        FILES+=("$i")
    elif [[ -d "$i" ]]; then
        DIRS+=("$i")
    else
        >&2 echo "Something wrong with $i, skipping."
    fi
done;

echo "Arguments parsed:"
echo "files: ${FILES[@}]}"
echo "dirs: ${DIRS[@]}"
if [[ $COVERAGE == "true" ]]; then
    echo "Run with making coverage plots."
fi
if [[ $TEST == "--dry-run" ]]; then
    echo "Only output what commands are run (dry run)."
fi
if [[ $KEEP == "true" ]]; then
    echo "Keep inbetween results."
fi

#########################
# CHECKING DEPENDENCIES #
#########################

program_missing () {
    $program=$1
    >&2 echo "$program not installed"
    >&2 echo "You can install with:"
    >&2 echo "conda -c bioconda install ${program%.*} (Anaconda)"
    >&2 echo "or:"
    >&2 echo "brew install ${program%.*} (Homebrew)"
    >&2 echo "Provided that one of these package managers is installed."
}

echo "Checking dependencies."
echo "Checking core dependencies."
for program in free ps printf awk perl sed; do
    if [[ ! `command -v $program` ]]; then
        >&2 echo "$program not installed."
        >&2 echo "This is a core OS program."
        >&2 echo "You should ask you sysadmin to install this."
        >&2 echo "Exiting"
        exit 1
    fi
done
if [[ ! `samtools --version` ]]; then
    >&2 echo "You have installed an too old version for samtools."
    >&2 echo "Please install a version above 1.0."
    >&2 echo "This program is only needed for the coverage plots."
    >&2 echo "If you do not ask for coverage plots with the --nocov option,"
    >&2 echo "you can still assemble."
    >&2 echo "Exiting"
    exit 1
fi
for program in parallel qcat filtlong flye medaka racon; do
    if [[ ! `command -v $program` ]]; then
        program_missing $program
        >&2 echo "Exiting"
        exit 1

    fi
done
if [[ ! `command -v NanoFilt` ]]; then
    program_missing nanofilt
    >&2 echo "Exiting"
    exit 1
fi
if [[ $COVERAGE == "true" ]]; then
    echo "Checking dependencies for coverage plots."
    echo "If these dependencies fail, you can still assemble with the --nocov option."
    for program in gnuplot; do
        if [[ ! `command -v $program` ]]; then
            program_missing $program
            >&2 echo "Exiting"
            exit 1
        fi
    done
fi

echo "Dependencies met."

#########################################
# GNU Parallel helper functions/strings #
#########################################

# I am not writing/using a full blown load manager,
# but GNU Parallel can do its best.
CALCFREEPROC () {
    # Calculate the amount of free processing powers left (in # of processors)
    ps -eo pcpu | awk -v P=`nproc` 'NR!=1{S+=$1}END{printf "%.2f",P-S/100}'
}
CALCFREEMEM () {
    # Calculate the amount of RAM left (in GB)
    free --giga | awk 'NR==2{print $NF}'
}
# $TEST can be the --dry-run option.
PARALLEL () {
    # PARALLEL (as the/this function) should work equivalent as parallel (as the command)
    # Calculating free processors and free memory at call time.
    # However, I am not a 100 % percent sure, as using parallel often encounters
    # complicated quoting rules, and I don't know if bash correctly solves them like this.
    #
    # If the first argument matches -j*, than the number after -j is taken as
    # the amount of processors.
    # The free mem is the amount of free memory divided by the amount of
    # processors.
    if [[ $1 =~ "-j*" ]]; then
        PROC=${1#-j*}
        shift
    else
        PROC=$(CALCFREEPROC)
    fi
    MEM="$(printf '%.f' $(bc -l <<<"$(CALCFREEMEM) / $PROC"))G"
    PROC=$(printf '%.f' $PROC)

    parallel $TEST -j $PROC --memfree $MEM --load 100% --delay 30 "$@"
}
# GNU Parallel --rpl strings.
# Perl regex to find the first directory:
FIRSTDIRECTORY='{m} s:.*/(?=.*/)::;s:/.*::;'

# move

# Change this here so that user can just serve the whole file instead.

# TODO: Remove inbetween results.

INPUT=ont.fastq
cat "${FILES[@]}" >> $INPUT
find "${DIRS[@]}" -type f | xargs cat >> $INPUT
# demultiplex
mkdir -p demultiplex
>&2 echo "Using kit NBD103/NBD104, change script to change."
qcat -f $INPUT -b demultiplex --trim -k NBD103/NBD104 --detect-middle
# filter good bar codes
mkdir -p demultiplex_filter
find demultiplex ! -name none.fastq -type f -size +1M | xarp cp -t demultiplex_filter

# Filter on quality
mkdir -p nanofilt
PARALLEL NanoFilt -q 7 -l 500 {} ">" nanofilt/{/} ::: demultiplex_filter/*
mkdir -p filtlong
PARALLEL filtlong -t 500000000 {} ">" filtlong/{/} ::: nanofilt/*

# assemble
mkdir -p flye
>&2 echo "Assembling for phage target length (40k bp)."
PARALLEL -j1 flye -g 40k -m 3000 -i 4 --meta -t 8 --nano-raw {} -o flye/{/.} ::: filtlong/*

# polishing
mkdir -p racon{0..4}
mkdir -p mapped{0..3}
parallel --dry-run cp {}/assembly.fasta racon0/{/}.fa ::: flye/* | sh
PARALLEL \
    "for i in {0..3}; do minimap2 -ax map-ont racon\$i/{/.}.fa filtlong/{/.}.fastq > mapped\$i/{/.}.sam && racon -t 4 -m 8 -x -6 -g -8 -w 500 filtlong/{/.}.fastq mapped\$i/{/.}.sam racon\$i/{/.}.fa > racon\$((i+1))/{/.}.fa; done" ::: racon0/*
>&2 echo "Assuming r941_flip235 base call model."
PARALLEL -j2 medaka_consensus -i filtlong/{/.}.fastq -d {} -o medaka/{/.} -t 4 -m r941_flip235 ::: racon4/*


# Coverage plots.
if [[ $COVERAGE == "true" ]]; then
    mkdir -p {medaka_mapped,stats,figs}

    echo "Mapping to raw reads."
    PARALLEL \
        test -s {}.consensus.fasta "&&"
        minimap2 -ax map-ont {} filtlong/{/}.fastq "|" \
        awk -F\\t -v OFS=\\t "'{\$1=substr(\$1,1,251)}1'" \
        ">" medaka_mapped/{/}.sam \
        ::: medaka/*

    echo "Sorting and indexing of sam files."
    PARALLEL --plus 'samtools sort {} > {.}.bam && samtools index {.}.bam && samtools calmd {.}.bam medaka/{/.}/consensus.fasta > {.}.md.sam && samtools sort {.}.md.sam > {.}.md.bam && samtools index {.}.md.bam' \
        ::: medaka_mapped/*.sam

    echo "Calculating depth."
    parallel --dry-run mkdir -p stats/{/.} ::: mapped/*.sam | sh
    shopt -s extglob
    PARALLEL 'samtools depth -a {} | awk "{print $3 > \"stats/{\.}/\"\$1}"' \
        ::: medaka_mapped/!(*.md).bam
    shopt -u extglob

    COVPLOT=$(cat <<-'EOF'
    set term png;
    set output figs.''/'.sample.'_'.contig;
    set title sample.'-'.contig noenhanced;
    set xlabel 'bp';
    set ylabel 'coverage';
    unset key;
    stats file nooutput;
    plot [:STATS_records] file with dots;
EOF
    )

    echo "Plotting coverage depth."
    PARALLEL -q --rpl $FIRSTDIRECTORY \
        gnuplot -e "sample='{m}'; contig='{/}'; file='{}'; figs='figs'"$COVPLOT \
        ::: stats/*/*
fi

# Variance calling / sequence contribution?