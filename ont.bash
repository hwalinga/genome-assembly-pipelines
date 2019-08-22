#!/bin/bash

###################
# PARSE ARGUMENTS #
###################

HELP=$(cat << EOF
Usage:
$0 [--options] "FolderPath" OR/AND "FastqFiles"
--test:
    Show commands to run, but not execute them.
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
    (Currently not implemented.)
--nocov:
    Do not plot the coverage plots.
-i [FILE]
    Instead of supplying as an argument you can pass the fastq file
    with the -i option. Note that this way the file will not be copied.
    (Copying might be desirable if it is on an unstable filesystem.)
--help,-h
    Plot this help and exit.
EOF
)

echo "========"
echo "Parsing arguments"
echo "========"

KEEP=false
TEST=""
COVERAGE=true
ERRORLOG="error.log"
INPUT=""

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
            echo "$HELP"
            >&2 echo "Exiting"
            exit 1
            ;;
        -i)
            shift
            INPUT="$1"
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

echo "========"
echo "Arguments parsed:"
echo "files: ${FILES[@]}"
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
echo "========"

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

echo "========"
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
    if [[ -z "$(gnuplot --version | awk -F' |\\.' '$2>=5')" ]]; then
        >&2 echo "gnuplot not installed or a version installed lower than 5."
        program_missing $pgram
        >&2 echo "You can still assemble with the --nocov option."
        >&2 echo "Exiting"
        exit 1
    fi
fi

echo "========"
echo "Dependencies met."
echo "========"

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

    parallel --plus $TEST -j $PROC --memfree $MEM --load 100% --delay 30 "$@"
}
if [[ $TEST == "--dry-run" ]]; then
    MKDIR="echo mkdir -p"
    CP="echo cp"
    KEEP=true
else
    MKDIR="mkdir -p"
    CP="cp"
fi

# GNU Parallel --rpl strings.
# Perl regex to find the first directory:
FIRSTDIRECTORY='{m} s:.*/(?=.*/)::;s:/.*::;'

# MOVE

# TODO: Remove inbetween results.

echo "$INPUT"
if [[ -z "$INPUT" ]]; then
    INPUT=ont.fastq
    if [[ $TEST == "--dry-run" ]]; then
        echo "FILES"
        echo "${FILES[@]}"
        echo "DIRS"
        echo "${DIRS[@]}"
    else
        echo "We are collecting the fastq files"
        if [[ ! -z "$FILES" ]]; then
            cat "${FILES[@]}" >> $INPUT
        fi
        if [[ ! -z "$DIRS" ]]; then
            find "${DIRS[@]}" -type f -print0 | xargs -0 cat >> $INPUT
        fi
    fi
fi

if [[ ! -s "$INPUT" ]]; then
    rm -f "$INPUT"
    >&2 echo "Could not find any fastq reads in the input arguments you provided."
    >&2 echo "Exiting"
    exit 1
fi

# demultiplex
$MKDIR demultiplex
>&2 echo "Using kit NBD103/NBD104, change script to change."
if [[ $TEST == "--dry-run" ]]; then
    echo "qcat -f $INPUT -b demultiplex --trim -k NBD103/NBD104 --detect-middle"
else
    qcat -f $INPUT -b demultiplex --trim -k NBD103/NBD104 --detect-middle
fi
# filter good bar codes
$MKDIR demultiplex_filter
if [[ $TEST == "--dry-run" ]]; then
    echo "find demultiplex ! -name none.fastq -type f -size +1M | xarp cp -t demultiplex_filter"
else
    find demultiplex ! -name none.fastq -type f -size +1M | xargs cp -t demultiplex_filter
fi

# Filter on quality
$MKDIR nanofilt
PARALLEL NanoFilt -q 7 -l 500 {} ">" nanofilt/{/} ::: demultiplex_filter/*
$MKDIR filtlong
PARALLEL filtlong -t 500000000 {} ">" filtlong/{/} ::: nanofilt/*

# assemble
$MKDIR flye
>&2 echo "Assembling for phage target length (40k bp)."
PARALLEL -j1 flye -g 40k -m 3000 -i 4 --meta -t 8 --nano-raw {} -o flye/{/.} ::: filtlong/*

# polishing
$MKDIR racon{0..4}
$MKDIR mapped{0..3}
parallel --dry-run cp {}/assembly.fasta racon0/{/}.fa ::: flye/* | sh
PARALLEL \
    "for i in {0..3}; do minimap2 -ax map-ont racon\$i/{/.}.fa filtlong/{/.}.fastq > mapped\$i/{/.}.sam && racon -t 4 -m 8 -x -6 -g -8 -w 500 filtlong/{/.}.fastq mapped\$i/{/.}.sam racon\$i/{/.}.fa > racon\$((i+1))/{/.}.fa; done" ::: racon0/*
>&2 echo "Assuming r941_flip235 base call model."
PARALLEL -j2 medaka_consensus -i filtlong/{/.}.fastq -d {} -o medaka/{/.} -t 4 -m r941_flip235 ::: racon4/*


shopt -s extglob
echo "========"
echo "Making coverage plots"
echo "========"
# Coverage plots.
if [[ $COVERAGE == "true" ]]; then
    $MKDIR {medaka_mapped,stats,figs}

    echo "Mapping to raw reads."
    PARALLEL \
        test -s {}/consensus.fasta "&&" \
        minimap2 -ax map-ont {}/consensus.fasta filtlong/{/}.fastq "|" \
        awk -F\\\\t -v OFS=\\\\t "'{\$1=substr(\$1,1,251)}1'" \
        ">" medaka_mapped/{/}.sam \
        ::: medaka/*

    echo "Sorting and indexing of sam files."
    PARALLEL 'samtools sort {} > {.}.bam && samtools index {.}.bam && samtools calmd {.}.bam medaka/{/.}/consensus.fasta > {.}.md.sam && samtools sort {.}.md.sam > {.}.md.bam && samtools index {.}.md.bam' \
        ::: medaka_mapped/!(*.md).sam

    echo "Calculating depth."
    parallel --dry-run $MKDIR stats/{/.} ::: medaka_mapped/!(*.md).sam | sh
    PARALLEL 'samtools depth -a {} | awk "{gsub(\":\",\"_\",\$1); print \$3 > \"stats/{/.}/\"\$1}"' ::: medaka_mapped/!(*.md).bam

    COVPLOT=$(cat <<-'EOF'
    set term png;
    set output figs.'/'.sample.'_'.contig.'.png';
    set title sample.'-'.contig noenhanced;
    set xlabel 'bp';
    set ylabel 'coverage';
    unset key;
    stats file nooutput;
    plot [:STATS_records] file with dots;
EOF
    )

    echo "Plotting coverage depth."
    parallel -q --rpl "$FIRSTDIRECTORY" \
        gnuplot -e "sample='{m}'; contig='{/}'; file='{}'; figs='figs';$COVPLOT" \
        ::: stats/*/*
fi
shopt -u extglob

# Variance calling / sequence contribution?
echo "Done"
date
echo "Assemblies in medaka folder"
exit 0
