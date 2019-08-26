#!/bin/bash

###################
# PARSE ARGUMENTS #
###################

HELP=$(cat << EOF
Usage:
ngs.bash [--options] "SOURCE1_glob" "SOURCE2_glob"
Provide the source1 and source2 interleaved glob patterns quoted, so that
they will not be expanded before the program can read them.
--test:
    Show commands to run, but not execute them.
--nomove:
    If source files are on a different parition than the current directory,
    they will first will be moved to this partition.
    If you have a reliable connection between the two drives, you can disable
    that default behavior with this flag.
    (Results still end up on the current directory regardless of this option.)
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
    (NB. However, at the moment there are no checks in place that check
    if a part of the pipeline has already been completed in a previous run.)
--bacteria,--phage,--both:
    Set one of this option if you have bacteria or phage DNA.
    In this case the difference with a phage assembly is that it samples the
    data, so that the assembly is faster. There isn't any difference further
    at the moment.
--nocov:
    Do not plot the coverage plots.
--help,-h
    Plot this help and exit.
EOF
)

echo "========"
echo "Parsing arguments"
echo "========"

NOMOVE=false
KEEP=false
BACTERIA=true
PHAGE=true
TEST=""
COVERAGE=true
PDF=true
ERRORLOG="error.log"

while [[ -n "$1" ]]; do
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
        --bacteria)
            PHAGE=false
            BACTERIA=true
            ;;
        --phage)
            PHAGE=true
            BACTERIA=false
            ;;
        --both)
            PHAGE=true
            BACTERIA=true
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
        --)
            SOURCE1="$2"
            SOURCE2="$3"
            break
            ;;
        -*)
            >&2 echo "Unrecognized option"
            >&2 echo "Exiting."
            exit 1
            ;;
        *)
            if [[ -z "$SOURCE1" ]]; then
                SOURCE1="$1"
            else
                SOURCE2="$1"
            fi
            ;;
    esac
    shift
done

if [[ -z "$SOURCE1" ]] || [[ -z "$SOURCE2" ]]; then
    >&2 echo "Have not set both source files"
    echo "source1 is:"
    echo "$SOURCE1"
    echo "source2 is:"
    echo "$SOURCE2"
    >&2 "Exiting"
    exit 1
fi

echo "========"
echo "Arguments parsed:"
echo "SOURCE1 files: $SOURCE1"
echo "SOURCE2 files: $SOURCE2"
if [[ $COVERAGE == "true" ]]; then
    echo "Run with making coverage plots."
fi
if [[ $TEST == "--dry-run" ]]; then
    echo "Only output what commands are run (dry run)."
fi
if [[ $KEEP == "true" ]]; then
    echo "Keep inbetween results."
fi

# Create symlink if there are spaces in filename problems.
echo "========"
echo "Checking spaces in filename problems."
if [[ "$SOURCE1" =~ .*\ .* ]] || [[ "$SOURCE2" =~ .*\ .* ]]; then
    if [[ "$PWD" =~ .*\ .* ]]; then  # Path contains spaces.
        >&2 "The provided paths and the current directory contains spaces"
        >&2 "There is no workaround I can apply here."
        >&2 "Change the provided paths or the current directory so that"
        >&2 "at least one of these does not contain any spaces."
        >&2 "Exiting"
        exit 1
    fi
    BASENAME1="${SOURCE1##*/}"
    DIRNAME1="${SOURCE1%/*}"
    BASENAME2="${SOURCE2##*/}"
    DIRNAME2="${SOURCE2%/*}"
    if [[ "$DIRNAME1" =~ .*\*.* ]] || [[ "$DIRNAME2" =~ .*\*.* ]] || [[ "$BASENAME1" =~ .*\ .* ]] || [[ "$BASENAME2" =~ .*\ .* ]]; then
        >&2 "Combining glob patterns with spaces is not possible."
        >&2 "Rename or move your files."
        >&2 "Exiting"
        exit 1
    fi
    ln -s "$DIRNAME1" symlink_files_1
    ln -s "$DIRNAME2" symlink_files_2
    SOURCE1="symlink_files_1/$BASENAME1"
    SOURCE2="symlink_files_2/$BASENAME2"
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
if [[ ! `command -v wkhtmltopdf` ]]; then
    >&2 echo "wkhtmltopdf not installed."
    >&2 echo "Really not a problem, only used for html to pdf conversion."
    >&2 echo "You can just open the html instead."
    PDF=false
fi
echo "Checking dependencies for assembly."
if [[ ! `command -v SOAPnuke` ]]; then
    >&2 echo "SOAPnuke not installed"
    >&2 echo "Download SOAPnuke from github and compile and add the executable to your \$PATH"
    >&2 echo "Exiting"
    exit 1
fi
for program in fastp spades.py parallel; do
    if [[ ! `command -v $program` ]]; then
        program_missing $program
        >&2 echo "Exiting"
        exit 1

    fi
done
if [[ $PHAGE == "true" ]] && [[ ! `command -v seqtk` ]]; then
    >&2 echo "Need seqtk for downsampling data, when assembling phages"
    >&2 echo "There is not a strict requirement for this, but recommended."
    >&2 echo "But you can also assemble as if you assemble bacteria."
    >&2 echo ""
    program_missing $program
    >&2 echo "Exiting"
    exit 1
fi
if [[ $COVERAGE == "true" ]]; then
    echo "Checking dependencies for coverage plots."
    echo "If these dependencies fail, you can still assemble with the --nocov option."
    for program in samtools gnuplot minimap2; do
        if [[ ! `command -v $program` ]]; then
            program_missing $program
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

# I can inject some Perl inside GNU Parallel,
# that's what these are doing:
# (No bioinformatics project is complete without some unreadable Perl.)
# What is found here is the common prefix of two files
# that are having the paired end reads.
# The second line will find the common prefix including the first directory,
# replacing the directory seperator (/) with an underscore (_).
COMMONPREFIX='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*[^_-]).*\0.*/\1`;$_=$1;'
COMMONPREFIXWITHDIR='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*/.*[^_-]).*\0.*/\1`;$_=$1;s:/:_:'
# Define common prefix if first directory is needed:
WITHDIRECTORYREGEX='s:.*/(.*)/(.*):\1_\2:'
# Perl regex to find the first directory:
# FIRSTDIRECTORY='{m} s:*?/::;s:/.*::;'
FIRSTDIRECTORY='{m} s:.*/(?=.*/)::;s:/.*::;'

#################
# PERFORM CHECK #
#################

# Check if glob pattern expands.
if [[ ! `stat -t $SOURCE1` ]]; then
    >&2 echo "Glob pattern of first files not found ($SOURCE1)."
    >&2 echo "Exiting"
    exit 1
fi
if [[ ! `stat -t $SOURCE2` ]]; then
    >&2 echo "Glob pattern of second files not found ($SOURCE2)."
    >&2 echo "Exiting"
    exit 1
fi

# Check if the glob pattern expand to unique files.
WITHDIRECTORY=false  # This variable notes if we need the directory name to preserve uniqueness.
if [[ `ls $SOURCE1 $SOURCE2 | sed 's:.*/::' | sort | uniq -D` ]]; then
    # not unique, maybe if we add one directory more.
    if [[ `ls $SOURCE1 $SOURCE2 | perl -pe 's:.*/(?=.*/)::' | sort | uniq -D` ]]; then
        >&2 echo "Files found (including first directory) are not unique to each other"
        >&2 echo "Exiting"
        exit 1
    fi
    WITHDIRECTORY=true
fi

# Check if first and second files have a common substring.
if [[ $WITHDIRECTORY == "false" ]]; then
    if [[ ! $(parallel --rpl '$COMMONPREFIX' echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
        >&2 echo "There is no common substring between the two lists of sources files"
        >&2 echo "Exiting"
        exit 1
    fi
else
    if [[ ! $(parallel --rpl '$COMMONPREFIXWITHDIR' echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
        >&2 echo "There is no common substring between the two lists of the source files."
        >&2 echo "Exiting"
        exit 1
    fi
fi

echo "Checks done."
echo "========"

# Check if source files are on a different disk, if so,
# it is better to move them to a folder on the local first,
# unless --nomove is specified.
MOVED=false
if [[ $NOMOVE == "false" ]] && [[ `stat --printf '%d' $(echo $SOURCE1 | head -1)` -ne `stat --printf '%d' ./` ]]; then
    MOVED=true
    echo "Copying files from remote to local 'raw'"
    $MKDIR raw
    if [[ WITHDIRECTORY == "true" ]]; then
        echo "Copying files"
        for i in $SOURCE1 $SOURCE2; do
            echo "Copying file $i"
            $CP $i raw/$(perl -pe '$WITHDIRECTORYREGEX'<<<$i)
        done
        echo "Copying done"
    else
        echo "Copying first files"
        $CP $SOURCE1 raw
        echo "Copying second files"
        $CP $SOURCE2 raw
        echo "Copying done"
    fi
fi

if [[ WITHDIRECTORY == "true" ]]; then
    BASESOURCE1=$(perl -pe '$WITHDIRECTORYREGEX'<<<"$SOURCE1")
    BASESOURCE2=$(perl -pe '$WITHDIRECTORYREGEX'<<<"$SOURCE2")
else
    BASESOURCE1=$(sed 's:.*/::'<<<"$SOURCE1")
    BASESOURCE2=$(sed 's:.*/::'<<<"$SOURCE2")
fi

##############
# PROCESSING #
##############

echo "Starting fastp"
echo "========"

# Figuring out where the raw files are
$MKDIR fastp html json
if [[ $MOVED == "true" ]]; then
    RAWSOURCE1=raw/$BASESOURCE1
    RAWSOURCE2=raw/$BASESOURCE2
else
    RAWSOURCE1=$SOURCE1
    RAWSOURCE2=$SOURCE2
fi

if [[ $WITHDIRECTORY == "true" ]]; then
    COMMONPREFIXRAW=$COMMONPREFIXWITHDIR
    TARGET="{t} $WITHDIRECTORYREGEX"
else
    COMMONPREFIXRAW=$COMMONPREFIX
    TARGET='{t} s:.*/::'
fi

# Now we can run fastp in parallel
PARALLEL --rpl "$COMMONPREFIXRAW" --rpl "$TARGET" \
    fastp -i {1} -o fastp/{1t} -I {2} -O fastp/{2t} \
    -5 -3 --correction --qualified_quality_phred 20 --length_required 30 \
    -j json/{1cp}.json -h html/{1cp}.html --report_title {1cp} \
    ::: $RAWSOURCE1 :::+ $RAWSOURCE2

echo "Finished fastp"
echo "========"
if [[ $PDF == "true" ]]; then
    echo "Making pdf about quality, this is the same as the HTML, just not interactive."
    $MKDIR pdf
    PARALLEL wkhtmltopdf {} pdf/{/.}.pdf ::: html/*
fi

if [[ $MOVED == "true" ]] && [[ $KEEP == "false" ]]; then
    echo "Removing raw"
    rm -rf raw
fi

echo "Starting SOAPNuke"
echo "========"

$MKDIR soapnuke
PARALLEL --rpl "$COMMONPREFIX" \
    SOAPnuke filter -1 {1} -2 {2} \
    -C {1/} -D {2/} -o soapnuke/{1cp} \
    ::: fastp/$BASESOURCE1 :::+ fastp/$BASESOURCE2

if [[ $KEEP == "false" ]]; then
    echo "Removing fastp (not the quality output)."
    rm -rf fastp
fi

if [[ $PHAGE == "true" ]]; then
    echo "Starting seqtk for phage analysis."
    $MKDIR seqtk
    PARALLEL seqtk sample {} 25000 ">" seqtk/{/.} ::: soapnuke/*/*.fq.gz
    if [[ $BACTERIA == "false" ]] && [[ $KEEP == "false" ]]; then
        echo "Removing SOAPnuke files"
        rm -rf soapnuke
    fi
fi

echo "Assembly"
echo "========"
phage_suffix="_phage"
bac_suffix="_bac"

targets=""
if [[ $BACTERIA == "true" ]]; then
    targets+=" $bac_suffix"
fi
if [[ $PHAGE == "true" ]]; then
    targets+=" $phage_suffix"
fi
for t in $targets; do

    if [[ $t == $phage_suffix ]]; then
        # phage corrected raw sources
        corrawsource1=seqtk/${BASESOURCE1%.gz}
        corrawsource2=seqtk/${BASESOURCE2%.gz}
    else
        # bacteria corrected raw sources
        corrawsource1=soapnuke/*/${BASESOURCE1%.gz}.gz
        corrawsource2=soapnuke/*/${BASESOURCE2%.gz}.gz
    fi

    echo "Running spades.py for ${t#_}."
    $MKDIR spades$t
    PARALLEL --rpl "$COMMONPREFIX" \
        mkdir -p spades$t/{1cp} "&&" \
        spades.py -1 {1} -2 {2} --careful -o spades$t/{1cp} \
        ::: $corrawsource1 :::+ $corrawsource2

    # If assembly failed, it is easier for the implementation to just
    # have an empty file instead of no file:
    touch -a spades$t/*/contigs.fasta
    for assm in spades$t/*/contigs.fasta; do
        if [[ ! -s $assm ]]; then
            echo "Seems that the assembly of $assm has failed." |
                tee -a $ERRORLOG | cat 1>&2
            echo "Here are the last 40 lines of spades.log" >> $ERRORLOG
            tail -40 ${ass%/*}/spades.log >> $ERRORLOG
        fi
    done

    if [[ $COVERAGE == "false" ]] && [[ $KEEP == "false" ]]; then
        if [[ $t == $phage_suffix ]]; then
            echo "Removing seqtk files (for phage assembly)"
            rm -rf seqtk
        else
            echo "Removing SOAPnuke files (for bacteria assembly)"
            rm -rf soapnuke
        fi
    fi
done

# Coverage plots.
if [[ $COVERAGE == "true" ]]; then
    echo "Making assembly plots."
    targets=""
    if [[ $BACTERIA == "true" ]]; then
        targets+=" $bac_suffix"
    fi
    if [[ $PHAGE == "true" ]]; then
        targets+=" $phage_suffix"
    fi
    for t in $targets; do

        if [[ $t == $phage_suffix ]]; then
            # phage corrected raw sources
            corrawsource1=seqtk/${BASESOURCE1%.gz}
            corrawsource2=seqtk/${BASESOURCE2%.gz}
        else
            # bacteria corrected raw sources
            corrawsource1=soapnuke/*/${BASESOURCE1%.gz}.gz
            corrawsource2=soapnuke/*/${BASESOURCE2%.gz}.gz
        fi

        echo "========"
        echo "Working on ${t#_}"
        $MKDIR {mapped,stats,figs}$t

        echo "Mapping to raw reads."
        PARALLEL --rpl "$FIRSTDIRECTORY" \
            test -s {1} "&&" \
            minimap2 -ax sr {1} {2} {3} "|" \
            awk -F\\t -v OFS=\\t "'{\$1=substr(\$1,1,251)}1'" \
            ">" mapped$t/{1m}.sam \
            ::: spades$t/*/contigs.fasta :::+ $corrawsource1 :::+ $corrawsource2

        echo "Sorting and indexing of sam files."
        PARALLEL --plus 'samtools sort {} > {.}.bam && samtools index {.}.bam' \
            ::: mapped$t/*.sam

        echo "Calculating depth."
        parallel --dry-run $MKDIR stats$t/{/.} ::: mapped$t/*.sam | sh
        PARALLEL 'samtools depth -a {} | awk "{print \$3 > \"stats'$t'/{/.}/\"\$1}"' \
            ::: mapped$t/*.bam

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
        gnuplot -e "sample='{m}'; contig='{/}'; file='{}'; figs='figs$t';$COVPLOT" \
        ::: stats$t/*/*
    done
fi

echo "========"
echo "Done"
echo "========"
exit 0
