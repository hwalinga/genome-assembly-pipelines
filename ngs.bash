#!/bin/bash

###################
# PARSE ARGUMENTS #
###################

HELP=$(cat << 'EOF'
Usage:
--test:
    Show commands to run
--nomove:
    If source files are on a different parition than the current directory,
    they will first will be moved to this partition.
    If you have a reliable connection between the two drives, you can disable
    that default behavior with this flag.
    (Results still end up on the current directory regardless of this option.)
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly
--bacteria,--phage,--both:
    Set one of this option if you have bacteria or phage DNA.
    In this case the difference with a phage assembly is that it samples the
    data, so that the assembly is faster. There isn't any
EOF
)

echo "Parsing arguments"

NOMOVE=false
KEEP=false
BACTERIA=true
PHAGE=true
TEST=""
COVERAGE=false

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
        --coverage)
            COVERAGE=true
            ;;
        --)
            SOURCE1="$2"
            SOURCE2="$3"
            break
            ;;
        -*)
            >&2 echo "Unrecognized option ($1), but we are moving on."
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
    >&2 echo "Have not set both source files, exiting"
    echo "source1 is:"
    echo "$SOURCE1"
    echo "source2 is:"
    echo "$SOURCE2"
    exit 1
fi

echo "Arguments parsed:"
echo "SOURCE1 files: $SOURCE1"
echo "SOURCE2 files: $SOURCE2"

#########################
# CHECKING DEPENDENCIES #
#########################

program_missing () {
    $program=$1
    >&2 echo "$program not installed"
    >&2 echo "You can install with:"
    >&2 echo "conda -c bioconda install ${program%.*}"
    >&2 echo "or:"
    >&2 echo "brew install ${program%.*}"
    >&2 echo "Provided that one of these package managers are installed."
    exit 1
}

echo "Checking dependencies."
if [[ ! `command -v SOAPnuke` ]]; then
    >&2 echo "SOAPnuke not installed"
    >&2 echo "Download SOAPnuke from github and compile and add the executable to your \$PATH"
    exit 1
fi
for program in fastp spades.py parallel; do
    if [[ ! `command -v $program` ]]; then
        program_missing $program
    fi
done
if [[ $PHAGE == "true" ]] && [[ ! `command -v seqtk` ]]; then
    >&2 echo "Need seqtk for downsampling data."
    >&2 echo "There is not a strict requirement for this,"
    >&2 echo "So you can also assemble as if you assemble bacteria."
    program_missing $program
fi
if [[ $COVERAGE == "true" ]]; then
    echo "Checking "

    for program in fastp spades.py parallel; do
        if [[ ! `command -v $program` ]]; then
            program_missing $program
        fi
    done



echo "Dependencies met."

#########################################
# GNU Parallel helper functions/strings #
#########################################

# I am not writing a full blown load manager,
# but GNU Parallel can do its best.
CALCFREEPROC () {
    # Calculate the amount of free processing powers left (in # of processors)
    ps -eo pcpu | awk -v P=`nproc` 'NR!=1{S+=$1}END{printf "%.2f",P-S/100}'
}
CALCFREEPROCABS () {
    # report the integer value of the amount of processort available.
    printf '%.f' $(CALCFREEPROC)
}
CALCFREEMEM () {
    # Calculate the amount of RAM left (in GB)
    free --giga | awk 'NR==2{print $NF}'
}
RECOMFREEMEM () {
    # For each processor there should be the same amount of RAM.
    bc -l <<<"$(CALCFREEMEM) / $(CALCFREEPROC)"
}
# $TEST can be the --dry-run option.
PARALLEL () {
    parallel $TEST -j $(CALCFREEPROCABS) --memfree $(RECOMFREEMEM) --load 100% $*
}

# I can insert some Perl inside GNU Parallel,
# that's what these are doing:
# (No bioinformatics project is complete without some unreadable Perl.)
# What is found here is the common prefix of two files
# having the paired end reads.
COMMONPREFIX='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*[^_-]).*\0.*/\1`;$_=$1;'
COMMONPREFIXWITHDIR='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*/.*[^_-]).*\0.*/\1`;$_=$1;s:/:_:'

#################
# PERFORM CHECK #
#################

# Check if glob pattern expands.
if [[ `stat -t $SOURCE1` ]]; then
    >&2 echo "Glob pattern of first files not found ($SOURCE1)."
    exit 1
fi
if [[ `stat -t $SOURCE2` ]]; then
    >&2 echo "Glob pattern of second files not found ($SOURCE2)."
    exit 1
fi

# Check if the glob pattern expand to unique files.
WITHDIRECTORY=false  # This variable notes if we need the directory name to preserve uniqueness.
if [[ `ls $SOURCE1 $SOURCE2 | sed 's:.*/::' | sort | uniq -D` ]]; then
    # not unique, maybe if we add one directory more.
    if [[ `ls $SOURCE1 $SOURCE2 | perl -pe 's:.*/(?=.*/)::' | sort | uniq -D` ]]; then
        >&2 echo "Files found (including first directory) are not unique to each other"
        exit 1
    fi
    WITHDIRECTORY=true
fi

# Check if first and second files have a common substring.
if [[ $WITHDIRECTORY == "false" ]]; then
    if [[ ! $(parallel --rpl $COMMONPREFIX echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
        >&2 echo "There is no common substring between the two lists of sources files"
        exit 1
    fi
else
    if [[ ! $(parallel --rpl $COMMONPREFIXWITHDIR echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
        >&2 echo "There is no common substring between the two lists of the source files."
        exit 1
    fi
fi

echo Checks done.

# Check if source files are on a different disk, if so,
# it is better to move them to a folder on the local first,
# unless --nomove is specified.
MOVED=false
if [[ $NOMOVE == "false" ]] && [[ `stat --printf '%d' $(echo $SOURCE1 | head -1)` -ne `stat --printf '%d' ./` ]]; then
    MOVED=true
    echo "Moving files from remote to local 'raw'"
    mkdir -p raw
    if [[ WITHDIRECTORY == "true" ]]; then
        echo "Moving files"
        for i in $SOURCE1 $SOURCE2; do
            echo "Moving file $i"
            cp $i raw/$(perl -pe 's:.*/(.*)/(.*):\1_\2:'<<<$i)
        done
        BASESOURCE1=$(perl -pe 's:.*/(.*)/(.*):\1_\2:'<<<"$SOURCE1")
        BASESOURCE2=$(perl -pe 's:.*/(.*)/(.*):\1_\2:'<<<"$SOURCE2")
        echo "Moving done"
    else
        echo "Moving first files"
        cp $SOURCE1 raw
        echo "Moving second files"
        cp $SOURCE2 raw
        echo "Moving done"
        BASESOURCE1=$(sed 's:.*/::'<<<"$SOURCE1")
        BASESOURCE2=$(sed 's:.*/::'<<<"$SOURCE2")
    fi
fi

##############
# PROCESSING #
##############

echo "Starting fastp"

# Figuring out where the raw files are
mkdir -p fastp html json
if [[ $MOVED == "true" ]]; then
    RAWSOURCE1=raw/$BASESOURCE1
    RAWSOURCE2=raw/$BASESOURCE2
else
    RAWSOURCE1=$SOURCE1
    RAWSOURCE2=$SOURCE2
fi

if [[ $WITHDIRECTORY == "true" ]]; then
    COMMONPREFIXRAW=$COMMONPREFIXWITHDIR
    TARGET='{t} s:.*/(.*)/(.*):\1_\2:'
else
    COMMONPREFIXRAW=$COMMONPREFIX
    TARGET='{t} s:.*/::'
fi

# Now we can run fastp in parallel
PARALLEL --rpl $COMMONPREFIXRAW --rpl $TARGET \
    fastp -i {1} -o fastp/{1t} -I {2} -O fastp/{2t} \
    -5 -3 --correction -qualified_quality_phred 20 --lenght_required 30 \
    -j json/{1cp}.json -h html/{1cp}.html --report_title {1cp} \
    ::: $RAWSOURCE1 :::+ $RAWSOURCE2

echo "Finished fastp"
echo "Making pdf about quality."
mkdir -p pdf
PARALLEL wkhtmltopdf {} pdf/{/.} ::: html/*

if [[ $MOVED == "true" ]] && [[ $KEEP == "false" ]]; then
    echo "Removing raw"
    rm -rf raw
fi

echo "Starting SOAPNuke"

mkdir -p soapnuke
PARALLEL --rpl $COMMONPREFIX \
    SOAPnuke filter -1 {1} -2 {2} \
    -C {1/} -D {2/} -o soapnuke/{1cp} \
    ::: fastp/$BASESOURCE1 :::+ fastp/$BASESOURCE2

if [[ $KEEP == "false" ]]; then
    echo "Removing fastp (not the quality output)"
    rm -rf fastp
fi

if [[ $PHAGE == "true" ]]; then
    echo "Starting seqtk for phage analysis."
    mkdir -p seqtk
    PARALLEL seqtk sample {} 25000 ">" seqtk/{/.} ::: soapnuke/*/*.fq.gz
    if [[ $BACTERIA == "false" ]] && [[ $KEEP == "false" ]]; then
        echo "Removing SOAPnuke files"
        rm -rf soapnuke
    fi
fi

echo "Assembly"

if [[ $PHAGE == "true" ]]; then
    echo "Running spades.py for phages."
    mkdir -p spades_phage
    PARALLEL --rpl $COMMONPREFIX \
        mkdir -p spades_phage/{1cp} "&&" \
        spades.py -1 {1} -2 {2} --carefull -o spades_phage/{1cp} \
        ::: seqtk/$BASESOURCE1 :::+ seqtk/$BASESOURCE2
    if [[ $KEEP == "false" ]]; then
        echo "Removing seqtk files"
        rm -rf seqtk
    fi
fi

if [[ $BACTERIA == "true" ]]; then
    echo "Running spades.py for bacteria."
    mkdir -p spades_bac
    PARALLEL --rpl $COMMONPREFIX \
        mkdir -p spades_bac/{1cp} "&&" \
        spades.py -1 {1} -2 {2} --carefull -o spades_bac/{1cp} \
        ::: soapnuke/*/$BASESOURCE1 :::+ soapnuke/*/$BASESOURCE2
    if [[ $KEEP == "false" ]]; then
        echo "Removing SOAPnuke files"
        rm -rf soapnuke
    fi
fi

echo "Done"
exit 0
