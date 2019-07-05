#!/bin/bash

# parse arguments
#
# move to local preference (--nomove)
# keep for easy resume (--keep)
# --bacteria or --phage, or --both (last one applies).
# two globs for 1 and 2, (positional)
#
#
# check deps (SOAPNuke, seqtk, fastp, spades.py, parallel)
#
# check if on a remote (prefer move)
# check if filenames have uniqueness
#
# run parallel

# SOURCE1
# TEST (--dry-run or "")
# SOURCE2
# NOMOVE
# TYPE

echo "Checking dependencies."
if [[ ! `command -v SOAPnuke` ]]; then
    echo "SOAPnuke not installed"
    echo "Download SOAPnuke from github and compile and add the executable to your \$PATH"
for program in seqtk fastp spades.py parallel; do
    if [[ ! `command -v $program` ]]; then
        echo "$program not installed"
        echo "You can install with:"
        echo "conda -c bioconda install ${program%.*}"
        echo "or:"
        echo "brew install ${program%.*}"
        echo "Provided that one of these package managers are installed."
    fi
done

CALCFREEPROC="ps -eo pcpu | awk -v P=`nproc` 'NR!=1{S+=\$1}END{printf \"%.2f\",P-S/100}'"
CALCFREEMEM="free --giga --total | tail -1 | awk '{print \$4}'"
RECOMFREEMEM="printf '%.f' \$(( \$($CALCFREEMEM) / \$($CALCFREEPROC)))"
PARALLEL="parallel $TEST -j \$($CALCFREEPROC) --memfree \$($RECOMFREEMEM)"
COMMONPREFIX='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*[^_-]).*\0.*/\1`;$_=$1;'
COMMONPREFIXWITHDIR='{cp} "$arg[1]\0$arg[2]"=~m`^.*/(.*/.*[^_-]).*\0.*/\1`;$_=$1;'

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
$PARALLEL --rpl $COMMONPREFIXRAW --rpl $TARGET \
    fastp -i {1} -o fastp/{1t} -I {2} -O fastp/{2t} \
    -5 -3 --correction -qualified_quality_phred 20 --lenght_required 30 \
    -j json/{1cp}.json -h html/{1cp}.html --report_title {1cp} \
    ::: $RAWSOURCE1 :::+ $RAWSOURCE2

echo "Finished fastp"
echo "Making pdf about quality."
mkdir -p pdf
$PARALLEL wkhtmltopdf {} pdf/{/.} ::: html/*

if [[ $MOVED == "true" ]] && [[ $KEEP == "false" ]]; then
    echo "Removing raw"
    rm -rf raw
fi

echo "Starting SOAPNuke"

mkdir -p soapnuke
$PARALLEL --rpl $COMMONPREFIX \
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
    $PARALLEL seqtk sample {} 25000 ">" seqtk/{/.} ::: soapnuke/*/*.fq.gz
    if [[ $BACTERIA == "false" ]] && [[ $KEEP == "false" ]]; then
        echo "Removing SOAPnuke files"
        rm -rf soapnuke
    fi
fi

echo "Assembly"

if [[ $PHAGE == "true" ]]; then
    echo "Running spades.py for phages."
    mkdir -p spades_phage
    $PARALLEL --rpl $COMMONPREFIX \
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
    $PARALLEL --rpl $COMMONPREFIX \
        mkdir -p spades_bac/{1cp} "&&" \
        spades.py -1 {1} -2 {2} --carefull -o spades_bac/{1cp} \
        ::: soapnuke/*/$BASESOURCE1 :::+ soapnuke/*/$BASESOURCE2
    if [[ $KEEP == "false" ]]; then
        echo "Removing SOAPnuke files"
        rm -rf soapnuke
    fi
fi

echo "Done"
