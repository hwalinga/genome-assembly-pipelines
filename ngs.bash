#!/bin/bash

# parse arguments
#
# move to local preference (--nomove)
# keep for easy resume (--keep)
# --bacteria or --phage, or --both (last one applies).
# two globs for 1 and 2, (positional)
#
#
# check deps (SOAPNuke, seqtk, fastp, spades.py)
#
# check if on a remote (prefer move)
# check if filenames have uniqueness
#
# run parallel

SOURCE1
SOURCE2
NOMOVE
TYPE

CALCFREEMEM="free --giga --total | tail -1 | awk '{print $4}'"
CALCFREEPROC="ps -eo pcpu | awk -v P=`nproc` 'NR!=1{S+=$1}END{printf \"%.2f\",P-S/100}'"

# This perl program checks if the files expanded by the glob patterns are unique.
# The first argument is a regex that removes any leading directories.
# Following arguments are the files.
UNIQUETEST=$(cat<<'EOF'
my @files = ((glob $ARGV[1]), (glob $ARGV[2]));
map{s:$ARGV[0]::}@files;
my @unique = keys{map{$_=>1}@files};
print $#files == $#unique ? "true" : "false";
EOF
)

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
    if [[ ! $(parallel --rpl '{cp} "@arg[1]\0@arg[2]"=~m`^.*/(.*).*\0.*/\1`;$_=$1;' echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
        >&2 echo "There is no common substring between the two lists of sources files"
        exit 1
    fi
else
    if [[ ! $(parallel --rpl '{cp} "@arg[1]\0@arg[2]"=~m`.*/(.*/.*).*\0.*/\1`;$_=$1;' echo {1cp} ::: $SOURCE1 :::+ $SOURCE2) ]]; then
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
    mkdir raw
    if [[ WITHDIRECTORY == "true" ]]; then
        echo "Moving files"
        for i in $SOURCE1 $SOURCE2; do
            echo "Moving file $i"
            cp $i raw/$(perl -pe 's:.*/(.*)/(.*):\1_\2:'<<<$i)
        done
    else
        echo "Moving first files"
        cp $SOURCE1 raw
        echo "Moving second files"
        cp $SOURCE2 raw
        echo "Moving done"
    fi
fi



FREEMEM=$($FREEMEMCALC)

fastp
