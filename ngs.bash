#!/bin/bash

# parse arguments
#
# move to local preference (--nomove)
# keep for easy resume (--keep)
# --bacteria or --phage
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

WITHDIRECTORY=false  # This variable notes if we need the directory name to preserve uniqueness.
if [[ $(perl <<<$UNIQUETEST "" "$SOURCE1" "$SOURCE2") == "false" ]]; then
    # not unique, maybe if we add one directory
    if [[ $(perl <<<$UNIQUETEST "" "$SOURCE1" "$SOURCE2") == "false" ]]; then
        >&2 echo "Files found are not unique to each other"
        exit 1
    fi
    WITHDIRECTORY=true
fi

# Check if source files are on a different disk, if so, it is better to move
# them to a folder on the local first, unless --nomove is specified.
SOURCETEST=$(echo $SOURCE1 | head -1)
if [[ $NOMOVE == "false" ]] && [[ `stat --printf '%d' $SOURCETEST` -ne `stat --printf '%d' ./` ]]; then
    echo hi
fi

FREEMEM=$($FREEMEMCALC)

fastp
