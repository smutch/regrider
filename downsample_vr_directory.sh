#!/usr/bin/env bash

if [[ ! -x ./build/bin/regrider ]]; then
    echo "You need to compile the code first!"
    exit -1
fi

usage() {
    echo "
$0 usage:

    -b <dir>    Bin directory of hdf5 (location of h5ls and h5copy)
    -s <dir>    Source directory containing VR grids
    -d <dir>    Destination directory to place VR grids (will be created if necessary)
    -n <int>    New grid dimension (n x n x n)

    All options must be present.
"
}

while getopts ":hb:s:d:n:" opt; do
    case $opt in
        b)
            h5bin=${OPTARG}
            ;;
        s)
            source=${OPTARG}
            ;;
        d)
            destination=${OPTARG}
            ;;
        n)
            dim=${OPTARG}
            ;;
        h | *)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

h5ls=$h5bin/h5ls
h5copy=$h5bin/h5copy
downsample=./build/bin/regrider

if [[ ! -d $destination ]]; then
    mkdir -p $destination
fi

for file_in in $source/snap_*.hdf5; do
    file_out=$destination/$(basename $file_in)
    echo "$file_in -> $file_out"

    groups=$($h5ls $file_in | cut -d ' ' -f 1)
    for group in $groups; do
        if [[ $group != "PartType1" ]]; then
            $h5copy -i $file_in -o $file_out -s $group -d $group
        fi
    done

    $downsample -v $file_in -d $dim -o $file_out
done
