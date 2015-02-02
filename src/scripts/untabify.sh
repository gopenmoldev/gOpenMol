#! /bin/sh
##############################################################################
#                           Copyright (c) 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################
# usage:
#     untabify.sh [file1 [directory2 [file3 [...]]]]
#
#     Directory refers to directory and to all subdirectries.
#     Convert all tabulars to 4 spaces in all C and C++ files
case $# in
0) set . ;;
esac

# Expand directories
for arg
do
    if test -d "$arg"
    then
	ifs="$IFS"
	# Split on new line.
	IFS='
'
	set x "$@" `find "$arg" -type f '(' \
	    -name '*.[ch]' -o \
	    -name '*.[ch]pp' -o \
	    -name '*.[ch][ch]' \
	    ')' -follow`
	IFS="$ifs"
    else
	set x "$@" "$arg"
    fi
    shift # x
    shift # $arg
done

# Untabify files.
for file
do
    if tr -dc '\t' <"$file" | grep . >/dev/null
    then
	echo -n "untabifying \`${file}' ..."
	expand -t 4 <"$file" >"$file.tmp" &&
	    diff -b "$file" "$file.tmp" &&
	    cat "$file.tmp" >"$file" &&
	    rm "$file.tmp" &&
	    echo "done"
    fi
done
