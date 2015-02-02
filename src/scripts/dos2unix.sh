#! /bin/sh
##############################################################################
#			    Copyright (c) 2003 by:
#	   Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#		       Confidential unpublished property of 
#			       Leif Laaksonen	
#			     All rights reserved
#
#	   Coded by: Eero Häkkinen
##############################################################################
# usage:
#     dos2unix.sh file1 [file2 [file3 [...]]]
#
#     Directory refers to directory and to all subdirectries.
#     Removes all carrier returns from text files.
#
#     Converts all files in all subdirectories mathing one of patterns:
#	*.c, *.cc, *.cpp,			(C and C++ source)
#	*.h, *.hh, *.hpp,			(C and C++ headers)
#	*.f, *.for, *.fpp, *.f77, *.f90, *.f95	(Fortran)
#	*.tcl, *.sh, *.csh, *.bat,		(script files)
#	Makefile, *.ac, *.mk, *.m4		(makefiles etc.)
#	*.html, *.htm, *.rtf, *.tex, *.txt,
#	*.terms					(text files)
#	*.dic, *.out, *.xmol, *.xyz, *.crd,
#	*.gom					(ASCII data files)
case $# in
0) set . ;;
esac

# Expand directories.
for arg
do
    if test -d "$arg"
    then
	ifs="$IFS"
	# Split on new line.
	IFS='
'
	set x "$@" `find "\$arg" -type f '(' \
	    -name '*.[cfh]' -o \
	    -name '*.[cfh]pp' -o \
	    -name '*.[ch][ch]' -o \
	    -name '*.for' -o \
	    -name '*.f[0-9][0-9]' -o \
	    -name '*.tcl' -o \
	    -name '*.sh' -o \
	    -name '*.csh' -o \
	    -name '*.bat' -o \
	    -name 'Makefile' -o \
	    -name '*.ac' -o \
	    -name '*.mk' -o \
	    -name '*.m4' -o \
	    -name '*.html' -o \
	    -name '*.htm' -o \
	    -name '*.rtf' -o \
	    -name '*.tex' -o \
	    -name '*.txt' -o \
	    -name '*.terms' -o \
	    -name '*.dic' -o \
	    -name '*.out' -o \
	    -name '*.xmol' -o \
	    -name '*.xyz' -o \
	    -name '*.crd' -o \
	    -name '*.gom' \
	    ')' -follow`
    else
	set x "$@" "$arg"
    fi
    shift # x
    shift # $arg
done

# Convert to Unix text format.
for file
do
    if tr -dc '\r' <"$file" | grep . >/dev/null
    then
	echo -n "converting \`${file} ..."
	tr -d '\r' <"$file" >"$file.tmp" &&
	    diff -bu "$file" "$file.tmp" &&
	    cat "$file.tmp" >"$file" &&
	    rm "$file.tmp" &&
	    echo "done"
    fi
done
