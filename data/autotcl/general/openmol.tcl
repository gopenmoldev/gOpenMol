############################################################################
# Filter to read a cli__line_input.data file 
#
# Copyright Leif Laaksonen 1997
#
# EXAMPLE:
#
#: label
#H2O
#: gaussian_basis_set_charge_centers
#O1  8.0  000.000000  000.000000    0.009000 /O_10.9s5p-4s2p
#H1  1.0     .000000    1.515263    1.058898 /H_2.4s-2s
#H2  1.0     .000000   -1.515263    1.058898 /H_2.4s-2s
#
# The OpenMol coordinates are in Atom Units!!!!!
#
proc lulOpenMoldata2XYZ { FileName } {

    set ConvUnits 0.52917715

    set File_p [open $FileName r]

    if {$File_p == 0} {
       puts "Can't open file '$FileName'"
       return
    }

    set Loop  0

    while { [gets $File_p line] >= 0 } {

# hook on ': label'
     if {[string match ": label*" $line]} {
# hook on ': gaussian_basis_set_charge_centers'
     } elseif {[string match ": gaussian_basis_set_charge_centers*" $line]} {
# start looping for the atom records 
            while { [gets $File_p line] >= 0 } {

		if {[string match ":*" $line]} break

            set AtomName    [lindex $line 0]
            set AtomCharge  [lindex $line 1]
            set AtomXC      [expr $ConvUnits * [lindex $line 2]]
            set AtomYC      [expr $ConvUnits * [lindex $line 3]]
            set AtomZC      [expr $ConvUnits * [lindex $line 4]]
            set AtomGbasis  [lindex $line 5]

            set XYZinfo($Loop) "$AtomName [format "%f %f %f" $AtomXC $AtomYC $AtomZC]"
            incr Loop
	    }
	}
    }

    close $File_p

    if {$Loop < 1} {
        lulErrorDialog "Could not find any atoms in the file '$FileName'"
        return;
    }

#    set FileSplit [split $FileName "."]
#    set OutFileName "[lindex $FileSplit 0].xyz"

     set OutFileName "$FileName.xyz"

# check if there is a file with the same name 

    if {[file exists $OutFileName]} {
        lulErrorDialog "The file with name '$OutFileName' already exists. Can't continue"
        return;
    }

    if {[file isdirectory $OutFileName]} {
            lulErrorDialog "ERROR - output is a directory"
            return
    }

    set Out_p [open $OutFileName w]
    if {$Out_p == 0} {
        lulErrorDialog "Can't open output file '$OutFileName'"
        return
    }

# first number of atoms
        puts $Out_p $Loop
    for {set i 0} { $i < $Loop } {incr i} {
# then the atoms
        puts $Out_p $XYZinfo($i)
    }

    close $Out_p

    return $OutFileName
}


