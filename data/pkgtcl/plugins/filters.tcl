###########################################################################
# PROC
namespace eval gom::FileFilters::AtomCoordinates::Chem3d {
proc Read {FileName Append} {
    puts "Reading file...."
    set NumAtoms [Chem3DImport $FileName]
    if {$NumAtoms > 0} {
	define structure "$FileName" $NumAtoms $Append
	set StrucNum [show molstructures]
	for {set i 1} {$i<=$NumAtoms} {incr i} {
	    set Xc [getChem3DVec $i 0]
	    set Yc [getChem3DVec $i 1]
	    set Zc [getChem3DVec $i 2]
	    set Atm [getChem3DName $StrucNum $i]
	    set Seg "USER"
	    set Res [getChem3DName $StrucNum $i]
	    set Ch "0.0"
	    set ResN 1
	    define atom    label    $Atm        $i $StrucNum
	    define residue label    $Res        $i $StrucNum
	    define segment label    $Seg        $i $StrucNum
	    define atom coordinates $Xc $Yc $Zc $i $StrucNum
	    define atom charge      $Ch         $i $StrucNum
	    define atom resnumber   $ResN       $i $StrucNum
	}
    }
    puts "Releasing space."
    releasec3d $NumAtoms
    puts Done!
}
}

###########################################################################
# PROC
namespace eval gom::FileFilters::AtomCoordinates::Spartan {
proc Read {FileName Append} {
    puts "Reading file...."
    set NumAtoms [SpartanImport $FileName]
    if {$NumAtoms > 0} {
	define structure "$FileName" $NumAtoms $Append
	set StrucNum [show molstructures]
	for {set i 1} {$i<=$NumAtoms} {incr i} {
	    set Xc [getChem3DVec $i 0]
	    set Yc [getChem3DVec $i 1]
	    set Zc [getChem3DVec $i 2]
	    set Atm [getChem3DName $StrucNum $i]
	    set Seg "USER"
	    set Res [getChem3DName $StrucNum $i]
	    set Ch "0.0"
	    set ResN 1
	    define atom    label    $Atm        $i $StrucNum
	    define residue label    $Res        $i $StrucNum
	    define segment label    $Seg        $i $StrucNum
	    define atom coordinates $Xc $Yc $Zc $i $StrucNum
	    define atom charge      $Ch         $i $StrucNum
	    define atom resnumber   $ResN       $i $StrucNum
	}
    }
    puts "Releasing space."
    releasec3d $NumAtoms
    puts Done!
}
}

namespace eval FiltersPlugin {

proc Enable {} {
    array set ::gom::FileFilters::AtomCoordinates::coordTypes {
	chem3d  {{1 0} Chem3D  {Chem 3D files} {*}}
	spartan {{1 0} Spartan {Spartan files} {*}}
    }
}

gom::RegisterPlugin "filters" "Chem3D/Spartan coordinate filter"

}
