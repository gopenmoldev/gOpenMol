namespace eval FCHKPlugin {

proc gridsave { id fn} {
	savegrid [show addr cont $id] $fn
}

proc loadmolecule {fn} {
	import coord gaussian $fn
	readfchk $fn
}


proc fchk {} {
	global fchkLoaded
	global fchkGridCalc
	global fchkOrb
	global fchkSpin
	global fchkFilename
	toplevel .fchkw
	frame .fchkw.filestuff
	frame .fchkw.gridstuff
	frame .fchkw.gridstuff.limits
	frame .fchkw.gridstuff.mesh
	frame .fchkw.gridoptions
	frame .fchkw.orbcmd
	frame .fchkw.orbcmd.orbitals
	frame .fchkw.orbcmd.commands

	entry .fchkw.filestuff.filename -width 32
	entry .fchkw.gridstuff.limits.xmin -width 10
	entry .fchkw.gridstuff.limits.xmax -width 10
	entry .fchkw.gridstuff.limits.ymin -width 10
	entry .fchkw.gridstuff.limits.ymax -width 10
	entry .fchkw.gridstuff.limits.zmin -width 10
	entry .fchkw.gridstuff.limits.zmax -width 10
	entry .fchkw.gridstuff.mesh.nx -width 5
	entry .fchkw.gridstuff.mesh.ny -width 5
	entry .fchkw.gridstuff.mesh.nz -width 5
	entry .fchkw.gridstuff.limits.scalf -width 5
	label .fchkw.gridstuff.grid -text "Grid Parameters"
	label .fchkw.gridstuff.limits.scallab -text "Scale box by"
	label .fchkw.gridstuff.mesh.nxlab -text "X points"
	label .fchkw.gridstuff.mesh.nylab -text "Y points"
	label .fchkw.gridstuff.mesh.nzlab -text "Z points"
	label .fchkw.gridstuff.limits.xsize -text "X Extents"
	label .fchkw.gridstuff.limits.ysize -text "Y Extents"
	label .fchkw.gridstuff.limits.zsize -text "Z Extents"
	label .fchkw.filestuff.filenamelabel -text "Current FCHK file"
	radiobutton .fchkw.gridoptions.orb -variable orb -value 1
	label .fchkw.gridoptions.orblab -text "Use selected orbital"
	radiobutton .fchkw.gridoptions.dens -variable orb -value 2
	label .fchkw.gridoptions.denlab -text "Use Density"
	radiobutton .fchkw.gridoptions.lplc -variable orb -value 3
	label .fchkw.gridoptions.lplclab -text "Use Laplacian of Density"
	radiobutton .fchkw.gridoptions.pspin -variable orb -value 4
	label .fchkw.gridoptions.spinlab -text "Use Spin Density Difference"
	radiobutton .fchkw.gridoptions.spind -variable orb -value 5
	label .fchkw.gridoptions.spindlab -text "Use Spin Density"
	radiobutton .fchkw.gridoptions.mesp -variable orb -value 6
	label .fchkw.gridoptions.mesplab -text "Use Electrostatic Potential"
	label .fchkw.gridoptions.spin -text "Spin to Use"
	label .fchkw.gridoptions.spina -text "Alpha" 
	radiobutton .fchkw.gridoptions.alpha -variable spin -value 0 -command {
		set	OrbEList [listfchkenergies $spin]
	}
	label .fchkw.gridoptions.spinb -text "Beta"
	radiobutton .fchkw.gridoptions.beta -variable spin -value 1 -command {
		set OrbEList [listfchkenergies $spin]
	}
	label .fchkw.orbcmd.orbitals.energylab -text "Orbital Energies"
	listbox .fchkw.orbcmd.orbitals.orblist -yscrollcommand {.fchkw.orbcmd.orbitals.yscroll set} -listvariable OrbEList
	scrollbar .fchkw.orbcmd.orbitals.yscroll -orient vertical -command {.fchkw.orbcmd.orbitals.orblist yview}
	button .fchkw.filestuff.getfile -text Browse -command {set fn [tk_getOpenFile -filetypes {{"Formatted checkpoint files" {.fchk .fch}}}] ;  .fchkw.filestuff.filename insert 0 "$fn" }
	button .fchkw.filestuff.openfile -text Load -command {
		set fn [.fchkw.filestuff.filename get]
		if {(([string compare $fn ""] != 0) && ($fchkLoaded == 0))} then {
			FCHKPlugin::loadmolecule $fn
			set tlist [getparms 0]
			.fchkw.gridstuff.mesh.nx insert 0 $tlist
			set tlist [getparms 1]
			.fchkw.gridstuff.mesh.ny insert 0 $tlist
			set tlist [getparms 2]
			.fchkw.gridstuff.mesh.nz insert 0 $tlist
			set tlist [getparms 3]
			.fchkw.gridstuff.limits.xmin insert 0 [lindex $tlist 0]
			.fchkw.gridstuff.limits.ymin insert 0 [lindex $tlist 1]
			.fchkw.gridstuff.limits.zmin insert 0 [lindex $tlist 2]
			set tlist [getparms 4]
			.fchkw.gridstuff.limits.xmax insert 0 [lindex $tlist 0]
			set tlist [getparms 5]
			.fchkw.gridstuff.limits.ymax insert 0 [lindex $tlist 1]
			set tlist [getparms 6]
			.fchkw.gridstuff.limits.zmax insert 0 [lindex $tlist 2]
			set OrbEList [listfchkenergies $spin]
			set orb fchkOrb
			set fchkLoaded 1
		}
	}
	button .fchkw.orbcmd.commands.calcgrid -text "Calculate Grid" -command {
		if {$orb == 1} {
			if {$spin == 1} {
				newfillgrid [.fchkw.orbcmd.orbitals.orblist curselection] "" b
			} else {
				newfillgrid [.fchkw.orbcmd.orbitals.orblist curselection]
			}
		} else {
			if {$orb == 3} {
				newfillgrid -1 "" a l
			} else {
				if {$orb == 4} {
					newfillgrid -1 "" d s
				} else {
					if {$orb == 5} {
						if {$spin == 1} {
							newfillgrid -1 "" b s
						} else {
							newfillgrid -1 "" a s
						}
					} else {
						if {$orb == 6 } {
							newfillgrid -1 "" a m
						} else {
							newfillgrid -1
						}
					}
				}
			}
		}
		set fchkGridCalc 1
	}
	button .fchkw.orbcmd.commands.savegrid -text "Save Grid" -command {
		if {$fchkGridCalc == 1} {
			toplevel .fchksavecontour
			entry .fchksavecontour.filename -width 32
			label .fchksavecontour.filetext -text "Filename"
			button .fchksavecontour.browse -text "Browse" -command {
				set fn [tk_getSaveFile -defaultextension .plt -filetypes {{"Grid files" {.plt .PLT}}}]
				.fchksavecontour.filename insert 0 $fn
			}
			entry .fchksavecontour.contourid 
			label .fchksavecontour.contidtext -text "Contour to save"
			.fchksavecontour.contourid insert 0 [show contour defined]
			button .fchksavecontour.save -text "Save Contour" -command {
				set fn [.fchksavecontour.filename get]
				set ffn [lindex $fn 0]
				if {[string compare $ffn ""]!=0} then {
					set fn [.fchksavecontour.contourid get]
					set cid [lindex $fn 0]
					FCHKPlugin::gridsave  $cid $ffn
					destroy .fchksavecontour
				}
			}
			button .fchksavecontour.abort -text "Cancel" -command {destroy .fchksavecontour}
			grid .fchksavecontour.filetext -row 1 -column 1
			grid .fchksavecontour.filename -row 1 -column 2
			grid .fchksavecontour.browse -row 1 -column 3
			grid .fchksavecontour.contidtext -row 2 -column 1
			grid .fchksavecontour.contourid -row 2 -column 2
			grid .fchksavecontour.save -row 3 -column 2
			grid .fchksavecontour.abort -row 3 -column 3
		}
	}
	button .fchkw.filestuff.freemem -text "Free FCHK Space" -command {
		if {$fchkLoaded == 1} {
			clearmem
			set fchkLoaded 0
			set fchkOrb 2
			set fchkSpin 0
			set spin $fchkSpin
			set orb $fchkOrb
			set OrbEList ""
			.fchkw.gridstuff.mesh.nx delete 0 end
			.fchkw.gridstuff.mesh.ny delete 0 end
			.fchkw.gridstuff.mesh.nz delete 0 end
			.fchkw.gridstuff.limits.xmin delete 0 end
			.fchkw.gridstuff.limits.xmax delete 0 end
			.fchkw.gridstuff.limits.ymin delete 0 end
			.fchkw.gridstuff.limits.ymax delete 0 end
			.fchkw.gridstuff.limits.zmin delete 0 end
			.fchkw.gridstuff.limits.zmax delete 0 end
			.fchkw.filestuff.filename delete 0 end
		}
	}
	button .fchkw.gridstuff.mesh.grain -text "Reset Grain" -command {
		if {$fchkLoaded == 1} {
			setgrain [.fchkw.gridstuff.mesh.nx get] [.fchkw.gridstuff.mesh.ny get] [.fchkw.gridstuff.mesh.nz get]
			.fchkw.gridstuff.mesh.nx delete 0 end
			.fchkw.gridstuff.mesh.ny delete 0 end
			.fchkw.gridstuff.mesh.nz delete 0 end
			set tlist [getparms 0]
			.fchkw.gridstuff.mesh.nx insert 0 $tlist
			set tlist [getparms 1]
			.fchkw.gridstuff.mesh.ny insert 0 $tlist
			set tlist [getparms 2]
			.fchkw.gridstuff.mesh.nz insert 0 $tlist
		}
	}
	button .fchkw.gridstuff.limits.scalebox -text "Scale Box" -command {
		if {$fchkLoaded == 1} {
			set cv [.fchkw.gridstuff.limits.scalf get]
			puts $cv
			scalebox $cv
			.fchkw.gridstuff.limits.xmin delete 0 end
			.fchkw.gridstuff.limits.xmax delete 0 end
			.fchkw.gridstuff.limits.ymin delete 0 end
			.fchkw.gridstuff.limits.ymax delete 0 end
			.fchkw.gridstuff.limits.zmin delete 0 end
			.fchkw.gridstuff.limits.zmax delete 0 end
			set tlist [getparms 3]
			.fchkw.gridstuff.limits.xmin insert 0 [lindex $tlist 0]
			.fchkw.gridstuff.limits.ymin insert 0 [lindex $tlist 1]
			.fchkw.gridstuff.limits.zmin insert 0 [lindex $tlist 2]
			set tlist [getparms 4]
			.fchkw.gridstuff.limits.xmax insert 0 [lindex $tlist 0]
			set tlist [getparms 5]
			.fchkw.gridstuff.limits.ymax insert 0 [lindex $tlist 1]
			set tlist [getparms 6]
			.fchkw.gridstuff.limits.zmax insert 0 [lindex $tlist 2]
		}
	}
	button .fchkw.dismiss -text Dismiss -command {set fchkSpin $spin; set fchkOrb $orb; set fchkFilename [.fchkw.filestuff.filename get]; destroy .fchkw}
	pack .fchkw.filestuff
	grid .fchkw.filestuff.filenamelabel -row 1 -column 1 -in .fchkw.filestuff
	grid .fchkw.filestuff.filename -row 1 -column 2 -in .fchkw.filestuff
	grid .fchkw.filestuff.getfile -row 1 -column 3 -in .fchkw.filestuff
	grid .fchkw.filestuff.openfile -row 2 -column 2 -in .fchkw.filestuff
	grid .fchkw.filestuff.freemem -row 2 -column 3 -in .fchkw.filestuff
	pack .fchkw.gridstuff 
	grid .fchkw.gridstuff.grid -row 1 -column 1 -in .fchkw.gridstuff 
	grid .fchkw.gridstuff.mesh -row 2 -column 1 -in .fchkw.gridstuff
	grid .fchkw.gridstuff.mesh.nxlab -row 1 -column 1 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.nx -row 1 -column 2 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.nylab -row 2 -column 1 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.ny -row 2 -column 2 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.nzlab -row 3 -column 1 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.nz -row 3 -column 2 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.mesh.grain -row 2 -column 4 -in .fchkw.gridstuff.mesh
	grid .fchkw.gridstuff.limits.xsize -row 1 -column 1 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.xmin -row 1 -column 2 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.xmax -row 1 -column 3 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.ysize -row 2 -column 1 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.ymin -row 2 -column 2 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.ymax -row 2 -column 3 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.zsize -row 3 -column 1 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.zmin -row 3 -column 2 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.zmax -row 3 -column 3 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.scallab -row 4 -column 1 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.scalf -row 4 -column 2 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits.scalebox -row 4 -column 3 -in .fchkw.gridstuff.limits
	grid .fchkw.gridstuff.limits -row 2 -column 2 -in .fchkw.gridstuff
	pack .fchkw.gridoptions
	grid .fchkw.gridoptions.denlab -row 1 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.dens -row 1 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.lplclab -row 2 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.lplc -row 2 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spinlab -row 3 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.pspin -row 3 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spindlab -row 4 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spind -row 4 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.orblab -row 5 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.orb -row 5 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.mesplab -row 6 -column 1 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.mesp -row 6 -column 2 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spin -row 1 -column 3 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spina -row 2 -column 3 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.alpha -row 2 -column 4 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.spinb -row 3 -column 3 -in .fchkw.gridoptions
	grid .fchkw.gridoptions.beta -row 3 -column 4 -in .fchkw.gridoptions
	pack .fchkw.orbcmd
	grid .fchkw.orbcmd.orbitals.energylab -row 1 -column 1 -in .fchkw.orbcmd.orbitals
	grid .fchkw.orbcmd.orbitals.orblist -row 2 -column 1 -in .fchkw.orbcmd.orbitals
	grid .fchkw.orbcmd.orbitals.yscroll  -row 2 -column 2 -in .fchkw.orbcmd.orbitals
	grid .fchkw.orbcmd.commands.calcgrid -row 1 -column 1 -in .fchkw.orbcmd.commands
	grid .fchkw.orbcmd.commands.savegrid -row 2 -column 1 -in .fchkw.orbcmd.commands
	grid .fchkw.orbcmd.commands -row 1 -column 1 -in .fchkw.orbcmd
	grid .fchkw.orbcmd.orbitals -row 1 -column 2 -in .fchkw.orbcmd
	pack .fchkw.dismiss
	set spin $fchkSpin
	set orb $fchkOrb
	if {$spin == 0} then {
		.fchkw.gridoptions.alpha select
	} else {
		.fchkw.gridoptions.beta select
	}
	if {$orb == 1} then {
		.fchkw.gridoptions.orb select
	} elseif {$orb == 2} then {
		.fchkw.gridoptions.dens select
	} elseif {$orb == 3} then {
		.fchkw.gridoptions.lplc select
	} elseif {$orb == 4} then {
		.fchkw.gridoptions.pspin select
	} elseif {$orb == 5} then {
		.fchkw.gridoptions.spind select
	} elseif {$orb == 6} then {
		.fchkw.gridoptions.mesp select
	}
	if {[string compare $fchkFilename ""]!=0} then {
		.fchkw.filestuff.filename insert 0 $fchkFilename
	}
	if {$fchkLoaded==1} then  {
		set tlist [getparms 0]
		.fchkw.gridstuff.mesh.nx insert 0 $tlist
		set tlist [getparms 1]
		.fchkw.gridstuff.mesh.ny insert 0 $tlist
		set tlist [getparms 2]
		.fchkw.gridstuff.mesh.nz insert 0 $tlist
		set tlist [getparms 3]
		.fchkw.gridstuff.limits.xmin insert 0 [lindex $tlist 0]
		.fchkw.gridstuff.limits.ymin insert 0 [lindex $tlist 1]
		.fchkw.gridstuff.limits.zmin insert 0 [lindex $tlist 2]
		set tlist [getparms 4]
		.fchkw.gridstuff.limits.xmax insert 0 [lindex $tlist 0]
		set tlist [getparms 5]
		.fchkw.gridstuff.limits.ymax insert 0 [lindex $tlist 1]
		set tlist [getparms 6]
		.fchkw.gridstuff.limits.zmax insert 0 [lindex $tlist 2]
		set OrbEList [listfchkenergies $spin]
	}
}

proc InitFCHK {} {
	lulLoadPlugin "fchk"
	global fchkLoaded
	global fchkGridCalc
	global fchkSpin
	global fchkOrb
	global fchkFilename
	set fchkLoaded 0
	set fchkGridCalc 0
	set fchkSpin 0
	set fchkOrb 2
	set fchkFilename ""
    set ::gomPluginsMenuItems(fchk) [list "FCHK Handling" {FCHKPlugin::fchk}]
    ::lulUpdatePluginsMenu
	fchk
}
}

if {[show graphics]} {
    set ::gomPluginsMenuItems(fchk) [list "FCHK Handling" {FCHKPlugin::InitFCHK}]
    ::lulUpdatePluginsMenu
}

