#######################################################################
#
# Make atom tree
#
# Leif Laaksonen 2000
#
# Based on code from the BWidget demos
#
# Modified and improved by
#
# Eero Häkkinen 2002/2003
#

package require BWidget

namespace eval lulAtomTree {

    variable dblclick
    variable numstruct
    variable tname
    variable lname

}

########################################################################
# PROC
proc lulAtomTree::create {} {
    global gomHelpFile
    global gomControlFont

    variable numstruct
    variable tname
    variable lname

    set  numstruct [show molstruct]

    if {$numstruct < 1} return

    set nb    .atomtree
    catch {destroy $nb}
    toplevel $nb 
    wm title    $nb "Atom tree"
    wm iconname $nb "Atom tree"

    set  frame [frame $nb.flat -borderwidth 2 -relief raised]
    set pw    [PanedWindow $frame.pw -side top -weights extra]
    set pane  [$pw add -weight 0]
    set title [TitleFrame $pane.lf -text "Atom tree"]
    set sw    [ScrolledWindow [$title getframe].sw \
                  -relief sunken -borderwidth 2]
    set tree  [Tree $sw.tree \
                   -relief flat -borderwidth 0 -highlightthickness 0 \
                   -height 25 -width 25 \
		   -redraw 0 -dropenabled 1 -dragenabled 1 \
                   -dragevent 3 \
                   -droptypes {
                       TREE_NODE    {copy {} move {} link {}}
                       LISTBOX_ITEM {copy {} move {} link {}}
                   } \
                   -opencmd   "lulAtomTree::modstructure 1 $sw.tree" \
                   -closecmd  "lulAtomTree::modstructure 0 $sw.tree"]
    $sw setwidget $tree
    set tname $tree
    pack $sw    -side top -fill both -expand 1

    # This has to be done before creating the next pane.
    # Otherwise PanedWindow won't calculate the desired height and will
    # create a very tiny window.
    lulAtomTree::initTree $tree
    ::update

    set pane [$pw add -weight 1]
    set lf   [TitleFrame $pane.lf -text "Content"]
    set sw   [ScrolledWindow [$lf getframe].sw \
                  -scrollbar horizontal -auto none -relief sunken -borderwidth 2]
    set list [ListBox::create $sw.lb -bg white \
                  -relief flat -borderwidth 0 \
                  -dragevent 3 \
                  -dropenabled 1 -dragenabled 1 \
                  -width 30 -highlightthickness 0 \
		  -multicolumn true \
                  -redraw 0 -dragenabled 1 \
                  -droptypes {
                      TREE_NODE    {copy {} move {} link {}}
                      LISTBOX_ITEM {copy {} move {} link {}}}]
    $sw setwidget $list
    set lname $list
    pack $sw $lf          -fill both -expand 1
    pack $pw              -fill both -expand 1
    pack $title -side top -fill both -expand 1
    pack $frame -side top -fill both -expand 1

    # The tree is created. So select the first item and update the list.
    lulAtomTree::selectFirstItem $tree $list

    $tree bindText  <ButtonPress-1>        "lulAtomTree::select tree 1 $tree $list"
    $tree bindText  <Double-ButtonPress-1> "lulAtomTree::select tree 2 $tree $list"
    $list bindText  <ButtonPress-1>        "lulAtomTree::select list 1 $tree $list"
    $list bindText  <Double-ButtonPress-1> "lulAtomTree::select list 2 $tree $list"

    frame  $nb.buttons -borderwidth 2 -relief raised
    pack   $nb.buttons -side bottom -fill x -pady 2m
    button $nb.buttons.update -text Update -font "$gomControlFont" \
           -command "lulAtomTree::update $tree $list"
    button $nb.buttons.dismiss -text Dismiss -font "$gomControlFont" \
           -command "destroy $nb"
    button $nb.buttons.apply   -text "Reset list"   -font "$gomControlFont" \
           -command "picking reset;lulInvalidateDisplay"
    button $nb.buttons.help    -text Help    -font "$gomControlFont" \
           -command \
           "htmlShowHelp $gomHelpFile(atomtree)"
    pack   $nb.buttons.update $nb.buttons.dismiss \
	$nb.buttons.apply $nb.buttons.help -side left -padx 5m -expand 1
}


########################################################################
## PROC
proc lulAtomTree::initTree { tree } {
    getAtoms $tree root
    $tree configure -redraw 1
}

########################################################################
## PROC
proc lulAtomTree::selectFirstItem { tree list } {
    set nodes [$tree nodes root]
    if {[llength $nodes] > 0} {
	lulAtomTree::select_node $tree $list [lindex $nodes 0]
    }
    $list configure -redraw 1
}
########################################################################
# PROC
proc lulAtomTree::addAtoms { tree node Wstr min max } {
    global gomControlFont

    for {set i $min} {$i <= $max} {incr i} {
	set SegName   [show atom segment $i $Wstr]
	set ResName   [show atom residue $i $Wstr]
	set AtmName   [show atom atom    $i $Wstr]
	set ResNum    [show atom resnum  $i $Wstr]
	set Mass      [show atom amass   $i $Wstr]
	set Hbond     [show atom hbond   $i $Wstr]
	set Connect   [show atom connect $i $Wstr]
	set NucCharge [show atom nuclear $i $Wstr]
	set APCharge  [show atom charge  $i $Wstr]
	scan [show atom colour  $i $Wstr] "%f %f %f" red green blue
	scan [show atom coord   $i $Wstr] "%f %f %f" xc yc zc
	set Atype     [show atom type $i $Wstr]
	set Vdw       [show atom vdw  $i $Wstr]

	set Colour [lulColourFloat2Hex $red $green $blue]

	set DataText "\{m${Wstr}a$i\} \{Type: $Atype\}
		      \{Atom coords x,y,z (A): $xc $yc $zc\}
		      \{Atom colour: $red $green $blue\}
		      \{Vdw radius (A): $Vdw\}
		      \{Mass: $Mass\}
		      \{Nuclear charge: $NucCharge\}
		      \{Atom partial charge: $APCharge\}
		      \{Hbond: $Hbond\}
		      \{Connectivity: $Connect\}"

	if {$Colour == "#ffffff"} {set Colour #000000}

	$tree insert end $node m${Wstr}a$i \
	    -text      "$SegName:$ResName\($ResNum\):$AtmName\($i\)" \
	    -fill      $Colour \
	    -image     [Bitmap::get file] \
	    -drawcross never \
	    -data      $DataText \
	    -font      $gomControlFont
    }
}

########################################################################
# PROC
proc lulAtomTree::addResidues { tree node Wstr min max force } {
    global gomControlFont

    set i $min
    while {$i <= $max} {
        set segName   [show atom segment   $i $Wstr]
        set resName   [show atom residue   $i $Wstr]
        set resNumber [show atom resnumber $i $Wstr]
	set first $i
	incr i
	while {$i <= $max && [show atom resnumber $i $Wstr] == $resNumber } {
	    incr i
	}
	set last [expr ($i-1)]
	if {$force == 0} {
	    if {$first > 1 || $last < $max} {
		# multiple residues
		return 1
	    }
	    set force 1
	}
	if {$force == 2 || ($first < $last)} {
	    $tree insert end $node m${Wstr}r$first-$last \
		-text      $segName:${resName}(${resNumber}) \
		-image     [Bitmap::get folder] \
		-drawcross allways \
		-font      $gomControlFont
	} else {
	    lulAtomTree::addAtoms $tree $node $Wstr $first $last
	}
    }
    return 0
}

########################################################################
# PROC
proc lulAtomTree::addSegments { tree node Wstr } {
    global gomControlFont

    set count [show numatoms $Wstr]
    set i 1
    while {$i <= $count} {
        set segName [show atom segment $i $Wstr]
	set first $i
	incr i
	while {$i <= $count && [show atom segment $i $Wstr] == $segName } {
	    incr i
	}
	set last [expr ($i-1)]
	if {$first == 1 && $last == $count} {
	    if {[lulAtomTree::addResidues $tree $node $Wstr 1 $count 0] == 0} {
		# only one residue in a segment
		return
	    }
	}
	$tree insert end $node m${Wstr}s$first-$last \
	    -text      $segName \
	    -image     [Bitmap::get folder] \
	    -drawcross allways \
	    -font      $gomControlFont
    }
}

########################################################################
# PROC
proc lulAtomTree::addStructures { tree node } {
    global gomStructureNames
    global gomControlFont

    set count [show molstruct]
    for {set i 1} {$i <= $count} {incr i} {
	$tree insert end $node m$i \
	    -text      $gomStructureNames($i) \
	    -image     [Bitmap::get openfold] \
	    -drawcross auto \
	    -font      $gomControlFont
	getAtoms $tree m$i
	$tree opentree m$i 0
    }
}

########################################################################
# PROC
proc lulAtomTree::getAtoms { tree node } {
    scan $node "m%d%*1s%d-%d" Wstr first last
    switch -glob $node {
    m*r* {lulAtomTree::addAtoms      $tree $node $Wstr $first $last}
    m*s* {lulAtomTree::addResidues   $tree $node $Wstr $first $last 1}
    m*   {lulAtomTree::addSegments   $tree $node $Wstr}
    root {lulAtomTree::addStructures $tree $node}
    }
    catch {$tree itemconfigure $node -drawcross auto}
}

########################################################################
# PROC
proc lulAtomTree::modstructure { open tree node } {
    if { $open } {
	if {[$tree itemcget $node -drawcross] == "allways" } {
	    getAtoms $tree $node
	}
	if { [llength [$tree nodes $node]] } {
	    $tree itemconfigure $node -image [Bitmap::get openfold]
	} else {
	    $tree itemconfigure $node -image [Bitmap::get file]
	}
    } else {
	$tree itemconfigure $node -image [Bitmap::get folder]
    }
}


########################################################################
# PROC
proc lulAtomTree::select { where num tree list node } {
    variable dblclick

    set dblclick 1
    switch $num {
    1 {
	if { $where == "tree" && [lsearch [$tree selection get] $node] != -1 } {
            unset dblclick
            after 500 "lulAtomTree::edit tree $tree $list $node"
            return
        }
        if { $where == "list" && [lsearch [$list selection get] $node] != -1 } {
            unset dblclick
            after 500 "lulAtomTree::edit list $tree $list $node"
            return
        }
        if { $where == "tree" } {
            select_node $tree $list $node
        } else {
            $list selection set $node
        }
    }
    2 {
	if { $where == "list" && [$tree exists $node]} {
	    if {![string match $node "m*:*"]} {
	        set parent [$tree parent $node]
		modstructure 1 $tree $parent
		$tree opentree $parent 0

		$tree see $node
		$tree selection set $node
		select_node $tree $list $node
	    }
	}
    } }
}


########################################################################
# PROC
proc lulAtomTree::select_node { tree list node } {
    global gomControlFont

    $tree selection set $node
    ::update
    eval {$list} delete [$list item 0 end]

    set dir [$tree itemcget $node -data]
    if { [$tree itemcget $node -drawcross] == "allways" } {
        getAtoms $tree $node
        set dir [$tree itemcget $node -data]
    }

    if {"$dir" == ""} {
	foreach subnode [$tree nodes $node] {
	    $list insert end $subnode \
		-text  [$tree itemcget $subnode -text] \
		-image [Bitmap::get folder]
	}
    } else {
	set temp [split [lindex $dir 0] :]

	scan $temp "m%da%d" Wstr i
	scan [show atom colour  $i $Wstr] "%f %f %f" red green blue
	set Colour [lulColourFloat2Hex $red $green $blue]

	picking atom $i $Wstr
	display
	::update
       
	if {$Colour == "#ffffff"} {set Colour #000000}
#	$list configure -bg [::lulGetVisibleForegroundColour $Colour]

	set loop 0
	foreach f $dir {
	    if {$loop > 0} {
		$list insert end m${Wstr}a$i:$loop \
		    -text  $f \
		    -font  $gomControlFont \
		    -fill  $Colour \
		    -image [Bitmap::get file]
	    }
	    incr loop
	}
    }
}

########################################################################
# PROC
proc lulAtomTree::pick_node { Wstr Atom } {
    variable tname
    variable lname

    set tree $tname
    set list $lname

    set parent root

    if {[$tree exists m$Wstr]} {set parent m$Wstr}

    set node $parent

    while {"$parent" == "$node"} {
	if {"$parent" != "root"} {
	    modstructure 1 $tree $parent
	    $tree opentree $parent 0
	}
	foreach node [$tree nodes $parent] {
	    if {[scan $node "m%d%*1s%d-%d" Wstr2 first last] < 3} {
		set last $first
	    }
	    if {$Wstr2 != $Wstr} return
	    if {$Atom >= $first && $Atom <= $last} {
		set parent $node
		break
	    }
	}

	if {"$node" == "m${Wstr}a$Atom"} {
	    $tree see $node
	    $tree selection set $node
	    select_node $tree $list $node
	    return
	}
    }
}

########################################################################
# PROC
proc lulAtomTree::edit { where tree list node } {
    global gomStructureNames
    variable dblclick

    if { [info exists dblclick] } {
        return
    }

    switch $where {
    list {set cntrl $list}
    tree {set cntrl $tree}
    }

    if { [lsearch [$cntrl selection get] $node] == -1 } return

    if { $where == "list" && [string match m*a*:* $node] } {

	scan $node "m%da%d" Wstr i
	set text  [$list itemcget $node -text]
	regsub {: .*$} $text {} title
	regsub {^.*: } $text {} text

	switch -glob $title {
	    "Type" {
		# Only capital letters and at least one letter.
		# Title can be removed.
		regsub {^.*: *} [$list edit $node "$title: $text" \
		    {regexp {^(.*:)? *[A-Za-z]{1,4} *$}}] \
		    {} newText
		set newText [string trim $newText]
		if {$newText == "" || $newText == $text} return
		set label [show atom atomname $i $Wstr]
		define atom label $newText $i $Wstr
		fill   atom                $i $Wstr
		define atom label $label   $i $Wstr
	    }
	    "Atom coords*" {
		# Three floats.
		# Only digits (at least one) and optionally one point.
		# Title can be removed.
		regsub {^.*: *} [$list edit $node "$title: $text" \
		    {regexp {^(.*:)? *-?([0-9]+\.|\.?[0-9])[0-9]*( +-?([0-9]+\.|\.?[0-9])[0-9]*){2,2} *$}}] \
		    {} newText
		set newText [string trim $newText]
		if {$newText == "" || $newText == $text} return
		foreach var {x y z} coord $newText trans [show translation array] {
		    set $var [expr $coord-$trans]
		}
		define atom coordinates $x $y $z $i $Wstr
	    }
	    "Atom colour*" {
		# Three float in range 0...1 or colour name.
		# Title can be removed.
		regsub {^.*: *} [$list edit $node "$title: $text" \
		    {regexp {^(.*:)? *(0?\.[0-9]+|[01]\.0*)( +(0?\.[0-9]+|[01]\.0*)){2,2}|[A-Za-z]+ *$}}] \
		    {} newText
		set newText [string trim $newText]
		if {$newText == "" || $newText == $text} return
		define atom colour $newText $i $Wstr
	    }
	    "Vdw radius*" {
		# Only digits (at least one) and optionally one point.
		# Title can be removed.
		regsub {^.*: *} [$list edit $node "$title: $text" \
		    {regexp {^(.*:)? *-?([0-9]+\.|\.?[0-9])[0-9]* *$}}] \
		    {} newText
		set newText [string trim $newText]
		if {$newText == "" || $newText == $text} return
		define atom vdw $newText $i $Wstr
	    }
	    "Atom partial charge*" {
		# Only digits (at least one) and optionally one point.
		# Title can be removed.
		regsub {^.*: *} [$list edit $node "$title: $text" \
		    {regexp {^(.*:)? *-?([0-9]+\.|\.?[0-9])[0-9]* *$}}] \
		    {} newText
		set newText [string trim $newText]
		if {$newText == "" || $newText == $text} return
		define atom charge $newText $i $Wstr
	    }
	    default {return}
	}

	lulInvalidateDisplay

	set parent [$tree parent m${Wstr}a$i]

	$tree closetree $parent
	$tree selection clear
	eval {$tree} delete [$tree nodes $parent]
	getAtoms $tree $parent
	modstructure 1 $tree $parent
	$tree opentree $parent 0

	$tree see m${Wstr}a$i
	$tree selection set m${Wstr}a$i
	select_node $tree $list m${Wstr}a$i
	$list selection set $node

	return
    }

    set text [$cntrl itemcget $node -text]
    switch -glob $node {
	m*a* {set res [$cntrl edit $node $text {regexp {^[^:()]+:[^:()]+\([0-9]+\):[^:()]+(\([0-9]+\))?$}}]}
	m*r* {set res [$cntrl edit $node $text {regexp {^[^:()]+:[^:()]+\([0-9]+\)$}}]}
	m*s* {set res [$cntrl edit $node $text]}
	m*   {set res [$cntrl edit $node $text]}
    }
    if {$res == ""} return
    if {[regexp {^m[0-9]+$} $node]} {
	scan $node "m%d" Wstr
	set gomStructureNames($Wstr) $res
	lulPushName2NameStack ""
	for {set i 1} {$i <= [show molstructures]} {incr i} {
	    lulPushName2NameStack "$gomStructureNames($i)"
	}
        $cntrl itemconfigure $node -text $res
    } else {
	regsub -all {[:()]} "$text" " " text
	regsub -all {[:()]} "$res"  " " res
	set count [scan "$text" "%s %s %d %s" oldSeg oldRes oldResN oldAtm]
	           scan "$res"  "%s %s %d %s" newSeg newRes newResN newAtm
	scan $node "m%d%1s%d-%d" Wstr type first last
	if {$type == "a"} {set last $first}
	for {set i $first} {$i <= $last} {incr i} {
	    if {$count >= 1 && $newSeg  != $oldSeg } {
		define segment label  $newSeg  $i $Wstr}
	    if {$count >= 2 && $newRes  != $oldRes } {
		define residue label  $newRes  $i $Wstr}
	    if {$count >= 3 && $newResN != $oldResN} {
		define atom resnumber $newResN $i $Wstr}
	    if {$count >= 4 && $newAtm  != $oldAtm } {
		define atom    label  $newAtm  $i $Wstr
		fill   atom                    $i $Wstr}
	}

	set parent [$tree parent $node]

	$tree closetree $parent
	$tree selection clear
	eval {$tree} delete [$tree nodes $parent]
	getAtoms $tree $parent
	modstructure 1 $tree $parent
	$tree opentree $parent 0

	if {![$tree exists $node]} {set node $parent}
	$tree see $node
	$tree selection set $node
	select_node $tree $list $node
    }
}

########################################################################
# PROC
proc lulAtomTree::update { tree list } {
    set treesel [$tree selection get]

    if {[llength $treesel] > 0} {
	set node [lindex $treesel 0]
	set treesel [list $node]
	set parent [$tree parent $node]
	while {"$parent" != "root"} {
	    set treesel [linsert $treesel 0 $parent]
	    set node $parent
	    set parent [$tree parent $node]
	}
    }

    $tree selection clear
    eval {$tree} delete [$tree nodes root]

    lulAtomTree::initTree $tree
    lulAtomTree::selectFirstItem $tree $list

    if {[llength $treesel] > 0} {
	scan [lindex $treesel 0] "m%d" Wstr
	set node ""
	foreach i [concat m$Wstr $treesel] {
	    if {[$tree exists $i]} {
		if {$node != ""} {$tree opentree $node 0}
		set node $i
		modstructure 1 $tree $node
	    }
	}
	if {$node != ""} {
	    $tree see $node
	    $tree selection set $node
	    select_node $tree $list $node
	}
    }
}
