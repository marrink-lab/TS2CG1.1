  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #   --- DISCLAIMER (by Melanie König m.konig@rug.nl):  	              	  	#
  #                     	                                                      	#
  #   This script can read triangulated surface (TS) files from DTS 	      		#
  #   like .tsi and .q files as well as 3D mesh files.	 				#
  #                    							               	#
  #                                                     	                     	#
  #   As always, you can modify, redistribute and make everything you 			#
  #   want with these few lines of code; if you write major improvement, 		#
  #   please let me know/test it!                                                  	#
  #                                                                          	 	#
  #											#
  #   TCL Script to visualize and modify triangulated surface    	          # # # #
  #   files like .tsi, .q from DTS simualtions                      	          #    #
  #   or 3D mesh files like .STL, .PLY    			          	  #   #
  #        								          # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### ----------------------------------------------------------------------------------------------------------------------------- USAGE
proc ts2vmd_usage {} {
	puts ""
	puts " USAGE"
  	puts "-------"
  	puts ""
  	puts "These few lines are given by the \"ts2vmd_usage\" command."
  	puts ""
  	puts "Loading TS or 3D mesh files:"
  	puts ""
  	puts " 		load_ts \[OPTIONS\]"
  	puts ""
  	puts "Writing TS or 3D mesh files:"
  	puts ""
  	puts "		write_ts \[OPTIONS\]"
  	puts ""
  	puts "Modifying vertices by:"
  	puts ""
  	puts "	adding domains"
  	puts " 		add_domain \[OPTIONS\]"
  	puts "	adding inclusions"
  	puts " 		add_inc \[OPTIONS\]"
  	puts "	adding exclusions"
  	puts " 		add_exc \[OPTIONS\]"
  	puts "	removing domains/inclusions/exclusions"
  	puts " 		clean_ts \[OPTIONS\]"
  	puts ""
  	puts "Show/hide periodic box:"
  	puts ""
  	puts "		show_box / hide_box"
  	puts ""
  	puts "Manully resize periodic box:"
  	puts ""
  	puts "		resize_box \[OPTIONS\]"
  	puts ""
  	puts "Show surface:"
  	puts ""
  	puts "		show_surf \[OPTIONS\]"
  	puts ""
}

ts2vmd_usage

#change basic vmd settings
display shadows on
display ambientocclusion on

### ---------------------------------------------------------------------------------------------------------------------------- UTILS

# check if file exists
proc file_exists { file } { if { [file exists $file] } { return -code 0 } else { return -code 1 "\nError: file $file does not exist.\n" } }

# check the if the number of entries matches the number from the keyword lines
proc check_entries {} {
	if {$::nver == [llength $::vertices]} { return -code 0
	} else { return -code 1 "\nError: Number vertices doesn't match.\n" }
	if {$::ntri == [llength $::triangles]} { return -code 0
	} else { return -code 1 "\nError: Number triangles doesn't match.\n" }
	if {$::ninc == [llength $::inclusions]} { return -code 0
	} else { return -code 1 "\nError: Number inclusions doesn't match.\n" }
	if {$::nexc == [llength $::exclusions]} { return -code 0
	} else { return -code 1 "\nError: Number exclusions doesn't match.\n" }
}

# process tsi file
proc load_ts { args } {

	set args [join $args]
	set args [split $args]

	# default values
	set mesh "off"
	set color 15
	set pbc "on"
	set scaling 1
	set convert_in 10

	# Arg parsing
	if { [llength $args] < 2 } {
		load_ts_usage
		return -code 1
	}
	foreach { n m } $args {
		if { $n == "-in" } {
			set file_in $m
			file_exists $file_in
			set format [ lindex [split $file_in "."] 1 ]
		}
		if { $n == "-mesh" } { set mesh $m }
		if { $n == "-color" } { set color $m }
		if { $n == "-box" } { set pbc $m }
                if { $n == "-unit" } {
                        switch $m {
                                "A" { set convert_in 1 }
                                "nm" { set convert_in 10 }
                        }
                }
		if { $n == "-rescale" } { set scaling $m }
	}

	# global variables defining the TS file
	set ::vertices {}
	set ::triangles {}
	set ::inclusions {}
	set ::exclusions {}
	set ::domains {}
	set ::ntri 0
	set ::nver 0
	set ::ninc 0
	set ::nexc 0

	switch $format {
		tsi {
			load_tsi $file_in $scaling $convert_in
			check_entries
			get_domains
		}
		q {
			load_q $file_in $scaling $convert_in
			check_entries
			get_domains
		}
		ply { load_ply $file_in $scaling $convert_in }
		xyz { load_xyz $file_in $scaling $convert_in }
		stl { load_stl $file_in $scaling $convert_in }
	}

	gen_coords $file_in

	# if inclusions and/or exclusions are present in the input file
	# they will be colored accordingly
	if { [llength $::domains] > 0 }  { show_dom }
	if { [llength $::inclusions] > 0 } { show_inc }
	if { [llength $::exclusions] > 0} { show_exc }

	# show box if wanted
	if { $pbc == "on" } { show_box }

	# show surface mesh if wanted
	if { $mesh != "off" } { show_surf [ format "-mesh %s -color" $mesh $color]}

	puts [format "\nYour triangulated surface %s consits of \n%s vertices connected by %s triangles\nwith %s inclusions and %s exclusions \n" $file_in $::nver $::ntri $::ninc $::nexc]
}

proc load_tsi { file_in scaling convert_in } {

	set tsi [open $file_in r]

	# read the different sections of a .tsi file
	while {[gets $tsi line] > 0} {
		if {[string match box* $line]} {
			set ::box "[expr [lindex $line 1] * $convert_in] [expr [lindex $line 2] * $convert_in] [expr [lindex $line 3] * $convert_in]"
      		} elseif {[string match vertex* $line]} {
			set ::nver "[lindex $line 1]"
			while {[gets $tsi line] >= 0 } {
        			if {[string match triangle* $line]} {
					set ::ntri "[lindex $line 1]"
					while {[gets $tsi line] >= 0} {
						if {[string match inclusion* $line]} {
							set ::ninc "[lindex $line 1]"
							while {[gets $tsi line] >= 0} {
	                                if {[string match exclusion* $line]} {
									set ::nexc "[lindex $line 1]"
									while {[gets $tsi line] >= 0} {
										if {[llength $line] > 0} {lappend ::exclusions $line }
									}
								}
								if {[llength $line] > 0} {lappend ::inclusions $line }
							}
						}
						if {[llength $line] > 0} {lappend ::triangles $line }
					}
				}
				if {[llength $line] > 0} {
					set x [expr [lindex $line 1] * $scaling * $convert_in ]
                        		set y [expr [lindex $line 2] * $scaling * $convert_in ]
                        		set z [expr [lindex $line 3] * $scaling * $convert_in ]
                        		lappend ::vertices [list [lindex $line 0] $x $y $z [lindex $line 4]]
				}
    			}
		}
	}
}


proc load_q { file_in scaling convert_in } {

	set q [open $file_in r]

	set read_ver "false"
	set read_tri "false"
	set i 0

	# read the different sections of a .tsi file
	while {[gets $q line] > 0} {
		if { $i == 0 } { 
			set ::box "[expr [lindex $line 0] * $convert_in] [expr [lindex $line 1] * $convert_in] [expr [lindex $line 2] * $convert_in]"
		} elseif { $i == 1 } {
			set ::nver $line
			set read_ver "true"
			set index 0
		} elseif { $read_ver == "true"} {
			set x [expr [lindex $line 1] * $scaling * $convert_in ]
			set y [expr [lindex $line 2] * $scaling * $convert_in ]
			set z [expr [lindex $line 3] * $scaling * $convert_in ]
			lappend ::vertices [list $index $x $y $z [lindex $line 4]]
			incr index
			if { $index == $::nver } {
				set read_ver "false"
				set read_tri "true"
				set index -1
			}
		} elseif { $read_tri == "true"} {
			if { $index == -1 } { set ::ntri $line
			} else { lappend ::triangles [join [list $index [lrange $line 1 end-1]]] }
			incr index
			if { $index == $::ntri } { set read_tri "false" }
		}
		incr i
	}

}

proc load_xyz { file_in } {

	mol new $file_in
	mol delrep 0 top
	mol representation vdw 1.200000 100.000000
	mol addrep 22
}

proc load_stl { file_in scaling } {

	set stl [open $file_in r]

	set ::facet_normal {}
	set ::box { 100 100 100 }

	# read the different sections of an ascii .stl file
	while {[gets $stl line] > 0} {
		if {[string match solid* $line]} {
			if { [lindex $line 1] != ""} { set name "[lindex $line 1]" }
		} elseif {[string match facet* $line]} {
			lappend ::facet_normal [lrange $line 1 end]
			set triangle {$::ntri}
		}  elseif {[string match vertex* $line]} {
			 if { [lsearch $::vertices [lrange $line 1 end]] == -1 } {
			 	set index [llength $::vertices]
			 	#lappend ::vertices [join [list $index [lrange $line 1 end]]]
                        	set x [expr [lindex $line 1] * $scaling]
                        	set y [expr [lindex $line 2] * $scaling]
                        	set z [expr [lindex $line 3] * $scaling]
                        	lappend ::vertices [list $index $x $y $z]
			 	incr ::nver
				lappend triangle $index
			 }  else { lappend triangle [lsearch $::vertices [lrange $line 1 end]] }
		} elseif {[string match endface* $line]} {
			incr ::ntri
			lappend ::triangles $triangle
		}
	}
}

proc load_ply { file_in scaling } {

	set ply [open $file_in r]

	set ::box { 100 100 100 }
	set read_ver "false"
	set read_tri "false"

	# read the different sections of an ascii .stl file
	while {[gets $ply line] > 0} {
		if {[string match {element vertex*} $line]} {
			set ::nver "[lindex $line 2]"
			#puts $::nver
		} elseif {[string match {element face*} $line]} {
			set ::ntri "[lindex $line 2]"
			#puts $::ntri
		} elseif {[string match end_header $line]} {
			set read_ver "true"
			set index 0
		} elseif { $read_ver == "true"} {
			set x [expr [lindex $line 0] * $scaling]
			set y [expr [lindex $line 1] * $scaling]
			set z [expr [lindex $line 2] * $scaling]
			lappend ::vertices [list $index $x $y $z]
			incr index
			if { $index == $::nver } {
				set read_ver "false"
				set read_tri "true"
				set index 0
			}
		} elseif { $read_tri == "true"} {
			lappend ::triangles [join [list $index [lrange $line 1 end]]]
			incr index
			if { $index == $::ntri } { set read_tri "false" }
		}
	}
}

proc get_domains {} {

	set ::domains {}

	foreach vertex $::vertices {
			if { [lindex $vertex 4] != ""} { lappend ::domains [list [lindex $vertex 0] [lindex $vertex 4]] }
	}
}

proc load_ts_usage {} {
	puts ""
	puts " USAGE"
  	puts "-------"
  	puts ""
  	puts "Loading TS or 3D mesh files:"
  	puts ""
  	puts "   load_ts \[OPTIONS\]"
  	puts ""
  	puts "Options and default values:"
  	puts "   -in        TS.tsi	      input file (*.tsi | *.q | *.stl | *.ply) "
  	puts "   -mesh      \"off\"           draw edges AND/OR faces (off | edge | face | both)"
  	puts "   -color     \"red\"           color (color name or VMD-defined ID) for mesh"
  	puts "   -box       \"on\"            show periodic box ( on | off )"
	puts "   -unit      \"nm\"            unit of coordinates ( A | nm )"
        puts "   -rescale   \"1\"             rescale input coordinates"
  	puts ""
}


#---------------------------------------------------------------------------------------------------------------------------------------- 
# VMD representation
# max. 4 different representations will be generated
# 1) All vertices
# 2) Domains
# 3) Inclusions
# 4) Exclusions
#----------------------------------------------------------------------------------------------------------------------------------------  

proc gen_coords { name } {

	# genearte dummy coordinated {0 0 0}
	set new [mol new atoms $::nver]
	animate dup $new
	set sel [atomselect $new all]

	set vertex 0
  	set newcoords {}

	# modify the dummy coordinates
  	foreach coord [$sel get {x y z}] {
		set offset "[lindex [lindex $::vertices $vertex] 1] [lindex [lindex $::vertices $vertex] 2] [lindex [lindex $::vertices $vertex] 3]"
		lappend newcoords [vecadd $coord $offset]
		incr vertex
  	}

  	$sel set {x y z} $newcoords
	set molid [molinfo top]
	mol rename $molid $name
	mol delrep $molid top
	mol representation VDW 1.200000 100.000000
	mol addrep $molid
	display resetview
	mol modcolor 0 $molid ColorID 9

}

proc show_dom {} {

	set molid [molinfo top]
	set repid [molinfo $molid get numreps]

	# add domains using different names and IDs
	foreach domain $::domains {
		set sel [atomselect top [format "index %s" [lindex $domain 0]]]
		#if {[lindex $domain 1] != 0} {$sel set name Domain}
		$sel set resname Domain
		$sel set type [lindex $domain 1]
	}

	mol addrep $molid
	#display resetview
	mol representation VDW 1.201000 100.000000
	mol modcolor $repid $molid Type
	mol modselect $repid $molid resname Domain
}

proc show_inc {} {

	set molid [molinfo top]
	set repid [molinfo $molid get numreps]

	# add inclusions using different names
	foreach inc $::inclusions {
		set sel [atomselect top [format "index %s" [lindex $inc 2]]]
		$sel set name Inc
		$sel set resid [lindex $inc 1]
	}

	mol addrep $molid
	#display resetview
	mol representation VDW 1.205000 100.000000
	mol modcolor $repid $molid Resid
	mol modselect $repid $molid name Inc

}

proc show_exc {} {

	set molid [molinfo top]
	set repid [molinfo $molid get numreps]

	# add exclusions using different names
	foreach exc $::exclusions {
		set sel [atomselect top [format "index %s" [lindex $exc 1]]]
		$sel set name Exc
	}

	mol addrep $molid
	#display resetview
	mol representation VDW 1.205000 100.000000
	mol modcolor $repid $molid ColorID 8
	mol modselect $repid $molid name Exc

}

#---------------------------------------------------------------------------------------------------------------------------------------- 
# Modify triangulated surface
# You can add domains, inclusions and exclusions or
# remove all modifications which will give you the blank vertices
#---------------------------------------------------------------------------------------------------------------------------------------- 

proc vertex_avail {vertex} {

	foreach inc $::inclusions exc $::exclusions {
		if {[lindex $inc 2] == $vertex || [lindex $exc 1] == $vertex} {return -code 1 "\nError: The vertex is already occupied by another inclusion or exclusion.\n"}
	}
}

proc add_domain { args } {

	set args [join $args]
	set args [split $args]

	# default values
	set domain_radius 0
	set overwrite "no"

	if {[ llength $args ] < 4 } {
		add_domain_usage
		return -code 1
	}
	foreach { n m } $args {
		if { $n == "-id" } { set id $m }
		if { $n == "-vertex" } { set vertex $m }
		if { $n == "-within" } { set domain_radius $m }
		if { $n == "-ignore" } { set overwrite $m }
	}

	# here we do not check if the vertex is available
	# a vertex can contain both an inclusion/exclusion and a domain

	set sel [atomselect top [format "within %s of index %s" $domain_radius $vertex]]
	$sel set resname Domain

	if { [llength $::domains]  == 0 } {
		foreach vertex [$sel get index] { lappend ::domains [list $vertex $id] }
		$sel set type $id
		show_dom
	} else {
		foreach vertex [$sel get index] {
			set sel2 [atomselect top [format "index %s"  $vertex]]
			set id_used [lsearch -index 0 $::domains $vertex]
			if { $id_used == -1 } {
				lappend ::domains [list $vertex $id]
				$sel2 set type $id
			} else { if { $overwrite == "yes" } {
				set ::domains [ lreplace $::domains $id_used $id_used [list $vertex $id] ]
				$sel2 set type $id
				}
			}
		}
	}
}

proc add_domain_usage {} {
	puts ""
	puts " USAGE"
  	puts "-------"
  	puts ""
  	puts "Adding vertex domains:"
  	puts ""
  	puts "   add_domain \[OPTIONS\]"
  	puts ""
  	puts "Options and default values:"
  	puts "   -id         1       domain id "
  	puts "   -vetex      0       vertex index"
  	puts "   -ignore    \"no\"     overwrite existing domains"
  	puts "   -within     10      included vertices within radius (in A)"
  	puts ""
}


proc add_inc { args } {

	set args [join $args]
	set args [split $args]

        if {[ llength $args ] < 4 } {
                add_inc_usage
                return -code 1
	}
	foreach { n m } $args {
		if { $n == "-type" } { set type $m }
		if { $n == "-vertex" } { set vertex $m }
	}
	# check if vertex is free
	vertex_avail $vertex

	set sel [atomselect top [format "index %s" $vertex]]
	$sel set name Inc
	$sel set resid [lindex $type]

	if { [llength $::inclusions]  == 0 } {
		lappend ::inclusions [format "%s  %s  %s     0    1" [expr [lindex $::inclusions end 0] + 1 ] $type $vertex]
		show_inc
	} else { lappend ::inclusions [format "%s  %s  %s     0    1" [expr [lindex $::inclusions end 0] + 1 ] $type $vertex] }

}

proc add_inc_usage {} {
        puts ""
        puts " USAGE"
        puts "-------"
        puts ""
        puts "Adding inclusions:"
        puts ""
        puts "   add_inc \[OPTIONS\]"
        puts ""
        puts "Options and default values:"
        puts "   -type     1      inclusion type"
        puts "   -vetex    0      vertex index"
        puts ""
}



proc add_exc { args } {

	set args [join $args]
	set args [split $args]

        if {[ llength $args ] < 4 } {
                add_exc_usage
                return -code 1
	}
	foreach { n m } $args {
		if { $n == "-radius" } { set radius $m }
		if { $n == "-vertex" } { set vertex $m }
	}
	# check if vertex is free
	vertex_avail $vertex

	set sel [atomselect top [format "index %s" $vertex]]
	$sel set name Exc

	if { [llength $::exclusions]  == 0 } {
		lappend ::exclusions [format "%s   %s   %s" [expr [lindex $::exclusions end 0] + 1 ] $vertex $radius]
			show_exc
	} else { lappend ::exclusions [format "%s   %s   %s" [expr [lindex $::exclusions end 0] + 1 ] $vertex $radius] }
}

proc add_exc_usage {} {
        puts ""
        puts " USAGE"
        puts "-------"
        puts ""
        puts "Adding exclusions:"
        puts ""
        puts "   add_exc \[OPTIONS\]"
        puts ""
        puts "Options and default values:"
        puts "   -radius      1      exclusion size"
        puts "   -vetex       0      vertex index"
        puts ""
}


#---------------------------------------------------------------------------------------------------------------------------------------- 
# clean vertices
#---------------------------------------------------------------------------------------------------------------------------------------- 

proc clean_ts { args } {

	set args [join $args]
	set args [split $args]
	set molid [molinfo top]

        if {[ llength $args ] < 4 } {
                clean_ts_usage
                return -code 1
        }

	foreach { n m } $args {
		if { $n == "-vertex" } { set vertex $m }
		if { $n == "-domain" } { set dom_id $m}
		if { $n == "-inclusion" } { set type_id $m}
		if { $n == "-exclusion" } { set exc_id $m }
		}

		if { [info exists vertex] } {
			if { $vertex == "all" } {
				set ::domains {}
				set ::inclusions {}
				set ::exclusions {}
				# remove reprensentations

				for {set i 1} {$i < [molinfo $molid get numreps] } {incr i}  {mol delrep 1 $molid}
				# update entries
				set sel [atomselect top "all"]
				$sel set name "X"
				$sel set resname "UNK"
				$sel set resid 0
				$sel set type "X"
			} else {
				# remove entries and update presentation
				clean_dom $vertex 0
				clean_inc $vertex 2
				clean_exc $vertex 1
			}
		}

		if { [info exists dom_id] } {
			if { $dom_id == "all" } {
				set ::domains {}
				# remove reprensentations
				for {set i 1} {$i < [molinfo $molid get numreps] } {incr i}  {
					if { [molinfo top get "{selection $i}"] == "{resname Domain}"} { mol delrep $i $molid }
				}
				# update entries
				set sel [atomselect top "resname Domain"]
				$sel set resname "UNK"
				$sel set type "X"
			} else { clean_dom $dom_id 1 }
		}

		if { [info exists type_id] } {
			if { $type_id == "all" } {
				set ::inclusions {}
				# remove reprensentations
				for {set i 1} {$i < [molinfo $molid get numreps] } {incr i}  {
					if { [molinfo top get "{selection $i}"] == "{name Inc}"} { mol delrep $i $molid }
				}
				# update entries
				set sel [atomselect top "name Inc"]
				$sel set name "X"
				$sel set resid 0
			} else { clean_inc $type_id 1 }
		}

		if { [info exists exc_id] } {
			if { $exc_id == "all" } {
				set ::exclusions {}
				# remove reprensentations
				for {set i 1} {$i < [molinfo $molid get numreps] } {incr i}  {
					if { [molinfo top get "{selection $i}"] == "{name Exc}"} { mol delrep $i $molid }
				}
				# update entries
				set sel [atomselect top "name Exc"]
				$sel set name "X"
			} else { clean_exc $exc_id 1 }
		}

}

proc clean_dom { val col } {

	set new_dom {}
	foreach domain $::domains { if { [lindex $domain $col ] != $val } {
		lappend new_dom $domain
		} else {
			set sel [atomselect top [format "index %s" [lindex $domain 0 ] ] ]
			$sel set type "X"
		}
	}
	set ::domains $new_dom
}

proc clean_inc { val col } {

	set new_inc {}
	foreach inc $::inclusions { if { [lindex $inc $col] != $val }  {
		lappend new_inc $inc
		} else {
			set sel [atomselect top [format "index %s" [lindex $inc 2 ] ] ]
			$sel set name "X"
			$sel set resid 0
		}
	}
	set ::inclusions $new_inc
}

proc clean_exc { val col } {

	set new_exc {}
	foreach exc $::exclusions { if { [lindex $exc $col]  != $val }  {
		lappend new_exc $exc
		} else {
			set sel [atomselect top [format "index %s" [lindex $exc 1 ] ] ]
			$sel set name "X"
		}
	}
	set ::exclusions $new_exc
}

proc clean_ts_usage {} {
        puts ""
        puts " USAGE"
        puts "-------"
        puts ""
        puts "Remove domains, inclusions and/or exclusions from TS file:"
        puts ""
        puts "   clean_ts \[OPTIONS\]"
        puts ""
        puts "Options and default values:"
        puts "   -domain      1      domain id (or all)"
        puts "   -vetex       0      vertex index (or all)"
	puts "   -inclusion   1      inclusion type (or all)"
	puts "   -exclusion   1      exclusion type (or all)"
        puts ""
}



#---------------------------------------------------------------------------------------------------------------------------------------- 
# Write to file
#---------------------------------------------------------------------------------------------------------------------------------------- 

proc write_ts { args } {

	set args [join $args]
	set args [split $args]

	#default unit is A
	set convert_out 0.1

        # Arg parsing
        if { [llength $args] < 2 } {
                write_ts_usage
                return -code 1
        }
	foreach { n m } $args {
		if { $n == "-out" } {
			set file_out $m
			set format [ lindex [split $m "."] 1 ]
		}
		if { $n == "-unit" } {
			switch $m {
				"A" { set convert_out 1 }
				"nm" { set convert-out 0.1 }
			}
		}

	}

	switch $format {
		tsi { write_tsi $file_out $convert_out }
		xyz { write_xyz $file_out }
	}

}

proc write_ts_usage {} {
        puts ""
        puts " USAGE"
        puts "-------"
        puts ""
        puts "Writing TS file:"
        puts ""
        puts "   write_ts \[OPTIONS\]"
        puts ""
        puts "Options and default values:"
        puts "   -out     out.tsi      output file ( *.tsi | *.xyz )"
	puts "   -unit    \"nm\"          coordinate units ( A | nm )"
        puts ""

}

proc write_tsi { file_out convert } {

	set tsi [open $file_out w]
	puts $tsi "version 1.1"
	set box "[expr [lindex $::box 0] * $convert] [expr [lindex $::box 1] * $convert] [expr [lindex $::box 2] * $convert]"
	puts $tsi [format "box	%s" $box ]
	puts $tsi [format "vertex	%s" [llength $::vertices]]
	foreach vertex $::vertices {
		set domain_id [lsearch -index 0 $::domains [lindex $vertex 0]]
		if { $domain_id == -1 } { set domain 0
    	} else { set domain [lindex $::domains $domain_id 1]}
		set x "[format "%.10f" [expr [lindex $vertex 1] * $convert]]"
  		set y "[format "%.10f" [expr [lindex $vertex 2] * $convert]]"
  		set z "[format "%.10f" [expr [lindex $vertex 3] * $convert]]"

    		puts $tsi [format "%s	%s	%s	%s	%s" [lindex $vertex 0] $x $y $z $domain ]
    	}
 	puts $tsi [format "triangle	%s" [llength $::triangles]]
 	foreach triangle $::triangles { puts $tsi [format "%s" $triangle]}
 	puts $tsi [format "inclusion	%s" [llength $::inclusions]]
 	foreach inclusion $::inclusions { puts $tsi [format "%s" $inclusion]}
 	puts $tsi [format "exclusion	%s" [llength $::exclusions]]
 	foreach exclusion $::exclusions { puts $tsi [format "%s" $exclusion]}

	close $tsi
}

# converts tsi vertices to xyz file format
proc write_xyz { file_out } {

	set xyz [open $file_out w]
	puts $xyz [format "%s\n" $::nver]

	set elements "C O N P S He Li Be B F Ne Na Mg Al Si Cl Ar"

	foreach vertex $::vertices {
    	if {[llength $vertex] != 0} {
			if {[lindex $vertex 4] == ""} { set domain 0
			} else {set domain [lindex $vertex 4]}
                        puts $xyz [format "%s  %.10f  %.10f  %.10f" [lindex $elements $domain] [expr [lindex $vertex 1] * $convert] [expr [lindex $vertex 2] * $convert] [expr [lindex $vertex 3] * $convert]]
                }
        }
        close $xyz
}

#---------------------------------------------------------------------------------------------------------------------------------------- 
# Draw pbc box
#---------------------------------------------------------------------------------------------------------------------------------------- 

proc show_box {} {
	pbc set $::box -all -molid top
	pbc box
}

proc hide_box {} { pbc box -off	}

proc resize_box { args } {

	set args [join $args]
	set args [split $args]

	if {[ llength $args ] < 6 } {
		resize_box_usage
		return -code 1
	}

	foreach { n m } $args {
		if { $n == "-x" } { set x $m }
		if { $n == "-y" } { set y $m }
		if { $n == "-z" } { set z $m }
	}

	# set new box
	remove
	set ::box [list $x $y $z]
	set center_box [list [expr $x/2] [expr $y/2] [expr $z/2]]

	# move the coordinates
	set sel [atomselect top all]
	set minus_com [vecsub $center_box [measure center $sel]]
	$sel moveby $minus_com
	show_box
	display resetview

	# uppate coordinates of vertices
	set vertex 0
  	foreach coord [$sel get {x y z}] {
  		set x "[format "%.10f" [lindex $coord 0]]"
  		set y "[format "%.10f" [lindex $coord 1]]"
  		set z "[format "%.10f" [lindex $coord 2]]"
		set ::vertices [ lreplace $::vertices $vertex $vertex "[format "%s	%s	%s	%s" $vertex $x $y $z]" ]
		incr vertex
  	}
}

proc resize_box_usage {} {
	puts ""
	puts " USAGE"
  	puts "-------"
  	puts ""
  	puts "Resizing the periodic box (only rectangular boxes):"
  	puts ""
  	puts "   resize_box \[OPTIONS\]"
  	puts ""
  	puts "Options and default values:"
  	puts "   -x	10	 size of box in x dimension"
  	puts "   -y	10	 size of box in y dimension"
  	puts "   -z	10	 size of box in z dimension"
  	puts ""
}


#---------------------------------------------------------------------------------------------------------------------------------------- 
# Draw surface and/or bonds of triangles
#---------------------------------------------------------------------------------------------------------------------------------------- 

proc tri_coords {vertices triangle} {
	set vert1 "[lindex $triangle 1]"
	set vert2 "[lindex $triangle 2]"
	set vert3 "[lindex $triangle 3]"
	set points {}
	if {[llength $triangle] != 0} {
		lappend points "[lindex [lindex $vertices $vert1] 1] [lindex [lindex $vertices $vert1] 2] [lindex [lindex $vertices $vert1] 3]"
		lappend points "[lindex [lindex $vertices $vert2] 1] [lindex [lindex $vertices $vert2] 2] [lindex [lindex $vertices $vert2] 3]"
		lappend points "[lindex [lindex $vertices $vert3] 1] [lindex [lindex $vertices $vert3] 2] [lindex [lindex $vertices $vert3] 3]"
	}
	return $points
}

proc show_surf { args } {

	set args [join $args]
	set args [split $args]

	# default values
	set color 15
	set material AOEdgy

        if {[ llength $args ] < 2 } {
                show_surf_usage
                return -code 1
        }

	foreach { n m } $args {
		if { $n == "-mesh" } { set mesh $m }
		if { $n == "-color" } { set color $m }
		if { $n == "-material" } { set material $m }
	}

	draw color $color
	draw material $material

	foreach triangle $::triangles {
		set points [tri_coords $::vertices $triangle]

		# show edges and or faces
		switch $mesh {
			edge { draw_lines $points }
			face { draw_triangles $points  }
			both {
				draw_lines $points
				draw_triangles $points
			}
		}
	}
}

# draw triangle faces between vertices
proc draw_triangles { points } {

	set cross_pbc "no"

	for {set i 0} {$i < 3} {incr i} {
		for {set j [expr $i + 1] } {$j < 3} {incr j} {
			for {set k 0} {$k < 3} {incr k} {
				set coord [expr abs([lindex [vecsub [lindex $points $i] [lindex $points $j]] $k])]
				set half_box [expr [lindex $::box $k]/2]
				if { $coord > $half_box } {
					set cross_pbc "yes"
					break
					}
			}
		}
	}
	if { $cross_pbc == "no" } { draw triangle  [lindex $points 0] [lindex $points 1]  [lindex $points 2] }
}

# draw triangles edges between vertices
proc draw_lines { points } {

	set box_x [lindex $::box 0]
	set box_y [lindex $::box 1]
	set box_z [lindex $::box 2]

	for {set i 0} {$i < 3} {incr i} {
		for {set j [expr $i + 1] } {$j < 3} {incr j} {
			set x [expr abs([lindex [vecsub [lindex $points $i] [lindex $points $j]] 0])]
			set y [expr abs([lindex [vecsub [lindex $points $i] [lindex $points $j]] 1])]
			set z [expr abs([lindex [vecsub [lindex $points $i] [lindex $points $j]] 2])]
			if { $x < [expr $box_x/2] && $y < [expr $box_y/2]  && $z < [expr $box_z/2] } {
				draw cylinder [lindex $points $i] [lindex $points $j] radius 0.1
				} else { continue }
		}
	}
}


proc show_surf_usage {} {
        puts ""
        puts " USAGE"
        puts "-------"
        puts ""
        puts "Show triangulated surface:"
        puts ""
        puts "   show_surf \[OPTIONS\]"
        puts ""
        puts "Options and default values:"
        puts "   -mesh       edge        draw edges AND/OR faces ( edge | face | both )"
        puts "   -color      \"15\"        edge/face color"
        puts "   -material   \"AOEdgy\"    edge/face material"
        puts ""
}

proc remove {} {draw delete all}

