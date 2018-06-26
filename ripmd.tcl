package provide ripmd 1.0

namespace eval RIPMD:: {
	namespace export ripmd
	# window handles
	variable w                                          ;# handle to main window
	variable resultFolder
	##########################################
	##
	## variables for general options section
	##
	##########################################
	
	#checkbuttons variables
	variable calpha 1zx
	variable hbonds 1
	variable salt 1
	variable disulphide 1
	variable capi 1
	variable pipi 1
	variable ArgArg 1
	variable vdw 0
	variable coulomb 0
	variable pearson 0
	
	#output path variable
	variable dirName
	
	#final graph format
	variable format "GML"
	
	#selection to use
	
	variable proteinSelection "protein"
	
	#############################################
	##
	##variables for interaction options section
	##
	#############################################
	
	#for C alpha
   	variable calphaDist "8"
   	
   	#for hbond
	variable hbondDist "3"	
	variable hbondAngle "120"
	
	#for salt
	variable saltDist "6"
 
 	#pi - pi
	variable pipiDist "6"

 	#Arg - Arg
	variable ArgArgDist "5"
		
	#coulomb
	variable simulatedPermittivity "1"
	variable reactionFieldPermittivity "82"
	variable reactionField 1
	variable temp "298.5"
	variable inverseDebye "0"
	variable coulombThreshold "12"
	variable coulombCovalentDist "3"

	
	#cartion-pi
   	variable capiDist "7"
	variable capiAngle1 "0"
   	variable capiAngle2 "60"
   	variable capiAngle3 "120"
   	variable capiAngle4 "180"	

	#disulphide
	variable disulphideDist "3"
	variable disulphideAng1 "60"
	variable disulphideAng2 "90"
 
 	#vdw
	variable covalentDist "3"
  	variable distVDW1 "-0.1"
	variable distVDW2 "3"
	
	##########################################
	##
	## variables for MD options section
	##
	##########################################
	
	#MD time
	variable percentaje "75"
	
	##########################################
	##
	## variables for PDB options section
	##
	##########################################	
	variable missing 0
	#ph
	variable pH "7.0"
	
	#################################
	##
	## for RIP-MD execution
	##
	#################################
	
	#files on VMD
	variable pathToFiles
	variable pdb "*"
	variable dcd "*"
	variable psf "*"
	variable command 
	variable outputFolder "*"
	variable proc "1"
	
	variable errors
	
	
	################################
	##
	## for VMD representation and
	## displaying
	##
	###############################
	variable displayCalpha 1
	variable displayHbonds 1
	variable displaySalt 1
	variable displayDisulphide 1
	variable displayCapi 1
	variable displayPipi 1
	variable displayArgArg 1
	variable displayVdw 0
	variable displayCoulomb 0
	
	#frame to use
	variable frame "0"
	variable frame_start "0"
	variable frame_end "-1"
	variable frame_separation "0"
	
	#to show pearson
	variable pearson1 "All"
	variable pearson2 "All"
	
	#general variable to draw bonds
	variable residIndexDict
	variable indexResidDict
	variable matrixCalpha
	variable matrixHbond
	variable matrixSalt
	variable matrixDisulphide
	variable matrixCationPi
	variable matrixPiPi
	variable matrixArgArg
	variable matrixCoulomb
	variable matrixVdw
	variable nodes ""
	variable selection ""
	
	variable toUse "segid"
	

	variable forceField "RIP_MD/dat/top_all22_prot.rtf"
	variable parameterFile "RIP_MD/dat/par_all22_prot.prm"	
	variable hist "All"
}


#
# function to draw the target residues of each interaction
#
proc RIPMD::drawNotSelectedResidues {indices} {
	if {[llength $indices] == 0} {
		#if we dont have target residues we will return
		return
	}
	set chain {}
	foreach num $indices {
		set aux1 "[dict get $RIPMD::indexResidDict $num]"
		set aux [split $aux1 ":"]
		#if chain is not in the chain variable we will append it
		if {[string match *[lindex $aux 1]* $chain] != 1} {
			append chain "[lindex $aux 1] "
		}
	}
###
    #now we will append to a list of list the different residues
    set resids {}
    set resids [lrepeat [llength $chain] {}]
    set Selection ""
    foreach index $indices {
		#setting a list of residues to append in order of the chain
		set aux [split [dict get $RIPMD::indexResidDict $index] ":"]
		#now we will look for the actual chain to append the residue
		set position [lsearch $chain [lindex $aux 1]]
		set auxList [lindex $resids $position]
		append auxList "[lindex $aux 2] "
		set resids [lreplace $resids $position $position $auxList]
	}
	#at this point the selection is configured
	set countChain 0
	foreach item $chain {
		
		if {$countChain != [expr [llength $chain] - 1]} {
			append Selection "($RIPMD::toUse $item and resid [lindex $resids $countChain] ) or "
		} else {
			append Selection "($RIPMD::toUse $item and resid [lindex $resids $countChain] )"
		}

		incr countChain
	}
	#creating new representations and setting protein style

	set numbersOfRep [molinfo top get numreps] 
	for { set i 0 } { $i < $numbersOfRep} {incr i} {
		mol delrep end top
	}
	mol addrep top
	mol modselect 0 top protein
	mol modstyle 0 top trace 0.2
	mol modmaterial 0 top Transparent
	mol modcolor 0 top colorid 8

	mol addrep top
	mol modselect 1 top alpha
	mol modstyle 1 top vdw 0.5
	mol modcolor 1 top colorid 1	
	
	mol addrep top
	mol modselect 2 top "$RIPMD::selection"
	mol modstyle 2 top Bonds 0.2	
	
	mol addrep top
	mol modselect 3 top "$Selection"
	mol modstyle 3 top Bonds 0.2
}

#
# function to draw bonds
#
proc RIPMD::drawBonds {} {
	set indices {}
	#first we will delete all draws
	foreach ID [graphics top list] {
		graphics top delete $ID
	}
	#and now we will look for interactions to draw
	if {$RIPMD::nodes != ""} {
		set indices {}
		#we are sure that  anode has been selected
		foreach index [$RIPMD::nodes curselection] {
	
			#######################################
			##
			##To display Calpha interactions
			##
			######################################
			if {$RIPMD::displayCalpha == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				#toUse is = to chain or segid
				set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and alpha"]
				set coord1 [$res1 get {x y z}]	
				foreach item [lindex $RIPMD::matrixCalpha $index] {
					if {$item != 0} {
						#we will set the second coordinate
						set aux [split [$RIPMD::nodes get $cont] ":"]
						set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and alpha"]
						set coord2 [$res2 get {x y z}]	
						
						#now we will draw the dashed line
						set a [lindex $coord1 0]
						set b [lindex $coord2 0]
						#to select yellow color for contacts
						graphics top color 0
						graphics top line $a $b width 5 style dashed
						
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							
							lappend indices $cont
						}
						
					}
					incr cont
				}
			}

			#######################################
			##
			##To display HBonds
			##
			######################################
			if {$RIPMD::displayHbonds == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				foreach item [lindex $RIPMD::matrixHbond $index] {
					if {$item != {}} {
						foreach list $item {
							set resids [split [lindex [lindex $item 0] 0] " "]
							set firstAtom [lindex [split [lindex $resids 0] ":"] 3]
							set secondAtom [lindex [split [lindex $resids 2] ":"] 3]
							set aux [split [lindex $resids 0] ":"]
							set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $firstAtom"]
							set coord1 [$res1 get {x y z}]
							set aux [split [lindex $resids 2] ":"]
							set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $secondAtom"]
							set coord2 [$res2 get {x y z}]							
							#now we will draw the dashed line
							set a [lindex $coord1 0]
							set b [lindex $coord2 0]
							#to select yellow color for contacts
							graphics top color 8
							graphics top line $a $b width 5 style dashed
						}
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							lappend indices $cont
						}
					}
					incr cont
				}
			}

			#######################################
			##
			##To display Salt Bridges
			##
			######################################
			if {$RIPMD::displaySalt == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				foreach item [lindex $RIPMD::matrixSalt $index] {
					if {$item != {}} {
						foreach list $item {
							set resids [split [lindex [lindex $item 0] 0] " "]
							set firstAtom [lindex [split [lindex $resids 0] ":"] 3]
							set secondAtom [lindex [split [lindex $resids 2] ":"] 3]
							set aux [split [lindex $resids 0] ":"]
							set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $firstAtom"]
							set coord1 [$res1 get {x y z}]
							set aux [split [lindex $resids 2] ":"]
							set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $secondAtom"]
							set coord2 [$res2 get {x y z}]							
							#now we will draw the dashed line
							set a [lindex $coord1 0]
							set b [lindex $coord2 0]
							#to select yellow color for contacts
							graphics top color 3
							graphics top line $a $b width 5 style dashed
						}
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							lappend indices $cont
						}
					}
					incr cont
				}
			}
						
			#######################################
			##
			##To display Disulphide interactions
			##
			######################################
			if {$RIPMD::displayDisulphide == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type SG"]
				set coord1 [$res1 get {x y z}]	
				foreach item [lindex $RIPMD::matrixDisulphide $index] {
					if {$item != 0} {
						#we will set the second coordinate
						set aux [split [$RIPMD::nodes get $cont] ":"]
						set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type SG"]
						set coord2 [$res2 get {x y z}]	
						
						#now we will draw the dashed line
						set a [lindex $coord1 0]
						set b [lindex $coord2 0]
						#to select blue color for contacts
						graphics top color 4
						graphics top line $a $b width 5 style dashed
						
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							lappend indices $cont
						}
						
					}
					incr cont
				}
			}		

			#######################################
			##
			##To display cation pi interaction
			##
			######################################
			if {$RIPMD::displayCapi == 1} {
				set cont 0
				set aux [split [$RIPMD::nodes get $index] ":"]
				foreach item [lindex $RIPMD::matrixCationPi $index] {
					if {$item != {}} {
						foreach list $item {
							set firstRes [lindex [split [lindex $list 0] " "] 0]
							set secondRes [lindex [split [lindex $list 0] " "] 3]
							#getting the coord of the cation
							set cation [lindex [split $firstRes ":"] 0]
							set res1 ""
							set chain [lindex [split $firstRes ":"] 1]
							set resid [lindex [split $firstRes ":"] 2]
							if {$cation == "LYS"} {
								set res1 [atomselect top "$RIPMD::toUse $chain and resid $resid and type NZ"]
							} elseif {$cation == "HIS"} {
								set res1 [atomselect top "$RIPMD::toUse $chain and resid $resid and type NE2"]						
							} elseif {$cation == "ARG"} {
								set res1 [atomselect top "$RIPMD::toUse $chain and resid $resid and type CZ"]								
							}
							set coord1 [$res1 get {x y z}]
							
							#now we will extract the coord for the pi system
							set pi [lindex [split $secondRes ":"] 0]
							set res2 ""
							set chain [lindex [split $secondRes ":"] 1]
							set resid [lindex [split $secondRes ":"] 2]	
							set coord2 ""						

							if  {$pi == "PHE" || $pi == "TYR"} {
								set cont2 0	
								set coord2 [veczero]
								set res2 [atomselect top "$RIPMD::toUse $chain and resid $resid and type CG CD1 CD2 CE1 CE2 CZ"]
								set coords [$res2 get {x y z}]	
								
								foreach coord $coords {
									set coord2 [vecadd $coord2 $coord]
									incr cont2
								}
								set coord2 [vecscale [expr 1.0 / $cont2] $coord2]
								
							} elseif  {$pi == "TRP"} {
								set cont2 0	
								set coord2 [veczero]
								set res2 [atomselect top "$RIPMD::toUse $chain and resid $resid and type CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"]
								set coords [$res2 get {x y z}]	
								
								foreach coord $coords {
									set coord2 [vecadd $coord2 $coord]
									incr cont2
								}
								set coord2 [vecscale [expr 1.0 / $cont2] $coord2]
								
							} elseif  {$pi == "HIS"} {
								set cont2 0	
								set coord2 [veczero]
								set res2 [atomselect top "chain $chain and resid $resid and type CG ND1 CE1 NE2 CD2"]
								set coords [$res2 get {x y z}]	
								
								foreach coord $coords {
									set coord2 [vecadd $coord2 $coord]
									incr cont2
								}
								set coord2 [vecscale [expr 1.0 / $cont2] $coord2]
							}
							#to select pink color for contacts
							graphics top color 7
							graphics top line [lindex $coord1 0] $coord2 width 5 style dashed
						}
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							lappend indices $cont
						}
					}
					incr cont
				}
			}
						
			###################################
			##
			##To display Pi Pi interactions
			##
			###################################			
			if {$RIPMD::displayPipi == 1} {
				set cont 0
				set coord1 [veczero]
				set coord2 [veczero]
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				if  {[lindex $aux 0] == "PHE" || [lindex $aux 0] == "TYR"} {
					set cont 0	
					set coord1 [veczero]
					set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG CD1 CD2 CE1 CE2 CZ"]
					set coords [$res1 get {x y z}]	
					
					foreach coord $coords {
						set coord1 [vecadd $coord1 $coord]
						incr cont
					}
					set coord1 [vecscale [expr 1.0 / $cont] $coord1]
					
				} elseif  {[lindex $aux 0] == "TRP"} {
					set cont 0	
					set coord1 [veczero]
					set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"]
					set coords [$res1 get {x y z}]	
					
					foreach coord $coords {
						set coord1 [vecadd $coord1 $coord]
						incr cont
					}
					set coord1 [vecscale [expr 1.0 / $cont] $coord1]
					
				} elseif  {[lindex $aux 0] == "HIS" || [lindex $aux 0] == "HSD" || [lindex $aux 0] == "HSE" || [lindex $aux 0] == "HSP"} {
					set cont 0	
					set coord1 [veczero]
					set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG ND1 CE1 NE2 CD2"]
					set coords [$res1 get {x y z}]	
					foreach coord $coords {
						set coord1 [vecadd $coord1 $coord]
						incr cont
					}
					set coord1 [vecscale [expr 1.0 / $cont] $coord1]
				}		
				set cont2 0		
				
				foreach item [lindex $RIPMD::matrixPiPi $index] {
					if {$item != 0} {
						set coord2 [veczero]
						set aux [split [$RIPMD::nodes get $cont2] ":"]
						if  {[lindex $aux 0] == "PHE" || [lindex $aux 0] == "TYR"} {
	
							set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG CD1 CD2 CE1 CE2 CZ"]
	
							
						} elseif  {[lindex $aux 0] == "TRP"} {
							set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"]
							
						} elseif  {[lindex $aux 0] == "HIS" || [lindex $aux 0] == "HSD" || [lindex $aux 0] == "HSE" || [lindex $aux 0] == "HSP"} {
							set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CG ND1 CE1 NE2 CD2"]				
						}
						set coords [$res2 get {x y z}]	
						
						set cont 0			
						foreach coord $coords {
							set coord2 [vecadd $coord2 $coord]
							incr cont
						}
						set coord2 [vecscale [expr 1.0 / $cont] $coord2]

						#now we will draw the dashed line
						set a [lindex $coord1 0]
						set b [lindex $coord2 0]
	
			
						#to select pink color for contacts
						graphics top color 9
						graphics top line $coord1 $coord2 width 5 style dashed
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont2] == -1} {
							lappend indices $cont2
						}			
					}							
					incr cont2
	
				}
	
			}
			
			#######################################
			##
			##To display Arg - Arg interactions
			##
			######################################
			if {$RIPMD::displayArgArg == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CZ"]
				set coord1 [$res1 get {x y z}]	
				foreach item [lindex $RIPMD::matrixArgArg $index] {
					if {$item != 0} {
						#we will set the second coordinate
						set aux [split [$RIPMD::nodes get $cont] ":"]
						set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type CZ"]
						set coord2 [$res2 get {x y z}]	
						
						#now we will draw the dashed line
						set a [lindex $coord1 0]
						set b [lindex $coord2 0]
						#to select yellow color for contacts
						graphics top color 1
						graphics top line $a $b width 5 style dashed
						
						#append to indices list to then draw the new residues
						if {[lsearch $indices $cont] == -1} {
							lappend indices $cont
						}
						
					}
					incr cont
				}
			}
			#######################################
			##
			##To display vdW
			##
			######################################
			if {$RIPMD::displayVdw == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				foreach item [lindex $RIPMD::matrixVdw $index] {
					if {$item != {}} {
						foreach list $item {
							set resids [split [lindex [lindex $item 0] 0] " "]
							set firstAtom [lindex [split [lindex $resids 0] ":"] 3]
							set secondAtom [lindex [split [lindex $resids 2] ":"] 3]
							set aux [split [lindex $resids 0] ":"]
							set res1 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $firstAtom"]
							set coord1 [$res1 get {x y z}]
							if {$coord1 != ""} {
								set aux [split [lindex $resids 2] ":"]
								set res2 [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $secondAtom"]
								set coord2 [$res2 get {x y z}]
								if {$coord2 != ""} {							
									#now we will draw the dashed line
									set a [lindex $coord1 0]
									set b [lindex $coord2 0]
									#to select yellow color for contacts
									graphics top color 11
									graphics top line $a $b width 5 style dashed
								}
								#append to indices list to then draw the new residues
								if {[lsearch $indices $cont] == -1} {
									lappend indices $cont
								}								
							}
						}
					}
					incr cont
				}
			}
##########################################################################################################################################			
			###################################
			##
			##To display Coulomb interactions
			##
			###################################			
			if {$RIPMD::displayCoulomb == 1} {
				set cont 0
				#we will set the first coordinate
				set aux [split [$RIPMD::nodes get $index] ":"]
				foreach item [lindex $RIPMD::matrixCoulomb $index] {
					set coord1 [veczero]
					set coord2 [veczero]
					if {$item != {}} {
						foreach list $item {
							set resids [split [lindex [lindex $item 0] 0] " "]
							set firstAtoms [string map {"," " "} [split [lindex [split [lindex $resids 0] ":"] 3]]]
							set secondAtoms [string map {"," " "} [split [lindex [split [lindex $resids 2] ":"] 3]]]
							#now we have two list of atoms groups forming the interaction
					
							set aux [split [lindex $resids 0] ":"]
							
							#to set the first coord
							set coord1 [veczero]
							set cont1 0
							foreach atom $firstAtoms {
								set coords ""
								set sel [atomselect top "$RIPMD::toUse [lindex $aux 1] and resid [lindex $aux 2] and type $atom"]
								set coords [$sel get {x y z}]	
								if {$coords != ""} {
									foreach coord $coords {
										set coord1 [vecadd $coord1 $coord]
										incr cont1
									}
								}
							}
							if {$cont1 != 0} {
								#getting the average of coord1
								set coord1 [vecscale [expr 1.0 / $cont1] $coord1]
								#to set coord2
								set coord2 [veczero]
								set cont1 0
								set aux2 [split [lindex $resids 2] ":"]
								foreach atom $secondAtoms {
									set coords ""
									set sel [atomselect top "$RIPMD::toUse [lindex $aux2 1] and resid [lindex $aux2 2] and type $atom"]
									set coords [$sel get {x y z}]	
									if {$coords != ""} {
										foreach coord $coords {
											set coord2 [vecadd $coord2 $coord]
											incr cont1
										}
									}
								}
								if {$cont1 != 0} {
									set coord2 [vecscale [expr 1.0 / $cont1] $coord2]
									#to select tan color for contacts
									graphics top color 5
									graphics top line $coord1 $coord2 width 5 style dashed	
																		
								}
								#append to indices list to then draw the new residues
								if {[lsearch $indices $cont] == -1} {
									lappend indices $cont
								}																
							}
						}
					}
					incr cont
				}
			}
###################################################################################################################														
		}
	}
	RIPMD::drawNotSelectedResidues $indices
	
}

#
# Function to select nodes and make bonds
#
proc RIPMD::selectionMade {w} {
	
	set chain {}
    # --- loop through each selected element
    foreach index [$w curselection] {
		set numbersOfRep [molinfo top get numreps] 
		for { set i 0 } { $i < $numbersOfRep} {incr i} {
			mol delrep end top
		}	
		#creating new representations and setting protein style
		mol addrep top
		mol modselect 0 top protein
		mol modstyle 0 top trace 0.2
		mol modmaterial 0 top Transparent
		mol modcolor 0 top colorid 8

		mol addrep top
		mol modselect 1 top alpha
		mol modstyle 1 top vdw 0.5
		mol modcolor 1 top colorid 1
		set aux [split [$w get $index] ":"]
		#if chain is not in the chain variable we will append it
		if {[string match *[lindex $aux 1]* $chain] != 1} {
			append chain "[lindex $aux 1] "
		}		
    }
    #now we will append to a list of list the different residues
    set resids {}
    set resids [lrepeat [llength $chain] {}]
    set selection ""
    foreach index [$w curselection] {
		#setting a list of residues to append in order of the chain
		set aux [split [$w get $index] ":"]
		#now we will look for the actual chain to append the residue
		set position [lsearch $chain [lindex $aux 1]]
		set auxList [lindex $resids $position]
		append auxList "[lindex $aux 2] "
		set resids [lreplace $resids $position $position $auxList]
	}
	#at this point the selection is configured
	set countChain 0
	#to use is to know if we have chain or segid

	foreach item $chain {
		
		if {$countChain != [expr [llength $chain] - 1]} {
			append selection "($RIPMD::toUse $item and resid [lindex $resids $countChain] ) or "
		} else {
			append selection "($RIPMD::toUse $item and resid [lindex $resids $countChain] )"
		}

		incr countChain
	}
	set sel [atomselect top $selection]
	if { [llength [$sel getbonds]] == 0 } {
		set RIPMD::toUse "chain"
		set selection [string map {segid chain} $selection]
	}	
	set RIPMD::selection "$selection"
	mol addrep top
	mol modselect 2 top "$selection"
	mol modstyle 2 top Bonds 0.2	
	
    set RIPMD::nodes $w
    RIPMD::drawBonds
}

#
#function to fill the matrices
#
proc RIPMD::fillMatrices {} {
	
	set fp [open "$RIPMD::resultFolder/RIP-MD_Results/Graphs/consensus_as_list" r]
	set file_data [read $fp]
	close $fp			
	##  Process data file
	set data2 [split $file_data "\n"]
	foreach line2 $data2 {
		set line3 [split $line2 " "]
		if { $line3 != "" } {
			#adding interaction to the matrices
			set matrixCoord1 [dict get $RIPMD::residIndexDict [lindex $line3 0]]
			set matrixCoord2 [dict get $RIPMD::residIndexDict [lindex $line3 1]]

			#alpha carbon			
			if {[string match *-(Ca)-* $line2] == 1} {
				set vector [lindex $RIPMD::matrixCalpha $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]

				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [incr actualValue]]
				#two matrixcoord1 because is the start and end of the position   
				set RIPMD::matrixCalpha [lreplace $RIPMD::matrixCalpha $matrixCoord1 $matrixCoord1 $vector]
	        
				set vector [lindex $RIPMD::matrixCalpha $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]

				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [incr actualValue]]
				set RIPMD::matrixCalpha [lreplace $RIPMD::matrixCalpha $matrixCoord2 $matrixCoord2 $vector]				
			} 
			
			#HBond
			if {[string match *-(HB)-* $line2] == 1} {
				#getting the interacting atoms in the HBond
				set interaction [split $line2 "'"]
				if {[string match *time* $line2] == 1} {
					set interactionName [lindex $interaction [expr [llength $interaction] - 6]]
					set interaction [split $interactionName "-"]
				} else {
					set interactionName [lindex $interaction [expr [llength $interaction] - 2]]
					set interaction [split $interactionName "-"]
				}
				set atomList {}
				lappend atomList $interaction

				
				#filling the NxNxM matrix
				set vector [lindex $RIPMD::matrixHbond $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [lappend actualValue $atomList]]
				set RIPMD::matrixHbond [lreplace $RIPMD::matrixHbond $matrixCoord1 $matrixCoord1 $vector]			
				
				set vector [lindex $RIPMD::matrixHbond $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [lappend actualValue $atomList]]
				set RIPMD::matrixHbond [lreplace $RIPMD::matrixHbond $matrixCoord2 $matrixCoord2 $vector]						
			}
 
 			#Salt Bridges
			if {[string match *-(salt)-* $line2] == 1} {
				#getting the interacting atoms in the Salt bridge
				set interaction [split $line2 "'"]
				if {[string match *time* $line2] == 1} {
					set interactionName [lindex $interaction [expr [llength $interaction] - 6]]
					set interaction [split $interactionName "-"]					
				} else {
					set interactionName [lindex $interaction [expr [llength $interaction] - 2]]
					set interaction [split $interactionName "-"]
				}
				set atomList {}
				lappend atomList $interaction

				
				#filling the NxNxM matrix
				set vector [lindex $RIPMD::matrixSalt $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [lappend actualValue $atomList]]
				set RIPMD::matrixSalt [lreplace $RIPMD::matrixSalt $matrixCoord1 $matrixCoord1 $vector]			
				
				set vector [lindex $RIPMD::matrixSalt $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [lappend actualValue $atomList]]
				set RIPMD::matrixSalt [lreplace $RIPMD::matrixSalt $matrixCoord2 $matrixCoord2 $vector]						
			}           
			
			#disulphide bridges
			if {[string match *-(SS)-* $line2] == 1} {
				set vector [lindex $RIPMD::matrixDisulphide $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [incr actualValue]]
				set RIPMD::matrixDisulphide [lreplace $RIPMD::matrixDisulphide $matrixCoord1 $matrixCoord1 $vector]
	
				set vector [lindex $RIPMD::matrixDisulphide $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [incr actualValue]]
				set RIPMD::matrixDisulphide [lreplace $RIPMD::matrixDisulphide $matrixCoord2 $matrixCoord2 $vector]				
			}	
			
			#cation pi interaction
			if {[string match *-(cation-pi)-* $line2] == 1} {
				
				#getting the interacting atoms in the cation pi interaction
				set interaction [split $line2 "'"]
				if {[string match *time* $line2] == 1} {
					set interactionName [lindex $interaction [expr [llength $interaction] - 6]]
					set interaction [split $interactionName "-"]				
				} else {		
					set interactionName [lindex $interaction [expr [llength $interaction] - 6]]
					set interaction [split $interactionName "-"]
				}
				set atomList {}
				lappend atomList $interaction
                
				
				#filling the NxNxM matrix
				set vector [lindex $RIPMD::matrixCationPi $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [lappend actualValue $atomList]]
				set RIPMD::matrixCationPi [lreplace $RIPMD::matrixCationPi $matrixCoord1 $matrixCoord1 $vector]			
				
				set vector [lindex $RIPMD::matrixCationPi $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [lappend actualValue $atomList]]
				set RIPMD::matrixCationPi [lreplace $RIPMD::matrixCationPi $matrixCoord2 $matrixCoord2 $vector]		
			}
			
			#pi-pi interaction
			if {[string match *-(pi-pi)-* $line2] == 1} {
				set vector [lindex $RIPMD::matrixPiPi $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [incr actualValue]]
				set RIPMD::matrixPiPi [lreplace $RIPMD::matrixPiPi $matrixCoord1 $matrixCoord1 $vector]
				
				set vector [lindex $RIPMD::matrixPiPi $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [incr actualValue]]
				set RIPMD::matrixPiPi [lreplace $RIPMD::matrixPiPi $matrixCoord2 $matrixCoord2 $vector]
				
			}
			

 			#Arg - Arg interaction			
			if {[string match *-(Arg-Arg)-* $line2] == 1} {
				set vector [lindex $RIPMD::matrixArgArg $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [incr actualValue]]
				#two matrixcoord1 because is the start and end of the position   
				set RIPMD::matrixArgArg [lreplace $RIPMD::matrixArgArg $matrixCoord1 $matrixCoord1 $vector]
	        
				set vector [lindex $RIPMD::matrixArgArg $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [incr actualValue]]
				set RIPMD::matrixArgArg [lreplace $RIPMD::matrixArgArg $matrixCoord2 $matrixCoord2 $vector]				
			}    
            
 			#Coulomb
			if {[string match *-(coulomb)-* $line2] == 1} {
				
				set interaction [split $line2 "'"]
				
				if {[string match *time* $line2] == 1} {
					set interactionName [lindex $interaction [expr [llength $interaction] - 10]]
					set interaction [split $interactionName "-"]	
					
				} else {
					set interactionName [lindex $interaction [expr [llength $interaction] - 2]]
					set interaction [split $interactionName "-"]
				}
				#set interactionName [lindex $interaction [expr [llength $interaction] - 2]]
				#set interaction [split $interactionName "-"]
				set atomList {}
				lappend atomList $interaction
				#filling the NxNxM matrix
				set vector [lindex $RIPMD::matrixCoulomb $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [lappend actualValue $atomList]]
				set RIPMD::matrixCoulomb [lreplace $RIPMD::matrixCoulomb $matrixCoord1 $matrixCoord1 $vector]			
				
				set vector [lindex $RIPMD::matrixCoulomb $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [lappend actualValue $atomList]]
				set RIPMD::matrixCoulomb [lreplace $RIPMD::matrixCoulomb $matrixCoord2 $matrixCoord2 $vector]						
			} 
			
 			#vdW
			if {[string match *-(vdW)-* $line2] == 1} {
				#getting the interacting atoms in the vdW
				set interaction [split $line2 "'"]
				if {[string match *time* $line2] == 1} {
					set interactionName [lindex $interaction [expr [llength $interaction] - 10]]
					set interaction [split $interactionName "-"]
				} else {
					set interactionName [lindex $interaction [expr [llength $interaction] - 2]]
					set interaction [split $interactionName "-"]
				}
				set atomList {}
				lappend atomList $interaction
            
				
				#filling the NxNxM matrix
				set vector [lindex $RIPMD::matrixVdw $matrixCoord1]
				set actualValue [lindex $vector $matrixCoord2]
				set vector [lreplace $vector $matrixCoord2 $matrixCoord2 [lappend actualValue $atomList]]
				set RIPMD::matrixVdw [lreplace $RIPMD::matrixVdw $matrixCoord1 $matrixCoord1 $vector]			
				
				set vector [lindex $RIPMD::matrixVdw $matrixCoord2]
				set actualValue [lindex $vector $matrixCoord1]
				set vector [lreplace $vector $matrixCoord1 $matrixCoord1 [lappend actualValue $atomList]]
				set RIPMD::matrixVdw [lreplace $RIPMD::matrixVdw $matrixCoord2 $matrixCoord2 $vector]						
			} 			                
							
		}
	}
	
}


#
# function to initialice each matrix filled with 0 or with a list; residIndexDict is a dict  with node and id values
#
proc RIPMD::initialMatrices {} {
	set length [dict size $RIPMD::residIndexDict]
	set RIPMD::matrixCalpha [lrepeat $length [lrepeat $length 0]]
	set RIPMD::matrixHbond [lrepeat $length [lrepeat $length {}]]
	set RIPMD::matrixSalt [lrepeat $length [lrepeat $length {}]]
	set RIPMD::matrixDisulphide [lrepeat $length [lrepeat $length 0]]
	set RIPMD::matrixCationPi [lrepeat $length [lrepeat $length {}]]
	set RIPMD::matrixPiPi [lrepeat $length [lrepeat $length 0]]
	set RIPMD::matrixArgArg [lrepeat $length [lrepeat $length 0]]
	set RIPMD::matrixCoulomb [lrepeat $length [lrepeat $length {}]]
	set RIPMD::matrixVdw	[lrepeat $length [lrepeat $length {}]]
}

#
# Function to load results
#
proc RIPMD::loadResults {} {
	#we will look for molecules loaded, if there are molecules we will delete it
	if {[catch {set RIPMD::pathToFiles [molinfo top get filename]} errmsg]} {
		set a ""
	} else {
		set molLoaded [mol list]
		foreach item $molLoaded {
			mol delete top
		}
	}
	
	#now we load the pqr stucture (or last frame). the reason is to manage statics coords
	mol load pdb "$RIPMD::resultFolder/RIP-MD_Results/representative.pdb"
	
	
	#enabling checkbuttons to display different types of bonds	
	$RIPMD::w.n.f5.calpha configure -state normal
	$RIPMD::w.n.f5.hbond configure -state normal
	$RIPMD::w.n.f5.salt configure -state normal
	$RIPMD::w.n.f5.disulphide configure -state normal
	$RIPMD::w.n.f5.capi configure -state normal
	$RIPMD::w.n.f5.pipi configure -state normal
	$RIPMD::w.n.f5.argarg configure -state normal
	$RIPMD::w.n.f5.vdw configure -state normal
	$RIPMD::w.n.f5.coulomb configure -state normal
	$RIPMD::w.n.f5.showPearson configure -state normal
	$RIPMD::w.n.f5.showHist configure -state normal
	
	#deleting actual representations
	$RIPMD::w.n.f5.nodesBox.l delete 0 end	
	$RIPMD::w.n.f1.textBox.l insert end "Results loaded, please see the Result Tab"
	$RIPMD::w.n.f1.textBox.l insert end "In that tab you can select residues to display their interactions"
	$RIPMD::w.n.f1.textBox.l insert end "In the same tab you will also be able to generate correlation plots"
	$RIPMD::w.n.f1.textBox.l insert end ""
	$RIPMD::w.n.f1.textBox.l insert end "For more information, please read the README file in $RIPMD::resultFolder"
	
	#we will delete actual representations
	if {[catch {set numbersOfRep [molinfo top get numreps]} errmsg]} {
		set a ""
		
	} else {
		for { set i 0 } { $i < $numbersOfRep} {incr i} {
			mol delrep end top
		}	

	}
#	set numbersOfRep [molinfo top get numreps] 

	
	#creating new representations and setting protein style
	mol addrep top
	mol modselect 0 top protein
	mol modstyle 0 top trace 0.2
	mol modcolor 0 top colorid 8
	mol modmaterial 0 top Transparent
	mol addrep top
	mol modselect 1 top alpha
	mol modstyle 1 top vdw 0.5
	mol modcolor 1 top colorid 1	
	
	
	#settin the dict in white to fill it (to do not have problems if i read others result later)
	set RIPMD::residIndexDict ""
	#reading graph info
	#putting nodes
	set fp2 [open "$RIPMD::resultFolder/RIP-MD_Results/nodes" r]			
	set file_data2 [read $fp2]
	close $fp2				
	##  Process data file
	set data2 [split $file_data2 "\n"]
	set count 0
	set i 0
	foreach line2 $data2 {
		if {$count!=0} {
			set line3 [split $line2 "\t"]
			if { $line3 != "" } {
				#also we will put the information in a dictionary {node:id}
				#dict append RIPMD::residIndexDict [lindex $line3 end] [lindex $line3 0] 
				#dict append RIPMD::indexResidDict [lindex $line3 0] [lindex $line3 end] 
				dict append RIPMD::residIndexDict [lindex $line3 end] $i 
				dict append RIPMD::indexResidDict $i [lindex $line3 end] 
				set i [expr $i + 1]
				$RIPMD::w.n.f5.nodesBox.l insert end "[lindex $line3 end]"
			}
		} else {
			incr count
		}
	}
	
	#creating diffrerent matrices for saving each type of interaction
	#we will call a function to start each matrix filled with 0
	RIPMD::initialMatrices
	
	#now we will call a function to fill each matrix
	RIPMD::fillMatrices	
}


#
# Creating a function to check if some parameters are or not numbers
#
proc RIPMD::checkNumbers {} {
	
	#try to divide numbers by 1, if we catch an error we will append the error
	if {[catch {expr $RIPMD::calphaDist / 1} errmsg]} {
		append RIPMD::errors "* Distance between alpha C's is not numerical.\n"
	}
	if {[catch {expr $RIPMD::hbondDist / 1} errmsg]} {
		append RIPMD::errors "* Distance between H bonds is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::hbondAngle / 1} errmsg]} {
		append RIPMD::errors "* Angle value for H bonds is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::saltDist / 1} errmsg]} {
		append RIPMD::errors "* Distance value for atoms forming a salt bridge is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::pipiDist / 1} errmsg]} {
		append RIPMD::errors "* Distance value for atoms forming a pi - pi interaction is not numerical.\n"
	}
	if {[catch {expr $RIPMD::ArgArgDist / 1} errmsg]} {
		append RIPMD::errors "* Distance value for atoms forming a Arg - Arg interaction is not numerical.\n"
	}		
	if {[catch {expr $RIPMD::simulatedPermittivity / 1} errmsg]} {
		append RIPMD::errors "* Value for simulated permittivity is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::capiDist / 1} errmsg]} {
		append RIPMD::errors "* Distance value for atoms forming a cation - pi interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::capiAngle1 / 1} errmsg]} {
		append RIPMD::errors "* First angle for cation - pi interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::capiAngle2 / 1} errmsg]} {
		append RIPMD::errors "* Second angle for cation - pi interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::capiAngle3 / 1} errmsg]} {
		append RIPMD::errors "* Third angle for cation - pi interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::capiAngle4 / 1} errmsg]} {
		append RIPMD::errors "* Fourth angle for cation - pi interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::disulphideDist / 1} errmsg]} {
		append RIPMD::errors "* Distance value for atoms forming a disulfide bridge interaction is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::disulphideAng1 / 1} errmsg]} {
		append RIPMD::errors "* First angle for disulfide bridge is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::disulphideAng2 / 1} errmsg]} {
		append RIPMD::errors "* Second angle for disulfide bridge is not numerical.\n"
	}
	if {[catch {expr $RIPMD::covalentDist / 1} errmsg]} {
		append RIPMD::errors "* Covalent bond threshold for VdW is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::distVDW1 / 1} errmsg]} {
		append RIPMD::errors "* Left dstance between atoms to compute VdW contacts is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::distVDW2 / 1} errmsg]} {
		append RIPMD::errors "* Right dstance between atoms to compute VdW contacts is not numerical.\n"
	}		
	if {[catch {expr $RIPMD::percentaje / 1} errmsg]} {
		append RIPMD::errors "* Percentaje of time is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::frame_start / 1} errmsg]} {
		append RIPMD::errors "* Frame to start RIN computation is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::frame_end / 1} errmsg]} {
		append RIPMD::errors "* Frame to stop RIN computation is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::frame_separation / 1} errmsg]} {
		append RIPMD::errors "* Number of separation frames is not numerical.\n"
	}	
	if {[catch {expr $RIPMD::frame_separation / 1} errmsg]} {
		append RIPMD::errors "* Number of separation frames is not numerical.\n"
	}
	if {[catch {expr $RIPMD::reactionFieldPermittivity / 1} errmsg]} {
		append RIPMD::errors "* Reaction Field Permittivity is not numerical.\n"
	}
	if {[catch {expr $RIPMD::temp / 1} errmsg]} {
		append RIPMD::errors "* Coulomb temperature fo thermal noise is not numerical.\n"
	}
	if {[catch {expr $RIPMD::inverseDebye / 1} errmsg]} {
		append RIPMD::errors "* Inverse debye screening is not numerical.\n"
	}
	if {[catch {expr $RIPMD::coulombThreshold / 1} errmsg]} {
		append RIPMD::errors "* Coulomb Distance Threshold is not numerical.\n"
	}
	if {[catch {expr $RIPMD::coulombCovalentDist / 1} errmsg]} {
		append RIPMD::errors "* Coulomb covalent bonds distance cutoff is not numerical.\n"
	}

	if {[catch {expr $RIPMD::pH / 1} errmsg]} {
		append RIPMD::errors "* pH is not numerical.\n"
	}	
	
	if {[catch {expr $RIPMD::proc / 1} errmsg]} {
		append RIPMD::errors "* Number of processor is not numerical.\n"
	} else {
		if {$RIPMD::proc < 1} {		
			append RIPMD::errors "* Number of processors to use can not be less than 1.\n"
		}
		
	}	
	return
}


#
# Create the window and initialize data structures
#
proc RIPMD::ripmd {} {
	variable w
	# If already initialized, just turn on
	if { [winfo exists .ripmd] } {
		wm deiconify $w
		return
	}
	##############################
	##
	## CHARMM variables
	##
	#############################
	


	if {[catch {set RIPMD::forceField "$::env(RIP_MD)/dat/top_all22_prot.rtf"} errmsg]} {
		set a [tk_dialog .myDialog "RIP-MD" "Can not locate RIP-MD system variable. Have you installed the RIP-MD standalone core?" error 0 "Ok"]
	} else {
		set RIPMD::forceField "$::env(RIP_MD)/dat/top_all22_prot.rtf"
		set RIPMD::parameterFile "$::env(RIP_MD)/dat/par_all22_prot.prm"		
	}


	
	set w [toplevel ".ripmd"]
	#
	wm title $w "RIP-MD - Residue Interactions in Protein Molecular Dynamic"
	wm maxsize $w 1024 600
	wm minsize $w 1024 600
	wm resizable $w 0 0
	
	## make the menu bar
	frame $w.menubar -relief raised -bd 2 ;#frame for menu bar
	pack $w.menubar -padx 1 -fill x
	
	#file button
	menubutton $w.menubar.file -text File -underline 0 -menu $w.menubar.file.menu
	
	#help button
	menubutton $w.menubar.help -text Help -underline 0 -menu $w.menubar.help.menu
	
	##
	## file menu
	##
	
	menu $w.menubar.file.menu -tearoff no
	$w.menubar.file.menu add command -label "Open results" -command {
		
		set a [tk_dialog .myDialog "RIP-MD" "Please, select a folder with all result files" warning 0 "Ok"]
		set aux [tk_chooseDirectory]
		if {$aux eq ""} {
			set a [tk_dialog .myDialog "RIP-MD" "No directory selected" error 0 "Ok"]
		} else {
			set RIPMD::resultFolder "$aux"
			RIPMD::loadResults
		}

	}
	
	
	#$w.menubar.help.menu add command -label "About RIP-MD" -command "set a [tk_dialog .myDialog 'About RIP-MD' 'asdasdasd' questhead 0 'Ok']"
	# XXX - set menubutton width to avoid truncation in OS X
	$w.menubar.file config -width 5
	
	##
	## help menu
	##
	menu $w.menubar.help.menu -tearoff no
	$w.menubar.help.menu add command -label "Help" -command {set aux [tk_dialog .myDialog "Help" "For detailed information, please refer to the RIP-MD user manual. You can find it at http://www.dlab.cl/ripmd" questhead 0 "Ok"]}
	#$w.menubar.help.menu add command -label "About RIP-MD" -command "set a [tk_dialog .myDialog 'About RIP-MD' 'asdasdasd' questhead 0 'Ok']"
	# XXX - set menubutton width to avoid truncation in OS X
	$w.menubar.help config -width 5
	
	
	#option tabs
	ttk::notebook $w.n -width 1300 -height 500
	ttk::frame $w.n.f1;
	ttk::frame $w.n.f2;
	ttk::frame $w.n.f3;
	ttk::frame $w.n.f5;

	$w.n add $w.n.f1 -text "General Options"
	$w.n add $w.n.f2 -text "Interaction Options"
	$w.n add $w.n.f3 -text "Advanced Options"
	
	#####################################
	##
	## General Options
	##
	#####################################
	grid [label $w.n.f1.interactionText -text "Select interactions to compute"]
	place $w.n.f1.interactionText -x 10 -y 10

	grid  [checkbutton $w.n.f1.calpha -text "Compute Contact Map (C Alpha)" -variable RIPMD::calpha] -rowspan 200
	place $w.n.f1.calpha -x 20 -y 40
	grid [checkbutton $w.n.f1.hbond -text "Compute H Bonds" -variable RIPMD::hbonds] -rowspan 50
	place $w.n.f1.hbond -x 20 -y 60
	grid [checkbutton $w.n.f1.salt -text "Compute Salt Bridges" -variable RIPMD::salt] -rowspan 50
	place $w.n.f1.salt -x 20 -y 80 ;#I dont know why but this put in a correctly position all checkbuttons
	grid [checkbutton $w.n.f1.disulphide -text "Compute Disulfide Bridges" -variable RIPMD::disulphide] -rowspan 50
	place $w.n.f1.disulphide -x 20 -y 100 
	grid [checkbutton $w.n.f1.capi -text "Compute Cation - pi interactions" -variable RIPMD::capi] -rowspan 50
	place $w.n.f1.capi -x 20 -y 120
	#\u3c0
	grid [checkbutton $w.n.f1.pipi -text "Compute pi - pi interactions" -variable RIPMD::pipi] -rowspan 50
	place $w.n.f1.pipi -x 20 -y 140
	grid  [checkbutton $w.n.f1.argarg -text "Compute Arg - Arg interactions" -variable RIPMD::ArgArg] -rowspan 200
	place $w.n.f1.argarg -x 20 -y 160	
	grid [checkbutton $w.n.f1.vdw -text "Compute VdW contacts *" -variable RIPMD::vdw] -rowspan 50
	place $w.n.f1.vdw -x 20 -y 180
	grid [checkbutton $w.n.f1.coulomb -text "Compute Coulomb interactions *" -variable RIPMD::coulomb] -rowspan 50
	place $w.n.f1.coulomb -x 20 -y 200
	grid [checkbutton $w.n.f1.pearson -text "Compute Pearson correlation" -variable RIPMD::pearson] -rowspan 50
	place $w.n.f1.pearson -x 20 -y 220	
	grid [label $w.n.f1.pearson2 -text "(only for MD)" ] -rowspan 50
	place $w.n.f1.pearson2 -x 40 -y 240	
	
	
	grid [label $w.n.f1.selectionToComputeText -text "Selection to compute"]
	place $w.n.f1.selectionToComputeText -x 10 -y 280
	grid  [entry $w.n.f1.selToCompute -width 15 -justify left -textvariable RIPMD::proteinSelection]
	place $w.n.f1.selToCompute -x 180 -y 280
	grid [label $w.n.f1.selectionToComputeText2 -text "(It use a syntax very similar to CHARMM."]
	place $w.n.f1.selectionToComputeText2 -x 10 -y 300	
	grid [label $w.n.f1.selectionToComputeText3 -text "For more information please read the User manual )"]
	place $w.n.f1.selectionToComputeText3 -x 10 -y 315		
	
	grid [label $w.n.f1.warningText1 -text "* Coulomb potential and VdW contacts could take long time to be computed."] 
	place $w.n.f1.warningText1 -x 0 -y 465	
	grid [label $w.n.f1.warningText2 -text "  These options could use a high amount of RAM memory, depending on the number of atoms of your structure."]
	place $w.n.f1.warningText2 -x 0 -y 480

	#button to select the output folder
	pack [button $w.n.f1.selectButton -text "Browse..." -width 10 -command {
		set RIPMD::dirName [tk_chooseDirectory]
		set RIPMD::outputFolder $RIPMD::dirName
		set position [llength [split $RIPMD::dirName /] ]
		set aux [split $RIPMD::dirName / ]
		#puts $aux
		if {$position==0} {
			set RIPMD::dirName "Select output folder for results..."
			set RIPMD::outputFolder "*"
		} else {
			set RIPMD::dirName "Folder:  [lindex $aux [expr $position - 1]]"
			
		}
	}]
	place $w.n.f1.selectButton -x 310 -y 80
	pack [label $w.n.f1.dirPath -text "Select output folder for results..." -textvariable RIPMD::dirName]
	place $w.n.f1.dirPath -x 310 -y 50
	

	#Output options
	grid [label $w.n.f1.outputText -text "Output options"]
	place $w.n.f1.outputText -x 300 -y 10
	
	grid [label $w.n.f1.outputFormatText -text "Output format"]
	place $w.n.f1.outputFormatText -x 310 -y 120
	
	ttk::menubutton $w.n.f1.outFormatMenu -menu $w.n.f1.outFormatMenu.menu  -textvariable RIPMD::format -width 20
	menu $w.n.f1.outFormatMenu.menu
	$w.n.f1.outFormatMenu.menu add command -label "GML" -command { set RIPMD::format "GML" }
	$w.n.f1.outFormatMenu.menu add command -label "EdgeList" -command { set RIPMD::format "edgelist" }
	$w.n.f1.outFormatMenu.menu add command -label "GraphML" -command { set RIPMD::format "GraphML" }
	$w.n.f1.outFormatMenu.menu add command -label "Pajek" -command { set RIPMD::format "Pajek" }
	
	pack $w.n.f1.outFormatMenu
	place $w.n.f1.outFormatMenu -x 310 -y 150
	
	
	#####################################
	##
	## log widget
	##
	####################################
	
	
	frame $w.n.f1.textBox
	grid [tk::listbox $w.n.f1.textBox.l -yscrollcommand "$w.n.f1.textBox.s1 set" -xscrollcommand "$w.n.f1.textBox.s2 set" -height 20 -width 50] -column 0 -row 0 -sticky nwes
	grid [ttk::scrollbar $w.n.f1.textBox.s1 -command "$w.n.f1.textBox.l yview" -orient vertical   ] -column 1 -row 0 -sticky ns
	grid [ttk::scrollbar $w.n.f1.textBox.s2 -command "$w.n.f1.textBox.l xview" -orient horizontal ] -column 0 -row 1 -sticky ew

	$RIPMD::w.n.f1.textBox.l insert end "When you run RIP-MD, the log file will be placed in the"		
	$RIPMD::w.n.f1.textBox.l insert end "Output folder."		
	$RIPMD::w.n.f1.textBox.l insert end "In addition you will be able to see it in the VMD console."
	$RIPMD::w.n.f1.textBox.l insert end ""		
	
	

	pack $w.n.f1.textBox
	place $w.n.f1.textBox -x 600 -y 30
	
	######################################
	##
	## Interaction options
	##
	######################################

	#calpha

	
	grid  [label $w.n.f2.calphaOpt -text "Options for C alpha contacts"] -rowspan 50
	place $w.n.f2.calphaOpt -x 0 -y 20
	grid  [label $w.n.f2.calpha -text "Distance cuttoff between C alpha"] -rowspan 50
	place $w.n.f2.calpha -x 20 -y 40
	grid  [entry $w.n.f2.calphaDist -width 8 -justify right -textvariable RIPMD::calphaDist]
	place $w.n.f2.calphaDist -x 250 -y 40
	grid  [label $w.n.f2.calpha2 -text "\u212b"] -rowspan 50
	place $w.n.f2.calpha2 -x 315 -y 40
	
	#hbonds
	grid  [label $w.n.f2.hBondOpt -text "Options for HBonds"] -rowspan 50
	place $w.n.f2.hBondOpt -x 0 -y 80
	grid  [label $w.n.f2.hbondDistText -text "Distance between donnor/acceptor"] -rowspan 50
	place $w.n.f2.hbondDistText -x 20 -y 100
	grid  [label $w.n.f2.hbondDistText2 -text "atoms"] -rowspan 50
	place $w.n.f2.hbondDistText2 -x 20 -y 115	
	grid  [entry $w.n.f2.hBondDist -width 8 -justify right -textvariable RIPMD::hbondDist]
	place $w.n.f2.hBondDist -x 250 -y 100
	grid  [label $w.n.f2.hbondText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.hbondText2 -x 315 -y 100
	
	grid  [label $w.n.f2.hbondAngleText -text "Angle between donnor/acceptor"] -rowspan 50
	place $w.n.f2.hbondAngleText -x 20 -y 140
	grid  [label $w.n.f2.hbondAngleText2 -text "vector"] -rowspan 50
	place $w.n.f2.hbondAngleText2 -x 20 -y 155
	grid  [entry $w.n.f2.hBondAngle -width 8 -justify right -textvariable RIPMD::hbondAngle]
	place $w.n.f2.hBondAngle -x 250 -y 140
	grid  [label $w.n.f2.hbondText3 -text "°"] -rowspan 50
	place $w.n.f2.hbondText3 -x 315 -y 140

	#salt bridges
	
	grid  [label $w.n.f2.saltOpt -text "Options for Salt Bridges"] -rowspan 50
	place $w.n.f2.saltOpt -x 0 -y 180
	grid  [label $w.n.f2.salt -text "Distance cuttoff between atoms"] -rowspan 50
	place $w.n.f2.salt -x 20 -y 200	
	grid  [entry $w.n.f2.saltDist -width 8 -justify right -textvariable RIPMD::saltDist]
	place $w.n.f2.saltDist -x 250 -y 200	
	grid  [label $w.n.f2.saltText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.saltText2 -x 315 -y 200	

	#pi pi interaction

	grid  [label $w.n.f2.pipiOpt -text "Options for pi - pi interactions"] -rowspan 50
	place $w.n.f2.pipiOpt -x 0 -y 240
	grid  [label $w.n.f2.pipi -text "Distance between ring center"] -rowspan 50
	place $w.n.f2.pipi -x 20 -y 260
	grid  [label $w.n.f2.pipi2 -text "of mass "] -rowspan 50
	place $w.n.f2.pipi2 -x 20 -y 275
	grid  [entry $w.n.f2.pipiDistance -width 8 -justify right -textvariable RIPMD::pipiDist]
	place $w.n.f2.pipiDistance -x 250 -y 260
	grid  [label $w.n.f2.pipiText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.pipiText2 -x 315 -y 260

	#arg arg interaction

	grid  [label $w.n.f2.argargOpt -text "Options for Arg - Arg interactions"] -rowspan 50
	place $w.n.f2.argargOpt -x 0 -y 300
	grid  [label $w.n.f2.argarg -text "Distance between Cz's"] -rowspan 50
	place $w.n.f2.argarg -x 20 -y 320
	grid  [entry $w.n.f2.argargDistance -width 8 -justify right -textvariable RIPMD::ArgArgDist]
	place $w.n.f2.argargDistance -x 250 -y 320
	grid  [label $w.n.f2.argargText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.argargText2 -x 315 -y 320

	#coulomb
	grid  [label $w.n.f2.couOpt -text "Options for Coulomb Potential contacts"] -rowspan 50
	place $w.n.f2.couOpt -x 0 -y 360
	grid  [checkbutton $w.n.f2.react -text "Use reaction field" -variable RIPMD::reactionField] -rowspan 200
	place $w.n.f2.react -x 0 -y 380

    
	grid  [label $w.n.f2.couOpt3 -text "Reaction field permittivity"] -rowspan 50
	place $w.n.f2.couOpt3 -x 20 -y 400
	grid  [entry $w.n.f2.dielectric2  -width 8 -justify right -textvariable RIPMD::reactionFieldPermittivity]
	place $w.n.f2.dielectric2 -x 250 -y 400
	
	grid  [label $w.n.f2.couT -text "Temperature for thermal noise"] -rowspan 50
	place $w.n.f2.couT -x 20 -y 420
	grid  [label $w.n.f2.couT2 -text "cutoff"] -rowspan 50
	place $w.n.f2.couT2 -x 20 -y 435
	grid  [entry $w.n.f2.couT3  -width 8 -justify right -textvariable RIPMD::temp]
	place $w.n.f2.couT3 -x 250 -y 420
	grid  [label $w.n.f2.couT4 -text "K"] -rowspan 50
	place $w.n.f2.couT4 -x 315 -y 420
		

	###
	###
	### entre items separacion de 40, entre opciones esde 20
	###
	###
	###
	grid  [label $w.n.f2.couOpt2 -text "Simulated permittivity"] -rowspan 50
	place $w.n.f2.couOpt2 -x 390 -y 20
	grid  [entry $w.n.f2.dielectric  -width 8 -justify right -textvariable RIPMD::simulatedPermittivity]
	place $w.n.f2.dielectric -x 645 -y 20
	
	grid  [label $w.n.f2.couInv -text "Inverse debye screening"] -rowspan 50
	place $w.n.f2.couInv -x 390 -y 40
	grid  [entry $w.n.f2.couInv2  -width 8 -justify right -textvariable RIPMD::inverseDebye]
	place $w.n.f2.couInv2 -x 645 -y 40
		
	grid  [label $w.n.f2.couThr -text "Distance threshold"] -rowspan 50
	place $w.n.f2.couThr -x 390 -y 60
	grid  [entry $w.n.f2.couThr2  -width 8 -justify right -textvariable RIPMD::coulombThreshold]
	place $w.n.f2.couThr2 -x 645 -y 60
	grid  [label $w.n.f2.couThr3 -text "\u212b"] -rowspan 50
	place $w.n.f2.couThr3 -x 710 -y 60		

	grid  [label $w.n.f2.couBonds -text "Covalent bonds distance cutoff"] -rowspan 50
	place $w.n.f2.couBonds -x 390 -y 80	
	grid  [entry $w.n.f2.couBonds2  -width 8 -justify right -textvariable RIPMD::coulombCovalentDist]
	place $w.n.f2.couBonds2 -x 645 -y 80
	
	#cation pi interaction 

	grid  [label $w.n.f2.capiOpt -text "Option for Cation - pi interactions"] -rowspan 50
	place $w.n.f2.capiOpt -x 370 -y 120
	grid  [label $w.n.f2.capi -text "Distance between cation and ring "] -rowspan 50
	place $w.n.f2.capi -x 390 -y 140
	grid  [label $w.n.f2.capi2 -text "center of mass"] -rowspan 50
	place $w.n.f2.capi2 -x 390 -y 155
	grid  [entry $w.n.f2.capiDist -width 8 -justify right -textvariable RIPMD::capiDist]
	place $w.n.f2.capiDist -x 645 -y 140
	grid  [label $w.n.f2.capiText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.capiText2 -x 710 -y 140

	grid  [label $w.n.f2.capi3 -text "Angle between center of mass and"] -rowspan 50
	place $w.n.f2.capi3 -x 390 -y 175
	grid  [label $w.n.f2.capi4 -text "normal vector of aromatic ring"] -rowspan 50
	place $w.n.f2.capi4 -x 390 -y 190	
	#grid  [label $w.n.f2.capi5 -text "aromatic ring."] -rowspan 50
	#place $w.n.f2.capi5 -x 390 -y 205	
	
	grid  [entry $w.n.f2.capiAngle1 -width 8 -justify right -textvariable RIPMD::capiAngle1]
	place $w.n.f2.capiAngle1 -x 645 -y 175		
	
	grid  [label $w.n.f2.capi6 -text "°\tand"] -rowspan 50
	place $w.n.f2.capi6 -x 710 -y 175

	grid  [entry $w.n.f2.capiAngle2 -width 8 -justify right -textvariable RIPMD::capiAngle2]
	place $w.n.f2.capiAngle2 -x 825 -y 175
	
	grid  [label $w.n.f2.capi7 -text "°"] -rowspan 50
	place $w.n.f2.capi7 -x 890 -y 175
	
	grid  [label $w.n.f2.capi8 -text "\tor"] -rowspan 50
	place $w.n.f2.capi8 -x 716 -y 195

	grid  [entry $w.n.f2.capiAngle3 -width 8 -justify right -textvariable RIPMD::capiAngle3]
	place $w.n.f2.capiAngle3 -x 645 -y 215		
	
	grid  [label $w.n.f2.capi9 -text "°\tand"] -rowspan 50
	place $w.n.f2.capi9 -x 710 -y 215

	grid  [entry $w.n.f2.capiAngle4 -width 8 -justify right -textvariable RIPMD::capiAngle4]
	place $w.n.f2.capiAngle4 -x 825 -y 215
	grid  [label $w.n.f2.capi10 -text "°"] -rowspan 50
	place $w.n.f2.capi10 -x 890 -y 215
				
	#disulphide bridge
	grid  [label $w.n.f2.disOpt -text "Options for Disulfide Bridges"] -rowspan 50
	place $w.n.f2.disOpt -x 370 -y 255
	grid  [label $w.n.f2.dis -text "Distance cuttoff between S atoms"] -rowspan 50
	place $w.n.f2.dis -x 390 -y 275


	grid  [entry $w.n.f2.disDist -width 8 -justify right -textvariable RIPMD::disulphideDist]
	place $w.n.f2.disDist -x 645 -y 275
	
	grid  [label $w.n.f2.disText2 -text "\u212b"] -rowspan 50
	place $w.n.f2.disText2 -x 710 -y 275
		
	
	grid  [label $w.n.f2.textDisAng -text "Range of Angles between C-S vectors"] -rowspan 50
	place $w.n.f2.textDisAng -x 390 -y 295
	
	grid  [entry $w.n.f2.disAng1 -width 8 -justify right -textvariable RIPMD::disulphideAng1]
	place $w.n.f2.disAng1 -x 645 -y 295
	
	grid  [label $w.n.f2.textDisAng3 -text "°\tand"] -rowspan 50
	place $w.n.f2.textDisAng3 -x 710 -y 295
	
	grid  [entry $w.n.f2.disAng2 -width 8 -justify right -textvariable RIPMD::disulphideAng2]
	place $w.n.f2.disAng2 -x 825 -y 295	

	grid  [label $w.n.f2.textDisAng4 -text "°"] -rowspan 50
	place $w.n.f2.textDisAng4 -x 890 -y 295		

	#vdw
	grid  [label $w.n.f2.vdwOpt -text "Options for VdW contacts"] -rowspan 50
	place $w.n.f2.vdwOpt -x 370 -y 335
	grid  [label $w.n.f2.vdwOpt2 -text "Covalent bonds distance cutoff"] -rowspan 50
	place $w.n.f2.vdwOpt2 -x 390 -y 355
	
	grid  [entry $w.n.f2.covDist -width 8 -justify right -textvariable RIPMD::covalentDist]
	place $w.n.f2.covDist -x 645 -y 355	


	grid  [label $w.n.f2.vdwOpt3 -text "Range of distance between two atoms"] -rowspan 50
	place $w.n.f2.vdwOpt3 -x 390 -y 375
		
	grid  [entry $w.n.f2.distVDW1 -width 8 -justify right -textvariable RIPMD::distVDW1]
	place $w.n.f2.distVDW1 -x 645 -y 375			
	
	grid  [label $w.n.f2.vdwOpt4 -text "\u212b\tand"] -rowspan 50
	place $w.n.f2.vdwOpt4 -x 710 -y 375	

	grid  [entry $w.n.f2.distVDW2 -width 8 -justify right -textvariable RIPMD::distVDW2]
	place $w.n.f2.distVDW2 -x 825 -y 375	

	grid  [label $w.n.f2.vdwOpt5 -text "\u212b"] -rowspan 50
	place $w.n.f2.vdwOpt5 -x 890 -y 375	

	#######################################
	##
	## advanced  options
	##
	#######################################
	grid [label $w.n.f3.pdbText -text "PDB options"]
	place $w.n.f3.pdbText -x 0 -y 20
	grid  [checkbutton $w.n.f3.missing -text "Use PDB2PQR to add missing atoms" -variable RIPMD::missing] -rowspan 200
	place $w.n.f3.missing -x 28 -y 45
	grid [label $w.n.f3.phText -text "PDB2PQR pH"]
	place $w.n.f3.phText -x 35 -y 95
	grid  [entry $w.n.f3.pH -width 8 -justify right -textvariable pH -textvariable RIPMD::pH]
	place $w.n.f3.pH -x 300 -y 95
	
	grid [label $w.n.f3.otherText -text "Other options"]
	place $w.n.f3.otherText -x 0 -y 250	
	grid [label $w.n.f3.procText -text "Number of processors to use"]
	place $w.n.f3.procText -x 35 -y 275
	#grid [label $w.n.f3.procText2 -text "in program execution"]
	#place $w.n.f3.procText2 -x 35 -y 290
	grid  [entry $w.n.f3.proc -width 8 -justify right -textvariable nProc -textvariable RIPMD::proc]
	set nProc "1"
	place $w.n.f3.proc -x 300 -y 275
	
	grid [label $w.n.f3.fff -text "Path to CHARMM force field file"]
	place $w.n.f3.fff -x 35 -y 350
	
	pack [label $w.n.f3.fff2 -textvariable RIPMD::forceField]
	place $w.n.f3.fff2 -x 35 -y 375
	
	pack [button $w.n.f3.fffSelectButton -text "Browse..." -width 10 -command {
		set RIPMD::fffPath [tk_getOpenFile]
		if { $RIPMD::fffPath != "" } {
			set RIPMD::forceField $RIPMD::fffPath

		
		} else {
			set a [tk_dialog .myDialog "RIP-MD Error" "You need to provide a valid file" error 0 "Ok"]		
		}
	}]
	place $w.n.f3.fffSelectButton -x 300 -y 345


	grid [label $w.n.f3.pf -text "Path to CHARMM parameters file"]
	place $w.n.f3.pf -x 35 -y 425
	pack [label $w.n.f3.pf2 -textvariable RIPMD::parameterFile]
	place $w.n.f3.pf2 -x 35 -y 445
	
	pack [button $w.n.f3.pfSelectButton -text "Browse..." -width 10 -command {
		set RIPMD::pfPath [tk_getOpenFile]
		if { $RIPMD::pfPath != "" } {
			set RIPMD::parameterFile $RIPMD::pfPath

		
		} else {
			set a [tk_dialog .myDialog "RIP-MD Error" "You need to provide a valid file" error 0 "Ok"]		
		}
	}]
	place $w.n.f3.pfSelectButton -x 300 -y 420

	
	grid [label $w.n.f3.mdText -text "MD options"]
	place $w.n.f3.mdText -x 500 -y 20
	grid  [label $w.n.f3.percentText1 -text "Percentaje of time to consider"] -rowspan 50
	place $w.n.f3.percentText1 -x 528 -y 45
	grid  [label $w.n.f3.percentText2 -text "interactions in the consensus graph"] -rowspan 50
	place $w.n.f3.percentText2 -x 528 -y 60
	grid  [entry $w.n.f3.percentajeTime -width 8 -justify right -textvariable RIPMD::percentaje]
	place $w.n.f3.percentajeTime -x 793 -y 45
	grid  [label $w.n.f3.percentText3 -text "%"] -rowspan 50
	place $w.n.f3.percentText3 -x 860 -y 45
    
	grid  [label $w.n.f3.frameText1 -text "Frame to use as a representative"] -rowspan 50
	place $w.n.f3.frameText1 -x 528 -y 120
	grid  [label $w.n.f3.frameText2 -text "structure to display interactions"] -rowspan 50
	place $w.n.f3.frameText2 -x 528 -y 135
	grid  [entry $w.n.f3.frameToUse -width 8 -justify right -textvariable RIPMD::frame]
	place $w.n.f3.frameToUse -x 793 -y 120
    
	
	grid  [label $w.n.f3.win1 -text "Frame to start RIP-MD computation"] -rowspan 50
	place $w.n.f3.win1 -x 528 -y 180
	grid  [entry $w.n.f3.timespace -width 8 -justify right -textvariable RIPMD::frame_start]
	place $w.n.f3.timespace -x 793 -y 180

	grid  [label $w.n.f3.win2 -text "Frame to stop RIP-MD computation."] -rowspan 50
	place $w.n.f3.win2 -x 528 -y 240
	grid  [label $w.n.f3.winO1 -text "Negative numbers indicate frames to the end"] -rowspan 50
	place $w.n.f3.winO1 -x 528 -y 255
	grid  [label $w.n.f3.winO2 -text "of the simulation, being -1 the last frame"] -rowspan 50
	place $w.n.f3.winO2 -x 528 -y 270 	
	#grid  [label $w.n.f3.winO3 -text "Example -1 is the last frame)"] -rowspan 50
	#place $w.n.f3.winO3 -x 528 -y 285	
	grid  [entry $w.n.f3.timespace2 -width 8 -justify right -textvariable RIPMD::frame_end]
	place $w.n.f3.timespace2 -x 793 -y 240


	grid  [label $w.n.f3.win3 -text "Number of separation between frames"] -rowspan 50
	place $w.n.f3.win3 -x 528 -y 390
	grid  [entry $w.n.f3.timespace3 -width 8 -justify right -textvariable RIPMD::frame_separation]
	place $w.n.f3.timespace3 -x 793 -y 390

	###########################
	##
	## frame for button
	##
	###########################

	frame $w.buttonFrameCommand -height 30 -width 100
	#pack $w.buttonFrame
	pack [button $w.buttonFrameCommand.ripmdButtonCommand -text "Command" -command {
			#try to set RIPMD::pathToFiles with loaded molecules, if can not read it we will put an error message
			if {[catch {set RIPMD::pathToFiles [molinfo top get filename]} errmsg]} {
				set a [tk_dialog .myDialog "RIP-MD" "There are no molecules loaded" error 0 "Ok"]
				return
			}
			set RIPMD::errors "You present the following errors:\n\n" 
			#change directory
			#setting initial values
			set RIPMD::pdb "*"
			set RIPMD::dcd "*"
			set RIPMD::psf "*"
			set RIPMD::command "python $::env(RIP_MD)/main.py"	
					
			#files will have a array with one list, so files 2 will have the list and we will loop over this 
			#to recover filenames and type of file
			
			set files [lindex $RIPMD::pathToFiles 0]
			set files2 [lindex $files 0]	

			#now we will set  some counters to know how many files of each type of file we have loaded
			set countPDB 0
			set countPSF 0
			set countDCD 0
			
			
			
			for {set index 0} {$index < [llength $files]} {incr index } {
				set file [lindex $files $index]
				set splittedFile [split $file .]
				set format [lindex $splittedFile [ expr [llength $splittedFile] - 1 ] ]
				if { $format == "pdb" || $format == "PDB" } {
					set RIPMD::pdb $file
					incr countPDB
				} elseif { $format == "dcd" || $format == "DCD"} {
					set RIPMD::dcd $file
					incr countDCD
				} elseif { $format == "psf" || $format == "PSF"} {
					set RIPMD::psf $file
					incr countPSF
				} 
			}
			
			#now we will look for some conditions, like we can not compute networks if user
			#have a DCD and a PDB and others
			if { $RIPMD::pdb != "*" && $RIPMD::dcd != "*"} {
				append RIPMD::errors "* There are a PDB and a DCD file, please remove one.\n"
			} 
			if { $RIPMD::dcd == "*" && $RIPMD::pdb == "*" && $RIPMD::psf!= "*"} {
				append RIPMD::errors "* There is only a PSF file loaded, please load a PDB or a DCD file.\n"
			}
			if { $RIPMD::dcd != "*" && $RIPMD::psf == "*" } {
				append RIPMD::errors  "* There is only a DCD file loaded, please load a PSF file.\n"
			}
			if { $countPDB > 1 } {
				append RIPMD::errors  "* There are more than one PDB file loaded.\n"
			}						                                                               
			if { $countDCD > 1 } {                                                                 
				append RIPMD::errors  "* There are more than one DCD file loaded.\n"
			}                                                                                    
			if { $countPSF > 1 } {                                                              
				append RIPMD::errors  "* There are more than one PSF file loaded.\n"
			}
			
			#
			# last check, if numbers are really numbers
			#
			RIPMD::checkNumbers

			############################################################
			##
			## now we will format the command line to execute RIP-MD
			##
			############################################################
			
			set RIPMD::command "python $::env(RIP_MD)/main.py"
			
			#setting files
			if {$RIPMD::pdb != "*"} {
				append RIPMD::command " --pdb  $RIPMD::pdb"
				if {$RIPMD::psf!= "*"} {
					append RIPMD::command " --psf " $RIPMD::psf
				}				
			
			} elseif {$RIPMD::dcd !="*" && $RIPMD::psf != "*"} {
	
				append RIPMD::command " --dcd  $RIPMD::dcd --psf $RIPMD::psf --reference_frame $RIPMD::frame --separation_frame $RIPMD::frame_separation --frame_start $RIPMD::frame_start --frame_end $RIPMD::frame_end "
				
				#percentaje of time
				append RIPMD::command " --time " $RIPMD::percentaje
				
				#for pearson correlation
				if {$RIPMD::pearson == 1} {
					append RIPMD::command " --pearson_corr "
					append RIPMD::command " --plot_pearson "
				}
				
			} else {
				if { $RIPMD::dcd == "*" && $RIPMD::pdb == "*" && $RIPMD::psf!= "*"} {
					#it is for not repeat the error of only psf is charged
					set aux 1
				} else {
					set a [tk_dialog .myDialog "RIP-MD Error" "You are not using correct input files. Aborting" error 0 "Ok"]
					return
				}
			}
			
			#setting output folder
			if {$RIPMD::outputFolder != "*"} {
	
				append RIPMD::command " --output $RIPMD::outputFolder"
			} else {
				append RIPMD::errors "* You must select an output folder.\n"
			}
			
			#setting selection to analyze
			append RIPMD::command " --selection " $RIPMD::proteinSelection

			#setting graph output format
			append RIPMD::command " --gformat " $RIPMD::format
			
			#setting number of processors and pH
			
			append RIPMD::command " --nproc " $RIPMD::proc
			append RIPMD::command " --pH " $RIPMD::pH
			if {$RIPMD::missing == 1} {
				append RIPMD::command " --missing_atoms "
			}
			
			#setting for C alpha
			if {$RIPMD::calpha == 1} {
				append RIPMD::command " --calpha"
				append RIPMD::command " --ca_dist " $RIPMD::calphaDist		
				
			}
			
			#settings for hbond
			if {$RIPMD::hbonds == 1} {
				append RIPMD::command " --hbond"
				append RIPMD::command " --h_dist " $RIPMD::hbondDist
				append RIPMD::command " --h_angle " $RIPMD::hbondAngle
			}
			
			#settings for salt bridges
			if {$RIPMD::salt ==1} {
				append RIPMD::command " --salt"
				append RIPMD::command " --s_distance " $RIPMD::saltDist
			}
			
			#settings for disulphide bridges
			if {$RIPMD::disulphide==1} {
				append RIPMD::command " --disulfide"
				append RIPMD::command " --ss_distance " $RIPMD::disulphideDist
				append RIPMD::command " --ss_angle \[" $RIPMD::disulphideAng1 "," $RIPMD::disulphideAng2 "\]"
			}
			
			#settings for cation pi interactions
			if {$RIPMD::capi == 1} {
				append RIPMD::command " --cation_pi"
				append RIPMD::command " --cation_pi_distance " $RIPMD::capiDist
				append RIPMD::command " --cp_angle1 \[" $RIPMD::capiAngle1 "," $RIPMD::capiAngle2 "," $RIPMD::capiAngle3 "," $RIPMD::capiAngle4 "\]"
			}
			
			#settings for pi pi interactions
			if {$RIPMD::pipi == 1} {
				append RIPMD::command " --pi_pi"
				append RIPMD::command " --pi_pi_distance " $RIPMD::pipiDist
			}
	
			#settings for arg arg interactions
			if {$RIPMD::ArgArg == 1} {
				append RIPMD::command " --arg_arg"
				append RIPMD::command " --arg_arg_distance " $RIPMD::ArgArgDist
			}			
			
			#settings for vdw contact
			if {$RIPMD::vdw == 1} {
				append RIPMD::command " --vdw"
				append RIPMD::command " --vdw_excluded " $RIPMD::covalentDist
				append RIPMD::command " --vdw_range \[" $RIPMD::distVDW1 "," $RIPMD::distVDW2 "\]"
			}
			
			#settings for coulomb
			if {$RIPMD::coulomb == 1} {
				append RIPMD::command " --coulomb"
				if {$RIPMD::reactionField == 1} {
					append RIPMD::command " --reaction_field "
				}
				append RIPMD::command " --simulated_permittivity " $RIPMD::simulatedPermittivity
				append RIPMD::command " --RF_permittivity " $RIPMD::reactionFieldPermittivity
				append RIPMD::command " --temperature " $RIPMD::temp
				append RIPMD::command " --inverse_debye_screening " $RIPMD::inverseDebye
				append RIPMD::command " --distance_threshold " $RIPMD::coulombThreshold
				append RIPMD::command " --coulomb_excluded " $RIPMD::coulombCovalentDist				
				
			}
							
			#appending charmm force field file and parameter file
			append RIPMD::command " --force_field  $RIPMD::forceField  --parameter_file  $RIPMD::parameterFile"

			############################################################
			##
			## if any error was detected we will show errors and then 
			## block the execution of RIP-MD
			##	
			############################################################
			
			if { $RIPMD::errors != "You present the following errors:\n\n" } {
				set a [tk_dialog .myDialog "RIP-MD Error" "$RIPMD::errors\nPlease fix error/s before continueing (for more information, please read the help section)" error 0 "Ok"]
				return
			}
			
			
			############################################################
			##
			## saving command
			##
			############################################################
			set file [tk_getSaveFile -title "Saving command" -parent .]
			if { $file == "" } {
				return; # they clicked cancel
				}
			set x [catch {set fid [open $file w+]}]
			set y [catch {puts $fid "$RIPMD::command"}]
			set z [catch {close $fid}]
			#puts $y
			if { $x || $y || $z || ![file exists $file] || ![file isfile $file] || ![file readable $file] } {
			tk_messageBox -parent . -icon error \
							-message "An error occurred while saving to \"$file\""
				} else {
			tk_messageBox -parent . -icon info \
							-message "Save successful, do not forget change the necessary paths"
				}			
		}]
	place $w.buttonFrameCommand -x 860 -y 560


###############
	frame $w.buttonFrame -height 30 -width 100
	#pack $w.buttonFrame
	pack [button $w.buttonFrame.ripmdButton -text "Run" -command {
			#try to set RIPMD::pathToFiles with loaded molecules, if can not read it we will put an error message
			if {[catch {set RIPMD::pathToFiles [molinfo top get filename]} errmsg]} {
				set a [tk_dialog .myDialog "RIP-MD" "There are no molecules loaded" error 0 "Ok"]
				return
			}
			set RIPMD::errors "You present the following errors:\n\n" 
			#change directory
			#cd "$env(VMDDIR)/plugins/noarch/tcl/RIP-MD/python"
			#pwd
			#setting initial values
			set RIPMD::pdb "*"
			set RIPMD::dcd "*"
			set RIPMD::psf "*"
			#set RIPMD::command "python ./python/main.py "
			set RIPMD::command "python $::env(RIP_MD)/main.py "	
					
			#files will have a array with one list, so files 2 will have the list and we will loop over this 
			#to recover filenames and type of file
			
			set files [lindex $RIPMD::pathToFiles 0]
			set files2 [lindex $files 0]	

			#now we will set  some counters to know how many files of each type of file we have loaded
			set countPDB 0
			set countPSF 0
			set countDCD 0
			
			
			
			for {set index 0} {$index < [llength $files]} {incr index } {
				set file [lindex $files $index]
				set splittedFile [split $file .]
				set format [lindex $splittedFile [ expr [llength $splittedFile] - 1 ] ]
				if { $format == "pdb" || $format == "PDB" } {
					set RIPMD::pdb $file
					incr countPDB
				} elseif { $format == "dcd" || $format == "DCD"} {
					set RIPMD::dcd $file
					incr countDCD
				} elseif { $format == "psf" || $format == "PSF"} {
					set RIPMD::psf $file
					incr countPSF
				} 
			}
			
			#now we will look for some conditions, like we can not compute networks if user
			#have a DCD and a PDB and others
			if { $RIPMD::pdb != "*" && $RIPMD::dcd != "*"} {
				append RIPMD::errors "* There are a PDB and a DCD file.\n"
			} 
			if { $RIPMD::dcd == "*" && $RIPMD::pdb == "*" && $RIPMD::psf!= "*"} {
				append RIPMD::errors "* There is only a PSF file loaded.\n"
			}
			if { $RIPMD::dcd != "*" && $RIPMD::psf == "*" } {
				append RIPMD::errors  "* There is only a DCD file loaded.\n"
			}
			if { $countPDB > 1 } {
				append RIPMD::errors  "* There are more than one PDB file loaded.\n"
			}						                                                               
			if { $countDCD > 1 } {                                                                 
				append RIPMD::errors  "* There are more than one DCD file loaded.\n"
			}                                                                                    
			if { $countPSF > 1 } {                                                              
				append RIPMD::errors  "* There are more than one PSF file loaded.\n"
			}
			
			#
			# last check, if numbers are really numbers
			#
			RIPMD::checkNumbers

			############################################################
			##
			## now we will format the command line to execute RIPMD
			##
			############################################################
			
			set RIPMD::command "python $::env(RIP_MD)/main.py"
			
			#setting files
			if {$RIPMD::pdb != "*"} {	
				append RIPMD::command " --pdb " $RIPMD::pdb
				if {$RIPMD::psf!= "*"} {
					append RIPMD::command " --psf " $RIPMD::psf
				}
			
			} elseif {$RIPMD::dcd !="*" && $RIPMD::psf != "*"} {
				set answer [tk_dialog .dialog1 "RIP-MD" "Using RIP-MD when a MD is charged on VMD may cause that your pc turn freeze\n\nDo you want to continue?" question 1 "Yes" "No"]
				if {$answer == 1} {
					return
				}
							
				append RIPMD::command " --dcd  $RIPMD::dcd --psf $RIPMD::psf --reference_frame $RIPMD::frame --separation_frame $RIPMD::frame_separation --frame_start $RIPMD::frame_start --frame_end $RIPMD::frame_end "
				
				#percentaje of time
				append RIPMD::command " --time " $RIPMD::percentaje
				
				#for pearson correlation
				if {$RIPMD::pearson == 1} {
					append RIPMD::command " --pearson_corr "
					append RIPMD::command " --plot_pearson "
				}
				
			} else {
				if { $RIPMD::dcd == "*" && $RIPMD::pdb == "*" && $RIPMD::psf!= "*"} {
					#it is for not repeat the error of only psf is charged
					set aux 1
				} else {
					set a [tk_dialog .myDialog "RIP-MD Error" "Unexpected Error. Aborting" error 0 "Ok"]
					return
				}
			}
			
			#setting output folder
			if {$RIPMD::outputFolder != "*"} {
			
				append RIPMD::command " --output $RIPMD::outputFolder"
			} else {
				append RIPMD::errors "* You must select an output folder.\n"
			}

			#setting selection to analyze
			append RIPMD::command " --selection " $RIPMD::proteinSelection
			
			#setting graph output format
			append RIPMD::command " --gformat " $RIPMD::format
			
			#setting number of processors and pH
			
			append RIPMD::command " --nproc " $RIPMD::proc
			append RIPMD::command " --pH " $RIPMD::pH
			if {$RIPMD::missing == 1} {
				append RIPMD::command " --missing_atoms "
			}			
			#setting for C alpha
			if {$RIPMD::calpha == 1} {
				append RIPMD::command " --calpha"
				append RIPMD::command " --ca_dist " $RIPMD::calphaDist		
				
			}
			
			#settings for hbond
			if {$RIPMD::hbonds == 1} {
				append RIPMD::command " --hbond"
				append RIPMD::command " --h_dist " $RIPMD::hbondDist
				append RIPMD::command " --h_angle " $RIPMD::hbondAngle
			}
			
			#settings for salt bridges
			if {$RIPMD::salt ==1} {
				append RIPMD::command " --salt"
				append RIPMD::command " --s_distance " $RIPMD::saltDist
			}
			
			#settings for disulphide bridges
			if {$RIPMD::disulphide==1} {
				append RIPMD::command " --disulfide"
				append RIPMD::command " --ss_distance " $RIPMD::disulphideDist
				append RIPMD::command " --ss_angle \[" $RIPMD::disulphideAng1 "," $RIPMD::disulphideAng2 "\]"
			}
			
			#settings for cation pi interactions
			if {$RIPMD::capi == 1} {
				append RIPMD::command " --cation_pi"
				append RIPMD::command " --cation_pi_distance " $RIPMD::capiDist
				append RIPMD::command " --cp_angle1 \[" $RIPMD::capiAngle1 "," $RIPMD::capiAngle2 "," $RIPMD::capiAngle3 "," $RIPMD::capiAngle4 "\]"
			}
			
			#settings for pi pi interactions
			if {$RIPMD::pipi == 1} {
				append RIPMD::command " --pi_pi"
				append RIPMD::command " --pi_pi_distance " $RIPMD::pipiDist
			}
	
			#settings for arg arg interactions
			if {$RIPMD::ArgArg == 1} {
				append RIPMD::command " --arg_arg"
				append RIPMD::command " --arg_arg_distance " $RIPMD::ArgArgDist
			}			
						
			#settings for vdw contact
			if {$RIPMD::vdw == 1} {
				append RIPMD::command " --vdw"
				append RIPMD::command " --vdw_excluded " $RIPMD::covalentDist
				append RIPMD::command " --vdw_range \[" $RIPMD::distVDW1 "," $RIPMD::distVDW2 "\]"
			}
			
			#settings for coulomb
			if {$RIPMD::coulomb == 1} {
				append RIPMD::command " --coulomb"
				
				if {$RIPMD::reactionField == 1} {
					append RIPMD::command " --reaction_field "
				}
				append RIPMD::command " --simulated_permittivity " $RIPMD::simulatedPermittivity
				append RIPMD::command " --RF_permittivity " $RIPMD::reactionFieldPermittivity
				append RIPMD::command " --temperature " $RIPMD::temp
				append RIPMD::command " --inverse_debye_screening " $RIPMD::inverseDebye
				append RIPMD::command " --distance_threshold " $RIPMD::coulombThreshold
				append RIPMD::command " --coulomb_excluded " $RIPMD::coulombCovalentDist
				
			}
			#for charmm topology file and charmm parameter file
					
			append RIPMD::command " --force_field  $RIPMD::forceField  --parameter_file  $RIPMD::parameterFile"
			
			
			
			
			############################################################
			##
			## if any error was detected we will show errors and then 
			## block the execution of RIP-MD
			##	
			############################################################
			
			if { $RIPMD::errors != "You present the following errors:\n\n" } {
				set a [tk_dialog .myDialog "RIP-MD Error" "$RIPMD::errors\nPlease fix error/s before continueing (for more information, please read the help section)" error 0 "Ok"]
				return
			}
			
			
			############################################################
			##
			## executing RIP-MD !!!
			##
			############################################################
			open |$RIPMD::command
	
			
			set end 0
			#we will look if the output log exist
			while { $end == 0} {
				after 1000
				if { [file exists "$RIPMD::outputFolder/output.log"] == 1 } {
					set end 1
				}
			}		
			set cont 0	
			#when file exists we will read it and puts the new lines
			while { $end == 1 } {
				after 5000
				set aux 0
				set fp [open "$RIPMD::outputFolder/output.log" r]
				set file_data [read $fp]
				close $fp
				
				set data [split $file_data "\n"]
				
				foreach line $data {
					
					if { $aux >= $cont } {
						#puts to print in console the log status
						puts "[lindex $data $aux]"
						
						if { [string match *error* $line] == 1 || [string match *Error* $line] == 1 || [string match *ERROR* $line] == 1 } {
							set end 3
						}
						if { [lindex $data $aux] == "END" } {
							set end 2							
						}
					}
					incr aux
					
				}
				set cont $aux
			}
			
			if { $end == 3 } {
				$RIPMD::w.n.f1.textBox.l insert end ""	
				$RIPMD::w.n.f1.textBox.l insert end "RIP-MD has been finished with errors, for more information please see the log file"	
				
			}
			if { $end == 2 } {
				$RIPMD::w.n.f1.textBox.l insert end ""	
				$RIPMD::w.n.f1.textBox.l insert end "RIP-MD has been finished without errors"
				$RIPMD::w.n.f1.textBox.l insert end "Results files are in $RIPMD::outputFolder"					
				$RIPMD::w.n.f1.textBox.l insert end "for more information, see the README file"
				$RIPMD::w.n.f1.textBox.l insert end "Loading Nodes..."	
				$RIPMD::w.n.f1.textBox.l insert end ""	
				set RIPMD::resultFolder "$RIPMD::outputFolder"
				RIPMD::loadResults
			}
		}]
	
	####################################################################
	##
	## Result Tab
	##
	####################################################################
	$RIPMD::w.n add $RIPMD::w.n.f5 -text "Results"
	#putting a label
	grid  [label $RIPMD::w.n.f5.resLab -text "Network nodes"] -rowspan 50
	place $RIPMD::w.n.f5.resLab -x 10 -y 10
	
	frame $RIPMD::w.n.f5.nodesBox
	grid [listbox $RIPMD::w.n.f5.nodesBox.l -selectmode extended -yscrollcommand "$RIPMD::w.n.f5.nodesBox.s1 set" -xscrollcommand "$RIPMD::w.n.f5.nodesBox.s2 set" -height 25 -width 20 ] -column 0 -row 0 -sticky nwes
	grid [ttk::scrollbar $RIPMD::w.n.f5.nodesBox.s1 -command "$RIPMD::w.n.f5.nodesBox.l yview" -orient vertical   ] -column 1 -row 0 -sticky ns
	grid [ttk::scrollbar $RIPMD::w.n.f5.nodesBox.s2 -command "$RIPMD::w.n.f5.nodesBox.l xview" -orient horizontal ] -column 0 -row 1 -sticky ew
	pack $RIPMD::w.n.f5.nodesBox
	place $RIPMD::w.n.f5.nodesBox -x 55 -y 40	
	
	#event to know nodes selected
	bind $RIPMD::w.n.f5.nodesBox.l <<ListboxSelect>> { 
		RIPMD::selectionMade %W
	}
	#interactions to display
	grid  [label $RIPMD::w.n.f5.disLab -text "Interactions to Display"] -rowspan 50
	place $RIPMD::w.n.f5.disLab -x 280 -y 10	
	grid  [checkbutton $RIPMD::w.n.f5.calpha -text "C alpha contacts (Blue)" -variable RIPMD::displayCalpha -state disable -command { RIPMD::drawBonds }] -rowspan 200
	place $RIPMD::w.n.f5.calpha -x 300 -y 40
	grid [checkbutton $RIPMD::w.n.f5.hbond -text "H Bonds (White)" -variable RIPMD::displayHbonds -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.hbond -x 300 -y 60
	grid [checkbutton $RIPMD::w.n.f5.salt -text "Salt Bridges (Orange)" -variable RIPMD::displaySalt -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.salt -x 300 -y 80 ;#I dont know why but this put in a correctly position all checkbuttons
	grid [checkbutton $RIPMD::w.n.f5.disulphide -text "Disulfide Bridges (Yellow)" -variable RIPMD::displayDisulphide -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.disulphide -x 300 -y 100 
	grid [checkbutton $RIPMD::w.n.f5.capi -text "Cation - pi interactions (Green)" -variable RIPMD::displayCapi -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.capi -x 300 -y 120
	#\u3c0
	grid [checkbutton $RIPMD::w.n.f5.pipi -text "pi - pi interactions (Pink)" -variable RIPMD::displayPipi -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.pipi -x 300 -y 140
	grid [checkbutton $RIPMD::w.n.f5.argarg -text "Arg - Arg interactions (Red)" -variable RIPMD::displayArgArg -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.argarg -x 300 -y 160	
	grid [checkbutton $RIPMD::w.n.f5.vdw -text "VdW contacts (Purple)" -variable RIPMD::displayVdw -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.vdw -x 300 -y 180
	grid [checkbutton $RIPMD::w.n.f5.coulomb -text "Coulomb interactions (Tan)" -variable RIPMD::displayCoulomb -state disable -command { RIPMD::drawBonds }] -rowspan 50
	place $RIPMD::w.n.f5.coulomb -x 300 -y 200		
	
	#button to select the pearson correlation
	
	grid [label $w.n.f5.selectText -text "Select interactions to display correlation map"]
	place $w.n.f5.selectText -x 600 -y 10
	
	ttk::menubutton $w.n.f5.pearson1 -menu $w.n.f5.pearson1.type  -textvariable RIPMD::pearson1 -width 10
	menu $w.n.f5.pearson1.type
	$w.n.f5.pearson1.type add command -label "All" -command { set RIPMD::pearson1 "All" }
	$w.n.f5.pearson1.type add command -label "C Alpha" -command { set RIPMD::pearson1 "C Alpha" }
	$w.n.f5.pearson1.type add command -label "H Bonds" -command { set RIPMD::pearson1 "H Bonds" }
	$w.n.f5.pearson1.type add command -label "Salt Bridges" -command { set RIPMD::pearson1 "Salt Bridges" }	
	$w.n.f5.pearson1.type add command -label "Cation - pi interaction" -command { set RIPMD::pearson1 "Cation - pi interaction" }
	$w.n.f5.pearson1.type add command -label "pi - pi interaction" -command { set RIPMD::pearson1 "pi - pi interaction" }
	$w.n.f5.pearson1.type add command -label "Arg - Arg interaction" -command { set RIPMD::pearson1 "Arg - Arg interaction" }
	$w.n.f5.pearson1.type add command -label "Coulomb interaction" -command { set RIPMD::pearson1 "Coulomb interaction" }
	$w.n.f5.pearson1.type add command -label "VdW contacts" -command { set RIPMD::pearson1 "VdW contacts" }
	
	pack $w.n.f5.pearson1
	place $w.n.f5.pearson1 -x 620 -y 50


	ttk::menubutton $w.n.f5.pearson2 -menu $w.n.f5.pearson2.type  -textvariable RIPMD::pearson2 -width 10
	menu $w.n.f5.pearson2.type
	$w.n.f5.pearson2.type add command -label "All" -command { set RIPMD::pearson2 "All" }
	$w.n.f5.pearson2.type add command -label "C Alpha" -command { set RIPMD::pearson2 "C Alpha" }
	$w.n.f5.pearson2.type add command -label "H Bonds" -command { set RIPMD::pearson2 "H Bonds" }
	$w.n.f5.pearson2.type add command -label "Salt Bridges" -command { set RIPMD::format "Salt Bridges" }	
	$w.n.f5.pearson2.type add command -label "Cation - pi interaction" -command { set RIPMD::pearson2 "Cation - pi interaction" }
	$w.n.f5.pearson2.type add command -label "pi - pi interaction" -command { set RIPMD::pearson2 "pi - pi interaction" }
	$w.n.f5.pearson2.type add command -label "Arg - Arg interaction" -command { set RIPMD::pearson2 "Arg - Arg interaction" }
	$w.n.f5.pearson2.type add command -label "Coulomb interaction" -command { set RIPMD::pearson2 "Coulomb interaction" }
	$w.n.f5.pearson2.type add command -label "VdW contacts" -command { set RIPMD::pearson2 "VdW contacts" }
	
	pack $w.n.f5.pearson2
	place $w.n.f5.pearson2 -x 800 -y 50
	
	pack [button $w.n.f5.showPearson -text "Display" -width 10 -state disable -command {
		
		if { [catch { exec python $::env(RIP_MD)/libs/plotCorr.py $RIPMD::pearson1 $RIPMD::pearson2 $RIPMD::resultFolder } msg] } {
			set a [tk_dialog .myDialog "RIP-MD Error" "Failed to open Pearson Plot. If you compute these interaction pair (and did you select the compute pearson correlation option), please execute this command in a terminal\n\npython $::env(RIP_MD)/libs/plotCorr.py $RIPMD::pearson1 $RIPMD::pearson2 $RIPMD::resultFolder" error 0 "Ok"]
			return
		}

    
	}]
	place $w.n.f5.showPearson -x 823 -y 100


	grid [label $w.n.f5.histogramText -text "Select interaction to compute interaction presence Histograms"]
	place $w.n.f5.histogramText -x 600 -y 200
	
	ttk::menubutton $w.n.f5.histogram -menu $w.n.f5.histogram.type  -textvariable RIPMD::hist -width 41
	menu $w.n.f5.histogram.type
	$w.n.f5.histogram.type add command -label "All" -command { set RIPMD::hist "All interactions" }
	$w.n.f5.histogram.type add command -label "C Alpha" -command { set RIPMD::hist "C Alpha" }
	$w.n.f5.histogram.type add command -label "H Bonds" -command { set RIPMD::hist "H Bonds" }
	$w.n.f5.histogram.type add command -label "Salt Bridges" -command { set RIPMD::hist "Salt Bridges" }	
	$w.n.f5.histogram.type add command -label "Disulfide Bridges" -command { set RIPMD::hist "Disulfide Bridges" }	
	$w.n.f5.histogram.type add command -label "Cation - pi interaction" -command { set RIPMD::hist "Cation - pi interaction" }
	$w.n.f5.histogram.type add command -label "pi - pi interaction" -command { set RIPMD::hist "pi - pi interaction" }
	$w.n.f5.histogram.type add command -label "Arg - Arg interaction" -command { set RIPMD::hist "Arg - Arg interaction" }
	$w.n.f5.histogram.type add command -label "Coulomb interaction" -command { set RIPMD::hist "Coulomb interaction" }
	$w.n.f5.histogram.type add command -label "VdW contacts" -command { set RIPMD::hist "VdW contacts" }
	
	pack $w.n.f5.histogram
	place $w.n.f5.histogram -x 600 -y 240
	
################################################################################################

	pack [button $w.n.f5.showHist -text "Display" -width 10 -state disable -command {
		
		set a [tk_dialog .myDialog "RIP-MD" "Computing Histograms, this could take several time" warning 0 "Ok"]
		if { [catch { exec python $::env(RIP_MD)/libs/plotHist.py $RIPMD::hist $RIPMD::resultFolder } msg] } {
			set a [tk_dialog .myDialog "RIP-MD Error" "Failed to open histogram Plot. Please execute this command in a terminal\n\npython $::env(RIP_MD)/libs/plotHist.py $RIPMD::pearson1 $RIPMD::pearson2 $RIPMD::resultFolder" error 0 "Ok"]
			return
		}    
	}]
	place $w.n.f5.showHist -x 840 -y 290

################################################################################################	
	
	##### packing some graphical packages		
	place $w.buttonFrame -x 950 -y 560

	pack $w.menubar.file -side left	
	pack $w.menubar.help -side left
	pack $w.n

}

proc ripmd_tk {} {
	RIPMD::ripmd
	return $RIPMD::w
}
