wipe
# ============================================================================
# Cores Information
# ============================================================================
set regcores       5
set pmlcores       3
set drmcores       1
set structurecores 1
set AnalysisType   PML
# ============================================================================
# Model Information
# ============================================================================
set StructureType     "STEEL"
set NStory            10
set NBay              4
set NBayZ             4
set StartNodeX       -9.0
set StartNodeY       -9.0
set StartNodeZ       -1.5
set LCol              3.0
set LBeam             4.5
set LGird             4.5
set SectionType       Elastic
set HaveStructure     "YES"
# ============================================================================
# Soil Information
# ============================================================================
set dx                2.5
set dy                2.5
set dz                2.5
set llx               60.0
set lly               60.0
set llz               40.0
set drmthicknessx     2.5
set drmthicknessy     2.5
set drmthicknessz     2.5
set numdrmlayers      2
set lx                50.0
set ly                50.0
set lz                35.0
set nx                20.0
set ny                20.0
set nz                14.0
# ============================================================================
# PML information
# ============================================================================
set AbsorbingElements "PML"
set numpmllayers      2
set pmlthicknessx     2.5
set pmlthicknessy     2.5
set pmlthicknessz     2.5
set pmltotalthickness 5.0
set HaveAbsorbingElements "NO"
# ============================================================================
# General information
# ============================================================================
set meshdir           "OpenSeesmesh/PML"
set outputdir         "Results/PML"
set DRMFile           "/home/amnp95/Projects/OpenSeesProjects/Structure/DRMload/SurfaceWave.h5drm"
# ============================================================================
# Embedding foundation
# ============================================================================
set EmbeddingFoundation "YES"
set EmbeddedFoundation [dict create xmax 10.0 xmin -10.0 ymax 10.0 ymin -10.0 zmax 0.0 zmin -5.0]
# ============================================================================
# Fondation information
# ============================================================================
set HaveFoundation "YES"
set foundationBlocks {}
lappend foundationBlocks [dict create matTag 1 xmax 10.0 xmin -10.0 ymax 10.0 ymin -10.0 zmax -1.5 zmin -4.5 Xmeshsize 1.0 Ymeshsize 1.0 Zmeshsize 1.0]
# ============================================================================
# Piles information
# ============================================================================
set HavePiles "YES"
set pilelist {}
lappend pilelist [dict create xtop -7.0 ytop -7.0 ztop -3.0 xbottom -7.0 ybottom -7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -3.5 ytop -7.0 ztop -3.0 xbottom -3.5 ybottom -7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 0.0 ytop -7.0 ztop -3.0 xbottom 0.0 ybottom -7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 3.5 ytop -7.0 ztop -3.0 xbottom 3.5 ybottom -7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 7.0 ytop -7.0 ztop -3.0 xbottom 7.0 ybottom -7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -7.0 ytop -3.5 ztop -3.0 xbottom -7.0 ybottom -3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -3.5 ytop -3.5 ztop -3.0 xbottom -3.5 ybottom -3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 0.0 ytop -3.5 ztop -3.0 xbottom 0.0 ybottom -3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 3.5 ytop -3.5 ztop -3.0 xbottom 3.5 ybottom -3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 7.0 ytop -3.5 ztop -3.0 xbottom 7.0 ybottom -3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -7.0 ytop 0.0 ztop -3.0 xbottom -7.0 ybottom 0.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -3.5 ytop 0.0 ztop -3.0 xbottom -3.5 ybottom 0.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 0.0 ytop 0.0 ztop -3.0 xbottom 0.0 ybottom 0.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 3.5 ytop 0.0 ztop -3.0 xbottom 3.5 ybottom 0.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 7.0 ytop 0.0 ztop -3.0 xbottom 7.0 ybottom 0.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -7.0 ytop 3.5 ztop -3.0 xbottom -7.0 ybottom 3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -3.5 ytop 3.5 ztop -3.0 xbottom -3.5 ybottom 3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 0.0 ytop 3.5 ztop -3.0 xbottom 0.0 ybottom 3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 3.5 ytop 3.5 ztop -3.0 xbottom 3.5 ybottom 3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 7.0 ytop 3.5 ztop -3.0 xbottom 7.0 ybottom 3.5 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -7.0 ytop 7.0 ztop -3.0 xbottom -7.0 ybottom 7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop -3.5 ytop 7.0 ztop -3.0 xbottom -3.5 ybottom 7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 0.0 ytop 7.0 ztop -3.0 xbottom 0.0 ybottom 7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 3.5 ytop 7.0 ztop -3.0 xbottom 3.5 ybottom 7.0 zbottom -10.0 numberofElements 6]
lappend pilelist [dict create xtop 7.0 ytop 7.0 ztop -3.0 xbottom 7.0 ybottom 7.0 zbottom -10.0 numberofElements 6]
# ============================================================================
# cells and nodes information
# ============================================================================
set soilfoundation_num_cells 15184
set soilfoundation_num_nodes 19806
# ===================================================================== #
# 3D test model for the pml element modeling the plane strain field     #
# University of Washington, Department of Civil and Environmental Eng   #
# Geotechnical Eng Group, A. Pakzad, P. Arduino - Jun 2023              #
# Basic units are m, Ton(metric), s										#
# ===================================================================== #

# ============================================================================
# define helper functions
# ============================================================================
proc writeNodesinFile {filename nodes} {
    # check if the file exists
    if {[file exists $filename] == 1} { file delete $filename }
    set f [open $filename "w"]
    foreach node $nodes {
        puts $f "$node [nodeCoord $node]"
    }
    close $f
}

proc writeElesinFile {filename eles} {
    # check if the file exists
    if {[file exists $filename] == 1} { file delete $filename }
    set f [open $filename "w"]
    foreach ele $eles {
        puts $f "$ele [eleNodes $ele]"
    }
    close $f
}


proc getMaxNodeTag {np pid} {
    set maxNodeTag -1
    foreach node [getNodeTags] {if {$node > $maxNodeTag} {set maxNodeTag $node}}
    barrier
    # puts "maxNodeTag: $maxNodeTag pid: $pid np: $np"
    if {$pid == 0} {set maxNodeTaglist {}}
    if {$pid == 0} {
        for {set i 0} {$i < $np} {incr i} {lappend maxNodeTaglist -1}
        # set the first element of the list to the maxNodeTag
        set maxNodeTaglist [lreplace $maxNodeTaglist 0 0 $maxNodeTag]
    }
    barrier
    if {$pid > 0} {
        send -pid 0 $maxNodeTag
    } else {
        for {set i 1} {$i < $np} {incr i} {
            recv -pid $i maxNodeTag
            set maxNodeTaglist [lreplace $maxNodeTaglist $i $i $maxNodeTag]
        }
    }
    barrier
    if {$pid == 0} {
        set maxNodeTag -1
        foreach node $maxNodeTaglist {if {$node > $maxNodeTag} {set maxNodeTag $node}}
        # puts "maximum: $maxNodeTag"
    }

    barrier
    # return the maxNodeTag
    if {$pid == 0} {
        for {set i 1} {$i < $np} {incr i} {
            send -pid $i $maxNodeTag
        }
    } else {
        recv -pid 0 maxNodeTag
    }
    barrier
    if {$maxNodeTag == -1} {set maxNodeTag 0}
    return $maxNodeTag
}


proc getMaxEleTag {np pid} {
    set maxEleTag -1
    foreach ele [getEleTags] {if {$ele > $maxEleTag} {set maxEleTag $ele}}
    barrier
    # puts "maxEleTag: $maxEleTag pid: $pid np: $np"
    if {$pid == 0} {set maxEleTaglist {}}
    if {$pid == 0} {
        for {set i 0} {$i < $np} {incr i} {lappend maxEleTaglist -1}
        # set the first element of the list to the maxEleTag
        set maxEleTaglist [lreplace $maxEleTaglist 0 0 $maxEleTag]
    }
    barrier
    if {$pid > 0} {
        send -pid 0 $maxEleTag
    } else {
        for {set i 1} {$i < $np} {incr i} {
            recv -pid $i maxEleTag
            set maxEleTaglist [lreplace $maxEleTaglist $i $i $maxEleTag]
        }
    }
    barrier
    if {$pid == 0} {
        set maxEleTag -1
        foreach ele $maxEleTaglist {if {$ele > $maxEleTag} {set maxEleTag $ele}}
        # puts "maximum: $maxEleTag"
    }

    barrier
    # return the maxEleTag
    if {$pid == 0} {
        for {set i 1} {$i < $np} {incr i} {
            send -pid $i $maxEleTag
        }
    } else {
        recv -pid 0 maxEleTag
    }
    barrier
    if {$maxEleTag == -1} {set maxEleTag 0}
    return $maxEleTag
}

proc addVartoModelInfoFile {fileName varName varValue pid writerpid} {
    if {$pid != $writerpid} {return}
    set f [open "$fileName" "r"]
    set lines [split [read $f] "\n"]
    close $f
    set f [open "$fileName" "w+"]
    set j 1
    foreach line $lines {
        set nextline [lindex $lines $j]
        if {$line == "\}"} {
            puts $f "\t\"$varName\": \"$varValue\""
            puts $f "\}"
            break
        } 
        if {$nextline == "\}" && $j > 1} {
            puts $f "$line,"
        } else {
            puts $f "$line"
        }
        incr j
    }
    close $f
}
     

proc initializModelinfoFile {fileName pid writerpid } {
    if {$pid != $writerpid} {return}
    if {[file exists $fileName] == 1} { file delete $fileName }
    set f [open "$fileName" "w+"]
    puts $f "\{"
    puts $f "\}"
    close $f
}



# ============================================================================
#  get the number of processors
# ============================================================================

set pid [getPID]
set np  [getNP]
set ModelInfoFile "$outputdir/ModelInfo.json"
initializModelinfoFile $ModelInfoFile $pid 0
addVartoModelInfoFile $ModelInfoFile "numberOfProcessors" $np $pid 0



# ============================================================================
# create the structure model
# ============================================================================
#  structure model will be built on the pid=0 processor all 
#  other processors will be used for the for fundation, soil, DRM, and PML elements

if {$pid < $structurecores} {
    if {$StructureType == "STEEL"} {
        source StructresFiles/SteelStructure.tcl
    }
    if {$StructureType == "CONCRETE"} {
        source ConcreteStructure.tcl
    }
    if {$StructureType == "Custom"} {
        source CustomStructure.tcl
    }
    puts "Structure model is created"
}
barrier

# ============================================================================
# update the maxNodeTag and maxEleTag
# ============================================================================
set maxNodeTag [getMaxNodeTag $np $pid]
set maxEleTag  [getMaxEleTag $np $pid]
barrier

# ============================================================================
#  Setting the maxNodeTag and maxEleTag for the structure model
# ============================================================================
# this is really important to set the maxNodeTag and maxEleTag for the structure model
set StructureMaxNodeTag $maxNodeTag
set StructureMaxEleTag  $maxEleTag
barrier
addVartoModelInfoFile $ModelInfoFile "StructureMaxNodeTag" $StructureMaxNodeTag $pid 0
addVartoModelInfoFile $ModelInfoFile  "StructureMaxEleTag"  $StructureMaxEleTag $pid 0


# ============================================================================
# bulding regular elements
# ============================================================================
model BasicBuilder -ndm 3 -ndf 3
set Vs          200.0
set nu          0.25         ;# --- Poisson's Ratio
set rho         2300.0                 ;# --- Density KG/M^3
set G           [expr $Vs*$Vs*$rho]
set E           [expr 2*$G*(1+$nu)]
nDMaterial ElasticIsotropic 1 $E $nu $rho;
nDMaterial ElasticIsotropic 2 [expr $E *100] $nu [expr $rho*4.0];
set SoilmatTag1 "1 0.0 0.0 -9.81";
set FoundationmatTag1 "2 0.0 0.0 -9.81";


# ============================================================================
# create regular nodes and elements
# ============================================================================
if {$HaveStructure == "YES" } {
    puts "StructureMaxNodeTag : $StructureMaxNodeTag\n"
    if {$pid>= $structurecores && $pid < [expr $regcores+$structurecores]} {
        model BasicBuilder -ndm 3 -ndf 3
        eval "source $meshdir/Nodes$pid.tcl"

        set matTag1 "1 0.0 0.0 -9.81";
        set elementType "stdBrick"
        eval "source $meshdir/Elements$pid.tcl"
    
        puts "Regular elements are created"
        set recordList [getNodeTags]
        set elerecordList [getEleTags]
    }
}
barrier

# ============================================================================
# creating DRM elements
# ============================================================================
if {$pid >= [expr $regcores +$structurecores] && $pid < [expr $regcores + $drmcores + $structurecores]} {
    model BasicBuilder -ndm 3 -ndf 3
    eval "source $meshdir/Nodes$pid.tcl"
    set matTag1 "1 0.0 0.0 -9.81";
    set elementType "stdBrick"
    eval "source $meshdir/Elements$pid.tcl"

    set elelist [getEleTags]
    # region 1 -ele $elelist -rayleigh $Damp_alpha $Damp_beta 0.0 0.0
    puts "DRM elements are created"
}
barrier

# ============================================================================
#  Adding pile elements
# ============================================================================
if {$HavePiles == "YES"} {
    if {$pid == $structurecores} {
        model BasicBuilder -ndm 3 -ndf 6
        set pileElements {}
        set pileNodes    {}
        # creating pile elements

        set secTag        1
        set transfTag     1
        set diameter      1. ; # pile diameter (m)
        set radius        [expr $diameter/2.]
        set pi            3.141593
        set Epile         1e10
        set nu            0.3
        set Gpile         [expr $Epile/(2*(1+$nu))]
        set Area          [expr ($diameter**2)*$pi/2.]
        set Iy            [expr ($diameter**4)*$pi/64.]
        set Iz            [expr ($diameter**4)*$pi/64.]
        set J             [expr ($diameter**4)*$pi/32.]
        set transfType    "PDelta"; # PDelta, Linear, Corotational
        
        # section Elastic $secTag $E $A $Iz $Iy $G $J 
        section Elastic $secTag $Epile $Area $Iz $Iy $Gpile $J

        set numpiles [llength $pilelist]
        # puts "Number of piles: $numpiles"
        set j 0
        foreach pile $pilelist {
            set xbottom [dict get $pile xbottom]
            set ybottom [dict get $pile ybottom]
            set zbottom [dict get $pile zbottom]
            set xtop    [dict get $pile xtop]
            set ytop    [dict get $pile ytop]
            set ztop    [dict get $pile ztop]
            set pilenumelements [dict get $pile numberofElements]
            set nPeri         8
            set nLong         8
            set numIntgrPts   5


            # creating pile nodes
            set pilenumnodes [expr $pilenumelements + 1]
            for {set i 1} {$i <= $pilenumnodes} {incr i} {
                set xcoord [expr $xbottom + ($xtop-$xbottom)*($i-1)/$pilenumelements]
                set ycoord [expr $ybottom + ($ytop-$ybottom)*($i-1)/$pilenumelements]
                set zcoord [expr $zbottom + ($ztop-$zbottom)*($i-1)/$pilenumelements]
                set Nodeid [expr $soilfoundation_num_nodes + $maxNodeTag + $j*$pilenumnodes + $i]
                node $Nodeid $xcoord $ycoord $zcoord
                lappend pileNodes $Nodeid
            }

            set P1x $xtop
            set P1y $ytop
            set P1z $ztop



            # normal vector to the pile
            set normalX [expr $xtop - $xbottom]
            set normalY [expr $ytop - $ybottom]
            set normalZ [expr $ztop - $zbottom]
            set norm    [expr sqrt($normalX**2 + $normalY**2 + $normalZ**2)]
            set normalX [expr $normalX/$norm]
            set normalY [expr $normalY/$norm]
            set normalZ [expr $normalZ/$norm]

            # ax + by + cz = d
            set d         [expr $normalX*$xtop + $normalY*$ytop + $normalZ*$ztop]

            # find another point on the plane
            set P2x [expr $xtop + 1]
            set P2y [expr $ytop + 1]
            set P2z [expr ($d - $normalX*$P2x - $normalY*$P2y)/$normalZ]

            set VecX_x [expr $P2x - $P1x]
            set VecX_y [expr $P2y - $P1y]
            set VecX_z [expr $P2z - $P1z]
            set norm   [expr sqrt($VecX_x**2 + $VecX_y**2 + $VecX_z**2)]
            set VecX_x [expr $VecX_x/$norm]
            set VecX_y [expr $VecX_y/$norm]
            set VecX_z [expr $VecX_z/$norm]





            set transfTag [expr $soilfoundation_num_cells + $maxEleTag + $j*$pilenumelements + 1]
            eval "geomTransf $transfType $transfTag $VecX_x $VecX_y $VecX_z"
            # creating pile elements
            for {set i 1} {$i < $pilenumnodes} {incr i} {
                set node1  [expr $soilfoundation_num_nodes + $maxNodeTag + $j*$pilenumnodes + $i]
                set node2  [expr $soilfoundation_num_nodes + $maxNodeTag + $j*$pilenumnodes + $i + 1]
                set eleTag [expr $soilfoundation_num_cells + $maxEleTag + $j*$pilenumelements + $i]
                element dispBeamColumn $eleTag $node1 $node2 $numIntgrPts $secTag $transfTag 
                lappend pileElements   $eleTag
            }

            # creating soil-pile interface
            set interfaceElemsFile "$outputdir/PileInterfaceElements$j.dat"
            set connecFile        "$outputdir/PileInterfaceConnections$j.dat"
            if {[file exists $interfaceElemsFile] == 1} { file delete $interfaceElemsFile }
            if {[file exists $connecFile] == 1} { file delete $connecFile }

            set interfaceElems {}
            set first_beam_ele [expr $soilfoundation_num_cells + $maxEleTag + $j*$pilenumelements + 1]
            set last_beam_ele  [expr $soilfoundation_num_cells + $maxEleTag + $j*$pilenumelements + $pilenumelements]
            set interfaceElems [generateInterfacePoints -beamEleRange $first_beam_ele $last_beam_ele   -gPenalty -shape circle -nP $nPeri -nL $nLong -crdTransf $transfTag -radius $radius -penaltyParam 1.0e12 -file $interfaceElemsFile -connectivity $connecFile  ]
            puts "interfaceElems: $interfaceElems"



            incr j
        }







        # puts "writing pile nodes and elements to a file"
        # write pile nodes to a file
        set pileNodesFile "$outputdir/PileNodes.dat"
        writeNodesinFile $pileNodesFile $pileNodes


        # write beam elements to a file
        set beamElemsFile "$outputdir/PileElements.dat"
        writeElesinFile $beamElemsFile $pileElements

        puts "Pile elements are created"
    }
}
barrier
set maxNodeTag [getMaxNodeTag $np $pid]
set maxEleTag  [getMaxEleTag $np $pid]

# ============================================================================
# Creating column and foundation elements
# ============================================================================
if {$HaveStructure == "YES" } {
    if {$pid == $structurecores && $structurecores > 0} {
        
        model BasicBuilder -ndm 3 -ndf 6
        set StrucFoundConnecElems {}
        set StrucFoundConnecNodes {}
        
        if {$StructureType == "STEEL" || $StructureType == "CONCRETE"} {set BaseColumnsFile "$outputdir/BaseColumnsNodes.dat"}
        
        # open BaseColumnsNodes.dat file
        set f [open $BaseColumnsFile "r"]
        set lines [split [read $f] "\n"]
        close $f

        set EmbeddingDepth 0.75
        set j 0
        foreach line $lines {
            # check if the line is empty
            if {[string length $line] == 0} {continue}
            incr j
            set nodeTag1 [lindex $line 0]
            set x        [lindex $line 1]
            set y        [lindex $line 2]
            set z        [lindex $line 3]
            node $nodeTag1 $x $y $z
            set nodeTag2 [expr $maxNodeTag + $j]
            set x        [expr $x]
            set y        [expr $y]
            set z        [expr $z - $EmbeddingDepth]
            node $nodeTag2 $x $y $z
            set numIntgrPts 3
            lappend StrucFoundConnecNodes $nodeTag1
            lappend StrucFoundConnecNodes $nodeTag2
            
            set nPeri         5
            set nLong         3
            set secTag        [expr $maxNodeTag + $j]
            set transfTag     [expr $maxNodeTag + $j]
            set diameter      0.25 ; 
            set radius        [expr $diameter/2.]
            set pi            3.141593
            set Epile         1e10
            set nu            0.3
            set Gpile         [expr $Epile/(2*(1+$nu))]
            set Area          [expr ($diameter**2)*$pi/2.]
            set Iy            [expr ($diameter**4)*$pi/64.]
            set Iz            [expr ($diameter**4)*$pi/64.]
            set J             [expr ($diameter**4)*$pi/32.]
            set transfType    "Linear"; # PDelta, Linear, Corotational


            section Elastic $secTag $Epile $Area $Iz $Iy $Gpile $J
            geomTransf $transfType $transfTag 1 0 0


            set eleTag [expr $maxEleTag + $j]
            element dispBeamColumn $eleTag $nodeTag2 $nodeTag1 $numIntgrPts $secTag $transfTag
            lappend StrucFoundConnecElems $eleTag

            # creating soil-pile interface
            set num [expr $j-1]
            set interfaceElemsFile "$outputdir/StructureFoundationInterfaceElements$num.dat"
            set connecFile         "$outputdir/StructureFoundationInterfaceConnections$num.dat"
            if {[file exists $interfaceElemsFile] == 1} { file delete $interfaceElemsFile }
            if {[file exists $connecFile] == 1} { file delete $connecFile }

            set interfaceElems {}
            set interfaceElems [generateInterfacePoints -beamEleRange $eleTag $eleTag -gPenalty -shape circle -nP $nPeri -nL $nLong -crdTransf $transfTag -radius $radius -penaltyParam 1.0e12 -file $interfaceElemsFile -connectivity $connecFile  ]
            set maxEleTag $interfaceElems
        }
        puts "Base columns and foundation elements are attached"

        # write pile nodes to a file
        set StrucFoundConnecNodesFile "$outputdir/StructureFoundationBeamNodes.dat"
        writeNodesinFile $StrucFoundConnecNodesFile $StrucFoundConnecNodes

        # write beam elements to a file
        set StrucFoundConnecElemsFile "$outputdir/StructureFoundationBeamElements.dat"
        writeElesinFile $StrucFoundConnecElemsFile $StrucFoundConnecElems
    }
}
barrier
# ============================================================================
# bulding PML layer
# ============================================================================

if {$HaveAbsorbingElements == "YES" && $pid >= [expr  $structurecores + $regcores + $drmcores] } {
    # create PML material
    set gamma           0.5                    ;# --- Coefficient gamma, newmark gamma = 0.5
    set beta            0.25                   ;# --- Coefficient beta,  newmark beta  = 0.25
    set eta             [expr 1.0/12.]         ;# --- Coefficient eta,   newmark eta   = 1/12 
    set E               $E                     ;# --- Young's modulus
    set nu              $nu                    ;# --- Poisson's Ratio
    set rho             $rho                   ;# --- Density
    set EleType         6                      ;# --- Element type, See line
    set PML_L           $pmltotalthickness     ;# --- Thickness of the PML
    set afp             2.0                    ;# --- Coefficient m, typically m = 2
    set PML_Rcoef       1.0e-8                ;# --- Coefficient R, typically R = 1e-8
    set RD_half_width_x [expr $llx/2.]         ;# --- Halfwidth of the regular domain in
    set RD_half_width_y [expr $lly/2.]         ;# --- Halfwidth of the regular domain in
    set RD_depth        [expr $llz/1.]         ;# --- Depth of the regular domain
    set pi              3.141593               ;# --- pi 

    set Damp_alpha      0.0    ;# --- Rayleigh damping coefficient alpha
    set Damp_beta       0.0    ;# --- Rayleigh damping coefficient beta

    
    if {$AbsorbingElements == "PML"} {
        model BasicBuilder -ndm 3 -ndf 9;
        set PMLmatTag1 "$eta $beta $gamma $E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"
        set elementType "PMLVISCOUS"
        # set elementType "PML"
        eval "source $meshdir/Nodes$pid.tcl"
        eval "source $meshdir/Elements$pid.tcl"

        # tie pml nodes to the regular nodes
        model BasicBuilder -ndm 3 -ndf 3;
        eval "source $meshdir/Boundary$pid.tcl"
    }
    if {$AbsorbingElements == "ASDA"} {
        model BasicBuilder -ndm 3 -ndf 3;
        set matTag1 "$G $nu $rho";
        set elementType "ASDAbsorbingBoundary3D"
        # set elementType "PML"
        eval "source $meshdir/Nodes$pid.tcl"
        eval "source $meshdir/Elements$pid.tcl"
    }
}
barrier


# ============================================================================
# creating fixities
# ============================================================================
if {$HaveAbsorbingElements == "YES"} {
    set totalThickness [expr $numdrmlayers*$drmthicknessx + $numpmllayers*$pmlthicknessx]
    fixX [expr -$lx/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixX [expr  $lx/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr -$ly/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr  $ly/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixZ [expr -$lz/1. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    puts "Boundary conditions are applied to the absorbing Nodes"
} else {
    set totalThickness [expr $numdrmlayers*$drmthicknessx]
    fixX [expr -$lx/2. - $totalThickness] 1 1 1;
    fixX [expr  $lx/2. + $totalThickness] 1 1 1;
    fixY [expr -$ly/2. - $totalThickness] 1 1 1;
    fixY [expr  $ly/2. + $totalThickness] 1 1 1;
    fixZ [expr -$lz/1. - $totalThickness] 1 1 1;
    puts "Boundary conditions are applied to the DRM Nodes"
}
barrier


# ============================================================================
# Gravity analysis
# ============================================================================
domainChange
constraints      Plain
numberer         ParallelRCM
system           Mumps -ICNTL14 400
test             NormDispIncr 1e-4 10 2
# algorithm        Linear -factorOnce 
algorithm        ModifiedNewton -factoronce
# algorithm        Newton
integrator       Newmark 0.5 0.25
analysis         Transient

for {set i 0} { $i < 100 } { incr i 1 } {
    if {$pid ==0 } {puts "Time step: $i";}
    set ok [analyze 1 0.01]
    if {$ok != 0} {
        puts "Analysis failed change the algorithm to ModifiedNewton"
        algorithm         ModifiedNewton
        analyze 1 0.01
    }
    if {$ok != 0} {
        puts "Analysis failed change the algorithm to Newton"
        algorithm         Newton
        analyze 1 0.01
    }
    if {$ok != 0} {puts "Analysis failed at time step: $i"; exit}
    algorithm             ModifiedNewton -factoronce     
}
# for {set i 0} { $i < 10 } { incr i 1 } {
#     if {$pid ==0 } {puts "Time step: $i";}
#     analyze 1 0.1
# }
# for {set i 0} { $i < 10 } { incr i 1 } {
#     if {$pid ==0 } {puts "Time step: $i";}
#     analyze 1 0.2
# }
# for {set i 0} { $i < 10 } { incr i 1 } {
#     if {$pid ==0 } {puts "Time step: $i";}
#     analyze 1 0.4
# }
# for {set i 0} { $i < 10 } { incr i 1 } {
#     if {$pid ==0 } {puts "Time step: $i";}
#     analyze 1 0.8
# }
# for {set i 0} { $i < 10 } { incr i 1 } {
#     if {$pid ==0 } {puts "Time step: $i";}
#     analyze 1 1.0
# }

barrier
loadConst -time 0.0
wipeAnalysis
puts "Gravity analysis is done"
# ============================================================================
# Gravirty recorders
# ============================================================================
if {$pid >= $structurecores && $pid < [expr $regcores + $structurecores]} {
    eval "recorder Node -file $outputdir/GravityNodeDisp$pid.out  -node $recordList -dof 1 2 3 disp"
    eval "recorder Node -file $outputdir/GravityNodeAccl$pid.out  -node $recordList -dof 1 2 3 accel"
    eval "recorder Node -file $outputdir/GravityNodeVelo$pid.out  -node $recordList -dof 1 2 3 vel"
    eval "recorder Element -file $outputdir/GravityElementStress$pid.out -ele $elerecordList stresses"
    eval "recorder Element -file $outputdir/GravityElementStrain$pid.out -ele $elerecordList strains"

    # print recordlist and elerecordlist to a file
    set f [open "$outputdir/nodeOuputTags$pid.out" "w+"]
    puts $f "$recordList"
    close $f
    set f [open "$outputdir/eleOuputTags$pid.out" "w+"]
    puts $f "$elerecordList"
    close $f
}
record
record
barrier
remove recorders
exit


 









# ============================================================================
# Post gravity settings
# ============================================================================

# ASDA Absorbing layer settings
if {$pid >= [expr $regcores + $drmcores] } {
    set abs_elements [getEleTags]
    eval "setParameter -val 1 -ele $abs_elements stage"
}

# set the initial displacements to zero
foreach node [getNodeTags] {
    foreach dof {1 2 3} {
        setNodeAccel $node $dof 0.0 -commit
        setNodeVel   $node $dof 0.0 -commit
        setNodeDisp  $node $dof 0.0 -commit
    }
}
# ============================================================================
# loading 
# ============================================================================

# if {$pid < $regcores} {  
if {$pid>=$regcores && $pid < [expr $regcores + $drmcores] } {
    set prop1 "1.0 1.0 0.001 1"; # factor crd_scale distance_tolerance do_coordinate_transformation
    set prop2 "1. 0. 0."; # T00, T01, T02
    set prop3 "0. 1. 0."; # T10, T11, T12
    set prop4 "0. 0. 1."; # T20, T21, T22
    set prop5 "0. 0. 0."; # x00, x01, x02
    set DRMtag  1
    eval "pattern H5DRM $DRMtag $DRMFile $prop1 $prop2 $prop3 $prop4 $prop5"
}

# # ============================================================================
# # recorders
# # ============================================================================
# source mesh/record.tcl
set deltaT 0.001
if {$pid < $regcores} {
    
    if {$pid == 0 } {set n [lindex $recordList 0] ;eval "recorder Node -file results/$modeltype/Time.out -time -dT $deltaT -node $n -dof 1 disp"}
    eval "recorder Node -file results/$modeltype/NodeDisp$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 disp"
    eval "recorder Node -file results/$modeltype/NodeAccl$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 accel"
    eval "recorder Node -file results/$modeltype/NodeVelo$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 vel"
    


    eval "recorder Element -file results/$modeltype/ElementStress$pid.out -dT $deltaT -ele $elerecordList stresses"
    eval "recorder Element -file results/$modeltype/ElementStrain$pid.out -dT $deltaT -ele $elerecordList strains"



    # print recordlist and elerecordlist to a file
    set f [open "results/$modeltype/nodeOuputTags$pid.out" "w+"]
    puts $f "$recordList"
    close $f
    set f [open "results/$modeltype/eleOuputTags$pid.out" "w+"]
    puts $f "$elerecordList"
    close $f


    # record the pile response
    if {$pid == 0} {
        eval "recorder Element -file results/$modeltype/Interfacepoints$pid.out -dT $deltaT -ele  $interfaceElems  displacement"
        eval "recorder Node    -file results/$modeltype/BeamDisp$pid.out        -dT $deltaT -node $beamNodes      -dof 1 2 3 disp"
        eval "recorder Node    -file results/$modeltype/BeamAccl$pid.out        -dT $deltaT -node $beamNodes      -dof 1 2 3 accel"
        eval "recorder Element -file results/$modeltype/BeamForce$pid.out       -dT $deltaT -ele  $beamElems       force"
    }
}

# if {$pid >= $regcores && $pid < [expr $regcores + $drmcores]} {
#     source "$meshdir/DRMNodes.tcl"
#     if {$DOPML=="YES" } {set modeltype "PML"} else {set modeltype "FIXED"}
#     eval "recorder Node -file results/$modeltype/NodeDispInternal$pid.out -dT $deltaT -node $internalNodes -dof 1 2 3 disp"
#     eval "recorder Node -file results/$modeltype/NodeAcclInternal$pid.out -dT $deltaT -node $internalNodes -dof 1 2 3 accel"
#     eval "recorder Node -file results/$modeltype/NodeVeloInternal$pid.out -dT $deltaT -node $internalNodes -dof 1 2 3 vel"
#     eval "recorder Node -file results/$modeltype/NodeDispExternal$pid.out -dT $deltaT -node $externalNodes -dof 1 2 3 disp"
#     eval "recorder Node -file results/$modeltype/NodeAcclExternal$pid.out -dT $deltaT -node $externalNodes -dof 1 2 3 accel"
#     eval "recorder Node -file results/$modeltype/NodeVeloExternal$pid.out -dT $deltaT -node $externalNodes -dof 1 2 3 vel"
# }
# recorde the initial state of the system
record

# ============================================================================
# Analysis 
# ============================================================================
domainChange

# constraints      Plain
constraints      Transformation
numberer         ParallelRCM
system           Mumps -ICNTL14 400
# system           ParallelProfileSPD
test             NormDispIncr 1e-4 10 2
algorithm        Linear -factorOnce 
# algorithm        ModifiedNewton -factoronce
integrator       Newmark 0.5 0.25
# integrator       Explicitdifference
analysis         Transient
set pi              3.141593              ;# --- pi 
set damp            0.05                  ;# --- Damping ratio
set omega1          [expr 2*$pi*0.2]      ; # lower frequency
set omega2          [expr 2*$pi*20]       ; # upper frequency
set Damp_alpha      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set Damp_beta       [expr 2*$damp/($omega1 + $omega2)]

# rayleigh $Damp_alpha $Damp_beta $Damp_alpha $Damp_beta
set startTime [clock milliseconds]

set dT 0.001
while { [getTime] < 20.00 } {
    if {$pid ==0 } {puts "Time: [getTime]";}
    analyze 1 $dT
}

set endTime [clock milliseconds]
set elapsedTime [expr {$endTime - $startTime}]
puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"



# if {$pid == 0} {
#     puts "natural frequency of the pile: $f Hz (assuming cantilever beam)"
#     puts "wavelegth: [expr $Vs/$f] m"
#     puts "Vs: $Vs"
# }
wipeAnalysis
remove recorders
remove loadPattern 2