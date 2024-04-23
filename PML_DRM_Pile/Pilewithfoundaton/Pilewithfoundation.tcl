# ===================================================================== #
# 3D test model for the pml element modeling the plane strain field     #
# University of Washington, Department of Civil and Environmental Eng   #
# Geotechnical Eng Group, A. Pakzad, P. Arduino - Jun 2023              #
# Basic units are m, Ton(metric), s										#
# ===================================================================== #
set pid [getPID]
set np  [getNP]

# get DOPML from command line
if {$argc >0} {
    set DOPML [lindex $argv 0]
} else {
    set DOPML "NO"
}
if {$pid==0} {
    puts "pid: $pid"
    puts "np: $np"
    puts "DOPML: $DOPML"
    file mkdir "results"
}
# ============================================================================
# define geometry and meshing parameters
# ============================================================================
wipe 
set llx           16.0;
set lly           6.0;
set llz           8.0;
set dy            0.5;
set dx            0.5;
set dz            0.5;
set drmthicknessx $dx
set drmthicknessy $dy
set drmthicknessz $dz
set numdrmlayers  1
set lx            [expr $llx - 2*$numdrmlayers*$drmthicknessx]
set ly            [expr $lly - 2*$numdrmlayers*$drmthicknessy]
set lz            [expr $llz - 1*$numdrmlayers*$drmthicknessz]
set nx            [expr $lx/$dx ]
set ny            [expr $ly/$dy ]
set nz            [expr $lz/$dz ]
set pmlthicknessx $dx
set pmlthicknessy $dy
set pmlthicknessz $dz
set numpmllayers  2
set regcores      1
set pmlcores      4
set drmcores      1
set meshdir       "OpenSeesmesh"
set outputdir     "results"
set pmltotalthickness [expr $numpmllayers*$pmlthicknessx]
barrier
# ============================================================================
#  run the mesh generator
# ============================================================================
# find that if python is exisiting in the system
if {$pid == 0} {
    set pythonexec "python3"
    if { [catch {exec python3 -V} python3_version] } {
        if { [catch {exec python -V} python_version] } {
            puts "Python is not installed in the system"
            exit
        } else {
            set pythonexec "python"
        }
    }
    puts "pythonexec: $pythonexec"

    # passing the arguments to the python script: 
    catch {eval "exec $pythonexec Pilewithfoundation.py $lx $ly $lz $dx $dy $dz $pmlthicknessx $pmlthicknessy $pmlthicknessz $numpmllayers $drmthicknessx $drmthicknessy $drmthicknessz $numdrmlayers $regcores $drmcores $pmlcores $meshdir $outputdir"} result 
    puts "result: $result"
    
}
# wait for the mesh generator to finish 
barrier

# ============================================================================
# bulding regular elements
# ============================================================================
model BasicBuilder -ndm 3 -ndf 3
set E           200000000.0            ;# --- Young's modulus Pa
set nu          0.30                   ;# --- Poisson's Ratio
set rho         2100.0                 ;# --- Density KG/M^3
set Vs          [expr {sqrt($E / (2.0 * (1.0 + $nu) * $rho))}]
if {$pid == 0} {
    puts "Vs: $Vs"
}
nDMaterial ElasticIsotropic 1 $E $nu $rho;
nDMaterial ElasticIsotropic 2 25000000000.0 0.22 2400.0; # --- concrete material
# ============================================================================
# create regular nodes and elements
# ============================================================================
if {$pid < $regcores} {
    model BasicBuilder -ndm 3 -ndf 3
    eval "source $meshdir/Nodes$pid.tcl"


    set matTag1 1;
    set matTag2 2;
    set elementType "SSPbrick"
    set elementType "stdBrick"
    eval "source $meshdir/Elements$pid.tcl"
    # set pi              3.141593              ;# --- pi 
    # set damp            0.02                  ;# --- Damping ratio
    # set omega1          [expr 2*$pi*0.2]      ; # lower frequency
    # set omega2          [expr 2*$pi*20]       ; # upper frequency
    # set Damp_alpha      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
    # set Damp_beta       [expr 2*$damp/($omega1 + $omega2)]
    # rayleigh $Damp_alpha $Damp_beta 0.0 0.0
    puts "Regular elements are created"
    set recordList [getNodeTags]
    set elerecordList [getEleTags]
}
barrier

# ============================================================================
# creating DRM elements
# ============================================================================
if {$pid >= $regcores && $pid < [expr $regcores + $drmcores]} {
    model BasicBuilder -ndm 3 -ndf 3
    eval "source $meshdir/Nodes$pid.tcl"
    set matTag1 1;
    set elementType "SSPbrick"
    set elementType "stdBrick"
    eval "source $meshdir/Elements$pid.tcl"
    puts "DRM elements are created"
}
barrier
# ============================================================================
#  Adding pile elements and a lumped mass at the top
# ============================================================================

if {$pid <$regcores} {
    set L1 2.0 ; # length of pile head (above ground surface) (m)
    set L2 0.75; # length of embedded pile (below ground surface) (m)
    set nElePile 11  ;# number of pile elements
    set eleSize [expr ($L1+$L2)/$nElePile]; # pile element length 
    set nNodePile [expr 1 + $nElePile]; # number of total pile nodes
    set foundationLevel 1.0; # the elevation of where the surface of the foundation is located


    # creating pile nodes
    model BasicBuilder -ndm 3 -ndf 6
    set beamNodes {}
    # creating nodes
    for {set i 1} {$i <= $nNodePile} {incr i} {
        set zCoord [expr $eleSize*(-$nNodePile + $i)+$L1 + $foundationLevel] 
        node [expr $i + 1000000]   0.00   0.00   $zCoord
        lappend beamNodes [expr $i + 1000000]
    }



    # creating pile elements
    set numIntgrPts   5
    set secTag        1
    set transfTag     1
    set diameter      0.8 ; # pile diameter (m)
    set radius        [expr $diameter/2.]
    set pi            3.141593
    set Epile         1e10
    set nu            0.3
    set Gpile         [expr $Epile/(2*(1+$nu))]
    set Area          [expr ($diameter**2)*$pi/2.]
    set Iy            [expr ($diameter**4)*$pi/64.]
    set Iz            [expr ($diameter**4)*$pi/64.]
    set J             [expr ($diameter**4)*$pi/32.]

    # geomTransf PDelta $transfTag $vecxzX $vecxzY $vecxzZ
    geomTransf   PDelta $transfTag  -1.0   0.0     0.0

    # section Elastic $secTag $E $A $Iz $Iy $G $J 
    section Elastic $secTag $Epile $Area $Iz $Iy $Gpile $J

  

    for {set i 1} {$i < $nNodePile} {incr i} {
        set node1  [expr $i + 1000000]
        set node2  [expr $i + 1000001]
        set eleTag [expr $i + 1000000]
        element dispBeamColumn $eleTag $node1 $node2 $numIntgrPts $secTag $transfTag 
    }

    # create a lumped mass at the top of the pile
    set lumpedMass 50000.0; # lumped mass at the top of the pile kg
    mass [expr $nNodePile + 1000000] $lumpedMass 0.0 0.0 0.0 0.0 0.0

    # natural frequency of the pile
    set K [expr ($Epile*$Iz*3)/($L1**3)] 
    set omega [expr sqrt($K/$lumpedMass)]
    set f [expr $omega/(2*$pi)]
    puts "natural frequency of the pile: $f Hz (assuming cantilever beam)"
    puts "wavelegth: [expr $Vs/$f] m"

    # ============================================================================
    # creating soil-pile interface
    # ============================================================================

    set nPeri 8              ;# number of points in the primeter of beam cross section
    set nLong 2              ;# number of points in the longitudinal direction of beam 

    set interfaceElemsFile "results/interfaceInfo.dat"
    set connecFile "results/connectivity.dat"
    if {[file exists $interfaceElemsFile] == 1} { file delete $interfaceElemsFile }
    if {[file exists $connecFile] == 1} { file delete $connecFile }

    set interfaceElems {}
    set first_beam_ele [expr 1000000 + 1]
    set last_beam_ele  [expr 1000000 + $nNodePile - 1]
    set interfaceElems [generateInterfacePoints -beamEleRange $first_beam_ele $last_beam_ele   -gPenalty -shape circle -nP $nPeri -nL $nLong -crdTransf $transfTag -radius $radius -penaltyParam 1.0e12 -file $interfaceElemsFile -connectivity $connecFile  ]
    puts "interfaceElems: $interfaceElems"
    

    # write beam nodes to a file
    set beamNodesFile "results/beamNodes.dat"
    if {[file exists $beamNodesFile] == 1} { file delete $beamNodesFile }
    set beamNodesID [open $beamNodesFile "w"]
    foreach node $beamNodes {
        puts $beamNodesID "$node [nodeCoord $node]"
    }
    close $beamNodesID

    # write beam elements to a file
    set beamElemsFile "results/beamElems.dat"
    if {[file exists $beamElemsFile] == 1} { file delete $beamElemsFile }
    set beamElemsID [open $beamElemsFile "w"]
    for {set i 1} {$i < $nNodePile} {incr i} {
        puts $beamElemsID "$i $i [expr $i+1]"
    }
    close $beamElemsID
}
if {$pid == 0} {puts "pile elements are created"}
barrier

# ============================================================================
# bulding PML layer
# ============================================================================
#create PML nodes and elements
if {$DOPML == "YES" && $pid >= [expr $regcores + $drmcores] } {
    # create PML material
    set gamma           0.5                   ;# --- Coefficient gamma, newmark gamma = 0.5
    set beta            0.25                   ;# --- Coefficient beta,  newmark beta  = 0.25
    set eta             [expr 1.0/12.]        ;# --- Coefficient eta,   newmark eta   = 1/12 
    set E               $E                    ;# --- Young's modulus
    set nu              $nu                   ;# --- Poisson's Ratio
    set rho             $rho                  ;# --- Density
    set EleType         6                     ;# --- Element type, See line
    set PML_L           $pmltotalthickness    ;# --- Thickness of the PML
    set afp             2.0                   ;# --- Coefficient m, typically m = 2
    set PML_Rcoef       1.0e-8                ;# --- Coefficient R, typically R = 1e-8
    set RD_half_width_x [expr $lx/2.]         ;# --- Halfwidth of the regular domain in
    set RD_half_width_y [expr $ly/2.]         ;# --- Halfwidth of the regular domain in
    set RD_depth        [expr $lz/1.]         ;# --- Depth of the regular domain
    set pi              3.141593              ;# --- pi 
    set damp            0.02                  ;# --- Damping ratio
    set omega1          [expr 2*$pi*0.2]      ; # lower frequency
    set omega2          [expr 2*$pi*20]       ; # upper frequency
    set Damp_alpha      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
    set Damp_beta       [expr 2*$damp/($omega1 + $omega2)]
    puts "Damp_alpha: $Damp_alpha"
    puts "Damp_beta: $Damp_beta"
    # set Damp_alpha      0.06220975551662957   ;# --- Rayleigh damping coefficient alpha
    # set Damp_beta       0.00157579151576134   ;# --- Rayleigh damping coefficient beta


    
    
    model BasicBuilder -ndm 3 -ndf 9;
    set matTag1 "$eta $beta $gamma $E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"
    set elementType "PMLVISCOUS"
    eval "source $meshdir/Nodes$pid.tcl"
    eval "source $meshdir/Elements$pid.tcl"


    # tie pml nodes to the regular nodes
    model BasicBuilder -ndm 3 -ndf 3;
    eval "source $meshdir/Boundary$pid.tcl"
}

barrier

# ============================================================================
# creating fixities
# ============================================================================
if {$DOPML == "YES"} {
    set totalThickness [expr $numdrmlayers*$drmthicknessx + $numpmllayers*$pmlthicknessx]
    fixX [expr -$lx/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixX [expr  $lx/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr -$ly/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr  $ly/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixZ [expr -$lz/1. - $totalThickness] 1 1 1 1 1 1 1 1 1;
} else {
    set totalThickness [expr $numdrmlayers*$drmthicknessx]
    fixX [expr -$lx/2. - $totalThickness] 1 1 1;
    fixX [expr  $lx/2. + $totalThickness] 1 1 1;
    fixY [expr -$ly/2. - $totalThickness] 1 1 1;
    fixY [expr  $ly/2. + $totalThickness] 1 1 1;
    fixZ [expr -$lz/1. - $totalThickness] 1 1 1;
}

# ============================================================================
# eigen analysis
# ============================================================================
# if {$DOPML == "NO"} {
#     set lambda [eigen  5];
#     set omega {}
#     set ff {}
#     set T {}
#     set pi 3.141593
#     puts "lambda: $lambda"

#     foreach lam $lambda {
#         lappend omega [expr sqrt($lam)]
#         lappend ff [expr sqrt($lam)/(2*$pi)]
#         lappend T [expr (2*$pi)/sqrt($lam)]
#     }
#     puts "lam: $lambda"
#     puts "omega: $omega"
#     puts "f: $ff"
#     puts "T: $T"
#     # exit
#     # Prompt the user
#     puts "Please press Enter to continue..."

#     # Wait for the user to press Enter
#     gets stdin
# }
# ============================================================================
# gravity analysis
# ============================================================================
constraints      Plain
numberer         ParallelRCM
system           Mumps
test             NormDispIncr 1e-4 10 0
algorithm        ModifiedNewton -factoronce
integrator       Newmark 0.6 0.25
analysis         Transient

for {set i 0} { $i < 100 } { incr i 1 } {
    if {$pid ==0 } {puts "Time step: $i";}
    analyze 1 0.01
}
for {set i 0} { $i < 10 } { incr i 1 } {
    if {$pid ==0 } {puts "Time step: $i";}
    analyze 1 0.1
}
for {set i 0} { $i < 10 } { incr i 1 } {
    if {$pid ==0 } {puts "Time step: $i";}
    analyze 1 1.0
}

barrier
loadConst -time 0.0
wipeAnalysis
puts "gravity analysis is done"
# ============================================================================
# loading 
# ============================================================================
set dT 0.001

if {$pid>=$regcores && $pid < [expr $regcores + $drmcores] } {
    # factor: multiply DRM dataset displacements and accelerations by this value
    # crd_scale: multiply (x, y, z) coordinates of dataset points by this value (typically for unit change between dataset and FE model)
    # distance_tolerance: tolerance for DRM point to FE mesh matching. Set to a size similar to the node distance near the DRM boundary.
    # do_coordinate_transformation: =1 apply a coordinate transformation to the (x, y, z) coordinates of dataset points before matching
    # T00, T01, T02, T10, T11, T12, T20, T21, T22: components of the transformation matrix
    # x00, x01, x02: x, y, z location of the top-center node of the DRM Box after transformation
    # pattern H5DRM $patternTag "/path/to/H5DRM/dataset" $factor $crd_scale $distance_tolerance $do_coordinate_transformation $T00 $T01 $T02 $T10 $T11 $T12 $T20 $T21 $T22 $x00 $x01 $x02
    # pattern H5DRM 2 "/home/amnp95/Projects/OpenSeesProjects/PML_DRM_Pile/EmbeddedInterfaceElement/MP/rickterDRM.h5drm" 1.0 1.0 0.001 0
    pattern H5DRM 2 "rickterDRM0_5.h5drm" 1.0 1.0 0.001 0
    # pattern H5DRM 2 "SurfaceWave.h5drm" 1.0 1.0 0.001
}

# # ============================================================================
# # recorders
# # ============================================================================
# source mesh/record.tcl
set deltaT [expr 2*$dT]
if {$pid < $regcores} {
    if {$DOPML=="YES" } {
        eval "recorder Node -file results/NodeDispPML$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file results/NodeAcclPML$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file results/NodeVeloPML$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 vel"
        eval "recorder Element -file results/ElementStressPML$pid.out   -time -dT $deltaT -ele $elerecordList  stress"
        eval "recorder Element -file results/InterfacepointsPML$pid.out -time -dT $deltaT -ele $interfaceElems  displacement"
        eval "recorder Node -file results/BeamDispPML$pid.out           -time -dT $deltaT -node $beamNodes  -dof 1 2 3 disp"
    } else {
        eval "recorder Node -file results/NodeDisp$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file results/NodeAccl$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file results/NodeVelo$pid.out -time -dT $deltaT -node $recordList  -dof 1 2 3 vel"
        eval "recorder Element -file results/ElementStress$pid.out   -time -dT $deltaT -ele $elerecordList  stress"
        eval "recorder Element -file results/Interfacepoints$pid.out -time -dT $deltaT -ele $interfaceElems  displacement"
        eval "recorder Node -file results/BeamDisp$pid.out           -time -dT $deltaT -node $beamNodes  -dof 1 2 3 disp"
    }
    # print recordlist and elerecordlist to a file
    set f [open "results/ouputTags$pid.out" "w+"]
    puts $f "recordList: $recordList"
    puts $f "elerecordList: $elerecordList"
    close $f

}

# ============================================================================
# Analysis 
# ============================================================================
domainChange
if {$DOPML == "YES"} {
    constraints      Plain
    numberer         ParallelRCM
    system           Mumps
    # system           BandGeneral
    test             NormDispIncr 1e-4 10 2
    # algorithm        Linear -factorOnce 
    algorithm        ModifiedNewton -factoronce
    integrator       Newmark 0.5 0.25
    # integrator       HHT 0.67
    analysis         Transient
    # rayleigh 0.06220975551662957 0.00157579151576134 0.0 0.0
    set startTime [clock milliseconds]
    for {set i 0} { $i < 20000 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i";}
        analyze 1 $dT
    }
    # set dT [expr 2*$dT]
    # for {set i 0} { $i < 4000 } { incr i 1 } {
    #     if {$pid ==0 } {puts "Time step: $i";}
    #     analyze 1 $dT
    # }

    set endTime [clock milliseconds]
    set elapsedTime [expr {$endTime - $startTime}]
    puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"

} else {
    constraints      Plain
    numberer         ParallelRCM
    system           Mumps
    # system           SparseSYM
    test             NormDispIncr 1e-4 3 2
    # algorithm        ModifiedNewton -factoronce
    algorithm        Linear -factorOnce 
    integrator       Newmark 0.5 0.25
    analysis         Transient

    # rayleigh 0.06220975551662957 0.00157579151576134 0.0 0.0
    set startTime [clock milliseconds]
    for {set i 0} { $i < 8000 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i";}
        analyze 1 $dT
    }
    set dT [expr 2*$dT]
    # for {set i 0} { $i < 4000 } { incr i 1 } {
    #     if {$pid ==0 } {puts "Time step: $i";}
    #     analyze 1 $dT
    # }
    set endTime [clock milliseconds]
    set elapsedTime [expr {$endTime - $startTime}]
    puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"

    
}

# if {$pid == 0} {
#     puts "natural frequency of the pile: $f Hz (assuming cantilever beam)"
#     puts "wavelegth: [expr $Vs/$f] m"
#     puts "Vs: $Vs"
# }
wipeAnalysis
remove recorders
remove loadPattern 2