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
set pmlcores      6
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
    catch {eval "exec $pythonexec ReducedDomain.py $lx $ly $lz $dx $dy $dz $pmlthicknessx $pmlthicknessy $pmlthicknessz $numpmllayers $drmthicknessx $drmthicknessy $drmthicknessz $numdrmlayers $regcores $drmcores $pmlcores $meshdir $outputdir"} result 
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
    set damp            0.00                  ;# --- Damping ratio
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
    file mkdir "results/PML"
    set totalThickness [expr $numdrmlayers*$drmthicknessx + $numpmllayers*$pmlthicknessx]
    fixX [expr -$lx/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixX [expr  $lx/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr -$ly/2. - $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr  $ly/2. + $totalThickness] 1 1 1 1 1 1 1 1 1;
    fixZ [expr -$lz/1. - $totalThickness] 1 1 1 1 1 1 1 1 1;
} else {
    file mkdir "results/FIXED"
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
        if {$pid == 0 } {set n [lindex $recordList 0] ;eval "recorder Node -file results/PML/Time.out -time -dT $deltaT -node $n -dof 1 disp"}
        eval "recorder Node -file results/PML/NodeDisp$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file results/PML/NodeAccl$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file results/PML/NodeVelo$pid.out  -dT $deltaT -node $recordList  -dof 1 2 3 vel"
    } else {
        if {$pid == 0 } {set n [lindex $recordList 0] ;eval "recorder Node -file results/FIXED/Time.out -time -dT $deltaT -node $n -dof 1 disp"}
        eval "recorder Node -file results/FIXED/NodeDisp$pid.out -dT $deltaT -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file results/FIXED/NodeAccl$pid.out -dT $deltaT -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file results/FIXED/NodeVelo$pid.out -dT $deltaT -node $recordList  -dof 1 2 3 vel"
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
    algorithm        Linear -factorOnce 
    # algorithm        ModifiedNewton -factoronce
    integrator       Newmark 0.5 0.25
    # integrator       HHT 0.67
    analysis         Transient
    # rayleigh 0.06220975551662957 0.00157579151576134 0.0 0.0
    set startTime [clock milliseconds]
    for {set i 0} { $i < 2000 } { incr i 1 } {
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
    for {set i 0} { $i < 2000 } { incr i 1 } {
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