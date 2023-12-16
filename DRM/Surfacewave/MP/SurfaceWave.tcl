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
    if {$argc >1} {
        set ignore [lindex $argv 1]
    } else {
        set ignore "NO"
    }
} else {
    set DOPML "NO"
    set ignore "NO"
}

if {$pid==0} {
    puts "pid: $pid"
    puts "np: $np"
    puts "DOPML: $DOPML"
    puts "ignore: $ignore"
}


# ============================================================================
# define geometry and meshing parameters
# ============================================================================
wipe 
set lx           185.0;
set ly           10.0;
set lz           42.5;
set dy           2.5;
set dx           2.5;
set dz           2.5;
set nx           [expr $lx/$dx ]
set ny           [expr $ly/$dy ]
set nz           [expr $lz/$dz ]
set pmlthickness 5.0
set regcores     1
set pmlcores     1


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


    # run the 3D_2DfieldMESH.py to generate the mesh and check if it is finished using catch

    # passing the arguments to the python script: lx ly lz dx dy dz pmlthickness
    if {$ignore == "NO"} {
        catch {eval "exec $pythonexec SurfaceWaveMesh.py $regcores $pmlcores $lx $ly $lz $dx $dy $dz $pmlthickness"} result 
        puts "result: $result"
    }
}

# wait for the mesh generator to finish 
barrier

# ============================================================================
# bulding regular elements
# ============================================================================
set E           12937500000.0          ;# --- Young's modulus
set nu          0.25                   ;# --- Poisson's Ratio
set rho         2300.0                 ;# --- Density
set Vs          [expr {sqrt($E / (2.0 * (1.0 + $nu) * $rho))}]
puts "Vs: $Vs"

# create nodes and elements
model BasicBuilder -ndm 3 -ndf 3
set materialTag 1;

if {$pid < $regcores} {
    nDMaterial ElasticIsotropic 1 $E $nu $rho;
    eval "source nodes$pid.tcl"
    eval "source elements$pid.tcl"
}
# ============================================================================
# bulding PML layer
# ============================================================================
#create PML nodes and elements
if {$DOPML == "YES" && $pid >= $regcores} {
    puts "Creating PML elements"
    model BasicBuilder -ndm 3 -ndf 9;
    # create PML material
    set gamma           0.5                   ;# --- Coefficient gamma, newmark gamma = 0.5
    set beta            0.25                   ;# --- Coefficient beta,  newmark beta  = 0.25
    set eta             [expr 1.0/12.]        ;# --- Coefficient eta,   newmark eta   = 1/12 
    set E               $E                    ;# --- Young's modulus
    set nu              $nu                   ;# --- Poisson's Ratio
    set rho             $rho                  ;# --- Density
    set EleType         6                     ;# --- Element type, See line
    set PML_L           $pmlthickness         ;# --- Thickness of the PML
    set afp             2.0                   ;# --- Coefficient m, typically m = 2
    set PML_Rcoef       1.0e-8                ;# --- Coefficient R, typically R = 1e-8
    set RD_half_width_x [expr $lx/2.]         ;# --- Halfwidth of the regular domain in
    set RD_half_width_y [expr $ly/2.]         ;# --- Halfwidth of the regular domain in
    set RD_depth        [expr $lz/1.]         ;# --- Depth of the regular domain
    set pi              3.141593              ;# --- pi 
    set damp            0.05                  ;# --- Damping ratio
    set omega1          [expr 2*$pi*0.2]      ; # lower frequency
    set omega2          [expr 2*$pi*20]       ; # upper frequency
    set Damp_alpha      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
    set Damp_beta       [expr 2*$damp/($omega1 + $omega2)]
    puts "Damp_alpha: $Damp_alpha"
    puts "Damp_beta: $Damp_beta"
    # set Damp_alpha      0.06220975551662957   ;# --- Rayleigh damping coefficient alpha
    # set Damp_beta       0.00157579151576134   ;# --- Rayleigh damping coefficient beta


    set PMLMaterial "$eta $beta $gamma $E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"
    # set PMLMaterial "$E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"

    puts "PMLMaterial: $PMLMaterial"
    eval "source pmlnodes$pid.tcl"
    eval "source pmlelements$pid.tcl"

    # tie pml nodes to the regular nodes
    model BasicBuilder -ndm 3 -ndf 3;
    eval "source boundary$pid.tcl"
}

barrier

# ============================================================================
# creating fixities
# ============================================================================
if {$DOPML == "YES"} {

    fixX [expr -$lx/2. - $pmlthickness] 1 1 1 1 1 1 1 1 1;
    fixX [expr  $lx/2. + $pmlthickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr -$ly/2. - $pmlthickness] 1 1 1 1 1 1 1 1 1;
    fixY [expr  $ly/2. + $pmlthickness] 1 1 1 1 1 1 1 1 1;
    fixZ [expr -$lz/1. - $pmlthickness] 1 1 1 1 1 1 1 1 1;

} else {
    fixX [expr -$lx/2.] 1 1 1;
    fixX [expr  $lx/2.] 1 1 1;
    fixY [expr -$ly/2.] 1 1 1;
    fixY [expr  $ly/2.] 1 1 1;
    fixZ [expr -$lz/1.] 1 1 1;
}

# ============================================================================
# eigen analysis
# ============================================================================
# set lambda [eigen  8];
# set omega {}
# set f {}
# set T {}
# set pi 3.141593

# foreach lam $lambda {
# 	lappend omega [expr sqrt($lam)]
# 	lappend f [expr sqrt($lam)/(2*$pi)]
# 	lappend T [expr (2*$pi)/sqrt($lam)]
# }
# puts "lam: $lambda"
# puts "omega: $omega"
# puts "f: $f"
# puts "T: $T"
# exit

# ============================================================================
# loading 
# ============================================================================
set dT 0.001
if {$pid<$regcores} {
    pattern H5DRM 2 "SurfaceWave.h5drm" 1.0 1.0 0.001 1   1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0  0.0 0.0 0.0
    # pattern H5DRM 2 "SurfaceWave.h5drm" 1.0 1.0 0.001
}


# set dT 0.0005
# timeSeries Path 1 -dt $dT -filePath force.dat -factor -1.0
# pattern Plain 1 1 {
#     source load.tcl
# }


# # ============================================================================
# # DRM LOADING
# # ============================================================================
# source load.tcl
# set dT 0.01
# foreach tag $loadList {
#     set  num1 $tag
#     set  num2 [expr $tag + 10000]
#     set  num3 [expr $tag + 20000]
#     eval "timeSeries Path $num1 -dt 0.01 -filePath DRMLOAD/nodeX$tag.dat -factor 1.0"
#     eval "timeSeries Path $num2 -dt 0.01 -filePath DRMLOAD/nodeY$tag.dat -factor 1.0"
#     eval "timeSeries Path $num3 -dt 0.01 -filePath DRMLOAD/nodeZ$tag.dat -factor 1.0" 

#     pattern Plain $num1 $num1 {load $tag 1.0 0.0 0.0}
#     pattern Plain $num2 $num2 {load $tag 0.0 1.0 0.0}
#     pattern Plain $num3 $num3 {load $tag 0.0 0.0 1.0}
# }
# setTime 8.0

# ============================================================================
# recorders
# ============================================================================
source record.tcl
if {$pid < $regcores} {
    if {$DOPML=="YES" } {
        eval "recorder Node -file NodeDisp_PML.out -time -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file NodeAccl_PML.out -time -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file NodeVelo_PML.out -time -node $recordList  -dof 1 2 3 vel"
    } else {
        eval "recorder Node -file NodeDisp.out -time -node $recordList  -dof 1 2 3 disp"
        eval "recorder Node -file NodeAccl.out -time -node $recordList  -dof 1 2 3 accel"
        eval "recorder Node -file NodeVelo.out -time -node $recordList  -dof 1 2 3 vel"
    }
}
# ============================================================================
# Analysis 
# ============================================================================
# # Analysis 
# print "DRM_PML_3D.info" 
# domainChange
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
    for {set i 0} { $i < 500 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i";}
        analyze 1 $dT
    }
    set endTime [clock milliseconds]
    set elapsedTime [expr {$endTime - $startTime}]
    puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"
} else {
    constraints      Plain
    numberer         RCM
    # system           Mumps
    system           SparseSYM
    test             NormDispIncr 1e-4 3 2
    # algorithm        ModifiedNewton -factoronce
    algorithm        Linear -factorOnce 
    integrator       Newmark 0.5 0.25
    analysis         Transient

    # rayleigh 0.06220975551662957 0.00157579151576134 0.0 0.0
    set startTime [clock milliseconds]
    for {set i 0} { $i < 500 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i";}
        analyze 1 $dT
    }
    set endTime [clock milliseconds]
    set elapsedTime [expr {$endTime - $startTime}]
    puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"

    
}
puts "Vs: $Vs"
wipeAnalysis
remove recorders
remove loadPattern 2