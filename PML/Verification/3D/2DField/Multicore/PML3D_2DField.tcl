# ===================================================================== #
# 3D test model for the pml element modeling the plane strain field     #
# University of Washington, Department of Civil and Environmental Eng   #
# Geotechnical Eng Group, A. Pakzad, P. Arduino - Jun 2023              #
# Basic units are m, Ton(metric), s										#
# ===================================================================== #
set pid [getPID]
set np  [getNP]

# get DOPML from command line
if {$argc > 0} {
    set DOPML [lindex $argv 0]
} else {
    set DOPML "NO"
}

if {$pid==0} {
    puts "pid: $pid"
    puts "np: $np"
    puts "DOPML: $DOPML"
}

# ============================================================================
# define geometry and meshing parameters
# ============================================================================
wipe 
set lx           20.0;
set ly           1.0;
set lz           10.0;
set dy           1.0;
set dx           1.0;
set dz           1.0;
set nx           [expr $lx/$dx ]
set ny           [expr $ly/$dy ]
set nz           [expr $lz/$dz ]
set pmlthickness 2.0
set regcores     1
set pmlcores     1

barrier
# ============================================================================
#  run the mesh generator
# ============================================================================
if {$pid==0} {
    # find that if python is exisiting in the system
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
    catch {eval "exec $pythonexec 3D_2DfieldMESh.py $regcores $pmlcores $lx $ly $lz $dx $dy $dz $pmlthickness"} result 
    puts "result: $result"

}
# wait for the mesh generator to finish 
barrier

# ============================================================================
# bulding regular elements
# ============================================================================
if {$pid < $regcores} {
    # create nodes and elements
    model BasicBuilder -ndm 3 -ndf 3
    set materialTag 1;
    nDMaterial ElasticIsotropic 1 2.08e8 0.3 2000.0
    eval "source nodes$pid.tcl"
    eval "source fixity$pid.tcl"
    eval "source elements$pid.tcl"
}
barrier
# ============================================================================
# bulding PML layer
# ============================================================================
#create PML nodes and elements
if {$DOPML == "YES" && $pid >= $regcores} {

    model BasicBuilder -ndm 3 -ndf 18;
    # create PML material
    set E               2.08e+08              ;# --- Young's modulus
    set nu              0.3                   ;# --- Poisson's Ratio
    set rho             2000.0                ;# --- Density
    set EleType         6                     ;# --- Element type, See line
    set PML_L           $pmlthickness         ;# --- Thickness of the PML
    set afp             2.0                   ;# --- Coefficient m, typically m = 2
    set PML_Rcoef       1.0e-8                ;# --- Coefficient R, typically R = 1e-8
    set RD_half_width_x [expr $lx/2.]         ;# --- Halfwidth of the regular domain in
    set RD_half_width_y [expr $ly/2.]         ;# --- Halfwidth of the regular domain in
    set RD_depth        [expr $lz/1.]         ;# --- Depth of the regular domain
    set Damp_alpha      0.0                   ;# --- Rayleigh damping coefficient alpha
    set Damp_beta       0.0                   ;# --- Rayleigh damping coefficient beta 
    set PMLMaterial "$E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"

    eval "source pmlnodes$pid.tcl"
    eval "source pmlfixity$pid.tcl"
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
    if {$pid >=$regcores} {
        fixX [expr -$lx/2. - $pmlthickness] 1 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0;
        fixX [expr  $lx/2. + $pmlthickness] 1 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0;
        fixZ [expr -$lz/1. - $pmlthickness] 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
    }
} else {
    if {$pid < $regcores} {
        fixX [expr -$lx/2.] 1 0 1;
        fixX [expr  $lx/2.] 1 0 1;
        fixZ [expr  $lz/1.] 0 0 1;
    }
}

# ============================================================================
# loading 
# ============================================================================
set dT 0.001
timeSeries Path 1 -dt 0.001 -filePath force.dat -factor 1.0
pattern Plain 1 1 {
    source load.tcl
}


# ============================================================================
# recorders
# ============================================================================
if {$pid == 0 } {
    eval "recorder Node -file NodeDisp.out -time -node $recordList  -dof 3 disp"
}
# ============================================================================
# Analysis 
# ============================================================================
# Analysis 
print "PML3D_1DExample2MP1_pid$pid.info" 
domainChange
if {$DOPML == "YES"} {
    constraints      Plain
    numberer         ParallelRCM
    system           Mumps -ICNTL14 200
    test             NormDispIncr 1e-3 3 0
    algorithm        Linear -factorOnce 
    # algorithm        ModifiedNewton -factoronce 
    integrator       Newmark 0.5 0.25
    analysis         Transient
    set startTime [clock milliseconds]
    for {set i 0} { $i < 1000 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i";}
        analyze 1 $dT
    }
    set endTime [clock milliseconds]
    set elapsedTime [expr {$endTime - $startTime}]
    puts "Elapsed time: $elapsedTime milliseconds in $pid"
} else {
    if {$pid ==0 } {
        numberer         RCM
        system           BandGEN
        constraints      Plain
        test             NormDispIncr 1e-3 20 1
        algorithm        Linear -factorOnce
        integrator       Newmark 0.5 0.25
        analysis         Transient
    
        for {set i 0} { $i < 1000 } { incr i 1 } {
            puts "Time step: $i"
            analyze 1 $dT
        }
    }
}