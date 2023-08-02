# ===================================================================== #
# 3D test model for the pml element modeling the plane strain field     #
# University of Washington, Department of Civil and Environmental Eng   #
# Geotechnical Eng Group, A. Pakzad, P. Arduino - Jun 2023              #
# Basic units are m, Ton(metric), s										#
# ===================================================================== #
# erase the boundary.tcl file if it exists
if {[file exists boundary.tcl]} {file delete boundary.tcl}
if {[file exists elements.tcl]} {file delete elements.tcl}
if {[file exists nodes.tcl]} {file delete nodes.tcl}
if {[file exists pmlelements.tcl]} {file delete pmlelements.tcl}
if {[file exists pmlnodes.tcl]} {file delete pmlnodes.tcl}
if {[file exists fixity.tcl.tcl]} {file delete fixity.tcl.tcl}
if {[file exists pmlfixity.tcl]} {file delete pmlfixity.tcl}
if {[file exists load.tcl]} {file delete load.tcl}

# get DOPML from command line
if {$argc > 0} {
    set DOPML [lindex $argv 0]
} else {
    set DOPML "NO"
}

# ============================================================================
# define geometry and meshing parameters
# ============================================================================
wipe 
set lx           10.0;
set ly           1.0;
set lz           10.0;
set dy           1.0;
set dx           1.0;
set dz           1.0;
set nx           [expr $lx/$dx ]
set ny           [expr $ly/$dy ]
set nz           [expr $lz/$dz ]
set pmlthickness 2.0


# ============================================================================
#  run the mesh generator
# ============================================================================
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
catch {eval "exec $pythonexec 3D_2DfieldMESh.py $lx $ly $lz $dx $dy $dz $pmlthickness"} result 
puts "result: $result"


#run 
# ============================================================================
# bulding regular elements
# ============================================================================
model BasicBuilder -ndm 3 -ndf 3

# create material
set materialTag 1;
nDMaterial ElasticIsotropic 1 2.08e8 0.3 2000.0
source nodes.tcl
source fixity.tcl
source elements.tcl


# ============================================================================
# bulding PML layer
# ============================================================================
#create PML nodes and elements
if {$DOPML == "YES"} {
    model BasicBuilder -ndm 3 -ndf 9;
    # create PML material
    set gamma           0.5                   ;# --- Coefficient gamma, newmark gamma = 0.5
    set beta            0.25                  ;# --- Coefficient beta,  newmark beta  = 0.25
    set eta             [expr 1.0/12.]        ;# --- Coefficient eta,   newmark eta   = 1/12 
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
    set PMLMaterial "$eta $beta $gamma $E $nu $rho $EleType $PML_L $afp $PML_Rcoef $RD_half_width_x $RD_half_width_y $RD_depth $Damp_alpha $Damp_beta"
    puts "PMLMaterial: $PMLMaterial"
    # 2.08e+08 0.3 2000.0  6. 5.0 2.0 1.0e-8 25.0 25.0 25.0 0.0 0.0;"

    source pmlnodes.tcl
    source pmlfixity.tcl
    source pmlelements.tcl


    # # delete the corner elements and nodes of the pml layer
    # remove element 101
    # remove element 122
    
    # remove sp 243 2
    # remove sp 243 11
    # remove sp 255 2
    # remove sp 255 11

    # remove node 243
    # remove node 243
    # remove node 255
    # remove node 255


    # remove sp 351 2
    # remove sp 351 11
    # remove sp 363 2
    # remove sp 363 11

    # remove node 351
    # remove node 351




    # tie pml nodes to the regular nodes
    model BasicBuilder -ndm 3 -ndf 3;
    source boundary.tcl
}
printGID "mesh.msh"

# ============================================================================
# creating fixities
# ============================================================================
if {$DOPML == "YES"} {
    fixX [expr -$lx/2. - $pmlthickness] 1 0 1 0 0 0 0 0 0;
    fixX [expr  $lx/2. + $pmlthickness] 1 0 1 0 0 0 0 0 0;
    fixZ [expr -$lz/1. - $pmlthickness] 1 0 1 0 0 0 0 0 0;
} else {
    fixX [expr -$lx/2.] 1 0 1;
    fixX [expr  $lx/2.] 1 0 1;
<<<<<<< HEAD
    fixZ [expr  $lz/1.] 1 0 1;
=======
    fixZ [expr -$lz/1.] 1 0 1;
>>>>>>> main
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
eval "recorder Node -file NodeDisp.out -time -node $recordList  -dof 3 disp"
eval "recorder Node -file NodeDispx.out -time -node $recordList  -dof 1 disp"
# ============================================================================
# Analysis 
# ============================================================================
print "PML3D_2DField.info" 
constraints   Plain
numberer      RCM
integrator    Newmark 0.5 0.25
# integrator    HHT 1.0
system        BandGEN
test          NormDispIncr 1.0e-5 20 1
# test          EnergyIncr 1.0e-5 20 2
# algorithm     Linear -factorOnce
# algorithm     ModifiedNewton -FactorOnce
algorithm     Newton
analysis      Transient


# timing the analysis
set start_time [clock seconds]
for {set i 0} { $i < 1000 } { incr i 1 } {
    puts "Time step: $i"
    analyze 1 $dT
}
<<<<<<< HEAD
set end_time [clock seconds]
puts "Total time: [expr $end_time - $start_time] seconds"
=======




>>>>>>> main
