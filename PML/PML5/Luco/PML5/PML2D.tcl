# ===================================================================== #
# 2D test model for the pml element modeling the plane strain field     #
# University of Washington, Department of Civil and Environmental Eng   #
# Geotechnical Eng Group, A. Pakzad, P. Arduino - Jun 2023              #
# Basic units are m, kg(metric), s										#
# ===================================================================== #
# erase the boundary.tcl file if it exists
set pid [getPID]
set np  [getNP]
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
set lx              10.0;
set ly              4.0;
set dy              1.0;
set dx              1.0;
set nx              [expr $lx/$dx ]
set ny              [expr $ly/$dy ]
set pmlthickness    8.0
set pmlthicknessX   $pmlthickness
set pmlthicknessY   $pmlthickness
set dxpml           $dx
set dypml           $dy
set regcores        1
set pmlcores        0


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


# run the PML2D.py to generate the mesh and check if it is finished using catch
# passing the arguments to the python script: lx ly dx dy pmlthicknessX pmlthicknessY  dxpml dypml
catch {eval "exec $pythonexec PML2D.py $lx $ly $dx $dy $pmlthicknessX  $pmlthicknessY $dxpml $dypml $regcores $pmlcores"} result 
# puts "result: $result"

# ============================================================================
# bulding Soil elements
# ============================================================================
model BasicBuilder -ndm 2 -ndf 2
set pid [getPID]
puts "pid: $pid"

# create soil material
set Thick                  1.0
set Type                   "PlaneStrain"
set bx                     0.0
set by                     0.0
set materialTag            1;
set VSSoil                 229;
set nuSoil                 0.25;
set rhoSoil                1800;
set ESoil                  [expr 2 * $rhoSoil * $VSSoil * $VSSoil * (1 + $nuSoil) ];
nDMaterial ElasticIsotropic $materialTag $ESoil $nuSoil $rhoSoil
set Soilmaterial "$Thick $Type $materialTag 0.0 0.0 $bx $by"


# create structural material
set StructureMaterialTag  2;
set VSStructure                 229000;
set nuStructure                 0.25;
set rhoStructure                1;
set EStructure                  [expr 2 * $rhoStructure * $VSStructure * $VSStructure * (1 + $nuStructure)];
nDMaterial ElasticIsotropic $StructureMaterialTag  $EStructure $nuStructure $rhoStructure
set Structurematerial "$Thick $Type $StructureMaterialTag 0.0 0.0 $bx $by"

eval "source nodes$pid.tcl"
eval "source elements$pid.tcl"
eval "source structurenodes$pid.tcl"
eval "source structureelements$pid.tcl"
eval "source interface.tcl"
# ============================================================================
# bulding PML layer
# ============================================================================
#create PML nodes and elements
if {$DOPML == "YES"} {
    # create PML material
    set E               $ESoil             ;# --- Young's modulus
    set nu              $nuSoil                   ;# --- Poisson's Ratio
    set rho             $rhoSoil                ;# --- Density
    set RD_half_width_x [expr $lx/2.]         ;# --- Halfwidth of the regular domain in
    set Depth [expr $ly]            ;# --- Halfwidth of the regular domain in
    set PMLmaterial "$E $nu $rho $pmlthicknessX $pmlthicknessY $RD_half_width_x $Depth"
    
    puts "PMLMaterial: $PMLmaterial"

    source pmlnodes$pid.tcl
    model BasicBuilder -ndm 2 -ndf 5;
    source pmlcenternodes$pid.tcl
    model BasicBuilder -ndm 2 -ndf 2;
    source pmlelements$pid.tcl

    # tie pml nodes to the regular nodes
    model BasicBuilder -ndm 2 -ndf 2;
    source boundary$pid.tcl
}

# ============================================================================
# creating fixities
# ============================================================================
if {$DOPML == "YES"} {
    fixX [expr -$lx/2. - $pmlthickness] 1 1 0 0 0;
    fixX [expr  $lx/2. + $pmlthickness] 1 1 0 0 0;
    fixY [expr -$ly/1. - $pmlthickness] 1 1 0 0 0;
} else {
    fixX [expr -$lx/2.] 1 1;
    fixX [expr  $lx/2.] 1 1;
    fixY [expr -$ly/1.] 1 1;
}
# ============================================================================
# Eigen analysis
# ============================================================================
# perform eigen analysis
if {$DOPML == "Eigen"} {
    set lamda [eigen 1]
    set omega [expr sqrt($lamda)]
    set f     [expr $omega/(2.*3.14159265359)]
    puts "lamda: $lamda"
    puts "omega: $omega"
    puts "f: $f"
    exit
}
# ============================================================================
# loading 
# ============================================================================
set dT 0.001
timeSeries Path 1 -dt 0.001 -filePath "RichterHF.dat" -factor -1000.0
pattern Plain 1 1 {
    source load.tcl
}

# ============================================================================
# recorders
# ============================================================================

eval "recorder Node -file NodeDisp$pid.out -time -node $recordList  -dof 2 disp"
# ============================================================================
# Analysis
# ============================================================================
constraints      Plain
numberer         RCM
system           BandGeneral
test             NormDispIncr 1e-4 10 0
algorithm        Linear -factorOnce 
# algorithm        ModifiedNewton -factoronce 
# algorithm        ModifiedNewton
# algorithm        Newton 
integrator       Newmark 0.5 0.25
# integrator       HHT 1.0
analysis         Transient
analyze 1 $dT
set startTime [clock milliseconds]
for {set i 0} { $i < 1500 } { incr i 1 } {
    if {$pid ==0 } {puts "Time step: $i";}
    analyze 1 $dT
}
set endTime [clock milliseconds]
set elapsedTime [expr {$endTime - $startTime}]
puts "Elapsed time: [expr $elapsedTime/1000.] seconds in $pid"
