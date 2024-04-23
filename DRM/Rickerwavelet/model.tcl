wipe 


set Xmeshsize 2.5
set Ymeshsize 2.5

set numElementsInX 1
set numElementsInY 16

set height   [expr $numElementsInY * $Ymeshsize]
set width    [expr $numElementsInX * $Xmeshsize]

set nomnodesX 2
set numnodesY [expr $numElementsInY + 1]


set E 3000000000.0
set poiss 0.3
set rho 2100.0
set thickness 1.0
set G [expr $E/(2.0*(1.0+$poiss))]



# builder
model Basic -ndm 2 -ndf 2


# nodes 
# ====================
set nodetag 0
for {set i 0} {$i < $numnodesY} {incr i} {
    incr nodetag
    set y [expr $i * $Ymeshsize]
    node $nodetag 0.0 $y
    puts "node $nodetag 0.0 $y"
    incr nodetag
    node $nodetag $width $y
    puts "node $nodetag $width $y"
    if {$i>1} {
        # laminar boundaries
        equalDOF $nodetag [expr $nodetag - 1] 1
    }
} 

# fix nodeTag (ndf constrValues)
fix 1 1 1
fix 2 1 1


# materials
# ====================
set matTag 1
nDMaterial ElasticIsotropic $matTag $E $poiss $rho


# elements
# ====================
set eleTag 0
set type "PlaneStrain"
set thick 1.0
set b1 0.0
set b2 -9.81
for {set i 0} {$i < $numElementsInY} {incr i} {
    incr eleTag
    set iNode [expr $i * 2 + 1]
    set jNode [expr $i * 2 + 2]
    set kNode [expr $i * 2 + 4]
    set lNode [expr $i * 2 + 3]
    element SSPquad $eleTag $iNode $jNode $kNode $lNode $matTag $type $thick $b1 $b2
}







# Static analysis (or quasti static)
# The absorbing boundaries now are in STAGE 0, so they act as constraints
constraints Transformation
numberer    RCM
system      UmfPack
test        NormUnbalance 0.0001 10 1
algorithm   Newton
integrator  LoadControl 1.0
analysis    Static
set ok [analyze 10]
if {$ok != 0} {
    error "Gravity analysis failed"
}
loadConst -time 0.0
wipeAnalysis


# recorders
# ====================  
set nodes {}
for {set i 0} {$i < $numnodesY} {incr i} {
    lappend nodes [expr $i * 2 + 1]
}
puts "nodes: $nodes"

eval "recorder Node -file \"acceleration.txt\" -time -node $nodes -dof 1 2 accel"
eval "recorder Node -file \"velocity.txt\" -time -node $nodes -dof 1 2 vel"
eval "recorder Node -file \"displacement.txt\" -time -node $nodes -dof 1 2 disp"


# time series
# ====================
set timeSeriesTag 1
set gmTag 1
timeSeries Path $timeSeriesTag -fileTime "ricker.dis" -filePath "ricker.time" -factor 10.0
pattern MultipleSupport 1 {
    groundMotion $gmTag Plain -disp $timeSeriesTag
    imposedSupportMotion 1 1 $gmTag
    imposedSupportMotion 2 1 $gmTag
}



# calculte rayleigh damping coefficients
# ====================
set pi 3.14159265358979323846
set damp 0.05
set omega1 [expr 2.0 * $pi * 1.0]
set omega2 [expr 2.0 * $pi * 20.0]

set Damp_alpha      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set Damp_beta       [expr 2*$damp/($omega1 + $omega2)]


# rayleigh $Damp_alpha $Damp_beta 0.0 0.0
# Dynamic analysis
# The absorbing boundaries now are in STAGE 0, so they act as constraints
constraints Transformation
numberer RCM
system UmfPack
test NormUnbalance 0.0001 10 1
algorithm Newton
integrator TRBDF2
analysis Transient
set nsteps 1000
set dt 0.001
set ok [analyze $nsteps $dt]
if {$ok != 0} {
    error "Dynamic analysis failed"
}