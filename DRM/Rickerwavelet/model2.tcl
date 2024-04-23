# ===================================================
# User parameters
# ===================================================
# material parameters
set E 200000000.0
set poiss 0.30
set rho 2100.0
set thickness 5.0
set G [expr $E/(2.0*(1.0+$poiss))]
# domain size
set Lx 80.0
set Ly 140.0
# mesh size
set hx 2.5
set hy 0.5
# time increment
set dt 0.001
# predominant frequency of the Ricker Wavelet
set freq 10.0
# total duration of the dynamic analysis
set duration 2.0

# builder
model Basic -ndm 2 -ndf 2

# time series
# we want to apply a Ricker Wavelet with predominant frequency = 10 Hz.
# It should be applied as velocity
set pi [expr acos(-1.0)]
set wl [expr sqrt(3.0/2.0)/$pi/$freq*10.0]
set ndiv [expr int($wl/$dt)]
set dt [expr $wl/$ndiv.0]
set ts_vals {}
for {set i 0} {$i < $ndiv} {incr i} {
    set ix [expr $i.0*$dt-$wl/2.0]
    set iy [expr $ix*exp(-$pi*$pi*$freq*$freq*$ix*$ix)]
    lappend ts_vals $iy
}
set tsX 1
timeSeries Path $tsX -dt $dt -values $ts_vals  -factor 9.806

# material
set matTag 1
nDMaterial ElasticIsotropic $matTag $E $poiss $rho

# Define nodes on a regular grid with sizes hx-hy.
# For a more clear visualization we set the size of the absorbing elements larger.
# (note: the size of this element does not influence the results. The only constraint is that it
# should have a non-zero size!)
set nodes {}
set ndivx [expr int($Lx/$hx) + 2]; # add 2 layers of absorbing elements (left and right)
set ndivy [expr int($Ly/$hy) + 1]; # add 1 layer of absorbing elements (bottom)
set abs_h [expr $hx*2.0]
for {set j 0} {$j <= $ndivy} {incr j} {
    if {$j == 0} {set y [expr -$abs_h]} else {set y [expr ($j-1) * $hy]}
    for {set i 0} {$i <= [expr $ndivx]} {incr i} {
        if {$i == 0} {set x [expr -$abs_h]} elseif {$i == [expr $ndivx]} {set x [expr $Lx+$abs_h]} else {set x [expr ($i-1) * $hx]}
        node [expr $j*($ndivx+1)+$i+1] [expr $x-$Lx/2.0] $y
        if {[expr abs($x-$Lx/2.0) < 1e-6] && [expr $y > -1e-6]} {
            lappend nodes [expr $j*($ndivx+1)+$i+1]
        }
    }
}
# write corrdinates to file
set f [open "coords.txt" "w"]
foreach node $nodes {
    set coord [nodeCoord $node]
    set x [lindex $coord 0]
    set z [lindex $coord 1]
    puts $f "$x 0.0 $z"
}







# set outernodes  {}
# set innernodes {}

# set tol 1e-6
# foreach node [getNodeTags] {
#     set coord [nodeCoord $node]
#     set x [lindex $coord 0]
#     set z [lindex $coord 1]
#     set a [expr -$Lx/2.]
#     set b [expr $Lx/2.]
#     if {[expr abs($x -$a)] < $tol || [expr abs($x-$b)] < $tol} {
#         if {$z > -0.001} {
#             lappend outernodes $node
#         }
#     }
#     set a  [expr $a + $tol]
#     set b  [expr $b - $tol]
#     if {$x > $a && $x < $b && [expr abs($z)] < $tol} {
#         lappend outernodes $node
#     }

#     # inner nodes 
#     set a  [expr -$Lx/2. + $hx]
#     set b  [expr  $Lx/2. - $hx]
#     set c  [expr  0 + $hy -$tol]
#     if { [expr abs($x-$a)] < $tol || [expr abs($x-$b)] < $tol} {
#         if {$z > $c} {
#             lappend innernodes $node
#         }
#     }
#     set a  [expr $a + $tol]
#     set b  [expr $b - $tol]
#     set c  [expr 0 + $hy]
#     if {$x > $a && $x < $b && [expr abs($z-$c)] < $tol} {
#         lappend innernodes $node
#     }
# }

# set nodes {}
# set f [open "coords.txt" "w"]
# foreach node $outernodes {
#     set coord [nodeCoord $node]
#     set x [lindex $coord 0]
#     set z [lindex $coord 1]
#     puts $f "$x 0.0 $z 0"
#     lappend nodes $node
# }

# foreach node $innernodes {
#     set coord [nodeCoord $node]
#     set x [lindex $coord 0]
#     set z [lindex $coord 1]
#     puts $f "$x 0.0 $z 1"
#     lappend nodes $node
# }
close $f



# puts "Nodes on the base of the soil: $nodes"
# Define elements.
# Save absorbing elements tags in a list
set abs_elements {}
for {set j 0} {$j < $ndivy} {incr j} {
    # Yflag
    if {$j == 0} {set Yflag "B"} else {set Yflag ""}
    for {set i 0} {$i < [expr $ndivx]} {incr i} {
        # Tags
        set Etag [expr $j*($ndivx)+$i+1]
        set N1 [expr $j*($ndivx+1)+$i+1]
        set N2 [expr $N1+1]
        set N4 [expr ($j+1)*($ndivx+1)+$i+1]
        set N3 [expr $N4+1]
        # Xflag
        if {$i == 0} {set Xflag "L"} elseif {$i == [expr $ndivx-1]} {set Xflag "R"} else {set Xflag ""}
        set btype "$Xflag$Yflag"
        if {$btype != ""} {
            # absorbing element
            lappend abs_elements $Etag
            if {$Yflag != ""} {
                # bottom element
                element ASDAbsorbingBoundary2D $Etag $N1 $N2 $N3 $N4 $G $poiss $rho $thickness $btype -fx $tsX
            } else {
                # vertical element
                element ASDAbsorbingBoundary2D $Etag $N1 $N2 $N3 $N4 $G $poiss $rho $thickness $btype
            }
        } else {
            # soil element
            element quad $Etag $N1 $N2 $N3 $N4 $thickness PlaneStrain $matTag 0.0 0.0 0.0 [expr -9.806*$rho]
        }
    }
}




# Static analysis (or quasti static)
# The absorbing boundaries now are in STAGE 0, so they act as constraints
constraints Transformation
numberer RCM
system UmfPack
test NormUnbalance 0.0002 10 1
algorithm Newton
integrator LoadControl 0.05
analysis Static
set ok [analyze 20]
if {$ok != 0} {
    error "Gravity analysis failed"
}
loadConst -time 0.0
wipeAnalysis

# update absorbing elements to STAGE 1 (absorbing)
setParameter -val 1 -ele {*}$abs_elements stage

# recorders
# set soil_base [expr 1*($ndivx+1)+int($ndivx/2)+1]
# set soil_top [expr $ndivy*($ndivx+1)+int($ndivx/2)+1]
# recorder Node -file "soil_base.txt" -time -node $soil_base -dof 1 accel
# recorder Node -file "soil_top.txt" -time -node $soil_top -dof 1 accel

eval "recorder Node -file \"acceleration.txt\" -time -node $nodes -dof 1 accel"
eval "recorder Node -file \"velocity.txt\" -time -node $nodes -dof 1 vel"
eval "recorder Node -file \"displacement.txt\" -time -node $nodes -dof 1 disp"

# Dynamic analysis
# The absorbing boundaries now are in STAGE 0, so they act as constraints
constraints Transformation
numberer RCM
system UmfPack
test NormUnbalance 0.0001 10 1
algorithm Newton
integrator TRBDF2
analysis Transient
set dt 0.001
set nsteps [expr int($duration/$dt)]
set dt [expr $duration/$nsteps.0]
set ok [analyze $nsteps $dt]
if {$ok != 0} {
    error "Dynamic analysis failed"
}