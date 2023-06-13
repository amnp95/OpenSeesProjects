# PML test 2D field example
# Written: Amin Pakzad, 2023


# get DOPML from command line
if {$argc > 0} {
    set DOPML [lindex $argv 0]
} else {
    set DOPML "NO"
}

# building nodes and elements
wipe 
model BasicBuilder -ndm 3 -ndf 3

set lx      25.0;
set ly      2.0;
set lz      2.0;
set dy      1.0;
set dx      1.0;
set dz      1.0;
set nx      [expr $lx/$dx ]
set ny      [expr $ly/$dy ]
set nz      [expr $lz/$dz ]


set xlist       {}
set ylist       {}
set zlist       {}
set Doflist     {}
set Loadinglist {}


for {set i 0} { $i <= $nx } { incr i} {lappend xlist [expr $dx*$i];}
for {set i 0} { $i <= $ny } { incr i} {lappend ylist [expr $dy*$i];}
for {set i 0} { $i <= $nz } { incr i} {lappend zlist [expr $dz*$i];}



set nodeTag    1;
set elementTag 1;




# creating nodes
set count 1;
foreach x $xlist {
    foreach y $ylist {
        foreach z $zlist {
            node  $nodeTag $x $y $z;
            fix $nodeTag 0 1 1
            puts "node $nodeTag $x $y $z;"
            if {$count == 1} {lappend Loadinglist [expr $nodeTag];}
            if {$count == [expr $nx+1]} {lappend Doflist [expr $nodeTag];}
            incr nodeTag;
        } 
    } 
    incr count;
}



# create material
set materialTag 1;
nDMaterial ElasticIsotropic 1 2.08e8 0.3 2000.0



# creating elements for 4 first nodes in x-y plane and then other in z direction (counterclock wise)
for {set x 0} {$x < $nx} {incr x 1} {
    for {set y 0} {$y < $ny} {incr y 1} {
        for {set z 0} {$z < $nz} {incr z 1} {
            set node1 [expr int($x    *($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1)];
            set node2 [expr int(($x+1)*($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1)];
            set node3 [expr int(($x+1)*($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1)];
            set node4 [expr int($x    *($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1)];
            set node5 [expr $node1 + 1];
            set node6 [expr $node2 + 1];
            set node7 [expr $node3 + 1];
            set node8 [expr $node4 + 1];
            puts "element stdBrick $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 $materialTag;"

            element stdBrick $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 $materialTag;
            incr elementTag;
        }
    }
}