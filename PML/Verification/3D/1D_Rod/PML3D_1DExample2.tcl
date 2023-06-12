# PML test 1D rod
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





#create PML nodes and elements
if {$DOPML == "YES"} {
    model BasicBuilder -ndm 3 -ndf 18;
    
    set dxPML 1.0;
    set dyPML $dy;
    set dzPML $dz;
    set lxPML 5.0;
    set lyPML $ly;
    set lzPML $lz;
    set nxPML [expr $lxPML/$dxPML ]
    set nyPML [expr $lyPML/$dyPML ]
    set nzPML [expr $lzPML/$dzPML ]

    set PMLxlist {}
    set PMLylist {}
    set PMLzlist {}
    set PMLDoflist {}

    for {set i 0} {$i<=$nxPML} {incr i} {lappend PMLxlist [expr $dxPML*$i + $lx];}
    for {set i 0} {$i<=$nyPML} {incr i} {lappend PMLylist [expr $dyPML*$i];}
    for {set i 0} {$i<=$nzPML} {incr i} {lappend PMLzlist [expr $dzPML*$i];}


    set count 1;
    foreach x $PMLxlist {
        foreach y $PMLylist {
            foreach z $PMLzlist {
                node  $nodeTag $x $y $z;
                fix $nodeTag 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
                if {$count == 1} {lappend PMLDoflist [expr $nodeTag];}
                puts "node $nodeTag $x $y $z;"
                incr nodeTag;
            } 
        } 
        incr count;
    }



    # creating elements
    for {set x 0} { $x < $nxPML } { incr x 1 } {
        for {set y 0} { $y < $nyPML } { incr y 1 } {
            for {set z 0} { $z < $nzPML } { incr z 1 } {
                set node1 [expr int($x    *($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1))];
                set node2 [expr int(($x+1)*($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1))];
                set node3 [expr int(($x+1)*($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1))];
                set node4 [expr int($x    *($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1))];
                set node5 [expr $node1 + 1];
                set node6 [expr $node2 + 1];
                set node7 [expr $node3 + 1];
                set node8 [expr $node4 + 1];
                element PML $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 2.08e+08 0.3 2000.0  6. 5.0 2.0 1.0e-8 25.0 25.0 25.0 0.0 0.0;
                puts "element PML $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 2.08e+08 0.3 2000.0  6. 5.0 2.0 1.0e-8 25.0 25.0 25.0 0.0 0.0;"
                incr elementTag;
            }
        }
    }


    # tie PML nodes to the main nodes
    for {set i 0} { $i < [llength $PMLDoflist] } { incr i 1 } {
        equalDOF [lindex $Doflist $i] [lindex $PMLDoflist $i] 1;
        puts "equalDOF [lindex $Doflist $i] [lindex $PMLDoflist $i] 1;"
    }
}




# creating fixities
if {$DOPML == "YES"} {
    fixX [expr $lxPML+$lx] 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
} else {
    fixX $lx 1 0 0 0;
}


# loading 
set dT 0.001
timeSeries Path 1 -dt 0.001 -filePath force.dat -factor 1.0
pattern Plain 1 1 {
    foreach node $Loadinglist {
        load $node 1.0 0.0 0.0
    }
}


# recorders
eval "recorder Node -file node_disp2.out -time -node $Doflist -dof 1 disp"

print "PML3D_1DExample2.info" 

# Analysis 
constraints   Plain
numberer      RCM
integrator    Newmark 0.5 0.25
system        BandGEN
test          NormDispIncr 1.0e-6 10 1
algorithm     Linear
analysis      Transient


for {set i 0} { $i < 1000 } { incr i 1 } {
    puts "Time step: $i"
    analyze 1 $dT
}







