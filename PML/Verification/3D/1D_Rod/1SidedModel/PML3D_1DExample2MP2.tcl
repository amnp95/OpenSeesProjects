# PML test 1D rod
# Written: Amin Pakzad, 2023


# get DOPML from command line
if {$argc > 0} {
    set DOPML [lindex $argv 0]
} else {
    set DOPML "NO"
}

# get PID from system
set pid [getPID]
set np  [getNP]


# printing some info
puts "DOPML = $DOPML"
puts "pid = $pid"
puts "np = $np"



# building nodes and elements
wipe 
set lx      25.0;
set ly      2.0;
set lz      2.0;
set dy      1.0;
set dx      1.0;
set dz      1.0;
set nx      [expr $lx/$dx ]
set ny      [expr $ly/$dy ]
set nz      [expr $lz/$dz ]


if {$pid==0} {
    model BasicBuilder -ndm 3 -ndf 3


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
                # puts "node $nodeTag $x $y $z;"
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
                # puts "element stdBrick $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 $materialTag;"

                element stdBrick $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 $materialTag;
                incr elementTag;
            }
        }
    }
}

if {$DOPML == "YES"} {
    # sending DOFlist to other processors
    if {$pid==0} {
        # sending data to the other processors
        send -pid 1 $Doflist
    } else {
        if {$pid==1} {
            # receiving data from the other processors
            recv -pid 0 Doflist
            # puts "Doflist = $Doflist"
        }
    }
    barrier
}



if {$DOPML == "YES"} {
    
    
    set PMLTotalLx 2.0

    # first partition of PML Elements
    if {$pid==1} {

        #create PML nodes and elements
        model BasicBuilder -ndm 3 -ndf 18;
        set dxPML   1.0;
        set dyPML   $dy;
        set dzPML   $dz;
        set lxPML   1.0;
        set lyPML   $ly;
        set lzPML   $lz;
        set nxPML   [expr $lxPML/$dxPML ]
        set nyPML   [expr $lyPML/$dyPML ]
        set nzPML   [expr $lzPML/$dzPML ]

        set PMLxlist {}
        set PMLylist {}
        set PMLzlist {}
        set PMLDoflist {}

        for {set i 0} {$i<=$nxPML} {incr i} {lappend PMLxlist [expr $dxPML*$i + $lx];}
        for {set i 0} {$i<=$nyPML} {incr i} {lappend PMLylist [expr $dyPML*$i];}
        for {set i 0} {$i<=$nzPML} {incr i} {lappend PMLzlist [expr $dzPML*$i];}


        set count 1;
        set nodeTag [expr int(($nx+1)*($ny+1)*($nz+1)+1)];
        set elementTag [expr int(($nx)*($ny)*($nz)+1)];


        foreach x $PMLxlist {
            foreach y $PMLylist {
                foreach z $PMLzlist {

                    node  $nodeTag $x $y $z;
                    fix $nodeTag 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;

                    if {$count == 1} {

                        lappend PMLDoflist [expr $nodeTag];

                        model BasicBuilder -ndm 3 -ndf 3
                        set tag [lindex $Doflist [expr int([llength $PMLDoflist] - 1)]]
                        node $tag $x $y $z;
                        # puts "$tag $x $y $z"


                        model BasicBuilder -ndm 3 -ndf 18
                    }
                    # puts "node $nodeTag $x $y $z;"
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
                    # puts "element PML $elementTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 2.08e+08 0.3 2000.0  6. 5.0 2.0 1.0e-8 25.0 25.0 25.0 0.0 0.0;"
                    incr elementTag;
                }
            }
        }



        # tie PML nodes to the main nodes
        model BasicBuilder -ndm 3 -ndf 18;
        for {set i 0} { $i < [llength $PMLDoflist] } { incr i 1 } {
            equalDOF [lindex $Doflist $i] [lindex $PMLDoflist $i] 1;
            # puts "equalDOF [lindex $Doflist $i] [lindex $PMLDoflist $i] 1;"
        }
    }

    barrier
    # Second partion of PML elements
    if {$pid==2} {
        
        #create PML nodes and elements
        model BasicBuilder -ndm 3 -ndf 18;
        
        set dxPML 1.0;
        set dyPML $dy;
        set dzPML $dz;
        set lxPML 1.0;
        set lyPML $ly;
        set lzPML $lz;
        set nxPML [expr $lxPML/$dxPML ]
        set nyPML [expr $lyPML/$dyPML ]
        set nzPML [expr $lzPML/$dzPML ]

        set PMLxlist {}
        set PMLylist {}
        set PMLzlist {}
        set PMLDoflist {}

        for {set i 0} {$i<=$nxPML} {incr i} {lappend PMLxlist [expr $dxPML*$i + $lx + ($PMLTotalLx - $lxPML)];}
        for {set i 0} {$i<=$nyPML} {incr i} {lappend PMLylist [expr $dyPML*$i];}
        for {set i 0} {$i<=$nzPML} {incr i} {lappend PMLzlist [expr $dzPML*$i];}


        set count 1;
        set nodeTag [expr int(($nx+1)*($ny+1)*($nz+1)+([expr ($PMLTotalLx - $lxPML)/$dxPML]+1)*($nyPML+1)*($nzPML+1)+1)];
        set elementTag [expr int(($nx)*($ny)*($nz)+([expr ($PMLTotalLx - $lxPML)/$dxPML])*($nyPML)*($nzPML)+1)];


        foreach x $PMLxlist {
            foreach y $PMLylist {
                foreach z $PMLzlist {
                    node  $nodeTag $x $y $z;
                    fix $nodeTag 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
                    if {$count == 1} {lappend PMLDoflist [expr $nodeTag]}
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
                    set node1 [expr int($x    *($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1) + ([expr ($PMLTotalLx - $lxPML)/$dxPML]+1)*($nyPML+1)*($nzPML+1))];
                    set node2 [expr int(($x+1)*($ny+1)*($nz+1) + $y    *($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1) + ([expr ($PMLTotalLx - $lxPML)/$dxPML]+1)*($nyPML+1)*($nzPML+1))];
                    set node3 [expr int(($x+1)*($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1) + ([expr ($PMLTotalLx - $lxPML)/$dxPML]+1)*($nyPML+1)*($nzPML+1))];
                    set node4 [expr int($x    *($ny+1)*($nz+1) + ($y+1)*($nz+1) + $z + 1 + ($nx+1)*($ny+1)*($nz+1) + ([expr ($PMLTotalLx - $lxPML)/$dxPML]+1)*($nyPML+1)*($nzPML+1))];
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



        # tie PML nodes to the pml nodes
        model BasicBuilder -ndm 3 -ndf 18;
        set j 1
        for {set i  [expr [llength $PMLDoflist]-1]} { $i >= 0 } { incr i -1 } {

            set tag [expr [lindex $PMLDoflist 0]-$j] 
            set coord [nodeCoord [lindex $PMLDoflist $i]];
            incr j;

            eval "node $tag $coord";
            puts "node $tag $coord;"

            equalDOF $tag [lindex $PMLDoflist $i] 1;
            puts "equalDOF $tag [lindex $PMLDoflist $i] 1;"
        }
    }

}

barrier
# STOP

# Creating fixities
if {$DOPML == "YES"} {
    if {$pid == 2} { fixX [expr $PMLTotalLx+$lx] 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; }
} else {
    if {$pid==0} {fixX $lx 1 0 0}
}



# loading 
set dT 0.001
if {$pid == 0} {
    timeSeries Path 1 -dt 0.001 -filePath force.dat -factor 1.0
    pattern Plain 1 1 {
        foreach node $Loadinglist {
            load $node 1.0 0.0 0.0
        }
    }
}





# recorders
if {$pid == 0} {
    eval "recorder Node -file NodeDispHead.out -time -node $Loadinglist  -dof 1 disp"
    eval "recorder Node -file NodeDispEnd.out  -time -node $Doflist      -dof 1 disp"
}


print "PML3D_1DExample2MP2_pid$pid.info" 
domainChange

 
# Analysis 
if {$DOPML == "YES"} {
    constraints      Plain
    numberer         ParallelRCM
    system           Mumps -ICNTL14 200
    test             NormDispIncr 1e-3 20 1
    algorithm        Linear -factorOnce
    integrator       Newmark 0.5 0.25
    analysis         Transient
    for {set i 0} { $i < 1000 } { incr i 1 } {
        if {$pid ==0 } {puts "Time step: $i"}
        analyze 1 $dT
    }
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



barrier
puts "Analysis completed by processor $pid"
barrier



