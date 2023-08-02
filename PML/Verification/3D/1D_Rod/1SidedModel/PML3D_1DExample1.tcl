foreach doPML {"YES"} {

    wipe
    
    set eleSize 1.0
    set h 1.0
    set b -1.0

    #
    # first add regular nodes and elements (2dof nodes & quads)
    #

    model Basic -ndm 3 -ndf 3

    
    for {set i 0; set nodeTag 1} {$i <= 25} {incr i 1; incr nodeTag 4} {
	set xCrd [expr $i * $eleSize - 25.0]
    	node $nodeTag           0. 0. $xCrd
    	node [expr $nodeTag+1]  $h 0. $xCrd
        node [expr $nodeTag+2]  $h $h $xCrd 
        node [expr $nodeTag+3]  0. $h $xCrd 
		puts "node $nodeTag 0. 0. $xCrd"
    	puts "node [expr $nodeTag+1]  $h 0. $xCrd"
        puts "node [expr $nodeTag+2]  $h $h $xCrd" 
        puts "node [expr $nodeTag+3]  0. $h $xCrd" 
        fix $nodeTag 1 1 0
        fix [expr $nodeTag+1] 1 1 0
        fix [expr $nodeTag+2] 1 1 0
        fix [expr $nodeTag+3] 1 1 0
    }
    
    nDMaterial ElasticIsotropic 1 2.08e8 0.3 2000.0
    for {set i 1; set iNode 1;} {$i <= 25} {incr i 1} {
    element stdBrick $i $iNode [expr $iNode+1] [expr $iNode+2] [expr $iNode+3] [expr $iNode+4] [expr $iNode+5] [expr $iNode+6] [expr $iNode+7] 1
	incr iNode 4
    }
    
    if {$doPML == "YES"} {

	#
	# add nodes and elements for PML (5dof nodes and PML ele)
	#   note: tie nodes that overlap with quad elements
	#

	model Basic -ndm 3 -ndf 18

    for {set i 0; set nodeTag 105} {$i <= 5} {incr i 1; incr nodeTag 4} {
    set xCrd [expr $i * $eleSize - 30.0]
        node $nodeTag           0. 0. $xCrd
        node [expr $nodeTag+1]  $h 0. $xCrd
        node [expr $nodeTag+2]  $h $h $xCrd 
        node [expr $nodeTag+3]  0. $h $xCrd 
		puts "node $nodeTag           0. 0. $xCrd"
        puts "node [expr $nodeTag+1]  $h 0. $xCrd"
        puts "node [expr $nodeTag+2]  $h $h $xCrd" 
        puts "node [expr $nodeTag+3]  0. $h $xCrd" 
        fix $nodeTag 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
        fix [expr $nodeTag+1] 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
        fix [expr $nodeTag+2] 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
        fix [expr $nodeTag+3] 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
    }
	
	for {set i 26; set iNode 105;} {$i <= 30} {incr i 1} {
		puts "$iNode"
	    element PML $i $iNode [expr $iNode+1] [expr $iNode+2] [expr $iNode+3] [expr $iNode+4] [expr $iNode+5] [expr $iNode+6] [expr $iNode+7] 2.08e+08 0.3 2000.0  6. 5.0 2.0 1.0e-8 25.0 25.0 25.0 0.0 0.0
	    incr iNode 4
	}


    # puts i am fixing PML
	
	fix 105 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
	fix 106 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
    fix 107 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
    fix 108 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 

    equalDOF 1 125 3
    equalDOF 2 126 3
    equalDOF 3 127 3
    equalDOF 4 128 3


    } else {

	fix 1 0 0 1
	fix 2 0 0 1
    fix 3 0 0 1
    fix 4 0 0 1

    }

    #
    # add load at end
    #

    set dT 0.001
    timeSeries Path 1 -dt 0.001 -filePath force.dat -factor 1.0
    pattern Plain 1 1 {
	load 101 0.0 0.0 1.0
	load 102 0.0 0.0 1.0
    load 103 0.0 0.0 1.0
    load 104 0.0 0.0 1.0
    }

    #
    # output node response at tip
    #

    recorder Node -file node_disp.out -time -node 101 125 -dof 3 disp

    #
    # perform analysis
    #
    print "PML3D_1DExample1.info" 
    constraints Plain
    numberer Plain
    integrator Newmark 0.5 0.25
    system BandGEN
    test NormDispIncr 1.0e-6 10 1
    algorithm Linear
    analysis Transient
    
    set ok [analyze 1000 0.001]

    #
    # spit out success message
    #

    if {$ok == 0} {
	puts "doPML? $doPML SUCCESS"
    } else {
	puts "doPML? $doPML FAIL"
    }
}
wipe