
####################################################################
#  This file was made using the UWquad2D problem type in GiD.      #
#  Created by:  Chris McGann and Pedro Arduino, 2010               #
#              -University of Washington-                          #
####################################################################

wipe

#-----------------------------------------------------------------------------------------
#  1. CREATE SOIL NODES AND FIXITIES
#-----------------------------------------------------------------------------------------
model BasicBuilder -ndm 2 -ndf 2

# define soil nodes
source nodeInfo.dat
source PMLnodeInfo.dat
puts "Finished creating all -ndf 2 and 5 nodes..."

# define fixities for soil nodes

# define equal degrees of freedom for free-field columns
equalDOF 91 88  1 2
equalDOF 93 92  1 2
equalDOF 96 97  1 2
equalDOF 102 101  1 2
equalDOF 106 105  1 2
equalDOF 117 116  1 2
equalDOF 125 124  1 2
equalDOF 140 139  1 2
equalDOF 150 149  1 2
equalDOF 167 166  1 2
equalDOF 181 182  1 2
equalDOF 199 200  1 2
equalDOF 220 219  1 2
equalDOF 243 242  1 2
equalDOF 267 266  1 2
equalDOF 288 287  1 2
equalDOF 315 314  1 2
equalDOF 342 343  1 2
equalDOF 371 372  1 2
equalDOF 405 404  1 2
equalDOF 436 437  1 2
equalDOF 469 471  1 2
equalDOF 505 506  1 2
equalDOF 545 544  1 2
equalDOF 579 580  1 2
equalDOF 620 621  1 2
equalDOF 665 663  1 2
equalDOF 705 707  1 2
equalDOF 753 752  1 2
equalDOF 802 801  1 2
equalDOF 848 847  1 2
equalDOF 898 899  1 2
equalDOF 948 949  1 2
equalDOF 1004 1001  1 2
equalDOF 1058 1059  1 2
equalDOF 1116 1115  1 2
equalDOF 1174 1172  1 2
equalDOF 1235 1234  1 2
equalDOF 1292 1291  1 2
equalDOF 1355 1354  1 2
equalDOF 1420 1421  1 2
equalDOF 1486 1487  1 2
equalDOF 1550 1552  1 2
equalDOF 1622 1621  1 2
equalDOF 1689 1688  1 2
equalDOF 1763 1764  1 2
equalDOF 1836 1835  1 2
equalDOF 1909 1908  1 2
equalDOF 1982 1983  1 2
equalDOF 2065 2064  1 2
equalDOF 2138 2139  1 2
equalDOF 2223 2221  1 2
equalDOF 2305 2304  1 2
equalDOF 2392 2393  1 2
equalDOF 2474 2475  1 2
equalDOF 2561 2560  1 2
equalDOF 2655 2656  1 2
equalDOF 2743 2742  1 2
equalDOF 2838 2837  1 2
equalDOF 2930 2929  1 2
equalDOF 3024 3023  1 2
equalDOF 3123 3124  1 2
equalDOF 3221 3222  1 2
equalDOF 3322 3321  1 2
equalDOF 3421 3422  1 2
equalDOF 3526 3524  1 2
equalDOF 3636 3637  1 2
equalDOF 3741 3742  1 2
equalDOF 3849 3850  1 2
equalDOF 3956 3957  1 2
equalDOF 4063 4064  1 2
equalDOF 4186 4185  1 2
equalDOF 4293 4294  1 2
equalDOF 4412 4411  1 2
equalDOF 4529 4527  1 2
equalDOF 4646 4645  1 2
equalDOF 4771 4769  1 2
equalDOF 4889 4890  1 2
equalDOF 5015 5014  1 2
equalDOF 5141 5140  1 2
equalDOF 5268 5264  1 2
equalDOF 5402 5401  1 2
equalDOF 5523 5522  1 2
equalDOF 5657 5656  1 2
equalDOF 5791 5790  1 2
equalDOF 5924 5923  1 2
equalDOF 6063 6062  1 2
equalDOF 6194 6195  1 2
equalDOF 6339 6340  1 2
equalDOF 6484 6483  1 2
equalDOF 6625 6624  1 2
equalDOF 6767 6768  1 2
equalDOF 6909 6908  1 2
equalDOF 7065 7063  1 2
equalDOF 7207 7206  1 2
equalDOF 7362 7360  1 2
equalDOF 7510 7509  1 2
equalDOF 7664 7665  1 2
equalDOF 7829 7830  1 2
equalDOF 7981 7980  1 2
equalDOF 8138 8139  1 2
equalDOF 8150 8149  1 2
equalDOF 8175 8173  1 2
equalDOF 8186 8187  1 2
equalDOF 8206 8205  1 2
equalDOF 8237 8236  1 2
equalDOF 8275 8273  1 2
equalDOF 8295 8294  1 2
equalDOF 8321 8323  1 2
equalDOF 8341 8342  1 2
equalDOF 8377 8374  1 2
equalDOF 8422 8421  1 2
equalDOF 8444 8445  1 2
equalDOF 8489 8488  1 2
equalDOF 8520 8521  1 2
equalDOF 8560 8561  1 2
equalDOF 8608 8607  1 2
equalDOF 8646 8647  1 2
equalDOF 8692 8691  1 2
equalDOF 8741 8742  1 2
equalDOF 8794 8793  1 2
equalDOF 8837 8836  1 2
equalDOF 8876 8879  1 2
equalDOF 8937 8935  1 2
equalDOF 8994 8995  1 2
equalDOF 9039 9038  1 2
equalDOF 9111 9109  1 2
equalDOF 9168 9166  1 2
equalDOF 9223 9221  1 2
equalDOF 9287 9286  1 2
equalDOF 9349 9350  1 2
equalDOF 9416 9417  1 2
equalDOF 9484 9482  1 2
equalDOF 9548 9545  1 2
equalDOF 9625 9626  1 2
equalDOF 9697 9692  1 2
equalDOF 9766 9765  1 2
equalDOF 9846 9844  1 2
equalDOF 9906 9905  1 2
equalDOF 9983 9984  1 2
equalDOF 10049 10053  1 2
equalDOF 10122 10123  1 2
equalDOF 10202 10201  1 2
equalDOF 10273 10274  1 2
equalDOF 10348 10346  1 2
equalDOF 10423 10422  1 2
equalDOF 10496 10497  1 2
equalDOF 10562 10563  1 2
equalDOF 10642 10643  1 2
equalDOF 10721 10722  1 2
equalDOF 10790 10789  1 2
equalDOF 10866 10867  1 2
equalDOF 10950 10949  1 2
equalDOF 11028 11027  1 2
equalDOF 11104 11103  1 2
equalDOF 11179 11181  1 2
equalDOF 11254 11251  1 2
equalDOF 11342 11341  1 2
equalDOF 11411 11410  1 2
equalDOF 11497 11495  1 2
equalDOF 11574 11576  1 2
equalDOF 11659 11661  1 2
equalDOF 11742 11743  1 2
equalDOF 11819 11818  1 2
equalDOF 11901 11900  1 2
equalDOF 11981 11980  1 2
equalDOF 12066 12067  1 2
equalDOF 12153 12152  1 2
equalDOF 12233 12235  1 2
equalDOF 12320 12321  1 2
equalDOF 12402 12403  1 2
equalDOF 12480 12479  1 2
equalDOF 12571 12574  1 2
equalDOF 12650 12651  1 2
equalDOF 12740 12738  1 2
equalDOF 12834 12836  1 2
equalDOF 12913 12914  1 2
equalDOF 13004 13005  1 2
equalDOF 13093 13092  1 2
equalDOF 13179 13180  1 2
equalDOF 13270 13268  1 2
equalDOF 13358 13360  1 2
equalDOF 13440 13441  1 2
equalDOF 13534 13532  1 2
equalDOF 13622 13620  1 2
equalDOF 13711 13712  1 2
equalDOF 13799 13798  1 2
equalDOF 13894 13895  1 2
equalDOF 13976 13977  1 2
equalDOF 14065 14066  1 2
equalDOF 14161 14162  1 2
equalDOF 14251 14252  1 2
equalDOF 14340 14341  1 2
equalDOF 14432 14434  1 2
equalDOF 14531 14530  1 2
equalDOF 14616 14615  1 2
equalDOF 14712 14711  1 2
equalDOF 14807 14808  1 2
equalDOF 14902 14901  1 2
equalDOF 14991 14990  1 2
equalDOF 15089 15088  1 2
equalDOF 15180 15182  1 2
equalDOF 15275 15276  1 2
equalDOF 15368 15369  1 2
equalDOF 15463 15462  1 2
equalDOF 15556 15555  1 2
equalDOF 15652 15651  1 2
equalDOF 15742 15741  1 2
equalDOF 15838 15837  1 2
equalDOF 15938 15939  1 2
equalDOF 16033 16032  1 2
equalDOF 16120 16121  1 2
equalDOF 16222 16221  1 2
equalDOF 16324 16323  1 2
equalDOF 16418 16419  1 2
equalDOF 16510 16511  1 2
equalDOF 16607 16608  1 2
equalDOF 16711 16708  1 2
equalDOF 16806 16807  1 2
equalDOF 16900 16899  1 2
equalDOF 16996 16998  1 2
equalDOF 17097 17096  1 2
equalDOF 17194 17193  1 2
equalDOF 17293 17294  1 2
equalDOF 17390 17389  1 2
equalDOF 17485 17484  1 2
equalDOF 17584 17583  1 2
equalDOF 17685 17686  1 2
equalDOF 17784 17782  1 2
equalDOF 17884 17883  1 2
equalDOF 17983 17981  1 2
equalDOF 18078 18079  1 2
equalDOF 18183 18184  1 2
equalDOF 18282 18281  1 2
equalDOF 18379 18378  1 2
equalDOF 18478 18477  1 2
equalDOF 18574 18573  1 2
equalDOF 18675 18676  1 2
equalDOF 18776 18777  1 2
equalDOF 18879 18881  1 2
equalDOF 18981 18979  1 2
equalDOF 19077 19078  1 2
equalDOF 19173 19174  1 2
equalDOF 19275 19274  1 2
equalDOF 19382 19381  1 2
equalDOF 19484 19486  1 2
equalDOF 19582 19583  1 2
equalDOF 19689 19690  1 2
equalDOF 19783 19782  1 2
equalDOF 19883 19882  1 2
equalDOF 19985 19986  1 2
equalDOF 20092 20091  1 2
equalDOF 20190 20191  1 2
equalDOF 20297 20296  1 2
equalDOF 20397 20396  1 2
equalDOF 20496 20494  1 2
equalDOF 20597 20598  1 2
equalDOF 20703 20702  1 2
equalDOF 20806 20807  1 2
equalDOF 20908 20909  1 2
equalDOF 21010 21009  1 2
equalDOF 21109 21110  1 2
equalDOF 21209 21210  1 2
equalDOF 21312 21311  1 2
equalDOF 21418 21419  1 2
equalDOF 21524 21522  1 2
equalDOF 21628 21629  1 2
equalDOF 21729 21728  1 2
equalDOF 21834 21833  1 2
equalDOF 21936 21937  1 2
equalDOF 22035 22036  1 2
equalDOF 22141 22140  1 2
equalDOF 22247 22246  1 2
equalDOF 22350 22349  1 2
equalDOF 22457 22456  1 2
equalDOF 22517 22516  1 2
equalDOF 22519 22518  1 2
equalDOF 22522 22520  1 2
equalDOF 22526 22525  1 2
equalDOF 22527 22528  1 2
equalDOF 22534 22533  1 2
equalDOF 22536 22537  1 2
equalDOF 22545 22544  1 2
equalDOF 22547 22546  1 2
equalDOF 22553 22552  1 2
equalDOF 22559 22560  1 2
equalDOF 22568 22567  1 2
equalDOF 22572 22571  1 2
equalDOF 22578 22577  1 2
equalDOF 22586 22585  1 2
equalDOF 22594 22593  1 2
equalDOF 22599 22600  1 2
equalDOF 22606 22605  1 2
equalDOF 22621 22620  1 2
equalDOF 22627 22628  1 2
equalDOF 22640 22639  1 2
equalDOF 22648 22647  1 2
equalDOF 22665 22666  1 2
equalDOF 22681 22680  1 2
equalDOF 22685 22686  1 2
equalDOF 22695 22696  1 2
equalDOF 22710 22709  1 2
equalDOF 22721 22720  1 2
equalDOF 22738 22737  1 2
equalDOF 22752 22751  1 2
equalDOF 22766 22765  1 2
equalDOF 22782 22781  1 2
equalDOF 22796 22795  1 2
equalDOF 22806 22804  1 2
equalDOF 22825 22826  1 2
equalDOF 22837 22836  1 2
equalDOF 22859 22860  1 2
equalDOF 22878 22877  1 2
equalDOF 22895 22894  1 2
equalDOF 22901 22899  1 2
equalDOF 22921 22922  1 2
equalDOF 22941 22942  1 2
equalDOF 22963 22962  1 2
equalDOF 22985 22986  1 2
equalDOF 23004 23005  1 2
equalDOF 23012 23011  1 2
equalDOF 23028 23027  1 2
equalDOF 23052 23053  1 2
equalDOF 23079 23080  1 2
equalDOF 23104 23103  1 2
equalDOF 23122 23121  1 2
equalDOF 23131 23130  1 2
equalDOF 23150 23151  1 2
equalDOF 23184 23183  1 2
equalDOF 23211 23210  1 2
equalDOF 23230 23229  1 2
equalDOF 23235 23236  1 2
equalDOF 23260 23262  1 2
equalDOF 23294 23293  1 2
equalDOF 23323 23322  1 2
equalDOF 23339 23338  1 2
equalDOF 23350 23351  1 2
equalDOF 23376 23375  1 2
equalDOF 23409 23410  1 2
equalDOF 23442 23443  1 2
equalDOF 23450 23449  1 2
equalDOF 23476 23475  1 2
equalDOF 23504 23503  1 2
equalDOF 23539 23538  1 2
equalDOF 23554 23555  1 2
equalDOF 23574 23575  1 2
equalDOF 23604 23602  1 2
equalDOF 23635 23634  1 2
equalDOF 23661 23660  1 2
equalDOF 23677 23676  1 2
equalDOF 23707 23708  1 2
equalDOF 23744 23742  1 2
equalDOF 23765 23763  1 2
equalDOF 23782 23783  1 2
equalDOF 23815 23814  1 2
equalDOF 23851 23849  1 2
equalDOF 23871 23870  1 2
equalDOF 23884 23885  1 2
equalDOF 23913 23914  1 2
equalDOF 23946 23947  1 2
equalDOF 23959 23960  1 2
equalDOF 23979 23980  1 2
equalDOF 24011 24010  1 2
equalDOF 24038 24037  1 2
equalDOF 24044 24043  1 2
equalDOF 24070 24068  1 2
equalDOF 24097 24098  1 2
equalDOF 24117 24118  1 2
equalDOF 24129 24128  1 2
equalDOF 24152 24153  1 2
equalDOF 24181 24182  1 2
equalDOF 24186 24187  1 2
equalDOF 24208 24209  1 2
equalDOF 24235 24236  1 2
equalDOF 24254 24255  1 2
equalDOF 24261 24260  1 2
equalDOF 24284 24285  1 2
equalDOF 24313 24312  1 2
equalDOF 24319 24320  1 2
equalDOF 24340 24339  1 2
equalDOF 24364 24363  1 2
equalDOF 24378 24379  1 2
equalDOF 24388 24389  1 2
equalDOF 24410 24411  1 2
equalDOF 24434 24433  1 2
equalDOF 24435 24437  1 2
equalDOF 24458 24457  1 2
equalDOF 24480 24479  1 2
equalDOF 24487 24488  1 2
equalDOF 24501 24502  1 2
equalDOF 24526 24525  1 2
equalDOF 24533 24534  1 2
equalDOF 24547 24548  1 2
equalDOF 24569 24568  1 2
equalDOF 24582 24581  1 2
equalDOF 24592 24591  1 2
equalDOF 24610 24611  1 2
equalDOF 24626 24625  1 2
equalDOF 24629 24628  1 2
equalDOF 24650 24649  1 2
equalDOF 24667 24668  1 2
puts "Finished creating equalDOF constraints for soil columns..."

puts "Finished creating equalDOF for base..."

#-----------------------------------------------------------------------------------------
#  2. CREATE SOIL MATERIALS
#-----------------------------------------------------------------------------------------

source material.tcl
puts "Finished creating all soil materials..."

#-----------------------------------------------------------------------------------------
#  3. CREATE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------

source elementInfo.dat
source PMLelementInfo.dat

set dT 0.001
timeSeries Path 1 -dt 0.001 -filePath force.dat -factor 0.004
pattern Plain 1 1 {
load 9741 0.0 -1.0
}

#-----------------------------------------------------------------------------------------
#  4. CREATE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------

# record nodal displacement
# recorder Node -file displacement.out -time -nodeRange 1 24932 -dof 1 2 disp
# recorder Node -file acceleration.out -time -nodeRange 1 24932 -dof 1 2 accel

recorder Node -file displacement_select.out -time -node 9741 15088 22516 24668 -dof 1 2 disp

puts "Finished creating recorders..."


constraints Plain
numberer Plain
integrator Newmark 0.5 0.25
system UmfPack
test NormDispIncr 1.0e-6 10 1
algorithm ModifiedNewton -factoronce
analysis Transient

set ok [analyze 2000 0.001]
wipe

