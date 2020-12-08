package require psfgen 
topology top_all36_carb.rtf
set name  hemi 


resetpsf
segment  he {
 residue 1 BXYL
 residue 2 BXYL
 residue 3 BXYL
 residue 4 BXYL
 residue 5 BXYL 
 residue 6 BXYL
 residue 7 BMAN
 residue 8 BXYL 
 residue 33 ALCA ;#
 residue 9 BXYL 
 residue 10 BXYL
 residue 11 BXYL
 residue 12 BXYL
 residue 34 ALCA ;#
 residue 13 BXYL
 residue 14 BMAN
 residue 15 BXYL 
 residue 16 BXYL
 residue 17 BXYL
 residue 18 BXYL 
 residue 19 BXYL 
 residue 20 BXYL
 residue 35 ALCA ;#
 residue 21 BXYL
 residue 22 BXYL
 residue 23 BXYL
 residue 24 BXYL
 residue 25 BXYL 
 residue 26 BXYL
 residue 27 BXYL
 residue 28 BMAN 
 residue 29 BXYL
 residue 36 ALCA ;#
 residue 30 BXYL
 residue 31 BMAN
 residue 32 BXYL
}
patch 14bb he:1 he:2
patch 14bb he:2 he:3
patch 14bb he:3 he:4
patch 14bb he:4 he:5
patch 14bb he:5 he:6
patch 14bb he:6 he:7
patch 14bb he:7 he:8 
patch 12ab he:8  he:33 
patch OME  he:33
patch 14bb he:8 he:9 
patch 14bb he:9 he:10 
patch 14bb he:10 he:11
patch 14bb he:11 he:12
patch 12ab he:12 he:34
patch OME  he:34
patch 14bb he:12 he:13
patch 14bb he:13 he:14
patch 14bb he:14 he:15
patch 14bb he:15 he:16
patch 14bb he:16 he:17
patch 14bb he:17 he:18 
patch 14bb he:18 he:19 
patch 14bb he:19 he:20
patch 12ab he:20 he:35
patch OME  he:35
patch 14bb he:20 he:21
patch 14bb he:21 he:22
patch 14bb he:22 he:23
patch 14bb he:23 he:24
patch 14bb he:24 he:25
patch 14bb he:25 he:26
patch 14bb he:26 he:27
patch 14bb he:27 he:28 
patch 14bb he:28 he:29 
patch 12ab he:29 he:36
patch OME  he:36 
patch 14bb he:29 he:30 
patch 14bb he:30 he:31
patch 14bb he:31 he:32


regenerate angles dihedrals
coordpdb hemi.res1.pdb
guesscoord 	 
writepdb $name.pdb	 
writepsf $name.psf

exit
