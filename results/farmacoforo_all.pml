load structures/processed/7AAA_chainA.pdb, protein
hide everything
show cartoon, protein
color gray90, protein
set cartoon_transparency, 0.65
bg_color white
select pocket, protein and resi 863+889+896+907
show sticks, pocket
color gray70, pocket and elem C
color blue, pocket and elem N
color red, pocket and elem O
color yellow, pocket and elem S
label pocket and name CA, resn + resi
set label_size, 18
set label_color, black
set label_font_id, 7
set label_position, (2,1,1)
pseudoatom feature_0, pos=[14.31,44.11,5.27]
hide nonbonded, feature_0
show spheres, feature_0
color blue, feature_0
alter feature_0, vdw=1.21
set sphere_transparency, 0.15, feature_0
pseudoatom feature_1, pos=[15.24,42.21,7.69]
hide nonbonded, feature_1
show spheres, feature_1
color forest, feature_1
alter feature_1, vdw=1.46
set sphere_transparency, 0.15, feature_1
pseudoatom feature_2, pos=[14.68,41.42,9.2]
hide nonbonded, feature_2
show spheres, feature_2
color yellow, feature_2
alter feature_2, vdw=2.62
set sphere_transparency, 0.35, feature_2
pseudoatom feature_3, pos=[14.9,41.23,7.8]
hide nonbonded, feature_3
show spheres, feature_3
color yellow, feature_3
alter feature_3, vdw=2.67
set sphere_transparency, 0.35, feature_3
pseudoatom feature_4, pos=[13.35,39.25,11.27]
hide nonbonded, feature_4
show spheres, feature_4
color yellow, feature_4
alter feature_4, vdw=1.34
set sphere_transparency, 0.35, feature_4
pseudoatom feature_5, pos=[17.86,43.37,8.32]
hide nonbonded, feature_5
show spheres, feature_5
color yellow, feature_5
alter feature_5, vdw=0.84
set sphere_transparency, 0.65, feature_5
pseudoatom feature_6, pos=[13.46,39.37,10.17]
hide nonbonded, feature_6
show spheres, feature_6
color forest, feature_6
alter feature_6, vdw=1.46
set sphere_transparency, 0.65, feature_6
pseudoatom feature_7, pos=[16.92,44.43,6.73]
hide nonbonded, feature_7
show spheres, feature_7
color yellow, feature_7
alter feature_7, vdw=0.14
set sphere_transparency, 0.65, feature_7
pseudoatom feature_8, pos=[11.69,35.79,11.32]
hide nonbonded, feature_8
show spheres, feature_8
color blue, feature_8
alter feature_8, vdw=1.0
set sphere_transparency, 0.65, feature_8
pseudoatom feature_9, pos=[11.09,44.08,9.96]
hide nonbonded, feature_9
show spheres, feature_9
color blue, feature_9
alter feature_9, vdw=1.0
set sphere_transparency, 0.65, feature_9
pseudoatom feature_10, pos=[3.57,40.79,11.59]
hide nonbonded, feature_10
show spheres, feature_10
color yellow, feature_10
alter feature_10, vdw=1.0
set sphere_transparency, 0.65, feature_10
pseudoatom feature_11, pos=[4.62,43.12,10.44]
hide nonbonded, feature_11
show spheres, feature_11
color blue, feature_11
alter feature_11, vdw=1.0
set sphere_transparency, 0.65, feature_11
pseudoatom feature_12, pos=[4.81,40.79,10.83]
hide nonbonded, feature_12
show spheres, feature_12
color yellow, feature_12
alter feature_12, vdw=1.0
set sphere_transparency, 0.65, feature_12
pseudoatom feature_13, pos=[14.22,38.2,10.9]
hide nonbonded, feature_13
show spheres, feature_13
color forest, feature_13
alter feature_13, vdw=1.0
set sphere_transparency, 0.65, feature_13
rebuild
set_view (\
0.295335442,    0.495653808,    0.816764057,\
0.644749045,    0.527470946,   -0.553234816,\
-0.705034792,    0.689999223,   -0.163792834,\
0.000000000,    0.000000000,  -80.127052307,\
13.514476776,   40.709243774,    9.326225281,\
-8592.375000000, 8752.620117188,  -20.000000000 )
set internal_gui, 0