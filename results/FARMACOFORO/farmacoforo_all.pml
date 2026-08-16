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
pseudoatom feature_0, pos=[13.588,43.338,5.593]
hide nonbonded, feature_0
show spheres, feature_0
color blue, feature_0
alter feature_0, vdw=0.445
set sphere_transparency, 0.15, feature_0
pseudoatom feature_1, pos=[14.976,44.86,4.775]
hide nonbonded, feature_1
show spheres, feature_1
color red, feature_1
alter feature_1, vdw=0.297
set sphere_transparency, 0.15, feature_1
pseudoatom feature_2, pos=[14.671,42.225,7.373]
hide nonbonded, feature_2
show spheres, feature_2
color forest, feature_2
alter feature_2, vdw=1.774
set sphere_transparency, 0.15, feature_2
pseudoatom feature_3, pos=[14.707,41.61,8.927]
hide nonbonded, feature_3
show spheres, feature_3
color yellow, feature_3
alter feature_3, vdw=2.545
set sphere_transparency, 0.35, feature_3
pseudoatom feature_4, pos=[15.076,41.54,7.671]
hide nonbonded, feature_4
show spheres, feature_4
color yellow, feature_4
alter feature_4, vdw=2.788
set sphere_transparency, 0.35, feature_4
pseudoatom feature_5, pos=[13.348,39.246,11.268]
hide nonbonded, feature_5
show spheres, feature_5
color yellow, feature_5
alter feature_5, vdw=1.338
set sphere_transparency, 0.35, feature_5
pseudoatom feature_6, pos=[16.847,44.389,6.672]
hide nonbonded, feature_6
show spheres, feature_6
color yellow, feature_6
alter feature_6, vdw=0.204
set sphere_transparency, 0.65, feature_6
pseudoatom feature_7, pos=[13.457,39.371,10.175]
hide nonbonded, feature_7
show spheres, feature_7
color forest, feature_7
alter feature_7, vdw=1.464
set sphere_transparency, 0.65, feature_7
pseudoatom feature_8, pos=[17.858,43.368,8.318]
hide nonbonded, feature_8
show spheres, feature_8
color yellow, feature_8
alter feature_8, vdw=0.839
set sphere_transparency, 0.65, feature_8
pseudoatom feature_9, pos=[2.519,42.647,10.828]
hide nonbonded, feature_9
show spheres, feature_9
color red, feature_9
alter feature_9, vdw=3.1
set sphere_transparency, 0.65, feature_9
pseudoatom feature_10, pos=[3.108,41.242,10.56]
hide nonbonded, feature_10
show spheres, feature_10
color yellow, feature_10
alter feature_10, vdw=1.72
set sphere_transparency, 0.65, feature_10
pseudoatom feature_11, pos=[11.691,35.786,11.321]
hide nonbonded, feature_11
show spheres, feature_11
color blue, feature_11
alter feature_11, vdw=1.0
set sphere_transparency, 0.65, feature_11
pseudoatom feature_12, pos=[4.808,40.793,10.83]
hide nonbonded, feature_12
show spheres, feature_12
color yellow, feature_12
alter feature_12, vdw=1.0
set sphere_transparency, 0.65, feature_12
pseudoatom feature_13, pos=[11.094,44.085,9.958]
hide nonbonded, feature_13
show spheres, feature_13
color red, feature_13
alter feature_13, vdw=1.0
set sphere_transparency, 0.65, feature_13
pseudoatom feature_14, pos=[1.733,42.209,13.053]
hide nonbonded, feature_14
show spheres, feature_14
color blue, feature_14
alter feature_14, vdw=1.0
set sphere_transparency, 0.65, feature_14
pseudoatom feature_15, pos=[14.215,38.202,10.901]
hide nonbonded, feature_15
show spheres, feature_15
color forest, feature_15
alter feature_15, vdw=1.0
set sphere_transparency, 0.65, feature_15
rebuild
set_view (\
0.295335442,    0.495653808,    0.816764057,\
0.644749045,    0.527470946,   -0.553234816,\
-0.705034792,    0.689999223,   -0.163792834,\
0.000000000,    0.000000000,  -80.127052307,\
13.514476776,   40.709243774,    9.326225281,\
-8592.375000000, 8752.620117188,  -20.000000000 )
set internal_gui, 0