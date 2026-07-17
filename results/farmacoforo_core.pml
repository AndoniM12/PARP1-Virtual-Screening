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
pseudoatom feature_0, pos=[14.312,44.108,5.27]
hide nonbonded, feature_0
show spheres, feature_0
color blue, feature_0
alter feature_0, vdw=1.209
set sphere_transparency, 0.15, feature_0
pseudoatom feature_1, pos=[15.245,42.208,7.686]
hide nonbonded, feature_1
show spheres, feature_1
color forest, feature_1
alter feature_1, vdw=1.456
set sphere_transparency, 0.15, feature_1
rebuild
set_view (\
0.295335442,    0.495653808,    0.816764057,\
0.644749045,    0.527470946,   -0.553234816,\
-0.705034792,    0.689999223,   -0.163792834,\
0.000000000,    0.000000000,  -80.127052307,\
13.514476776,   40.709243774,    9.326225281,\
-8592.375000000, 8752.620117188,  -20.000000000 )
set internal_gui, 0