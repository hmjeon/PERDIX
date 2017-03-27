     12
      1   12.8559  -18.4078    4.5927
      2    2.6353  -18.2908  -18.0086
      3   19.4301   -1.3006  -15.5467
      4   20.6532    2.9721   11.8687
      5  -19.4518    1.5734   10.8306
      6  -21.3999    1.9363  -13.0532
      7    0.3631  -12.8147   19.8932
      8    2.7443    8.5393   18.8170
      9   -1.9827   14.4621  -20.2903
     10  -14.9666  -19.2470    2.2430
     11   11.3689   20.2643   -0.3120
     12  -12.2499   20.3134   -1.0344

     20
      3               4      1      3
      3               3     11      4
      3               3      2      9
      3               3      9     11
      3               3      1      2
      3               7      4      8
      3               1      4      7
      3               4     11      8
      3               6     10      5
      3               5     12      6
      3               7      8      5
      3               5     10      7
      3               8     12      5
      3               9      2      6
      3               6     12      9
      3               6      2     10
      3               7     10      1
      3               8     11     12
      3               9     12     11
      3               2      1     10

     30
       2       4       1
       2       1       3
       2       3       4
       2       3      11
       2      11       4
       2       3       2
       2       2       9
       2       9       3
       2       9      11
       2       1       2
       2       7       4
       2       4       8
       2       8       7
       2       7       1
       2      11       8
       2       6      10
       2      10       5
       2       5       6
       2       5      12
       2      12       6
       2       8       5
       2       5       7
       2      10       7
       2       8      12
       2       2       6
       2       6       9
       2      12       9
       2       2      10
       2      10       1
       2      11      12


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =   12
    geom.n_face =   20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [    12.8559d0,   -18.4078d0,     4.5927d0 ]
    geom.iniP(      2).pos(1:3) = [     2.6353d0,   -18.2908d0,   -18.0086d0 ]
    geom.iniP(      3).pos(1:3) = [    19.4301d0,    -1.3006d0,   -15.5467d0 ]
    geom.iniP(      4).pos(1:3) = [    20.6532d0,     2.9721d0,    11.8687d0 ]
    geom.iniP(      5).pos(1:3) = [   -19.4518d0,     1.5734d0,    10.8306d0 ]
    geom.iniP(      6).pos(1:3) = [   -21.3999d0,     1.9363d0,   -13.0532d0 ]
    geom.iniP(      7).pos(1:3) = [     0.3631d0,   -12.8147d0,    19.8932d0 ]
    geom.iniP(      8).pos(1:3) = [     2.7443d0,     8.5393d0,    18.8170d0 ]
    geom.iniP(      9).pos(1:3) = [    -1.9827d0,    14.4621d0,   -20.2903d0 ]
    geom.iniP(     10).pos(1:3) = [   -14.9666d0,   -19.2470d0,     2.2430d0 ]
    geom.iniP(     11).pos(1:3) = [    11.3689d0,    20.2643d0,    -0.3120d0 ]
    geom.iniP(     12).pos(1:3) = [   -12.2499d0,    20.3134d0,    -1.0344d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [      4,       1,       3 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [      3,      11,       4 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [      3,       2,       9 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [      3,       9,      11 ]
    geom.face(      5).n_poi =       3; allocate(geom.face(      5).poi(      3)); geom.face(      5).poi(1:      3) = [      3,       1,       2 ]
    geom.face(      6).n_poi =       3; allocate(geom.face(      6).poi(      3)); geom.face(      6).poi(1:      3) = [      7,       4,       8 ]
    geom.face(      7).n_poi =       3; allocate(geom.face(      7).poi(      3)); geom.face(      7).poi(1:      3) = [      1,       4,       7 ]
    geom.face(      8).n_poi =       3; allocate(geom.face(      8).poi(      3)); geom.face(      8).poi(1:      3) = [      4,      11,       8 ]
    geom.face(      9).n_poi =       3; allocate(geom.face(      9).poi(      3)); geom.face(      9).poi(1:      3) = [      6,      10,       5 ]
    geom.face(     10).n_poi =       3; allocate(geom.face(     10).poi(      3)); geom.face(     10).poi(1:      3) = [      5,      12,       6 ]
    geom.face(     11).n_poi =       3; allocate(geom.face(     11).poi(      3)); geom.face(     11).poi(1:      3) = [      7,       8,       5 ]
    geom.face(     12).n_poi =       3; allocate(geom.face(     12).poi(      3)); geom.face(     12).poi(1:      3) = [      5,      10,       7 ]
    geom.face(     13).n_poi =       3; allocate(geom.face(     13).poi(      3)); geom.face(     13).poi(1:      3) = [      8,      12,       5 ]
    geom.face(     14).n_poi =       3; allocate(geom.face(     14).poi(      3)); geom.face(     14).poi(1:      3) = [      9,       2,       6 ]
    geom.face(     15).n_poi =       3; allocate(geom.face(     15).poi(      3)); geom.face(     15).poi(1:      3) = [      6,      12,       9 ]
    geom.face(     16).n_poi =       3; allocate(geom.face(     16).poi(      3)); geom.face(     16).poi(1:      3) = [      6,       2,      10 ]
    geom.face(     17).n_poi =       3; allocate(geom.face(     17).poi(      3)); geom.face(     17).poi(1:      3) = [      7,      10,       1 ]
    geom.face(     18).n_poi =       3; allocate(geom.face(     18).poi(      3)); geom.face(     18).poi(1:      3) = [      8,      11,      12 ]
    geom.face(     19).n_poi =       3; allocate(geom.face(     19).poi(      3)); geom.face(     19).poi(1:      3) = [      9,      12,      11 ]
    geom.face(     20).n_poi =       3; allocate(geom.face(     20).poi(      3)); geom.face(     20).poi(1:      3) = [      2,       1,      10 ]

--------------------------------------------------

     30
line(      1, 1:2) = [       4,       1 ]
line(      2, 1:2) = [       1,       3 ]
line(      3, 1:2) = [       3,       4 ]
line(      4, 1:2) = [       3,      11 ]
line(      5, 1:2) = [      11,       4 ]
line(      6, 1:2) = [       3,       2 ]
line(      7, 1:2) = [       2,       9 ]
line(      8, 1:2) = [       9,       3 ]
line(      9, 1:2) = [       9,      11 ]
line(     10, 1:2) = [       1,       2 ]
line(     11, 1:2) = [       7,       4 ]
line(     12, 1:2) = [       4,       8 ]
line(     13, 1:2) = [       8,       7 ]
line(     14, 1:2) = [       7,       1 ]
line(     15, 1:2) = [      11,       8 ]
line(     16, 1:2) = [       6,      10 ]
line(     17, 1:2) = [      10,       5 ]
line(     18, 1:2) = [       5,       6 ]
line(     19, 1:2) = [       5,      12 ]
line(     20, 1:2) = [      12,       6 ]
line(     21, 1:2) = [       8,       5 ]
line(     22, 1:2) = [       5,       7 ]
line(     23, 1:2) = [      10,       7 ]
line(     24, 1:2) = [       8,      12 ]
line(     25, 1:2) = [       2,       6 ]
line(     26, 1:2) = [       6,       9 ]
line(     27, 1:2) = [      12,       9 ]
line(     28, 1:2) = [       2,      10 ]
line(     29, 1:2) = [      10,       1 ]
line(     30, 1:2) = [      11,      12 ]

