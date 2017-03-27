      6
      1    1.4393   -1.8152   18.6786
      2   -2.5337  -15.4490    0.6888
      3   15.9284   -0.7617   -1.4565
      4    1.1919   -0.7591  -19.0232
      5   -0.9576   16.1710    2.6237
      6  -15.0682    2.6140   -1.5113

      8
      3               1      3      5
      3               1      2      3
      3               1      6      2
      3               1      5      6
      3               4      5      3
      3               4      3      2
      3               4      2      6
      3               4      6      5

     12
       2       1       3
       2       3       5
       2       5       1
       2       1       2
       2       2       3
       2       1       6
       2       6       2
       2       5       6
       2       4       5
       2       3       4
       2       2       4
       2       6       4


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =    6
    geom.n_face =    8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [     1.4393d0,    -1.8152d0,    18.6786d0 ]
    geom.iniP(      2).pos(1:3) = [    -2.5337d0,   -15.4490d0,     0.6888d0 ]
    geom.iniP(      3).pos(1:3) = [    15.9284d0,    -0.7617d0,    -1.4565d0 ]
    geom.iniP(      4).pos(1:3) = [     1.1919d0,    -0.7591d0,   -19.0232d0 ]
    geom.iniP(      5).pos(1:3) = [    -0.9576d0,    16.1710d0,     2.6237d0 ]
    geom.iniP(      6).pos(1:3) = [   -15.0682d0,     2.6140d0,    -1.5113d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [      1,       3,       5 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [      1,       2,       3 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [      1,       6,       2 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [      1,       5,       6 ]
    geom.face(      5).n_poi =       3; allocate(geom.face(      5).poi(      3)); geom.face(      5).poi(1:      3) = [      4,       5,       3 ]
    geom.face(      6).n_poi =       3; allocate(geom.face(      6).poi(      3)); geom.face(      6).poi(1:      3) = [      4,       3,       2 ]
    geom.face(      7).n_poi =       3; allocate(geom.face(      7).poi(      3)); geom.face(      7).poi(1:      3) = [      4,       2,       6 ]
    geom.face(      8).n_poi =       3; allocate(geom.face(      8).poi(      3)); geom.face(      8).poi(1:      3) = [      4,       6,       5 ]

--------------------------------------------------

     12
line(      1, 1:2) = [       1,       3 ]
line(      2, 1:2) = [       3,       5 ]
line(      3, 1:2) = [       5,       1 ]
line(      4, 1:2) = [       1,       2 ]
line(      5, 1:2) = [       2,       3 ]
line(      6, 1:2) = [       1,       6 ]
line(      7, 1:2) = [       6,       2 ]
line(      8, 1:2) = [       5,       6 ]
line(      9, 1:2) = [       4,       5 ]
line(     10, 1:2) = [       3,       4 ]
line(     11, 1:2) = [       2,       4 ]
line(     12, 1:2) = [       6,       4 ]

