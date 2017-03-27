      8
      1  -10.8699   -9.0403   12.8253
      2   13.1783    9.3523    9.2211
      3   -9.4579   13.0216   10.4285
      4  -10.6283  -11.0596   -8.5752
      5   10.3701   11.8547  -10.4220
      6    8.6643  -13.3699    9.2543
      7   11.1136  -12.4990  -11.8562
      8  -12.3703   11.7402  -10.8758

     12
      3               6      2      1
      3               1      2      3
      3               2      5      3
      3               3      5      8
      3               8      5      4
      3               4      5      7
      3               7      6      4
      3               4      6      1
      3               6      7      2
      3               2      7      5
      3               1      3      4
      3               4      3      8

     18
       2       6       2
       2       2       1
       2       1       6
       2       2       3
       2       3       1
       2       2       5
       2       5       3
       2       5       8
       2       8       3
       2       5       4
       2       4       8
       2       5       7
       2       7       4
       2       7       6
       2       6       4
       2       1       4
       2       7       2
       2       3       4


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =    8
    geom.n_face =   12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [   -10.8699d0,    -9.0403d0,    12.8253d0 ]
    geom.iniP(      2).pos(1:3) = [    13.1783d0,     9.3523d0,     9.2211d0 ]
    geom.iniP(      3).pos(1:3) = [    -9.4579d0,    13.0216d0,    10.4285d0 ]
    geom.iniP(      4).pos(1:3) = [   -10.6283d0,   -11.0596d0,    -8.5752d0 ]
    geom.iniP(      5).pos(1:3) = [    10.3701d0,    11.8547d0,   -10.4220d0 ]
    geom.iniP(      6).pos(1:3) = [     8.6643d0,   -13.3699d0,     9.2543d0 ]
    geom.iniP(      7).pos(1:3) = [    11.1136d0,   -12.4990d0,   -11.8562d0 ]
    geom.iniP(      8).pos(1:3) = [   -12.3703d0,    11.7402d0,   -10.8758d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [      6,       2,       1 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [      1,       2,       3 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [      2,       5,       3 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [      3,       5,       8 ]
    geom.face(      5).n_poi =       3; allocate(geom.face(      5).poi(      3)); geom.face(      5).poi(1:      3) = [      8,       5,       4 ]
    geom.face(      6).n_poi =       3; allocate(geom.face(      6).poi(      3)); geom.face(      6).poi(1:      3) = [      4,       5,       7 ]
    geom.face(      7).n_poi =       3; allocate(geom.face(      7).poi(      3)); geom.face(      7).poi(1:      3) = [      7,       6,       4 ]
    geom.face(      8).n_poi =       3; allocate(geom.face(      8).poi(      3)); geom.face(      8).poi(1:      3) = [      4,       6,       1 ]
    geom.face(      9).n_poi =       3; allocate(geom.face(      9).poi(      3)); geom.face(      9).poi(1:      3) = [      6,       7,       2 ]
    geom.face(     10).n_poi =       3; allocate(geom.face(     10).poi(      3)); geom.face(     10).poi(1:      3) = [      2,       7,       5 ]
    geom.face(     11).n_poi =       3; allocate(geom.face(     11).poi(      3)); geom.face(     11).poi(1:      3) = [      1,       3,       4 ]
    geom.face(     12).n_poi =       3; allocate(geom.face(     12).poi(      3)); geom.face(     12).poi(1:      3) = [      4,       3,       8 ]

--------------------------------------------------

     18
line(      1, 1:2) = [       6,       2 ]
line(      2, 1:2) = [       2,       1 ]
line(      3, 1:2) = [       1,       6 ]
line(      4, 1:2) = [       2,       3 ]
line(      5, 1:2) = [       3,       1 ]
line(      6, 1:2) = [       2,       5 ]
line(      7, 1:2) = [       5,       3 ]
line(      8, 1:2) = [       5,       8 ]
line(      9, 1:2) = [       8,       3 ]
line(     10, 1:2) = [       5,       4 ]
line(     11, 1:2) = [       4,       8 ]
line(     12, 1:2) = [       5,       7 ]
line(     13, 1:2) = [       7,       4 ]
line(     14, 1:2) = [       7,       6 ]
line(     15, 1:2) = [       6,       4 ]
line(     16, 1:2) = [       1,       4 ]
line(     17, 1:2) = [       7,       2 ]
line(     18, 1:2) = [       3,       4 ]

