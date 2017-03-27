      4
      1   -5.9812   12.5978   -4.5636
      2    0.8677   -0.9451   12.6382
      3   11.6103    1.0457   -4.1139
      4   -6.4967  -12.6984   -3.9607

      4
      3               2      1      4
      3               2      4      3
      3               2      3      1
      3               4      1      3

      6
       2       2       1
       2       1       4
       2       4       2
       2       4       3
       2       3       2
       2       3       1


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =    4
    geom.n_face =    4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [    -5.9812d0,    12.5978d0,    -4.5636d0 ]
    geom.iniP(      2).pos(1:3) = [     0.8677d0,    -0.9451d0,    12.6382d0 ]
    geom.iniP(      3).pos(1:3) = [    11.6103d0,     1.0457d0,    -4.1139d0 ]
    geom.iniP(      4).pos(1:3) = [    -6.4967d0,   -12.6984d0,    -3.9607d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [      2,       1,       4 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [      2,       4,       3 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [      2,       3,       1 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [      4,       1,       3 ]

--------------------------------------------------

      6
line(      1, 1:2) = [       2,       1 ]
line(      2, 1:2) = [       1,       4 ]
line(      3, 1:2) = [       4,       2 ]
line(      4, 1:2) = [       4,       3 ]
line(      5, 1:2) = [       3,       2 ]
line(      6, 1:2) = [       3,       1 ]

