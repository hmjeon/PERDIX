      12
       1   -0.0000    0.0000   19.0212
       2   17.0130    0.0000    8.5065
       3    5.2573   16.1804    8.5065
       4  -13.7638   10.0000    8.5065
       5  -13.7638  -10.0000    8.5065
       6    5.2573  -16.1804    8.5065
       7   13.7638   10.0000   -8.5065
       8   13.7638  -10.0000   -8.5065
       9   -5.2573   16.1804   -8.5065
      10  -17.0130    0.0000   -8.5065
      11   -5.2573  -16.1804   -8.5065
      12   -0.0000    0.0000  -19.0212

      20
       3               1       2       3
       3               1       3       4
       3               1       4       5
       3               1       5       6
       3               1       6       2
       3               2       6       8
       3               2       8       7
       3               2       7       3
       3               3       7       9
       3               3       9       4
       3               4       9      10
       3               4      10       5
       3               5      10      11
       3               5      11       6
       3               6      11       8
       3               7       8      12
       3               7      12       9
       3               8      11      12
       3               9      12      10
       3              10      12      11

      30
       2       1       2
       2       2       3
       2       3       1
       2       3       4
       2       4       1
       2       4       5
       2       5       1
       2       5       6
       2       6       1
       2       6       2
       2       6       8
       2       8       2
       2       8       7
       2       7       2
       2       7       3
       2       7       9
       2       9       3
       2       9       4
       2       9      10
       2      10       4
       2      10       5
       2      10      11
       2      11       5
       2      11       6
       2      11       8
       2       8      12
       2      12       7
       2      12       9
       2      11      12
       2      12      10


-------------------------------------------------------------------------------

      12
point(    1).pos(1:3) = [    -0.0000d0,     0.0000d0,    19.0212d0 ]
point(    2).pos(1:3) = [    17.0130d0,     0.0000d0,     8.5065d0 ]
point(    3).pos(1:3) = [     5.2573d0,    16.1804d0,     8.5065d0 ]
point(    4).pos(1:3) = [   -13.7638d0,    10.0000d0,     8.5065d0 ]
point(    5).pos(1:3) = [   -13.7638d0,   -10.0000d0,     8.5065d0 ]
point(    6).pos(1:3) = [     5.2573d0,   -16.1804d0,     8.5065d0 ]
point(    7).pos(1:3) = [    13.7638d0,    10.0000d0,    -8.5065d0 ]
point(    8).pos(1:3) = [    13.7638d0,   -10.0000d0,    -8.5065d0 ]
point(    9).pos(1:3) = [    -5.2573d0,    16.1804d0,    -8.5065d0 ]
point(   10).pos(1:3) = [   -17.0130d0,     0.0000d0,    -8.5065d0 ]
point(   11).pos(1:3) = [    -5.2573d0,   -16.1804d0,    -8.5065d0 ]
point(   12).pos(1:3) = [    -0.0000d0,     0.0000d0,   -19.0212d0 ]


      30
line(    1, 1:2) = [     1,     2 ]
line(    2, 1:2) = [     2,     3 ]
line(    3, 1:2) = [     3,     1 ]
line(    4, 1:2) = [     3,     4 ]
line(    5, 1:2) = [     4,     1 ]
line(    6, 1:2) = [     4,     5 ]
line(    7, 1:2) = [     5,     1 ]
line(    8, 1:2) = [     5,     6 ]
line(    9, 1:2) = [     6,     1 ]
line(   10, 1:2) = [     6,     2 ]
line(   11, 1:2) = [     6,     8 ]
line(   12, 1:2) = [     8,     2 ]
line(   13, 1:2) = [     8,     7 ]
line(   14, 1:2) = [     7,     2 ]
line(   15, 1:2) = [     7,     3 ]
line(   16, 1:2) = [     7,     9 ]
line(   17, 1:2) = [     9,     3 ]
line(   18, 1:2) = [     9,     4 ]
line(   19, 1:2) = [     9,    10 ]
line(   20, 1:2) = [    10,     4 ]
line(   21, 1:2) = [    10,     5 ]
line(   22, 1:2) = [    10,    11 ]
line(   23, 1:2) = [    11,     5 ]
line(   24, 1:2) = [    11,     6 ]
line(   25, 1:2) = [    11,     8 ]
line(   26, 1:2) = [     8,    12 ]
line(   27, 1:2) = [    12,     7 ]
line(   28, 1:2) = [    12,     9 ]
line(   29, 1:2) = [    11,    12 ]
line(   30, 1:2) = [    12,    10 ]


          20
geom.face(    1).n_poi =     3; allocate(geom.face(    1).poi(    3)); geom.face(    1).poi(1:    3) = [    1,     2,     3 ]
geom.face(    2).n_poi =     3; allocate(geom.face(    2).poi(    3)); geom.face(    2).poi(1:    3) = [    1,     3,     4 ]
geom.face(    3).n_poi =     3; allocate(geom.face(    3).poi(    3)); geom.face(    3).poi(1:    3) = [    1,     4,     5 ]
geom.face(    4).n_poi =     3; allocate(geom.face(    4).poi(    3)); geom.face(    4).poi(1:    3) = [    1,     5,     6 ]
geom.face(    5).n_poi =     3; allocate(geom.face(    5).poi(    3)); geom.face(    5).poi(1:    3) = [    1,     6,     2 ]
geom.face(    6).n_poi =     3; allocate(geom.face(    6).poi(    3)); geom.face(    6).poi(1:    3) = [    2,     6,     8 ]
geom.face(    7).n_poi =     3; allocate(geom.face(    7).poi(    3)); geom.face(    7).poi(1:    3) = [    2,     8,     7 ]
geom.face(    8).n_poi =     3; allocate(geom.face(    8).poi(    3)); geom.face(    8).poi(1:    3) = [    2,     7,     3 ]
geom.face(    9).n_poi =     3; allocate(geom.face(    9).poi(    3)); geom.face(    9).poi(1:    3) = [    3,     7,     9 ]
geom.face(   10).n_poi =     3; allocate(geom.face(   10).poi(    3)); geom.face(   10).poi(1:    3) = [    3,     9,     4 ]
geom.face(   11).n_poi =     3; allocate(geom.face(   11).poi(    3)); geom.face(   11).poi(1:    3) = [    4,     9,    10 ]
geom.face(   12).n_poi =     3; allocate(geom.face(   12).poi(    3)); geom.face(   12).poi(1:    3) = [    4,    10,     5 ]
geom.face(   13).n_poi =     3; allocate(geom.face(   13).poi(    3)); geom.face(   13).poi(1:    3) = [    5,    10,    11 ]
geom.face(   14).n_poi =     3; allocate(geom.face(   14).poi(    3)); geom.face(   14).poi(1:    3) = [    5,    11,     6 ]
geom.face(   15).n_poi =     3; allocate(geom.face(   15).poi(    3)); geom.face(   15).poi(1:    3) = [    6,    11,     8 ]
geom.face(   16).n_poi =     3; allocate(geom.face(   16).poi(    3)); geom.face(   16).poi(1:    3) = [    7,     8,    12 ]
geom.face(   17).n_poi =     3; allocate(geom.face(   17).poi(    3)); geom.face(   17).poi(1:    3) = [    7,    12,     9 ]
geom.face(   18).n_poi =     3; allocate(geom.face(   18).poi(    3)); geom.face(   18).poi(1:    3) = [    8,    11,    12 ]
geom.face(   19).n_poi =     3; allocate(geom.face(   19).poi(    3)); geom.face(   19).poi(1:    3) = [    9,    12,    10 ]
geom.face(   20).n_poi =     3; allocate(geom.face(   20).poi(    3)); geom.face(   20).poi(1:    3) = [   10,    12,    11 ]
