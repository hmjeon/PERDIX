     27
      1  -79.9882   -1.6440    0.0000
      2  -70.6419   19.5864    0.0000
      3  -67.5988  -20.9006    0.0000
      4  -54.9956    3.7753    0.0000
      5  -49.5614   32.1046    0.0000
      6  -46.2169  -32.1627    0.0000
      7  -40.6266  -12.5319    0.0000
      8  -32.9305   13.4195    0.0000
      9  -24.7469   38.7198    0.0000
     10  -23.2681  -37.7112    0.0000
     11  -19.9458  -11.3786    0.0000
     12  -10.6472   14.5326    0.0000
     13   -0.0000  -39.4368    0.0000
     14   -0.0000   40.6682    0.0000
     15   -0.0000  -12.8497    0.0000
     16   10.6472   14.5326    0.0000
     17   19.9458  -11.3786    0.0000
     18   23.2681  -37.7112    0.0000
     19   24.7469   38.7198    0.0000
     20   32.9305   13.4195    0.0000
     21   40.6265  -12.5319    0.0000
     22   46.2168  -32.1627    0.0000
     23   49.5614   32.1046    0.0000
     24   54.9956    3.7753    0.0000
     25   67.5988  -20.9006    0.0000
     26   70.6420   19.5864    0.0000
     27   79.9882   -1.6440    0.0000

     36
      3              17     21     20
      3              20     21     24
      3              27     26     24
      3              24     25     27
      3              21     25     24
      3              23     19     20
      3              20     24     23
      3              23     24     26
      3               1      3      4
      3               4      2      1
      3               2      4      5
      3              22     25     21
      3              13     15     10
      3              10     15     11
      3              11     15     12
      3              14      9     12
      3               7     10     11
      3               6     10      7
      3               3      6      7
      3               7      4      3
      3              18     15     13
      3              17     15     18
      3              21     17     18
      3              18     22     21
      3              16     15     17
      3              16     12     15
      3              16     17     20
      3              20     19     16
      3              16     19     14
      3              14     12     16
      3               8      5      4
      3               4      7      8
      3               9      5      8
      3               8     12      9
      3              11     12      8
      3               8      7     11

     62
       2      17      21
       2      21      20
       2      20      17
       2      21      24
       2      24      20
       2      27      26
       2      26      24
       2      24      27
       2      24      25
       2      25      27
       2      21      25
       2      23      19
       2      19      20
       2      20      23
       2      24      23
       2      26      23
       2       1       3
       2       3       4
       2       4       1
       2       4       2
       2       2       1
       2       4       5
       2       5       2
       2      22      25
       2      21      22
       2      13      15
       2      15      10
       2      10      13
       2      15      11
       2      11      10
       2      15      12
       2      12      11
       2      14       9
       2       9      12
       2      12      14
       2       7      10
       2      11       7
       2       6      10
       2       7       6
       2       3       6
       2       7       3
       2       7       4
       2      18      15
       2      13      18
       2      17      15
       2      18      17
       2      18      21
       2      18      22
       2      16      15
       2      17      16
       2      16      12
       2      20      16
       2      19      16
       2      19      14
       2      14      16
       2       8       5
       2       4       8
       2       7       8
       2       9       5
       2       8       9
       2       8      12
       2       8      11


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =   27
    geom.n_face =   36

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [   -79.9882d0,    -1.6440d0,     0.0000d0 ]
    geom.iniP(      2).pos(1:3) = [   -70.6419d0,    19.5864d0,     0.0000d0 ]
    geom.iniP(      3).pos(1:3) = [   -67.5988d0,   -20.9006d0,     0.0000d0 ]
    geom.iniP(      4).pos(1:3) = [   -54.9956d0,     3.7753d0,     0.0000d0 ]
    geom.iniP(      5).pos(1:3) = [   -49.5614d0,    32.1046d0,     0.0000d0 ]
    geom.iniP(      6).pos(1:3) = [   -46.2169d0,   -32.1627d0,     0.0000d0 ]
    geom.iniP(      7).pos(1:3) = [   -40.6266d0,   -12.5319d0,     0.0000d0 ]
    geom.iniP(      8).pos(1:3) = [   -32.9305d0,    13.4195d0,     0.0000d0 ]
    geom.iniP(      9).pos(1:3) = [   -24.7469d0,    38.7198d0,     0.0000d0 ]
    geom.iniP(     10).pos(1:3) = [   -23.2681d0,   -37.7112d0,     0.0000d0 ]
    geom.iniP(     11).pos(1:3) = [   -19.9458d0,   -11.3786d0,     0.0000d0 ]
    geom.iniP(     12).pos(1:3) = [   -10.6472d0,    14.5326d0,     0.0000d0 ]
    geom.iniP(     13).pos(1:3) = [    -0.0000d0,   -39.4368d0,     0.0000d0 ]
    geom.iniP(     14).pos(1:3) = [    -0.0000d0,    40.6682d0,     0.0000d0 ]
    geom.iniP(     15).pos(1:3) = [    -0.0000d0,   -12.8497d0,     0.0000d0 ]
    geom.iniP(     16).pos(1:3) = [    10.6472d0,    14.5326d0,     0.0000d0 ]
    geom.iniP(     17).pos(1:3) = [    19.9458d0,   -11.3786d0,     0.0000d0 ]
    geom.iniP(     18).pos(1:3) = [    23.2681d0,   -37.7112d0,     0.0000d0 ]
    geom.iniP(     19).pos(1:3) = [    24.7469d0,    38.7198d0,     0.0000d0 ]
    geom.iniP(     20).pos(1:3) = [    32.9305d0,    13.4195d0,     0.0000d0 ]
    geom.iniP(     21).pos(1:3) = [    40.6265d0,   -12.5319d0,     0.0000d0 ]
    geom.iniP(     22).pos(1:3) = [    46.2168d0,   -32.1627d0,     0.0000d0 ]
    geom.iniP(     23).pos(1:3) = [    49.5614d0,    32.1046d0,     0.0000d0 ]
    geom.iniP(     24).pos(1:3) = [    54.9956d0,     3.7753d0,     0.0000d0 ]
    geom.iniP(     25).pos(1:3) = [    67.5988d0,   -20.9006d0,     0.0000d0 ]
    geom.iniP(     26).pos(1:3) = [    70.6420d0,    19.5864d0,     0.0000d0 ]
    geom.iniP(     27).pos(1:3) = [    79.9882d0,    -1.6440d0,     0.0000d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [     17,      21,      20 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [     20,      21,      24 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [     27,      26,      24 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [     24,      25,      27 ]
    geom.face(      5).n_poi =       3; allocate(geom.face(      5).poi(      3)); geom.face(      5).poi(1:      3) = [     21,      25,      24 ]
    geom.face(      6).n_poi =       3; allocate(geom.face(      6).poi(      3)); geom.face(      6).poi(1:      3) = [     23,      19,      20 ]
    geom.face(      7).n_poi =       3; allocate(geom.face(      7).poi(      3)); geom.face(      7).poi(1:      3) = [     20,      24,      23 ]
    geom.face(      8).n_poi =       3; allocate(geom.face(      8).poi(      3)); geom.face(      8).poi(1:      3) = [     23,      24,      26 ]
    geom.face(      9).n_poi =       3; allocate(geom.face(      9).poi(      3)); geom.face(      9).poi(1:      3) = [      1,       3,       4 ]
    geom.face(     10).n_poi =       3; allocate(geom.face(     10).poi(      3)); geom.face(     10).poi(1:      3) = [      4,       2,       1 ]
    geom.face(     11).n_poi =       3; allocate(geom.face(     11).poi(      3)); geom.face(     11).poi(1:      3) = [      2,       4,       5 ]
    geom.face(     12).n_poi =       3; allocate(geom.face(     12).poi(      3)); geom.face(     12).poi(1:      3) = [     22,      25,      21 ]
    geom.face(     13).n_poi =       3; allocate(geom.face(     13).poi(      3)); geom.face(     13).poi(1:      3) = [     13,      15,      10 ]
    geom.face(     14).n_poi =       3; allocate(geom.face(     14).poi(      3)); geom.face(     14).poi(1:      3) = [     10,      15,      11 ]
    geom.face(     15).n_poi =       3; allocate(geom.face(     15).poi(      3)); geom.face(     15).poi(1:      3) = [     11,      15,      12 ]
    geom.face(     16).n_poi =       3; allocate(geom.face(     16).poi(      3)); geom.face(     16).poi(1:      3) = [     14,       9,      12 ]
    geom.face(     17).n_poi =       3; allocate(geom.face(     17).poi(      3)); geom.face(     17).poi(1:      3) = [      7,      10,      11 ]
    geom.face(     18).n_poi =       3; allocate(geom.face(     18).poi(      3)); geom.face(     18).poi(1:      3) = [      6,      10,       7 ]
    geom.face(     19).n_poi =       3; allocate(geom.face(     19).poi(      3)); geom.face(     19).poi(1:      3) = [      3,       6,       7 ]
    geom.face(     20).n_poi =       3; allocate(geom.face(     20).poi(      3)); geom.face(     20).poi(1:      3) = [      7,       4,       3 ]
    geom.face(     21).n_poi =       3; allocate(geom.face(     21).poi(      3)); geom.face(     21).poi(1:      3) = [     18,      15,      13 ]
    geom.face(     22).n_poi =       3; allocate(geom.face(     22).poi(      3)); geom.face(     22).poi(1:      3) = [     17,      15,      18 ]
    geom.face(     23).n_poi =       3; allocate(geom.face(     23).poi(      3)); geom.face(     23).poi(1:      3) = [     21,      17,      18 ]
    geom.face(     24).n_poi =       3; allocate(geom.face(     24).poi(      3)); geom.face(     24).poi(1:      3) = [     18,      22,      21 ]
    geom.face(     25).n_poi =       3; allocate(geom.face(     25).poi(      3)); geom.face(     25).poi(1:      3) = [     16,      15,      17 ]
    geom.face(     26).n_poi =       3; allocate(geom.face(     26).poi(      3)); geom.face(     26).poi(1:      3) = [     16,      12,      15 ]
    geom.face(     27).n_poi =       3; allocate(geom.face(     27).poi(      3)); geom.face(     27).poi(1:      3) = [     16,      17,      20 ]
    geom.face(     28).n_poi =       3; allocate(geom.face(     28).poi(      3)); geom.face(     28).poi(1:      3) = [     20,      19,      16 ]
    geom.face(     29).n_poi =       3; allocate(geom.face(     29).poi(      3)); geom.face(     29).poi(1:      3) = [     16,      19,      14 ]
    geom.face(     30).n_poi =       3; allocate(geom.face(     30).poi(      3)); geom.face(     30).poi(1:      3) = [     14,      12,      16 ]
    geom.face(     31).n_poi =       3; allocate(geom.face(     31).poi(      3)); geom.face(     31).poi(1:      3) = [      8,       5,       4 ]
    geom.face(     32).n_poi =       3; allocate(geom.face(     32).poi(      3)); geom.face(     32).poi(1:      3) = [      4,       7,       8 ]
    geom.face(     33).n_poi =       3; allocate(geom.face(     33).poi(      3)); geom.face(     33).poi(1:      3) = [      9,       5,       8 ]
    geom.face(     34).n_poi =       3; allocate(geom.face(     34).poi(      3)); geom.face(     34).poi(1:      3) = [      8,      12,       9 ]
    geom.face(     35).n_poi =       3; allocate(geom.face(     35).poi(      3)); geom.face(     35).poi(1:      3) = [     11,      12,       8 ]
    geom.face(     36).n_poi =       3; allocate(geom.face(     36).poi(      3)); geom.face(     36).poi(1:      3) = [      8,       7,      11 ]

--------------------------------------------------

     62
line(      1, 1:2) = [      17,      21 ]
line(      2, 1:2) = [      21,      20 ]
line(      3, 1:2) = [      20,      17 ]
line(      4, 1:2) = [      21,      24 ]
line(      5, 1:2) = [      24,      20 ]
line(      6, 1:2) = [      27,      26 ]
line(      7, 1:2) = [      26,      24 ]
line(      8, 1:2) = [      24,      27 ]
line(      9, 1:2) = [      24,      25 ]
line(     10, 1:2) = [      25,      27 ]
line(     11, 1:2) = [      21,      25 ]
line(     12, 1:2) = [      23,      19 ]
line(     13, 1:2) = [      19,      20 ]
line(     14, 1:2) = [      20,      23 ]
line(     15, 1:2) = [      24,      23 ]
line(     16, 1:2) = [      26,      23 ]
line(     17, 1:2) = [       1,       3 ]
line(     18, 1:2) = [       3,       4 ]
line(     19, 1:2) = [       4,       1 ]
line(     20, 1:2) = [       4,       2 ]
line(     21, 1:2) = [       2,       1 ]
line(     22, 1:2) = [       4,       5 ]
line(     23, 1:2) = [       5,       2 ]
line(     24, 1:2) = [      22,      25 ]
line(     25, 1:2) = [      21,      22 ]
line(     26, 1:2) = [      13,      15 ]
line(     27, 1:2) = [      15,      10 ]
line(     28, 1:2) = [      10,      13 ]
line(     29, 1:2) = [      15,      11 ]
line(     30, 1:2) = [      11,      10 ]
line(     31, 1:2) = [      15,      12 ]
line(     32, 1:2) = [      12,      11 ]
line(     33, 1:2) = [      14,       9 ]
line(     34, 1:2) = [       9,      12 ]
line(     35, 1:2) = [      12,      14 ]
line(     36, 1:2) = [       7,      10 ]
line(     37, 1:2) = [      11,       7 ]
line(     38, 1:2) = [       6,      10 ]
line(     39, 1:2) = [       7,       6 ]
line(     40, 1:2) = [       3,       6 ]
line(     41, 1:2) = [       7,       3 ]
line(     42, 1:2) = [       7,       4 ]
line(     43, 1:2) = [      18,      15 ]
line(     44, 1:2) = [      13,      18 ]
line(     45, 1:2) = [      17,      15 ]
line(     46, 1:2) = [      18,      17 ]
line(     47, 1:2) = [      18,      21 ]
line(     48, 1:2) = [      18,      22 ]
line(     49, 1:2) = [      16,      15 ]
line(     50, 1:2) = [      17,      16 ]
line(     51, 1:2) = [      16,      12 ]
line(     52, 1:2) = [      20,      16 ]
line(     53, 1:2) = [      19,      16 ]
line(     54, 1:2) = [      19,      14 ]
line(     55, 1:2) = [      14,      16 ]
line(     56, 1:2) = [       8,       5 ]
line(     57, 1:2) = [       4,       8 ]
line(     58, 1:2) = [       7,       8 ]
line(     59, 1:2) = [       9,       5 ]
line(     60, 1:2) = [       8,       9 ]
line(     61, 1:2) = [       8,      12 ]
line(     62, 1:2) = [       8,      11 ]

