     20
      1   17.1453   19.3422  -21.2961
      2   15.7068  -20.0261  -20.6164
      3   15.2401    4.7228   28.5272
      4   34.1515   12.9296    0.6245
      5   34.1515  -13.1214    0.6245
      6  -11.7329   -1.1709   32.9814
      7  -19.4756  -21.5229   19.2166
      8    1.5612  -34.1971   13.6500
      9   20.0907  -21.6051   18.2737
     10    0.5849  -31.8230  -14.9453
     11  -18.9211  -20.8913  -21.1337
     12   13.0759   -0.0959  -33.4768
     13  -19.2511   23.0880  -16.9362
     14   20.8505   22.2340   16.9087
     15    4.0958   32.8048   -8.2061
     16  -20.1670   15.2966   22.9382
     17   -1.8085   33.5410   12.2026
     18  -34.0509   12.9296    0.6245
     19  -17.1963    0.6866  -30.5855
     20  -34.0509  -13.1214    0.6245

     36
      3              12      1      2
      3               1      4      2
      3               5      2      4
      3              14      3      4
      3               3      9      4
      3               9      5      4
      3               8     10      9
      3               9     10      5
      3               5     10      2
      3               3      6      9
      3               6      7      9
      3               7      8      9
      3               8      7     10
      3               7     20     10
      3              20     11     10
      3              10     11      2
      3               2     11     12
      3              11     19     12
      3              12     19      1
      3              19     13      1
      3              13     15      1
      3              17     14     15
      3              14      4     15
      3               4      1     15
      3               6      3     16
      3              16      3     17
      3              14     17      3
      3               6     16      7
      3               7     16     20
      3              20     16     18
      3              15     13     17
      3              13     18     17
      3              18     16     17
      3              13     19     18
      3              18     19     20
      3              19     11     20

     54
       2      12       1
       2       1       2
       2       2      12
       2       1       4
       2       4       2
       2       5       2
       2       4       5
       2      14       3
       2       3       4
       2       4      14
       2       3       9
       2       9       4
       2       9       5
       2       8      10
       2      10       9
       2       9       8
       2      10       5
       2      10       2
       2       3       6
       2       6       9
       2       6       7
       2       7       9
       2       7       8
       2       7      10
       2       7      20
       2      20      10
       2      20      11
       2      11      10
       2      11       2
       2      11      12
       2      11      19
       2      19      12
       2      19       1
       2      19      13
       2      13       1
       2      13      15
       2      15       1
       2      17      14
       2      14      15
       2      15      17
       2       4      15
       2       3      16
       2      16       6
       2       3      17
       2      17      16
       2      16       7
       2      16      20
       2      16      18
       2      18      20
       2      13      17
       2      13      18
       2      18      17
       2      19      18
       2      19      20


--------------------------------------------------

    ! Allocate point and face structure
    geom.n_iniP =   20
    geom.n_face =   36

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [    17.1453d0,    19.3422d0,   -21.2961d0 ]
    geom.iniP(      2).pos(1:3) = [    15.7068d0,   -20.0261d0,   -20.6164d0 ]
    geom.iniP(      3).pos(1:3) = [    15.2401d0,     4.7228d0,    28.5272d0 ]
    geom.iniP(      4).pos(1:3) = [    34.1515d0,    12.9296d0,     0.6245d0 ]
    geom.iniP(      5).pos(1:3) = [    34.1515d0,   -13.1214d0,     0.6245d0 ]
    geom.iniP(      6).pos(1:3) = [   -11.7329d0,    -1.1709d0,    32.9814d0 ]
    geom.iniP(      7).pos(1:3) = [   -19.4756d0,   -21.5229d0,    19.2166d0 ]
    geom.iniP(      8).pos(1:3) = [     1.5612d0,   -34.1971d0,    13.6500d0 ]
    geom.iniP(      9).pos(1:3) = [    20.0907d0,   -21.6051d0,    18.2737d0 ]
    geom.iniP(     10).pos(1:3) = [     0.5849d0,   -31.8230d0,   -14.9453d0 ]
    geom.iniP(     11).pos(1:3) = [   -18.9211d0,   -20.8913d0,   -21.1337d0 ]
    geom.iniP(     12).pos(1:3) = [    13.0759d0,    -0.0959d0,   -33.4768d0 ]
    geom.iniP(     13).pos(1:3) = [   -19.2511d0,    23.0880d0,   -16.9362d0 ]
    geom.iniP(     14).pos(1:3) = [    20.8505d0,    22.2340d0,    16.9087d0 ]
    geom.iniP(     15).pos(1:3) = [     4.0958d0,    32.8048d0,    -8.2061d0 ]
    geom.iniP(     16).pos(1:3) = [   -20.1670d0,    15.2966d0,    22.9382d0 ]
    geom.iniP(     17).pos(1:3) = [    -1.8085d0,    33.5410d0,    12.2026d0 ]
    geom.iniP(     18).pos(1:3) = [   -34.0509d0,    12.9296d0,     0.6245d0 ]
    geom.iniP(     19).pos(1:3) = [   -17.1963d0,     0.6866d0,   -30.5855d0 ]
    geom.iniP(     20).pos(1:3) = [   -34.0509d0,   -13.1214d0,     0.6245d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       3; allocate(geom.face(      1).poi(      3)); geom.face(      1).poi(1:      3) = [     12,       1,       2 ]
    geom.face(      2).n_poi =       3; allocate(geom.face(      2).poi(      3)); geom.face(      2).poi(1:      3) = [      1,       4,       2 ]
    geom.face(      3).n_poi =       3; allocate(geom.face(      3).poi(      3)); geom.face(      3).poi(1:      3) = [      5,       2,       4 ]
    geom.face(      4).n_poi =       3; allocate(geom.face(      4).poi(      3)); geom.face(      4).poi(1:      3) = [     14,       3,       4 ]
    geom.face(      5).n_poi =       3; allocate(geom.face(      5).poi(      3)); geom.face(      5).poi(1:      3) = [      3,       9,       4 ]
    geom.face(      6).n_poi =       3; allocate(geom.face(      6).poi(      3)); geom.face(      6).poi(1:      3) = [      9,       5,       4 ]
    geom.face(      7).n_poi =       3; allocate(geom.face(      7).poi(      3)); geom.face(      7).poi(1:      3) = [      8,      10,       9 ]
    geom.face(      8).n_poi =       3; allocate(geom.face(      8).poi(      3)); geom.face(      8).poi(1:      3) = [      9,      10,       5 ]
    geom.face(      9).n_poi =       3; allocate(geom.face(      9).poi(      3)); geom.face(      9).poi(1:      3) = [      5,      10,       2 ]
    geom.face(     10).n_poi =       3; allocate(geom.face(     10).poi(      3)); geom.face(     10).poi(1:      3) = [      3,       6,       9 ]
    geom.face(     11).n_poi =       3; allocate(geom.face(     11).poi(      3)); geom.face(     11).poi(1:      3) = [      6,       7,       9 ]
    geom.face(     12).n_poi =       3; allocate(geom.face(     12).poi(      3)); geom.face(     12).poi(1:      3) = [      7,       8,       9 ]
    geom.face(     13).n_poi =       3; allocate(geom.face(     13).poi(      3)); geom.face(     13).poi(1:      3) = [      8,       7,      10 ]
    geom.face(     14).n_poi =       3; allocate(geom.face(     14).poi(      3)); geom.face(     14).poi(1:      3) = [      7,      20,      10 ]
    geom.face(     15).n_poi =       3; allocate(geom.face(     15).poi(      3)); geom.face(     15).poi(1:      3) = [     20,      11,      10 ]
    geom.face(     16).n_poi =       3; allocate(geom.face(     16).poi(      3)); geom.face(     16).poi(1:      3) = [     10,      11,       2 ]
    geom.face(     17).n_poi =       3; allocate(geom.face(     17).poi(      3)); geom.face(     17).poi(1:      3) = [      2,      11,      12 ]
    geom.face(     18).n_poi =       3; allocate(geom.face(     18).poi(      3)); geom.face(     18).poi(1:      3) = [     11,      19,      12 ]
    geom.face(     19).n_poi =       3; allocate(geom.face(     19).poi(      3)); geom.face(     19).poi(1:      3) = [     12,      19,       1 ]
    geom.face(     20).n_poi =       3; allocate(geom.face(     20).poi(      3)); geom.face(     20).poi(1:      3) = [     19,      13,       1 ]
    geom.face(     21).n_poi =       3; allocate(geom.face(     21).poi(      3)); geom.face(     21).poi(1:      3) = [     13,      15,       1 ]
    geom.face(     22).n_poi =       3; allocate(geom.face(     22).poi(      3)); geom.face(     22).poi(1:      3) = [     17,      14,      15 ]
    geom.face(     23).n_poi =       3; allocate(geom.face(     23).poi(      3)); geom.face(     23).poi(1:      3) = [     14,       4,      15 ]
    geom.face(     24).n_poi =       3; allocate(geom.face(     24).poi(      3)); geom.face(     24).poi(1:      3) = [      4,       1,      15 ]
    geom.face(     25).n_poi =       3; allocate(geom.face(     25).poi(      3)); geom.face(     25).poi(1:      3) = [      6,       3,      16 ]
    geom.face(     26).n_poi =       3; allocate(geom.face(     26).poi(      3)); geom.face(     26).poi(1:      3) = [     16,       3,      17 ]
    geom.face(     27).n_poi =       3; allocate(geom.face(     27).poi(      3)); geom.face(     27).poi(1:      3) = [     14,      17,       3 ]
    geom.face(     28).n_poi =       3; allocate(geom.face(     28).poi(      3)); geom.face(     28).poi(1:      3) = [      6,      16,       7 ]
    geom.face(     29).n_poi =       3; allocate(geom.face(     29).poi(      3)); geom.face(     29).poi(1:      3) = [      7,      16,      20 ]
    geom.face(     30).n_poi =       3; allocate(geom.face(     30).poi(      3)); geom.face(     30).poi(1:      3) = [     20,      16,      18 ]
    geom.face(     31).n_poi =       3; allocate(geom.face(     31).poi(      3)); geom.face(     31).poi(1:      3) = [     15,      13,      17 ]
    geom.face(     32).n_poi =       3; allocate(geom.face(     32).poi(      3)); geom.face(     32).poi(1:      3) = [     13,      18,      17 ]
    geom.face(     33).n_poi =       3; allocate(geom.face(     33).poi(      3)); geom.face(     33).poi(1:      3) = [     18,      16,      17 ]
    geom.face(     34).n_poi =       3; allocate(geom.face(     34).poi(      3)); geom.face(     34).poi(1:      3) = [     13,      19,      18 ]
    geom.face(     35).n_poi =       3; allocate(geom.face(     35).poi(      3)); geom.face(     35).poi(1:      3) = [     18,      19,      20 ]
    geom.face(     36).n_poi =       3; allocate(geom.face(     36).poi(      3)); geom.face(     36).poi(1:      3) = [     19,      11,      20 ]

--------------------------------------------------

     54
line(      1, 1:2) = [      12,       1 ]
line(      2, 1:2) = [       1,       2 ]
line(      3, 1:2) = [       2,      12 ]
line(      4, 1:2) = [       1,       4 ]
line(      5, 1:2) = [       4,       2 ]
line(      6, 1:2) = [       5,       2 ]
line(      7, 1:2) = [       4,       5 ]
line(      8, 1:2) = [      14,       3 ]
line(      9, 1:2) = [       3,       4 ]
line(     10, 1:2) = [       4,      14 ]
line(     11, 1:2) = [       3,       9 ]
line(     12, 1:2) = [       9,       4 ]
line(     13, 1:2) = [       9,       5 ]
line(     14, 1:2) = [       8,      10 ]
line(     15, 1:2) = [      10,       9 ]
line(     16, 1:2) = [       9,       8 ]
line(     17, 1:2) = [      10,       5 ]
line(     18, 1:2) = [      10,       2 ]
line(     19, 1:2) = [       3,       6 ]
line(     20, 1:2) = [       6,       9 ]
line(     21, 1:2) = [       6,       7 ]
line(     22, 1:2) = [       7,       9 ]
line(     23, 1:2) = [       7,       8 ]
line(     24, 1:2) = [       7,      10 ]
line(     25, 1:2) = [       7,      20 ]
line(     26, 1:2) = [      20,      10 ]
line(     27, 1:2) = [      20,      11 ]
line(     28, 1:2) = [      11,      10 ]
line(     29, 1:2) = [      11,       2 ]
line(     30, 1:2) = [      11,      12 ]
line(     31, 1:2) = [      11,      19 ]
line(     32, 1:2) = [      19,      12 ]
line(     33, 1:2) = [      19,       1 ]
line(     34, 1:2) = [      19,      13 ]
line(     35, 1:2) = [      13,       1 ]
line(     36, 1:2) = [      13,      15 ]
line(     37, 1:2) = [      15,       1 ]
line(     38, 1:2) = [      17,      14 ]
line(     39, 1:2) = [      14,      15 ]
line(     40, 1:2) = [      15,      17 ]
line(     41, 1:2) = [       4,      15 ]
line(     42, 1:2) = [       3,      16 ]
line(     43, 1:2) = [      16,       6 ]
line(     44, 1:2) = [       3,      17 ]
line(     45, 1:2) = [      17,      16 ]
line(     46, 1:2) = [      16,       7 ]
line(     47, 1:2) = [      16,      20 ]
line(     48, 1:2) = [      16,      18 ]
line(     49, 1:2) = [      18,      20 ]
line(     50, 1:2) = [      13,      17 ]
line(     51, 1:2) = [      13,      18 ]
line(     52, 1:2) = [      18,      17 ]
line(     53, 1:2) = [      19,      18 ]
line(     54, 1:2) = [      19,      20 ]

