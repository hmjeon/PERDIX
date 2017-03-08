!
! ---------------------------------------------------------------------------------------
!
!                             Module for Exam_Miscellaneous
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2016/04/06
!                                                            Last  modified  : 2016/07/14
!
! ---------------------------------------------------------------------------------------
!
module Exam_Miscellaneous

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Miscellaneous_Heptagonal_Bipyramid        ! 36. V=9,   E=21,  F=14
    public Exam_Miscellaneous_Enneagonal_Trapezohedron    ! 37. V=20,  E=36,  F=18
    public Exam_Miscellaneous_Small_Stell_Dodecahedron    ! 38. V=32,  E=90,  F=60
    public Exam_Miscellaneous_Rhombic_Hexecontahedron     ! 39. V=62,  E=120, F=60
    public Exam_Miscellaneous_Goldberg_Dk5dgD             ! 40. V=140, E=210, F=72
    public Exam_Miscellaneous_Double_Helix                ! 41. V=36,  E=100, F=66
    public Exam_Miscellaneous_Nested_Cube                 ! 42. V=16,  E=32,  F=16
    public Exam_Miscellaneous_Nested_Octahedron           ! 43. V=12,  E=30,  F=18
    public Exam_Miscellaneous_Torus                       ! 44. V=36,  E=108, F=72
    public Exam_Miscellaneous_Double_Torus                ! 45. V=44,  E=92,  F=46

contains

! ---------------------------------------------------------------------------------------

! Example of Heptagonal Bipyramid
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Heptagonal_Bipyramid(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "36_Hepta_Bipyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Heptagonal bipyramid"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.9d0      ! Cylindrical model
    prob.move_x   =-0.5d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP =  9
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  8.2664d0,  10.0011d0,  19.0513d0 ]
    geom.iniP( 2).pos(1:3) = [-11.3790d0,   6.2349d0,  19.0513d0 ]
    geom.iniP( 3).pos(1:3) = [  8.2664d0, -43.1226d0,  19.0513d0 ]
    geom.iniP( 4).pos(1:3) = [ 21.6853d0,   6.2349d0,   4.7043d0 ]
    geom.iniP( 5).pos(1:3) = [ -8.2668d0,  43.1235d0, -19.0513d0 ]
    geom.iniP( 6).pos(1:3) = [-22.4538d0,  -2.2257d0,   4.7043d0 ]
    geom.iniP( 7).pos(1:3) = [ 18.7751d0,  -2.2257d0, -13.1829d0 ]
    geom.iniP( 8).pos(1:3) = [-16.6214d0,  -9.0102d0, -13.1829d0 ]
    geom.iniP( 9).pos(1:3) = [  1.7279d0,  -9.0102d0, -21.1455d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 3, 1, 2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 3, 4, 1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1, 5, 2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2, 6, 3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 4, 5 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4, 3, 7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 2, 5, 6 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 3, 6, 8 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 4, 7, 5 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 3, 9, 7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 5, 8, 6 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 8, 9, 3 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 5, 7, 9 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 8, 5, 9 ]
end subroutine Exam_Miscellaneous_Heptagonal_Bipyramid

! ---------------------------------------------------------------------------------------

! Example of Enneagonal Trapezohedron
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Enneagonal_Trapezohedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "37_Ennea_Trapezohedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Enneagonal trapezohedron"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.86d0     ! Cylindrical model
    prob.move_x   =-3.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 20
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ 16.0321d0,   10.0013d0,  46.1619d0 ]
    geom.iniP( 2).pos(1:3) = [ -1.9342d0,   18.7964d0,  46.1619d0 ]
    geom.iniP( 3).pos(1:3) = [-17.7323d0,    6.5288d0,  46.1619d0 ]
    geom.iniP( 4).pos(1:3) = [ 16.0321d0, -155.8380d0,  46.1619d0 ]
    geom.iniP( 5).pos(1:3) = [ 42.5275d0,    6.5288d0,  25.2332d0 ]
    geom.iniP( 6).pos(1:3) = [ 30.1299d0,   18.7964d0,  35.0265d0 ]
    geom.iniP( 7).pos(1:3) = [-16.0321d0,  155.8380d0, -46.1619d0 ]
    geom.iniP( 8).pos(1:3) = [-33.3263d0,   12.2676d0,  35.0265d0 ]
    geom.iniP( 9).pos(1:3) = [-42.9655d0,   -2.2663d0,  25.2332d0 ]
    geom.iniP(10).pos(1:3) = [ 49.3563d0,   -2.2663d0,  -6.8289d0 ]
    geom.iniP(11).pos(1:3) = [ 47.8621d0,   12.2676d0,   6.8289d0 ]
    geom.iniP(12).pos(1:3) = [-49.3563d0,    2.2663d0,   6.8289d0 ]
    geom.iniP(13).pos(1:3) = [-47.8621d0,  -12.2676d0,  -6.8289d0 ]
    geom.iniP(14).pos(1:3) = [ 33.3263d0,  -12.2676d0, -35.0265d0 ]
    geom.iniP(15).pos(1:3) = [ 42.9655d0,    2.2663d0, -25.2332d0 ]
    geom.iniP(16).pos(1:3) = [-42.5275d0,   -6.5288d0, -25.2332d0 ]
    geom.iniP(17).pos(1:3) = [-30.1299d0,  -18.7964d0, -35.0265d0 ]
    geom.iniP(18).pos(1:3) = [  1.9342d0,  -18.7964d0, -46.1619d0 ]
    geom.iniP(19).pos(1:3) = [ 17.7323d0,   -6.5288d0, -46.1619d0 ]
    geom.iniP(20).pos(1:3) = [-16.0321d0,  -10.0013d0, -46.1619d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  1,  2,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  1,  4,  5,  6 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  2,  1,  6,  7 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  3,  2,  7,  8 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  4,  3,  8,  9 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  5,  4, 10, 11 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  6,  5, 11,  7 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  8,  7, 12,  9 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  4,  9, 12, 13 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 10,  4, 14, 15 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 11, 10, 15,  7 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 12,  7, 16, 13 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  4, 13, 16, 17 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 14,  4, 18, 19 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 15, 14, 19,  7 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 16,  7, 20, 17 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  4, 17, 20, 18 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 19, 18, 20,  7 ]
end subroutine Exam_Miscellaneous_Enneagonal_Trapezohedron

! ---------------------------------------------------------------------------------------

! Example of Small Stell Dodecahedron
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Small_Stell_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "38_Small_Stell_Dodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Small stell dodecahedron"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.86d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 32
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.0000d0,   0.0000d0,  28.0264d0 ]
    geom.iniP( 2).pos(1:3) = [ 18.6841d0,   0.0000d0,  20.8896d0 ]
    geom.iniP( 3).pos(1:3) = [ -9.3420d0,  16.1801d0,  20.8896d0 ]
    geom.iniP( 4).pos(1:3) = [ -9.3420d0, -16.1800d0,  20.8896d0 ]
    geom.iniP( 5).pos(1:3) = [ 20.8896d0,  16.1808d0,   9.3421d0 ]
    geom.iniP( 6).pos(1:3) = [ 20.8896d0, -16.1800d0,   9.3421d0 ]
    geom.iniP( 7).pos(1:3) = [-24.4577d0,  10.0005d0,   9.3421d0 ]
    geom.iniP( 8).pos(1:3) = [  3.5683d0,  26.1813d0,   9.3421d0 ]
    geom.iniP( 9).pos(1:3) = [  3.5683d0, -26.1813d0,   9.3421d0 ]
    geom.iniP(10).pos(1:3) = [-24.4577d0, -10.0004d0,   9.3421d0 ]
    geom.iniP(11).pos(1:3) = [ 24.4578d0,  10.0005d0,  -9.3419d0 ]
    geom.iniP(12).pos(1:3) = [ 24.4578d0, -10.0004d0,  -9.3419d0 ]
    geom.iniP(13).pos(1:3) = [-20.8895d0,  16.1801d0,  -9.3419d0 ]
    geom.iniP(14).pos(1:3) = [ -3.5682d0,  26.1813d0,  -9.3419d0 ]
    geom.iniP(15).pos(1:3) = [ -3.5682d0, -26.1813d0,  -9.3419d0 ]
    geom.iniP(16).pos(1:3) = [-20.8895d0, -16.1808d0,  -9.3419d0 ]
    geom.iniP(17).pos(1:3) = [  9.3421d0,  16.1801d0, -20.8894d0 ]
    geom.iniP(18).pos(1:3) = [  9.3421d0, -16.1800d0, -20.8894d0 ]
    geom.iniP(19).pos(1:3) = [-18.6840d0,   0.0000d0, -20.8894d0 ]
    geom.iniP(20).pos(1:3) = [  0.0000d0,   0.0000d0, -28.0262d0 ]
    geom.iniP(21).pos(1:3) = [ 15.1158d0,  26.1811d0,  39.5736d0 ]
    geom.iniP(22).pos(1:3) = [-30.2310d0,   0.0001d0,  39.5729d0 ]
    geom.iniP(23).pos(1:3) = [ 15.1146d0, -26.1807d0,  39.5730d0 ]
    geom.iniP(24).pos(1:3) = [ 48.9148d0,   0.0004d0,   9.3420d0 ]
    geom.iniP(25).pos(1:3) = [-24.4591d0,  42.3617d0,   9.3425d0 ]
    geom.iniP(26).pos(1:3) = [-24.4565d0, -42.3623d0,   9.3423d0 ]
    geom.iniP(27).pos(1:3) = [ 24.4567d0,  42.3624d0,  -9.3422d0 ]
    geom.iniP(28).pos(1:3) = [ 24.4590d0, -42.3620d0,  -9.3423d0 ]
    geom.iniP(29).pos(1:3) = [-48.9150d0,  -0.0004d0,  -9.3418d0 ]
    geom.iniP(30).pos(1:3) = [ 30.2318d0,   0.0000d0, -39.5737d0 ]
    geom.iniP(31).pos(1:3) = [-15.1162d0,  26.1806d0, -39.5746d0 ]
    geom.iniP(32).pos(1:3) = [-15.1159d0, -26.1813d0, -39.5737d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2, 21 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  2,  5, 21 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  5,  8, 21 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8,  3, 21 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  3,  1, 21 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  1,  3, 22 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  3,  7, 22 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  7, 10, 22 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 10,  4, 22 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4,  1, 22 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  1,  4, 23 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  4,  9, 23 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  9,  6, 23 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  6,  2, 23 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  2,  1, 23 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  2,  6, 24 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  6, 12, 24 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 12, 11, 24 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 11,  5, 24 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  5,  2, 24 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [  3,  8, 25 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [  8, 14, 25 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 14, 13, 25 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 13,  7, 25 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [  7,  3, 25 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [  4, 10, 26 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 10, 16, 26 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 16, 15, 26 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 15,  9, 26 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [  9,  4, 26 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  5, 11, 27 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 11, 17, 27 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 17, 14, 27 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 14,  8, 27 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [  8,  5, 27 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [  6,  9, 28 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [  9, 15, 28 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 15, 18, 28 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 18, 12, 28 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 12,  6, 28 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [  7, 13, 29 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 13, 19, 29 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 19, 16, 29 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 16, 10, 29 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 10,  7, 29 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 11, 12, 30 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 12, 18, 30 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 18, 20, 30 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 20, 17, 30 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 17, 11, 30 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 13, 14, 31 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 14, 17, 31 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 17, 20, 31 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 20, 19, 31 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 19, 13, 31 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 15, 16, 32 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 16, 19, 32 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 19, 20, 32 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 20, 18, 32 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 18, 15, 32 ]
end subroutine Exam_Miscellaneous_Small_Stell_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Rhombic Hexecontahedron
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Rhombic_Hexecontahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "39_Rhombic_Hexeconta"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Rhombic hexecontahedron"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.86d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 62
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -0.0000d0,   0.0000d0,  47.6814d0 ]
    geom.iniP( 2).pos(1:3) = [ 31.7873d0,   0.0000d0,  35.5395d0 ]
    geom.iniP( 3).pos(1:3) = [-15.8936d0,  27.5273d0,  35.5395d0 ]
    geom.iniP( 4).pos(1:3) = [-15.8936d0, -27.5273d0,  35.5395d0 ]
    geom.iniP( 5).pos(1:3) = [ 35.5395d0,  27.5286d0,  15.8936d0 ]
    geom.iniP( 6).pos(1:3) = [ 35.5395d0, -27.5273d0,  15.8936d0 ]
    geom.iniP( 7).pos(1:3) = [-41.6102d0,  17.0139d0,  15.8936d0 ]
    geom.iniP( 8).pos(1:3) = [  6.0707d0,  44.5425d0,  15.8936d0 ]
    geom.iniP( 9).pos(1:3) = [  6.0707d0, -44.5425d0,  15.8936d0 ]
    geom.iniP(10).pos(1:3) = [-41.6102d0, -17.0139d0,  15.8936d0 ]
    geom.iniP(11).pos(1:3) = [ 41.6102d0,  17.0139d0, -15.8936d0 ]
    geom.iniP(12).pos(1:3) = [ 41.6102d0, -17.0139d0, -15.8936d0 ]
    geom.iniP(13).pos(1:3) = [-35.5395d0,  27.5273d0, -15.8936d0 ]
    geom.iniP(14).pos(1:3) = [ -6.0707d0,  44.5425d0, -15.8936d0 ]
    geom.iniP(15).pos(1:3) = [ -6.0707d0, -44.5425d0, -15.8936d0 ]
    geom.iniP(16).pos(1:3) = [-35.5395d0, -27.5286d0, -15.8936d0 ]
    geom.iniP(17).pos(1:3) = [ 15.8936d0,  27.5273d0, -35.5395d0 ]
    geom.iniP(18).pos(1:3) = [ 15.8936d0, -27.5273d0, -35.5395d0 ]
    geom.iniP(19).pos(1:3) = [-31.7873d0,   0.0000d0, -35.5395d0 ]
    geom.iniP(20).pos(1:3) = [ -0.0000d0,   0.0000d0, -47.6814d0 ]
    geom.iniP(21).pos(1:3) = [ 12.1417d0,   0.0000d0,  31.7876d0 ]
    geom.iniP(22).pos(1:3) = [ 25.7166d0,  10.5150d0,  19.6457d0 ]
    geom.iniP(23).pos(1:3) = [ 15.8937d0,  27.5287d0,  12.1417d0 ]
    geom.iniP(24).pos(1:3) = [ -6.0708d0,  10.5145d0,  31.7876d0 ]
    geom.iniP(25).pos(1:3) = [-21.9645d0,  17.0132d0,  19.6457d0 ]
    geom.iniP(26).pos(1:3) = [-31.7874d0,   0.0000d0,  12.1417d0 ]
    geom.iniP(27).pos(1:3) = [ -6.0708d0, -10.5145d0,  31.7876d0 ]
    geom.iniP(28).pos(1:3) = [ -3.7520d0, -27.5282d0,  19.6457d0 ]
    geom.iniP(29).pos(1:3) = [ 25.7166d0, -10.5145d0,  19.6457d0 ]
    geom.iniP(30).pos(1:3) = [ 29.4686d0, -17.0132d0,   0.0000d0 ]
    geom.iniP(31).pos(1:3) = [ -3.7520d0,  27.5282d0,  19.6457d0 ]
    geom.iniP(32).pos(1:3) = [ -0.0000d0,  34.0274d0,   0.0000d0 ]
    geom.iniP(33).pos(1:3) = [-21.9645d0, -17.0132d0,  19.6457d0 ]
    geom.iniP(34).pos(1:3) = [-29.4686d0, -17.0137d0,   0.0000d0 ]
    geom.iniP(35).pos(1:3) = [ 29.4686d0,  17.0137d0,   0.0000d0 ]
    geom.iniP(36).pos(1:3) = [ 21.9645d0,  17.0132d0, -19.6457d0 ]
    geom.iniP(37).pos(1:3) = [ 15.8937d0, -27.5282d0,  12.1417d0 ]
    geom.iniP(38).pos(1:3) = [ -0.0000d0, -34.0274d0,   0.0000d0 ]
    geom.iniP(39).pos(1:3) = [  3.7520d0, -27.5282d0, -19.6457d0 ]
    geom.iniP(40).pos(1:3) = [-29.4686d0,  17.0132d0,   0.0000d0 ]
    geom.iniP(41).pos(1:3) = [-25.7166d0,  10.5145d0, -19.6457d0 ]
    geom.iniP(42).pos(1:3) = [ 31.7874d0,   0.0000d0, -12.1417d0 ]
    geom.iniP(43).pos(1:3) = [ 21.9645d0, -17.0132d0, -19.6457d0 ]
    geom.iniP(44).pos(1:3) = [  6.0708d0, -10.5145d0, -31.7876d0 ]
    geom.iniP(45).pos(1:3) = [-15.8937d0,  27.5282d0, -12.1417d0 ]
    geom.iniP(46).pos(1:3) = [  3.7520d0,  27.5282d0, -19.6457d0 ]
    geom.iniP(47).pos(1:3) = [  6.0708d0,  10.5145d0, -31.7876d0 ]
    geom.iniP(48).pos(1:3) = [-15.8937d0, -27.5287d0, -12.1417d0 ]
    geom.iniP(49).pos(1:3) = [-25.7166d0, -10.5150d0, -19.6457d0 ]
    geom.iniP(50).pos(1:3) = [-12.1417d0,   0.0000d0, -31.7876d0 ]
    geom.iniP(51).pos(1:3) = [  6.0708d0,  10.5149d0,  15.8938d0 ]
    geom.iniP(52).pos(1:3) = [-12.1417d0,   0.0000d0,  15.8938d0 ]
    geom.iniP(53).pos(1:3) = [  6.0708d0, -10.5147d0,  15.8938d0 ]
    geom.iniP(54).pos(1:3) = [ 19.6457d0,   0.0001d0,   3.7520d0 ]
    geom.iniP(55).pos(1:3) = [ -9.8229d0,  17.0134d0,   3.7520d0 ]
    geom.iniP(56).pos(1:3) = [ -9.8229d0, -17.0135d0,   3.7520d0 ]
    geom.iniP(57).pos(1:3) = [  9.8229d0,  17.0135d0,  -3.7520d0 ]
    geom.iniP(58).pos(1:3) = [  9.8229d0, -17.0134d0,  -3.7520d0 ]
    geom.iniP(59).pos(1:3) = [-19.6457d0,  -0.0001d0,  -3.7520d0 ]
    geom.iniP(60).pos(1:3) = [ 12.1417d0,   0.0000d0, -15.8938d0 ]
    geom.iniP(61).pos(1:3) = [ -6.0708d0,  10.5147d0, -15.8938d0 ]
    geom.iniP(62).pos(1:3) = [ -6.0708d0, -10.5149d0, -15.8938d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 51, 21,  2, 22 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 51, 22,  5, 23 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 51, 23,  8, 31 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 51, 31,  3, 24 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 51, 24,  1, 21 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 52, 24,  3, 25 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 52, 25,  7, 26 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 52, 26, 10, 33 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 52, 33,  4, 27 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 52, 27,  1, 24 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 53, 27,  4, 28 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 53, 28,  9, 37 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 53, 37,  6, 29 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 53, 29,  2, 21 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 53, 21,  1, 27 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 54, 29,  6, 30 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 54, 30, 12, 42 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 54, 42, 11, 35 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 54, 35,  5, 22 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 54, 22,  2, 29 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 55, 31,  8, 32 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 55, 32, 14, 45 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 55, 45, 13, 40 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 55, 40,  7, 25 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 55, 25,  3, 31 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 56, 33, 10, 34 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 56, 34, 16, 48 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 56, 48, 15, 38 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 56, 38,  9, 28 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 56, 28,  4, 33 ]
    geom.face(31).n_poi = 4; allocate(geom.face(31).poi(4)); geom.face(31).poi(1:4) = [ 57, 35, 11, 36 ]
    geom.face(32).n_poi = 4; allocate(geom.face(32).poi(4)); geom.face(32).poi(1:4) = [ 57, 36, 17, 46 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [ 57, 46, 14, 32 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [ 57, 32,  8, 23 ]
    geom.face(35).n_poi = 4; allocate(geom.face(35).poi(4)); geom.face(35).poi(1:4) = [ 57, 23,  5, 35 ]
    geom.face(36).n_poi = 4; allocate(geom.face(36).poi(4)); geom.face(36).poi(1:4) = [ 58, 37,  9, 38 ]
    geom.face(37).n_poi = 4; allocate(geom.face(37).poi(4)); geom.face(37).poi(1:4) = [ 58, 38, 15, 39 ]
    geom.face(38).n_poi = 4; allocate(geom.face(38).poi(4)); geom.face(38).poi(1:4) = [ 58, 39, 18, 43 ]
    geom.face(39).n_poi = 4; allocate(geom.face(39).poi(4)); geom.face(39).poi(1:4) = [ 58, 43, 12, 30 ]
    geom.face(40).n_poi = 4; allocate(geom.face(40).poi(4)); geom.face(40).poi(1:4) = [ 58, 30,  6, 37 ]
    geom.face(41).n_poi = 4; allocate(geom.face(41).poi(4)); geom.face(41).poi(1:4) = [ 59, 40, 13, 41 ]
    geom.face(42).n_poi = 4; allocate(geom.face(42).poi(4)); geom.face(42).poi(1:4) = [ 59, 41, 19, 49 ]
    geom.face(43).n_poi = 4; allocate(geom.face(43).poi(4)); geom.face(43).poi(1:4) = [ 59, 49, 16, 34 ]
    geom.face(44).n_poi = 4; allocate(geom.face(44).poi(4)); geom.face(44).poi(1:4) = [ 59, 34, 10, 26 ]
    geom.face(45).n_poi = 4; allocate(geom.face(45).poi(4)); geom.face(45).poi(1:4) = [ 59, 26,  7, 40 ]
    geom.face(46).n_poi = 4; allocate(geom.face(46).poi(4)); geom.face(46).poi(1:4) = [ 60, 42, 12, 43 ]
    geom.face(47).n_poi = 4; allocate(geom.face(47).poi(4)); geom.face(47).poi(1:4) = [ 60, 43, 18, 44 ]
    geom.face(48).n_poi = 4; allocate(geom.face(48).poi(4)); geom.face(48).poi(1:4) = [ 60, 44, 20, 47 ]
    geom.face(49).n_poi = 4; allocate(geom.face(49).poi(4)); geom.face(49).poi(1:4) = [ 60, 47, 17, 36 ]
    geom.face(50).n_poi = 4; allocate(geom.face(50).poi(4)); geom.face(50).poi(1:4) = [ 60, 36, 11, 42 ]
    geom.face(51).n_poi = 4; allocate(geom.face(51).poi(4)); geom.face(51).poi(1:4) = [ 61, 45, 14, 46 ]
    geom.face(52).n_poi = 4; allocate(geom.face(52).poi(4)); geom.face(52).poi(1:4) = [ 61, 46, 17, 47 ]
    geom.face(53).n_poi = 4; allocate(geom.face(53).poi(4)); geom.face(53).poi(1:4) = [ 61, 47, 20, 50 ]
    geom.face(54).n_poi = 4; allocate(geom.face(54).poi(4)); geom.face(54).poi(1:4) = [ 61, 50, 19, 41 ]
    geom.face(55).n_poi = 4; allocate(geom.face(55).poi(4)); geom.face(55).poi(1:4) = [ 61, 41, 13, 45 ]
    geom.face(56).n_poi = 4; allocate(geom.face(56).poi(4)); geom.face(56).poi(1:4) = [ 62, 48, 16, 49 ]
    geom.face(57).n_poi = 4; allocate(geom.face(57).poi(4)); geom.face(57).poi(1:4) = [ 62, 49, 19, 50 ]
    geom.face(58).n_poi = 4; allocate(geom.face(58).poi(4)); geom.face(58).poi(1:4) = [ 62, 50, 20, 44 ]
    geom.face(59).n_poi = 4; allocate(geom.face(59).poi(4)); geom.face(59).poi(1:4) = [ 62, 44, 18, 39 ]
    geom.face(60).n_poi = 4; allocate(geom.face(60).poi(4)); geom.face(60).poi(1:4) = [ 62, 39, 15, 48 ]
end subroutine Exam_Miscellaneous_Rhombic_Hexecontahedron

! ---------------------------------------------------------------------------------------

! Example of Goldberg dk5dgD
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Goldberg_Dk5dgD(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "40_Goldberg_dk5dgD"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Goldberg dk5dgD"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.91d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP =  140
    geom.n_face = 72

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(  1).pos(1:3) = [  0.0000d0,   0.0000d0,  80.1233d0 ]
    geom.iniP(  2).pos(1:3) = [ 53.4158d0,   0.0000d0,  59.7205d0 ]
    geom.iniP(  3).pos(1:3) = [-26.7075d0,  46.2589d0,  59.7205d0 ]
    geom.iniP(  4).pos(1:3) = [-26.7075d0, -46.2588d0,  59.7205d0 ]
    geom.iniP(  5).pos(1:3) = [ 59.7205d0,  46.2589d0,  26.7075d0 ]
    geom.iniP(  6).pos(1:3) = [ 59.7205d0, -46.2588d0,  26.7075d0 ]
    geom.iniP(  7).pos(1:3) = [-69.9215d0,  28.5895d0,  26.7075d0 ]
    geom.iniP(  8).pos(1:3) = [ 10.2010d0,  74.8484d0,  26.7075d0 ]
    geom.iniP(  9).pos(1:3) = [ 10.2018d0, -74.8483d0,  26.7083d0 ]
    geom.iniP( 10).pos(1:3) = [-69.9215d0, -28.5895d0,  26.7075d0 ]
    geom.iniP( 11).pos(1:3) = [ 69.9215d0,  28.5895d0, -26.7075d0 ]
    geom.iniP( 12).pos(1:3) = [ 69.9215d0, -28.5895d0, -26.7075d0 ]
    geom.iniP( 13).pos(1:3) = [-59.7205d0,  46.2589d0, -26.7075d0 ]
    geom.iniP( 14).pos(1:3) = [-10.2018d0,  74.8484d0, -26.7083d0 ]
    geom.iniP( 15).pos(1:3) = [-10.2010d0, -74.8483d0, -26.7075d0 ]
    geom.iniP( 16).pos(1:3) = [-59.7205d0, -46.2588d0, -26.7075d0 ]
    geom.iniP( 17).pos(1:3) = [ 26.7075d0,  46.2589d0, -59.7205d0 ]
    geom.iniP( 18).pos(1:3) = [ 26.7075d0, -46.2588d0, -59.7205d0 ]
    geom.iniP( 19).pos(1:3) = [-53.4158d0,   0.0000d0, -59.7205d0 ]
    geom.iniP( 20).pos(1:3) = [  0.0000d0,   0.0000d0, -80.1233d0 ]
    geom.iniP( 21).pos(1:3) = [  2.1082d0,  66.0404d0,  45.2832d0 ]
    geom.iniP( 22).pos(1:3) = [-11.6553d0,  34.8457d0,  71.1783d0 ]
    geom.iniP( 23).pos(1:3) = [ 20.6154d0,   7.3298d0,  77.0561d0 ]
    geom.iniP( 24).pos(1:3) = [ 54.3221d0,  21.5186d0,  54.7944d0 ]
    geom.iniP( 25).pos(1:3) = [ 42.8843d0,  57.8035d0,  35.1571d0 ]
    geom.iniP( 26).pos(1:3) = [-58.2471d0, -31.1947d0,  45.2832d0 ]
    geom.iniP( 27).pos(1:3) = [-24.3500d0, -27.5159d0,  71.1775d0 ]
    geom.iniP( 28).pos(1:3) = [-16.6554d0,  14.1881d0,  77.0561d0 ]
    geom.iniP( 29).pos(1:3) = [-45.7969d0,  36.2849d0,  54.7944d0 ]
    geom.iniP( 30).pos(1:3) = [-71.5016d0,   8.2369d0,  35.1571d0 ]
    geom.iniP( 31).pos(1:3) = [ 56.1389d0, -34.8456d0,  45.2832d0 ]
    geom.iniP( 32).pos(1:3) = [ 36.0053d0,  -7.3289d0,  71.1775d0 ]
    geom.iniP( 33).pos(1:3) = [ -3.9600d0, -21.5178d0,  77.0561d0 ]
    geom.iniP( 34).pos(1:3) = [ -8.5260d0, -57.8034d0,  54.7944d0 ]
    geom.iniP( 35).pos(1:3) = [ 28.6174d0, -66.0403d0,  35.1579d0 ]
    geom.iniP( 36).pos(1:3) = [ 73.6035d0,  31.1947d0,  -5.0798d0 ]
    geom.iniP( 37).pos(1:3) = [ 65.6015d0,  27.5167d0,  36.8193d0 ]
    geom.iniP( 38).pos(1:3) = [ 63.7848d0, -14.1880d0,  46.3313d0 ]
    geom.iniP( 39).pos(1:3) = [ 70.6646d0, -36.2848d0,  10.3101d0 ]
    geom.iniP( 40).pos(1:3) = [ 76.7328d0,  -8.2369d0, -21.4628d0 ]
    geom.iniP( 41).pos(1:3) = [-63.8175d0,  48.1449d0,  -5.0798d0 ]
    geom.iniP( 42).pos(1:3) = [-56.6310d0,  43.0539d0,  36.8193d0 ]
    geom.iniP( 43).pos(1:3) = [-19.6047d0,  62.3329d0,  46.3313d0 ]
    geom.iniP( 44).pos(1:3) = [ -3.9082d0,  79.3396d0,  10.3101d0 ]
    geom.iniP( 45).pos(1:3) = [-31.2329d0,  70.5706d0, -21.4636d0 ]
    geom.iniP( 46).pos(1:3) = [ -9.7860d0, -79.3395d0,  -5.0798d0 ]
    geom.iniP( 47).pos(1:3) = [ -8.9705d0, -70.5706d0,  36.8201d0 ]
    geom.iniP( 48).pos(1:3) = [-44.1801d0, -48.1448d0,  46.3313d0 ]
    geom.iniP( 49).pos(1:3) = [-66.7564d0, -43.0539d0,  10.3101d0 ]
    geom.iniP( 50).pos(1:3) = [-45.4998d0, -62.3336d0, -21.4628d0 ]
    geom.iniP( 51).pos(1:3) = [  3.9074d0,  79.3396d0, -10.3109d0 ]
    geom.iniP( 52).pos(1:3) = [ 31.2321d0,  70.5706d0,  21.4628d0 ]
    geom.iniP( 53).pos(1:3) = [ 63.8175d0,  48.1449d0,   5.0798d0 ]
    geom.iniP( 54).pos(1:3) = [ 56.6311d0,  43.0539d0, -36.8193d0 ]
    geom.iniP( 55).pos(1:3) = [ 19.6047d0,  62.3329d0, -46.3313d0 ]
    geom.iniP( 56).pos(1:3) = [ 44.1801d0, -48.1448d0, -46.3313d0 ]
    geom.iniP( 57).pos(1:3) = [ 66.7564d0, -43.0539d0, -10.3101d0 ]
    geom.iniP( 58).pos(1:3) = [ 45.4998d0, -62.3328d0,  21.4636d0 ]
    geom.iniP( 59).pos(1:3) = [  9.7868d0, -79.3395d0,   5.0798d0 ]
    geom.iniP( 60).pos(1:3) = [  8.9713d0, -70.5706d0, -36.8193d0 ]
    geom.iniP( 61).pos(1:3) = [-70.6646d0, -36.2856d0, -10.3101d0 ]
    geom.iniP( 62).pos(1:3) = [-76.7319d0,  -8.2377d0,  21.4628d0 ]
    geom.iniP( 63).pos(1:3) = [-73.6035d0,  31.1939d0,   5.0798d0 ]
    geom.iniP( 64).pos(1:3) = [-65.6015d0,  27.5159d0, -36.8193d0 ]
    geom.iniP( 65).pos(1:3) = [-63.7848d0, -14.1888d0, -46.3305d0 ]
    geom.iniP( 66).pos(1:3) = [ 16.6554d0,  14.1881d0, -77.0561d0 ]
    geom.iniP( 67).pos(1:3) = [ 45.7961d0,  36.2857d0, -54.7944d0 ]
    geom.iniP( 68).pos(1:3) = [ 71.5016d0,   8.2377d0, -35.1571d0 ]
    geom.iniP( 69).pos(1:3) = [ 58.2471d0, -31.1939d0, -45.2832d0 ]
    geom.iniP( 70).pos(1:3) = [ 24.3500d0, -27.5167d0, -71.1775d0 ]
    geom.iniP( 71).pos(1:3) = [-20.6154d0,   7.3290d0, -77.0561d0 ]
    geom.iniP( 72).pos(1:3) = [-54.3229d0,  21.5178d0, -54.7944d0 ]
    geom.iniP( 73).pos(1:3) = [-42.8851d0,  57.8035d0, -35.1579d0 ]
    geom.iniP( 74).pos(1:3) = [ -2.1082d0,  66.0404d0, -45.2832d0 ]
    geom.iniP( 75).pos(1:3) = [ 11.6553d0,  34.8457d0, -71.1783d0 ]
    geom.iniP( 76).pos(1:3) = [  3.9600d0, -21.5178d0, -77.0561d0 ]
    geom.iniP( 77).pos(1:3) = [  8.5260d0, -57.8034d0, -54.7944d0 ]
    geom.iniP( 78).pos(1:3) = [-28.6174d0, -66.0411d0, -35.1571d0 ]
    geom.iniP( 79).pos(1:3) = [-56.1388d0, -34.8464d0, -45.2832d0 ]
    geom.iniP( 80).pos(1:3) = [-36.0053d0,  -7.3297d0, -71.1775d0 ]
    geom.iniP( 81).pos(1:3) = [  7.6133d0,  39.7064d0,  68.9124d0 ]
    geom.iniP( 82).pos(1:3) = [ 23.1481d0,  26.4471d0,  71.7509d0 ]
    geom.iniP( 83).pos(1:3) = [ 39.3894d0,  33.2703d0,  61.0331d0 ]
    geom.iniP( 84).pos(1:3) = [ 33.8923d0,  50.7477d0,  51.5704d0 ]
    geom.iniP( 85).pos(1:3) = [ 14.2533d0,  54.7252d0,  56.4399d0 ]
    geom.iniP( 86).pos(1:3) = [-38.1931d0, -13.2593d0,  68.9124d0 ]
    geom.iniP( 87).pos(1:3) = [-34.4777d0,   6.8240d0,  71.7509d0 ]
    geom.iniP( 88).pos(1:3) = [-48.5080d0,  17.4766d0,  61.0331d0 ]
    geom.iniP( 89).pos(1:3) = [-60.8945d0,   3.9775d0,  51.5696d0 ]
    geom.iniP( 90).pos(1:3) = [-54.5196d0, -15.0187d0,  56.4399d0 ]
    geom.iniP( 91).pos(1:3) = [ 30.5799d0, -26.4462d0,  68.9124d0 ]
    geom.iniP( 92).pos(1:3) = [ 11.3296d0, -33.2702d0,  71.7509d0 ]
    geom.iniP( 93).pos(1:3) = [  9.1186d0, -50.7468d0,  61.0331d0 ]
    geom.iniP( 94).pos(1:3) = [ 27.0030d0, -54.7243d0,  51.5704d0 ]
    geom.iniP( 95).pos(1:3) = [ 40.2671d0, -39.7056d0,  56.4399d0 ]
    geom.iniP( 96).pos(1:3) = [ 74.4095d0,  13.2602d0,  25.9015d0 ]
    geom.iniP( 97).pos(1:3) = [ 73.5318d0,  -6.8232d0,  30.4946d0 ]
    geom.iniP( 98).pos(1:3) = [ 76.8443d0, -17.4766d0,  13.1526d0 ]
    geom.iniP( 99).pos(1:3) = [ 79.7680d0,  -3.9775d0,  -2.1584d0 ]
    geom.iniP(100).pos(1:3) = [ 78.2635d0,  15.0188d0,   5.7209d0 ]
    geom.iniP(101).pos(1:3) = [-48.6880d0,  57.8098d0,  25.9015d0 ]
    geom.iniP(102).pos(1:3) = [-30.8570d0,  67.0925d0,  30.4946d0 ]
    geom.iniP(103).pos(1:3) = [-23.2867d0,  75.2872d0,  13.1526d0 ]
    geom.iniP(104).pos(1:3) = [-36.4401d0,  71.0700d0,  -2.1584d0 ]
    geom.iniP(105).pos(1:3) = [-52.1390d0,  60.2685d0,   5.7209d0 ]
    geom.iniP(106).pos(1:3) = [-25.7215d0, -71.0699d0,  25.9023d0 ]
    geom.iniP(107).pos(1:3) = [-42.6756d0, -60.2684d0,  30.4946d0 ]
    geom.iniP(108).pos(1:3) = [-53.5575d0, -57.8106d0,  13.1526d0 ]
    geom.iniP(109).pos(1:3) = [-43.3287d0, -67.0925d0,  -2.1584d0 ]
    geom.iniP(110).pos(1:3) = [-26.1253d0, -75.2872d0,   5.7209d0 ]
    geom.iniP(111).pos(1:3) = [ 36.4394d0,  71.0700d0,   2.1584d0 ]
    geom.iniP(112).pos(1:3) = [ 52.1383d0,  60.2685d0,  -5.7209d0 ]
    geom.iniP(113).pos(1:3) = [ 48.6880d0,  57.8106d0, -25.9023d0 ]
    geom.iniP(114).pos(1:3) = [ 30.8562d0,  67.0925d0, -30.4946d0 ]
    geom.iniP(115).pos(1:3) = [ 23.2867d0,  75.2872d0, -13.1526d0 ]
    geom.iniP(116).pos(1:3) = [ 53.5583d0, -57.8098d0, -13.1526d0 ]
    geom.iniP(117).pos(1:3) = [ 43.3295d0, -67.0925d0,   2.1584d0 ]
    geom.iniP(118).pos(1:3) = [ 26.1253d0, -75.2872d0,  -5.7209d0 ]
    geom.iniP(119).pos(1:3) = [ 25.7215d0, -71.0699d0, -25.9023d0 ]
    geom.iniP(120).pos(1:3) = [ 42.6756d0, -60.2684d0, -30.4946d0 ]
    geom.iniP(121).pos(1:3) = [-79.7680d0,  -3.9783d0,   2.1584d0 ]
    geom.iniP(122).pos(1:3) = [-78.2635d0,  15.0180d0,  -5.7209d0 ]
    geom.iniP(123).pos(1:3) = [-74.4095d0,  13.2594d0, -25.9015d0 ]
    geom.iniP(124).pos(1:3) = [-73.5318d0,  -6.8240d0, -30.4946d0 ]
    geom.iniP(125).pos(1:3) = [-76.8442d0, -17.4774d0, -13.1526d0 ]
    geom.iniP(126).pos(1:3) = [ 48.5080d0,  17.4766d0, -61.0330d0 ]
    geom.iniP(127).pos(1:3) = [ 60.8945d0,   3.9775d0, -51.5704d0 ]
    geom.iniP(128).pos(1:3) = [ 54.5197d0, -15.0187d0, -56.4399d0 ]
    geom.iniP(129).pos(1:3) = [ 38.1932d0, -13.2593d0, -68.9124d0 ]
    geom.iniP(130).pos(1:3) = [ 34.4777d0,   6.8240d0, -71.7509d0 ]
    geom.iniP(131).pos(1:3) = [-39.3894d0,  33.2703d0, -61.0330d0 ]
    geom.iniP(132).pos(1:3) = [-33.8923d0,  50.7469d0, -51.5704d0 ]
    geom.iniP(133).pos(1:3) = [-14.2533d0,  54.7244d0, -56.4399d0 ]
    geom.iniP(134).pos(1:3) = [ -7.6133d0,  39.7056d0, -68.9124d0 ]
    geom.iniP(135).pos(1:3) = [-23.1481d0,  26.4463d0, -71.7509d0 ]
    geom.iniP(136).pos(1:3) = [ -9.1186d0, -50.7476d0, -61.0330d0 ]
    geom.iniP(137).pos(1:3) = [-27.0022d0, -54.7251d0, -51.5704d0 ]
    geom.iniP(138).pos(1:3) = [-40.2663d0, -39.7064d0, -56.4399d0 ]
    geom.iniP(139).pos(1:3) = [-30.5798d0, -26.4470d0, -68.9124d0 ]
    geom.iniP(140).pos(1:3) = [-11.3296d0, -33.2702d0, -71.7509d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 6; allocate(geom.face( 1).poi(6)); geom.face( 1).poi(1:6) = [  22,  81,  85,  21,  43,  3 ]
    geom.face( 2).n_poi = 6; allocate(geom.face( 2).poi(6)); geom.face( 2).poi(1:6) = [  23,  82,  81,  22,  28,  1 ]
    geom.face( 3).n_poi = 6; allocate(geom.face( 3).poi(6)); geom.face( 3).poi(1:6) = [  24,  83,  82,  23,  32,  2 ]
    geom.face( 4).n_poi = 6; allocate(geom.face( 4).poi(6)); geom.face( 4).poi(1:6) = [  25,  84,  83,  24,  37,  5 ]
    geom.face( 5).n_poi = 6; allocate(geom.face( 5).poi(6)); geom.face( 5).poi(1:6) = [  21,  85,  84,  25,  52,  8 ]
    geom.face( 6).n_poi = 6; allocate(geom.face( 6).poi(6)); geom.face( 6).poi(1:6) = [  27,  86,  90,  26,  48,  4 ]
    geom.face( 7).n_poi = 6; allocate(geom.face( 7).poi(6)); geom.face( 7).poi(1:6) = [  28,  87,  86,  27,  33,  1 ]
    geom.face( 8).n_poi = 6; allocate(geom.face( 8).poi(6)); geom.face( 8).poi(1:6) = [  29,  88,  87,  28,  22,  3 ]
    geom.face( 9).n_poi = 6; allocate(geom.face( 9).poi(6)); geom.face( 9).poi(1:6) = [  30,  89,  88,  29,  42,  7 ]
    geom.face(10).n_poi = 6; allocate(geom.face(10).poi(6)); geom.face(10).poi(1:6) = [  26,  90,  89,  30,  62, 10 ]
    geom.face(11).n_poi = 6; allocate(geom.face(11).poi(6)); geom.face(11).poi(1:6) = [  32,  91,  95,  31,  38,  2 ]
    geom.face(12).n_poi = 6; allocate(geom.face(12).poi(6)); geom.face(12).poi(1:6) = [  33,  92,  91,  32,  23,  1 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [  34,  93,  92,  33,  27,  4 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [  35,  94,  93,  34,  47,  9 ]
    geom.face(15).n_poi = 6; allocate(geom.face(15).poi(6)); geom.face(15).poi(1:6) = [  31,  95,  94,  35,  58,  6 ]
    geom.face(16).n_poi = 6; allocate(geom.face(16).poi(6)); geom.face(16).poi(1:6) = [  37,  96, 100,  36,  53,  5 ]
    geom.face(17).n_poi = 6; allocate(geom.face(17).poi(6)); geom.face(17).poi(1:6) = [  38,  97,  96,  37,  24,  2 ]
    geom.face(18).n_poi = 6; allocate(geom.face(18).poi(6)); geom.face(18).poi(1:6) = [  39,  98,  97,  38,  31,  6 ]
    geom.face(19).n_poi = 6; allocate(geom.face(19).poi(6)); geom.face(19).poi(1:6) = [  40,  99,  98,  39,  57, 12 ]
    geom.face(20).n_poi = 6; allocate(geom.face(20).poi(6)); geom.face(20).poi(1:6) = [  36, 100,  99,  40,  68, 11 ]
    geom.face(21).n_poi = 6; allocate(geom.face(21).poi(6)); geom.face(21).poi(1:6) = [  42, 101, 105,  41,  63,  7 ]
    geom.face(22).n_poi = 6; allocate(geom.face(22).poi(6)); geom.face(22).poi(1:6) = [  43, 102, 101,  42,  29,  3 ]
    geom.face(23).n_poi = 6; allocate(geom.face(23).poi(6)); geom.face(23).poi(1:6) = [  44, 103, 102,  43,  21,  8 ]
    geom.face(24).n_poi = 6; allocate(geom.face(24).poi(6)); geom.face(24).poi(1:6) = [  45, 104, 103,  44,  51, 14 ]
    geom.face(25).n_poi = 6; allocate(geom.face(25).poi(6)); geom.face(25).poi(1:6) = [  41, 105, 104,  45,  73, 13 ]
    geom.face(26).n_poi = 6; allocate(geom.face(26).poi(6)); geom.face(26).poi(1:6) = [  47, 106, 110,  46,  59,  9 ]
    geom.face(27).n_poi = 6; allocate(geom.face(27).poi(6)); geom.face(27).poi(1:6) = [  48, 107, 106,  47,  34,  4 ]
    geom.face(28).n_poi = 6; allocate(geom.face(28).poi(6)); geom.face(28).poi(1:6) = [  49, 108, 107,  48,  26, 10 ]
    geom.face(29).n_poi = 6; allocate(geom.face(29).poi(6)); geom.face(29).poi(1:6) = [  50, 109, 108,  49,  61, 16 ]
    geom.face(30).n_poi = 6; allocate(geom.face(30).poi(6)); geom.face(30).poi(1:6) = [  46, 110, 109,  50,  78, 15 ]
    geom.face(31).n_poi = 6; allocate(geom.face(31).poi(6)); geom.face(31).poi(1:6) = [  52, 111, 115,  51,  44,  8 ]
    geom.face(32).n_poi = 6; allocate(geom.face(32).poi(6)); geom.face(32).poi(1:6) = [  53, 112, 111,  52,  25,  5 ]
    geom.face(33).n_poi = 6; allocate(geom.face(33).poi(6)); geom.face(33).poi(1:6) = [  54, 113, 112,  53,  36, 11 ]
    geom.face(34).n_poi = 6; allocate(geom.face(34).poi(6)); geom.face(34).poi(1:6) = [  55, 114, 113,  54,  67, 17 ]
    geom.face(35).n_poi = 6; allocate(geom.face(35).poi(6)); geom.face(35).poi(1:6) = [  51, 115, 114,  55,  74, 14 ]
    geom.face(36).n_poi = 6; allocate(geom.face(36).poi(6)); geom.face(36).poi(1:6) = [  57, 116, 120,  56,  69, 12 ]
    geom.face(37).n_poi = 6; allocate(geom.face(37).poi(6)); geom.face(37).poi(1:6) = [  58, 117, 116,  57,  39,  6 ]
    geom.face(38).n_poi = 6; allocate(geom.face(38).poi(6)); geom.face(38).poi(1:6) = [  59, 118, 117,  58,  35,  9 ]
    geom.face(39).n_poi = 6; allocate(geom.face(39).poi(6)); geom.face(39).poi(1:6) = [  60, 119, 118,  59,  46, 15 ]
    geom.face(40).n_poi = 6; allocate(geom.face(40).poi(6)); geom.face(40).poi(1:6) = [  56, 120, 119,  60,  77, 18 ]
    geom.face(41).n_poi = 6; allocate(geom.face(41).poi(6)); geom.face(41).poi(1:6) = [  62, 121, 125,  61,  49, 10 ]
    geom.face(42).n_poi = 6; allocate(geom.face(42).poi(6)); geom.face(42).poi(1:6) = [  63, 122, 121,  62,  30,  7 ]
    geom.face(43).n_poi = 6; allocate(geom.face(43).poi(6)); geom.face(43).poi(1:6) = [  64, 123, 122,  63,  41, 13 ]
    geom.face(44).n_poi = 6; allocate(geom.face(44).poi(6)); geom.face(44).poi(1:6) = [  65, 124, 123,  64,  72, 19 ]
    geom.face(45).n_poi = 6; allocate(geom.face(45).poi(6)); geom.face(45).poi(1:6) = [  61, 125, 124,  65,  79, 16 ]
    geom.face(46).n_poi = 6; allocate(geom.face(46).poi(6)); geom.face(46).poi(1:6) = [  67, 126, 130,  66,  75, 17 ]
    geom.face(47).n_poi = 6; allocate(geom.face(47).poi(6)); geom.face(47).poi(1:6) = [  68, 127, 126,  67,  54, 11 ]
    geom.face(48).n_poi = 6; allocate(geom.face(48).poi(6)); geom.face(48).poi(1:6) = [  69, 128, 127,  68,  40, 12 ]
    geom.face(49).n_poi = 6; allocate(geom.face(49).poi(6)); geom.face(49).poi(1:6) = [  70, 129, 128,  69,  56, 18 ]
    geom.face(50).n_poi = 6; allocate(geom.face(50).poi(6)); geom.face(50).poi(1:6) = [  66, 130, 129,  70,  76, 20 ]
    geom.face(51).n_poi = 6; allocate(geom.face(51).poi(6)); geom.face(51).poi(1:6) = [  72, 131, 135,  71,  80, 19 ]
    geom.face(52).n_poi = 6; allocate(geom.face(52).poi(6)); geom.face(52).poi(1:6) = [  73, 132, 131,  72,  64, 13 ]
    geom.face(53).n_poi = 6; allocate(geom.face(53).poi(6)); geom.face(53).poi(1:6) = [  74, 133, 132,  73,  45, 14 ]
    geom.face(54).n_poi = 6; allocate(geom.face(54).poi(6)); geom.face(54).poi(1:6) = [  75, 134, 133,  74,  55, 17 ]
    geom.face(55).n_poi = 6; allocate(geom.face(55).poi(6)); geom.face(55).poi(1:6) = [  71, 135, 134,  75,  66, 20 ]
    geom.face(56).n_poi = 6; allocate(geom.face(56).poi(6)); geom.face(56).poi(1:6) = [  77, 136, 140,  76,  70, 18 ]
    geom.face(57).n_poi = 6; allocate(geom.face(57).poi(6)); geom.face(57).poi(1:6) = [  78, 137, 136,  77,  60, 15 ]
    geom.face(58).n_poi = 6; allocate(geom.face(58).poi(6)); geom.face(58).poi(1:6) = [  79, 138, 137,  78,  50, 16 ]
    geom.face(59).n_poi = 6; allocate(geom.face(59).poi(6)); geom.face(59).poi(1:6) = [  80, 139, 138,  79,  65, 19 ]
    geom.face(60).n_poi = 6; allocate(geom.face(60).poi(6)); geom.face(60).poi(1:6) = [  76, 140, 139,  80,  71, 20 ]
    geom.face(61).n_poi = 5; allocate(geom.face(61).poi(5)); geom.face(61).poi(1:5) = [  82,  83,  84,  85,  81 ]
    geom.face(62).n_poi = 5; allocate(geom.face(62).poi(5)); geom.face(62).poi(1:5) = [  87,  88,  89,  90,  86 ]
    geom.face(63).n_poi = 5; allocate(geom.face(63).poi(5)); geom.face(63).poi(1:5) = [  92,  93,  94,  95,  91 ]
    geom.face(64).n_poi = 5; allocate(geom.face(64).poi(5)); geom.face(64).poi(1:5) = [  97,  98,  99, 100,  96 ]
    geom.face(65).n_poi = 5; allocate(geom.face(65).poi(5)); geom.face(65).poi(1:5) = [ 102, 103, 104, 105, 101 ]
    geom.face(66).n_poi = 5; allocate(geom.face(66).poi(5)); geom.face(66).poi(1:5) = [ 107, 108, 109, 110, 106 ]
    geom.face(67).n_poi = 5; allocate(geom.face(67).poi(5)); geom.face(67).poi(1:5) = [ 112, 113, 114, 115, 111 ]
    geom.face(68).n_poi = 5; allocate(geom.face(68).poi(5)); geom.face(68).poi(1:5) = [ 117, 118, 119, 120, 116 ]
    geom.face(69).n_poi = 5; allocate(geom.face(69).poi(5)); geom.face(69).poi(1:5) = [ 122, 123, 124, 125, 121 ]
    geom.face(70).n_poi = 5; allocate(geom.face(70).poi(5)); geom.face(70).poi(1:5) = [ 127, 128, 129, 130, 126 ]
    geom.face(71).n_poi = 5; allocate(geom.face(71).poi(5)); geom.face(71).poi(1:5) = [ 132, 133, 134, 135, 131 ]
    geom.face(72).n_poi = 5; allocate(geom.face(72).poi(5)); geom.face(72).poi(1:5) = [ 137, 138, 139, 140, 136 ]
end subroutine Exam_Miscellaneous_Goldberg_Dk5dgD

! ---------------------------------------------------------------------------------------

! Example of Double Helix
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Double_Helix(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "41_Double_Helix"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Double helix"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.88d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 36
    geom.n_face = 66

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-80.0000d0,  -0.0000d0, -26.1313d0 ]
    geom.iniP( 2).pos(1:3) = [-80.0000d0,  18.4776d0, -18.4776d0 ]
    geom.iniP( 3).pos(1:3) = [-80.0000d0,  -0.0000d0,  26.1313d0 ]
    geom.iniP( 4).pos(1:3) = [-80.0000d0, -18.4776d0,  18.4776d0 ]
    geom.iniP( 5).pos(1:3) = [-60.0000d0,  18.4776d0, -18.4776d0 ]
    geom.iniP( 6).pos(1:3) = [-60.0000d0,  26.1313d0,   0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [-60.0000d0, -18.4776d0,  18.4776d0 ]
    geom.iniP( 8).pos(1:3) = [-60.0000d0, -26.1313d0,   0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [-40.0000d0,  26.1313d0,   0.0000d0 ]
    geom.iniP(10).pos(1:3) = [-40.0000d0,  18.4776d0,  18.4776d0 ]
    geom.iniP(11).pos(1:3) = [-40.0000d0, -26.1313d0,   0.0000d0 ]
    geom.iniP(12).pos(1:3) = [-40.0000d0, -18.4776d0, -18.4776d0 ]
    geom.iniP(13).pos(1:3) = [-20.0000d0,  18.4776d0,  18.4776d0 ]
    geom.iniP(14).pos(1:3) = [-20.0000d0,  -0.0000d0,  26.1313d0 ]
    geom.iniP(15).pos(1:3) = [-20.0000d0, -18.4776d0, -18.4776d0 ]
    geom.iniP(16).pos(1:3) = [-20.0000d0,  -0.0000d0, -26.1313d0 ]
    geom.iniP(17).pos(1:3) = [  0.0000d0,  -0.0000d0,  26.1313d0 ]
    geom.iniP(18).pos(1:3) = [  0.0000d0, -18.4776d0,  18.4776d0 ]
    geom.iniP(19).pos(1:3) = [  0.0000d0,  -0.0000d0, -26.1313d0 ]
    geom.iniP(20).pos(1:3) = [  0.0000d0,  18.4776d0, -18.4776d0 ]
    geom.iniP(21).pos(1:3) = [ 20.0000d0, -18.4776d0,  18.4776d0 ]
    geom.iniP(22).pos(1:3) = [ 20.0000d0, -26.1313d0,   0.0000d0 ]
    geom.iniP(23).pos(1:3) = [ 20.0000d0,  18.4776d0, -18.4776d0 ]
    geom.iniP(24).pos(1:3) = [ 20.0000d0,  26.1313d0,   0.0000d0 ]
    geom.iniP(25).pos(1:3) = [ 40.0000d0, -26.1313d0,   0.0000d0 ]
    geom.iniP(26).pos(1:3) = [ 40.0000d0, -18.4776d0, -18.4776d0 ]
    geom.iniP(27).pos(1:3) = [ 40.0000d0,  26.1313d0,   0.0000d0 ]
    geom.iniP(28).pos(1:3) = [ 40.0000d0,  18.4776d0,  18.4776d0 ]
    geom.iniP(29).pos(1:3) = [ 60.0000d0, -18.4776d0, -18.4776d0 ]
    geom.iniP(30).pos(1:3) = [ 60.0000d0,  -0.0000d0, -26.1313d0 ]
    geom.iniP(31).pos(1:3) = [ 60.0000d0,  18.4776d0,  18.4776d0 ]
    geom.iniP(32).pos(1:3) = [ 60.0000d0,  -0.0000d0,  26.1313d0 ]
    geom.iniP(33).pos(1:3) = [ 80.0000d0,  -0.0000d0, -26.1313d0 ]
    geom.iniP(34).pos(1:3) = [ 80.0000d0,  18.4776d0, -18.4776d0 ]
    geom.iniP(35).pos(1:3) = [ 80.0000d0,  -0.0000d0,  26.1313d0 ]
    geom.iniP(36).pos(1:3) = [ 80.0000d0, -18.4776d0,  18.4776d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  5,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  5,  2,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  6,  2,  3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  3,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  7,  3,  4 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  7,  4,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  8,  4,  5 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  4,  1,  5 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  9,  5,  6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  9,  6, 10 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 10,  6,  7 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 10,  7, 11 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 11,  7,  8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 11,  8, 12 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 12,  8,  9 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  8,  5,  9 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 13,  9, 10 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 13, 10, 14 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 14, 10, 11 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 14, 11, 15 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 15, 11, 12 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 15, 12, 16 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 16, 12, 13 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 12,  9, 13 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 17, 13, 14 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 17, 14, 18 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 18, 14, 15 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 18, 15, 19 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 19, 15, 16 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 19, 16, 20 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 20, 16, 17 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 16, 13, 17 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 21, 17, 18 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 21, 18, 22 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 22, 18, 19 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 22, 19, 23 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 23, 19, 20 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 23, 20, 24 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 24, 20, 21 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 20, 17, 21 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 25, 21, 22 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 25, 22, 26 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 26, 22, 23 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 26, 23, 27 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 27, 23, 24 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 27, 24, 28 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 28, 24, 25 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 24, 21, 25 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 29, 25, 26 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 29, 26, 30 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 30, 26, 27 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 30, 27, 31 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 31, 27, 28 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 31, 28, 32 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 32, 28, 29 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 28, 25, 29 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 33, 29, 30 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 33, 30, 34 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 34, 30, 31 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 34, 31, 35 ]
    geom.face(61).n_poi = 3; allocate(geom.face(61).poi(3)); geom.face(61).poi(1:3) = [ 35, 31, 32 ]
    geom.face(62).n_poi = 3; allocate(geom.face(62).poi(3)); geom.face(62).poi(1:3) = [ 35, 32, 36 ]
    geom.face(63).n_poi = 3; allocate(geom.face(63).poi(3)); geom.face(63).poi(1:3) = [ 36, 32, 33 ]
    geom.face(64).n_poi = 3; allocate(geom.face(64).poi(3)); geom.face(64).poi(1:3) = [ 32, 29, 33 ]
    geom.face(65).n_poi = 4; allocate(geom.face(65).poi(4)); geom.face(65).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face(66).n_poi = 4; allocate(geom.face(66).poi(4)); geom.face(66).poi(1:4) = [ 33, 34, 35, 36 ]
end subroutine Exam_Miscellaneous_Double_Helix

! ---------------------------------------------------------------------------------------

! Example of Nested Cube
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Nested_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "42_Nested_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Nested cube"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 0.7d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 16
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ 11.5470d0,  11.5470d0,  11.5470d0 ]
    geom.iniP( 2).pos(1:3) = [-11.5470d0,  11.5470d0,  11.5470d0 ]
    geom.iniP( 3).pos(1:3) = [-11.5470d0, -11.5470d0,  11.5470d0 ]
    geom.iniP( 4).pos(1:3) = [ 11.5470d0, -11.5470d0,  11.5470d0 ]
    geom.iniP( 5).pos(1:3) = [ 11.5470d0, -11.5470d0, -11.5470d0 ]
    geom.iniP( 6).pos(1:3) = [ 11.5470d0,  11.5470d0, -11.5470d0 ]
    geom.iniP( 7).pos(1:3) = [-11.5470d0,  11.5470d0, -11.5470d0 ]
    geom.iniP( 8).pos(1:3) = [-11.5470d0, -11.5470d0, -11.5470d0 ]
    geom.iniP( 9).pos(1:3) = [ 23.0940d0,  23.0940d0,  23.0940d0 ]
    geom.iniP(10).pos(1:3) = [-23.0940d0,  23.0940d0,  23.0940d0 ]
    geom.iniP(11).pos(1:3) = [-23.0940d0, -23.0940d0,  23.0940d0 ]
    geom.iniP(12).pos(1:3) = [ 23.0940d0, -23.0940d0,  23.0940d0 ]
    geom.iniP(13).pos(1:3) = [ 23.0940d0, -23.0940d0, -23.0940d0 ]
    geom.iniP(14).pos(1:3) = [ 23.0940d0,  23.0940d0, -23.0940d0 ]
    geom.iniP(15).pos(1:3) = [-23.0940d0,  23.0940d0, -23.0940d0 ]
    geom.iniP(16).pos(1:3) = [-23.0940d0, -23.0940d0, -23.0940d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  1,  2,  7,  6 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  2,  3,  8,  7 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  3,  4,  5,  8 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  4,  1,  6,  5 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  9, 14, 15, 10 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 10, 15, 16, 11 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 11, 16, 13, 12 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 12, 13, 14,  9 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  2,  1,  9, 10 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  3,  2, 10, 11 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  4,  3, 11, 12 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  1,  4, 12,  9 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  5,  6, 14, 13 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  6,  7, 15, 14 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  7,  8, 16, 15 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  8,  5, 13, 16 ]
end subroutine Exam_Miscellaneous_Nested_Cube

! ---------------------------------------------------------------------------------------

! Example of Nested Octahedron
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Nested_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "43_Nested_Octa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Nested octahedron"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.05d0     ! Cylindrical model
    prob.move_x   =-1.5d0      ! Cylindrical model
    prob.move_y   =-1.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 12
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.0000d0,   0.0000d0,  20.0000d0 ]
    geom.iniP( 2).pos(1:3) = [  0.0000d0,   0.0000d0, -20.0000d0 ]
    geom.iniP( 3).pos(1:3) = [  0.0000d0,  20.0000d0,   0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [  0.0000d0, -20.0000d0,   0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ 20.0000d0,   0.0000d0,   0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [-20.0000d0,   0.0000d0,   0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [  0.0000d0,   0.0000d0,  40.0000d0 ]
    geom.iniP( 8).pos(1:3) = [  0.0000d0,   0.0000d0, -40.0000d0 ]
    geom.iniP( 9).pos(1:3) = [  0.0000d0,  40.0000d0,   0.0000d0 ]
    geom.iniP(10).pos(1:3) = [  0.0000d0, -40.0000d0,   0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ 40.0000d0,   0.0000d0,   0.0000d0 ]
    geom.iniP(12).pos(1:3) = [-40.0000d0,   0.0000d0,   0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  5,  4 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  4,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  6,  3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  5,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  4,  5 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 2,  3,  6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 7, 10, 11 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 7, 12, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 7,  9, 12 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8,  9, 11 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 8, 11, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 8, 12,  9 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 1,  3,  9,  7 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 5,  1,  7, 11 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 3,  5, 11,  9 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 2,  6, 12,  8 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 4,  2,  8, 10 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 6,  4, 10, 12 ]
end subroutine Exam_Miscellaneous_Nested_Octahedron

! ---------------------------------------------------------------------------------------

! Example of Torus
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Torus(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "44_Torus"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Torus"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.96d0     ! Cylindrical model
    prob.move_x   = 15.0d0     ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 36
    geom.n_face = 72

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-10.0000d0,  17.3206d0,  -0.0001d0 ]
    geom.iniP( 2).pos(1:3) = [-20.0001d0,   0.0000d0,  -0.0001d0 ]
    geom.iniP( 3).pos(1:3) = [-30.0001d0,   0.0000d0,  17.3208d0 ]
    geom.iniP( 4).pos(1:3) = [ 15.0001d0, -25.9809d0,  17.3208d0 ]
    geom.iniP( 5).pos(1:3) = [ 20.0001d0,  -0.0000d0,  -0.0001d0 ]
    geom.iniP( 6).pos(1:3) = [ 30.0001d0,  -0.0000d0,  17.3208d0 ]
    geom.iniP( 7).pos(1:3) = [ 15.0001d0,  25.9809d0,  17.3208d0 ]
    geom.iniP( 8).pos(1:3) = [-15.0001d0,  25.9809d0,  17.3208d0 ]
    geom.iniP( 9).pos(1:3) = [-15.0001d0, -25.9809d0,  17.3208d0 ]
    geom.iniP(10).pos(1:3) = [ 50.0002d0,  -0.0000d0,  17.3208d0 ]
    geom.iniP(11).pos(1:3) = [ 25.0001d0,  43.3014d0,  17.3208d0 ]
    geom.iniP(12).pos(1:3) = [-60.0003d0,   0.0000d0,  -0.0001d0 ]
    geom.iniP(13).pos(1:3) = [-25.0001d0,  43.3014d0,  17.3208d0 ]
    geom.iniP(14).pos(1:3) = [-30.0001d0,  51.9619d0,  -0.0001d0 ]
    geom.iniP(15).pos(1:3) = [-50.0002d0,   0.0000d0,  17.3208d0 ]
    geom.iniP(16).pos(1:3) = [-25.0001d0, -43.3014d0,  17.3208d0 ]
    geom.iniP(17).pos(1:3) = [ 25.0001d0, -43.3014d0,  17.3208d0 ]
    geom.iniP(18).pos(1:3) = [ 30.0001d0,  51.9619d0,  -0.0001d0 ]
    geom.iniP(19).pos(1:3) = [-50.0002d0,   0.0000d0, -17.3206d0 ]
    geom.iniP(20).pos(1:3) = [ 30.0001d0, -51.9619d0,  -0.0001d0 ]
    geom.iniP(21).pos(1:3) = [ 25.0001d0, -43.3014d0, -17.3206d0 ]
    geom.iniP(22).pos(1:3) = [-30.0001d0, -51.9619d0,  -0.0001d0 ]
    geom.iniP(23).pos(1:3) = [-25.0001d0, -43.3014d0, -17.3206d0 ]
    geom.iniP(24).pos(1:3) = [ 50.0002d0,  -0.0000d0, -17.3206d0 ]
    geom.iniP(25).pos(1:3) = [ 25.0001d0,  43.3014d0, -17.3206d0 ]
    geom.iniP(26).pos(1:3) = [ 60.0003d0,  -0.0000d0,  -0.0001d0 ]
    geom.iniP(27).pos(1:3) = [-25.0001d0,  43.3014d0, -17.3206d0 ]
    geom.iniP(28).pos(1:3) = [-30.0001d0,   0.0000d0, -17.3206d0 ]
    geom.iniP(29).pos(1:3) = [ 30.0001d0,  -0.0000d0, -17.3206d0 ]
    geom.iniP(30).pos(1:3) = [ 15.0001d0,  25.9809d0, -17.3206d0 ]
    geom.iniP(31).pos(1:3) = [-15.0001d0,  25.9809d0, -17.3206d0 ]
    geom.iniP(32).pos(1:3) = [-15.0001d0, -25.9809d0, -17.3206d0 ]
    geom.iniP(33).pos(1:3) = [ 10.0000d0, -17.3206d0,  -0.0001d0 ]
    geom.iniP(34).pos(1:3) = [-10.0000d0, -17.3206d0,  -0.0001d0 ]
    geom.iniP(35).pos(1:3) = [ 15.0001d0, -25.9809d0, -17.3206d0 ]
    geom.iniP(36).pos(1:3) = [ 10.0000d0,  17.3206d0,  -0.0001d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1, 36,  8 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  8, 36,  7 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  1,  3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  3,  1,  8 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 34,  2,  9 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  9,  2,  3 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 33, 34,  4 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  4, 34,  9 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  5, 33,  6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  6, 33,  4 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 36,  5,  7 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  7,  5,  6 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  8,  7, 13 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 13,  7, 11 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  3,  8, 15 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 15,  8, 13 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  9,  3, 16 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 16,  3, 15 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  4,  9, 17 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 17,  9, 16 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [  6,  4, 10 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 10,  4, 17 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  7,  6, 11 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 11,  6, 10 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 13, 11, 14 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 14, 11, 18 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 15, 13, 12 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 12, 13, 14 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 16, 15, 22 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 22, 15, 12 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 17, 16, 20 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 20, 16, 22 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 10, 17, 26 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 26, 17, 20 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 11, 10, 18 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 18, 10, 26 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 14, 18, 27 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 27, 18, 25 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 12, 14, 19 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 19, 14, 27 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 22, 12, 23 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 23, 12, 19 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 20, 22, 21 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 21, 22, 23 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 26, 20, 24 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 24, 20, 21 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 18, 26, 25 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 25, 26, 24 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 27, 25, 31 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 31, 25, 30 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 19, 27, 28 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 28, 27, 31 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 23, 19, 32 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 32, 19, 28 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 21, 23, 35 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 35, 23, 32 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 24, 21, 29 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 29, 21, 35 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 25, 24, 30 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 30, 24, 29 ]
    geom.face(61).n_poi = 3; allocate(geom.face(61).poi(3)); geom.face(61).poi(1:3) = [ 31, 30,  1 ]
    geom.face(62).n_poi = 3; allocate(geom.face(62).poi(3)); geom.face(62).poi(1:3) = [  1, 30, 36 ]
    geom.face(63).n_poi = 3; allocate(geom.face(63).poi(3)); geom.face(63).poi(1:3) = [ 28, 31,  2 ]
    geom.face(64).n_poi = 3; allocate(geom.face(64).poi(3)); geom.face(64).poi(1:3) = [  2, 31,  1 ]
    geom.face(65).n_poi = 3; allocate(geom.face(65).poi(3)); geom.face(65).poi(1:3) = [ 32, 28, 34 ]
    geom.face(66).n_poi = 3; allocate(geom.face(66).poi(3)); geom.face(66).poi(1:3) = [ 34, 28,  2 ]
    geom.face(67).n_poi = 3; allocate(geom.face(67).poi(3)); geom.face(67).poi(1:3) = [ 35, 32, 33 ]
    geom.face(68).n_poi = 3; allocate(geom.face(68).poi(3)); geom.face(68).poi(1:3) = [ 33, 32, 34 ]
    geom.face(69).n_poi = 3; allocate(geom.face(69).poi(3)); geom.face(69).poi(1:3) = [ 29, 35,  5 ]
    geom.face(70).n_poi = 3; allocate(geom.face(70).poi(3)); geom.face(70).poi(1:3) = [  5, 35, 33 ]
    geom.face(71).n_poi = 3; allocate(geom.face(71).poi(3)); geom.face(71).poi(1:3) = [ 30, 29, 36 ]
    geom.face(72).n_poi = 3; allocate(geom.face(72).poi(3)); geom.face(72).poi(1:3) = [ 36, 29,  5 ]
end subroutine Exam_Miscellaneous_Torus

! ---------------------------------------------------------------------------------------

! Example of Double Torus
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Miscellaneous_Double_Torus(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "45_Double_Torus"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Double torus"

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.9d0      ! Cylindrical model
    prob.move_x   = 2.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 44
    geom.n_face = 46

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ 51.9615d0,  10.0000d0,  20.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ 34.6410d0,  20.0000d0,  20.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ 17.3205d0,  10.0000d0,  20.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ 17.3205d0, -10.0000d0,  20.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ 34.6410d0, -20.0000d0,  20.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ 51.9615d0, -10.0000d0,  20.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ 51.9615d0,  10.0000d0, -20.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ 34.6410d0,  20.0000d0, -20.0000d0 ]
    geom.iniP( 9).pos(1:3) = [ 17.3205d0,  10.0000d0, -20.0000d0 ]
    geom.iniP(10).pos(1:3) = [ 17.3205d0, -10.0000d0, -20.0000d0 ]
    geom.iniP(11).pos(1:3) = [ 34.6410d0, -20.0000d0, -20.0000d0 ]
    geom.iniP(12).pos(1:3) = [ 51.9615d0, -10.0000d0, -20.0000d0 ]
    geom.iniP(13).pos(1:3) = [ 69.2821d0,  20.0000d0,  20.0000d0 ]
    geom.iniP(14).pos(1:3) = [ 34.6410d0,  40.0000d0,  20.0000d0 ]
    geom.iniP(15).pos(1:3) = [  0.0000d0,  20.0000d0,  20.0000d0 ]
    geom.iniP(16).pos(1:3) = [  0.0000d0, -20.0000d0,  20.0000d0 ]
    geom.iniP(17).pos(1:3) = [ 34.6410d0, -40.0000d0,  20.0000d0 ]
    geom.iniP(18).pos(1:3) = [ 69.2821d0, -20.0000d0,  20.0000d0 ]
    geom.iniP(19).pos(1:3) = [ 69.2821d0,  20.0000d0, -20.0000d0 ]
    geom.iniP(20).pos(1:3) = [ 34.6410d0,  40.0000d0, -20.0000d0 ]
    geom.iniP(21).pos(1:3) = [  0.0000d0,  20.0000d0, -20.0000d0 ]
    geom.iniP(22).pos(1:3) = [  0.0000d0, -20.0000d0, -20.0000d0 ]
    geom.iniP(23).pos(1:3) = [ 34.6410d0, -40.0000d0, -20.0000d0 ]
    geom.iniP(24).pos(1:3) = [ 69.2821d0, -20.0000d0, -20.0000d0 ]
    geom.iniP(25).pos(1:3) = [-17.3205d0,  20.0000d0, -10.0000d0 ]
    geom.iniP(26).pos(1:3) = [-34.6410d0,  20.0000d0, -20.0000d0 ]
    geom.iniP(27).pos(1:3) = [-51.9615d0,  20.0000d0, -10.0000d0 ]
    geom.iniP(28).pos(1:3) = [-51.9615d0,  20.0000d0,  10.0000d0 ]
    geom.iniP(29).pos(1:3) = [-34.6410d0,  20.0000d0,  20.0000d0 ]
    geom.iniP(30).pos(1:3) = [-17.3205d0,  20.0000d0,  10.0000d0 ]
    geom.iniP(31).pos(1:3) = [-17.3205d0, -20.0000d0, -10.0000d0 ]
    geom.iniP(32).pos(1:3) = [-34.6410d0, -20.0000d0, -20.0000d0 ]
    geom.iniP(33).pos(1:3) = [-51.9615d0, -20.0000d0, -10.0000d0 ]
    geom.iniP(34).pos(1:3) = [-51.9615d0, -20.0000d0,  10.0000d0 ]
    geom.iniP(35).pos(1:3) = [-34.6410d0, -20.0000d0,  20.0000d0 ]
    geom.iniP(36).pos(1:3) = [-17.3205d0, -20.0000d0,  10.0000d0 ]
    geom.iniP(37).pos(1:3) = [-34.6410d0,  20.0000d0, -40.0000d0 ]
    geom.iniP(38).pos(1:3) = [-69.2821d0,  20.0000d0, -20.0000d0 ]
    geom.iniP(39).pos(1:3) = [-69.2821d0,  20.0000d0,  20.0000d0 ]
    geom.iniP(40).pos(1:3) = [-34.6410d0,  20.0000d0,  40.0000d0 ]
    geom.iniP(41).pos(1:3) = [-34.6410d0, -20.0000d0, -40.0000d0 ]
    geom.iniP(42).pos(1:3) = [-69.2821d0, -20.0000d0, -20.0000d0 ]
    geom.iniP(43).pos(1:3) = [-69.2821d0, -20.0000d0,  20.0000d0 ]
    geom.iniP(44).pos(1:3) = [-34.6410d0, -20.0000d0,  40.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 13, 14,  2,  1 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 14, 15,  3,  2 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 15, 16,  4,  3 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 16, 17,  5,  4 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 17, 18,  6,  5 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 18, 13,  1,  6 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  7,  8, 20, 19 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  8,  9, 21, 20 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  9, 10, 22, 21 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 10, 11, 23, 22 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 11, 12, 24, 23 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 12,  7, 19, 24 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  1,  2,  8,  7 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  2,  3,  9,  8 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  3,  4, 10,  9 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  4,  5, 11, 10 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  5,  6, 12, 11 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  6,  1,  7, 12 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 14, 13, 19, 20 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 15, 14, 20, 21 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 17, 16, 22, 23 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 18, 17, 23, 24 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 13, 18, 24, 19 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 21, 37, 26, 25 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 37, 38, 27, 26 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 38, 39, 28, 27 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 39, 40, 29, 28 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 40, 15, 30, 29 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 15, 21, 25, 30 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 31, 32, 41, 22 ]
    geom.face(31).n_poi = 4; allocate(geom.face(31).poi(4)); geom.face(31).poi(1:4) = [ 32, 33, 42, 41 ]
    geom.face(32).n_poi = 4; allocate(geom.face(32).poi(4)); geom.face(32).poi(1:4) = [ 33, 34, 43, 42 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [ 34, 35, 44, 43 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [ 35, 36, 16, 44 ]
    geom.face(35).n_poi = 4; allocate(geom.face(35).poi(4)); geom.face(35).poi(1:4) = [ 36, 31, 22, 16 ]
    geom.face(36).n_poi = 4; allocate(geom.face(36).poi(4)); geom.face(36).poi(1:4) = [ 25, 26, 32, 31 ]
    geom.face(37).n_poi = 4; allocate(geom.face(37).poi(4)); geom.face(37).poi(1:4) = [ 26, 27, 33, 32 ]
    geom.face(38).n_poi = 4; allocate(geom.face(38).poi(4)); geom.face(38).poi(1:4) = [ 27, 28, 34, 33 ]
    geom.face(39).n_poi = 4; allocate(geom.face(39).poi(4)); geom.face(39).poi(1:4) = [ 28, 29, 35, 34 ]
    geom.face(40).n_poi = 4; allocate(geom.face(40).poi(4)); geom.face(40).poi(1:4) = [ 29, 30, 36, 35 ]
    geom.face(41).n_poi = 4; allocate(geom.face(41).poi(4)); geom.face(41).poi(1:4) = [ 30, 25, 31, 36 ]
    geom.face(42).n_poi = 4; allocate(geom.face(42).poi(4)); geom.face(42).poi(1:4) = [ 37, 21, 22, 41 ]
    geom.face(43).n_poi = 4; allocate(geom.face(43).poi(4)); geom.face(43).poi(1:4) = [ 38, 37, 41, 42 ]
    geom.face(44).n_poi = 4; allocate(geom.face(44).poi(4)); geom.face(44).poi(1:4) = [ 39, 38, 42, 43 ]
    geom.face(45).n_poi = 4; allocate(geom.face(45).poi(4)); geom.face(45).poi(1:4) = [ 40, 39, 43, 44 ]
    geom.face(46).n_poi = 4; allocate(geom.face(46).poi(4)); geom.face(46).poi(1:4) = [ 15, 40, 44, 16 ]
end subroutine Exam_Miscellaneous_Double_Torus

! ---------------------------------------------------------------------------------------

end module Exam_Miscellaneous