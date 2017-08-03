!
! ---------------------------------------------------------------------------------------
!
!                                       Exam_Prism
!
!                                                                    Updated : 2017/03/27
!
! Comments: Prism and antiprism
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
! Reference: http://www.georgehart.com/virtual-polyhedra/prisms-index.html
!
module Exam_Prism

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    ! Prism
    public Exam_Prism_Triangular          ! V= 6, E= 9, F= 5
    public Exam_Prism_Pentagonal          ! V=10, E=15, F= 7
    public Exam_Prism_Hexagonal           ! V=12, E=18, F= 8
    public Exam_Prism_Heptagonal          ! V=14, E=21, F= 9
    public Exam_Prism_Octagonal           ! V=16, E=24, F=10
    public Exam_Prism_Enneagonal          ! V=18, E=27, F=11
    public Exam_Prism_Decagonal           ! V=20, E=30, F=12

    ! Antiprism
    public Exam_Antiprism_Square          ! V= 8, E=16, F=10
    public Exam_Antiprism_Pentagonal      ! V=10, E=20, F=12
    public Exam_Antiprism_Hexagonal       ! V=12, E=24, F=14
    public Exam_Antiprism_Heptagonal      ! V=14, E=28, F=16
    public Exam_Antiprism_Octagonal       ! V=16, E=32, F=18
    public Exam_Antiprism_Enneagonal      ! V=18, E=36, F=20
    public Exam_Antiprism_Decagonal       ! V=20, E=40, F=22

contains

! ---------------------------------------------------------------------------------------

! Example of triangular prism
subroutine Exam_Prism_Triangular(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Triangular"
    prob.name_file = "Prism Triangular"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 6
    geom.n_face = 5

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  -0.0000d0,   0.0000d0,  15.2753d0 ]
    geom.iniP(2).pos(1:3) = [  15.1186d0,   0.0000d0,   2.1822d0 ]
    geom.iniP(3).pos(1:3) = [ -11.3389d0,  10.0000d0,   2.1822d0 ]
    geom.iniP(4).pos(1:3) = [   1.8898d0, -15.0000d0,   2.1822d0 ]
    geom.iniP(5).pos(1:3) = [   3.7796d0,  10.0000d0, -10.9109d0 ]
    geom.iniP(6).pos(1:3) = [  -9.4491d0,  -5.0000d0, -10.9109d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 4, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 3, 5, 6 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 1, 2, 5, 3 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 3, 6, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 2, 4, 6, 5 ]
end subroutine Exam_Prism_Triangular

! ---------------------------------------------------------------------------------------

! Example of pentagonal prism
subroutine Exam_Prism_Pentagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Pentagonal"
    prob.name_file = "Prism_Pentagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 10
    geom.n_face = 7

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  19.7343d0 ]
    geom.iniP( 2).pos(1:3) = [  17.2421d0,  -0.0000d0,   9.5997d0 ]
    geom.iniP( 3).pos(1:3) = [  -5.9570d0,  16.1804d0,   9.5997d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.1260d0, -11.1804d0,   9.5997d0 ]
    geom.iniP( 5).pos(1:3) = [  11.2851d0,  16.1804d0,  -0.5350d0 ]
    geom.iniP( 6).pos(1:3) = [  14.7724d0, -11.1804d0,  -6.7985d0 ]
    geom.iniP( 7).pos(1:3) = [ -19.0829d0,   5.0000d0,  -0.5350d0 ]
    geom.iniP( 8).pos(1:3) = [  -3.9961d0, -18.0902d0,  -6.7985d0 ]
    geom.iniP( 9).pos(1:3) = [   8.8154d0,   5.0000d0, -16.9332d0 ]
    geom.iniP(10).pos(1:3) = [  -9.9531d0,  -1.9098d0, -16.9332d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1, 2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1, 3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4, 7, 10,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 8, 10,  9 ]
    geom.face(6).n_poi = 5; allocate(geom.face(6).poi(5)); geom.face(6).poi(1:5) = [ 1, 4,  8,  6, 2 ]
    geom.face(7).n_poi = 5; allocate(geom.face(7).poi(5)); geom.face(7).poi(1:5) = [ 3, 5,  9, 10, 7 ]
end subroutine Exam_Prism_Pentagonal

! ---------------------------------------------------------------------------------------

! Example of hexagonal prism
subroutine Exam_Prism_Hexagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Hexagonal"
    prob.name_file = "Prism_Hexagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  22.3606d0 ]
    geom.iniP( 2).pos(1:3) = [  17.8886d0,  -0.0000d0,  13.4164d0 ]
    geom.iniP( 3).pos(1:3) = [  -4.4721d0,  17.3205d0,  13.4164d0 ]
    geom.iniP( 4).pos(1:3) = [ -15.6525d0,  -8.6603d0,  13.4164d0 ]
    geom.iniP( 5).pos(1:3) = [  13.4164d0,  17.3205d0,   4.4721d0 ]
    geom.iniP( 6).pos(1:3) = [  20.1246d0,  -8.6603d0,  -4.4721d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.1246d0,   8.6603d0,   4.4721d0 ]
    geom.iniP( 8).pos(1:3) = [ -13.4164d0, -17.3205d0,  -4.4721d0 ]
    geom.iniP( 9).pos(1:3) = [  15.6525d0,   8.6603d0, -13.4164d0 ]
    geom.iniP(10).pos(1:3) = [   4.4721d0, -17.3205d0, -13.4164d0 ]
    geom.iniP(11).pos(1:3) = [ -17.8886d0,  -0.0000d0, -13.4164d0 ]
    geom.iniP(12).pos(1:3) = [   0.0000d0,  -0.0000d0, -22.3606d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 10, 12,  9 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 8, 11, 12, 10 ]
    geom.face(7).n_poi = 6; allocate(geom.face(7).poi(6)); geom.face(7).poi(1:6) = [ 1,  4,  8, 10,  6, 2 ]
    geom.face(8).n_poi = 6; allocate(geom.face(8).poi(6)); geom.face(8).poi(1:6) = [ 3,  5,  9, 12, 11, 7 ]
end subroutine Exam_Prism_Hexagonal

! ---------------------------------------------------------------------------------------

! Example of hetagonal prism
subroutine Exam_Prism_Heptagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Heptagonal"
    prob.name_file = "Prism_Heptagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 9

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  25.1236d0 ]
    geom.iniP( 2).pos(1:3) = [  18.3475d0,  -0.0000d0,  17.1630d0 ]
    geom.iniP( 3).pos(1:3) = [  -3.4540d0,  18.0194d0,  17.1630d0 ]
    geom.iniP( 4).pos(1:3) = [ -17.0470d0,  -6.7845d0,  17.1630d0 ]
    geom.iniP( 5).pos(1:3) = [  14.8935d0,  18.0194d0,   9.2023d0 ]
    geom.iniP( 6).pos(1:3) = [  24.1794d0,  -6.7845d0,  -0.7245d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.5010d0,  11.2349d0,   9.2023d0 ]
    geom.iniP( 8).pos(1:3) = [ -19.9568d0, -15.2446d0,  -0.7245d0 ]
    geom.iniP( 9).pos(1:3) = [  20.7254d0,  11.2349d0,  -8.6852d0 ]
    geom.iniP(10).pos(1:3) = [  13.1042d0, -15.2446d0, -15.0691d0 ]
    geom.iniP(11).pos(1:3) = [ -23.4107d0,   2.7748d0,  -8.6852d0 ]
    geom.iniP(12).pos(1:3) = [  -6.5383d0, -19.0097d0, -15.0691d0 ]
    geom.iniP(13).pos(1:3) = [   9.6502d0,   2.7748d0, -23.0298d0 ]
    geom.iniP(14).pos(1:3) = [  -9.9923d0,  -0.9903d0, -23.0298d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 8, 11, 14, 12 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [10, 12, 14, 13 ]
    geom.face(8).n_poi = 7; allocate(geom.face(8).poi(7)); geom.face(8).poi(1:7) = [ 1,  4,  8, 12, 10,  6, 2 ]
    geom.face(9).n_poi = 7; allocate(geom.face(9).poi(7)); geom.face(9).poi(1:7) = [ 3,  5,  9, 13, 14, 11, 7 ]
end subroutine Exam_Prism_Heptagonal

! ---------------------------------------------------------------------------------------

! Example of octagonal prism
subroutine Exam_Prism_Octagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Octagonal"
    prob.name_file = "16_Prism_Octagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xyz2", 1.0d0, 1.02d0, 0.0d0, 1.0d0)

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  27.9793d0 ]
    geom.iniP( 2).pos(1:3) = [  18.6790d0,  -0.0000d0,  20.8312d0 ]
    geom.iniP( 3).pos(1:3) = [  -2.7355d0,  18.4776d0,  20.8312d0 ]
    geom.iniP( 4).pos(1:3) = [ -17.8778d0,  -5.4120d0,  20.8312d0 ]
    geom.iniP( 5).pos(1:3) = [  15.9435d0,  18.4776d0,  13.6831d0 ]
    geom.iniP( 6).pos(1:3) = [  27.2173d0,  -5.4120d0,   3.5741d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.6133d0,  13.0657d0,  13.6831d0 ]
    geom.iniP( 8).pos(1:3) = [ -24.4818d0, -13.0657d0,   3.5741d0 ]
    geom.iniP( 9).pos(1:3) = [  24.4818d0,  13.0657d0,  -3.5741d0 ]
    geom.iniP(10).pos(1:3) = [  20.6133d0, -13.0657d0, -13.6831d0 ]
    geom.iniP(11).pos(1:3) = [ -27.2173d0,   5.4120d0,  -3.5741d0 ]
    geom.iniP(12).pos(1:3) = [ -15.9435d0, -18.4776d0, -13.6831d0 ]
    geom.iniP(13).pos(1:3) = [  17.8778d0,   5.4120d0, -20.8312d0 ]
    geom.iniP(14).pos(1:3) = [   2.7355d0, -18.4776d0, -20.8312d0 ]
    geom.iniP(15).pos(1:3) = [ -18.6790d0,  -0.0000d0, -20.8312d0 ]
    geom.iniP(16).pos(1:3) = [   0.0000d0,  -0.0000d0, -27.9793d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [10, 14, 16, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [12, 15, 16, 14 ]
    geom.face( 9).n_poi = 8; allocate(geom.face( 9).poi(8)); geom.face( 9).poi(1:8) = [ 1,  4,  8, 12, 14, 10,  6, 2 ]
    geom.face(10).n_poi = 8; allocate(geom.face(10).poi(8)); geom.face(10).poi(1:8) = [ 3,  5,  9, 13, 16, 15, 11, 7 ]
end subroutine Exam_Prism_Octagonal

! ---------------------------------------------------------------------------------------

! Example of enneagonal prism
subroutine Exam_Prism_Enneagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Enneagonal"
    prob.name_file = "17_Prism_Enneagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xyz2", 1.0d0, 1.0d0, -1.0d0, 1.0d0)

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 11

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  30.9008d0 ]
    geom.iniP( 2).pos(1:3) = [  18.9238d0,  -0.0000d0,  24.4286d0 ]
    geom.iniP( 3).pos(1:3) = [  -2.2137d0,  18.7939d0,  24.4286d0 ]
    geom.iniP( 4).pos(1:3) = [ -18.4059d0,  -4.3969d0,  24.4286d0 ]
    geom.iniP( 5).pos(1:3) = [  16.7101d0,  18.7939d0,  17.9562d0 ]
    geom.iniP( 6).pos(1:3) = [  29.5109d0,  -4.3969d0,   8.0401d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.6196d0,  14.3969d0,  17.9562d0 ]
    geom.iniP( 8).pos(1:3) = [ -27.6816d0, -11.1334d0,   8.0401d0 ]
    geom.iniP( 9).pos(1:3) = [  27.2972d0,  14.3969d0,   1.5678d0 ]
    geom.iniP(10).pos(1:3) = [  26.8073d0, -11.1334d0, -10.5962d0 ]
    geom.iniP(11).pos(1:3) = [ -29.8953d0,   7.6605d0,   1.5678d0 ]
    geom.iniP(12).pos(1:3) = [ -23.4868d0, -17.0574d0, -10.5962d0 ]
    geom.iniP(13).pos(1:3) = [  24.5937d0,   7.6605d0, -17.0685d0 ]
    geom.iniP(14).pos(1:3) = [  12.0783d0, -17.0574d0, -22.7602d0 ]
    geom.iniP(15).pos(1:3) = [ -25.7005d0,   1.7365d0, -17.0685d0 ]
    geom.iniP(16).pos(1:3) = [  -7.7844d0, -19.3970d0, -22.7602d0 ]
    geom.iniP(17).pos(1:3) = [   9.8646d0,   1.7365d0, -29.2325d0 ]
    geom.iniP(18).pos(1:3) = [  -9.9981d0,  -0.6031d0, -29.2325d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [10, 14, 17, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [12, 15, 18, 16 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [14, 16, 18, 17 ]
    geom.face(10).n_poi = 9; allocate(geom.face(10).poi(9)); geom.face(10).poi(1:9) = [ 1,  4,  8, 12, 16, 14, 10,  6, 2 ]
    geom.face(11).n_poi = 9; allocate(geom.face(11).poi(9)); geom.face(11).poi(1:9) = [ 3,  5,  9, 13, 17, 18, 15, 11, 7 ]
end subroutine Exam_Prism_Enneagonal

! ---------------------------------------------------------------------------------------

! Example of decagonal prism
subroutine Exam_Prism_Decagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Prism Decagonal"
    prob.name_file = "Prism_Decagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,   0.0000d0,  33.8707d0 ]
    geom.iniP( 2).pos(1:3) = [  19.1085d0,   0.0000d0,  27.9657d0 ]
    geom.iniP( 3).pos(1:3) = [  -1.8247d0,  19.0211d0,  27.9657d0 ]
    geom.iniP( 4).pos(1:3) = [ -18.7600d0,  -3.6327d0,  27.9657d0 ]
    geom.iniP( 5).pos(1:3) = [  17.2838d0,  19.0211d0,  22.0609d0 ]
    geom.iniP( 6).pos(1:3) = [  31.2666d0,  -3.6327d0,  12.5066d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.5847d0,  15.3884d0,  22.0609d0 ]
    geom.iniP( 8).pos(1:3) = [ -30.0058d0,  -9.5106d0,  12.5066d0 ]
    geom.iniP( 9).pos(1:3) = [  29.4419d0,  15.3884d0,   6.6018d0 ]
    geom.iniP(10).pos(1:3) = [  31.8305d0,  -9.5106d0,  -6.6018d0 ]
    geom.iniP(11).pos(1:3) = [ -31.8305d0,   9.5106d0,   6.6018d0 ]
    geom.iniP(12).pos(1:3) = [ -29.4419d0, -15.3884d0,  -6.6018d0 ]
    geom.iniP(13).pos(1:3) = [  30.0058d0,   9.5106d0, -12.5066d0 ]
    geom.iniP(14).pos(1:3) = [  20.5847d0, -15.3884d0, -22.0609d0 ]
    geom.iniP(15).pos(1:3) = [ -31.2666d0,   3.6327d0, -12.5066d0 ]
    geom.iniP(16).pos(1:3) = [ -17.2838d0, -19.0211d0, -22.0609d0 ]
    geom.iniP(17).pos(1:3) = [  18.7600d0,   3.6327d0, -27.9657d0 ]
    geom.iniP(18).pos(1:3) = [   1.8247d0, -19.0211d0, -27.9657d0 ]
    geom.iniP(19).pos(1:3) = [ -19.1085d0,   0.0000d0, -27.9657d0 ]
    geom.iniP(20).pos(1:3) = [   0.0000d0,   0.0000d0, -33.8707d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  4; allocate(geom.face( 1).poi( 4)); geom.face( 1).poi(1: 4) = [ 1,  2,  5,  3 ]
    geom.face( 2).n_poi =  4; allocate(geom.face( 2).poi( 4)); geom.face( 2).poi(1: 4) = [ 1,  3,  7,  4 ]
    geom.face( 3).n_poi =  4; allocate(geom.face( 3).poi( 4)); geom.face( 3).poi(1: 4) = [ 2,  6,  9,  5 ]
    geom.face( 4).n_poi =  4; allocate(geom.face( 4).poi( 4)); geom.face( 4).poi(1: 4) = [ 4,  7, 11,  8 ]
    geom.face( 5).n_poi =  4; allocate(geom.face( 5).poi( 4)); geom.face( 5).poi(1: 4) = [ 6, 10, 13,  9 ]
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi( 4)); geom.face( 6).poi(1: 4) = [ 8, 11, 15, 12 ]
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi( 4)); geom.face( 7).poi(1: 4) = [10, 14, 17, 13 ]
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi( 4)); geom.face( 8).poi(1: 4) = [12, 15, 19, 16 ]
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi( 4)); geom.face( 9).poi(1: 4) = [14, 18, 20, 17 ]
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi( 4)); geom.face(10).poi(1: 4) = [16, 19, 20, 18 ]
    geom.face(11).n_poi = 10; allocate(geom.face(11).poi(10)); geom.face(11).poi(1:10) = [ 1,  4,  8, 12, 16, 18, 14, 10,  6, 2 ]
    geom.face(12).n_poi = 10; allocate(geom.face(12).poi(10)); geom.face(12).poi(1:10) = [ 3,  5,  9, 13, 17, 20, 19, 15, 11, 7 ]
end subroutine Exam_Prism_Decagonal

! ---------------------------------------------------------------------------------------

! Example of square antiprism
subroutine Exam_Antiprism_Square(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Square"
    prob.name_file = "Antiprism_Square"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 8
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  16.4533d0 ]
    geom.iniP( 2).pos(1:3) = [  15.8821d0,  -0.0000d0,   4.2977d0 ]
    geom.iniP( 3).pos(1:3) = [   3.2893d0,  15.5378d0,   4.2977d0 ]
    geom.iniP( 4).pos(1:3) = [ -14.5196d0,   6.4360d0,   4.2977d0 ]
    geom.iniP( 5).pos(1:3) = [  -9.3035d0, -12.8719d0,   4.2977d0 ]
    geom.iniP( 6).pos(1:3) = [   6.5786d0, -12.8719d0,  -7.8580d0 ]
    geom.iniP( 7).pos(1:3) = [   7.9411d0,   6.4360d0, -12.8930d0 ]
    geom.iniP( 8).pos(1:3) = [  -9.8679d0,  -2.6659d0, -12.8930d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1, 2, 3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1, 3, 4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1, 4, 5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2, 6, 7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2, 7, 3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4, 8, 5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5, 8, 6 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 8, 7 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 1, 5, 6, 2 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 3, 7, 8, 4 ]
end subroutine Exam_Antiprism_Square

! ---------------------------------------------------------------------------------------

! Example of pentagonal antiprism
subroutine Exam_Antiprism_Pentagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Pentagonal"
    prob.name_file = "18_Antiprism_Pentagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xyz2", 1.0d0, 1.15d0, 0.0d0, 0.5d0)

    ! The number of points and faces
    geom.n_iniP = 10
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,   0.0000d0,  19.0212d0 ]
    geom.iniP( 2).pos(1:3) = [  17.0130d0,   0.0000d0,   8.5065d0 ]
    geom.iniP( 3).pos(1:3) = [   5.2573d0,  16.1804d0,   8.5065d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.7638d0,  10.0000d0,   8.5065d0 ]
    geom.iniP( 5).pos(1:3) = [ -13.7638d0, -10.0000d0,   8.5065d0 ]
    geom.iniP( 6).pos(1:3) = [  13.7638d0, -10.0000d0,  -8.5065d0 ]
    geom.iniP( 7).pos(1:3) = [  13.7638d0,  10.0000d0,  -8.5065d0 ]
    geom.iniP( 8).pos(1:3) = [ -17.0130d0,   0.0000d0,  -8.5065d0 ]
    geom.iniP( 9).pos(1:3) = [  -5.2573d0, -16.1804d0,  -8.5065d0 ]
    geom.iniP(10).pos(1:3) = [  -0.0000d0,   0.0000d0, -19.0212d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6,  9, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 10,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 10,  9 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 1,  5,  9, 6, 2 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 3,  7, 10, 8, 4 ]
end subroutine Exam_Antiprism_Pentagonal

! ---------------------------------------------------------------------------------------

! Example of hexagonal antiprism
subroutine Exam_Antiprism_Hexagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Hexagonal"
    prob.name_file = "19_Antiprism_Hexagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xyz2", 1.0d0, 1.1d0, 1.5d0, 1.0d0)

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  21.7533d0 ]
    geom.iniP( 2).pos(1:3) = [  17.7615d0,  -0.0000d0,  12.5593d0 ]
    geom.iniP( 3).pos(1:3) = [   6.5012d0,  16.5290d0,  12.5593d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.0023d0,  12.1000d0,  12.5593d0 ]
    geom.iniP( 5).pos(1:3) = [ -16.0195d0,  -7.6711d0,  12.5593d0 ]
    geom.iniP( 6).pos(1:3) = [  19.5035d0,  -7.6711d0,  -5.8288d0 ]
    geom.iniP( 7).pos(1:3) = [  17.7615d0,  12.1000d0,  -3.3652d0 ]
    geom.iniP( 8).pos(1:3) = [ -21.2456d0,   3.2422d0,  -3.3652d0 ]
    geom.iniP( 9).pos(1:3) = [ -14.2775d0, -15.3422d0,  -5.8288d0 ]
    geom.iniP(10).pos(1:3) = [   3.4840d0, -15.3422d0, -15.0228d0 ]
    geom.iniP(11).pos(1:3) = [   9.5184d0,   3.2422d0, -19.2898d0 ]
    geom.iniP(12).pos(1:3) = [  -9.9851d0,  -1.1867d0, -19.2898d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 12, 11 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [ 1,  5,  9, 10, 6, 2 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [ 3,  7, 11, 12, 8, 4 ]
end subroutine Exam_Antiprism_Hexagonal

! ---------------------------------------------------------------------------------------

! Example of heptagonal antiprism
subroutine Exam_Antiprism_Heptagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Heptagonal"
    prob.name_file = "20_Antiprism_Heptagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xyz2", 1.0d0, 0.98d0, 0.0d0, 0.5d0)

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -0.0000d0,  -0.0000d0,  24.5946d0 ]
    geom.iniP( 2).pos(1:3) = [  18.2722d0,  -0.0000d0,  16.4627d0 ]
    geom.iniP( 3).pos(1:3) = [   7.3266d0,  16.7390d0,  16.4627d0 ]
    geom.iniP( 4).pos(1:3) = [ -12.3967d0,  13.4237d0,  16.4627d0 ]
    geom.iniP( 5).pos(1:3) = [ -17.2680d0,  -5.9741d0,  16.4627d0 ]
    geom.iniP( 6).pos(1:3) = [  23.7892d0,  -5.9741d0,  -1.8095d0 ]
    geom.iniP( 7).pos(1:3) = [  20.5286d0,  13.4237d0,   1.8095d0 ]
    geom.iniP( 8).pos(1:3) = [ -23.7892d0,   5.9741d0,   1.8095d0 ]
    geom.iniP( 9).pos(1:3) = [ -20.5286d0, -13.4237d0,  -1.8095d0 ]
    geom.iniP(10).pos(1:3) = [  12.3967d0, -13.4237d0, -16.4627d0 ]
    geom.iniP(11).pos(1:3) = [  17.2680d0,   5.9741d0, -16.4627d0 ]
    geom.iniP(12).pos(1:3) = [ -18.2722d0,  -0.0000d0, -16.4627d0 ]
    geom.iniP(13).pos(1:3) = [  -7.3266d0, -16.7390d0, -16.4627d0 ]
    geom.iniP(14).pos(1:3) = [  -0.0000d0,  -0.0000d0, -24.5946d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 13, 14 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 14, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 14, 13 ]
    geom.face(15).n_poi = 7; allocate(geom.face(15).poi(7)); geom.face(15).poi(1:7) = [ 1,  5,  9, 13, 10, 6, 2 ]
    geom.face(16).n_poi = 7; allocate(geom.face(16).poi(7)); geom.face(16).poi(1:7) = [ 3,  7, 11, 14, 12, 8, 4 ]
end subroutine Exam_Antiprism_Heptagonal

! ---------------------------------------------------------------------------------------

! Example of octagonal antiprism
subroutine Exam_Antiprism_Octagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Octagonal"
    prob.name_file = "Antiprism_Octagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,   0.0000d0,  27.5111d0 ]
    geom.iniP( 2).pos(1:3) = [  18.6320d0,   0.0000d0,  20.2411d0 ]
    geom.iniP( 3).pos(1:3) = [   7.8977d0,  16.8753d0,  20.2411d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.9366d0,  14.3062d0,  20.2411d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.0171d0,  -4.7471d0,  20.2411d0 ]
    geom.iniP( 6).pos(1:3) = [  26.9645d0,  -4.7471d0,   2.6902d0 ]
    geom.iniP( 7).pos(1:3) = [  22.4908d0,  14.3062d0,   6.8083d0 ]
    geom.iniP( 8).pos(1:3) = [ -25.3935d0,   8.1038d0,   6.8083d0 ]
    geom.iniP( 9).pos(1:3) = [ -24.8651d0, -11.4605d0,   2.6902d0 ]
    geom.iniP(10).pos(1:3) = [  20.1164d0, -11.4605d0, -14.8607d0 ]
    geom.iniP(11).pos(1:3) = [  23.2942d0,   8.1038d0, -12.1887d0 ]
    geom.iniP(12).pos(1:3) = [ -24.5901d0,   1.9014d0, -12.1887d0 ]
    geom.iniP(13).pos(1:3) = [ -16.5326d0, -16.2076d0, -14.8607d0 ]
    geom.iniP(14).pos(1:3) = [   2.0994d0, -16.2076d0, -22.1306d0 ]
    geom.iniP(15).pos(1:3) = [   9.8372d0,   1.9014d0, -25.6216d0 ]
    geom.iniP(16).pos(1:3) = [  -9.9971d0,  -0.6677d0, -25.6216d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [13, 16, 14 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [14, 16, 15 ]
    geom.face(17).n_poi = 8; allocate(geom.face(17).poi(8)); geom.face(17).poi(1:8) = [ 1,  5,  9, 13, 14, 10, 6, 2 ]
    geom.face(18).n_poi = 8; allocate(geom.face(18).poi(8)); geom.face(18).poi(1:8) = [ 3,  7, 11, 15, 16, 12, 8, 4 ]
end subroutine Exam_Antiprism_Octagonal

! ---------------------------------------------------------------------------------------

! Example of enneagonal antiprism
subroutine Exam_Antiprism_Enneagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Enneagonal"
    prob.name_file = "Antiprism_Enneagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  30.4809d0 ]
    geom.iniP( 2).pos(1:3) = [  18.8930d0,  -0.0000d0,  23.9195d0 ]
    geom.iniP( 3).pos(1:3) = [   8.3072d0,  16.9688d0,  23.9195d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.5879d0,  14.9221d0,  23.9195d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.4974d0,  -3.8465d0,  23.9195d0 ]
    geom.iniP( 6).pos(1:3) = [  29.3415d0,  -3.8465d0,   7.3052d0 ]
    geom.iniP( 7).pos(1:3) = [  23.9195d0,  14.9221d0,  11.5879d0 ]
    geom.iniP( 8).pos(1:3) = [ -26.4565d0,   9.7397d0,  11.5879d0 ]
    geom.iniP( 9).pos(1:3) = [ -27.9439d0,  -9.7397d0,   7.3052d0 ]
    geom.iniP(10).pos(1:3) = [  26.4565d0,  -9.7397d0, -11.5879d0 ]
    geom.iniP(11).pos(1:3) = [  27.9439d0,   9.7397d0,  -7.3052d0 ]
    geom.iniP(12).pos(1:3) = [ -29.3415d0,   3.8465d0,  -7.3052d0 ]
    geom.iniP(13).pos(1:3) = [ -23.9195d0, -14.9221d0, -11.5879d0 ]
    geom.iniP(14).pos(1:3) = [  11.5879d0, -14.9221d0, -23.9195d0 ]
    geom.iniP(15).pos(1:3) = [  18.4974d0,   3.8465d0, -23.9195d0 ]
    geom.iniP(16).pos(1:3) = [ -18.8930d0,  -0.0000d0, -23.9195d0 ]
    geom.iniP(17).pos(1:3) = [  -8.3072d0, -16.9688d0, -23.9195d0 ]
    geom.iniP(18).pos(1:3) = [   0.0000d0,  -0.0000d0, -30.4809d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [13, 16, 17 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [14, 17, 18 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [14, 18, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [16, 18, 17 ]
    geom.face(19).n_poi = 9; allocate(geom.face(19).poi(9)); geom.face(19).poi(1:9) = [ 1,  5,  9, 13, 17, 14, 10, 6, 2 ]
    geom.face(20).n_poi = 9; allocate(geom.face(20).poi(9)); geom.face(20).poi(1:9) = [ 3,  7, 11, 15, 18, 16, 12, 8, 4 ]
end subroutine Exam_Antiprism_Enneagonal

! ---------------------------------------------------------------------------------------

! Example of decagonal antiprism
subroutine Exam_Antiprism_Decagonal(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Antiprism Decagonal"
    prob.name_file = "Antiprism_Decagonal"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [77, 175, 74], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 22

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.0000d0,  -0.0000d0,  33.4901d0 ]
    geom.iniP( 2).pos(1:3) = [  19.0876d0,  -0.0000d0,  27.5182d0 ]
    geom.iniP( 3).pos(1:3) = [   8.6096d0,  17.0356d0,  27.5182d0 ]
    geom.iniP( 4).pos(1:3) = [ -11.3208d0,  15.3681d0,  27.5182d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.8222d0,  -3.1719d0,  27.5182d0 ]
    geom.iniP( 6).pos(1:3) = [  31.1498d0,  -3.1719d0,  11.8835d0 ]
    geom.iniP( 7).pos(1:3) = [  24.9860d0,  15.3681d0,  16.1589d0 ]
    geom.iniP( 8).pos(1:3) = [ -27.1924d0,  11.0023d0,  16.1589d0 ]
    geom.iniP( 9).pos(1:3) = [ -30.1896d0,  -8.3041d0,  11.8835d0 ]
    geom.iniP(10).pos(1:3) = [  31.5792d0,  -8.3041d0,  -7.4420d0 ]
    geom.iniP(11).pos(1:3) = [  31.5533d0,  11.0023d0,  -2.2207d0 ]
    geom.iniP(12).pos(1:3) = [ -32.9429d0,   5.6060d0,  -2.2207d0 ]
    geom.iniP(13).pos(1:3) = [ -29.7602d0, -13.4364d0,  -7.4420d0 ]
    geom.iniP(14).pos(1:3) = [  20.2118d0, -13.4364d0, -23.0767d0 ]
    geom.iniP(15).pos(1:3) = [  25.8028d0,   5.6060d0, -20.6004d0 ]
    geom.iniP(16).pos(1:3) = [ -26.3756d0,   1.2402d0, -20.6004d0 ]
    geom.iniP(17).pos(1:3) = [ -17.6980d0, -16.6083d0, -23.0767d0 ]
    geom.iniP(18).pos(1:3) = [   1.3896d0, -16.6083d0, -29.0487d0 ]
    geom.iniP(19).pos(1:3) = [   9.9312d0,   1.2402d0, -31.9597d0 ]
    geom.iniP(20).pos(1:3) = [  -9.9992d0,  -0.4274d0, -31.9597d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [ 1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 2,  6,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 2,  7,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [ 4,  8,  5 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [ 5,  8,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [ 6, 10, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [ 6, 11,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [ 8, 12,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi( 3)); geom.face(11).poi(1: 3) = [ 9, 12, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi( 3)); geom.face(12).poi(1: 3) = [10, 14, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi( 3)); geom.face(13).poi(1: 3) = [10, 15, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi( 3)); geom.face(14).poi(1: 3) = [12, 16, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi( 3)); geom.face(15).poi(1: 3) = [13, 16, 17 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi( 3)); geom.face(16).poi(1: 3) = [14, 18, 19 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi( 3)); geom.face(17).poi(1: 3) = [14, 19, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi( 3)); geom.face(18).poi(1: 3) = [16, 20, 17 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi( 3)); geom.face(19).poi(1: 3) = [17, 20, 18 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi( 3)); geom.face(20).poi(1: 3) = [18, 20, 19 ]
    geom.face(21).n_poi =10; allocate(geom.face(21).poi(10)); geom.face(21).poi(1:10) = [ 1,  5,  9, 13, 17, 18, 14, 10, 6, 2 ]
    geom.face(22).n_poi =10; allocate(geom.face(22).poi(10)); geom.face(22).poi(1:10) = [ 3,  7, 11, 15, 19, 20, 16, 12, 8, 4 ]
end subroutine Exam_Antiprism_Decagonal

! ---------------------------------------------------------------------------------------

end module Exam_Prism