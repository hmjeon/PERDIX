!
! ---------------------------------------------------------------------------------------
!
!                               Biscribed Chiral Solids
!
!                                                                    Updated : 2017/05/03
!
! Comments: Biscribed chiral solids.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Exam_Chiral

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Chiral_Biscribed_Propello_Tetrahedron
    public Exam_Chiral_Biscribed_Propello_Cube
    public Exam_Chiral_Biscribed_Propello_Octahedron
    public Exam_Chiral_Biscribed_Snub_Cube
    public Exam_Chiral_Biscribed_Pentagonal_Icositetrahedron

    public Exam_Chiral_Asym_Tetrahedron        ! V=4,  E=6,  F=4
    public Exam_Chiral_Asym_Cube               ! V=8,  E=12, F=6
    public Exam_Chiral_Asym_Octahedron         ! V=6,  E=12, F=8
    public Exam_Chiral_Asym_Dodecahedron       ! V=20, E=30, F=12
    public Exam_Chiral_Asym_Icosahedron        ! V=12, E=30, F=20

contains

! ---------------------------------------------------------------------------------------

! Example of chiral biscribed propello tetrahedron
subroutine Exam_Chiral_Biscribed_Propello_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2, c3
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Bi_Pro_Tetrahedron"
    prob.name_file = "Chiral_Bi_Pro_Tetrahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    c0 = 0.041226573367477290238151003637d0
    c1 = 0.447594291856692697052238850001d0
    c2 = 0.577350269189625764509148780502d0
    c3 = 0.893285911422363030309232997113d0

    geom.iniP( 1).pos(1:3) = [  c1,  c0,  c3 ]
    geom.iniP( 2).pos(1:3) = [  c1, -c0, -c3 ]
    geom.iniP( 3).pos(1:3) = [ -c1, -c0,  c3 ]
    geom.iniP( 4).pos(1:3) = [ -c1,  c0, -c3 ]
    geom.iniP( 5).pos(1:3) = [  c3,  c1,  c0 ]
    geom.iniP( 6).pos(1:3) = [  c3, -c1, -c0 ]
    geom.iniP( 7).pos(1:3) = [ -c3, -c1,  c0 ]
    geom.iniP( 8).pos(1:3) = [ -c3,  c1, -c0 ]
    geom.iniP( 9).pos(1:3) = [  c0,  c3,  c1 ]
    geom.iniP(10).pos(1:3) = [  c0, -c3, -c1 ]
    geom.iniP(11).pos(1:3) = [ -c0, -c3,  c1 ]
    geom.iniP(12).pos(1:3) = [ -c0,  c3, -c1 ]
    geom.iniP(13).pos(1:3) = [  c2, -c2,  c2 ]
    geom.iniP(14).pos(1:3) = [  c2,  c2, -c2 ]
    geom.iniP(15).pos(1:3) = [ -c2,  c2,  c2 ]
    geom.iniP(16).pos(1:3) = [ -c2, -c2, -c2 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  13,  1,  3, 11 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  13, 11, 10,  6 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  13,  6,  5,  1 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  14,  2,  4, 12 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  14, 12,  9,  5 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  14,  5,  6,  2 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  15,  3,  1,  9 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  15,  9, 12,  8 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  15,  8,  7,  3 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  16,  4,  2, 10 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  16, 10, 11,  7 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  16,  7,  8,  4 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [   1,  5,  9 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [   2,  6, 10 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [   3,  7, 11 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [   4,  8, 12 ]
end subroutine Exam_Chiral_Biscribed_Propello_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of chiral biscribed propello cube
subroutine Exam_Chiral_Biscribed_Propello_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Bi_Pro_Cube"
    prob.name_file = "Chiral_Bi_Pro_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 32
    geom.n_face = 30

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    c0 = 0.268318503889892044924746743411d0
    c1 = 0.437902382065032587008255021493d0
    c2 = 0.649038600739057169667625426214d0

    geom.iniP( 1).pos(1:3) = [     c1,     c0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [     c1,    -c0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [    -c1,    -c0,  1.0d0 ]
    geom.iniP( 4).pos(1:3) = [    -c1,     c0, -1.0d0 ]
    geom.iniP( 5).pos(1:3) = [  1.0d0,     c1,     c0 ]
    geom.iniP( 6).pos(1:3) = [  1.0d0,    -c1,    -c0 ]
    geom.iniP( 7).pos(1:3) = [ -1.0d0,    -c1,     c0 ]
    geom.iniP( 8).pos(1:3) = [ -1.0d0,     c1,    -c0 ]
    geom.iniP( 9).pos(1:3) = [     c0,  1.0d0,     c1 ]
    geom.iniP(10).pos(1:3) = [     c0, -1.0d0,    -c1 ]
    geom.iniP(11).pos(1:3) = [    -c0, -1.0d0,     c1 ]
    geom.iniP(12).pos(1:3) = [    -c0,  1.0d0,    -c1 ]
    geom.iniP(13).pos(1:3) = [     c0,    -c1,  1.0d0 ]
    geom.iniP(14).pos(1:3) = [     c0,     c1, -1.0d0 ]
    geom.iniP(15).pos(1:3) = [    -c0,     c1,  1.0d0 ]
    geom.iniP(16).pos(1:3) = [    -c0,    -c1, -1.0d0 ]
    geom.iniP(17).pos(1:3) = [  1.0d0,    -c0,     c1 ]
    geom.iniP(18).pos(1:3) = [  1.0d0,     c0,    -c1 ]
    geom.iniP(19).pos(1:3) = [ -1.0d0,     c0,     c1 ]
    geom.iniP(20).pos(1:3) = [ -1.0d0,    -c0,    -c1 ]
    geom.iniP(21).pos(1:3) = [     c1, -1.0d0,     c0 ]
    geom.iniP(22).pos(1:3) = [     c1,  1.0d0,    -c0 ]
    geom.iniP(23).pos(1:3) = [    -c1,  1.0d0,     c0 ]
    geom.iniP(24).pos(1:3) = [    -c1, -1.0d0,    -c0 ]
    geom.iniP(25).pos(1:3) = [     c2,     c2,     c2 ]
    geom.iniP(26).pos(1:3) = [     c2,     c2,    -c2 ]
    geom.iniP(27).pos(1:3) = [     c2,    -c2,     c2 ]
    geom.iniP(28).pos(1:3) = [     c2,    -c2,    -c2 ]
    geom.iniP(29).pos(1:3) = [    -c2,     c2,     c2 ]
    geom.iniP(30).pos(1:3) = [    -c2,     c2,    -c2 ]
    geom.iniP(31).pos(1:3) = [    -c2,    -c2,     c2 ]
    geom.iniP(32).pos(1:3) = [    -c2,    -c2,    -c2 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  2, 12,  0, 14 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3, 13,  1, 15 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  4, 16,  5, 17 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  7, 19,  6, 18 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  8, 21, 11, 22 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  9, 20, 10, 23 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 24,  0, 16,  4 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 24,  4, 21,  8 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 24,  8, 14,  0 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 25, 13, 11, 21 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 25, 21,  4, 17 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 25, 17,  1, 13 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 26, 12, 10, 20 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 26, 20,  5, 16 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 26, 16,  0, 12 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 27,  1, 17,  5 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 27,  5, 20,  9 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 27,  9, 15,  1 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 28, 14,  8, 22 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 28, 22,  7, 18 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 28, 18,  2, 14 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 29,  3, 19,  7 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 29,  7, 22, 11 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 29, 11, 13,  3 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 30,  2, 18,  6 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 30,  6, 23, 10 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 30, 10, 12,  2 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 31, 15,  9, 23 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 31, 23,  6, 19 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 31, 19,  3, 15 ]
end subroutine Exam_Chiral_Biscribed_Propello_Cube

! ---------------------------------------------------------------------------------------

! Example of chiral biscribed propello octahedron
subroutine Exam_Chiral_Biscribed_Propello_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Bi_Pro_Octahedron"
    prob.name_file = "Chiral_Bi_Pro_Octahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 30
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.149825404557484083971488566549d0
    c1 = 0.623938024555993280267936683079d0
    c2 = 0.766976981181541688825873625866d0

    geom.iniP( 1).pos(1:3) = [  0.0d0,  0.0d0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [  0.0d0,  0.0d0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [  1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 4).pos(1:3) = [ -1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  0.0d0,  1.0d0,  0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  0.0d0, -1.0d0,  0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   c1,  -c0,   c2 ]
    geom.iniP( 8).pos(1:3) = [   c1,   c0,  -c2 ]
    geom.iniP( 9).pos(1:3) = [  -c1,   c0,   c2 ]
    geom.iniP(10).pos(1:3) = [  -c1,  -c0,  -c2 ]
    geom.iniP(11).pos(1:3) = [   c2,  -c1,   c0 ]
    geom.iniP(12).pos(1:3) = [   c2,   c1,  -c0 ]
    geom.iniP(13).pos(1:3) = [  -c2,   c1,   c0 ]
    geom.iniP(14).pos(1:3) = [  -c2,  -c1,  -c0 ]
    geom.iniP(15).pos(1:3) = [   c0,  -c2,   c1 ]
    geom.iniP(16).pos(1:3) = [   c0,   c2,  -c1 ]
    geom.iniP(17).pos(1:3) = [  -c0,   c2,   c1 ]
    geom.iniP(18).pos(1:3) = [  -c0,  -c2,  -c1 ]
    geom.iniP(19).pos(1:3) = [   c0,   c1,   c2 ]
    geom.iniP(20).pos(1:3) = [   c0,  -c1,  -c2 ]
    geom.iniP(21).pos(1:3) = [  -c0,  -c1,   c2 ]
    geom.iniP(22).pos(1:3) = [  -c0,   c1,  -c2 ]
    geom.iniP(23).pos(1:3) = [   c2,   c0,   c1 ]
    geom.iniP(24).pos(1:3) = [   c2,  -c0,  -c1 ]
    geom.iniP(25).pos(1:3) = [  -c2,  -c0,   c1 ]
    geom.iniP(26).pos(1:3) = [  -c2,   c0,  -c1 ]
    geom.iniP(27).pos(1:3) = [   c1,   c2,   c0 ]
    geom.iniP(28).pos(1:3) = [   c1,  -c2,  -c0 ]
    geom.iniP(29).pos(1:3) = [  -c1,  -c2,   c0 ]
    geom.iniP(30).pos(1:3) = [  -c1,   c2,  -c0 ]

    ! Set connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  0,  6, 22, 18 ] + 1
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  0, 18, 16,  8 ] + 1
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  0,  8, 24, 20 ] + 1
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  0, 20, 14,  6 ] + 1
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  1,  7, 23, 19 ] + 1
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  1, 19, 17,  9 ] + 1
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  1,  9, 25, 21 ] + 1
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  1, 21, 15,  7 ] + 1
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  2, 10, 27, 23 ] + 1
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  2, 23,  7, 11 ] + 1
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  2, 11, 26, 22 ] + 1
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  2, 22,  6, 10 ] + 1
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  3, 12, 29, 25 ] + 1
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  3, 25,  9, 13 ] + 1
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  3, 13, 28, 24 ] + 1
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  3, 24,  8, 12 ] + 1
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  4, 15, 21, 29 ] + 1
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  4, 29, 12, 16 ] + 1
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [  4, 16, 18, 26 ] + 1
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [  4, 26, 11, 15 ] + 1
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [  5, 14, 20, 28 ] + 1
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [  5, 28, 13, 17 ] + 1
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [  5, 17, 19, 27 ] + 1
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [  5, 27, 10, 14 ] + 1
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 14, 10,  6 ] + 1
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 15, 11,  7 ] + 1
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 16, 12,  8 ] + 1
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 17, 13,  9 ] + 1
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 18, 22, 26 ] + 1
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 19, 23, 27 ] + 1
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 20, 24, 28 ] + 1
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 21, 25, 29 ] + 1
end subroutine Exam_Chiral_Biscribed_Propello_Octahedron

! ---------------------------------------------------------------------------------------

! Example of chiral biscribed snub cube
subroutine Exam_Chiral_Biscribed_Snub_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Bi_Snub_Cube"
    prob.name_file = "Chiral_Bi_Snub_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 24
    geom.n_face = 38

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.33775397381375235213753224516503d0
    c1 = 0.62122641055658531169250095449000d0
    c2 = 1.14261350892596209347948408672000d0

    geom.iniP( 1).pos(1:3) = [  c1,  c0,  c2 ]
    geom.iniP( 2).pos(1:3) = [  c1, -c0, -c2 ]
    geom.iniP( 3).pos(1:3) = [ -c1, -c0,  c2 ]
    geom.iniP( 4).pos(1:3) = [ -c1,  c0, -c2 ]
    geom.iniP( 5).pos(1:3) = [  c2,  c1,  c0 ]
    geom.iniP( 6).pos(1:3) = [  c2, -c1, -c0 ]
    geom.iniP( 7).pos(1:3) = [ -c2, -c1,  c0 ]
    geom.iniP( 8).pos(1:3) = [ -c2,  c1, -c0 ]
    geom.iniP( 9).pos(1:3) = [  c0,  c2,  c1 ]
    geom.iniP(10).pos(1:3) = [  c0, -c2, -c1 ]
    geom.iniP(11).pos(1:3) = [ -c0, -c2,  c1 ]
    geom.iniP(12).pos(1:3) = [ -c0,  c2, -c1 ]
    geom.iniP(13).pos(1:3) = [  c0, -c1,  c2 ]
    geom.iniP(14).pos(1:3) = [  c0,  c1, -c2 ]
    geom.iniP(15).pos(1:3) = [ -c0,  c1,  c2 ]
    geom.iniP(16).pos(1:3) = [ -c0, -c1, -c2 ]
    geom.iniP(17).pos(1:3) = [  c2, -c0,  c1 ]
    geom.iniP(18).pos(1:3) = [  c2,  c0, -c1 ]
    geom.iniP(19).pos(1:3) = [ -c2,  c0,  c1 ]
    geom.iniP(20).pos(1:3) = [ -c2, -c0, -c1 ]
    geom.iniP(21).pos(1:3) = [  c1, -c2,  c0 ]
    geom.iniP(22).pos(1:3) = [  c1,  c2, -c0 ]
    geom.iniP(23).pos(1:3) = [ -c1,  c2,  c0 ]
    geom.iniP(24).pos(1:3) = [ -c1, -c2, -c0 ]

    ! Set connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  2, 12,  0, 14 ] + 1
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3, 13,  1, 15 ] + 1
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  4, 16,  5, 17 ] + 1
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  7, 19,  6, 18 ] + 1
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  8, 21, 11, 22 ] + 1
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  9, 20, 10, 23 ] + 1
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  0,  8, 14 ] + 1
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  1,  9, 15 ] + 1
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  2, 10, 12 ] + 1
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  3, 11, 13 ] + 1
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4,  0, 16 ] + 1
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  5,  1, 17 ] + 1
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  6,  2, 18 ] + 1
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  7,  3, 19 ] + 1
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  8,  4, 21 ] + 1
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  9,  5, 20 ] + 1
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 10,  6, 23 ] + 1
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 11,  7, 22 ] + 1
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 12, 16,  0 ] + 1
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 13, 17,  1 ] + 1
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 14, 18,  2 ] + 1
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 15, 19,  3 ] + 1
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 16, 20,  5 ] + 1
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 17, 21,  4 ] + 1
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 18, 22,  7 ] + 1
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 19, 23,  6 ] + 1
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 20, 12, 10 ] + 1
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 21, 13, 11 ] + 1
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 22, 14,  8 ] + 1
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 23, 15,  9 ] + 1
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  8,  0,  4 ] + 1
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [  9,  1,  5 ] + 1
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 10,  2,  6 ] + 1
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 11,  3,  7 ] + 1
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 12, 20, 16 ] + 1
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 13, 21, 17 ] + 1
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 14, 22, 18 ] + 1
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 15, 23, 19 ] + 1
end subroutine Exam_Chiral_Biscribed_Snub_Cube

! ---------------------------------------------------------------------------------------

! Example of chiral biscribed pentagonal icositetrahedron
subroutine Exam_Chiral_Biscribed_Pentagonal_Icositetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: c0, c1, c2, c3
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Bi_Penta_Icositetrahedron"
    prob.name_file = "Chiral_Bi_Penta_Icositetrahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 38
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    c0 = 0.107540025570073586711065317103d0
    c1 = 0.577350269189625764509148780502d0
    c2 = 0.643438410432318999656172634091d0
    c3 = 0.757906428842451843370418995303d0

    geom.iniP( 1).pos(1:3) = [  0.0d0,  0.0d0,  1.0d0 ]
    geom.iniP( 2).pos(1:3) = [  0.0d0,  0.0d0, -1.0d0 ]
    geom.iniP( 3).pos(1:3) = [  1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 4).pos(1:3) = [ -1.0d0,  0.0d0,  0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  0.0d0,  1.0d0,  0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  0.0d0, -1.0d0,  0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   c2,  -c0,   c3 ]
    geom.iniP( 8).pos(1:3) = [   c2,   c0,  -c3 ]
    geom.iniP( 9).pos(1:3) = [  -c2,   c0,   c3 ]
    geom.iniP(10).pos(1:3) = [  -c2,  -c0,  -c3 ]
    geom.iniP(11).pos(1:3) = [   c3,  -c2,   c0 ]
    geom.iniP(12).pos(1:3) = [   c3,   c2,  -c0 ]
    geom.iniP(13).pos(1:3) = [  -c3,   c2,   c0 ]
    geom.iniP(14).pos(1:3) = [  -c3,  -c2,  -c0 ]
    geom.iniP(15).pos(1:3) = [   c0,  -c3,   c2 ]
    geom.iniP(16).pos(1:3) = [   c0,   c3,  -c2 ]
    geom.iniP(17).pos(1:3) = [  -c0,   c3,   c2 ]
    geom.iniP(18).pos(1:3) = [  -c0,  -c3,  -c2 ]
    geom.iniP(19).pos(1:3) = [   c0,   c2,   c3 ]
    geom.iniP(20).pos(1:3) = [   c0,  -c2,  -c3 ]
    geom.iniP(21).pos(1:3) = [  -c0,  -c2,   c3 ]
    geom.iniP(22).pos(1:3) = [  -c0,   c2,  -c3 ]
    geom.iniP(23).pos(1:3) = [   c3,   c0,   c2 ]
    geom.iniP(24).pos(1:3) = [   c3,  -c0,  -c2 ]
    geom.iniP(25).pos(1:3) = [  -c3,  -c0,   c2 ]
    geom.iniP(26).pos(1:3) = [  -c3,   c0,  -c2 ]
    geom.iniP(27).pos(1:3) = [   c2,   c3,   c0 ]
    geom.iniP(28).pos(1:3) = [   c2,  -c3,  -c0 ]
    geom.iniP(29).pos(1:3) = [  -c2,  -c3,   c0 ]
    geom.iniP(30).pos(1:3) = [  -c2,   c3,  -c0 ]
    geom.iniP(31).pos(1:3) = [   c1,   c1,   c1 ]
    geom.iniP(32).pos(1:3) = [   c1,   c1,  -c1 ]
    geom.iniP(33).pos(1:3) = [   c1,  -c1,   c1 ]
    geom.iniP(34).pos(1:3) = [   c1,  -c1,  -c1 ]
    geom.iniP(35).pos(1:3) = [  -c1,   c1,   c1 ]
    geom.iniP(36).pos(1:3) = [  -c1,   c1,  -c1 ]
    geom.iniP(37).pos(1:3) = [  -c1,  -c1,   c1 ]
    geom.iniP(38).pos(1:3) = [  -c1,  -c1,  -c1 ]

    ! Set connectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [ 0,  6, 22, 30, 18 ] + 1
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [ 0, 18, 16, 34,  8 ] + 1
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [ 0,  8, 24, 36, 20 ] + 1
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [ 0, 20, 14, 32,  6 ] + 1
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [ 1,  7, 23, 33, 19 ] + 1
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [ 1, 19, 17, 37,  9 ] + 1
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [ 1,  9, 25, 35, 21 ] + 1
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [ 1, 21, 15, 31,  7 ] + 1
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [ 2, 10, 27, 33, 23 ] + 1
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 2, 23,  7, 31, 11 ] + 1
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 2, 11, 26, 30, 22 ] + 1
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 2, 22,  6, 32, 10 ] + 1
    geom.face(13).n_poi = 5; allocate(geom.face(13).poi(5)); geom.face(13).poi(1:5) = [ 3, 12, 29, 35, 25 ] + 1
    geom.face(14).n_poi = 5; allocate(geom.face(14).poi(5)); geom.face(14).poi(1:5) = [ 3, 25,  9, 37, 13 ] + 1
    geom.face(15).n_poi = 5; allocate(geom.face(15).poi(5)); geom.face(15).poi(1:5) = [ 3, 13, 28, 36, 24 ] + 1
    geom.face(16).n_poi = 5; allocate(geom.face(16).poi(5)); geom.face(16).poi(1:5) = [ 3, 24,  8, 34, 12 ] + 1
    geom.face(17).n_poi = 5; allocate(geom.face(17).poi(5)); geom.face(17).poi(1:5) = [ 4, 15, 21, 35, 29 ] + 1
    geom.face(18).n_poi = 5; allocate(geom.face(18).poi(5)); geom.face(18).poi(1:5) = [ 4, 29, 12, 34, 16 ] + 1
    geom.face(19).n_poi = 5; allocate(geom.face(19).poi(5)); geom.face(19).poi(1:5) = [ 4, 16, 18, 30, 26 ] + 1
    geom.face(20).n_poi = 5; allocate(geom.face(20).poi(5)); geom.face(20).poi(1:5) = [ 4, 26, 11, 31, 15 ] + 1
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [ 5, 14, 20, 36, 28 ] + 1
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 5, 28, 13, 37, 17 ] + 1
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 5, 17, 19, 33, 27 ] + 1
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 5, 27, 10, 32, 14 ] + 1
end subroutine Exam_Chiral_Biscribed_Pentagonal_Icositetrahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric object
subroutine Exam_Chiral_Asym_Object(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Object"
    prob.name_file = "Chiral_Asym_Object"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! Preset parameters
    if(para_preset == "on") then
        if(para_vertex_design == "flat") then
            para_junc_ang        = "min"    ! Junctional gap
            para_unpaired_scaf   = "off"    ! Unpaired scaffold nucleotides
            para_n_base_tn       = 7        ! The number of nucleotides

            ! Folding conditions
            para_const_edge_mesh = "on"     ! Constant edge length
        else if(para_vertex_design == "beveled") then
            para_const_edge_mesh = "on"     ! Constant edge length
            para_junc_ang        = "opt"    ! Junctional gap
            para_unpaired_scaf   = "on"     ! Unpaired scaffold nucleotides
            para_n_base_tn       = -1       ! The number of nucleotides
        end if
    end if

    ! Allocate point and face structure
    geom.n_iniP = 9
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.000000d0,  0.000000d0,  1.000000d0 ]
    geom.iniP(2).pos(1:3) = [  0.000000d0,  0.000000d0, -1.000000d0 ]
    geom.iniP(3).pos(1:3) = [  0.000000d0,  1.000000d0,  0.000000d0 ]
    geom.iniP(4).pos(1:3) = [  0.000000d0, -1.000000d0,  0.000000d0 ]
    geom.iniP(5).pos(1:3) = [  1.000000d0,  0.000000d0,  0.000000d0 ]
    geom.iniP(6).pos(1:3) = [ -1.000000d0,  0.000000d0,  0.000000d0 ]
    geom.iniP(7).pos(1:3) = [ -1.000000d0,  1.000000d0, -1.000000d0 ]
    geom.iniP(8).pos(1:3) = [  0.333333d0,  1.333333d0, -1.333333d0 ]
    geom.iniP(9).pos(1:3) = [  1.000000d0,  1.000000d0,  1.000000d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 1, 5, 9 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 1, 9, 3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 9, 5, 3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 1, 4, 5 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 6, 4 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 1, 3, 6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 2, 3, 5 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 2, 5, 4 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 2, 4, 6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 2, 6, 7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 3, 7, 6 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 3, 8, 7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 3, 2, 8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 2, 7, 8 ]
end subroutine Exam_Chiral_Asym_Object

! ---------------------------------------------------------------------------------------

! Example of asymmetric tetrahedron
subroutine Exam_Chiral_Asym_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Tetrahedron"
    prob.name_file = "Chiral_Asym_Tetrahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    geom.iniP(1).pos(1:3) = [ -5.9812d0,  12.5978d0, -4.5636d0 ]
    geom.iniP(2).pos(1:3) = [  0.8677d0,  -0.9451d0, 12.6382d0 ]
    geom.iniP(3).pos(1:3) = [ 11.6103d0,   1.0457d0, -4.1139d0 ]
    geom.iniP(4).pos(1:3) = [ -6.4967d0, -12.6984d0, -3.9607d0 ]

    ! Set connectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 2, 1, 4 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 2, 4, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 2, 3, 1 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 4, 1, 3 ]
end subroutine Exam_Chiral_Asym_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric cube
subroutine Exam_Chiral_Asym_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Cube"
    prob.name_file = "Chiral_Asym_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xyz")

    ! The number of points and faces
    geom.n_iniP = 8
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vectors
    geom.iniP(1).pos(1:3) = [ -10.8699d0,  -9.0403d0,  12.8253d0 ]
    geom.iniP(2).pos(1:3) = [  13.1783d0,   9.3523d0,   9.2211d0 ]
    geom.iniP(3).pos(1:3) = [  -9.4579d0,  13.0216d0,  10.4285d0 ]
    geom.iniP(4).pos(1:3) = [ -10.6283d0, -11.0596d0,  -8.5752d0 ]
    geom.iniP(5).pos(1:3) = [  10.3701d0,  11.8547d0, -10.4220d0 ]
    geom.iniP(6).pos(1:3) = [   8.6643d0, -13.3699d0,   9.2543d0 ]
    geom.iniP(7).pos(1:3) = [  11.1136d0, -12.4990d0, -11.8562d0 ]
    geom.iniP(8).pos(1:3) = [ -12.3703d0,  11.7402d0, -10.8758d0 ]

    ! Set connectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 6, 2, 3, 1 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 2, 5, 8, 3 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 8, 5, 7, 4 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 7, 6, 1, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 6, 7, 5, 2 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 1, 3, 8, 4 ]
end subroutine Exam_Chiral_Asym_Cube

! ---------------------------------------------------------------------------------------

! Example of asymmetric octahedron
subroutine Exam_Chiral_Asym_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Octahedron"
    prob.name_file = "Chiral_Asym_Octahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 6
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [   1.4393d0,    -1.8152d0,    18.6786d0 ]
    geom.iniP(2).pos(1:3) = [  -2.5337d0,   -15.4490d0,     0.6888d0 ]
    geom.iniP(3).pos(1:3) = [  15.9284d0,    -0.7617d0,    -1.4565d0 ]
    geom.iniP(4).pos(1:3) = [   1.1919d0,    -0.7591d0,   -19.0232d0 ]
    geom.iniP(5).pos(1:3) = [  -0.9576d0,    16.1710d0,     2.6237d0 ]
    geom.iniP(6).pos(1:3) = [ -15.0682d0,     2.6140d0,    -1.5113d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face( 1).poi(1:3) = [ 1, 3, 5 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face( 2).poi(1:3) = [ 1, 2, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face( 3).poi(1:3) = [ 1, 6, 2 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face( 4).poi(1:3) = [ 1, 5, 6 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face( 5).poi(1:3) = [ 4, 5, 3 ]
    geom.face(6).n_poi = 3; allocate(geom.face(6).poi(3)); geom.face( 6).poi(1:3) = [ 4, 3, 2 ]
    geom.face(7).n_poi = 3; allocate(geom.face(7).poi(3)); geom.face( 7).poi(1:3) = [ 4, 2, 6 ]
    geom.face(8).n_poi = 3; allocate(geom.face(8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 6, 5 ]
end subroutine Exam_Chiral_Asym_Octahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric dodecahedron
subroutine Exam_Chiral_Asym_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp
    
    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Dodecahedron"
    prob.name_file = "Chiral_Asym_Dodecahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xz")

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  17.1453d0,  19.3422d0, -21.2961d0 ]
    geom.iniP( 2).pos(1:3) = [  15.7068d0, -20.0261d0, -20.6164d0 ]
    geom.iniP( 3).pos(1:3) = [  15.2401d0,   4.7228d0,  28.5272d0 ]
    geom.iniP( 4).pos(1:3) = [  34.1515d0,  12.9296d0,   0.6245d0 ]
    geom.iniP( 5).pos(1:3) = [  34.1515d0, -13.1214d0,   0.6245d0 ]
    geom.iniP( 6).pos(1:3) = [ -11.7329d0,  -1.1709d0,  32.9814d0 ]
    geom.iniP( 7).pos(1:3) = [ -19.4756d0, -21.5229d0,  19.2166d0 ]
    geom.iniP( 8).pos(1:3) = [   1.5612d0, -34.1971d0,  13.6500d0 ]
    geom.iniP( 9).pos(1:3) = [  20.0907d0, -21.6051d0,  18.2737d0 ]
    geom.iniP(10).pos(1:3) = [   0.5849d0, -31.8230d0, -14.9453d0 ]
    geom.iniP(11).pos(1:3) = [ -18.9211d0, -20.8913d0, -21.1337d0 ]
    geom.iniP(12).pos(1:3) = [  13.0759d0,  -0.0959d0, -33.4768d0 ]
    geom.iniP(13).pos(1:3) = [ -19.2511d0,  23.0880d0, -16.9362d0 ]
    geom.iniP(14).pos(1:3) = [  20.8505d0,  22.2340d0,  16.9087d0 ]
    geom.iniP(15).pos(1:3) = [   4.0958d0,  32.8048d0,  -8.2061d0 ]
    geom.iniP(16).pos(1:3) = [ -20.1670d0,  15.2966d0,  22.9382d0 ]
    geom.iniP(17).pos(1:3) = [  -1.8085d0,  33.5410d0,  12.2026d0 ]
    geom.iniP(18).pos(1:3) = [ -34.0509d0,  12.9296d0,   0.6245d0 ]
    geom.iniP(19).pos(1:3) = [ -17.1963d0,   0.6866d0, -30.5855d0 ]
    geom.iniP(20).pos(1:3) = [ -34.0509d0, -13.1214d0,   0.6245d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [ 12,  1,  4,  5,  2 ]
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [ 14,  3,  9,  5,  4 ]
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [  8, 10,  2,  5,  9 ]
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [  3,  6,  7,  8,  9 ]
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [  8,  7, 20, 11, 10 ]
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [ 10, 11, 19, 12,  2 ]
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [ 12, 19, 13, 15,  1 ]
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [ 17, 14,  4,  1, 15 ]
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [  6,  3, 14, 17, 16 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [  6, 16, 18, 20,  7 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 15, 13, 18, 16, 17 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 13, 19, 11, 20, 18 ]
end subroutine Exam_Chiral_Asym_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric icosahedron
subroutine Exam_Chiral_Asym_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_prob = "Chiral_Asym_Icosahedron"
    prob.name_file = "Chiral_Asym_Icosahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&
        "_"//trim(para_cut_stap_method)

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy")

    ! The number of points and faces
    geom.n_iniP = 12
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  12.8559d0, -18.4078d0,   4.5927d0 ]
    geom.iniP( 2).pos(1:3) = [   2.6353d0, -18.2908d0, -18.0086d0 ]
    geom.iniP( 3).pos(1:3) = [  19.4301d0,  -1.3006d0, -15.5467d0 ]
    geom.iniP( 4).pos(1:3) = [  20.6532d0,   2.9721d0,  11.8687d0 ]
    geom.iniP( 5).pos(1:3) = [ -19.4518d0,   1.5734d0,  10.8306d0 ]
    geom.iniP( 6).pos(1:3) = [ -21.3999d0,   1.9363d0, -13.0532d0 ]
    geom.iniP( 7).pos(1:3) = [   0.3631d0, -12.8147d0,  19.8932d0 ]
    geom.iniP( 8).pos(1:3) = [   2.7443d0,   8.5393d0,  18.8170d0 ]
    geom.iniP( 9).pos(1:3) = [  -1.9827d0,  14.4621d0, -20.2903d0 ]
    geom.iniP(10).pos(1:3) = [ -14.9666d0, -19.2470d0,   2.2430d0 ]
    geom.iniP(11).pos(1:3) = [  11.3689d0,  20.2643d0,  -0.3120d0 ]
    geom.iniP(12).pos(1:3) = [ -12.2499d0,  20.3134d0,  -1.0344d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 4,  1,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 3, 11,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 3,  2,  9 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 3,  9, 11 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 3,  1,  2 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 7,  4,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 1,  4,  7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 11,  8 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 10,  5 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 5, 12,  6 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 7,  8,  5 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 5, 10,  7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 8, 12,  5 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 9,  2,  6 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 6, 12,  9 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 6,  2, 10 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 7, 10,  1 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 8, 11, 12 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 9, 12, 11 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 2,  1, 10 ]
end subroutine Exam_Chiral_Asym_Icosahedron

! ---------------------------------------------------------------------------------------

end module Exam_Chiral