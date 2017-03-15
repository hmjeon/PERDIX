!
! ---------------------------------------------------------------------------------------
!
!                                Module for Exam_Johnson
!
!                                                             First modified : 2017/03/13
!                                                             Last  modified : 2017/03/13
!
! Comments: This module contains the geometric definition of 2D open structures
!
! by Hyungmin Jun (Hyungminjun@outlook.com), MIT, Bathe Lab, 2017
!
! Copyright 2017. Massachusetts Institute of Technology. Rights Reserved.
! M.I.T. hereby makes following copyrightable material available to the
! public under GNU General Public License, version 2 (GPL-2.0). A copy of
! this license is available at https://opensource.org/licenses/GPL-2.0
!
! ---------------------------------------------------------------------------------------
!
! Reference: http://www.georgehart.com/virtual-polyhedra/johnson-index.html
! Reference: https://en.wikipedia.org/wiki/List_of_Johnson_solids
!
module Exam_Johnson

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Johnson_Square_Pyramid_J1                       ! V= 5, E= 8, F= 5
    public Exam_Johnson_Pentagonal_Pyramid_J2                   ! V= 6, E=10, F= 6
    public Exam_Johnson_Triangular_Cupola_J3                    ! V= 9, E=15, F= 8
    public Exam_Johnson_Square_Cupola_J4                        ! V=12, E=20, F=10
    public Exam_Johnson_Pentagonal_Cupola_J5                    ! V=15, E=25, F=12
    public Exam_Johnson_Pentagonal_Rotunda_J6                   ! V=20, E=35, F=17
    public Exam_Johnson_Elongated_Triangular_Cupola_J18         ! V=15, E=27, F=14
    public Exam_Johnson_Elongated_Square_Cupola_J19             ! V=20, E=36, F=18
    public Exam_Johnson_Elongated_Pentagonal_Cupola_J20         ! V=25, E=45, F=22
    public Exam_Johnson_Elongated_Pentagonal_Rotunda_J21        ! V=30, E=55, F=27
    public Exam_Johnson_Gyroelongated_Triangular_Cupola_J22     ! V=15, E=33, F=20
    public Exam_Johnson_Gyroelongated_Square_Cupola_J23         ! V=20, E=44, F=26
    public Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24     ! V=25, E=55, F=32
    public Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25    ! V=30, E=65, F=37
    public Exam_Johnson_Gyrobifastigium_J26                     ! V= 8, E=14, F= 8

contains

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Square_Pyramid_J1(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "21_Square_Pyramid_J1"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Square Pyramid J1"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 5
    geom.n_face = 5

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -10.1105d0,  9.2855d0,   4.4223d0 ]
    geom.iniP(2).pos(1:3) = [  -9.0792d0, -4.0479d0, -10.4491d0 ]
    geom.iniP(3).pos(1:3) = [  -1.3014d0, -8.4125d0,   7.4522d0 ]
    geom.iniP(4).pos(1:3) = [   9.7299d0,  8.2541d0,   6.7229d0 ]
    geom.iniP(5).pos(1:3) = [  10.7612d0, -5.0792d0,  -8.1485d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 2, 5, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 2, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 4, 1, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 5, 4, 3 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 5, 2, 1, 4 ]
end subroutine Exam_Johnson_Square_Pyramid_J1

! ---------------------------------------------------------------------------------------

! Example of pentagonal pyramid - J2
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Pyramid_J2(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "22_Penta_Pyramid_J2"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Penta Pyramid J2"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 6
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -13.9168d0,  -1.6024d0,  9.8118d0 ]
    geom.iniP(2).pos(1:3) = [  -5.2771d0,  15.6346d0,  4.4974d0 ]
    geom.iniP(3).pos(1:3) = [  -4.2653d0,  -0.2210d0, -7.6508d0 ]
    geom.iniP(4).pos(1:3) = [  -2.1451d0, -16.5638d0,  3.6813d0 ]
    geom.iniP(5).pos(1:3) = [  11.8342d0,  11.3262d0, -4.9176d0 ]
    geom.iniP(6).pos(1:3) = [  13.7699d0,  -8.5736d0, -5.4220d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 4, 1, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 6, 4, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 5, 6, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 5, 3 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 1, 2, 3 ]
    geom.face(6).n_poi = 5; allocate(geom.face(6).poi(5)); geom.face(6).poi(1:5) = [ 1, 4, 6, 5, 2 ]
end subroutine Exam_Johnson_Pentagonal_Pyramid_J2

! ---------------------------------------------------------------------------------------

! Example of triangular cupola - J3
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Triangular_Cupola_J3(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Trir_Cupola_J3"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tri Cupola J3"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 9
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -17.5081d0,  10.0668d0,   4.6647d0 ]
    geom.iniP(2).pos(1:3) = [ -14.3927d0,   4.3854d0, -14.2567d0 ]
    geom.iniP(3).pos(1:3) = [ -13.0636d0,  -8.9941d0,   0.5497d0 ]
    geom.iniP(4).pos(1:3) = [  -2.2041d0,  10.8607d0,  17.5163d0 ]
    geom.iniP(5).pos(1:3) = [   2.2403d0,  -8.2002d0,  13.4012d0 ]
    geom.iniP(6).pos(1:3) = [   4.0267d0,  -0.5021d0, -20.3265d0 ]
    geom.iniP(7).pos(1:3) = [   5.3557d0, -13.8816d0,  -5.5202d0 ]
    geom.iniP(8).pos(1:3) = [  16.2152d0,   5.9732d0,  11.4464d0 ]
    geom.iniP(9).pos(1:3) = [  19.3306d0,   0.2918d0,  -7.4750d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 3, 7, 5 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 7, 6, 9 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 5, 8, 4 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 3, 1, 2 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 7, 3, 2, 6 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 5, 7, 9, 8 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [ 3, 5, 4, 1 ]
    geom.face(8).n_poi = 6; allocate(geom.face(8).poi(6)); geom.face(8).poi(1:6) = [ 1, 4, 8, 9, 6, 2 ]
end subroutine Exam_Johnson_Triangular_Cupola_J3

! ---------------------------------------------------------------------------------------

! Example of square cupola - J4
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Square_Cupola_J4(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Square_Cupola_J4"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Square Cupola J4"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 12
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -14.2820d0,   9.4779d0, -20.2796d0 ]
    geom.iniP( 2).pos(1:3) = [ -13.9347d0, -10.5173d0, -20.0069d0 ]
    geom.iniP( 3).pos(1:3) = [ -13.9144d0,   9.7569d0,  -0.2849d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.5672d0, -10.2382d0,  -0.0122d0 ]
    geom.iniP( 5).pos(1:3) = [  -4.3469d0,  23.9272d0, -10.6612d0 ]
    geom.iniP( 6).pos(1:3) = [  -3.5085d0, -24.3455d0, -10.0030d0 ]
    geom.iniP( 7).pos(1:3) = [   0.4831d0,  10.1962d0,  13.5902d0 ]
    geom.iniP( 8).pos(1:3) = [   0.8304d0,  -9.7990d0,  13.8629d0 ]
    geom.iniP( 9).pos(1:3) = [  10.0506d0,  24.3665d0,   3.2140d0 ]
    geom.iniP(10).pos(1:3) = [  10.8890d0, -23.9062d0,   3.8722d0 ]
    geom.iniP(11).pos(1:3) = [  20.4767d0,  10.5383d0,  13.2179d0 ]
    geom.iniP(12).pos(1:3) = [  20.8240d0,  -9.4568d0,  13.4906d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 4,  2,  6 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 8, 10, 12 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 7, 11,  9 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 3,  5,  1 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 3,  4,  8,  7 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 4,  3,  1,  2 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 8,  4,  6, 10 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 7,  8, 12, 11 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 3,  7,  9,  5 ]
    geom.face(10).n_poi = 8; allocate(geom.face(10).poi(8)); geom.face(10).poi(1:8) = [ 5,  9, 11, 12, 10, 6, 2, 1 ]
end subroutine Exam_Johnson_Square_Cupola_J4

! ---------------------------------------------------------------------------------------

! Example of pentagonal cupola - J5
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Cupola_J5(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "23_Penta_Cupola_J5"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Penta Cupola J5"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 15
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -27.8521d0,   3.4402d0, -16.4903d0 ]
    geom.iniP( 2).pos(1:3) = [ -24.1621d0, -16.1327d0, -14.6775d0 ]
    geom.iniP( 3).pos(1:3) = [ -20.3511d0,  21.6903d0, -13.2237d0 ]
    geom.iniP( 4).pos(1:3) = [ -17.0151d0,   7.0047d0,  -0.0630d0 ]
    geom.iniP( 5).pos(1:3) = [ -13.3251d0, -12.5683d0,   1.7497d0 ]
    geom.iniP( 6).pos(1:3) = [ -10.6906d0, -29.5523d0,  -8.4779d0 ]
    geom.iniP( 7).pos(1:3) = [  -4.5242d0,  31.6466d0,  -6.1256d0 ]
    geom.iniP( 8).pos(1:3) = [  -1.1882d0,  16.9610d0,   7.0351d0 ]
    geom.iniP( 9).pos(1:3) = [   4.7823d0, -14.7087d0,   9.9682d0 ]
    geom.iniP(10).pos(1:3) = [   7.4168d0, -31.6927d0,  -0.2594d0 ]
    geom.iniP(11).pos(1:3) = [  12.2833d0,   3.5414d0,  13.2348d0 ]
    geom.iniP(12).pos(1:3) = [  13.5831d0,  29.5063d0,   2.0929d0 ]
    geom.iniP(13).pos(1:3) = [  23.2436d0, -21.7363d0,   6.8388d0 ]
    geom.iniP(14).pos(1:3) = [  27.0546d0,  16.0867d0,   8.2926d0 ]
    geom.iniP(15).pos(1:3) = [  30.7447d0,  -3.4862d0,  10.1053d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [  5,  2,  6 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [  9, 10, 13 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 11, 15, 14 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [  8, 12,  7 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [  4,  3,  1 ]
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi( 4)); geom.face( 6).poi(1: 4) = [  5,  4,  1,  2 ]
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi( 4)); geom.face( 7).poi(1: 4) = [  9,  5,  6, 10 ]
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi( 4)); geom.face( 8).poi(1: 4) = [ 11,  9, 13, 15 ]
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi( 4)); geom.face( 9).poi(1: 4) = [  8, 11, 14, 12 ]
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi( 4)); geom.face(10).poi(1: 4) = [  4,  8,  7,  3 ]
    geom.face(11).n_poi =  5; allocate(geom.face(11).poi( 5)); geom.face(11).poi(1: 5) = [  4,  5,  9, 11,  8 ]
    geom.face(12).n_poi = 10; allocate(geom.face(12).poi(10)); geom.face(12).poi(1:10) = [  3,  7, 12, 14, 15, 13, 10, 6, 2, 1 ]
end subroutine Exam_Johnson_Pentagonal_Cupola_J5

! ---------------------------------------------------------------------------------------

! Example of pentagonal rotunda - J6
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Rotunda_J6(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Penta_Rotunda_J6"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Penta Rotunda J6"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 20
    geom.n_face = 17

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -27.5191d0, -12.0355d0, -16.4029d0 ]
    geom.iniP( 2).pos(1:3) = [ -26.8439d0,  -7.8632d0,   3.1454d0 ]
    geom.iniP( 3).pos(1:3) = [ -21.8689d0, -26.1084d0,  -3.3634d0 ]
    geom.iniP( 4).pos(1:3) = [ -21.3804d0,   4.1651d0, -26.3958d0 ]
    geom.iniP( 5).pos(1:3) = [ -20.2880d0,  10.9160d0,   5.2341d0 ]
    geom.iniP( 6).pos(1:3) = [ -16.9114d0,  18.3499d0, -13.0234d0 ]
    geom.iniP( 7).pos(1:3) = [ -14.6378d0,  -3.1570d0,  18.2736d0 ]
    geom.iniP( 8).pos(1:3) = [  -6.5880d0, -32.6783d0,   7.7421d0 ]
    geom.iniP( 9).pos(1:3) = [  -5.7976d0,  16.3054d0, -29.5251d0 ]
    geom.iniP(10).pos(1:3) = [  -4.7053d0,  23.0562d0,   2.1048d0 ]
    geom.iniP(11).pos(1:3) = [  -2.1190d0, -18.4935d0,  21.1146d0 ]
    geom.iniP(12).pos(1:3) = [   4.4370d0,   0.2858d0,  23.2032d0 ]
    geom.iniP(13).pos(1:3) = [  10.5756d0,  16.4864d0,  13.2104d0 ]
    geom.iniP(14).pos(1:3) = [  12.4868d0, -29.2355d0,  12.6718d0 ]
    geom.iniP(15).pos(1:3) = [  13.2772d0,  19.7481d0, -24.5955d0 ]
    geom.iniP(16).pos(1:3) = [  13.9523d0,  23.9203d0,  -5.0472d0 ]
    geom.iniP(17).pos(1:3) = [  23.0945d0,   1.1499d0,  16.0513d0 ]
    geom.iniP(18).pos(1:3) = [  28.0696d0, -17.0953d0,   9.5425d0 ]
    geom.iniP(19).pos(1:3) = [  28.5581d0,  13.1783d0, -13.4899d0 ]
    geom.iniP(20).pos(1:3) = [  34.2082d0,  -0.8946d0,  -0.4504d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [ 12, 17, 13 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 17, 18, 20 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 13, 16, 10 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 16, 19, 15 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 10,  6,  5 ]
    geom.face( 6).n_poi =  3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [  6,  9,  4 ]
    geom.face( 7).n_poi =  3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [  5,  2,  7 ]
    geom.face( 8).n_poi =  3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [  2,  1,  3 ]
    geom.face( 9).n_poi =  3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [  7, 11, 12 ]
    geom.face(10).n_poi =  3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [ 11,  8, 14 ]
    geom.face(11).n_poi =  5; allocate(geom.face(11).poi( 5)); geom.face(11).poi(1: 5) = [ 12, 13, 10,  5,  7 ]
    geom.face(12).n_poi =  5; allocate(geom.face(12).poi( 5)); geom.face(12).poi(1: 5) = [ 12, 11, 14, 18, 17 ]
    geom.face(13).n_poi =  5; allocate(geom.face(13).poi( 5)); geom.face(13).poi(1: 5) = [ 13, 17, 20, 19, 16 ]
    geom.face(14).n_poi =  5; allocate(geom.face(14).poi( 5)); geom.face(14).poi(1: 5) = [ 10, 16, 15,  9,  6 ]
    geom.face(15).n_poi =  5; allocate(geom.face(15).poi( 5)); geom.face(15).poi(1: 5) = [  5,  6,  4,  1,  2 ]
    geom.face(16).n_poi =  5; allocate(geom.face(16).poi( 5)); geom.face(16).poi(1: 5) = [  7,  2,  3,  8, 11 ]
    geom.face(17).n_poi = 10; allocate(geom.face(17).poi(10)); geom.face(17).poi(1:10) = [  3,  1,  4,  9, 15, 19, 20, 18, 14, 8 ]
end subroutine Exam_Johnson_Pentagonal_Rotunda_J6

! ---------------------------------------------------------------------------------------

! Example of elongated triangular cupola - J18
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Triangular_Cupola_J18(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Elon_Tri_Cupola_J18"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elong Tri Cupola J18"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 15
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -19.3986d0,   1.1770d0,   6.6873d0 ]
    geom.iniP( 2).pos(1:3) = [ -15.9208d0,   0.3831d0, -12.9920d0 ]
    geom.iniP( 3).pos(1:3) = [ -13.6126d0, -17.8838d0,   8.4788d0 ]
    geom.iniP( 4).pos(1:3) = [ -13.2535d0,  20.2069d0,   7.0056d0 ]
    geom.iniP( 5).pos(1:3) = [ -10.1349d0, -18.6777d0, -11.2005d0 ]
    geom.iniP( 6).pos(1:3) = [  -9.7757d0,  19.4130d0, -12.6737d0 ]
    geom.iniP( 7).pos(1:3) = [  -4.9323d0,  -3.7104d0,  19.6040d0 ]
    geom.iniP( 8).pos(1:3) = [   1.2128d0,  15.3194d0,  19.9222d0 ]
    geom.iniP( 9).pos(1:3) = [   2.0232d0,  -5.2983d0, -19.7546d0 ]
    geom.iniP(10).pos(1:3) = [   4.3314d0, -23.5653d0,   1.7161d0 ]
    geom.iniP(11).pos(1:3) = [   8.1683d0,  13.7316d0, -19.4364d0 ]
    geom.iniP(12).pos(1:3) = [  13.0117d0,  -9.3918d0,  12.8413d0 ]
    geom.iniP(13).pos(1:3) = [  16.4895d0, -10.1858d0,  -6.8380d0 ]
    geom.iniP(14).pos(1:3) = [  19.1569d0,   9.6381d0,  13.1596d0 ]
    geom.iniP(15).pos(1:3) = [  22.6346d0,   8.8441d0,  -6.5197d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  5, 10,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 10, 13, 12 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  3,  7,  1 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  5,  2,  9 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  1,  4,  6,  2 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  2,  6, 11,  9 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  9, 11, 15, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 13, 15, 14, 12 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 12, 14,  8,  7 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  7,  8,  4,  1 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 10,  5,  9, 13 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  3, 10, 12,  7 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  5,  3,  1,  2 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [ 15, 11,  6,  4, 8, 14 ]
end subroutine Exam_Johnson_Elongated_Triangular_Cupola_J18

! ---------------------------------------------------------------------------------------

! Example of elongated square cupola - J19
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Square_Cupola_J19(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Elong_Square_Cupola_J19"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elon Square Cupola J19"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 20
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -24.5202d0,   3.1911d0,  -9.9080d0 ]
    geom.iniP( 2).pos(1:3) = [ -21.8375d0,  -6.3764d0,   7.4491d0 ]
    geom.iniP( 3).pos(1:3) = [ -21.8162d0,  13.6173d0,   6.9439d0 ]
    geom.iniP( 4).pos(1:3) = [ -14.3984d0, -11.2065d0, -19.4086d0 ]
    geom.iniP( 5).pos(1:3) = [ -14.3683d0,  17.0689d0, -20.1231d0 ]
    geom.iniP( 6).pos(1:3) = [ -11.7157d0, -20.7740d0,  -2.0515d0 ]
    geom.iniP( 7).pos(1:3) = [ -11.6642d0,  27.4950d0,  -3.2712d0 ]
    geom.iniP( 8).pos(1:3) = [  -7.8915d0,  -6.0291d0,  21.7806d0 ]
    geom.iniP( 9).pos(1:3) = [  -7.8702d0,  13.9645d0,  21.2754d0 ]
    geom.iniP(10).pos(1:3) = [  -4.2465d0,   2.6713d0, -29.6236d0 ]
    geom.iniP(11).pos(1:3) = [   2.2303d0, -20.4267d0,  12.2800d0 ]
    geom.iniP(12).pos(1:3) = [   2.2817d0,  27.8424d0,  11.0603d0 ]
    geom.iniP(13).pos(1:3) = [   2.6201d0, -21.1415d0, -15.9926d0 ]
    geom.iniP(14).pos(1:3) = [   9.1482d0,   4.0295d0,  24.6914d0 ]
    geom.iniP(15).pos(1:3) = [  12.7720d0,  -7.2638d0, -26.2077d0 ]
    geom.iniP(16).pos(1:3) = [  16.5660d0, -20.7942d0,  -1.6611d0 ]
    geom.iniP(17).pos(1:3) = [  19.2700d0, -10.3681d0,  15.1908d0 ]
    geom.iniP(18).pos(1:3) = [  19.3002d0,  17.9072d0,  14.4763d0 ]
    geom.iniP(19).pos(1:3) = [  26.7179d0,  -6.9165d0, -11.8762d0 ]
    geom.iniP(20).pos(1:3) = [  29.4221d0,   3.5097d0,   4.9757d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 11, 16, 17 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  8, 14,  9 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  3,  1 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  4, 13 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  3,  7,  5,  1 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  1,  5, 10,  4 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  4, 10, 15, 13 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 13, 15, 19, 16 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 16, 19, 20, 17 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 17, 20, 18, 14 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 14, 18, 12,  9 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  9, 12,  7,  3 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  6, 11,  8,  2 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 11,  6, 13, 16 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  8, 11, 17, 14 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  2,  8,  9,  3 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  6,  2,  1,  4 ]
    geom.face(18).n_poi = 8; allocate(geom.face(18).poi(8)); geom.face(18).poi(1:8) = [ 19, 15, 10,  5, 7, 12, 18, 20 ]
end subroutine Exam_Johnson_Elongated_Square_Cupola_J19

! ---------------------------------------------------------------------------------------

! Example of elongated pentagonal cupola - J20
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Cupola_J20(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Elon_Penta_Cupola_J20"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elon Penta Cupola J20"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 25
    geom.n_face = 22

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -30.1850d0,   9.7035d0,  -8.7581d0 ]
    geom.iniP( 2).pos(1:3) = [ -27.0859d0,  -8.4039d0, -16.6650d0 ]
    geom.iniP( 3).pos(1:3) = [ -22.9724d0,  23.1749d0,   4.1454d0 ]
    geom.iniP( 4).pos(1:3) = [ -22.9406d0,  -5.0679d0,   2.6142d0 ]
    geom.iniP( 5).pos(1:3) = [ -19.3708d0,  17.9564d0, -23.4190d0 ]
    geom.iniP( 6).pos(1:3) = [ -16.2717d0,  -0.1510d0, -31.3260d0 ]
    geom.iniP( 7).pos(1:3) = [ -15.7280d0,   8.4036d0,  15.5177d0 ]
    geom.iniP( 8).pos(1:3) = [ -14.8588d0, -24.2307d0, -16.5554d0 ]
    geom.iniP( 9).pos(1:3) = [ -12.1582d0,  31.4279d0, -10.5156d0 ]
    geom.iniP(10).pos(1:3) = [ -10.7135d0, -20.8947d0,   2.7239d0 ]
    geom.iniP(11).pos(1:3) = [  -8.2031d0,  26.8649d0,  17.1166d0 ]
    geom.iniP(12).pos(1:3) = [  -4.0446d0, -15.9778d0, -31.2164d0 ]
    geom.iniP(13).pos(1:3) = [   0.9567d0,   0.9026d0,  23.6021d0 ]
    geom.iniP(14).pos(1:3) = [   1.8259d0, -31.7317d0,  -8.4709d0 ]
    geom.iniP(15).pos(1:3) = [   2.6111d0,  35.1178d0,   2.4556d0 ]
    geom.iniP(16).pos(1:3) = [   4.0558d0, -17.2047d0,  15.6951d0 ]
    geom.iniP(17).pos(1:3) = [   8.4816d0,  19.3639d0,  25.2011d0 ]
    geom.iniP(18).pos(1:3) = [  12.6400d0, -23.4788d0, -23.1319d0 ]
    geom.iniP(19).pos(1:3) = [  16.5952d0, -28.0418d0,   4.5003d0 ]
    geom.iniP(20).pos(1:3) = [  19.2958d0,  27.6168d0,  10.5401d0 ]
    geom.iniP(21).pos(1:3) = [  20.7086d0,   3.5371d0,  25.3107d0 ]
    geom.iniP(22).pos(1:3) = [  23.8078d0, -14.5703d0,  17.4037d0 ]
    geom.iniP(23).pos(1:3) = [  27.4093d0, -19.7888d0, -10.1607d0 ]
    geom.iniP(24).pos(1:3) = [  31.5228d0,  11.7900d0,  10.6498d0 ]
    geom.iniP(25).pos(1:3) = [  34.6221d0,  -6.3173d0,   2.7428d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [ 16, 19, 22 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 13, 21, 17 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [  7, 11,  3 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [  4,  1,  2 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 10,  8, 14 ]
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi( 4)); geom.face( 6).poi(1: 4) = [  3,  9,  5,  1 ]
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi( 4)); geom.face( 7).poi(1: 4) = [  1,  5,  6,  2 ]
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi( 4)); geom.face( 8).poi(1: 4) = [  2,  6, 12,  8 ]
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi( 4)); geom.face( 9).poi(1: 4) = [  8, 12, 18, 14 ]
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi( 4)); geom.face(10).poi(1: 4) = [ 14, 18, 23, 19 ]
    geom.face(11).n_poi =  4; allocate(geom.face(11).poi( 4)); geom.face(11).poi(1: 4) = [ 19, 23, 25, 22 ]
    geom.face(12).n_poi =  4; allocate(geom.face(12).poi( 4)); geom.face(12).poi(1: 4) = [ 22, 25, 24, 21 ]
    geom.face(13).n_poi =  4; allocate(geom.face(13).poi( 4)); geom.face(13).poi(1: 4) = [ 21, 24, 20, 17 ]
    geom.face(14).n_poi =  4; allocate(geom.face(14).poi( 4)); geom.face(14).poi(1: 4) = [ 17, 20, 15, 11 ]
    geom.face(15).n_poi =  4; allocate(geom.face(15).poi( 4)); geom.face(15).poi(1: 4) = [ 11, 15,  9,  3 ]
    geom.face(16).n_poi =  4; allocate(geom.face(16).poi( 4)); geom.face(16).poi(1: 4) = [ 16, 10, 14, 19 ]
    geom.face(17).n_poi =  4; allocate(geom.face(17).poi( 4)); geom.face(17).poi(1: 4) = [ 13, 16, 22, 21 ]
    geom.face(18).n_poi =  4; allocate(geom.face(18).poi( 4)); geom.face(18).poi(1: 4) = [  7, 13, 17, 11 ]
    geom.face(19).n_poi =  4; allocate(geom.face(19).poi( 4)); geom.face(19).poi(1: 4) = [  4,  7,  3,  1 ]
    geom.face(20).n_poi =  4; allocate(geom.face(20).poi( 4)); geom.face(20).poi(1: 4) = [ 10,  4,  2,  8 ]
    geom.face(21).n_poi =  5; allocate(geom.face(21).poi( 5)); geom.face(21).poi(1: 5) = [ 10, 16, 13,  7,  4 ]
    geom.face(22).n_poi = 10; allocate(geom.face(22).poi(10)); geom.face(22).poi(1:10) = [ 23, 18, 12,  6,  5, 9, 15, 20, 24, 25 ]
end subroutine Exam_Johnson_Elongated_Pentagonal_Cupola_J20

! ---------------------------------------------------------------------------------------

! Example of elongated pentagonal rotunda - J21
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Rotunda_J21(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Elon_Penta_Rotunda_J21"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elon Penta Rotunda J21"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 30
    geom.n_face = 27

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -31.3933d0,   4.7766d0,  -3.6992d0 ]
    geom.iniP( 2).pos(1:3) = [ -27.5261d0,   1.6602d0,  15.6743d0 ]
    geom.iniP( 3).pos(1:3) = [ -26.7983d0, -11.9322d0, -13.6844d0 ]
    geom.iniP( 4).pos(1:3) = [ -23.8423d0,  19.5337d0,   7.4906d0 ]
    geom.iniP( 5).pos(1:3) = [ -23.1145d0,   5.9412d0, -21.8681d0 ]
    geom.iniP( 6).pos(1:3) = [ -20.5410d0, -16.9746d0,  17.6626d0 ]
    geom.iniP( 7).pos(1:3) = [ -20.0912d0, -25.3752d0,  -0.4820d0 ]
    geom.iniP( 8).pos(1:3) = [ -16.0837d0, -20.7384d0, -28.0944d0 ]
    geom.iniP( 9).pos(1:3) = [ -12.9900d0,  -2.2175d0,  28.8525d0 ]
    geom.iniP(10).pos(1:3) = [ -12.3998d0,  -2.8650d0, -36.2782d0 ]
    geom.iniP(11).pos(1:3) = [ -10.8966d0,  29.8186d0,  -3.7625d0 ]
    geom.iniP(12).pos(1:3) = [ -10.4468d0,  21.4180d0, -21.9071d0 ]
    geom.iniP(13).pos(1:3) = [  -9.3766d0, -34.1814d0, -14.8920d0 ]
    geom.iniP(14).pos(1:3) = [  -7.0294d0,  26.7022d0,  15.6111d0 ]
    geom.iniP(15).pos(1:3) = [  -5.5552d0, -29.2530d0,  12.6962d0 ]
    geom.iniP(16).pos(1:3) = [  -0.3224d0,  13.2592d0,  28.8135d0 ]
    geom.iniP(17).pos(1:3) = [   0.2678d0,  12.6117d0, -36.3171d0 ]
    geom.iniP(18).pos(1:3) = [   5.1595d0, -38.0593d0,  -1.7138d0 ]
    geom.iniP(19).pos(1:3) = [   6.3660d0,  28.5865d0, -13.7867d0 ]
    geom.iniP(20).pos(1:3) = [   6.6627d0,  -5.3756d0,  30.8018d0 ]
    geom.iniP(21).pos(1:3) = [  11.2577d0, -22.0844d0,  20.8167d0 ]
    geom.iniP(22).pos(1:3) = [  12.6233d0,  23.5442d0,  17.5603d0 ]
    geom.iniP(23).pos(1:3) = [  17.0807d0,  19.7803d0, -28.1967d0 ]
    geom.iniP(24).pos(1:3) = [  20.9021d0,  24.7088d0,  -0.6085d0 ]
    geom.iniP(25).pos(1:3) = [  21.9723d0, -30.8907d0,   6.4067d0 ]
    geom.iniP(26).pos(1:3) = [  23.9253d0,  -6.6076d0,  20.7776d0 ]
    geom.iniP(27).pos(1:3) = [  27.6092d0,  11.2657d0,  12.5939d0 ]
    geom.iniP(28).pos(1:3) = [  31.6168d0,  15.9025d0, -15.0185d0 ]
    geom.iniP(29).pos(1:3) = [  34.6400d0, -15.4139d0,   6.3676d0 ]
    geom.iniP(30).pos(1:3) = [  38.3238d0,   2.4595d0,  -1.8161d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [  9, 20, 16 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 20, 21, 26 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 16, 22, 14 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 22, 27, 24 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 14, 11,  4 ]
    geom.face( 6).n_poi =  3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [ 11, 19, 12 ]
    geom.face( 7).n_poi =  3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [  4,  1,  2 ]
    geom.face( 8).n_poi =  3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [  1,  5,  3 ]
    geom.face( 9).n_poi =  3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [  2,  6,  9 ]
    geom.face(10).n_poi =  3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [  6,  7, 15 ]
    geom.face(11).n_poi =  4; allocate(geom.face(11).poi( 4)); geom.face(11).poi(1: 4) = [ 12, 17, 10,  5 ]
    geom.face(12).n_poi =  4; allocate(geom.face(12).poi( 4)); geom.face(12).poi(1: 4) = [  5, 10,  8,  3 ]
    geom.face(13).n_poi =  4; allocate(geom.face(13).poi( 4)); geom.face(13).poi(1: 4) = [  3,  8, 13,  7 ]
    geom.face(14).n_poi =  4; allocate(geom.face(14).poi( 4)); geom.face(14).poi(1: 4) = [  7, 13, 18, 15 ]
    geom.face(15).n_poi =  4; allocate(geom.face(15).poi( 4)); geom.face(15).poi(1: 4) = [ 15, 18, 25, 21 ]
    geom.face(16).n_poi =  4; allocate(geom.face(16).poi( 4)); geom.face(16).poi(1: 4) = [ 21, 25, 29, 26 ]
    geom.face(17).n_poi =  4; allocate(geom.face(17).poi( 4)); geom.face(17).poi(1: 4) = [ 26, 29, 30, 27 ]
    geom.face(18).n_poi =  4; allocate(geom.face(18).poi( 4)); geom.face(18).poi(1: 4) = [ 27, 30, 28, 24 ]
    geom.face(19).n_poi =  4; allocate(geom.face(19).poi( 4)); geom.face(19).poi(1: 4) = [ 24, 28, 23, 19 ]
    geom.face(20).n_poi =  4; allocate(geom.face(20).poi( 4)); geom.face(20).poi(1: 4) = [ 19, 23, 17, 12 ]
    geom.face(21).n_poi =  5; allocate(geom.face(21).poi( 5)); geom.face(21).poi(1: 5) = [  9, 16, 14,  4,  2 ]
    geom.face(22).n_poi =  5; allocate(geom.face(22).poi( 5)); geom.face(22).poi(1: 5) = [  9,  6, 15, 21, 20 ]
    geom.face(23).n_poi =  5; allocate(geom.face(23).poi( 5)); geom.face(23).poi(1: 5) = [ 16, 20, 26, 27, 22 ]
    geom.face(24).n_poi =  5; allocate(geom.face(24).poi( 5)); geom.face(24).poi(1: 5) = [ 14, 22, 24, 19, 11 ]
    geom.face(25).n_poi =  5; allocate(geom.face(25).poi( 5)); geom.face(25).poi(1: 5) = [  4, 11, 12,  5,  1 ]
    geom.face(26).n_poi =  5; allocate(geom.face(26).poi( 5)); geom.face(26).poi(1: 5) = [  2,  1,  3,  7,  6 ]
    geom.face(27).n_poi = 10; allocate(geom.face(27).poi(10)); geom.face(27).poi(1:10) = [ 25, 18, 13,  8, 10, 17, 23, 28, 30, 29 ]
end subroutine Exam_Johnson_Elongated_Pentagonal_Rotunda_J21

! ---------------------------------------------------------------------------------------

! Example of gyroelongated triangular cupola - J22
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Triangular_Cupola_J22(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Gyro_Tri_Cupola_J22"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyro Tri Cupola J22"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 15
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -19.0184d0,   1.4823d0,   6.9937d0 ]
    geom.iniP( 2).pos(1:3) = [ -17.2045d0,  15.2402d0,  -7.4088d0 ]
    geom.iniP( 3).pos(1:3) = [ -15.9030d0,  -4.1991d0, -11.9277d0 ]
    geom.iniP( 4).pos(1:3) = [ -14.5740d0, -17.5785d0,   2.8786d0 ]
    geom.iniP( 5).pos(1:3) = [ -10.1675d0,  18.9787d0,  10.9353d0 ]
    geom.iniP( 6).pos(1:3) = [  -4.7715d0,   9.1382d0, -21.8374d0 ]
    geom.iniP( 7).pos(1:3) = [  -3.7145d0,   2.2762d0,  19.8452d0 ]
    geom.iniP( 8).pos(1:3) = [   0.7299d0, -16.7846d0,  15.7301d0 ]
    geom.iniP( 9).pos(1:3) = [   2.5163d0,  -9.0866d0, -17.9975d0 ]
    geom.iniP(10).pos(1:3) = [   3.8453d0, -22.4661d0,  -3.1912d0 ]
    geom.iniP(11).pos(1:3) = [   9.3026d0,  16.6153d0,  14.8507d0 ]
    geom.iniP(12).pos(1:3) = [  14.6986d0,   6.7748d0, -17.9221d0 ]
    geom.iniP(13).pos(1:3) = [  14.7048d0,  -2.6113d0,  13.7753d0 ]
    geom.iniP(14).pos(1:3) = [  17.8202d0,  -8.2927d0,  -5.1460d0 ]
    geom.iniP(15).pos(1:3) = [  21.7357d0,  10.5133d0,   0.4220d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  2,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  3,  6,  9 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  9,  6, 12 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  9, 12, 14 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 14, 12, 15 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  4, 10,  8 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10,  9, 14 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  8, 13,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4,  1,  3 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 14, 15, 13 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 13, 15, 11 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 13, 11,  7 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  7, 11,  5 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  7,  5,  1 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  1,  5,  2 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 10,  4,  3,  9 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  8, 10, 14, 13 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [  4,  8,  7,  1 ]
    geom.face(20).n_poi = 6; allocate(geom.face(20).poi(6)); geom.face(20).poi(1:6) = [ 12,  6,  2,  5, 11, 15 ]
end subroutine Exam_Johnson_Gyroelongated_Triangular_Cupola_J22

! ---------------------------------------------------------------------------------------

! Example of gyroelongated square cupola - J23
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Square_Cupola_J23(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "24_Gyro_Square_Cupola_J23"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyro Square Cupola J23"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 20
    geom.n_face = 26

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -26.1543d0,   8.6723d0,  -9.8268d0 ]
    geom.iniP( 2).pos(1:3) = [ -24.3469d0,   3.9666d0,   9.5276d0 ]
    geom.iniP( 3).pos(1:3) = [ -23.9077d0, -10.4309d0,  -4.3475d0 ]
    geom.iniP( 4).pos(1:3) = [ -18.9083d0,  22.1068d0,   3.0965d0 ]
    geom.iniP( 5).pos(1:3) = [ -18.0967d0,  -4.4964d0, -22.5414d0 ]
    geom.iniP( 6).pos(1:3) = [ -10.5188d0,  14.3928d0,  19.5316d0 ]
    geom.iniP( 7).pos(1:3) = [ -10.1766d0,  -5.6009d0,  19.9039d0 ]
    geom.iniP( 8).pos(1:3) = [  -9.7374d0, -19.9985d0,   6.0288d0 ]
    geom.iniP( 9).pos(1:3) = [  -9.4583d0, -20.3660d0, -13.9660d0 ]
    geom.iniP(10).pos(1:3) = [  -0.6033d0,  27.9373d0,   8.6582d0 ]
    geom.iniP(11).pos(1:3) = [   0.5446d0,  -9.6853d0, -27.5992d0 ]
    geom.iniP(12).pos(1:3) = [   9.4764d0,  14.7400d0,  19.8043d0 ]
    geom.iniP(13).pos(1:3) = [   9.8186d0,  -5.2536d0,  20.1766d0 ]
    geom.iniP(14).pos(1:3) = [  10.2578d0, -19.6512d0,   6.3014d0 ]
    geom.iniP(15).pos(1:3) = [  10.5369d0, -20.0187d0, -13.6933d0 ]
    geom.iniP(16).pos(1:3) = [  18.0380d0,  22.7484d0,   3.6003d0 ]
    geom.iniP(17).pos(1:3) = [  18.8496d0,  -3.8548d0, -22.0376d0 ]
    geom.iniP(18).pos(1:3) = [  23.9258d0,   4.8050d0,  10.1859d0 ]
    geom.iniP(19).pos(1:3) = [  24.3650d0,  -9.5926d0,  -3.6893d0 ]
    geom.iniP(20).pos(1:3) = [  26.0956d0,   9.5797d0,  -9.1143d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  2,  1,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  1,  5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  3,  5,  9 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  9,  5, 11 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  9, 11, 15 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 15, 11, 17 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 15, 17, 19 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 19, 17, 20 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 14, 15, 19 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 13, 18, 12 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7,  6,  2 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  8,  3,  9 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 19, 20, 18 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 20, 16, 18 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 18, 16, 12 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 12, 16, 10 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 12, 10,  6 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  6, 10,  4 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  6,  4,  2 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  2,  4,  1 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [  8, 14, 13,  7 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 14,  8,  9, 15 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 13, 14, 19, 18 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [  7, 13, 12,  6 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [  8,  7,  2,  3 ]
    geom.face(26).n_poi = 8; allocate(geom.face(26).poi(8)); geom.face(26).poi(1:8) = [ 11,  5,  1,  4, 10, 16, 20, 17 ]
end subroutine Exam_Johnson_Gyroelongated_Square_Cupola_J23

! ---------------------------------------------------------------------------------------

! Example of gyroelongated pentagonal cupola - J24
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "25_Gyro_Penta_Cupola_J24"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyro Penta Cupola J24"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 25
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -32.0925d0,   8.3800d0, -10.1052d0 ]
    geom.iniP( 2).pos(1:3) = [ -31.7012d0,  -7.9497d0,   1.4357d0 ]
    geom.iniP( 3).pos(1:3) = [ -29.5608d0,  10.1577d0,   9.6542d0 ]
    geom.iniP( 4).pos(1:3) = [ -27.9833d0,  -9.4603d0, -18.1577d0 ]
    geom.iniP( 5).pos(1:3) = [ -23.9121d0,  24.9820d0,  -2.5251d0 ]
    geom.iniP( 6).pos(1:3) = [ -21.7448d0, -23.7765d0,  -5.6625d0 ]
    geom.iniP( 7).pos(1:3) = [ -17.0155d0,  -4.6136d0,  14.5964d0 ]
    geom.iniP( 8).pos(1:3) = [ -16.1412d0,  23.6292d0,  15.8538d0 ]
    geom.iniP( 9).pos(1:3) = [ -13.1543d0, -21.7245d0, -23.6067d0 ]
    geom.iniP(10).pos(1:3) = [  -7.0592d0, -20.4405d0,   7.4982d0 ]
    geom.iniP(11).pos(1:3) = [  -6.5669d0,  34.0042d0,   1.6872d0 ]
    geom.iniP(12).pos(1:3) = [  -3.5959d0,   8.8579d0,  20.7961d0 ]
    geom.iniP(13).pos(1:3) = [  -3.4947d0, -31.2775d0,  -8.9291d0 ]
    geom.iniP(14).pos(1:3) = [   3.4317d0,  27.3192d0,  17.6666d0 ]
    geom.iniP(15).pos(1:3) = [   6.7304d0, -23.7280d0, -24.3710d0 ]
    geom.iniP(16).pos(1:3) = [  12.5138d0, -16.7505d0,   9.3110d0 ]
    geom.iniP(17).pos(1:3) = [  13.3179d0,  32.0008d0,   0.9229d0 ]
    geom.iniP(18).pos(1:3) = [  14.6541d0,   1.3569d0,  17.5295d0 ]
    geom.iniP(19).pos(1:3) = [  16.0782d0, -27.5875d0,  -7.1163d0 ]
    geom.iniP(20).pos(1:3) = [  21.6818d0,  19.8182d0,  14.4000d0 ]
    geom.iniP(21).pos(1:3) = [  24.0757d0, -14.7057d0, -20.1586d0 ]
    geom.iniP(22).pos(1:3) = [  28.1469d0,  19.7366d0,  -4.5261d0 ]
    geom.iniP(23).pos(1:3) = [  29.4978d0, -14.1160d0,  -0.9166d0 ]
    geom.iniP(24).pos(1:3) = [  31.6382d0,   3.9914d0,   7.3019d0 ]
    geom.iniP(25).pos(1:3) = [  32.2559d0,   1.8963d0, -12.5786d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [  2,  4,  6 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [  6,  4,  9 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [  6,  9, 13 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 13,  9, 15 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 13, 15, 19 ]
    geom.face( 6).n_poi =  3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [ 19, 15, 21 ]
    geom.face( 7).n_poi =  3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [ 19, 21, 23 ]
    geom.face( 8).n_poi =  3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [ 23, 21, 25 ]
    geom.face( 9).n_poi =  3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [ 23, 25, 24 ]
    geom.face(10).n_poi =  3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [ 24, 25, 22 ]
    geom.face(11).n_poi =  3; allocate(geom.face(11).poi( 3)); geom.face(11).poi(1: 3) = [ 16, 19, 23 ]
    geom.face(12).n_poi =  3; allocate(geom.face(12).poi( 3)); geom.face(12).poi(1: 3) = [ 18, 24, 20 ]
    geom.face(13).n_poi =  3; allocate(geom.face(13).poi( 3)); geom.face(13).poi(1: 3) = [ 12, 14,  8 ]
    geom.face(14).n_poi =  3; allocate(geom.face(14).poi( 3)); geom.face(14).poi(1: 3) = [  7,  3,  2 ]
    geom.face(15).n_poi =  3; allocate(geom.face(15).poi( 3)); geom.face(15).poi(1: 3) = [ 10,  6, 13 ]
    geom.face(16).n_poi =  3; allocate(geom.face(16).poi( 3)); geom.face(16).poi(1: 3) = [ 24, 22, 20 ]
    geom.face(17).n_poi =  3; allocate(geom.face(17).poi( 3)); geom.face(17).poi(1: 3) = [ 20, 22, 17 ]
    geom.face(18).n_poi =  3; allocate(geom.face(18).poi( 3)); geom.face(18).poi(1: 3) = [ 20, 17, 14 ]
    geom.face(19).n_poi =  3; allocate(geom.face(19).poi( 3)); geom.face(19).poi(1: 3) = [ 14, 17, 11 ]
    geom.face(20).n_poi =  3; allocate(geom.face(20).poi( 3)); geom.face(20).poi(1: 3) = [ 14, 11,  8 ]
    geom.face(21).n_poi =  3; allocate(geom.face(21).poi( 3)); geom.face(21).poi(1: 3) = [  8, 11,  5 ]
    geom.face(22).n_poi =  3; allocate(geom.face(22).poi( 3)); geom.face(22).poi(1: 3) = [  8,  5,  3 ]
    geom.face(23).n_poi =  3; allocate(geom.face(23).poi( 3)); geom.face(23).poi(1: 3) = [  3,  5,  1 ]
    geom.face(24).n_poi =  3; allocate(geom.face(24).poi( 3)); geom.face(24).poi(1: 3) = [  3,  1,  2 ]
    geom.face(25).n_poi =  3; allocate(geom.face(25).poi( 3)); geom.face(25).poi(1: 3) = [  2,  1,  4 ]
    geom.face(26).n_poi =  4; allocate(geom.face(26).poi( 4)); geom.face(26).poi(1: 4) = [ 16, 10, 13, 19 ]
    geom.face(27).n_poi =  4; allocate(geom.face(27).poi( 4)); geom.face(27).poi(1: 4) = [ 18, 16, 23, 24 ]
    geom.face(28).n_poi =  4; allocate(geom.face(28).poi( 4)); geom.face(28).poi(1: 4) = [ 12, 18, 20, 14 ]
    geom.face(29).n_poi =  4; allocate(geom.face(29).poi( 4)); geom.face(29).poi(1: 4) = [  7, 12,  8,  3 ]
    geom.face(30).n_poi =  4; allocate(geom.face(30).poi( 4)); geom.face(30).poi(1: 4) = [ 10,  7,  2,  6 ]
    geom.face(31).n_poi =  5; allocate(geom.face(31).poi( 5)); geom.face(31).poi(1: 5) = [ 10, 16, 18, 12,  7 ]
    geom.face(32).n_poi = 10; allocate(geom.face(32).poi(10)); geom.face(32).poi(1:10) = [ 21, 15,  9,  4,  1, 5, 11, 17, 22, 25 ]
end subroutine Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24

! ---------------------------------------------------------------------------------------

! Example of gyroelongated pentagonal rotunda - J25
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Gyro_Penta_Rotunda_J25"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyro Penta Rotunda J25"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 30
    geom.n_face = 37

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -30.3609d0,  -6.5425d0,  -9.2432d0 ]
    geom.iniP( 2).pos(1:3) = [ -29.6858d0,  -2.3702d0,  10.3052d0 ]
    geom.iniP( 3).pos(1:3) = [ -24.7107d0, -20.6154d0,   3.7964d0 ]
    geom.iniP( 4).pos(1:3) = [ -24.2223d0,   9.6582d0, -19.2362d0 ]
    geom.iniP( 5).pos(1:3) = [ -23.7981d0, -24.2419d0, -15.8510d0 ]
    geom.iniP( 6).pos(1:3) = [ -23.5413d0,  -8.3261d0, -27.9599d0 ]
    geom.iniP( 7).pos(1:3) = [ -23.1299d0,  16.4091d0,  12.3939d0 ]
    geom.iniP( 8).pos(1:3) = [ -19.7532d0,  23.8431d0,  -5.8637d0 ]
    geom.iniP( 9).pos(1:3) = [ -17.4797d0,   2.3361d0,  25.4335d0 ]
    geom.iniP(10).pos(1:3) = [ -12.7939d0, -35.0946d0,  -3.1572d0 ]
    geom.iniP(11).pos(1:3) = [ -12.1216d0,   6.5736d0, -34.8585d0 ]
    geom.iniP(12).pos(1:3) = [  -9.4298d0, -27.1853d0,  14.9019d0 ]
    geom.iniP(13).pos(1:3) = [  -8.6394d0,  21.7985d0, -22.3654d0 ]
    geom.iniP(14).pos(1:3) = [  -7.5470d0,  28.5494d0,   9.2646d0 ]
    geom.iniP(15).pos(1:3) = [  -4.9608d0, -13.0004d0,  28.2744d0 ]
    geom.iniP(16).pos(1:3) = [   1.5952d0,   5.7788d0,  30.3631d0 ]
    geom.iniP(17).pos(1:3) = [   5.2680d0, -36.7384d0,   5.2731d0 ]
    geom.iniP(18).pos(1:3) = [   6.0991d0,  14.7661d0, -33.9123d0 ]
    geom.iniP(19).pos(1:3) = [   7.7339d0,  21.9795d0,  20.3702d0 ]
    geom.iniP(20).pos(1:3) = [   9.6451d0, -23.7426d0,  19.8316d0 ]
    geom.iniP(21).pos(1:3) = [  10.4355d0,  25.2413d0, -17.4358d0 ]
    geom.iniP(22).pos(1:3) = [  11.1106d0,  29.4135d0,   2.1126d0 ]
    geom.iniP(23).pos(1:3) = [  20.2528d0,   6.6430d0,  23.2112d0 ]
    geom.iniP(24).pos(1:3) = [  23.4887d0, -28.5460d0,   6.2195d0 ]
    geom.iniP(25).pos(1:3) = [  24.1610d0,  13.1221d0, -25.4820d0 ]
    geom.iniP(26).pos(1:3) = [  25.2279d0, -11.6022d0,  16.7023d0 ]
    geom.iniP(27).pos(1:3) = [  25.7164d0,  18.6714d0,  -6.3302d0 ]
    geom.iniP(28).pos(1:3) = [  31.3667d0,   4.5985d0,   6.7094d0 ]
    geom.iniP(29).pos(1:3) = [  34.9082d0, -13.6463d0,  -0.6792d0 ]
    geom.iniP(30).pos(1:3) = [  35.1653d0,   2.2695d0, -12.7881d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi( 3)); geom.face( 1).poi(1: 3) = [  3, 10, 12 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi( 3)); geom.face( 2).poi(1: 3) = [ 12, 10, 17 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi( 3)); geom.face( 3).poi(1: 3) = [ 12, 17, 20 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi( 3)); geom.face( 4).poi(1: 3) = [ 20, 17, 24 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi( 3)); geom.face( 5).poi(1: 3) = [ 20, 24, 26 ]
    geom.face( 6).n_poi =  3; allocate(geom.face( 6).poi( 3)); geom.face( 6).poi(1: 3) = [ 26, 24, 29 ]
    geom.face( 7).n_poi =  3; allocate(geom.face( 7).poi( 3)); geom.face( 7).poi(1: 3) = [ 26, 29, 28 ]
    geom.face( 8).n_poi =  3; allocate(geom.face( 8).poi( 3)); geom.face( 8).poi(1: 3) = [ 28, 29, 30 ]
    geom.face( 9).n_poi =  3; allocate(geom.face( 9).poi( 3)); geom.face( 9).poi(1: 3) = [ 28, 30, 27 ]
    geom.face(10).n_poi =  3; allocate(geom.face(10).poi( 3)); geom.face(10).poi(1: 3) = [ 27, 30, 25 ]
    geom.face(11).n_poi =  3; allocate(geom.face(11).poi( 3)); geom.face(11).poi(1: 3) = [ 16, 23, 19 ]
    geom.face(12).n_poi =  3; allocate(geom.face(12).poi( 3)); geom.face(12).poi(1: 3) = [ 23, 26, 28 ]
    geom.face(13).n_poi =  3; allocate(geom.face(13).poi( 3)); geom.face(13).poi(1: 3) = [ 19, 22, 14 ]
    geom.face(14).n_poi =  3; allocate(geom.face(14).poi( 3)); geom.face(14).poi(1: 3) = [ 22, 27, 21 ]
    geom.face(15).n_poi =  3; allocate(geom.face(15).poi( 3)); geom.face(15).poi(1: 3) = [ 14,  8,  7 ]
    geom.face(16).n_poi =  3; allocate(geom.face(16).poi( 3)); geom.face(16).poi(1: 3) = [  8, 13,  4 ]
    geom.face(17).n_poi =  3; allocate(geom.face(17).poi( 3)); geom.face(17).poi(1: 3) = [  7,  2,  9 ]
    geom.face(18).n_poi =  3; allocate(geom.face(18).poi( 3)); geom.face(18).poi(1: 3) = [  2,  1,  3 ]
    geom.face(19).n_poi =  3; allocate(geom.face(19).poi( 3)); geom.face(19).poi(1: 3) = [  9, 15, 16 ]
    geom.face(20).n_poi =  3; allocate(geom.face(20).poi( 3)); geom.face(20).poi(1: 3) = [ 15, 12, 20 ]
    geom.face(21).n_poi =  3; allocate(geom.face(21).poi( 3)); geom.face(21).poi(1: 3) = [ 27, 25, 21 ]
    geom.face(22).n_poi =  3; allocate(geom.face(22).poi( 3)); geom.face(22).poi(1: 3) = [ 21, 25, 18 ]
    geom.face(23).n_poi =  3; allocate(geom.face(23).poi( 3)); geom.face(23).poi(1: 3) = [ 21, 18, 13 ]
    geom.face(24).n_poi =  3; allocate(geom.face(24).poi( 3)); geom.face(24).poi(1: 3) = [ 13, 18, 11 ]
    geom.face(25).n_poi =  3; allocate(geom.face(25).poi( 3)); geom.face(25).poi(1: 3) = [ 13, 11,  4 ]
    geom.face(26).n_poi =  3; allocate(geom.face(26).poi( 3)); geom.face(26).poi(1: 3) = [  4, 11,  6 ]
    geom.face(27).n_poi =  3; allocate(geom.face(27).poi( 3)); geom.face(27).poi(1: 3) = [  4,  6,  1 ]
    geom.face(28).n_poi =  3; allocate(geom.face(28).poi( 3)); geom.face(28).poi(1: 3) = [  1,  6,  5 ]
    geom.face(29).n_poi =  3; allocate(geom.face(29).poi( 3)); geom.face(29).poi(1: 3) = [  1,  5,  3 ]
    geom.face(30).n_poi =  3; allocate(geom.face(30).poi( 3)); geom.face(30).poi(1: 3) = [  3,  5, 10 ]
    geom.face(31).n_poi =  5; allocate(geom.face(31).poi( 5)); geom.face(31).poi(1: 5) = [ 16, 19, 14,  7,  9 ]
    geom.face(32).n_poi =  5; allocate(geom.face(32).poi( 5)); geom.face(32).poi(1: 5) = [ 16, 15, 20, 26, 23 ]
    geom.face(33).n_poi =  5; allocate(geom.face(33).poi( 5)); geom.face(33).poi(1: 5) = [ 19, 23, 28, 27, 22 ]
    geom.face(34).n_poi =  5; allocate(geom.face(34).poi( 5)); geom.face(34).poi(1: 5) = [ 14, 22, 21, 13,  8 ]
    geom.face(35).n_poi =  5; allocate(geom.face(35).poi( 5)); geom.face(35).poi(1: 5) = [  7,  8,  4,  1,  2 ]
    geom.face(36).n_poi =  5; allocate(geom.face(36).poi( 5)); geom.face(36).poi(1: 5) = [  9,  2,  3, 12, 15 ]
    geom.face(37).n_poi = 10; allocate(geom.face(37).poi(10)); geom.face(37).poi(1:10) = [ 29, 24, 17, 10,  5, 6, 11, 18, 25, 30 ]
end subroutine Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25

! ---------------------------------------------------------------------------------------

! Example of gyrobifastigium - J26
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyrobifastigium_J26(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Gyrobifastigium_J26"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyrobifastigium J26"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set geometric type and view
    geom.n_iniP = 8
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -10.0000d0,  10.0000d0,   0.0000d0 ]
    geom.iniP(2).pos(1:3) = [  10.0000d0,  10.0000d0,   0.0000d0 ]
    geom.iniP(3).pos(1:3) = [  10.0000d0, -10.0000d0,   0.0000d0 ]
    geom.iniP(4).pos(1:3) = [ -10.0000d0, -10.0000d0,   0.0000d0 ]
    geom.iniP(5).pos(1:3) = [   0.0000d0,  10.0000d0,  17.3205d0 ]
    geom.iniP(6).pos(1:3) = [   0.0000d0, -10.0000d0,  17.3205d0 ]
    geom.iniP(7).pos(1:3) = [ -10.0000d0,   0.0000d0, -17.3205d0 ]
    geom.iniP(8).pos(1:3) = [  10.0000d0,   0.0000d0, -17.3205d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 5, 1, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 6, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 7, 1, 4 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 8, 3, 2 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 1, 5, 6, 4 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 3, 6, 5, 2 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [ 1, 7, 8, 2 ]
    geom.face(8).n_poi = 4; allocate(geom.face(8).poi(4)); geom.face(8).poi(1:4) = [ 3, 8, 7, 4 ]
end subroutine Exam_Johnson_Gyrobifastigium_J26

! ---------------------------------------------------------------------------------------

end module Exam_Johnson