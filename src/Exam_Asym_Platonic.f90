!
! ---------------------------------------------------------------------------------------
!
!                               Module - Exam_Asym_Platonic
!
!                                                                    Updated : 2017/03/18
!
! Comments: Asymmetric platonic solids
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Exam_Asym_Platonic

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Asym_Tetrahedron        ! 26. V=4,  E=6,  F=4
    public Exam_Asym_Cube               ! 27. V=8,  E=12, F=6
    public Exam_Asym_Octahedron         ! 28. V=6,  E=12, F=8
    public Exam_Asym_Dodecahedron       ! 29. V=20, E=30, F=12
    public Exam_Asym_Icosahedron        ! 30. V=12, E=30, F=20

contains

! ---------------------------------------------------------------------------------------

! Example of asymmetric tetrahedron
! Last updated on Mon 27 Mar 2017 by Hyungmin
subroutine Exam_Asym_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "26_Asym_Tetrahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Asym Tetrahedron"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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
end subroutine Exam_Asym_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric cube
! Last updated on Mon 27 Mar 2017 by Hyungmin
subroutine Exam_Asym_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "27_Asym_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Asym Cube"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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
end subroutine Exam_Asym_Cube

! ---------------------------------------------------------------------------------------

! Example of asymmetric octahedron
! Last updated on Mon 27 Mar 2017 by Hyungmin
subroutine Exam_Asym_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "28_Asym_Octahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Asym Octahedron"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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
end subroutine Exam_Asym_Octahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric dodecahedron
! Last updated on Mon 27 Mar 2017 by Hyungmin
subroutine Exam_Asym_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp
    
    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "29_Asym_Dodecahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Asym Dodecahedron"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! Allocate point and face structure
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
end subroutine Exam_Asym_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of asymmetric icosahedron
! Last updated on Mon 27 Mar 2017 by Hyungmin
subroutine Exam_Asym_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "30_Asym_Icosahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Asym Icosahedron"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [150, 58, 228], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! Allocate point and face structure
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
end subroutine Exam_Asym_Icosahedron

! ---------------------------------------------------------------------------------------

end module Exam_Asym_Platonic