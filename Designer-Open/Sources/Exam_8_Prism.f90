!
! ---------------------------------------------------------------------------------------
!
!                                Module for Exam_Prism
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2017/01/09
!                                                            Last  modified  : 2017/01/09
!
! ---------------------------------------------------------------------------------------
!
module Exam_Prism

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    ! http://www.georgehart.com/virtual-polyhedra/prisms-index.html

    ! Prism
    public Exam_Prism_Triangular_Prism          ! 64. V= 6,  E= 9,  F= 5
    public Exam_Prism_Pentagonal_Prism          ! 65. V= 7,  E=15,  F=10
    public Exam_Prism_Hexagonal_Prism           ! 28. V=26,  E=48,  F=24
    !public Exam_Prism_Heptagonal_Prism          ! 29. V=38,  E=60,  F=24
    !public Exam_Prism_Octagonal_Prism           ! 30. V=14,  E=36,  F=24
    !public Exam_Prism_Enneagonal_Prism          ! 31. V=26,  E=72,  F=48
    !public Exam_Prism_Decagonal_Prism           ! 32. V=32,  E=90,  F=60

    ! Antiprism
    !public Exam_Prism_Square_Antiprism          ! 33. V=32,  E=90,  F=60
    !public Exam_Prism_Pentagonal_Antiprism      ! 34. V=14,  E=36,  F=24
    !public Exam_Prism_Hexagonal_Antiprism       ! 35. V=8,   E=18,  F=12
    !public Exam_Prism_Heptagonal_Antiprism      ! 56. V=62,  E=180, F=120
    !public Exam_Prism_Octagonal_Antiprism       ! 57. V=62,  E=120, F=60
    !public Exam_Prism_Enneagonal_Antiprism      ! 58. V=92,  E=150, F=60
    !public Exam_Prism_Decagonal_Antiprism       ! 58. V=92,  E=150, F=60
contains

! ---------------------------------------------------------------------------------------

! Example of triangular prism
! Last updated on Tuesday 18 Octaober 2016 by Hyungmin
subroutine Exam_Prism_Triangular_Prism(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp
    double precision :: height

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "64_Triangular_Prism"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triangular Prism"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 8
        height               = 1.0d0
        !height               = 1.0d0 * 63.0d0 / 42.0d0      ! Long height
        height               = 1.0d0 * 42.0d0 / 63.0d0     ! Short height
    end if

    ! Set geometric type and view
    prob.color    = [150, 58, 228]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.03d0     ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   =-0.3d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 6
    geom.n_face = 5

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.5d0,  height/2.0d0, 0.0d0              ]
    geom.iniP(2).pos(1:3) = [  0.0d0,  height/2.0d0, dsqrt(3.0d0)/2.0d0 ]
    geom.iniP(3).pos(1:3) = [ -0.5d0,  height/2.0d0, 0.0d0              ]
    geom.iniP(4).pos(1:3) = [  0.5d0, -height/2.0d0, 0.0d0              ]
    geom.iniP(5).pos(1:3) = [  0.0d0, -height/2.0d0, dsqrt(3.0d0)/2.0d0 ]
    geom.iniP(6).pos(1:3) = [ -0.5d0, -height/2.0d0, 0.0d0              ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 3, 2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 4, 5, 6 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 1, 2, 5, 4 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 2, 3, 6, 5 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 3, 1, 4, 6 ]
end subroutine Exam_Prism_Triangular_Prism

! ---------------------------------------------------------------------------------------

! Example of pentagonal prism
! Last updated on Monday 9 January 2017 by Hyungmin
subroutine Exam_Prism_Pentagonal_Prism(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "65_Pentagonal_Prism"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triangular Prism"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.03d0     ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   =-0.3d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 10
    geom.n_face = 7

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.0d0,  0.0d0,  1.159953d0 ]
	geom.iniP( 2).pos(1:3) = [  1.013464d0,  0.0d0, 0.5642542d0 ]
	geom.iniP( 3).pos(1:3) = [  -0.3501431d0, 0.9510565d0, 0.5642542d0 ]
	geom.iniP( 4).pos(1:3) = [  -0.7715208d0, -0.6571639d0, 0.5642542d0 ]
	geom.iniP( 5).pos(1:3) = [  0.6633206d0, 0.9510565d0, -0.03144481d0 ]
	geom.iniP( 6).pos(1:3) = [  0.8682979d0, -0.6571639d0, -0.3996071d0 ]
	geom.iniP( 7).pos(1:3) = [  -1.121664d0, 0.2938926d0, -0.03144481d0 ]
	geom.iniP( 8).pos(1:3) = [  -0.2348831d0, -1.063314d0, -0.3996071d0 ]
	geom.iniP( 9).pos(1:3) = [  0.5181548d0, 0.2938926d0, -0.9953061d0 ]
	geom.iniP(10).pos(1:3) = [  -0.5850262d0, -0.112257d0, -0.9953061d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 5; allocate(geom.face(1).poi(5)); geom.face(1).poi(1:5) = [ 1,4,8,6,2 ]
    geom.face(2).n_poi = 5; allocate(geom.face(2).poi(5)); geom.face(2).poi(1:5) = [ 3,5,9,10,7 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 1,2,5,3 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1,3,7,4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 2,6,9,5 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 4,7,10,8 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [ 6,8,10,9 ]
end subroutine Exam_Prism_Pentagonal_Prism

! ---------------------------------------------------------------------------------------

! Example of hexagonal prism
! Last updated on Monday 9 January 2017 by Hyungmin
subroutine Exam_Prism_Hexagonal_Prism(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "65_Pentagonal_Prism"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triangular Prism"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.03d0     ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   =-0.3d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP =   12
    geom.n_face =    8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(      1).pos(1:3) = [     0.0000d0,    -0.0000d0,    22.3606d0 ]
    geom.iniP(      2).pos(1:3) = [    17.8886d0,    -0.0000d0,    13.4164d0 ]
    geom.iniP(      3).pos(1:3) = [    -4.4721d0,    17.3205d0,    13.4164d0 ]
    geom.iniP(      4).pos(1:3) = [   -15.6525d0,    -8.6603d0,    13.4164d0 ]
    geom.iniP(      5).pos(1:3) = [    13.4164d0,    17.3205d0,     4.4721d0 ]
    geom.iniP(      6).pos(1:3) = [    20.1246d0,    -8.6603d0,    -4.4721d0 ]
    geom.iniP(      7).pos(1:3) = [   -20.1246d0,     8.6603d0,     4.4721d0 ]
    geom.iniP(      8).pos(1:3) = [   -13.4164d0,   -17.3205d0,    -4.4721d0 ]
    geom.iniP(      9).pos(1:3) = [    15.6525d0,     8.6603d0,   -13.4164d0 ]
    geom.iniP(     10).pos(1:3) = [     4.4721d0,   -17.3205d0,   -13.4164d0 ]
    geom.iniP(     11).pos(1:3) = [   -17.8886d0,    -0.0000d0,   -13.4164d0 ]
    geom.iniP(     12).pos(1:3) = [     0.0000d0,    -0.0000d0,   -22.3606d0 ]

    ! Set point position vectors
    geom.face(      1).n_poi =       4; allocate(geom.face(      1).poi(      4)); geom.face(      1).poi(1:      4) = [      1,       2,       5,       3 ]
    geom.face(      2).n_poi =       4; allocate(geom.face(      2).poi(      4)); geom.face(      2).poi(1:      4) = [      1,       3,       7,       4 ]
    geom.face(      3).n_poi =       4; allocate(geom.face(      3).poi(      4)); geom.face(      3).poi(1:      4) = [      2,       6,       9,       5 ]
    geom.face(      4).n_poi =       4; allocate(geom.face(      4).poi(      4)); geom.face(      4).poi(1:      4) = [      4,       7,      11,       8 ]
    geom.face(      5).n_poi =       4; allocate(geom.face(      5).poi(      4)); geom.face(      5).poi(1:      4) = [      6,      10,      12,       9 ]
    geom.face(      6).n_poi =       4; allocate(geom.face(      6).poi(      4)); geom.face(      6).poi(1:      4) = [      8,      11,      12,      10 ]
    geom.face(      7).n_poi =       6; allocate(geom.face(      7).poi(      6)); geom.face(      7).poi(1:      6) = [      1,       4,       8,      10,       6,       2 ]
    geom.face(      8).n_poi =       6; allocate(geom.face(      8).poi(      6)); geom.face(      8).poi(1:      6) = [      3,       5,       9,      12,      11,       7 ]

end subroutine Exam_Prism_Hexagonal_Prism

! ---------------------------------------------------------------------------------------

end module Exam_Prism