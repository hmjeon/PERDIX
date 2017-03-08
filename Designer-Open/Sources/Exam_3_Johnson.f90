!
! ---------------------------------------------------------------------------------------
!
!                                Module for Exam_Johnson
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2016/04/07
!                                                            Last  modified  : 2016/07/14
!
! ---------------------------------------------------------------------------------------
!
module Exam_Johnson

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Johnson_Gyroelongated_Pentagonal_Pyramid_J11        ! 16. V=11, E=25, F=16
    public Exam_Johnson_Triangular_Bipyramid_J12                    ! 17. V=5,  E=9,  F=6
    public Exam_Johnson_Pentagonal_Bipyramid_J13                    ! 18. V=7,  E=15, F=10
    public Exam_Johnson_Gyroelongated_Square_Bipyramid_J17          ! 19. V=10, E=24, F=16
    public Exam_Johnson_Square_Gyrobicupola_J29                     ! 20. V=16, E=32, F=18
    public Exam_Johnson_Pentagonal_Orthocupolarotunda_J32           ! 21. V=25, E=50, F=27
    public Exam_Johnson_Pentagonal_Orthobirotunda_J34               ! 22. V=30, E=60, F=32
    public Exam_Johnson_Elongated_Pentagonal_Gyrobicupola_J39       ! 23. V=30, E=60, F=32
    public Exam_Johnson_Elongated_Pentagonal_Gyrobirotunda_J43      ! 24. V=40, E=80, F=42
    public Exam_Johnson_Gyroelongated_Square_Bicupola_J45           ! 25. V=24, E=56, F=34

    public Exam_Johnson_Pentagonal_Pyramid_J2                       ! 62. V=6,  E=10, F=6
    public Exam_Johnson_Elongated_Square_Bipyramid_J15              ! 63. V=10, E=20, F=12

contains

! ---------------------------------------------------------------------------------------

! Example of gyroelongated pentagonal pyramid(J11)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Pentagonal_Pyramid_J11(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "16_Gyroelon_Penta_Pyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyroelongated pentagonal pyramid"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.03d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "YZ"

    ! Allocate point and face structure
    geom.n_iniP = 11
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-13.6909d0,  -8.0685d0,  11.9047d0 ]
    geom.iniP( 2).pos(1:3) = [-12.6788d0,  11.7887d0,   9.7246d0 ]
    geom.iniP( 3).pos(1:3) = [ -9.5106d0, -16.4470d0,  -5.7704d0 ]
    geom.iniP( 4).pos(1:3) = [ -7.8705d0,  15.6810d0,  -9.2966d0 ]
    geom.iniP( 5).pos(1:3) = [ -5.9124d0,  -1.7701d0, -18.8732d0 ]
    geom.iniP( 6).pos(1:3) = [  2.5562d0,   1.8501d0,  18.0411d0 ]
    geom.iniP( 7).pos(1:3) = [  4.5143d0, -15.6010d0,   8.4645d0 ]
    geom.iniP( 8).pos(1:3) = [  6.1544d0,  16.5270d0,   4.9383d0 ]
    geom.iniP( 9).pos(1:3) = [  9.3226d0, -11.7087d0, -10.5567d0 ]
    geom.iniP(10).pos(1:3) = [ 10.3346d0,   8.1485d0, -12.7368d0 ]
    geom.iniP(11).pos(1:3) = [ 16.7811d0,  -0.4000d0,   4.1603d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  7,  3,  9 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  9,  3,  5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  9,  5, 10 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 10,  5,  4 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 10,  4,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  8,  4,  2 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  8,  2,  6 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  6,  2,  1 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  6,  1,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  7,  1,  3 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7,  9, 11 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  9, 10, 11 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 10,  8, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  8,  6, 11 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  6,  7, 11 ]
    geom.face(16).n_poi = 5; allocate(geom.face(16).poi(5)); geom.face(16).poi(1:5) = [  2,  4,  5,  3,  1 ]
end subroutine Exam_Johnson_Gyroelongated_Pentagonal_Pyramid_J11

! ---------------------------------------------------------------------------------------

! Example of triangular bipyramid(J12)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Triangular_Bipyramid_J12(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "17_Triang_Bipyramid"  //&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triangular bipyramid"

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
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 1.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "YZ"

    ! Allocate point and face structure
    geom.n_iniP =  5
    geom.n_face =  6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -8.3400d0,   3.3336d0,   7.2580d0 ]
    geom.iniP( 2).pos(1:3) = [ -2.5660d0,  -6.6664d0,  -9.0720d0 ]
    geom.iniP( 3).pos(1:3) = [ -2.5660d0,  13.3336d0,  -9.0720d0 ]
    geom.iniP( 4).pos(1:3) = [  2.5660d0, -13.3344d0,   9.0720d0 ]
    geom.iniP( 5).pos(1:3) = [ 10.9060d0,   3.3336d0,   1.8140d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 2, 4, 1 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 4, 5, 1 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 4, 2, 5 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 1, 3, 2 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 1, 5, 3 ]
    geom.face(6).n_poi = 3; allocate(geom.face(6).poi(3)); geom.face(6).poi(1:3) = [ 3, 5, 2 ]
end subroutine Exam_Johnson_Triangular_Bipyramid_J12

! ---------------------------------------------------------------------------------------

! Example of pentagonal bipyramid(J13)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Bipyramid_J13(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "18_Penta_Bipyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal bipyramid"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        !para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.95d0     ! Cylindrical model
    prob.move_x   = -2.5d0     ! Cylindrical model
    prob.move_y   = 0.5d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP =  7
    geom.n_face = 10

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [-15.8832d0,   6.0527d0,  -0.7538d0 ]
    geom.iniP(2).pos(1:3) = [ -9.8887d0,  -9.9765d0,   9.6010d0 ]
    geom.iniP(3).pos(1:3) = [ -1.9321d0,  -6.1082d0,  -8.3383d0 ]
    geom.iniP(4).pos(1:3) = [  0.0720d0,  13.7173d0, -10.0665d0 ]
    geom.iniP(5).pos(1:3) = [  1.9321d0,   6.1087d0,   8.3389d0 ]
    geom.iniP(6).pos(1:3) = [  9.7727d0, -12.2186d0,   6.6868d0 ]
    geom.iniP(7).pos(1:3) = [ 15.9272d0,   2.4245d0,  -5.4681d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 4, 3, 1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 3, 2, 1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 3, 6, 2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 1, 5, 4 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 2, 5 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 5, 2, 6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 3, 4, 7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 5, 7 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 6, 3, 7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 5, 6, 7 ]
end subroutine Exam_Johnson_Pentagonal_Bipyramid_J13

! ---------------------------------------------------------------------------------------

! Example of gyroelongated square bipyramid(J17)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Square_Bipyramid_J17(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "19_Gyroelong_Square_Bipyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyroelongated square bipyramid"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
    end if

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.9d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 10
    geom.n_face = 16

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-13.8671d0,   8.6628d0,   1.8381d0 ]
    geom.iniP( 2).pos(1:3) = [-12.0490d0, -10.0885d0,  -4.8763d0 ]
    geom.iniP( 3).pos(1:3) = [ -6.7767d0,  -5.6322d0,  13.8969d0 ]
    geom.iniP( 4).pos(1:3) = [ -3.9585d0,   5.0425d0, -15.1550d0 ]
    geom.iniP( 5).pos(1:3) = [ -0.6182d0,  21.9737d0,  -5.0423d0 ]
    geom.iniP( 6).pos(1:3) = [  0.6178d0, -21.9733d0,   5.0423d0 ]
    geom.iniP( 7).pos(1:3) = [  3.4980d0,  11.3450d0,  11.3928d0 ]
    geom.iniP( 8).pos(1:3) = [  7.2363d0, -10.7545d0, -10.1347d0 ]
    geom.iniP( 9).pos(1:3) = [ 12.5107d0,  -6.2982d0,   8.6366d0 ]
    geom.iniP(10).pos(1:3) = [ 13.4067d0,   7.7227d0,  -5.5984d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  7,  9, 10 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 10,  9,  8 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 10,  8,  4 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  4,  8,  2 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  4,  2,  1 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  1,  2,  3 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  1,  3,  7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  7,  3,  9 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  7, 10,  5 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 10,  4,  5 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4,  1,  5 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  1,  7,  5 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  8,  9,  6 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  2,  8,  6 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  3,  2,  6 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  9,  3,  6 ]
end subroutine Exam_Johnson_Gyroelongated_Square_Bipyramid_J17

! ---------------------------------------------------------------------------------------

! Example of square gyrobicupola(J29)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Square_Gyrobicupola_J29(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "20_Square_Gyrobicupola"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Square gyrobicupola"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 16
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-25.7139d0,  -1.8023d0,  -4.3023d0 ]
    geom.iniP( 2).pos(1:3) = [-20.0835d0,  16.7031d0,   0.7821d0 ]
    geom.iniP( 3).pos(1:3) = [-16.2812d0, -19.2535d0,  -6.8665d0 ]
    geom.iniP( 4).pos(1:3) = [-14.3630d0,  -9.2868d0,  10.3688d0 ]
    geom.iniP( 5).pos(1:3) = [-11.3508d0,   7.4844d0, -14.6711d0 ]
    geom.iniP( 6).pos(1:3) = [ -8.7326d0,   9.2205d0,  15.4531d0 ]
    geom.iniP( 7).pos(1:3) = [ -2.6882d0,  25.4257d0,   5.4084d0 ]
    geom.iniP( 8).pos(1:3) = [ -1.9181d0,  -9.9668d0, -17.2353d0 ]
    geom.iniP( 9).pos(1:3) = [  2.6882d0, -25.4260d0,  -5.4084d0 ]
    geom.iniP(10).pos(1:3) = [  4.6063d0, -15.4592d0,  11.8269d0 ]
    geom.iniP(11).pos(1:3) = [  6.0444d0,  16.2051d0, -10.0447d0 ]
    geom.iniP(12).pos(1:3) = [ 10.2367d0,   3.0481d0,  16.9112d0 ]
    geom.iniP(13).pos(1:3) = [ 15.4771d0,  -1.2442d0, -12.6089d0 ]
    geom.iniP(14).pos(1:3) = [ 16.2812d0,  19.2533d0,   6.8665d0 ]
    geom.iniP(15).pos(1:3) = [ 20.0835d0, -16.7033d0,  -0.7821d0 ]
    geom.iniP(16).pos(1:3) = [ 25.7139d0,   1.8020d0,   4.3023d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  5,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 11,  7, 14 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 13, 16, 15 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8,  9,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  4,  1,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 10,  9, 15 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 12, 16, 14 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  6,  7,  2 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  8,  5, 11, 13 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  5,  8,  3,  1 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 11,  5,  2,  7 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 13, 11, 14, 16 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  8, 13, 15,  9 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  6,  4, 10, 12 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  4,  6,  2,  1 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 10,  4,  3,  9 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 12, 10, 15, 16 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  6, 12, 14,  7 ]
end subroutine Exam_Johnson_Square_Gyrobicupola_J29

! ---------------------------------------------------------------------------------------

! Example of pentagonal orthocupolarotunda(J32)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Orthocupolarotunda_J32(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "21_Penta_Orthocupolarotunda"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal orthocupolarotunda"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.96d0     ! Cylindrical model
    prob.move_x   = 2.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 25
    geom.n_face = 27

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-32.0843d0,   7.9926d0,  -0.6560d0 ]
    geom.iniP( 2).pos(1:3) = [-28.0900d0,   0.4821d0,  17.4473d0 ]
    geom.iniP( 3).pos(1:3) = [-24.9218d0,  10.1988d0, -19.1992d0 ]
    geom.iniP( 4).pos(1:3) = [-21.4836d0,  -6.7184d0,  -9.0985d0 ]
    geom.iniP( 5).pos(1:3) = [-20.0255d0,  17.9073d0,  11.8489d0 ]
    geom.iniP( 6).pos(1:3) = [-17.4913d0, -14.2309d0,   9.0047d0 ]
    geom.iniP( 7).pos(1:3) = [-14.4651d0,  -9.4666d0,  28.1920d0 ]
    geom.iniP( 8).pos(1:3) = [ -9.3367d0,   6.2565d0, -31.1021d0 ]
    geom.iniP( 9).pos(1:3) = [ -8.4347d0,  21.4776d0, -18.1572d0 ]
    geom.iniP(10).pos(1:3) = [ -5.9005d0, -10.6607d0, -21.0014d0 ]
    geom.iniP(11).pos(1:3) = [ -5.4105d0,  26.2419d0,   1.0322d0 ]
    geom.iniP(12).pos(1:3) = [ -1.4162d0,  18.7314d0,  19.1334d0 ]
    geom.iniP(13).pos(1:3) = [  0.5600d0, -22.8135d0,   8.2887d0 ]
    geom.iniP(14).pos(1:3) = [  2.0201d0,   1.8122d0,  29.2361d0 ]
    geom.iniP(15).pos(1:3) = [  3.5862d0, -18.0492d0,  27.4780d0 ]
    geom.iniP(16).pos(1:3) = [  7.7245d0, -20.6073d0, -10.2546d0 ]
    geom.iniP(17).pos(1:3) = [  8.7145d0,  -2.3261d0, -31.8161d0 ]
    geom.iniP(18).pos(1:3) = [ 10.1726d0,  22.2996d0, -10.8707d0 ]
    geom.iniP(19).pos(1:3) = [ 16.6351d0,  10.1468d0,  18.4194d0 ]
    geom.iniP(20).pos(1:3) = [ 19.1692d0, -21.9914d0,  15.5752d0 ]
    geom.iniP(21).pos(1:3) = [ 20.7734d0,   7.5886d0, -19.3133d0 ]
    geom.iniP(22).pos(1:3) = [ 22.3395d0, -12.2748d0, -21.0714d0 ]
    geom.iniP(23).pos(1:3) = [ 23.7976d0,  12.3529d0,  -0.1239d0 ]
    geom.iniP(24).pos(1:3) = [ 26.3337d0, -19.7853d0,  -2.9681d0 ]
    geom.iniP(25).pos(1:3) = [ 27.2338d0,  -4.5642d0,   9.9768d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  4,  1,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 10,  8, 17 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 16, 22, 24 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 13, 20, 15 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  6,  7,  2 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 19, 25, 23 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 25, 20, 24 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 23, 21, 18 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 21, 22, 17 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 18,  9, 11 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  9,  8,  3 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 11,  5, 12 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  5,  1,  2 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 12, 14, 19 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 14,  7, 15 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  4,  6,  2,  1 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 10,  4,  3,  8 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 16, 10, 17, 22 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 13, 16, 24, 20 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [  6, 13, 15,  7 ]
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [  6,  4, 10, 16, 13 ]
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 19, 23, 18, 11, 12 ]
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 19, 14, 15, 20, 25 ]
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 23, 25, 24, 22, 21 ]
    geom.face(25).n_poi = 5; allocate(geom.face(25).poi(5)); geom.face(25).poi(1:5) = [ 18, 21, 17,  8,  9 ]
    geom.face(26).n_poi = 5; allocate(geom.face(26).poi(5)); geom.face(26).poi(1:5) = [ 11,  9,  3,  1,  5 ]
    geom.face(27).n_poi = 5; allocate(geom.face(27).poi(5)); geom.face(27).poi(1:5) = [ 12,  5,  2,  7, 14 ]
end subroutine Exam_Johnson_Pentagonal_Orthocupolarotunda_J32

! ---------------------------------------------------------------------------------------

! Example of pentagonal orthobirotunda(J34)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Orthobirotunda_J34(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "22_Penta_Orthobirotunda"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal orthobirotunda"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.99d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 30
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-31.5870d0,   0.6861d0,   7.0107d0 ]
    geom.iniP( 2).pos(1:3) = [-29.0827d0, -10.9010d0,  -9.0989d0 ]
    geom.iniP( 3).pos(1:3) = [-25.9184d0,  19.2559d0,   2.2082d0 ]
    geom.iniP( 4).pos(1:3) = [-25.1924d0, -18.1476d0,   9.1349d0 ]
    geom.iniP( 5).pos(1:3) = [-22.0281d0,  12.0112d0,  20.4419d0 ]
    geom.iniP( 6).pos(1:3) = [-21.8641d0,   0.5081d0, -23.8562d0 ]
    geom.iniP( 7).pos(1:3) = [-19.9079d0,  19.1459d0, -16.8676d0 ]
    geom.iniP( 8).pos(1:3) = [-16.9576d0, -26.6684d0,  -6.9787d0 ]
    geom.iniP( 9).pos(1:3) = [-11.6811d0, -18.4597d0,  23.8782d0 ]
    geom.iniP(10).pos(1:3) = [-10.3470d0,  30.4709d0,  -3.4363d0 ]
    geom.iniP(11).pos(1:3) = [ -9.7249d0,   0.1781d0,  30.8669d0 ]
    geom.iniP(12).pos(1:3) = [ -6.2926d0,  11.7232d0, -29.5028d0 ]
    geom.iniP(13).pos(1:3) = [ -5.2765d0,  -8.2067d0, -30.8589d0 ]
    geom.iniP(14).pos(1:3) = [ -4.0544d0,  18.7498d0,  26.0664d0 ]
    geom.iniP(15).pos(1:3) = [ -2.2442d0, -25.0023d0, -20.4259d0 ]
    geom.iniP(16).pos(1:3) = [  1.6442d0, -32.2470d0,  -2.1942d0 ]
    geom.iniP(17).pos(1:3) = [  3.1643d0,  30.1589d0,  11.3091d0 ]
    geom.iniP(18).pos(1:3) = [  4.9065d0, -27.1745d0,  16.8776d0 ]
    geom.iniP(19).pos(1:3) = [  6.2926d0, -11.7230d0,  29.5028d0 ]
    geom.iniP(20).pos(1:3) = [  9.1749d0,  30.0489d0,  -7.7687d0 ]
    geom.iniP(21).pos(1:3) = [ 11.6811d0,  18.4598d0, -23.8782d0 ]
    geom.iniP(22).pos(1:3) = [ 13.3252d0, -13.7872d0, -26.0724d0 ]
    geom.iniP(23).pos(1:3) = [ 15.4675d0,  18.3258d0,  21.7340d0 ]
    geom.iniP(24).pos(1:3) = [ 19.6198d0, -25.5103d0,   3.4303d0 ]
    geom.iniP(25).pos(1:3) = [ 21.8641d0,  -0.5080d0,  23.8562d0 ]
    geom.iniP(26).pos(1:3) = [ 23.8062d0,   2.6943d0, -21.7580d0 ]
    geom.iniP(27).pos(1:3) = [ 25.1924d0,  18.1478d0,  -9.1349d0 ]
    geom.iniP(28).pos(1:3) = [ 26.8385d0, -14.1013d0, -11.3271d0 ]
    geom.iniP(29).pos(1:3) = [ 29.0827d0,  10.9011d0,   9.0989d0 ]
    geom.iniP(30).pos(1:3) = [ 30.0988d0,  -9.0288d0,   7.7427d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 16,  8, 15 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  8,  4,  2 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 15, 13, 22 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 13,  6, 12 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 22, 26, 28 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 26, 21, 27 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 28, 30, 24 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 30, 29, 25 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 24, 18, 16 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 18, 19,  9 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  5, 11, 14 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 11,  9, 19 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 14, 23, 17 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 23, 25, 29 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 17, 20, 10 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 20, 27, 21 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 10,  7,  3 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  7, 12,  6 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  3,  1,  5 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  1,  2,  4 ]
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [ 16, 15, 22, 28, 24 ]
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 16, 18,  9,  4,  8 ]
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 15,  8,  2,  6, 13 ]
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 22, 13, 12, 21, 26 ]
    geom.face(25).n_poi = 5; allocate(geom.face(25).poi(5)); geom.face(25).poi(1:5) = [ 28, 26, 27, 29, 30 ]
    geom.face(26).n_poi = 5; allocate(geom.face(26).poi(5)); geom.face(26).poi(1:5) = [ 24, 30, 25, 19, 18 ]
    geom.face(27).n_poi = 5; allocate(geom.face(27).poi(5)); geom.face(27).poi(1:5) = [  5, 14, 17, 10,  3 ]
    geom.face(28).n_poi = 5; allocate(geom.face(28).poi(5)); geom.face(28).poi(1:5) = [  5,  1,  4,  9, 11 ]
    geom.face(29).n_poi = 5; allocate(geom.face(29).poi(5)); geom.face(29).poi(1:5) = [ 14, 11, 19, 25, 23 ]
    geom.face(30).n_poi = 5; allocate(geom.face(30).poi(5)); geom.face(30).poi(1:5) = [ 17, 23, 29, 27, 20 ]
    geom.face(31).n_poi = 5; allocate(geom.face(31).poi(5)); geom.face(31).poi(1:5) = [ 10, 20, 21, 12,  7 ]
    geom.face(32).n_poi = 5; allocate(geom.face(32).poi(5)); geom.face(32).poi(1:5) = [  3,  7,  6,  2,  1 ]
end subroutine Exam_Johnson_Pentagonal_Orthobirotunda_J34

! ---------------------------------------------------------------------------------------

! Example of elongated pentagonal gyrobicupola(J39)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Gyrobicupola_J39(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "23_Elong_Penta_Gyrobicupola"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elongated pentagonal gyrobicupola"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.97d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 30
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-31.8656d0,   6.8743d0,  -9.1965d0 ]
    geom.iniP( 2).pos(1:3) = [-31.3436d0,   6.9564d0,  10.7965d0 ]
    geom.iniP( 3).pos(1:3) = [-29.8915d0, -13.0287d0,  -9.1685d0 ]
    geom.iniP( 4).pos(1:3) = [-29.3695d0, -12.9467d0,  10.8266d0 ]
    geom.iniP( 5).pos(1:3) = [-21.7691d0,  24.1372d0,  -9.5325d0 ]
    geom.iniP( 6).pos(1:3) = [-21.2451d0,  24.2172d0,  10.4625d0 ]
    geom.iniP( 7).pos(1:3) = [-17.4629d0,  -1.7621d0, -20.0590d0 ]
    geom.iniP( 8).pos(1:3) = [-16.6008d0, -27.9714d0,  -9.4545d0 ]
    geom.iniP( 9).pos(1:3) = [-16.0768d0, -27.8894d0,  10.5385d0 ]
    geom.iniP(10).pos(1:3) = [-14.1427d0,   8.6764d0,  20.8571d0 ]
    geom.iniP(11).pos(1:3) = [-12.1686d0, -11.2266d0,  20.8871d0 ]
    geom.iniP(12).pos(1:3) = [ -7.3644d0,  15.5008d0, -20.3930d0 ]
    geom.iniP(13).pos(1:3) = [ -4.1702d0, -16.7049d0, -20.3450d0 ]
    geom.iniP(14).pos(1:3) = [ -3.4562d0,  32.1636d0, -10.0445d0 ]
    geom.iniP(15).pos(1:3) = [ -2.9321d0,  32.2456d0,   9.9505d0 ]
    geom.iniP(16).pos(1:3) = [  2.9321d0, -32.2456d0,  -9.9505d0 ]
    geom.iniP(17).pos(1:3) = [  3.4562d0, -32.1636d0,  10.0445d0 ]
    geom.iniP(18).pos(1:3) = [  4.1702d0,  16.7049d0,  20.3450d0 ]
    geom.iniP(19).pos(1:3) = [  7.3644d0, -15.5008d0,  20.3930d0 ]
    geom.iniP(20).pos(1:3) = [ 12.1686d0,  11.2266d0, -20.8871d0 ]
    geom.iniP(21).pos(1:3) = [ 14.1427d0,  -8.6764d0, -20.8571d0 ]
    geom.iniP(22).pos(1:3) = [ 16.0768d0,  27.8894d0, -10.5385d0 ]
    geom.iniP(23).pos(1:3) = [ 16.6008d0,  27.9714d0,   9.4545d0 ]
    geom.iniP(24).pos(1:3) = [ 17.4629d0,   1.7621d0,  20.0590d0 ]
    geom.iniP(25).pos(1:3) = [ 21.2451d0, -24.2172d0, -10.4625d0 ]
    geom.iniP(26).pos(1:3) = [ 21.7691d0, -24.1372d0,   9.5325d0 ]
    geom.iniP(27).pos(1:3) = [ 29.3695d0,  12.9467d0, -10.8266d0 ]
    geom.iniP(28).pos(1:3) = [ 29.8915d0,  13.0287d0,   9.1685d0 ]
    geom.iniP(29).pos(1:3) = [ 31.3436d0,  -6.9564d0, -10.7965d0 ]
    geom.iniP(30).pos(1:3) = [ 31.8656d0,  -6.8743d0,   9.1965d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 12,  5, 14 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 20, 22, 27 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 21, 29, 25 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 13, 16,  8 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  7,  3,  1 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 10,  6,  2 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 11,  4,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 19, 17, 26 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 24, 30, 28 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 18, 23, 15 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 12,  7,  1,  5 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 20, 12, 14, 22 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 21, 20, 27, 29 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 13, 21, 25, 16 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  7, 13,  8,  3 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 10, 18, 15,  6 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 11, 10,  2,  4 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 19, 11,  9, 17 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 24, 19, 26, 30 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 18, 24, 28, 23 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 25, 26, 17, 16 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 16, 17,  9,  8 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [  8,  9,  4,  3 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [  3,  4,  2,  1 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [  1,  2,  6,  5 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [  5,  6, 15, 14 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 14, 15, 23, 22 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 22, 23, 28, 27 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 27, 28, 30, 29 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 29, 30, 26, 25 ]
    geom.face(31).n_poi = 5; allocate(geom.face(31).poi(5)); geom.face(31).poi(1:5) = [  7, 12, 20, 21, 13 ]
    geom.face(32).n_poi = 5; allocate(geom.face(32).poi(5)); geom.face(32).poi(1:5) = [ 18, 10, 11, 19, 24 ]
end subroutine Exam_Johnson_Elongated_Pentagonal_Gyrobicupola_J39

! ---------------------------------------------------------------------------------------

! Example of elongated pentagonal gyrobirotunda(J43)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Gyrobirotunda_J43(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "24_Elong_Penta_Gyrobirotunda"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elongated pentagonal gyrobirotunda"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.9d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 40
    geom.n_face = 42

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-40.7127d0,  -6.3204d0,  -0.6760d0 ]
    geom.iniP( 2).pos(1:3) = [-37.5985d0,   6.8305d0, -15.4190d0 ]
    geom.iniP( 3).pos(1:3) = [-36.2404d0,  12.6869d0,   3.6582d0 ]
    geom.iniP( 4).pos(1:3) = [-32.9802d0, -23.8816d0,  -6.3204d0 ]
    geom.iniP( 5).pos(1:3) = [-31.6241d0, -18.0252d0,  12.7569d0 ]
    geom.iniP( 6).pos(1:3) = [-27.9399d0,  -2.6022d0, -30.1760d0 ]
    geom.iniP( 7).pos(1:3) = [-25.0857d0, -21.5835d0, -24.5537d0 ]
    geom.iniP( 8).pos(1:3) = [-24.3856d0,  12.7289d0,  19.7673d0 ]
    geom.iniP( 9).pos(1:3) = [-23.4656d0,  16.4051d0, -25.8437d0 ]
    geom.iniP(10).pos(1:3) = [-21.5315d0,  -6.2524d0,  25.3897d0 ]
    geom.iniP(11).pos(1:3) = [-21.2694d0,  25.8797d0,   5.0223d0 ]
    geom.iniP(12).pos(1:3) = [-15.9951d0, -33.2882d0, -11.1207d0 ]
    geom.iniP(13).pos(1:3) = [-13.7989d0, -23.8136d0,  19.7453d0 ]
    geom.iniP(14).pos(1:3) = [-13.3749d0,  28.1779d0, -13.2109d0 ]
    geom.iniP(15).pos(1:3) = [-10.9527d0, -12.0088d0, -34.9784d0 ]
    geom.iniP(16).pos(1:3) = [ -6.8705d0,  17.7972d0,  27.9899d0 ]
    geom.iniP(17).pos(1:3) = [ -4.1403d0, -33.2462d0,   4.9883d0 ]
    geom.iniP(18).pos(1:3) = [ -4.0183d0,  -1.1841d0,  33.6123d0 ]
    geom.iniP(19).pos(1:3) = [ -3.7543d0,  30.9481d0,  13.2449d0 ]
    geom.iniP(20).pos(1:3) = [ -3.7163d0,  18.7453d0, -27.9679d0 ]
    geom.iniP(21).pos(1:3) = [  3.7163d0, -18.7453d0,  27.9679d0 ]
    geom.iniP(22).pos(1:3) = [  3.7543d0, -30.9481d0, -13.2449d0 ]
    geom.iniP(23).pos(1:3) = [  4.0183d0,   1.1841d0, -33.6123d0 ]
    geom.iniP(24).pos(1:3) = [  4.1403d0,  33.2462d0,  -4.9883d0 ]
    geom.iniP(25).pos(1:3) = [  6.8705d0, -17.7972d0, -27.9899d0 ]
    geom.iniP(26).pos(1:3) = [ 10.9527d0,  12.0088d0,  34.9784d0 ]
    geom.iniP(27).pos(1:3) = [ 13.3749d0, -28.1779d0,  13.2109d0 ]
    geom.iniP(28).pos(1:3) = [ 13.7989d0,  23.8136d0, -19.7453d0 ]
    geom.iniP(29).pos(1:3) = [ 15.9951d0,  33.2882d0,  11.1207d0 ]
    geom.iniP(30).pos(1:3) = [ 21.2694d0, -25.8797d0,  -5.0223d0 ]
    geom.iniP(31).pos(1:3) = [ 21.5315d0,   6.2524d0, -25.3897d0 ]
    geom.iniP(32).pos(1:3) = [ 23.4656d0, -16.4051d0,  25.8437d0 ]
    geom.iniP(33).pos(1:3) = [ 24.3856d0, -12.7289d0, -19.7673d0 ]
    geom.iniP(34).pos(1:3) = [ 25.0857d0,  21.5835d0,  24.5537d0 ]
    geom.iniP(35).pos(1:3) = [ 27.9399d0,   2.6022d0,  30.1760d0 ]
    geom.iniP(36).pos(1:3) = [ 31.6241d0,  18.0252d0, -12.7569d0 ]
    geom.iniP(37).pos(1:3) = [ 32.9802d0,  23.8816d0,   6.3204d0 ]
    geom.iniP(38).pos(1:3) = [ 36.2404d0, -12.6869d0,  -3.6582d0 ]
    geom.iniP(39).pos(1:3) = [ 37.5985d0,  -6.8305d0,  15.4190d0 ]
    geom.iniP(40).pos(1:3) = [ 40.7127d0,   6.3204d0,   0.6760d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  3,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  8, 11 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  9,  6 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  9, 14, 20 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  6, 15,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 15, 23, 25 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  7, 12,  4 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 12, 22, 17 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4,  5,  1 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  5, 13, 10 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 34, 26, 35 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 26, 16, 18 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 35, 32, 39 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 32, 21, 27 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 39, 38, 40 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 38, 30, 33 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 40, 36, 37 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 36, 31, 28 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 37, 29, 34 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 29, 24, 19 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 25, 33, 30, 22 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 22, 30, 27, 17 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 17, 27, 21, 13 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 13, 21, 18, 10 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 10, 18, 16,  8 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [  8, 16, 19, 11 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 11, 19, 24, 14 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 14, 24, 28, 20 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 20, 28, 31, 23 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 23, 31, 33, 25 ]
    geom.face(31).n_poi = 5; allocate(geom.face(31).poi(5)); geom.face(31).poi(1:5) = [  1,  2,  6,  7,  4 ]
    geom.face(32).n_poi = 5; allocate(geom.face(32).poi(5)); geom.face(32).poi(1:5) = [  1,  5, 10,  8,  3 ]
    geom.face(33).n_poi = 5; allocate(geom.face(33).poi(5)); geom.face(33).poi(1:5) = [  2,  3, 11, 14,  9 ]
    geom.face(34).n_poi = 5; allocate(geom.face(34).poi(5)); geom.face(34).poi(1:5) = [  6,  9, 20, 23, 15 ]
    geom.face(35).n_poi = 5; allocate(geom.face(35).poi(5)); geom.face(35).poi(1:5) = [  7, 15, 25, 22, 12 ]
    geom.face(36).n_poi = 5; allocate(geom.face(36).poi(5)); geom.face(36).poi(1:5) = [  4, 12, 17, 13,  5 ]
    geom.face(37).n_poi = 5; allocate(geom.face(37).poi(5)); geom.face(37).poi(1:5) = [ 34, 35, 39, 40, 37 ]
    geom.face(38).n_poi = 5; allocate(geom.face(38).poi(5)); geom.face(38).poi(1:5) = [ 34, 29, 19, 16, 26 ]
    geom.face(39).n_poi = 5; allocate(geom.face(39).poi(5)); geom.face(39).poi(1:5) = [ 35, 26, 18, 21, 32 ]
    geom.face(40).n_poi = 5; allocate(geom.face(40).poi(5)); geom.face(40).poi(1:5) = [ 39, 32, 27, 30, 38 ]
    geom.face(41).n_poi = 5; allocate(geom.face(41).poi(5)); geom.face(41).poi(1:5) = [ 40, 38, 33, 31, 36 ]
    geom.face(42).n_poi = 5; allocate(geom.face(42).poi(5)); geom.face(42).poi(1:5) = [ 37, 36, 28, 24, 29 ]
end subroutine Exam_Johnson_Elongated_Pentagonal_Gyrobirotunda_J43

! ---------------------------------------------------------------------------------------

! Example of gyroelongated square bicupola(J45)
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Square_Bicupola_J45(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "25_Gyroelong_Square_Bicupola"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Gyroelongated square bicupola"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 1.0d0      ! Cylindrical model
    prob.move_y   =-1.2d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 24
    geom.n_face = 34

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [-26.8535d0,  -5.8744d0,   1.1680d0 ]
    geom.iniP( 2).pos(1:3) = [-22.7752d0,  11.3726d0,  10.4366d0 ]
    geom.iniP( 3).pos(1:3) = [-21.1811d0,   2.1460d0, -16.2551d0 ]
    geom.iniP( 4).pos(1:3) = [-18.5810d0,  -6.0205d0,  19.3771d0 ]
    geom.iniP( 5).pos(1:3) = [-17.5289d0, -17.1511d0, -12.4689d0 ]
    geom.iniP( 6).pos(1:3) = [-17.1029d0,  19.3931d0,  -6.9865d0 ]
    geom.iniP( 7).pos(1:3) = [-15.7408d0, -21.4574d0,   6.9803d0 ]
    geom.iniP( 8).pos(1:3) = [ -8.2063d0,  10.4146d0,  24.1074d0 ]
    geom.iniP( 9).pos(1:3) = [ -7.6843d0,  24.4894d0,   9.9065d0 ]
    geom.iniP(10).pos(1:3) = [ -3.9161d0,   3.4461d0, -26.2717d0 ]
    geom.iniP(11).pos(1:3) = [ -1.8719d0, -16.9371d0,  20.6652d0 ]
    geom.iniP(12).pos(1:3) = [ -1.3519d0, -26.8557d0,  -5.8204d0 ]
    geom.iniP(13).pos(1:3) = [ -0.2658d0, -15.8511d0, -22.4855d0 ]
    geom.iniP(14).pos(1:3) = [  0.1602d0,  20.6952d0, -17.0031d0 ]
    geom.iniP(15).pos(1:3) = [  8.5027d0,  -0.5041d0,  25.3955d0 ]
    geom.iniP(16).pos(1:3) = [  9.3027d0,  18.2170d0,  18.3991d0 ]
    geom.iniP(17).pos(1:3) = [  9.5788d0,  25.7915d0,  -0.1101d0 ]
    geom.iniP(18).pos(1:3) = [ 12.5169d0, -22.3375d0,   7.8644d0 ]
    geom.iniP(19).pos(1:3) = [ 14.8251d0,  -2.7343d0, -23.0155d0 ]
    geom.iniP(20).pos(1:3) = [ 16.1572d0, -19.0533d0, -11.5288d0 ]
    geom.iniP(21).pos(1:3) = [ 18.9033d0,  14.5148d0, -13.7469d0 ]
    geom.iniP(22).pos(1:3) = [ 22.8916d0,  -5.9024d0,  12.5927d0 ]
    geom.iniP(23).pos(1:3) = [ 23.6916d0,  12.8187d0,   5.5983d0 ]
    geom.iniP(24).pos(1:3) = [ 26.5298d0,  -2.6202d0,  -6.7985d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 14, 17, 21 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 10, 19, 13 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  3,  5,  1 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  2,  9 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 15, 16,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 11,  4,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 18, 12, 20 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 22, 24, 23 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 13, 12,  5 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  5, 12,  7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  5,  7,  1 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  1,  7,  4 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  1,  4,  2 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  2,  4,  8 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  2,  8,  9 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  9,  8, 16 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  9, 16, 17 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 17, 16, 23 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 17, 23, 21 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 21, 23, 24 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 21, 24, 19 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 19, 24, 20 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 19, 20, 13 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 13, 20, 12 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [  6, 14, 10,  3 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 14,  6,  9, 17 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 10, 14, 21, 19 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [  3, 10, 13,  5 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [  6,  3,  1,  2 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 22, 15, 11, 18 ]
    geom.face(31).n_poi = 4; allocate(geom.face(31).poi(4)); geom.face(31).poi(1:4) = [ 15, 22, 23, 16 ]
    geom.face(32).n_poi = 4; allocate(geom.face(32).poi(4)); geom.face(32).poi(1:4) = [ 11, 15,  8,  4 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [ 18, 11,  7, 12 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [ 22, 18, 20, 24 ]
end subroutine Exam_Johnson_Gyroelongated_Square_Bicupola_J45

! ---------------------------------------------------------------------------------------

! Example of Pentagonal pyramid(J2)
! Last updated on Wednesday 12 Octaober 2016 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Pyramid_J2(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "62_Pentagonal_Pyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal Pyramid"

    ! Set geometric type and view
    prob.color    = [77, 175, 74]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   =-3.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 6
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -13.9168d0,  -1.6021d0,  9.8125d0 ]
    geom.iniP(2).pos(1:3) = [  -5.2781d0,  15.6353d0,  4.4980d0 ]
    geom.iniP(3).pos(1:3) = [  -4.2660d0,  -0.2200d0, -7.6510d0 ]
    geom.iniP(4).pos(1:3) = [  -2.1458d0, -16.5654d0,  3.6820d0 ]
    geom.iniP(5).pos(1:3) = [  11.8353d0,  11.3269d0, -4.9187d0 ]
    geom.iniP(6).pos(1:3) = [  13.7715d0,  -8.5747d0, -5.4228d0 ]

    ! Set point position vectors
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 4, 1, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 6, 4, 3 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 5, 6, 3 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 5, 3 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 1, 2, 3 ]
    geom.face(6).n_poi = 5; allocate(geom.face(6).poi(5)); geom.face(6).poi(1:5) = [ 1, 4, 6, 5, 2 ]
end subroutine Exam_Johnson_Pentagonal_Pyramid_J2

! ---------------------------------------------------------------------------------------

! Example of Elongated square bipyramid(J15)
! Last updated on Tuesday 18 Octaober 2016 by Hyungmin
subroutine Exam_Johnson_Elongated_Square_Bipyramid_J15(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "63_Elongated_Square_Bipyramid"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Elongated Square Bipyramid"

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
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 10
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -12.6491d0,   6.3246d0, -10.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -12.6491d0,   6.3246d0,  10.0000d0 ]
    geom.iniP( 3).pos(1:3) = [  -7.6344d0,  22.9032d0,  -0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [  -6.3246d0, -12.6491d0, -10.0000d0 ]
    geom.iniP( 5).pos(1:3) = [  -6.3246d0, -12.6491d0,  10.0000d0 ]
    geom.iniP( 6).pos(1:3) = [   6.3246d0,  12.6491d0, -10.0000d0 ]
    geom.iniP( 7).pos(1:3) = [   6.3246d0,  12.6491d0,  10.0000d0 ]
    geom.iniP( 8).pos(1:3) = [   7.6344d0, -22.9032d0,  -0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [  12.6491d0,  -6.3246d0, -10.0000d0 ]
    geom.iniP(10).pos(1:3) = [  12.6491d0,  -6.3246d0,  10.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  9, 10,  8 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  7,  6,  3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  4,  9,  8 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  1,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  5,  4,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  1,  2,  3 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 10,  5,  8 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  2,  7,  3 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 10,  9,  6, 7 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  9,  4,  1, 6 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  4,  5,  2, 1 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  5, 10,  7, 2 ]
end subroutine Exam_Johnson_Elongated_Square_Bipyramid_J15

! ---------------------------------------------------------------------------------------

end module Exam_Johnson