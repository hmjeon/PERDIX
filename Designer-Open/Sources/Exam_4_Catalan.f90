!
! ---------------------------------------------------------------------------------------
!
!                                Module for Exam_Catalan 
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/08/03
!                                                            Last  modified  : 2015/07/14
!
! ---------------------------------------------------------------------------------------
!
module Exam_Catalan

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Catalan_Rhombic_Dodecahedron           ! 26. V=14,  E=24,  F=12
    public Exam_Catalan_Rhombic_Triacontahedron        ! 27. V=32,  E=60,  F=30
    public Exam_Catalan_Deltoidal_Icositetrahedron     ! 28. V=26,  E=48,  F=24  (Trapezoidal Icositetrahedron, Tetragonal Icosikaitetrahedron, Strombic Icositetrahedron)
    public Exam_Catalan_Pentagonal_Icositetrahedron    ! 29. V=38,  E=60,  F=24  (Pentagonal Icosikaitetrahedron)
    public Exam_Catalan_Triakis_Octahedron             ! 30. V=14,  E=36,  F=24  (Kisoctahedron)
    public Exam_Catalan_Disdyakis_Dodecahedron         ! 31. V=26,  E=72,  F=48  (Hexakis Octahedron, Kisrhombic Dodecahedron)
    public Exam_Catalan_Triakis_Icosahedron            ! 32. V=32,  E=90,  F=60  (Kisicosahedron)
    public Exam_Catalan_Pentakis_Dodecahedron          ! 33. V=32,  E=90,  F=60  (Kisdodecahedron)
    public Exam_Catalan_Tetrakis_Hexahedron            ! 34. V=14,  E=36,  F=24  (Tetrahexahedron, Kiscube, Jetty snowflake)
    public Exam_Catalan_Triakis_Tetrahedron            ! 35. V=8,   E=18,  F=12  (Kistetrahedron)

    public Exam_Catalan_Disdyakis_Triacontahedron      ! 56. V=62,  E=180, F=120 (Hexakis Icosahedron, Kisrhombic Triacontahedron)    
    public Exam_Catalan_Deltoidal_Hexecontahedron      ! 57. V=62,  E=120, F=60  (Trapezoidal Hexecontahedron, Strombic Hexecontahedron, Tetragonal Hexacontahedron)
    public Exam_Catalan_Pentagonal_Hexecontahedron     ! 58. V=92,  E=150, F=60

contains

! ---------------------------------------------------------------------------------------

! Example of Rhombic Dodecahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Rhombic_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "26_Rhom_Dodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Rhombic dodecahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        !para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
    end if

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.04d0     ! Cylindrical model
    prob.move_x   = 1.2d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 14
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  4.71405d0,  3.33333d0,  8.16497d0 ]; geom.iniP( 2).pos(1:3) = [ -4.71405d0,   6.66667d0,  8.16497d0 ]
    geom.iniP( 3).pos(1:3) = [ -4.71405d0, -3.33333d0,  8.16497d0 ]; geom.iniP( 4).pos(1:3) = [  4.71405d0,  -6.66667d0,  8.16497d0 ]
    geom.iniP( 5).pos(1:3) = [  9.42809d0, -3.33333d0,  0.00000d0 ]; geom.iniP( 6).pos(1:3) = [  9.42809d0,   6.66667d0,  0.00000d0 ]
    geom.iniP( 7).pos(1:3) = [  0.00000d0, 10.00000d0,  0.00000d0 ]; geom.iniP( 8).pos(1:3) = [ -9.42809d0,   3.33333d0,  0.00000d0 ]
    geom.iniP( 9).pos(1:3) = [ -9.42809d0, -6.66667d0,  0.00000d0 ]; geom.iniP(10).pos(1:3) = [  0.00000d0, -10.00000d0,  0.00000d0 ]
    geom.iniP(11).pos(1:3) = [  4.71405d0, -6.66667d0, -8.16497d0 ]; geom.iniP(12).pos(1:3) = [  4.71405d0,   3.33333d0, -8.16497d0 ]
    geom.iniP(13).pos(1:3) = [ -4.71405d0,  6.66667d0, -8.16497d0 ]; geom.iniP(14).pos(1:3) = [ -4.71405d0,  -3.33333d0, -8.16497d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  1,  2,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  4,  5,  6,  1 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  1,  6,  7,  2 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  2,  8,  9,  3 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  3,  9, 10,  4 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  6,  5, 11, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  5,  4, 10, 11 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  2,  7, 13,  8 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  7,  6, 12, 13 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  9,  8, 13, 14 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 10,  9, 14, 11 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 11, 14, 13, 12 ]
end subroutine Exam_Catalan_Rhombic_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Rhombic Triacontahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Rhombic_Triacontahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "27_Rhom_Triaconta"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Rhombic triacontahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.96d0     ! Cylindrical model
    prob.move_x   = 1.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 32
    geom.n_face = 30

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   4.47215d0,   2.76394d0,  13.76385d0 ]; geom.iniP( 2).pos(1:3) = [  -4.47215d0,   7.23609d0,  13.76385d0 ]
    geom.iniP( 3).pos(1:3) = [  -4.47215d0,  -2.76394d0,  13.76385d0 ]; geom.iniP( 4).pos(1:3) = [   4.47215d0,  -7.23609d0,  13.76385d0 ]
    geom.iniP( 5).pos(1:3) = [  11.70824d0,  -2.76394d0,   8.50653d0 ]; geom.iniP( 6).pos(1:3) = [  11.70824d0,   7.23609d0,   8.50653d0 ]
    geom.iniP( 7).pos(1:3) = [   2.76394d0,  11.70824d0,   8.50653d0 ]; geom.iniP( 8).pos(1:3) = [ -11.70824d0,   2.76394d0,   8.50653d0 ]
    geom.iniP( 9).pos(1:3) = [ -11.70824d0,  -7.23609d0,   8.50653d0 ]; geom.iniP(10).pos(1:3) = [  -2.76394d0, -11.70824d0,   8.50653d0 ]
    geom.iniP(11).pos(1:3) = [  14.47218d0,  -7.23609d0,   0.00000d0 ]; geom.iniP(12).pos(1:3) = [  14.47218d0,   2.76394d0,   0.00000d0 ]
    geom.iniP(13).pos(1:3) = [   7.23609d0, -11.70824d0,   5.25732d0 ]; geom.iniP(14).pos(1:3) = [   0.00000d0,  16.18033d0,   0.00000d0 ]
    geom.iniP(15).pos(1:3) = [  -7.23609d0,  11.70824d0,   5.25732d0 ]; geom.iniP(16).pos(1:3) = [   8.94430d0,  11.70824d0,   0.00000d0 ]
    geom.iniP(17).pos(1:3) = [ -14.47218d0,   7.23609d0,   0.00000d0 ]; geom.iniP(18).pos(1:3) = [ -14.47218d0,  -2.76394d0,   0.00000d0 ]
    geom.iniP(19).pos(1:3) = [   0.00000d0, -16.18033d0,   0.00000d0 ]; geom.iniP(20).pos(1:3) = [  -8.94430d0, -11.70824d0,   0.00000d0 ]
    geom.iniP(21).pos(1:3) = [  11.70824d0,  -2.76394d0,  -8.50653d0 ]; geom.iniP(22).pos(1:3) = [  11.70824d0,   7.23609d0,  -8.50653d0 ]
    geom.iniP(23).pos(1:3) = [   7.23609d0, -11.70824d0,  -5.25732d0 ]; geom.iniP(24).pos(1:3) = [  -7.23609d0,  11.70824d0,  -5.25732d0 ]
    geom.iniP(25).pos(1:3) = [   2.76394d0,  11.70824d0,  -8.50653d0 ]; geom.iniP(26).pos(1:3) = [ -11.70824d0,   2.76394d0,  -8.50653d0 ]
    geom.iniP(27).pos(1:3) = [ -11.70824d0,  -7.23609d0,  -8.50653d0 ]; geom.iniP(28).pos(1:3) = [  -2.76394d0, -11.70824d0,  -8.50653d0 ]
    geom.iniP(29).pos(1:3) = [   4.47215d0,  -7.23609d0, -13.76385d0 ]; geom.iniP(30).pos(1:3) = [   4.47215d0,   2.76394d0, -13.76385d0 ]
    geom.iniP(31).pos(1:3) = [  -4.47215d0,   7.23609d0, -13.76385d0 ]; geom.iniP(32).pos(1:3) = [  -4.47215d0,  -2.76394d0, -13.76385d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  1,  2,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  4,  5,  6,  1 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  1,  6,  7,  2 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  2,  8,  9,  3 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  3,  9, 10,  4 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  6,  5, 11, 12 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  5,  4, 13, 11 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  2,  7, 14, 15 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  7,  6, 16, 14 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  9,  8, 17, 18 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  8,  2, 15, 17 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  4, 10, 19, 13 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 10,  9, 20, 19 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 11, 21, 22, 12 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 12, 22, 16,  6 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 13, 19, 23, 11 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 14, 24, 17, 15 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [ 16, 22, 25, 14 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 17, 26, 27, 18 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 18, 27, 20,  9 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 20, 27, 28, 19 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 22, 21, 29, 30 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 21, 11, 23, 29 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 23, 19, 28, 29 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 17, 24, 31, 26 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 24, 14, 25, 31 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 25, 22, 30, 31 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 27, 26, 31, 32 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 28, 27, 32, 29 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 29, 32, 31, 30 ]
end subroutine Exam_Catalan_Rhombic_Triacontahedron

! ---------------------------------------------------------------------------------------

! Example of Deltoidal Icositetrahedron (Trapezoidal Icositetrahedron, Tetragonal Icosikaitetrahedron, Strombic Icositetrahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Deltoidal_Icositetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "28_Deltoi_Icositetra"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Deltoidal icositetrahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.96d0     ! Cylindrical model
    prob.move_x   = 2.5d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 26
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   5.57791d0,   6.46448d0,  14.57578d0 ]; geom.iniP( 2).pos(1:3) = [  -4.31428d0,   5.00000d0,  14.57578d0 ]
    geom.iniP( 3).pos(1:3) = [  -7.21164d0,  -4.57107d0,  14.57578d0 ]; geom.iniP( 4).pos(1:3) = [   5.57791d0,  -6.46448d0,  14.57578d0 ]
    geom.iniP( 5).pos(1:3) = [  15.09999d0,  -4.57107d0,   6.03748d0 ]; geom.iniP( 6).pos(1:3) = [  12.94286d0,   5.00000d0,   7.97176d0 ]
    geom.iniP( 7).pos(1:3) = [   2.31045d0,  15.60661d0,   6.03748d0 ]; geom.iniP( 8).pos(1:3) = [  -9.52209d0,  11.03555d0,   8.53830d0 ]
    geom.iniP( 9).pos(1:3) = [ -15.77673d0,   0.00000d0,   6.03748d0 ]; geom.iniP(10).pos(1:3) = [  -6.84153d0, -12.07108d0,   7.97176d0 ]
    geom.iniP(11).pos(1:3) = [   2.31045d0, -15.60661d0,   6.03748d0 ]; geom.iniP(12).pos(1:3) = [  12.78954d0,  11.03555d0,   0.00000d0 ]
    geom.iniP(13).pos(1:3) = [  15.77673d0,   0.00000d0,  -6.03748d0 ]; geom.iniP(14).pos(1:3) = [  10.41560d0, -12.07108d0,   1.36774d0 ]
    geom.iniP(15).pos(1:3) = [  -2.31045d0,  15.60661d0,  -6.03748d0 ]; geom.iniP(16).pos(1:3) = [ -10.41560d0,  12.07108d0,  -1.36774d0 ]
    geom.iniP(17).pos(1:3) = [ -12.78954d0, -11.03555d0,   0.00000d0 ]; geom.iniP(18).pos(1:3) = [ -15.09999d0,   4.57107d0,  -6.03748d0 ]
    geom.iniP(19).pos(1:3) = [  -2.31045d0, -15.60661d0,  -6.03748d0 ]; geom.iniP(20).pos(1:3) = [   6.84153d0,  12.07108d0,  -7.97176d0 ]
    geom.iniP(21).pos(1:3) = [   9.52209d0, -11.03555d0,  -8.53830d0 ]; geom.iniP(22).pos(1:3) = [   7.21164d0,   4.57107d0, -14.57578d0 ]
    geom.iniP(23).pos(1:3) = [  -5.57791d0,   6.46448d0, -14.57578d0 ]; geom.iniP(24).pos(1:3) = [ -12.94286d0,  -5.00000d0,  -7.97176d0 ]
    geom.iniP(25).pos(1:3) = [  -5.57791d0,  -6.46448d0, -14.57578d0 ]; geom.iniP(26).pos(1:3) = [   4.31428d0,  -5.00000d0, -14.57578d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  1,  2,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  4,  5,  6,  1 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  1,  7,  8,  2 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  2,  8,  9,  3 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  3, 10, 11,  4 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  1,  6, 12,  7 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  6,  5, 13, 12 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  5,  4, 11, 14 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  8,  7, 15, 16 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  3,  9, 17, 10 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  9,  8, 16, 18 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 11, 10, 17, 19 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 12, 20, 15,  7 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  5, 14, 21, 13 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [ 13, 22, 20, 12 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 11, 19, 21, 14 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 15, 23, 18, 16 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  9, 18, 24, 17 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 17, 24, 25, 19 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 15, 20, 22, 23 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 13, 21, 26, 22 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 21, 19, 25, 26 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 18, 23, 25, 24 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 22, 26, 25, 23 ]
end subroutine Exam_Catalan_Deltoidal_Icositetrahedron

! ---------------------------------------------------------------------------------------

! Example of Pentagonal Icositetrahedron (Pentagonal Icosikaitetrahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Pentagonal_Icositetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "29_Penta_Icositetra"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal icositetrahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 3.0d0      ! Cylindrical model
    prob.move_y   = 1.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 38
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   7.82010d0,   5.00001d0,  19.50683d0 ]; geom.iniP( 2).pos(1:3) = [  -1.25679d0,   9.19644d0,  19.50683d0 ]
    geom.iniP( 3).pos(1:3) = [  -8.87492d0,   2.71846d0,  19.50683d0 ]; geom.iniP( 4).pos(1:3) = [  -6.19181d0,  -6.91489d0,  19.50683d0 ]
    geom.iniP( 5).pos(1:3) = [   7.82010d0,  -9.19644d0,  19.50683d0 ]; geom.iniP( 6).pos(1:3) = [  17.95182d0,  -6.91489d0,   9.82789d0 ]
    geom.iniP( 7).pos(1:3) = [  19.89193d0,   2.71846d0,   7.97449d0 ]; geom.iniP( 8).pos(1:3) = [  14.38341d0,   9.19644d0,  13.23682d0 ]
    geom.iniP( 9).pos(1:3) = [  10.81503d0,  16.91491d0,   7.97449d0 ]; geom.iniP(10).pos(1:3) = [  -2.75425d0,  20.15389d0,  10.60564d0 ]
    geom.iniP(11).pos(1:3) = [ -14.95692d0,  13.39288d0,   7.97449d0 ]; geom.iniP(12).pos(1:3) = [ -16.32352d0,   5.00001d0,  13.23682d0 ]
    geom.iniP(13).pos(1:3) = [ -21.38942d0,  -5.95744d0,   5.76617d0 ]; geom.iniP(14).pos(1:3) = [ -11.38851d0, -12.71847d0,  13.23682d0 ]
    geom.iniP(15).pos(1:3) = [  -5.88002d0, -19.19647d0,   7.97449d0 ]; geom.iniP(16).pos(1:3) = [   3.93990d0, -18.82978d0,   9.82789d0 ]
    geom.iniP(17).pos(1:3) = [   9.76021d0, -19.19647d0,   1.70447d0 ]; geom.iniP(18).pos(1:3) = [  17.37833d0, -12.71847d0,   1.70447d0 ]
    geom.iniP(19).pos(1:3) = [  18.83711d0,  -7.95600d0,  -6.96681d0 ]; geom.iniP(20).pos(1:3) = [  21.38942d0,   5.95744d0,  -5.76617d0 ]
    geom.iniP(21).pos(1:3) = [  13.32862d0,  16.91491d0,  -1.70447d0 ]; geom.iniP(22).pos(1:3) = [   5.88001d0,  19.19646d0,  -7.97449d0 ]
    geom.iniP(23).pos(1:3) = [  -2.88509d0,  21.11138d0,  -3.55788d0 ]; geom.iniP(24).pos(1:3) = [ -10.81503d0,  16.91491d0,  -7.97449d0 ]
    geom.iniP(25).pos(1:3) = [ -17.37833d0,  12.71846d0,  -1.70447d0 ]; geom.iniP(26).pos(1:3) = [ -20.77720d0,   3.75954d0,  -4.56555d0 ]
    geom.iniP(27).pos(1:3) = [ -15.84219d0, -13.95891d0,  -4.56555d0 ]; geom.iniP(28).pos(1:3) = [  -8.30143d0, -19.87089d0,  -1.70447d0 ]
    geom.iniP(29).pos(1:3) = [   2.75425d0, -20.15389d0, -10.60564d0 ]; geom.iniP(30).pos(1:3) = [  12.44332d0, -10.43691d0, -14.24450d0 ]
    geom.iniP(31).pos(1:3) = [   8.87492d0,  -2.71846d0, -19.50683d0 ]; geom.iniP(32).pos(1:3) = [  12.27381d0,   6.24046d0, -16.64576d0 ]
    geom.iniP(33).pos(1:3) = [   5.30651d0,  13.39288d0, -16.09789d0 ]; geom.iniP(34).pos(1:3) = [  -7.82010d0,   9.19644d0, -19.50683d0 ]
    geom.iniP(35).pos(1:3) = [ -17.06654d0,   0.43689d0, -13.23682d0 ]; geom.iniP(36).pos(1:3) = [ -14.38341d0,  -9.19644d0, -13.23682d0 ]
    geom.iniP(37).pos(1:3) = [  -5.61832d0, -11.11133d0, -17.65343d0 ]; geom.iniP(38).pos(1:3) = [  -0.94500d0,  -3.08512d0, -21.36021d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [  5,  1,  2,  3,  4 ]
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [  1,  5,  6,  7,  8 ]
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [  2,  1,  8,  9, 10 ]
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [  3,  2, 10, 11, 12 ]
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [  4,  3, 12, 13, 14 ]
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [  5,  4, 14, 15, 16 ]
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [  6,  5, 16, 17, 18 ]
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [  7,  6, 18, 19, 20 ]
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [  8,  7, 20, 21,  9 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 10,  9, 21, 22, 23 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 11, 10, 23, 24, 25 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 12, 11, 25, 26, 13 ]
    geom.face(13).n_poi = 5; allocate(geom.face(13).poi(5)); geom.face(13).poi(1:5) = [ 14, 13, 27, 28, 15 ]
    geom.face(14).n_poi = 5; allocate(geom.face(14).poi(5)); geom.face(14).poi(1:5) = [ 16, 15, 28, 29, 17 ]
    geom.face(15).n_poi = 5; allocate(geom.face(15).poi(5)); geom.face(15).poi(1:5) = [ 18, 17, 29, 30, 19 ]
    geom.face(16).n_poi = 5; allocate(geom.face(16).poi(5)); geom.face(16).poi(1:5) = [ 20, 19, 30, 31, 32 ]
    geom.face(17).n_poi = 5; allocate(geom.face(17).poi(5)); geom.face(17).poi(1:5) = [ 21, 20, 32, 33, 22 ]
    geom.face(18).n_poi = 5; allocate(geom.face(18).poi(5)); geom.face(18).poi(1:5) = [ 23, 22, 33, 34, 24 ]
    geom.face(19).n_poi = 5; allocate(geom.face(19).poi(5)); geom.face(19).poi(1:5) = [ 25, 24, 34, 35, 26 ]
    geom.face(20).n_poi = 5; allocate(geom.face(20).poi(5)); geom.face(20).poi(1:5) = [ 13, 26, 35, 36, 27 ]
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [ 28, 27, 36, 37, 29 ]
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 30, 29, 37, 38, 31 ]
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 32, 31, 38, 34, 33 ]
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 35, 34, 38, 37, 36 ]
end subroutine Exam_Catalan_Pentagonal_Icositetrahedron

! ---------------------------------------------------------------------------------------

! Example of Triakis Octahedron (Kisoctahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Triakis_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "30_Tria_Octa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triakis octahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 2.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 14
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  2.39921d0,   8.53554d0,   8.19141d0 ]; geom.iniP( 2).pos(1:3) = [ -6.49491d0,  -6.03554d0,   8.19141d0 ]
    geom.iniP( 3).pos(1:3) = [  2.39921d0,  -1.46447d0,   8.19141d0 ]; geom.iniP( 4).pos(1:3) = [  9.88791d0,  -6.03554d0,   3.39300d0 ]
    geom.iniP( 5).pos(1:3) = [ -5.79220d0,   3.53554d0,   5.38056d0 ]; geom.iniP( 6).pos(1:3) = [  7.77977d0,   3.53554d0,   1.40543d0 ]
    geom.iniP( 7).pos(1:3) = [ -9.88791d0,   6.03554d0,  -3.39300d0 ]; geom.iniP( 8).pos(1:3) = [  0.41164d0,  -8.53554d0,   1.40543d0 ]
    geom.iniP( 9).pos(1:3) = [  6.49491d0,   6.03554d0,  -8.19141d0 ]; geom.iniP(10).pos(1:3) = [ -7.77977d0,  -3.53554d0,  -1.40543d0 ]
    geom.iniP(11).pos(1:3) = [ -0.41164d0,   8.53554d0,  -1.40543d0 ]; geom.iniP(12).pos(1:3) = [ -2.39921d0,  -8.53554d0,  -8.19141d0 ]
    geom.iniP(13).pos(1:3) = [  5.79220d0,  -3.53554d0,  -5.38056d0 ]; geom.iniP(14).pos(1:3) = [ -2.39921d0,   1.46447d0,  -8.19141d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  5,  2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  2,  4,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  6 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  2,  5,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  5,  1,  7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  4,  2,  8 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4,  9,  6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  6,  9,  1 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7, 10,  2 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  1, 11,  7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  2, 12,  8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  8, 12,  4 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  9,  4, 13 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  1,  9, 11 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  2, 10, 12 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 10,  7, 12 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  7, 11,  9 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  4, 12, 13 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 13, 12,  9 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [  7, 14, 12 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  9, 14,  7 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [  9, 12, 14 ]
end subroutine Exam_Catalan_Triakis_Octahedron

! ---------------------------------------------------------------------------------------

! Example of Disdyakis Dodecahedron (Hexakis Octahedron, Kisrhombic Dodecahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Disdyakis_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "31_Disdya_Dodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Disdyakis dodecahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.96d0     ! Cylindrical model
    prob.move_x   =-1.5d0      ! Cylindrical model
    prob.move_y   =-1.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 26
    geom.n_face = 48

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   3.36696d0,   3.53554d0,  15.23910d0 ]; geom.iniP( 2).pos(1:3) = [  -6.62112d0,   3.04738d0,  15.23910d0 ]
    geom.iniP( 3).pos(1:3) = [   3.36696d0,  -9.84156d0,  15.23910d0 ]; geom.iniP( 4).pos(1:3) = [  12.42529d0,   3.04738d0,  11.03093d0 ]
    geom.iniP( 5).pos(1:3) = [   2.12322d0,  15.60662d0,   9.60980d0 ]; geom.iniP( 6).pos(1:3) = [  -8.98385d0,  -6.03555d0,  11.78682d0 ]
    geom.iniP( 7).pos(1:3) = [  13.11357d0,  -6.03555d0,   6.90455d0 ]; geom.iniP( 8).pos(1:3) = [  -9.74661d0,   9.57108d0,   8.33454d0 ]
    geom.iniP( 9).pos(1:3) = [  -8.86576d0, -13.45179d0,   5.07970d0 ]; geom.iniP(10).pos(1:3) = [ -18.01612d0,   0.00000d0,   3.98051d0 ]
    geom.iniP(11).pos(1:3) = [  12.35083d0,   9.57108d0,   3.45228d0 ]; geom.iniP(12).pos(1:3) = [  10.18065d0, -13.45179d0,   0.87154d0 ]
    geom.iniP(13).pos(1:3) = [  18.01612d0,   0.00000d0,  -3.98051d0 ]; geom.iniP(14).pos(1:3) = [ -10.18065d0,  13.45179d0,  -0.87154d0 ]
    geom.iniP(15).pos(1:3) = [   0.76275d0, -15.60662d0,   3.45228d0 ]; geom.iniP(16).pos(1:3) = [   8.86576d0,  13.45179d0,  -5.07970d0 ]
    geom.iniP(17).pos(1:3) = [  -0.76275d0,  15.60662d0,  -3.45228d0 ]; geom.iniP(18).pos(1:3) = [ -12.35083d0,  -9.57108d0,  -3.45228d0 ]
    geom.iniP(19).pos(1:3) = [  -2.12322d0, -15.60662d0,  -9.60980d0 ]; geom.iniP(20).pos(1:3) = [   9.74661d0,  -9.57108d0,  -8.33454d0 ]
    geom.iniP(21).pos(1:3) = [ -13.11357d0,   6.03555d0,  -6.90455d0 ]; geom.iniP(22).pos(1:3) = [  -3.36696d0,   9.84156d0, -15.23910d0 ]
    geom.iniP(23).pos(1:3) = [ -12.42529d0,  -3.04738d0, -11.03093d0 ]; geom.iniP(24).pos(1:3) = [   8.98385d0,   6.03555d0, -11.78682d0 ]
    geom.iniP(25).pos(1:3) = [   6.62112d0,  -3.04738d0, -15.23910d0 ]; geom.iniP(26).pos(1:3) = [  -3.36696d0,  -3.53554d0, -15.23910d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  5,  2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  2,  6,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  5 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  4,  3,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  2,  5,  8 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  3,  6,  9 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  6,  2, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4, 11,  5 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  3, 12,  7 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  7, 13,  4 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  5, 14,  8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  8, 10,  2 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  6, 10,  9 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  9, 15,  3 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  5, 11, 16 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 11,  4, 13 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  7, 12, 13 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 12,  3, 15 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [  8, 14, 10 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 14,  5, 17 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  9, 10, 18 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 15,  9, 19 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 11, 13, 16 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 16, 17,  5 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 12, 20, 13 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 15, 19, 12 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 14, 21, 10 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 17, 22, 14 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 10, 23, 18 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 18, 19,  9 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 16, 13, 24 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 17, 16, 22 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 13, 20, 25 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 20, 12, 19 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 10, 21, 23 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 21, 14, 22 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 18, 23, 19 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 13, 25, 24 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 24, 22, 16 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 20, 19, 25 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 21, 22, 23 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 23, 26, 19 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 24, 25, 22 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 25, 19, 26 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 23, 22, 26 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 25, 26, 22 ]
end subroutine Exam_Catalan_Disdyakis_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Triakis Icosahedron (Kisicosahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Triakis_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "32_Tria_Icosa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triakis icosahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.97d0     ! Cylindrical model
    prob.move_x   = 2.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! Allocate point and face structure
    geom.n_iniP = 32
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   2.34795d0,   8.61805d0,  13.74519d0 ]; geom.iniP( 2).pos(1:3) = [  -6.39492d0,  -6.23608d0,  13.74519d0 ]
    geom.iniP( 3).pos(1:3) = [   2.34795d0,  -1.38197d0,  13.74519d0 ]; geom.iniP( 4).pos(1:3) = [  10.59508d0,  -6.23608d0,  10.84295d0 ]
    geom.iniP( 5).pos(1:3) = [  -6.14703d0,   3.61804d0,  12.06138d0 ]; geom.iniP( 6).pos(1:3) = [   9.80231d0,   3.61804d0,   9.33690d0 ]
    geom.iniP( 7).pos(1:3) = [ -13.09623d0,   7.70821d0,   6.14703d0 ]; geom.iniP( 8).pos(1:3) = [   1.59493d0, -10.32625d0,   9.33690d0 ]
    geom.iniP( 9).pos(1:3) = [  14.39413d0,   7.70821d0,   1.45111d0 ]; geom.iniP(10).pos(1:3) = [ -12.15026d0,  -2.23608d0,   6.61244d0 ]
    geom.iniP(11).pos(1:3) = [  -3.94289d0,  11.70823d0,   6.61244d0 ]; geom.iniP(12).pos(1:3) = [   0.24788d0, -16.32625d0,   1.45111d0 ]
    geom.iniP(13).pos(1:3) = [  13.65631d0,  -2.23608d0,   2.20415d0 ]; geom.iniP(14).pos(1:3) = [   5.91435d0,  11.70823d0,   4.92862d0 ]
    geom.iniP(15).pos(1:3) = [ -14.39413d0,  -7.70821d0,  -1.45111d0 ]; geom.iniP(16).pos(1:3) = [  -0.24788d0,  16.32625d0,  -1.45111d0 ]
    geom.iniP(17).pos(1:3) = [  -7.36546d0, -10.85411d0,   4.92862d0 ]; geom.iniP(18).pos(1:3) = [   8.58388d0, -10.85411d0,   2.20415d0 ]
    geom.iniP(19).pos(1:3) = [  13.09623d0,  -7.70821d0,  -6.14703d0 ]; geom.iniP(20).pos(1:3) = [ -13.65631d0,   2.23608d0,  -2.20415d0 ]
    geom.iniP(21).pos(1:3) = [  -8.58388d0,  10.85411d0,  -2.20415d0 ]; geom.iniP(22).pos(1:3) = [  12.15026d0,   2.23608d0,  -6.61244d0 ]
    geom.iniP(23).pos(1:3) = [   7.36546d0,  10.85411d0,  -4.92862d0 ]; geom.iniP(24).pos(1:3) = [ -10.59508d0,   6.23608d0, -10.84295d0 ]
    geom.iniP(25).pos(1:3) = [  -5.91435d0, -11.70823d0,  -4.92862d0 ]; geom.iniP(26).pos(1:3) = [   3.94289d0, -11.70823d0,  -6.61244d0 ]
    geom.iniP(27).pos(1:3) = [   6.39492d0,   6.23608d0, -13.74519d0 ]; geom.iniP(28).pos(1:3) = [  -9.80231d0,  -3.61804d0,  -9.33690d0 ]
    geom.iniP(29).pos(1:3) = [  -1.59493d0,  10.32625d0,  -9.33690d0 ]; geom.iniP(30).pos(1:3) = [  -2.34795d0,  -8.61805d0, -13.74519d0 ]
    geom.iniP(31).pos(1:3) = [   6.14703d0,  -3.61804d0, -12.06138d0 ]; geom.iniP(32).pos(1:3) = [  -2.34795d0,   1.38197d0, -13.74519d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  5,  2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  2,  4,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  6 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  2,  5,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  5,  1,  7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  4,  2,  8 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4,  9,  6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  6,  9,  1 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7, 10,  2 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  1, 11,  7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  2, 12,  8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  8, 12,  4 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  9,  4, 13 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  1,  9, 14 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  2, 10, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 10,  7, 15 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  7, 11, 16 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 11,  1, 16 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 12,  2, 17 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [  4, 12, 18 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  4, 19, 13 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 13, 19,  9 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [  9, 16, 14 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 14, 16,  1 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 15, 17,  2 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [  7, 20, 15 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 16, 21,  7 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 17, 15, 12 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 12, 19, 18 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 18, 19,  4 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [  9, 19, 22 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 16,  9, 23 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 15, 20, 24 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 20,  7, 24 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [  7, 21, 24 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 21, 16, 24 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 12, 15, 25 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 19, 12, 26 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 19, 27, 22 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 22, 27,  9 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [  9, 27, 23 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 23, 27, 16 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 24, 28, 15 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 16, 29, 24 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 15, 30, 25 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 25, 30, 12 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 12, 30, 26 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 26, 30, 19 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 27, 19, 31 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 16, 27, 29 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 15, 28, 30 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 28, 24, 30 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 24, 29, 27 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 19, 30, 31 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 31, 30, 27 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 24, 32, 30 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 27, 32, 24 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 27, 30, 32 ]
end subroutine Exam_Catalan_Triakis_Icosahedron

! ---------------------------------------------------------------------------------------

! Example of Pentakis Dodecahedron (Kisdodecahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Pentakis_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "33_Penta_Dodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentakis dodecahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 0.8d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "YZ"

    ! Allocate point and face structure
    geom.n_iniP = 32
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   2.97755d0,   5.63662d0,  14.45333d0 ]; geom.iniP( 2).pos(1:3) = [  -6.33420d0,  -0.71767d0,  14.45333d0 ]
    geom.iniP( 3).pos(1:3) = [   2.97755d0,  -4.36339d0,  14.45333d0 ]; geom.iniP( 4).pos(1:3) = [  11.53109d0,  -0.71767d0,  10.77288d0 ]
    geom.iniP( 5).pos(1:3) = [  -5.95510d0,   8.72679d0,  11.18856d0 ]; geom.iniP( 6).pos(1:3) = [  -3.53565d0, -10.99911d0,  10.77288d0 ]
    geom.iniP( 7).pos(1:3) = [   9.89246d0,   8.72679d0,   7.92378d0 ]; geom.iniP( 8).pos(1:3) = [   7.50572d0, -10.99911d0,   8.49823d0 ]
    geom.iniP( 9).pos(1:3) = [ -13.60562d0,   3.75774d0,   7.09242d0 ]; geom.iniP(10).pos(1:3) = [   1.46112d0,  14.03918d0,   7.09242d0 ]
    geom.iniP(11).pos(1:3) = [ -11.73268d0,  -6.03006d0,   7.92378d0 ]; geom.iniP(12).pos(1:3) = [  15.30108d0,   3.75774d0,   1.13732d0 ]
    geom.iniP(13).pos(1:3) = [  13.90920d0,  -6.03006d0,   2.64125d0 ]; geom.iniP(14).pos(1:3) = [  -8.78785d0,  12.87798d0,   2.54314d0 ]
    geom.iniP(15).pos(1:3) = [   0.54413d0, -15.15027d0,   2.64125d0 ]; geom.iniP(16).pos(1:3) = [  -9.07746d0, -12.87798d0,   1.13732d0 ]
    geom.iniP(17).pos(1:3) = [   9.07746d0,  12.87798d0,  -1.13732d0 ]; geom.iniP(18).pos(1:3) = [   8.78785d0, -12.87798d0,  -2.54314d0 ]
    geom.iniP(19).pos(1:3) = [ -13.90920d0,   6.03006d0,  -2.64125d0 ]; geom.iniP(20).pos(1:3) = [ -15.30108d0,  -3.75774d0,  -1.13732d0 ]
    geom.iniP(21).pos(1:3) = [  -0.54413d0,  15.15027d0,  -2.64125d0 ]; geom.iniP(22).pos(1:3) = [  11.73268d0,   6.03006d0,  -7.92378d0 ]
    geom.iniP(23).pos(1:3) = [  13.60562d0,  -3.75774d0,  -7.09242d0 ]; geom.iniP(24).pos(1:3) = [  -7.50572d0,  10.99911d0,  -8.49823d0 ]
    geom.iniP(25).pos(1:3) = [  -1.46112d0, -14.03918d0,  -7.09242d0 ]; geom.iniP(26).pos(1:3) = [  -9.89246d0,  -8.72679d0,  -7.92378d0 ]
    geom.iniP(27).pos(1:3) = [   3.53565d0,  10.99911d0, -10.77288d0 ]; geom.iniP(28).pos(1:3) = [   5.95510d0,  -8.72679d0, -11.18856d0 ]
    geom.iniP(29).pos(1:3) = [ -11.53109d0,   0.71767d0, -10.77288d0 ]; geom.iniP(30).pos(1:3) = [   6.33420d0,   0.71767d0, -14.45333d0 ]
    geom.iniP(31).pos(1:3) = [  -2.97755d0,   4.36339d0, -14.45333d0 ]; geom.iniP(32).pos(1:3) = [  -2.97755d0,  -5.63662d0, -14.45333d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  5,  2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  2,  6,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  4,  3,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  2,  5,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  5,  1, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  3,  6,  8 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  6,  2, 11 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4, 12,  7 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  7, 10,  1 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  8, 13,  4 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  5, 14,  9 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  9, 11,  2 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 10, 14,  5 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  6, 15,  8 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 11, 16,  6 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  7, 12, 17 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 12,  4, 13 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 10,  7, 17 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 13,  8, 18 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  9, 14, 19 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 11,  9, 20 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 14, 10, 21 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [  8, 15, 18 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 15,  6, 16 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 16, 11, 20 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 12, 22, 17 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 13, 23, 12 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 17, 21, 10 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 18, 23, 13 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 14, 24, 19 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 19, 20,  9 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 21, 24, 14 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 15, 25, 18 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 16, 25, 15 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 20, 26, 16 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 17, 22, 27 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 22, 12, 23 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 21, 17, 27 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 23, 18, 28 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 19, 24, 29 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 20, 19, 29 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 24, 21, 27 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 18, 25, 28 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 25, 16, 26 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 26, 20, 29 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 22, 30, 27 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 23, 30, 22 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 28, 30, 23 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 24, 31, 29 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 27, 31, 24 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 25, 32, 28 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 26, 32, 25 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 29, 32, 26 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 27, 30, 31 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 30, 28, 32 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 29, 31, 32 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 30, 32, 31 ]
end subroutine Exam_Catalan_Pentakis_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Tetrakis Hexahedron (Tetrahexahedron, Kiscube, Jetty Snowflake)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Tetrakis_Hexahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "34_Tetra_Hexa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tetrakis hexahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.95d0     ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   =-0.5d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 14
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  2.98143d0,   6.66668d0,   8.94429d0 ]; geom.iniP( 2).pos(1:3) = [  -6.95667d0,  -2.22222d0,   8.94429d0 ]
    geom.iniP( 3).pos(1:3) = [  2.98143d0,  -3.33334d0,   8.94429d0 ]; geom.iniP( 4).pos(1:3) = [  10.93189d0,  -2.22222d0,   2.98143d0 ]
    geom.iniP( 5).pos(1:3) = [ -5.96286d0,   6.66668d0,   4.47215d0 ]; geom.iniP( 6).pos(1:3) = [   0.99381d0, -11.11112d0,   2.98143d0 ]
    geom.iniP( 7).pos(1:3) = [  7.45357d0,   6.66668d0,   0.00000d0 ]; geom.iniP( 8).pos(1:3) = [ -10.93189d0,   2.22222d0,  -2.98143d0 ]
    geom.iniP( 9).pos(1:3) = [ -0.99381d0,  11.11112d0,  -2.98143d0 ]; geom.iniP(10).pos(1:3) = [  -7.45357d0,  -6.66668d0,   0.00000d0 ]
    geom.iniP(11).pos(1:3) = [  6.95667d0,   2.22222d0,  -8.94429d0 ]; geom.iniP(12).pos(1:3) = [   5.96286d0,  -6.66668d0,  -4.47215d0 ]
    geom.iniP(13).pos(1:3) = [ -2.98143d0,  -6.66668d0,  -8.94429d0 ]; geom.iniP(14).pos(1:3) = [  -2.98143d0,   3.33334d0,  -8.94429d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  1,  2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  5,  2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  2,  6,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  4,  3,  6 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  2,  5,  8 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  5,  1,  9 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  6,  2, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4, 11,  7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7,  9,  1 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  6, 12,  4 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  5,  9,  8 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  8, 10,  2 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 10, 13,  6 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  7, 11,  9 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 11,  4, 12 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 12,  6, 13 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  8,  9, 14 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 10,  8, 13 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 11, 14,  9 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 12, 13, 11 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 14, 13,  8 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 14, 11, 13 ]
end subroutine Exam_Catalan_Tetrakis_Hexahedron

! ---------------------------------------------------------------------------------------

! Example of Triakis Tetrahedron (Kistetrahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Triakis_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "35_Tria_Tetra"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Triakis tetrahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 0.95d0     ! Atomic model
    prob.size     = 1.04d0     ! Cylindrical model
    prob.move_x   = 2.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 8
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  2.51260d0,  8.33333d0,  5.33003d0 ]
    geom.iniP(2).pos(1:3) = [ -6.70026d0, -5.55556d0,  5.33003d0 ]
    geom.iniP(3).pos(1:3) = [  2.51260d0, -1.66667d0,  5.33003d0 ]
    geom.iniP(4).pos(1:3) = [  8.37534d0, -5.55556d0, -1.77667d0 ]
    geom.iniP(5).pos(1:3) = [ -5.02520d0,  3.33334d0,  1.06601d0 ]
    geom.iniP(6).pos(1:3) = [  4.02016d0,  3.33334d0, -3.19801d0 ]
    geom.iniP(7).pos(1:3) = [ -4.18767d0,  2.77778d0, -8.88340d0 ]
    geom.iniP(8).pos(1:3) = [ -1.50756d0, -5.00001d0, -3.19801d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 3, 1, 2 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 3, 4, 1 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 1, 5, 2 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 2, 4, 3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 1, 4, 6 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 2, 5, 7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 5, 1, 7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 4, 2, 8 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 4, 7, 6 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 6, 7, 1 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 7, 8, 2 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 8, 7, 4 ]
end subroutine Exam_Catalan_Triakis_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of Disdyakis Triacontahedron (Hexakis Icosahedron, Kisrhombic Triacontahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Disdyakis_Triacontahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "56_Disdy_Triacontahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Disdyakis Triacontahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   =-1.9d0      ! Cylindrical model
    prob.move_y   =-1.5d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! Allocate point and face structure
    geom.n_iniP = 62
    geom.n_face = 120

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   3.55493d0,   3.61805d0,  26.79975d0 ]; geom.iniP( 2).pos(1:3) = [  -6.44354d0,   3.44209d0,  26.79975d0 ]
    geom.iniP( 3).pos(1:3) = [   3.55493d0, -12.09020d0,  26.79975d0 ]; geom.iniP( 4).pos(1:3) = [  13.20762d0,   3.44209d0,  24.19306d0 ]
    geom.iniP( 5).pos(1:3) = [   3.01178d0,  18.77348d0,  22.70511d0 ]; geom.iniP( 6).pos(1:3) = [ -10.49633d0,  -5.42707d0,  24.58309d0 ]
    geom.iniP( 7).pos(1:3) = [  16.54235d0,  -5.42707d0,  20.99647d0 ]; geom.iniP( 8).pos(1:3) = [ -10.79036d0,  11.28117d0,  22.36645d0 ]
    geom.iniP( 9).pos(1:3) = [ -13.52817d0, -13.76834d0,  19.97535d0 ]; geom.iniP(10).pos(1:3) = [ -22.94388d0,   2.06525d0,  18.61046d0 ]
    geom.iniP(11).pos(1:3) = [  16.24832d0,  11.28117d0,  18.77981d0 ]; geom.iniP(12).pos(1:3) = [  18.26807d0, -13.76834d0,  15.75765d0 ]
    geom.iniP(13).pos(1:3) = [  27.00235d0,   2.06525d0,  11.98520d0 ]; geom.iniP(14).pos(1:3) = [ -14.08765d0,  18.02300d0,  15.75765d0 ]
    geom.iniP(15).pos(1:3) = [  -6.19306d0, -20.06235d0,  17.40983d0 ]; geom.iniP(16).pos(1:3) = [  17.70862d0,  18.02300d0,  11.53994d0 ]
    geom.iniP(17).pos(1:3) = [  10.51777d0, -20.06235d0,  15.19320d0 ]; geom.iniP(18).pos(1:3) = [  -6.96284d0,  23.68038d0,  11.60656d0 ]
    geom.iniP(19).pos(1:3) = [ -20.53835d0, -12.39921d0,  12.97653d0 ]; geom.iniP(20).pos(1:3) = [   1.74445d0, -24.40498d0,  13.15096d0 ]
    geom.iniP(21).pos(1:3) = [ -14.99482d0, -24.96925d0,   5.35994d0 ]; geom.iniP(22).pos(1:3) = [   9.74799d0,  23.68038d0,   9.38992d0 ]
    geom.iniP(23).pos(1:3) = [  23.21117d0, -12.39921d0,   7.17325d0 ]; geom.iniP(24).pos(1:3) = [  15.87365d0, -24.96925d0,   1.26531d0 ]
    geom.iniP(25).pos(1:3) = [ -21.01410d0,  14.63528d0,   9.38992d0 ]; geom.iniP(26).pos(1:3) = [   0.83920d0,  27.03449d0,   6.32656d0 ]
    geom.iniP(27).pos(1:3) = [ -15.87365d0,  24.96925d0,  -1.26531d0 ]; geom.iniP(28).pos(1:3) = [ -25.55084d0,  -9.82406d0,   4.71555d0 ]
    geom.iniP(29).pos(1:3) = [  22.73541d0,  14.63528d0,   3.58661d0 ]; geom.iniP(30).pos(1:3) = [  14.99482d0,  24.96925d0,  -5.35994d0 ]
    geom.iniP(31).pos(1:3) = [  25.89661d0,  -9.82406d0,  -2.10885d0 ]; geom.iniP(32).pos(1:3) = [ -25.89661d0,   9.82406d0,   2.10885d0 ]
    geom.iniP(33).pos(1:3) = [ -27.03855d0,   0.00000d0,   3.58661d0 ]; geom.iniP(34).pos(1:3) = [   0.47576d0, -27.03449d0,   3.58661d0 ]
    geom.iniP(35).pos(1:3) = [  25.55084d0,   9.82406d0,  -4.71555d0 ]; geom.iniP(36).pos(1:3) = [  27.03855d0,   0.00000d0,  -3.58661d0 ]
    geom.iniP(37).pos(1:3) = [  -0.47576d0,  27.03449d0,  -3.58661d0 ]; geom.iniP(38).pos(1:3) = [ -22.73541d0, -14.63528d0,  -3.58661d0 ]
    geom.iniP(39).pos(1:3) = [ -27.00235d0,  -2.06525d0, -11.98520d0 ]; geom.iniP(40).pos(1:3) = [  -0.83920d0, -27.03449d0,  -6.32656d0 ]
    geom.iniP(41).pos(1:3) = [  21.01410d0, -14.63528d0,  -9.38992d0 ]; geom.iniP(42).pos(1:3) = [  22.94388d0,  -2.06525d0, -18.61046d0 ]
    geom.iniP(43).pos(1:3) = [ -23.21117d0,  12.39921d0,  -7.17325d0 ]; geom.iniP(44).pos(1:3) = [  -1.74445d0,  24.40498d0, -13.15096d0 ]
    geom.iniP(45).pos(1:3) = [ -17.70862d0, -18.02300d0, -11.53994d0 ]; geom.iniP(46).pos(1:3) = [  -9.74799d0, -23.68038d0,  -9.38992d0 ]
    geom.iniP(47).pos(1:3) = [  20.53835d0,  12.39921d0, -12.97653d0 ]; geom.iniP(48).pos(1:3) = [  14.08765d0, -18.02300d0, -15.75765d0 ]
    geom.iniP(49).pos(1:3) = [   6.96284d0, -23.68038d0, -11.60656d0 ]; geom.iniP(50).pos(1:3) = [ -18.26807d0,  13.76834d0, -15.75765d0 ]
    geom.iniP(51).pos(1:3) = [ -10.51777d0,  20.06235d0, -15.19320d0 ]; geom.iniP(52).pos(1:3) = [  -3.01178d0, -18.77348d0, -22.70511d0 ]
    geom.iniP(53).pos(1:3) = [  13.52817d0,  13.76834d0, -19.97535d0 ]; geom.iniP(54).pos(1:3) = [   6.19306d0,  20.06235d0, -17.40983d0 ]
    geom.iniP(55).pos(1:3) = [  -3.55493d0,  12.09020d0, -26.79975d0 ]; geom.iniP(56).pos(1:3) = [ -16.24832d0, -11.28117d0, -18.77981d0 ]
    geom.iniP(57).pos(1:3) = [  10.79036d0, -11.28117d0, -22.36645d0 ]; geom.iniP(58).pos(1:3) = [ -16.54235d0,   5.42707d0, -20.99647d0 ]
    geom.iniP(59).pos(1:3) = [ -13.20762d0,  -3.44209d0, -24.19306d0 ]; geom.iniP(60).pos(1:3) = [  10.49633d0,   5.42707d0, -24.58309d0 ]
    geom.iniP(61).pos(1:3) = [   6.44354d0,  -3.44209d0, -26.79975d0 ]; geom.iniP(62).pos(1:3) = [  -3.55493d0,  -3.61805d0, -26.79975d0 ]

    ! Set face connnectivity
    geom.face(  1).n_poi = 3; allocate(geom.face(  1).poi(3)); geom.face(  1).poi(1:3) = [  3,  1,  2 ]
    geom.face(  2).n_poi = 3; allocate(geom.face(  2).poi(3)); geom.face(  2).poi(1:3) = [  3,  4,  1 ]
    geom.face(  3).n_poi = 3; allocate(geom.face(  3).poi(3)); geom.face(  3).poi(1:3) = [  1,  5,  2 ]
    geom.face(  4).n_poi = 3; allocate(geom.face(  4).poi(3)); geom.face(  4).poi(1:3) = [  2,  6,  3 ]
    geom.face(  5).n_poi = 3; allocate(geom.face(  5).poi(3)); geom.face(  5).poi(1:3) = [  1,  4,  5 ]
    geom.face(  6).n_poi = 3; allocate(geom.face(  6).poi(3)); geom.face(  6).poi(1:3) = [  4,  3,  7 ]
    geom.face(  7).n_poi = 3; allocate(geom.face(  7).poi(3)); geom.face(  7).poi(1:3) = [  2,  5,  8 ]
    geom.face(  8).n_poi = 3; allocate(geom.face(  8).poi(3)); geom.face(  8).poi(1:3) = [  3,  6,  9 ]
    geom.face(  9).n_poi = 3; allocate(geom.face(  9).poi(3)); geom.face(  9).poi(1:3) = [  6,  2, 10 ]
    geom.face( 10).n_poi = 3; allocate(geom.face( 10).poi(3)); geom.face( 10).poi(1:3) = [  4, 11,  5 ]
    geom.face( 11).n_poi = 3; allocate(geom.face( 11).poi(3)); geom.face( 11).poi(1:3) = [  3, 12,  7 ]
    geom.face( 12).n_poi = 3; allocate(geom.face( 12).poi(3)); geom.face( 12).poi(1:3) = [  7, 13,  4 ]
    geom.face( 13).n_poi = 3; allocate(geom.face( 13).poi(3)); geom.face( 13).poi(1:3) = [  5, 14,  8 ]
    geom.face( 14).n_poi = 3; allocate(geom.face( 14).poi(3)); geom.face( 14).poi(1:3) = [  8, 10,  2 ]
    geom.face( 15).n_poi = 3; allocate(geom.face( 15).poi(3)); geom.face( 15).poi(1:3) = [  6, 10,  9 ]
    geom.face( 16).n_poi = 3; allocate(geom.face( 16).poi(3)); geom.face( 16).poi(1:3) = [  9, 15,  3 ]
    geom.face( 17).n_poi = 3; allocate(geom.face( 17).poi(3)); geom.face( 17).poi(1:3) = [  5, 11, 16 ]
    geom.face( 18).n_poi = 3; allocate(geom.face( 18).poi(3)); geom.face( 18).poi(1:3) = [ 11,  4, 13 ]
    geom.face( 19).n_poi = 3; allocate(geom.face( 19).poi(3)); geom.face( 19).poi(1:3) = [  7, 12, 13 ]
    geom.face( 20).n_poi = 3; allocate(geom.face( 20).poi(3)); geom.face( 20).poi(1:3) = [ 12,  3, 17 ]
    geom.face( 21).n_poi = 3; allocate(geom.face( 21).poi(3)); geom.face( 21).poi(1:3) = [  8, 14, 10 ]
    geom.face( 22).n_poi = 3; allocate(geom.face( 22).poi(3)); geom.face( 22).poi(1:3) = [ 14,  5, 18 ]
    geom.face( 23).n_poi = 3; allocate(geom.face( 23).poi(3)); geom.face( 23).poi(1:3) = [  9, 10, 19 ]
    geom.face( 24).n_poi = 3; allocate(geom.face( 24).poi(3)); geom.face( 24).poi(1:3) = [  3, 15, 20 ]
    geom.face( 25).n_poi = 3; allocate(geom.face( 25).poi(3)); geom.face( 25).poi(1:3) = [ 15,  9, 21 ]
    geom.face( 26).n_poi = 3; allocate(geom.face( 26).poi(3)); geom.face( 26).poi(1:3) = [ 11, 13, 16 ]
    geom.face( 27).n_poi = 3; allocate(geom.face( 27).poi(3)); geom.face( 27).poi(1:3) = [ 16, 22,  5 ]
    geom.face( 28).n_poi = 3; allocate(geom.face( 28).poi(3)); geom.face( 28).poi(1:3) = [ 12, 23, 13 ]
    geom.face( 29).n_poi = 3; allocate(geom.face( 29).poi(3)); geom.face( 29).poi(1:3) = [  3, 20, 17 ]
    geom.face( 30).n_poi = 3; allocate(geom.face( 30).poi(3)); geom.face( 30).poi(1:3) = [ 17, 24, 12 ]
    geom.face( 31).n_poi = 3; allocate(geom.face( 31).poi(3)); geom.face( 31).poi(1:3) = [ 14, 25, 10 ]
    geom.face( 32).n_poi = 3; allocate(geom.face( 32).poi(3)); geom.face( 32).poi(1:3) = [  5, 26, 18 ]
    geom.face( 33).n_poi = 3; allocate(geom.face( 33).poi(3)); geom.face( 33).poi(1:3) = [ 18, 27, 14 ]
    geom.face( 34).n_poi = 3; allocate(geom.face( 34).poi(3)); geom.face( 34).poi(1:3) = [ 10, 28, 19 ]
    geom.face( 35).n_poi = 3; allocate(geom.face( 35).poi(3)); geom.face( 35).poi(1:3) = [ 19, 21,  9 ]
    geom.face( 36).n_poi = 3; allocate(geom.face( 36).poi(3)); geom.face( 36).poi(1:3) = [ 15, 21, 20 ]
    geom.face( 37).n_poi = 3; allocate(geom.face( 37).poi(3)); geom.face( 37).poi(1:3) = [ 16, 13, 29 ]
    geom.face( 38).n_poi = 3; allocate(geom.face( 38).poi(3)); geom.face( 38).poi(1:3) = [  5, 22, 26 ]
    geom.face( 39).n_poi = 3; allocate(geom.face( 39).poi(3)); geom.face( 39).poi(1:3) = [ 22, 16, 30 ]
    geom.face( 40).n_poi = 3; allocate(geom.face( 40).poi(3)); geom.face( 40).poi(1:3) = [ 13, 23, 31 ]
    geom.face( 41).n_poi = 3; allocate(geom.face( 41).poi(3)); geom.face( 41).poi(1:3) = [ 23, 12, 24 ]
    geom.face( 42).n_poi = 3; allocate(geom.face( 42).poi(3)); geom.face( 42).poi(1:3) = [ 17, 20, 24 ]
    geom.face( 43).n_poi = 3; allocate(geom.face( 43).poi(3)); geom.face( 43).poi(1:3) = [ 10, 25, 32 ]
    geom.face( 44).n_poi = 3; allocate(geom.face( 44).poi(3)); geom.face( 44).poi(1:3) = [ 25, 14, 27 ]
    geom.face( 45).n_poi = 3; allocate(geom.face( 45).poi(3)); geom.face( 45).poi(1:3) = [ 18, 26, 27 ]
    geom.face( 46).n_poi = 3; allocate(geom.face( 46).poi(3)); geom.face( 46).poi(1:3) = [ 19, 28, 21 ]
    geom.face( 47).n_poi = 3; allocate(geom.face( 47).poi(3)); geom.face( 47).poi(1:3) = [ 28, 10, 33 ]
    geom.face( 48).n_poi = 3; allocate(geom.face( 48).poi(3)); geom.face( 48).poi(1:3) = [ 20, 21, 34 ]
    geom.face( 49).n_poi = 3; allocate(geom.face( 49).poi(3)); geom.face( 49).poi(1:3) = [ 13, 35, 29 ]
    geom.face( 50).n_poi = 3; allocate(geom.face( 50).poi(3)); geom.face( 50).poi(1:3) = [ 29, 30, 16 ]
    geom.face( 51).n_poi = 3; allocate(geom.face( 51).poi(3)); geom.face( 51).poi(1:3) = [ 22, 30, 26 ]
    geom.face( 52).n_poi = 3; allocate(geom.face( 52).poi(3)); geom.face( 52).poi(1:3) = [ 23, 24, 31 ]
    geom.face( 53).n_poi = 3; allocate(geom.face( 53).poi(3)); geom.face( 53).poi(1:3) = [ 31, 36, 13 ]
    geom.face( 54).n_poi = 3; allocate(geom.face( 54).poi(3)); geom.face( 54).poi(1:3) = [ 20, 34, 24 ]
    geom.face( 55).n_poi = 3; allocate(geom.face( 55).poi(3)); geom.face( 55).poi(1:3) = [ 25, 27, 32 ]
    geom.face( 56).n_poi = 3; allocate(geom.face( 56).poi(3)); geom.face( 56).poi(1:3) = [ 32, 33, 10 ]
    geom.face( 57).n_poi = 3; allocate(geom.face( 57).poi(3)); geom.face( 57).poi(1:3) = [ 26, 37, 27 ]
    geom.face( 58).n_poi = 3; allocate(geom.face( 58).poi(3)); geom.face( 58).poi(1:3) = [ 28, 38, 21 ]
    geom.face( 59).n_poi = 3; allocate(geom.face( 59).poi(3)); geom.face( 59).poi(1:3) = [ 33, 39, 28 ]
    geom.face( 60).n_poi = 3; allocate(geom.face( 60).poi(3)); geom.face( 60).poi(1:3) = [ 21, 40, 34 ]
    geom.face( 61).n_poi = 3; allocate(geom.face( 61).poi(3)); geom.face( 61).poi(1:3) = [ 29, 35, 30 ]
    geom.face( 62).n_poi = 3; allocate(geom.face( 62).poi(3)); geom.face( 62).poi(1:3) = [ 35, 13, 36 ]
    geom.face( 63).n_poi = 3; allocate(geom.face( 63).poi(3)); geom.face( 63).poi(1:3) = [ 26, 30, 37 ]
    geom.face( 64).n_poi = 3; allocate(geom.face( 64).poi(3)); geom.face( 64).poi(1:3) = [ 31, 24, 41 ]
    geom.face( 65).n_poi = 3; allocate(geom.face( 65).poi(3)); geom.face( 65).poi(1:3) = [ 36, 31, 42 ]
    geom.face( 66).n_poi = 3; allocate(geom.face( 66).poi(3)); geom.face( 66).poi(1:3) = [ 24, 34, 40 ]
    geom.face( 67).n_poi = 3; allocate(geom.face( 67).poi(3)); geom.face( 67).poi(1:3) = [ 32, 27, 43 ]
    geom.face( 68).n_poi = 3; allocate(geom.face( 68).poi(3)); geom.face( 68).poi(1:3) = [ 33, 32, 39 ]
    geom.face( 69).n_poi = 3; allocate(geom.face( 69).poi(3)); geom.face( 69).poi(1:3) = [ 27, 37, 44 ]
    geom.face( 70).n_poi = 3; allocate(geom.face( 70).poi(3)); geom.face( 70).poi(1:3) = [ 21, 38, 45 ]
    geom.face( 71).n_poi = 3; allocate(geom.face( 71).poi(3)); geom.face( 71).poi(1:3) = [ 38, 28, 39 ]
    geom.face( 72).n_poi = 3; allocate(geom.face( 72).poi(3)); geom.face( 72).poi(1:3) = [ 40, 21, 46 ]
    geom.face( 73).n_poi = 3; allocate(geom.face( 73).poi(3)); geom.face( 73).poi(1:3) = [ 35, 47, 30 ]
    geom.face( 74).n_poi = 3; allocate(geom.face( 74).poi(3)); geom.face( 74).poi(1:3) = [ 36, 42, 35 ]
    geom.face( 75).n_poi = 3; allocate(geom.face( 75).poi(3)); geom.face( 75).poi(1:3) = [ 30, 44, 37 ]
    geom.face( 76).n_poi = 3; allocate(geom.face( 76).poi(3)); geom.face( 76).poi(1:3) = [ 24, 48, 41 ]
    geom.face( 77).n_poi = 3; allocate(geom.face( 77).poi(3)); geom.face( 77).poi(1:3) = [ 41, 42, 31 ]
    geom.face( 78).n_poi = 3; allocate(geom.face( 78).poi(3)); geom.face( 78).poi(1:3) = [ 40, 49, 24 ]
    geom.face( 79).n_poi = 3; allocate(geom.face( 79).poi(3)); geom.face( 79).poi(1:3) = [ 27, 50, 43 ]
    geom.face( 80).n_poi = 3; allocate(geom.face( 80).poi(3)); geom.face( 80).poi(1:3) = [ 43, 39, 32 ]
    geom.face( 81).n_poi = 3; allocate(geom.face( 81).poi(3)); geom.face( 81).poi(1:3) = [ 44, 51, 27 ]
    geom.face( 82).n_poi = 3; allocate(geom.face( 82).poi(3)); geom.face( 82).poi(1:3) = [ 38, 39, 45 ]
    geom.face( 83).n_poi = 3; allocate(geom.face( 83).poi(3)); geom.face( 83).poi(1:3) = [ 45, 46, 21 ]
    geom.face( 84).n_poi = 3; allocate(geom.face( 84).poi(3)); geom.face( 84).poi(1:3) = [ 46, 52, 40 ]
    geom.face( 85).n_poi = 3; allocate(geom.face( 85).poi(3)); geom.face( 85).poi(1:3) = [ 30, 47, 53 ]
    geom.face( 86).n_poi = 3; allocate(geom.face( 86).poi(3)); geom.face( 86).poi(1:3) = [ 47, 35, 42 ]
    geom.face( 87).n_poi = 3; allocate(geom.face( 87).poi(3)); geom.face( 87).poi(1:3) = [ 44, 30, 54 ]
    geom.face( 88).n_poi = 3; allocate(geom.face( 88).poi(3)); geom.face( 88).poi(1:3) = [ 41, 48, 42 ]
    geom.face( 89).n_poi = 3; allocate(geom.face( 89).poi(3)); geom.face( 89).poi(1:3) = [ 48, 24, 49 ]
    geom.face( 90).n_poi = 3; allocate(geom.face( 90).poi(3)); geom.face( 90).poi(1:3) = [ 49, 40, 52 ]
    geom.face( 91).n_poi = 3; allocate(geom.face( 91).poi(3)); geom.face( 91).poi(1:3) = [ 43, 50, 39 ]
    geom.face( 92).n_poi = 3; allocate(geom.face( 92).poi(3)); geom.face( 92).poi(1:3) = [ 50, 27, 51 ]
    geom.face( 93).n_poi = 3; allocate(geom.face( 93).poi(3)); geom.face( 93).poi(1:3) = [ 51, 44, 55 ]
    geom.face( 94).n_poi = 3; allocate(geom.face( 94).poi(3)); geom.face( 94).poi(1:3) = [ 45, 39, 56 ]
    geom.face( 95).n_poi = 3; allocate(geom.face( 95).poi(3)); geom.face( 95).poi(1:3) = [ 46, 45, 52 ]
    geom.face( 96).n_poi = 3; allocate(geom.face( 96).poi(3)); geom.face( 96).poi(1:3) = [ 47, 42, 53 ]
    geom.face( 97).n_poi = 3; allocate(geom.face( 97).poi(3)); geom.face( 97).poi(1:3) = [ 53, 54, 30 ]
    geom.face( 98).n_poi = 3; allocate(geom.face( 98).poi(3)); geom.face( 98).poi(1:3) = [ 54, 55, 44 ]
    geom.face( 99).n_poi = 3; allocate(geom.face( 99).poi(3)); geom.face( 99).poi(1:3) = [ 48, 57, 42 ]
    geom.face(100).n_poi = 3; allocate(geom.face(100).poi(3)); geom.face(100).poi(1:3) = [ 49, 52, 48 ]
    geom.face(101).n_poi = 3; allocate(geom.face(101).poi(3)); geom.face(101).poi(1:3) = [ 50, 58, 39 ]
    geom.face(102).n_poi = 3; allocate(geom.face(102).poi(3)); geom.face(102).poi(1:3) = [ 51, 55, 50 ]
    geom.face(103).n_poi = 3; allocate(geom.face(103).poi(3)); geom.face(103).poi(1:3) = [ 39, 59, 56 ]
    geom.face(104).n_poi = 3; allocate(geom.face(104).poi(3)); geom.face(104).poi(1:3) = [ 56, 52, 45 ]
    geom.face(105).n_poi = 3; allocate(geom.face(105).poi(3)); geom.face(105).poi(1:3) = [ 53, 42, 60 ]
    geom.face(106).n_poi = 3; allocate(geom.face(106).poi(3)); geom.face(106).poi(1:3) = [ 54, 53, 55 ]
    geom.face(107).n_poi = 3; allocate(geom.face(107).poi(3)); geom.face(107).poi(1:3) = [ 42, 57, 61 ]
    geom.face(108).n_poi = 3; allocate(geom.face(108).poi(3)); geom.face(108).poi(1:3) = [ 57, 48, 52 ]
    geom.face(109).n_poi = 3; allocate(geom.face(109).poi(3)); geom.face(109).poi(1:3) = [ 39, 58, 59 ]
    geom.face(110).n_poi = 3; allocate(geom.face(110).poi(3)); geom.face(110).poi(1:3) = [ 58, 50, 55 ]
    geom.face(111).n_poi = 3; allocate(geom.face(111).poi(3)); geom.face(111).poi(1:3) = [ 56, 59, 52 ]
    geom.face(112).n_poi = 3; allocate(geom.face(112).poi(3)); geom.face(112).poi(1:3) = [ 42, 61, 60 ]
    geom.face(113).n_poi = 3; allocate(geom.face(113).poi(3)); geom.face(113).poi(1:3) = [ 60, 55, 53 ]
    geom.face(114).n_poi = 3; allocate(geom.face(114).poi(3)); geom.face(114).poi(1:3) = [ 57, 52, 61 ]
    geom.face(115).n_poi = 3; allocate(geom.face(115).poi(3)); geom.face(115).poi(1:3) = [ 58, 55, 59 ]
    geom.face(116).n_poi = 3; allocate(geom.face(116).poi(3)); geom.face(116).poi(1:3) = [ 59, 62, 52 ]
    geom.face(117).n_poi = 3; allocate(geom.face(117).poi(3)); geom.face(117).poi(1:3) = [ 60, 61, 55 ]
    geom.face(118).n_poi = 3; allocate(geom.face(118).poi(3)); geom.face(118).poi(1:3) = [ 61, 52, 62 ]
    geom.face(119).n_poi = 3; allocate(geom.face(119).poi(3)); geom.face(119).poi(1:3) = [ 59, 55, 62 ]
    geom.face(120).n_poi = 3; allocate(geom.face(120).poi(3)); geom.face(120).poi(1:3) = [ 61, 62, 55 ]
end subroutine Exam_Catalan_Disdyakis_Triacontahedron

! ---------------------------------------------------------------------------------------

! Example of Deltoidal Hexecontahedron (Trapezoidal Hexecontahedron, Strombic Hexecontahedron, Tetragonal Hexacontahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Deltoidal_Hexecontahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "57_Deltoi_Hexecontahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method
    prob.name_prob = "Deltoidal Hexecontahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "YZ"

    ! Allocate point and face structure
    geom.n_iniP = 62
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   6.05353d0,   6.38197d0,  26.34802d0 ]; geom.iniP( 2).pos(1:3) = [  -3.93254d0,   5.85410d0,  26.34802d0 ]
    geom.iniP( 3).pos(1:3) = [  -8.19710d0,  -3.19099d0,  26.34802d0 ]; geom.iniP( 4).pos(1:3) = [   6.05353d0,  -9.01148d0,  26.34802d0 ]
    geom.iniP( 5).pos(1:3) = [  18.87512d0,  -3.19099d0,  20.12809d0 ]; geom.iniP( 6).pos(1:3) = [  15.03820d0,   5.85410d0,  21.98944d0 ]
    geom.iniP( 7).pos(1:3) = [   4.51203d0,  20.15031d0,  19.63865d0 ]; geom.iniP( 8).pos(1:3) = [  -9.08029d0,  13.51724d0,  22.50388d0 ]
    geom.iniP( 9).pos(1:3) = [ -20.36039d0,   3.44208d0,  19.63865d0 ]; geom.iniP(10).pos(1:3) = [ -11.10565d0, -11.70823d0,  21.98944d0 ]
    geom.iniP(11).pos(1:3) = [  -4.18288d0, -18.68036d0,  20.12809d0 ]; geom.iniP(12).pos(1:3) = [  17.99193d0,  13.51724d0,  16.28397d0 ]
    geom.iniP(13).pos(1:3) = [  26.89027d0,   3.44208d0,   8.78266d0 ]; geom.iniP(14).pos(1:3) = [  12.54867d0, -18.68036d0,  16.28397d0 ]
    geom.iniP(15).pos(1:3) = [  19.58964d0, -11.70823d0,  14.93709d0 ]; geom.iniP(16).pos(1:3) = [  -6.49514d0,  25.06234d0,  10.06404d0 ]
    geom.iniP(17).pos(1:3) = [ -12.72597d0,  18.94430d0,  14.93709d0 ]; geom.iniP(18).pos(1:3) = [ -19.31673d0, -11.54510d0,  16.28397d0 ]
    geom.iniP(19).pos(1:3) = [ -20.74577d0,  15.48938d0,  10.06404d0 ]; geom.iniP(20).pos(1:3) = [   3.43184d0, -22.56232d0,  14.93709d0 ]
    geom.iniP(21).pos(1:3) = [ -13.35417d0, -23.59239d0,   8.78266d0 ]; geom.iniP(22).pos(1:3) = [  17.96936d0,  18.94430d0,   7.88477d0 ]
    geom.iniP(23).pos(1:3) = [  10.23641d0,  25.06234d0,   6.21993d0 ]; geom.iniP(24).pos(1:3) = [  24.48704d0, -11.54510d0,   6.21993d0 ]
    geom.iniP(25).pos(1:3) = [  23.05800d0,  15.48938d0,   0.00000d0 ]; geom.iniP(26).pos(1:3) = [  15.84834d0, -23.59239d0,   2.07331d0 ]
    geom.iniP(27).pos(1:3) = [   0.81015d0,  27.03447d0,   3.52616d0 ]; geom.iniP(28).pos(1:3) = [ -15.84834d0,  23.59239d0,  -2.07331d0 ]
    geom.iniP(29).pos(1:3) = [ -27.07232d0,   0.00000d0,   6.21993d0 ]; geom.iniP(30).pos(1:3) = [ -24.33232d0,  -9.47215d0,   7.88477d0 ]
    geom.iniP(31).pos(1:3) = [ -25.33374d0,   9.47215d0,   3.52616d0 ]; geom.iniP(32).pos(1:3) = [   1.42905d0, -27.03447d0,   6.21993d0 ]
    geom.iniP(33).pos(1:3) = [ -23.05800d0, -15.48938d0,   0.00000d0 ]; geom.iniP(34).pos(1:3) = [  13.35417d0,  23.59239d0,  -8.78266d0 ]
    geom.iniP(35).pos(1:3) = [  25.33374d0,  -9.47215d0,  -3.52616d0 ]; geom.iniP(36).pos(1:3) = [  27.07232d0,   0.00000d0,  -6.21993d0 ]
    geom.iniP(37).pos(1:3) = [  24.33232d0,   9.47215d0,  -7.88477d0 ]; geom.iniP(38).pos(1:3) = [  20.74577d0, -15.48938d0, -10.06404d0 ]
    geom.iniP(39).pos(1:3) = [  -1.42905d0,  27.03447d0,  -6.21993d0 ]; geom.iniP(40).pos(1:3) = [ -24.48704d0,  11.54510d0,  -6.21993d0 ]
    geom.iniP(41).pos(1:3) = [ -26.89027d0,  -3.44208d0,  -8.78266d0 ]; geom.iniP(42).pos(1:3) = [ -10.23641d0, -25.06234d0,  -6.21993d0 ]
    geom.iniP(43).pos(1:3) = [  -0.81015d0, -27.03447d0,  -3.52616d0 ]; geom.iniP(44).pos(1:3) = [   6.49514d0, -25.06234d0, -10.06404d0 ]
    geom.iniP(45).pos(1:3) = [ -17.96936d0, -18.94430d0,  -7.88477d0 ]; geom.iniP(46).pos(1:3) = [  19.31673d0,  11.54510d0, -16.28397d0 ]
    geom.iniP(47).pos(1:3) = [  20.36039d0,  -3.44208d0, -19.63865d0 ]; geom.iniP(48).pos(1:3) = [  12.72597d0, -18.94430d0, -14.93709d0 ]
    geom.iniP(49).pos(1:3) = [  -3.43184d0,  22.56232d0, -14.93709d0 ]; geom.iniP(50).pos(1:3) = [ -12.54867d0,  18.68036d0, -16.28397d0 ]
    geom.iniP(51).pos(1:3) = [   4.18288d0,  18.68036d0, -20.12809d0 ]; geom.iniP(52).pos(1:3) = [ -19.58964d0,  11.70823d0, -14.93709d0 ]
    geom.iniP(53).pos(1:3) = [ -17.99193d0, -13.51724d0, -16.28397d0 ]; geom.iniP(54).pos(1:3) = [ -18.87512d0,   3.19099d0, -20.12809d0 ]
    geom.iniP(55).pos(1:3) = [  -4.51203d0, -20.15031d0, -19.63865d0 ]; geom.iniP(56).pos(1:3) = [  11.10565d0,  11.70823d0, -21.98944d0 ]
    geom.iniP(57).pos(1:3) = [   9.08029d0, -13.51724d0, -22.50388d0 ]; geom.iniP(58).pos(1:3) = [   8.19710d0,   3.19099d0, -26.34802d0 ]
    geom.iniP(59).pos(1:3) = [  -6.05353d0,   9.01148d0, -26.34802d0 ]; geom.iniP(60).pos(1:3) = [ -15.03820d0,  -5.85410d0, -21.98944d0 ]
    geom.iniP(61).pos(1:3) = [  -6.05353d0,  -6.38197d0, -26.34802d0 ]; geom.iniP(62).pos(1:3) = [   3.93254d0,  -5.85410d0, -26.34802d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  1,  2,  3 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  4,  5,  6,  1 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  1,  7,  8,  2 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  2,  8,  9,  3 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  3, 10, 11,  4 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  1,  6, 12,  7 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [  6,  5, 13, 12 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [  5,  4, 14, 15 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  8,  7, 16, 17 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  3,  9, 18, 10 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  9,  8, 17, 19 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  4, 11, 20, 14 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [ 11, 10, 18, 21 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [ 12, 22, 23,  7 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  5, 15, 24, 13 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [ 13, 25, 22, 12 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [ 14, 26, 24, 15 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  7, 23, 27, 16 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [ 16, 28, 19, 17 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [  9, 29, 30, 18 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 19, 31, 29,  9 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 11, 21, 32, 20 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 20, 32, 26, 14 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 18, 30, 33, 21 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 23, 22, 25, 34 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 13, 24, 35, 36 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [ 25, 13, 36, 37 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [ 24, 26, 38, 35 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 16, 27, 39, 28 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 27, 23, 34, 39 ]
    geom.face(31).n_poi = 4; allocate(geom.face(31).poi(4)); geom.face(31).poi(1:4) = [ 19, 28, 40, 31 ]
    geom.face(32).n_poi = 4; allocate(geom.face(32).poi(4)); geom.face(32).poi(1:4) = [ 30, 29, 41, 33 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [ 29, 31, 40, 41 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [ 32, 21, 42, 43 ]
    geom.face(35).n_poi = 4; allocate(geom.face(35).poi(4)); geom.face(35).poi(1:4) = [ 26, 32, 43, 44 ]
    geom.face(36).n_poi = 4; allocate(geom.face(36).poi(4)); geom.face(36).poi(1:4) = [ 21, 33, 45, 42 ]
    geom.face(37).n_poi = 4; allocate(geom.face(37).poi(4)); geom.face(37).poi(1:4) = [ 25, 37, 46, 34 ]
    geom.face(38).n_poi = 4; allocate(geom.face(38).poi(4)); geom.face(38).poi(1:4) = [ 35, 38, 47, 36 ]
    geom.face(39).n_poi = 4; allocate(geom.face(39).poi(4)); geom.face(39).poi(1:4) = [ 36, 47, 46, 37 ]
    geom.face(40).n_poi = 4; allocate(geom.face(40).poi(4)); geom.face(40).poi(1:4) = [ 26, 44, 48, 38 ]
    geom.face(41).n_poi = 4; allocate(geom.face(41).poi(4)); geom.face(41).poi(1:4) = [ 39, 49, 50, 28 ]
    geom.face(42).n_poi = 4; allocate(geom.face(42).poi(4)); geom.face(42).poi(1:4) = [ 34, 51, 49, 39 ]
    geom.face(43).n_poi = 4; allocate(geom.face(43).poi(4)); geom.face(43).poi(1:4) = [ 28, 50, 52, 40 ]
    geom.face(44).n_poi = 4; allocate(geom.face(44).poi(4)); geom.face(44).poi(1:4) = [ 41, 53, 45, 33 ]
    geom.face(45).n_poi = 4; allocate(geom.face(45).poi(4)); geom.face(45).poi(1:4) = [ 40, 52, 54, 41 ]
    geom.face(46).n_poi = 4; allocate(geom.face(46).poi(4)); geom.face(46).poi(1:4) = [ 42, 55, 44, 43 ]
    geom.face(47).n_poi = 4; allocate(geom.face(47).poi(4)); geom.face(47).poi(1:4) = [ 45, 53, 55, 42 ]
    geom.face(48).n_poi = 4; allocate(geom.face(48).poi(4)); geom.face(48).poi(1:4) = [ 34, 46, 56, 51 ]
    geom.face(49).n_poi = 4; allocate(geom.face(49).poi(4)); geom.face(49).poi(1:4) = [ 47, 38, 48, 57 ]
    geom.face(50).n_poi = 4; allocate(geom.face(50).poi(4)); geom.face(50).poi(1:4) = [ 46, 47, 58, 56 ]
    geom.face(51).n_poi = 4; allocate(geom.face(51).poi(4)); geom.face(51).poi(1:4) = [ 48, 44, 55, 57 ]
    geom.face(52).n_poi = 4; allocate(geom.face(52).poi(4)); geom.face(52).poi(1:4) = [ 50, 49, 51, 59 ]
    geom.face(53).n_poi = 4; allocate(geom.face(53).poi(4)); geom.face(53).poi(1:4) = [ 52, 50, 59, 54 ]
    geom.face(54).n_poi = 4; allocate(geom.face(54).poi(4)); geom.face(54).poi(1:4) = [ 53, 41, 54, 60 ]
    geom.face(55).n_poi = 4; allocate(geom.face(55).poi(4)); geom.face(55).poi(1:4) = [ 55, 53, 60, 61 ]
    geom.face(56).n_poi = 4; allocate(geom.face(56).poi(4)); geom.face(56).poi(1:4) = [ 56, 58, 59, 51 ]
    geom.face(57).n_poi = 4; allocate(geom.face(57).poi(4)); geom.face(57).poi(1:4) = [ 57, 62, 58, 47 ]
    geom.face(58).n_poi = 4; allocate(geom.face(58).poi(4)); geom.face(58).poi(1:4) = [ 55, 61, 62, 57 ]
    geom.face(59).n_poi = 4; allocate(geom.face(59).poi(4)); geom.face(59).poi(1:4) = [ 59, 61, 60, 54 ]
    geom.face(60).n_poi = 4; allocate(geom.face(60).poi(4)); geom.face(60).poi(1:4) = [ 59, 58, 62, 61 ]
end subroutine Exam_Catalan_Deltoidal_Hexecontahedron

! ---------------------------------------------------------------------------------------

! Example of Pentagonal Hexecontahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Catalan_Pentagonal_Hexecontahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "58_Penta_Hexecontahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Pentagonal Hexecontahedron"

    ! Set geometric type and view
    prob.color    = [247, 147, 30]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 3.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 92
    geom.n_face = 60

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   8.34392d0,   5.00003d0,  34.99539d0 ]; geom.iniP( 2).pos(1:3) = [  -0.47434d0,   9.71578d0,  34.99539d0 ]
    geom.iniP( 3).pos(1:3) = [  -8.79133d0,   4.16344d0,  34.99539d0 ]; geom.iniP( 4).pos(1:3) = [  -7.81719d0,  -5.78902d0,  34.99539d0 ]
    geom.iniP( 5).pos(1:3) = [   8.34392d0, -12.49858d0,  34.99539d0 ]; geom.iniP( 6).pos(1:3) = [  22.76645d0,  -5.78902d0,  27.70334d0 ]
    geom.iniP( 7).pos(1:3) = [  23.63575d0,   4.16344d0,  27.26382d0 ]; geom.iniP( 8).pos(1:3) = [  16.21354d0,   9.71578d0,  31.01651d0 ]
    geom.iniP( 9).pos(1:3) = [  14.81747d0,  18.87924d0,  27.26382d0 ]; geom.iniP(10).pos(1:3) = [  -1.25719d0,  25.75047d0,  28.03292d0 ]
    geom.iniP(11).pos(1:3) = [ -16.58607d0,  17.34607d0,  27.26382d0 ]; geom.iniP(12).pos(1:3) = [ -17.08287d0,   8.09019d0,  31.01651d0 ]
    geom.iniP(13).pos(1:3) = [ -29.06245d0,  -2.84454d0,  24.44985d0 ]; geom.iniP(14).pos(1:3) = [ -15.19001d0, -11.24898d0,  31.01651d0 ]
    geom.iniP(15).pos(1:3) = [ -12.90796d0, -20.23285d0,  27.26382d0 ]; geom.iniP(16).pos(1:3) = [  -3.38283d0, -23.24608d0,  27.70334d0 ]
    geom.iniP(17).pos(1:3) = [  15.51890d0, -23.24608d0,  23.19658d0 ]; geom.iniP(18).pos(1:3) = [  23.82095d0, -20.23285d0,  18.50652d0 ]
    geom.iniP(19).pos(1:3) = [  27.55076d0, -11.24898d0,  20.82581d0 ]; geom.iniP(20).pos(1:3) = [  32.93240d0,  -6.44602d0,  13.89975d0 ]
    geom.iniP(21).pos(1:3) = [  33.44232d0,  11.03467d0,  14.50380d0 ]; geom.iniP(22).pos(1:3) = [  20.92299d0,  21.96943d0,  19.97177d0 ]
    geom.iniP(23).pos(1:3) = [  17.49512d0,  28.81069d0,  13.53377d0 ]; geom.iniP(24).pos(1:3) = [   8.15663d0,  32.17287d0,  14.75383d0 ]
    geom.iniP(25).pos(1:3) = [ -11.25182d0,  31.22528d0,  14.75383d0 ]; geom.iniP(26).pos(1:3) = [ -20.21842d0,  26.96943d0,  13.53377d0 ]
    geom.iniP(27).pos(1:3) = [ -22.96356d0,  19.82679d0,  19.97177d0 ]; geom.iniP(28).pos(1:3) = [ -29.47532d0,  12.91062d0,  16.84694d0 ]
    geom.iniP(29).pos(1:3) = [ -26.41262d0, -18.38078d0,  16.84694d0 ]; geom.iniP(30).pos(1:3) = [ -18.68342d0, -23.90308d0,  19.97177d0 ]
    geom.iniP(31).pos(1:3) = [ -11.54755d0, -35.23304d0,   8.70624d0 ]; geom.iniP(32).pos(1:3) = [  -0.17468d0, -29.75824d0,  20.82581d0 ]
    geom.iniP(33).pos(1:3) = [   9.55269d0, -29.75824d0,  18.50652d0 ]; geom.iniP(34).pos(1:3) = [  12.22772d0, -32.88695d0,   9.39303d0 ]
    geom.iniP(35).pos(1:3) = [  27.08250d0, -26.65520d0,   2.55926d0 ]; geom.iniP(36).pos(1:3) = [  34.27834d0, -10.90004d0,   5.04817d0 ]
    geom.iniP(37).pos(1:3) = [  35.83914d0,  -5.11098d0,  -2.95507d0 ]; geom.iniP(38).pos(1:3) = [  35.96527d0,   4.80292d0,  -1.65170d0 ]
    geom.iniP(39).pos(1:3) = [  28.54306d0,  22.36463d0,   2.10099d0 ]; geom.iniP(40).pos(1:3) = [  21.41665d0,  29.01407d0,   4.33699d0 ]
    geom.iniP(41).pos(1:3) = [  11.54755d0,  35.23305d0,  -8.70624d0 ]; geom.iniP(42).pos(1:3) = [   3.27049d0,  35.54727d0,   6.70776d0 ]
    geom.iniP(43).pos(1:3) = [  -6.71765d0,  35.05964d0,   6.70776d0 ]; geom.iniP(44).pos(1:3) = [ -11.40781d0,  34.42013d0,  -2.10099d0 ]
    geom.iniP(45).pos(1:3) = [ -27.08250d0,  26.65521d0,  -2.55926d0 ]; geom.iniP(46).pos(1:3) = [ -32.87174d0,  13.53028d0,   7.46178d0 ]
    geom.iniP(47).pos(1:3) = [ -35.83914d0,   5.11098d0,   2.95507d0 ]; geom.iniP(48).pos(1:3) = [ -35.24143d0,  -3.44932d0,   8.08968d0 ]
    geom.iniP(49).pos(1:3) = [ -34.16797d0, -11.96307d0,   2.95507d0 ]; geom.iniP(50).pos(1:3) = [ -29.62444d0, -19.64722d0,   7.46178d0 ]
    geom.iniP(51).pos(1:3) = [ -24.92445d0, -26.36397d0,   1.73500d0 ]; geom.iniP(52).pos(1:3) = [   5.02331d0, -35.83781d0,   3.11692d0 ]
    geom.iniP(53).pos(1:3) = [   6.71765d0, -35.05964d0,  -6.70776d0 ]; geom.iniP(54).pos(1:3) = [  15.52006d0, -31.37480d0,  -9.69779d0 ]
    geom.iniP(55).pos(1:3) = [  29.14808d0, -17.78613d0, -12.38306d0 ]; geom.iniP(56).pos(1:3) = [  33.19898d0,  -8.65479d0, -11.92566d0 ]
    geom.iniP(57).pos(1:3) = [  29.06246d0,   2.84455d0, -24.44985d0 ]; geom.iniP(58).pos(1:3) = [  33.44409d0,  10.60947d0,  -9.39303d0 ]
    geom.iniP(59).pos(1:3) = [  29.62444d0,  19.64723d0,  -7.46178d0 ]; geom.iniP(60).pos(1:3) = [  23.51789d0,  23.73371d0, -14.24495d0 ]
    geom.iniP(61).pos(1:3) = [  -5.84321d0,  34.30464d0, -10.40897d0 ]; geom.iniP(62).pos(1:3) = [  -9.55269d0,  29.75824d0, -18.50652d0 ]
    geom.iniP(63).pos(1:3) = [ -18.61586d0,  25.58577d0, -17.83578d0 ]; geom.iniP(64).pos(1:3) = [ -31.88131d0,  12.67516d0, -11.92566d0 ]
    geom.iniP(65).pos(1:3) = [ -35.32943d0,   4.67092d0,  -7.02226d0 ]; geom.iniP(66).pos(1:3) = [ -33.44232d0, -11.03467d0, -14.50380d0 ]
    geom.iniP(67).pos(1:3) = [ -25.03523d0, -25.01475d0,  -8.17296d0 ]; geom.iniP(68).pos(1:3) = [ -17.49512d0, -28.81069d0, -13.53376d0 ]
    geom.iniP(69).pos(1:3) = [ -10.27290d0, -33.74006d0,  -8.68185d0 ]; geom.iniP(70).pos(1:3) = [  -1.15421d0, -33.98006d0, -12.77974d0 ]
    geom.iniP(71).pos(1:3) = [   1.25719d0, -25.75047d0, -28.03292d0 ]; geom.iniP(72).pos(1:3) = [  15.95019d0, -26.81990d0, -18.58985d0 ]
    geom.iniP(73).pos(1:3) = [  22.96356d0, -19.82678d0, -19.97177d0 ]; geom.iniP(74).pos(1:3) = [  21.18151d0, -12.62007d0, -26.67161d0 ]
    geom.iniP(75).pos(1:3) = [  21.57815d0,  18.55013d0, -22.57372d0 ]; geom.iniP(76).pos(1:3) = [  12.90796d0,  20.23285d0, -27.26382d0 ]
    geom.iniP(77).pos(1:3) = [   6.67038d0,  27.00346d0, -23.35848d0 ]; geom.iniP(78).pos(1:3) = [  -3.11285d0,  26.00083d0, -25.17071d0 ]
    geom.iniP(79).pos(1:3) = [  -8.34392d0,  12.49858d0, -34.99539d0 ]; geom.iniP(80).pos(1:3) = [ -20.72400d0,  17.89313d0, -23.86736d0 ]
    geom.iniP(81).pos(1:3) = [ -27.55076d0,  11.24898d0, -20.82581d0 ]; geom.iniP(82).pos(1:3) = [ -26.91452d0,   1.89963d0, -24.31663d0 ]
    geom.iniP(83).pos(1:3) = [ -20.55236d0, -16.44713d0, -25.02781d0 ]; geom.iniP(84).pos(1:3) = [ -15.18810d0, -24.40154d0, -22.20774d0 ]
    geom.iniP(85).pos(1:3) = [  12.48739d0, -12.81618d0, -31.60868d0 ]; geom.iniP(86).pos(1:3) = [   8.79133d0,  -4.16343d0, -34.99539d0 ]
    geom.iniP(87).pos(1:3) = [  13.99944d0,   4.19348d0, -33.25248d0 ]; geom.iniP(88).pos(1:3) = [   9.00776d0,  12.84450d0, -32.75938d0 ]
    geom.iniP(89).pos(1:3) = [ -19.48767d0,  -0.27403d0, -30.65052d0 ]; geom.iniP(90).pos(1:3) = [ -16.21354d0,  -9.71577d0, -31.01651d0 ]
    geom.iniP(91).pos(1:3) = [  -6.75715d0, -11.32154d0, -33.84468d0 ]; geom.iniP(92).pos(1:3) = [  -1.11246d0,  -3.39426d0, -36.14616d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [  5,  1,  2,  3,  4 ]
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [  1,  5,  6,  7,  8 ]
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [  2,  1,  8,  9, 10 ]
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [  3,  2, 10, 11, 12 ]
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [  4,  3, 12, 13, 14 ]
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [  5,  4, 14, 15, 16 ]
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [  6,  5, 17, 18, 19 ]
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [  7,  6, 19, 20, 21 ]
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [  8,  7, 21, 22,  9 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 10,  9, 22, 23, 24 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 11, 10, 25, 26, 27 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 12, 11, 27, 28, 13 ]
    geom.face(13).n_poi = 5; allocate(geom.face(13).poi(5)); geom.face(13).poi(1:5) = [ 14, 13, 29, 30, 15 ]
    geom.face(14).n_poi = 5; allocate(geom.face(14).poi(5)); geom.face(14).poi(1:5) = [ 16, 15, 30, 31, 32 ]
    geom.face(15).n_poi = 5; allocate(geom.face(15).poi(5)); geom.face(15).poi(1:5) = [  5, 16, 32, 33, 17 ]
    geom.face(16).n_poi = 5; allocate(geom.face(16).poi(5)); geom.face(16).poi(1:5) = [ 18, 17, 33, 34, 35 ]
    geom.face(17).n_poi = 5; allocate(geom.face(17).poi(5)); geom.face(17).poi(1:5) = [ 19, 18, 35, 36, 20 ]
    geom.face(18).n_poi = 5; allocate(geom.face(18).poi(5)); geom.face(18).poi(1:5) = [ 21, 20, 36, 37, 38 ]
    geom.face(19).n_poi = 5; allocate(geom.face(19).poi(5)); geom.face(19).poi(1:5) = [ 22, 21, 39, 40, 23 ]
    geom.face(20).n_poi = 5; allocate(geom.face(20).poi(5)); geom.face(20).poi(1:5) = [ 24, 23, 40, 41, 42 ]
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [ 10, 24, 42, 43, 25 ]
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [ 26, 25, 43, 44, 45 ]
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [ 27, 26, 45, 46, 28 ]
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [ 13, 28, 46, 47, 48 ]
    geom.face(25).n_poi = 5; allocate(geom.face(25).poi(5)); geom.face(25).poi(1:5) = [ 29, 13, 48, 49, 50 ]
    geom.face(26).n_poi = 5; allocate(geom.face(26).poi(5)); geom.face(26).poi(1:5) = [ 30, 29, 50, 51, 31 ]
    geom.face(27).n_poi = 5; allocate(geom.face(27).poi(5)); geom.face(27).poi(1:5) = [ 32, 31, 52, 34, 33 ]
    geom.face(28).n_poi = 5; allocate(geom.face(28).poi(5)); geom.face(28).poi(1:5) = [ 35, 34, 52, 53, 54 ]
    geom.face(29).n_poi = 5; allocate(geom.face(29).poi(5)); geom.face(29).poi(1:5) = [ 36, 35, 55, 56, 37 ]
    geom.face(30).n_poi = 5; allocate(geom.face(30).poi(5)); geom.face(30).poi(1:5) = [ 38, 37, 56, 57, 58 ]
    geom.face(31).n_poi = 5; allocate(geom.face(31).poi(5)); geom.face(31).poi(1:5) = [ 21, 38, 58, 59, 39 ]
    geom.face(32).n_poi = 5; allocate(geom.face(32).poi(5)); geom.face(32).poi(1:5) = [ 40, 39, 59, 60, 41 ]
    geom.face(33).n_poi = 5; allocate(geom.face(33).poi(5)); geom.face(33).poi(1:5) = [ 42, 41, 61, 44, 43 ]
    geom.face(34).n_poi = 5; allocate(geom.face(34).poi(5)); geom.face(34).poi(1:5) = [ 45, 44, 61, 62, 63 ]
    geom.face(35).n_poi = 5; allocate(geom.face(35).poi(5)); geom.face(35).poi(1:5) = [ 46, 45, 64, 65, 47 ]
    geom.face(36).n_poi = 5; allocate(geom.face(36).poi(5)); geom.face(36).poi(1:5) = [ 48, 47, 65, 66, 49 ]
    geom.face(37).n_poi = 5; allocate(geom.face(37).poi(5)); geom.face(37).poi(1:5) = [ 50, 49, 66, 67, 51 ]
    geom.face(38).n_poi = 5; allocate(geom.face(38).poi(5)); geom.face(38).poi(1:5) = [ 31, 51, 67, 68, 69 ]
    geom.face(39).n_poi = 5; allocate(geom.face(39).poi(5)); geom.face(39).poi(1:5) = [ 52, 31, 69, 70, 53 ]
    geom.face(40).n_poi = 5; allocate(geom.face(40).poi(5)); geom.face(40).poi(1:5) = [ 54, 53, 70, 71, 72 ]
    geom.face(41).n_poi = 5; allocate(geom.face(41).poi(5)); geom.face(41).poi(1:5) = [ 35, 54, 72, 73, 55 ]
    geom.face(42).n_poi = 5; allocate(geom.face(42).poi(5)); geom.face(42).poi(1:5) = [ 56, 55, 73, 74, 57 ]
    geom.face(43).n_poi = 5; allocate(geom.face(43).poi(5)); geom.face(43).poi(1:5) = [ 58, 57, 75, 60, 59 ]
    geom.face(44).n_poi = 5; allocate(geom.face(44).poi(5)); geom.face(44).poi(1:5) = [ 41, 60, 75, 76, 77 ]
    geom.face(45).n_poi = 5; allocate(geom.face(45).poi(5)); geom.face(45).poi(1:5) = [ 61, 41, 77, 78, 62 ]
    geom.face(46).n_poi = 5; allocate(geom.face(46).poi(5)); geom.face(46).poi(1:5) = [ 63, 62, 78, 79, 80 ]
    geom.face(47).n_poi = 5; allocate(geom.face(47).poi(5)); geom.face(47).poi(1:5) = [ 45, 63, 80, 81, 64 ]
    geom.face(48).n_poi = 5; allocate(geom.face(48).poi(5)); geom.face(48).poi(1:5) = [ 65, 64, 81, 82, 66 ]
    geom.face(49).n_poi = 5; allocate(geom.face(49).poi(5)); geom.face(49).poi(1:5) = [ 67, 66, 83, 84, 68 ]
    geom.face(50).n_poi = 5; allocate(geom.face(50).poi(5)); geom.face(50).poi(1:5) = [ 69, 68, 84, 71, 70 ]
    geom.face(51).n_poi = 5; allocate(geom.face(51).poi(5)); geom.face(51).poi(1:5) = [ 72, 71, 85, 74, 73 ]
    geom.face(52).n_poi = 5; allocate(geom.face(52).poi(5)); geom.face(52).poi(1:5) = [ 57, 74, 85, 86, 87 ]
    geom.face(53).n_poi = 5; allocate(geom.face(53).poi(5)); geom.face(53).poi(1:5) = [ 75, 57, 87, 88, 76 ]
    geom.face(54).n_poi = 5; allocate(geom.face(54).poi(5)); geom.face(54).poi(1:5) = [ 77, 76, 88, 79, 78 ]
    geom.face(55).n_poi = 5; allocate(geom.face(55).poi(5)); geom.face(55).poi(1:5) = [ 80, 79, 89, 82, 81 ]
    geom.face(56).n_poi = 5; allocate(geom.face(56).poi(5)); geom.face(56).poi(1:5) = [ 66, 82, 89, 90, 83 ]
    geom.face(57).n_poi = 5; allocate(geom.face(57).poi(5)); geom.face(57).poi(1:5) = [ 84, 83, 90, 91, 71 ]
    geom.face(58).n_poi = 5; allocate(geom.face(58).poi(5)); geom.face(58).poi(1:5) = [ 85, 71, 91, 92, 86 ]
    geom.face(59).n_poi = 5; allocate(geom.face(59).poi(5)); geom.face(59).poi(1:5) = [ 87, 86, 92, 79, 88 ]
    geom.face(60).n_poi = 5; allocate(geom.face(60).poi(5)); geom.face(60).poi(1:5) = [ 89, 79, 92, 91, 90 ]
end subroutine Exam_Catalan_Pentagonal_Hexecontahedron

! ---------------------------------------------------------------------------------------

end module Exam_Catalan