!
! ---------------------------------------------------------------------------------------
!
!                                 Module for Exam_Archi
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/08/03
!                                                            Last  modified  : 2016/07/14
!
! ---------------------------------------------------------------------------------------
!
module Exam_Archi

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Archi_Cubeoctahedron                 !  6. V=12,  E=24,  F=14  (Rhombitetratetrahedron)
    public Exam_Archi_Icosidodecahedron              !  7. V=30,  E=60,  F=32
    public Exam_Archi_Rhombicuboctahedron            !  8. V=24,  E=48,  F=26
    public Exam_Archi_Snub_Cube                      !  9. V=24,  E=60,  F=38  (Snub Cuboctahedron)
    public Exam_Archi_Truncated_Cube                 ! 10. V=24,  E=36,  F=14  (Truncated Hexahedron)
    public Exam_Archi_Truncated_Cuboctahedron        ! 11. V=48,  E=72,  F=26  (Truncated Cuboctahedron)
    public Exam_Archi_Truncated_Dodecahedron         ! 12. V=60,  E=90,  F=32
    public Exam_Archi_Truncated_Icosahedron          ! 13. V=60,  E=90,  F=32  (Bucky Ball)
    public Exam_Archi_Truncated_Octahedron           ! 14. V=24,  E=36,  F=14  (Truncated Tetratetrahedron)
    public Exam_Archi_Truncated_Tetrahedron          ! 15. V=12,  E=18,  F=8

    public Exam_Archi_Truncated_Icosidodecahedron    ! 53. V=120, E=240, F=122 (Truncated Icosidodecahedron)
    public Exam_Archi_Rhombicosidodecahedron         ! 54. V=60,  E=120, F=62
    public Exam_Archi_Snub_Dodecahedron              ! 55. V=60,  E=150, F=92  (Snub Icosidodecahedron)

contains

! ---------------------------------------------------------------------------------------

! Example of Cubeoctahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Cubeoctahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "06_Cubeocta"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Cubeoctahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! allocate point, line and face structure
    geom.n_iniP = 12
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.00000d0,  0.00000d0, 10.00000d0 ]; geom.iniP( 2).pos(1:3) = [  8.66026d0,  0.00000d0,   5.00000d0 ]
    geom.iniP( 3).pos(1:3) = [  2.88675d0,  8.16497d0,  5.00000d0 ]; geom.iniP( 4).pos(1:3) = [ -8.66026d0,  0.00000d0,   5.00000d0 ]
    geom.iniP( 5).pos(1:3) = [ -2.88675d0, -8.16497d0,  5.00000d0 ]; geom.iniP( 6).pos(1:3) = [  8.66026d0,  0.00000d0,  -5.00000d0 ]
    geom.iniP( 7).pos(1:3) = [  5.77351d0, -8.16497d0,  0.00000d0 ]; geom.iniP( 8).pos(1:3) = [ -5.77351d0,  8.16497d0,   0.00000d0 ]
    geom.iniP( 9).pos(1:3) = [  2.88675d0,  8.16497d0, -5.00000d0 ]; geom.iniP(10).pos(1:3) = [ -8.66026d0,  0.00000d0,  -5.00000d0 ]
    geom.iniP(11).pos(1:3) = [ -2.88675d0, -8.16497d0, -5.00000d0 ]; geom.iniP(12).pos(1:3) = [  0.00000d0,  0.00000d0, -10.00000d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  4,  5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  7,  6 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  3,  9,  8 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  4,  8, 10 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  5, 11,  7 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  6, 12,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10, 12, 11 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  1,  3,  8,  4 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  1,  5,  7,  2 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  2,  6,  9,  3 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  4, 10, 11,  5 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  6,  7, 11, 12 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  8,  9, 12, 10 ]
end subroutine Exam_Archi_Cubeoctahedron

! ---------------------------------------------------------------------------------------

! Example of Icosidodecahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Icosidodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "07_Icosidodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Icosidodecahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.97d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! allocate point, line and face structure
    geom.n_iniP = 30
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  16.18033d0 ]; geom.iniP( 2).pos(1:3) = [   9.51058d0,   0.00000d0,  13.09020d0 ]
    geom.iniP( 3).pos(1:3) = [   4.25326d0,   8.50652d0,  13.09020d0 ]; geom.iniP( 4).pos(1:3) = [  -9.51058d0,   0.00000d0,  13.09020d0 ]
    geom.iniP( 5).pos(1:3) = [  -4.25326d0,  -8.50652d0,  13.09020d0 ]; geom.iniP( 6).pos(1:3) = [  15.38845d0,   0.00000d0,   5.00001d0 ]
    geom.iniP( 7).pos(1:3) = [  11.13519d0,  -8.50652d0,   8.09018d0 ]; geom.iniP( 8).pos(1:3) = [  -2.62865d0,  13.76384d0,   8.09018d0 ]
    geom.iniP( 9).pos(1:3) = [   6.88193d0,  13.76384d0,   5.00001d0 ]; geom.iniP(10).pos(1:3) = [ -15.38845d0,   0.00000d0,   5.00001d0 ]
    geom.iniP(11).pos(1:3) = [ -11.13519d0,   8.50652d0,   8.09018d0 ]; geom.iniP(12).pos(1:3) = [   2.62865d0, -13.76384d0,   8.09018d0 ]
    geom.iniP(13).pos(1:3) = [  -6.88193d0, -13.76384d0,   5.00001d0 ]; geom.iniP(14).pos(1:3) = [  15.38845d0,   0.00000d0,  -5.00001d0 ]
    geom.iniP(15).pos(1:3) = [  13.76384d0,   8.50652d0,   0.00000d0 ]; geom.iniP(16).pos(1:3) = [   8.50652d0, -13.76384d0,   0.00000d0 ]
    geom.iniP(17).pos(1:3) = [  -8.50652d0,  13.76384d0,   0.00000d0 ]; geom.iniP(18).pos(1:3) = [   6.88193d0,  13.76384d0,  -5.00001d0 ]
    geom.iniP(19).pos(1:3) = [ -15.38845d0,   0.00000d0,  -5.00001d0 ]; geom.iniP(20).pos(1:3) = [ -13.76384d0,  -8.50652d0,   0.00000d0 ]
    geom.iniP(21).pos(1:3) = [  -6.88193d0, -13.76384d0,  -5.00001d0 ]; geom.iniP(22).pos(1:3) = [   9.51058d0,   0.00000d0, -13.09020d0 ]
    geom.iniP(23).pos(1:3) = [  11.13519d0,  -8.50652d0,  -8.09018d0 ]; geom.iniP(24).pos(1:3) = [   2.62865d0, -13.76384d0,  -8.09018d0 ]
    geom.iniP(25).pos(1:3) = [ -11.13519d0,   8.50652d0,  -8.09018d0 ]; geom.iniP(26).pos(1:3) = [  -2.62865d0,  13.76384d0,  -8.09018d0 ]
    geom.iniP(27).pos(1:3) = [   4.25326d0,   8.50652d0, -13.09020d0 ]; geom.iniP(28).pos(1:3) = [  -9.51058d0,   0.00000d0, -13.09020d0 ]
    geom.iniP(29).pos(1:3) = [  -4.25326d0,  -8.50652d0, -13.09020d0 ]; geom.iniP(30).pos(1:3) = [   0.00000d0,   0.00000d0, -16.18033d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  4,  5 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  7,  6 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  3,  9,  8 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  4, 11, 10 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  5, 13, 12 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  6, 14, 15 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  7, 12, 16 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  8, 17, 11 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  9, 15, 18 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 10, 19, 20 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 13, 20, 21 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 14, 23, 22 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 16, 24, 23 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 17, 26, 25 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 18, 27, 26 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 19, 25, 28 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 21, 29, 24 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 22, 30, 27 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 28, 30, 29 ]
    geom.face(21).n_poi = 5; allocate(geom.face(21).poi(5)); geom.face(21).poi(1:5) = [  1,  3,  8, 11,  4 ]
    geom.face(22).n_poi = 5; allocate(geom.face(22).poi(5)); geom.face(22).poi(1:5) = [  1,  5, 12,  7,  2 ]
    geom.face(23).n_poi = 5; allocate(geom.face(23).poi(5)); geom.face(23).poi(1:5) = [  2,  6, 15,  9,  3 ]
    geom.face(24).n_poi = 5; allocate(geom.face(24).poi(5)); geom.face(24).poi(1:5) = [  4, 10, 20, 13,  5 ]
    geom.face(25).n_poi = 5; allocate(geom.face(25).poi(5)); geom.face(25).poi(1:5) = [  6,  7, 16, 23, 14 ]
    geom.face(26).n_poi = 5; allocate(geom.face(26).poi(5)); geom.face(26).poi(1:5) = [  8,  9, 18, 26, 17 ]
    geom.face(27).n_poi = 5; allocate(geom.face(27).poi(5)); geom.face(27).poi(1:5) = [ 10, 11, 17, 25, 19 ]
    geom.face(28).n_poi = 5; allocate(geom.face(28).poi(5)); geom.face(28).poi(1:5) = [ 12, 13, 21, 24, 16 ]
    geom.face(29).n_poi = 5; allocate(geom.face(29).poi(5)); geom.face(29).poi(1:5) = [ 14, 22, 27, 18, 15 ]
    geom.face(30).n_poi = 5; allocate(geom.face(30).poi(5)); geom.face(30).poi(1:5) = [ 19, 28, 29, 21, 20 ]
    geom.face(31).n_poi = 5; allocate(geom.face(31).poi(5)); geom.face(31).poi(1:5) = [ 22, 23, 24, 29, 30 ]
    geom.face(32).n_poi = 5; allocate(geom.face(32).poi(5)); geom.face(32).poi(1:5) = [ 25, 26, 27, 30, 28 ]
end subroutine Exam_Archi_Icosidodecahedron

! ---------------------------------------------------------------------------------------

! Example of Rhombicuboctahedron (Rhombicuboctahedron or Small Rhombicuboctahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Rhombicuboctahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "08_Rhombicubocta"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Rhombicuboctahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 24
    geom.n_face = 26

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  13.98965d0 ]; geom.iniP( 2).pos(1:3) = [   9.33950d0,   0.00000d0,  10.41561d0 ]
    geom.iniP( 3).pos(1:3) = [  -1.36774d0,   9.23881d0,  10.41561d0 ]; geom.iniP( 4).pos(1:3) = [  -8.93890d0,   2.70599d0,  10.41561d0 ]
    geom.iniP( 5).pos(1:3) = [  -1.36774d0,  -9.23881d0,  10.41561d0 ]; geom.iniP( 6).pos(1:3) = [   7.97177d0,   9.23881d0,   6.84154d0 ]
    geom.iniP( 7).pos(1:3) = [  13.60866d0,   2.70599d0,   1.78704d0 ]; geom.iniP( 8).pos(1:3) = [   7.97177d0,  -9.23881d0,   6.84154d0 ]
    geom.iniP( 9).pos(1:3) = [  -4.66976d0,  13.06565d0,   1.78704d0 ]; geom.iniP(10).pos(1:3) = [ -10.30664d0,  -6.53283d0,   6.84154d0 ]
    geom.iniP(11).pos(1:3) = [ -12.24092d0,   6.53283d0,   1.78704d0 ]; geom.iniP(12).pos(1:3) = [  -4.66976d0, -13.06565d0,   1.78704d0 ]
    geom.iniP(13).pos(1:3) = [   4.66976d0,  13.06565d0,  -1.78704d0 ]; geom.iniP(14).pos(1:3) = [  12.24092d0,  -6.53283d0,  -1.78704d0 ]
    geom.iniP(15).pos(1:3) = [  10.30664d0,   6.53283d0,  -6.84154d0 ]; geom.iniP(16).pos(1:3) = [   4.66976d0, -13.06565d0,  -1.78704d0 ]
    geom.iniP(17).pos(1:3) = [  -7.97177d0,   9.23881d0,  -6.84154d0 ]; geom.iniP(18).pos(1:3) = [ -13.60866d0,  -2.70599d0,  -1.78704d0 ]
    geom.iniP(19).pos(1:3) = [  -7.97177d0,  -9.23881d0,  -6.84154d0 ]; geom.iniP(20).pos(1:3) = [   1.36774d0,   9.23881d0, -10.41561d0 ]
    geom.iniP(21).pos(1:3) = [   8.93890d0,  -2.70599d0, -10.41561d0 ]; geom.iniP(22).pos(1:3) = [   1.36774d0,  -9.23881d0, -10.41561d0 ]
    geom.iniP(23).pos(1:3) = [  -9.33950d0,   0.00000d0, -10.41561d0 ]; geom.iniP(24).pos(1:3) = [   0.00000d0,   0.00000d0, -13.98965d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  3,  4 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  2,  7,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  5, 10, 12 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8, 16, 14 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  9, 17, 11 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 13, 15, 20 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 18, 23, 19 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 21, 22, 24 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [  1,  2,  6,  3 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  1,  4, 10,  5 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [  1,  5,  8,  2 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [  2,  8, 14,  7 ]
    geom.face(13).n_poi = 4; allocate(geom.face(13).poi(4)); geom.face(13).poi(1:4) = [  3,  6, 13,  9 ]
    geom.face(14).n_poi = 4; allocate(geom.face(14).poi(4)); geom.face(14).poi(1:4) = [  3,  9, 11,  4 ]
    geom.face(15).n_poi = 4; allocate(geom.face(15).poi(4)); geom.face(15).poi(1:4) = [  4, 11, 18, 10 ]
    geom.face(16).n_poi = 4; allocate(geom.face(16).poi(4)); geom.face(16).poi(1:4) = [  5, 12, 16,  8 ]
    geom.face(17).n_poi = 4; allocate(geom.face(17).poi(4)); geom.face(17).poi(1:4) = [  6,  7, 15, 13 ]
    geom.face(18).n_poi = 4; allocate(geom.face(18).poi(4)); geom.face(18).poi(1:4) = [  7, 14, 21, 15 ]
    geom.face(19).n_poi = 4; allocate(geom.face(19).poi(4)); geom.face(19).poi(1:4) = [  9, 13, 20, 17 ]
    geom.face(20).n_poi = 4; allocate(geom.face(20).poi(4)); geom.face(20).poi(1:4) = [ 10, 18, 19, 12 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [ 11, 17, 23, 18 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [ 12, 19, 22, 16 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [ 14, 16, 22, 21 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [ 15, 21, 24, 20 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [ 17, 20, 24, 23 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [ 19, 23, 24, 22 ]
end subroutine Exam_Archi_Rhombicuboctahedron

! ---------------------------------------------------------------------------------------

! Example of Snub Cube (Snub Cuboctahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Snub_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "09_Snub_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Snub Cube"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 24
    geom.n_face = 38

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  13.43713d0 ]; geom.iniP( 2).pos(1:3) = [   9.28194d0,   0.00000d0,   9.71614d0 ]
    geom.iniP( 3).pos(1:3) = [   3.89511d0,   8.42512d0,   9.71614d0 ]; geom.iniP( 4).pos(1:3) = [  -6.01283d0,   7.07109d0,   9.71614d0 ]
    geom.iniP( 5).pos(1:3) = [  -8.94159d0,  -2.49045d0,   9.71614d0 ]; geom.iniP( 6).pos(1:3) = [  -1.49173d0,  -9.16128d0,   9.71614d0 ]
    geom.iniP( 7).pos(1:3) = [   7.79021d0,  -9.16128d0,   5.99509d0 ]; geom.iniP( 8).pos(1:3) = [  13.17708d0,  -2.49045d0,   0.84898d0 ]
    geom.iniP( 9).pos(1:3) = [  11.05933d0,   7.07109d0,   2.87207d0 ]; geom.iniP(10).pos(1:3) = [   3.26911d0,  13.00571d0,   0.84898d0 ]
    geom.iniP(11).pos(1:3) = [  -6.63883d0,  11.65173d0,   0.84898d0 ]; geom.iniP(12).pos(1:3) = [ -12.55110d0,   3.84448d0,   2.87207d0 ]
    geom.iniP(13).pos(1:3) = [  -8.75656d0,  -9.77915d0,   2.87207d0 ]; geom.iniP(14).pos(1:3) = [   0.34035d0, -13.40595d0,   0.84898d0 ]
    geom.iniP(15).pos(1:3) = [   8.31560d0,  -9.77915d0,  -3.97199d0 ]; geom.iniP(16).pos(1:3) = [  10.24828d0,  -1.13642d0,  -8.61620d0 ]
    geom.iniP(17).pos(1:3) = [   8.13056d0,   8.42512d0,  -6.59310d0 ]; geom.iniP(18).pos(1:3) = [  -1.49173d0,  10.91557d0,  -7.69304d0 ]
    geom.iniP(19).pos(1:3) = [ -10.09298d0,   5.93467d0,  -6.59310d0 ]; geom.iniP(20).pos(1:3) = [ -12.36601d0,  -3.44423d0,  -3.97199d0 ]
    geom.iniP(21).pos(1:3) = [  -5.57187d0, -10.29771d0,  -6.59310d0 ]; geom.iniP(22).pos(1:3) = [   2.40337d0,  -6.67085d0, -11.41408d0 ]
    geom.iniP(23).pos(1:3) = [   2.30277d0,   3.22661d0, -12.83920d0 ]; geom.iniP(24).pos(1:3) = [  -6.29848d0,  -1.75428d0, -11.73922d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  1,  5,  6 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  2,  7,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  2,  8,  9 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  2,  9,  3 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  3,  9, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4, 11, 12 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4, 12,  5 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  5, 13,  6 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  6, 13, 14 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  6, 14,  7 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  7, 14, 15 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  7, 15,  8 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  8, 15, 16 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  9, 17, 10 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 10, 17, 18 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 10, 18, 11 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 11, 18, 19 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 11, 19, 12 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 12, 19, 20 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 13, 20, 21 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 13, 21, 14 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 15, 22, 16 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 16, 22, 23 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 16, 23, 17 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 17, 23, 18 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 19, 24, 20 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 20, 24, 21 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 21, 24, 22 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 22, 24, 23 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [  1,  6,  7,  2 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [  3, 10, 11,  4 ]
    geom.face(35).n_poi = 4; allocate(geom.face(35).poi(4)); geom.face(35).poi(1:4) = [  5, 12, 20, 13 ]
    geom.face(36).n_poi = 4; allocate(geom.face(36).poi(4)); geom.face(36).poi(1:4) = [  8, 16, 17,  9 ]
    geom.face(37).n_poi = 4; allocate(geom.face(37).poi(4)); geom.face(37).poi(1:4) = [ 14, 21, 22, 15 ]
    geom.face(38).n_poi = 4; allocate(geom.face(38).poi(4)); geom.face(38).poi(1:4) = [ 18, 23, 24, 19 ]
end subroutine Exam_Archi_Snub_Cube

! ---------------------------------------------------------------------------------------

! Example of Truncated Cube (Truncated Hexahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "10_Trunc_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated cube"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "opt"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "off"    ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "on"     ! [on, off], Unpaired scaffold nucleotides
    end if

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.95d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XZ"

    ! allocate point, line and face structure
    geom.n_iniP = 24
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  17.78826d0 ]; geom.iniP( 2).pos(1:3) = [   9.59685d0,   0.00000d0,  14.97742d0 ]
    geom.iniP( 3).pos(1:3) = [  -8.19143d0,   5.00001d0,  14.97742d0 ]; geom.iniP( 4).pos(1:3) = [   4.38678d0,  -8.53555d0,  14.97742d0 ]
    geom.iniP( 5).pos(1:3) = [  14.97742d0,   5.00001d0,   8.19143d0 ]; geom.iniP( 6).pos(1:3) = [ -15.38905d0,   3.53555d0,   8.19143d0 ]
    geom.iniP( 7).pos(1:3) = [ -10.17900d0,  12.07110d0,   8.19143d0 ]; geom.iniP( 8).pos(1:3) = [   2.39921d0, -15.60664d0,   8.19143d0 ]
    geom.iniP( 9).pos(1:3) = [  17.37668d0,   3.53555d0,  -1.40543d0 ]; geom.iniP(10).pos(1:3) = [  12.98985d0,  12.07110d0,   1.40543d0 ]
    geom.iniP(11).pos(1:3) = [ -17.37668d0,  -3.53555d0,   1.40543d0 ]; geom.iniP(12).pos(1:3) = [  -4.79843d0,  17.07110d0,   1.40543d0 ]
    geom.iniP(13).pos(1:3) = [  -4.79843d0, -17.07110d0,   1.40543d0 ]; geom.iniP(14).pos(1:3) = [   4.79843d0, -17.07110d0,  -1.40543d0 ]
    geom.iniP(15).pos(1:3) = [  15.38905d0,  -3.53555d0,  -8.19143d0 ]; geom.iniP(16).pos(1:3) = [   4.79843d0,  17.07110d0,  -1.40543d0 ]
    geom.iniP(17).pos(1:3) = [ -12.98985d0, -12.07110d0,  -1.40543d0 ]; geom.iniP(18).pos(1:3) = [ -14.97742d0,  -5.00001d0,  -8.19143d0 ]
    geom.iniP(19).pos(1:3) = [  -2.39921d0,  15.60664d0,  -8.19143d0 ]; geom.iniP(20).pos(1:3) = [  10.17900d0, -12.07110d0,  -8.19143d0 ]
    geom.iniP(21).pos(1:3) = [   8.19143d0,  -5.00001d0, -14.97742d0 ]; geom.iniP(22).pos(1:3) = [  -9.59685d0,   0.00000d0, -14.97742d0 ]
    geom.iniP(23).pos(1:3) = [  -4.38678d0,   8.53555d0, -14.97742d0 ]; geom.iniP(24).pos(1:3) = [   0.00000d0,   0.00000d0, -17.78826d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 3; allocate(geom.face(1).poi(3));  geom.face( 1).poi(1:3) = [  1,  4,  2 ]
    geom.face(2).n_poi  = 3; allocate(geom.face(2).poi(3));  geom.face( 2).poi(1:3) = [  3,  7,  6 ]
    geom.face(3).n_poi  = 3; allocate(geom.face(3).poi(3));  geom.face( 3).poi(1:3) = [  5,  9, 10 ]
    geom.face(4).n_poi  = 3; allocate(geom.face(4).poi(3));  geom.face( 4).poi(1:3) = [  8, 13, 14 ]
    geom.face(5).n_poi  = 3; allocate(geom.face(5).poi(3));  geom.face( 5).poi(1:3) = [ 11, 18, 17 ]
    geom.face(6).n_poi  = 3; allocate(geom.face(6).poi(3));  geom.face( 6).poi(1:3) = [ 12, 16, 19 ]
    geom.face(7).n_poi  = 3; allocate(geom.face(7).poi(3));  geom.face( 7).poi(1:3) = [ 15, 20, 21 ]
    geom.face(8).n_poi  = 3; allocate(geom.face(8).poi(3));  geom.face( 8).poi(1:3) = [ 22, 23, 24 ]
    geom.face(9).n_poi  = 8; allocate(geom.face(9).poi(8));  geom.face( 9).poi(1:8) = [  1,  2,  5, 10, 16, 12,  7,  3 ]
    geom.face(10).n_poi = 8; allocate(geom.face(10).poi(8)); geom.face(10).poi(1:8) = [  1,  3,  6, 11, 17, 13,  8,  4 ]
    geom.face(11).n_poi = 8; allocate(geom.face(11).poi(8)); geom.face(11).poi(1:8) = [  2,  4,  8, 14, 20, 15,  9,  5 ]
    geom.face(12).n_poi = 8; allocate(geom.face(12).poi(8)); geom.face(12).poi(1:8) = [  6,  7, 12, 19, 23, 22, 18, 11 ]
    geom.face(13).n_poi = 8; allocate(geom.face(13).poi(8)); geom.face(13).poi(1:8) = [  9, 15, 21, 24, 23, 19, 16, 10 ]
    geom.face(14).n_poi = 8; allocate(geom.face(14).poi(8)); geom.face(14).poi(1:8) = [ 13, 17, 18, 22, 24, 21, 20, 14 ]
end subroutine Exam_Archi_Truncated_Cube

! ---------------------------------------------------------------------------------------

! Example of Truncated Cuboctahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Cuboctahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "11_Trunc_Cubocta"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated cuboctahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.97d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XYZ"

    ! allocate point, line and face structure
    geom.n_iniP = 48
    geom.n_face = 26

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  23.17622d0 ]; geom.iniP( 2).pos(1:3) = [   9.76453d0,   0.00000d0,  21.01875d0 ]
    geom.iniP( 3).pos(1:3) = [  -0.47666d0,   9.75288d0,  21.01875d0 ]; geom.iniP( 4).pos(1:3) = [  -7.71827d0,  -5.98116d0,  21.01875d0 ]
    geom.iniP( 5).pos(1:3) = [   9.28787d0,   9.75288d0,  18.86135d0 ]; geom.iniP( 6).pos(1:3) = [  15.85538d0,  -5.98116d0,  15.81033d0 ]
    geom.iniP( 7).pos(1:3) = [  -8.67160d0,  13.52461d0,  16.70395d0 ]; geom.iniP( 8).pos(1:3) = [  -8.86902d0, -14.43979d0,  15.81033d0 ]
    geom.iniP( 9).pos(1:3) = [ -15.91320d0,  -2.20943d0,  16.70395d0 ]; geom.iniP(10).pos(1:3) = [  14.90206d0,  13.52461d0,  11.49553d0 ]
    geom.iniP(11).pos(1:3) = [  14.70463d0, -14.43979d0,  10.60192d0 ]; geom.iniP(12).pos(1:3) = [  21.46957d0,  -2.20943d0,   8.44452d0 ]
    geom.iniP(13).pos(1:3) = [ -10.49647d0,  18.85864d0,   8.44452d0 ]; geom.iniP(14).pos(1:3) = [ -16.38987d0,   7.54345d0,  14.54656d0 ]
    geom.iniP(15).pos(1:3) = [ -17.06398d0, -10.66805d0,  11.49553d0 ]; geom.iniP(16).pos(1:3) = [  -2.77817d0, -20.42095d0,  10.60192d0 ]
    geom.iniP(17).pos(1:3) = [  13.07719d0,  18.85864d0,   3.23610d0 ]; geom.iniP(18).pos(1:3) = [  20.99291d0,   7.54345d0,   6.28713d0 ]
    geom.iniP(19).pos(1:3) = [  20.31882d0, -10.66805d0,   3.23610d0 ]; geom.iniP(20).pos(1:3) = [   6.98634d0, -20.42095d0,   8.44452d0 ]
    geom.iniP(21).pos(1:3) = [ -18.21474d0,  12.87750d0,   6.28713d0 ]; geom.iniP(22).pos(1:3) = [  -4.88225d0,  22.63037d0,   1.07870d0 ]
    geom.iniP(23).pos(1:3) = [ -19.16804d0, -12.87750d0,   1.97232d0 ]; geom.iniP(24).pos(1:3) = [  -4.88225d0, -22.63037d0,   1.07870d0 ]
    geom.iniP(25).pos(1:3) = [  19.16804d0,  12.87750d0,  -1.97232d0 ]; geom.iniP(26).pos(1:3) = [   4.88225d0,  22.63037d0,  -1.07870d0 ]
    geom.iniP(27).pos(1:3) = [  18.21474d0, -12.87750d0,  -6.28713d0 ]; geom.iniP(28).pos(1:3) = [   4.88225d0, -22.63037d0,  -1.07870d0 ]
    geom.iniP(29).pos(1:3) = [ -20.31882d0,  10.66805d0,  -3.23610d0 ]; geom.iniP(30).pos(1:3) = [  -6.98634d0,  20.42095d0,  -8.44452d0 ]
    geom.iniP(31).pos(1:3) = [ -20.99291d0,  -7.54345d0,  -6.28713d0 ]; geom.iniP(32).pos(1:3) = [ -13.07719d0, -18.85864d0,  -3.23610d0 ]
    geom.iniP(33).pos(1:3) = [  17.06398d0,  10.66805d0, -11.49553d0 ]; geom.iniP(34).pos(1:3) = [   2.77817d0,  20.42095d0, -10.60192d0 ]
    geom.iniP(35).pos(1:3) = [  16.38987d0,  -7.54345d0, -14.54656d0 ]; geom.iniP(36).pos(1:3) = [  10.49647d0, -18.85864d0,  -8.44452d0 ]
    geom.iniP(37).pos(1:3) = [ -21.46957d0,   2.20943d0,  -8.44452d0 ]; geom.iniP(38).pos(1:3) = [ -14.70463d0,  14.43979d0, -10.60192d0 ]
    geom.iniP(39).pos(1:3) = [ -14.90206d0, -13.52461d0, -11.49553d0 ]; geom.iniP(40).pos(1:3) = [  15.91320d0,   2.20943d0, -16.70395d0 ]
    geom.iniP(41).pos(1:3) = [   8.86902d0,  14.43979d0, -15.81033d0 ]; geom.iniP(42).pos(1:3) = [   8.67160d0, -13.52461d0, -16.70395d0 ]
    geom.iniP(43).pos(1:3) = [ -15.85538d0,   5.98116d0, -15.81033d0 ]; geom.iniP(44).pos(1:3) = [  -9.28787d0,  -9.75288d0, -18.86135d0 ]
    geom.iniP(45).pos(1:3) = [   7.71827d0,   5.98116d0, -21.01875d0 ]; geom.iniP(46).pos(1:3) = [   0.47666d0,  -9.75288d0, -21.01875d0 ]
    geom.iniP(47).pos(1:3) = [  -9.76453d0,   0.00000d0, -21.01875d0 ]; geom.iniP(48).pos(1:3) = [   0.00000d0,   0.00000d0, -23.17622d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 4; allocate(geom.face(1).poi(4));  geom.face( 1).poi(1:4) = [  1,  2,  5,  3 ]
    geom.face(2).n_poi  = 4; allocate(geom.face(2).poi(4));  geom.face( 2).poi(1:4) = [  4,  9, 15,  8 ]
    geom.face(3).n_poi  = 4; allocate(geom.face(3).poi(4));  geom.face( 3).poi(1:4) = [  6, 11, 19, 12 ]
    geom.face(4).n_poi  = 4; allocate(geom.face(4).poi(4));  geom.face( 4).poi(1:4) = [  7, 13, 21, 14 ]
    geom.face(5).n_poi  = 4; allocate(geom.face(5).poi(4));  geom.face( 5).poi(1:4) = [ 10, 18, 25, 17 ]
    geom.face(6).n_poi  = 4; allocate(geom.face(6).poi(4));  geom.face( 6).poi(1:4) = [ 16, 24, 28, 20 ]
    geom.face(7).n_poi  = 4; allocate(geom.face(7).poi(4));  geom.face( 7).poi(1:4) = [ 22, 26, 34, 30 ]
    geom.face(8).n_poi  = 4; allocate(geom.face(8).poi(4));  geom.face( 8).poi(1:4) = [ 23, 31, 39, 32 ]
    geom.face(9).n_poi  = 4; allocate(geom.face(9).poi(4));  geom.face( 9).poi(1:4) = [ 27, 36, 42, 35 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 29, 38, 43, 37 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 33, 40, 45, 41 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 44, 47, 48, 46 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [  1,  3,  7, 14,  9,  4 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [  2,  6, 12, 18, 10,  5 ]
    geom.face(15).n_poi = 6; allocate(geom.face(15).poi(6)); geom.face(15).poi(1:6) = [  8, 15, 23, 32, 24, 16 ]
    geom.face(16).n_poi = 6; allocate(geom.face(16).poi(6)); geom.face(16).poi(1:6) = [ 11, 20, 28, 36, 27, 19 ]
    geom.face(17).n_poi = 6; allocate(geom.face(17).poi(6)); geom.face(17).poi(1:6) = [ 13, 22, 30, 38, 29, 21 ]
    geom.face(18).n_poi = 6; allocate(geom.face(18).poi(6)); geom.face(18).poi(1:6) = [ 17, 25, 33, 41, 34, 26 ]
    geom.face(19).n_poi = 6; allocate(geom.face(19).poi(6)); geom.face(19).poi(1:6) = [ 31, 37, 43, 47, 44, 39 ]
    geom.face(20).n_poi = 6; allocate(geom.face(20).poi(6)); geom.face(20).poi(1:6) = [ 35, 42, 46, 48, 45, 40 ]
    geom.face(21).n_poi = 8; allocate(geom.face(21).poi(8)); geom.face(21).poi(1:8) = [  1,  4,  8, 16, 20, 11,  6,  2 ]
    geom.face(22).n_poi = 8; allocate(geom.face(22).poi(8)); geom.face(22).poi(1:8) = [  3,  5, 10, 17, 26, 22, 13,  7 ]
    geom.face(23).n_poi = 8; allocate(geom.face(23).poi(8)); geom.face(23).poi(1:8) = [  9, 14, 21, 29, 37, 31, 23, 15 ]
    geom.face(24).n_poi = 8; allocate(geom.face(24).poi(8)); geom.face(24).poi(1:8) = [ 12, 19, 27, 35, 40, 33, 25, 18 ]
    geom.face(25).n_poi = 8; allocate(geom.face(25).poi(8)); geom.face(25).poi(1:8) = [ 24, 32, 39, 44, 46, 42, 36, 28 ]
    geom.face(26).n_poi = 8; allocate(geom.face(26).poi(8)); geom.face(26).poi(1:8) = [ 30, 34, 41, 45, 48, 47, 43, 38 ]
end subroutine Exam_Archi_truncated_Cuboctahedron

! ---------------------------------------------------------------------------------------

! Example of Truncated Dodecahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "12_Trunc_Dodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated dodecahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.9d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 60
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  29.69470d0 ]; geom.iniP( 2).pos(1:3) = [   9.85725d0,   0.00000d0,  28.01074d0 ]
    geom.iniP( 3).pos(1:3) = [  -8.49499d0,   5.00000d0,  28.01074d0 ]; geom.iniP( 4).pos(1:3) = [   4.78480d0,  -8.61805d0,  28.01074d0 ]
    geom.iniP( 5).pos(1:3) = [  17.31158d0,   5.00000d0,  23.60244d0 ]; geom.iniP( 6).pos(1:3) = [ -17.45539d0,   4.47213d0,  23.60244d0 ]
    geom.iniP( 7).pos(1:3) = [ -12.38297d0,  13.09021d0,  23.60244d0 ]; geom.iniP( 8).pos(1:3) = [   4.03179d0, -17.56234d0,  23.60244d0 ]
    geom.iniP( 9).pos(1:3) = [  24.30055d0,   4.47213d0,  16.46967d0 ]; geom.iniP(10).pos(1:3) = [  19.51574d0,  13.09021d0,  18.15349d0 ]
    geom.iniP(11).pos(1:3) = [ -23.45864d0,  -1.38197d0,  18.15349d0 ]; geom.iniP(12).pos(1:3) = [ -10.17881d0,  21.18039d0,  18.15349d0 ]
    geom.iniP(13).pos(1:3) = [  -1.97145d0, -23.41646d0,  18.15349d0 ]; geom.iniP(14).pos(1:3) = [   7.88578d0, -23.41646d0,  16.46967d0 ]
    geom.iniP(15).pos(1:3) = [  28.15454d0,  -1.38197d0,   9.33690d0 ]; geom.iniP(16).pos(1:3) = [  15.62776d0,  21.18039d0,  13.74520d0 ]
    geom.iniP(17).pos(1:3) = [ -24.21165d0, -10.32625d0,  13.74520d0 ]; geom.iniP(18).pos(1:3) = [ -28.09963d0,  -2.23607d0,   9.33690d0 ]
    geom.iniP(19).pos(1:3) = [ -11.68487d0,  25.65252d0,   9.33690d0 ]; geom.iniP(20).pos(1:3) = [  -2.72447d0,  26.18039d0,  13.74520d0 ]
    geom.iniP(21).pos(1:3) = [ -10.93186d0, -23.94432d0,  13.74520d0 ]; geom.iniP(22).pos(1:3) = [  14.87475d0, -23.94432d0,   9.33690d0 ]
    geom.iniP(23).pos(1:3) = [  27.40153d0, -10.32625d0,   4.92861d0 ]; geom.iniP(24).pos(1:3) = [  29.60572d0,  -2.23607d0,  -0.52033d0 ]
    geom.iniP(25).pos(1:3) = [  14.12170d0,  25.65252d0,   4.92861d0 ]; geom.iniP(26).pos(1:3) = [   7.13277d0,  26.18039d0,  12.06138d0 ]
    geom.iniP(27).pos(1:3) = [ -19.42685d0, -18.94432d0,  12.06138d0 ]; geom.iniP(28).pos(1:3) = [ -29.60572d0,   2.23607d0,   0.52033d0 ]
    geom.iniP(29).pos(1:3) = [ -16.32587d0,  24.79844d0,   0.52033d0 ]; geom.iniP(30).pos(1:3) = [ -15.57285d0, -24.79844d0,   4.92861d0 ]
    geom.iniP(31).pos(1:3) = [  16.32587d0, -24.79844d0,  -0.52033d0 ]; geom.iniP(32).pos(1:3) = [  22.32908d0, -18.94432d0,   4.92861d0 ]
    geom.iniP(33).pos(1:3) = [  28.09963d0,   2.23607d0,  -9.33690d0 ]; geom.iniP(34).pos(1:3) = [  15.57285d0,  24.79844d0,  -4.92861d0 ]
    geom.iniP(35).pos(1:3) = [ -28.15454d0,   1.38197d0,  -9.33690d0 ]; geom.iniP(36).pos(1:3) = [ -27.40153d0,  10.32625d0,  -4.92861d0 ]
    geom.iniP(37).pos(1:3) = [ -22.32908d0,  18.94432d0,  -4.92861d0 ]; geom.iniP(38).pos(1:3) = [ -14.87475d0,  23.94432d0,  -9.33690d0 ]
    geom.iniP(39).pos(1:3) = [ -14.12170d0, -25.65252d0,  -4.92861d0 ]; geom.iniP(40).pos(1:3) = [  11.68487d0, -25.65252d0,  -9.33690d0 ]
    geom.iniP(41).pos(1:3) = [  23.45864d0,   1.38197d0, -18.15349d0 ]; geom.iniP(42).pos(1:3) = [  24.21165d0,  10.32625d0, -13.74520d0 ]
    geom.iniP(43).pos(1:3) = [  19.42685d0,  18.94432d0, -12.06138d0 ]; geom.iniP(44).pos(1:3) = [  10.93186d0,  23.94432d0, -13.74520d0 ]
    geom.iniP(45).pos(1:3) = [ -24.30055d0,  -4.47213d0, -16.46967d0 ]; geom.iniP(46).pos(1:3) = [  -7.88578d0,  23.41646d0, -16.46967d0 ]
    geom.iniP(47).pos(1:3) = [ -15.62776d0, -21.18039d0, -13.74520d0 ]; geom.iniP(48).pos(1:3) = [  -7.13277d0, -26.18039d0, -12.06138d0 ]
    geom.iniP(49).pos(1:3) = [   2.72447d0, -26.18039d0, -13.74520d0 ]; geom.iniP(50).pos(1:3) = [  10.17881d0, -21.18039d0, -18.15349d0 ]
    geom.iniP(51).pos(1:3) = [  17.45539d0,  -4.47213d0, -23.60244d0 ]; geom.iniP(52).pos(1:3) = [   1.97145d0,  23.41646d0, -18.15349d0 ]
    geom.iniP(53).pos(1:3) = [ -19.51574d0, -13.09021d0, -18.15349d0 ]; geom.iniP(54).pos(1:3) = [ -17.31158d0,  -5.00000d0, -23.60244d0 ]
    geom.iniP(55).pos(1:3) = [  -4.03179d0,  17.56234d0, -23.60244d0 ]; geom.iniP(56).pos(1:3) = [  12.38297d0, -13.09021d0, -23.60244d0 ]
    geom.iniP(57).pos(1:3) = [   8.49499d0,  -5.00000d0, -28.01074d0 ]; geom.iniP(58).pos(1:3) = [  -9.85725d0,   0.00000d0, -28.01074d0 ]
    geom.iniP(59).pos(1:3) = [  -4.78480d0,   8.61805d0, -28.01074d0 ]; geom.iniP(60).pos(1:3) = [   0.00000d0,   0.00000d0, -29.69470d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi =  3; allocate(geom.face( 1).poi(3));  geom.face( 1).poi(1:3)  = [  1,  4,  2 ]
    geom.face( 2).n_poi =  3; allocate(geom.face( 2).poi(3));  geom.face( 2).poi(1:3)  = [  3,  7,  6 ]
    geom.face( 3).n_poi =  3; allocate(geom.face( 3).poi(3));  geom.face( 3).poi(1:3)  = [  5,  9, 10 ]
    geom.face( 4).n_poi =  3; allocate(geom.face( 4).poi(3));  geom.face( 4).poi(1:3)  = [  8, 13, 14 ]
    geom.face( 5).n_poi =  3; allocate(geom.face( 5).poi(3));  geom.face( 5).poi(1:3)  = [ 11, 18, 17 ]
    geom.face( 6).n_poi =  3; allocate(geom.face( 6).poi(3));  geom.face( 6).poi(1:3)  = [ 12, 20, 19 ]
    geom.face( 7).n_poi =  3; allocate(geom.face( 7).poi(3));  geom.face( 7).poi(1:3)  = [ 15, 23, 24 ]
    geom.face( 8).n_poi =  3; allocate(geom.face( 8).poi(3));  geom.face( 8).poi(1:3)  = [ 16, 25, 26 ]
    geom.face( 9).n_poi =  3; allocate(geom.face( 9).poi(3));  geom.face( 9).poi(1:3)  = [ 21, 27, 30 ]
    geom.face(10).n_poi =  3; allocate(geom.face(10).poi(3));  geom.face(10).poi(1:3)  = [ 22, 31, 32 ]
    geom.face(11).n_poi =  3; allocate(geom.face(11).poi(3));  geom.face(11).poi(1:3)  = [ 28, 36, 35 ]
    geom.face(12).n_poi =  3; allocate(geom.face(12).poi(3));  geom.face(12).poi(1:3)  = [ 29, 38, 37 ]
    geom.face(13).n_poi =  3; allocate(geom.face(13).poi(3));  geom.face(13).poi(1:3)  = [ 33, 41, 42 ]
    geom.face(14).n_poi =  3; allocate(geom.face(14).poi(3));  geom.face(14).poi(1:3)  = [ 34, 43, 44 ]
    geom.face(15).n_poi =  3; allocate(geom.face(15).poi(3));  geom.face(15).poi(1:3)  = [ 39, 47, 48 ]
    geom.face(16).n_poi =  3; allocate(geom.face(16).poi(3));  geom.face(16).poi(1:3)  = [ 40, 49, 50 ]
    geom.face(17).n_poi =  3; allocate(geom.face(17).poi(3));  geom.face(17).poi(1:3)  = [ 45, 54, 53 ]
    geom.face(18).n_poi =  3; allocate(geom.face(18).poi(3));  geom.face(18).poi(1:3)  = [ 46, 52, 55 ]
    geom.face(19).n_poi =  3; allocate(geom.face(19).poi(3));  geom.face(19).poi(1:3)  = [ 51, 56, 57 ]
    geom.face(20).n_poi =  3; allocate(geom.face(20).poi(3));  geom.face(20).poi(1:3)  = [ 58, 59, 60 ]
    geom.face(21).n_poi = 10; allocate(geom.face(21).poi(10)); geom.face(21).poi(1:10) = [  1,  2,  5, 10, 16, 26, 20, 12,  7,  3 ]
    geom.face(22).n_poi = 10; allocate(geom.face(22).poi(10)); geom.face(22).poi(1:10) = [  1,  3,  6, 11, 17, 27, 21, 13,  8,  4 ]
    geom.face(23).n_poi = 10; allocate(geom.face(23).poi(10)); geom.face(23).poi(1:10) = [  2,  4,  8, 14, 22, 32, 23, 15,  9,  5 ]
    geom.face(24).n_poi = 10; allocate(geom.face(24).poi(10)); geom.face(24).poi(1:10) = [  6,  7, 12, 19, 29, 37, 36, 28, 18, 11 ]
    geom.face(25).n_poi = 10; allocate(geom.face(25).poi(10)); geom.face(25).poi(1:10) = [  9, 15, 24, 33, 42, 43, 34, 25, 16, 10 ]
    geom.face(26).n_poi = 10; allocate(geom.face(26).poi(10)); geom.face(26).poi(1:10) = [ 13, 21, 30, 39, 48, 49, 40, 31, 22, 14 ]
    geom.face(27).n_poi = 10; allocate(geom.face(27).poi(10)); geom.face(27).poi(1:10) = [ 17, 18, 28, 35, 45, 53, 47, 39, 30, 27 ]
    geom.face(28).n_poi = 10; allocate(geom.face(28).poi(10)); geom.face(28).poi(1:10) = [ 19, 20, 26, 25, 34, 44, 52, 46, 38, 29 ]
    geom.face(29).n_poi = 10; allocate(geom.face(29).poi(10)); geom.face(29).poi(1:10) = [ 23, 32, 31, 40, 50, 56, 51, 41, 33, 24 ]
    geom.face(30).n_poi = 10; allocate(geom.face(30).poi(10)); geom.face(30).poi(1:10) = [ 35, 36, 37, 38, 46, 55, 59, 58, 54, 45 ]
    geom.face(31).n_poi = 10; allocate(geom.face(31).poi(10)); geom.face(31).poi(1:10) = [ 41, 51, 57, 60, 59, 55, 52, 44, 43, 42 ]
    geom.face(32).n_poi = 10; allocate(geom.face(32).poi(10)); geom.face(32).poi(1:10) = [ 47, 53, 54, 58, 60, 57, 56, 50, 49, 48 ]
end subroutine Exam_Archi_Truncated_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Truncated Icosahedron (bucky ball)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "13_Trunc_Icosa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated icosahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.97d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 60
    geom.n_face = 32

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  24.78026d0 ]; geom.iniP( 2).pos(1:3) = [   9.79434d0,   0.00000d0,  22.76250d0 ]
    geom.iniP( 3).pos(1:3) = [  -5.52068d0,   8.09019d0,  22.76250d0 ]; geom.iniP( 4).pos(1:3) = [  -3.57076d0,  -9.12025d0,  22.76250d0 ]
    geom.iniP( 5).pos(1:3) = [  14.06800d0,   8.09019d0,  18.72701d0 ]; geom.iniP( 6).pos(1:3) = [  12.27684d0,  -9.12025d0,  19.49772d0 ]
    geom.iniP( 7).pos(1:3) = [ -14.61213d0,   7.06014d0,  18.72701d0 ]; geom.iniP( 8).pos(1:3) = [  -1.24704d0,  16.18039d0,  18.72701d0 ]
    geom.iniP( 9).pos(1:3) = [   4.01676d0, -14.75688d0,  19.49772d0 ]; geom.iniP(10).pos(1:3) = [ -12.66220d0, -10.15032d0,  18.72701d0 ]
    geom.iniP(11).pos(1:3) = [  20.82416d0,   7.06014d0,  11.42673d0 ]; geom.iniP(12).pos(1:3) = [   8.54732d0,  16.18039d0,  16.70927d0 ]
    geom.iniP(13).pos(1:3) = [  19.03299d0, -10.15032d0,  12.19745d0 ]; geom.iniP(14).pos(1:3) = [ -15.95730d0,  14.51371d0,  12.19745d0 ]
    geom.iniP(15).pos(1:3) = [ -18.18289d0,  -2.06012d0,  16.70927d0 ]; geom.iniP(16).pos(1:3) = [  -7.69722d0,  20.15033d0,  12.19745d0 ]
    geom.iniP(17).pos(1:3) = [   2.51281d0, -21.42357d0,  12.19745d0 ]; geom.iniP(18).pos(1:3) = [ -14.16613d0, -16.81699d0,  11.42673d0 ]
    geom.iniP(19).pos(1:3) = [  19.47898d0,  14.51371d0,   4.89717d0 ]; geom.iniP(20).pos(1:3) = [  23.30665d0,  -2.06012d0,   8.16195d0 ]
    geom.iniP(21).pos(1:3) = [  11.89149d0,  20.15033d0,   8.16195d0 ]; geom.iniP(22).pos(1:3) = [  17.52904d0, -16.81699d0,   4.89717d0 ]
    geom.iniP(23).pos(1:3) = [ -20.87321d0,  12.84705d0,   3.65015d0 ]; geom.iniP(24).pos(1:3) = [ -23.09880d0,  -3.72680d0,   8.16195d0 ]
    geom.iniP(25).pos(1:3) = [  -4.35305d0,  24.12030d0,   3.65015d0 ]; geom.iniP(26).pos(1:3) = [   9.26896d0, -22.45364d0,   4.89717d0 ]
    geom.iniP(27).pos(1:3) = [  -6.57864d0, -22.45364d0,   8.16195d0 ]; geom.iniP(28).pos(1:3) = [ -20.61630d0, -12.84705d0,   4.89717d0 ]
    geom.iniP(29).pos(1:3) = [  20.61630d0,  12.84705d0,  -4.89717d0 ]; geom.iniP(30).pos(1:3) = [  24.44387d0,  -3.72680d0,  -1.63239d0 ]
    geom.iniP(31).pos(1:3) = [   5.44129d0,  24.12030d0,   1.63239d0 ]; geom.iniP(32).pos(1:3) = [  20.87321d0, -12.84705d0,  -3.65015d0 ]
    geom.iniP(33).pos(1:3) = [ -17.52904d0,  16.81699d0,  -4.89717d0 ]; geom.iniP(34).pos(1:3) = [ -24.44387d0,   3.72680d0,   1.63239d0 ]
    geom.iniP(35).pos(1:3) = [  -9.26896d0,  22.45364d0,  -4.89717d0 ]; geom.iniP(36).pos(1:3) = [   4.35305d0, -24.12030d0,  -3.65015d0 ]
    geom.iniP(37).pos(1:3) = [  -5.44129d0, -24.12030d0,  -1.63239d0 ]; geom.iniP(38).pos(1:3) = [ -19.47898d0, -14.51371d0,  -4.89717d0 ]
    geom.iniP(39).pos(1:3) = [  14.16613d0,  16.81699d0, -11.42673d0 ]; geom.iniP(40).pos(1:3) = [  23.09880d0,   3.72680d0,  -8.16195d0 ]
    geom.iniP(41).pos(1:3) = [   6.57864d0,  22.45364d0,  -8.16195d0 ]; geom.iniP(42).pos(1:3) = [  15.95730d0, -14.51371d0, -12.19745d0 ]
    geom.iniP(43).pos(1:3) = [ -19.03299d0,  10.15032d0, -12.19745d0 ]; geom.iniP(44).pos(1:3) = [ -23.30665d0,   2.06012d0,  -8.16195d0 ]
    geom.iniP(45).pos(1:3) = [  -2.51281d0,  21.42357d0, -12.19745d0 ]; geom.iniP(46).pos(1:3) = [   7.69722d0, -20.15033d0, -12.19745d0 ]
    geom.iniP(47).pos(1:3) = [ -11.89149d0, -20.15033d0,  -8.16195d0 ]; geom.iniP(48).pos(1:3) = [ -20.82416d0,  -7.06014d0, -11.42673d0 ]
    geom.iniP(49).pos(1:3) = [  12.66220d0,  10.15032d0, -18.72701d0 ]; geom.iniP(50).pos(1:3) = [  18.18289d0,   2.06012d0, -16.70927d0 ]
    geom.iniP(51).pos(1:3) = [  14.61213d0,  -7.06014d0, -18.72701d0 ]; geom.iniP(52).pos(1:3) = [ -12.27684d0,   9.12025d0, -19.49772d0 ]
    geom.iniP(53).pos(1:3) = [  -4.01676d0,  14.75688d0, -19.49772d0 ]; geom.iniP(54).pos(1:3) = [   1.24704d0, -16.18039d0, -18.72701d0 ]
    geom.iniP(55).pos(1:3) = [  -8.54732d0, -16.18039d0, -16.70927d0 ]; geom.iniP(56).pos(1:3) = [ -14.06800d0,  -8.09019d0, -18.72701d0 ]
    geom.iniP(57).pos(1:3) = [   3.57076d0,   9.12025d0, -22.76250d0 ]; geom.iniP(58).pos(1:3) = [   5.52068d0,  -8.09019d0, -22.76250d0 ]
    geom.iniP(59).pos(1:3) = [  -9.79434d0,   0.00000d0, -22.76250d0 ]; geom.iniP(60).pos(1:3) = [   0.00000d0,   0.00000d0, -24.78026d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 5; allocate(geom.face( 1).poi(5)); geom.face( 1).poi(1:5) = [  1,  4,  9,  6,  2 ]
    geom.face( 2).n_poi = 5; allocate(geom.face( 2).poi(5)); geom.face( 2).poi(1:5) = [  3,  8, 16, 14,  7 ]
    geom.face( 3).n_poi = 5; allocate(geom.face( 3).poi(5)); geom.face( 3).poi(1:5) = [  5, 11, 19, 21, 12 ]
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [ 10, 15, 24, 28, 18 ]
    geom.face( 5).n_poi = 5; allocate(geom.face( 5).poi(5)); geom.face( 5).poi(1:5) = [ 13, 22, 32, 30, 20 ]
    geom.face( 6).n_poi = 5; allocate(geom.face( 6).poi(5)); geom.face( 6).poi(1:5) = [ 17, 27, 37, 36, 26 ]
    geom.face( 7).n_poi = 5; allocate(geom.face( 7).poi(5)); geom.face( 7).poi(1:5) = [ 23, 33, 43, 44, 34 ]
    geom.face( 8).n_poi = 5; allocate(geom.face( 8).poi(5)); geom.face( 8).poi(1:5) = [ 25, 31, 41, 45, 35 ]
    geom.face( 9).n_poi = 5; allocate(geom.face( 9).poi(5)); geom.face( 9).poi(1:5) = [ 29, 40, 50, 49, 39 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 38, 48, 56, 55, 47 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 42, 46, 54, 58, 51 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 52, 53, 57, 60, 59 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [  1,  2,  5, 12,  8,  3 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [  1,  3,  7, 15, 10,  4 ]
    geom.face(15).n_poi = 6; allocate(geom.face(15).poi(6)); geom.face(15).poi(1:6) = [  2,  6, 13, 20, 11,  5 ]
    geom.face(16).n_poi = 6; allocate(geom.face(16).poi(6)); geom.face(16).poi(1:6) = [  4, 10, 18, 27, 17,  9 ]
    geom.face(17).n_poi = 6; allocate(geom.face(17).poi(6)); geom.face(17).poi(1:6) = [  6,  9, 17, 26, 22, 13 ]
    geom.face(18).n_poi = 6; allocate(geom.face(18).poi(6)); geom.face(18).poi(1:6) = [  7, 14, 23, 34, 24, 15 ]
    geom.face(19).n_poi = 6; allocate(geom.face(19).poi(6)); geom.face(19).poi(1:6) = [  8, 12, 21, 31, 25, 16 ]
    geom.face(20).n_poi = 6; allocate(geom.face(20).poi(6)); geom.face(20).poi(1:6) = [ 11, 20, 30, 40, 29, 19 ]
    geom.face(21).n_poi = 6; allocate(geom.face(21).poi(6)); geom.face(21).poi(1:6) = [ 14, 16, 25, 35, 33, 23 ]
    geom.face(22).n_poi = 6; allocate(geom.face(22).poi(6)); geom.face(22).poi(1:6) = [ 18, 28, 38, 47, 37, 27 ]
    geom.face(23).n_poi = 6; allocate(geom.face(23).poi(6)); geom.face(23).poi(1:6) = [ 19, 29, 39, 41, 31, 21 ]
    geom.face(24).n_poi = 6; allocate(geom.face(24).poi(6)); geom.face(24).poi(1:6) = [ 22, 26, 36, 46, 42, 32 ]
    geom.face(25).n_poi = 6; allocate(geom.face(25).poi(6)); geom.face(25).poi(1:6) = [ 24, 34, 44, 48, 38, 28 ]
    geom.face(26).n_poi = 6; allocate(geom.face(26).poi(6)); geom.face(26).poi(1:6) = [ 30, 32, 42, 51, 50, 40 ]
    geom.face(27).n_poi = 6; allocate(geom.face(27).poi(6)); geom.face(27).poi(1:6) = [ 33, 35, 45, 53, 52, 43 ]
    geom.face(28).n_poi = 6; allocate(geom.face(28).poi(6)); geom.face(28).poi(1:6) = [ 36, 37, 47, 55, 54, 46 ]
    geom.face(29).n_poi = 6; allocate(geom.face(29).poi(6)); geom.face(29).poi(1:6) = [ 39, 49, 57, 53, 45, 41 ]
    geom.face(30).n_poi = 6; allocate(geom.face(30).poi(6)); geom.face(30).poi(1:6) = [ 43, 52, 59, 56, 48, 44 ]
    geom.face(31).n_poi = 6; allocate(geom.face(31).poi(6)); geom.face(31).poi(1:6) = [ 49, 50, 51, 58, 60, 57 ]
    geom.face(32).n_poi = 6; allocate(geom.face(32).poi(6)); geom.face(32).poi(1:6) = [ 54, 55, 56, 59, 60, 58 ]
end subroutine Exam_Archi_Truncated_Icosahedron

! ---------------------------------------------------------------------------------------

! Example of Truncated Octahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "14_Trunc_Octa"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated octahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   =-0.5d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 24
    geom.n_face = 14

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  15.81138d0 ]; geom.iniP( 2).pos(1:3) = [   9.48686d0,   0.00000d0,  12.64913d0 ]
    geom.iniP( 3).pos(1:3) = [  -6.32457d0,   7.07107d0,  12.64913d0 ]; geom.iniP( 4).pos(1:3) = [  -1.05409d0,  -9.42810d0,  12.64913d0 ]
    geom.iniP( 5).pos(1:3) = [  12.64913d0,   7.07107d0,   6.32457d0 ]; geom.iniP( 6).pos(1:3) = [   8.43276d0,  -9.42810d0,   9.48686d0 ]
    geom.iniP( 7).pos(1:3) = [ -13.70323d0,   4.71406d0,   6.32457d0 ]; geom.iniP( 8).pos(1:3) = [  -3.16228d0,  14.14216d0,   6.32457d0 ]
    geom.iniP( 9).pos(1:3) = [  -8.43276d0, -11.78513d0,   6.32457d0 ]; geom.iniP(10).pos(1:3) = [  14.75733d0,   4.71406d0,  -3.16228d0 ]
    geom.iniP(11).pos(1:3) = [   6.32457d0,  14.14216d0,   3.16228d0 ]; geom.iniP(12).pos(1:3) = [  10.54094d0, -11.78513d0,   0.00000d0 ]
    geom.iniP(13).pos(1:3) = [ -10.54094d0,  11.78513d0,   0.00000d0 ]; geom.iniP(14).pos(1:3) = [ -14.75733d0,  -4.71406d0,   3.16228d0 ]
    geom.iniP(15).pos(1:3) = [  -6.32457d0, -14.14216d0,  -3.16228d0 ]; geom.iniP(16).pos(1:3) = [   8.43276d0,  11.78513d0,  -6.32457d0 ]
    geom.iniP(17).pos(1:3) = [  13.70323d0,  -4.71406d0,  -6.32457d0 ]; geom.iniP(18).pos(1:3) = [   3.16228d0, -14.14216d0,  -6.32457d0 ]
    geom.iniP(19).pos(1:3) = [  -8.43276d0,   9.42810d0,  -9.48686d0 ]; geom.iniP(20).pos(1:3) = [ -12.64913d0,  -7.07107d0,  -6.32457d0 ]
    geom.iniP(21).pos(1:3) = [   1.05409d0,   9.42810d0, -12.64913d0 ]; geom.iniP(22).pos(1:3) = [   6.32457d0,  -7.07107d0, -12.64913d0 ]
    geom.iniP(23).pos(1:3) = [  -9.48686d0,   0.00000d0, -12.64913d0 ]; geom.iniP(24).pos(1:3) = [   0.00000d0,   0.00000d0, -15.81138d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  1,  4,  6,  2 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3,  8, 13,  7 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  5, 10, 16, 11 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  9, 14, 20, 15 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 12, 18, 22, 17 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [ 19, 21, 24, 23 ]
    geom.face( 7).n_poi = 6; allocate(geom.face( 7).poi(6)); geom.face( 7).poi(1:6) = [  1,  2,  5, 11,  8,  3 ]
    geom.face( 8).n_poi = 6; allocate(geom.face( 8).poi(6)); geom.face( 8).poi(1:6) = [  1,  3,  7, 14,  9,  4 ]
    geom.face( 9).n_poi = 6; allocate(geom.face( 9).poi(6)); geom.face( 9).poi(1:6) = [  2,  6, 12, 17, 10,  5 ]
    geom.face(10).n_poi = 6; allocate(geom.face(10).poi(6)); geom.face(10).poi(1:6) = [  4,  9, 15, 18, 12,  6 ]
    geom.face(11).n_poi = 6; allocate(geom.face(11).poi(6)); geom.face(11).poi(1:6) = [  7, 13, 19, 23, 20, 14 ]
    geom.face(12).n_poi = 6; allocate(geom.face(12).poi(6)); geom.face(12).poi(1:6) = [  8, 11, 16, 21, 19, 13 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [ 10, 17, 22, 24, 21, 16 ]
    geom.face(14).n_poi = 6; allocate(geom.face(14).poi(6)); geom.face(14).poi(1:6) = [ 15, 20, 23, 24, 22, 18 ]
end subroutine Exam_Archi_Truncated_Octahedron

! ---------------------------------------------------------------------------------------

! Example of Truncated Tetrahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    prob.name_file = "15_Trunc_Tetra"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated tetrahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.5d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 12
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.00000d0,  0.00000d0, 11.72605d0 ]; geom.iniP( 2).pos(1:3) = [   9.04536d0,   0.00000d0,   7.46204d0 ]
    geom.iniP( 3).pos(1:3) = [ -7.53780d0,  5.00000d0,  7.46204d0 ]; geom.iniP( 4).pos(1:3) = [   3.51764d0,  -8.33335d0,   7.46204d0 ]
    geom.iniP( 5).pos(1:3) = [ 10.55292d0,  5.00000d0, -1.06600d0 ]; geom.iniP( 6).pos(1:3) = [ -11.55794d0,   1.66667d0,  -1.06600d0 ]
    geom.iniP( 7).pos(1:3) = [ -6.03024d0, 10.00002d0, -1.06600d0 ]; geom.iniP( 8).pos(1:3) = [  -0.50252d0, -11.66665d0,  -1.06600d0 ]
    geom.iniP( 9).pos(1:3) = [  6.53275d0,  1.66667d0, -9.59405d0 ]; geom.iniP(10).pos(1:3) = [   3.01512d0,  10.00002d0,  -5.33003d0 ]
    geom.iniP(11).pos(1:3) = [ -8.04032d0, -6.66668d0, -5.33003d0 ]; geom.iniP(12).pos(1:3) = [   1.00504d0,  -6.66668d0,  -9.59405d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1,  4,  2 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 3,  7,  6 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 5,  9, 10 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 8, 11, 12 ]
    geom.face(5).n_poi = 6; allocate(geom.face(5).poi(6)); geom.face(5).poi(1:6) = [ 1,  2,  5, 10,  7,  3 ]
    geom.face(6).n_poi = 6; allocate(geom.face(6).poi(6)); geom.face(6).poi(1:6) = [ 1,  3,  6, 11,  8,  4 ]
    geom.face(7).n_poi = 6; allocate(geom.face(7).poi(6)); geom.face(7).poi(1:6) = [ 2,  4,  8, 12,  9,  5 ]
    geom.face(8).n_poi = 6; allocate(geom.face(8).poi(6)); geom.face(8).poi(1:6) = [ 6,  7, 10,  9, 12, 11 ]
end subroutine Exam_Archi_Truncated_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of Truncated Icosidodecahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Truncated_Icosidodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "53_Trunc_Icosidodeca"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Truncated Icosidodecahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.93d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 120
    geom.n_face = 62

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(  1).pos(1:3) = [   0.00000d0,   0.00000d0,  38.02410d0 ]; geom.iniP(  2).pos(1:3) = [   9.91319d0,   0.00000d0,  36.70911d0 ]
    geom.iniP(  3).pos(1:3) = [  -0.17443d0,   9.91165d0,  36.70911d0 ]; geom.iniP(  4).pos(1:3) = [  -8.33548d0,  -5.36573d0,  36.70911d0 ]
    geom.iniP(  5).pos(1:3) = [   9.73878d0,   9.91165d0,  35.39416d0 ]; geom.iniP(  6).pos(1:3) = [  17.61760d0,  -5.36573d0,  33.26649d0 ]
    geom.iniP(  7).pos(1:3) = [  -8.68434d0,  14.45761d0,  34.07917d0 ]; geom.iniP(  8).pos(1:3) = [ -11.90939d0, -14.04769d0,  33.26649d0 ]
    geom.iniP(  9).pos(1:3) = [ -16.84540d0,  -0.81981d0,  34.07917d0 ]; geom.iniP( 10).pos(1:3) = [  17.26874d0,  14.45761d0,  30.63655d0 ]
    geom.iniP( 11).pos(1:3) = [  20.17038d0, -14.04769d0,  29.01115d0 ]; geom.iniP( 12).pos(1:3) = [  25.14760d0,  -0.81981d0,  28.50888d0 ]
    geom.iniP( 13).pos(1:3) = [ -12.54050d0,  21.81305d0,  28.50888d0 ]; geom.iniP( 14).pos(1:3) = [ -17.01985d0,   9.09184d0,  32.76422d0 ]
    geom.iniP( 15).pos(1:3) = [ -20.41931d0,  -9.50177d0,  30.63655d0 ]; geom.iniP( 16).pos(1:3) = [  -9.35661d0, -22.72965d0,  29.01115d0 ]
    geom.iniP( 17).pos(1:3) = [  19.53931d0,  21.81305d0,  24.25358d0 ]; geom.iniP( 18).pos(1:3) = [  24.97315d0,   9.09184d0,  27.19393d0 ]
    geom.iniP( 19).pos(1:3) = [  27.70038d0,  -9.50177d0,  24.25358d0 ]; geom.iniP( 20).pos(1:3) = [  16.59651d0, -22.72965d0,  25.56856d0 ]
    geom.iniP( 21).pos(1:3) = [ -20.87598d0,  16.44732d0,  27.19393d0 ]; geom.iniP( 22).pos(1:3) = [ -10.26993d0,  29.16852d0,  22.12591d0 ]
    geom.iniP( 23).pos(1:3) = [ -26.37646d0, -13.63777d0,  23.75130d0 ]; geom.iniP( 24).pos(1:3) = [  -1.65220d0, -28.09538d0,  25.56856d0 ]
    geom.iniP( 25).pos(1:3) = [ -15.31375d0, -26.86565d0,  22.12591d0 ]; geom.iniP( 26).pos(1:3) = [  27.24372d0,  16.44732d0,  20.81095d0 ]
    geom.iniP( 27).pos(1:3) = [  15.68319d0,  29.16852d0,  18.68328d0 ]; geom.iniP( 28).pos(1:3) = [  31.65643d0, -13.63777d0,  16.05338d0 ]
    geom.iniP( 29).pos(1:3) = [   8.26100d0, -28.09538d0,  24.25358d0 ]; geom.iniP( 30).pos(1:3) = [  20.55256d0, -26.86565d0,  17.36833d0 ]
    geom.iniP( 31).pos(1:3) = [ -26.94089d0,  18.43703d0,  19.49600d0 ]; geom.iniP( 32).pos(1:3) = [  -2.73994d0,  33.71444d0,  17.36833d0 ]
    geom.iniP( 33).pos(1:3) = [ -16.33484d0,  31.15823d0,  14.42798d0 ]; geom.iniP( 34).pos(1:3) = [ -32.44137d0, -11.64806d0,  16.05338d0 ]
    geom.iniP( 35).pos(1:3) = [ -23.82364d0, -22.31973d0,  19.49600d0 ]; geom.iniP( 36).pos(1:3) = [  -7.60934d0, -32.23141d0,  18.68328d0 ]
    geom.iniP( 37).pos(1:3) = [  31.09200d0,  18.43703d0,  11.79804d0 ]; geom.iniP( 38).pos(1:3) = [   7.17326d0,  33.71444d0,  16.05338d0 ]
    geom.iniP( 39).pos(1:3) = [  19.53143d0,  31.15823d0,   9.67037d0 ]; geom.iniP( 40).pos(1:3) = [  35.50472d0, -11.64806d0,   7.04047d0 ]
    geom.iniP( 41).pos(1:3) = [  28.08256d0, -22.31973d0,  12.61076d0 ]; geom.iniP( 42).pos(1:3) = [  12.21708d0, -32.23141d0,  16.05338d0 ]
    geom.iniP( 43).pos(1:3) = [ -32.89803d0,  14.30103d0,  12.61076d0 ]; geom.iniP( 44).pos(1:3) = [ -24.67036d0,  25.79250d0,  13.11303d0 ]
    geom.iniP( 45).pos(1:3) = [  -8.80488d0,  35.70419d0,   9.67037d0 ]; geom.iniP( 46).pos(1:3) = [ -29.88859d0, -20.33002d0,  11.79804d0 ]
    geom.iniP( 47).pos(1:3) = [ -36.29749d0,  -4.29258d0,  10.48309d0 ]; geom.iniP( 48).pos(1:3) = [  -3.65326d0, -36.36742d0,  10.48309d0 ]
    geom.iniP( 49).pos(1:3) = [  35.04805d0,  14.30103d0,   3.59783d0 ]; geom.iniP( 50).pos(1:3) = [  27.23584d0,  25.79250d0,   6.22775d0 ]
    geom.iniP( 51).pos(1:3) = [  11.02155d0,  35.70419d0,   7.04047d0 ]; geom.iniP( 52).pos(1:3) = [  31.93081d0, -20.33002d0,   3.59783d0 ]
    geom.iniP( 53).pos(1:3) = [  37.77532d0,  -4.29258d0,   0.65748d0 ]; geom.iniP( 54).pos(1:3) = [   6.25994d0, -36.36742d0,   9.16810d0 ]
    geom.iniP( 55).pos(1:3) = [ -30.62746d0,  21.65650d0,   6.22775d0 ]; geom.iniP( 56).pos(1:3) = [ -36.47194d0,   5.61907d0,   9.16810d0 ]
    geom.iniP( 57).pos(1:3) = [  -4.95660d0,  37.69390d0,   0.65748d0 ]; geom.iniP( 58).pos(1:3) = [ -31.19193d0, -21.65650d0,   1.97245d0 ]
    geom.iniP( 59).pos(1:3) = [ -37.60084d0,  -5.61907d0,   0.65748d0 ]; geom.iniP( 60).pos(1:3) = [  -4.95660d0, -37.69390d0,   0.65748d0 ]
    geom.iniP( 61).pos(1:3) = [  31.19193d0,  21.65650d0,  -1.97245d0 ]; geom.iniP( 62).pos(1:3) = [  37.60084d0,   5.61907d0,  -0.65748d0 ]
    geom.iniP( 63).pos(1:3) = [   4.95660d0,  37.69390d0,  -0.65748d0 ]; geom.iniP( 64).pos(1:3) = [  30.62746d0, -21.65650d0,  -6.22775d0 ]
    geom.iniP( 65).pos(1:3) = [  36.47194d0,  -5.61907d0,  -9.16810d0 ]; geom.iniP( 66).pos(1:3) = [   4.95660d0, -37.69390d0,  -0.65748d0 ]
    geom.iniP( 67).pos(1:3) = [ -31.93081d0,  20.33002d0,  -3.59783d0 ]; geom.iniP( 68).pos(1:3) = [ -37.77532d0,   4.29258d0,  -0.65748d0 ]
    geom.iniP( 69).pos(1:3) = [  -6.25994d0,  36.36742d0,  -9.16810d0 ]; geom.iniP( 70).pos(1:3) = [ -27.23584d0, -25.79250d0,  -6.22775d0 ]
    geom.iniP( 71).pos(1:3) = [ -35.04805d0, -14.30103d0,  -3.59783d0 ]; geom.iniP( 72).pos(1:3) = [ -11.02155d0, -35.70419d0,  -7.04047d0 ]
    geom.iniP( 73).pos(1:3) = [  29.88859d0,  20.33002d0, -11.79804d0 ]; geom.iniP( 74).pos(1:3) = [  36.29749d0,   4.29258d0, -10.48309d0 ]
    geom.iniP( 75).pos(1:3) = [   3.65326d0,  36.36742d0, -10.48309d0 ]; geom.iniP( 76).pos(1:3) = [  24.67036d0, -25.79250d0, -13.11303d0 ]
    geom.iniP( 77).pos(1:3) = [  32.89803d0, -14.30103d0, -12.61076d0 ]; geom.iniP( 78).pos(1:3) = [   8.80488d0, -35.70419d0,  -9.67037d0 ]
    geom.iniP( 79).pos(1:3) = [ -28.08256d0,  22.31973d0, -12.61076d0 ]; geom.iniP( 80).pos(1:3) = [ -35.50472d0,  11.64806d0,  -7.04047d0 ]
    geom.iniP( 81).pos(1:3) = [ -12.21708d0,  32.23141d0, -16.05338d0 ]; geom.iniP( 82).pos(1:3) = [ -31.09200d0, -18.43703d0, -11.79804d0 ]
    geom.iniP( 83).pos(1:3) = [ -19.53143d0, -31.15823d0,  -9.67037d0 ]; geom.iniP( 84).pos(1:3) = [  -7.17326d0, -33.71444d0, -16.05338d0 ]
    geom.iniP( 85).pos(1:3) = [  23.82364d0,  22.31973d0, -19.49600d0 ]; geom.iniP( 86).pos(1:3) = [  32.44137d0,  11.64806d0, -16.05338d0 ]
    geom.iniP( 87).pos(1:3) = [   7.60934d0,  32.23141d0, -18.68328d0 ]; geom.iniP( 88).pos(1:3) = [  26.94089d0, -18.43703d0, -19.49600d0 ]
    geom.iniP( 89).pos(1:3) = [  16.33484d0, -31.15823d0, -14.42798d0 ]; geom.iniP( 90).pos(1:3) = [   2.73994d0, -33.71444d0, -17.36833d0 ]
    geom.iniP( 91).pos(1:3) = [ -31.65643d0,  13.63777d0, -16.05338d0 ]; geom.iniP( 92).pos(1:3) = [ -20.55256d0,  26.86565d0, -17.36833d0 ]
    geom.iniP( 93).pos(1:3) = [  -8.26100d0,  28.09538d0, -24.25358d0 ]; geom.iniP( 94).pos(1:3) = [ -27.24372d0, -16.44732d0, -20.81095d0 ]
    geom.iniP( 95).pos(1:3) = [ -15.68319d0, -29.16852d0, -18.68328d0 ]; geom.iniP( 96).pos(1:3) = [  26.37646d0,  13.63777d0, -23.75130d0 ]
    geom.iniP( 97).pos(1:3) = [  15.31375d0,  26.86565d0, -22.12591d0 ]; geom.iniP( 98).pos(1:3) = [   1.65220d0,  28.09538d0, -25.56856d0 ]
    geom.iniP( 99).pos(1:3) = [  20.87598d0, -16.44732d0, -27.19393d0 ]; geom.iniP(100).pos(1:3) = [  10.26993d0, -29.16852d0, -22.12591d0 ]
    geom.iniP(101).pos(1:3) = [ -27.70038d0,   9.50177d0, -24.25358d0 ]; geom.iniP(102).pos(1:3) = [ -16.59651d0,  22.72965d0, -25.56856d0 ]
    geom.iniP(103).pos(1:3) = [ -24.97315d0,  -9.09184d0, -27.19393d0 ]; geom.iniP(104).pos(1:3) = [ -19.53931d0, -21.81305d0, -24.25358d0 ]
    geom.iniP(105).pos(1:3) = [  20.41931d0,   9.50177d0, -30.63655d0 ]; geom.iniP(106).pos(1:3) = [   9.35661d0,  22.72965d0, -29.01115d0 ]
    geom.iniP(107).pos(1:3) = [  17.01985d0,  -9.09184d0, -32.76422d0 ]; geom.iniP(108).pos(1:3) = [  12.54050d0, -21.81305d0, -28.50888d0 ]
    geom.iniP(109).pos(1:3) = [ -25.14760d0,   0.81981d0, -28.50888d0 ]; geom.iniP(110).pos(1:3) = [ -20.17038d0,  14.04769d0, -29.01115d0 ]
    geom.iniP(111).pos(1:3) = [ -17.26874d0, -14.45761d0, -30.63655d0 ]; geom.iniP(112).pos(1:3) = [  16.84540d0,   0.81981d0, -34.07917d0 ]
    geom.iniP(113).pos(1:3) = [  11.90939d0,  14.04769d0, -33.26649d0 ]; geom.iniP(114).pos(1:3) = [   8.68434d0, -14.45761d0, -34.07917d0 ]
    geom.iniP(115).pos(1:3) = [ -17.61760d0,   5.36573d0, -33.26649d0 ]; geom.iniP(116).pos(1:3) = [  -9.73878d0,  -9.91165d0, -35.39416d0 ]
    geom.iniP(117).pos(1:3) = [   8.33548d0,   5.36573d0, -36.70911d0 ]; geom.iniP(118).pos(1:3) = [   0.17443d0,  -9.91165d0, -36.70911d0 ]
    geom.iniP(119).pos(1:3) = [  -9.91319d0,   0.00000d0, -36.70911d0 ]; geom.iniP(120).pos(1:3) = [   0.00000d0,   0.00000d0, -38.02410d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi =  4; allocate(geom.face( 1).poi(4));  geom.face( 1).poi(1: 4) = [   1,  2,  5,  3 ]
    geom.face( 2).n_poi =  4; allocate(geom.face( 2).poi(4));  geom.face( 2).poi(1: 4) = [   4,  9, 15,  8 ]
    geom.face( 3).n_poi =  4; allocate(geom.face( 3).poi(4));  geom.face( 3).poi(1: 4) = [   6, 11, 19, 12 ]
    geom.face( 4).n_poi =  4; allocate(geom.face( 4).poi(4));  geom.face( 4).poi(1: 4) = [   7, 13, 21, 14 ]
    geom.face( 5).n_poi =  4; allocate(geom.face( 5).poi(4));  geom.face( 5).poi(1: 4) = [  10, 18, 26, 17 ]
    geom.face( 6).n_poi =  4; allocate(geom.face( 6).poi(4));  geom.face( 6).poi(1: 4) = [  16, 25, 36, 24 ]
    geom.face( 7).n_poi =  4; allocate(geom.face( 7).poi(4));  geom.face( 7).poi(1: 4) = [  20, 29, 42, 30 ]
    geom.face( 8).n_poi =  4; allocate(geom.face( 8).poi(4));  geom.face( 8).poi(1: 4) = [  22, 32, 45, 33 ]
    geom.face( 9).n_poi =  4; allocate(geom.face( 9).poi(4));  geom.face( 9).poi(1: 4) = [  23, 34, 46, 35 ]
    geom.face(10).n_poi =  4; allocate(geom.face(10).poi(4));  geom.face(10).poi(1: 4) = [  27, 39, 51, 38 ]
    geom.face(11).n_poi =  4; allocate(geom.face(11).poi(4));  geom.face(11).poi(1: 4) = [  28, 41, 52, 40 ]
    geom.face(12).n_poi =  4; allocate(geom.face(12).poi(4));  geom.face(12).poi(1: 4) = [  31, 44, 55, 43 ]
    geom.face(13).n_poi =  4; allocate(geom.face(13).poi(4));  geom.face(13).poi(1: 4) = [  37, 49, 61, 50 ]
    geom.face(14).n_poi =  4; allocate(geom.face(14).poi(4));  geom.face(14).poi(1: 4) = [  47, 56, 68, 59 ]
    geom.face(15).n_poi =  4; allocate(geom.face(15).poi(4));  geom.face(15).poi(1: 4) = [  48, 60, 66, 54 ]
    geom.face(16).n_poi =  4; allocate(geom.face(16).poi(4));  geom.face(16).poi(1: 4) = [  53, 65, 74, 62 ]
    geom.face(17).n_poi =  4; allocate(geom.face(17).poi(4));  geom.face(17).poi(1: 4) = [  57, 63, 75, 69 ]
    geom.face(18).n_poi =  4; allocate(geom.face(18).poi(4));  geom.face(18).poi(1: 4) = [  58, 71, 82, 70 ]
    geom.face(19).n_poi =  4; allocate(geom.face(19).poi(4));  geom.face(19).poi(1: 4) = [  64, 76, 88, 77 ]
    geom.face(20).n_poi =  4; allocate(geom.face(20).poi(4));  geom.face(20).poi(1: 4) = [  67, 79, 91, 80 ]
    geom.face(21).n_poi =  4; allocate(geom.face(21).poi(4));  geom.face(21).poi(1: 4) = [  72, 83, 95, 84 ]
    geom.face(22).n_poi =  4; allocate(geom.face(22).poi(4));  geom.face(22).poi(1: 4) = [  73, 86, 96, 85 ]
    geom.face(23).n_poi =  4; allocate(geom.face(23).poi(4));  geom.face(23).poi(1: 4) = [  78, 90,100, 89 ]
    geom.face(24).n_poi =  4; allocate(geom.face(24).poi(4));  geom.face(24).poi(1: 4) = [  81, 93,102, 92 ]
    geom.face(25).n_poi =  4; allocate(geom.face(25).poi(4));  geom.face(25).poi(1: 4) = [  87, 97,106, 98 ]
    geom.face(26).n_poi =  4; allocate(geom.face(26).poi(4));  geom.face(26).poi(1: 4) = [  94,103,111,104 ]
    geom.face(27).n_poi =  4; allocate(geom.face(27).poi(4));  geom.face(27).poi(1: 4) = [  99,108,114,107 ]
    geom.face(28).n_poi =  4; allocate(geom.face(28).poi(4));  geom.face(28).poi(1: 4) = [ 101,110,115,109 ]
    geom.face(29).n_poi =  4; allocate(geom.face(29).poi(4));  geom.face(29).poi(1: 4) = [ 105,112,117,113 ]
    geom.face(30).n_poi =  4; allocate(geom.face(30).poi(4));  geom.face(30).poi(1: 4) = [ 116,119,120,118 ]
    geom.face(31).n_poi =  6; allocate(geom.face(31).poi(6));  geom.face(31).poi(1: 6) = [   1,  3,  7, 14,  9,  4 ]
    geom.face(32).n_poi =  6; allocate(geom.face(32).poi(6));  geom.face(32).poi(1: 6) = [   2,  6, 12, 18, 10,  5 ]
    geom.face(33).n_poi =  6; allocate(geom.face(33).poi(6));  geom.face(33).poi(1: 6) = [   8, 15, 23, 35, 25, 16 ]
    geom.face(34).n_poi =  6; allocate(geom.face(34).poi(6));  geom.face(34).poi(1: 6) = [  11, 20, 30, 41, 28, 19 ]
    geom.face(35).n_poi =  6; allocate(geom.face(35).poi(6));  geom.face(35).poi(1: 6) = [  13, 22, 33, 44, 31, 21 ]
    geom.face(36).n_poi =  6; allocate(geom.face(36).poi(6));  geom.face(36).poi(1: 6) = [  17, 26, 37, 50, 39, 27 ]
    geom.face(37).n_poi =  6; allocate(geom.face(37).poi(6));  geom.face(37).poi(1: 6) = [  24, 36, 48, 54, 42, 29 ]
    geom.face(38).n_poi =  6; allocate(geom.face(38).poi(6));  geom.face(38).poi(1: 6) = [  32, 38, 51, 63, 57, 45 ]
    geom.face(39).n_poi =  6; allocate(geom.face(39).poi(6));  geom.face(39).poi(1: 6) = [  34, 47, 59, 71, 58, 46 ]
    geom.face(40).n_poi =  6; allocate(geom.face(40).poi(6));  geom.face(40).poi(1: 6) = [  40, 52, 64, 77, 65, 53 ]
    geom.face(41).n_poi =  6; allocate(geom.face(41).poi(6));  geom.face(41).poi(1: 6) = [  43, 55, 67, 80, 68, 56 ]
    geom.face(42).n_poi =  6; allocate(geom.face(42).poi(6));  geom.face(42).poi(1: 6) = [  49, 62, 74, 86, 73, 61 ]
    geom.face(43).n_poi =  6; allocate(geom.face(43).poi(6));  geom.face(43).poi(1: 6) = [  60, 72, 84, 90, 78, 66 ]
    geom.face(44).n_poi =  6; allocate(geom.face(44).poi(6));  geom.face(44).poi(1: 6) = [  69, 75, 87, 98, 93, 81 ]
    geom.face(45).n_poi =  6; allocate(geom.face(45).poi(6));  geom.face(45).poi(1: 6) = [  70, 82, 94,104, 95, 83 ]
    geom.face(46).n_poi =  6; allocate(geom.face(46).poi(6));  geom.face(46).poi(1: 6) = [  76, 89,100,108, 99, 88 ]
    geom.face(47).n_poi =  6; allocate(geom.face(47).poi(6));  geom.face(47).poi(1: 6) = [  79, 92,102,110,101, 91 ]
    geom.face(48).n_poi =  6; allocate(geom.face(48).poi(6));  geom.face(48).poi(1: 6) = [  85, 96,105,113,106, 97 ]
    geom.face(49).n_poi =  6; allocate(geom.face(49).poi(6));  geom.face(49).poi(1: 6) = [ 103,109,115,119,116,111 ]
    geom.face(50).n_poi =  6; allocate(geom.face(50).poi(6));  geom.face(50).poi(1: 6) = [ 107,114,118,120,117,112 ]
    geom.face(51).n_poi = 10; allocate(geom.face(51).poi(10)); geom.face(51).poi(1:10) = [   1,  4,  8, 16, 24, 29, 20, 11,  6,  2 ]
    geom.face(52).n_poi = 10; allocate(geom.face(52).poi(10)); geom.face(52).poi(1:10) = [   3,  5, 10, 17, 27, 38, 32, 22, 13,  7 ]
    geom.face(53).n_poi = 10; allocate(geom.face(53).poi(10)); geom.face(53).poi(1:10) = [   9, 14, 21, 31, 43, 56, 47, 34, 23, 15 ]
    geom.face(54).n_poi = 10; allocate(geom.face(54).poi(10)); geom.face(54).poi(1:10) = [  12, 19, 28, 40, 53, 62, 49, 37, 26, 18 ]
    geom.face(55).n_poi = 10; allocate(geom.face(55).poi(10)); geom.face(55).poi(1:10) = [  25, 35, 46, 58, 70, 83, 72, 60, 48, 36 ]
    geom.face(56).n_poi = 10; allocate(geom.face(56).poi(10)); geom.face(56).poi(1:10) = [  30, 42, 54, 66, 78, 89, 76, 64, 52, 41 ]
    geom.face(57).n_poi = 10; allocate(geom.face(57).poi(10)); geom.face(57).poi(1:10) = [  33, 45, 57, 69, 81, 92, 79, 67, 55, 44 ]
    geom.face(58).n_poi = 10; allocate(geom.face(58).poi(10)); geom.face(58).poi(1:10) = [  39, 50, 61, 73, 85, 97, 87, 75, 63, 51 ]
    geom.face(59).n_poi = 10; allocate(geom.face(59).poi(10)); geom.face(59).poi(1:10) = [  59, 68, 80, 91,101,109,103, 94, 82, 71 ]
    geom.face(60).n_poi = 10; allocate(geom.face(60).poi(10)); geom.face(60).poi(1:10) = [  65, 77, 88, 99,107,112,105, 96, 86, 74 ]
    geom.face(61).n_poi = 10; allocate(geom.face(61).poi(10)); geom.face(61).poi(1:10) = [  84, 95,104,111,116,118,114,108,100, 90 ]
    geom.face(62).n_poi = 10; allocate(geom.face(62).poi(10)); geom.face(62).poi(1:10) = [  93, 98,106,113,117,120,119,115,110,102 ]
end subroutine Exam_Archi_Truncated_Icosidodecahedron

! ---------------------------------------------------------------------------------------

! Example of Rhombicosidodecahedron (Rhombicosidodecahedron, Small Rhombicosidodecahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Rhombicosidodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "54_Rhombicosidodecahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Rhombicosidodecahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 60
    geom.n_face = 62

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  22.32948d0 ]; geom.iniP( 2).pos(1:3) = [   9.74610d0,   0.00000d0,  20.09037d0 ]
    geom.iniP( 3).pos(1:3) = [  -0.51446d0,   9.73252d0,  20.09037d0 ]; geom.iniP( 4).pos(1:3) = [  -8.81543d0,   4.15628d0,  20.09037d0 ]
    geom.iniP( 5).pos(1:3) = [  -3.68514d0,  -9.02254d0,  20.09037d0 ]; geom.iniP( 6).pos(1:3) = [   9.23164d0,   9.73252d0,  17.85118d0 ]
    geom.iniP( 7).pos(1:3) = [  16.70020d0,   4.15628d0,  14.22808d0 ]; geom.iniP( 8).pos(1:3) = [  12.08438d0,  -9.02254d0,  16.46727d0 ]
    geom.iniP( 9).pos(1:3) = [  -5.03203d0,  16.45752d0,  14.22808d0 ]; geom.iniP(10).pos(1:3) = [ -12.50059d0,  -4.86626d0,  17.85118d0 ]
    geom.iniP(11).pos(1:3) = [ -13.33301d0,  10.88128d0,  14.22808d0 ]; geom.iniP(12).pos(1:3) = [   3.78340d0, -14.59876d0,  16.46727d0 ]
    geom.iniP(13).pos(1:3) = [ -10.16231d0, -13.88880d0,  14.22808d0 ]; geom.iniP(14).pos(1:3) = [  10.73750d0,  16.45752d0,  10.60498d0 ]
    geom.iniP(15).pos(1:3) = [  19.03846d0,  -4.86626d0,  10.60498d0 ]; geom.iniP(16).pos(1:3) = [  18.20606d0,  10.88128d0,   6.98189d0 ]
    geom.iniP(17).pos(1:3) = [  15.35331d0, -13.88880d0,   8.36579d0 ]; geom.iniP(18).pos(1:3) = [   1.92206d0,  20.61380d0,   8.36579d0 ]
    geom.iniP(19).pos(1:3) = [ -11.82715d0,  17.60628d0,   6.98189d0 ]; geom.iniP(20).pos(1:3) = [ -19.29569d0,  -3.71748d0,  10.60498d0 ]
    geom.iniP(21).pos(1:3) = [ -19.81016d0,   6.01502d0,   8.36579d0 ]; geom.iniP(22).pos(1:3) = [  -2.69375d0, -19.46502d0,  10.60498d0 ]
    geom.iniP(23).pos(1:3) = [   7.05235d0, -19.46502d0,   8.36579d0 ]; geom.iniP(24).pos(1:3) = [ -16.95743d0, -12.74002d0,   6.98189d0 ]
    geom.iniP(25).pos(1:3) = [  13.68848d0,  17.60628d0,   1.11960d0 ]; geom.iniP(26).pos(1:3) = [  21.98955d0,  -3.71748d0,   1.11960d0 ]
    geom.iniP(27).pos(1:3) = [  21.47499d0,   6.01502d0,  -1.11960d0 ]; geom.iniP(28).pos(1:3) = [  18.30430d0, -12.74002d0,  -1.11960d0 ]
    geom.iniP(29).pos(1:3) = [  -4.87305d0,  21.76256d0,   1.11960d0 ]; geom.iniP(30).pos(1:3) = [   4.87305d0,  21.76256d0,  -1.11960d0 ]
    geom.iniP(31).pos(1:3) = [ -18.30430d0,  12.74002d0,   1.11960d0 ]; geom.iniP(32).pos(1:3) = [ -21.47499d0,  -6.01502d0,   1.11960d0 ]
    geom.iniP(33).pos(1:3) = [ -21.98955d0,   3.71748d0,  -1.11960d0 ]; geom.iniP(34).pos(1:3) = [  -4.87305d0, -21.76256d0,   1.11960d0 ]
    geom.iniP(35).pos(1:3) = [   4.87305d0, -21.76256d0,  -1.11960d0 ]; geom.iniP(36).pos(1:3) = [ -13.68848d0, -17.60628d0,  -1.11960d0 ]
    geom.iniP(37).pos(1:3) = [  16.95743d0,  12.74002d0,  -6.98189d0 ]; geom.iniP(38).pos(1:3) = [  19.81016d0,  -6.01502d0,  -8.36579d0 ]
    geom.iniP(39).pos(1:3) = [  19.29569d0,   3.71748d0, -10.60498d0 ]; geom.iniP(40).pos(1:3) = [  11.82715d0, -17.60628d0,  -6.98189d0 ]
    geom.iniP(41).pos(1:3) = [  -7.05235d0,  19.46502d0,  -8.36579d0 ]; geom.iniP(42).pos(1:3) = [   2.69375d0,  19.46502d0, -10.60498d0 ]
    geom.iniP(43).pos(1:3) = [ -15.35331d0,  13.88880d0,  -8.36579d0 ]; geom.iniP(44).pos(1:3) = [ -18.20606d0, -10.88128d0,  -6.98189d0 ]
    geom.iniP(45).pos(1:3) = [ -19.03846d0,   4.86626d0, -10.60498d0 ]; geom.iniP(46).pos(1:3) = [  -1.92206d0, -20.61380d0,  -8.36579d0 ]
    geom.iniP(47).pos(1:3) = [ -10.73750d0, -16.45752d0, -10.60498d0 ]; geom.iniP(48).pos(1:3) = [  10.16231d0,  13.88880d0, -14.22808d0 ]
    geom.iniP(49).pos(1:3) = [  13.33301d0, -10.88128d0, -14.22808d0 ]; geom.iniP(50).pos(1:3) = [  12.50059d0,   4.86626d0, -17.85118d0 ]
    geom.iniP(51).pos(1:3) = [   5.03203d0, -16.45752d0, -14.22808d0 ]; geom.iniP(52).pos(1:3) = [  -3.78340d0,  14.59876d0, -16.46727d0 ]
    geom.iniP(53).pos(1:3) = [ -12.08438d0,   9.02254d0, -16.46727d0 ]; geom.iniP(54).pos(1:3) = [ -16.70020d0,  -4.15628d0, -14.22808d0 ]
    geom.iniP(55).pos(1:3) = [  -9.23164d0,  -9.73252d0, -17.85118d0 ]; geom.iniP(56).pos(1:3) = [   3.68514d0,   9.02254d0, -20.09037d0 ]
    geom.iniP(57).pos(1:3) = [   8.81543d0,  -4.15628d0, -20.09037d0 ]; geom.iniP(58).pos(1:3) = [   0.51446d0,  -9.73252d0, -20.09037d0 ]
    geom.iniP(59).pos(1:3) = [  -9.74610d0,   0.00000d0, -20.09037d0 ]; geom.iniP(60).pos(1:3) = [   0.00000d0,   0.00000d0, -22.32948d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  3,  4 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  2,  7,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  5, 10, 13 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8, 17, 15 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  9, 19, 11 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 12, 22, 23 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 14, 16, 25 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 18, 30, 29 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 20, 32, 24 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 21, 31, 33 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 26, 28, 38 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 27, 39, 37 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 34, 46, 35 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 36, 44, 47 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 40, 51, 49 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 41, 42, 52 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 43, 53, 45 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 48, 50, 56 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 54, 59, 55 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 57, 58, 60 ]
    geom.face(21).n_poi = 4; allocate(geom.face(21).poi(4)); geom.face(21).poi(1:4) = [  1,  2,  6,  3 ]
    geom.face(22).n_poi = 4; allocate(geom.face(22).poi(4)); geom.face(22).poi(1:4) = [  1,  4, 10,  5 ]
    geom.face(23).n_poi = 4; allocate(geom.face(23).poi(4)); geom.face(23).poi(1:4) = [  2,  8, 15,  7 ]
    geom.face(24).n_poi = 4; allocate(geom.face(24).poi(4)); geom.face(24).poi(1:4) = [  3,  9, 11,  4 ]
    geom.face(25).n_poi = 4; allocate(geom.face(25).poi(4)); geom.face(25).poi(1:4) = [  5, 13, 22, 12 ]
    geom.face(26).n_poi = 4; allocate(geom.face(26).poi(4)); geom.face(26).poi(1:4) = [  6,  7, 16, 14 ]
    geom.face(27).n_poi = 4; allocate(geom.face(27).poi(4)); geom.face(27).poi(1:4) = [  8, 12, 23, 17 ]
    geom.face(28).n_poi = 4; allocate(geom.face(28).poi(4)); geom.face(28).poi(1:4) = [  9, 18, 29, 19 ]
    geom.face(29).n_poi = 4; allocate(geom.face(29).poi(4)); geom.face(29).poi(1:4) = [ 10, 20, 24, 13 ]
    geom.face(30).n_poi = 4; allocate(geom.face(30).poi(4)); geom.face(30).poi(1:4) = [ 11, 19, 31, 21 ]
    geom.face(31).n_poi = 4; allocate(geom.face(31).poi(4)); geom.face(31).poi(1:4) = [ 14, 25, 30, 18 ]
    geom.face(32).n_poi = 4; allocate(geom.face(32).poi(4)); geom.face(32).poi(1:4) = [ 15, 17, 28, 26 ]
    geom.face(33).n_poi = 4; allocate(geom.face(33).poi(4)); geom.face(33).poi(1:4) = [ 16, 27, 37, 25 ]
    geom.face(34).n_poi = 4; allocate(geom.face(34).poi(4)); geom.face(34).poi(1:4) = [ 20, 21, 33, 32 ]
    geom.face(35).n_poi = 4; allocate(geom.face(35).poi(4)); geom.face(35).poi(1:4) = [ 22, 34, 35, 23 ]
    geom.face(36).n_poi = 4; allocate(geom.face(36).poi(4)); geom.face(36).poi(1:4) = [ 24, 32, 44, 36 ]
    geom.face(37).n_poi = 4; allocate(geom.face(37).poi(4)); geom.face(37).poi(1:4) = [ 26, 38, 39, 27 ]
    geom.face(38).n_poi = 4; allocate(geom.face(38).poi(4)); geom.face(38).poi(1:4) = [ 28, 40, 49, 38 ]
    geom.face(39).n_poi = 4; allocate(geom.face(39).poi(4)); geom.face(39).poi(1:4) = [ 29, 30, 42, 41 ]
    geom.face(40).n_poi = 4; allocate(geom.face(40).poi(4)); geom.face(40).poi(1:4) = [ 31, 43, 45, 33 ]
    geom.face(41).n_poi = 4; allocate(geom.face(41).poi(4)); geom.face(41).poi(1:4) = [ 34, 36, 47, 46 ]
    geom.face(42).n_poi = 4; allocate(geom.face(42).poi(4)); geom.face(42).poi(1:4) = [ 35, 46, 51, 40 ]
    geom.face(43).n_poi = 4; allocate(geom.face(43).poi(4)); geom.face(43).poi(1:4) = [ 37, 39, 50, 48 ]
    geom.face(44).n_poi = 4; allocate(geom.face(44).poi(4)); geom.face(44).poi(1:4) = [ 41, 52, 53, 43 ]
    geom.face(45).n_poi = 4; allocate(geom.face(45).poi(4)); geom.face(45).poi(1:4) = [ 42, 48, 56, 52 ]
    geom.face(46).n_poi = 4; allocate(geom.face(46).poi(4)); geom.face(46).poi(1:4) = [ 44, 54, 55, 47 ]
    geom.face(47).n_poi = 4; allocate(geom.face(47).poi(4)); geom.face(47).poi(1:4) = [ 45, 53, 59, 54 ]
    geom.face(48).n_poi = 4; allocate(geom.face(48).poi(4)); geom.face(48).poi(1:4) = [ 49, 51, 58, 57 ]
    geom.face(49).n_poi = 4; allocate(geom.face(49).poi(4)); geom.face(49).poi(1:4) = [ 50, 57, 60, 56 ]
    geom.face(50).n_poi = 4; allocate(geom.face(50).poi(4)); geom.face(50).poi(1:4) = [ 55, 59, 60, 58 ]
    geom.face(51).n_poi = 5; allocate(geom.face(51).poi(5)); geom.face(51).poi(1:5) = [  1,  5, 12,  8,  2 ]
    geom.face(52).n_poi = 5; allocate(geom.face(52).poi(5)); geom.face(52).poi(1:5) = [  3,  6, 14, 18,  9 ]
    geom.face(53).n_poi = 5; allocate(geom.face(53).poi(5)); geom.face(53).poi(1:5) = [  4, 11, 21, 20, 10 ]
    geom.face(54).n_poi = 5; allocate(geom.face(54).poi(5)); geom.face(54).poi(1:5) = [  7, 15, 26, 27, 16 ]
    geom.face(55).n_poi = 5; allocate(geom.face(55).poi(5)); geom.face(55).poi(1:5) = [ 13, 24, 36, 34, 22 ]
    geom.face(56).n_poi = 5; allocate(geom.face(56).poi(5)); geom.face(56).poi(1:5) = [ 17, 23, 35, 40, 28 ]
    geom.face(57).n_poi = 5; allocate(geom.face(57).poi(5)); geom.face(57).poi(1:5) = [ 19, 29, 41, 43, 31 ]
    geom.face(58).n_poi = 5; allocate(geom.face(58).poi(5)); geom.face(58).poi(1:5) = [ 25, 37, 48, 42, 30 ]
    geom.face(59).n_poi = 5; allocate(geom.face(59).poi(5)); geom.face(59).poi(1:5) = [ 32, 33, 45, 54, 44 ]
    geom.face(60).n_poi = 5; allocate(geom.face(60).poi(5)); geom.face(60).poi(1:5) = [ 38, 49, 57, 50, 39 ]
    geom.face(61).n_poi = 5; allocate(geom.face(61).poi(5)); geom.face(61).poi(1:5) = [ 46, 47, 55, 58, 51 ]
    geom.face(62).n_poi = 5; allocate(geom.face(62).poi(5)); geom.face(62).poi(1:5) = [ 52, 56, 60, 59, 53 ]
end subroutine Exam_Archi_Rhombicosidodecahedron

! ---------------------------------------------------------------------------------------

! Example of Snub Dodecahedron (Snub Icosidodecahedron)
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Archi_Snub_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "55_Snub_Dodecahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Snub Dodecahedron"

    ! Set geometric type and view
    prob.color    = [231, 76, 60]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.94d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 60
    geom.n_face = 92

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  21.55838d0 ]; geom.iniP( 2).pos(1:3) = [   9.72735d0,   0.00000d0,  19.23912d0 ]
    geom.iniP( 3).pos(1:3) = [   4.58719d0,   8.57783d0,  19.23912d0 ]; geom.iniP( 4).pos(1:3) = [  -5.40095d0,   8.09018d0,  19.23912d0 ]
    geom.iniP( 5).pos(1:3) = [  -9.68109d0,  -0.94756d0,  19.23912d0 ]; geom.iniP( 6).pos(1:3) = [  -3.72978d0,  -8.98388d0,  19.23912d0 ]
    geom.iniP( 7).pos(1:3) = [  12.00939d0,  -8.98388d0,  15.48644d0 ]; geom.iniP( 8).pos(1:3) = [  17.32044d0,  -0.94756d0,  12.80115d0 ]
    geom.iniP( 9).pos(1:3) = [  13.50075d0,   8.09018d0,  14.73240d0 ]; geom.iniP(10).pos(1:3) = [   7.33232d0,  15.72045d0,  12.80115d0 ]
    geom.iniP(11).pos(1:3) = [  -8.82880d0,  14.93143d0,  12.80115d0 ]; geom.iniP(12).pos(1:3) = [ -14.22463d0,   6.73659d0,  14.73240d0 ]
    geom.iniP(13).pos(1:3) = [ -12.64849d0,  -9.36684d0,  14.73240d0 ]; geom.iniP(14).pos(1:3) = [  -5.76609d0, -16.35994d0,  12.80115d0 ]
    geom.iniP(15).pos(1:3) = [   3.69244d0, -14.53622d0,  15.48644d0 ]; geom.iniP(16).pos(1:3) = [  10.92177d0, -16.35994d0,   8.82228d0 ]
    geom.iniP(17).pos(1:3) = [  17.93513d0,  -9.36684d0,   7.44036d0 ]; geom.iniP(18).pos(1:3) = [  21.24194d0,  -0.74416d0,   3.60438d0 ]
    geom.iniP(19).pos(1:3) = [  15.06157d0,  13.87921d0,   6.72919d0 ]; geom.iniP(20).pos(1:3) = [   7.62006d0,  19.82678d0,   3.68768d0 ]
    geom.iniP(21).pos(1:3) = [  -0.95921d0,  19.64721d0,   8.82228d0 ]; geom.iniP(22).pos(1:3) = [  -9.51517d0,  18.99020d0,   3.68768d0 ]
    geom.iniP(23).pos(1:3) = [ -16.34193d0,  12.34603d0,   6.72919d0 ]; geom.iniP(24).pos(1:3) = [ -20.00010d0,   3.06636d0,   7.44036d0 ]
    geom.iniP(25).pos(1:3) = [ -19.02597d0,  -6.88611d0,   7.44036d0 ]; geom.iniP(26).pos(1:3) = [ -13.63794d0, -15.28038d0,   6.72919d0 ]
    geom.iniP(27).pos(1:3) = [   1.77400d0, -20.15587d0,   7.44036d0 ]; geom.iniP(28).pos(1:3) = [   8.28160d0, -19.90372d0,  -0.14832d0 ]
    geom.iniP(29).pos(1:3) = [  19.62945d0,  -8.58865d0,  -2.38431d0 ]; geom.iniP(30).pos(1:3) = [  20.61278d0,   0.56459d0,  -6.28964d0 ]
    geom.iniP(31).pos(1:3) = [  19.84588d0,   8.41929d0,  -0.14832d0 ]; geom.iniP(32).pos(1:3) = [  14.05989d0,  16.06940d0,  -2.97650d0 ]
    geom.iniP(33).pos(1:3) = [  -1.05013d0,  21.50952d0,  -1.00239d0 ]; geom.iniP(34).pos(1:3) = [  -9.00548d0,  18.55011d0,  -6.28964d0 ]
    geom.iniP(35).pos(1:3) = [ -20.05139d0,   7.79963d0,  -1.36838d0 ]; geom.iniP(36).pos(1:3) = [ -21.43277d0,  -2.09777d0,  -1.00239d0 ]
    geom.iniP(37).pos(1:3) = [ -18.15855d0, -11.53950d0,  -1.36838d0 ]; geom.iniP(38).pos(1:3) = [ -10.96291d0, -18.40910d0,  -2.38431d0 ]
    geom.iniP(39).pos(1:3) = [  -1.43781d0, -21.42228d0,  -1.94478d0 ]; geom.iniP(40).pos(1:3) = [   4.68074d0, -18.82084d0,  -9.41447d0 ]
    geom.iniP(41).pos(1:3) = [  13.66326d0, -15.10081d0,  -7.07436d0 ]; geom.iniP(42).pos(1:3) = [  16.71257d0,  -6.82374d0, -11.78523d0 ]
    geom.iniP(43).pos(1:3) = [  16.91669d0,   9.21730d0,  -9.67635d0 ]; geom.iniP(44).pos(1:3) = [   9.36974d0,  15.42990d0, -11.78523d0 ]
    geom.iniP(45).pos(1:3) = [   0.03122d0,  18.79205d0, -10.56517d0 ]; geom.iniP(46).pos(1:3) = [  -7.41394d0,  13.70985d0, -14.89425d0 ]
    geom.iniP(47).pos(1:3) = [ -15.51722d0,  11.63397d0,  -9.41447d0 ]; geom.iniP(48).pos(1:3) = [ -19.12567d0,   2.31138d0,  -9.67635d0 ]
    geom.iniP(49).pos(1:3) = [ -13.82800d0, -12.96569d0, -10.26854d0 ]; geom.iniP(50).pos(1:3) = [  -5.22304d0, -18.05166d0, -10.56517d0 ]
    geom.iniP(51).pos(1:3) = [   9.01955d0, -12.08952d0, -15.40314d0 ]; geom.iniP(52).pos(1:3) = [  10.60604d0,  -2.73725d0, -18.56838d0 ]
    geom.iniP(53).pos(1:3) = [  10.73220d0,   7.17666d0, -17.26501d0 ]; geom.iniP(54).pos(1:3) = [   2.04245d0,  12.10409d0, -17.72242d0 ]
    geom.iniP(55).pos(1:3) = [ -12.62208d0,   5.35293d0, -16.63713d0 ]; geom.iniP(56).pos(1:3) = [ -14.42570d0,  -4.40537d0, -15.40314d0 ]
    geom.iniP(57).pos(1:3) = [  -7.00509d0, -10.84495d0, -17.26501d0 ]; geom.iniP(58).pos(1:3) = [   1.79731d0,  -7.16015d0, -20.25506d0 ]
    geom.iniP(59).pos(1:3) = [   2.67868d0,   2.75476d0, -21.21321d0 ]; geom.iniP(60).pos(1:3) = [  -6.38449d0, -1.41768d0,  -20.54248d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  3,  4 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  1,  4,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  1,  5,  6 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  2,  7,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  2,  8,  9 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  2,  9,  3 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  3,  9, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4, 11, 12 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4, 12,  5 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  5, 13,  6 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  6, 13, 14 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  6, 14, 15 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  7, 15, 16 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  7, 16, 17 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  7, 17,  8 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  8, 17, 18 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  9, 19, 10 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 10, 19, 20 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 10, 20, 21 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 11, 21, 22 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 11, 22, 23 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 11, 23, 12 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 12, 23, 24 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 13, 25, 26 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 13, 26, 14 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 14, 27, 15 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 15, 27, 16 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 16, 27, 28 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 17, 29, 18 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 18, 29, 30 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 18, 30, 31 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 19, 31, 32 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 19, 32, 20 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 20, 33, 21 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 21, 33, 22 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 22, 33, 34 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 23, 35, 24 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 24, 35, 36 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 24, 36, 25 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 25, 36, 37 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 25, 37, 26 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 26, 37, 38 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 27, 39, 28 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 28, 39, 40 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 28, 40, 41 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 29, 41, 42 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 29, 42, 30 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 30, 43, 31 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 31, 43, 32 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 32, 43, 44 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 33, 45, 34 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 34, 45, 46 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 34, 46, 47 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 35, 47, 48 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 35, 48, 36 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 37, 49, 38 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 38, 49, 50 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 38, 50, 39 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 39, 50, 40 ]
    geom.face(61).n_poi = 3; allocate(geom.face(61).poi(3)); geom.face(61).poi(1:3) = [ 40, 51, 41 ]
    geom.face(62).n_poi = 3; allocate(geom.face(62).poi(3)); geom.face(62).poi(1:3) = [ 41, 51, 42 ]
    geom.face(63).n_poi = 3; allocate(geom.face(63).poi(3)); geom.face(63).poi(1:3) = [ 42, 51, 52 ]
    geom.face(64).n_poi = 3; allocate(geom.face(64).poi(3)); geom.face(64).poi(1:3) = [ 43, 53, 44 ]
    geom.face(65).n_poi = 3; allocate(geom.face(65).poi(3)); geom.face(65).poi(1:3) = [ 44, 53, 54 ]
    geom.face(66).n_poi = 3; allocate(geom.face(66).poi(3)); geom.face(66).poi(1:3) = [ 44, 54, 45 ]
    geom.face(67).n_poi = 3; allocate(geom.face(67).poi(3)); geom.face(67).poi(1:3) = [ 45, 54, 46 ]
    geom.face(68).n_poi = 3; allocate(geom.face(68).poi(3)); geom.face(68).poi(1:3) = [ 46, 55, 47 ]
    geom.face(69).n_poi = 3; allocate(geom.face(69).poi(3)); geom.face(69).poi(1:3) = [ 47, 55, 48 ]
    geom.face(70).n_poi = 3; allocate(geom.face(70).poi(3)); geom.face(70).poi(1:3) = [ 48, 55, 56 ]
    geom.face(71).n_poi = 3; allocate(geom.face(71).poi(3)); geom.face(71).poi(1:3) = [ 49, 56, 57 ]
    geom.face(72).n_poi = 3; allocate(geom.face(72).poi(3)); geom.face(72).poi(1:3) = [ 49, 57, 50 ]
    geom.face(73).n_poi = 3; allocate(geom.face(73).poi(3)); geom.face(73).poi(1:3) = [ 51, 58, 52 ]
    geom.face(74).n_poi = 3; allocate(geom.face(74).poi(3)); geom.face(74).poi(1:3) = [ 52, 58, 59 ]
    geom.face(75).n_poi = 3; allocate(geom.face(75).poi(3)); geom.face(75).poi(1:3) = [ 52, 59, 53 ]
    geom.face(76).n_poi = 3; allocate(geom.face(76).poi(3)); geom.face(76).poi(1:3) = [ 53, 59, 54 ]
    geom.face(77).n_poi = 3; allocate(geom.face(77).poi(3)); geom.face(77).poi(1:3) = [ 55, 60, 56 ]
    geom.face(78).n_poi = 3; allocate(geom.face(78).poi(3)); geom.face(78).poi(1:3) = [ 56, 60, 57 ]
    geom.face(79).n_poi = 3; allocate(geom.face(79).poi(3)); geom.face(79).poi(1:3) = [ 57, 60, 58 ]
    geom.face(80).n_poi = 3; allocate(geom.face(80).poi(3)); geom.face(80).poi(1:3) = [ 58, 60, 59 ]
    geom.face(81).n_poi = 5; allocate(geom.face(81).poi(5)); geom.face(81).poi(1:5) = [  1,  6, 15,  7,  2 ]
    geom.face(82).n_poi = 5; allocate(geom.face(82).poi(5)); geom.face(82).poi(1:5) = [  3, 10, 21, 11,  4 ]
    geom.face(83).n_poi = 5; allocate(geom.face(83).poi(5)); geom.face(83).poi(1:5) = [  5, 12, 24, 25, 13 ]
    geom.face(84).n_poi = 5; allocate(geom.face(84).poi(5)); geom.face(84).poi(1:5) = [  8, 18, 31, 19,  9 ]
    geom.face(85).n_poi = 5; allocate(geom.face(85).poi(5)); geom.face(85).poi(1:5) = [ 14, 26, 38, 39, 27 ]
    geom.face(86).n_poi = 5; allocate(geom.face(86).poi(5)); geom.face(86).poi(1:5) = [ 16, 28, 41, 29, 17 ]
    geom.face(87).n_poi = 5; allocate(geom.face(87).poi(5)); geom.face(87).poi(1:5) = [ 20, 32, 44, 45, 33 ]
    geom.face(88).n_poi = 5; allocate(geom.face(88).poi(5)); geom.face(88).poi(1:5) = [ 22, 34, 47, 35, 23 ]
    geom.face(89).n_poi = 5; allocate(geom.face(89).poi(5)); geom.face(89).poi(1:5) = [ 30, 42, 52, 53, 43 ]
    geom.face(90).n_poi = 5; allocate(geom.face(90).poi(5)); geom.face(90).poi(1:5) = [ 36, 48, 56, 49, 37 ]
    geom.face(91).n_poi = 5; allocate(geom.face(91).poi(5)); geom.face(91).poi(1:5) = [ 40, 50, 57, 58, 51 ]
    geom.face(92).n_poi = 5; allocate(geom.face(92).poi(5)); geom.face(92).poi(1:5) = [ 46, 54, 59, 60, 55 ]
end subroutine Exam_Archi_Snub_Dodecahedron

! ---------------------------------------------------------------------------------------

end module Exam_Archi