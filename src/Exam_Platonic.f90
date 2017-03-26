!
! ---------------------------------------------------------------------------------------
!
!                               Module for Exam_Platonic
!
!                                                             First modified : 2015/08/03
!                                                             Last  modified : 2016/07/14
!
! Comments: This module contains the geometric definition of platonic solids
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
module Exam_Platonic

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Platonic_Tetrahedron        ! 1. V=4,  E=6,  F=4
    public Exam_Platonic_Cube               ! 2. V=8,  E=12, F=6
    public Exam_Platonic_Octahedron         ! 3. V=6,  E=12, F=8
    public Exam_Platonic_Dodecahedron       ! 4. V=20, E=30, F=12
    public Exam_Platonic_Icosahedron        ! 5. V=12, E=30, F=20

contains

! ---------------------------------------------------------------------------------------

! Example of Tetrahedron
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Exam_Platonic_Tetrahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "01_Tetrahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tetrahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 4
    geom.n_face = 4

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.00000d0,  0.00000d0,  6.12374d0 ]
    geom.iniP(2).pos(1:3) = [  5.77351d0,  0.00000d0, -2.04125d0 ]
    geom.iniP(3).pos(1:3) = [ -2.88676d0,  5.00000d0, -2.04125d0 ]
    geom.iniP(4).pos(1:3) = [ -2.88676d0, -5.00000d0, -2.04125d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 4, 2 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 2, 4, 3 ]
end subroutine Exam_Platonic_Tetrahedron

! ---------------------------------------------------------------------------------------

! Example of Cube
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Platonic_Cube(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "02_Cube"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Cube"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 0.7d0      ! Atomic model
    prob.size     = 1.00d0     ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    para_fig_view = "XYZ"

    ! allocate point, line and face structure
    geom.n_iniP = 8
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [ -1.0d0, -1.0d0, -1.0d0 ]; geom.iniP(2).pos(1:3) = [  1.0d0, -1.0d0, -1.0d0 ]
    geom.iniP(3).pos(1:3) = [  1.0d0,  1.0d0, -1.0d0 ]; geom.iniP(4).pos(1:3) = [ -1.0d0,  1.0d0, -1.0d0 ]
    geom.iniP(5).pos(1:3) = [ -1.0d0, -1.0d0,  1.0d0 ]; geom.iniP(6).pos(1:3) = [  1.0d0, -1.0d0,  1.0d0 ]
    geom.iniP(7).pos(1:3) = [  1.0d0,  1.0d0,  1.0d0 ]; geom.iniP(8).pos(1:3) = [ -1.0d0,  1.0d0,  1.0d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 1, 4, 3, 2 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 5, 6, 7, 8 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 2, 3, 7, 6 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 1, 5, 8, 4 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 1, 2, 6, 5 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 3, 4, 8, 7 ]
end subroutine Exam_Platonic_Cube

! ---------------------------------------------------------------------------------------

! Example of Octahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Platonic_Octahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "03_Octahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Octahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = -1.0d0     ! Cylindrical model
    prob.move_y   = -1.0d0     ! Cylindrical model
    para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 6
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  0.00000d0,  0.00000d0,  7.07107d0 ]
    geom.iniP(2).pos(1:3) = [  7.07107d0,  0.00000d0,  0.00000d0 ]
    geom.iniP(3).pos(1:3) = [  0.00000d0,  7.07107d0,  0.00000d0 ]
    geom.iniP(4).pos(1:3) = [ -7.07107d0,  0.00000d0,  0.00000d0 ]
    geom.iniP(5).pos(1:3) = [  0.00000d0, -7.07107d0,  0.00000d0 ]
    geom.iniP(6).pos(1:3) = [  0.00000d0,  0.00000d0, -7.07107d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi = 3; allocate(geom.face(1).poi(3)); geom.face(1).poi(1:3) = [ 1, 2, 3 ]
    geom.face(2).n_poi = 3; allocate(geom.face(2).poi(3)); geom.face(2).poi(1:3) = [ 1, 3, 4 ]
    geom.face(3).n_poi = 3; allocate(geom.face(3).poi(3)); geom.face(3).poi(1:3) = [ 1, 4, 5 ]
    geom.face(4).n_poi = 3; allocate(geom.face(4).poi(3)); geom.face(4).poi(1:3) = [ 1, 5, 2 ]
    geom.face(5).n_poi = 3; allocate(geom.face(5).poi(3)); geom.face(5).poi(1:3) = [ 2, 5, 6 ]
    geom.face(6).n_poi = 3; allocate(geom.face(6).poi(3)); geom.face(6).poi(1:3) = [ 2, 6, 3 ]
    geom.face(7).n_poi = 3; allocate(geom.face(7).poi(3)); geom.face(7).poi(1:3) = [ 3, 6, 4 ]
    geom.face(8).n_poi = 3; allocate(geom.face(8).poi(3)); geom.face(8).poi(1:3) = [ 4, 6, 5 ]
end subroutine Exam_Platonic_Octahedron

! ---------------------------------------------------------------------------------------

! Example of Dodecahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Platonic_Dodecahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp
    
    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "04_Dodecahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Dodecahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 0.95d0     ! Cylindrical model
    prob.move_x   = 3.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 20
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   0.00000d0,   0.00000d0,  14.01264d0 ]; geom.iniP( 2).pos(1:3) = [   9.34173d0,   0.00000d0,  10.44437d0 ]
    geom.iniP( 3).pos(1:3) = [  -4.67086d0,   8.09018d0,  10.44437d0 ]; geom.iniP( 4).pos(1:3) = [  -4.67086d0,  -8.09018d0,  10.44437d0 ]
    geom.iniP( 5).pos(1:3) = [  10.44437d0,   8.09018d0,   4.67086d0 ]; geom.iniP( 6).pos(1:3) = [  10.44437d0,  -8.09018d0,   4.67086d0 ]
    geom.iniP( 7).pos(1:3) = [ -12.22848d0,   5.00000d0,   4.67086d0 ]; geom.iniP( 8).pos(1:3) = [   1.78411d0,  13.09018d0,   4.67086d0 ]
    geom.iniP( 9).pos(1:3) = [   1.78411d0, -13.09018d0,   4.67086d0 ]; geom.iniP(10).pos(1:3) = [ -12.22848d0,  -5.00000d0,   4.67086d0 ]
    geom.iniP(11).pos(1:3) = [  12.22848d0,   5.00000d0,  -4.67086d0 ]; geom.iniP(12).pos(1:3) = [  12.22848d0,  -5.00000d0,  -4.67086d0 ]
    geom.iniP(13).pos(1:3) = [ -10.44437d0,   8.09018d0,  -4.67086d0 ]; geom.iniP(14).pos(1:3) = [  -1.78411d0,  13.09018d0,  -4.67086d0 ]
    geom.iniP(15).pos(1:3) = [  -1.78411d0, -13.09018d0,  -4.67086d0 ]; geom.iniP(16).pos(1:3) = [ -10.44437d0,  -8.09018d0,  -4.67086d0 ]
    geom.iniP(17).pos(1:3) = [   4.67086d0,   8.09018d0, -10.44437d0 ]; geom.iniP(18).pos(1:3) = [   4.67086d0,  -8.09018d0, -10.44437d0 ]
    geom.iniP(19).pos(1:3) = [  -9.34173d0,   0.00000d0, -10.44437d0 ]; geom.iniP(20).pos(1:3) = [   0.00000d0,   0.00000d0, -14.01264d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 5; allocate(geom.face(1).poi(5));  geom.face( 1).poi(1:5) = [  1,  2,  5,  8,  3 ]
    geom.face(2).n_poi  = 5; allocate(geom.face(2).poi(5));  geom.face( 2).poi(1:5) = [  1,  3,  7, 10,  4 ]
    geom.face(3).n_poi  = 5; allocate(geom.face(3).poi(5));  geom.face( 3).poi(1:5) = [  1,  4,  9,  6,  2 ]
    geom.face(4).n_poi  = 5; allocate(geom.face(4).poi(5));  geom.face( 4).poi(1:5) = [  2,  6, 12, 11,  5 ]
    geom.face(5).n_poi  = 5; allocate(geom.face(5).poi(5));  geom.face( 5).poi(1:5) = [  3,  8, 14, 13,  7 ]
    geom.face(6).n_poi  = 5; allocate(geom.face(6).poi(5));  geom.face( 6).poi(1:5) = [  4, 10, 16, 15,  9 ]
    geom.face(7).n_poi  = 5; allocate(geom.face(7).poi(5));  geom.face( 7).poi(1:5) = [  5, 11, 17, 14,  8 ]
    geom.face(8).n_poi  = 5; allocate(geom.face(8).poi(5));  geom.face( 8).poi(1:5) = [  6,  9, 15, 18, 12 ]
    geom.face(9).n_poi  = 5; allocate(geom.face(9).poi(5));  geom.face( 9).poi(1:5) = [  7, 13, 19, 16, 10 ]
    geom.face(10).n_poi = 5; allocate(geom.face(10).poi(5)); geom.face(10).poi(1:5) = [ 11, 12, 18, 20, 17 ]
    geom.face(11).n_poi = 5; allocate(geom.face(11).poi(5)); geom.face(11).poi(1:5) = [ 13, 14, 17, 20, 19 ]
    geom.face(12).n_poi = 5; allocate(geom.face(12).poi(5)); geom.face(12).poi(1:5) = [ 15, 16, 19, 20, 18 ]
end subroutine Exam_Platonic_Dodecahedron

! ---------------------------------------------------------------------------------------

! Example of Icosahedron
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Exam_Platonic_Icosahedron(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "05_Icosahedron"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Icosahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7        ! [-1   ],   The number of nucleotides in poly T loop, -1 : depending on distance
        para_set_seq_scaf    = 1        ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from seq.txt, 2 - random
        para_set_start_scaf  = 1        ! [7217, 4141], Starting nucleotide position of scaffold strand
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.05d0     ! Cylindrical model
    prob.move_x   =-0.5d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 12
    geom.n_face = 20

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  0.00000d0,  0.00000d0,  9.51058d0 ]; geom.iniP( 2).pos(1:3) = [  8.50650d0,  0.00000d0,  4.25326d0 ]
    geom.iniP( 3).pos(1:3) = [  2.62866d0,  8.09018d0,  4.25326d0 ]; geom.iniP( 4).pos(1:3) = [ -6.88192d0,  5.00001d0,  4.25326d0 ]
    geom.iniP( 5).pos(1:3) = [ -6.88192d0, -5.00001d0,  4.25326d0 ]; geom.iniP( 6).pos(1:3) = [  2.62866d0, -8.09018d0,  4.25326d0 ]
    geom.iniP( 7).pos(1:3) = [  6.88192d0,  5.00001d0, -4.25326d0 ]; geom.iniP( 8).pos(1:3) = [  6.88192d0, -5.00001d0, -4.25326d0 ]
    geom.iniP( 9).pos(1:3) = [ -2.62866d0,  8.09018d0, -4.25326d0 ]; geom.iniP(10).pos(1:3) = [ -8.50650d0,  0.00000d0, -4.25326d0 ]
    geom.iniP(11).pos(1:3) = [ -2.62866d0, -8.09018d0, -4.25326d0 ]; geom.iniP(12).pos(1:3) = [  0.00000d0,  0.00000d0, -9.51058d0 ]

    ! Set face connnectivity
    geom.face(1).n_poi  = 3; allocate(geom.face(1).poi(3));  geom.face( 1).poi(1:3) = [  1,  2,  3 ]
    geom.face(2).n_poi  = 3; allocate(geom.face(2).poi(3));  geom.face( 2).poi(1:3) = [  1,  3,  4 ]
    geom.face(3).n_poi  = 3; allocate(geom.face(3).poi(3));  geom.face( 3).poi(1:3) = [  1,  4,  5 ]
    geom.face(4).n_poi  = 3; allocate(geom.face(4).poi(3));  geom.face( 4).poi(1:3) = [  1,  5,  6 ]
    geom.face(5).n_poi  = 3; allocate(geom.face(5).poi(3));  geom.face( 5).poi(1:3) = [  1,  6,  2 ]
    geom.face(6).n_poi  = 3; allocate(geom.face(6).poi(3));  geom.face( 6).poi(1:3) = [  2,  6,  8 ]
    geom.face(7).n_poi  = 3; allocate(geom.face(7).poi(3));  geom.face( 7).poi(1:3) = [  2,  8,  7 ]
    geom.face(8).n_poi  = 3; allocate(geom.face(8).poi(3));  geom.face( 8).poi(1:3) = [  2,  7,  3 ]
    geom.face(9).n_poi  = 3; allocate(geom.face(9).poi(3));  geom.face( 9).poi(1:3) = [  3,  7,  9 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  3,  9,  4 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  4,  9, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  4, 10,  5 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  5, 10, 11 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  5, 11,  6 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  6, 11,  8 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  7,  8, 12 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  7, 12,  9 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  8, 11, 12 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  9, 12, 10 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 10, 12, 11 ]
end subroutine Exam_Platonic_Icosahedron

! ---------------------------------------------------------------------------------------

end module Exam_Platonic