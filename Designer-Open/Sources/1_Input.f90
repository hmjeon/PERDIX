!
! ---------------------------------------------------------------------------------------
!
!                                   Module for Input
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/08/03
!                                                            Last  modified  : 2016/10/12
!
! ---------------------------------------------------------------------------------------
!
module Input

    use Importer

    use Exam_Platonic
    use Exam_Archi
    use Exam_Catalan
    use Exam_Johnson
    use Exam_Miscellaneous
    use Exam_Asymmetric
    use Exam_OpenGeo
    use Exam_Prism

    use Section

    use Para
    use Math
    use List

    implicit none

    public  Input_Initialization
    public  Input_Initialization_Autorun

    private Input_Read_Parameter
    private Input_Reset_Parameter
    private Input_Set_Command
    private Input_Print_Problem
    private Input_Print_Section
    private Input_Print_Num_BP_Edge
    private Input_Set_Problem
    private Input_Select_Problem
    private Input_Select_File
    private Input_Set_Section
    private Input_Find_Max_Min_Section
    private Input_Set_Section_Connectivity
    private Input_Set_Num_BP_Edge
    private Input_Convert_Face_To_Line
    private Input_Scale_Init_Geometry
    private Input_Set_Path
    private Input_Set_Workplace
    private Input_Write_GEO_File
    private Input_Chimera_Init_Geometry
    private Input_Tecplot_Init_Geometry
    private Input_Generate_Schlegel_Diagram

contains

! ---------------------------------------------------------------------------------------

! Initialize parameters and inputs
! Last updated on Friday 15 Apr 2016 by Hyungmin
subroutine Input_Initialization(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: arg, i, j, k, n_problem, n_section, n_edge_length
    character(10) :: char_problem, char_section, char_edge_length
    character(10) :: char_vertex_design, char_vertex_modify, char_cut_stap
    logical :: results

    ! Read parameters from env.dat
    call Input_Read_Parameter

    ! Set command environment
    call Input_Set_Command

    if(iargc() == 0) then

        ! ==================================================
        !
        ! Select geometry and cross-section from screen input
        !
        ! ==================================================
        ! Print pre-defined problems
        call Input_Print_Problem; read(*, *) prob.sel_prob
        if(prob.sel_prob <= 0) stop

        ! Clean the screen
        results = SYSTEMQQ("cls")

        ! Print pre-defined cross-sections
        call Input_Print_Section; read(*, *) prob.sel_sec
        if(prob.sel_sec <= 0 .or. prob.sel_sec >= 6) stop

        ! Print pre-defined edge length(bps)
        if(prob.sel_prob <= 58 .or. prob.sel_prob >= 62) then
            call Input_Print_Num_BP_Edge(prob)
            read(*, *) prob.sel_bp_edge
            if(prob.sel_bp_edge <= 0) stop
        end if
    else

        ! ==================================================
        !
        ! Input system by using main argument
        !
        ! ==================================================
        arg = 1; call getarg(arg, char_problem)         ! 1st argument, problem
        arg = 2; call getarg(arg, char_section)         ! 2nd argument, section
        arg = 3; call getarg(arg, char_edge_length)     ! 3rd argument, edge length
        arg = 4; call getarg(arg, char_vertex_design)   ! 4th argument, vertex design
        arg = 5; call getarg(arg, char_vertex_modify)   ! 5th argument, vertex modify
        arg = 6; call getarg(arg, char_cut_stap)        ! 6th argument, cutting method

        read(char_problem,     *), n_problem
        read(char_section,     *), n_section
        read(char_edge_length, *), n_edge_length

        ! Set parameters of problem, section, edge length, and so on
        prob.sel_prob        = n_problem
        prob.sel_sec         = n_section
        prob.sel_bp_edge     = n_edge_length
        para_vertex_design   = trim(char_vertex_design)
        para_vertex_modify   = trim(char_vertex_modify)
        para_cut_stap_method = trim(char_cut_stap)
        para_start_bp_ID     = -1
    end if

    ! ==================================================
    !
    ! Set problem, cross-section and the number of basepairs
    !
    ! ==================================================
    ! Set cross-section based on square or honeycomb lattices
    call Input_Set_Section(prob, geom)

    ! Set the number of basepair on edge that has minimum length
    call Input_Set_Num_BP_Edge(prob, geom)

    ! Set problem to be solved
    call Input_Set_Problem(prob, geom)

    ! External running
    if(iargc() /= 0 .and. para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "max"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        !para_const_edge_mesh = "off"    ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! ==================================================
    !
    ! Prepair geometry - line generation and scaling
    !
    ! ==================================================
    ! Convert surface to line connectivity
    call Input_Convert_Face_To_Line(geom)

    ! Set geometric scale with initial minimum length
    call Input_Scale_Init_Geometry(geom)

    ! ==================================================
    !
    ! Set environment and write initial geometry
    !
    ! ==================================================
    ! Set working and Chimera path
    call Input_Set_Path(prob)

    ! Remove previous working directory and make new one
    call Input_Set_Workplace(prob)

    ! Write *.geo file
    call Input_Write_GEO_File(prob, geom)

    ! Write initial geometry
    call Input_Chimera_Init_Geometry(prob, geom)

    ! Write initial geometry for Tecplot
    call Input_Tecplot_Init_Geometry(prob, geom)

    ! Generate Schlegel diagram
    call Input_Generate_Schlegel_Diagram(prob, geom)

    ! Open output progress file (unit 11 is used for global output file)
    open(unit=11, file=trim(prob.path_work1)//"DNAcs.txt", form="formatted")

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a )"), "   +--------------------------------------------------------------------+"
        write(i, "(a )"), "   |                                                                    |"
        write(i, "(a )"), "   |          1. Inputs - geometry, cross-section, edge length          |"
        write(i, "(a )"), "   |                                                                    |"
        write(i, "(a )"), "   +--------------------------------------------------------------------+"
        write(i, "(a )")
        call Space(i, 6)
        write(i, "(a)"), "1.1. Geometry"
        call Space(i, 11)
        write(i, "(a)"), "* Geometric name                    : "//trim(prob.name_prob)
        call Space(i, 11)
        write(i, "(a)"), "* Geometric file type               : "//trim(prob.type_file)
        call Space(i, 11)
        write(i, "(a)"), "* Surface type                      : "//trim(prob.type_geo)//" geometry"
        call Space(i, 11)
        write(i, "(a)"), "* The number of faces               : "//trim(adjustl(Int2Str(geom.n_face)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of points              : "//trim(adjustl(Int2Str(geom.n_iniP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of edges               : "//trim(adjustl(Int2Str(geom.n_iniL)))
        write(i, "(a)")

        call Space(i, 6)
        write(i, "(a)"), "1.2. Cross-section information"
        call Space(i, 11)
        write(i, "(a)"), "* Section type                      : "//trim(geom.sec.types)//" lattice"
        call Space(i, 11)
        write(i, "(a)"), "* The number of duplexes            : "//trim(adjustl(Int2Str(geom.n_sec)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of rows                : "//trim(adjustl(Int2Str(geom.sec.maxR-geom.sec.minR+1)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of columns             : "//trim(adjustl(Int2Str(geom.sec.maxC-geom.sec.minC+1)))
        call Space(i, 11)
        write(i, "(a)"), "* Reference row                     : "//trim(adjustl(Int2Str(geom.sec.ref_row)))
        call Space(i, 11)
        write(i, "(a)"), "* Reference min/max column          : "&
            //trim(adjustl(Int2Str(geom.sec.ref_minC)))//" / "//trim(adjustl(Int2Str(geom.sec.ref_maxC)))
        write(i, "(a)")

        call Space(i, 6)
        write(i, "(a)"), "1.3. Edge length"
        call Space(i, 11)
        write(i, "(a)"), "* The minimum edge length           : "//trim(adjustl(Int2Str(prob.n_bp_edge)))
        write(i, "(a)")
        call Space(i, 6)

        write(i, "(a)"), "1.4. Design parameters"
        call Space(i, 11)
        write(i, "(a)"), "* Junction modification             : "//trim(para_junc_ang)
        call Space(i, 11)
        write(i, "(a)"), "* Vertex design                     : "//trim(para_vertex_design)//" vertex"
        call Space(i, 11)
        write(i, "(a)"), "* Vertex modification               : "//trim(para_vertex_modify)
        call Space(i, 11)
        write(i, "(a)"), "* Sticky-end for self connection    : "//trim(para_sticky_self)
        call Space(i, 11)
        write(i, "(a)"), "* Unpaired scffold nucleotides      : "//trim(para_unpaired_scaf)
        call Space(i, 11)
        write(i, "(a)"), "* Unpaired nucleotides at vertex    : "//trim(para_unpaired_square)
        call Space(i, 11)
        write(i, "(a$)"), "* The number of bases in Tn         : "
        if(para_n_base_tn == -1) then
            write(i, "(a)"), "depending on distance"
        else
            write(i, "(a)"), trim(adjustl(Int2Str(para_n_base_tn)))
        end if
        call Space(i, 11)
        write(i, "(a)"), "* Distance btw phosphate groups     : "//trim(adjustl(Dble2Str(para_dist_pp)))
        call Space(i, 11)
        write(i, "(a)"), "* Starting base pair ID             : "//trim(adjustl(Int2Str(para_start_bp_ID)))
        call Space(i, 11)
        write(i, "(a)"), "* Axial rise distance [nm]          : "//trim(adjustl(Dble2Str(para_dist_bp)))
        call Space(i, 11)
        write(i, "(a)"), "* Radius of helix [nm]              : "//trim(adjustl(Dble2Str(para_rad_helix)))
        call Space(i, 11)
        write(i, "(a)"), "* The Gap between helixes           : "//trim(adjustl(Dble2Str(para_gap_helix)))
        call Space(i, 11)
        write(i, "(a)"), "* Angle of minor groove             : "//trim(adjustl(Dble2Str(para_ang_minor)))
        call Space(i, 11)
        write(i, "(a)"), "* Angle correction factor           : "//trim(adjustl(Dble2Str(para_ang_correct)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap between two scaf xovers       : "//trim(adjustl(Int2Str(para_gap_xover_two_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap between xover(stap) and bound : "//trim(adjustl(Int2Str(para_gap_xover_bound_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap between stap and scaf xovers  : "//trim(adjustl(Int2Str(para_gap_xover_two)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap between xover and first nick  : "//trim(adjustl(Int2Str(para_gap_xover_nick1)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap between xover and nick        : "//trim(adjustl(Int2Str(para_gap_xover_nick)))
        call Space(i, 11)
        write(i, "(a)"), "* Staple cutting method             : "//trim(para_cut_stap_method)
        call Space(i, 11)
        write(i, "(a)"), "* Non-circular stap by single xover : "//trim(para_set_stap_sxover)
        call Space(i, 11)
        write(i, "(a$)"), "* Minimum # of bases in scaf strand : "
        if(para_max_cut_scaf == 0) then
            write(i, "(a)"), "infinite"
        else
            write(i, "(a)"), trim(adjustl(Int2Str(para_max_cut_scaf)))
        end if
        call Space(i, 11)
        write(i, "(a)"), "* Minimum # of bases in stap strand : "//trim(adjustl(Int2Str(para_min_cut_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Midium # of bases in stap strand  : "//trim(adjustl(Int2Str(para_mid_cut_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Maximum # of bases in stap strand : "//trim(adjustl(Int2Str(para_max_cut_stap)))
        call Space(i, 11)
        write(i, "(a$)"), "* Scaffold sequence                 : "
        if(para_set_seq_scaf == 0) then
            write(i, "(a)"), "M13mp18 sequence"
        else if(para_set_seq_scaf == 1) then
            write(i, "(a)"), "user-defined sequence"
        else if(para_set_seq_scaf == 2) then
            write(i, "(a)"), "random sequence"
        end if
        write(i, "(a)")
    end do
end subroutine Input_Initialization

! ---------------------------------------------------------------------------------------

! Initialize all inputs for autorun
! Last updated on Friday 11 November 2016 by Hyungmin
subroutine Input_Initialization_Autorun(prob, geom, ii, sec, edge, vertex, char_cut, char_junc)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom
    integer,        intent(inout) :: sec
    integer,        intent(in)    :: ii, edge, vertex
    character(10),  intent(in)    :: char_cut, char_junc

    integer :: i

    ! Reset parameters
    call Input_Reset_Parameter
    para_external = "on"

    ! Reset data structures
    prob.n_cng_min_stap = 0
    prob.n_cng_max_stap = 0

    ! Set command environment
    call Input_Set_Command

    ! Set staple cutting and junction modification
    para_cut_stap_method = char_cut
    para_junc_ang        = char_junc

    ! Set parameters of problem
    prob.sel_prob = ii

    ! Set vertex design
    if(vertex == 1) para_vertex_design = "flat"
    if(vertex == 2) para_vertex_design = "beveled"

    ! Set junction modification
    para_vertex_modify = "const"

    ! For honeycomb section
    if(sec == 2 .or. sec == 3 .or. sec == 4) then

        ! Tetrahedron, cube, octahedron, triangular bipyramid, double helix, nested octahedron, torus
        if(ii == 1 .or. ii == 2 .or. ii == 3 .or. ii == 17 .or. ii == 18 .or. ii == 41 .or. ii == 43 .or. ii == 44 .or. ii ==59 .or. ii == 60 .or. ii == 61 .or. ii == 64) then
            sec = 2
            para_vertex_modify = "const"
        !else if(ii == 39) then
        !    sec = 3
        !    para_vertex_modify = "mod"
        else
            if(sec == 3) sec = 3
            if(sec == 4) sec = 4
            para_vertex_modify = "const"
        end if
    end if

    prob.sel_sec     = sec
    prob.sel_bp_edge = edge
    para_start_bp_ID = -1

    ! ==================================================
    !
    ! Set problem, cross-section and the number of basepairs
    !
    ! ==================================================
    ! Set cross-section based on square or honeycomb lattices
    call Input_Set_Section(prob, geom)

    ! Set the number of basepair on edge that has minimum length
    call Input_Set_Num_BP_Edge(prob, geom)

    ! Set problem to be solved
    call Input_Set_Problem(prob, geom)

    ! ==================================================
    !
    ! Prepair geometry - line generation and scaling
    !
    ! ==================================================
    ! Convert surface to line connectivity
    call Input_Convert_Face_To_Line(geom)

    ! Set geometric scale with initial minimum length
    call Input_Scale_Init_Geometry(geom)

    ! ==================================================
    !
    ! Set environment and write initial geometry
    !
    ! ==================================================
    ! Set working and Chimera path
    call Input_Set_Path(prob)

    ! Remove previous working directory and make new one
    call Input_Set_Workplace(prob)

    ! Write *.geo file
    call Input_Write_GEO_File(prob, geom)

    ! Write initial geometry
    call Input_Chimera_Init_Geometry(prob, geom)

    ! Generate Schlegel diagram
    call Input_Generate_Schlegel_Diagram(prob, geom)

    ! Open output progress file (unit 11 is used for global output file)
    open(unit=11, file=trim(prob.path_work1)//"DNAcs.txt", form="formatted")

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a )"), "   ----------------------------------------------------------------------"
        write(i, "(a )"), "   |                                                                    |"
        write(i, "(a )"), "   |              1. Inputs of Geometry and Cross-Section               |"
        write(i, "(a )"), "   |                                                                    |"
        write(i, "(a )"), "   ----------------------------------------------------------------------"
        write(i, "(a )")
        call Space(i, 6)
        write(i, "(a )"), "1.1. Geometric data"
        call Space(i, 11)
        write(i, "(a )"), "* Geometric name             :   "//trim(prob.name_file)
        call Space(i, 11)
        write(i, "(a )"), "* Geometric type             :   "//trim(prob.type_file)
        call Space(i, 11)
        write(i, "(a$)"), "* The number of BPs on edge  :  "
        write(i, "(i7)"), prob.n_bp_edge
        write(i, "(a )")

        call Space(i, 6)
        write(i, "(a )"), "1.2. Geometric information"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of points       : "
        write(i, "(i7)"), geom.n_iniP
        call Space(i, 11)
        write(i, "(a$)"), "* The number of edges        : "
        write(i, "(i7)"), geom.n_iniL
        call Space(i, 11)
        write(i, "(a$)"), "* The number of faces        : "
        write(i, "(i7)"), geom.n_face
        write(i, "(a )")

        call Space(i, 6)
        write(i, "(a )"), "1.3. Cross-section information"
        call Space(i, 11)
        write(i, "(a )"), "* Section type               :   "//trim(geom.sec.types)
        call Space(i, 11)
        write(i, "(a$)"), "* The number of helixes      : "
        write(i, "(i7)"), geom.n_sec
        call Space(i, 11)
        write(i, "(a$)"), "* The number of rows         : "
        write(i, "(i7)"), geom.sec.maxR
        call Space(i, 11)
        write(i, "(a$)"), "* The number of columns      : "
        write(i, "(i7)"), geom.sec.maxC
        write(i, "(a )")

        call Space(i, 6)
        write(i, "(a )"), "1.4. Design type"
        call Space(i, 11)
        write(i, "(a )"), "* Geometry type              :   "//trim(prob.type_geo)
        call Space(i, 11)
        write(i, "(a )"), "* Vertex design              :   "//trim(para_vertex_design)//" vertex"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of bases in Tn  : "
        write(i, "(i7)"), para_n_base_tn
        call Space(i, 11)
        write(i, "(a$)"), "* No crossover region (stap) : "
        write(i, "(i7)"), para_gap_xover_bound_stap
        call Space(i, 11)
        write(i, "(a$)"), "* Starting base pair ID      : "
        write(i, "(i7)"), para_start_bp_ID
        write(i, "(a )")
    end do
end subroutine Input_Initialization_Autorun

! ---------------------------------------------------------------------------------------

! Read parameters from external txt file, env.txt
! Last updated on Thursday 14 June 2016 by Hyungmin
subroutine Input_Read_Parameter
    character(200) :: ctemp
    integer :: i

    ! Open file
    open(unit=1, file="env.txt", form="formatted")

    ! External parameter loading, this should be always ".true."
    read(1, *), ctemp, para_external

    ! If the external mode is on
    if(para_external == "on") then

        ! Program parameters
        read(1, *), ctemp, para_preset
        read(1, *), ctemp, para_output_Tecplot
        read(1, *), ctemp, para_cmd_Tecplot
        read(1, *), ctemp, para_cmd_Chimera
        read(1, *), ctemp, para_fig_output
        read(1, *), ctemp, para_fig_route_step
        read(1, *), ctemp, para_fig_bgcolor
        read(1, *), ctemp, para_fig_view
        read(1, *), ctemp, para_n_route_step
        read(1, *), ctemp, para_type_cndo
        read(1, *), ctemp
        read(1, "(a100)"), para_path_Chimera

        ! Parameters for junction modification
        read(1, *), ctemp, para_junc_ang
        read(1, *), ctemp, para_const_edge_mesh
        read(1, *), ctemp, para_sticky_self
        read(1, *), ctemp, para_unpaired_scaf
        read(1, *), ctemp, para_unpaired_square
        read(1, *), ctemp, para_vertex_modify
        read(1, *), ctemp, para_vertex_design

        ! Paramters for B-from DNA generation
        read(1, *), ctemp, para_dist_pp
        read(1, *), ctemp, para_dist_bp
        read(1, *), ctemp, para_rad_helix
        read(1, *), ctemp, para_gap_helix
        read(1, *), ctemp, para_ang_minor
        read(1, *), ctemp, para_ang_correct
        read(1, *), ctemp, para_n_base_tn
        read(1, *), ctemp, para_start_bp_ID

        ! Paramters for scaffold route
        read(1, *), ctemp, para_weight_edge
        read(1, *), ctemp, para_method_MST
        read(1, *), ctemp, para_method_sort
        read(1, *), ctemp, para_adjacent_list
        read(1, *), ctemp, para_all_spanning

        ! Parameter for sequence design
        read(1, *), ctemp, para_cut_stap_method
        read(1, *), ctemp, para_set_stap_sxover
        read(1, *), ctemp, para_output_design
        read(1, *), ctemp, para_set_xover_scaf
        read(1, *), ctemp, para_gap_xover_two_scaf
        read(1, *), ctemp, para_gap_xover_bound_scaf
        read(1, *), ctemp, para_gap_xover_bound_stap
        read(1, *), ctemp, para_gap_xover_two
        read(1, *), ctemp, para_gap_xover_nick1
        read(1, *), ctemp, para_gap_xover_nick
        read(1, *), ctemp, para_max_cut_scaf
        read(1, *), ctemp, para_min_cut_stap
        read(1, *), ctemp, para_mid_cut_stap
        read(1, *), ctemp, para_max_cut_stap
        read(1, *), ctemp, para_set_seq_scaf
        read(1, *), ctemp, para_set_start_scaf

        ! UCSF Chimera output control
        read(1, *), ctemp, para_write_101
        read(1, *), ctemp, para_write_102
        read(1, *), ctemp, para_write_103
        read(1, *), ctemp, para_write_104
        read(1, *), ctemp, para_write_301
        read(1, *), ctemp, para_write_302
        read(1, *), ctemp, para_write_303
        read(1, *), ctemp, para_write_401
        read(1, *), ctemp, para_write_501
        read(1, *), ctemp, para_write_502
        read(1, *), ctemp, para_write_503
        read(1, *), ctemp, para_write_504
        read(1, *), ctemp, para_write_505
        read(1, *), ctemp, para_write_601_1
        read(1, *), ctemp, para_write_601_2
        read(1, *), ctemp, para_write_601_3
        read(1, *), ctemp, para_write_601_4
        read(1, *), ctemp, para_write_601_5
        read(1, *), ctemp, para_write_606
        read(1, *), ctemp, para_write_607
        read(1, *), ctemp, para_write_608
        read(1, *), ctemp, para_write_609
        read(1, *), ctemp, para_write_610
        read(1, *), ctemp, para_write_701
        read(1, *), ctemp, para_write_711
        read(1, *), ctemp, para_write_702
        read(1, *), ctemp, para_write_703
        read(1, *), ctemp, para_write_705
        read(1, *), ctemp, para_write_706
        read(1, *), ctemp, para_write_710
        read(1, *), ctemp, para_write_801
        read(1, *), ctemp, para_write_802
        read(1, *), ctemp, para_write_803
        read(1, *), ctemp, para_write_804
        read(1, *), ctemp, para_write_805
        read(1, *), ctemp, para_write_808

        read(1, *), ctemp, para_chimera_axis
        read(1, *), ctemp, para_chimera_102_info
        read(1, *), ctemp, para_chimera_301_info
        read(1, *), ctemp, para_chimera_302_info
        read(1, *), ctemp, para_chimera_303_info
        read(1, *), ctemp, para_chimera_401_info
        read(1, *), ctemp, para_chimera_502_ori
        read(1, *), ctemp, para_chimera_503_mod
        read(1, *), ctemp, para_chimera_504_info
        read(1, *), ctemp, para_chimera_601_dir
        read(1, *), ctemp, para_chimera_609_cyl
        read(1, *), ctemp, para_chimera_609_dir
    end if

    ! For auto generation of output figures
    if(para_fig_output == "on" .or. para_cmd_Chimera == "on") then
        para_fig_bgcolor      = "white"     ! [black, white, all], Background colo
        para_chimera_axis     = .false.     ! Plot with axis at the ceneter of geometry (*.bild)
        para_chimera_102_info = .false.     ! Plot with edge and point number (_init_geo.bild)
        para_chimera_301_info = .false.     ! Plot with edge and point number (_check_geo.bild)
        para_chimera_302_info = .false.     ! Plot with edge and point number (_init_geo_local.bild)
        para_chimera_303_info = .false.     ! Plot with edge and point number (_mod_geo.bild)
        para_chimera_401_info = .false.     ! Plot with edge and point number (_cross_geo.bild)
        para_chimera_502_ori  = .false.     ! Plot with helix z-direction (_line.bild / _node.bild)
        para_chimera_503_mod  = .false.     ! Plot with modified edges (_mesh.bild)
        para_chimera_504_info = .false.     ! Plot with edge and point number (_cross_geo_mod.bild)
        para_chimera_601_dir  = .false.     ! Plot with strand direction (_scaf.bild / _stap.bild)
        para_chimera_609_cyl  = .false.     ! Plot with cylinderical representation (_atom.bild)
        para_chimera_609_dir  = .false.     ! Plot with strand direction (_atom.bild)

        ! UCSF Chimera output control
        if(para_fig_output == "on") then
            para_write_101   = .false.      ! GEO file,                                Input_Write_GEO_File,             ".geo"
            para_write_102   = .false.      ! Initial geometry,                        Input_Chimera_Init_Geometry,      "init_geo.bild"
            para_write_103   = .false.      ! Faced initial geometry,                  Input_Tecplot_Init_Geometry,      "init_geo_face.dat"
            para_write_104   = .false.      ! Schlegel diagram,                        Input_Chimera_Schlegel_Diagram,   "_schlegel.bild"
            para_write_301   = .false.      ! Initial geometry with face orientation,  ModGeo_Chimera_Check_Geometry,    "_check_geo.bild"
            para_write_302   = .false.      ! Initial geometry with local vector,      ModGeo_Chimera_Init_Geometry_L,   "_init_geo_local.bild"
            para_write_303   = .false.      ! Modified geometry seperated from vertex, ModGeo_Chimera_Mod_Geometry,      "_mod_geo.bild"
            para_write_401   = .false.      ! Cross-sectional geometry,                Section_Chimera_Cross_Geometry,   "_cross_geo.bild"
            para_write_501   = .false.      ! Cylindrical model with orientation,      Basepair_Chimera_Cylinder_Ori,    "_cyl_ori1.bild"
            para_write_502   = .true.       ! Cylindrical model,                       Basepair_Chimera_Cylinder,        "_cyl1.bild", "_cyl2.bild"
            para_write_503   = .false.      ! Basepair model,                          Basepair_Chimera_Mesh,            "_mesh.bild"
            para_write_504   = .false.      ! Cross-sectional geometry,                Basepair_Chimera_Cross_Geometry,  "_cross_geo_mod.bild"
            para_write_505   = .true.       ! Txt file on edge length,                 Basepair_Write_Edge_Length,       "_edge_length_1.txt", "_edge_length_2.txt"
            para_write_601_1 = .false.      ! Route 1, seperated edges,                Route_Chimera_Route, step 1,      "_route1_scaf.bild", "_route1_stap.bild"
            para_write_601_2 = .false.      ! Route 2, contruction closed loop,        Route_Chimera_Route, step 2,      "_route2_scaf.bild", "_route2_stap.bild"
            para_write_601_3 = .false.      ! Route 3, centered crossovers             Route_Chimera_Route, step 3,      "_route3_scaf.bild", "_route3_stap.bild"
            para_write_601_4 = .false.      ! Route 4, modified centered crossovers,   Route_Chimera_Route, step 4,      "_route4_scaf.bild", "_route4_stap.bild"
            para_write_601_5 = .false.      ! Route 5, scaffold route,                 Route_Chimera_Route, step 5,      "_route5_scaf.bild", "_route5_stap.bild"
            para_write_606   = .false.      ! Sapnning tree for dual-graph,            Route_Graph_Chimera_Spanning_Tre, "_spantree.bild"
            para_write_607   = .false.      ! Crossovers based on basepair model,      Route_Chimera_Crossovers,         "_crossovers.bild"
            para_write_608   = .false.      ! 3-orientation vectors,                   Route_Chimera_Orientation,        "_orientation.bild"
            para_write_609   = .false.      ! Atomic model without sequence design,    Route_Chimera_Atom,               "_atom.bild"
            para_write_610   = .false.      ! Possible centered scaffold crossovers,   Route_Write_Centered_Scaf_Xover,  "_scaf_xover.txt"
            para_write_701   = .false.      ! Txt on sequence design data,             SeqDesign_Write_Strand,           "strand.txt"
            para_write_711   = .false.      ! Csv file for sequence data,              SeqDesign_Write_Strand,           "sequence.csv"
            para_write_702   = .false.      ! Atomic model with sequence design,       SeqDesign_Chimera_Atom,           "_atom_nick.bild"
            para_write_703   = .false.      ! Route 6, strand route with nick,         SeqDesign_Chimera_Route,          "_route6_scaf.bild", "_route6_stap.bild"
            para_write_705   = .true.       ! Sequence model,                          SeqDesign_Chimera_Sequence,       "_sequence_design.bild"
            para_write_706   = .false.      ! Atomic model bases on strands/sequence,  SeqDesign_Chimera_Strand,         "_strand.bild", "_sequence.bild"
            para_write_710   = .false.      ! Edge-based sequence design,              SeqDesign_Write_Graphical_Output, "_design_edgeX"
            para_write_801   = .false.      ! Txt on basepair based data,              Output_Write_Basepair,            "_basepair.txt"
            para_write_802   = .false.      ! Txt on nucleotide based data,            Output_Write_Base,                "_base.txt"
            para_write_803   = .false.      ! CanDo input file,                        Output_Write_CanDo,               ".cndo"
            para_write_804   = .false.      ! Tecplot input file,                      Output_Write_TecPlot,             "_tecplot.dat"
            para_write_805   = .false.      ! ADINA input file,                        Output_Write_ADINA,               "_adina.in"
            para_write_808   = .false.      ! Txt on sectional edges based sequence,   Output_Write_Sequence_CroL,       "_seq_line.txt"
        end if
    end if

    close(unit=1)
end subroutine Input_Read_Parameter

! ---------------------------------------------------------------------------------------

! Reset paramters as default values
! Last updated on Thursday 15 September 2016 by Hyungmin
subroutine Input_Reset_Parameter

    ! Program parameters
    para_output_Tecplot  = "off"      ! [off, on], Output files for Tecplot(http://www.tecplot.com/) to draw vector image
    para_cmd_Tecplot     = "off"      ! [off, on], Command file to run TecPlot automatically
    para_cmd_Chimera     = "off"      ! [off, on], Command file to run UCSF Chimera(https://www.cgl.ucsf.edu/chimera/) automatically
    para_fig_output      = "off"      ! [off, on], Automatic figure generation from outputs
    para_fig_route_step  = "off"      ! [off, on], Automatic figure generation from route steps
    para_fig_bgcolor     = "black"    ! [black, white, all], Background color for figures from UCSF Chimera
    para_fig_view        = "preset"   ! [preset, xy, xz, xyz, all], Viewpoint for figures from UCSF Chimera
    para_n_route_step    = 5          ! [5], The number of steps in routing progress
    para_type_cndo       = 2          ! [2, 1], CanDo file option, 1 : original format, 2 : updated format
    para_path_Chimera    = "C:\Program Files\Chimera 1.10.2\bin\chimera.exe"

    ! Parameters for junction modification
    para_junc_ang        = "opt"      ! [opt, max, ave, min], Junction gap modification for different arm angle
    para_const_edge_mesh = "off"      ! [off, on], Constant edge length from polyhedra mesh
    para_sticky_self     = "off"      ! [off, on], Sticky-end for self connection on henycomb cross-section
    para_unpaired_scaf   = "on"       ! [on, off], Unpaired scffold nucleotides
    para_unpaired_square = "on"       ! [on, off], Unpaired scffold and staple nucleotides at the vertex
    para_vertex_modify   = "const"    ! [mod, const], Vertex modification to avoid clash
    para_vertex_design   = "flat"     ! [flat, beveled], Vertex design

    ! Paramters for B-from DNA generation
    para_dist_pp         = 0.42d0     ! [0.42, 0.6], distance between adjacent phosphate groups, nm
    para_dist_bp         = 0.34d0     ! [0.34 ], Axial rise distance, nm
    para_rad_helix       = 1.0d0      ! [1.0  ], The radius of the DNA helix, nm
    para_gap_helix       = 0.25d0     ! [0.25 ], The Gap between two helixes, nm
    para_ang_minor       = 150.0d0    ! [150.0], An angle of minor groove, degree
    para_ang_correct     = 0.0d0      ! [0.0  ], Correction factor to adjust orientation, degree
    para_n_base_tn       = -1         ! [-1   ], The number of nucleotides in poly Tn loop, -1 : depending on distance
    para_start_bp_ID     = -1         ! [-1   ], Starting base pair ID for the reference, -1 : pre-defined starting BP

    ! Paramters for scaffold route
    para_weight_edge     = "on"       ! [on, off], Assign weight factor into edges of dual graph
    para_method_MST      = "prim"     ! [prim, kruskal, greedy], Minimum spanning tree algorithm
    para_method_sort     = "quick"    ! [none, quick, shell], Sorting algorithm to find MST for Prim or Kruskal
    para_adjacent_list   = "off"      ! [off, on], Output for adjacent list for Prim or Kruskal
    para_all_spanning    = "on"       ! [on, off], Possible all spanning tree when # of edges is less than 12 for Prim or Kruskal

    para_cut_stap_method = "max"      ! [max, mix, opt, min, mid], Cutting method to make short staple strand, opt - 14nt seeds
    para_set_stap_sxover = "off"      ! [off, on], To make non-circular staple by single crossover (when para_set_stap_sxover is "on")
    para_output_design   = "arrow"    ! [arrow, seq, strand], Graphical output type for sequence design
    para_set_xover_scaf  = "split"    ! [split, center], Setting possible scaffold strand

    para_gap_xover_two_scaf   = 3     ! [3 ], The minimum gap between two scaffold crossovers
    para_gap_xover_bound_scaf = 7     ! [7 ], The mimimum gap between scaffold crossover and vertex boundary
    para_gap_xover_bound_stap = 6     ! [6 ], The mimimum gap between staple crossover and vertex boundary
    para_gap_xover_two        = 6     ! [6 ], The minimum gap between scaffold and staple crossovers
    para_gap_xover_nick1      = 10    ! [10], The minimum gap between xover(scaf/stap)/Tn and first nick
    para_gap_xover_nick       = 3     ! [3 ], The minimum gap between xover and nick, if staple length exceeds 60, redesign with num - 1

    para_max_cut_scaf         = 0     ! [0, -1], The maximum number of nucleotides for one scaffold strand (maximum : 7249nt)
    para_min_cut_stap         = 20    ! [20], The minimum number of nucleotides for one staple strand
    para_mid_cut_stap         = 40    ! [40], The optimal number of nucleotides for one staple strand
    para_max_cut_stap         = 60    ! [60], The maximum number of nucleotides for one staple strand
    para_set_seq_scaf         = 0     ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from seq.txt, 2 - random
    para_set_start_scaf       = 1     ! [1], Starting nucleotide position of scaffold strand
end subroutine Input_Reset_Parameter

! ---------------------------------------------------------------------------------------

! Set command environment
! Last updated on Saturday 20 Feb 2016 by Hyungmin
subroutine Input_Set_Command
    logical :: results

    ! Set command environments
    results = SYSTEMQQ('title DNAcs')                       ! cmd title
    results = SYSTEMQQ('mode con: cols=135 lines=6000')     ! cmd size
    results = SYSTEMQQ('color')                             ! convert color, 02, f0, f1, f2
    results = SYSTEMQQ('date /t')                           ! display time
    !results = SYSTEMQQ('hostname')                          ! display hostname of the computer
    !results = SYSTEMQQ('ver')                               ! display version information
end subroutine Input_Set_Command

! ---------------------------------------------------------------------------------------

! Print pre-defined problems
! Last updated on Wednesday 24 Feb 2016 by Hyungmin
subroutine Input_Print_Problem
    write(0, "(a)")
    write(0, "(a)"), "       +=====================================================================================+"
    write(0, "(a)"), "       |                                                                                     |"
    write(0, "(a)"), "       |              DNAcs - DNA Nanostructure with Arbitrary Cross-Section                 |"
    write(0, "(a)"), "       |                               Programmed by Hyungmin Jun (hyeongminjeon@gmail.com)  |"
    write(0, "(a)"), "       |                                                                                     |"
    write(0, "(a)"), "       +=====================================================================================+"
    write(0, "(a)")
    write(0, "(a)"), "   A. First input - Geometry discretized by surface mesh [# of different edge length, problem size nomalized by Tetrahedron]"
    write(0, "(a)")
    write(0, "(a)"), "      I - Pre-defined geometries: Platonic solids"
    write(0, "(a)"), "          ---------------------------------------"
    write(0, "(a)"), "        *1. Tetrahedron[1, 1],   *2. *Cube[1, 2],   *3. Octahedron[1, 2],   4. Dodecahedron[1, 5],   5. Icosahedron[1, 5]"
    write(0, "(a)")
    write(0, "(a)"), "     II - Pre-defined geometries: Archimedean solids"
    write(0, "(a)"), "          ------------------------------------------"
    write(0, "(a)"), "         6. Cubeoctahedron[1, 5],             7. Icosidodecahedron[1, 12],        8. Rhombicuboctahedron[1, 10]"
    write(0, "(a)"), "         9. Snub Cube[1, 11],                10. Truncated Cube[1, 8],           11. Truncated Cuboctahedron[1, 14]"
    write(0, "(a)"), "        12. Truncated Dodecahedron[1, 21],   13. Truncated Icosahedron[1, 16],   14. Truncated Octahedron[1, 7]"
    write(0, "(a)"), "        15. Truncated Tetrahedron[1, 4]"
    write(0, "(a)")
    write(0, "(a)"), "    III - Pre-defined geometries: Johnson solids"
    write(0, "(a)"), "          --------------------------------------"
    write(0, "(a)"), "        16. Gyroelongated Pentagonal Pyramid (J11)[1, 4],     *17. Triangular Bipyramid (J12)[2, 2]"
    write(0, "(a)"), "       *18. Pentagonal Bipyramid (J13)[2, 3],                  19. Gyroelongated Square Bipyramid (J17)[2, 4]"
    write(0, "(a)"), "        20. Square Gyrobicupola (J29)[2, 6],                   21. Pentagonal Orthocupolarotunda (J32)[2, 10]"
    write(0, "(a)"), "        22. Pentagonal Orthobirotunda (J34)[1, 12],            23. Elongated Pentagonal Gyrobicupola (J39)[1, 13]"
    write(0, "(a)"), "        24  Elongated Pentagonal Gyrobirotunda (J43)[1, 17],   25. Gyroelongated Square Bicupola (J45)[1, 10]"
    write(0, "(a)")
    write(0, "(a)"), "     IV - Pre-defined geometries: Catalan solids"
    write(0, "(a)"), "          --------------------------------------"
    write(0, "(a)"), "        26. Rhombic Dodecahedron[1, 4],           27. Rhombic Triacontahedron[1, 10],   28. Deltoidal Icositetrahedron[2, 9]"
    write(0, "(a)"), "        29. Pentagonal Icositetrahedron[2, 12],   30. Triakis Octahedron[2, 7],         31. Disdyakis Dodecahedron[3, 16]"
    write(0, "(a)"), "        32. Triakis Icosahedron[2, 18],           33. Pentakis Dodecahedron[2, 16],     34. Tetrakis Hexahedron[2, 7]"
    write(0, "(a)"), "        35. Triakis Tetrahedron[2, 4]"
    write(0, "(a)")
    write(0, "(a)"), "      V - Pre-defined geometries: Miscellaneous polyhedra"
    write(0, "(a)"), "          -----------------------------------------------"
    write(0, "(a)"), "        36. Heptagonal Bipyramid[2, 7],       37. Enneagonal Trapezohedron[2, 26],   38. Small Stell Dodecahedron[2, 26]"
    write(0, "(a)"), "       #39. Rhombic Hexecontahedron[2, 23],   40. Goldberg dk5dgD[2, 38],           *41. Double Helix[15, 54]"
    write(0, "(a)"), "        42. Nested Cube[3, 12],              *43. Nested Octahedron[3, 16],         *44. Torus[10, 90]"
    write(0, "(a)"), "        45. Double Torus[X, 32]"
    write(0, "(a)")
    write(0, "(a)"), "       VI - Pre-defined geometries: Asymmetric and non-convex objects"
    write(0, "(a)"), "            ---------------------------------------------------------"
    write(0, "(a)"), "        46. Ball[7, 31],      47. Nickedtorus[33, 144],   48. Helix[23, 108],   49. Rod[6, 91],   50. Stickman[35, 91]"
    write(0, "(a)"), "        51. Bottle[51, 95],   52. Bunny[212, 218]"
    write(0, "(a)")
    write(0, "(a)"), "      VII - Pre-defined geometries: Ect"
    write(0, "(a)"), "            ---------------------------"
    write(0, "(a)"), "        53. Truncated Icosidodecahedron[1, 30],  54. Rhombicosidodecahedron[1, 20],     55. Snub Dodecahedron[1, 25]"
    write(0, "(a)"), "        56. Disdyakis Triacontahedron[3, 43],    57. Deltoidal Hexecontahedron[2, 25],  58. Pentagonal Hexecontahedron[2, 32]"
    write(0, "(a)"), "       *59. Asym Tetra I[63-63-63-73-52-42],    *60. Asym Tetra II[63-63-63-63-52-42], *61. Asym Tetra III[63-63-63-63-63-42]"
    write(0, "(a)"), "        62. Pentagonal Pyramid(J2)               63. Elongated Square Bipyramid(J15)"
    write(0, "(a)")
    write(0, "(a)"), "     VIII - Pre-defined geometries: Prism"
    write(0, "(a)"), "            ---------------------------"
    write(0, "(a)"), "        64. Triangular Prism,                    65. Pentagonal Prism,                  66. Hexagonal Prism"
    write(0, "(a)")
    !write(0, "(a)"), "      VII - Pre-defined geometries: Open surface mesh - under construction"
    !write(0, "(a)"), "        53. Plate Uniform Quad,     54. Plate Distorted Quad,           55. Plate Uniform Tri,     56. Plate Distorted Tri"
    !write(0, "(a)"), "        57. Circular Plate Quad,    58. Circular Plate Tri,             59. Annular Plate Quad"
    !write(0, "(a)"), "        60. Annular Plate Tri,      61. Hyperbolic Paraboloid Quad,     62. Hyperbolic Paraboloid Tri"
    !write(0, "(a)")
    write(0, "(a)"), "       ======================================================================================= "
    write(0, "(a)")
    write(0, "(a)"), "       100. Input from file (*.ply, *.geo, *.stl, *.wrl)"
    write(0, "(a)")
    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Problem

! ---------------------------------------------------------------------------------------

! Print pre-defined cross-sections
! Last updated on Monday 12 December 2016 by Hyungmin
subroutine Input_Print_Section
    write(0, "(a)")
    write(0, "(a)"), "   B. Second input - Pre-defined cross-sections"
    write(0, "(a)")
    write(0, "(a)"), "    - : crossover, = reference axis"
    write(0, "(a)")
    write(0, "(a)"), "    [Section defined on honeycomb lattice]"
    write(0, "(a)"), "     ================="
    write(0, "(a)")
    write(0, "(a)")
    write(0, "(a)"), "             [sec ID]                  [sec ID] "
    write(0, "(a)"), "      1.                   2.  @ @      4 3     "
    write(0, "(a)"), "        =@ @=  =0 1=          @   @    5   2    "
    write(0, "(a)"), "                              =@ @=    =0 1=    "
    write(0, "(a)"), "                                                "
    write(0, "(a)"), "                               [type 1, CW]     "
    write(0, "(a)"), "         [1 by 2]              [1 honeycomb]    "
    write(0, "(a)")
    write(0, "(a)")
    write(0, "(a)"), "                [sec ID]                [sec ID]                [sec ID]"
    write(0, "(a)"), "      3.  @-@      5 4        4.  @-@      5 4        5.  @-@      4 3  "
    write(0, "(a)"), "        =@   @=  =0   3=        =@   @=  =0   3=        =@   @=  =5   2="
    write(0, "(a)"), "          @-@      1 2            @-@      1 2            @-@      0 1  "
    write(0, "(a)"), "                                                                        "
    write(0, "(a)"), "         [type 2, CW]            [type 3, CCW]           [type 3, CW]   "
    write(0, "(a)"), "         [1 honeycomb]           [1 honeycomb]           [1 honeycomb]  "
    write(0, "(a)")
    write(0, "(a)")
    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Section

! ---------------------------------------------------------------------------------------

! Print pre-defined cross-sections
! Last updated on Thuesday 22 Mar 2016 by Hyungmin
subroutine Input_Print_Section_Old
    write(0, "(a)")
    write(0, "(a)"), "   B. Second input - Pre-defined cross-sections"
    write(0, "(a)")
    write(0, "(a)"), "    - : crossover, = reference axis"
    write(0, "(a)")
    write(0, "(a)"), "    [Section defined on square lattice]"
    write(0, "(a)"), "     =============="
    write(0, "(a)")
    write(0, "(a)"), "      1.   @-@     2.   @ @     3.   @-@     4. =@ @-@ @=   5.  @-@ @-@    6.  @ @-@ @ "
    write(0, "(a)"), "          =@ @=         @-@          @ @                       =@ @-@ @=       @-@ @-@ "
    write(0, "(a)"), "                       =@ @=         @-@                                      =@ @-@ @="
    write(0, "(a)"), "                                    =@ @=                                              "
    write(0, "(a)"), "        [2 by 2]     [3 by 2]     [4 by 2]      [1 by 4]       [2 by 4]       [3 by 4] "
    write(0, "(a)")
    write(0, "(a)"), "      7.  @-@ @-@      8. =@ @-@ @-@ @=     9.  @-@ @-@ @-@      10.  @ @-@ @-@ @ "
    write(0, "(a)"), "          @ @-@ @                              =@ @-@ @-@ @=          @-@ @-@ @-@ "
    write(0, "(a)"), "          @-@ @-@                                                    =@ @-@ @-@ @="
    write(0, "(a)"), "         =@ @-@ @=                                                                "
    write(0, "(a)"), "          [4 by 4]          [1 by 6]             [2 by 6]               [3 by 6]  "
    write(0, "(a)")
    write(0, "(a)"), "     11. @-@ @-@ @-@   12.  @ @-@ @   13.  @-@ @-@   14.    @ @     15.    @ @    "
    write(0, "(a)"), "         @ @-@ @-@ @        @     @        @     @       =@ @-@ @=       @-@ @-@  "
    write(0, "(a)"), "         @-@ @-@ @-@       =@ @-@ @=       @     @          @ @         =@ @-@ @= "
    write(0, "(a)"), "        =@ @-@ @-@ @=                     =@ @-@ @=                        @ @    "
    write(0, "(a)"), "          [4 by 6]        [3-4 hole]     [4-4 hole]     [3-4 cross]    [4-4 cross]"
    write(0, "(a)")
    write(0, "(a)"), "    [Section defined on honeycomb lattice]"
    write(0, "(a)"), "     ================="
    write(0, "(a)"), "                           [sec ID]             [sec ID]                                                        "
    write(0, "(a)"), "     16.       17.  @-@       5 4     18.  @ @     4 3    19.  @-@   @-@     20.     @-@        21.  @-@   @-@  "
    write(0, "(a)"), "        =@ @=     =@   @=   =0   3=       @   @   5   2      =@   @-@   @=        @-@   @-@         @   @-@   @ "
    write(0, "(a)"), "                    @-@       1 2         =@ @=   =0 1=        @-@   @-@        =@   @-@   @=        @-@   @-@  "
    write(0, "(a)"), "                                                                                  @-@   @-@        =@   @-@   @="
    write(0, "(a)"), "                  [type 1, const]         [type 2, mod]                              @-@             @-@   @-@  "
    write(0, "(a)"), "      [1 by 2]     [1 honeycomb]          [1 honeycomb]      [2 honeycomb]      [4 honeycomb]      [5 honeycomb]"
    write(0, "(a)")
    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Section_Old

! ---------------------------------------------------------------------------------------

! Print pre-defined the number of base pair on edges
! Last updated on Wednesday 9 Mar 2016 by Hyungmin
subroutine Input_Print_Num_BP_Edge(prob)
    type(ProbType), intent(in) :: prob

    ! --------------------------------------------------
    ! Read edge length that is predefined
    ! The dege that has minimum length is corresponding to the edge lenth
    ! --------------------------------------------------
    if(prob.sel_sec <= para_n_square_lattice) then
        write(0, "(a)")
        write(0, "(a)"), "   C. Third input - Pre-defined minimum edge length"
        write(0, "(a)")
        write(0, "(a)"), "   [Square lattice]"
        write(0, "(a)")
        write(0, "(a)"), "    *  1.  32 BPs :  3 turn ->  3 [turn] * 10.67 [BPs/turn] =  32.00 [10.88nm]"
        write(0, "(a)"), "       2.  43 BPs :  4 turn ->  4 [turn] * 10.67 [BPs/turn] =  42.67 [14.51nm]"
        write(0, "(a)"), "       3.  53 BPs :  5 turn ->  5 [turn] * 10.67 [BPs/turn] =  53.34 [18.14nm]"
        write(0, "(a)"), "    *  4.  64 BPs :  6 turn ->  6 [turn] * 10.67 [BPs/turn] =  64.00 [21.76nm]"
        write(0, "(a)"), "       5.  75 BPs :  7 turn ->  7 [turn] * 10.67 [BPs/turn] =  74.67 [25.39nm]"
        write(0, "(a)"), "       6.  85 BPs :  8 turn ->  8 [turn] * 10.67 [BPs/turn] =  85.34 [29.02nm]"
        write(0, "(a)"), "    *  7.  96 BPs :  9 turn ->  9 [turn] * 10.67 [BPs/turn] =  96.00 [32.64nm]"
        write(0, "(a)"), "       8. 107 BPs : 10 turn -> 10 [turn] * 10.67 [BPs/turn] = 106.67 [36.27nm]"
        write(0, "(a)"), "       9. 117 BPs : 11 turn -> 11 [turn] * 10.67 [BPs/turn] = 117.34 [39.90nm]"
        write(0, "(a)"), "    * 10. 128 BPs : 12 turn -> 12 [turn] * 10.67 [BPs/turn] = 128.00 [43.52nm]"
        write(0, "(a)")
    else
        write(0, "(a)")
        write(0, "(a)"), "   C. Third input - Pre-defined minimum edge length"
        write(0, "(a)")
        write(0, "(a)"), "   [Honeycomb lattice]"
        write(0, "(a)")
        write(0, "(a)"), "       1.  31 BPs :  3 turn ->  3 [turn] * 10.5 [BPs/turn] =  31.5 [10.71 nm]"
        write(0, "(a)"), "    *  2.  42 BPs :  4 turn ->  4 [turn] * 10.5 [BPs/turn] =  42.0 [14.28 nm]"
        write(0, "(a)"), "       3.  52 BPs :  5 turn ->  5 [turn] * 10.5 [BPs/turn] =  52.5 [17.85 nm]"
        write(0, "(a)"), "    *  4.  63 BPs :  6 turn ->  6 [turn] * 10.5 [BPs/turn] =  63.0 [21.42 nm]"
        write(0, "(a)"), "       5.  73 BPs :  7 turn ->  7 [turn] * 10.5 [BPs/turn] =  73.5 [24.99 nm]"
        write(0, "(a)"), "    *  6.  84 BPs :  8 turn ->  8 [turn] * 10.5 [BPs/turn] =  84.0 [28.56 nm]"
        write(0, "(a)"), "       7.  94 BPs :  9 turn ->  9 [turn] * 10.5 [BPs/turn] =  94.5 [32.13 nm]"
        write(0, "(a)"), "    *  8. 105 BPs : 10 turn -> 10 [turn] * 10.5 [BPs/turn] = 105.0 [35.70 nm]"
        write(0, "(a)"), "       9. 115 BPs : 11 turn -> 11 [turn] * 10.5 [BPs/turn] = 115.5 [39.27 nm]"
        write(0, "(a)"), "    * 10. 126 BPs : 12 turn -> 12 [turn] * 10.5 [BPs/turn] = 126.0 [42.84 nm]"
        write(0, "(a)")
    end if

    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Num_BP_Edge

! ---------------------------------------------------------------------------------------

! Set problem to be solved
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Input_Set_Problem(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    if(prob.sel_prob == 100) then

        ! File selector
        call Input_Select_File(prob, geom)
    else

        ! Select the pre-defined example from user input
        call Input_Select_Problem(prob, geom)
    end if
end subroutine Input_Set_Problem

! ---------------------------------------------------------------------------------------

! File selector
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Input_Select_File(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, len_char
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    ! Read geometric file
    write(0, "(a)")
    write(0, "(a)"), " Write the file name (*.ply, *.geo, *.stl, *.wrl), [Enter] : "
    read(*, *),  prob.name_file

    len_char       = LEN_TRIM(prob.name_file)
    prob.type_file = prob.name_file(len_char-2:len_char)
    prob.name_file = prob.name_file(1:len_char-4)

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Select file type
    if(prob.type_file == "ply") then

        ! Import ply format
        call Importer_PLY(prob, geom)
    else if(prob.type_file == "stl") then

        ! Import stl format -> ply format
        call Importer_STL(prob)
        call Importer_PLY(prob, geom)

    else if(prob.type_file == "wrl") then

        ! Import wrl format -> ply format
        call Importer_WRL(prob)
        call Importer_PLY(prob, geom)
    else if(prob.type_file == "geo") then

        ! Import geo format
        call Importer_GEO(prob, geom)
    else

        write(0, "(a$)"), "Error - Not defined file type : "
        write(0, "(a )"), "Input_Select_File"
        stop
    end if

    prob.name_prob = prob.name_file
    prob.name_file = trim(prob.name_file)//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    ! Print filename and type
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* File name : "//trim(prob.name_file)//"."//trim(prob.type_file)
        write(i, "(a)")
    end do
end subroutine Input_Select_File

! ---------------------------------------------------------------------------------------

! Select the pre-defined example from user input
! Last updated on Wednesday 9 Mar 2016 by Hyungmin
subroutine Input_Select_Problem(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set file type as primitive
    prob.type_file = "primitive"

    ! Select problem
    select case (prob.sel_prob)

        ! Pre-defined geometries: Platonic solids
        case (1);  call Exam_Platonic_Tetrahedron  (prob, geom)
        case (2);  call Exam_Platonic_Cube         (prob, geom)
        case (3);  call Exam_Platonic_Octahedron   (prob, geom)
        case (4);  call Exam_Platonic_Dodecahedron (prob, geom)
        case (5);  call Exam_Platonic_Icosahedron  (prob, geom)

        ! Pre-defined geometries: Archimedean solids
        case (6);  call Exam_Archi_Cubeoctahedron          (prob, geom)
        case (7);  call Exam_Archi_Icosidodecahedron       (prob, geom)
        case (8);  call Exam_Archi_Rhombicuboctahedron     (prob, geom)
        case (9);  call Exam_Archi_Snub_Cube               (prob, geom)
        case (10); call Exam_Archi_Truncated_Cube          (prob, geom)
        case (11); call Exam_Archi_Truncated_Cuboctahedron (prob, geom)
        case (12); call Exam_Archi_Truncated_Dodecahedron  (prob, geom)
        case (13); call Exam_Archi_Truncated_Icosahedron   (prob, geom)
        case (14); call Exam_Archi_Truncated_Octahedron    (prob, geom)
        case (15); call Exam_Archi_Truncated_Tetrahedron   (prob, geom)

        ! Pre-defined geometries: Johnson solids
        case (16); call Exam_Johnson_Gyroelongated_Pentagonal_Pyramid_J11   (prob, geom)
        case (17); call Exam_Johnson_Triangular_Bipyramid_J12               (prob, geom)
        case (18); call Exam_Johnson_Pentagonal_Bipyramid_J13               (prob, geom)
        case (19); call Exam_Johnson_Gyroelongated_Square_Bipyramid_J17     (prob, geom)
        case (20); call Exam_Johnson_Square_Gyrobicupola_J29                (prob, geom)
        case (21); call Exam_Johnson_Pentagonal_Orthocupolarotunda_J32      (prob, geom)
        case (22); call Exam_Johnson_Pentagonal_Orthobirotunda_J34          (prob, geom)
        case (23); call Exam_Johnson_Elongated_Pentagonal_Gyrobicupola_J39  (prob, geom)
        case (24); call Exam_Johnson_Elongated_Pentagonal_Gyrobirotunda_J43 (prob, geom)
        case (25); call Exam_Johnson_Gyroelongated_Square_Bicupola_J45      (prob, geom)

        ! Pre-defined geometries: Catalan solids
        case (26); call Exam_Catalan_Rhombic_Dodecahedron        (prob, geom)
        case (27); call Exam_Catalan_Rhombic_Triacontahedron     (prob, geom)
        case (28); call Exam_Catalan_Deltoidal_Icositetrahedron  (prob, geom)
        case (29); call Exam_Catalan_Pentagonal_Icositetrahedron (prob, geom)
        case (30); call Exam_Catalan_Triakis_Octahedron          (prob, geom)
        case (31); call Exam_Catalan_Disdyakis_Dodecahedron      (prob, geom)
        case (32); call Exam_Catalan_Triakis_Icosahedron         (prob, geom)
        case (33); call Exam_Catalan_Pentakis_Dodecahedron       (prob, geom)
        case (34); call Exam_Catalan_Tetrakis_Hexahedron         (prob, geom)
        case (35); call Exam_Catalan_Triakis_Tetrahedron         (prob, geom)

        ! Pre-defined geometries: Miscellaneous polyhedra
        case (36); call Exam_Miscellaneous_Heptagonal_Bipyramid     (prob, geom)
        case (37); call Exam_Miscellaneous_Enneagonal_Trapezohedron (prob, geom)
        case (38); call Exam_Miscellaneous_Small_Stell_Dodecahedron (prob, geom)
        case (39); call Exam_Miscellaneous_Rhombic_Hexecontahedron  (prob, geom)
        case (40); call Exam_Miscellaneous_Goldberg_Dk5dgD          (prob, geom)
        case (41); call Exam_Miscellaneous_Double_Helix             (prob, geom)
        case (42); call Exam_Miscellaneous_Nested_Cube              (prob, geom)
        case (43); call Exam_Miscellaneous_Nested_Octahedron        (prob, geom)
        case (44); call Exam_Miscellaneous_Torus                    (prob, geom)
        case (45); call Exam_Miscellaneous_Double_Torus             (prob, geom)

        ! Pre-defined geometries: Asymmetric and non-convex objects
        case (46); call Exam_Asymmetric_Ball        (prob, geom)
        case (47); call Exam_Asymmetric_Nickedtorus (prob, geom)
        case (48); call Exam_Asymmetric_Helix       (prob, geom)
        case (49); call Exam_Asymmetric_Rod         (prob, geom)
        case (50); call Exam_Asymmetric_Stickman    (prob, geom)
        case (51); call Exam_Asymmetric_Bottle      (prob, geom)
        case (52); call Exam_Asymmetric_Bunny       (prob, geom)

        ! Pre-defined geometries: Ect (3 Asymmetric tetrahedrons, 3 Archimedean solids, 3 Catalan solids)
        case (53); call Exam_Archi_Truncated_Icosidodecahedron      (prob, geom)
        case (54); call Exam_Archi_Rhombicosidodecahedron           (prob, geom)
        case (55); call Exam_Archi_Snub_Dodecahedron                (prob, geom)
        case (56); call Exam_Catalan_Disdyakis_Triacontahedron      (prob, geom)
        case (57); call Exam_Catalan_Deltoidal_Hexecontahedron      (prob, geom)
        case (58); call Exam_Catalan_Pentagonal_Hexecontahedron     (prob, geom)
        case (59); call Exam_Asym_Tetra_I_63_73_52_42               (prob, geom)
        case (60); call Exam_Asym_Tetra_II_63_52_42                 (prob, geom)
        case (61); call Exam_Asym_Tetra_III_63_42                   (prob, geom)
        case (62); call Exam_Johnson_Pentagonal_Pyramid_J2          (prob, geom)
        case (63); call Exam_Johnson_Elongated_Square_Bipyramid_J15 (prob, geom)

        ! Pre-defined geometries: prism
        case (64); call Exam_Prism_Triangular_Prism     (prob, geom)
        case (65); call Exam_Prism_Pentagonal_Prism     (prob, geom)
        case (66); call Exam_Prism_hexagonal_Prism      (prob, geom)

        ! Primitive examples 7 - open geometry
        !case (53); call Exam_OpenGeo_Plate_Uniform_Quad         (prob, geom)
        !case (54); call Exam_OpenGeo_Plate_Distorted_Quad       (prob, geom)
        !case (55); call Exam_OpenGeo_Plate_Uniform_Tri          (prob, geom)
        !case (56); call Exam_OpenGeo_Plate_Distorted_Tri        (prob, geom)
        !case (57); call Exam_OpenGeo_Circular_Plate_Quad        (prob, geom)
        !case (58); call Exam_OpenGeo_Circular_Plate_Tri         (prob, geom)
        !case (59); call Exam_OpenGeo_Annular_Plate_Quad         (prob, geom)
        !case (60); call Exam_OpenGeo_Annular_Plate_Tri          (prob, geom)
        !case (61); call Exam_OpenGeo_Hyperbolic_Paraboloid_Quad (prob, geom)
        !case (62); call Exam_OpenGeo_Hyperbolic_Paraboloid_Tri  (prob, geom)

        case default
            write(0, "(a$)"), "Error - Not defined problem : "
            write(0, "(a )"), "Input_Select_Problem"
        stop
    end select
end subroutine Input_Select_Problem

! ---------------------------------------------------------------------------------------

! Set cross-section based on square or honeycomb lattices
! Last updated on Monday 12 December 2016 by Hyungmin
subroutine Input_Set_Section(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: bp_id

    ! Cross-section definition on local coordinate - t3-t2
    !        t2
    !        ¡è
    !        |
    !     ---|-------¡æ t3
    !        |
    ! The number of columns should be even
    if(prob.sel_sec == 1) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 3 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 3, 13
        !     ¡Ü¡Ü     .00   01.  <------- reference axis
        geom.sec.dir      = 90
        geom.n_sec        = 2
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2

    else if(prob.sel_sec == 2) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 13 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 3, 4 / 13, 14       | caDNAno   02    (CW)
        !      ¡Ü¡Ü        04 03              |         03  01
        !     ¡Ü  ¡Ü      05   02             |         04  00
        !      ¡Ü¡Ü       .00 01.   <--- ref  |           05
        geom.sec.dir      = 90
        geom.n_sec        = 6
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 3; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 1

        ! Increase or decrease edge length except for reference row and column section
        para_vertex_modify = "mod2"

    else if(prob.sel_sec == 3) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 11 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 1 / 11              | caDNAno   02    (CW)
        !      ¡Ü¡Ü        05=04              |         03  01
        !     ¡Ü  ¡Ü     .00   03.  <--- ref  |         04  00
        !      ¡Ü¡Ü        01=02              |           05
        geom.sec.dir      = 150
        geom.n_sec        = 6
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 2; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 1
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 2
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 1

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"
    else if(prob.sel_sec == 4) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 18 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9 / 18, 19       | caDNAno   00    (CCW)
        !      ¡Ü¡Ü        05=04              |         01  05
        !     ¡Ü  ¡Ü     .00   03.  <--- ref  |         02  04
        !      ¡Ü¡Ü        01=02              |           03
        geom.sec.dir      = -90
        geom.n_sec        = 6
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 2; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 1
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 2
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 1

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"
    else if(prob.sel_sec == 5) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 18 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9 / 18, 19       | caDNAno   02    (CW)
        !      ¡Ü¡Ü        04 03              |         03  01
        !     ¡Ü  ¡Ü     .05   02.  <--- ref  |         04  00
        !      ¡Ü¡Ü        00 01              |           05
        geom.sec.dir      = 90
        geom.n_sec        = 6
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 3; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 1

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"
    else

        write(0, "(a$)"), "Error - Not defined cross-section : "
        write(0, "(a )"), "Input_Set_Section"
        stop
    end if

    !! Set section connectivity in the defined initial section
    call Input_Set_Section_Connectivity(prob, geom)

    ! Find maximum and minimum sectional row and column
    call Input_Find_Max_Min_Section(geom)
end subroutine Input_Set_Section

! ---------------------------------------------------------------------------------------

! Set cross-section based on square or honeycomb lattices
! Last updated on Friday 05 August 2016 by Hyungmin
subroutine Input_Set_Section_Old(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: bp_id

    ! Cross-section definition on local coordinate - t3-t2
    !        t2
    !        ¡è
    !        |
    !     ---|-------¡æ t3
    !        |
    ! The number of columns should be even
    if(prob.sel_sec == 1) then
        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü      03 = 02
        !     ¡Ü¡Ü     .00   01.  <------- reference axis
        geom.n_sec        = 4
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 1

    else if(prob.sel_sec == 2) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü     .04   05.
        !     ¡Ü¡Ü      03 = 02
        !     ¡Ü¡Ü     .00   01.  <------- reference axis
        geom.n_sec        = 6
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 1
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 2

    else if(prob.sel_sec == 3) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü      07 = 06
        !     ¡Ü¡Ü     .04   05.
        !     ¡Ü¡Ü      03 = 02
        !     ¡Ü¡Ü     .00   01.  <------- reference axis
        geom.n_sec        = 8
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 1
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 2
        geom.sec.id(7) = 6; geom.sec.posR(7) = 4; geom.sec.posC(7) = 2
        geom.sec.id(8) = 7; geom.sec.posR(8) = 4; geom.sec.posC(8) = 1

    else if(prob.sel_sec == 4) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 4
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 3
        geom.sec.id(4) = 3; geom.sec.posR(4) = 1; geom.sec.posC(4) = 4

    else if(prob.sel_sec == 5) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü      07 = 06   05 = 04
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 8
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 3
        geom.sec.id(4) = 3; geom.sec.posR(4) = 1; geom.sec.posC(4) = 4
        geom.sec.id(5) = 4; geom.sec.posR(5) = 2; geom.sec.posC(5) = 4
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 3
        geom.sec.id(7) = 6; geom.sec.posR(7) = 2; geom.sec.posC(7) = 2
        geom.sec.id(8) = 7; geom.sec.posR(8) = 2; geom.sec.posC(8) = 1

    else if(prob.sel_sec == 6) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü     .08   09 = 10   11.
        !     ¡Ü¡Ü¡Ü¡Ü      07 = 06   05 = 04
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 12
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 2; geom.sec.posC(5)  = 4
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 2; geom.sec.posC(6)  = 3
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 2
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2; geom.sec.posC(8)  = 1
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 3; geom.sec.posC(9)  = 1
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 3; geom.sec.posC(10) = 2
        geom.sec.id(11) = 10; geom.sec.posR(11) = 3; geom.sec.posC(11) = 3
        geom.sec.id(12) = 11; geom.sec.posR(12) = 3; geom.sec.posC(12) = 4

    else if(prob.sel_sec == 7) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü      15 = 14   13 = 12
        !     ¡Ü¡Ü¡Ü¡Ü     .08   09 = 10   11.
        !     ¡Ü¡Ü¡Ü¡Ü      07 = 06   05 = 04
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 16
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 2; geom.sec.posC(5)  = 4
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 2; geom.sec.posC(6)  = 3
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 2
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2; geom.sec.posC(8)  = 1
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 3; geom.sec.posC(9)  = 1
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 3; geom.sec.posC(10) = 2
        geom.sec.id(11) = 10; geom.sec.posR(11) = 3; geom.sec.posC(11) = 3
        geom.sec.id(12) = 11; geom.sec.posR(12) = 3; geom.sec.posC(12) = 4
        geom.sec.id(13) = 12; geom.sec.posR(13) = 4; geom.sec.posC(13) = 4
        geom.sec.id(14) = 13; geom.sec.posR(14) = 4; geom.sec.posC(14) = 3
        geom.sec.id(15) = 14; geom.sec.posR(15) = 4; geom.sec.posC(15) = 2
        geom.sec.id(16) = 15; geom.sec.posR(16) = 4; geom.sec.posC(16) = 1

    else if(prob.sel_sec == 8) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03 = 04   05.  <------- reference axis
        geom.n_sec        = 6
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 6

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 1; geom.sec.posC(5)  = 5
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 1; geom.sec.posC(6)  = 6

    else if(prob.sel_sec == 9) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü      11 = 10   09 = 08   07 = 06
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03 = 04   05.  <------- reference axis
        geom.n_sec        = 12
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 6

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 1; geom.sec.posC(5)  = 5
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 1; geom.sec.posC(6)  = 6
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 6
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2; geom.sec.posC(8)  = 5
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 2; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 2; geom.sec.posC(10) = 3
        geom.sec.id(11) = 10; geom.sec.posR(11) = 2; geom.sec.posC(11) = 2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 2; geom.sec.posC(12) = 1

    else if(prob.sel_sec == 10) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü     .12   13 = 14   15 = 16   17.
        !     ¡Ü¡Ü¡Ü¡Ü      11 = 10   09 = 08   07 = 06
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03 = 04   05.  <------- reference axis
        geom.n_sec        = 18
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 6

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 1; geom.sec.posC(5)  = 5
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 1; geom.sec.posC(6)  = 6
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 6
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2; geom.sec.posC(8)  = 5
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 2; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 2; geom.sec.posC(10) = 3
        geom.sec.id(11) = 10; geom.sec.posR(11) = 2; geom.sec.posC(11) = 2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 2; geom.sec.posC(12) = 1
        geom.sec.id(13) = 12; geom.sec.posR(13) = 3; geom.sec.posC(13) = 1
        geom.sec.id(14) = 13; geom.sec.posR(14) = 3; geom.sec.posC(14) = 2
        geom.sec.id(15) = 14; geom.sec.posR(15) = 3; geom.sec.posC(15) = 3
        geom.sec.id(16) = 15; geom.sec.posR(16) = 3; geom.sec.posC(16) = 4
        geom.sec.id(17) = 16; geom.sec.posR(17) = 3; geom.sec.posC(17) = 5
        geom.sec.id(18) = 17; geom.sec.posR(18) = 3; geom.sec.posC(18) = 6

    else if(prob.sel_sec == 11) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü      23 = 22   21 = 20   19 = 18
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü     .12   13 = 14   15 = 16   17.
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü      11 = 10   09 = 08   07 = 06
        !     ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03 = 04   05.  <------- reference axis
        geom.n_sec        = 24
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 6

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 1; geom.sec.posC(5)  = 5
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 1; geom.sec.posC(6)  = 6
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 6
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2; geom.sec.posC(8)  = 5
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 2; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 2; geom.sec.posC(10) = 3
        geom.sec.id(11) = 10; geom.sec.posR(11) = 2; geom.sec.posC(11) = 2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 2; geom.sec.posC(12) = 1
        geom.sec.id(13) = 12; geom.sec.posR(13) = 3; geom.sec.posC(13) = 1
        geom.sec.id(14) = 13; geom.sec.posR(14) = 3; geom.sec.posC(14) = 2
        geom.sec.id(15) = 14; geom.sec.posR(15) = 3; geom.sec.posC(15) = 3
        geom.sec.id(16) = 15; geom.sec.posR(16) = 3; geom.sec.posC(16) = 4
        geom.sec.id(17) = 16; geom.sec.posR(17) = 3; geom.sec.posC(17) = 5
        geom.sec.id(18) = 17; geom.sec.posR(18) = 3; geom.sec.posC(18) = 6
        geom.sec.id(19) = 18; geom.sec.posR(19) = 4; geom.sec.posC(19) = 6
        geom.sec.id(20) = 19; geom.sec.posR(20) = 4; geom.sec.posC(20) = 5
        geom.sec.id(21) = 20; geom.sec.posR(21) = 4; geom.sec.posC(21) = 4
        geom.sec.id(22) = 21; geom.sec.posR(22) = 4; geom.sec.posC(22) = 3
        geom.sec.id(23) = 22; geom.sec.posR(23) = 4; geom.sec.posC(23) = 2
        geom.sec.id(24) = 23; geom.sec.posR(24) = 4; geom.sec.posC(24) = 1

    else if(prob.sel_sec == 12) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü     .06   07 = 08   09.
        !     ¡Ü  ¡Ü     .05             04.
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 10
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0; geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1; geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2; geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3; geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4; geom.sec.posR(5)  = 2; geom.sec.posC(5)  = 4
        geom.sec.id(6)  = 5; geom.sec.posR(6)  = 2; geom.sec.posC(6)  = 1
        geom.sec.id(7)  = 6; geom.sec.posR(7)  = 3; geom.sec.posC(7)  = 1
        geom.sec.id(8)  = 7; geom.sec.posR(8)  = 3; geom.sec.posC(8)  = 2
        geom.sec.id(9)  = 8; geom.sec.posR(9)  = 3; geom.sec.posC(9)  = 3
        geom.sec.id(10) = 9; geom.sec.posR(10) = 3; geom.sec.posC(10) = 4

    else if(prob.sel_sec == 13) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !     ¡Ü¡Ü¡Ü¡Ü      11 = 10   09 = 08
        !     ¡Ü  ¡Ü     .06             07.
        !     ¡Ü  ¡Ü     .05             04.
        !     ¡Ü¡Ü¡Ü¡Ü     .00   01 = 02   03.  <------- reference axis
        geom.n_sec        = 12
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 2
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 1; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 2; geom.sec.posC(5)  = 4
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 2; geom.sec.posC(6)  = 1
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 3; geom.sec.posC(7)  = 1
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 3; geom.sec.posC(8)  = 4
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 4; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 4; geom.sec.posC(10) = 3
        geom.sec.id(11) = 10; geom.sec.posR(11) = 4; geom.sec.posC(11) = 2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 4; geom.sec.posC(12) = 1

    else if(prob.sel_sec == 14) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !      ¡Ü¡Ü           .06   07.
        !     ¡Ü¡Ü¡Ü¡Ü     .04   05 = 02   03.  <------- reference axis
        !      ¡Ü¡Ü           .00   01.
        geom.n_sec        = 8
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 2
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 3
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 3
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 4
        geom.sec.id(5) = 4; geom.sec.posR(5) = 2; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 2
        geom.sec.id(7) = 6; geom.sec.posR(7) = 3; geom.sec.posC(7) = 2
        geom.sec.id(8) = 7; geom.sec.posR(8) = 3; geom.sec.posC(8) = 3

    else if(prob.sel_sec == 15) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 0
        bp_id = mod(para_start_bp_ID, 32)

        ! Starting BP - 0, 10, 11, 20, 21
        !      ¡Ü¡Ü            11 = 10
        !     ¡Ü¡Ü¡Ü¡Ü      07 = 06   09 = 08
        !     ¡Ü¡Ü¡Ü¡Ü     .04   05 = 02   03.  <------- reference axis
        !      ¡Ü¡Ü           .00   01.
        geom.n_sec        = 12
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "square")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 1; geom.sec.posC(1)  = 2
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 3
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 2; geom.sec.posC(3)  = 3
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 2; geom.sec.posC(4)  = 4
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 2; geom.sec.posC(5)  = 1
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 2; geom.sec.posC(6)  = 2
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 3; geom.sec.posC(7)  = 2
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 3; geom.sec.posC(8)  = 1
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 3; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 3; geom.sec.posC(10) = 3
        geom.sec.id(11) = 10; geom.sec.posR(11) = 4; geom.sec.posC(11) = 3
        geom.sec.id(12) = 11; geom.sec.posR(12) = 4; geom.sec.posC(12) = 2

    else if(prob.sel_sec == 16) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 3 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 3, 13
        !     ¡Ü¡Ü     .00   01.  <------- reference axis
        geom.n_sec        = 2
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2

    else if(prob.sel_sec == 17) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 18 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9, 18, 19
        !      ¡Ü¡Ü        05=04
        !     ¡Ü  ¡Ü     .00   03.  <------- reference axis
        !      ¡Ü¡Ü        01=02
        geom.n_sec        = 6
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 2; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 1
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 2
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 1

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"

    else if(prob.sel_sec == 18) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 13 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 13, 14
        !      ¡Ü¡Ü        04 03
        !     ¡Ü  ¡Ü      05   02
        !      ¡Ü¡Ü       .00 01.     <------- reference axis
        geom.n_sec        = 6
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 3; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 1

        ! Increase or decrease edge length except for reference row and column section
        para_vertex_modify = "mod2"

    else if(prob.sel_sec == 19) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 8 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9, 18, 19
        !      ¡Ü¡Ü  ¡Ü¡Ü       05=04   11=10
        !     ¡Ü  ¡Ü¡Ü  ¡Ü    .00   03=06   09.  <------- reference axis
        !      ¡Ü¡Ü  ¡Ü¡Ü       01=02   07=08
        geom.n_sec        = 12
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 2; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 1
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 2
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 2; geom.sec.posC(4)  = 2
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 3; geom.sec.posC(5)  = 2
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 3; geom.sec.posC(6)  = 1
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 2; geom.sec.posC(7)  = 3
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 1; geom.sec.posC(8)  = 3
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 1; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 2; geom.sec.posC(10) = 4
        geom.sec.id(11) = 10; geom.sec.posR(11) = 3; geom.sec.posC(11) = 4
        geom.sec.id(12) = 11; geom.sec.posR(12) = 3; geom.sec.posC(12) = 3

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"

    else if(prob.sel_sec == 20) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 8 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9, 18, 19
        !        ¡Ü¡Ü              15=14
        !      ¡Ü¡Ü  ¡Ü¡Ü        05=04   11=10
        !     ¡Ü  ¡Ü¡Ü  ¡Ü     .00   03=06   09.  <------- reference axis
        !      ¡Ü¡Ü  ¡Ü¡Ü        01=02   07=08
        !        ¡Ü¡Ü              13=12
        geom.n_sec        = 16
        geom.sec.ref_row  = 3
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 3+1; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 2+1; geom.sec.posC(2)  = 1
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 2+1; geom.sec.posC(3)  = 2
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 3+1; geom.sec.posC(4)  = 2
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 4+1; geom.sec.posC(5)  = 2
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 4+1; geom.sec.posC(6)  = 1
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 3+1; geom.sec.posC(7)  = 3
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 2+1; geom.sec.posC(8)  = 3
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 2+1; geom.sec.posC(9)  = 4
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 3+1; geom.sec.posC(10) = 4
        geom.sec.id(11) = 10; geom.sec.posR(11) = 4+1; geom.sec.posC(11) = 4
        geom.sec.id(12) = 11; geom.sec.posR(12) = 4+1; geom.sec.posC(12) = 3
        geom.sec.id(13) = 12; geom.sec.posR(13) = 1+1; geom.sec.posC(13) = 3
        geom.sec.id(14) = 13; geom.sec.posR(14) = 1+1; geom.sec.posC(14) = 2
        geom.sec.id(15) = 14; geom.sec.posR(15) = 5+1; geom.sec.posC(15) = 3
        geom.sec.id(16) = 15; geom.sec.posR(16) = 5+1; geom.sec.posC(16) = 2

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"

    else if(prob.sel_sec == 21) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 8 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 8, 9, 18, 19
        !      ¡Ü¡Ü  ¡Ü¡Ü        07=06   17=16
        !     ¡Ü  ¡Ü¡Ü  ¡Ü     .08   05=18   15.
        !      ¡Ü¡Ü  ¡Ü¡Ü        09=04   19=14
        !     ¡Ü  ¡Ü¡Ü  ¡Ü     .00   03=10   13.  <------- reference axis
        !      ¡Ü¡Ü  ¡Ü¡Ü        01=02   11=12
        geom.n_sec        = 20
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1)  = 0;  geom.sec.posR(1)  = 2; geom.sec.posC(1)  = 1
        geom.sec.id(2)  = 1;  geom.sec.posR(2)  = 1; geom.sec.posC(2)  = 1
        geom.sec.id(3)  = 2;  geom.sec.posR(3)  = 1; geom.sec.posC(3)  = 2
        geom.sec.id(4)  = 3;  geom.sec.posR(4)  = 2; geom.sec.posC(4)  = 2
        geom.sec.id(5)  = 4;  geom.sec.posR(5)  = 3; geom.sec.posC(5)  = 2
        geom.sec.id(6)  = 5;  geom.sec.posR(6)  = 4; geom.sec.posC(6)  = 2
        geom.sec.id(7)  = 6;  geom.sec.posR(7)  = 5; geom.sec.posC(7)  = 2
        geom.sec.id(8)  = 7;  geom.sec.posR(8)  = 5; geom.sec.posC(8)  = 1
        geom.sec.id(9)  = 8;  geom.sec.posR(9)  = 4; geom.sec.posC(9)  = 1
        geom.sec.id(10) = 9;  geom.sec.posR(10) = 3; geom.sec.posC(10) = 1
        geom.sec.id(11) = 10; geom.sec.posR(11) = 2; geom.sec.posC(11) = 3
        geom.sec.id(12) = 11; geom.sec.posR(12) = 1; geom.sec.posC(12) = 3
        geom.sec.id(13) = 12; geom.sec.posR(13) = 1; geom.sec.posC(13) = 4
        geom.sec.id(14) = 13; geom.sec.posR(14) = 2; geom.sec.posC(14) = 4
        geom.sec.id(15) = 14; geom.sec.posR(15) = 3; geom.sec.posC(15) = 4
        geom.sec.id(16) = 15; geom.sec.posR(16) = 4; geom.sec.posC(16) = 4
        geom.sec.id(17) = 16; geom.sec.posR(17) = 5; geom.sec.posC(17) = 4
        geom.sec.id(18) = 17; geom.sec.posR(18) = 5; geom.sec.posC(18) = 3
        geom.sec.id(19) = 18; geom.sec.posR(19) = 4; geom.sec.posC(19) = 3
        geom.sec.id(20) = 19; geom.sec.posR(20) = 3; geom.sec.posC(20) = 3

        ! Decrease edge length at the bottom
        if(para_vertex_modify == "mod") para_vertex_modify = "mod1"
    else

        write(0, "(a$)"), "Error - Not defined cross-section : "
        write(0, "(a )"), "Input_Set_Section"
        stop
    end if

    !! Set section connectivity in the defined initial section
    call Input_Set_Section_Connectivity(prob, geom)

    ! Find maximum and minimum sectional row and column
    call Input_Find_Max_Min_Section(geom)
end subroutine Input_Set_Section_Old

! ---------------------------------------------------------------------------------------

! Find maximum and minimum sectional row and column
! Last updated on Wednesday 24 Feb 2016 by Hyungmin
subroutine Input_Find_Max_Min_Section(geom)
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! sec.maxR ~ minC were initialized from Mani_Init_SecType
    do i = 1, geom.n_sec
        ! Find max and min of row
        if(geom.sec.posR(i) > geom.sec.maxR) geom.sec.maxR = geom.sec.posR(i)
        if(geom.sec.posR(i) < geom.sec.minR) geom.sec.minR = geom.sec.posR(i)

        ! Find max and min of col
        if(geom.sec.posC(i) > geom.sec.maxC) geom.sec.maxC = geom.sec.posC(i)
        if(geom.sec.posC(i) < geom.sec.minC) geom.sec.minC = geom.sec.posC(i)
    end do

    ! Find the size of columns and rows
    geom.sec.n_row = geom.sec.maxR - geom.sec.minR + 1
    geom.sec.n_col = geom.sec.maxC - geom.sec.minC + 1

    ! Check even number
    if(mod(geom.sec.n_col, 2) /= 0) then
        write(0, "(a$)"), "Error - the column number should be even : "
        write(0, "(a )"), "Input_Find_Max_Min_Section"
        stop
    end if
end subroutine Input_Find_Max_Min_Section

! ---------------------------------------------------------------------------------------

! Set section connectivity in the defined initial section
! Last updated on Wednesday 03 August 2016 by Hyungmin
subroutine Input_Set_Section_Connectivity(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: sec_cur, sec_com, row_cur, row_com
    integer :: i, j, count, id_bp
    logical :: b_connect

    ! If self connection route turns on, find connection to link each other
    ! Section connection was initialized as -1
    ! Loop for current section
    do i = 1, geom.n_sec
        sec_cur = geom.sec.id(i)
        row_cur = geom.sec.posR(i)

        ! Loop for comparing section
        do j = 1, geom.n_sec
            sec_com = geom.sec.id(j)
            row_com = geom.sec.posR(j)

            ! If section numbers are the same
            if(sec_cur == sec_com) cycle

            ! Determine the section connection for scaffold strand
            id_bp     = 1
            b_connect = Section_Connection_Scaf(geom, sec_cur, sec_com, id_bp)

            ! Set section connectivity
            if( (para_vertex_design == "flat"    .and. b_connect == .true.) .or. &
                (para_vertex_design == "beveled" .and. para_vertex_modify == "mod1" .and. &
                 row_cur < geom.sec.ref_row .and. row_cur == row_com ) ) then
                geom.sec.conn(i) = sec_com
                exit
            end if
        end do
    end do

    ! For alternative 6-helice bundle
    if(prob.sel_sec == 2 .and. para_vertex_design == "flat" .and. para_vertex_modify == "mod2") then
        geom.sec.conn(1) = -1       ! Sec 0 -> neighbor
        geom.sec.conn(2) = -1       ! Sec 1 -> neighbor
        geom.sec.conn(3) = 3        ! Sec 2 -> 3
        geom.sec.conn(4) = 2        ! Sec 3 -> 2
        geom.sec.conn(5) = 5        ! Sec 4 -> 5
        geom.sec.conn(6) = 4        ! Sec 5 -> 4
    end if

    ! Print information for self connection route
    count = 0
    write(0, "(a)"), "   --------------------------------------------------"
    write(0, "(a)"), "     Sectional connection to determine route method  "
    write(0, "(a)"), "   --------------------------------------------------"
    do i = 1, geom.n_sec

        write(0, "(i10, a$)"), geom.sec.id(i), " th section  ->"

        if(geom.sec.conn(i) /= -1) then
            write(0, "(i7, a)"), geom.sec.conn(i), "    th section"
        else
            write(0, "(a)"), "  Neighbor connection"
            count = count + 1
        end if
    end do
    write(0, "(a)"), "   --------------------------------------------------"
    write(0, "(a)")

    ! Check neighboring connection
    ! There should be at least two neighboring connection
    ! Also, these number should be even
    if(mod(count, 2) /= 0 .or. count == 0) then
        write(0, "(a$)"), "Error - The neighboring connection should be "
        write(0, "(a$)"), "even and larger than 1 : "
        write(0, "(a )"), "Input_Set_Section_Connectivity"
        stop
    end if
end subroutine Input_Set_Section_Connectivity

! ---------------------------------------------------------------------------------------

! Set length of edge
! Last updated on Thursday 25 Feb 2016 by Hyungmin
subroutine Input_Set_Num_BP_Edge(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(in)    :: geom

    ! One base pair will be added at the end of nodes after FE meshing
    if(geom.sec.types == "square") then

        if(prob.sel_bp_edge == 1)  prob.n_bp_edge = 32       ! 32bp * 1
        if(prob.sel_bp_edge == 2)  prob.n_bp_edge = 43
        if(prob.sel_bp_edge == 3)  prob.n_bp_edge = 53
        if(prob.sel_bp_edge == 4)  prob.n_bp_edge = 64       ! 32bp * 2
        if(prob.sel_bp_edge == 5)  prob.n_bp_edge = 75
        if(prob.sel_bp_edge == 6)  prob.n_bp_edge = 85
        if(prob.sel_bp_edge == 7)  prob.n_bp_edge = 96       ! 32bp * 3
        if(prob.sel_bp_edge == 8)  prob.n_bp_edge = 107
        if(prob.sel_bp_edge == 9)  prob.n_bp_edge = 117
        if(prob.sel_bp_edge == 10) prob.n_bp_edge = 128      ! 32bp * 4

    else if(geom.sec.types == "honeycomb") then

        if(prob.sel_bp_edge == 1)  prob.n_bp_edge = 31
        if(prob.sel_bp_edge == 2)  prob.n_bp_edge = 42       ! 21bp * 2
        if(prob.sel_bp_edge == 3)  prob.n_bp_edge = 52
        if(prob.sel_bp_edge == 4)  prob.n_bp_edge = 63       ! 21bp * 3
        if(prob.sel_bp_edge == 5)  prob.n_bp_edge = 73
        if(prob.sel_bp_edge == 6)  prob.n_bp_edge = 84       ! 21bp * 4
        if(prob.sel_bp_edge == 7)  prob.n_bp_edge = 94
        if(prob.sel_bp_edge == 8)  prob.n_bp_edge = 105      ! 21bp * 5
        if(prob.sel_bp_edge == 9)  prob.n_bp_edge = 115
        if(prob.sel_bp_edge == 10) prob.n_bp_edge = 126      ! 21bp * 6

        ! Unpaired scaffold and staple nucleotides at the vertex when square lattice
        para_unpaired_square = "off"
    end if

    ! Increase the number of basepair for self-connection
    !if(para_vertex_design == 1) then
    !    prob.n_bp_edge = prob.n_bp_edge + 1
    !end if
end subroutine Input_Set_Num_BP_Edge

! ---------------------------------------------------------------------------------------

! Convert surface to line connectivity
! Last updated on Thursday 25 Feb 2016 by Hyungmin
subroutine Input_Convert_Face_To_Line(geom)
    type(GeomType), intent(inout) :: geom
    
    integer :: i, j, k, f_zero, flag
    integer :: point_1, point_2, point_1_com, point_2_com

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)   ! Maximum connectivity
    end type MeshType

    type(MeshType), allocatable, dimension (:) :: Basepair_con  ! 1st: # of meshes, 2nd: points
    type(ListConn), pointer :: straight_con
    type(ListConn), pointer :: ptr, ptr1

    allocate(Basepair_con(geom.n_face))

    ! Nullify the linked list for the junction data
    nullify(straight_con)
    nullify(ptr)
    nullify(ptr1)

    ! Initialze variable
    f_zero      = 0
    geom.n_iniL = 0

    ! Read mesh
    do i = 1, geom.n_face

        ! Read the number of vectices in the mesh and connectivity
        Basepair_con(i).cn(1:geom.face(i).n_poi) = geom.face(i).poi(1:geom.face(i).n_poi)

        ! If there is zero value in the connectivity
        do j = 1, geom.face(i).n_poi
            if(Basepair_con(i).cn(j) == 0) f_zero = 1
        end do

        do j = 1, geom.face(i).n_poi

            flag = 1    ! Flag 1 means that there is no connectivity in the existing array

            if(j == geom.face(i).n_poi) then
                point_1 = Basepair_con(i).cn(j)
                point_2 = Basepair_con(i).cn(1)
            else
                point_1 = Basepair_con(i).cn(j)
                point_2 = Basepair_con(i).cn(j+1)
            end if

            if(i == 1) then     ! First connectivities are always saved
                flag = 1
            else
                allocate(ptr1)
                ptr1 => straight_con

                do k = 1, geom.n_iniL

                    point_1_com = ptr1%point(1)
                    point_2_com = ptr1%point(2)

                    ! Check where there is the same connectivity in the existing array
                    if((point_1_com == point_1 .and. point_2_com == point_2) ) then
                        flag = 0
                        exit
                    else if((point_1_com == point_2 .and. point_2_com == point_1) ) then
                        flag = 0
                        exit
                    end if

                    ! Pointer to move next list
                    ptr1 => ptr1%next
                end do
            end if

            ! Set connectivity adding this into list
            if(flag == 1) then

                allocate(ptr)
                geom.n_iniL  = geom.n_iniL + 1
                ptr%point(1) = point_1
                ptr%point(2) = point_2
                straight_con => List_Insert_Conn(straight_con, ptr)

                if(i == 1 .and. geom.face(i).n_poi == 2) exit
            end if

        end do
    end do

    ! Allocate stragiht line data structure
    allocate(geom.iniL(geom.n_iniL))

    ! Initialize line data type
    call Mani_Init_LineType(geom.iniL, geom.n_iniL)

    ! Set straight line information
    ptr1 => straight_con
    do i = 1, geom.n_iniL

        ! If there is zero in the connectivity, add one
        if(f_zero == 1) then
            geom.iniL(geom.n_iniL+1-i).poi(1:2) = ptr1%point(1:2) + 1
        else
            geom.iniL(geom.n_iniL+1-i).poi(1:2) = ptr1%point(1:2)
        end if

        ! Pointer to move next linked list
        ptr1 => ptr1%next
    end do

    if(f_zero == 1) then
        do i = 1, geom.n_face
            do j = 1, geom.face(i).n_poi
                geom.face(i).poi(j) = geom.face(i).poi(j) + 1
            end do
        end do
    end if

    ! Copy from poi to iniP that always save initial points
    do i = 1, geom.n_iniL
        geom.iniL(i).iniP(1:2) = geom.iniL(i).poi(1:2)
    end do

    ! Delete linked list allocated
    call List_Delete_Conn(straight_con)
    !call List_Delete_Conn(ptr)
    !call List_Delete_Conn(ptr1)

    ! Deallocate mesh connectivity data
    deallocate(Basepair_con)
end subroutine Input_Convert_Face_To_Line

! ---------------------------------------------------------------------------------------

! Set geometric scale, initial minimum length of edge is set up 20 nm
! Last updated on Thuesday 12 Apr 2016 by Hyungmin
subroutine Input_Scale_Init_Geometry(geom)
    type(GeomType), intent(inout) :: geom

    double precision :: pos_c(3), length, min_len
    integer :: i, poi_1, poi_2

    ! Find center position of the structure
    pos_c(1:3) = 0.0d0
    do i = 1, geom.n_iniP
        pos_c(1:3) = pos_c + geom.iniP(i).pos
    end do
    pos_c(1:3) = pos_c / dble(geom.n_iniP)

    ! Move whole geometry to center position
    do i = 1, geom.n_iniP
        geom.iniP(i).pos(1:3) = geom.iniP(i).pos(1:3) - pos_c(1:3)
    end do

    ! Find minimum length of line
    poi_1 = geom.iniL(1).poi(1)
    poi_2 = geom.iniL(1).poi(2)
    min_len = Size_Vector(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

    do i = 2, geom.n_iniL

        poi_1  = geom.iniL(i).poi(1)
        poi_2  = geom.iniL(i).poi(2)
        length = Size_Vector(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

        ! If short edge exists
        if(length < min_len) then
            min_len = length
        end if
    end do

    ! Recalucate edge length
    do i = 1, geom.n_iniP
        geom.iniP(i).pos(1:3) = geom.iniP(i).pos / min_len * para_init_scale
    end do

    ! Find edge that has minimum length 
    !do i = 1, geom.n_iniL
    !    poi_1  = geom.iniL(i).poi(1)
    !    poi_2  = geom.iniL(i).poi(2)
    !    length = Size_Vector(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)
    !    write(0, "(i, f15.5)"), i, length
    !end do
end subroutine Input_Scale_Init_Geometry

! ---------------------------------------------------------------------------------------

! Set working and Chimera path
! Last updated on Wednesday 16 Mar 2016 by Hyungmin
subroutine Input_Set_Path(prob)
    type(ProbType), intent(inout) :: prob

    ! Set working directory
    if(para_external == "on") then
        prob.path_work1 = "Output\"//trim(prob.name_file)//"\"
        prob.path_work2 = "Output/"//trim(prob.name_file)//"/"
    else
        prob.path_work1 = "Output\"
        prob.path_work2 = "Output/"
    end if

    ! Set Chimera path
    prob.path_chimera = trim(para_path_chimera)
end subroutine Input_Set_Path

! ---------------------------------------------------------------------------------------

! Remove previous working directory and make new one
! Last updated on Saturday 20 Feb 2016 by Hyungmin
subroutine Input_Set_Workplace(prob)
    type(ProbType), intent(in) :: prob

    logical :: results

    ! Remove the directory and files
    results = SYSTEMQQ("rd "//trim(prob.path_work1)//' /s /q')

    ! Make new working directory
    results = SYSTEMQQ("md "//trim(prob.path_work1))

    ! Directory for Tecplot output
    if(para_output_Tecplot == "on") then
        results = SYSTEMQQ("md "//trim(prob.path_work1)//"\Tecplot")
    end if

    write(0, "(a)"), "  ...Removed the existing working directory"
    write(0, "(a)")
end subroutine Input_Set_Workplace

! ---------------------------------------------------------------------------------------

! Write geo file file
! Last updated on Friday 4 Mar 2016 by Hyungmin
subroutine Input_Write_GEO_File(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    character(200) :: path
    integer :: i, j

    ! Exception
    if(para_write_101 == .false.) return

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=101, file=trim(path)//".geo", form="formatted")

    ! Write points
    write(101, "(i7)"), geom.n_iniP
    do i = 1, geom.n_iniP
        write(101,"(i7, 3f10.4)"), i, geom.iniP(i).pos(1:3)
    end do
    write(101, "(a)")

    ! Write faces
    write(101, "(i7)"), geom.n_face
    do i = 1, geom.n_face
        write(101,"(i7, i16$)"), geom.face(i).n_poi, geom.face(i).poi(1)
        do j = 2, geom.face(i).n_poi - 1
            write(101,"(i7$)"), geom.face(i).poi(j)
        end do
        write(101,"(i7)"), geom.face(i).poi(geom.face(i).n_poi)
    end do
    write(101, "(a)")

    ! Write edges
    write(101, "(i7)"), geom.n_iniL
    do i = 1, geom.n_iniL
        write(101, "(3i8)"), 2, geom.iniL(i).poi(1:2)
    end do
    write(101, "(a)")

    ! Write geometric for input module
    write(101, "(a)")
    write(101, "(a)"), "--------------------------------------------------"
    write(101, "(a)")
    call Space(101, 4)
    write(101, "(a)") "! Allocate point and face structure"
    call Space(101, 4)
    write(101, "(a, i4)") "geom.n_iniP = ", geom.n_iniP
    call Space(101, 4)
    write(101, "(a, i4)") "geom.n_face = ", geom.n_face
    write(101, "(a)")
    call Space(101, 4)
    write(101, "(a)") "allocate(geom.iniP(geom.n_iniP))"
    call Space(101, 4)
    write(101, "(a)") "allocate(geom.face(geom.n_face))"
    write(101, "(a)")

    ! For points
    call Space(101, 4)
    write(101, "(a)") "! Set point position vectors"
    do i = 1, geom.n_iniP
        call Space(101, 4)
        write(101, "(a, i7, a$)"), "geom.iniP(", i, ").pos(1:3) = [ "
        write(101, "(f10.4, a$)"), geom.iniP(i).pos(1), "d0, "
        write(101, "(f10.4, a$)"), geom.iniP(i).pos(2), "d0, "
        write(101, "(f10.4, a )"), geom.iniP(i).pos(3), "d0 ]"
    end do
    write(101, "(a)")

    ! For faces
    call Space(101, 4)
    write(101, "(a)") "! Set point position vectors"
    do i = 1, geom.n_face
        call Space(101, 4)
        write(101, "(a, i7$)"), "geom.face(",            i
        write(101, "(a, i7$)"), ").n_poi = ",            geom.face(i).n_poi
        write(101, "(a, i7$)"), "; allocate(geom.face(", i
        write(101, "(a, i7$)"), ").poi(",                geom.face(i).n_poi
        write(101, "(a, i7$)"), ")); geom.face(",        i
        write(101, "(a, i7$)"), ").poi(1:",              geom.face(i).n_poi
        write(101, "(a$    )"), ") = ["

        do j = 1, geom.face(i).n_poi - 1
            write(101, "(i7, a$)"), geom.face(i).poi(j), ", "
        end do
        write(101, "(i7, a)"), geom.face(i).poi(geom.face(i).n_poi), " ]"
    end do

    write(101, "(a)")
    write(101, "(a)"), "--------------------------------------------------"
    write(101, "(a)")

    ! For edges
    write(101, "(i7)"), geom.n_iniL
    do i = 1, geom.n_iniL
        write(101, "(a, i7, a$)"), "line(", i, ", 1:2) = [ "
        write(101, "(i7, a$   )"), geom.iniL(i).poi(1), ", "
        write(101, "(i7, a    )"), geom.iniL(i).poi(2), " ]"
    end do
    write(101, "(a)")

    close(unit=101)
end subroutine Input_Write_GEO_File

! ---------------------------------------------------------------------------------------

! Write initial geometry for Chimera
! Last updated on Monday 12 September 2016 by Hyungmin
subroutine Input_Chimera_Init_Geometry(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    integer :: i, j
    logical :: f_info, f_axis
    character(200) :: path

    if(para_write_102 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_102_info

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=102, file=trim(path)//"_init_geo.bild", form="formatted")

    ! Write initial points
    write(102, "(a)"), ".color red"
    do i = 1, geom.n_iniP
        write(102, "(a$    )"), ".sphere "
        write(102, "(3f9.2$)"), geom.iniP(i).pos(1:3)
        write(102, "(1f9.2 )"), 0.75d0
    end do

    ! Write initial edges
    write(102, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL

        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)

        write(102,"(a$    )"), ".cylinder "
        write(102,"(3f9.2$)"), pos_1(1:3)
        write(102,"(3f9.2$)"), pos_2(1:3)
        write(102,"(1f9.2 )"), 0.3d0
    end do

    ! Information on initial geometry
    if(f_info == .true.) then

        ! For points
        do i = 1, geom.n_iniP
            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), geom.iniP(i).pos + 1.0d0
            write(102, "(a    )"), ".color red"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i    )"), i
        end do

        ! For edges
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0
            length     = Size_Vector(pos_2 - pos_1)

            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), pos_c(1:3) + 0.5d0
            write(102, "(a    )"), ".color dark green"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i,    a$)"), i, "("
            write(102, "(f5.1, a )"), length, ")"
        end do

        ! For faces
        do i = 1, geom.n_face

            ! Find center position in mesh
            pos_c(1:3) = 0.0d0
            do j = 1, geom.face(i).n_poi
                pos_c(1:3) = pos_c + geom.iniP(geom.face(i).poi(j)).pos
            end do
            pos_c(1:3) = pos_c / dble(geom.face(i).n_poi)

            ! Write face number
            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), pos_c(1:3) + 1.0d0
            write(102, "(a    )"), ".color black"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i7   )"), i
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(102, "(a)"), ".translate 0.0 0.0 0.0"
        write(102, "(a)"), ".scale 0.5"
        write(102, "(a)"), ".color grey"
        write(102, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(102, "(a)"), ".color red"             ! x-axis
        write(102, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(102, "(a)"), ".color blue"            ! y-axis
        write(102, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(102, "(a)"), ".color yellow"          ! z-axis
        write(102, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=102)

    ! ==================================================
    !
    ! Write the file for Tecplot
    !
    ! ==================================================
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=102, file=trim(path)//"_init_geo.dat", form="formatted")

    write(102, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(102, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(102, "(a$)"), 'ZONE F = FEPOINT'
    write(102, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniP)))
    write(102, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(102, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_iniP
        write(102, "(3f9.3$)"), geom.iniP(i).pos(1:3)
        write(102, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(102, "(1i7$)"), geom.iniL(i).poi(1)
        write(102, "(1i7 )"), geom.iniL(i).poi(2)
    end do

    close(unit=102)
end subroutine Input_Chimera_Init_Geometry

! ---------------------------------------------------------------------------------------

! Write initial geometry for Tecplot
! Last updated on Monday 12 September 2016 by Hyungmin
subroutine Input_Tecplot_Init_Geometry(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    logical :: f_info, f_axis
    integer :: i, j, nline
    character(200) :: path

    if(para_write_103 == .false.) return

    ! Open file
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=102, file=trim(path)//"_init_geo_face.dat", form="formatted")

    ! Find the number of lines
    nline = 0
    do i = 1, geom.n_face
        nline = nline + geom.face(i).n_poi
    end do

    ! For Tecplot output
    write(102, "(a)"), 'VARIABLES = "X", "Y", "Z"'
    write(102, "(a)"), 'ZONE T    = "'//trim(prob.name_file)//'"'
    write(102, "(a)"), 'ZONETYPE  = FEPOLYGON'
    write(102, "(a)"), 'NODES     = '//trim(adjustl(Int2Str(geom.n_iniP)))
    write(102, "(a)"), 'ELEMENTS  = '//trim(adjustl(Int2Str(geom.n_face)))
    write(102, "(a)"), 'FACES     = '//trim(adjustl(Int2Str(nline)))
    write(102, "(a)"), 'NumConnectedBoundaryFaces   = 0'
    write(102, "(a)"), 'TotalNumBoundaryConnections = 0'
    write(102, "(a)")

    ! Write x-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(1)
    end do
    write(102, "(a)")

    ! Write y-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(2)
    end do
    write(102, "(a)")

    ! Write z-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(3)
    end do
    write(102, "(a)"); write(102, "(a)")

    ! Write line connectivity
    write(102, "(a)"), "# Face Nodes List"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            if(j == geom.face(i).n_poi) then
                write(102, "(2i8)"), geom.face(i).poi(j), geom.face(i).poi(1)
            else
                write(102, "(2i8)"), geom.face(i).poi(j), geom.face(i).poi(j+1)
            end if
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    ! # Face Left Elements (In this case, they all point to the single element in this zone)
    write(102, "(a)"), "# Face Left Elements"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            write(102, "(i7$)"), i
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    ! # Face Right elements (0 means no boundary, -n means use the nth boundary connection)
    write(102, "(a)"), "# Face Right elements"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            write(102, "(i7$)"), 0
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    close(unit=102)
end subroutine Input_Tecplot_Init_Geometry

! ---------------------------------------------------------------------------------------

! Generate Schlegel diagram
! Last updated on Sunday 30 May 2016 by Hyungmin
subroutine Input_Generate_Schlegel_Diagram(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision, allocatable :: pos_xy(:,:), pos_nei(:,:)
    integer, allocatable :: face(:), vert_row(:), vert_col(:), nei(:)

    double precision :: angle, pos_mid(2)
    integer :: i, j, k, n_vert, n_rept, n_nei, max_n_vert, max_face, nbr
    logical :: b_continue

    ! Number of vertices and iterations to calculate pos_xy
    n_vert = geom.n_iniP
    n_rept = 10*n_vert

    ! Choose face that has the maximum number of vertices
    max_n_vert = geom.face(1).n_poi
    max_face   = 1
    do i = 2, geom.n_face
        if(geom.face(i).n_poi > max_n_vert) then
            max_n_vert = geom.face(i).n_poi
            max_face   = i
        end if
    end do

    ! Identify vertices associated with biggest face
    allocate(face(max_n_vert))
    do i = 1, max_n_vert
        face(i) = geom.face(max_face).poi(i)
    end do

    ! Initialize pos_xy
    allocate(pos_xy(n_vert, 2))
    pos_xy(1:n_vert, 1:2) = 0.0d0

    ! Set big face on unit circle
    angle = 2.0d0*pi / dble(max_n_vert)

    ! For each vertex in big face
    do i = 1, max_n_vert
        pos_xy(face(i), 1) = dcos(angle*dble(i-1))
        pos_xy(face(i), 2) = dsin(angle*dble(i-1))
    end do

    ! Calculate positions of other vertices through iterative process
    do i = 1, n_rept
        do j = 1, n_vert

            b_continue = .true.
            do k = 1, max_n_vert
                if(face(k) == j) then
                    b_continue = .false.
                    exit
                end if
            end do

            if(b_continue == .true.) then

                ! Find j th vert
                allocate(vert_row(geom.n_iniL))
                allocate(vert_col(geom.n_iniL))

                n_nei = 0
                do k = 1, geom.n_iniL
                    if(geom.iniL(k).poi(1) == j) then
                        n_nei           = n_nei + 1
                        vert_row(n_nei) = k
                        vert_col(n_nei) = 1
                    else if(geom.iniL(k).poi(2) == j) then
                        n_nei           = n_nei + 1
                        vert_row(n_nei) = k
                        vert_col(n_nei) = 2
                    end if
                end do

                ! Find neighbors
                allocate(nei(n_nei))
                nei(1:n_nei) = 0

                do k = 1, n_nei
                    nbr    = 3 - vert_col(k)
                    nei(k) = geom.iniL(vert_row(k)).poi(nbr)
                end do

                ! Get neighbor position
                allocate(pos_nei(n_nei,2))
                pos_nei(:, 1:2) = pos_xy(nei, 1:2)

                ! Find mid position
                pos_mid(1:2) = 0.0d0
                do k = 1, n_nei
                    pos_mid(1:2) = pos_mid(1:2) + pos_nei(k, 1:2)
                end do
                pos_mid(1:2) = pos_mid(1:2)/dble(n_nei)

                ! Store as new pos_xy for jth vertex
                pos_xy(j, 1:2) = pos_mid(1:2)

                ! Deallocate memory
                deallocate(vert_row, vert_col, pos_nei, nei)
            end if
        end do
    end do

    ! Write Schlegel diagram
    call Input_Chimera_Schlegel_Diagram(prob, geom, pos_xy)

    ! Deallocate memory
    deallocate(face)
    deallocate(pos_xy)
end subroutine Input_Generate_Schlegel_Diagram

! ---------------------------------------------------------------------------------------

! Write Schlegel diagram
! Last updated on Sunday 30 May 2016 by Hyungmin
subroutine Input_Chimera_Schlegel_Diagram(prob, geom, pos_xy)
    type(ProbType),   intent(in) :: prob
    type(GeomType),   intent(in) :: geom
    double precision, intent(in) :: pos_xy(:,:)

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j
    logical :: f_axis
    character(200) :: path

    if(para_write_104 == .false.) return

    ! Set option
    f_axis = para_chimera_axis

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=104, file=trim(path)//"_schlegel.bild", form="formatted")

    ! Write initial points
    write(104, "(a)"), ".color red"
    do i = 1, geom.n_iniP
        write(104, "(a$    )"), ".sphere "
        write(104, "(3f9.2$)"), pos_xy(i,1:2) * 30.0d0, 0.0d0
        write(104, "(1f9.2 )"), 0.75d0
    end do

    ! Write initial edges
    write(104, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL
        pos_1(3)   = 0.0d0
        pos_2(3)   = 0.0d0
        pos_1(1:2) = pos_xy(geom.iniL(i).poi(1),1:2) * 30.0d0
        pos_2(1:2) = pos_xy(geom.iniL(i).poi(2),1:2) * 30.0d0

        write(104,"(a$    )"), ".cylinder "
        write(104,"(3f9.2$)"), pos_1(1:3)
        write(104,"(3f9.2$)"), pos_2(1:3)
        write(104,"(1f9.2 )"), 0.3d0
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(104, "(a)"), ".translate 0.0 0.0 0.0"
        write(104, "(a)"), ".scale 0.5"
        write(104, "(a)"), ".color grey"
        write(104, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(104, "(a)"), ".color red"             ! x-axis
        write(104, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(104, "(a)"), ".color blue"            ! y-axis
        write(104, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(104, "(a)"), ".color yellow"          ! z-axis
        write(104, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=104)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=104, file=trim(path)//"_schlegel.dat", form="formatted")

    write(104, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(104, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(104, "(a$)"), 'ZONE F = FEPOINT'
    write(104, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniP)))
    write(104, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(104, "(a )"), ', ET=LINESEG'

    ! Write vertices
    do i = 1, geom.n_iniP
        write(104, "(3f9.3$)"), pos_xy(i, 1:2), 0.0d0
        write(104, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(104, "(1i7$)"), geom.iniL(i).poi(1)
        write(104, "(1i7 )"), geom.iniL(i).poi(2)
    end do

    close(unit=104)
end subroutine Input_Chimera_Schlegel_Diagram

! ---------------------------------------------------------------------------------------

end module Input