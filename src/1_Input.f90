!
! =============================================================================
!
! Module - Input
! Last Updated : 04/10/2018, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX-2L, which allows scientists to build and solve
! the sequence design of complex DNAnanostructures.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! PERDIX-2L is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! PERDIX-2L is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Input

    use Importer

    use Exam_2D_Open

    use Section

    use Para
    use Math
    use List

    implicit none

    public  Input_Initialize
    public  Input_Initialize_Report

    private Input_Print_Parameters
    private Input_Read_Parameter
    private Input_Reset_Para_Report
    private Input_Set_Command
    private Input_Print_Problem
    private Input_Print_Num_BP_Edge
    private Input_Set_Problem
    private Input_Set_Vertex_Design
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

! -----------------------------------------------------------------------------

! Initialize parameters and inputs
subroutine Input_Initialize(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: arg, i, j, n_section, len_char, ppos
    character(10) :: c_sec, c_edge_len, c_edge
    character(100) :: c_prob
    logical :: results, here

    ! Read parameters from env.dat
    call Input_Read_Parameter

    ! Set command environment
    call Input_Set_Command

    if(iargc() == 0) then

        ! ==================================================
        ! Running from Win32 console interface menu
        ! ==================================================
        ! Print pre-defined problems
        call Input_Print_Problem
        read(*, *), c_prob

        ppos = scan(trim(c_prob), ".", BACK = .true.)
        if(ppos > 0) then
            len_char       = len_trim(c_prob)
            prob.name_file = trim(adjustl(c_prob(1:ppos-1)))
            prob.type_file = trim(adjustl(c_prob(ppos+1:len_char)))

            ! Check file
            inquire(file="input/"//trim(prob.name_file)//"."//trim(prob.type_file), EXIST=here)
            if(here == .false.) then
                write(0, "(a)")
                write(0, "(a)"), "   +============================ E R R O R =============================+"
                write(0, "(a)"), "   |                                                                    |"
                write(0, "(a)"), "   |   The file does not exist.                                         |"
                write(0, "(a)"), "   |                                                                    |"
                write(0, "(a)"), "   +====================================================================+"
                write(0, "(a)")
                if(para_platform == "win") pause
                stop
            end if
        else
            read(c_prob, *), prob.sel_prob

            ! The negative value terminate the program
            if(prob.sel_prob <= 0) stop
        end if

        ! Print pre-defined edge length(bps)
        call Input_Print_Num_BP_Edge(prob)
        read(*, *) prob.sel_bp_edge

        ! Reference edge - shortest edge
        prob.sel_edge = 0

        ! Choose specific edge number and take edge-length as an input
        if(prob.sel_bp_edge == 0) then
            write(0, "(a)"), "   Type the specific edge ID [Enter] : "
            write(0, "(a)")
            read(*, *) prob.sel_edge
            write(0, "(a)"), "   Type the minimum edge length [Enter] : "
            write(0, "(a)")
            read(*, *) prob.sel_bp_edge
        end if

        ! The negative value terminate the program
        if(prob.sel_bp_edge < 0) stop
    else

        ! ==================================================
        ! Running from a command shell with options
        ! ==================================================
        arg = 1; call getarg(arg, c_prob)       ! 1st argument, problem
        arg = 2; call getarg(arg, c_edge)       ! 2nd argument, edge
        arg = 3; call getarg(arg, c_edge_len)   ! 3rd argument, edge length

        ppos = scan(trim(c_prob), ".", BACK = .true.)
        if(ppos > 0) then
            len_char       = len_trim(c_prob)
            prob.name_file = trim(adjustl(c_prob(1:ppos-1)))
            prob.type_file = trim(adjustl(c_prob(ppos+1:len_char)))

            ! Check file
            inquire(file="input/"//trim(prob.name_file)//"."//trim(prob.type_file), EXIST=here)
            if(here == .false.) then
                write(0, "(a)")
                write(0, "(a)"), "   +============================ E R R O R =============================+"
                write(0, "(a)"), "   |                                                                    |"
                write(0, "(a)"), "   |   The file does not exist.                                         |"
                write(0, "(a)"), "   |                                                                    |"
                write(0, "(a)"), "   +====================================================================+"
                write(0, "(a)")
                if(para_platform == "win") pause
                stop
            end if
        else
            read(c_prob, *), prob.sel_prob

            ! The negative value terminate the program
            if(prob.sel_prob <= 0) stop
        end if

        ! Reference edge
        read(c_edge, *), prob.sel_edge

        ! Edge length
        read(c_edge_len, *), prob.sel_bp_edge
    end if

    ! Set section - DX tile, vertex - mitered
    prob.sel_sec    = 1
    prob.sel_vertex = 2

    ! ==================================================
    ! Set problem, section, edge length and vertex design
    ! ==================================================
    ! Set vertex design
    call Input_Set_Vertex_Design(prob)

    ! Set cross-section
    call Input_Set_Section(prob, geom)

    ! Set the minimum edge length
    call Input_Set_Num_BP_Edge(prob, geom)

    ! Set problem
    call Input_Set_Problem(prob, geom)

    ! ==================================================
    ! Prepair geometry - line generation and scaling
    ! ==================================================
    ! Convert surface to line connectivity
    call Input_Convert_Face_To_Line(geom)

    ! Set geometric scale with initial minimum length
    call Input_Scale_Init_Geometry(geom)

    ! ==================================================
    ! Set environment and write initial geometry
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

    ! Print progress
    call Input_Print_Parameters(prob, geom)
end subroutine Input_Initialize

! -----------------------------------------------------------------------------

! Initialize all parameters
subroutine Input_Initialize_Report(prob, geom, mesh, i, sec, edge, char_vert, char_cut)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom
    type(MeshType), intent(inout) :: mesh
    integer,        intent(inout) :: sec
    integer,        intent(in)    :: i, edge
    character(10),  intent(in)    :: char_vert, char_cut

    ! Reset parameters
    call Input_Reset_Para_Report
    para_platform        = "dev"
    para_vertex_design   = char_vert    ! Flat or mitered vertex
    para_cut_stap_method = char_cut     ! Staple-break

    ! Reset data structures
    prob.n_cng_min_stap = 0
    prob.n_cng_max_stap = 0
    mesh.n_node         = 0
    mesh.n_ele          = 0

    ! Set command environment
    call Input_Set_Command

    ! Set parameters of problem
    prob.sel_prob    = i
    prob.sel_sec     = sec
    prob.sel_bp_edge = edge

    ! Set UCSF Chimera output control
    para_write_101   = .false.       !  GEO file
    para_write_102   = .true.        ! *Initial geometry
    para_write_103   = .false.       !  Faced initial geometry
    para_write_104   = .false.       !  Schlegel diagram
    para_write_301   = .false.       !  Initial geometry with face orientation
    para_write_302   = .true.        ! *Initial geometry with local vector
    para_write_303   = .true.        ! *Modified geometry seperated from vertex
    para_write_401   = .false.       !  Cross-sectional geometry
    para_write_501   = .false.       !  Cylindrical model with orientation
    para_write_502   = .true.        ! *Cylindrical model
    para_write_503   = .false.       !  Basepair mode,
    para_write_504   = .true.        ! *Multiple lines
    para_write_505   = .true.        ! *Txt file on edge length
    para_write_601_1 = .false.       !  Route 1, seperated edges
    para_write_601_2 = .false.       !  Route 2, contruction closed loop
    para_write_601_3 = .false.       !  Route 3, centered crossovers
    para_write_601_4 = .false.       !  Route 4, modified centered crossovers
    para_write_601_5 = .false.       !  Route 5, scaffold route
    para_write_606   = .true.        ! *Sapnning tree for dual-graph
    para_write_607   = .true.        ! *Crossovers based on basepair model
    para_write_608   = .false.       !  3-orientation vectors
    para_write_609   = .false.       !  Atomic model without sequence design
    para_write_610   = .false.       !  Possible centered scaffold crossovers
    para_write_701   = .true.        ! *Txt on sequence design data
    para_write_711   = .false.       !  Csv file for sequence data
    para_write_702   = .true.        ! *Atomic model with sequence design
    para_write_703   = .true.        ! *Route 6, strand route with nick
    para_write_705   = .true.        ! *Route design
    para_write_706   = .false.       !  Atomic model bases on strands/sequence
    para_write_710   = .false.       !  Edge-based sequence design
    para_write_801   = .false.       !  Txt on basepair based data
    para_write_802   = .false.       !  Txt on nucleotide based data
    para_write_803   = .true.        ! *CanDo input file
    para_write_804   = .false.       !  Tecplot input file
    para_write_805   = .false.       !  ADINA input file
    para_write_808   = .false.       !  Txt on sectional edges based sequence

    ! UCSF Chimera output option
    para_chimera_axis     = .false.  !  Plot with axis at the ceneter of geometry
    para_chimera_102_info = .true.   ! *Plot with edge and point number
    para_chimera_301_info = .false.  !  Plot with edge and point number
    para_chimera_302_info = .true.   ! *Plot with edge and point number
    para_chimera_303_info = .true.   ! *Plot with edge and point number
    para_chimera_401_info = .false.  !  Plot with edge and point number
    para_chimera_502_ori  = .false.  !  Plot with helix z-direction
    para_chimera_503_mod  = .false.  !  Plot with modified edges
    para_chimera_504_info = .true.   ! *Plot with edge and point number
    para_chimera_601_dir  = .false.  !  Plot with strand direction
    para_chimera_609_cyl  = .false.  !  Plot with cylinderical representation)
    para_chimera_609_dir  = .false.  !  Plot with strand direction

    ! ==================================================
    ! Set problem, cross-section and edge length
    ! ==================================================
    ! Set cross-section
    call Input_Set_Section(prob, geom)

    ! Set the minimum edge length
    call Input_Set_Num_BP_Edge(prob, geom)

    ! Set problem
    call Input_Set_Problem(prob, geom)

    ! ==================================================
    ! Prepair geometry - line generation and scaling
    ! ==================================================
    ! Convert surface to line connectivity
    call Input_Convert_Face_To_Line(geom)

    ! Set geometric scale with initial minimum length
    call Input_Scale_Init_Geometry(geom)

    ! ==================================================
    ! Set environment and write initial geometry
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

    ! Print progress
    call Input_Print_Parameters(prob, geom)
end subroutine Input_Initialize_Report

! -----------------------------------------------------------------------------

! Print progress
subroutine Input_Print_Parameters(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Open output progress file (unit 11 is used for global output file)
    open(unit=11, file=trim(prob.path_work)//"/TXT_PERDIX_2L.txt", form="formatted")

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
end subroutine Input_Print_Parameters

! -----------------------------------------------------------------------------

! Read parameters from external txt file, env.txt
subroutine Input_Read_Parameter
    character(200) :: ctemp
    integer :: i

    ! Open file
    open(unit=1, file="env.txt", form="formatted")

    read(1, *), ctemp, para_platform

    ! If the mode is not dev
    if(para_platform /= "dev") then
        read(1, *), ctemp, para_set_seq_scaf
    end if

    close(unit=1)
end subroutine Input_Read_Parameter

! -----------------------------------------------------------------------------

! Reset paramters as default values
subroutine Input_Reset_Para_Report

    ! Program parameters
    para_preset          = "on"       ! [on, off], Preset parameter defined in pre-defined examples
    para_output_Tecplot  = "on"       ! [on, off], Output files for Tecplot(http://www.tecplot.com/) to draw vector image
    para_fig_view        = "xy"       ! [xy, xz, xyz, all], Viewpoint for figures from UCSF Chimera
    para_type_cndo       = 2          ! [1, 2], CanDo file option, 1 : original format, 2 : updated format

    ! Parameters for junction modification
    para_junc_ang        = "opt"      ! [opt, max, ave, min], Junction gap modification for different arm angle
    para_const_edge_mesh = "off"      ! [off, on], Constant edge length from polyhedra mesh
    para_sticky_self     = "off"      ! [off, on], Sticky-end for self connection on henycomb cross-section
    para_unpaired_scaf   = "on"       ! [on, off], Unpaired scaffold nucleotides
    para_vertex_modify   = "const"    ! [const, mod], Vertex modification to avoid clash
    !para_vertex_design   = "flat"     ! [flat, mitered], Vertex design

    ! Paramters for B-from DNA generation
    para_dist_pp       = 0.42d0     ! [0.42, 0.6], distance between adjacent phosphate groups, nm
    para_dist_bp       = 0.34d0     ! [0.34 ], Axial rise distance, nm
    para_rad_helix     = 1.0d0      ! [1.0  ], The radius of the DNA helix, nm
    para_gap_helix     = 0.25d0     ! [0.25 ], The gap between two helixes, nm
    para_ang_minor     = 150.0d0    ! [150.0], An angle of minor groove, degree
    para_ang_correct   = 0.0d0      ! [0.0  ], Correction factor to adjust orientation, degree
    para_n_base_tn     = -1         ! [-1   ], The number of nucleotides in poly T loop, -1 : depending on distance
    para_start_bp_ID   = -1         ! [-1   ], Starting base pair ID for the reference, -1 : pre-defined starting BP

    ! Paramters for scaffold route
    para_weight_edge   = "on"       ! [on, off], Assign weight factor into edges of dual graph
    para_method_MST    = "prim"     ! [prim, kruskal, greedy], Minimum spanning tree algorithm
    para_method_sort   = "quick"    ! [none, quick, shell], Sorting algorithm to find MST for Prim or Kruskal
    para_adjacent_list = "off"      ! [off, on], Output for adjacent list for Prim or Kruskal
    para_all_spanning  = "off"      ! [off, on], All possible spanning trees when # of edges is less than 12 for Prim or Kruskal

    ! Parameter for sequence design
    !para_cut_stap_method  = "max"   ! [max, mix, opt, min, mid], Cutting method to make short staple strand, opt - 14nt seeds
    para_set_stap_sxover  = "off"   ! [off, on], To make non-circular staple by single crossover (when para_set_stap_sxover is "on")
    para_output_design    = "arrow" ! [arrow, seq, strand], Graphical output type for sequence design
    para_set_xover_scaf   = "split" ! [split, center], Setting possible scaffold strand

    para_gap_xover_two_scaf   = 3   ! [3 ], The minimum gap between two scaffold crossovers
    para_gap_xover_bound_scaf = 7   ! [7 ], The mimimum gap between scaffold crossover and vertex boundary
    para_gap_xover_bound_stap = 6   ! [6 ], The mimimum gap between staple crossover and vertex boundary
    para_gap_xover_two        = 6   ! [6 ], The minimum gap between scaffold and staple crossovers
    para_gap_xover_nick1      = 10  ! [10], The minimum gap between xover(scaf/stap)/Tn and first nick
    para_gap_xover_nick       = 3   ! [3 ], The minimum gap between xover and nick, if staple length exceeds 60, redesign with num - 1

    para_max_cut_scaf         = 0   ! [0, 7249], Scaffold break - 0 : not breaking, num : breaking over num
    para_min_cut_stap         = 20  ! [20], The minimum number of nucleotides for one staple strand
    para_mid_cut_stap         = 40  ! [40], The optimal number of nucleotides for one staple strand
    para_max_cut_stap         = 60  ! [60], The maximum number of nucleotides for one staple strand
    para_set_seq_scaf         = 0   ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from env.txt, 2 - random
    para_set_start_scaf       = 1   ! [1], Starting nucleotide position of scaffold strand
end subroutine Input_Reset_Para_Report

! -----------------------------------------------------------------------------

! Set command environment
subroutine Input_Set_Command
    integer :: results

    ! Set command environments
    if(para_platform == "dev" .or. para_platform == "win") then
        results = systemqq('title PERDIX-2L')                   ! cmd title
        results = systemqq('mode con: cols=135 lines=6000')     ! cmd size
        results = systemqq('color')                             ! convert color, 02, f0, f1, f2
        results = systemqq('date /t')                           ! display time
        !results = systemqq('hostname')                          ! display hostname of the computer
        !results = systemqq('ver')                               ! display version information
    end if
end subroutine Input_Set_Command

! -----------------------------------------------------------------------------

! Print pre-defined problems
subroutine Input_Print_Problem
    write(0, "(a)")
    write(0, "(a)"), "   +============================================================================+"
    write(0, "(a)"), "   |                                                                            |"
    write(0, "(a)"), "   |    PERDIX-2L by Hyungmin Jun (hyungminjun@outlook.com), Bathe Lab, MIT     |"
    write(0, "(a)"), "   |                                                                            |"
    write(0, "(a)"), "   +============================================================================+"
    write(0, "(a)")
    write(0, "(a)"), "   A. First input - Pre-defined 2D target geometries"
    write(0, "(a)"), "   ================================================="
    write(0, "(a)")
    write(0, "(a)"), "    [ Triangular-mesh objects ]"
    write(0, "(a)"), "       1. Square,               2. Honeycomb"
    write(0, "(a)"), "       3. Circle,               4. Wheel,                    5. Ellipse"
    write(0, "(a)")
    write(0, "(a)"), "    [ Quadrilateral-mesh objects ]"
    write(0, "(a)"), "       6. Rhombic Tiling,       7. Quarter Circle"
    write(0, "(a)"), "       8. Cross,                9. Arrowhead,               10. Annulus"
    write(0, "(a)")
    write(0, "(a)"), "    [ N-polygon-mesh objects ]"
    write(0, "(a)"), "      11. Cairo Penta Tiling,  12. Lotus"
    write(0, "(a)"), "      13. Hexagonal Tiling,    14. Prismatic Penta Tiling,  15. Hepta Penta Tiling"
    write(0, "(a)")
    write(0, "(a)"), "    [ Variable vertex-number, edge-length, and internal mesh ]"
    write(0, "(a)"), "      16. 4-Sided Polygon,     17. 5-Sided Polygon,         18. 6-Sided Polygon"
    write(0, "(a)"), "      19. L-Shape [42-bp],     20. L-Shape [63-bp],         21. L-Shape [84-bp]"
    write(0, "(a)"), "      22. Curved Arm [Quad],   23. Curved Arm [Tri],        24. Curved Arm [Mixed]"
    write(0, "(a)")
    write(0, "(a)"), "   Select the number or type geometry file (*.geo, *.igs, *.iges, *.ply) [Enter]: "
    write(0, "(a)")
end subroutine Input_Print_Problem

! -----------------------------------------------------------------------------

! Print pre-defined minimum edge lengths
subroutine Input_Print_Num_BP_Edge(prob)
    type(ProbType), intent(in) :: prob

    ! The minimum edge lengths pre-defined
    write(0, "(a)")
    write(0, "(a)"), "   B. Second input - Pre-defined minimum edge lengths"
    write(0, "(a)"), "   =================================================="
    write(0, "(a)")
    write(0, "(a)"), "   * 1.  42 bp =  4 turn * 10.5 bp/turn ->  42 bp * 0.34nm/bp = 14.28nm"
    write(0, "(a)"), "     2.  52 bp =  5 turn * 10.5 bp/turn ->  52 bp * 0.34nm/bp = 17.85nm"
    write(0, "(a)"), "   * 3.  63 bp =  6 turn * 10.5 bp/turn ->  63 bp * 0.34nm/bp = 21.42nm"
    write(0, "(a)"), "     4.  73 bp =  7 turn * 10.5 bp/turn ->  73 bp * 0.34nm/bp = 24.99nm"
    write(0, "(a)"), "   * 5.  84 bp =  8 turn * 10.5 bp/turn ->  84 bp * 0.34nm/bp = 28.56nm"
    write(0, "(a)"), "     6.  94 bp =  9 turn * 10.5 bp/turn ->  94 bp * 0.34nm/bp = 32.13nm"
    write(0, "(a)"), "   * 7. 105 bp = 10 turn * 10.5 bp/turn -> 105 bp * 0.34nm/bp = 35.70nm"
    write(0, "(a)"), "     8. 115 bp = 11 turn * 10.5 bp/turn -> 115 bp * 0.34nm/bp = 39.27nm"
    write(0, "(a)"), "   * 9. 126 bp = 12 turn * 10.5 bp/turn -> 126 bp * 0.34nm/bp = 42.84nm"
    write(0, "(a)")
    write(0, "(a)"), "   Select the number or type the minimum edge length [Enter] : "
    write(0, "(a)")
end subroutine Input_Print_Num_BP_Edge

! -----------------------------------------------------------------------------

! Set problem
subroutine Input_Set_Problem(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    if(prob.sel_prob == 0) call Input_Select_File(prob, geom)
    if(prob.sel_prob /= 0) call Input_Select_Problem(prob, geom)
end subroutine Input_Set_Problem

! -----------------------------------------------------------------------------

! Set vertex design
subroutine Input_Set_Vertex_Design(prob)
    type(ProbType), intent(inout) :: prob

    if(prob.sel_vertex == 1) then
        para_vertex_design = "flat"
    else
        para_vertex_design = "mitered"
    end if

    !print *, para_vertex_design
end subroutine Input_Set_Vertex_Design

! -----------------------------------------------------------------------------

! Select geometry file
subroutine Input_Select_File(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, len_char

    ! Read geometric file
    !write(0, "(a)")
    !write(0, "(a)"), " Write the file name (*.PLY, *.IGES, *.GEO), [Enter] : "
    !read(*, *), prob.name_file

    !len_char = LEN_TRIM(prob.name_file)
    !if(prob.name_file(len_char-3:len_char-3) == '.') then
    !
    !    prob.type_file = prob.name_file(len_char-2:len_char)
    !    prob.name_file = prob.name_file(1:len_char-4)
    !
    !else if(prob.name_file(len_char-4:len_char-4) == '.') then
    !    prob.type_file = prob.name_file(len_char-3:len_char)
    !    prob.name_file = prob.name_file(1:len_char-5)
    !else
    !    write(0, "(a)"), "Wrong file type"
    !    stop
    !end if

    ! Select file format
    if(prob.type_file == "ply") then
        call Importer_PLY(prob, geom)
    else if(prob.type_file == "geo" .or. prob.type_file == "igs" .or. prob.type_file == "iges") then
        call Importer_GEO(prob, geom)
    else
        write(0, "(a)")
        write(0, "(a)"), "   +============================ E R R O R =============================+"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   |   The file extension is not supported.                             |"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   +====================================================================+"
        write(0, "(a)")
        if(para_platform == "win") pause
        stop
    end if

    prob.name_prob = prob.name_file
    call Mani_Set_Problem(prob, [52, 152, 219], "xy")

    ! Print filename and type
    call Space(0, 11)
    write(0, "(a)"), "* File name : "//trim(prob.name_file)//"."//trim(prob.type_file)
    write(0, "(a)")
end subroutine Input_Select_File

! -----------------------------------------------------------------------------

! Select the pre-defined geometry
subroutine Input_Select_Problem(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set file type as primitive
    prob.type_file = "primitive"

    ! Select problem
    select case (prob.sel_prob)

        ! Triangular-mesh objects
        case ( 1); call Exam_Open2D_Square          (prob, geom)
        case ( 2); call Exam_Open2D_Honeycomb       (prob, geom)
        case ( 3); call Exam_Open2D_Circle          (prob, geom)
        case ( 4); call Exam_Open2D_Wheel           (prob, geom)
        case ( 5); call Exam_Open2D_Ellipse         (prob, geom)

        ! Quadrilateral-mesh objects
        case ( 6); call Exam_Open2D_Rhombic_Tiling  (prob, geom)
        case ( 7); call Exam_Open2D_Quarter_Circle  (prob, geom)
        case ( 8); call Exam_Open2D_Cross           (prob, geom)
        case ( 9); call Exam_Open2D_Arrowhead       (prob, geom)
        case (10); call Exam_Open2D_Annulus         (prob, geom)

        ! 5-, 6-, 7-, 8-polygon-mesh object
        case (11); call Exam_Open2D_Cairo_Penta_Tiling     (prob, geom)
        case (12); call Exam_Open2D_Lotus                  (prob, geom)
        case (13); call Exam_Open2D_Hexagonal_Tiling       (prob, geom)
        case (14); call Exam_Open2D_Prismatic_Penta_Tiling (prob, geom)
        case (15); call Exam_Open2D_Hepta_Penta_Tiling     (prob, geom)

        ! Variable vertex-number
        case (16); call Exam_Open2D_4_Sided_Polygon (prob, geom)
        case (17); call Exam_Open2D_5_Sided_Polygon (prob, geom)
        case (18); call Exam_Open2D_6_Sided_Polygon (prob, geom)

        ! Variable edge-length
        case (19); call Exam_Open2D_L_Shape_42bp (prob, geom)
        case (20); call Exam_Open2D_L_Shape_63bp (prob, geom)
        case (21); call Exam_Open2D_L_Shape_84bp (prob, geom)

        ! Variable internal mesh
        case (22); call Exam_Open2D_Curved_Arm_Quad (prob, geom)
        case (23); call Exam_Open2D_Curved_Arm_Tri  (prob, geom)
        case (24); call Exam_Open2D_Curved_Arm_Mix  (prob, geom)

        case default
        write(0, "(a)")
        write(0, "(a)"), "   +============================ E R R O R =============================+"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   |   The pre-defined geometries does not exist.                       |"
        write(0, "(a)"), "   |   Select the number from 1 to 24, or type geometry file            |"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   +====================================================================+"
        write(0, "(a)")
        if(para_platform == "win") pause
        stop
    end select
end subroutine Input_Select_Problem

! -----------------------------------------------------------------------------

! Set cross-section
subroutine Input_Set_Section(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: bp_id

    ! The cross-section is defined on the local coordinate, t3-t2
    !        t2
    !        ¡è
    !        |
    !     ---|------¡æ t3
    !        |
    ! The number of columns of the crosssection should be even
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

! -----------------------------------------------------------------------------

! Find maximum and minimum sectional row and column
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

! -----------------------------------------------------------------------------

! Set section connectivity in the defined initial section
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
                (para_vertex_design == "mitered" .and. para_vertex_modify == "mod1" .and. &
                 row_cur < geom.sec.ref_row .and. row_cur == row_com ) ) then
                geom.sec.conn(i) = sec_com
                exit
            end if
        end do
    end do

    ! Print information for self connection route
    count = 0
    write(0, "(a)"), "   --------------------------------------------------"
    write(0, "(a)"), "     Sectional connection to determine route method  "
    write(0, "(a)"), "   --------------------------------------------------"
    do i = 1, geom.n_sec

        write(0, "(i10, a$)"), geom.sec.id(i), " section  ->"

        if(geom.sec.conn(i) /= -1) then
            write(0, "(i7, a)"), geom.sec.conn(i), "    section"
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

! -----------------------------------------------------------------------------

! Set the minimum edge length
subroutine Input_Set_Num_BP_Edge(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(in)    :: geom

    if(prob.sel_bp_edge == 1) prob.n_bp_edge =  42     ! 10.5bp/turn *  4 turn
    if(prob.sel_bp_edge == 2) prob.n_bp_edge =  52     ! 10.5bp/turn *  5 turn
    if(prob.sel_bp_edge == 3) prob.n_bp_edge =  63     ! 10.5bp/turn *  6 turn
    if(prob.sel_bp_edge == 4) prob.n_bp_edge =  73     ! 10.5bp/turn *  7 turn
    if(prob.sel_bp_edge == 5) prob.n_bp_edge =  84     ! 10.5bp/turn *  8 turn
    if(prob.sel_bp_edge == 6) prob.n_bp_edge =  94     ! 10.5bp/turn *  9 turn
    if(prob.sel_bp_edge == 7) prob.n_bp_edge = 105     ! 10.5bp/turn * 10 turn
    if(prob.sel_bp_edge == 8) prob.n_bp_edge = 115     ! 10.5bp/turn * 11 turn
    if(prob.sel_bp_edge == 9) prob.n_bp_edge = 126     ! 10.5bp/turn * 12 turn

    if(prob.sel_bp_edge >= 10 .and. prob.sel_bp_edge <= 36) then
        write(0, "(a)")
        write(0, "(a)"), "   +============================ E R R O R =============================+"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   |  The minimum edge length should be over 37-bp.                     |"
        write(0, "(a)"), "   |                                                                    |"
        write(0, "(a)"), "   +====================================================================+"
        write(0, "(a)")
        if(para_platform == "win") pause
        stop
    end if

    if(prob.sel_bp_edge >= 21) then
        prob.n_bp_edge = prob.sel_bp_edge
    end if
end subroutine Input_Set_Num_BP_Edge

! -----------------------------------------------------------------------------

! Convert surface to line connectivity
subroutine Input_Convert_Face_To_Line(geom)
    type(GeomType), intent(inout) :: geom
    
    integer :: i, j, k, f_zero, flag
    integer :: point_1, point_2, point_1_com, point_2_com

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)   ! Maximum connectivity
    end type MeshType

    type(MeshType), allocatable :: Basepair_con(:)  ! 1st: # of meshes, 2nd: points
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

! -----------------------------------------------------------------------------

! Set geometric scale, initial minimum length of edge is set up 20 nm
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
    min_len = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

    do i = 2, geom.n_iniL

        poi_1  = geom.iniL(i).poi(1)
        poi_2  = geom.iniL(i).poi(2)
        length = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

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
    !    length = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)
    !    write(0, "(i, f15.5)"), i, length
    !end do
end subroutine Input_Scale_Init_Geometry

! -----------------------------------------------------------------------------

! Set working and Chimera path
subroutine Input_Set_Path(prob)
    type(ProbType), intent(inout) :: prob

    ! Set working directory
    if(para_platform == "dev") then
        prob.path_work = "output"
    else if(para_platform == "win") then
        prob.path_work = "output\"//trim(prob.name_file)
    else
        prob.path_work = "output/"//trim(prob.name_file)
    end if
end subroutine Input_Set_Path

! -----------------------------------------------------------------------------

! Remove previous working directory and make new one
subroutine Input_Set_Workplace(prob)
    type(ProbType), intent(in) :: prob

    logical :: results

    ! Remove the directory and files
    if(para_platform == "dev") then
        results = systemqq("rd "//trim(prob.path_work)//' /s /q')   ! Windows
    else if(para_platform == "win") then
        results = systemqq("rd "//trim(prob.path_work)//' /s /q')   ! Windows
    else
        results = systemqq("rm -r "//trim(prob.path_work))          ! Linux & Mac
    end if

    ! Make new working directory
    if(para_platform == "dev") then
        results = systemqq("md "//trim(prob.path_work))         ! Dev
    else if(para_platform == "win") then
        results = systemqq("md "//trim(prob.path_work))         ! Windows
    else
        results = systemqq("mkdir -p "//trim(prob.path_work))   ! Linux & Mac
    end if

    ! Directory for Tecplot output
    if(para_output_Tecplot == "on") then
        if(para_platform == "dev" )then
            results = systemqq("md "//trim(prob.path_work)//"\tecplot")         ! Dev
        else if(para_platform == "win") then
            results = systemqq("md "//trim(prob.path_work)//"\tecplot")         ! Windows
        else
            results = systemqq("mkdir -p "//trim(prob.path_work)//"/tecplot")   ! Linux & Mac
        end if
    end if

    write(0, "(a)"), "  ...Removed the existing working directory"
    write(0, "(a)")
end subroutine Input_Set_Workplace

! -----------------------------------------------------------------------------

! Write geo file file
subroutine Input_Write_GEO_File(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    character(200) :: path
    integer :: i, j

    ! Exception
    if(para_write_101 == .false.) return

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
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
    write(101, "(a)") "! The number of points and faces"
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
    write(101, "(a)") "! Set face connectivity"
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

! -----------------------------------------------------------------------------

! Write initial geometry for Chimera
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

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit=102, file=trim(path)//"_01_target_geometry.bild", form="formatted")

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
            length     = Norm(pos_2 - pos_1)

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
    ! Write the FE format output
    ! ==================================================
    !open(unit=102, file=trim(path)//"_19_FE_Format.txt", form="formatted")
    !write(102, "(i)"), geom.n_iniP
    !do i = 1, geom.n_iniP
    !    write(102, "(4f10.3, 2i10)"), geom.iniP(i).pos(1:2), 0.0d0, 0.0d0, 1, 1
    !end do

    !write(102, "(i)"), geom.n_iniL
    !do i = 1, geom.n_iniL
    !    write(102, "(2i10)"), geom.iniL(i).poi(1), geom.iniL(i).poi(2)
    !end do
    !close(unit=102)

    ! ==================================================
    ! Write the file for Tecplot
    ! ==================================================
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit=102, file=trim(path)//"_01_init_geo.dat", form="formatted")

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

! -----------------------------------------------------------------------------

! Write initial geometry for Tecplot
subroutine Input_Tecplot_Init_Geometry(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    logical :: f_info, f_axis
    integer :: i, j, nline
    character(200) :: path

    if(para_write_103 == .false.) return

    ! Open file
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
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

! -----------------------------------------------------------------------------

! Generate Schlegel diagram
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

! -----------------------------------------------------------------------------

! Write Schlegel diagram
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

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
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

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
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

! -----------------------------------------------------------------------------

end module Input