!
! ---------------------------------------------------------------------------------------
!
!                                   Designer-Open
!
!                                                                    Updated : 2017/03/18
!
! Comments: 
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com).
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
program Designer_Open

    use Ifport

    use Data_Prob       ! Data structure for problem description
    use Data_Geom       ! Data strucutre for face, edge, point and section
    use Data_Mesh       ! Data structure for basepair
    use Data_DNA        ! Data structure for B-form DNA

    use Para            ! Parameters used

    use Input           ! 1st step : Input for geometry and cross-section
    use ModGeo          ! 2nd step : Seperated line from the vertex
    use Section         ! 3rd step : Multiple lines seperated from vertex
    use Basepair        ! 4th step : Basepair design from cross-sectional geometry
    use Route           ! 5th step : B-form DNA genertion and scaffold routing
    use SeqDesign       ! 6th step : Sequence design
    use Output          ! 7th step : Outputs for Chimera and Tecplot

    implicit none

    call Main           ! Main module
    !call Autorun        ! Main module for auto run
    !call Autorun_CHK    ! Main module for auto run

    contains

! ---------------------------------------------------------------------------------------

! Main subroutine
! Last updated on Sun 12 Mar 2017 by Hyungmin
subroutine Main()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem data
    type(GeomType)  :: geom     ! Geometry data
    type(BoundType) :: bound    ! Boundary data
    type(MeshType)  :: mesh     ! Basepaire data
    type(DNAType)   :: dna      ! B-form DNA data
    real            :: time     ! Computational time

    ! 1st step : Initialize input
    call Input_Initialization(prob, geom)

    ! Check time
    call cpu_time(time)

    ! 2nd step : Modified geometry seperated from vertex
    call ModGeo_Modification(prob, geom, bound)

    ! 3rd step : Cross-sectional geometry
    call Section_Generation(prob, geom, bound)

    ! 4th step : Basepair generation from cross-sectional geometry
    call Basepair_Discretize(prob, geom, bound, mesh)

    ! 5th step : B-form DNA and scaffold route
    call Route_Generation(prob, geom, bound, mesh, dna)

    ! 6th step : Sequence design
    call SeqDesign_Design(prob, geom, mesh, dna)

    ! 7th step : Generate outputs and run post-processing tools
    call Output_Generation(prob, mesh, dna)

    ! Print information
    call Print_Information(prob, geom, bound, mesh, dna)

    ! Deallocate global dynamic array
    call Deallocate_Variables(geom, bound, mesh, dna)

    ! Check time consuming
    call Print_TimeConsuming(time)

    !if(para_external == .true.) pause
end subroutine Main

! ---------------------------------------------------------------------------------------

! Autorun to calculate base length in scaffold and staple strands
! Last updated on Friday 11 November 2016 by Hyungmin
subroutine Autorun()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem description
    type(GeomType)  :: geom     ! Geometric data(section, point, edge, face)
    type(BoundType) :: bound    ! Boundary data(outer, junction)
    type(MeshType)  :: mesh     ! Fintie element node data
    type(DNAType)   :: dna      ! B-form DNA data

    double precision :: time, srt_time, end_time
    integer :: arg, i, sec, edge, vertex, max_stap, min_stap
    logical :: results, b_external
    character(10) :: char_sec, char_edge, char_cut, char_junc, char_vert

    b_external = .true.    ! Flag for external run

    if(b_external == .true.) then

        ! Input system by using main argument
        arg = 1; call getarg(arg, char_sec)     ! Argument, section
        arg = 2; call getarg(arg, char_edge)    ! Argument, edge length
        arg = 3; call getarg(arg, char_cut)     ! Argument, cutting method
        arg = 4; call getarg(arg, char_junc)    ! Argument, junction
        arg = 5; call getarg(arg, char_vert)    ! Argument, # of vertex designs

        read(char_sec,  *), sec
        read(char_edge, *), edge
        read(char_vert, *), vertex

    else if(b_external == .false.) then

        sec       = 2       ! Section number
        edge      = 3       ! Edge length
        vertex    = 1       ! Non-beveled and beveled vertex design
        char_cut  = "mix"
        char_junc = "opt"
    end if

    ! Open file
    open(unit=90, file="Autorun_"//trim(adjustl(Int2Str(sec)))//&
        "_"//trim(adjustl(Int2Str(edge)))//"_"//trim(char_cut)//"_"//trim(char_junc)//&
        "_"//trim(adjustl(Int2Str(vertex)))//".txt", form="formatted")

    ! Remove the directory and files
    results = SYSTEMQQ("rd "//trim("Output")//' /s /q')

    ! Check time
    call cpu_time(srt_time)

    ! Section infomation
    write(90, "(a)"), "================================="
    if(edge == 2) write(90, "(a)"), "Edge length           : 42bp"
    if(edge == 3) write(90, "(a)"), "Edge length           : 52bp"
    if(edge == 4) write(90, "(a)"), "Edge length           : 63bp"
    write(90, "(a)"), "Staple cutting        : "//trim(char_cut)
    write(90, "(a)"), "Junction modification : "//trim(char_junc)
    if(sec == 1)  write(90, "(a)"), "Cross-section         : square"
    if(sec == 16) write(90, "(a)"), "Cross-section         : DX tile"
    if(sec == 17) write(90, "(a)"), "Cross-section         : honeycomb"
    if(sec == 18) write(90, "(a)"), "Cross-section         : honeycomb"
    write(90, "(a)"), "================================="
    write(90, "(a)")

    write(90, "(a)"), "L_Scaf - Required scaffold length(nt)"
    write(90, "(a)"), "L_Stap - Required total staple length(nt)"
    write(90, "(a)"), "Len_BP - Base paire length(nt)"
    write(90, "(a)"), "MaxE   - Maximum edge length"
    write(90, "(a)"), "MinE   - Minimum edge length"
    write(90, "(a)"), "Sec    - Cross-section number"
    write(90, "(a)"), "#Stap  - # of staples"
    write(90, "(a)"), "Min    - Minimum length of staples"
    write(90, "(a)"), "Max    - Maximum length of staples"
    write(90, "(a)"), "14nt   - # of 14nt seeds"
    write(90, "(a)"), "S14    - # of secondary 14nt seeds"
    write(90, "(a)"), "4nt    - # of 4nt windows"
    write(90, "(a)"), "R14    - # of 14nt seeds / # of staples"
    write(90, "(a)"), "RS14   - # of secondary 14nt seeds / # of staples"
    write(90, "(a)"), "RN14   - # of no 14nt seeds / # of staples"
    write(90, "(a)"), "#TR    - # of total regions"
    write(90, "(a)"), "TR14   - # of total 14nt seeds"
    write(90, "(a)"), "TR4    - # of total 4nt regions"
    write(90, "(a)"), "a1     - # of paramter changing for minimum staple length"
    write(90, "(a)"), "a2     - # of paramter changing for maximum staple length"
    write(90, "(a)"), "b1     - # of scaffold crossovers"
    write(90, "(a)"), "b2     - # of staple crossovers"
    write(90, "(a)"), "b3     - # of staple single crossovers"
    write(90, "(a)"), "c1     - # of unpaired scaffold nucleotides"
    write(90, "(a)"), "c1     - # of unpaired staple nucleotides"
    write(90, "(a)")

    call space (90, 30)
    write(90, "(a)"), "L_Scaf| L_Stap| Len_BP| MaxE| MinE| Sec| #Stap| Min| Max| 14nt| S14| 4nt|  R14| RS14| RN14|  #TR| TR14|  TR4| a1| a2|   b1|    b2|  b3|   c1|   c2"

    ! Print border
    write(90, "(a$)"), "---------------------------- "
    write(90, "(a )"), " =================================================================================================================================================="

    min_stap = 10000
    max_stap =-10000

    ! Loop for problem
    do i = 1, 62

        ! Problem section
        if( i/=1  .and. i/=2  .and. i/=3  .and. i/=4  .and. i/=5  .and. i/=6  .and. i/=8  .and. i/=10 .and. i/=14 .and. i/=15 .and. &
            i/=16 .and. i/=17 .and. i/=18 .and. i/=19 .and. i/=20 .and. i/=26 .and. i/=34 .and. i/=35 .and. i/=59 ) then

            !write(90, "(a2, a$)"), trim(adjustl(Int2Str(i))), ". ------------------------, "
            !write(90, "(a     )"), " -----, ------, ------, ----, ----, ---, -----, ---, ---, ----, ---, ---, ----, ----, ----, ----, ----, ----, --, --, ----, -----, ---, ----, ----"
            !cycle
        end if

        ! For main problems
        if( i >= 46 ) then
            !write(90, "(a2, a$)"), trim(adjustl(Int2Str(i))), ". ------------------------, "
            !write(90, "(a     )"), " -----, ------, ------, ----, ----, ---, -----, ---, ---, ----, ---, ---, ----, ----, ----, ----, ----, ----, --, --, ----, -----, ---, ----, ----"
            !cycle
        end if

        ! Initialize input
        call Input_Initialization_Autorun(prob, geom, i, sec, edge, vertex, char_cut, char_junc)

        ! Set parameters
        para_write_101   = .false.; para_write_102   = .false.; para_write_103   = .false.
        para_write_104   = .false.; para_write_301   = .false.; para_write_302   = .false.
        para_write_303   = .false.; para_write_401   = .false.; para_write_501   = .false.
        para_write_502   = .false.; para_write_503   = .false.; para_write_504   = .false.
        para_write_505   = .true. ; para_write_601_1 = .false.; para_write_601_2 = .false.
        para_write_601_3 = .false.; para_write_601_4 = .false.; para_write_601_5 = .false.
        para_write_606   = .false.; para_write_607   = .false.; para_write_608   = .false.
        para_write_609   = .false.; para_write_610   = .false.; para_write_701   = .true.
        para_write_702   = .false.; para_write_703   = .false.; para_write_705   = .false.
        para_write_706   = .false.; para_write_710   = .false.; para_write_801   = .false.
        para_write_802   = .false.; para_write_803   = .false.; para_write_804   = .false.
        para_write_805   = .false.; para_write_808   = .false.

        ! 3rd step : Modifed geometry seperated from vertex
        call ModGeo_Modification(prob, geom, bound)

        ! 4th step : Cross-sectional geometry
        call Section_Generation(prob, geom, bound)

        ! 5th step : Basepair generation from cross-sectional geometry
        call Basepair_Discretize(prob, geom, bound, mesh)

        ! 6th step : B-form DNA generation and scaffold route
        call Route_Generation(prob, geom, bound, mesh, dna)

        ! 7th step : Sequence design
        call SeqDesign_Design(prob, geom, mesh, dna)

        ! Generate outputs and run post-processing tools
        call Output_Generation(prob, mesh, dna)

        ! Print information
        call Print_Information(prob, geom, bound, mesh, dna)

        ! Deallocate global dynamic array
        call Deallocate_Variables(geom, bound, mesh, dna)

        ! Print problem
        write(90, "(a4, a24, a$)"), trim(adjustl(Int2Str(i)))//". ", trim(prob.name_prob), ","

        ! Print information
        write(90, "(a8$)"), trim(adjustl(Int2Str(dna.n_base_scaf)))//","
        write(90, "(a8$)"), trim(adjustl(Int2Str(dna.n_base_stap)))//","
        write(90, "(a8$)"), trim(adjustl(Int2Str(mesh.n_node)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(geom.max_edge_length)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(geom.min_edge_length)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(prob.sel_sec)))//","
        write(90, "(a7$)"), trim(adjustl(Int2Str(dna.n_stap)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.len_min_stap)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.len_max_stap)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_14nt)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_s14nt)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_4nt)))//","
        write(90, "(a6$)"), trim(adjustl(Dble2Str2(dble(dna.n_14nt)/dble(dna.n_stap))))//","
        write(90, "(a6$)"), trim(adjustl(Dble2Str2(dble(dna.n_s14nt)/dble(dna.n_stap))))//","
        write(90, "(a6$)"), trim(adjustl(Dble2Str2(dble(dna.n_stap-dna.n_14nt)/dble(dna.n_stap))))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_tot_region)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_tot_14nt)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_tot_4nt)))//","
        write(90, "(a4$)"), trim(adjustl(Int2Str(prob.n_cng_min_stap)))//","
        write(90, "(a4$)"), trim(adjustl(Int2Str(prob.n_cng_max_stap)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_xover_scaf)))//","
        write(90, "(a7$)"), trim(adjustl(Int2Str(dna.n_xover_stap)))//","
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_sxover_stap)))//","
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_nt_unpaired_scaf)))//","
        write(90, "(a5 )"), trim(adjustl(Int2Str(dna.n_nt_unpaired_stap)))

        if(max_stap < dna.len_max_stap) max_stap = dna.len_max_stap
        if(min_stap > dna.len_min_stap) min_stap = dna.len_min_stap
    end do

    write(90, "(a)")
    write(90, "(a)"), "Minimum staple length   : "//trim(adjustl(Int2Str(min_stap)))
    write(90, "(a)"), "Maximum staple length   : "//trim(adjustl(Int2Str(max_stap)))

    ! Set end time
    call cpu_time(end_time)
    time = end_time - srt_time
    write(90, "(a)"), "Time consuming for test : "// trim(adjustl(Dble2Str(time/60.0d0)))//" [min]"
    write(90, "(a)"); write(90, "(a)")

    ! Close file
    close(unit=90)
end subroutine Autorun

! ---------------------------------------------------------------------------------------

! Autorun to calculate base length in scaffold and staple strands
! Last updated on Monday 12 December 2016 by Hyungmin
subroutine Autorun_CHK()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem description
    type(GeomType)  :: geom     ! Geometric data(section, point, edge, face)
    type(BoundType) :: bound    ! Boundary data(outer, junction)
    type(MeshType)  :: mesh     ! Fintie element node data
    type(DNAType)   :: dna      ! B-form DNA data

    integer :: i, j, sec
    character(10) :: char_cut, char_junc

    ! Open file
    open(unit=90, file="Autorun_Output.txt", form="formatted")

    ! Loop for problem
    do i = 1, 45
        if(i == 5 .or. i == 6) then

            ! Edge length
            do j = 2, 6, 2

                sec       = 3
                char_cut  = "max"
                char_junc = "opt"

                ! Initialize input
                call Input_Initialization_Autorun(prob, geom, i, sec, j, 1, char_cut, char_junc)

                ! Set parameters
                para_write_101   = .false.; para_write_102   = .false.; para_write_103   = .false.
                para_write_104   = .false.; para_write_301   = .false.; para_write_302   = .false.
                para_write_303   = .false.; para_write_401   = .false.; para_write_501   = .false.
                para_write_502   = .false.; para_write_503   = .false.; para_write_504   = .false.
                para_write_505   = .true. ; para_write_601_1 = .false.; para_write_601_2 = .false.
                para_write_601_3 = .false.; para_write_601_4 = .false.; para_write_601_5 = .false.
                para_write_606   = .false.; para_write_607   = .false.; para_write_608   = .false.
                para_write_609   = .false.; para_write_610   = .false.; para_write_701   = .true.
                para_write_702   = .false.; para_write_703   = .false.; para_write_705   = .false.
                para_write_706   = .false.; para_write_710   = .false.; para_write_801   = .false.
                para_write_802   = .false.; para_write_803   = .false.; para_write_804   = .false.
                para_write_805   = .false.; para_write_808   = .false.

                ! 3rd step : Modifed geometry seperated from vertex
                call ModGeo_Modification(prob, geom, bound)

                ! 4th step : Cross-sectional geometry
                call Section_Generation(prob, geom, bound)

                ! 5th step : Basepair generation from cross-sectional geometry
                call Basepair_Discretize(prob, geom, bound, mesh)

                ! 6th step : B-form DNA generation and scaffold route
                call Route_Generation(prob, geom, bound, mesh, dna)

                ! 7th step : Sequence design
                call SeqDesign_Design(prob, geom, mesh, dna)

                ! Generate outputs and run post-processing tools
                call Output_Generation(prob, mesh, dna)

                ! Print information
                call Print_Information(prob, geom, bound, mesh, dna)

                ! Deallocate global dynamic array
                call Deallocate_Variables(geom, bound, mesh, dna)

                ! Print problem
                write(90, "(i6, i6, a)"), i, j, "   "//trim(adjustl(Dble2Str(dble(dna.n_14nt)/dble(dna.n_stap)*100.0d0)))
                write(90, "(i6, i6, a)"), i, j, "   "//trim(adjustl(Dble2Str(dble(dna.n_nt_14nt)/dble(dna.n_base_stap)*100.0d0)))
                write(90, "(a)")
            end do
        end if
    end do

    ! Close file
    close(unit=90)
end subroutine Autorun_CHK

! ---------------------------------------------------------------------------------------

! Print information of the DNAcs
! Last updated on Thursday 14 June 2016 by Hyungmin
subroutine Print_Information(prob, geom, bound, mesh, dna)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound
    type(MeshType),  intent(in) :: mesh
    type(DNAType),   intent(in) :: dna

    integer :: i

    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |   8. Summary of DNA Nanostructures with Arbitrary Cross-section    |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")

        call Space(i, 6)
        write(i, "(a)"), "8.1. Geometry, cross-section and problem description"
        call Space(i, 11)
        write(i, "(a)"), "* Full geometric name               : "//trim(prob.name_file)
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial faces       : "//trim(adjustl(Int2Str(geom.n_face)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial points      : "//trim(adjustl(Int2Str(geom.n_iniP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial edges       : "//trim(adjustl(Int2Str(geom.n_iniL)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified points     : "//trim(adjustl(Int2Str(geom.n_modP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified edges      : "//trim(adjustl(Int2Str(geom.n_iniL)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of sectional points    : "//trim(adjustl(Int2Str(geom.n_croP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of sectional edges     : "//trim(adjustl(Int2Str(geom.n_croL)))
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
        call Space(i, 11)
        write(i, "(a)"), "* The minimum edge-length           : "//trim(adjustl(Int2Str(prob.n_bp_edge)))// "bp"
        call Space(i, 11)
        write(i, "(a)"), "* Junction modification             : "//trim(para_junc_ang)
        call Space(i, 11)
        write(i, "(a)"), "* Constant edge length from mesh    : "//trim(para_const_edge_mesh)
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
        write(i, "(a$)"), "* Maximum # of bases in scaf strand : "
        if(para_max_cut_scaf == 0) then
            write(i, "(a)"), "infinite"
        else
            write(i, "(a)"), trim(adjustl(Int2Str(para_max_cut_scaf)))
        end if
        write(i, "(a)")

        call Space(i, 6)
        write(i, "(a)"), "8.2. Scaffold routing and sequence design"
        call Space(i, 11)
        write(i, "(a$)"), "* Scaffold sequence                 : "
        if(para_set_seq_scaf == 0) then
            write(i, "(a)"), "M13mp18(7249nt) sequence"
        else if(para_set_seq_scaf == 1) then
            write(i, "(a)"), "user-defined sequence from seq.txt"
        else if(para_set_seq_scaf == 2) then
            write(i, "(a)"), "random sequence"
        end if
        call Space(i, 11)
        write(i, "(a)"), "* Vertex design method              : "//trim(para_vertex_design)//" vertex"
        call Space(i, 11)
        write(i, "(a)"), "* Vertex design to avoid clash      : "//trim(para_vertex_modify)
        call Space(i, 11)
        write(i, "(a)"), "* Cutting method for short staples  : "//trim(para_cut_stap_method)
        call Space(i, 11)
        write(i, "(a)"), "* Non-circular stap by single xover : "//trim(para_set_stap_sxover)
        call Space(i, 11)
        write(i, "(a)"), "* Unpaired scaffold nucleotides     : "//trim(para_unpaired_scaf)
        call Space(i, 11)
        write(i, "(a)"), "* Start position of scaffold strand : "//trim(adjustl(Int2Str(para_set_start_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap btw xover and nick [staple]   : "//trim(adjustl(Int2Str(para_gap_xover_nick)))
        call Space(i, 11)
        write(i, "(a)"), "* Gap btw xover and vertex [staple] : "//trim(adjustl(Int2Str(para_gap_xover_bound_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of basepairs                    : "//trim(adjustl(Int2Str(mesh.n_node)))
        call Space(i, 11)
        write(i, "(a)"), "* # of nucleotides in scaffold      : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of nucleotides in staple        : "//trim(adjustl(Int2Str(dna.n_base_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of scaffold crossovers          : "//trim(adjustl(Int2Str(dna.n_xover_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of all staple crossovers        : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of single staple crossovers     : "//trim(adjustl(Int2Str(dna.n_sxover_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of nucleotides in scaffold      : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of unpaired scaffolds           : "//trim(adjustl(Int2Str(dna.n_unpaired_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of unpaired nt in scaffolds     : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of unpaired staples             : "//trim(adjustl(Int2Str(dna.n_unpaired_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of unpaired nt in staples       : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Minimum staple length             : "//trim(adjustl(Int2Str(dna.len_min_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Maximum staple length             : "//trim(adjustl(Int2Str(dna.len_max_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Average staple length             : "//trim(adjustl(Dble2Str(dna.len_ave_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of changing for min staple      : "//trim(adjustl(Int2Str(prob.n_cng_min_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of changing for max staple      : "//trim(adjustl(Int2Str(prob.n_cng_max_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of scaffold strands             : "//trim(adjustl(Int2Str(dna.n_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* # of staple strands               : "//trim(adjustl(Int2Str(dna.n_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of strands with 4nt seeds       : "//trim(adjustl(Int2Str(dna.n_4nt)))
        call Space(i, 11)
        write(i, "(a)"), "* # of strands with 14nt seeds      : "//trim(adjustl(Int2Str(dna.n_14nt)))
        call Space(i, 11)
        write(i, "(a)"), "* The ratio - #strand_4nts/#staps   : "//trim(adjustl(Dble2Str(dble(dna.n_4nt)/dble(dna.n_stap)*100.0d0)))//" %"
        call Space(i, 11)
        write(i, "(a)"), "* The ratio - #strand_14nts/#staps  : "//trim(adjustl(Dble2Str(dble(dna.n_14nt)/dble(dna.n_stap)*100.0d0)))//" %"
        call Space(i, 11)
        write(i, "(a)"), "* # of total regions in staples     : "//trim(adjustl(Int2Str(dna.n_tot_region)))
        call Space(i, 11)
        write(i, "(a)"), "* # of total 14nt seeds             : "//trim(adjustl(Int2Str(dna.n_tot_14nt)))
        call Space(i, 11)
        write(i, "(a)"), "* # of total 4nt seeds              : "//trim(adjustl(Int2Str(dna.n_tot_4nt)))
        call Space(i, 11)
        write(i, "(a)"), "* # of total nucleotides in 14nt    : "//trim(adjustl(Int2Str(dna.n_nt_14nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_nt_14nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"
        call Space(i, 11)
        write(i, "(a)"), "* # of total nucleotides in 4nt     : "//trim(adjustl(Int2Str(dna.n_nt_4nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_nt_4nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"
        call Space(i, 11)
        write(i, "(a)"), "* Minimum edge length               : "//trim(adjustl(Int2Str(geom.min_edge_length)))//" bp"
        call Space(i, 11)
        write(i, "(a)"), "* Maximum edge length               : "//trim(adjustl(Int2Str(geom.max_edge_length)))//" bp"
        write(i, "(a)")
        write(i, "(a)"), "   --------------------------------------------------------------------------------"
        write(i, "(a)")
        write(i, "(a)"), " Completed"
    end do

    close(unit=11)
end subroutine Print_Information

! ---------------------------------------------------------------------------------------

! Deallocate global dynamic array
! Last updated on Wednesday 24 Feb 2016 by Hyungmin
subroutine Deallocate_Variables(geom, bound, mesh, dna)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    type(MeshType),  intent(inout) :: mesh
    type(DNAType),   intent(inout) :: dna

    integer :: i

    ! Deallocate geom array
    do i = 1, geom.n_face
        deallocate(geom.face(i).poi)
    end do

    if(allocated(geom.sec.id))   deallocate(geom.sec.id)
    if(allocated(geom.sec.posR)) deallocate(geom.sec.posR)
    if(allocated(geom.sec.posC)) deallocate(geom.sec.posC)
    if(allocated(geom.sec.conn)) deallocate(geom.sec.conn)
    if(allocated(geom.iniP))     deallocate(geom.iniP)
    if(allocated(geom.modP))     deallocate(geom.modP)
    if(allocated(geom.croP))     deallocate(geom.croP)
    if(allocated(geom.iniL))     deallocate(geom.iniL)
    if(allocated(geom.croL))     deallocate(geom.croL)
    if(allocated(geom.face))     deallocate(geom.face)

    ! Deallocate bound array
    do i = 1, bound.n_outer
        deallocate(bound.outer(i).neiP)
        deallocate(bound.outer(i).newP)
    end do

    do i = 1, bound.n_junc
        deallocate(bound.junc(i).iniL)
        deallocate(bound.junc(i).modP)
        deallocate(bound.junc(i).croP)
        deallocate(bound.junc(i).node)
        deallocate(bound.junc(i).conn)
        deallocate(bound.junc(i).type_conn)
    end do

    if(allocated(bound.outer)) deallocate(bound.outer)
    if(allocated(bound.junc))  deallocate(bound.junc)

    ! Deallocate mesh arrays
    if(allocated(mesh.node)) deallocate(mesh.node)
    if(allocated(mesh.ele))  deallocate(mesh.ele)

    ! Deallocate dna arrays
    do i = 1, dna.n_strand
        if(allocated(dna.strand(i).base)) deallocate(dna.strand(i).base)
    end do

    if(allocated(dna.base_scaf))  deallocate(dna.base_scaf)
    if(allocated(dna.base_stap))  deallocate(dna.base_stap)
    if(allocated(dna.top))        deallocate(dna.top)
    if(allocated(dna.strand))     deallocate(dna.strand)
    if(allocated(dna.order_stap)) deallocate(dna.order_stap)
end subroutine Deallocate_Variables

! ---------------------------------------------------------------------------------------

! Print time consuming
! Last updated on Thursday 18 Feb 2016 by Hyungmin
subroutine Print_TimeConsuming(time_start)
    real, intent(in) :: time_start

    real :: time, time_end

    ! Set end time
    call cpu_time(time_end)

    time = time_end - time_start

    write(0, "(a       )")
    write(0, "(a       )"), " --------------------------------------------------"
    write(0, "(a$      )"), "   Time consuming : "
    write(0, "(f6.2, a$)"), time, " [sec], "
    write(0, "(f6.2, a )"), time/60.0d0, " [min]"
    write(0, "(a       )"), " --------------------------------------------------"

    write(0, "(a)")
end subroutine Print_TimeConsuming

! ---------------------------------------------------------------------------------------

end program Designer_Open