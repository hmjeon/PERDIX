!
! ---------------------------------------------------------------------------------------
!
!                                   PERDIX-OPEN
!
!                                                                    Updated : 2017/08/02
!
! Comments: PERDIX is an open-source Fortran library, which allows scientists
! to build and solve the sequence design of complex DNA nanostructures.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
program PERDIX_OPEN

    use Ifport

    use Data_Prob       ! Data structure for problem definition
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

    !call Main           ! Main module
    call Report         ! Main module for auto run

contains

! ---------------------------------------------------------------------------------------

! Main subroutine
subroutine Main()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem data
    type(GeomType)  :: geom     ! Geometry data
    type(BoundType) :: bound    ! Boundary data
    type(MeshType)  :: mesh     ! Basepaire data
    type(DNAType)   :: dna      ! B-form DNA data
    real            :: time     ! Computational time

    ! 1st step : Initialize input
    call Input_Initialize(prob, geom)

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

    ! Verify solution
    call Verify_Solution(mesh, dna)

    ! Deallocate global dynamic array
    call Deallocate_Variables(geom, bound, mesh, dna)

    ! Check time consuming
    call Print_TimeConsuming(time)

    !if(para_external == .true.) pause
end subroutine Main

! ---------------------------------------------------------------------------------------

! Autorun to calculate design parameters
subroutine Report()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem data
    type(GeomType)  :: geom     ! Geometry data
    type(BoundType) :: bound    ! Boundary data
    type(MeshType)  :: mesh     ! Basepaire data
    type(DNAType)   :: dna      ! B-form DNA data

    character(10) :: char_sec, char_edge, char_cut, char_vert
    integer :: i, sec, edge_in, edge, max_stap, min_stap
    logical :: results

    sec       = 1           ! Section number
    edge_in   = 2           ! Edge length
    char_vert = "beveled"   ! Flat or beveled vertex
    char_cut  = "max"       ! Staple-break rule

    ! Open file
    open(unit = 90, file = "Report_2D_Flat_"//trim(char_cut)//".txt", form="formatted")

    ! Remove the directory and files
    results = SYSTEMQQ("rd "//trim("output")//' /s /q')

    ! Infomation
    write(90, "(a)"), "==========================================="
    write(90, "(a)"), "Sec 1: DX tile defined on honeycomb lattice"
    write(90, "(a)"), "Edge 1: 31bp, 2: 42bp , 3: 52bp, 4: 63bp"
    write(90, "(a)"), "Staple cutting: "//trim(adjustl(char_cut))
    write(90, "(a)"), "==========================================="
    write(90, "(a)")

    call space (90, 30)
    write(90, "(a)"), "|=========================================================================================================================================================================|"
    call space (90, 30)
    write(90, "(a)"), "|=========|==================== TOTAL LENGTH ==================|======== STAPLE =======|===== SEED ====|== STRAND RATIO =|= NU RATIO=|= PARA=|=== CROSSOVER ===|=UNPAIRED=|"
    call space (90, 30)
    write(90, "(a)"), " Sec| Edge| L_Scaf| L_Stap|   L_BP|L_Beveled(ratio)| MaxE| MinE| nStap| Min| Max|   ave| 14nt| S14| 4nt| 14nt|  S14|  4nt| 14nt|  4nt| p1| p2| scaf|  stap| one| scaf| stap"
    write(90, "(a$)"), "----------------------------  "
    write(90, "(a )"), "----|-----|-------|-------|-------|----------------|-----|-----|------|----|----|------|-----|----|----|-----|-----|-----|-----|-----|---|---|-----|------|----|-----|----|"

    min_stap = 10000
    max_stap =-10000

    ! Problem
    do i = 1, 24

        ! Edge length
        if(i == 5 .or. i == 10 .or. i == 15) then
            edge = edge_in - 1
        else if(i == 23) then
            edge = edge_in + 2
        else if(i == 24) then
            edge = edge_in + 4
        else
            edge = edge_in
        end if

        ! Initialize input
        call Input_Initialize_Report(prob, geom, mesh, i, sec, edge, char_vert, char_cut)

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

        ! Print problem
        write(90, "(a4, a24, a$)"), trim(adjustl(Int2Str(i)))//". ", trim(prob.name_prob), ","

        ! Print information
        write(90, "(a6$)"), trim(adjustl(Int2Str(sec)))//"|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(edge)))//"|"
        write(90, "(a8$)"), trim(adjustl(Int2Str(dna.n_base_scaf)))//"|"
        write(90, "(a8$)"), trim(adjustl(Int2Str(dna.n_base_stap)))//"|"
        write(90, "(a8$)"), trim(adjustl(Int2Str(mesh.n_node)))//"|"
        write(90, "(a8$)"), trim(adjustl(Int2Str(mesh.n_beveled)))//"("
        write(90, "(a9$)"), trim(adjustl(Dble2Str(dble(mesh.n_beveled)/dble(mesh.n_node)*100.0d0)))//"%)|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(geom.max_edge_length)))//"|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(geom.min_edge_length)))//"|"
        write(90, "(a7$)"), trim(adjustl(Int2Str(dna.n_stap)))//"|"
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.len_min_stap)))//"|"
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.len_max_stap)))//"|"
        write(90, "(a7$)"), trim(adjustl(Dble2Str2(dna.len_ave_stap)))//"|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_14nt)))//"|"
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_s14nt)))//"|"
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_4nt)))//"|"
        write(90, "(a6$)"), trim(adjustl(Dble2Str(dble(dna.n_14nt)/dble(dna.n_stap))))//"|"
        write(90, "(a6$)"), trim(adjustl(Dble2Str(dble(dna.n_s14nt)/dble(dna.n_stap))))//"|"
        write(90, "(a6$)"), trim(adjustl(Dble2Str(dble(dna.n_4nt)/dble(dna.n_stap))))//"|"
        write(90, "(a6$)"), trim(adjustl(Dble2Str(dble(dna.n_nt_14nt)/dble(dna.n_base_stap))))//"|"
        write(90, "(a6$)"), trim(adjustl(Dble2Str(dble(dna.n_nt_4nt)/dble(dna.n_base_stap))))//"|"
        write(90, "(a4$)"), trim(adjustl(Int2Str(prob.n_cng_min_stap)))//"|"
        write(90, "(a4$)"), trim(adjustl(Int2Str(prob.n_cng_max_stap)))//"|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_xover_scaf)))//"|"
        write(90, "(a7$)"), trim(adjustl(Int2Str(dna.n_xover_stap)))//"|"
        write(90, "(a5$)"), trim(adjustl(Int2Str(dna.n_sxover_stap)))//"|"
        write(90, "(a6$)"), trim(adjustl(Int2Str(dna.n_nt_unpaired_scaf)))//"|"
        write(90, "(a5 )"), trim(adjustl(Int2Str(dna.n_nt_unpaired_stap)))

        if(max_stap < dna.len_max_stap) max_stap = dna.len_max_stap
        if(min_stap > dna.len_min_stap) min_stap = dna.len_min_stap
    end do

    write(90, "(a)")
    write(90, "(a)"), "Minimum staple length   : "//trim(adjustl(Int2Str(min_stap)))
    write(90, "(a)"), "Maximum staple length   : "//trim(adjustl(Int2Str(max_stap)))

    ! Close file
    close(unit = 90)
end subroutine Report

! ---------------------------------------------------------------------------------------

! Print information of the PERDIX
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
        write(i, "(a)"), "* # of changing for min staple      : "//trim(adjustl(Int2Str(prob.n_cng_min_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* # of changing for max staple      : "//trim(adjustl(Int2Str(prob.n_cng_max_stap)))
        write(i, "(a)")

        ! ============================================================
        ! Base pair information
        ! ============================================================
        call Space(i, 11)
        write(i, "(a)"), "[ BASE PAIR ]"
        call Space(i, 16)
        write(i, "(a)"), "* # of basepairs            : "//trim(adjustl(Int2Str(mesh.n_node)))
        call Space(i, 16)
        write(i, "(a)"), "* # of beveled nucleotides  : "//&
            trim(adjustl(Int2Str(mesh.n_beveled)))//" ["//&
            trim(adjustl(Dble2Str(dble(mesh.n_beveled)/dble(mesh.n_node)*100.0d0)))//" %]"
        call Space(i, 16)
        write(i, "(a)"), "* Edge length [ min - max ] : ["//&
            trim(adjustl(Int2Str(geom.min_edge_length)))//" - "// &
            trim(adjustl(Int2Str(geom.max_edge_length)))//"]"
        write(i, "(a)")

        ! ============================================================
        ! Scaffold information
        ! ============================================================
        call Space(i, 11)
        write(i, "(a)"), "[ SCAFFOLD ]"
        call Space(i, 16)
        write(i, "(a)"), "* # of scaffold strands     : "//trim(adjustl(Int2Str(dna.n_scaf)))
        call Space(i, 16)
        write(i, "(a)"), "* # of total nucleotides    : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
        !call Space(i, 16)
        !write(i, "(a)"), "* # of unpaired regions     : "//trim(adjustl(Int2Str(dna.n_unpaired_scaf)))
        call Space(i, 16)
        write(i, "(a)"), "* # of unpaired nucleotides : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_scaf)))
        call Space(i, 16)
        write(i, "(a)"), "* # of double-crossovers    : "//trim(adjustl(Int2Str(dna.n_xover_scaf/2)))
        write(i, "(a)")

        ! ============================================================
        ! Staple information
        ! ============================================================
        call Space(i, 11)
        write(i, "(a)"), "[ STAPLE ]"
        call Space(i, 16)
        write(i, "(a)"), "* # of staples              : "//trim(adjustl(Int2Str(dna.n_stap)))
        call Space(i, 25)
        write(i, "(a)"), "@ with the 4nt dsDNA domain  - "//&
            trim(adjustl(Int2Str(dna.n_4nt)))//" ["//&
            trim(adjustl(Dble2Str(dble(dna.n_4nt)/dble(dna.n_stap)*100.0d0)))//" %]"
        call Space(i, 25)
        write(i, "(a)"), "@ with the 14nt dsDNA domain - "//&
            trim(adjustl(Int2Str(dna.n_14nt)))//" ["//&
            trim(adjustl(Dble2Str(dble(dna.n_14nt)/dble(dna.n_stap)*100.0d0)))//" %]"

        call Space(i, 16)
        write(i, "(a)"), "* # of nucleotides          : "//trim(adjustl(Int2Str(dna.n_base_stap)))
        call Space(i, 25)
        write(i, "(a)"), "@ in 4nt dsDNA domains       - "//&
            trim(adjustl(Int2Str(dna.n_nt_4nt)))//" ["//&
            trim(adjustl(Dble2Str(dble(dna.n_nt_4nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"
        call Space(i, 25)
        write(i, "(a)"), "@ in 14nt dsDNA domains      - "//&
            trim(adjustl(Int2Str(dna.n_nt_14nt)))//" ["//&
            trim(adjustl(Dble2Str(dble(dna.n_nt_14nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"

        !call Space(i, 16)
        !write(i, "(a)"), "* # of unpaired regions     : "//trim(adjustl(Int2Str(dna.n_unpaired_stap)))
        call Space(i, 16)
        write(i, "(a)"), "* # of unpaired nucleotides : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_stap)))
        call Space(i, 16)
        write(i, "(a)"), "* # of total crossovers     : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
        call Space(i, 16)
        write(i, "(a)"), "* # of single-crossovers    : "//trim(adjustl(Int2Str(dna.n_sxover_stap)))
        call Space(i, 16)
        write(i, "(a)"), "* # of double-crossovers    : "//trim(adjustl(Int2Str((dna.n_xover_stap-dna.n_sxover_stap)/2)))
        call Space(i, 16)
        write(i, "(a)"), "* Length [min - max- ave]   : ["//&
            trim(adjustl(Int2Str(dna.len_min_stap)))//" - "//&
            trim(adjustl(Int2Str(dna.len_max_stap)))//" - "//&
            trim(adjustl(Dble2Str(dna.len_ave_stap)))//"]"
        !call Space(i, 16)
        !write(i, "(a)"), "* # of total regions            : "//trim(adjustl(Int2Str(dna.n_tot_region)))
        !call Space(i, 16)
        !write(i, "(a)"), "* # of total 14nt dsDNA domains : "//trim(adjustl(Int2Str(dna.n_tot_14nt)))
        !call Space(i, 16)
        !write(i, "(a)"), "* # of total 4nt dsDNA domains  : "//trim(adjustl(Int2Str(dna.n_tot_4nt)))

        write(i, "(a)")
        write(i, "(a)"), "   --------------------------------------------------------------------------------"
        write(i, "(a)")
        write(i, "(a)"), " Completed"
    end do

    close(unit=11)
end subroutine Print_Information

! ---------------------------------------------------------------------------------------

! Deallocate global dynamic array
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
    do i = 1, bound.n_junc
        deallocate(bound.junc(i).iniL)
        deallocate(bound.junc(i).modP)
        deallocate(bound.junc(i).croP)
        deallocate(bound.junc(i).node)
        deallocate(bound.junc(i).conn)
        deallocate(bound.junc(i).type_conn)
    end do

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

! Verify solution
subroutine Verify_Solution(mesh, dna)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    double precision :: verify
    integer :: i

    verify = 0.0d0
    do i = 1, dna.n_top
        if( dna.top(i).xover == -1 .or. &
            dna.top(i).up    == -1 .or. &
            dna.top(i).dn    == -1 .or. &
            dna.top(i).node  == -1 ) then
            verify = verify + dna.top(i).pos(1)
            verify = verify + dna.top(i).pos(2)
            verify = verify + dna.top(i).pos(3)
            verify = verify + dble(dna.top(i).strand) + dble(dna.top(i).address)
         end if

        verify = verify + dble(dna.top(i).strand) + dble(dna.top(i).address)
    end do

    verify = verify + dble(dna.n_nt_14nt)/dble(dna.n_base_stap)*100.0d0
    verify = verify + dble(dna.n_nt_4nt) /dble(dna.n_base_stap)*100.0d0
    verify = verify + dble(dna.len_ave_stap + dna.n_14nt + dna.n_s14nt + dna.n_4nt + dna.n_only_4nt)
    verify = verify + dble(dna.n_nt_14nt + dna.n_nt_4nt + dna.n_tot_region + dna.n_tot_14nt + dna.n_tot_4nt)
    verify = verify + dble(dna.len_min_stap + dna.len_max_stap + dna.n_unpaired_scaf + dna.n_nt_unpaired_scaf)
    verify = verify + dble(dna.n_unpaired_stap + dna.n_nt_unpaired_stap)

    do i = 1, mesh.n_node
        verify = verify + mesh.node(i).pos(1)
        verify = verify + mesh.node(i).pos(2)
        verify = verify + mesh.node(i).pos(3)
    end do

    write(0, "(a)")
    write(0, "(a)"), "[ONLY DEBUG MODE]"
    write(0, "(a25, a)"), " 3.41299052999253E+07"," - Reference: 1 - 2"
    write(0, "(a25, a)"), " 3.88022588879539E+07"," - Reference: 2 - 2"
    write(0, "(a25, a)"), " 2.61079451760243E+07"," - Reference: 9 - 2"
    write(0, "(es25.14)"), verify
end subroutine Verify_Solution

! ---------------------------------------------------------------------------------------

end program PERDIX_OPEN