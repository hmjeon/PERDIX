!
! ---------------------------------------------------------------------------------------
!
!                                   Module for SeqDesign
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/11/11
!                                                            Last  modified  : 2016/11/14
!
! ---------------------------------------------------------------------------------------
!
module SeqDesign

    use Ifport

    use Data_Prob
    use Data_Mesh
    use Data_DNA

    use Section

    use Para
    use Mani
    use List
    use Math

    implicit none

    public  SeqDesign_Design

    private SeqDesign_Build_dnaTop
    private SeqDesign_Get_Rand_Sequence
    private SeqDesign_Get_Comp_Sequence
    private SeqDesign_Build_Strand
    private SeqDesign_Make_Noncir_Stap_Single_Xover
    private SeqDesign_Make_Noncir_Stap_Nick
    private SeqDesign_Make_Nick_Scaf
    private SeqDesign_Make_Short_Scaf
    private SeqDesign_Build_Sequence_Design_Max
    private SeqDesign_Build_Sequence_Design_Opt
    private SeqDesign_Build_Sequence_Design_Mix
    private SeqDesign_Build_Sequence_Design
    private SeqDesign_Avoid_Barrier
    private SeqDesign_Make_Short_Strand
    private SeqDesign_Count_Remainder
    private SeqDesign_Rebuild_Strand
    private SeqDesign_Build_Region_Staple
    private SeqDesign_Order_Staple
    private SeqDesign_Print_14nt_Region
    private SeqDesign_CirGraph_Count_Edge
    private SeqDesign_CirGraph_Init_Variable
    private SeqDesign_Assign_Sequence
    private SeqDesign_Set_M13mp18
    private SeqDesign_Get_M13mp18
    private SeqDesign_Import_Sequence
    private SeqDesign_Set_Rand_Sequence
    private SeqDesign_Write_Strand
    private SeqDesign_Write_Graphical_Output
    private SeqDesign_Chimera_Atom
    private SeqDesign_Chimera_Route
    private SeqDesign_Chimera_Sequence_Design
    private SeqDesign_Chimera_Strand

    type :: RegionType
        integer :: types,    length                 ! 1-vertex, 2-edge, region length
        integer :: sta_pos,  cen_pos,  end_pos      ! Position
        integer :: sta_base, cen_base, end_base     ! Base ID
        integer :: n_14nt,   n_4nt                  ! # of 14nt seeds
    end type RegionType

    ! GraphType data strucutre is corresponding to node
    type GraphType
        ! node2base(i) : From base to index
        ! base2node(i) : From index to base
        ! node(i)      : 0-nomarl, 1-first 14nt, 2-secondary 14nt, 3-4nt
        ! edge(i,1)    : [1]->2 : source of connectivity
        ! edge(i,2)    : 1->[2] : target of connectivity
        ! edge(i,3)    : Weight factor
        integer, allocatable :: node2base(:)
        integer, allocatable :: base2node(:)
        integer, allocatable :: node(:)
        integer, allocatable :: edge(:,:)
    end type GraphType
contains

! ---------------------------------------------------------------------------------------

! Design topology
! Last updated on Saturday 17 September 2016 by Hyungmin
subroutine SeqDesign_Design(prob, geom, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom
    type(MeshType), intent(inout) :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Build dnaTop data from dna base data
    call SeqDesign_Build_dnaTop(dna)

    ! Build strand data from dnaTop
    call SeqDesign_Build_Strand(dna)

    ! Make non-circular staple strand
    call SeqDesign_Make_Noncir_Stap_Nick(mesh, dna)



    ! Build sequence design with non-circular staple strands
    ! ++++++++++++++++++++++++++
    ! max에 문제가 발생됨
    ! ++++++++++++++++++++++++++
    if(para_cut_stap_method == "max") call SeqDesign_Build_Sequence_Design_Max(prob, mesh, dna)
    if(para_cut_stap_method == "opt") call SeqDesign_Build_Sequence_Design_Opt(prob, mesh, dna)

    ! ==================================================
    !
    ! Rebuild strand data from dnaTop
    !call SeqDesign_Rebuild_Strand(dna)

    ! Chimera sequence design
    !call SeqDesign_Chimera_Sequence_Design(prob, geom, mesh, dna)
    !stop
    !
    ! ==================================================

    ! Make nick in scaffold strand
    call SeqDesign_Make_Nick_Scaf(geom, mesh, dna)

    ! Make short scaffold strand
    call SeqDesign_Make_Short_Scaf(mesh, dna)

    ! Rebuild strand data from dnaTop
    call SeqDesign_Rebuild_Strand(dna)

    ! List in long length order of the staple
    call SeqDesign_Order_Staple(dna)

    ! Print 14nt region with various representations
    if(para_max_cut_scaf == 0) call SeqDesign_Print_14nt_Region(prob, geom, mesh, dna)

    ! Assign DNA sequence according to para_set_seq_scaf
    call SeqDesign_Assign_Sequence(dna)

    ! Write strand data
    call SeqDesign_Write_Strand(prob, geom, mesh, dna)

    ! Write atom model by dnaTop and strand data
    call SeqDesign_Chimera_Atom(prob, dna)

    ! Chimera topology route
    call SeqDesign_Chimera_Route(prob, mesh, dna)

    ! Chimera sequence design
    call SeqDesign_Chimera_Sequence_Design(prob, geom, mesh, dna)

    ! Write Chimera file for strand and sequence
    call SeqDesign_Chimera_Strand(prob, dna)
end subroutine SeqDesign_Design

! ---------------------------------------------------------------------------------------

! Reset possible staple crossovers
! Last updated on Sunday 21 August 2016 by Hyungmin
subroutine SeqDesign_Reset_Possible_Stap_Xover(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Crossover based on cross-sectional edges
    type :: CroLType
        integer :: max_bp, min_bp
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, k, croL_cur, croL_com, sec_cur, sec_com, id_bp
    integer :: up_scaf1, dn_scaf1, up_scaf2, dn_scaf2
    integer :: up_cur, up_com, dn_cur, dn_com
    logical :: b_nei_up, b_nei_dn, b_scaf

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Reset possible staple crossovers"
    end do

    ! Reset connectivity for bases in staple at crossovers
    dna.n_xover_stap  = 0
    dna.n_sxover_stap = 0
    do i = 1, mesh.n_node

        ! Reset connection at crossovers
        if(dna.base_stap(i).xover /= -1) then
            if(dna.base_scaf(i).up == -1) then
                dna.base_stap(i).up = dna.base_scaf(i).dn
                dna.base_stap(i).dn = -1
            else if(dna.base_scaf(i).dn == -1) then
                dna.base_stap(i).up = -1
                dna.base_stap(i).dn = dna.base_scaf(i).up
            else
                dna.base_stap(i).up = dna.base_scaf(i).dn
                dna.base_stap(i).dn = dna.base_scaf(i).up
            end if
        end if

        dna.base_stap(i).xover = -1
    end do

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    croL(1:geom.n_croL).max_bp = -999999
    croL(1:geom.n_croL).min_bp =  999999

    ! Find maximum and minimum basepair ID in cross-sectional edges
    do i = 1, mesh.n_node
        croL_cur = mesh.node(i).croL
        id_bp    = mesh.node(i).bp

        ! Set maximum and minimum base ID
        if(croL(croL_cur).max_bp < id_bp) croL(croL_cur).max_bp = id_bp
        if(croL(croL_cur).min_bp > id_bp) croL(croL_cur).min_bp = id_bp
    end do

    ! Find the possible staple double crossovers
    dna.n_xover_stap  = 0
    dna.n_sxover_stap = 0
    do i = 1, mesh.n_node       ! Loop for current node

        ! Print progress bar
        call Mani_Progress_Bar(i, mesh.n_node)

        ! Loop for comparing node
        do j = i + 1, mesh.n_node

            ! Exception for the pre-constructed crossovers (due to double crossover)
            if(dna.base_stap(i).xover /= -1 .and. dna.base_stap(j).xover /= -1) then
                cycle
            end if

            ! It should be skipped when condition below
            ! Basepair ID and iniL shoud be the same and croL and section ID should be different
            if(mesh.node(i).bp   /= mesh.node(j).bp  ) cycle
            if(mesh.node(i).iniL /= mesh.node(j).iniL) cycle
            if(mesh.node(i).croL == mesh.node(j).croL) cycle
            if(mesh.node(i).sec  == mesh.node(j).sec ) cycle

            ! Find section ID
            sec_cur  = mesh.node(i).sec
            sec_com  = mesh.node(j).sec
            croL_cur = mesh.node(i).croL
            croL_com = mesh.node(j).croL
            id_bp    = mesh.node(i).bp

            ! To eliminate boundary staple crossovers
            if(croL(croL_cur).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_cur).max_bp - para_gap_xover_bound_stap < id_bp) cycle
            if(croL(croL_com).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_com).max_bp - para_gap_xover_bound_stap < id_bp) cycle

            ! Determine whether the node has crossover or not
            if(Section_Connection_Stap(geom, sec_cur, sec_com, id_bp) == .true.) then

                ! To eliminate crossover if there is neighboring scaffold crossover
                !
                !     dn_scaf1   i   up_scaf1
                !    *---*---*---*---*---*---*-->  : node i
                !            |       |
                ! <--*---*---*---*---*---*---*     : node j
                !     up_scaf2   j   dn_scaf2
                b_scaf   = .false.
                up_scaf1 = mesh.node(i).id
                dn_scaf1 = mesh.node(i).id
                up_scaf2 = mesh.node(j).id
                dn_scaf2 = mesh.node(j).id

                ! Check neighbor scaffold crossovers
                do k = 1, para_gap_xover_two
                    if( (dna.base_scaf(mesh.node(up_scaf1).id).xover == dna.base_scaf(mesh.node(dn_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(dn_scaf2).id).xover == dna.base_scaf(mesh.node(up_scaf1).id).id ) .or.  &
                        (dna.base_scaf(mesh.node(dn_scaf1).id).xover == dna.base_scaf(mesh.node(up_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(up_scaf2).id).xover == dna.base_scaf(mesh.node(dn_scaf1).id).id ) ) then
                        b_scaf = .true.
                        exit
                    else
                        up_scaf1 = mesh.node(up_scaf1).up
                        dn_scaf1 = mesh.node(dn_scaf1).dn
                        up_scaf2 = mesh.node(up_scaf2).up
                        dn_scaf2 = mesh.node(dn_scaf2).dn
                    end if
                end do

                if(b_scaf == .true.) cycle

                ! Find upper or downward neighboring crossovers
                ! Node numbering is opposite to staple ID
                up_cur   = mesh.node(i).dn
                up_com   = mesh.node(j).up
                id_bp    = mesh.node(up_cur).bp
                b_nei_up = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                dn_cur   = mesh.node(i).up
                dn_com   = mesh.node(j).dn
                id_bp    = mesh.node(dn_cur).bp
                b_nei_dn = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                ! Set current and previous or next crossovers (double crossover)
                ! For current crossover
                dna.n_xover_stap       = dna.n_xover_stap + 2
                dna.base_stap(i).xover = dna.base_stap(j).id
                dna.base_stap(j).xover = dna.base_stap(i).id

                ! For neighboring crossover
                if(b_nei_up == .true.) then
                    dna.base_stap(up_cur).xover = dna.base_stap(up_com).id
                    dna.base_stap(up_com).xover = dna.base_stap(up_cur).id

                    ! Set connectivity
                    dna.base_stap(i).up      = dna.base_stap(j).id
                    dna.base_stap(j).dn      = dna.base_stap(i).id
                    dna.base_stap(up_cur).dn = dna.base_stap(up_com).id
                    dna.base_stap(up_com).up = dna.base_stap(up_cur).id

                else if(b_nei_dn == .true.) then
                    dna.base_stap(dn_cur).xover = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).xover = dna.base_stap(dn_cur).id

                    ! Set connectivity
                    dna.base_stap(i).dn      = dna.base_stap(j).id
                    dna.base_stap(j).up      = dna.base_stap(i).id
                    dna.base_stap(dn_cur).up = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).dn = dna.base_stap(dn_cur).id
                end if
            end if
        end do
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The number of possible staple crossovers    : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
        write(i, "(a)")
    end do

    ! Deallocate memory
    deallocate(croL)
end subroutine SeqDesign_Reset_Possible_Stap_Xover

! ---------------------------------------------------------------------------------------

! Build dnaTop data from dna base data
! Last updated on Thursday 12 May 2016 by Hyungmin
subroutine SeqDesign_Build_dnaTop(dna)
    type(DNAType), intent(inout) :: dna

    integer   :: i, j, n_jump, id_up, id_down
    character :: seq_across

    ! Set the number of dna top
    dna.n_top = dna.n_base_scaf + dna.n_base_stap

    ! Allocate dna top data
    allocate(dna.top(dna.n_top))

    ! Copy data from bases in scaffold strand
    do i = 1, dna.n_base_scaf

        ! ID and node ID
        dna.top(i).id   = dna.base_scaf(i).id
        dna.top(i).node = dna.base_scaf(i).node

        ! Connectivity
        dna.top(i).up = dna.base_scaf(i).up
        dna.top(i).dn = dna.base_scaf(i).dn

        ! Crossover and across ID
        dna.top(i).xover = dna.base_scaf(i).xover

        if(dna.base_scaf(i).across == -1) then
            dna.top(i).across = -1
        else
            dna.top(i).across = dna.base_scaf(i).across + dna.n_base_scaf
        end if

        ! Strand, residue ID and sequence
        dna.top(i).strand  = -1
        dna.top(i).address = -1
        dna.top(i).b_14nt  = .false.
        dna.top(i).seq     = "N"

        ! Assign position vector
        dna.top(i).pos(1:3) = dna.base_scaf(i).pos(1:3)
    end do

    ! Copy data from bases in staple strand
    n_jump = dna.n_base_scaf
    do i = 1, dna.n_base_stap

        ! ID and node ID
        dna.top(i+n_jump).id   = dna.base_stap(i).id + n_jump
        dna.top(i+n_jump).node = dna.base_stap(i).node

        ! Connectivity
        if(dna.base_stap(i).up == -1) then
            dna.top(i+n_jump).up = -1
        else
            dna.top(i+n_jump).up = dna.base_stap(i).up + n_jump
        end if

        if(dna.base_stap(i).dn == -1) then
            dna.top(i+n_jump).dn = -1
        else
            dna.top(i+n_jump).dn = dna.base_stap(i).dn + n_jump
        end if

        ! Crossover and across ID
        if(dna.base_stap(i).xover == -1) then
            dna.top(i+n_jump).xover = -1
        else
            dna.top(i+n_jump).xover = dna.base_stap(i).xover + n_jump
        end if
        dna.top(i+n_jump).across = dna.base_stap(i).across

        ! Strand, residue ID and sequence
        dna.top(i+n_jump).strand  = -1
        dna.top(i+n_jump).address = -1
        dna.top(i+n_jump).b_14nt  = .false.
        dna.top(i+n_jump).seq     = "N"

        ! Sequence for poly Tn loop
        if(dna.base_stap(i).across == -1) dna.top(i+n_jump).seq = "T"

        ! Assign position vector
        dna.top(i+n_jump).pos(1:3) = dna.base_stap(i).pos(1:3)
    end do

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |                       6. Sequence design                           |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
        call Space(i, 6)
        write(i, "(a)"), "6.1. Build DNA nucleotide based data"
        call Space(i, 11)
        write(i, "(a)"), "* The total number of bases                : "//trim(adjustl(Int2Str(dna.n_top)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of bases in scaffold strand   : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of bases in staple strand     : "//trim(adjustl(Int2Str(dna.n_base_stap)))
    end do

    write(0, "(a)"); write(11, "(a)")
end subroutine SeqDesign_Build_dnaTop

! ---------------------------------------------------------------------------------------

! Get random sequence
! Last updated on Sunday 28 Feb 2016 by Hyungmin
function SeqDesign_Get_Rand_Sequence result(seq)
    character :: seq
    integer :: random

    ! Random number generation (range from 1 to 4)
    random = int(4*rand()) + 1

    if(random == 1) then
        seq = "A"
    else if(random == 2) then
        seq = "T"
    else if(random == 3) then
        seq = "G"
    else if(random == 4) then
        seq = "C"
    end if
end function SeqDesign_Get_Rand_Sequence

! ---------------------------------------------------------------------------------------

! Get complementary sequence
! Last updated on Sunday 27 Mar 2016 by Hyungmin
function SeqDesign_Get_Comp_Sequence(seq) result(com_seq)
    character, intent(in) :: seq

    character :: com_seq

    if(seq == "A") then
        com_seq = "T"
    else if(seq == "T") then
        com_seq = "A"
    else if(seq == "C") then
        com_seq = "G"
    else if(seq == "G") then
        com_seq = "C"
    else
        write(0, "(a$)"), "Error - Not defined sequnce : "
        write(0, "(a )"), "SeqDesign_Get_Comp_Sequence"
        stop
    end if
end function SeqDesign_Get_Comp_Sequence

! ---------------------------------------------------------------------------------------

! Build strand data from dnaTop
! Last updated on Saturday 26 Feb 2016 by Hyungmin
subroutine SeqDesign_Build_Strand(dna)
    type(DNAType), intent(inout) :: dna

    type(StrandType), allocatable, dimension(:) :: strand
    logical,          allocatable, dimension(:) :: b_visit
    type(ListBase),   pointer :: list_base, ptr_base
    type(TopType) :: cur_base

    integer, parameter :: max_strand = 2000
    integer :: i, j, n_strand, n_base, int_base
    logical :: b_end

    ! Allocate and initialize strand data
    allocate(strand(max_strand))
    call Mani_Init_StrandType(strand, max_strand)

    ! Nullify the linked lists
    nullify(list_base, ptr_base)

    ! Allocate and initilize the b_visit data
    allocate(b_visit(dna.n_top))
    b_visit(1:dna.n_top) = .false.

    n_strand = 0
    do
        ! Check b_visit if there is 0 (0 means not visiting base)
        b_end = .true.
        do i = 1, dna.n_top
            if(b_visit(dna.top(i).id) == .false.) then
                b_end = .false.
                exit
            end if
        end do
        if(b_end == .true.) exit

        ! Increase the number of strands
        n_strand = n_strand + 1
        cur_base = dna.top(i)
        int_base = cur_base.id

        if(n_strand > max_strand) then
            write(0, "(a$)"), "Error - Excess of maximum number of strands : "
            write(0, "(a )"), "SeqDesign_Build_Strand"
            write(0, "(a$)"), "The current number of strands : "
            write(0, "(i7)"), n_strand
            stop
        end if

        ! Find the first base in the current strand
        do
            if(cur_base.dn == -1) then
                ! cur_base is at the 3'-end of the strand
                strand(n_strand).b_circular = .false.
                exit
            else if(cur_base.dn == int_base) then
                ! cur_base goes back to the starting point
                strand(n_strand).b_circular = .true.
                cur_base = dna.top(int_base)
                exit
            end if

            ! Update current base
            cur_base = dna.top(cur_base.dn)

            if(b_visit(cur_base.id) == .true.) then
                write(0, "(a$)"), "Error - Reached a visited base : "
                write(0, "(a )"), "SeqDesign_Build_Strand"
                stop
            end if
        end do

        ! Walk through the current strand
        n_base = 1

        ! Insert data to linked list
        allocate(ptr_base)
        ptr_base%id = cur_base.id
        list_base => List_Insert_Base(list_base, ptr_base)

        dna.top(cur_base.id).strand  = n_strand
        dna.top(cur_base.id).address = n_base
        dna.top(cur_base.id).b_14nt  = .false.
        b_visit(cur_base.id)         = .true.

        ! Loop to add a new base into current strand
        do
            ! Check for going out loop
            if(strand(n_strand).b_circular == .false.) then
                ! If non-circular and upper ID is equal to -1
                if(cur_base.up == -1) exit
            else
                ! If circular and base ID is equal to init base ID
                if(cur_base.up == int_base) exit
            end if

            cur_base = dna.top(cur_base.up)

            if(b_visit(cur_base.id) == .true.) then
                write(0, "(a$)"), "Error - Reached a visited base : "
                write(0, "(a )"), "SeqDesign_Build_Strand"
                stop
            end if

            n_base = n_base + 1

            ! Insert data to linked list
            allocate(ptr_base)
            ptr_base%id = cur_base.id
            list_base => List_Insert_Base(list_base, ptr_base)

            dna.top(cur_base.id).strand  = n_strand
            dna.top(cur_base.id).address = n_base
            dna.top(cur_base.id).b_14nt  = .false.
            b_visit(cur_base.id)         = .true.
        end do

        strand(n_strand).n_base = n_base

        ! Allocate array using the linked list
        allocate(strand(n_strand).base(n_base))

        ! Put strand data from linked list
        ptr_base => list_base
        do i = 1, n_base
            strand(n_strand).base(n_base+1-i) = ptr_base%id
            ptr_base => ptr_base%next
        end do
    end do

    ! Deallocate b_visit data
    deallocate(b_visit)

    ! Set the number of strands and copy data
    dna.n_strand = n_strand
    dna.n_scaf   = 1
    dna.n_stap   = dna.n_strand - 1

    allocate(dna.strand(dna.n_strand))
    do i = 1, dna.n_strand

        dna.strand(i).n_base     = strand(i).n_base
        dna.strand(i).b_circular = strand(i).b_circular

        ! Copy data from base to dna.base
        allocate(dna.strand(i).base(dna.strand(i).n_base))
        do j = 1, dna.strand(i).n_base
            dna.strand(i).base(j) = strand(i).base(j)
        end do

        ! Set strand type
        dna.strand(i).types = "stap"
        if(i == 1) dna.strand(i).types = "scaf"

        ! Deallocate strand base data
        deallocate(strand(i).base)
    end do

    ! Deallocate strand 
    deallocate(strand)

    ! Delete linked list allocated
    call List_Delete_Base(list_base)
    !call List_Delete_Base(ptr_base)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.2. Build strand based data"
        call Space(i, 11)
        write(i, "(a)"), "* The total number of strands              : "//trim(adjustl(Int2Str(dna.n_strand)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of scaffold strands           : "//trim(adjustl(Int2Str(1)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of staple strands             : "//trim(adjustl(Int2Str(dna.n_strand - 1)))
        call Space(i, 11)
        write(i, "(a)"), "* Detailed DNA strand data builded"
    end do

    ! Print progress in detail
    do i = 1, dna.n_strand
        write(11, "(i20, a$)"), i, " th strand, type : "
        if(dna.strand(i).types == "scaf") then
            write(11, "(a$)"), "scaffold, "
        else if(dna.strand(i).types == "stap") then
            write(11, "(a$)"), "staple, "
        end if
        write(11, "(a, i7$ )"), "# of bases : ", dna.strand(i).n_base
        write(11, "(a, l, a)"), ", circular : ", dna.strand(i).b_circular

        ! Print bases numbering
        write(11, "(a27, i7, a$)"), "Bases : ", dna.strand(i).base(1), " ->"
        do j = 2, dna.strand(i).n_base - 1
            if(mod(j, 100)== 0) then
                write(11, "(i7, a$ )"), dna.strand(i).base(j), " ->"
            else if(mod(j, 100)== 1) then
                write(11, "(i34, a$)"), dna.strand(i).base(j), " ->"
            else
                write(11, "(i7, a$ )"), dna.strand(i).base(j), " ->"
            end if
        end do
        write(11, "(i7)"), dna.strand(i).base(dna.strand(i).n_base)
        write(11, "(a )")
    end do

    write(0, "(a)"); write(11, "(a)")
end subroutine SeqDesign_Build_Strand

! ---------------------------------------------------------------------------------------

! Make non-circular staple strand by single crossover
! Last updated on Saturday 17 September 2016 by Hyungmin
subroutine SeqDesign_Make_Noncir_Stap_Single_Xover(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, base, init_base, across, node, dn_base

    ! Loop to make non-circular staple strand
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! Find non-circular strand
        if(dna.strand(i).b_circular == .true.) then

            ! Find first base that has double crossover
            base = dna.strand(i).base(1)
            do j = 1, dna.strand(i).n_base
                if( dna.top(base).xover /= -1 .and. &
                    dna.top(dna.top(base).dn).xover /= -1 ) then
                    exit
                end if
                base = dna.top(base).dn
            end do

            ! Find downward base
            dn_base = dna.top(base).dn

            ! Disconnect between current and downward bases
            dna.top(base).dn    = -1
            dna.top(dn_base).up = -1

            ! If this point is crossover
            if( dna.top(base).xover    == dna.top(dn_base).id .and. &
                dna.top(dn_base).xover == dna.top(base).id ) then

                ! Delete single crossover
                dna.n_xover_stap  = dna.n_xover_stap  - 1
                dna.n_sxover_stap = dna.n_sxover_stap + 1

                dna.top(base).xover    = -1
                dna.top(dn_base).xover = -1
            end if
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.3. Make non-circular strand by single crossover"
        call Space(i, 11)
        write(i, "(a)"), "* The number of staple strands             : "//trim(adjustl(Int2Str(dna.n_base_stap)))
        write(i, "(a)")
    end do
end subroutine SeqDesign_Make_Noncir_Stap_Single_Xover

! ---------------------------------------------------------------------------------------

! Make non-circular staple strand
! Last updated on Wednesday 09 August 2016 by Hyungmin
subroutine SeqDesign_Make_Noncir_Stap_Nick(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, base, init_base, across, node, dn_base
    integer :: count, max_count, start_base, max_start_base

    ! Loop to make non-circular staples
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! Find non-circular strand
        if(dna.strand(i).b_circular == .true.) then

            ! Find first base position in which the double crossover exists
            base = dna.strand(i).base(1)
            do j = 1, dna.strand(i).n_base
                if( dna.top(base).xover /= -1 .and. &
                    dna.top(dna.top(base).dn).xover /= -1 ) then
                    exit
                end if
                base = dna.top(base).dn
            end do

            ! In case of the small or vertex dsDNA region
            init_base = base

            ! Find maximum length region
            count      = 1
            max_count  = 0
            start_base = base
            do j = 1, dna.strand(i).n_base

                base   = dna.top(base).id
                across = dna.top(base).across

                if(across == -1) then
                    ! For base in Tn-loop
                    if(max_count < count) then
                        max_count      = count
                        max_start_base = start_base
                    end if
                    count      = 0
                    start_base = base
                else
                    ! Check scaffold and staple crossovers
                    if(dna.top(base).xover /= -1 .or. dna.top(across).xover /= -1) then
                        if(max_count < count) then
                            max_count      = count
                            max_start_base = start_base
                        end if
                        count      = 0
                        start_base = base
                    end if
                end if

                count = count + 1
                base  = dna.top(base).up
            end do

            ! Exclude starting and ending bases
            max_count = max_count - 1

            ! The minimum gap between xover/Tn and first nick, para_gap_xover_nick1(default is 10)
            if( max_count > para_gap_xover_nick1 - 2 ) then

                ! Set base at the center of the maximum length region
                base = dna.top(max_start_base).id
                do j = 1, max_count / 2 + 1
                    base = dna.top(base).up
                end do
            else
                ! For Exception in case of small maximum region
                base = init_base
            end if

            dn_base = dna.top(base).dn

            ! Disconnect between current and downward bases
            dna.top(base).dn    = -1
            dna.top(dn_base).up = -1

            ! If this is the crossover position
            ! In case of "max_count > para_gap_xover_nick1 - 2", there are bug in making single xover
            if( dna.top(base).xover    == dna.top(dn_base).id .and. &
                dna.top(dn_base).xover == dna.top(base).id ) then

                ! Delete single crossover
                dna.n_xover_stap  = dna.n_xover_stap  - 1
                dna.n_sxover_stap = dna.n_sxover_stap + 1

                dna.top(base).xover    = -1
                dna.top(dn_base).xover = -1
            end if
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.3. Make non-circular strand"
        call Space(i, 11)
        write(i, "(a)"), "* The number of staple strands             : "//trim(adjustl(Int2Str(dna.n_base_stap)))
        write(i, "(a)")
    end do
end subroutine SeqDesign_Make_Noncir_Stap_Nick

! ---------------------------------------------------------------------------------------

! Build sequence design combined maximum and optimal cutting
! Last updated on Thursday 10 November 2016 by Hyungmin
subroutine SeqDesign_Build_Sequence_Design_Mix(prob, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, cen_base, pre_base, up_base, pre_reg
    integer :: n_region, cn_tn, length, final_length, cen_pos, pre_pos, bgn_pos
    integer :: cng_para, up, xover, upxover
    integer :: cng_gap_nick
    logical :: b_ext, b_cut, b_14nt

    cng_gap_nick = 8

    ! Loop to make short staple based on mixed method
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! For short staple, para_max_cut_stap <= 60
        if(dna.strand(i).n_base <= para_max_cut_stap) cycle

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! For vertex in case of DX tile design
        !if(cn_tn > 0 .and. prob.sel_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! ==================================================
        !
        ! Build region of staple strands
        !
        ! ==================================================
        ! Tn loop |<----->| : 7 nt poly T loop
        !     ~~~~~=======*===========*=====*~~~~~~~
        !                 |<--------->|<--->|
        !           11 and 5 nt region not including crossover and nick
        !
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Print information on region of staple strands
        do k = 0, 11, 11
            write(k, "(i$    )"), i
            write(k, "(a, i5$)"), " - th strand, # of total bases : ", dna.strand(i).n_base
            write(k, "(a, i5 )"), ", # of bases in Tn : ", cn_tn

            do j = 1, n_region
                write(k, "(i20, a$)"), j, "  -  th region"
                write(k, "(a,  i3$)"), ", type : ",          region(j).types
                write(k, "(a,  i4$)"), ", region length : ", region(j).length
                write(k, "(a,  i4$)"), ", start pos : ",     region(j).sta_pos
                write(k, "(a,  i4$)"), ", center pos : ",    region(j).cen_pos
                write(k, "(a,  i4 )"), ", end pos : ",       region(j).end_pos
            end do
            write(k, "(a)"); call Space(k, 19)
            write(k, "(a)"), "----- Make nick position -----------------------------------------------------------------------------"
        end do

        ! ==================================================
        !
        ! Cut staples to make short multi staples with 14nt seeds
        !
        ! ==================================================
        bgn_pos = 0
        pre_reg = 0
        b_cut   = .false.
        b_14nt  = .false.
        j       = 0

        do
            ! Check loop
            j = j + 1
            if(j == n_region + 1) exit

            ! Skip when the region length is smaller than 5
            if(j /= n_region .and. region(j).length < 5) cycle

            ! Cutted length from center to beginning
            length = region(j).cen_pos - bgn_pos
            if(j == n_region) length = dna.strand(i).n_base - bgn_pos + 1

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_cut_stap) then
                b_cut = .true.
            end if

            ! Skip this current loop to prevent cutting this region
            if(j /= n_region .and. length <= para_max_cut_stap) then
                if( (b_14nt == .false. .and. region(j).length + 1 >= 14 .and. region(j).types == 1) .or. &
                    (b_14nt == .false. .and. region(j).length + 2 >= 14 .and. region(j).types == 2) ) then
                    b_14nt = .true.
                    cycle
                end if
            end if

            if(b_cut == .true. .and. b_14nt == .true. .and. length <= para_max_cut_stap) then
                !
                ! ==================================================
                !
                ! Optimal cutting to have 14nt seeds
                !
                ! ==================================================
                !
                ! Centered base ID and position
                cen_base = region(j).cen_base
                cen_pos  = region(j).cen_pos

                ! To avoid short region cutting
                if(region(j).length < cng_gap_nick) cycle

                ! If the final length is smaller than minimum length
                if(dna.strand(i).n_base - cen_pos < para_min_cut_stap) cycle

                ! Print progress
                do k = 0, 11, 11
                    write(k, "(i20, a$)"), j, " -> cut region"
                    write(k, "(a,  i4$)"), ", cutted length : ",   cen_pos - bgn_pos
                    write(k, "(a,  i4$)"), ", cutted pos : ",      cen_pos
                    write(k, "(a,  i4$)"), ", remained length : ", dna.strand(i).n_base - cen_pos
                    write(k, "(a      )"), " --> 14nt cutting"
                end do

                ! Increase the number of staples
                dna.n_stap = dna.n_stap + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(cen_base).up
                dna.top(cen_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update beginning position and flags
                bgn_pos = cen_pos
                b_cut   = .false.
                b_14nt  = .false.
            else if(b_cut == .true. .and. length > para_max_cut_stap) then
                !
                ! ==================================================
                !
                ! If the length exceeds the maximum cutting length
                !
                ! ==================================================
                !
                ! To avoid small region in previous region
                jj       = j
                b_ext    = .false.
                cng_para = para_gap_xover_nick

                ! Find previous region
                do
                    ! If the region length is longer than paramter value (default : 8)
                    if( (region(jj).types == 1 .and. region(jj).length >= cng_para * 2 + 2 + 1) .or. &
                        (region(jj).types == 2 .and. region(jj).length >= cng_para * 2 + 2) ) then
                        pre_base = region(jj).cen_base
                        pre_pos  = region(jj).cen_pos

                        ! Check cutted and remained length
                        if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                            pre_pos - bgn_pos <= para_max_cut_stap .and. &
                            dna.strand(i).n_base - pre_pos >= para_min_cut_stap ) exit
                    end if

                    ! Go back previous region
                    jj = jj - 1

                    ! Exception that the region cutted already
                    if(jj == pre_reg .or. jj == 0) then

                        ! Reset parameter and region index
                        jj       = j
                        cng_para = cng_para - 1

                        if(cng_para == 0) then
                            b_ext    = .true.
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos
                            exit
                        end if

                        if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                        if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                        if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                        do k = 0, 11, 11
                            call space(k, 16)
                            write(k, "(a)"), "** Adjusting parameter, para_gap_xover_nick : From "//&
                                trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                                trim(adjustl(Int2Str(cng_para)))//" - "//&
                                trim(adjustl(Int2Str(cng_para * 2 + 2)))
                        end do
                    end if
                end do

                ! No cutting due to exception, which may exceeds 60
                if(b_ext == .true.) then
                    do k = 0, 11, 11
                        call space(k, 19)
                        write(k, "(a$)"), "|-last region with exception, remaining length : "
                        write(k, "(i4)"), dna.strand(i).n_base - bgn_pos
                    end do
                    if(j == n_region) cycle
                end if

                ! If the cutted length is smaller than para_min_cut_stap
                if(pre_pos - bgn_pos < para_min_cut_stap) then
                    cycle
                end if

                ! Check 14nt seed cutting
                up      = dna.top(region(jj).end_base).up
                xover   = dna.top(up).xover
                upxover = dna.top(dna.top(up).up).xover

                if( region(jj).end_pos + 1 - bgn_pos >= para_min_cut_stap .and. &
                    region(jj).end_pos + 1 - bgn_pos <= para_max_cut_stap .and. &
                    dna.strand(i).n_base - region(jj).end_pos + 1 >= para_min_cut_stap .and. &
                    region(jj).length >= 12 .and. xover /= -1 .and. upxover /= -1 .and. para_set_stap_sxover == "on") then

                    ! Print progress
                    do k = 0, 11, 11
                        write(k, "(i20, a$)"), jj, " -> cut region"
                        write(k, "(a, i4$ )"), ", cutted length : ",   region(jj).end_pos + 1 - bgn_pos
                        write(k, "(a, i4$ )"), ", cutted pos : ",      region(jj).end_pos + 1
                        write(k, "(a, i4$ )"), ", remained length : ", dna.strand(i).n_base - region(jj).end_pos + 1
                        write(k, "(a      )"), " --> max cutting with single xover"
                    end do

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    if(dna.top(up).up == xover) dna.top(up).up = -1
                    if(dna.top(up).dn == xover) dna.top(up).dn = -1

                    if(dna.top(xover).up == up) dna.top(xover).up = -1
                    if(dna.top(xover).dn == up) dna.top(xover).dn = -1

                    dna.top(up).xover    = -1
                    dna.top(xover).xover = -1

                    dna.n_xover_stap  = dna.n_xover_stap  - 1
                    dna.n_sxover_stap = dna.n_sxover_stap + 1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = region(jj).end_pos + 1
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                else

                    ! Print progress
                    do k = 0, 11, 11
                        write(k, "(i20, a$)"), jj, " -> cut region"
                        write(k, "(a,  i4$)"), ", cutted length : ",   pre_pos - bgn_pos
                        write(k, "(a,  i4$)"), ", cutted pos : ",      pre_pos
                        write(k, "(a,  i4$)"), ", remained length : ", dna.strand(i).n_base - pre_pos
                        write(k, "(a      )"), " --> max cutting with nick"
                    end do

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = pre_pos
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                end if
            end if

            ! If remained length is smaller than para_max_cut_stap, exit this loop
            if(dna.strand(i).n_base - bgn_pos < para_max_cut_stap) then
                do k = 0, 11, 11
                    call space(k, 19)
                    write(k, "(a$ )"), "|-->last region, cutted length : "
                    write(k, "(i4$)"), dna.strand(i).n_base - bgn_pos
                    write(k, "(a$ )"), ", cutted pos :   ++"
                    write(k, "(a$ )"), ", remained length :   ++"
                    write(k, "(a  )"), " --> 14nt cutting"
                end do
                exit
            end if

            ! Check final staple and print information
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - bgn_pos

                ! If the last staple exceeds maximum length
                if(final_length > para_max_cut_stap) then

                    jj       = j
                    cng_para = para_gap_xover_nick

                    ! Find previous region
                    do
                        ! If the region length is longer than paramter value (default : 8)
                        if(region(jj).length >= cng_para * 2 + 2) then
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos

                            ! Check cutted and remained length
                            if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                                pre_pos - bgn_pos <= para_max_cut_stap .and. &
                                dna.strand(i).n_base - pre_pos >= para_min_cut_stap ) exit
                        end if

                        ! Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_reg .or. jj == 0) then

                            jj       = j
                            cng_para = cng_para - 1

                            if(cng_para == 0) then
                                deallocate(region)
                                return
                            end if
                        end if
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(0, "(a)"), " WARNING : This strand length exceeds para_max_cut_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_cut_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_cut_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_cut_stap) ) then
                        do k = 0, 11, 11
                            write(k, "(i20, a, i4)"), j, "-final region3, cutted length : ", &
                                dna.strand(i).n_base - bgn_pos
                        end do
                        cycle
                    end if

                    ! Add # of staple
                    up_base    = dna.top(pre_base).up
                    dna.n_stap = dna.n_stap + 1
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    do k = 0, 11, 11
                        write(k, "(i20, a, i4)"), jj, "-final region4, cutted length : ", &
                            final_length - (dna.strand(i).n_base - pre_pos)
                        write(k, "(i20, a, i4)"), jj, "-final region5, cutted length : ", &
                            dna.strand(i).n_base - pre_pos
                    end do
                else
                    do k = 0, 11, 11
                        write(k, "(i20, a, i4)"), j, "-final region6, cutted length : ", &
                            dna.strand(i).n_base - bgn_pos
                    end do
                end if
            end if
        end do
        write(0, "(a)"); write(11, "(a)")

        ! Deallocate memory
        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design_Mix

! ---------------------------------------------------------------------------------------

! Build sequence design with non-circular staple strands with 14nt seeds
! Last updated on Thursday 10 November 2016 by Hyungmin
subroutine SeqDesign_Build_Sequence_Design_Opt(prob, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, cen_base, pre_base, up_base, pre_reg
    integer :: n_region, cn_tn, length, final_length, cen_pos, pre_pos, bgn_pos
    integer :: cng_para, up, upup, xover, upxover
    logical :: b_ext, b_cut, b_14nt

    ! Loop to make short staples with 14bt seed as much as possible
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! For short staple, para_max_cut_stap <= 60
        if(dna.strand(i).n_base <= para_max_cut_stap) cycle

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! Skip for vertex slightly long staple < 70
        if(cn_tn > 0 .and. prob.sel_sec == 3 .and. dna.strand(i).n_base < 70) cycle

        ! For vertex in case of DX tile design
        !if(cn_tn > 0 .and. prob.sel_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! ==================================================
        !
        ! Build region of staple strands
        !
        ! ==================================================
        ! Tn loop |<----->| : 7 nt poly T loop
        !     ~~~~~=======*===========*=====*~~~~~~~
        !                 |<--------->|<--->|
        !           11 and 5 nt region not including crossover and nick
        !
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Print information on region of staple strands
        do k = 0, 11, 11
            write(k, "(i$    )"), i
            write(k, "(a, i5$)"), " - th strand, # of total bases : ", dna.strand(i).n_base
            write(k, "(a, i5 )"), ", # of bases in Tn : ", cn_tn

            do j = 1, n_region
                write(k, "(i20, a$)"), j, "  -  th region"
                write(k, "(a,  i3$)"), ", type : ",          region(j).types
                write(k, "(a,  i4$)"), ", region length : ", region(j).length
                write(k, "(a,  i4$)"), ", start pos : ",     region(j).sta_pos
                write(k, "(a,  i4$)"), ", center pos : ",    region(j).cen_pos
                write(k, "(a,  i4 )"), ", end pos : ",       region(j).end_pos
            end do
            write(k, "(a)"); call Space(k, 19)
            write(k, "(a)"), "----- Make nick position -----------------------------------------------------------------------------"
        end do

        ! ==================================================
        !
        ! Cut staples to make short multi staples with 14nt seeds
        !
        ! ==================================================
        bgn_pos = 0
        pre_reg = 0
        b_cut   = .false.
        b_14nt  = .false.
        j       = 0

        do
            ! Check loop
            j = j + 1
            if(j == n_region + 1) exit

            ! Skip when the region length is smaller than 5
            if(j /= n_region .and. region(j).length < 5) cycle

            ! Cutted length from center to beginning
            length = region(j).cen_pos - bgn_pos
            if(j == n_region) length = dna.strand(i).n_base - bgn_pos + 1

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_cut_stap) then
                b_cut = .true.
            end if

            ! Skip this current loop to prevent cutting this region
            if(j /= n_region .and. length <= para_max_cut_stap) then
                if( (b_14nt == .false. .and. region(j).length + 1 >= 14 .and. region(j).types == 1) .or. &
                    (b_14nt == .false. .and. region(j).length + 2 >= 14 .and. region(j).types == 2) ) then
                    b_14nt = .true.
                    cycle
                end if
            end if

            if(b_cut == .true. .and. b_14nt == .true. .and. length <= para_max_cut_stap) then
                !
                ! ==================================================
                !
                ! Optimal cutting to have 14nt seeds
                !
                ! ==================================================
                !
                ! Centered base ID and position
                cen_base = region(j).cen_base
                cen_pos  = region(j).cen_pos

                ! Option 1 - Secondary 14nt seeds will be introduced
                !if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                !    (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then
                !    cycle
                !end if

                ! Option 2 - Make strand long length
                do
                    exit
                    cen_base = region(j).cen_base
                    cen_pos  = region(j).cen_pos

                    ! If the cutting length is larger than maximum length
                    if(cen_pos - bgn_pos > para_max_cut_stap) then
                        j = j - 1
                        cen_base = region(j).cen_base
                        cen_pos  = region(j).cen_pos
                        exit
                    end if

                    if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                        (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then
                        if(j == n_region) exit
                        j = j + 1
                    else
                        j = j + 1
                        !exit
                    end if
                end do

                ! Option 3 - To avoid short region cutting
                !If(region(j).length < 10) cycle

                ! If the final length is smaller than minimum length
                if(dna.strand(i).n_base - cen_pos < para_min_cut_stap) cycle

                ! Print progress
                do k = 0, 11, 11
                    write(k, "(i20, a$)"), j, " -> cut region"
                    write(k, "(a,  i4$)"), ", cutted length : ",   cen_pos - bgn_pos
                    write(k, "(a,  i4$)"), ", cutted pos : ",      cen_pos
                    write(k, "(a,  i4$)"), ", remained length : ", dna.strand(i).n_base - cen_pos
                    write(k, "(a      )"), " --> 14nt cutting"
                end do

                ! Increase the number of staples
                dna.n_stap = dna.n_stap + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(cen_base).up
                dna.top(cen_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update beginning position and flags
                bgn_pos = cen_pos
                b_cut   = .false.
                b_14nt  = .false.
            else if(b_cut == .true. .and. length > para_max_cut_stap) then
                !
                ! ==================================================
                !
                ! If the length exceeds the maximum cutting length
                !
                ! ==================================================
                !
                ! To avoid small region in previous region
                jj       = j
                b_ext    = .false.
                cng_para = para_gap_xover_nick

                ! Find previous region
                do
                    ! If the region length is longer than paramter value (default : 8)
                    if(region(jj).length >= cng_para * 2 + 2) then
                        pre_base = region(jj).cen_base
                        pre_pos  = region(jj).cen_pos

                        ! Check cutted and remained length
                        if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                            pre_pos - bgn_pos <= para_max_cut_stap .and. &
                            dna.strand(i).n_base - pre_pos >= para_min_cut_stap ) exit
                    end if

                    ! Go back previous region
                    jj = jj - 1

                    ! Exception that the region cutted already
                    if(jj == pre_reg .or. jj == 0) then

                        ! Reset parameter and region index
                        jj       = j
                        cng_para = cng_para - 1

                        if(cng_para == 0) then
                            b_ext    = .true.
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos
                            exit
                        end if

                        if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                        if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                        if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                        do k = 0, 11, 11
                            call space(k, 16)
                            write(k, "(a)"), "** Adjusting parameter, para_gap_xover_nick : From "//&
                                trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                                trim(adjustl(Int2Str(cng_para)))//" - "//&
                                trim(adjustl(Int2Str(cng_para * 2 + 2)))
                        end do
                    end if
                end do

                ! No cutting due to exception, which may exceeds 60
                if(b_ext == .true.) then
                    do k = 0, 11, 11
                        call space(k, 19)
                        write(k, "(a$)"), "|-last region with exception, remaining length : "
                        write(k, "(i4)"), dna.strand(i).n_base - bgn_pos
                    end do
                    if(j == n_region) cycle
                end if

                ! If the cutted length is smaller than para_min_cut_stap
                if(pre_pos - bgn_pos < para_min_cut_stap) then
                    cycle
                end if

                ! Check 14nt seed cutting
                up      = dna.top(region(jj).end_base).up
                xover   = dna.top(up).xover
                upxover = dna.top(dna.top(up).up).xover

                if( region(jj).end_pos + 1 - bgn_pos >= para_min_cut_stap .and. &
                    region(jj).end_pos + 1 - bgn_pos <= para_max_cut_stap .and. &
                    dna.strand(i).n_base - region(jj).end_pos + 1 >= para_min_cut_stap .and. &
                    region(jj).length >= 12 .and. xover /= -1 .and. upxover /= -1 .and. para_set_stap_sxover == "on") then

                    ! Print progress
                    do k = 0, 11, 11
                        write(k, "(i20, a$)"), jj, " -> cut region"
                        write(k, "(a, i4$ )"), ", cutted length : ",   region(jj).end_pos + 1 - bgn_pos
                        write(k, "(a, i4$ )"), ", cutted pos : ",      region(jj).end_pos + 1
                        write(k, "(a, i4$ )"), ", remained length : ", dna.strand(i).n_base - region(jj).end_pos + 1
                        write(k, "(a      )"), " --> max cutting with single xover"
                    end do

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    if(dna.top(up).up == xover) dna.top(up).up = -1
                    if(dna.top(up).dn == xover) dna.top(up).dn = -1

                    if(dna.top(xover).up == up) dna.top(xover).up = -1
                    if(dna.top(xover).dn == up) dna.top(xover).dn = -1

                    dna.top(up).xover    = -1
                    dna.top(xover).xover = -1

                    dna.n_xover_stap  = dna.n_xover_stap  - 1
                    dna.n_sxover_stap = dna.n_sxover_stap + 1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = region(jj).end_pos + 1
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                else

                    ! Print progress
                    do k = 0, 11, 11
                        write(k, "(i20, a$)"), jj, " -> cut region"
                        write(k, "(a,  i4$)"), ", cutted length : ",   pre_pos - bgn_pos
                        write(k, "(a,  i4$)"), ", cutted pos : ",      pre_pos
                        write(k, "(a,  i4$)"), ", remained length : ", dna.strand(i).n_base - pre_pos
                        write(k, "(a      )"), " --> max cutting with nick"
                    end do

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = pre_pos
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                end if
            end if

            ! If remained length is smaller than para_max_cut_stap, exit this loop
            if(dna.strand(i).n_base - bgn_pos < para_max_cut_stap) then
                do k = 0, 11, 11
                    call space(k, 19)
                    write(k, "(a$ )"), "|-->last region, cutted length : "
                    write(k, "(i4$)"), dna.strand(i).n_base - bgn_pos
                    write(k, "(a$ )"), ", cutted pos :   ++"
                    write(k, "(a$ )"), ", remained length :   ++"
                    write(k, "(a  )"), " --> 14nt cutting"
                end do
                exit
            end if

            ! Check final staple and print information
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - bgn_pos

                ! If the last staple exceeds maximum length
                if(final_length > para_max_cut_stap) then

                    jj       = j
                    cng_para = para_gap_xover_nick

                    ! Find previous region
                    do
                        ! If the region length is longer than paramter value (default : 8)
                        if(region(jj).length >= cng_para * 2 + 2) then
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos

                            ! Check cutted and remained length
                            if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                                pre_pos - bgn_pos <= para_max_cut_stap .and. &
                                dna.strand(i).n_base - pre_pos >= para_min_cut_stap ) exit
                        end if

                        ! Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_reg .or. jj == 0) then

                            jj       = j
                            cng_para = cng_para - 1

                            if(cng_para == 0) then
                                deallocate(region)
                                return
                            end if
                        end if
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(0, "(a)"), " WARNING : This strand length exceeds para_max_cut_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_cut_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_cut_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_cut_stap) ) then
                        do k = 0, 11, 11
                            write(k, "(i20, a, i4)"), j, "-final region3, cutted length : ", &
                                dna.strand(i).n_base - bgn_pos
                        end do
                        cycle
                    end if

                    ! Add # of staple
                    up_base    = dna.top(pre_base).up
                    dna.n_stap = dna.n_stap + 1
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    do k = 0, 11, 11
                        write(k, "(i20, a, i4)"), jj, "-final region4, cutted length : ", &
                            final_length - (dna.strand(i).n_base - pre_pos)
                        write(k, "(i20, a, i4)"), jj, "-final region5, cutted length : ", &
                            dna.strand(i).n_base - pre_pos
                    end do
                else
                    do k = 0, 11, 11
                        write(k, "(i20, a, i4)"), j, "-final region6, cutted length : ", &
                            dna.strand(i).n_base - bgn_pos
                    end do
                end if
            end if
        end do
        write(0, "(a)"); write(11, "(a)")

        ! Deallocate memory
        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design_Opt

! ---------------------------------------------------------------------------------------

! Build sequence design with maximum cutting
! Last updated on Wednesday 2 November 2016 by Hyungmin
subroutine SeqDesign_Build_Sequence_Design_Max(prob, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, cen_base, pre_base, up_base, pre_region, cng_para
    integer :: n_region, cn_tn, length, final_length, cen_pos, pre_pos, bgn_pos, base1, base2
    logical :: b_cut, b_ext

    ! Loop for sequence design
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! For short staple, para_max_cut_stap <= 60
        if(dna.strand(i).n_base <= para_max_cut_stap) cycle

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! For vertex in case of DX tile design
        !if(cn_tn > 0 .and. prob.sel_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! Skip for vertex slightly long staple < 70
        if(cn_tn > 0 .and. prob.sel_sec == 3 .and. dna.strand(i).n_base < 70) cycle

        ! ==================================================
        !
        ! Build region of staple strands
        !
        ! ==================================================
        ! Tn loop |<----->| : 7 nt poly T loop
        !     ~~~~~=======*===========*=====*~~~~~~~
        !                 |<--------->|<--->|
        !           11 and 5 nt region not including crossover and nick
        !
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Print information on staple region
        do k = 0, 11, 11
            write(k, "(i$    )"), i
            write(k, "(a, i5$)"), " - th strand, # of total bases : ", dna.strand(i).n_base
            write(k, "(a, i5 )"), ", # of bases in Tn : ", cn_tn
            do j = 1, n_region
                write(k, "(i20, a$)"), j, "  -  th region"
                write(k, "(a,  i3$ )"), ", type : ",          region(j).types
                write(k, "(a,  i4$ )"), ", region length : ", region(j).length
                write(k, "(a,  i4$ )"), ", start pos : ",     region(j).sta_pos
                write(k, "(a,  i4$ )"), ", center pos : ",    region(j).cen_pos
                write(k, "(a,  i4  )"), ", end pos : ",       region(j).end_pos
            end do
            write(k, "(a)"); call Space(k, 19)
            write(k, "(a)"), "----- Make nick position -----------------------------------------------------------------------------"
        end do

        ! ==================================================
        !
        ! Cut staple strand to make short multi staples
        !
        ! ==================================================
        bgn_pos    = 0
        pre_region = 0
        b_cut      = .false.

        do j = 1, n_region

            ! Find centered base and position
            cen_base = region(j).cen_base
            cen_pos  = region(j).cen_pos
            length   = cen_pos - bgn_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_cut_stap) b_cut = .true.

            ! Staple cutting with maximum length
            if(b_cut == .true. .and. length >= para_max_cut_stap) then

                ! To avoid small region that will be 4nt region
                jj       = j - 1
                b_ext    = .false.
                cng_para = para_gap_xover_nick

                ! Find previous large enough region
                do
                    ! If the region length exceed, para_gap_xover_nick = 3 (default : 8)
                    if( (region(jj).types == 1 .and. region(jj).length >= cng_para * 2 + 2 + 1) .or. &
                        (region(jj).types == 2 .and. region(jj).length >= cng_para * 2 + 2) ) then
                        pre_base = region(jj).cen_base
                        pre_pos  = region(jj).cen_pos

                        ! Check remained staple length
                        if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                            pre_pos - bgn_pos <= para_max_cut_stap .and. &
                            dna.strand(i).n_base - pre_pos >= para_min_cut_stap) exit
                    end if

                    !Go back previous region
                    jj = jj - 1

                    ! Exception that the region cutted already
                    if(jj == pre_region .or. jj == 0) then

                        ! Reset parameter and region index
                        jj       = j - 1
                        cng_para = cng_para - 1
                        if(cng_para == 0) then
                            b_ext    = .true.
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos
                            exit
                        end if

                        ! Print information on changing parameter
                        if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                        if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                        if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                        do k = 0, 11, 11
                            call space(k, 16)
                            write(k, "(a)"), "** Adjusting parameter, para_gap_xover_nick : From "//&
                                trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                                trim(adjustl(Int2Str(cng_para)))//" - "//&
                                trim(adjustl(Int2Str(cng_para * 2 + 2)))
                        end do
                    end if
                end do

                ! No cutting due to exception, which might exceed 60 edge length
                if(b_ext == .true.) then
                    do k = 0, 11, 11
                        call space(k, 19)
                        write(k, "(a$)"), "|-last region with exception, remaining length : "
                        write(k, "(i4)"), dna.strand(i).n_base - bgn_pos
                    end do
                    if(j == n_region) cycle
                end if

                ! Print progress
                do k = 0, 11, 11
                    write(k, "(i20, a$)"), jj, " -> cut region"
                    write(k, "(a,  i4$)"), ", cutted length : ",   pre_pos - bgn_pos
                    write(k, "(a,  i4$)"), ", cutted pos : ",      pre_pos
                    write(k, "(a,  i4$)"), ", remained length : ", dna.strand(i).n_base - pre_pos
                    write(k, "(a      )"), " --> max cutting with nick"
                end do

                ! 1 : with single crossover
                if(1 .and. region(jj).types == 1 .and. region(jj).length <= 9) then

                    ! Make single crossover
                    if(dna.top(dna.top(region(jj).end_base).up).across == -1) then
                        ! Back crossover
                        base1      = dna.top(region(jj).sta_base).dn
                        base2      = dna.top(base1).dn
                        pre_pos    = region(jj-1).end_pos + 1
                        pre_region = jj - 1
                    else
                        ! Front crossover
                        base1      = dna.top(region(jj).end_base).up
                        base2      = dna.top(base1).up
                        pre_pos    = region(jj).end_pos + 1
                        pre_region = jj
                    end if

                    ! Make single crossover
                    do k = 0, 11, 11
                        call space(k, 16)
                        write(k, "(a$)"), "***** Single Xover"
                        write(k, "(a,  i4$)"), ", cutted length : ",   pre_pos - bgn_pos
                        write(k, "(a,  i4$)"), ", cutted pos : ",      pre_pos
                        write(k, "(a,  i4 )"), ", remained length : ", dna.strand(i).n_base - pre_pos
                    end do

                    ! Cut crossover and make new connectivity
                    dna.top(base1).xover = -1
                    dna.top(base2).xover = -1
                    dna.top(base1).dn    = -1
                    dna.top(base2).up    = -1
                    dna.n_stap           = dna.n_stap + 1
                    dna.n_xover_stap     = dna.n_xover_stap - 1
                    dna.n_sxover_stap    = dna.n_sxover_stap + 1

                    ! Update flag
                    bgn_pos = pre_pos
                    b_cut   = .false.
                else

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    ! Update starting position and flag
                    pre_region = jj
                    bgn_pos    = pre_pos
                    b_cut      = .false.
                end if
            end if

            ! Check remained length whether it is smaller than para_max_cut_stap
            if(dna.strand(i).n_base - bgn_pos <= para_max_cut_stap) then
                do k = 0, 11, 11
                    write(k, "(i20, a, i4)"), j, "-final region2, cutted length : ", &
                        dna.strand(i).n_base - bgn_pos
                end do
                exit
            end if

            ! Check if it is the final region
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - bgn_pos

                ! If the last staple exceeds maximum length
                if(final_length >= para_max_cut_stap) then

                    jj       = j
                    cng_para = para_gap_xover_nick
                    do
                        ! Exception
                        if(jj == 0) then
                            deallocate(region)
                            return
                        end if

                        if(region(jj).length >= para_gap_xover_nick * 2 + 2) then
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos

                            ! Check remained staple length
                            if( pre_pos - bgn_pos >= para_min_cut_stap .and. &
                                pre_pos - bgn_pos <= para_max_cut_stap .and. &
                                dna.strand(i).n_base - pre_pos >= para_min_cut_stap) exit
                        end if

                        !Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_region .or. jj == 0) then

                            ! Reset parameter and region index
                            jj       = j - 1
                            cng_para = cng_para - 1
                            if(cng_para == 0) then
                                b_ext    = .true.
                                pre_base = region(jj).cen_base
                                pre_pos  = region(jj).cen_pos
                                exit
                            end if

                            ! Print information on changing parameter
                            if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                            if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                            if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                            do k = 0, 11, 11
                                call space(k, 16)
                                write(k, "(a)"), "** Adjusting parameter, para_gap_xover_nick : From "//&
                                    trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                                    trim(adjustl(Int2Str(cng_para)))//" - "//&
                                    trim(adjustl(Int2Str(cng_para * 2 + 2)))
                            end do
                        end if
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(0, "(a)"), " WARNING : This strand length exceeds para_max_cut_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_cut_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_cut_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_cut_stap) ) then
                        do k = 0, 11, 11
                            write(k, "(i20, a, i5)"), j, "-final region3, cutted length : ", &
                                dna.strand(i).n_base - bgn_pos
                        end do
                        cycle
                    end if

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    do k = 0, 11, 11
                        write(k, "(i20, a, i5)"), jj, "-final region4, cutted length : ", &
                            final_length - (dna.strand(i).n_base - pre_pos)
                        write(k, "(i20, a, i5)"), jj, "-final region5, cutted length : ", &
                            dna.strand(i).n_base - pre_pos
                    end do
                else
                    do k = 0, 11, 11
                        write(k, "(i20, a, i5)"), j, "-final region6, cutted length : ", &
                            dna.strand(i).n_base - bgn_pos
                    end do
                end if
            end if

        end do
        write(0, "(a)"); write(11, "(a)")

        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design_Max

! ---------------------------------------------------------------------------------------

! Build staple region without crossovers and unpaired nucleotides
! Last updated on Wednesday 2 November 2016 by Hyungmin
subroutine SeqDesign_Build_Region_Staple(dna, i, region, n_region)
    type(DNAType),    intent(inout) :: dna
    type(RegionType), intent(inout) :: region(:)
    integer,          intent(in)    :: i
    integer,          intent(inout) :: n_region

    integer :: j, k, base, across, sta_base, cen_base, end_base, cen_pos, len_cen
    logical :: b_region, b_vertex

    ! Find the starting point (down == -1)
    base = Mani_Go_Start_Base(dna, i)

    n_region = 0
    b_region = .false.
    b_vertex = .false.

    do j = 1, dna.strand(i).n_base

        ! Find across ID
        across = dna.top(base).across

        if(across == -1) then

            ! For Tn loop
            b_region = .true.
            b_vertex = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
        else
            ! Check dn, up, xover and across's over
            ! Starting base always has -1 (down = -1)
            if( dna.top(base).dn      == -1 .or. &
                dna.top(base).up      == -1 .or. &
                dna.top(base).xover   /= -1 .or. &
                dna.top(across).xover /= -1 ) then

            b_region = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
            end if
        end if

        if(b_region == .true.) then

            ! Make new region
            b_region = .false.
            n_region = n_region + 1

            ! Set region data
            region(n_region).sta_base = base
            region(n_region).length   = 1
            region(n_region).sta_pos  = j
            region(n_region).end_pos  = j

            ! Set vertex region
            if(b_vertex == .true.) then
                region(n_region).types   = 1
                region(n_region-1).types = 1
                b_vertex = .false.
            else
                region(n_region).types = 2
            end if
        else

            ! Increase length
            region(n_region).length  = region(n_region).length + 1
            region(n_region).end_pos = j
        end if

        ! Update base to go up
        base = dna.top(base).up
    end do

    ! Set center/end base and position
    do j = 1, n_region

        len_cen = (region(j).length+1)/2 - 1

        ! If the vetrex region
        if(region(j).types == 1) then
            if(dna.top(dna.top(region(j).sta_base).dn).across /= -1) then
                len_cen = len_cen - 1
            else
                len_cen = len_cen + 1
            end if
        end if

        ! Find center base
        cen_base = region(j).sta_base
        do k = 1, len_cen
            cen_base = dna.top(cen_base).up
        end do
        region(j).cen_base = cen_base

        ! Set center position
        region(j).cen_pos = region(j).sta_pos + len_cen

        ! Find end base
        end_base = region(j).sta_base
        do k = 1, region(j).length - 1
            end_base = dna.top(end_base).up
        end do
        region(j).end_base = end_base

        ! Set end position
        region(j).end_pos = region(j).sta_pos + (region(j).length - 1)

        ! Set nucleotide of 14nt seeds
        if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. & 
            (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then

            sta_base = region(j).sta_base
            dna.top(dna.top(sta_base).dn).b_14nt = .true.
            dna.top(sta_base).b_14nt             = .true.

            do k = 1, region(j).length - 1
                sta_base = dna.top(sta_base).up
                dna.top(sta_base).b_14nt = .true.
            end do

            if(dna.top(dna.top(sta_base).up).across /= -1) then
                dna.top(dna.top(sta_base).up).b_14nt = .true.
            end if
        end if
    end do
end subroutine SeqDesign_Build_Region_Staple

! ---------------------------------------------------------------------------------------

! Build sequence design with non-circular staple strands
! Last updated on Thursday 7 June 2016 by Hyungmin
subroutine SeqDesign_Build_Sequence_Design(prob, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type :: RegionType
        integer :: types        ! 1 - vertex, 2 - edge
        integer :: length       ! total length at the region not including boundary

        integer :: strt_pos     ! start position
        integer :: end_pos      ! end position
        integer :: cntr_pos     ! center position

        integer :: strt_base    ! start base
        integer :: cntr_base    ! center base
    end type RegionType

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, across, cntr_base, pre_base, up_base, pre_region
    integer :: n_region, cn_tn, length, final_length, cntr_pos, pre_pos, begin_pos
    logical :: b_region, b_cut, b_ext, b_vertex

    ! Loop to make short staple strand
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! For short staple (less than para_max_cut_stap, 60)
        if(dna.strand(i).n_base < para_max_cut_stap) cycle

        ! Count bases if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! ==================================================
        !
        ! For vertex for DX tile
        !
        ! ==================================================
        !if(cn_tn > 0 .and. prob.sel_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)

        ! --------------------------------------------------
        !
        ! Build region array
        !
        ! --------------------------------------------------
        ! Tn loop |<---------->| 12 bp
        !     ~~~~~============*===========*=====*  : base pair
        !                      |<--------->|
        !               region 11 not including crossover and nick region
        n_region = 0
        allocate(region(dna.strand(i).n_base))

        do j = 1, dna.strand(i).n_base

            ! Find across ID
            across = dna.top(base).across

            if(across == -1) then

                ! For Tn loop
                b_region = .true.
                b_vertex = .true.

                ! Update base to go up
                base = dna.top(base).up
                cycle
            else
                ! Check down, up, xover and across
                ! Starting base always has -1 (down = -1)
                if( dna.top(base).dn      == -1 .or. &
                    dna.top(base).up      == -1 .or. &
                    dna.top(base).xover   /= -1 .or. &
                    dna.top(across).xover /= -1 ) then

                    b_region = .true.

                    ! Update base to go up
                    base = dna.top(base).up
                    cycle
                end if
            end if

            if(b_region == .true.) then

                ! Make new region
                b_region = .false.
                n_region = n_region + 1

                ! Set region data
                region(n_region).strt_base = base
                region(n_region).length    = 1
                region(n_region).strt_pos  = j
                region(n_region).end_pos   = j
                region(n_region).types     = 2

                ! Set vertex region
                if(b_vertex == .true.) then
                    region(n_region).types   = 1
                    region(n_region-1).types = 1
                    b_vertex = .false.
                end if
            else

                ! Increase length
                region(n_region).length  = region(n_region).length + 1
                region(n_region).end_pos = j

                ! Set vertex region
                if(b_vertex == .true.) then
                    region(n_region).types   = 1
                    region(n_region-1).types = 1
                    b_vertex = .false.
                end if
            end if

            ! Update base to go up
            base = dna.top(base).up
        end do

        ! Set centered base and position
        do j = 1, n_region

            ! Find centered base
            cntr_base = region(j).strt_base
            do k = 1, (region(j).length + 1) / 2 - 1
                cntr_base = dna.top(cntr_base).up
            end do
            region(j).cntr_base = cntr_base

            ! Set centered position
            region(j).cntr_pos = region(j).strt_pos + ((region(j).length+1)/2-1)
        end do

        ! Print information on strand and region
        do k = 0, 11, 11
            write(k, "(i, 2(a, i5))"), i, &
                " - th strand, # of total bases : ", dna.strand(i).n_base, &
                ", # of bases in Tn : ", cn_tn
            do j = 1, n_region
                write(k, "(i20, a$)"), j, "  -  th region"
                write(k, "(a,  i3$ )"), ", type : ",          region(j).types
                write(k, "(a,  i4$ )"), ", region length : ", region(j).length
                write(k, "(a,  i4$ )"), ", start pos : ",     region(j).strt_pos
                write(k, "(a,  i4$ )"), ", center pos : ",    region(j).cntr_pos
                write(k, "(a,  i4  )"), ", end pos : ",       region(j).end_pos
            end do
            write(k, "(a)")
            call Space(k, 19)
            write(k, "(a)"), "----- Make nick position -----------------------------------------------------------------------------"
        end do

        ! ==================================================
        !
        ! Cut staple strand to make it short
        !
        ! ==================================================
        begin_pos  = 0
        pre_region = 0
        b_cut      = .false.
        do j = 1, n_region

            if( para_cut_stap_method == "min" .and. &
                region(j).length < para_gap_xover_nick*2+2 ) cycle

            ! Find centered base and pos
            cntr_base = region(j).cntr_base
            cntr_pos  = region(j).cntr_pos
            length    = cntr_pos - begin_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_cut_stap) then
                b_cut = .true.
            end if

            ! Cutting method # 1 (minimum length staples)
            ! If minimum cutting length exceed, cut strand at the centered region
            if(para_cut_stap_method == "min" .and. b_cut == .true.) then

                pre_base = region(j).cntr_base
                pre_pos  = region(j).cntr_pos

                ! Skip the end iteration
                if(dna.strand(i).n_base - pre_pos < para_min_cut_stap) then
                    do k = 0, 11, 11
                        write(k, "(i20, a, i5)"), j, " -final region, cutted length : ", &
                            dna.strand(i).n_base - begin_pos
                    end do
                    exit
                end if

                do k = 0, 11, 11
                    write(k, "(i20, a, 2(a, i5))"), j, " -> cut region", &
                        ", cutted length : ",   pre_pos-begin_pos,       &
                        ", remained length : ", dna.strand(i).n_base - pre_pos
                end do

                ! Add # of staple
                dna.n_stap = dna.n_stap + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(pre_base).up
                dna.top(pre_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update new count and flag
                begin_pos = pre_pos
                b_cut     = .false.
            end if

            ! Cutting method # 2 and 3 (maximum or optimal length staples)
            if((para_cut_stap_method == "max" .or. para_cut_stap_method == "mid") .and. b_cut == .true.) then

                if( (para_cut_stap_method == "max" .and. length >= para_max_cut_stap) .or. &
                    (para_cut_stap_method == "mid" .and. length >= para_mid_cut_stap) ) then

                    ! To avoid small region in previous region
                    jj    = j - 1
                    b_ext = .false.

                    ! Find previous region
                    do
                        ! If the region length is longer than paramter value (default : 8)
                        if(region(jj).length >= para_gap_xover_nick*2+2) then
                            pre_base = region(jj).cntr_base
                            pre_pos  = region(jj).cntr_pos

                            ! Check remained length for the end of the strand
                            if(dna.strand(i).n_base - pre_pos >= para_min_cut_stap) exit
                        end if

                        !Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_region) then
                            b_ext = .true.
                            exit
                        end if
                    end do

                    ! Make staple strand that exceeds para_max_cut_stap/para_mid_cut_stap
                    if(b_ext == .true.) then
                        if(j == n_region) then
                            do k = 0, 11, 11
                                write(k, "(i20, a, i5)"), j, "-final region1, cutted length : ", &
                                    dna.strand(i).n_base - begin_pos
                            end do
                        end if
                        cycle
                    end if

                    ! If the cutted length is smaller than para_min_cut_stap
                    if(pre_pos-begin_pos < para_min_cut_stap) cycle

                    do k = 0, 11, 11
                        write(k, "(i20, a, 3(a, i5))"), jj, " -> cut region", &
                            ", cutted length : ",   pre_pos-begin_pos,        &
                            ", cutted pos : ",      pre_pos,                  &
                            ", remained length : ", dna.strand(i).n_base - pre_pos
                    end do

                    ! Add # of staple
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1
                    pre_region           = jj

                    ! Update starting position and flag
                    begin_pos = pre_pos
                    b_cut     = .false.
                end if
            end if

            ! If remained length is smaller than para_max_cut_stap, exit this loop
            if( (para_cut_stap_method == "max" .and. dna.strand(i).n_base-begin_pos < para_max_cut_stap) .or. &
                (para_cut_stap_method == "mid" .and. dna.strand(i).n_base-begin_pos < para_mid_cut_stap) ) then
                do k = 0, 11, 11
                    write(k, "(i20, a, i5)"), j, "-final region2, cutted length : ", &
                        dna.strand(i).n_base - begin_pos
                end do
                exit
            end if

            ! Check final staple and print information (for only mid/max)
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - begin_pos

                ! If the last staple exceeds maximum length
                if( (para_cut_stap_method == "max" .and. final_length >= para_max_cut_stap) .or. &
                    (para_cut_stap_method == "mid" .and. final_length >= para_mid_cut_stap) ) then

                    jj = j
                    do
                        ! Exception
                        if(jj == 0) then
                            deallocate(region)
                            return
                        end if

                        if(region(jj).length >= para_gap_xover_nick*2+2) then
                            pre_base = region(jj).cntr_base
                            pre_pos  = region(jj).cntr_pos
                            if(dna.strand(i).n_base - pre_pos >= para_min_cut_stap) exit
                        end if
                        jj = jj - 1
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(0, "(a)"), " WARNING : This strand length exceeds para_max_cut_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_cut_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_cut_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_cut_stap) ) then
                        do k = 0, 11, 11
                            write(k, "(i20, a, i5)"), j, "-final region3, cutted length : ", &
                                dna.strand(i).n_base - begin_pos
                        end do
                        cycle
                    end if

                    ! Add # of staple
                    up_base = dna.top(pre_base).up
                    dna.n_stap = dna.n_stap + 1
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    do k = 0, 11, 11
                        write(k, "(i20, a, i5)"), jj, "-final region4, cutted length : ", &
                            final_length - (dna.strand(i).n_base - pre_pos)
                        write(k, "(i20, a, i5)"), jj, "-final region5, cutted length : ", &
                            dna.strand(i).n_base - pre_pos
                    end do
                else
                    do k = 0, 11, 11
                        write(k, "(i20, a, i5)"), j, "-final region6, cutted length : ", &
                            dna.strand(i).n_base - begin_pos
                    end do
                end if
            end if

        end do
        write(0, "(a)"); write(11, "(a)")

        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design

! ---------------------------------------------------------------------------------------

! Make nick in scaffold strand
! Last updated on Thuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Make_Nick_Scaf(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Count # of crossovers based on cross-sectional edges
    type :: CroLType
        integer :: n_xover
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, node, sec, across, edge, base, dn_base, strt_base
    integer :: max_strt_base, max_count, count, max_edge

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    do i = 1, geom.n_croL
        croL(i).n_xover = 0
    end do

    ! Find the number of crossovers and bps in terms of croL
    do i = 1, dna.strand(1).n_base
        base = dna.strand(1).base(i)
        node = dna.top(base).node

        ! For unpaired nucleotide
        if(node == -1) cycle

        edge = mesh.node(node).croL

        ! If the base is crossover
        if(dna.top(base).xover /= -1) then
            croL(edge).n_xover = croL(edge).n_xover + 1
        end if
    end do

    !do i = 1, geom.n_croL
    !    print *, i, "-th : ", croL(i).n_xover
    !end do

    ! Loop to make nick position in scaffold strand
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).types == "stap") cycle

        ! Find first base(going down) that has single crossover
        base = dna.strand(i).base(1)
        do j = 1, dna.strand(i).n_base
            node = dna.top(base).node
            if(mesh.node(node).dn == -1) exit
            base = dna.top(base).dn
        end do

        ! Find maximum region in only non-crossover edges
        max_strt_base = base
        max_count     = 0
        count         = 0
        do j = 1, dna.strand(i).n_base

            ! To avoid node ID is negative
            do
                if(dna.top(base).node /= -1) then

                    node = dna.top(base).node
                    sec  = mesh.node(node).sec

                    ! -------------------------------------
                    ! To replece nick at the bottom section
                    ! -------------------------------------
                    if(geom.sec.posR(sec+1) == 1) exit
                end if
                base = dna.top(base).up
            end do

            node = dna.top(base).node
            edge = mesh.node(node).croL

            ! Find max region only in non-crossover edges
            if(croL(edge).n_xover == 0) then
                node   = dna.top(base).node
                across = dna.top(base).across

                if( dna.top(base).xover   /= -1 .or. &  ! If there is crossover
                    dna.top(across).xover /= -1 .or. &
                    dna.top(across).up    == -1 .or. &  ! If there is staple nick
                    dna.top(across).dn    == -1 .or. &
                    mesh.node(node).up    == -1 .or. &  ! If there is single crossover
                    mesh.node(node).dn    == -1 ) then

                    ! If the region exceeds max_count
                    if(max_count < count) then
                        max_count     = count
                        max_edge      = edge
                        max_strt_base = strt_base
                    end if

                    ! Set new starting base
                    count     = 0
                    strt_base = base
                else
                    count = count + 1
                end if
            end if

            ! Go to the upper base
            base = dna.top(base).up
        end do

        !max_count = max_count + 1

        ! Set base in the middle at the maximum region
        base = dna.top(max_strt_base).id
        do j = 1, max_count / 2 + 1
            base = dna.top(base).up
        end do
        dn_base = dna.top(base).dn

        ! Disconnect between current and downward bases
        dna.top(base).dn    = -1
        dna.top(dn_base).up = -1
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.4. Set nick position in scaffold strand"
        call Space(i, 11)
        write(i, "(a)"), "* The edge number at nick position         : "//trim(adjustl(Int2Str(max_edge)))
        write(i, "(a)")
    end do

    ! Deallocate memory
    deallocate(croL)
end subroutine SeqDesign_Make_Nick_Scaf

! ---------------------------------------------------------------------------------------

! Make short scaffold strand
! Last updated on Thuesday 28 June 2016 by Hyungmin
subroutine SeqDesign_Make_Short_Scaf(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type :: RegionType
        integer :: length       ! total length at the region not including boundary

        integer :: strt_pos     ! start position
        integer :: end_pos      ! end position
        integer :: cntr_pos     ! center position

        integer :: strt_base    ! start base
        integer :: cntr_base    ! center base
    end type RegionType

    type(RegionType), allocatable :: region(:)

    integer :: i, j, k, base, across, cntr_base, pre_base, up_base
    integer :: n_region, length, cntr_pos, pre_pos, begin_pos
    logical :: b_region, b_cut

    ! Exception current subroutine
    if(para_max_cut_scaf == 0) then
        return
    else if(para_max_cut_scaf == -1) then
        if(dna.strand(1).n_base > 7249) then
            para_max_cut_scaf = 3000
        else
            return
        end if
    end if

    ! Loop to make short scaffold strand
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).types == "stap") cycle

        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)

        ! --------------------------------------------------
        !
        ! Build region array
        !
        ! --------------------------------------------------
        ! Tn loop |<---------->| 12 bp
        !     ~~~~~============*===========*=====*  : base pair
        !                      |<--------->|
        !               region 11 not including crossover and nick region
        n_region = 0
        allocate(region(dna.strand(i).n_base))

        do j = 1, dna.strand(i).n_base

            ! Find across ID
            across = dna.top(base).across

            if(across == -1) then
                ! For Tn loop
                b_region = .true.

                ! Update base to go up
                base = dna.top(base).up
                cycle
            else
                ! Check down, up, xover and across
                ! Starting base always has -1 (down = -1)
                if( dna.top(base).dn      == -1 .or. &
                    dna.top(base).up      == -1 .or. &
                    dna.top(base).xover   /= -1 .or. &
                    dna.top(across).xover /= -1 ) then

                    b_region = .true.

                    ! Update base to go up
                    base = dna.top(base).up
                    cycle
                end if
            end if

            if(b_region == .true.) then
                ! Make new region
                b_region = .false.
                n_region = n_region + 1

                ! Set region data
                region(n_region).strt_base = base
                region(n_region).length    = 1
                region(n_region).strt_pos  = j
                region(n_region).end_pos   = j
            else
                ! Increase length
                region(n_region).length  = region(n_region).length + 1
                region(n_region).end_pos = j
            end if

            ! Update base to go up
            base = dna.top(base).up
        end do

        ! Set centered base and position
        do j = 1, n_region

            ! Find centered base
            cntr_base = region(j).strt_base
            do k = 1, (region(j).length+1)/2-1
                cntr_base = dna.top(cntr_base).up
            end do
            region(j).cntr_base = cntr_base

            ! Set centered position
            region(j).cntr_pos = region(j).strt_pos + ((region(j).length+1)/2-1)
        end do

        ! Print information on strand and region
        !write(0, "(i, a,i5)"), i, &
        !    " - th strand, # of total bases : ", dna.strand(i).n_base
        !write(0, "(a)")

        !do j = 1, n_region
        !    write(0, "(i20, a, 4(a, i5))"), j, " - th region", &
        !        ", region length : ", region(j).length,   &
        !        ", start pos : ",     region(j).strt_pos, &
        !        ", centerd pos : ",   region(j).cntr_pos, &
        !        ", end pos : ",       region(j).end_pos
        !end do
        !write(0, "(a)")
        !call Space(0, 19)
        !write(0, "(a)"), "--------------------------------------------------------------------------------"
        !write(0, "(a)")

        ! --------------------------------------------------
        !
        ! Cut strand to make it short
        !
        ! --------------------------------------------------
        begin_pos  = 0
        b_cut      = .false.
        do j = 1, n_region

            ! Find centered base and pos
            cntr_base = region(j).cntr_base
            cntr_pos  = region(j).cntr_pos
            length    = cntr_pos - begin_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length > para_max_cut_scaf) then
                b_cut = .true.
            end if

            ! Cutting method # 1 (minimum length staples)
            ! If minimum cutting length exceed, cut strand at the centered region
            if(b_cut == .true.) then

                pre_base = region(j).cntr_base
                pre_pos  = region(j).cntr_pos

                ! Skip the end iteration
                if(j == n_region) cycle

                write(0, "(i20, a, 2(a, i5))"), j, " -> cut region", &
                    ", cutted length : ",   pre_pos-begin_pos,       &
                    ", remained length : ", dna.strand(i).n_base - pre_pos

                ! Add # of staple
                dna.n_scaf = dna.n_scaf + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(pre_base).up
                dna.top(pre_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update new count and flag
                begin_pos = pre_pos
                b_cut     = .false.
            end if

            ! Print final strand
            if(j == n_region) then
                write(0, "(i20, a, i4)"), j, " -1 cut region, cutted length : ", &
                    dna.strand(i).n_base - pre_pos
            end if
        end do
        write(0, "(a)")

        deallocate(region)
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.5. Make short scaffold strand"
        call Space(i, 11)
        write(i, "(a)"), "* The maximum number of bases per scaffold : "//trim(adjustl(Int2Str(para_max_cut_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified scaffold strands  : "//trim(adjustl(Int2Str(dna.n_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified scaffold strands  : "//trim(adjustl(Int2Str(dna.n_scaf)))
        write(i, "(a)")
    end do
end subroutine SeqDesign_Make_Short_Scaf

! ---------------------------------------------------------------------------------------

! Move cur_base to avoid node without ID and crossovers
! Last updated on Tuesday 24 May 2016 by Hyungmin
function SeqDesign_Avoid_Barrier(mesh, dna, base, gap) result(cur_base)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: base
    integer,        intent(in) :: gap

    integer :: i, cur_base, chk_base, node, xover
    logical :: nofront, noback

    cur_base = base

    ! --------------------------------------------------
    ! Move cur_base to avoid node without ID
    ! --------------------------------------------------
    do
        ! Check front region if there is node without ID
        nofront  = .true.
        chk_base = cur_base

        do i = 1, gap
            node = dna.top(chk_base).node

            ! For staple strand
            if(node == -1) then
                noback = .false.
                exit
            end if

            ! For scaffold strand
            if(mesh.node(node).up == -1 .or. mesh.node(node).dn == -1) then
                nofront = .false.
                exit
            end if
            chk_base = dna.top(chk_base).up
        end do

        ! Check backward region if there is node without ID
        noback   = .true.
        chk_base = cur_base

        do i = 1, gap
            node = dna.top(chk_base).node

            ! For staple strand
            if(node == -1) then
                noback = .false.
                exit
            end if

            ! For scaffold strand
            if(mesh.node(node).up == -1 .or. mesh.node(node).dn == -1) then
                noback = .false.
                exit
            end if
            chk_base = dna.top(chk_base).dn
        end do

        ! If two critera was satified
        if(nofront == .true. .and. noback == .true.) exit

        ! Update the cur_base
        cur_base = dna.top(cur_base).dn
    end do

    ! --------------------------------------------------
    ! Move cur_base to crossover if there is crossover nearby
    ! --------------------------------------------------
    do
        ! Check front region whether there is crossovers
        nofront  = .true.
        chk_base = cur_base

        do i = 1, gap
            xover = dna.top(chk_base).xover
            if(xover /= -1) then
                exit
            end if
            chk_base = dna.top(chk_base).up
        end do

        if(xover /= -1) then
            cur_base = dna.top(chk_base).up
            if(dna.top(cur_base).xover == -1) then
                cur_base = dna.top(cur_base).dn
            end if
            exit
        end if

        ! Check backward region whether there is crossovers
        noback   = .true.
        chk_base = cur_base

        do i = 1, gap
            xover = dna.top(chk_base).xover
            if(xover /= -1) then
                exit
            end if
            chk_base = dna.top(chk_base).dn
        end do

        if(xover /= -1) then
            cur_base = chk_base
            exit
        end if

        ! If two critera was satified
        if(nofront == .true. .and. noback == .true.) exit

        ! Update the cur_base
        cur_base = dna.top(cur_base).dn
    end do
end function SeqDesign_Avoid_Barrier

! ---------------------------------------------------------------------------------------

! Make short strand
! Last updated on Monday 4 Apr 2016 by Hyungmin
subroutine SeqDesign_Make_Short_Strand(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, max_strand, rem_count, cur_base, down_base
    logical :: b_last

    ! Loop for all strand
    do i = 1, dna.n_strand

        ! Set the maximum number of bases per strand
        if(dna.strand(i).types == "scaf") then
            max_strand = para_max_cut_scaf
        else if(dna.strand(i).types == "stap") then
            max_strand = para_max_cut_stap
        end if

        if(max_strand == 0) cycle

        ! If it exceeds certain length, cut strand
        if(dna.strand(i).n_base > max_strand) then

            ! Find starting point
            cur_base = dna.strand(i).base(1)
            do
                if(dna.top(cur_base).dn == -1) exit
                cur_base = dna.top(cur_base).dn
            end do

            ! Cut strand
            do
                ! Count the remaining number of bases
                rem_count = SeqDesign_Count_Remainder(dna, cur_base)

                ! Go cur_base up
                if(rem_count > 2 * max_strand) then
                    b_last = .false.
                    do j = 1, max_strand
                        cur_base = dna.top(cur_base).up
                    end do
                else
                    b_last = .true.
                    do j = 1, rem_count / 2
                        cur_base = dna.top(cur_base).up
                    end do
                end if

                ! Move cur_base to avoid node without ID and crossovers
                cur_base  = SeqDesign_Avoid_Barrier(mesh, dna, cur_base, 5)
                down_base = dna.top(cur_base).dn

                ! Disconnect between current and downward bases
                dna.top(cur_base).dn  = -1
                dna.top(down_base).up = -1

                ! If the cur_base is crossover, disconnect the crossover
                if(dna.top(cur_base).xover == dna.top(down_base).id) then
                    dna.top(cur_base).xover  = -1
                    dna.top(down_base).xover = -1
                end if

                if(dna.strand(i).types == "scaf") then
                    dna.n_scaf = dna.n_scaf + 1
                else if(dna.strand(i).types == "stap") then
                    dna.n_stap = dna.n_stap + 1
                end if

                if(b_last == .true.) exit
            end do
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "6.4. Make short staple strand"
        call Space(i, 11)
        write(i, "(a$)"), "* The maximum number of bases per scaffold : "
        write(i, "(i7)"), para_max_cut_scaf
        call Space(i, 11)
        write(i, "(a$)"), "* The maximum number of bases per staple   : "
        write(i, "(i7)"), para_max_cut_stap
        call Space(i, 11)
        write(i, "(a$)"), "* The number of modified scaffold strands  : "
        write(i, "(i7)"), dna.n_scaf
        call Space(i, 11)
        write(i, "(a$)"), "* The number of modified staple strands    : "
        write(i, "(i7)"), dna.n_stap
        write(i, "(a )")
    end do
end subroutine SeqDesign_Make_Short_Strand

! ---------------------------------------------------------------------------------------

! Count the number of remainder bases, strand should be non-circular
! Last updated on Monday 4 Apr 2016 by Hyungmin
function SeqDesign_Count_Remainder(dna, cur_base) result(count)
    type(DNAType), intent(in) :: dna
    integer, intent(in) :: cur_base

    integer :: i, count, base

    count = 1
    base  = cur_base

    ! Find the number of remaining bases
    do
        base = dna.top(base).up
        if(base == -1) then
            exit
        else
            count = count + 1
        end if
    end do
end function SeqDesign_Count_Remainder

! ---------------------------------------------------------------------------------------

! Rebuild strand data from dnaTop data
! Last updated on Wednesday 30 Mar 2016 by Hyungmin
subroutine SeqDesign_Rebuild_Strand(dna)
    type(DNAType), intent(inout) :: dna

    logical, allocatable, dimension(:) :: b_visit
    type(ListBase), pointer :: list_base
    type(ListBase), pointer :: ptr_base
    type(TopType) :: cur_base

    integer :: i, j, n_strand, n_base, int_base
    logical :: b_end

    ! Deallocate previous strand data
    deallocate(dna.strand)

    ! Allocate and initialize strand data
    dna.n_strand = dna.n_scaf + dna.n_stap
    allocate(dna.strand(dna.n_strand))
    call Mani_Init_StrandType(dna.strand, dna.n_strand)

    ! Set strand type
    do i = 1, dna.n_strand
        if(i <= dna.n_scaf) then
            dna.strand(i).types = "scaf"
        else
            dna.strand(i).types = "stap"
        end if
    end do

    ! Nullify the linked lists
    nullify(list_base, ptr_base)

    ! Allocate and initilize the b_visit data
    allocate(b_visit(dna.n_top))
    do i = 1, dna.n_top
        b_visit(i) = .false.
    end do
    
    n_strand = 0
    do
        ! Check b_visit if there is 0 (0 means not visiting base)
        b_end = .true.
        do i = 1, dna.n_top
            if(b_visit(dna.top(i).id) == .false.) then
                ! All visit of all bases
                b_end = .false.
                exit
            end if
        end do
        if(b_end == .true.) exit

        ! Increase the number of strands
        n_strand = n_strand + 1
        cur_base = dna.top(i)
        int_base = cur_base.id

        ! Find the first base in the current strand
        do
            if(cur_base.dn == -1) then
                ! cur_base is at the 3'-end of the strand
                dna.strand(n_strand).b_circular = .false.
                exit
            else if(cur_base.dn == int_base) then
                ! cur_base goes back to the starting point
                dna.strand(n_strand).b_circular = .true.
                cur_base = dna.top(int_base)
                exit
            end if

            cur_base = dna.top(cur_base.dn)

            if(b_visit(cur_base.id)) then
                write(0, "(a$)"), "Error - Reached a visited base : "
                write(0, "(a )"), "SeqDesign_Rebuild_Strand"
                stop
            end if
        end do

        ! Walk through the current strand
        n_base = 1

        ! Insert data to linked list
        allocate(ptr_base)
        ptr_base%id = cur_base.id
        list_base => List_Insert_Base(list_base, ptr_base)

        dna.top(cur_base.id).strand  = n_strand
        dna.top(cur_base.id).address = n_base
        dna.top(cur_base.id).b_14nt  = .false.
        b_visit(cur_base.id)         = .true.

        ! Loop to add a new base into current strand
        do
            ! Check for going out loop
            if(dna.strand(n_strand).b_circular == .false.) then
                ! If non-circular and upper ID is equal to -1
                if(cur_base.up == -1) exit
            else
                ! If circular and base ID is equal to init base ID
                if(cur_base.up == int_base) exit
            end if

            cur_base = dna.top(cur_base.up)

            if(b_visit(cur_base.id)) then
                write(0, "(a$)"), "Error - Reached a visited base : "
                write(0, "(a )"), "SeqDesign_Rebuild_Strand"
                stop
            end if

            n_base = n_base + 1

            ! Insert data to linked list
            allocate(ptr_base)
            ptr_base%id = cur_base.id
            list_base => List_Insert_Base(list_base, ptr_base)

            dna.top(cur_base.id).strand  = n_strand
            dna.top(cur_base.id).address = n_base
            dna.top(cur_base.id).b_14nt  = .false.
            b_visit(cur_base.id)         = .true.
        end do

        dna.strand(n_strand).n_base = n_base

        ! Allocate array using the linked list
        allocate(dna.strand(n_strand).base(n_base))

        ! Put strand data from linked list
        ptr_base => list_base
        do i = 1, n_base
            dna.strand(n_strand).base(n_base+1-i) = ptr_base%id
            ptr_base => ptr_base%next
        end do

    end do

    ! Deallocate b_visit data
    deallocate(b_visit)

    ! Delete linked list allocated
    call List_Delete_Base(list_base)
    !call List_Delete_Base(ptr_base)

    ! Write sequence data
    dna.len_min_stap =  10000
    dna.len_max_stap = -10000
    do i = 1, dna.n_strand
        if(dna.strand(i).types == "stap") then
            if(dna.strand(i).n_base < dna.len_min_stap) dna.len_min_stap = dna.strand(i).n_base
            if(dna.strand(i).n_base > dna.len_max_stap) dna.len_max_stap = dna.strand(i).n_base
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.6. Rebuild strand based data"
        call Space(i, 11)
        write(i, "(a)"), "* Total number of strands                  : "&
            //trim(adjustl(Int2Str(dna.n_strand)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of scaffold strands           : "&
            //trim(adjustl(Int2Str(dna.n_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of staple strands             : "&
            //trim(adjustl(Int2Str(dna.n_stap)))
                call Space(i, 11)
        write(i, "(a)"), "* The minimum staple length                : "&
            //trim(adjustl(Int2Str(dna.len_min_stap)))
                call Space(i, 11)
        write(i, "(a)"), "* The maximum staple length                : "&
            //trim(adjustl(Int2Str(dna.len_max_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* Detailed DNA strand data"
    end do

    ! Print progress in detail
    do i = 1, dna.n_strand
        write(11, "(i20, a$)"), i, " th strand -> "
        write(11, "(a, i7$ )"), "# of bases : ", dna.strand(i).n_base
        write(11, "(a, l, a)"), ", circular : ", dna.strand(i).b_circular, ", bases : "

        ! Print bases numbering
        write(11, "(i30, a$)"), dna.strand(i).base(1), "->"
        do j = 2, dna.strand(i).n_base - 1
            if(mod(j, 100)== 0) then
                write(11, "(i7, a  )"), dna.strand(i).base(j), "->"
            else if(mod(j, 100)== 1) then
                write(11, "(i30, a$)"), dna.strand(i).base(j), "->"
            else
                write(11, "(i7, a$ )"), dna.strand(i).base(j), "->"
            end if
        end do
        write(11, "(i7)"), dna.strand(i).base(dna.strand(i).n_base)
        write(11, "(a )")
    end do

    write(0, "(a)"); write(11, "(a)")
end subroutine SeqDesign_Rebuild_Strand

! ---------------------------------------------------------------------------------------

! List in long length order of the staple
! Last updated on Wednesday 23 November 2016 by Hyungmin
subroutine SeqDesign_Order_Staple(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j

    ! dna.n_stap(i, 1) - original ID
    ! dna.n_stap(i, 2) - short order ID
    allocate(dna.order_stap(dna.n_stap, 2))

    do i = 1, dna.n_stap
        dna.order_stap(i, 1) = i + 1
        dna.order_stap(i, 2) = dna.strand(i+1).n_base
    end do

    ! Sort
    call Sort2(dna.order_stap(:, 2), dna.order_stap(:, 1), dna.n_stap)

    ! Print progress
    do j = 0, 11, 11
        call Space(j, 15)
        write(j, "(a)"), "Short length order of staples"
        call Space(j, 15)
        write(j, "(a)"), "#ID        # original ID     staple length"
        call Space(j, 15)
        write(j, "(a)"), "---        -------------     -------------"
        do i = 1, dna.n_stap
            write(j, "(i18$)"), i
            write(j, "(i18$)"), dna.order_stap(i, 1)
            write(j, "(i18 )"), dna.order_stap(i, 2)
        end do
        write(j, "(a)")
    end do
end subroutine SeqDesign_Order_Staple

! ---------------------------------------------------------------------------------------

! Print 14nt region with various representations
! Last updated on Wednesday 7 November 2016 by Hyungmin
subroutine SeqDesign_Print_14nt_Region(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)
    integer,          allocatable :: region_14nt(:)
    integer,          allocatable :: region_4nt(:)
    integer,          allocatable :: scaf_base2idx(:)
    type(GraphType)               :: graph

    integer :: n_nt_14nt, n_nt_4nt, cn_tn, n_win, n_14nt, n_4nt, n_only_4nt, n_sec_14nt, n_sec_4nt
    integer :: len_ave_stap, one_14nt, two_14nt, three_14nt, other_14nt, n_count(30)
    integer :: i, j, k, base, n_region, n_edge, tot_14nt, tot_4nt
    integer :: across1, across2, node, iniL
    integer :: cir_rgb(3), cir_size
    logical :: b_14nt, b_4nt
    character(200) :: path

    ! ==================================================
    !
    ! Intialize and setup the data structures
    !
    ! ==================================================
    !
    n_14nt   = 0; n_4nt    = 0; n_only_4nt = 0
    one_14nt = 0; two_14nt = 0; three_14nt = 0; other_14nt = 0

    n_nt_14nt = 0; n_nt_4nt = 0
    len_ave_stap = 0

    ! Initialize the variable for couting regions
    allocate(region_14nt(dna.n_strand))
    allocate(region_4nt(dna.n_strand))
    region_14nt(:) = 0
    region_4nt(:)  = 0
    n_count(1:30)  = 0

    ! Count total edges for circular graph
    n_edge = SeqDesign_CirGraph_Count_Edge(dna)

    ! Allocate and initialize variables for the circular graph
    call SeqDesign_CirGraph_Init_Variable(dna, graph, n_edge)

    ! Build scaffold index for bar and circular graph
    allocate(scaf_base2idx(dna.n_base_scaf))
    base = Mani_Go_Start_Base(dna, 1)
    do i = 1, dna.n_base_scaf

        scaf_base2idx(base)   = i       ! From base to index
        graph.base2node(base) = i       ! From base to index
        graph.node2base(i)    = base    ! From index to base

        ! Node type depending on edge number
        if(dna.top(base).node == -1) then
            graph.node(base) = 0
        else
            graph.node(base) = mesh.node(dna.top(base).node).iniL
        end if

        base = dna.top(base).up
    end do

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=791, file=trim(path)//"_count_staple.dat",  form="formatted")
    open(unit=792, file=trim(path)//"_count_region.txt",  form="formatted")
    open(unit=793, file=trim(path)//"_graph_bar.txt",     form="formatted")
    open(unit=794, file=trim(path)//"_graph_cir.graphml", form="formatted")

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "6.7. Print 14nt seed with various representations"
        call Space(i, 11)
        write(i, "(a)"), "* The number of staple strands             : "&
            //trim(adjustl(Int2Str(dna.n_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of total regions              : "&
            //trim(adjustl(Int2Str(ubound(graph.edge, 1))))
        call Space(i, 11)
        write(i, "(a)"), "* Detailed strand region data"
        write(i, "(a)")
    end do

    ! ==================================================
    !
    ! Calculate regions and # of staples
    !
    ! ==================================================
    n_edge = 0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") then
            ! ==================================================
            ! For scaffold strand in bar graph - Graph #3
            ! ==================================================
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                node = dna.top(base).node
                if(node /= -1) then
                    iniL = mesh.node(node).iniL
                else
                    iniL = 0
                end if
                write(793, "(i$)"), i - 1
                write(793, "(f$)"), dble(scaf_base2idx(base)) - 0.5d0
                write(793, "(f$)"), dble(scaf_base2idx(base)) + 0.5d0
                write(793, "( i)"), iniL + 2
            end do
            cycle
        end if

        ! Count staple length
        len_ave_stap = len_ave_stap + dna.strand(i).n_base

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! Build region of staple strands
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! ==================================================
        ! Count strand length - Graph #1
        ! ==================================================
        write(791, "(3i$)"), i - 1, dna.strand(i).n_base, i - 1
        if(n_region > 20) stop
        do j = 1, 20
            if(j <= n_region) then
                if(j /= 1) then
                    if(region(j).sta_pos - region(j-1).end_pos /= 3) then
                        if(region(j).types == 1 .and. region(j-1).types == 1) then
                            write(791, "(i$)"), region(j).sta_pos - region(j-1).end_pos - 1
                        else
                            write(791, "(i$)"), region(j).sta_pos - region(j-1).end_pos - 3
                        end if
                    end if
                end if
                if(region(j).types == 1) write(791, "(i$)"), region(j).length + 1
                if(region(j).types == 2) write(791, "(i$)"), region(j).length + 2
            else
                write(791, "(a$)"), " -- "
            end if
        end do
        write(791, "(a$)")

        ! ==================================================
        ! Count 14nt regions - Graph #2
        ! ==================================================
        n_sec_14nt = 0;   n_sec_4nt  = 0
        b_14nt = .false.; b_4nt= .false.

        do j = 1, n_region
            if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then

                ! Just check one 14nt seeds
                if(b_14nt == .false.) then
                    n_14nt = n_14nt + 1
                    b_14nt = .true.
                end if

                ! All 14nt seeds in the region
                n_sec_14nt = n_sec_14nt + 1
                
                ! Count # of nucleotides in 14nt seed domains
                if(region(j).types == 1) n_nt_14nt = n_nt_14nt + region(j).length + 1
                if(region(j).types == 2) n_nt_14nt = n_nt_14nt + region(j).length + 2
            end if

            ! Global 14nt region couting
            if(region(j).length + 2 <= 30) then     ! There are no region longer than 30bp length
                if(region(j).types == 1) n_count(region(j).length + 1) = n_count(region(j).length + 1) + 1
                if(region(j).types == 2) n_count(region(j).length + 2) = n_count(region(j).length + 2) + 1

                ! Find two crossovers region, 2nt region
                if(j == 1) cycle 
                if(dna.top(dna.top(region(j  ).sta_base).dn).across == -1) cycle
                if(dna.top(dna.top(region(j-1).end_base).up).across == -1) cycle

                n_win = (region(j).sta_pos -1) - (region(j-1).end_pos + 1) - 1
                if(n_win >= 1) n_count(n_win) = n_count(n_win) + 1
            end if
        end do

        ! Count 4nt regions
        do j = 1, n_region
            if( (region(j).types == 1 .and. region(j).length + 1 <= 4) .or. &
                (region(j).types == 2 .and. region(j).length + 2 <= 4) ) then

                ! Just check one 4nt seeds
                if(b_4nt == .false.) then
                    n_4nt = n_4nt + 1
                    b_4nt = .true.
                end if

                ! # of 4nt seeds in the region
                n_sec_4nt = n_sec_4nt + 1

                ! Count # of nucleotides in 4nt seed domains
                if(region(j).types == 1) n_nt_4nt = n_nt_4nt + region(j).length + 1
                if(region(j).types == 2) n_nt_4nt = n_nt_4nt + region(j).length + 2
            end if
        end do

        if(n_sec_14nt == 1) region_14nt(i) = 1
        if(n_sec_14nt == 2) region_14nt(i) = 2
        if(n_sec_14nt == 3) region_14nt(i) = 3
        if(n_sec_14nt == 4) region_14nt(i) = 4

        if(n_sec_4nt  == 1) region_4nt(i)  = 1
        if(n_sec_4nt  == 2) region_4nt(i)  = 2
        if(n_sec_4nt  == 3) region_4nt(i)  = 3
        if(n_sec_4nt  == 4) region_4nt(i)  = 4

        ! ==================================================
        ! For staple strand in bar graph - Graph #3
        ! ==================================================
        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            across1 = dna.top(base).across

            if(across1 /= -1) then
                write(793, "(i$)"), i - 1

                ! To change order
                !if(i - 1 == 96) write(793, "(i$)"), 1
                !if(i - 1 == 95) write(793, "(i$)"), 2
                !if(i - 1 == 90) write(793, "(i$)"), 3

                write(793, "(f$)"), dble(scaf_base2idx(across1)) - 0.5d0
                write(793, "(f$)"), dble(scaf_base2idx(across1)) + 0.5d0

                if(dna.top(base).b_14nt == .false.) write(793, "(i )"), 1
                if(dna.top(base).b_14nt == .true.)  write(793, "(i )"), 2
            end if

            base = dna.top(base).up
        end do

        ! ==================================================
        ! Build node and edge for circular graph - Graph #4
        ! ==================================================
        !
        n_sec_14nt = 0
        do j = 1, n_region

            ! Set node type with different color
            base = dna.top(region(j).cen_base).across

            ! Set node type depending on region size
            if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then

                n_sec_14nt = n_sec_14nt + 1
                !if(n_sec_14nt == 1) graph.node(base) = 1
                !if(n_sec_14nt == 2) graph.node(base) = 2
            else if( (region(j).types == 1 .and. region(j).length + 1 <= 4) .or. &
                     (region(j).types == 2 .and. region(j).length + 2 <= 4) ) then
                !graph.node(base) = 3
            else
                !graph.node(base) = 0
            end if

            ! Set Edge data
            if(j /= n_region) then

                n_edge = n_edge + 1

                ! Set edge in graph
                graph.edge(n_edge, 1) = dna.top(region(j  ).cen_base).across
                graph.edge(n_edge, 2) = dna.top(region(j+1).cen_base).across
            end if
        end do

        ! Update variables
        dna.strand(i).n_14nt = n_sec_14nt
        dna.strand(i).n_4nt  = n_sec_4nt
        dna.len_ave_stap     = dble(len_ave_stap)/dble(dna.n_stap)

        ! Check 14nt only region
        if(b_14nt == .false. .and. b_4nt == .true.) then
            n_only_4nt = n_only_4nt + 1
        end if

        ! Print information on region of staple strands
        do k = 0, 11, 11

            write(k, "(i$    )"), i
            write(k, "(a, i4$)"), " - th staple, # of total bases : ", dna.strand(i).n_base
            write(k, "(a, i2$)"), ", # of bases in poly T : ", cn_tn
            write(k, "(a, l$ )"), ", 14nt region("//trim(adjustl(Int2Str(n_sec_14nt)))//") - ", b_14nt
            write(k, "(a, l  )"), ", 4nt region("//trim(adjustl(Int2Str(n_sec_4nt)))//") - ", b_4nt

            do j = 1, n_region
                write(k, "(i20, a$)"), j, " - th region"

                if(region(j).types == 1) write(k, "(a$)"), " [vertex], "
                if(region(j).types == 2) write(k, "(a$)"), " [ edge ], "

                if(region(j).types == 1) write(k, "(a, 2(i2, a)$)"), "length : ", region(j).length, "(", region(j).length+1, ")"
                if(region(j).types == 2) write(k, "(a, 2(i2, a)$)"), "length : ", region(j).length, "(", region(j).length+2, ")"

                write(k, "(a, i3$)"), ", start pos : ", region(j).sta_pos
                write(k, "(a, i3 )"), ", end pos : ",   region(j).end_pos
            end do
            write(k, "(a)")
        end do

        ! Deallocate memory
        deallocate(region)
    end do

    ! Update global nucleotides in 14nt and 4nt dsDNA domains
    dna.n_nt_14nt = n_nt_14nt
    dna.n_nt_4nt  = n_nt_4nt

    ! ==================================================
    ! Print 14nt regions - Graph #2
    ! ==================================================
    !
    tot_4nt  = 0
    tot_14nt = 0
    do i = 1, dna.n_strand
        if(region_14nt(i) == 1) one_14nt   = one_14nt   + 1
        if(region_14nt(i) == 2) two_14nt   = two_14nt   + 1
        if(region_14nt(i) == 3) three_14nt = three_14nt + 1
        if(region_14nt(i) == 4) other_14nt = other_14nt + 1

        tot_14nt = tot_14nt + region_14nt(i)
        tot_4nt  = tot_4nt  + region_4nt(i)
    end do

    ! Update n_14nt and n_4nt seeds in dna data
    dna.n_14nt       = n_14nt
    dna.n_s14nt      = two_14nt
    dna.n_4nt        = n_4nt
    dna.n_only_4nt   = n_only_4nt
    dna.n_tot_region = n_edge
    dna.n_tot_14nt   = tot_14nt
    dna.n_tot_4nt    = tot_4nt

    ! Print progress
    do i = 0, 11, 11

        call Space(i, 11)
        write(i, "(a)"), "* The number of staples with 14nt regions  : "&
            //trim(adjustl(Int2Str(n_14nt)))
        call Space(i, 16)
        write(i, "(a)"), "* One   14nt seeds : "//trim(adjustl(Int2Str(one_14nt)))
        call Space(i, 16)
        write(i, "(a)"), "* Two   14nt seeds : "//trim(adjustl(Int2Str(two_14nt)))
        call Space(i, 16)
        write(i, "(a)"), "* Three 14nt seeds : "//trim(adjustl(Int2Str(three_14nt)))
        call Space(i, 16)
        write(i, "(a)"), "* More  14nt seeds : "//trim(adjustl(Int2Str(other_14nt)))
        call Space(i, 11)
        write(i, "(a)")

        call Space(i, 11)
        write(i, "(a)"), "* Length of dsDNA domain formed with scaffold (bp)"
        do k = 1, 25
            call Space(i, 16)
            write(i, "(a, i3, a)"), "* ", k, &
                " bp dsDNA domain : "//trim(adjustl(Int2Str(n_count(k))))
        end do
        write(i, "(a)")
    end do

    ! Write text file
    do k = 1, 25
        if(k == 1) then
            write(792, "(6i)"), k, n_count(k), 4, tot_4nt, 14, one_14nt
        else if(k == 2) then
            write(792, "(6i)"), k, n_count(k), 14, tot_14nt, 1414, two_14nt
        else if(k == 3) then
            write(792, "(2i$)"), k, n_count(k)
            write(792, "(2i$)"), 0, n_count(1) + n_count(2) + n_edge - (tot_4nt+tot_14nt)
            write(792, "(2i )"), 0, dna.n_stap - (one_14nt+two_14nt+three_14nt+other_14nt)
        else
            write(792, "(6i)"), k, n_count(k), 0, 0, 0, 0
        end if
    end do

    ! ==================================================
    ! Print circle graph - Graph #4
    ! ==================================================
    !
    ! Write Graphml Format
    write(794, "(a)"), '<?xml version="1.0" encoding="UTF-8"?>'
    write(794, "(a)"), '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">'
    write(794, "(a)"), '<key attr.name="label" attr.type="string" for="node" id="label"/>'
    write(794, "(a)"), '<key attr.name="Edge Label" attr.type="string" for="edge" id="edgelabel"/>'
    write(794, "(a)"), '<key attr.name="weight" attr.type="double" for="edge" id="weight"/>'
    write(794, "(a)"), '<key attr.name="r" attr.type="int" for="node" id="r"/>'
    write(794, "(a)"), '<key attr.name="g" attr.type="int" for="node" id="g"/>'
    write(794, "(a)"), '<key attr.name="b" attr.type="int" for="node" id="b"/>'
    write(794, "(a)"), '<key attr.name="x" attr.type="float" for="node" id="x"/>'
    write(794, "(a)"), '<key attr.name="y" attr.type="float" for="node" id="y"/>'
    write(794, "(a)"), '<key attr.name="size" attr.type="float" for="node" id="size"/>'
    write(794, "(a)"), '<graph edgedefault="undirected">'

    ! Write node with label
    do i = 1, dna.n_base_scaf

        ! Find base from index number
        base = graph.node2base(i)

        write(794, "(a)"), '<node id="'//trim(adjustl(Int2Str(base)))//'">'

        ! Write label corresponding to index number
        if(i < 10) then
            write(794, "(a$)"), '<data key="label">'//'0000'//trim(adjustl(Int2Str(i)))//'</data>'
        else if(i < 100) then
            write(794, "(a$)"), '<data key="label">'//'000'//trim(adjustl(Int2Str(i)))//'</data>'
        else if(i < 1000) then
            write(794, "(a$)"), '<data key="label">'//'00'//trim(adjustl(Int2Str(i)))//'</data>'
        else if(i < 10000) then
            write(794, "(a$)"), '<data key="label">'//'0'//trim(adjustl(Int2Str(i)))//'</data>'
        end if

        ! Change color and size of nodes
        !if(graph.node(base) == 0) then
        !    cir_size = 10; cir_r = 0; cir_g = 0; cir_b = 0
        !else if(graph.node(base) == 1)then
        !    cir_size = 300; cir_r = 255; cir_g = 0; cir_b = 0
        !else if(graph.node(base) == 2) then
        !    cir_size = 400; cir_r = 0; cir_g = 0; cir_b = 255
        !else if(graph.node(base) == 3) then
        !    cir_size = 200; cir_r = 0; cir_g = 192; cir_b = 0
        !end if

        cir_rgb(:) = [64, 64, 64]
        if(graph.node(base) ==  0) cir_rgb(:) = [  64,  64,  64 ]
        if(graph.node(base) ==  1) cir_rgb(:) = [  28, 117, 188 ]
        if(graph.node(base) ==  2) cir_rgb(:) = [ 247, 148,  30 ]
        if(graph.node(base) ==  3) cir_rgb(:) = [ 190,  30,  45 ]
        if(graph.node(base) ==  4) cir_rgb(:) = [   0, 104,  56 ]
        if(graph.node(base) ==  5) cir_rgb(:) = [ 102,  45, 145 ]
        if(graph.node(base) ==  6) cir_rgb(:) = [ 255, 206,   7 ]
        if(graph.node(base) ==  7) cir_rgb(:) = [   1,   1,   1 ]
        if(graph.node(base) ==  8) cir_rgb(:) = [ 141, 198,  63 ]
        if(graph.node(base) ==  9) cir_rgb(:) = [ 194, 181, 155 ]
        if(graph.node(base) == 10) cir_rgb(:) = [  39, 170, 225 ]
        if(graph.node(base) == 11) cir_rgb(:) = [  46,  49, 146 ]
        if(graph.node(base) == 12) cir_rgb(:) = [ 236,   0, 140 ]
        if(graph.node(base) == 13) cir_rgb(:) = [ 147, 149, 152 ]
        if(graph.node(base) == 14) cir_rgb(:) = [ 139,  94,  60 ]
        if(graph.node(base) == 15) cir_rgb(:) = [   0, 167, 157 ]

        if(graph.node(base) == 16) cir_rgb(:) = [  38,  34,  98 ]
        if(graph.node(base) == 17) cir_rgb(:) = [ 117,  76,  41 ] 
        if(graph.node(base) == 18) cir_rgb(:) = [ 241,  90,  41 ]
        if(graph.node(base) == 19) cir_rgb(:) = [  43, 182, 115 ]
        if(graph.node(base) == 20) cir_rgb(:) = [ 146,  39, 143 ]
        if(graph.node(base) == 21) cir_rgb(:) = [ 215, 223,  35 ]
        if(graph.node(base) == 22) cir_rgb(:) = [  65,  64,  66 ]
        if(graph.node(base) == 23) cir_rgb(:) = [   0, 148,  68 ]
        if(graph.node(base) == 24) cir_rgb(:) = [ 155, 133, 121 ]
        if(graph.node(base) == 25) cir_rgb(:) = [  33,  64, 154 ]
        if(graph.node(base) == 26) cir_rgb(:) = [ 237,  28,  36 ]
        if(graph.node(base) == 27) cir_rgb(:) = [  89,  74,  66 ]
        if(graph.node(base) == 28) cir_rgb(:) = [ 188, 190, 192 ]
        if(graph.node(base) == 29) cir_rgb(:) = [ 169, 124,  80 ]
        if(graph.node(base) == 30) cir_rgb(:) = [ 239,  65,  54 ]

        cir_size = 100
        !if(dna.top(base).node == -1) then
        !    cir_size = 1
        !else
        !    cir_size = (mesh.node(dna.top(base).node).sec + 1)
        !end if

        write(794, "(a)"), '<data key="size">'//trim(adjustl(Int2Str(cir_size)))//'.0</data>'
        write(794, "(a)"), '<data key="r">'   //trim(adjustl(Int2Str(cir_rgb(1))))//'</data>'
        write(794, "(a)"), '<data key="g">'   //trim(adjustl(Int2Str(cir_rgb(2))))//'</data>'
        write(794, "(a)"), '<data key="b">'   //trim(adjustl(Int2Str(cir_rgb(3))))//'</data>'
        write(794, "(a)"), '<data key="x">0.0</data>'
        write(794, "(a)"), '<data key="y">0.0</data>'
        write(794, "(a)"), '</node>'
    end do

    ! Write edge
    do i = 1, n_edge
        write(794, "(a$)"), '<edge id="'//trim(adjustl(Int2Str(i-1)))

        if( graph.base2node(graph.edge(i, 1)) > graph.base2node(graph.edge(i, 2))) then
            ! If the source number is larger than target number
            if( graph.base2node(graph.edge(i, 1)) - graph.base2node(graph.edge(i, 2)) > ubound(graph.node, 1) / 2) then
                write(794, "(a$)"), '" source="'//trim(adjustl(Int2Str(graph.edge(i, 2))))
                write(794, "(a )"), '" target="'//trim(adjustl(Int2Str(graph.edge(i, 1))))//'">'
            else
                write(794, "(a$)"), '" source="'//trim(adjustl(Int2Str(graph.edge(i, 1))))
                write(794, "(a )"), '" target="'//trim(adjustl(Int2Str(graph.edge(i, 2))))//'">'
            end if
        else
            ! Normal direction, source < target
            if( graph.base2node(graph.edge(i, 2)) - graph.base2node(graph.edge(i, 1)) > ubound(graph.node, 1) / 2) then
                write(794, "(a$)"), '" source="'//trim(adjustl(Int2Str(graph.edge(i, 1))))
                write(794, "(a )"), '" target="'//trim(adjustl(Int2Str(graph.edge(i, 2))))//'">'
            else
                write(794, "(a$)"), '" source="'//trim(adjustl(Int2Str(graph.edge(i, 2))))
                write(794, "(a )"), '" target="'//trim(adjustl(Int2Str(graph.edge(i, 1))))//'">'
            end if
        end if

        ! Check whether it has 14nt or 4nt and change the thickness of edges
        if(graph.node(graph.edge(i, 1)) == 1 .or. graph.node(graph.edge(i, 2)) == 1) then
            write(794, "(a)"), '<data key="weight">3.0</data>'
        else if(graph.node(graph.edge(i, 1)) == 2 .or. graph.node(graph.edge(i, 2)) == 2) then
            write(794, "(a)"), '<data key="weight">4.0</data>'
        else if(graph.node(graph.edge(i, 1)) == 3 .or. graph.node(graph.edge(i, 2)) == 3) then
            write(794, "(a)"), '<data key="weight">2.0</data>'
        else
            write(794, "(a)"), '<data key="weight">1.0</data>'
        end if
        write(794, "(a)"), '</edge>'
    end do

    write(794, "(a)"), '</graph>'
    write(794, "(a)"), '</graphml>'

    ! Close the file stream
    close(unit=791)
    close(unit=792)
    close(unit=793)
    close(unit=794)

    ! Deallocate memory
    deallocate(region_14nt)
    deallocate(region_4nt)
    deallocate(scaf_base2idx)
    deallocate(graph.node2base)
    deallocate(graph.base2node)
    deallocate(graph.node)
    deallocate(graph.edge)
end subroutine SeqDesign_Print_14nt_Region

! ---------------------------------------------------------------------------------------

! Count the total number of regions for all staples
! Last updated on Tuesday 8 November 2016 by Hyungmin
function SeqDesign_CirGraph_Count_Edge(dna) result(n_edge)
    type(DNAType), intent(inout) :: dna

    type(RegionType), allocatable :: region(:)
    integer :: i, n_region, n_edge

    ! Count the number of edges
    n_edge = 0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types == "scaf") cycle

        ! Find region
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Increse edge number
        n_edge = n_edge + (n_region - 1)
        deallocate(region)
    end do
end function SeqDesign_CirGraph_Count_Edge

! ---------------------------------------------------------------------------------------

! Allocate and initialize variables for the circular graph
! Last updated on Tuesday 8 November 2016 by Hyungmin
subroutine SeqDesign_CirGraph_Init_Variable(dna, graph, n_edge)
    type(DNAType),   intent(in)    :: dna
    type(GraphType), intent(inout) :: graph
    integer,         intent(in)    :: n_edge

    allocate(graph.node2base(dna.n_base_scaf))
    allocate(graph.base2node(dna.n_base_scaf))
    allocate(graph.node(dna.n_base_scaf))
    allocate(graph.edge(n_edge, 3))

    graph.node2base(1:dna.n_base_scaf) = 0
    graph.base2node(1:dna.n_base_scaf) = 0
    graph.node(1:dna.n_base_scaf)      = 0
    graph.edge(1:n_edge, 1) = 0
    graph.edge(1:n_edge, 2) = 0
    graph.edge(1:n_edge, 3) = 0
end subroutine SeqDesign_CirGraph_Init_Variable

! ---------------------------------------------------------------------------------------

! Assign DNA sequence according to para_set_seq_scaf
! Last updated on Tuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Assign_Sequence(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "6.7. Set staple sequence"
        call Space(i, 11)
        write(i, "(a$)"), "* Scaffold sequence                        : "
        if(para_set_seq_scaf == 0) then
            write(i, "(a)"), "M13mp18(7249nt) sequence"
        else if(para_set_seq_scaf == 1) then
            write(i, "(a)"), "user-defined sequence"
        else if(para_set_seq_scaf == 2) then
            write(i, "(a)"), "random sequence"
        end if
        call Space(i, 11)
        write(i, "(a)"), "* The number of scaffold strands           : "&
            //trim(adjustl(Int2Str(dna.n_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of nucleotides in scaffold    : "&
            //trim(adjustl(Int2Str(dna.n_base_scaf)))
    end do

    if(para_set_seq_scaf == 0) then

        ! Set M13mp18 DNA sequence
        call SeqDesign_Set_M13mp18(dna)
    else if(para_set_seq_scaf == 1) then

        ! Import sequence from txt file
        call SeqDesign_Import_Sequence(dna)
    else if(para_set_seq_scaf == 2) then

        ! Set random sequence
        call SeqDesign_Set_Rand_Sequence(dna)
    else
        write(0, "(a$)"), "Error - Not assign para_set_seq_scaf properly : "
        write(0, "(a )"), "SeqDesign_Assign_Sequence"
        stop
    end if

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Detailed bases data builded"
    end do

    ! Print detail information
    do i = 1, dna.n_top
        write(11, "(i20, a$)"), dna.top(i).id, " th base ->"
        write(11, "(a,  i7$)"), "  up : ",         dna.top(i).up
        write(11, "(a,  i7$)"), ", down : ",       dna.top(i).dn
        write(11, "(a,  i7$)"), ", across : ",     dna.top(i).across
        write(11, "(a,  i7$)"), ", xover : ",      dna.top(i).xover
        write(11, "(a,  i7$)"), ", node : ",       dna.top(i).node
        write(11, "(a,  a8$)"), ", seq : ", dna.top(i).seq
        write(11, "(a,  i7$)"), ", strand : ",     dna.top(i).strand
        write(11, "(a,  i7 )"), ", address : ",    dna.top(i).address
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine SeqDesign_Assign_Sequence

! ---------------------------------------------------------------------------------------

! Set M13mp18 DNA sequence
! Last updated on Tuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Set_M13mp18(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j, count, n_base_scaf, base, across
    integer, parameter :: len_M13 = 7249
    character(len_M13) :: M13_seq

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(0, "(a$)"), "Error - The number of bases in scaffolds are not consistent : "
        write(0, "(a )"), "SeqDesign_Set_M13mp18"
        stop
    end if

    ! Set M13mp18 sequence
    M13_seq = SeqDesign_Get_M13mp18(len_M13)

    ! Set scaffold sequence as full M13 scaffold
    do i = 1, dna.n_strand

        ! Set sequence for scaffold
        if(dna.strand(i).types /= "scaf") cycle

        ! Find the starting point (down == -1)
        base   = Mani_Go_Start_Base(dna, i)
        across = dna.top(base).across

        do j = 1, dna.strand(i).n_base

            count = j+para_set_start_scaf-1
            if(count >= len_M13) then
                count = mod(count, len_M13)
                if(count == 0) count = len_M13
            end if

            ! Assign M13mp18 sequence
            dna.top(base).seq = M13_seq(count:count)

            ! Set complementary sequence
            if(across /= -1) then
                dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
            end if

            ! Update base
            if(j /= dna.strand(i).n_base) then
                base   = dna.top(base).up
                across = dna.top(base).across
            end if
        end do
    end do
end subroutine SeqDesign_Set_M13mp18

! ---------------------------------------------------------------------------------------

! Get M13mp18 DNA sequence
! Last updated on Monday 28 Mar 2016 by Hyungmin
function SeqDesign_Get_M13mp18(len_M13) result(M13_seq)
    integer, intent(in) :: len_M13

    character(len_M13) :: M13_seq, M13_seq_old
    integer :: i, count

    ! Initialize M13 sequence data as "N"
    do i = 1, len_M13
        M13_seq(i:i)     = "N"
        M13_seq_old(i:i) = "N"
    end do

    ! Set scaffold sequence as full M13 scaffold, 7249 nucleotides
    ! M13mp18 [length=7249] [version=09-MAY-2008] [topology=circular] Cloning vector M13mp18
    M13_seq = &
        "AATGCTACTACTATTAGTAGAATTGATGCCACCTTTTCAGCTCGCGCCCC"//&
        "AAATGAAAATATAGCTAAACAGGTTATTGACCATTTGCGAAATGTATCTA"//&   !    1 -  100
        "ATGGTCAAACTAAATCTACTCGTTCGCAGAATTGGGAATCAACTGTTATA"//&
        "TGGAATGAAACTTCCAGACACCGTACTTTAGTTGCATATTTAAAACATGT"//&   !  101 -  200
        "TGAGCTACAGCATTATATTCAGCAATTAAGCTCTAAGCCATCCGCAAAAA"//&
        "TGACCTCTTATCAAAAGGAGCAATTAAAGGTACTCTCTAATCCTGACCTG"//&   !  201 -  300
        "TTGGAGTTTGCTTCCGGTCTGGTTCGCTTTGAAGCTCGAATTAAAACGCG"//&
        "ATATTTGAAGTCTTTCGGGCTTCCTCTTAATCTTTTTGATGCAATCCGCT"//&   !  301 -  400

        "TTGCTTCTGACTATAATAGTCAGGGTAAAGACCTGATTTTTGATTTATGG"//&
        "TCATTCTCGTTTTCTGAACTGTTTAAAGCATTTGAGGGGGATTCAATGAA"//&   !  401 -  500
        "TATTTATGACGATTCCGCAGTATTGGACGCTATCCAGTCTAAACATTTTA"//&
        "CTATTACCCCCTCTGGCAAAACTTCTTTTGCAAAAGCCTCTCGCTATTTT"//&   !  501 -  600
        "GGTTTTTATCGTCGTCTGGTAAACGAGGGTTATGATAGTGTTGCTCTTAC"//&
        "TATGCCTCGTAATTCCTTTTGGCGTTATGTATCTGCATTAGTTGAATGTG"//&   !  601 -  700
        "GTATTCCTAAATCTCAACTGATGAATCTTTCTACCTGTAATAATGTTGTT"//&
        "CCGTTAGTTCGTTTTATTAACGTAGATTTTTCTTCCCAACGTCCTGACTG"//&   !  701 -  800

        "GTATAATGAGCCAGTTCTTAAAATCGCATAAGGTAATTCACAATGATTAA"//&
        "AGTTGAAATTAAACCATCTCAAGCCCAATTTACTACTCGTTCTGGTGTTT"//&   !  801 -  900
        "CTCGTCAGGGCAAGCCTTATTCACTGAATGAGCAGCTTTGTTACGTTGAT"//&
        "TTGGGTAATGAATATCCGGTTCTTGTCAAGATTACTCTTGATGAAGGTCA"//&   !  901 - 1000
        "GCCAGCCTATGCGCCTGGTCTGTACACCGTTCATCTGTCCTCTTTCAAAG"//&
        "TTGGTCAGTTCGGTTCCCTTATGATTGACCGTCTGCGCCTCGTTCCGGCT"//&   ! 1001 - 1100
        "AAGTAACATGGAGCAGGTCGCGGATTTCGACACAATTTATCAGGCGATGA"//&
        "TACAAATCTCCGTTGTACTTTGTTTCGCGCTTGGTATAATCGCTGGGGGT"//&   ! 1101 - 1200

        "CAAAGATGAGTGTTTTAGTGTATTCTTTTGCCTCTTTCGTTTTAGGTTGG"//&
        "TGCCTTCGTAGTGGCATTACGTATTTTACCCGTTTAATGGAAACTTCCTC"//&   ! 1201 - 1300
        "ATGAAAAAGTCTTTAGTCCTCAAAGCCTCTGTAGCCGTTGCTACCCTCGT"//&
        "TCCGATGCTGTCTTTCGCTGCTGAGGGTGACGATCCCGCAAAAGCGGCCT"//&   ! 1301 - 1400
        "TTAACTCCCTGCAAGCCTCAGCGACCGAATATATCGGTTATGCGTGGGCG"//&
        "ATGGTTGTTGTCATTGTCGGCGCAACTATCGGTATCAAGCTGTTTAAGAA"//&   ! 1401 - 1500
        "ATTCACCTCGAAAGCAAGCTGATAAACCGATACAATTAAAGGCTCCTTTT"//&
        "GGAGCCTTTTTTTTGGAGATTTTCAACGTGAAAAAATTATTATTCGCAAT"//&   ! 1501 - 1600

        "TCCTTTAGTTGTTCCTTTCTATTCTCACTCCGCTGAAACTGTTGAAAGTT"//&
        "GTTTAGCAAAATCCCATACAGAAAATTCATTTACTAACGTCTGGAAAGAC"//&   ! 1601 - 1700
        "GACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGC"//&
        "TACAGGCGTTGTAGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACAT"//&   ! 1701 - 1800
        "GGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAGGGTGGTGGCTCTGAG"//&
        "GGTGGCGGTTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTACTAAACCTCC"//&   ! 1801 - 1900
        "TGAGTACGGTGATACACCTATTCCGGGCTATACTTATATCAACCCTCTCG"//&
        "ACGGCACTTATCCGCCTGGTACTGAGCAAAACCCCGCTAATCCTAATCCT"//&   ! 1901 - 2000

        "TCTCTTGAGGAGTCTCAGCCTCTTAATACTTTCATGTTTCAGAATAATAG"//&
        "GTTCCGAAATAGGCAGGGGGCATTAACTGTTTATACGGGCACTGTTACTC"//&   ! 2001 - 2100
        "AAGGCACTGACCCCGTTAAAACTTATTACCAGTACACTCCTGTATCATCA"//&
        "AAAGCCATGTATGACGCTTACTGGAACGGTAAATTCAGAGACTGCGCTTT"//&   ! 2101 - 2200
        "CCATTCTGGCTTTAATGAGGATTTATTTGTTTGTGAATATCAAGGCCAAT"//&
        "CGTCTGACCTGCCTCAACCTCCTGTCAATGCTGGCGGCGGCTCTGGTGGT"//&   ! 2201 - 2300
        "GGTTCTGGTGGCGGCTCTGAGGGTGGTGGCTCTGAGGGTGGCGGTTCTGA"//&
        "GGGTGGCGGCTCTGAGGGAGGCGGTTCCGGTGGTGGCTCTGGTTCCGGTG"//&   ! 2301 - 2400

        "ATTTTGATTATGAAAAGATGGCAAACGCTAATAAGGGGGCTATGACCGAA"//&
        "AATGCCGATGAAAACGCGCTACAGTCTGACGCTAAAGGCAAACTTGATTC"//&   ! 2401 - 2500
        "TGTCGCTACTGATTACGGTGCTGCTATCGATGGTTTCATTGGTGACGTTT"//&
        "CCGGCCTTGCTAATGGTAATGGTGCTACTGGTGATTTTGCTGGCTCTAAT"//&   ! 2501 - 2600
        "TCCCAAATGGCTCAAGTCGGTGACGGTGATAATTCACCTTTAATGAATAA"//&
        "TTTCCGTCAATATTTACCTTCCCTCCCTCAATCGGTTGAATGTCGCCCTT"//&   ! 2601 - 2700
        "TTGTCTTTGGCGCTGGTAAACCATATGAATTTTCTATTGATTGTGACAAA"//&
        "ATAAACTTATTCCGTGGTGTCTTTGCGTTTCTTTTATATGTTGCCACCTT"//&   ! 2701 - 2800

        "TATGTATGTATTTTCTACGTTTGCTAACATACTGCGTAATAAGGAGTCTT"//&
        "AATCATGCCAGTTCTTTTGGGTATTCCGTTATTATTGCGTTTCCTCGGTT"//&   ! 2801 - 2900
        "TCCTTCTGGTAACTTTGTTCGGCTATCTGCTTACTTTTCTTAAAAAGGGC"//&
        "TTCGGTAAGATAGCTATTGCTATTTCATTGTTTCTTGCTCTTATTATTGG"//&   ! 2901 - 3000
        "GCTTAACTCAATTCTTGTGGGTTATCTCTCTGATATTAGCGCTCAATTAC"//&
        "CCTCTGACTTTGTTCAGGGTGTTCAGTTAATTCTCCCGTCTAATGCGCTT"//&   ! 3001 - 3100
        "CCCTGTTTTTATGTTATTCTCTCTGTAAAGGCTGCTATTTTCATTTTTGA"//&
        "CGTTAAACAAAAAATCGTTTCTTATTTGGATTGGGATAAATAATATGGCT"//&   ! 3101 - 3200

        "GTTTATTTTGTAACTGGCAAATTAGGCTCTGGAAAGACGCTCGTTAGCGT"//&
        "TGGTAAGATTCAGGATAAAATTGTAGCTGGGTGCAAAATAGCAACTAATC"//&   ! 3201 - 3300
        "TTGATTTAAGGCTTCAAAACCTCCCGCAAGTCGGGAGGTTCGCTAAAACG"//&
        "CCTCGCGTTCTTAGAATACCGGATAAGCCTTCTATATCTGATTTGCTTGC"//&   ! 3301 - 3400
        "TATTGGGCGCGGTAATGATTCCTACGATGAAAATAAAAACGGCTTGCTTG"//&
        "TTCTCGATGAGTGCGGTACTTGGTTTAATACCCGTTCTTGGAATGATAAG"//&   ! 3401 - 3500
        "GAAAGACAGCCGATTATTGATTGGTTTCTACATGCTCGTAAATTAGGATG"//&
        "GGATATTATTTTTCTTGTTCAGGACTTATCTATTGTTGATAAACAGGCGC"//&   ! 3501 - 3600

        "GTTCTGCATTAGCTGAACATGTTGTTTATTGTCGTCGTCTGGACAGAATT"//&
        "ACTTTACCTTTTGTCGGTACTTTATATTCTCTTATTACTGGCTCGAAAAT"//&   ! 3601 - 3700
        "GCCTCTGCCTAAATTACATGTTGGCGTTGTTAAATATGGCGATTCTCAAT"//&
        "TAAGCCCTACTGTTGAGCGTTGGCTTTATACTGGTAAGAATTTGTATAAC"//&   ! 3701 - 3800
        "GCATATGATACTAAACAGGCTTTTTCTAGTAATTATGATTCCGGTGTTTA"//&
        "TTCTTATTTAACGCCTTATTTATCACACGGTCGGTATTTCAAACCATTAA"//&   ! 3801 - 3900
        "ATTTAGGTCAGAAGATGAAATTAACTAAAATATATTTGAAAAAGTTTTCT"//&
        "CGCGTTCTTTGTCTTGCGATTGGATTTGCATCAGCATTTACATATAGTTA"//&   ! 3901 - 4000

        "TATAACCCAACCTAAGCCGGAGGTTAAAAAGGTAGTCTCTCAGACCTATG"//&
        "ATTTTGATAAATTCACTATTGACTCTTCTCAGCGTCTTAATCTAAGCTAT"//&   ! 4001 - 4100
        "CGCTATGTTTTCAAGGATTCTAAGGGAAAATTAATTAATAGCGACGATTT"//&
        "ACAGAAGCAAGGTTATTCACTCACATATATTGATTTATGTACTGTTTCCA"//&   ! 4101 - 4200
        "TTAAAAAAGGTAATTCAAATGAAATTGTTAAATGTAATTAATTTTGTTTT"//&
        "CTTGATGTTTGTTTCATCATCTTCTTTTGCTCAGGTAATTGAAATGAATA"//&   ! 4201 - 4300
        "ATTCGCCTCTGCGCGATTTTGTAACTTGGTATTCAAAGCAATCAGGCGAA"//&
        "TCCGTTATTGTTTCTCCCGATGTAAAAGGTACTGTTACTGTATATTCATC"//&   ! 4301 - 4400

        "TGACGTTAAACCTGAAAATCTACGCAATTTCTTTATTTCTGTTTTACGTG"//&
        "CAAATAATTTTGATATGGTAGGTTCTAACCCTTCCATTATTCAGAAGTAT"//&   ! 4401 - 4500
        "AATCCAAACAATCAGGATTATATTGATGAATTGCCATCATCTGATAATCA"//&
        "GGAATATGATGATAATTCCGCTCCTTCTGGTGGTTTCTTTGTTCCGCAAA"//&   ! 4501 - 4600
        "ATGATAATGTTACTCAAACTTTTAAAATTAATAACGTTCGGGCAAAGGAT"//&
        "TTAATACGAGTTGTCGAATTGTTTGTAAAGTCTAATACTTCTAAATCCTC"//&   ! 4601 - 4700
        "AAATGTATTATCTATTGACGGCTCTAATCTATTAGTTGTTAGTGCTCCTA"//&
        "AAGATATTTTAGATAACCTTCCTCAATTCCTTTCAACTGTTGATTTGCCA"//&   ! 4701 - 4800

        "ACTGACCAGATATTGATTGAGGGTTTGATATTTGAGGTTCAGCAAGGTGA"//&
        "TGCTTTAGATTTTTCATTTGCTGCTGGCTCTCAGCGTGGCACTGTTGCAG"//&   ! 4801 - 4900
        "GCGGTGTTAATACTGACCGCCTCACCTCTGTTTTATCTTCTGCTGGTGGT"//&
        "TCGTTCGGTATTTTTAATGGCGATGTTTTAGGGCTATCAGTTCGCGCATT"//&   ! 4901 - 5000
        "AAAGACTAATAGCCATTCAAAAATATTGTCTGTGCCACGTATTCTTACGC"//&
        "TTTCAGGTCAGAAGGGTTCTATCTCTGTTGGCCAGAATGTCCCTTTTATT"//&   ! 5001 - 5100
        "ACTGGTCGTGTGACTGGTGAATCTGCCAATGTAAATAATCCATTTCAGAC"//&
        "GATTGAGCGTCAAAATGTAGGTATTTCCATGAGCGTTTTTCCTGTTGCAA"//&   ! 5101 - 5200

        "TGGCTGGCGGTAATATTGTTCTGGATATTACCAGCAAGGCCGATAGTTTG"//&
        "AGTTCTTCTACTCAGGCAAGTGATGTTATTACTAATCAAAGAAGTATTGC"//&   ! 5201 - 5300
        "TACAACGGTTAATTTGCGTGATGGACAGACTCTTTTACTCGGTGGCCTCA"//&
        "CTGATTATAAAAACACTTCTCAGGATTCTGGCGTACCGTTCCTGTCTAAA"//&   ! 5301 - 5400
        "ATCCCTTTAATCGGCCTCCTGTTTAGCTCCCGCTCTGATTCTAACGAGGA"//&
        "AAGCACGTTATACGTGCTCGTCAAAGCAACCATAGTACGCGCCCTGTAGC"//&   ! 5401 - 5500
        "GGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTAC"//&
        "ACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTC"//&   ! 5501 - 5600

        "TCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCT"//&
        "TTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGA"//&   ! 5601 - 5700
        "TTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTC"//&
        "GCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAA"//&   ! 5701 - 5800
        "ACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGG"//&
        "GATTTTGCCGATTTCGGAACCACCATCAAACAGGATTTTCGCCTGCTGGG"//&   ! 5801 - 5900
        "GCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGA"//&
        "AGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTG"//&   ! 5901 - 6000

        "GCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAAT"//&
        "GCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAAC"//&   ! 6001 - 6100
        "GCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTT"//&
        "TATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTC"//&   ! 6101 - 6200
        "ACACAGGAAACAGCTATGACCATGATTACGAATTCGAGCTCGGTACCCGG"//&
        "GGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCG"//&   ! 6201 - 6300
        "TTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGC"//&
        "CTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCG"//&   ! 6301 - 6400

        "CACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCT"//&
        "TTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGC"//&   ! 6401 - 6500
        "GATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCA"//&
        "CGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCA"//&   ! 6501 - 6600
        "ATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACA"//&
        "TTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTT"//&   ! 6601 - 6700
        "TGATGGCGTTCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAA"//&
        "TGCGAATTTTAACAAAATATTAACGTTTACAATTTAAATATTTGCTTATA"//&   ! 6701 - 6800

        "CAATCTTCCTGTTTTTGGGGCTTTTCTGATTATCAACCGGGGTACATATG"//&
        "ATTGACATGCTAGTTTTACGATTACCGTTCATCGATTCTCTTGTTTGCTC"//&   ! 6801 - 6900
        "CAGACTCTCAGGCAATGACCTGATAGCCTTTGTAGATCTCTCAAAAATAG"//&
        "CTACCCTCTCCGGCATTAATTTATCAGCTAGAACGGTTGAATATCATATT"//&   ! 6901 - 7000
        "GATGGTGATTTGACTGTCTCCGGCCTTTCTCACCCTTTTGAATCTTTACC"//&
        "TACACATTACTCAGGCATTGCATTTAAAATATATGAGGGTTCTAAAAATT"//&   ! 7001 - 7100
        "TTTATCCTTGCGTTGAAATAAAGGCTTCTCCCGCAAAAGTATTACAGGGT"//&
        "CATAATGTTTTTGGTACAACCGATTTAGCTTTATGCTCTGAGGCTTTATT"//&   ! 7101 - 7200
        "GCTTAATTTTGCTAATTCTTTGCCTTGCCTGTATGATTTATTGGATGTT"       ! 7200 - 7249

    M13_seq_old(   1:  50) = Mani_To_Upper("aatgctactactattagtagaattgatgccaccttttcagctcgcgcccc")
    M13_seq_old(  51: 100) = Mani_To_Upper("aaatgaaaatatagctaaacaggttattgaccatttgcgaaatgtatcta")
    M13_seq_old( 101: 150) = Mani_To_Upper("atggtcaaactaaatctactcgttcgcagaattgggaatcaactgttaca")
    M13_seq_old( 151: 200) = Mani_To_Upper("tggaatgaaacttccagacaccgtactttagttgcatatttaaaacatgt")
    M13_seq_old( 201: 250) = Mani_To_Upper("tgagctacagcaccagattcagcaattaagctctaagccatccgcaaaaa")
    M13_seq_old( 251: 300) = Mani_To_Upper("tgacctcttatcaaaaggagcaattaaaggtactctctaatcctgacctg")
    M13_seq_old( 301: 350) = Mani_To_Upper("ttggagtttgcttccggtctggttcgctttgaagctcgaattaaaacgcg")
    M13_seq_old( 351: 400) = Mani_To_Upper("atatttgaagtctttcgggcttcctcttaatctttttgatgcaatccgct")
    M13_seq_old( 401: 450) = Mani_To_Upper("ttgcttctgactataatagtcagggtaaagacctgatttttgatttatgg")
    M13_seq_old( 451: 500) = Mani_To_Upper("tcattctcgttttctgaactgtttaaagcatttgagggggattcaatgaa")
    M13_seq_old( 501: 550) = Mani_To_Upper("tatttatgacgattccgcagtattggacgctatccagtctaaacatttta")
    M13_seq_old( 551: 600) = Mani_To_Upper("ctattaccccctctggcaaaacttcttttgcaaaagcctctcgctatttt")
    M13_seq_old( 601: 650) = Mani_To_Upper("ggtttttatcgtcgtctggtaaacgagggttatgatagtgttgctcttac")
    M13_seq_old( 651: 700) = Mani_To_Upper("tatgcctcgtaattccttttggcgttatgtatctgcattagttgaatgtg")
    M13_seq_old( 701: 750) = Mani_To_Upper("gtattcctaaatctcaactgatgaatctttctacctgtaataatgttgtt")
    M13_seq_old( 751: 800) = Mani_To_Upper("ccgttagttcgttttattaacgtagatttttcttcccaacgtcctgactg")
    M13_seq_old( 801: 850) = Mani_To_Upper("gtataatgagccagttcttaaaatcgcataaggtaattcacaatgattaa")
    M13_seq_old( 851: 900) = Mani_To_Upper("agttgaaattaaaccatctcaagcccaatttactactcgttctggtgttc")
    M13_seq_old( 901: 950) = Mani_To_Upper("tcgtcagggcaagccttattcactgaatgagcagctttgttacgttgatt")
    M13_seq_old( 951:1000) = Mani_To_Upper("tgggtaatgaatatccggttcttgtcaagattactcttgatgaaggtcag")
    M13_seq_old(1001:1150) = Mani_To_Upper("ccagcctatgcgcctggtctgtacaccgttcatctgtcctctttcaaagt")
    M13_seq_old(1051:1200) = Mani_To_Upper("tggtcagttcggttcccttatgattgaccgtctgcgcctcgttccggcta")
    M13_seq_old(1101:1250) = Mani_To_Upper("agtaacatggagcaggtcgcggatttcgacacaatttatcaggcgatgat")
    M13_seq_old(1151:1200) = Mani_To_Upper("acaaatctccgttgtactttgtttcgcgcttggtataatcgctgggggtc")
    M13_seq_old(1201:1250) = Mani_To_Upper("aaagatgagtgttttagtgtattctttcgcctctttcgttttaggttggt")
    M13_seq_old(1251:1300) = Mani_To_Upper("gccttcgtagtggcattacgtattttacccgtttaatggaaacttcctca")
    M13_seq_old(1301:1350) = Mani_To_Upper("tgaaaaagtctttagtcctcaaagcctctgtagccgttgctaccctcgtt")
    M13_seq_old(1351:1400) = Mani_To_Upper("ccgatgctgtctttcgctgctgagggtgacgatcccgcaaaagcggcctt")
    M13_seq_old(1401:1450) = Mani_To_Upper("taactccctgcaagcctcagcgaccgaatatatcggttatgcgtgggcga")
    M13_seq_old(1451:1500) = Mani_To_Upper("tggttgttgtcattgtcggcgcaactatcggtatcaagctgtttaagaaa")
    M13_seq_old(1501:1550) = Mani_To_Upper("ttcacctcgaaagcaagctgataaaccgatacaattaaaggctccttttg")
    M13_seq_old(1551:1600) = Mani_To_Upper("gagcctttttttttggagattttcaacgtgaaaaaattattattcgcaat")
    M13_seq_old(1601:1650) = Mani_To_Upper("tcctttagttgttcctttctattctcactccgctgaaactgttgaaagtt")
    M13_seq_old(1651:1700) = Mani_To_Upper("gtttagcaaaaccccatacagaaaattcatttactaacgtctggaaagac")
    M13_seq_old(1701:1750) = Mani_To_Upper("gacaaaactttagatcgttacgctaactatgagggttgtctgtggaatgc")
    M13_seq_old(1751:1800) = Mani_To_Upper("tacaggcgttgtagtttgtactggtgacgaaactcagtgttacggtacat")
    M13_seq_old(1801:1850) = Mani_To_Upper("gggttcctattgggcttgctatccctgaaaatgagggtggtggctctgag")
    M13_seq_old(1851:1900) = Mani_To_Upper("ggtggcggttctgagggtggcggttctgagggtggcggtactaaacctcc")
    M13_seq_old(1901:1950) = Mani_To_Upper("tgagtacggtgatacacctattccgggctatacttatatcaaccctctcg")
    M13_seq_old(1951:2000) = Mani_To_Upper("acggcacttatccgcctggtactgagcaaaaccccgctaatcctaatcct")
    M13_seq_old(2001:2050) = Mani_To_Upper("tctcttgaggagtctcagcctcttaatactttcatgtttcagaataatag")
    M13_seq_old(2051:2100) = Mani_To_Upper("gttccgaaataggcagggggcattaactgtttatacgggcactgttactc")
    M13_seq_old(2101:2150) = Mani_To_Upper("aaggcactgaccccgttaaaacttattaccagtacactcctgtatcatca")
    M13_seq_old(2151:2200) = Mani_To_Upper("aaagccatgtatgacgcttactggaacggtaaattcagagactgcgcttt")
    M13_seq_old(2201:2250) = Mani_To_Upper("ccattctggctttaatgaagatccattcgtttgtgaatatcaaggccaat")
    M13_seq_old(2251:2300) = Mani_To_Upper("cgtctgacctgcctcaacctcctgtcaatgctggcggcggctctggtggt")
    M13_seq_old(2301:2350) = Mani_To_Upper("ggttctggtggcggctctgagggtggtggctctgagggtggcggttctga")
    M13_seq_old(2351:2400) = Mani_To_Upper("gggtggcggctctgagggaggcggttccggtggtggctctggttccggtg")
    M13_seq_old(2401:2450) = Mani_To_Upper("attttgattatgaaaagatggcaaacgctaataagggggctatgaccgaa")
    M13_seq_old(2451:2500) = Mani_To_Upper("aatgccgatgaaaacgcgctacagtctgacgctaaaggcaaacttgattc")
    M13_seq_old(2501:2550) = Mani_To_Upper("tgtcgctactgattacggtgctgctatcgatggtttcattggtgacgttt")
    M13_seq_old(2551:2600) = Mani_To_Upper("ccggccttgctaatggtaatggtgctactggtgattttgctggctctaat")
    M13_seq_old(2601:2650) = Mani_To_Upper("tcccaaatggctcaagtcggtgacggtgataattcacctttaatgaataa")
    M13_seq_old(2651:2700) = Mani_To_Upper("tttccgtcaatatttaccttccctccctcaatcggttgaatgtcgccctt")
    M13_seq_old(2701:2750) = Mani_To_Upper("ttgtctttagcgctggtaaaccatatgaattttctattgattgtgacaaa")
    M13_seq_old(2751:2800) = Mani_To_Upper("ataaacttattccgtggtgtctttgcgtttcttttatatgttgccacctt")
    M13_seq_old(2801:2850) = Mani_To_Upper("tatgtatgtattttctacgtttgctaacatactgcgtaataaggagtctt")
    M13_seq_old(2851:2900) = Mani_To_Upper("aatcatgccagttcttttgggtattccgttattattgcgtttcctcggtt")
    M13_seq_old(2901:2950) = Mani_To_Upper("tccttctggtaactttgttcggctatctgcttacttttcttaaaaagggc")
    M13_seq_old(2951:3000) = Mani_To_Upper("ttcggtaagatagctattgctatttcattgtttcttgctcttattattgg")
    M13_seq_old(3001:3050) = Mani_To_Upper("gcttaactcaattcttgtgggttatctctctgatattagcgctcaattac")
    M13_seq_old(3051:3100) = Mani_To_Upper("cctctgactttgttcagggtgttcagttaattctcccgtctaatgcgctt")
    M13_seq_old(3101:3150) = Mani_To_Upper("ccctgtttttatgttattctctctgtaaaggctgctattttcatttttga")
    M13_seq_old(3151:3200) = Mani_To_Upper("cgttaaacaaaaaatcgtttcttatttggattgggataaataatatggct")
    M13_seq_old(3201:3250) = Mani_To_Upper("gtttattttgtaactggcaaattaggctctggaaagacgctcgttagcgt")
    M13_seq_old(3251:3300) = Mani_To_Upper("tggtaagattcaggataaaattgtagctgggtgcaaaatagcaactaatc")
    M13_seq_old(3301:3350) = Mani_To_Upper("ttgatttaaggcttcaaaacctcccgcaagtcgggaggttcgctaaaacg")
    M13_seq_old(3351:3400) = Mani_To_Upper("cctcgcgttcttagaataccggataagccttctatatctgatttgcttgc")
    M13_seq_old(3401:3450) = Mani_To_Upper("tattgggcgcggtaatgattcctacgatgaaaataaaaacggcttgcttg")
    M13_seq_old(3451:3500) = Mani_To_Upper("ttctcgatgagtgcggtacttggtttaatacccgttcttggaatgataag")
    M13_seq_old(3501:3550) = Mani_To_Upper("gaaagacagccgattattgattggtttctacatgctcgtaaattaggatg")
    M13_seq_old(3551:3600) = Mani_To_Upper("ggatattatttttcttgttcaggacttatctattgttgataaacaggcgc")
    M13_seq_old(3601:3650) = Mani_To_Upper("gttctgcattagctgaacatgttgtttattgtcgtcgtctggacagaatt")
    M13_seq_old(3651:3700) = Mani_To_Upper("actttaccttttgtcggtactttatattctcttattactggctcgaaaat")
    M13_seq_old(3701:3750) = Mani_To_Upper("gcctctgcctaaattacatgttggcgttgttaaatatggcgattctcaat")
    M13_seq_old(3751:3800) = Mani_To_Upper("taagccctactgttgagcgttggctttatactggtaagaatttgtataac")
    M13_seq_old(3801:3850) = Mani_To_Upper("gcatatgatactaaacaggctttttctagtaattatgattccggtgttta")
    M13_seq_old(3851:3900) = Mani_To_Upper("ttcttatttaacgccttatttatcacacggtcggtatttcaaaccattaa")
    M13_seq_old(3901:3950) = Mani_To_Upper("atttaggtcagaagatgaaattaactaaaatatatttgaaaaagttttct")
    M13_seq_old(3951:4000) = Mani_To_Upper("cgcgttctttgtcttgcgattggatttgcatcagcatttacatatagtta")
    M13_seq_old(4001:4050) = Mani_To_Upper("tataacccaacctaagccggaggttaaaaaggtagtctctcagacctatg")
    M13_seq_old(4051:4100) = Mani_To_Upper("attttgataaattcactattgactcttctcagcgtcttaatctaagctat")
    M13_seq_old(4101:4150) = Mani_To_Upper("cgctatgttttcaaggattctaagggaaaattaattaatagcgacgattt")
    M13_seq_old(4151:4200) = Mani_To_Upper("acagaagcaaggttattcactcacatatattgatttatgtactgtttcca")
    M13_seq_old(4201:4250) = Mani_To_Upper("ttaaaaaaggtaattcaaatgaaattgttaaatgtaattaattttgtttt")
    M13_seq_old(4251:4300) = Mani_To_Upper("cttgatgtttgtttcatcatcttcttttgctcaggtaattgaaatgaata")
    M13_seq_old(4301:4350) = Mani_To_Upper("attcgcctctgcgcgattttgtaacttggtattcaaagcaatcaggcgaa")
    M13_seq_old(4351:4400) = Mani_To_Upper("tccgttattgtttctcccgatgtaaaaggtactgttactgtatattcatc")
    M13_seq_old(4401:4450) = Mani_To_Upper("tgacgttaaacctgaaaatctacgcaatttctttatttctgttttacgtg")
    M13_seq_old(4451:4500) = Mani_To_Upper("ctaataattttgatatggttggttcaattccttccataattcagaagtat")
    M13_seq_old(4501:4550) = Mani_To_Upper("aatccaaacaatcaggattatattgatgaattgccatcatctgataatca")
    M13_seq_old(4551:4600) = Mani_To_Upper("ggaatatgatgataattccgctccttctggtggtttctttgttccgcaaa")
    M13_seq_old(4601:4650) = Mani_To_Upper("atgataatgttactcaaacttttaaaattaataacgttcgggcaaaggat")
    M13_seq_old(4651:4700) = Mani_To_Upper("ttaatacgagttgtcgaattgtttgtaaagtctaatacttctaaatcctc")
    M13_seq_old(4701:4750) = Mani_To_Upper("aaatgtattatctattgacggctctaatctattagttgttagtgcaccta")
    M13_seq_old(4751:4800) = Mani_To_Upper("aagatattttagataaccttcctcaattcctttctactgttgatttgcca")
    M13_seq_old(4801:4850) = Mani_To_Upper("actgaccagatattgattgagggtttgatatttgaggttcagcaaggtga")
    M13_seq_old(4851:4900) = Mani_To_Upper("tgctttagatttttcatttgctgctggctctcagcgtggcactgttgcag")
    M13_seq_old(4901:4950) = Mani_To_Upper("gcggtgttaatactgaccgcctcacctctgttttatcttctgctggtggt")
    M13_seq_old(4951:5000) = Mani_To_Upper("tcgttcggtatttttaatggcgatgttttagggctatcagttcgcgcatt")
    M13_seq_old(5001:5050) = Mani_To_Upper("aaagactaatagccattcaaaaatattgtctgtgccacgtattcttacgc")
    M13_seq_old(5051:5100) = Mani_To_Upper("tttcaggtcagaagggttctatctctgttggccagaatgtcccttttatt")
    M13_seq_old(5101:5150) = Mani_To_Upper("actggtcgtgtgactggtgaatctgccaatgtaaataatccatttcagac")
    M13_seq_old(5151:5200) = Mani_To_Upper("gattgagcgtcaaaatgtaggtatttccatgagcgtttttcctgttgcaa")
    M13_seq_old(5201:5250) = Mani_To_Upper("tggctggcggtaatattgttctggatattaccagcaaggccgatagtttg")
    M13_seq_old(5251:5300) = Mani_To_Upper("agttcttctactcaggcaagtgatgttattactaatcaaagaagtattgc")
    M13_seq_old(5301:5350) = Mani_To_Upper("tacaacggttaatttgcgtgatggacagactcttttactcggtggcctca")
    M13_seq_old(5351:5400) = Mani_To_Upper("ctgattataaaaacacttctcaagattctggcgtaccgttcctgtctaaa")
    M13_seq_old(5401:5450) = Mani_To_Upper("atccctttaatcggcctcctgtttagctcccgctctgattccaacgagga")
    M13_seq_old(5451:5500) = Mani_To_Upper("aagcacgttatacgtgctcgtcaaagcaaccatagtacgcgccctgtagc")
    M13_seq_old(5501:5550) = Mani_To_Upper("ggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctac")
    M13_seq_old(5551:5600) = Mani_To_Upper("acttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttc")
    M13_seq_old(5601:5650) = Mani_To_Upper("tcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccct")
    M13_seq_old(5651:5700) = Mani_To_Upper("ttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttga")
    M13_seq_old(5701:5750) = Mani_To_Upper("tttgggtgatggttcacgtagtgggccatcgccctgatagacggtttttc")
    M13_seq_old(5751:5800) = Mani_To_Upper("gccctttgacgttggagtccacgttctttaatagtggactcttgttccaa")
    M13_seq_old(5801:5850) = Mani_To_Upper("actggaacaacactcaaccctatctcgggctattcttttgatttataagg")
    M13_seq_old(5851:5900) = Mani_To_Upper("gattttgccgatttcggaaccaccatcaaacaggattttcgcctgctggg")
    M13_seq_old(5901:5950) = Mani_To_Upper("gcaaaccagcgtggaccgcttgctgcaactctctcagggccaggcggtga")
    M13_seq_old(5951:6000) = Mani_To_Upper("agggcaatcagctgttgcccgtctcgctggtgaaaagaaaaaccaccctg")
    M13_seq_old(6001:6050) = Mani_To_Upper("gcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaat")
    M13_seq_old(6051:6100) = Mani_To_Upper("gcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaac")
    M13_seq_old(6101:6150) = Mani_To_Upper("gcaattaatgtgagttagctcactcattaggcaccccaggctttacactt")
    M13_seq_old(6151:6200) = Mani_To_Upper("tatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttc")
    M13_seq_old(6201:6250) = Mani_To_Upper("acacaggaaacagctatgaccatgattacgaattcgagctcggtacccgg")
    M13_seq_old(6251:6300) = Mani_To_Upper("ggatcctctagagtcgacctgcaggcatgcaagcttggcactggccgtcg")
    M13_seq_old(6301:6350) = Mani_To_Upper("ttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgc")
    M13_seq_old(6351:6400) = Mani_To_Upper("cttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccg")
    M13_seq_old(6401:6450) = Mani_To_Upper("caccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgct")
    M13_seq_old(6451:6500) = Mani_To_Upper("ttgcctggtttccggcaccagaagcggtgccggaaagctggctggagtgc")
    M13_seq_old(6501:6550) = Mani_To_Upper("gatcttcctgaggccgatacggtcgtcgtcccctcaaactggcagatgca")
    M13_seq_old(6551:6600) = Mani_To_Upper("cggttacgatgcgcccatctacaccaacgtaacctatcccattacggtca")
    M13_seq_old(6601:6650) = Mani_To_Upper("atccgccgtttgttcccacggagaatccgacgggttgttactcgctcaca")
    M13_seq_old(6651:6700) = Mani_To_Upper("tttaatgttgatgaaagctggctacaggaaggccagacgcgaattatttt")
    M13_seq_old(6701:6750) = Mani_To_Upper("tgatggcgttcctattggttaaaaaatgagctgatttaacaaaaatttaa")
    M13_seq_old(6751:6800) = Mani_To_Upper("cgcgaattttaacaaaatattaacgtttacaatttaaatatttgcttata")
    M13_seq_old(6801:6850) = Mani_To_Upper("caatcttcctgtttttggggcttttctgattatcaaccggggtacatatg")
    M13_seq_old(6851:6900) = Mani_To_Upper("attgacatgctagttttacgattaccgttcatcgattctcttgtttgctc")
    M13_seq_old(6901:6950) = Mani_To_Upper("cagactctcaggcaatgacctgatagcctttgtagatctctcaaaaatag")
    M13_seq_old(6951:7000) = Mani_To_Upper("ctaccctctccggcattaatttatcagctagaacggttgaatatcatatt")
    M13_seq_old(7001:7050) = Mani_To_Upper("gatggtgatttgactgtctccggcctttctcacccttttgaatctttacc")
    M13_seq_old(7051:7100) = Mani_To_Upper("tacacattactcaggcattgcatttaaaatatatgagggttctaaaaatt")
    M13_seq_old(7101:7150) = Mani_To_Upper("tttatccttgcgttgaaataaaggcttctcccgcaaaagtattacagggt")
    M13_seq_old(7151:7200) = Mani_To_Upper("cataatgtttttggtacaaccgatttagctttatgctctgaggctttatt")
    M13_seq_old(7201:7249) = Mani_To_Upper("gcttaattttgctaattctttgccttgcctgtatgatttattggatgtt" )

    ! Compare M13 with old version
    !count = 0
    !do i = 1, len_M13
    !    if(M13_seq(i:i) /= M13_seq_old(i:i)) then
    !        count = count + 1
    !        write(0, "(a)"), trim(adjustl(Int2Str(count)))//" th, pos : "&
    !            //trim(adjustl(Int2Str(i)))//", current : "&
    !            //M13_seq(i:i)//", previous : "//M13_seq_old(i:i)
    !    end if
    !end do
    !stop

    ! Check M13 sequence data
    do i = 1, len_M13
        if(M13_seq(i:i) == "N") then
            write(0, "(a$)"), "Error - Wrong assigned sequence in M13mp18 : "
            write(0, "(a )"), "SeqDesign_Get_M13mp18"
            stop
        end if
    end do
end function SeqDesign_Get_M13mp18

! ---------------------------------------------------------------------------------------

! Import sequence from txt file
! Last updated on Saturday 9 Apr 2016 by Hyungmin
subroutine SeqDesign_Import_Sequence(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j, count, len_seq, n_base_scaf, base, across
    integer, parameter :: max_len = 20000
    character(max_len) :: seq

    ! Encoding should be UTF-8
    open(unit=701, file="seq.txt")

    ! Read sequence and set the sequence lenghth
    read(701, "(a)"), seq
    read(701, "(a)"), seq
    len_seq = len(trim(seq))

    ! Change upper case if low case
    seq = Mani_To_Upper(seq)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "6.6. Set sequence from file"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of sequence from file         : "
        write(i, "(i7)"), len_seq
        call Space(i, 11)
        write(i, "(a$)"), "* The number of scaffold strands           : "
        write(i, "(i7)"), dna.n_scaf
        call Space(i, 11)
        write(i, "(a$)"), "* The number of bases in scaffold strands  : "
        write(i, "(i7)"), dna.n_base_scaf
        write(i, "(a )")
    end do

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(0, "(a$)"), "Error - The number of bases in scaffolds are not consistent : "
        write(0, "(a )"), "SeqDesign_Set_M13mp18"
        stop
    end if

    ! Set scaffold sequence from file
    do i = 1, dna.n_strand

        ! Set sequence for scaffold
        if(dna.strand(i).types /= "scaf") cycle

        ! Find the starting point (down == -1)
        base   = Mani_Go_Start_Base(dna, i)
        across = dna.top(base).across

        do j = 1, dna.strand(i).n_base

            count = j+para_set_start_scaf-1
            if(count >= len_seq) then
                count = mod(count, len_seq)
                if(count == 0) count = len_seq
            end if

            ! Assign sequence from file data
            dna.top(base).seq = seq(count:count)

            ! Set complementary sequence
            if(across /= -1) then
                dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
            end if

            ! Update base
            if(j /= dna.strand(i).n_base) then
                base   = dna.top(base).up
                across = dna.top(base).across
            end if
        end do
    end do
end subroutine SeqDesign_Import_Sequence

! ---------------------------------------------------------------------------------------

! Set random sequence
! Last updated on Saturday 9 Apr 2016 by Hyungmin
subroutine SeqDesign_Set_Rand_Sequence(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j, len_seq, n_base_scaf, base, across

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "6.6. Set random sequence"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of scaffold strands           : "
        write(i, "(i7)"), dna.n_scaf
        call Space(i, 11)
        write(i, "(a$)"), "* The number of bases in scaffold strands  : "
        write(i, "(i7)"), dna.n_base_scaf
        write(i, "(a )")
    end do

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(0, "(a$)"), "Error - The number of bases in scaffolds are not consistent : "
        write(0, "(a )"), "SeqDesign_Set_M13mp18"
        stop
    end if

    ! Set scaffold sequence from file
    do i = 1, dna.n_strand

        ! For scaffold strand
        if(dna.strand(i).types /= "scaf") cycle

        do j = 1, dna.strand(i).n_base

            base   = dna.strand(i).base(j)
            across = dna.top(base).across

            ! Assign random sequence
            dna.top(base).seq = SeqDesign_Get_Rand_Sequence()

            ! Set complementary sequence
            dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
        end do
    end do
end subroutine SeqDesign_Set_Rand_Sequence

! ---------------------------------------------------------------------------------------

! Write txt file for strand data
! Last updated on Thuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Write_Strand(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer, allocatable :: length_stap(:)

    double precision :: length
    integer :: start_iniL, end_iniL, start_sec, end_sec
    integer :: i, j, base, strt_base, end_base, n_base, count
    character(4)   :: types
    character(100) :: seq
    character(200) :: path

    if(para_write_701 == .false.) return

    path = trim(prob.path_work1)
    open(unit=701, file=trim(path)//"strand.txt",   form="formatted")

    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                      1. Sequence data                          +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")

    ! DNA strand data
    write(701, "(a, i7)"), " The total number of strands : ", dna.n_strand
    write(701, "(a)")

    dna.len_min_stap =  10000
    dna.len_max_stap = -10000

    ! Write sequence data
    do i = 1, dna.n_strand
        write(701, "(i10, a$)"), i, " th strand "

        if(dna.strand(i).types == "scaf") then
            write(701, "(a$)"), "[[scaf]]"
        else if(dna.strand(i).types == "stap") then
            write(701, "(a$)"), "( stap )"

            if(dna.strand(i).n_base < dna.len_min_stap) dna.len_min_stap = dna.strand(i).n_base
            if(dna.strand(i).n_base > dna.len_max_stap) dna.len_max_stap = dna.strand(i).n_base

            if(0 .and. dna.strand(i).n_base > 70) then

                do j = 0, 701, 701
                    call space(j, 5)
                    write(j, "(a     )"), "=================================================="
                    call space(j, 5)
                    write(j, "(a     )"), "There are long staple strand"
                    call space(j, 5)
                    write(j, "(a, i4$)"), "Strand #", i
                    call space(j, 5)
                    write(j, "(a, i4 )"), ", # of nucleotides : ", dna.strand(i).n_base
                    call space(j, 5)
                    write(j, "(a, i4 )"), "The current parameter, para-gap_xover_nick : ", para_gap_xover_nick
                    call space(j, 5)
                    write(j, "(a, i4 )"), "The recommended parameter value is ", para_gap_xover_nick - 1
                    call space(j, 5)
                    write(j, "(a     )"), "=================================================="
                    stop
                end do
            end if
        end if

        write(701, "(a, i4$)"), ", # of 14nts :", dna.strand(i).n_14nt
        write(701, "(a, i4$)"), ", start edge # :", mesh.node(dna.top(dna.strand(i).base(1)).node).iniL
        write(701, "(a$    )"), ", # of bases :"
        write(701, "(i6, a$)"), dna.strand(i).n_base, " -> "

        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            write(701, "(a$)"), dna.top(base).seq
        end do
        write(701, "(a)")
    end do
    write(701, "(a)")

    ! Write csv file for squence data
    if(para_write_711 == .true.) then
        open(unit=702, file=trim(path)//"sequence.csv", form="formatted")

        ! DNA strand data
        write(702, "(a)"), " The total number of strands : "//trim(adjustl(Int2Str(dna.n_strand)))
        write(702, "(a)"), ","
        write(702, "(a)"), " Type, num, # nucleo, sequence"

        ! Write sequence data
        do i = 1, dna.n_strand
            if(dna.strand(i).types == "scaf") then
                write(702, "(a$)"), "scaf,"//trim(adjustl(Int2Str(i)))//","
            else if(dna.strand(i).types == "stap") then
                write(702, "(a$)"), "stap,"//trim(adjustl(Int2Str(i)))//","
            end if
            write(702, "(a$)"), trim(adjustl(Int2Str(dna.strand(i).n_base)))//","

            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                write(702, "(a$)"), dna.top(base).seq
            end do
            write(702, "(a)")
        end do
        close(unit=702)
    end if

    ! Write graphical output
    call SeqDesign_Write_Graphical_Output(prob, geom, mesh, dna, 701)

    ! Write unpaired nucleotide and poly Tn loop information
    write(701, "(a)"), " ======================================================== "
    write(701, "(a)")
    write(701, "(a)"), "    [Unpaired nucleotides and poly Tn loop information]   "
    write(701, "(a)")
    write(701, "(a)"), " ======================================================== "
    write(701, "(a)")

    dna.n_nt_unpaired_scaf = 0
    dna.n_nt_unpaired_stap = 0
    dna.n_unpaired_scaf    = 0
    dna.n_unpaired_stap    = 0

    n_base = 0
    do i = 1, dna.n_top

        ! Find the end bases in basepairs
        if(dna.top(i).up /= -1 .and. dna.top(i).dn /= -1) then
            !
            !                i
            ! *--*--*--*--*->*--*--*--*-->
            ! |  |  |  |  |  |
            ! *<-*--*--*--*--*
            if( dna.top(i).across /= -1 .and. &
                dna.top(dna.top(i).up).node == -1 .and. &
                dna.top(dna.top(i).dn).node /= -1 ) then

                base       = dna.top(i).id
                types      = dna.strand(dna.top(i).strand).types
                start_iniL = mesh.node(dna.top(dna.top(base).dn).node).iniL
                start_sec  = mesh.node(dna.top(dna.top(base).dn).node).sec
                count      = 0
                n_base     = n_base + 1
                strt_base  = base

                if(types == "scaf") dna.n_unpaired_scaf = dna.n_unpaired_scaf + 1
                if(types == "stap") dna.n_unpaired_stap = dna.n_unpaired_stap + 1

                ! Count the number of bases
                do
                    base = dna.top(base).up
                    if(dna.top(base).across /= -1) exit
                    count = count + 1
                    seq(count:count) = dna.top(base).seq

                    if(types == "scaf") dna.n_nt_unpaired_scaf = dna.n_nt_unpaired_scaf + 1
                    if(types == "stap") dna.n_nt_unpaired_stap = dna.n_nt_unpaired_stap + 1
                end do

                end_iniL = mesh.node(dna.top(base).node).iniL
                end_sec  = mesh.node(dna.top(base).node).sec
                end_base = dna.top(base).id
                length   = Size_Vector(dna.top(strt_base).pos - dna.top(end_base).pos)

                write(701, "(i10, a$        )"), n_base, "th "//trim(types)
                write(701, "(a, i3, a, f5.2$)"), ", # of bases : ", count, ", Total length : ", length
                !write(701, "(a, f5.2, a$    )"), ", Divided length : ", length/dble(count+1), ", Edge(Section)) : "
                write(701, "(a, f5.2, a$    )"), ", Divided length : ", length/dble(para_dist_pp) - 1.0d0, ", Edge(Section)) : "
                write(701, "(i3, a, i3, a$  )"), start_iniL, "(", start_sec, ") -> "
                write(701, "(i3, a, i3, a$  )"), end_iniL,   "(", end_sec,   "), Sequence : "

                do j = 1, count
                    write(701, "(a$)"), seq(j:j)
                end do
                write(701, "(a)")
            end if
        end if
    end do
    write(701, "(a)")

    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                     3. Strand information                      +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")

    ! Strand information
    do i = 1, dna.n_strand
        write(701, "(i10, a$)"), i, " th strand -> type : "//trim(dna.strand(i).types)
        write(701, "(a, i6$ )"), ", # of bases : ", dna.strand(i).n_base
        write(701, "(a, l, a)"), ", circular : ",   dna.strand(i).b_circular, ", bases : "

        ! Print bases numbering
        write(701, "(i20, a$)"), dna.strand(i).base(1), " -> "
        do j = 2, dna.strand(i).n_base - 1
            if(mod(j, 100)== 0) then
                write(701, "(a9)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//" -> "
            else if(mod(j, 100)== 1) then
                call space(701, 15)
                write(701, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//" -> "
            else
                write(701, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//" -> "
            end if
        end do
        write(701, "(i7)"), dna.strand(i).base(dna.strand(i).n_base)
        write(701, "(a )")
    end do

    ! DNA base information
    write(701, "(a )")
    write(701, "(a )"), " DNA base information in detail"
    write(701, "(a$)"), " The total number of bases : "
    write(701, "(i7)"), dna.n_top

    do i = 1, dna.n_top
        write(701, "(i10, a$)"), dna.top(i).id, " th base ->"
        write(701, "(a$     )"), " type : "//trim(dna.strand(dna.top(i).strand).types)
        write(701, "(a$     )"), ", seq : "//trim(dna.top(i).seq)
        write(701, "(a, i6$ )"), ", up : ",       dna.top(i).up
        write(701, "(a, i6$ )"), ", down : ",     dna.top(i).dn
        write(701, "(a, i6$ )"), ", across : ",   dna.top(i).across
        write(701, "(a, i6$ )"), ", basepair : ", dna.top(i).node
        write(701, "(a, i4$ )"), ", strand : ",   dna.top(i).strand
        write(701, "(a, i5  )"), ", address : ",  dna.top(i).address
    end do

    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                       4. Staple length                         +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")

    ! Allocate and initialize memory
    allocate(length_stap(dna.len_min_stap:dna.len_max_stap))
    length_stap(dna.len_min_stap:dna.len_max_stap) = 0

    ! Find staple length
    do i = 1, dna.n_strand
        if(dna.strand(i).types == "stap") then
            length_stap(dna.strand(i).n_base) = length_stap(dna.strand(i).n_base) + 1
        end if
    end do

    ! Write staple length
    do i = dna.len_min_stap, dna.len_max_stap
        if(length_stap(i) /= 0) then
            write(701, "(2i6)"), i, length_stap(i)
        end if
    end do

    ! Deallocate memory
    deallocate(length_stap)

    close(unit=701)
end subroutine SeqDesign_Write_Strand

! ---------------------------------------------------------------------------------------

! Write graphical output
! Last updated on Tuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Write_Graphical_Output(prob, geom, mesh, dna, unit)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: unit

    ! Section type data
    type :: SecType
        integer :: start_bp, end_bp
        integer :: start_nei_edge, start_nei_sec
        integer :: end_nei_sec, end_nei_edge

        integer,   allocatable :: bp(:)
        integer,   allocatable :: node(:)
        integer,   allocatable :: strnd_scaf(:), strnd_stap(:)
        integer,   allocatable :: xover_scaf(:), xover_stap(:)
        integer,   allocatable :: nick_scaf(:),  nick_stap(:)
        character, allocatable :: seq_scaf(:),   seq_stap(:)
    end type SecType

    ! Initial line data
    type :: EdgeType
        integer :: n_sec
        type(SecType), allocatable :: sec(:)
    end type EdgeType

    ! Tecplot drawing data type
    type :: TecType
        integer :: types
        double precision :: size
        double precision :: color
        double precision :: pos_cen(3)
        double precision :: pos_node1(3)
        double precision :: pos_node2(3)
    end type TecType

    type(EdgeType), allocatable :: edge(:)
    type(TecType),  allocatable :: tec(:)
    integer,        allocatable :: conn_scaf(:,:), conn_stap(:,:)

    integer :: i, j, k, m, iniL, sec, strnd, factor, pos, type_tec(6), tec_factor
    integer :: node, dn_node, up_node, bp, max_bp, min_bp, mid_bp, mid_sec
    integer :: n_conn_scaf, n_conn_stap, n_edge, n_tec
    logical :: b_conn_scaf, b_conn_stap, b_sec
    character(200) :: path
    character(20)  :: col_list(20)

    col_list(1:4)   = ["magenta",         "orange red",   "deep pink",     "turquoise"      ]
    col_list(5:8)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(9:12)  = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(13:16) = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(17:20) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Find maximum and minimum base pair ID
    max_bp = mesh.node(1).bp
    min_bp = mesh.node(1).bp
    do i = 1, mesh.n_node
        if(mesh.node(i).bp > max_bp) max_bp = mesh.node(i).bp
        if(mesh.node(i).bp < min_bp) min_bp = mesh.node(i).bp
    end do
    min_bp = min_bp - 2
    max_bp = max_bp + 2

    ! Allocate and initialize edge data
    n_edge = geom.n_iniL
    allocate(edge(n_edge))
    do i = 1, n_edge
        edge(i).n_sec = geom.n_sec
        allocate(edge(i).sec(edge(i).n_sec))

        do j = 1, edge(i).n_sec
            allocate(edge(i).sec(j).bp        (min_bp:max_bp))
            allocate(edge(i).sec(j).node      (min_bp:max_bp))
            allocate(edge(i).sec(j).strnd_scaf(min_bp:max_bp))
            allocate(edge(i).sec(j).strnd_stap(min_bp:max_bp))
            allocate(edge(i).sec(j).xover_scaf(min_bp:max_bp))
            allocate(edge(i).sec(j).xover_stap(min_bp:max_bp))
            allocate(edge(i).sec(j).nick_scaf (min_bp:max_bp))
            allocate(edge(i).sec(j).nick_stap (min_bp:max_bp))
            allocate(edge(i).sec(j).seq_scaf  (min_bp:max_bp))
            allocate(edge(i).sec(j).seq_stap  (min_bp:max_bp))

            do k = min_bp, max_bp
                edge(i).sec(j).bp(k)         = k
                edge(i).sec(j).node(k)       = -1
                edge(i).sec(j).strnd_scaf(k) = -1
                edge(i).sec(j).strnd_stap(k) = -1
                edge(i).sec(j).xover_scaf(k) = -1
                edge(i).sec(j).xover_stap(k) = -1
                edge(i).sec(j).nick_scaf(k)  = -1
                edge(i).sec(j).nick_stap(k)  = -1
                edge(i).sec(j).seq_scaf(k)   = "N"
                edge(i).sec(j).seq_stap(k)   = "N"
            end do
        end do
    end do

    ! Allocate conn data
    allocate(conn_scaf(dna.n_xover_scaf, 3))
    allocate(conn_stap(dna.n_xover_stap, 3))

    ! Build node information based on initial edges with cross-section
    do i = 1, mesh.n_node

        ! Find starting node depending on cross-section ID
        node = mesh.node(i).id
        sec  = mesh.node(i).sec
        iniL = mesh.node(node).iniL

        ! It depends on cross-section ID and direction
        if(mesh.node(i).dn == -1 .and. mod(mesh.node(i).sec, 2) == 0) then

            ! If section ID is even and negative z-direction
            ! Find neiboring edge and section
            if(dna.top(node).dn /= -1) then

                ! To avoid unpaired nucleotide
                dn_node = dna.top(node).dn
                do
                    if(dna.top(dn_node).across /= -1) exit
                    dn_node = dna.top(dn_node).dn
                end do

                edge(iniL).sec(sec+1).start_nei_edge = mesh.node(dn_node).iniL
                edge(iniL).sec(sec+1).start_nei_sec  = mesh.node(dn_node).sec
            else
                edge(iniL).sec(sec+1).start_nei_edge = -1
                edge(iniL).sec(sec+1).start_nei_sec  = -1
            end if

            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).up
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp

            ! Find neiboring edge and section
            if(dna.top(edge(iniL).sec(sec+1).node(bp)).up /= -1) then

                ! To avoid unpaired nucleotide
                up_node = dna.top(edge(iniL).sec(sec+1).node(bp)).up
                do
                    if(dna.top(up_node).across /= -1) exit
                    up_node = dna.top(up_node).up
                end do

                edge(iniL).sec(sec+1).end_nei_edge = mesh.node(up_node).iniL
                edge(iniL).sec(sec+1).end_nei_sec  = mesh.node(up_node).sec
            else
                edge(iniL).sec(sec+1).end_nei_edge = -1
                edge(iniL).sec(sec+1).end_nei_sec  = -1
            end if
        else if(mesh.node(i).up == -1 .and. mod(mesh.node(i).sec, 2) == 1) then

            ! If section ID is odd and positive z-direction
            ! Find neiboring edge and section
            if(dna.top(node).up /= -1) then
                
                ! To avoid unpaired nucleotide
                up_node = dna.top(node).up
                do
                    if(dna.top(up_node).across /= -1) exit
                    up_node = dna.top(up_node).up
                end do
                
                edge(iniL).sec(sec+1).start_nei_edge = mesh.node(up_node).iniL
                edge(iniL).sec(sec+1).start_nei_sec  = mesh.node(up_node).sec
            else
                edge(iniL).sec(sec+1).start_nei_edge = -1
                edge(iniL).sec(sec+1).start_nei_sec  = -1
            end if

            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).dn
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp

            ! Find neiboring edge and section
            if(dna.top(edge(iniL).sec(sec+1).node(bp)).dn /= -1) then

                ! To avoid unpaired nucleotide
                dn_node = dna.top(edge(iniL).sec(sec+1).node(bp)).dn
                do
                    if(dna.top(dn_node).across /= -1) exit
                    dn_node = dna.top(dn_node).dn
                end do

                edge(iniL).sec(sec+1).end_nei_edge = mesh.node(dn_node).iniL
                edge(iniL).sec(sec+1).end_nei_sec  = mesh.node(dn_node).sec
            else
                edge(iniL).sec(sec+1).end_nei_edge = -1
                edge(iniL).sec(sec+1).end_nei_sec  = -1
            end if
        end if
    end do

    ! Build additional data fields
    do i = 1, dna.n_top
        node = dna.top(i).node

        ! Skip the loop if poly Tn loop and unpaired nucleotides
        if(node == -1) cycle

        strnd = dna.top(i).strand
        sec   = mesh.node(node).sec
        iniL  = mesh.node(node).iniL
        bp    = mesh.node(node).bp

        if(dna.strand(strnd).types == "scaf") then

            ! For scaffold strand
            edge(iniL).sec(sec+1).seq_scaf(bp) = dna.top(i).seq
            if(dna.top(i).xover /= -1) then
                edge(iniL).sec(sec+1).xover_scaf(bp) = &
                    mesh.node(dna.top(dna.top(i).xover).node).sec
            end if
            if(dna.top(i).up == -1) then
                if(mod(sec, 2) == 0) then
                    edge(iniL).sec(sec+1).nick_scaf(bp) = 2
                else
                    edge(iniL).sec(sec+1).nick_scaf(bp) = 3
                end if
            end if
            if(dna.top(i).dn == -1) edge(iniL).sec(sec+1).nick_scaf(bp) = 1
            edge(iniL).sec(sec+1).strnd_scaf(bp) = strnd
        else if(dna.strand(strnd).types == "stap") then

            ! for staple strand
            edge(iniL).sec(sec+1).seq_stap(bp) = dna.top(i).seq
            if(dna.top(i).xover /= -1) then
                edge(iniL).sec(sec+1).xover_stap(bp) = &
                    mesh.node(dna.top(dna.top(i).xover).node).sec
            end if
            if(dna.top(i).up == -1) then
                if(mod(sec, 2) == 0) then
                    edge(iniL).sec(sec+1).nick_stap(bp) = 3
                else
                    edge(iniL).sec(sec+1).nick_stap(bp) = 2
                end if
            end if
            if(dna.top(i).dn == -1) edge(iniL).sec(sec+1).nick_stap(bp) = 1
            edge(iniL).sec(sec+1).strnd_stap(bp) = strnd
        end if
    end do

    ! Open file stream
    if(para_write_710 == .true.) then
        path = trim(prob.path_work1)//trim(prob.name_file)
        do i = 1, n_edge
            open(unit = 710+i, file=trim(path)//"_design_edge"&
                //trim(adjustl(Int2Str(i)))//".bild", form="formatted")
        end do

        ! Tecplot drawing
        if(para_output_Tecplot == "on") then

            ! 1 - normal base, 2 - xover base, 3, xover, 4 - nick, 5 - right arrow, 6- left arrow
            allocate(tec(2*dna.n_top))

            do i = 1, 2*dna.n_top
                tec(i).types          = 0
                tec(i).pos_cen(1:3)   = 0.0d0
                tec(i).pos_node1(1:3) = 0.0d0
                tec(i).pos_node2(1:3) = 0.0d0
                tec(i).size           = 0.0d0
                tec(i).color          = 0.0d0
            end do
            type_tec(1:6) = 0
            n_tec         = 0
        end if
    end if

    ! Print information based on edges
    do i = 1, n_edge

        ! Print base pair ID
        if(i == 1) then
            !do k = min_bp, max_bp
            !    write(unit, "(i5$)"), edge(i).sec(j).bp(k)
            !end do
            !write(unit, "(a)"); write(unit, "(a)")

            !do k = min_bp, max_bp
            !    write(unit, "(a5$)"), "-----"
            !end do
            !write(unit, "(a)"); write(unit, "(a)")

            write(unit, "(a)")
            write(unit, "(a)"), "--------------------------------------------------------------------"
            write(unit, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
            write(unit, "(a)"), "|+                                                                +|"
            write(unit, "(a)"), "|+                      2. Graphical output                       +|"
            write(unit, "(a)"), "|+                                                                +|"
            write(unit, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
            write(unit, "(a)"), "--------------------------------------------------------------------"
            write(unit, "(a)")
            write(unit, "(a)"), "1. [-], [num], [.], [>], [<] : Nucleotide"
            write(unit, "(a)"), "2. [|]                       : Crossover"
            write(unit, "(a)"), "3. [A], [T], [G], [C]        : Sequence"
            write(unit, "(a)"), "4. [ a: b], c                : From starting base ID(a) to ending base ID(b), total # of bases (c)"
            write(unit, "(a)"), "5. [5'], [3']                : Strand direction"
            write(unit, "(a)"); write(unit, "(a)")
        end if

        ! --------------------------------------------------
        !
        ! Scaffold strand
        !
        ! --------------------------------------------------
        n_conn_scaf = 0
        do j = 1, edge(i).n_sec

            ! Factor : 0,0,1,1,2,2,3,3,...
            mid_sec    = int((6*edge(i).n_sec-5+1)/2)
            factor     = 3*(int(j+1)/2-1) + mid_sec
            tec_factor = (i-1)*(5+6*edge(i).n_sec-5-3*(int(edge(i).n_sec+1)/2-1))

            ! Print edge and point information, Edge 1 - (point 1 -> point 2)
            if(j == 1) then
                write(unit, "(a)") " ==================================================="
                call Space(unit, 10)
                write(unit, "(a)"), " [Edge "//trim(adjustl(Int2Str(i)))//&
                    " : point "//trim(adjustl(Int2Str(geom.iniL(i).poi(1))))//&
                    " -> point "//trim(adjustl(Int2Str(geom.iniL(i).poi(2))))//"]"
                write(unit, "(a)") " ==================================================="
            end if

            if(j == 1) then
                write(unit, "(a)")
                write(unit, "(a)"), "       [[ SCAFFOLD STRAND ]]"
                write(unit, "(a)")
            end if

            ! Print section ID and base length, sec xx -->> [xxx:xxx], xxx : - 26 length
            write(unit, "(a5$)"), " sec "
            write(unit, "(a2$)"), trim(adjustl(Int2Str(j-1)))
            write(unit, "(a4$)"), " - ["
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).start_bp + para_start_bp_ID - 1)))
            write(unit, "(a1$)"), ":"
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp + para_start_bp_ID - 1)))
            write(unit, "(a3$)"), "], "
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp - edge(i).sec(j).start_bp + 1)))
            write(unit, "(a2$)"), " :"
            call Space(unit, 4)

            ! Graphic representation, bases with crossover and nick
            do k = min_bp, max_bp

                mid_bp = int((max_bp + min_bp) / 2)

                if(edge(i).sec(j).node(k) == -1) then
                    if(edge(i).sec(j).start_bp - 2 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "5"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "3"
                    else if(edge(i).sec(j).end_bp + 1 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "3"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "5"
                    else if(edge(i).sec(j).start_bp - 1 == k .or. edge(i).sec(j).end_bp + 2 == k) then
                        write(unit, "(a$)"), "'"
                    else
                        write(unit, "(a$)"), " "
                    end if
                else

                    if(edge(i).sec(j).xover_scaf(k) /= -1) then
                        !write(unit, "(a$)"), "+"

                        if(para_output_design == "arrow") then
                            if(edge(i).sec(j).xover_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).xover_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).xover_scaf(k))))
                            end if
                        end if

                        if(para_output_design == "seq") write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)

                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        n_conn_scaf = n_conn_scaf + 1
                        conn_scaf(n_conn_scaf, 1) = edge(i).sec(j).xover_scaf(k)
                        conn_scaf(n_conn_scaf, 2) = k
                        conn_scaf(n_conn_scaf, 3) = edge(i).sec(j).strnd_scaf(k)

                        ! Draw sphere (crossover point)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_output_Tecplot == "on") then
                                n_tec                 = n_tec + 1
                                type_tec(3)           = type_tec(3) + 1
                                tec(n_tec).types      = 3
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 1) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "."
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw sphere (nick point)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_output_Tecplot == "on") then
                                n_tec                 = n_tec + 1
                                type_tec(4)           = type_tec(4) + 1
                                tec(n_tec).types      = 4
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 2) then

                        if(para_output_design == "arrow")  write(unit, "(a$)"), ">"
                        if(para_output_design == "seq")    write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw arrow (right arrow, >)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_output_Tecplot == "on") then
                                n_tec                   = n_tec + 1
                                type_tec(5)             = type_tec(5) + 1
                                tec(n_tec).types        = 5
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp) - 0.05d0
                                tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).pos_node1(1) = 0.9d0
                                tec(n_tec).color        = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 3) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "<"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw arrow (left arrow, <)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_output_Tecplot == "on") then
                                n_tec                   = n_tec + 1
                                type_tec(6)             = type_tec(6) + 1
                                tec(n_tec).types        = 6
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp) + 0.05d0
                                tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).pos_node1(1) =-0.9d0
                                tec(n_tec).color        = dble(1)
                            end if
                        end if
                    else

                        if(para_output_design == "arrow") write(unit, "(a$)"), "-"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw cylinder (normal base)
                        if(para_write_710 == .true.) then
                            strnd = edge(i).sec(j).strnd_scaf(k)

                            if(mesh.node(edge(i).sec(j).node(k)).up == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 0) then
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                   = n_tec + 1
                                    type_tec(5)             = type_tec(5) + 1
                                    tec(n_tec).types        = 5
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                    tec(n_tec).pos_node1(1) = 1.0d0
                                    tec(n_tec).color        = dble(1)
                                end if
                            else if(mesh.node(edge(i).sec(j).node(k)).up == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 1) then
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                   = n_tec + 1
                                    type_tec(6)             = type_tec(6) + 1
                                    tec(n_tec).types        = 6
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                    tec(n_tec).pos_node1(1) =-1.0d0
                                    tec(n_tec).color        = dble(1)
                                end if
                            else
                                write(710+i, "(a     )"), ".color steel blue"
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(1)          = type_tec(1) + 1
                                    tec(n_tec).types     = 1
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor) - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor) - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                    end if
                end if
            end do

            ! For Neighbor connection status, [ START : 3( 3), END : 5( 0)]
            call Space(unit, 5)
            write(unit, "(a, i2$)"), "[ START : ", edge(i).sec(j).start_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).start_nei_sec
            write(unit, "(a, i2$)"), "), END : ",  edge(i).sec(j).end_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).end_nei_sec
            write(unit, "(a$    )"), ") ]"
            write(unit, "(a     )")

            ! Draw crossovers
            call Space(unit, 30)
            do k = min_bp, max_bp
                if(edge(i).sec(j).node(k) == -1) then
                    write(unit, "(a$)"), " "
                else
                    b_conn_scaf = .false.
                    do m = 1, n_conn_scaf
                        if(conn_scaf(m, 1) >= j .and. conn_scaf(m, 2) == k) then
                            b_conn_scaf = .true.
                            exit
                        end if
                    end do

                    if(b_conn_scaf == .true.) then
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(para_write_710 == .true.) then
                            if(mod(j, 2) == 1) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            else
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                    else
                        if( edge(i).sec(j).start_nei_edge == i .and. &
                            edge(i).sec(j).start_bp       == k .and. &
                            edge(i).sec(j).start_nei_sec  == j ) then

                            ! Neighbor connection for starting bp
                            write(unit, "(a$)"), "|"

                            ! Draw cylinder
                            if(para_write_710 == .true.) then
                                if(mod(j, 2) == 1) then
                                    write(710+i, "(a$    )"), ".cylinder "
                                    write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                    write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                    write(710+i, "(1f9.2 )"), 0.1d0

                                    if(para_output_Tecplot == "on") then
                                        n_tec                = n_tec + 1
                                        type_tec(2)          = type_tec(2) + 1
                                        tec(n_tec).types     = 2
                                        tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                        tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                        tec(n_tec).color     = dble(1)
                                    end if
                                else
                                    write(710+i, "(a$    )"), ".cylinder "
                                    write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                    write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                    write(710+i, "(1f9.2 )"), 0.1d0

                                    if(para_output_Tecplot == "on") then
                                        n_tec                = n_tec + 1
                                        type_tec(2)          = type_tec(2) + 1
                                        tec(n_tec).types     = 2
                                        tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                        tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                        tec(n_tec).color     = dble(1)
                                    end if
                                end if
                            end if
                        else if( edge(i).sec(j).end_nei_edge == i .and. &
                            edge(i).sec(j).end_bp       == k .and. &
                            edge(i).sec(j).end_nei_sec  == j ) then

                        ! Neighbor connection for ending bp
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(para_write_710 == .true.) then
                            if(mod(j, 2) == 1) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            else
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                        else
                            write(unit, "(a$)"), " "
                        end if
                    end if
                end if
            end do
            write(unit, "(a)")
        end do

        ! ==================================================
        !
        ! Staple strand
        !
        ! ==================================================
        n_conn_stap = 0
        b_sec       = .true.
        do j = 1, edge(i).n_sec

            ! Factor : 0,0,1,1,2,2,3,3,...
            mid_sec    = int((6*edge(i).n_sec-5+1)/2)
            factor     = 3*(int(j + 1)/2 - 1) + mid_sec
            tec_factor = (i-1)*(5+6*edge(i).n_sec-5-3*(int(edge(i).n_sec+1)/2-1))

            if(j == 1) then
                write(unit, "(a)")
                write(unit, "(a)"), "       [[ STAPLE STRAND ]]"
                write(unit, "(a)")
            end if

            ! Print section ID and base length, sec xx -->> [xxx:xxx], xxx : - 26 length
            write(unit, "(a5$)"), " sec "
            write(unit, "(a2$)"), trim(adjustl(Int2Str(j-1)))
            write(unit, "(a4$)"), " - ["
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).start_bp + para_start_bp_ID - 1)))
            write(unit, "(a1$)"), ":"
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp + para_start_bp_ID - 1)))
            write(unit, "(a3$)"), "], "
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp - edge(i).sec(j).start_bp + 1)))
            write(unit, "(a2$)"), " :"
            call Space(unit, 4)

            ! Graphic representation, bases with crossover and nick
            do k = min_bp, max_bp
                mid_bp = int((min_bp+max_bp)/2)
                if(edge(i).sec(j).node(k) == -1) then
                    if(edge(i).sec(j).start_bp - 2 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "3"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "5"
                    else if(edge(i).sec(j).end_bp + 1 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "5"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "3"
                    else if(edge(i).sec(j).start_bp - 1 == k .or. edge(i).sec(j).end_bp + 2 == k) then
                        write(unit, "(a$)"), "'"
                    else
                        write(unit, "(a$)"), " "
                    end if
                else

                    if(mod(j, 2) == 1) then
                        pos = 6*int((j+1)/2)-5 + 1
                    else
                        pos = 6*int((j+1)/2)-5 + 2
                    end if
                    pos = 2*pos-1-factor

                    if(edge(i).sec(j).xover_stap(k) /= -1) then
                        !write(unit, "(a$)"), "+"

                        if(para_output_design == "arrow") then
                            if(edge(i).sec(j).xover_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).xover_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).xover_stap(k))))
                            end if
                        end if

                        if(para_output_design == "seq") write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)

                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        n_conn_stap = n_conn_stap + 1
                        conn_stap(n_conn_stap, 1) = edge(i).sec(j).xover_stap(k)
                        conn_stap(n_conn_stap, 2) = k
                        conn_stap(n_conn_stap, 3) = edge(i).sec(j).strnd_stap(k)

                        ! Draw sphere (crossover point)
                        if(para_write_710 == .true.) then
                            strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(pos), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_output_Tecplot == "on") then
                                n_tec                 = n_tec + 1
                                type_tec(3)           = type_tec(3) + 1
                                tec(n_tec).types      = 3
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(pos) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 1) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "."
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw sphere (nick point)
                        if(para_write_710 == .true.) then
                            strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(pos), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_output_Tecplot == "on") then
                                n_tec                 = n_tec + 1
                                type_tec(4)           = type_tec(4) + 1
                                tec(n_tec).types      = 4
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(pos) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 2) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), ">"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw arrow (right arrow, >)
                        if(para_write_710 == .true.) then
                            strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_output_Tecplot == "on") then
                                n_tec                   = n_tec + 1
                                type_tec(5)             = type_tec(5) + 1
                                tec(n_tec).types        = 5
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                tec(n_tec).pos_node1(1) = 0.8d0
                                tec(n_tec).color        = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 3) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "<"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw arrow (left arrow, <)
                        if(para_write_710 == .true.) then
                            strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_output_Tecplot == "on") then
                                n_tec                   = n_tec + 1
                                type_tec(6)             = type_tec(6) + 1
                                tec(n_tec).types        = 6
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                tec(n_tec).pos_node1(1) =-0.8d0
                                tec(n_tec).color        = dble(strnd)
                            end if
                        end if
                    else

                        if(para_output_design == "arrow") write(unit, "(a$)"), "-"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Write Chimera
                        if(para_write_710 == .true.) then
                            strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                            if(mesh.node(edge(i).sec(j).node(k)).dn == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 1) then
                                ! Draw arrow (right arrow, >)
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.09d0, 0.4d0, 0.55d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                   = n_tec + 1
                                    type_tec(5)             = type_tec(5) + 1
                                    tec(n_tec).types        = 5
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                    tec(n_tec).pos_node1(1) = 1.0d0
                                    tec(n_tec).color        = dble(strnd)
                                end if
                            else if(mesh.node(edge(i).sec(j).node(k)).dn == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 0) then
                                ! Draw arrow (left arrow, <)
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.09d0, 0.4d0, 0.55d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                   = n_tec + 1
                                    type_tec(6)             = type_tec(6) + 1
                                    tec(n_tec).types        = 6
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                    tec(n_tec).pos_node1(1) =-1.0d0
                                    tec(n_tec).color        = dble(strnd)
                                end if
                            else
                                ! Draw cylinder, standard base
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(1)          = type_tec(1) + 1
                                    tec(n_tec).types     = 1
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp)-0.4d0, -dble(pos) - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp)+0.4d0, -dble(pos) - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            end if
                        end if
                    end if

                    if(para_write_710 == .true.) then
                        ! Draw sequence
                        if(mod(j, 2) == 0) then
                            write(710+i, "(a, 3f9.3)"), ".cmov ", dble(k-mid_bp-0.3d0), -dble(pos)-1.0d0, 0.0d0
                        else
                            write(710+i, "(a, 3f9.3)"), ".cmov ", dble(k-mid_bp-0.3d0), -dble(pos)+0.4d0, 0.0d0
                        end if
                        write(710+i, "(a)"), ".font arial 20 bold"
                        write(710+i, "(a)"), dna.top(dna.top(edge(i).sec(j).node(k)).across).seq

                        ! Draw section ID
                        if(b_sec == .true.) then
                            write(710+i, "(a)"), ".color red"
                            if(mod(j, 2) == 0) then
                                write(710+i, "(a, 3f9.3)"), ".cmov ", dble(min_bp-mid_bp-5.0d0), -dble(pos)-1.0d0, 0.0d0
                            else
                                write(710+i, "(a, 3f9.3)"), ".cmov ", dble(min_bp-mid_bp-5.0d0), -dble(pos)+0.4d0, 0.0d0
                            end if
                            write(710+i, "(a)"), ".font arial 20 bold"
                            write(710+i, "(a)"), "sec "//trim(adjustl(Int2Str(j-1)))
                            b_sec = .false.
                        end if
                    end if
                end if
            end do
            b_sec = .true.

            ! For Neighbor connection status, [ START : 3( 3), END : 5( 0)]
            call Space(unit, 5)
            do k = 1, edge(i).n_sec
                if(edge(i).sec(k).start_nei_edge == i) then
                    edge(i).sec(k).start_nei_edge = -1
                    edge(i).sec(k).start_nei_sec  = -1
                end if
                if(edge(i).sec(k).end_nei_edge == i) then
                    edge(i).sec(k).end_nei_edge = -1
                    edge(i).sec(k).end_nei_sec  = -1
                end if
            end do
            write(unit, "(a, i2$)"), "[ START : ", edge(i).sec(j).start_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).start_nei_sec
            write(unit, "(a, i2$)"), "), END : ",  edge(i).sec(j).end_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).end_nei_sec
            write(unit, "(a$    )"), ") ]"
            write(unit, "(a     )")

            ! Draw crossovers
            call Space(unit, 30)
            do k = min_bp, max_bp
                if(edge(i).sec(j).node(k) == -1) then
                    write(unit, "(a$)"), " "
                else
                    b_conn_stap = .false.
                    do m = 1, n_conn_stap
                        if(conn_stap(m, 1) >= j .and. conn_stap(m, 2) == k) then
                            b_conn_stap = .true.
                            exit
                        end if
                    end do

                    if(b_conn_stap == .true.) then
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(para_write_710 == .true.) then
                            strnd = mod(conn_stap(m, 3), 20) + 1
                            if(mod(j, 2) == 0) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5+5-factor)+0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-5-factor)-0.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-5+5-factor)+0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-5-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            else
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2+1-factor)+0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-1-factor)-0.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_output_Tecplot == "on") then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2+1-factor)+0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-1-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            end if
                        end if
                    else
                        write(unit, "(a$)"), " "
                    end if
                end if
            end do
            write(unit, "(a)")
        end do
        write(unit, "(a)")
    end do

    ! Write tecplot output
    if(para_output_Tecplot == "on") then

        ! Open file stream
        path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
        open(unit=709, file=trim(path)//"_design_edge.dat", form="formatted")
        do i = 1, 6

            if(type_tec(i) == 0) cycle

            write(709, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
            write(709, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3", "c"'
            write(709, "(a$)"), 'ZONE F = FEPOINT'

            if(i == 1 .or. i == 2) then

                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(2*type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal connectivity
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_node1(1:3), 1.0d0, 1.0d0, 1.0d0, tec(j).color
                        write(709, "(7f8.2)"), tec(j).pos_node2(1:3), 1.0d0, 1.0d0, 1.0d0, tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j*2-1, j*2
                end do
            else if(type_tec(i) > 0 .and. (i == 3 .or. i == 4)) then

                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal position
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_cen(1:3), tec(j).size, 1.0d0, 1.0d0, tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j, j
                end do
            else if(type_tec(i) > 0 .and. (i == 5 .or. i == 6)) then
                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal position
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_cen(1:3), tec(j).pos_node1(1:3), tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j, j
                end do
            end if
        end do
    end if

    ! Deallocate memory
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            deallocate(edge(i).sec(j).bp        )
            deallocate(edge(i).sec(j).node      )
            deallocate(edge(i).sec(j).strnd_scaf)
            deallocate(edge(i).sec(j).strnd_stap)
            deallocate(edge(i).sec(j).xover_scaf)
            deallocate(edge(i).sec(j).xover_stap)
            deallocate(edge(i).sec(j).nick_scaf )
            deallocate(edge(i).sec(j).nick_stap )
            deallocate(edge(i).sec(j).seq_scaf  )
            deallocate(edge(i).sec(j).seq_stap  )
        end do
        deallocate(edge(i).sec)
    end do
    deallocate(edge)
    deallocate(conn_scaf)
    deallocate(conn_stap)

    do i = 1, n_edge
        close(unit = 710+i)
    end do

    ! Deallocate memory and close file for Tecplot output
    if(para_output_Tecplot == "on") then
        deallocate(tec)
        close(unit=709)
    end if
end subroutine SeqDesign_Write_Graphical_Output

! ---------------------------------------------------------------------------------------

! Write atom model by dnaTop and strand data
! Last updated on Wednesday 25 May 2016 by Hyungmin
subroutine SeqDesign_Chimera_Atom(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j, base, up, xover, across
    logical :: f_axis
    character(200) :: path

    if(para_write_702 == .false.) return

    f_axis = para_chimera_axis

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=702, file=trim(path)//"_atom_nick.bild", form="formatted")

    ! For all bases
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base

            if(dna.strand(i).types == "scaf") write(702, "(a)"), ".color steel blue"
            if(dna.strand(i).types == "stap") write(702, "(a)"), ".color orange"

            ! Draw bases
            base = dna.strand(i).base(j)
            write(702, "(a$    )"), ".sphere "
            write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(702, "(1f9.3 )"), 0.15d0

            ! Draw backbones
            up = dna.top(base).up
            if(up /= -1) then
                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
                write(702, "(3f9.3$)"), dna.top(up).pos(1:3)
                write(702, "(1f9.3 )"), 0.05d0
            end if

            ! Draw crossovers
            xover = dna.top(base).xover
            if(xover /= -1 .and. base < xover) then

                if(dna.strand(i).types == "scaf") write(702, "(a)"), ".color blue"
                if(dna.strand(i).types == "stap") write(702, "(a)"), ".color red"

                pos_1(:) = dna.top(base).pos(1:3)
                pos_2(:) = dna.top(xover).pos(1:3)

                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), pos_1(1:3)
                write(702, "(3f9.3$)"), pos_2(1:3)
                write(702, "(1f9.3 )"), 0.08d0
            end if

            ! Draw the Watson-Crick connections
            across = dna.top(base).across
            if(across > 0) then
                write(702, "(a     )"), ".color light gray"
                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
                write(702, "(3f9.3$)"), dna.top(across).pos(1:3)
                write(702, "(1f9.3 )"), 0.025d0
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(702, "(a)"), ".translate 0.0 0.0 0.0"
        write(702, "(a)"), ".scale 0.5"
        write(702, "(a)"), ".color grey"
        write(702, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(702, "(a)"), ".color red"             ! x-axis
        write(702, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(702, "(a)"), ".color blue"            ! y-axis
        write(702, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(702, "(a)"), ".color yellow"          ! z-axis
        write(702, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=702)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=702, file=trim(path)//"_atom_nick.dat", form="formatted")

    write(702, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'

    ! For bases in scaffold strands
    do i = 1, dna.n_strand

        write(702, "(a )"), 'VARIABLES = "X", "Y", "Z", "Weight"'
        write(702, "(a$)"), 'ZONE F = FEPOINT'
        write(702, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.strand(i).n_base)))
        write(702, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.strand(i).n_base - 1)))
        write(702, "(a )"), ', ET=LINESEG'

        ! Draw bases
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            write(702, "(4f9.3)"), dna.top(base).pos(1:3), 1.0d0
        end do

        ! Write elements
        do j = 1, dna.strand(i).n_base - 1
            write(702, "(2i7)"), j, j + 1
        end do
    end do

    write(702, "(a )"), 'VARIABLES = "X", "Y", "Z", "Weight"'
    write(702, "(a$)"), 'ZONE F = FEPOINT'
    write(702, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.n_base_scaf*2)))
    write(702, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.n_base_scaf)))
    write(702, "(a )"), ', ET=LINESEG'

    ! For bases in scaffold strands
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base
            if(dna.strand(i).types == "scaf") then
                base   = dna.strand(i).base(j)
                across = dna.top(base).across

                if(across > 0) then
                    write(702, "(4f9.3)"), dna.top(base).pos(1:3),   1.0d0
                    write(702, "(4f9.3)"), dna.top(across).pos(1:3), 1.0d0
                end if
            end if
        end do
    end do

    ! Write elements
    do i = 1, dna.n_base_scaf
        write(702, "(2i7)"), 2*i-1, 2*i
    end do

    close(unit=702)
end subroutine SeqDesign_Chimera_Atom

! ---------------------------------------------------------------------------------------

! Write scaffold or staple route for Chimera
! Last updated on Tuesday 9 August 2016 by Hyungmin
subroutine SeqDesign_Chimera_Route(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    double precision, allocatable :: base_scaf(:,:), base_stap(:,:)

    double precision :: pos_1(3), pos_2(3), vec(3), radius
    integer :: i, j, k, cur_base, up_base, down_base, cur_node, up_node
    integer :: n_base_scaf, n_base_stap
    logical :: f_axis
    character(200) :: path

    if(para_write_703 == .false.) return

    ! Exception for Tecplot drawing
    if(para_output_Tecplot == "on") then
        allocate(base_scaf (dna.n_base_scaf *2, 3))
        allocate(base_stap (dna.n_base_stap *2, 3))
        n_base_scaf = 0
        n_base_stap = 0
    end if

    ! Set option
    f_axis = para_chimera_axis

    ! File open for route step
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=703, file=trim(path)//"_route6_scaf.bild", form="formatted")
    open(unit=704, file=trim(path)//"_route6_stap.bild", form="formatted")

    ! --------------------------------------------------
    !
    ! For bases in the scaffold strand
    !
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! For scaffold strand
        if(dna.strand(i).types /= "scaf") cycle

        ! Loop for strand
        do j = 1, dna.strand(i).n_base

            ! Find the base
            cur_base  = dna.strand(i).base(j)
            up_base   = dna.top(cur_base).up
            down_base = dna.top(cur_base).dn

            ! Find the node for cur_base and up_base
            cur_node = dna.top(cur_base).node

            if(cur_node == -1) cycle

            if(up_base /= -1) then

                ! Find the node for up_base
                up_node = dna.top(up_base).node
            else

                ! Draw point and skip if there is no upward base
                !write(703, "(a     )"), ".color dark green"
                !write(703, "(a$    )"), ".sphere "
                !write(703, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                !write(703, "(1f9.3 )"), 0.2d0
                cycle
            end if

            ! Draw starting point as the point and arrow
            if(down_base == -1 .and. cur_node /= -1 .and. up_node /= -1) then
                pos_1(1:3) = mesh.node(cur_node).pos(1:3)
                pos_2(1:3) = mesh.node(up_node).pos(1:3)
                vec(1:3)   = Normalize_Vector(pos_2 - pos_1)

                !write(703, "(a     )"), ".color dark green"
                !write(703, "(a$    )"), ".sphere "
                !write(703, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                !write(703, "(1f9.3 )"), 0.2d0

                !write(703, "(a$    )"), ".arrow "
                !write(703, "(3f8.2$)"), pos_1(1:3)
                !write(703, "(3f8.2$)"), pos_1(1:3) + vec(1:3)*1.3d0
                !write(703, "(2f8.2 )"), 0.21d0, 0.4d0
            end if

            ! Color section for Tn loop and crossover
            if(cur_node == -1 .or. up_node == -1) then

                ! Check Tn loop
                if(cur_node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base = dna.top(up_base).up
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
                write(703, "(a)"), ".color orange"
                radius = 0.1d0
            else if(dna.top(cur_base).xover == up_base) then

                ! Crossovers color
                write(703, "(a)"), ".color red"
                radius = 0.1d0
            else

                ! Scaffold strand color
                write(703, "(a)"), ".color steel blue"
                radius = 0.1d0
            end if

            ! Draw route path
            write(703, "(a$    )"), ".cylinder "
            write(703, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
            write(703, "(3f9.3$)"), mesh.node(up_node).pos(1:3)
            write(703, "(1f9.3 )"), radius

            if(para_output_Tecplot == "on") then
                base_scaf(n_base_scaf + 1, 1:3) = mesh.node(cur_node).pos(1:3)
                base_scaf(n_base_scaf + 2, 1:3) = mesh.node(up_node).pos(1:3)
                n_base_scaf = n_base_scaf + 2
            end if
        end do
    end do

    ! --------------------------------------------------
    ! For bases in the staple strand
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! For staple strand
        if(dna.strand(i).types /= "stap") cycle

        ! Loop for strand
        do j = 1, dna.strand(i).n_base

            ! Find the base
            cur_base  = dna.strand(i).base(j)
            up_base   = dna.top(cur_base).up
            down_base = dna.top(cur_base).dn

            ! Find the node for cur_base ann up_base
            cur_node = dna.top(cur_base).node

            if(cur_node == -1) cycle

            if(up_base /= -1) then

                ! Find the node for up_base
                up_node = dna.top(up_base).node
            else

                ! Draw point and skip if there is no upward base
                !write(704, "(a     )"), ".color dark green"
                !write(704, "(a$    )"), ".sphere "
                !write(704, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                !write(704, "(1f9.3 )"), 0.2d0
                cycle
            end if

            ! Draw starting point as the point and arrow
            if(down_base == -1 .and. cur_node /= -1 .and. up_node /= -1) then
                pos_1(1:3) = mesh.node(cur_node).pos(1:3)
                pos_2(1:3) = mesh.node(up_node).pos(1:3)
                vec(1:3)   = Normalize_Vector(pos_2 - pos_1)

                !write(704, "(a     )"), ".color dark green"
                !write(704, "(a$    )"), ".sphere "
                !write(704, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                !write(704, "(1f9.3 )"), 0.2d0

                !write(704, "(a$    )"), ".arrow "
                !write(704, "(3f8.2$)"), pos_1(1:3)
                !write(704, "(3f8.2$)"), pos_1(1:3) + vec(1:3)*1.5d0
                !write(704, "(2f8.2 )"), 0.31d0, 0.55d0
            end if

            ! Color section for Tn loop and crossover
            if(cur_node == -1 .or. up_node == -1) then

                ! Check Tn loop
                if(cur_node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base = dna.top(up_base).up
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
                write(704, "(a)"), ".color orange"
                radius = 0.1d0
            else if(dna.top(cur_base).xover == up_base) then

                ! Crossovers color
                write(704, "(a)"), ".color red"
                radius = 0.1d0
            else

                ! Staple strand color
                write(704, "(a)"), ".color orange"
                radius = 0.1d0
            end if

            ! Draw route path
            write(704, "(a$    )"), ".cylinder "
            write(704, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
            write(704, "(3f9.3$)"), mesh.node(up_node).pos(1:3)
            write(704, "(1f9.3 )"), radius

            if(para_output_Tecplot == "on") then
                base_stap(n_base_stap + 1, 1:3) = mesh.node(cur_node).pos(1:3)
                base_stap(n_base_stap + 2, 1:3) = mesh.node(up_node).pos(1:3)
                n_base_stap = n_base_stap + 2
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) then
        do i = 0, 1
            write(703+i, "(a)"), ".translate 0.0 0.0 0.0"
            write(703+i, "(a)"), ".scale 0.5"
            write(703+i, "(a)"), ".color grey"
            write(703+i, "(a)"), ".sphere 0 0 0 0.5"      ! Center
            write(703+i, "(a)"), ".color red"             ! x-axis
            write(703+i, "(a)"), ".arrow 0 0 0 4 0 0 "
            write(703+i, "(a)"), ".color blue"            ! y-axis
            write(703+i, "(a)"), ".arrow 0 0 0 0 4 0 "
            write(703+i, "(a)"), ".color yellow"          ! z-axis
            write(703+i, "(a)"), ".arrow 0 0 0 0 0 4 "
        end do
    end if

    close(unit=703)
    close(unit=704)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=703, file=trim(path)//"_route6_scaf.dat", form="formatted")
    open(unit=704, file=trim(path)//"_route6_stap.dat", form="formatted")

    ! For scaffold bases
    write(703, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(703, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(703, "(a$)"), 'ZONE F = FEPOINT'
    write(703, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_scaf)))
    write(703, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_scaf/2)))
    write(703, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_scaf
        write(703, "(3f9.3$)"), base_scaf(i, 1:3)
        write(703, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_base_scaf/2
        write(703, "(2i7)"), 2*i, 2*i-1
    end do

    ! For crossover for scaffold strand
    write(704, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(704, "(a$)"), 'ZONE F = FEPOINT'
    write(704, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_stap)))
    write(704, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_stap/2)))
    write(704, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_stap
        write(704, "(3f9.3$)"), base_stap(i, 1:3)
        write(704, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_base_stap/2
        write(704, "(2i7)"), 2*i, 2*i-1
    end do

    if(allocated(base_scaf))  deallocate(base_scaf)
    if(allocated(base_stap))  deallocate(base_stap)

    close(unit=703)
    close(unit=704)
end subroutine SeqDesign_Chimera_Route

! ---------------------------------------------------------------------------------------

! Chimera sequence design
! Last updated on Thursday 23 June 2016 by Hyungmin
subroutine SeqDesign_Chimera_Sequence_Design(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(in)    :: dna

    double precision, allocatable :: base_scaf(:,:), base_stap(:,:)
    integer,          allocatable :: col_stap(:)

    double precision :: vec_n(3), vec_jn(3), vec_jt1(3), vec_jt2(3)
    double precision :: pos_1(3), pos_2(3), vec(3), RGB(3)
    integer :: i, j, k, base, up_base, node, up_node
    integer :: croL, iniL, count_scaf, count_stap
    logical :: f_axis
    character(15)  :: col_list(16)
    character(200) :: path

    if(para_write_705 == .false.) return

    col_list(1:4)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(5:8)   = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(9:12)  = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(13:16) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Set option
    f_axis = para_chimera_axis

    ! File open for sequence design
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=705, file=trim(path)//"_sequence_design.bild", form="formatted")

    if(para_output_Tecplot == "on") then
        allocate(base_scaf(dna.n_base_scaf*2, 3))
        allocate(base_stap(dna.n_base_stap*2, 4))
        allocate(col_stap(dna.n_base_stap*2))
        count_scaf = 0
        count_stap = 0
    end if

    ! --------------------------------------------------
    !
    ! Nucleotides of staples
    !
    ! --------------------------------------------------
    do i = 1, geom.n_iniL
        geom.iniL(i).n_xover = 0
    end do

    do i = 1, dna.n_base_scaf
        if(dna.top(i).xover /= -1) then
            geom.iniL(mesh.node(dna.top(i).node).iniL).n_xover &
                = geom.iniL(mesh.node(dna.top(i).node).iniL).n_xover + 1
        end if
    end do

    ! --------------------------------------------------
    !
    ! For bases of the scaffold strand
    !
    ! --------------------------------------------------
    write(705, "(a)"), ".color steel blue"
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).types /= "scaf") cycle

        do j = 1, dna.strand(i).n_base

            ! Find the current and up base
            base    = dna.strand(i).base(j)
            up_base = dna.top(base).up

            ! Skip if there is no upward base
            if(up_base == -1) cycle

            ! Find the node
            node    = dna.top(base).node
            up_node = dna.top(up_base).node

            ! Draw bases in unpaired nucleotides
            if(node == -1 .or. up_node == -1) then
                if(node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base = dna.top(up_base).up
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
            end if

            ! Get cross-sectional line
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line
            !iniL  = mesh.node(node).iniL
            !pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            !pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            !vec(1:3) = mesh.node(node).pos(1:3) - pos_2(1:3)
            !vec(1:3) = Normalize_Vector(vec(1:3))

            ! Vector projection
            !vec_n(1:3)   = Normalize_Vector(pos_2 - pos_1)
            !vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            !vec_jt1(1:3) = Normalize_Vector(vec - vec_jn)

            ! Get cross-sectional line #1
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(up_node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line #2
            !iniL  = mesh.node(up_node).iniL
            !pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            !pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            !vec(1:3) = mesh.node(up_node).pos(1:3) - pos_2(1:3)
            !vec(1:3) = Normalize_Vector(vec(1:3))

            ! Vector projection
            !vec_n(1:3)   = Normalize_Vector(pos_2 - pos_1)
            !vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            !vec_jt2(1:3) = Normalize_Vector(vec - vec_jn)

            ! Draw route path
            write(705, "(a$    )"), ".cylinder "
            write(705, "(3f9.3$)"), mesh.node(node).pos    !+ vec_jt1*0.5d0
            write(705, "(3f9.3$)"), mesh.node(up_node).pos !+ vec_jt2*0.5d0
            write(705, "(1f9.3 )"), 0.1d0

            ! Draw starting point
            if(dna.top(base).dn == -1) then
                write(705, "(a     )"), ".color red"
                write(705, "(a$    )"), ".sphere "
                write(705, "(3f9.3$)"), mesh.node(node).pos !+ vec_jt1*0.5d0
                write(705, "(1f9.3 )"), 0.2d0
                write(705, "(a     )"), ".color steel blue"
            end if

            if(para_output_Tecplot == "on") then
                base_scaf(count_scaf + 1, 1:3) = mesh.node(node).pos    !+ vec_jt1*0.5d0
                base_scaf(count_scaf + 2, 1:3) = mesh.node(up_node).pos !+ vec_jt2*0.5d0
                count_scaf = count_scaf + 2
            end if
        end do
    end do

    ! --------------------------------------------------
    !
    ! For bases of the staple strand
    !
    ! --------------------------------------------------
    vec_jt1(1:3) = 0.0d0
    vec_jt2(1:3) = 0.0d0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).types /= "stap") cycle

        !RGB(1:3) = [imod(irand(), 256), imod(irand(), 256), imod(irand(), 256)]
        !write(705, "(a, 3f9.4)"), ".color ", dble(RGB(1:3))/255.0d0
        write(705, "(a)"), ".color "//trim(col_list(mod(i-1, 16) + 1))

        !if(i /= 2) cycle
        do j = 1, dna.strand(i).n_base

            ! Find the base
            base    = dna.strand(i).base(j)
            up_base = dna.top(base).up

            ! Ending point as arrow
            if(up_base == -1) then
                pos_1 = mesh.node(dna.top(dna.top(base).dn).node).pos
                pos_2 = mesh.node(dna.top(base).node).pos
                vec   = pos_1 - pos_2
                write(705, "(a$    )"), ".arrow "
                write(705, "(3f9.3$)"), pos_1 + 1.5d0*vec - vec_jt1*0.5d0
                write(705, "(3f9.3$)"), pos_2 - 0.9d0*vec - vec_jt2*0.5d0
                write(705, "(3f9.3 )"), 0.1d0, 0.3d0, 0.3d0
                cycle
            end if

            ! Find node
            node    = dna.top(base).node
            up_node = dna.top(up_base).node

            ! Draw bases in Tn loop
            if(node == -1 .or. up_node == -1) then
                if(node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base = dna.top(up_base).up
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
            end if

            ! Get cross-sectional line
            pos_1(1:3) = 0.0d0
            pos_2(1:3) = 0.0d0
            do k = 1, geom.n_sec
                croL = mesh.node(node).croL
                pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
                pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            end do

            ! Get cross-sectional line
            iniL  = mesh.node(node).iniL
            pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            vec(1:3) = mesh.node(node).pos(1:3) - pos_2(1:3)
            vec(1:3) = Normalize_Vector(vec(1:3))

            ! Vector projection
            vec_n(1:3)   = Normalize_Vector(pos_2 - pos_1)
            vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            vec_jt1(1:3) = Normalize_Vector(vec - vec_jn)

            ! Get cross-sectional line #1
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(up_node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line #2
            iniL  = mesh.node(up_node).iniL
            pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            vec(1:3) = mesh.node(up_node).pos(1:3) - pos_2(1:3)
            vec(1:3) = Normalize_Vector(vec(1:3))

            ! Vector projection
            vec_n(1:3)   = Normalize_Vector(pos_2 - pos_1)
            vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            vec_jt2(1:3) = Normalize_Vector(vec - vec_jn)

            ! Draw route path
            write(705, "(a$    )"), ".cylinder "
            write(705, "(3f9.3$)"), mesh.node(node).pos(1:3)    - vec_jt1*0.5d0
            write(705, "(3f9.3$)"), mesh.node(up_node).pos(1:3) - vec_jt2*0.5d0
            write(705, "(1f9.3 )"), 0.1d0

            if(para_output_Tecplot == "on") then
                base_stap(count_stap + 1, 1:3) = mesh.node(node).pos(1:3)    - vec_jt1*0.5d0
                base_stap(count_stap + 2, 1:3) = mesh.node(up_node).pos(1:3) - vec_jt2*0.5d0

                if(mesh.node(node).iniL == mesh.node(up_node).iniL) then
                    col_stap(count_stap + 1) = geom.iniL(mesh.node(node).iniL).n_xover
                    col_stap(count_stap + 2) = geom.iniL(mesh.node(node).iniL).n_xover
                else
                    col_stap(count_stap + 1) = 20
                    col_stap(count_stap + 2) = 20
                end if

                count_stap = count_stap + 2
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(705, "(a)"), ".translate 0.0 0.0 0.0"
        write(705, "(a)"), ".scale 0.5"
        write(705, "(a)"), ".color grey"
        write(705, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(705, "(a)"), ".color red"             ! x-axis
        write(705, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(705, "(a)"), ".color blue"            ! y-axis
        write(705, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(705, "(a)"), ".color yellow"          ! z-axis
        write(705, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if

    close(unit=705)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=705, file=trim(path)//"_sequence_design.dat", form="formatted")
    open(unit=706, file=trim(path)//"_init_geo_col.dat",    form="formatted")

    ! For scaffold bases
    write(705, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(705, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(705, "(a$)"), 'ZONE F = FEPOINT'
    write(705, "(a$)"), ', N='//trim(adjustl(Int2Str(count_scaf)))
    write(705, "(a$)"), ', E='//trim(adjustl(Int2Str(count_scaf/2)))
    write(705, "(a )"), ', ET=LINESEG'

    do i = 1, count_scaf
        write(705, "(3f9.3$)"), base_scaf(i, 1:3)
        write(705, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, count_scaf/2
        write(705, "(2i7)"), 2*i, 2*i-1
    end do

    ! For staple bases
    write(705, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(705, "(a$)"), 'ZONE F = FEPOINT'
    write(705, "(a$)"), ', N='//trim(adjustl(Int2Str(count_stap)))
    write(705, "(a$)"), ', E='//trim(adjustl(Int2Str(count_stap/2)))
    write(705, "(a )"), ', ET=LINESEG'

    do i = 1, count_stap
        write(705, "(3f9.3$)"), base_stap(i, 1:3)
        write(705, "(1f9.3 )"), dble(col_stap(i))
    end do
    do i = 1, count_stap/2
        write(705, "(2i7)"), 2*i, 2*i-1
    end do

    ! ---------------------------------------------
    !
    ! New initial geometry
    !
    ! ---------------------------------------------
    write(706, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(706, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(706, "(a$)"), 'ZONE F = FEPOINT'
    write(706, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_modP)))
    write(706, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(706, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_iniL
        write(706, "(3f9.3$)"), geom.modP(geom.iniL(i).poi(1)).pos(1:3)
        write(706, "(1i7 )"), geom.iniL(i).n_xover
        write(706, "(3f9.3$)"), geom.modP(geom.iniL(i).poi(2)).pos(1:3)
        write(706, "(1i7 )"), geom.iniL(i).n_xover
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(706, "(1i7$)"), geom.iniL(i).poi(1)
        write(706, "(1i7 )"), geom.iniL(i).poi(2)
    end do

    ! Deallocate memory
    if(allocated(base_scaf)) deallocate(base_scaf)
    if(allocated(base_stap)) deallocate(base_stap)
    if(allocated(col_stap))  deallocate(col_stap)

    close(unit=705)
    close(unit=706)
end subroutine SeqDesign_Chimera_Sequence_Design

! ---------------------------------------------------------------------------------------

! Write Chimera file for strand and sequence
! Last updated on Thuesday 1 Mar 2016 by Hyungmin
subroutine SeqDesign_Chimera_Strand(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    double precision :: vec(3), azure(3), tweetybird(3), green(3), carmine(3)
    integer :: i, j, k, base, down, up, across
    logical :: f_axis
    character(15)  :: col_list(16)
    character(200) :: path

    if(para_write_706 == .false.) return

    f_axis = para_chimera_axis

    col_list(1:4)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(5:8)   = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(9:12)  = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(13:16) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Colors for four bases, reference: http://www.umass.edu/molvis/tutorials/dna/atgc.htm
    azure      = [0.000d0, 127.0d0, 255.0d0] / 255.0d0       ! color for A
    tweetybird = [253.0d0, 245.0d0, 1.000d0] / 255.0d0       ! color for T
    green      = [0.000d0, 255.0d0, 0.000d0] / 255.0d0       ! color for G
    carmine    = [150.0d0, 0.000d0, 24.00d0] / 255.0d0       ! color for C

    ! Write the file
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=706, file=trim(path)//"_strand.bild",   form="formatted")
    open(unit=707, file=trim(path)//"_sequence.bild", form="formatted")

    ! Draw the structure using strand data
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base

            ! Find the base
            base = dna.strand(i).base(j)

            ! Color for strand
            write(706, "(a)"), ".color " // trim(col_list(mod(i-1, 16) + 1))

            ! Color for sequence
            if(dna.top(base).seq == "A") then
                write(707, "(a, 3f9.3)"), ".color ", azure(1:3)
            else if(dna.top(base).seq == "T") then
                write(707, "(a, 3f9.3)"), ".color ", tweetybird(1:3)
            else if(dna.top(base).seq == "G") then
                write(707, "(a, 3f9.3)"), ".color ", green(1:3)
            else if(dna.top(base).seq == "C") then
                write(707, "(a, 3f9.3)"), ".color ", carmine(1:3)
            else
                write(0, "(a$)"), "Error - Not assigned sequence : " 
                write(0, "(a )"), "SeqDesign_Chimera_Strand" 
                stop
            end if

            write(706, "(a$    )"), ".sphere "
            write(706, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(706, "(1f9.3 )"), 0.15d0

            write(707, "(a$    )"), ".sphere "
            write(707, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(707, "(1f9.3 )"), 0.15d0

            do k = 0, 1
                ! Draw the backbones
                down = dna.top(base).dn
                if(down >= 0) then
                    write(706+k, "(a     )"), ".color light gray"
                    write(706+k, "(a$    )"), ".cylinder "
                    write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
                    write(706+k, "(3f9.3$)"), dna.top(down).pos(1:3)
                    write(706+k, "(1f9.3 )"), 0.05d0
                end if

                ! Draw the Watson-Crick connections
                across = dna.top(base).across
                if(across >= 0) then
                    write(706+k, "(a     )"), ".color light gray"
                    write(706+k, "(a$    )"), ".cylinder "
                    write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
                    write(706+k, "(3f9.3$)"), dna.top(across).pos(1:3)
                    write(706+k, "(1f9.3 )"), 0.025d0
                end if
            end do
        end do

        ! Draw the strand directionality
        if(dna.strand(i).n_base > 1) then
            ! Get last base
            base = dna.strand(i).base(dna.strand(i).n_base)

            if(dna.strand(i).b_circular == .true.) then
                down = dna.strand(i).base(1)
                vec  = Normalize_Vector(dna.top(down).pos - dna.top(base).pos)
            else
                up  = dna.strand(i).base(dna.strand(i).n_base - 1)
                vec = Normalize_Vector(dna.top(base).pos - dna.top(up).pos)
            end if
        else
            vec = [1.0d0, 0.0d0, 0.0d0]
        end if

        ! Color for strand
        write(706, "(a)"), ".color " // trim(col_list(mod(i-1, 16) + 1))

        ! Color for sequence
        if(dna.top(base).seq == "A") then
            write(707, "(a, 3f9.3)"), ".color ", azure(1:3)
        else if(dna.top(base).seq == "T") then
            write(707, "(a, 3f9.3)"), ".color ", tweetybird(1:3)
        else if(dna.top(base).seq == "G") then
            write(707, "(a, 3f9.3)"), ".color ", green(1:3)
        else if(dna.top(base).seq == "C") then
            write(707, "(a, 3f9.3)"), ".color ", carmine(1:3)
        else
            write(707, "(a)"), ".color light gray"
        end if

        do k = 0, 1
            write(706+k, "(a$    )"), ".cone "
            write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3) + vec(1:3) * 0.36d0
            write(706+k, "(1f9.3 )"), 0.18d0
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) then
        do i = 0, 1
            write(706+i, "(a)"), ".translate 0.0 0.0 0.0"
            write(706+i, "(a)"), ".scale 0.5"
            write(706+i, "(a)"), ".color grey"
            write(706+i, "(a)"), ".sphere 0 0 0 0.5"      ! Center
            write(706+i, "(a)"), ".color red"             ! x-axis
            write(706+i, "(a)"), ".arrow 0 0 0 4 0 0 "
            write(706+i, "(a)"), ".color blue"            ! y-axis
            write(706+i, "(a)"), ".arrow 0 0 0 0 4 0 "
            write(706+i, "(a)"), ".color yellow"          ! z-axis
            write(706+i, "(a)"), ".arrow 0 0 0 0 0 4 "
        end do
    end if

    close(unit=706)
    close(unit=707)
end subroutine SeqDesign_Chimera_Strand

! ---------------------------------------------------------------------------------------

end module SeqDesign