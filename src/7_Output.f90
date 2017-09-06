!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Output
!
!                                                                    Updated : 2017/09/05
!
! Comments: This module is for outputs.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Output

    use Ifport
    use Data_Mesh
    use Data_DNA
    use Data_Bound

    use Mani
    use Para
    use Math

    use Chimera
    use Tecplot

    implicit none

    public  Output_Generation

    private Output_Write_Outputs
    private Output_Write_Out_Sequences
    private Output_Write_Out_Unpaired
    private Output_Write_Out_Graphics
    private Output_Write_Out_Strand_Base
    private Output_Write_Out_Guide_JSON
    private Output_Write_Out_JSON

    private Output_Check_Output
    private Output_Write_Basepair
    private Output_Write_Base
    private Output_Write_CanDo
    private Output_Write_CanDo_New
    private Output_Write_DNA_Info
    private Output_Write_TecPlot
    private Output_Write_ADINA
    private Output_Make_Route_Step
    private Output_Chimera_Route_Step
    private Output_Figure_Route_Step
    private Output_Write_Sequence_CroL

contains

! ---------------------------------------------------------------------------------------

! Generate output files and run postprocessing tools
subroutine Output_Generation(prob, geom, bound, mesh, dna)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound
    type(MeshType),  intent(in) :: mesh
    type(DNAType),   intent(inout) :: dna

    integer :: i

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |                            7. Outputs                              |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
        call Space(i, 6)
        write(i, "(a)"), "7.1. Write output files"
    end do

    ! Check connectivity in dna data
    !call Output_Check_Output(dna)

    ! Write outputs
    call Output_Write_Outputs(prob, geom, bound, mesh, dna)

    ! Write information related with basepair
    call Output_Write_Basepair(prob, mesh, dna)

    ! Write information related with base
    call Output_Write_Base(prob, dna)

    ! Write cndo file for PDB atom generation and CanDo simulation
    call Output_Write_CanDo(prob, mesh, dna)
    !call Output_Write_CanDo_New(prob, mesh, dna)

    !call Output_Write_DNA_Info(prob, dna)

    ! Write TecPlot input file
    call Output_Write_TecPlot(prob, mesh)

    ! Write ADINA input file
    call Output_Write_ADINA(prob, mesh)

    ! Command for Tecplot and Chimera, and figures for output and route step
    if(para_cmd_Tecplot    == "on") call TecPlot_Command_Orientation(prob)
    if(para_cmd_Chimera    == "on") call Chimera_Command_Output(prob)
    if(para_fig_output     == "on") call Chimera_Figure_Output(prob)
    if(para_fig_route_step == "on") call Output_Make_Route_Step(prob, mesh, dna)

    ! Write sequence based on cross-sectional line
    call Output_Write_Sequence_CroL(prob, mesh, dna)
end subroutine Output_Generation

! ---------------------------------------------------------------------------------------

! Check connectivity in dna data
subroutine Output_Check_Output(dna)
    type(DNAType), intent(in) :: dna

    integer :: i, j, count, cur_base, across, xover
    integer, allocatable :: base(:)

    ! --------------------------------------------------
    ! Check the # of dnatop and bases in all strands
    ! --------------------------------------------------
    count = 0
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base
            count = count + 1
        end do
    end do

    ! Check the number of bases in all strands
    if(count /= dna.n_top) then
        write(0, "(a$)"), "Error - # of bases in dna top and strand are not the same : "
        write(0, "(a )"), "Output_Check_Output"
        stop
    end if

    ! --------------------------------------------------
    ! There should be non-circular strand and check # of bases
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Find starting point
        cur_base = dna.strand(i).base(1)
        do
            if(dna.top(cur_base).dn == -1) exit
            cur_base = dna.top(cur_base).dn
        end do

        count = 1
        do
            if(dna.top(cur_base).up == -1) exit
            cur_base = dna.top(cur_base).up
            count = count + 1
        end do

        ! Check the number of bases in dna strand
        if(count /= dna.strand(i).n_base) then
            write(0, "(a$)"), "Error - # of bases in dna strand are not consistent : "
            write(0, "(a )"), "Output_Check_Output"
            stop
        end if
    end do

    ! --------------------------------------------------
    ! Check across
    ! --------------------------------------------------
    do i = 1, dna.n_top
        if(dna.top(i).across /= -1) then
            across = dna.top(i).across
            if(dna.top(across).across /= dna.top(i).id) then
                write(0, "(a$)"), "Error - The across ID is not consistent : "
                write(0, "(a )"), "Output_Check_Output"
                stop
            end if
        end if
    end do

    ! --------------------------------------------------
    ! Check crossover
    ! --------------------------------------------------
    do i = 1, dna.n_top
        if(dna.top(i).xover /= -1) then
            xover = dna.top(i).xover
            if(dna.top(xover).xover /= dna.top(i).id) then
                write(0, "(a$)"), "Error - The crossover ID is not consistent : "
                write(0, "(a )"), "Output_Check_Output"
                stop
            end if
        end if
    end do

    ! --------------------------------------------------
    ! Check whether strand has consistent direction (non-circular strand)
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Allocate base data to store entity
        allocate(base(dna.strand(i).n_base))

        ! Find starting point
        cur_base = dna.strand(i).base(1)
        do
            if(dna.top(cur_base).dn == -1) exit
            cur_base = dna.top(cur_base).dn
        end do

        ! Go base up and save base number in base data
        base(1) = cur_base
        do j = 2, dna.strand(i).n_base
            cur_base = dna.top(cur_base).up
            base(j)  = cur_base
        end do

        ! Go base down and check whether the same bases
        do j = 1, dna.strand(i).n_base
            if(base(dna.strand(i).n_base - j + 1) == cur_base) then
                cur_base = dna.top(cur_base).dn
            else
                write(0, "(a$)"), "Error - The strand has not consistent direction : "
                write(0, "(a )"), "Output_Check_Output"
                stop
            end if
        end do

        ! Deallocate base memory
        deallocate(base)
    end do
end subroutine Output_Check_Output

! ---------------------------------------------------------------------------------------

! Write output file for sequence and json
subroutine Output_Write_Outputs(prob, geom, bound, mesh, dna)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound
    type(MeshType),  intent(in) :: mesh
    type(DNAType),   intent(inout) :: dna

    integer :: max_unpaired

    if(para_write_701 == .false.) return
    open(unit = 701, file=trim(prob.path_work1)//"TXT_Sequence.txt", form="formatted")

    ! Write staple sequence results
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                      1. Sequence data                          +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Sequences(prob, mesh, dna, 701)

    ! Write graphical routing
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                     2. Graphical routing                       +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    write(701, "(a)"), "1. [-], [num], [.], [>], [<] : Nucleotide"
    write(701, "(a)"), "2. [|]                       : Crossover"
    write(701, "(a)"), "3. [A], [T], [G], [C]        : Sequence"
    write(701, "(a)"), "4. [ a: b], c                : From starting base ID(a) to ending base ID(b), total # of bases (c)"
    write(701, "(a)"), "5. [5'], [3']                : Strand direction"
    write(701, "(a)")
    call Output_Write_Out_Graphics(prob, geom, mesh, dna, 701)

    ! Write information on unpaired nucleotides and poly Tn loop
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+    3. Information on unpaired nucleotides and poly Tn loop     +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    max_unpaired = Output_Write_Out_Unpaired(mesh, dna, 701)

    ! Outputs based on strands and nucleotides
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+             4. Strand and nucleotide based outputs             +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Strand_Base(mesh, dna, 701)

    ! Output of staple length
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                       5. Staple length                         +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Staple_Length(dna, 701)

    ! JSON output
    call Output_Write_Out_Guide_JSON(prob, geom, bound, mesh)
    call Output_Write_Out_JSON(prob, geom, mesh, dna, max_unpaired)

    close(unit=701)
end subroutine Output_Write_Outputs

! ---------------------------------------------------------------------------------------

! Write output of sequence data
subroutine Output_Write_Out_Sequences(prob, mesh, dna, unit)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(inout) :: dna
    integer, intent(in) :: unit

    integer :: i, j, base, n_scaf, n_stap, dozen, add
    character(200) :: path

    n_scaf = 0
    n_stap = 0
    dna.len_min_stap =  10000
    dna.len_max_stap = -10000

    ! Write sequence data
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            n_scaf = n_scaf + 1
            write(unit, "(i8, a$)"), n_scaf, " - [scaf]"
        else if(dna.strand(i).type1 == "stap") then
            n_stap = n_stap + 1
            write(unit, "(i8, a$)"), n_stap, " - (stap)"
        end if

        if(dna.strand(i).n_base < dna.len_min_stap) dna.len_min_stap = dna.strand(i).n_base
        if(dna.strand(i).n_base > dna.len_max_stap) dna.len_max_stap = dna.strand(i).n_base

        if(0 .and. dna.strand(i).n_base > 80) then
            do j = 0, unit, unit
                call space(j,5); write(j, "(a     )"), "=================================================="
                call space(j,5); write(j, "(a     )"), "There are long staple strand"
                call space(j,5); write(j, "(a, i4$)"), "Strand #", i
                call space(j,5); write(j, "(a, i4 )"), ", # of nucleotides : ", dna.strand(i).n_base
                call space(j,5); write(j, "(a, i4 )"), "The current parameter, para-gap_xover_nick : ", para_gap_xover_nick
                call space(j,5); write(j, "(a, i4 )"), "The recommended parameter value is ", para_gap_xover_nick - 1
                call space(j,5); write(j, "(a     )"), "=================================================="
                stop
            end do
        end if

        if(dna.strand(i).type1 == "stap") then
            write(unit, "(a$)"), "-"//dna.strand(i).type2
            write(unit, "(a$)"), ", # of 14nt seeds: "//trim(adjustl(Int2Str(dna.strand(i).n_14nt)))
        end if
        write(unit, "(a$)"), ", # of nts: "//trim(adjustl(Int2Str(dna.strand(i).n_base)))//" -> "

        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            write(unit, "(a$)"), dna.top(base).seq
            base = dna.Top(base).up
        end do
        write(unit, "(a)")
        if(dna.strand(i).type1 == "scaf") write(unit, "(a)")
    end do
    write(unit, "(a)")

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=702, file=trim(path)//"_16_sequence.csv", form = "formatted")

    ! DNA sequence data according to strands
    write(702, "(a)")
    write(702, "(a)"), "No, Type1, Type2, ABCD, Strand name, # of nts, Sequence"
    write(702, "(a)")

    ! Write sequence data
    n_scaf = 0
    n_stap = 0
    add = 0
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            n_scaf = n_scaf + 1
            write(702, "(a$)"), trim(adjustl(Int2Str(n_scaf)))//", "
            write(702, "(a$)"), "Scaffold, -, -, -, "
        else if(dna.strand(i).type1 == "stap") then
            n_stap = n_stap + 1
            write(702, "(a$)"), trim(adjustl(Int2Str(n_stap)))//", "
            write(702, "(a$)"), "Staple, "//dna.strand(i).type2//", "

            dozen = mod(n_stap, 12)
            if(dozen == 0) dozen = 12
            if(dozen == 1) add = add + 1
            write(702, "(a$)"), achar(65 + add - 1)//trim(adjustl(Int2Str(dozen)))//", "

            write(702, "(a$)"), trim(adjustl(prob.name_prob))//"_"
            if(para_vertex_design == "flat")    write(702, "(a$)"), "F_"
            if(para_vertex_design == "beveled") write(702, "(a$)"), "M_"
            write(702, "(a$)"), trim(adjustl(Int2Str(prob.n_bp_edge)))//"bp_"
            write(702, "(a$)"), trim(adjustl(Int2Str(n_stap)))//", "
        end if
        write(702, "(a$)"), trim(adjustl(Int2Str(dna.strand(i).n_base)))//", "

        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            write(702, "(a$)"), dna.top(base).seq
            base = dna.Top(base).up
        end do

        if(dna.strand(i).type1 == "scaf") write(702, "(a)")
        write(702, "(a)")
    end do
    close(unit=702)
end subroutine Output_Write_Out_Sequences

! ---------------------------------------------------------------------------------------

! Write graphical output
subroutine Output_Write_Out_Graphics(prob, geom, mesh, dna, unit)
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

        if(dna.strand(strnd).type1 == "scaf") then

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
        else if(dna.strand(strnd).type1 == "stap") then

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
end subroutine Output_Write_Out_Graphics

! ---------------------------------------------------------------------------------------

! Write output of unpaired nucleotides
function Output_Write_Out_Unpaired(mesh, dna, unit) result(max_base)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: unit

    double precision :: length
    integer :: s_base, e_base, s_iniL, e_iniL, s_sec, e_sec, max_base
    integer :: i, j, n_base, base, count
    character(100) :: seq
    character(4) :: types

    dna.n_nt_unpaired_scaf = 0
    dna.n_nt_unpaired_stap = 0
    dna.n_unpaired_scaf    = 0
    dna.n_unpaired_stap    = 0

    max_base = 0
    n_base   = 0
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

                base   = dna.top(i).id
                types  = dna.strand(dna.top(i).strand).type1
                s_iniL = mesh.node(dna.top(dna.top(base).dn).node).iniL
                s_sec  = mesh.node(dna.top(dna.top(base).dn).node).sec
                count  = 0
                n_base = n_base + 1
                s_base = base

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

                e_iniL = mesh.node(dna.top(base).node).iniL
                e_sec  = mesh.node(dna.top(base).node).sec
                e_base = dna.top(base).id
                length = Norm(dna.top(s_base).pos - dna.top(e_base).pos)

                write(unit, "(i10, a$ )"), n_base, " "//trim(types)
                write(unit, "(a$      )"), ", # of unpaired nts : "//trim(adjustl(Int2Str(count)))
                write(unit, "(a, f5.2$)"), " <-- Total length : ", length
                write(unit, "(a$      )"), ", Edge(sec)) : "
                write(unit, "(a$      )"), trim(adjustl(Int2Str(s_iniL)))//"("//trim(adjustl(Int2Str(s_sec)))//") -> "
                write(unit, "(a$      )"), trim(adjustl(Int2Str(e_iniL)))//"("//trim(adjustl(Int2Str(e_sec)))//"), Sequence : "

                ! Check maximum number of unpaired nucleotides
                if(max_base < count) max_base = count

                do j = 1, count
                    write(unit, "(a$)"), seq(j:j)
                end do
                write(unit, "(a)")
            end if
        end if
    end do
    write(unit, "(a)")
end function Output_Write_Out_Unpaired

! ---------------------------------------------------------------------------------------

! Outputs based on strands and nucleotides
subroutine Output_Write_Out_Strand_Base(mesh, dna, unit)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer, intent(in) :: unit

    integer :: i, j

    write(unit, "(a)"), " Output based on strands"
    write(unit, "(a)"), " The total number of strands : "//trim(adjustl(Int2Str(dna.n_strand)))
    write(unit, "(a)")

    ! Strand based output
    do i = 1, dna.n_strand
        write(unit, "(i10,a$)"), i, " "//trim(dna.strand(i).type1)//" ==>"
        write(unit, "(a, i6$)"), " # of nts : ", dna.strand(i).n_base
        write(unit, "(a, l  )"), ", circular : ", dna.strand(i).b_circular

        ! Print bases numbering
        call space(unit, 15)
        write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(1))))//"->"
        do j = 2, dna.strand(i).n_base - 1
            if(mod(j, 10) == 0) then
                write(unit, "(a9)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            else if(mod(j, 10) == 1) then
                call space(unit, 15)
                write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            else
                write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            end if
        end do
        if(mod(j, 10) == 1) call space(unit, 15)
        write(unit, "(a7)"), trim(adjustl(Int2Str(dna.strand(i).base(dna.strand(i).n_base))))
        write(unit, "(a )")
    end do

    ! DNA base information
    write(unit, "(a)"), " Output based on nucleotides"
    write(unit, "(a)"), " The total number of nucleotides : "//trim(adjustl(Int2Str(dna.n_top)))
    write(unit, "(a)")

    ! node-bp-InitL-sec-up-dn-ac-xo
    do i = 1, dna.n_top
        write(unit, "(i10, a$)"), dna.top(i).id, " nt ==>"
        write(unit, "(a$     )"), trim(dna.strand(dna.top(i).strand).type1)
        write(unit, "(a$     )"), ", seq: "//trim(dna.top(i).seq)
        write(unit, "(a, i6$ )"), ", nde: ", dna.top(i).node

        if(dna.top(i).node /= -1) then
            write(unit, "(a, i6$)"), ", bp: ",   mesh.node(dna.top(i).node).bp
            write(unit, "(a, i3$)"), ", iniL: ", mesh.node(dna.top(i).node).iniL
            write(unit, "(a, i3$)"), ", sec: ",  mesh.node(dna.top(i).node).sec
        else
            write(unit, "(a$)"), ",      UNPAIRED NUCLEOTIDES      "
        end if

        write(unit, "(a, i6$)"), ", up: ",     dna.top(i).up
        write(unit, "(a, i6$)"), ", dn: ",     dna.top(i).dn
        write(unit, "(a, i6$)"), ", xover: ",  dna.top(i).xover
        write(unit, "(a, i6$)"), ", across: ", dna.top(i).across
        write(unit, "(a, i3$)"), ", strand: ", dna.top(i).strand
        !write(unit, "(a, i5)"), ", addr: ", dna.top(i).address
        write(unit, "(a)")
    end do
    write(unit, "(a)")
end subroutine Output_Write_Out_Strand_Base

! ---------------------------------------------------------------------------------------

! Output about staple length
subroutine Output_Write_Out_Staple_Length(dna, unit)
    type(DNAType), intent(in) :: dna
    integer,       intent(in) :: unit

    integer, allocatable :: length_stap(:)
    integer :: i

    ! Allocate and initialize memory
    allocate(length_stap(dna.len_min_stap:dna.len_max_stap))
    length_stap(dna.len_min_stap:dna.len_max_stap) = 0

    ! Find staple length
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "stap") then
            length_stap(dna.strand(i).n_base) = length_stap(dna.strand(i).n_base) + 1
        end if
    end do

    ! Write staple length
    write(unit, "(a)"), "   Edge length     Number"
    write(unit, "(a)"), "   ----------------------"
    do i = dna.len_min_stap, dna.len_max_stap
        if(length_stap(i) /= 0) then
            write(unit, "(2i11)"), i, length_stap(i)
        end if
    end do

    ! Deallocate memory
    deallocate(length_stap)
end subroutine Output_Write_Out_Staple_Length

! ---------------------------------------------------------------------------------------

! Write JSON guide model
subroutine Output_Write_Out_Guide_JSON(prob, geom, bound, mesh)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound
    type(MeshType),  intent(in) :: mesh

    double precision :: length, pos_1(3), pos_2(3), pos_c(3), t1(3)
    integer :: i, j, iter, nbp, add
    logical :: f_axis, f_info
    character(200) :: path

    ! Set option
    f_axis = para_chimera_axis

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=998, file=trim(path)//"_13_json_guide.bild", form="formatted")

    ! Write points for multi lines
    write(998, "(a, 3f6.2)"), ".color ", 0.0d0/255.0d0, 114.0d0/255.0d0, 178.0d0/255.0d0
    do i = 1, geom.n_croP

        write(998, "(a$    )"), ".sphere "
        write(998, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(998, "(1f9.3 )"), 0.3d0
    end do

    ! Write multi lines
    write(998, "(a, 3f6.2)"), ".color ", 0.0d0/255.0d0, 114.0d0/255.0d0, 178.0d0/255.0d0
    do i = 1, geom.n_croL

        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)

        write(998, "(a$   )"), ".cylinder "
        write(998, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.15d0
    end do

    ! Write junctional connection
    write(998, "(a, 3f6.2)"), ".color ", 213.0d0/255.0d0, 94.0d0/255.0d0, 0.0d0/255.0d0
    do i = 1, bound.n_junc
        do j = 1, geom.n_sec*bound.junc(i).n_arm

            pos_1(1:3) = mesh.node(bound.junc(i).conn(j,1)).pos(1:3)
            pos_2(1:3) = mesh.node(bound.junc(i).conn(j,2)).pos(1:3)

            write(998, "(a$   )"), ".cylinder "
            write(998, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.05d0
        end do
    end do

    ! Write scaffold direction
    write(998, "(a)"), ".color red"
    do i = 1, geom.n_croL
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        t1(1:3) = geom.croL(i).t(1, 1:3)

        write(998, "(a$    )"), ".arrow "
        write(998, "(3f8.2$)"), pos_c(1:3)
        write(998, "(3f8.2$)"), pos_c(1:3) + t1(1:3) * 1.5d0
        write(998, "(3f8.2 )"), 0.2d0, 0.6d0, 0.6d0
    end do

    ! Write edge number
    add = 0
    write(998, "(a)"), ".color dark green"
    do i = 1, geom.n_croL
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        add = i / (geom.n_croL / geom.n_sec)

        write(998, "(a$   )"), ".cmov "
        write(998, "(3f9.3)"), pos_c(1:2) + 0.6d0, pos_c(3) - 0.6d0
        write(998, "(a    )"), ".font Helvetica 12 bold"
        write(998, "(i7)"), mod(2*i-1, geom.n_croL) + add - 1
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(998, "(a)"), ".translate 0.0 0.0 0.0"
        write(998, "(a)"), ".scale 0.5"
        write(998, "(a)"), ".color grey"
        write(998, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(998, "(a)"), ".color red"             ! x-axis
        write(998, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(998, "(a)"), ".color blue"            ! y-axis
        write(998, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(998, "(a)"), ".color yellow"          ! z-axis
        write(998, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=998)

    ! ---------------------------------------------
    ! Write the file for Tecplot
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=998, file=trim(path)//"_13_json_guide.dat", form="formatted")

    write(998, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(998, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(998, "(a$)"), 'ZONE F = FEPOINT'
    write(998, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_croP)))
    write(998, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_croL)))
    write(998, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_croP
        write(998, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(998, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_croL
        write(998, "(1i7$)"), geom.croL(i).poi(1)
        write(998, "(1i7 )"), geom.croL(i).poi(2)
    end do

    close(unit=998)
end subroutine Output_Write_Out_Guide_JSON

! ---------------------------------------------------------------------------------------

! Write JSON output
subroutine Output_Write_Out_JSON(prob, geom, mesh, dna, max_unpaired)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: max_unpaired

    ! Conn type data
    type :: ConnType
        integer :: conn(4)  ! [ down sec | down bp | up sec | up bp ]
    end type ConnType

    ! Section type data - Global section => (E-1)*iniE + sec
    type :: SecType
        type(ConnType), allocatable :: scaf(:)
        type(ConnType), allocatable :: stap(:)
        integer :: n_stap_col
        integer :: stap_col(20,2)
    end type SecType

    ! Initial line data
    type :: EdgeType
        integer :: n_sec
        type(SecType), allocatable :: sec(:)
    end type EdgeType

    type(EdgeType), allocatable :: edge(:)
    integer :: i, j, k, min_bp, max_bp, n_edge, width, num, shift
    integer :: c_base, c_node, c_edge, c_sec, c_bp
    integer :: dn_base, dn_node, dn_edge, dn_sec, dn_bp
    integer :: up_base, up_node, up_edge, up_sec, up_bp
    integer :: cup_bp, cdn_bp, col, row, col_shift, row_shift, color(12)
    character(200) :: path

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=999, file=trim(path)//"_14_json.json", form="formatted")

    ! Hex color code
    color( 1) = 13369344    ! #cc0000
    color( 2) = 16204552    ! #f74308
    color( 3) = 16225054    ! #f7931e
    color( 4) = 11184640    ! #aaaa00
    color( 5) = 5749504     ! #57bb00
    color( 6) = 29184       ! #007200
    color( 7) = 243362      ! #03b6a2
    color( 8) = 1507550     ! #1700de
    color( 9) = 7536862     ! #7300de
    color(10) = 12060012    ! #b8056c
    color(11) = 3355443     ! #333333
    color(12) = 8947848     ! #888888

    ! Find maximum and minimum bp ID
    max_bp = mesh.node(1).bp
    min_bp = mesh.node(1).bp
    do i = 2, mesh.n_node
        if(mesh.node(i).bp > max_bp) max_bp = mesh.node(i).bp
        if(mesh.node(i).bp < min_bp) min_bp = mesh.node(i).bp
    end do

    min_bp = min_bp + para_start_bp_ID - 1
    max_bp = max_bp + para_start_bp_ID - 1

    shift = para_start_bp_ID + 21

    ! Possible maximum edge length = max_bp - min_bp + 1 + max_unpaired
    width = (max_bp - min_bp + 1 + max_unpaired) + 2 * para_start_bp_ID + 21
    width = (width / 21) * 21 + 21

    !print *, min_bp, max_bp, para_start_bp_ID, max_unpaired, max_bp - min_bp + 1 + max_unpaired

    ! Allocate and initialize edge data
    n_edge = geom.n_iniL
    allocate(edge(n_edge))

    do i = 1, n_edge
        edge(i).n_sec = geom.n_sec
        allocate(edge(i).sec(edge(i).n_sec))

        do j = 1, edge(i).n_sec
            allocate(edge(i).sec(j).scaf(width))
            allocate(edge(i).sec(j).stap(width))

            edge(i).sec(j).n_stap_col = 0
            edge(i).sec(j).stap_col   = 0

            do k = 1, width
                edge(i).sec(j).scaf(k).conn(1:4) = -1
                edge(i).sec(j).stap(k).conn(1:4) = -1
            end do
        end do
    end do

    ! Strand based loop
    do i = 1, dna.n_strand

        c_base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base

            ! Current base
            c_node = dna.top(c_base).node
            if(c_node /= -1) then
                c_edge = mesh.node(c_node).iniL
                c_sec  = mesh.node(c_node).sec
                c_bp   = mesh.node(c_node).bp + shift
            else
                if(mod(c_sec, 2) == 0) then
                    if(dna.strand(i).type1 == "scaf") c_bp = c_bp + 1
                    if(dna.strand(i).type1 == "stap") c_bp = c_bp - 1
                else
                    if(dna.strand(i).type1 == "scaf") c_bp = c_bp - 1
                    if(dna.strand(i).type1 == "stap") c_bp = c_bp + 1
                end if
            end if

            ! Downward base
            dn_base = dna.top(c_base).dn
            if(dn_base /= -1) then
                dn_node = dna.top(dn_base).node
                if(dn_node /= -1) then
                    dn_edge = mesh.node(dn_node).iniL
                    dn_sec  = mesh.node(dn_node).sec
                    dn_bp   = mesh.node(dn_node).bp + shift
                else
                    if(mod(c_sec, 2) == 0) then
                        if(dna.strand(i).type1 == "scaf") dn_bp = dn_bp + 1
                        if(dna.strand(i).type1 == "stap") dn_bp = dn_bp - 1
                    else
                        if(dna.strand(i).type1 == "scaf") dn_bp = dn_bp - 1
                        if(dna.strand(i).type1 == "stap") dn_bp = dn_bp + 1
                    end if
                end if
                if(dna.strand(i).type1 == "scaf") then
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(1) = (dn_edge - 1) * edge(1).n_sec + dn_sec
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(2) = dn_bp - 1
                else
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(1) = (dn_edge - 1) * edge(1).n_sec + dn_sec
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(2) = dn_bp - 1
                end if
            else if(dn_base == -1 .and. dna.strand(i).type1 == "stap") then

                ! Set color
                edge(c_edge).sec(c_sec+1).n_stap_col = edge(c_edge).sec(c_sec+1).n_stap_col + 1
                num = edge(c_edge).sec(c_sec+1).n_stap_col
                if(num >= 21) then
                    write(0, "(a)"), "Error: Check the number of crossovers per edge"
                    stop
                end if
                edge(c_edge).sec(c_sec+1).stap_col(num,1) = c_bp - 1
                edge(c_edge).sec(c_sec+1).stap_col(num,2) = color(mod(i, 12)+1)
            end if

            ! Upper base
            up_base = dna.top(c_base).up
            if(up_base /= -1) then
                up_node = dna.top(up_base).node
                if(up_node /= -1) then
                    up_edge = mesh.node(up_node).iniL
                    up_sec  = mesh.node(up_node).sec
                    up_bp   = mesh.node(up_node).bp + shift
                else
                    if(mod(c_sec, 2) == 0) then
                        if(dna.strand(i).type1 == "scaf") up_bp = up_bp + 1
                        if(dna.strand(i).type1 == "stap") up_bp = up_bp - 1
                    else
                        if(dna.strand(i).type1 == "scaf") up_bp = up_bp - 1
                        if(dna.strand(i).type1 == "stap") up_bp = up_bp + 1
                    end if
                end if
                if(dna.strand(i).type1 == "scaf") then
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(3) = (up_edge - 1) * edge(1).n_sec + up_sec
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(4) = up_bp - 1
                else
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(3) = (up_edge - 1) * edge(1).n_sec + up_sec
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(4) = up_bp - 1
                end if
            end if

            ! Update base
            c_base = dna.Top(c_base).up
        end do
    end do

    ! Print JSON-style data structure
    write(999, "(a)"), "{"
    write(999, "(a)"), '"vstrands":'
    write(999, "(a)"), '['

    row_shift = 0
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            write(999, "(a)"), '{'

            ! Skip
            write(999, "(a$)"), '"skip":['
            do k = 1, width - 1
                write(999, "(a$)"), "0,"
            end do
            write(999, "(a)"), "0],"

            ! Stap
            write(999, "(a$)"), '"stap":['
            do k = 1, width
                write(999, "(a$)"), "["//&
                    trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(1))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(2))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(3))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(4))))

                if(k /= width) write(999, "(a$)"), "]"
                if(k == width) write(999, "(a )"), "]],"
            end do

            ! Scaffold loop
            write(999, "(a)"), '"scafLoop":[],'

            ! Staple colors
            write(999, "(a$)"), '"stap_colors":['
            do k = 1, edge(i).sec(j).n_stap_col
                write(999, "(a$)"), '['//trim(adjustl(Int2Str(edge(i).sec(j).stap_col(k,1))))
                write(999, "(a$)"), ','//trim(adjustl(Int2Str(edge(i).sec(j).stap_col(k,2))))
                
                if(k /= edge(i).sec(j).n_stap_col) write(999, "(a$)"), '],'
                if(k == edge(i).sec(j).n_stap_col) write(999, "(a$)"), ']'
            end do
            write(999, "(a )"), '],'

            ! Scaffold loop
            write(999, "(a)"), '"stapLoop":[],'

            ! Scaffold routing
            write(999, "(a$)"), '"scaf":['
            do k = 1, width
                write(999, "(a$)"), "["//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(1))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(2))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(3))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(4))))

                if(k /= width) write(999, "(a$)"), "]"
                if(k == width) write(999, "(a )"), "]],"
            end do

            ! Cross-section position
            col_shift = mod(i, 15)
            if(col_shift == 1 .and. j == 1) row_shift = row_shift + 2
            if(col_shift == 0) col_shift = 15
            col_shift = 2*col_shift

            num = (i - 1) * edge(i).n_sec + j
            col = col_shift - geom.sec.posR(j)
            row = row_shift - geom.sec.posC(j)

            write(999, "(a)"), '"num":'//trim(adjustl(Int2Str(num - 1)))//","
            write(999, "(a)"), '"col":'//trim(adjustl(Int2Str(col)))//","
            write(999, "(a)"), '"row":'//trim(adjustl(Int2Str(row)))//","

            ! Loop
            write(999, "(a$)"), '"loop":['
            do k = 1, width - 1
                write(999, "(a$)"), "0,"
            end do
            write(999, "(a)"), "0]"

            if(i == n_edge .and. j == edge(i).n_sec) then
                write(999, "(a)"), '}'
            else
                write(999, "(a)"), '},'
            end if
        end do
    end do

    write(999, "(a)"), '],'
    write(999, "(a)"), '"name":"2D lattice design by PERDIX-OPEN"'
    write(999, "(a)"), "}"

    ! Deallocate edge data
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            deallocate(edge(i).sec(j).scaf)
            deallocate(edge(i).sec(j).stap)
        end do
        deallocate(edge(i).sec)
    end do
    deallocate(edge)
end subroutine Output_Write_Out_JSON

! ---------------------------------------------------------------------------------------

! Write information related with basepair
subroutine Output_Write_Basepair(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    character(200) :: path
    integer :: i

    if(para_write_801 == .false.) return

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=801, file=trim(path)//"_basepair.txt", form="formatted")

    ! The number of base pairs
    write(801, "(i7)"), mesh.n_node

    ! Write information related with base pair
    do i = 1, mesh.n_node
        write(801, "(i7$   )"), i
        write(801, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(1, 1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(2, 1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(3, 1:3)
        write(801, "(i7$   )"), dna.base_scaf(i).xover
        write(801, "(i7    )"), dna.base_stap(i).xover
    end do
    write(801, "(a)")

    ! The number connectivity of base pairs
    write(801, "(i7)"), mesh.n_ele

    ! Write the connectivity of basepair
    do i = 1, mesh.n_ele
        write(801, "(i7, 2i8)"), i, mesh.ele(i).cn(1:2)
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Write information related with basepair"
    end do

    close(unit=801)
end subroutine Output_Write_Basepair

! ---------------------------------------------------------------------------------------

! Write information related with base
subroutine Output_Write_Base(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    character(200) :: path
    integer :: i, id, id_up, id_down

    if(para_write_802 == .false.) return

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=802, file=trim(path)//"_base.txt", form="formatted")

    ! The total number of bases
    write (802, "(i7)"), dna.n_base_scaf + dna.n_base_stap

    ! For scaffold strand
    do i = 1, dna.n_top

        ! ID, up ID, down ID, across ID, position vector
        write(802, "(i7$   )"), dna.top(i).id
        write(802, "(i7$   )"), dna.top(i).up
        write(802, "(i7$   )"), dna.top(i).dn
        write(802, "(i7$   )"), dna.top(i).across
        write(802, "(3f10.4)"), dna.top(i).pos(1:3)
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Write information related with bases"
    end do

    close(unit=802)
end subroutine Output_Write_Base

! ---------------------------------------------------------------------------------------

! Write cndo file for PDB atom generation and CanDo simulation
subroutine Output_Write_CanDo(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    integer, allocatable :: node_nt(:,:)

    double precision :: min_pos(3), pos(3), fscale
    integer :: i, j, base, count, strand
    character(200) :: path

    if(para_write_803 == .false.) return

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=803, file=trim(path)//"_15_cndo.cndo", form="formatted")

    write(803, "(a)"), '"CanDo (.cndo) file format version 1.0"'
    write(803, "(a)")

    ! Set viewpoint and RGB color
    if(para_type_cndo /= 1) then
        write(803, "(a)"), trim(para_fig_view)
        write(803, "(a)")

        write(803, "(f8.2)"), prob.scale
        write(803, "(a)")

        write(803, "(3i)"), prob.color(1:3)
        write(803, "(a)")
    end if

    ! For dnatop data that is defined by bases
    if(para_type_cndo == 1) then

        ! For cndo 1.0
        write(803, "(a)"), 'dnaTop,id,up,down,across,seq'
        do i = 1, dna.n_top

            ! ID, up, down, across and sequence
            write(803, "(a$)"), &
                trim(adjustl(Int2Str(i)))                 //","//&
                trim(adjustl(Int2Str(dna.top(i).id)))     //","//&
                trim(adjustl(Int2Str(dna.top(i).dn)))     //","//&      ! cndo convention is opposite
                trim(adjustl(Int2Str(dna.top(i).up)))     //","//&      ! cndo convention is opposite
                trim(adjustl(Int2Str(dna.top(i).across))) //","

            write(803, "(a)"), trim(dna.top(i).seq)
        end do
    else

        ! For updated cndo with strand type
        write(803, "(a)"), 'dnaTop,id,up,down,across,seq,strand'
        do i = 1, dna.n_top

            ! ID, up, down, across and sequence
            write(803, "(a$)"), &
                trim(adjustl(Int2Str(i)))                 //","//&
                trim(adjustl(Int2Str(dna.top(i).id)))     //","//&
                trim(adjustl(Int2Str(dna.top(i).dn)))     //","//&      ! cndo convention is opposite
                trim(adjustl(Int2Str(dna.top(i).up)))     //","//&      ! cndo convention is opposite
                trim(adjustl(Int2Str(dna.top(i).across))) //","//&
                trim(dna.top(i).seq)                      //","

            strand = dna.top(i).strand
            if(dna.strand(strand).type1 == "scaf") then
                write(803, "(a)"), trim("0")
            else if(dna.strand(strand).type1 == "stap") then
                write(803, "(a)"), trim("1")
            end if
        end do
    end if
    write(803, "(a)")

    ! Pre-calculation
    min_pos(1:3) = mesh.node(1).pos(1:3)
    do i = 1, mesh.n_node
        if(min_pos(1) > mesh.node(i).pos(1)) min_pos(1) = mesh.node(i).pos(1)
        if(min_pos(2) > mesh.node(i).pos(2)) min_pos(2) = mesh.node(i).pos(2)
        if(min_pos(3) > mesh.node(i).pos(3)) min_pos(3) = mesh.node(i).pos(3)
    end do

    fscale = 10.0d0
    if(min_pos(1) < 0 .and. min_pos(1) < -99.0d0) then
        min_pos(1) = (min_pos(1) + 80.0d0)
    else
        min_pos(1) = 0.0d0
    end if

    if(min_pos(2) < 0 .and. min_pos(2) < -99.0d0) then
        min_pos(2) = (min_pos(2) + 80.0d0)
    else
        min_pos(2) = 0.0d0
    end if

    if(min_pos(3) < 0 .and. min_pos(3) < -99.0d0) then
        min_pos(3) = (min_pos(3) + 80.0d0)
    else
        min_pos(3) = 0.0d0
    end if

    ! The number of base pairs
    write(803, "(a)"), 'dNode,"e0(1)","e0(2)","e0(3)"'
    do i = 1, mesh.n_node

        pos(:) = (mesh.node(i).pos(:)-min_pos(:))*fscale

        write(803, "(a)"), &
            trim(adjustl(Int2Str(i)))//","//&
            trim(adjustl(Dble2Str(pos(1)))) //","//&
            trim(adjustl(Dble2Str(pos(2)))) //","//&
            trim(adjustl(Dble2Str(pos(3))))

        if( (pos(1) > 9999.0d0 .or. pos(1) < -999.0d0) .or. &
            (pos(2) > 9999.0d0 .or. pos(2) < -999.0d0) .or. &
            (pos(3) > 9999.0d0 .or. pos(3) < -999.0d0) ) then

            ! Print error for generating atomic model
            write(0, "(a, f10.3)"), "x-pos : ", pos(1)
            write(0, "(a, f10.3)"), "y-pos : ", pos(2)
            write(0, "(a, f10.3)"), "z-pos : ", pos(3)
            stop
        end if
    end do
    write(803, "(a)")

    ! Write orientations
    write(803, "(a$)"), 'triad,"e1(1)","e1(2)","e1(3)",'
    write(803, "(a$)"), '"e2(1)","e2(2)","e2(3)",'
    write(803, "(a )"), '"e3(1)","e3(2)","e3(3)"'

    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i)))                       //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 3)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 3)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 3))))
    end do
    write(803, "(a)")

    ! nt1, and nt2 based on basepair
    allocate(node_nt(mesh.n_node, 2))
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                if(dna.top(base).across /= -1) then
                    node_nt(dna.top(base).node, 1) = dna.top(base).id
                    node_nt(dna.top(base).node, 2) = dna.top(base).across
                end if
            end do
        end if
    end do

    ! Set base ID
    write(803, "(a)"), 'id_nt,id1,id2'
    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i))) //","//&
            trim(adjustl(Int2Str(node_nt(i, 1))))//","//&
            trim(adjustl(Int2Str(node_nt(i, 2))))
    end do
    deallocate(node_nt)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Write CanDo file (*.cndo)"
    end do

    close(unit=803)
end subroutine Output_Write_CanDo

! ---------------------------------------------------------------------------------------

! Write cndo file for PDB atom generation and CanDo simulation
subroutine Output_Write_CanDo_New(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    integer, allocatable :: node_nt(:,:)
    integer :: i, j, base, strand, types
    character(200) :: path

    if(para_write_803 == .false.) return

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=803, file=trim(path)//"_15_cndo.cndo", form="formatted")

    ! For dnatoop data that is defined by bases
    write(803, "(i10)"), dna.n_top
    do i = 1, dna.n_top

        ! ID, up ID, down ID, across ID and sequence
        write(803, "(i10$)"), i
        write(803, "(i10$)"), dna.top(i).id
        write(803, "(i10$)"), dna.top(i).up
        write(803, "(i10$)"), dna.top(i).dn
        write(803, "(i10$)"), dna.top(i).across
        write(803, "(a10$)"), dna.top(i).seq

        strand = dna.top(i).strand
        if(dna.strand(strand).type1 == "scaf") then
            write(803, "(i10)"), 0
        else if(dna.strand(strand).type1 == "stap") then
            write(803, "(i10)"), 1
        end if
    end do
    write(803, "(a)")

    ! The number of base pairs
    write(803, "(i10)"), mesh.n_node
    do i = 1, mesh.n_node
        write(803, "(i10$  )"), i
        write(803, "(f15.6$)"), mesh.node(i).pos(1)*10.0d0
        write(803, "(f15.6$)"), mesh.node(i).pos(2)*10.0d0
        write(803, "(f15.6 )"), mesh.node(i).pos(3)*10.0d0
    end do
    write(803, "(a)")

    ! e1 : the mojor groove
    ! e2 : the preferred nucleotide
    ! e3 : along the duplex axis towards the 3'-direction of
    !      the strand with the preferred nucleotide
    write(803, "(i10)"), mesh.n_node
    do i = 1, mesh.n_node
        write(803, "(i10$  )"), i
        write(803, "(3f15.6$)"), mesh.node(i).ori(1, 1:3)
        write(803, "(3f15.6$)"), mesh.node(i).ori(2, 1:3)
        write(803, "(3f15.6 )"), mesh.node(i).ori(3, 1:3)
    end do
    write(803, "(a)")

    ! nt1, and nt2 based on basepair
    allocate(node_nt(mesh.n_node, 2))
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                if(dna.top(base).across /= -1) then
                    node_nt(dna.top(base).node, 1) = dna.top(base).id
                    node_nt(dna.top(base).node, 2) = dna.top(base).across
                end if
            end do
        end if
    end do

    ! Set base ID
    write(803, "(a)"), 'id_nt,id1,id2'
    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i))) //","//&
            trim(adjustl(Int2Str(node_nt(i, 1))))//","//&
            trim(adjustl(Int2Str(node_nt(i, 2))))
    end do
    deallocate(node_nt)

    close(unit=803)
end subroutine Output_Write_CanDo_New

! ---------------------------------------------------------------------------------------

! Write basepair and nucleotide information
subroutine Output_Write_DNA_Info(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    integer :: i, id, id_up, id_down
    character(200) :: path

    ! open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=100, file=trim(path)//"_dnaInfo.dat", form="formatted")

    ! the total number of nucleotide
    write (100, "(i8)"), dna.n_top

    ! for scaffold strand
    do i = 1, dna.n_top

        ! id, up id, down id, across id, position vector
        write(100, "(i8$   )"), i
        write(100, "(i8$   )"), dna.top(i).id
        write(100, "(i8$   )"), dna.top(i).up
        write(100, "(i8$   )"), dna.top(i).dn
        write(100, "(i8$   )"), dna.top(i).across
        write(100, "(3f10.4)"), dna.top(i).pos(1:3)
    end do

    close(unit=100)
end subroutine Output_Write_DNA_Info

! ---------------------------------------------------------------------------------------

! Write TecPlot input file
subroutine Output_Write_TecPlot(prob, mesh)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh

    character(200) :: path
    integer :: i, k

    if(para_write_804 == .false.) return

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=804, file=trim(path)//"_tecplot.dat", form="formatted")

    ! Loop for 3 orientation vectors
    do k = 1, 3

        ! For Tecplot output
        write(804, "(a     )"), "variables = 'X', 'Y', 'Z', 'e3', 'e2', 'e1' "
        write(804, "(a, i7$)"), "ZONE F=FEPOINT, N=", mesh.n_node
        write(804, "(a, i7$)"), ", E=",               mesh.n_ele
        write(804, "(a     )"), ", ET=QUADRILATERAL"

        ! Write nodal position vector
        do i = 1, mesh.n_node
            write(804, "(3f10.4$)"), mesh.node(i).pos(1:3)
            write(804, "(3f10.4 )"), mesh.node(i).ori(k, 1:3)
        end do

        ! Element connectivity
        do i = 1, mesh.n_ele
            write(804, "(2i6$)"), mesh.ele(i).cn(1:2)
            write(804, "(2i6 )"), mesh.ele(i).cn(1:2)
        end do
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Write TecPlot input file"
    end do

    close(unit=804)
end subroutine Output_Write_TecPlot

! ---------------------------------------------------------------------------------------

! Write ADINA input file
subroutine Output_Write_ADINA(prob, mesh)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh

    character(200) :: path
    integer :: i, time(8)

    if(para_write_805 == .false.) return

    call date_and_time(VALUES=time)

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=805, file=trim(path)//"_adina.in", form="formatted")

    ! Write of ADINA infile
    write(805, "(a     )"), "*"
    write(805, "(a     )"), "* Command file created from PERDIX for command import"
    write(805, "(a     )"), "*"
    write(805, "(a$    )"), "*--- Command file created "
    write(805, "(i4, a$)"), time(2), "/"
    write(805, "(i4, a$)"), time(3), "/"
    write(805, "(i4, a$)"), time(1), ", "
    write(805, "(i2, a$)"), time(5), ":"
    write(805, "(i2, a$)"), time(6), ":"
    write(805, "(i2, a$)"), time(7), " ---*"
    write(805, "(a     )"), "*--- for ADINA: AUI version 9.0.6 ---*"
    write(805, "(a     )"), "*"
    write(805, "(a     )"), "DATABASE NEW SAVE=NO PROMPT=NO"    ! That creates a new database
    write(805, "(a     )"), "FEPROGRAM ADINA"                   ! For displacement and stress analysis
    write(805, "(a     )"), "COORDINATES POINT SYSTEM=0"        ! Definition of coordinates for geometry points
    write(805, "(a     )"), "@CLEAR"
    write(805, "(a     )")

    ! Write point information
    do i = 1, mesh.n_node
        write(805, "(i, 3f13.5, i)"), i, mesh.node(i).pos(:), 0
    end do
    write(805, "(a)") "@"

    ! Write line information
    do i = 1, mesh.n_ele
        write(805, "(a, i$)"), "LINE STRAIGHT NAME=", i
        write(805, "(a, i$)"), " P1=", mesh.ele(i).cn(1)
        write(805, "(a, i )"), " P2=", mesh.ele(i).cn(2)
        write(805, "(a    )"), "*"
    end do

    ! Set material properties
    write(805, "(a)"), "*"
    write(805, "(a)"), "EGROUP BEAM NAME=1 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=1 RINT=5,"
    write(805, "(a)"), "     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,"
    write(805, "(a)"), "     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,"
    write(805, "(a)"), "     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,"
    write(805, "(a)"), "     BOLT-TOL=0.00000000000000 DESCRIPT='NONE' SECTION=1,"
    write(805, "(a)"), "     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,"
    write(805, "(a)"), "     TDEATH=0.00000000000000 SPOINT=2 BOLTFORC=0.00000000000000,"
    write(805, "(a)"), "     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,"
    write(805, "(a)"), "     WARP=NO" ! ENDRELEA=ACCURATE"

    ! Geometry lines are assigned a number of subdivisions
    write(805, "(a)"), "*"
    write(805, "(a)"), "SUBDIVIDE LINE NAME=1 MODE=DIVISIONS NDIV=1 RATIO=1.00000000000000,"
    write(805, "(a)"), "      PROGRESS=GEOMETRIC CBIAS=NO"
    write(805, "(a)"), "@CLEAR"

    ! Element connectivity
    do i = 1, mesh.n_ele
        write(805, "(i7)") i
    end do
    write(805, "(a)") "@"

    ! Mesh generation
    write(805, "(a)"), "*"
    write(805, "(a)"), "GLINE NODES=2 AUXPOINT=0 NCOINCID=ENDS NCENDS=12,"
    write(805, "(a)"), "     NCTOLERA=1.00000000000000E-05 SUBSTRUC=0 GROUP=1 MIDNODES=CURVED,"
    write(805, "(a)"), "     XO=0.00000000000000 YO=1.00000000000000 ZO=0.00000000000000,"
    write(805, "(a)"), "     XYZOSYST=SKEW"
    write(805, "(a)"), "@CLEAR"
    do i = 1, mesh.n_ele
        write(805, "(i7)") i
    end do
    write(805, "(a)") "@"

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Write ADINA input file (*.in)"
        write(i, "(a)")
    end do

    close(unit = 805)
end subroutine Output_Write_ADINA

! ---------------------------------------------------------------------------------------

! Make route step using Chimera
subroutine Output_Make_Route_Step(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    integer :: i, j, id, step
    logical :: results
    character(10) :: str

    ! Make folder to save route step Chimera files
    results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_bild\")

    ! Go to the first base of the strand
    id = dna.strand(1).base(1)
    do
        if(dna.top(id).dn == -1) exit
        id = dna.top(id).dn
    end do

    step = dna.n_base_scaf / para_n_route_step
    do i = 1, para_n_route_step - 1

        ! Go to upper base
        do j = 1, step
            id = dna.base_scaf(id).up
        end do

        write(str, "(i0)"), i
        call Output_Chimera_Route_Step(prob, mesh, dna, id, str)
    end do

    ! go to the end of the strand
    do
        if(dna.top(id).up == -1) exit
        id = dna.top(id).up
    end do

    ! The last route step
    write(str, "(i0)"), para_n_route_step
    call Output_Chimera_Route_Step(prob, mesh, dna, id, str)

    ! Make figures from route step
    call Output_Figure_Route_Step(prob, para_n_route_step)
end subroutine Output_Make_Route_Step

! ---------------------------------------------------------------------------------------

! Write Chimera for route progress
subroutine Output_Chimera_Route_Step(prob, mesh, dna, stop_base, str)
    type(ProbType),  intent(in) :: prob
    type(MeshType),  intent(in) :: mesh
    type(DNAType),   intent(in) :: dna
    integer,         intent(in) :: stop_base
    character(10),   intent(in) :: str

    double precision :: pos_1(3), pos_2(3)
    integer :: i, id, up, base_1, base_2, node_1, node_2
    logical :: f_axis
    character(200) :: path

    ! Set option
    f_axis = para_chimera_axis

    ! Open file
    path = trim(prob.path_work1)//"step_bild\"//trim(prob.name_file)
    open(unit=806, file=trim(path)//"_rstep_"//trim(str)//".bild", form="formatted")

    ! Go to the first base of the strand
    id = dna.strand(1).base(1)
    do
        if(dna.top(id).dn == -1) exit
        id = dna.top(id).dn
    end do

    write(806, "(a)"), ".color steel blue"
    do
        ! If the ID meets the stop base number, terminate the current subroutine
        if(id == stop_base) then
            !exit
            write(806, "(a)"), ".color grey"
        end if

        ! If the upper ID is negative, terminate loop
        if(dna.top(id).up == -1) then
            exit
        end if

        ! Set base ID
        base_1 = dna.top(id).id
        base_2 = dna.top(id).up

        ! Set node number
        node_1 = dna.top(base_1).node
        node_2 = dna.top(base_2).node

        ! Set position vector
        pos_1(1:3) = mesh.node(node_1).pos(1:3)
        pos_2(1:3) = mesh.node(node_2).pos(1:3)

        write(806, "(a$    )"), ".cylinder "
        write(806, "(3f9.3$)"), pos_1(1:3)
        write(806, "(3f9.3$)"), pos_2(1:3)
        write(806, "(1f9.3 )"), 0.15d0

        ! Update ID number to upward base
        id = dna.top(id).up
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(806, "(a)"), ".translate 0.0 0.0 0.0"
        write(806, "(a)"), ".scale 0.5"
        write(806, "(a)"), ".color grey"
        write(806, "(a)"), ".sphere 0 0 0 0.5"      ! center
        write(806, "(a)"), ".color red"             ! x-axis
        write(806, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(806, "(a)"), ".color blue"            ! y-axis
        write(806, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(806, "(a)"), ".color yellow"          ! z-axis
        write(806, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if

    close(unit=806)
end subroutine Output_Chimera_Route_Step

! ---------------------------------------------------------------------------------------

! Make figures from route step
subroutine Output_Figure_Route_Step(prob, n_progress)
    type(ProbType), intent(in) :: prob
    integer,  intent(in) :: n_progress

    character(200) :: path, str
    logical :: results
    integer :: i

    ! Make step directories for PNG files
    if(para_fig_bgcolor == "white") then
        ! --------------------------------------------------
        ! White background
        ! --------------------------------------------------
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_YZ\")
        else if(para_fig_view == "xyz" .or. para_fig_view == "all" ) then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XYZ\")
        end if
    else if(para_fig_bgcolor == "black") then
        ! --------------------------------------------------
        ! Black background
        ! --------------------------------------------------
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_YZ\")
        else if(para_fig_view == "xyz" .or. para_fig_view == "all" ) then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XYZ\")
        end if
    else if(para_fig_bgcolor == "all") then
        ! --------------------------------------------------
        ! White and black background
        ! --------------------------------------------------
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_YZ\")
        else if(para_fig_view == "xyz" .or. para_fig_view == "all" ) then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_white_XYZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"step_fig_black_XYZ\")
        end if
    end if

    ! Make batch file to make figures
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=807, file=trim(path)//"_rstep.bat", form="formatted")

    write(807, "(a)"), "@echo off"
    do i = 1, n_progress

        write(str, "(i0)"), i
        str = trim(prob.name_file)//trim("_rstep_"//trim(str))

        ! Make Python script to save imanges
        if(i == 1) then
            ! First Python script
            call Chimera_Figure_Route_Step(prob, str, 1)
        else if(i == n_progress) then
            ! Last Python script
            call Chimera_Figure_Route_Step(prob, str, 2)
        else
            call Chimera_Figure_Route_Step(prob, str, 0)
        end if
    end do

    ! Write python file into batch script
    path = " "//trim(prob.path_work1)//trim(prob.name_file)//"_rstep.py"
    write(807, "(a)"), '"'//trim(prob.path_chimera)//'" --script'//path

    ! Run batch file
    path = trim(prob.path_work1)//trim(prob.name_file)
    results = SYSTEMQQ("start cmd /C call "//trim(path)//"_rstep.bat")

    close (unit=807)
end subroutine Output_Figure_Route_Step

! ---------------------------------------------------------------------------------------

! Write sequence based on cross-sectional line
subroutine Output_Write_Sequence_CroL(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    character, allocatable :: seq(:)
    integer :: i, j, base, node
    character(200) :: path

    if(para_write_808 == .false.) return

    ! Allocate memory
    allocate(seq(mesh.n_node))

    ! Find sequence from top data
    do i = 1, dna.n_strand

        ! For staple strand
        if(dna.strand(i).type1 /= "stap") cycle
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            node = dna.top(base).node

            if(node /= -1) then
                seq(node) = dna.top(base).seq
            end if
        end do
    end do

    ! Open files
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=808, file=trim(path)//"_seq_line.txt", form="formatted")

    ! Write information on sequence
    do i = 1, mesh.n_node
        write(808, "(i5$)"), mesh.node(i).croL
        write(808, "(a$ )"), " multi line -> node ID : "
        write(808, "(i5$)"), mesh.node(i).id
        write(808, "(a$ )"), ", BP ID : "
        write(808, "(i5$)"), mesh.node(i).bp
        write(808, "(a$ )"), ", sequence for staple : "
        write(808, "(a  )"), seq(i)
    end do

    ! Deallocate memory
    deallocate(seq)

    close(unit=808)
end subroutine Output_Write_Sequence_CroL

! ---------------------------------------------------------------------------------------

end module Output