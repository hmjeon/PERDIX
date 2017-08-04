!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Output
!
!                                                                    Updated : 2017/03/27
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

    use Mani
    use Para
    use Math

    use Chimera
    use Tecplot

    implicit none

    public  Output_Generation

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
subroutine Output_Generation(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

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
    open(unit=803, file=trim(path)//"_13_cndo.cndo", form="formatted")

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
            if(dna.strand(strand).types == "scaf") then
                write(803, "(a)"), trim("0")
            else if(dna.strand(strand).types == "stap") then
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
        if(dna.strand(i).types == "scaf") then
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
    open(unit=803, file=trim(path)//"_13_cndo.cndo", form="formatted")

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
        if(dna.strand(strand).types == "scaf") then
            write(803, "(i10)"), 0
        else if(dna.strand(strand).types == "stap") then
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
        if(dna.strand(i).types == "scaf") then
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
        if(dna.strand(i).types /= "stap") cycle
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
        write(808, "(a$ )"), " th sectional line -> node ID : "
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