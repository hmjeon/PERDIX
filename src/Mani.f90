!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Mani
!
!                                                                    Updated : 2017/03/27
!
! Comments: This module is for manipulation of data structures used.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Mani

    use Data_Geom
    use Data_Mesh
    use Data_DNA
    use Data_Prob

    use Para

    implicit none

    public Mani_Set_View_Color
    public Space
    public Mani_To_Upper
    public Mani_Progress_Bar

    ! Manipulate data structure
    public Mani_Init_LineType
    public Mani_Allocate_SecType
    public Mani_Init_SecType
    public Mani_Init_MeshType
    public Mani_Init_Node
    public Mani_Init_Ele
    public Mani_Copy_NodeType
    public Mani_Copy_EleType
    public Mani_Init_BaseType
    public Mani_Init_StrandType
    public Mani_Copy_BaseType
    public Mani_Go_Start_Base
    public Mani_Go_Base

contains

! ---------------------------------------------------------------------------------------

! Set geometric type and view
subroutine Mani_Set_View_Color(prob, color, view, scale, size, move_x, move_y)
    type(ProbType), intent(inout) :: prob
    double precision, intent(in)  :: scale, size, move_x, move_y
    character(len=*), intent(in)  :: view
    integer, intent(in) :: color(3)

    prob.color  = color
    prob.scale  = scale     ! Atomic model
    prob.size   = size      ! Cylindrical model
    prob.move_x = move_x    ! Cylindrical model
    prob.move_y = move_y    ! Cylindrical model

    ! Set view points
    para_fig_view = trim(adjustl(view))
end subroutine Mani_Set_View_Color

! ---------------------------------------------------------------------------------------

! Print space
subroutine Space(n_unit, num)
    integer, intent(in) :: n_unit, num

    integer :: i

    do i = 1, num
    write(n_unit, "(a$)"), " "
    end do
end subroutine Space

! ---------------------------------------------------------------------------------------

! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (1 Mar 2016)
function Mani_To_Upper(strIn) result(strOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut

    integer :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do
end function Mani_To_Upper

! ---------------------------------------------------------------------------------------   
    
! Print progress bar
subroutine Mani_Progress_Bar(index, max)
    integer, intent(in) :: index, max

    integer :: i, step

    character(len=56) :: bar = "           * Progressing.... [                    ] ???%"
    step = index * 10 / max

    write(unit=bar(53:55), fmt = "(i3)"), 10*step

    ! Initialize progress bar
    if(index == 1) then
        do i = 1, 10
            bar(30+2*i  :30+2*i  ) = " "
            bar(30+2*i-1:30+2*i-1) = " "
        end do
    end if

    ! Print star
    do i = 1, step
        if(i == step .and. i /= 10) then
            bar(30+2*i  :30+2*i  ) = ">"
            bar(30+2*i-1:30+2*i-1) = "-"
        else
            bar(30+2*i  :30+2*i  ) = "-"
            bar(30+2*i-1:30+2*i-1) = "-"
        end if
    end do

    ! Print the progress bar, char(13) - "\n"
    write(unit=0, fmt="(a1, a56, $)"), char(13), bar

    if(step /= 10) then
        flush(unit=0)
    else
        write(unit=0, fmt=*)
    end if
end subroutine Mani_Progress_Bar

! ---------------------------------------------------------------------------

! Initialize line data type
subroutine Mani_Init_LineType(line, n_line)
    type(LineType), allocatable, intent(inout) :: line(:)
    integer, intent(in) :: n_line

    integer :: i, j

    do i = 1, n_line
        line(i).iniL = -1   ! Inital line number
        line(i).sec  = -1   ! Cross-section ID

        do j = 1, 2
            line(i).poi(j)       = -1   ! Point connectivity
            line(i).neiP(j, 1:2) = -1   ! Neighbor points
            line(i).neiL(j, 1:2) = -1   ! Neighbor lines
        end do

        do j = 1, 3
            line(i).t(j, 1:3) = 0.0d0   ! Local vector at the center
        end do
    end do
end subroutine Mani_Init_LineType

! ---------------------------------------------------------------------------

! Allocate section type data
subroutine Mani_Allocate_SecType(sec, n_sec)
    type(SecType), intent(inout) :: sec
    integer,       intent(in)    :: n_sec

    allocate(sec.id(n_sec))     ! Section ID, start from 0
    allocate(sec.posR(n_sec))   ! Row position number, start from 1
    allocate(sec.posC(n_sec))   ! Column position number, start from 1
    allocate(sec.conn(n_sec))   ! Connetivity for self connection route
end subroutine Mani_Allocate_SecType

! ---------------------------------------------------------------------------

! Initialize section type data
subroutine Mani_Init_SecType(sec, n_sec, types)
    type(SecType), intent(inout) :: sec
    integer,       intent(in)    :: n_sec
    character(*),  intent(in)    :: types

    integer :: i, m_value

    m_value = 9999999

    sec.types    = types
    sec.maxR     = -m_value     ! Maximum row number
    sec.minR     =  m_value     ! Minimum row number
    sec.maxC     = -m_value     ! Maximum column number
    sec.minC     =  m_value     ! Minimum column number
    sec.n_row    = -1           ! Size of row
    sec.n_col    = -1           ! Size of column
    !sec.ref_row  = -1           ! Reference row to set t-axis
    !sec.ref_maxC = -1           ! Maximum column in reference row
    !sec.ref_minC = -1           ! Minimum column in reference row

    do i = 1, n_sec
        sec.id(i)   = -1        ! Section ID
        sec.posR(i) = -1        ! Row position number
        sec.posC(i) = -1        ! Column position number
        sec.conn(i) = -1        ! Connetivity for self connection route
    end do
end subroutine Mani_Init_SecType

! ---------------------------------------------------------------------------

! Initialize mesh type data
subroutine Mani_Init_MeshType(mesh)
    type(MeshType), intent(inout) :: mesh

    integer :: i, j

    do i = 1, mesh.n_node
        mesh.node(i).id      = -1   ! Node ID
        mesh.node(i).bp      = -1   ! Base pair ID
        mesh.node(i).up      = -1   ! Upward ID
        mesh.node(i).dn      = -1   ! Downward ID
        mesh.node(i).sec     = -1   ! Section ID
        mesh.node(i).iniL    = -1   ! Initial line
        mesh.node(i).croL    = -1   ! Cross-section line
        mesh.node(i).beveled = -1   ! Beveled node
        mesh.node(i).conn    = -1   ! -1 - no-connection, 1 - neighbor, 2 - self, 3 - modified neighbor, 4- modified self
        mesh.node(i).ghost   = -1   ! Ghost node

        ! Position and orientation vector
        do j = 1, 3
            mesh.node(i).pos(j)      = 0.0d0
            mesh.node(i).ori(j, 1:3) = 0.0d0
        end do
    end do

    ! Initialize connectivity
    do i = 1, mesh.n_ele
        mesh.ele(i).cn(1:2) = -1
    end do
end subroutine Mani_Init_MeshType

! ---------------------------------------------------------------------------

! Initialize node type data
subroutine Mani_Init_Node(node, n_node)
    type(NodeType), intent(inout) :: node(:)
    integer,        intent(in)    :: n_node

    integer :: i, j

    do i = 1, n_node
        node(i).id      = -1    ! Node ID
        node(i).bp      = -1    ! Base pair ID
        node(i).up      = -1    ! Upward ID
        node(i).dn      = -1    ! Downward ID
        node(i).sec     = -1    ! Section ID
        node(i).iniL    = -1    ! Initial line
        node(i).croL    = -1    ! Cross-section line
        node(i).beveled = -1    ! Beveled node
        node(i).conn    = -1    ! -1 - no-connection, 1 - neighbor, 2 - self, 3 - modified neighbor, 4 - modified self
        node(i).ghost   = -1    ! Ghost node

        ! Position and orientation vector
        do j = 1, 3
            node(i).pos(j)      = 0.0d0
            node(i).ori(j, 1:3) = 0.0d0
        end do
    end do
end subroutine Mani_Init_Node

! ---------------------------------------------------------------------------

! Initialize mesh type data
subroutine Mani_Init_Ele(ele, n_ele)
    type(EleType), intent(inout) :: ele(:)
    integer,       intent(in)    :: n_ele

    integer :: i

    ! Initialize connectivity
    do i = 1, n_ele
        ele(i).cn(1:2) = -1
    end do
end subroutine Mani_Init_Ele

! ---------------------------------------------------------------------------

! Copy node data from ori to copy with size num
subroutine Mani_Copy_NodeType(ori, copy, num)
    type(NodeType), intent(in)    :: ori(:)
    type(NodeType), intent(inout) :: copy(:)
    integer,        intent(in)    :: num

    integer :: i, j

    ! Copy node data from ori to copy with size num
    do i = 1, num
        copy(i).id      = ori(i).id         ! Node ID
        copy(i).bp      = ori(i).bp         ! Base pair ID
        copy(i).up      = ori(i).up         ! Upward ID
        copy(i).dn      = ori(i).dn         ! Downward ID
        copy(i).sec     = ori(i).sec        ! Cross-section ID
        copy(i).iniL    = ori(i).iniL       ! Initial line
        copy(i).croL    = ori(i).croL       ! Crossectional line
        copy(i).beveled = ori(i).beveled    ! Beveled node
        copy(i).conn    = ori(i).conn       ! Connection type
        copy(i).ghost   = ori(i).ghost      ! Ghost node

        ! Position and orientation vector
        do j = 1, 3
            copy(i).pos(j)   = ori(i).pos(j)
            copy(i).ori(j,:) = ori(i).ori(j,:)
        end do
    end do
end subroutine

! ---------------------------------------------------------------------------

! Copy element data from ori to copy with size num
subroutine Mani_Copy_EleType(ori, copy, num)
    type(EleType), intent(in)    :: ori(:)
    type(EleType), intent(inout) :: copy(:)
    integer,       intent(in)    :: num

    integer :: i

    ! Copy connectivity connectivity
    do i = 1, num
        copy(i).cn(:) = ori(i).cn(:)
    end do
end subroutine

! ---------------------------------------------------------------------------

! Initialize base type
subroutine Mani_Init_BaseType(base, n_base)
    type(BaseType), intent(inout) :: base(:)
    integer,        intent(in)    :: n_base

    integer :: i

    ! Initialize base data
    do i = 1, n_base
        base(i).id       = -1       ! Base ID
        base(i).node     = -1       ! Heritage node ID
        base(i).up       = -1       ! Upward base ID
        base(i).dn       = -1       ! Downward base ID
        base(i).xover    = -1       ! Possible crossover ID
        base(i).across   = -1       ! Complementary base ID
        base(i).strand   = -1       ! Strand ID
        base(i).pos(1:3) = 0.0d0    ! Position vector
    end do
end subroutine Mani_Init_BaseType

! ---------------------------------------------------------------------------

! Initialize strand type data
subroutine Mani_Init_StrandType(strand, n_strand)
    type(StrandType), allocatable, intent(inout) :: strand(:)
    integer, intent(in) :: n_strand

    integer :: i

    do i = 1, n_strand
        strand(i).n_base     = 0
        strand(i).b_circular = .false.
        strand(i).type1      = "NNNN"
        strand(i).type2      = "edge"
    end do
end subroutine Mani_Init_StrandType

! ---------------------------------------------------------------------------

! Copy base type
subroutine Mani_Copy_BaseType(ori, copy, num)
    type(BaseType), intent(in)    :: ori(:)
    type(BaseType), intent(inout) :: copy(:)
    integer,        intent(in)    :: num

    integer :: i

    ! Copy node data from ori to copy with size num
    do i = 1, num
        copy(i).id     = ori(i).id      ! Base ID
        copy(i).node   = ori(i).node    ! Heritage node ID
        copy(i).up     = ori(i).up      ! Upward base ID
        copy(i).dn     = ori(i).dn      ! Downward base ID
        copy(i).xover  = ori(i).xover   ! Possible crossover ID
        copy(i).across = ori(i).across  ! Complementary base ID
        copy(i).strand = ori(i).strand  ! Strand ID
        copy(i).pos(:) = ori(i).pos(:)  ! Position vector
    end do
end subroutine Mani_Copy_BaseType

! ---------------------------------------------------------------------------

! Go to starting position (down = -1)
function Mani_Go_Start_Base(dna, strand) result(base)
    type(DNAType), intent(in) :: dna
    integer,       intent(in) :: strand

    integer :: i, base

    base = dna.strand(strand).base(1)
    do i = 1, dna.strand(strand).n_base
        if(dna.top(base).dn == -1) exit
        base = dna.top(base).dn
    end do
end function Mani_Go_Start_Base

! ---------------------------------------------------------------------------

! Move to base far away n from start position
function Mani_Go_Base(dna, strand, n) result(base)
    type(DNAType), intent(in) :: dna
    integer,       intent(in) :: strand
    integer,       intent(in) :: n

    integer :: i, base

    base = Mani_Go_Start_Base(dna, strand)

    do i = 1, n
        base = dna.top(base).up
    end do
end function Mani_Go_Base

! ---------------------------------------------------------------------------

end module Mani