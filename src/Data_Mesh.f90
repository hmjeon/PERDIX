!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Data_Mesh
!
!                                                                    Updated : 2017/03/27
!
! Comments: The data structure of the basepair model
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2018. All rights reserved.
!
module Data_Mesh

! ---------------------------------------------------------------------------------------

    ! Node(base pair) data type structure
    type :: NodeType
        integer :: id       ! Node ID
        integer :: bp       ! Base pair ID
        integer :: up       ! Upward ID
        integer :: dn       ! Downward ID
        integer :: sec      ! Section ID
        integer :: iniL     ! Initial line ID
        integer :: croL     ! Cross-section line ID
        integer :: beveled  ! Beveled node, 1: beveled node

        ! Nodal connectivity at the section
        ! -1 - no-connection, 1 - neighbor, 2 - self, 3 - modified neighbor, 4 - modified self
        integer :: conn = -1

        ! Ghost node : -1 : normal node, 1 : ghost node to be deleted
        integer :: ghost

        double precision :: pos(3)      ! Position vector
        double precision :: ori(3, 3)   ! Orientation vector
    end type NodeType

    ! ElementType structure
    type :: EleType
        integer :: cn(2)    ! Connectivity
    end type EleType

! ---------------------------------------------------------------------------------------

    ! MeshType data structure
    type :: MeshType
        integer :: n_node = 0   ! The number of nodes
        integer :: n_ele  = 0   ! The number of elements
        integer :: n_beveled    ! The number of beveled nodes

        type(NodeType), allocatable :: node(:)   ! Node array
        type(EleType),  allocatable :: ele(:)    ! Element array
    end type

! ---------------------------------------------------------------------------------------

end module Data_Mesh