!
! ---------------------------------------------------------------------------------------
!
!                                  Module for Data_Mesh
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/04/21
!                                                            Last  modified  : 2016/07/28
!
! ---------------------------------------------------------------------------------------
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

        type(NodeType), allocatable, dimension(:) :: node   ! Node array
        type(EleType),  allocatable, dimension(:) :: ele    ! Element array
    end type

! ---------------------------------------------------------------------------------------

end module Data_Mesh