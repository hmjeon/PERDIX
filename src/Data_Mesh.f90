!
! =============================================================================
!
! Module - Data_Mesh
! Last Updated : 04/10/2018, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX, which allows scientists to build and solve
! the sequence design of complex DNAnanostructures.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! PERDIX is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! PERDIX is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Data_Mesh

! -----------------------------------------------------------------------------

    ! Node(base pair) data type structure
    type :: NodeType
        integer :: id       ! Node ID
        integer :: bp       ! Base pair ID
        integer :: up       ! Upward ID
        integer :: dn       ! Downward ID
        integer :: sec      ! Section ID
        integer :: iniL     ! Initial line ID
        integer :: croL     ! Cross-section line ID
        integer :: mitered  ! Mitered node, 1: mitered node

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

! -----------------------------------------------------------------------------

    ! MeshType data structure
    type :: MeshType
        integer :: n_node = 0   ! The number of nodes
        integer :: n_ele  = 0   ! The number of elements
        integer :: n_mitered    ! The number of mitered nodes

        type(NodeType), allocatable :: node(:)   ! Node array
        type(EleType),  allocatable :: ele(:)    ! Element array
    end type

! -----------------------------------------------------------------------------

end module Data_Mesh