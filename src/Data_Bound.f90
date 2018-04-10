!
! =============================================================================
!
! Module - Data_Bound
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
module Data_Bound

! -----------------------------------------------------------------------------

    ! Junction data structure
    type :: JuncType
        integer :: n_arm                ! The number of arms
        integer :: poi_c                ! Center position

        integer :: n_un_scaf = 0        ! # of unpaired nucleotide in scaf
        integer :: n_un_stap = 0        ! # of unpaired nucleotide in stap

        double precision :: ref_ang     ! Reference angle between two neighboring edges
        double precision :: tot_ang     ! Total angle at the junction
        double precision :: gap         ! Gap distance between junction and end of edges

        integer, allocatable :: iniL(:)         ! Initial line
        integer, allocatable :: modP(:)         ! Modified point
        integer, allocatable :: croP(:,:)       ! Sectional point (# of arms, # of sections)
        integer, allocatable :: node(:,:)       ! Nodes (# of arms, # of sections)
        integer, allocatable :: conn(:,:)       ! Node connectivity (node #, node # to be connected))
        integer, allocatable :: type_conn(:)    ! Section connection type, negihbor(1) and self(2)
    end type JuncType

! -----------------------------------------------------------------------------

    ! Boundary data structure
    type :: BoundType
        integer :: n_outer      ! The number of outers
        integer :: n_junc       ! The number of juncs

        type(JuncType), allocatable :: junc(:)
    end type BoundType

! -----------------------------------------------------------------------------

end module Data_Bound