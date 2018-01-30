!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Data_Bound
!
!                                                                    Updated : 2017/03/27
!
! Comments: Data structure of the junction
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Data_Bound

! ---------------------------------------------------------------------------------------

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

! ---------------------------------------------------------------------------------------

    ! Boundary data structure
    type :: BoundType
        integer :: n_outer      ! The number of outers
        integer :: n_junc       ! The number of juncs

        type(JuncType), allocatable :: junc(:)
    end type BoundType

! ---------------------------------------------------------------------------------------

end module Data_Bound