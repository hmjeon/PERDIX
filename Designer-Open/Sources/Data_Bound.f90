!
! ---------------------------------------------------------------------------------------
!
!                                 Module for Data_Bound
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/04/21
!                                                            Last  modified  : 2016/07/14
!
! ---------------------------------------------------------------------------------------
!
module Data_Bound

! ---------------------------------------------------------------------------------------

    ! Outer data structure
    type :: OuterType
        integer :: typeP    ! Point type, -1:internal, 1:outer
        integer :: n_neiP   ! The number of neighbor points
        integer :: n_newP   ! The number of new points

        integer, allocatable, dimension(:)   :: neiP    ! Neighbor point connectivity sorted
        integer, allocatable, dimension(:,:) :: newP    ! New point connectivity(index, heritage)
    end type OuterType

! ---------------------------------------------------------------------------------------

    ! Junction data structure
    type :: JuncType
        integer :: n_arm                ! The number of arms
        integer :: poi_c                ! Center position

        double precision :: ref_ang     ! Reference angle between two neighboring edges
        double precision :: tot_ang     ! Total angle at the junction
        double precision :: gap         ! Gap distance between junction and end of edges

        integer, allocatable, dimension(:)   :: iniL        ! Initial line
        integer, allocatable, dimension(:)   :: modP        ! Modified point
        integer, allocatable, dimension(:,:) :: croP        ! Sectional point (# of arms, # of sections)
        integer, allocatable, dimension(:,:) :: node        ! Nodes (# of arms, # of sections)
        integer, allocatable, dimension(:,:) :: conn        ! Node connectivity (node #, node # to be connected))
        integer, allocatable, dimension(:)   :: type_conn   ! Section connection type, negihbor(1) and self(2)
    end type JuncType

! ---------------------------------------------------------------------------------------

    ! Boundary data structure
    type :: BoundType
        integer :: n_outer              ! The number of outers
        integer :: n_junc               ! The number of juncs
        integer :: max_gap, min_gap     ! Maximum and minimum gap from vertex

        type(OuterType), allocatable, dimension(:) :: outer
        type(JuncType),  allocatable, dimension(:) :: junc
    end type BoundType

! ---------------------------------------------------------------------------------------

end module Data_Bound