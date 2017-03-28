!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Data_Geom
!
!                                                                    Updated : 2017/03/27
!
! Comments: The data structure of the geometry and cross-section
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
module Data_Geom

! ---------------------------------------------------------------------------------------

    ! Section data for arbitrary cross-section
    type :: SecType
        character(len=10) :: types                  ! Lattice type, squre or honeycomb

        integer :: dir                              ! Section convection of caDNAnano
        integer :: maxR, minR                       ! Maximum and minimum row
        integer :: maxC, minC                       ! Maximum and minimum column
        integer :: n_row, n_col                     ! Size of row and column
        integer :: ref_row                          ! Reference row to set t-axis
        integer :: ref_maxC, ref_minC               ! Maximum and minimum column in reference row

        integer, allocatable, dimension(:) :: id    ! Section ID
        integer, allocatable, dimension(:) :: posR  ! Row position number
        integer, allocatable, dimension(:) :: posC  ! Column position number

        ! Connetivity for self connection route
        ! -1 : neighbor connection, num : section number to be connected for self-connection
        integer, allocatable, dimension(:) :: conn
    end type SecType

! ---------------------------------------------------------------------------------------

    ! Point type data structure
    type :: PointType
        double precision :: pos(3)      ! Position vector
    end type PointType

! ---------------------------------------------------------------------------------------

    ! Line type data structure
    type :: LineType
        integer :: iniL             ! Inital line number
        integer :: sec              ! Cross-section ID
        integer :: iniP(2)          ! Initial point connectivity
        integer :: poi(2)           ! Point connectivity
        integer :: neiF(2)          ! Neighbor faces  (direction)
        integer :: neiP(2, 2)       ! Neighbor points (point number, direction)
        integer :: neiL(2, 2)       ! Neighbor lines  (point number, direction)
                                    ! Point number : 1 - starting, 2 - ending
                                    ! Direction    : 1 - left, 2 - right
        integer :: n_xover          ! The number of scaffold crossovers
        double precision :: t(3, 3) ! Local vector at the center
    end type LineType

! ---------------------------------------------------------------------------------------

    ! Face type data structure
    type :: FaceType
        integer :: n_poi                            ! The number of points
        integer, allocatable, dimension (:) :: poi  ! Connectivity
    end type FaceType

! ---------------------------------------------------------------------------------------

    ! Geometry data type to manage section, point, line and face data
    type :: Geomtype
        integer :: n_sec                    ! The number of sections
        integer :: n_iniP, n_modP, n_croP   ! The number of initial, modified and sectional points
        integer :: n_iniL, n_croL           ! The number of initial, sectional lines
        integer :: n_face                   ! The number of faces
        integer :: min_edge_length
        integer :: max_edge_length

        type(SecType) :: sec                                            ! section
        type(PointType), allocatable, dimension(:) :: iniP, modP, croP  ! points
        type(LineType),  allocatable, dimension(:) :: iniL, croL        ! lines
        type(FaceType),  allocatable, dimension(:) :: face              ! face
    end type Geomtype

! ---------------------------------------------------------------------------------------

end module Data_Geom