!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Data_Prob
!
!                                                                    Updated : 2017/03/27
!
! Comments: The data structure of the problem definition
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2018. All rights reserved.
!
module Data_Prob

! ---------------------------------------------------------------------------------------

    ! Problem type data structure
    type :: ProbType
        integer :: sel_prob             ! Number for pre-defined problem
        integer :: sel_vertex           ! Beveled or non-beveled
        integer :: sel_sec              ! Number for pre-defined cross-section
        integer :: sel_bp_edge          ! Number for pre-defined # of base pairs on edges
        integer :: sel_edge             ! Specific edge number to set as reference
        integer :: n_bp_edge            ! The number of bps each edge

        integer :: color(3)             ! Problem color, [52, 152, 219]
        integer :: n_cng_min_stap = 0   ! The number of changing parameter for minimum staple length
        integer :: n_cng_max_stap = 0   ! The number of changing parameter for maximum staple length
        double precision :: scale       ! Problem scale for post-processing(atomic model)
        double precision :: size        ! Problem size for post-processing(cylindrical model)
        double precision :: move_x      ! To adjust center position
        double precision :: move_y      ! To adjust conter position

        character(200) :: name_file             ! File name
        character(200) :: name_prob             ! Problem name
        character(10)  :: type_file             ! File type
        character(10)  :: type_geo = "closed"   ! Geometric type, open or closed
        character(200) :: path_work1            ! Working directory path
        character(200) :: path_work2            ! Working directory path
        character(200) :: path_chimera          ! Chimera path
    end type ProbType

! ---------------------------------------------------------------------------------------

end module Data_Prob