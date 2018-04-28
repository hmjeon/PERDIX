!
! =============================================================================
!
! Module - Data_Prob
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
module Data_Prob

! -----------------------------------------------------------------------------

    ! Problem type data structure
    type :: ProbType
        integer :: sel_prob             ! Number for pre-defined problem
        integer :: sel_vertex           ! Mitered or non-mitered
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
        character(200) :: path_work             ! Working directory path
    end type ProbType

! -----------------------------------------------------------------------------

end module Data_Prob