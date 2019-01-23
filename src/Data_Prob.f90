!
! =============================================================================
!
! Module - Data_Prob
! Last Updated : 01/09/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX, which allows scientists to build and solve
! the sequence design of complex DNA nanostructures.
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
module Data_Prob

! -----------------------------------------------------------------------------

    ! Problem type data structure
    type :: ProbType
        integer :: sel_prob             ! Number for pre-defined problem
        integer :: sel_vertex           ! Mitered or non-mitered
        integer :: sel_edge_sec         ! Edge section
        integer :: sel_edge_len         ! Edge length
        integer :: sel_edge_ref = 0     ! Reference edge
        integer :: n_edge_len           ! The number of bps each edge

        double precision :: p_mesh      ! Mesing parameter

        character(30000) :: scaf_seq    ! User-defined scaffold sequence

        integer :: color(3)             ! Problem color, [52, 152, 219]
        integer :: n_cng_min_stap = 0   ! The number of changing parameter for minimum staple length
        integer :: n_cng_max_stap = 0   ! The number of changing parameter for maximum staple length

        character(200) :: name_file             ! File name
        character(200) :: name_prob             ! Problem name
        character(10)  :: type_file             ! File type
        character(10)  :: type_geo = "closed"   ! Geometric type, open or closed
        character(200) :: path_work             ! Path for the working directory
        character(200) :: path_input            ! Path for the input geometry
    end type ProbType

! -----------------------------------------------------------------------------

end module Data_Prob