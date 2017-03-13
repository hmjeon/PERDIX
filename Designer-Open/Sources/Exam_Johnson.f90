!
! ---------------------------------------------------------------------------------------
!
!                                Module for Exam_Johnson
!
!                                                             First modified : 2017/03/13
!                                                             Last  modified : 2017/03/13
!
! Comments: This module contains the geometric definition of 2D open structures
!
! by Hyungmin Jun (Hyungminjun@outlook.com), MIT, Bathe Lab, 2017
!
! Copyright 2017. Massachusetts Institute of Technology. Rights Reserved.
! M.I.T. hereby makes following copyrightable material available to the
! public under GNU General Public License, version 2 (GPL-2.0). A copy of
! this license is available at https://opensource.org/licenses/GPL-2.0
!
! ---------------------------------------------------------------------------------------
!
! Reference: http://www.georgehart.com/virtual-polyhedra/johnson-index.html
! Reference: https://en.wikipedia.org/wiki/List_of_Johnson_solids
!
module Exam_Johnson

    use Data_Prob
    use Data_Geom

    use Mani

    implicit none

    public Exam_Johnson_Square_Pyramid_J1                       ! V= 5, E= 8, F= 5
    public Exam_Johnson_Pentagonal_Pyramid_J2                   ! V= 6, E=10, F= 6
    public Exam_Johnson_Triangular_Cupola_J3                    ! V= 9, E=15, F= 8
    public Exam_Johnson_Square_Cupola_J4                        ! V=12, E=20, F=10
    public Exam_Johnson_Pentagonal_Cupola_J5                    ! V=15, E=25, F=12
    public Exam_Johnson_Pentagonal_Rotunda_J6                   ! V=20, E=35, F=17
    public Exam_Johnson_Elongated_Triangular_Cupola_J18         ! V=15, E=27, F=14
    public Exam_Johnson_Elongated_Square_Cupola_J19             ! V=20, E=36, F=18
    public Exam_Johnson_Elongated_Pentagonal_Cupola_J20         ! V=25, E=45, F=22
    public Exam_Johnson_Elongated_Pentagonal_Rotunda_J21        ! V=30, E=55, F=27
    public Exam_Johnson_Gyroelongated_Triangular_Cupola_J22     ! V=15, E=33, F=20
    public Exam_Johnson_Gyroelongated_Square_Cupola_J23         ! V=20, E=44, F=26
    public Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24     ! V=25, E=55, F=32
    public Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25    ! V=30, E=65, F=37
    public Exam_Johnson_Gyrobifastigium_J26                     ! V= 8, E=14, F= 8

    contains

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Square_Pyramid_J1(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Square_Pyramid_J1

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Pyramid_J2(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Pentagonal_Pyramid_J2

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Triangular_Cupola_J3(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Triangular_Cupola_J3

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Square_Cupola_J4(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Square_Cupola_J4

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Cupola_J5(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Pentagonal_Cupola_J5

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Pentagonal_Rotunda_J6(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Pentagonal_Rotunda_J6

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Triangular_Cupola_J18(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Elongated_Triangular_Cupola_J18

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Square_Cupola_J19(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Elongated_Square_Cupola_J19

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Cupola_J20(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Elongated_Pentagonal_Cupola_J20

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Elongated_Pentagonal_Rotunda_J21(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Elongated_Pentagonal_Rotunda_J21

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Triangular_Cupola_J22(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Gyroelongated_Triangular_Cupola_J22

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Square_Cupola_J23(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Gyroelongated_Square_Cupola_J23

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Gyroelongated_Pentagonal_Cupola_J24

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Gyroelongated_Pentagonal_Rotunda_J25

! ---------------------------------------------------------------------------------------

! Example of square pyramid - J1
! Last updated on Mon 13 Mar 2017 by Hyungmin
subroutine Exam_Johnson_Gyrobifastigium_J26(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "test"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "test"

    ! Set geometric type and view
end subroutine Exam_Johnson_Gyrobifastigium_J26

! ---------------------------------------------------------------------------------------

end module Exam_Johnson