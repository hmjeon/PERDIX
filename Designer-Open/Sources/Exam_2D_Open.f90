!
! ---------------------------------------------------------------------------------------
!
!                               Module for Exam_2D Open
!
!                                                             First modified : 2017/03/10
!                                                             Last  modified : 2017/03/10
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
module Exam_2D_Open

    use Data_Prob
    use Data_Geom

    use Para
    use Mani
    use Math

    implicit none

    public Exam_Open2D_Plate_Uniform_Quad           !  1. Plate with uniform quad mesh
    public Exam_Open2D_Plate_Uniform_Tri            !  2. Plate with uniform tri mesh
    public Exam_Open2D_Plate_Distorted_Quad         !  3. Plate with distorted quad mesh
    public Exam_Open2D_Plate_Distorted_Tri          !  4. Plate with distorted tri mesh
    public Exam_Open2D_Circular_Plate_Quad          !  5. Circular plate with quad mesh
    public Exam_Open2D_Circular_Plate_Tri           !  6. Circular plate with tri mesh
    public Exam_Open2D_Annular_Plate_Quad           !  7. Annular plate with quad mesh
    public Exam_Open2D_Annular_Plate_Tri            !  8. Annular plate with tri mesh
    public Exam_Open2D_Hyperbolic_Paraboloid_Quad   !  9. Hyperbolic paraboloid with quad mesh
    public Exam_Open2D_Hyperbolic_Paraboloid_Tri    ! 10. Hyperbolic paraboloid with tri mesh
    public Exam_Open2D_Circle_Tri_Coarse            ! 11. Circle with tri coarse mesh
    public Exam_Open2D_Circle_Tri_Fine              ! 12. Circle with tri fine mesh
    public Exam_Open2D_Ellipse_Tri_Coarse           ! 13. Ellipse with tri coarse mesh
    public Exam_Open2D_Ellipse_Tri_Fine             ! 14. Ellipse with tri fine mesh

    private Exam_Open2D_Cross_Point
    private Exam_Open2D_Merge_Point_Face_Quad
    private Exam_Open2D_Merge_Point_Face_Tri
    private Exam_Open2D_Comp_XYZ

contains

! ---------------------------------------------------------------------------------------

! Example of plate with uniform quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Plate_Uniform_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: i, j, index, n_i_poi, n_j_poi, n, nx, ny
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "1_Plate_Uniform_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Uniform Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    ! Set options
    n  = 3
    nx = n
    ny = n

    n_i_poi = nx + 1
    n_j_poi = ny + 1
    x_width = dble(nx)
    y_width = dble(ny)
    del_x   = x_width / dble(n_i_poi - 1)
    del_y   = y_width / dble(n_j_poi - 1)

    geom.n_iniP = n_i_poi * n_j_poi
    geom.n_face = nx * ny

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    do j = 1, n_j_poi
        do i = 1, n_i_poi
            index = n_i_poi * (j - 1) + i
            geom.iniP(index).pos(1) = del_x * dble(i - 1)
            geom.iniP(index).pos(2) = del_y * dble(j - 1)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, ny
        do i = 1, nx
            index = nx * (j - 1) + i
            geom.face(index).poi(1) = n_i_poi * (j - 1) + i
            geom.face(index).poi(2) = n_i_poi * (j - 1) + i + 1
            geom.face(index).poi(3) = n_i_poi * j + i + 1
            geom.face(index).poi(4) = n_i_poi * j + i
        end do
    end do
end subroutine Exam_Open2D_Plate_Uniform_Quad

! ---------------------------------------------------------------------------------------

! Example of plate with uniform tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Plate_Uniform_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: i, j, index, n_i_poi, n_j_poi, n, nx, ny
    character :: pn
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "02_Plate_Uniform_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Uniform Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    ! Set options
    n  = 2
    nx = n * 2
    ny = n
    pn = "\"

    n_i_poi = nx + 1
    n_j_poi = ny + 1
    x_width = dble(nx)
    y_width = dble(ny)
    del_x   = x_width / dble(n_i_poi - 1)
    del_y   = y_width / dble(n_j_poi - 1)

    geom.n_iniP = n_i_poi * n_j_poi
    geom.n_face = nx  * ny * 2
    
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    do j = 1, n_j_poi
        do i = 1, n_i_poi
            index = n_i_poi * (j - 1) + i
            geom.iniP(index).pos(1) = del_x * dble(i - 1)
            geom.iniP(index).pos(2) = del_y * dble(j - 1)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, ny
        do i = 1, nx
            if(pn == "\") then          ! Mesh pattern, \
                index = 2*(nx * (j-1) + i) - 1
                geom.face(index).poi(1) = n_i_poi * (j-1) + i
                geom.face(index).poi(2) = n_i_poi * (j-1) + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i

                index = 2*(nx * (j-1) + i)
                geom.face(index).poi(1) = n_i_poi * (j-1) + i + 1
                geom.face(index).poi(2) = n_i_poi * j + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i
            else if(pn == "/") then     ! Mesh pattern, /
                index = 2*(nx * (j-1) + i) - 1
                geom.face(index).poi(1) = n_i_poi * (j - 1) + i
                geom.face(index).poi(2) = n_i_poi * j + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i

                index = 2*(nx * (j-1) + i)
                geom.face(index).poi(1) = n_i_poi * (j - 1) + i
                geom.face(index).poi(2) = n_i_poi * (j - 1) + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i + 1
            end if
        end do
    end do
end subroutine Exam_Open2D_Plate_Uniform_Tri

! ---------------------------------------------------------------------------------------

! Example of Plate with distorted quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Plate_Distorted_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision, allocatable :: edge(:,:,:)
    double precision :: x_width, y_width, del_x, del_y, ff, dd, xy(2)
    integer :: i, j, dt, index, n_i_poi, n_j_poi, n, nx, ny
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "3_Plate_Distorted_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    ! Set options
    n  = 3
    dt = 1

    nx      = n
    ny      = n
    n_i_poi = nx + 1
    n_j_poi = ny + 1
    x_width = 1.0d0
    y_width = 1.0d0
    del_x   = x_width / dble(n_i_poi - 1)
    del_y   = y_width / dble(n_j_poi - 1)

    geom.n_iniP = n_i_poi * n_j_poi
    geom.n_face = nx  * ny

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    allocate(edge(4, n + 1, 2))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Routines for distorted mesh & calculate nodes along adges
    ff = 0.0d0
    do i = 1, n
        ff = ff + i
    end do

    dd = 1.0d0 / ff
    do i = 1, n + 1
        ff = 0.0d0
        if(i > 1) then
            do j = 1, i - 1
                ff = ff + j
            end do
        end if
        edge(1, i, 1) = ff * dd               ! Bottom
        edge(1, i, 2) = 0.0d0                 ! Bottom
        edge(2, n+2-i, 1) = 1.0d0 - ff * dd   ! Top
        edge(2, i, 2) = 1.0d0                 ! Top
        edge(3, i, 1) = 0.0d0                 ! Left
        edge(3, i, 2) = ff * dd               ! Left
        edge(4, i, 1) = 1.0d0                 ! Right
        edge(4, n+2-i, 2) = 1.0d0 - ff * dd   ! Right
    end do

    ! Set position vector
    do j = 1, n_j_poi
        do i = 1, n_i_poi
            index = n_i_poi * (j - 1) + i

            call Exam_Open2D_Cross_Point(edge(1,i,:), edge(2,i,:), edge(3,j,:), edge(4,j,:), xy)

            geom.iniP(index).pos(1) = xy(1)
            geom.iniP(index).pos(2) = xy(2)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, ny
        do i = 1, nx
            index = nx * (j - 1) + i
            geom.face(index).poi(1) = n_i_poi * (j - 1) + i
            geom.face(index).poi(2) = n_i_poi * (j - 1) + i + 1
            geom.face(index).poi(3) = n_i_poi * j + i + 1
            geom.face(index).poi(4) = n_i_poi * j + i
        end do
    end do

    ! Deallocate edge data
    deallocate(edge)
end subroutine Exam_Open2D_Plate_Distorted_Quad

! ---------------------------------------------------------------------------------------

! Example of Plate with distorted tri mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Plate_Distorted_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision, allocatable :: edge(:,:,:)
    double precision :: x_width, y_width, del_x, del_y, ff, dd, xy(2)
    integer :: i, j, n, dt, index, n_i_poi, n_j_poi, nx, ny
    character :: pn
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "4_Plate_Distorted_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    n  = 3
    pn = "\"
    dt = 1

    nx      = n
    ny      = n
    n_i_poi = nx + 1
    n_j_poi = ny + 1
    x_width = 1.0d0
    y_width = 1.0d0
    del_x   = x_width / dble(n_i_poi - 1)
    del_y   = y_width / dble(n_j_poi - 1)

    geom.n_iniP = n_i_poi * n_j_poi
    geom.n_face = nx  * ny * 2

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    allocate(edge(4, n + 1, 2))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! routines for distorted mesh & calculate nodes along adges
    ff = 0.0d0
    do i = 1, n
        ff = ff + i
    end do

    dd = 1.0d0 / ff
    do i = 1, n + 1
        ff = 0.0d0
        if(i > 1) then
            do j = 1, i - 1
                ff = ff + j
            end do
        end if
        edge(1,i,1)     = ff * dd               ! bottom
        edge(1,i,2)     = 0.0d0                 ! bottom
        edge(2,n+2-i,1) = 1.0d0 - ff * dd       ! top
        edge(2,i,2)     = 1.0d0                 ! top
        edge(3,i,1)     = 0.0d0                 ! left
        edge(3,i,2)     = ff * dd               ! left
        edge(4,i,1)     = 1.0d0                 ! right
        edge(4,n+2-i,2) = 1.0d0 - ff * dd       ! right
    end do

    ! Set position vector
    do j = 1, n_j_poi
        do i = 1, n_i_poi
            index = n_i_poi * (j - 1) + i

            call Exam_Open2D_Cross_Point(edge(1, i, :), edge(2, i, :), edge(3, j, :), edge(4, j, :), xy)

            geom.iniP(index).pos(1) = xy(1)
            geom.iniP(index).pos(2) = xy(2)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, ny
        do i = 1, nx
            if(pn == "\") then          ! mesh pattern, \
                index = 2*(nx * (j-1) + i) - 1
                geom.face(index).poi(1) = n_i_poi * (j-1) + i
                geom.face(index).poi(2) = n_i_poi * (j-1) + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i

                index = 2*(nx * (j-1) + i)
                geom.face(index).poi(1) = n_i_poi * (j-1) + i + 1
                geom.face(index).poi(2) = n_i_poi * j + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i
            else if(pn == "/") then     ! mesh pattern, /
                index = 2*(nx * (j-1) + i) - 1
                geom.face(index).poi(1) = n_i_poi * (j - 1) + i
                geom.face(index).poi(2) = n_i_poi * j + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i

                index = 2*(nx * (j-1) + i)
                geom.face(index).poi(1) = n_i_poi * (j - 1) + i
                geom.face(index).poi(2) = n_i_poi * (j - 1) + i + 1
                geom.face(index).poi(3) = n_i_poi * j + i + 1
            end if
        end do
    end do

    ! deallocate edge data
    deallocate(edge)
end subroutine Exam_Open2D_Plate_Distorted_Tri

! ---------------------------------------------------------------------------------------

! Example of circular plate with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Circular_Plate_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision,allocatable  :: edge(:,:,:), joint(:,:)
    integer, allocatable :: conn(:,:)

    double precision :: dx, dy
    integer :: i, j, n, count_n, count_e, count_n_t, count_e_t
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "5_Circular_Plate_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Circular Plate Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    n           = 2
    geom.n_iniP = (n + 1)*(n + 1 ) + n*(n + 1) + n*n
    geom.n_face = 3*n*n

    ! allocate data structures
    allocate(edge(2, n+1, 2), conn(n*n, 4), joint((n+1)**2, 3))
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    
    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! First mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1, i, 1) = dble(i-1) * 0.5d0 / dble(n)             ! x1
        edge(1, i, 2) = 0.0d0                                   ! y1
        ! Second edge
        edge(2, i, 1) = dble(i-1) * 1.0d0 / 3.0d0 / dble(n)     ! x2
        edge(2, i, 2) = 0.5d0 - 0.5d0 * edge(2, i, 1)           ! y2
    end do

    ! Generation of joint
    count_n = 0
    do i = 1, n + 1
        dx = (edge(2,i,1) - edge(1,i,1)) / dble(n)
        dy = (edge(2,i,2) - edge(1,i,2)) / dble(n)

        do j = 1, n + 1
            count_n = count_n + 1
            geom.iniP(count_n).pos(1) = edge(1,i,1) + dble(j-1) * dx
            geom.iniP(count_n).pos(2) = edge(1,i,2) + dble(j-1) * dy
            geom.iniP(count_n).pos(3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e = 0
    do i = 1, n
        do j = 1, n
            count_e = count_e + 1
            geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 1
            geom.face(count_e).poi(2) = (n+1)*(i+0) + (j-1) + 2
            geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 2
            geom.face(count_e).poi(4) = (n+1)*(i-1) + (j-1) + 1
        end do
    end do

    ! Second mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1, n+2-i, 1) = dble(i-1) * 1.0d0 / 3.0d0 / dble(n)     ! x1
        edge(1, n+2-i, 2) = 0.5d0 - 0.5d0 * edge(1, n+2-i, 1)       ! y1
        ! Second edge
        edge(2, n+2-i, 1) = dsin(dble(i-1)*(pi/180.0d0) * 45.d0 / dble(n))    ! x2
        edge(2, n+2-i, 2) = dcos(dble(i-1)*(pi/180.0d0) * 45.d0 / dble(n))    ! y2
    end do

    ! Generation of joint
    count_n_t = 0
    do i = 1, n + 1
        dx = (edge(2, i, 1) - edge(1, i, 1)) / dble(n)
        dy = (edge(2, i, 2) - edge(1, i, 2)) / dble(n)

        do j = 1, n + 1
            count_n_t = count_n_t + 1
            joint(count_n_t, 1) = edge(1,i,1) + dble(j-1) * dx
            joint(count_n_t, 2) = edge(1,i,2) + dble(j-1) * dy
            joint(count_n_t, 3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e_t = 0
    do i = 1, n
        do j = 1, n
            count_e_t = count_e_t + 1

            conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 1
            conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
            conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 2
            conn(count_e_t, 4) = (n+1)*(i+0) + (j-1) + 1
        end do
    end do

    call Exam_Open2D_Merge_Point_Face_Quad(n, joint, conn, count_n, count_e, geom)

    ! Third mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1,i,1) = 1.0d0/3.d0 + dble(i-1) * (0.5d0 - 1.0d0 / 3.0d0) / dble(n)        ! x1
        edge(1,i,2) = 1.0d0 - 2.0d0 * edge(1,i,1)                                       ! y1
        ! Second edge
        edge(2,i,1) = dsin(pi/180.0d0 * 45.d0 + dble(i-1) * (pi/180.0d0) * 45.d0 / dble(n))   ! x2
        edge(2,i,2) = dcos(pi/180.0d0 * 45.d0 + dble(i-1) * (pi/180.0d0) * 45.d0 / dble(n))   ! y2
    end do

    ! Generation of joint
    count_n_t = 0
    do i = 1, n + 1
        dx = (edge(2, i, 1) - edge(1, i, 1)) / dble(n)
        dy = (edge(2, i, 2) - edge(1, i, 2)) / dble(n)

        do j = 1, n + 1
            count_n_t = count_n_t + 1
            joint(count_n_t, 1) = edge(1, i, 1) + dble(j-1) * dx
            joint(count_n_t, 2) = edge(1, i, 2) + dble(j-1) * dy
            joint(count_n_t, 3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e_t = 0
    do i = 1, n
        do j = 1, n
            count_e_t = count_e_t + 1
            conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
            conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
            conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 2
            conn(count_e_t, 4) = (n+1)*(i-1) + (j-1) + 1
        end do
    end do

    call Exam_Open2D_Merge_Point_Face_Quad(n, joint, conn, count_n, count_e, geom)

    ! Deallocate data structure
    deallocate(edge, joint, conn)
end subroutine Exam_Open2D_Circular_Plate_Quad

! ---------------------------------------------------------------------------------------

! Example of circular plate with tri mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Circular_Plate_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision,allocatable  :: edge(:,:,:), joint(:,:)
	integer, allocatable :: conn(:,:)

	double precision :: dx, dy
    integer :: i, j, n, count_n, count_e, count_n_t, count_e_t
    character :: pn
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "6_Circular_Plate_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Circular Plate Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    n  = 2
    pn = "\"

    geom.n_iniP = (n + 1)*(n + 1) + n*(n + 1) + n*n
    geom.n_face = 3*n*n*2

    ! Allocate data structures
    allocate(edge(2, n+1, 2), conn(n*n*2, 4), joint((n+1)**2, 3))
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! First mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1, i, 1) = dble(i-1) * 0.5d0 / dble(n)             ! x1
        edge(1, i, 2) = 0.0d0                                   ! y1
        ! Second edge
        edge(2, i, 1) = dble(i-1) * 1.0d0 / 3.0d0 / dble(n)     ! x2
        edge(2, i, 2) = 0.5d0 - 0.5d0 * edge(2, i, 1)           ! y2
    end do

    ! Generation of joint
    count_n = 0
    do i = 1, n + 1
        dx = (edge(2, i, 1) - edge(1, i, 1)) / dble(n)
        dy = (edge(2, i, 2) - edge(1, i, 2)) / dble(n)

        do j = 1, n + 1
            count_n = count_n + 1
            geom.iniP(count_n).pos(1) = edge(1, i, 1) + dble(j-1) * dx
            geom.iniP(count_n).pos(2) = edge(1, i, 2) + dble(j-1) * dy
            geom.iniP(count_n).pos(3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e = 0
    do i = 1, n
        do j = 1, n
            if(pn == "\") then
                count_e = count_e + 1
                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 1
                geom.face(count_e).poi(2) = (n+1)*(i+0) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 1

                count_e = count_e + 1
                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 2
                geom.face(count_e).poi(2) = (n+1)*(i-1) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 1
            else if(pn == "/") then
                count_e = count_e + 1
                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 1
                geom.face(count_e).poi(2) = (n+1)*(i-1) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 1

                count_e = count_e + 1
                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 1
                geom.face(count_e).poi(2) = (n+1)*(i+0) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 2
            end if
        end do
    end do

    ! Second mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1, n+2-i, 1) = dble(i-1) * 1.0d0 / 3.0d0 / dble(n)     ! x1
        edge(1, n+2-i, 2) = 0.5d0 - 0.5d0 * edge(1, n+2-i, 1)       ! y1
        ! Second edge
        edge(2, n+2-i, 1) = dsin(dble(i-1)*(pi/180.0d0) * 45.d0 / dble(n))    ! x2
        edge(2, n+2-i, 2) = dcos(dble(i-1)*(pi/180.0d0) * 45.d0 / dble(n))    ! y2
    end do

    ! Generation of joint
    count_n_t = 0
    do i = 1, n + 1
        dx = (edge(2, i, 1) - edge(1, i, 1)) / dble(n)
        dy = (edge(2, i, 2) - edge(1, i, 2)) / dble(n)

        do j = 1, n + 1
            count_n_t = count_n_t + 1
            joint(count_n_t, 1) = edge(1,i,1) + dble(j-1) * dx
            joint(count_n_t, 2) = edge(1,i,2) + dble(j-1) * dy
            joint(count_n_t, 3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e_t = 0
    do i = 1, n
        do j = 1, n
            if(pn == "\") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 1
            else if(pn == "/") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 2
            end if
        end do
    end do

    call Exam_Open2D_Merge_Point_Face_Tri(n, joint, conn, count_n, count_e, geom)

    ! Third mesh, generate positions of edges
    do i = 1, n + 1
        ! First edge
        edge(1,i,1) = 1.0d0/3.d0 + dble(i-1) * (0.5d0 - 1.0d0 / 3.0d0) / dble(n)        ! x1
        edge(1,i,2) = 1.0d0 - 2.0d0 * edge(1,i,1)                                       ! y1
        ! Second edge
        edge(2,i,1) = dsin(pi/180.0d0 * 45.d0 + dble(i-1) * (pi/180.0d0) * 45.d0 / dble(n))   ! x2
        edge(2,i,2) = dcos(pi/180.0d0 * 45.d0 + dble(i-1) * (pi/180.0d0) * 45.d0 / dble(n))   ! y2
    end do

    ! Generation of joint
    count_n_t = 0
    do i = 1, n + 1
        dx = (edge(2, i, 1) - edge(1, i, 1)) / dble(n)
        dy = (edge(2, i, 2) - edge(1, i, 2)) / dble(n)

        do j = 1, n + 1
            count_n_t = count_n_t + 1
            joint(count_n_t, 1) = edge(1, i, 1) + dble(j-1) * dx
            joint(count_n_t, 2) = edge(1, i, 2) + dble(j-1) * dy
            joint(count_n_t, 3) = 0.0d0
        end do
    end do

    ! Generation of connectivities for quadrilaterlar elements
    count_e_t = 0
    do i = 1, n
        do j = 1, n
            if(pn == "\") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 2
            else if(pn == "/") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 1
            end if
        end do
    end do

    call Exam_Open2D_Merge_Point_Face_Tri(n, joint, conn, count_n, count_e, geom)

    ! Deallocate data structure
    deallocate(edge, joint, conn)
end subroutine Exam_Open2D_Circular_Plate_Tri

! ---------------------------------------------------------------------------------------

! Example of Annular plate with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Annular_Plate_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: i_radius, o_radius, angle, radius, force
    integer :: i, j, count, n, nx, nr
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "7_Annular_Plate_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Annular Plate Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"

    n        = 2
    o_radius = 4.0d0        ! Outer radius
    i_radius = 2.4d0        ! Internal radius
    nx       = n
    nr       = 5*n

    geom.n_iniP = (nx + 1) * nr
    geom.n_face = nx * nr

    ! Allocate point structure and set position vector
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    count = 0
    do j = 1, nr
        do i = 1, nx + 1
            count  = count + 1
            angle  = 2.0d0 * pi / dble(nr) * dble(j - 1)
            radius = ((o_radius - i_radius) / dble(nx) * dble(i - 1)) + i_radius

            geom.iniP(count).pos(1) = radius * dcos(angle)
            geom.iniP(count).pos(2) = radius * dsin(angle)
            geom.iniP(count).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    count = 0
    do i = 1, nr
        do j = 1, nx
            if(i /= nr) then
                count = count + 1
                geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 2
                geom.face(count).poi(4) = (nx+1)*(i-0) + (j-1) + 1
            else if(i == nr) then
                count = count + 1
                geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 2
                geom.face(count).poi(4) = (nx+1)*(1-1) + (j-1) + 1
            end if
        end do
    end do
end subroutine Exam_Open2D_Annular_Plate_Quad

! ---------------------------------------------------------------------------------------

! Example of annular plate with triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Annular_Plate_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: i_radius, o_radius, angle, radius, force
    integer :: i, j, count, n, nx, nr
    character :: pattern
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "60_Annular_Plate_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Annular Plate Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n        = 3
    pattern  = "\"
    o_radius = 1.0d0        ! Outer radius
    i_radius = 0.6d0        ! Internal radius
    nx       = n
    nr       = 5 * n

    geom.n_iniP = (nx + 1)*nr
    geom.n_face = nx*nr*2

    ! Set position vector
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    count = 0
    do j = 1, nr
        do i = 1, nx + 1
            count  = count + 1
            angle  = 2.0d0 * pi / dble(nr) * dble(j - 1)
            radius = ((o_radius - i_radius) / dble(nx) * dble(i - 1)) + i_radius

            geom.iniP(count).pos(1) = radius * dcos(angle)
            geom.iniP(count).pos(2) = radius * dsin(angle)
            geom.iniP(count).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    count = 0
    do i = 1, nr
        do j = 1, nx
            if(i /= nr) then
                if(pattern == "\") then
                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 1

                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(2) = (nx+1)*(i-0) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 1
                else if(pattern == "\") then
                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(i-0) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 1

                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 2
                end if
            else if(i == nr) then
                if(pattern == "\") then
                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 1

                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(2) = (nx+1)*(1-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 1
                else if(pattern == "\") then

                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(1-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 1

                    count = count + 1
                    geom.face(count).n_poi = 3
                    allocate(geom.face(count).poi(3))
                    geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 2
                end if
            end if
        end do
    end do
end subroutine Exam_Open2D_Annular_Plate_Tri

! ---------------------------------------------------------------------------------------

! Example of hyperbolic paraboloid with quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Hyperbolic_Paraboloid_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, j, count, n, nx, ny
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "61_Hyperbolic_Paraboloid_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hyperbolic Paraboloid Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n  = 3
    nx = n
    ny = n

    geom.n_iniP = (nx + 1)*(ny + 1)
    geom.n_face = nx*ny

    ! Allocate point structure and set position vector
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    count = 0
    do i = 1, nx + 1
        do j = 1, ny + 1
            count = count + 1

            geom.iniP(count).pos(1) = (dble(i) - 1.0d0) * (1.0d0 / dble(nx)) - 0.5d0
            geom.iniP(count).pos(2) = (dble(j) - 1.0d0) * (1.0d0 / dble(ny)) - 0.5d0
            geom.iniP(count).pos(3) = geom.iniP(count).pos(1)**2.0d0 - geom.iniP(count).pos(2)**2.0d0
        end do
    end do

    ! Allocate face structure and set connectivity
    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    count = 0
    do i = 1, nx
        do j = 1, ny
            count = count + 1
            geom.face(count).poi(1) = (ny+1)*(i+0) + (j-1) + 1
            geom.face(count).poi(2) = (ny+1)*(i+0) + (j-1) + 2
            geom.face(count).poi(3) = (ny+1)*(i-1) + (j-1) + 2
            geom.face(count).poi(4) = (ny+1)*(i-1) + (j-1) + 1
        end do
    end do
end subroutine Exam_Open2D_Hyperbolic_Paraboloid_Quad

! ---------------------------------------------------------------------------------------

! Example of hyperbolic paraboloid with triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Hyperbolic_Paraboloid_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, j, count, n, nx, ny
    character :: pattern
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "62_Hyperbolic_Paraboloid_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hyperbolic Paraboloid Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n       = 3
    pattern = "\"
    nx      = n
    ny      = n

    geom.n_iniP = (nx + 1)*(ny + 1)
    geom.n_face = nx*ny*2

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    ! Allocate face structure and set connectivity
    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do
    ! Set position vector
    count = 0
    do i = 1, nx + 1
        do j = 1, ny + 1
            count = count + 1

            geom.iniP(count).pos(1) = (dble(i) - 1.0d0) * (1.0d0 / dble(nx)) - 0.5d0
            geom.iniP(count).pos(2) = (dble(j) - 1.0d0) * (1.0d0 / dble(ny)) - 0.5d0
            geom.iniP(count).pos(3) = geom.iniP(count).pos(1)**2.0d0 - geom.iniP(count).pos(2)**2.0d0
        end do
    end do

    ! Set connectivity
    count = 0
    do i = 1, nx
        do j = 1, ny
            if(pattern == "\") then
                count = count + 1
                geom.face(count).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(count).poi(2) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(count).poi(3) = (ny+1)*(i-1) + (j-1) + 1

                count = count + 1
                geom.face(count).poi(1) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(count).poi(2) = (ny+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (ny+1)*(i-1) + (j-1) + 1
            else if(pattern == "/") then
                count = count + 1
                geom.face(count).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(count).poi(2) = (ny+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (ny+1)*(i-1) + (j-1) + 1

                count = count + 1
                geom.face(count).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(count).poi(2) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(count).poi(3) = (ny+1)*(i-1) + (j-1) + 2
            end if
        end do
    end do
end subroutine Exam_Open2D_Hyperbolic_Paraboloid_Tri

! ---------------------------------------------------------------------------------------

! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Cross_Point(AP1, AP2, BP1, BP2, CP)
    double precision, intent(in) :: AP1(2), AP2(2), BP1(2), BP2(2)
    double precision, intent(out) :: CP(2)

    double precision :: t, buf

    buf = (BP2(2)-BP1(2))*(AP2(1)-AP1(1))-(BP2(1)-BP1(1))*(AP2(2)-AP1(2))
    t   = (BP2(1)-BP1(1))*(AP1(2)-BP1(2))-(BP2(2)-BP1(2))*(AP1(1)-BP1(1))
    t   = t/buf

    CP(1) = AP1(1) + t * (AP2(1)-AP1(1))
    CP(2) = AP1(2) + t * (AP2(2)-AP1(2))
end subroutine Exam_Open2D_Cross_Point

! ---------------------------------------------------------------------------------------

! Merge points and faces
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Merge_Point_Face_Quad(n, joint, conn, n_node, n_element, geom) 
    type(GeomType), intent(inout) :: geom
    double precision, intent(in) :: joint((n+1)**2, 3)
    integer, intent(in) :: n, conn(n*n, 4)
    integer, intent(inout) :: n_node, n_element

    integer :: i, j, k, m, flag

    ! Merge points
    m = n_node
    do i = 1, (n+1)**2
        ! Find a coninside node
        flag = 0
        do j=1, n_node
            if( Exam_Open2D_Comp_XYZ( joint(i,:), geom.iniP(j).pos(:) ) == 1) then
                flag = j; exit
            end if
        end do

        ! Copy points
        if(flag == 0) then
            m = m + 1
            geom.iniP(m).pos(:) = joint(i,:)
        end if
    end do

    n_node = m

    ! Merge faces
    do i = 1, n * n
        n_element = n_element + 1

        do j = 1, 4
            do k = 1, n_node
                if(Exam_Open2D_Comp_XYZ(joint(conn(i,j), :), geom.iniP(k).pos(:)) == 1) exit
            end do
            
            geom.face(n_element).poi(j) = k
        end do
    end do
end subroutine Exam_Open2D_Merge_Point_Face_Quad

! ---------------------------------------------------------------------------------------

! Merge points and faces
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Merge_Point_Face_Tri(n, joint, conn, n_node, n_element, geom)  
    type(GeomType), intent(inout) :: geom
    double precision, intent(in) :: joint((n+1)**2, 3)
    integer, intent(in) :: n, conn(n*n*2, 3)
    integer, intent(inout) :: n_node, n_element

    integer :: i, j, k, m, flag

    ! Merge points
    m = n_node
    do i = 1, (n+1)**2
        ! Find a coninside node
        flag = 0
        do j=1, n_node
            if( Exam_Open2D_Comp_XYZ( joint(i,:), geom.iniP(j).pos(:) ) == 1) then
                flag = j; exit
            end if
        end do

        ! Copy points
        if(flag == 0) then
            m = m + 1
            geom.iniP(m).pos(:) = joint(i,:)
        end if
    end do

    n_node = m

    ! Merge faces
    do i = 1, n * n * 2
        n_element = n_element + 1
        do j = 1, 3
            do k = 1, n_node
                if(Exam_Open2D_Comp_XYZ(joint(conn(i,j), :), geom.iniP(k).pos(:)) == 1) exit
            end do
            geom.face(n_element).poi(j) = k
        end do
    end do
end subroutine Exam_Open2D_Merge_Point_Face_Tri

! ---------------------------------------------------------------------------------------

! Last updated on Tuesday 30 August 2016 by Hyungmin
integer function Exam_Open2D_Comp_XYZ(a, b)
   double precision, intent(in) :: a(3), b(3) 
   double precision :: d

   Exam_Open2D_Comp_XYZ = 0
   d = dsqrt((a(1)-b(1))**2.0d0 + (a(2)-b(2))**2.0d0 + (a(3)-b(3))**2.0d0)
   if(d<1.0d-5) Exam_Open2D_Comp_XYZ = 1
end function Exam_Open2D_Comp_XYZ

! ---------------------------------------------------------------------------------------

! Example of circle generated by Distmesh using 0.4 factor
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_circle_Tri_Coarse(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Exam_Open2D_circle_Tri_Coarse"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tetrahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! allocate point, line and face structure
    geom.n_iniP = 19
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -42.8419d0,   0.3155d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -37.2731d0, -21.2527d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -36.9329d0,  21.5902d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -22.7947d0,  -0.0741d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ -21.6133d0, -37.1190d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -21.2310d0,  37.0880d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -11.4451d0, -19.8967d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -11.3495d0,  19.5828d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [   0.0000d0,  42.7183d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [   0.0000d0,  -0.2165d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [   0.0000d0, -42.9701d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [  11.3495d0,  19.5828d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [  11.4451d0, -19.8967d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [  21.2310d0,  37.0880d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [  21.6133d0, -37.1190d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [  22.7947d0,  -0.0741d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [  36.9329d0,  21.5902d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [  37.2731d0, -21.2527d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [  42.8419d0,   0.3155d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 12,  8, 10 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  9,  8, 12 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 11,  7,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 10,  7, 13 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 15, 18, 13 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 11, 15, 13 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 13,  7, 11 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10,  8,  4 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4,  7, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  1,  4,  3 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  3,  4,  8 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 16, 12, 10 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 10, 13, 16 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 16, 18, 19 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 16, 13, 18 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 17, 16, 19 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 12, 16, 17 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  2,  5,  7 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  2,  4,  1 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  7,  4,  2 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [  6,  8,  9 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [  6,  3,  8 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  9, 12, 14 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 12, 17, 14 ]
end subroutine Exam_Open2D_circle_Tri_Coarse

! ---------------------------------------------------------------------------------------

! Example of circle generated by Distmesh using 0.3 factor
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_circle_Tri_Fine(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Exam_Open2D_circle_Tri_Fine"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tetrahedron"

    ! Problem specified preset parameters
    if(para_vertex_design == "flat" .and. para_preset == "on") then
        para_junc_ang        = "min"    ! [opt, max, ave, min], Junction gap modification for different arm angle
        para_const_edge_mesh = "on"     ! [off, on], Constant edge length from polyhedra mesh
        para_unpaired_scaf   = "off"    ! [on, off], Unpaired scaffold nucleotides
        para_n_base_tn       = 7
    end if

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   =-1.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP =   41
    geom.n_face =   62

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -73.3614d0, -11.8370d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -72.5666d0,  16.2298d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -64.7977d0, -36.1929d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -62.3795d0,  40.4627d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ -54.2818d0,   1.2099d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -49.1997d0, -55.4151d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -46.3176d0, -20.6773d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -45.7579d0,  22.8010d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [ -45.0771d0,  59.0648d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [ -34.8025d0, -40.7945d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -32.8146d0,  42.8006d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -30.0051d0,   0.3843d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [ -27.6648d0, -68.6009d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [ -22.5483d0,  70.6993d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [ -19.5813d0,  22.4158d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [ -19.1358d0, -21.9646d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [ -12.3871d0, -48.0049d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [ -10.4760d0,  49.1958d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [  -4.8450d0,  -0.4156d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [  -2.9056d0, -73.7752d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [   2.8966d0,  25.1169d0, 0.0000d0 ]
    geom.iniP(22).pos(1:3) = [   3.2060d0,  74.0134d0, 0.0000d0 ]
    geom.iniP(23).pos(1:3) = [   4.8069d0, -25.7542d0, 0.0000d0 ]
    geom.iniP(24).pos(1:3) = [  10.0466d0, -50.2788d0, 0.0000d0 ]
    geom.iniP(25).pos(1:3) = [  11.9852d0,  48.6052d0, 0.0000d0 ]
    geom.iniP(26).pos(1:3) = [  18.1354d0,   0.8058d0, 0.0000d0 ]
    geom.iniP(27).pos(1:3) = [  22.6068d0, -70.1560d0, 0.0000d0 ]
    geom.iniP(28).pos(1:3) = [  25.4837d0,  21.0371d0, 0.0000d0 ]
    geom.iniP(29).pos(1:3) = [  26.5683d0, -18.3038d0, 0.0000d0 ]
    geom.iniP(30).pos(1:3) = [  28.8291d0,  68.0887d0, 0.0000d0 ]
    geom.iniP(31).pos(1:3) = [  30.8778d0, -41.2878d0, 0.0000d0 ]
    geom.iniP(32).pos(1:3) = [  33.2332d0,  42.4044d0, 0.0000d0 ]
    geom.iniP(33).pos(1:3) = [  44.1880d0, -58.8829d0, 0.0000d0 ]
    geom.iniP(34).pos(1:3) = [  46.7318d0,   0.3023d0, 0.0000d0 ]
    geom.iniP(35).pos(1:3) = [  48.8206d0,  22.5392d0, 0.0000d0 ]
    geom.iniP(36).pos(1:3) = [  49.2164d0, -22.9020d0, 0.0000d0 ]
    geom.iniP(37).pos(1:3) = [  51.4691d0,  52.8839d0, 0.0000d0 ]
    geom.iniP(38).pos(1:3) = [  60.9300d0, -41.2265d0, 0.0000d0 ]
    geom.iniP(39).pos(1:3) = [  65.9409d0,  32.8782d0, 0.0000d0 ]
    geom.iniP(40).pos(1:3) = [  71.8189d0, -15.9070d0, 0.0000d0 ]
    geom.iniP(41).pos(1:3) = [  73.1140d0,   8.4381d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 17, 13, 20 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 10, 13, 17 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  6, 13, 10 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 10, 17, 16 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 16,  7, 10 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 14,  9, 11 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  3,  7,  1 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10,  7,  3 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  3,  6, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 26, 21, 19 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 28, 21, 26 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 31, 38, 36 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 36, 29, 31 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 38, 40, 36 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  2,  8,  4 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 11,  9,  4 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  4,  8, 11 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  5,  2,  1 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  5,  8,  2 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  1,  7,  5 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 11,  8, 15 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 19, 21, 15 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 33, 38, 31 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 31, 27, 33 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 24, 27, 31 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 24, 17, 20 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 20, 27, 24 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 18, 14, 11 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 11, 15, 18 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 18, 15, 21 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 22, 14, 18 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 35, 32, 28 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 41, 39, 35 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 30, 32, 37 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 37, 35, 39 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 32, 35, 37 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 34, 26, 29 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 29, 36, 34 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 28, 26, 34 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 34, 35, 28 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 41, 35, 34 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 34, 40, 41 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 34, 36, 40 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [  8,  5, 12 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 12, 15,  8 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [  7, 16, 12 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 12,  5,  7 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 12, 16, 19 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 19, 15, 12 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 19, 16, 23 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 31, 29, 23 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 23, 24, 31 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 23, 16, 17 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 17, 24, 23 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 23, 26, 19 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 29, 26, 23 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 25, 30, 22 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 22, 18, 25 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 25, 32, 30 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 25, 18, 21 ]
    geom.face(61).n_poi = 3; allocate(geom.face(61).poi(3)); geom.face(61).poi(1:3) = [ 25, 21, 28 ]
    geom.face(62).n_poi = 3; allocate(geom.face(62).poi(3)); geom.face(62).poi(1:3) = [ 28, 32, 25 ]
end subroutine Exam_Open2D_circle_Tri_Fine

! ---------------------------------------------------------------------------------------

! Example of circle generated by Distmesh using 0.5 factor
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Ellipse_Tri_Coarse(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "13_Ellipse_Tri_Coarse"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Ellipse Tri Coarse"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   =-1.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 27
    geom.n_face = 36

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -79.9881d0,  -1.6440d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -70.6418d0,  19.5864d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -67.5987d0, -20.9006d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -54.9955d0,   3.7753d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ -49.5613d0,  32.1045d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -46.2168d0, -32.1627d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -40.6265d0, -12.5319d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -32.9305d0,  13.4195d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [ -24.7469d0,  38.7197d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [ -23.2681d0, -37.7112d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -19.9458d0, -11.3786d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -10.6472d0,  14.5326d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [  -0.0000d0, -39.4368d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [   0.0000d0,  40.6681d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [   0.0000d0, -12.8497d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [  10.6472d0,  14.5326d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [  19.9458d0, -11.3786d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [  23.2681d0, -37.7112d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [  24.7469d0,  38.7197d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [  32.9305d0,  13.4195d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [  40.6265d0, -12.5319d0, 0.0000d0 ]
    geom.iniP(22).pos(1:3) = [  46.2168d0, -32.1627d0, 0.0000d0 ]
    geom.iniP(23).pos(1:3) = [  49.5613d0,  32.1045d0, 0.0000d0 ]
    geom.iniP(24).pos(1:3) = [  54.9955d0,   3.7753d0, 0.0000d0 ]
    geom.iniP(25).pos(1:3) = [  67.5987d0, -20.9006d0, 0.0000d0 ]
    geom.iniP(26).pos(1:3) = [  70.6419d0,  19.5864d0, 0.0000d0 ]
    geom.iniP(27).pos(1:3) = [  79.9881d0,  -1.6440d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 17, 21, 20 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 20, 21, 24 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 27, 26, 24 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 24, 25, 27 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 21, 25, 24 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 23, 19, 20 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 20, 24, 23 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 23, 24, 26 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  1,  3,  4 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4,  2,  1 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  2,  4,  5 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 22, 25, 21 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 13, 15, 10 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 10, 15, 11 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 11, 15, 12 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 14,  9, 12 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  7, 10, 11 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  6, 10,  7 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  3,  6,  7 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  7,  4,  3 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 18, 15, 13 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 17, 15, 18 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 21, 17, 18 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 18, 22, 21 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 16, 15, 17 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 16, 12, 15 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 16, 17, 20 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 20, 19, 16 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 16, 19, 14 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 14, 12, 16 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  8,  5,  4 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [  4,  7,  8 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [  9,  5,  8 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [  8, 12,  9 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 11, 12,  8 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [  8,  7, 11 ]
end subroutine Exam_Open2D_Ellipse_Tri_Coarse

! ---------------------------------------------------------------------------------------

! Example of ellipse generated by Distmesh using 0.4 factor
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Ellipse_Tri_Fine(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "14_Ellipse_Tri_Fine"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Ellipse Tri Fine"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    prob.type_geo = "closed"
    if(para_fig_view == "preset") para_fig_view = "XY"

    ! Allocate point and face structure
    geom.n_iniP = 43
    geom.n_face = 65

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -96.9318d0,   7.5528d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -93.8328d0, -14.7675d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -83.7804d0,  25.3542d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -76.8633d0, -30.8500d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ -75.5667d0,  -4.1628d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -65.9050d0,  13.3487d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -62.7335d0,  37.5336d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -59.7035d0, -20.1866d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [ -52.7392d0, -41.7165d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [ -48.9984d0,  -3.0574d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -46.5389d0,  23.3181d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -38.3842d0,  44.9491d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [ -37.2164d0, -26.1347d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [ -32.1970d0,   8.4553d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [ -27.1676d0, -47.4715d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [ -23.6173d0, -11.0471d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [ -22.9875d0,  28.9760d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [ -12.6915d0,  48.4467d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [ -12.5057d0, -31.1852d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [ -10.5958d0,   8.6903d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [  -0.0000d0,  29.4515d0, 0.0000d0 ]
    geom.iniP(22).pos(1:3) = [   0.0000d0, -49.3949d0, 0.0000d0 ]
    geom.iniP(23).pos(1:3) = [   0.0000d0, -12.1473d0, 0.0000d0 ]
    geom.iniP(24).pos(1:3) = [  10.5958d0,   8.6903d0, 0.0000d0 ]
    geom.iniP(25).pos(1:3) = [  12.5057d0, -31.1852d0, 0.0000d0 ]
    geom.iniP(26).pos(1:3) = [  12.6915d0,  48.4467d0, 0.0000d0 ]
    geom.iniP(27).pos(1:3) = [  22.9875d0,  28.9760d0, 0.0000d0 ]
    geom.iniP(28).pos(1:3) = [  23.6173d0, -11.0471d0, 0.0000d0 ]
    geom.iniP(29).pos(1:3) = [  27.1676d0, -47.4715d0, 0.0000d0 ]
    geom.iniP(30).pos(1:3) = [  32.1970d0,   8.4553d0, 0.0000d0 ]
    geom.iniP(31).pos(1:3) = [  37.2164d0, -26.1347d0, 0.0000d0 ]
    geom.iniP(32).pos(1:3) = [  38.3842d0,  44.9491d0, 0.0000d0 ]
    geom.iniP(33).pos(1:3) = [  46.5389d0,  23.3181d0, 0.0000d0 ]
    geom.iniP(34).pos(1:3) = [  48.9984d0,  -3.0574d0, 0.0000d0 ]
    geom.iniP(35).pos(1:3) = [  52.7392d0, -41.7165d0, 0.0000d0 ]
    geom.iniP(36).pos(1:3) = [  59.7035d0, -20.1866d0, 0.0000d0 ]
    geom.iniP(37).pos(1:3) = [  62.7335d0,  37.5336d0, 0.0000d0 ]
    geom.iniP(38).pos(1:3) = [  65.9050d0,  13.3487d0, 0.0000d0 ]
    geom.iniP(39).pos(1:3) = [  75.5667d0,  -4.1628d0, 0.0000d0 ]
    geom.iniP(40).pos(1:3) = [  76.8633d0, -30.8500d0, 0.0000d0 ]
    geom.iniP(41).pos(1:3) = [  83.7804d0,  25.3542d0, 0.0000d0 ]
    geom.iniP(42).pos(1:3) = [  93.8328d0, -14.7675d0, 0.0000d0 ]
    geom.iniP(43).pos(1:3) = [  96.9318d0,   7.5528d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 40, 42, 39 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [ 43, 41, 39 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 39, 42, 43 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 26, 18, 21 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 28, 24, 23 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 40, 39, 36 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 34, 36, 39 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 17, 21, 18 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 18, 12, 17 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 17, 11, 14 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 17, 12, 11 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 24, 21, 20 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 20, 23, 24 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 20, 17, 14 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 21, 17, 20 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 28, 23, 25 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 22, 29, 25 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 13, 15, 19 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 19, 15, 22 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 22, 25, 19 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 19, 25, 23 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 40, 36, 35 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 10,  8, 13 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 14, 11, 10 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 10, 11,  6 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [  1,  2,  5 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [  2,  4,  5 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [  4,  8,  5 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [  5, 10,  6 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [  8, 10,  5 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  5,  3,  1 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [  6,  3,  5 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 38, 34, 39 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 41, 37, 38 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 38, 39, 41 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 27, 21, 24 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 26, 21, 27 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 27, 32, 26 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 11, 12,  7 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [  6, 11,  7 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [  7,  3,  6 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 13, 19, 16 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 16, 19, 23 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 14, 10, 16 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 16, 10, 13 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 16, 20, 14 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 23, 20, 16 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [  9, 15, 13 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 13,  8,  9 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [  9,  8,  4 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 31, 35, 36 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 31, 34, 28 ]
    geom.face(53).n_poi = 3; allocate(geom.face(53).poi(3)); geom.face(53).poi(1:3) = [ 36, 34, 31 ]
    geom.face(54).n_poi = 3; allocate(geom.face(54).poi(3)); geom.face(54).poi(1:3) = [ 29, 35, 31 ]
    geom.face(55).n_poi = 3; allocate(geom.face(55).poi(3)); geom.face(55).poi(1:3) = [ 28, 25, 31 ]
    geom.face(56).n_poi = 3; allocate(geom.face(56).poi(3)); geom.face(56).poi(1:3) = [ 31, 25, 29 ]
    geom.face(57).n_poi = 3; allocate(geom.face(57).poi(3)); geom.face(57).poi(1:3) = [ 30, 27, 24 ]
    geom.face(58).n_poi = 3; allocate(geom.face(58).poi(3)); geom.face(58).poi(1:3) = [ 30, 24, 28 ]
    geom.face(59).n_poi = 3; allocate(geom.face(59).poi(3)); geom.face(59).poi(1:3) = [ 28, 34, 30 ]
    geom.face(60).n_poi = 3; allocate(geom.face(60).poi(3)); geom.face(60).poi(1:3) = [ 34, 38, 33 ]
    geom.face(61).n_poi = 3; allocate(geom.face(61).poi(3)); geom.face(61).poi(1:3) = [ 33, 30, 34 ]
    geom.face(62).n_poi = 3; allocate(geom.face(62).poi(3)); geom.face(62).poi(1:3) = [ 33, 38, 37 ]
    geom.face(63).n_poi = 3; allocate(geom.face(63).poi(3)); geom.face(63).poi(1:3) = [ 37, 32, 33 ]
    geom.face(64).n_poi = 3; allocate(geom.face(64).poi(3)); geom.face(64).poi(1:3) = [ 32, 27, 33 ]
    geom.face(65).n_poi = 3; allocate(geom.face(65).poi(3)); geom.face(65).poi(1:3) = [ 27, 30, 33 ]
end subroutine Exam_Open2D_Ellipse_Tri_Fine

! ---------------------------------------------------------------------------------------

end module Exam_2D_Open