!
! ---------------------------------------------------------------------------------------
!
!                                Module - Exam_2D_Open
!
!                                                                    Updated : 2017/03/27
!
! Comments: 2D open geometry
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
! Reference: https://people.sc.fsu.edu/~jburkardt/m_src/distmesh/distmesh.html
!
module Exam_2D_Open

    use Data_Prob
    use Data_Geom

    use Para
    use Math
    use Mani

    implicit none

    public Exam_Open2D_Plate_Uniform_Quad       !  1. Plate with uniform quad mesh
    public Exam_Open2D_Plate_Uniform_Tri        !  2. Plate with uniform tri mesh
    public Exam_Open2D_Quarter_Circle_Quad      !  3. Quarter circle with quad mesh
    public Exam_Open2D_Disk_Quad                !  4. Disk with quad mesh
    public Exam_Open2D_Circle_Tri_Coarse        !  5. Circle with tri coarse mesh
    public Exam_Open2D_Ellipse_Tri_Coarse       !  6. Ellipse with tri coarse mesh
    public Exam_Open2D_L_Shape_Regular          !  7. L-shape with regular mesh
    public Exam_Open2D_Hexagonal_Mesh           !  8. Hexagonal mesh
    public Exam_Open2D_Honeycomb                !  9. Honeycomb
    public Exam_Open2D_Stickman                 ! 10. 2D stickman

    public Exam_Open2D_Quarter_Circle_Tri           ! Quarter circle with tri mesh
    public Exam_Open2D_Disk_Tri                     ! Disk with tri mesh
    public Exam_Open2D_Circle_Tri_Fine              ! Circle with tri fine mesh
    public Exam_Open2D_Ellipse_Tri_Fine             ! Ellipse with tri fine mesh
    public Exam_Open2D_L_Shape_Irregular            ! L-shape with irregular mesh

    !public Exam_Open2D_Plate_Distorted_Quad         ! Plate with distorted quad mesh
    !public Exam_Open2D_Plate_Distorted_Tri          ! Plate with distorted tri mesh
    !public Exam_Open2D_Hyperbolic_Paraboloid_Quad   ! Hyperbolic paraboloid with quad mesh
    !public Exam_Open2D_Hyperbolic_Paraboloid_Tri    ! Hyperbolic paraboloid with tri mesh

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

    prob.name_file = "01_Plate_Uniform_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Uniform Quad"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! Set options
    !n  = 2
    nx = 4
    ny = 3

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! Set options
    !n  = 2
    nx = 3
    ny = 4
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

! Example of quarter circle with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Quarter_Circle_Quad(prob, geom)
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

    prob.name_file = "03_Quarter_Circle_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Quarter Circle Quad"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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
end subroutine Exam_Open2D_Quarter_Circle_Quad

! ---------------------------------------------------------------------------------------

! Example of disk with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Disk_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: i_rad, o_rad, ang, rad
    integer :: i, j, index, n, nx, nr
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "04_Disk_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Disk Quad"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    n  = 2
    nx = n
    nr = n * 5

    o_rad = 4.0d0       ! Outer radius
    i_rad = 1.8d0       ! Internal radius

    geom.n_iniP = (nx + 1) * nr
    geom.n_face = nx * nr

    ! Allocate point structure and set position vector
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    index = 0
    do j = 1, nr
        do i = 1, nx + 1
            index  = index + 1
            ang = 2.0d0 * pi / dble(nr) * dble(j - 1)
            rad = ((o_rad - i_rad) / dble(nx) * dble(i - 1)) + i_rad

            geom.iniP(index).pos(1) = rad * dcos(ang)
            geom.iniP(index).pos(2) = rad * dsin(ang)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nr
        do j = 1, nx
            if(i /= nr) then
                index = index + 1
                geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(index).poi(3) = (nx+1)*(i-0) + (j-1) + 2
                geom.face(index).poi(4) = (nx+1)*(i-0) + (j-1) + 1
            else if(i == nr) then
                index = index + 1
                geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(index).poi(3) = (nx+1)*(1-1) + (j-1) + 2
                geom.face(index).poi(4) = (nx+1)*(1-1) + (j-1) + 1
            end if
        end do
    end do
end subroutine Exam_Open2D_Disk_Quad

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

    prob.name_file = "05_Circle_Tri_Coarse"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Circle Tri Coarse"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
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

! Example of circle generated by Distmesh using 0.5 factor
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open2D_Ellipse_Tri_Coarse(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "06_Ellipse_Tri_Coarse"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Ellipse Tri Coarse"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
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

! Example of L-shape with regular mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_L_Shape_Regular(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "07_L_Shape_Regular"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "L Shape Regular"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    geom.n_iniP = 21
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -34.2857d0, -34.2857d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -14.2857d0, -34.2857d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [   5.7143d0, -34.2857d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [  25.7143d0, -34.2857d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [  45.7143d0, -34.2857d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -34.2857d0, -14.2857d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -14.2857d0, -14.2857d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [   5.7143d0, -14.2857d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [  25.7143d0, -14.2857d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [  45.7143d0, -14.2857d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -34.2857d0,   5.7143d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -14.2857d0,   5.7143d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [   5.7143d0,   5.7143d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [  25.7143d0,   5.7143d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [  45.7143d0,   5.7143d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [ -34.2857d0,  25.7143d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [ -14.2857d0,  25.7143d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [   5.7143d0,  25.7143d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [ -34.2857d0,  45.7143d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [ -14.2857d0,  45.7143d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [   5.7143d0,  45.7143d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  6 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  7,  6,  2 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  3,  7 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8,  7,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  3,  4,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  9,  8,  4 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  4,  5,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10,  9,  5 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  6,  7, 11 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 12, 11,  7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7,  8, 12 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 13, 12,  8 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  8,  9, 13 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 14, 13,  9 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [  9, 10, 14 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 15, 14, 10 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 11, 12, 16 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 17, 16, 12 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 12, 13, 17 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 18, 17, 13 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 16, 17, 19 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 20, 19, 17 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 17, 18, 20 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 21, 20, 18 ]
end subroutine Exam_Open2D_L_Shape_Regular

! ---------------------------------------------------------------------------------------

! Example of hexagonal mesh
! Last updated on Thu 23 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Hexagonal_Mesh(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "08_Hexagonal_Mesh"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hexagonal Mesh"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    geom.n_iniP = 38
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set the position vector
    geom.iniP( 1).pos(1:3) = [ -10.0000d0, -82.0445d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [  10.0000d0, -82.0445d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -40.0000d0, -64.7240d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -20.0000d0, -64.7240d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [  20.0000d0, -64.7240d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [  40.0000d0, -64.7240d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -50.0000d0, -47.4035d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -10.0000d0, -47.4035d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [  10.0000d0, -47.4035d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [  50.0000d0, -47.4035d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -40.0000d0, -30.0830d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -20.0000d0, -30.0830d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [  20.0000d0, -30.0830d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [  40.0000d0, -30.0830d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [ -50.0000d0, -12.7625d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [ -10.0000d0, -12.7625d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [  10.0000d0, -12.7625d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [  50.0000d0, -12.7625d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [ -40.0000d0,   4.5580d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [ -20.0000d0,   4.5580d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [  20.0000d0,   4.5580d0, 0.0000d0 ]
    geom.iniP(22).pos(1:3) = [  40.0000d0,   4.5580d0, 0.0000d0 ]
    geom.iniP(23).pos(1:3) = [ -50.0000d0,  21.8785d0, 0.0000d0 ]
    geom.iniP(24).pos(1:3) = [ -10.0000d0,  21.8785d0, 0.0000d0 ]
    geom.iniP(25).pos(1:3) = [  10.0000d0,  21.8785d0, 0.0000d0 ]
    geom.iniP(26).pos(1:3) = [  50.0000d0,  21.8785d0, 0.0000d0 ]
    geom.iniP(27).pos(1:3) = [ -40.0000d0,  39.1991d0, 0.0000d0 ]
    geom.iniP(28).pos(1:3) = [ -20.0000d0,  39.1991d0, 0.0000d0 ]
    geom.iniP(29).pos(1:3) = [  20.0000d0,  39.1991d0, 0.0000d0 ]
    geom.iniP(30).pos(1:3) = [  40.0000d0,  39.1991d0, 0.0000d0 ]
    geom.iniP(31).pos(1:3) = [ -50.0000d0,  56.5196d0, 0.0000d0 ]
    geom.iniP(32).pos(1:3) = [ -10.0000d0,  56.5196d0, 0.0000d0 ]
    geom.iniP(33).pos(1:3) = [  10.0000d0,  56.5196d0, 0.0000d0 ]
    geom.iniP(34).pos(1:3) = [  50.0000d0,  56.5196d0, 0.0000d0 ]
    geom.iniP(35).pos(1:3) = [ -40.0000d0,  73.8401d0, 0.0000d0 ]
    geom.iniP(36).pos(1:3) = [ -20.0000d0,  73.8401d0, 0.0000d0 ]
    geom.iniP(37).pos(1:3) = [  20.0000d0,  73.8401d0, 0.0000d0 ]
    geom.iniP(38).pos(1:3) = [  40.0000d0,  73.8401d0, 0.0000d0 ]

    ! Set connectivity
    geom.face( 1).n_poi = 6; allocate(geom.face( 1).poi(6)); geom.face( 1).poi(1:6) = [  1,  2,  5,  9,  8,  4 ]
    geom.face( 2).n_poi = 6; allocate(geom.face( 2).poi(6)); geom.face( 2).poi(1:6) = [  3,  4,  8, 12, 11,  7 ]
    geom.face( 3).n_poi = 6; allocate(geom.face( 3).poi(6)); geom.face( 3).poi(1:6) = [  5,  6, 10, 14, 13,  9 ]
    geom.face( 4).n_poi = 6; allocate(geom.face( 4).poi(6)); geom.face( 4).poi(1:6) = [  8,  9, 13, 17, 16, 12 ]
    geom.face( 5).n_poi = 6; allocate(geom.face( 5).poi(6)); geom.face( 5).poi(1:6) = [ 11, 12, 16, 20, 19, 15 ]
    geom.face( 6).n_poi = 6; allocate(geom.face( 6).poi(6)); geom.face( 6).poi(1:6) = [ 13, 14, 18, 22, 21, 17 ]
    geom.face( 7).n_poi = 6; allocate(geom.face( 7).poi(6)); geom.face( 7).poi(1:6) = [ 16, 17, 21, 25, 24, 20 ]
    geom.face( 8).n_poi = 6; allocate(geom.face( 8).poi(6)); geom.face( 8).poi(1:6) = [ 19, 20, 24, 28, 27, 23 ]
    geom.face( 9).n_poi = 6; allocate(geom.face( 9).poi(6)); geom.face( 9).poi(1:6) = [ 21, 22, 26, 30, 29, 25 ]
    geom.face(10).n_poi = 6; allocate(geom.face(10).poi(6)); geom.face(10).poi(1:6) = [ 24, 25, 29, 33, 32, 28 ]
    geom.face(11).n_poi = 6; allocate(geom.face(11).poi(6)); geom.face(11).poi(1:6) = [ 27, 28, 32, 36, 35, 31 ]
    geom.face(12).n_poi = 6; allocate(geom.face(12).poi(6)); geom.face(12).poi(1:6) = [ 29, 30, 34, 38, 37, 33 ]
end subroutine Exam_Open2D_Hexagonal_Mesh

! ---------------------------------------------------------------------------------------

! Example of honeycomb
! Last updated on Thu 23 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Honeycomb(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "09_Honeycomb"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Honeycomb"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! Allocate point and face structure
    geom.n_iniP =   30
    geom.n_face =  36

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -51.9618d0,  10.0000d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -17.3206d0, -10.0000d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -34.6411d0,   0.0000d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -51.9618d0, -10.0000d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [   0.0000d0, -20.0001d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -34.6411d0, -40.0002d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -17.3206d0, -30.0001d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -34.6411d0, -20.0001d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [ -17.3206d0, -50.0002d0, 0.0000d0 ]
    geom.iniP(10).pos(1:3) = [ -17.3206d0,  10.0000d0, 0.0000d0 ]
    geom.iniP(11).pos(1:3) = [ -34.6411d0,  40.0002d0, 0.0000d0 ]
    geom.iniP(12).pos(1:3) = [ -34.6411d0,  20.0001d0, 0.0000d0 ]
    geom.iniP(13).pos(1:3) = [ -17.3206d0,  30.0001d0, 0.0000d0 ]
    geom.iniP(14).pos(1:3) = [   0.0000d0,  20.0001d0, 0.0000d0 ]
    geom.iniP(15).pos(1:3) = [ -17.3206d0,  50.0002d0, 0.0000d0 ]
    geom.iniP(16).pos(1:3) = [  17.3205d0, -50.0002d0, 0.0000d0 ]
    geom.iniP(17).pos(1:3) = [  34.6412d0, -40.0002d0, 0.0000d0 ]
    geom.iniP(18).pos(1:3) = [  34.6412d0, -20.0001d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [   0.0000d0, -40.0002d0, 0.0000d0 ]
    geom.iniP(20).pos(1:3) = [  17.3205d0, -30.0001d0, 0.0000d0 ]
    geom.iniP(21).pos(1:3) = [  51.9618d0, -10.0000d0, 0.0000d0 ]
    geom.iniP(22).pos(1:3) = [  51.9618d0,  10.0000d0, 0.0000d0 ]
    geom.iniP(23).pos(1:3) = [  34.6412d0,  20.0001d0, 0.0000d0 ]
    geom.iniP(24).pos(1:3) = [  17.3205d0, -10.0000d0, 0.0000d0 ]
    geom.iniP(25).pos(1:3) = [  34.6412d0,   0.0000d0, 0.0000d0 ]
    geom.iniP(26).pos(1:3) = [  34.6412d0,  40.0002d0, 0.0000d0 ]
    geom.iniP(27).pos(1:3) = [  17.3205d0,  50.0002d0, 0.0000d0 ]
    geom.iniP(28).pos(1:3) = [   0.0000d0,  40.0002d0, 0.0000d0 ]
    geom.iniP(29).pos(1:3) = [  17.3205d0,  10.0000d0, 0.0000d0 ]
    geom.iniP(30).pos(1:3) = [  17.3205d0,  30.0001d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  8,  2,  3 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  2, 10,  3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 10, 12,  3 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 12,  1,  3 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  1,  4,  3 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  4,  8,  3 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  9, 19,  7 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 19,  5,  7 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  5,  2,  7 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  2,  8,  7 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  8,  6,  7 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  6,  9,  7 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 10, 14, 13 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 14, 28, 13 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 28, 15, 13 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 15, 11, 13 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 11, 12, 13 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 12, 10, 13 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 16, 17, 20 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 17, 18, 20 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 18, 24, 20 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 24,  5, 20 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  5, 19, 20 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 19, 16, 20 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [ 18, 21, 25 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 21, 22, 25 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 22, 23, 25 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 23, 29, 25 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 29, 24, 25 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 24, 18, 25 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [ 29, 23, 30 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [ 23, 26, 30 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 26, 27, 30 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 27, 28, 30 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 28, 14, 30 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 14, 29, 30 ]
end subroutine Exam_Open2D_Honeycomb

! ---------------------------------------------------------------------------------------

! Example of stickman
! Last updated on Thu 23 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Stickman(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "10_Stickman"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Stickman"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    geom.n_iniP = 48
    geom.n_face = 52

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set the position vector
    geom.iniP( 1).pos(1:3) = [ -16.42769d0, -18.33558d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  -5.47590d0, -18.33558d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [   5.47590d0, -18.33558d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  16.42769d0, -18.33558d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [ -16.42769d0,  -7.38379d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -5.47590d0,  -7.38379d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   5.47590d0,  -7.38379d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  16.42769d0,  -7.38379d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -16.42769d0,   3.56800d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  -5.47590d0,   3.56800d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [   5.47590d0,   3.56800d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [  16.42769d0,   3.56800d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -16.42769d0,  14.51980d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  -5.47590d0,  14.51980d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [   5.47590d0,  14.51980d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  16.42769d0,  14.51980d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [ -16.42769d0,  25.47159d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  -5.47590d0,  25.47159d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [   5.47590d0,  25.47159d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  16.42769d0,  25.47159d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [  -7.38151d0, -29.44070d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [  -9.31998d0, -40.54582d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [ -11.26940d0, -51.63999d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [ -18.32235d0, -28.80550d0, 0.0d0 ]
    geom.iniP(25).pos(1:3) = [ -20.19511d0, -39.27541d0, 0.0d0 ]
    geom.iniP(26).pos(1:3) = [ -22.04596d0, -49.74533d0, 0.0d0 ]
    geom.iniP(27).pos(1:3) = [  18.32235d0, -28.80550d0, 0.0d0 ]
    geom.iniP(28).pos(1:3) = [  20.19511d0, -39.27541d0, 0.0d0 ]
    geom.iniP(29).pos(1:3) = [  22.04596d0, -49.74533d0, 0.0d0 ]
    geom.iniP(30).pos(1:3) = [   7.38151d0, -29.44070d0, 0.0d0 ]
    geom.iniP(31).pos(1:3) = [   9.31998d0, -40.54582d0, 0.0d0 ]
    geom.iniP(32).pos(1:3) = [  11.26940d0, -51.63999d0, 0.0d0 ]
    geom.iniP(33).pos(1:3) = [  26.86475d0,  20.08331d0, 0.0d0 ]
    geom.iniP(34).pos(1:3) = [  37.26895d0,  25.80014d0, 0.0d0 ]
    geom.iniP(35).pos(1:3) = [  47.61840d0,  31.68126d0, 0.0d0 ]
    geom.iniP(36).pos(1:3) = [  24.95914d0,  30.85987d0, 0.0d0 ]
    geom.iniP(37).pos(1:3) = [  33.52344d0,  36.09483d0, 0.0d0 ]
    geom.iniP(38).pos(1:3) = [  42.14250d0,  41.16551d0, 0.0d0 ]
    geom.iniP(39).pos(1:3) = [ -24.95914d0,   9.13151d0, 0.0d0 ]
    geom.iniP(40).pos(1:3) = [ -33.52344d0,   3.89656d0, 0.0d0 ]
    geom.iniP(41).pos(1:3) = [ -42.14250d0,  -1.17412d0, 0.0d0 ]
    geom.iniP(42).pos(1:3) = [ -26.86475d0,  19.90808d0, 0.0d0 ]
    geom.iniP(43).pos(1:3) = [ -37.26895d0,  14.19124d0, 0.0d0 ]
    geom.iniP(44).pos(1:3) = [ -47.61840d0,   8.31013d0, 0.0d0 ]
    geom.iniP(45).pos(1:3) = [  -6.57108d0,  36.42338d0, 0.0d0 ]
    geom.iniP(46).pos(1:3) = [   6.57108d0,  36.42338d0, 0.0d0 ]
    geom.iniP(47).pos(1:3) = [  -5.25686d0,  47.37518d0, 0.0d0 ]
    geom.iniP(48).pos(1:3) = [   5.25686d0,  47.37518d0, 0.0d0 ]

    ! Set face connnectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  1,  2,  5 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  5,  2,  6 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  2,  3,  6 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  3,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  3,  4,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  7,  4,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [  5,  6,  9 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [  9,  6, 10 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  6,  7, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 10,  7, 11 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7,  8, 11 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 11,  8, 12 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  9, 10, 13 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 13, 10, 14 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 10, 11, 14 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 14, 11, 15 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 11, 12, 15 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [ 15, 12, 16 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 13, 14, 17 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [ 17, 14, 18 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [ 14, 15, 18 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [ 18, 15, 19 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [ 15, 16, 19 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 19, 16, 20 ]
    geom.face(25).n_poi = 3; allocate(geom.face(25).poi(3)); geom.face(25).poi(1:3) = [  2,  1, 21 ]
    geom.face(26).n_poi = 3; allocate(geom.face(26).poi(3)); geom.face(26).poi(1:3) = [ 21,  1, 24 ]
    geom.face(27).n_poi = 3; allocate(geom.face(27).poi(3)); geom.face(27).poi(1:3) = [ 21, 24, 22 ]
    geom.face(28).n_poi = 3; allocate(geom.face(28).poi(3)); geom.face(28).poi(1:3) = [ 22, 24, 25 ]
    geom.face(29).n_poi = 3; allocate(geom.face(29).poi(3)); geom.face(29).poi(1:3) = [ 22, 25, 23 ]
    geom.face(30).n_poi = 3; allocate(geom.face(30).poi(3)); geom.face(30).poi(1:3) = [ 23, 25, 26 ]
    geom.face(31).n_poi = 3; allocate(geom.face(31).poi(3)); geom.face(31).poi(1:3) = [  3, 30,  4 ]
    geom.face(32).n_poi = 3; allocate(geom.face(32).poi(3)); geom.face(32).poi(1:3) = [  4, 30, 27 ]
    geom.face(33).n_poi = 3; allocate(geom.face(33).poi(3)); geom.face(33).poi(1:3) = [ 30, 31, 27 ]
    geom.face(34).n_poi = 3; allocate(geom.face(34).poi(3)); geom.face(34).poi(1:3) = [ 27, 31, 28 ]
    geom.face(35).n_poi = 3; allocate(geom.face(35).poi(3)); geom.face(35).poi(1:3) = [ 31, 32, 28 ]
    geom.face(36).n_poi = 3; allocate(geom.face(36).poi(3)); geom.face(36).poi(1:3) = [ 28, 32, 29 ]
    geom.face(37).n_poi = 3; allocate(geom.face(37).poi(3)); geom.face(37).poi(1:3) = [ 16, 33, 20 ]
    geom.face(38).n_poi = 3; allocate(geom.face(38).poi(3)); geom.face(38).poi(1:3) = [ 20, 33, 36 ]
    geom.face(39).n_poi = 3; allocate(geom.face(39).poi(3)); geom.face(39).poi(1:3) = [ 33, 34, 36 ]
    geom.face(40).n_poi = 3; allocate(geom.face(40).poi(3)); geom.face(40).poi(1:3) = [ 36, 34, 37 ]
    geom.face(41).n_poi = 3; allocate(geom.face(41).poi(3)); geom.face(41).poi(1:3) = [ 34, 35, 37 ]
    geom.face(42).n_poi = 3; allocate(geom.face(42).poi(3)); geom.face(42).poi(1:3) = [ 37, 35, 38 ]
    geom.face(43).n_poi = 3; allocate(geom.face(43).poi(3)); geom.face(43).poi(1:3) = [ 17, 42, 13 ]
    geom.face(44).n_poi = 3; allocate(geom.face(44).poi(3)); geom.face(44).poi(1:3) = [ 13, 42, 39 ]
    geom.face(45).n_poi = 3; allocate(geom.face(45).poi(3)); geom.face(45).poi(1:3) = [ 42, 43, 39 ]
    geom.face(46).n_poi = 3; allocate(geom.face(46).poi(3)); geom.face(46).poi(1:3) = [ 39, 43, 40 ]
    geom.face(47).n_poi = 3; allocate(geom.face(47).poi(3)); geom.face(47).poi(1:3) = [ 43, 44, 40 ]
    geom.face(48).n_poi = 3; allocate(geom.face(48).poi(3)); geom.face(48).poi(1:3) = [ 40, 44, 41 ]
    geom.face(49).n_poi = 3; allocate(geom.face(49).poi(3)); geom.face(49).poi(1:3) = [ 18, 19, 45 ]
    geom.face(50).n_poi = 3; allocate(geom.face(50).poi(3)); geom.face(50).poi(1:3) = [ 45, 19, 46 ]
    geom.face(51).n_poi = 3; allocate(geom.face(51).poi(3)); geom.face(51).poi(1:3) = [ 45, 46, 47 ]
    geom.face(52).n_poi = 3; allocate(geom.face(52).poi(3)); geom.face(52).poi(1:3) = [ 47, 46, 48 ]
end subroutine Exam_Open2D_Stickman

! ---------------------------------------------------------------------------------------

! Example of quarter circle with tri mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Quarter_Circle_Tri(prob, geom)
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

    prob.name_file = "04_Quarter_Circle_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Quarter Circle Tri"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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
end subroutine Exam_Open2D_Quarter_Circle_Tri

! ---------------------------------------------------------------------------------------

! Example of disk with tri mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Disk_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: i_rad, o_rad, ang, rad
    integer :: i, j, index, n, nx, nr
    character :: pn
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "06_Disk_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Disk Tri"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    n  = 2
    nx = n
    nr = n * 5
    pn = "\"

    o_rad = 4.0d0       ! Outer radius
    i_rad = 1.8d0       ! Internal radius

    geom.n_iniP = (nx + 1) * nr
    geom.n_face = nx * nr * 2

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    index = 0
    do j = 1, nr
        do i = 1, nx + 1
            index = index + 1

            ang = 2.0d0 * pi / dble(nr) * dble(j - 1)
            rad = ((o_rad - i_rad) / dble(nx) * dble(i - 1)) + i_rad

            geom.iniP(index).pos(1) = rad * dcos(ang)
            geom.iniP(index).pos(2) = rad * dsin(ang)
            geom.iniP(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nr
        do j = 1, nx
            if(i /= nr) then
                if(pn == "\") then

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(i-0) + (j-1) + 1

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(2) = (nx+1)*(i-0) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(i-0) + (j-1) + 1
                else if(pn == "\") then

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(i-0) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(i-0) + (j-1) + 1

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(i-0) + (j-1) + 2
                end if
            else if(i == nr) then
                if(pn == "\") then

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(1-1) + (j-1) + 1

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(2) = (nx+1)*(1-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(1-1) + (j-1) + 1
                else if(pn == "\") then

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(1-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(1-1) + (j-1) + 1

                    index = index + 1
                    geom.face(index).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                    geom.face(index).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                    geom.face(index).poi(3) = (nx+1)*(1-1) + (j-1) + 2
                end if
            end if
        end do
    end do
end subroutine Exam_Open2D_Disk_Tri

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

    prob.name_file = "08_Circle_Tri_Fine"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Circle Tri Fine"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
    geom.n_iniP = 41
    geom.n_face = 62

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

! Example of ellipse generated by Distmesh using 0.4 factor
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Ellipse_Tri_Fine(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "10_Ellipse_Tri_Fine"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Ellipse Tri Fine"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    ! The number of points and faces
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

! Example of L-shape with irregular
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open2D_L_Shape_Irregular(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "12_L_Shape_Irregular"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "L Shape Irregular"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    geom.n_iniP =  100
    geom.n_face =  156

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(  1).pos(1:3) = [ -104.0849d0, -107.3907d0, 0.0000d0 ]
    geom.iniP(  2).pos(1:3) = [  146.3091d0, -107.3907d0, 0.0000d0 ]
    geom.iniP(  3).pos(1:3) = [  146.3091d0,   17.8063d0, 0.0000d0 ]
    geom.iniP(  4).pos(1:3) = [   21.1121d0,   17.8063d0, 0.0000d0 ]
    geom.iniP(  5).pos(1:3) = [   21.1121d0,  143.0032d0, 0.0000d0 ]
    geom.iniP(  6).pos(1:3) = [ -104.0849d0,  143.0032d0, 0.0000d0 ]
    geom.iniP(  7).pos(1:3) = [ -104.0849d0,  -85.3976d0, 0.0000d0 ]
    geom.iniP(  8).pos(1:3) = [ -104.0849d0,  -62.6196d0, 0.0000d0 ]
    geom.iniP(  9).pos(1:3) = [ -104.0849d0,  -40.1192d0, 0.0000d0 ]
    geom.iniP( 10).pos(1:3) = [ -104.0849d0,  -17.0284d0, 0.0000d0 ]
    geom.iniP( 11).pos(1:3) = [ -104.0849d0,    5.4728d0, 0.0000d0 ]
    geom.iniP( 12).pos(1:3) = [ -104.0849d0,   28.5741d0, 0.0000d0 ]
    geom.iniP( 13).pos(1:3) = [ -104.0849d0,   51.3316d0, 0.0000d0 ]
    geom.iniP( 14).pos(1:3) = [ -104.0849d0,   74.3921d0, 0.0000d0 ]
    geom.iniP( 15).pos(1:3) = [ -104.0849d0,   97.1146d0, 0.0000d0 ]
    geom.iniP( 16).pos(1:3) = [ -104.0849d0,  119.9825d0, 0.0000d0 ]
    geom.iniP( 17).pos(1:3) = [  -83.4878d0,  143.0032d0, 0.0000d0 ]
    geom.iniP( 18).pos(1:3) = [  -81.1939d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 19).pos(1:3) = [  -73.6894d0,  -85.1848d0, 0.0000d0 ]
    geom.iniP( 20).pos(1:3) = [  -82.9189d0,  -63.1748d0, 0.0000d0 ]
    geom.iniP( 21).pos(1:3) = [  -74.3278d0,  -40.8403d0, 0.0000d0 ]
    geom.iniP( 22).pos(1:3) = [  -83.0265d0,  -18.1901d0, 0.0000d0 ]
    geom.iniP( 23).pos(1:3) = [  -75.5905d0,    4.6463d0, 0.0000d0 ]
    geom.iniP( 24).pos(1:3) = [  -83.0822d0,   27.8965d0, 0.0000d0 ]
    geom.iniP( 25).pos(1:3) = [  -75.6186d0,   50.8833d0, 0.0000d0 ]
    geom.iniP( 26).pos(1:3) = [  -82.6556d0,   73.6899d0, 0.0000d0 ]
    geom.iniP( 27).pos(1:3) = [  -74.9044d0,   95.7833d0, 0.0000d0 ]
    geom.iniP( 28).pos(1:3) = [  -80.9939d0,  118.2975d0, 0.0000d0 ]
    geom.iniP( 29).pos(1:3) = [  -62.8829d0,  143.0032d0, 0.0000d0 ]
    geom.iniP( 30).pos(1:3) = [  -55.7317d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 31).pos(1:3) = [  -44.6483d0,  -84.7518d0, 0.0000d0 ]
    geom.iniP( 32).pos(1:3) = [  -57.0548d0,  -62.7116d0, 0.0000d0 ]
    geom.iniP( 33).pos(1:3) = [  -45.1985d0,  -40.2681d0, 0.0000d0 ]
    geom.iniP( 34).pos(1:3) = [  -57.9642d0,  -18.4281d0, 0.0000d0 ]
    geom.iniP( 35).pos(1:3) = [  -48.1252d0,    4.3575d0, 0.0000d0 ]
    geom.iniP( 36).pos(1:3) = [  -58.5145d0,   27.6753d0, 0.0000d0 ]
    geom.iniP( 37).pos(1:3) = [  -48.2792d0,   50.8161d0, 0.0000d0 ]
    geom.iniP( 38).pos(1:3) = [  -56.8227d0,   73.1672d0, 0.0000d0 ]
    geom.iniP( 39).pos(1:3) = [  -47.1925d0,   95.2388d0, 0.0000d0 ]
    geom.iniP( 40).pos(1:3) = [  -55.9556d0,  118.1338d0, 0.0000d0 ]
    geom.iniP( 41).pos(1:3) = [  -41.4100d0,  143.0032d0, 0.0000d0 ]
    geom.iniP( 42).pos(1:3) = [  -29.4498d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 43).pos(1:3) = [  -16.5712d0,  -83.9986d0, 0.0000d0 ]
    geom.iniP( 44).pos(1:3) = [  -29.7684d0,  -62.1370d0, 0.0000d0 ]
    geom.iniP( 45).pos(1:3) = [  -16.8814d0,  -39.4998d0, 0.0000d0 ]
    geom.iniP( 46).pos(1:3) = [  -30.9884d0,  -17.3621d0, 0.0000d0 ]
    geom.iniP( 47).pos(1:3) = [  -20.9735d0,    6.1689d0, 0.0000d0 ]
    geom.iniP( 48).pos(1:3) = [  -32.4620d0,   26.9436d0, 0.0000d0 ]
    geom.iniP( 49).pos(1:3) = [  -23.1061d0,   48.8409d0, 0.0000d0 ]
    geom.iniP( 50).pos(1:3) = [  -28.6495d0,   72.7932d0, 0.0000d0 ]
    geom.iniP( 51).pos(1:3) = [  -22.2007d0,   96.9358d0, 0.0000d0 ]
    geom.iniP( 52).pos(1:3) = [  -30.0711d0,  119.4122d0, 0.0000d0 ]
    geom.iniP( 53).pos(1:3) = [  -19.0191d0,  143.0032d0, 0.0000d0 ]
    geom.iniP( 54).pos(1:3) = [   -3.2699d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 55).pos(1:3) = [   10.2373d0,  -83.2952d0, 0.0000d0 ]
    geom.iniP( 56).pos(1:3) = [   -2.0004d0,  -60.1991d0, 0.0000d0 ]
    geom.iniP( 57).pos(1:3) = [   12.1414d0,  -35.8825d0, 0.0000d0 ]
    geom.iniP( 58).pos(1:3) = [   -3.7155d0,  -17.3483d0, 0.0000d0 ]
    geom.iniP( 59).pos(1:3) = [    2.4635d0,    5.2952d0, 0.0000d0 ]
    geom.iniP( 60).pos(1:3) = [   -3.9081d0,   30.9332d0, 0.0000d0 ]
    geom.iniP( 61).pos(1:3) = [   -1.8629d0,   59.3932d0, 0.0000d0 ]
    geom.iniP( 62).pos(1:3) = [   -1.8254d0,   85.6454d0, 0.0000d0 ]
    geom.iniP( 63).pos(1:3) = [   21.1121d0,   95.9976d0, 0.0000d0 ]
    geom.iniP( 64).pos(1:3) = [   -3.1287d0,  115.5400d0, 0.0000d0 ]
    geom.iniP( 65).pos(1:3) = [    0.9809d0,  143.0032d0, 0.0000d0 ]
    geom.iniP( 66).pos(1:3) = [   22.5006d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 67).pos(1:3) = [   36.1306d0,  -82.6261d0, 0.0000d0 ]
    geom.iniP( 68).pos(1:3) = [   24.5357d0,  -58.6468d0, 0.0000d0 ]
    geom.iniP( 69).pos(1:3) = [   38.5154d0,  -32.9882d0, 0.0000d0 ]
    geom.iniP( 70).pos(1:3) = [   24.8134d0,   -9.7106d0, 0.0000d0 ]
    geom.iniP( 71).pos(1:3) = [   43.4690d0,   17.8063d0, 0.0000d0 ]
    geom.iniP( 72).pos(1:3) = [   21.1121d0,   43.6398d0, 0.0000d0 ]
    geom.iniP( 73).pos(1:3) = [   21.1121d0,   70.0289d0, 0.0000d0 ]
    geom.iniP( 74).pos(1:3) = [   21.1121d0,  119.6866d0, 0.0000d0 ]
    geom.iniP( 75).pos(1:3) = [   47.7733d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 76).pos(1:3) = [   61.2144d0,  -81.9198d0, 0.0000d0 ]
    geom.iniP( 77).pos(1:3) = [   50.0518d0,  -57.1298d0, 0.0000d0 ]
    geom.iniP( 78).pos(1:3) = [   65.0537d0,  -32.1778d0, 0.0000d0 ]
    geom.iniP( 79).pos(1:3) = [   50.8326d0,   -7.4510d0, 0.0000d0 ]
    geom.iniP( 80).pos(1:3) = [   67.4113d0,   17.8063d0, 0.0000d0 ]
    geom.iniP( 81).pos(1:3) = [   72.5103d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 82).pos(1:3) = [   85.8801d0,  -82.6271d0, 0.0000d0 ]
    geom.iniP( 83).pos(1:3) = [   75.3803d0,  -56.6101d0, 0.0000d0 ]
    geom.iniP( 84).pos(1:3) = [   92.8008d0,  -33.3313d0, 0.0000d0 ]
    geom.iniP( 85).pos(1:3) = [   76.4936d0,   -8.1051d0, 0.0000d0 ]
    geom.iniP( 86).pos(1:3) = [   91.9000d0,   17.8063d0, 0.0000d0 ]
    geom.iniP( 87).pos(1:3) = [   97.3376d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 88).pos(1:3) = [  111.6736d0,  -86.2066d0, 0.0000d0 ]
    geom.iniP( 89).pos(1:3) = [  100.0002d0,  -60.1166d0, 0.0000d0 ]
    geom.iniP( 90).pos(1:3) = [  120.3833d0,  -45.9497d0, 0.0000d0 ]
    geom.iniP( 91).pos(1:3) = [  100.3391d0,   -6.8079d0, 0.0000d0 ]
    geom.iniP( 92).pos(1:3) = [  117.4091d0,   17.8063d0, 0.0000d0 ]
    geom.iniP( 93).pos(1:3) = [  123.7122d0, -107.3907d0, 0.0000d0 ]
    geom.iniP( 94).pos(1:3) = [  146.3091d0,  -84.4810d0, 0.0000d0 ]
    geom.iniP( 95).pos(1:3) = [  125.8913d0,  -70.5457d0, 0.0000d0 ]
    geom.iniP( 96).pos(1:3) = [  146.3091d0,  -33.3054d0, 0.0000d0 ]
    geom.iniP( 97).pos(1:3) = [  119.9165d0,  -24.1763d0, 0.0000d0 ]
    geom.iniP( 98).pos(1:3) = [  125.8604d0,   -1.2392d0, 0.0000d0 ]
    geom.iniP( 99).pos(1:3) = [  146.3091d0,  -57.5196d0, 0.0000d0 ]
    geom.iniP(100).pos(1:3) = [  146.3091d0,   -9.1617d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face(  1).n_poi = 3; allocate(geom.face(  1).poi(3)); geom.face(  1).poi(1:3) = [  15,  28,  16 ]
    geom.face(  2).n_poi = 3; allocate(geom.face(  2).poi(3)); geom.face(  2).poi(1:3) = [  24,  13,  12 ]
    geom.face(  3).n_poi = 3; allocate(geom.face(  3).poi(3)); geom.face(  3).poi(1:3) = [  24,  23,  36 ]
    geom.face(  4).n_poi = 3; allocate(geom.face(  4).poi(3)); geom.face(  4).poi(1:3) = [  28,  29,  17 ]
    geom.face(  5).n_poi = 3; allocate(geom.face(  5).poi(3)); geom.face(  5).poi(1:3) = [   6,  16,  17 ]
    geom.face(  6).n_poi = 3; allocate(geom.face(  6).poi(3)); geom.face(  6).poi(1:3) = [  17,  16,  28 ]
    geom.face(  7).n_poi = 3; allocate(geom.face(  7).poi(3)); geom.face(  7).poi(1:3) = [  27,  28,  15 ]
    geom.face(  8).n_poi = 3; allocate(geom.face(  8).poi(3)); geom.face(  8).poi(1:3) = [  15,  26,  27 ]
    geom.face(  9).n_poi = 3; allocate(geom.face(  9).poi(3)); geom.face(  9).poi(1:3) = [  14,  26,  15 ]
    geom.face( 10).n_poi = 3; allocate(geom.face( 10).poi(3)); geom.face( 10).poi(1:3) = [  13,  26,  14 ]
    geom.face( 11).n_poi = 3; allocate(geom.face( 11).poi(3)); geom.face( 11).poi(1:3) = [  99,  96,  90 ]
    geom.face( 12).n_poi = 3; allocate(geom.face( 12).poi(3)); geom.face( 12).poi(1:3) = [  84,  90,  97 ]
    geom.face( 13).n_poi = 3; allocate(geom.face( 13).poi(3)); geom.face( 13).poi(1:3) = [  97,  96, 100 ]
    geom.face( 14).n_poi = 3; allocate(geom.face( 14).poi(3)); geom.face( 14).poi(1:3) = [  97,  90,  96 ]
    geom.face( 15).n_poi = 3; allocate(geom.face( 15).poi(3)); geom.face( 15).poi(1:3) = [  21,   9,  20 ]
    geom.face( 16).n_poi = 3; allocate(geom.face( 16).poi(3)); geom.face( 16).poi(1:3) = [  20,  32,  21 ]
    geom.face( 17).n_poi = 3; allocate(geom.face( 17).poi(3)); geom.face( 17).poi(1:3) = [  11,  24,  12 ]
    geom.face( 18).n_poi = 3; allocate(geom.face( 18).poi(3)); geom.face( 18).poi(1:3) = [  23,  24,  11 ]
    geom.face( 19).n_poi = 3; allocate(geom.face( 19).poi(3)); geom.face( 19).poi(1:3) = [  22,   9,  21 ]
    geom.face( 20).n_poi = 3; allocate(geom.face( 20).poi(3)); geom.face( 20).poi(1:3) = [  22,  10,   9 ]
    geom.face( 21).n_poi = 3; allocate(geom.face( 21).poi(3)); geom.face( 21).poi(1:3) = [  23,  11,  22 ]
    geom.face( 22).n_poi = 3; allocate(geom.face( 22).poi(3)); geom.face( 22).poi(1:3) = [  22,  11,  10 ]
    geom.face( 23).n_poi = 3; allocate(geom.face( 23).poi(3)); geom.face( 23).poi(1:3) = [  25,  26,  13 ]
    geom.face( 24).n_poi = 3; allocate(geom.face( 24).poi(3)); geom.face( 24).poi(1:3) = [  13,  24,  25 ]
    geom.face( 25).n_poi = 3; allocate(geom.face( 25).poi(3)); geom.face( 25).poi(1:3) = [  25,  24,  36 ]
    geom.face( 26).n_poi = 3; allocate(geom.face( 26).poi(3)); geom.face( 26).poi(1:3) = [  36,  23,  35 ]
    geom.face( 27).n_poi = 3; allocate(geom.face( 27).poi(3)); geom.face( 27).poi(1:3) = [  21,  32,  33 ]
    geom.face( 28).n_poi = 3; allocate(geom.face( 28).poi(3)); geom.face( 28).poi(1:3) = [  28,  27,  40 ]
    geom.face( 29).n_poi = 3; allocate(geom.face( 29).poi(3)); geom.face( 29).poi(1:3) = [  41,  29,  40 ]
    geom.face( 30).n_poi = 3; allocate(geom.face( 30).poi(3)); geom.face( 30).poi(1:3) = [  40,  29,  28 ]
    geom.face( 31).n_poi = 3; allocate(geom.face( 31).poi(3)); geom.face( 31).poi(1:3) = [  98,  97, 100 ]
    geom.face( 32).n_poi = 3; allocate(geom.face( 32).poi(3)); geom.face( 32).poi(1:3) = [ 100,   3,  98 ]
    geom.face( 33).n_poi = 3; allocate(geom.face( 33).poi(3)); geom.face( 33).poi(1:3) = [  98,   3,  92 ]
    geom.face( 34).n_poi = 3; allocate(geom.face( 34).poi(3)); geom.face( 34).poi(1:3) = [  91,  98,  92 ]
    geom.face( 35).n_poi = 3; allocate(geom.face( 35).poi(3)); geom.face( 35).poi(1:3) = [  84,  97,  91 ]
    geom.face( 36).n_poi = 3; allocate(geom.face( 36).poi(3)); geom.face( 36).poi(1:3) = [  97,  98,  91 ]
    geom.face( 37).n_poi = 3; allocate(geom.face( 37).poi(3)); geom.face( 37).poi(1:3) = [  80,  71,  79 ]
    geom.face( 38).n_poi = 3; allocate(geom.face( 38).poi(3)); geom.face( 38).poi(1:3) = [  67,  76,  77 ]
    geom.face( 39).n_poi = 3; allocate(geom.face( 39).poi(3)); geom.face( 39).poi(1:3) = [  42,  31,  30 ]
    geom.face( 40).n_poi = 3; allocate(geom.face( 40).poi(3)); geom.face( 40).poi(1:3) = [   8,  20,   9 ]
    geom.face( 41).n_poi = 3; allocate(geom.face( 41).poi(3)); geom.face( 41).poi(1:3) = [   1,  18,   7 ]
    geom.face( 42).n_poi = 3; allocate(geom.face( 42).poi(3)); geom.face( 42).poi(1:3) = [  20,   8,   7 ]
    geom.face( 43).n_poi = 3; allocate(geom.face( 43).poi(3)); geom.face( 43).poi(1:3) = [  47,  35,  46 ]
    geom.face( 44).n_poi = 3; allocate(geom.face( 44).poi(3)); geom.face( 44).poi(1:3) = [  46,  58,  47 ]
    geom.face( 45).n_poi = 3; allocate(geom.face( 45).poi(3)); geom.face( 45).poi(1:3) = [  23,  22,  34 ]
    geom.face( 46).n_poi = 3; allocate(geom.face( 46).poi(3)); geom.face( 46).poi(1:3) = [  34,  35,  23 ]
    geom.face( 47).n_poi = 3; allocate(geom.face( 47).poi(3)); geom.face( 47).poi(1:3) = [  34,  22,  21 ]
    geom.face( 48).n_poi = 3; allocate(geom.face( 48).poi(3)); geom.face( 48).poi(1:3) = [  21,  33,  34 ]
    geom.face( 49).n_poi = 3; allocate(geom.face( 49).poi(3)); geom.face( 49).poi(1:3) = [  46,  35,  34 ]
    geom.face( 50).n_poi = 3; allocate(geom.face( 50).poi(3)); geom.face( 50).poi(1:3) = [  34,  33,  46 ]
    geom.face( 51).n_poi = 3; allocate(geom.face( 51).poi(3)); geom.face( 51).poi(1:3) = [  46,  33,  45 ]
    geom.face( 52).n_poi = 3; allocate(geom.face( 52).poi(3)); geom.face( 52).poi(1:3) = [  45,  58,  46 ]
    geom.face( 53).n_poi = 3; allocate(geom.face( 53).poi(3)); geom.face( 53).poi(1:3) = [  45,  56,  57 ]
    geom.face( 54).n_poi = 3; allocate(geom.face( 54).poi(3)); geom.face( 54).poi(1:3) = [  57,  58,  45 ]
    geom.face( 55).n_poi = 3; allocate(geom.face( 55).poi(3)); geom.face( 55).poi(1:3) = [  65,  53,  64 ]
    geom.face( 56).n_poi = 3; allocate(geom.face( 56).poi(3)); geom.face( 56).poi(1:3) = [  64,  63,  74 ]
    geom.face( 57).n_poi = 3; allocate(geom.face( 57).poi(3)); geom.face( 57).poi(1:3) = [   5,  65,  74 ]
    geom.face( 58).n_poi = 3; allocate(geom.face( 58).poi(3)); geom.face( 58).poi(1:3) = [  74,  65,  64 ]
    geom.face( 59).n_poi = 3; allocate(geom.face( 59).poi(3)); geom.face( 59).poi(1:3) = [  52,  53,  41 ]
    geom.face( 60).n_poi = 3; allocate(geom.face( 60).poi(3)); geom.face( 60).poi(1:3) = [  41,  40,  52 ]
    geom.face( 61).n_poi = 3; allocate(geom.face( 61).poi(3)); geom.face( 61).poi(1:3) = [  64,  53,  52 ]
    geom.face( 62).n_poi = 3; allocate(geom.face( 62).poi(3)); geom.face( 62).poi(1:3) = [  61,  49,  60 ]
    geom.face( 63).n_poi = 3; allocate(geom.face( 63).poi(3)); geom.face( 63).poi(1:3) = [  26,  25,  38 ]
    geom.face( 64).n_poi = 3; allocate(geom.face( 64).poi(3)); geom.face( 64).poi(1:3) = [  38,  27,  26 ]
    geom.face( 65).n_poi = 3; allocate(geom.face( 65).poi(3)); geom.face( 65).poi(1:3) = [  86,  91,  92 ]
    geom.face( 66).n_poi = 3; allocate(geom.face( 66).poi(3)); geom.face( 66).poi(1:3) = [  60,  47,  59 ]
    geom.face( 67).n_poi = 3; allocate(geom.face( 67).poi(3)); geom.face( 67).poi(1:3) = [  59,  47,  58 ]
    geom.face( 68).n_poi = 3; allocate(geom.face( 68).poi(3)); geom.face( 68).poi(1:3) = [  73,  61,  72 ]
    geom.face( 69).n_poi = 3; allocate(geom.face( 69).poi(3)); geom.face( 69).poi(1:3) = [  72,  61,  60 ]
    geom.face( 70).n_poi = 3; allocate(geom.face( 70).poi(3)); geom.face( 70).poi(1:3) = [  80,  79,  85 ]
    geom.face( 71).n_poi = 3; allocate(geom.face( 71).poi(3)); geom.face( 71).poi(1:3) = [  85,  86,  80 ]
    geom.face( 72).n_poi = 3; allocate(geom.face( 72).poi(3)); geom.face( 72).poi(1:3) = [  91,  86,  85 ]
    geom.face( 73).n_poi = 3; allocate(geom.face( 73).poi(3)); geom.face( 73).poi(1:3) = [  79,  78,  85 ]
    geom.face( 74).n_poi = 3; allocate(geom.face( 74).poi(3)); geom.face( 74).poi(1:3) = [  84,  91,  85 ]
    geom.face( 75).n_poi = 3; allocate(geom.face( 75).poi(3)); geom.face( 75).poi(1:3) = [  85,  78,  84 ]
    geom.face( 76).n_poi = 3; allocate(geom.face( 76).poi(3)); geom.face( 76).poi(1:3) = [  94,  93,   2 ]
    geom.face( 77).n_poi = 3; allocate(geom.face( 77).poi(3)); geom.face( 77).poi(1:3) = [  75,  67,  66 ]
    geom.face( 78).n_poi = 3; allocate(geom.face( 78).poi(3)); geom.face( 78).poi(1:3) = [  76,  67,  75 ]
    geom.face( 79).n_poi = 3; allocate(geom.face( 79).poi(3)); geom.face( 79).poi(1:3) = [  89,  90,  84 ]
    geom.face( 80).n_poi = 3; allocate(geom.face( 80).poi(3)); geom.face( 80).poi(1:3) = [  68,  67,  77 ]
    geom.face( 81).n_poi = 3; allocate(geom.face( 81).poi(3)); geom.face( 81).poi(1:3) = [  57,  56,  68 ]
    geom.face( 82).n_poi = 3; allocate(geom.face( 82).poi(3)); geom.face( 82).poi(1:3) = [  19,   7,  18 ]
    geom.face( 83).n_poi = 3; allocate(geom.face( 83).poi(3)); geom.face( 83).poi(1:3) = [  20,   7,  19 ]
    geom.face( 84).n_poi = 3; allocate(geom.face( 84).poi(3)); geom.face( 84).poi(1:3) = [  18,  30,  19 ]
    geom.face( 85).n_poi = 3; allocate(geom.face( 85).poi(3)); geom.face( 85).poi(1:3) = [  19,  30,  31 ]
    geom.face( 86).n_poi = 3; allocate(geom.face( 86).poi(3)); geom.face( 86).poi(1:3) = [  32,  20,  19 ]
    geom.face( 87).n_poi = 3; allocate(geom.face( 87).poi(3)); geom.face( 87).poi(1:3) = [  19,  31,  32 ]
    geom.face( 88).n_poi = 3; allocate(geom.face( 88).poi(3)); geom.face( 88).poi(1:3) = [  43,  31,  42 ]
    geom.face( 89).n_poi = 3; allocate(geom.face( 89).poi(3)); geom.face( 89).poi(1:3) = [  42,  54,  43 ]
    geom.face( 90).n_poi = 3; allocate(geom.face( 90).poi(3)); geom.face( 90).poi(1:3) = [  60,  49,  48 ]
    geom.face( 91).n_poi = 3; allocate(geom.face( 91).poi(3)); geom.face( 91).poi(1:3) = [  48,  47,  60 ]
    geom.face( 92).n_poi = 3; allocate(geom.face( 92).poi(3)); geom.face( 92).poi(1:3) = [  36,  35,  48 ]
    geom.face( 93).n_poi = 3; allocate(geom.face( 93).poi(3)); geom.face( 93).poi(1:3) = [  35,  47,  48 ]
    geom.face( 94).n_poi = 3; allocate(geom.face( 94).poi(3)); geom.face( 94).poi(1:3) = [  73,  63,  62 ]
    geom.face( 95).n_poi = 3; allocate(geom.face( 95).poi(3)); geom.face( 95).poi(1:3) = [  62,  61,  73 ]
    geom.face( 96).n_poi = 3; allocate(geom.face( 96).poi(3)); geom.face( 96).poi(1:3) = [  62,  63,  64 ]
    geom.face( 97).n_poi = 3; allocate(geom.face( 97).poi(3)); geom.face( 97).poi(1:3) = [  37,  25,  36 ]
    geom.face( 98).n_poi = 3; allocate(geom.face( 98).poi(3)); geom.face( 98).poi(1:3) = [  37,  38,  25 ]
    geom.face( 99).n_poi = 3; allocate(geom.face( 99).poi(3)); geom.face( 99).poi(1:3) = [  36,  48,  37 ]
    geom.face(100).n_poi = 3; allocate(geom.face(100).poi(3)); geom.face(100).poi(1:3) = [  37,  48,  49 ]
    geom.face(101).n_poi = 3; allocate(geom.face(101).poi(3)); geom.face(101).poi(1:3) = [  39,  40,  27 ]
    geom.face(102).n_poi = 3; allocate(geom.face(102).poi(3)); geom.face(102).poi(1:3) = [  27,  38,  39 ]
    geom.face(103).n_poi = 3; allocate(geom.face(103).poi(3)); geom.face(103).poi(1:3) = [  39,  52,  40 ]
    geom.face(104).n_poi = 3; allocate(geom.face(104).poi(3)); geom.face(104).poi(1:3) = [  60,  59,   4 ]
    geom.face(105).n_poi = 3; allocate(geom.face(105).poi(3)); geom.face(105).poi(1:3) = [   4,  72,  60 ]
    geom.face(106).n_poi = 3; allocate(geom.face(106).poi(3)); geom.face(106).poi(1:3) = [  77,  78,  69 ]
    geom.face(107).n_poi = 3; allocate(geom.face(107).poi(3)); geom.face(107).poi(1:3) = [  57,  68,  69 ]
    geom.face(108).n_poi = 3; allocate(geom.face(108).poi(3)); geom.face(108).poi(1:3) = [  69,  68,  77 ]
    geom.face(109).n_poi = 3; allocate(geom.face(109).poi(3)); geom.face(109).poi(1:3) = [  69,  78,  79 ]
    geom.face(110).n_poi = 3; allocate(geom.face(110).poi(3)); geom.face(110).poi(1:3) = [  76,  75,  81 ]
    geom.face(111).n_poi = 3; allocate(geom.face(111).poi(3)); geom.face(111).poi(1:3) = [  77,  76,  83 ]
    geom.face(112).n_poi = 3; allocate(geom.face(112).poi(3)); geom.face(112).poi(1:3) = [  83,  78,  77 ]
    geom.face(113).n_poi = 3; allocate(geom.face(113).poi(3)); geom.face(113).poi(1:3) = [  84,  78,  83 ]
    geom.face(114).n_poi = 3; allocate(geom.face(114).poi(3)); geom.face(114).poi(1:3) = [  83,  89,  84 ]
    geom.face(115).n_poi = 3; allocate(geom.face(115).poi(3)); geom.face(115).poi(1:3) = [  82,  81,  87 ]
    geom.face(116).n_poi = 3; allocate(geom.face(116).poi(3)); geom.face(116).poi(1:3) = [  76,  81,  82 ]
    geom.face(117).n_poi = 3; allocate(geom.face(117).poi(3)); geom.face(117).poi(1:3) = [  82,  83,  76 ]
    geom.face(118).n_poi = 3; allocate(geom.face(118).poi(3)); geom.face(118).poi(1:3) = [  89,  83,  82 ]
    geom.face(119).n_poi = 3; allocate(geom.face(119).poi(3)); geom.face(119).poi(1:3) = [  93,  94,  88 ]
    geom.face(120).n_poi = 3; allocate(geom.face(120).poi(3)); geom.face(120).poi(1:3) = [  89,  82,  88 ]
    geom.face(121).n_poi = 3; allocate(geom.face(121).poi(3)); geom.face(121).poi(1:3) = [  87,  93,  88 ]
    geom.face(122).n_poi = 3; allocate(geom.face(122).poi(3)); geom.face(122).poi(1:3) = [  88,  82,  87 ]
    geom.face(123).n_poi = 3; allocate(geom.face(123).poi(3)); geom.face(123).poi(1:3) = [  55,  54,  66 ]
    geom.face(124).n_poi = 3; allocate(geom.face(124).poi(3)); geom.face(124).poi(1:3) = [  55,  43,  54 ]
    geom.face(125).n_poi = 3; allocate(geom.face(125).poi(3)); geom.face(125).poi(1:3) = [  66,  67,  55 ]
    geom.face(126).n_poi = 3; allocate(geom.face(126).poi(3)); geom.face(126).poi(1:3) = [  56,  43,  55 ]
    geom.face(127).n_poi = 3; allocate(geom.face(127).poi(3)); geom.face(127).poi(1:3) = [  55,  68,  56 ]
    geom.face(128).n_poi = 3; allocate(geom.face(128).poi(3)); geom.face(128).poi(1:3) = [  67,  68,  55 ]
    geom.face(129).n_poi = 3; allocate(geom.face(129).poi(3)); geom.face(129).poi(1:3) = [  32,  31,  44 ]
    geom.face(130).n_poi = 3; allocate(geom.face(130).poi(3)); geom.face(130).poi(1:3) = [  31,  43,  44 ]
    geom.face(131).n_poi = 3; allocate(geom.face(131).poi(3)); geom.face(131).poi(1:3) = [  44,  33,  32 ]
    geom.face(132).n_poi = 3; allocate(geom.face(132).poi(3)); geom.face(132).poi(1:3) = [  44,  45,  33 ]
    geom.face(133).n_poi = 3; allocate(geom.face(133).poi(3)); geom.face(133).poi(1:3) = [  56,  45,  44 ]
    geom.face(134).n_poi = 3; allocate(geom.face(134).poi(3)); geom.face(134).poi(1:3) = [  44,  43,  56 ]
    geom.face(135).n_poi = 3; allocate(geom.face(135).poi(3)); geom.face(135).poi(1:3) = [  52,  39,  51 ]
    geom.face(136).n_poi = 3; allocate(geom.face(136).poi(3)); geom.face(136).poi(1:3) = [  64,  52,  51 ]
    geom.face(137).n_poi = 3; allocate(geom.face(137).poi(3)); geom.face(137).poi(1:3) = [  51,  62,  64 ]
    geom.face(138).n_poi = 3; allocate(geom.face(138).poi(3)); geom.face(138).poi(1:3) = [  70,  69,  79 ]
    geom.face(139).n_poi = 3; allocate(geom.face(139).poi(3)); geom.face(139).poi(1:3) = [  70,   4,  59 ]
    geom.face(140).n_poi = 3; allocate(geom.face(140).poi(3)); geom.face(140).poi(1:3) = [  57,  69,  70 ]
    geom.face(141).n_poi = 3; allocate(geom.face(141).poi(3)); geom.face(141).poi(1:3) = [  70,  79,  71 ]
    geom.face(142).n_poi = 3; allocate(geom.face(142).poi(3)); geom.face(142).poi(1:3) = [  71,   4,  70 ]
    geom.face(143).n_poi = 3; allocate(geom.face(143).poi(3)); geom.face(143).poi(1:3) = [  70,  58,  57 ]
    geom.face(144).n_poi = 3; allocate(geom.face(144).poi(3)); geom.face(144).poi(1:3) = [  70,  59,  58 ]
    geom.face(145).n_poi = 3; allocate(geom.face(145).poi(3)); geom.face(145).poi(1:3) = [  90,  89,  95 ]
    geom.face(146).n_poi = 3; allocate(geom.face(146).poi(3)); geom.face(146).poi(1:3) = [  89,  88,  95 ]
    geom.face(147).n_poi = 3; allocate(geom.face(147).poi(3)); geom.face(147).poi(1:3) = [  99,  90,  95 ]
    geom.face(148).n_poi = 3; allocate(geom.face(148).poi(3)); geom.face(148).poi(1:3) = [  95,  94,  99 ]
    geom.face(149).n_poi = 3; allocate(geom.face(149).poi(3)); geom.face(149).poi(1:3) = [  95,  88,  94 ]
    geom.face(150).n_poi = 3; allocate(geom.face(150).poi(3)); geom.face(150).poi(1:3) = [  62,  51,  50 ]
    geom.face(151).n_poi = 3; allocate(geom.face(151).poi(3)); geom.face(151).poi(1:3) = [  49,  61,  50 ]
    geom.face(152).n_poi = 3; allocate(geom.face(152).poi(3)); geom.face(152).poi(1:3) = [  61,  62,  50 ]
    geom.face(153).n_poi = 3; allocate(geom.face(153).poi(3)); geom.face(153).poi(1:3) = [  50,  37,  49 ]
    geom.face(154).n_poi = 3; allocate(geom.face(154).poi(3)); geom.face(154).poi(1:3) = [  38,  37,  50 ]
    geom.face(155).n_poi = 3; allocate(geom.face(155).poi(3)); geom.face(155).poi(1:3) = [  50,  39,  38 ]
    geom.face(156).n_poi = 3; allocate(geom.face(156).poi(3)); geom.face(156).poi(1:3) = [  50,  51,  39 ]
end subroutine Exam_Open2D_L_Shape_Irregular

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

    prob.name_file = "Plate_Distorted_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Quad"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    prob.name_file = "Plate_Distorted_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Tri"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

! Example of hyperbolic paraboloid with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Hyperbolic_Paraboloid_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width
    integer :: i, j, index, n, nx, ny
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Hyperbolic_Paraboloid_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hyperbolic Paraboloid Quad"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    n  = 3
    nx = n
    ny = n

    geom.n_iniP = (nx + 1) * (ny + 1)
    geom.n_face = nx * ny

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    index = 0
    do i = 1, nx + 1
        do j = 1, ny + 1
            index = index + 1
            geom.iniP(index).pos(1) = (dble(i) - 1.0d0) * (1.0d0 / dble(nx)) - 0.5d0
            geom.iniP(index).pos(2) = (dble(j) - 1.0d0) * (1.0d0 / dble(ny)) - 0.5d0
            geom.iniP(index).pos(3) = geom.iniP(index).pos(1)**2.0d0 - geom.iniP(index).pos(2)**2.0d0
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nx
        do j = 1, ny
            index = index + 1
            geom.face(index).poi(1) = (ny+1)*(i+0) + (j-1) + 1
            geom.face(index).poi(2) = (ny+1)*(i+0) + (j-1) + 2
            geom.face(index).poi(3) = (ny+1)*(i-1) + (j-1) + 2
            geom.face(index).poi(4) = (ny+1)*(i-1) + (j-1) + 1
        end do
    end do
end subroutine Exam_Open2D_Hyperbolic_Paraboloid_Quad

! ---------------------------------------------------------------------------------------

! Example of hyperbolic paraboloid with quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Hyperbolic_Paraboloid_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, j, index, n, nx, ny
    character :: pn
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "Hyperbolic_Paraboloid_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hyperbolic Paraboloid Tri"

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [52, 152, 219], "XY", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    n  = 3
    nx = n
    ny = n
    pn = "\"

    geom.n_iniP = (nx + 1) * (ny + 1)
    geom.n_face = nx * ny * 2

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    index = 0
    do i = 1, nx + 1
        do j = 1, ny + 1
            index = index + 1
            geom.iniP(index).pos(1) = (dble(i) - 1.0d0) * (1.0d0 / dble(nx)) - 0.5d0
            geom.iniP(index).pos(2) = (dble(j) - 1.0d0) * (1.0d0 / dble(ny)) - 0.5d0
            geom.iniP(index).pos(3) = geom.iniP(index).pos(1)**2.0d0 - geom.iniP(index).pos(2)**2.0d0
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nx
        do j = 1, ny
            if(pn == "\") then
                index = index + 1
                geom.face(index).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(index).poi(2) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(index).poi(3) = (ny+1)*(i-1) + (j-1) + 1

                index = index + 1
                geom.face(index).poi(1) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(index).poi(2) = (ny+1)*(i-1) + (j-1) + 2
                geom.face(index).poi(3) = (ny+1)*(i-1) + (j-1) + 1
            else if(pn == "/") then
                index = index + 1
                geom.face(index).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(index).poi(2) = (ny+1)*(i-1) + (j-1) + 2
                geom.face(index).poi(3) = (ny+1)*(i-1) + (j-1) + 1

                index = index + 1
                geom.face(index).poi(1) = (ny+1)*(i+0) + (j-1) + 1
                geom.face(index).poi(2) = (ny+1)*(i+0) + (j-1) + 2
                geom.face(index).poi(3) = (ny+1)*(i-1) + (j-1) + 2
            end if
        end do
    end do
end subroutine Exam_Open2D_Hyperbolic_Paraboloid_Tri

! ---------------------------------------------------------------------------------------

! Return cross point
! Last updated on Sat 11 Mar 2017 by Hyungmin
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

! Merge points and faces for quad mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Merge_Point_Face_Quad(n, joint, conn, n_node, n_element, geom) 
    type(GeomType),   intent(inout) :: geom
    double precision, intent(in)    :: joint((n+1)**2, 3)
    integer, intent(in)    :: n, conn(n*n, 4)
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

! Merge points and faces for tri mesh
! Last updated on Sat 11 Mar 2017 by Hyungmin
subroutine Exam_Open2D_Merge_Point_Face_Tri(n, joint, conn, n_node, n_element, geom)  
    type(GeomType), intent(inout) :: geom
    double precision, intent(in)  :: joint((n+1)**2, 3)
    integer, intent(in)    :: n, conn(n*n*2, 3)
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

! Return 0 or 1
! Last updated on Sat 11 Mar 2017 by Hyungmin
integer function Exam_Open2D_Comp_XYZ(a, b)
   double precision, intent(in) :: a(3), b(3) 
   double precision :: d

   Exam_Open2D_Comp_XYZ = 0
   d = dsqrt((a(1)-b(1))**2.0d0 + (a(2)-b(2))**2.0d0 + (a(3)-b(3))**2.0d0)
   if(d < 1.0d-5) Exam_Open2D_Comp_XYZ = 1
end function Exam_Open2D_Comp_XYZ

! ---------------------------------------------------------------------------------------

end module Exam_2D_Open