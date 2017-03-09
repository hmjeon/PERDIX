!
! ---------------------------------------------------------------------------------------
!
!                               Module for Exam_OpenGeo
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/10/20
!                                                            Last  modified  : 2016/08/30
!
! ---------------------------------------------------------------------------------------
!
module Exam_OpenGeo

    use Data_Prob
    use Data_Geom

    use Para
    use Mani
    use Math

    implicit none

    public  Exam_OpenGeo_Plate_Uniform_Quad             ! 53. Plate with uniform quadrilater mesh
    public  Exam_OpenGeo_Plate_Distorted_Quad           ! 54. Plate with distorted quadrilater mesh
    public  Exam_OpenGeo_Plate_Uniform_Tri              ! 55. Plate with uniform triangular mesh
    public  Exam_OpenGeo_Plate_Distorted_Tri            ! 56. Plate with distorted triangular mesh
    public  Exam_OpenGeo_Circular_Plate_Quad            ! 57. Circular plate with quadrilater mesh
    public  Exam_OpenGeo_Circular_Plate_Tri             ! 58. Circular plate with triangular mesh
    public  Exam_OpenGeo_Annular_Plate_Quad             ! 59. Annular plate with quadrilater mesh
    public  Exam_OpenGeo_Annular_Plate_Tri              ! 60. Annular plate with triangular mesh
    public  Exam_OpenGeo_Hyperbolic_Paraboloid_Quad     ! 61. Hyperbolic paraboloid with quadrilater mesh
    public  Exam_OpenGeo_Hyperbolic_Paraboloid_Tri      ! 62. Hyperbolic paraboloid with triangular mesh

    private Exam_OpenGeo_Cross_Point
    private Exam_OpenGeo_Merge_Point_Face_Quad
    private Exam_OpenGeo_Merge_Point_Face_Tri
    private Exam_OpenGeo_Comp_XYZ

    contains

! ---------------------------------------------------------------------------------------

! Example of retagular plate geometry with the uniform mesh of quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Plate_Uniform_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: n_i_point, n_j_point, n_i_face, n_j_face
    integer :: i, j, n, numbering
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
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n         = 2
    x_width   = 1.0d0           ! x length
    y_width   = 1.0d0           ! y length
    n_i_point = n + 1           ! # of points in x-direction
    n_j_point = n + 1           ! # of points in y-direction
    n_i_face  = n_i_point - 1
    n_j_face  = n_j_point - 1
    del_x     = x_width / dble(n_i_point - 1)
    del_y     = y_width / dble(n_j_point - 1)

    geom.n_iniP = n_i_point * n_j_point
    geom.n_face = n_i_face  * n_j_face

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    
    ! Set position vector
    do j = 1, n_j_point
        do i = 1, n_i_point
            numbering = n_i_point * (j - 1) + i

            geom.iniP(numbering).pos(1) = del_x * dble(i - 1)
            geom.iniP(numbering).pos(2) = del_y * dble(j - 1)
            geom.iniP(numbering).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, n_j_face
        do i = 1, n_i_face
            numbering = n_i_face * (j - 1) + i

            geom.face(numbering).n_poi = 4
            allocate(geom.face(numbering).poi(4))

            geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
            geom.face(numbering).poi(2) = n_i_point * (j - 1) + i + 1
            geom.face(numbering).poi(3) = n_i_point * j + i + 1
            geom.face(numbering).poi(4) = n_i_point * j + i
        end do
    end do
end subroutine Exam_OpenGeo_Plate_Uniform_Quad

! ---------------------------------------------------------------------------------------

! Example of retagular plate geometry with the distorted mesh of quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Plate_Distorted_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision, allocatable :: edge(:,:,:)
    double precision :: ff, dd, xy(2)
    double precision :: x_width, y_width, del_x, del_y
    integer :: n_i_point, n_j_point, n_i_face, n_j_face
    integer :: i, j, n, distort, numbering
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "54_Plate_Distorted_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n         = 2
    distort   = 1
    x_width   = 1.0d0           ! x length
    y_width   = 1.0d0           ! y length
    n_i_point = n + 1           ! # of points in x-direction
    n_j_point = n + 1           ! # of points in y-direction
    n_i_face  = n_i_point - 1
    n_j_face  = n_j_point - 1
    del_x     = x_width / dble(n_i_point - 1)
    del_y     = y_width / dble(n_j_point - 1)

    geom.n_iniP = n_i_point * n_j_point
    geom.n_face = n_i_face  * n_j_face

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    allocate(edge(4, n + 1, 2))

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
    do j = 1, n_j_point
        do i = 1, n_i_point
            numbering = n_i_point * (j - 1) + i

            call Exam_OpenGeo_Cross_Point(edge(1,i,:), edge(2,i,:), edge(3,j,:), edge(4,j,:), xy)

            geom.iniP(numbering).pos(1) = xy(1)
            geom.iniP(numbering).pos(2) = xy(2)
            geom.iniP(numbering).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, n_j_face
        do i = 1, n_i_face
            numbering = n_i_face * (j - 1) + i

            geom.face(numbering).n_poi = 4
            allocate(geom.face(numbering).poi(4))

            geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
            geom.face(numbering).poi(2) = n_i_point * (j - 1) + i + 1
            geom.face(numbering).poi(3) = n_i_point * j + i + 1
            geom.face(numbering).poi(4) = n_i_point * j + i
        end do
    end do

    ! Deallocate edge data
    deallocate(edge)
end subroutine Exam_OpenGeo_Plate_Distorted_Quad

! ---------------------------------------------------------------------------------------

! Example of plate with the uniform mesh of triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Plate_Uniform_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: n_i_point, n_j_point, n_i_face, n_j_face
    integer :: i, j, n, numbering
    character :: pattern
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "55_Plate_Uniform_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "55_Plate Uniform Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n         = 2
    pattern   = "\"
    x_width   = 1.0d0           ! x length
    y_width   = 1.0d0           ! y length
    n_i_point = n + 1           ! # of points in x-direction
    n_j_point = n + 1           ! # of points in y-direction
    n_i_face  = n_i_point - 1
    n_j_face  = n_j_point - 1
    del_x     = x_width / dble(n_i_point - 1)
    del_y     = y_width / dble(n_j_point - 1)

    geom.n_iniP = n_i_point * n_j_point
    geom.n_face = n_i_face  * n_j_face * 2
    
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vector
    do j = 1, n_j_point
        do i = 1, n_i_point
            numbering = n_i_point * (j - 1) + i

            geom.iniP(numbering).pos(1) = del_x * dble(i - 1)
            geom.iniP(numbering).pos(2) = del_y * dble(j - 1)
            geom.iniP(numbering).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, n_j_face
        do i = 1, n_i_face
            if(pattern == "\") then             ! Mesh pattern, \
                numbering = 2*(n_i_face * (j-1) + i) - 1

                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j-1) + i
                geom.face(numbering).poi(2) = n_i_point * (j-1) + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i

                numbering = 2*(n_i_face * (j-1) + i)
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j-1) + i + 1
                geom.face(numbering).poi(2) = n_i_point * j + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i
            else if(pattern == "/") then        ! Mesh pattern, /
                numbering = 2*(n_i_face * (j-1) + i) - 1
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
                geom.face(numbering).poi(2) = n_i_point * j + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i

                numbering = 2*(n_i_face * (j-1) + i)
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
                geom.face(numbering).poi(2) = n_i_point * (j - 1) + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i + 1
            end if
        end do
    end do
end subroutine Exam_OpenGeo_Plate_Uniform_Tri

! ---------------------------------------------------------------------------------------

! Example of plate with the distorted mesh of triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Plate_Distorted_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision, allocatable :: edge(:,:,:)
    double precision :: ff, dd, xy(2)
    double precision :: x_width, y_width, del_x, del_y
    integer :: n_i_point, n_j_point, n_i_face, n_j_face
    integer :: i, j, n, distort, numbering 
    character :: pattern
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "56_Plate_Distorted_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Distorted Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n         = 2
    pattern   = "\"
    distort   = 3
    x_width   = 1.0d0           ! x length
    y_width   = 1.0d0           ! y length
    n_i_point = n + 1           ! # of points in x-direction
    n_j_point = n + 1           ! # of points in y-direction
    n_i_face  = n_i_point - 1
    n_j_face  = n_j_point - 1
    del_x     = x_width / dble(n_i_point - 1)
    del_y     = y_width / dble(n_j_point - 1)

    geom.n_iniP = n_i_point * n_j_point
    geom.n_face = n_i_face  * n_j_face * 2

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))
    allocate(edge(4, n + 1, 2))

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
    do j = 1, n_j_point
        do i = 1, n_i_point
            numbering = n_i_point * (j - 1) + i

            call Exam_OpenGeo_Cross_Point(edge(1, i, :), edge(2, i, :), edge(3, j, :), edge(4, j, :), xy)

            geom.iniP(numbering).pos(1) = xy(1)
            geom.iniP(numbering).pos(2) = xy(2)
            geom.iniP(numbering).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do j = 1, n_j_face
        do i = 1, n_i_face
            if(pattern == "\") then             ! mesh pattern, \
                numbering = 2*(n_i_face * (j-1) + i) - 1
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j-1) + i
                geom.face(numbering).poi(2) = n_i_point * (j-1) + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i

                numbering = 2*(n_i_face * (j-1) + i)
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j-1) + i + 1
                geom.face(numbering).poi(2) = n_i_point * j + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i
            else if(pattern == "/") then        ! mesh pattern, /
                numbering = 2*(n_i_face * (j-1) + i) - 1
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
                geom.face(numbering).poi(2) = n_i_point * j + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i

                numbering = 2*(n_i_face * (j-1) + i)
                
                geom.face(numbering).n_poi = 3
                allocate(geom.face(numbering).poi(3))
                geom.face(numbering).poi(1) = n_i_point * (j - 1) + i
                geom.face(numbering).poi(2) = n_i_point * (j - 1) + i + 1
                geom.face(numbering).poi(3) = n_i_point * j + i + 1
            end if
        end do
    end do

    ! deallocate edge data
    deallocate(edge)
end subroutine Exam_OpenGeo_Plate_Distorted_Tri

! ---------------------------------------------------------------------------------------

! Example of circular plate with quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Circular_Plate_Quad(prob, geom)
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

    prob.name_file = "57_Circular_Plate_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Circular Plate Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

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

    call Exam_OpenGeo_Merge_Point_Face_Quad(n, joint, conn, count_n, count_e, geom)

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

    call Exam_OpenGeo_Merge_Point_Face_Quad(n, joint, conn, count_n, count_e, geom)

    ! Deallocate data structure
    deallocate(edge, joint, conn)
end subroutine Exam_OpenGeo_Circular_Plate_Quad

! ---------------------------------------------------------------------------------------

! Example of circular plate with triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Circular_Plate_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision,allocatable  :: edge(:,:,:), joint(:,:)
	integer, allocatable :: conn(:,:)

	double precision :: dx, dy
    integer :: i, j, n, count_n, count_e, count_n_t, count_e_t
    character :: pattern
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "58_Circular_Plate_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Tetrahedron"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n           = 3
    pattern     = "\"
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
            if(pattern == "\") then
                count_e = count_e + 1
                

                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 1
                geom.face(count_e).poi(2) = (n+1)*(i+0) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 1

                count_e = count_e + 1
                

                geom.face(count_e).poi(1) = (n+1)*(i+0) + (j-1) + 2
                geom.face(count_e).poi(2) = (n+1)*(i-1) + (j-1) + 2
                geom.face(count_e).poi(3) = (n+1)*(i-1) + (j-1) + 1
            else if(pattern == "/") then
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
            if(pattern == "\") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i+0) + (j-1) + 1
            else if(pattern == "/") then
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

    call Exam_OpenGeo_Merge_Point_Face_Tri(n, joint, conn, count_n, count_e, geom)

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
            if(pattern == "\") then
                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i-1) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 1

                count_e_t = count_e_t + 1
                conn(count_e_t, 1) = (n+1)*(i+0) + (j-1) + 1
                conn(count_e_t, 2) = (n+1)*(i+0) + (j-1) + 2
                conn(count_e_t, 3) = (n+1)*(i-1) + (j-1) + 2
            else if(pattern == "/") then
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

    call Exam_OpenGeo_Merge_Point_Face_Tri(n, joint, conn, count_n, count_e, geom)

    ! Deallocate data structure
    deallocate(edge, joint, conn)
end subroutine Exam_OpenGeo_Circular_Plate_Tri

! ---------------------------------------------------------------------------------------

! Example of annular plate with quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Annular_Plate_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: i_radius, o_radius, angle, radius, force
    integer :: i, j, count, n, nx, nr
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "59_Annular_Plate_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Annular Plate Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "PRESET" .or. para_fig_view == "preset") para_fig_view = "XY"

    n        = 3
    o_radius = 4.0d0        ! Outer radius
    i_radius = 2.4d0        ! Internal radius
    nx       = n
    nr       = 10 * n

    geom.n_iniP = (nx + 1) * nr
    geom.n_face = nx * nr

    ! Allocate point structure and set position vector
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
                count = count + 1
                
                geom.face(count).n_poi = 4
                allocate(geom.face(count).poi(4))
                geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (nx+1)*(i-0) + (j-1) + 2
                geom.face(count).poi(4) = (nx+1)*(i-0) + (j-1) + 1
            else if(i == nr) then
                count = count + 1
                geom.face(count).n_poi = 4
                allocate(geom.face(count).poi(4))
                geom.face(count).poi(1) = (nx+1)*(i-1) + (j-1) + 1
                geom.face(count).poi(2) = (nx+1)*(i-1) + (j-1) + 2
                geom.face(count).poi(3) = (nx+1)*(1-1) + (j-1) + 2
                geom.face(count).poi(4) = (nx+1)*(1-1) + (j-1) + 1
            end if
        end do
    end do
end subroutine Exam_OpenGeo_Annular_Plate_Quad

! ---------------------------------------------------------------------------------------

! Example of annular plate with triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Annular_Plate_Tri(prob, geom)
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
end subroutine Exam_OpenGeo_Annular_Plate_Tri

! ---------------------------------------------------------------------------------------

! Example of hyperbolic paraboloid with quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Hyperbolic_Paraboloid_Quad(prob, geom)
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
end subroutine Exam_OpenGeo_Hyperbolic_Paraboloid_Quad

! ---------------------------------------------------------------------------------------

! Example of hyperbolic paraboloid with triangles
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Hyperbolic_Paraboloid_Tri(prob, geom)
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
end subroutine Exam_OpenGeo_Hyperbolic_Paraboloid_Tri

! ---------------------------------------------------------------------------------------

! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Cross_Point(AP1, AP2, BP1, BP2, CP)
    double precision, intent(in) :: AP1(2), AP2(2), BP1(2), BP2(2)
    double precision, intent(out) :: CP(2)

    double precision :: t, buf

    buf = (BP2(2)-BP1(2))*(AP2(1)-AP1(1))-(BP2(1)-BP1(1))*(AP2(2)-AP1(2))
    t   = (BP2(1)-BP1(1))*(AP1(2)-BP1(2))-(BP2(2)-BP1(2))*(AP1(1)-BP1(1))
    t   = t/buf

    CP(1) = AP1(1) + t * (AP2(1)-AP1(1))
    CP(2) = AP1(2) + t * (AP2(2)-AP1(2))
end subroutine Exam_OpenGeo_Cross_Point

! ---------------------------------------------------------------------------------------

! Merge points and faces
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Merge_Point_Face_Quad(n, joint, conn, n_node, n_element, geom) 
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
            if( Exam_OpenGeo_Comp_XYZ( joint(i,:), geom.iniP(j).pos(:) ) == 1) then
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
                if(Exam_OpenGeo_Comp_XYZ(joint(conn(i,j), :), geom.iniP(k).pos(:)) == 1) exit
            end do
            
            geom.face(n_element).poi(j) = k
        end do
    end do
end subroutine Exam_OpenGeo_Merge_Point_Face_Quad

! ---------------------------------------------------------------------------------------

! Merge points and faces
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_OpenGeo_Merge_Point_Face_Tri(n, joint, conn, n_node, n_element, geom)  
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
            if( Exam_OpenGeo_Comp_XYZ( joint(i,:), geom.iniP(j).pos(:) ) == 1) then
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
                if(Exam_OpenGeo_Comp_XYZ(joint(conn(i,j), :), geom.iniP(k).pos(:)) == 1) exit
            end do
            geom.face(n_element).poi(j) = k
        end do
    end do
end subroutine Exam_OpenGeo_Merge_Point_Face_Tri

! ---------------------------------------------------------------------------------------

! Last updated on Tuesday 30 August 2016 by Hyungmin
integer function Exam_OpenGeo_Comp_XYZ(a, b)
   double precision, intent(in) :: a(3), b(3) 
   double precision :: d

   Exam_OpenGeo_Comp_XYZ = 0
   d = dsqrt((a(1)-b(1))**2.0d0 + (a(2)-b(2))**2.0d0 + (a(3)-b(3))**2.0d0)
   if(d<1.0d-5) Exam_OpenGeo_Comp_XYZ = 1
end function Exam_OpenGeo_Comp_XYZ

! ---------------------------------------------------------------------------------------

end module Exam_OpenGeo