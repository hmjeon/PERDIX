!
! ---------------------------------------------------------------------------------------
!
!                               Module - Exam_3D_Open
!
!                                                                    Updated : 2017/03/27
!
! Comments: 3D open geometry
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Exam_3D_Open

    use Data_Prob
    use Data_Geom

    use Para
    use Mani
    use Math

    implicit none

    ! Quad mesh
    public Exam_Open3D_End_Triangular_Prism_Quad    ! 11. Open end triangular prism with quad mesh
    public Exam_Open3D_End_Cube_Quad                ! 12. Open end cube with quad mesh
    public Exam_Open3D_End_Pentagonal_Prism_Quad    ! 13. Open end pentagonal prism with quad mesh
    public Exam_Open3D_End_Cylinder_Quad            ! 14. Open end cylinder with quad mesh
    public Exam_Open3D_Hemisphere_Quad              ! 15. Hemisphere with quad mesh

    ! Tri mesh
    public Exam_Open3D_End_Triangular_Prism_Tri     ! 11. Open end triangular prism with tri mesh
    public Exam_Open3D_End_Cube_Tri                 ! 12. Open end cube with tri mesh
    public Exam_Open3D_End_Pentagonal_Prism_Tri     ! 13. Open end pentagonal prism with tri mesh
    public Exam_Open3D_End_Cylinder_Tri             ! 14. Open end cylinder with tri mesh
    public Exam_Open3D_Hemisphere_Tri               ! 15. Hemisphere with tri mesh

contains

! ---------------------------------------------------------------------------------------

! Example of open end triangular prism with quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Triangular_Prism_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "11_End_Triangular_Prism_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Triangular Prism Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 2
    nz = n
    nr = 3

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Triangular_Prism_Quad

! ---------------------------------------------------------------------------------------

! Example of open end cube with quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Cube_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "12_End_Cube_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Cube Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = 4

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Cube_Quad

! ---------------------------------------------------------------------------------------

! Example of open end pentagonal prism with quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Pentagonal_Prism_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "13_End_Pentagonal_Prism_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Pentagonal Prism Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = 5

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Pentagonal_Prism_Quad

! ---------------------------------------------------------------------------------------

! Example of open end cylinder with quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Cylinder_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "14_End_Cylinder_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Cylinder Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 4
    nz = n
    nr = 8

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)  
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set face connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j + 1
            geom.face(n_face).poi(4) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + 1
        geom.face(n_face).poi(4) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Cylinder_Quad

! ---------------------------------------------------------------------------------------

! Example of hemisphere with quad mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_Hemisphere_Quad(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang1, ang2, length
    integer :: i, j, n, nz, nr, index
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "15_Hemisphere_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hemisphere Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = n * 3

    ! Set dimension
    rad =  10.0d0

    geom.n_face = nz * nr
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    ! Set position vector
    index = 0
    do i = 1, nz + 1
        ang2 = (72.0d0-dble(i-1)*(72.0d0/dble(nz)))*(pi/180.0d0)

        do j = 1, nr
            index = index + 1
            ang1 = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(index).pos(1) = rad*cos(ang2)*cos(ang1)
            geom.iniP(index).pos(2) = rad*cos(ang2)*sin(ang1)
            geom.iniP(index).pos(3) = rad*sin(ang2)
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nz
        do j = 1, nr - 1
            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j
            geom.face(index).poi(2) = (i - 1) * nr + j + 1
            geom.face(index).poi(3) = i * nr + j + 1
            geom.face(index).poi(4) = i * nr + j
        end do

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + j
        geom.face(index).poi(2) = (i - 1) * nr + 1
        geom.face(index).poi(3) = i * nr + 1
        geom.face(index).poi(4) = i * nr + j
    end do
end subroutine Exam_Open3D_Hemisphere_Quad

! ---------------------------------------------------------------------------------------

! Example of open end triangular prism with tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Triangular_Prism_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "11_End_Triangular_Prism_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Triangular Prism Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = 3

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Triangular_Prism_Tri

! ---------------------------------------------------------------------------------------

! Example of open end cube with tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Cube_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "12_End_Cube_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Cube Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = 4

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Cube_Tri

! ---------------------------------------------------------------------------------------

! Example of open end pentagonal prism with tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Pentagonal_Prism_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "13_End_Pentagonal_Prism_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Pentagonal Prism Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = 5

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Pentagonal_Prism_Tri

! ---------------------------------------------------------------------------------------

! Example of open end cylinder with tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_End_Cylinder_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang, length
    integer :: i, j, n, nz, nr, n_poi, n_face
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "14_End_Cylinder_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "End Cylinder Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 5
    nz = n
    nr = 8

    ! Set dimension
    rad    = 1.0d0
    length = 2.0d0 * rad * dsin(pi/dble(nr)) * dble(nz)

    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    n_poi = 0
    do i = 1, nz + 1
        do j = 1, nr
            n_poi = n_poi + 1
            ang   = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(n_poi).pos(1) = rad*dcos(ang)
            geom.iniP(n_poi).pos(2) = -(i-1)*(length/nz)
            geom.iniP(n_poi).pos(3) = rad*dsin(ang)
        end do
    end do

    ! Set connectivity
    n_face = 0
    do i = 1, nz
        do j = 1, nr - 1
            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j
            geom.face(n_face).poi(2) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j

            n_face = n_face + 1
            geom.face(n_face).poi(1) = (i - 1) * nr + j + 1
            geom.face(n_face).poi(2) = i * nr + j + 1
            geom.face(n_face).poi(3) = i * nr + j
        end do

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + j
        geom.face(n_face).poi(2) = (i - 1) * nr + 1
        geom.face(n_face).poi(3) = i * nr + j

        n_face = n_face + 1
        geom.face(n_face).poi(1) = (i - 1) * nr + 1
        geom.face(n_face).poi(2) = i * nr + 1
        geom.face(n_face).poi(3) = i * nr + j
    end do
end subroutine Exam_Open3D_End_Cylinder_Tri

! ---------------------------------------------------------------------------------------

! Example of hemisphere with tri mesh
! Last updated on Fri 10 Mar 2017 by Hyungmin
subroutine Exam_Open3D_Hemisphere_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: rad, ang1, ang2, length
    integer :: i, j, n, nz, nr, index
    character(10) :: char_sec, char_bp, char_start_bp

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "15_Hemisphere_Tri"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Hemisphere Tri"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0           ! Atomic model
    prob.size     = 1.0d0           ! Cylindrical model
    prob.move_x   = 0.0d0           ! Cylindrical model
    prob.move_y   = 0.0d0           ! Cylindrical model
    para_fig_view = "XY"

    ! Set mesh
    n  = 3
    nz = n
    nr = n * 3

    ! Set dimension
    rad =  10.0d0

    geom.n_face = nz * nr * 2
    geom.n_iniP = (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    index = 0
    do i = 1, nz + 1
        ang2 = (72.0d0-dble(i-1)*(72.0d0/dble(nz)))*(pi/180.0d0)

        do j = 1, nr
            index = index + 1
            ang1 = (360.0d0-dble(j-1)*(360.0d0/dble(nr)))*(pi/180.0d0)

            geom.iniP(index).pos(1) = rad*cos(ang2)*cos(ang1)
            geom.iniP(index).pos(2) = rad*cos(ang2)*sin(ang1)
            geom.iniP(index).pos(3) = rad*sin(ang2)
        end do
    end do

    ! Set connectivity
    index = 0
    do i = 1, nz
        do j = 1, nr - 1
            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j
            geom.face(index).poi(2) = (i - 1) * nr + j + 1
            geom.face(index).poi(3) = i * nr + j

            index = index + 1
            geom.face(index).poi(1) = (i - 1) * nr + j + 1
            geom.face(index).poi(2) = i * nr + j + 1
            geom.face(index).poi(3) = i * nr + j
        end do

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + j
        geom.face(index).poi(2) = (i - 1) * nr + 1
        geom.face(index).poi(3) = i * nr + j

        index = index + 1
        geom.face(index).poi(1) = (i - 1) * nr + 1
        geom.face(index).poi(2) = i * nr + 1
        geom.face(index).poi(3) = i * nr + j
    end do
end subroutine Exam_Open3D_Hemisphere_Tri

! ---------------------------------------------------------------------------------------

end module Exam_3D_Open