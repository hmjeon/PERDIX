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
    !public Exam_Open3D_Cooling_Tower_Quad           ! Cooling tower with quad mesh

    ! Tri mesh
    public Exam_Open3D_End_Triangular_Prism_Tri     ! 11. Open end triangular prism with tri mesh
    public Exam_Open3D_End_Cube_Tri                 ! 12. Open end cube with tri mesh
    public Exam_Open3D_End_Pentagonal_Prism_Tri     ! 13. Open end pentagonal prism with tri mesh
    public Exam_Open3D_End_Cylinder_Tri             ! 14. Open end cylinder with tri mesh
    public Exam_Open3D_Hemisphere_Tri               ! 15. Hemisphere with tri mesh
    public Exam_Open3D_Cooling_Tower_Tri            ! Cooling tower with tri mesh

contains

! ---------------------------------------------------------------------------------------

! Example of open end triangular prism with quad mesh
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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

    ! Set geometric type and view (atom, cylinder size, move_x, move_y)
    call Mani_Set_View_Color(prob, [231, 76, 60], "xyz", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

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

! ! Example of cooling tower with tri mesh
subroutine Exam_Open3D_Cooling_Tower_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: i, j, index, n_i_poi, n_j_poi, n, nx, ny
    character(10) :: char_sec, char_bp, char_start_bp

    double precision, allocatable :: jnt(:,:)
    integer, allocatable :: con(:,:), out_con(:,:)

    double precision :: length, angle, radius
    integer :: jj(6), n_jnt, n_con

    !jj(:) = [1, 2, 3, 4, 3, 2]      ! "/" mesh
    jj(:) = [2, 4, 1, 3, 1, 4]      ! "\" mesh

    write(0, "(a)"), "[COOLING TOWER UNIFORM]"
    write(0, "(a)")

    n       = 10

    ! Allocate memories
    allocate(con(n*n, 4), jnt((n+1)*n, 3))

    ! Generation of joints
    n_jnt = 0
    do i = 1, n + 1
        do j = 1, n

            n_jnt = n_jnt + 1
            angle = (360.0d0 - dble(j-1)*(360.0d0/dble(n)))*(pi/180.0d0)

            jnt(n_jnt,2) = 1.0d0 - dble(i-1)*(2.0d0/dble(n))

            radius = dsqrt(0.5d0 + jnt(n_jnt,2)*jnt(n_jnt,2))
            jnt(n_jnt,1) = radius * dcos(angle)
            jnt(n_jnt,3) = radius * dsin(angle)
        end do
    end do

    ! Connectivities
    con(:,:) = 0
    n_con = 0
    do i = 1, n
        do j = 1, n - 1

            n_con = n_con + 1

            con(n_con, 1) = (i-1)*(n) + (j-1) + 1
            con(n_con, 2) = (i-1)*(n) + (j-1) + 2
            con(n_con, 3) = (i-1)*(n) + (j-1) + (n) + 1
            con(n_con, 4) = (i-1)*(n) + (j-1) + (n) + 2
        end do

        n_con = n_con + 1

        con(n_con, 1) = (i-1)*(n) + (n-2) + 2
        con(n_con, 2) = (i-1)*(n) + (1-1) + 1
        con(n_con, 3) = (i-1)*(n) + (n-2) + (n) + 2
        con(n_con, 4) = (i-1)*(n) + (1-1) + (n) + 1
    end do

    allocate(out_con(2*n_con, 3))
    call Prob_Elem_Mesh(jj, con, n_con, out_con, 3)


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
    call Mani_Set_View_Color(prob, [52, 152, 219], "xy", 1.0d0, 1.0d0, 0.0d0, 0.0d0)

    geom.n_iniP = n_jnt
    geom.n_face = n_con
    
    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    ! Set position vector
    do i = 1, n_jnt
        geom.iniP(i).pos(1:3) = jnt(i,:)
    end do

    ! Set connectivity
    do i = 1, n_con
        geom.face(i).poi(1:3) = out_con(i, 1:3)
    end do

    ! Deallocate
    deallocate(jnt, con, out_con)
end subroutine Exam_Open3D_Cooling_Tower_Tri

! ---------------------------------------------------------------------------------------

subroutine Prob_Elem_Mesh(con, element, n_con, out_con, nNPE)
    integer, intent(in)    :: con(:)
    integer, intent(in)    :: element(:,:), nNPE
    integer, intent(inout) :: n_con
    integer, intent(out)   :: out_con(:,:)

    integer :: i, j, k

    j = 0
    do i = 1, n_con
        j = j + 1
        do k = 1, nNPE
            out_con(j,k) = element(i, con(k))
        end do
        j = j + 1
        do k = 1, nNPE
            out_con(j,k) = element(i, con(nNPE+k))
        end do
    end do
    n_con = j
end subroutine Prob_Elem_Mesh

! ---------------------------------------------------------------------------------------

end module Exam_3D_Open