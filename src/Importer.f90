!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Importer
!
!                                                                    Updated : 2017/03/27
!
! Comments: This module is for the import specific file type.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Importer

    use Ifport

    use Data_Prob
    use Data_Geom

    use Math

    implicit none

    public Importer_PLY
    public Importer_STL
    public Importer_WRL
    public Importer_GEO

contains

! ---------------------------------------------------------------------------------------

! import PLY format that has surface mesh information
subroutine Importer_PLY(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, status, npoint
    character(200) :: ctemp, path

    ! mesh data structure
    type :: MeshType
        integer :: cn(100)   ! maximum connectivity
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable, dimension (:) :: Basepair_con

    path = "..\examples\"//trim(prob.name_file)//"."//trim(prob.type_file)
    open(unit=1001, file=path, form="formatted")

    do
        read(1001, "(a100)", IOSTAT=status), ctemp

        ! negative value, means the end-of-file (EOF) mark was read
        if(status < 0) exit

        if(index(ctemp, "format")) then

            ! read the number of points
            read(1001, *), ctemp, ctemp, geom.n_iniP
            allocate(geom.iniP(geom.n_iniP))

        else if(index(ctemp, "property float32 z") .or. index(ctemp, "property float z")) then

            ! read the number of faces
            read(1001, *), ctemp, ctemp, geom.n_face

        else if(index(ctemp,"end_header")) then

            ! read point
            do i = 1, geom.n_iniP
                read(1001, *), geom.iniP(i).pos(1:3)
            end do

            allocate(geom.face(geom.n_face))
            allocate(Basepair_con(geom.n_face))

            ! read # of vectices in the mesh and connectivity
            do i = 1, geom.n_face
                read(1001, *), npoint, Basepair_con(i).cn(1:npoint)

                geom.face(i).n_poi = npoint
                allocate(geom.face(i).poi(npoint))

                geom.face(i).poi(1:npoint) = Basepair_con(i).cn(1:npoint)
            end do

        end if
    end do

    deallocate(Basepair_con)
    close(unit=1001)
end subroutine Importer_PLY

! ---------------------------------------------------------------------------------------

! import STL format using meshconv
subroutine Importer_STL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! run meshconv to generate *.PLY fileformat
    results = SYSTEMQQ(trim("Library\meshconv -c ply Input\")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    ! change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_STL

! ---------------------------------------------------------------------------------------

! import WRL format using meshconv
subroutine Importer_WRL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! run meshconv to generate *.PLY fileformat
    results = SYSTEMQQ(trim("Library\meshconv -c ply Input\")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    ! change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_WRL

! ---------------------------------------------------------------------------------------

! Import .geo format and convert face meshes
subroutine Importer_GEO(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, j, n_poi, n_line, temp
    logical :: results, b_face
    character(200) :: path, file

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)   ! Maximum connectivity
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable :: face_con(:)

    b_face = .true.

    file = trim(prob.name_file)//"."//trim(prob.type_file)
    path = "input\"//file
    open(unit=1002, file=path, form="formatted")

    ! Read number of points and faces
    if(prob.type_file == 'geo') read(1002, *), geom.n_iniP, n_line, geom.n_face

    if(prob.type_file == 'iges' .or. prob.type_file == 'igs' .or. n_line /= 0) then

        b_face = .false.

        close(unit=1002)
        write(0, "(a)"), "Converting geometry with faced mesh"

        ! Run convertorv to generate the geometry with face meshes
        results = SYSTEMQQ(trim("Convertor")//" input\"//trim(file))
        !results = SYSTEMQQ(trim("copy ")//trim(prob.name_file)//".tmp input\")
        !results = SYSTEMQQ(trim("del geometry.tmp"))

        path = trim(prob.name_file)//".tmp"
        open(unit=1002, file="input\"//path, form="formatted")

        ! Read number of points and faces
        read(1002, *), geom.n_iniP, n_line, geom.n_face
    end if

    ! Read point data
    allocate(geom.iniP(geom.n_iniP))
    do i = 1, geom.n_iniP
        read(1002, *), temp, geom.iniP(i).pos(1:2)
        geom.iniP(i).pos(3) = 0.0d0

        if(b_face == .true.) geom.iniP(i).pos(2) = -geom.iniP(i).pos(2)
    end do

    ! Read face
    allocate(geom.face(geom.n_face))
    allocate(face_con(geom.n_face))
    do i = 1, geom.n_face

        read(1002, *), temp, n_poi, face_con(i).cn(1:n_poi)

        geom.face(i).n_poi = n_poi
        allocate(geom.face(i).poi(n_poi))

        do j = 1, n_poi
            geom.face(i).poi(j) = face_con(i).cn(n_poi-j+1)
            !geom.face(i).poi(j) = face_con(i).cn(j)
        end do
    end do

    ! Deallocate memory and close file
    deallocate(face_con)
    close(unit=1002)

    ! Delete temp file
    results = SYSTEMQQ(trim("del input\")//trim(prob.name_file)//".tmp")
end subroutine Importer_GEO

! ---------------------------------------------------------------------------------------

end module Importer