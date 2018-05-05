!
! =============================================================================
!
! Module - Importer
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
module Importer

    use Ifport

    use Para
    use Data_Prob
    use Data_Geom

    use Math

    implicit none

    public Importer_PLY
    public Importer_STL
    public Importer_WRL
    public Importer_GEO

contains

! -----------------------------------------------------------------------------

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

    path = "input/"//trim(prob.name_file)//"."//trim(prob.type_file)
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

! -----------------------------------------------------------------------------

! import STL format using meshconv
subroutine Importer_STL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    ! change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_STL

! -----------------------------------------------------------------------------

! import WRL format using meshconv
subroutine Importer_WRL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    ! change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_WRL

! -----------------------------------------------------------------------------

! Import .geo, iges and igs format to convert polygon meshes
subroutine Importer_GEO(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: p_mesh
    integer :: i, j, n_poi, n_line, temp
    logical :: results
    character(200) :: path, fullname

    ! Data structure for meshing
    type :: MeshType
        integer :: cn(100)   ! Maximum connectivity
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable :: face_con(:)

    fullname = trim(prob.name_file)//"."//trim(prob.type_file)
    path     = "input/"//fullname

    ! Read number of points and faces
    if(prob.type_file == 'geo') then
        open(unit=1002, file=path, form="formatted")
        read(1002, *), geom.n_iniP, n_line, geom.n_face
    end if

    ! Boundary & internal mesh design
    if(geom.n_face == 0) then
        if(prob.type_file == 'geo') close(unit=1002)

        write(0, "(a)"), "Converting geometry with faced mesh"

        ! Convert to face meshes from lines
        if(para_platform == "dev") results = systemqq(trim("python tools/Shapely/Shapely.py")//" input/"//trim(fullname))
        !if(para_platform == "win") results = systemqq(trim("python tools/Shapely/Shapely.py")//" input/"//trim(fullname))
        if(para_platform == "win") results = systemqq(trim("tools\Shapely\Shapely.exe")//" input/"//trim(fullname))
        if(para_platform == "mac") results = systemqq(trim("python tools/Shapely/Shapely.py")//" input/"//trim(fullname))

        fullname = trim(prob.name_file)//trim("_shapely.geo")
        open(unit=1002, file="input/"//trim(fullname), form="formatted")

        ! Read number of points and faces
        read(1002, *), geom.n_iniP, n_line, geom.n_face
    end if

    ! Boundary design
    if(geom.n_face == 1) then
        write(0, "(a)")
        write(0, "(a)"), "   Type the mesh spacing parameter (0.0 ~ 1.0) [Enter] : "
        read(*, *), p_mesh

        if(p_mesh > 0.0d0 .and. p_mesh < 1.0) then
            close(unit=1002)

            if(para_platform == "dev") then
                ! MATLAB - DistMesh
                ! matlab -wait -nodisplay -nosplash -nodesktop -r
                ! "addpath tools/DistMesh/src; addpath tools/DistMesh;
                ! DistMesh('input/ex_des1_shapely.geo', 0.3)
                results = systemqq("matlab -wait -nodisplay -nosplash -nodesktop -r "//&
                    '"addpath tools/DistMesh/src; addpath tools/DistMesh; DistMesh('//&
                    "'input/"//trim(fullname)//"',"//trim(Dble2Str(p_mesh))//')"')

                ! Python - DistMesh
                ! python tools/PyDistMesh/PyDistMesh.py input/ex_des1_shapely.geo 0.3
                !results = systemqq(&
                !trim("python tools/PyDistMesh/PyDistMesh.py input/")&
                !    //trim(fullname)//' '//trim(Dble2Str(p_mesh)))
            else if(para_platform == "win") then

                ! MATLAB exe - DistMesh
                ! tools\DistMesh\DistMesh.exe input\ex_des1_shapely.geo 0.3
                results = systemqq('tools\DistMesh\DistMesh.exe input\'&
                    //trim(fullname)//' '//trim(Dble2Str(p_mesh)))
            else if(para_platform == "mac") then

                ! Python - DistMesh
                ! python tools/PyDistMesh/PyDistMesh.py input/ex_des1_shapely.geo 0.3
                results = systemqq(&
                    trim("python tools/PyDistMesh/PyDistMesh.py input/")&
                    //trim(fullname)//' '//trim(Dble2Str(p_mesh)))
            end if

            fullname = trim(prob.name_file)//trim("_shapely_distmesh.geo")
            open(unit=1002, file="input/"//trim(fullname), form="formatted")
            read(1002, *), geom.n_iniP, n_line, geom.n_face
        end if
    end if

    ! Read point data
    allocate(geom.iniP(geom.n_iniP))
    do i = 1, geom.n_iniP
        read(1002, *), temp, geom.iniP(i).pos(1:2)
        geom.iniP(i).pos(3) = 0.0d0

        geom.iniP(i).pos(2) = geom.iniP(i).pos(2)
    end do

    ! Read face
    allocate(geom.face(geom.n_face))
    allocate(face_con(geom.n_face))
    do i = 1, geom.n_face

        face_con(i).cn(:) = 0
        read(1002, *), temp, n_poi, face_con(i).cn(1:n_poi)

        geom.face(i).n_poi = n_poi
        allocate(geom.face(i).poi(n_poi))

        do j = 1, n_poi
            !geom.face(i).poi(j) = face_con(i).cn(n_poi-j+1)
            geom.face(i).poi(j) = face_con(i).cn(j)
        end do
    end do

    ! Deallocate memory and close file
    deallocate(face_con)
    close(unit=1002)

    ! Delete temp file
    if(para_platform == "win" .or. para_platform == "dev") then
        results = systemqq(trim("del input\")//trim(prob.name_file)//trim("_shapely.geo"))
    else
        results = systemqq(trim("rm input/")//trim(prob.name_file)//trim("_shapely.geo"))
    end if

    if(para_platform == "win" .or. para_platform == "dev") then
        results = systemqq(trim("del input\")//trim(prob.name_file)//trim("_shapely_distmesh.geo"))
    else
        results = systemqq(trim("rm input/")//trim(prob.name_file)//trim("_shapely_distmesh.geo"))
    end if
end subroutine Importer_GEO

! -----------------------------------------------------------------------------

end module Importer