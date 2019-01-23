!
! =============================================================================
!
! Module - Importer
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
    public Importer_SVG

contains

! -----------------------------------------------------------------------------

! Import PLY file
subroutine Importer_PLY(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(200) :: ctemp, path
    integer :: i, status, npoint, ioerr
    logical :: here

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable, dimension (:) :: Basepair_con

    path = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file)

    inquire(file = trim(path), exist=here)
    if(here == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | No ply file.                                     |"
        write(p_redir, "(a)"), " +==================================================+"
        write(p_redir, "(a)")
        stop
    end if
    open(unit = 1001, file = trim(path), form = "formatted")

    do
        read(1001, "(a100)", iostat = status), ctemp

        ! Negative value, means the end-of-file (EOF) mark was read
        if(status < 0) exit

        if(index(ctemp, "format")) then

            ! Read the number of points
            read(1001, *, iostat = ioerr), ctemp, ctemp, geom.n_iniP
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The format for # is not correct in GEO file.     |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
            allocate(geom.iniP(geom.n_iniP))

        else if(index(ctemp, "property float32 z") .or. index(ctemp, "property float z")) then

            ! Read the number of faces
            read(1001, *, iostat = ioerr), ctemp, ctemp, geom.n_face
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The format for faces is not correct in PLY file. |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

        else if(index(ctemp,"end_header")) then

            ! Read point
            do i = 1, geom.n_iniP
                read(1001, *, iostat = ioerr), geom.iniP(i).pos(1:3)
                if(ioerr /= 0) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error ========================================+"
                    write(p_redir, "(a)"), " | The format for points is not correct in PLY file.|"
                    write(p_redir, "(a)"), " +==================================================+"
                    stop
                end if
            end do

            allocate(geom.face(geom.n_face))
            allocate(Basepair_con(geom.n_face))

            ! Read face connectivity
            do i = 1, geom.n_face
                read(1001, *, iostat = ioerr), npoint, Basepair_con(i).cn(1:npoint)
                if(ioerr /= 0) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error ========================================+"
                    write(p_redir, "(a)"), " | The format for conns is not correct in PLY file. |"
                    write(p_redir, "(a)"), " +==================================================+"
                    stop
                end if
                geom.face(i).n_poi = npoint
                allocate(geom.face(i).poi(npoint))

                geom.face(i).poi(1:npoint) = Basepair_con(i).cn(1:npoint)
            end do

        end if
    end do

    deallocate(Basepair_con)
    close(unit = 1001)
end subroutine Importer_PLY

! -----------------------------------------------------------------------------

! Import STL format using meshconv
subroutine Importer_STL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! Run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    if(results == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Failed to load the meshconv.                     |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_STL

! -----------------------------------------------------------------------------

! Import WRL format using meshconv
subroutine Importer_WRL(prob)
    type(ProbType), intent(inout) :: prob

    logical :: results

    ! Run meshconv to generate *.PLY fileformat
    results = systemqq(trim("Library/meshconv -c ply Input/")// &
        trim(prob.name_file)//"."//trim(prob.type_file)//trim(" -ascii"))

    if(results == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Failed to load the meshconv.                     |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Change file type to PLY
    prob.type_file = "ply"
end subroutine Importer_WRL

! -----------------------------------------------------------------------------

! Import .geo, and igs file format
subroutine Importer_GEO(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: p_mesh
    integer :: i, j, n_poi, n_line, temp, ioerr
    logical :: results, here
    character(200) :: path, fullname, fullname1, fullname2

    ! Data structure for meshing
    type :: MeshType
        integer :: cn(100)
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable :: face_con(:)

    fullname = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file)

    inquire(file = trim(fullname), exist = here)
    if(here == .false.) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | No geo file.                                     |"
        write(p_redir, "(a)"), " +==================================================+"
        write(p_redir, "(a)")
        stop
    end if

    ! Read number of points and faces
    if(prob.type_file == 'geo') then

        inquire(file = trim(fullname), exist = here)
        if(here == .false.) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | No geo file.                                     |"
            write(p_redir, "(a)"), " +==================================================+"
            write(p_redir, "(a)")
            stop
        end if
        open(unit = 1002, file = trim(fullname), form = "formatted")

        read(1002, *, iostat = ioerr), geom.n_iniP, n_line, geom.n_face
        if(ioerr /= 0) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | The format1 for # is not correct in GEO file.    |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if
    end if

    ! Boundary and internal mesh design
    if(geom.n_face == 0) then
        if(prob.type_file == 'geo') close(unit = 1002)

        write(0, "(a)"), "   * Converting lines to meshes by Shapely"

        ! Convert to face meshes from lines
        !DEC$ IF DEFINED(_WIN32)
        results = systemqq("python tools\Shapely\Shapely.py "//trim(fullname))
        !results = systemqq("tools\Shapely\Shapely.exe "//trim(fullname))
        !DEC$ ELSE
        results = systemqq("python /home/polyhedra/PERDIX/make/tools/Shapely/Shapely.py "//trim(fullname))
        !results = systemqq("python tools/Shapely/Shapely.py "//trim(fullname))
        !DEC$ ENDIF
        if(results == .false.) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | Failed to load Python Shapely.                   |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if

        fullname = trim(prob.path_input)//trim(prob.name_file)//trim("_shapely.geo")

        ! Check file
        inquire(file = trim(fullname), exist = here)
        if(here == .false.) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | Fail to convert lines to meshes by Shapely       |"
            write(p_redir, "(a)"), " +==================================================+"
            write(p_redir, "(a)")
            stop
        end if
        open(unit = 1002, file = trim(fullname), form = "formatted")

        ! Read number of points and faces
        read(1002, *, iostat = ioerr), geom.n_iniP, n_line, geom.n_face
        if(ioerr /= 0) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | The format2 for # is not correct in GEO file.    |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if
    end if

    ! Exception module
    if(geom.n_face == 0) then
        close(1002)
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Some faces are open.                             |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Boundary design
    if(geom.n_face == 1) then
        if(iargc() == 0) then

            write(0, "(a)")
            write(0, "(a)"), "   Type the value (0.0 ~ 1.0) for the mesh spacing parameter [Enter]"
            write(0, "(a)"), "   * The small value generates finer meshes"
            write(0, "(a)"), "   * '1' for the structure without mesh"
            read(*, *), p_mesh
        else

            p_mesh = prob.p_mesh
        end if

        if(p_mesh < 0.0d0 .or. p_mesh > 1.0d0) then
            close(1002)
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | The spacing parameter should be from 0.0 to 1.0. |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if

        if(p_mesh >= 0.0d0 .and. p_mesh < 1.0) then
            ! Interpolation of the parameter value
            ! 0.0          x          1.0
            ! 0.2          y          0.7
            p_mesh = 0.2d0 + (0.7d0 - 0.2d0) * 1.0d0 * (p_mesh)
            close(unit = 1002)

            write(0, "(a)")
            write(0, "(a)"), "   * Generating the internal mesh by DistMesh"
            write(0, "(a)")

            !DEC$ IF DEFINED(_WIN32)
            !! === DEV
            !! MATLAB - DistMesh
            !! matlab -wait -nodisplay -nosplash -nodesktop -r
            !! "addpath tools/DistMesh/src; addpath tools/DistMesh;
            !! DistMesh('ex_des1_shapely.geo', 0.3)
            !results = systemqq("matlab -wait -nodisplay -nosplash -nodesktop -r "//&
            !    '"addpath tools/DistMesh/src; addpath tools/DistMesh; DistMesh('//&
            !    trim(fullname)//"',"//trim(Dble2Str(p_mesh))//')"')
            !
            !! Python - DistMesh
            !! python tools/PyDistMesh/PyDistMesh.py ex_des1_shapely.geo 0.3
            !!results = systemqq(trim("python tools/PyDistMesh/PyDistMesh.py ")//trim(fullname)//' '//trim(Dble2Str(p_mesh)))
            !! === DEV

            ! MATLAB exe - DistMesh
            ! tools\DistMesh\DistMesh.exe ex_des1_shapely.geo 0.3
            !results = systemqq('tools\DistMesh\DistMesh.exe '//trim(fullname)//' '//trim(Dble2Str(p_mesh)))
            results = systemqq("matlab -wait -nodisplay -nosplash -nodesktop -r "//'"addpath tools/DistMesh/src; addpath tools/DistMesh; DistMesh('//"'"//trim(fullname)//"',"//trim(Dble2Str(p_mesh))//')"')
            !DEC$ ELSEIF DEFINED(__APPLE__)

            ! Python - DistMesh
            ! python tools/PyDistMesh/PyDistMesh.py ex_des1_shapely.geo 0.3
            results = systemqq(trim("python tools/PyDistMesh/PyDistMesh.py ")//trim(fullname)//" "//trim(Dble2Str(p_mesh)))
            !DEC$ ELSEIF DEFINED(__linux)

            ! MATLAB - DistMesh
            results = systemqq("/usr/local/MATLAB/R2018b/bin/matlab wait -nodisplay -nosplash -nodesktop -noawt -r "//&
                '"addpath /home/polyhedra/PERDIX/make/tools/DistMesh/src; addpath /home/polyhedra/PERDIX/make/tools/DistMesh; DistMesh_linux('&
                //"'"//trim(fullname)//"',"//trim(Dble2Str(p_mesh))//')"')
            !results = systemqq("matlab wait -nodisplay -nosplash -nodesktop -noawt -r "//&
            !    '"addpath /home/polyhedra/PERDIX/make/tools/DistMesh/src; addpath /home/polyhedra/PERDIX/make/tools/DistMesh; DistMesh_linux('&
            !    //"'"//trim(fullname)//"',"//trim(Dble2Str(p_mesh))//')"')
            !DEC$ ENDIF
            if(results == .false.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Failed to load DistMesh due to the curvature.    |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            fullname = trim(prob.name_file)//trim("_shapely_distmesh.geo")

            ! Check file
            inquire(file = fullname, exist = here)
            if(here == .false.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Fail to generate the mesh by DistMesh            |"
                write(p_redir, "(a)"), " +==================================================+"
                write(p_redir, "(a)")
                stop
            end if
            open(unit = 1002, file = trim(prob.path_input)//trim(fullname), form = "formatted")

            read(1002, *, iostat = ioerr), geom.n_iniP, n_line, geom.n_face
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The format3 for # is not correct in GEO file.    |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
        end if
    end if

    ! Read point data
    allocate(geom.iniP(geom.n_iniP))
    do i = 1, geom.n_iniP
        read(1002, *, iostat = ioerr), temp, geom.iniP(i).pos(1:2)
        if(ioerr /= 0) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | The format for points is not correct in GEO file.|"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if
        geom.iniP(i).pos(3) = 0.0d0
        geom.iniP(i).pos(2) = geom.iniP(i).pos(2)
    end do

    ! Read face
    allocate(geom.face(geom.n_face))
    allocate(face_con(geom.n_face))
    do i = 1, geom.n_face

        face_con(i).cn(:) = 0
        read(1002, *, iostat = ioerr), temp, n_poi, face_con(i).cn(1:n_poi)
        if(ioerr /= 0) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | The format for faces is not correct in GEO file. |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if

        geom.face(i).n_poi = n_poi
        allocate(geom.face(i).poi(n_poi))

        do j = 1, n_poi
            !geom.face(i).poi(j) = face_con(i).cn(n_poi-j+1)
            geom.face(i).poi(j) = face_con(i).cn(j)
        end do
    end do

    ! Deallocate memory and close file
    deallocate(face_con)
    close(unit = 1002)

    ! Delete temp files
    fullname1 = trim(prob.path_input)//trim(prob.name_file)//trim("_shapely.geo")
    fullname2 = trim(prob.path_input)//trim(prob.name_file)//trim("_shapely_distmesh.geo")

    !DEC$ IF DEFINED(_WIN32)
    inquire(file = trim(fullname1), exist=here)
    if(here == .true.) results = systemqq(trim("del")//" "//trim(fullname1))
    inquire(file = trim(fullname2), exist=here)
    if(here == .true.) results = systemqq(trim("del")//" "//trim(fullname2))
    !DEC$ ELSE
    inquire(file = trim(fullname1), exist=here)
    if(here == .true.) results = systemqq(trim("rm -rf")//" "//trim(fullname1))
    inquire(file = trim(fullname2), exist=here)
    if(here == .true.) results = systemqq(trim("rm -rf")//" "//trim(fullname2))
    !DEC$ ENDIF
end subroutine Importer_GEO

! -----------------------------------------------------------------------------

! Import .svg format to convert .geo file
subroutine Importer_SVG(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x1, y1, x2, y2
    integer :: i, j, n_line, temp, length, ioerr
    character(len=200) :: fullname
    character(len=100) :: text, word, xx, cx1, cx2, cy1, cy2
    character(len=*), parameter :: search_str = "<line "
    logical :: here

    fullname = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file)

    ! Read number of points and faces
    n_line = 0
    if(prob.type_file == 'svg') then

        inquire(file = trim(fullname), exist = here)
        if(here == .false.) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | No svg file.                                     |"
            write(p_redir, "(a)"), " +==================================================+"
            write(p_redir, "(a)")
            stop
        end if
        open(unit = 1002, file = trim(fullname), action="read")

        ! Check the number of lines
        do
            ! Read line into character variable
            read(1002, "(a)", iostat = ioerr), text
            if(ioerr /= 0) exit

            ! Read word line
            read(text, *), word

            ! Found search string at beginning of line
            if(word == search_str) then
                n_line = n_line + 1
            end if
        end do
        close(unit = 1002)

        open(unit = 1003, file = trim(prob.path_input)//trim(prob.name_file)//".geo", form = "formatted")
        write(1003, "(3i)"), n_line * 2, n_line, 0

        n_line = 0
        open(unit = 1002, file = trim(fullname), action="read")
        do
            ! Read line into character variable
            read(1002, "(a)", iostat = ioerr), text
            if(ioerr /= 0) exit

            ! Read word line
            read(text, *), word

            ! Found search string at beginning of line
            if(word == search_str) then
                n_line = n_line + 1
                read(text, *, iostat = ioerr), xx, xx, cx1, cy1, cx2, cy2
                if(ioerr /= 0) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error ========================================+"
                    write(p_redir, "(a)"), " | The format is not correct in SVG file.           |"
                    write(p_redir, "(a)"), " +==================================================+"
                    stop
                end if

                length = len_trim(cx1); xx = cx1(5:length-1); read(xx,*), x1
                length = len_trim(cy1); xx = cy1(5:length-1); read(xx,*), y1
                length = len_trim(cx2); xx = cx2(5:length-1); read(xx,*), x2
                length = len_trim(cy2); xx = cy2(5:length-1); read(xx,*), y2

                write(1003, "(i, 2f)"), 2*n_line-1, x1, y1
                write(1003, "(i, 2f)"), 2*n_line,   x2, y2
            end if
        end do
    end if

    ! Construct the line connectivity
    do i = 1, n_line
        write(1003, "(3i)"), i, 2*i-1, 2*i
    end do

    close(unit = 1002)
    close(unit = 1003)

    prob.type_file = "geo"
    call Importer_GEO(prob, geom)
end subroutine Importer_SVG

! -----------------------------------------------------------------------------

end module Importer