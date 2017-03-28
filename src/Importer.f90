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
! Last updated on Monday 7 Mar 2016 by Hyungmin
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
! Last updated on Monday 7 Mar 2016 by Hyungmin
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
! Last updated on Monday 7 Mar 2016 by Hyungmin
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

! import geom format that is standard for PERDIX
! Last updated on Monday 7 Mar 2016 by Hyungmin
subroutine Importer_GEO(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i, number, npoint
    character(200) :: path

    ! mesh data structure
    type :: MeshType
        integer :: cn(100)   ! maximum connectivity
    end type MeshType

    ! 1st: # of meshes, 2nd: points
    type(MeshType), allocatable, dimension (:) :: Basepair_con

    path = "Input\"//trim(prob.name_file)//"."//trim(prob.type_file)
    open(unit=1002, file=path, form="formatted")

    ! read points
    read(1002, *), geom.n_iniP
    allocate(geom.iniP(geom.n_iniP))

    do i = 1, geom.n_iniP
        read(1002, *), number, geom.iniP(i).pos(1:3)
    end do

    ! read face
    read(1002, *), geom.n_face
    allocate(geom.face(geom.n_face))
    allocate(Basepair_con(geom.n_face))

    ! read # of vectices in the mesh and connectivity
    do i = 1, geom.n_face
        read(1002, *), npoint, Basepair_con(i).cn(1:npoint)

        geom.face(i).n_poi = npoint
        allocate(geom.face(i).poi(npoint))

        geom.face(i).poi(1:npoint) = Basepair_con(i).cn(1:npoint)
    end do

    deallocate(Basepair_con)
end subroutine Importer_GEO

! ---------------------------------------------------------------------------------------

end module Importer