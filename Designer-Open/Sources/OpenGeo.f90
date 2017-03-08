!
! ---------------------------------------------------------------------------------------
!
!                                   Module for OpenGeo
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/10/28
!                                                            Last  modified  : 2016/07/15
!
! ---------------------------------------------------------------------------------------
!
module OpenGeo

    use Data_Prob
    use Data_Geom
    use Data_Bound

    use Mani
    use Math

    implicit none

    public  OpenGeo_Check

    private OpenGeo_Check_Geometry_Type
    private OpenGeo_Set_Nei_Boundary
    private OpenGeo_Count_New_Line_Point
    private OpenGeo_Resize_Line_Point
    private OpenGeo_Add_New_Line_Point
    private OpenGeo_Add_Corner_2Arm
    private OpenGeo_Add_Corner_3Arm
    private OpenGeo_Add_Boundary
    private OpenGeo_Check_Quad_Triangle
    private OpenGeo_Reverse_Face_Connectivity
    private OpenGeo_Write_Init_Open_Geometry

contains

! ---------------------------------------------------------------------------------------

! check geometry type and modify point and line at boundary if open-geometry
subroutine OpenGeo_Check(prob, geom, bound)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    integer :: i, n_new_line, n_new_point

    ! print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |                  2. Open-geometry modification                     |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
    end do

    ! set reverse face connectivity
    !call OpenGeo_Reverse_Face_Connectivity(face)

    ! check open-geometry and set boundary data for open-geometry
    !call OpenGeo_Check_Geometry_Type(prob, geom, bound)

    if(prob.type_geo == "open") then

        ! set neighbor points at boundary in open-geometry
        call OpenGeo_Set_Nei_Boundary(geom, bound)

        ! count the number of new lines and points at the boundary
        call OpenGeo_Count_New_Line_Point(geom, bound, n_new_line, n_new_point)

        ! resize line and point data
        call OpenGeo_Resize_Line_Point(geom, bound, n_new_line, n_new_point)

        ! add new lines and points at the boundary
        call OpenGeo_Add_New_Line_Point(geom, bound, n_new_line, n_new_point)

        ! write initial open-geometry
        call OpenGeo_Write_Init_Open_Geometry(prob, geom, bound)

    else

        ! print progress
        do i = 0, 11, 11
            call Space(i, 11)
            write(i, "(a)"), "-> It is closed geometric problem"
            write(i, "(a)")
        end do
    end if
end subroutine OpenGeo_Check

! ---------------------------------------------------------------------------------------

! reverse face connectivity
subroutine OpenGeo_Reverse_Face_Connectivity(geom)
    type(GeomType), intent(inout) :: geom

    integer, allocatable, dimension(:) :: conn_face
    integer :: i, j

    do i = 1, geom.n_face

        ! allocate data structure to save face connectivity temporarily
        allocate(conn_face(geom.face(i).n_poi))

        ! copy the face connectivity from original data to temporal data
        do j = 1, geom.face(i).n_poi
            conn_face(j) = geom.face(i).poi(j)
        enddo

        ! reverse face connectivity
        do j = 1, geom.face(i).n_poi
            geom.face(i).poi(j) = conn_face(geom.face(i).n_poi+1-j)
        enddo

        ! deallocate face data
        deallocate(conn_face)

    end do
end subroutine OpenGeo_Reverse_Face_Connectivity

! ---------------------------------------------------------------------------------------

! check open-geometry and set boundary data for open-geometry
subroutine OpenGeo_Check_Geometry_Type(prob, geom, bound)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(inout) :: bound

    integer :: point1_ref, point2_ref, point1_com, point2_com
    integer :: i, j, m, n, count

    ! initialize type of geometry
    prob.type_geo = "closed"

    ! allocate and initialize bounday data
    allocate(bound.outer(geom.n_iniP))
    do i = 1, geom.n_iniP
        ! - 1: internal, 1: outer, 0: newly generated point
        bound.outer(i).typeP = -1

        ! initalize data to zero
        bound.outer(i).n_neiP = 0
        bound.outer(i).n_newP = 0
    end do

    ! if line shares two faces it is internal line
    ! reference loop to check open or closed geometry
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi

            count = 0

            if(j /= geom.face(i).n_poi) then
                point1_ref = geom.face(i).poi(j)
                point2_ref = geom.face(i).poi(j + 1)
            else
                point1_ref = geom.face(i).poi(j)
                point2_ref = geom.face(i).poi(1)
            end if

            ! comparision loop
            do m = 1, geom.n_face
                if(i == m) cycle
                do n = 1, geom.face(m).n_poi

                    if(n /= geom.face(m).n_poi) then
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(n + 1)
                    else
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(1)
                    end if

                    ! check boundary line
                    if( (point1_ref == point1_com .and. point2_ref == point2_com) .or. &
                        (point1_ref == point2_com .and. point2_ref == point1_com) ) then
                        count = count + 1
                    end if

                end do
            end do

            ! if open-side edge, boundary point : 1, internal point : -1
            if(count == 0) then
                prob.type_geo = "open"
                bound.outer(point1_ref).typeP  = 1
                bound.outer(point2_ref).typeP  = 1
                bound.outer(point1_ref).n_neiP = 1
                bound.outer(point2_ref).n_neiP = 1
            end if

        end do
    end do
end subroutine OpenGeo_Check_Geometry_Type
    
! ---------------------------------------------------------------------------------------

! set neighbor points at boundary in open-geometry
subroutine OpenGeo_Set_Nei_Boundary(geom, bound)
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(inout) :: bound

    integer, allocatable, dimension(:) :: conn
    integer :: point1_ref, point2_ref, point1_com, point2_com
    integer :: i, j, m, n, count, sort_count, count_front, count_back

    ! count neighbor line for open points
    do i = 1, geom.n_iniP

        ! skip the loop if internal point
        if(bound.outer(i).typeP <= 0) cycle

        count = 1
        point1_ref = i

        do j = 1, geom.n_face
            do m = 1, geom.face(j).n_poi
                point1_com = geom.face(j).poi(m)
                if(point1_ref == point1_com) then
                    count = count + 1
                end if
            end do
        end do

        bound.outer(i).n_neiP = count
        allocate(bound.outer(i).neiP(count))

    end do

    ! set line connectivity
    do i = 1, geom.n_iniP

        ! skip the loop if internal point
        if(bound.outer(i).typeP <= 0) cycle

        point1_ref = i
        count = 0

        do j = 1, geom.n_iniL
            point1_com = geom.iniL(j).poi(1)
            point2_com = geom.iniL(j).poi(2)

            if(point1_ref == point1_com) then
                count = count + 1
                bound.outer(i).neiP(count) = point2_com
            else if(point1_ref == point2_com) then
                count = count + 1
                bound.outer(i).neiP(count) = point1_com
            end if
        end do

        ! sort point in bound connectivity
        count_front = 1     ! first two neighbor is used to determine line neighbor
        count_back  = 3

        ! allocate conn data with # of neighbor points
        allocate(conn(bound.outer(i).n_neiP))

        do j = 1, bound.outer(i).n_neiP

            sort_count = 0
            point1_ref = i
            point2_ref = bound.outer(i).neiP(j)

            do m = 1, geom.n_face
                do n = 1, geom.face(m).n_poi

                    if(n /= geom.face(m).n_poi) then
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(n + 1)
                    else
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(1)
                    end if

                    ! check boundary points
                    if( (point1_ref == point1_com .and. point2_ref == point2_com) .or. &
                        (point1_ref == point2_com .and. point2_ref == point1_com) ) then
                        sort_count = sort_count + 1
                    end if

                end do
            end do

            ! swap to last index of array
            if(sort_count == 1) then        ! outer line
                conn(count_front) = bound.outer(i).neiP(j)
                count_front = count_front + 1
            else if(sort_count == 2) then   ! internal line
                conn(count_back) = bound.outer(i).neiP(j)
                count_back = count_back + 1
            end if

        end do

        ! copy from conn to original data
        do j = 1, bound.outer(i).n_neiP
            bound.outer(i).neiP(j) = conn(j)
        end do

        ! deallocate conn data
        deallocate(conn)

    end do
end subroutine OpenGeo_Set_Nei_Boundary

! ---------------------------------------------------------------------------------------

! count the number of newly added lines and points at boundary
subroutine OpenGeo_Count_New_Line_Point(geom, bound, n_addL, n_addP)
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(inout) :: n_addL, n_addP

    integer :: point_ref, point1_com, point2_com, nface
    integer :: i, j, m, n

    ! count the number of additional points and lines
    n_addP = 0

    do i = 1, geom.n_iniP

        ! skip the loop if internal point
        if(bound.outer(i).typeP <= 0) cycle

        if(bound.outer(i).n_neiP == 2) then

            ! corner at boundary with 4-node surface element
            n_addP = n_addP + 2

            ! allocate new point data in boudary data
            bound.outer(i).n_newP = 2
            allocate(bound.outer(i).newP(2, 2))
            bound.outer(i).newP(1, 1:2) = -1
            bound.outer(i).newP(2, 1:2) = -1

        else if(bound.outer(i).n_neiP >= 3) then

            ! check types of surface
            point_ref = i
            do m = 1, geom.n_face
                do n = 1, geom.face(m).n_poi
                    if(n /= geom.face(m).n_poi) then
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(n + 1)
                    else
                        point1_com = geom.face(m).poi(n)
                        point2_com = geom.face(m).poi(1)
                    end if

                    ! check boundary points
                    if( point_ref == point1_com .or. point_ref == point2_com ) then
                        nface = geom.face(m).n_poi
                    end if

                end do
            end do

            if(nface == 3 .and. bound.outer(i).n_neiP == 3) then

                ! corner with triangular element
                n_addP = n_addP + 3

                ! allocate new point data in boudary data
                bound.outer(i).n_newP = 3
                allocate(bound.outer(i).newP(3, 2))
                do j = 1, 3
                    bound.outer(i).newP(j, 1:2) = -1
                end do

            else

                ! general case of side line
                n_addP = n_addP + bound.outer(i).n_neiP - 2

                ! allocate new point data in boudary data
                bound.outer(i).n_newP = bound.outer(i).n_neiP - 2
                allocate(bound.outer(i).newP(bound.outer(i).n_neiP - 2, 2))
                do j = 1, bound.outer(i).n_neiP - 2
                    bound.outer(i).newP(j, 1:2) = -1
                end do

            end if

        end if
    end do
end subroutine OpenGeo_Count_New_Line_Point

! ---------------------------------------------------------------------------------------

! resize line and point data
subroutine OpenGeo_Resize_Line_Point(geom, bound, n_addL, n_addP)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(inout) :: n_addL, n_addP

    type(OuterType), allocatable, dimension(:) :: temp_bound
    double precision, allocatable, dimension(:, :) :: pos_temp
    integer, allocatable, dimension(:,:)  :: line_temp
    integer :: i, j

    ! # of additional line is the same to the one of additional points
    n_addL = n_addP

    ! copy to temporary line, point and bound data
    ! temporary line data
    allocate(line_temp(geom.n_iniL, 2))
    do i = 1, geom.n_iniL
        line_temp(i, 1:2) = geom.iniL(i).poi(1:2)
    end do

    ! temporary point data
    allocate(pos_temp(geom.n_iniP, 3))
    do i = 1, geom.n_iniP
        pos_temp(i, 1:3) = geom.iniP(i).pos(1:3)
    end do

    ! temporary bound data
    allocate(temp_bound(geom.n_iniP))
    do i = 1, geom.n_iniP
        temp_bound(i).typeP  = bound.outer(i).typeP
        temp_bound(i).n_neiP = bound.outer(i).n_neiP
        temp_bound(i).n_newP = bound.outer(i).n_newP

        ! skip the loop if internal point
        if(bound.outer(i).typeP <= 0) cycle

        allocate(temp_bound(i).neiP(bound.outer(i).n_neiP))
        do j = 1, bound.outer(i).n_neiP
            temp_bound(i).neiP(j) = bound.outer(i).neiP(j)
        end do
        
        allocate(temp_bound(i).newP(bound.outer(i).n_newP, 2))
        do j = 1, bound.outer(i).n_newP
            temp_bound(i).newP(j, 1:2) = bound.outer(i).newP(j, 1:2)
        end do
    end do

    ! deallocate original data
    deallocate(geom.iniL)
    deallocate(geom.iniP)
    
    do i = 1, geom.n_iniP
        if(bound.outer(i).typeP /= -1) &
            deallocate(bound.outer(i).neiP, bound.outer(i).newP)
    end do
    deallocate(bound.outer)

    ! increase # of lines and points due to boundary lines and points
    geom.n_iniP = geom.n_iniP + n_addP
    geom.n_iniL = geom.n_iniL + n_addL

    ! reallocate original data with increased # of lines and points
    allocate(geom.iniL(geom.n_iniL))
    allocate(geom.iniP(geom.n_iniP))
    allocate(bound.outer(geom.n_iniP))

    ! copy from temporary data to updated original data
    ! line data
    do i = 1, geom.n_iniL
        if(i <= geom.n_iniL - n_addL) then
            geom.iniL(i).poi(1:2) = line_temp(i, 1:2)
        else
            geom.iniL(i).poi(1) = -1
            geom.iniL(i).poi(2) = -1
        end if

        geom.iniL(i).sec  = -1
        geom.iniL(i).iniL = -1
        geom.iniL(i).neiP(1, 1:2) = -1
        geom.iniL(i).neiP(2, 1:2) = -1
    end do

    ! point data
    do i = 1, geom.n_iniP
        if(i <= geom.n_iniP - n_addP) then
            geom.iniP(i).pos(1:3) = pos_temp(i, 1:3)
        else
            geom.iniP(i).pos(1:3) = 0.0d0
        end if
    end do

    ! bound data
    do i = 1, geom.n_iniP
        if(i <= geom.n_iniP - n_addP) then
            
            bound.outer(i).typeP  = temp_bound(i).typeP
            bound.outer(i).n_neiP = temp_bound(i).n_neiP
            bound.outer(i).n_newP = temp_bound(i).n_newP

            ! skip the loop if internal point
            if(temp_bound(i).typeP <= 0) cycle

            allocate(bound.outer(i).neiP(temp_bound(i).n_neiP))
            do j = 1, bound.outer(i).n_neiP
                bound.outer(i).neiP(j) = temp_bound(i).neiP(j)
            end do

            allocate(bound.outer(i).newP(temp_bound(i).n_newP, 2))
            do j = 1, bound.outer(i).n_newP
                bound.outer(i).newP(j, 1:2) = temp_bound(i).newP(j, 1:2)
            end do
        else
            ! type of new generated point should be zero
            bound.outer(i).typeP  = 0
            bound.outer(i).n_neiP = 0
            bound.outer(i).n_newP = 0
        end if
    end do

    ! deallocate temporary data
    deallocate(line_temp)
    deallocate(pos_temp)
    do i = 1, geom.n_iniP - n_addP
        if(temp_bound(i).typeP /= -1) &
            deallocate(temp_bound(i).neiP, temp_bound(i).newP)
    end do
    deallocate(temp_bound)
end subroutine OpenGeo_Resize_Line_Point

! ---------------------------------------------------------------------------------------

! add new lines and points at the boundary
subroutine OpenGeo_Add_New_Line_Point(geom, bound, n_addL, n_addP)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(in)    :: n_addL, n_addP

    integer :: i, j, chk_face, count, iline, ipoint

    ! set boundary points and lines
    count = 0
    do i = 1, geom.n_iniP

        ! skip the loop if internal point
        if(bound.outer(i).typeP <= 0) cycle

        ! for 2-point boundary arm for quad surface element
        if(bound.outer(i).n_neiP == 2) then

            ! add two lines and points
            do j = 1, 2
                count  = count + 1
                iline  = geom.n_iniL - n_addL + count
                ipoint = geom.n_iniP - n_addP + count

                ! add lines and points at corner with 2 arm (qaud)
                call OpenGeo_Add_Corner_2Arm(geom, bound, i, j, iline, ipoint)
            end do

        else if(bound.outer(i).n_neiP >= 3) then  ! > 3-points boundary arm

            ! check triangular or quad element
            call OpenGeo_Check_Quad_Triangle(geom, i, chk_face)

            ! for triangular mesh with 3-arm at the boundary
            if(chk_face == 3 .and. bound.outer(i).n_neiP == 3) then

                do j = 1, 3
                    count  = count + 1
                    iline  = geom.n_iniL - n_addL + count
                    ipoint = geom.n_iniP - n_addP + count

                    ! add lines and points at corner with 3 arm (triangle)
                    call OpenGeo_Add_Corner_3Arm(geom, bound, i, j, iline, ipoint)
                end do

            else    ! for general cases

                do j = 1, bound.outer(i).n_neiP - 2
                    count  = count + 1
                    iline  = geom.n_iniL - n_addL + count
                    ipoint = geom.n_iniP - n_addP + count

                    ! add lines and points at boundary
                    call OpenGeo_Add_Boundary(geom, bound, i, j, iline, ipoint)
                end do
            end if
        end if
    end do

    ! print information
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), " 1.1. Open-geometry modification (outer points))"
        write(i, "(a)")
    end do

    do j = 1, geom.n_iniP
        ! skip the loop if internal point
        if(bound.outer(j).typeP <= 0) cycle

        write(11, "(i10, a, 25i4)"), j, " th point : neighbor(sorted) point index : ", &
            bound.outer(j).neiP(1:bound.outer(j).n_neiP)
        write(11, "(a, 25i4)"), "                      added point index            : ", &
            bound.outer(j).newP(1:bound.outer(j).n_newP, 1)
        write(11, "(a, 25i4)"), "                      heritage point index         : ", &
            bound.outer(j).newP(1:bound.outer(j).n_newP, 2)
    end do
end subroutine OpenGeo_Add_New_Line_Point

! ---------------------------------------------------------------------------------------

! add lines and points at corner with 2 arm(quad)
subroutine OpenGeo_Add_Corner_2Arm(geom, bound, ibound, count, iline, ipoint)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(in)    :: ibound, count, iline, ipoint

    integer :: i, point1_com, point2_com
    double precision :: vec(3), length

    length = 0.5d0

    ! set new point position
    vec(1:3) = geom.iniP(ibound).pos(1:3) - geom.iniP(bound.outer(ibound).neiP(count)).pos(1:3)
    geom.iniP(ipoint).pos(1:3) = length * vec(1:3) + geom.iniP(ibound).pos(1:3)

    do i = 1, geom.n_iniL
        point1_com = geom.iniL(i).poi(1)
        point2_com = geom.iniL(i).poi(2)

        ! outward to the boundary
        if(point1_com == ibound .and. point2_com == bound.outer(ibound).neiP(count)) then

            ! new point direction is outward
            geom.iniL(iline).poi(1) = ibound
            geom.iniL(iline).poi(2) = ipoint

            ! neighbor points related to point 2 are free
            geom.iniL(iline).neiP(2, 1) = -1
            geom.iniL(iline).neiP(2, 2) = -1

            ! set line neighbor connectivity
            if(count == 1) then
                geom.iniL(iline).neiP(1, 1) = ipoint + 1
                geom.iniL(iline).neiP(1, 2) = bound.outer(ibound).neiP(2)
            else if(count == 2) then
                geom.iniL(iline).neiP(1, 1) = ipoint - 1
                geom.iniL(iline).neiP(1, 2) = bound.outer(ibound).neiP(1)
            end if

            exit

        ! inward to the boundary
        else if(point2_com == ibound .and. point1_com == bound.outer(ibound).neiP(count)) then

            ! new point direction is outward
            geom.iniL(iline).poi(1) = ibound
            geom.iniL(iline).poi(2) = ipoint

            ! neighbor points related to point 2 are free
            geom.iniL(iline).neiP(2, 1) = -1
            geom.iniL(iline).neiP(2, 2) = -1

            ! set line neighbor connectivity
            if(count == 1) then
                geom.iniL(iline).neiP(1, 1) = bound.outer(ibound).neiP(2)
                geom.iniL(iline).neiP(1, 2) = ipoint + 1
            else if(count == 2) then
                geom.iniL(iline).neiP(1, 2) = ipoint - 1
                geom.iniL(iline).neiP(1, 1) = bound.outer(ibound).neiP(1)
            end if

            exit

        end if
    end do

    ! set index and heritage index for new point
    bound.outer(ibound).newP(count, 1) = ipoint
    bound.outer(ibound).newP(count, 2) = bound.outer(ibound).neiP(count)
end subroutine OpenGeo_Add_Corner_2Arm

! ---------------------------------------------------------------------------------------

! add lines and points at corner with 3 arm (triangle)
subroutine OpenGeo_Add_Corner_3Arm(geom, bound, ibound, count, iline, ipoint)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(in)    :: ibound, count, iline, ipoint

    integer :: k, point1_com, point2_com
    double precision :: vec(3), length

    length = 0.5d0

    ! set point position
    vec(1:3) = geom.iniP(ibound).pos(1:3) - geom.iniP(bound.outer(ibound).neiP(count)).pos(1:3)
    geom.iniP(ipoint).pos(1:3) = length * vec(1:3) + geom.iniP(ibound).pos(1:3)

    ! set line connectivity
    geom.iniL(iline).poi(1) = ibound
    geom.iniL(iline).poi(2) = ipoint

    ! set new point index, heritage
    bound.outer(ibound).newP(count, 1) = ipoint
    bound.outer(ibound).newP(count, 2) = bound.outer(ibound).neiP(count)
end subroutine OpenGeo_Add_Corner_3Arm

! ---------------------------------------------------------------------------------------

! add lines and points at boundary
subroutine OpenGeo_Add_Boundary(geom, bound, ibound, count, iline, ipoint)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound
    integer,         intent(in)    :: ibound, count, iline, ipoint

    integer :: i, point1_com, point2_com
    double precision :: vec(3), length

    length = 0.5d0

    ! set new point position
    vec(1:3) = geom.iniP(ibound).pos(1:3) - geom.iniP(bound.outer(ibound).neiP(count + 2)).pos(1:3)
    geom.iniP(ipoint).pos(1:3) = length * vec(1:3) + geom.iniP(ibound).pos(1:3)

    ! new point direction is outward
    geom.iniL(iline).poi(1) = ibound
    geom.iniL(iline).poi(2) = ipoint

    ! set neighbor points related to point 2
    geom.iniL(iline).neiP(2, 1) = -1
    geom.iniL(iline).neiP(2, 2) = -1

    ! find direction to find neighbor points
    do i = 1, geom.n_iniL
        point1_com = geom.iniL(i).poi(1)
        point2_com = geom.iniL(i).poi(2)

        ! outward to the boundary
        if(point1_com == ibound .and. point2_com == bound.outer(ibound).neiP(count)) then

            ! set line neighbor connectivity
            geom.iniL(iline).neiP(1, 1) = bound.outer(ibound).neiP(1)
            geom.iniL(iline).neiP(1, 2) = bound.outer(ibound).neiP(2)

            exit

        ! inward to the boundary
        else if(point2_com == ibound .and. point1_com == bound.outer(ibound).neiP(count)) then

            ! set line neighbor connectivity
            geom.iniL(iline).neiP(1, 1) = bound.outer(ibound).neiP(2)
            geom.iniL(iline).neiP(1, 2) = bound.outer(ibound).neiP(1)

            exit

        end if
    end do

    ! set index and heritage index for new point
    bound.outer(ibound).newP(count, 1) = ipoint
    bound.outer(ibound).newP(count, 2) = bound.outer(ibound).neiP(count + 2)
end subroutine OpenGeo_Add_Boundary

! ---------------------------------------------------------------------------------------

! check triangular or quad element
subroutine OpenGeo_Check_Quad_Triangle(geom, idxP, chk_face)
    type(GeomType), intent(in)    :: geom
    integer,        intent(in)    :: idxP
    integer,        intent(inout) :: chk_face

    integer :: m, n, point_ref, point1_com, point2_com

    point_ref = idxP
    do m = 1, geom.n_face
        do n = 1, geom.face(m).n_poi
            if(n /= geom.face(m).n_poi) then
                point1_com = geom.face(m).poi(n)
                point2_com = geom.face(m).poi(n + 1)
            else
                point1_com = geom.face(m).poi(n)
                point2_com = geom.face(m).poi(1)
            end if

            ! check boundary points
            if( point_ref == point1_com .or. point_ref == point2_com ) then
                chk_face = geom.face(m).n_poi
            end if
        end do
    end do
end subroutine OpenGeo_Check_Quad_Triangle

! ---------------------------------------------------------------------------------------

! write initial open geometry
subroutine OpenGeo_Write_Init_Open_Geometry(prob, geom, bound)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound

    double precision :: length, pos_1(3), pos_2(3), pos_c(3), vec(3)
    logical :: f_info, f_axis
    integer :: i, j

    f_axis = .false.
    f_info = .false.

    open(unit=102, file=trim(prob.path_work1)//trim(prob.name_file)//&
        "_init_open_geo.bild", form="formatted")

    ! write initial line for open geometry
    write(102, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
        write(102,"(a10, 7f9.3)"), ".cylinder ", pos_1(1:3), pos_2(1:3), 0.2d0
    end do

    ! write initial point for open-geometry
    do i = 1, geom.n_iniP
        if(prob.type_geo == "closed") then
            write(102, "(a)"), ".color 53"
        else if(prob.type_geo == "open") then
            if(bound.outer(i).typeP == -1) then         ! Internal points
                write(102, "(a)"), ".color red"
            else if(bound.outer(i).typeP == 0) then     ! New generated points
                write(102, "(a)"), ".color black"
            else                                        ! Outer points
                write(102, "(a)"), ".color 53"
            end if
        end if

        write(102, "(a, 4f9.3)"), ".sphere ", geom.iniP(i).pos(1:3), 0.5d0

    end do

    ! information on initial geometry
    if(f_info == .true.) then

        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0
            vec(1:3)   = Normalize_Vector(pos_2(1:3) - pos_1(1:3))
            length     = Size_Vector(pos_2 - pos_1)

            write(102, "(a, 3f9.3)"), ".cmov ", pos_c(1:3) + 0.5d0
            write(102, "(a)"), ".color black"
            write(102, "(a)"), ".font Helvetica 12 bold"
            write(102, "(i, a, f5.2, a)"), i, "(", length, ")"
            write(102, "(a)"), ".color blue"
            write(102, "(a, 8f8.2)"), ".arrow ", pos_c(1:3), pos_c(1:3) + 1.5d0 * vec(1:3), 0.25d0, 0.5d0
        end do

        do i = 1, geom.n_iniP
            write(102, "(a, 3f9.3)"), ".cmov ", geom.iniP(i).pos(1:3) + 1.0d0
            write(102, "(a)"), ".color red"
            write(102, "(a)"), ".font Helvetica 12 bold"
            write(102, "(i)"), i
        end do

        ! for faces
        do i = 1, geom.n_face

            ! find center position in mesh
            pos_c(1:3) = 0.0d0
            do j = 1, geom.face(i).n_poi
                pos_c(1:3) = pos_c(1:3) + geom.iniP(geom.face(i).poi(j)).pos(1:3)
            end do
            pos_c(1:3) = pos_c(1:3) / dble(geom.face(i).n_poi)

            ! write face number
            write(102, "(a, 3f9.3)"), ".cmov ", pos_c(1:3) + 1.0d0
            write(102, "(i7)"), i
        end do
        
    end if

    ! write global axis
    if(f_axis == .true.) then

        write(102, "(a)"), ".translate 0.0 0.0 0.0"
        write(102, "(a)"), ".scale 0.5"
        write(102, "(a)"), ".color grey"
        write(102, "(a)"), ".sphere 0 0 0 0.5"      ! center
        write(102, "(a)"), ".color red"             ! x-axis
        write(102, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(102, "(a)"), ".color blue"            ! y-axis
        write(102, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(102, "(a)"), ".color yellow"          ! z-axis
        write(102, "(a)"), ".arrow 0 0 0 0 0 4 "

    end if

    close(unit=102)
end subroutine OpenGeo_Write_Init_Open_Geometry

! ---------------------------------------------------------------------------------------

end module OpenGeo