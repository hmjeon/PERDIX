!
! ---------------------------------------------------------------------------------------
!
!                                   Module - ModGeo
!
!                                                                    Updated : 2017/03/27
!
! Comments: This module is for geometric modification to generate the seperated lines.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module ModGeo

    use Data_Prob
    use Data_Geom
    use Data_Bound

    use Para
    use List
    use Math

    use Mani

    implicit none

    public  ModGeo_Modification

    private ModGeo_Set_Neighbor_Point
    private ModGeo_Find_Neighbor_Face
    private ModGeo_Find_Neighbor_Points
    private ModGeo_Set_Neighbor_Line
    private ModGeo_Find_Neighbor_Line
    private ModGeo_Chimera_Check_Geometry
    private ModGeo_Set_Junction_Data
    private ModGeo_Set_Local_Coorindate
    private ModGeo_Set_Local_Vectors
    private ModGeo_Chimera_Init_Geometry_Local
    private ModGeo_Seperate_Line
    private ModGeo_Find_Scale_Factor
    private ModGeo_Set_Angle_Junction
    private ModGeo_Set_Width_Section
    private ModGeo_Set_Scale_Geometry
    private ModGeo_Set_Gap_Junction
    private ModGeo_Chimera_Sep_Geometry
    private ModGeo_Set_Const_Geometric_Ratio

contains

! ---------------------------------------------------------------------------------------

! Geometry modification
subroutine ModGeo_Modification(prob, geom, bound)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    double precision :: scale
    integer :: i

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |                2. Build modified edges and points                  |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
    end do

    ! ==================================================
    !
    ! Set local coordinate system
    !
    ! ==================================================
    ! Set neighbor point
    call ModGeo_Set_Neighbor_Point(prob, geom)

    ! Set neighbor line
    call ModGeo_Set_Neighbor_Line(prob, geom, bound)

    ! Chimera check geometry
    call ModGeo_Chimera_Check_Geometry(prob, geom)

    ! Set arm junction data
    call ModGeo_Set_Junction_Data(geom, bound)

    ! Set local coordinate system on each line
    call ModGeo_Set_Local_Coorindate(geom)

    ! Write initial geometry with local coordinate system
    call ModGeo_Chimera_Init_Geometry_Local(prob, geom)

    ! ==================================================
    !
    ! Seperated lines from the vertex
    !
    ! ==================================================
    ! Seperate the line from the vertex without off-set distance
    call ModGeo_Seperate_Line(geom, bound)

    ! Find the factor to scale the geometry
    scale = ModGeo_Find_Scale_Factor(prob, geom, bound)

    ! Scale the geometry with the factor
    call ModGeo_Set_Scale_Geometry(geom, scale)

    ! Write seperated geometry
    !call ModGeo_Chimera_Sep_Geometry(prob, geom, "sep0")

    ! Set the gap distance between the junction and end of edges
    call ModGeo_Set_Gap_Junction(geom, bound)

    ! Set constant modified edge ratio based on original geometry
    if(para_const_edge_mesh == "on") call ModGeo_Set_Const_Geometric_Ratio(geom)

    ! Write seperated geometry
    call ModGeo_Chimera_Sep_Geometry(prob, geom, "sep")
end subroutine ModGeo_Modification

! ---------------------------------------------------------------------------------------

! Set neighbor points from two faces sharing with the line
subroutine ModGeo_Set_Neighbor_Point(prob, geom)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(inout) :: geom

    integer :: i, nei_face(2), nei_poi(2, 2)

    ! nei_poi(a, b), geom.iniL(i).neiP(a, b)
    ! a - point index from the connectivies of the line
    ! b - 1: left and positive sign, 2: right and negative sign

    ! For total number of lines
    do i = 1, geom.n_iniL

        ! Initialize neighbor face and point
        nei_face(1:2)   = -1
        nei_poi(1, 1:2) = -1
        nei_poi(2, 1:2) = -1

        ! Find neighbor face sharing with the i th line
        nei_face(1:2) = ModGeo_Find_Neighbor_Face(geom, i)

        ! Set open geometry
        if(nei_face(1) == -1 .or. nei_face(2) == -1) prob.type_geo = "open"

        ! Set nighbor face, 1: left and positive, 2: right and negative
        geom.iniL(i).neiF(1) = nei_face(1)
        geom.iniL(i).neiF(2) = nei_face(2)

        ! Find neighbor points from two faces
        nei_poi = ModGeo_Find_Neighbor_Points(geom, i, nei_face)

        ! Update
        geom.iniL(i).neiP(1, 1:2) = nei_poi(1, 1:2)
        geom.iniL(i).neiP(2, 1:2) = nei_poi(2, 1:2)
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.1. Neighbor points from two faces sharing with the line"
    end do
    do i = 1, geom.n_iniL
        write(11, "(i20, a, i3, a, i3, a$)"), i, " th edge (", geom.iniL(i).poi(1), "->", geom.iniL(i).poi(2), ") "
        write(11, "(a,   i3$)"), "| (Poi 1) -> left(+) : ", geom.iniL(i).neiP(1, 1)
        write(11, "(a,   i3$)"), ", right(-) : ",           geom.iniL(i).neiP(1, 2)
        write(11, "(a,   i3$)"), ", (Poi 2) -> left(+) : ", geom.iniL(i).neiP(2, 1)
        write(11, "(a,   i3 )"), ", right(-) : ",           geom.iniL(i).neiP(2, 2)
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Set_Neighbor_Point

! ---------------------------------------------------------------------------------------

! Find faces sharing with the current line
function ModGeo_Find_Neighbor_Face(geom, line) result(nei_face)
    type(GeomType), intent(in) :: geom
    integer,        intent(in) :: line
    
    integer :: nei_face(2), poi_1, poi_2, poi_a, poi_b
    integer :: i, j

    ! nei_face(1:2) : face ID sharing with current line
    nei_face(1) = -1
    nei_face(2) = -1

    poi_1 = geom.iniL(line).poi(1)
    poi_2 = geom.iniL(line).poi(2)

    ! For the total number of faces
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi

            ! Set face orientation
            if(j /= geom.face(i).n_poi) then
                poi_a = geom.face(i).poi(j)
                poi_b = geom.face(i).poi(j + 1)
            else
                poi_a = geom.face(i).poi(j)
                poi_b = geom.face(i).poi(1)
            end if

            ! If poi_12 and poi_ab are the same direction
            if(poi_a == poi_1 .and. poi_b == poi_2) then
                nei_face(1) = i
            end if

            ! If poi_12 and poi_ab are the different direction
            if(poi_b == poi_1 .and. poi_a == poi_2) then
                nei_face(2) = i
            end if
        end do
    end do
end function ModGeo_Find_Neighbor_Face

! ---------------------------------------------------------------------------------------

! Find neighbor points from the face sharing the line
function ModGeo_Find_Neighbor_Points(geom, line, face) result(nei_poi)
    type(GeomType), intent(in) :: geom
    integer,        intent(in) :: line
    integer,        intent(in) :: face(2)

    integer :: nei_poi(2,2), i, j, poi_1, poi_2, num_face, npoi_face
    integer :: poi_a, poi_b, poi_first, poi_last

    ! Points from the line connectivities
    poi_1 = geom.iniL(line).poi(1)
    poi_2 = geom.iniL(line).poi(2)

    ! Initialize the return value before calculating
    nei_poi(1, 1:2) = -1
    nei_poi(2, 1:2) = -1

    ! Check the number of faces sharing line
    do i = 1, 2

        ! Check whether it is boundary or not
        if(face(i) /= -1) then
            num_face  = face(i)
            npoi_face = geom.face(num_face).n_poi
        else
            ! For bounday faces
            cycle
        end if

        do j = 1, npoi_face - 1

            poi_first = geom.face(num_face).poi(1)
            poi_last  = geom.face(num_face).poi(npoi_face)
            poi_a     = geom.face(num_face).poi(j)
            poi_b     = geom.face(num_face).poi(j + 1)

            ! If the direction of the line is the same with the face orientation
            if(poi_1 == poi_a .and. poi_2 == poi_b) then

                ! The next line is assigned to left
                if(j + 1 == npoi_face) then
                    nei_poi(2, 1) = poi_first
                else
                    nei_poi(2, 1) = geom.face(num_face).poi(j + 2)
                end if

                ! The previous line is assigned to left
                if(j == 1) then
                    nei_poi(1, 1) = poi_last
                else
                    nei_poi(1, 1) = geom.face(num_face).poi(j - 1)
                end if
            else if(poi_2 == poi_first .and. poi_1 == poi_last) then

                ! At point 2, next line is left
                nei_poi(2, 1) = geom.face(num_face).poi(2)

                ! At point 1, previous line is left
                nei_poi(1, 1) = geom.face(num_face).poi(npoi_face - 1)
            end if

            ! If the direction of the line is the opposite to the face orientation
            if(poi_2 == poi_a .and. poi_1 == poi_b) then

                ! Previous line is right
                if(j == 1) then
                    nei_poi(2, 2) = poi_last
                else
                    nei_poi(2, 2) = geom.face(num_face).poi(j - 1)
                end if

                ! Next line is left
                if(j + 1 == npoi_face) then
                    nei_poi(1, 2) = poi_first
                else
                    nei_poi(1, 2) = geom.face(num_face).poi(j + 2)
                end if
            else if(poi_1 == poi_first .and. poi_2 == poi_last) then

                ! At point 2, previous line is right
                nei_poi(2, 2) = geom.face(num_face).poi(npoi_face - 1)

                ! At point 1, next line is right
                nei_poi(1, 2) = geom.face(num_face).poi(2)
            end if
        end do
    end do
end function ModGeo_Find_Neighbor_Points

! ---------------------------------------------------------------------------------------

! Set neighbor lines from neighbor points
subroutine ModGeo_Set_Neighbor_Line(prob, geom, bound)
    type(ProbType),  intent(in)    :: prob
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    integer :: i, nei_line(2, 2)

    ! nei_line(a, b)
    ! a : point index from the connectivity
    ! b : 1- left and positive, 2- right and negative

    ! Find the neighbor line from neighbor points
    do i = 1, geom.n_iniL

        call ModGeo_Find_Neighbor_Line(geom, i, nei_line)

        ! Update
        geom.iniL(i).neiL(1, 1:2) = nei_line(1, 1:2)
        geom.iniL(i).neiL(2, 1:2) = nei_line(2, 1:2)
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.2. Find neighboring edges from initial geometry"
    end do
    do i = 1, geom.n_iniL
        write(11, "(i20, a$ )"), i, " th line"
        write(11, "(a,   i4$)"), ", (point 1) -> left : ", geom.iniL(i).neiL(1, 1)
        write(11, "(a,   i4$)"), ", right : ",             geom.iniL(i).neiL(1, 2)
        write(11, "(a,   i4$)"), ", (point 2) -> left : ", geom.iniL(i).neiL(2, 1)
        write(11, "(a,   i4 )"), ", right : ",             geom.iniL(i).neiL(2, 2)
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Set_Neighbor_Line

! ---------------------------------------------------------------------------------------

! Find neighbor line from neighbor points
subroutine ModGeo_Find_Neighbor_Line(geom, line, nei_line)
    type(GeomType), intent(in)  :: geom
    integer,        intent(in)  :: line
    integer,        intent(out) :: nei_line(2, 2)

    integer :: i, poi_1, poi_2

    nei_line(1, 1:2) = -1
    nei_line(2, 1:2) = -1

    poi_1 = geom.iniL(line).poi(1)
    poi_2 = geom.iniL(line).poi(2)

    ! Set neighboring lines
    do i = 1, geom.n_iniL

        if( (geom.iniL(line).neiP(1, 1) == geom.iniL(i).poi(1) .and. poi_1 == geom.iniL(i).poi(2)) .or. &
            (geom.iniL(line).neiP(1, 1) == geom.iniL(i).poi(2) .and. poi_1 == geom.iniL(i).poi(1)) ) then
            nei_line(1,1) = i
        end if

        if( (geom.iniL(line).neiP(2, 1) == geom.iniL(i).poi(1) .and. poi_2 == geom.iniL(i).poi(2)) .or. &
            (geom.iniL(line).neiP(2, 1) == geom.iniL(i).poi(2) .and. poi_2 == geom.iniL(i).poi(1)) ) then
            nei_line(2,1) = i
        end if

        if( (geom.iniL(line).neiP(1, 2) == geom.iniL(i).poi(1) .and. poi_1 == geom.iniL(i).poi(2)) .or. &
            (geom.iniL(line).neiP(1, 2) == geom.iniL(i).poi(2) .and. poi_1 == geom.iniL(i).poi(1)) ) then
            nei_line(1,2) = i
        end if

        if( (geom.iniL(line).neiP(2, 2) == geom.iniL(i).poi(1) .and. poi_2 == geom.iniL(i).poi(2)) .or. &
            (geom.iniL(line).neiP(2, 2) == geom.iniL(i).poi(2) .and. poi_2 == geom.iniL(i).poi(1)) ) then
            nei_line(2,2) = i
        end if
    end do
end subroutine ModGeo_Find_Neighbor_Line

! ---------------------------------------------------------------------------------------

! Write Chimera check geometry
subroutine ModGeo_Chimera_Check_Geometry(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: pos_1(3), pos_2(3), pos_c(3), pos_c1(3), pos_c2(3)
    double precision :: vec_a(3), vec_b(3), vec(3)
    integer :: i, j, point_1, point_2
    logical :: f_axis, f_info
    character(200) :: path

    if(para_write_301 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_301_info

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=301, file=trim(path)//"_check_geo.bild", form="formatted")

    ! Write initial lines
    do i = 1, geom.n_iniL
        point_1    = geom.iniL(i).poi(1)
        point_2    = geom.iniL(i).poi(2)
        pos_1(1:3) = geom.iniP(point_1).pos(1:3)
        pos_2(1:3) = geom.iniP(point_2).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        ! Line and orientation vector
        write(301, "(a     )"), ".color dark green"
        write(301, "(a$    )"), ".cylinder "
        write(301, "(3f9.3$)"), pos_1(1:3)
        write(301, "(3f9.3$)"), pos_2(1:3)
        write(301, "(1f9.3 )"), 0.2d0
    end do

    ! Write initial geometric points
    do i = 1, geom.n_iniP
        write(301,"(a     )"), ".color red"
        write(301,"(a$    )"), ".sphere "
        write(301,"(3f9.3$)"), geom.iniP(i).pos(1:3)
        write(301,"(1f9.3 )"), 0.5d0
    end do

    ! Write outward vector
    do i = 1, geom.n_face

        ! Find center position in mesh
        pos_c(1:3) = 0.0d0
        do j = 1, geom.face(i).n_poi
            pos_c(1:3) = pos_c + geom.iniP(geom.face(i).poi(j)).pos
        end do
        pos_c(1:3) = pos_c / dble(geom.face(i).n_poi)

        ! Find orientation
        vec_a(1:3) = geom.iniP(geom.face(i).poi(1)).pos - pos_c
        vec_b(1:3) = geom.iniP(geom.face(i).poi(2)).pos - pos_c
        vec(1:3)   = Normalize(Cross(vec_a, vec_b))

        ! Draw the orientation of the outward vector of the face
        write(301, "(a     )"), ".color salmon"
        write(301, "(a$    )"), ".arrow "
        write(301, "(3f8.2$)"), pos_c(1:3)
        write(301, "(3f8.2$)"), pos_c(1:3) + 3.0d0*vec(1:3)
        write(301, "(3f8.2 )"), 0.2d0, 0.5d0, 0.6d0
    end do

    ! ==================================================
    !
    ! Write information on line and neighbor numbers
    !
    ! ==================================================
    if (f_info == .true.) then

        ! For lines
        do i = 1, geom.n_iniL
            point_1     = geom.iniL(i).poi(1)
            point_2     = geom.iniL(i).poi(2)
            pos_1(1:3)  = geom.iniP(point_1).pos(1:3)
            pos_2(1:3)  = geom.iniP(point_2).pos(1:3)
            pos_c(1:3)  = (pos_1(1:3) + pos_2(1:3)) / 2.0d0
            vec(1:3)    = Normalize(pos_2(1:3) - pos_1(1:3))
            pos_c2(1:3) = (6.0d0*pos_2 + 4.0d0*pos_1) / (6.0d0 + 4.0d0)
            pos_c1(1:3) = (6.0d0*pos_1 + 4.0d0*pos_2) / (6.0d0 + 4.0d0)

            ! Write line index in blue
            write(301, "(a$   )"), ".cmov "
            write(301, "(3f9.3)"), pos_c(1:3) + 0.4d0
            write(301, "(a    )"), ".color blue"
            write(301, "(a    )"), ".font Helvetica 12 bold"
            write(301, "(i7   )"), i

            ! Draw the direction of line connectivity in blue
            write(301, "(a$    )"), ".arrow "
            write(301, "(3f8.2$)"), pos_c(1:3)
            write(301, "(3f8.2$)"), pos_c(1:3) + 1.5d0*vec(1:3)
            write(301, "(2f8.2 )"), 0.25d0, 0.5d0

            ! Write the positive(left) sign at point 1 from line connectivity
            if(geom.iniL(i).neiL(1, 1) /= -1) then
                vec(1:3) = geom.iniP(geom.iniL(i).neiP(1, 1)).pos - pos_c1
                vec(1:3) = Normalize(vec)
                write(301, "(a$   )"), ".cmov "
                write(301, "(3f9.3)"), pos_c1(1:3) + 2.0d0*vec(1:3)
                write(301, "(a    )"), ".color red"
                write(301, "(i7   )"), geom.iniL(i).neiL(1, 1)
            end if

            ! Write the negative(right) sign at point 1 from line connectivity
            if(geom.iniL(i).neiL(1, 2) /= -1) then
                vec(1:3) = geom.iniP(geom.iniL(i).neiP(1, 2)).pos - pos_c1
                vec(1:3) = Normalize(vec)
                write(301, "(a$   )"), ".cmov "
                write(301, "(3f9.3)"), pos_c1(1:3) + 2.0d0*vec(1:3)
                write(301, "(a    )"), ".color red"
                write(301, "(i7   )"), geom.iniL(i).neiL(1, 2)
            end if

            ! Write the positive(left) sign at point 2 from line connectivity
            if(geom.iniL(i).neiL(2, 1) /= -1) then
                vec(1:3) = geom.iniP(geom.iniL(i).neiP(2, 1)).pos - pos_c2
                vec(1:3) = Normalize(vec)
                write(301, "(a$   )"), ".cmov "
                write(301, "(3f9.3)"), pos_c2(1:3) + 2.0d0*vec(1:3)
                write(301, "(a    )"), ".color red"
                write(301, "(i7   )"), geom.iniL(i).neiL(2, 1)
            end if

            ! Write the negative(right) sign at point 2 from line connectivity
            if(geom.iniL(i).neiL(2, 2) /= -1) then
                vec(1:3) = geom.iniP(geom.iniL(i).neiP(2, 2)).pos - pos_c2
                vec(1:3) = Normalize(vec)
                write(301, "(a$   )"), ".cmov "
                write(301, "(3f9.3)"), pos_c2(1:3) + 2.0d0*vec(1:3)
                write(301, "(a    )"), ".color red"
                write(301, "(i7   )"), geom.iniL(i).neiL(2, 2)
            end if
        end do

        ! Write point index
        do i = 1, geom.n_iniP
            write(301, "(a$   )"), ".cmov "
            write(301, "(3f9.3)"), geom.iniP(i).pos(1:3) + 1.0d0
            write(301, "(a    )"), ".color black"
            write(301, "(i7   )"), i
        end do

        ! Write face index
        do i = 1, geom.n_face
            pos_c(1:3) = 0.0d0
            do j = 1, geom.face(i).n_poi
                pos_c(1:3) = pos_c + geom.iniP(geom.face(i).poi(j)).pos
            end do
            pos_c(1:3) = pos_c / dble(geom.face(i).n_poi)

            write(301, "(a$   )"), ".cmov "
            write(301, "(3f9.3)"), pos_c(1:3) + 1.0d0
            write(301, "(a    )"), ".color dark green"
            write(301, "(i7   )"), i
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(301, "(a)"), ".translate 0.0 0.0 0.0"
        write(301, "(a)"), ".scale 0.5"
        write(301, "(a)"), ".color grey"
        write(301, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(301, "(a)"), ".color red"             ! x-axis
        write(301, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(301, "(a)"), ".color blue"            ! y-axis
        write(301, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(301, "(a)"), ".color yellow"          ! z-axis
        write(301, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if

    close(unit=301)
end subroutine ModGeo_Chimera_Check_Geometry

! ---------------------------------------------------------------------------------------

! Set arm junction data
subroutine ModGeo_Set_Junction_Data(geom, bound)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    type(ListJunc), pointer :: junc, ptr
    integer :: i, j, k, count, poi_cur, poi_com

    ! Nullify linked list for junction data
    nullify(junc)
    nullify(ptr)

    ! --------------------------------------------------
    ! Find junctions and allocate linked list
    ! --------------------------------------------------
    do i = 1, geom.n_iniP

        count   = 0
        poi_cur = i

        ! Allocate ptr linked list
        allocate(ptr)

        ! Loop for comparing line
        do j = 1, geom.n_iniL
            do k = 1, 2

                poi_com = geom.iniL(j).poi(k)
                if(poi_cur == poi_com) then
                    count          = count + 1
                    ptr%n_arm      = count
                    ptr%cnL(count) = j
                    ptr%poi_c      = poi_cur
                end if
            end do
        end do

        ! If junction is connected to more than 2 arms
        if(count >= 2) then
            junc => List_Insert_Junc(junc, ptr)
        else
            write(0, "(a$)"), "Error - Not closed geometry : "
            write(0, "(a )"), "ModGeo_Set_Junction_Data"
            stop
        end if
    end do

    ! Allocate global junc data structure
    bound.n_junc = List_Count_Junc(junc)
    allocate(bound.junc(bound.n_junc))

    ! Copy from linked list to global vector
    ptr => junc
    do i = 1, bound.n_junc

        bound.junc(bound.n_junc-i+1).n_arm = ptr%n_arm
        bound.junc(bound.n_junc-i+1).poi_c = ptr%poi_c

        allocate(bound.junc(bound.n_junc-i+1).iniL(ptr%n_arm))
        allocate(bound.junc(bound.n_junc-i+1).modP(ptr%n_arm))
        allocate(bound.junc(bound.n_junc-i+1).croP(ptr%n_arm, geom.n_sec))
        allocate(bound.junc(bound.n_junc-i+1).node(ptr%n_arm, geom.n_sec))
        allocate(bound.junc(bound.n_junc-i+1).conn(ptr%n_arm*geom.n_sec, 2))
        allocate(bound.junc(bound.n_junc-i+1).type_conn(ptr%n_arm*geom.n_sec))

        do j = 1, ptr%n_arm
            bound.junc(bound.n_junc-i+1).iniL(j)    = ptr%cnL(j)    ! Initial line
            bound.junc(bound.n_junc-i+1).modP(j)    = -1            ! Modified point
            bound.junc(bound.n_junc-i+1).croP(j, :) = -1            ! Sectional point
            bound.junc(bound.n_junc-i+1).node(j, :) = -1            ! Node
        end do

        do j = 1, geom.n_sec*ptr%n_arm
            bound.junc(bound.n_junc-i+1).conn(j, :)   = -1      ! Sectional connection
            bound.junc(bound.n_junc-i+1).type_conn(j) = -1      ! Type of sectional connection
        end do

        ptr => ptr%next
    end do

    ! Delete linked list allocated
    call List_Delete_Junc(junc)
    !call List_Delete_Junc(ptr)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.3. Construct arm junction data"
    end do
    do i = 1, bound.n_junc
        ! 1 th junction : 3-arms -> center point # 1, lines # : (1, 2, 3)
        write(11, "(i20, a$)"), i, " th junc : "
        write(11, "(i4,  a$)"), bound.junc(i).n_arm, "-arms -> center point # : "
        write(11, "(i4,  a$)"), bound.junc(i).poi_c, ", lines # : ("
        do j = 1, bound.junc(i).n_arm - 1
            write(11, "(i4, a$)"), bound.junc(i).iniL(j), ", "
        end do
        write(11, "(i4, a)"), bound.junc(i).iniL(bound.junc(i).n_arm), ")"
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Set_Junction_Data

! ---------------------------------------------------------------------------------------

! Set local coordinate system on each edge
subroutine ModGeo_Set_Local_Coorindate(geom)
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.4. Set local coordinate system on initial edges"
    end do

    ! Set local vecotrs on each line
    do i = 1, geom.n_iniL

        ! Set local coordinate system
        geom.iniL(i).t(:,:) = ModGeo_Set_Local_Vectors(geom, i)

        ! Print detailed information
        ! 1 th line : t1 -> [**, **, **], t2 -> [**, **, **], t3 -> [**, **, **]
        write(11, "(i20,  a$)"), i, " th line : t1 -> [ "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(1, 1), ", "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(1, 2), ", "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(1, 3), " ], t2 -> [ "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(2, 1), ", "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(2, 2), ", "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(2, 3), " ], t3 -> [ "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(3, 1), ", "
        write(11, "(f5.2, a$)"), geom.iniL(i).t(3, 2), ", "
        write(11, "(f5.2, a )"), geom.iniL(i).t(3, 3), " ]"
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Set_Local_Coorindate

! ---------------------------------------------------------------------------------------

! Set local vecotrs on each edge
function ModGeo_Set_Local_Vectors(geom, line) result(local)
    type(GeomType), intent(in) :: geom
    integer,        intent(in) :: line

    double precision :: local(3, 3), pos_1(3), pos_2(3), vec_face1(3), vec_face2(3)
    double precision :: vec_a(3), vec_b(3), pos_c(3)
    integer :: i, poi_1, poi_2, face1, face2

    ! ==================================================
    !
    ! Set first local vector, t1 - along the direction of the connectivity
    !
    ! ==================================================
    poi_1 = geom.iniL(line).poi(1)
    poi_2 = geom.iniL(line).poi(2)
    pos_1(1:3) = geom.iniP(poi_1).pos(1:3)
    pos_2(1:3) = geom.iniP(poi_2).pos(1:3)
    local(1,:) = Normalize(pos_2 - pos_1)

    ! ==================================================
    !
    ! Set second local vector, t2
    !
    ! ==================================================
    face1 = geom.iniL(line).neiF(1)
    if(face1 /= -1) then

        ! Find outward vector of the face 1
        vec_a(1:3) = geom.iniP(geom.face(face1).poi(2)).pos - geom.iniP(geom.face(face1).poi(1)).pos
        vec_b(1:3) = geom.iniP(geom.face(face1).poi(3)).pos - geom.iniP(geom.face(face1).poi(2)).pos
        vec_face1(1:3) = Normalize(Cross(vec_a, vec_b))
    else

        ! If the neighbor face 1 is boundary
        vec_face1(1:3) = 0.0d0
    end if

    face2 = geom.iniL(line).neiF(2)
    if(face2 /= -1) then

        ! Find outward vector of the face 2
        vec_a(1:3) = geom.iniP(geom.face(face2).poi(2)).pos - geom.iniP(geom.face(face2).poi(1)).pos
        vec_b(1:3) = geom.iniP(geom.face(face2).poi(3)).pos - geom.iniP(geom.face(face2).poi(2)).pos
        vec_face2(1:3) = Normalize(Cross(vec_a, vec_b))
    else

        ! If the neighbor face 2 is boundary
        vec_face2(1:3) = 0.0d0
    end if

    ! Set second local vector, t2
    if(face1 == -1 .or. face2 == -1) then
        local(2,:) = (vec_face1 + vec_face2)
    else
        local(2,:) = 0.5d0*(vec_face1 + vec_face2)
    end if
    local(2,:) = Normalize(local(2,:))

    ! ==================================================
    !
    ! Set third local vector, t3
    !
    ! ==================================================
    local(3,:) = Cross(local(1,:), local(2,:))
    local(3,:) = Normalize(local(3,:))

    ! Check local vector
    if(dabs(Norm(local(1, 1:3)) - 1.0d0) > eps) then
        write(0, "(a$)"), "Error - The t1 local vector was not defined : "
        write(0, "(a$)"), "ModGeo_Set_Local_Vectors"
        stop
    else if(dabs(Norm(local(2, 1:3)) - 1.0d0) > eps) then
        write(0, "(a$)"), "Error - The t2 local vector was not defined : "
        write(0, "(a$)"), "ModGeo_Set_Local_Vectors"
        stop
    else if(dabs(Norm(local(3, 1:3)) - 1.0d0) > eps) then
        write(0, "(a$)"), "Error - The t3 local vector was not defined : "
        write(0, "(a$)"), "ModGeo_Set_Local_Vectors"
        stop
    end if
end function ModGeo_Set_Local_Vectors

! ---------------------------------------------------------------------------------------

! Write initial geometry with the local coordinate system
subroutine ModGeo_Chimera_Init_Geometry_Local(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: pos_1(3), pos_2(3), pos_c(3)
    logical :: f_axis, f_info
    integer :: i
    character(200) :: path

    if(para_write_302 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_302_info

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=302, file=trim(path)//"_init_geo_local.bild", form="formatted")

    ! Write edges
    write(302, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos
        write(302, "(a$    )"), ".cylinder "
        write(302, "(3f9.3$)"), pos_1(1:3)
        write(302, "(3f9.3$)"), pos_2(1:3)
        write(302, "(1f9.3 )"), 0.2d0
    end do

    ! Write point as sphere
    write(302, "(a)"), ".color red"
    do i = 1, geom.n_iniP
        write(302, "(a$    )"), ".sphere "
        write(302, "(3f9.3$)"), geom.iniP(i).pos(1:3)
        write(302, "(1f9.3 )"), 0.5d0
    end do

    ! Draw local vectors
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(302, "(a     )"), ".color red"     ! first vector
        write(302, "(a$    )"), ".arrow "
        write(302, "(3f8.2$)"), pos_c(1:3)
        write(302, "(3f8.2$)"), pos_c + geom.iniL(i).t(1,:) * 1.8d0
        write(302, "(3f8.2 )"), 0.22d0, 0.44d0, 0.6d0

        write(302, "(a     )"), ".color blue"    ! second vector
        write(302, "(a$    )"), ".arrow "
        write(302, "(3f8.2$)"), pos_c(1:3)
        write(302, "(3f8.2$)"), pos_c + geom.iniL(i).t(2,:) * 1.8d0
        write(302, "(3f8.2 )"), 0.2d0, 0.4d0, 0.6d0

        write(302, "(a     )"), ".color yellow"  ! third vector
        write(302, "(a$    )"), ".arrow "
        write(302, "(3f8.2$)"), pos_c(1:3)
        write(302, "(3f8.2$)"), pos_c + geom.iniL(i).t(3,:) * 1.8d0
        write(302, "(3f8.2 )"), 0.2d0, 0.4d0, 0.6d0
    end do

    ! Information on the index of the each line and point
    if(f_info == .true.) then

        ! For edge index
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos
            pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

            write(302, "(a$   )"), ".cmov "
            write(302, "(3f9.3)"), pos_c(1:3) + 0.4d0
            write(302, "(a    )"), ".color black"
            write(302, "(a    )"), ".font Helvetica 12 bold"
            write(302, "(i7   )"), i
        end do

        ! For point index
        do i = 1, geom.n_iniP
            write(302, "(a$   )"), ".cmov "
            write(302, "(3f9.3)"), geom.iniP(i).pos(1:3) + 0.4d0
            write(302, "(a    )"), ".color red"
            write(302, "(a    )"), ".font Helvetica 12 bold"
            write(302, "(i7   )"), i
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(302, "(a)"), ".translate 0.0 0.0 0.0"
        write(302, "(a)"), ".scale 0.5"
        write(302, "(a)"), ".color grey"
        write(302, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(302, "(a)"), ".color red"             ! x-axis
        write(302, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(302, "(a)"), ".color blue"            ! y-axis
        write(302, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(302, "(a)"), ".color yellow"          ! z-axis
        write(302, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=302)

    ! ==================================================
    !
    ! Write the file for Tecplot
    !
    ! ==================================================
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=302, file=trim(path)//"_init_geo_local.dat", form="formatted")

    write(302, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(302, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(302, "(a$)"), 'ZONE F = FEPOINT'
    write(302, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniP)))
    write(302, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a )"), ', ET=LINESEG'

    ! Write vertices
    do i = 1, geom.n_iniP
        write(302, "(3f9.3$)"), geom.iniP(i).pos(1:3)
        write(302, "(3f9.3 )"), 1.0d0, 1.0d0, 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(302, "(2i7)"), geom.iniL(i).poi(1), geom.iniL(i).poi(2)
    end do

    write(302, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(302, "(a$)"), 'ZONE F = FEPOINT'
    write(302, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(302, "(3f9.3$)"), pos_c(1:3)
        write(302, "(3f9.3 )"), geom.iniL(i).t(1, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(302, "(2i7)"), i, i
    end do

    write(302, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(302, "(a$)"), 'ZONE F = FEPOINT'
    write(302, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(302, "(3f9.3$)"), pos_c(1:3)
        write(302, "(3f9.3 )"), geom.iniL(i).t(2, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(302, "(2i7)"), i, i
    end do

    write(302, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(302, "(a$)"), 'ZONE F = FEPOINT'
    write(302, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(302, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(302, "(3f9.3$)"), pos_c(1:3)
        write(302, "(3f9.3 )"), geom.iniL(i).t(3, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(302, "(2i7)"), i, i
    end do

    close(unit=302)
end subroutine ModGeo_Chimera_Init_Geometry_Local

! ---------------------------------------------------------------------------------------

! Seperate the line from the vertex without off-set distance
subroutine ModGeo_Seperate_Line(geom, bound)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    integer, allocatable, dimension(:) :: count_arm
    integer :: i, j, poi_1, poi_2, poi_c

    ! Set the number of the seperated points and lines
    geom.n_modP = 2 * geom.n_iniL
    geom.n_iniL = geom.n_iniL

    ! Reallocate global point data
    allocate(geom.modP(geom.n_modP))

    ! Initialize the number of the arm junctions
    allocate(count_arm(bound.n_junc))
    count_arm(1:bound.n_junc) = 0

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.5. Seperate edges from the vertex"
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial vertex   : "//trim(adjustl(Int2Str(geom.n_iniP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial edges    : "//trim(adjustl(Int2Str(geom.n_iniL)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of seperated points : "//trim(adjustl(Int2Str(geom.n_modP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of seperated lines  : "//trim(adjustl(Int2Str(geom.n_iniL)))
        write(i, "(a)")
    end do

    ! Set connectivity of the each seperate line
    do i = 1, geom.n_iniL

        ! Point of initial line
        poi_1 = geom.iniL(i).poi(1)
        poi_2 = geom.iniL(i).poi(2)

        ! Set the position vector for seperated points
        geom.modP(2*i-1).pos(1:3) = geom.iniP(poi_1).pos(1:3)
        geom.modP(2*i-0).pos(1:3) = geom.iniP(poi_2).pos(1:3)

        ! Set new connectivity for new points
        geom.iniL(i).poi(1) = 2 * i - 1
        geom.iniL(i).poi(2) = 2 * i

        ! Update junction data with the seperated points
        do j = 1, bound.n_junc

            ! The original point at the junction
            poi_c = bound.junc(j).poi_c

            if(poi_1 == poi_c) then
                count_arm(j) = count_arm(j) + 1
                bound.junc(j).modP(count_arm(j)) = geom.iniL(i).poi(1)
            else if(poi_2 == poi_c) then
                count_arm(j) = count_arm(j) + 1
                bound.junc(j).modP(count_arm(j)) = geom.iniL(i).poi(2)
            end if
        end do
    end do

    ! Deallocation
    deallocate(count_arm)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.6. Update junction data with the seperated points"
    end do

    do i = 1, bound.n_junc
        ! 1 th junction : 3-arms -> mod points # : (1, 2, 3)
        write(11, "(i20, a$)"), i, " th junc : "
        write(11, "(i3,  a$)"), bound.junc(i).n_arm, "-arms -> seperated points - ("
        do j = 1, bound.junc(i).n_arm - 1
            write(11, "(i4, a$)"), bound.junc(i).modP(j), ", "
        end do
        write(11, "(i4, a)"), bound.junc(i).modP(j), ")"
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Seperate_Line

! ---------------------------------------------------------------------------------------

! Find the factor to scale the geometry
function ModGeo_Find_Scale_Factor(prob, geom, bound) result(scale)
    type(ProbType),  intent(in)    :: prob
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    double precision, allocatable, dimension (:,:) :: pos_modP
    double precision :: tot_ang, ref_ang, ang, width, factor
    double precision :: scale, length, min_length, diff, cur_length
    double precision :: pos_cur(3), pos_opp(3), vec_a(3), vec_b(3)
    integer :: i, j, poi_cur, poi_1, poi_2

    ! Scale up to avoid small edge length after modification
    !call ModGeo_Set_Scale_Geometry(geom, 100.0d0/para_init_scale)

    ! To save the position of the seperated points with off-set distance
    allocate(pos_modP(geom.n_modP, 3))

    ! Set the angle between two adjacent lines sharing the same faces
    call ModGeo_Set_Angle_Junction(geom, bound)

    ! Set the width of cross-section
    width = ModGeo_Set_Width_Section(geom)

    ! Loop for all junction to pre-calculate the edge length
    do i = 1, bound.n_junc

        ! Find reference and total angle
        ref_ang = bound.junc(i).ref_ang
        tot_ang = bound.junc(i).tot_ang

        if(prob.type_geo == "open") tot_ang = 2.0d0 * pi

        ! Total angle : Reference angle = 2 pi : ang
        ang = ref_ang * (2.0d0 * pi / tot_ang)

        ! 0.2    -------- 60
        ! factor -------- ang-ref_ang
        ! 0.0    -------- 0
        if(geom.sec.n_col == 2) then
            factor = (0.20d0-0.0d0)*(ang-ref_ang)/Deg2Rad(60.0d0) + 0.0d0
        else
            factor = 0.0d0
        end if

        if(tot_ang <= 2.0d0 * pi) then
            ! Total angle : Inradius = 2 pi : bound.junc(i).gap
            bound.junc(i).gap = (width/2.0d0/dtan(ang/2.0d0) + factor) * (2.0d0 * pi / tot_ang)

            ! Find the apothem of a regular polygon at the junction
            ! a = 0.5*s/tan(180/n), https://en.wikipedia.org/wiki/Apothem
            !bound.junc(i).gap = (0.5d0 * width / dtan(pi/dble(bound.junc(i).n_arm))) * (ang/ref_ang)
        else
            bound.junc(i).gap = width/2.0d0/dtan(ang/2.0d0)
        end if

        ! Find modified position due to off-set distance
        do j = 1, bound.junc(i).n_arm

            ! The seperated points at the junction
            poi_cur = bound.junc(i).modP(j)
            pos_cur = geom.modP(poi_cur).pos(1:3)

            poi_1 = geom.iniL(bound.junc(i).iniL(j)).poi(1)
            poi_2 = geom.iniL(bound.junc(i).iniL(j)).poi(2)

            ! Find point that is opposite to poi_cur
            if(poi_1 == poi_cur) then
                pos_opp(1:3) = geom.modP(poi_2).pos(1:3)
            else if(poi_2 == poi_cur) then
                pos_opp(1:3) = geom.modP(poi_1).pos(1:3)
            end if

            ! Total length
            length = Norm(pos_opp(1:3) - pos_cur(1:3))

            ! Find modified position due to off-set distance
            vec_a(1:3) = bound.junc(i).gap * pos_opp(1:3)
            vec_b(1:3) = (length-bound.junc(i).gap) * pos_cur(1:3)
            pos_modP(poi_cur, 1:3) = (vec_a + vec_b)/length
        end do
    end do

    ! Find edge with minimum length
    do i = 1, geom.n_iniL

        ! Find modified point from the seperated line
        poi_1 = geom.iniL(i).poi(1)
        poi_2 = geom.iniL(i).poi(2)

        ! Length of the modified edge
        length = Norm(pos_modP(poi_1,:) - pos_modP(poi_2,:))

        ! Find modified edge with minimum length
        if(i == 1 .or. min_length > length) then
            min_length = length
            cur_length = Norm(geom.modP(poi_1).pos - geom.modP(poi_2).pos)
            diff       = cur_length - min_length
        end if
    end do

    ! Desired length
    if(geom.sec.types == "square") then
        length = diff + para_dist_bp * (dble(prob.n_bp_edge-1))
    else if(geom.sec.types == "honeycomb") then
        length = diff + para_dist_bp * (dble(prob.n_bp_edge-1-1))
    end if
    scale = length / cur_length

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.7. Find the scale factor to adjust polyhedra size"
        call Space(i, 11)
        write(i, "(a)"), "* The minumum edge length       : "//trim(adjustl(Int2Str(prob.n_bp_edge)))
        call Space(i, 11)
        write(i, "(a)"), "* Scale factor to adjust size   : "//trim(adjustl(Dble2Str(scale)))
        write(i, "(a)")
    end do

    ! Deallocate
    deallocate(pos_modP)
end function ModGeo_Find_Scale_Factor

! ---------------------------------------------------------------------------------------

! Set the angle between two adjacent lines sharing the same faces
subroutine ModGeo_Set_Angle_Junction(geom, bound)
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(inout) :: bound

    type :: JuncType
        integer :: n_arm    ! The number of arms
        double precision, allocatable, dimension(:) :: ang
    end type JuncType

    type(JuncType), allocatable, dimension(:) :: junc
    double precision :: pos_cur(3), pos_pre(3), pos_next(3)
    double precision :: tot_ang, ref_ang, max_ang, min_ang, ang
    double precision :: vec_1(3), vec_2(3)
    integer :: i, j, poi, poi_c

    ! Allocate and initialize junc memory
    allocate(junc(bound.n_junc))
    do i = 1, bound.n_junc
        junc(i).n_arm = bound.junc(i).n_arm
        allocate(junc(i).ang(junc(i).n_arm))

        ! Initialize memory
        junc(i).n_arm = 0
        do j = 1, junc(i).n_arm
            junc(i).ang(j) = 0.0d0
        end do
    end do

    ! Find angle at vertex
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi

            ! Find current point in the i-th face
            poi_c = geom.face(i).poi(j)
            pos_cur(:) = geom.iniP(poi_c).pos(1:3)

            ! Find previous point in the i-th face
            if(j == 1) then
                poi = geom.face(i).poi(geom.face(i).n_poi)
                pos_pre(1:3) = geom.iniP(poi).pos(1:3)
            else
                poi = geom.face(i).poi(j - 1)
                pos_pre(1:3) = geom.iniP(poi).pos(1:3)
            end if

            ! Find next point in the i-th face
            if(j == geom.face(i).n_poi) then
                poi = geom.face(i).poi(1)
                pos_next(1:3) = geom.iniP(poi).pos(1:3)
            else
                poi = geom.face(i).poi(j + 1)
                pos_next(1:3) = geom.iniP(poi).pos(1:3)
            end if

            vec_1(1:3) = pos_pre(1:3)  - pos_cur(1:3)
            vec_2(1:3) = pos_next(1:3) - pos_cur(1:3)

            ! Find angle between two vectors in 3D
            ang = datan2(Norm(Cross(vec_1, vec_2)), dot_product(vec_1, vec_2))

            junc(poi_c).n_arm = junc(poi_c).n_arm + 1
            junc(poi_c).ang(junc(poi_c).n_arm) = ang
        end do
    end do

    ! Print angle between two lines at the junction
    !do i = 1, bound.n_junc
    !    write(0, "(2i$)"), i, bound.junc(i).n_arm
    !    do j = 1, junc(i).n_arm
    !        write(0, "(f$)"), junc(i).ang(j) * 180.0d0 / pi
    !    end do
    !    write(0, "(a)")
    !end do

    ! Set total, maximum and minimum angle at the junction
    do i = 1, bound.n_junc

        ! Find total angle
        tot_ang = 0.0d0
        do j = 1, junc(i).n_arm
            tot_ang = tot_ang + junc(i).ang(j)
        end do
        bound.junc(i).tot_ang = tot_ang

        ! Find maximum and minimum angle
        max_ang = junc(i).ang(1)
        min_ang = junc(i).ang(1)
        do j = 2, junc(i).n_arm
            if(max_ang < junc(i).ang(j)) max_ang = junc(i).ang(j)
            if(min_ang > junc(i).ang(j)) min_ang = junc(i).ang(j)
        end do
        !write(0, "(i, 3f12.4)"), i, tot_ang*180.0d0/pi, max_ang*180.0d0/pi, min_ang*180.0d0/pi

        ! Set the reference anle in junction data
        if(para_junc_ang == "min" .or. para_junc_ang == "opt") then
            bound.junc(i).ref_ang = min_ang
        else if(para_junc_ang == "max") then
            bound.junc(i).ref_ang = max_ang
        else if(para_junc_ang == "ave") then
            bound.junc(i).ref_ang = tot_ang / dble(junc(i).n_arm)
        end if
    end do

    ! Deallocate memory
    do i = 1, bound.n_junc
        deallocate(junc(i).ang)
    end do
    deallocate(junc)
end subroutine ModGeo_Set_Angle_Junction

! ---------------------------------------------------------------------------------------

! Set width cross-section
function ModGeo_Set_Width_Section(geom) result(width)
    type(GeomType), intent(in) :: geom

    double precision :: width
    integer :: n_column

    ! Find the number of columns in terms of reference row
    if(geom.sec.types == "square") then
        n_column = geom.sec.ref_maxC - geom.sec.ref_minC + 1
    else if(geom.sec.types == "honeycomb") then
        if(geom.sec.ref_row == 1) then
            n_column = 2
        else
            ! 2 -> 3*1, 4 -> 3*2, 6 -> 3*3
            n_column = 3 * (geom.sec.ref_maxC - geom.sec.ref_minC + 1) / 2
        end if
    end if

    ! Find width of cross-section
    width = (2.0d0*para_rad_helix + para_gap_helix) * dble(n_column)
end function ModGeo_Set_Width_Section

! ---------------------------------------------------------------------------------------

! Set the geometric size with the scale factor
subroutine ModGeo_Set_Scale_Geometry(geom, scale)
    type(GeomType),   intent(inout) :: geom
    double precision, intent(in)    :: scale

    integer :: i

    ! Scale seperated lines with the input
    ! para_init_scale is initial scale factor from input
    do i = 1, geom.n_modP
        geom.modP(i).pos(1:3) = geom.modP(i).pos(1:3) * scale
    end do

    ! Rescale initial points
    do i = 1, geom.n_iniP
        geom.iniP(i).pos(1:3) = geom.iniP(i).pos(1:3) * scale
    end do
end subroutine ModGeo_Set_Scale_Geometry

! ---------------------------------------------------------------------------------------

! Write seperated geometry
subroutine ModGeo_Chimera_Sep_Geometry(prob, geom, mode)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    character(*),   intent(in) :: mode

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    integer :: i
    logical :: f_axis, f_info
    character(200) :: path

    if(para_write_303 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_303_info

    path = trim(prob.path_work1)//trim(prob.name_file)//"_"
    open(unit=303, file=trim(path)//trim(mode)//"_geo.bild", form="formatted")

    ! Write seperated points
    write(303, "(a)"), ".color red"
    do i = 1, geom.n_modP
        write(303, "(a$   )"), ".sphere "
        write(303, "(4f9.3)"), geom.modP(i).pos(1:3), 0.35d0
    end do

    ! Write seperated lines
    write(303, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL

        pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)

        write(303, "(a$   )"), ".cylinder "
        write(303, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.15d0
    end do

    ! Information on the index of points and lines
    if(f_info == .true.) then

        ! For point index
        do i = 1, geom.n_modP
            write(303, "(a$   )"), ".cmov "
            write(303, "(3f9.3)"), geom.modP(i).pos(1:3) + 0.4d0
            write(303, "(a    )"), ".color red"
            write(303, "(a    )"), ".font Helvetica 12 bold"
            write(303, "(i7   )"), i
        end do

        ! For line index
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0
            length     = Norm(pos_2 - pos_1)

            write(303, "(a$   )"), ".cmov "
            write(303, "(3f9.3)"), pos_c(1:3) + 0.4d0
            write(303, "(a    )"), ".color dark green"
            write(303, "(a    )"), ".font Helvetica 12 bold"
            !write(303, "(i7, a, f5.1, a)"), i, "(", length, ")"
            write(303, "(i7, a, f5.1, i, a)"), i, "("//&
                trim(adjustl(Dble2Str1(length)))//", "//&
                trim(adjustl(Int2Str(nint(length/para_dist_bp)+1)))//")"
        end do

        ! Draw three local vectors on each line
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

            write(303, "(a     )"), ".color red"     ! first vector
            write(303, "(a$    )"), ".arrow "
            write(303, "(3f8.2$)"), pos_c(1:3)
            write(303, "(3f8.2$)"), pos_c(1:3) + geom.iniL(i).t(1,1:3)*1.5d0
            write(303, "(3f8.2 )"), 0.18d0, 0.36d0, 0.6d0

            write(303, "(a     )"), ".color blue"    ! section vector
            write(303, "(a$    )"), ".arrow "
            write(303, "(3f8.2$)"), pos_c(1:3)
            write(303, "(3f8.2$)"), pos_c(1:3) + geom.iniL(i).t(2,1:3)*1.2d0
            write(303, "(3f8.2 )"), 0.18d0, 0.36d0, 0.6d0

            write(303, "(a     )"), ".color yellow"  ! third vector
            write(303, "(a$    )"), ".arrow "
            write(303, "(3f8.2$)"), pos_c(1:3)
            write(303, "(3f8.2$)"), pos_c(1:3) + geom.iniL(i).t(3,1:3)*1.2d0
            write(303, "(3f8.2 )"), 0.18d0, 0.36d0, 0.6d0
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(303, "(a)"), ".translate 0.0 0.0 0.0"
        write(303, "(a)"), ".scale 0.5"
        write(303, "(a)"), ".color grey"
        write(303, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(303, "(a)"), ".color red"             ! x-axis
        write(303, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(303, "(a)"), ".color blue"            ! y-axis
        write(303, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(303, "(a)"), ".color yellow"          ! z-axis
        write(303, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=303)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=303, file=trim(path)//"_sep_geo.dat", form="formatted")

    write(303, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(303, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(303, "(a$)"), 'ZONE F = FEPOINT'
    write(303, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_modP)))
    write(303, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a )"), ', ET=LINESEG'

    ! Write vertices
    do i = 1, geom.n_modP
        write(303, "(3f9.3$)"), geom.modP(i).pos(1:3)
        write(303, "(3f9.3 )"), 1.0d0, 1.0d0, 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(303, "(2i7)"), geom.iniL(i).poi(1), geom.iniL(i).poi(2)
    end do

    write(303, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(303, "(a$)"), 'ZONE F = FEPOINT'
    write(303, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(303, "(3f9.3$)"), pos_c(1:3)
        write(303, "(3f9.3 )"), geom.iniL(i).t(1, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(303, "(2i7)"), i, i
    end do

    write(303, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(303, "(a$)"), 'ZONE F = FEPOINT'
    write(303, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(303, "(3f9.3$)"), pos_c(1:3)
        write(303, "(3f9.3 )"), geom.iniL(i).t(2, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(303, "(2i7)"), i, i
    end do

    write(303, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(303, "(a$)"), 'ZONE F = FEPOINT'
    write(303, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(303, "(a )"), ', ET=LINESEG'

    ! Local coordinate system on edges
    do i = 1, geom.n_iniL
        pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(303, "(3f9.3$)"), pos_c(1:3)
        write(303, "(3f9.3 )"), geom.iniL(i).t(3, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(303, "(2i7)"), i, i
    end do

    close(unit=303)
end subroutine ModGeo_Chimera_Sep_Geometry

! ---------------------------------------------------------------------------------------

! Set the gap distance between the junction and end of edges
subroutine ModGeo_Set_Gap_Junction(geom, bound)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    double precision :: dist_gap, length, remain
    double precision :: pos_cur(3), pos_opp(3), vec_a(3), vec_b(3)
    integer :: i, j, poi_cur, poi_1, poi_2, n_bp

    double precision :: pos_c(3), pos_1(3), pos_2(3)
    double precision :: ref_length, len_ini, len_mod, ratio
    integer :: ref_edge, iniP1, iniP2, modP1, modP2

    ! Loop for junction to modify gap distance
    do i = 1, bound.n_junc

        ! Find gap distance
        dist_gap = bound.junc(i).gap

        ! Modify the edge length
        do j = 1, bound.junc(i).n_arm

            ! Position vector of arm junction point
            poi_cur = bound.junc(i).modP(j)
            pos_cur = geom.modP(poi_cur).pos(1:3)

            poi_1 = geom.iniL(bound.junc(i).iniL(j)).poi(1)
            poi_2 = geom.iniL(bound.junc(i).iniL(j)).poi(2)

            ! Find point that is opposite to poi_cur
            if(poi_1 == poi_cur) then
                pos_opp(1:3) = geom.modP(poi_2).pos(1:3)
            else if(poi_2 == poi_cur) then
                pos_opp(1:3) = geom.modP(poi_1).pos(1:3)
            end if

            ! Original edge length
            length = Norm(pos_opp(1:3) - pos_cur(1:3))

            ! Find modified position vector
            vec_a(1:3) = dist_gap*pos_opp(1:3)
            vec_b(1:3) = (length-dist_gap)*pos_cur(1:3)
            geom.modP(poi_cur).pos(1:3) = (vec_a + vec_b)/length
        end do
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "2.8. Update edge length with the scale factor"
    end do

    do i = 1, geom.n_iniL
        poi_1  = geom.iniL(i).poi(1)
        poi_2  = geom.iniL(i).poi(2)
        length = Norm(geom.modP(poi_1).pos - geom.modP(poi_2).pos)
        n_bp   = nint(length / para_dist_bp) + 1
        remain = dmod(length, para_dist_bp)

        if(dmod(length, para_dist_bp) >= para_dist_bp/2.0d0) then
            remain = 0.0d0
        end if

        write(11, "(i20,  a$)"), i,      " th edge -> edge length : "
        write(11, "(f8.3, a$)"), length, " [nm] -> "
        write(11, "(i4,   a$)"), n_bp,   " BP, remaining length : "
        write(11, "(f8.3, a )"), remain, " [nm]"
    end do
    write(0, "(a)"); write(11, "(a)")
end subroutine ModGeo_Set_Gap_Junction

! ---------------------------------------------------------------------------------------

! Set constant modified edge ratio based on original geometry
subroutine ModGeo_Set_Const_Geometric_Ratio(geom)
    type(GeomType), intent(inout) :: geom

    double precision :: length, vec_a(3), vec_b(3), pos_c(3), pos_1(3), pos_2(3)
    double precision :: ref_length, len_ini, len_mod, len_new, len_ref, ratio, magic
    integer :: i, ref_edge, iniP1, iniP2, modP1, modP2

    ! Calculate magic depending on types of section
    if(geom.sec.types == "square")    magic = 0.0d0
    if(geom.sec.types == "honeycomb") magic = 0.34d0

    ! Find reference edge
    do i = 1, geom.n_iniL
        iniP1  = geom.iniL(i).poi(1)
        iniP2  = geom.iniL(i).poi(2)
        length = Norm(geom.modP(iniP1).pos - geom.modP(iniP2).pos)

        ! Find reference edge
        if(i == 1 .or. ref_length > length) then
            ref_length = length
            ref_edge   = i
        end if
    end do

    ! Find ratio divied by initial edge length
    modP1   = geom.iniL(ref_edge).poi(1)
    modP2   = geom.iniL(ref_edge).poi(2)
    iniP1   = geom.iniL(ref_edge).iniP(1)
    iniP2   = geom.iniL(ref_edge).iniP(2)
    len_ini = Norm(geom.iniP(iniP1).pos - geom.iniP(iniP2).pos)
    len_mod = Norm(geom.modP(modP1).pos - geom.modP(modP2).pos)
    len_ref = len_ini
    ratio   = (len_mod + magic) / len_ini

    ! Rescale the edge with ratio
    do i = 1, geom.n_iniL
        iniP1 = geom.iniL(i).iniP(1)
        iniP2 = geom.iniL(i).iniP(2)
        modP1 = geom.iniL(i).poi(1)
        modP2 = geom.iniL(i).poi(2)

        pos_1(1:3) = geom.iniP(iniP1).pos(1:3)
        pos_2(1:3) = geom.iniP(iniP2).pos(1:3)
        pos_c(1:3) = 0.5d0 * (pos_1 + pos_2)
        len_ini    = Norm(geom.iniP(iniP1).pos - geom.iniP(iniP2).pos)
        len_mod    = Norm(geom.modP(modP1).pos - geom.modP(modP2).pos)

        vec_a(1:3) = Normalize(pos_1(1:3) - pos_c(1:3))
        vec_b(1:3) = Normalize(pos_2(1:3) - pos_c(1:3))

        ! Recalculate point position
        geom.modP(modP1).pos(1:3) = pos_c(1:3) + vec_a(1:3) * (ratio * len_ini / 2.0d0 - magic/2.0d0 + 0.05d0)
        geom.modP(modP2).pos(1:3) = pos_c(1:3) + vec_b(1:3) * (ratio * len_ini / 2.0d0 - magic/2.0d0 + 0.05d0)

        len_new = Norm(geom.modP(modP1).pos - geom.modP(modP2).pos)

        ! Print progress
        call space(0, 11)
        write(0, "(a$)"), trim(adjustl(Int2Str(i)))//" - th edge"
        write(0, "(a$)"), ", len(mesh) : "//trim(adjustl(Dble2Str(len_ini)))
        write(0, "(a$)"), " ["//trim(adjustl(Dble2Str(42.0d0*(len_ini/len_ref))))//"] "
        write(0, "(a$)"), ", len(mod) : "//trim(adjustl(Dble2Str(len_mod)))
        write(0, "(a$)"), " ["//trim(adjustl(Int2Str(nint(len_mod/0.34d0) + 1)))//"] "
        write(0, "(a$)"), ", len(new) : "//trim(adjustl(Dble2Str(len_new)))
        write(0, "(a )"), " ["//trim(adjustl(Int2Str(nint(len_new/0.34d0) + 1)))//"] "
    end do
end subroutine ModGeo_Set_Const_Geometric_Ratio

! ---------------------------------------------------------------------------------------

end module ModGeo