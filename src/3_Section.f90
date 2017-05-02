!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Section
!
!                                                                    Updated : 2017/03/27
!
! Comments: This module is to generate multiple lines.
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Section

    use Data_Prob
    use Data_Geom
    use Data_Bound

    use Para
    use Math
    use Mani

    implicit none

    public  Section_Generation
    public  Section_Connection_Scaf
    public  Section_Connection_Stap

    private Section_Set_Sectional_Data
    private Section_Generate_Section_Geometry
    private Section_Generate_Square
    private Section_Generate_Honeycomb
    private Section_Get_Position
    private Section_Get_Shape_Function
    private Section_Get_Director
    private Section_Get_Parameter
    private Section_Reset_Local_Coordinate
    private Section_Chimera_Cro_Geometry

contains

! ---------------------------------------------------------------------------------------

! Cross-section generation
subroutine Section_Generation(prob, geom, bound)
    type(ProbType),  intent(in)    :: prob
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    double precision :: i

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |                 3. Build cross-sectional edges                     |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
    end do

    ! Set sectional data, croP and croL
    call Section_Set_Sectional_Data(geom, bound)

    ! Generate cross-sectional geometry
    call Section_Generate_Section_Geometry(geom)

    ! Write cross-sectional geometry
    call Section_Chimera_Cro_Geometry(prob, geom)
end subroutine Section_Generation

! ---------------------------------------------------------------------------

! Determine the section connection for scaffold strand
function Section_Connection_Scaf(geom, sec_cur, sec_com, bp_id) result(b_connect)
    type(GeomType), intent(in) :: geom
    integer, intent(in) :: sec_cur
    integer, intent(in) :: sec_com
    integer, intent(in) :: bp_id

    integer :: i, j, bp, row_cur, row_com, col_cur, col_com
    logical :: b_connect

    if(geom.sec.types == "square") then
        bp = bp_id + para_start_bp_ID - 1
        if(bp < 0) bp = 32 + bp
        bp = mod(bp, 32)
    else if(geom.sec.types == "honeycomb") then
        bp = bp_id + para_start_bp_ID - 1
        if(bp < 0) bp = 21 + bp
        bp = mod(bp, 21)
    end if

    ! initialize return value
    b_connect = .false.

    ! find section position
    row_cur = geom.sec.posR(sec_cur + 1)
    col_cur = geom.sec.posC(sec_cur + 1)
    row_com = geom.sec.posR(sec_com + 1)
    col_com = geom.sec.posC(sec_com + 1)

    !write(0, "(i10, a, 2i)"), sec_cur, " th : section(current),   row and col # : ", row_cur, col_cur
    !write(0, "(i10, a, 2i)"), sec_com, " th : section(comparing), row and col # : ", row_com, col_com
    !write(0, "(a)")

    ! determine whether the section connects or not
    if(geom.sec.types == "square") then

        ! --------------------------------------------------
        !
        ! for square lattice, possible crossovers in scaffold
        !
        ! --------------------------------------------------
        if(row_cur == row_com .and. col_cur == col_com - 1) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then

                ! left(cur) to right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16 .or. bp==26 .or. bp==27) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Left  ---> Right"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then

                ! left(cur) to right(com) connection
                !  1  ------->  0        ===>  row : same
                ! cur(odd)    com(even)        col : cur < com
                if(bp==0 .or. bp==10 .or. bp==11 .or. bp==20 .or. bp==21 .or. bp==31) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Left  ---> Right"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com .and. col_cur == col_com + 1) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then

                ! right(cur) to left(com) connection
                !  1  <-------  0       ===>  row : same
                ! com(odd)     cur(even)      col : cur > com
                if(bp==0 .or. bp==10 .or. bp==11 .or. bp==20 .or. bp==21 .or. bp==31) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Right ---> Left"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then

                ! right(cur) to left(com) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16 .or. bp==26 .or. bp==27) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Right ---> Left"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com + 1 .and. col_cur == col_com) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then

                ! up(cur) to bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp==7 .or. bp==8 .or. bp==18 .or. bp==19 .or. bp==28 .or. bp==29) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Top   ---> Bottom"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then

                ! up(cur) to bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp==2 .or. bp==3 .or. bp==12 .or. bp==13 .or. bp==23 .or. bp==24) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Top   ---> Bottom"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com - 1 .and. col_cur == col_com) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then

                ! bottom(cur) to up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp==2 .or. bp==3 .or. bp==12 .or. bp==13 .or. bp==23 .or. bp==24) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Bottom ---> Top"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then

                ! bottom(cur) to up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp==7 .or. bp==8 .or. bp==18 .or. bp==19 .or. bp==28 .or. bp==29) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Bottom ---> Top"
                    b_connect = .true.
                end if
            end if
        end if

    else if(geom.sec.types == "honeycomb") then

        ! --------------------------------------------------
        !
        ! For honeycomb lattice, possible crossovers in scaffold
        !
        ! --------------------------------------------------
        if(row_cur == row_com .and. col_cur == col_com - 1) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(cur) -> right(com) connection
                !  1  ------->  0        ===>  row : same
                ! cur(odd)    com(even)        col : cur < com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(cur) -> right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(cur) -> right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com .and. col_cur == col_com + 1) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(com) <- right(cur) connection
                !  1  <-------  0       ===>  row : same
                ! com(odd)     cur(even)      col : cur > com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(com) <- right(cur) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(com) <- right(cur) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com + 1 .and. col_cur == col_com) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com - 1 .and. col_cur == col_com) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp==1 .or. bp==2 .or. bp==11 .or. bp==12) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp==8 .or. bp==9 .or. bp==18 .or. bp==19) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp==4 .or. bp==5 .or. bp==15 .or. bp==16) then
                    b_connect = .true.
                end if
            end if
        end if
    end if
end function Section_Connection_Scaf

! ---------------------------------------------------------------------------

! determine the section connection for staple strand
function Section_Connection_Stap(geom, sec_cur, sec_com, bp_id) result(b_connect)
    type(GeomType), intent(in) :: geom
    integer,        intent(in) :: sec_cur, sec_com, bp_id

    integer :: i, j, bp, row_cur, row_com, col_cur, col_com
    logical :: b_connect

    if(geom.sec.types == "square") then
        bp = bp_id + para_start_bp_ID - 1
        if(bp < 0) bp = 32 + bp
        bp = mod(bp, 32)
    else if(geom.sec.types == "honeycomb") then
        bp = bp_id + para_start_bp_ID - 1
        if(bp < 0) bp = 21 + bp
        bp = mod(bp, 21)
    end if

    ! initialize boolean return variable
    b_connect = .false.

    ! find section position
    row_cur = geom.sec.posR(sec_cur + 1)
    col_cur = geom.sec.posC(sec_cur + 1)
    row_com = geom.sec.posR(sec_com + 1)
    col_com = geom.sec.posC(sec_com + 1)

    !write(0, "(i10, a, 2i)"), sec_cur, " th : section(current),   row and col # : ", row_cur, col_cur
    !write(0, "(i10, a, 2i)"), sec_com, " th : section(comparing), row and col # : ", row_com, col_com
    !write(0, "(a)")

    ! --------------------------------------------------
    ! determine whether the section connects or not
    ! --------------------------------------------------
    if(geom.sec.types == "square") then

        ! --------------------------------------------------
        !
        ! for square lattice, possible crossovers in staple
        !
        ! --------------------------------------------------
        if(row_cur == row_com .and. col_cur == col_com - 1) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(cur) to right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp == 0 .or. bp == 31) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Left  ---> Right"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(cur) to right(com) connection
                !  1  ------->  0        ===>  row : same
                ! cur(odd)    com(even)        col : cur < com
                if(bp == 15 .or. bp == 16) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Left  ---> Right"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com .and. col_cur == col_com + 1) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! right(cur) to left(com) connection
                !  1  <-------  0       ===>  row : same
                ! com(odd)     cur(even)      col : cur > com
                if(bp == 15 .or. bp == 16) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Right ---> Left"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! right(cur) to left(com) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp == 0 .or. bp == 31) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Right ---> Left"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com + 1 .and. col_cur == col_com) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) to bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp == 23 .or. bp == 24) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Top   ---> Bottom"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) to bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp == 7 .or. bp == 8) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Top   ---> Bottom"
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com - 1 .and. col_cur == col_com) then

            if(mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) to up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp == 7 .or. bp == 8) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Bottom ---> Top"
                    b_connect = .true.
                end if
            else if(mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) to up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp == 23 .or. bp == 24) then
                    !write(0, "(2i4, a)"), sec_cur, sec_com, ",   Bottom ---> Top"
                    b_connect = .true.
                end if
            end if
        end if

    else if(geom.sec.types == "honeycomb") then

        ! --------------------------------------------------
        !
        ! for honeycomb lattice, possible crossovers in staple
        !
        ! --------------------------------------------------
        if(row_cur == row_com .and. col_cur == col_com - 1) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(cur) -> right(com) connection
                !  1  ------->  0        ===>  row : same
                ! cur(odd)    com(even)        col : cur < com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(cur) -> right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(cur) -> right(com) connection
                !  0  ------->  1        ===>  row : same
                ! cur(even)    com(odd)        col : cur < com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com .and. col_cur == col_com + 1) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(com) <- right(cur) connection
                !  1  <-------  0       ===>  row : same
                ! com(odd)     cur(even)      col : cur > com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! left(com) <- right(cur) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! left(com) <- right(cur) connection
                !  0  <-------  1       ===>  row : same
                ! com(even)     cur(odd)      col : cur > com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com + 1 .and. col_cur == col_com) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! up(cur) -> bottom(com) connection
                !  1  cur(odd)
                !  ก้                ===>  col : same
                !  0  com(even)           row : cur > com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! up(cur) -> bottom(com) connection
                !  0  cur(even)
                !  ก้                ===>  col : same
                !  1  com(odd)            row : cur > com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            end if

        else if(row_cur == row_com - 1 .and. col_cur == col_com) then

            if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp == 6 .or. bp == 7) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp == 13 .or. bp == 14) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == -90 .and. mod(sec_cur, 2) == 1 .and. mod(sec_com, 2) == 0) then
                ! bottom(cur) -> up(com) connection
                !  0  com(even)
                !  ก่                ===>  col : same
                !  1  cur(odd)            row : cur < com
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 90 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            else if(geom.sec.dir == 150 .and. mod(sec_cur, 2) == 0 .and. mod(sec_com, 2) == 1) then
                ! bottom(cur) -> up(com) connection
                !  1  com(odd)
                !  ก่                ===>  col : same
                !  0  cur(even)           row : cur < com
                if(bp == 0 .or. bp == 20) then
                    b_connect = .true.
                end if
            end if
        end if
    end if
end function Section_Connection_Stap

! ---------------------------------------------------------------------------------------

! Set crosssectional data, croP and croL and update junction data
subroutine Section_Set_Sectional_Data(geom, bound)
    type(GeomType),  intent(inout) :: geom
    type(BoundType), intent(inout) :: bound

    integer :: i, j, m, n, point

    ! set total cross-sectional points and lines
    geom.n_croP = geom.n_sec * geom.n_modP
    geom.n_croL = geom.n_sec * geom.n_iniL

    ! allocate sectional data
    allocate(geom.croP(geom.n_croP))
    allocate(geom.croL(geom.n_croL))

    ! copy information from modified one to sectional one, loop for cross-section
    ! the sectional point, croP is derived from modified point, modP
    ! the sectional line, croL is derived from initial line, iniL that is seperated from the junction
    do i = 1, geom.n_sec

        ! for sectional point, croP
        ! initially, new positioin vector is assigned as the same to modified points
        do j = 1, geom.n_modP
            geom.croP(geom.n_modP*(i-1)+j).pos(1:3) = geom.modP(j).pos(1:3)

            ! update junction data
            do m = 1, bound.n_junc
                do n = 1, bound.junc(m).n_arm

                    ! if junction data is equal to modified point number, j
                    if(bound.junc(m).modP(n) == j) then
                        bound.junc(m).croP(n, i) = geom.n_modP*(i-1)+j
                    end if
                end do
            end do
        end do

        ! for sectional line, croL
        do j = 1, geom.n_iniL

            ! point connectivity
            geom.croL(geom.n_iniL*(i-1)+j).poi(1:2) = (i-1)*geom.n_modP + geom.iniL(j).poi(1:2)

            ! ID of section and initial line
            geom.croL(geom.n_iniL*(i-1)+j).sec  = geom.sec.id(i)
            geom.croL(geom.n_iniL*(i-1)+j).iniL = j

            do m = 1, 3
                geom.croL(geom.n_iniL*(i-1)+j).t(m, 1:3) = geom.iniL(j).t(m, 1:3)
            end do
        end do
    end do

    ! print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "3.1. Cross-section information"
        call Space(i, 11)
        write(i, "(a)"), "* Section type                   : "//trim(geom.sec.types)//" lattice"
        call Space(i, 11)
        write(i, "(a)"), "* The number of duplexes         : "//trim(adjustl(Int2Str(geom.n_sec)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of rows             : "//trim(adjustl(Int2Str(geom.sec.maxR-geom.sec.minR+1)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of columns          : "//trim(adjustl(Int2Str(geom.sec.maxC-geom.sec.minC+1)))
        call Space(i, 11)
        write(i, "(a)"), "* Reference row                  : "//trim(adjustl(Int2Str(geom.sec.ref_row)))
        call Space(i, 11)
        write(i, "(a)"), "* Reference min/max column       : "&
            //trim(adjustl(Int2Str(geom.sec.ref_minC)))//" / "//trim(adjustl(Int2Str(geom.sec.ref_maxC)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified points  : "//trim(adjustl(Int2Str(geom.n_modP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of modified edges   : "//trim(adjustl(Int2Str(geom.n_iniL)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of sectional points : "//trim(adjustl(Int2Str(geom.n_croP)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of sectional edges  : "//trim(adjustl(Int2Str(geom.n_croL)))
        write(i, "(a)")
        call Space(i, 6)
        write(i, "(a)"), "3.2. Update junctional data with cross-section infomation"
    end do

    ! print modified junction connectivity
    do m = 1, bound.n_junc

        write(11, "(i20$)"), m
        write(11, "(a$  )"), " th junc, the number of arms : "
        write(11, "(i5  )"), bound.junc(m).n_arm

        do n = 1, bound.junc(m).n_arm

            write(11, "(a40$)"), "Mod point # : "
            write(11, "(i5$ )"), bound.junc(m).modP(n)
            write(11, "(a$  )"), " -> cross point # : ( "

            do i = 1, geom.n_sec - 1
                write(11, "(i5, a$)"), bound.junc(m).croP(n, i), ", "
            end do

            write(11, "(i5, a)"), bound.junc(m).croP(n, geom.n_sec), " )"
        end do
    end do

    write(0, "(a)"); write(11, "(a)")
end subroutine Section_Set_Sectional_Data

! ---------------------------------------------------------------------------------------

! Generate cross-sectional geometry
subroutine Section_Generate_Section_Geometry(geom)
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Generate cross-sectional geometry
    if(geom.sec.types == "square") then

        ! For square lattice
        call Section_Generate_Square(geom)
    else if(geom.sec.types == "honeycomb") then

        ! For honeycomb lattice
        call Section_Generate_Honeycomb(geom)
    end if

    ! Reset local coordinate system according to section ID(direction)
    call Section_Reset_Local_Coordinate(geom)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "3.3. Generate cross-sectional geometry"
        write(i, "(a)")
        call Space(i, 6)
        write(i, "(a)"), "3.4. Update local vectors with the sectional ID"
        write(i, "(a)")
    end do
end subroutine Section_Generate_Section_Geometry

! ---------------------------------------------------------------------------------------

! Generate cross-sectional geometry defined by square lattice
subroutine Section_Generate_Square(geom)
    type(GeomType), intent(inout) :: geom

    double precision :: r, s, t, t_mod, pos(3)
    integer :: i, j, node, pos_row, pos_col

    ! Find cross-sectional position vector
    do i = 1, geom.n_sec
        do j = 1, geom.n_iniL

            ! Set row and column position
            pos_row = geom.sec.posR(i) - geom.sec.minR + 1
            pos_col = geom.sec.posC(i) - geom.sec.minC + 1

            ! Interpolation using the ratio of the natural coordinate system
            !
            !  1-----pos_col-----geom.sec.n_col,  1-----pos_row-----geom.sec.n_row
            ! -1-----   s   -----    1           -1-----   t   -----      1
            !
            !      1-(-1)             (s-(-1))
            ! ------------------ = --------------   --> s = 2*(pos_col-1)/(geom.sec.n_col-1)-1
            ! (geom.sec.n_col-1)    (pos_col-1)
            !
            !      1-(-1)             (t-(-1))
            ! ------------------ = --------------   --> t = 2*(pos_row-1)/(geom.sec.n_row-1)-1
            ! (geom.sec.n_row-1)    (pos_row-1)
            !
            r = 1.0d0
            s = 2.0d0 * dble(pos_col-1) / dble(geom.sec.n_col-1) - 1.0d0
            t = 2.0d0 * dble(pos_row-1) / dble(geom.sec.n_row-1) - 1.0d0

            ! To translate natural coordinate origin to reference row
            t_mod = 2.0d0 * dble(geom.sec.ref_row-1) / dble(geom.sec.n_row-1) - 1.0d0
            t = t - t_mod

            ! Exception handing when there is only one duplex in column or row
            if(geom.sec.n_row == 1) t = 0.0d0

            ! Set cross-sectional position vector at starting part
            pos(1:3) = Section_Get_Position(geom, -r, s, t, j, geom.sec.n_row, geom.sec.n_col)
            node     = geom.iniL(j).poi(1) + (i-1)*geom.n_modP
            geom.croP(node).pos(1:3) = pos(1:3)

            ! Set cross-sectional position vector at endding part
            pos(1:3) = Section_Get_Position(geom, r, s, t, j, geom.sec.n_row, geom.sec.n_col)
            node     = geom.iniL(j).poi(2) + (i-1)*geom.n_modP
            geom.croP(node).pos(1:3) = pos(1:3)

        end do
    end do
end subroutine Section_Generate_Square

! ---------------------------------------------------------------------------------------

! Generate cross-sectional geometry defined by honeycomb lattice
subroutine Section_Generate_Honeycomb(geom)
    type(GeomType), intent(inout) :: geom

    double precision :: r, s, t, t_mod, pos(3)
    integer :: i, j, node, pos_row, pos_col

    ! Find cross-sectional position vector
    do i = 1, geom.n_sec
        do j = 1, geom.n_iniL

            ! Set row and column position
            pos_row = geom.sec.posR(i) - geom.sec.minR + 1
            pos_col = geom.sec.posC(i) - geom.sec.minC + 1

            ! Interpolation using the ratio of the natural coordinate system
            !
            !  1-----3pos_col-----3geom.sec.n_col+1  |   1-----3pos_col-1-----3geom.sec.n_col+1
            ! -1-----   s    -----        1          |  -1-----    s     -----        1
            !                                        |
            !          2               (s+1)         |          2                (s+1)
            ! ------------------ = --------------    |  ------------------ = --------------
            ! (3geom.sec.n_col)     (3pos_col-1)     |  (3geom.sec.n_col)     (3pos_col-2)
            !                                        |
            ! s = 2*(3pos_col-1)/(3geom.sec.n_col)-1 |  s = 2*(3pos_col-2)/(3geom.sec.n_col)-1
            !
            ! ################################################################################
            !
            !  1-----pos_row-----geom.sec.n_row
            ! -1-----   t   -----      1
            !
            !      1-(-1)             (t-(-1))
            ! ------------------ = --------------  --> t = 2*(pos_row-1)/(geom.sec.n_row-1)-1
            ! (geom.sec.n_row-1)    (pos_row-1)
            !

            r = 1.0d0

            ! 3-line generation at each cross-section
            if( (mod(geom.sec.posC(i), 2) == 0 .and. mod(geom.sec.posR(i), 2) == 0) .or. &
                (mod(geom.sec.posC(i), 2) /= 0 .and. mod(geom.sec.posR(i), 2) /= 0) ) then
                s = 2.0d0 * dble(3*pos_col-1) / dble(3*geom.sec.n_col) - 1.0d0
            else
                s = 2.0d0 * dble(3*pos_col-2) / dble(3*geom.sec.n_col) - 1.0d0
            end if

            t = 2.0d0 * dble(pos_row-1) / dble(geom.sec.n_row-1) - 1.0d0

            ! To translate natural coordinate origin to reference row
            t_mod = 2.0d0 * dble(geom.sec.ref_row-1) / dble(geom.sec.n_row-1) - 1.0d0
            t = t - t_mod

            ! Exception handing when there is only one duplex in column or row
            if(geom.sec.n_row == 1) t = 0.0d0

            ! Set cross-sectional position vector at starting part
            pos(1:3) = Section_Get_Position(geom, -r, s, t, j, geom.sec.n_row, geom.sec.n_col)
            node     = geom.iniL(j).poi(1) + (i-1)*geom.n_modP
            geom.croP(node).pos(1:3) = pos(1:3)

            ! Set cross-sectional position vector at endding part
            pos(1:3) = Section_Get_Position(geom,  r, s, t, j, geom.sec.n_row, geom.sec.n_col)
            node     = geom.iniL(j).poi(2) + (i-1)*geom.n_modP
            geom.croP(node).pos(1:3) = pos(1:3)

        end do
    end do
end subroutine Section_Generate_Honeycomb

! ---------------------------------------------------------------------------------------

! Get position vector using natural coordinate system
function Section_Get_Position(geom, r, s, t, num, n_row, n_col) result(pos)
    type(GeomType),   intent(in) :: geom
    double precision, intent(in) :: r, s, t
    integer,          intent(in) :: num, n_row, n_col

    double precision :: pos(3), yz_bar(2), Vs(3), Vt(3)
    double precision :: hr(2), hst(4)

    ! Get shape function
    call Section_Get_Shape_Function(r, s, t, hr, hst)

    ! Director vector
    call Section_Get_Director(geom, num, Vs, Vt)

    ! Get parameter y and z bar
    call Section_Get_Parameter(geom.sec.types, hst, n_row, n_col, yz_bar)

    ! Set position vector from interpolation function
    pos(1) = hr(1) * geom.modP(geom.iniL(num).poi(1)).pos(1)       + &
             hr(2) * geom.modP(geom.iniL(num).poi(2)).pos(1)       + &
             hr(1) * yz_bar(1) * Vt(1) + hr(2) * yz_bar(1) * Vt(1) + &
             hr(1) * yz_bar(2) * Vs(1) + hr(2) * yz_bar(2) * Vs(1)

    pos(2) = hr(1) * geom.modP(geom.iniL(num).poi(1)).pos(2)       + &
             hr(2) * geom.modP(geom.iniL(num).poi(2)).pos(2)       + &
             hr(1) * yz_bar(1) * Vt(2) + hr(2) * yz_bar(1) * Vt(2) + &
             hr(1) * yz_bar(2) * Vs(2) + hr(2) * yz_bar(2) * Vs(2)

    pos(3) = hr(1) * geom.modP(geom.iniL(num).poi(1)).pos(3)       + &
             hr(2) * geom.modP(geom.iniL(num).poi(2)).pos(3)       + &
             hr(1) * yz_bar(1) * Vt(3) + hr(2) * yz_bar(1) * Vt(3) + &
             hr(1) * yz_bar(2) * Vs(3) + hr(2) * yz_bar(2) * Vs(3)
end function Section_Get_Position

! ---------------------------------------------------------------------------------------

! Get 1D and 2D finite element shape function
subroutine Section_Get_Shape_Function(r, s, t, hr, hst)
    double precision, intent(in)  :: r, s, t
    double precision, intent(out) :: hr(2), hst(4)

    ! 1D shape function
    hr(1) = 0.5d0 * (1.0d0 - r)
    hr(2) = 0.5d0 * (1.0d0 + r)

    ! 2D shape function
    hst(1) = 0.25d0 * (1.0d0 - s) * (1.0d0 + t)
    hst(2) = 0.25d0 * (1.0d0 + s) * (1.0d0 + t)
    hst(3) = 0.25d0 * (1.0d0 - s) * (1.0d0 - t)
    hst(4) = 0.25d0 * (1.0d0 + s) * (1.0d0 - t)
end subroutine Section_Get_Shape_Function

! ---------------------------------------------------------------------------------------

! Get director vectors, local coordinate system
subroutine Section_Get_Director(geom, line, Vs, Vt)
    type(GeomType),   intent(in)  :: geom
    integer,          intent(in)  :: line
    double precision, intent(out) :: Vs(3), Vt(3)

    Vs(1:3) = geom.croL(line).t(2, 1:3)
    Vt(1:3) = geom.croL(line).t(3, 1:3)
end subroutine Section_Get_Director

! ---------------------------------------------------------------------------------------

! Get parameter y and z bar
subroutine Section_Get_Parameter(type_sec, hst, n_row, n_col, yz_bar)
    character(10),    intent(in)  :: type_sec
    double precision, intent(in)  :: hst(4)
    integer,          intent(in)  :: n_row, n_col
    double precision, intent(out) :: yz_bar(2)

    double precision :: rad_cylinder, y_scale, z_scale
    double precision :: y(4), z(4)

    !
    ! ^       ^   base    axis
    ! *------>|   |<------*
    ! * 1.0   |   |       *     Radius of helix : 1 nm
    ! *       |   |<------*     Gap between two helixes : 0.25 nm
    ! *       |   |   1.0 *     Cylinder diameter : 2.25 nm
    !       ->|   !<-  gap = 0.25
    ! axis   base   
    !                 z
    !                 ^
    !                 |
    !   (y1,z1) ************* (y2,z2)
    !           *     |     *
    !           *     |     *
    !           *     |-----*-------> y
    !           *           *
    !           *           *
    !   (y3,z3) ************* (y4,z4)
    !
    ! Radius of cylinder
    rad_cylinder = para_rad_helix + para_gap_helix / 2.0d0

    ! Find scale factor
    if(type_sec == "square") then
        !->--->--->
        ! * | * | * |
        ! 1   3   5         : 2->1r, 4->3r, 6->5r
        y_scale = rad_cylinder * dble(n_col - 1)
        z_scale = rad_cylinder * dble(n_row - 1)

    else if(type_sec == "honeycomb") then

        y_scale = 0.5d0 * rad_cylinder *      3.0d0  * dble(n_col)
        z_scale = 0.5d0 * rad_cylinder * sqrt(3.0d0) * dble(n_row - 1)

    end if

    ! Position each corner
    y(1) = -1.0d0 * y_scale
    y(2) =  1.0d0 * y_scale
    y(3) = -1.0d0 * y_scale
    y(4) =  1.0d0 * y_scale

    z(1) =  1.0d0 * z_scale
    z(2) =  1.0d0 * z_scale
    z(3) = -1.0d0 * z_scale
    z(4) = -1.0d0 * z_scale

    ! Return parameter value
    yz_bar(1) = hst(1)*y(1) + hst(2)*y(2) + hst(3)*y(3) + hst(4)*y(4)
    yz_bar(2) = hst(1)*z(1) + hst(2)*z(2) + hst(3)*z(3) + hst(4)*z(4)
end subroutine Section_Get_Parameter

! ---------------------------------------------------------------------------------------

! Reset local coordinates according to section ID
subroutine Section_Reset_Local_Coordinate(geom)
    type(GeomType), intent(inout) :: geom

    double precision :: pos_1(3), pos_2(3), t(3,3)
    integer :: i, point_1, point_2

    ! Reset local coordinate system according to section ID
    do i = 1, geom.n_croL

        point_1  = geom.croL(i).poi(1)
        point_2  = geom.croL(i).poi(2)
        pos_1(:) = geom.croP(geom.croL(i).poi(1)).pos
        pos_2(:) = geom.croP(geom.croL(i).poi(2)).pos

        ! Reset local coordinate t1, t2, t3
        if(mod(geom.croL(i).sec, 2) == 0) then
            ! + z-direction
            t(1,:) = Normalize(pos_2 - pos_1)
            t(2,:) = geom.croL(i).t(2,:)
            t(3,:) = Cross(t(1,:), t(2,:))
        else
            ! - z-direction
            t(1,:) = Normalize(pos_1 - pos_2)
            t(2,:) = geom.croL(i).t(2,:)
            t(3,:) = Cross(t(1,:), t(2,:))
        end if

        ! Update local coordinate system
        geom.croL(i).t(1,:) = t(1,:)
        geom.croL(i).t(2,:) = t(2,:)
        geom.croL(i).t(3,:) = t(3,:)
    end do
end subroutine Section_Reset_Local_Coordinate

! ---------------------------------------------------------------------------------------

! Write cross-sectional geometry
subroutine Section_Chimera_Cro_Geometry(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: pos_1(3), pos_2(3), pos_c(3), pos_a(3), pos_b(3)
    double precision :: length, t1(3), t2(3), t3(3)
    integer :: i, j, iter, nbp
    logical :: f_axis, f_info
    character(200) :: path

    if(para_write_401 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_401_info

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=401, file=trim(path)//"_cro_geo.bild", form="formatted")

    ! Write cross-sectional points
    write(401, "(a)"), ".color red"
    do i = 1, geom.n_croP

        write(401, "(a$    )"), ".sphere "
        write(401, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(401, "(1f9.3 )"), 0.35d0
    end do

    ! Write cross-sectional edges
    write(401, "(a)"), ".color dark green"
    do i = 1, geom.n_croL

        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)

        write(401, "(a$   )"), ".cylinder "
        write(401, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.15d0
    end do

    ! Write information
    if(f_info == .true.) then

        ! Write modified points
        write(401, "(a)"), ".color dark green"
        do i = 1, geom.n_modP
            write(401, "(a$    )"), ".sphere "
            write(401, "(3f9.3$)"), geom.modP(i).pos(1:3)
            write(401, "(1f9.3 )"), 0.3d0
        end do

        ! Write modified edges
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.modP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.modP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = Normalize(pos_2(1:3) - pos_1(1:3))

            length = Norm(pos_2(1:3) - pos_1(1:3))
            iter   = length / 2

            do j = 1, iter * 2
                if(mod(j, 2) == 1) then
                    pos_a(1:3) = pos_1(1:3) + dble(j-1) * pos_c(1:3)
                    pos_b(1:3) = pos_1(1:3) + dble(j)   * pos_c(1:3)

                    write(401, "(a$    )"), ".cylinder "
                    write(401, "(3f9.3$)"), pos_a(1:3)
                    write(401, "(3f9.3$)"), pos_b(1:3)
                    write(401, "(1f9.3 )"), 0.1d0
                end if
            end do
        end do

        ! Information on edge number
        do i = 1, geom.n_croL
            pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

            write(401, "(a$    )"), ".cmov "
            write(401, "(3f9.3 )"), pos_c(1:3) + 0.4d0
            write(401, "(a     )"), ".color dark green"
            write(401, "(a     )"), ".font Helvetica 12 bold"
            write(401, "(i7, a$)"), i,                "("
            write(401, "(i3, a )"), geom.croL(i).sec, ")"
        end do

        ! Local coordinate system on cross-sectional edges
        do i = 1, geom.n_croL
            pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

            t1(1:3) = geom.croL(i).t(1, 1:3)
            t2(1:3) = geom.croL(i).t(2, 1:3)
            t3(1:3) = geom.croL(i).t(3, 1:3)

            write(401, "(a     )"), ".color red"    ! x-axis
            write(401, "(a$    )"), ".arrow "
            write(401, "(3f8.2$)"), pos_c(1:3)
            write(401, "(3f8.2$)"), pos_c(1:3) + t1(1:3) * 1.5d0
            write(401, "(3f8.2 )"), 0.18d0, 0.36d0, 0.6d0

            write(401, "(a     )"), ".color blue"   ! y-axis
            write(401, "(a$    )"), ".arrow "
            write(401, "(3f8.2$)"), pos_c(1:3)
            write(401, "(3f8.2$)"), pos_c(1:3) + t2(1:3) * 1.2d0
            write(401, "(3f8.2 )"), 0.15d0, 0.3d0, 0.6d0

            write(401, "(a     )"), ".color yellow" ! z-axis
            write(401, "(a$    )"), ".arrow "
            write(401, "(3f8.2$)"), pos_c(1:3)
            write(401, "(3f8.2$)"), pos_c(1:3) + t3(1:3) * 1.2d0
            write(401, "(3f8.2 )"), 0.15d0, 0.3d0, 0.6d0
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(401, "(a)"), ".translate 0.0 0.0 0.0"
        write(401, "(a)"), ".scale 0.5"
        write(401, "(a)"), ".color grey"
        write(401, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(401, "(a)"), ".color red"             ! x-axis
        write(401, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(401, "(a)"), ".color blue"            ! y-axis
        write(401, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(401, "(a)"), ".color yellow"          ! z-axis
        write(401, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=401)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=401, file=trim(path)//"_cross_geo.dat", form="formatted")

    write(401, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(401, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(401, "(a$)"), 'ZONE F = FEPOINT'
    write(401, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_croP)))
    write(401, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_croL)))
    write(401, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_croP
        write(401, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(401, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_croL
        write(401, "(1i7$)"), geom.croL(i).poi(1)
        write(401, "(1i7 )"), geom.croL(i).poi(2)
    end do

    close(unit=401)
end subroutine Section_Chimera_Cro_Geometry

! ---------------------------------------------------------------------------------------

end module Section