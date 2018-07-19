!
! =============================================================================
!
! Module - Math
! Last Updated : 04/10/2018, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX, which allows scientists to build and solve
! the sequence design of complex DNAnanostructures.
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
module Math

    implicit none

    double precision, parameter :: pi  = 3.141592653589793d0
    double precision, parameter :: eps = 0.0000001d0

    public Find_Intersection
    public Check_Intersection
    public Find_Cloest_Point
    public Is_Same_Vector
    public Rotate_Vector
    public Normalize
    public Norm
    public Cross
    public Inverse_22

    public Find_Minimum
    public Swap
    public Sort
    public Sort2

    public Rad2Deg
    public Deg2Rad
    public Int2Str
    public Dble2Str
    public Dble2Str1
    public Dble2Str2

    interface Reallocate
        module procedure Reallocate_Int_1D
        module procedure Reallocate_Int_2D
    end interface

    public Math_Determinant
    public Math_Triangularization
    public Math_Get_Low_Tri
    public Math_Set_Entity_One

contains

! -----------------------------------------------------------------------------

! Find intersection point in two 3D vectors
! http://mathforum.org/library/drmath/view/63719.html
function Find_Intersection(a1, a2, b1, b2) result(pos)
    double precision, intent(in) :: a1(3)
    double precision, intent(in) :: a2(3)
    double precision, intent(in) :: b1(3)
    double precision, intent(in) :: b2(3)

    double precision :: pos(3), mat(2,2), invmat(2,2), t1, t2
    logical :: ok

    ok = Check_Intersection(a1, a2, b1, b2)
    print *, ok

    mat(1,:) = [a2(1) - a1(1), -(b2(1)-b1(1))]
    mat(2,:) = [a2(2) - a1(2), -(b2(2)-b1(2))]

    invmat = Inverse_22(mat)

    t1 = invmat(1,1) * (b1(1) - a1(1)) + invmat(1,2) * (b1(2) - a1(2))
    t2 = invmat(2,1) * (b1(1) - a1(1)) + invmat(2,2) * (b1(2) - a1(2))

    pos(1) = a1(1) + (a2(1) - a1(1)) * t1
    pos(2) = a1(2) + (a2(2) - a1(2)) * t1
    pos(3) = a1(3) + (a2(3) - a1(3)) * t1
end function Find_Intersection

! -----------------------------------------------------------------------------

! Check intersection between two 3D vectors
! http://stackoverflow.com/questions/2316490/the-algorithm-to-find-the-point-of-intersection-of-two-3d-line-segment
function Check_Intersection(a1, a2, b1, b2) result(ok)
    double precision, intent(in) :: a1(3)
    double precision, intent(in) :: a2(3)
    double precision, intent(in) :: b1(3)
    double precision, intent(in) :: b2(3)

    double precision :: da(3), db(3), dc(3)
    logical :: ok

    da = a2 - a1
    db = b2 - b1
    dc = b1 - a1

    if(dabs(dot_product(dc, Cross(da, db))) > 0.00001d0) then
        ok = .true.
    else
        ok = .false.
    end if
end function Check_Intersection

! -----------------------------------------------------------------------------

! Closet point between lines
! http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
function Find_Cloest_Point(a1, a2, b1, b2, check) result(pos)
    double precision, intent(in) :: a1(3)
    double precision, intent(in) :: a2(3)
    double precision, intent(in) :: b1(3)
    double precision, intent(in) :: b2(3)
    logical, intent(out) :: check

    double precision, parameter :: SMALL_NUM = 0.00000001d0
    double precision :: pos_u(3), pos_v(3), pos(3), u(3), v(3), w(3), dp(3)
    double precision :: a, b, c, d, e, dd, sc, tc, dist

    u = a2(1:3) - a1(1:3)
    v = b2(1:3) - b1(1:3)
    w = a1(1:3) - b1(1:3)
    a = dot_product(u, u)   ! Always >= 0
    b = dot_product(u, v)
    c = dot_product(v, v)   ! Always >= 0
    d = dot_product(u, w)
    e = dot_product(v, w)

    ! Always >= 0
    dd = a*c - b*b

    ! Compute the line parameters of the two closest points
    if(dd < SMALL_NUM) then

        ! Lines are almost parallel
        sc = 0.0d0

        ! Use the largest denominator
        if(b > c) then
            tc = d / b
        else
            tc = e / c
        end if
        check = .false.
    else

        sc = (b*e - c*d) / dd
        tc = (a*e - b*d) / dd
        check = .true.
    end if

    ! Get the difference of the two closest points
    dP = w + (sc * u) - (tc * v)

    ! The closest distance
    dist = Norm(dP)

    ! The position at the closet point
    pos_u = a1 + u * sc
    pos_v = b1 + v * tc
    pos   = 0.5d0 * (pos_u + pos_v)
end function Find_Cloest_Point

! -----------------------------------------------------------------------------

! Compare two vectors
function Is_Same_Vector(pos_1, pos_2) result(flag)
    double precision, intent(in) :: pos_1(3)
    double precision, intent(in) :: pos_2(3)

    logical :: flag

    flag = .false.

    if( dabs(dabs(pos_1(1)) - dabs(pos_2(1))) < 0.00001d0 .and. &
        dabs(dabs(pos_1(2)) - dabs(pos_2(2))) < 0.00001d0 .and. &
        dabs(dabs(pos_1(3)) - dabs(pos_2(3))) < 0.00001d0 .and. &
        pos_1(1) * pos_2(1) >= 0.0 .and. &
        pos_1(2) * pos_2(2) >= 0.0 .and. &
        pos_1(3) * pos_2(3) >= 0.0 ) then

    flag = .true.
    end if
end function Is_Same_Vector

! -----------------------------------------------------------------------------

! Rotate V2 orientation vector using 3D finite rotation
subroutine Rotate_Vector(vec, pseudo, angle)
    double precision, intent(inout) :: vec(3)
    double precision, intent(in)    :: pseudo(3), angle

    double precision :: theta_x, theta_y, theta_z, gamma
    double precision :: vec_new(3), R(3, 3), S(3, 3)

    ! 3D finite rotate formulation
    theta_x = pseudo(1) * angle
    theta_y = pseudo(2) * angle
    theta_z = pseudo(3) * angle

    gamma = dsqrt(theta_x**2.0d0 + theta_y**2.0d0 + theta_z**2.0d0)

    R(:,:) = 0.0d0
    R(1,1) = 1.0d0
    R(2,2) = 1.0d0
    R(3,3) = 1.0d0

    if(dabs(gamma) /= 0.0d0) then
        S(1,2) = -theta_z
        S(1,3) =  theta_y
        S(2,1) =  theta_z
        S(2,3) = -theta_x
        S(3,1) = -theta_y
        S(3,2) =  theta_x
        R      =  R + dsin(gamma)/gamma * S + 0.5d0 * &
                  (dsin(gamma/2.0d0)/(gamma/2.0d0))**2.0d0 * matmul(S,S)
    end if

    vec_new(1) = R(1,1)*vec(1) + R(1,2)*vec(2) + R(1,3)*vec(3)
    vec_new(2) = R(2,1)*vec(1) + R(2,2)*vec(2) + R(2,3)*vec(3)
    vec_new(3) = R(3,1)*vec(1) + R(3,2)*vec(2) + R(3,3)*vec(3)

    vec_new = Normalize(vec_new)

    vec(1:3) = vec_new(1:3)
end subroutine Rotate_Vector

! -----------------------------------------------------------------------------

! Size vector
function Norm(vec) result(size)
    double precision, intent(in) :: vec(3)

    double precision :: size

    size = dsqrt(vec(1)**2.0d0 + vec(2)**2.0d0 + vec(3)**2.0d0)
end function Norm

! -----------------------------------------------------------------------------

! Normalize vector
function Normalize(vec) result(vector)
    double precision, intent(in) :: vec(3)

    double precision :: length, vector(3)

    length = dsqrt(vec(1)**2.0d0 + vec(2)**2.0d0 + vec(3)**2.0d0)
    vector = vec(1:3) / length
end function Normalize

! -----------------------------------------------------------------------------

! Cross product
function Cross(vec1, vec2) result(vec_cross)
    double precision, intent(in) :: vec1(3), vec2(3)

    double precision :: vec_cross(3)

    vec_cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    vec_cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    vec_cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end function Cross

! -----------------------------------------------------------------------------

! Inverse matrix
function Inverse_22(mat) result (inv_mat)
    double precision, intent(in)  :: mat(2,2)

    double precision :: inv_mat(2,2), cofactor(2,2), det

    det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

    cofactor(1,1) = +mat(2,2)
    cofactor(1,2) = -mat(2,1)
    cofactor(2,1) = -mat(1,2)
    cofactor(2,2) = +mat(1,1)

    inv_mat = transpose(cofactor) / det
end function Inverse_22

! -----------------------------------------------------------------------------

! Return the location of the minimum in the section between start and end
function Find_Minimum(array, start, end) result(value)
    integer, intent(in) :: array(:)
    integer, intent(in) :: start
    integer, intent(in) :: end

    integer :: i, min, loc, value

    min = array(start)  ! Assume the first is the min
    loc = start         ! Record its position

    ! Start with next elements
    do i = start+1, end

        ! If array(i) less than the min?
        if(array(i) < min) then

            ! Record its position
            min = array(i)
            loc = i
        end if
    end do

    ! Return the position
    value = loc
end function Find_Minimum

! -----------------------------------------------------------------------------

! Swap the values of its two formal arguments
subroutine Swap(entity_a, entity_b)
    integer, intent(inout) :: entity_a
    integer, intent(inout) :: entity_b

    integer :: entity_t

    entity_t = entity_a
    entity_a = entity_b
    entity_b = entity_t
end subroutine Swap

! -----------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
subroutine Sort(array, size)
    integer, intent(inout) :: array(:)
    integer, intent(in)    :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array, i, size)

        ! Swap this and the minimum
        call  Swap(array(i), array(loc))
    end do
end subroutine Sort

! -----------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
subroutine Sort2(array1, array2, size)
    integer, intent(inout) :: array1(:)
    integer, intent(inout) :: array2(:)
    integer, intent(in)    :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array1, i, size)

        ! Swap this and the minimum
        call Swap(array1(i), array1(loc))
        call Swap(array2(i), array2(loc))
    end do
end subroutine Sort2

! -----------------------------------------------------------------------------

! Radian to degree
function Rad2Deg(rad) result(deg)
    double precision, intent(in) :: rad
    
    double precision :: deg
    deg = rad * 180.0d0 / pi
end function Rad2Deg

! -----------------------------------------------------------------------------

! Degree to radian
function Deg2Rad(deg) result(rad)
    double precision, intent(in) :: deg
    
    double precision :: rad
    rad = deg * pi / 180.0d0
end function Deg2Rad

! -----------------------------------------------------------------------------

! Integer to string
function Int2Str(int) result(str)
    integer, intent(in) :: int

    character(len=100) :: str

    write(str, "(i)"), int
end function Int2Str

! -----------------------------------------------------------------------------

! Double precision to string
function Dble2Str(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.3)"), dbl
end function Dble2Str

! -----------------------------------------------------------------------------

! Double precision to string
function Dble2Str1(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.1)"), dbl
end function Dble2Str1

! -----------------------------------------------------------------------------

! Double precision to string
function Dble2Str2(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.2)"), dbl
end function Dble2Str2

! -----------------------------------------------------------------------------

! Reallocate integer 1D array
subroutine Reallocate_Int_1D(array, num_new)
    integer, allocatable, intent(inout) :: array(:)
    integer, intent(in) :: num_new

    integer, allocatable :: temp(:)
    integer :: num_old

    num_old = SIZE(array)

    allocate(temp(num_new))

    temp(1:num_old) = array

    call move_alloc(temp, array)
end subroutine Reallocate_Int_1D

! -----------------------------------------------------------------------------

! Reallocate integer 2D array
subroutine Reallocate_Int_2D(array, num_i_new, num_j_new)
    integer, allocatable, intent(inout) :: array(:,:)
    integer, intent(in) :: num_i_new
    integer, intent(in) :: num_j_new

    integer, allocatable :: temp(:,:)
    integer :: num_i_old, num_j_old

    num_i_old = ubound(array, 1)
    num_j_old = ubound(array, 2)

    allocate(temp(num_i_new, num_j_new))

    temp(1:num_i_old,1:num_j_old) = array

    call move_alloc(temp, array)
end subroutine Reallocate_Int_2D
    
! -----------------------------------------------------------------------------

! The determinant of a real square matrix mat by Gauss method with full pivoting
function Math_Determinant(mat, eps) result(det)
    double precision, intent(in) :: mat(:,:)
    double precision, intent(in) :: eps

    double precision, allocatable :: diag(:,:)
    integer, allocatable :: kp(:)
    integer, allocatable :: lp(:)
    double precision :: det
    integer :: i, count, dim, flag

    dim = ubound(mat, 1)

    ! Allocate local matrix diag and vectors kp, lp
    allocate(diag(dim, dim))
    allocate(Kp(dim))
    allocate(Lp(dim))

    ! Triangularization subroutine
    flag = Math_Triangularization(mat, eps, diag, kp, lp)

    if(flag == 0) then
        ! If the matrix singular, det = 0
        det = 0.0d0
    else
        ! If the matrix regular
        det = 1.0d0
        do i = 1, dim
            det = det * diag(i,i)
        end do

        count = 0
        do i = 1, dim - 1
            if(lp(i) /= i) count = count + 1
            if(kp(i) /= i) count = count + 1
        end do

        ! If count is odd
        if(mod(count, 2) /= 0) det = -det
    end if

    ! Dellocate memory
    deallocate(diag)
    deallocate(Kp)
    deallocate(Lp)
end function Math_Determinant

! -----------------------------------------------------------------------------

! The upper triangularization algorithm of Gauss method with full pivoting
function Math_Triangularization(mat, eps, diag, kp, lp) result(flag)
    double precision, intent(in)    :: mat(:,:)
    double precision, intent(in)    :: eps
    double precision, intent(inout) :: diag(:,:)
    integer,          intent(inout) :: kp(:)
    integer,          intent(inout) :: lp(:)

    double precision :: po, t0
    integer :: i, j, k, lo, ko, flag, n

    n = ubound(mat, 1)

    diag = mat
    flag = 1
    k    = 1

    do while(flag == 1 .and. k < n)
        po = diag(k,k)
        lo = k
        ko = k

        do i = k, n
            do j = k, n
                if(dabs(diag(i,j)) > dabs(po)) then
                    po = diag(i,j)
                    lo = i
                    ko = j
                end if
            end do
        end do
        Lp(k) = lo
        Kp(k) = ko

        if(dabs(po) < eps) then
            flag = 0
        else
            if(lo /= k) then
                do j = k, n
                    t0      = diag(k,j)
                    diag(k,j)  = diag(lo,j)
                    diag(lo,j) = t0
                end do
            end if
            if(ko /= k) then
                do i = 1, n
                    t0      = diag(i,k)
                    diag(i,k)  = diag(i,ko)
                    diag(i,ko) = t0
                end do
            end if
            do i = k + 1, n
                diag(i,k) = diag(i,k) / po
                do j = k + 1, n
                    diag(i,j) = diag(i,j) - diag(i,k) * diag(k,j)
                end do
            end do
            k = k + 1
        end if
    end do

    if(flag == 1 .and. dabs(diag(n,n)) < eps) flag = 0
end function Math_Triangularization

! -----------------------------------------------------------------------------

! Get lower triangular matrix
subroutine Math_Get_Low_Tri(mat)
    integer, intent(inout) :: mat(:,:)

    integer :: i, j

    ! Get lower triangular part of matrix
    do i = 1, size(mat, 1)
        do j = 1, size(mat, 2)
            if(i < j) then
                mat(i, j) = 0
            end if
        end do
    end do
end subroutine Math_Get_Low_Tri

! -----------------------------------------------------------------------------

! Set matrix entity one
subroutine Math_Set_Entity_One(adj)
    integer, allocatable, intent(inout) :: adj(:,:)

    integer :: i, j, isize, jsize

    isize = ubound(adj, 1)
    jsize = ubound(adj, 2)

    do i = 1, isize
        do j = 1, jsize
            if(adj(i,j) /= 0) then
                adj(i,j) = 1
            end if
        end do
    end do
end subroutine Math_Set_Entity_One

! -----------------------------------------------------------------------------

end module Math