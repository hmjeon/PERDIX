!
! ---------------------------------------------------------------------------------------
!
!                                   Module for Math
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/04/29
!                                                            Last  modified  : 2016/07/14
!
! ---------------------------------------------------------------------------------------
!
module Math

    implicit none

    public Rotate_Vector
    public Normalize_Vector
    public Size_Vector
    public Cross_Product_General
    public Cross_Product_Exception
    public Cross_Product
    public Cross_Product_Sub
    public Inverse_Matrix
    public Determinant_Matrix

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

    double precision, parameter :: pi      = 3.1415926535898d0  ! pi
    double precision, parameter :: epsilon = 0.00000001d0       ! floating point comparison

contains

! ---------------------------------------------------------------------------------------

! Rotate V2 orientation vector using 3D finite rotation
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
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

    vec_new = Normalize_Vector(vec_new)

    vec(1:3) = vec_new(1:3)
end subroutine Rotate_Vector

! ---------------------------------------------------------------------------------------

! Size vector
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Size_Vector(vec) result(size)
    double precision, intent(in) :: vec(3)

    double precision :: size

    size = dsqrt(vec(1)**2.0d0 + vec(2)**2.0d0 + vec(3)**2.0d0)
end function Size_Vector

! ---------------------------------------------------------------------------------------

! Normalize vector
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Normalize_Vector(vec) result(vector)
    double precision, intent(in) :: vec(3)

    double precision :: length, vector(3)

    length = dsqrt(vec(1)**2.0d0 + vec(2)**2.0d0 + vec(3)**2.0d0)
    vector = vec(1:3) / length
end function Normalize_Vector

! ---------------------------------------------------------------------------------------

! Cross product
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
subroutine Cross_Product_General(vec1, vec2, vec_cross)
    double precision, intent(in)  :: vec1(3), vec2(3)
    double precision, intent(out) :: vec_cross(3)

    vec_cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    vec_cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    vec_cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end subroutine Cross_Product_General

! ---------------------------------------------------------------------------------------

! Cross product
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Cross_Product_Exception(vec1, vec2) result(cross)
    double precision, intent(in) :: vec1(3), vec2(3)

    double precision :: cross(3)

    ! If vec1 is equal to ey
    if( dabs(vec1(1)) == dabs(vec2(1)) .and. dabs(vec1(2)) == dabs(vec2(2)) .and. dabs(vec1(3)) == dabs(vec2(3)) ) then
        cross(1) = 0.0d0
        cross(2) = 0.0d0
        cross(3) = 1.0d0
    else
        cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end if
end function Cross_Product_Exception

! ---------------------------------------------------------------------------------------

! Cross product
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Cross_Product(vec1, vec2) result(cross)
    double precision, intent(in)  :: vec1(3), vec2(3)

    double precision :: cross(3)

    cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end function Cross_Product

! ---------------------------------------------------------------------------------------

! Cross product
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
subroutine Cross_Product_Sub(vec1, vec2, vec3)
    double precision, intent(in)    :: vec1(3), vec2(3)
    double precision, intent(inout) :: vec3(3)

    vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end subroutine Cross_Product_Sub

! ---------------------------------------------------------------------------------------

! Inverse matrix
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Inverse_Matrix(mat) result(inv_mat)
    double precision, intent(in) :: mat(3, 3)

    double precision :: det, inv_mat(3, 3)

    inv_mat(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
    inv_mat(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
    inv_mat(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
    inv_mat(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
    inv_mat(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
    inv_mat(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
    inv_mat(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
    inv_mat(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
    inv_mat(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

    det = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
        - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
        + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

    inv_mat(:,:) = inv_mat(:,:) / det
end function Inverse_Matrix

! ---------------------------------------------------------------------------------------

! Determinant
! Last updated on Wednesday 7 Apr 2016 by Hyungmin
function Determinant_Matrix(mat) result(deter)
    double precision, intent(in)  :: mat(3,3)

    double precision :: deter

    deter &
        = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
        - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
        + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
end function Determinant_Matrix

! ---------------------------------------------------------------------------------------

! Return the location of the minimum in the section between start and end
! Last updated on Tuesday 26 Apr 2016 by Hyungmin
function Find_Minimum(array, start, end) result(value)
    integer, dimension(1:), intent(in) :: array
    integer,                intent(in) :: start
    integer,                intent(in) :: end

    integer :: i, min, loc, value

    min = array(start)  ! assume the first is the min
    loc = start			! record its position

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

! ---------------------------------------------------------------------------------------

! Swap the values of its two formal arguments
! Last updated on Tuesday 26 Apr 2016 by Hyungmin
subroutine Swap(entity_a, entity_b)
    integer, intent(inout) :: entity_a
    integer, intent(inout) :: entity_b

    integer :: entity_t

    entity_t = entity_a
    entity_a = entity_b
    entity_b = entity_t
end subroutine Swap

! ---------------------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
! Last updated on Tuesday 26 Apr 2016 by Hyungmin
subroutine Sort(array, size)
    integer, dimension(1:), intent(inout) :: array
    integer,                intent(in)    :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array, i, size)

        ! Swap this and the minimum
        call  Swap(array(i), array(loc))
    end do
end subroutine Sort

! ---------------------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
! Last updated on Tuesday 26 Apr 2016 by Hyungmin
subroutine Sort2(array1, array2, size)
    integer, dimension(1:), intent(inout) :: array1
    integer, dimension(1:), intent(inout) :: array2
    integer,                intent(in)    :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array1, i, size)

        ! Swap this and the minimum
        call  Swap(array1(i), array1(loc))
        call  Swap(array2(i), array2(loc))
    end do
end subroutine Sort2

! ---------------------------------------------------------------------------------------

! Radian to degree
! Last updated on Friday 29 Apr 2016 by Hyungmin
function Rad2Deg(rad) result(deg)
    double precision, intent(in) :: rad
    
    double precision :: deg
    deg = rad * 180.0d0 / pi
end function Rad2Deg

! ---------------------------------------------------------------------------------------

! Degree to radian
! Last updated on Monday 2 May 2016 by Hyungmin
function Deg2Rad(deg) result(rad)
    double precision, intent(in) :: deg
    
    double precision :: rad
    rad = deg * pi / 180.0d0
end function Deg2Rad

! ---------------------------------------------------------------------------------------

! Integer to string
! Last updated on Wendesday 4 May 2016 by Hyungmin
function Int2Str(int) result(str)
    integer, intent(in) :: int

    character(len=100) :: str

    write(str, "(i)"), int
end function Int2Str

! ---------------------------------------------------------------------------------------

! Double precision to string
! Last updated on Wendesday 4 May 2016 by Hyungmin
function Dble2Str(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.3)"), dbl
end function Dble2Str

! ---------------------------------------------------------------------------------------

! Double precision to string
! Last updated on Thursday 20 October 2016 by Hyungmin
function Dble2Str1(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.1)"), dbl
end function Dble2Str1

! ---------------------------------------------------------------------------------------

! Double precision to string
! Last updated on Thursday 20 October 2016 by Hyungmin
function Dble2Str2(dbl) result(str)
    double precision, intent(in) :: dbl

    character(len=100) :: str

    write(str, "(f0.2)"), dbl
end function Dble2Str2

! ---------------------------------------------------------------------------------------

! Reallocate integer 1D array
! Last updated on Wendesday 4 May 2016 by Hyungmin
subroutine Reallocate_Int_1D(array, num_new)
    integer, allocatable, intent(inout) :: array(:)
    integer, intent(in) :: num_new

    integer, allocatable, dimension(:) :: temp
    integer :: num_old

    num_old = SIZE(array)

    allocate(temp(num_new))

    temp(1:num_old) = array

    call move_alloc(temp, array)
end subroutine Reallocate_Int_1D

! ---------------------------------------------------------------------------------------

! Reallocate integer 2D array
! Last updated on Wendesday 4 May 2016 by Hyungmin
subroutine Reallocate_Int_2D(array, num_i_new, num_j_new)
    integer, allocatable, intent(inout) :: array(:,:)
    integer, intent(in) :: num_i_new
    integer, intent(in) :: num_j_new

    integer, allocatable, dimension(:,:) :: temp
    integer :: num_i_old, num_j_old

    num_i_old = ubound(array, 1)
    num_j_old = ubound(array, 2)

    allocate(temp(num_i_new, num_j_new))

    temp(1:num_i_old,1:num_j_old) = array

    call move_alloc(temp, array)
end subroutine Reallocate_Int_2D
    
! ---------------------------------------------------------------------------------------

! The determinant of a real square matrix mat by Gauss method with full pivoting
! Last updated on Sunday 8 May 2016 by Hyungmin
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

! ---------------------------------------------------------------------------------------

! The upper triangularization algorithm of Gauss method with full pivoting
! Last updated on Sunday 8 May 2016 by Hyungmin
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

! ---------------------------------------------------------------------------------------

! Get lower triangular matrix
! Last updated on Sunday 8 May 2016 by Hyungmin
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

! ---------------------------------------------------------------------------------------

! Set matrix entity one
! Last updated on Thursday 5 May 2016 by Hyungmin
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

! ---------------------------------------------------------------------------------------

end module Math