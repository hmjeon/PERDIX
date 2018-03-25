!
! ---------------------------------------------------------------------------------------
!
!                                       Module - List
!
!                                                                    Updated : 2017/03/27
!
! Comments: The linked list used
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2018. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module List

    implicit none

    public List_Insert_Conn
    public List_Insert_Junc
    public List_Insert_Base
    public List_Insert_Scaf
    public List_Delete_Conn
    public List_Delete_Junc
    public List_Delete_Base
    public List_Delete_Scaf
    public List_Count_Junc
    public List_Count_Base
    public List_Count_Scaf

    type :: ListConn            ! Linked list data for line connectivity
        integer :: point(2)
        type(ListConn), pointer :: next
    end type ListConn

    type :: ListJunc            ! Linked list data for junctions
        integer :: n_arm
        integer :: cnL(100)
        integer :: poi_c
        type(ListJunc), pointer :: next
    end type ListJunc

    type :: ListBase            ! Linked list data for base
        integer :: id
        type(ListBase), pointer :: next
    end type ListBase

    type :: ListScaf            ! Linked list data for scaffold strand
        integer :: id
        type(ListScaf), pointer :: next
    end type ListScaf

contains

! ---------------------------------------------------------------------------------------

! Insert list for connection
function List_Insert_Conn(head, elem) result(list)
    type(ListConn), pointer, intent(in) :: head
    type(ListConn), pointer, intent(inout) :: elem

    type(ListConn), pointer :: list

    ! Insert in the head position
    elem%next => head

    ! Return list
    list => elem
end function List_Insert_Conn

! ---------------------------------------------------------------------------------------

! Insert list for junction
function List_Insert_Junc(head, elem) result(list)
    type(ListJunc), pointer, intent(in) :: head
    type(ListJunc), pointer, intent(in) :: elem

    type(ListJunc), pointer :: list

    ! Insert in the head position
    elem%next => head

    ! Return list
    list => elem
end function List_Insert_Junc

! ---------------------------------------------------------------------------------------

! Insert list for base
function List_Insert_Base(head, elem) result(list)
    type(ListBase), pointer, intent(in) :: head
    type(ListBase), pointer, intent(in) :: elem

    type(ListBase), pointer :: list

    ! Insert in the head position
    elem%next => head

    ! Return list
    list => elem
end function List_Insert_Base

! ---------------------------------------------------------------------------------------

! Insert list for scaffold strand
function List_Insert_Scaf(head, elem) result(list)
    type(ListScaf), pointer, intent(in) :: head
    type(ListScaf), pointer, intent(in) :: elem

    type(ListScaf), pointer :: list

    ! Insert in the head position
    elem%next => head

    ! Return list
    list => elem
end function List_Insert_Scaf

! ---------------------------------------------------------------------------------------

! Delete list for connection
subroutine List_Delete_Conn(self)
    type(ListConn), pointer :: self
    type(ListConn), pointer :: current
    type(ListConn), pointer :: next

    current => self
    do while (associated(current))

        next => current%next
        if(associated(current)) then
            deallocate(current)
            nullify(self)
        end if
        deallocate(current)
        nullify(current)
        current => next
    end do
end subroutine List_Delete_Conn

! ---------------------------------------------------------------------------------------

! Delete list for junction
subroutine List_Delete_Junc(self)
    type(ListJunc), pointer :: self
    type(ListJunc), pointer :: current
    type(ListJunc), pointer :: next

    current => self
    do while (associated(current))

        next => current%next
        if(associated(current)) then
            deallocate(current)
            nullify(self)
        end if
        deallocate(current)
        nullify(current)
        current => next
    end do
end subroutine List_Delete_Junc

! ---------------------------------------------------------------------------------------

! Delete list for base
subroutine List_Delete_Base(self)
    type(ListBase), pointer :: self
    type(ListBase), pointer :: current
    type(ListBase), pointer :: next

    current => self
    do while (associated(current))
        next => current%next
        if(associated(current)) then
            deallocate(current)
            nullify(self)
        end if
        deallocate(current)
        nullify(current)
        current => next
    end do
end subroutine List_Delete_Base

! ---------------------------------------------------------------------------------------

! Delete list for scaffold
subroutine List_Delete_Scaf(self)
    type(ListScaf), pointer :: self
    type(ListScaf), pointer :: current
    type(ListScaf), pointer :: next

    current => self
    do while (associated(current))

        next => current%next
        if(associated(current)) then
            deallocate(current)
            nullify(self)
        end if
        deallocate(current)
        nullify(current)
        current => next
    end do
end subroutine List_Delete_Scaf

! ---------------------------------------------------------------------------------------

! Count junction list
function List_Count_Junc(head) result(count)
    type(ListJunc), pointer, intent(in) :: head

    type(ListJunc), pointer :: ptr
    integer :: count

    ptr => head

    ! Loop for counting
    count = 0
    do while(associated(ptr))
        ptr => ptr%next
        count = count + 1
    end do
end function List_Count_Junc

! ---------------------------------------------------------------------------------------

! Count base list
function List_Count_Base(head) result(count)
    type(ListBase), pointer, intent(in) :: head

    type(ListBase), pointer :: ptr
    integer :: count

    ptr => head

    ! Loop for counting
    count = 0
    do while(associated(ptr))
        ptr => ptr%next
        count = count + 1
    end do
end function List_Count_Base

! ---------------------------------------------------------------------------------------

! Count scaffold list
function List_Count_Scaf(head) result(count)
    type(ListScaf), pointer, intent(in) :: head

    type(ListScaf), pointer :: ptr
    integer :: count

    ptr => head

    ! Loop for counting
    count = 0
    do while(associated(ptr))
        ptr => ptr%next
        count = count + 1
    end do
end function List_Count_Scaf

! ---------------------------------------------------------------------------------------

end module List