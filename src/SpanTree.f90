!
! =============================================================================
!
! Module - SpanTree
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
module SpanTree

    use Math
    use Mani

    implicit none

    ! Counting spanning trees
    public  SpanTree_Count_Spanning_Trees
    public  SpanTree_Count_Spanning_Trees2
    public  SpanTree_Generate_Spanning_Trees

    ! MST algorithm
    public  SpanTree_Kruskal_Algorithm
    public  SpanTree_Prim_Algorithm_1
    public  SpanTree_Prim_Algorithm_2

    ! Manipulation
    public  SpanTree_List2Adj
    public  SpanTree_Adj2List
    public  SpanTree_Print_All_Trees
    public  SpanTree_Print_Matrix
    public  SpanTree_Print_Vector

    ! --------------------------------------------------
    ! Private subroutine definition
    ! --------------------------------------------------
    private SpanTree_Check_Undirected_Graph
    private SpanTree_Generate_Graph_Data
    private SpanTree_Find_Index
    private SpanTree_Delete
    private SpanTree_Undelete
    private SpanTree_Get_Spanning_Tree
    private SpanTree_Is_Bridge

    private SpanTree_Count_Degree_Node
    private SpanTree_Factor_Matrix
    private SpanTree_Swap
    private SpanTree_Get_Determinant

    private SpanTree_Interface_1_Sort
    private SpanTree_Interface_2_Sort
    private SpanTree_Quick_Sort
    private SpanTree_Shell_Sort
    private SpanTree_Remove

contains

! -----------------------------------------------------------------------------

! The number of spanning trees using Kirchhoff's matrix tree theorem
function SpanTree_Count_Spanning_Trees(adj) result(count)
    integer, allocatable, intent(in) :: adj(:,:)

    double precision, allocatable :: lap(:,:)
    integer, allocatable :: deg(:,:), mat(:,:)

    double precision :: det
    integer :: i, j, dim, entity
    integer(kind=8) :: count

    dim = ubound(adj, 1)
    allocate(mat(dim, dim))
    mat = adj

    ! Check adjacent matrix
    call SpanTree_Check_Undirected_Graph(mat)

    ! Get lower triangular matrix
    call Math_Get_Low_Tri(mat)

    ! Set matrix entity one
    call Math_Set_Entity_One(mat)

    mat = mat + transpose(mat)

    ! Degree and Laplacian matrix
    allocate(deg(dim, dim))
    allocate(lap(dim, dim))

    deg(1:dim, 1:dim) = 0
    lap(1:dim, 1:dim) = 0.0d0

    do i = 1, dim
        entity = 0
        do j = 1, dim
            entity = entity + mat(i,j)
        end do
        deg(i,i) = entity
    end do

    lap   = dble(deg) - dble(mat)
    det   = Math_Determinant(lap(2:dim,2:dim), 0.0000000001d0)
    count = nint(det, kind=8)

    ! Deallocate memory
    deallocate(deg, lap)
end function SpanTree_Count_Spanning_Trees

! -----------------------------------------------------------------------------

! Count the number of the spanning trees of a graph
function SpanTree_Count_Spanning_Trees2(adj) result(count)
    integer, allocatable, intent(in) :: adj(:,:)

    double precision, allocatable :: mat(:,:)
    integer, allocatable :: deg(:), piv(:), inp(:,:)

    double precision :: det
    integer :: i, info, dim
    integer(kind = 8) count

    dim = ubound(adj, 1)

    allocate(inp(dim, dim))
    allocate(mat(dim, dim))
    allocate(deg(dim))
    allocate(piv(dim))
    
    inp = adj

    ! Set matrix entity one
    call Math_Set_Entity_One(inp)

    ! Count the degree of each node
    call SpanTree_Count_Degree_Node(inp, deg)

    mat(1:dim,1:dim) = -dble(inp(1:dim,1:dim))

    do i = 1, dim
        mat(i,i) = mat(i,i) + dble(deg(i))
    end do

    ! Factor the dim-1 order matrix
    info = SpanTree_Factor_Matrix(mat, piv, dim-1)

    if(info /= 0) then
        count = 0
    else
        ! Get the determinant
        det   = SpanTree_Get_Determinant(mat, piv, dim, dim-1)
        count = nint(det, kind=8)
    end if

    ! Deallocate memory
    deallocate(inp)
    deallocate(mat)
    deallocate(deg)
    deallocate(piv)
end function SpanTree_Count_Spanning_Trees2

! -----------------------------------------------------------------------------

! Generate all spanning trees for an undirected graph
subroutine SpanTree_Generate_Spanning_Trees(adj, idx, src, dst)
    integer, allocatable, intent(inout) :: adj(:,:)
    integer, allocatable, intent(inout) :: idx(:,:)
    integer, allocatable, intent(inout) :: src(:)
    integer, allocatable, intent(inout) :: dst(:)

    integer, allocatable :: d(:), t(:), n(:), p(:), a(:), s(:), l(:)
    integer, allocatable :: mate(:), idx_map(:)

    integer :: i, j, k, u, e, f, g, x, vv
    integer :: n_vert, n_edge, level, temp
    logical :: skip_s2_to_s4
    integer(kind=8) :: total_num_st

    ! Generate data and define operations on it
    call SpanTree_Generate_Graph_Data(adj, n_vert, n_edge, d, t, n, p, mate, src, dst, idx_map)

    ! Initialization
    skip_s2_to_s4 = .false.
    i             = 0
    total_num_st  = SpanTree_Count_Spanning_Trees(adj)

    ! Edge index matrix for spanning trees
    allocate(idx(n_vert-1, total_num_st))
    do j = 1, n_vert-1
        idx(j, 1:total_num_st) = 0
    end do

    call SpanTree_Get_Spanning_Tree(n_vert, t, n, a)
    !call SpanTree_Print_Vector(a, "a")

    ! Stack for temporarily removed arcs
    allocate(s(n_vert-2))
    s(1:n_vert-2) = 0

    ! Links to successors
    allocate(l(n_vert+2*n_edge))
    l(1:n_vert+2*n_edge) = 0

    ! Current level (Called 'l' in TAOCP)
    ! Arc node that was replaced in last ST to obtain current ST
    level = 1
    x     = 0

    if(n_vert == 2) then
        vv = 1
        e  = n(vv)
        skip_s2_to_s4 = .true.
    end if

99  continue
    do
        if(skip_s2_to_s4 ==.false.) then

            ! Choose edge that should be part of the ST and set u and v to the adjacent vertices
            e  = a(level + 1)
            u  = t(e)
            vv = t(mate(e))

            ! Ensure a higher (or equal) degree at the tail vertex
            if(d(u) > d(vv)) then
                temp = vv
                vv   = u
                u    = temp
                e    = mate(e)
            end if

            ! Contract edge 'e' (v -> u) and combine vertex u into vertex v
            ! Degree of the combined vertex
            k = d(u) + d(vv)
            f = n(u)
            g = 0

            ! Go through the arc node list of u
            ! If the edge leads to v, remove it (not part of reduced graph)
            do while(t(f) /= 0)
                if(t(f) == vv) then
                    call SpanTree_Delete(f, p, n)
                    call SpanTree_Delete(mate(f), p, n)
                    k    = k - 2
                    l(f) = g
                    g    = f
                else
                    ! Mate edge now leads to v (instead of u)
                    t(mate(f)) = vv
                end if
                f = n(f)
            end do

            ! Remember stack of arcs removed for contraction of 'e'
            l(e)  = g
            d(vv) = k

            ! Insert the reduced node arc list of u at the beginning of the node arc list of v
            n(p(u))  = n(vv)
            p(n(vv)) = p(u)
            p(n(u))  = vv
            n(vv)    = n(u)

            ! Edge e is part of the spanning tree (on this level)
            a(level) = e

            ! Move down one level
            level = level + 1

            ! Until the reduced graph has only 2 vertices:
            if(level < n_vert - 1) then

                ! Clear 'removed arc'-stack of new level
                s(level) = 0
                goto 99
            end if
            e = n(vv)
        end if

        ! Reduced graph comprises two vertices now, visit all its spanning trees
        do
            ! Update edge on highest level of the spanning tree
            a(n_vert - 1) = e

            ! New spanning tree in 'a' - save the result
            ! Increase spanning tree counter
            i = i + 1
            do j = 1, n_vert-1
                ! Convert arc node indices to edge indices and save result
                idx(j, i) = idx_map(a(j))
            end do

            ! Save the previously processed edge
            x = e
            e = n(e)
            if(t(e) == 0) then
                exit
            end if
        end do

        ! Loop
        do

            ! Move up one level
            level = level - 1
            if(level == 0) then
                return
            end if
            e  = a(level)
            u  = t(e)
            vv = t(mate(e))

            ! Revert contraction of edge 'e' (vv -> u)
            ! Remove the part of the arc node list of v that was inherited from u during edge contraction
            ! and let both ends of the extracted part point to the header node of u
            n(vv)    = n(p(u))
            p(n(vv)) = vv
            n(p(u))  = u
            p(n(u))  = u
            f        = p(u)

            ! Go through the arc node list of u
            do while(t(f) /= 0)

                ! Mate edge leads to u again (instead of v)
                t(mate(f)) = u
                f = p(f)
            end do
            f = l(e)
            k = d(vv)

            ! Go through stack of arcs removed during contraction of 'e'
            do while(f /= 0)
                k = k + 2

                ! Restore the removed edges
                call SpanTree_Undelete(mate(f), p, n)
                call SpanTree_Undelete(f, p, n)
                f = l(f)
            end do
            d(vv) = k - d(u)

            ! If 'e' is not a bridge, remove it from the graph and repeat S2
            if(SpanTree_Is_Bridge(mate(e), u, vv, n_vert, t, n) == .false.) then
                x = e

                ! Put 'e' on top of the 'removed arc'-stack
                l(e)     = s(level)
                s(level) = e

                ! Remove it from the graph
                call SpanTree_Delete(e, p, n)
                call SpanTree_Delete(mate(e), p, n)
                d(u)  = d(u) - 1
                d(vv) = d(vv) - 1
                exit
            end if

            ! Undo edge deletions of level 'level'
            e = s(level)

            ! Go through the 'removed arc'-stack for this level and recover the edges
            do while(e /= 0)
                u     = t(e)
                vv    = t(mate(e))
                d(u)  = d(u) + 1
                d(vv) = d(vv) + 1

                call SpanTree_Undelete(mate(e), p, n)
                call SpanTree_Undelete(e, p, n)

                e = l(e)
            end do
        end do
    end do

    ! Deallocate memory
    deallocate(d, t, n, p, a, s, l)
    deallocate(mate, idx_map)
end subroutine SpanTree_Generate_Spanning_Trees

! -----------------------------------------------------------------------------

! Kruskal's Algorithm to find minimal spanning tree
subroutine SpanTree_Kruskal_Algorithm(tail, head, cost, tree, n_node, n_edge, length, mode)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer, allocatable, intent(inout) :: tree(:)
    integer,              intent(in)    :: n_node
    integer,              intent(in)    :: n_edge
    integer,              intent(inout) :: length
    character(*),         intent(in)    :: mode

    integer, allocatable :: rot(:), list(:), inic(:), fim(:)

    integer :: es, k, va, an, su, po, ro, no

    ! Allocate memory
    allocate(tree(n_node-1))
    allocate(rot(n_node))
    allocate(list(n_edge))
    allocate(inic(n_node))
    allocate(fim(n_node))

    ! Initialize data
    rot (:) = 0
    list(:) = 0
    inic(:) = 0
    fim (:) = 0
    
    ! Data is sorted by arc costs
    if(mode == "quick") then
        call SpanTree_Quick_Sort(tail, head, cost, n_edge)
    else if(mode == "shell") then
        call SpanTree_Shell_Sort(tail, head, cost, n_edge)
    else
        continue
    end if

    length = 0
    es     = 0
    k      = 1
    va     = 0

    do while(va < n_node - 1)

        an = tail(k)
        su = head(k)
        if(rot(an) + rot(su) == 0) then

            es       = es + 1
            fim(es)  = k
            inic(es) = k
            rot(an)  = es
            rot(su)  = es
            va       = va + 1
            tree(va) = k
            length   = length + cost(k)
        else if(rot(an) == 0) then

            po       = rot(su)
            rot(an)  = po
            list(k)  = inic(po)
            inic(po) = k
            va       = va + 1
            tree(va) = k
            length   = length + cost(k)
        else if(rot(su) == 0) then

            po       = rot(an)
            rot(su)  = po
            list(k)  = inic(po)
            inic(po) = k
            va       = va + 1
            tree(va) = k
            length   = length + cost(k)
        else if(rot(su) /= rot(an)) then

            if(rot(su) < rot(an)) then
                po = rot(su)
                ro = rot(an)
            else
                po = rot(an)
                ro = rot(su)
            end if

            no = inic(po)

10          rot(head(no)) = ro
            rot(tail(no)) = ro

            if(no /= fim(po)) then
                no = list(no)
                goto 10
            end if

            list(fim(rot(an))) = inic(po)
            fim(rot(an))       = fim(po)
            va                 = va + 1
            tree(va)           = k
            length             = length + cost(k)
        end if

        k = k + 1
    end do

    ! Deallocate memory
    deallocate(rot)
    deallocate(list)
    deallocate(inic)
    deallocate(fim)
end subroutine SpanTree_Kruskal_Algorithm

! -----------------------------------------------------------------------------

! Prim's Algorithm, version 1 with quick or shell sort, or without sort
subroutine SpanTree_Prim_Algorithm_1(tail, head, cost, tree, n_node, n_edge, length, mode)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer, allocatable, intent(inout) :: tree(:)
    integer,              intent(in)    :: n_node
    integer,              intent(in)    :: n_edge
    integer,              intent(inout) :: length
    character(*),         intent(in)    :: mode

    integer, allocatable :: up(:), down(:), lsup(:), lsdn(:), node(:), posi(:)
    integer, parameter   :: maxint = 2**31-1

    integer :: h, i, j, k, x, es, fi, me, wo, am, nm, su

    ! Allocate memory
    allocate(tree(n_node-1))
    allocate(up(n_node))
    allocate(down(n_node))
    allocate(lsup(-n_edge:n_edge))
    allocate(lsdn(-n_edge:n_edge))
    allocate(node(n_node))
    allocate(posi(n_node))

    ! Initialize data
    up(:)   = 0
    down(:) = 0
    lsup(:) = 0
    lsdn(:) = 0
    node(:) = 0
    posi(:) = 0

    ! Sorting's interface
    call SpanTree_Interface_1_Sort(tail, head, cost, up, down, lsup, lsdn, mode)

    length   = 0
    fi       = 1
    node(fi) = 1
    posi(1)  = 1
    es       = 0

    do while(es < n_node - 1)
        me = maxint
        do j = 1, fi
            i  = node(j)
            x  = up(i)
            wo = cost(iabs(x))

            if(wo < me) then
                me = wo
                am = x
                h  = iabs(am)
                if(head(h) /= i) then
                    nm = head(h)
                else
                    nm = tail(h)
                end if
            end if
        end do

        es       = es + 1
        k        = iabs(am)
        tree(es) = k
        length   = length + cost(k)
        fi       = fi + 1
        node(fi) = nm
        posi(nm) = fi
        x        = up(nm)

120     if(x > 0) then
            su = head(x)
            if(posi(su) > 0) then
                call SpanTree_Remove(nm,  x, up, down, lsup, lsdn, node, posi, fi)
                call SpanTree_Remove(su, -x, up, down, lsup, lsdn, node, posi, fi)
            end if
        else
            su = tail(-x)
            if(posi(su) > 0) then
                call SpanTree_Remove(su, -x, up, down, lsup, lsdn, node, posi, fi)
                call SpanTree_Remove(nm,  x, up, down, lsup, lsdn, node, posi, fi)
            end if
        end if

        x = lsup(x)
        if(x /= 0) goto 120
    end do

    ! Deallocate memory
    deallocate(up)
    deallocate(down)
    deallocate(lsup)
    deallocate(lsdn)
    deallocate(node)
    deallocate(posi)
end subroutine SpanTree_Prim_Algorithm_1

! -----------------------------------------------------------------------------

! Prim's Algorithm Version 2 with quick or shell sort, or without sort
subroutine SpanTree_Prim_Algorithm_2(tail, head, cost, tree, n_node, n_edge, length, mode)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer, allocatable, intent(inout) :: tree(:)
    integer,              intent(in)    :: n_node
    integer,              intent(in)    :: n_edge
    integer,              intent(inout) :: length
    character(*),         intent(in)    :: mode

    integer, allocatable :: up(:), node(:), list(:)
    logical, allocatable :: marc(:)
    integer, parameter   :: maxint = 2**31-1

    integer :: i, j, x, xx, fi, es, me, wo, no, am

    ! Allocate memory
    allocate(tree(n_node-1))
    allocate(up(n_node))
    allocate(node(n_node))
    allocate(marc(n_node))
    allocate(list(-n_edge:n_edge))

    ! Initialize data
    up(:)   = 0
    node(:) = 0
    marc(:) = 0
    list(:) = 0

    ! Sorting's interface
    call SpanTree_Interface_2_Sort(tail, head, cost, up, list, mode)

    length  = 0
    fi      = 1
    node(1) = 1
    marc(1) = .true.
    es      = 0

    do while(es < n_node - 1)

        me = maxint
        j  = 1
55      i  = node(j)
        if(i == 0 .or. fi == 0) then
            print *, "There is no solution. ", es, " Arcs in tree"
            stop
        end if

        x = up(i)

        if(x == 0) then
            if(fi > j) node(j) = node(fi)
            fi = fi - 1
            j   = j - 1
            goto 65
        end if

60      xx = iabs(x)
        if(marc(tail(xx)) .and. marc(head(xx))) then
            x     = list(x)
            up(i) = x

            if(x == 0) then
                if(fi > j) node(j) = node(fi)
                fi = fi - 1
                j  = j - 1
                goto 65
            else
                goto 60
            end if
        else
            wo = cost(xx)
            if(wo < me) then
                me = wo
                no = i
                am = xx
            end if
        end if
65      j = j + 1
        if(j <= fi) goto 55

        up(no)   = list(up(no))
        es       = es + 1
        tree(es) = am
        length   = length + cost(am)

        fi = fi + 1
        if(head(am) == no) then
            node(fi)       = tail(am)
            marc(tail(am)) = .true.
        else
            if(tail(am) /= no) then
                print *, " Ame t h no ", am, tail(am), head(am), no
                print *, "GAITA"
            end if
            node(fi)       = head(am)
            marc(head(am)) = .true.
        end if
    end do

    ! Deallocate memory
    deallocate(up)
    deallocate(node)
    deallocate(marc)
    deallocate(list)
end subroutine SpanTree_Prim_Algorithm_2

! -----------------------------------------------------------------------------

! Convert from list to adjacent matrix
subroutine SpanTree_List2Adj(adj, tail, head, cost, n_node, n_edge)
    integer, allocatable, intent(inout) :: adj(:,:)
    integer, allocatable, intent(in)    :: tail(:)
    integer, allocatable, intent(in)    :: head(:)
    integer, allocatable, intent(in)    :: cost(:)
    integer,              intent(in)    :: n_node
    integer,              intent(in)    :: n_edge

    integer :: i

    ! If the memory has been allocated
    if(allocated(adj) == .true.) then
        deallocate(adj)
    end if

    ! Allocate and initialize memory
    allocate(adj(n_node, n_node))
    adj(1:n_node,1:n_node) = 0

    ! Set adjacent matrix
    do i = 1, n_edge
        adj(tail(i), head(i)) = cost(i)
        adj(head(i), tail(i)) = cost(i)
    end do
end subroutine SpanTree_List2Adj

! -----------------------------------------------------------------------------

! Convert from adjacent matrix to list
subroutine SpanTree_Adj2List(adj, tail, head, cost, n_node, n_edge)
    integer, allocatable, intent(in)    :: adj(:,:)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer,              intent(in)    :: n_node
    integer,              intent(in)    :: n_edge

    integer :: i, j, count

    ! If the memory has been allocated
    if(allocated(tail) == .true.) then
        deallocate(tail)
    end if

    if(allocated(head) == .true.) then
        deallocate(head)
    end if

    if(allocated(cost) == .true.) then
        deallocate(cost)
    end if

    ! Allocate and initialize memory
    allocate(tail(n_edge)); tail(1:n_edge) = 0
    allocate(head(n_edge)); head(1:n_edge) = 0
    allocate(cost(n_edge)); cost(1:n_edge) = 0

    ! Only upper triangular entity
    count = 0
    do i = 1, n_node
        do j = i + 1, n_node
            if(adj(i,j) /= 0) then
                count       = count + 1
                cost(count) = adj(i,j)
                tail(count) = i
                head(count) = j
            end if
        end do
    end do
end subroutine SpanTree_Adj2List

! -----------------------------------------------------------------------------

! Print all spanning tree
subroutine SpanTree_Print_All_Trees(idx, src, dst)
    integer, allocatable, intent(in) :: idx(:,:)
    integer, allocatable, intent(in) :: src(:)
    integer, allocatable, intent(in) :: dst(:)

    integer :: i, j, n_tree, n_bran

    n_tree = ubound(idx, 2)
    n_bran = ubound(idx, 1)

    call Space(0, 15)
    write(0, "(a)"), "- The number of spanning trees : "&
        //trim(adjustl(Int2Str(n_tree)))
    call Space(0, 15)
    write(0, "(a)"), "- The number of branches       : "&
        //trim(adjustl(Int2Str(n_bran)))

    do i = 1, n_tree
        call Space(0, 18)
        write(0, "(i3, a$)"), i, " spanning tree : "

        do j = 1, n_bran
            write(0, "(a, i2, a, i2, a$)"), "(", &
                src(idx(j,i)), ",", dst(idx(j,i)), ")"

            if(j /= n_bran) then
                write(0, "(a$)"), "  "
            end if
        end do
        write(0, "(a)")
    end do
end subroutine SpanTree_Print_All_Trees

! -----------------------------------------------------------------------------

! Print matrix
subroutine SpanTree_Print_Matrix(adj, str)
    integer, allocatable, intent(in) :: adj(:,:)
    character(*), intent(in) :: str

    integer :: i, j, num_i, num_j

    num_i = ubound(adj, 1)
    num_j = ubound(adj, 2)

    write(0, "(a, 2i5)"), "   Matrix - "//trim(str)//", size : ", num_i, num_j
    do i = 1, num_i
        do j = 1, num_j
            write(0, "(i4$)"), adj(i,j)
        end do
        write(0, "(a)")
    end do
    write(0, "(a)")
end subroutine SpanTree_Print_Matrix

! -----------------------------------------------------------------------------

! Print vector
subroutine SpanTree_Print_Vector(adj, str)
    integer, allocatable, intent(in) :: adj(:)
    character(*), intent(in) :: str

    integer :: i, j, num_i

    num_i = size(adj)

    write(0, "(a, i5)"), "   Vector - "//trim(str)//", size : ", num_i
    do i = 1, num_i
        write(0, "(i4$)"), adj(i)
    end do
    write(0, "(a)")
    write(0, "(a)")
end subroutine SpanTree_Print_Vector

! -----------------------------------------------------------------------------

! Check if the adjacency matrix
subroutine SpanTree_Check_Undirected_Graph(adj)
    integer, allocatable, intent(in) :: adj(:,:)

    integer :: i, j, isize, jsize

    ! 1. Adjacency matrix must be square
    isize = ubound(adj, 1)
    jsize = ubound(adj, 2)
    if(isize /= jsize) then
        write(0, "(a$)"), "Error - Adjacency matrix must be square : "
        write(0, "(a )"), "SpanTree_Check_Undirected_Graph"
        stop
    end if

    ! 2. Graph must not contain self-loops
    do i = 1, isize
        if(adj(i,i) /= 0) then
            write(0, "(a$)"), "Error - Graph must not contain self-loops : "
            write(0, "(a )"), "SpanTree_Check_Undirected_Graph"
            stop
        end if
    end do

    ! 3. Graph must be undirected(Adjacency matrix is symmetric or lower triangular)
    ! lower triangular matrix if a diagonal matrix
    do i = 1, isize
        do j = 1, jsize
            if(i == j) cycle
            if(adj(i,j) /= adj(j,i)) then
                write(0, "(a$)"), "Error - Adjacency matrix is symmetric : "
                write(0, "(a )"), "SpanTree_Check_Undirected_Graph"
                stop
            end if
        end do
    end do

    ! 4. Graph must have two or more vertices
    if(isize <= 1 .or. jsize <= 1) then
        write(0, "(a$)"), "Error - Graph must have two or more vertices : "
        write(0, "(a )"), "SpanTree_Check_Undirected_Graph"
        stop
    end if
end subroutine SpanTree_Check_Undirected_Graph

! -----------------------------------------------------------------------------

! This function generates the data for the undirected simple graph
subroutine SpanTree_Generate_Graph_Data(adj, n_vert, n_edge, d, t, n, p, mate, src, dst, idx_map)
    integer, allocatable, intent(in)    :: adj(:,:)
    integer,              intent(inout) :: n_vert, n_edge
    integer, allocatable, intent(inout) :: d(:), t(:), n(:), p(:)
    integer, allocatable, intent(inout) :: mate(:), src(:), dst(:), idx_map(:)

    integer, allocatable :: mat(:,:)
    integer, allocatable :: d1(:), d2(:), a(:)

    integer :: i, j, ch, dim, sum, count

    allocate(mat(size(adj, 1), size(adj, 2)))
    mat = adj
    call SpanTree_Check_Undirected_Graph(mat)
    
    ! Get lower triangular matrix
    call Math_Get_Low_Tri(mat)

    ! Set matrix entity one
    call Math_Set_Entity_One(mat)

    ! Find index if there is entity
    call SpanTree_Find_Index(mat, dst, src)

    n_vert = size(mat, 1)  ! Number of vertices (in TAOCP: n)
    n_edge = size(src)      ! Number of edges
    !write(0, "(2(a, i5))"), " # of vertices : ", n_vert, ", # of edges : ", n_edge
    !write(0, "(a)")

    ! Degree vector
    dim = ubound(mat, 1)
    allocate(d(dim))
    allocate(d1(dim))
    allocate(d2(dim))

    do i = 1, dim
        sum = 0
        do j = 1, dim
            sum = sum + mat(i,j)
        end do
        d1(i) = sum

        sum = 0
        do j = 1, dim
            sum = sum + mat(j,i)
        end do
        d2(i) = sum
    end do

    d = d1 + d2

    !call SpanTree_Print_Vector(d1, "d1")
    !call SpanTree_Print_Vector(d2, "d2")
    !call SpanTree_Print_Vector(d,  "d")

    ! Tip vector
    allocate(t(2*n_edge+n_vert))
    t(:) = 0
    do i = 1, n_edge
        t(2*i-1 + n_vert) = src(i)
        t(2*i   + n_vert) = dst(i)
    end do

    !call SpanTree_Print_Vector(t, "Tip vector, t")

    ! Next element in arc list
    allocate(n(2*n_edge+n_vert))
    n(1:2*n_edge+n_vert) = 0
    
    ! Previous element in arc list
    allocate(p(2*n_edge+n_vert))
    p(1:2*n_edge+n_vert) = 0

    do i = 1, n_vert
        count = 0
        do j = 1, n_edge
            if(dst(j) == i .or. src(j) == i) count = count + 1
        end do

        allocate(a(count))
        count = 0
        do j = 1, n_edge
            if(dst(j) == i) then
                count = count + 1
                a(count) = n_vert + 2*(j - 1) + 1
            end if
        end do
        do j = 1, n_edge
            if(src(j) == i) then
                count = count + 1
                a(count) = n_vert + 2*j
            end if
        end do
        !call SpanTree_Print_Vector(a, "a")

        n(i)    = a(1)
        p(a(1)) = i
        do j = 1, count
            if(j == count) then
                n(a(j)) = i
                p(i)     = a(j)
            else
                n(a(j))   = a(j+1)
                p(a(j+1)) = a(j)
            end if
        end do

        deallocate(a)
    end do
    !call SpanTree_Print_Vector(n, "n")
    !call SpanTree_Print_Vector(p, "p")

    ! The arc node index of its "mate" edge
    allocate(mate(2*n_edge+n_vert))
    mate(1:2*n_edge+n_vert) = 0

    ! Maps from arc node index of an edge to
    do i = 1, n_edge
        mate(2*i-1 + n_vert) = 2*i + n_vert
        mate(2*i   + n_vert) = 2*(i-1)+1 + n_vert
    end do
    !call SpanTree_Print_Vector(mate, "mate")

    ! Maps from arc node index to edge index in src/dst
    allocate(idx_map(2*n_edge+n_vert))
    idx_map(1:2*n_edge+n_vert) = 0
    do i = 1, n_edge
        idx_map(2*i-1 + n_vert) = i
        idx_map(2*i   + n_vert) = i
    end do
    !call SpanTree_Print_Vector(idx_map, "idx_map")

    deallocate(mat)
    deallocate(d1)
    deallocate(d2)
end subroutine SpanTree_Generate_Graph_Data

! -----------------------------------------------------------------------------

! Find index if there is entity
subroutine SpanTree_Find_Index(adj, dst, src)
    integer, allocatable, intent(in)    :: adj(:,:)
    integer, allocatable, intent(inout) :: dst(:)
    integer, allocatable, intent(inout) :: src(:)

    integer :: i, j, size, count

    size = ubound(adj, 1)

    count = 0
    do i = 1, size
        do j = 1, size
            if(adj(i,j) /= 0) count = count + 1
        end do
    end do

    allocate(dst(count))
    allocate(src(count))

    count = 0
    do i = 1, size
        do j = 1, size
            if(adj(j,i) /= 0) then
                count = count + 1
                src(count) = i
                dst(count) = j
            end if
        end do
    end do

    ! Print progress
    !call SpanTree_Print_Matrix(adj, "adj")
    !call SpanTree_Print_Vector(dst, "dst")
    !call SpanTree_Print_Vector(src, "src")
end subroutine SpanTree_Find_Index

! -----------------------------------------------------------------------------

! Delete
subroutine SpanTree_Delete(a, p, n)
    integer,              intent(in)    :: a
    integer, allocatable, intent(inout) :: n(:), p(:)
    
    n(p(a)) = n(a)
    p(n(a)) = p(a)
end subroutine SpanTree_Delete

! -----------------------------------------------------------------------------

! Undelete
subroutine SpanTree_Undelete(a, p, n)
    integer,              intent(in)    :: a
    integer, allocatable, intent(inout) :: n(:), p(:)

    p(n(a)) = a
    n(p(a)) = a
end subroutine SpanTree_Undelete

! -----------------------------------------------------------------------------

! Generates a spanning tree
subroutine SpanTree_Get_Spanning_Tree(n_vert, t, n, a)
    integer,              intent(in)    :: n_vert
    integer, allocatable, intent(in)    :: t(:), n(:)
    integer, allocatable, intent(inout) :: a(:)
    
    integer, allocatable :: b(:)
    integer :: vv, w, k, e, u

    ! Stores the arc pointers for the spanning tree
    allocate(a(n_vert-1))
    allocate(b(n_vert))
    a(1:n_vert-1) = 0
    b(1:n_vert)   = 0

    vv   = 1            ! Start at (header node of) vertex 1
    w    = 1            ! Last visited vertex: Vertex 1
    b(1) = 1            ! Indicate vertex 1 was visited
    k    = n_vert - 1   ! Need to visit n_vert - 1 more vertices

    do

        ! Get first edge in vertex v's node arc list
        e = n(vv)

        ! Continue until end of node arc list
        do while(t(e) /= 0)

            ! Get tip vertex of this edge
            u = t(e)

            ! If vertex u was not visited yet:
            if(b(u) == 0) then

                ! Indicate u was visited & save previously visited vertex
                ! Update the last visited vertex
                ! Add this edge to the spanning tree
                ! One vertex less on our 'to-visit-list'
                b(u) = w
                w    = u
                a(k) = e
                k    = k - 1

                ! All vertices were visited?
                if(k == 0) then
                    goto 100
                end if
            end if

            ! Process the next arc in v's node arc list
            e = n(e)
        end do

        ! We were not able to visit all n_vert vertices by following edges
        if(w == 1) then
            write(0, "(a)"), "Graph must be connected."
            stop
        end if

        ! Move on to (header node of) vertex that was added last
        ! If w doesn't lead to unvisited vertices, continue with previous vertex
        vv = w
        w  = b(w)
    end do
100 continue

    ! Deallocate memory
    deallocate(b)
end subroutine SpanTree_Get_Spanning_Tree

! -----------------------------------------------------------------------------

! The bridge test for edge 'e' (v -> u), for which 'e_mate' is its mate
function SpanTree_Is_Bridge(e_mate, u, vv, n_vert, t, n) result(flag)
    integer,              intent(in) :: e_mate, u, vv, n_vert
    integer, allocatable, intent(in) :: t(:), n(:)

    integer, allocatable :: b(:)

    integer :: w, f, v_prime, uu
    logical :: flag

    uu = u

    allocate(b(n_vert))

    ! Stores the "vertex check list", which contains the vertices
    ! that are reachable from u but have not yet been checked if they lead to v
    b(1:n_vert) = 0
    w           = uu
    b(w)        = vv    ! Indicate the end of the "vertex check list"

    do
        f = n(uu)
        do while(t(f) /= 0)                 ! Go through the arc node list of u

            v_prime = t(f)
            if(b(v_prime) == 0) then        ! If vertex v_prime was not visited yet:

                if(v_prime /= vv) then      ! Add v_prime to the "vertex check list"
                    b(v_prime) = vv
                    b(w) = v_prime
                    w = v_prime
                else if(f /= e_mate) then   ! Edge is not e_mate? We found a path from u to v!
                    flag = .false.          ! => Edge e/e_mate is not a bridge
                    goto 100
                end if
            end if

            f = n(f)
        end do

        uu = b(uu)                          ! Process the next vertex of the "vertex check list"
        if(uu == vv) then                   ! End of "vertex check list"?
            flag = .true.                   ! then v is only reachable from u via via e_mate!
            goto 100
        end if
    end do
100 continue

    ! Deallocate memory
    deallocate(b)
end function SpanTree_Is_Bridge

! -----------------------------------------------------------------------------

! Count the degree of each node
subroutine SpanTree_Count_Degree_Node(adj, deg)
    integer, intent(in)    :: adj(:,:)
    integer, intent(inout) :: deg(:)

    integer :: i, j, n_node

    ! The number of nodes
    n_node = ubound(adj, 1)

    deg(1:n_node) = 0

    do i = 1, n_node
        do j = 1, n_node
            if(adj(i,j) /= 0) then
                deg(i) = deg(i) + adj(i,j)
            end if
        end do
    end do
end subroutine SpanTree_Count_Degree_Node

! -----------------------------------------------------------------------------

! Factor a general matrix
! n : the order of the matrix, piv : a vector of pivot indices
function SpanTree_Factor_Matrix(mat, piv, n) result(info)
    double precision, intent(inout) :: mat(:,:)
    integer,          intent(inout) :: piv(:)
    integer,          intent(in)    :: n

    integer :: i, j, k, l, info

    ! info, singularity flag, 0 - no singular matrix
    info = 0

    do k = 1, n - 1

        ! Find l, the index of the pivot row
        l = k
        do i = k + 1, n
            if(abs(mat(l,k)) < abs(mat(i,k))) then
                l = i
            end if
        end do

        piv(k) = l

        ! If the pivot index is zero, the algorithm has failed
        if(mat(l,k) == 0.0d0) then
            info = k
            write(0, "(a)")
            write(0, "(a, i8)"), "Error - Zero pivot on step ", info
            stop
        end if

        ! Interchange rows l and k if necessary
        if(l /= k) call SpanTree_Swap(mat(l,k), mat(k,k))

        ! Normalize the values that lie below the pivot entry
        mat(k+1:n,k) = -mat(k+1:n,k) / mat(k,k)

        ! Row elimination with column indexing
        do j = k + 1, n
            if(l /= k) call SpanTree_Swap(mat(l,j), mat(k,j))
            mat(k+1:n,j) = mat(k+1:n,j) + mat(k+1:n,k) * mat(k,j)
        end do
    end do

    piv(n) = n

    if(mat(n,n) == 0.0d0) then
        info = n
        write(0, "(a)")
        write(0, "(a,i8)"), "Error : Zero pivot on step ", info
        stop
    end if
end function SpanTree_Factor_Matrix

! -----------------------------------------------------------------------------

! Swaps two double precision values
subroutine SpanTree_Swap(x, y)
    double precision, intent(inout) :: x
    double precision, intent(inout) :: y

    double precision :: temp
    temp = x
    x    = y
    y    = temp
end subroutine SpanTree_Swap

! -----------------------------------------------------------------------------

! Get the determinant of a matrix factored
function SpanTree_Get_Determinant(mat, piv, lda, n) result(det)
    integer (kind = 4) lda
    integer (kind = 4) n
    real    (kind = 8) mat(lda,n)
    integer (kind = 4) piv(n)

    double precision :: det
    integer :: i

    det = 1.0d0

    do i = 1, n
        det = det * mat(i,i)
    end do

    do i = 1, n
        if(piv(i) /= i) then
            det = -det
        end if
    end do
end function SpanTree_Get_Determinant

! -----------------------------------------------------------------------------

! Sorting's interface 1
subroutine SpanTree_Interface_1_Sort(tail, head, cost, up, down, listup, listdn, mode)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer, allocatable, intent(inout) :: up(:)
    integer, allocatable, intent(inout) :: down(:)
    integer, allocatable, intent(inout) :: listup(:)
    integer, allocatable, intent(inout) :: listdn(:)
    character(*),         intent(in)    :: mode

    integer :: i, k, x, kk, n_edge

    n_edge = size(tail)

    ! Data is sorted by arc costs
    if(mode == "quick") then
        call SpanTree_Quick_Sort(tail, head, cost, n_edge)
    else if(mode == "shell") then
        call SpanTree_Shell_Sort(tail, head, cost, n_edge)
    else
        continue
    end if

    do k = 1, n_edge
        i = tail(k)
        if(up(i) /= 0) then
            x         = down(i)
            listup(x) = k
            listdn(k) = x
            down(i)   = k
        else
            up(i)   = k
            down(i) = k
        end if

        kk = -k
        i  = head(k)
        if(up(i) /= 0) then
            x          = down(i)
            listup(x)  = kk
            listdn(kk) = x
            down(i)    = kk
        else
            up(i)   = kk
            down(i) = kk
        end if
    end do
end subroutine SpanTree_Interface_1_Sort

! -----------------------------------------------------------------------------

! Sorting's interface 2
! The data structure used in the code is dynamic; however, arcs that forms a circuit 
! with the (actual) spanning tree are marked to be removed as soon as "they are the cheapest arc"
subroutine SpanTree_Interface_2_Sort(tail, head, cost, up, list, mode)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer, allocatable, intent(inout) :: up(:)
    integer, allocatable, intent(inout) :: list(:)
    character(*)                        :: mode

    integer :: n_edge, i, k, kk

    n_edge = size(tail)

    ! Data is sorted by arc costs
    if(mode == "quick") then
        call SpanTree_Quick_Sort(tail, head, cost, n_edge)
    else if(mode == "shell") then
        call SpanTree_Shell_Sort(tail, head, cost, n_edge)
    else
        continue
    end if

    do 40 k = 1, n_edge
        i        = tail(k)
        list(k)  = up(i)
        up(i)    = k
        i        = head(k)
        kk       = -k
        list(kk) = up(i)
        up(i)    = kk
40  continue
end subroutine SpanTree_Interface_2_Sort

! -----------------------------------------------------------------------------

! Iterative quick sort algorithm
subroutine SpanTree_Quick_Sort(tail, head, cost, n_edge)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer,              intent(in)    :: n_edge

    integer, allocatable :: sl(:), sr(:)

    integer :: i, j, s, l, r, ha, me, temp

    allocate(sl(-n_edge:n_edge))
    allocate(sr(-n_edge:n_edge))

    s     = 1
    sr(1) = n_edge
    sl(1) = 1

    do while(s /= 0)
        l  = sl(s)
        r  = sr(s)
        s  = s - 1
10      i  = l
        j  = r
        ha = (l + r) / 2
        me = cost(ha)

20      if(cost(i) < me) then
            i = i + 1
            goto 20
        end if

       do while(cost(j) > me)
            j = j - 1
        end do

        if(i <= j) then
            temp    = head(i)
            head(i) = head(j)
            head(j) = temp
            temp    = tail(i)
            tail(i) = tail(j)
            tail(j) = temp
            temp    = cost(i)
            cost(i) = cost(j)
            cost(j) = temp

            i = i + 1
            j = j - 1
        end if

        if(i <= j) goto 20

        if(i < r) then
            s     = s + 1
            sl(s) = i
            sr(s) = r
        end if

        r = j
        if(l < r) goto 10
    end do

    ! Deallocate memory
    deallocate(sl)
    deallocate(sr)
end subroutine SpanTree_Quick_Sort

! -----------------------------------------------------------------------------

! Shell sort algorithm
subroutine SpanTree_Shell_Sort(tail, head, cost, n_edge)
    integer, allocatable, intent(inout) :: tail(:)
    integer, allocatable, intent(inout) :: head(:)
    integer, allocatable, intent(inout) :: cost(:)
    integer,              intent(in)    :: n_edge

    integer :: i, k, x, h, temp

    h = n_edge - 1

    do while (h >= 1)
        do i = 1, n_edge - h
            k = i
20          if(k >= 1) then
                x = k + h
                if(cost(x) < cost(k)) then
                    temp    = tail(x)
                    tail(x) = tail(k)
                    tail(k) = temp
                    temp    = head(x)
                    head(x) = head(k)
                    head(k) = temp
                    temp    = cost(x)
                    cost(x) = cost(k)
                    cost(k) = temp
                    k = k - h
                    goto 20
                end if
            end if
        end do
        h = h / 2
    end do
end subroutine SpanTree_Shell_Sort

! -----------------------------------------------------------------------------

subroutine SpanTree_Remove(ant, arc, up, down, listup, listd, node, posicao, fim)
    integer, intent(in) :: ant, arc
    integer, allocatable, intent(inout) :: up(:)
    integer, allocatable, intent(inout) :: down(:)
    integer, allocatable, intent(inout) :: node(:)
    integer, allocatable, intent(inout) :: posicao(:)
    integer, allocatable, intent(inout) :: listup(:)
    integer, allocatable, intent(inout) :: listd(:)
    integer,              intent(inout) :: fim

    integer :: x, y, z

    if(up(ant) /= arc) then
        y = listd(arc)
        if(listup(arc) /= 0) then
            z         = listup(arc)
            listup(y) = z
            listd(z)  = y
        else
            listup(y) = 0
            down(ant) = y
        end if
    else
        if(down(ant) /= arc) then
            x       = listup(arc)
            up(ant) = x
        else
            up(ant)    = 0
            down(ant)  = 0
            x          = posicao(ant)
            y          = node(fim)
            node(x)    = y
            posicao(y) = x
            fim        = fim - 1
        end if
    end if
end subroutine SpanTree_Remove

! -----------------------------------------------------------------------------

end module SpanTree