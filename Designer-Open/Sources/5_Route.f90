!
! ---------------------------------------------------------------------------------------
!
!                                    Module for Route
!
!                                                            First programed : 2015/04/21
!                                                            Last  modified  : 2016/03/21
!
! Comments: The module is to generate the dual graph and compute spanning tree.
!
! by Hyungmin Jun (Hyungminjun@outlook.com), MIT, Bathe Lab, 2017
!
! Copyright 2017. Massachusetts Institute of Technology. Rights Reserved.
! M.I.T. hereby makes following copyrightable material available to the
! public under GNU General Public License, version 2 (GPL-2.0). A copy of
! this license is available at https://opensource.org/licenses/GPL-2.0
!
! ---------------------------------------------------------------------------------------
module Route

    use Ifport

    use Data_Prob
    use Data_Geom
    use Data_Bound
    use Data_Mesh
    use Data_DNA
    use SpanTree

    use Section

    use Para
    use List
    use Mani
    use Math

    implicit none

    public  Route_Generation

    private Route_Init_Base_Connectivity
    private Route_Set_Base_Position
    private Route_Chimera_Route
    private Route_Connect_Strand_Junc
    private Route_Connect_Scaf
    private Route_Connect_Stap
    private Route_Add_Nucleotide
    private Route_Set_Strand_ID_Scaf
    private Route_Find_Centered_Scaf_Xover
    private Route_Check_Nei_Xover
    private Route_Split_Centered_Scaf_Xover
    private Route_Write_Centered_Scaf_Xover
    private Route_Modify_Scaf_Xover

    ! Function for graph theory
    private Route_Graph_Build_Origami
    private Route_Graph_Allocate_Data
    private Route_Graph_Set_Data
    private Route_Graph_Write_Adjacent
    private Route_Graph_Chimera_All_Spanning_Tree
    private Route_Graph_Write_List
    private Route_Graph_Chimera_Spanning_Tree
    private Route_Graph_Delete_Scaf_Xover

    private Route_Make_Scaf_Origami
    private Route_Find_Possible_Stap_Xover
    private Route_Set_Stap_Crossover
    private Route_Set_Orientation
    private Route_Chimera_Crossovers
    private Route_Chimera_Orientation
    private Route_Chimera_Atom

    private Route_Delete_End_Pos_Xover_Scaf
    private Route_Delete_End_Pos_Xover_Stap
    private Route_Center_Scaf_Crossover

contains

! ---------------------------------------------------------------------------------------

! Generate B-form DNA (set dna.base_scaf and dna.base_stap)
subroutine Route_Generation(prob, geom, bound, mesh, dna)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(in)    :: bound
    type(MeshType),  intent(inout) :: mesh
    type(DNAType),   intent(inout) :: dna

    integer :: i

    ! Print progress
    do i = 0, 11, 11
        write(i, "(a)")
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   |              5. Build B-form DNA and scaffold route                |"
        write(i, "(a)"), "   |                                                                    |"
        write(i, "(a)"), "   +--------------------------------------------------------------------+"
        write(i, "(a)")
    end do

    ! Initialize base connectivity
    call Route_Init_Base_Connectivity(mesh, dna)

    ! Set base position vector
    call Route_Set_Base_Position(geom, mesh, dna)

    ! Write scaffold route, step 1 - no connection at the junction
    call Route_Chimera_Route(prob, geom, mesh, dna, "route1")

    ! Connect strands at the junction by unpaired nucleotides
    call Route_Connect_Strand_Junc(geom, bound, mesh, dna)

    ! Set the strand ID in scaffold
    call Route_Set_Strand_ID_Scaf(dna)

    ! Write scaffold route, step 2 - connected junction
    call Route_Chimera_Route(prob, geom, mesh, dna, "route2")

    ! Find centered scaffold crossovers
    call Route_Find_Centered_Scaf_Xover(prob, geom, mesh, dna)

    ! Write centered scaffold crossovers
    call Route_Write_Centered_Scaf_Xover(prob, geom, mesh, dna)

    ! Modify scaffold crossovers to avoid duplicated indication
    call Route_Modify_Scaf_Xover(geom, mesh, dna)

    ! Write scaffold route, step 3 - strand with possible crossovers
    call Route_Chimera_Route(prob, geom, mesh, dna, "route3")

    ! Build scaffold DNA origami
    if(para_method_MST == "prim" .or. para_method_MST == "kruskal") then

        ! Prim or Kruskal algorithm
        call Route_Graph_Build_Origami(prob, mesh, dna)
    else if(para_method_MST == "greedy") then

        ! Greedy algorithm
        call Route_Make_Scaf_Origami(prob, mesh, dna)
    end if

    ! Write scaffold route, step 4 - one crossover at each edge in the section
    call Route_Chimera_Route(prob, geom, mesh, dna, "route4")

    ! Find possible staple crossovers
    call Route_Find_Possible_Stap_Xover(geom, mesh, dna)

    ! Set staple crossover
    !call Route_Set_Stap_Crossover(prob, dna)

    ! Write scaffold route, step 4 - continuous scaffold strand
    call Route_Chimera_Route(prob, geom, mesh, dna, "route5")

    ! Set three orientation vectors based on 3DNA convention
    call Route_Set_Orientation(mesh, dna)

    ! Write crossovers based on base pairs
    call Route_Chimera_Crossovers(prob, geom, bound, mesh, dna)

    ! Write 3 orientation vectors based on 3DNA convention
    call Route_Chimera_Orientation(prob, mesh, dna)

    ! Write atomic model
    call Route_Chimera_Atom(prob, geom, dna)
end subroutine Route_Generation

! ---------------------------------------------------------------------------------------

! Initialize base connectivity from mesh data
! Last updated on Wednesday 12 May 2016 by Hyungmin
subroutine Route_Init_Base_Connectivity(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i

    ! Set the initial number of scaf / stap strands from # of nodes
    dna.n_base_scaf = mesh.n_node
    dna.n_base_stap = mesh.n_node

    ! Allocate scaffold and staple strand data
    allocate(dna.base_scaf(dna.n_base_scaf))
    allocate(dna.base_stap(dna.n_base_stap))

    ! Set ID and across of bases in scaffold and staple strands
    do i = 1, mesh.n_node

        ! The ID of strand comes from the one of the bp
        dna.base_scaf(i).id   = mesh.node(i).id      ! ID for bases in scaffold
        dna.base_stap(i).id   = mesh.node(i).id      ! ID for bases in staple
        dna.base_scaf(i).node = mesh.node(i).id      ! Heritage node id
        dna.base_stap(i).node = mesh.node(i).id      ! Heritage node id

        ! Set connectivity for bases in scaffold, which is the same to nodal connectivity
        dna.base_scaf(i).up = mesh.node(i).up
        dna.base_scaf(i).dn = mesh.node(i).dn

        ! Set connectivity for bases in staple, which is the opposite to the nodal connectivity
        if(dna.base_scaf(i).up == -1) then
            dna.base_stap(i).up = dna.base_scaf(i).dn
            dna.base_stap(i).dn = -1
        else if(dna.base_scaf(i).dn == -1) then
            dna.base_stap(i).up = -1
            dna.base_stap(i).dn = dna.base_scaf(i).up
        else
            dna.base_stap(i).up = dna.base_scaf(i).dn
            dna.base_stap(i).dn = dna.base_scaf(i).up
        end if

        ! Initialize crossover ID as negative value
        dna.base_scaf(i).xover = -1
        dna.base_stap(i).xover = -1

        ! Set across ID
        dna.base_scaf(i).across = dna.base_stap(i).id
        dna.base_stap(i).across = dna.base_scaf(i).id

        ! Set strand ID
        dna.base_scaf(i).strand = -1
        dna.base_stap(i).strand = -1

        ! Set position vector
        dna.base_scaf(i).pos(:) = 0.0d0
        dna.base_stap(i).pos(:) = 0.0d0
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.1. Initialize DNA data structure derived from basepairs"
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial bases in scaffold     : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of initial bases in stapple      : "//trim(adjustl(Int2Str(dna.n_base_stap)))
        write(i, "(a)")
    end do
end subroutine Route_Init_Base_Connectivity

! ---------------------------------------------------------------------------------------

! Set base position vector
! Last updated on Friday 19 Feb 2016 by Hyungmin
subroutine Route_Set_Base_Position(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    double precision :: ang_BP, ang_init, ang_start, ang_scaf, ang_stap
    double precision :: pos_scaf(3), pos_stap(3)
    double precision :: t(3,3), e(3,3), rot(3,3)
    integer :: i, j

    ! para_rad_helix   = 1.0d0,   Radius of DNA helices (nm)
    ! para_ang_minor   = 150.0d0, Angle of the minor groove (degree)
    ! para_ang_correct = 0.0d0,   Angle to be correct
    ! Angle rotated between neighboring base-pairs
    if(geom.sec.types == "square") then
        ang_BP = 360.0d0*3.0d0/32.0d0     ! 360/10.67 = 33.75
    else if(geom.sec.types == "honeycomb") then
        ang_BP = 360.0d0*2.0d0/21.0d0     ! 360/10.5  = 34.28
    end if

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.2. Set base position vectors using by parameters"
        call Space(i, 11)
        write(i, "(a)"), "* Section type                                : "//trim(geom.sec.types)//" lattice"
        call Space(i, 11)
        write(i, "(a)"), "* Starting base pair ID                       : "//trim(adjustl(Int2Str(para_start_bp_ID)))
        call Space(i, 11)
        write(i, "(a)"), "* Axial rise distance [nm]                    : "//trim(adjustl(Dble2Str(para_dist_bp)))
        call Space(i, 11)
        write(i, "(a)"), "* Radius of helix [nm]                        : "//trim(adjustl(Dble2Str(para_rad_helix)))
        call Space(i, 11)
        write(i, "(a)"), "* The Gap between helixes                     : "//trim(adjustl(Dble2Str(para_gap_helix)))
        call Space(i, 11)
        write(i, "(a)"), "* Angle of minor groove                       : "//trim(adjustl(Dble2Str(para_ang_minor)))
        call Space(i, 11)
        write(i, "(a)"), "* Angle correction factor                     : "//trim(adjustl(Dble2Str(para_ang_correct)))
    end do

    ! Loop for all finte nodes that are equal to base pairs
    do i = 1, mesh.n_node
        !
        ! Initialize initial angle that was referred with CANDO
        ! Find angle for scaffold and staple bases corresponding to basepair id
        !
        if(mod(mesh.node(i).sec, 2) == 0) then

            ! Positive z-direction
            if(geom.sec.types == "square") then
                ang_init = 180.0d0 - ang_BP / 2.0d0 + ang_BP
                ang_init = 180.0d0 + ang_BP / 2.0d0
            else if(geom.sec.types == "honeycomb") then
                if(geom.sec.dir ==  90) ang_init = 90.0d0
                if(geom.sec.dir == 150) ang_init = 90.0d0+60.0d0
                if(geom.sec.dir == -90) ang_init = 270.0d0
            end if

            ang_start = dmod(ang_init  + ang_BP*dble(para_start_bp_ID),  360.0d0)
            ang_start = dmod(ang_start + ang_BP*dble(mesh.node(i).bp-1), 360.0d0)
            ang_scaf  = dmod(ang_start + para_ang_correct,               360.0d0)
            ang_stap  = dmod(ang_scaf  + para_ang_minor,                 360.0d0)

        else if(mod(mesh.node(i).sec, 2) == 1) then

            ! Negative z-direction
            if(geom.sec.types == "square") then
                ang_init = 0.0d0 - ang_BP / 2.0d0 + ang_BP
                ang_init = 0.0d0 + ang_BP / 2.0d0
            else if(geom.sec.types == "honeycomb") then
                if(geom.sec.dir ==  90) ang_init = 270.0d0
                if(geom.sec.dir == 150) ang_init = 270.0d0+60.0d0
                if(geom.sec.dir == -90) ang_init = 90.0d0
            end if

            ang_start = dmod(ang_init  + ang_BP*dble(para_start_bp_ID),  360.0d0)
            ang_start = dmod(ang_start + ang_BP*dble(mesh.node(i).bp-1), 360.0d0)
            ang_scaf  = dmod(ang_start - para_ang_correct,               360.0d0)
            ang_stap  = dmod(ang_scaf  - para_ang_minor,                 360.0d0)

        end if

        ! Make negative angles to positive angles
        if(ang_scaf < 0.0d0 .or. ang_stap < 0.0d0) then
            ang_scaf = dmod(360.0d0 + ang_scaf, 360.0d0)
            ang_stap = dmod(360.0d0 + ang_stap, 360.0d0)
        end if

        ! Global Cartesian coordinate
        e(1,:) = 0.0d0; e(2,:) = 0.0d0; e(3,:) = 0.0d0
        e(1,1) = 1.0d0; e(2,2) = 1.0d0; e(3,3) = 1.0d0

        ! Local coordinate system
        t(1, 1:3) = geom.iniL(mesh.node(i).iniL).t(1, 1:3)
        t(2, 1:3) = geom.iniL(mesh.node(i).iniL).t(2, 1:3)
        t(3, 1:3) = geom.iniL(mesh.node(i).iniL).t(3, 1:3)

        ! The direction cosines of the coordinate transformation
        rot(1,1) = dot_product(e(1,:), t(1,:))
        rot(1,2) = dot_product(e(1,:), t(2,:))
        rot(1,3) = dot_product(e(1,:), t(3,:))
        rot(2,1) = dot_product(e(2,:), t(1,:))
        rot(2,2) = dot_product(e(2,:), t(2,:))
        rot(2,3) = dot_product(e(2,:), t(3,:))
        rot(3,1) = dot_product(e(3,:), t(1,:))
        rot(3,2) = dot_product(e(3,:), t(2,:))
        rot(3,3) = dot_product(e(3,:), t(3,:))

        ! Set base position vector for the scaffold strand
        pos_scaf(1) = para_rad_helix * 0.0d0
        pos_scaf(2) = para_rad_helix * dsin(-ang_scaf * pi/180.0d0)
        pos_scaf(3) = para_rad_helix * dcos(-ang_scaf * pi/180.0d0)
        dna.base_scaf(i).pos(1:3) = mesh.node(i).pos + matmul(rot, pos_scaf)

        ! Set base position vector for the staple strand
        pos_stap(1) = para_rad_helix * 0.0d0
        pos_stap(2) = para_rad_helix * dsin(-ang_stap * pi/180.0d0)
        pos_stap(3) = para_rad_helix * dcos(-ang_stap * pi/180.0d0)
        dna.base_stap(i).pos(1:3) = mesh.node(i).pos + matmul(rot, pos_stap)

        write(11, "(i20$)"), mesh.node(i).id
        write(11, "(a, a3$)"), " th node -> base ID : ", trim(adjustl(Int2Str(mesh.node(i).bp)))
        write(11, "(a, a6$)"), ", up : ",   trim(adjustl(Int2Str(mesh.node(i).up)))
        write(11, "(a, a6$)"), ", down : ", trim(adjustl(Int2Str(mesh.node(i).dn)))
        write(11, "(a, a2$)"), ", sec : ",  trim(adjustl(Int2Str(mesh.node(i).sec)))
        write(11, "(a, a3$)"), ", iniL : ", trim(adjustl(Int2Str(mesh.node(i).iniL)))
        write(11, "(a, a3$)"), ", croL : ", trim(adjustl(Int2Str(mesh.node(i).croL)))
        write(11, "(a, a6$)"), ", start ang : ", trim(adjustl(Dble2Str(ang_start)))
        write(11, "(a, a6$)"), ", scaf ang  : ", trim(adjustl(Dble2Str(ang_scaf)))
        write(11, "(a, a6 )"), ", stap ang  : ", trim(adjustl(Dble2Str(ang_stap)))
    end do

    ! Print progress
    write(0, "(a)"); write(11, "(a)")
end subroutine Route_Set_Base_Position

! ---------------------------------------------------------------------------------------

! Write scaffold or staple route for Chimera
! Last updated on Saturday 16 July 2016 by Hyungmin
subroutine Route_Chimera_Route(prob, geom, mesh, dna, step_route)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    character(6),   intent(in) :: step_route

    double precision, allocatable :: base_scaf(:,:),  base_stap(:,:)
    double precision, allocatable :: xover_scaf(:,:), xover_stap(:,:)
    double precision :: pos_1(3), pos_2(3), pos_x(3), pos_c(3), vec(3)
    integer :: i, j, base_1, base_2, base_x, node_1, node_2, node_x
    integer :: n_base_scaf, n_base_stap, n_xover_scaf, n_xover_stap
    logical :: f_axis, f_dir
    character(200) :: path

    if( (step_route == "route1" .and. para_write_601_1 == .true.) .or. &
        (step_route == "route2" .and. para_write_601_2 == .true.) .or. &
        (step_route == "route3" .and. para_write_601_3 == .true.) .or. &
        (step_route == "route4" .and. para_write_601_4 == .true.) .or. &
        (step_route == "route5" .and. para_write_601_5 == .true.) ) then
    else
        return
    end if

    ! Exception for Tecplot drawing
    if(step_route == "route3" .and. para_output_Tecplot == "on") then
        allocate(base_scaf (dna.n_base_scaf *2, 3))
        allocate(base_stap (dna.n_base_stap *2, 3))
        allocate(xover_scaf(dna.n_xover_scaf*4, 3))
        allocate(xover_stap(dna.n_xover_stap*4, 3))
        n_base_scaf  = 0
        n_base_stap  = 0
        n_xover_scaf = 0
        n_xover_stap = 0
    end if

    ! Set option
    f_axis = para_chimera_axis
    f_dir  = para_chimera_601_dir

    ! File open for route step
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=601, file=trim(path)//"_"//step_route//"_scaf.bild", form="formatted")
    open(unit=602, file=trim(path)//"_"//step_route//"_stap.bild", form="formatted")

    ! --------------------------------------------------
    ! For each base of the scaffold strand
    ! --------------------------------------------------
    do i = 1, dna.n_base_scaf
        if(dna.base_scaf(i).up == -1) cycle
        ! If upper ID exists
        ! current node
        !     *-------->*-------->*-------->*
        !  base_1    base_2
        base_1 = dna.base_scaf(i).id
        base_2 = dna.base_scaf(i).up

        ! Set node number
        node_1 = dna.base_scaf(base_1).node
        node_2 = dna.base_scaf(base_2).node

        ! Determine whether bases is in single nucleotide
        if(node_1 == -1 .or. node_2 == -1) then

            if(node_1 /= -1 .and. node_2 == -1) then
                ! Find the next base avoiding unpaired nucleotide
                do
                    if(dna.base_scaf(base_2).node /= -1) exit
                    base_2 = dna.base_scaf(base_2).up
                end do
            else
                cycle
            end if

            node_2 = dna.base_scaf(base_2).node

            pos_1(1:3) = mesh.node(node_1).pos(1:3)
            pos_2(1:3) = mesh.node(node_2).pos(1:3)

            write(601, "(a     )"), ".color red"
            write(601, "(a$    )"), ".cylinder "
            write(601, "(3f9.3$)"), pos_1(1:3)
            write(601, "(3f9.3$)"), pos_2(1:3)
            write(601, "(1f9.3 )"), 0.1d0

            if(step_route == "route3" .and. para_output_Tecplot == "on") then
                base_scaf(n_base_scaf + 1, 1:3) = pos_1(1:3)
                base_scaf(n_base_scaf + 2, 1:3) = pos_2(1:3)
                n_base_scaf = n_base_scaf + 2
            end if
            cycle
        end if

        pos_1(1:3) = mesh.node(node_1).pos(1:3)
        pos_2(1:3) = mesh.node(node_2).pos(1:3)

        ! Draw crossovers
        if(dna.base_scaf(base_1).xover /= -1) then

            base_x = dna.base_scaf(base_1).xover
            node_x = dna.base_scaf(base_x).node
            pos_x(1:3) = mesh.node(node_x).pos(1:3)

            write(601, "(a     )"), ".color tan"
            write(601, "(a$    )"), ".cylinder "
            write(601, "(3f9.3$)"), pos_1(1:3)
            write(601, "(3f9.3$)"), pos_x(1:3)
            write(601, "(1f9.3 )"), 0.1d0

            if(step_route == "route3" .and. para_output_Tecplot == "on") then
                xover_scaf(n_xover_scaf + 1, 1:3) = pos_1(1:3)
                xover_scaf(n_xover_scaf + 2, 1:3) = pos_x(1:3)
                n_xover_scaf = n_xover_scaf + 2
            end if

            if(dna.base_scaf(base_1).xover == dna.base_scaf(base_1).up) cycle
        end if

        ! Draw route
        write(601, "(a     )"), ".color steel blue"
        write(601, "(a$    )"), ".cylinder "
        write(601, "(3f9.3$)"), pos_1(1:3)
        write(601, "(3f9.3$)"), pos_2(1:3)
        write(601, "(1f9.3 )"), 0.1d0

        if(step_route == "route3" .and. para_output_Tecplot == "on") then
            base_scaf(n_base_scaf + 1, 1:3) = pos_1(1:3)
            base_scaf(n_base_scaf + 2, 1:3) = pos_2(1:3)
            n_base_scaf = n_base_scaf + 2
        end if
    end do

    ! --------------------------------------------------
    ! For each base of the staple strand
    ! --------------------------------------------------
    do i = 1, dna.n_base_stap
        if(dna.base_stap(i).up == -1) cycle
        ! If upper ID exists
        ! current node
        !     *-------->*-------->*-------->*
        !  base_1    base_2
        base_1 = dna.base_stap(i).id
        base_2 = dna.base_stap(i).up

        ! Set node number
        node_1 = dna.base_stap(base_1).node
        node_2 = dna.base_stap(base_2).node

        ! Determine whether bases is in Tn loop
        if(node_1 == -1 .or. node_2 == -1) then

            if(node_1 /= -1 .and. node_2 == -1) then
                ! Find the next base avoiding Tn loop
                do
                    if(dna.base_stap(base_2).node /= -1) exit
                    base_2 = dna.base_stap(base_2).up
                end do
            else
                cycle
            end if

            node_2 = dna.base_stap(base_2).node

            pos_1(1:3) = mesh.node(node_1).pos(1:3)
            pos_2(1:3) = mesh.node(node_2).pos(1:3)

            write(602, "(a     )"), ".color red"
            write(602, "(a$    )"), ".cylinder "
            write(602, "(3f9.3$)"), pos_1(1:3)
            write(602, "(3f9.3$)"), pos_2(1:3)
            write(602, "(1f9.3 )"), 0.1d0

            if(step_route == "route3" .and. para_output_Tecplot == "on") then
                base_stap(n_base_stap + 1, 1:3) = pos_1(1:3)
                base_stap(n_base_stap + 2, 1:3) = pos_2(1:3)
                n_base_stap = n_base_stap + 2
            end if
            cycle
        end if

        pos_1(1:3) = mesh.node(node_1).pos(1:3)
        pos_2(1:3) = mesh.node(node_2).pos(1:3)

        ! Draw crossovers
        if(dna.base_stap(base_1).xover /= -1) then

            base_x = dna.base_stap(base_1).xover
            node_x = dna.base_stap(base_x).node
            pos_x(1:3) = mesh.node(node_x).pos(1:3)

            write(602, "(a     )"), ".color tan"
            write(602, "(a$    )"), ".cylinder "
            write(602, "(3f9.3$)"), pos_1(1:3)
            write(602, "(3f9.3$)"), pos_x(1:3)
            write(602, "(1f9.3 )"), 0.1d0

            if(step_route == "route3" .and. para_output_Tecplot == "on") then
                xover_stap(n_xover_stap + 1, 1:3) = pos_1(1:3)
                xover_stap(n_xover_stap + 2, 1:3) = pos_x(1:3)
                n_xover_stap = n_xover_stap + 2
            end if

            if(dna.base_stap(base_1).xover == dna.base_stap(base_1).up) cycle
        end if

        write(602, "(a     )"), ".color orange"
        write(602, "(a$    )"), ".cylinder "
        write(602, "(3f9.3$)"), pos_1(1:3)
        write(602, "(3f9.3$)"), pos_2(1:3)
        write(602, "(1f9.3 )"), 0.1d0
        
        if(step_route == "route3" .and. para_output_Tecplot == "on") then
            base_stap(n_base_stap + 1, 1:3) = pos_1(1:3)
            base_stap(n_base_stap + 2, 1:3) = pos_2(1:3)
            n_base_stap = n_base_stap + 2
        end if
    end do

    ! Write direction arrow on end of edges
    if(f_dir == .true.) then
        do i = 1, geom.n_croL

            ! Even number of cross-section is the same direction of the crosL
            if(mod(geom.croL(i).sec, 2) == 0) then
                ! +z-direction
                pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
                pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
            else
                ! -z-direction
                pos_2(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
                pos_1(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
            end if

            ! Find direction vector for scaffolds
            vec(1:3) = Normalize_Vector(pos_2 - pos_1)

            write(601, "(a     )"), ".color red"
            write(601, "(a$    )"), ".arrow "
            write(601, "(3f8.2$)"), pos_1(1:3)
            write(601, "(3f8.2$)"), pos_1(1:3) + vec(1:3)*1.3d0
            write(601, "(3f8.2 )"), 0.12d0, 0.3d0, 0.6d0

            write(601, "(a$    )"), ".arrow "
            write(601, "(3f8.2$)"), pos_2(1:3) - vec(1:3)*1.3d0
            write(601, "(3f8.2$)"), pos_2
            write(601, "(3f8.2 )"), 0.12d0, 0.3d0, 0.6d0

            ! The opposite direction for staple strands
            pos_c(1:3) = pos_1(1:3)
            pos_1(1:3) = pos_2(1:3)
            pos_2(1:3) = pos_c(1:3)
            vec(1:3)   = Normalize_Vector(pos_2 - pos_1)

            write(602, "(a     )"), ".color red"
            write(602, "(a$    )"), ".arrow "
            write(602, "(3f8.2$)"), pos_1(1:3)
            write(602, "(3f8.2$)"), pos_1(1:3) + vec(1:3)*1.3d0
            write(602, "(3f8.2 )"), 0.12d0, 0.3d0, 0.6d0

            write(602, "(a$    )"), ".arrow "
            write(602, "(3f8.2$)"), pos_2(1:3) - vec(1:3)*1.3d0
            write(602, "(3f8.2$)"), pos_2
            write(602, "(3f8.2 )"), 0.12d0, 0.3d0, 0.6d0
        end do
    end if

    ! --------------------------------------------------
    ! Write starting point for scaffold route
    ! --------------------------------------------------
    !write(601, "(a     )"), ".color yellow"
    !write(601, "(a$    )"), ".sphere "
    !write(601, "(3f9.3$)"), mesh.node(1).pos(1:3)
    !write(601, "(1f9.3 )"), 0.25d0

    !pos_1(1:3) = mesh.node(1).pos(1:3)
    !pos_2(1:3) = mesh.node(mesh.node(1).up).pos(1:3)
    !vec(1:3)   = Normalize_Vector(pos_2 - pos_1)

    !write(601, "(a$    )"), ".arrow "
    !write(601, "(3f8.2$)"), pos_1(1:3)
    !write(601, "(3f8.2$)"), pos_1(1:3) + vec(1:3)*1.3d0
    !write(601, "(2f8.2 )"), 0.16d0, 0.4d0

    ! Write global axis
    if(f_axis == .true.) then
        do i = 0, 1
            write(601+i, "(a)"), ".translate 0.0 0.0 0.0"
            write(601+i, "(a)"), ".scale 0.5"
            write(601+i, "(a)"), ".color grey"
            write(601+i, "(a)"), ".sphere 0 0 0 0.5"      ! Center
            write(601+i, "(a)"), ".color red"             ! x-axis
            write(601+i, "(a)"), ".arrow 0 0 0 4 0 0 "
            write(601+i, "(a)"), ".color blue"            ! y-axis
            write(601+i, "(a)"), ".arrow 0 0 0 0 4 0 "
            write(601+i, "(a)"), ".color yellow"          ! z-axis
            write(601+i, "(a)"), ".arrow 0 0 0 0 0 4 "
        end do
    end if
    close(unit=601)
    close(unit=602)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(step_route /= "route3" .or. para_output_Tecplot /= "on") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=601, file=trim(path)//"_"//step_route//"_scaf.dat", form="formatted")
    open(unit=602, file=trim(path)//"_"//step_route//"_stap.dat", form="formatted")

    ! For scaffold bases
    write(601, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(601, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(601, "(a$)"), 'ZONE F = FEPOINT'
    write(601, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_scaf)))
    write(601, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_scaf/2)))
    write(601, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_scaf
        write(601, "(3f9.3$)"), base_scaf(i, 1:3)
        write(601, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_base_scaf/2
        write(601, "(2i7)"), 2*i, 2*i-1
    end do

    ! For crossover for scaffold strand
    write(601, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(601, "(a$)"), 'ZONE F = FEPOINT'
    write(601, "(a$)"), ', N='//trim(adjustl(Int2Str(n_xover_scaf)))
    write(601, "(a$)"), ', E='//trim(adjustl(Int2Str(n_xover_scaf/2)))
    write(601, "(a )"), ', ET=LINESEG'

    do i = 1, n_xover_scaf
        write(601, "(3f9.3$)"), xover_scaf(i, 1:3)
        write(601, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_xover_scaf/2
        write(601, "(2i7)"), 2*i, 2*i-1
    end do

    ! For staple bases
    write(602, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(602, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(602, "(a$)"), 'ZONE F = FEPOINT'
    write(602, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_stap)))
    write(602, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_stap/2)))
    write(602, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_stap
        write(602, "(3f9.3$)"), base_stap(i, 1:3)
        write(602, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_base_stap/2
        write(602, "(2i7)"), 2*i, 2*i-1
    end do

    ! For crossover for staple strand
    write(602, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(602, "(a$)"), 'ZONE F = FEPOINT'
    write(602, "(a$)"), ', N='//trim(adjustl(Int2Str(n_xover_stap)))
    write(602, "(a$)"), ', E='//trim(adjustl(Int2Str(n_xover_stap/2)))
    write(602, "(a )"), ', ET=LINESEG'

    do i = 1, n_xover_stap
        write(602, "(3f9.3$)"), xover_stap(i, 1:3)
        write(602, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, n_xover_stap/2
        write(602, "(2i7)"), 2*i, 2*i-1
    end do

    ! Deallocate memory
    if(allocated(base_scaf))  deallocate(base_scaf)
    if(allocated(base_stap))  deallocate(base_stap)
    if(allocated(xover_scaf)) deallocate(xover_scaf)
    if(allocated(xover_stap)) deallocate(xover_stap)

    close(unit=601)
    close(unit=602)
end subroutine Route_Chimera_Route

! ---------------------------------------------------------------------------------------

! Connect strands at the junction by unpaired nucleotides
! Last updated on Tue 21 Mar 2017 by Hyungmin
subroutine Route_Connect_Strand_Junc(geom, bound, mesh, dna)
    type(GeomType),  intent(in)    :: geom
    type(BoundType), intent(in)    :: bound
    type(MeshType),  intent(in)    :: mesh
    type(DNAType),   intent(inout) :: dna

    integer, allocatable, dimension(:,:) :: conn
    integer, allocatable, dimension(:,:) :: join

    integer :: node_cur, iniL_cur, croL_cur, sec_cur
    integer :: node_com, iniL_com, croL_com, sec_com
    integer :: n_col, start, count, n_scaf, n_stap, n_face, n_conn
    integer :: i, j, k, m, n, n_add_base_scaf, n_add_base_stap
    logical :: b_conn
    character(10) :: dir_cur, dir_com

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.3. Connect strands at the junction"
        call Space(i, 11)
        write(i, "(a)"), "* Vertex design method                        : "//trim(para_vertex_design)// " vertex"
        call Space(i, 11)
        write(i, "(a)"), "* The number of junctions                     : "//trim(adjustl(Int2Str(bound.n_junc)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of rows on cross-section         : "//trim(adjustl(Int2Str(geom.sec.maxR-geom.sec.minR+1)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of cols on cross-section         : "//trim(adjustl(Int2Str(geom.sec.maxC-geom.sec.minC+1)))
        call Space(i, 11)
        write(i, "(a)"), "* Detailed information on connection"
    end do

    ! The number of added unpaired nucleotides
    n_add_base_scaf = 0
    n_add_base_stap = 0

    ! Loop for junction data
    do i = 1, bound.n_junc

        ! Print progress bar
        call Mani_Progress_Bar(i, bound.n_junc)

        ! Allocate conn arrary to store node to be connected
        ! conn(j, 1) : current node
        ! conn(j, 2) : node to be joined with conn(j, 1)
        allocate(conn(geom.n_sec*bound.junc(i).n_arm, 2))

        conn(:, 1) = bound.junc(i).conn(:, 1)
        conn(:, 2) = bound.junc(i).conn(:, 2)

        ! Print progress
        write(11, "(i20$)"), i
        write(11, "(a$  )"), " th junc -> # of nodes to be connected : "
        write(11, "(i7  )"), geom.n_sec*bound.junc(i).n_arm

        do j = 1, geom.n_sec*bound.junc(i).n_arm
            write(11, "(i30, a$)"), conn(j, 1), " th node --> "
            write(11, "(i7, a  )"), conn(j, 2), " th node"
        end do

        ! ==================================================
        ! Build join without any duplicate data
        ! ==================================================
        ! Allocate join
        allocate(join(geom.n_sec*bound.junc(i).n_arm/2, 2))

        ! Find duplicate data
        n_conn = 0
        do j = 1, geom.n_sec*bound.junc(i).n_arm

            ! Compare with previous storaged data
            b_conn = .false.
            do k = 1, n_conn
                if(join(k, 1) == conn(j, 1)) then
                    if(join(k, 2) == conn(j, 2)) then
                        b_conn = .true.
                    end if
                else if(join(k, 1) == conn(j, 2)) then
                    if(join(k, 2) == conn(j, 1)) then
                        b_conn = .true.
                    end if
                end if
            end do

            ! If there is no entity in existing array
            if(b_conn == .false.) then
                n_conn = n_conn + 1
                join(n_conn, 1:2) = conn(j, 1:2)
            end if
        end do
        deallocate(conn)

        ! ==================================================
        !
        ! Connect nodes between join(j,1) and join(j,2)
        ! Size j : geom.n_sec * bound.junc(i).n_arm / 2
        !
        ! ==================================================
        ! Connect node from join1 to join2
        do j = 1, geom.n_sec * bound.junc(i).n_arm / 2
            !
            ! *--------->*         or    *<---------*
            ! node_cur   node_com        node_cur   node_com
            !
            ! Set current and comparing node to be connected
            node_cur = join(j, 1)
            node_com = join(j, 2)
//  여기 부터
            ! ==================================================
            !
            ! Neighbor connection
            ! For scaffold - connect with/without unpaired nucleotides
            ! For staple - connect with poly Tn loop
            !
            ! ==================================================
            if(geom.sec.conn(mesh.node(node_cur).sec + 1) == -1) then
                ! Connect scaffold with/without unpaired nucleotides
                n_add_base_scaf = n_add_base_scaf + &
                    Route_Connect_Scaf(mesh, dna, node_cur, node_com)

                ! Connect staple with poly Tn loop
                n_add_base_stap = n_add_base_stap + &
                    Route_Connect_Stap(mesh, dna, node_cur, node_com)
            end if
        end do

        deallocate(join)
    end do

    

    
    stop
    
    
    ! ==================================================
    !
    ! Self connection
    ! Scaffold strand - connection in the same section
    !
    ! ==================================================
    ! Loop for junction
    do i = 1, bound.n_junc
        do j = 1, bound.junc(i).n_arm

            do m = 1, geom.n_sec
                node_cur = bound.junc(i).node(j, m)
                sec_cur  = mesh.node(node_cur).sec

                do n = m + 1, geom.n_sec
                    node_com = bound.junc(i).node(j, n)
                    sec_com  = mesh.node(node_com).sec

                    if(sec_cur == sec_com) cycle
                    if(geom.sec.conn(sec_cur+1) /= sec_com) cycle

                    ! Connect scaffold strand in the same section without unpaired nucleotides
                    n_add_base_scaf = n_add_base_scaf + &
                        Route_Connect_Scaf(mesh, dna, node_cur, node_com)
                end do
            end do
        end do
    end do

    ! Check the closed scaffold strand
    do i = 1, dna.n_base_scaf
        if(dna.base_scaf(i).up == -1 .or. dna.base_scaf(i).dn == -1) then
            write(0, "(a$)"), "Error - Not circular scaffold : "
            write(0, "(a )"), "Route_Reconnect_Junction"
            stop
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The number of bases in scaffold strand      : "&
            //trim(adjustl(Int2Str(dna.n_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of bases in staple strand        : "&
            //trim(adjustl(Int2Str(dna.n_base_stap)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of unpaired scaffold nucleotides : "&
            //trim(adjustl(Int2Str(n_add_base_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of added bases in poly Tn loop   : "&
            //trim(adjustl(Int2Str(n_add_base_stap)))
    end do

    do j = 1, dna.n_base_scaf
        write(11, "(i20$  )"), dna.base_scaf(j).id
        write(11, "(a, i7$)"), " th scaf -> node # : ", dna.base_scaf(j).node
        write(11, "(a, i7$)"), ", up # : ",             dna.base_scaf(j).up
        write(11, "(a, i7$)"), ", down # : ",           dna.base_scaf(j).dn
        write(11, "(a, i7$)"), ", xover # : ",          dna.base_scaf(j).xover
        write(11, "(a, i7 )"), ", across # : ",         dna.base_scaf(j).across
    end do

    do j = 1, dna.n_base_stap
        write(11, "(i20$  )"), dna.base_stap(j).id
        write(11, "(a, i7$)"), " th stap -> node # : ", dna.base_stap(j).node
        write(11, "(a, i7$)"), ", up # : ",             dna.base_stap(j).up
        write(11, "(a, i7$)"), ", down # : ",           dna.base_stap(j).dn
        write(11, "(a, i7$)"), ", xover # : ",          dna.base_stap(j).xover
        write(11, "(a, i7 )"), ", across # : ",         dna.base_stap(j).across
    end do

    write(0, "(a)"); write(11, "(a)")

end subroutine Route_Connect_Strand_Junc

! ---------------------------------------------------------------------------------------

! Connect scaffold with/without unpaired nucleotides
! Last updated on Saturday 13 August 2016 by Hyungmin
function Route_Connect_Scaf(mesh, dna, node_cur, node_com) result(n_add_base)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: node_cur
    integer,        intent(in)    :: node_com

    double precision :: pos_cur(3), pos_com(3), pos(3), length
    integer :: i, cur, com, ttt, n_add_base

    ! Set current and comparing nodes
    cur = node_cur
    com = node_com

    ! Find position vector of the curent and comparing bases
    pos_cur(1:3) = dna.base_scaf(cur).pos(1:3)
    pos_com(1:3) = dna.base_scaf(com).pos(1:3)

    ! Only for modified neighbor connection
    !if(para_unpaired_scaf == "on") then
    if(para_unpaired_scaf == "on" .and. mesh.node(cur).conn == 3 .and. mesh.node(com).conn == 3) then

        ! mesh.node(cur).conn = 3 means modified neighbor connection
        length     = Size_Vector(pos_cur - pos_com)
        n_add_base = idnint(length/para_dist_pp) - 1

        if(para_unpaired_square == "on") then
            n_add_base = n_add_base + 2
        end if
    else if(para_unpaired_square == "on" .and. mesh.node(cur).conn == 1 .and. mesh.node(com).conn == 1) then

        ! mesh.node(cur).conn = 2 means original neighbor connection
        n_add_base = 2
    else

        n_add_base = 0
    end if

    if(dna.base_scaf(cur).up == -1 .and. dna.base_scaf(com).dn == -1) then
        !
        ! --------------------------------------------------
        !
        ! If the current node is inward to the junction
        ! ------->*====+====>*-------->
        !        cur  junc  com
        !
        ! --------------------------------------------------
        ! Set base connectivity of the scaffold strand
        dna.base_scaf(cur).up = dna.base_scaf(com).id
        dna.base_scaf(com).dn = dna.base_scaf(cur).id

        ! Add bases for depending on the distance between two bases
        do i = 1, n_add_base

            ! Find additional base position
            ! *---#---#---#---#---*     ex)4 bases
            pos(1:3) = (dble(i)*pos_com(1:3) + dble(n_add_base+1-i)*pos_cur(1:3)) / dble(n_add_base+1)

            ! Add one neucleotide
            call Route_Add_Nucleotide(dna.base_scaf, dna.n_base_scaf, pos(1:3))

            ! Connect scaffold strand
            ttt = dna.n_base_scaf
            dna.base_scaf(ttt).up = dna.base_scaf(com).id
            dna.base_scaf(ttt).dn = dna.base_scaf(cur).id
            dna.base_scaf(cur).up = dna.base_scaf(ttt).id
            dna.base_scaf(com).dn = dna.base_scaf(ttt).id

            ! Update current node
            cur = ttt
        end do
    else if(dna.base_scaf(cur).dn == -1 .and. dna.base_scaf(com).up == -1) then
        !
        !
        ! --------------------------------------------------
        !
        ! If the current node is inward to the junction
        ! <-------*<====+====*<--------
        !        cur  junc   com
        !
        ! --------------------------------------------------
        ! Set base connectivity of the scaffold strand
        dna.base_scaf(cur).dn = dna.base_scaf(com).id
        dna.base_scaf(com).up = dna.base_scaf(cur).id

        ! Add bases depending on the distance between two bases
        do i = 1, n_add_base

            ! Find additional base position
            pos(1:3) = (dble(i)*pos_com(1:3) + dble(n_add_base+1-i)*pos_cur(1:3)) / dble(n_add_base+1)

            ! Add one neucleotide
            call Route_Add_Nucleotide(dna.base_scaf, dna.n_base_scaf, pos(1:3))

            ! connect scaffold strand
            ttt = dna.n_base_scaf
            dna.base_scaf(ttt).up = dna.base_scaf(cur).id
            dna.base_scaf(ttt).dn = dna.base_scaf(com).id
            dna.base_scaf(cur).dn = dna.base_scaf(ttt).id
            dna.base_scaf(com).up = dna.base_scaf(ttt).id

            ! Update current node
            cur = ttt
        end do
    else
        write(0, "(a$)"), "Error - Wrong current ID for connection : "
        write(0, "(a )"), "Route_Connect_Scaf"
        stop
    end if
end function Route_Connect_Scaf

! ---------------------------------------------------------------------------------------

! Connect staple nucleotides by buiding poly T loop
! Last updated on Monday 08 August 2016 by Hyungmin
function Route_Connect_Stap(mesh, dna, node_cur, node_com) result(n_poly_Tn)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: node_cur
    integer,        intent(in)    :: node_com

    double precision :: pos_cur(3), pos_com(3), pos(3), length
    integer :: i, cur, com, ttt, n_poly_Tn

    ! Set current and comparing nodes
    cur = node_cur
    com = node_com

    ! Find position vector of the curent and comparing bases
    pos_cur(1:3) = dna.base_stap(cur).pos(1:3)
    pos_com(1:3) = dna.base_stap(com).pos(1:3)

    ! Add bases for Tn loop depending on distance between two bases
    if(para_n_base_tn == -1) then

        length    = Size_Vector(pos_cur - pos_com)
        n_poly_Tn = idnint(length/para_dist_pp) - 1

        ! To make region for sequence design
        if(n_poly_Tn == 0 ) n_poly_Tn = 1
        !if(n_poly_Tn == 6 ) n_poly_Tn = 7
    else
        n_poly_Tn = para_n_base_tn
    end if

    if(para_unpaired_square == "on") n_poly_Tn = n_poly_Tn + 2

    ! For reference connection without change
    if(mesh.node(cur).conn == 1 .and. mesh.node(com).conn == 1) then
        !n_poly_Tn = 7
    end if

    if(dna.base_stap(cur).up == -1 .and. dna.base_stap(com).dn == -1) then
        !
        ! --------------------------------------------------
        !
        ! If the current node is inward to the junction
        ! ------->*====+====>*-------->
        !        cur  junc  com
        !
        ! --------------------------------------------------
        ! Set base connectivity of the staple strand
        dna.base_stap(cur).up = dna.base_stap(com).id
        dna.base_stap(com).dn = dna.base_stap(cur).id

        ! Add bases for Tn loop depending on the distance between two bases
        do i = 1, n_poly_Tn

            ! Find additional base position
            ! *---#---#---#---#---*     ex)4 bases for Tn loop
            pos(1:3) = (dble(i)*pos_com(1:3) + dble(n_poly_Tn+1-i)*pos_cur(1:3)) / dble(n_poly_Tn+1)

            ! Add one neucleotide to make poly Tn loop(staple)
            call Route_Add_Nucleotide(dna.base_stap, dna.n_base_stap, pos(1:3))

            ! Connect staple strand
            ttt = dna.n_base_stap
            dna.base_stap(ttt).up = dna.base_stap(com).id
            dna.base_stap(ttt).dn = dna.base_stap(cur).id
            dna.base_stap(cur).up = dna.base_stap(ttt).id
            dna.base_stap(com).dn = dna.base_stap(ttt).id

            ! Update current node
            cur = ttt
        end do
    else if(dna.base_stap(cur).dn == -1 .and. dna.base_stap(com).up == -1) then
        !
        !
        ! --------------------------------------------------
        !
        ! If the current node is inward to the junction
        ! <-------*<====+====*<--------
        !        cur  junc   com
        !
        ! --------------------------------------------------
        ! Set base connectivity of the staple strand
        dna.base_stap(cur).dn = dna.base_stap(com).id
        dna.base_stap(com).up = dna.base_stap(cur).id

        ! Add bases for Tn loop depending on the distance between two bases
        do i = 1, n_poly_Tn

            ! Find additional base position
            pos(1:3) = (dble(i)*pos_com(1:3) + dble(n_poly_Tn+1-i)*pos_cur(1:3)) / dble(n_poly_Tn+1)

            ! Add one neucleotide to make poly Tn loop(staple)
            call Route_Add_Nucleotide(dna.base_stap, dna.n_base_stap, pos(1:3))

            ! connect staple strand
            ttt = dna.n_base_stap
            dna.base_stap(ttt).up = dna.base_stap(cur).id
            dna.base_stap(ttt).dn = dna.base_stap(com).id
            dna.base_stap(cur).dn = dna.base_stap(ttt).id
            dna.base_stap(com).up = dna.base_stap(ttt).id

            ! Update current node
            cur = ttt
        end do
    else
        write(0, "(a$)"), "Error - Wrong current ID for connection : "
        write(0, "(a )"), "Route_Connect_Staple"
        stop
    end if
end function Route_Connect_Stap

! ---------------------------------------------------------------------------------------

! Add one neucleotide to make single strand(scaffold) or poly Tn loop(staple)
! Last updated on Thuesday 9 August 2016 by Hyungmin
subroutine Route_Add_Nucleotide(base, n_base, pos)
    type(BaseType), allocatable, intent(inout) :: base(:)
    integer,                     intent(inout) :: n_base
    double precision,            intent(in)    :: pos(3)

    type(BaseType), allocatable, dimension(:) :: t_base

    integer :: i

    ! Copy allocation dynamic array of base to t_base
    call move_alloc(base, t_base)

    ! Increase the number of nucleotides
    n_base = n_base + 1

    ! Reallocate global base
    allocate(base(n_base))

    ! Copy data from temporary array to original one
    call Mani_Copy_BaseType(t_base, base, n_base - 1)

    ! For new added base
    base(n_base).id     = n_base
    base(n_base).node   = -1
    base(n_base).up     = -1
    base(n_base).dn     = -1
    base(n_base).xover  = -1
    base(n_base).across = -1
    base(n_base).strand = -1
    base(n_base).pos(:) = pos(1:3)

    ! Deallocate
    deallocate(t_base)
end subroutine Route_Add_Nucleotide

! ---------------------------------------------------------------------------------------

! Set the strand ID in scaffold
! Last updated on Wednesday 3 May 2016 by Hyungmin
subroutine Route_Set_Strand_ID_Scaf(dna)
    type(DNAType), intent(inout) :: dna

    logical, allocatable, dimension(:) :: visit
    integer :: i, n_strand, start, search, count

    ! Print information
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.4. Set strand ID on bases of the scaffold strand"
        call Space(i, 11)
        write(i, "(a)"), "* Set strand ID into dna.base_scaf data"
    end do

    ! Allocate and initialize variable to check visiting
    allocate(visit(dna.n_base_scaf))
    do i = 1, dna.n_base_scaf
        visit(i) = .false.
    end do

    ! Set strand number to dna.base_scaf(i).strand
    n_strand = 0
    do i = 1, dna.n_base_scaf

        ! If the base is visited, skip the loop
        if(visit(i) == .true.) cycle

        n_strand = n_strand + 1
        start    = dna.base_scaf(i).id
        search   = start
        visit(i) = .true.

        ! Check strand number
        dna.base_scaf(i).strand = n_strand

        count = 0
        ! Infinite loop to check original base id
        do
            search = dna.base_scaf(search).up
            visit(search) = .true.

            ! If the original base is found
            if(start == search) exit

            ! Check strand number
            dna.base_scaf(search).strand = n_strand

            count = count + 1
            if(count > dna.n_base_scaf) stop
        end do
    end do

    dna.n_scaf = n_strand

    ! Deallocate array
    deallocate(visit)

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The number of scaffold strands              : "//trim(adjustl(Int2Str(dna.n_scaf)))
        write(i, "(a)")
    end do
end subroutine Route_Set_Strand_ID_Scaf

! ---------------------------------------------------------------------------------------

! Find center scaffold crossovers with splitted algorithm
! Last updated on Saturday 9 June 2016 by Hyungmin
subroutine Route_Find_Centered_Scaf_Xover(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Crossover based on cross-sectional edges
    type :: CroLType
        integer :: max_bp, min_bp
        logical, allocatable :: b_xover(:)
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, k, ave_bp, sec_cur, sec_com, croL_cur, croL_com, id_bp
    integer :: up_cur, up_com, dn_cur, dn_com, node_cur, node_com
    integer :: n_gap, max_gap, max_cur, max_com, min_bp1, max_bp1, min_bp2, max_bp2
    integer :: iniL, pre_iniL, split_cur, split_com, step
    logical :: b_nei, b_nei_up, b_nei_dn, b_fail

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.5. Find centered scaffold crossovers"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of cross-sectional edges         : "
        write(i, "(a )"), trim(adjustl(Int2Str(geom.n_croL)))
        call Space(i, 11)
        write(i, "(a$)"), "* The number of sections                      : "
        write(i, "(a )"), trim(adjustl(Int2Str(geom.n_sec)))
    end do

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    do i = 1, geom.n_croL
        allocate(croL(i).b_xover(geom.n_sec))
        croL(i).b_xover(1:geom.n_sec) = .false.
        croL(i).max_bp = -999999
        croL(i).min_bp =  999999
    end do

    ! Set maximum and minimum base ID in cross-sectional edges
    do i = 1, mesh.n_node
        croL_cur = mesh.node(i).croL
        id_bp    = mesh.node(i).bp

        ! Set maximum base ID
        if(croL(croL_cur).max_bp < id_bp) then
            croL(croL_cur).max_bp = id_bp
        end if

        ! Set minimum base ID
        if(croL(croL_cur).min_bp > id_bp) then
            croL(croL_cur).min_bp = id_bp
        end if
    end do

    ! Find centered scaffold crossovers with spiting algorithm
    dna.n_xover_scaf = 0
    pre_iniL         = 0
    do i = 1, mesh.n_node

        ! Print progress bar
        call Mani_Progress_Bar(i, mesh.n_node)

        ! Comparing node
        do j = i + 1, mesh.n_node

            ! It should be skipped if satisfied
            if(mesh.node(i).bp   /= mesh.node(j).bp  ) cycle
            if(mesh.node(i).iniL /= mesh.node(j).iniL) cycle
            if(mesh.node(i).sec  == mesh.node(j).sec ) cycle

            ! Find section, cross line and section ID
            sec_cur  = mesh.node(i).sec
            sec_com  = mesh.node(j).sec
            croL_cur = mesh.node(i).croL
            croL_com = mesh.node(j).croL
            id_bp    = mesh.node(i).bp
            min_bp1  = croL(croL_cur).min_bp
            max_bp1  = croL(croL_cur).max_bp
            min_bp2  = croL(croL_com).min_bp
            max_bp2  = croL(croL_com).max_bp
            ave_bp   = (min_bp1 + max_bp1) / 2 - 2

            ! To put double crossover at center
            if(id_bp < ave_bp) cycle
            if(id_bp > max_bp1 - para_gap_xover_bound_scaf) cycle
            if(croL(croL_cur).b_xover(sec_com+1) == .true.) cycle
            if(croL(croL_com).b_xover(sec_cur+1) == .true.) cycle

            ! Check crossover
            b_nei = Section_Connection_Scaf(geom, sec_cur, sec_com, id_bp)

            ! Check neighbor crossover position
            if(b_nei == .false.) cycle

            node_cur = mesh.node(i).id
            node_com = mesh.node(j).id

            ! ==================================================
            !
            ! Build spllited crossover
            !
            ! ==================================================
            max_gap   = n_gap
            max_cur   = node_cur
            max_com   = node_com
            split_cur = node_cur
            split_com = node_com
            iniL      = mesh.node(i).iniL

            !print *, sec_cur, sec_com, id_bp, id_bp + para_start_bp_ID - 1

            if((pre_iniL == 0 .or. iniL /= pre_iniL) .and. sec_cur == 0) then

                ! Consistent pre-defined scaffold crossover for first connection, 0-5
                if((prob.sel_sec /= 1) .and. (sec_cur /= 0 .or. sec_com /= 5)) cycle

                !
                ! --------------------------------------------------
                !
                ! First crossover at each edge
                !
                ! --------------------------------------------------
                pre_iniL = iniL

                ! Set centered crossovers when DX tile, otherwise, spliting
                if(geom.sec.maxR == 1 .and. geom.sec.maxC == 2) then
                    ! Centered crossover for DX tile
                    if(para_set_xover_scaf == "center") step = 1
                    if(para_set_xover_scaf == "split")  step = 5
                else
                    ! Other cross-sections depend on parameter
                    if(para_set_xover_scaf == "center") step = 1
                    if(para_set_xover_scaf == "split")  step = 3
                end if

                ! Split crossover to avoid existing one, return crossover position
                ! 1-going down further, 2-going down, 3-going up further, 4-going up, 5-center
                b_fail = Route_Split_Centered_Scaf_Xover&
                    (geom, mesh, split_cur, split_com, min_bp1, max_bp1, min_bp2, max_bp2, step)

                ! Check the gap to the vertex boundary
                if(step == 3 .and. b_fail == .true.) then
                    b_fail = Route_Split_Centered_Scaf_Xover&
                        (geom, mesh, split_cur, split_com, min_bp1, max_bp1, min_bp2, max_bp2, step + 1)
                end if

                ! Make centered crossover
                node_cur = split_cur
                node_com = split_com
            else
                !
                ! --------------------------------------------------
                !
                ! Split crossover from center
                !
                ! --------------------------------------------------
                do k = 1, 5

                    split_cur = node_cur
                    split_com = node_com

                    ! Split crossover to avoid existing one, return crossover position
                    ! 1-going down further, 2-going down, 3-going up further, 4-going up, 5-center
                    b_fail = Route_Split_Centered_Scaf_Xover&
                        (geom, mesh, split_cur, split_com, min_bp1, max_bp1, min_bp2, max_bp2, k)

                    ! Check gap distance
                    n_gap = Route_Check_Nei_Xover(geom, mesh, dna, split_cur, split_com)

                    ! Check the criteria
                    if(b_fail == .false. .and. n_gap == para_gap_xover_two_scaf) then

                        ! Exit loop
                        exit
                    else if(b_fail == .false. .and. n_gap /= para_gap_xover_two_scaf) then

                        ! Check the gap distance
                        if(n_gap >= max_gap) then
                            max_gap = n_gap
                            max_cur = split_cur
                            max_com = split_com
                        end if
                    end if

                    ! For max gap crossover for last iteration
                    if(k == 5) then
                        split_cur = max_cur
                        split_com = max_com
                    end if
                end do

                node_cur = split_cur
                node_com = split_com
            end if

            ! Find upper and downward crossovers
            croL(croL_cur).b_xover(sec_com+1) = .true.

            ! Exception module by 9/27/2016 to treat error : 51_1_1_beveled
            if(node_cur == -1 .or. node_com == -1) cycle

            up_cur = mesh.node(node_cur).up
            up_com = mesh.node(node_com).dn

            ! Exception module by 9/27/2016 to treat error : 51_1_1_beveled
            if(up_cur == -1 .or. up_com == -1) cycle

            sec_cur  = mesh.node(up_cur).sec
            sec_com  = mesh.node(up_com).sec
            id_bp    = mesh.node(up_cur).bp
            b_nei_up = Section_Connection_Scaf(geom, sec_cur, sec_com, id_bp)

            dn_cur = mesh.node(node_cur).dn
            dn_com = mesh.node(node_com).up

            ! Exception module by 9/27/2016 to treat error : 51_1_1_beveled
            if(dn_cur == -1 .or. dn_com == -1) cycle

            sec_cur  = mesh.node(dn_cur).sec
            sec_com  = mesh.node(dn_com).sec
            id_bp    = mesh.node(dn_cur).bp
            b_nei_dn = Section_Connection_Scaf(geom, sec_cur, sec_com, id_bp)

            ! Set current and previous or next crossovers
            dna.n_xover_scaf = dna.n_xover_scaf + 2
            dna.base_scaf(node_cur).xover = dna.base_scaf(node_com).id
            dna.base_scaf(node_com).xover = dna.base_scaf(node_cur).id

            if(b_nei_up == .true.) then
                dna.base_scaf(up_cur).xover = dna.base_scaf(up_com).id
                dna.base_scaf(up_com).xover = dna.base_scaf(up_cur).id
            else if(b_nei_dn == .true.) then
                dna.base_scaf(dn_cur).xover = dna.base_scaf(dn_com).id
                dna.base_scaf(dn_com).xover = dna.base_scaf(dn_cur).id
            end if
        end do
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The number of centered scaffold crossovers  : "&
            //trim(adjustl(Int2Str(dna.n_xover_scaf)))
        write(i, "(a)")
    end do

    ! Deallocate memory
    do i = 1, geom.n_croL
        deallocate(croL(i).b_xover)
    end do
    deallocate(croL)
end subroutine Route_Find_Centered_Scaf_Xover

! ---------------------------------------------------------------------------------------

! Check whether nighbor crossover exists or not
! Last updated on Friday 8 June 2016 by Hyungmin
function Route_Check_Nei_Xover(geom, mesh, dna, cur, com) result(n_gap)
    type(GeomType), intent(in) :: geom    
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: cur
    integer,        intent(in) :: com

    integer :: i, n_gap, sec1, sec2, up_node1, up_node2, dn_node1, dn_node2
    logical :: b_nei_up, b_nei_dn

    sec1 = mesh.node(cur).sec
    sec2 = mesh.node(com).sec

    ! Exception
    if(mesh.node(cur).up == -1 .or. mesh.node(cur).dn == -1) then
        n_gap = 0
        return
    end if

    ! Check neighbor crossover position
    b_nei_up = Section_Connection_Scaf(geom, sec1, sec2, mesh.node(mesh.node(cur).up).bp)
    b_nei_dn = Section_Connection_Scaf(geom, sec1, sec2, mesh.node(mesh.node(cur).dn).bp)

    if(b_nei_up == .true.) then
        up_node1 = mesh.node(cur).up
        up_node2 = mesh.node(com).dn
        dn_node1 = cur
        dn_node2 = com
    else if(b_nei_dn == .true.) then
        up_node1 = cur
        up_node2 = com
        dn_node1 = mesh.node(cur).dn
        dn_node2 = mesh.node(com).up
    end if

    ! Check neighboring crossovers
    n_gap = 0
    do i = 1, para_gap_xover_two_scaf

        ! If the node approches the boundary
        if( up_node1 == -1 .or. up_node2 == -1 .or. &
            dn_node1 == -1 .or. dn_node2 == -1 ) then
            exit
        end if

        up_node1 = mesh.node(up_node1).up
        up_node2 = mesh.node(up_node2).dn

        dn_node1 = mesh.node(dn_node1).dn
        dn_node2 = mesh.node(dn_node2).up

        ! If the node approches the boundary
        if( up_node1 == -1 .or. up_node2 == -1 .or. &
            dn_node1 == -1 .or. dn_node2 == -1 ) then
            exit
        end if

        ! If there is another crossover
        if( dna.base_scaf(up_node1).xover /= -1 .or. &
            dna.base_scaf(up_node2).xover /= -1 .or. &
            dna.base_scaf(dn_node1).xover /= -1 .or. &
            dna.base_scaf(dn_node2).xover /= -1 ) then
            exit
        end if

        ! Update
        n_gap = n_gap + 1
    end do
end function Route_Check_Nei_Xover

! ---------------------------------------------------------------------------------------

! Split crossover to avoid existing one, return crossover position
! Last updated on Monday 15 August 2016 by Hyungmin
function Route_Split_Centered_Scaf_Xover(geom, mesh, cur, com, min_bp1, max_bp1, min_bp2, max_bp2, step) result(b_fail)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    integer,        intent(inout) :: cur
    integer,        intent(inout) :: com
    integer,        intent(in)    :: min_bp1, max_bp1
    integer,        intent(in)    :: min_bp2, max_bp2
    integer,        intent(in)    :: step

    integer :: node1, node2, sec1, sec2
    integer :: id_bp1, id_bp2, n_skip, n_cross
    logical :: b_fail
    character(10) :: direction

    ! node1, sec1 - current node and its section ID
    ! node2, sec2 - comparing node and its section ID
    b_fail = .false.
    node1  = cur
    node2  = com
    sec1   = mesh.node(node1).sec
    sec2   = mesh.node(node2).sec

    ! Set direction and number of crossovers to skip
    if(para_set_xover_scaf == "split") then

        if(step == 1) then
            n_cross = 3; direction = "down"   ! Going down further
        else if(step == 2) then
            n_cross = 1; direction = "down"   ! Going down
        else if(step == 3) then
            n_cross = 3; direction = "up"     ! Going up further
        else if(step == 4) then
            n_cross = 1; direction = "up"     ! Going up
        else if(step == 5) then
            n_cross = 0; direction = "center" ! Remain at center
        end if
    else if(para_set_xover_scaf == "center") then

        if(step == 1) then
            n_cross = 0; direction = "center" ! Remain at center
        else if(step == 2) then
            n_cross = 1; direction = "down"   ! Going down
        else if(step == 3) then
            n_cross = 1; direction = "up"     ! Going up
        else if(step == 4) then
            n_cross = 3; direction = "down"   ! Going down further
        else if(step == 5) then
            n_cross = 3; direction = "up"     ! Going up further
        end if
    end if

    ! 1-going down further, 2-going down, 3-going up further, 4-going up, 5-center
    select case (direction)

    case ("down")
        ! ==================================================
        !
        ! Going down
        !
        ! ==================================================
        ! To avoid current crossover

        ! Exception
        if(mesh.node(node1).dn == -1 .or. mesh.node(node2).up == -1) then
            b_fail = .true.
            return
        end if

        node1  = mesh.node(mesh.node(node1).dn).dn
        node2  = mesh.node(mesh.node(node2).up).up

        ! Exception
        if(node1 == -1 .or. node2 == -1) then
            b_fail = .true.
            return
        end if

        id_bp1 = mesh.node(node1).bp
        id_bp2 = mesh.node(node2).bp

        ! Find other downward crossovers
        n_skip = 0
        do
            ! Exception module for crossovers near vertex
            if( id_bp1 < min_bp1 + para_gap_xover_bound_scaf .or. &
                id_bp1 > max_bp1 - para_gap_xover_bound_scaf .or. &
                id_bp2 < min_bp2 + para_gap_xover_bound_scaf .or. &
                id_bp2 > max_bp2 - para_gap_xover_bound_scaf ) then
                b_fail = .true.
                exit
            end if

            ! Find crossover
            if(Section_Connection_Scaf(geom, sec1, sec2, id_bp1) == .true.) then
                if(n_skip == n_cross) exit
                n_skip = n_skip + 1
            end if

            ! Update node and base pair ID
            node1  = mesh.node(node1).dn
            node2  = mesh.node(node2).up
            id_bp1 = mesh.node(node1).bp
            id_bp2 = mesh.node(node2).bp
        end do
    case ("up")
        ! ==================================================
        !
        ! Going up
        !
        ! ==================================================
        ! To avoid current crossover

        ! Exception
        if(mesh.node(node1).up == -1 .or. mesh.node(node2).dn == -1) then
            b_fail = .true.
            return
        end if

        node1  = mesh.node(mesh.node(node1).up).up
        node2  = mesh.node(mesh.node(node2).dn).dn

        ! Exception
        if(node1 == -1 .or. node2 == -1) then
            b_fail = .true.
            return
        end if

        id_bp1 = mesh.node(node1).bp
        id_bp2 = mesh.node(node2).bp

        ! Find another upper crossover
        n_skip = 0
        do
            ! Exception module for crossovers near vertex
            if( id_bp1 < min_bp1 + para_gap_xover_bound_scaf .or. &
                id_bp1 > max_bp1 - para_gap_xover_bound_scaf .or. &
                id_bp2 < min_bp2 + para_gap_xover_bound_scaf .or. &
                id_bp2 > max_bp2 - para_gap_xover_bound_scaf ) then
                b_fail = .true.
                exit
            end if

            ! Find another crossover
            if(Section_Connection_Scaf(geom, sec1, sec2, id_bp1) == .true.) then
                if(n_skip == n_cross) exit
                n_skip = n_skip + 1
            end if

            ! Update node and base pair ID
            node1  = mesh.node(node1).up
            node2  = mesh.node(node2).dn
            id_bp1 = mesh.node(node1).bp
            id_bp2 = mesh.node(node2).bp
        end do
    case ("center")
        ! ==================================================
        !
        ! Reamin at center
        !
        ! ==================================================
        b_fail = .false.
    end select

    ! Update return value
    if(b_fail == .false.) then
        cur = node1
        com = node2
    end if
end function Route_Split_Centered_Scaf_Xover

! ---------------------------------------------------------------------------------------

! Write centered scaffold crossovers
! Last updated on Friday 15 July 2016 by Hyungmin
subroutine Route_Write_Centered_Scaf_Xover(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    ! Section type data
    type :: SecType
        integer :: start_bp, end_bp
        integer,   allocatable :: bp(:)
        integer,   allocatable :: node(:)
        integer,   allocatable :: xover(:)
    end type SecType

    ! Initial line data
    type :: EdgeType
        integer :: n_sec
        type(SecType), allocatable :: sec(:)
    end type EdgeType

    type(EdgeType), allocatable :: edge(:)
    integer,        allocatable :: conn_scaf(:,:)

    integer :: i, j, k, m, iniL, node, sec, bp, max_bp, min_bp
    integer :: n_conn_scaf, n_edge
    logical :: b_conn_scaf
    character(200) :: path

    if(para_write_610 == .false.) return

    path = trim(prob.path_work1)
    open(unit=610, file=trim(path)//"scaf_xover.txt", form="formatted")

    ! Find maximum and minimum base pair ID
    max_bp = mesh.node(1).bp
    min_bp = mesh.node(1).bp
    do i = 1, mesh.n_node
        if(mesh.node(i).bp > max_bp) max_bp = mesh.node(i).bp
        if(mesh.node(i).bp < min_bp) min_bp = mesh.node(i).bp
    end do
    min_bp = min_bp - 2
    max_bp = max_bp + 2

    ! Allocate and initialize edge data
    n_edge = geom.n_iniL
    allocate(edge(n_edge))
    do i = 1, n_edge
        edge(i).n_sec = geom.n_sec
        allocate(edge(i).sec(edge(i).n_sec))

        do j = 1, edge(i).n_sec
            allocate(edge(i).sec(j).bp   (min_bp:max_bp))
            allocate(edge(i).sec(j).node (min_bp:max_bp))
            allocate(edge(i).sec(j).xover(min_bp:max_bp))

            do k = min_bp, max_bp
                edge(i).sec(j).bp(k)    = k
                edge(i).sec(j).node(k)  = -1
                edge(i).sec(j).xover(k) = -1
            end do
        end do
    end do

    ! Allocate conn data
    allocate(conn_scaf(dna.n_xover_scaf, 2))

    ! Build node information based on initial edges
    do i = 1, mesh.n_node

        ! Find starting node depending on cross-section ID
        node = mesh.node(i).id
        sec  = mesh.node(i).sec
        iniL = mesh.node(node).iniL

        ! It depends on cross-section ID and direction
        if(mesh.node(i).dn == -1 .and. mod(mesh.node(i).sec, 2) == 0) then

            ! If section ID is even and negative z-direction
            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).up
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp
        else if(mesh.node(i).up == -1 .and. mod(mesh.node(i).sec, 2) == 1) then

            ! If section ID is odd and positive z-direction
            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).dn
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp
        end if
    end do

    ! Build additional data fields
    do i = 1, dna.n_base_scaf
        node = dna.base_scaf(i).node

        ! Skip the loop if Tn
        if(node == -1) cycle

        sec  = mesh.node(node).sec
        iniL = mesh.node(node).iniL
        bp   = mesh.node(node).bp

        ! For scaffold strand
        if(dna.base_scaf(i).xover /= -1) then
            edge(iniL).sec(sec+1).xover(bp) = mesh.node(dna.base_scaf(dna.base_scaf(i).xover).node).sec
        end if
    end do

    ! Print information based on edges
    do i = 1, n_edge

        ! Print base pair ID
        if(i == 1) then
            write(610, "(a)")
            write(610, "(a)"), "+------------------------------------------------------------------+"
            write(610, "(a)"), "|                                                                  |"
            write(610, "(a)"), "|              Possible centered scaffold crossovers               |"
            write(610, "(a)"), "|                                                                  |"
            write(610, "(a)"), "+------------------------------------------------------------------+"
            write(610, "(a)")
            write(610, "(a)"), "1. [-], [num] : Nucleotide"
            write(610, "(a)"), "2. [|]        : Crossover"
            write(610, "(a)"), "4. [ a: b], c : From starting base ID(a) to ending base ID(b), total # of bases (c)"
            write(610, "(a)"), "5. [5'], [3'] : Strand direction"
            write(610, "(a)"); write(610, "(a)")
        end if

        ! Scaffold strand
        n_conn_scaf = 0
        do j = 1, edge(i).n_sec

            ! Print edge and point information, Edge 1 - (point 1 -> point 2)
            if(j == 1) then
                write(610, "(a)") " ==================================================="
                call Space(610, 10)
                write(610, "(a)"), " [Edge "//trim(adjustl(Int2Str(i)))//&
                    " : point "//trim(adjustl(Int2Str(geom.iniL(i).poi(1))))//&
                    " -> point "//trim(adjustl(Int2Str(geom.iniL(i).poi(2))))//"]"
                write(610, "(a)") " ==================================================="
                write(610, "(a)")
            end if

            ! Print section ID and base length, sec xx -->> [xxx:xxx], xxx : - 26 length
            write(610, "(a5$)"), " sec "
            write(610, "(a2$)"), trim(adjustl(Int2Str(j-1)))
            write(610, "(a4$)"), " - ["
            write(610, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).start_bp + para_start_bp_ID - 1)))
            write(610, "(a1$)"), ":"
            write(610, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp + para_start_bp_ID - 1)))
            write(610, "(a3$)"), "], "
            write(610, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp - edge(i).sec(j).start_bp + 1)))
            write(610, "(a2$)"), " :"
            call Space(610, 4)

            ! Graphic representation, bases with crossover and nick
            do k = min_bp, max_bp

                if(edge(i).sec(j).node(k) == -1) then
                    if(edge(i).sec(j).start_bp - 2 == k) then
                        if(mod(j-1, 2) == 0) write(610, "(a$)"), "5"
                        if(mod(j-1, 2) == 1) write(610, "(a$)"), "3"
                    else if(edge(i).sec(j).end_bp + 1 == k) then
                        if(mod(j-1, 2) == 0) write(610, "(a$)"), "3"
                        if(mod(j-1, 2) == 1) write(610, "(a$)"), "5"
                    else if(edge(i).sec(j).start_bp - 1 == k .or. edge(i).sec(j).end_bp + 2 == k) then
                        write(610, "(a$)"), "'"
                    else
                        write(610, "(a$)"), " "
                    end if
                else
                    if(edge(i).sec(j).xover(k) /= -1) then
                        !write(610, "(a$)"), "+"

                        if(edge(i).sec(j).xover(k) >= 10) then
                            ! ANSI character code, 65 -> A
                            write(610, "(a$)"), achar(55+edge(i).sec(j).xover(k))
                        else
                            write(610, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).xover(k))))
                        end if

                        n_conn_scaf = n_conn_scaf + 1
                        conn_scaf(n_conn_scaf, 1) = edge(i).sec(j).xover(k)
                        conn_scaf(n_conn_scaf, 2) = k
                    else
                        write(610, "(a$)"), "-"
                    end if
                end if
            end do
            write(610, "(a)")

            ! Draw crossovers
            call Space(610, 30)
            do k = min_bp, max_bp
                if(edge(i).sec(j).node(k) == -1) then
                    write(610, "(a$)"), " "
                else
                    b_conn_scaf = .false.
                    do m = 1, n_conn_scaf
                        if(conn_scaf(m, 1) >= j .and. conn_scaf(m, 2) == k) then
                            b_conn_scaf = .true.
                            exit
                        end if
                    end do

                    if(b_conn_scaf == .true.) then
                        write(610, "(a$)"), "|"
                    else
                        write(610, "(a$)"), " "
                    end if
                end if
            end do
            write(610, "(a)")
        end do
        write(610, "(a)")
    end do

    ! Deallocate memory
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            deallocate(edge(i).sec(j).bp   )
            deallocate(edge(i).sec(j).node )
            deallocate(edge(i).sec(j).xover)
        end do
        deallocate(edge(i).sec)
    end do
    deallocate(edge)
    deallocate(conn_scaf)

    ! Close file
    close(unit=610)
end subroutine Route_Write_Centered_Scaf_Xover

! ---------------------------------------------------------------------------------------

! Modify scaffold crossovers to avoid duplicated indication
! Last updated on Thuesday 12 June 2016 by Hyungmin
subroutine Route_Modify_Scaf_Xover(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type :: CheckType
        integer :: strd(2)
        integer :: base1(2)
        integer :: base2(2)
    end type CheckType

    type(CheckType), allocatable :: check(:)

    integer :: i, j, count, strd_cur, strd_com, n_xover1, n_xover2
    integer :: base, up_base, dn_base, xover, up_xover, dn_xover
    logical :: b_delete, b_skip

    ! Print information
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.6.  Modify centered crossover to avoid duplicate one"
    end do

    ! 1. Delete crossovers indicating the same strand number
    n_xover1 = dna.n_xover_scaf
    do i = 1, dna.n_base_scaf

        ! Check crossover
        xover = dna.base_scaf(i).xover

        ! If there if crossover
        if(xover /= -1) then

            ! Delete crossover if the strand IDs are identical
            if(dna.base_scaf(i).strand == dna.base_scaf(xover).strand) then
                dna.n_xover_scaf           = dna.n_xover_scaf - 1
                dna.base_scaf(i).xover     = -1
                dna.base_scaf(xover).xover = -1
            end if
        end if
    end do

    ! 2. Delete crossovers indicating the same strand number
    n_xover2 = dna.n_xover_scaf
    allocate(check(dna.n_xover_scaf))

    count = 0
    do i = 1, dna.n_base_scaf

        base  = dna.base_scaf(i).id
        xover = dna.base_scaf(base).xover

        up_base  = dna.base_scaf(base).up
        dn_base  = dna.base_scaf(base).dn
        up_xover = dna.base_scaf(up_base).xover
        dn_xover = dna.base_scaf(dn_base).xover

        ! Check double crossovers
        if(xover /= -1 .and. (up_xover /= -1 .or. dn_xover /= -1)) then

            ! Find strand ID
            strd_cur = dna.base_scaf(base).strand
            strd_com = dna.base_scaf(xover).strand

            ! Check whether the crossover exists or not in check array
            b_delete = .false.
            do j = 1, count
                if( (check(j).strd(1) == strd_cur .and. check(j).strd(2) == strd_com) .or. &
                    (check(j).strd(2) == strd_cur .and. check(j).strd(1) == strd_com) ) then
                    b_delete = .true.
                    exit
                end if
            end do

            if(b_delete == .false.) then

                ! Add new crossover into check array
                count = count + 1
                check(count).base1(1) = base
                check(count).base1(2) = xover
                check(count).strd(1)  = dna.base_scaf(base).strand
                check(count).strd(2)  = dna.base_scaf(xover).strand

                if(up_xover /= -1) then
                    check(count).base2(1) = up_base
                    check(count).base2(2) = up_xover
                else if(dn_xover /= -1) then
                    check(count).base2(1) = dn_base
                    check(count).base2(2) = dn_xover
                end if

            else if(b_delete == .true.) then

                ! Check double crossover to be deleted
                b_skip = .false.
                do j = 1, count
                    if( (base == check(j).base1(1) .and. xover == check(j).base1(2)) .or. &
                        (base == check(j).base1(2) .and. xover == check(j).base1(1)) .or. &
                        (base == check(j).base2(1) .and. xover == check(j).base2(2)) .or. &
                        (base == check(j).base2(2) .and. xover == check(j).base2(1)) &
                        ) then
                        b_skip = .true.
                        exit
                    end if
                end do

                ! Delete double crossovers
                if(b_skip == .false.) then
                    dna.n_xover_scaf = dna.n_xover_scaf - 2
                    dna.base_scaf(base).xover  = -1
                    dna.base_scaf(xover).xover = -1

                    if(up_xover /= -1) then
                        dna.base_scaf(up_base).xover  = -1
                        dna.base_scaf(up_xover).xover = -1
                    else if(dn_xover /= -1) then
                        dna.base_scaf(dn_base).xover  = -1
                        dna.base_scaf(dn_xover).xover = -1
                    end if
                end if
            end if
        end if
    end do

    ! Deallocate memory
    deallocate(check)

    ! Print information
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The deleted self connecting crossovers      : "&
            //trim(adjustl(Int2Str(n_xover1 - n_xover2)))
        call Space(i, 11)
        write(i, "(a)"), "* The deleted duplicate indexing crossovers   : "&
            //trim(adjustl(Int2Str(n_xover2 - dna.n_xover_scaf)))
        call Space(i, 11)
        write(i, "(a)"), "* The number of scaffold crossovers           : "&
            //trim(adjustl(Int2Str(dna.n_xover_scaf)))
        write(i, "(a)")
    end do
end subroutine Route_Modify_Scaf_Xover

! ---------------------------------------------------------------------------------------

! Build scaffold DNA origami
! Last updated on Tuesday 10 Mar 2016 by Hyungmin
subroutine Route_Graph_Build_Origami(prob, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(Meshtype), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    double precision, allocatable :: pos_node(:,:)
    integer, allocatable :: adj(:,:), tail(:), head(:), cost(:)
    integer, allocatable :: idx(:,:), src(:),  dst(:),  tree(:)

    integer :: i, j, n_node, n_edge, length

    ! Set the number of node and edge that should be even
    n_node = dna.n_scaf
    n_edge = dna.n_xover_scaf / 2

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.7. Make scaffold origami"
        call Space(i, 11)
        write(i, "(a)"), "* Build dual-graph based on closed scaffolds"
        call Space(i, 15)
        write(i, "(a)"), "- The number of nodes in dual-graph       : "&
            //trim(adjustl(Int2Str(n_node)))
        call Space(i, 15)
        write(i, "(a)"), "- The number of edges in dual-graph       : "&
            //trim(adjustl(Int2Str(n_edge)))
        call Space(i, 11)
        write(i, "(a)"), "* Generate spanning trees"
    end do

    ! Allocate graph data
    call Route_Graph_Allocate_Data(pos_node, tail, head, cost, n_node, n_edge)

    ! Set node and edge for dual graph
    call Route_Graph_Set_Data(prob, mesh, dna, pos_node, tail, head, cost)

    ! Possible all spanning tree when # of edges is less than 12 for Prim or Kruskal
    if(para_all_spanning == "on" .and. n_edge <= 12) then

        ! Convert from edge list to adjacent matrix
        call SpanTree_List2Adj(adj, tail, head, cost, n_node, n_edge)

        ! Write Adjacent matrix
        !call Route_Graph_Write_Adjacent(prob, adj)

        ! Print the number of all spanning trees
        do i = 0, 11, 11
            call Space(i, 15)
            write(i, "(a$)"), "- The number of all spanning trees        : "
            write(i, "(i7)"), SpanTree_Count_Spanning_Trees(adj)
        end do

        ! Generate spanning trees
        call SpanTree_Generate_Spanning_Trees(adj, idx, src, dst)

        ! Print all spanning trees
        call SpanTree_Print_All_Trees(idx, src, dst)

        ! Write all spanning trees
        call Route_Graph_Chimera_All_Spanning_Tree(prob, pos_node, tail, head, idx, src, dst)
    end if

    dna.graph_node = n_node
    dna.graph_edge = n_edge

    ! Find minimal spanning tree
    if(para_method_MST == "prim") then

        ! Prim's algorithm (none, quick, shell sorting)
        call SpanTree_Prim_Algorithm_1(tail, head, cost, tree, n_node, n_edge, length, trim(para_method_sort))
        !call SpanTree_Prim_Algorithm_2(tail, head, cost, tree, n_node, n_edge, length, trim(para_method_sort))

    else if(para_method_MST == "kruskal") then

        ! Kruskal's algorithm (none, quick, shell sorting)
        call SpanTree_Kruskal_Algorithm(tail, head, cost, tree, n_node, n_edge, length, trim(para_method_sort))

    end if

    ! Print edges with minimal spanning tree
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* Find the minimal spanning tree (all edge weight - 1)"
        call Space(i, 15)
        write(i, "(a)"), "- The toal cost for this path             : "//trim(adjustl(Int2Str(length)))
        call Space(i, 15)
        write(i, "(a)"), "- The minimal spanning tree"
        call Space(i, 15)
        write(i, "(a)"), "  Arc    Tail node Head node   Cost"
        call Space(i, 15)
        write(i, "(a)"), "  ---    --------- ---------   ----"
        do j = 1, size(tree)
            call Space(i, 9)
            write(i, "(4i10)"), tree(j), tail(tree(j)), head(tree(j)), cost(tree(j))
            if(cost(tree(j)) == 3) stop
        end do
        write(i, "(a)")
    end do

    ! Write adjacent list for dual graph
    if(para_adjacent_list == "on") call Route_Graph_Write_List(prob, tail, head, tree)

    ! Write Chimera for the spanning tree
    call Route_Graph_Chimera_Spanning_Tree(prob, pos_node, tail, head, tree)

    ! Delete non-spanning scaffold crossover
    call Route_Graph_Delete_Scaf_Xover(dna, tail, head, tree)

    ! Deallocate graph memory
    if(allocated(pos_node)) deallocate(pos_node)
    if(allocated(adj) )     deallocate(adj)
    if(allocated(tail))     deallocate(tail)
    if(allocated(head))     deallocate(head)
    if(allocated(cost))     deallocate(cost)
    if(allocated(idx) )     deallocate(idx)
    if(allocated(src) )     deallocate(src)
    if(allocated(dst) )     deallocate(dst)
    if(allocated(tree))     deallocate(tree)
end subroutine Route_Graph_Build_Origami

! ---------------------------------------------------------------------------------------

! Allocate graph data
! Last updated on Tuesday 10 May 2016 by Hyungmin
subroutine Route_Graph_Allocate_Data(pos_node, tail, head, cost, n_node, n_edge)
    double precision, allocatable, intent(inout) :: pos_node(:,:)
    integer,          allocatable, intent(inout) :: tail(:)
    integer,          allocatable, intent(inout) :: head(:)
    integer,          allocatable, intent(inout) :: cost(:)
    integer, intent(in) :: n_node
    integer, intent(in) :: n_edge

    ! Allocate graph memory
    allocate(pos_node(n_node, 3))
    allocate(tail(n_edge))
    allocate(head(n_edge))
    allocate(cost(n_edge))

    ! Initialize node and edge graph data
    pos_node(1:n_node, :) = 0.0d0
    tail(1:n_edge)        = 0
    head(1:n_edge)        = 0
    cost(1:n_edge)        = 0
end subroutine Route_Graph_Allocate_Data

! ---------------------------------------------------------------------------------------
    
! Set node and edge for dual graph
! Last updated on Thursday 10 November 2016 by Hyungmin
subroutine Route_Graph_Set_Data(prob, mesh, dna, pos_node, tail, head, cost)
    type(ProbType),   intent(in)    :: prob
    type(MeshType),   intent(in)    :: mesh
    type(DNAType),    intent(in)    :: dna
    double precision, intent(inout) :: pos_node(:,:)
    integer,          intent(inout) :: tail(:)
    integer,          intent(inout) :: head(:)
    integer,          intent(inout) :: cost(:)

    integer, allocatable :: n_base(:)

    !integer :: pre_con1_m(2), pre_con2_m(2), pre_con1_b(2), pre_con2_b(2), pre_conm_m(2), pre_conm_b(2)
    integer :: con_pri1(2), con_pri2(2), con_span(2)
    integer :: i, j, strand, xover, strand_xover, n_node, n_edge
    integer :: node1, sec1, node2, sec2
    logical :: b_add

    ! count the number of nodes and edges for dual graph
    n_node = ubound(pos_node,1)
    n_edge = size(tail)

    ! Array for counting base in terms of closed strand loops
    allocate(n_base(n_node))
    n_base(1:n_node) = 0

    ! Build node and edge data
    n_edge = 0
    do i = 1, dna.n_base_scaf
        strand = dna.base_scaf(i).strand
        xover  = dna.base_scaf(i).xover

        n_base(strand)     = n_base(strand) + 1
        pos_node(strand,:) = pos_node(strand,:) + dna.base_scaf(i).pos(:)

        ! If there is crossover
        if(xover /= -1) then

            ! Strand ID of the base's crossover
            strand_xover = dna.base_scaf(xover).strand

            ! Check previous queue
            b_add = .true.
            do j = 1, n_edge
                if( (tail(j) == strand .and. head(j) == strand_xover) .or. &
                    (head(j) == strand .and. tail(j) == strand_xover) ) then
                    b_add = .false.
                    exit
                end if
            end do

            ! If there is no dupilcated edge
            if(b_add == .true.) then

                ! The number of edges
                n_edge = n_edge + 1

                ! Make edge list
                tail(n_edge) = strand
                head(n_edge) = strand_xover
                cost(n_edge) = 1

                ! Set weight factor depending on cross-sectional ID
                node1 = dna.base_scaf(i).id
                node2 = dna.base_scaf(i).xover
                sec1  = mesh.node(node1).sec
                sec2  = mesh.node(node2).sec

                ! Set priority
                if(prob.sel_sec == 1) then

                    ! 6HB bottom connection
                    con_pri1(:) = [0, 1]    ! Connection priority 1
                    con_pri2(:) = [0, 1]    ! Connection priority 2
                    con_span(:) = [0, 1]    ! Connection for spanning tree
                else if(prob.sel_sec == 2) then

                    ! 6HB bottom connection
                    con_pri1(:) = [0, 5]    ! Connection priority 1
                    con_pri2(:) = [3, 4]    ! Connection priority 2
                    con_span(:) = [0, 1]    ! Connection for spanning tree
                
                else if(prob.sel_sec == 3) then

                    ! Middle reference section, pattern #13
                    con_pri1(:) = [2, 3]    ! Connection priority 1
                    con_pri2(:) = [3, 4]    ! Connection priority 2
                    con_span(:) = [0, 5]    ! Connection for spanning tree

                    ! Middle reference section, pattern #23
                    !con_pri1(:) = [0, 5]    ! Connection priority 1
                    !con_pri2(:) = [2, 3]    ! Connection priority 2
                    !con_span(:) = [3, 4]    ! Connection for spanning tree
                else if(prob.sel_sec == 4) then

                    ! 6HB middle connection
                    con_pri1(:) = [2, 3]    ! Connection priority 1
                    con_pri2(:) = [3, 4]    ! Connection priority 2
                    con_span(:) = [0, 5]    ! Connection for spanning tree

                else if(prob.sel_sec == 5) then

                    ! 6HB middle connection
                    con_pri1(:) = [0, 5]    ! Connection priority 1
                    con_pri2(:) = [4, 5]    ! Connection priority 2
                    con_span(:) = [1, 2]    ! Connection for spanning tree
                end if

                if(para_weight_edge == "on") then
                    if( (sec1 == con_pri1(1) .and. sec2 == con_pri1(2)) .or. (sec1 == con_pri1(2) .and. sec2 == con_pri1(1)) ) then
                        cost(n_edge) = 1
                    else if( (sec1 == con_pri2(1) .and. sec2 == con_pri2(2)) .or. (sec1 == con_pri2(2) .and. sec2 == con_pri2(1)) ) then
                        cost(n_edge) = 1
                    else if( (sec1 == con_span(1) .and. sec2 == con_span(2)) .or. (sec1 == con_span(2) .and. sec2 == con_span(1)) ) then
                        cost(n_edge) = 2
                    else
                        cost(n_edge) = 3
                    end if
                end if
            end if
        end if
    end do

    ! Find center position of close loop strands (= node)
    do i = 1, n_node
        pos_node(i,:) = pos_node(i,:) / n_base(i)
    end do

    ! Deallocate memory
    deallocate(n_base)
end subroutine Route_Graph_Set_Data

! ---------------------------------------------------------------------------------------

! Write Adjacent matrix
! Last updated on Wednesday 4 May 2016 by Hyungmin
subroutine Route_Graph_Write_Adjacent(prob, adj)
    type(ProbType),   intent(in) :: prob
    integer,          intent(in) :: adj(:,:)

    character(200) :: path
    integer :: i, j, dim

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=603, file=trim(path)//"_adjacent.dl", form="formatted")

    ! Find the number of array size
    dim = ubound(adj, 1)

    ! Write UCINET DL Format
    ! UCINET DL format is the most common file format used by UCINET package
    write(603, "(a)"), "dl n="//trim(adjustl(Int2Str(dim)))
    write(603, "(a)"), "format = fullmatrix"

    ! Print label as edge number
    write(603, "(a)"), "labels:"
    do i = 1, dim
        if(i == dim) then
            write(603, "(a)"), trim(adjustl(Int2Str(i)))
        else
            write(603, "(a$)"), trim(adjustl(Int2Str(i)))//","
        end if
    end do

    ! Write full matrix
    write(603, "(a)"), "data:"
    do i = 1, dim
        do j = 1, dim
            if(j /= dim) then
                write(603, "(i7$)"), adj(i, j)
            else
                write(603, "(i7)"), adj(i, j)
            end if
        end do
    end do

    ! Closed unit
    close(unit=603)
end subroutine Route_Graph_Write_Adjacent

! ---------------------------------------------------------------------------------------

! Write all spanning trees
! Last updated on Friday 22 July 2016 by Hyungmin
subroutine Route_Graph_Chimera_All_Spanning_Tree(prob, pos_node, tail, head, idx, src, dst)
    type(ProbType),   intent(in) :: prob
    double precision, intent(in) :: pos_node(:,:)
    integer,          intent(in) :: tail(:)
    integer,          intent(in) :: head(:)
    integer,          intent(in) :: idx(:,:)
    integer,          intent(in) :: src(:)
    integer,          intent(in) :: dst(:)

    double precision :: pos_1(3), pos_2(3), radius
    logical :: f_axis, b_branch, results
    integer :: i, j, k, n_tree, n_node, n_edge
    character(200) :: path

    f_axis = para_chimera_axis

    n_tree  = ubound(idx, 2)
    n_node  = size(pos_node, 1)
    n_edge  = size(tail)
    results = SYSTEMQQ("md "//trim(prob.path_work1)//"Spantree\")

    do k = 1, n_tree

        ! File open
        path = trim(prob.path_work1)//"Spantree\"//trim(prob.name_file)
        open(unit=604, file=trim(path)//"_graph"//trim(adjustl(Int2Str(k)))//".bild", form="formatted")

        ! Write node as sphere
        write(604, "(a)"), ".color tan"
        do i = 1, ubound(pos_node, 1)
            write(604, "(a$    )"), ".sphere "
            write(604, "(3f9.3$)"), pos_node(i, 1:3)
            write(604, "(1f9.3 )"), 0.4d0
        end do

        ! Write edge as cylinder
        do i = 1, n_edge
            pos_1(:) = pos_node(tail(i), 1:3)
            pos_2(:) = pos_node(head(i), 1:3)

            b_branch = .false.
            do j = 1, size(idx, 1)
                if( (src(idx(j,k)) == tail(i) .and. dst(idx(j,k)) == head(i)) .or. &
                    (src(idx(j,k)) == head(i) .and. dst(idx(j,k)) == tail(i)) ) then
                    b_branch = .true.
                    exit
                end if
            end do

            if(b_branch == .true.) then
                write(604, "(a)"), ".color red"
                radius = 0.15d0
            else
                write(604, "(a)"), ".color steel blue"
                radius = 0.15d0
            end if

            write(604, "(a$    )"), ".cylinder "
            write(604, "(3f9.3$)"), pos_1(1:3)
            write(604, "(3f9.3$)"), pos_2(1:3)
            write(604, "(1f9.3 )"), radius
        end do

        ! Write global axis
        if(f_axis == .true.) then
            write(604, "(a)"), ".translate 0.0 0.0 0.0"
            write(604, "(a)"), ".scale 0.5"
            write(604, "(a)"), ".color grey"
            write(604, "(a)"), ".sphere 0 0 0 0.5"      ! Center
            write(604, "(a)"), ".color red"             ! x-axis
            write(604, "(a)"), ".arrow 0 0 0 4 0 0 "
            write(604, "(a)"), ".color blue"            ! y-axis
            write(604, "(a)"), ".arrow 0 0 0 0 4 0 "
            write(604, "(a)"), ".color yellow"          ! z-axis
            write(604, "(a)"), ".arrow 0 0 0 0 0 4 "
        end if
        close(unit=604)

        ! ---------------------------------------------
        !
        ! Write the file for Tecplot
        !
        ! ---------------------------------------------
        if(para_output_Tecplot == "on") then

            path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
            open(unit=604, file=trim(path)//"_graph"//trim(adjustl(Int2Str(k)))//".dat", form="formatted")

            write(604, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
            write(604, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
            write(604, "(a$)"), 'ZONE F = FEPOINT'
            write(604, "(a$)"), ', N='//trim(adjustl(Int2Str(n_node)))
            write(604, "(a$)"), ', E='//trim(adjustl(Int2Str(n_edge - size(idx, 1))))
            write(604, "(a )"), ', ET=LINESEG'

            ! Write nodes
            do i = 1, n_node
                write(604, "(3f9.3$)"), pos_node(i, 1:3)
                write(604, "(1f9.3 )"), 1.0d0
            end do

            ! Write edges
            do i = 1, n_edge

                b_branch = .false.
                do j = 1, size(idx, 1)
                    if( (src(idx(j,k)) == tail(i) .and. dst(idx(j,k)) == head(i)) .or. &
                        (src(idx(j,k)) == head(i) .and. dst(idx(j,k)) == tail(i)) ) then
                    b_branch = .true.
                    exit
                    end if
                end do

                if(b_branch == .false.) then
                    write(604, "(1i7$)"), tail(i)
                    write(604, "(1i7 )"), head(i)
                end if
            end do

            write(604, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
            write(604, "(a$)"), 'ZONE F = FEPOINT'
            write(604, "(a$)"), ', N='//trim(adjustl(Int2Str(n_node)))
            write(604, "(a$)"), ', E='//trim(adjustl(Int2Str(size(idx, 1))))
            write(604, "(a )"), ', ET=LINESEG'

            ! Write nodes
            do i = 1, n_node
                write(604, "(3f9.3$)"), pos_node(i, 1:3)
                write(604, "(1f9.3 )"), 1.0d0
            end do

            ! Write edges
            do i = 1, n_edge

                b_branch = .false.
                do j = 1, size(idx, 1)
                    if( (src(idx(j,k)) == tail(i) .and. dst(idx(j,k)) == head(i)) .or. &
                        (src(idx(j,k)) == head(i) .and. dst(idx(j,k)) == tail(i)) ) then
                    b_branch = .true.
                    exit
                    end if
                end do

                if(b_branch == .true.) then
                    write(604, "(1i7$)"), tail(i)
                    write(604, "(1i7 )"), head(i)
                end if
            end do
            close(unit=604)
        end if
    end do
end subroutine Route_Graph_Chimera_All_Spanning_Tree

! ---------------------------------------------------------------------------------------

! Write edge list
! Last updated on Tuesday 10 May 2016 by Hyungmin
subroutine Route_Graph_Write_List(prob, tail, head, tree)
    type(ProbType), intent(in) :: prob
    integer,        intent(in) :: tail(:)
    integer,        intent(in) :: head(:)
    integer,        intent(in) :: tree(:)

    double precision :: weight
    integer :: i, j, n_node
    logical :: b_branch
    character(200) :: path, str

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=605, file=trim(path)//"_list.dl", form="formatted")

    ! Find the number of nodes
    n_node = size(tree) + 1

    ! Write UCINET DL Format
    ! UCINET DL format is the most common file format used by UCINET package
    write(605, "(a)"), "dl n="//trim(adjustl(Int2Str(n_node)))
    write(605, "(a)"), "format = edgelist1"

    ! Print label as edge number
    write(605, "(a)"), "labels:"
    do i = 1, n_node
        if(i == n_node) then
            write(605, "(a)"), trim(adjustl(Int2Str(i)))
        else
            write(605, "(a$)"), trim(adjustl(Int2Str(i)))//","
        end if
    end do

    ! Write edge list
    write(605, "(a)"), "data:"
    do i = 1, size(tail)

        weight   = 1.0d0
        b_branch = .false.
        do j = 1, size(tree)
            if( (tail(tree(j)) == tail(i) .and. head(tree(j)) == head(i)) .or. &
                (tail(tree(j)) == head(i) .and. head(tree(j)) == tail(i)) ) then
                b_branch = .true.
                exit
            end if
        end do

        if(b_branch == .true.) weight = 2.0d0

        write(605, "(2i6, f8.3)"), tail(i), head(i), weight
        write(605, "(2i6, f8.3)"), head(i), tail(i), weight
    end do

    ! Closed unit
    close(unit=605)
end subroutine Route_Graph_Write_List

! ---------------------------------------------------------------------------------------

! Write Chimera for the spanning tree
! Last updated on Saturday 16 July 2016 by Hyungmin
subroutine Route_Graph_Chimera_Spanning_Tree(prob, pos_node, tail, head, tree)
    type(ProbType),   intent(in) :: prob
    double precision, intent(in) :: pos_node(:,:)
    integer,          intent(in) :: tail(:)
    integer,          intent(in) :: head(:)
    integer,          intent(in) :: tree(:)

    double precision :: pos_1(3), pos_2(3), radius
    logical :: f_axis, b_branch
    integer :: i, j, n_node, n_edge
    character(200) :: path

    if(para_write_606 == .false.) return

    f_axis = para_chimera_axis

    n_node = ubound(pos_node, 1)
    n_edge = size(tail)

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=606, file=trim(path)//"_spantree.bild", form="formatted")

    ! Write node as sphere
    write(606, "(a)"), ".color steel blue"
    do i = 1, n_node
        write(606, "(a$    )"), ".sphere "
        write(606, "(3f9.3$)"), pos_node(i, 1:3)
        write(606, "(1f9.3 )"), 0.35d0
    end do

    ! Write edge as cylinder
    do i = 1, n_edge
        pos_1(:) = pos_node(tail(i), 1:3)
        pos_2(:) = pos_node(head(i), 1:3)

        b_branch = .false.
        do j = 1, size(tree)
            if( (tail(tree(j)) == tail(i) .and. head(tree(j)) == head(i)) .or. &
                (tail(tree(j)) == head(i) .and. head(tree(j)) == tail(i)) ) then
                b_branch = .true.
                exit
            end if
        end do

        if(b_branch == .true.) then
            write(606, "(a)"), ".color red"
            radius = 0.15d0
        else
            write(606, "(a)"), ".color tan"
            radius = 0.1d0
        end if

        write(606, "(a$    )"), ".cylinder "
        write(606, "(3f9.3$)"), pos_1(1:3)
        write(606, "(3f9.3$)"), pos_2(1:3)
        write(606, "(1f9.3 )"), radius
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(606, "(a)"), ".translate 0.0 0.0 0.0"
        write(606, "(a)"), ".scale 0.5"
        write(606, "(a)"), ".color grey"
        write(606, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(606, "(a)"), ".color red"             ! x-axis
        write(606, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(606, "(a)"), ".color blue"            ! y-axis
        write(606, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(606, "(a)"), ".color yellow"          ! z-axis
        write(606, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if

    close(unit=606)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=606, file=trim(path)//"_spantree.dat", form="formatted")

    write(606, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(606, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(606, "(a$)"), 'ZONE F = FEPOINT'
    write(606, "(a$)"), ', N='//trim(adjustl(Int2Str(n_node)))
    write(606, "(a$)"), ', E='//trim(adjustl(Int2Str(n_edge - size(tree))))
    write(606, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, n_node
        write(606, "(3f9.3$)"), pos_node(i, 1:3)
        write(606, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges - no spanning branch
    do i = 1, n_edge

        b_branch = .false.
        do j = 1, size(tree)
            if( (tail(tree(j)) == tail(i) .and. head(tree(j)) == head(i)) .or. &
                (tail(tree(j)) == head(i) .and. head(tree(j)) == tail(i)) ) then
                b_branch = .true.
                exit
            end if
        end do

        if(b_branch == .false.) then
        write(606, "(1i7$)"), tail(i)
        write(606, "(1i7 )"), head(i)
        end if
    end do

    write(606, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(606, "(a$)"), 'ZONE F = FEPOINT'
    write(606, "(a$)"), ', N='//trim(adjustl(Int2Str(n_node)))
    write(606, "(a$)"), ', E='//trim(adjustl(Int2Str(size(tree))))
    write(606, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, n_node
        write(606, "(3f9.3$)"), pos_node(i, 1:3)
        write(606, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges - spanning branch
    do i = 1, n_edge

        b_branch = .false.
        do j = 1, size(tree)
            if( (tail(tree(j)) == tail(i) .and. head(tree(j)) == head(i)) .or. &
                (tail(tree(j)) == head(i) .and. head(tree(j)) == tail(i)) ) then
                b_branch = .true.
                exit
            end if
        end do

        if(b_branch == .true.) then
        write(606, "(1i7$)"), tail(i)
        write(606, "(1i7 )"), head(i)
        end if
    end do

    close(unit = 606)
end subroutine Route_Graph_Chimera_Spanning_Tree

! ---------------------------------------------------------------------------------------

! Delete non-spanning scaffold crossover
! Last updated on Tuesday 10 Mar 2016 by Hyungmin
subroutine Route_Graph_Delete_Scaf_Xover(dna, tail, head, tree)
    type(DNAType), intent(inout) :: dna
    integer,       intent(in)    :: tail(:)
    integer,       intent(in)    :: head(:)
    integer,       intent(in)    :: tree(:)

    integer :: i, j, n_node, strand, xover, strand_xover, up, up_xover
    logical :: b_spanning

    n_node = size(tree) + 1

    ! Build node and edge data
    do i = 1, dna.n_base_scaf
        xover = dna.base_scaf(i).xover

        ! If there is crossover
        if(xover /= -1) then

            ! Strand ID of the base's crossover
            strand       = dna.base_scaf(i).strand
            strand_xover = dna.base_scaf(xover).strand

            ! Find spanning or non-spanning branch
            b_spanning = .false.
            do j = 1, size(tree)
                if( (tail(tree(j)) == strand       .and. head(tree(j)) == strand_xover) .or. &
                    (tail(tree(j)) == strand_xover .and. head(tree(j)) == strand)     ) then
                    b_spanning = .true.
                    exit
                end if
            end do

            ! Delete or reconstruct crossover
            if(b_spanning == .false.) then

                ! If non-spanning branch, delete
                dna.base_scaf(i).xover     = -1
                dna.base_scaf(xover).xover = -1
                dna.n_xover_scaf           = dna.n_xover_scaf - 1
            else

                ! If spanning branch, construct connections
                up       = dna.base_scaf(i).up
                up_xover = dna.base_scaf(up).xover

                ! If there is crossover
                if(up_xover /= -1) then

                    ! Crossover strand should be identical
                    if(dna.base_scaf(xover).strand == dna.base_scaf(up_xover).strand) then
                        dna.base_scaf(i).up    = xover
                        dna.base_scaf(up).dn = up_xover
                    end if
                end if
            end if
        end if
    end do
end subroutine Route_Graph_Delete_Scaf_Xover

! ---------------------------------------------------------------------------------------

! Make scaffold strand origami
! Last updated on Thursday 24 Mar 2016 by Hyungmin
subroutine Route_Make_Scaf_Origami(prob, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer, allocatable, dimension(:) :: id_entry
    logical, allocatable, dimension(:) :: b_visit

    integer :: i, j, count, n_xover, base
    integer :: strand, up, start, xover

    ! Allocate and initialize array
    allocate(b_visit(dna.n_scaf))
    allocate(id_entry(dna.n_scaf))

    do i = 1, dna.n_scaf
        b_visit(i)  = .false.
        id_entry(i) = -1
    end do

    ! Make discontinuous scaffold continuous
    ! Starting base ID is 1
    base   = dna.base_scaf(1).id
    strand = dna.base_scaf(base).strand

    b_visit(strand)  = .true.   ! Check visiting
    id_entry(strand) = base

    ! Loop for all bases in scaffold strands
    n_xover = 0
    do i = 2, dna.n_base_scaf

        xover = dna.base_scaf(base).xover
        start = id_entry(strand)
        up    = dna.base_scaf(base).up

        if(xover == -1) then
            ! Go to upper base
            base   = dna.base_scaf(base).up
            strand = dna.base_scaf(base).strand

        else if(xover /= -1 .and. b_visit(dna.base_scaf(xover).strand) == .false.) then
            ! --------------------------------------------------
            ! If there is crossovers and never visted before
            ! going in, set connectivity of two bases
            ! --------------------------------------------------
            dna.base_scaf(base).up  = xover
            dna.base_scaf(xover).dn = base

            ! Go into new strand
            base   = dna.base_scaf(base).xover
            strand = dna.base_scaf(base).strand

            ! Set entry base corresponding to strand
            id_entry(strand) = base

        else if(xover /= -1 .and. start == up) then
            ! --------------------------------------------------
            ! If all visited in the certain scaffold
            ! going out, set connectivity of two bases
            ! --------------------------------------------------
            dna.base_scaf(base).up  = xover
            dna.base_scaf(xover).dn = base

            ! Go out to the previous strand
            base   = dna.base_scaf(base).xover
            strand = dna.base_scaf(base).strand

        else
            ! Delete possible crossovers
            if(xover /= -1 .and. &
                xover /= dna.base_scaf(base).up .and. &
                xover /= dna.base_scaf(base).dn) then

                dna.base_scaf(xover).xover = -1
                dna.base_scaf(base).xover  = -1
                n_xover = n_xover + 1
            end if

            ! Go to upper base
            base   = dna.base_scaf(base).up
            strand = dna.base_scaf(base).strand
        end if

        b_visit(strand) = .true.
    end do

    ! Update the number of crossovers in scaffold strand
    dna.n_xover_scaf = dna.n_xover_scaf - n_xover

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.12. Continuous scaffold strand"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of deleted possible crossovers   : "
        write(i, "(i7)"), n_xover
        call Space(i, 11)
        write(i, "(a$)"), "* The number of possible scaffold crossovers  : "
        write(i, "(i7)"), dna.n_xover_scaf
        call Space(i, 11)
        write(i, "(a )"), "* Detailed information on continuous scaffold"
    end do

    ! Print progress in detail
    count = 1
    start = dna.base_scaf(1).id
    up    = dna.base_scaf(1).id

    do
        if(mod(count, 100) == 0) then
            write(11, "(i7 , a )"), up,    "->"
        else if(mod(count, 100) == 1) then
            write(11, "(i20, a$)"), count, " th base - "
            write(11, "(i7 , a$)"), up,    "->"
        else
            write(11, "(i7 , a$)"), up,    "->"
        end if

        ! Update bases
        up = dna.base_scaf(up).up
        if(up == start) exit

        count = count + 1
        if(count > dna.n_base_scaf) then
            write(0, "(a$)"), "Error - Not single circular scaffold : "
            write(0, "(a )"), "Route_Make_Scaf_Origami"
            stop
        end if
    end do
    write(0, "(a)"); write(11, "(a)")

    ! Check all visit
    if(count /= dna.n_base_scaf) then
        call Space(0, 11)
        write(0, "(a$)"), "Error - Not visited bases or strands : "
        write(0, "(a )"), "Route_Make_Scaf_Origami"
        write(0, "(a )")
        call Space(0, 6)
        write(0, "(a )"), " -----------------------------------------------------------------------------------"
        call Space(0, 6)
        write(0, "(a )"), " It cannot make one continuous scaffold strand with the possible scaffold crossovers"
        call Space(0, 6)
        write(0, "(a )")
        call Space(0, 6)
        write(0, "(a )"), "                       Check the possible scaffold crossovers       "
        call Space(0, 6)
        write(0, "(a )"), " -----------------------------------------------------------------------------------"
        write(0, "(a )")
        stop
    end if

    ! Deallocate memory
    deallocate(b_visit)
    deallocate(id_entry)
end subroutine Route_Make_Scaf_Origami

! ---------------------------------------------------------------------------------------

! Find possible staple crossovers
! Last updated on Friday 15 June 2016 by Hyungmin
subroutine Route_Find_Possible_Stap_Xover(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Crossover based on cross-sectional edges
    type :: CroLType
        integer :: max_bp, min_bp
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, k, croL_cur, croL_com, sec_cur, sec_com, id_bp
    integer :: up_scaf1, dn_scaf1, up_scaf2, dn_scaf2
    integer :: up_cur, up_com, dn_cur, dn_com
    logical :: b_nei_up, b_nei_dn, b_scaf

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.8. Set possible staple crossovers"
    end do

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    croL(1:geom.n_croL).max_bp = -999999
    croL(1:geom.n_croL).min_bp =  999999

    ! Find maximum and minimum basepair ID in cross-sectional edges
    do i = 1, mesh.n_node
        croL_cur = mesh.node(i).croL
        id_bp    = mesh.node(i).bp

        ! Set maximum and minimum base ID
        if(croL(croL_cur).max_bp < id_bp) croL(croL_cur).max_bp = id_bp
        if(croL(croL_cur).min_bp > id_bp) croL(croL_cur).min_bp = id_bp
    end do

    ! Find the possible staple double crossovers
    dna.n_xover_stap  = 0
    dna.n_sxover_stap = 0
    do i = 1, mesh.n_node       ! Loop for current node

        ! Print progress bar
        call Mani_Progress_Bar(i, mesh.n_node)

        ! Loop for comparing node
        do j = i + 1, mesh.n_node

            ! Exception for the pre-constructed crossovers (due to double crossover)
            if(dna.base_stap(i).xover /= -1 .and. dna.base_stap(j).xover /= -1) then
                cycle
            end if

            ! It should be skipped when condition below
            ! Basepair ID and iniL shoud be the same and croL and section ID should be different
            if(mesh.node(i).bp   /= mesh.node(j).bp  ) cycle
            if(mesh.node(i).iniL /= mesh.node(j).iniL) cycle
            if(mesh.node(i).croL == mesh.node(j).croL) cycle
            if(mesh.node(i).sec  == mesh.node(j).sec ) cycle

            ! Find section ID
            sec_cur  = mesh.node(i).sec
            sec_com  = mesh.node(j).sec
            croL_cur = mesh.node(i).croL
            croL_com = mesh.node(j).croL
            id_bp    = mesh.node(i).bp

            ! To eliminate boundary staple crossovers
            if(croL(croL_cur).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_cur).max_bp - para_gap_xover_bound_stap < id_bp) cycle
            if(croL(croL_com).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_com).max_bp - para_gap_xover_bound_stap < id_bp) cycle

            ! Determine whether the node has crossover or not
            if(Section_Connection_Stap(geom, sec_cur, sec_com, id_bp) == .true.) then

                ! To eliminate crossover if there is neighboring scaffold crossover
                !
                !     dn_scaf1   i   up_scaf1
                !    *---*---*---*---*---*---*-->  : node i
                !            |       |
                ! <--*---*---*---*---*---*---*     : node j
                !     up_scaf2   j   dn_scaf2
                b_scaf   = .false.
                up_scaf1 = mesh.node(i).id
                dn_scaf1 = mesh.node(i).id
                up_scaf2 = mesh.node(j).id
                dn_scaf2 = mesh.node(j).id

                ! Check neighbor scaffold crossovers
                do k = 1, para_gap_xover_two
                    if( (dna.base_scaf(mesh.node(up_scaf1).id).xover == dna.base_scaf(mesh.node(dn_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(dn_scaf2).id).xover == dna.base_scaf(mesh.node(up_scaf1).id).id ) .or.  &
                        (dna.base_scaf(mesh.node(dn_scaf1).id).xover == dna.base_scaf(mesh.node(up_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(up_scaf2).id).xover == dna.base_scaf(mesh.node(dn_scaf1).id).id ) ) then
                        b_scaf = .true.
                        exit
                    else
                        up_scaf1 = mesh.node(up_scaf1).up
                        dn_scaf1 = mesh.node(dn_scaf1).dn
                        up_scaf2 = mesh.node(up_scaf2).up
                        dn_scaf2 = mesh.node(dn_scaf2).dn
                    end if
                end do

                if(b_scaf == .true.) cycle

                ! Find upper or downward neighboring crossovers
                ! Node numbering is opposite to staple ID
                up_cur   = mesh.node(i).dn
                up_com   = mesh.node(j).up
                id_bp    = mesh.node(up_cur).bp
                b_nei_up = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                dn_cur   = mesh.node(i).up
                dn_com   = mesh.node(j).dn
                id_bp    = mesh.node(dn_cur).bp
                b_nei_dn = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                ! Set current and previous or next crossovers (double crossover)
                ! For current crossover
                dna.n_xover_stap       = dna.n_xover_stap + 2
                dna.base_stap(i).xover = dna.base_stap(j).id
                dna.base_stap(j).xover = dna.base_stap(i).id

                ! For neighboring crossover
                if(b_nei_up == .true.) then
                    dna.base_stap(up_cur).xover = dna.base_stap(up_com).id
                    dna.base_stap(up_com).xover = dna.base_stap(up_cur).id

                    ! Set connectivity
                    dna.base_stap(i).up      = dna.base_stap(j).id
                    dna.base_stap(j).dn      = dna.base_stap(i).id
                    dna.base_stap(up_cur).dn = dna.base_stap(up_com).id
                    dna.base_stap(up_com).up = dna.base_stap(up_cur).id

                else if(b_nei_dn == .true.) then
                    dna.base_stap(dn_cur).xover = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).xover = dna.base_stap(dn_cur).id

                    ! Set connectivity
                    dna.base_stap(i).dn      = dna.base_stap(j).id
                    dna.base_stap(j).up      = dna.base_stap(i).id
                    dna.base_stap(dn_cur).up = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).dn = dna.base_stap(dn_cur).id
                end if
            end if
        end do
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a)"), "* The number of possible staple crossovers    : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
        write(i, "(a)")
    end do

    ! Deallocate memory
    deallocate(croL)
end subroutine Route_Find_Possible_Stap_Xover

! ---------------------------------------------------------------------------------------

! Set staple crossover
! Last updated on Thursday 18 Feb 2016 by Hyungmin
subroutine Route_Set_Stap_Crossover(prob, dna)
    type(ProbType), intent(in)    :: prob
    type(DNAType),  intent(inout) :: dna

    integer :: i, cur, up, down
    integer :: cur_xover, up_xover, down_xover

    ! Print information
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.12. Set staple cross-overs"
    end do

    ! Loop for all bases of staple strand
    do i = 1, dna.n_base_stap

        ! Print progress bar
        call Mani_Progress_Bar(i, dna.n_base_stap)

        ! Current, up and down id
        cur  = dna.base_stap(i).id
        up   = dna.base_stap(cur).up
        down = dna.base_stap(cur).dn

        if(up == -1 .or. down == -1) cycle

        ! Crossover id
        cur_xover  = dna.base_stap(cur).xover
        up_xover   = dna.base_stap(up).xover
        down_xover = dna.base_stap(down).xover

        if(cur_xover /= -1 .and. (up /= cur_xover .or. down /= cur_xover)) then
            if(up_xover /= -1) then
                !     cur    up
                ! ---->*---->*---->
                !      |     |
                ! <----*<----*<----
                ! cur_xover  up_xover
                dna.base_stap(cur      ).up = cur_xover
                dna.base_stap(cur_xover).dn = cur
                dna.base_stap(up       ).dn = up_xover
                dna.base_stap(up_xover ).up = up
            else if(down_xover /= -1) then
                !     cur   down
                ! <----*<----*<----
                !      |     |
                ! ---->*---->*---->
                ! cur_xover  down_xover
                dna.base_stap(cur       ).dn = cur_xover
                dna.base_stap(cur_xover ).up = cur
                dna.base_stap(down      ).up = down_xover
                dna.base_stap(down_xover).dn = down
            else
                write(0, "(a$)"), "Error - The number of crossover should be even : "
                write(0, "(a )"), "Route_Set_Stap_Crossover"
                stop
            end if
        end if
    end do
    
    write(0, "(a)"); write(11, "(a)")
end subroutine Route_Set_Stap_Crossover

! ---------------------------------------------------------------------------------------

! Set three orientation vectors based on 3DNA convention
! Last updated on Thursday 26 May 2016 by Hyungmin
subroutine Route_Set_Orientation(mesh, dna)
    type(MeshType), intent(inout) :: mesh
    type(DNAType),  intent(in)    :: dna

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j, up, down, id

    ! For all base pairs
    do i = 1, mesh.n_node

        ! Check the base connectivity
        if(mesh.node(i).up /= -1) then
            pos_1(1:3) = mesh.node(i).pos(1:3)
            pos_2(1:3) = mesh.node(mesh.node(i).up).pos(1:3)
        else
            ! If the end base does not exist, using downward base
            pos_1(1:3) = mesh.node(mesh.node(i).dn).pos(1:3)
            pos_2(1:3) = mesh.node(i).pos(1:3)
        end if

        ! Third orientation vector, e3
        ! Along the duplex axis towards the 3'-direction of the strand with the preferred nucleotide
        mesh.node(i).ori(3,:) = Normalize_Vector(pos_2 - pos_1)

        ! Second orientation vector, e2 which indicates the preferred nucleotide
        mesh.node(i).ori(2,:) = dna.base_scaf(i).pos - mesh.node(i).pos
        mesh.node(i).ori(2,:) = Normalize_Vector(mesh.node(i).ori(2,:))

        ! First orientation vector, e1 which indicates the mojor groove
        mesh.node(i).ori(1,:) = Cross_Product(mesh.node(i).ori(2,:), mesh.node(i).ori(3,:))
        mesh.node(i).ori(1,:) = Normalize_Vector(mesh.node(i).ori(1,:))
    end do

    return

    ! For all base pairs
    do i = 1, mesh.n_node
        if(mesh.node(i).up == -1 .or. mesh.node(i).dn == -1) then

            ! Check the base connectivity
            if(mesh.node(i).up /= -1) then
                pos_1(1:3) = mesh.node(i).pos(1:3)
                pos_2(1:3) = mesh.node(mesh.node(i).up).pos(1:3)
            else
                ! If the end base does not exist, using downward base
                pos_1(1:3) = mesh.node(mesh.node(i).dn).pos(1:3)
                pos_2(1:3) = mesh.node(i).pos(1:3)
            end if

            ! Third orientation vector, e3
            ! Along the duplex axis towards the 3'-direction of the strand with the preferred nucleotide
            mesh.node(i).ori(3,:) = Normalize_Vector(pos_2 - pos_1)

            ! Second orientation vector, e2 which indicates the preferred nucleotide
            !mesh.node(i).ori(2,:) = dna.base_scaf(i).pos - mesh.node(i).pos
            !mesh.node(i).ori(2,:) = Normalize_Vector(mesh.node(i).ori(2,:))

            id = mesh.node(i).id
            if(mesh.node(i).dn == -1) then
                do j = 1, iabs(mesh.node(i).bp) - 1
                    id = mesh.node(id).up
                end do
            else if(mesh.node(i).up == -1) then
                do j = 1, iabs(mesh.node(i).bp) - 1
                    id = mesh.node(id).dn
                end do
            end if

            mesh.node(i).ori(2,:) = mesh.node(id).ori(2,:)

            ! First orientation vector, e1 which indicates the mojor groove
            mesh.node(i).ori(1,:) = Cross_Product(mesh.node(i).ori(2,:), mesh.node(i).ori(3,:))
            mesh.node(i).ori(1,:) = Normalize_Vector(mesh.node(i).ori(1,:))
        end if
    end do
end subroutine Route_Set_Orientation

! ---------------------------------------------------------------------------------------

! Write crossovers based on base pairs
! Last updated on Friday 05 August 2016 by Hyungmin
subroutine Route_Chimera_Crossovers(prob, geom, bound, mesh, dna)
    type(ProbType),  intent(in) :: prob
    type(GeomType),  intent(in) :: geom
    type(BoundType), intent(in) :: bound
    type(MeshType),  intent(in) :: mesh
    type(DNAType),   intent(in) :: dna

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j, xover
    logical :: f_axis
    character(200) :: path

    if(para_write_607 == .false.) return

    f_axis = para_chimera_axis
    path   = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=607, file=trim(path)//"_crossovers.bild", form="formatted")

    ! write vertex connection
    do i = 1, bound.n_junc
        do j = 1, geom.n_sec*bound.junc(i).n_arm

            ! Cylinder color depending on vertex connection
            if( mesh.node(bound.junc(i).conn(j, 1)).conn == 2 .or. &
                mesh.node(bound.junc(i).conn(j, 1)).conn == 4 ) then
                write(607, "(a)"), ".color red"
            else if( mesh.node(bound.junc(i).conn(j, 1)).conn == 1 .or. &
                     mesh.node(bound.junc(i).conn(j, 1)).conn == 3 ) then
                write(607, "(a)"), ".color blue"
            end if

            pos_1(:) = mesh.node(bound.junc(i).conn(j, 1)).pos
            pos_2(:) = mesh.node(bound.junc(i).conn(j, 2)).pos

            write(607, "(a$    )"), ".cylinder "
            write(607, "(3f9.3$)"), pos_1(1:3)
            write(607, "(3f9.3$)"), pos_2(1:3)
            write(607, "(1f9.3 )"), 0.06d0
        end do
    end do

    ! Write base pair
    write(607, "(a)"), ".color tan"
    do i = 1, mesh.n_node
        write(607, "(a$    )"), ".sphere "
        write(607, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(607, "(1f9.3 )"), 0.08d0
    end do

    ! First orientation vector, e3
    ! Along the duplex axis towards the 3'-direction of the strand with the preferred nucleotide
    write(607, "(a)"), ".color steel blue"
    do i = 1, mesh.n_node
        write(607, "(a$    )"), ".arrow "
        write(607, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(607, "(3f8.2$)"), mesh.node(i).pos(1:3) + mesh.node(i).ori(3,:)*0.25d0
        write(607, "(3f8.2 )"), 0.04d0, 0.12d0, 0.6d0
    end do

    ! Write centered scaffold crossovers
    write(607, "(a)"), ".color steel blue"
    do i = 1, mesh.n_node

        xover = dna.base_scaf(i).xover
        if(xover /= -1 .and. i < xover) then

            pos_1(:) = mesh.node(i).pos
            pos_2(:) = mesh.node(xover).pos

            write(607, "(a$    )"), ".cylinder "
            write(607, "(3f9.3$)"), pos_1(1:3)
            write(607, "(3f9.3$)"), pos_2(1:3)
            write(607, "(1f9.3 )"), 0.06d0
        end if
    end do

    ! Write possible staple crossovers
    write(607, "(a)"), ".color orange"
    do i = 1, mesh.n_node

        xover = dna.base_stap(i).xover
        if(xover /= -1 .and. i < xover) then

            pos_1(:) = mesh.node(i).pos(1:3)
            pos_2(:) = mesh.node(xover).pos(1:3)

            write(607, "(a$    )"), ".cylinder "
            write(607, "(3f9.3$)"), pos_1(1:3)
            write(607, "(3f9.3$)"), pos_2(1:3)
            write(607, "(1f9.3 )"), 0.06d0
        end if
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(607, "(a)"), ".translate 0.0 0.0 0.0"
        write(607, "(a)"), ".scale 0.5"
        write(607, "(a)"), ".color grey"
        write(607, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(607, "(a)"), ".color red"             ! x-axis
        write(607, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(607, "(a)"), ".color blue"            ! y-axis
        write(607, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(607, "(a)"), ".color yellow"          ! z-axis
        write(607, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=607)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=607, file=trim(path)//"_crossovers.dat", form="formatted")

    write(607, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(607, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(607, "(a$)"), 'ZONE F = FEPOINT'
    write(607, "(a$)"), ', N='//trim(adjustl(Int2Str(mesh.n_node)))
    write(607, "(a$)"), ', E='//trim(adjustl(Int2Str(mesh.n_ele)))
    write(607, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, mesh.n_node
        write(607, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(607, "(1f9.3 )"), 1.0d0
    end do

    ! Write elements
    do i = 1, mesh.n_ele
        write(607, "(2i7)"), mesh.ele(i).cn(1), mesh.ele(i).cn(2)
    end do

    ! Write scaffold crossovers
    write(607, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(607, "(a$)"), 'ZONE F = FEPOINT'
    write(607, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.n_xover_scaf*2)))
    write(607, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.n_xover_scaf)))
    write(607, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, mesh.n_node

        xover = dna.base_scaf(i).xover
        if(xover /= -1 .and. i < xover) then

            pos_1(:) = mesh.node(i).pos
            pos_2(:) = mesh.node(xover).pos

            write(607, "(4f9.3)"), pos_1(1:3), 1.0d0
            write(607, "(4f9.3)"), pos_2(1:3), 1.0d0
        end if
    end do

    ! Write elements
    do i = 1, dna.n_xover_scaf
        write(607, "(2i7)"), 2*i-1, 2*i
    end do

    ! Write staple crossovers
    write(607, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(607, "(a$)"), 'ZONE F = FEPOINT'
    write(607, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.n_xover_stap*2)))
    write(607, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.n_xover_stap)))
    write(607, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, mesh.n_node

        xover = dna.base_stap(i).xover
        if(xover /= -1 .and. i < xover) then

            pos_1(:) = mesh.node(i).pos
            pos_2(:) = mesh.node(xover).pos

            write(607, "(4f9.3)"), pos_1(1:3), 1.0d0
            write(607, "(4f9.3)"), pos_2(1:3), 1.0d0
        end if
    end do

    ! Write elements
    do i = 1, dna.n_xover_stap
        write(607, "(2i7)"), 2*i-1, 2*i
    end do

    close(unit=607)
end subroutine Route_Chimera_Crossovers

! ---------------------------------------------------------------------------------------

! Write 3 orientation vectors based on 3DNA convention
! Last updated on Saturday 16 July 2016 by Hyungmin
subroutine Route_Chimera_Orientation(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    character(200) :: path
    logical :: f_axis
    integer :: i

    if(para_write_608 == .false.) return

    f_axis = para_chimera_axis

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=608, file=trim(path)//"_orientation.bild", form="formatted")

    ! Write base pair
    write(608, "(a)"), ".color tan"
    do i = 1, mesh.n_node
        write(608, "(a$    )"), ".sphere "
        write(608, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(608, "(1f9.3 )"), 0.08d0
    end do

    ! First orientation vector, e1 which indicates the mojor groove
    write(608, "(a)"), ".color salmon"
    do i = 1, mesh.n_node
        write(608, "(a$    )"), ".arrow "
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3) + mesh.node(i).ori(1,:)*0.85d0
        write(608, "(2f8.2 )"), 0.03d0, 0.09d0
    end do

    ! Second orientation vector, e2 which indicates the preferred nucleotide
    write(608, "(a)"), ".color sea green"
    do i = 1, mesh.n_node
        write(608, "(a$    )"), ".arrow "
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3) + mesh.node(i).ori(2,:)*0.85d0
        write(608, "(2f8.2 )"), 0.03d0, 0.09d0
    end do

    ! First orientation vector, e3
    ! Along the duplex axis towards the 3'-direction of the strand with the preferred nucleotide
    write(608, "(a)"), ".color steel blue"
    do i = 1, mesh.n_node
        write(608, "(a$    )"), ".arrow "
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(608, "(3f8.2$)"), mesh.node(i).pos(1:3) + mesh.node(i).ori(3,:)*0.25d0
        write(608, "(3f8.2 )"), 0.04d0, 0.1d0, 0.6d0
    end do

    ! Write global axis
    if(f_axis == .true.) then
        write(608, "(a)"), ".translate 0.0 0.0 0.0"
        write(608, "(a)"), ".scale 0.5"
        write(608, "(a)"), ".color grey"
        write(608, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(608, "(a)"), ".color red"             ! x-axis
        write(608, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(608, "(a)"), ".color blue"            ! y-axis
        write(608, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(608, "(a)"), ".color yellow"          ! z-axis
        write(608, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if
    close(unit=608)

    ! ---------------------------------------------
    !
    ! Write the file for Tecplot
    !
    ! ---------------------------------------------
    if(para_output_Tecplot == "off") return

    path = trim(prob.path_work1)//"Tecplot\"//trim(prob.name_file)
    open(unit=608, file=trim(path)//"_orientation.dat", form="formatted")

    write(608, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(608, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(608, "(a$)"), 'ZONE F = FEPOINT'
    write(608, "(a$)"), ', N='//trim(adjustl(Int2Str(mesh.n_node)))
    write(608, "(a$)"), ', E='//trim(adjustl(Int2Str(mesh.n_ele)))
    write(608, "(a )"), ', ET=LINESEG'

    ! Write nodes
    do i = 1, mesh.n_node
        write(608, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(608, "(3f9.3 )"), 1.0d0, 1.0d0, 1.0d0
    end do

    ! Write elements
    do i = 1, mesh.n_ele
        write(608, "(2i7)"), mesh.ele(i).cn(1), mesh.ele(i).cn(2)
    end do

    ! Write staple crossovers
    write(608, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(608, "(a$)"), 'ZONE F = FEPOINT'
    write(608, "(a$)"), ', N='//trim(adjustl(Int2Str(mesh.n_node)))
    write(608, "(a$)"), ', E='//trim(adjustl(Int2Str(mesh.n_ele)))
    write(608, "(a )"), ', ET=LINESEG'

    ! Write nodes with orientation vector 1
    do i = 1, mesh.n_node
        write(608, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(608, "(3f9.3 )"), mesh.node(i).ori(1, 1:3)
    end do

    ! Write elements
    do i = 1, mesh.n_ele
        write(608, "(2i7)"), mesh.ele(i).cn(1), mesh.ele(i).cn(2)
    end do

    ! Write staple crossovers
    write(608, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(608, "(a$)"), 'ZONE F = FEPOINT'
    write(608, "(a$)"), ', N='//trim(adjustl(Int2Str(mesh.n_node)))
    write(608, "(a$)"), ', E='//trim(adjustl(Int2Str(mesh.n_ele)))
    write(608, "(a )"), ', ET=LINESEG'

    ! Write nodes with orientation vector 2
    do i = 1, mesh.n_node
        write(608, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(608, "(3f9.3 )"), mesh.node(i).ori(2, 1:3)
    end do

    ! Write elements
    do i = 1, mesh.n_ele
        write(608, "(2i7)"), mesh.ele(i).cn(1), mesh.ele(i).cn(2)
    end do

    ! Write staple crossovers
    write(608, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(608, "(a$)"), 'ZONE F = FEPOINT'
    write(608, "(a$)"), ', N='//trim(adjustl(Int2Str(mesh.n_node)))
    write(608, "(a$)"), ', E='//trim(adjustl(Int2Str(mesh.n_ele)))
    write(608, "(a )"), ', ET=LINESEG'

    ! Write nodes with orientation vector 3
    do i = 1, mesh.n_node
        write(608, "(3f9.3$)"), mesh.node(i).pos(1:3)
        write(608, "(3f9.3 )"), mesh.node(i).ori(3, 1:3)
    end do

    ! Write elements
    do i = 1, mesh.n_ele
        write(608, "(2i7)"), mesh.ele(i).cn(1), mesh.ele(i).cn(2)
    end do

    close(unit=608)
end subroutine Route_Chimera_Orientation

! ---------------------------------------------------------------------------------------

! Write atomic model
! Last updated on Thuesday 9 August 2016 by Hyungmin
subroutine Route_Chimera_Atom(prob, geom, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(DNAType),  intent(in) :: dna

    double precision :: pos_1(3), pos_2(3)
    integer :: i, across
    logical :: f_axis, f_dir, f_cyn
    character(200) :: path

    if(para_write_609 == .false.) return

    f_axis = para_chimera_axis
    f_cyn  = para_chimera_609_cyl
    f_dir  = para_chimera_609_dir

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=609, file=trim(path)//"_atom.bild", form="formatted")

    ! --------------------------------------------------
    !
    ! For the nucleotide of scaffold strand
    !
    ! --------------------------------------------------
    write(609, "(a)"), ".color steel blue"
    do i = 1, dna.n_base_scaf

        ! Write bases
        write(609, "(a$    )"), ".sphere "
        write(609, "(3f9.3$)"), dna.base_scaf(i).pos(1:3)
        write(609, "(1f9.3 )"), 0.15d0

        ! Write backbone
        if(dna.base_scaf(i).up /= -1) then
            pos_1(1:3) = dna.base_scaf(i).pos(1:3)
            pos_2(1:3) = dna.base_scaf(dna.base_scaf(i).up).pos

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), pos_1(1:3)
            write(609, "(3f9.3$)"), pos_2(1:3)
            write(609, "(1f9.3 )"), 0.05d0
        end if
    end do

    ! --------------------------------------------------
    !
    ! For the nucleotide of staple strand
    !
    ! --------------------------------------------------
    write(609, "(a)"), ".color orange"
    do i = 1, dna.n_base_stap

        ! Write bases
        write(609, "(a$    )"), ".sphere "
        write(609, "(3f9.3$)"), dna.base_stap(i).pos(1:3)
        write(609, "(1f9.3 )"), 0.15d0

        ! Write backbone
        if(dna.base_stap(i).up /= -1) then
            pos_1(1:3) = dna.base_stap(i).pos(1:3)
            pos_2(1:3) = dna.base_stap(dna.base_stap(i).up).pos

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), pos_1(1:3)
            write(609, "(3f9.3$)"), pos_2(1:3)
            write(609, "(1f9.3 )"), 0.05d0
        end if
    end do

    ! --------------------------------------------------
    !
    ! Write crossovers of the scaffold and staple strand
    !
    ! --------------------------------------------------
    write(609, "(a)"), ".color blue"
    do i = 1, dna.n_base_scaf
        if(dna.base_scaf(i).xover /= -1 .and. i < dna.base_scaf(i).xover) then
            pos_1(:) = dna.base_scaf(i).pos
            pos_2(:) = dna.base_scaf(dna.base_scaf(i).xover).pos

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), pos_1(1:3)
            write(609, "(3f9.3$)"), pos_2(1:3)
            write(609, "(1f9.3 )"), 0.08d0
        end if
    end do

    write(609, "(a)"), ".color red"
    do i = 1, dna.n_base_stap
        if(dna.base_stap(i).xover /= -1 .and. i < dna.base_stap(i).xover) then
            pos_1(:) = dna.base_stap(i).pos
            pos_2(:) = dna.base_stap(dna.base_stap(i).xover).pos

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), pos_1(1:3)
            write(609, "(3f9.3$)"), pos_2(1:3)
            write(609, "(1f9.3 )"), 0.08d0
        end if
    end do

    ! --------------------------------------------------
    !
    ! Write cylinder or Watson-Crick base pair connections
    !
    ! --------------------------------------------------
    if(f_cyn == .true.) then

        ! for cylinder
        write(609, "(a)"), ".color grey"
        do i = 1, geom.n_croL
            pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos
            pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), pos_1(1:3)
            write(609, "(3f9.3$)"), pos_2(1:3)
            write(609, "(1f9.3 )"), 0.85d0
        end do
    else

        ! For Watson-Crick base pair
        write(609, "(a)"), ".color grey"
        do i = 1, dna.n_base_scaf

            across = dna.base_scaf(i).across
            if(across == -1) cycle

            write(609, "(a$    )"), ".cylinder "
            write(609, "(3f9.3$)"), dna.base_scaf(i).pos
            write(609, "(3f9.3$)"), dna.base_stap(across).pos
            write(609, "(1f9.3 )"), 0.03d0
        end do
    end if

    ! Write direction
    if(f_dir == .true.) then
        do i = 1, geom.n_croL
            if(mod(geom.croL(i).sec, 2) == 0) then
                write(609, "(a)"), ".color purple"
                pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos
                pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos
            else
                write(609, "(a)"), ".color dark green"
                pos_1(1:3) = geom.croP(geom.croL(i).poi(2)).pos
                pos_2(1:3) = geom.croP(geom.croL(i).poi(1)).pos
            end if
            write(609, "(a$    )"), ".arrow "
            write(609, "(3f8.2$)"), pos_1(1:3)
            write(609, "(3f8.2$)"), pos_2(1:3)
            write(609, "(1f8.2$)"), 0.08d0
            write(609, "(1f8.2 )"), 0.25d0
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) then
        write(609, "(a)"), ".translate 0.0 0.0 0.0"
        write(609, "(a)"), ".scale 0.5"
        write(609, "(a)"), ".color grey"
        write(609, "(a)"), ".sphere 0 0 0 0.5"      ! Center
        write(609, "(a)"), ".color red"             ! x-axis
        write(609, "(a)"), ".arrow 0 0 0 4 0 0 "
        write(609, "(a)"), ".color blue"            ! y-axis
        write(609, "(a)"), ".arrow 0 0 0 0 4 0 "
        write(609, "(a)"), ".color yellow"          ! z-axis
        write(609, "(a)"), ".arrow 0 0 0 0 0 4 "
    end if

    close(unit=609)
end subroutine Route_Chimera_Atom

! ---------------------------------------------------------------------------------------

! Delete possible scaffold crossovers near the junction
! Last updated on Wednesday 2 Mar 2016 by Hyungmin
subroutine Route_Delete_End_Pos_Xover_Scaf(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, id, xover, sec_1, sec_2, n_xover_scaf_pre
    integer, parameter :: para_avoid_gap_scaf = 7

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.8. Deleting possible crossovers in scaffold strand"
        call Space(i, 11)
        write(i, "(a$)"), "* Parameter for deleting xovers near junction : "
        write(i, "(i7)"), para_avoid_gap_scaf
    end do

    ! Original # of possible scaf crossovers
    n_xover_scaf_pre = dna.n_xover_scaf

    ! Loop to find negative upIDvalue
    do i = 1, mesh.n_node

        ! Search for end of edges
        if(mesh.node(i).dn == -1) then

            ! If searching node has -1 in id_down
            ! Go up and delete up to para_avoid_gap_scaf
            ! *----->*----->*----->*----->*
            ! ↑
            ! node(i)
            id = mesh.node(i).id

            do j = 1, para_avoid_gap_scaf

                xover = dna.base_scaf(id).xover

                ! If there are crossovers, delete crossovers
                if(xover /= -1) then
                    dna.base_scaf(id).xover    = -1
                    dna.base_scaf(xover).xover = -1
                    dna.n_xover_scaf = dna.n_xover_scaf - 1

                    ! Delete one crossover when crossover exists at para_avoid_gap_scaf
                    if(j == para_avoid_gap_scaf) then
                        sec_1 = mesh.node(xover).sec
                        id    = mesh.node(id).up
                        xover = dna.base_scaf(id).xover

                        if(xover == -1) cycle
                        sec_2 = mesh.node(xover).sec

                        ! If there are crossovers at the same section, delete crossovers
                        if(xover /= -1 .and. sec_1 == sec_2) then
                            dna.base_scaf(id).xover    = -1
                            dna.base_scaf(xover).xover = -1
                            dna.n_xover_scaf = dna.n_xover_scaf - 1
                        end if
                    end if
                end if

                ! Update ID for next node
                id = mesh.node(id).up
            end do

        else if(mesh.node(i).up == -1) then
            !
            ! If searching node has -1 in id_down
            ! Go up and delete up to para_avoid_gap_scaf
            ! *<------*<-----*<-----*<-----*
            ! ↑
            ! node(i)
            id = mesh.node(i).id

            do j = 1, para_avoid_gap_scaf

                xover = dna.base_scaf(id).xover

                ! If there are crossovers, delete crossovers
                if(xover /= -1) then
                    dna.base_scaf(id).xover    = -1
                    dna.base_scaf(xover).xover = -1
                    dna.n_xover_scaf = dna.n_xover_scaf - 1

                    ! Delete one crossover when crossover exists at para_avoid_gap_scaf
                    if(j == para_avoid_gap_scaf) then
                        sec_1 = mesh.node(xover).sec
                        id    = mesh.node(id).dn
                        xover = dna.base_scaf(id).xover

                        if(xover == -1) cycle
                        sec_2 = mesh.node(xover).sec

                        ! If there are crossovers at the same section, delete crossovers
                        if(xover /= -1 .and. sec_1 == sec_2) then
                            dna.base_scaf(id).xover    = -1
                            dna.base_scaf(xover).xover = -1
                            dna.n_xover_scaf = dna.n_xover_scaf - 1
                        end if
                    end if
                end if

                ! Update ID for next node
                id = mesh.node(id).dn
            end do
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a$)"), "* The number of original possible crossovers  : "
        write(i, "(i7)"), n_xover_scaf_pre
        call Space(i, 11)
        write(i, "(a$)"), "* The number of deleted possible crossovers   : "
        write(i, "(i7)"), n_xover_scaf_pre - dna.n_xover_scaf
        call Space(i, 11)
        write(i, "(a$)"), "* The number of current possible crossovers   : "
        write(i, "(i7)"), dna.n_xover_scaf
        write(i, "(a )")
    end do
end subroutine Route_Delete_End_Pos_Xover_Scaf

! ---------------------------------------------------------------------------------------

! Delete possible staple crossovers near the junction
! Last updated on Wednesday 2 Mar 2016 by Hyungmin
subroutine Route_Delete_End_Pos_Xover_Stap(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, id, xover, sec_1, sec_2, n_xover_stap_pre

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a )"), "5.9. Deleting possible crossovers in staple strand"
        call Space(i, 11)
        write(i, "(a$)"), "* Parameter for deleting xovers near junction : "
        write(i, "(i7)"), para_gap_xover_bound_stap
    end do

    ! Original # of possible stap crossovers
    n_xover_stap_pre = dna.n_xover_stap

    ! Find up or down has negative value
    do i = 1, mesh.n_node

        ! Search for end of edges
        if(mesh.node(i).dn == -1) then
            !
            ! If searching node has -1 in id_down
            ! Go up and delete up to para_gap_xover_bound_stap
            ! *----->*----->*----->*----->*
            ! ↑
            ! node(i)
            id = mesh.node(i).id

            do j = 1, para_gap_xover_bound_stap

                xover = dna.base_stap(id).xover

                ! If there are crossovers, delete crossovers
                if(xover /= -1) then
                    dna.base_stap(id).xover    = -1
                    dna.base_stap(xover).xover = -1
                    dna.n_xover_stap  = dna.n_xover_stap  - 1
                    dna.n_sxover_stap = dna.n_sxover_stap + 1

                    ! Delete one crossover when crossover exists at para_gap_xover_bound_stap
                    if(j == para_gap_xover_bound_stap) then
                        sec_1 = mesh.node(xover).sec
                        id    = mesh.node(id).up
                        xover = dna.base_stap(id).xover

                        if(xover == -1) cycle
                        sec_2 = mesh.node(xover).sec

                        ! If there are crossovers, delete crossovers
                        if(xover /= -1 .and. sec_1 == sec_2) then
                            dna.base_stap(id).xover    = -1
                            dna.base_stap(xover).xover = -1
                            dna.n_xover_stap  = dna.n_xover_stap  - 1
                            dna.n_sxover_stap = dna.n_sxover_stap + 1
                        end if
                    end if
                end if

                id = mesh.node(id).up
            end do

        else if(mesh.node(i).up == -1) then
            !
            ! If searching node has -1 in id_up
            ! Go up and delete up to para_gap_xover_bound_stap
            ! *<-----*<-----*<-----*<-----*
            ! ↑
            ! node(i)
            id = mesh.node(i).id

            do j = 1, para_gap_xover_bound_stap

                xover = dna.base_stap(id).xover

                ! If there are crossovers, delete crossovers
                if(xover /= -1) then
                    dna.base_stap(id).xover    = -1
                    dna.base_stap(xover).xover = -1
                    dna.n_xover_stap  = dna.n_xover_stap  - 1
                    dna.n_sxover_stap = dna.n_sxover_stap + 1

                    ! Delete one crossover when crossover exists at para_gap_xover_bound_stap
                    if(j == para_gap_xover_bound_stap) then
                        sec_1 = mesh.node(xover).sec
                        id    = mesh.node(id).dn
                        xover = dna.base_stap(id).xover

                        if(xover == -1) cycle
                        sec_2 = mesh.node(xover).sec

                        ! If there are crossovers, delete crossovers
                        if(xover /= -1 .and. sec_1 == sec_2) then
                            dna.base_stap(id).xover    = -1
                            dna.base_stap(xover).xover = -1
                            dna.n_xover_stap  = dna.n_xover_stap  - 1
                            dna.n_sxover_stap = dna.n_sxover_stap + 1
                        end if
                    end if
                end if

                id = mesh.node(id).dn
            end do
        end if
    end do

    ! Print progress
    do i = 0, 11, 11
        call Space(i, 11)
        write(i, "(a$)"), "* The number of original possible crossovers  : "
        write(i, "(i7)"), n_xover_stap_pre
        call Space(i, 11)
        write(i, "(a$)"), "* The number of deleted possible crossovers   : "
        write(i, "(i7)"), n_xover_stap_pre - dna.n_xover_stap
        call Space(i, 11)
        write(i, "(a$)"), "* The number of current staple crossovers     : "
        write(i, "(i7)"), dna.n_xover_stap
        write(i, "(a$)"), "* The number of single staple crossovers      : "
        write(i, "(i7)"), dna.n_sxover_stap
        write(i, "(a )")
    end do
end subroutine Route_Delete_End_Pos_Xover_Stap

! ---------------------------------------------------------------------------------------
    
! Centered scaffold crossover
! Last updated on Monday 11 Apr 2016 by Hyungmin
subroutine Route_Center_Scaf_Crossover(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    
    ! Cross-section information on crossovers
    type :: CroLType
        integer, allocatable, dimension(:) :: cn        ! The section connectivity
        integer, allocatable, dimension(:) :: count     ! Counter
    end type CroLType

    type(CroLType), allocatable :: croL(:)
    integer :: cur_xover, cur_node, cur_sec, cur_croL
    integer :: i, j, num_xover, com_node, com_sec, value
    logical :: b_new, cur_strd, com_strd

    ! Print information
    do i = 0, 11, 11
        call Space(i, 6)
        write(i, "(a)"), "5.11. Delete boundary crossovers in scaffolds"
        call Space(i, 11)
        write(i, "(a$)"), "* The number of original crossovers           : "
        write(i, "(i7)"), dna.n_xover_scaf
    end do

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    do i = 1, geom.n_croL
        allocate(croL(i).cn(geom.n_sec))
        allocate(croL(i).count(geom.n_sec))
        croL(i).cn(1:geom.n_sec)    = 0
        croL(i).count(1:geom.n_sec) = 0
    end do

    ! Count the number of crossovers for each edge
    do i = 1, dna.n_base_scaf

        cur_xover = dna.base_scaf(i).xover
        if(cur_xover /= -1) then
            cur_node = dna.base_scaf(i).node
            cur_croL = mesh.node(cur_node).croL

            com_node = dna.base_scaf(cur_xover).node
            com_sec  = mesh.node(com_node).sec

            croL(cur_croL).cn(com_sec + 1) &
                = croL(cur_croL).cn(com_sec + 1) + 1
        end if
    end do

    !do i = 1, geom.n_croL
    !    write(0, "(i$)"), i
    !    do j = 1, geom.n_sec
    !        if(croL(i).cn(j) /= 0) then
    !            value =croL(i).cn(j) / 2
    !            if(mod(value, 2) == 1) then
    !                croL(i).cn(j) = value
    !            else
    !                croL(i).cn(j) = value - 1
    !            end if
    !        end if
    !        write(0, "(i6$)"), croL(i).cn(j)
    !    end do
    !    write(0, "(a)")
    !end do

    ! Delete boundary crossovers to remain centered crossovers
    do i = 1, dna.n_base_scaf
        cur_node  = dna.base_scaf(i).node
        cur_croL  = mesh.node(cur_node).croL
        cur_sec   = mesh.node(cur_node).sec
        cur_xover = dna.base_scaf(i).xover

        ! For positive z-direction and base with crossover
        if(mod(cur_sec, 2) /= 0 .or. cur_xover == -1) cycle

        com_node = dna.base_scaf(cur_xover).node
        com_sec  = mesh.node(com_node).sec

        croL(cur_croL).count(com_sec + 1) &
            = croL(cur_croL).count(com_sec + 1) + 1

        ! Find center position
        num_xover = croL(cur_croL).cn(com_sec + 1) / 2
        if(num_xover >= 2 .and. mod(num_xover, 2) == 0) then
            num_xover = num_xover - 1
        end if

        ! Delete crossovers
        if( num_xover     /= croL(cur_croL).count(com_sec + 1) .and. &
            num_xover + 1 /= croL(cur_croL).count(com_sec + 1)) then
            dna.base_scaf(i).xover         = -1
            dna.base_scaf(cur_xover).xover = -1
            dna.n_xover_scaf = dna.n_xover_scaf - 1
        end if
    end do

    ! Deallocate memory
    do i = 1, geom.n_croL
        deallocate(croL(i).cn)
        deallocate(croL(i).count)
    end do
    deallocate(croL)
end subroutine Route_Center_Scaf_Crossover

! ---------------------------------------------------------------------------------------

end module Route