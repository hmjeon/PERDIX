!
! ---------------------------------------------------------------------------------------
!
!                               Module for Exam_OpenGeo
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/10/20
!                                                            Last  modified  : 2016/08/30
!
! ---------------------------------------------------------------------------------------
!
module Exam_3D_Open

    use Data_Prob
    use Data_Geom

    use Para
    use Mani
    use Math

    implicit none

    public  Exam_Open_End_Cylinder             ! 53. Plate with uniform quadrilater mesh

    contains

! ---------------------------------------------------------------------------------------

! Example of retagular plate geometry with the uniform mesh of quadrilaterals
! Last updated on Tuesday 30 August 2016 by Hyungmin
subroutine Exam_Open_End_Cylinder(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: n_i_point, n_j_point, n_i_face, n_j_face
    integer :: i, j, n, numbering
    character(10) :: char_sec, char_bp, char_start_bp

    double precision :: young, possion, thick
    double precision :: radius, length
    double precision :: angle, layer_length
    integer :: domainsize, nz, nr, jointN, elementN

    write(unit=char_sec,      fmt = "(i10)"), prob.sel_sec
    write(unit=char_bp,       fmt = "(i10)"), prob.n_bp_edge
    write(unit=char_start_bp, fmt = "(i10)"), para_start_bp_ID

    prob.name_file = "1_Plate_Uniform_Quad"//&
        "_"//trim(adjustl(trim(char_sec)))//"cs"//&     ! Cross-section
        "_"//trim(adjustl(trim(char_bp)))//"bp"//&      ! Edge length
        "_"//trim(para_vertex_design)//&                ! Vertex design
        "_"//trim(para_vertex_modify)//&                ! Vertex modification
        "_"//trim(para_cut_stap_method)                 ! Cutting method

    prob.name_prob = "Plate Uniform Quad"

    ! Set geometric type and view
    prob.color    = [52, 152, 219]
    prob.scale    = 1.0d0      ! Atomic model
    prob.size     = 1.0d0      ! Cylindrical model
    prob.move_x   = 0.0d0      ! Cylindrical model
    prob.move_y   = 0.0d0      ! Cylindrical model
    prob.type_geo = "open"
    if(para_fig_view == "preset") para_fig_view = "XY"
    
    domainsize   =  6
    radius       =  4.953d0                 ! radius
    layer_length =  10.35d0 * 2.0d0
    nz           =  domainsize
    nr           =  2 * domainsize
    geom.n_face  =  nz * nr * 2
    geom.n_iniP  =  (nz + 1) * nr

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! set nodal position vector
    jointN = 0
    do i = 1, nz + 1
	    do j = 1, nr
            jointN = jointN + 1
            angle  = (360.0d0-(j-1)*(360.0d0/nr))*(3.141592653589793d0/180.0d0)
		    
            geom.iniP(jointN).pos(1) =  radius*dcos(angle)           ! x-coordinate
            geom.iniP(jointN).pos(2) = -(i-1)*(layer_length/nz)      ! y-coordinate
            geom.iniP(jointN).pos(3) =  radius*dsin(angle)           ! z-coordinate

	    end do
    end do
    
    ! set connectivity
    !elementN = 0
    !do i = 1, nz
    !    do j = 1, nr - 1
    !        elementN = elementN + 1
    !        geom.face(elementN).n_poi = 4
    !        allocate(geom.face(elementN).poi(4))
    !        
    !        geom.face(elementN).poi(1) = (i - 1) * nr + j
    !        geom.face(elementN).poi(2) = (i - 1) * nr + j + 1
    !        geom.face(elementN).poi(3) = i * nr + j + 1
    !        geom.face(elementN).poi(4) = i * nr + j
    !        
    !        print *, geom.face(elementN).poi(1:4)
    !    end do

    !    elementN = elementN + 1
    !    geom.face(elementN).n_poi = 4
    !    allocate(geom.face(elementN).poi(4))
    !    geom.face(elementN).poi(1) = (i - 1) * nr + j
	!    geom.face(elementN).poi(2) = (i - 1) * nr + 1
	!    geom.face(elementN).poi(3) = i * nr + 1
    !    geom.face(elementN).poi(4) = i * nr + j
    !end do
    
    
    
    
    ! set connectivity
    elementN = 0
    do i = 1, nz
        do j = 1, nr - 1
            elementN = elementN + 1
                    geom.face(elementN).n_poi = 3
        allocate(geom.face(elementN).poi(3))
            geom.face(elementN).poi(1) = (i - 1) * nr + j
            geom.face(elementN).poi(2) = (i - 1) * nr + j + 1
            geom.face(elementN).poi(3) = i * nr + j

            elementN = elementN + 1
                    geom.face(elementN).n_poi = 3
        allocate(geom.face(elementN).poi(3))
            geom.face(elementN).poi(1) = (i - 1) * nr + j + 1
            geom.face(elementN).poi(2) = i * nr + j + 1
            geom.face(elementN).poi(3) = i * nr + j
        end do

        elementN = elementN + 1
                geom.face(elementN).n_poi = 3
        allocate(geom.face(elementN).poi(3))
        geom.face(elementN).poi(1) = (i - 1) * nr + j
	    geom.face(elementN).poi(2) = (i - 1) * nr + 1
        geom.face(elementN).poi(3) = i * nr + j

        elementN = elementN + 1
                geom.face(elementN).n_poi = 3
        allocate(geom.face(elementN).poi(3))
	    geom.face(elementN).poi(1) = (i - 1) * nr + 1
	    geom.face(elementN).poi(2) = i * nr + 1
        geom.face(elementN).poi(3) = i * nr + j
    end do
    
    
    
!stop
end subroutine Exam_Open_End_Cylinder

! ---------------------------------------------------------------------------------------

end module Exam_3D_Open