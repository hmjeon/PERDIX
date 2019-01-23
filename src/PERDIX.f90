!
! =============================================================================
!
! PERDIX v2.0
! Last Updated : 01/22/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! PERDIX is an open-source software, which allows scientists to build and
! solve the sequence design of complex 2D DNA wireframe with 6HB edges.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
program PERDIX

    use Para            ! Parameters

    use Data_Prob       ! Problem data
    use Data_Geom       ! Geometry data
    use Data_Mesh       ! Basepair data
    use Data_DNA        ! B-form DNA data

    use Input           ! Input
    use ModGeo          ! Geometry modification
    use Section         ! Edge cross-section
    use Basepair        ! Basepair design
    use Route           ! Scaffold routing
    use SeqDesign       ! Sequence design
    use Output          ! Outputs

    implicit none

    call Main

contains

! -----------------------------------------------------------------------------

! Main subroutine
subroutine Main()

    ! Declare variables
    type(ProbType)  :: prob     ! Problem data
    type(GeomType)  :: geom     ! Geometry data
    type(MeshType)  :: mesh     ! Basepair data
    type(DNAType)   :: dna      ! B-form DNA data
    real            :: time     ! Computational time

    ! Initialize input
    call Input_Initialize(prob, geom)

    ! Check time
    call cpu_time(time)

    ! Modify the geometry
    call ModGeo_Modification(prob, geom)

    ! Set the edge cross-section
    call Section_Generation(prob, geom)

    ! Build the basepair
    call Basepair_Discretize(prob, geom, mesh)

    ! Solve the scaffold routing
    call Route_Generation(prob, geom, mesh, dna)

    ! Design the staple path and sequence
    call SeqDesign_Design(prob, geom, mesh, dna)

    ! Generate outputs
    call Output_Generation(prob, geom, mesh, dna)

    ! Print information
    call Print_Summary(prob, geom, mesh, dna)

    ! Verify the solution and check time
    !DEC$ IF DEFINED(_WIN32)
    if(iargc() == 0) then
        call Verify_Solution(mesh, dna)
        call Print_TimeConsuming(time)
    end if
    !DEC$ ENDIF

    ! Deallocate variables
    call Deallocate_Variables(geom, mesh, dna)

    if(p_redir /= 0) close(unit = p_redir)
end subroutine Main

! -----------------------------------------------------------------------------

! Print summary
subroutine Print_Summary(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    write(p_redir, "(a)")
    write(p_redir, "(a)"), " +==================================================+"
    write(p_redir, "(a)"), " | 8. Summary                                       |"
    write(p_redir, "(a)"), " +==================================================+"
    write(p_redir, "(a)")

    write(p_redir, "(a)"), " 8.1. Design parameters"
    write(p_redir, "(a)"), "   * Scaffold sequence              : "//trim(para_scaf_seq)
    write(p_redir, "(a)"), "   * Vertex design                  : "//trim(para_vertex_design)
    write(p_redir, "(a)"), "   * Vertex angle                   : "//trim(para_vertex_angle)
    write(p_redir, "(a)"), "   * Vertex clash                   : "//trim(para_vertex_crash)
    write(p_redir, "(a)"), "   * Constant edge length from mesh : "//trim(para_const_edge_mesh)
    write(p_redir, "(a)"), "   * Gap b/w two scaf xovers        : "//trim(adjustl(Int2Str(para_gap_xover_two_scaf)))
    write(p_redir, "(a)"), "   * Gap b/w xover(stap) and bound  : "//trim(adjustl(Int2Str(para_gap_xover_bound_stap)))
    write(p_redir, "(a)"), "   * Gap b/w stap and scaf xovers   : "//trim(adjustl(Int2Str(para_gap_xover_two)))
    write(p_redir, "(a)"), "   * Gap b/w xover and first nick   : "//trim(adjustl(Int2Str(para_gap_xover_nick1)))
    write(p_redir, "(a)"), "   * Gap b/w xover and nick         : "//trim(adjustl(Int2Str(para_gap_xover_nick)))
    write(p_redir, "(a)"), "   * Break rule for staples         : "//trim(para_cut_stap_method)
    write(p_redir, "(a)"), "   * Unpaired scaffold nucleotides  : "//trim(para_unpaired_scaf)
    write(p_redir, "(a)"), "   * # of changing for min staple   : "//trim(adjustl(Int2Str(prob.n_cng_min_stap)))
    write(p_redir, "(a)"), "   * # of changing for max staple   : "//trim(adjustl(Int2Str(prob.n_cng_max_stap)))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " 8.2. Target geometry"
    write(p_redir, "(a)"), "   * File name            : "//trim(prob.name_file)
    write(p_redir, "(a)"), "   * The number of faces  : "//trim(adjustl(Int2Str(geom.n_face)))
    write(p_redir, "(a)"), "   * The number of points : "//trim(adjustl(Int2Str(geom.n_iniP)))
    write(p_redir, "(a)"), "   * The number of edges  : "//trim(adjustl(Int2Str(geom.n_iniL)))
    write(p_redir, "(a)"), "   * Minimum edge-length  : "//trim(adjustl(Int2Str(prob.n_edge_len)))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " 8.3. Basepair"
    write(p_redir, "(a)"), "   * # of basepairs            : "//trim(adjustl(Int2Str(mesh.n_node)))
    write(p_redir, "(a)"), "   * # of mitered nts          : "//trim(adjustl(Int2Str(mesh.n_mitered)))//" ["//trim(adjustl(Dble2Str(dble(mesh.n_mitered)/dble(mesh.n_node)*100.0d0)))//" %]"
    write(p_redir, "(a)"), "   * Edge length [ min - max ] : ["//trim(adjustl(Int2Str(geom.min_edge_length)))//" - "//trim(adjustl(Int2Str(geom.max_edge_length)))//"]"
    write(p_redir, "(a)"), "   * Min # xovers [scaf, stap] : "//trim(adjustl(Int2Str(dna.min_xover_scaf+dna.min_xover_stap)))//" ["//trim(adjustl(Int2Str(dna.min_xover_scaf)))//", "//trim(adjustl(Int2Str(dna.min_xover_stap)))//"]"
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " 8.4. Scaffold"
    write(p_redir, "(a)"), "   * # of the scaffold  : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)"), "   * # of total nts     : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
    write(p_redir, "(a)"), "   * # of unpaired nts  : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_scaf)))
    write(p_redir, "(a)"), "   * # of double Xovers : "//trim(adjustl(Int2Str(dna.n_xover_scaf/2)))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " 8.5. Staple"
    write(p_redir, "(a)"), "   * # of staples         : "//trim(adjustl(Int2Str(dna.n_stap)))
    write(p_redir, "(a)"), "    @ with the 4nt dsDNA  : "//trim(adjustl(Int2Str(dna.n_4nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_4nt)/dble(dna.n_stap)*100.0d0)))//" %]"
    write(p_redir, "(a)"), "    @ with the 14nt dsDNA : "//trim(adjustl(Int2Str(dna.n_14nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_14nt)/dble(dna.n_stap)*100.0d0)))//" %]"
    write(p_redir, "(a)"), "   * # of nts             : "//trim(adjustl(Int2Str(dna.n_base_stap)))
    write(p_redir, "(a)"), "    @ in 4nt dsDNA        : "//trim(adjustl(Int2Str(dna.n_nt_4nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_nt_4nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"
    write(p_redir, "(a)"), "    @ in 14nt dsDNA       : "//trim(adjustl(Int2Str(dna.n_nt_14nt)))//" ["//trim(adjustl(Dble2Str(dble(dna.n_nt_14nt)/dble(dna.n_base_stap)*100.0d0)))//" %]"
    write(p_redir, "(a)"), "   * # of unpaired nts    : "//trim(adjustl(Int2Str(dna.n_nt_unpaired_stap)))
    write(p_redir, "(a)"), "   * # of total Xover     : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
    write(p_redir, "(a)"), "   * # of single-Xovers   : "//trim(adjustl(Int2Str(dna.n_sxover_stap)))
    write(p_redir, "(a)"), "   * # of double-Xovers   : "//trim(adjustl(Int2Str((dna.n_xover_stap-dna.n_sxover_stap)/2)))
    write(p_redir, "(a)"), "   * Length [min-max-ave] : ["//trim(adjustl(Int2Str(dna.len_min_stap)))//" - "//trim(adjustl(Int2Str(dna.len_max_stap)))//" - "//trim(adjustl(Dble2Str(dna.len_ave_stap)))//"]"
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " +=== completed =====================================+"
    write(p_redir, "(a)"), " | PERDIX generated output files.                    |"
    write(p_redir, "(a)"), " +===================================================+"
end subroutine Print_Summary

! -----------------------------------------------------------------------------

! Deallocate dynamic variables
subroutine Deallocate_Variables(geom, mesh, dna)
    type(GeomType), intent(inout) :: geom
    type(MeshType), intent(inout) :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i

    ! Deallocate the geom array
    do i = 1, geom.n_face
        deallocate(geom.face(i).poi)
    end do

    if(allocated(geom.sec.id))   deallocate(geom.sec.id)
    if(allocated(geom.sec.posR)) deallocate(geom.sec.posR)
    if(allocated(geom.sec.posC)) deallocate(geom.sec.posC)
    if(allocated(geom.sec.conn)) deallocate(geom.sec.conn)
    if(allocated(geom.iniP))     deallocate(geom.iniP)
    if(allocated(geom.modP))     deallocate(geom.modP)
    if(allocated(geom.croP))     deallocate(geom.croP)
    if(allocated(geom.iniL))     deallocate(geom.iniL)
    if(allocated(geom.croL))     deallocate(geom.croL)
    if(allocated(geom.face))     deallocate(geom.face)

    ! Deallocate the bound array
    do i = 1, geom.n_junc
        deallocate(geom.junc(i).iniL)
        deallocate(geom.junc(i).modP)
        deallocate(geom.junc(i).croP)
        deallocate(geom.junc(i).node)
        deallocate(geom.junc(i).conn)
        deallocate(geom.junc(i).type_conn)
    end do

    if(allocated(geom.junc)) deallocate(geom.junc)

    ! Deallocate mesh arrays
    if(allocated(mesh.node)) deallocate(mesh.node)
    if(allocated(mesh.ele))  deallocate(mesh.ele)

    ! Deallocate dna arrays
    do i = 1, dna.n_strand
        if(allocated(dna.strand(i).base)) deallocate(dna.strand(i).base)
    end do

    if(allocated(dna.base_scaf))  deallocate(dna.base_scaf)
    if(allocated(dna.base_stap))  deallocate(dna.base_stap)
    if(allocated(dna.top))        deallocate(dna.top)
    if(allocated(dna.strand))     deallocate(dna.strand)
    if(allocated(dna.order_stap)) deallocate(dna.order_stap)
end subroutine Deallocate_Variables

! -----------------------------------------------------------------------------

! Print computational time
subroutine Print_TimeConsuming(time_start)
    real, intent(in) :: time_start

    real :: time, time_end

    call cpu_time(time_end)

    time = time_end - time_start

    write(0, "(a       )")
    write(0, "(a       )"), "   --------------------------------------------------"
    write(0, "(a$      )"), "   Time consuming : "
    write(0, "(f6.2, a$)"), time, " [sec], "
    write(0, "(f6.2, a )"), time/60.0d0, " [min]"
    write(0, "(a       )"), "   --------------------------------------------------"
    write(0, "(a       )")
end subroutine Print_TimeConsuming

! -----------------------------------------------------------------------------

! Verify solution
subroutine Verify_Solution(mesh, dna)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    double precision :: verify
    integer :: i

    verify = 0.0d0
    do i = 1, dna.n_top
        if( dna.top(i).xover == -1 .or. &
            dna.top(i).up    == -1 .or. &
            dna.top(i).dn    == -1 .or. &
            dna.top(i).node  == -1 ) then
            verify = verify + dna.top(i).pos(1)
            verify = verify + dna.top(i).pos(2)
            verify = verify + dna.top(i).pos(3)
            verify = verify + (dble(dna.top(i).strand) + dble(dna.top(i).address))/dble(dna.n_top)
        end if
        verify = verify + (dble(dna.top(i).strand) + dble(dna.top(i).address))/dble(dna.n_top)

        if(dna.top(i).seq == "A") verify = verify + 1.0d0
        if(dna.top(i).seq == "T") verify = verify + 2.0d0
        if(dna.top(i).seq == "C") verify = verify + 3.0d0
        if(dna.top(i).seq == "G") verify = verify + 4.0d0
    end do

    verify = verify + dble(dna.n_nt_14nt)/dble(dna.n_base_stap)
    verify = verify + dble(dna.n_nt_4nt) /dble(dna.n_base_stap)
    verify = verify + dble(dna.len_ave_stap + dna.n_14nt + dna.n_s14nt + dna.n_4nt + dna.n_only_4nt)
    verify = verify + dble(dna.n_nt_14nt + dna.n_nt_4nt + dna.n_tot_region + dna.n_tot_14nt + dna.n_tot_4nt)
    verify = verify + dble(dna.len_min_stap + dna.len_max_stap + dna.n_unpaired_scaf + dna.n_nt_unpaired_scaf)
    verify = verify + dble(dna.n_unpaired_stap + dna.n_nt_unpaired_stap)

    do i = 1, mesh.n_node
        verify = verify + mesh.node(i).pos(1)
        verify = verify + mesh.node(i).pos(2)
        verify = verify + mesh.node(i).pos(3)
    end do

    write(0, "(a)")
    write(0, "(a)"), "   [ Check the solution - Debug mode only ]"
    write(0, "(a)"), "    3.41395107452898E+04 - Ref.: 1-1-1"
    write(0, "(a)"), "    3.77170400414344E+04 - Ref.: 2-1-1"
    write(0, "(a)"), "    1.35568076299989E+05 - Ref.: 7-1-1"
    write(0, "(es24.14)"), verify
end subroutine Verify_Solution

! -----------------------------------------------------------------------------

end program PERDIX