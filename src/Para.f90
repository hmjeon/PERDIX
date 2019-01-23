!
! =============================================================================
!
! Module - Para
! Last Updated : 01/14/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX, which allows scientists to build and solve
! the sequence design of complex DNA nanostructures.
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
module Para

    implicit none

    integer, parameter :: p_redir  = 0          ! 0:screen, 11:log file
    logical, parameter :: p_detail = .false.    ! .true.: detail, .false.: simple

    ! Parameters in the 'seq.txt' file
    character(10) :: para_scaf_seq

    ! Design parameters for the vertex
    character(10)  :: para_vertex_design   = "mitered"  ! [mitered, flat] Vertex design
    character(10)  :: para_vertex_angle    = "opt"      ! [opt, max, ave, min], Vertex angle
    character(10)  :: para_vertex_crash    = "const"    ! [const, mod], Vertex crash
    character(10)  :: para_const_edge_mesh = "off"      ! [off, on, round], Const. edge length design
    character(10)  :: para_sticky_self     = "off"      ! [off, on], Sticky-end (10-11) for flat vertex
    character(10)  :: para_unpaired_scaf   = "on"       ! [on, off], Unpaired scaffold nts

    ! Design paramters for the B-from DNA
    double precision :: para_dist_pp       = 0.42d0     ! [0.42, 0.6], distance b/w two nts, nm
    double precision :: para_dist_bp       = 0.34d0     ! [0.34 ], Axial rise distance, nm
    double precision :: para_rad_helix     = 1.0d0      ! [1.0  ], Radius of the DNA helix, nm
    double precision :: para_gap_helix     = 0.25d0     ! [0.25 ], Gap b/w two helixes, nm
    double precision :: para_ang_minor     = 150.0d0    ! [150.0], Angle of minor groove, degree
    double precision :: para_ang_correct   = 0.0d0      ! [0.0  ], Factor for adjusting ori., degree
    integer          :: para_n_base_tn     = -1         ! [-1   ], # of nts in Tn, -1: by distance
    integer          :: para_start_bp_ID   = -1         ! [-1   ], Starting bp ID as a reference, -1: pre-defined

    ! Design paramters for the scaffold route
    character(10)    :: para_weight_edge   = "on"       ! [on, off], Assign the weight factor
    character(10)    :: para_method_MST    = "prim"     ! [prim, kruskal, greedy], Minimum spanning tree
    character(10)    :: para_method_sort   = "quick"    ! [none, quick, shell], Sorting for MST
    character(10)    :: para_adjacent_list = "off"      ! [off, on], Adjacent list for Prim or Kruskal
    character(10)    :: para_all_spanning  = "off"      ! [off, on], All possible spanning trees

    ! Parameter for sequence design
    character(10) :: para_cut_stap_method  = "max"      ! [max, mix, opt, min, mid], Staple break rule
    character(10) :: para_set_stap_sxover  = "off"      ! [off, on], Single Xover
    character(10) :: para_output_design    = "arrow"    ! [arrow, seq, strand], Graphical outputs
    character(10) :: para_set_xover_scaf   = "split"    ! [split, center], Pattern for scaffold xovers

    integer   :: para_gap_xover_two_scaf   = 3          ! [3 ], Min. gap b/w two scaffold Xovers
    integer   :: para_gap_xover_bound_scaf = 7          ! [7 ], Min. gap b/w scaffold Xover and end-edge
    integer   :: para_gap_xover_bound_stap = 6          ! [6 ], Min. gap b/w staple Xover and end-edge
    integer   :: para_gap_xover_two        = 6          ! [6 ], Min. gap b/w scaffold and staple Xovers
    integer   :: para_gap_xover_nick1      = 10         ! [10], Min. gap b/w Xover(scaf/stap)/Tn and first nick
    integer   :: para_gap_xover_nick       = 3          ! [3 ], Min. gap b/w Xover and nick

    integer   :: para_max_cut_scaf         = 0          ! [0, 7249], Scaffold break
    integer   :: para_min_cut_stap         = 20         ! [20], Min. # of nts for one staple
    integer   :: para_mid_cut_stap         = 40         ! [40], Opt. # of nts for one staple
    integer   :: para_max_cut_stap         = 60         ! [60], Max. # of nts for one staple
    integer   :: para_set_start_scaf       = 1          ! [1], Starting nt position of the scaffold

    ! Output parameter
    logical :: para_tecplot     = .false.   !  Tecplot output
    logical :: para_write_101   = .false.   !  geo format file
    logical :: para_write_102   = .true.    ! *01_target_geometry.bild
    logical :: para_write_103   = .false.   !  init_geo_face.dat
    logical :: para_write_104   = .false.   !  schlegel.bild
    logical :: para_write_301   = .false.   !  check_geo.bild
    logical :: para_write_302   = .true.    ! *02_geometry_local.bild
    logical :: para_write_303   = .true.    ! *03_separate_lines.bild
    logical :: para_write_401   = .false.   !  cro_geo.bild
    logical :: para_write_501   = .false.   !  cyl_ori1.bild
    logical :: para_write_502   = .true.    ! *05_cylinder_prior / 06_cylinder_final.bild
    logical :: para_write_503   = .false.   !  mesh.bild
    logical :: para_write_504   = .true.    ! *04_section_lines.bild
    logical :: para_write_505   = .true.    ! *TXT_Edge_Length.txt
    logical :: para_write_601_1 = .false.   !  route1_scaf.bild / route1_stap.bild
    logical :: para_write_601_2 = .false.   !  route2_scaf.bild / route2_stap.bild
    logical :: para_write_601_3 = .false.   !  route3_scaf.bild / route3_stap.bild
    logical :: para_write_601_4 = .false.   !  route4_scaf.bild / route4_stap.bild
    logical :: para_write_601_5 = .false.   !  route5_scaf.bild / route5_stap.bild
    logical :: para_write_606   = .true.    ! *07_spantree.bild
    logical :: para_write_607   = .true.    ! *08_crossovers.bild
    logical :: para_write_608   = .false.   !  orientation.bild
    logical :: para_write_609   = .false.   !  atom.bild
    logical :: para_write_610   = .false.   !  scaf_xover.txt
    logical :: para_write_701   = .true.    ! *TXT_Sequence.txt
    logical :: para_write_711   = .false.   !  sequence.csv
    logical :: para_write_702   = .true.    ! *09_atomic_model.bild
    logical :: para_write_703   = .true.    ! *10_routing_scaf.bild / 11_routing_stap.bild
    logical :: para_write_705   = .true.    ! *12_routing_all.bild
    logical :: para_write_706   = .false.   !  strand.bild / sequence.bild
    logical :: para_write_710   = .false.   !  design_edgeX
    logical :: para_write_801   = .false.   !  basepair.txt
    logical :: para_write_802   = .false.   !  base.txt
    logical :: para_write_803   = .true.    ! *CanDo.cndo
    logical :: para_write_804   = .false.   !  tecplot.dat
    logical :: para_write_805   = .false.   !  adina.in
    logical :: para_write_808   = .false.   !  seq_line.tx

    ! BILD option
    logical :: para_chimera_axis     = .false.  !  Plot with the axis
    logical :: para_chimera_102_info = .true.   ! *Plot with edge and point #
    logical :: para_chimera_301_info = .false.  !  Plot with edge and point #
    logical :: para_chimera_302_info = .true.   ! *Plot with edge and point #
    logical :: para_chimera_303_info = .true.   ! *Plot with edge and point #
    logical :: para_chimera_401_info = .false.  !  Plot with edge and point #
    logical :: para_chimera_502_ori  = .false.  !  Plot with helix z-directions
    logical :: para_chimera_503_mod  = .false.  !  Plot with modified edges
    logical :: para_chimera_504_info = .true.   ! *Plot with edge and point #
    logical :: para_chimera_601_dir  = .false.  !  Plot with strand directions
    logical :: para_chimera_609_cyl  = .false.  !  Plot with cylinders
    logical :: para_chimera_609_dir  = .false.  !  Plot with strand directions

contains

! -----------------------------------------------------------------------------

end module Para