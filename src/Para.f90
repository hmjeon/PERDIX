!
! =============================================================================
!
! Module - Para
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
module Para

    implicit none

    ! Geometric paprameters
    double precision, parameter :: para_init_scale = 20.0d0     ! Initial geometric scale

    ! Program parameters
    character(10)  :: para_platform        = "dev"      ! [dev, win, mac, linux], External parameter loading flag
    character(10)  :: para_preset          = "on"       ! [on, off], Preset parameter defined in pre-defined examples
    character(10)  :: para_output_Tecplot  = "off"      ! [off, on], Output files for Tecplot(http://www.tecplot.com/) to draw vector image
    character(10)  :: para_fig_view        = "xy"       ! [xy, xz, xyz1, xyz2, xyz, all], Viewpoint for figures from UCSF Chimera
    integer        :: para_type_cndo       = 2          ! [1, 2], CanDo file option, 1 : original format, 2 : updated format

    ! Parameters for junction modification
    character(10)  :: para_junc_ang        = "opt"      ! [opt, max, ave, min], Junction gap modification for different arm angle
    character(10)  :: para_const_edge_mesh = "off"      ! [off, on, round], Constant edge length from polyhedra mesh
    character(10)  :: para_sticky_self     = "off"      ! [off, on], Sticky-end for self connection on henycomb cross-section
    character(10)  :: para_unpaired_scaf   = "on"       ! [on, off], Unpaired scaffold nucleotides
    character(10)  :: para_vertex_modify   = "const"    ! [const, mod], Vertex modification to avoid clash
    character(10)  :: para_vertex_design   = "mitered"  ! [mitered, flat] Vertex design

    ! Paramters for B-from DNA generation
    double precision :: para_dist_pp       = 0.42d0     ! [0.42, 0.6], distance between adjacent phosphate groups, nm
    double precision :: para_dist_bp       = 0.34d0     ! [0.34 ], Axial rise distance, nm
    double precision :: para_rad_helix     = 1.0d0      ! [1.0  ], The radius of the DNA helix, nm
    double precision :: para_gap_helix     = 0.25d0     ! [0.25 ], The gap between two helixes, nm
    double precision :: para_ang_minor     = 150.0d0    ! [150.0], An angle of minor groove, degree
    double precision :: para_ang_correct   = 0.0d0      ! [0.0  ], Correction factor to adjust orientation, degree
    integer          :: para_n_base_tn     = -1         ! [-1   ], The number of nucleotides in poly T loop, -1 : depending on distance
    integer          :: para_start_bp_ID   = -1         ! [-1   ], Starting base pair ID for the reference, -1 : pre-defined starting BP

    ! Paramters for scaffold route
    character(10)    :: para_weight_edge   = "on"       ! [on, off], Assign weight factor into edges of dual graph
    character(10)    :: para_method_MST    = "prim"     ! [prim, kruskal, greedy], Minimum spanning tree algorithm
    character(10)    :: para_method_sort   = "quick"    ! [none, quick, shell], Sorting algorithm to find MST for Prim or Kruskal
    character(10)    :: para_adjacent_list = "off"      ! [off, on], Output for adjacent list for Prim or Kruskal
    character(10)    :: para_all_spanning  = "off"      ! [off, on], All possible spanning trees when # of edges is less than 12 for Prim or Kruskal

    ! Parameter for sequence design
    character(10) :: para_cut_stap_method  = "max"      ! [max, mix, opt, min, mid], Cutting method to make short staple strand, opt - 14nt seeds
    character(10) :: para_set_stap_sxover  = "off"      ! [off, on], To make non-circular staple by single crossover (when para_set_stap_sxover is "on")
    character(10) :: para_output_design    = "arrow"    ! [arrow, seq, strand], Graphical output type for sequence design
    character(10) :: para_set_xover_scaf   = "split"    ! [split, center], Setting possible scaffold strand

    integer   :: para_gap_xover_two_scaf   = 3          ! [3 ], The minimum gap between two scaffold crossovers
    integer   :: para_gap_xover_bound_scaf = 7          ! [7 ], The mimimum gap between scaffold crossover and vertex boundary
    integer   :: para_gap_xover_bound_stap = 6          ! [6 ], The mimimum gap between staple crossover and vertex boundary
    integer   :: para_gap_xover_two        = 6          ! [6 ], The minimum gap between scaffold and staple crossovers
    integer   :: para_gap_xover_nick1      = 10         ! [10], The minimum gap between xover(scaf/stap)/Tn and first nick
    integer   :: para_gap_xover_nick       = 3          ! [3 ], The minimum gap between xover and nick, if staple length exceeds 60, redesign with num - 1

    integer   :: para_max_cut_scaf         = 0          ! [0, 7249], Scaffold break - 0 : not breaking, num : breaking over num
    integer   :: para_min_cut_stap         = 20         ! [20], The minimum number of nucleotides for one staple strand
    integer   :: para_mid_cut_stap         = 40         ! [40], The optimal number of nucleotides for one staple strand
    integer   :: para_max_cut_stap         = 60         ! [60], The maximum number of nucleotides for one staple strand
    integer   :: para_set_seq_scaf         = 0          ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from env.txt, 2 - random
    integer   :: para_set_start_scaf       = 1          ! [1], Starting nucleotide position of scaffold strand

    ! UCSF Chimera output control
    logical :: para_write_101   = .false.       !  GEO file,                               ".geo"
    logical :: para_write_102   = .true.        ! *Initial geometry,                       "_01_target_geometry.bild"
    logical :: para_write_103   = .false.       !  Faced initial geometry,                 "init_geo_face.bild"
    logical :: para_write_104   = .false.       !  Schlegel diagram,                       "_schlegel.bild"
    logical :: para_write_301   = .true.       !  Initial geometry with face orientation, "_check_geo.bild"
    logical :: para_write_302   = .true.        ! *Initial geometry with local vector,     "_02_target_geometry.bild"
    logical :: para_write_303   = .true.        ! *Seperated lines from vertex,            "_03_seperated_lines.bild"
    logical :: para_write_401   = .false.       !  Cross-sectional geometry,               "_cro_geo.bild"
    logical :: para_write_501   = .false.       !  Cylindrical model with orientation,     "_cyl_ori1.bild"
    logical :: para_write_502   = .true.        ! *Cylindrical model,                      "05_cylindrical_model_1.bild", "06_cylindrical_model_2.bild"
    logical :: para_write_503   = .false.       !  Basepair model,                         "_mesh.bild"
    logical :: para_write_504   = .true.        ! *Multiple lines,                         "_04_doubled_lines.bild"
    logical :: para_write_505   = .true.        ! *Txt file on edge length,                "TXT_Edge_Length.txt"
    logical :: para_write_601_1 = .false.       !  Route 1, seperated edges,               "_route1_scaf.bild", "_route1_stap.bild"
    logical :: para_write_601_2 = .false.       !  Route 2, contruction closed loop,       "_route2_scaf.bild", "_route2_stap.bild"
    logical :: para_write_601_3 = .false.       !  Route 3, centered crossovers            "_route3_scaf.bild", "_route3_stap.bild"
    logical :: para_write_601_4 = .false.       !  Route 4, modified centered crossovers,  "_route4_scaf.bild", "_route4_stap.bild"
    logical :: para_write_601_5 = .false.       !  Route 5, scaffold route,                "_route5_scaf.bild", "_route5_stap.bild"
    logical :: para_write_606   = .true.        ! *Sapnning tree for dual-graph,           "_07_spantree.bild"
    logical :: para_write_607   = .true.        ! *Crossovers based on basepair model,     "_08_crossovers.bild"
    logical :: para_write_608   = .false.       !  3-orientation vectors,                  "_orientation.bild"
    logical :: para_write_609   = .false.       !  Atomic model without sequence design,   "_atom.bild"
    logical :: para_write_610   = .false.       !  Possible centered scaffold crossovers,  "_scaf_xover.txt"
    logical :: para_write_701   = .true.        ! *Txt on sequence design data,            "TXT_Sequence.txt"
    logical :: para_write_711   = .false.       !  Csv file for sequence data,             "sequence.csv"
    logical :: para_write_702   = .true.        ! *Atomic model with sequence design,      "_09_atomic_model.bild"
    logical :: para_write_703   = .true.        ! *Route 6, strand route with nick,        "_10_routing_scaf.bild", "_11_routing_stap.bild"
    logical :: para_write_705   = .true.        ! *Sequence model,                         "_12_routing_all.bild"
    logical :: para_write_706   = .false.       !  Atomic model bases on strands/sequence, "_strand.bild", "_sequence.bild"
    logical :: para_write_710   = .false.       !  Edge-based sequence design,             "_design_edgeX"
    logical :: para_write_801   = .false.       !  Txt on basepair based data,             "_basepair.txt"
    logical :: para_write_802   = .false.       !  Txt on nucleotide based data,           "_base.txt"
    logical :: para_write_803   = .true.        ! *CanDo input file,                       "_16_cndo_foramt.cndo"
    logical :: para_write_804   = .false.       !  Tecplot input file,                     "_tecplot.dat"
    logical :: para_write_805   = .false.       !  ADINA input file,                       "_adina.in"
    logical :: para_write_808   = .false.       !  Txt on sectional edges based sequence,  "_seq_line.txt"

    ! UCSF Chimera output option
    logical :: para_chimera_axis     = .false.  !  Plot with axis at the ceneter of geometry
    logical :: para_chimera_102_info = .true.   ! *Plot with edge and point number
    logical :: para_chimera_301_info = .false.  !  Plot with edge and point number
    logical :: para_chimera_302_info = .true.   ! *Plot with edge and point number
    logical :: para_chimera_303_info = .true.   ! *Plot with edge and point number
    logical :: para_chimera_401_info = .false.  !  Plot with edge and point number
    logical :: para_chimera_502_ori  = .false.  !  Plot with helix z-direction
    logical :: para_chimera_503_mod  = .false.  !  Plot with modified edges
    logical :: para_chimera_504_info = .true.   ! *Plot with edge and point number
    logical :: para_chimera_601_dir  = .false.  !  Plot with strand direction
    logical :: para_chimera_609_cyl  = .false.  !  Plot with cylinderical representation
    logical :: para_chimera_609_dir  = .false.  !  Plot with strand direction

contains

! -----------------------------------------------------------------------------

end module Para