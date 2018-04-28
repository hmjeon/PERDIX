!
! =============================================================================
!
! Module - Para
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
    integer   :: para_set_seq_scaf         = 0          ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from seq.txt, 2 - random
    integer   :: para_set_start_scaf       = 1          ! [1], Starting nucleotide position of scaffold strand

    ! UCSF Chimera output control
    logical :: para_write_101   = .false.       !  GEO file,                                Input_Write_GEO_File,             ".geo"
    logical :: para_write_102   = .true.        ! *Initial geometry,                        Input_Chimera_Init_Geometry,      "_01_init_geo.bild"
    logical :: para_write_103   = .false.       !  Faced initial geometry,                  Input_Tecplot_Init_Geometry,      "init_geo_face.bild"
    logical :: para_write_104   = .false.       !  Schlegel diagram,                        Input_Chimera_Schlegel_Diagram,   "_schlegel.bild"
    logical :: para_write_301   = .false.       !  Initial geometry with face orientation,  ModGeo_Chimera_Check_Geometry,    "_check_geo.bild"
    logical :: para_write_302   = .true.        ! *Initial geometry with local vector,      ModGeo_Chimera_Init_Geometry_L,   "_02_init_geo_local.bild"
    logical :: para_write_303   = .true.        ! *Seperated lines from vertex,             ModGeo_Chimera_Mod_Geometry,      "_03_sep_line.bild"
    logical :: para_write_401   = .false.       !  Cross-sectional geometry,                Section_Chimera_Cross_Geometry,   "_cro_geo.bild"
    logical :: para_write_501   = .false.       !  Cylindrical model with orientation,      Basepair_Chimera_Cylinder_Ori,    "_cyl_ori1.bild"
    logical :: para_write_502   = .true.        ! *Cylindrical model,                       Basepair_Chimera_Cylinder,        "04_cylinder_1.bild", "05_cylinder_2.bild"
    logical :: para_write_503   = .false.       !  Basepair model,                          Basepair_Chimera_Mesh,            "_mesh.bild"
    logical :: para_write_504   = .true.        ! *Multiple lines,                          Basepair_Chimera_Cross_Geometry,  "_06_multi_line.bild"
    logical :: para_write_505   = .true.        ! *Txt file on edge length,                 Basepair_Write_Edge_Length,       "TXT_Edge_Length.txt"
    logical :: para_write_601_1 = .false.       !  Route 1, seperated edges,                Route_Chimera_Route, step 1,      "_route1_scaf.bild", "_route1_stap.bild"
    logical :: para_write_601_2 = .false.       !  Route 2, contruction closed loop,        Route_Chimera_Route, step 2,      "_route2_scaf.bild", "_route2_stap.bild"
    logical :: para_write_601_3 = .false.       !  Route 3, centered crossovers             Route_Chimera_Route, step 3,      "_route3_scaf.bild", "_route3_stap.bild"
    logical :: para_write_601_4 = .false.       !  Route 4, modified centered crossovers,   Route_Chimera_Route, step 4,      "_route4_scaf.bild", "_route4_stap.bild"
    logical :: para_write_601_5 = .false.       !  Route 5, scaffold route,                 Route_Chimera_Route, step 5,      "_route5_scaf.bild", "_route5_stap.bild"
    logical :: para_write_606   = .true.        ! *Sapnning tree for dual-graph,            Route_Graph_Chimera_Spanning_Tre, "_07_spantree.bild"
    logical :: para_write_607   = .true.        ! *Crossovers based on basepair model,      Route_Chimera_Crossovers,         "_08_xovers.bild"
    logical :: para_write_608   = .false.       !  3-orientation vectors,                   Route_Chimera_Orientation,        "_orientation.bild"
    logical :: para_write_609   = .false.       !  Atomic model without sequence design,    Route_Chimera_Atom,               "_atom.bild"
    logical :: para_write_610   = .false.       !  Possible centered scaffold crossovers,   Route_Write_Centered_Scaf_Xover,  "_scaf_xover.txt"
    logical :: para_write_701   = .true.        ! *Txt on sequence design data,             SeqDesign_Write_Outputs,          "TXT_Sequence.txt"
    logical :: para_write_711   = .false.       !  Csv file for sequence data,              SeqDesign_Write_Outputs,          "sequence.csv"
    logical :: para_write_702   = .true.        ! *Atomic model with sequence design,       SeqDesign_Chimera_Atom,           "_09_atomic_model.bild"
    logical :: para_write_703   = .true.        ! *Route 6, strand route with nick,         SeqDesign_Chimera_Route,          "_10_route_scaf.bild", "_11_route_stap.bild"
    logical :: para_write_705   = .true.        ! *Sequence model,                          SeqDesign_Chimera_Sequence,       "_12_route_all.bild"
    logical :: para_write_706   = .false.       !  Atomic model bases on strands/sequence,  SeqDesign_Chimera_Strand,         "_strand.bild", "_sequence.bild"
    logical :: para_write_710   = .false.       !  Edge-based sequence design,              SeqDesign_Write_Graphical_Output, "_design_edgeX"
    logical :: para_write_801   = .false.       !  Txt on basepair based data,              Output_Write_Basepair,            "_basepair.txt"
    logical :: para_write_802   = .false.       !  Txt on nucleotide based data,            Output_Write_Base,                "_base.txt"
    logical :: para_write_803   = .true.        ! *CanDo input file,                        Output_Write_CanDo,               "_16_cndo.cndo"
    logical :: para_write_804   = .false.       !  Tecplot input file,                      Output_Write_TecPlot,             "_tecplot.dat"
    logical :: para_write_805   = .false.       !  ADINA input file,                        Output_Write_ADINA,               "_adina.in"
    logical :: para_write_808   = .false.       !  Txt on sectional edges based sequence,   Output_Write_Sequence_CroL,       "_seq_line.txt"

    ! UCSF Chimera output option
    logical :: para_chimera_axis     = .false.  !  Plot with axis at the ceneter of geometry (*.bild)
    logical :: para_chimera_102_info = .true.   ! *Plot with edge and point number (_01_init_geo.bild)
    logical :: para_chimera_301_info = .false.  !  Plot with edge and point number (_check_geo.bild)
    logical :: para_chimera_302_info = .true.   ! *Plot with edge and point number (_02_init_geo_local.bild)
    logical :: para_chimera_303_info = .true.   ! *Plot with edge and point number (_03_sep_line.bild)
    logical :: para_chimera_401_info = .false.  !  Plot with edge and point number (_cro_geo.bild)
    logical :: para_chimera_502_ori  = .false.  !  Plot with helix z-direction (_line.bild / _node.bild)
    logical :: para_chimera_503_mod  = .false.  !  Plot with modified edges (_mesh.bild)
    logical :: para_chimera_504_info = .true.   ! *Plot with edge and point number (_06_multi_line.bild)
    logical :: para_chimera_601_dir  = .false.  !  Plot with strand direction (_scaf.bild / _stap.bild)
    logical :: para_chimera_609_cyl  = .false.  !  Plot with cylinderical representation (_atom.bild)
    logical :: para_chimera_609_dir  = .false.  !  Plot with strand direction (_atom.bild)

contains

! -----------------------------------------------------------------------------

end module Para