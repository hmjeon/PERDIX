!
! ---------------------------------------------------------------------------------------
!
!                                    Module for Para
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2016/02/03
!                                                            Last  modified  : 2016/11/21
!
! ---------------------------------------------------------------------------------------
!
module Para

    implicit none

    ! Geometric paprameters
    double precision, parameter :: para_init_scale       = 20.0d0       ! Initial geometric scale
    integer         , parameter :: para_n_square_lattice = 0            ! # of pre-defined square lattice cross-section

    ! Program parameters
    character(10)  :: para_external                     ! [off, on], External parameter loading flag
    character(10)  :: para_preset          = "on"       ! [on, off], Preset parameter defined in pre-defined examples
    character(10)  :: para_output_Tecplot  = "off"      ! [off, on], Output files for Tecplot(http://www.tecplot.com/) to draw vector image
    character(10)  :: para_cmd_Tecplot     = "off"      ! [off, on], Command file to run TecPlot automatically
    character(10)  :: para_cmd_Chimera     = "off"      ! [off, on], Command file to run UCSF Chimera(https://www.cgl.ucsf.edu/chimera/) automatically
    character(10)  :: para_fig_output      = "off"      ! [off, on], Automatic figure generation from outputs
    character(10)  :: para_fig_route_step  = "off"      ! [off, on], Automatic figure generation from route steps
    character(10)  :: para_fig_bgcolor     = "black"    ! [black, white, all], Background color for figures from UCSF Chimera
    character(10)  :: para_fig_view        = "preset"   ! [preset, xy, xz, xyz, all], Viewpoint for figures from UCSF Chimera
    integer        :: para_n_route_step    = 5          ! [5], The number of steps in routing progress
    integer        :: para_type_cndo       = 2          ! [2, 1], CanDo file option, 1 : original format, 2 : updated format
    character(200) :: para_path_Chimera    = &          ! UCSF Chimera program path
        "C:\Program Files\Chimera 1.10.2\bin\chimera.exe"

    ! Parameters for junction modification
    character(10)  :: para_junc_ang        = "opt"      ! [opt, max, ave, min], Junction gap modification for different arm angle
    character(10)  :: para_const_edge_mesh = "off"      ! [off, on], Constant edge length from polyhedra mesh
    character(10)  :: para_sticky_self     = "off"      ! [off, on], Sticky-end for self connection on henycomb cross-section
    character(10)  :: para_unpaired_scaf   = "on"       ! [on, off], Unpaired scaffold nucleotides
    character(10)  :: para_unpaired_square = "on"       ! [on, off], Unpaired scaffold and staple nucleotides at the vertex when square lattice
    character(10)  :: para_vertex_modify   = "const"    ! [const, mod], Vertex modification to avoid clash
    character(10)  :: para_vertex_design   = "flat"     ! [flat, beveled], Vertex design

    ! Paramters for B-from DNA generation
    double precision :: para_dist_pp       = 0.42d0     ! [0.42, 0.6], distance between adjacent phosphate groups, nm
    double precision :: para_dist_bp       = 0.34d0     ! [0.34 ], Axial rise distance, nm
    double precision :: para_rad_helix     = 1.0d0      ! [1.0  ], The radius of the DNA helix, nm
    double precision :: para_gap_helix     = 0.25d0     ! [0.25 ], The Gap between two helixes, nm
    double precision :: para_ang_minor     = 150.0d0    ! [150.0], An angle of minor groove, degree
    double precision :: para_ang_correct   = 0.0d0      ! [0.0  ], Correction factor to adjust orientation, degree
    integer          :: para_n_base_tn     = -1         ! [-1   ], The number of nucleotides in poly T loop, -1 : depending on distance
    integer          :: para_start_bp_ID   = -1         ! [-1   ], Starting base pair ID for the reference, -1 : pre-defined starting BP

    ! Paramters for scaffold route
    character(10)    :: para_weight_edge   = "on"       ! [on, off], Assign weight factor into edges of dual graph
    character(10)    :: para_method_MST    = "prim"     ! [prim, kruskal, greedy], Minimum spanning tree algorithm
    character(10)    :: para_method_sort   = "quick"    ! [none, quick, shell], Sorting algorithm to find MST for Prim or Kruskal
    character(10)    :: para_adjacent_list = "off"      ! [off, on], Output for adjacent list for Prim or Kruskal
    character(10)    :: para_all_spanning  = "on"       ! [on, off], All possible spanning trees when # of edges is less than 12 for Prim or Kruskal

    ! Parameter for sequence design
    character(10) :: para_cut_stap_method  = "opt"      ! [max, mix, opt, min, mid], Cutting method to make short staple strand, opt - 14nt seeds
    character(10) :: para_set_stap_sxover  = "off"      ! [off, on], To make non-circular staple by single crossover (when para_set_stap_sxover is "on")
    character(10) :: para_output_design    = "arrow"    ! [arrow, seq, strand], Graphical output type for sequence design
    character(10) :: para_set_xover_scaf   = "split"    ! [split, center], Setting possible scaffold strand

    integer   :: para_gap_xover_two_scaf   = 3          ! [3 ], The minimum gap between two scaffold crossovers
    integer   :: para_gap_xover_bound_scaf = 7          ! [7 ], The mimimum gap between scaffold crossover and vertex boundary
    integer   :: para_gap_xover_bound_stap = 6          ! [6 ], The mimimum gap between staple crossover and vertex boundary
    integer   :: para_gap_xover_two        = 6          ! [6 ], The minimum gap between scaffold and staple crossovers
    integer   :: para_gap_xover_nick1      = 10         ! [10], The minimum gap between xover(scaf/stap)/Tn and first nick
    integer   :: para_gap_xover_nick       = 3          ! [3 ], The minimum gap between xover and nick, if staple length exceeds 60, redesign with num - 1

    integer   :: para_max_cut_scaf         = 0          ! [0, -1], The maximum number of nucleotides for one scaffold strand, 0 : not cutting, -1 : cutting over 7249
    integer   :: para_min_cut_stap         = 20         ! [20], The minimum number of nucleotides for one staple strand
    integer   :: para_mid_cut_stap         = 40         ! [40], The optimal number of nucleotides for one staple strand
    integer   :: para_max_cut_stap         = 60         ! [60], The maximum number of nucleotides for one staple strand
    integer   :: para_set_seq_scaf         = 0          ! [0, 1, 2], Scaffold sequence, 0 - M13mp18(7249nt), 1 - import sequence from seq.txt, 2 - random
    integer   :: para_set_start_scaf       = 7217       ! [7217, 4141], Starting nucleotide position of scaffold strand

    ! UCSF Chimera output control
    logical :: para_write_101   = .false.       !  GEO file,                                Input_Write_GEO_File,             ".geo"
    logical :: para_write_102   = .true.        ! *Initial geometry,                        Input_Chimera_Init_Geometry,      "init_geo.bild"
    logical :: para_write_103   = .false.       !  Faced initial geometry,                  Input_Tecplot_Init_Geometry,      "init_geo_face.bild"
    logical :: para_write_104   = .false.       !  Schlegel diagram,                        Input_Chimera_Schlegel_Diagram,   "_schlegel.bild"
    logical :: para_write_301   = .true.        !  Initial geometry with face orientation,  ModGeo_Chimera_Check_Geometry,    "_check_geo.bild"
    logical :: para_write_302   = .true.        !  Initial geometry with local vector,      ModGeo_Chimera_Init_Geometry_L,   "_init_geo_local.bild"
    logical :: para_write_303   = .true.        ! *Modified geometry seperated from vertex, ModGeo_Chimera_Mod_Geometry,      "_mod_geo.bild"
    logical :: para_write_401   = .true.        ! *Cross-sectional geometry,                Section_Chimera_Cross_Geometry,   "_cross_geo.bild"
    logical :: para_write_501   = .true.       !  Cylindrical model with orientation,      Basepair_Chimera_Cylinder_Ori,    "_cyl_ori1.bild"
    logical :: para_write_502   = .true.        ! *Cylindrical model,                       Basepair_Chimera_Cylinder,        "_cyl1.bild", "_cyl2.bild"
    logical :: para_write_503   = .true.       !  Basepair model,                          Basepair_Chimera_Mesh,            "_mesh.bild"
    logical :: para_write_504   = .true.        ! *Cross-sectional geometry,                Basepair_Chimera_Cross_Geometry,  "_cross_geo_mod.bild"
    logical :: para_write_505   = .true.        ! *Txt file on edge length,                 Basepair_Write_Edge_Length,       "_edge_length.txt"
    logical :: para_write_601_1 = .true.       !  Route 1, seperated edges,                Route_Chimera_Route, step 1,      "_route1_scaf.bild", "_route1_stap.bild"
    logical :: para_write_601_2 = .true.       !  Route 2, contruction closed loop,        Route_Chimera_Route, step 2,      "_route2_scaf.bild", "_route2_stap.bild"
    logical :: para_write_601_3 = .true.        !  Route 3, centered crossovers             Route_Chimera_Route, step 3,      "_route3_scaf.bild", "_route3_stap.bild"
    logical :: para_write_601_4 = .true.        !  Route 4, modified centered crossovers,   Route_Chimera_Route, step 4,      "_route4_scaf.bild", "_route4_stap.bild"
    logical :: para_write_601_5 = .true.        ! *Route 5, scaffold route,                 Route_Chimera_Route, step 5,      "_route5_scaf.bild", "_route5_stap.bild"
    logical :: para_write_606   = .true.        ! *Sapnning tree for dual-graph,            Route_Graph_Chimera_Spanning_Tre, "_spantree.bild"
    logical :: para_write_607   = .true.        ! *Crossovers based on basepair model,      Route_Chimera_Crossovers,         "_crossovers.bild"
    logical :: para_write_608   = .true.        ! *3-orientation vectors,                   Route_Chimera_Orientation,        "_orientation.bild"
    logical :: para_write_609   = .true.       !  Atomic model without sequence design,    Route_Chimera_Atom,               "_atom.bild"
    logical :: para_write_610   = .true.        ! *Possible centered scaffold crossovers,   Route_Write_Centered_Scaf_Xover,  "_scaf_xover.txt"
    logical :: para_write_701   = .true.        ! *Txt on sequence design data,             SeqDesign_Write_Strand,           "strand.txt"
    logical :: para_write_711   = .true.       !  Csv file for sequence data,              SeqDesign_Write_Strand,           "sequence.csv"
    logical :: para_write_702   = .true.        ! *Atomic model with sequence design,       SeqDesign_Chimera_Atom,           "_atom_nick.bild"
    logical :: para_write_703   = .true.        ! *Route 6, strand route with nick,         SeqDesign_Chimera_Route,          "_route6_scaf.bild", "_route6_stap.bild"
    logical :: para_write_705   = .true.        ! *Sequence model,                          SeqDesign_Chimera_Sequence,       "_sequence_design.bild"
    logical :: para_write_706   = .true.       !  Atomic model bases on strands/sequence,  SeqDesign_Chimera_Strand,         "_strand.bild", "_sequence.bild"
    logical :: para_write_710   = .true.        ! *Edge-based sequence design,              SeqDesign_Write_Graphical_Output, "_design_edgeX"
    logical :: para_write_801   = .true.       !  Txt on basepair based data,              Output_Write_Basepair,            "_basepair.txt"
    logical :: para_write_802   = .true.       !  Txt on nucleotide based data,            Output_Write_Base,                "_base.txt"
    logical :: para_write_803   = .true.        ! *CanDo input file,                        Output_Write_CanDo,               ".cndo"
    logical :: para_write_804   = .true.       !  Tecplot input file,                      Output_Write_TecPlot,             "_tecplot.dat"
    logical :: para_write_805   = .true.       !  ADINA input file,                        Output_Write_ADINA,               "_adina.in"
    logical :: para_write_808   = .true.       !  Txt on sectional edges based sequence,   Output_Write_Sequence_CroL,       "_seq_line.txt"

    ! UCSF Chimera output option
    logical :: para_chimera_axis     = .true.   ! *Plot with axis at the ceneter of geometry (*.bild)
    logical :: para_chimera_102_info = .true.   ! *Plot with edge and point number (_init_geo.bild)
    logical :: para_chimera_301_info = .true.   !  Plot with edge and point number (_check_geo.bild)
    logical :: para_chimera_302_info = .true.   !  Plot with edge and point number (_init_geo_local.bild)
    logical :: para_chimera_303_info = .true.   ! *Plot with edge and point number (_mod_geo.bild)
    logical :: para_chimera_401_info = .true.   ! *Plot with edge and point number (_cross_geo.bild)
    logical :: para_chimera_502_ori  = .false.  !  Plot with helix z-direction (_line.bild / _node.bild)
    logical :: para_chimera_503_mod  = .false.  !  Plot with modified edges (_mesh.bild)
    logical :: para_chimera_504_info = .false.  !  Plot with edge and point number (_cross_geo_mod.bild)
    logical :: para_chimera_601_dir  = .false.  !  Plot with strand direction (_scaf.bild / _stap.bild)
    logical :: para_chimera_609_cyl  = .false.  !  Plot with cylinderical representation (_atom.bild)
    logical :: para_chimera_609_dir  = .false.  !  Plot with strand direction (_atom.bild)
contains

! ---------------------------------------------------------------------------------------

end module Para