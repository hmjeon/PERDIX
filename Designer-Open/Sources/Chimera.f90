!
! ---------------------------------------------------------------------------------------
!
!                                    Module for Chimera
!
!                                             Programmed by Hyungmin Jun (hmjeon@mit.edu)
!                                                   Massachusetts Institute of Technology
!                                                    Department of Biological Engineering
!                                         Laboratory for computational Biology & Biophics
!                                                            First programed : 2015/08/03
!                                                            Last  modified  : 2016/04/15
!
! ---------------------------------------------------------------------------------------
!
module Chimera

    use Ifport
    use Data_Prob

    use Math
    use Para

    implicit none

    public  Chimera_Command_Output
    public  Chimera_Figure_Output
    public  Chimera_Figure_Route_Step

    private Chimera_Init_ChimeraType
    private Chimera_Make_Python_Figure

! ---------------------------------------------------------------------------------------

    ! Chimera data type
    type :: ChimeraType

        logical :: run_init_geo         !  1. Input_Chimera_Init_Geometry,     "_init_geo.bild"
        logical :: run_check_geo        !  2. ModGeo_Chimera_Check_Geometry,   "_check_geo.bild"
        logical :: run_init_geo_local   !  3. ModGeo_Chimera_Init_Geometry_L,  "_init_geo_local.bild"
        logical :: run_mod_geo          !  4. ModGeo_Chimera_Mod_Geometry,     "_mod_geo.bild"
        logical :: run_cross_geo        !  5. Section_Chimera_Cross_Geometry,  "_cross_geo.bild"
        logical :: run_cylinder1_ori    !  6. Basepair_Chimera_Cylinder_Ori,   "_cyl1_ori.bild"
        logical :: run_cylinder1        !  7. Basepair_Chimera_Cylinder,       "_cyl1.bild"
        logical :: run_cylinder2        !  8. Basepair_Chimera_Cylinder,       "_cyl2.bild"
        logical :: run_mesh             !  9. Basepair_Chimera_Mesh,           "_mesh.bild"
        logical :: run_cross_geo_mod    ! 10. Basepair_Chimera_Cross_Geometry, "_cross_geo_mod.bild"
        logical :: run_route1_scaf      ! 11. Route_Chimera_Route, step 1,     "_route1_scaf.bild"
        logical :: run_route1_stap      ! 12. Route_Chimera_Route, step 1,     "_route1_stap.bild"
        logical :: run_route2_scaf      ! 13. Route_Chimera_Route, step 2,     "_route2_scaf.bild"
        logical :: run_route2_stap      ! 14. Route_Chimera_Route, step 2,     "_route2_stap.bild"
        logical :: run_route3_scaf      ! 15. Route_Chimera_Route, step 3,     "_route3_scaf.bild"
        logical :: run_route3_stap      ! 16. Route_Chimera_Route, step 3,     "_route3_stap.bild"
        logical :: run_route4_scaf      ! 17. Route_Chimera_Route, step 4,     "_route4_scaf.bild"
        logical :: run_route4_stap      ! 18. Route_Chimera_Route, step 4,     "_route4_stap.bild"
        logical :: run_route5_scaf      ! 19. Route_Chimera_Route, step 5,     "_route5_scaf.bild"
        logical :: run_route5_stap      ! 20. Route_Chimera_Route, step 5,     "_route5_stap.bild"
        logical :: run_crossovers       ! 21. Route_Chimera_Crossovers,        "_crossovers.bild"
        logical :: run_orientation      ! 22. Route_Chimera_Orientation,       "_orientation.bild"
        logical :: run_atom             ! 23. Route_Chimera_Atom,              "_atom.bild"
        logical :: run_atom_nick        ! 24. SeqDesign_Chimera_Atom,          "_atom_nick.bild"
        logical :: run_route6_scaf      ! 25. SeqDesign_Chimera_Route,         "_route6_scaf.bild"
        logical :: run_route6_stap      ! 26. SeqDesign_Chimera_Route,         "_route6_stap.bild"
        logical :: run_seq_design       ! 27. SeqDesign_Chimera_Sequence,      "_sequence_design.bild"
        logical :: run_strand           ! 28. SeqDesign_Chimera_Strand,        "_strand.bild",
        logical :: run_sequence         ! 29. SeqDesign_Chimera_Strand,        "_sequence.bild"

        double precision :: scale

        ! Control parameter
        character(10)    :: view                    ! XY, XZ, YZ, XYZ
        integer          :: size_x, size_y          ! Window size
        double precision :: turn_x, turn_y, turn_z  ! Turn
        double precision :: move_x, move_y, move_z  ! Move
        double precision :: zoom                    ! Zoom

        ! Light parameter
        character(15)    :: projection_mode     ! Projection mode
        character(15)    :: light_mode          ! Light mode
        double precision :: light_brightness    ! Brightness
        double precision :: light_contrast      ! Contrast
        double precision :: light_ratio         ! Ratio
    end type ChimeraType

contains
    
! ---------------------------------------------------------------------------------------

! Make command file to run outputs by Chimera
! Last updated on Friday 15 Apr 2016 by Hyungmin
subroutine Chimera_Command_Output(prob)
    type(ProbType), intent(in) :: prob

    type(ChimeraType) :: Chimera
    character(200) :: path
    logical :: results
    integer :: i

    ! Initialzie Chimera data
    call Chimera_Init_ChimeraType(Chimera)

    Chimera.scale = 1.0d0
    Chimera.view  = para_fig_view
    
    Chimera.run_init_geo        = para_write_102    ! 1.  Input_Chimera_Init_Geometry,     "_init_geo.bild"
    Chimera.run_check_geo       = para_write_301    ! 2.  ModGeo_Chimera_Check_Geometry,   "_check_geo.bild"
    Chimera.run_init_geo_local  = para_write_302    ! 3.  ModGeo_Chimera_Init_Geometry_L,  "_init_geo_local.bild"
    Chimera.run_mod_geo         = para_write_303    ! 4.  ModGeo_Chimera_Mod_Geometry,     "_mod_geo.bild"
    Chimera.run_cross_geo       = para_write_401    ! 5.  Section_Chimera_Cross_Geometry,  "_cross_geo.bild"
    Chimera.run_cylinder1_ori   = para_write_501    ! 6.  Basepair_Chimera_Cylinder_Ori,   "_cyl1_ori.bild"
    Chimera.run_cylinder1       = para_write_502    ! 7.  Basepair_Chimera_Cylinder,       "_cyl1.bild"
    Chimera.run_cylinder2       = para_write_502    ! 8.  Basepair_Chimera_Cylinder,       "_cyl2.bild"
    Chimera.run_mesh            = para_write_503    ! 9.  Basepair_Chimera_Mesh,           "_mesh.bild"
    Chimera.run_cross_geo_mod   = para_write_504    ! 10. Basepair_Chimera_Cross_Geometry, "_cross_geo_mod.bild"
    Chimera.run_route1_scaf     = para_write_601_1  ! 11. Route_Chimera_Route, step 1,     "_route1_scaf.bild"
    Chimera.run_route1_stap     = para_write_601_1  ! 12. Route_Chimera_Route, step 1,     "_route1_stap.bild"
    Chimera.run_route2_scaf     = para_write_601_2  ! 13. Route_Chimera_Route, step 2,     "_route2_scaf.bild"
    Chimera.run_route2_stap     = para_write_601_2  ! 14. Route_Chimera_Route, step 2,     "_route2_stap.bild"
    Chimera.run_route3_scaf     = para_write_601_3  ! 15. Route_Chimera_Route, step 3,     "_route3_scaf.bild"
    Chimera.run_route3_stap     = para_write_601_3  ! 16. Route_Chimera_Route, step 3,     "_route3_stap.bild"
    Chimera.run_route4_scaf     = para_write_601_4  ! 17. Route_Chimera_Route, step 4,     "_route4_scaf.bild"
    Chimera.run_route4_stap     = para_write_601_4  ! 18. Route_Chimera_Route, step 4,     "_route4_stap.bild"
    Chimera.run_route5_scaf     = para_write_601_5  ! 19. Route_Chimera_Route, step 5,     "_route5_scaf.bild"
    Chimera.run_route5_stap     = para_write_601_5  ! 20. Route_Chimera_Route, step 5,     "_route5_stap.bild"
    Chimera.run_crossovers      = para_write_607    ! 21. Route_Chimera_Crossovers,        "_crossovers.bild"
    Chimera.run_orientation     = para_write_608    ! 22. Route_Chimera_Orientation,       "_orientation.bild"
    Chimera.run_atom            = para_write_609    ! 23. Route_Chimera_Atom,              "_atom.bild"
    Chimera.run_atom_nick       = para_write_702    ! 24. SeqDesign_Chimera_Atom,          "_atom_nick.bild"
    Chimera.run_route6_scaf     = para_write_703    ! 25. SeqDesign_Chimera_Route,         "_route6_scaf.bild"
    Chimera.run_route6_stap     = para_write_703    ! 26. SeqDesign_Chimera_Route,         "_route6_stap.bild"
    Chimera.run_seq_design      = para_write_705    ! 27. SeqDesign_Chimera_Sequence,      "_sequence_design.bild"
    Chimera.run_strand          = para_write_706    ! 28. SeqDesign_Chimera_Strand,        "_strand.bild"
    Chimera.run_sequence        = para_write_706    ! 29. SeqDesign_Chimera_Strand,        "_sequence.bild"

    ! Make step directories for PNG files
    results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_py\")

    ! Python script for geometric data
    path = trim(prob.path_work1)//"output_py\"//trim(prob.name_file)

    open(unit=101, file=trim(path)//"_init_geo.py",        form="formatted")
    open(unit=102, file=trim(path)//"_check_geo.py",       form="formatted")
    open(unit=103, file=trim(path)//"_init_geo_local.py",  form="formatted")
    open(unit=104, file=trim(path)//"_mod_geo.py",         form="formatted")
    open(unit=105, file=trim(path)//"_cross_geo.py",       form="formatted")
    open(unit=106, file=trim(path)//"_cyl_ori1.py",        form="formatted")
    open(unit=107, file=trim(path)//"_cyl1.py",            form="formatted")
    open(unit=108, file=trim(path)//"_cyl2.py",            form="formatted")
    open(unit=109, file=trim(path)//"_mesh.py",            form="formatted")
    open(unit=110, file=trim(path)//"_cross_geo_mod.py",   form="formatted")
    open(unit=111, file=trim(path)//"_route1_scaf.py",     form="formatted")
    open(unit=112, file=trim(path)//"_route1_stap.py",     form="formatted")
    open(unit=113, file=trim(path)//"_route2_scaf.py",     form="formatted")
    open(unit=114, file=trim(path)//"_route2_stap.py",     form="formatted")
    open(unit=115, file=trim(path)//"_route3_scaf.py",     form="formatted")
    open(unit=116, file=trim(path)//"_route3_stap.py",     form="formatted")
    open(unit=117, file=trim(path)//"_route4_scaf.py",     form="formatted")
    open(unit=118, file=trim(path)//"_route4_stap.py",     form="formatted")
    open(unit=119, file=trim(path)//"_route5_scaf.py",     form="formatted")
    open(unit=120, file=trim(path)//"_route5_stap.py",     form="formatted")
    open(unit=121, file=trim(path)//"_crossovers.py",      form="formatted")
    open(unit=122, file=trim(path)//"_orientation.py",     form="formatted")
    open(unit=123, file=trim(path)//"_atom.py",            form="formatted")
    open(unit=124, file=trim(path)//"_atom_nick.py",       form="formatted")
    open(unit=125, file=trim(path)//"_route6_scaf.py",     form="formatted")
    open(unit=126, file=trim(path)//"_route6_stap.py",     form="formatted")
    open(unit=127, file=trim(path)//"_sequence_design.py", form="formatted")
    open(unit=128, file=trim(path)//"_strand.py",          form="formatted")
    open(unit=129, file=trim(path)//"_sequence.py",        form="formatted")

    ! Python script for geometric data
    do i = 101, 129

        ! http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/framecommand.html

        write(i, "(a)"), "from chimera import runCommand"

        path = trim(prob.name_file)
        if(i == 101) write(i, "(a)"), "runCommand('open "//trim(path)//"_init_geo.bild')"
        if(i == 102) write(i, "(a)"), "runCommand('open "//trim(path)//"_check_geo.bild')"
        if(i == 103) write(i, "(a)"), "runCommand('open "//trim(path)//"_init_geo_local.bild')"
        if(i == 104) write(i, "(a)"), "runCommand('open "//trim(path)//"_mod_geo.bild')"
        if(i == 105) write(i, "(a)"), "runCommand('open "//trim(path)//"_cross_geo.bild')"
        if(i == 106) write(i, "(a)"), "runCommand('open "//trim(path)//"_cyl_ori1.bild')"
        if(i == 107) write(i, "(a)"), "runCommand('open "//trim(path)//"_cyl1.bild')"
        if(i == 108) write(i, "(a)"), "runCommand('open "//trim(path)//"_cyl2.bild')"
        if(i == 109) write(i, "(a)"), "runCommand('open "//trim(path)//"_mesh.bild')"
        if(i == 110) write(i, "(a)"), "runCommand('open "//trim(path)//"_cross_geo_mod.bild')"
        if(i == 111) write(i, "(a)"), "runCommand('open "//trim(path)//"_route1_scaf.bild')"
        if(i == 112) write(i, "(a)"), "runCommand('open "//trim(path)//"_route1_stap.bild')"
        if(i == 113) write(i, "(a)"), "runCommand('open "//trim(path)//"_route2_scaf.bild')"
        if(i == 114) write(i, "(a)"), "runCommand('open "//trim(path)//"_route2_stap.bild')"
        if(i == 115) write(i, "(a)"), "runCommand('open "//trim(path)//"_route3_scaf.bild')"
        if(i == 116) write(i, "(a)"), "runCommand('open "//trim(path)//"_route3_stap.bild')"
        if(i == 117) write(i, "(a)"), "runCommand('open "//trim(path)//"_route4_scaf.bild')"
        if(i == 118) write(i, "(a)"), "runCommand('open "//trim(path)//"_route4_stap.bild')"
        if(i == 119) write(i, "(a)"), "runCommand('open "//trim(path)//"_route5_scaf.bild')"
        if(i == 120) write(i, "(a)"), "runCommand('open "//trim(path)//"_route5_stap.bild')"
        if(i == 121) write(i, "(a)"), "runCommand('open "//trim(path)//"_crossovers.bild')"
        if(i == 122) write(i, "(a)"), "runCommand('open "//trim(path)//"_orientation.bild')"
        if(i == 123) write(i, "(a)"), "runCommand('open "//trim(path)//"_atom.bild')"
        if(i == 124) write(i, "(a)"), "runCommand('open "//trim(path)//"_atom_nick.bild')"
        if(i == 125) write(i, "(a)"), "runCommand('open "//trim(path)//"_route6_scaf.bild')"
        if(i == 126) write(i, "(a)"), "runCommand('open "//trim(path)//"_route6_stap.bild')"
        if(i == 127) write(i, "(a)"), "runCommand('open "//trim(path)//"_sequence_design.bild')"
        if(i == 128) write(i, "(a)"), "runCommand('open "//trim(path)//"_strand.bild')"
        if(i == 129) write(i, "(a)"), "runCommand('open "//trim(path)//"_sequence.bild')"

        write(i, "(a)"), "runCommand('windowsize "//trim(adjustl(Int2Str(Chimera.size_x)))//" "//trim(adjustl(Int2Str(Chimera.size_y)))//"')"
        write(i, "(a)"), "runCommand('set projection "//trim(Chimera.projection_mode)//"')"
        write(i, "(a)"), "runCommand('lighting mode " //trim(Chimera.light_mode)//"')"

        if(i == 127) then
            write(i, "(a)"), "runCommand('lighting brightness 0.95')"
            write(i, "(a)"), "runCommand('lighting contrast 0.0')"
            write(i, "(a)"), "runCommand('lighting ratio 1.0')"
        else
            write(i, "(a)"), "runCommand('lighting brightness "//trim(adjustl(Dble2Str(Chimera.light_brightness)))//"')"
            write(i, "(a)"), "runCommand('lighting contrast "//trim(adjustl(Dble2Str(Chimera.light_contrast)))//"')"
            write(i, "(a)"), "runCommand('lighting ratio "//trim(adjustl(Dble2Str(Chimera.light_ratio)))//"')"
        end if

        write(i, "(a)"), "runCommand('~set depthCue')"
        write(i, "(a)"), "runCommand('preset apply publication 3')"

        ! Set background color
        if(para_fig_bgcolor == "black") then
            write(i, "(a)"), "runCommand('unset bgTransparency')"
            write(i, "(a)"), "runCommand('set bgcolor black')"
        else
            write(i, "(a)"), "runCommand('set bgTransparency')"
            write(i, "(a)"), "runCommand('set bgcolor white')"
        end if
        write(i, "(a)"), "runCommand('window')"

        ! Set x and y movement
        if(i == 127) then
            write(i, "(a)"), "runCommand('scale 0.9')"
        else
            write(i, "(a)"), "runCommand('move x "//trim(adjustl(Dble2Str(prob.move_x)))//"')"
            write(i, "(a)"), "runCommand('move y "//trim(adjustl(Dble2Str(prob.move_y)))//"')"
            write(i, "(a)"), "runCommand('scale "//trim(adjustl(Dble2Str(prob.size)))//"')"
        end if

        ! Set view
        if(Chimera.view == "XZ") then
            write(i, "(a)"), "runCommand('turn x -90')"
        else if(Chimera.view == "YZ") then
            write(i, "(a)"), "runCommand('turn x -90')"
            write(i, "(a)"), "runCommand('turn y -90')"
        else if(Chimera.view == "XYZ") then
            write(i, "(a)"), "runCommand('turn x -90')"
            write(i, "(a)"), "runCommand('turn y -120')"
            write(i, "(a)"), "runCommand('turn x 35')"
        end if

        ! Close files
        close(unit=i)
    end do

    ! Make batch file to run based on Python script
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=1001, file=trim(path)//"_output_run.bat", form="formatted")

    write(1001, "(a )"), "@echo off"

    path = trim('"'//trim(prob.path_chimera)//'" --script')&
        //" "//"output_py\"//trim(prob.name_file)

    if(Chimera.run_init_geo        == .true.) write(1001, "(a)"), trim(path)//"_init_geo.py"
    if(Chimera.run_check_geo       == .true.) write(1001, "(a)"), trim(path)//"_check_geo.py"
    if(Chimera.run_init_geo_local  == .true.) write(1001, "(a)"), trim(path)//"_init_geo_local.py"
    if(Chimera.run_mod_geo         == .true.) write(1001, "(a)"), trim(path)//"_mod_geo.py"
    if(Chimera.run_cross_geo       == .true.) write(1001, "(a)"), trim(path)//"_cross_geo.py"
    if(Chimera.run_cylinder1_ori   == .true.) write(1001, "(a)"), trim(path)//"_cyl1_ori.py"
    if(Chimera.run_cylinder1       == .true.) write(1001, "(a)"), trim(path)//"_cyl1.py"
    if(Chimera.run_cylinder2       == .true.) write(1001, "(a)"), trim(path)//"_cyl2.py"
    if(Chimera.run_mesh            == .true.) write(1001, "(a)"), trim(path)//"_mesh.py"
    if(Chimera.run_cross_geo_mod   == .true.) write(1001, "(a)"), trim(path)//"_cross_geo_mod.py"
    if(Chimera.run_route1_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route1_scaf.py"
    if(Chimera.run_route1_stap     == .true.) write(1001, "(a)"), trim(path)//"_route1_stap.py"
    if(Chimera.run_route2_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route2_scaf.py"
    if(Chimera.run_route2_stap     == .true.) write(1001, "(a)"), trim(path)//"_route2_stap.py"
    if(Chimera.run_route3_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route3_scaf.py"
    if(Chimera.run_route3_stap     == .true.) write(1001, "(a)"), trim(path)//"_route3_stap.py"
    if(Chimera.run_route4_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route4_scaf.py"
    if(Chimera.run_route4_stap     == .true.) write(1001, "(a)"), trim(path)//"_route4_stap.py"
    if(Chimera.run_route5_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route5_scaf.py"
    if(Chimera.run_route5_stap     == .true.) write(1001, "(a)"), trim(path)//"_route5_stap.py"
    if(Chimera.run_crossovers      == .true.) write(1001, "(a)"), trim(path)//"_crossovers.py"
    if(Chimera.run_orientation     == .true.) write(1001, "(a)"), trim(path)//"_orientation.py"
    if(Chimera.run_atom            == .true.) write(1001, "(a)"), trim(path)//"_atom.py"
    if(Chimera.run_atom_nick       == .true.) write(1001, "(a)"), trim(path)//"_atom_nick.py"
    if(Chimera.run_route6_scaf     == .true.) write(1001, "(a)"), trim(path)//"_route6_scaf.py"
    if(Chimera.run_route6_stap     == .true.) write(1001, "(a)"), trim(path)//"_route6_stap.py"
    if(Chimera.run_seq_design      == .true.) write(1001, "(a)"), trim(path)//"_sequence_design.py"
    if(Chimera.run_strand          == .true.) write(1001, "(a)"), trim(path)//"_strand.py"
    if(Chimera.run_sequence        == .true.) write(1001, "(a)"), trim(path)//"_sequence.py"

    close(unit=1001)
end subroutine Chimera_Command_Output

! ---------------------------------------------------------------------------------------

! Make figures from outputs
! Last updated on Friday 15 Apr 2016 by Hyungmin
subroutine Chimera_Figure_Output(prob)
    type(ProbType), intent(in) :: prob

    type(ChimeraType) :: Chimera

    ! Initialzie Chimera data
    call Chimera_Init_ChimeraType(Chimera)

    Chimera.size_x = 800
    Chimera.size_y = 800
    Chimera.scale  = 1.0d0
    Chimera.view   = para_fig_view

    Chimera.run_init_geo        = para_write_102    ! 1.  Input_Chimera_Init_Geometry,     "_init_geo.bild"
    Chimera.run_check_geo       = para_write_301    ! 2.  ModGeo_Chimera_Check_Geometry,   "_check_geo.bild"
    Chimera.run_init_geo_local  = para_write_302    ! 3.  ModGeo_Chimera_Init_Geometry_L,  "_init_geo_local.bild"
    Chimera.run_mod_geo         = para_write_303    ! 4.  ModGeo_Chimera_Mod_Geometry,     "_mod_geo.bild"
    Chimera.run_cross_geo       = para_write_401    ! 5.  Section_Chimera_Cross_Geometry,  "_cross_geo.bild"
    Chimera.run_cylinder1_ori   = para_write_501    ! 6.  Basepair_Chimera_Cylinder_Ori,   "_cyl1_ori.bild"
    Chimera.run_cylinder1       = para_write_502    ! 7.  Basepair_Chimera_Cylinder,       "_cyl1.bild"
    Chimera.run_cylinder2       = para_write_502    ! 8.  Basepair_Chimera_Cylinder,       "_cyl2.bild"
    Chimera.run_mesh            = para_write_503    ! 9.  Basepair_Chimera_Mesh,           "_mesh.bild"
    Chimera.run_cross_geo_mod   = para_write_504    ! 10. Basepair_Chimera_Cross_Geometry, "_cross_geo_mod.bild"
    Chimera.run_route1_scaf     = para_write_601_1  ! 11. Route_Chimera_Route, step 1,     "_route1_scaf.bild"
    Chimera.run_route1_stap     = para_write_601_1  ! 12. Route_Chimera_Route, step 1,     "_route1_stap.bild"
    Chimera.run_route2_scaf     = para_write_601_2  ! 13. Route_Chimera_Route, step 2,     "_route2_scaf.bild"
    Chimera.run_route2_stap     = para_write_601_2  ! 14. Route_Chimera_Route, step 2,     "_route2_stap.bild"
    Chimera.run_route3_scaf     = para_write_601_3  ! 15. Route_Chimera_Route, step 3,     "_route3_scaf.bild"
    Chimera.run_route3_stap     = para_write_601_3  ! 16. Route_Chimera_Route, step 3,     "_route3_stap.bild"
    Chimera.run_route4_scaf     = para_write_601_4  ! 17. Route_Chimera_Route, step 4,     "_route4_scaf.bild"
    Chimera.run_route4_stap     = para_write_601_4  ! 18. Route_Chimera_Route, step 4,     "_route4_stap.bild"
    Chimera.run_route5_scaf     = para_write_601_5  ! 19. Route_Chimera_Route, step 5,     "_route5_scaf.bild"
    Chimera.run_route5_stap     = para_write_601_5  ! 20. Route_Chimera_Route, step 5,     "_route5_stap.bild"
    Chimera.run_crossovers      = para_write_607    ! 21. Route_Chimera_Crossovers,        "_crossovers.bild"
    Chimera.run_orientation     = para_write_608    ! 22. Route_Chimera_Orientation,       "_orientation.bild"
    Chimera.run_atom            = para_write_609    ! 23. Route_Chimera_Atom,              "_atom.bild"
    Chimera.run_atom_nick       = para_write_702    ! 24. SeqDesign_Chimera_Atom,          "_atom_nick.bild"
    Chimera.run_route6_scaf     = para_write_703    ! 25. SeqDesign_Chimera_Route,         "_route6_scaf.bild"
    Chimera.run_route6_stap     = para_write_703    ! 26. SeqDesign_Chimera_Route,         "_route6_stap.bild"
    Chimera.run_seq_design      = para_write_705    ! 27. SeqDesign_Chimera_Sequence,      "_sequence_design.bild"
    Chimera.run_strand          = para_write_706    ! 28. SeqDesign_Chimera_Strand,        "_strand.bild"
    Chimera.run_sequence        = para_write_706    ! 29. SeqDesign_Chimera_Strand,        "_sequence.bild"

    ! Make Python script to make figures
    call Chimera_Make_Python_Figure(prob, Chimera)
end subroutine Chimera_Figure_Output

! ---------------------------------------------------------------------------------------

! Make Figures from route step
! Last updated on Monday 14 Mar 2016 by Hyungmin
subroutine Chimera_Figure_Route_Step(prob, name_file, nnn)
    type(ProbType), intent(in) :: prob
    character(200), intent(in) :: name_file
    integer,        intent(in) :: nnn

    type(ChimeraType) :: Chimera
    character(200) :: path
    logical :: results

    ! Initialzie Chimera data
    call Chimera_Init_ChimeraType(Chimera)

    Chimera.scale = 1.0d0
    Chimera.view  = para_fig_view

    ! Batch file to be excuted
    if(nnn == 1) then
        path = trim(prob.path_work1)//trim(prob.name_file)
        open(unit=1002, file=trim(path)//"_rstep.py", form="formatted")
    end if

    ! Python script to make the figure
    if(nnn == 1) then
        write(1002, "(a)"), "from chimera import runCommand"
    end if

    write(1002, "(a$)"), "runCommand('open "//trim(prob.path_work2)
    write(1002, "(a )"), "step_bild/"//trim(name_file)//".bild')"

    if(nnn == 1) then
        write(1002, "(a, 2i7, a )"), "runCommand('windowsize ", Chimera.size_x, Chimera.size_y, "')"
        write(1002, "(a         )"), "runCommand('set projection "//trim(Chimera.projection_mode)//"')"
        write(1002, "(a         )"), "runCommand('lighting mode " //trim(Chimera.light_mode)//"')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting brightness ", Chimera.light_brightness, "')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting contrast ", Chimera.light_contrast, "')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting ratio ", Chimera.light_ratio, "')"
        write(1002, "(a         )"), "runCommand('~set depthCue')"
        write(1002, "(a         )"), "runCommand('preset apply publication 3')"
    end if

    write(1002, "(a         )"), "runCommand('window')"
    write(1002, "(a, f8.2, a)"), "runCommand('scale ", Chimera.scale, "')"

    ! Set background color
    if(para_fig_bgcolor == "white") then
        ! --------------------------------------------------
        ! White background
        ! --------------------------------------------------
        write(1002, "(a)"), "runCommand('set bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor white')"

        ! Set view
        if(Chimera.view == "XY" .or. Chimera.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(Chimera.view == "XZ" .or. Chimera.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "YZ" .or. Chimera.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "XYZ" .or. Chimera.view == "xyz" .or. &
                Chimera.view == "ALL" .or. Chimera.view == "all" ) then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -120')"
            write(1002, "(a )"), "runCommand('turn x 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -35')"
            write(1002, "(a )"), "runCommand('turn y 120')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        end if
    else if(para_fig_bgcolor == "black") then
        ! --------------------------------------------------
        ! Black background
        ! --------------------------------------------------
        write(1002, "(a)"), "runCommand('unset bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor black')"

        ! Set view
        if(Chimera.view == "XY" .or. Chimera.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(Chimera.view == "XZ" .or. Chimera.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "YZ" .or. Chimera.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "XYZ" .or. Chimera.view == "xyz" .or. &
                Chimera.view == "ALL" .or. Chimera.view == "all" ) then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -120')"
            write(1002, "(a )"), "runCommand('turn x 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -35')"
            write(1002, "(a )"), "runCommand('turn y 120')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        end if
    else if(para_fig_bgcolor == "all") then
        ! --------------------------------------------------
        ! White and black background (white)
        ! --------------------------------------------------
        write(1002, "(a)"), "runCommand('set bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor white')"

        ! Set view
        if(Chimera.view == "XY" .or. Chimera.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(Chimera.view == "XZ" .or. Chimera.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "YZ" .or. Chimera.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "XYZ" .or. Chimera.view == "xyz" .or. &
                Chimera.view == "ALL" .or. Chimera.view == "all" ) then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -120')"
            write(1002, "(a )"), "runCommand('turn x 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -35')"
            write(1002, "(a )"), "runCommand('turn y 120')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        end if

        ! --------------------------------------------------
        ! White and black background (black)
        ! --------------------------------------------------
        write(1002, "(a)"), "runCommand('unset bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor black')"

        ! Set view
        if(Chimera.view == "XY" .or. Chimera.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(Chimera.view == "XZ" .or. Chimera.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "YZ" .or. Chimera.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(Chimera.view == "XYZ" .or. Chimera.view == "xyz" .or. &
                Chimera.view == "ALL" .or. Chimera.view == "all" ) then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -120')"
            write(1002, "(a )"), "runCommand('turn x 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -35')"
            write(1002, "(a )"), "runCommand('turn y 120')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        end if
    end if
    write(1002, "(a)"), "runCommand('close all')"

    if(nnn == 2) then
        write(1002, "(a)"), "runCommand('stop yes')"
        close(unit=1002)
    end if
end subroutine Chimera_Figure_Route_Step

! ---------------------------------------------------------------------------------------

! Initialize ChimeraType
! Last updated on Wendesday 16 Mar 2016 by Hyungmin
subroutine Chimera_Init_ChimeraType(Chimera)
    type(ChimeraType), intent(inout) :: Chimera

    Chimera.size_x = 800
    Chimera.size_y = 800
    Chimera.zoom   = 1.0d0
    Chimera.scale  = 1.0d0

    Chimera.turn_x = 0.0d0
    Chimera.turn_y = 0.0d0
    Chimera.turn_z = 0.0d0

    Chimera.move_x = 0.0d0
    Chimera.move_y = 0.0d0
    Chimera.move_z = 0.0d0

    ! Set mode
    Chimera.projection_mode  = "perspective"    ! [perspective / orthographic]
    Chimera.light_mode       = "two-point"      ! [ambient / single / two-point / three-point]

    ! Set light
    Chimera.light_brightness = 1.16d0           ! From 0 to 2
    Chimera.light_contrast   = 0.83d0           ! From 0 to 1
    Chimera.light_ratio      = 1.25d0           ! From 1
end subroutine Chimera_Init_ChimeraType

! ---------------------------------------------------------------------------------------

! Make Python script to make figures
! Last updated on Wendesday 16 Mar 2016 by Hyungmin
subroutine Chimera_Make_Python_Figure(prob, chi)
    type(ProbType),    intent(in)    :: prob
    type(ChimeraType), intent(inout) :: chi

    character(200) :: path, des
    logical :: results
    integer :: i

    ! Make output directories for TIF images
    if(para_fig_bgcolor == "white") then
        ! ==================================================
        !
        ! White background
        !
        ! ==================================================
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_YZ\")
        else if(para_fig_view == "xyz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XYZ\")
        else if(para_fig_view == "all") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XYZ\")
        end if
    else if(para_fig_bgcolor == "black") then
        ! ==================================================
        !
        ! Black background
        !
        ! ==================================================
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_YZ\")
        else if(para_fig_view == "xyz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XYZ\")
        else if(para_fig_view == "all") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XYZ\")
        end if
    else if(para_fig_bgcolor == "all") then
        ! ==================================================
        !
        ! White and black background
        !
        ! ==================================================
        if(para_fig_view == "xy") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XY\")
        else if(para_fig_view == "xz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XZ\")
        else if(para_fig_view == "yz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_YZ\")
        else if(para_fig_view == "xyz") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XYZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XYZ\")
        else if(para_fig_view == "all") then
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_white_XYZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XY\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_YZ\")
            results = SYSTEMQQ("md "//trim(prob.path_work1)//"output_figure_black_XYZ\")
        end if
    end if

    ! Python script for geometric data
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=101, file=trim(path)//"_output.py", form="formatted")

    write(101, "(a)"), "from chimera import runCommand"
    do i = 1, 29

        if(     i ==  1 .and. chi.run_init_geo       == .true.) then; des = "_init_geo"
        else if(i ==  2 .and. chi.run_check_geo      == .true.) then; des = "_check_geo"
        else if(i ==  3 .and. chi.run_init_geo_local == .true.) then; des = "_init_geo_local"
        else if(i ==  4 .and. chi.run_mod_geo        == .true.) then; des = "_mod_geo"
        else if(i ==  5 .and. chi.run_cross_geo      == .true.) then; des = "_cross_geo"
        else if(i ==  6 .and. chi.run_cylinder1_ori  == .true.) then; des = "_cyl1_ori"
        else if(i ==  7 .and. chi.run_cylinder1      == .true.) then; des = "_cyl1"
        else if(i ==  8 .and. chi.run_cylinder2      == .true.) then; des = "_cyl2"
        else if(i ==  9 .and. chi.run_mesh           == .true.) then; des = "_mesh"
        else if(i == 10 .and. chi.run_cross_geo_mod  == .true.) then; des = "_cross_geo_mod"
        else if(i == 11 .and. chi.run_route1_scaf    == .true.) then; des = "_route1_scaf"
        else if(i == 12 .and. chi.run_route1_stap    == .true.) then; des = "_route1_stap"
        else if(i == 13 .and. chi.run_route2_scaf    == .true.) then; des = "_route2_scaf"
        else if(i == 14 .and. chi.run_route2_stap    == .true.) then; des = "_route2_stap"
        else if(i == 15 .and. chi.run_route3_scaf    == .true.) then; des = "_route3_scaf"
        else if(i == 16 .and. chi.run_route3_stap    == .true.) then; des = "_route3_stap"
        else if(i == 17 .and. chi.run_route4_scaf    == .true.) then; des = "_route4_scaf"
        else if(i == 18 .and. chi.run_route4_stap    == .true.) then; des = "_route4_stap"
        else if(i == 19 .and. chi.run_route5_scaf    == .true.) then; des = "_route5_scaf"
        else if(i == 20 .and. chi.run_route5_stap    == .true.) then; des = "_route5_stap"
        else if(i == 21 .and. chi.run_crossovers     == .true.) then; des = "_crossovers"
        else if(i == 22 .and. chi.run_orientation    == .true.) then; des = "_orientation"
        else if(i == 23 .and. chi.run_atom           == .true.) then; des = "_atom"
        else if(i == 24 .and. chi.run_atom_nick      == .true.) then; des = "_atom_nick"
        else if(i == 25 .and. chi.run_route6_scaf    == .true.) then; des = "_route6_scaf"
        else if(i == 26 .and. chi.run_route6_stap    == .true.) then; des = "_route6_stap"
        else if(i == 27 .and. chi.run_seq_design     == .true.) then; des = "_sequence_design"
        else if(i == 28 .and. chi.run_strand         == .true.) then; des = "_strand"
        else if(i == 29 .and. chi.run_sequence       == .true.) then; des = "_sequence"
        else
            cycle
        end if

        write(101, "(a)"), "runCommand('open "//trim(prob.path_work2)//trim(prob.name_file)//trim(des)//".bild')"
        write(101, "(a)"), "runCommand('windowsize "//trim(adjustl(Int2Str(chi.size_x)))//" "//trim(adjustl(Int2Str(chi.size_y)))//"')"
        write(101, "(a)"), "runCommand('set projection "//trim(chi.projection_mode)//"')"
        write(101, "(a)"), "runCommand('lighting mode " //trim(chi.light_mode)//"')"

        if(i == 27) then
            write(101, "(a)"), "runCommand('lighting brightness 0.85')"
            write(101, "(a)"), "runCommand('lighting contrast 0.0')"
            write(101, "(a)"), "runCommand('lighting ratio 1.0')"
        else
            write(101, "(a)"), "runCommand('lighting brightness "//trim(adjustl(Dble2Str(chi.light_brightness)))//"')"
            write(101, "(a)"), "runCommand('lighting contrast "//trim(adjustl(Dble2Str(chi.light_contrast)))//"')"
            write(101, "(a)"), "runCommand('lighting ratio "//trim(adjustl(Dble2Str(chi.light_ratio)))//"')"
        end if

        write(101, "(a)"), "runCommand('~set depthCue')"
        write(101, "(a)"), "runCommand('preset apply publication 3')"
        write(101, "(a)"), "runCommand('window')"

        ! Set x and y movement
        if(i == 27) then
            write(101, "(a)"), "runCommand('scale 0.9')"
        else
            write(101, "(a)"), "runCommand('move x "//trim(adjustl(Dble2Str(prob.move_x)))//"')"
            write(101, "(a)"), "runCommand('move y "//trim(adjustl(Dble2Str(prob.move_y)))//"')"
            write(101, "(a)"), "runCommand('scale "//trim(adjustl(Dble2Str(prob.size)))//"')"
        end if

        ! Set background color
        if(para_fig_bgcolor == "white") then
            ! ==================================================
            !
            ! White background
            !
            ! ==================================================
            write(101, "(a)"), "runCommand('set bgTransparency')"
            write(101, "(a)"), "runCommand('set bgcolor white')"

            ! Set view
            if(chi.view == "XY" .or. chi.view == "xy") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"
            else if(chi.view == "XZ" .or. chi.view == "xz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "YZ" .or. chi.view == "yz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn y 90')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "XYZ" .or. chi.view == "xyz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -120')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "ALL" .or. chi.view == "all") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -30')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            end if

        else if(para_fig_bgcolor == "black") then
            ! ==================================================
            !
            ! Black background
            !
            ! ==================================================
            write(101, "(a)"), "runCommand('unset bgTransparency')"
            write(101, "(a)"), "runCommand('set bgcolor black')"

            ! Set view
            if(chi.view == "XY" .or. chi.view == "xy") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"
            else if(chi.view == "XZ" .or. chi.view == "xz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "YZ" .or. chi.view == "yz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn y 90')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "XYZ" .or. chi.view == "xyz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -120')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "ALL" .or. chi.view == "all") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -30')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            end if

        else if(para_fig_bgcolor == "all") then
            ! ==================================================
            !
            ! White and black background (white)
            !
            ! ==================================================
            write(101, "(a)"), "runCommand('set bgTransparency')"
            write(101, "(a)"), "runCommand('set bgcolor white')"

            ! Set view
            if(chi.view == "XY" .or. chi.view == "xy") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"
            else if(chi.view == "XZ" .or. chi.view == "xz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "YZ" .or. chi.view == "yz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn y 90')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "XYZ" .or. chi.view == "xyz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -120')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "ALL" .or. chi.view == "all") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -30')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_white_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            end if

            ! ==================================================
            !
            ! White and black background (black)
            !
            ! ==================================================
            write(101, "(a)"), "runCommand('unset bgTransparency')"
            write(101, "(a)"), "runCommand('set bgcolor black')"

            ! Set view
            if(chi.view == "XY" .or. chi.view == "xy") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"
            else if(chi.view == "XZ" .or. chi.view == "xz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "YZ" .or. chi.view == "yz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn y 90')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "XYZ" .or. chi.view == "xyz") then
                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('turn y -120')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            else if(chi.view == "ALL" .or. chi.view == "all") then
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XY/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XY.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn x -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -90')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_YZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_YZ.tif supersample 3"//"')"

                write(101, "(a )"), "runCommand('turn y -30')"
                write(101, "(a )"), "runCommand('turn x 35')"
                write(101, "(a )"), "runCommand('wait')"
                write(101, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
                write(101, "(a$)"), "output_figure_black_XYZ/"//trim(prob.name_file)
                write(101, "(a )"), trim(des)//"_XYZ.tif supersample 3"//"')"
                write(101, "(a )"), "runCommand('turn x -35')"
                write(101, "(a )"), "runCommand('turn y 120')"
                write(101, "(a )"), "runCommand('turn x 90')"
            end if
        end if

        write(101, "(a)"), "runCommand('close all')"
    end do

    ! Close files
    write(101, "(a)"), "runCommand('stop yes')"
    close(unit=101)

    ! Make batch file to make TIF images
    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=1003, file=trim(path)//"_output.bat", form="formatted")

    path = trim('"'//trim(prob.path_chimera)//'" --script')//" "&
        //trim(prob.path_work1)//trim(prob.name_file)
    write(1003, "(a )"), "@echo off"
    write(1003, "(a)"), trim(path)//"_output.py"

    path = trim(prob.path_work1)//trim(prob.name_file)
    results = SYSTEMQQ("start cmd /C call "//trim(path)//"_output.bat")

    close (unit=1003)
end subroutine Chimera_Make_Python_Figure

! ---------------------------------------------------------------------------------------

end module Chimera