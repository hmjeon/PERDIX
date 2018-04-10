!
! =============================================================================
!
! Module - Chimera
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
module Chimera

    use Ifport
    use Data_Prob

    use Math
    use Para

    implicit none

    public  Chimera_Figure_Route_Step
    private Chimera_Init_ChimType

! -----------------------------------------------------------------------------

    ! Chimera data type
    type :: ChimType

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
        character(15)    :: preset              ! Preset
        logical          :: depthcue            ! Depth cue
        double precision :: light_brightness    ! Brightness
        double precision :: light_contrast      ! Contrast
        double precision :: light_ratio         ! Ratio
    end type ChimType

contains

! -----------------------------------------------------------------------------

! Make Figures from route step
subroutine Chimera_Figure_Route_Step(prob, name_file, nnn)
    type(ProbType), intent(in) :: prob
    character(200), intent(in) :: name_file
    integer,        intent(in) :: nnn

    type(ChimType) :: chim
    character(200) :: path
    logical :: results

    ! Initialzie Chimera data
    call Chimera_Init_ChimType(chim)

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
        write(1002, "(a, 2i7, a )"), "runCommand('windowsize ", chim.size_x, chim.size_y, "')"
        write(1002, "(a         )"), "runCommand('set projection "//trim(chim.projection_mode)//"')"
        write(1002, "(a         )"), "runCommand('lighting mode " //trim(chim.light_mode)//"')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting brightness ", chim.light_brightness, "')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting contrast ", chim.light_contrast, "')"
        write(1002, "(a, f8.2, a)"), "runCommand('lighting ratio ", chim.light_ratio, "')"
        write(1002, "(a         )"), "runCommand('preset apply "//trim(chim.preset)//"')"

        if(chim.depthcue == .false.) write(1002, "(a)"), "runCommand('~set depthCue')"
        if(chim.depthcue == .true. ) write(1002, "(a)"), "runCommand('set depthCue')"
    end if

    write(1002, "(a         )"), "runCommand('window')"
    write(1002, "(a, f8.2, a)"), "runCommand('scale ", chim.scale, "')"

    ! Set background color
    if(para_fig_bgcolor == "white") then
        ! --------------------------------------------------
        ! White background
        ! --------------------------------------------------
        write(1002, "(a)"), "runCommand('set bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor white')"

        ! Set view
        if(chim.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(chim.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz1") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -45')"
            write(1002, "(a )"), "runCommand('turn z 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ1/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn z -35')"
            write(1002, "(a )"), "runCommand('turn y 45')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz2") then
            write(1002, "(a )"), "runCommand('turn x 60')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ2/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -60')"
        else if(chim.view == "xyz" .or. chim.view == "all" ) then
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
        ! ==================================================
        ! Black background
        ! ==================================================
        write(1002, "(a)"), "runCommand('unset bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor black')"

        ! Set view
        if(chim.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(chim.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz1") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -45')"
            write(1002, "(a )"), "runCommand('turn z 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ1/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn z -35')"
            write(1002, "(a )"), "runCommand('turn y 45')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz2") then
            write(1002, "(a )"), "runCommand('turn x 60')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ2/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -60')"
        else if(chim.view == "xyz" .or. chim.view == "all" ) then
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
        ! ==================================================
        ! White and black background (white)
        ! ==================================================
        write(1002, "(a)"), "runCommand('set bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor white')"

        ! Set view
        if(chim.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(chim.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz1") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -45')"
            write(1002, "(a )"), "runCommand('turn z 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ1/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn z -35')"
            write(1002, "(a )"), "runCommand('turn y 45')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz2") then
            write(1002, "(a )"), "runCommand('turn x 60')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_white_XYZ2/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -60')"
        else if(chim.view == "xyz" .or. chim.view == "all" ) then
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

        ! ==================================================
        ! White and black background (black)
        ! ==================================================
        write(1002, "(a)"), "runCommand('unset bgTransparency')"
        write(1002, "(a)"), "runCommand('set bgcolor black')"

        ! Set view
        if(chim.view == "xy") then
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XY/"//trim(name_file)//".tif supersample 3"//"')"
        else if(chim.view == "xz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "yz") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -90')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_YZ/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn y 90')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz1") then
            write(1002, "(a )"), "runCommand('turn x -90')"
            write(1002, "(a )"), "runCommand('turn y -45')"
            write(1002, "(a )"), "runCommand('turn z 35')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ1/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn z -35')"
            write(1002, "(a )"), "runCommand('turn y 45')"
            write(1002, "(a )"), "runCommand('turn x 90')"
        else if(chim.view == "xyz2") then
            write(1002, "(a )"), "runCommand('turn x 60')"
            write(1002, "(a )"), "runCommand('wait')"
            write(1002, "(a$)"), "runCommand('copy file "//trim(prob.path_work2)
            write(1002, "(a )"), "step_fig_black_XYZ2/"//trim(name_file)//".tif supersample 3"//"')"
            write(1002, "(a )"), "runCommand('turn x -60')"
        else if(chim.view == "xyz" .or. chim.view == "all" ) then
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

! -----------------------------------------------------------------------------

! Initialize ChimType
subroutine Chimera_Init_ChimType(chim)
    type(ChimType), intent(out) :: chim

    chim.size_x = 800
    chim.size_y = 800
    chim.zoom   = 1.0d0
    chim.scale  = 1.0d0

    chim.turn_x = 0.0d0
    chim.turn_y = 0.0d0
    chim.turn_z = 0.0d0

    chim.move_x = 0.0d0
    chim.move_y = 0.0d0
    chim.move_z = 0.0d0

    ! Set mode
    chim.view            = para_fig_view
    chim.projection_mode = "orthographic"   ! [perspective | orthographic]
    chim.light_mode      = "two-point"      ! [ambient | single | two-point | three-point]
    chim.preset          = "publication 3"  ! Preset
    chim.depthcue        = .true.           ! ~set depthCue

    ! Set light as default
    chim.light_brightness = 1.16d0      ! From 0 to 2
    chim.light_contrast   = 0.83d0      ! From 0 to 1
    chim.light_ratio      = 1.25d0      ! From 1
end subroutine Chimera_Init_ChimType

! -----------------------------------------------------------------------------

end module Chimera