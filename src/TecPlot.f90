!
! ---------------------------------------------------------------------------------------
!
!                                   Module - TecPlot
!
!                                                                    Updated : 2017/03/27
!
! Comments: This module is to write outputs by Tecplot. 
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2018. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module TecPlot

    use Ifport
    use Data_Prob

    implicit none

    public  TecPlot_Command_Orientation

    private Tecplot_Write_Layout

! ---------------------------------------------------------------------------------------

    ! TecPlot data type
    type :: TecPlotType
        character(10) :: type_viewpoint     ! XY, XZ, YZ, XYZ
        character(10) :: show_scatter       ! YES, NO
        character(10) :: show_vector1       ! YES, NO
        character(10) :: show_vector2       ! YES, NO
        character(10) :: show_vector3       ! YES, NO

        double precision :: size_scatter
        double precision :: Norm
        double precision :: scale
    end type TecPlotType

contains

! ---------------------------------------------------------------------------------------

! Print 3 orientation vectors with bp
subroutine TecPlot_Command_Orientation(prob)
    type(ProbType), intent(in) :: prob

    type(TecPlotType) :: TecPlot

    ! Set TecPlot environment
    TecPlot.type_viewpoint = "ZXY"      ! XY, ZY, ZXY
    TecPlot.show_scatter   = "YES"
    TecPlot.show_vector1   = "YES"
    TecPlot.show_vector2   = "YES"
    TecPlot.show_vector3   = "YES"
    TecPlot.size_scatter   = 0.5d0
    TecPlot.Norm    = 1.0d0
    TecPlot.scale          = 50.0d0

    ! Write layout file in TecPlot
    call Tecplot_Write_Layout(prob, TecPlot)
end subroutine TecPlot_Command_Orientation

! ---------------------------------------------------------------------------------------

! Write layout in TecPlot
subroutine Tecplot_Write_Layout(prob, TecPlot)
    type(ProbType),    intent(in) :: prob
    type(TecPlotType), intent(in) :: TecPlot

    character(200) :: path

    path = trim(prob.path_work1)//trim(prob.name_file)
    open(unit=2001, file=trim(path)//".lay", form="formatted")

    write(2001, "(a)"), "#!MC 1410"
    write(2001, "(a)"), "$!VarSet |LFDSFN1| = '"//trim(prob.name_file)//"_tecplot.dat"//"'"
    write(2001, "(a)"), "$!VarSet |LFDSVL1| = '"//'"'//"'x'"//'"'//' "'//"'y'"&
        //'"'//' "'//"'z'"//'"'//' "'//"'V1x'"//'"'&
        //' "'//"'V1y'"//'"'//' "'//"'V1z'"//'"'//"'"

    write(2001, "(a)"), "$!READDATASET  '|LFDSFN1|'"

    ! Setvalue command that sets the position, border, and background attributes for the active
    write(2001, "(a)"), "$!FRAMELAYOUT"
    write(2001, "(a)"), "   SHOWHEADER = NO"
    write(2001, "(a)"), "   SHOWBORDER = NO"
    write(2001, "(a)"), "   HEADERCOLOR = BLACK"
    write(2001, "(a)"), "   WIDTH = 8"              ! Value is in inches
    write(2001, "(a)"), "   HEIGHT = 8"             ! Value is in inches
    write(2001, "(a)"), "   XYPOS"                  ! Position of upper left corner of the frame in inches from left and top edge of the paper
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       X = 1.5"
    write(2001, "(a)"), "       Y = 0.25"
    write(2001, "(a)"), "   }"

    ! Setvalue command that assigns attributes for axis in a 3D frame
    write(2001, "(a)"), "$!THREEDAXIS"
    write(2001, "(a)"), "   XDETAIL"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOWAXIS = NO"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   YDETAIL"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOWAXIS = NO"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   ZDETAIL"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOWAXIS = NO"
    write(2001, "(a)"), "   }"

    ! Changes plot types between valid Tecplot 360 modes such as XYLine and Cartesian2D
    write(2001, "(a)"), "$!PLOTTYPE  = CARTESIAN3D"

    ! Set the name for the active frame
    write(2001, "(a)"), "$!FRAMENAME  = '3DNA orientation vectors'"

    ! Setvalue command that changes global attributes associated with 3D vector plots
    write(2001, "(a)"), "$!GLOBALTHREEDVECTOR"
    write(2001, "(a)"), "   UVAR = 4"
    write(2001, "(a)"), "   VVAR = 5"
    write(2001, "(a)"), "   WVAR = 6"
    write(2001, "(a, e)"), "   RELATIVELENGTH = ", TecPlot.Norm

    ! Setvalue command that assigns zone attributes for field plots
    write(2001, "(a)"), "$!FIELDMAP  [1]"
    write(2001, "(a)"), "   MESH"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   VECTOR"
    write(2001, "(a)"), "   {"

    if(TecPlot.show_vector1 == "YES") then
        write(2001, "(a)"), "       SHOW = YES"
    else
        write(2001, "(a)"), "       SHOW = NO"
    end if

    write(2001, "(a)"), "       ARROWHEADSTYLE = FILLED"
    write(2001, "(a)"), "       COLOR = RED"
    write(2001, "(a)"), "       PATTERNLENGTH = 0.01"
    write(2001, "(a)"), "       LINETHICKNESS = 0.02"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   SCATTER"
    write(2001, "(a)"), "   {"

    if(TecPlot.show_scatter == "YES") then
        write(2001, "(a)"), "       SHOW = YES"
    else
        write(2001, "(a)"), "       SHOW = NO"
    end if

    write(2001, "(a)"), "       SYMBOLSHAPE"
    write(2001, "(a)"), "       {"
    write(2001, "(a)"), "           GEOMSHAPE = SPHERE"
    write(2001, "(a)"), "       }"
    write(2001, "(a)"), "   COLOR = YELLOW"
    write(2001, "(a)"), "   FILLCOLOR = YELLOW"
    write(2001, "(a, e)"), "    FRAMESIZE = ", TecPlot.size_scatter
    write(2001, "(a)"), "   }"

    write(2001, "(a)"), "$!FIELDMAP  [2]"
    write(2001, "(a)"), "   MESH"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   VECTOR"
    write(2001, "(a)"), "   {"

    if(TecPlot.show_vector2 == "YES") then
        write(2001, "(a)"), "       SHOW = YES"
    else
        write(2001, "(a)"), "       SHOW = NO"
    end if

    write(2001, "(a)"), "       ARROWHEADSTYLE = FILLED"
    write(2001, "(a)"), "       COLOR = BLUE"
    write(2001, "(a)"), "       PATTERNLENGTH = 1.5"
    write(2001, "(a)"), "       LINETHICKNESS = 0.02"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   SCATTER"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "       SYMBOLSHAPE"
    write(2001, "(a)"), "       {"
    write(2001, "(a)"), "           GEOMSHAPE = SPHERE"
    write(2001, "(a)"), "       }"
    write(2001, "(a)"), "   COLOR = YELLOW"
    write(2001, "(a)"), "   FILLCOLOR = YELLOW"
    write(2001, "(a)"), "   FRAMESIZE = 0.8"
    write(2001, "(a)"), "   }"

    write(2001, "(a)"), "$!FIELDMAP  [3]"
    write(2001, "(a)"), "   MESH"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   VECTOR"
    write(2001, "(a)"), "   {"

    if(TecPlot.show_vector3 == "YES") then
        write(2001, "(a)"), "       SHOW = YES"
    else
        write(2001, "(a)"), "       SHOW = NO"
    end if

    write(2001, "(a)"), "       ARROWHEADSTYLE = FILLED"
    write(2001, "(a)"), "       COLOR = CUSTOM3"        ! Orange
    write(2001, "(a)"), "       PATTERNLENGTH = 1.5"
    write(2001, "(a)"), "       LINETHICKNESS = 0.02"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   SCATTER"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "       SYMBOLSHAPE"
    write(2001, "(a)"), "       {"
    write(2001, "(a)"), "           GEOMSHAPE = SPHERE"
    write(2001, "(a)"), "       }"
    write(2001, "(a)"), "   COLOR = YELLOW"
    write(2001, "(a)"), "   FILLCOLOR = YELLOW"
    write(2001, "(a)"), "   FRAMESIZE = 0.8"
    write(2001, "(a)"), "   }"

    write(2001, "(a)"), "$!FIELDMAP  [4]"
    write(2001, "(a)"), "   MESH"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       COLOR = CUSTOM17"
    write(2001, "(a)"), "       LINETHICKNESS = 0.6"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   VECTOR"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "       ARROWHEADSTYLE = FILLED"
    write(2001, "(a)"), "       COLOR = BLUE"
    write(2001, "(a)"), "       PATTERNLENGTH = 1.5"
    write(2001, "(a)"), "       LINETHICKNESS = 0.02"
    write(2001, "(a)"), "   }"
    write(2001, "(a)"), "   SCATTER"
    write(2001, "(a)"), "   {"
    write(2001, "(a)"), "       SHOW = NO"
    write(2001, "(a)"), "       SYMBOLSHAPE"
    write(2001, "(a)"), "       {"
    write(2001, "(a)"), "           GEOMSHAPE = SPHERE"
    write(2001, "(a)"), "       }"
    write(2001, "(a)"), "   COLOR = GREEN"
    write(2001, "(a)"), "   FILLCOLOR = GREEN"
    write(2001, "(a)"), "   FRAMESIZE = 0.8"
    write(2001, "(a)"), "   }"

    ! Setvalue command that changes global attributes associated with the 3D view
    write(2001, "(a)"), "$!THREEDVIEW"

    if(TecPlot.type_viewpoint == "XY") then

        write(2001, "(a)"), "   PSIANGLE = 0.0"
        write(2001, "(a)"), "   THETAANGLE = 0.0"
        write(2001, "(a)"), "   ALPHAANGLE = 0.0"
        write(2001, "(a)"), "   VIEWERPOSITION"
        write(2001, "(a)"), "   {"
        write(2001, "(a)"), "       X = 0.0"
        write(2001, "(a)"), "       Y = 0.0"
        write(2001, "(a)"), "       Z = 200.0"
        write(2001, "(a)"), "   }"

    else if(TecPlot.type_viewpoint == "ZY") then

        write(2001, "(a)"), "   PSIANGLE = 90.0"
        write(2001, "(a)"), "   THETAANGLE = -90.0"
        write(2001, "(a)"), "   ALPHAANGLE = 0.0"
        write(2001, "(a)"), "   VIEWERPOSITION"
        write(2001, "(a)"), "   {"
        write(2001, "(a)"), "       X = 200.0"
        write(2001, "(a)"), "       Y = 0.0"
        write(2001, "(a)"), "       Z = 0.0"
        write(2001, "(a)"), "   }"

    else if(TecPlot.type_viewpoint == "ZXY") then

        write(2001, "(a)"), "   PSIANGLE = 48.36"
        write(2001, "(a)"), "   THETAANGLE = 117.24"
        write(2001, "(a)"), "   ALPHAANGLE = -108.88"
        write(2001, "(a)"), "   VIEWERPOSITION"
        write(2001, "(a)"), "   {"
        write(2001, "(a)"), "       X = -132.893"
        write(2001, "(a)"), "       Y = 68.404"
        write(2001, "(a)"), "       Z = 132.893"
        write(2001, "(a)"), "   }"

    end if

    write(2001, "(a, e)"), "   VIEWWIDTH = ", TecPlot.scale     ! Zoom in and out

    ! Setvalue command that turns field plot layers on or off, or sets the 2D draw order
    write(2001, "(a)"), "$!FIELDLAYERS"
    write(2001, "(a)"), "   SHOWVECTOR = YES"
    write(2001, "(a)"), "   SHOWSCATTER = YES"
    write(2001, "(a)"), "   SHOWEDGE = NO"
    write(2001, "(a)"), "   SHOWMESH = YES"

    close(unit=2001)
end subroutine Tecplot_Write_Layout

! ---------------------------------------------------------------------------------------

end module TecPlot