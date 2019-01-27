#
# Makefiles for PERDIX, written by Hyungmin Jun
#
#	-fast
#	Maximizes speed across the entire program.
#	Description:
#	This option maximizes speed across the entire program.
#
#	It sets the following options:
#	- On Linux* systems: -ipo, -O3, -no-prec-div, -static, and -xHost
#	- On OS X* systems: -ipo, -mdynamic-no-pic, -O3, -no-prec-div, and -xHost
#	- On Windows*  systems:  /O3,  /Qipo,  /Qprec-div-,  and /QxHost
#
#	NOTE: Option fast sets some aggressive optimizations that may not be
#   appropriate for all applications. The resulting executable may not run on
#   processor types different from the one on which you compile. You should
#   make sure that you understand the individual optimization options
#   that are enabled by option fast.

objects = Data_DNA.o Data_Geom.o Data_Mesh.o Data_Prob.o \
	Para.o List.o Math.o Mani.o SpanTree.o Importer.o \
	Exam_PERDIX.o Exam_PERDIX.o \
	Section.o Input.o ModGeo.o Basepair.o Route.o SeqDesign.o Output.o PERDIX.o

f90comp = ifort

# Makefile
PERDIX: $(objects)
	$(f90comp) -o PERDIX $(objects)

Data_DNA.mod: Data_DNA.o Data_DNA.f90
	$(f90comp) -c Data_DNA.f90
Data_DNA.o: Data_DNA.f90
	$(f90comp) -c Data_DNA.f90

Data_Geom.mod: Data_Geom.o Data_Geom.f90
	$(f90comp) -c Data_Geom.f90
Data_Geom.o: Data_Geom.f90
	$(f90comp) -c Data_Geom.f90

Data_Mesh.mod: Data_Mesh.o Data_Mesh.f90
	$(f90comp) -c Data_Mesh.f90
Data_Mesh.o: Data_Mesh.f90
	$(f90comp) -c Data_Mesh.f90

Data_Prob.mod: Data_Prob.o Data_Prob.f90
	$(f90comp) -c Data_Prob.f90
Data_Prob.o: Data_Prob.f90
	$(f90comp) -c Data_Prob.f90

Para.mod: Para.o Para.f90
	$(f90comp) -c Para.f90
Para.o: Para.f90
	$(f90comp) -c Para.f90

List.mod: List.o List.f90
	$(f90comp) -c List.f90
List.o: List.f90
	$(f90comp) -c List.f90

Math.mod: Math.o Math.f90
	$(f90comp) -c Math.f90
Math.o: Math.f90
	$(f90comp) -c Math.f90

Mani.mod: Mani.o Mani.f90
	$(f90comp) -c Mani.f90
Mani.o: Mani.f90
	$(f90comp) -c Mani.f90

SpanTree.mod: SpanTree.o SpanTree.f90
	$(f90comp) -c SpanTree.f90
SpanTree.o: SpanTree.f90
	$(f90comp) -c SpanTree.f90

Importer.mod: Importer.o Importer.f90
	$(f90comp) -c Importer.f90
Importer.o: Importer.f90
	$(f90comp) -c Importer.f90

Exam_PERDIX.mod: Exam_PERDIX.o Exam_PERDIX.f90
	$(f90comp) -c Exam_PERDIX.f90
Exam_PERDIX.o: Exam_PERDIX.f90
	$(f90comp) -c Exam_PERDIX.f90

Exam_PERDIX.mod: Exam_PERDIX.o Exam_PERDIX.f90
	$(f90comp) -c Exam_PERDIX.f90
Exam_PERDIX.o: Exam_PERDIX.f90
	$(f90comp) -c Exam_PERDIX.f90

Section.mod : Section.o Section.f90
	$(f90comp) -c Section.f90
Section.o : Section.f90
	$(f90comp) -c Section.f90

Input.mod : Input.o Input.f90
	$(f90comp) -c Input.f90
Input.o : Input.f90
	$(f90comp) -c Input.f90

ModGeo.mod : ModGeo.o ModGeo.f90
	$(f90comp) -c ModGeo.f90
ModGeo.o : ModGeo.f90
	$(f90comp) -c ModGeo.f90

Basepair.mod : Basepair.o Basepair.f90
	$(f90comp) -c Basepair.f90
Basepair.o : Basepair.f90
	$(f90comp) -c Basepair.f90

Route.mod : Route.o Route.f90
	$(f90comp) -c Route.f90
Route.o : Route.f90
	$(f90comp) -c Route.f90

SeqDesign.mod : SeqDesign.o SeqDesign.f90
	$(f90comp) -c SeqDesign.f90
SeqDesign.o : SeqDesign.f90
	$(f90comp) -c SeqDesign.f90

Output.mod : Output.o Output.f90
	$(f90comp) -c Output.f90
Output.o : Output.f90
	$(f90comp) -c Output.f90

PERDIX.mod : PERDIX.o PERDIX.f90
	$(f90comp) -c PERDIX.f90
PERDIX.o : PERDIX.f90
	$(f90comp) -c PERDIX.f90

# Cleaning everything
clean :
	rm *.o
	rm *.mod
# End of the makefile
