# Makefiles for PERDIX-2L
objects = Data_Bound.o Data_DNA.o Data_Geom.o Data_Mesh.o Data_Prob.o \
	Para.o List.o Mani.o Math.o SpanTree.o Importer.o \
	Exam_2D_Open.o \
	3_Section.o 1_Input.o 2_ModGeo.o 4_Basepair.o 5_Route.o \
	6_SeqDesign.o 7_Output.o PERDIX-2L.o

f90comp = ifort

# Makefile
PERDIX-2L: $(objects)
	$(f90comp) -o PERDIX-2L $(objects)

Data_Bound.mod: Data_Bound.o Data_Bound.f90
	$(f90comp) -c Data_Bound.f90
Data_Bound.o: Data_Bound.f90
	$(f90comp) -c Data_Bound.f90

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

Exam_2D_Open.mod: Exam_2D_Open.o Exam_2D_Open.f90
	$(f90comp) -c Exam_2D_Open.f90
Exam_2D_Open.o: Exam_2D_Open.f90
	$(f90comp) -c Exam_2D_Open.f90

3_Section.mod : 3_Section.o 3_Section.f90
	$(f90comp) -c 3_Section.f90
3_Section.o : 3_Section.f90
	$(f90comp) -c 3_Section.f90

1_Input.mod : 1_Input.o 1_Input.f90
	$(f90comp) -c 1_Input.f90
1_Input.o : 1_Input.f90
	$(f90comp) -c 1_Input.f90

2_ModGeo.mod : 2_ModGeo.o 2_ModGeo.f90
	$(f90comp) -c 2_ModGeo.f90
2_ModGeo.o : 2_ModGeo.f90
	$(f90comp) -c 2_ModGeo.f90

4_Basepair.mod : 4_Basepair.o 4_Basepair.f90
	$(f90comp) -c 4_Basepair.f90
4_Basepair.o : 4_Basepair.f90
	$(f90comp) -c 4_Basepair.f90

5_Route.mod : 5_Route.o 5_Route.f90
	$(f90comp) -c 5_Route.f90
5_Route.o : 5_Route.f90
	$(f90comp) -c 5_Route.f90

6_SeqDesign.mod : 6_SeqDesign.o 6_SeqDesign.f90
	$(f90comp) -c 6_SeqDesign.f90
6_SeqDesign.o : 6_SeqDesign.f90
	$(f90comp) -c 6_SeqDesign.f90

7_Output.mod : 7_Output.o 7_Output.f90
	$(f90comp) -c 7_Output.f90
7_Output.o : 7_Output.f90
	$(f90comp) -c 7_Output.f90

PERDIX-2L.mod : PERDIX-2L.o PERDIX-2L.f90
	$(f90comp) -c PERDIX-2L.f90
PERDIX-2L.o : PERDIX-2L.f90
	$(f90comp) -c PERDIX-2L.f90

# Cleaning everything
clean :
	-rm *.o $(objects)
	-rm *.mod $(objects)
# End of the makefile