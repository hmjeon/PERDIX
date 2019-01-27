Three ways to launch PERDIX

note:
It needs the 'seq.txt' file if you want to use user-defined scaffold sequence.
Otherwise, PERDIX runs based on the M13mp18 scaffold sequence.

1. Clicking the icon on Window
2. Run on command without parameters

Note that two inputs are needed to run PERDIX, geometry file and scaffold sequence data file.
2.1. gemetry file
	The user defiend geometry file should be located at the folder named "input".
2.2. Sequence file
	The sequence data file (seq.txt) should be located at the same folder where PERDIX.exe exists.
	If there is no file 'seq.txt', the m13mp18 sequence is used.
	To use user's defined sequence, write 'user' at the first row in the 'seq.txt' file and write scaffold sequences in the second row.

3. Run on command with seven parameters

3.1.Command run
	PERDIX [geometry file] [sequence file] [edge] [edge length] [mesh spacing] [run mode]
	ex) PERDIX ./test.ply ./seq.txt 1 42 0.0 s

3.1.1. Option 1 - Goemetry file
	integer  - PERDIX uses the pre-defined geometry
	filename - PERDIX reads the user's geometry file
	
3.1.2. Option 2 - Sequence file
	m13     - PERDIX will use m13mp18 scaffold sequence without importing the sequence data
	seq.txt - PERDIX will import sequence data from the txt file named seq.txt
	
3.1.3. Option 3 - Edge
	0       - PERDIX find the minimum edge length and the length option will be appiled to min. edge
	integer - The user can point edge number and the following option 5 will be appiled to this edge

3.1.4. Option 4 - Edge length
	integer - The length will be applied to the edge choosen from the option 4.
	
3.1.5. Option 5 - Mesh spacing
	0.0 - 1.0 - The mesh density will be defined with this option when using boundary only design.

3.1.6. Option 6 - Run mode
	s - Single run: Design files are generated in the folder named 'outputs'
	m - Multiple run: Design files are generated in the folder names 'outputs/design_name'

3.2. Examples of running PERDIX with parameters
mkdir results
..\make\x64\Debug\PERDIX.exe  1 m13 0 63  0.0 m > outputs\01_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  2 m13 0 42  0.0 m > outputs\02_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  3 m13 0 42  0.0 m > outputs\03_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  4 m13 0 42  0.0 m > outputs\04_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  5 m13 0 42  0.0 m > outputs\05_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  6 m13 0 42  0.0 m > outputs\06_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  7 m13 0 42  0.0 m > outputs\07_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  8 m13 0 42  0.0 m > outputs\08_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe  9 m13 0 42  0.0 m > outputs\09_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 10 m13 0 42  0.0 m > outputs\10_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 11 m13 0 42  0.0 m > outputs\11_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 12 m13 0 42  0.0 m > outputs\12_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 13 m13 0 42  0.0 m > outputs\13_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 14 m13 0 42  0.0 m > outputs\14_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 15 m13 0 42  0.0 m > outputs\15_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 16 m13 0 42  0.0 m > outputs\16_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 17 m13 0 84  0.0 m > outputs\17_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 18 m13 0 84  0.0 m > outputs\18_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 19 m13 0 84  0.0 m > outputs\19_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 20 m13 0 128 0.0 m > outputs\20_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 21 m13 0 105 0.0 m > outputs\21_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 22 m13 0 74  0.0 m > outputs\22_PERDIX.log 2>&1
..\make\x64\Debug\PERDIX.exe 23 m13 0 57  0.0 m > outputs\23_PERDIX.log 2>&1