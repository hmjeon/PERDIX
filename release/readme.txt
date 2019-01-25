Three ways to launch METIS

note:
It needs the 'seq.txt' file if you want to use user-defined scaffold sequence.
Otherwise, METIS runs based on the M13mp18 scaffold sequence.

1. Clicking the icon on Window
2. Run on command without parameters

Note that two inputs are needed to run METIS, geometry file and scaffold sequence data file.
2.1. gemetry file
	The user defiend geometry file should be located at the folder named "input".
2.2. Sequence file
	The sequence data file (seq.txt) should be located at the same folder where METIS.exe exists.
	If there is no file 'seq.txt', the m13mp18 sequence is used.
	To use user's defined sequence, write 'user' at the first row in the 'seq.txt' file and write scaffold sequences in the second row.

3. Run on command with seven parameters

3.1.Command run
	METIS [geometry file] [sequence file] [edge-section] [edge] [edge length] [mesh spacing] [run mode]
	ex) METIS ./test.ply ./seq.txt 2 1 42 0.0 s

3.1.1. Option 1 - Goemetry file
	integer  - METIS uses the pre-defined geometry
	filename - METIS reads the user's geometry file
	
3.1.2. Option 2 - Sequence file
	m13     - METIS will use m13mp18 scaffold sequence without importing the sequence data
	seq.txt - METIS will import sequence data from the txt file named seq.txt
	
3.1.3. Option 3 - Edge section, METIS should be 2 as edge section
	1 - DX tile edge
	2 - 6HB edge
	
3.1.4. Option 4 - Edge
	0       - METIS find the minimum edge length and the length option will be appiled to min. edge
	integer - The user can point edge number and the following option 5 will be appiled to this edge

3.1.5. Option 5 - Edge length
	integer - The length will be applied to the edge choosen from the option 4.
	
3.1.6. Option 6 - Mesh spacing
	0.0 - 1.0 - The mesh density will be defined with this option when using boundary only design.

3.1.7. Option 7 - Run mode
	s - Single run: Design files are generated in the folder named 'outputs'
	m - Multiple run: Design files are generated in the folder names 'outputs/design_name'

3.2. Examples of running METIS with parameters
mkdir results
..\make\x64\Debug\METIS.exe  1 m13 2 0 63  0.0 m > outputs\01_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  2 m13 2 0 42  0.0 m > outputs\02_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  3 m13 2 0 42  0.0 m > outputs\03_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  4 m13 2 0 42  0.0 m > outputs\04_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  5 m13 2 0 42  0.0 m > outputs\05_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  6 m13 2 0 42  0.0 m > outputs\06_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  7 m13 2 0 42  0.0 m > outputs\07_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  8 m13 2 0 42  0.0 m > outputs\08_METIS.log 2>&1
..\make\x64\Debug\METIS.exe  9 m13 2 0 42  0.0 m > outputs\09_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 10 m13 2 0 42  0.0 m > outputs\10_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 11 m13 2 0 42  0.0 m > outputs\11_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 12 m13 2 0 42  0.0 m > outputs\12_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 13 m13 2 0 42  0.0 m > outputs\13_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 14 m13 2 0 42  0.0 m > outputs\14_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 15 m13 2 0 42  0.0 m > outputs\15_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 16 m13 2 0 42  0.0 m > outputs\16_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 17 m13 2 0 84  0.0 m > outputs\17_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 18 m13 2 0 84  0.0 m > outputs\18_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 19 m13 2 0 84  0.0 m > outputs\19_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 20 m13 2 0 128 0.0 m > outputs\20_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 21 m13 2 0 105 0.0 m > outputs\21_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 22 m13 2 0 74  0.0 m > outputs\22_METIS.log 2>&1
..\make\x64\Debug\METIS.exe 23 m13 2 0 57  0.0 m > outputs\23_METIS.log 2>&1