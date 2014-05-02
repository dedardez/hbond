Halogen Interaction Angle Calculator
------------------------------------
*How to analyze PDB files*
1. Create a folder on your computer that contains only the pdb files that you wish to analyze
2. Open hbond2.py
3. Select the folder where you placed all of your pdb files
4. Select the folder in which you wish the output text file to be saved in
5. Click on the black console window and type in your preferred name for the output file
<Note: Do not follow your name with .txt,.py,etc. Simply type a name consisting of letter, spaces, and/or numbers>
6.Press Enter
7.The program will begin to analyze the given files(This may take some times...it depends on how many files you are analyzing)
8.A prompt will appear that confirms analysis completion
9. You may now close the black console window
10. You should see a .txt file that has been created within the folder, and with the file name that you have specified

*How to read the Output*

Line            Meaning
****            ************************************************
1-4             ID associated with a given PDB file
6-11            atomic serial number for the ligand carbon
14-16           full name of the ligand carbon
18-20           name of the ligand residue
22-23           atom type(In this column, it will always be C)
25-30           atomic serial number of the ligand halogen
32-35           full name of ligand halogen
37-39           name of ligand residue
41-42           atom type(In this column, this wil be a halogen)
47-52           atomic serial number of proteinacious element
55-57           full name of protein atom
59-61           amino acid residue
63-66           amino acid sequence number
68-69           protein atom type
71-76           bond angle of a C-X--D system(in degrees)
80-83           distance between X and D(in angstroms)
85-90           bond angle of a C-X--pi system(where applicable)

*Special Character*

<Note: If the symbol 'R' appears in line 55-57, it means that the halogen interacting with a proteinacious
       benzene ring.>


Example of C-X--D output:

         |10       |20       |30       |40       |50       |60       |70       |80       |90              
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
PDB  Het#    Het Res Ty Het#   Het  Res Ty  | Pro#    Atm Res R#   Ty C-X--D   Dis  C-X--ã
115D     46  C5  BRU  C     50 BR   BRU BR  |     31  N9   DG    2  N  81.87   3.92
115D     46  C5  BRU  C     50 BR   BRU BR  |     32  C8   DG    2  C  98.04   3.81
115D     46  C5  BRU  C     50 BR   BRU BR  |     33  N7   DG    2  N  93.06   3.77

Example of C-X--pi output:
         |10       |20       |30       |40       |50       |60       |70       |80       |90              
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
PDB  Het#    Het Res Ty Het#   Het  Res Ty  | Pro#    Atm Res R#   Ty C-X--D   Dis  C-X--ã
1WBG   2383  C11 L03  C   2384 CL12 L03 CL  |          R  TYR  228    169.47   4.09  65.54
1YL1   1010  C1  ETF  C   1008  F2  ETF  F  |          R  TRP  108     86.98   5.69 -47.89
2PIT   2079  C9  4HY  C   2087  I3  4HY  I  |          R  TYR  834    164.21   3.93  55.39







