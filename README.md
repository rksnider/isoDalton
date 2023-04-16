# isoDalton

The Matlab files implement isoDalton, which computes the isotopic mass distributions
of a given molecule.

This software is associated with the following paper:
Snider, R.K. Efficient Calculation of Exact Mass Isotopic Distributions
J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.
The digital object identifier (DOI) link to paper:  http://dx.doi.org/10.1016/j.jasms.2007.05.016

-------------------------------------------------------------------------------
Matlab Installation:
-------------------------------------------------------------------------------
Create a directory named isoDalton, copy all the files in the Matlab directory
into \isoDalton, and then set Matlab's path to this directory.



A C version of isoDalton has been developed.  isoDalton can generate the fine isotopic structure given any molecular formula.  It is much faster than the Matlab version. The table below gives the speedup for bovine insulin.

![image](https://user-images.githubusercontent.com/5913180/232347394-7e8b52d2-05a9-4d35-b98c-d4f921f36720.png)


IsoDalton run times for Bovine Insulin [C254 H378 N65 O75 S6]
States
(most probable masses)	Matlab
(time in seconds)	C
(time in seconds)	Speedup
100	4.6499	0.0780	59.60
1,000	28.9499	0.5780	50.08
10,000	239.7523	7.0470	34.02
Matlab (2009b), C program using Visual Studio 2005, run on a Windows XP computer (Dell Precision WorkStation 390, Core 2 CPU 6300 @ 1.862 GHz, 4 GB RAM)

-------------------------------------------------------------------------------
C Installation:
-------------------------------------------------------------------------------
See document C_installation in C directory


-------------------------------------------------------------------------------
There are three demos:
-------------------------------------------------------------------------------
demo_glycine1.m          
demo_glycine2.m
demo_bovine_insuline.m



-------------------------------------------------------------------------------
Usage (exact mass):
isoDalton_exact_mass(molecule_string,maxstates);
-------------------------------------------------------------------------------

 Input:  molecule_string :
         String of elements <Ex> (any of the first 112 on the periodic chart of the elements)
         paired with the number of atoms for each element <Cx>, i.e. ['E1 C1 E2 C2 E3 C3 etc.'].  
         Spaces are not required, i.e. all the following cases are allowed, although the top case is 
         the most readable.:
         molecule_string = 'C254 H377 N65 O75 S6 Fe2';
         molecule_string = 'C254H377N65O75S6Fe2';
         molecule_string = 'C 254 H 377 N 65 O 75 S 6 Fe 2';

         maxstates : maximum number of mass states
                     set maxstates = realmax; if one wishes to keep *ALL* the states.  Note that this
                     is practical only for small molecules, since for large molecules this will slow
                     to a catatonic crawl and run out of memory.

---------------------------------------------------------------------------------------------------

 Output:  	returns a two column matrix where the first column contains the exact masses
           (most probable ones) and the second column are the probabilities (pruned).
---------------------------------------------------------------------------------------------------
Example:
-------------------------------------------------------------------------------
molecule_string = 'C2 H5 N1 O2';  % glycine
maxstates = 1000;          % maximum number of states
states = isoDalton_exact_mass(molecule_string,maxstates);





-------------------------------------------------------------------------------
Isotopic Compositions
-------------------------------------------------------------------------------
To change isotopic compositions, modify file isoDalton_modify_isotope_composition.m
and uncomment the example below the line: 
-----------     User Input Below this Line -------------------
and modify the example.



-------------------------------------------------------------------------------
Files:  (not complete yet....)
-------------------------------------------------------------------------------
isoDalton_NIST_all_isotopes.txt      
            list of all the isotopes for elements 1:112  
            from   http://physics.nist.gov/PhysRefData/Compositions/index.html
read by :   isoDalton_NIST_isotopes_read()
------------
isoDalton_NIST_isotopes_read() 
            reads the file isoDalton_NIST_all_isotopes.txt to get isotope information
parent   :  isoDalton_get_isotope_info()
children :  none
------------
isoDalton_get_isotope_info()
            puts the isotope information in a more usable/structured form
parent   :  isoDalton_get_molecule_isotopes()
children :  isoDalton_NIST_isotopes_read()
            isoDalton_element_symbols_read()
            isoDalton_modify_isotope_composition()
------------
isoDalton_element_symbols_read()
            reads the file isoDalton_element_symbols.txt to get element symbols
parent   :  isoDalton_get_isotope_info()
children :  none
------------
isoDalton_element_symbols.txt
            list of symbols for all the elements
            from http://www.chem.qmul.ac.uk/iupac/AtWt/
read by  :  isoDalton_element_symbols_read()
------------
isoDalton_modify_isotope_composition()
            allows the user to change isotopic compositions for the elements
parent   :  isoDalton_get_isotope_info()
children :  none
------------
isoDalton_get_molecule_isotopes()
            gets the isotope information for elements in the molecule string
parent   :  isoDalton_exact_mass()
children :  isoDalton_get_isotope_info()
            isoDalton_element_symbols_read()
            isoDalton_element_sym2num()
------------
isoDalton_element_sym2num()
            returns the atomic number for a given element symbol
parent   :  isoDalton_get_molecule_isotopes()
children :  none
------------
isoDalton_exact_mass()





