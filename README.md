# isoDalton
isoDalton is a program to compute the exact masses found in an isotopic distribution.  The software in this repository is associated with the following paper:
Snider, R.K. Efficient Calculation of Exact Mass Isotopic Distributions
J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.
The digital object identifier (DOI) link to paper:  http://dx.doi.org/10.1016/j.jasms.2007.05.016

Paper Abstract:
This paper presents a new method for efficiently calculating the exact masses in an isotopic distribution using a dynamic programming approach. The resulting program, isoDalton, can generate extremely high isotopic resolutions as demonstrated by a FWHM resolution of 2 × 1011. This resolution allows very fine mass structures in isotopic distributions to be seen, even for large molecules. Since the number of exact masses grows exponentially with molecular size, only the most probable exact masses are kept, the number of which is user specified.

------------------------------------------------------------------------------

A C version of isoDalton has been developed.  isoDalton can generate the fine isotopic structure given any molecular formula.  It is much faster than the Matlab version. The table below gives the speedup for bovine insulin.

![image](https://user-images.githubusercontent.com/5913180/232347394-7e8b52d2-05a9-4d35-b98c-d4f921f36720.png)

-------------------------------------------------------------------------------
C Installation:
-------------------------------------------------------------------------------
See document c_installation.pdf in the \C directory

-------------------------------------------------------------------------------
Matlab Installation:
-------------------------------------------------------------------------------
Create a directory named isoDalton, copy all the files in the Matlab directory
into \isoDalton, and then set Matlab's path to this directory.

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

 Output:  	returns a two column matrix where the first column contains the exact masses
           (most probable ones) and the second column are the probabilities (pruned).
---------------------------------------------------------------------------------------------------
Example:

molecule_string = 'C2 H5 N1 O2';  % glycine   
maxstates = 1000;          % maximum number of states   
states = isoDalton_exact_mass(molecule_string,maxstates);   





-------------------------------------------------------------------------------
Modifying Isotopic Compositions
Matlab Example
-------------------------------------------------------------------------------
To change isotopic compositions, modify file isoDalton_modify_isotope_composition.m
and uncomment the example below the line: 
-----------     User Input Below this Line -------------------
and modify the example.


