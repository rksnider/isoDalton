/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  isotopes.h                                              */
/*               Defines data structures to contain isotope information  */
/*               of the chemical elements.                               */
/*-----------------------------------------------------------------------*/
/* This software is associated with the following paper:                 */
/* Snider,R.K. Efficient Calculation of Exact Mass Isotopic Distributions*/
/* J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.                  */
/* The digital object identifier (DOI) link to the paper is:             */
/* http://dx.doi.org/10.1016/j.jasms.2007.05.016                         */
/*-----------------------------------------------------------------------*/
/* Author:       Ross Snider                                             */
/* Company:      Montana State University                                */
/* Create Date:  June 2010                                               */
/* Revision:     1.0                                                     */
/* License:      GPL-2.0 or MIT  (opensource.org/licenses/MIT)           */
/*-----------------------------------------------------------------------*/ 

#ifndef ISOTOPE_STRUCTS
     
#define ISOTOPE_STRUCTS  
// There are 118 elements + the electron at index 0
//---------------------------------------------------------------------------------------------
#define ELEMENT_TOTAL 119
//---------------------------------------------------------------------------------------------
// Structure to contain Isotope Information
//---------------------------------------------------------------------------------------------
struct isotope_info {
	//---------------------------------------------------------------------------------------------
	// Name of Isotope if there is a name. Set to None if there is no name.
	//---------------------------------------------------------------------------------------------
	char *Name;
	char *Symbol;
	//---------------------------------------------------------------------------------------------
	// http://physics.nist.gov/PhysRefData/Compositions/notes.html
	// The mass number, also called atomic mass number or nucleon number, 
	// is the total number of protons and neutrons in an atomic nucleus.
	//---------------------------------------------------------------------------------------------
	int MassNumber;   
	//---------------------------------------------------------------------------------------------
	// Relative Atomic Mass (of the isotope):  Ar(X), where X is an isotope
	// formerly called atomic weight (see Standard Atomic Weight)
        // These values are scaled to Ar(12C) = 12, where 12C is a neutral atom 
	// in its nuclear and electronic ground state. 
	// Thus, the relative atomic mass of entity X is given by: 
	// Ar(X) = m(X) / [m(12C) / 12] . 
	// The relative atomic mass of an element is derived by averaging the 
	// relative atomic masses of the isotopes of that element. 
	// Brackets [ ] indicate the mass number of the most stable isotope. 
	//---------------------------------------------------------------------------------------------
	double AtomicMass;
	//---------------------------------------------------------------------------------------------
	// Representative Isotopic Composition (%):  Mole fraction of the various isotopes
        // In the opinion of the Subcommittee for Isotopic Abundance Measurements (SIAM), 
	// these values represent the isotopic composition of the chemicals and/or materials 
	// most commonly encountered in the laboratory. They may not, therefore, 
	// correspond to the most abundant natural material. 
	//---------------------------------------------------------------------------------------------
	double CompositionFraction;
};
//---------------------------------------------------------------------------------------------
// Structure to contain information of a chemical element
//---------------------------------------------------------------------------------------------
struct element_info {
	//---------------------------------------------------------------------------------------------
	// The Atomic Number (also known as the proton number) is the number of protons  
	// found in the nucleus of an atom and therefore identical to the charge number of the nucleus
	//---------------------------------------------------------------------------------------------
	int AtomicNumber;
	char *Name;
	char *Symbol;
	//---------------------------------------------------------------------------------------------
	// AverageMass = Standard Atomic Weight (common usage):  Ar(X), where X is an element
	// derived by averaging the relative atomic masses of the isotopes of that element
	// (weighted by their composition percentages)
	//---------------------------------------------------------------------------------------------
	double AverageMass;
	//---------------------------------------------------------------------------------------------
	// Isotopes are different types of atoms (nuclides) of the same chemical element, 
	// each having a different number of neutrons.
	//---------------------------------------------------------------------------------------------
	struct isotope_info **Isotope;
	int                   IsotopeTotal;
	int MostCommonIsotopeIndex;  // index to isotope with greatest composition fraction (most common)
	int *NonzeroIsotopeIndex;    // index values of isotopes with nonzero composition fractions
	int  NonzeroIsotopeTotal;    // total number of isoteopes with nonzero composition fractions
	//---------------------------------------------------------------------------------------------
	// Element Stability 
	// Stable = True if stable, False if unstable
	// If Stable = False, then AverageMass is set to the mass of the most stable isotope
	//---------------------------------------------------------------------------------------------
	bool Stable;  // True of stable, false if unstable
};
//---------------------------------------------------------------------------------------------
// Structure to contain information of all the chemical elements
//---------------------------------------------------------------------------------------------
struct element_list {
	//------------------------------------------------------------
	// List of the Elements and their isotopes
	// Note : Element zero is an electron and elements are
	//        indexed by their atomic number (proton number)
	//------------------------------------------------------------
	struct element_info Element[ELEMENT_TOTAL];
	int                 Element_Total;
};
//---------------------------------------------------------------------------------------------
// Structure to contain information of a molecule
//---------------------------------------------------------------------------------------------
struct molecule_info {
	char *Formula;
	int   ElementTotal;
	int  *AtomCount;     // Number of atoms of each element in the molecule
	int  *AtomicNumber;  // index values of the elements in the element list Elements
};

#endif
