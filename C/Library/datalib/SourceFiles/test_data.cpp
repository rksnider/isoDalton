/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  test_data.cpp                                           */
/*               Source code for testing the functions that read         */
/*               the data files                                          */
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

#include "data.h"

int main(int argc, char **argv)
{
	char ch;
	char *SpectraSimPath;
	char *DataPath;
	char *DataPathUser;
	char *FractFilename;
	size_t Nchar;
	struct element_list   Elements;
	struct element_list *pElements;

    FractFilename = (char *)malloc(50*sizeof(char));
	//--------------------------------------------------------------------------
	// Create the paths to the data directories
	// Warning : Hard coded paths.
	//--------------------------------------------------------------------------
	Nchar = strlen("C:\\users\\ross\\MSU\\isoDalton\\C")+1;
    SpectraSimPath = (char *)malloc(Nchar*sizeof(char));
	strcpy(SpectraSimPath,"C:\\users\\ross\\MSU\\isoDalton\\C");
    DataPath = (char *)malloc((strlen(SpectraSimPath)+15)*sizeof(char));
	strcpy(DataPath,SpectraSimPath);
	strcat(DataPath,"\\DataFiles");
    DataPathUser = (char *)malloc((strlen(SpectraSimPath)+15)*sizeof(char));
	strcpy(DataPathUser,SpectraSimPath);
	strcat(DataPathUser,"\\DataFilesUser");

	//-----------------------------------------------------------------------
	// Read in element data from NIST
	// Note 1: This function should be called before data_read_AtomTabl as
	//         it allocates memory for all possible isotopes whereas
	//         the file AtomTabl only lists the isotopes with nonzero
	//         fractions and we want to list all possible isotopes since
	//         the user can change the composition fractions and memory
	//         needs to be allocated for this possibility.
	// Note 2: Any discrepancies between NIST and AtomTabl.xml will be
	//         noted and will default to AtomTabl.xml
	//-----------------------------------------------------------------------
	pElements = &Elements;
	data_read_NIST(DataPath, pElements);
	data_read_AtomTabl(DataPath, pElements);

	//-----------------------------------------------------------------------
	// Create a default User Isotopic file
	//-----------------------------------------------------------------------
	strcpy(FractFilename, "UserIsotopesALL.xml");
    data_write_UserIsotopes(DataPathUser, FractFilename, pElements);

	//-----------------------------------------------------------------------
	// Read a User Isotopic Composition Fraction file
	//-----------------------------------------------------------------------
	strcpy(FractFilename, "UserIsotopesNIST_HCNOS.xml");
    data_read_UserIsotopes(DataPathUser, FractFilename, pElements);
	
	//-----------------------------------------------------------------------
	//Make sure the Isotopic Composition Fractions sum to 1.0
	//-----------------------------------------------------------------------
	data_normalize_fractions(pElements);



	printf("\n press the key e to exit \n");
	while(1){
		ch = _getch();
		if( 'e' == ch ){
			break;
		}
	}

    return 0;
}
