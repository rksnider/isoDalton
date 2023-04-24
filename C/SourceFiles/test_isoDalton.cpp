/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  test_isoDalton.cpp                                      */
/*               Source code for testing the isoDalton.cpp code          */
/*               The bovine insulin example is run with 1000 states.     */
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

#include "isoDalton.h"

int main(int argc, char **argv)
{

	FILE * pFile;
	char *pathfilename;
	char *SpectraSimPath;
	char *DataPath;
	char *DataPathUser;
	char *UserCompFilename;
	size_t Nchar;
	struct element_list   Elements;
	struct element_list *pElements;
	struct molecule_info Molecule;
	struct molecule_info *pMolecule;
	struct istates_info isostates;
	struct istates_info *pisostates;
	int state_index;
	char  formula[50];
	char ch;
	int print_this;
	int element_index;
	int Nstates;
	int log10flag;


	//--------------------------------------------------------------------------
	// Create the paths to the data directories
    // Note: Paths are hard coded and need to be changed.
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
    UserCompFilename = (char *)malloc(50*sizeof(char));
	strcpy(UserCompFilename,"UserIsotopesNIST_HCNOS.xml");


	pElements = &Elements;
	isoDalton_get_isotopes(DataPath, DataPathUser, UserCompFilename, pElements);


	pMolecule = &Molecule;
	strcpy(formula,"C 254 H 378 N 65 O 75 S 6");
    isoDalton_parse_molecular_formula(formula, pMolecule, pElements);

	print_this=1;
	if(print_this==1){
		printf("-----------------------------------------------------------\n");
		printf("Molecular Formula = %s\n",pMolecule->Formula);
		printf("    Element Total = %d\n",pMolecule->ElementTotal);
		for(element_index=0; element_index<pMolecule->ElementTotal; element_index++){
			printf("    Element %s has %3d atoms and %2d nonzero isotopes\n",pElements->Element[pMolecule->AtomicNumber[element_index]].Symbol,pMolecule->AtomCount[element_index],pElements->Element[pMolecule->AtomicNumber[element_index]].NonzeroIsotopeTotal);
		}
		printf("-----------------------------------------------------------\n");
	}

 	Nstates    = 10000;
	log10flag  = 0;  // set to one to carry terms in the log10 domain (terms added), set to zero if probabilities are to be multiplied
    pisostates = &isostates;
	pisostates->StateTotal = Nstates;
    pisostates->mass       = (double *)malloc(Nstates*sizeof(double));
	pisostates->prob       = (double *)malloc(Nstates*sizeof(double));
    isoDalton_exact_mass(pMolecule, pElements, Nstates, pisostates, log10flag);

	//printf("-----------------------------------------------------------\n");
	//printf("Most probable masses:\n");
	//for(state_index=Nstates-1; state_index>=0; state_index--){
	//	printf("%4d mass=%17.15f prob=%17.15f\n",state_index,pisostates->mass[state_index],pisostates->prob[state_index]);
	//}
	//printf("-----------------------------------------------------------\n");

	pathfilename = (char *)malloc((strlen(DataPathUser)+15)*sizeof(char));
	strcpy(pathfilename,DataPathUser);
	strcpy(pathfilename,"\\");
	strcat(pathfilename,"test.txt");
	pFile = fopen (pathfilename,"w");
	if (pFile!=NULL){
		for(state_index=Nstates-1; state_index>=0; state_index--){
			fprintf(pFile,"%20.15f %20.15f\n",pisostates->mass[state_index],pisostates->prob[state_index]);
		}
	}
	

	printf("\n press the key e to exit \n");
	while(1){
		ch = _getch();
		if( 'e' == ch ){
			break;
		}
	}

    return 0;
}
