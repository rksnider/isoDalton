/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  isoDalton.cpp                                           */
/*               Source code for isoDalton, which has been translated    */
/*               from the original Matlab version                        */
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
#include "sort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <conio.h>  // for testing
#include <math.h>
#include <float.h>
#include <time.h>

//--------------------------------------------------------
// Create the isotope information
//--------------------------------------------------------
void isoDalton_get_isotopes(char *datapath, char *userdatapath, char *usercompfilename, struct element_list *pElements){

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
	data_read_NIST(datapath, pElements);
	data_read_AtomTabl(datapath, pElements);

	//-----------------------------------------------------------------------
	// Create a default User Isotopic Compostion Fraction file
	//-----------------------------------------------------------------------
	//strcpy(FractFilename, "UserIsoFraction1.xml");
    //data_write_UserCompFraction(DataPathUser, FractFilename, pElements);

	//-----------------------------------------------------------------------
	// Read a User Isotopic Composition Fraction file
	//-----------------------------------------------------------------------
    data_read_UserIsotopes(userdatapath, usercompfilename, pElements);
	
	//-----------------------------------------------------------------------
	//Make sure the Isotopic Composition Fractions sum to 1.0
	//-----------------------------------------------------------------------
	data_normalize_fractions(pElements);





}

void isoDalton_parse_molecular_formula(char *molecular_formula, struct molecule_info *pMolecule, struct element_list *pElements){
	int i;
	char *pch,*pchstart,*pchend;
	char symbol[10];
	int char_count;
	int element_index;
	int element_index_saved;
	int found_flag;
	int char_index;
	int Nch;

	pMolecule->Formula = (char *)malloc((strlen(molecular_formula)+1)*sizeof(char));
	strcpy(pMolecule->Formula,molecular_formula);
	pMolecule->ElementTotal = 0;
	pMolecule->AtomicNumber = NULL;
	pMolecule->AtomCount    = NULL;
	//printf("%s\n",pMolecule->Formula);

	pch = molecular_formula;
	if(isspace(pch[0])){
		while(isspace(*pch)){
			pch++;// advance to first letter
		}
	}
	Nch = (int)strlen(pch);
	for(char_index=0; char_index<Nch; ){ // note char_index is incremented below
		//---------------------------------------
		// Find element symbol
		//---------------------------------------
		pchstart = pchend = pch;
		char_count=0;
		while(isalpha(*pchend)){
			pchend++;
			char_index++;
			char_count += 1;
		}
		for(i=0; i<char_count; i++){
			symbol[i]=*(pchstart+i);
		}
		symbol[char_count]='\0';
		//printf("symbol=%s\n",symbol);
		pMolecule->ElementTotal += 1;
		if( 1 == pMolecule->ElementTotal ){
			pMolecule->AtomicNumber = (int *)malloc(pMolecule->ElementTotal*sizeof(int));
			pMolecule->AtomCount    = (int *)malloc(pMolecule->ElementTotal*sizeof(int));
		}else{
			pMolecule->AtomicNumber = (int *)realloc(pMolecule->AtomicNumber,pMolecule->ElementTotal*sizeof(int));
			pMolecule->AtomCount    = (int *)realloc(pMolecule->AtomCount,pMolecule->ElementTotal*sizeof(int));
		}
		//---------------------------------------
		// Find element symbol in element list
		//---------------------------------------
		found_flag = 0;
		for(element_index=1; element_index<ELEMENT_TOTAL; element_index++){  // start at 1 since 0=electron
			if( (element_index != 117) && (0 == strcmp(symbol,pElements->Element[element_index].Symbol)) ){
				found_flag          = 1;
				element_index_saved = element_index;
				break;
			}
		}
		//printf("element_index_saved = %d\n",element_index_saved);
		//printf("index=%d  symbol=%s\n",element_index_saved,pElements->Element[element_index_saved].Symbol);
		pMolecule->AtomicNumber[pMolecule->ElementTotal-1]=element_index_saved;
		//---------------------------------------
		// Find element count
		//---------------------------------------
		pchstart = pchend;
		while(isspace(*pchstart)){
			pchstart++;// advance to next number
			char_index++;
		}
		pchend = pchstart;
		char_count=0;
		while(isdigit(*pchend)){
			pchend++;
			char_index++;
			char_count += 1;
		}
		for(i=0; i<char_count; i++){
			symbol[i]=*(pchstart+i);
		}
		symbol[char_count]='\0';
		//printf("symbol=%s\n",symbol);
		//printf("number=%d\n",atoi(symbol));
		pMolecule->AtomCount[pMolecule->ElementTotal-1]=atoi(symbol);

		pch = pchend;
		while(isspace(*pch)){
			pch++;// advance to first letter
			char_index++;
		}
		//printf("char_index=%d\n",char_index);
	}
}

void isoDalton_combine_masses(int *Nelements, double *mass, double *prob, int log10flag) {
	// Note: this function assumes the mass values have been sorted in ascending order
	int mindex;
	int start_index;
	int insert_index;
	int stop_index;
	double start_mass;
	double msum,psum;
	int Nelements2;
	int Celements;
	double mass_threshold;

	start_index  = 0;
	stop_index   = 1;
	insert_index = 0;
	Nelements2   = *Nelements;

	while( stop_index <= *Nelements ){
		start_mass = mass[start_index];
		mass_threshold = start_mass/pow(10.0,15.0);  // double precision eps
		if( stop_index < *Nelements ){
			while( fabs(mass[stop_index]-start_mass) <= mass_threshold ){
				//printf("     start_index = %d stop_index = %d\n",start_index,stop_index);
				//printf("     mass[stop_index] = %f start_mass = %f\n",mass[stop_index], start_mass);
				//printf("     difference = %f\n",(mass[stop_index]-start_mass));
				stop_index++;
				if(stop_index >= *Nelements){
					break;
				}
			}
		}
		Celements = stop_index - start_index; // the number of elements being combined
		//printf("Celements=%d insert_index=%d\n",Celements,insert_index);
		if( Celements > 1 ){ // condense mass block
			// average masses
			msum = 0;
			for(mindex=start_index; mindex<stop_index; mindex++){
				msum += mass[mindex];
			}
			mass[insert_index] = msum/(double)Celements;
			// add probs that are in log10 domain
			psum = 0;
			if(1 == log10flag ){
				for(mindex=start_index; mindex<stop_index; mindex++){
					psum += pow(10,prob[mindex]);
				}
				prob[insert_index] = log10(psum);
			}else{
				for(mindex=start_index; mindex<stop_index; mindex++){
					psum += prob[mindex];
				}
				prob[insert_index] = psum;
			}
			// update element count
			Nelements2 = Nelements2 - Celements + 1;
			//printf("       Inserting mass=%f prob=%f into location %d\n",mass[insert_index],prob[insert_index],insert_index);
			//printf("       Nelements2=%d\n",Nelements2);
		}else{
			mass[insert_index] = mass[start_index];
			prob[insert_index] = prob[start_index];
		}
		insert_index++;
		start_index = stop_index;
		stop_index  = start_index+1;
		if( stop_index > *Nelements ){
			mass[insert_index] = mass[start_index];
			prob[insert_index] = prob[start_index];
			break;
		}
	}

    *Nelements = Nelements2;
}


void isoDalton_exact_mass(struct molecule_info *pMolecule, struct element_list *pElements, int Mstates, struct istates_info *pisostates, int log10flag){

	int Nelements;
	int Natoms;
	int Nisotopes,maxNisotopes;
	int isotope_index;
	int index,index1,index2,index3;
	int *NonzeroIsotopeTotal;
	int *Eindex,*Mindex;
	double mass_min;
	double mass_max;
	double prob_min;
	double prob_max;
	double prob_isotope;
	double mass_isotope;
	double term_lightest;
	double term_heaviest;
	double term_most_probable;
	double term_least_probable;
	double distribution_span;
	double *state1_mass,*state2_mass,*state3_mass;
	double *state1_prob,*state2_prob,*state3_prob;
	double *average_mass1,*average_mass2;
	int state1_index,state2_index;
	int Nstate1,Nstate2;
	int state_index;
	clock_t time0,time1;

	Nelements=pMolecule->ElementTotal;
	NonzeroIsotopeTotal = (int *)malloc(Nelements*sizeof(int));
	Eindex              = (int *)malloc(Nelements*sizeof(int));
	Mindex              = (int *)malloc(Nelements*sizeof(int));
	average_mass1       = (double *)malloc(Nelements*sizeof(double));
	average_mass2       = (double *)malloc(Nelements*sizeof(double));
	for(index1=0; index1<Nelements; index1++){
		NonzeroIsotopeTotal[index1] = pElements->Element[pMolecule->AtomicNumber[index1]].NonzeroIsotopeTotal;
		Eindex[index1]              = pMolecule->AtomicNumber[index1];
		Mindex[index1]              = index1;
		average_mass1[index1]       = pElements->Element[pMolecule->AtomicNumber[index1]].AverageMass;
	}

	//--------------------------------------------------------------
	// Get index that can sort based in number of nonzero isotopes
	//--------------------------------------------------------------
	heapsort_3int_up(Nelements, NonzeroIsotopeTotal, Eindex, Mindex);
	for(index1=0; index1<Nelements; index1++){
		average_mass2[index1] = average_mass1[Mindex[index1]];
		//printf("  %f\n",average_mass2[index1]);
	}
	heapsort_1dbl_2int_up(Nelements, average_mass2, Eindex, Mindex);  // do a secondary sort on increasing mass

	printf("The molecular elements sorted by increasing isotope numbers (and then by mass)\n");
	for(index1=0; index1<Nelements; index1++){
		printf("Element %10s has %2d nonzero isotopes.\n",pElements->Element[Eindex[index1]].Name, pElements->Element[Eindex[index1]].NonzeroIsotopeTotal);
	}
	printf("-----------------------------------------------------------\n");

	//-----------------------------------------
	// Get mass and probability spanning info
	//-----------------------------------------
	term_lightest       = 0;
	term_heaviest       = 0;
	term_most_probable  = 0;
	term_least_probable = 0;
	for(index1=0; index1<Nelements; index1++){
		Natoms      = pMolecule->AtomCount[Mindex[index1]];
		Nisotopes   = pElements->Element[Eindex[index1]].NonzeroIsotopeTotal;
		mass_min =  DBL_MAX;
		mass_max = -DBL_MAX;
		prob_min =  DBL_MAX;
		prob_max = -DBL_MAX;
		for(index2=0; index2<Nisotopes; index2++){
			index3       = pElements->Element[Eindex[index1]].NonzeroIsotopeIndex[index2];
			mass_isotope = pElements->Element[Eindex[index1]].Isotope[index3]->AtomicMass;
			prob_isotope = pElements->Element[Eindex[index1]].Isotope[index3]->CompositionFraction;
			printf("%10s(%2d) isotope(%2d) mass %17.15f fraction %17.10f \n",pElements->Element[Eindex[index1]].Name,pElements->Element[Eindex[index1]].AtomicNumber,pElements->Element[Eindex[index1]].Isotope[index3]->MassNumber,pElements->Element[Eindex[index1]].Isotope[index3]->AtomicMass,pElements->Element[Eindex[index1]].Isotope[index3]->CompositionFraction);
			if( mass_min > mass_isotope ){
				mass_min = mass_isotope;
			}
			if( mass_max < mass_isotope ){
				mass_max = mass_isotope;
			}
			if( prob_min > prob_isotope ){
				prob_min = prob_isotope;
			}
			if( prob_max < prob_isotope ){
				prob_max = prob_isotope;
			}
		}
		//printf("      mass_min = %f\n",mass_min);
		//printf("      mass_max = %f\n",mass_max);
		//printf("      prob_min = %f\n",prob_min);
		//printf("      prob_max = %f\n",prob_max);
		//printf("      Natoms = %d  (%f)\n",Natoms,(double)Natoms);
		term_lightest       = term_lightest +(double)Natoms*mass_min;
		term_heaviest       = term_heaviest + (double)Natoms*mass_max;
		term_most_probable  = term_most_probable  + (double)Natoms*log10(prob_max);
		term_least_probable = term_least_probable + (double)Natoms*log10(prob_min);
		//printf("      term_lightest       = %f\n",term_lightest);
		//printf("      term_lightest       = %f\n",term_heaviest);
		//printf("      term_most_probable  = %f\n",term_most_probable);
		//printf("      term_least_probable= %f\n",term_least_probable);
	}
	printf("-----------------------------------------------------------\n");
	printf("Information regarding molecule [%s]\n",pMolecule->Formula);
	printf("The lightest mass term = %f daltons\n",term_lightest);
	printf("The heaviest mass term = %f daltons\n",term_heaviest);
	distribution_span = term_heaviest - term_lightest;
	printf("The isotopic distribution spans = %f daltons\n",distribution_span);
	printf("The most  probable term (log10) = %f\n",term_most_probable);
	printf("The least probable term (log10) = %f\n",term_least_probable);
	printf("-----------------------------------------------------------\n");


	//---------------------------------------------------------
	// Setup state vectors
	// The max number of states to preallocate will be 
	// Mstates*maxNisotopes where maxNisotopes is maximum 
	// number of isotopes over all elements in the molecule.
	//---------------------------------------------------------
	maxNisotopes = 0;
	for(index1=0; index1<Nelements; index1++){
		Nisotopes = pElements->Element[Eindex[index1]].NonzeroIsotopeTotal;
		if( Nisotopes > maxNisotopes){
			 maxNisotopes = Nisotopes;
		}
	}
	printf("maxNisotope = %d\n",maxNisotopes);
	state1_mass = (double *)malloc(Mstates*maxNisotopes*sizeof(double));
	state1_prob = (double *)malloc(Mstates*maxNisotopes*sizeof(double));
	state2_mass = (double *)malloc(Mstates*maxNisotopes*sizeof(double));  //****findme**** free memory
	state2_prob = (double *)malloc(Mstates*maxNisotopes*sizeof(double));

	time0 = clock();
	//---------------------------------------------------------
	// Load initial states
	//---------------------------------------------------------
	Nisotopes   = pElements->Element[Eindex[0]].NonzeroIsotopeTotal;
	for(index=0; index<Nisotopes; index++){
		index2            = pElements->Element[Eindex[0]].NonzeroIsotopeIndex[index];
		state1_mass[index] = pElements->Element[Eindex[0]].Isotope[index2]->AtomicMass;
		if(1 == log10flag ){
		    state1_prob[index] = log10(pElements->Element[Eindex[0]].Isotope[index2]->CompositionFraction);
		}else{
			state1_prob[index] = pElements->Element[Eindex[0]].Isotope[index2]->CompositionFraction;
		}
	}
	printf("Initialized lattice with Element %s\n",pElements->Element[Eindex[0]].Name);
	Nstate1=Nisotopes;
	Nstate2=0;
	for(index=0; index<Nisotopes; index++){
		if(1 == log10flag ){
		    printf("state %d mass=%f log10(prob)=%f\n",index,state1_mass[index],state1_prob[index]);
		}else{
		    printf("state %d mass=%f prob=%f\n",index,state1_mass[index],state1_prob[index]);
		}
	}

	//---------------------------------------------------------
	// Trellis
	//---------------------------------------------------------
	for(index1=0; index1<Nelements; index1++){
		Natoms      = pMolecule->AtomCount[Mindex[index1]];
		Nisotopes   = pElements->Element[Eindex[index1]].NonzeroIsotopeTotal;
		for(index2=0; index2<Natoms; index2++){
			if( !(index1==0 && index2==0) ){ // skip the first entry since this was preloaded in the trellis
				//printf("Nstate1=%d  Nisotopes=%d\n",Nstate1,Nisotopes);
				//---------------------------------------------------------
				// Expand state1 by Nisotopes
				//---------------------------------------------------------
				state2_index=0;
				for(state1_index=0; state1_index<Nstate1; state1_index++){
					for(isotope_index=0; isotope_index<Nisotopes; isotope_index++){
						index3                    = pElements->Element[Eindex[index1]].NonzeroIsotopeIndex[isotope_index];
						state2_mass[state2_index] = state1_mass[state1_index] + pElements->Element[Eindex[index1]].Isotope[index3]->AtomicMass;
						if(1 == log10flag ){
							state2_prob[state2_index] = log10(  pow(10,state1_prob[state1_index]) * pElements->Element[Eindex[index1]].Isotope[index3]->CompositionFraction  );
						}else{
							state2_prob[state2_index] = state1_prob[state1_index] * pElements->Element[Eindex[index1]].Isotope[index3]->CompositionFraction;
						}
						state2_index++;
					}
				}
				Nstate2 = state2_index;

				//printf("Nstates2 = %d\n",Nstate2);
				//for(state2_index=0; state2_index<Nstate2; state2_index++){
				//	printf("%3d %17.15f %17.15f\n",state2_index,state2_mass[state2_index],state2_prob[state2_index]);
				//}

				//---------------------------------------------------------
				// sort state2 by ascending masses
				//---------------------------------------------------------
				heapsort_2dbl_up(Nstate2, state2_mass, state2_prob);

				//printf("Nstates2 = %d   sorted by mass\n",Nstate2);
				//for(state2_index=0; state2_index<Nstate2; state2_index++){
				//	printf("%3d %f %f\n",state2_index,state2_mass[state2_index],state2_prob[state2_index]);
				//}

				//------------------------------------------------------------
				// combine mass states that are closer than a mass threshold
				//------------------------------------------------------------
				isoDalton_combine_masses(&Nstate2, state2_mass, state2_prob, log10flag);

				//printf("Nstates2 = %d   mass combined\n",Nstate2);
				//for(state2_index=0; state2_index<Nstate2; state2_index++){
				//	printf("%3d %f %f\n",state2_index,state2_mass[state2_index],state2_prob[state2_index]);
				//}

				//---------------------------------------------------------
				// sort state2 by decending probability
				//---------------------------------------------------------
				heapsort_2dbl_down(Nstate2, state2_prob, state2_mass);

				//printf("sorted by probability\n");
				//for(state2_index=0; state2_index<Nstate2; state2_index++){
				//	printf("%3d %17.15f %17.15f\n",state2_index,state2_mass[state2_index],state2_prob[state2_index]);
				//}


				//printf("Element %10s Atom Count %d\n",pElements->Element[Eindex[index1]].Name, index2);


			//printf("\n\n\n press the key c to continue \n");
			//while(1){
			//	ch = _getch();
			//	if( 'c' == ch ){
			//		break;
			//	}
			//}

			    //----------------------------------------------
			    // Swap pointers so state2 becomes state1
				// and vice versa and cap number of states
			    //----------------------------------------------
			    state3_mass = state1_mass;
			    state3_prob = state1_prob;
			    state1_mass = state2_mass;
			    state1_prob = state2_prob;
			    state2_mass = state3_mass;
			    state2_prob = state3_prob;
				if( Nstate2 > Mstates ){
                    Nstate1 = Mstates;
				}else{
					Nstate1 = Nstate2;
				}
			}
		}
	}
	time1 = clock();
	printf("clock rate of timer = %d\n",CLOCKS_PER_SEC);
	printf("It took %8.4f seconds to compute isotope spectra\n",(double)(time1-time0)/(double)(CLOCKS_PER_SEC));
	printf("Number of States = %d\n",Mstates);


	for(state_index=0; state_index<Mstates; state_index++){
		pisostates->mass[state_index] = state1_mass[state_index];
		pisostates->prob[state_index] = state1_prob[state_index];
	}

}			   
			   



