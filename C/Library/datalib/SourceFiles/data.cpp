/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  data.cpp                                                */
/*               Source code for reading all the data files              */
/*-----------------------------------------------------------------------*/
/* Author:       Ross Snider                                             */
/* Company:      Montana State University                                */
/* Create Date:  June 2010                                               */
/* Revision:     1.0                                                     */
/* License:      GPL-2.0 or MIT  (opensource.org/licenses/MIT)           */  
/*-----------------------------------------------------------------------*/

#include "data.h"

void data_normalize_fractions(struct element_list *pElements){
	int Nentries;
	int Nisotopes;
	int element_index;
	int iso_index;
	double fraction_sum;

	Nentries = pElements->Element_Total;
	for(element_index=1; element_index<Nentries; element_index++){
		if( pElements->Element[element_index].AtomicNumber > 0 ){
			Nisotopes = pElements->Element[element_index].IsotopeTotal;
			if( 0 < Nisotopes ){
				fraction_sum = 0;
				for(iso_index=0; iso_index<Nisotopes; iso_index++){
					fraction_sum += pElements->Element[element_index].Isotope[iso_index]->CompositionFraction;
				}
				if( fraction_sum != 1.0 ){
					//printf("Warning: The composition fraction sum for %s (Atomic Number %d) doesn't sum to 1.0 (%f) - It now does\n",pElements->Element[element_index].Name,pElements->Element[element_index].AtomicNumber,fraction_sum);
				}
			}
		}
	}
}


void data_read_UserIsotopes(char *path, char *filename, struct element_list *pElements){
	char *pathfilename;
	char name[30];
	XMLNode xMainNode,xNode,xNode2;
	int i;
	int Mnumber;
	int Nentries;
	int Eindex;
	int Nisotopes;
	int iso_index;
	int AtomicNumber;
	int MassNumber;
	int Nmass,Nfraction;
	double mass,imass;
	double fraction,ifraction;
	clock_t time0,time1;

	//---------------------------------------------------
	// open user isotopic composition xml file 
	//---------------------------------------------------
	pathfilename = (char *)malloc((strlen(path)+strlen(filename)+30)*sizeof(char));
	strcpy(pathfilename,path);
	strcat(pathfilename,"\\");
	strcat(pathfilename,filename);
	printf("Reading pathfile: %s\n",pathfilename);
	printf("Reading file: %s\n",filename);

	//---------------------------------------------------
    // Open and parse the XML file:
	//---------------------------------------------------
	time0 = clock();
    xMainNode=XMLNode::openFileHelper(pathfilename);
	time1 = clock();

	//---------------------------------------------------
	// Count how many Element entries there are
	//---------------------------------------------------
	Nentries = xMainNode.getChildNode("user_isotopes").nChildNode("element");
	//printf("There are %d element nodes\n", Nentries);
	//---------------------------------------------------
    // Read in each entry
	//---------------------------------------------------
	for(Eindex=0; Eindex<Nentries; Eindex++){
		//--------------------------------------------------------------------
		// Get Atomic Number
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("atomic_number");
		AtomicNumber = atoi(xNode.getText());
		//printf("Atomic Number = %d\n",AtomicNumber);
		//--------------------------------------------------------------------
		// Get Element Name
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("name");
		strcpy(name,xNode.getText());
		//printf("    Name   = %s\n",name);
		//--------------------------------------------------------------------
		// Get Element Symbol
		//--------------------------------------------------------------------
		//xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("symbol");
		//printf("    Symbol   = %s\n",xNode.getText());
		//--------------------------------------------------------------------
		// Get Number of Isotopes
		//--------------------------------------------------------------------
		 Nisotopes = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).nChildNode("isotope");
		 for(iso_index=0; iso_index<Nisotopes; iso_index++){
			//--------------------------------------------------------------------
			// Get Mass Number
			//--------------------------------------------------------------------
			xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("isotope",iso_index).getChildNode("mass_number");
			MassNumber = atoi(xNode.getText());
			//printf("    MassNumber = %d\n",MassNumber);
			//--------------------------------------------------------------------
			// Get mass if present
			//--------------------------------------------------------------------
			Nmass = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("isotope",iso_index).nChildNode("mass");
			if(1 == Nmass){
			    xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("isotope",iso_index).getChildNode("mass");
			    imass = atof(xNode.getText());
			}
			//--------------------------------------------------------------------
			// Get Fraction if present
			//--------------------------------------------------------------------
			Nfraction = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("isotope",iso_index).nChildNode("fraction");
			if(1 == Nfraction){
				xNode = xMainNode.getChildNode("user_isotopes").getChildNode("element",Eindex).getChildNode("isotope",iso_index).getChildNode("fraction");
				ifraction = atof(xNode.getText());
				//printf("    fraction = %f\n",ifraction);
			}
			//--------------------------------------------------------------------
			// Check against default values loaded from AtomTabl.xml
			//--------------------------------------------------------------------
			for(i=0; i<pElements->Element[AtomicNumber].IsotopeTotal; i++){
				Mnumber = pElements->Element[AtomicNumber].Isotope[i]->MassNumber;
				if(Mnumber == MassNumber){
					if( 1 == Nfraction){
						fraction = pElements->Element[AtomicNumber].Isotope[i]->CompositionFraction;
						if( fraction != ifraction ){
							printf("The composition fraction for %s (atomic number %d, mass number %d) has been changed from %f to %f\n",name,AtomicNumber,MassNumber,pElements->Element[AtomicNumber].Isotope[i]->CompositionFraction,ifraction);
							pElements->Element[AtomicNumber].Isotope[i]->CompositionFraction = ifraction;
						}
					}
					if(1 == Nmass){
						mass = pElements->Element[AtomicNumber].Isotope[i]->AtomicMass;
						if( mass != imass ){
							printf("   The mass for %s (atomic number %d, mass number %d) has been changed from %f to %f\n",name,AtomicNumber,MassNumber,pElements->Element[AtomicNumber].Isotope[i]->AtomicMass,imass);
							pElements->Element[AtomicNumber].Isotope[i]->AtomicMass = imass;
						}
					}
					break;
				 }
			 }
		 }
	}

}


void data_write_UserIsotopes(char *path, char *filename, struct element_list *pElements){

	FILE * pFile;
	char *pathfilename;
    int  element_index;
    int  iso_index;


	pathfilename = (char *)malloc((strlen(path)+strlen(filename)+1)*sizeof(char));
	strcpy(pathfilename,path);
	strcat(pathfilename,"\\");
	strcat(pathfilename,filename);
	pFile = fopen (pathfilename,"w");
	if (pFile!=NULL)
	{
		printf("Writing user isotopic composition file: %s\n",filename);
		fprintf(pFile,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		fprintf(pFile,"<!--For isotope entries the <mass_number> must be present. \n");
		fprintf(pFile,"       The <mass> entry may or may not be present. \n");
		fprintf(pFile,"       The <fraction> entry may or may not be present.  -->\n");
		fprintf(pFile,"<user_isotopes>\n");
		for(element_index=1; element_index<ELEMENT_TOTAL; element_index++){
			if( pElements->Element[element_index].AtomicNumber > 0 ){
				fprintf(pFile,"   <element><atomic_number>%d</atomic_number><name>%s</name><symbol>%s</symbol>\n",pElements->Element[element_index].AtomicNumber,pElements->Element[element_index].Name,pElements->Element[element_index].Symbol);
				for(iso_index=0; iso_index<pElements->Element[element_index].IsotopeTotal; iso_index++){
					if( 0 < pElements->Element[element_index].Isotope[iso_index]->MassNumber ){
						fprintf(pFile,"      <isotope><mass_number>%d</mass_number><mass>%20.15f</mass><fraction>%20.15f</fraction></isotope>\n",pElements->Element[element_index].Isotope[iso_index]->MassNumber,pElements->Element[element_index].Isotope[iso_index]->AtomicMass,pElements->Element[element_index].Isotope[iso_index]->CompositionFraction);
					}
				}
				fprintf(pFile,"   </element>\n");
			}
		}	
		fprintf(pFile,"</user_isotopes>\n");
		fclose(pFile);
	}else{
		printf("Error : writing file: %s\n",filename);
	}

}


//----------------------------------------------------------------------
// Create the element list
// Input data_read_NIST(path, pElements) 
//		where path is the directory path to the file NIST_isotopes.txt
//      and pElements is a pointer to the element_list structure
//----------------------------------------------------------------------
void data_read_NIST(char *path, struct element_list *pElements){

  FILE * pFile;
  char *filename;
  char line[100];
  char *pline;
  char *pch,*pch1;
  int  AtomicNumber,AtomicNumber_count[ELEMENT_TOTAL];
  int  MassNumber;
  int  element_index;
  double AtomicMass,AtomicWeight;
  double Composition;

	//---------------------------------------------------
	// Allocate memory for element list and set NULL
	// pointers initially
	//---------------------------------------------------
	for(element_index=0; element_index<ELEMENT_TOTAL; element_index++){
		pElements->Element[element_index].Name                = NULL;
		pElements->Element[element_index].Symbol              = NULL;
		pElements->Element[element_index].AtomicNumber        = 0;
		pElements->Element[element_index].Isotope             = NULL;
		pElements->Element[element_index].IsotopeTotal        = 0;
		pElements->Element[element_index].NonzeroIsotopeIndex = NULL;
		pElements->Element[element_index].NonzeroIsotopeTotal = 0;
		AtomicNumber_count[element_index]                     = 0;
	}



  //---------------------------------------------------
  // open file NIST_isotopes.txt and read in data
  //---------------------------------------------------
  filename = (char *)malloc((strlen(path)+30)*sizeof(char));
  strcpy(filename,path);
  strcat(filename,"\\NIST_isotopes.txt");
  pFile = fopen (filename,"r");
  if (pFile!=NULL)
  {
    printf("Reading file NIST_isotopes.txt\n");
	while (!feof(pFile)) {
		fgets(line , 100 , pFile);
		pline = line;
		if( NULL != line ){
			//printf("line=%s",line);
			//-------------------------------------------------------------
			// Look for "Atomic Number" as this is the first line of entry
			//-------------------------------------------------------------
			pch = strstr(line,"Atomic Number");
			if( NULL != pch ){
				pch=strchr(line,'=');
				AtomicNumber = atoi(pch+1);
				AtomicNumber_count[AtomicNumber] += 1;
				if( 1 == AtomicNumber_count[AtomicNumber] ){
				    pElements->Element[AtomicNumber].AtomicNumber = AtomicNumber;
				}
				//printf("AtomicNumber=%d  count = %d\n",AtomicNumber,AtomicNumber_count[AtomicNumber]);
				//-------------------------------------------------------------
				// Allocate memory for isotope
				//-------------------------------------------------------------
				if( NULL == pElements->Element[AtomicNumber].Isotope ){
					pElements->Element[AtomicNumber].Isotope       = (struct isotope_info **)malloc(sizeof(struct isotope_info *));
					pElements->Element[AtomicNumber].Isotope[0]    = (struct isotope_info *) malloc(sizeof(struct isotope_info));
					pElements->Element[AtomicNumber].IsotopeTotal = 1;
				}else{
					pElements->Element[AtomicNumber].Isotope       = (struct isotope_info **)realloc(pElements->Element[AtomicNumber].Isotope,AtomicNumber_count[AtomicNumber]*sizeof(struct isotope_info *));
					pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1] = (struct isotope_info *) malloc(sizeof(struct isotope_info));
					pElements->Element[AtomicNumber].IsotopeTotal = AtomicNumber_count[AtomicNumber];
				}
				//-------------------------------------------------------------
				// initialize isotope structure
				//-------------------------------------------------------------
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->AtomicMass			= 0;
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->CompositionFraction	= 0;
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->MassNumber			= 0;
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->Name					= NULL;
			}
			//-------------------------------------------------------------
			// Look for "Atomic Symbol" 
			//-------------------------------------------------------------
			pch = strstr(line,"Atomic Symbol");
			if( NULL != pch ){
				if( 1 == AtomicNumber_count[AtomicNumber]){ // only get symbol for first occurrance
					//printf("  s1=%s\n",line);
					pch=strchr(line,'=');					
					pch++;
					while(isspace(*pch)){pch++;}  // go to next non space character
					pch1=pch;
					while(isalpha(*pch1)){pch1++;} // go to next non letter
					*pch1='\0';  // insert null character to end string 
					//printf("  s=%s\n",pch);
					pElements->Element[AtomicNumber].Symbol = (char *)malloc((strlen(pch)+1)*sizeof(char));
					strcpy(pElements->Element[AtomicNumber].Symbol,pch);
					//printf("   symbol=%s\n",pElements->Element[AtomicNumber].Symbol);
				}
			}
			//-------------------------------------------------------------
			// Look for "Mass Number" 
			//-------------------------------------------------------------
			pch = strstr(line,"Mass Number");
			if( NULL != pch ){
				pch=strchr(line,'=');
				MassNumber = atoi(pch+1);
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->MassNumber = MassNumber;
				//printf("    MassNumber=%d\n",pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->MassNumber);
			}
			//-------------------------------------------------------------
			// Look for "Relative Atomic Mass" 
			//-------------------------------------------------------------
			pch = strstr(line,"Relative Atomic Mass");
			if( NULL != pch ){
				pch=strchr(line,'=');
				pch1=strchr(line,'(');  // stop at the first (
				if( NULL != pch1 ){
					*pch1 = '\0';
				}
				AtomicMass = atof(pch+1);
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->AtomicMass = AtomicMass;
				//printf("   AtomicMass=%f\n",pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->AtomicMass);
			}
			//-------------------------------------------------------------
			// Look for "Isotopic Composition" 
			//-------------------------------------------------------------
			pch = strstr(line,"Isotopic Composition");
			if( NULL != pch ){
				pch=strchr(line,'=');
				pch1=strchr(line,'(');  // stop at the first (
				if( NULL != pch1 ){
					*pch1 = '\0';
				}
				Composition = atof(pch+1);
				pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->CompositionFraction = Composition/100;
				//printf("   Composition=%f\n",pElements->Element[AtomicNumber].Isotope[AtomicNumber_count[AtomicNumber]-1]->CompositionFraction);
			}
			//-------------------------------------------------------------
			// Look for "Standard Atomic Weight" 
			//-------------------------------------------------------------
			pch = strstr(line,"Standard Atomic Weight");
			if( NULL != pch ){
				pch=strchr(line,'=');
				pch1=strchr(line,'(');  // stop at the first (
				if( NULL != pch1 ){
					*pch1 = '\0';
				}
				if( 1 == AtomicNumber_count[AtomicNumber]){ // only need to do this once for each element
					AtomicWeight = atof(pch+1);
					pElements->Element[AtomicNumber].AverageMass = AtomicWeight;
					//printf("   AtomicWeight=%f\n",pElements->Element[AtomicNumber].AtomicWeight);
				}
			}
		}
	}
    fclose (pFile);
  }else{
    printf("Error opening file NIST_isotopes.txt\n");
	printf("filename = %s\n",filename);
  }
}


void data_read_AtomTabl(char *path, struct element_list *pElements)
{
	char *filename;
	char symbol[10];
	char Mstring[30];
	char iname[30];
	char isymbol[30];
	char *pstring;
	char *pch;
	int i;
	XMLNode xMainNode,xNode,xNode2;
	XMLNode NodeEntry;
	int Nentries;
	int Nchildren;
	int Eindex;
	size_t Nchar;
	int AtomicNumber;
	int Inumber;
	int iso_index;
	int Nisotopes;
	int  element_index;
	int max_index;
	int nonzero_count;
	double max_fraction;
	double AverageMass;
	int MassNumber;
	int Mnumber;
	double imass,ifraction;
	int foundflag;
	int pElement_isotope_index;
	clock_t time0,time1;

	//---------------------------------------------------
	// open file AtomTabl.xml and read in data
	//---------------------------------------------------
	filename = (char *)malloc((strlen(path)+30)*sizeof(char));
	strcpy(filename,path);
	strcat(filename,"\\AtomTabl.XML");
    printf("Reading file AtomTabl.XML\n");

	//---------------------------------------------------
    // Open and parse the XML file:
	//---------------------------------------------------
	time0 = clock();
    xMainNode=XMLNode::openFileHelper(filename);
	time1 = clock();

	//---------------------------------------------------
	// Count how many Element entries there are
	//---------------------------------------------------
	Nentries = xMainNode.getChildNode("isotope_table").nChildNode("element");
	printf("There are %d element nodes\n", Nentries);
	pElements->Element_Total = Nentries;

	//---------------------------------------------------
    // Read in each entry
	//---------------------------------------------------
	for(Eindex=0; Eindex<Nentries; Eindex++){
		//--------------------------------------------------------------------
		// Get Atomic Number
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.atom_no");
		AtomicNumber = atoi(xNode.getText());
		if( 112< AtomicNumber ){
		    pElements->Element[AtomicNumber].AtomicNumber = AtomicNumber;  // the NIST file only goes to 112
		}
		//printf("Atomic Number = %d %d\n",pElements->Element[AtomicNumber].AtomicNumber,AtomicNumber);
		//--------------------------------------------------------------------
		// Set default values
		//--------------------------------------------------------------------
		pElements->Element[AtomicNumber].Stable = true;
		//--------------------------------------------------------------------
		// Get Element Name
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.name");
		Nchar = strlen(xNode.getText()) + 1;
		pElements->Element[AtomicNumber].Name = (char *)malloc(Nchar*sizeof(char));
		strcpy(pElements->Element[AtomicNumber].Name, xNode.getText());
		//printf("    Name   = %s (%d)\n",pElements->Element[AtomicNumber].Name,AtomicNumber);
		//--------------------------------------------------------------------
		// Get Element Symbol
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.symbol");
		strcpy(symbol,xNode.getText());
		if( (0 == AtomicNumber) || (112 < AtomicNumber)){  // NIST doesn't list the electron as element zero
			Nchar = strlen(symbol) + 1;
			pElements->Element[AtomicNumber].Symbol = (char *)malloc(Nchar*sizeof(char));
			strcpy(pElements->Element[AtomicNumber].Symbol,symbol);
		}
		if( 0 != strcmp(pElements->Element[AtomicNumber].Symbol,symbol)){ // check symbols across files
			printf("   Warning : Symbols don't match for Atomic Number %d  (%s != %s)  Using %s.\n",AtomicNumber,pElements->Element[AtomicNumber].Symbol,symbol,symbol);
			// Elements 110 and 110 don't match so defer to the AtomTabl listing
			strcpy(pElements->Element[AtomicNumber].Symbol,symbol);
		}
		//printf("    Symbol = %s\n",pElements->Element[AtomicNumber].Symbol);
		//--------------------------------------------------------------------
		// Get Average Mass
		//--------------------------------------------------------------------
		xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.avg_mass");
		strcpy(Mstring,xNode.getText());
		pch=strchr(Mstring,'(');  // stop at the first opening parenthesis
		if( NULL != pch ){
			*pch = '\0';
		}
		AverageMass = atof(Mstring);
		if( 0 == AtomicNumber ){
			pElements->Element[AtomicNumber].AverageMass = AverageMass*(-1.0);  // mass of electron 
			//printf("Electron mass = %18.15f\n",pElements->Element[AtomicNumber].AverageMass);
		}
		if( 112< AtomicNumber ){
		    pElements->Element[AtomicNumber].AverageMass = AverageMass;  // the NIST file only goes to 112
		}

		if( (AverageMass != pElements->Element[AtomicNumber].AverageMass) && (0 != AtomicNumber) ){
			printf("   Warning : Average Masses don't match for Atomic Number %d  \n   (%f != %f)  Using AtomTabl mass %f\n",AtomicNumber,pElements->Element[AtomicNumber].AverageMass,AverageMass,AverageMass);
			// Use AtomTabl average mass
			pElements->Element[AtomicNumber].AverageMass = AverageMass;
			//strcpy(pElements->Element[AtomicNumber].Symbol,symbol);
		}
		if(0 == pElements->Element[AtomicNumber].AverageMass){
			pstring=strchr(Mstring,'[');  // if average mass is enclosed in square brackets, it means it is an unstable element and the mass is set to the most stable isotope
			if( NULL != pstring ){
				//---------------------------------
				// we have an unstable element
				//---------------------------------
				pElements->Element[AtomicNumber].Stable = false;
				pch=strchr(Mstring,']');
				*pch = '\0';  // get rid of the closing square bracket
				pstring++; //get rid of starting square bracket
				Inumber = atoi(pstring);
				//printf("    Unstable Element  Mass=%f  Isonumber=%d  %s\n",pElements->Element[AtomicNumber].AverageMass,Inumber,Mstring);
				//printf("    Isotope Count = %d\n",pElements->Element[AtomicNumber].IsotopeTotal);
				foundflag = 0;
				for(iso_index=0; iso_index<pElements->Element[AtomicNumber].IsotopeTotal; iso_index++){
					if( pElements->Element[AtomicNumber].Isotope[iso_index]->MassNumber == Inumber ){
						pElements->Element[AtomicNumber].AverageMass = pElements->Element[AtomicNumber].Isotope[iso_index]->AtomicMass;
						pElements->Element[AtomicNumber].Isotope[iso_index]->CompositionFraction = 1.0;
						foundflag = 1;
						break;
					}
				}
				if( 0 == foundflag ){ // the MassNumber wasn't found so just pick one
					pElements->Element[AtomicNumber].AverageMass = pElements->Element[AtomicNumber].Isotope[0]->AtomicMass;
					pElements->Element[AtomicNumber].Isotope[0]->CompositionFraction = 1.0;
				}
				//printf("    **********  Mass=%f\n",pElements->Element[AtomicNumber].AverageMass);
			}
		}
		//--------------------------------------------------------------------
		// Get AtomTabl.xml Isotope info
		//--------------------------------------------------------------------
		Nisotopes = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).nChildNode("element.isotope");
		for(iso_index=0; iso_index<Nisotopes; iso_index++){
			//--------------------------------------------------------------------
			// Get Mass Number
			//--------------------------------------------------------------------
			xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).getChildNode("element.isotope.mass_no");
			MassNumber = atoi(xNode.getText());
			//--------------------------------------------------------------------
			// Get pElement isotope index
			//--------------------------------------------------------------------
			for(i=0; i<pElements->Element[AtomicNumber].IsotopeTotal; i++){
				Mnumber = pElements->Element[AtomicNumber].Isotope[i]->MassNumber;
				if(Mnumber == MassNumber){
					pElement_isotope_index = i;
					break;
				}
			 }
			//--------------------------------------------------------------------
			// Get Isotope Name
			//--------------------------------------------------------------------
			Nchildren = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).nChildNode("element.isotope.name");
			if( 0 < Nchildren ){
				xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).getChildNode("element.isotope.name");
				strcpy(iname,xNode.getText());
			}else{
				strcpy(iname,"None");
			}
			pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->Name = (char *)malloc((strlen(iname)+1)*sizeof(char));
			strcpy(pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->Name,iname);
			//--------------------------------------------------------------------
			// Get Isotope Symbol
			//--------------------------------------------------------------------
			Nchildren = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).nChildNode("element.isotope.symbol");
			if( 0 < Nchildren ){
				xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).getChildNode("element.isotope.symbol");
				strcpy(isymbol,xNode.getText());
			}else{
				strcpy(isymbol,"None");
			}
			pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->Symbol = (char *)malloc((strlen(isymbol)+1)*sizeof(char));
			strcpy(pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->Symbol,isymbol);
			//--------------------------------------------------------------------
			// Get Isotope Mass
			//--------------------------------------------------------------------
			xNode = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).getChildNode("element.isotope.mass");
			imass = atof(xNode.getText());
			//--------------------------------------------------------------------
			// Get Isotope Fraction
			//--------------------------------------------------------------------
			xNode     = xMainNode.getChildNode("isotope_table").getChildNode("element",Eindex).getChildNode("element.isotope",iso_index).getChildNode("element.isotope.fract");
			ifraction = atof(xNode.getText());

			//printf("         %d %s %s %f %f\n",MassNumber,iname,isymbol,imass,ifraction);
			//--------------------------------------------------------------------
			// Compare to NIST and default to AtomTabl.xml
			//--------------------------------------------------------------------
			if( pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->AtomicMass != imass ){
				printf("             mass difference:     %f != %f  (%d %d)\n",pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->AtomicMass,imass,AtomicNumber,MassNumber);
				pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->AtomicMass = imass;
			}
			if(  pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->CompositionFraction != ifraction ){
				printf("             fraction difference: %f != %f  (%d %d)\n",pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->CompositionFraction,ifraction,AtomicNumber,MassNumber);
				pElements->Element[AtomicNumber].Isotope[pElement_isotope_index]->CompositionFraction = ifraction;
			}
		}
	}
	//---------------------------------------------------------------
	// Find isotope with greatest composition percentage
	//---------------------------------------------------------------
	for(element_index=1; element_index<ELEMENT_TOTAL; element_index++){
		max_index = 0;
		max_fraction = 0;
		for(iso_index=0; iso_index<pElements->Element[element_index].IsotopeTotal; iso_index++){
			if( pElements->Element[element_index].Isotope[iso_index]->CompositionFraction > max_fraction ){
				max_fraction = pElements->Element[element_index].Isotope[iso_index]->CompositionFraction;
				max_index      = iso_index;
			}
		}
		pElements->Element[element_index].MostCommonIsotopeIndex = max_index;
	}
	//---------------------------------------------------------------------
	// Find indicies of the isotopes with nonzero composition percentages
	//---------------------------------------------------------------------
     for(element_index=1; element_index<ELEMENT_TOTAL; element_index++){
		nonzero_count = 0;
		for(iso_index=0; iso_index<pElements->Element[element_index].IsotopeTotal; iso_index++){
			if( 0 < pElements->Element[element_index].Isotope[iso_index]->CompositionFraction ){
				nonzero_count += 1;
			}
		}
		if( nonzero_count == 0 ){
			printf("Element %d has no isotopes\n",element_index);
		}
		pElements->Element[element_index].NonzeroIsotopeTotal    = nonzero_count;
		pElements->Element[element_index].NonzeroIsotopeIndex = (int *)malloc(nonzero_count*sizeof(int));
		nonzero_count = 0;
		for(iso_index=0; iso_index<pElements->Element[element_index].IsotopeTotal; iso_index++){
			if( 0 < pElements->Element[element_index].Isotope[iso_index]->CompositionFraction ){
				pElements->Element[element_index].NonzeroIsotopeIndex[nonzero_count] = iso_index;
				nonzero_count += 1;
			}
		}
	}
}


