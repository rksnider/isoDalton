/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  data.h                                                  */
/*               Header file for data.cpp, which contains source code    */ 
/*               for reading all the data files                          */
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

// The following is defined to suppress warnings about depracated functions strcpy and strcat etc 
// Warnings such as "This function or variable may be unsafe. Consider using strcat_s instead..." 
#ifdef WIN32
	#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <string.h>
#include <stdio.h>
#include <conio.h>  // for _kbhit()
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "isotopes.h"  // isotopes data structures
#include "xmlParser.h"


void data_read_RESID(char *, struct RESID_info *);
void data_read_NIST(char *, struct element_list *);
void data_read_AtomTabl(char *, struct element_list *);
void data_write_UserIsotopes(char *, char *, struct element_list *);
void data_read_UserIsotopes(char *, char *, struct element_list *);
void data_normalize_fractions(struct element_list *);



