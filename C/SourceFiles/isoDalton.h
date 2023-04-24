/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  isoDalton.h                                             */
/*               Header file for isoDalton.cpp                           */
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

struct istates_info {
	int StateTotal;
	double *mass;
	double *prob;
};


void isoDalton_get_isotopes(char *, char *, char *, struct element_list *);
void isoDalton_parse_molecular_formula(char *, struct molecule_info *, struct element_list *);
void isoDalton_combine_masses(int* , double *, double *, int);
void isoDalton_exact_mass(struct molecule_info *, struct element_list *, int, struct istates_info *, int);


