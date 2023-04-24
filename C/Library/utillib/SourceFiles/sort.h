/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  sort.h                                                  */
/*               Header code for sort.cpp                                */
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

void heapsort_dbl_up(int , double *);
void heapsort_2dbl_up(int , double *, double *);
void heapsort_2dbl_down(int , double *, double *);
void heapsort_3dbl_up(int , double *, double *, double *);
void heapsort_3int_up(int , int *, int *, int *);
void heapsort_1dbl_2int_up(int , double *, int *, int *);
