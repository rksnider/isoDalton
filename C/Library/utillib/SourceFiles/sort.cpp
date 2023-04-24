/* SPDX-License-Identifier: GPL-2.0 or MIT                               */
/*-----------------------------------------------------------------------*/
/* Description:  sort.cpp                                                */
/*               Sorting code, adapted from numerical recipies in C code.*/
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

#include "sort.h"
#include <stdio.h>

//-----------------------------------------------------
// heap sort  ascending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
//-----------------------------------------------------
void heapsort_dbl_up(int Nelements, double *array1) {
	int i,j,k,m;
	double tmp1;
	double *array1offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
		}else{
			tmp1      = array1offset[m];
			array1offset[m] = array1offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] < array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 < array1offset[j]) {
				array1offset[i] = array1offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
	}
}

//-----------------------------------------------------
// heap sort (Two vectors) ascending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
// This is version 2 where the sorting is determined by
// array1 and the second array (array2) gets the
// same rearrangement.
//-----------------------------------------------------
void heapsort_2dbl_up(int Nelements, double *array1, double *array2) {
	int i,j,k,m;
	double tmp1,tmp2;
	double *array1offset;
	double *array2offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 
	array2offset = array2-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
			tmp2 = array2offset[k];
		}else{
			tmp1      = array1offset[m];
			tmp2      = array2offset[m];
			array1offset[m] = array1offset[1];
			array2offset[m] = array2offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				array2offset[1] = tmp2;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] < array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 < array1offset[j]) {
				array1offset[i] = array1offset[j];
				array2offset[i] = array2offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
		array2offset[i] = tmp2;
	}
}

//-----------------------------------------------------
// heap sort (Two vectors) decending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
// This is version 2 where the sorting is determined by
// array1 and the second array (array2) gets the
// same rearrangement.
//-----------------------------------------------------
void heapsort_2dbl_down(int Nelements, double *array1, double *array2) {
	int i,j,k,m;
	double tmp1,tmp2;
	double *array1offset;
	double *array2offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 
	array2offset = array2-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
			tmp2 = array2offset[k];
		}else{
			tmp1      = array1offset[m];
			tmp2      = array2offset[m];
			array1offset[m] = array1offset[1];
			array2offset[m] = array2offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				array2offset[1] = tmp2;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] > array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 > array1offset[j]) {
				array1offset[i] = array1offset[j];
				array2offset[i] = array2offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
		array2offset[i] = tmp2;
	}
}


//-----------------------------------------------------
// heap sort (Three vectors) ascending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
// This is version 2 where the sorting is determined by
// array1 and the second and thirde arrays get the
// same rearrangement.
//-----------------------------------------------------
void heapsort_3dbl_up(int Nelements, double *array1, double *array2, double *array3) {
	int i,j,k,m;
	double tmp1,tmp2,tmp3;
	double *array1offset;
	double *array2offset;
	double *array3offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 
	array2offset = array2-1; 
	array3offset = array3-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
			tmp2 = array2offset[k];
			tmp3 = array3offset[k];
		}else{
			tmp1      = array1offset[m];
			tmp2      = array2offset[m];
			tmp3      = array3offset[m];
			array1offset[m] = array1offset[1];
			array2offset[m] = array2offset[1];
			array3offset[m] = array3offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				array2offset[1] = tmp2;
				array3offset[1] = tmp3;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] < array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 < array1offset[j]) {
				array1offset[i] = array1offset[j];
				array2offset[i] = array2offset[j];
				array3offset[i] = array3offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
		array2offset[i] = tmp2;
		array3offset[i] = tmp3;
	}
}

//-----------------------------------------------------
// heap sort (Three vectors) ascending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
// This is version 2 where the sorting is determined by
// array1 and the second and thirde arrays get the
// same rearrangement.
//-----------------------------------------------------
void heapsort_3int_up(int Nelements, int *array1, int *array2, int *array3) {
	int i,j,k,m;
	int tmp1,tmp2,tmp3;
	int *array1offset;
	int *array2offset;
	int *array3offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 
	array2offset = array2-1; 
	array3offset = array3-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
			tmp2 = array2offset[k];
			tmp3 = array3offset[k];
		}else{
			tmp1      = array1offset[m];
			tmp2      = array2offset[m];
			tmp3      = array3offset[m];
			array1offset[m] = array1offset[1];
			array2offset[m] = array2offset[1];
			array3offset[m] = array3offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				array2offset[1] = tmp2;
				array3offset[1] = tmp3;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] < array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 < array1offset[j]) {
				array1offset[i] = array1offset[j];
				array2offset[i] = array2offset[j];
				array3offset[i] = array3offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
		array2offset[i] = tmp2;
		array3offset[i] = tmp3;
	}
}

//-----------------------------------------------------
// heap sort (Three vectors) ascending
// This has a O(Nlog2N) search time and the
// worst case is also O(Nlog2N) unlike quicksort
// This is version 2 where the sorting is determined by
// array1 and the second and thirde arrays get the
// same rearrangement.
//-----------------------------------------------------
void heapsort_1dbl_2int_up(int Nelements, double *array1, int *array2, int *array3) {
	int i,j,k,m;
	double tmp1;
	int tmp2,tmp3;
	double *array1offset;
	int *array2offset;
	int *array3offset;

	// This came from Numerical routines which expects the vectors to have offset of array1[1..Nelements] 
	// so we need to shift the array by 1 so when it references array[1], it is really array[0].
	array1offset = array1-1; 
	array2offset = array2-1; 
	array3offset = array3-1; 

	m = Nelements;
	k = (Nelements >> 1) + 1;
	while(1) {
		if (k > 1) {
			tmp1 = array1offset[--k];
			tmp2 = array2offset[k];
			tmp3 = array3offset[k];
		}else{
			tmp1      = array1offset[m];
			tmp2      = array2offset[m];
			tmp3      = array3offset[m];
			array1offset[m] = array1offset[1];
			array2offset[m] = array2offset[1];
			array3offset[m] = array3offset[1];
			if (--m == 1) {
				array1offset[1] = tmp1;
				array2offset[1] = tmp2;
				array3offset[1] = tmp3;
				return;
			}
		}
		i = k;
		j = k << 1;
		while (j <= m) {
			if (j < m && array1offset[j] < array1offset[j+1]){ 
				++j; 
			}
			if (tmp1 < array1offset[j]) {
				array1offset[i] = array1offset[j];
				array2offset[i] = array2offset[j];
				array3offset[i] = array3offset[j];
				j += (i = j);
			}else{
				j = m + 1;
			}
		}
		array1offset[i] = tmp1;
		array2offset[i] = tmp2;
		array3offset[i] = tmp3;
	}
}
