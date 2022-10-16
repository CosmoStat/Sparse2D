/*******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Jean-Luc Starck
 **
 **    Date: 19/09/08
 **    
 **    File:  MatMask.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Class for manipulating mask on the sphere   
 **    ----------- 
 **                 
 ******************************************************************************/

#ifndef _MATMASK_H_
#define _MATMASK_H_


#include "HealpixClass.h"

#define  DEF_NBR_MASTER_ITER 20

arr<double> mrs_wigner3j2( double il1, double il2, arr<double> il3, arr<double> lnwpro );
dblarray mrs_mmake_mll( arr<double> ell, arr<double> well, long lmax );

dblarray mrs_matmask( Hmap<double> Mask, int & lmax );
// Return the matrix associated to the mask retalive to the Healpix map Mask

void matmask_mult_cl(PowSpec & PS, dblarray & Mat, PowSpec & PSout, bool transpose=false);
// Multiply the powerspectrum PS with the Matrix Mat. 
// if transpose == true, then multiply the powerspectrum PS with the transpose of the matrix with 
void matmask_mult_tabcl(dblarray & PS, dblarray & Mat, dblarray & PSout, bool transpose=false);

void iter_master_decconv(PowSpec & ClData, PowSpec & ClSol, Hdmap & MaskMap, dblarray & MatMask, int Niter, int Lmax_Tot, bool Verbose=false);


#endif

