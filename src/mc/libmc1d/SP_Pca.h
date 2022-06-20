/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: %V%
**
**    Author: 10/30/02
**
**    Date:  2.1
**    
**    File:  SP_Pca.h
**
*******************************************************************************
**
**    DESCRIPTION  do a principal component analysis on images  
**    ----------- 
**
***************************************************************************/

#ifndef _SPPCA_H_
#define _SPPCA_H_

const int MAX_NBR_PCA_IMA = 50;

class CorrelMatAna1D {
 
public:
   int      _NbrSpectre;  // Number of spectres to analyse
   CPCA     _MatCor;      // class for principal component analysis
   fltarray _CorrelMat;   // correlation matrix (2D)
   fltarray _TabMean;     // mean value of each image (1D array)
   
   CorrelMatAna1D(int NbSpectre); // class initialixation

   void print(); // print the result of the analysis
                 // mean value per spectre
                 // correlation matrix
                 // eigen vectors
                 // eigen values
  
   void subtract_mean(fltarray *TabSpectre);
   // calculate the mean of each image in TabSpectre,
   // store it in TabMean and subtract it
   
   void add_mean(fltarray *TabSpectre);
   // add to each spectre k the mean value TabMean[k]
   
   // calulate the correlation matrix and diogonalize it
   // find the eigen values and the eigen vectors   
   void compute(fltarray *TabSpectre, float ValMin=0);
         
   void transform (fltarray *TabSpectreIn, fltarray *TabSpectreOut);
   // transform an spectre array into its principal component
   // (TabSpectreOut can be set to the same pointer as TabSpectreIn, then
   // (the input data will be overwritten, but memory will be safed).
  
   void  invtransform (fltarray *TabSpectreIn, fltarray *TabSpectreOut);
   // apply an inverse transform: the images are reconstructed from
   // the eigen vector images.
   // (TabImaOut can be set to the same pointer as TabImaIn, then
   // (the input data will be overwritten, but memory will be safed).
  
   void  invsubtransform (fltarray *TabSpectreIn, fltarray *TabSpectreOut, 
                         int NbrEigen, int *TabKill=NULL);
   // apply an inverse transform: the images are reconstructed from
   // a sub set of the eigen vector images.
  
  ~CorrelMatAna1D(){_NbrSpectre=0;};
};
 
float correlation (const fltarray &Spec1, const fltarray &Spec2, float ValMin=0.0);         


#endif
