/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  98/02/03
**    
**    File:  PCA.cc
**
*******************************************************************************
**
**    DESCRIPTION  do a principal component analysis on images  
**    ----------- 
**
***************************************************************************/

#ifndef _IMPCA_H_
#define _IMPCA_H_

const int MAX_NBR_PCA_IMA = 50;

class CorrelMatAna2D {
 
public:
   int Nbr_Image;        // Number of images to analyse
   CPCA MatCor;          // class for principal component analysis
   fltarray CorrelMat;   // correlation matrix (2D)
   fltarray TabMean;     // mean value of each image (1D array)
   
   CorrelMatAna2D(int Nima); // class initialixation

  void print(); // print the result of the analysis
                // mean value per image
                // correlation matrix
                // eigen vectors
                // eigen values
  
   void subtract_mean(Ifloat *TabIma);
   // calculate the mean of each image in TabIma,
   // store it in TabMean and subtract it
   
   void add_mean(Ifloat *TabIma);
   // add to each image k the mean value TabMean[k]
   
   // calulate the correlation matrix and diogonalize it
   // find the eigen values and the eigen vectors   
   void compute(Ifloat *TabIma, float ValMin=0);
         
  void transform (Ifloat *TabImaIn, Ifloat *TabImaOut);
  // transform an image array into its principal component
  // (TabImaOut can be set to the same pointer as TabImaIn, then
  // (the input data will be overwritten, but memory will be safed).
  
  void  invtransform ( Ifloat *TabImaIn, Ifloat *TabImaOut);
  // apply an inverse transform: the images are reconstructed from
  // the eigen vector images.
  // (TabImaOut can be set to the same pointer as TabImaIn, then
  // (the input data will be overwritten, but memory will be safed).
  
  void  invsubtransform (Ifloat *TabImaIn, Ifloat *TabImaOut, 
                         int NbrEigen, int *TabKill=NULL);
  // apply an inverse transform: the images are reconstructed from
  // a sub set of the eigen vector images.
  
  ~CorrelMatAna2D(){Nbr_Image=0;};
};
          
#endif
