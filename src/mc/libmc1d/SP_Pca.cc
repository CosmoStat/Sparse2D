/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: %V%
**
**    Author: 08/23/99
**
**    Date:  1.1
**    
**    File:  SP_Pca.cc
**
*******************************************************************************
**
**    DESCRIPTION  class for principal component analysis for images  
**    ----------- 
**                 
******************************************************************************/
 
#include "IM_Obj.h"
// #include "IM_Errors.h"
#include "NR.h"
#include "CPca.h"
#include "SP_Pca.h"

/************************************************************/

CorrelMatAna1D::CorrelMatAna1D(int NbSpectre) {

     _NbrSpectre = NbSpectre;
     _MatCor.alloc(NbSpectre);
     _TabMean.alloc(NbSpectre);
     _CorrelMat.alloc(NbSpectre, NbSpectre);
     _TabMean.init();
}

/************************************************************/

void CorrelMatAna1D::print() {

    for(int k=0; k<_NbrSpectre; k++)
         cout  << "spectre " << k+1 <<  ": mean = " << _TabMean(k) << endl;
    _MatCor.print();
}

/************************************************************/

void CorrelMatAna1D::subtract_mean(fltarray *TabSpectre) {

    int i,k;
    int NbPoint = TabSpectre[0].n_elem();
        
    for(k=0; k<_NbrSpectre ; k++) {
       _TabMean(k) = TabSpectre[k].mean();   
       for (i=0; i<NbPoint; i++) TabSpectre[k](i) -= _TabMean(k);
    }
}

/************************************************************/

void CorrelMatAna1D::add_mean(fltarray *TabSpectre) {

    int i,k;
    int NbPoint = TabSpectre[0].n_elem();

    for(k=0; k<_NbrSpectre ; k++) {
        for (i=0; i<NbPoint; i++) TabSpectre[k](i) += _TabMean(k);
    }
}

/************************************************************/

void CorrelMatAna1D::compute(fltarray *TabSpectre, float ValMin) {

    int i,j;
        
    // compute the correlation matrix
    _CorrelMat.init();
        
    // diagonal elements are equal to 1
    for(i=0; i<_NbrSpectre ; i++) _CorrelMat(i,i) = 1.;
       
    for(i=0; i<_NbrSpectre ; i++)
    for(j=i+1; j<_NbrSpectre ; j++)
    { 
       if (ValMin > 0)  
 	    _CorrelMat(i,j) = correlation(TabSpectre[i], TabSpectre[j], ValMin);
       else _CorrelMat(i,j) = correlation(TabSpectre[i], TabSpectre[j]);
 	  
       // the matrix is symetric  
       _CorrelMat(j,i) =  _CorrelMat(i,j);
    }
       
       // compute the eigen values
       _MatCor.compute_eigen(_CorrelMat);
}      

/************************************************************/

void CorrelMatAna1D::transform (fltarray *TabSpectreIn, fltarray *TabSpectreOut){

   int i,k;    
   int NbPoint = TabSpectreIn[0].n_elem();

   fltarray Vect(_NbrSpectre);
   fltarray Vectout(_NbrSpectre);     
     
   for (i=0; i<NbPoint ; i++) {
      for (k=0; k<_NbrSpectre; k++) Vect(k) = (TabSpectreIn[k])(i);
      _MatCor.transform(Vect, Vectout);
      for (k=0; k<_NbrSpectre; k++) (TabSpectreOut[k])(i) = Vectout(k);
   } 
}

/************************************************************/

void CorrelMatAna1D::invtransform (fltarray *TabSpectreIn, 
                                   fltarray *TabSpectreOut) {
   int i,k;   
   int NbPoint = TabSpectreIn[0].n_elem();

   fltarray Vect(_NbrSpectre);
   fltarray Vectout(_NbrSpectre);
    
    for (i=0; i<NbPoint; i++) {
       for (k=0; k<_NbrSpectre; k++) Vect(k) = TabSpectreIn[k](i);
       _MatCor.invtransform(Vect, Vectout);
       for (k=0; k<_NbrSpectre; k++) TabSpectreOut[k](i) = Vectout(k);
    }   
}

/************************************************************/

void CorrelMatAna1D::invsubtransform (fltarray *TabSpectreIn, 
                                      fltarray *TabSpectreOut, 
                                      int NbrEigen, int *TabKill) {
				       
    int i,k;   
    int NbPoint = TabSpectreIn[0].n_elem();

    fltarray Vect(_NbrSpectre);
    fltarray Vectout(_NbrSpectre);
    int N = NbrEigen;
    
    if ((N < 1) || (N > _NbrSpectre)) N = _NbrSpectre;
    
    for (i=0; i<NbPoint; i++) {
        for (k=0; k<N; k++) {
          if ((TabKill != NULL) && (TabKill[k] == 1)) Vect(k) = 0.;
          else Vect(k) = TabSpectreIn[k](i);
        }
        _MatCor.invtransform(Vect, Vectout, N);
        for (k=0; k<_NbrSpectre; k++) TabSpectreOut[k](i) = Vectout(k);
    }   
}

/************************************************************/



float correlation (const fltarray &Spec1, const fltarray &Spec2, float ValMin) {
    double Sum_X2,Sum_Y2,Sum_XY;
    int i;
    int NbPoint = Spec1.n_elem();
    double Coef;
    float ValRet;

   // test spectre size
   if (Spec1.n_elem() != Spec2.n_elem())
   {
      cerr << "Error in correllation routine: spectres have different sizes ..." << endl ; 
      cerr << "   spectre 1: " <<  Spec1.n_elem() << endl ;
      cerr << "   spectre 2: " <<  Spec2.n_elem() << endl ;
      exit(-1);
   }
       
    Sum_X2 = Sum_Y2 = Sum_XY = 0.;
    for (i = 0; i<NbPoint; i++) {
        if (ABS (Spec1 (i)) > ValMin) {
           Sum_X2 += Spec1 (i) * Spec1 (i);
           Sum_Y2 += Spec2 (i) * Spec2 (i);
           Sum_XY += Spec1 (i) * Spec2 (i);
        }
    }
    Coef = sqrt (Sum_X2 * Sum_Y2);
    if (Coef > 0.) ValRet = (float)(Sum_XY / Coef);
    else ValRet = 0;
    return (ValRet);
}




