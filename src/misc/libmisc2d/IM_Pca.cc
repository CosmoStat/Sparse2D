/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12
**    
**    File:  IM_Pca.cc
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
#include "IM_Pca.h"

/************************************************************/

CorrelMatAna2D::CorrelMatAna2D(int Nima)
{
     Nbr_Image = Nima;
     MatCor.alloc(Nima);
     TabMean.alloc(Nima);
     CorrelMat.alloc(Nbr_Image, Nbr_Image);
     TabMean.init();
}

/************************************************************/

void CorrelMatAna2D::print()
{
    for(int k=0; k < Nbr_Image; k++)
         cout  << "image " << k+1 <<  ": mean = " << TabMean(k) << endl;
    MatCor.print();
}

/************************************************************/

void CorrelMatAna2D::subtract_mean(Ifloat *TabIma)
{
    int i,k;
    int Nl = TabIma[0].nl();
    int Nc = TabIma[0].nc();
            
    for(k=0; k < Nbr_Image; k++)
    {
       TabMean(k) = average (TabIma[k]);   
       for (i=0; i < Nl*Nc; i++) TabIma[k](i) -= TabMean(k);
    }
}

/************************************************************/

void CorrelMatAna2D::add_mean(Ifloat *TabIma)
{
    int i,k;
    int Nl = TabIma[0].nl();
    int Nc = TabIma[0].nc();
    for(k=0; k < Nbr_Image; k++)
    {
        for (i=0; i < Nl*Nc; i++) TabIma[k](i) += TabMean(k);
    }
}

/************************************************************/

void CorrelMatAna2D::compute(Ifloat *TabIma, float ValMin)
{
    int i,j;
        
    // compute the correlation matrix
    CorrelMat.init();
        
    // diagonal elements are equal to 1
    for(i=0; i < Nbr_Image; i++) CorrelMat(i,i) = 1.;
       
    for(i=0; i < Nbr_Image; i++)
    for(j=i+1; j < Nbr_Image; j++)
    { 
       if (ValMin > 0)  
 	       CorrelMat(i,j) = correlation(TabIma[i], TabIma[j], ValMin);
       else CorrelMat(i,j) = correlation(TabIma[i], TabIma[j]);
 	  
       // the matrix is symetric  
       CorrelMat(j,i) =  CorrelMat(i,j);
    }
       
       // compute the eigen values
       MatCor.compute_eigen(CorrelMat);
}      

/************************************************************/

void CorrelMatAna2D::transform (Ifloat *TabImaIn, Ifloat *TabImaOut)
{
   int i,j,k;
   fltarray Vect(Nbr_Image);
   fltarray Vectout(Nbr_Image);     
     
   for (i=0; i < TabImaIn[0].nl(); i++)
   for (j=0; j < TabImaIn[0].nc(); j++)
   {
      for (k = 0; k < Nbr_Image; k++) Vect(k) = (TabImaIn[k])(i,j);
      MatCor.transform(Vect, Vectout);
      for (k = 0; k < Nbr_Image; k++) (TabImaOut[k])(i,j) = Vectout(k);
   } 
}

/************************************************************/

void CorrelMatAna2D::invtransform ( Ifloat *TabImaIn, Ifloat *TabImaOut)
{
   int i,j,k;
   fltarray Vect(Nbr_Image);
   fltarray Vectout(Nbr_Image);
    
    for (i=0; i < TabImaIn[0].nl(); i++)
    for (j=0; j < TabImaIn[0].nc(); j++)
    {
       for (k = 0; k < Nbr_Image; k++) Vect(k) = TabImaIn[k](i,j);
       MatCor.invtransform(Vect, Vectout);
       for (k = 0; k < Nbr_Image; k++) TabImaOut[k](i,j) = Vectout(k);
    }   
}

/************************************************************/

void CorrelMatAna2D::invsubtransform (Ifloat *TabImaIn, Ifloat *TabImaOut, 
                         int NbrEigen, int *TabKill)
{
    int i,j,k;
    fltarray Vect(Nbr_Image);
    fltarray Vectout(Nbr_Image);
    int N = NbrEigen;
    
    if ((N < 1) || (N > Nbr_Image)) N = Nbr_Image;
    
    for (i=0; i < TabImaIn[0].nl(); i++)
    for (j=0; j < TabImaIn[0].nc(); j++)
    {
        for (k = 0; k < N; k++)
        {
          if ((TabKill != NULL) && (TabKill[k] == 1)) Vect(k) = 0.;
          else Vect(k) = TabImaIn[k](i,j);
        }
        MatCor.invtransform(Vect, Vectout, N);
        for (k = 0; k < Nbr_Image; k++) TabImaOut[k](i,j) = Vectout(k);
    }   
}

/************************************************************/








