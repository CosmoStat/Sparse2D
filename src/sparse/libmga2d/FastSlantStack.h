/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/09/02 
**    
**    File:  FastSlantStack.h
**
*******************************************************************************
**
**    DESCRIPTION:   Fast Slant Stack radon transform of an image
**    ----------- 
**
******************************************************************************/

#ifndef __IM_FSS__
#define __IM_FSS__

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_Line.h"
#include "FFTN_1D.h"
#include "FFTN_2D.h"

#define DEF_RADON_FSS_NBR_ITER  10

/**********************************************************************/

class slant_stack_radon 
{
   FFTN_2D FFT2D; // FTT2D class
   FFTN_1D FFT1D; // FTT1D class     
   int N;  // number of lines and columns of the input image
   int Nr; // number of lines and columns of the radon transform
   cfarray TF1DDataX; // Buffer used for the Radon transform
   cfarray TF1DDataY; // ditto
   cfarray FTLineX;   // ditto
   cfarray FTLineY;   // ditto
   public:
     Bool Verbose;
     Bool UseLeviCode;
     Bool OptIter;   // Applying an optimal step gradient method for the iterative reconstruction
     Bool Positive; // If true, the iterative reconstruction enforces the
                  // positivity of the solution. Default is false.
     void alloc(int N);  // Allocate the CLASS for an image of size NxN
     int ima_nl() { return N;} // return the number of lines of the input
                               // image
     int ima_nc() { return N;} // return the number of columns of the input
                               // image
     int nlr () { return Nr;}  // return the number of lines of radon 
                               // transform
     int ncr () { return Nr;}  // return the number of columns of radon 
                               // transform
     slant_stack_radon(){Verbose=False;UseLeviCode=True;OptIter=True;}
     ~slant_stack_radon(){}
     void transform(Ifloat &Data, Ifloat & Radon); // Radon transform
     void inverse(Ifloat & Radon, Ifloat &Data); // Backprojection
     void inverse_fft(Ifloat & Radon, Ifloat &Data);
           // Backprojection in the Fourier domain
     void iter_inverse(Ifloat & Radon, Ifloat &Data, int Niter=DEF_RADON_FSS_NBR_ITER);
           // Iterative Backprojection
};

/**********************************************************************/
// Fast Slant Stack radon transform of an image in the Fourier domain
/**********************************************************************/

class fft_slant_stack_radon 
{
   FFTN_2D FFT2D;
   FFTN_1D FFT1D;      
   int Nl,Nc;   // number of lines and columns of the input image
   int Nl1,Nc1; // number of lines and columns of the first buffer
                // Nl1 = 2Nl and Nc1 = Nc
   int Nl2,Nc2; // number of lines and columns of the second buffer
                // Nl2 = Nl and Nc2 = 2Nc
   int Nlr,Ncr; // number of lines and columns of the radon transform
                // Nlr= 2Nl and Ncr = 2Nc
   int Nl12, Nc12;
   int Nl22, Nc22;
   Icomplex_f TF_Buff1; // Buffer used for the Radon transform
   Icomplex_f TF_Buff2; // ditto
   Icomplex_f TF_Rec;   // ditto
   intarray TX0,TY0;   // ditto
   intarray TX1,TY1;   // ditto
   intarray TX2,TY2;   // ditto
   Ifloat Buff1, Buff2;  // ditto

   void fft_radon(Icomplex_f & TFBuff1,  Icomplex_f & TFBuff2, Ifloat & Radon);
   void inv_fft_radon(Ifloat & Radon, Icomplex_f & TFBuff);
   public:
     Bool Positive; // If true, the iterative reconstruction enforces the
                  // positivity of the solution
     void alloc(int N);   // Allocate the CLASS for an image of size NxN
                          // WARNING: N must be an odd number
     int ima_nl() { return Nl;}// return the number of lines of the input
                               // image
     int ima_nc() { return Nc;}// return the number of columns of the input
                               // image
     int nlr () { return Nlr;}// return the number of lines of radon 
                               // transform
     int ncr () { return Ncr;}// return the number of columns of radon 
                               // transform
     fft_slant_stack_radon(){}
     ~fft_slant_stack_radon(){}
     void transform(Ifloat &Data, Ifloat & Radon);// Radon transform
     void inverse(Ifloat & Radon, Ifloat &Data);// inverse Radon transform
     void iter_inverse(Ifloat & Radon, Ifloat &Data, int Niter=10);
     // iterative inverse Radon transform
};

 
#endif
