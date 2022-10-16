/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  IM_Radon.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _RADON_H_
#define _RADON_H_

#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "FastSlantStack.h"

/***********************************************************************/

#define NBR_ALL_RADON 5
#define NBR_RADON 4
enum type_radon {RADON_PROJECT_BACK,RADON_PROJECT_FFT,RADON_FFT, RADON_FINITE, RADON_FSS, RADON_FFT_2};

#define DEF_RADON  RADON_FFT
#define DEF_RADON_FILTER_WIDTH  100
#define DEF_RADON_SIGMA  10.

inline const char * StringRadon (type_radon type)
{
   switch (type)
   {
      case RADON_FFT_2:
	      return ("Slant Stack Radon transform in Fourier domain");
      case RADON_PROJECT_BACK:
	      return ("Radon transform (resp. backprojection) in spatial domain"); 
      case RADON_PROJECT_FFT: 
              return ("Radon transform in spatial domain and reconstruction in Fourier domain "); 
      case RADON_FFT: 
              return ("Radon transform (resp. reconstruction) in Fourier space (Linogram)"); 
      case RADON_FINITE: 
              return ("Finite Radon Transform"); 
      case RADON_FSS:
              return ("Slant Stack Radon transform"); 
   }
   return NULL;
}

class Radon {
   Bool Alloc;
   int Nl,Nc;
   int ProjectNbr;   // Number of project
   int ResolNbr;     // Number of pixels per projection (resolution)
   type_radon RadonMethod; // Radon method to use
   void make_table();  // make sin and cos table
   void reset_param(); // reinitialize the field
   double *SinTab, *CosTab; // table used for direct and backprojection
   double *TabAngle; // table used for transform and
   double *TabStep;  // reconstruction in Fourier space
   float NormRadonFFT; // Normalization parameter (only used with RADON_FFT)
   void project(Ifloat &Data1, Ifloat &RadonTrans);
   // Radon transform from project in the direct space
   // Data1 = input image
   // RadonTrans = output transform

   void back_project(Ifloat &RadonTrans, Ifloat &ImagRec);
   // Reconstruction by backprojection
   // RadonTrans = input Radon transform
   // ImagRec= output image reconstruction

   void invfft_project(Ifloat & RadonImage, Ifloat &ImageOut, Bool Permut=True);
   // inverse radon transform, for a radon calculated by the project method

   void fft_trans(Ifloat &Data1, Ifloat &Data2);
   void fft_rec(Ifloat &Data1, Ifloat &Data2, Bool Permut=False);
   void mr_io_fill_header(fitsfile *fptr);

  public:
   int NbrIter; // Number of iteration for the Slant Stack Radon reconstruction;
                // Default is 1.
   slant_stack_radon SSR;
   type_format FormatInputImag;   // data format of the input image
   fft_slant_stack_radon FSSR;
   Radon() {extern type_format Format_Imag;
            SinTab=NULL;CosTab=NULL;TabAngle=NULL;TabStep=NULL;reset_param();
            FormatInputImag = Format_Imag;NbrIter=DEF_RADON_FSS_NBR_ITER;}
   Bool Verbose;
   int FilterWidth;   // width of the filter 
   float FilterSigma; // standard deviation of the filter
   type_radon radon_method() { return RadonMethod;}
   int radon_nl()   { return ProjectNbr;}
   int radon_nc()   { return ResolNbr;}
   int imag_nl()   { return Nl;}
   int imag_nc()   { return Nc;}

   void alloc(int Nl_Imag, int Nc_Imag, type_radon RM=DEF_RADON, int RadonNl=0, int RadonNc=0);
   // allocate the image dimension 
   // per default, we fixe the number of project equal to the image number of lines
   //              and the resolution number to the number of column

   void filter(Ifloat &RadonImage);
   // Make a filtering of the Radon image in the Fourier space
   // This step is usually perform before backprojection.

   void transform(Ifloat &Image, Ifloat &RadonImage);
   // Radon transform of an image
   //              RadonImage  should be first allocated before calling
   //              if it is not alloc is called with default parameter
   //              i.e. number of projection = Image.nl*2
   //                   resolution per projection = Image.nc

   void recons(Ifloat &RadonImage, Ifloat &ImageRec);
   // Radon inverse transform of an image
   // if the allocation has not been done, it calls the alloc routine
   // with the following parameter:
   //      Nl_Imag = RadonImage.nc
   //      Nc_Imag = RadonImage.nc
   //      RM      =  DEF_RADON  
   //      RadonNl = RadonImage.nl
   //      RadonNc = RadonImage.nc

   void fft_trans_cf(Ifloat &Image, Icomplex_f & Data2);
   void invfft_line(Icomplex_f & Data2, Ifloat &RadonImage);
   void fft_line(Ifloat &RadonImage, Icomplex_f & Data1, Bool Permu=True);
   void fft_rec_cf(Icomplex_f & Data1, Ifloat &Image);
   void finite_radon(Ifloat &Image, Ifloat & Radon);
   void inv_finite_radon(Ifloat & Radon, Ifloat &Image);
   ~Radon(){reset_param();}

    void read  (char *Name, Ifloat &Ima);
    // Read a radon transform and initialize the class
    
    void write (char *Name, Ifloat &Ima);
    // write a radon transform in the FITS format
};


#endif
