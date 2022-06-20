/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/10/2001 
**    
**    File:  IM3D_Radon.h
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

#ifndef _RADON3D_H_
#define _RADON3D_H_

#include "Array.h"
#include "IM_Obj.h"
#include "FFTN_3D.h"
#include "IM3D_IO.h"

/***********************************************************************/

#define NBR_RADON3D 1
enum type_radon3d {RADON3D_FFT};

#define DEF_RADON3D  RADON3D_FFT

inline const char * StringRadon3d (type_radon3d type)
{
   switch (type)
   {
      case RADON3D_FFT: 
              return ("Radon transform and reconstruction in Fourier space"); 
   }
   return NULL;
}

class Radon3D {
   FFTN_3D CFFT3D;     // FFT Class
   Bool Alloc;
   int Nx,Ny,Nz;
   int ProjectNbr2D; // Number of projections on a plane
   int ProjectNbr;   // Number of project
   int ResolNbr;     // Number of pixels per projection (resolution)
   type_radon3d RadonMethod; // Radon method to use
   void make_table();  // make sin and cos table
   void reset_param(); // reinitialize the field
   double *SinTab, *CosTab; // table used for direct and backprojection
   double *TabAngle; // table used for transform and
   double *TabStep;  // reconstruction in Fourier space
   double *TabAngle2;
   float NormRadonFFT; // Normalization parameter (only used with RADON_FFT)
 
   void fft_trans(fltarray &Cube, Ifloat &Radon3DImage);
   void fft_rec(Ifloat &Radon3DImage, fltarray &cube);

  public:
   Radon3D() {SinTab=(double*)NULL;CosTab=(double*)NULL;
              TabAngle=(double*)NULL;TabStep=(double*)NULL;
              TabAngle2=(double*)NULL;TabStep=(double*)NULL;
	      reset_param();}
   Bool Verbose;
   type_radon3d radon_method() { return RadonMethod;}
   int radon_nl()   { return ProjectNbr;}
   int radon_nc()   { return ResolNbr;}
   int cube_nx()   { return Nx;}
   int cube_ny()   { return Ny;}
   int cube_nz()   { return Nz;}
   // void num_project(int p, int & u, int & v, int & w, float &Angle1, float &Angle2);
   void alloc(int Nx_Dat, int Ny_Dat, int Nz_Dat, type_radon3d RM=DEF_RADON3D);
   // allocate the image dimension 
 
   void transform(fltarray &Cube, Ifloat &RadonImage);
   // Radon transform of a cube
   //              RadonImage  should be first allocated before calling
   //              if it is not alloc is called with default parameter
   //              i.e. number of projection = cube.nx*2
   //                   resolution per projection = Image.nx

   void recons(Ifloat &RadonImage, fltarray &CubeRec);
   // Radon inverse transform of an image
   // if the allocation has not been done, it calls the alloc routine
   // with the following parameter:
   //      Nx_Cube = RadonImage.nc
   //      Ny_Cube = RadonImage.nc
   //      Nz_Cube = RadonImage.nc
   //      RM      =  DEF_RADON  
   //      RadonNl = RadonImage.nl
   //      RadonNc = RadonImage.nc

   void fft_trans_cf(fltarray &Cube,  cfarray &  Data2);
   void invfft_line(cfarray & Data2, Ifloat &RadonImage);
 
   void fft_line(Ifloat &RadonImage, cfarray &Data1);
   void fft_rec_cf(cfarray &Data1, fltarray &Cube);
   ~Radon3D(){reset_param();}
};

#endif
