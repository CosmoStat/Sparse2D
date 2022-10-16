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
**    Date:  20/03/01
**    
**    File:  RID_STF.h
**
*******************************************************************************
**
**    DESCRIPTION   
**    -----------   
**                 
******************************************************************************/

#ifndef _RID_STF_H_
#define _RID_STF_H_

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "FFTN_1D.h"
#include "Ridgelet.h"

class STF_RID {
   int Np;              // Input data size
   int Nl,Nc;           // STF size
   int NlRid,NcRid;     // Ridgelet size 
   Icomplex_f STF_Buff; // Buffer used for STF transform
   Icomplex_f Freq0;    // Zero frequency buffer 
   int RidBlockSize;    // Ridgelet block size
   public:
   ST_FFTN STF;         // Short term Fourier transform class
   Ridgelet RG;         // Ridgelet class
   Ifloat DataSpec;     // Buffer for the spectrogram calculation
   Bool Verbose;
   STF_RID() {Verbose=False;RidBlockSize=Nl=Nc=Np=NlRid=NcRid=0;}
   void alloc(int Nx, type_std_win TW, float WinParam, int WindowSize, 
              int Step, int Nbr_Plan, Bool BlockOverlap, int BlockSize);
         // Nx = IN: number of pixel of the input 1D data
	 // TW = IN: type of window to be used for the Short Term Fourier 
	 //          Transform
	 // WinParam = IN: parameter of the window
	 // WindowSize = IN: window size
	 // Step = IN: temporal step between two Fourier analysis
	 // Nbr_Plan = IN: number of scale used in the ridgelet transform
	 // BlockOverlap = IN: block overlapping in the ridgelet transform
	 // BlockSize= IN: block size in the ridgelet transform

   void transform(fltarray &Data, Ifloat & Ridtrans);
         // Apply the Chirplet Transform = STF + Rdigelet transform
	 // Data = IN: input 1D data
	 // Ridtrans= OUT: output 2D data
	 // The spectrogram is also saved in DataSpec
	 
   void recons(Ifloat & Ridtrans, fltarray &DataRec);
         // Apply the inverse Chirplet Transform 
	 // Ridtrans= IN: input 2D data
	 // Data = OUT: output 1D data
	 // The reconstructed spectrogram is also saved in DataSpec
	 
   void thresholding(Ifloat & ImaRid, Ifloat & ImaSigma, float N_Sigma);
        // Threshold the Chirplet Transform using a noise map

   void thresholding(Ifloat & ImaRid, Ifloat & Support);
        // Threshold the Chirplet Transform using a support.
	// if (Support(i,j) == 0) then ImaRid(i,j) = 0

   void transform_cf(fltarray &Data, Ifloat & TransDat_re, Ifloat & TransDat_im);
   void recons_cf(Ifloat & TransDat_re, Ifloat & TransDat_im, fltarray &DataRec);


};

#endif

