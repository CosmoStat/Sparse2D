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
**    Date:  04/09/99
**    
**    File:  ColRest.cc
**
*******************************************************************************
**
**    DESCRIPTION   Class for color image restoration
**    -----------                 
**
******************************************************************************/

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_Sigma.h"
#include "ColRest.h"
#include "FFTN_2D.h"
#include "IM_Math.h"
#include "MR_Contrast.h"

/****************************************************************************/

// static void histo_equal(fltarray &Data, int Band, int NbrBin=1024);
// static void histo_equal(fltarray &Data, int Band, int NbrBin)
// {
//     Contrast CO;
//     Ifloat ImaFrame;
//     float *Ptr = Data.buffer();
//     int Nx = Data.nx();
//     int Ny = Data.ny();
//     Ptr += Band*Nx*Ny;
//     ImaFrame.alloc(Ptr, Ny, Nx);
//     CO.histo_equal(Data, NbrBin);
// }

/*****************************************************/

void ColorRestore::rescale_0_255(fltarray &Data)
{
   int b,i,j;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   float MinT = Data.min();
   float MaxT =  Data.max();
   float ScaleT = 255. / (MaxT - MinT);
   
   if (Verbose == True)
     cout << "DMin = " << Data.min() << " DMax = " << Data.max() << endl;        

   for (b = 0; b < Nz; b++) 
   for (i=0; i < Ny; i++)
   for (j=0; j < Nx; j++) 
   {
      float Val = Data(j,i ,b);
      Data(j,i,b) = (Val - MinT)*ScaleT;
   }
   // cout << "Min = " << MinT << " Max = " << MaxT  << endl; 
}
   
/***************************************************************************/

void ColorRestore::sature_clipping(fltarray &Data)
{
   int i,j,k;
   int Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   float MeanClip, SigmaClip;
   
   Data.sigma_clip (MeanClip, SigmaClip, 3);
   float MinVal = MeanClip - ClipVal*SigmaClip;
   float MaxVal = MeanClip + ClipVal*SigmaClip;
   for (k = 0; k < Nz; k++)
   {
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) 
      {
          if (Data(j,i,k) > MaxVal) Data(j,i,k) = MaxVal;
          else if (Data(j,i,k) < MinVal ) Data(j,i,k) = MinVal;
      }
   }
}

/***************************************************************************/

void ColorRestore::sature_clipping(Ifloat &Data)
{
   int i,j;
   int Nc = Data.nc();
   int Nl = Data.nl();
   float MeanClip, SigmaClip;
   
   sigma_clip (Data, MeanClip, SigmaClip, 3);
   float MinVal = MeanClip - ClipVal*SigmaClip;
   float MaxVal = MeanClip + ClipVal*SigmaClip;
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) 
   {
      if (Data(i,j) > MaxVal) Data(i,j) = MaxVal;
      else if (Data(i,j) < MinVal ) Data(i,j) = MinVal;
   }
}

/***************************************************************************/

void ColorRestore::info(fltarray &Data)
{
   int Nx = Data.nx();
   int Ny = Data.ny();
   // int Nz = Data.nz();
   Ifloat ImaFrame;
   float *Ptr = Data.buffer();
   Ptr = Data.buffer();
   ImaFrame.alloc(Ptr, Ny, Nx);
   INFO_X(ImaFrame,"Red");
   Ptr += Nx*Ny;  
   ImaFrame.alloc(Ptr, Ny, Nx);
   INFO_X(ImaFrame,"Green");
   Ptr += Nx*Ny;  
   ImaFrame.alloc(Ptr, Ny, Nx);
   INFO_X(ImaFrame,"Blue");
}

/***************************************************************************/

// void ColorRestore::retinex(fltarray &Data)
// { 
//    int i,j,k;
//    int Nc = Data.nx();
//    int Nl = Data.ny();
//    int Nz = Data.nz();
//    fltarray LuvDat(Nc,Nl,Nz);
//  
//  
//    
//    sature_clipping(Data);
// }

/***************************************************************************/

void ColorRestore::retinex(fltarray &Data)
{ 
   int i,j,k;
   int Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   float Sigma = 80.;
   Ifloat Gaussian(Nl,Nc,"gauss");
   fltarray Ret(Nc,Nl,Nz);
   FFTN_2D CFFT;
      
   float *PtrData = Data.buffer();
   float *PtrRet = Ret.buffer();
   Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
 
   for (k = 0; k < 3; k++)
   {
      ImaFrame.alloc(PtrData, Nl, Nc);
      ImaRet.alloc(PtrRet, Nl, Nc);
      Gaussian.init();
      make_gaussian(Gaussian, (float) (Sigma/sqrt(2.)));
      CFFT.convolve(ImaFrame, Gaussian, Buff);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) 
      {
           if ((ImaFrame(i,j) > 0) && (Buff(i,j) > 0))
	        ImaRet(i,j) = log(ImaFrame(i,j)) - log( Buff(i,j));
      }
      PtrData += Nl*Nc;
      PtrRet += Nl*Nc;
   }
   Data = Ret;
}
   
/***************************************************************************/

void ColorRestore::multiscale_retinex(fltarray &Data)
{ 
   int i,j,b,k;
   int Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   int N = 3;
   float TabSigma[3] = {15.,80.,250.};
   // float TabSigma[3] = {50.,2.,2.};
   float G = 192;
   float B =  -30.;
   float Alpha = 125.;
   float Beta = 46.;
   Ifloat Gaussian(Nl,Nc,"gauss");
   fltarray Ret(Nc,Nl,Nz);
   fltarray Ci(Nc,Nl,Nz);
   FFTN_2D CFFT;

   // Ci calculation
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) 
   {
      float Sum=0.;
      for (k = 0; k < 3; k++) Sum += Data(j,i,k);
      for (k = 0; k < 3; k++) 
         if (Data(j,i,k) > 0) Ci(j,i,k) = Beta*(log(Alpha*Data(j,i,k)/Sum));
   }
   
   
      
   float *PtrData = Data.buffer();
   float *PtrRet = Ret.buffer();
   Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
 
   for (k = 0; k < 3; k++)
   {
      ImaFrame.alloc(PtrData, Nl, Nc);
      ImaRet.alloc(PtrRet, Nl, Nc);
      for (b = 0; b < N; b++)
      {
         Gaussian.init();
         make_gaussian(Gaussian, (float) (TabSigma[b]/sqrt(2.)));
         CFFT.convolve(ImaFrame, Gaussian, Buff);
         for (i = 0; i < Nl; i++)
         for (j = 0; j < Nc; j++) 
	 {
            if ((ImaFrame(i,j) > 0) && (Buff(i,j) > 0))
 	     ImaRet(i,j) += 
 	        G*(Ci(j,i,k)*(log(ImaFrame(i,j)) - log( Buff(i,j))) + B) / (float) N;
            // ImaRet(i,j) += (log(ImaFrame(i,j)) - log( Buff(i,j))) / (float) N;
         }
      }
      PtrData += Nl*Nc;
      PtrRet += Nl*Nc;
   }
   Data = Ret; 
   
//    for (k = 0; k < 3; k++)
//    {
//       for (i = 0; i < Nl; i++)
//       for (j = 0; j < Nc; j++) Data(j,i,k) = Ret(j,i,k);
//    }
}
   
/***************************************************************************/

void ColorRestore::equalize_histo_luminance(fltarray &Data)
{ 
   Contrast CO;
   
   int i,j,Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   fltarray LuvDat(Nc,Nl,Nz);
   Ifloat KHisto(Nl,Nc,"Khisto");
   LuvDat = Data;
   rgb_to_luv(LuvDat);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++) KHisto(i,j) = LuvDat(j,i,0);
   CO.histo_equal(KHisto);
 
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
      float Scale = (LuvDat(j,i,0) != 0) ? (KHisto(i,j) / LuvDat(j,i,0)): 1.;
      Data(j,i,0) *= Scale;
      Data(j,i,1) *= Scale;
      Data(j,i,2) *= Scale;
   }  
}

/***************************************************************************/

void ColorRestore::equalize_histo(fltarray &Data, int BandNumber)
{ 
   Contrast CO;
   int Nx = Data.nx();
   int Ny = Data.ny();
   // int Nz = Data.nz();
   Ifloat KHisto;
   float *Ptr = Data.buffer() + Nx*Ny*BandNumber;
   KHisto.alloc(Ptr, Ny, Nx);
   rgb_to_hsv(Data);
   CO.histo_equal(KHisto);
   hsv_to_rgb(Data);
}

/***************************************************************************/

void ColorRestore::col_filter_band (Ifloat & Band, float  Noise, float NSigma)
{
    int Nl = Band.nl();
    int Nc = Band.nc();
    int i,j;
    float Level= NSigma* Noise;

    if (HardThreshold == True)
    {
         for (i=0;i < Nl; i++)
         for (j=0;j < Nc; j++)
                 if (ABS(Band(i,j)) < Level) Band(i,j) = 0.;
    }
    else
    {
       Level /= 2;
       for (i=0;i < Nl; i++)
       for (j=0;j < Nc; j++)
       if (ABS(Band(i,j)) < Level) Band(i,j) = 0.;
       else
       {
	       if (Band(i,j) > 0) Band(i,j) -= Level;
	       else Band(i,j) += Level;
       }
    }
    // Kill isolated pixel  
    for (i = 1; i < Nl-1; i++)
    for (j = 1; j < Nc-1; j++)
       if ((Band(i,j) != 0) && (
           ((Band(i-1,j) == 0)  && (Band(i+1,j) == 0)  &&
           (Band(i,j+1) == 0)  && (Band(i,j-1) == 0))))  Band(i,j) = 0;
}

/***************************************************************************/

void ColorRestore::col_filter(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima)
{
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   int b,i,j,k;
   Ifloat *TabBand;
   int Nbr_UndecScale = (NbrUndec < 0) ? Nbr_Plan: NbrUndec;
   FilterAnaSynt SelectFilter(F_MALLAT_7_9);
   SubBandFilter SB1D(SelectFilter, NORM_L2);
   SB1D.Border = I_MIRROR;
   Bool *TabDec = new Bool [Nbr_Plan]; 
   for (i = 0; i < MIN(Nbr_UndecScale, Nbr_Plan); i++) TabDec[i] = False;
   for (i = Nbr_UndecScale; i < Nbr_Plan; i++) TabDec[i] = True;  
   HALF_DECIMATED_2D_WT WT(SB1D);
   int NbrBand = WT.alloc(TabBand, Ny, Nx, Nbr_Plan, TabDec);   
   float SigmaNoise=0.;
   Ifloat Frame;
   Bool PositivContraint = True;
   float *Ptr = Data.buffer();
   rgb_to_yuv(Data);
    
    for (b = 0; b < Nz; b++)
    {
       Frame.alloc(Ptr, Ny, Nx);
       if (b != Nz-1) Ptr += Ny*Nx;
       if (MAD == False)
       {
          if (Noise_Ima < FLOAT_EPSILON) 
                SigmaNoise = detect_noise_from_med (Frame);
          else  SigmaNoise = Noise_Ima;
        
          if (Verbose == True) 
          {     
            cout << "YUV band " << b+1;
            cout << ": Nl = " << Frame.nl() << " Nc = " << Frame.nc() << endl;
            cout << "  Estimated noise standard deviation = " <<  SigmaNoise   << endl;
          }
       }
       if (Verbose == True) cout << "   WT transform: " << Nbr_UndecScale <<   endl;
       WT.transform(Frame, TabBand, Nbr_Plan, TabDec);
       if (Verbose == True) cout << "   Filtering ... " << endl;
       for (k=0; k < NbrBand-1; k++) 
       {       
          int N1 = (TabBand[k]).nl();
	  int N2 = (TabBand[k]).nc();
          if (MAD == True)
 	  {
              fltarray Buff(N1*N2);
              int Ind = 0;             
	      for (i=0;i<N1;i++)
              for (j=0;j<N2;j++) Buff(Ind++) = ABS((TabBand[k])(i,j));
              SigmaNoise = get_median(Buff.buffer(), N1*N2) / 0.6745;
	  }       
	  if (Verbose == True)
	  {
	     float T = SigmaNoise*N_Sigma;
	     if (HardThreshold == False) T /= 2.;
	     cout << "Band " << k+1 << " Nl = " << N1 << " Nc = " << N2 << " T = " << T << endl;
	  }
	  // if (b < 3) col_filter_band(TabBand[k], SigmaNoise, N_Sigma+1);
	  // else 
	     col_filter_band(TabBand[k], SigmaNoise, N_Sigma);
       }
       
       // Kill isolated pixel
       
       if (Verbose == True) cout << "   Reconstruction ... " << endl;
       WT.recons(TabBand, Frame, Nbr_Plan, TabDec);
       for (i=0; i < Ny; i++)
       for (j=0; j < Nx; j++) Data(j,i,b) =  Frame(i,j);
    }
    yuv_to_rgb(Data);
    
    if (PositivContraint == True)
    {
      for (b = 0; b < Nz; b++)
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++)
        if (Data(j,i,b) < 0) Data(j,i,b) = 0.;
    }
    WT.free(TabBand, Nbr_Plan);
    delete [] TabDec;
}

/***************************************************************************/

// void ColorRestore::atrou_retinex(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima)
// {
//    int i,j,b,k;
//    int Nc = Data.nx();
//    int Nl = Data.ny();
//    int Nz = Data.nz();
//    type_transform Transform = TO_PAVE_BSPLINE; // TO_UNDECIMATED_MALLAT; 
//    MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
//    fltarray Ret(Nc,Nl,Nz);
//    float SigmaNoise, *PtrData = Data.buffer();
//    float *PtrRet = Ret.buffer();
//    Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
// 
//    cout << "Ret " << endl;
//    Ret = Data;
//    rgb_to_yuv(Ret);
//  
//    ImaRet.alloc(PtrRet, Nl, Nc);
//    MR_Data.transform(ImaRet);
//     
//    PtrRet = Ret.buffer();
//    ImaFrame.alloc(PtrRet, Nl, Nc);
//    Buff = ImaFrame;
//    if (Noise_Ima < FLOAT_EPSILON) 
//              SigmaNoise = detect_noise_from_med (ImaFrame);
//    else  SigmaNoise = Noise_Ima;
//           
//    for (b = MR_Data.nbr_band()-2; b >= 0; b--)
//    {
//       // float Level = SigmaNoise*MR_Data.band_norm(b)*N_Sigma;
//       for (i = 0; i < Nl; i++)
//       for (j = 0; j < Nc; j++) 
//       {
//      	 // float W = ABS(MR_Data(b,i,j)) / Level;
// 	 // if (W > 1) W = 1;
// 	 // else W=0;
// 	 
//          if (MR_Data(b,i,j) >= 0)
//               MR_Data(b,i,j) =   log(1. + MR_Data(b,i,j));
// 	 else MR_Data(b,i,j) = - log(1. + ABS(MR_Data(b,i,j)));
//       }
//    }
//    MR_Data.recons(ImaFrame);
//    for (i=0; i < Nl; i++)
//    for (j=0; j < Nc; j++)
//    {
//       float Scale = (Buff(i,j) != 0) ? (ImaFrame(i,j) / Buff(i,j)): 1.;
//       Data(j,i,0) *= Scale;
//       Data(j,i,1) *= Scale;
//       Data(j,i,2) *= Scale;
//    }  
// }

/***************************************************************************/

void ColorRestore::atrou_retinex(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima)
{
   int i,j,b,k;
   int Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   type_transform Transform = TO_PAVE_BSPLINE;
   MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
   fltarray Ret(Nc,Nl,Nz);
   float *PtrData = Data.buffer();
   float *PtrRet = Ret.buffer();
   Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
   for (k = 0; k < 3; k++)
   {
      ImaFrame.alloc(PtrData, Nl, Nc);
      ImaRet.alloc(PtrRet, Nl, Nc);
      MR_Data.transform(ImaFrame);
      Buff = MR_Data.band(Nbr_Plan-1);
      for (b = Nbr_Plan-2; b >= 0; b--)
      {
         for (i = 0; i < Nl; i++)
         for (j = 0; j < Nc; j++) 
	 {
            if ((ImaFrame(i,j) > 0) && (Buff(i,j) > 0))
              ImaRet(i,j) += (log(ImaFrame(i,j)) - log( Buff(i,j))) / (float) (Nbr_Plan-1);
         }
	 Buff += MR_Data.band(b);
      }
      PtrData += Nl*Nc;
      PtrRet += Nl*Nc;
   }
   Data = Ret;    
}

/***************************************************************************/

void ColorRestore::atrou_log(fltarray  &Data, int Nbr_Plan, 
                             Bool Sign, float Noise_Ima, float LogParam)
{
   int i,j,b,k;
   int Nc = Data.nx();
   int Nl = Data.ny();
   int Nz = Data.nz();
   type_transform Transform = TO_PAVE_BSPLINE;
   MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
   float *PtrData = Data.buffer();
   Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
   fltarray BuffData(Nc,Nl,Nz);
   
   BuffData = Data;
   rgb_to_luv(Data);
   //Nz=1;
   cout << "atrou_log: Noise_Ima = " << Noise_Ima << endl;
   for (k = 0; k < 1; k++)
   {
      ImaFrame.alloc(PtrData, Nl, Nc);
      MR_Data.transform(ImaFrame);
      ImaFrame.init();
      for (b = Nbr_Plan-2; b >= 0; b--)
      {
         float Noise = Noise_Ima*MR_Data.band_norm(b);
         if (Noise < FLOAT_EPSILON) Noise = 1.;
         for (i = 0; i < Nl; i++)
         for (j = 0; j < Nc; j++)  
         {
            float Coef = MR_Data(b,i,j);
	    float Gamma = contrast_function(ABS(Coef));
            MR_Data(b,i,j) *= Gamma;
	 }
      }
      PtrData += Nl*Nc;
      MR_Data.recons(Buff);
   }

   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
      if (ImaFrame(i,i) != 0) Buff(i,j) /= ImaFrame(i,i); 
	 
   for (k = 0; k < Nz; k++)
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
   {
      Data(j,i,k) = BuffData(j,i,k)*Buff(i,j);
   }
   //yuv_to_rgb(Data);
}

/*************************************************************************/
// produce poor result ==> not used
// void ColorRestore::atrou_clipping(fltarray &Data, int Nbr_Plan)
// {
//    int b,k;
//    int Nc = Data.nx();
//    int Nl = Data.ny();
//    int Nz = Data.nz();
//    type_transform Transform = TO_PAVE_BSPLINE;
//    cout << "ATROUS CLIP" << endl;
//    rgb_to_luv(Data);
//    MultiResol MR_Data(Nl, Nc, Nbr_Plan, Transform, "diadic");
//    fltarray Ret(Nc,Nl,Nz);
//    float *PtrData = Data.buffer();
//    Ifloat ImaFrame, ImaRet, Buff(Nl,Nc, "buff");
//    
//    for (k = 0; k < 3; k++)
//    {
//       ImaFrame.alloc(PtrData, Nl, Nc);
//       MR_Data.transform(ImaFrame);
//       Buff = MR_Data.band(Nbr_Plan-1);
//       for (b = Nbr_Plan-1; b >= 0; b--)
//       {   
//          sature_clipping(MR_Data.band(b));
// //          for (i = 0; i < Nl; i++)
// //          for (j = 0; j < Nc; j++) 
// // 	 {
// //             if ((ImaFrame(i,j) > 0) && (Buff(i,j) > 0))
// //               ImaRet(i,j) += (log(ImaFrame(i,j)) - log( Buff(i,j))) / (float) (Nbr_Plan-1);
// //          }
// // 	 Buff += MR_Data.band(b);
//       }
//       MR_Data.recons(ImaFrame);
//       PtrData += Nl*Nc;
//    }
//    luv_to_rgb(Data);
//    // Data = Ret;    
// }

/***************************************************************************/

void ColorRestore::enhance(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima)
{
   int Nx = Data.nx();
   int Ny = Data.ny();
   // int Nz = Data.nz();
   int b,k,i,j;
   float MeanU, MeanV, MeanTU, MeanTV;
   fltarray TabMin(3), TabMax(3);
   fltarray TabTMin(3), TabTMax(3);
   
   // RGB to LUV transformation
    rgb_to_luv(Data);
    
    type_transform Transform = TO_DIADIC_MALLAT;
    MultiResol MR_Data(Ny, Nx, Nbr_Plan, Transform, "diadic");
    MultiResol MR_Data1(Ny, Nx, Nbr_Plan, Transform, "diadic");
    MultiResol MR_Data2(Ny, Nx, Nbr_Plan, Transform, "diadic");
    double Q100 = pow((double) 100., Contrast_Q_Param);

    // Non linear mapping of the luminance
    if (Contrast_Q_Param != 0)
    for (i=0; i < Ny; i++)
    for (j=0; j < Nx; j++)
 	 Data(j,i,0) = (float) (pow((double) Data(j,i,0), 1.-Contrast_Q_Param)*Q100);
 
    // Wavelet Transform
    Ifloat ImaFrame;
    float Gamma,  *Ptr = Data.buffer();
    ImaFrame.alloc(Ptr, Ny, Nx);
    TabMin(0) = min(ImaFrame);
    TabMax(0) = max(ImaFrame);
    if (Verbose == True) cout <<"Transform Band 1 " <<endl;
    MR_Data.transform(ImaFrame);
    Ptr += Nx*Ny;
    if (Verbose == True) cout <<"Transform Band 2 " <<endl;
    ImaFrame.alloc(Ptr, Ny, Nx);
    MeanU = average(ImaFrame);
    TabMin(1) = min(ImaFrame);
    TabMax(1) = max(ImaFrame);
    MR_Data1.transform(ImaFrame);
    Ptr += Nx*Ny;
    if (Verbose == True) cout <<"Transform Band 3 " <<endl;
    ImaFrame.alloc(Ptr, Ny, Nx);
    MeanV = average(ImaFrame);
    TabMin(2) = min(ImaFrame);
    TabMax(2) = max(ImaFrame);
    MeanV = average(ImaFrame);
    MR_Data2.transform(ImaFrame);
    if (Verbose == True) cout <<"Contrast function calculation " <<endl;
    
    // Wavelet coefficients correction
    for (k=0; k < MR_Data.nbr_scale()-1; k++) 
    {   
         // Ifloat Frame(Ny,Nx,"Frame");    
         // Gradient calculation  
	if (Verbose == True) cout <<"Scale " << k+1 << endl;
	for (i=0; i < Ny; i++)
        for (j=0; j < Nx; j++) 
	{
 	   Gamma = sqr(MR_Data(2*k,i,j)/4.)+ sqr(MR_Data(2*k+1,i,j)/4.);
	   Gamma += sqr(MR_Data1(2*k,i,j)/4.)+ sqr(MR_Data1(2*k+1,i,j)/4.);
	   Gamma += sqr(MR_Data2(2*k,i,j)/4.)+ sqr(MR_Data2(2*k+1,i,j)/4.);
	   Gamma = sqrt(Gamma);
 	   Gamma = contrast_function(Gamma);
	   // Frame(i,j) = Gamma;
	   // Gamma = 1;
	   MR_Data(2*k,i,j) *= Gamma;
	   MR_Data(2*k+1,i,j) *= Gamma;
	   MR_Data1(2*k,i,j) *= Gamma;
	   MR_Data1(2*k+1,i,j) *= Gamma;
	   MR_Data2(2*k,i,j) *= Gamma;
	   MR_Data2(2*k+1,i,j) *= Gamma;
       }
       // io_write_ima_float("xx_frame.fits", Frame);
    }
    
    // Wavelet Reconstruction
    Ptr = Data.buffer();
    ImaFrame.alloc(Ptr, Ny, Nx);
    if (Verbose == True) cout <<"Recons Band 1 " <<endl;
    MR_Data.recons(ImaFrame);
    Ptr += Nx*Ny;
    ImaFrame.alloc(Ptr, Ny, Nx);
    TabTMin(0) = min(ImaFrame);
    TabTMax(0) = max(ImaFrame);
    MR_Data1.recons(ImaFrame);
    TabTMin(1) = min(ImaFrame);
    TabTMax(1) = max(ImaFrame);
    Ptr += Nx*Ny;
    ImaFrame.alloc(Ptr, Ny, Nx);
    MR_Data2.recons(ImaFrame);
    TabTMin(2) = min(ImaFrame);
    TabTMax(2) = max(ImaFrame);
    
    // Rescale L,u,v
    for (b=0; b < 3; b++)
    for (i=0; i < Ny; i++)
    for (j=0; j < Nx; j++)
    {
       float Scale = (TabMax(b) - TabMin(b)) / (TabTMax(b) - TabTMin(b));
       Data(j,i,b) =  (Data(j,i,b) - TabTMin(b)) * Scale + TabMin(b);
    }
     
    // inverse mapping of the Luminance
    if (Contrast_Q_Param != 0)
    {
       double  Q1 = 1. / (1-Contrast_Q_Param);
       Q100 = pow((double) 100., -Contrast_Q_Param*Q1);
       for (i=0; i < Ny; i++)
       for (j=0; j < Nx; j++)
 	  Data(j,i,0) = (float) (pow((double) Data(j,i,0), Q1)*Q100);
    }
    
    // Shift the component u and v
    Ptr = Data.buffer() + Nx*Ny;
    ImaFrame.alloc(Ptr, Ny, Nx);
    MeanTU = average(ImaFrame);
    Ptr += Nx*Ny;
    ImaFrame.alloc(Ptr, Ny, Nx);
    MeanTV = average(ImaFrame);
    for (i=0; i < Ny; i++)
    for (j=0; j < Nx; j++)
    {
       Data(j,i,1) += MeanU - MeanTU;
       Data(j,i,2) += MeanV - MeanTV;
    }
    
    // LUV to RGB
    luv_to_rgb(Data);
}

/***************************************************************************/

