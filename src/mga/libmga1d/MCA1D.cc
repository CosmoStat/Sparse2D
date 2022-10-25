/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Hubert Druesne, Philippe Querre
**
**    Date:  02/06/2003 
**    
**    File:  MCA1D.cc
**
*******************************************************************************
**
**    DESCRIPTION  Decomposition of a signal on multiple bases
**    ----------- 
**                 
******************************************************************************/


#include "GlobalInc.h"
 #include "Usage.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Sigma.h"
#include "MR1D_Obj.h"
#include "MR_SoftRegul.h" 
#include "IM1D_Dct.h"
#include "IM1D_ALDCT.h"
#include "MDCT1D.h"
#include "MCA1D.h" 
#include <vector>
#include <algorithm>
 
/*********************************************************************/
 
inline float MCA1D::mca1d_update (float CoefSol, float Threshold, 
                                  float SoftLevel) 
{
// L2 norm minimization with L1 norm penalt
   
   float NewCoef = hard_threshold(CoefSol, Threshold);
   if (UseNormL1 == True) {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return NewCoef; 
}

/***********************************************************************/

inline bool MCA1D::mca1d_is_on_support (float CoefSol, float Threshold, 
                                        float SoftLevel) 
{
   float NewCoef = hard_threshold(CoefSol, Threshold);
   if (UseNormL1 == True) 
   {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return (!(NewCoef == 0));  
} 
  
/***********************************************************************/

void MCA1D::reconstruction (fltarray& Result) 
{
   int b;
   Result = TabSigRec[0];
   for (b=1; b < NbrBase; b++)  Result += TabSigRec[b];
   // if (PositivRecSig == True) threshold(Result);
}

/****************************************************************************/

void MCA1D::make_residual (fltarray& Data, fltarray& Residual) 
{  
   // Get the new solution
   reconstruction (Residual);
   
   // Get the Residual
   if (UseMask == False)
    {
		for (int i=0; i<Nx; i++) Residual(i) = Data(i) - Residual(i);
		}
  else
	  {
		for (int i=0; i<Nx; i++) Residual(i) = (Data(i) - Residual(i))*MaskedData(i);
    }
}

/****************************************************************************/

void MCA1D::reset() 
{ // Reset all internal variables 

    Bord = I_MIRROR;
    NbrBase = 0;
    LastSoftThreshold =  DEF_MCA1D_LAST_SOFT_THRESHOLD;
    Nbr_Iter = DEF_MCA1D_NBR_ITER;
    PositivRecSig = False;
    DataSigmaNoise = 0.;
    Verbose = False;
    RemoveSmoothPlane = False;
    Write = True;
     
    AT_NbrScale1D = 4;
    AT_PositivRecSig = PositivRecSig;
    AT_KillLastScale = False;
    AT_Trans = NULL;    
    AT_OnlyPositivDetect = False;
    AT_SuppressIsolatedPixel = False;
    
    WT_NbrScale1D = 4;
    WT_SelectFilter = NULL;
    WT_SB1D = NULL;
    WT = NULL;
    WT_PositivRecSig = PositivRecSig;
    WT_OnlyPositivDetect  = False;
    WT_SuppressIsolatedPixel = False;
    WT_NonOrthoFilterBank= False;
    
    UseNormL1= False;
    LambdaTV = 0;
    
    COS_PositivRecSig = False;
    COS_BlockSize = 256;
    COS_Overlapping=False;
    COS_Sensibility = 1.;
    COSMin = 0.;
    CosWeightFirst = False;
    COS_SuppressDctCompCont = False;
      
    MCOS_PositivRecSig = False;
    MCOS_FirstBlockSize = 8;
    MCOS_Overlapping = False;
    MCOS_FirstDetectScale = 0;    
    
    TabALDCT = NULL;
    ALDCT_InfoCost=0;
    
    NSigmaSoft = 1.5;
    TotalVariation = False;
    Linear = True;
    
    UseMadForLambdaEstimation=False;
    CFAR=CFDR=False;
    CFAR_NsigmaMad=5.;
    CFDR_Qparam=0.25;
}
 
/****************************************************************************/

void MCA1D::alloc(int NSig) // , fltarray * Tab) 
{
   Nx = NSig;
   // int s;
   // for (int i=0; i < NbrBase; i++) Tab[i].alloc(Nx, "alloc");

   cout << "ALLO ALLO Signal size: Nx = " << Nx << endl;
   Resi.alloc(Nx,"resi");
   for (int i =0; i < NbrBase; i++) 
   { 
       TabSigRec[i].alloc(Nx,StringMCA1DBase(TabSelect[i]));
       cout << "ALLOC: Selection:" << StringMCA1DBase ((type_mca1dbase) TabSelect[i])  << endl;
       
       switch (TabSelect[i]) {
          
             case MCA1D_ATROU:
 
                 AWT.alloc (AT_Trans, Nx, AT_NbrScale1D);
                 AWT.Verbose = Verbose;
                 AWT.Write = Write;
                 AWT.Bord = I_MIRROR;
                 AWT.OnlyPositivDetect = AT_OnlyPositivDetect;
                 AWT.SuppressIsolatedPixel = AT_SuppressIsolatedPixel;
                 AWT.SigmaNoise = DataSigmaNoise;
                 AWT.UseNormL1 = UseNormL1;
                 AWT.NSigmaSoft = NSigmaSoft;                 
                 break;
             case MCA1D_WT:
                 cout << "ALLOC WT " << endl;
                 if (WT_NbrScale1D < 2) {
                     cout << "Error: bad number of scales ... " << endl;
                     exit(-1);
                 }
		 
		 if (WT_NonOrthoFilterBank == True)
		 {
 		    WT_USB1D = new  UndecSubBandFilter(U_B3SPLINE_2);
                    WT = new  PAVE_1D_WT(*WT_USB1D);
 		 }
		 else
		 {
                    WT_SelectFilter = new FilterAnaSynt(F_LEMARIE_3);
                    WT_SB1D = new SubBandFilter(*WT_SelectFilter, NORM_L2);
                    WT = new  PAVE_1D_WT(*WT_SB1D);
		 }
		 WT_Trans.alloc(Nx, WT_NbrScale1D);
                 if (Verbose == True)
                 {
                    cout << "WT: Nbr scales = " << WT_NbrScale1D << endl;
                 }   
                 cout << "ALLOC WT " << endl;
                 break;
              case MCA1D_COS:
		 LDCT.alloc(Nx, COS_BlockSize, COS_Overlapping, CosWeightFirst);
                 LDCT.Verbose = Verbose;
                 LDCT.Write = Write;
                 LDCT.SigmaNoise = DataSigmaNoise;
                 LDCT.UseNormL1 = UseNormL1;
                 LDCT.NSigmaSoft = NSigmaSoft;
                 LDCT.SigmaNoise  =  DataSigmaNoise;
                 LDCT.COSSupportMin = COSMin;
                 LDCT.SuppressDctCompCont = COS_SuppressDctCompCont;  
                 LDCT.COS_Sensibility =  COS_Sensibility;            
                 break;
                 
	     case MCA1D_MCOS:
		 M_DCT.Verbose=Verbose;
	         M_DCT.BlockOverlap = MCOS_Overlapping;
	         M_DCT.alloc(Nx, MCOS_NbrScale1D, MCOS_FirstBlockSize);
                 M_DCT.SigmaNoise = DataSigmaNoise;
                 M_DCT.COSMin = COSMin;
                 M_DCT.COS_Sensibility =  COS_Sensibility;            
                 M_DCT.Write = Write;
                 M_DCT.Verbose = Verbose;
                 M_DCT.UseNormL1 = UseNormL1;
                 M_DCT.NSigmaSoft = NSigmaSoft;                 
                 M_DCT.OnlyPositivDetect = FALSE;
                 break;
                 
             case MCA1D_ALDCT:
                 AL_DCT.alloc(TabALDCT, Nx, ALDCT_NbrScale1D);
                 AL_DCT.SigmaNoise = DataSigmaNoise;
                 AL_DCT.Sensibility = ALDCT_Sensibility;
                 AL_DCT.Write = Write;
                 AL_DCT.Verbose = Verbose;
                 AL_DCT.UseNormL1 = UseNormL1;
                 AL_DCT.NSigmaSoft = NSigmaSoft;
                 break;
                 
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
        }
   }

   if (NbrBase  < 1) {
       cout << "Error: no selected transform ... " << endl;
       exit(-1); 
   }
   
   //if (Verbose == True)
   //         cout << "Number of selected bases = " << NbrBase << endl;
}


/****************************************************************************/

void MCA1D::free() 
{
   for (int b=0; b<NbrBase; b++) 
   {
      switch(TabSelect[b]) 
      {
          case MCA1D_ATROU: AWT.free(AT_Trans, AT_NbrScale1D); break;
           case MCA1D_WT:
                 if  (WT_NbrScale1D > 2) {
                    delete WT_SelectFilter;
                    delete WT_SB1D;
                    delete WT;
                    WT_NbrScale1D = 0;
                    WT_SelectFilter = NULL;
                    WT_SB1D = NULL;
                     WT = NULL;
                 }
               break;
            case MCA1D_COS:
	       break;
            case MCA1D_MCOS:
	         M_DCT.free();
	         break;
	       break;
            case MCA1D_ALDCT:
	       break;
	    default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
      }    
   }
   NbrBase = 0;
}   
   
/****************************************************************************/

void MCA1D::remove_smooth_plane (fltarray& Signal, int NbSmoothRemove) 
{
   Smooth.alloc(Nx);

   for (int i=0;i<NbSmoothRemove;i++) 
   {
      AWT.Bord = I_MIRROR;
      AWT.transform (Signal, AT_Trans);
      Smooth += AT_Trans[AT_NbrScale1D-1];
      Signal -= AT_Trans[AT_NbrScale1D-1];
   }

   if (Write) io_1d_write_data ("Smooth", Smooth);
}

/****************************************************************************/

float MCA1D::compute_first_lambda (fltarray& Signal) 
{
   vector<float> VectMax;
   float M;
   for (int b=0; b < NbrBase; b++) 
   {
      switch(TabSelect[b]) 
      {
           case MCA1D_ATROU:
                AWT.transform (Signal, AT_Trans);
		M = AWT.getAbsMaxTransf (AT_Trans);
		if (Verbose == True) 
		  cout << "Atrou: first threshold = " << M << endl;
                VectMax.push_back(M);
                break;
            case MCA1D_WT:
                WT->transform (Signal, WT_Trans, WT_NbrScale1D);
		M = WT->getAbsMaxTransf (WT_Trans, DataSigmaNoise, WT_OnlyPositivDetect);
		if (Verbose == True) 
                   cout << "Undecimated WT: first threshold = " << M << endl;
                VectMax.push_back(M);
	        break;
 	    case MCA1D_COS:
                LDCT.transform(Signal);
		M = LDCT.getAbsMaxTransf (Signal);
		if (Verbose == True) 
                   cout << "Local DCT: first threshold = " << M << endl;
 		VectMax.push_back(M);
                break;
 	     case MCA1D_MCOS:      
                M_DCT.transform(Signal);
		M = M_DCT.getAbsMaxTransf ();
		if (Verbose == True) 
                   cout << "MultiScale DCT: first threshold = " << M << endl;
		VectMax.push_back(M);
                break;
           case MCA1D_ALDCT:
                AL_DCT.transform (Signal, TabALDCT);
		M = AL_DCT.getAbsMaxTransf (TabALDCT);
		if (Verbose == True) 
                   cout << "Adapted Local DCT: first threshold = " << M << endl;
                VectMax.push_back(M);
                break;
	   default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
      }
   }
   if (NbrBase == 1) return VectMax[0];
   else 
   {
      sort (VectMax.begin(), VectMax.end());
      M = 1.01*VectMax[NbrBase-2];
      if (Verbose == True) cout << "Fist Threshold = " << M << endl;
      return M;
   }
   return 1.;    
}
  
/****************************************************************************/

void MCA1D::decomposition (fltarray& Signal) 
{ 
   int i,b;
   // RegulIma RIM;
   MR_Regul RIM;
   RIM.NbrScale = 5;
   RIM.ExpDecreasingLambda = False;
   
   if ((CFAR == True) || (CFDR == True)) UseMadForLambdaEstimation = True;
   
   Lambda = FirstSoftThreshold;
   if (Nbr_Iter == 1) Lambda = LastSoftThreshold;
   StepL = (FirstSoftThreshold - LastSoftThreshold ) / (float) (Nbr_Iter-1);

   float DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
   // if (Linear == False) StepL = DeltaThreshold / (float) (Nbr_Iter-1);

   // int NbrScaleTV = 5;
   // RIM.GradOperType = OPER_ENERGY;
   // RIM.GradOperType = OPER_LAPLACIAN;
   fltarray TVSig;
   
   if (Signal.nx() != Nx)  
   {  
       cout << "Error: MCA class not initialized with this signal size" 
            << Data.nx() << endl;
       cout << "       Expected size is " << Nx << endl;
       exit(-1);
   }
   Resi = Signal;
   if ((TotalVariation == True) && (LambdaTV > 0)) TVSig.alloc(Nx, "TV");
       
    // if (LambdaTV > 0)   LambdaTV *= sqrt((double) 2.) / Nbr_Iter;
  
   for (int Iter=0; Iter < Nbr_Iter; Iter++) 
   {
       if (Iter>0) 
       {
          if (Linear == True) Lambda -= StepL; 
          else Lambda =  LastSoftThreshold 
                          + DeltaThreshold  *(1.-erf(2.8*Iter/ Nbr_Iter));
          if (Lambda < LastSoftThreshold) Lambda = LastSoftThreshold;
       }
        
       if (Verbose == True)
           cout << endl << "Iteration " << Iter+1 << " Lambda = " << Lambda << endl;

       
       for (b=0; b < NbrBase; b++) 
       {
           cout << " Base " << b+1 << ": Before" << endl;
           TabSigRec[b].info("REC");
           Resi.info("RESI");
           
          TVSig = TabSigRec[b];
          if (Write)
	      {
	         char Name[256];
             sprintf(Name, "TVSig_%d_%d", b, Iter);
             if (Iter % 1 == 0) io_1d_write_data (Name,TVSig );
	      }
	      TabSigRec[b] += Resi;
          if (Write == False) TabSigRec[b].info();

          switch(TabSelect[b]) 
	      {
 	          case MCA1D_ATROU:
                if (Verbose == True) cout << "Atrou proj " << endl;
                atrou_proj(TabSigRec[b], Lambda, Iter);
                break;
              case MCA1D_WT:
                if (Verbose == True) cout << "WT proj " << endl;
                wt_proj(TabSigRec[b],  Lambda, Iter);                    
                break;
 	          case MCA1D_COS:
                if (Verbose == True) cout << "Local DCT proj " << endl;
                cos_proj (TabSigRec[b], Lambda, Iter);   
                break;
 	          case MCA1D_MCOS:
                if (Verbose == True) cout << "Local MDCT proj " << endl;
                mcos_proj(TabSigRec[b], Lambda, Iter);                    
                break;
 	          case MCA1D_ALDCT:
                if (Verbose == True) cout << "ALDCT proj " << endl;
		        aldct_proj(TabSigRec[b], Lambda, Iter);
		        break;
 	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }
          
// 	  if (UseNormL1 == True) {
	  
	     for (i=0; i<Nx; i++)
	          Resi(i) = (TabSigRec[b])(i) - TVSig(i); 
          
//	  }
// 	  else {
// 	     for (i = 0; i < Nx; i++)
// 	       {
// 	          Resi(i) = TabSigRec[b](i);
// 	          TabSigRec[b](i) += TVSig(i); 
// 	       }
//  	     }
	  
         // New residual calculation
         if (    (TotalVariation == True) 
             && (LambdaTV > 0)
             && (TabSelect[b] != MCA1D_COS))
         {
//                RIM.obj_regul(TVIma, TabSigRec[b], LambdaTV);
//	          regul_tv_haarwt(TabSigRec[b], TabSigRec[b], LambdaTV, NbrScaleTV);
// 	          RIM.obj_regul(TVIma, TabSigRec[b], LambdaTV*Noise_Ima);
//	          RIM.im_soft_threshold (TabSigRec[b], TabSigRec[b], 
//                                      LambdaTV*Noise_Sig);
	    }
	  
	 // Positivity constraint
	 switch(TabSelect[b]) 
	 {
             case MCA1D_ATROU: 
	        if (AT_PositivRecSig == True) (TabSigRec[b]).inf_threshold(0.0);
                break;
             case MCA1D_WT:
                if (WT_PositivRecSig == True) (TabSigRec[b]).inf_threshold(0.0);
                break;
	     case MCA1D_MCOS:
	        break;
	     case MCA1D_COS: 
 	        break;
	     case MCA1D_ALDCT:
	        break;                 
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }     
	  make_residual (Signal, Resi);
         if (Write == True)
	 {
	    char Name[256];
	    sprintf(Name, "Resi_%d_%d", b, Iter);
            if (Iter % 1 == 0) io_1d_write_data (Name,Resi );
	 }
	  
	  // if (Verbose == True)
	  {
             cout << " Base " << b+1 << ": Reconstruction" << endl;
             TabSigRec[b].info();
             Resi.info();
	  }
	  
      } // ENDFOR (b=0 ...)
      if (Write && Iter % 1 == 0) write_allima("xx", Iter);
   } // ENDFOR (Iter=...)
   
} 




/****************************************************************************/
/*         PROJECTION                                                       */
/****************************************************************************/

void MCA1D::atrou_proj (fltarray& SigAtrou, float Lambda_Sigma,  int IterNumber) 
{
    fltarray TabLevel(AT_NbrScale1D);
    fltarray TabNoise(AT_NbrScale1D);
    
    if (UseMadForLambdaEstimation == True)
    {  
         AWT.transform (Resi, AT_Trans);
         for (int s=AT_FirstDetectScale; s<AT_NbrScale1D-1; s++)
	 { 
            float *Ptr = (AT_Trans[s]).buffer();
            int Npix = Nx;
            float SigmaMad = get_sigma_mad(Ptr, Npix);
	    if (SigmaMad < DataSigmaNoise*AWT.norm_band(s)) SigmaMad =  DataSigmaNoise*AWT.norm_band(s);
 	    // float SigmaMad = (AT_Trans[s]).sigma();
            TabNoise(s) = SigmaMad;
            if (CFAR == True) 
            {
	      TabLevel(s) = CFAR_NsigmaMad*SigmaMad;
	      if (Verbose == True) cout << "   AT: Sigma Mad = " << SigmaMad << " NSig = " << CFAR_NsigmaMad << " Lambda = " << CFAR_NsigmaMad*SigmaMad << endl;
            }
            else
            {
              double  NSig,PDet=fdr_gauss_threshold(Ptr, Npix, CFDR_Qparam, SigmaMad);
	      if (PDet < 1.e-07) NSig = 5.2;
	      else
	      {
	         NSig = ABS(xerfc(0.5+(1-PDet)/2.));
	         if ((NSig > 5.2) ||(NSig < 0)) NSig = 5.2;
 	      }
              if (NSig < LastSoftThreshold) NSig = LastSoftThreshold;
	      TabLevel(s) = NSig*SigmaMad;
 	      if (Verbose == True) cout << "   AT: Sigma Mad = " << SigmaMad << " NSigFDR = " << NSig << " Lambda = " << NSig*SigmaMad << endl;
	    } 
	 }        
    }
      
   AWT.transform (SigAtrou, AT_Trans);
   AWT.KillScaleNotUsed (AT_Trans, AT_FirstDetectScale);
   // AWT.Threshold (AT_Trans, Lamba_Sigma, IterNumber);
   
   int NumberDetectedCoef=0;
   for (int s=AT_FirstDetectScale; s<AT_NbrScale1D-1; s++) 
   {
      float NSig = (s == 0) ? Lambda_Sigma + 1: Lambda_Sigma;
      if (Lambda_Sigma == 0) NSig = 0;
      float Norm = (s == AT_NbrScale1D-1) ? AWT.norm_band(s-1): AWT.norm_band(s);
      float Noise =  DataSigmaNoise*Norm;
      float Level = NSig*Noise;
      if (UseMadForLambdaEstimation == True) 
      {
         Level = TabLevel(s);
	 Noise = TabNoise(s);
      }
      for (int i=0; i< Nx; i++) 
      {
          float Coef = (AT_Trans[s])(i);  
	  (AT_Trans[s])(i) = mca1d_update(Coef, Level, Noise);
          if (AT_OnlyPositivDetect == True)  if ((AT_Trans[s])(i) < 0) (AT_Trans[s])(i) = 0;
          if (fabs((AT_Trans[s])(i)) > 0) NumberDetectedCoef++;
      }
       
      // Support 
//       for (int i=1; i< Nx-1; i++) 
//       {
//             float Pixm =  (AT_Trans[s])(i-1);
//             float Pixp =  (AT_Trans[s])(i+1);
// 	    if ((Pixm == 0) && (Pixp == 0))  (AT_Trans[s])(i) = 0.;
//             if (fabs((AT_Trans[s])(i)) > 0) NumberDetectedCoef++;
//       }
   }
   if (Verbose)  cout << "   AT: Number detected coef : " << NumberDetectedCoef << endl;
      
   if (AT_KillLastScale == True)  AWT.KillLastScale (AT_Trans);
   AWT.recons(AT_Trans, SigAtrou);
}

/****************************************************************************/

void MCA1D::cos_proj (fltarray& SigRec, float Lamba_Sigma,  int IterNumber) 
{
   float NSig=Lamba_Sigma;
   
   LDCT.transform (Resi);
   io_1d_write_data("xx_resi_dct.fits", LDCT.trans_Sig());  
   if (UseMadForLambdaEstimation == True)
   {
      float *Ptr = LDCT.trans_Sig().buffer();
      int NpixDCT = LDCT.nx();
      float SigmaMad = get_sigma_mad(Ptr, NpixDCT);
      if (SigmaMad < DataSigmaNoise) SigmaMad = DataSigmaNoise;
      // float SigmaMad =  LDCT.trans_Sig().sigma();
      LDCT.SigmaNoise = SigmaMad;
      if (CFAR == True) 
      {
         NSig = CFAR_NsigmaMad;
	 if (Verbose == True) cout << "   DCT: Sigma Mad = " << SigmaMad << " NSig = " << CFAR_NsigmaMad << " Lambda = " << CFAR_NsigmaMad*SigmaMad << endl;
      }
      else
      {
         float PDet=fdr_gauss_threshold(Ptr, NpixDCT, CFDR_Qparam, SigmaMad);
	 if (PDet < 1.e-07) NSig = 5.2;
         else  NSig = ABS(xerfc(0.5+(1-PDet)/2.));
	 if ((NSig > 5.2) ||(NSig < 0)) NSig = 5.2;
 	 if (Verbose == True) cout << "   DCT: Sigma Mad = " << SigmaMad << " NSigFDR = " << NSig << " Lambda = " << NSig*SigmaMad << endl;
      }
   }
   
   LDCT.transform (SigRec); 
   LDCT.threshold (NSig);
   LDCT.recons (SigRec);
}

/****************************************************************************/
           
void MCA1D::wt_proj(fltarray& Sig, float Lamba_Sigma, int IterNumber) 
{
   float NSigma = Lamba_Sigma;
   fltarray WTab1DSignifLevel(MAX_SCALE_1D);
   WTab1DSignifLevel(0) = 0.723490;
   WTab1DSignifLevel(1) = 0.285450;
   WTab1DSignifLevel(2) = 0.177948 ;
   WTab1DSignifLevel(3) = 0.122223;
   WTab1DSignifLevel(4) = 0.0858113;
   WTab1DSignifLevel(5) = 0.0605703;
   WTab1DSignifLevel(6) = 0.0428107;
   WTab1DSignifLevel(7) = 0.023;
   WTab1DSignifLevel(8) = 0.017;
	     
  for (int i = 9; i < MAX_SCALE_1D; i++) WTab1DSignifLevel(i) = WTab1DSignifLevel(i-1) / 2;
	        
//cout << "TRans  " << endl;
   WT->transform (Sig, WT_Trans, WT_NbrScale1D);
//   cout << "TRans OK" << endl;
   
   // WT->threshold (WT_Trans, Lamba_Sigma, IterNumber);
    for (int s=0; s < WT_Trans.ny()-1; s++)
   {      
      float NSig = (s < 1) ? NSigma + 1: NSigma;
      float Noise = DataSigmaNoise;
      float Level = NSig*Noise;
      if (WT_NonOrthoFilterBank == True) Level *= WTab1DSignifLevel(s); 
      for (int i=0; i <  WT_Trans.nx(); i++)
      {
         float Coef = WT_Trans(i,s); 
         WT_Trans(i,s) = mca1d_update(Coef, Level, Noise);
	 if (s < WT_FirstDetectScale) WT_Trans(i,s) = 0;
         else if (WT_OnlyPositivDetect == True)  if (WT_Trans(i,s) < 0)  WT_Trans(i,s) = 0;
      }
    }
   
  // cout << "TRES OK" << endl;
   
   if (AT_KillLastScale == True)  
         for (int i=0; i <  WT_Trans.nx(); i++) WT_Trans(i,WT_Trans.ny()-1) = 0;
   
   WT->recons(WT_Trans, Sig, WT_NbrScale1D); 
  // cout << "REC OK" << endl;  
}

/****************************************************************************/

void MCA1D::mcos_proj(fltarray& SigRec, float Lamba_Sigma, int IterNumber) 
{
   M_DCT.transform(SigRec);
   M_DCT.KillScaleNotUsed (MCOS_FirstDetectScale);
   M_DCT.threshold (Lamba_Sigma, IterNumber);
   M_DCT.recons(SigRec);
}

/****************************************************************************/

void MCA1D::aldct_proj(fltarray& SigRec, float Lamba_Sigma, int IterNumber)
{
   AL_DCT.transform (SigRec, TabALDCT);
   AL_DCT.threshold (TabALDCT, Lamba_Sigma, IterNumber);
   fltarray BestBasis;
   AL_DCT.getBestBasis (TabALDCT, BestBasis, ALDCT_InfoCost);
   AL_DCT.recons(SigRec);
}

/****************************************************************************/

void MCA1D::write_allima(char *FileName, fitsstruct & Header) 
{
   char Pre[MAXCHAR];
   char FName[MAXCHAR];
   
   io_strcpy_prefix(Pre, FileName);
   for (int b=0; b < NbrBase; b++)
   {
      switch(TabSelect[b])
      {
          case MCA1D_ATROU:    
	       sprintf(FName,"%s_atrou", Pre);
	       break;
          case MCA1D_WT:
	       sprintf(FName,"%s_wt", Pre);
	       break;
	  case MCA1D_COS:
	       sprintf(FName,"%s_cos", Pre);
	       break;
	  case MCA1D_MCOS:
	       sprintf(FName,"%s_mcos", Pre);
	       break;
	  case MCA1D_ALDCT:
	       sprintf(FName,"%s_aldct", Pre);
	       break;
	  default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
       }
       io_1d_write_data(FName, TabSigRec[b], &Header);  
       // fits_write_fltarr (FName, TabSigRec[b], &Header);    
   }
   sprintf(FName,"%s_resi", Pre);
   io_1d_write_data(FName, Resi, &Header);
}

/****************************************************************************/

void MCA1D::write_allima(char *FileName, int It) 
{
   char Pre[MAXCHAR];
   char FName[MAXCHAR];
   io_strcpy_prefix(Pre, FileName);
   for (int b=0; b < NbrBase; b++)
   {
      switch(TabSelect[b])
      {
 	  case MCA1D_ATROU:    
	       sprintf(FName,"%s_atrou", Pre);
	       break;
          case MCA1D_WT:
	       sprintf(FName,"%s_wt", Pre);
	       break;
	  case MCA1D_COS:
	       sprintf(FName,"%s_cos", Pre);
	       break;
          case MCA1D_MCOS:
	       sprintf(FName,"%s_mcos", Pre);
	       break;
          case MCA1D_ALDCT:
	       sprintf(FName,"%s_aldct", Pre);
	       break;
          default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
       }
       sprintf(FName, "%s_%d", (char *) FName, It);
       io_1d_write_data (FName, TabSigRec[b]);    
   }
   sprintf(FName,"%s_resi", Pre);
   sprintf(FName, "%s_%d", FName, It);
   io_1d_write_data (FName, Resi);
}

/****************************************************************************/

void MCA1D::infoVerbose (char *nameSignalIn, char *nameSignalOut) 
{
  if (Verbose == True) 
    { 
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << nameSignalIn << endl;
    cout << "File Name Out = " << nameSignalOut << endl; 
    cout << "Nbr iterations = " << Nbr_Iter << endl;
    cout << "First threshold = " << FirstSoftThreshold << endl;
    cout << "Last threshold = " << LastSoftThreshold << endl;
    if (RemoveSmoothPlane)
      cout << "Remove smooth plane" << endl;
    if (Linear) cout << "Linear step descent" << endl;
    else cout << "Non linear step descent" << endl;
    cout << "Estimated noise in signal : " << DataSigmaNoise << endl;
    if (UseNormL1) 
      cout << "Used L1 norm in threshold" << endl;
      cout << "Selected basis : " << endl;
      for (int i =0; i < NbrBase; i++) 
        cout << " ==> " 
        << StringMCA1DBase ((type_mca1dbase) TabSelect[i]) 
        << endl; 
            
      for (int i =0; i < NbrBase; i++) 
        {
        if (TabSelect[i] == MCA1D_ATROU) 
          {
          cout << "A Trous transform" << endl;
          cout << "   Nbr scale = " << AT_NbrScale1D << endl;
          if (AT_KillLastScale)
            cout << "   Kill last scale" << endl;
          if (AT_PositivRecSig)
            cout << "   Positiv constraint on signal" << endl;
            cout << "   Scale of first detection : " << AT_FirstDetectScale << endl;
          if (AT_SuppressIsolatedPixel) 
            cout << "   Suppres isolated pixel on support" << endl;
          }
        if (TabSelect[i] == MCA1D_WT) 
          {   
          cout << "HDWT transform" << endl;
          cout << "   nb scale : " << WT_NbrScale1D << endl;
          if (AT_PositivRecSig)
            cout << "   Positivity constraint on signal" << endl;
            cout << "   Scale of first detection : " << AT_FirstDetectScale
            << endl;
          if (WT_SuppressIsolatedPixel) 
            cout << "   Suppres isolated pixel on support" << endl;
          }
        if (TabSelect[i] == MCA1D_COS) 
          {             
          cout << "DCT transform" << endl;
          cout << "   COS_BlockSize : " << COS_BlockSize << endl;
          if (COS_Overlapping)
            cout << "   COS Overlapping" << endl;
            cout << "   Cos sensibility :  " << COS_Sensibility << endl;
          if (COS_SuppressDctCompCont)
            cout << "   Suppres continue comp (first coeff)" << endl;
            cout << "   Cos support min : " << COSMin << endl;
          if (COS_PositivRecSig)
            cout << "   Positivity constraint on signal" << endl;  
         }
       if (TabSelect[i] == MCA1D_MCOS) 
         {  
         cout << "MDCT transform" << endl;
         cout << "   nb scale : " << MCOS_NbrScale1D  << endl;
         }
       if (TabSelect[i] == MCA1D_ALDCT) 
         {  
         cout << "ALDCT transform" << endl;
         cout << "   ALDCT sensibility : " << ALDCT_Sensibility << endl;
         }
       }
    }   
}     
    

/***********************************************************************/
