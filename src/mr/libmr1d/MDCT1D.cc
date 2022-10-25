/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Hubert Druesne, Philippe Querre
**
**    Date:  24/07/03 
**    
**    File:  MDCT1D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale DCT
**    -----------  
**                 
******************************************************************************/

#include "MDCT1D.h"
  
/****************************************************************************/
//                  Multiscale DCT
/****************************************************************************/

void MDCT1D::alloc(int Nx, int Nbr_Plan, int FirstBlockSize){
    int s;
    NbrScale = Nbr_Plan;
    TabBlockSize.alloc(Nbr_Plan);
    TabBlockSize(0) =  FirstBlockSize;
    for (s=1; s < Nbr_Plan; s++) 
    {
        if (2*TabBlockSize(s-1) <= Nx) TabBlockSize(s) = 2*TabBlockSize(s-1);
	else TabBlockSize(s) = TabBlockSize(s-1);
    }
    TabDCT = new LOCAL_DCT1D [Nbr_Plan];
    // TabTransMDCT = new fltarray [Nbr_Plan];
    
    if (Write == True) cout << " Multiscale DCT allocation: Nbr_Plan =  " << Nbr_Plan << endl;
    for (s=0; s < Nbr_Plan; s++)
    {
        TabDCT[s].alloc(Nx, TabBlockSize(s), BlockOverlap);
	// TabTransMDCT [s] = & (TabDCT[s]._DCTSig);
	// if (Verbose == True) 
	   printf("   Band %d, BlockSize = %d, Nx = %d \n", 
	          s+1, TabBlockSize(s), nx(s));
    }
    AWT.alloc(TabTransWT, Nx, Nbr_Plan);
}

/****************************************************************************/

void MDCT1D::transform(fltarray &Signal)
{
   int s;
   AWT.transform(Signal, TabTransWT);
   for (s = 0; s < NbrScale; s++)
   {
       if (Write == True) cout << "  Local DCT, band " << s + 1 << endl;
       TabDCT[s].transform(TabTransWT[s]);
       if (Verbose == True)
       {
          TabDCT[s]._DCTSig.info("DCT band");
 	  if (Write == True)
             cout << "  Norm band  = " << norm_band(s) << endl;
       }
   }
}


/****************************************************************************/

void MDCT1D::KillScaleNotUsed (int FirstDetScale) 
{
   FirstDetectScale = FirstDetScale;
   if (FirstDetectScale != 0) {   
      for (int j=0; j < FirstDetectScale; j++)
         TabDCT[j]._DCTSig.init();

   }
}

void MDCT1D::KillLastScale () 
{
   TabDCT[NbrScale-1]._DCTSig.init();
}

/****************************************************************************/

void MDCT1D::threshold (float NSigma, int IterNumber) {

  for (int s=FirstDetectScale; s<nbr_scale()-1; s++) {
 
      int BS = TabDCT[s]._B1DTrans.block_size();
      fltarray BlockCosSig(BS, "blocksig");
      fltarray MeanBlockCosIma(BS, "blocksig");
      float COSSupportMin = COSMin;
      if ((COSSupportMin < 0) || (COSSupportMin >= 0.5)) COSSupportMin = 0;
      COSSupportMin *= BS;
      
      float Noise = SigmaNoise * norm_band(s);   
      float Level = NSigma*Noise*COS_Sensibility;
      
      if (COSSupportMin == 0) {
     
         for (int i=0; i<nx(s); i++){
            float Coef = (*this)(s,i);  
            (*this)(s,i) = update(Coef, Level, Noise);
         }
      } 
      else {
   
         for (int bx=0; bx<TabDCT[s]._B1DTrans.nbr_block_nx(); bx++) {
        
            TabDCT[s]._B1DTrans.get_block_sig ( bx, TabDCT[s]._DCTSig, 
                                                BlockCosSig);
         
            for (int i=0; i<BS; i++) {
	       int Indi = (i < BS/2) ? i+BS/2: i - BS/2;
               int Rad = Indi - BS/2;
               BlockCosSig(i) =  update (BlockCosSig(i), Level, Noise);
	     
               if (Rad < COSSupportMin) BlockCosSig(i) = 0; 
 	    }
	    TabDCT[s]._B1DTrans.put_block_sig ( bx, TabDCT[s]._DCTSig, 
                                                BlockCosSig);
         }
      }
      
      if (Write){
         char Name[256];
         sprintf(Name,"MDCT__Thres_sc%d__iter%d", s+1, IterNumber);
         fits_write_fltarr (Name, TabDCT[s]._DCTSig);
      }
   
   }
   if (Write){
      char Name[256];
      sprintf(Name,"MDCT__Thres_sc%d__iter%d", nbr_scale(), IterNumber);
      fits_write_fltarr (Name, band(nbr_scale()-1));
   }
}

/****************************************************************************/

void MDCT1D::recons(fltarray &Signal){
   int s;
   for (s = 0; s < NbrScale; s++){
        if (Write == True) cout << "  Local DCT recons, band " << s + 1 << endl;
        TabDCT[s].recons(TabTransWT[s]);
   }
   if (Write == True) cout << "Isotropic Wavelet recons " << endl;
   AWT.recons(TabTransWT, Signal);
} 
   
/****************************************************************************/

float MDCT1D::norm_band(int s){
   static double TN[10] = 
    {0.721618,  0.286209, 0.178281, 0.121010, 0.083104, 0.0588543, 0.0446171,
     0.0289142, 0.0177307, 0.};
   if ((s < 0) || (s >= 10)) return 0;
   else return TN[s];
}

 
/****************************************************************************/

void MDCT1D::threshold(float NoiseLevel, float NSigma){
   int s,i;
   
   for (s=0; s < nbr_scale(); s++){
      float Noise = NoiseLevel * norm_band(s);   
      if (Write == True)
         cout << "  Local DCT, band " << s + 1 << " Level = " << Noise << endl;

      for (i=0; i < nx(s); i++){
         float Level = NSigma*Noise;
	 if (ABS((TabDCT[s])(i)) < Level) (TabDCT[s])(i) = 0.;
      }
   }
}

/****************************************************************************/

inline float MDCT1D::update (float CoefSol, float Threshold, 
                             float SoftLevel) {
					 
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

/****************************************************************************/

float MDCT1D::getAbsMaxTransf () {

   float absMaxLevel=0;
   
   for (int s=0; s<nbr_scale()-1; s++) {

      if (OnlyPositivDetect) (band(s)).inf_threshold(0.0);
      float MaxLevel = fabs(band(s).maxfabs());
      float NormMaxLevel = MaxLevel/SigmaNoise/norm_band(s)/COS_Sensibility;
                
      if (Write)
         cout << "Lambda M_Dct_"<< s+1 << "  : " <<  NormMaxLevel
              << ", Max amplitude_"<< s+1 << " : " << MaxLevel << endl;
      
      if (absMaxLevel <= NormMaxLevel) absMaxLevel = NormMaxLevel;
   }
   if (Verbose) cout << "Lambda MDCT1D : " << absMaxLevel << endl;

   return absMaxLevel;   
}

