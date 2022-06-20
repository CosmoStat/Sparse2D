/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.2
**
**    Author: 09/20/99
**
**    Date:  01/03/02
**    
**    File:  wth_MC2D.cc
**
*******************************************************************************/

#include "wth_MC2D.h"
#include "MR1D_Obj.h"
#include "SB_Filter.h"

// for test
int xPosTrace = 20;
int yPosTrace = 31;
int bNumPlan  = 0;


void WTH_MC2D::filt_Alloc (int Length) {

   i_Length = Length;
   o_wthMc2d.alloc(i_Length);
}


WTH_MC2D::~WTH_MC2D () {}


void WTH_MC2D::filt_Transform (MultiResol  *TabMRin,  MultiResol *TabMRout,
                               MRNoiseModel *TabNoiseIn, MRNoiseModel *TabNoiseOut) { 
	 		  
   if (TabNoiseIn != (MRNoiseModel*)NULL) {
 
      filt_TransfSignal (TabMRin, TabMRout, NULL);
      filt_ComputeNoise (TabMRin, TabNoiseIn, TabNoiseOut); 
      filt_Threshold (TabMRout, TabNoiseOut);
   }
 
   else filt_TransfSignal (TabMRin, TabMRout);
}

 
void WTH_MC2D::filt_TransfSignal (MultiResol  *TabMRin,  MultiResol *TabMRout,
                                  MRNoiseModel *TabNoiseIn) {
				  
   int Nbr_Plan = iround(log((float)i_Length/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   MR_1D Mr1d_Data;
   Mr1d_Data.alloc (i_Length, TO1_MALLAT, Nbr_Plan, &o_FAS, NORM_L2, False);
   //Mr1d_Data.SB_Filter = F_HAAR;
   //Mr1d_Data.Border =  I_CONT;
   
   for (int b=0; b<TabMRin[0].nbr_band()-1; b++)
      for (int i=0; i<TabMRin[0].size_band_nl(b); i++)
         for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
	 
	    int k=0;
	    for (k=0; k<i_Length; k++)
               o_wthMc2d(k) = (TabMRin[k])(b,i,j);
	       
	    // if (b==bNumPlan && i==xPosTrace && j==yPosTrace) fits_write_fltarr("testIN.fits", o_wthMc2d);
	    Mr1d_Data.transform (o_wthMc2d);
	        
	    k=0;
	    for (int s=0;s<Mr1d_Data.nbr_band();s++)
	       for (int p=0;p<Mr1d_Data.size_scale_np(s);p++) {
                  (TabMRout[k])(b,i,j) = Mr1d_Data(s,p);
		  k++;
	       }
	       
	    if (b==bNumPlan && i==xPosTrace && j==yPosTrace) {
	       for (k=0; k<i_Length; k++) o_wthMc2d(k) = (TabMRin[k])(b,i,j);
	       // fits_write_fltarr("testOUT.fits", o_wthMc2d);
	    }
	    
	 } 
} 


void WTH_MC2D::filt_InvTransform (MultiResol *TabMRin , MultiResol *TabMRout, Bool KillLastScale) {
   
   int Nbr_Plan = iround(log((float)i_Length/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   MR_1D Mr1d_Data;
   Mr1d_Data.alloc (i_Length, TO1_MALLAT, Nbr_Plan, &o_FAS, NORM_L2, False);   
   //MR_1D Mr1d_Data (i_Length, WP1_MALLAT, "MC1D_Transform", Nbr_Plan, False);
   //Mr1d_Data.SB_Filter = F_HAAR;   
   //Mr1d_Data.Border =  I_CONT;
   
   for (int b=0; b<TabMRin[0].nbr_band()-1; b++)
      for (int i=0; i<TabMRin[0].size_band_nl(b); i++)
         for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
	 
	    int k=0;
	    for (int s=0;s<Mr1d_Data.nbr_band();s++)
	       for (int p=0;p<Mr1d_Data.size_scale_np(s);p++) {
                  Mr1d_Data(s,p) = (TabMRin[k])(b,i,j);
		  k++;
		}
		
	    //if (b==bNumPlan && i==xPosTrace && j==yPosTrace) fits_write_fltarr("testIN.fits", o_wthMc2d);
	    if (KillLastScale == True) Mr1d_Data.rec_adjoint (o_wthMc2d, False);
	    else Mr1d_Data.rec_adjoint (o_wthMc2d);
	    // if (b==bNumPlan && i==xPosTrace && j==yPosTrace) fits_write_fltarr("testREC.fits", o_wthMc2d);
	      
	    for (k=0; k<i_Length; k++) {
 //if (o_wthMc2d(k) < 0) o_wthMc2d(k) = 0.0;
               (TabMRout[k])(b,i,j) = o_wthMc2d(k);
	    }	    
	 } 
}


void mallat_1d_noise_transform (fltarray &Data, fltarray &Mallat, 
                                int Nbr_Plan, SubBand1D &SBF) {
   int i,s;
   int Np_2,Nps;
   float *ImagLow, *ImagHigh, *DataResol;
   int Np = Data.nx();
      
   // cout << "mallat : " << Nbr_Plan << " np = " << Np << endl;
    
   if (Np !=  Mallat.n_elem()) Mallat.alloc(Np);
    
   /* Allocation  */
   Nps = size_resol(1, Np);
   DataResol    = new float[Np];
   ImagHigh     = new float[Nps];
   ImagLow      = new float[Nps];
    
   for (i=0; i< Np; i++) DataResol[i] = Data(i);
    
   /* Compute the wavelet coefficients */
   for (s = 0; s < Nbr_Plan-1; s++) {  
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
    	
        SBF.noise_transform (Nps,DataResol,ImagLow,ImagHigh);
        for (i=0; i < Nps/2; i++)  Mallat(Np_2+i) = ImagHigh[i];
    
    	if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++)  
    	               DataResol[i] = ImagLow[i];
   }
   Np_2 = size_resol(Nbr_Plan-1, Np);
   for (i=0; i< Np_2; i++) Mallat(i) = ImagLow[i];
   delete [] ImagHigh;
   delete [] ImagLow;
   delete [] DataResol;
}

void WTH_MC2D::filt_ComputeNoise (MultiResol *TabMRin, MRNoiseModel *TabNoiseIn,
                                  MRNoiseModel *TabNoiseOut) {
   
   int Nbr_Plan = iround(log((float)i_Length/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   int k;
   
   type_noise noise = TabNoiseIn[0].which_noise();

   for (int b=0; b<TabMRin[0].nbr_band()-1; b++) {
   
      if (!one_level_per_pos_2d(noise)) {
       	       
         //cout << "sigma IN (scale:" << b << ") :"; 
 	 //for (k=0; k<i_Length; k++) cout << TabNoiseIn[k].sigma(b,0,0)  << ",";
 	 //cout << endl;
	 
 	 for (k=0; k<i_Length; k++)
             o_wthMc2d(k) = pow ((double) TabNoiseIn[k].sigma(b,0,0), 2.); 
	          
 	 fltarray ao_wthNoise(i_Length);
 	 SubBandFilter SBF(F_HAAR, NORM_L2);
 	 SBF.Border   = I_CONT;
 	 mallat_1d_noise_transform (o_wthMc2d, ao_wthNoise, Nbr_Plan, SBF); 
	 
         for (k=0; k<i_Length; k++)
             TabNoiseOut[k].sigma(b,0,0) = sqrt (ao_wthNoise(k));	      
         
	 //cout << "sigma OUT (scale:" << b << ") :";
 	 //for (k=0; k<i_Length; k++) cout << TabNoiseIn[k].sigma(b,0,0)  << ",";
 	 //cout << endl;   
	    
      } else {
         cout << "not implemented yet" << endl;
	 exit (-1);
      }
   }
   
//       for (int i=0; i<TabMRin[0].size_band_nl(b); i++)
//          for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
// 	 
// 	    int k=0;
// 	    for (k=0; k<i_Length; k++)
//                o_wthMc2d(k) = pow (TabNoiseIn[k].sigma(b,i,j), 2.);
// 	       
// 	    if (i==0 && j==0) {
// 	       cout << "sigma IN (scale:" << b << ") :"; 
// 	       for (k=0; k<i_Length; k++) cout << TabNoiseIn[k].sigma(b,i,j) 
// 	                                       << ",";
// 	       cout << endl;
// 	    }
// 	    if (b==bNumPlan && i==xPosTrace && j==yPosTrace) fits_write_fltarr("testnoiseIN.fits", o_wthMc2d);
// 	    
// 	    fltarray ao_wthNoise(i_Length);
// 	    SubBandFilter SBF(F_HAAR, NORM_L1);
// 	    SBF.Border   = I_CONT;
// 	    mallat_1d_noise_transform (o_wthMc2d, ao_wthNoise,
// 	                               Nbr_Plan, SBF);
// 				       
// 				       	   
// 	    if (b==bNumPlan && i==xPosTrace && j==yPosTrace) fits_write_fltarr("testnoiseOUT.fits", ao_wthNoise);
// 	    
//             for (k=0; k<i_Length; k++)
//                TabNoiseOut[k].sigma(b,i,j) = sqrt (ao_wthNoise(k));
// 	       
// 	    if (i==0 && j==0) {
// 	       cout << "sigma: OUT (scale:" << b << ") :"; 
// 	       for (k=0; k<i_Length; k++) cout << TabNoiseOut[k].sigma(b,i,j) 
// 	                                       << ",";
// 	      cout << endl;
// 	    }
// 	 }
//    } 
}


void WTH_MC2D::filt_Threshold (MultiResol  *TabMRin_out,
                               MRNoiseModel *TabNoiseOut) {
			       
   for (int k=0; k<i_Length; k++) {
   
   
      TabNoiseOut[k].threshold(TabMRin_out[k]);
 
//       if (int b==bNumPlan && i==xPosTrace && j==yPosTrace) {
// 	 fltarray o_test (i_Length);
// 	 for (int k=0; k<i_Length; k++) o_test(k) = (TabMRin_out[k])(bNumPlan,xPosTrace,yPosTrace);
// 	 fits_write_fltarr ("afterthre.fits", o_test);
//       }     
   }
}			       
			       
/*void WTH_MC2D::filt_Compute (MultiResol *TabMR, MRNoiseModel *TabMRNoiseModel=NULL) {

void WTH_MC2D::filt_Print () {}*/













