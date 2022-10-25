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
**    Date:  10/04/02 
**    
**    File:  WT_Mirror.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION: Mirror Wavelet Decomposition  
**    ----------- 
**
******************************************************************************/

#ifndef _WT_MIRROR_H_
#define _WT_MIRROR_H_

#include "SB_Filter.h"


class MIRROR_2D_WT {
        int NlIma, NcIma;
        Ifloat *TabTransWT; // Pointer to the transformation  
	// Ifloat *TabMirrorBand;      
        SubBand1D *Ptr_SB1D; // Pointer to the 1D subband decomposition
        SubBand2D *Ptr_SBT;  // Pointer to the 2D subband decomposition
        LineCol *Ptr_LC;     // Pointer to the 2D Line-Column decomposition
	HALF_DECIMATED_2D_WT *Ptr_HD; // Pointer to the 2D partially decimated WT
	Ortho_2D_WT *Ptr_OWT; // Pointer to the orthogonal WT
	int NbrTotBand;   // Total number of bands
	int NbrBand;      // number of wavelet bands
	int Nbr_Plan;      // Number of dyadic scales
	intarray TabNl;   // TabNl(0..NbrTotBand) size of the bands
	intarray TabNc;   // TabNl(0..NbrTotBand) size of the bands
	intarray TabDepi;  
	intarray TabDepj;
	intarray TabIndBand;
	intarray TabResolBand;
	          // TabDepi, TabDepj, TabIndBand allows us to knows the
		  // position of a coefficents of scale b and position i,j
		  // in the tranformation TabTransWT
		  //   ind = TabIndBand(b)
		  //   Depi = TabDepi(b)
		  //   Depj = TabDepj(b)
		  //   Coef = TabTransWT[Ind](i+Depi,j+Depj)
		  
	Bool *TabDec; // TabDec[0..NbrPlan] 
	              // TabDec[i] == True if the scale i is decimated
		      // in the wavelet decomposition
		      	     
        void  set_tabdec(int NumUndec, int NbrPlan);
	// initialize TabDec
	void  make_alloc (int Nl, int Nc, int NbrPlan);
        // initialize all arrays and variables

       public:
          MIRROR_2D_WT (SubBand1D &SB1D) 
	    {  NlIma = NcIma = 0;
	       Ptr_SB1D = &SB1D;
	      Ptr_LC = new LineCol(SB1D);
	      Ptr_LC->Mirror=True;
	      Ptr_HD = new HALF_DECIMATED_2D_WT(SB1D);
	      Ptr_OWT = new Ortho_2D_WT(SB1D);
	      Ptr_SBT = new SubBand2D(SB1D);
	      TabTransWT = NULL;
	      TabDec = NULL;
	    }
	  void alloc(int Nl, int Nc, int NbrPlan, Bool *TabDec);
	  // Allocate for the class for images of size NlxNc and for
	  // NbrPlan scales. TabDec indicates which scale must be decimated
	  
	  void alloc(int Nl, int Nc, int NbrPlan, int NbrUndecimatedScale=-1);
	  // Allocate for the class for images of size NlxNc and for
	  // NbrPlan scales. The first NbrUndecimatedScale wavelet scales are 
	  // not decimated. By default, all scale are not decimated.
	  
	  void free();
	  // deallocate the CLASS
	  
 	  void transform (Ifloat &Imag);
	  // Apply the mirror wavelet transform
	  
 	  void recons(Ifloat &Imag);
	  // Apply the mirror wavelet reconstruction
	  
          int nbr_band() {return NbrTotBand;}
	  // return the total number of bands in the mirror wavelet transform
	  
	  int nbr_wtband() {return NbrBand;}
	  // return the number of wavelet bands
	  
	  int nbr_scale() {return Nbr_Plan;}
	  // return the number dyadic scales
	  
          int size_band_nl(int b) {return TabNl(b);}
	  // return the number of lines of the band b
	  
	  int size_band_nc(int b) {return TabNc(b);}
	  // return the number of columns of the band b
	  
	  void get_band(Ifloat &Band, int b);
	  // extract the band b
	  
	  int get_resol_band(int b) {return TabResolBand(b);}
	  // return the resolution level of a given scale
	  //        it is in [-NbrPlan,+NbrPlan]
	  //        -NbrPlan corresponds to the highest frequency.
	  //        NbrPlan corresponds to the lowest frequency.
	  
	  int get_indwt_band(int b) {return TabIndBand(b);}
	  // Return the index in the wavelet representation of
	  //        a given band
	  
	  void put_band(Ifloat &Band, int b);
	  // put the band b
	  
          float & operator() (int b, int i, int j);
	  // reference to a mirror wavelet coefficient

	  void get_tab_mirrorband(Ifloat * & TabMirrorBand);
          // creates an array of band, and copy the transformed data
	  // in this array:  TabMirrorBand[0..NbrTotBand-1]
	  // The size of a given and is given by  size_band_nl
	  // and size_band_nc
	  
	  void put_tab_mirrorband(Ifloat * & TabMirrorBand);
	  // set the array TabTransWT with the transformed data
	  // TabMirrorBand[0..NbrTotBand-1]
	  
	  // ================================================
          // The following routines allows the manipulation
	  // of a data buffer given by the calling the routine	  
	  Ifloat * get_trans() {return TabTransWT;}
	  // get a pointer to the data: TabTransWT[0..NbrBand-1]
	  
	  void get_band(Ifloat *TabTrans, Ifloat &Band, int b);
	  // extract the band b from TabTrans
	  
          void put_band(Ifloat *TabTrans, Ifloat &Band, int b);
	  // put the band b in TabTrans
	  
	  float get_val(Ifloat *TabTrans, int b, int i, int j);
	  // get coeff value from TabTrans
	  
	  void put_val(Ifloat *TabTrans, int b, int i, int j, float Val);
	  // put a coeff value  in TabTrans
	  
          void get_tab_mirrorband(Ifloat *TabBand, Ifloat * & TabMirrorBand);
	  // creates an array of band from TabBand and copy the transformed data
	  // in this array:  TabMirrorBand[0..NbrTotBand-1]
	  // The size of a given and is given by  size_band_nl
	  // and size_band_nc
	  
          void put_tab_mirrorband(Ifloat *TabBand, Ifloat * & TabMirrorBand);
	  // set the  array  TabBand[0..NbrBand-1] with the transformed data
	  // in this array TabMirrorBand[0..NbrTotBand-1]
 	  
          ~MIRROR_2D_WT(){Ptr_SB1D = NULL; delete Ptr_LC;
	                 delete Ptr_HD; delete Ptr_OWT;}
};

/****************************************************************************/
 
#endif
