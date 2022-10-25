

#ifndef _MR3DFEWEVENT_H_
#define _MR3DFEWEVENT_H_

#include "MR1D_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "IM_Noise.h" 
#include <iomanip>


// for wavelet psi
#define BIN3D 150 	
#define BIN3D_2 2*BIN3D+1

//sigma wavelet 3D
#define SIGMA_BSPLINE_WAVELET 0.0102225

// length max of autoconvolution
#define MAXHISTONBPOINT  2048+1

// length of first histo (construct without convolution)
#define LENGTH_FIRST_HISTO 1024+1

//  output names
#define Name_Bspline_   "Aba_bspline"
#define Name_Wavelet_   "Aba_wavelet"
#define Name_HistoConv  "_Aba_histo"
#define Name_Threshold "_threshold"
#define Name_HistoBin "_histobin"
#define Name_HistoBound "_histobound"
#define Name_HistoDistrib "_histodistrib"
#define Name_Param "_param"


#define	CUBE_FABS(x)	(fabs(x) * fabs(x) * fabs(x))
#define CUBE(x)		((x) * (x) * (x))
#define DEBUG_FEW 1     


//definition of differents values of the support
#define VAL_SupNull 0
#define VAL_SupOK 1
#define VAL_SupDill 2
#define VAL_SupEdge 3
#define VAL_SupLastOK 9
#define VAL_SupKill 10
#define VAL_SupMinEv 11
#define VAL_SupFirstScale 12  


class FewEventPoisson {

private:
   float _Epsilon;	     //epsilon for threshold
   int  _NbAutoConv;         // number max of autoconv used
   bool _InitOk;             // initialisation flag
   
   // Computes the one dimensional bspline  
   // Computes the histogram of the wavelet in histo(0,*)
   void bspline_histo_3D (Bool WriteAllInfo=False, Bool Verbose=False);  
   
   // performs auto-convolution of the first histogram. 
   void histo_convolution (Bool WriteAllInfo=False, Bool Verbose=False, int param=0);
   
   // Normalizes the histograms 
   //void histo_normalisation (Bool WriteAllInfo=False, Bool Verbose=False);
   
   // Computes the repartition function 
   void histo_distribution (Bool WriteAllInfo=False, Bool Verbose=False);
   
   void shape_signal (fltarray &Signal1d);
   
   void show_param (char* text, float nbconv, float min, float max, 
                    float bin, float nbechant);
     
   int get_ev(int i,int j,int k,intarray & nb_event, type_border Border);
   
public:
   fltarray _Param;	      // _Param(0) =value of rytme
   			      // _Param(1) =value of _NbAutoConv
   fltarray _HistoBound;      // min and max reduced coefficient
   fltarray _HistoBin;        // bin and number of points  
   fltarray _HistoConv;       // Histogram information on autoconv
   fltarray _HistoDistrib;    // Distribution of Histograms
   Ifloat   _Threshold;      // Threshold(p, 0) = Min value threshold
                             // Threshold(p, 1) = Max value threshold
  
   // init attribute
   FewEventPoisson(Bool InitOk, int AutoConv,float epsilon);
  

   // Computes the histogram autoconvolutions, and the repartition function.
   void compute_distribution (Bool WriteAllInfo, Bool WriteHisto, Bool
   		Verbose, int param);

   // set the Threshold array for a confidence interval
   // given by Epsilon  
   void find_threshold (Bool WriteAllInfo,Bool WriteHisto,
   			Bool Verbose, int param);
   
   // set the Threshold array for a confidence interval
   // given by Epsilon  
   void histo_threshold (Bool WriteAllInfo=False, Bool
   	Verbose=False);

   void get_threshold(int NbrEvent, int CurrentScale, float & SeuilMin, float & SeuilMax);
   // Return the threshold min and max corresponding the number of events
   // NbrEvent

   void event_set_support(	fltarray *&	Data_In,
		        int             CurrentScale,
			int		FirstDectectScale, 
			int		NEGFirstDectectScale, 
			int 		MinEvent,
			Bool		OnlyPosVal,
		        type_border     Border,
		        intarray *& 	Nb_Event,
		        intarray *& 	Support,
		        Bool            WriteAllInfo);
   
   
   // give the value of a pixel depending to the border
   int get_pix(int i,int j,int k, intarray & data, type_border Border);
 
   // compute nb_event with no approximation
   void event_scale_3D (intarray & data_event, intarray *& nb_event,
	int Current_Scale, type_border Border, Bool WriteAllInfo) ;
		
   // compute nb_event with no approximation
   void event_scale_3D_Approx (	intarray & last_nb_event, 
				intarray & nb_event,
				int Current_Scale, 
				type_border Border,
				Bool WriteAllInfo) ;

   void get_sigma_band( fltarray & Nb_Event_In_Boxes,fltarray & tab_sigma,
   		float N_Sigma,int Nb_Scale);
};


#endif
