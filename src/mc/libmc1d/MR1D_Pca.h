/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: %V%
**
**    Author: 08/23/99
**
**    Date:  1.1
**    
**    File:  MR1D_Pca.h
**
*******************************************************************************
**
**    DESCRIPTION  do a principal component analysis on scales of a MR transform  
**    ----------- 
**
***************************************************************************/

#ifndef _MR1DPCA_H_
#define _MR1DPCA_H_


#include "MR1D_NoiseModel.h"
#include "CPca.h"
#include "mc1d_com.h"



#define MAX_NBR_PCA_SPEC 500
#define MAX_NB_BAND 8

const int MAX_NBR_MR_PCA_SPEC = 500;
typedef int Tab [MAX_NBR_PCA_SPEC][MAX_NB_BAND]; 

#define MAX_CORREL_COMP 4
#define MAX_CORREL_NOISE 4

			  

			  
class PCA_MR1D {
   int _NbSpectre;       // Number of images to analyse
   int _NbBand;          // Number of bands 
   Bool _TraceMatCor;    // write correl Mat on disk
   Bool _FirstTraceOk;   // trace only the first time  
public:
   int nbr_spectre() {return _NbSpectre;}
   int nbr_band()  {return  _NbBand;}
   
   CPCA         *MatCor;        // class for principal component analysis (one per band)
   fltarray     *CorrelMat;     // correlation matrix (2D) (one per band)			  
   fltarray     TabMean;       // mean value of each image (2D array) (one per band) 			  
  
   PCA_MR1D (); 
   PCA_MR1D (int NbSpectre, int NbBand, Bool TraceMatCor=False); 
   void alloc (int NbSpectre, int NbBand, Bool TraceMatCor=False);
   
   void print(); // print the analysis results (correl. matrix,eigen vectors, eigen values)  
   void print(int Band); // print the analysis results for one band
 
private:   
   void subtract_mean(MR_1D *TabMR);   
   void add_mean(MR_1D *TabMR);    
   // calculate the mean of each image in TabIma,
   // store it in TabMean and subtract it   
   // add to each image k the mean value TabMean 

   float Correl_Without_Normalisation (const fltarray & Sp1, 
                                       const fltarray& Sp2, 
				       float ValMin);
	       
	   
public:
   void compute(MR_1D *TabMR, MR1DNoiseModel *TabMRNoiseModel=NULL, 
                Bool NormCorrelMat = False,
		CorrelComputeType1d ComputeType = E_LOCAL1D, 
		CorrelNoiseType1d NoiseType = E_WITHOUT1D,
		Ifloat * ImportCorrelMat = NULL, float ValMin=0);
   void transform (MR_1D *TabMRin, MR_1D *TabMRout, 
                   MR1DNoiseModel *TabNoiseIn=NULL, MR1DNoiseModel *TabNoiseOut=NULL);
   
   //void invtransform (MR_1D *TabMRin , MR_1D *TabMRout);
   void invsubtransform (MR_1D *TabMRin , MR_1D *TabMRout,
                         int *NbrEigen=NULL, 
                         int TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND]=NULL);
   ~PCA_MR1D ();
  
public:
   void pca_Thresold (MR_1D *TabMRin_out, MR1DNoiseModel *TabNoiseOut);
   void pca_TransfSignal (MR_1D *TabMRin , MR_1D *TabMRout, MR1DNoiseModel *TabNoiseIn = NULL);
   void pca_ComputeNoise  (MR_1D *TabMRin, MR1DNoiseModel *TabNoiseIn, MR1DNoiseModel *TabNoiseOut);
};


  

          
#endif
