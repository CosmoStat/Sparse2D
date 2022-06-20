/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  27/07/98
**    
**    File:  MR_PCA.h
**
*******************************************************************************
**
**    DESCRIPTION  do a principal component analysis on scales of a MR transform  
**    ----------- 
**
***************************************************************************/

#ifndef _MRPCA_H_
#define _MRPCA_H_

#include "MR_NoiseModel.h"
#include "CPca.h"


#define MAX_NBR_PCA_IMA 500
#define MAX_NB_BAND 8
 
const int MAX_NBR_MR_PCA_IMA = 500;
typedef int Tab [MAX_NBR_PCA_IMA][MAX_NB_BAND]; 

#define MAX_CORREL_COMP 4
#define MAX_CORREL_NOISE 4
enum CorrelComputeType   {E_LOCAL=0,
			  E_GLOBAL=1,
			  E_GLOBAL_WITHOUT_LAST=2,
			  E_IMPORT=3};
enum CorrelNoiseType     {E_WITHOUT=0,
                          E_THRESHOLD=1,
                          E_LINEAR=2,
			  E_PROBA=3};
inline const char* CorCompTransform (CorrelComputeType CorrelComp) 
{
   switch (CorrelComp) 
   {
   case E_LOCAL  : return ("One correlation matrix per band"); break;
   case E_GLOBAL : return ("One correlation matrix for all bands"); break;
   case E_GLOBAL_WITHOUT_LAST : 
      return ("One correlation matrix for all band, except the last scale"); break;
   case E_IMPORT : return ("Import correlation matrix"); break;
   default : return ("correl type unknown !!!"); break;
   }
};
inline const char* CorNoiseTransform (CorrelNoiseType CorrelNoise) 
{
   switch (CorrelNoise) 
   {
   case E_WITHOUT : return ("Compute the correl. matrix with all wavelet coefficients"); break;
   case E_THRESHOLD : return ("Compute the correl. matrix only with significant wavelet coefficients"); break;
   case E_LINEAR : return ("Compute the correl. matrix with weighted wavelet coefficients"); break;
   case E_PROBA : return ("Compute the correl. matrix with weighted wavelet coefficients"); break;
   default : return ("Correl noise unknown !!!"); break;
   }
};   

class PCA_MR2D {
   int _NbrImage;        // Number of images to analyse
   int _NbrBand;         // Number of bands 
   Bool _TraceMatCor;    // write correl Mat on disk
   Bool _FirstTraceOk;   // trace only the first time
public:
   int nbr_image() {return _NbrImage;}
   int nbr_band()  {return  _NbrBand;}
   
   CPCA         *MatCor;        // class for principal component analysis (one per band)
   fltarray     *CorrelMat;     // correlation matrix (2D) (one per band)
   fltarray     TabMean;       // mean value of each image (2D array) (one per band) 			  
  
   PCA_MR2D (); 
   PCA_MR2D (int Nima, int Nband, Bool TraceMatCor=False); 
   void alloc (int Nima, int Nband, Bool TraceMatCor=False);
   
   void print(); // print the analysis results (correl. matrix,eigen vectors, eigen values)  
   void print(int Band); // print the analysis results for one band
   
private:
   void subtract_mean(MultiResol *TabMR);   
   void add_mean(MultiResol *TabMR);    
   // calculate the mean of each image in TabIma,
   // store it in TabMean and subtract it   
   // add to each image k the mean value TabMean 

   float Correl_Without_Normalisation (const Ifloat & Im1, 
                                                 const Ifloat & Im2, 
						 float ValMin);
   // compute correlation coef without normalisation

public:     
   void compute(MultiResol *TabMR, MRNoiseModel *TabMRNoiseModel=NULL, 
                Bool NormCorrelMat = False,
		CorrelComputeType ComputeType = E_LOCAL, 
		CorrelNoiseType NoiseType = E_WITHOUT,
		Ifloat * ImportCorrelMat = NULL, float ValMin=0);
   void transform (MultiResol  *TabMRin,  MultiResol *TabMRout, 
                   MRNoiseModel *TabNoiseIn=NULL, MRNoiseModel *TabNoiseOut=NULL);
   void invtransform (MultiResol *TabMRin , MultiResol *TabMRout);
   void invsubtransform (MultiResol *TabMRin , MultiResol *TabMRout, 
                         int *NbrEigen=NULL, 
                         int TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND]=NULL);
   ~PCA_MR2D ();
  
public:
   void pca_Thresold (MultiResol  *TabMRin_out,MRNoiseModel *TabNoiseOut);
   void pca_TransfSignal (MultiResol *TabMRin , MultiResol *TabMRout, MRNoiseModel *TabNoiseIn = NULL);
   void pca_ComputeNoise  (MultiResol *TabMRin, MRNoiseModel *TabNoiseIn, MRNoiseModel *TabNoiseOut); 
};


  
  
  
  
  


/*  void compute(MultiResol *TabMR, 
                MRNoiseModel *TabMRNoiseModel=NULL, 
                float ValMin=0);
    // For each scale:
    //    calulate the correlation matrix and diagonalize it
    //    find the eigen values and the eigen vectors  
    //    If UseNoiseModel == True, the correlation matrix is 
    //        calculated only from the significant wavelet coefficients
    //        TabMR[i](b,x,y) is significant if 
    //            | TabMR[i](b,x,y) | >  TabNsigma(i,b) * Sigma
    //           with Sigma =  TabNoise(i) * TabMR[i].band_norm(b)
    //    EigenNoise is correctly initialized.
    //
    // before calling this routine, TabNoise(*) must be initiazied!
    
  //void transform (MultiResol  *TabMRin,  MultiResol *TabMRout);
  void transform (MultiResol  *TabMRin,  MultiResol *TabMRout,
                  MRNoiseModel *TabNoiseIn=NULL);
  // transform an MR object array into its principal component
  // (TabMRout  can be set to the same pointer as TabMRin, then
  // (the input data will be overwritten, but memory will be safed).
  
  void  invtransform (MultiResol *TabMRin , MultiResol *TabMRout);
  // apply an inverse transform: the MR objects are reconstructed from
  // the eigen vector.
  // (TabMROut  can be set to the same pointer as TabMRin, then
  // (the input data will be overwritten, but memory will be safed).
  
  void  invsubtransform (MultiResol *TabMRin , MultiResol *TabMRout, 
                         int *NbrEigen=NULL, 
                         int TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND]=NULL);
  // apply an inverse transform: the images are reconstructed from
  // a sub set of the eigen vector images.
  // NbrEigen[b] : input = number of eigen values to keep at scale b
  // TabKill[0..NbrImage-1][0..NbrBand-1] = 1 or 0
  //     if TabKill[i][b] = 0, eigen vector number i, at scale b is not 
  //     taken into account for the reconstruction

  void thresold_eigen (MultiResol  *TabMRin_out);
  // threshold all values in TabMRin_out such that:
  //   | TabMRin_out[i](b,x,y) | < TabNsigma(i,b) * EigenNoise(i,b) 
*/

          
#endif
