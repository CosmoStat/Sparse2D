/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.2
**
**    Author: 12/07/98
**
**    Date:  98/12/07
**    
**    File:  MR_PCA.cc
**
*******************************************************************************
**
**    DESCRIPTION  class for principal component analysis for images  
**    ----------- 
**                 
******************************************************************************/
 
#include "IM_Obj.h"
#include "MR_Obj.h"
#include "NR.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "IM_Noise.h"
#include "MR_NoiseModel.h"

#include "MR_PCA.h"

#define DEF_TYPENOISE NOISE_GAUSSIAN



PCA_MR2D::PCA_MR2D () {}

PCA_MR2D::PCA_MR2D (int Nima, int Nband, Bool TraceMatCor) {
        _TraceMatCor = TraceMatCor;
	alloc (Nima, Nband, _TraceMatCor);
}

void PCA_MR2D::alloc (int Nima, int Nband, Bool TraceMatCor) {
  
  int b;
  _NbrImage = Nima;		      // number of images to analyse
  _NbrBand = Nband;                    // number of band
  _TraceMatCor = TraceMatCor;
  _FirstTraceOk = False;
  MatCor = new CPCA[_NbrBand];         // class for principal component analysis
  for (b=0;b<_NbrBand;b++) {
     MatCor[b].init();
     MatCor[b].alloc(_NbrImage);
  }
  CorrelMat = new fltarray[_NbrBand];  // correlation matrix (2D)
  for (b=0;b<_NbrBand;b++) CorrelMat[b].alloc(_NbrImage,_NbrImage);
  TabMean.alloc(_NbrBand,_NbrImage);
  
}


void PCA_MR2D::subtract_mean(MultiResol *TabMR) {
    
   for(int k=0; k <_NbrImage; k++) {  
      for (int b=0;b<_NbrBand;b++) {    
       
         Ifloat ao_Band = TabMR[k].extract_band(b);
         TabMean(b,k) = average (ao_Band);
	 //cout << "!!!!! Mean(b:"<<b<<",i:"<<k<<"): " << TabMean(b,k) << endl;
         for (int l=0; l<ao_Band.nl(); l++) {
            for (int n=0; n<ao_Band.nc(); n++) {
	       (TabMR[k])(b,l,n) -= TabMean(b,k);
	    }
	 }
      }
   }
}

void PCA_MR2D::add_mean(MultiResol *TabMR) {

    for(int k=0; k <_NbrImage; k++) {    
       for (int b=0;b<_NbrBand;b++) {    
       
         Ifloat ao_Band = TabMR[k].extract_band(b);
         for (int l=0; l<ao_Band.nl(); l++) {
            for (int n=0; n<ao_Band.nc(); n++) {
	       (TabMR[k])(b,l,n) += TabMean(b,k);
	    }
	 }
      }
   }
}


/*void PCA_MR2D::compute(MultiResol *TabMR, float ValMin) {

	    // For each scale:
    //    calulate the correlation matrix and diagonalize it
    //    find the eigen values and the eigen vectors  
    //    If UseNoiseModel == True, the correlation matrix is 
    
  if (ValMin < 0.) ValMin = 0.;
  for (int b=0;b<_NbrBand;b++) {

    // compute the correlation matrix
    CorrelMat[b].init();
        
    // diagonal elements are equal to 1
    for(int i=0; i < _NbrImage; i++) CorrelMat[b](i,i) = 1.;
       
    for(int i=0; i<_NbrImage; i++) {
      for(int j=i+1; j<_NbrImage; j++) { 
         
          CorrelMat[b](i,j) = correlation(TabMR[i].extract_band(b), 
                                          TabMR[j].extract_band(b), ValMin);
 	  
        // the matrix is symetric  
        CorrelMat[b](j,i) =  CorrelMat[b](i,j);
      }
    }
       
    // compute the eigen values
    MatCor[b].compute_eigen(CorrelMat[b]);
  }
}  */          



void PCA_MR2D::compute (MultiResol *TabMR, MRNoiseModel *TabMRNoise,
                        Bool NormCorrelMat,
                        CorrelComputeType ComputeType, 
			CorrelNoiseType NoiseType, 
			Ifloat* ImportCorrelMat, float ValMin) {

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
    
  subtract_mean (TabMR);
    
  if (ValMin < 0.) ValMin = 0.;
  
  int GlobalSize = 0;
  
  switch (ComputeType) {
  
  case E_LOCAL: 
  
    for (int b=0;b<_NbrBand;b++) {

      // compute the correlation matrix
      CorrelMat[b].init();
     
      int i;   
      // diagonal elements are equal to 1
//      for(i=0; i < _NbrImage; i++) CorrelMat[b](i,i) = 1.;
      
      for(i=0; i<_NbrImage; i++) {
//        for(int j=i+1; j<_NbrImage; j++) { 
         for(int j=i; j<_NbrImage; j++) {       // on calcule le coef diagonal

	  Ifloat ao_BandI = TabMR[i].extract_band(b);
          Ifloat ao_BandJ = TabMR[j].extract_band(b);
	  
	  switch (NoiseType) {
	  
	  case E_WITHOUT: break;
	  
	  case E_LINEAR:
             if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
                for (int l=0; l<ao_BandI.nl(); l++) {
                   for (int n=0; n<ao_BandI.nc(); n++) {
	              float af_Coef, af_ksigma;
		      af_ksigma= TabMRNoise[i].NSigma[b] * TabMRNoise[i].sigma(b,l,n);
		      af_Coef = (ao_BandI(l,n)/af_ksigma >= 1 ? 1 : ao_BandI(l,n)/af_ksigma);
		      ao_BandI(l,n) *= af_Coef;
		      af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l,n);
		      af_Coef = (ao_BandJ(l,n)/af_ksigma >= 1 ? 1 : ao_BandJ(l,n)/af_ksigma);
	              ao_BandJ(l,n) *= af_Coef;
		   }
		}
             }
	     break;
	  
	  case E_PROBA:
             if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
                for (int l=0; l<ao_BandI.nl(); l++) {
                   for (int n=0; n<ao_BandI.nc(); n++) {
		      ao_BandI(l,n) *= (1-TabMRNoise[i].prob_noise(ao_BandI(l,n),b,l,n));
	              ao_BandJ(l,n) *= (1-TabMRNoise[j].prob_noise(ao_BandJ(l,n),b,l,n));
		   }
		}
             }
	     break;
	     		      
	  case E_THRESHOLD:      	      
             if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
                for (int l=0; l<ao_BandI.nl(); l++) {
                   for (int n=0; n<ao_BandI.nc(); n++) {		      
	              if (TabMRNoise[i](b,l,n) != True) ao_BandI(l,n) = 0.;
	              if (TabMRNoise[j](b,l,n) != True) ao_BandJ(l,n) = 0.;
                   }
                }
             }
	     break;
          }

          if (!NormCorrelMat) CorrelMat[b](i,j) = Correl_Without_Normalisation (ao_BandI, ao_BandJ, ValMin);
          else CorrelMat[b](i,j) = correlation (ao_BandI, ao_BandJ, ValMin);

          // the matrix is symetric  
          CorrelMat[b](j,i) =  CorrelMat[b](i,j);
        }
      }
       
      // compute the eigen values
      MatCor[b].compute_eigen(CorrelMat[b]);
    }
    break;
  
  case E_GLOBAL:
  
     for (int b=0;b<_NbrBand;b++) {
        GlobalSize +=  TabMR[0].nbr_coeff_in_band(b); 
	        
        // compute the correlation matrix
        CorrelMat[b].init();
	 
        // diagonal elements are equal to 1
//        for(int i=0; i < _NbrImage; i++) CorrelMat[b](i,i) = 1.;	 
     }
      
     for(int i=0; i<_NbrImage; i++) {
//        for(int j=i+1; j<_NbrImage; j++) { 
         for(int j=i; j<_NbrImage; j++) {       // on calcule le coef diagonal

	   //Ifloat ao_BandI = TabMR[i].extract_band(b);
           //Ifloat ao_BandJ = TabMR[j].extract_band(b);
	   Ifloat ao_GlobalI (1, GlobalSize, "");
	   Ifloat ao_GlobalJ (1, GlobalSize, "");
	  	  
           int CurrentInd=0;
	   for (int b=0;b<_NbrBand;b++) {
	  
	      for (int l=0; l<TabMR[0].size_band_nl(b); l++) {
	         for (int n=0; n<TabMR[0].size_band_nc(b); n++) {
	
	            ao_GlobalI (0, CurrentInd) = (TabMR[i])(b,l,n);
		    ao_GlobalJ (0, CurrentInd) = (TabMR[j])(b,l,n);    
	            
		    switch (NoiseType) {
		  
	            case E_WITHOUT: break;
	  
	            case E_LINEAR:
		  
	       	       if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
		      
		          float af_Coef, af_ksigma;
		          af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l,n);
		          af_Coef = (ao_GlobalI(0, CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalI(0, CurrentInd)/af_ksigma);
		          ao_GlobalI(0, CurrentInd) *= af_Coef;
		          af_Coef = (ao_GlobalJ(0, CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalJ(0, CurrentInd)/af_ksigma);
	                  ao_GlobalJ(0, CurrentInd) *= af_Coef;
	               } 
		       break; 
	  
	            case E_PROBA:
                      
		      if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
		      
		         float Val = ao_GlobalI(0, CurrentInd);
		         ao_GlobalI(0, CurrentInd) *= (1-TabMRNoise[i].prob_noise(Val,b,l,n));
	                 Val = ao_GlobalJ(0, CurrentInd);
			 ao_GlobalJ(0, CurrentInd) *= (1-TabMRNoise[j].prob_noise(Val,b,l,n));
		      }
	              break;
	     			          
	            case E_THRESHOLD:
	       	     
		       if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
		
	                  if (TabMRNoise[i](b,l,n) != True) ao_GlobalI(0, CurrentInd) = 0.;
	                  if (TabMRNoise[j](b,l,n) != True) ao_GlobalJ(0, CurrentInd) = 0.;
                       }
		       break;
		  
		    }
		    CurrentInd++;
	         }
              }
	   }
	    
  
           if (!NormCorrelMat) CorrelMat[0](i,j) = Correl_Without_Normalisation (ao_GlobalI, ao_GlobalJ, ValMin);
           else CorrelMat[0](i,j) = correlation (ao_GlobalI, ao_GlobalJ, ValMin);

           // the matrix is symetric  
           CorrelMat[0](j,i) =  CorrelMat[0](i,j);
        }
     } 
      
     // compute the eigen values
     MatCor[0].compute_eigen(CorrelMat[0]);       
      
     for (int b=1;b<_NbrBand;b++) {
     
        for(int i=0; i<_NbrImage; i++) {
//           for(int j=i+1; j<_NbrImage; j++) {   
           for(int j=i; j<_NbrImage; j++) {       // on calcule le coef diagonal
    
              CorrelMat[b](i,i) =  CorrelMat[0](i,i);
	      CorrelMat[b](i,j) =  CorrelMat[0](i,j);
	      CorrelMat[b](j,i) =  CorrelMat[b](i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }
     break;
     
  case E_GLOBAL_WITHOUT_LAST:
  
     for (int b=0;b<_NbrBand-1;b++) { // without last scale
        GlobalSize +=  TabMR[0].nbr_coeff_in_band(b); 
	        
        // compute the correlation matrix
        CorrelMat[b].init();
	 
        // diagonal elements are equal to 1
//        for(int i=0; i < _NbrImage; i++) CorrelMat[b](i,i) = 1.;	 
     }
      
     for(int i=0; i<_NbrImage; i++) {
//        for(int j=i+1; j<_NbrImage; j++) { 
        for(int j=i; j<_NbrImage; j++) {       // on calcule le coef diagonal

	   //Ifloat ao_BandI = TabMR[i].extract_band(b);
           //Ifloat ao_BandJ = TabMR[j].extract_band(b);
	   Ifloat ao_GlobalI (1, GlobalSize, "");
	   Ifloat ao_GlobalJ (1, GlobalSize, "");
	  	  
           int CurrentInd=0;
	   for (int b=0;b<_NbrBand-1;b++) { // without last scale
	  
	      for (int l=0; l<TabMR[0].size_band_nl(b); l++) {
	         for (int n=0; n<TabMR[0].size_band_nc(b); n++) {
	
	            ao_GlobalI (0, CurrentInd) = (TabMR[i])(b,l,n);
		    ao_GlobalJ (0, CurrentInd) = (TabMR[j])(b,l,n);    
	            
		    switch (NoiseType) {
		  
	            case E_WITHOUT: break;
	  
	            case E_LINEAR:
		  
	       	       if (TabMRNoise != (MRNoiseModel*)NULL) {
		      
		          float af_Coef, af_ksigma;
		          af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l,n);
		          af_Coef = (ao_GlobalI(0, CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalI(0, CurrentInd)/af_ksigma);
		          ao_GlobalI(0, CurrentInd) *= af_Coef;
		          af_Coef = (ao_GlobalJ(0, CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalJ(0, CurrentInd)/af_ksigma);
	                  ao_GlobalJ(0, CurrentInd) *= af_Coef;
	               } 
		       break; 
	  
	            case E_PROBA:
                      
		      if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbrBand-1) {
		      
		         float Val = ao_GlobalI(0, CurrentInd);
		         ao_GlobalI(0, CurrentInd) *= (1-TabMRNoise[i].prob_noise(Val,b,l,n));
	                 Val = ao_GlobalJ(0, CurrentInd);
			 ao_GlobalJ(0, CurrentInd) *= (1-TabMRNoise[j].prob_noise(Val,b,l,n));
		      }
	              break;
	     			          		          
	            case E_THRESHOLD:
	       	     
		       if (TabMRNoise != (MRNoiseModel*)NULL) {
		
	                  if (TabMRNoise[i](b,l,n) != True) ao_GlobalI(0, CurrentInd) = 0.;
	                  if (TabMRNoise[j](b,l,n) != True) ao_GlobalJ(0, CurrentInd) = 0.;
                       }
		       break;
		  
		    }
		    CurrentInd++;
	         }
              }
	   }
	    
  
           if (!NormCorrelMat) CorrelMat[0](i,j) = Correl_Without_Normalisation (ao_GlobalI, ao_GlobalJ, ValMin);
           else CorrelMat[0](i,j) = correlation (ao_GlobalI, ao_GlobalJ, ValMin);

           // the matrix is symetric  
           CorrelMat[0](j,i) =  CorrelMat[0](i,j);
        }
     } 
      
     // compute the eigen values
     MatCor[0].compute_eigen(CorrelMat[0]);       
      
     for (int b=1;b<_NbrBand;b++) {
     
        for(int i=0; i<_NbrImage; i++) {
//           for(int j=i+1; j<_NbrImage; j++) {       
           for(int j=i; j<_NbrImage; j++) {       // on calcule le coef diagonal
              
	      CorrelMat[b](i,i) =  CorrelMat[0](i,i);
	      CorrelMat[b](i,j) =  CorrelMat[0](i,j);
	      CorrelMat[b](j,i) =  CorrelMat[b](i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }
     break;
       
  case E_IMPORT:
   
     for (int b=0;b<_NbrBand;b++) {
     
        for(int i=0; i<_NbrImage; i++) {
           for(int j=0; j<_NbrImage; j++) {       
              CorrelMat[b](i,j) =  (*ImportCorrelMat)(i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }  
  
     break;
  }
  if (_TraceMatCor) {
     fltarray CM = fltarray (_NbrImage,_NbrImage,_NbrBand);
     for (int b=0;b<_NbrBand;b++) {
        for(int i=0; i<_NbrImage; i++) {
	   for(int j=0; j<_NbrImage; j++) {
	      CM(i,j,b) = CorrelMat[b](i,j);
	   }
	} 
     }
     fits_write_fltarr ("Correl_Matrix.fits", CM);
  }
} 


void PCA_MR2D::transform (MultiResol  *TabMRin,  MultiResol *TabMRout,
                          MRNoiseModel *TabNoiseIn, MRNoiseModel *TabNoiseOut) {

   if (TabNoiseIn != (MRNoiseModel*)NULL) {
 
      pca_TransfSignal (TabMRin, TabMRout, NULL);
      pca_ComputeNoise (TabMRin, TabNoiseIn, TabNoiseOut); 
      pca_Thresold (TabMRout, TabNoiseOut);
   }
 
   else pca_TransfSignal (TabMRin, TabMRout);
}




void PCA_MR2D::invtransform (MultiResol *TabMRin , MultiResol *TabMRout) {

  fltarray Vect(_NbrImage);
  fltarray Vectout(_NbrImage);
  
  for (int b=0;b<_NbrBand-1;b++) {
    for (int i=0; i<TabMRin[0].size_band_nl(b); i++) {
      for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
      
        int k;
        for (k=0; k<_NbrImage; k++) Vect(k) = (TabMRin[k])(b,i,j);
        MatCor[b].transform(Vect, Vectout);
        for (k = 0; k<_NbrImage; k++) (TabMRout[k])(b,i,j) = Vectout(k);
      }
    }
  } 
  add_mean (TabMRout);
}


void  PCA_MR2D::invsubtransform (MultiResol *TabMRin , MultiResol *TabMRout,
                                 int *NbrEigen, 
                                 int TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND]){
 
  fltarray Vect(_NbrImage);
  fltarray Vectout(_NbrImage);
  
  int Ne=0;
  //cout << "Number of eigen vector Ne = " 
  //     << NbrEigen[0] << " (same for all band)" << endl;
  char list[80]="";
  for (int l=0;l<NbrEigen[0];l++)
     if (TabKill[l][0] == 1)
        sprintf (list, "%s %d", list, l+1);
     else
        Ne++;
  if (!_FirstTraceOk) {
     cout << "Number of eigen vector not used :" << list << endl;
     cout << "Number of eigen vector used for the reconstruction: Ne = " 
          << Ne << " (same for all band)" << endl; 
     _FirstTraceOk = True;
  } 
  
  
  for (int b=0;b<_NbrBand-1;b++) {

    int N = NbrEigen[b];
    if ((N < 1) || (N > _NbrImage)) N = _NbrImage;
    
    for (int i=0; i < TabMRin[0].size_band_nl(b); i++) {
      for (int j=0; j < TabMRin[0].size_band_nc(b); j++) {

        int k;
        for (k = 0; k < N; k++) {
          if ((TabKill != NULL) && (TabKill[k][b] == 1)) Vect(k) = 0.;
          else Vect(k) = TabMRin[k](b,i,j);
        }	
        MatCor[b].invtransform(Vect, Vectout, N);
        for (k = 0; k < _NbrImage; k++) TabMRout[k](b,i,j) = Vectout(k);
      }
    }   
  }
  add_mean (TabMRout);
}

  

PCA_MR2D::~PCA_MR2D () {

  for (int b=0;b<_NbrBand;b++) CorrelMat[b].free();
  delete [] MatCor;       
  delete [] CorrelMat;
}

   
void PCA_MR2D::print(int Band) {
    cout << "band " << Band+1 << " :" << endl;
    MatCor[Band].print();
    cout << endl;
}
 

void PCA_MR2D::print() {
  for (int b=0;b<_NbrBand;b++) print(b);
}


void PCA_MR2D::pca_TransfSignal (MultiResol *TabMRin , MultiResol *TabMRout,
                                 MRNoiseModel *TabNoiseIn) {

  fltarray Vect(_NbrImage);
  fltarray Vectout(_NbrImage);     
     
  for (int b=0;b<_NbrBand-1;b++) {
  
    for (int i=0; i<TabMRin[0].size_band_nl(b); i++) {
      for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
 
        int k;
        for (k=0; k<_NbrImage; k++) Vect(k) = (TabMRin[k])(b,i,j);
	   /*if (TabNoiseIn != (MRNoiseModel*)NULL && b != _NbrBand-1 )
              Vect(k) = (((TabNoiseIn[k])(b,i,j)==true) ?  
                                                   (TabMRin[k])(b,i,j) : 0.0);
           else Vect(k) = (TabMRin[k])(b,i,j);*/
	   
  
        MatCor[b].transform(Vect, Vectout);
        for (k=0; k<_NbrImage; k++) (TabMRout[k])(b,i,j) = Vectout(k);
      }
    }
  } 
} 


void PCA_MR2D::pca_ComputeNoise  (MultiResol *TabMRin,
                                  MRNoiseModel *TabNoiseIn,
				  MRNoiseModel *TabNoiseOut) {
       
       

  type_noise ae_TypeNoise = TabNoiseIn[0].which_noise();
  	
  if (! one_level_per_pos_2d (ae_TypeNoise)) {

    fltarray Vect(_NbrImage);
    fltarray Vectout(_NbrImage);     

    for (int b=0;b<_NbrBand-1;b++) {

      int k;
      for (k=0; k<_NbrImage; k++) Vect(k) = pow((double)  TabNoiseIn[k].sigma(b,0,0),2.);
      //for (int k=0; k<_NbrImage; k++) cout <<"sigma("<<b<<","<<k<<") in:"<<TabNoiseIn[k].sigma(b)<<endl;
      MatCor[b].transform_noise(Vect, Vectout);
      for (k=0; k<_NbrImage; k++) TabNoiseOut[k].sigma(b,0,0) = sqrt(Vectout(k));
      //for (int k=0; k<_NbrImage; k++) cout <<"sigma("<<b<<","<<k<<") out:"<<TabNoiseIn[k].sigma(b)<<endl;      

    }

  } else { // one sigma at (b,i,j)

    fltarray Vect(_NbrImage);
    fltarray Vectout(_NbrImage);     
     
    for (int b=0;b<_NbrBand-1;b++) {
      for (int i=0; i<TabMRin[0].size_band_nl(b); i++) {
        for (int j=0; j<TabMRin[0].size_band_nc(b); j++) {
 
          int k;
          for (k=0; k<_NbrImage; k++) Vect(k) = pow ((double) TabNoiseIn[k].sigma(b,i,j),2.);
          MatCor[b].transform_noise(Vect, Vectout);
          for (k=0; k<_NbrImage; k++) TabNoiseOut[k].sigma(b,i,j) = sqrt (Vectout(k));
        }
      }
    } 
  }
}


void PCA_MR2D::pca_Thresold (MultiResol  *TabMRin_out, MRNoiseModel *TabNoiseOut) {

   for (int k=0; k<_NbrImage; k++) 
      TabNoiseOut[k].threshold(TabMRin_out[k]);

}
          

float PCA_MR2D::Correl_Without_Normalisation (const Ifloat & Im1, const Ifloat & Im2, float ValMin)
 {
     double Sum_XY;
     int i,j;
     int Nl = Im1.nl();
     int Nc = Im1.nc();
 
    // test image size
    if ((Im1.nl() != Im2.nl()) || ( Im1.nc() != Im2.nc()))
    {
       cerr << "Error in correlation routine: images have different sizes ..." << endl ; 
       cerr << "   image 1: " <<  Im1.nl() << "X"  <<  Im1.nc() << endl ;
       cerr << "   image 2: " <<  Im2.nl() << "X"  <<  Im2.nc() << endl ;
       exit(-1);
    }
    
     Sum_XY = 0.;
     for (i = 0; i < Nl; i++) {
        for (j = 0; j < Nc; j++) {
	   if (ABS(Im1 (i,j)) > ValMin && ABS(Im2 (i,j)) > ValMin)
              Sum_XY += Im1 (i,j) * Im2 (i,j);
        }
     }
     
     return (Sum_XY/(Nl*Nc));
 }
