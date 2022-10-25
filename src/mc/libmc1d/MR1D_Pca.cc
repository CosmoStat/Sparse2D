/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: %V%
**
**    Author: %G$
**
**    Date:  1.2
**    
**    File:  MR1D_Pca.cc
**
*******************************************************************************
**
**    DESCRIPTION  class for principal component analysis for images  
**    ----------- 
**                 
******************************************************************************/
 
#include "MR1D_Obj.h"
#include "CPca.h"
#include "SP_Pca.h"
#include "IM_Noise.h"
#include "MR1D_NoiseModel.h"

#include "MR1D_Pca.h"

#define DEF_TYPENOISE NOISE_GAUSSIAN

Bool one_level_per_pos(type_noise TNoise)
// return True if we need only one level per position
// in the multiresolution space
{
    Bool ValRet= False;

    if ((TNoise == NOISE_EVENT_POISSON) || (TNoise == NOISE_NON_UNI_ADD)
         || (TNoise == NOISE_NON_UNI_MULT)  || (TNoise == NOISE_UNDEFINED))  ValRet = True;

    return ValRet;
}

PCA_MR1D::PCA_MR1D () {}

PCA_MR1D::PCA_MR1D (int NbSpectre, int NbBand, Bool TraceMatCor) {
        _TraceMatCor = TraceMatCor;
	alloc (NbSpectre, NbBand, _TraceMatCor);
}


void PCA_MR1D::alloc (int NbSpectre, int NbBand, Bool TraceMatCor) {
  
  int b;
  _NbSpectre = NbSpectre;	       // number of spectres to analyse
  _NbBand = NbBand;                    // number of band
  _TraceMatCor = TraceMatCor;
  _FirstTraceOk = False; 
  MatCor = new CPCA[_NbBand];          // class for principal component analysis
  for (b=0;b<_NbBand;b++) MatCor[b].alloc(_NbSpectre);
  CorrelMat = new fltarray[_NbBand];  // correlation matrix (2D)
  for (b=0;b<_NbBand;b++) CorrelMat[b].alloc(_NbSpectre,_NbSpectre);
  TabMean.alloc(_NbBand,_NbSpectre);
}

void PCA_MR1D::subtract_mean(MR_1D *TabMR) {
    
   for(int k=0; k <_NbSpectre; k++) {  
      for (int b=0;b<_NbBand;b++) {    
       
         fltarray ao_Band;
	 TabMR[k].scale(ao_Band, b);        
         TabMean(b,k) = ao_Band.mean();
	 //cout << "!!!!! Mean(b:"<<b<<",i:"<<k<<"): " << TabMean(b,k) << endl;
         for (int l=0; l<ao_Band.n_elem(); l++) {
	     (TabMR[k])(b,l) -= TabMean(b,k);
	 }
      }
   }
}

void PCA_MR1D::add_mean(MR_1D *TabMR) {

    for(int k=0; k <_NbSpectre; k++) {    
       for (int b=0;b<_NbBand;b++) {    
          
	 fltarray ao_Band;
	 TabMR[k].scale(ao_Band, b);           
         for (int l=0; l<ao_Band.n_elem(); l++) {
	    (TabMR[k])(b,l) += TabMean(b,k);
	 }
      }
   }
}




void PCA_MR1D::compute (MR_1D *TabMR, MR1DNoiseModel *TabMRNoise,
                        Bool NormCorrelMat,
		        CorrelComputeType1d ComputeType, 
		        CorrelNoiseType1d NoiseType,
		        Ifloat * ImportCorrelMat, float ValMin) {

  subtract_mean (TabMR);
    
  if (ValMin < 0.) ValMin = 0.;
  
  int GlobalSize = 0;
  
  switch (ComputeType) {
  
  case E_LOCAL1D: 
  
     for (int b=0;b<_NbBand;b++) {
  
        // compute the correlation matrix
        CorrelMat[b].init();
     
        int i;   
        // diagonal elements are equal to 1
        for(i=0; i<_NbSpectre; i++) 
           for(int j=i; j<_NbSpectre; j++) {       // on calcule le coef diagonal
  
              fltarray ao_BandI;
	      TabMR[i].scale(ao_BandI, b);
	      fltarray ao_BandJ;
	      TabMR[j].scale(ao_BandJ, b);

 	      switch (NoiseType) {
	  
	      case E_WITHOUT1D: break;
	  
	      case E_LINEAR1D:
                 
		 if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
                    for (int l=0; l<ao_BandI.n_elem(); l++) {
	               float af_Coef, af_ksigma;
		       af_ksigma= TabMRNoise[i].NSigma[b] * TabMRNoise[i].sigma(b,l);
		       af_Coef = (ao_BandI(l)/af_ksigma >= 1 ? 1 : ao_BandI(l)/af_ksigma);
		       ao_BandI(l) *= af_Coef;
		       af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l);
		       af_Coef = (ao_BandJ(l)/af_ksigma >= 1 ? 1 : ao_BandJ(l)/af_ksigma);
	               ao_BandJ(l) *= af_Coef;
		    }
		 }
	         break;
	  
	      //case E_PROBA1D:
             
	      //   if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
              //      for (int l=0; l<ao_BandI.n_elem(); l++) {
		//      ao_BandI(l) *= (1-TabMRNoise[i].prob_noise(ao_BandI(l),b,l));
	      //        ao_BandJ(l) *= (1-TabMRNoise[j].prob_noise(ao_BandJ(l),b,l));
		//    }
		// }
	        // break;
	     		      
	      case E_THRESHOLD1D:      	      
                 if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
                    for (int l=0; l<ao_BandI.n_elem(); l++) {
	               if (TabMRNoise[i](b,l) != True) ao_BandI(l) = 0.;
	               if (TabMRNoise[j](b,l) != True) ao_BandJ(l) = 0.;
                    }
                 }
	         break;
              } 
	      
           if (!NormCorrelMat) 
	      CorrelMat[b](i,j) = Correl_Without_Normalisation (ao_BandI, ao_BandJ, ValMin);
           else CorrelMat[b](i,j) = correlation (ao_BandI, ao_BandJ, ValMin);

           // the matrix is symetric  
           CorrelMat[b](j,i) =  CorrelMat[b](i,j);
     }
      
     // compute the eigen values
     MatCor[b].compute_eigen(CorrelMat[b]);
  }
  break; 
  
  case E_GLOBAL1D:
  
     for (int b=0;b<_NbBand;b++) {
        GlobalSize +=  TabMR[0].size_scale_np(b); 
	        
        // compute the correlation matrix
        CorrelMat[b].init();
	 
        // diagonal elements are equal to 1
//        for(int i=0; i < _NbrImage; i++) CorrelMat[b](i,i) = 1.;	 
     }
      
     for(int i=0; i<_NbSpectre; i++) {
        for(int j=i; j<_NbSpectre; j++) {       // on calcule le coef diagonal
           
           fltarray ao_GlobalI (GlobalSize);
	   fltarray ao_GlobalJ (GlobalSize);
           int CurrentInd=0;
	   
	   for (int b=0;b<_NbBand;b++) {
	  
	      for (int l=0; l<TabMR[0].size_scale_np(b); l++) {
	
	            ao_GlobalI (CurrentInd) = (TabMR[i])(b,l);
		    ao_GlobalJ (CurrentInd) = (TabMR[j])(b,l);    
	            
		    switch (NoiseType) {
		  
	            case E_WITHOUT1D: break;
	  
	            case E_LINEAR1D:
		  
	       	       if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
		      
		          float af_Coef, af_ksigma;
		          af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l);
		          af_Coef = (ao_GlobalI(CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalI(CurrentInd)/af_ksigma);
		          ao_GlobalI(CurrentInd) *= af_Coef;
		          af_Coef = (ao_GlobalJ(CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalJ(CurrentInd)/af_ksigma);
	                  ao_GlobalJ(CurrentInd) *= af_Coef;
	               } 
		       break; 
	  
	           // case E_PROBA1D:
                      
		   //    if (TabMRNoise != (MRNoiseModel*)NULL && b < _NbBand-1) {
		    //  
		   //       float Val = ao_GlobalI(CurrentInd);
		   //       ao_GlobalI(CurrentInd) *= (1-TabMRNoise[i].prob_noise(Val,b,l));
	           //       Val = ao_GlobalJ(CurrentInd);
		//	  ao_GlobalJ(CurrentInd) *= (1-TabMRNoise[j].prob_noise(Val,b,l));
		 //      }
	           //    break;
	     			          
	            case E_THRESHOLD1D:
	       	     
		       if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
		
	                  if (TabMRNoise[i](b,l) != True) ao_GlobalI(CurrentInd) = 0.;
	                  if (TabMRNoise[j](b,l) != True) ao_GlobalJ(CurrentInd) = 0.;
                       }
		       break;
		  
		    }
		    CurrentInd++;
              }
	   }
	    
  
           if (!NormCorrelMat) 
	      CorrelMat[0](i,j) = Correl_Without_Normalisation (ao_GlobalI, ao_GlobalJ, ValMin);
           else 
	      CorrelMat[0](i,j) = correlation (ao_GlobalI, ao_GlobalJ, ValMin);

           // the matrix is symetric  
           CorrelMat[0](j,i) =  CorrelMat[0](i,j);
        }
     } 
      
     // compute the eigen values
     MatCor[0].compute_eigen(CorrelMat[0]);       
      
     for (int b=1;b<_NbBand;b++) {
     
        for(int i=0; i<_NbSpectre; i++) {
           for(int j=i; j<_NbSpectre; j++) {       // on calcule le coef diagonal
    
              CorrelMat[b](i,i) =  CorrelMat[0](i,i);
	      CorrelMat[b](i,j) =  CorrelMat[0](i,j);
	      CorrelMat[b](j,i) =  CorrelMat[b](i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }
     break;
  
  case E_GLOBAL_WITHOUT_LAST1D:
  
     for (int b=0;b<_NbBand-1;b++) { // without last scale
        GlobalSize +=  TabMR[0].size_scale_np(b); 
	        
        // compute the correlation matrix
        CorrelMat[b].init();
	 
        // diagonal elements are equal to 1
//        for(int i=0; i < _NbSpectre; i++) CorrelMat[b](i,i) = 1.;	 
     }
      
     for(int i=0; i<_NbSpectre; i++) {
        for(int j=i; j<_NbSpectre; j++) {       // on calcule le coef diagonal
           
	   fltarray ao_GlobalI (GlobalSize);
	   fltarray ao_GlobalJ (GlobalSize);
           int CurrentInd=0;
	  	  
	   for (int b=0;b<_NbBand-1;b++) { // without last scale
	  
	      for (int l=0; l<TabMR[0].size_scale_np(b); l++) {
	
	         ao_GlobalI (CurrentInd) = (TabMR[i])(b,l);
		 ao_GlobalJ (CurrentInd) = (TabMR[j])(b,l);    
	            
		 switch (NoiseType) {
		  
	         case E_WITHOUT1D: break;
	  
	         case E_LINEAR1D:
		  
	       	    if (TabMRNoise != (MR1DNoiseModel*)NULL) {
		      
		       float af_Coef, af_ksigma;
		       af_ksigma= TabMRNoise[j].NSigma[b] * TabMRNoise[j].sigma(b,l);
		       af_Coef = (ao_GlobalI(CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalI(CurrentInd)/af_ksigma);
		       ao_GlobalI(CurrentInd) *= af_Coef;
		       af_Coef = (ao_GlobalJ(CurrentInd)/af_ksigma >= 1 ? 1 : ao_GlobalJ(CurrentInd)/af_ksigma);
	               ao_GlobalJ(CurrentInd) *= af_Coef;
	            } 
		    break; 
	  
	          //  case E_PROBA1D:
                      
		   //   if (TabMRNoise != (MR1DNoiseModel*)NULL && b < _NbBand-1) {
		   //   
		   //      float Val = ao_GlobalI(CurrentInd);
		   //      ao_GlobalI(CurrentInd) *= (1-TabMRNoise[i].prob_noise(Val,b,l));
	           //      Val = ao_GlobalJ(CurrentInd);
		//	 ao_GlobalJ(CurrentInd) *= (1-TabMRNoise[j].prob_noise(Val,b,l));
		   //   }
	          //    break;
	     			          		          
	            case E_THRESHOLD1D:
	       	     
		       if (TabMRNoise != (MR1DNoiseModel*)NULL) {
		
	                  if (TabMRNoise[i](b,l) != True) ao_GlobalI(CurrentInd) = 0.;
	                  if (TabMRNoise[j](b,l) != True) ao_GlobalJ(CurrentInd) = 0.;
                       }
		       break;
		  
		    }
		    CurrentInd++;
              }
	   }
	    
  
           if (!NormCorrelMat) 
	      CorrelMat[0](i,j) = Correl_Without_Normalisation (ao_GlobalI, ao_GlobalJ, ValMin);
           else 
	      CorrelMat[0](i,j) = correlation (ao_GlobalI, ao_GlobalJ, ValMin);

           // the matrix is symetric  
           CorrelMat[0](j,i) =  CorrelMat[0](i,j);
        }
     } 
      
     // compute the eigen values
     MatCor[0].compute_eigen(CorrelMat[0]);       
      
     for (int b=1;b<_NbBand;b++) {
     
        for(int i=0; i<_NbSpectre; i++) {
           for(int j=i; j<_NbSpectre; j++) {       // on calcule le coef diagonal
              
	      CorrelMat[b](i,i) =  CorrelMat[0](i,i);
	      CorrelMat[b](i,j) =  CorrelMat[0](i,j);
	      CorrelMat[b](j,i) =  CorrelMat[b](i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }
     break; 
  
  case E_IMPORT1D:
   
     for (int b=0;b<_NbBand;b++) {
     
        for(int i=0; i<_NbSpectre; i++) {
           for(int j=0; j<_NbSpectre; j++) {       
              CorrelMat[b](i,j) =  (*ImportCorrelMat)(i,j);
	   }
        }
       
        // compute the eigen values
        MatCor[b].compute_eigen(CorrelMat[b]);
     }  
  
     break;
  }
  if (_TraceMatCor) {
     fltarray CM = fltarray (_NbSpectre, _NbSpectre, _NbBand);
     for (int b=0;b<_NbBand;b++) {
        for(int i=0; i<_NbSpectre; i++) {
	   for(int j=0; j<_NbSpectre; j++) {
	      CM(i,j,b) = CorrelMat[b](i,j);
	   }
	} 
     }
     fits_write_fltarr ("Correl_Matrix.fits", CM);
  }
} 
 
  


void PCA_MR1D::transform (MR_1D *TabMRin, MR_1D *TabMRout,
                          MR1DNoiseModel *TabNoiseIn, MR1DNoiseModel *TabNoiseOut) {

   if (TabNoiseIn != (MR1DNoiseModel*)NULL) {
 
      pca_TransfSignal (TabMRin, TabMRout, NULL);
      pca_ComputeNoise (TabMRin, TabNoiseIn, TabNoiseOut); 
      pca_Thresold (TabMRout, TabNoiseOut);
   }
 
   else pca_TransfSignal (TabMRin, TabMRout);
}




/*void PCA_MR1D::invtransform (MR_1D *TabMRin , MR_1D *TabMRout) {

  fltarray Vect(_NbSpectre);
  fltarray Vectout(_NbSpectre);
  
  for (int b=0;b<_NbBand;b++) {
    for (int i=0; i<TabMRin[0].size_scale_np(b); i++) {
      
      int k;
      for (k=0; k<_NbSpectre; k++) Vect(k) = (TabMRin[k])(b,i);
      MatCor[b].transform(Vect, Vectout);
      for (k = 0; k<_NbSpectre; k++) (TabMRout[k])(b,i) = Vectout(k);
    }
  } 
  add_mean (TabMRout);
}*/


void  PCA_MR1D::invsubtransform (MR_1D *TabMRin , MR_1D *TabMRout,
                                 int *NbrEigen, 
                                 int TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND]){
 
  fltarray Vect(_NbSpectre);
  fltarray Vectout(_NbSpectre);
  
  // number of eigen used ...
  //-------------------------
  int Ne=0;
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

  for (int b=0;b<_NbBand;b++) {

     int N = NbrEigen[b];
     if ((N < 1) || (N > _NbSpectre)) N = _NbSpectre;
    
     for (int i=0; i<TabMRin[0].size_scale_np(b); i++) {
	   
        for (int k=0; k<N; k++) {
            if ((TabKill != NULL) && (TabKill[k][b] == 1)) Vect(k) = 0.;
            else Vect(k) = TabMRin[k](b,i);
        }	
        MatCor[b].invtransform(Vect, Vectout, N);
        for (int k=0; k<_NbSpectre; k++) TabMRout[k](b,i) = Vectout(k);
     }   
  } 
  
  add_mean (TabMRout);
}






PCA_MR1D::~PCA_MR1D () {

  for (int b=0;b<_NbBand;b++) CorrelMat[b].free();
  free (MatCor);       
  free (CorrelMat);
}

   
void PCA_MR1D::print(int Band) {
    cout << "band " << Band+1 << " :" << endl;
    MatCor[Band].print();
    cout << endl;
}
 

void PCA_MR1D::print() {
  for (int b=0;b<_NbBand;b++) print(b);
}


void PCA_MR1D::pca_TransfSignal (MR_1D *TabMRin , MR_1D *TabMRout,
                                 MR1DNoiseModel *TabNoiseIn) {

  fltarray Vect(_NbSpectre);
  fltarray Vectout(_NbSpectre);     
     
  for (int b=0;b<_NbBand;b++) {
    for (int i=0; i<TabMRin[0].size_scale_np(b); i++) {
 
      int k;
      for (k=0; k<_NbSpectre; k++) Vect(k) = (TabMRin[k])(b,i);
	 /*if (TabNoiseIn != (MRNoiseModel*)NULL && b != NbrBand-1 )
            Vect(k) = (((TabNoiseIn[k])(b,i,j)==true) ?  
                                                 (TabMRin[k])(b,i,j) : 0.0);
         else Vect(k) = (TabMRin[k])(b,i,j);*/
      MatCor[b].transform(Vect, Vectout);
      for (k=0; k<_NbSpectre; k++) (TabMRout[k])(b,i) = Vectout(k);
    }
  } 
} 


void PCA_MR1D::pca_ComputeNoise  (MR_1D *TabMRin,
                                  MR1DNoiseModel *TabNoiseIn,
				  MR1DNoiseModel *TabNoiseOut) {
       
       

  type_noise ae_TypeNoise = TabNoiseIn[0].which_noise();
  	
  if (! one_level_per_pos (ae_TypeNoise)) {

    fltarray Vect(_NbSpectre);
    fltarray Vectout(_NbSpectre);     

    for (int b=0;b<_NbBand-1;b++) {

      int k;
      for (k=0; k<_NbSpectre; k++) Vect(k) = pow((double) TabNoiseIn[k].sigma(b,0),2.);
      //for (int k=0; k<_NbSpectre; k++) cout <<"sigma("<<b<<","<<k<<") in:"<<TabNoiseIn[k].sigma(b)<<endl;
      MatCor[b].transform_noise(Vect, Vectout);
      for (k=0; k<_NbSpectre; k++) TabNoiseOut[k].sigma(b,0) = sqrt(Vectout(k));
      //for (int k=0; k<_NbSpectre; k++) cout <<"sigma("<<b<<","<<k<<") out:"<<TabNoiseIn[k].sigma(b)<<endl;      

    }

  } else { // one sigma at (b,i,j)

    fltarray Vect(_NbSpectre);
    fltarray Vectout(_NbSpectre);     
     
    for (int b=0;b<_NbBand-1;b++) {
      for (int i=0; i<TabMRin[0].size_scale_np(b); i++) {
 
        int k;
        for (k=0; k<_NbSpectre; k++) Vect(k) = pow ((double) TabNoiseIn[k].sigma(b,i),2.);
        MatCor[b].transform_noise(Vect, Vectout);
        for (k=0; k<_NbSpectre; k++) TabNoiseOut[k].sigma(b,i) = sqrt (Vectout(k));
      }
    } 
  }
}


void PCA_MR1D::pca_Thresold (MR_1D *TabMRin_out, MR1DNoiseModel *TabNoiseOut) {

   for (int k=0; k<_NbSpectre; k++) 
      TabNoiseOut[k].threshold(TabMRin_out[k]);

}


float PCA_MR1D::Correl_Without_Normalisation (const fltarray & Sp1, 
                                              const fltarray & Sp2, 
					      float ValMin) {
   double Sum_XY;
 
   // test image size
   if (Sp1.n_elem() != Sp2.n_elem()) {
       
      cerr << "Error in correlation routine: images have different sizes ..." << endl ; 
      cerr << "   spectre 1: " <<  Sp1.n_elem() << endl ;
      cerr << "   spectre 2: " <<  Sp2.n_elem() << endl ;
      exit(-1);
   }
    
   Sum_XY = 0.;
   int Nl=Sp1.n_elem();
   for (int i=0;i<Nl;i++) {
      if (ABS(Sp1 (i)) > ValMin && ABS(Sp2 (i)) > ValMin)
         Sum_XY += Sp1 (i) * Sp2 (i);
   }
     
   return (Sum_XY/Nl);
}

