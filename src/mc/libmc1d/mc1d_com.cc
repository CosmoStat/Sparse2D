/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.3
**
**    Author: 02/29/00
**
**    Date:  01/02/09
**    
**    File:  mc2d_com.cc
**
*******************************************************************************/

#include "mc1d_com.h"
 
const char* CorCompTransform1d (CorrelComputeType1d CorrelComp) {
   switch (CorrelComp) {
   case E_LOCAL1D  : return ("One correlation matrix per band"); break;
   case E_GLOBAL1D : return ("One correlation matrix for all bands"); break;
   case E_GLOBAL_WITHOUT_LAST1D : 
      return ("One correlation matrix for all band, except the last scale"); break;
   case E_IMPORT1D : return ("Import correlation matrix"); break;
   default : return ("correl type unknown !!!"); break;
   }
};

const char* CorNoiseTransform1d (CorrelNoiseType1d CorrelNoise) {
   switch (CorrelNoise) {
   case E_WITHOUT1D : return ("Compute the correl. matrix with all wavelet coefficients"); break;
   case E_THRESHOLD1D : return ("Compute the correl. matrix only with significant wavelet coefficients"); break;
   case E_LINEAR1D : return ("Compute the correl. matrix with weighted wavelet coefficients"); break;
   //case E_PROBA1D : return ("Compute the correl. matrix with weighted wavelet coefficients"); break;
   default : return ("Correl noise unknown !!!"); break;
   }
};   
			  
			  
void correl_compute1d (CorrelComputeType1d CorrelComp) {
    fprintf(OUTMAN, "         [-x CorrelMat_Method]\n");
    fprintf(OUTMAN, "              0 : %s \n", CorCompTransform1d(E_LOCAL1D));
    fprintf(OUTMAN, "              1 : %s \n", CorCompTransform1d(E_GLOBAL1D));
    fprintf(OUTMAN, "              2 : %s \n", CorCompTransform1d(E_GLOBAL_WITHOUT_LAST1D));   
    fprintf(OUTMAN, "              3 : %s \n", CorCompTransform1d(E_IMPORT1D));    
    fprintf(OUTMAN, "              default is %s.\n",CorCompTransform1d(CorrelComp));  
} 

void correl_noise1d (CorrelNoiseType1d CorrelNoise) {
    fprintf(OUTMAN, "         [-y Noise correlation type]\n");
    fprintf(OUTMAN, "              0 : %s \n", CorNoiseTransform1d(E_WITHOUT1D));
    fprintf(OUTMAN, "              1 : %s \n", CorNoiseTransform1d(E_THRESHOLD1D));
    fprintf(OUTMAN, "              2 : %s \n", CorNoiseTransform1d(E_LINEAR1D));
    //fprintf(OUTMAN, "              3 : %s \n", CorNoiseTransform1d(E_PROBA1D));       
    fprintf(OUTMAN, "              default is %s.\n",CorNoiseTransform1d(CorrelNoise));  
}

/******************************************************************************
classe to_Iter
******************************************************************************/


void to_Iter1d::iter_Begin () {

  //cout << " ==> BEGIN" << endl;

  o_2dResidu.alloc (iter_getpResult()->res_GetpParam()->i_NbPts, 
                    iter_getpResult()->res_GetpParam()->i_NbSpectre);  
  o_2dInit.alloc (iter_getpResult()->res_GetpParam()->i_NbPts, 
                  iter_getpResult()->res_GetpParam()->i_NbSpectre); 
  NbIter = 0;

  for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbSpectre;k++) {
    for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbPts;i++) {
        o_2dInit(i,k) = iter_getpResult()->res_GetpParam()->po_TabSpectre[k](i);   // Einit
	o_2dResidu(i,k) = 0.0;                                            // I(0) = 0
      }
  }
  
  NbIter++;
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_Residu();
}
   
Bool to_Iter1d::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > iter_getpResult()->res_GetpParam()->i_NbIter) { 
   
      for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbSpectre;k++) 
         for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbPts;i++)
            iter_getpResult()->res_GetpParam()->po_TabSpectre[k](i) = o_2dResidu(i,k); 
	            
      return False;
   } else return True;

}

void to_Iter1d::operator++ () {(*this)++;}
void to_Iter1d::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_Residu ();
}

void to_Iter1d::iter_Residu () {

   for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbSpectre;k++)
      for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbPts;i++) {
	 
            // I(i) = I(i) + E(i)
            o_2dResidu(i,k) += iter_getpResult()->res_GetpParam()->po_TabSpectre[k](i);
	                
	    // E(i) = Einit - I(i)
	    iter_getpResult()->res_GetpParam()->po_TabSpectre[k](i) = o_2dInit(i,k) - o_2dResidu(i,k);
	    
   }    
}
 
 
 
 
/******************************************************************************
classe to_SoftIter
******************************************************************************/

 
void to_SoftIter1d::iter_Begin () {

  //cout << " ==> BEGIN" << endl;

  //o_2dCurrentSol.alloc (iter_getpResult()->res_GetpParam()->i_NbPts, 
  //                      iter_getpResult()->res_GetpParam()->i_NbSpectre);  
  NbIter = 0;

  //for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbSpectre;k++)
  //  for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbPts;i++) 
  //	o_2dCurrentSol(i,k) = 0.0;                                            
   
  NbIter++;
  
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_ActionOnSol();
}


Bool to_SoftIter1d::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > iter_getpResult()->res_GetpParam()->i_NbIter) { 
   
      //for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbSpectre;k++) 
      //   for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbPts;i++)
      //      iter_getpResult()->res_GetpParam()->po_TabSpectre[k](i) 
      //        = o_2dCurrentSol(i,k); 
	            
      return False;
   } else return True;
}
 
void to_SoftIter1d::operator++ () {(*this)++;}
void to_SoftIter1d::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_ActionOnSol ();
}

 
 
 
 
 
 
 
Bool control_mr1d (int NbMultRes1d, MR_1D* TabMr1d) {


  // same number of band
  int NbBand = TabMr1d[0].nbr_band();
  for (int k=1; k<NbMultRes1d; k++) {
    if (TabMr1d[k].nbr_band() != NbBand)  {
      cerr << "Error : image have different number of band" << endl;
      exit (-1);
    }
  }

  // read the data, and verify that all images have the same size
  for (int b=0; b<NbBand; b++) {
    for (int k=1; k<NbMultRes1d; k++) {

      if (   TabMr1d[k].size_scale_np(b) != TabMr1d[0].size_scale_np(b)) {
        cerr << "Error : spectre of different size ..." << endl ;
        cerr << "   spectre 1: " <<  TabMr1d[0].size_scale_np(b) << endl ;
        cerr << "   spetcre " << k+1 << ": " << TabMr1d[k].size_scale_np(b) << endl ;
        exit(-1);
      }
    }
  }
  return (True);
}


Bool control_spectre (int NbSpectre, fltarray* TabSpectre) {

  if (NbSpectre == 0) return(False);
  if (TabSpectre[0].ny() != 0) return(False);

  // read the data, and verify that all images have the same size
  for (int k=1; k<NbSpectre; k++) {

    if (   TabSpectre[k].n_elem() != TabSpectre[0].n_elem()) {
      cerr << "Error : spectre of different size ..." << endl ;
      cerr << "   spetcre 1: " <<  TabSpectre[0].n_elem() << endl;
      cerr << "   spectre " << k+1 << ": " << TabSpectre[k].n_elem() << endl;
      exit(-1);
    }
    if (TabSpectre[k].ny() != 0) return(False);
  }
  return (True);
}
