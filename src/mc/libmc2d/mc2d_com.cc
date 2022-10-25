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
#include "mc2d_com.h"








/******************************************************************************
classe to_Iter
******************************************************************************/


void to_Iter::iter_Begin () {

  //cout << " ==> BEGIN" << endl;

  o_3dResidu.alloc (iter_getpResult()->res_GetpParam()->i_NbLin, 
                    iter_getpResult()->res_GetpParam()->i_NbCol, 
                    iter_getpResult()->res_GetpParam()->i_NbImage);  
  o_3dInit.alloc (iter_getpResult()->res_GetpParam()->i_NbLin, 
                  iter_getpResult()->res_GetpParam()->i_NbCol, 
                  iter_getpResult()->res_GetpParam()->i_NbImage); 
  NbIter = 0;

  for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbImage;k++) {
    for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbLin;i++)
      for (int j=0;j<iter_getpResult()->res_GetpParam()->i_NbCol;j++) {
        o_3dInit(i,j,k) = iter_getpResult()->res_GetpParam()->po_TabIma[k](i,j);   // Einit
	o_3dResidu(i,j,k) = 0.0;                                            // I(0) = 0
      }
  }
  
  NbIter++;
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_Residu();
}
   
Bool to_Iter::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > iter_getpResult()->res_GetpParam()->i_NbIter) { 
   
      for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbImage;k++) 
         for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbLin;i++)
            for (int j=0;j<iter_getpResult()->res_GetpParam()->i_NbCol;j++)
               iter_getpResult()->res_GetpParam()->po_TabIma[k](i,j) = o_3dResidu(i,j,k); 
	            
      return False;
   } else return True;

}

void to_Iter::operator++ () {(*this)++;}
void to_Iter::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_Residu ();
}

void to_Iter::iter_Residu () {

   for (int k=0;k<iter_getpResult()->res_GetpParam()->i_NbImage;k++)
      for (int i=0;i<iter_getpResult()->res_GetpParam()->i_NbLin;i++)
         for (int j=0;j<iter_getpResult()->res_GetpParam()->i_NbCol;j++) {
	 
            // I(i) = I(i) + E(i)
            o_3dResidu(i,j,k) += iter_getpResult()->res_GetpParam()->po_TabIma[k](i,j);
	                
	    // E(i) = Einit - I(i)
	    iter_getpResult()->res_GetpParam()->po_TabIma[k](i,j) = o_3dInit(i,j,k) - o_3dResidu(i,j,k);
	    
   }    
}
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
Bool control_mr2d (int NbMultRes2d, MultiResol* TabMr2d) {


  // same number of band
  int NbBand = TabMr2d[0].nbr_band();
  for (int k=1; k<NbMultRes2d; k++) {
    if (TabMr2d[k].nbr_band() != NbBand)  {
      cerr << "Error : image have different number of band" << endl;
      exit (-1);
    }
  }

  // read the data, and verify that all images have the same size
  for (int b=0; b<NbBand; b++) {
    for (int k=1; k<NbMultRes2d; k++) {

      if (   TabMr2d[k].size_band_nl(b) != TabMr2d[0].size_band_nl(b) 
          || TabMr2d[k].size_band_nc(b) != TabMr2d[0].size_band_nc(b)) {
        cerr << "Error : image of different size ..." << endl ;
        cerr << "   image 1: " <<  TabMr2d[0].size_band_nl(b) ;
        cerr << "X"  <<  TabMr2d[0].size_band_nc(b) << endl ;
        cerr << "   image " << k+1 << ": " << TabMr2d[k].size_band_nl(b) ;
        cerr << "X"  <<  TabMr2d[k].size_band_nc(b) << endl ;
        exit(-1);
      }
    }
  }
  return (True);
}


Bool control_image (int NbImage, Ifloat* TabImage) {

  if (NbImage == 0) return(False);

  // read the data, and verify that all images have the same size
  for (int k=1; k<NbImage; k++) {

    if (   TabImage[k].nl() != TabImage[0].nl() 
        || TabImage[k].nc() != TabImage[0].nc()) {
      cerr << "Error : image of different size ..." << endl ;
      cerr << "   image 1: " <<  TabImage[0].nl();
      cerr << "X"  <<  TabImage[0].nc() << endl;
      cerr << "   image " << k+1 << ": " << TabImage[k].nl();
      cerr << "X"  <<  TabImage[k].nc() << endl;
      exit(-1);
    }
  }
  return (True);
}
