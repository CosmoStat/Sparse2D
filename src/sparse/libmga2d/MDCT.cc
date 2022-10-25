/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  21/04/02
**
**    File:  MDCT.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale DCT
**    -----------
**
******************************************************************************/

#include "MDCT.h"

/****************************************************************************/
//                  Multiscale DCT
/****************************************************************************/

void MDCT::alloc(int Nl, int Nc, int Nbr_Plan, int FirstBlockSize)
{
    int s;
    NbrScale = Nbr_Plan;
    TabBlockSize.alloc(Nbr_Plan);
    TabBlockSize(0) =  FirstBlockSize;
    for (s=1; s < Nbr_Plan; s++) TabBlockSize(s) = 2*TabBlockSize(s-1);
    TabDCT = new LOCAL_DCT2D [Nbr_Plan];
    // TabTransMDCT = new Ifloat [Nbr_Plan];

    if (Verbose == True) cout << " Multiscale DCT allocation: Nbr_Plan =  " << Nbr_Plan << endl;
    for (s=0; s < Nbr_Plan; s++)
    {
        TabDCT[s].alloc(Nl, Nc,TabBlockSize(s), BlockOverlap);
	// TabTransMDCT [s] = & (TabDCT[s].DCTIma);
	if (Verbose == True)
	{
	   printf("   Band %d, BlockSize = %d, Nl = %d, Nc = %d \n",
	          s+1, TabBlockSize(s), nl(s), nc(s));
	}
    }
    AWT.alloc(TabTransWT, Nl, Nc, Nbr_Plan);
}

/****************************************************************************/

void MDCT::transform(Ifloat &Image)
{
   int s,i,j;
   if (Verbose == True) cout << "Isotropic WT " << endl;
   AWT.transform(Image, TabTransWT, NbrScale);
   for (s = 0; s < NbrScale; s++)
   {
       if (Verbose == True) cout << "  Local DCT, band " << s + 1 << endl;
       TabDCT[s].transform(TabTransWT[s]);
       if (Verbose == True)
       {
          INFO_X(TabDCT[s].DCTIma, "DCT band");
	  cout << "  Norm band  = " << norm_band(s) << endl;
       }
   }
}

/****************************************************************************/

void MDCT::recons(Ifloat &Image)
{
   int s;
   for (s = 0; s < NbrScale; s++)
   {
        if (Verbose == True) cout << "  Local DCT recons, band " << s + 1 << endl;
        TabDCT[s].recons(TabTransWT[s]);
   }
   if (Verbose == True) cout << "Isotropic Wavelet recons " << endl;
   AWT.recons(TabTransWT, Image, NbrScale);
}

/****************************************************************************/

float MDCT::norm_band(int s)
{
   static double TN[10] =
    {0.892954,  0.201288, 0.0856941, 0.0419776, 0.0201896, 0.0106804, 0.00813682,
     0., 0., 0.};
   if ((s < 0) || (s >= 10)) return 0;
   else return TN[s];
}


/****************************************************************************/

void MDCT::threshold(float NoiseLevel, float NSigma)
{
   int s,i,j;

   for (s=0; s < nbr_scale(); s++)
   {
      float Noise = NoiseLevel * norm_band(s);
      cout << "  Local DCT, band " << s + 1 << " Level = " << Noise << endl;

      for (i=0; i < nl(s); i++)
      for (j=0; j < nc(s); j++)
      {
         float Level = NSigma*Noise;
         float Coef = (*this)(s,i,j);
	 if (ABS((TabDCT[s])(i,j)) < Level) (TabDCT[s])(i,j) = 0.;
      }
   }
}
/****************************************************************************/
