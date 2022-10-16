/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  13.08.99
**    
**    File:  MR_HaarPoissonFilter.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise management
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "FFTN_2D.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_HaarPoisson.h"
#include "MR_HaarPoissonFilter.h"
  
/****************************************************************************/
#define DEBUG_HH 1

// void HFilter::mr_ortho_regul_ima_rec(MultiResol &MR_Step, int s, Ifloat &Result, int Step)
// {
//    int Nl = Result.nl();
//    int Nc = Result.nc();
//    int b,i,j,w,Iter=0;
//    Ifloat Lap(Nl,Nc, "Dat");
//    MultiResol MR_Iter (Nl, Nc, 2, MR_Step.Type_Transform, "MR_Iter");
//    MR_Iter.SBFilter = MR_Step.SBFilter;
//    MR_Iter.TypeNorm =  MR_Step.TypeNorm;
//    int Nw = MR_Iter.nbr_mr_coeff();
//    unsigned char *TabSup = new unsigned char [Nw];
//    float F;
//    float F_old=0; 
//    float Delta_F=10;;
//    float Level;
//    
//    type_border Bord=I_CONT;
//    float *TabSave = new float [Nw];
//    float Delta_w,LevelQ;
//    double LambdaCoef, NormVal = POW2(s+1);
// 
//    w =0;
//    for (b=0; b < MR_Iter.nbr_band()-1; b++)
//    for (i=0; i < MR_Iter.size_band_nl(b); i++)
//    for (j=0; j < MR_Iter.size_band_nc(b); j++)
//    {
//       // Save initial value of dequantized coefficients for range constraint
//       TabSave[w]=MR_Step(b,i,j);
//       // compute the multiresolution support
//       TabSup[w] = (ABS(MR_Step(b,i,j)) > FLOAT_EPSILON) ? 1 : 0; 
//       w++;
//    }
// 	  
//    while (Iter < MaxIterPerScaleRec)
//    // while (((Delta_F*Delta_F) > .001 )&&(Iter<MaxIter))
//    {
//        MR_Step.recons(Result);
//        
//        // compute the Laplacian
//        for (i=0;i<Nl;i++)
//        for (j=0;j<Nc;j++)
// 	   Lap(i,j) =  4*Result(i,j) - (Result(i-Step,j,Bord) +
// 					Result(i,j-Step, Bord) + 
// 					Result(i,j+Step, Bord) +
// 					Result(i+Step,j, Bord));
// 
//        // positivity constraint on reconstructed image
//        threshold(Result);
//        
//        // compute the multiresolution transform of the laplacian 
//        MR_Iter.transform(Lap); 
//        F=0.;   
//        // compute the  new solution
//        w=0;
//        if (TypeBGR == HPF_BGR_FLAT)  Level = (float) get_level(Lambda, s);
//        for (b=0; b < MR_Step.nbr_band()-1; b++)
//        for (i=0; i < MR_Step.size_band_nl(b); i++)
//        for (j=0; j < MR_Step.size_band_nc(b); j++)
//        {
//            switch (TypeBGR)
//            {
//              case HPF_BGR_FLAT: 
//                  // Level = (float) get_level(Lambda, s);
//                break;
//              case HPF_BGR_MR:    
//                  // Get lambda from baground model image
// 	         // Lambda value must be normalized in count per pixel
// 	         //              ==> we multiply it by the  normalization value
//                  LambdaCoef = MR_Step(3,i,j) / NormVal;
//                  Level = get_level(LambdaCoef,s);
//                break;
//               case HPF_BGR_IMA:
//               case HPF_BGR_ITER:
// 	         LambdaCoef = (TabImaModel[s])(i,j);
//                  Level = get_level(LambdaCoef, s);
//               break;
//            }
// 
//  	  //main constraint
// 	  MR_Step(b,i,j) -= 0.2 * MR_Iter(b,i,j);
// 	  F += MR_Iter(b,i,j);
// 
// 	  // range constraint on thresholded-restored coefficients
//           if (TabSup[w] == 0)
//           {
// 	     if (ABS(MR_Step(b,i,j))>Level)   
// 	             MR_Step(b,i,j) = (MR_Step(b,i,j)>0.0)? Level:-Level;
// 	  }
//           else
//           {
//  	      LevelQ =  Level / 4.;
//   	      Delta_w = MR_Step(b,i,j) - TabSave[w];
// 	      if (ABS(Delta_w) > LevelQ) 
// 	       MR_Step(b,i,j) = (Delta_w>0.0)? TabSave[w]+LevelQ:TabSave[w]-LevelQ;
//           }
//           w++;
//        }
//        if (MR_Step.nbr_mr_coeff()!=0)  
// 	 Delta_F = (F-F_old)/ MR_Step.nbr_mr_coeff();
//        else  Delta_F = 0.;
// #if DEBUG_HH
//        cerr << "Iter = " << Iter+1<< " F = " << F << "  Delta_F =" <<  Delta_F*Delta_F <<  endl;
// #endif
//        F_old=F;
//        Iter++;
//    }
//    MR_Step.recons(Result);
// 
// 
//    delete [] TabSup;
//    delete [] TabSave;  
// }
//  
/****************************************************************************/
 

/*****************************************************************************/

void HFilter::laplacien_recons(MultiResol &MR_Step, int s, Ifloat &Result)
{
   float *TabLevel;
   int SizeTabLevel;
   double LambdaCoef, NormVal = POW2(s+1);
   int i,j,b,w,Nw = MR_Step.nbr_mr_coeff();

   if (TypeBGR == HPF_BGR_FLAT)
   {
        TabLevel = new float [1];
        TabLevel[0] = (float) get_level(Lambda, s);
        SizeTabLevel = 1;
   }
   else
   {
        TabLevel = new float [Nw];
        SizeTabLevel = Nw;
        w = 0;  
        for (b=0; b < MR_Step.nbr_band()-1; b++)
        for (i=0; i < MR_Step.size_band_nl(b); i++)
        for (j=0; j < MR_Step.size_band_nc(b); j++)
        {
            if (TypeBGR == HPF_BGR_MR) LambdaCoef = MR_Step(3,i,j) / NormVal;
            else LambdaCoef = (TabImaModel[s])(i,j);
            TabLevel[w++] = get_level(LambdaCoef, s);
        }
    }
   mr_ortho_regul_ima_rec(MR_Step,  Result, TabLevel,
                          MaxIterPerScaleRec,SizeTabLevel);
   delete [] TabLevel; 
}

/****************************************************************************/

void HFilter::threshold_scale(MultiResol &MR_Data, Ifloat & LowResol, int s)
{
   int i,j,b;
   float Level;
   int Cpt=0;
   double LambdaCoef, NormVal = POW2(s+1);

   // Search the Lambda parameter per scale  
   if (TypeBGR == HPF_BGR_FLAT)
   {
      Level = (float) get_level(Lambda, s);
      if (Verbose == True)  
      {
 	  cout << "Scale " << s +1 << "  Level = " << Level << endl; //  << " Unnorm. Level = " << Level * Norm << endl;
	// cout << " coef = " << NormCoef1 << endl;
      }
   }

  for (b=0; b <= 2; b++)
  {
    int BB = 3*s+b;
    for (i=0; i < MR_Data.size_band_nl(BB); i++)
    for (j=0; j < MR_Data.size_band_nc(BB); j++)
    {
       float HaarCoeff = MR_Data(BB,i,j);
       switch (TypeBGR)
       {
          case HPF_BGR_FLAT: 
               // Level = (float) get_level(Lambda, s);
               break;
          case HPF_BGR_MR:    
               // Get lambda from baground model image
	       // Lambda value must be normalized in count per pixel
	       //              ==> we multiply it by the  normalization value
               LambdaCoef =  LowResol(i,j) / NormVal;
               Level = get_level(LambdaCoef,s);
               break;
          case HPF_BGR_IMA:
          case HPF_BGR_ITER:
	      LambdaCoef = (TabImaModel[s])(i,j);
              Level = get_level(LambdaCoef, s);
              break;
        }

       if ((s < FirstScale) || 
              ((KeepPositivSup == True) && (HaarCoeff < 0))) HaarCoeff  = 0.;  
       else
       {
	   if (SoftThreshold == False)
	   {
	        if (ABS(HaarCoeff) < Level) HaarCoeff  = 0.;
           }
	   else
	   {
	        if (HaarCoeff>0) HaarCoeff = MAX(0,HaarCoeff-Level);
		else HaarCoeff = MIN(0,HaarCoeff+Level);
           }
        }
	if (HaarCoeff != 0) Cpt ++;
 	MR_Data(BB,i,j) = HaarCoeff;
   }
   NbDetect += Cpt;
  }
   //if (Verbose == True)
   //   cout << "Scale "<< s+1 << " Number of detect. coef = " << Cpt << endl;
}

/****************************************************************************/
 
void HFilter::threshold_all_scale(MultiResol &MR_Data)
{  
   int i,j,b,s;
   float Level;
   double LambdaCoef;
    // Search the Lambda parameter per scale  


   for (s = NbrPlan-2; s >= 0; s--)
   {
      // Variable used for threshold calculation
      NormCoef1 = 1. / (double) POW2(s+2);
      NormCoef2 = (double) POW2(2*(s+1));

      if (TypeBGR == HPF_BGR_FLAT)
      {
         Level = (float) get_level(Lambda, s);
         if (Verbose == True)  
         {
 	    cout << "Scale " << s+1 << " Level = " << Level << endl; //  << " Unnorm. Level = " << Level * Norm << endl;
	   // cout << " coef = " << NormCoef1 << endl;
         }
      }

      for (b=0; b <= 2; b++)
      {
         int Band = 3*s+b;
         for (i=0; i < MR_Data.size_band_nl(Band); i++)
         for (j=0; j < MR_Data.size_band_nc(Band); j++)
         {
            float HaarCoeff = MR_Data(Band,i,j);
            switch (TypeBGR)
            {
               case HPF_BGR_FLAT: 
               // Level = (float) get_level(Lambda, s);
               break;
               case HPF_BGR_MR:    
                 cerr << "Error: HPF_BGR_MR cannot be used in this method ... " << endl;
                 exit(-1);
               break;
              case HPF_BGR_IMA:
              case HPF_BGR_ITER:
	        LambdaCoef = (TabImaModel[s])(i,j);
                Level = get_level(LambdaCoef, s);
              break;
            }
 
           if ((s < FirstScale) || 
              ((KeepPositivSup == True) && (HaarCoeff < 0))) HaarCoeff  = 0.;  
           else
           {
	     if (SoftThreshold == False)
	          HaarCoeff = hard_threshold(HaarCoeff, Level);
 	     else HaarCoeff = soft_threshold(HaarCoeff, Level);
           }
  	  MR_Data(Band,i,j) = HaarCoeff;
       }
    }
  }
}

/****************************************************************************/


double HFilter::get_level(float LambdaVal, int Scale)
// return the threshold level for a given scale and a given lambda
// and for EpsilonPoisson
{
    double Level;

    if ((TypeThreshold == HPF_LAMBDA_K_1) || (TypeThreshold == HPF_LAMBDA_K_2)) 
    {
        double z,z2,z4,Lambdaj = NormCoef2*LambdaVal;
	
	z = (TypeThreshold == HPF_LAMBDA_K_1) ? N_Sigma: sqrt(2.*log((double) Np));
	z2 = z*z;
        z4 = z2*z2;
        Level = NormCoef1*(z2+sqrt(z4+4.*Lambdaj*z2));
    }
    else  
    {
       int Norm = haar_norm_coef(3*Scale, NbrBand);
       double L = LambdaVal * POW2(2*(Scale+1));
       
       Level = get_harr_poisson_threshold(L, EpsilonPoisson) / Norm;
    }
    return Level;
}

/****************************************************************************/

// void HFilter::get_scale_level(Ifloat &LowResol, int s, float *TabLevel)
// // get the threshold level using the low resolution scale as a background 
// // model: MR_Step.band(3) = background model
// // Norm = normalization coefficient when L2 normalization is used to
// //        get an unirmalized coefficient
// // s = scale 
// // TabLevel : out = level table 
// {
//     int w=0,i,j;
//     double NormVal = POW2(s+1);
//     for (i=0; i < LowResol.nl(); i++)
//     for (j=0; j < LowResol._nc(); j++)
//     {
//         // Get lambda from baground model image
// 	// Lambda value must be normalized in count per pixel
// 	//              ==> we multiply it by the  normalization value
//        double LambdaCoef = MR_Step(3,i,j) / NormVal;
//        TabLevel[w] = get_level(LambdaCoef,s);
//        w++;
//     }
// }     
// 
// /****************************************************************************/
// 
// void HFilter::get_model_scale_level(int Scale, Ifloat & ImaModel, MultiResol &MR_Step, float *TabLevel)
// // get the threshold level using an image model given by ImaModel
// // ImaModel = model image
// // MR_Step: in = 2 scales (4 bands) Haar transform at a given resolution 
// {
//     int w=0,b,i,j;
//      
//     for (b=0; b < MR_Step.nbr_band()-1; b++)
//     for (i=0; i < MR_Step.size_band_nl(b); i++)
//     for (j=0; j < MR_Step.size_band_nc(b); j++)
//     {
//         // get the threshold value
// 	double LambdaCoef = ImaModel(i,j);
//        TabLevel[w] = get_level(LambdaCoef, Scale);
//        w++;
//     }
// }

/****************************************************************************/

void HFilter::get_mr_ima(Ifloat &Ima, int NbrScale)
// set TabImaModel, the background model at each scale, from the
// the background image Ima
// Ima: in = background image model
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int Nls = (DecimatedTransform == True) ? (Nl+1) / 2: Nl;
   int Ncs = (DecimatedTransform == True) ? (Nc+1) / 2: Nc;
   int s,i,j;
   int Step = 1;
   for (s = 0; s < NbrScale-1; s++)
   {
      TabImaModel[s].alloc(Nls,Ncs,"ima model");
      for (i=0; i < Nls; i++)
      for (j=0; j < Ncs; j++)
      {
         if (DecimatedTransform == True) 
         {
	    if (s == 0)
 	      TabImaModel[s](i,j) = 0.25*(Ima(2*i,2*j,I_MIRROR) + Ima(2*i,2*j+1,I_MIRROR)
		                     + Ima(2*i+1,2*j,I_MIRROR) + Ima(2*i+1,2*j+1,I_MIRROR));
 	    else
 	    TabImaModel[s](i,j) =   0.25*(TabImaModel[s-1](2*i,2*j,I_MIRROR)   + TabImaModel[s-1](2*i,2*j+1,I_MIRROR)
		                     + TabImaModel[s-1](2*i+1,2*j,I_MIRROR) + TabImaModel[s-1](2*i+1,2*j+1,I_MIRROR));
         }
         else
         {
	    if (s == 0)
 	      TabImaModel[s](i,j) = 0.25*(Ima(i,j,I_MIRROR) + Ima(i,j+Step,I_MIRROR)
		                     + Ima(i+Step,j,I_MIRROR) + Ima(i+Step,j+Step,I_MIRROR));
 	    else
 	    TabImaModel[s](i,j) =   0.25*(TabImaModel[s-1](i,j,I_MIRROR)   + TabImaModel[s-1](i,j+Step,I_MIRROR)
		                     + TabImaModel[s-1](i+Step,j,I_MIRROR) + TabImaModel[s-1](i+Step,j+Step,I_MIRROR));
         }
      }
      if (DecimatedTransform == True) 
      {
         Nls = (Nls+1) / 2;
         Ncs = (Ncs+1) / 2;
      }
      else Step *= 2;
   }
}
/****************************************************************************/

void HFilter::alloc(type_background_model TBGR, 
                    type_haar_poisson_threshold TThreshold, 
		    MultiResol &MR_Data, char *Name_BGR_In)
{
    FirstScale=0;
    KeepPositivSup=False;
    SoftThreshold=False;
    TabImaModel = NULL;
    MaxFilterIter=4;
    DecimatedTransform = (MR_Data.Type_Transform == TO_MALLAT) ? True: False; 
    if (DecimatedTransform == True)  MaxIterPerScaleRec=6;
    else MaxIterPerScaleRec=1;
    TypeBGR=TBGR;
    TypeThreshold=TThreshold;
    NbrPlan = MR_Data.nbr_scale();
    NbrBand = MR_Data.nbr_band();
    Nlima = MR_Data.size_ima_nl();
    Ncima = MR_Data.size_ima_nc();
    Np = Nlima*Ncima;
    
    // For optimal PRESS filtering, we need the same information,
    // i.e. image at different resolution level
    // as with the automatic background estimation from a coarser scale
    if (TypeThreshold==HPF_PRESS)
    {
        TypeBGR = HPF_BGR_IMA;
	TabImaModel = new Ifloat [NbrPlan];
    }
    else
    {
       if (TypeBGR == HPF_BGR_ITER) BGR_Ima.alloc(Nlima,  Ncima, "BGR Ima"); 
       if (TypeBGR == HPF_BGR_IMA)
       {
          io_read_ima_float(Name_BGR_In, BGR_Ima);
          if ((BGR_Ima.nl() != Nlima) || (BGR_Ima.nc() != Ncima))
          {
             cerr << "Background image must have the same size as the input image ... " << endl;
	     exit(-1);
          }
       }
    
       if ((TypeBGR == HPF_BGR_IMA) || (TypeBGR == HPF_BGR_ITER) || 
         (TypeThreshold==HPF_PRESS))
       {
          TabImaModel = new Ifloat [NbrPlan];
          if (TypeBGR == HPF_BGR_IMA) get_mr_ima(BGR_Ima, NbrPlan);
       }
    }    
}

/****************************************************************************/

void HFilter::filter_noiter (MultiResol &MR_Data, Ifloat &Result)
{
   Ifloat Imag;
   int Nb = NbrBand;
   int s,Nls,Ncs;
   int Step=1;
   SubBandFilter WT1D(MR_Data.SBFilter, MR_Data.TypeNorm);

   MultiResol MR_Step;
   s = NbrPlan-2;
   Nls = (DecimatedTransform == True) ? MR_Data.size_scale_nl(s, I_SMOOTH): Nlima;
   Ncs = (DecimatedTransform == True) ? MR_Data.size_scale_nc(s, I_SMOOTH): Ncima; 
   // cout << "Nls = " << Nls << " Ncs = " << Ncs << endl;
   if (DecimatedTransform == True) Result.resize(Nls,Ncs);
   if ((Nls != MR_Data.band(Nb-1).nl()) || (Ncs != MR_Data.band(Nb-1).nc()))
   {
       cerr << "Error: band size are not the same in filter_noiter ... " << endl;
       cerr  << "   Nls = " << Nls << " Ncs = " << Ncs << endl;
       exit(-1);
   }   
   Result = MR_Data.band(Nb-1);
   // io_write_ima_float("xx_last.fits", Result);

   // Filtering scale by scale
   NbDetect = 0;
   if (DecimatedTransform == False)
                  for (s = 0; s < NbrPlan-2; s++) Step *= 2;

   for (s = NbrPlan-2; s >= 0; s--)
   {
      // Variable used for threshold calculation
      NormCoef1 = 1. / (double) POW2(s+2);
      NormCoef2 = (double) POW2(2*(s+1));

      threshold_scale(MR_Data, Result, s);

      // Reconstruct the solution
      if (DecimatedTransform == True) 
      {
         // Create a MR object with two scales
         Nls = MR_Data.size_scale_nl(s-1, I_SMOOTH);
         Ncs = MR_Data.size_scale_nc(s-1, I_SMOOTH);
         MR_Step.alloc(Nls, Ncs, 2, MR_Data.Type_Transform, "MR_Step");   
         MR_Step.SBFilter = MR_Data.SBFilter;
         MR_Step.TypeNorm =  MR_Data.TypeNorm; 
         MR_Step.band(0) = MR_Data.band(3*s);
         MR_Step.band(1) = MR_Data.band(3*s+1);
         MR_Step.band(2) = MR_Data.band(3*s+2);
         MR_Step.band(3) = Result;
         Result.resize(Nls,Ncs);
         if (MaxIterPerScaleRec > 1) laplacien_recons(MR_Step, s, Result);
         // mr_ortho_regul_ima_rec(MR_Step, s, Result);
         else MR_Step.recons(Result);
         MR_Step.free();
      }
      else 
      {
	 PAVE_2D_WT P2WT(WT1D);
         Ifloat *TabBand = MR_Data.tabband();
         P2WT.recons2d(TabBand[3*s], TabBand[3*s+1], TabBand[3*s+2],
	               Result, Result, Step);
         Step /= 2;
      }
      threshold(Result);
      // INFO(Result, "result");
   }
   // if (Verbose == True)
   //            cout <<  " Number of detect. coef = " <<  NbDetect << endl;
 
}

/****************************************************************************/

void HFilter::press_filter (Ifloat & Data, MultiResol &MR_Data, Ifloat &Result)
{
    int b,i,j;
    get_mr_ima(Data, NbrPlan);
    for (b = 0; b < MR_Data.nbr_band()-1; b++)
    {
       int s = b / 3;
       for (i=0; i < MR_Data.size_band_nl(b); i++)
       for (j=0; j < MR_Data.size_band_nc(b); j++)
       {
           float Coef = MR_Data(b,i,j);
	   float Coef2 = Coef*Coef;
	   float Sigma = TabImaModel[s](i,j);
           float Hpress = (Coef2 > FLOAT_EPSILON) ? (Coef2 - Sigma) / Coef2: 0.;
	   MR_Data(b,i,j) *= Hpress;
       }
    }
    MR_Data.recons(Result);
    threshold(Result);
}

/****************************************************************************/

void HFilter::filter (Ifloat & Data, MultiResol &MR_Data, Ifloat &Result)
{    
   if (TypeThreshold == HPF_PRESS) press_filter(Data, MR_Data, Result);
   else 
   {
      if (TypeBGR != HPF_BGR_ITER) filter_noiter (MR_Data, Result);
      else 
      {
        int Iter,i,j;
        Ifloat Resi;
        float Mean, Delta;
      
        if (Verbose == True) Resi.alloc(Nlima,Ncima,"resi");
        Mean = average(Data);
        BGR_Ima.init(Mean);
        for (Iter=0; Iter < MaxFilterIter; Iter++)
        {
          get_mr_ima(BGR_Ima, NbrPlan);
          filter_noiter(MR_Data, Result);
          if (Verbose == True) 
	  {
             for (i=0; i < Nlima; i++)
             for (j=0; j < Ncima; j++) Resi(i,j) = BGR_Ima(i,j) - Result(i,j);
 	     Delta = sigma(Resi);
             cout << "Filter. Iter " << Iter+1 << " Delta = " << Delta << endl << endl;
          }
	  BGR_Ima = Result;
        }
     }
   }   
}

/****************************************************************************/

void HFilter::test(MultiResol &MR_Data)
{
    int b,Nb = NbrBand; 
    fltarray Tab_Proba(Nb);
    // Tab_Proba.init(0.00269);
    for (b=0; b < Nb-1; b++) Tab_Proba(b) = EpsilonPoisson; // 0.00269;
    fltarray ThresholdMax(Nb);
    //StatPoissonHaar StatHaar(Lambda,  MR_Data.nbr_band()); 
    //StatHaar.write();
    //fltarray ThresholdMin(Nb);
    // StatHaar.find_threshold(MR_Data.nbr_band(), Tab_Proba, ThresholdMin, ThresholdMax);
    cout << "TEST routine ... " << endl;
    cout << "   Lambda = " << Lambda << endl;
    cout << "   EpsilonPoisson = " << EpsilonPoisson << endl;
    for (b=0; b < MR_Data.nbr_band()-1; b+=3) 
    { 
       int s = b / 3;
       double L = Lambda * POW2(2*s+1);
       cout << "nb = " << Nb << "   L = " << L << " Eps = " <<  Tab_Proba(b) << endl;
       ThresholdMax(b) = (float) get_harr_poisson_threshold(L, (double) Tab_Proba(b));
       double NCoef1 = (double) haar_norm_coef(b, NbrBand);
       //for (int i=0; i < MR_Data.size_band_nl(b); i++)  
       //for (int j=0; j < MR_Data.size_band_nc(b); j++)  MR_Data(b,i,j) *= NCoef1;
       ThresholdMax(b) /= NCoef1;
       cout << "Band " << b+1 << " Threshold = " << ThresholdMax(b)*NCoef1 << " Norm. Threshold = " << ThresholdMax(b) << endl;
       cout << "   Lambda = " << L << " NormCoef = " << NCoef1 << endl;
       cout << "   Kol = " << get_level(Lambda, s) << endl;;
    } 
    //MR_Data.write("xx_h");
}

