/******************************************************************************
**                   Copyright (C) 1997 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_BudgetComp.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image compression/decompression using wavelet
**    -----------  
**       
*******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Noise.h"
// #include "MR_Filter.h"
#include "IM_CompTool.h"
#include "MR_Sigma.h"
#include "MR_Calloc.h"
#include "IM_Comp.h"
#include "MR_Comp.h"

 
/*************************************************************************/

static float getentrop(int *qdata, int NbrCoef, int Minq, int Maxq)
{
   int i;
   int Nh = Maxq - Minq + 1;
   int *Histo = new int [Nh]; 
   float TotEntr=0.;
   double Log2 = log(2.);
   
   for (i=0; i < Nh; i++) Histo[i]=0;
   for (i=0; i < NbrCoef; i++)
   {
       int Ind = qdata[i] - Minq;
       if ((Ind < 0) || (Ind >= Nh))
       {
          cerr << "Error in getentrop: Ind = " << Ind << " Nh = " << Nh << endl;
          cerr << "    Minq = " << Minq << " Maxq  = " << Maxq << endl;
          cerr << "    qdata = " << qdata[i] << " i = " << i << endl;
          exit(-1);
       }
       Histo[Ind] ++;
   }
   // cout << "end histo " << endl;
   for (i=0; i < Nh; i++) 
   {
      double Prob = Histo[i] / (double) NbrCoef;
      if (Histo[i] != 0) TotEntr -= Prob * log(Prob)/Log2;
   }
   delete [] Histo;
   return TotEntr;
}

/*************************************************************************/

static void getratedist(MultiResol &MR_Data,BandRateDist & BRD,float *TabSigma)
{
   int NbrBand = BRD.NbrBand;
   int NbrQuant = BRD.NbrQuant;
   int i,j,k;
   float QuantifPar;
   float Min,Max,Mean,Coef,Decod,Val;
   
   // cout << "GRD " << endl;
   
   for (i = 0; i < NbrBand; i++)
   {
      int NbrCoef = MR_Data.nbr_coeff_in_band(i);
      int *qdata = i_alloc(NbrCoef);
      Min = Max = MR_Data(i,0);
      Mean=0.;
      for (j = 1; j < NbrCoef; j++) 
      {
          Coef = MR_Data(i,j);
          if (Min >  Coef) Min =  Coef;
          else if (Max <  Coef) Max = Coef;
          Mean +=  Coef;
      }
      if ( ABS(Min) > Max) Max = ABS(Min);
      
      if (i ==  NbrBand-1) Mean /= (float) NbrCoef;
      else Mean = 0;
      BRD.TabMin(i) = Min;
      BRD.TabMax(i) = Max;
      BRD.TabMean(i) = Mean;
      
      // k = 0, we do not keep the data
      BRD.Rate(i,0) = 0;
      BRD.Dist(i,0) = 0;
      BRD.QuantCoeff(i,0) = 0;
      for (j = 0; j < NbrCoef; j++) 
               BRD.Dist(i,0) += (MR_Data(i,j)-Mean)*(MR_Data(i,j)-Mean);
      BRD.Dist(i,0) /= (float) NbrCoef;
      
      for (k = 1; k < NbrQuant; k++) 
      {
        int Minq,Maxq;
         
        int  Prec = (1 << k)-1;
        QuantifPar =  Max / (float) Prec * TabSigma[i];
        BRD.QuantCoeff(i,k) = QuantifPar;
        Minq = quant(Min, QuantifPar);
        Maxq = quant(Max, QuantifPar);
        BRD.Dist(i,k)=0.;
        for (j = 0; j < NbrCoef; j++)
        {
             qdata[j] = quant( MR_Data(i,j), QuantifPar);
             Decod = qdata[j] * QuantifPar;
             // Val = (MR_Data(i,j) - Decod) / TabSigma[i];
             Val = (MR_Data(i,j) - Decod);
             BRD.Dist(i,k) += Val*Val;
        }
        BRD.Dist(i,k) /= (float) NbrCoef;
        // cout << "getentrop" << endl;
        BRD.Rate(i,k) = getentrop(qdata, NbrCoef, Minq, Maxq) * NbrCoef;
        // cout << "Rate " << BRD.Rate(i,k) << endl;
      }
      i_free(qdata);
   }
} 

/*************************************************************************/

void mr_getquant(MultiResol & MR_Data,float *TabSigma, int budget, 
               float *TabQuantif, float *TabMean, float *TabRate,
	       Bool Verbose)
{  
   int i;
   int Budget = budget;
   int nSets = MR_Data.nbr_band();
   float *weight = new float[nSets];  
                  // perceptual weights for coefficient sets
   for (i = 0; i < nSets; i++) weight[i] = 1.0;
   if (MR_Data.Type_Transform == TO_MALLAT) 
      for (i = 0; i < nSets; i++) 
      {
	  if (i==2) weight[i] /= 2;
	  else if (i % 3 == 2) weight[i] /= 1.5;
	  if (i < 3) weight[i] /= 1.5;	    
      }
   
   int nQuant =  NBR_QUANT;          
   
   // if (Verbose == True) cout << "nset = " << nSets << endl;
   // for (i = 0; i < nSets; i++) TabSigma[i] = 0.25;
   
  BandRateDist BRD(nSets, nQuant);
  getratedist(MR_Data, BRD, TabSigma);
  
   /*
   for (i = 0; i < nSets; i++)
   for (int j = 0; j <  nQuant; j++)
   {
      cout << "Set " << i+1 << " Quant " << j+1 << " Rate = " 
         << BRD.Rate(i,j) << " Dist = " << BRD.Dist(i,j) << 
            " Quantif = " << BRD.QuantCoeff(i,j) << endl;
   }
   */
   
   Allocator *allocator = new Allocator ();
  // Use rate/distortion information for each subband to find bit
  //    allocation that minimizes total (weighted) distortion subject
  //    to a byte budget
  Budget -= nSets * 12;  // subtract off approximate size of header info
  allocator->optimalAllocate (BRD, nSets, Budget, True, weight);
  if (Verbose == True) printf ("Target rate = %d bytes\n", Budget);
    
  for (i = nSets-1; i >= 0; i--)
  {      
     int prec  = allocator->precision[i];
     TabQuantif[i] = BRD.QuantCoeff(i,prec);
     TabMean[i] = BRD.TabMean(i);
     if (TabRate != NULL) TabRate[i] = BRD.Rate(i,prec);
  } 
  
  delete [] weight;
}

/*********************************************************************/


static void mr_iter_comp(Ifloat &Imag, MultiResol &MR_Data, float *TabQuantif, 
                  int MaxIter, MRNoiseModel *ModelData)
{
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int i,j,NbrCoef,Iter;
   Ifloat ImagRec(Nl,Nc,"ImagRec");
   Ifloat Result(Nl,Nc,"Result");
   ImagRec = Imag;
   Result.init();
   int NbrIterRec=3;
   float Quantif;
   
   for (Iter = 0; Iter < MaxIter; Iter ++)
   {  
      if (ModelData != NULL) ModelData->threshold(MR_Data);
      
      for (i = MR_Data.nbr_band()-1; i >= 0; i--)
      {      
          NbrCoef = MR_Data.nbr_coeff_in_band(i);
          int qdata; 
          Quantif =  TabQuantif[i];
         if (Quantif > 0) 
         {
            for (j = 0; j < NbrCoef; j++)
            {
               qdata  = quant( MR_Data(i,j), Quantif);
               MR_Data(i,j) = (float) qdata * Quantif;
            }
         }
         else for (j = 0; j < NbrCoef; j++) MR_Data(i,j) = 0.;
      }
      mr_ortho_regul_rec (MR_Data,  ImagRec,  TabQuantif, NbrIterRec);
      Result += ImagRec;
      ImagRec = Imag - Result;
      MR_Data.transform(ImagRec);
   }
   MR_Data.transform(Result);
   if (ModelData != NULL) ModelData->threshold(MR_Data);
}

/*********************************************************************/

float get_dwt_quant_param(int s)
{
   float Q = 7. / (float) DEFAULT_SIGNAL_QUANTIF;
   if (s==2) Q *= 2;
   else if (s % 3 == 2) Q *= 1.5;
   if (s < 3) Q *= 1.5;
   return Q;
}    

/*********************************************************************/

void mr_budcompress(Ifloat &Imag, type_transform Transform, 
                    int Nbr_Plan, int Bud, FILE *FileDes, Ifloat &ImagRec,
                    Bool KillNoise, float &Noise, float NSigma, 
                    float SignalQuant, Bool SupIsol, Bool SupDil, Bool SupNeg,
                    int MedianWindowSize, int MaxIter, Bool Verbose, Bool RecImag)
{ 
    int i,j,s;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    // MultiResol MR_Data (Nl, Nc, Nbr_Plan, Transform, "MR_Imag");
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    sb_type_norm Norm = NORM_L1;
    int NbrUndec = -1;
    type_sb_filter SB_Filter = F_MALLAT_7_9;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MultiResol MR_Data;
    MR_Data.alloc (Nl, Nc, Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);

    MR_Data.ExactPyrRec = True;
    MR_Data.MedianWindowSize = MedianWindowSize;
    int Budget = Bud;
    noise_compute (MR_Data);
    float *TabSigma = new float [MR_Data.nbr_band()];
    float *TabParamQuant = new float [MR_Data.nbr_band()];
    float *TabRate = new float [MR_Data.nbr_band()];
    
    // initialize TabSigma and TabParamQuant
    for (i=0; i < MR_Data.nbr_band(); i++)
    {
        //float Diff;
        //TabParamQuant[i] =  SignalQuant;
        //Diff = SignalQuant - DEFAULT_SIGNAL_QUANTIF;
        //if (Diff > FLOAT_EPSILON)
        //        TabParamQuant[i] *= (1. +  Diff /  POW2(s));
	if (Transform != TO_MALLAT)
          TabParamQuant[i] = get_quant_param(i, MR_Data.nbr_band(), 
                                          SignalQuant, DEFAULT_SIGNAL_QUANTIF);
	else TabParamQuant[i] = SignalQuant*get_dwt_quant_param(i);
        TabSigma[i]= MR_Data.band_norm(i);
        if (TabSigma[i] < FLOAT_EPSILON) TabSigma[i] = 0.05;   
    }  
    
    // estimate the quantification parameters
    float *TabQuantif = new float [MR_Data.nbr_band()];
    float *TabMean = new float [MR_Data.nbr_band()];
    Budget -= 12;  // Nl,Nc, and Nbr_Plan space
                   // but header is not taken into account
 
  // cout << "Budget = " << Budget << endl;                  
    if (KillNoise == False)
    {
       MR_Data.transform(Imag);   
       if (Budget > 1) mr_getquant(MR_Data, TabSigma, Budget, TabQuantif, TabMean, TabRate);
       else 
       for (i=0; i < MR_Data.nbr_band(); i++)
       { 
          TabQuantif[i] =  TabParamQuant[i]*TabSigma[i];
          TabMean[i] = 0.;
       }
       if ((MR_Data.Set_Transform == TRANSF_MALLAT) && (MaxIter > 1))
            mr_iter_comp(Imag, MR_Data, TabQuantif, MaxIter, NULL);
           
// for (i=0; i < MR_Data.nbr_band(); i++) 
//   cout << i+1 << " Tq = " << TabQuantif[i]/TabSigma[i] << " Par = " << TabParamQuant[i] << " Sig = " << TabSigma[i] << endl;
    
    }
    else
    {
   
       // if (verbose) cout << "Data thresholding at : " << NSigma << endl;
       // MRNoiseModel ModelData(NOISE_GAUSSIAN, Nl,Nc, Nbr_Plan, Transform);    
       MRNoiseModel ModelData;
       ModelData.alloc(NOISE_GAUSSIAN,  Nl,Nc, 
                    Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);

       if (Noise > FLOAT_EPSILON) ModelData.SigmaNoise = Noise;
        for (s=0; s < MR_Data.nbr_band(); s++) ModelData.NSigma[s]=(s==0)?NSigma+1:NSigma;
       ModelData.OnlyPositivDetect = SupNeg;
       ModelData.SigmaDetectionMethod = SIGMA_BSPLINE;
       ModelData.model(Imag, MR_Data);
// cout << "Nl = " << Nl << " Nc = " << Nc << " Noise = " << ModelData.SigmaNoise << " NSigma " << NSigma << endl;
//       INFO( Imag, "data ");

//       INFO(MR_Data.band(0), "wavelet scale 1 ");
//       INFO(MR_Data.band(4), "wavelet scale 4 ");

        if (SupIsol == True) ModelData.kill_isol(0);
       Noise = ModelData.SigmaNoise;
       if (SupDil == True)   
                  for (s = 0; s <  Nbr_Plan-1; s++) ModelData.dilate_support(s);
                        
       if (Budget > 1)
       {
           if (Noise < 1.) 
                    for (i=0; i < MR_Data.nbr_band(); i++) TabSigma[i] *= Noise;
           mr_getquant(MR_Data, TabSigma, Budget, TabQuantif, TabMean, TabRate);
       }
       else 
       {
          // cout << "TabQuantif " << endl;
          for (i=0; i < MR_Data.nbr_band(); i++)
          { 
          TabQuantif[i] = Noise*TabParamQuant[i]*TabSigma[i];
          TabMean[i] = 0.;
          }
          TabQuantif[MR_Data.nbr_band()-1] = TabQuantif[MR_Data.nbr_band()-2];
//    for (i=0; i < MR_Data.nbr_band(); i++) 
//                  cout << i+1 << " Tq = " << TabQuantif[i] << " Par = " << TabParamQuant[i] << " Sig = " << TabSigma[i] << endl;
       }
       if ((MaxIter > 1) && (MR_Data.Set_Transform == TRANSF_MALLAT))
                   mr_iter_comp(Imag, MR_Data, TabQuantif, MaxIter, &ModelData);
       else ModelData.threshold(MR_Data);
     }
    // INFO(MR_Data.band(0), "wavelet scale 1 ");
    // INFO(MR_Data.band(4), "wavelet scale 4 ");
    
    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 
    // if (Verbose == True) cout << "encode" << endl;
     
    // encode the data
    // extern long size_enc;
    int NbrCoef,Nls,Ncs;
    float Quantif;
    int *qdata;
    for (i = MR_Data.nbr_band()-1; i >= 0; i--)
    {      
       // INFO(MR_Data.band(i), "wavelet scale");
       NbrCoef = MR_Data.nbr_coeff_in_band(i);

       qdata = i_alloc(NbrCoef); 
       Quantif =  TabQuantif[i];
       Nls = MR_Data.size_band_nl (i);
       Ncs = MR_Data.size_band_nc (i);
       if (Quantif > 0) 
       {
         for (j = 0; j < NbrCoef; j++)
         {
           qdata[j] = quant( MR_Data(i,j), Quantif);
           MR_Data(i,j) = (float) qdata[j] * Quantif;
         }
         encode(FileDes,  qdata,  Nls,  Ncs, Quantif);
       }
       else 
       { 
          if (i !=  MR_Data.nbr_band()-1) qdata[0] = 0;
          else qdata[0] =  quant(TabMean[i], 1.);
          encode( FileDes,  qdata, 1,1, (float) 0.);
         //for (j = 0; j < NbrCoef; j++) qdata[j] = (int) qmean;
         //encode( FileDes,  qdata,   MR_Data.size_band_nl (i), 
         //       MR_Data.size_band_nc(i), 1.);
          for (j = 0; j <  NbrCoef; j++)  MR_Data(i,j) = qdata[0];
       }
       if (Verbose == True)
       {
           //extern long size_enc;
           //cerr << "Set " << i+1 << "  NbrCoef = " <<  NbrCoef 
           //     << " Quant = " << Quantif <<  " Size = " << size_enc << endl;
           //if (Budget > 1)
           //  cout << "  Entropy = " << TabRate[i]/8.  << " SizeEnc = " << size_enc << endl;
        }
        i_free(qdata);
     }
    // MR_Data.write("xx_rec.mr");
   if (RecImag == True) MR_Data.recons(ImagRec);
   
   delete [] TabSigma;
   delete [] TabParamQuant;
   delete [] TabRate;
   delete [] TabMean;
   delete [] TabQuantif;

}
 
/*************************************
 
       DECOMPRESSION PART

**************************************/

void mr_buddecompress (Ifloat &Imag, type_transform Transform, 
                       int Resol, FILE *FileDes, int &Nli, int &Nci, int NbrIterRec,
		       Bool Verbose)
{
    int i,j;
    int Nl,Nc;
    int Nbr_Plan;
 
    /* read the original image size */
    Nl =readint(FileDes); /* x size of image */
    Nc =readint(FileDes); /* y size of image */
    Nli = Nl;
    Nci = Nc;
    Imag.resize (Nl, Nc);

    /* read the number of scales */
    Nbr_Plan=readint(FileDes); 
     
    if (Verbose == True)
    {
       fprintf (stderr,"DECODING wavelet transform: nbr_scale = %d\n", Nbr_Plan);
       fprintf (stderr,"Original image: Nl = %d, Nc = %d\n", Nl, Nc);
     }

    /* Resol must be < Nbr_Plan and > 0 */
    if ((Resol < 0) || (Resol >= Nbr_Plan))
    {
      fprintf (stderr, "Error: bad parameter Resol: Resol = %d\n", Resol);
      fprintf (stderr, "       Resol must verify 0 <= Resol < %d\n", Nbr_Plan);
      exit (0);
    }

    /* if we don't want the full resolution */
    if (Resol > 0)
    {
       for (int s=0; s < Resol; s++)
       {
          Nl = (Nl+1)/2;
          Nc = (Nc+1)/2;
       }
       Nbr_Plan -= Resol;
    }     
    if (Verbose == True) printf ("\n Create an image of size (%d, %d)\n", Nl, Nc);
    Imag.resize (Nl, Nc);
    MultiResol MR_Data (Nl, Nc, Nbr_Plan, Transform, "MR_Imag");
    Ifloat Buffer(Nl, Nc, "Buffer mr_decomp");
    extern long size_enc;
    int NbrBand = MR_Data.nbr_band();
    float *LevelQ = new float[NbrBand];
    for (i =  NbrBand -1 ; i  >= 0; i--) 
    {
       int NbrCoef;
       NbrCoef = MR_Data.nbr_coeff_in_band(i);
       int *data;
       float Level;
       int Nls,Ncs;
       decode (FileDes, &data, &Nls, &Ncs, &Level);
       LevelQ[i] = Level;
       int SizeScale = size_enc;
       if (Verbose == True)
             cerr << "Set " << i+1 << " " << Nls << "x" << Ncs << 
                 " Quantif = " << Level  <<  " Size = " << SizeScale << endl;
       // details which_detail;
//        int s;
//        MR_Data.band_to_scale(i, s, which_detail);
//        if (verbose)
//        {
//           cerr << "Scale  " << s+1;
//           if (which_detail == D_DIAGONAL) cerr << " D=diagonal " << endl;
//           else if (which_detail == D_VERTICAL) cerr << " D=vertical " << endl;
//           else if (which_detail == D_HORIZONTAL) cerr << " D=horizontal " << endl;
//           else if (which_detail == I_SMOOTH) cerr << " D=smooth " << endl;
//           else cerr << " D=-1 " << endl;
//           cerr << "sizescale nc = " << MR_Data.size_scale_nc(s,which_detail) << endl;
//           cerr << "sizeband nc = " << MR_Data.size_band_nc(i) << endl;
//       }
       
       
       if (Level > FLOAT_EPSILON)
       {
          for (j = 0; j <  NbrCoef; j++)
          {
             if (NbrBand > 1) MR_Data(i,j) = ((float) data[j]) * Level;
             else Imag(j) = ((float) data[j]) * Level;
          }
       }
       else
       {
          if (NbrBand > 1)
                for (j = 0; j <  NbrCoef; j++)  MR_Data(i,j) =  data[0];
          else  for (j = 0; j <  NbrCoef; j++)  Imag(j) = data[0];
       }
       i_free(data);
    }        

    if (NbrBand > 1)
    {       
       if ((NbrIterRec > 1) && (MR_Data.Set_Transform == TRANSF_MALLAT))
          mr_ortho_regul_rec (MR_Data, Imag, LevelQ, NbrIterRec);
       else MR_Data.recons (Imag);
    }
    delete [] LevelQ;
}
    
     

/****************************************************************************/

 
