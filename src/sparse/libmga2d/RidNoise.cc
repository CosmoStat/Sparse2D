/******************************************************************************
**                   Copyright (C) 200 by CEA + Stanford university
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  06/04/00
**    
**    File:  RidNoise.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for ridgelet noise management
**    ----------- 
**                 
******************************************************************************/

#include "RidNoise.h"
#include "NR.h"

/****************************************************************************/

void StatRidNoise::init(Ifloat &RidTrans, int Nscale, Ridgelet &Rid)
{
  int i,j,k,s;
  fltarray ThresholdMin, ThresholdMax ;
  fltarray Buff;
 
  if (RidTrans.n_elem() == 0)
  {
     cout << "Error in StatRidNoise::alloc: incorrect RidTrans parameter ... " << endl;
     exit(-1);
  }
  if (Nscale != NbrBand)
  {
     NbrBand = Nscale;
     StepHisto.alloc(NbrBand-1);  
     MinHisto.alloc(NbrBand-1);  
     MaxHisto.alloc(NbrBand-1);       
     TabSigma.alloc(NbrBand-1);
     Histo.alloc(NbrBand-1, RID_SIZE_HISTO);
     Repart.alloc(NbrBand-1, RID_SIZE_HISTO);
     Tab_Bin.alloc(NbrBand-1, RID_SIZE_HISTO);
  }
 
  for (s=0; s < NbrBand-1; s++)
  {
     int NFirst = Rid.rid_pos(s);
     int NLast = Rid.rid_pos(s) + Rid.rid_size(s);
     int Nl = RidTrans.nl();
     int Nelem = Nl*Rid.rid_size(s);
     int Ind=0;
     Buff.reform(Nelem);
     for (i=0; i < Nl; i++)
     for (j=NFirst; j < NLast; j++) Buff(Ind++) =  RidTrans(i,j);
 
      // Compute min max
      MinHisto(s) = Buff.min();
      MaxHisto(s)= Buff.max();
      TabSigma(s)= Buff.sigma();
      if (Verbose == True) cout << "band " << s << ":  Histogram [min,max] = [" << MinHisto(s) << "," << MaxHisto(s) << "]" << endl;
      
      // Compute histogram bin size  
      StepHisto(s) = (MaxHisto(s) - MinHisto(s))/(float)(RID_SIZE_HISTO-1);

      // initialization
      for (k=0;k<RID_SIZE_HISTO;k++)
      {
	  Tab_Bin(s,k) = MinHisto(s) + (float) k * StepHisto(s);
	  Histo(s,k)=0.;
	  Repart(s,k)= 0.;
      }
 
      // compute histogram 
      for(i=0; i< Nelem;i++)
      { 
	  k = (int)  (( Buff(i)-MinHisto(s))/StepHisto(s)) ;
	  if ((k < 0) || (k >= RID_SIZE_HISTO))
	  {
	      cout << "Error: k = " << k << " in band " << s << endl;
              cout << "min = " << MinHisto(s) << endl;
	      cout << "max = " << MaxHisto(s) << endl;
	      cout << "step = " << StepHisto(s) << endl;
	      exit(-1);
	  }
	  Histo(s,k)+=1; 
      }

      // histogram normalisation and repartition function
      for (k=0; k < RID_SIZE_HISTO; k++)
      { 
	  Histo(s,k)/=(Nelem);
	  if (k==0) Repart(s,k)=Histo(s,k);
	  else Repart(s,k) = Repart(s,k-1)+Histo(s,k);
      }
   }
}

/****************************************************************************/

float StatRidNoise::prob(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if ((k>=RID_SIZE_HISTO)|| (k< 0))  return 0.;
  else return Histo(Band,k);
}

/****************************************************************************/

float StatRidNoise::repartition(int Band, float Val)
{
  int k = (int)  ( (Val-Tab_Bin(Band,0)) / (Tab_Bin(Band,1)-Tab_Bin(Band,0)) );
  if (k>=RID_SIZE_HISTO)  return 1.;
  else if (k<0)  return 0.;
  else return Repart(Band,k);
}

/****************************************************************************/ 

void StatRidNoise::find_threshold(int NbrMaxBand, fltarray &Tab_Proba, 
                       fltarray &ThresholdMin, fltarray &ThresholdMax)
{
  int k;

  for (int s = 0; s < NbrMaxBand-1; s++)
  {
    k=0;
    while (Repart(s,k) < Tab_Proba(s)) k++;
    ThresholdMin(s) = Tab_Bin(s,k);
    while (Repart(s,k) < (1.-Tab_Proba(s))) k++;
    ThresholdMax(s) = Tab_Bin(s,k);
  }
}

/****************************************************************************/ 

void StatRidNoise::find_nsigma_threshold(int NbrMaxBand, fltarray &Tab_NSigma, fltarray &ThresholdMin, fltarray &ThresholdMax)
{
  int k;

  for (int s = 0; s < NbrMaxBand-1; s++)
  {          
    float Proba = (1. - erff((double) Tab_NSigma(s) / sqrt((double) 2.)));;
    k=0;
    while (Repart(s,k) < Proba) k++;
    ThresholdMin(s) = Tab_Bin(s,k);
    while (Repart(s,k) < (1.- Proba)) k++;
    ThresholdMax(s) = Tab_Bin(s,k);
  }
}

/****************************************************************************/

void StatRidNoise::write(char *Name)
{
  char FileName[255];

  sprintf(FileName, "%s_histo", Name);
  fits_write_fltarr(FileName, Histo);
  sprintf(FileName, "%s_repart", Name);
  fits_write_fltarr(FileName, Repart);
  sprintf(FileName, "%s_bin", Name);
  fits_write_fltarr(FileName, Tab_Bin);
}

/****************************************************************************/

