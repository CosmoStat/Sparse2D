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
**    Date:  1.03.99
**    
**    File:  MR_CorrNoise.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise management
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "RadNoise.h"
#include "NR.h"
#include "IM_IO.h"

/****************************************************************************/

RadNoiseMap::RadNoiseMap()
{
    Verbose = False;
    Histo.alloc(SIZE_CRAD_HISTO);
    Repart.alloc(SIZE_CRAD_HISTO);
    Tab_Bin.alloc(SIZE_CRAD_HISTO);
    Sigma = 0.;
    StepHisto=MinHisto=MaxHisto=0;
} 

/****************************************************************************/
 
void RadNoiseMap::set(Ifloat &NoiseIma)
{
  int i,j,k;
   
  if (NoiseIma.n_elem() == 0)
  {
     cout << "Error in RadNoiseMap::alloc: incorrect NoiseIma parameter ... " << endl;
     exit(-1);
  }
  int Nlb =  NoiseIma.nl();
  int Ncb =  NoiseIma.nc();
  MinHisto=min(NoiseIma);
  MaxHisto=max(NoiseIma);
 
  if (Verbose == True) cout << "Ima :  Histogram [min,max] = [" << MinHisto << "," << MaxHisto << "]" << endl;
      
  // Compute histogram bin size  
  StepHisto = (MaxHisto - MinHisto)/(float)(SIZE_CRAD_HISTO-1);

  // initialization
  for (k=0;k<SIZE_CRAD_HISTO;k++)
  {
      Tab_Bin(k) = MinHisto + (float) k * StepHisto;
      Histo(k)=0.;
      Repart(k)= 0.;
  }
      
  // compute histogram 
  for(i=0; i< Nlb;i++)
  for(j=0; j< Ncb;j++)
  { 
     k = (int)  ((NoiseIma(i,j)-MinHisto)/StepHisto) ;
     if ((k < 0) || (k >= SIZE_CRAD_HISTO))
     {
	      cout << "Error: k = " << k << " in band " <<  endl;
              cout << "min = " << MinHisto << endl;
	      cout << "max = " << MaxHisto << endl;
	      cout << "step = " << StepHisto << endl;
	      cout << "val = " <<  NoiseIma(i,j) << endl;
	      exit(-1);
      }
      Histo(k)+=1; 
  }
 
  // histogram normalisation and repartition function
  for (k=0; k < SIZE_CRAD_HISTO; k++)
  { 
     Histo(k)/=(Nlb*Ncb);
     if (k==0) Repart(k)=Histo(k);
      else Repart(k) = Repart(k-1)+Histo(k);
  }
}

/****************************************************************************/

float RadNoiseMap::prob(float Val)
{
  int k = (int)  ( (Val-Tab_Bin(0)) / (Tab_Bin(1)-Tab_Bin(0)) );
  if ((k>=SIZE_CRAD_HISTO)|| (k< 0))  return 0.;
  else return Histo(k);
}

/****************************************************************************/

float RadNoiseMap::repartition(float Val)
{
  int k = (int)  ( (Val-Tab_Bin(0)) / (Tab_Bin(1)-Tab_Bin(0)) );
  if (k>=SIZE_CRAD_HISTO)  return 1.;
  else if (k<0)  return 0.;
  else return Repart(k);
}

/****************************************************************************/ 

void RadNoiseMap::find_threshold(float Proba, 
                       float &ThresholdMin, float &ThresholdMax)
{
   int k=0;

   while (Repart(k) < Proba) k++;
   ThresholdMin = Tab_Bin(k);
   while (Repart(k) < (1.-Proba)) k++;
   ThresholdMax = Tab_Bin(k);
}

/****************************************************************************/ 

void RadNoiseMap::find_gthreshold(float N_Sigma, 
                       float &ThresholdMin, float &ThresholdMax)
{
   float Proba = (1. - erff((double) N_Sigma / sqrt((double) 2.)));;
   find_threshold(Proba, ThresholdMin, ThresholdMax);
}

/****************************************************************************/

void RadNoiseMap::write()
{
  fits_write_fltarr((char*)"Tab_Histo", Histo);
  fits_write_fltarr((char*)"Tab_Repart", Repart);
  fits_write_fltarr((char*)"Tab_Bin", Tab_Bin);
}

/****************************************************************************/

