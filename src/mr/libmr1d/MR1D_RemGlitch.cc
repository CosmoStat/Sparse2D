/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 14:46:37
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  MR1D_RemGlitch.cc
**
**
****************************************************************************/

// static char sccsid[] = "@(#)MR1D_RemGlitch.cc 3.1 96/05/02 CEA 1996 @(#)";
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// #include "IM1D_IO.h"
#include "GlobalInc.h"
#include "MR1D_Obj.h"
#include "MR1D_Sigma.h"
// #include "MR1D_Filter.h"
#include "MR1D_RemGlitch.h"

float TabCoefPaveMedian[8] = {0.967920, 0.296151, 0.212961, 0.197451, 
                      0.171384, 0.144828, 0.120119, 0.130501};
float TabCoefPyrMedian[8] = {0.918146, 0.373589, 0.260029, 0.203192, 
                      0.135948, 0.112718, 0.087961, 0.157093};

/************************************************************************/

void mr1d_set_support(MR_1D &MR_Data, MR_1D& Sup1d, float Noise_Ima, float N_Sigma, float *TabCoef)
{
    int s;
    float Level;
    int i;
    int Nbr_Plan = MR_Data.nbr_scale()-1;

    for (s = 0; s < Nbr_Plan; s++)
    {
       Level = Noise_Ima * N_Sigma * TabCoef[s];
       for (i = 0; i < MR_Data.size_scale_np(s); i++)
       {
          if (ABS (MR_Data(s,i)) >= Level) Sup1d(s,i) = 1.;
          else Sup1d(s,i) = 0.;
       }
     }
}

/************************************************************************/

void mr1d_threshold_support(MR_1D &MR_Data, MR_1D& Sup1d, float Noise_Ima, float N_Sigma,  float *TabCoef)
{
    int s,i;
    int Nbr_Plan = MR_Data.nbr_scale()-1;
    float Level;

    for (s = 0; s < Nbr_Plan; s++)
    {
       Level = Noise_Ima * N_Sigma * TabCoef[s];
       for (i = 0; i < MR_Data.size_scale_np(s); i++)
          if (Sup1d (s,i) < FLOAT_EPSILON) 
          {
             if (ABS (MR_Data(s,i)) >= Level) Sup1d(s,i) = 1.;
             else   MR_Data(s,i) = 0.;
          }
    }
}

/************************************************************************/

float noise_tmp_cube(fltarray & Data, int Nit)
{
    int It, x,y,z;
    float Mean,Sigma=0.,S0,S1,S2,Sm=0,Val;
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    float c1=-1./sqrt(6.);
    float c2=2./sqrt(6.);

    Mean = 0.;
    for (It = 0; It < Nit; It++)
    {
       S0 = S1 = S2 = 0.;
       for (x=0; x < Nx; x++)
       for (y=0; y < Ny; y++)
       for (z = 1; z < Nz-1; z++)
       {
           Val = c2*Data(x,y,z) + c1*(Data(x,y,z+1)+Data(x,y,z-1));
	   if ((It == 0) || (ABS(Val - Mean) < Sm))
	   { 
	       S0 ++;
	       S1 += Val;
	       S2 += Val*Val;
	   }
       }
       Mean = S1 / S0;
       Sigma = sqrt(S2/S0- Mean * Mean);
       Sm = 3. * Sigma;
    }
    return Sigma;
}

/************************************************************************/

void mr1d_threshold_glitch(MR_1D &MR_Data, fltarray &Mask, float Sigma_Noise, float N_Sigma, int Last_Scale)
{
    int s,i;
    float Level, Level1;

    for (s = Last_Scale-2; s >= 0; s--)
    {
       Level = Sigma_Noise * N_Sigma * TabCoefPaveMedian[s];
       Level1 = Sigma_Noise * TabCoefPaveMedian[s] * 1.;
       for (i = 0; i < MR_Data.size_scale_np(s); i++)
       {
           if (ABS(MR_Data(s,i)) > Level) 
           {
               MR_Data(s,i) = 0.;
               Mask(i) = 1.;
           }
           else if ((ABS(MR_Data(s,i)) > Level1) && (Mask(i) == 1))
           {
               MR_Data(s,i) = 0.;
           }
           else if ((i > 0) && (ABS(MR_Data(s,i)) > Level1) 
                                && (Mask(i-1) == 1))
           {
               MR_Data(s,i) = 0.;
               Mask(i) = 1.;
           }
        }
    }
// cout << "Threshold glitch : " << "Sigma_Noise = " << Sigma_Noise << endl;
// cout << "N_Sigma = " << N_Sigma << endl;
// cout << "Last_Scale = " << Last_Scale << endl;
// cout << "FIND = " << total(Mask) << endl;

}

/************************************************************************/

void mr1d_deglitch(fltarray &Data, fltarray &Ima_SigmaNoise, fltarray &Ima_MeanClip, Bool Deglitch, Bool Filter, int NbrScale, int NscaleFilter, float N_Sigma, int MaxIter_Deglitch, int MaxIter_Filter, float Epsilon)
{
   int i, Nx, Ny, Nz, x, y, z;
   float Sigma;
   int TotalMask=0;

   // new declarations and initializations
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   fltarray Vector(Nz);
   fltarray Mask(Nz);
   fltarray FilterVector(Nz);
   fltarray Resi(Nz);

   MR_1D MR_Data_Pave (Nz, TM1_PAVE_MEDIAN, "MR_Deglitch", NbrScale);
   MR_1D MR_Data_Pyr (Nz, TM1_PYR_MEDIAN, "MR_Filtering", NscaleFilter);
   MR_1D Sup1D (Nz, TM1_PYR_MEDIAN, "MR_Filtering", NscaleFilter);

   for (x=0; x < Nx; x++)
   for (y=0; y < Ny; y++)
   {
       Mask.init();
       for (z=0; z < Nz; z++) Vector(z) = Data(x,y,z);
       Ima_SigmaNoise(x,y) = detect1d_noise_from_med(Vector);
       Sigma =  Ima_SigmaNoise(x,y);

       // Data Filtering
       if (Deglitch ==  True) 
       {
          for (i =0; i < MaxIter_Deglitch; i ++)
          {
             MR_Data_Pave.transform (Vector);
             mr1d_threshold_glitch (MR_Data_Pave, Mask, Sigma,
                                    N_Sigma, MR_Data_Pave.nbr_scale());
             MR_Data_Pave.recons (Vector);
          }
          TotalMask += (int) Mask.total();
       }

       // Data Filtering 
       if (Filter == True)  
       {
          float Delta, Old_Sigma, Sigmav, Norm;
          int Iter = 0;

          MR_Data_Pyr.transform (Vector);
          mr1d_set_support(MR_Data_Pyr, Sup1D, Sigma,N_Sigma, TabCoefPyrMedian);
          mr1d_threshold_support(MR_Data_Pyr, Sup1D,
                                 Sigma, N_Sigma, TabCoefPyrMedian);
          MR_Data_Pyr.recons (FilterVector);
          Norm = Vector.sigma();
          Sigmav = Norm;
          while (True)
          {
             Old_Sigma = Sigmav; 
             Resi = Vector - FilterVector;
             Sigmav = Resi.sigma();
             MR_Data_Pyr.transform (Resi);
             mr1d_threshold_support(MR_Data_Pyr, Sup1D, Sigma,
                                    N_Sigma, TabCoefPyrMedian);
             MR_Data_Pyr.recons (Resi);
             Delta = (Old_Sigma - Sigmav)  / Norm;
             FilterVector += Resi;
             Iter ++;
             if ((Iter >= MaxIter_Filter) || ((Delta < Epsilon) && (Iter > 2))) break;
          }
       }

       if (Filter == True) for (z=0; z < Nz; z++) Vector(z) = FilterVector(z);
       for (z=0; z < Nz; z++) Data(x,y,z) = Vector(z);
       Vector.sigma_clip (Ima_MeanClip(x,y), Sigma);
   }
   cout << "Number of masked values: " << TotalMask << endl;
}

