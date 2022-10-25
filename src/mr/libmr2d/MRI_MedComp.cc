/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.5
**
**    Author: Jean-Luc Starck
**
**    Date:  95/08/28
**    
**    File:  mr_med_comp.c
**
*******************************************************************************
**
**    DESCRIPTION routines which allows to code a pyramidal multiresolution
**    ----------- transform and to decode it.
**
*******************************************************************************
** 
** int mr_level_signif (int s, float Noise, float Nsigma)
**
**  returns Noise * Nsigma * TabSignificantLevel[s]
**
*******************************************************************************
**
** void mri_scaling (Median, Nbr_Plan, Tab_Level)
** 
** reduces the dynamic of range on the multiresolution transform
**
*******************************************************************************
**
** void mr_comp_trans  (int *Imag, int **Median, int *Residual, int Nl, 
**                      int Nc, int Nbr_Plan, float **Noise, 
**                      float Nsigma, int MaxIter, float *Tab_Level)
**
** computes the multiresolution coefficients will allow to represent
** an image with the minimum of coefficient
** Median: OUT = multiresolution transform
** Residual: OUT = Error
** Imag: IN= input image
** Nl,Nc: IN = size of the image
** Nbr_Plan: IN = number of scales
** Noise: OUT = estimated noise
** Nsigma: IN = threshold at Nsigma*Noise at each scale
** MaxIter: IN = number of iterations
** Tab_Level: IN = significant level per scale
**
*******************************************************************************
**
** void encode_med_pyr (int *Imag, int Nl, int Nc, int Nbr_Plan, 
**                      float *Noise, int NSigma, 
**                      int MaxIter, FILE *FileDes, int *Residual)
**
** ecrit sur un fichier la compression de l'image, et l'erreur est
** mis dans residual
**    . transformee multiresolution
**    . quantification iterative
**    . quadtree sur chaque echelle
**    . code de Huffmann sur chaque echelle
**    . l'erreur (image - image decompressee) est mise dans Residual
**
*******************************************************************************
**
** void decode_med_pyr (int **Imag, int *Nl, int *Nc, int *Nbr_Plan, 
**                      int Resol, FILE *FileDes)
**
** read a coded pyramid from file, and decodes it.
**
** Imag: OUT = decoded image
** Nl, Nc: OUT = image size
** Nbr_Plan: OUT = number of scales
** Resol: IN = resolution we want.
**             if (Resol = 0, we want the image at full resolution
** FileDes: IN = file where we read the coded pyramid
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_med_comp.c 1.5 95/08/28 CEA 1994 @(#)";

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Noise.h"
// #include "MR_Filter.h"
#include "IM_CompTool.h"
#include "MR_Sigma.h"
#include "IM_Comp.h"
#include "MR_Comp.h"


int Tab_Nl[MAX_SCALE], Tab_Col[MAX_SCALE], Tab_Pos[MAX_SCALE];

#define DEBUG 0

/****************************************************************************/

static void mri_threshold (int *Median, int Nbr_Plan, float *Tab_Level, float *Tab_Quant, Bool SupIsol, Bool SupNeg)
{
    int s,i,j,ind;
    int *Plan;
    int Nl,Nc,Nls,Ncs;
    float Level;

    // Thresholding and support dilation
    for (s = 0; s < Nbr_Plan-1; s++)
    {
        Plan = Median + Tab_Pos[s];
        Level = Tab_Level[s];
        Nls = Tab_Nl[s];
        Ncs = Tab_Col[s];
        ind=0;
        for (i = 0; i < Nls; i++)
        for (j = 0; j < Ncs; j++)
        {
           if ((s == 0) || (i == 0) || (j == 0) || (i == Nls-1) || (j==Ncs-1))
           {
              if (ABS(Plan[ind]) < Level) Plan[ind] = 0;
           }
           else
           {
                if (        (ABS(Plan[ind]) < Level)  
                        && ( ABS(Plan[(i-1)*Ncs+j]) < Level)
                        && (ABS(Plan[(i-1)*Ncs+j-1]) < Level)
                        && (ABS(Plan[(i-1)*Ncs+j+1]) < Level)
                        && (ABS(Plan[ind-1]) < Level)
                        && (ABS(Plan[ind+1]) < Level)
                        && (ABS(Plan[(i+1)*Ncs+j-1]) < Level)
                        && (ABS(Plan[(i+1)*Ncs+j])   < Level)
                        && (ABS(Plan[(i+1)*Ncs+j+1]) < Level))  Plan[ind] = 0;
           }
           ind++;
        }
    }

   /* kill isolated pixels in the first scale */
   if (SupIsol)
   {
      s = 0;
      Nl = Tab_Nl[s];
      Nc = Tab_Col[s];
      Plan = Median + Tab_Pos[s];

      for (i =1; i < Nl-1; i++)
      for (j =1; j < Nc-1; j++)
      if  (Plan[i*Nc+j] != 0)
      {
              if (    !Plan[(i-1)*Nc+j] && !Plan[(i-1)*Nc+j-1] 
                        && !Plan[(i-1)*Nc+j+1] && !Plan[i*Nc+j-1] 
                        && Plan[i*Nc+j+1] && !Plan[(i+1)*Nc+j-1]
                        && !Plan[(i+1)*Nc+j] && !Plan[(i+1)*Nc+j+1])
                    Plan[i*Nc+j] = 0;
      }
   }


    /* quantization */
    for (s = 0; s < Nbr_Plan; s++)
    {
        Plan = Median + Tab_Pos[s];
        for (i =0; i < Tab_Nl[s]*Tab_Col[s]; i++)
        {
           if (Plan[i] > FLOAT_EPSILON) 
                          Plan[i] = (int)((float) Plan[i] / Tab_Quant[s] + 0.5);
           else if ((SupNeg == False) && (Plan[i] < - FLOAT_EPSILON)) 
                          Plan[i] = (int)((float) Plan[i] / Tab_Quant[s] - 0.5);
                 else Plan[i] = 0;
        }
    }
}

/****************************************************************************/
/*
static void mri_scaling (int *Median, int Nbr_Plan, float *Tab_Level, float *Tab_Quant, int *Support, Bool SupIsol)
{
    int s,i,j,Nl,Nc;
    int *Plan, *PlanSup;

    for (s = 0; s < Nbr_Plan-1; s++)
    {
        Plan = Median + Tab_Pos[s];
        PlanSup = Support + Tab_Pos[s];

        for (i =0; i < Tab_Nl[s]*Tab_Col[s]; i++)
        {
           if ((Plan[i] > Tab_Level[s]) || (PlanSup[i] != 0))
                  PlanSup[i] += (int) ((float) Plan[i] / Tab_Quant[s] + 0.5);
        }
    }

    s = Nbr_Plan-1;
    Plan = Median + Tab_Pos[s];
    PlanSup = Support + Tab_Pos[s];
    for (i =0; i < Tab_Nl[s]*Tab_Col[s]; i++) 
                     PlanSup[i] += (int) ((float)Plan[i] / Tab_Quant[s] + 0.5);


   // kill isolated pixels in the first scale 
   if (SupIsol)
   {
      s = 0;
      Nl = Tab_Nl[s];
      Nc = Tab_Col[s];
      Plan = Support + Tab_Pos[s];

      if (Tab_Quant[s] < Tab_Level[s])
      {
          for (i =1; i < Nl-1; i++)
          for (j =1; j < Nc-1; j++)
          if  (Plan[i*Nc+j] != 0)
          {
              if (    !Plan[(i-1)*Nc+j] && !Plan[(i-1)*Nc+j-1] 
                        && !Plan[(i-1)*Nc+j+1] && !Plan[i*Nc+j-1] 
                        && Plan[i*Nc+j+1] && !Plan[(i+1)*Nc+j-1]
                        && !Plan[(i+1)*Nc+j] && !Plan[(i+1)*Nc+j+1])
                    Plan[i*Nc+j] = 0;
          }
      }
   }

}
 */
/****************************************************************************/

static void mri_comp_trans  (int *Imag, int **Median, int *Result, 
                     int Nl, int Nc, int Nbr_Plan, 
                     float *Noise, float Nsigma, float SignalQuantif,
                     float *Tab_Quant, Bool SupIsol, 
                     Bool NoiseInData, Bool Rec, int MedianWindowSize, Bool SupNeg)
{
    int  s,  Size;
    float Tab_Level[MAX_SCALE], NoiseDat;
    Size = mri_size_medpyr (Nl, Nc, Nbr_Plan);
 
    /* computes the indice position */
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);
    /* computes the multiresolution transform */
    mri_pyrmedian_transform (Imag, Nl, Nc, Median, Nbr_Plan, MedianWindowSize);

    /* noise estimation in data */
    if (NoiseInData == True)
    {
       if (*Noise < FLOAT_EPSILON)
             *Noise = sigma_clip_int (*Median, Nl, Nc, 0) / 0.970291;
       NoiseDat = *Noise;
    }
    else NoiseDat=1.;
    
    /* search the signification level at each scale */
    for (s = 0; s < Nbr_Plan-1; s++)
    {
        if (s == 0) Tab_Level[s] =  NoiseDat * (Nsigma+1) * mr_tab_noise(s);
        else Tab_Level[s] =  NoiseDat * Nsigma * mr_tab_noise(s);
        if (NoiseInData == False) Tab_Level[s] = 0.;
        Tab_Quant[s] = mr_tab_noise(s) * NoiseDat *  
             get_quant_param(s, Nbr_Plan,SignalQuantif, DEFAULT_SIGNAL_QUANTIF);
        if (Tab_Quant[s] < 1.) Tab_Quant[s] = 1.;
    }
    Tab_Quant[Nbr_Plan-1] =  Tab_Quant[Nbr_Plan-2];
 
    for (s = Nbr_Plan; s < MAX_SCALE; s++) Tab_Level[s] = 1;

    /* decrease the dynamic */
    mri_threshold  (*Median, Nbr_Plan, Tab_Level, Tab_Quant, SupIsol, SupNeg);

    /* increase the dynamic */
    if (Rec == True) mri_pyrmedian_rec (*Median, Result, Nl, Nc, Nbr_Plan, Tab_Quant);
}

/****************************************************************************/

void mri_compress (int *Imag, int Nl, int Nc, int Nbr_Plan, float *Noise,
                    float SignalQuantif, float  NSigma,
                    FILE *FileDes, int *ImagRec, Bool SupIsol,
                    Bool NoiseInData, Bool Rec, int MedianWindowSize, 
		    Bool SupNeg, Bool Verbose)
{
    int *Median=NULL;
    float *Tab_Quant;
    int Tab_Nl[MAX_SCALE], Tab_Col[MAX_SCALE], Tab_Pos[MAX_SCALE];
    int s, *Plan;
    float Sigmaf;
    
    Tab_Quant = new float [MAX_SCALE];
    /* calcules index positions */
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);

    /* calcules  the multiresolution transform */ 
    mri_comp_trans  (Imag, &Median, ImagRec, Nl, Nc, Nbr_Plan, 
                     Noise, NSigma, SignalQuantif, Tab_Quant, 
                     SupIsol, NoiseInData, Rec, MedianWindowSize, SupNeg);

    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 

    /* codes each scale */
    for (s = Nbr_Plan-1; s >= 0; s--)
    {
        Plan = Median + Tab_Pos[s];
        Sigmaf = Tab_Quant[s];
        encode(FileDes, Plan, Tab_Nl[s], Tab_Col[s], Sigmaf);       
        if (Verbose == True)
        {
           cout << "Set " << s+1 << "  NbrCoef = " 
                <<    Tab_Nl[s]*Tab_Col[s] 
                << " Quant = " <<   Sigmaf <<  endl;
        }
    }
    i_free(Median);
    delete [] Tab_Quant;
}

/****************************************************************************/
// decode_med_pyr
void mri_decompress (int **Imag, int *Nl, int *Nc, int *Nbr_Plan, 
                     int Resol, FILE *FileDes, int &Nli, int &Nci, Bool Verbose)
{
    int *Median;
    int Tab_Nl[MAX_SCALE], Tab_Col[MAX_SCALE], Tab_Pos[MAX_SCALE];
    int i,s, *Plan, *PlanMed, Size;
    int Nls, Ncs;
    float Level;
    float Tab_Level[MAX_SCALE];

    /* read the original image size */
    *Nl =readint(FileDes); /* x size of image */
    *Nc =readint(FileDes); /* y size of image */
    Nli = *Nl;
    Nci = *Nc;
    
    /* read the number of scales */
    *Nbr_Plan=readint(FileDes); 

    if (Verbose == True)
    {
       printf ("DECODING pyramid: nbr_scale = %d\n", *Nbr_Plan);
       printf ("Originale image: Nl = %d, Nc = %d\n", *Nl, *Nc);
    }

    /* index position */
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, *Nl, *Nc, *Nbr_Plan);

    /* Resol must be < Nbr_Plan and > 0 */
    if ((Resol < 0) || (Resol >= *Nbr_Plan))
    {
      fprintf (stderr, "Error: bad parameter Resol: Resol = %d\n", Resol);
      fprintf (stderr, "       Resol must verify 0 <= Resol < %d\n", *Nbr_Plan);
      exit (0);
    }

    if (Resol > 0)
    {
        *Nl = Tab_Nl[Resol];
        *Nc = Tab_Col[Resol];
        *Nbr_Plan -= Resol;
         mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, *Nl, *Nc, *Nbr_Plan);
    }
    if (Verbose == True) printf ("\n Create an image of size (%d, %d)\n", *Nl, *Nc);

    Size = mri_size_medpyr (*Nl, *Nc, *Nbr_Plan);
    // *Imag = i_vector_alloc (*Nl * *Nc);
    // Median = i_vector_alloc (Size);
    *Imag = i_alloc(*Nl * *Nc);
    Median = i_alloc(Size);

    /* decodes each scale */
    for (s = *Nbr_Plan - 1; s >= 0; s--)
    {
        PlanMed = Median + Tab_Pos[s];
        decode(FileDes, &Plan, &Nls, &Ncs, &Level);
        Tab_Level[s] = Level;
        if (Verbose == True) 
             printf ("\nScale %d: Nls = %d, Ncs = %d, Scaling = %5.2f\n", 
                  s+1, Nls, Ncs, Level);
        /* print_info_pict_i (Plan, Nls, Ncs);*/

        Tab_Level[s] = Level;
        for (i = 0; i < Nls*Ncs; i++) PlanMed[i] = Plan[i];
        i_free (Plan);
    }

    /* image reconstruction */
    mri_pyrmedian_rec (Median, *Imag, *Nl, *Nc, *Nbr_Plan, Tab_Level);
    i_free(Median);
}

/****************************************************************************/
