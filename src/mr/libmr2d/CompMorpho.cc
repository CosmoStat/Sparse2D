/******************************************************************************
**                   Copyright (C) 1997 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Simona Mei && Jean-Luc Starck
**
**    Date:  97/07/14 
**    
**    File:  CompMorpho.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image compression using mathematical morphology
**    -----------  
**       
*******************************************************************************
**
** void morpho_compress (Ifloat &Imag, FILE *FileDes, int Nbr_Plan,
**                       float &Noise_Ima, float N_Sigma  
**                       int Elstr_Size, int Elstr_Shape)
**
** Imag : IN image to code
** Nbr_Plan: IN number of scales used for the coding
** Noise_Ima: IN-OUT : noise estimation
** NSigma: IN Detection at NSigma
** FileDes: file descriptor
** Elstr_Size: Size of the structural element
** Elstr_Shape: Shape of the structural element (square or cercle) 
**
*******************************************************************************
**
** void morphoi_compress (int *Imag, int Nl, int Nc, int Nbr_Plan, FILE *FileDes, 
**                        float & Noise, float NSigma, float N_SigmaNoise, 
**                        int Elstr_Size, int Elstr_Shape, int MedianWindowSize)                          **
** Imag : IN image to code
** Nl, Nc: Image dimensions
** Nbr_Plan: IN number of scales used for the coding
** Noise: IN-OUT : noise estimation
** NSigma: IN Detection at NSigma
** FileDes: file descriptor
** Elstr_Size: Size of the structural element
** Elstr_Shape: Shape of the structural element (square or cercle)
** MedianWindowSize: Median window for background estimation by MR
**
*******************************************************************************
**
** void morpho_decompress (Ifloat &Imag, int N_Iter,  FILE *FileDes, 
**                        int Elstr_Size, int Elstr_Shape)
**
** Imag: decompressed image
** FileDes: file descriptor
** Elstr_Size: Size of the structural element
** Elstr_Shape: Shape of the structural element (square or cercle)
** N_Iter: Morphological transformation iteration number
**
*******************************************************************************
**
** void morphoi_decompress (Int *Imag, int N_Iter, int Nc, int Nl, FILE *FileDes, 
**                        int Elstr_Size, int Elstr_Shape)
**
** Imag: decompressed image
** FileDes: file descriptor
** Elstr_Size: Size of the structural element
** Elstr_Shape: Shape of the structural element (square or cercle)
** N_Iter: Morphological transformation iteration number
**
******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Noise.h"
// #include "MR_Filter.h"
#include "IM_CompTool.h"
#include "MR_Sigma.h"
#include "IM_Sigma.h"

#define PRINT_DATA 1
#define WRITE_DATA 0

//Morphomate part
static  int MaxBuffer = 100;
 
/****************************************************************************/

void morpho_compress (Ifloat &Imag, FILE *FileDes, 
                     float &Noise_Ima, float N_Sigma, 
                     int ElStr_Size, int Elstr_Shape, 
                     Bool UseMR_for_BGR, int NPixBgr, Bool Verbose)
{
    int i,j,k,l;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int N_Iter;
    Ifloat Imag_Er(Nl,Nc,"er");
    Ifloat Imag_Dil(Nl,Nc,"dil");
    Ifloat Germe(Nl,Nc,"ger");
    int *PlanMorpho = i_alloc(Nc*Nl);
    float ThresholdLevel;
 
    writeint(FileDes, Nl);
    writeint(FileDes, Nc);

    /* Background and Noise calculation*/
    if (UseMR_for_BGR == True)
    {
       writeint(FileDes, 1);
       int Nbr_Plan = 1;
       type_transform Transform = TM_PYR_MEDIAN;  
       int Ns = MIN(Nl,Nc);
       while (NPixBgr < Ns) 
       {
          Nbr_Plan ++;
          Ns /= 2;
       }
       if (Nbr_Plan < 3) Nbr_Plan = 3;
       if (Verbose == True)
            cout << "Number of scales for bgr calculation = " << Nbr_Plan << endl;
       MultiResol MR_Data(Imag.nl(),Imag.nc(),Nbr_Plan,Transform,"MR_Transform");
       writeint(FileDes, Nbr_Plan);
       MR_Data.transform(Imag);           

       if (Noise_Ima < FLOAT_EPSILON) Noise_Ima = mr_noise_estimation (MR_Data);
       else noise_compute (MR_Data);
       writefloat(FileDes,  Noise_Ima);
      for (k=0;k<MR_Data.nbr_scale()-1;k++) MR_Data.scale(k).init();
      Ifloat Bkg(Imag.nl(),Imag.nc(),"Background");
      MR_Data.recons(Bkg);

      ThresholdLevel = (N_Sigma+1)*Noise_Ima;
      for (i=0; i < Nl*Nc;i++)
      {
           if (Imag(i) < Bkg(i)+ ThresholdLevel) Imag(i) = 0;
           else Imag(i) = (int) ((Imag(i) - Bkg(i)) / Noise_Ima);
      }
      
      // encode the backgound
      for (i=0; i < MR_Data.scale(Nbr_Plan-1).n_elem(); i++)
         PlanMorpho[i] = (int)((MR_Data.scale(Nbr_Plan-1))(i) / Noise_Ima);
      encode(FileDes, PlanMorpho, MR_Data.scale(Nbr_Plan-1).nl(),
             MR_Data.scale(Nbr_Plan-1).nc(), Noise_Ima);
    }
    else
    {
       writeint(FileDes, 0);
       // background and sigma calculation by iteration
       float Av_Im = 0;
       if (Noise_Ima < FLOAT_EPSILON) sigma_clip(Imag, Av_Im, Noise_Ima); 
       else 
       {
          float NSig;
          sigma_clip(Imag, Av_Im, NSig);
       }
       writefloat(FileDes, Av_Im);
       writefloat(FileDes, Noise_Ima); 
       ThresholdLevel = Av_Im+Noise_Ima*(N_Sigma);
       for (i = 0; i < Nl*Nc; i++)
       {
           if (Imag(i) < ThresholdLevel) Imag(i) = 0;
           else Imag(i) = (int) ((Imag(i) - Av_Im) / Noise_Ima);
       }
    }
    
    
    /**delete the isolated pixels***/

    for (i=0; i < Nl;i++)
    for (j=0; j < Nc;j++)
    {
        float Max=0;
        for (k = i - 1; k <= i + 1; k++)
        for (l = j - 1; l <= j + 1; l++)
             if (Imag(k,l, I_CONT) > 0) Max = 1;
        if (Max == 0) Imag(i,j)=0;            
    }
       
   /* Calculation of the planes with a grey level mathematical
      morphology transformation
   */
   N_Iter=1;
   while (max(Imag) > 0) 
   {   
       writeint(FileDes, N_Iter);
       Germe = Imag;
       if (Elstr_Shape == 1) /* Calculation with a cercle */
       {
       //cout << "cercle : " << ElStr_Size << endl;
          morpho_cercle_erosion(Imag, Imag_Er, ElStr_Size);
          morpho_cercle_dilation(Imag_Er,Imag_Dil,ElStr_Size);
       }
       else /*Calculation with a square*/
       {
       //cout << "square  : " << ElStr_Size << endl;
          morpho_erosion(Imag, Imag_Er, ElStr_Size);
          morpho_dilation(Imag_Er,Imag_Dil,ElStr_Size);
       }
       Germe-=Imag_Dil;   
       for (i=0;i<Nl;i++)
       for (j=0;j<Nc;j++)
              PlanMorpho[i*Nc+j] = iround(Germe(i,j));
       encode(FileDes, PlanMorpho, Nl, Nc, 1.);
       N_Iter++;
      Imag = Imag_Er;
   }
   
   writeint(FileDes, 0);
   if (Verbose == True) fprintf (stderr, "Iteration number: = %d \n", N_Iter-1);
   i_free(PlanMorpho);
}

// DECOMPRESSION MORPHOMATE
 /****************************************************************************/

void morpho_decompress (Ifloat &Imag, int & N_Iter, FILE *FileDes, 
                        int ElStr_Size, int Elstr_Shape, Bool Verbose)
{
    int Nls,Ncs, i,s, *Plan;
    int Nl,Nc;
    float Level, Noise, Avg=0.;
    N_Iter = 0;
     
    /* read the original image size */
    Nl =readint(FileDes); /* x size of image */
    Nc =readint(FileDes); /* y size of image */
        
    Imag.alloc (Nl, Nc, "decomp image");
    Imag.init();
    
    Ifloat *PlanBuffer;
    PlanBuffer = new Ifloat [MaxBuffer];
    Ifloat Imag_Dil(Nl,Nc,"Imag_Dil");
    Ifloat Bkg(Nl, Nc,"Background");

    if (Verbose == True)
    {
       fprintf (stderr, "Original image: Nl = %d, Nc = %d\n", Nl, Nc);
       fprintf (stderr, "\n Create an image of size (%d, %d)\n", Nl, Nc);
    }
     
    /* decodes each plan */

    // read the background
    int Bgr = readint(FileDes);
    if (Bgr == 1)
    {
       int DNbrPlan = readint(FileDes);       
       Noise = readfloat(FileDes);       
       decode(FileDes, &Plan, &Nls, &Ncs, &Level);
       MultiResol MR_Data(Nl,Nc,DNbrPlan, TM_PYR_MEDIAN,"MR_Transform");
       for (i=0; i < MR_Data.scale(DNbrPlan-1).n_elem(); i++)
          MR_Data.scale(DNbrPlan-1)(i) = (float) Plan[i] * Level;
       MR_Data.recons(Bkg);
       i_free(Plan);
    }
    else 
    {
       Avg=readfloat(FileDes);
       Noise =readfloat(FileDes);       
       Bkg.init(Avg); 
    }
    
    // Read the planes
    int ValNiter = readint(FileDes);
    while (ValNiter != 0)
    {
       N_Iter  = ValNiter;
       s = N_Iter -1;
       if (N_Iter >= MaxBuffer)
       {
          cerr << "Error: buffer for morpho compression too small ... " << endl;
          exit(-1);
       }
       PlanBuffer[s].alloc(Nl,Nc,"Plan");
       PlanBuffer[s].init();
       decode(FileDes, &Plan, &Nls, &Ncs, &Level);       
       if (Verbose == True) 
          fprintf (stderr, "Iter %d: Nl = %d, Nc = %d, Scaling = %5.2f\n", 
                                              N_Iter , Nls, Ncs, Level);

       for (i = 0; i < Nls*Ncs; i++)  PlanBuffer[s](i) = float(Plan[i]) * Level;
       ValNiter = readint(FileDes);
       i_free(Plan);
    }
    if (Verbose == True)
      fprintf (stderr, "Iteration number: = %d \n",  N_Iter);
      
    // Sum the planes to constract the final decompressed image  
    if (N_Iter > 0)
    {
    Imag = PlanBuffer[N_Iter-1];
    for (s = N_Iter-1 ;s > 0;s--)
    {   
       // Structurant element is a cercle
       if (Elstr_Shape ==1) morpho_cercle_dilation(Imag,Imag_Dil,ElStr_Size);
       //Structurant element is a square 
       else morpho_dilation(Imag,Imag_Dil,ElStr_Size); 
       Imag=Imag_Dil;
       
       Imag+=PlanBuffer[s-1];
    }
    }
    for (i=0; i <Nl*Nc; i++) Imag(i) = Imag(i)*Noise + Bkg(i);
    
    delete [] PlanBuffer;
    
}
 
/***************************************************************************/
 
void morphoi_quant_compress(int *Imag, int Nl, int Nc, FILE *FileDes, 
                            int ElStr_Size,  int Elstr_Shape, Bool Verbose)
{
    int i =0;  
    int N_Iter=0;  
    int *Imag_Er;
    int *Imag_Dil;
    int *Germe;
     
    Imag_Er = i_alloc(Nl*Nc);
    Imag_Dil = i_alloc(Nl*Nc);
    Germe =  i_alloc(Nl*Nc);
    Iint Dat(Nl, Nc, "test");
    for (i = 0; i < Nl*Nc; i++)
   { 
      Imag_Er[i]=0;
      Imag_Dil[i]=0;
      Germe[i]=0;
   }
   
   /* calcules the gray level morpho transform plans*/
   N_Iter=1;

   Bool ContLoop = True;
   while (ContLoop == True) 
   {     
     writeint(FileDes, N_Iter);
     if (Elstr_Shape == 1)
     {
         morphoi_cercle_erosion (Imag, Imag_Er, Nl, Nc,ElStr_Size);
         morphoi_cercle_dilation (Imag_Er,Imag_Dil,Nl, Nc,ElStr_Size);
     }
     else
     {
        morphoi_erosion (Imag, Imag_Er, Nl, Nc,ElStr_Size);
        morphoi_dilation (Imag_Er,Imag_Dil,Nl, Nc,ElStr_Size);
     }
     for (i=0; i<Nl*Nc; i++) Germe[i]=Imag[i]-Imag_Dil[i];
     
     encode(FileDes, Germe, Nl, Nc, 1.);
     N_Iter++;
     int count=0;
     for (i=0; i<Nl*Nc; i++)
     {
        Imag[i]=Imag_Er[i];
        if (Imag_Er[i]==0) count++; 
     }
     if (count == Nl*Nc) ContLoop = False;    
   }   
  
   writeint(FileDes, 0);
   if (Verbose == True) fprintf (stderr, "Iteration number: = %d \n", N_Iter-1); 
   
   i_free(Imag_Er);   
   i_free(Imag_Dil);
   i_free(Germe);
}

/***************************************************************************/

void morphoi_compress (int *Imag, int Nl, int Nc, FILE *FileDes, 
                       float & Noise, float N_SigmaNoise, int ElStr_Size, 
                       int Elstr_Shape, Bool UseMR_for_BGR, int NPixBgr,
		       Bool Verbose)
{
    int i =0;  
 
    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
   
   /* Background and Noise calculation*/
   if (UseMR_for_BGR == True)
   {
       writeint(FileDes, 1);
       int Nbr_Plan = 1; 
       int Tab_Nl[MAX_SCALE], Tab_Col[MAX_SCALE], Tab_Pos[MAX_SCALE];
       int *Bkg= i_alloc(Nl*Nc);   
       int *Median;
       int Ns = MIN(Nl,Nc);
       while (NPixBgr < Ns) 
       {
          Nbr_Plan ++;
          Ns /= 2;
       }
       if (Nbr_Plan < 3) Nbr_Plan = 3;
       if (Verbose == True)
            cout << "Number of scales for bgr calculation = " << Nbr_Plan << endl;
       mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);  
       mri_pyrmedian_transform(Imag, Nl, Nc, &Median, Nbr_Plan, 5);
       if (Noise < FLOAT_EPSILON) Noise = sigma_clip_int(Median,Nl,Nc,0)/0.970291;
       writeint(FileDes,Nbr_Plan);
       writefloat(FileDes,Noise);
          

      // encode the backgound
      int Nl_bgr = Tab_Nl[Nbr_Plan-1];
      int Nc_bgr = Tab_Col[Nbr_Plan-1];
      int *Bgr = Median + Tab_Pos[Nbr_Plan-1];
      for (i=0; i < Nl_bgr*Nc_bgr; i++)  Bkg[i] = (int) (Bgr[i]/Noise);
 
      encode(FileDes, Bkg, Nl_bgr, Nc_bgr, Noise);
      
      morphoi_noise_rec(Median,Bkg,Nl,Nc,Nbr_Plan);
      float ThresholdLevel =  Noise*(N_SigmaNoise);
      for (i = 0; i < Nl*Nc; i++)
      {
          if (Imag[i] < ThresholdLevel+Bkg[i]) Imag[i] = 0;
          else Imag[i] = (int) ((Imag[i] - Bkg[i]) / Noise);      
      }
      i_free(Bkg);
      i_free(Median);
     }
    else
    {
       writeint(FileDes, 0);
       // background and sigma calculation by iteration
       float Av_Im = 0;
       if (Noise  < FLOAT_EPSILON) sigma_clip_int(Imag, Av_Im, Noise,Nl,Nc); 
       else 
       {
          float NSig;
          sigma_clip_int(Imag, Av_Im, NSig,Nl,Nc);
       }
       writefloat(FileDes, Av_Im);
       writefloat(FileDes, Noise);
       float ThresholdLevel = Av_Im+Noise*(N_SigmaNoise);
       for (i = 0; i < Nl*Nc; i++)
       {
           if (Imag[i] < ThresholdLevel) Imag[i] = 0;
           else Imag[i] = (int) ((Imag[i] - Av_Im) / Noise);
       }
    }
    morphoi_quant_compress(Imag, Nl, Nc, FileDes, ElStr_Size, Elstr_Shape, Verbose);  
}
 
/***************************************************************************/

void morphoi_quant_decompress(int *Imag, int & N_Iter, int Nl, int Nc,
                         FILE *FileDes, int ElStr_Size, int Elstr_Shape, Bool Verbose)
{
    int i,*Plan;
    float Level;
 
    /* read the original image size */
    int *Imag_Dil;
    Iint *PlanBuffer;
    PlanBuffer= new Iint [MaxBuffer];    
    Imag_Dil= i_alloc(Nl*Nc);
  
    if (Verbose == True)
    {
       fprintf (stderr, "Original image: Nl = %d, Nc = %d\n", Nl, Nc);
       fprintf (stderr, "\n Create an image of size (%d, %d)\n", Nl, Nc);
    }
    
    /* decodes each plan */
    N_Iter = readint(FileDes) - 1;
    Bool Cont = True;
    if (N_Iter < 0) Cont = False;
    
    int NbrLoop = 0;
    while (Cont == True)
    {
       NbrLoop ++;
       int Nls,Ncs;
       if (N_Iter >= MaxBuffer)
       {
          cerr << "Error: buffer for morpho compression too small ... " << endl;
          exit(-1);
       }
       
       decode(FileDes, &Plan, &Nls, &Ncs, &Level);
       PlanBuffer[N_Iter].alloc(Nls,Ncs,"plan");
       PlanBuffer[N_Iter].init();     
       if (Verbose == True) 
          fprintf (stderr, "Iter %d: Nl = %d, Nc = %d, Scaling = %5.2f\n", 
                                              N_Iter +1, Nls, Ncs, Level);
                                     
       for (i = 0; i < Nls*Ncs; i++)  PlanBuffer[N_Iter](i) =  Plan[i];
       N_Iter = readint(FileDes) - 1;
       if (N_Iter < 0) Cont = False;
       i_free(Plan);
    }
        
    for (i=0;i<Nl*Nc;i++)  Imag[i]=PlanBuffer[NbrLoop-1](i);
    for (int s =  NbrLoop-1; s > 0;s--)
    {          
       if (Elstr_Shape ==1) //case cercle
           morphoi_cercle_dilation(Imag,Imag_Dil, Nl, Nc,ElStr_Size);       
       else //case square
          morphoi_dilation(Imag,Imag_Dil, Nl, Nc, ElStr_Size);
       for (i=0;i<Nl*Nc;i++)
          Imag[i]=PlanBuffer[s-1](i)+Imag_Dil[i];
    }  
    
    delete [] PlanBuffer; 
    i_free(Imag_Dil); 
}

/**************************************************************************/

void morphoi_decompress(int *Imag, int & N_Iter, int Nl, int Nc,
                         FILE *FileDes, int ElStr_Size, int Elstr_Shape,
			 Bool Verbose)
{
    int i,*Plan;
    float Level,Noise;
 

    // read the background
    int Bgr = readint(FileDes);
    if (Bgr == 1)
    {
       int Nls,Ncs,NbrPlan = readint(FileDes);
       Noise = readfloat(FileDes);
       decode(FileDes, &Plan, &Nls, &Ncs, &Level);
       Ifloat Bkg(Nl,Nc,"Bkg");   
       MultiResol MR_Data(Nl,Nc,NbrPlan, TM_PYR_MEDIAN,"MR_Transform");
       for (i=0; i < MR_Data.scale(NbrPlan-1).n_elem(); i++)
          MR_Data.scale(NbrPlan-1)(i) = (float) Plan[i] * Level;
       MR_Data.recons(Bkg);
       morphoi_quant_decompress(Imag,N_Iter,Nl,Nc,FileDes,
                                ElStr_Size,Elstr_Shape, Verbose);
       for (i=0;i <Nl*Nc; i++) Imag[i] = (int)(Imag[i]*Noise + Bkg(i));
       i_free(Plan);
    }
    else 
    {
       float Avg=readfloat(FileDes); 
       Noise = readfloat(FileDes);    
       morphoi_quant_decompress(Imag,N_Iter,Nl,Nc,FileDes,
                                ElStr_Size,Elstr_Shape,Verbose);
       for (i=0;i <Nl*Nc; i++) Imag[i] = (int) (Imag[i]*Noise + Avg); 
    }    
}
    
    
