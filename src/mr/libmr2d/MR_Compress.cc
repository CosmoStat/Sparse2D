/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_Compress.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image compression using multiresolution
**    -----------  
**       
*******************************************************************************
**
** void mr_compress (Ifloat &Imag, int Nbr_Plan, float &Noise, float NSigma, 
**              float SignalQuant, int MaxIter, FILE *FileDes, Ifloat &Residual)
**
** Imag : IN image to code
** Nbr_Plan: IN number of scales used for the coding
** Noise: IN-OUT : noise estimation
** NSigma: IN Detection at NSigma
** SignalQuant: IN Quantification of the coefficient w by:
**                  c = w / (Sigma_Noise*SignalQuant)
** MaxIter: number of iterations
** FileDes: file descriptor
** Residual: image resiudal
**
*******************************************************************************
**
** void mr_decompress (Ifloat &Imag,  int Resol, FILE *FileDes)
**
** Imag: decompressed image
** Resol: resolution we want to extract
** FileDes: file descriptor
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
**                        float Noise, float NSigma, float N_SigmaNoise, 
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

// static char sccsid[] = "@(#)MR_Compress.cc 3.2 96/06/13 CEA 1995 @(#)";

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Noise.h"
// #include "MR_Filter.h"
#include "IM_CompTool.h"
#include "MR_Sigma.h"
#include "IM_Comp.h"
#include "MR_Comp.h"
#include "OptMedian.h"

#define PRINT_DATA 1
#define WRITE_DATA 0

static Bool VerboseComp=False;

/***************************************************************************/

static void kill_isol_sup(Iint *Median, MultiResol &MR_Data, float *Tab_Level,
                          float Noise_Ima, float N_Sigma)
{
      int i,j,s = 0;
      int Nls = MR_Data.band(s).nl();
      int Ncs = MR_Data.band(s).nc();
      float Level = Noise_Ima * mr_level_noise (MR_Data, s, N_Sigma);
      if (Tab_Level[s] < Level)
      {
          for (i =1; i < Nls-1; i++)
          for (j =1; j < Ncs-1; j++)
          if  (Median[s](i,j) != 0)
          {
              if (    !Median[s](i-1,j) && !Median[s](i-1,j-1) 
                        && !Median[s](i-1,j+1) && !Median[s](i,j-1) 
                        && !Median[s](i,j+1) && !Median[s](i+1,j-1)
                        && !Median[s](i+1,j) && !Median[s](i+1,j+1))
                    Median[s](i,j) = 0;
          }
      }
}

/***************************************************************************/


static void dilate_sup (Iint *Median, MultiResol &MR_Data, 
                 int Nbr_Plan, float *Tab_Level, 
                 float Noise_Ima, float N_Sigma)  
{
    int s,i,j;
    for (s = 0; s < Nbr_Plan-1; s++)
    {
        int Nls = MR_Data.band(s).nl();
        int Ncs = MR_Data.band(s).nc();
        float Level = Noise_Ima * mr_level_noise (MR_Data, s, N_Sigma);

        if (Tab_Level[s] < Level)
        {
            for (i =1; i < Nls-1; i++)
            for (j =1; j < Ncs-1; j++)
            {
               if ( Median[s](i,j) == 0)
               {
                   if (   ( fabs(MR_Data(s,i-1,j)) > Level)
                        || (fabs(MR_Data(s,i-1,j-1)) > Level)
                        || (fabs(MR_Data(s,i-1,j+1)) > Level)
                        || (fabs(MR_Data(s,i,j-1)) > Level)
                        || (fabs(MR_Data(s,i,j+1)) > Level)
                        || (fabs(MR_Data(s,i+1,j-1)) > Level)
                        || (fabs(MR_Data(s,i+1,j))   > Level)
                        || (fabs(MR_Data(s,i+1,j+1)) > Level))
                    Median[s](i,j) += quant( MR_Data(s,i,j), Tab_Level[s]);
                }
            }
        }
    }
}
/***************************************************************************/

static void mr_quantif (Iint *Median, MultiResol &MR_Data, 
                 int Nbr_Plan, float *Tab_Level, 
                 float Noise_Ima, float N_Sigma, Bool SupNeg)
{ 
    int s,i,j;
    float Level;
    int *ptr;

    for (s = 0; s < Nbr_Plan; s++)
    {

        int Nls = MR_Data.band(s).nl();
        int Ncs = MR_Data.band(s).nc();
        Level = Noise_Ima * mr_level_noise (MR_Data, s, N_Sigma);

        ptr = Median[s].buffer();
         
        for (i =0; i < Nls; i++)
        for (j =0; j < Ncs; j++)
        {
            if (Tab_Level[s] > 0)
            {
            if (MR_Data(s,i,j) >= Level)
                  Median[s](i,j)=(int) (MR_Data(s,i,j) / Tab_Level[s] + 0.5);
            else if ((SupNeg == False) && (MR_Data(s,i,j) < -Level))
                 Median[s](i,j)= (int) (MR_Data(s,i,j) / Tab_Level[s] - 0.5);
            }
            else 
            {
               if (s != Nbr_Plan-1) Median[s](i,j) = 0;
               else 
               {
                  float Level1 =  Noise_Ima / 10.;
                  Median[s](i,j) = quant(MR_Data(s,i,j), Level1);
               }
            }
        }
    } 
}

/****************************************************************************/

static void mr_invers_quantif (Iint *Median, MultiResol &MR_Data, 
                 int Nbr_Plan, float *Tab_Level)
{
    int s,i,j;

    for (s = 0; s < Nbr_Plan; s++)
    {
        int Nls = MR_Data.scale(s).nl();
        int Ncs = MR_Data.scale(s).nc();

        for (i =0; i < Nls; i++)
        for (j =0; j < Ncs; j++)
              MR_Data(s,i,j) = Median[s](i,j) * Tab_Level[s];
    }
}

/****************************************************************************/

static void mr_median_tr(Ifloat &Imag, MultiResol &MR_Data)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int Nls, Ncs;
    int s = 0;
    int si, sj,ei,ej;
    int k,l,i,j,ind_fen;
    int Nbr_Plan = MR_Data.nbr_scale();
    int Window_Size=MR_Data.MedianWindowSize;
    int Window2 = (Window_Size - 1) / 2;
    extern float hmedian(float *ra, int n);
    float *fenetre;
    extern float select(unsigned long k, unsigned long n, float arr[]);
    int  Size = Window_Size*Window_Size;

    fenetre = new float [Size];
 
    MR_Data.band(0) = Imag;
    for (s = 0; s < Nbr_Plan -1; s++)
    {
        Nls = MR_Data.band(s+1).nl();
        Ncs = MR_Data.band(s+1).nc();
        Nl = MR_Data.band(s).nl();
        Nc = MR_Data.band(s).nc();

        for (i = 0; i < Nls; i++) 
        {
           si = 2*i - Window2;
           ei = 2*i + Window2;
           if (si < 0) si = 0;
           if (ei >= Nl) ei = Nl-1;
           for (j = 0; j < Ncs; j++) 
           {
              ind_fen = 0;
              sj = 2*j - Window2;
              ej = 2*j + Window2;
              if (sj < 0) sj = 0;
              if (ej >= Nc) ej = Nc-1;

              for (k = si; k <= ei; k ++)
              for (l = sj; l <= ej; l ++) 
	         fenetre[ind_fen++] = MR_Data(s,k ,l);
              // MR_Data(s+1, i,j) = hmedian(fenetre, ind_fen); 
              if (ind_fen == 9)  MR_Data(s+1, i,j) = opt_med9(fenetre);  
              else MR_Data(s+1, i,j) = get_median(fenetre, ind_fen);
	    }
        }
        im_increase_size_2 (MR_Data.band(s+1), Imag);

        for (i = 0; i < Nl; i++) 
        for (j = 0; j < Nc; j++) MR_Data(s,i,j) -= Imag(i,j);
        Imag.resize(Nls, Ncs);
    }
    delete [] fenetre;
}
    
/****************************************************************************/

static void mr_comp_trans  (Ifloat &Imag, Iint *PyrMedian, Ifloat &ImagRec,
                     int Nbr_Plan, float & Noise_Ima, float NSigma, 
                     float SignalQuant, int Budget, float Tab_Level[], 
                     Bool SupIsol, Bool NoiseInData,
                     int MedianWindowSize, Bool SupNeg, int & MeanLastScale)
{
    int s;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    type_transform Transform = TM_PYR_MEDIAN;
    MultiResol MR_Data (Nl, Nc, Nbr_Plan, Transform, "MR_Imag");
    extern float  mr_tab_noise(int s);
    int Nls, Ncs;
    MR_Data.MedianWindowSize = MedianWindowSize;
    float Noise;
     float *Tab_Sigma = new float[Nbr_Plan];
          
    /* computes the multiresolution transform */
    ImagRec = Imag;
    mr_median_tr(ImagRec, MR_Data);
 
    /* creates the  pyramid */
    for (s = 0; s < Nbr_Plan; s++)
    {
        Nls = MR_Data.scale(s).nl();
        Ncs = MR_Data.scale(s).nc();

        PyrMedian[s].alloc(Nls,Ncs,"pyr int");
        PyrMedian[s].init();
    }

    /* Noise estimation in the Data */
    noise_compute (MR_Data);
    if ((Noise_Ima < FLOAT_EPSILON) && (NoiseInData == True))
                                Noise_Ima = mr_noise_estimation (MR_Data);
    if (NoiseInData == False) Noise = 1.;
    else Noise = Noise_Ima;
    
    /* search the signification level at each scale */
    float *TabMean= new float[Nbr_Plan];
    MeanLastScale=0;
    
    // Sigma at each scale
    for (s = 0; s < Nbr_Plan; s++) Tab_Sigma[s] =  mr_tab_noise (s) * Noise;
    s = Nbr_Plan-1;
    Tab_Sigma[s] = Tab_Sigma[s-1];
    
    if (Budget < 1)
    {
       for (s = 0; s < Nbr_Plan; s++)
              Tab_Level[s] = Tab_Sigma[s]*
                        get_quant_param(s, MR_Data.nbr_scale(), 
                                         SignalQuant, DEFAULT_SIGNAL_QUANTIF);
    }
    else
    {
        if (Noise > 1.)
               for (s = 0; s < Nbr_Plan; s++) Tab_Sigma[s] /= Noise;
        mr_getquant(MR_Data, Tab_Sigma, Budget, Tab_Level, TabMean);
        MeanLastScale = quant(TabMean[Nbr_Plan-1],1.);
    }  
    
    if (VerboseComp == True) 
    {
       fprintf (stderr, "Estimated Noise Standard Deviation = %f\n", Noise_Ima);
       for (s = 0; s < Nbr_Plan; s++)
             printf ("Scale %d, Level = %5.5f\n", s+1, Tab_Level[s]);
    }
    
    /* quantification */
    if (NoiseInData == True)
        mr_quantif (PyrMedian, MR_Data, Nbr_Plan,Tab_Level,Noise_Ima, 
                    NSigma, SupNeg);
    else mr_quantif (PyrMedian, MR_Data, Nbr_Plan,Tab_Level, 0., 
                    NSigma, SupNeg);
                    
    /* kill isolated pixels in the first scale */
    if (SupIsol == True) kill_isol_sup(PyrMedian,  MR_Data,  Tab_Level,
                                     Noise_Ima, NSigma);
  
    /* Support dilatation */
    if (NoiseInData == True) 
       dilate_sup (PyrMedian, MR_Data, Nbr_Plan,  Tab_Level, Noise_Ima, NSigma);
        
    /* inverse transform */
    mr_invers_quantif (PyrMedian, MR_Data, Nbr_Plan, Tab_Level);
   
    /* Reconstruct the image */
    MR_Data.recons (ImagRec);       
    
    delete [] TabMean;
    delete [] Tab_Sigma;
}

/****************************************************************************/

void encode_null_block(FILE *outfile, int *a, int nx, int ny)
{
   long size_enc;

   writeint(outfile, 0);
   int Bef = ftell(outfile);        
   writeint(outfile, 1);	/* size of image */
   writeint(outfile, 1);
   writefloat(outfile, 0.);
   writefloat(outfile, *a); /* scale factor for digitization	*/
   // writeint(outfile, ny);
   // writeint(outfile, nx);
 
   int Aft = ftell(outfile);
   size_enc =  Aft - Bef;
   long Val = -size_enc-4;
   int nf = fseek(outfile, Val, SEEK_CUR);
   writeint(outfile, size_enc);
   nf = fseek(outfile, size_enc,SEEK_CUR);
}

/********************************************************************/

void mr_compress (Ifloat &Imag, int Nbr_Plan, float &Noise, float NSigma, 
                float SignalQuant, int Budget, FILE *FileDes, 
                Ifloat &Residual, Bool SupIsol, Bool NoiseInData,
                int MedianWindowSize, Bool SupNeg, Bool Verbose)
{
    float Tab_Level[MAX_SCALE];
    int s;
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    Iint *PyrMedian;
    int MeanLastScale;
    
    PyrMedian = new Iint [Nbr_Plan];
    VerboseComp = Verbose;
    /* calcules  the multiresolution transform */ 
    mr_comp_trans  (Imag, PyrMedian, Residual, Nbr_Plan, Noise, 
                    NSigma, SignalQuant, Budget, Tab_Level, 
                    SupIsol, NoiseInData, MedianWindowSize, 
                    SupNeg, MeanLastScale);

    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 

    /* codes each scale */
    for (s = Nbr_Plan-1; s >= 0; s--)
    {
       if (Tab_Level[s] > 0.)
           encode(FileDes, PyrMedian[s].buffer(), 
               PyrMedian[s].nl(), PyrMedian[s].nc(), Tab_Level[s]);
        else 
        {
           if (s != Nbr_Plan-1) MeanLastScale = 0;
           encode(FileDes, &MeanLastScale, 1,1, (float) 0.);
           //encode_null_block(FileDes, &MeanLastScale, 
	   //                  PyrMedian[s].nl(), PyrMedian[s].nc());
           // writeint(FileDes, PyrMedian[s].nl());
           // writeint(FileDes, PyrMedian[s].nc());
        }     
        if (Verbose == True)
        {
           extern long size_enc;
           cerr << "Set " << s+1 << "  NbrCoef = " 
                <<   PyrMedian[s].nl()*PyrMedian[s].nc() 
                << " Quant = " <<  Tab_Level[s] << " Size = " << size_enc << endl;
        }
     }
     
     delete []PyrMedian;
}

/****************************************************************************/

void mr_decompress (Ifloat &Imag,  int Resol, FILE *FileDes, int &Nli, int &Nci, Bool Verbose)
{
    int i,j,k,s, *Plan;
    int Nls, Ncs;
    float Level;
    int Nl, Nc, Nbr_Plan;
    // type_transform Transform = TM_PYR_MEDIAN;

    /* read the original image size */
    Nl =readint(FileDes); /* x size of image */
    Nc =readint(FileDes); /* y size of image */
    Nli = Nl;
    Nci = Nc;
    
    /* read the number of scales */
    Nbr_Plan=readint(FileDes); 

    //if (verbose)
    //{ 
    //   fprintf (stderr, "DECODING pyramid: nbr_scale = %d\n", Nbr_Plan);
    //   fprintf (stderr, "Original image: Nl = %d, Nc = %d\n", Nl, Nc);
    //}

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
        for (s = 0; s < Resol; s++)
        {
            Nl = (Nl+1) / 2;
            Nc = (Nc+1) / 2;
        }
        Nbr_Plan -= Resol;
    }
    //if (Verbose==True) fprintf (stderr,"\n Create an image of size (%d, %d)\n\n", Nl, Nc);
    Imag.alloc (Nl, Nc, "decomp image");
    // MultiResol MR_Data (Nl, Nc, Nbr_Plan, Transform, "MR_Imag");
    Ifloat Buffer(Nl, Nc, "Buffer mr_decomp");

    /* decodes each scale */
    extern long size_enc;
    for (s = Nbr_Plan - 1; s >= 0; s--)
    {
        decode(FileDes, &Plan, &Nls, &Ncs, &Level);
        if (Level < FLOAT_EPSILON)
        {
           // Nls = readint(FileDes);  
           // Ncs = readint(FileDes);
	   Nls = Nl;
	   Ncs = Nc;
	   for (k = 0; k < s; k++)
	   {
	      Nls = (Nls+1) / 2;
	      Ncs = (Ncs+1) / 2;
	   }
        }
        Buffer.resize (Nls, Ncs);
        if (s != Nbr_Plan - 1) im_increase_size_2 (Imag, Buffer);
        Imag.resize (Nls, Ncs);
        int SizeScale = size_enc;
        if (Verbose == True) 
          fprintf (stderr, "Scale %d: Nls = %d, Ncs = %d, Quantif = %5.2f, Size = %d\n", 
                  s+1, Nls, Ncs, Level, SizeScale);
        if (Level < FLOAT_EPSILON)
        { 
           for (i = 0; i < Nls; i++)
           for (j = 0; j < Ncs; j++) Imag(i,j)= Plan[0];
        }
        else
        {
        for (i = 0; i < Nls; i++) 
        for (j = 0; j < Ncs; j++) 
        {
                //    MR_Data(s,i,j) = Plan[i*Ncs+j]*Level;
            if (s != Nbr_Plan - 1) Imag(i,j)=Buffer(i,j)+ Plan[i*Ncs+j]*Level;
            else Imag(i,j) =  Plan[i*Ncs+j]*Level;
        }
        }
        i_free(Plan);
    }
    /* image reconstruction */
    // MR_Data.recons (Imag);
}

/****************************************************************************/

