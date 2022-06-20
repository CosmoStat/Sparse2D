/******************************************************************************
**                   Copyright (C) 2001 CEA
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
**    File:  LiftingComp.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image compression/decompression using lifting
**    -----------  
**       
*******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "IM_CompTool.h"
#include "MR_Sigma.h"
#include "MR_Calloc.h"
#include "IM_Comp.h"
#include "MR_Comp.h"

/***************************************************/

static void get_info_band(int Nl, int Nc, int NbrBand, 
                          int *TabBandNl, int *TabBandNc,
                          int *TabResolNl,int *TabResolNc)
{
   int s;
   int Nl_s=Nl;
   int Nc_s=Nc;

   for (s=0; s < NbrBand-1; s+=3)
   {
       TabBandNl[s] = (Nl_s+1)/2;
       TabBandNc[s] = Nc_s/2;
       TabBandNl[s+1] = Nl_s/2;
       TabBandNc[s+1] = (Nc_s+1)/2;
       TabBandNl[s+2] = Nl_s/2;
       TabBandNc[s+2] = Nc_s/2;
       Nl_s = (Nl_s+1)/2;
       Nc_s = (Nc_s+1)/2;
       TabResolNl[s] = Nl_s;
       TabResolNc[s] = Nc_s;
       TabResolNl[s+1] = Nl_s;
       TabResolNc[s+1] = Nc_s;
       TabResolNl[s+2] = Nl_s;
       TabResolNc[s+2] = Nc_s;
   }     
   s = NbrBand-1;
   TabBandNl[s] = Nl_s;
   TabBandNc[s] = Nc_s;
   TabResolNl[s] = Nl_s;
   TabResolNc[s] = Nc_s;
}

/***************************************************/

void lifting_compress(Ifloat &Imag, int Nbr_Plan, FILE *FileDes,  
                      float NSigma, float SignalQuantif, Bool UnifQuant,
                      Bool Verbose, Bool RecImag)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int i,j,s;
    int Nl_s=Nl,Nc_s=Nc;
    int NbrBand = (Nbr_Plan-1)*3 + 1;
    float *TabQuantif = new float [NbrBand];
    int *TabBandNl = new int [NbrBand];
    int *TabBandNc = new int [NbrBand];
    int *TabResolNl = new int [NbrBand];
    int *TabResolNc = new int [NbrBand];
    float NSig = (NSigma > 0) ? NSigma: 1.;
     
    // FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    // SubBandFilter SB1D(SelectFilter, NORM_L2);
    // SB1D.Border = I_MIRROR;
    // Ortho_2D_WT WT2D(SB1D);
    type_lift LiftingTrans=TL_F79;
    Lifting Clift1D(LiftingTrans);
    Clift1D.Border = I_MIRROR;
    Ortho_2D_WT WT2D(Clift1D);

    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 
    WT2D.transform(Imag, Nbr_Plan);
    // cout << "WT2D OK" << endl;
    // Band image size computation
 
    // else cout << "  noise = " << Noise << endl;

    get_info_band(Nl,Nc,NbrBand,TabBandNl,TabBandNc,TabResolNl,TabResolNc);
    // extern long size_enc;
    int NbrCoef;
    float Quantif;
    int *qdata;
    for (s = NbrBand-1; s >= 0; s--)
    {      
       int Ind=0, Depi=0, Depj=0;
       details Detail;
       int Scale;
       Nl_s = TabBandNl[s];
       Nc_s =  TabBandNc[s];
       TabQuantif[s] = SignalQuantif;
       NbrCoef =  Nl_s*Nc_s;
       qdata = i_alloc(NbrCoef);
       band2scale(s,TO_MALLAT, NbrBand, Scale, Detail);
       switch (Detail)
       {
        case D_HORIZONTAL:
               Depj += TabResolNc[s];
               break;
        case D_DIAGONAL:  
               Depi += TabResolNl[s];
               Depj += TabResolNc[s];
               if (Scale == 0) TabQuantif[s] *= 2;
	       else TabQuantif[s] *= 1.5;
              break;
        case D_VERTICAL:  
               Depi += TabResolNl[s];
               break;
        case I_SMOOTH: break;
        case D_NULL:
        default: cerr << "Error: unknown detail" << endl;
                 exit(0);
                 break;
       }         
       if (Scale == 0) TabQuantif[s] *= 1.5;
       if (UnifQuant == True) Quantif = SignalQuantif;
       else Quantif =  TabQuantif[s];

       if (Quantif > 0) 
       {
         float LevelDetect = SignalQuantif*NSig;

     // cerr << "Set " << s+1 << "  LevelDetect  = " <<  LevelDetect 
     //     << " Depi = " << Depi << "   Depj = " << Depj << endl;

         for (i = 0; i < Nl_s; i++)
         for (j = 0; j < Nc_s; j++)
         {
             float Val = Imag(i+Depi,j+Depj);
             if (ABS(Val) < LevelDetect)  qdata[Ind++] = 0;
             else qdata[Ind++] = quant(Val, Quantif);
         }
         encode(FileDes,  qdata,  Nl_s,  Nc_s, Quantif);

         if (RecImag == True)
         {
            Ind=0; 
            for (i = 0; i < Nl_s; i++)
            for (j = 0; j < Nc_s; j++) 
                  Imag(i+Depi,j+Depj) = (float) qdata[Ind++] * Quantif;
         }
       }
       else
       { 
          qdata[0] =  0;  // quant(TabMean[s], 1.);
          encode( FileDes,  qdata, 1,1, (float) 0.);
          if (RecImag == True)
          {
             for (i = 0; i < Nl_s; i++)
             for (j = 0; j < Nc_s; j++) 
                  Imag(i+Depi,j+Depj) = (float) qdata[0];
          }
       }
       if (Verbose == True)
       {
           //extern long size_enc;
           // cerr << "Set " << s+1 << "  NbrCoef = " <<  NbrCoef 
           //     << " Quant = " << Quantif <<  " Size = " << size_enc << endl;
       }
       i_free(qdata);
     }
   if (RecImag == True) WT2D.recons(Imag, Nbr_Plan);
    
   delete [] TabQuantif;
   delete [] TabBandNl;
   delete [] TabBandNc;
   delete [] TabResolNl;
   delete [] TabResolNc;
}

/*********************************************************************/

void lifting_decompress(Ifloat &Imag, int Resolution, FILE *FileDes, 
                        int &Nli, int &Nci, Bool Verbose)
{
    int i,j,s;
    int Nl,Nc,Nl_s,Nc_s;
    int Nbr_Plan;
    int Resol = MAX(0,Resolution);
     
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
       for (s=0; s < Resol; s++)
       {
          Nl = (Nl+1)/2;
          Nc = (Nc+1)/2;
       }
       Nbr_Plan -= Resol;
    }     
    if (Verbose == True) printf ("\n Create an image of size (%d, %d)\n", Nl, Nc);
    Imag.resize (Nl, Nc);
    int NbrBand = (Nbr_Plan-1)*3 + 1;

//     FilterAnaSynt SelectFilter(F_MALLAT_7_9);
//     SubBandFilter SB1D(SelectFilter, NORM_L2);
//     SB1D.Border = I_MIRROR;
//     Ortho_2D_WT WT2D(SB1D);
    type_lift LiftingTrans=TL_F79;
    Lifting Clift1D(LiftingTrans);
    Clift1D.Border = I_MIRROR;
    Ortho_2D_WT WT2D(Clift1D);
    int *TabBandNl = new int [NbrBand];
    int *TabBandNc = new int [NbrBand];
    int *TabResolNl = new int [NbrBand];
    int *TabResolNc = new int [NbrBand];
    get_info_band(Nl,Nc,NbrBand,TabBandNl,TabBandNc,TabResolNl,TabResolNc);

    //extern long size_enc;
    float *LevelQ = new float[NbrBand];
    for (s =  NbrBand -1 ; s  >= 0; s--) 
    {
       int NbrCoef;
       int Ind=0, Depi=0, Depj=0;
       Nl_s = TabBandNl[s];
       Nc_s =  TabBandNc[s];
       NbrCoef = Nl_s*Nc_s;
       int *data;
       float Level;
       int Nls,Ncs;
       decode (FileDes, &data, &Nls, &Ncs, &Level);
       LevelQ[s] = Level;
       //int SizeScale = size_enc;
       // if (Verbose == True)
       //      cerr << "Set " << s+1 << " " << Nls << "x" << Ncs << 
       //          " Quantif = " << Level  <<  " Size = " << SizeScale << endl;
       if (s != NbrBand -1)
       {
          details Detail;
          int Scale;
          band2scale(s,TO_MALLAT, NbrBand, Scale, Detail);
          switch (Detail)
          {
             case D_HORIZONTAL:
               Depj += TabResolNc[s];
               break;
             case D_DIAGONAL:  
               Depi += TabResolNl[s];
               Depj += TabResolNc[s];
               break;
             case D_VERTICAL:  
               Depi += TabResolNl[s];
               break;
             case I_SMOOTH: break;
             case D_NULL:
             default: cerr << "Error: unknown detail" << endl;
                 exit(0);
                 break;
           }         
       }
       if (Level > FLOAT_EPSILON)
       {
           if (NbrBand > 1)
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++)
                  Imag(i+Depi,j+Depj) = ((float) data[Ind++]) * Level;
           }
           else
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++) Imag(i,j) = ((float) data[Ind++]) * Level;
           }
       }
       else
       {
           if (NbrBand > 1)
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++) Imag(i+Depi,j+Depj) =  data[0];
           }
           else  
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++)  Imag(i,j) = data[0];
           }
       }
       i_free(data);
    }        

   if (NbrBand > 1) WT2D.recons(Imag, Nbr_Plan);
 
   // In case of sub resolution reconstruction, the reconstructued
   // image must be renormalized because we used a L2 normalization
   float NormCoef=1.;
   if (Resol > 0)
   {
      for (s=0; s < Resol; s++) NormCoef *=2;
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) Imag(i,j) /= NormCoef;
  
   }     

   delete [] LevelQ;
   delete [] TabBandNl;
   delete [] TabBandNc;
   delete [] TabResolNl;
   delete [] TabResolNc;
}
    

/****************************************************************************/
