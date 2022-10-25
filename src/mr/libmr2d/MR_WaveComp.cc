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
**    File:  MR_WaveComp.cc
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
#include "SB_Filter.h"

/****************************************************************************/
 
void mr_i_ortho_compress (Iint &Imag, int Nbr_Plan, FILE *FileDes, Bool Verbose)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc(); 
    int s,i,j,Nls,Ncs,ind_i,ind_j,Nl1,Nc1;   
    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 
    int *PlanInt = i_alloc(Nl*Nc);
    float Level=1.;
    
    
    // encode the last smooth array 
    s = Nbr_Plan-1;
    Nls = size_ima_resol(s,Nl);
    Ncs = size_ima_resol(s,Nc);
    for (i=0;i < Nls; i++)
    for (j=0;j < Ncs; j++)  
    {
       ind_orthog_transf(s, i, j, I_SMOOTH, Nl,  Nc,  ind_i,  ind_j);
       PlanInt[i*Ncs+j]  =  Imag(ind_i,ind_j);
    }
    encode(FileDes, PlanInt,Nls,Ncs,1.); 
   
    // encode the scales
    for (s = Nbr_Plan-2; s >= 0; s--)
    {
       // image size at resolution s
       Nl1 = size_ima_resol(s,Nl);
       Nc1 = size_ima_resol(s,Nc);
//cout << "Nl1 = " << Nl1 << " Nc1 = " << Nc1 << endl;

       // diagonal band      
       Nls = Nl1/2;
       Ncs = Nc1/2;
       for (i=0;i < Nls; i++)
       for (j=0;j < Ncs; j++) 
       {
          ind_orthog_transf(s, i, j,   D_DIAGONAL, Nl,  Nc,  ind_i,  ind_j); 
          PlanInt[i*Ncs+j]= Imag(ind_i,ind_j);
       }
       encode(FileDes, PlanInt, Nls, Ncs, Level); 
 
        // vertical band
       Nls = Nl1/2;
       Ncs = (Nc1+1) /2;       
       for (i=0;i < Nls; i++)
       for (j=0;j < Ncs; j++) 
       {
          ind_orthog_transf(s, i, j, D_VERTICAL, Nl,  Nc,  ind_i,  ind_j); 
          PlanInt[i*Ncs+j]= Imag(ind_i,ind_j);
       } 
       encode(FileDes, PlanInt, Nls, Ncs, Level);       
      
       // horizontal band   
       Nls = (Nl1+1)/2;
       Ncs = Nc1 /2;
       for (i=0;i < Nls; i++)
       for (j=0;j < Ncs; j++)  
       { 
          ind_orthog_transf(s, i, j,  D_HORIZONTAL, Nl,  Nc,  ind_i,  ind_j); 
          PlanInt[i*Ncs+j]= Imag(ind_i,ind_j);
       }
       encode(FileDes, PlanInt, Nls, Ncs, Level);
 
       if (Verbose == True)
       {
                 cout << "Encode: Scale " << s << " Nl = " << Nls << 
                          " Nc = " << Ncs  << endl;
       }
    }     
    i_free(PlanInt);
}


/****************************************************************************/
 
void mr_i_ortho_decompress(Iint &Imag, int Resol, int & Nbr_Plan, FILE *FileDes, int &Nli, int &Nci, Bool Verbose)
{
    int i,j,s, *Plan, ind_i,  ind_j;
    int Nls, Ncs;
    float Level;
    int Nl, Nc;
 
    /* read the original image size */
    Nl =readint(FileDes); /* x size of image */
    Nc =readint(FileDes); /* y size of image */
    Nli = Nl;
    Nci = Nc;
    
    /* read the number of scales */
    Nbr_Plan=readint(FileDes); 

    if (Verbose == True)
    {
       printf ("DECODING wavelet transform: nbr_scale = %d\n", Nbr_Plan);
       printf ("Original image: Nl = %d, Nc = %d\n", Nl, Nc);
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
       Nl = size_ima_resol(Resol,Nl);
       Nc = size_ima_resol(Resol,Nc);
       Nbr_Plan -= Resol;
    }     
    if (Verbose == True) printf ("\n Create an image of size (%d, %d)\n", Nl, Nc);
    Imag.alloc (Nl, Nc, "decomp image");
      
    s = Nbr_Plan-1;
    decode(FileDes, &Plan, &Nls, &Ncs, &Level);
    for (i=0;i < Nls; i++)
    for (j=0;j < Ncs; j++)  
    {      
        ind_orthog_transf(s, i, j, I_SMOOTH, Nl,  Nc,  ind_i,  ind_j);
        Imag(ind_i,ind_j) = Plan[i*Ncs+j];
    }
    i_free(Plan);
    
    if (Verbose == True)
    {
        cout << "Decode: Scale " << s << " Nl = " << Nls << 
                          " Nc = " << Ncs <<  " Level = " << Level << endl;
    }         
	    
    for (s = Nbr_Plan-2; s >= 0; s--)
    { 
        decode(FileDes, &Plan, &Nls, &Ncs, &Level);
        for (i=0;i < Nls; i++)
        for (j=0;j < Ncs; j++)
        {
           ind_orthog_transf(s, i, j,  D_DIAGONAL , Nl,  Nc,  ind_i,  ind_j); 
           Imag(ind_i,ind_j) = Plan[i*Ncs+j];
        }
        i_free(Plan);
        if (Verbose == True)
        {
            cout << "Decode: Scale " << s << " Nl = " << Nls << 
                          " Nc = " << Ncs <<  " Level = " << Level << endl;
        } 
        decode(FileDes, &Plan, &Nls, &Ncs, &Level);
        for (i=0;i < Nls; i++)
        for (j=0;j < Ncs; j++) 
        {
           ind_orthog_transf(s, i, j, D_VERTICAL, Nl,  Nc,  ind_i,  ind_j); 
           Imag(ind_i,ind_j) = Plan[i*Ncs+j];
        } 
        i_free(Plan);
        decode(FileDes, &Plan, &Nls, &Ncs, &Level);
        for (i=0;i < Nls; i++)
        for (j=0;j < Ncs; j++) 
        {
           ind_orthog_transf(s, i, j, D_HORIZONTAL, Nl,  Nc,  ind_i,  ind_j); 
           Imag(ind_i,ind_j) = Plan[i*Ncs+j];
        } 
        i_free(Plan);
    }             
}
