/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  CurContrast.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Curvelet contrast Enhancement
**    -----------  
**                 
******************************************************************************/
 
#include "Curvelet.h"
#include "IM_IO.h"
#include "CurContrast.h"

/****************************************************************************/
//
//                     Contrast enhancement
// 
/****************************************************************************/

void CurContrast::cur_enhance(Ifloat &Data, Curvelet & Cur, Bool PosIma)
{
   int i,j;
   fltarray Trans;
   double Q100 = pow((double) 100., Contrast_Q_Param);
   float Max=0.;
   // float TabMin, TabMax;
   // float TabTMin, TabTMax;    
      
   if (L_Coeff != 0)  Max = max(Data);
          
    if (Contrast_Q_Param != 0)
      for (i=0; i < Data.nl(); i++)
      for (j=0; j < Data.nc(); j++)
  	 Data(i,j) = (float) (pow((double) Data(i,j), 1.-Contrast_Q_Param)*Q100);

    // TabMin = min(Data);
    // TabMax = max(Data);
             
    //float N_Sigma=3.;
    //Bool FilterData = False;
    //  
    //Cur.enhance(Data, Data, Noise_Ima,    
    //            N_Sigma,   FilterData, 
    //	        Contrast_P_Param,Contrast_Q_Param, Contrast_M_Param);
 		      
    Cur.transform(Data,Trans);
    // fits_write_fltarr("xx_tr.fits", Trans);
    for (int s2d=0; s2d < Cur.NbrScale2D-1; s2d++)
    for (int s1d=0; s1d < Cur.nbr_rid_scale(s2d); s1d++)
    {
        float Norm = Cur.norm_band(s2d, s1d)*sqrt((float) Cur.TabBlockSize(s2d));
	int i,j,Posi = Cur.ipos(s2d, s1d);
	int Posj = Cur.jpos(s2d, s1d);
	int Nl = Cur.size_nl(s2d, s1d);
	int Nc = Cur.size_nc(s2d, s1d);
	float Noise = Noise_Ima*Norm;
 	Contrast_C_Param = Noise*NSigmaLow;
	float MaxCoef = 0.;

 	if ((Verbose == True) || (UseSigmaUp == False))
	{
 	   for (i = 0; i < Nl; i++)
	   for (j = 0; j < Nc; j++) 
	      if (MaxCoef < ABS(Trans(j+Posj, i+Posi, s2d)))
	                               MaxCoef = ABS(Trans(j+Posj, i+Posi, s2d));

        }
	if (UseSigmaUp == False) 
	     Contrast_M_Param = MaxCoef*Contrast_M_ParamCoef; //  / (float)(s2d+1);
	else Contrast_M_Param = NSigmaUp*Noise;
	// Contrast_M_Param = 100.;

        if (Verbose == True) 
	{
	  cout << "Scale " << s2d << " " << s1d << "  MaxCoef. = " << MaxCoef << endl;
	  cout  << " Noise = " << Noise << " Protect. Level (M_Param) = " << Contrast_M_Param << endl;
 	}
	for (i = 0; i < Nl; i++)
	for (j = 0; j < Nc; j++)
	{
	   float CoefCur,Gamma=1.;
	   CoefCur = Trans(j+Posj, i+Posi, s2d);
 	   if (UseVeldeFunction == True)
	        Gamma = cur_contrast_function_velde(ABS(CoefCur));
	   else Gamma = cur_contrast_function(ABS(CoefCur));
	   // Gamma /= Norm;
// 	   if (Gamma < 1) 
// 	   {
// 	       cout << "PB .... " << endl;
// 	       cout << "Noise = " << Noise << " Contrast_C_Param " << Contrast_C_Param << endl;
//                cout << "Contrast_M_Param = " << Contrast_M_Param << " Contrast_P_Param " << Contrast_P_Param << endl;
//                cout << "Contrast_S_Param = " << Contrast_S_Param << " CoefCur " << CoefCur << endl;
// 	       cout << "Gamma = " << Gamma << endl;
// 	       exit(-1);
// 	   } 
   	   Trans(j+Posj, i+Posi, s2d) *= Gamma;
 	}
     }
 
    Cur.recons(Trans,Data, False);		
    if (PosIma == True) threshold(Data);
//     TabTMin = min(Data);
//     TabTMax = max(Data);
//     
//     for (i=0; i < Data.nl(); i++)
//     for (j=0; j < Data.nc(); j++)
//     {
//        float Scale = (TabMax - TabMin) / (TabTMax - TabTMin);
//        Data(i,j) =  (Data(i,j) - TabTMin) * Scale + TabMin;
//     }
    
     // inverse mapping of the Luminance
    if (Contrast_Q_Param != 0)
    {
       if (Verbose == True) cout << "inverse mapping of the Luminance " <<  endl;
       double  Q1 = 1. / (1-Contrast_Q_Param);
       Q100 = pow((double) 100., -Contrast_Q_Param*Q1);
       for (i=0; i < Data.nl(); i++)
       for (j=0; j < Data.nc(); j++)
 	  Data(i,j) = (float) (pow((double) Data(i,j), Q1)*Q100);
    }
    
    // Sature the enhanced image at the max*L_Coeff of the restore image
    if (L_Coeff != 0)
    {
       Max *= L_Coeff;
       for (i = 0; i < Data.nl(); i++)
       for (j = 0; j < Data.nc(); j++)
       {
           // Data(i,j) =  (Data(i,j) - TMin) * Scale + Min;
           // if (Data(i,j) < Min) Data(i,j) = Min;
           if (Data(i,j) > Max) Data(i,j) = Max;
       }
    }     
}

/*************************************************************************/

inline void  noise_RGB_to_YUV(float R, float G, float B, float & Y, float &  Cr, float & Cb)
{
   Y  =  sqrt(sqr((double)(0.257 * R)) 
             + sqr((double)(0.504 * G))
	     + sqr((double)(0.098 * B)));
   Cr =  sqrt(sqr((double)(0.439 * R)) +
              sqr((double)(0.368 * G)) +
	      sqr((double)(0.071 * B)));
   Cb =  sqrt(sqr((double)(0.148 * R)) +
              sqr( (0.291 * G)) +
              sqr((double) (0.439 * B)));
}

void CurContrast::curcol_enhance(fltarray &Data, Curvelet & Cur, Bool UseYUV)
{
   int i,j,k;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   int Nl = Ny;
   int Nc = Nx;
   float *PtrData = Data.buffer();
   fltarray *TabTrans = new fltarray[Nz];
   float *TabNoise = new float [Nz];
   fltarray TabMean(Nz), TabTMean(Nz);
   fltarray TabMin(Nz), TabMax(Nz);
   fltarray TabTMin(Nz), TabTMax(Nz);
   
       
    if (Verbose == True)  printf("RGB --> LUV transformation.\n"); 
    if (Nz == 3) 
       if (UseYUV == True) rgb_to_yuv(Data);
       else rgb_to_luv(Data);
    
    Ifloat *TabFrame = new Ifloat [Nz];
    // Noise_Ima = 0.;
    if (Noise_Ima > FLOAT_EPSILON)
       noise_RGB_to_YUV(Noise_Ima, Noise_Ima, Noise_Ima,
                        TabNoise[0], TabNoise[1],TabNoise[2]);
    for (k = 0; k < Nz; k++)
    {
       TabFrame[k].alloc(PtrData, Nl, Nc);
       PtrData += Nl*Nc;
       TabMin(k) = min(TabFrame[k]);
       TabMax(k) = max(TabFrame[k]);
       TabMean(k) = average(TabFrame[k]);
       if (Noise_Ima < FLOAT_EPSILON)
          TabNoise[k] = detect_noise_from_bspline(TabFrame[k]);
       if (Verbose == True) printf("Frame %2d: Sigma Noise = %f\n", k+1,TabNoise[k]);
       if (Verbose == True) printf("           Curvelet transform.\n");
       Cur.transform(TabFrame[k],TabTrans[k]);
    }    

    Ifloat Grad;
    if (Verbose == True) 
    {
        printf("Curvelet coeff. enhancement.\n");
	cout << "   MParam = " << Contrast_M_ParamCoef << endl;
	cout << "   PParam = " << Contrast_P_Param   << endl;
	cout << "   TParam = " << Contrast_S_Param   << endl;
 	if (Contrast_Q_Param != 0) cout << "QParam = " << Contrast_Q_Param   << endl;
   }
   	
    for (int s2d=0; s2d < Cur.NbrScale2D-1; s2d++)
    for (int s1d=0; s1d < Cur.nbr_rid_scale(s2d); s1d++)
    {
        float Norm = Cur.norm_band(s2d, s1d)*sqrt((float) Cur.TabBlockSize(s2d));
	int i,j,Posi = Cur.ipos(s2d, s1d);
	int Posj = Cur.jpos(s2d, s1d);
	int Nlb = Cur.size_nl(s2d, s1d);
	int Ncb = Cur.size_nc(s2d, s1d);
  	float Noise = 0.;
	for (k = 0; k < Nz; k++)
	{
	   Noise += pow((double) TabNoise[k]*Norm, (double)4.);
	}
        // Noise is the noise related to the value sqrt(sum_k frame_k^2)
	Noise = pow((double) Noise, (double)(1./4.));
	
 	Contrast_C_Param = Noise*NSigmaLow;
	float MaxCoef = 0.;
	Grad.resize(Nlb, Ncb);
         
        for (i = 0; i < Nlb; i++)
        for (j = 0; j < Ncb; j++) 
	{
	   float Val = 0.;
	   for (k = 0; k < Nz; k++)
	         Val += (TabTrans[k])(j+Posj, i+Posi, s2d) * (TabTrans[k])(j+Posj, i+Posi, s2d);
	   Val = sqrt(Val);
	   Grad(i,j) = Val; 
	   if (MaxCoef < Val) MaxCoef = Val;
        }
	if (UseSigmaUp == False) Contrast_M_Param = MaxCoef*Contrast_M_ParamCoef;
	else  Contrast_M_Param = NSigmaUp*Noise;
        if (Verbose == True) 
	{
	  cout << "Scale " << s2d << " " << s1d << "  MaxCoef. = " << MaxCoef << "  M_ParamCoef " << Contrast_M_ParamCoef << endl;
	  cout  << " Noise = " << Noise << " Protect. Level = " << Contrast_M_Param << endl;
 	}
	
	for (i = 0; i < Nlb; i++)
        for (j = 0; j < Ncb; j++) 
	for (k = 0; k < Nz; k++)
	  (TabTrans[k])(j+Posj, i+Posi, s2d) *= cur_contrast_function(Grad(i,j));
    }
    
    for (k = 0; k < Nz; k++)
    {
       TabFrame[k].init();
       if (Verbose == True) cout << "Reconstruction frame " << k+1 << endl;
        Cur.recons(TabTrans[k],TabFrame[k]);		
       TabTMin(k) = min(TabFrame[k]);
       TabTMax(k) = max(TabFrame[k]);
       float Scale = (TabMax(k) - TabMin(k)) / (TabTMax(k) - TabTMin(k));
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
            (TabFrame[k])(i,j) =  
	        ((TabFrame[k])(i,j) - TabTMin(k)) * Scale + TabMin(k);
       TabTMean(k) = average(TabFrame[k]);
       if (k > 0)
       for (i=0; i < Nl; i++)
       for (j=0; j < Nc; j++)
       {
           Data(j,i,k) += TabMean(k)  - TabTMean(k);
       }
    }
    
    if (Nz == 3) 
      if  (UseYUV == True) yuv_to_rgb(Data);
      else luv_to_rgb(Data);
                
//     for (k = 0; k < Nz; k++)
//     for (i = 0; i < Nl; i++)
//     for (j = 0; j < Nc; j++)
//         if (Data(j,i,k) < 0.) Data(j,i,k) = 0.;

    delete [] TabFrame;
    delete [] TabTrans;
    delete [] TabNoise;
}

/*************************************************************************/

void CurContrast::curcol_enhance_luminance(fltarray &Data, Curvelet & Cur, Bool UseYUV)
{
   int i,j,k;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   int Nl = Ny;
   int Nc = Nx;
   float *PtrData = Data.buffer();
   fltarray Trans;
   fltarray BuffData(Nc,Nl,Nz);
   BuffData = Data;
     
     
    if (UseYUV == True) rgb_to_yuv(Data);
    else rgb_to_luv(Data);
    
    Ifloat ImaFrame, Buff(Nl,Nc, "buff");
    ImaFrame.alloc(PtrData, Nl, Nc);
    Buff = ImaFrame;
    
    if (Noise_Ima < FLOAT_EPSILON)
    {
       Noise_Ima = detect_noise_from_bspline(Buff);
             // detect_noise_from_med (Buff);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    cur_enhance(Buff, Cur);
    
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
      if (ImaFrame(i,j) != 0) Buff(i,j) /= ImaFrame(i,j); 
	 
   for (k = 0; k < Nz; k++)
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
   {
      Data(j,i,k) = BuffData(j,i,k)*Buff(i,j);
      if (Data(j,i,k) < 0.) Data(j,i,k) = 0.;
   }
}

/*************************************************************************/
