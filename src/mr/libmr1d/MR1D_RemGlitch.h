/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  MR1D_RemGlitch.h
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------    
**
**    this module contains the routines for temporal deglitching and filtering
**
**************************************************************************/ 

#ifndef __MR1D_DEGLITCH__
#define __MR1D_DEGLITCH__

#define DEF_NSCALE_FILTER  4
#define DEF_NSCALE_DEGLITCH  3
#define DEF_ITER_DEGLITCH  1
#define DEF_ITER_FILTER  5

float noise_tmp_cube(fltarray & Data, int Nit=3);
float noise_tmp_vect(fltarray & Data, int Nit=3);
void noise_tmp_vect(fltarray & Data, float Sigma, float &Mean, int Nit=3);

void filter_temp_cube(fltarray & Data, float N_Sigma, fltarray & ImaSigma_Noise, int Nbr_Scale, type_trans_1d Transform, float Epsilon, int Max_Iter);

void deglitch_temp_cube(fltarray & Data, fltarray & ImaSigma_Noise, float N_Sigma=3.,int Max_Iter=5,int Nbr_Scale=3);

void mr1d_deglitch(fltarray &Data, fltarray &Ima_SigmaNoise,
                fltarray & Ima_MeanClip,
                Bool Deglitch, Bool Filter,  
                int NbrScale=DEF_NSCALE_DEGLITCH, 
                int NscaleFilter=DEF_NSCALE_FILTER, 
                float N_Sigma=3., int MaxIter_Deglitch=DEF_ITER_DEGLITCH,
                int MaxIter=DEF_ITER_FILTER, float Epsilon=3);

#endif
