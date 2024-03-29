/******************************************************************************
**                   Copyright (C) 2014 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Fred Ngole
**
**    Date:  24/09/2014
**    
**    File:  sr_util.h
**
*******************************************************************************
**
**    DESCRIPTION   
**    ----------- 
**    Tools used in sprite.cc and useful beyond the scope of super-resolution             
**
**
******************************************************************************/
#ifndef SR_UTIL_H
#define SR_UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include "TempArray.h"
#include <cmath>
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Sigma.h"
#include "CoaddCorrel.h"
#include "IM_Deconv.h"
#include "OptMedian.h"
//#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


const double Pi = 3.14159265358979323846;
const int LR_DENOISING = 0;
const int MR_SUPP_DENOISING=1;
const int ONE_STEP_ANALYSIS=1;
const int ACCURATE_ANALYSIS=0;
const int MODEL_SINC_INTERP=0;
const int MODEL_NON_LIN_INTERP=1;
const int MODEL_GAUSS_FIT=2;
const int FISTA = 1;
const int GEN_FWBW = 0;
//gsl_rng *rng;
//const gsl_rng_type * T;
//T =  gsl_rng_default;
//rng = gsl_rng_alloc (T);

typedef struct{
  type_transform transf;
  int nb_sc;
  int nb_usc;
  sb_type_norm Norm;
  type_undec_filter U_Filter;
  type_sb_filter SB_Filter;
  type_border Bord;
  FilterAnaSynt *FAS ;
}mr_opt;
int compute_centroid(Ifloat *img,double sig,double *cent,int niter=10);
int compute_centroid_arr(fltarray *data,double sig,double **cent,int niter=10);
double gauss2d(double amp,double x,double y,double theta,double sigx,double sigy);
void gauss2darray(double amp,double theta,double sigx,double sigy, double xcen,double ycen,Ifloat &im_out);
void gauss2D_fitting_first_guess(Ifloat im,double*sigx,double*sigy,double*xcen,double*ycen,double*amp,int smooth_ker);
void gauss2d_fit(Ifloat im,double *xcen,double *ycen,double* amp,double sig0,double *sig,int nb_iter_alter=20,int nb_subiter=10,double *mse=NULL,int max_test=10);
void gauss2d_fit_2(Ifloat im,double *xcen,double *ycen,double *amp,double *sigx,double *sigy,int nb_iter_alter=20,int nb_subiter=10,double *mse=NULL,int max_test=10);
double gauss2d_der1sig(Ifloat res,Ifloat gauss_approx,Ifloat distances_map,double sig);
double gauss2d_der1sigx(Ifloat res, Ifloat gauss_approx, double sigx, double xcen);
double gauss2d_der1sigy(Ifloat res, Ifloat gauss_approx, double sigy, double ycen);
double gauss2d_der2sig(Ifloat res,Ifloat gauss_approx,Ifloat distances_map,double sig);
double gauss2d_der2sigx(Ifloat res, Ifloat gauss_approx, double sigx, double xcen);
double gauss2d_der2sigy(Ifloat res, Ifloat gauss_approx, double sigy, double ycen);
double gauss2d_der1xcen(Ifloat res, Ifloat gauss_approx, double sig, double xcen);
double gauss2d_der1ycen(Ifloat res, Ifloat gauss_approx, double sig, double ycen);
double gauss2d_der2xcen(Ifloat res, Ifloat gauss_approx, double sig,double xcen);
double gauss2d_der2ycen(Ifloat res, Ifloat gauss_approx, double sig,double ycen);
double gauss2d_der2xcenycen(Ifloat res, Ifloat gauss_approx, double sig, double xcen, double ycen);
void dist_map(Ifloat &map,double xcen,double ycen);
void thresholding(Ifloat *data_in,Ifloat *data_out,double thresh,bool thresh_type=0);
void thresholding(Ifloat *data,Ifloat *data_out,Ifloat *thresh,bool thresh_type=0);
void thresholding(fltarray *data_in,fltarray *data_out,fltarray *thresh,bool thresh_type=0);
void circ_thresh(Ifloat *data_in,Ifloat *data_out,double r,double *cen);
void circ_thresh(fltarray *data_in,fltarray *data_out,double r,double **cen);
int sign_num (double a);
int wl_trans(Ifloat *Dat,mr_opt opt,MultiResol* MR_Data);
void wl_thresholding(MultiResol*wav_coeff,fltarray *thresh,bool thresh_type=0);
void wl_filter(Ifloat *img,Ifloat*img_filt,mr_opt opt,double nsig=4,bool thresh_type=0,fltarray* noise_map=NULL,fltarray *coeff_thresh=NULL,Ifloat *coarse_scale=NULL);
void wl_filter_analysis(Ifloat img,Ifloat &img_filt,mr_opt opt,double nsig,bool thresh_type,fltarray std_map,MultiResol &coeff_init,int nb_iter=100,double mu=0.5,double *mse=NULL);
int wl_gnoise_est(MultiResol* wav_coeff,fltarray* noise_arr);
mr_opt mr_opt_init(int Nx,type_transform transf = TO_UNDECIMATED_MALLAT,int nb_sc=-1,int nb_usc=-1,sb_type_norm Norm=NORM_L1,type_undec_filter U_Filter=DEF_UNDER_FILTER,type_sb_filter SB_Filter=F_MALLAT_7_9,type_border Bord=I_CONT,Bool Verbose =False);
int noise_est(fltarray* data,fltarray* noise_arr);
void wl_arr_filter(fltarray *data,fltarray*data_filt,mr_opt opt,double nsig=4,bool thresh_type=0);
void coadd_frames(fltarray Data, Ifloat &Ima,Bool GetMedian=False);
void coadd_frames(fltarray Dat, Ifloat &Ima,float Zoom,
		      Bool GetMedian=True,int MaxDist=4, 
		      int Surface=10, type_interp_corr TypeInterp=ICF_TANH,Bool OptNoOffset=False,Bool MeanSub=False,Bool Verbose=False);
void flux_est(Ifloat *data, double r,double *cent,double*flux);
void flux_est(fltarray *data, double r,double **cent,double*flux,int ref_im_ind=0);
void shift_est(fltarray *data, double sig_gfit, double *sig_vec,double **cent,double **shift,int ref_im_ind=0,double nsig=4);
double sinc(double x);
void lanczos(double u,double*mask,int mask_rad=10);
void lanczos(double u[],double**mask,int mask_rad=10);
void lanczos_upsamp(Ifloat &im,int upfact,Ifloat &im_out, double *u=NULL,int kern_rad=5);
void lanczos_stacking(fltarray &data,Ifloat &im_upsamp,int upfact,double **u,int kern_rad=5,Bool mean_en=False);
void im_stacking(fltarray &data,Ifloat &im_stack,Bool mean_en=False);
void decim(Ifloat *img,Ifloat *img_filt,Ifloat *img_dec,int D,Ifloat* mask,Bool filt_en=True,Bool fft_en=True ); 
void rotate(Ifloat *input,Ifloat *output,int dir);
void transpose_decim(Ifloat*im_in,Ifloat*im_out,int D);
void power_meth(void(*opname)(double*,double*),int size_input,int *nb_iter,double *spec_rad,int nb_iter_max=100,double tol=0.01);
double norm2 (double * x,int siz);
void scale (double * x,double siz,double a);
double dot(Ifloat Im1,Ifloat Im2);
void hadamard(Ifloat Im1,Ifloat Im2,Ifloat &Im_out);
void ineq_cons(Ifloat &A, Ifloat B,double b=1);
void reverse(double*u,int length,double*ur);
void holes(double*h,int sc_nb,int length,double*hout);
void convol1D(double*in,double*out,double*h,int length,int filt_length);
void convol1D_cent(double*in,double*out,double*h,int length,int filt_length);
void correl2D(Ifloat im1,Ifloat im2, Ifloat &im_correl);
void sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg);
void transp_sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg);
void grad(Ifloat im,Ifloat &im_gradx,Ifloat &im_grady,double* hx,double* hy,int lhx,int lhy);
void grad_cons_noise_est(Ifloat im_in,Ifloat& im_out,Ifloat &res,double * sig,double gamma,int nb_iter=5,double * mse_iter=NULL,double mu=0.5, double *hx=NULL,double *hy=NULL,int lhx=-1,int lhy=-1);
void mad_res_noise_est(Ifloat im_in,Ifloat& im_out,Ifloat &res,double * sig,int nb_iter=5,double * mse_iter=NULL,double*sig_iter=NULL,double mu=1, int nsig=5,Bool sig_up=True,int thresh_type=0);
//void randomn(double *x,double sig=1,int length=1,double mean=0);
void randomngsl(double *x,double sig,int length,double mean=0);
//void randomn(Ifloat *mat,double sig=1,double mean=0);
void randomngsl(Ifloat *mat,double sig=1,double mean=0);
//void randomn(fltarray *mat,double *sig,double mean=0);
void randomngsl(fltarray *mat,double *sig,double mean=0);
//void randomn(fltarray *mat,double *sig,double*mean);
void randomngsl(fltarray *mat,double *sig,double*mean);
double var_iter(double x,int n,double *mean,double*M);
void check_ineq(Ifloat dat,double thresh,Bool abs_en,Ifloat &flag);
void mr_support_filt(Ifloat img,Ifloat &img_filt, mr_opt opt,double nsig=5,int nb_iter=1,double*mse=NULL,Bool Pos_coeff=True,double lambda=1,Bool coarse_cons=False,Bool pos_cons=True, Bool drop_coarse=True,Bool iso_cons=False, double sig=89.108911);

inline void mod_denoise_usage(void)
{
    fprintf(OUTMAN, "         [-m type_of_model_denoising]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Wavelet denoising of LR images before coaddition");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Iterative multiresolution support denoising after coaddition");
    fprintf(OUTMAN, "             default is %s.\n",  "Iterative multiresolution support denoising after coaddition");
}

inline void mod_comp_usage(void)
{
    fprintf(OUTMAN, "         [-M model_computing_method]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Sinc Interpolation");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Non linear interpolation");
    fprintf(OUTMAN, "              %d: %s \n",2,
                                            "Gaussian fitting after sinc interpolation (no denoising)");
    fprintf(OUTMAN, "             default is %s.\n",  "Sinc Interpolation");
}

inline void min_id_usage(void)
{
    fprintf(OUTMAN, "         [-A Minimisation algorithm]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Generalized forward-backard (analysis+positivity)");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "FISTA");
    fprintf(OUTMAN, "             default is %s.\n",  "Generalized forward-backard (analysis+positivity)");
}

inline void sparse_const_type_usage(void)
{
    fprintf(OUTMAN, "         [-D sparse_constraint_implementation (if FISTA is chosen for the minimisation)]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Accurate analysis constraint with subiteration");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "One step approximation of analysis constraint implementation");
    fprintf(OUTMAN, "             default is %s.\n",  "Accurate analysis constraint with subiteration");
}

inline void noise_est_usage(void)
{
    fprintf(OUTMAN, "         [-e Noise_est_meth]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Estimation on backprojected residual using a MAD");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Noise simulation");
    //fprintf(OUTMAN, "              %d: %s \n",2,
    //                                      "Direct calculation of noise standard deviation per pixel (option under development)");
	    
    fprintf(OUTMAN, "             default is %s.\n",  "Estimation on backprojected residual using a MAD");
}

inline void thresh_type_usage(void)
{
    fprintf(OUTMAN, "         [-c Thresh_type]\n");
    
    fprintf(OUTMAN, "              %d: %s \n",0,
                                            "Hard Thresholding");
    fprintf(OUTMAN, "              %d: %s \n",1,
                                            "Soft Thresholding");
    //fprintf(OUTMAN, "              %d: %s \n",2,
    //                                      "Direct calculation of noise standard deviation per pixel (option under development)");
	    
    fprintf(OUTMAN, "             default is %s.\n",  "Soft Thresholding");
}

void decim_conv_mat(Ifloat conv_mat,int im_siz[],int decim_fact,Ifloat &output,double flux=1,double sig=1);
void convmask_shift_extend(Ifloat mask,int im_siz[],int shift[],Ifloat &output);
void noise_map_comp(double **shifts, double *sig, double *flux,int nb_im,int lancrad, int im_siz[],mr_opt opt, int decim_fact,fltarray &output);
int renewSeed();
void mr_norm(mr_opt opt,double*nrm,int nb_band,int Nx,int Ny=-1);
void wl_gnoise_est_dat(Ifloat data,mr_opt opt,fltarray* noise_arr);
void shell(unsigned long n,double *a);
double median(double *arr,int n);
#endif


