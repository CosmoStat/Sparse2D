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
**    File:  sprite.h
**
*******************************************************************************
**
**    DESCRIPTION   
**    ----------- 
**    Sprite class definition             
**
**
******************************************************************************/

#ifndef SPRITE_H
#define SPRITE_H
#include <stdio.h>
#include <stdlib.h>
#include "TempArray.h"
#include <cmath>
#include "sr_util.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "MR_Sigma.h"
#include "CoaddCorrel.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



class Sprite
{
    public:
    
  Sprite(int upfact,fltarray *data,mr_opt mod_opt,mr_opt sr_opt,int nb_iter,int mod_denoising_type=MR_SUPP_DENOISING,double nsig=3,string output_dir="",string output_file_name="solution.fits",string res_file_name="resi.fits",int noise_est_meth=0,int sparse_cons_type=ACCURATE_ANALYSIS, int min_id=GEN_FWBW,Bool eq_flux=True,Bool eq_noise=True,int win=50,double tol=0.05,bool Verbose=False,char * shift_filename=NULL,Bool Sparse_en=True,char * dist_weights_filename=NULL,double spectral_rad=0,int model_type=MODEL_SINC_INTERP,Bool fft_en=True,int ref_im_ind=0,Bool pos_cons=True,Bool deconv_en=False,Ifloat*conv_mask=NULL,int lanc_rad=10,double w1=0.5,double lambdat=1.4);
  void mod_gen(int ref_im_ind,int denoising_type,float zoom=-1,double nsig=5,bool thresh_type=0);
  void mod_gen_2(int ref_im_ind,int denoising_type, int nsig=5,bool thresh_type=0);
  void mod_gen_gauss_fit(void);
  double lr_psf_size(void);
  Ifloat get_mod();
  double** get_cent();
  double** get_shift();
  double* get_flux();
  void set_spectral_rad(double spectral_rad);
  void sprite_op_mat(double *input=NULL,double *output=NULL);
  void sprite_op_mat(Ifloat input,fltarray &output); // For unitary tests
  void sprite_transp_op_mat(double *input=NULL,double *output=NULL);
  void sprite_transp_op_mat(fltarray input,Ifloat &output); // For unitary tests
  void sprite_global_op_mat(double *input,double *output);
  void spectral_rad_est(int *nb_iter,int nb_iter_max=100,double tol=0.001);
  void min1(Ifloat &delta,Bool thresh_type=True,int nb_rw=1);
  void min2(Ifloat &delta,Ifloat &z1,Ifloat &z2,Bool thresh_type=True,int nb_rw=1); 
  void write_outputs(Bool writeRes,Bool writeParam);
  void do_it(Bool wr_en,Bool writeRes,Bool writeParam,int noise_sim_iter=-1,Bool thresh_type=True,int nb_rw=2);
  void noise_simulator(int nb_monte);
  Ifloat get_sol();
  ~Sprite();
 private:
		 Bool deconv_en;
		 Bool fft_en;
		 Bool pos_cons;
		 Bool eq_flux;
		 Bool eq_noise;
		 bool Verbose;
		 Bool Sparse_en;
                 int Nx;
		 int Ny;
		 int Nz;
		 int ref_im_ind;
		 int lanc_rad;
		 int upfact;
		 int nb_iter;
		 int min_id; // Indicates the minimizer to use
		 int noise_est_meth; // Noise estimation method ; 0 = direct estimation on transformed gradient, 1 = noise simulation (the implemented methods assume a gaussian noise)
		 int cv_win;
		 int mod_denoising_type;
		 int nb_iter_cv;
		 int sparse_cons_type;
		 int model_type;
		 double spectral_rad;
		 double cv_tol;
		 double lambdat; // Generalized fw-bw parameter
		 double w1; // Generalized fw-bw parameter
		 double * iter_diff_buff; 
		 double * flux;
		 double * dist_weights;
		 double * nsig;
		 double * mse_steep;
		 double ** cent;
		 double ** shift;
		 double ** shift_cent;
		 double sig_max;
		 string output_dir;
		 string output_file_name; 
		 string resi_file_name; 
		 mr_opt mod_opt;
		 mr_opt sr_opt;
		 Ifloat model;
		 Ifloat sr_rec; // The solution
		 Ifloat gradient;
		 Ifloat noise_est_trace ;
		 Ifloat coarse_scale;
		 Ifloat * conv_mask;
		 fltarray *data;
		 fltarray resi;
		 fltarray resi_out;
		 fltarray sr_rec_lr;
		 fltarray op_mat_output;
		 fltarray noise_sig_map;
		 fltarray lanczos_ker;
		 fltarray data_nrm;
		 fltarray wvl_noise_map;
				 
};

#endif
