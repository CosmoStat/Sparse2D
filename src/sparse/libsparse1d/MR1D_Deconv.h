/******************************************************************************
**                   Copyright (C) 2003 by IRSN
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Philippe Querre
**
**    Date:  03/03/03
**
**    File:  MR_Deconv.h
**
*******************************************************************************
**
**    DESCRIPTION  Definition for multiresolution deconvolution
**    -----------
**
**    PARAMETRES
**    ----------
**
**    RESULTS
**    -------
**
**
******************************************************************************/

#ifndef __MR1D_DECONV__
#define __MR1D_DECONV__


#include "MR_Sigma.h"
#include "MR1D_Obj.h"
#include "MR1D_NoiseModel.h"
#include "FFTN_1D.h"
#include "MR1D_Regul.h"
#include "MR1D_Filter.h"


enum type_init_deconv { DEC_INIT_ZERO,
                        DEC_INIT_FLAT,
                        DEC_INIT_IMA,
                        DEC_INIT_GUESS};

class MR1D_Deconv {

protected:

   FFTN_1D FFT1D;

   // compute the residual
   virtual void compute_resi();

   // float find_optim(Ifloat &Gradient, Ifloat & Buff);
   void sig_grad(fltarray &gradient);

   void init_param();

   int Nx;            // working signal size
   Bool UseICF;       // if True, the ICF is given by the user
   fltarray Sig_ICF;  // ICF signal
   Bool KeepImagn;    // Keep the variable Signal_n
   Bool WaveFilterResi; // if true the residual is filtered in the wavelet
                        // space at each iteration using the multiresolution
			// support
   Bool ApplySoftRegul;      // Apply a soft thresholding at each iteration
   fltarray Mult;       // Multiplyer for MEM methods
   float  FluxImag;     // Image flux
   type_init_deconv TypeInit; // Initialisation type
   Bool AdjointRec;     // Use the adjoint for the reconstruction
   cfarray *TabCB_cf; // CLEAN BEAM WT in Fourier space (for MRC only)

   // optimization routine
   float find_optim_tikhonov (fltarray &Gradient);
   float find_optim_xi2 (fltarray &Gradient);
   float find_optim_poisson (fltarray &Gradient);
   float im_find_optim (fltarray &Gradient);

   // compute fonctional value
   float fonctional();

   // iterative deconvolution (init_deconv must have been called before)
   void sig_iter_deconv();

//   void init_mrc()      ; // initialization for MultiResolution CLEAN
//   void mrc(fltarray *ICF); // MultiResolution CLEAN deconvolution
//   void mrc_rec_obj(MultiResol &MR_SOL, cmparray *TabCB_cf);
                         // solution reconstruction in Obj
			 // using the dirac list in MR_SOL, and CLEAN BEAM
			 // table TabCB_cf
  void  vaguelet();
  void  va_inverse_data(fltarray &NoiseSimu);

public:


   fltarray Obj;      // solution
   fltarray Resi;     // residual
   fltarray Signal;   // input data
   fltarray Signal_n; // only used by Lucy method
   fltarray Psf;      // Point spread function
   cfarray Psf_cf;    // PSF Fourier transform

   Bool PositivConstraint; // if true, the solution must positive
   Bool NormFlux;     // if True, the flux is renormalized at each iteration
   float IterCvg;     // congergence paramter applied to the gradient
                      // at each iteration
   Bool OptimParam;   // Optimize the convergence parameter at each iteration
                      // if OptimParam == True then IterCvg is mot used
   float EpsCvg;      // congergence parameter for the final convergence
   type_deconv DecMethod; // type of deconvolution method

   float RegulParam;  // User regularization parameter
   int MaxIter;       // maximum number of iterations
   int MaxMRCleanIter;// maximum number of iterations for MR CLEAN method
   Bool UseMRCEnergy; // Use standard clean in MRC if UseMRCEnergy = False
                      // otherwise use the ENERGY map in MRC
   Bool CleanLastScale; // if True, CLEAN is also applied on the last scale
                        // otherwise apply a simple Van Citter method
   int CleanFirstScale;// First Scale to CLEAN
   Bool SupportDilate; // Dilate the multiresolution support
                       // (used by MRC method)

   Bool Verbose;      // verbose mode
   Bool PsfMaxShift;  // by default the PSF is shifted to place its max
                      // at the center.
		      // if PsfMaxShift == False, this is not done
   Bool GaussConv;    // If true, the solution has a limited resolution
                      // (i.e. it must a convolution product between
		      // a gaussian (fixed by Fwhm) and a hidden solution)
   float Fwhm;        // Full Width at Half Maximum of the Gaussian ICF
   type_noise StatNoise; // type of noise
   float Noise_Sig;   // Gaussian noise standard deviation
   float N_Sigma;     // CLEAN parameter: stop at N_Sigma above the noise
   fltarray MemModel; // Model image for mem restoration using Gulling entropy
   Bool UseModel;     // True when a model is used with MEM MODEL method
   type_border Border;// Border type used for multiscale transform calculation

   // Parameter for wavelet methods
   Bool KillLastScale;// Do not restore the last scale
   MR_1D MR1D_Data;
   MR1DNoiseModel *ModelData; // noise model class

   void convol_gauss(fltarray &Data, float FWHM);
   // convolve (using the FFT an image by a Gaussian of size FWHM

   MR1D_Deconv()  {init_param(); }

   // initialize the deconvolution
   void init_deconv(fltarray *FirstGuess=NULL, fltarray *ICF=NULL);

   // deconvolution
   void sig_deconv (fltarray *FirstGuess, fltarray *ICF);

   RegulSig RSig; // Regularization Class
   type_1d_filter Filter; // Type of filtering used in wavelet waguelet decomposition

   virtual ~MR1D_Deconv(){};
   };



void dec_clean (fltarray& Signal, fltarray& Psf, fltarray& Signal_Out,
                float Fwhm, float Noise, float Gamma=DEFAULT_CLEAN_GAMMA,
                int Niter=DEFAULT_CLEAN_NITER, int K_Aver=DEFAULT_CLEAN_KAVER,
                Bool AddResi=True, Bool SearchPositiv=True,
		Bool UseEnergy=False, Bool Verbose=False);
#endif
