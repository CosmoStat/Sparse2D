/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  03/03/03
**
**    File:  MR1D_Deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  signal deconvolution using multiresolution
**    -----------
**
**
***************************************************************************/

#include "MR1D_Deconv.h"
#include "DefFunc.h"



/**************************************************************************/
/**********                        ****************************************/
/**********     smooth_mediane     ****************************************/
/**********                        ****************************************/
/**************************************************************************/
void smooth_mediane (const fltarray& Sig1, fltarray& Sig2,
                     type_border Type, int Step_trou, int Window_Size) {

    int Nx = Sig1.nl();
    int k,i,ind_fen;
    int Step;
    float *fenetre;
    int  Window2 = (Window_Size - 1) / 2;
    int  Size = Window_Size;

    if (Step_trou > 0) Step = (int)(pow((double)2., (double) Step_trou) + 0.5);
    else Step = 1;

    fenetre = new float [Size];
    fenetre[0]=0.;
    for (i = 0; i < Nx; i++) {

       ind_fen = 0;
       for (k = i - Window2*Step; k <= i+Window2*Step; k += Step)
          fenetre[ind_fen++] = Sig1(k , Type);

       // Imag2(i,j) = hmedian(fenetre, ind_fen);
       if (Size == 9) Sig2(i) = opt_med9(fenetre);
       else Sig2(i) = get_median(fenetre, Size);
    }
    delete[] fenetre;
}


/**************************************************************************/
/**********                               *********************************/
/**********     detect_noise_from_med     *********************************/
/**********                               *********************************/
/**************************************************************************/
/* search for an estimation of sigma in noise with median transform       */
/**************************************************************************/
extern float BadPixalVal;
extern Bool BadPixel;
float detect_noise_from_med (const fltarray& Signal) {

    double Noise;
    int i;
    int N = Signal.nx();
    Bool Average_Non_Null = True;
    int Nit = 3;

    fltarray Buff(Signal.nx(), "Buff noise estimation");
    smooth_mediane (Signal, Buff, DEFAULT_BORDER, 0, 3);
    for (i=0; i < Signal.nx(); i++) Buff(i) = Signal(i) - Buff(i);
    // INFO(Buff, "diff");

    Noise=get_sigma_clip (Buff.buffer(), N, Nit,
                          Average_Non_Null, BadPixel, BadPixalVal);
    Noise /= 0.972463;
    // printf("Calc sigma(%d) = %f\n", N, (float) Noise);
    // if (Noise < FLOAT_EPSILON) Noise = detect_noise_from_mad(Buff);
    // INFO(Buff, "buff");
    return ((float) Noise);
}

/**************************************************************************/
/**********                     *******************************************/
/**********      dec_pos_max    *******************************************/
/**********                     *******************************************/
/**************************************************************************/
/* find max and ind_max of Tab in (only positive if SearchPositiv=True)   */
/**************************************************************************/
void dec_pos_max (fltarray& Tab, int &Ind_i, float &Val_Max,
                  Bool SearchPositiv=False);
void dec_pos_max (fltarray& Tab, int &Ind_i, float &Val_Max,
                  Bool SearchPositiv) {

   int i;
   float Val_Abs_Max;
   int Bande = 1;
   int Nx = Tab.nx();

   Val_Max = 0.;
   Val_Abs_Max = 0.;

   for (i = Bande; i < Nx-Bande; i++) {

      float Val_Abs = (SearchPositiv == False) ? ABS(Tab(i)) : Tab(i);
      if (Val_Abs > Val_Abs_Max) {

          Val_Abs_Max = Val_Abs;
          Val_Max = Tab(i);
          Ind_i = i;
      }
   }
}


/**************************************************************************/
/**********                        ****************************************/
/**********     dec_center_psf     ****************************************/
/**********                        ****************************************/
/**************************************************************************/
/* create a center psf of size Nx (signal in size), then normalize it     */
/**************************************************************************/
void dec_center_psf (fltarray& D_Beam, fltarray& Psf) {

   int Ind_i,i,Dep_i;
   float Val_Max;
   int Nx_Beam = D_Beam.nx();
   int Nx = Psf.nx();
   double Flux = 0.;

   dec_pos_max (D_Beam, Ind_i, Val_Max, True);

   Psf.init();
   for (i = 0; i < Nx_Beam; i++) {

      Dep_i = i - Ind_i  + Nx / 2;
      if ((Dep_i >= 0) && (Dep_i < Nx)) {

         Psf(Dep_i) = D_Beam(i);
	 Flux += Psf(Dep_i);
      }
    }
    for (i = 0; i < Nx; i++) Psf (i) = (float) (Psf(i) / Flux);
}




/**************************************************************************/
/**********                    ********************************************/
/**********      psf_extend    ********************************************/
/**********                    ********************************************/
/**************************************************************************/
/* create a psf of signal size (but not center on max)                    */
/**************************************************************************/
void psf_extend (const fltarray& Signal, fltarray& Signal_Out) {

    int Nx0 = Signal.nx();
    int Nx1 = Signal_Out.nx();
    int i1,i0,Depi;
    double Flux = 0.;

    Depi = (Nx1 - Nx0) / 2;

    for (i1 = 0; i1 < Nx1; i1++) {

       i0 = i1 - Depi;
       if ((i0 < 0) || (i0 >= Nx0))
            Signal_Out (i1) = 0.;
       else  {
          Signal_Out (i1) = Signal (i0);
	  Flux += Signal_Out (i1);
       }
    }
    for (i1 = 0; i1 < Nx1; i1++)
       Signal_Out (i1) = (float) (Signal_Out (i1) / Flux);
}


/**************************************************************************/
/**********                 ***********************************************/
/**********     psf_get     ***********************************************/
/**********                 ***********************************************/
/**************************************************************************/
/* compute the fft of the psf (signal size) in Psf_cf                     */
/**************************************************************************/
void psf_get (fltarray& InPsf, cfarray &Psf_cf, int SigNx, Bool PsfMaxShift)
{
   FFTN_1D FFT1D;
   int Nx = InPsf.nx();
   if (Nx < SigNx) Nx = SigNx;

   int Nx1=Nx;
   // Bool ModifSize=False;
   // dec_line_column (Nl, Nl1);
   // dec_line_column (Nc, Nc1);
   // if (Nl1 < Nc1) Nl1 = Nc1;
   // else Nc1 = Nl1;
   // if ((Nl != Nl1) || (Nc != Nc1)) ModifSize = True;

   fltarray OutPsf(Nx1, "Psf");
   if (PsfMaxShift == True) dec_center_psf (InPsf, OutPsf);
   else psf_extend (InPsf, OutPsf);
   // INFO(InPsf,"IN PSF");
   // INFO(OutPsf,"OUT PSF");

   cfarray OutPsf_cf(Nx1, "Psf_cf");
   for (int i=0;i<Nx1;i++) OutPsf_cf(i) = complex_f(OutPsf(i),0.);
   Psf_cf.resize(Nx1);
   // fft2d(OutPsf, Psf_cf);
   // fits_write_fltarr("xx_psf.fits", OutPsf);
   FFT1D.fftn1d (OutPsf_cf, Psf_cf);
}


/**************************************************************************/
/**********                    ********************************************/
/**********     test_ind       ********************************************/
/**********                    ********************************************/
/**************************************************************************/
int test_indice (int ind, int N)
{
    int Val;

    if (ind < 0) Val = - ind;
    else
    {
        if (ind >= N) Val = 2 * (N - 1) - ind;
        else Val = ind;
    }
    if ((Val >= N) || (Val < 0)) Val = -1;
    return (Val);
}

/**************************************************************************/
/**********                    ********************************************/
/**********     sig_extend     ********************************************/
/**********                    ********************************************/
/**************************************************************************/
/* extend signal with size of signal out (centered)                       */
/**************************************************************************/
void sig_extend (const fltarray& Signal, cfarray& Signal_Out) {

    int Nx0 = Signal.nx();
    int Nx1 = Signal_Out.nl();
    int i1,i0,Depi;
    complex_f Zero = complex_f(0., 0.);

    Depi = (Nx1 - Nx0) / 2;

    for (i1 = 0; i1 < Nx1; i1++) {

       i0 = test_indice (i1 - Depi, Nx0);
       if (i0 < 0) Signal_Out(i1) = Zero;
       // if ((i0 < 0) || (j0 < 0)) Signal_Out(i1) = Zero;
       else Signal_Out (i1) = complex_f(Signal(i0), 0.);
    }
}



/**************************************************************************/
/**********                    ********************************************/
/**********     sig_extract    ********************************************/
/**********                    ********************************************/
/**************************************************************************/
/* extract cenrtered signal of size Nx1 in Out Signal                     */
/**************************************************************************/
void sig_extract (const cfarray &Signal, fltarray &Signal_Out, float Norm)
{
    int Nx0 = Signal_Out.nx();
    int Nx1 = Signal.nx();
    int i1,i0,Depi;

    Depi = (Nx1 - Nx0) / 2;

    for (i0 = 0; i0 < Nx0; i0++) {

       i1 = i0 + Depi;
       Signal_Out (i0) = Signal (i1).real() / Norm;
    }
}



/**************************************************************************/
/**********                     *******************************************/
/**********     dec_inverse     *******************************************/
/**********                     *******************************************/
/**************************************************************************/
/* compute fft(sig)*conj(fft(psf)                                         */
/**************************************************************************/
void dec_inverse (fltarray& Signal, cfarray& Psf_cf, fltarray& Result,
                  float Eps=0.001, float Regul=0.);
void dec_inverse (fltarray& Signal, cfarray& Psf_cf, fltarray& Result,
                  float Eps, float Regul) {
   FFTN_1D FFT1D;
   int i;
   int Nx = Signal.nx();
   int Nx1 = Psf_cf.nx();
   cfarray O_cf (Nx1, "Buffer conv");
   // float FluxPsf = (FluxNorm == False) ? 1.: Psf_cf(Nl1/2,Nc1/2).real();
   float FluxPsf = 1. / Psf_cf(Nx1/2).real();

   if (Nx != Nx1) sig_extend (Signal, O_cf);
   else  {
       for (i=0;i<Nx;i++) O_cf(i) = complex_f(Signal(i),0);
   }
   FFT1D.fftn1d(O_cf);
   // fft2d (O_cf);
   float Den;

   for (i = 0; i < Nx1; i++) {

       Den = norm (Psf_cf(i)) + Regul;
       O_cf(i) *= conj(Psf_cf(i));
       if (Den > Eps) O_cf(i) /= Den;
       else O_cf(i) = complex_f(0.);
   }
   FFT1D.fftn1d(O_cf,True);
   // fft2d (O_cf, -1);
   if (Nx != Nx1) sig_extract (O_cf, Result, FluxPsf);
   else for (i = 0; i < Nx; i++) Result(i) = O_cf(i).real() / FluxPsf;
}



/**************************************************************************/
/**********                    ********************************************/
/**********     psf_convol     ********************************************/
/**********                    ********************************************/
/**************************************************************************/
/* compute convolution between signal and psf (in signal out)             */
/**************************************************************************/
void psf_convol (fltarray& Signal, cfarray& Psf_cf, fltarray& Signal_out,
                 Bool FluxNorm) {

   FFTN_1D FFT1D;
   int i,Nx = Signal.nx();
   int Nx1 = Psf_cf.nx();
   cfarray O_cf (Nx1, "Buffer conv");
   float FluxPsf = (FluxNorm == False) ? 1.: Psf_cf(Nx1/2).real();
   double Flux;

   if (Nx != Nx1) sig_extend (Signal, O_cf);
   else  {
       for (i = 0; i < Nx; i++) O_cf(i) = complex_f(Signal(i),0);
   }
   // fft2d (O_cf);
   FFT1D.fftn1d(O_cf);
   O_cf *= Psf_cf;
   // fft2d (O_cf, -1);
   FFT1D.fftn1d(O_cf,True);
   if (Nx != Nx1) {

      sig_extract (O_cf, Signal_out, FluxPsf);
      Flux = Signal.total() / Signal_out.total();
      for (i = 0; i < Nx; i++)
         Signal_out(i) = (float) (Signal_out(i) * Flux);
   }
   else
      for (i = 0; i < Nx; i++)
         Signal_out(i) = O_cf(i).real() / FluxPsf;

   // INFO(Imag,"IN ima");
   // INFO(Imag_out,"OUT ima");
}
void psf_convol (fltarray& Signal, cfarray& Psf_cf, Bool FluxNorm) {
   psf_convol (Signal, Psf_cf, Signal, FluxNorm);
}
void psf_convol (fltarray& Signal, fltarray& Psf, Bool FluxNorm) {
   cfarray Psf_cf;
   psf_get (Psf, Psf_cf, Signal.nx(), True);
   psf_convol (Signal, Psf_cf, Signal, FluxNorm);
}
void psf_convol (fltarray& Signal, fltarray& Psf, fltarray& Signal_out,
                  Bool FluxNorm) {
   cfarray Psf_cf;
   psf_get (Psf, Psf_cf, Signal.nx(), True);
   psf_convol (Signal, Psf_cf, Signal_out, FluxNorm);
}

/**************************************************************************/
/**********                          **************************************/
/**********      psf_convol_conj     **************************************/
/**********                          **************************************/
/**************************************************************************/
void psf_convol_conj (fltarray& Signal, cfarray& Psf_cf, fltarray& Signal_out) {

   FFTN_1D FFT1D;
   int Nx = Signal.nx();
   int Nx1 = Psf_cf.nx();
   cfarray O_cf (Nx1, "Buffer conv");
   float FluxPsf = Psf_cf(Nx1/2).real();

   if (Nx != Nx1) sig_extend (Signal, O_cf);
   else {
       for (int i=0; i<Nx; i++) O_cf(i) =  complex_f(Signal(i),0);
   }

   FFT1D.fftn1d(O_cf);
   for (int i=0; i<Nx1; i++) O_cf(i) *= conj(Psf_cf(i));

   FFT1D.fftn1d(O_cf, True);
   if (Nx != Nx1) sig_extract (O_cf, Signal_out, FluxPsf);
   else
      for (int i=0; i<Nx; i++)
         Signal_out(i) = O_cf(i).real() / FluxPsf;
}
void psf_convol_conj (fltarray& Signal, cfarray& Psf_cf) {
   psf_convol_conj (Signal, Psf_cf, Signal);
}






/**************************************************************************/
/**********                     *******************************************/
/**********      init_param     *******************************************/
/**********                     *******************************************/
/**************************************************************************/
void MR1D_Deconv::init_param()
{
    StatNoise = NOISE_GAUSSIAN;
    KeepImagn = False;
    Nx = 0;
    PsfMaxShift = True;
    GaussConv = False;
    Fwhm = 0.;
    Verbose = False;
    RegulParam = 0;
    MaxIter = DEFAULT_MAX_ITER_DECONV;
    DecMethod = DEFAULT_STD_DECONV;
    EpsCvg = 0.0000;
    IterCvg = 0.1;
    PositivConstraint = True;
    NormFlux = False;
    OptimParam = False;
    KillLastScale = False;
    UseICF = False;
    Noise_Sig = 0.;
    N_Sigma = 3.;
    FluxImag = 0.;
    TypeInit = DEC_INIT_FLAT;
    ModelData = NULL;
    WaveFilterResi = False;
    UseModel = False;
    Border = DEFAULT_BORDER;
    KeepImagn = True;
    AdjointRec = False;
    MaxMRCleanIter = DEFAULT_CLEAN_NITER;
    SupportDilate = False;
    UseMRCEnergy = True;
    CleanLastScale=True;
    CleanFirstScale=0;
    TabCB_cf=NULL;
    Filter = FIL_1D_THRESHOLD;
    ApplySoftRegul = False;
}



/**************************************************************************/
/**********                       *****************************************/
/**********     compute_resi      *****************************************/
/**********                       *****************************************/
/**************************************************************************/
void MR1D_Deconv::compute_resi() {

   int i;
   psf_convol (Obj, Psf_cf, Resi, True);
   if (KeepImagn == True) Signal_n = Resi;
   for (i=0; i< Nx; i++)  Resi(i) = Signal(i) - Resi(i);
}



/**************************************************************************/
/**********                   *********************************************/
/**********      sig_grad     *********************************************/
/**********                   *********************************************/
/**************************************************************************/
void MR1D_Deconv::sig_grad (fltarray& Gradient)
{

   int i;
   float Scale; // used by DEC_MEM
   double Logmin = 1e-9; // used by DEC_MEM

   // Regularization by residual filtering
   // INFO(Resi, "RESI");

   if ((WaveFilterResi == True) || (KillLastScale == True)) MR1D_Data.transform (Resi);
   if (WaveFilterResi == True)  ModelData->threshold(MR1D_Data);

   if (KillLastScale == True)
   {
      if (AdjointRec == False)
      {
         int s = MR1D_Data.nbr_band()-1;
         for(int w=0; w<MR1D_Data.size_scale_np(s); w++)  MR1D_Data(s,w) = 0;
	 MR1D_Data.recons (Resi);
      }
      else MR1D_Data.rec_adjoint (Resi, False);
   }
   else if (WaveFilterResi == True)
   {
      if (AdjointRec == False) MR1D_Data.recons (Resi);
      else MR1D_Data.rec_adjoint (Resi);
   }

   switch (DecMethod)
   {
      case DEC_MEM:
 	 Scale = Signal_n.total() / FluxImag;
	 for (i=0; i<Nx; i++)
	 {
            // Mult(i) += Imag(i)*Scale - Imag_n(i);
            float ValNum = Signal(i) * Scale;
            if (ValNum < Logmin) ValNum = Logmin;
	    float ValDen = Signal_n(i);
            if (ValDen < Logmin) ValDen = Logmin;
            Mult(i) += log(ValNum/ValDen);
	 }
         Gradient = Mult;
  	 psf_convol (Gradient, Psf_cf, True);
	 for (i=0; i<Nx; i++) Gradient(i) = exp(Gradient(i));
	 break;
       case DEC_MEM_MODEL:
	 Gradient = Resi;
         psf_convol_conj (Gradient, Psf_cf);
	 if (UseModel == True)
	 {
	    for (i=0; i<Nx; i++)
	       if (Obj(i) > MemModel(i))
	          Gradient(i) -= RegulParam*log(Obj(i)/MemModel(i));
	 }
	 else
	 {
	    for (i=0; i<Nx; i++)
               if (Obj(i) > FLOAT_EPSILON)
	          Gradient(i) -= RegulParam*(1+log(Obj(i)));
	 }
	 break;
      case DEC_MAP:
      case DEC_MR_MAP:
         for (i=0; i<Nx; i++)
            if (Signal_n(i) > FLOAT_EPSILON)
               Gradient(i) = Resi(i) / Signal_n(i);
            else
               Gradient(i) =  FLOAT_EPSILON;
 	 psf_convol_conj (Gradient, Psf_cf);
	 for (i=0; i<Nx; i++)
            Gradient(i) = exp( Gradient(i));
 	 break;
     case DEC_MR_LUCY:
     case DEC_LUCY:
     case DEC_MARKOV_LUCY:
         for (i=0; i<Nx; i++)
	 {
            if (Signal_n(i) > FLOAT_EPSILON)
               Gradient(i) = Resi(i) / Signal_n(i);
            else  Gradient(i) = 0.;
         }
 	 psf_convol_conj (Gradient, Psf_cf);
         Gradient *= Obj;
         // if ((RegulParam > 0) &&
         //   RSig.obj_regul (Obj, Gradient, RegulParam*Noise_Sig);
         break;
      case DEC_MARKOV:
         Gradient = Resi;
         psf_convol_conj (Gradient, Psf_cf);
         if (RegulParam > 0) RSig.obj_regul(Obj, Gradient, RegulParam*Noise_Sig);
         break;
      case DEC_TIKHONOV:
         Gradient = Resi;
         psf_convol_conj (Gradient, Psf_cf);
         if (RegulParam > 0) RSig.obj_regul(Obj, Gradient, RegulParam);
         break;
      case DEC_MR_GRADIENT:
      case DEC_GRADIENT:
         Gradient = Resi;
         psf_convol_conj (Gradient, Psf_cf);
         // if (RegulParam > 0) RSig.obj_regul(Obj, Gradient, RegulParam*Noise_Sig);
         break;
      case DEC_MR_CLEAN:
      case DEC_MR_CITTERT:
      case DEC_CITTERT:
  	 Gradient = Resi;
         // if (RegulParam > 0)
         //   RSig.obj_regul (Obj, Gradient, RegulParam*Noise_Sig);
         break;

      default:
         cerr << "mr_deconv: Not implemented in this procedure ... ";
         cerr << endl;
         exit (0);
         break;
   }
}


/**************************************************************************/
/**********                       *****************************************/
/**********      sig_gaussian     *****************************************/
/**********                       *****************************************/
/**************************************************************************/

fltarray sig_gaussian (int Nx, float Fwhm, int Indi=-1);
fltarray sig_gaussian (int Nx, float Fwhm, int Indi)
{
    int Delta_i;
    float sigma,sigma2,Energ = 0.;
    register int i;
    fltarray *Result = new fltarray(Nx, "gaussienne");
    int Ci=Indi;

    if (Ci < 0) Ci = Nx / 2;
    sigma = 0.5 * Fwhm / sqrt (2. * log ((double) 2.));
    sigma2 = -2. * sigma * sigma;

    for (i=0; i<Nx; i++)
    {
        Delta_i = (i - Ci) * (i - Ci);
	(*Result)(i) = exp((double)((float) Delta_i / sigma2));
	Energ +=  (*Result)(i);
    }
    return (*Result);
}

/**************************************************************************/
/**********                    ********************************************/
/**********      norm_flux     ********************************************/
/**********                    ********************************************/
/**************************************************************************/
void norm_flux (fltarray& Signal, float Val_Norm=1.);
void norm_flux (fltarray& Signal, float Val_Norm) {

    int i;
    float F = Signal.total() / Val_Norm;
    for (i=0; i<Signal.nx(); i++)
       Signal(i) /= F;
}



/**************************************************************************/
/**********                      ******************************************/
/**********      sigma_clip      ******************************************/
/**********                      ******************************************/
/**************************************************************************/
void sigma_clip (const fltarray& Signal, float &Moy, float &Sig, int Nit=3);
void sigma_clip (const fltarray& Signal, float &Moy, float &Sig, int Nit) {

    int It, i;
    double S0,S1,S2,Sm=0,x,Mean,Sigma=0.;
    int Nx = Signal.nx();

    Mean = 0.;
    for (It=0; It<Nit; It++) {

       S0 = S1 = S2 = 0.;
       for (i=0; i<Nx; i++)  {
          x = Signal(i);
          if ((BadPixel == False) || (ABS(x-BadPixalVal) > FLOAT_EPSILON)) {
	     if ((It == 0) || (ABS(x - Mean) < Sm)) {
	        S0 ++; S1 += x; S2 += x*x;
	     }
          }
       }
       if (S0 == 0) S0=1;
       Mean = S1 / S0;
       Sigma = S2/S0-  Mean * Mean;
       if ( Sigma >= 0.) Sigma = sqrt(Sigma);
       else Sigma = 0.;
       Sm = 3. * Sigma;
    }
    Moy= (float) Mean;
    Sig =(float) Sigma;
}




/**************************************************************************/
/**********                       *****************************************/
/**********      convol_gauss     *****************************************/
/**********                       *****************************************/
/**************************************************************************/

void MR1D_Deconv::convol_gauss (fltarray& Data, float FWHM)
{
   fltarray Gauss;
   Gauss = sig_gaussian (Data.nx(), FWHM);
   norm_flux (Gauss);
   psf_convol (Data, Gauss, True);
}


/**************************************************************************/
/**********                          **************************************/
/**********      kill_last_scale     **************************************/
/**********                          **************************************/
/**************************************************************************/
void kill_last_scale (MR_1D & MR1D_Data) {

   int s = MR1D_Data.nbr_band()-1;
   for (int w=0; w < MR1D_Data.size_scale_np(s); w++)
      MR1D_Data(s,w) = 0;
}


/**************************************************************************/
/**********                      ******************************************/
/**********      init_deconv     ******************************************/
/**********                      ******************************************/
/**************************************************************************/
void MR1D_Deconv::init_deconv(fltarray* FirstGuess, fltarray* ICF)
{
   float Flux;
   Nx = Signal.nx();
   cout << Nx << endl;
   Bool NonAddMethod = False;
   int i;

   // set NonAddMethod to True for LUCY and MAP
   if (    (DecMethod == DEC_LUCY) || (DecMethod == DEC_MR_LUCY)
        || (DecMethod == DEC_MAP)  || (DecMethod == DEC_MR_MAP))
	                                         NonAddMethod = True;

   // Multiplyer image allocation MEM in memthod
   if ((DecMethod == DEC_MEM) || (DecMethod == DEC_MEM_MODEL))
                                                 Mult.alloc(Nx,"mult");
   // if ((DecMethod == DEC_LUCY) ||
   //      ((OptimParam == True) && (StatNoise == NOISE_POISSON)))
   //         KeepImagn = True;
   if (KeepImagn == True) Signal_n.alloc(Nx,"Signal_n");

   // if an ICF is introduced, it can be either a Gaussian
   // or an image given by the user
   // the PSF must be convolved by the ICF
   if (GaussConv == True) convol_gauss(Psf, Fwhm);
   else if (ICF != NULL) {
      Sig_ICF.alloc(Psf_cf.nx(), "ICF");
      dec_center_psf (*ICF, Sig_ICF);
      UseICF = True;
      FFT1D.convolve(Psf, Sig_ICF, Psf);
   }

   // centering the PSF (normally the max of PSF is at the center of the image)
   // the size of Psf_cf is a power of 2, which is not necessary the size of
   // Psf, or the image size
   psf_get (Psf, Psf_cf, Nx, PsfMaxShift);

   // allocate Obj if FirstGuess==NULL
   if (FirstGuess != NULL)
      TypeInit = DEC_INIT_GUESS;
   else
      Obj.alloc(Nx, "object");

   // initialize Obj in function of TypeInit (with 0, flat, Signal or guess)
   switch (TypeInit)
   {
      case DEC_INIT_ZERO: Obj.init();  break;
      case DEC_INIT_FLAT:
             Flux = Signal.total() / (float) (Nx);
  	     Obj.init(Flux);
	     break;
      case DEC_INIT_IMA:
             Obj = Signal;
	     break;
      case DEC_INIT_GUESS:
             Obj = *FirstGuess;
             break;
   }

   // allocate resi
   Resi.alloc(Nx,"resi");

   /* Noise estimation in the Data */
   if (StatNoise == NOISE_GAUSSIAN)
   {
        if (Noise_Sig < FLOAT_EPSILON)
	{
	   Noise_Sig = detect_noise_from_med (Signal);
	   if (Verbose == True)
	      cout << "Sigma Noise = " << Noise_Sig << endl;
	}
   }
   if (Noise_Sig < FLOAT_EPSILON) Noise_Sig = 1.;

   if ((WaveFilterResi == True) || (KillLastScale == True))
   {

      WaveFilterResi = True;
      // MR1D_Data.alloc(Nl, Nc, ModelData->nbr_scale(),
      //              ModelData->type_trans(), "MR1D_Data");
      MR1D_Data.alloc (Nx, ModelData->type_trans(), ModelData->nbr_scale(),
                       ModelData->filter_bank(), ModelData->TypeNorm);

      MR1D_Data.Border = Border;
      ModelData->model(Signal, MR1D_Data);

      if (KillLastScale == True)
      {
         MR1D_Data.transform(Obj);
         if (TypeInit == DEC_INIT_IMA)  ModelData->threshold(MR1D_Data);
         kill_last_scale (MR1D_Data);
         MR1D_Data.recons(Obj);
      }
   }

   // Zero value in the first guess leads to problem with
   // Lucy and MAP methods
   if (NonAddMethod == True)
   {
      for (i=0; i < Nx; i++)
      {
 	  if (Obj(i) < FLOAT_EPSILON) Obj(i) = FLOAT_EPSILON;
	  if (Signal(i) < FLOAT_EPSILON) Signal(i) = FLOAT_EPSILON;
      }
   }
   switch (DecMethod)
   {
      case DEC_INVERSE:
      case DEC_MARKOV_LUCY:
      case DEC_MARKOV:
           RSig.GradOperType = OPER_MARKOV_2;
           break;
      case DEC_TIKHONOV:
           RSig.GradOperType = OPER_LAPLACIAN;
           break;
      case DEC_MR_MEM:
      case DEC_MR_MEM_NOISE:
      case DEC_MEM:
      case DEC_MEM_MODEL:
      case DEC_MAP:
           break;
      case DEC_CITTERT:
      case DEC_GRADIENT:
      case DEC_LUCY:
      case DEC_CLEAN:
      case DEC_MR_CITTERT:
      case DEC_MR_GRADIENT:
      case DEC_MR_LUCY:
      case DEC_MR_MAP:
           if (RegulParam > 0) ApplySoftRegul = True;
	   if (Verbose == True) cout << "Wavelet Regularization " << RegulParam << endl;
           break;
      case DEC_MR_CLEAN:
      case DEC_MR_VAGUELET:
      case DEC_MR_INVERSE:
           break;
   }
}


/**************************************************************************/
/**********                             ***********************************/
/**********      find_optim_poisson     ***********************************/
/**********                             ***********************************/
/**************************************************************************/
float MR1D_Deconv::find_optim_poisson (fltarray& Gradient) {

  float Cvgparam = 1., OldCvgparam;
  fltarray Delta(Nx,"Buff");
  int i,Iter=0;
  int Np = Nx;
  float MinVal=1., MaxVal=100.;
  float Val, Deriv1=0., Deriv2=0.;
  if (RegulParam > 0.5) MinVal = 1. / (2.*RegulParam);

  // sear the maximum value
  for (i=0; i < Np; i++)
    if (Gradient(i) < - FLOAT_EPSILON)
    {
       Val = - Obj(i) / Gradient(i);
       if (Val < MaxVal) MaxVal = Val;
    }
  // cout << "interval: Min = " << MinVal << " Max = " << MaxVal << endl;

  if (MaxVal > 0.) {

     // dec_convol (Gradient, Psf_cf, Delta);
     psf_convol (Gradient, Psf_cf, Delta, True);
     // Cvgparam = (MinVal + MaxVal) / 2.;
     do {

        Iter++;
        OldCvgparam = Cvgparam;
        // Calculate first and second derivative
        for (i=0; i<Np; i++) {

	   // First derivative
	   Val = Signal_n(i) + Cvgparam*Delta(i);
	   if (ABS(Val) > FLOAT_EPSILON)
	         Deriv1 += Signal(i)*Delta(i) / Val ;

           // Second derivative
	   Val *= Val;
	   if (ABS(Val) > FLOAT_EPSILON)
	          Deriv2 -= Signal(i)*Delta(i)*Delta(i) / Val;
	}
        Cvgparam -= Deriv1/ Deriv2;

	// cout << "Iter " << Iter << " Cvgparam = " << Cvgparam << endl;
     } while ((Iter < 10) && ((ABS(Cvgparam-OldCvgparam) > 0.1*OldCvgparam)));
     if (Cvgparam > MaxVal) Cvgparam = MaxVal;
  }
  if (Cvgparam < MinVal) Cvgparam = MinVal;
  return Cvgparam;
}



/**************************************************************************/
/**********                              **********************************/
/**********      find_optim_tikhonov     **********************************/
/**********                              **********************************/
/**************************************************************************/
float MR1D_Deconv::find_optim_tikhonov(fltarray& Gradient) {

   float Cvgparam = IterCvg;
   fltarray Buff(Nx,"Buff");
   float LapG,LapO, Num, Den;
   int i;

   if (RegulParam > 0.) {

      // dec_convol (Gradient, Psf_cf, Buff);
      psf_convol(Gradient, Psf_cf, Buff, True);
      Num = (Buff*Resi).total();
      Den = Buff.energy();
      for (i=0; i<Nx; i ++) {

         LapG = RSig.laplacian_val(Gradient,i);
	 LapO = RSig.laplacian_val(Obj,i);
	 Num -= RegulParam*LapG*LapO;
         Den += RegulParam*LapG*LapG;
      }
      Cvgparam =  Num / Den;
   } else
      Cvgparam = find_optim_xi2(Gradient);

   return Cvgparam;
}


/**************************************************************************/
/**********                         ***************************************/
/**********      find_optim_xi2     ***************************************/
/**********                         ***************************************/
/**************************************************************************/
float MR1D_Deconv::find_optim_xi2(fltarray& Gradient) {

   float Cvgparam = IterCvg;
   fltarray Buff(Nx,"Buff");
   float MinVal = 1.;
   if (RegulParam > 0.5) MinVal = 1. / (2.*RegulParam);
   // dec_convol (Gradient, Psf_cf, Buff);
   psf_convol(Gradient, Psf_cf, Buff, True);
   Cvgparam = (Buff*Resi).total()/Buff.energy();
   if (Cvgparam < MinVal) Cvgparam = MinVal;
   else if (Cvgparam > 10) Cvgparam = 10.;
   return Cvgparam;
}



/**************************************************************************/
/**********                        ****************************************/
/**********      im_find_optim     ****************************************/
/**********                        ****************************************/
/**************************************************************************/
float MR1D_Deconv::im_find_optim (fltarray& Gradient) {

   float Cvgparam = 1.;
   if (DecMethod == DEC_TIKHONOV) Cvgparam = find_optim_tikhonov(Gradient);
   else if (StatNoise == NOISE_GAUSSIAN)  Cvgparam = find_optim_xi2(Gradient);
   else if (StatNoise == NOISE_POISSON) Cvgparam = find_optim_poisson(Gradient);

   if (Verbose == True) {
      cout << "Optim Val = " << Cvgparam << endl;
   }

   if (Cvgparam < 0) Cvgparam = 0.;
   return Cvgparam;
}



/**************************************************************************/
/**********                     *******************************************/
/**********     fonctional      *******************************************/
/**********                     *******************************************/
/**************************************************************************/
float MR1D_Deconv::fonctional() {

   float Val,Func = 0.;
   int i;

   // Maxumum likehood part of the functional
   if (StatNoise == NOISE_GAUSSIAN) Func = Resi.energy();
   else { // POISSON noise

      Func = Signal.energy();  // we add a constant in order to have
                             // a positive functional
      for (i=0; i<Nx; i ++)
        if (Signal_n(i) > 0)
           Func -= Signal(i)*log(Signal_n(i));
   }

   // regularized part of the functional
   switch (DecMethod)
   {
      case DEC_TIKHONOV:
         for (i=0; i<Nx; i ++) {
            Val = RSig.laplacian_val(Obj,i);
            Func += Val*Val*RegulParam;
         }
         break;
      case DEC_MEM:
         for (i=0; i<Nx; i ++)
            if (Obj(i) > 0)
               Func += RegulParam*Obj(i)*log(Obj(i));
         break;
      case DEC_MEM_MODEL:
         if (UseModel == False) {
            for (i=0; i<Nx; i ++)
               if (Obj(i) > 0)
                  Func += RegulParam*Obj(i)*log(Obj(i));
         } else {
             for (i=0; i<Nx; i ++)
                if (Obj(i) > 0)
                   Func +=    RegulParam*((-Obj(i) + MemModel(i))
                            + Obj(i)*log(Obj(i)/MemModel(i)));
         }
         break;
      case DEC_MARKOV_LUCY:
      case DEC_MARKOV:
      case DEC_MR_MEM:
      case DEC_MR_MEM_NOISE:
      case DEC_MAP:
      case DEC_CITTERT:
      case DEC_GRADIENT:
      case DEC_INVERSE:
      case DEC_LUCY:
      case DEC_CLEAN:
      case DEC_MR_CITTERT:
      case DEC_MR_GRADIENT:
      case DEC_MR_LUCY:
      case DEC_MR_MAP:
      case DEC_MR_CLEAN:
      case DEC_MR_VAGUELET:
      case DEC_MR_INVERSE:
           break;
   }

   return Func;
}



/**************************************************************************/
/**********                         ***************************************/
/**********      im_iter_deconv     ***************************************/
/**********                         ***************************************/
/**************************************************************************/

void MR1D_Deconv::sig_iter_deconv()
{
    fltarray Gradient(Nx,"Info");
    float Func, OldFunc=1e10, Sigma;
    int i,Iter=0;
    Bool Stop = False;
    float Cvgparam = IterCvg;
    float Sigma2N = Noise_Sig*Noise_Sig*Nx;
    if (NormFlux == True) FluxImag = Signal.total();
    else FluxImag=0.;

    // Residual estimation
    compute_resi();
    // INFO(Resi, "RESI");
    // start iterations
    if ((Verbose == True) && (RegulParam > 0))
           cout << "Start iterate: Regul =  " <<  RegulParam << endl;
    do {

       // compute the gradient
       sig_grad(Gradient);

       if (OptimParam == True) Cvgparam = im_find_optim(Gradient);

       // calculate the next object estimation with positivity  constraint
       if ((DecMethod == DEC_MAP) || (DecMethod == DEC_MR_MAP)) {
          for (i = 0; i < Nx; i++)
             Obj(i) *= Gradient(i);
       } else if (DecMethod == DEC_MEM)
          Obj = Gradient;
       else {
          for (i=0; i<Nx; i++)
             Obj(i) += Cvgparam * Gradient(i);
       }

       if (ApplySoftRegul == True)
          RSig.mr1d_soft_threshold(Obj,  RegulParam, Noise_Sig);
       if (PositivConstraint == True)
          for (i=0; i<Nx; i++)
             if (Obj(i) < 0) Obj(i)=0;
       if ((NormFlux == True) && (KillLastScale == False)) {

          float FluxObj = Obj.total();
	  FluxObj = FluxImag / FluxObj;
	  for (i= 0; i < Nx; i++) Obj(i) *= FluxObj;
       }

       // next residual
       compute_resi();
       Sigma = Resi.energy() / Sigma2N;

       if ((Verbose == True) && ((Iter >= 0) && (Iter % 1 == 0))) {

          Func = fonctional() / (float) (Nx);
          printf("%d: Sigma = %f, Func = %f\n", Iter, Sigma,  Func);
       }

       Iter ++;
       if (Iter > MaxIter) Stop = True;
       else if (EpsCvg > FLOAT_EPSILON) {

          if (OldFunc < Sigma) Stop = True;
          if ((ABS(OldFunc - Sigma)) < EpsCvg) Stop = True;
       }
       OldFunc = Sigma;
    } while (Stop == False);

    if (Verbose == True)
    {
       Obj.info ("Solution: ");
       Resi.info ("Resi: ");
    }
}

/**************************************************************************/
/**********                     *******************************************/
/**********      sig_deconv     *******************************************/
/**********                     *******************************************/
/**************************************************************************/
void MR1D_Deconv::sig_deconv (fltarray* FirstGuess, fltarray* ICF) {


   float Mean,Sigma;                 // used in DEC_CLEAN
   float GammaClean = RegulParam;    // used in DEC_CLEAN
   Bool IterMethod = True;

   switch (DecMethod)
   {
      case DEC_INVERSE:
            IterMethod = False;
            StatNoise = NOISE_GAUSSIAN;
            init_deconv(FirstGuess,  ICF);
	    dec_inverse (Signal, Psf_cf, Obj, EpsCvg);
	    compute_resi();
            break;
      case DEC_TIKHONOV:
      case DEC_MARKOV:
      case DEC_GRADIENT:
            StatNoise = NOISE_GAUSSIAN;
 	    break;
      case DEC_MAP:
            TypeInit = DEC_INIT_FLAT;
            OptimParam = False;
	    StatNoise = NOISE_POISSON;
 	    break;
      case DEC_MARKOV_LUCY:
      case DEC_LUCY:
            TypeInit = DEC_INIT_FLAT;
            StatNoise = NOISE_POISSON;
 	    break;
      case DEC_CITTERT:
      case DEC_MEM:
      case DEC_MEM_MODEL:
            TypeInit = DEC_INIT_FLAT;
            NormFlux = True;
            OptimParam = False;
	    StatNoise = NOISE_GAUSSIAN;
 	    break;
      case DEC_CLEAN:
            IterMethod = False;
            StatNoise = NOISE_GAUSSIAN;
 	    Nx = Signal.nx();
	    Obj.alloc(Nx, "Obj");
	    Resi.alloc(Nx, "Resi");
            Resi = Signal;
            norm_flux (Psf, 1.);
            if (Noise_Sig < FLOAT_EPSILON)
               Noise_Sig = detect_noise_from_med (Resi);
            sigma_clip (Resi,Mean,Sigma);
	    cout << "Mean+N_Sigma*Noise_Sig = " << Mean+N_Sigma*Noise_Sig << endl;
	    cout << "GammaClean = " << GammaClean << endl;
	    if (GaussConv == True) cout << "GaussConv = " << Fwhm << endl;

            dec_clean (Resi, Psf, Obj, 0., Mean+N_Sigma*Noise_Sig,
                       GammaClean, MaxIter,
                       DEFAULT_CLEAN_KAVER, DEFAULT_CLEAN_ADDRESI,
                       True, False, True);
            if (ICF != NULL) {
               Sig_ICF.alloc(Nx, "Psf");
               dec_center_psf (*ICF, Sig_ICF);
	       UseICF = True;
	    }
	    break;
      case DEC_MR_VAGUELET:
           IterMethod = False;
	   WaveFilterResi = True;
	   init_deconv(FirstGuess,ICF);
	   vaguelet();
	   compute_resi();
           break;
      case DEC_MR_CITTERT:
      case DEC_MR_GRADIENT:
           WaveFilterResi = True;
	   // NormFlux = True;
           break;
      case DEC_MR_LUCY:
      case DEC_MR_MAP:
           TypeInit = DEC_INIT_FLAT;
           WaveFilterResi = True;
	   // NormFlux = True;
  	   break;
      case DEC_MR_CLEAN:
      case DEC_MR_MEM:
      case DEC_MR_MEM_NOISE:
      case DEC_MR_INVERSE:
      default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
   }
   if (IterMethod == True)
   {
      init_deconv(FirstGuess,ICF);
      sig_iter_deconv();
   }

   if (GaussConv == True)  convol_gauss(Obj, Fwhm);
   else if (UseICF == True) psf_convol(Obj, Sig_ICF, Obj, True);

}

/************************************************************************/

void MR1D_Deconv::va_inverse_data(fltarray  &NoiseSimu)
{
    if (StatNoise == NOISE_GAUSSIAN)  simu_gaussian_noise(NoiseSimu, Noise_Sig);
    else
   {
      // Find the noise map by filtering the data
      MR1DFiltering CFilter(*ModelData, FIL_1D_THRESHOLD);
      CFilter.Verbose = Verbose;
      CFilter.filter(Signal, Obj);

      for (int i=0; i < Nx; i++) NoiseSimu(i) = Signal(i) - Obj(i);
   }

   // inverse deconvolution
   // dec_std_inverse (Imag, Psf_cf, Obj);
   // dec_std_inverse (NoiseSimu, Psf_cf, NoiseSimu);
   dec_inverse(Signal, Psf_cf, Obj, EpsCvg,RegulParam);
   dec_inverse(NoiseSimu, Psf_cf, NoiseSimu, EpsCvg,RegulParam);
}

/************************************************************************/

void MR1D_Deconv::vaguelet()
{
   int b,i;
   int Nb = MR1D_Data.nbr_band()-1;
   fltarray ThresholdMax(Nb+1);
   fltarray NoiseSimu(Nx);
   fltarray Scale;

   // inverse deconvolution on both imag and noise map
   va_inverse_data(NoiseSimu);
   // io_write_ima_float("xx_inv.fits", Obj);
   // io_write_ima_float("xx_noise.fits", NoiseSimu);

   fltarray TabSigma(Nb);
   MR1D_Data.transform(NoiseSimu);
   for (b=0; b < Nb; b++)
   {
      MR1D_Data.scale (Scale, b);
      ThresholdMax(b) = Scale.sigma() * ModelData->NSigma[b];
   }
   for (b=0; b < Nb; b++) cout << "Band " << b << " sigma = " << TabSigma(b) << " T = " << ThresholdMax(b) << endl;

    // Update the multiresolution support
    MR1D_Data.transform(Obj);
    for (b=0; b < Nb; b++)
    for (i=0; i < MR1D_Data.size_scale_np(b); i++)
    {
        if (ModelData->OnlyPositivDetect == True)
	{
	  if (MR1D_Data(b,i) > ThresholdMax(b)) ModelData->support(b,i) = 1;
        }
        else
	{
	   if (ABS(MR1D_Data(b,i)) > ThresholdMax(b))  ModelData->support(b,i) = 1;
 	}
   }

   type_1d_filter Filter = FIL_1D_THRESHOLD;
   MR1DFiltering CFilter(*ModelData, Filter);
   NoiseSimu=Obj;
   CFilter.Max_Iter = MaxIter;
   CFilter.Epsilon = EpsCvg;
   CFilter.Sup_Set = False;
   CFilter.Verbose = Verbose;
   CFilter.KillLastScale = KillLastScale;
   CFilter.filter(NoiseSimu, Obj);   // mr_support_filter(NoiseSimu, MR1D_Data, Obj);
   if (PositivConstraint == True)
           for (i=0; i<Nx; i++)   if (Obj(i) < 0) Obj(i)=0;
}

/****************************************************************************/
