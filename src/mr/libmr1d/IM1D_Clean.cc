/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  08/04/03 
**    
**    File:  IM1D_Clean.cc
**
*******************************************************************************
**
**    DESCRIPTION  Signal Deconvolution: CLEAN METHOD
**    -----------  
**                 
**
*******************************************************************************
**
** void dec_clean (fltarray & Signal, fltarray & Psf, fltarray & Signal_Out, 
**                 float Fwhm, float Noise, float Gamma, int Niter, 
**                 int K_Aver, Bool AddResi, Bool SearchPositiv, Bool UseEnergy,
** 		   Bool Verbose)
**
** Deconvolution by the CLEAN method
** Signal = in: dirty map
** Psf = in: dirty beam
** Signal_Out = out: clean map
** Fwhm = in: FWHM of the clean beam (if Fwhm > 0) the clean map is convolved
**            by the cleab beam
** Noise: in: detection level
** Gamma: in: loop gain
** Niter: in: maximum number of iterations
** K_Aver: in: end loop parameter: stop if Max < K_Aver*mean(residual map)
** AddResi: in: if true, the residual map is added to the solution
** SearchPositiv: in: if true, only positive peaks are searched
** UseEnergy: in: if true, maxima are search on the energy map
** Verbose: in: verbose mode
** 
******************************************************************************/
		

#include "FFTN_1D.h"

fltarray sig_gaussian (int Nx, float Fwhm, int Indi=-1);
void dec_pos_max (fltarray& Tab, int &Ind_i, float &Val_Max, 
                  Bool SearchPositiv=False);
                  
/**************************************************************************/
/**********                            ************************************/
/**********      substract_signal      ************************************/
/**********                            ************************************/
/**************************************************************************/                  
void fft1d_conv (const fltarray& Sig_in1, const fltarray& Sig_in2, 
                 fltarray &Sig_out) {

    int Nx = Sig_in1.nx();
    cfarray  Sig_in1_cf(Nx, "Im_in1_cf");
    cfarray  Sig_in2_cf(Nx, "Im_in2_cf");
    FFTN_1D FFT1D;

    for (int i=0; i<Nx; i++) {
       Sig_in1_cf(i) = complex_f(Sig_in1(i),0.);
       Sig_in2_cf(i) = complex_f(Sig_in2(i),0.);     
    }
    FFT1D.fftn1d (Sig_in1_cf);
    FFT1D.fftn1d (Sig_in2_cf);
    Sig_in2_cf *= Sig_in1_cf;
    
    FFT1D.fftn1d (Sig_in2_cf, Sig_in1_cf, True);

    for (int i=0; i<Nx; i++)
       Sig_out(i) = Sig_in1_cf(i).real();
}      
       
       
       
                  
/**************************************************************************/
/**********                            ************************************/
/**********      substract_signal      ************************************/
/**********                            ************************************/
/**************************************************************************/
float substract_signal (fltarray& Dirty_Map, fltarray& Dirty_Beam, 
                        int Max_Beam_i, int Max_i, float  Coef_Mul) {
 
// Soustraction de la Dirty_Beam centree en Max_Beam_i, Max_Beam_j
//    a la dirty map. 
//    Dirty_Map[i,j] = Dirty_Map[i,j] - 
//      Coef_Mul*Dirty_Beam [ Max_Beam_i + i - Max_i, Max_Beam_j + j - Max_Beam_j ]
// 
//    La procedure renvoie la moyenne en valeur absolue de l'image .
// 

    register int i;
    float Aver_Abs = 0.;
    int Dep_Beam_i;
    int Nx = Dirty_Map.nx();
    
    for (i=0; i<Nx; i++) {
    
       Dep_Beam_i = i - Max_i + Max_Beam_i;
       Aver_Abs += ABS(Dirty_Map(i));  
       if ((Dep_Beam_i >= 0) && (Dep_Beam_i < Nx)) {
          Dirty_Map(i) -= Dirty_Beam(Dep_Beam_i)*Coef_Mul;
       }
    }
    Aver_Abs /= Nx;
    return (Aver_Abs);
}


/**************************************************************************/
/**********                      ******************************************/
/**********      dec_hogbom      ******************************************/
/**********                      ******************************************/
/**************************************************************************/
void dec_hogbom (fltarray& Dirty_Map, fltarray& Dirty_Beam, 
                 float Gamma, float Noise, fltarray& Clean_Map,
                 int Niter, int K_Aver, Bool SearchPositiv, 
                 Bool Verbose) {
                        
    int iteration=0;
    int Max_i;
    float Old_Max, Max_Sig = 0., Coef_Mul, Gain;
    float Aver_Abs = 0.;
    int Max_Beam_i;
    float Max_Beam;

    Clean_Map.init();

    dec_pos_max (Dirty_Beam, Max_Beam_i, Max_Beam, SearchPositiv);

    /* begin of the loop */
    Gain = Gamma / Max_Beam;
    
    /* search the max in the dirty map  */

    dec_pos_max (Dirty_Map, Max_i, Max_Sig, SearchPositiv);
    Old_Max = ABS(Max_Sig) + 1.;
    
cout << "Max_Sig:" << Max_Sig << endl;
cout << "K_Aver*Aver_Abs:" << K_Aver*Aver_Abs << endl;
cout << "Noise:" << Noise << endl;
cout << "iteration:" << iteration << endl;
cout << "Old_Max:" <<Old_Max  << endl;

    while (     (ABS(Max_Sig) >= K_Aver*Aver_Abs) 
             && (ABS(Max_Sig) >= Noise) 
             && (iteration < Niter)
             && ( (ABS(Old_Max) > ABS(Max_Sig)) || (Max_Sig*Old_Max > 0.))
          ) {
          
       Coef_Mul = Max_Sig * Gain;
       Clean_Map(Max_i) += Coef_Mul;
  
       // Subtract the dirty beam to dirty map. the dirty beam is centered
       // on Imax, and weighted by Gamma 
       Aver_Abs = substract_signal (Dirty_Map, Dirty_Beam,  Max_Beam_i,  
	                          Max_i, Coef_Mul);
		      
       iteration ++;
 
       if ((Verbose == True) && (iteration % 10 == 0))
          cout << iteration << " : Max_Ima = (" << Max_i << ") : "
               << Max_Sig << endl;
                
       // search the max in the dirty map 
       Old_Max = Max_Sig;
       dec_pos_max (Dirty_Map, Max_i, Max_Sig, SearchPositiv);
    }
    
    if (Verbose == True) {
    
       cout << "End loop: number of iterations = " << iteration 
            << ", Max_Ima = (" << Max_i << ") : " << Max_Sig << endl;
        
       cout << "end loop criterion is :" << endl;                                      
       if (ABS(Max_Sig) < K_Aver*Aver_Abs)
          cout << "  Max_Sig < " << K_Aver << " * MeanResidu" << endl;
       if (iteration >= Niter) 
          cout << "Nbr_of_iter = MAX_ITER"<< endl;
       if (ABS(Max_Sig) < Noise ) 
          cout << "Max_Sig < Noise Level"   << endl;
       if (ABS(Old_Max) < ABS(Max_Sig) )  
          cout << "Old_Max > New_Max" << endl;
    }
}



/**************************************************************************/
/**********                      ******************************************/
/**********      map_energy      ******************************************/
/**********                      ******************************************/
/**************************************************************************/
/* energy of a signal : E = som_i( I_i^2 + Der_x(I)_i^2)                  */
/**************************************************************************/
void map_energy (fltarray& Signal, fltarray& Signal_x, fltarray& Plan_Energ) {

    int i;
    int Nx = Signal.nx();
    
    for (i=0; i<Nx; i++)
       Plan_Energ(i)= Signal_x(i)*Signal_x(i) + Signal(i)*Signal(i);
}



/**************************************************************************/
/**********                            ************************************/
/**********      fft_derive_image      ************************************/
/**********                            ************************************/
/**************************************************************************/
/*  Calcul de la derive d'un signal dans la directions x par utilisation  */
/*  de la transormee de Fourier : Der_x (I) = TF-1 [ TF(I)*2*i*PI*u ]     */
/**************************************************************************/
void fft_derive_image (fltarray& Signal, fltarray& Signal_x) {

    int i;
    float u;
    int Nx = Signal.nx();
    cfarray Pict (Nx, "Buffer");
    cfarray Pict_x (Nx, "Buffer");
    FFTN_1D FFT1D;

    for (i=0; i<Nx; i++) 
       Pict(i) = complex_f(Signal(i),0.);
    FFT1D.fftn1d(Pict);
    
    // derivative calculation
    for (i=0; i<Nx; i++) {
        u = (float)(i - Nx/2) / (float) Nx;
	Pict_x(i) = Pict(i)*complex_f(0,2.*PI*u);
    }

    FFT1D.fftn1d(Pict_x, True);
 
    for (i=0; i<Nx; i++) 
        Signal_x(i) = Pict_x(i).real();
    
}


/**************************************************************************/
/**********                               *********************************/
/**********      clean_hogbom_energy      *********************************/
/**********                               *********************************/
/**************************************************************************/
void clean_hogbom_energy (fltarray& Dirty_Map, fltarray& Dirty_Beam, 
                 float Gamma, float Noise, fltarray& Clean_Map,
		 fltarray& Tab_Dirty_X, fltarray& Dirty_Beam_X,
                 int Niter, int K_Aver, Bool SearchPositiv, Bool Verbose) {

    int Nx = Dirty_Map.nx();
    int iteration=0;
    int Max_i;
    float Old_Max, Max_Sig = 0.;
    float Gain,Gain_X;
    float Coef_Mul;
    float Aver_Abs = 0., Aver_Abs_X;
    int Max_Beam_i;
    float Max_Beam, Max_Dirty_Map;
    fltarray Map_Energy(Nx,"map energy");

    Clean_Map.init();
  
    dec_pos_max (Dirty_Beam, Max_Beam_i, Max_Beam, SearchPositiv);
    Gain = Gamma / Max_Beam;

    // energy map creation
    map_energy (Dirty_Map, Tab_Dirty_X, Map_Energy);

    /* Detection du maximum de la map d'energie */
    dec_pos_max (Map_Energy, Max_i, Max_Sig, True);
    Max_Dirty_Map = Dirty_Map(Max_i);

    /* Debut de la boucle */
    Gain = Gamma / Max_Beam;
    Gain_X = Gamma / Max_Beam;
    Old_Max = ABS(Max_Sig) + 1;

    while (  (ABS(Max_Dirty_Map) >= Noise) 
             && (iteration < Niter)
             && (Old_Max - Max_Sig > FLOAT_EPSILON)
           ) {
           
       Coef_Mul = Dirty_Map(Max_i) * Gain;
       Clean_Map(Max_i)  += Coef_Mul;

       // Subtract the dirty beam to dirty map. the dirty beam is centered
       //   on Imax, and weighted by Gamma  
       Aver_Abs_X = substract_signal (Tab_Dirty_X, Dirty_Beam_X, Max_Beam_i, 
                                    Max_i, Coef_Mul);
       Aver_Abs = substract_signal (Dirty_Map, Dirty_Beam, Max_Beam_i, 
                                  Max_i, Coef_Mul);

       iteration ++;
       Old_Max = Max_Sig;
  	
       map_energy (Dirty_Map, Tab_Dirty_X, Map_Energy);
       dec_pos_max (Map_Energy, Max_i, Max_Sig, True);
       Max_Dirty_Map = Dirty_Map(Max_i);
    }

    if (Verbose == True) {
    
        cout << "End loop: number of iterations = " << iteration
             << ", Max_Ima = (" << Max_i << ") : " << Max_Dirty_Map << endl;
        cout << "end loop criterion is : " << endl;                                      
        if (iteration >= Niter) 
           cout << "Nbr_of_iter = MAX_ITER"<< endl;
        if (ABS(Max_Dirty_Map) < Noise ) 
           cout << "Max_Ima < Noise Level"   << endl;
        if (ABS(Old_Max) < ABS(Max_Sig) )  
           cout << "Old_Max > New_Max" << endl;
    }
}


/**************************************************************************/
/**********                     *******************************************/
/**********      dec_clean      *******************************************/
/**********                     *******************************************/
/**************************************************************************/
void dec_clean (fltarray& Signal, fltarray& Psf, fltarray& Signal_Out, 
                float Fwhm, float Noise, float Gamma, int Niter, 
                int K_Aver, Bool AddResi, Bool SearchPositiv, Bool UseEnergy,
		Bool Verbose) {
                
    int Nx = Signal.nx();
    fltarray Clean_Beam(Nx, "Clean_Beam");

    // find the dirac list position
    if (UseEnergy == False)
       dec_hogbom (Signal, Psf, Gamma, Noise, Signal_Out, 
                   Niter, K_Aver, SearchPositiv, Verbose);
    else {
    
      fltarray Signal_x(Nx,"Signal_x");
      fltarray Signal_y(Nx,"Signal_y");
      fltarray Dirty_Beam_X(Nx,"DirtBeam_x");
      fltarray Dirty_Beam_Y(Nx,"DirtBeam_y");
      fft_derive_image (Signal,  Signal_x);
      fft_derive_image (Psf,  Dirty_Beam_X);
      clean_hogbom_energy (Signal, Psf, Gamma, Noise, Signal_Out, Signal_x, 
		           Dirty_Beam_X, Niter, K_Aver, SearchPositiv, Verbose);
    }
    
    // create the clean beam 
    if (Fwhm > FLOAT_EPSILON) {
    
       Clean_Beam = sig_gaussian (Nx, Fwhm);

       // Convolution of the result by the clean beam 
       fft1d_conv (Signal_Out, Clean_Beam, Signal_Out);
    }

    if (AddResi == True) Signal_Out += Signal;
}

/***************************************************************************/

