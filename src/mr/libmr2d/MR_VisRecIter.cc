//----------------------------------------------------------
//		Copyright (C) 1996 OCA-CEA
//----------------------------------------------------------
//    UNIT
//
//    Version: 
//
//    Author:	Frederic Rue & Benoit VANDAME & Jean-Luc Starck
//
//    Date:	05/04/96 
//    
//    File:	MR_VisRecIter.C 
//
//----------------------------------------------------------

#include	"IM_Obj.h"
#include	"MR_Obj.h"
#include	"IM_VisTool.h"
#include	"MR_VisElem.h"
#include	"MR_VisTree.h"
#include	"MR_NoiseModel.h"
#include	"MR_ListObj.h"
#include	"MR_Deconv.h"

//----------------------------------------------------------
//	Calcul_masque
//		Calcul masque Min suivant W0  
//----------------------------------------------------------
static Icomplex_f Psf_cf;

void calcul_masque (MultiResol& W0, MultiResol & Min)
{
        int NScale = W0.nbr_band();
        int i,j,Nb_ech = NScale - 1;

	Min.band(Nb_ech).init (0.0);
	
	for (int n=0; n< Nb_ech; n++)
	{
           int Nls = W0.size_band_nl(n);
           int Ncs = W0.size_band_nc(n);

	   for (i=0 ; i<Nls ; i++)
	   for (j=0 ; j<Ncs ; j++)
           {
               if (W0(n,i,j)) Min(n,i,j) = 1.;
               else Min(n,i,j) = 0.;
           }
	}
}


//----------------------------------------------------------
//	Rec_iter_grad
//		reconstruction iterative utilisant la methode du gradient
//----------------------------------------------------------

void rec_iter_grad (MultiResol& W0, Ifloat &Im_rec, int Nb_iter,
				double& Erreur,float eps)
{
	float w;
	double	Erreur0=0., num, den;
	int s,i,j,n = 0;
        int Nl = W0.size_ima_nl();
        int Nc = W0.size_ima_nc();
        int NScale = W0.nbr_scale();
        type_transform Trans = W0.Type_Transform;
	MultiResol V (Nl, Nc, NScale, Trans, "V rec_iter_grad");
	MultiResol Min(Nl, Nc, NScale, Trans, "min rec_iter_grad");
        Ifloat I0(Nl, Nc, "grad I0");
        Ifloat temp(Nl, Nc, "grad I0");
        Ifloat Ir(Nl, Nc, "Ir I0");

	calcul_masque (W0, Min);
	W0.rec_adjoint (Im_rec, False);

        for (i=0; i< Nl; i++)
        for (j=0; j< Nc; j++)
        {
          Erreur0 += Im_rec(i,j)*Im_rec(i,j);
          temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
        }
	Erreur0 = sqrt (Erreur0);

	I0 = Im_rec;
                
	do {
		n++;

	        // TO IMAGE RECONSTRUITE
		V.transform (temp);
                for (s=0; s < V.nbr_band()-1; s++)
 		for (i=0; i< V.size_band_nl(s); i++) 
		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);		
		
	        //---- CALCUL IMAGE RESIDUE
		// rec_adjoint (V, temp);
		V.rec_adjoint(temp, False);
		Ir = I0 - temp;
         //  Erreur = sqrt(energy(Ir)) / Erreur0;
         //  cout << n << ": Resi = " << Erreur << endl;

		Erreur = 0.;
                for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) Erreur += Ir(i,j)*Ir(i,j);
	        Erreur = sqrt (Erreur) / Erreur0;

	        //---- CALCUL PARAMETRE DE CONVERGENCE
		V.transform (Ir);
                for (s=0; s < V.nbr_band()-1; s++)
                for (i=0; i< V.size_band_nl(s); i++) 
		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);
		V.rec_adjoint(temp, False);

		num = 0.;
		den = 0.;
                for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) 
                {
                   num += Ir(i,j)*temp(i,j);
                   den += temp(i,j)*temp(i,j);
                }
		if(!den) w = 1.0;
		else w = MAX (1.0, num/den);

	        //---- CALCUL NOUVELLE SOLUTION
		
	        for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) 
		{
		    Im_rec(i,j) +=  w * Ir(i,j);
                    temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
                }
	} while(n < Nb_iter && Erreur > eps);

	threshold(Im_rec);
}


//----------------------------------------------------------
//	Rec_iter_grad_conj
//		reconstruction iterative utilisant
//		la methode du gradient conjugue
//----------------------------------------------------------

 
void rec_iter_grad_conj (MultiResol& W0, Ifloat &Im_rec, int Nb_iter,
				double& Erreur, float eps)
{
	float		a, b;
	double		Erreur0, num, den;
	int		s,i,j,n = 0;
        int Nl = W0.size_ima_nl();
        int Nc = W0.size_ima_nc();
        int NScale = W0.nbr_scale();
        type_transform Trans = W0.Type_Transform;
	MultiResol V (Nl, Nc, NScale, Trans, "V rec_iter_grad");
	MultiResol Min(Nl, Nc, NScale, Trans, "min rec_iter_grad");
        Ifloat I0(Nl, Nc, "grad I0");
        Ifloat temp(Nl, Nc, "grad I0");
        Ifloat Ir(Nl, Nc, "Ir I0");
        Ifloat Ierr(Nl, Nc, "Ierr");
        Ifloat Iint1(Nl, Nc, "Iint1");
        Ifloat Iint2(Nl, Nc, "Iint2");

 	calcul_masque (W0, Min);
	W0.rec_adjoint (Im_rec,False);
	// rec_adjoint(W0, Im_rec);
	 
        Erreur0 = 0.;
	for (i=0; i< Nl; i++)
	for (j=0; j< Nc; j++)
        { 
            Erreur0 += Im_rec(i,j)*Im_rec(i,j);
            temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
        }
        Erreur0 = sqrt (Erreur0);

	I0 = Im_rec;
	Ir.init (0.0);
	a = b = 0.0;

	do {
		n++;

	//---- TO IMAGE RECONSTRUITE
		V.transform (temp);
 		for (s=0; s < V.nbr_band()-1; s++)
		for (i=0; i< V.size_band_nl(s); i++) 
		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);

	        //---- CALCUL IMAGE RESIDUE
	        V.rec_adjoint (temp, False);
		// rec_adjoint(V, temp);
		
		for (i=0; i< Nl; i++)
		for (j=0; j< Nc; j++) Ierr(i,j) = I0(i,j) - temp(i,j);
 
                V.transform(Ierr);
                for (s=0; s < V.nbr_band()-1; s++)
                for (i=0; i< V.size_band_nl(s); i++) 
		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);
 
		if (n > 1)
		{
			V.rec_adjoint (Iint1, False);
			// rec_adjoint(V, Iint1);
                        num = den = 0.;
                        for (i=0; i< Nl; i++)
                        for (j=0; j< Nc; j++) 
                        {
                           num += Iint1(i,j)*Iint2(i,j);
                           den += Iint2(i,j)*Iint2(i,j);
                        }
			if (!den) b = 0.0;
			else b = MIN (1, fabs (-num / den));
		}
                Erreur =0.;
                for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++)
                {
                   Erreur += Ierr(i,j)*Ierr(i,j);
                }
	        Erreur = sqrt (Erreur) / Erreur0;

	//---- CALCUL IMAGE RESIDUE
	        for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) 
		{
		    temp(i,j) = b * Ir(i,j);
                    Ir(i,j) = Ierr(i,j) + temp(i,j);
                }

	//---- CALCUL DE A
		V.transform (Ir);
                for (s=0; s < V.nbr_band()-1; s++)
                for (i=0; i< V.size_band_nl(s); i++) 
		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);
 
                V.rec_adjoint (Iint2, False);
                // rec_adjoint(V, Iint2);
                num = den = 0.;
                for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) 
                {
                    num += Ierr(i,j)*Iint2(i,j);
                    den += Iint2(i,j)*Iint2(i,j);
                }
		if (!den) a = 1.0;
		else a = MAX (1.0, num / den);
	
	//---- CALCUL NOUVELLE SOLUTION
	        for (i=0; i< Nl; i++)
                for (j=0; j< Nc; j++) 
		{
		    Im_rec(i,j) += a * Ir(i,j);
                    temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
                }


	//---- AFFICHAGE RESULTATS
	} while(n < Nb_iter && Erreur > eps);

        // cout << n << ": Error = " << Erreur << endl;
	threshold(Im_rec);
}

/***************************************************************/
 
void mr_sm_recons_obj (MultiResol& MR_Data, Ifloat &Imag, int Nb_iter,
		       double& Error, float eps)
{
   int i,j,s,Iter = 0;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int NbrScale = MR_Data.nbr_scale();
   double Error0;
   float Val;
   Ifloat Resi (Nl, Nc, "Resi");
   Ifloat Sol0 (Nl, Nc, "Sol0");
   MRNoiseModel NoiseModel(NOISE_GAUSSIAN, Nl, Nc, 
                           MR_Data.nbr_scale(), MR_Data.Type_Transform);

   // model noise initialization
   for (s = 0; s < NbrScale-1; s++)
   for (i = 0; i < MR_Data.size_scale_nl(s); i++)
   for (j = 0; j < MR_Data.size_scale_nc(s); j++)
   {
       if (ABS(MR_Data(s,i,j)) > FLOAT_EPSILON)  
              NoiseModel.support(s,i,j) = VAL_SupOK;
       else   NoiseModel.support(s,i,j) = VAL_SupNull;
   }
    MR_Data.rec_adjoint (Resi, False);
    // MR_Data.band(NbrScale-1).init();
    // MR_Data.recons (Resi);
    
    Error0 = sqrt(energy(Resi));
    Imag = Resi;
    Sol0 = Resi;
    threshold(Imag);
    do {
          MR_Data.transform (Imag);
          NoiseModel.threshold(MR_Data, False);
          MR_Data.rec_adjoint (Resi, False);
          // MR_Data.band(NbrScale-1).init();
          // MR_Data.recons (Resi);
          Error=0.;
          for (i = 0; i < Nl; i++)
          for (j = 0; j < Nc; j++)
          {
            Val =  Sol0(i,j) - Resi(i,j);
            Error += Val*Val;
            Imag(i,j) += Val;
            if (Imag(i,j) < 0.) Imag(i,j) = 0.;
          }
          Error = sqrt(Error) / Error0;
        Iter ++;
    } while ((Iter < Nb_iter) && (Error > eps));
    // cout << "Iter = " << Iter << "  Error = " << Error << endl;
}

/***************************************************************/

static void im_reduce(Ifloat &Imag, Ifloat &ImagOut, float Zoom)
{
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nlz = (Zoom > 1.) ? (int)(Nl/Zoom): Nl;
   int Ncz = (Zoom > 1.) ? (int)(Nc/Zoom): Nc;
   int i,j,k,l;
   
   if ( (ImagOut.nl() != Nlz) || (ImagOut.nc() != Ncz))
                                        ImagOut.resize(Nlz, Ncz);
   for (i=0;i<Nlz;i++)
   for (j=0;j<Ncz;j++)
   {
      int c=0;
      int Depi=(int)(i*Zoom);
      int Endi=(int)((i+1)*Zoom);
      int Depj=(int)(j*Zoom);
      int Endj=(int)((j+1)*Zoom);
      ImagOut(i,j) = 0.;
      for (k=Depi; k<Endi;k++)
      for (l=Depj; l<Endj;l++)
      {
         if ((k<Nl) && (l<Nc))
	 {
	    ImagOut(i,j) += Imag(k,l);
	    c++;
	 }
      }
      if (c != 0) ImagOut(i,j) /= (float) c;
   }
   // INFO(Imag,"Imag");
   // INFO(ImagOut, "Imag out");
}   
   
/******************************************************************/

static void im_reduce(Ifloat &Imag, float Zoom)
{
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nlz = (Zoom > 1.) ? (int)(Nl/Zoom): Nl;
   int Ncz = (Zoom > 1.) ? (int)(Nc/Zoom): Nc;
   Ifloat ImagOut(Nlz,Ncz,"zoom out");
   im_reduce(Imag, ImagOut, Zoom);
   Imag.resize(Nlz, Ncz);
   Imag = ImagOut;
}

/******************************************************************/

void mr_psf_recons_obj (MultiResol& MR_Data, Ifloat &Result, Ifloat &Psf,
                        int Nb_iter, double& Error, float eps, 
			float Fwhm, Ifloat &IMGauss, 
			Bool Deconv, float Zoom)
// Reconstruct a deconvolved object from its wavelet transform and the
// corresponding PSF.
// If Zoom > 1, the PSF is must be oversampled, when compared to thr data
// i.e. the pixel size in the PSF is smaller than in the data, and in
// this case, we reconstrct an object with the same pixel size as the PSF
{
   int i,j,s,Iter = 0;
   int NbrScale = MR_Data.nbr_band();
   int Nl = MR_Data.size_ima_nl(); // size of the reconstructed image
   int Nc = MR_Data.size_ima_nc(); 
   Ifloat Imag (Nl, Nc, "Imag"); // Imag = Object * PSF
   Ifloat Resi (Nl, Nc, "Residual");
   Ifloat RecIma (Nl, Nc, "RecIma");
 
   MultiResol MRAux(Nl, Nc, MR_Data.nbr_scale(), MR_Data.Type_Transform, "Aux");
   MRAux.Border = MR_Data.Border;
   int Nlo = Result.nl(); // size of the reconstructed object
   int Nco = Result.nc(); // Nlo = Nl * Zoom and Nco = Nc * Zoom
   if (Zoom > 1)
   {
      Nlo = (int) (Nlo*Zoom);
      Nco = (int) (Nco*Zoom);
      Result.resize(Nlo,Nco);
   } 
   Ifloat Temp (Nlo, Nco, "Temp");    
 
   MR_Data.recons(RecIma);
   // io_write_ima_float("xx_rec.fits", RecIma);
   
   // Calculate the Fourier transform of the PSF
   psf_get(Psf, Psf_cf, Nlo, Nco, True);

   // Calculate the starting point solution
   if (Zoom > 1) im_bilinear_interp (RecIma, Result);
   else Result = RecIma;
   float CoefFlux;
   double FluxData = total(RecIma);
   double FluxIma, Error0 = sqrt(energy(Result)); 
   do { 
          // solution in image space
	  if (Zoom > 1) 
	  {
	     psf_convol(Result, Psf_cf, Temp);
	     im_reduce(Temp, Imag, Zoom);
	  }
	  else psf_convol (Result, Psf_cf, Imag);
	  // INFO(Imag, "Imag");
	  // Flux Constraint
// 	  FluxIma=0.;
// 	  for (i = 0; i < Nl; i++)
//           for (j = 0; j < Nc; j++) 
// 	             if (RecIma(i,j) > FLOAT_EPSILON)  FluxIma += Imag(i,j);
//           CoefFlux = (float) (FluxData / FluxIma);
//    	  for (i = 0; i < Nl; i++)
// 	  for (j = 0; j < Nc; j++) Imag(i,j) *= CoefFlux;
// 	  for (i = 0; i < Nlo; i++)
// 	  for (j = 0; j < Nco; j++) Result(i,j) *= CoefFlux;
	  
 	  // Apply the WT to the rec. object reconvolved with the PSF
	  // and extract the difference with data wavelet coefficient
	  // MRAux.transform (Imag);
          MRAux.transform (Imag);
 	  if (MRAux.Type_Transform != TO_PAVE_BSPLINE)
	  {
	     MRAux.band(NbrScale-1).init();
	     MRAux.recons(Resi);
	     FluxIma = total(Resi);
	  }
	  else FluxIma = 0.;
 	  for (s = 0; s < NbrScale-1; s++)
	  for (i = 0; i < MR_Data.size_band_nl(s); i++)
	  for (j = 0; j < MR_Data.size_band_nc(s); j++)
          {
	      if (ABS(MR_Data(s,i,j)) > FLOAT_EPSILON) 
	      {
	            if (MRAux.Type_Transform == TO_PAVE_BSPLINE)
		         FluxIma += MRAux(s,i,j);
 	            MRAux(s,i,j) = MR_Data(s,i,j) - MRAux(s,i,j);
	      }
	      else  MRAux(s,i,j) = 0.;
              // if (ABS(MR_Data(s,i,j)) < FLOAT_EPSILON) MRAux(s,i,j) = 0.;
	      // else FluxIma += MRAux(s,i,j);
          }
	  MRAux.rec_adjoint(Resi, False);
     	  //MRAux.rec_adjoint(Imag, False);
 	  // MRAux.band(NbrScale-1).init();
	  // MRAux.recons(Resi);
	  
	  // Flux Constraint
	  if (Iter > 3)
	  {
              // CoefFlux = (float) (FluxData / total(Imag));
	      CoefFlux = (float) (FluxData / FluxIma);
    	      for (i = 0; i < Nl; i++)
	      for (j = 0; j < Nc; j++) Imag(i,j) *= CoefFlux;
 	      for (i = 0; i < Nlo; i++)
	      for (j = 0; j < Nco; j++) Result(i,j) *= CoefFlux;
           }
	   // for (i = 0; i < Nl; i++)
	   //for (j = 0; j < Nc; j++) Resi(i,j) = RecIma(i,j) - Imag(i,j);  
     	  if (Zoom > 1)
          {
	     im_bilinear_interp (Resi, Temp);
 	     Result += Temp;
	     threshold(Result);
	     for (i = 0; i < Nl; i++)
	     for (j = 0; j < Nc; j++) Error += Resi(i,j) * Resi(i,j);
  	  }
	  else 
	  {
  	      Error = 0.;
              for (i = 0; i < Nl; i++)
              for (j = 0; j < Nc; j++) 
	      {
	         Result(i,j) += Resi(i,j);
   	         if (Result(i,j) < 0.) Result(i,j) = 0.;
 		 Error += Resi(i,j) * Resi(i,j);
	      }
  	  }
          Error = sqrt(Error);
	  
          Iter ++;
//  	  cout << Iter << " Error = " << Error << " Flux = " << flux(Result) <<
//	      " Sigma = " << sigma(Resi) << " Min = " << min(Result) << " Max = " << max(Result) << " CoefFlux  = " << CoefFlux << endl;
       } while ((Iter < Nb_iter) && (Error/Error0 > eps));
       
//       io_write_ima_float("xx_data.fits", RecIma);
//       io_write_ima_float("xx_sol.fits", Result);
//       io_write_ima_float("xx_resi.fits", Resi); 
//       MR_Data.write("xx.mr"); 
//       MRAux.write("xx1.mr"); 

     psf_convol (Result, Psf_cf, Imag);
//     io_write_ima_float("xx_ima.fits", Imag);
//     io_write_ima_float("xx_psf.fits", Psf);
//     cout << "Nl psf = " << Psf_cf.nl() << "  Nc psf = " << Psf_cf.nc() << endl;
//     cout << "Nl Data " <<  Result.nl()  << "  Nc Data = " << Result.nc() << endl;
//     cout << "psf = " << Psf_cf(Psf_cf.nl()/2, Psf_cf.nc()/2).real() << "  " << Psf_cf(Psf_cf.nl()/2, Psf_cf.nc()/2).imag()  << endl;
//     cout << "psf = " << Psf_cf(Psf_cf.nl()/2-1, Psf_cf.nc()/2-1).real() << "  " << Psf_cf(Psf_cf.nl()/2, Psf_cf.nc()/2).imag()  << endl;
//     cout << "Iter = " << Iter << "  Error = " << Error << endl;

}



/***************************************************************/

#define SQR(a) (a)*(a)
void xmm_get_psf (float x, float y,    /* target x,y position in pixels */
	          Ifloat &Image)       /* PSF output image */
  /* calculates XMM off-axis PSF (normalized surface brightness)
   * The XMM QM PSF has been fit  with 4 Gaussians , and normalized, 
   * so that  Integral 2*PI*r*dr*f(r) from 0 to infinity = 1 [1/arcsec2]
   */
{
  int nx,ny,nqx,nqy;
  double sum=0.0,sig2, argum1,argum2, sig_blur;
  double th2, sscale, mscale, xc, yc, rr2, xqc, yqc; 
  double offsg=0.0324; /* PSF blurring */
    
    /* following are the norms after integrating each component to 1.00 */
  float norm[] = {0.4169, 0.3470, 0.1279, 0.1082};
  float sigma[] = {4.1450, 8.6218, 22.372, 66.698};
  
  /* Assuming 512x512 image and scale of 4"/px */
  /* NOTE! It is XMM-pn specific */
  nx = 512;
  ny = 512;

  sscale = 4.0;       /* scale in "/pixel */
  mscale = sscale/60.0; /* scale in '/pixel */
  
  /* centre of the field of view, assuming in the centre of the pixel */
  xc = nx/2.0;
  yc = ny/2.0;
  
  /* Off-axis angle in arcmin */
  th2 = mscale*mscale*(SQR(x-xc) + SQR(y-yc));
  
  sig_blur=offsg*th2;
  
  /* For PSF image we don't need more than 51x51 pixels */
  /* I should make it adaptive!!! */
  nqx = XMM_PSF_NL;
  nqy = XMM_PSF_NC;
  Image.resize(nqy,nqx);
  xqc = nqx/2;
  yqc = nqy/2;
  
  double Flux=0.;
  int i,j;
  for (i=0; i< nqx; i++)
  for (j=0; j< nqy; j++)
  {
     /* distance from the source position in arcsec */
     rr2 = sscale*sscale*(SQR(i-xqc) + SQR(j-yqc));
     sum = 0;
     for(int k=0; k<4; k++)
     {
	sig2=sigma[k]*sigma[k]+sig_blur*sig_blur;
	if (sig2 > 0) 
	{
	   argum1 = rr2/(2.0*sig2);
	   argum2 = norm[k]/(2.0*M_PI*sig2);
	   sum += argum2*exp(-argum1);
	}
     }
     if (sum < 1e-4) sum = 0.;
     Image(j,i) = (float) sum;
     Flux += sum;
   }
   for (i=0; i< nqx; i++)
   for (j=0; j< nqy; j++) Image(j,i) = (float) (Image(j,i) / Flux); 
   //cout << "PSF " << x << " " << y << " Flux = " << Flux << endl;
   //INFO(Image, "PSF");
   //io_write_ima_float("xx_psf.fits", Image);
}

/*********************************************************************/


void ListObj2D::recons_obj(MultiResol& W, Ifloat &Im_rec, double & Erreur, 
                          double & Flux, t_Gauss_2D & Gauss, 
			  float PosX, float PosY)
{
   double X,Y;
   
   if ((W.Type_Transform == TM_TO_PYR) || (W.Type_Transform == TM_TO_SEMI_PYR))
   {
      if (ReconsMethod == GRAD_PSF)
      {
         cerr << "Error: the PSF cannot be used for the reconstruction " << endl;
         cerr << "       with this transform ... " << endl;
         exit(-1);
      }
      // W.recons(Im_rec);
      Erreur = 0.;
   }
   else switch(ReconsMethod)
   {
	case GRAD_OPTI_STEP: 
             rec_iter_grad (W, Im_rec, Nb_iter_rec, Erreur, ErrorRec);
              break;	
	case GRAD_CONJUG: 
             rec_iter_grad_conj (W, Im_rec, Nb_iter_rec, Erreur, ErrorRec);
             break;
	case GRAD_FIX_STEP: 
             mr_sm_recons_obj(W, Im_rec, Nb_iter_rec, Erreur, ErrorRec);
             break;
	case GRAD_PSF:
	     if (Psf.n_elem() == 0)
	     { 
	        cerr << "Error: the PSF is not defined... " << endl;
	        exit(-1);
	     }
             mr_psf_recons_obj(W, Im_rec, Psf, 
			       Nb_iter_rec, Erreur, ErrorRec, 
			       Fwhm, IMGauss,  Deconv, ZoomPSF);
	    // if Fwhm > 0 the resolutiuon is limited
	    if (Fwhm > 0) psf_convol (Im_rec, IMGauss);
	    break;
	case GRAD_PSF_XMM:
	     if (Psf.n_elem() == 0) Psf.alloc(XMM_PSF_NL,XMM_PSF_NC,"XMM PSF");
 	     X =  PosX - Nc / 2.  + 256; 
	     Y =  PosY - Nl / 2.  + 256;  
	     // cout << "XMM PSF Pos = " << X << " " << Y << endl;
             xmm_get_psf (X, Y, Psf);	    
	     // io_write_ima_float("xx_psf.fits", Psf);
	     if ((IMGauss.nl() == Psf.nl()) || (IMGauss.nc() == Psf.nc()))
	     {
	        IMGauss = im_gaussian(Psf.nl(), Psf.nc(), Fwhm*ZoomPSF);
		norm_flux(IMGauss);
                psf_convol (Psf, IMGauss);
             }
 
 	     // rec_psf_iter_grad(W, Im_rec, Psf, Nb_iter_rec, Erreur, ErrorRec);
             mr_psf_recons_obj(W, Im_rec, Psf, 
			       Nb_iter_rec, Erreur, ErrorRec, 
			       Fwhm, IMGauss,  Deconv, ZoomPSF);
 	    // if Fwhm > 0 the resolutiuon is limited
	    if (Fwhm > 0) psf_convol (Im_rec, IMGauss);
	    break;
	default: 
             cerr << "Error: unknown reconstruction method ... "<<endl;
	     exit(-1);
	     break;
   }

   // find Gaussian parameters 
   Estime_param_gauss_2D(Gauss, Im_rec);
	
   if ((ReconsMethod == GRAD_PSF) || (ReconsMethod == GRAD_PSF_XMM))
   {
 	// Test if a deconvolution is wanted
        if (Deconv == False)
	{
	    // cout << "Reconv " << endl;
            // psf_convol (Im_rec, Psf);
	    psf_convol (Im_rec, Psf_cf);
            threshold(Im_rec);
	}
	// Test if the deconvolved image is oversampled
	if (ZoomPSF > 1)
	{
	    im_reduce(Im_rec, ZoomPSF);
	}
   }
   Gauss.sX /= ZoomPSF;
   Gauss.sY /= ZoomPSF;
   Gauss.mx /= ZoomPSF;
   Gauss.my /= ZoomPSF;
   // Gauss.amp = ValPixMax
	
   Flux = total(Im_rec);
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
// 	  // Flux Constraint
// 	  FluxIma=0.;
// 	  for (i = 0; i < Nl; i++)
//           for (j = 0; j < Nc; j++) 
// 	             if (RecIma(i,j) > FLOAT_EPSILON)  FluxIma += Imag(i,j);
//           CoefFlux = (float) (FluxData / FluxIma);
//    	  for (i = 0; i < Nl; i++)
// 	  for (j = 0; j < Nc; j++) Imag(i,j) *= CoefFlux;
// 	  for (i = 0; i < Nlo; i++)
// 	  for (j = 0; j < Nco; j++) Result(i,j) *= CoefFlux;
	  
	  
	  
	  
// void rec_psf_iter_grad (MultiResol& W0, Ifloat &Im_rec, Ifloat &Psf, int Nb_iter,
// 				double& Erreur,float eps)
// {
// 	float		w;
// 	double		Erreur0=0., num, den;
// 	int		s,i,j,n = 0;
//         int Nl = W0.size_ima_nl();
//         int Nc = W0.size_ima_nc();
//         int NScale = W0.nbr_scale();
//         type_transform Trans = W0.Type_Transform;
// 	MultiResol V (Nl, Nc, NScale, Trans, "V rec_iter_grad");
// 	MultiResol Min(Nl, Nc, NScale, Trans, "min rec_iter_grad");
//         Ifloat I0(Nl, Nc, "grad I0");
//         Ifloat temp(Nl, Nc, "grad I0");
//         Ifloat Ir(Nl, Nc, "Ir I0");
//      
//         Icomplex_f Psf_cf;
//         psf_get(Psf, Psf_cf, Nl, Nc, True);
//    
// 	calcul_masque (W0, Min);
// 	W0.rec_adjoint (Im_rec, False);
//  
//         for (i=0; i< Nl; i++)
//         for (j=0; j< Nc; j++)
//         {
//           Erreur0 += Im_rec(i,j)*Im_rec(i,j);
//           temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
//         }
// 	Erreur0 = sqrt (Erreur0);
// 
// 	I0 = Im_rec;
//                 
// 	do {
// 		n++;
// 
// 	        // TO IMAGE RECONSTRUITE
// 		psf_convol(temp, Psf_cf);
// 		V.transform (temp);
//                 for (s=0; s < V.nbr_band()-1; s++)
// 		for (i=0; i< V.size_band_nl(s); i++)
// 		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);
// 		
//  	        //---- CALCUL IMAGE RESIDUE
//  		V.rec_adjoint(temp, False);
// 		Ir = I0 - temp;
//  		
// 		Erreur = 0.;
//                 for (i=0; i< Nl; i++)
//                 for (j=0; j< Nc; j++) Erreur += Ir(i,j)*Ir(i,j);
// 	        Erreur = sqrt (Erreur) / Erreur0;
// 
// 	        //---- CALCUL PARAMETRE DE CONVERGENCE
//  		V.transform (Ir);
//                 for (s=0; s < V.nbr_band()-1; s++)
// 		for (i=0; i< V.size_band_nl(s); i++) 
// 		for (j=0; j< V.size_band_nc(s); j++) V(s,i,j) *= Min(s,i,j);
// 		V.rec_adjoint(temp, False);
//  		
// 		num = 0.;
// 		den = 0.;
//                 for (i=0; i< Nl; i++)
//                 for (j=0; j< Nc; j++) 
//                 {
//                    num += Ir(i,j)*temp(i,j);
//                    den += temp(i,j)*temp(i,j);
//                 }
// 		if(!den) w = 1.0;
// 		else w = MAX (1.0, num/den);
// 
// 	        //---- CALCUL NOUVELLE SOLUTION
// 		
// 	        for (i=0; i< Nl; i++)
//                 for (j=0; j< Nc; j++) 
// 		{
// 		    Im_rec(i,j) +=  w * Ir(i,j);
//                     temp(i,j) = (Im_rec(i,j) > 0.) ? Im_rec(i,j) : 0.;
//                 }
//   cout << n << " w = " << w <<  " Error = " <<  Erreur << " Flux = " << flux(temp) << " Max = " << max(temp) << endl;
// 
// 	} while(n < Nb_iter && Erreur > eps);
// 	
//   	threshold(Im_rec);
//        io_write_ima_float("xx_sol.fits", Im_rec);
// 
// }
 





//----------------------------------------------------------
//	Rec_iter
//		Reconstruction image a partir de sa T.O
//----------------------------------------------------------

// void mr_recons_obj (MultiResol& W, Ifloat &Im_rec, type_objrec_method Meth, 
//                     int Nb_iter, double& Erreur, float eps, Ifloat & Psf, 
// 		    Bool Deconv, float Alpha, float Zoom)
// {
//         if ((W.Type_Transform == TM_TO_PYR) || (W.Type_Transform == TM_TO_SEMI_PYR))
//         {
// 	    if (Meth == GRAD_PSF)
// 	    {
// 	       cerr << "Error: the PSF cannot be used for the reconstruction " << endl;
// 	       cerr << "       with this transform ... " << endl;
// 	       exit(-1);
// 	    }
//             W.recons(Im_rec);
//             Erreur = 0.;
//         }
//         else
// 	switch(Meth)
// 	{
// 		case GRAD_OPTI_STEP: 
//                         rec_iter_grad (W, Im_rec, Nb_iter, Erreur, eps);
//                         break;	
// 		case GRAD_CONJUG: 
//                         rec_iter_grad_conj (W, Im_rec, Nb_iter, Erreur, eps);
//                         break;
// 		case GRAD_FIX_STEP: 
//                         mr_sm_recons_obj(W, Im_rec, Nb_iter, Erreur, eps);
//                         break;
// 		case GRAD_PSF: 
// 		        if (Psf.n_elem() == 0)
// 			{
// 			   cerr << "Error: the PSF is not defined... " << endl;
// 	                   exit(-1);
// 	                }
// 			rec_psf_iter_grad_conj(W, Im_rec, Psf, Nb_iter, Erreur, eps);
//                         //mr_psf_recons_obj(W, Im_rec, Psf, 
// 			//                  Nb_iter, Erreur, eps, Deconv,
// 			//		  Alpha, Zoom);
//                         break;
// 		default: 
//                       cerr << "Error: unknown reconstruction method ... "<<endl;
// 		      exit(-1);
// 		      break;
//  	}
// }



/******************************************************************/


// void test_mr_psf_recons_obj (MultiResol& MR_Data, Ifloat &Result, Ifloat &Psf,
//                         int Nb_iter, double& Error, float eps, Bool Deconv,
// 			float Alpha, float Zoom)
// {
//    int i,j,s,Iter = 0;
//    int Nl = Result.nl();
//    int Nc = Result.nc();
//    int NbrScale = MR_Data.nbr_band();
//    double Error0;
//    int Nlz = (Zoom > 1.) ? (int)(Nl/Zoom): Nl;
//    int Ncz = (Zoom > 1.) ? (int)(Nc/Zoom): Nc;
//    Ifloat Resi (Nlz, Ncz, "Residual");
//    Ifloat Data (Nl, Nc, "Residual");
//    Ifloat Imag (Nl, Nc, "Residual");
//    Ifloat Temp (Nl, Nc, "Temp");
// 
//     cout << " TEST RECONS OBJ " << endl;
//     mr_sm_recons_obj(MR_Data, Data, 10, Error, eps);
//       
//     io_write_ima_float("xx_data.fits", Result);
//    MultiResol MRAux(Nlz, Ncz, MR_Data.nbr_scale(), MR_Data.Type_Transform, "Aux");
// 
//     Icomplex_f Psf_cf;
//     psf_get(Psf, Psf_cf, Nl, Nc, True);
//     //INFO(Psf,"PSF");
//     //cout << "PP = " << Psf_cf(Nl/2,Nc/2).real() << endl;
//     io_write_ima_float("xx_psf.fits", Psf);
//  
//     Result = Data;
//     float CoefFlux, FluxIma,FluxData = flux(Data);
//     Error0 = sqrt(energy(Result)); 
//     INFO(Result, "INIT");
//     MR_Data.write("xx.mr");
//     do { 
//           // solution in image space
//           psf_convol (Result, Psf_cf, Imag);
//  	  io_write_ima_float("xx_ima.fits", Result);
// 	  // Flux Constraint
// 	  FluxIma=0.;
// 	  for (i = 0; i < Nl; i++)
//           for (j = 0; j < Nc; j++) 
// 	             if (Data(i,j) > FLOAT_EPSILON)  FluxIma += Imag(i,j);
//           CoefFlux = FluxData / FluxIma;
//  	  for (i = 0; i < Nl; i++)
// 	  for (j = 0; j < Nc; j++)
// 	  {
// 	     Imag(i,j) *= CoefFlux;
// 	     Result(i,j) *= CoefFlux;
// 	  }
// 	  	  
//           MRAux.transform (Imag);
//   	  for (s = 0; s < NbrScale-1; s++)
//  	  for (i = 0; i < MR_Data.size_band_nl(s); i++)
//  	  for (j = 0; j < MR_Data.size_band_nc(s); j++)
//           {
// 	      if (ABS(MR_Data(s,i,j)) > FLOAT_EPSILON) 
//  	           MRAux(s,i,j) = MR_Data(s,i,j) - MRAux(s,i,j);
// 	      else  MRAux(s,i,j) = 0.;
//               // if (ABS(MR_Data(s,i,j)) < FLOAT_EPSILON) MRAux(s,i,j) = 0.;
// 	      
//           }
// 	  MRAux.band(NbrScale-1).init();
//  	  
//   	  MRAux.rec_adjoint(Resi);      
// 	  // MRAux.recons(Resi);   
//           //INFO(Data-Resi,"Resi");
//           // Resi = Imag;
// 
//  	  Error = 0.;
//           for (i = 0; i < Nl; i++)
//           for (j = 0; j < Nc; j++)
// 	  {
// 	      // Resi(i,j) = Data(i,j) - Imag(i,j);
// 	      Result(i,j) += Resi(i,j);
//  	      Error += Resi(i,j) * Resi(i,j);
// 	      // if (Data(i,j) == 0) Result(i,j) = 0.; 
//  	      if (Result(i,j) < 0.) Result(i,j) = 0.;
//  	  }
// 	  // io_write_ima_float("xx_resi.fits", Resi);
//  
//           Error = sqrt(Error) / Error0;
// 	  
//           Iter ++;
// 	  cout << Iter << " Error = " << Error << " Flux = " << flux(Result) <<
// 	      " Sigma = " << sigma(Resi) << " Min = " << min(Result) << " Max = " << max(Result) << endl;
//     } while ((Iter < Nb_iter) && (Error > eps));
//           
//     INFO(Result, "BEF CONV");
//     psf_convol (Result, Psf_cf, Imag);
//     INFO(Imag, "AF CONV");
//   
//     // cout << "Iter = " << Iter << "  Error = " << Error << endl;
//     io_write_ima_float("xx_sol.fits", Result);
//     io_write_ima_float("xx_resi.fits", Resi);
// 
// }
