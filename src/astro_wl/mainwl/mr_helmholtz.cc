/******************************************************************************
**                   Copyright (C) 2007 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Erwan Deriaz & Jean-Luc Starck
**
**    Date:  08/09/08
**    
**    File:  mr_helmholtz.cc
**
*******************************************************************************
**
**    DESCRIPTION  Helmholtz decomposition program  
**    ----------- 
**                 
******************************************************************************/

#include "MR_DivCurl.h"
#include "FFTN_2D.h"

char Name_Imag_In[512]; /* input data file   */
char Name_Imag_Out[512]; /* output data file */
char Name_Imag_Out_2[512]; /* output data file */
char Name_EB_FileOut[512];
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;   // Verbose mode
int NbrScale=0;       // Number of scales, default is automatically estimated

Bool V0Proj = True;
int NbrIter = DEF_ITER_HELMHOLTZ;  // Number of iterations for div/rot mode decomposition (i.e. helmholtz)
Bool EBMod = False;
type_helmholtz THelm = DEF_HELMHOLTZ;
Bool ReadCat = False;
float BinCat = 1;

type_constraints_helmholtz THCons = C_NO;
Bool UseConstraint_on_B = False;
Bool UseConstraint_on_E = False;
Bool WaveletBorder = False;
int ImageSize=0;
/***************************************/

const char * StringHelmholtzTransform (type_helmholtz TypeHel)
{
    switch (TypeHel)
    {
        case H_ANI_WT: 
			return ("Anisotropic Wavelet Decomposition");break;
        case H_ISO_WT: 
			return ("Isotropic Wavelet Decomposition");break;
        case H_FFT: 
			return ("Fourier Transform");break;
		default:
			return ("Undefined decomposition");
			break;
    }
}

const char * StringHelmholtzConstraint(type_constraints_helmholtz TypeHelConstraint)
{
    switch (TypeHelConstraint)
    {
        case C_NO: 
			return ("No constraint.");break;
        case C_HARD: 
			return ("l0 norm constraint.");break;
        case C_SOFT: 
			return ("l1 norm constraint.");break;
        case C_ZERO_BORDER: 
			return ("Zero border constraint.");break;
		case C_ZERO: 
			return ("Component set to zero.");break;
        default:
			return ("Undefined constraint");
			break;
    }
}

/***************************************/

void Helmholtz::alloc(type_helmholtz TypeHel, int Nx, int Ny)
{
	int NbrScale=0;
	int NbrUndec=0;
 	// cout << "Alloc in " << endl;
	MAD = NULL;
	MID = NULL;
	DC = NULL;
    switch (TypeHel)
    {
        case H_ANI_WT: 
			MAD = new MR_ANI_DIVCURL;
			// MAD->Verbose = Verbose;
			(MAD->alloc)(Ny, Nx, NbrScale);
            MAD->BorderWavelet=BorderWavelet;
			DC = MAD;
			break;
        case H_ISO_WT: 
            MID = new MR_ISO_DIVCURL;
			MID->Verbose = Verbose;
            (MID->alloc)(Ny, Nx, NbrUndec, NbrScale);
  			DC = MID;	
			break;
  		case H_FFT: 
			DC = NULL;
			break;
		default:
			cout << "Error: Undefined decomposition" << endl;
			exit(-1);
			break;
    }
	THelm=TypeHel;
 	// cout << "Alloc out  " << endl;
}

/***************************************/

void Helmholtz::alloc(type_helmholtz TypeHel, HelmCatalog & Cat)
{
    HCat = &Cat;
    UseCat = True;
    int Nx=HCat->nx();
    int Ny=HCat->ny();
    alloc(TypeHel, Nx, Ny);
}

/***************************************/

void infocf(string Mes, Icomplex_d &Image)
{
   string Mr = Mes + " Re: ";
   Ifloat Buff(Image.nl(), Image.nc(), "buff");
   real(Buff, Image);
   Buff.info(Mr);
   string Mi = Mes + " Im: ";
   imag(Buff, Image);
   Buff.info(Mi);
}

/***************************************/

void Helmholtz::fft_trans(fltarray &DataIn)
{
	// E =  1/ k^2 (k1^2 - k2^2) G1  + 1/ k^2 ( 2 k1 k2    ) G2
	// B =  1/ k^2 ( 2 k1 k2   ) G1  + 1/ k^2 ( k2^2 - k1^2) G2
	
	// (u1) = 1/ k^2 [  k2^2  -k1k2  ] (G1)  
	// (u2)          [  -k1k2  k1^2  ]  (G2)
	// 
	// (v1 ) = 1/ k^2 [ k1^2  k1k2  ]  (G1)  
	// (v2)           [ k1k2  k2^2  ]  (G2)
	fltarray Data;
	Data = DataIn;
	int Nl = Data.ny();
	int Nc = Data.nx();
	FFTN_2D FFT;
 	Icomplex_d G1,G2;
	Icomplex_d U1_cf, V1_cf, U2_cf, V2_cf;
	G1.alloc(Nl,Nc,"Buff1");
	G2.alloc(Nl,Nc,"Buff2");
	U1_cf.alloc(Nl,Nc,"Buff1");
	V1_cf.alloc(Nl,Nc,"Buff2");	
	U2_cf.alloc(Nl,Nc,"Buff1");
	V2_cf.alloc(Nl,Nc,"Buff2");
	
	Ifloat Imag;
    float *BuffData = Data.buffer();
    Imag.alloc (BuffData,Nl,Nc);
	FFT.fftn2d(Imag, G1);
	
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	FFT.fftn2d(Imag, G2);
	// Ifloat Bre(Nl,Nc,"rr");
	// Ifloat Bim(Nl,Nc,"ii");
	
	for (int i=0; i < Nl; i++)
    for (int j=0; j < Nc; j++)
	{
		if ((i == Nl/2) && (j == Nc/2))
		{	
			G1(i,j) = G2(i,j) = U1_cf(i,j) = U2_cf(i,j) = V1_cf(i,j) = V2_cf(i,j) = complex_f(0.,0.);
 		}
		else 
		{
		   // float kx = (float) j - Nc/2.;
		   // float ky = (float) i - Nl/2.;
		   double kx =   j - (int) (Nc/2);
		   double ky =   i - (int) (Nl/2);
		   double k2 = kx*kx + ky*ky;
		   double kx2 =  kx*kx / k2;
		   double ky2 =  ky*ky / k2;
		   double kxy =  kx*ky / k2;
		   double nkxy = -kxy;
 
 	        U1_cf(i,j) =  ky2 * G1(i,j) +  nkxy * G2(i,j);
		    U2_cf(i,j) =   nkxy * G1(i,j) +  kx2 * G2(i,j);	

		    V1_cf(i,j) =  kx2 * G1(i,j) +  kxy * G2(i,j);
 	        V2_cf(i,j) =   kxy * G1(i,j) +  ky2 * G2(i,j);
 			
			// G1(i,j) -= (U1_cf(i,j)+V1_cf(i,j));
 			// G2(i,j) -= (U2_cf(i,j)+V2_cf(i,j));
			// Bre(i,j) = G2(i,j).real();
			// Bim(i,j) = G2(i,j).imag();
		}
     }
	
 	// Bre.info();
	// io_write_ima_float("bre", Bre);
 	// Bim.info();
	// io_write_ima_cd("xx_u1", U1_cf);
	
    // G1 -= (U1_cf + V1_cf);
    // infocf("G1", G1);
	// G2 -= (U2_cf + V2_cf);
    // infocf("G2", G2);
	// io_write_ima_cd("xx_g2", G2);
	
	FFT.fftn2d (U1_cf, True);
    FFT.fftn2d (U2_cf, True);			
	FFT.fftn2d (V1_cf, True);
    FFT.fftn2d (V2_cf, True);
/*	io_write_ima_cd("xx_u1", U1_cf);
	io_write_ima_cd("xx_u2", U2_cf);
	io_write_ima_cd("xx_v1", V1_cf);
	io_write_ima_cd("xx_v2", V2_cf);
    infocf("U1_cf", U1_cf);
	infocf("U2_cf", U2_cf);
	infocf("V1_cf", V1_cf);
	infocf("V2_cf", V2_cf);

	for (int i=0; i < Nl; i++)
    for (int j=0; j < Nc; j++) 
	{		
		G1(i,j) = complex_d(Data(j,i,0)) - (U1_cf(i,j)+V1_cf(i,j));
		G2(i,j) = complex_d(Data(j,i,1)) - (U2_cf(i,j)+V2_cf(i,j));
	}	
	infocf("G1dir", G1);
    infocf("G2dir", G2);*/
							
    U.resize(Nc,Nl, 2);
    V.resize(Nc,Nl, 2);
	
 	for (int i=0; i < Nl; i++)
    for (int j=0; j < Nc; j++)
    {
		U(j,i,0) = (U1_cf(i,j)).real();
		U(j,i,1) = (U2_cf(i,j)).real();
		V(j,i,0) = (V1_cf(i,j)).real();
		V(j,i,1) = (V2_cf(i,j)).real();
		// Data(j,i,0) -= (U(j,i,0) + V(j,i,0));
		// Data(j,i,1) -= (U(j,i,1) + V(j,i,1));		
	}
	// U.info("U");
	// V.info("V");
 	// Data -= U;
	// Data -= V;
	// Data.info("datain");
}
									
/***************************************/

void set_to_0_border_wavelet(int N, float *tab)
{
    for (int j1=N;j1>0;j1=j1) 
    { 
        j1=j1/2;
        for (int j2=N;j2>0;j2=j2) 
        {
            j2=j2/2;
            for (int i1=0;i1<(j1+(j1<=1));i1++) 
            {
                for (int i2=0;i2<(j2+(j2<=1));i2++) 
                {
                    // if !((i1==0)||(i1+1==(j1+(j1<=1)))||(i2==0)||(i2+1==(j2+(j2<=1))))
                    if ((i1!=0)&& (i1+1!=(j1+(j1<=1)))&&(i2!=0)&&(i2+1!=(j2+(j2<=1))))
                        tab[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=0.;
                }
            }
        }
    }    
}

/***************************************/

void Helmholtz::iter_trans(fltarray &Data)
{
	fltarray R,Un;
	fltarray Vn;
	fltarray Resi;
    fltarray Trans;
	U.resize(Data.nx(), Data.ny(), Data.nz());
	V.resize(Data.nx(), Data.ny(), Data.nz());
	Un.resize(Data.nx(), Data.ny(), Data.nz());
	Vn.resize(Data.nx(), Data.ny(), Data.nz());
	R.resize(Data.nx(), Data.ny(), Data.nz());
	float LambdaMax,Lambda;
    Bool BordZeroConstraint = False;
    Bool Use_SoftThreshold = False;
    Bool Use_HardThreshold = False;
    Bool UseConstraint;         

    if (THCons != C_NO) UseConstraint = True;
    else UseConstraint = False;
    if ((BorderWavelet == True) && (THCons == C_ZERO_BORDER)) BordZeroConstraint = True;
    if (THCons == C_HARD) Use_HardThreshold = True;
    if (THCons == C_SOFT) Use_SoftThreshold = True;
    // Bool UseConstraint_on_B = True;   // if true,  All common modes in E mode
    // Bool UseConstraint_on_E =  (UseConstraint_on_B == True) ? False: True; // if true,  All common modes in B mode
 
       
 
	// (DC->RotDiv_Trans).alloc(Data.nx(), Data.ny(), Data.nz());
    if (Verbose == True)
    {
        if (BorderWavelet == True) cout << " ITER Wavelet Helmholtz with Border-wavelets " << endl;
        else cout << " ITER Wavelet Helmholtz " << endl;
        if (UseConstraint_on_B == True)  cout << "      : Use constraint on B mode." << endl;
        if (UseConstraint_on_E == True)  cout << "      : Use constraint on E mode." << endl;
        if (BordZeroConstraint == True)  cout << "      : Use Border zero contraints." << endl;
        if (Use_SoftThreshold == True)  cout << "      : l1 constraint (soft thresholding)." << endl;
        if (Use_HardThreshold == True)  cout << "      : l0 constraint (hard thresholding)." << endl;
    }
    
  /*  FFT initialization
    fft_trans(Data);
    Un = U;
    Vn = V;
    for (int x=0; x < Un.nx(); x++)
    for (int y=0; y < Un.ny(); y++)    Un(x,y,1) = Vn(x,y,1) = 0; */
             
  	for (int i=0; i < Niter; i++)
	{
     
		R = Data - (Un+Vn);
		
		// Divergence free decomposition
		DC->transform(R);
        
        if (i == 0) 
        {
           Trans.alloc( (DC->RotDiv_Trans).nx(),  (DC->RotDiv_Trans).ny(), 2);
           float Max=(DC->RotDiv_Trans).maxfabs();
           LambdaMax = ABS(Max) - 0.1 * ABS(Max);
           cout << " LambdaMax = " << LambdaMax << endl;
        }
        Lambda = LambdaMax * (1.  -  (float) ( i) / (float) ( Niter - 1));
        
        if (Verbose == True) cout << "Iter " << i+1  << ", Lamda = " << Lambda  << ", SigmaResi = " << R.sigma() << " " << (DC->RotDiv_Trans).nx() << " " << (DC->RotDiv_Trans).ny() << endl;

		// Set to zero the second component
		for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
		for (int y=0; y < (DC->RotDiv_Trans).ny(); y++)
        {
            (DC->RotDiv_Trans)(x,y,1) = 0;
             Trans (x,y,0) += (DC->RotDiv_Trans)(x,y,0);
        }

		// Reconstruction
		DC->recons(R);
		Un += R;
        if (UseConstraint == True)
        {
            // cout << "B mode cst " << endl;
            fltarray EB;
            EB = Vn;
            EB -= Un;
            // We keep only E mode to make the penalization term
            if ((UseConstraint_on_E == True) && (UseConstraint_on_B == False))
            {
               for (int x=0; x < EB.nx(); x++)
               for (int y=0; y < EB.ny(); y++) EB(x,y,1)= 0.;
            }
            // We keep only B mode to make the penalization term
            // U is div free, and should not contain any B
            if ((UseConstraint_on_B == True) && (UseConstraint_on_E == False))
            {
                // We put E to zero in order to not touch it
                for (int x=0; x < EB.nx(); x++)
                for (int y=0; y < EB.ny(); y++) EB(x,y,0)= 0.;
            }
 
            DC->transform(EB);
            // Div free WT of EB
            for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
            for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) 
            { 
                // We keep the part to subtract to the U mode
                if (Use_SoftThreshold == True) (DC->RotDiv_Trans)(x,y,0) = (DC->RotDiv_Trans)(x,y,0) - soft_threshold( (DC->RotDiv_Trans)(x,y,0),  Lambda);
                else  if (Use_HardThreshold == True) (DC->RotDiv_Trans)(x,y,0) = (DC->RotDiv_Trans)(x,y,0) - hard_threshold( (DC->RotDiv_Trans)(x,y,0),  Lambda);
                // We set to zero the complement part to obtain the divergence null field
                (DC->RotDiv_Trans)(x,y,1) = 0;
            }
            if (BordZeroConstraint == True)
            {
               float *tab = (DC->RotDiv_Trans).buffer();
               int N=(DC->RotDiv_Trans).nx()-2;
               set_to_0_border_wavelet(N,tab); // In fact, set to zero all coefficients not affected by the border
            }   
            DC->recons(R);
            if (Verbose == True) cout << "   DIV BMODE: sig = " << R.sigma()<< endl;
            for (int x=0; x < Un.nx(); x++)
                for (int y=0; y < Un.ny(); y++) 
                {
                    Un(x,y,0) += R(x,y,0);
                    Un(x,y,1) += R(x,y,1);
                }
       }

            
		// Curl free decomposition
		R = Data - (Un+Vn);
		DC->transform(R, False);
		for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
		for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) 
        {
           (DC->RotDiv_Trans)(x,y,0) = 0;
            Trans (x,y,1) += (DC->RotDiv_Trans)(x,y,1);
        }

		// Reconstruction
		DC->recons(R, False);		
		Vn += R;
        
        // Constraint on mode B
        if (UseConstraint == True)
        {
           // cout << "B mode cst " << endl;
            EB = Vn;
            EB -= Un;
            if ((UseConstraint_on_E == True) && (UseConstraint_on_B == False))
            {
               for (int x=0; x < EB.nx(); x++)
               for (int y=0; y < EB.ny(); y++) EB(x,y,1)= 0.;
            }
            if ((UseConstraint_on_B == True) && (UseConstraint_on_E == False))
            {
                for (int x=0; x < EB.nx(); x++)
                    for (int y=0; y < EB.ny(); y++) EB(x,y,0)= 0.;
            }
            
            DC->transform(EB,False);
            for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
            for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) 
            { 
                if ((x > 1) && (y  > 1))  
                {
                   // (DC->RotDiv_Trans)(x,y,1) = 0;
                }
                if (Use_SoftThreshold == True)         (DC->RotDiv_Trans)(x,y,1) = (DC->RotDiv_Trans)(x,y,1) - soft_threshold( (DC->RotDiv_Trans)(x,y,1),  Lambda);
                else  if (Use_HardThreshold == True)   (DC->RotDiv_Trans)(x,y,1) = (DC->RotDiv_Trans)(x,y,1) - hard_threshold( (DC->RotDiv_Trans)(x,y,1),  Lambda);
                (DC->RotDiv_Trans)(x,y,0) = 0;
            }
            
            if (BordZeroConstraint == True)
            {
                float *tab = (DC->RotDiv_Trans).buffer() + (DC->RotDiv_Trans).nx() * (DC->RotDiv_Trans).ny();
                int N=(DC->RotDiv_Trans).nx()-2;
                set_to_0_border_wavelet(N,tab);
            }   
                            
            DC->recons(R,False);
            if (Verbose == True) cout << "   CURL BMODE: sig = " << R.sigma() << endl;
            for (int x=0; x < Un.nx(); x++)
            for (int y=0; y < Un.ny(); y++) 
            {
                Vn(x,y,0) -= R(x,y,0);
                Vn(x,y,1) -= R(x,y,1);
            }
            
       }
}
    
/*    
    if (BorderWavelet == True)
    {
        for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
        for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) 
        {
           if ((x <= 1) || (y <=1))  (DC->RotDiv_Trans)(x,y,0) = 0;
           else (DC->RotDiv_Trans)(x,y,0) = Trans (x,y,0);
           (DC->RotDiv_Trans)(x,y,1) = 0;
        }
        DC->recons(U);
        
        
        for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
        for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) 
        {
            if ((x <= 1) || (y <=1))  (DC->RotDiv_Trans)(x,y,1) = 0;
            else (DC->RotDiv_Trans)(x,y,1) = Trans (x,y,1);
            (DC->RotDiv_Trans)(x,y,0) = 0;
        }
        DC->recons(V, False);
    }
    else 
    {
        U = Un;
        V = Vn;
    }
    */
    U = Un;
    V = Vn;

}

/***************************************/

void Helmholtz::cat_iter_trans(fltarray &Data)
{
    char FN[256];
    int Nx = HCat->nx();
    int Ny = HCat->ny();
    int Nz = 2;
	fltarray R,Un;
	fltarray Vn;
	fltarray Resi;
    Bool BorderWavelet = DC->BorderWavelet;
    Bool BordZeroConstraint = False;
    Bool Use_SoftThreshold = False;
    Bool Use_HardThreshold = False;
    Bool UseConstraint;         
    
    if (THCons != C_NO) UseConstraint = True;
    else UseConstraint = False;
    if ((BorderWavelet == True) && (THCons == C_ZERO_BORDER)) BordZeroConstraint = True;
    if (THCons == C_HARD) Use_HardThreshold = True;
    if (THCons == C_SOFT) Use_SoftThreshold = True;
    // Bool UseConstraint_on_B = True;   // if true,  All common modes in E mode
    // Bool UseConstraint_on_E =  (UseConstraint_on_B == True) ? False: True; // if true,  All common modes in B mode
    
    
    
	// (DC->RotDiv_Trans).alloc(Data.nx(), Data.ny(), Data.nz());
    if (Verbose == True)
    {
        if (BorderWavelet == True) cout << "CAT ITER Wavelet Helmholtz with Border-wavelets " << endl;
        else cout << "CAT ITER Wavelet Helmholtz " << endl;
        if (UseConstraint_on_B == True)  cout << "      : Use constraint on B mode." << endl;
        if (UseConstraint_on_E == True)  cout << "      : Use constraint on E mode." << endl;
        if (BordZeroConstraint == True)  cout << "      : Use Border zero contraints." << endl;
        if (Use_SoftThreshold == True)  cout << "      : l1 constraint (soft thresholding)." << endl;
        if (Use_HardThreshold == True)  cout << "      : l0 constraint (hard thresholding)." << endl;
    }


    
	U.resize(Nx,Ny,Nz);
	V.resize(Nx,Ny,Nz);
    
    Un.resize(HCat->np(), 2);
	Vn.resize(HCat->np(), 2);
	R.resize(HCat->np(), 2);
    fltarray Ima2D;
    
    (DC->RotDiv_Trans).alloc(Nx,Ny,Nz);
         
    if (UseCat == True) DC->V0_Proj = False;
    if (UseCat == True) DC->PointValueRec = False;
    int Debug=1;
  	for (int i=0; i < Niter; i++)
	{
    // div_quasi_interpol(fltarray & Shear_or_Pol1D, fltarray & Ima2D);
    // div_point_value(fltarray & Ima2D, fltarray & Shear_or_Pol1D);
    
		if (UseCat == True) R = HCat->DataG - (Un+Vn);
        else R = Data - (Un+Vn);
        
		if (Verbose == True) cout << "Iter " << i+1 << ", Sigma Resi = " << R.sigma() << endl;
		
		// Divergence free decomposition
        if (Debug == True) cout << " Divergence free decomposition " << endl;
        if (UseCat == True) 
        {
           if (Debug == True) cout << " div_quasi_interpol " << endl;
           if (BorderWavelet == False) (HCat->div_quasi_interpol)(R,Ima2D);
           else (HCat->div_bord_quasi_interpol)(R, Ima2D, Nx);
           
            if (Debug == True) 
            {
               sprintf(FN, "xx_resi_divquasi_%d.fits", i+1);
               fits_write_fltarr(FN, Ima2D);		
            }
		    DC->transform(Ima2D);
        }
        else DC->transform(R);
        
		// Set to zero the second component
		for (int x=0; x < (DC->RotDiv_Trans).nx(); x++)
        for (int y=0; y < (DC->RotDiv_Trans).ny(); y++) (DC->RotDiv_Trans)(x,y,1) = 0;
		
        // Reconstruction
        if (UseCat == True) 
        {
           if (Debug == True) cout << " Divergence free Rec " << endl;
		   DC->recons(Ima2D);
           if (Debug == True) cout << " div_point_value " << endl;
           if (BorderWavelet == False) (HCat->div_point_value)(Ima2D, R);
           else (HCat->div_bord_point_value)(Ima2D, R, Nx);
        }
        else DC->recons(R);
        
        if (Debug == True) cout << " Un update " << R.naxis() << " " << R.nx() << " " << R.ny() << endl;
        Un += R;
        if (Debug == True) cout << " Un " << Un.naxis() << " " << Un.nx() << " " << Un.ny() << endl;
		
        if (Debug == True) cout << " Residual calculation " << endl;
 		if (UseCat == True) R = HCat->DataG - (Un+Vn);
        else R = Data - (Un+Vn);

        // Curl free decomposition
        if (Debug == True) cout << " Curl free decomposition " << endl;
        if (UseCat == True) 
        {
            if (BorderWavelet == False) (HCat->curl_quasi_interpol)(R, Ima2D);
            else  (HCat->curl_bord_quasi_interpol)(R, Ima2D, Nx);
            if (Debug == True) 
            {
                sprintf(FN, "xx_resi_curlquasi_%d.fits", i+1);
                fits_write_fltarr(FN, Ima2D);		
            }
            DC->transform(Ima2D, False);
        }
        else  DC->transform(R, False);
        
        // Set to zero the fitst component
		for (int x=0; x < Ima2D.nx(); x++)
        for (int y=0; y < Ima2D.ny(); y++) (DC->RotDiv_Trans)(x,y,0) = 0;
        
		// Reconstruction
        if (Debug == True) cout << " Curl free Rec " << endl;
		if (UseCat == True) 
        {
            DC->recons(Ima2D, False);
            if (BorderWavelet == False) (HCat->curl_point_value)(Ima2D, R);
            else (HCat->curl_bord_point_value)(Ima2D, R, Nx);
        }
        else DC->recons(R, False);		
		Vn += R;
        if (Debug == True) 
        {
            HCat->cat2ima(Un, U);
            HCat->cat2ima(Vn, V);
            sprintf(FN, "xx_un_%d.fits", i+1);
            fits_write_fltarr(FN, U);
            sprintf(FN, "xx_vn_%d.fits", i+1);
            fits_write_fltarr(FN, V);
            V -= U;
            sprintf(FN, "xx_eb_%d.fits", i+1);
            fits_write_fltarr(FN, V);
        }
        
        if (Debug == True) cout << " Vn " << Vn.naxis() << " " << Vn.nx() << " " << Vn.ny() << endl;
	}
    
    if (UseCat == True) 
    {
       HCat->cat2ima(Un, U);
       HCat->cat2ima(Vn, V);
       HCat->cat2ima(Data);
    }
    else
    {
	   U = Un;
	   V = Vn;
    }
}

/***************************************/

void Helmholtz::decomposition(fltarray &Data)
{
	if (Verbose == True) cout << StringHelmholtzTransform(THelm) << endl;

    switch (THelm)
	{
      case H_ANI_WT: 
            if (UseCat ==False)
            {
			   MAD->PointValueRec=PointValueRec;
			   MAD->V0_Proj=V0_Proj; 
			   iter_trans(Data);
            }
            else cat_iter_trans(Data);
 		   break;
        case H_ISO_WT: 
            if (UseCat ==False)
            {
 			    MID->PointValueRec=PointValueRec;
			    MID->V0_Proj=V0_Proj;           
			   iter_trans(Data);
            }
            else cat_iter_trans(Data);
			break;
		case H_FFT: 
	        DC = NULL;
            if (UseCat ==False) fft_trans(Data);
            else
            {
                HCat->cat2ima(Data);
                Data.info();
                fft_trans(Data);
            }
			break;
	   default:
	        cout << "Error: Undefined decomposition" << endl;
			exit(-1);
 	        break;
	}
}

/***************************************/
											  
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
 
	fprintf(OUTMAN, "         [-t TypeHelmholtz_Decompposition]\n");
    fprintf(OUTMAN, "             1: Anisotropic Wavelet Decomposition.\n"); 
    fprintf(OUTMAN, "             2: Isotropic Wavelet Decomposition.\n"); 
    fprintf(OUTMAN, "             3: Fourier Transform.\n"); 
	fprintf(OUTMAN, "             Default is 1.\n"); 

 	fprintf(OUTMAN, "         [-c TypeHelmholtz_Constraint]\n");
    fprintf(OUTMAN, "             1: No constraint.\n"); 
    fprintf(OUTMAN, "             2: l0 norm constraint.\n"); 
    fprintf(OUTMAN, "             3: l1 norm constraint.\n"); 
    fprintf(OUTMAN, "             4: zero border constraint.\n");
    fprintf(OUTMAN, "             5: component to zero.\n");
	fprintf(OUTMAN, "             Default is 1.\n"); 
    
    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Number of iterations in the wavelet Helmholtz wavelet decomposition.\n"); 
    fprintf(OUTMAN, "             Default is %d\n", NbrIter); 
	
    fprintf(OUTMAN, "         [-W EB_FileName]\n");
    fprintf(OUTMAN, "             Write the EB mode components to the disk.\n");
	
	fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Do not apply the V0 projection or the Point Value interpolation.\n");	

	fprintf(OUTMAN, "         [-C ImageSize]\n");
    fprintf(OUTMAN, "             Input data file is a catalog instead of an image. \n");	
    fprintf(OUTMAN, "             ImageSize is the final image size (must be a power of 2). \n");	

    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Use Border Wavelets.\n");	
    
    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "             Use constraints on E mode.\n");	
    
    fprintf(OUTMAN, "         [-B]\n");
    fprintf(OUTMAN, "             Use constraints on B mode.\n");
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv, (char *) "c:EBbC:Pi:W:t:vzZ")) != -1) 
    {
	switch (c) 
        {
            case 'c':
				if (sscanf(OptArg,"%d",&c ) != 1) 
				{
					fprintf(OUTMAN, "ErrorL bad type of helmholtz decomposition: %s\n", OptArg);
					exit(-1);
 				}
                if ((c > 0) && (c <= 5)) 
                    THCons = (type_constraints_helmholtz) (c-1);
                else  
                {
					fprintf(OUTMAN, "Error: bad type of helmholtz decomposition: %s\n", OptArg);
					exit(-1);
				}
 				break;
            case 'E': UseConstraint_on_E = True; break;
            case 'B': UseConstraint_on_B = True; break;
            case 'b': WaveletBorder = True; break;
			case 't':
				if (sscanf(OptArg,"%d",&c ) != 1) 
				{
					fprintf(OUTMAN, "ErrorL bad type of helmholtz decomposition: %s\n", OptArg);
					exit(-1);
 				}
                if ((c > 0) && (c <= 3)) 
					 THelm = (type_helmholtz) (c-1);
                else  
                {
					fprintf(OUTMAN, "Error: bad type of helmholtz decomposition: %s\n", OptArg);
					exit(-1);
				}
 				break;
			case 'i':
				if (sscanf(OptArg,"%d",&NbrIter) != 1) 
                {
					fprintf(OUTMAN, "Error: bad number of iterations: %s\n", OptArg);
					exit(-1);
                    
				}
 				break;
 			case 'W':  
				if (sscanf(OptArg,"%s", Name_EB_FileOut) != 1) 
                {
					fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
					exit(-1);
				}
				EBMod = True; 
				break;		
			case 'C': 
                if (sscanf(OptArg,"%d",&ImageSize) != 1) 
				{
					fprintf(OUTMAN, "ErrorL bad image size: %s\n", OptArg);
					exit(-1);
 				}
                ReadCat = (ReadCat == True) ? False: True; 
             break;
			case 'P': V0Proj = (V0Proj == True) ? False: True; break;
		    case 'v': Verbose = True;break;
		    case '?': usage(argv); break;
	  default: usage(argv); break;
 		}
	} 
 
       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out_2, argv[OptInd++]);
	else usage(argv);
	
	
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}

/*********************************************************************/


int main(int argc, char *argv[])
{
    int k;
    
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

    if (Verbose == True)
    {
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out U = " << Name_Imag_Out << endl; 
		cout << "File Name Out V = " << Name_Imag_Out_2 << endl; 
		if (EBMod == True) cout << " EB FIle Name = " << Name_EB_FileOut<< endl;
		cout << "Transform: " << StringHelmholtzTransform(THelm) << endl;
		cout << "Constraint: " << StringHelmholtzConstraint(THCons) << endl;
        
     }
    if (THCons == C_ZERO_BORDER) WaveletBorder = True;

    fltarray Data;
	intarray IData;

	fits_read_fltarr(Name_Imag_In, Data);
	int Nx = Data.nx();
    int Ny = Data.ny();
	
    Helmholtz HD;
    HD.Verbose = Verbose;
    HelmCatalog Cat;


	if (Data.naxis() == 2) ReadCat = True;
  	if (ReadCat == True) 
	{
		// TabXYG(i,0) = x
		// TabXYG(i,1) = y
		// TabXYG(i,2) = G1
		// TabXYG(i,3) = G2
		// TabXYG(i,4) = SigmaX
        // TabXYG(i,5) = SigmaX
		
		fltarray TabXYG;
		// fits_read_fltarr(Name_Imag_In, TabXYG);
		TabXYG = Data;
		if (TabXYG.naxis() != 2) 
		{
			cout << "Error: catalog file must be a 2D file ... " << endl;
			exit(-1);
		}
		if (TabXYG.ny() != 6)
		{
			cout << "Error: catalog second dimension must be 5 (x,y,G1,G2,sigma) ... " << endl;
			exit(-1);
		}
        
        Cat.allocn(TabXYG, ImageSize);
        if (Verbose == True) Cat.TabNp.info();

        fltarray Ima2D, Shear;
        Shear = Cat.DataG;
        
        Cat.DataG.info("G");
 
        Cat.div_bord_quasi_interpol(Cat.DataG,  Ima2D, Cat.nx());
        Cat.div_bord_point_value(Ima2D, Shear, Cat.nx() );
          
        Shear -= Cat.DataG;
        Shear.info("xx_err");
        fits_write_fltarr("xx_err.fits", Shear);
        Cat.DataG = Shear;
        Cat.cat2ima(Ima2D);
        fits_write_fltarr(Name_Imag_Out,Ima2D);
        fits_write_fltarr(Name_Imag_Out_2,Ima2D);

        exit(0);


        
        //cout << "Cat2ima " << endl;
        //Cat.cat2ima(Data);
        //fits_write_fltarr("xx_data.fits", Data);
        // exit(0);
        /*
        // div_quasi_interpol(Data,Tabc);

        Data.info();
        fltarray ImaSig;
        cout << "Cat2sigma " << endl;
        Cat.cat2imasigma(ImaSig);
        fits_write_fltarr("xx_data.fits", Data);		
		fits_write_fltarr("xx_sigma.fits", ImaSig);	
        */
                
        // fits_write_intarr("xx_np.fits", Cat.TabNp);
        HD.BorderWavelet = WaveletBorder;
        HD.alloc(THelm, Cat);
        HD.V0_Proj = True;
        
        
    }	
    else
    {	 
       // float Min,Max;
       // float *Ptr = Data.buffer();  // pointer to the data
       if (Verbose == True)
       { 
       cout << "Nx = " << Nx <<  " Ny = " << Ny  <<  endl; 
       Data.info("Input stat: ");
       // cout << "Min = " << Data.min()  << "  Max = " << Data.max() <<  " Sigma = " << Data.sigma() << endl;
       }
       HD.BorderWavelet = WaveletBorder;
       HD.alloc(THelm, Nx,Ny);
	   HD.V0_Proj = V0Proj;
    }
    
    HD.PointValueRec = V0Proj;
    HD.Niter = NbrIter;
    HD.UseConstraint_on_B = UseConstraint_on_B;
    HD.UseConstraint_on_E = UseConstraint_on_E;
    HD.THCons=  THCons;
    
	   // if (Verbose == True) 
    cout << "Transform .... " << endl;
    HD.decomposition(Data);

	if (Verbose == True) cout << "Write the result .... " << endl;
	fits_write_fltarr(Name_Imag_Out, HD.U);
	fits_write_fltarr(Name_Imag_Out_2, HD.V);
  	if (EBMod == True) 
    { 
		if (Verbose == True) cout << "Write EB result .... " << endl;
 		HD.EB = HD.V;
		HD.EB -= HD.U;
		fits_write_fltarr(Name_EB_FileOut, HD.EB);
		
		// Residual computation: Data - (u+v)
		HD.U += HD.V;
		Data -= HD.U;
		Data.info("data - (u+v)");
		fits_write_fltarr("xx_resi.fits", Data);		
	}
    exit(0);
} 

