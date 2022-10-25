/******************************************************************************
**                   Copyright (C) 2010 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  Sept 2010
**    
**    File:  PoissonMeyerWT.cc
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  3D POISSON MEYER WAVELET TRANSFORM
**    ----------- 
******************************************************************************/

#include "PoissonMeyerWT3D.h"
#define _STAB_PBCUR3D_

extern bool Verbose;

/*
N=128
.run estim_moments_pmwt
       1.0000000       1.0000000       1.0000000
      0.35355342     0.046665615    0.0037008122
      0.12500000    0.0058331145   0.00016355428
     0.044194174   0.00071994093   7.2281575e-06
     0.015625000   7.9333781e-05   3.1944364e-07
       1.0000000
       1.6817928
       2.8284271
       4.7568284
       7.9999999
      0.37500000
     0.075839044
     0.026812336
    0.0092341374
    0.0024294028
       1.0000000       1.0000000       1.0000000
      0.99999999     0.046665610    0.0011558362
      0.35355341    0.0058331152   5.1081363e-05
      0.12500000   0.00071994074   2.2544096e-06
     0.044194178   7.9333789e-05   9.5744281e-08
       1.0000000
       1.0000000
       1.6817928
       2.8284271
       4.7568282
      0.37500000
     0.028448169
     0.010057657
    0.0034738939
   0.00096730185
   
.run estim_norm_pmwt
Std for       30 cubes of size      128 with lambda=     100
      0.48657167      0.28061025      0.28049698      0.27600660       610.00104

N=128
.run estim_moments_ppmwt
       1.0000000       1.0000000       1.0000000
      0.35355336     0.052734383    0.0053720602
      0.12500000    0.0065917976   0.00023742326
     0.044194178   0.00082397466   1.0485930e-05
     0.015625000   0.00010299682   5.4636389e-07
       1.0000000
       1.6817929
       2.8284271
       4.7568283
       8.0000000
      0.37500000
     0.079575870
     0.028133591
    0.0099508470
    0.0031154884
*/
// Poisson transform with direct reconstruction, but no reconstruction filters
	// Decimated Scales
	// b = sgn(T1)/sqrt(T1)
	// c = 7T2/8T1 - T3/2T2
	static float VST_D_b[5] = 
		{ 1.00000, 1.68179, 2.82843, 4.75683, 8.00000};
	static float VST_D_c[5] = 
		{ 0.375000, 0.0758391, 0.0268123, 0.00923414, 0.00242940};

	// Undecimated Scales
	static float VST_U_b[5] = 
		{ 1.00000, 1.00000, 1.68179, 2.82843, 4.75683};
	static float VST_U_c[5] = 
		{ 0.375000, 0.0284482, 0.0100577, 0.00347389, 0.000967302};

	// Variance in the detail coefficients bands, for the four finest scales
	#define VST_norm_PMWT3D_maxscale 4
	static float VST_norm[4] = 
		{0.48657167, 0.28061025, 0.28049698, 0.27600660};
	//	{1, 1, 1, 1, 1};

// Poisson transform Coherent with the gaussian transform
	static float VST_C_b[5] = 
		{ 1.0000000, 1.6817929, 2.8284271, 4.7568283, 8.0000000};
	static float VST_C_c[5] = 
		{ 0.37500000, 0.079575870, 0.028133591, 0.0099508470, 0.0031154884};

	// Variance in the detail coefficients bands, for the four finest scales
	#define VST_norm_CPMWT3D_maxscale 4
	static float VST_C_norm[4] = 
		{ 0.48823129, 0.30375488, 0.30472992, 0.30403518};

/*********************************************************************/
static inline void PrintError( int status)
{
    // ***************************************************** 
    // * Print out cfitsio error messages and exit program * 
    // ***************************************************** 

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        // get the error status description 
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  // get first message; null if stack is empty 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  // get remaining messages 
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );
}

/****************************************************************************/

static inline void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/

static int MIN3(int a, int b, int c)
{
	return (a>b?MIN(b,c):MIN(a,c));
}


POISSON_MWT3D::POISSON_MWT3D()
{
	FFT3D.CenterZeroFreq = True;
	Isotrop=False;
	NbrScale=0;
	Tabcf_WT_Band=NULL;
	Type_PMW3D = DEF_TYPE_PMW3D;
}

POISSON_MWT3D::~POISSON_MWT3D()
{
	if (Tabcf_WT_Band !=NULL) delete [] Tabcf_WT_Band;
}

void POISSON_MWT3D::get_hfilter(fltarray &H, double DNx, double DNy, double DNz)
{
//	cerr<<"get_hfilter"<<endl;
	int i,j,k;
	int Nx = H.nx();
	int Ny = H.ny();
	int Nz = H.nz();
	double Nx2 = (DNx > 0) ? DNx /2.: Nx / 2.;
	double Ny2 = (DNy > 0) ? DNy /2.: Ny / 2.;
	double Nz2 = (DNz > 0) ? DNz /2.: Nz / 2.;
	double Nx4 = Nx2 / 2.;
	double Ny4 = Ny2 / 2.;
	double Nz4 = Nz2 / 2.;

//	cerr<<" Nx DNX Nx2 Nx4 Nx/4. : "<<Nx<<"; "<<DNx<<"; "<<Nx2<<"; "<<Nx4<<"; "<<Nx/4.<<endl;
	if (Isotrop == False)
	{
		dblarray H1(Nx);
		dblarray H2(Ny);
		dblarray H3(Nz);

		H1.init(1.);
		H2.init(1.);
		H3.init(1.);
		for(i=0; i < Nx; i++)
		{
			double x = (i - int(Nx/2));
			double r = (ABS(x) - Nx4) / Nx4;
			if (r <= 0) H1(i) = 1.;
			else if (r < 1) H1(i) = lowpass_window(1. - r);
			else H1(i) = 0.;
		}
		for(j=0; j< Ny; j++)
		{
			double x = (j - int(Ny/2));
			double r = (ABS(x) - Ny4) / Ny4;
			if (r <= 0) H2(j) = 1.;
			else if (r < 1) H2(j) = lowpass_window(1. - r);
			else H2(j) = 0.;
		}
		for(k=0; k< Nz; k++)
		{
			double x = (k - int(Nz/2));
			double r = (ABS(x) - Nz4) / Nz4;
			if (r <= 0) H3(k) = 1.;
			else if (r < 1) H3(k) = lowpass_window(1. - r);
			else H3(k) = 0.;
		}
//		cerr<<" Recopie filtre"<<endl;
		for(i=0; i< Nx; i++)
			for(j=0; j< Ny; j++) 
				for(k=0; k< Nz; k++) 
					H(i,j,k) = H1(i)*H2(j)*H3(k);
	}
	else
	{
		// cout << "ISOTROP" << endl;
		int N = MIN3(Nx/4,Ny/4,Nz/4);
		H.init(1.);
		for(i=0; i< Nx; i++)
			for(j=0; j< Ny; j++)
				for(k=0; k< Nz; k++) 
				{
					double x = (i - Nx/2);
					double y = (j - Ny/2);
					double z = (k - Nz/2);
					double r = (sqrt(x*x+y*y+z*z) - N) / N;
					if (r <= 0) H(i,j,k) = 1.;
					else if (r < 1) H(i,j,k) = lowpass_window(1. - r);
					else H(i,j,k) = 0.;
				}
	}
	// io_write_cube_float("low_pass.fits", H);
	// exit(0);
//	cerr<<"end get_hfilter"<<endl;
}

/*********************************************************************/


void POISSON_MWT3D::init(int Nbr_Scale, int Nx, int Ny, int Nz, Bool IsotropWT, Bool WTNeedOddSize)
{
	bool LocVerbose = false && Verbose;
	if(Verbose) cerr<<"POISSON_MWT3D::init("<<Nbr_Scale<<","<<Nx<<","<<Ny<<","<<Nz<<","<<IsotropWT<<","<<WTNeedOddSize<<")"<<endl;

	NbrScale = Nbr_Scale;
	Nx_Cube = Nx;
	Ny_Cube = Ny;
	Nz_Cube = Nz;
	TabNx.alloc(NbrScale);
	TabNy.alloc(NbrScale);
	TabNz.alloc(NbrScale);
	Isotrop = IsotropWT;
	NeedOddSize = WTNeedOddSize;

	if (LocVerbose)
	{
		cout << " INIT WT: " << "CubeSize = " << Nx << " " << Ny << " " << Nz << " NbrScale = " << NbrScale << endl;
		if (Isotrop == True) cout << "   Use isotropic wavelets " << endl;
		else cout << "   Use Meyer's wavelets with cube extention (4/3) " << endl;
	}

	D_ExtNx = (float) Nx;
	D_ExtNy = (float) Ny;
	D_ExtNz = (float) Nz;
	if (LocVerbose) cerr<<" Nx "<<Nx<<","<<Ny<<","<<Nz<<endl;
	if (LocVerbose) cerr<<" D_ExtNx "<<D_ExtNx<<","<<D_ExtNy<<","<<D_ExtNz<<endl;
	
	double DNX=D_ExtNx;
	double DNY=D_ExtNy;
	double DNZ=D_ExtNz;
	
	if(LocVerbose) cerr<<"Need odd size = "<<NeedOddSize<<endl;	
	for (int s=0; s < Nbr_Scale; s++)
	{
		if (NeedOddSize == True)
		{
			TabNx(s) = 2 * int(floor(DNX/2)) + 1; 
			TabNy(s) = 2 * int(floor(DNY/2)) + 1;  
			TabNz(s) = 2 * int(floor(DNZ/2)) + 1;  
		}
		else 
		{
			TabNx(s) = int( DNX + 0.5);
			TabNy(s) = int( DNY + 0.5);
			TabNz(s) = int( DNZ + 0.5);
		}
		DNX /= 2.;
		DNY /= 2.;
		DNZ /= 2.;
		if (LocVerbose) cerr<<" TabNx("<<s<<") "<<TabNx(s)<<","<<TabNy(s)<<","<<TabNz(s)<<endl;
	}
	ExtNx = TabNx(0);
	ExtNy = TabNy(0);
	ExtNz = TabNz(0);
	if (LocVerbose) cerr<<" ExtNx=Tab(0) "<<ExtNx<<","<<ExtNy<<","<<ExtNz<<endl;
	
// Memory allocation
	Tabcf_WT_Band = new cfarray[NbrScale];
	for  (int s=0; s < NbrScale; s++)
		Tabcf_WT_Band[s].alloc(TabNx(s), TabNy(s), TabNz(s));
	TF_ExtData.alloc(ExtNx,ExtNy,ExtNz);
	H.alloc(ExtNx,ExtNy,ExtNz);
	
	if(Verbose) cerr<<"...end POISSON_MWT3D::init"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::get_extFourier(fltarray &Data, cfarray & TF_ExtData)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"WT::get_extFourier..."<<endl;
	
	int i,j,k;
	int Nx = Data.nx();
	int Ny = Data.ny();
	int Nz = Data.nz();

	if(LocVerbose) cerr<<" Nxyz "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
	if(LocVerbose) cerr<<" ext xyz "<<ExtNx<<" "<<ExtNy<<" "<<ExtNz<<endl;
	TF_ExtData.resize(ExtNx,ExtNy,ExtNz);
	if(LocVerbose) cerr<<" extDataxyz "<<TF_ExtData.nx()<<" "<<TF_ExtData.ny()<<" "<<TF_ExtData.nz()<<endl;
	
	FFT3D.fftn3d(Data, TF_ExtData, False);
	float  Norm = sqrt(float(Nx*Ny*Nz));
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				TF_ExtData(i,j,k) /= Norm;

	if(Verbose) cerr<<"...End WT::get_extFourier"<<endl;
}

/*********************************************************************/

float POISSON_MWT3D::get_norm(int s)
{
	// coef that takes into account the Fourier zero padding (rounded sizes and sometimes imposed odd size)
	// coef = sqrt( theoretical_redundancy / real_redundancy )
	double coef = sqrt( (float(ExtNx)/pow((double) 2,(double) s)*float(ExtNy)/pow((double) 2,(double) s)*float(ExtNz)/pow((double) 2,(double) s))/(TabNx(s)*TabNy(s)*TabNz(s)) );
	if( (s<NbrScale-1) && (s<VST_norm_PMWT3D_maxscale) )
	{
		if(Type_PMW3D == PMW3D_SUM)
			return VST_norm[s]*coef;
		else if(Type_PMW3D == PMW3D_COHERENT)
			return VST_C_norm[s]*coef;
		else 
			return coef;
	}
	else
		return coef;
}

/*********************************************************************/

void POISSON_MWT3D::extract_stat(fltarray *TabBand, char* Outname, bool normalize)
{
//	bool LocVerbose = true & Verbose;
	if(Verbose) cerr<<"Extract_stat... normalize="<<normalize<<endl;
	
// Output stat files
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_stat.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_maxcoef.dat");
	fstream ccoef;
	ccoef.open (MaxCoefname, fstream::out);
	
	for(int s=0;s<NbrScale;s++)
	{
		double Mean, Sig, Skew, Curt;
		float Mini, Maxi, norm=1;
		if(normalize)
			norm=get_norm(s);
		
		moment4(TabBand[s].buffer(), TabBand[s].n_elem(), Mean, Sig, Skew, Curt, Mini, Maxi);
		
	// Normalisation
		Maxi/=norm; Mini/=norm;
		Mean/=norm; Sig/=norm; Skew/=(norm*norm*norm); Curt/=(norm*norm*norm*norm);
		
	// Output
		if(Verbose)
		{
			cerr<<" Stat scale("<<s<<",norm="<<norm<<") : " <<Mean<<","<<Sig<<","<<Skew<<","<<Curt<<","<<endl;
			cerr << "  MinMaxCoef Values = "<< Mini << ":"<< Maxi<<endl;
		}
		cstat << Mean <<"\t"<< Sig <<"\t"<< Skew<<"\t"<<Curt<<"\t"<<endl;
		ccoef << s  << "\t" << Mini << "\t" << Maxi << endl;
	}
	
	cstat.close();
	ccoef.close();
	
	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::normalize(fltarray *TabBand, fltarray *TabBandNorm)
{
	if(Verbose) cerr<<"POISSON_MWT3D::normalize..."<<endl;
	for(int s=0;s<NbrScale;s++)
	{
		float norm = get_norm(s);
		for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
				for (int k=0; k < TabNz(s); k++)
					(TabBandNorm[s])(i,j,k)=(TabBand[s])(i,j,k)/norm;
	}
	if(Verbose) cerr<<"...End POISSON_MWT3D::normalize"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::normalize_self(fltarray *TabBand, bool inverse)
{
	if(Verbose) cerr<<"POISSON_MWT3D::normalize_self(.,"<<inverse<<")"<<endl;
	
	// Fine scales
	for(int s=0;s<NbrScale;s++)
	{
		float norm = get_norm(s);
		for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
				for (int k=0; k < TabNz(s); k++)
					if(inverse)
						(TabBand[s])(i,j,k) *= norm;
					else
						(TabBand[s])(i,j,k) /= norm;
	}
	
	if(Verbose) cerr<<"...End POISSON_MWT3D::normalize_self"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{
	if(Verbose) cerr<<"POISSON_MWT3D::threshold(.,"<<","<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;

	float lvl;
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		float cnt=0;
		int nx = TabBand[s].nx();
		int ny = TabBand[s].ny();
		int nz = TabBand[s].nz();
		int tot = nx*ny*nz;
		
		float norm = get_norm(s);
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
				for(int k = 0; k < nz; k++)
				{
					float Nsig=NSigma;
					if(s==0) Nsig=(force4 ? (NSigma+1) : NSigma);
					lvl = SigmaNoise * Nsig * norm;

					if( abs((TabBand[s])(i,j,k)) < lvl )
					{
						cnt++;
						(TabBand[s])(i,j,k)=0; // hard
					}
					else if(FilterType==FT_SOFT) (TabBand[s])(i,j,k) -= (2*int( (TabBand[s])(i,j,k)>0 )-1)*lvl;
				}
		if(Verbose) cerr<<" Scale "<<s<<", n#proportion non seuillee ("<<lvl<<")="<<tot-cnt<<"#"<<(tot-cnt)/tot<<endl;
	}

	// Coarse scale
	// Do not denoise the coarsest scale

	if(Verbose) cerr<<"...End POISSON_MWT3D::threshold"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::wiener(fltarray *TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"POISSON_MWT3D::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
//	bool LocVerbose = false & Verbose;
	
	// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		float norm = get_norm(s);
		int nx = nxs(s);
		int ny = nys(s);
		int nz = nzs(s);
		int Nx = nx/LocalBS + 1 ;
		int Ny = ny/LocalBS + 1 ;
		int Nz = nz/LocalBS + 1 ;
		fltarray coef_wiener(Nx,Ny,Nz);

		coef_wiener.init(-2);

		// Wiener blocks
		for(int kx=0 ; kx < Nx ; kx++)
			for(int ky=0 ; ky < Ny ; ky++)
				for(int kz = 0; kz < Nz; kz++)
				{
					double sigma2 = 0.0;
					float cnt=0;

				// Sigma calculation
					// Pixels in a wiener block
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								cnt+=1;
								// spatial position = position in the block + block_position
								int x = (bx + kx*LocalBS) % nx;
								int y = (by + ky*LocalBS) % ny;
								int z = (bz + kz*LocalBS) % nz;
								sigma2+=pow((double) (TabBand[s])(x,y,z) / norm,(double) 2);
							}
					float sig = sqrt(max( 0.0, sigma2/cnt - pow((double)  noise_lvl,(double) 2) ));
					float norm = sig / (sig+noise_lvl);
					coef_wiener(kx,ky,kz) = norm;
				}


	// Apply the wiener coefficient
		// Wiener blocks (int angle2 and space) for mainly vertical Ridgelets
		for(int kx=0 ; kx < Nx ; kx++)
			for(int ky=0 ; ky < Ny ; ky++)
				for(int kz = 0; kz < Nz; kz++)
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								int x = (bx + kx*LocalBS);
								int y = (by + ky*LocalBS);
								int z = (bz + kz*LocalBS);
								if( x<nx && y<ny && z<nz )
									(TabBand[s])(x,y,z) *= coef_wiener(kx,ky,kz);
							}
	}// end scale

	// Coarse scale
	// Do not denoise the coarsest scale

	if(Verbose) cerr<<"...End POISSON_MWT3D::wiener"<<endl;
}

/*********************************************************************/
static float sign(float a)
{
	return float((a>0)-(a<0));
}

static void stabilize(cfarray & A, int s, bool decimated)
{
#ifdef _STAB_PBCUR3D_
//	cerr<<"stab scale "<<s<<", "<<"decimated="<<decimated<<endl;
	float b,c;
	if(decimated)
	{
		b = VST_D_b[s];
		c = VST_D_c[s];
	}
	else
	{
		b = VST_U_b[s];
		c = VST_U_c[s];
	}
	
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = complex_f(b * sign(A(i,j,k).real()+c) * sqrt(abs(A(i,j,k).real()+c)),0);
//				A(i,j,k) = complex_f( A(i,j,k).real()*(s+1)+(s+1) ,0);
			}
#endif
}
static void unstabilize(cfarray & A, int s, bool decimated)
{
#ifdef _STAB_PBCUR3D_
//	cerr<<"unstab scale "<<s<<", "<<"decimated="<<decimated<<endl;
	float b,c;
	if(decimated)
	{
		b = VST_D_b[s];
		c = VST_D_c[s];
	}
	else
	{
		b = VST_U_b[s];
		c = VST_U_c[s];
	}
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = complex_f( sign(A(i,j,k).real())*A(i,j,k).real()*A(i,j,k).real()/(b*b) - c ,0);
//				A(i,j,k) = complex_f( A(i,j,k).real()/(s+1)+(s+1) ,0);
			}
#endif
}

static void stabilize(cfarray & A, int s)
{
#ifdef _STAB_PBCUR3D_
//	cerr<<"stab scale "<<s<<", "<<"decimated="<<decimated<<endl;
	float b,c;
	b = VST_C_b[s];
	c = VST_C_c[s];
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = complex_f(b * sign(A(i,j,k).real()+c) * sqrt(abs(A(i,j,k).real()+c)),0);
//				A(i,j,k) = complex_f( A(i,j,k).real()*(s+1)+(s+1) ,0);
			}
#endif
}
static void unstabilize(cfarray & A, int s)
{
#ifdef _STAB_PBCUR3D_
//	cerr<<"unstab scale "<<s<<", "<<"decimated="<<decimated<<endl;
	float b,c;
	b = VST_C_b[s];
	c = VST_C_c[s];
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = complex_f( sign(A(i,j,k).real())*A(i,j,k).real()*A(i,j,k).real()/(b*b) - c ,0);
//				A(i,j,k) = complex_f( A(i,j,k).real()/(s+1)+(s+1) ,0);
			}
#endif
}

void POISSON_MWT3D::transform(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc)
{
	if(Type_PMW3D == PMW3D_SUM)
		transform_sum(Data, Tab_Cube, Alloc);
	else if(Type_PMW3D == PMW3D_COHERENT)
		transform_coherent(Data, Tab_Cube, Alloc);
	else
	{
		cerr<<"Wrong Poisson Meyer Wavelet Type"<<endl;
		exit(0);
	}
}

void POISSON_MWT3D::transform_sum(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc)
{
	bool LocVerbose = true;//false & Verbose;
	if(Verbose) cerr<<"POISSON_MWT3D::transformifft..."<<endl;
	
	int i,j,k,s,Nxs,Nys,Nzs;
	int Nx = TF_ExtData.nx();
	int Ny = TF_ExtData.ny();
	int Nz = TF_ExtData.ny();
	double DNX = D_ExtNx;
	double DNY = D_ExtNy;
	double DNZ = D_ExtNz;
	float  Norm;
	
	if (Alloc == True) 
	{
		Tab_Cube = new fltarray[NbrScale];
		for  (int s=0; s < NbrScale-1; s++)
			Tab_Cube[s].alloc(TabNx(s), TabNy(s), TabNz(s));
		Tab_Cube[NbrScale-1].alloc(TabNx(NbrScale-1), TabNy(NbrScale-1), TabNz(NbrScale-1));
	}
	
	if(Verbose) cerr<<"POISSON_MWT3D::transform(.)..."<<endl;
	cfarray Tilde(Nx,Ny,Nz);
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	for(k=0; k<Nz; k++)
		Tilde(i,j,k) = complex_f(Data(i,j,k),0);
	cfarray TF_D;
	get_extFourier(Data, TF_ExtData); //FT
	
//	if(LocVerbose) 
//	{
//		char filename[64];
//		sprintf(filename,"%s_F.fits","out");
//		writecfarr(filename, TF_ExtData);
//	}

	// input : Tilde and its Fourier transform TF_ExtData
	for (s=0; s < NbrScale-1; s++) 
	{
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_Tilde%d.fits","out",s);
//		writecfarr(filename, Tilde);
//		}
		Nx = TabNx(s);// size of the current scale
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);// size of next scale
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		H.resize(Nxs,Nys,Nzs); // lowpass filter for next scale
		get_hfilter(H); // low pass filter for next scale

		TF_D.resize(Nxs,Nys,Nzs); // allocate to size of next scale
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
//			if(Verbose) cerr<<" s "<<s<<" ijk Indijk"<<i<<" "<<j<<" "<<k<<"   "<<Indi<<" "<<Indj<<" "<<Indk<<endl;
			
		// decimated
			TF_D(i,j,k) = H(i,j,k)*H(i,j,k) * TF_ExtData(Indi, Indj, Indk); // DWFT
		}
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_TFD%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}
		TF_ExtData.init(complex_f(0,0));
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
		// undecimated
			TF_ExtData(Indi, Indj, Indk) = TF_D(i,j,k); // WFT
		}
		
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_UTFD%d.fits","out",s);
		writecfarr(filename, TF_ExtData);
		}
		
	// Inverse fourier on the undecimated low scale
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) *= Norm;
		FFT3D.ifftn3d(TF_ExtData);// FWFT = a1
		
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_Uscale%d.fits","out",s+1);
		writecfarr(filename, TF_ExtData);
		}
		
	// w1 = tilde0 - a1    =>  TF_ExtData = tilde - TF_ExtData
		stabilize(Tilde,s,true);
		stabilize(TF_ExtData,s+1,false);
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_TUscale%d.fits","out",s);
		writecfarr(filename, TF_ExtData);
		}
		
		TF_ExtData = Tilde - TF_ExtData;
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_coef%d.fits","out",s);
//		writecfarr(filename, TF_ExtData);
//		}
		
	// export w1.real
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					Tab_Cube[s](i,j,k) = TF_ExtData(i,j,k).real();
		
	// initialize next FT
		TF_ExtData = TF_D;
		
	// Inverse fourier on the Decimated low scale = next Tilde
		Norm = sqrt(float(Nxs*Nys*Nzs));
		for(i=0; i<Nxs; i++)
			for(j=0; j<Nys; j++)
				for(k=0; k<Nzs; k++)
					TF_D(i,j,k) *= Norm;
		FFT3D.ifftn3d(TF_D);
		Tilde = TF_D;// Tilde = FDWFT
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_Dscale%d.fits","out",s+1);
		writecfarr(filename, Tilde);
		}
		
		DNX =  DNX/2.;
		DNY =  DNY/2.;
		DNZ =  DNZ/2.;
	}
	
// export last scale decimated : Tilde
	for(i=0; i<Nxs; i++)
		for(j=0; j<Nys; j++)
			for(k=0; k<Nzs; k++)
				Tab_Cube[NbrScale-1](i,j,k) = TF_D(i,j,k).real();
//	if(LocVerbose) 
//	{
//	char filename[64];
//	sprintf(filename,"%s_coef%d.fits","out",s);
//	writefltarr(filename, Tab_Cube[NbrScale-1]);
//	}
	
	if(Verbose) cerr<<"...End POISSON_MWT3D::transform"<<endl;
	
	if(Verbose) cerr<<"End POISSON_MWT3D::transformifft"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::transform_coherent(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc)
{
	bool LocVerbose = true;//false & Verbose;
	if(Verbose) cerr<<"POISSON_MWT3DBIS::transformifft..."<<endl;
	
	int i,j,k,s,Nxs,Nys,Nzs;
	int Nx = TF_ExtData.nx();
	int Ny = TF_ExtData.ny();
	int Nz = TF_ExtData.ny();
	double DNX = D_ExtNx;
	double DNY = D_ExtNy;
	double DNZ = D_ExtNz;
	float  Norm;
	
	if (Alloc == True) 
	{
		Tab_Cube = new fltarray[NbrScale];
		for  (int s=0; s < NbrScale; s++)
			Tab_Cube[s].alloc(TabNx(s), TabNy(s), TabNz(s));
	}
	
	if(Verbose) cerr<<"POISSON_MWT3DBIS::transform(.)..."<<endl;
	cfarray Tilde(Nx,Ny,Nz);
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	for(k=0; k<Nz; k++)
		Tilde(i,j,k) = complex_f(Data(i,j,k),0);
	cfarray TF_D;
	get_extFourier(Data, TF_ExtData); //FT
	
	if(LocVerbose) 
	{
		char filename[64];
		sprintf(filename,"%s_F.fits","out");
		writecfarr(filename, TF_ExtData);
	}

	// input : Tilde and its Fourier transform TF_ExtData
	for (s=0; s < NbrScale-1; s++) 
	{
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_Tilde%d.fits","out",s);
//		writecfarr(filename, Tilde);
//		}
		Nx = TabNx(s);// size of the current scale
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);// size of next scale
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		H.resize(Nxs,Nys,Nzs); // lowpass filter for next scale
		get_hfilter(H); // low pass filter for next scale

		TF_D.resize(Nxs,Nys,Nzs); // allocate to size of next scale
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
//			if(Verbose) cerr<<" s "<<s<<" ijk Indijk"<<i<<" "<<j<<" "<<k<<"   "<<Indi<<" "<<Indj<<" "<<Indk<<endl;
			
		// decimated
			TF_D(i,j,k) = H(i,j,k) * TF_ExtData(Indi, Indj, Indk); // DWFT
		}
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_TFD%d.fits","out",s);
		writecfarr(filename, TF_D);
		}

	// Stabilize the scale
		stabilize(Tilde,s);
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_stab%d.fits","out",s);
		writecfarr(filename, Tilde);
		}
		
	// Fourier of the stabilized scale
		TF_ExtData=Tilde;
		FFT3D.fftn3d(TF_ExtData);
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) /= Norm;
		
	// Apply the wavelet high-pass filter 
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
			TF_ExtData(Indi, Indj, Indk) *= sqrt(1.-H(i,j,k)*H(i,j,k));
		}

		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_WF%d.fits","out",s);
		writecfarr(filename, TF_ExtData);
		}

	// Inverse fourier to get the detail coefficients
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) *= Norm;
		FFT3D.ifftn3d(TF_ExtData);

		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_coef%d.fits","out",s);
		writecfarr(filename, TF_ExtData);
		}

	// Export w1.real
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					Tab_Cube[s](i,j,k) = TF_ExtData(i,j,k).real();

	// initialize next FT
		TF_ExtData = TF_D;
		
	// Inverse fourier on the Decimated low scale = next Tilde
		Tilde.alloc(Nxs,Nys,Nzs);
		Norm = sqrt(float(Nxs*Nys*Nzs));
		for(i=0; i<Nxs; i++)
			for(j=0; j<Nys; j++)
				for(k=0; k<Nzs; k++)
					Tilde(i,j,k) = TF_D(i,j,k) * Norm;
		FFT3D.ifftn3d(Tilde);
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_scale%d.fits","out",s+1);
		writecfarr(filename, Tilde);
		}
		
		DNX =  DNX/2.;
		DNY =  DNY/2.;
		DNZ =  DNZ/2.;
	}
	
// export last scale decimated : Tilde
	for(i=0; i<Nxs; i++)
		for(j=0; j<Nys; j++)
			for(k=0; k<Nzs; k++)
				Tab_Cube[NbrScale-1](i,j,k) = Tilde(i,j,k).real();
	
	if(LocVerbose) 
	{
	char filename[64];
	sprintf(filename,"%s_last%d.fits","out",s);
	writefltarr(filename, Tab_Cube[NbrScale-1]);
	}
	
	if(Verbose) cerr<<"...End POISSON_MWT3DBIS::transform"<<endl;
	
	if(Verbose) cerr<<"End POISSON_MWT3DBIS::transformifft"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{
	if(Type_PMW3D == PMW3D_SUM)
		recons_sum(Tab_Cube, Data, Alloc);
	else if(Type_PMW3D == PMW3D_COHERENT)
		recons_coherent(Tab_Cube, Data, Alloc);
	else
	{
		cerr<<"Wrong Poisson Meyer Wavelet Type"<<endl;
		exit(0);
	}
}

void POISSON_MWT3D::recons_sum(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"POISSON_MWT3D::recons(.)..."<<endl;
	
// Init
	int i,j,k,s,Nx,Ny,Nz,Nxs,Nys,Nzs;
	float  Norm;
	if (Alloc == True) 
		Data.resize(Nx_Cube, Ny_Cube, Nz_Cube);
	
// Previous approximation
	cfarray TF_D(TabNx(NbrScale-1),TabNy(NbrScale-1),TabNz(NbrScale-1));
	for(i=0; i<TabNx(NbrScale-1); i++)
		for(j=0; j<TabNy(NbrScale-1); j++)
			for(k=0; k<TabNz(NbrScale-1); k++)
				TF_D(i,j,k) = Tab_Cube[NbrScale-1](i,j,k);
	
	if(LocVerbose) 
	{
		char filename[64];
		sprintf(filename,"%s_F.fits","out");
		writecfarr(filename, TF_D);
	}

// input : Previous decimated approximation : TF_D
//			and fine scales in direct space : Tab_Cube[s]
	for (s = NbrScale-2; s > -1; s--)
	{
//		cerr<<" Scale "<<s<<endl;
		Nx = TabNx(s);// size of the current scale
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);// size of the approximation
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		
	// Fourier transform of the decimated scale
		FFT3D.fftn3d(TF_D);
		Norm = sqrt(float(Nxs*Nys*Nzs));
		for(i=0; i<Nxs; i++)
			for(j=0; j<Nys; j++)
				for(k=0; k<Nzs; k++)
					TF_D(i,j,k) /= Norm;
		
	// Enlarge the approximation in Fourier
		TF_ExtData.resize(Nx,Ny,Nz);
		TF_ExtData.init(complex_f(0,0));
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
//			if(Verbose) cerr<<" s "<<s<<" ijk Indijk"<<i<<" "<<j<<" "<<k<<"   "<<Indi<<" "<<Indj<<" "<<Indk<<endl;
			
		// undecimated approximation
			TF_ExtData(Indi, Indj, Indk) = TF_D(i,j,k);
		}

	// stabilize the approximation in direct space
	//  T(F-1(TF_ExtData))
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) *= Norm;
		FFT3D.ifftn3d(TF_ExtData);
		stabilize(TF_ExtData,s+1,false);
		
	
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_u%d.fits","out",s);
//		writecfarr(filename, TF_ExtData);
//		sprintf(filename,"%s_d%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}
		
	// Current scale
		TF_D.alloc(Nx,Ny,Nz);
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_D(i,j,k) = Tab_Cube[s](i,j,k);
		
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_wf%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}

	// Add previous scale
		for (i=0; i < Nx; i++) // loop on big scale
		for (j=0; j < Ny; j++)
		for (k=0; k < Nz; k++)
			TF_D(i,j,k) += TF_ExtData(i,j,k);
		
		unstabilize(TF_D,s,true);
		
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_af%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}

	}
	
//	if(LocVerbose) 
//	{
//	char filename[64];
//	sprintf(filename,"%s_a.fits","out");
//	writecfarr(filename, TF_D);
//	}

// Export the solution
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				Data(i,j,k) = TF_D(i,j,k).real();

	if(Verbose) cerr<<"...End POISSON_MWT3D::recons"<<endl;
	
	if(Verbose) cerr<<"End POISSON_MWT3D::reconsifft"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::recons_coherent(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{
	bool LocVerbose =  false & Verbose;
	if(Verbose) cerr<<"POISSON_MWT3DBIS::recons(.)..."<<endl;
	
// Init
	int i,j,k,s,Nx,Ny,Nz,Nxs,Nys,Nzs;
	s=NbrScale-1;
	Nx = TabNx(s);// size of the approximation
	Ny = TabNy(s);
	Nz = TabNz(s);
	float  Norm;
	if (Alloc == True) 
		Data.resize(Nx_Cube, Ny_Cube, Nz_Cube);
	
// Previous approximation
	cfarray TF_D(Nx,Ny,Nz);
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				TF_D(i,j,k) = Tab_Cube[s](i,j,k);
	
// Fourier transform of the previous approx
	FFT3D.fftn3d(TF_D);
	Norm = sqrt(float(Nx*Ny*Nz));
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				TF_D(i,j,k) /= Norm;
	
	if(LocVerbose) 
	{
		char filename[64];
		sprintf(filename,"%s_F.fits","out");
		writecfarr(filename, TF_D);
	}

// input : Previous decimated approximation, in Fourier: TF_ExtData
//			and fine scales in direct space : Tab_Cube[s]
	for (s = NbrScale-2; s > -1; s--)
	{
//		cerr<<" Scale "<<s<<endl;
		Nx = TabNx(s);// size of the current scale
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);// size of the approximation
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		H.resize(Nxs,Nys,Nzs); // lowpass filter for next scale
		get_hfilter(H); // low pass filter for next scale
		
	// Enlarge the approximation in Fourier
		TF_ExtData.resize(Nx,Ny,Nz);
		TF_ExtData.init(complex_f(0,0));
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;

		// undecimated approximation
			TF_ExtData(Indi, Indj, Indk) = TF_D(i,j,k);
		}

	// Current scale
		TF_D.alloc(Nx,Ny,Nz);
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_D(i,j,k) = Tab_Cube[s](i,j,k);
		
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_wf%d.fits","out",s);
		writecfarr(filename, TF_D);
		}

		unstabilize(TF_D,s);
		
	// Fourier transform of the scale
		FFT3D.fftn3d(TF_D);
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_D(i,j,k) /= Norm;
		
	// Mix with previous scale, and apply the filters
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
			TF_D(Indi,Indj,Indk) = sqrt(1-H(i,j,k)*H(i,j,k))*TF_D(Indi,Indj,Indk) + H(i,j,k)*TF_ExtData(Indi,Indj,Indk);
		}
		
		if(LocVerbose) 
		{
		char filename[64];
		sprintf(filename,"%s_af%d.fits","out",s);
		writecfarr(filename, TF_D);
		}

	}
	
// Inverse fourier of the solution
	Norm = sqrt(float(Nx*Ny*Nz));
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				TF_D(i,j,k) *= Norm;
	FFT3D.ifftn3d(TF_D);

	if(LocVerbose) 
	{
	char filename[64];
	sprintf(filename,"%s_a.fits","out");
	writecfarr(filename, TF_D);
	}

// Export the solution
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				Data(i,j,k) = TF_D(i,j,k).real();

	if(Verbose) cerr<<"...End POISSON_MWT3DBIS::recons"<<endl;
	
	if(Verbose) cerr<<"End POISSON_MWT3DBIS::reconsifft"<<endl;
}

/*********************************************************************/

/*
void POISSON_MWT3D::recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{// optimized for the gaussian case, but wrong in Poisson mode
	bool LocVerbose = true;// & Verbose;
	if(Verbose) cerr<<"POISSON_MWT3D::recons(.)..."<<endl;
	
// Init
	int i,j,k,s,Nx,Ny,Nz,Nxs,Nys,Nzs;
	float  Norm;
	if (Alloc == True) 
		Data.resize(Nx_Cube, Ny_Cube, Nz_Cube);
	
// Fourier transform of the approximation
	cfarray TF_D(TabNx(NbrScale-1),TabNy(NbrScale-1),TabNz(NbrScale-1));
	FFT3D.fftn3d(Tab_Cube[NbrScale-1], TF_D, False);
//	cerr<<"tabNx = "<<TabNx(NbrScale-1)<<", .nx()="<<TF_D.nx()<<endl;
	Norm = sqrt(float(TabNx(NbrScale-1)*TabNy(NbrScale-1)*TabNz(NbrScale-1)));
	for(i=0; i<TabNx(NbrScale-1); i++)
		for(j=0; j<TabNy(NbrScale-1); j++)
			for(k=0; k<TabNz(NbrScale-1); k++)
				TF_D(i,j,k) /= Norm;
	
//	if(LocVerbose) 
//	{
//		char filename[64];
//		sprintf(filename,"%s_F.fits","out");
//		writecfarr(filename, TF_D);
//	}

// input : Fourier of previous decimated approximation : TF_D
//			and fine scales in direct space : Tab_Cube[s]
	for (s = NbrScale-2; s > -1; s--)
	{
//		cerr<<" Scale "<<s<<endl;
		Nx = TabNx(s);// size of the current scale
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);// size of the approximation
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		
	// Enlarge the approximation in Fourier
		TF_ExtData.resize(Nx,Ny,Nz);
		TF_ExtData.init(complex_f(0,0));
		for (i=0; i < Nxs; i++) // loop on small scale
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2; // coordinates in big scale
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
//			if(Verbose) cerr<<" s "<<s<<" ijk Indijk"<<i<<" "<<j<<" "<<k<<"   "<<Indi<<" "<<Indj<<" "<<Indk<<endl;
			
		// undecimated approximation
			TF_ExtData(Indi, Indj, Indk) = TF_D(i,j,k);
		}
		if(s==NbrScale-2) // stabilize the approximation in direct space
		{ //  F(T(F-1(TF_ExtData)))
		// Inverse fourier of the approximation
			Norm = sqrt(float(Nx*Ny*Nz));
			for(i=0; i<Nx; i++)
				for(j=0; j<Ny; j++)
					for(k=0; k<Nz; k++)
						TF_ExtData(i,j,k) *= Norm;
			FFT3D.ifftn3d(TF_ExtData);
			stabilize(TF_ExtData,NbrScale-1);
			FFT3D.fftn3d(TF_ExtData);
			float  Norm = sqrt(float(Nx*Ny*Nz));
			for(i=0; i<Nx; i++)
				for(j=0; j<Ny; j++)
					for(k=0; k<Nz; k++)
						TF_ExtData(i,j,k) /= Norm;
		}
	
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_u%d.fits","out",s);
//		writecfarr(filename, TF_ExtData);
//		sprintf(filename,"%s_d%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}
		
	// Fourier transform of the current scale
		TF_D.alloc(Nx,Ny,Nz);
		FFT3D.fftn3d(Tab_Cube[s], TF_D, False);
//		cerr<<"tabNx = "<<TabNx(s)<<", .nx()="<<Tab_Cube[s].nx()<<", .nx()="<<TF_D.nx()<<endl;
		Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_D(i,j,k) /= Norm;
		
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_wf%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}

	// Add the current scale
//		cerr<<"Nx="<<Nx<<", approx="<<TF_ExtData.nx()<<", wf="<<TF_D.nx()<<endl;
		for (i=0; i < Nx; i++) // loop on big scale
		for (j=0; j < Ny; j++)
		for (k=0; k < Nz; k++)
			TF_D(i,j,k) += TF_ExtData(i,j,k);
		
//		if(LocVerbose) 
//		{
//		char filename[64];
//		sprintf(filename,"%s_af%d.fits","out",s);
//		writecfarr(filename, TF_D);
//		}

	}
	
// Inverse fourier of the solution
	Norm = sqrt(float(Nx*Ny*Nz));
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				TF_D(i,j,k) *= Norm;
	FFT3D.ifftn3d(TF_D);
	
	unstabilize(TF_D,0);
//	if(LocVerbose) 
//	{
//	char filename[64];
//	sprintf(filename,"%s_a.fits","out");
//	writecfarr(filename, TF_D);
//	}

// Export the solution
	for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
			for(k=0; k<Nz; k++)
				Data(i,j,k) = TF_D(i,j,k).real();

	if(Verbose) cerr<<"...End POISSON_MWT3D::recons"<<endl;
	
	if(Verbose) cerr<<"End POISSON_MWT3D::reconsifft"<<endl;
}
*/

/*********************************************************************/

void POISSON_MWT3D::test(fltarray *TabBand)
{
/*	
	for(int s3=0;s3<NbrScale;s3++)
	{
		cerr<<"s nx"<<s3<<" "<<TabBand[s3].nx()<<endl;
		for(int i=0;i<TabBand[s3].nx();i++)
		for(int j=0;j<TabBand[s3].ny();j++)
		for(int k=0;k<TabBand[s3].nz();k++)
		TabBand[s3](i,j,k)=0;
	}
	
	TabBand[2](0,0,0)=sqrt(Nx_Cube*Ny_Cube*Nz_Cube);
*/	
}

/*********************************************************************/

void POISSON_MWT3D::write (char *Name, fltarray * & data, bool Normalize)
{
	if(Verbose) cerr<<"POISSON_MWT3D::write_multi..."<<endl;
	
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1; 

	// .mr extention
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 
	
// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;

// write first header part (parameters)
	naxis=0;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"Meyer 3D Wavele Transform", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) Nx_Cube, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) Ny_Cube, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) Nz_Cube, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
	for (int s=0; s < nbr_scale(); s++)
	{
		naxis=3;
		naxes[0] = nxs(s);
		naxes[1] = nys(s);
		naxes[2] = nzs(s);

		if(ffcrhd(fptr,&status))
			PrintError( status );
		if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
			PrintError( status );
		
	// save the data
		if ( ffppre(fptr, group, 1, nxs(s)*nys(s)*nzs(s), (data[s]).buffer(), &status) )
			PrintError( status ); 
	}

// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end POISSON_MWT3D::write_multi"<<endl;
}

/*********************************************************************/

void POISSON_MWT3D::read(char *Name, fltarray * & data, bool *Normalize)
{
	if(Verbose) cerr<<"POISSON_MWT3D::read_multi..."<<endl;
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	char comment[FLEN_COMMENT];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );

// get number of Headers
	int nhead;
	fits_get_num_hdus(fptr, &nhead, &status);
	
// read primary header
	if ( ffmahd(fptr, 1, &hdutype, &status) ) PrintError( status );
	
	// Read params
	int _NbrScale,_Nx,_Ny,_Nz;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*Normalize = (bool)mon_long;
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
	init(_NbrScale, _Nx, _Ny, _Nz );
	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<NbrScale<<endl; exit(0); }
	
//	init the structure and the output vector
	data = new fltarray[NbrScale];
	for(int s=0; s < NbrScale; s++)
		data[s].alloc(TabNx(s), TabNy(s), TabNz(s));
	
// read data
	for(int s=0;s<_NbrScale;s++)
	{
		if (fits_movabs_hdu(fptr, s+2, NULL, &status)) PrintError( status );
		if (ffgpve(fptr, 1, 1, nxs(s)*nys(s)*nzs(s), nulval, (data[s]).buffer(), &anynul, &status)) PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end POISSON_MWT3D::read_multi"<<endl;
}



/*********************************************************************/

/*
void POISSON_MWT3D::write_mono (char *Name, fltarray * & data)
{
	if(Verbose) cerr<<"POISSON_MWT3D::write..."<<endl;
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1;
	int Nelem=0;

// .mr 
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 

// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;
	
	naxis = 1;
	Nelem=0;
	for (int s=0; s < nbr_scale(); s++)  Nelem += 3 + nxs(s)*nys(s)*nzs(s);
	naxes[0] = Nelem;

// write first header part (parameters)
	if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
		PrintError( status );          // call PrintError if error occurs 

// write optional keyword to the header 
	if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"Meyer 3D Wavele Transform", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of bands 3D", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) Nx_Cube, (char*)"x size of the original cube", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) Ny_Cube, (char*)"y size of the original cube", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) Nz_Cube, (char*)"z size of the original cube", &status))
		PrintError( status );  
	
// save the data
	int pos=1;
	for (int s=0; s < nbr_scale(); s++)
	{
		float* Nxyz = new float[3];
		Nxyz[0] = (float) nxs(s);
		Nxyz[1] = (float) nxs(s);
		Nxyz[2] = (float) nxs(s);
		if ( ffppre(fptr, group, pos, 3, Nxyz, &status) )
			PrintError( status );  
		pos += 3;
		if ( ffppre(fptr, group, pos, nxs(s)*nys(s)*nzs(s), (data[s]).buffer(), &status) )
			PrintError( status );  
		pos += nxs(s)*nys(s)*nzs(s);
	}

	// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end POISSON_MWT3D::write"<<endl;
}
void POISSON_MWT3D::read_mono (char *Name, fltarray * & data)
{
	if(Verbose) cerr<<"POISSON_MWT3D::read..."<<endl;
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	long hdunum;
	char comment[FLEN_COMMENT];
	int naxis;
	long naxes[3];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );
	
// read header
	hdunum = 1;  //read  table 
	if ( ffmahd(fptr, hdunum, &hdutype, &status) ) // move to the HDU 
		PrintError( status );
 
	int simple, bitpix, extend;
	long pcount, gcount;
	if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount, &gcount, &extend, &status)) // move to the HDU 
		PrintError( status );
	
	int _NbrScale,_Nx,_Ny,_Nz;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) 
		PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status))
		PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status))
		PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status))
		PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status))
		PrintError( status );
	_Nz = (int)mon_long;
	
	
//	init the structure and the output vector
	init(_NbrScale, _Nx, _Ny, _Nz );
	data = new fltarray[NbrScale];
	for  (int s=0; s < NbrScale; s++)
		data[s].alloc(TabNx(s), TabNy(s), TabNz(s));
	
// read the images
	int pos=1;
	for(int s=0;s<_NbrScale;s++)
	{
		float* Nxyz = new float[3];
		if (ffgpve(fptr, 1, pos,3, nulval, Nxyz, &anynul, &status))
			PrintError( status );
		pos+=3;
		if (ffgpve(fptr, 1, pos, (int)(Nxyz[0]*Nxyz[1]*Nxyz[2]), nulval, (data[s]).buffer(), &anynul, &status))
			PrintError( status );
		pos += (int)(Nxyz[0]*Nxyz[1]*Nxyz[2]);
	}

	// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end POISSON_MWT3D::read"<<endl;
}
*/

