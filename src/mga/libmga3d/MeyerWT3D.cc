/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  20/04/2008
**    
**    File:  MeyerWT.cc
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  3D MEYER WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/
 

#include "MeyerWT3D.h"
extern bool Verbose;

// Not extended Meyer wavelets :
//  Theoretical coefficients, for perfect decimation
//  first scale : sqrt(485/2)/16 = 0.973276
//  intermediate scales : sqrt(189/2)/16 = 0.607569
//  last scale : (3/4)^(3/2) = 0.649519
static float Tab_Meyer[3] = 
	{	0.973276, 0.607569, 0.649519 };

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


MEYER_WT3D::MEYER_WT3D()
{
	FFT3D.CenterZeroFreq = True;
	Extend=False;
	Isotrop=False;
	NbrScale=0;
	Tabcf_WT_Band=NULL;
}

MEYER_WT3D::~MEYER_WT3D()
{
	if (Tabcf_WT_Band !=NULL) delete [] Tabcf_WT_Band;
}

void MEYER_WT3D::get_hfilter(fltarray &H, double DNx, double DNy, double DNz)
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
		if (Extend == True)
		{
//			cerr<<" Extend"<<endl;
			for(i=0; i< Nx4; i++)
			{
				double x = (i -  int(Nx/2) + Nx2) / Nx4;
				H1(i) =  H1(Nx-i-1) =  lowpass_window(x);
			}
			for(j=0; j< Ny4; j++)
			{
				double x = (j - int(Ny/2) + Ny2) / Ny4;
				H2(j) =  H2(Ny-j-1) =  lowpass_window(x);
			}
			for(k=0; k< Nz4; k++)
			{
				double x = (k - int(Nz/2) + Nz2) / Nz4;
				H3(k) =  H3(Nz-k-1) =  lowpass_window(x);
			}
		}
		else
		{
//			cerr<<" NoExtend"<<endl;
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


void MEYER_WT3D::init(int Nbr_Scale, int Nx, int Ny, int Nz, Bool ExtendWT, Bool IsotropWT, Bool WTNeedOddSize)
{
	bool LocVerbose = true && Verbose;
	if(LocVerbose) cerr<<"MEYER_WT3D::init("<<Nbr_Scale<<","<<Nx<<","<<Ny<<","<<Nz<<","<<ExtendWT<<","<<IsotropWT<<","<<WTNeedOddSize<<")"<<endl;

	NbrScale = Nbr_Scale;
	Nx_Cube = Nx;
	Ny_Cube = Ny;
	Nz_Cube = Nz;
	TabNx.alloc(NbrScale);
	TabNy.alloc(NbrScale);
	TabNz.alloc(NbrScale);
	Extend = ExtendWT;
	Isotrop = IsotropWT;
	NeedOddSize = WTNeedOddSize;

	if (Isotrop)
		Extend = False;
	if (LocVerbose)
	{
		cout << " INIT WT: " << "CubeSize = " << Nx << " " << Ny << " " << Nz << " NbrScale = " << NbrScale << endl;
		if (Isotrop == True) cout << "   Use isotropic wavelets " << endl;
		else if (Extend == False) cout << "   Use Meyer's wavelets without cube extension  " << endl;
		else cout << "   Use Meyer's wavelets with cube extention (4/3) " << endl;
	}

	if (Extend)
	{
		D_ExtNx =(4. / 3. * (float) Nx);
		D_ExtNy =(4. / 3. * (float) Ny);
		D_ExtNz =(4. / 3. * (float) Nz);
	}
	else
	{
		D_ExtNx = (float) Nx;
		D_ExtNy = (float) Ny;
		D_ExtNz = (float) Nz;
	}
	if (LocVerbose) cerr<<" Nx "<<Nx<<","<<Ny<<","<<Nz<<endl;
	if (LocVerbose) cerr<<" D_ExtNx "<<D_ExtNx<<","<<D_ExtNy<<","<<D_ExtNz<<endl;
	
	double DNX=D_ExtNx;
	double DNY=D_ExtNy;
	double DNZ=D_ExtNz;
	
	if(LocVerbose) cerr<<"Need odd size = "<<NeedOddSize<<endl;	
	for (int s=0; s < Nbr_Scale; s++)
	{
		if ((Extend==True) || (NeedOddSize == True))
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
	
	if(Verbose) cerr<<"...end MEYER_WT3D::init"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::get_extFourier(fltarray &Data, cfarray & TF_ExtData)
{
	bool LocVerbose = true & Verbose;
	if(Verbose) cerr<<"WT::get_extFourier..."<<endl;
	
	int i,j,k;
	int Nx = Data.nx();
	int Ny = Data.ny();
	int Nz = Data.nz();
	int Nx2 = Nx/2; 
	int Ny2 = Ny/2;  
	int Nz2 = Nz/2;  

	if(LocVerbose) cerr<<" Nxyz "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
	if(LocVerbose) cerr<<" ext xyz "<<ExtNx<<" "<<ExtNy<<" "<<ExtNz<<endl;
	TF_ExtData.resize(ExtNx,ExtNy,ExtNz);
	if(LocVerbose) cerr<<" extDataxyz "<<TF_ExtData.nx()<<" "<<TF_ExtData.ny()<<" "<<TF_ExtData.nz()<<endl;
	
	if (Extend == True)
	{
		cfarray TF_Cube(Nx,Ny,Nz,(char*)"TFIMA");
		int ExtNx2 = ExtNx/2;
		int ExtNy2 = ExtNy/2;
		int ExtNz2 = ExtNz/2;
		if(LocVerbose) cerr<<" ExtNx2 "<<ExtNx2<<" "<<ExtNy2<<" "<<ExtNz2<<endl;
		if(LocVerbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
		if(LocVerbose) cerr<<" Nx2 "<<Nx2<<" "<<Ny2<<" "<<Nz2<<endl;
		FFT3D.fftn3d(Data, TF_Cube, False);
		
		float  Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_Cube(i,j,k) /= Norm;
		for(i=0; i<ExtNx; i++)
			for(j=0; j<ExtNy; j++) 
				for(k=0; k<ExtNz; k++) 
				{
//					if(LocVerbose) cerr<<" ijk "<<i<<" "<<j<<" "<<k<<endl;
					int u = i - ExtNx2 + Nx2;
					int v = j - ExtNy2 + Ny2;
					int w = k - ExtNz2 + Nz2;
					if (u < 0) u += Nx;
					else if (u >= Nx) u -= Nx;
					if (v < 0) v += Ny;
					else if (v >= Ny) v -= Ny;
					if (w < 0) w += Nz;
					else if (w >= Nz) w -= Nz;
					TF_ExtData(i,j,k) = TF_Cube(u,v,w);
				}
	}
	else 
	{
		FFT3D.fftn3d(Data, TF_ExtData, False);
		float  Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) /= Norm;
	}    
	if(Verbose) cerr<<"...End WT::get_extFourier"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::get_IFFT_Cube(cfarray & TF_ExtData, fltarray &Data)
{
	if(Verbose) cerr<<"Get_IFFT_Cube..."<<endl;
	int i,j,k;
	int Nx = Data.nx();
	int Ny = Data.ny(); 
	int Nz = Data.nz(); 

	// cout << "get_extFourier" << endl;

	if ( (ExtNx != TF_ExtData.nx()) || (ExtNy != TF_ExtData.ny()) || (ExtNz != TF_ExtData.nz()))
	{
		cout << "Error: bad cube size in get_IFFT_Cube ... " << endl;
		cout << "ExtNx = " << ExtNx << "ExtNy = " << ExtNy << endl;
		cout << "InExtNx = " << TF_ExtData.nx() << "InExtNy = " << TF_ExtData.ny() << "InExtNz = " << TF_ExtData.nz() << endl;
		exit(-1);
	}

	// cout << "alloc" << ExtNx << " " << ExtNy << endl;
	if (Extend == False)
	{
		float  Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					TF_ExtData(i,j,k) *= Norm;
		FFT3D.ifftn3d(TF_ExtData);
		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
				for(k=0; k<Nz; k++)
					Data(i,j,k) = TF_ExtData(i,j,k).real();    
	}
	else
	{
		cfarray TF_Cube(Nx,Ny,Nz,(char*)"TFCUBE"); 
		int ExtNx2 = ExtNx/2;
		int ExtNy2 = ExtNy/2;
		int ExtNz2 = ExtNz/2;
		int Nx2 = Nx/2; 
		int Ny2 = Ny/2;  
		int Nz2 = Nz/2;  
		for(i=0; i<ExtNx; i++)
			for(j=0; j<ExtNy; j++) 
				for(k=0; k<ExtNz; k++) 
				{
					int u = i - ExtNx2;
					int v = j - ExtNy2;
					int w = k - ExtNz2;
					if (u < -Nx2) u += Nx;
					else if (u >= Nx2) u -= Nx;
					if (v < -Ny2) v += Ny;
					else if (v >= Ny2) v -= Ny;
					if (w < -Nz2) w += Nz;
					else if (w >= Nz2) w -= Nz;
					TF_Cube(u+Nx2,v+Ny2,w+Nz2) += TF_ExtData(i,j,k);
				}
		// cout << "END get_extFourier" << endl;
		float  Norm = sqrt(float(Nx*Ny*Nz));
		for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
		for(k=0; k<Nz; k++)  TF_Cube(i,j,k) *= Norm;

		FFT3D.ifftn3d(TF_Cube);
		for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++)
		for(k=0; k<Nz; k++)  Data(i,j,k) = TF_Cube(i,j,k).real();
	}
	if(Verbose) cerr<<"...End Get_IFFT_Cube"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::transform_cf(cfarray & TF_ExtData,  cfarray * &TabWT)
{
	if(Verbose) cerr<<"MEYER_WT3D::transform_cf..."<<endl;
	int i,j,k,s,Nxs,Nys,Nzs;
	int Nx = TF_ExtData.nx();
	int Ny = TF_ExtData.ny();
	int Nz = TF_ExtData.ny();
	double DNX = D_ExtNx;
	double DNY = D_ExtNy;
	double DNZ = D_ExtNz;

	TabWT[0] = TF_ExtData;
	if (Extend == True)
	{
		get_hfilter(H, DNX, DNY, DNZ);
		for (i=0; i < Nx; i++)
		for (j=0; j < Ny; j++)
		for (k=0; k < Nz; k++) (TabWT[0])(i,j,k) *= H(i,j,k);
	}
	for (s=0; s < NbrScale-1; s++) 
	{
		Nx = TabNx(s);
		Ny = TabNy(s);   
		Nz = TabNz(s);   
		Nxs = TabNx(s+1);
		Nys = TabNy(s+1);
		Nzs = TabNz(s+1);
//		if(Verbose) cerr<<" Nx "<<Nx<<" "<<Ny<<" "<<Nz<<endl;
//		if(Verbose) cerr<<" Nxs "<<Nxs<<" "<<Nys<<" "<<Nzs<<endl;
		H.resize(Nxs,Nys,Nzs);
		TF_ExtData.resize(Nxs,Nys,Nzs);
		if (Extend == True) get_hfilter(H,DNX/2.,DNY/2.,DNZ/2.);
		else get_hfilter(H,DNX/2.,DNY/2.,DNZ/2.);
//		else get_hfilter(H); // change A.Woiselle sept 2010
		
//cerr<<"writing H filter of size "<<H.nx()<<endl;
//char filename[64];
//sprintf(filename,"%s_s%d.fits","out_filter",s);
//writefltarr(filename, H);
		
		for (i=0; i < Nxs; i++)
		for (j=0; j < Nys; j++)
		for (k=0; k < Nzs; k++)
		{
			int Indi = i - Nxs/2 + Nx/2;
			int Indj = j - Nys/2 + Ny/2;
			int Indk = k - Nzs/2 + Nz/2;
//			if(Verbose) cerr<<" s "<<s<<" ijk Indijk"<<i<<" "<<j<<" "<<k<<"   "<<Indi<<" "<<Indj<<" "<<Indk<<endl;
			TF_ExtData(i,j,k) = H(i,j,k) * (TabWT[s])(Indi, Indj, Indk);
			(TabWT[s])(Indi, Indj, Indk) *= sqrt(1-H(i,j,k)*H(i,j,k));
		}
		TabWT[s+1].resize(Nx,Ny,Nz);
		TabWT[s+1] = TF_ExtData; 
		DNX =  DNX/2.;
		DNY =  DNY/2.;
		DNZ =  DNZ/2.;
	}
	if(Verbose) cerr<<"...End MEYER_WT3D::transform_cf"<<endl;
}
/*********************************************************************/

void MEYER_WT3D::recons_cf(cfarray * &TabWT,  cfarray & TF_ExtData)
{
	bool LocVerbose = False && Verbose;
	if(Verbose) cerr<<"MEYER_WT3D::recons_cf..."<<endl;

	int i,j,k,s;
	int Nx = TabWT[0].nx();
	int Ny = TabWT[0].ny();
	int Nz = TabWT[0].nz();
	int Nxd2 = Nx/2;
	int Nyd2 = Ny/2;
	int Nzd2 = Nz/2;
	fltarray H(Nx,Ny,Nz,(char*)"LowPass");
	fltarray HB(Nx,Ny,Nz,(char*)"LowPass");

	// for  (s=0; s < NbrScale; s++) FFT2D.fftn3d(TabWT[s]);

	// cout << "REC ITER" << Nx << " " << Ny <<  endl;
	TF_ExtData.alloc(TabWT[0].nx(), TabWT[0].ny(), TabWT[0].nz(), (char*)"rec");
	// cout << "DATA NEW SIZE = " << TF_ExtData.nx() << "  " << TF_ExtData.ny() << "  " << TF_ExtData.nz() << endl;
	double DNX=D_ExtNx;
	double DNY=D_ExtNy;
	double DNZ=D_ExtNz;

	if (Extend == True)
	{
		get_hfilter(HB,D_ExtNx,D_ExtNy,D_ExtNz);
		for (s=0; s < NbrScale-1; s++)
		{
			if (LocVerbose == True) cout << " Rec WT Scale " << s+1 << " " << TabNx(s) << " " << TabNy(s) << " " << TabNz(s) << endl;
			Nx = TabNx(s);
			Ny = TabNy(s);
			Nz = TabNz(s);
			int Nxs = TabNx(s+1);
			int Nys = TabNy(s+1);
			int Nzs = TabNz(s+1);
			H.resize(Nxs, Nys,Nzs);
			get_hfilter(H,DNX/2.,DNY/2.,DNZ/2.);
			for (i=0; i < Nx; i++)
			for (j=0; j < Ny; j++)
			for (k=0; k < Nz; k++)   (TabWT[s])(i,j,k) *= HB(i,j,k);
			for (i=0; i < Nxs; i++)
			for (j=0; j < Nys; j++)
			for (k=0; k < Nzs; k++)  (TabWT[s])(i - Nxs/2 + Nx/2, j - Nys/2 + Ny/2, k - Nzs/2 + Nz/2) *= sqrt(1-H(i,j,k)*H(i,j,k));
			for (i=0; i < Nx; i++)
			for (j=0; j < Ny; j++)
			for (k=0; k < Nz; k++)   TF_ExtData(i - Nx/2 + Nxd2, j - Ny/2 + Nyd2, k - Nz/2 + Nzd2) += (TabWT[s])(i, j, k);

			// char Name[256];
			// sprintf(Name, "mrband_%d" , s+1);
			// io_write_cube_complex_f(Name, TabWT[s]);
			if (s != NbrScale-2) HB = H;
			DNX =  DNX/2.;
			DNY =  DNY/2.;
			DNZ =  DNZ/2.;
		}
		s = NbrScale-1;
		Nx = TabWT[s].nx();
		Ny = TabWT[s].ny();
		Nz = TabWT[s].nz();
		for (i=0; i < Nx; i++)
		for (j=0; j < Ny; j++)
		for (k=0; k < Nz; k++)  TF_ExtData(i - Nx/2 + Nxd2, j - Ny/2 + Nyd2, k - Nz/2 + Nzd2) += (TabWT[s])(i, j, k)*H(i,j,k);
	}
	else
	{
		s = NbrScale-1;
		Nx = TabNx(s);
		Ny = TabNy(s);
		Nz = TabNz(s);
		TF_ExtData.resize(Nx, Ny, Nz);
		TF_ExtData = TabWT[s];
		for (s=NbrScale-2; s >= 0; s--)
		{
			if (LocVerbose == True) cout << "Rec WT Scale " << s+1 << " " << TabNx(s) << " " << TabNy(s) << " " << TabNz(s) << endl;
			Nx = TabNx(s);
			Ny = TabNy(s);
			Nz = TabNz(s);
			int Nxs = TabNx(s+1);
			int Nys = TabNy(s+1);
			int Nzs = TabNz(s+1);
			DNX =  D_ExtNx/pow((double) 2.,(double) s);
			DNY =  D_ExtNy/pow((double) 2.,(double) s);
			DNZ =  D_ExtNz/pow((double) 2.,(double) s);
			H.resize(Nxs, Nys, Nzs);
//cerr<<"scale "<<s<<" : DNX/2.="<<DNX/2.<<" and H.nx="<<H.nx()<<endl;
			get_hfilter(H,DNX/2.,DNY/2.,DNZ/2.);
//			get_hfilter(H); // change A.Woiselle sept 2010

			for (i=0; i < Nxs; i++)
			for (j=0; j < Nys; j++)
			for (k=0; k < Nzs; k++)
			{
				int Indi = i - Nxs/2 + Nx/2;
				int Indj = j - Nys/2 + Ny/2;
				int Indk = k - Nzs/2 + Nz/2;
				(TabWT[s])(Indi, Indj, Indk) *= sqrt(1-H(i,j,k)*H(i,j,k));
				TF_ExtData(i,j,k) *= H(i,j,k);
				(TabWT[s])(Indi, Indj, Indk) += TF_ExtData(i,j,k);
			}   
			TF_ExtData.resize(Nx,Ny,Nz);
			TF_ExtData = TabWT[s];
		}
	}
	if(Verbose) cerr<<"...End MEYER_WT3D::recons_cf"<<endl;
}
 
/*********************************************************************/

void MEYER_WT3D::ifft_tabcube(cfarray * & TabCF_Cube, fltarray * & Tab_Cube, Bool Alloc)
{
	if(Verbose) cerr<<"MEYER_WT3D::ifft_tabcube, alloc="<<(bool)Alloc<<endl;
    if (Alloc == True) Tab_Cube = new fltarray[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm =  sqrt((float)(TabNx(s) * TabNy(s) * TabNz(s)));
       if (Alloc == True) Tab_Cube[s].alloc(TabNx(s), TabNy(s), TabNz(s));
       for (int i=0; i < TabNx(s); i++)
       for (int j=0; j < TabNy(s); j++)
       for (int k=0; k < TabNz(s); k++) (TabCF_Cube[s])(i,j,k) *= Norm;
       FFT3D.ifftn3d(TabCF_Cube[s]);
       for (int i=0; i < TabNx(s); i++)
       for (int j=0; j < TabNy(s); j++)
       for (int k=0; k < TabNz(s); k++) (Tab_Cube[s])(i,j,k) = (TabCF_Cube[s])(i,j,k).real();
       
//        {
//          char Name[256];
//          sprintf(Name, "iband_%d" , s+1);
//          cout << "Scale " << s+1 << " " << Name << endl;
//          io_write_cube_complex_f(Name, TabCF_Cube[s]);
//        }    
    }
	if(Verbose) cerr<<"...End MEYER_WT3D::ifft_tabcube"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::fft_tabcube(fltarray * & Tab_Cube, cfarray * & TabCF_Cube)
{
    // if (Alloc == True)  TabCF_Cube = new cfarray[NbrScale];
    for  (int s=0; s < NbrScale; s++)
    {
       float Norm = sqrt( (float) TabNx(s) * TabNy(s) * TabNz(s));
       // if (Alloc == True) TabCF_Cube[s].alloc(TabNx(s), TabNy(s), TabNz(s));
       FFT3D.fftn3d(Tab_Cube[s], TabCF_Cube[s]);
       for (int i=0; i < TabNx(s); i++)
       for (int j=0; j < TabNy(s); j++)
       for (int k=0; k < TabNz(s); k++) (TabCF_Cube[s])(i,j,k) /= Norm;
    }
}

/*********************************************************************/

float MEYER_WT3D::get_norm(int s)
{
	// coef that takes into account the Fourier zero padding (rounded sizes and sometimes imposed odd size)
	// coef = sqrt( theoretical_redundancy / real_redundancy )
	double coef = sqrt( (float(D_ExtNx)/pow((double) 2,(double) s)*float(D_ExtNy)/pow((double) 2,(double) s)*float(D_ExtNz)/pow((double) 2,(double) s))/(TabNx(s)*TabNy(s)*TabNz(s)) );
	if(s==0 && Extend==False) return Tab_Meyer[0]*coef;
	else if(s==NbrScale-1) return Tab_Meyer[2]*coef;
	else return Tab_Meyer[1]*coef;
}

/*********************************************************************/

// The same as extract_stat except that there is no normalisation 
void MEYER_WT3D::noise_calibration(fltarray *TabBand, char* Outname)
{
	bool LocVerbose = true & Verbose;
	if(Verbose) cerr<<"Noise_calibration..."<<endl;
	
// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_noise_calib.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	
	double mean, Sig;
	
	for(int s=0;s<NbrScale;s++)
	{
		mean=(TabBand[s]).mean();
		Sig=(TabBand[s]).sigma();
		if(LocVerbose)
			cerr<< "\tStat echelle "<<s<<" : mean,sigma = "<<mean<<" "<<Sig<<endl;
		cstat <<mean<<"\t"<<Sig<<"\t"<<endl;
	}
	
// coarse scale : N-1
	// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats

	cstat.close();
	if(Verbose) cerr<<"...End Noise_calibration"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::extract_stat(fltarray *TabBand, char* Outname, bool normalize)
{
//	bool LocVerbose = true & Verbose;
	if(Verbose) cerr<<"Extract_stat..."<<endl;
	
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
		if(normalize) norm=get_norm(s);
		
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

void MEYER_WT3D::normalize(fltarray *TabBand, fltarray *TabBandNorm)
{
	if(Verbose) cerr<<"MEYER_WT3D::normalize..."<<endl;
	for(int s=0;s<NbrScale;s++)
	{
		float norm = get_norm(s);
		for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
				for (int k=0; k < TabNz(s); k++)
					(TabBandNorm[s])(i,j,k)=(TabBand[s])(i,j,k)/norm;
	}
	if(Verbose) cerr<<"...End MEYER_WT3D::normalize"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::normalize_self(fltarray *TabBand, bool inverse)
{
	if(Verbose) cerr<<"MEYER_WT3D::normalize_self(.,"<<inverse<<")"<<endl;
	
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
	
	if(Verbose) cerr<<"...End MEYER_WT3D::normalize_self"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{
	if(Verbose) cerr<<"MEYER_WT3D::threshold(.,"<<","<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;

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

	if(Verbose) cerr<<"...End MEYER_WT3D::threshold"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::wiener(fltarray *TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"MEYER_WT3D::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
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
								sigma2+=pow((double)  (TabBand[s])(x,y,z) / norm, (double) 2);
							}
					float sig = sqrt(max( 0.0, sigma2/cnt - pow((double) noise_lvl,(double) 2) ));
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

	if(Verbose) cerr<<"...End MEYER_WT3D::wiener"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::values_at(fltarray *TabBand, char * filename, char* Outname)
{
	int s,x,y,z;
	float m;
	
	// Output stream
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_values_at.dat");
	fstream cval;
	cval.open (MaxCoefname, fstream::out);
	
	FILE* fic;
	fic=fopen(filename, "r");

	if(fic)
	{
		for(int s3=0;s3<NbrScale-1;s3++)
		{
			if(!fscanf(fic,"%d %f %d %d %d",&s,&m,&x,&y,&z)) cerr<<"Error while reading coordinates in "<<filename<<endl;
			cval<<s<<"\t"<<abs((TabBand[s])(x,y,z))<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
		}
		fclose(fic);
	}
	else cerr<<"Warnig: File "<<filename<<" not found, for use in BCurvelet3D::values_at"<<endl;
	cval.close();
}

/*********************************************************************/

void MEYER_WT3D::transform(fltarray &Data)
{
	if(Verbose) cerr<<"MEYER_WT3D::transform(.)..."<<endl;
	bool LocVerbose = false & Verbose;
	get_extFourier(Data, TF_ExtData);
	
	if(LocVerbose) 
	{
		char filename[64];
		fltarray TabWavelet;
		cfarray TabWaveletF;
		TabWaveletF=TF_ExtData;

		cerr<<"size TF = "<<TF_ExtData.nx()<<endl;
		TabWavelet.alloc(nxs(0), nys(0), nzs(0));

		// save cube in fourier space
		for (int i=0; i < nxs(0); i++)
		for (int j=0; j < nys(0); j++)
		for (int k=0; k < nzs(0); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real();
		sprintf(filename,"%s_F.fits","out");
		writefltarr(filename, TabWavelet);
	}

	transform_cf(TF_ExtData, Tabcf_WT_Band);
	
	if(LocVerbose) 
	{
		char filename[64];
		for  (int s=0; s < NbrScale; s++)
		{
			fltarray TabWavelet;
			cfarray TabWaveletF;
			TabWaveletF=Tabcf_WT_Band[s];

//			float Norm =  sqrt((float)(TabNx(s) * TabNy(s) * TabNz(s)));
			TabWavelet.alloc(TabNx(s), TabNy(s), TabNz(s));

			// save module in fourier space
			for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
			for (int k=0; k < TabNz(s); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real();
			sprintf(filename,"%s_FW%d.fits","out",s);
			writefltarr(filename, TabWavelet);
		}
	}
	if(Verbose) cerr<<"...End MEYER_WT3D::transform"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::transform(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc)
{
	if(Verbose) cerr<<"MEYER_WT3D::transformifft..."<<endl;
	transform(Data);
	ifft_tabcube(Tabcf_WT_Band, Tab_Cube, Alloc);
	if(Verbose) cerr<<"End MEYER_WT3D::transformifft"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc)
{
   if (Alloc == True) Data.resize(Nx_Cube, Ny_Cube, Nz_Cube);
   fft_tabcube(Tab_Cube, Tabcf_WT_Band);
   recons(Data);
}

/*********************************************************************/

void MEYER_WT3D::recons(fltarray &Data)
{
	recons_cf(Tabcf_WT_Band, TF_ExtData);

	if(Verbose) 
	{
		char filename[64];
		fltarray TabWavelet;
		cfarray TabWaveletF;
		TabWaveletF=TF_ExtData;

		cerr<<"size rTF = "<<TF_ExtData.nx()<<endl;
		TabWavelet.alloc(nxs(0), nys(0), nzs(0));

		// save module in fourier space
		for (int i=0; i < nxs(0); i++)
		for (int j=0; j < nys(0); j++)
		for (int k=0; k < nzs(0); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real();
		sprintf(filename,"%s_rF.fits","out");
		writefltarr(filename, TabWavelet);
	}

	get_IFFT_Cube(TF_ExtData, Data);
}

/*********************************************************************/

void MEYER_WT3D::test(fltarray *TabBand)
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

void MEYER_WT3D::write (char *Name, fltarray * & data, bool Normalize)
{
	if(Verbose) cerr<<"MEYER_WT3D::write_multi..."<<endl;
	
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
	
	if(Verbose) cerr<<"...end MEYER_WT3D::write_multi"<<endl;
}

/*********************************************************************/

void MEYER_WT3D::read(char *Name, fltarray * & data, bool *Normalize)
{
	if(Verbose) cerr<<"MEYER_WT3D::read_multi..."<<endl;
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
	
	if(Verbose) cerr<<"...end MEYER_WT3D::read_multi"<<endl;
}



/*********************************************************************/

/*
void MEYER_WT3D::write_mono (char *Name, fltarray * & data)
{
	if(Verbose) cerr<<"MEYER_WT3D::write..."<<endl;
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
	
	if(Verbose) cerr<<"...end MEYER_WT3D::write"<<endl;
}
void MEYER_WT3D::read_mono (char *Name, fltarray * & data)
{
	if(Verbose) cerr<<"MEYER_WT3D::read..."<<endl;
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
	
	if(Verbose) cerr<<"...end MEYER_WT3D::read"<<endl;
}
*/

