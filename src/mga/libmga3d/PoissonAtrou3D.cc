/*******************************************************************************
**
**    DESCRIPTION  Sub-Band decomposition with A trous algorithm
**    -----------  
**                 
******************************************************************************/
 

#include "PoissonAtrou3D.h"
#define _STAB_PATROU3D_

  
/****************************************************************************/
//                   A TROUS ALGORITHM 3D
/****************************************************************************/

extern Bool Verbose;

static float Tab_Atrou3D_Moment[6][3] = 
	{{ 1           , 1           , 1            },
	 { 1.000000e+00, 2.044439e-02, 6.027655e-04 },
	 { 1.000000e+00, 1.882471e-03, 5.212584e-06 },
	 { 1.000000e+00, 2.199011e-04, 7.140385e-08 },
	 { 1.000000e+00, 2.703688e-05, 1.080409e-09 },
	 { 1.000000e+00, 3.361847e-06, 1.673232e-11 }};
// b = sgn(T1)/sqrt(T1)  // print, 1/sqrt(tau(0,*))
static float VST_b[6] = 
	{ 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000};
// c = 7T2/8T1 - T3/2T2  // print, 7*tau(1,*)/(8*tau(0,*)) - tau(2,*)/(2*tau(1,*))
static float VST_c[6] = 
	{ 3.750000e-01, 3.147254e-03, 2.626563e-04, 3.005900e-05, 3.676988e-06, 4.530555e-07};

static float VST_norm[6]=
		{0.476018, 0.0601403, 0.0176514, 0.00607454, 0.00219340, 0.0675106};

POISSONATROUS_3D_WT::POISSONATROUS_3D_WT ()
{
	Bord=I_MIRROR; 
	no_fine=false;
}

void POISSONATROUS_3D_WT::b3spline_filtering(fltarray & Old, fltarray & New, int s)
{
    int Nx = Old.nx();
    int Ny = Old.ny();
    int Nz = Old.nz();    
    int i,j,k,Step;
    double Coeff_h0 = 3. / 8.;
    double Coeff_h1 = 1. / 4.;
    double  Coeff_h2 = 1. / 16.;
    fltarray Buff1(Nx,Ny,Nz); 
    fltarray Buff2(Nx,Ny,Nz);
     
    Step = (int)(pow((double)2., (double) s) + 0.5);
	//convolution in Z axe
    for (i = 0; i < Nx; i ++)
    for (j = 0; j < Ny; j ++)
    for (k = 0; k < Nz; k ++)
    {
       double Val = Coeff_h0 * (double)get_pix(Old,i,j,k) 
	 + Coeff_h1 * ((double)get_pix(Old,i,j,k-Step) 
	 + (double)get_pix(Old,i,j,k+Step))
	 + Coeff_h2 * ((double)get_pix(Old,i,j,k-2*Step) 
	 + (double)get_pix(Old,i,j,k+2*Step));
       Buff1(i,j,k) = (float) Val;
    }
 
	//convolution in Y axe
    for (i = 0; i < Nx; i ++)
    for (j = 0; j < Ny; j ++)
    for (k = 0; k < Nz; k ++)
    {
       double Val = Coeff_h0 * (double)get_pix(Buff1,i,j,k) 
	 + Coeff_h1 * ((double)get_pix(Buff1,i,j-Step,k) 
	 + (double)get_pix(Buff1,i,j+Step,k))
	 + Coeff_h2 * ((double)get_pix(Buff1,i,j-2*Step,k) 
	 + (double)get_pix(Buff1,i,j+2*Step,k));
       Buff2(i,j,k) = (float) Val;
    }


	//convolution in X axe
    for (i = 0; i < Nx; i ++)
    for (j = 0; j < Ny; j ++)
    for (k = 0; k < Nz; k ++)
    {
       New(i,j,k) = Coeff_h0 * (double)get_pix(Buff2,i,j,k) 
	 + Coeff_h1 * ((double)get_pix(Buff2,i-Step,j,k) 
	 + (double)get_pix(Buff2,i+Step,j,k))
	 + Coeff_h2 * ((double)get_pix(Buff2,i-2*Step,j,k) 
	 + (double)get_pix(Buff2,i+2*Step,j,k));
    }

}

// ***************************************************************************

static float sign(float a)
{
	return float((a>0)-(a<0));
}

static void stabilize(fltarray & A, int s)
{
#ifdef _STAB_PATROU3D_
//	cerr<<"stab scale "<<s<<endl;
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = VST_b[s] * sign(A(i,j,k)+VST_c[s]) * sqrt(abs(A(i,j,k)+VST_c[s]));
//				A(i,j,k) = A(i,j,k)*(s+1)+(s+1);
			}
#endif
}
static void unstabilize(fltarray & A, int s)
{
#ifdef _STAB_PATROU3D_
//	cerr<<"unstab scale "<<s<<endl;
	for(int i=0;i<A.nx();i++)
		for(int j=0;j<A.ny();j++)
			for(int k=0;k<A.nz();k++)
			{
				A(i,j,k) = sign(A(i,j,k))*A(i,j,k)*A(i,j,k)/(VST_b[s]*VST_b[s]) - VST_c[s];
//				A(i,j,k) = A(i,j,k).real()/(s+1)+(s+1);
			}
#endif
}

/****************************************************************************/
//	A trous transformation 
void POISSONATROUS_3D_WT::transform(fltarray &Cube, fltarray * & Wavelet_Coef, int NbrScale)
{
	_NbrScale = NbrScale;
	_Nx = Cube.nx();
	_Ny = Cube.ny();
	_Nz = Cube.nz();

	Wavelet_Coef[0] = Cube;
	stabilize(Wavelet_Coef[0],0);
	fltarray Data_out;

	//cout << "TRANS: " << Cube.min() << " " << Cube.max()<< endl;     
	for (int s = 0; s <  NbrScale-1; s++)
	{
	//compute the signal in Vs+1 space:  Signal(Vs+1)=h*Signal(Vs)
		unstabilize(Wavelet_Coef[s],s);
		b3spline_filtering(Wavelet_Coef[s], Wavelet_Coef[s+1], s);
		stabilize(Wavelet_Coef[s],s);
		stabilize(Wavelet_Coef[s+1],s+1);
		Wavelet_Coef[s] = Wavelet_Coef[s] - Wavelet_Coef[s+1];
	}
	
	unstabilize(Wavelet_Coef[NbrScale-1],NbrScale-1);
	if(no_fine)
		Wavelet_Coef[0].init(0.0);
}

/****************************************************************************/

void POISSONATROUS_3D_WT::recons(fltarray * & Wavelet_Coef, fltarray & Data_Out, 
                            int NbrScale, Bool AddLastScale)
{
//	if(Verbose) cerr<<"POISSONATROUS_3D_WT::recons(.,.,"<<NbrScale<<","<<AddLastScale<<")..."<<endl;
	
	int Nx = Wavelet_Coef[0].nx();
	int Ny = Wavelet_Coef[0].ny();
	int Nz = Wavelet_Coef[0].nz(); 
	int s;
	if(NbrScale>0) _NbrScale=NbrScale;
	
	Data_Out.alloc(Nx,Ny,Nz);
	
	int Last_Scale_Used= (AddLastScale == True) ? _NbrScale : _NbrScale-1;
	Data_Out = Wavelet_Coef[0];	
	for (s = 1; s < Last_Scale_Used; s++) Data_Out += Wavelet_Coef[s];
   
//	if(Verbose) cerr<<"POISSONATROUS_3D_WT::recons(.,.,"<<_NbrScale<<","<<AddLastScale<<")..."<<endl;
}

/****************************************************************************/
 
void POISSONATROUS_3D_WT::normalize_self(fltarray *TabBand, bool inverse)
{
	if(Verbose) cerr<<"POISSONATROUS_3D_WT::normalize_self(.,"<<inverse<<")..."<<endl;
	
	int Nx = (TabBand[0]).nx();
	int Ny = (TabBand[0]).ny();
	int Nz = (TabBand[0]).nz();
		
	for (int s=0; s < _NbrScale-1; s++)
	{
		float sig = VST_norm[s] ;
		for (int i=0; i < Nx; i++)
			for (int j=0; j < Ny; j++)
				for (int k=0; k < Nz; k++)
					if(inverse)
						(TabBand[s])(i,j,k) *= sig ;
					else
						(TabBand[s])(i,j,k) /= sig ;
	}
						
	if(Verbose) cerr<<"End POISSONATROUS_3D_WT::normalize_self"<<endl;
}

/****************************************************************************/
 
void POISSONATROUS_3D_WT::threshold(fltarray * & TabBandIn, float thresh, bool soft, bool normalized)
{
	if(Verbose) cerr<<"POISSONATROUS_3D_WT::threshold(.,"<<thresh<<", "<<soft<<")..."<<endl;
	
	int Nx = (TabBandIn[0]).nx();
	int Ny = (TabBandIn[0]).ny();
	int Nz = (TabBandIn[0]).nz();
	
	if(!soft) // hard
	{
		for (int s=0; s < _NbrScale-1; s++)
		{
			int cnt=0;
			float lvl = thresh;
			if(!normalized)
				lvl *= VST_norm[s] ;
			for (int i=0; i < Nx; i++)
				for (int j=0; j < Ny; j++)
					for (int k=0; k < Nz; k++)
						if(abs((TabBandIn[s])(i,j,k)) < lvl)
						{
							(TabBandIn[s])(i,j,k) = 0 ;
							cnt++;
						}
			if(Verbose) cerr<<" threshold scale "<<s<<"("<<lvl<<") : "<<cnt<<"#->"<<(float(Nx*Ny*Nz)-cnt)/float(Nx*Ny*Nz)<<endl;
		}
	}
	else // soft
	{
		for (int s=0; s < _NbrScale-1; s++)
		{
			int cnt=0;
			float lvl = thresh;
			if(!normalized)
				lvl *= VST_norm[s] ;
			for (int i=0; i < Nx; i++)
				for (int j=0; j < Ny; j++)
					for (int k=0; k < Nz; k++)
						if( abs((TabBandIn[s])(i,j,k)) < lvl )
						{
							(TabBandIn[s])(i,j,k)=0;
							cnt++;
						}
						else (TabBandIn[s])(i,j,k) -= (2*int( (TabBandIn[s])(i,j,k)>0 )-1)*lvl;
			if(Verbose) cerr<<" threshold scale "<<s<<"("<<lvl<<") : "<<cnt<<"#->"<<(float(Nx*Ny*Nz)-cnt)/float(Nx*Ny*Nz)<<endl;
		}
	}
	if(Verbose) cerr<<"End POISSONATROUS_3D_WT::threshold"<<endl;
}

// ********************************************************************

void POISSONATROUS_3D_WT::clean_single(fltarray * & TabBandIn, float Sigma=0.)
{
//	if(Verbose) cerr<<"atrou::clean_single(.)..."<<endl;
	
	float zero=Sigma/10.;
	for (int s=0; s < _NbrScale-1; s++)
		for (int i=0; i < _Nx; i++)
			for (int j=0; j < _Ny; j++)
				for (int k=0; k < _Nz; k++)
					if(abs((TabBandIn[s])(i,j,k)) > zero)
					{
						bool single=true;
						
						int ii=-1;
						
						while ((single==true) && (ii<2))
						{
							int jj=-1;
							while ((single==true) && (jj<2))
							{
								int kk=-1;
								while ((single==true) && (kk<2))
								{
									if(!(ii==0 && jj==0 && kk==0))
										if((i+ii>0)		&& (j+jj>0)		&& (k+kk>0)
										&& (i+ii<_Nx-1) && (j+jj<_Ny-1) && (k+kk<_Nz-1) )
											if( abs((TabBandIn[s])(i+ii,j+jj,k+kk)) > zero ) single=false;
									kk++;
								}
								jj++;
							}
							ii++;
						}
						if(single) (TabBandIn[s])(i,j,k) = 0;
					}
//	if(Verbose) cerr<<"...end atrou::clean_single"<<endl;
}

// ********************************************************************

void POISSONATROUS_3D_WT::extract_stat(fltarray * & TabBand, char *Outname)
{
//	if(Verbose) cerr<<"atrou::extract_stat(.,"<<Outname<<")..."<<endl;
	
	TabStat.alloc(_NbrScale,6);
	for(int s=0;s<_NbrScale;s++)
	{
		double m1=0,m2=0,m3=0,m4=0;
		float val;
		float Minni = (TabBand[s])(0,0,0);
		float Maxxi = Minni;
		
		for (int i=0; i < _Nx; i++)
			for (int j=0; j < _Ny; j++)
				for (int k=0; k < _Nz; k++)
				{
					val = (TabBand[s])(i,j,k) / VST_norm[s];
				// Moments
					m1 += val;
					m2 += pow( val , 2 );
					m3 += pow( val , 3 );
					m4 += pow( val , 4 );
				}
				
		int N = _Nx*_Ny*_Nz;

		TabStat(s,0) = m1/N;
		moment4_center(N, TabStat(s,0), m2/N, m3/N, m4/N, TabStat(s,1), TabStat(s,2), TabStat(s,3));
		TabStat(s,4) = Minni;
		TabStat(s,5) = Maxxi;
		
		if(Verbose)
			cerr << s <<"\t"<< TabStat(s,0) <<"\t"<< TabStat(s,1) <<"\t"<< TabStat(s,2) <<"\t"<<
					 TabStat(s,3) <<"\t"<< TabStat(s,4) << "\t" << TabStat(s,5) << endl;
	}
	
// Write the stat file
	if(Outname!=NULL) 
	{
		char filename[64];
		sprintf(filename,"%s_stat.fits",Outname);
		writefltarr(filename, TabStat);
	}

//	if(Verbose) cerr<<"...end atrou::extract_stat"<<endl;
}

/****************************************************************************/

static inline void PrintError( int status)
{
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
void POISSONATROUS_3D_WT::write(char *Name, fltarray * TabBand, bool Normalize)
{
//	if(Verbose) cerr<<"atrou::write("<<Name<<",.,)"<<endl;
	
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
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D A Trou wavelet transform", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NScale3D", (long) _NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) _Nx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) _Ny, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) _Nz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	for (int b=0; b < _NbrScale; b++)
	{
		naxis=3;
		naxes[0] = TabBand[b].nx();
		naxes[1] = TabBand[b].ny();
		naxes[2] = TabBand[b].nz();
//	cerr<<"save : "<<s<<" "<<naxes[0]<<" "<<naxes[1]<<" "<<naxes[2]<<endl;
		if(ffcrhd(fptr,&status))
			PrintError( status );
		if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
			PrintError( status );

	// save the data
		if ( ffppre(fptr, group, 1, naxes[0]*naxes[1]*naxes[2], (TabBand[b]).buffer(), &status) )
			PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
//	if(Verbose) cerr<<"...end atrou::write"<<endl;
}

// ********************************************************************

void POISSONATROUS_3D_WT::read(char *Name, fltarray * &TabBand, bool *NormalizeInv)
{
//	if(Verbose) cerr<<"atrou::read("<<Name<<",.,.,.)"<<endl;
	bool LocVerbose = false & Verbose;
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
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NScale3D", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*NormalizeInv = (bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
// TabBand allocation
	alloc(TabBand, _Nx, _Ny, _Nz, _NbrScale);
	
// Read data
	int cnt=1;
	for (int b=0; b < _NbrScale; b++)
	{
		if(LocVerbose) cerr<<" read : "<<b <<endl;
		int NX,NY,NZ;
		if (fits_movabs_hdu(fptr, ++cnt, NULL, &status)) PrintError( status );

		if (ffgkyj(fptr,(char*)"NAXIS1", &mon_long, comment, &status)) PrintError( status );
		NX = (int)mon_long;
		if (ffgkyj(fptr,(char*)"NAXIS2", &mon_long, comment, &status)) PrintError( status );
		NY = (int)mon_long;
		if (ffgkyj(fptr,(char*)"NAXIS3", &mon_long, comment, &status)) PrintError( status );
		NZ = (int)mon_long;

		if(LocVerbose) cerr<<" read TB("<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;

		TabBand[b].alloc(NX,NY,NZ);
		if (ffgpve(fptr, 1, 1, NX*NY*NZ, nulval, (TabBand[b]).buffer(), &anynul, &status)) PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
//	if(Verbose) cerr<<"...end atrou::read"<<endl;
}

// ********************************************************************
// ********************************************************************
//					Global functions - for MEX
// ********************************************************************
// ********************************************************************

void piwt3d_clear(vector< fltarray* > &C)
{
	// coresponds to POISSONATROUS_3D_WT::alloc 
	delete [] C[0];
}	

void piwt3d_transform(fltarray &Data, vector< fltarray* > &vTabBand, int _NbrScale)
{
//cerr<<"iwt3d_transform(fltarray &Data, vector< fltarray* > &vTabBand, "<<_NbrScale<<endl;
// iwt allocation
	fltarray* TabBand;
	POISSONATROUS_3D_WT *atrou = new POISSONATROUS_3D_WT();
	
	atrou->alloc(TabBand, Data.nx(), Data.ny(), Data.nz(), _NbrScale);
	
// Calculus
	atrou->transform(Data,TabBand,_NbrScale);
//	atrou->normalize_self(TabBand,false);
	
// vTabBand allocation
	vTabBand.resize(_NbrScale);
	for(int s=0;s<_NbrScale;s++)
		vTabBand[s] = &TabBand[s];
	
	delete atrou;
	return ;
}

void piwt3d_recons(vector< fltarray* > &vTabBand, fltarray &Data)
{
// Number of scales
	int _NbrScale = vTabBand.size();
	
// iwt allocation
	fltarray* TabBand;
	POISSONATROUS_3D_WT *atrou = new POISSONATROUS_3D_WT();
	atrou->set_nbr_sale(_NbrScale);
	int nx=vTabBand[0]->nx(), ny=vTabBand[0]->ny(), nz=vTabBand[0]->nz();

// TabBand allocation
	TabBand = new fltarray [_NbrScale];
	for(int s=0;s<_NbrScale;s++)
		TabBand[s].alloc(vTabBand[s]->buffer(), nx, ny, nz);
	
// Calculus
//	atrou->normalize_self(TabBand,true);
	atrou->recons(TabBand,Data,_NbrScale);
//	atrou->normalize_self(TabBand,false);
	
	delete [] TabBand;
	delete atrou;
}

void piwt3d_threshold(vector< fltarray* > &vTabBand, float threshold, filter_type FilterType)
{
// Number of scales
	int _NbrScale = vTabBand.size();
	
// iwt allocation
	fltarray* TabBand;
	POISSONATROUS_3D_WT *atrou = new POISSONATROUS_3D_WT();
	atrou->set_nbr_sale(_NbrScale);
	int nx=vTabBand[0]->nx(), ny=vTabBand[0]->ny(), nz=vTabBand[0]->nz();
	
// TabBand allocation
	TabBand = new fltarray [_NbrScale];
	for(int s=0;s<_NbrScale;s++)
		TabBand[s].alloc(vTabBand[s]->buffer(), nx, ny, nz);

// Calculus
	if(FilterType==FT_HARD || FilterType==FT_SOFT) atrou->threshold(TabBand, threshold, FilterType==FT_SOFT, true);// true=normalized
	else cerr<<"Unknown filtering method"<<endl;
	
	delete [] TabBand;
	delete atrou;
}

void piwt3d_filter(fltarray &Data, fltarray &Recons, int _NbrScale, float threshold, filter_type FilterType)
{
// iwt allocation
	fltarray* TabBand;
	POISSONATROUS_3D_WT *atrou = new POISSONATROUS_3D_WT();
	atrou->alloc(TabBand, Data.nx(), Data.ny(), Data.nz(), _NbrScale);

// Calculus
	atrou->transform(Data,TabBand,_NbrScale);
	if(FilterType==FT_HARD || FilterType==FT_SOFT) atrou->threshold(TabBand, threshold, FilterType==FT_SOFT, false);//false=!normalized
	else cerr<<"Unknown filtering method"<<endl;
	atrou->recons(TabBand,Recons,_NbrScale);
	
	delete [] TabBand;// allocated by atrou->alloc -> POISSONATROUS_3D_WT::alloc
	delete atrou;
}

