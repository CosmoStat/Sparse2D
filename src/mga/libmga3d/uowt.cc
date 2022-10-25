
#include "uowt.h"

extern bool Verbose;

static const double PI4 = M_PI/4.0;
static const double PI2 = M_PI/2.0;

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

UOWT::UOWT(SubBandFilter *SB) : PAVE_3D_WT(*SB)
{
// Filtering properties
	no_fine=false;
	no_coarse=false;
}

UOWT::~UOWT()
{
	
}

void UOWT::init(fltarray * &TabBand, int nx, int ny, int nz, int _NbrScale)
{
	DataNx = nx;
	DataNy = ny;
	DataNz = nz;
	NbrScale = _NbrScale;
	PAVE_3D_WT::alloc(TabBand, nx, ny, nz, NbrScale);
}

// ********************************************************************

void UOWT::apply_positivity(fltarray & Data)
{
	if(positivity)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) < 0 ) Data(i,j,k) = 0;
}

// ********************************************************************

void UOWT::extract_stat(fltarray * TabBand)
{
	for (int b=0;b<get_nbands();b++)
	{
		cerr<<"Scale "<<b/7<<", Band "<<b%7<<" : ";
		TabBand[b].info();
	}
}

float UOWT::get_maxabs(fltarray * TabBand)
{
	float Max = 0.;
	for (int b=0;b<get_nbands();b++)
	{
		for(int i=0;i<TabBand[b].nx();i++)
		for(int j=0;j<TabBand[b].ny();j++)
		for(int k=0;k<TabBand[b].nz();k++)
			Max = (Max > abs(TabBand[b](i,j,k))) ? Max : abs(TabBand[b](i,j,k)) ;
	}
	return Max;
}

// ********************************************************************

void UOWT::threshold(fltarray * TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"UOWT::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;
	
	float lvl = SigmaNoise * NSigma ;
	if(FilterType==FT_HARD)
		for (int b=0;b<7*(NbrScale-1);b++)
			for (int i=0; i < DataNx; i++)
				for (int j=0; j < DataNy; j++)
					for (int k=0; k < DataNz; k++)
						if(abs(TabBand[b](i,j,k)) < lvl) TabBand[b](i,j,k) = 0;
	if(FilterType==FT_SOFT)
		for (int b=0;b<7*(NbrScale-1);b++)
			for (int i=0; i < DataNx; i++)
				for (int j=0; j < DataNy; j++)
					for (int k=0; k < DataNz; k++)
						TabBand[b](i,j,k) = ( abs(TabBand[b](i,j,k)) < lvl ) ? 0 : TabBand[b](i,j,k) - (2*int(TabBand[b](i,j,k)>0)-1)*lvl ;
	
	if(no_coarse)
		TabBand[7*(NbrScale-1)].init(0.0);
	
	if(no_fine)
		if(FilterType==FT_HARD)
			for (int b=0;b<7;b++)
				for (int i=0; i < DataNx; i++)
					for (int j=0; j < DataNy; j++)
						for (int k=0; k < DataNz; k++)
							TabBand[b](i,j,k) = 0;
	
	if(Verbose) cerr<<"...End UOWT::threshold"<<endl;
}

/*********************************************************************/

void UOWT::threshold_on_mask(fltarray * TabBand, fltarray * TabBand_mask, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"UOWT::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;
	
	float lvl = SigmaNoise * NSigma ;
	if(FilterType==FT_HARD)
		for (int b=0;b<7*(NbrScale-1);b++)
			for (int i=0; i < DataNx; i++)
				for (int j=0; j < DataNy; j++)
					for (int k=0; k < DataNz; k++)
//						if(abs(TabBand_mask[b](i,j,k)) > 0.5 ) // mask -1,1,0
							if(abs(TabBand[b](i,j,k)) < lvl*abs(TabBand_mask[b](i,j,k))) TabBand[b](i,j,k) = 0;
	
	if(FilterType==FT_SOFT)
		for (int b=0;b<7*(NbrScale-1);b++)
			for (int i=0; i < DataNx; i++)
				for (int j=0; j < DataNy; j++)
					for (int k=0; k < DataNz; k++)
//						if(abs(TabBand_mask[b](i,j,k)) > 0.5 ) // mask -1,1,0
							TabBand[b](i,j,k) = ( abs(TabBand[b](i,j,k)) < lvl*abs(TabBand_mask[b](i,j,k)) ) ? 0 : TabBand[b](i,j,k) - (2*int(TabBand[b](i,j,k)>0)-1)*lvl*abs(TabBand_mask[b](i,j,k)) ;
	
	if(no_coarse)
		TabBand[7*(NbrScale-1)].init(0.);
	
	if(no_fine)
		if(FilterType==FT_HARD)
			for (int b=0;b<7;b++)
				TabBand[b].init(0.);
	
	if(Verbose) cerr<<"...End UOWT::threshold"<<endl;
}

/*********************************************************************/

void UOWT::wiener(fltarray * TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"UOWT::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
	
	float val;
	float noise2 = noise_lvl*noise_lvl;
	
	for(int b=0;b<get_nbands()-1;b++)
	{
		int Nx = ceil(float(DataNx)/float(LocalBS));
		int Ny = ceil(float(DataNy)/float(LocalBS));
		int Nz = ceil(float(DataNz)/float(LocalBS));

		fltarray coef_wiener(DataNx,DataNy,DataNz);
		coef_wiener.init(-2);

		// Wiener blocks : evaluate the wiener coef
		for(int kx=0 ; kx < Nx ; kx++)
			for(int ky=0 ; ky < Ny ; ky++)
				for(int kz = 0; kz < Nz; kz++)
				{
					double sigma2 = 0.0;
					float cnt=pow((float)LocalBS,3);

				// Sigma calculation
					// Pixels in a wiener block (angle)
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								//cnt+=1;
								int x = (bx + kx*LocalBS) % DataNx;
								int y = (by + ky*LocalBS) % DataNy;
								int z = (bz + kz*LocalBS) % DataNz;
								val = (TabBand[b])(x,y,z);
								sigma2+=pow(val,2);
							}

					float sig = max( 0.0, sigma2/cnt - noise2 );
					float norm = sig / (sig+noise2);


				// Store the coef in the table
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								int x = (bx + kx*LocalBS) % DataNx;
								int y = (by + ky*LocalBS) % DataNy;
								int z = (bz + kz*LocalBS) % DataNz;
								if( coef_wiener(x,y,z) < -1 )
									coef_wiener(x,y,z) = norm;
							}
				}

	// Wiener blocks : Apply the coefficients
		for(int kx=0 ; kx < Nx ; kx++)
			for(int ky=0 ; ky < Ny ; ky++)
				for(int kz = 0; kz < Nz; kz++)
				{
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								int x = (bx + kx*LocalBS) % DataNx;
								int y = (by + ky*LocalBS) % DataNy;
								int z = (bz + kz*LocalBS) % DataNz;
								(TabBand[b])(x,y,z) *= coef_wiener(x,y,z);
							}
				}
	}// end band
	
	if(Verbose) cerr<<"...End UOWT::wiener"<<endl;
}

/*********************************************************************/

void UOWT::fdr(fltarray * TabBand, float Alpha, float SigmaNoise)
{
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
//	bool LocVerbose = true & Verbose;
	

	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

/****************************************************************************/

void UOWT::stein_block_threshold(fltarray * TabBand, float SigmaNoise)
{
	if(Verbose) cerr<<"UOWT::stein_block_threshold(.)"<<endl;
	//bool LocVerbose = true & Verbose;

	if(Verbose) cerr<<"...End UOWT::stein_block_threshold"<<endl;
}

/*********************************************************************/

void UOWT::write(char *Name, fltarray * TabBand)
{
	if(Verbose) cerr<<"UOWT::write("<<Name<<",.,)"<<endl;
	
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
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D Undecimated orthogonal wavelet transform", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrBands", (long) get_nbands(), (char*)"Total number of bands", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) DataNx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) DataNy, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) DataNz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	for (int b=0; b < get_nbands(); b++)
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
	
	if(Verbose) cerr<<"...end UOWT::write"<<endl;
}

// ********************************************************************

void UOWT::read(char *Name, fltarray * &TabBand)
{
	if(Verbose) cerr<<"UOWT::read("<<Name<<",.,.,.)"<<endl;
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
	int _NbrScale, _Nx, _Ny, _Nz;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
//	init the structure
	init(TabBand, _Nx, _Ny, _Nz, _NbrScale);
	
// TabBand allocation
	TabBand = new fltarray[get_nbands()];
	
// Read data
	int cnt=1;
	for (int b=0; b < get_nbands(); b++)
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
		//TabCF_Band[b].alloc(NX,NY,NZ);
		if (ffgpve(fptr, 1, 1, NX*NY*NZ, nulval, (TabBand[b]).buffer(), &anynul, &status)) PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end UOWT::read"<<endl;
}


// ********************************************************************
// ********************************************************************
//					Floating functions - for Fista
// ********************************************************************
// ********************************************************************

void UOWT::alloc(fltarray& in, int _NbrScale)
{
	NbrScale=_NbrScale; 
	DataNx = in.nx(); 
	DataNy = in.ny(); 
	DataNz = in.nz();
	m_TabBand = new fltarray[get_nbands()];
}

void UOWT::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc)
		out = new float [get_nbands()*get_sizebands()];
	
	for(int i=0;i<get_nbands();i++)
		m_TabBand[i].alloc(out+i*get_sizebands(), DataNx, DataNy, DataNz);
	
	transform(in, m_TabBand);
}

void UOWT::recons(float* in, fltarray& out)
{
	for(int i=0;i<get_nbands();i++)
		m_TabBand[i].alloc(in+i*get_sizebands(), DataNx, DataNy, DataNz);
	
	recons(m_TabBand, out);
}

void UOWT::soft_threshold(float* in, float lvl, bool threshold_coarse)
{// threshold_coarse=true not implemented
	for(int i=0;i<get_nbands();i++)
		m_TabBand[i].alloc(in+i*get_sizebands(), DataNx, DataNy, DataNz);
	
	threshold(m_TabBand, lvl, 1, FT_SOFT, false);
}

void UOWT::hard_threshold(float* in, float lvl, bool threshold_coarse)
{// threshold_coarse=true not implemented
	for(int i=0;i<get_nbands();i++)
		m_TabBand[i].alloc(in+i*get_sizebands(), DataNx, DataNy, DataNz);
	
	threshold(m_TabBand, lvl, 1, FT_HARD, false);
}


// ********************************************************************
// ********************************************************************
//					Global functions - for MEX
// ********************************************************************
// ********************************************************************

void uwt3d_clear(vector< vector< fltarray* > > &C)
{
	// coresponds to PAVE_3D_WT::alloc 
	delete [] C[0][0];
}	

void uwt3d_transform(fltarray &Data, vector< vector< fltarray* > > &vTabBand, int _NbrScale, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
// uwt allocation
	fltarray* TabBand;
	UOWT *uwt = new UOWT(SB1D);
	uwt->init(TabBand, Data.nx(), Data.ny(), Data.nz(), _NbrScale);
	
// Calculus
	uwt->transform(Data, TabBand);
	
// vTabBand allocation
	vTabBand.resize(_NbrScale);
	for(int s=0;s<_NbrScale;s++)
	{
		int nb = (s==_NbrScale-1) ? 1:7;
		vTabBand[s].resize(nb);
		for(int b=0;b<nb;b++)
			vTabBand[s][b] = &TabBand[7*s+b];
	}
	
	delete uwt;
	return ;
}

void uwt3d_recons(vector< vector< fltarray* > > &vTabBand, fltarray &Data, int wavelet_type)
{
// Number of scales
	int _NbrScale = vTabBand.size();
	
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
// uwt allocation
	fltarray* TabBand;
	UOWT *uwt = new UOWT(SB1D);
	int nx=vTabBand[0][0]->nx(), ny=vTabBand[0][0]->ny(), nz=vTabBand[0][0]->nz();
	uwt->init(TabBand, nx, ny, nz, _NbrScale);

// TabBand allocation
	int NBands = 7*(_NbrScale-1)+1;
	for(int s=0;s<_NbrScale;s++)
	{
		int nb = (s==_NbrScale-1) ? 1:7;
		for(int b=0;b<nb;b++)
			TabBand[7*s+b].alloc(vTabBand[s][b]->buffer(), nx, ny, nz);
	}
	
// Calculus
	uwt->recons(TabBand,Data,_NbrScale);
	
	delete [] TabBand;
	delete uwt;
	delete SB1D;
}

void uwt3d_threshold(vector< vector< fltarray* > > &vTabBand, float threshold, filter_type FilterType)
{
// Number of scales
	int _NbrScale = vTabBand.size();
	
// Wavelet initialisation (unused but necessary)
	type_sb_filter wavelet_type = F_MALLAT_7_9;
	FilterAnaSynt SelectFilter(wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
// uwt allocation
	fltarray* TabBand;
	UOWT *uwt = new UOWT(SB1D);
	int nx=vTabBand[0][0]->nx(), ny=vTabBand[0][0]->ny(), nz=vTabBand[0][0]->nz();
	uwt->init(TabBand, nx, ny, nz, _NbrScale);

// TabBand allocation
	int NBands = 7*(_NbrScale-1)+1;
	TabBand = new fltarray[NBands];
	for(int s=0;s<_NbrScale;s++)
	{
		int nb = (s==_NbrScale-1) ? 1:7;
		for(int b=0;b<nb;b++)
			TabBand[7*s+b].alloc(vTabBand[s][b]->buffer(), nx, ny, nz);
	}
	
// Calculus
	if(FilterType==FT_HARD || FilterType==FT_SOFT) uwt->threshold(TabBand, threshold, 1, FilterType, false);
	else if(FilterType==FT_WIENER) uwt->wiener(TabBand, threshold, 3);
	else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet"<<endl;
	
	delete [] TabBand;
	delete uwt;
	delete SB1D;
}

void uwt3d_filter(fltarray &Data, fltarray &Recons, int _NbrScale, float threshold, filter_type FilterType, int wavelet_type)
{
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
// uwt allocation
	fltarray* TabBand;
	UOWT *uwt = new UOWT(SB1D);
	uwt->init(TabBand, Data.nx(), Data.ny(), Data.nz(), _NbrScale);

// Calculus
	uwt->transform(Data,TabBand,_NbrScale);
	if(FilterType==FT_HARD || FilterType==FT_SOFT) uwt->threshold(TabBand, threshold, 1, FilterType, false);
	else if(FilterType==FT_WIENER) uwt->wiener(TabBand, threshold, 3);
	else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet"<<endl;
	uwt ->recons(TabBand,Recons,_NbrScale);
	
	delete [] TabBand;// allocated by uwt->init -> PAVE_3D_WT::alloc
	delete uwt;
	delete SB1D;
}



