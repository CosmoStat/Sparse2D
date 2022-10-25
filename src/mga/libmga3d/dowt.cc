
#include "dowt.h"

extern bool Verbose;

static const double PI4 = M_PI/4.0;
static const double PI2 = M_PI/2.0;

DOWT::DOWT(SubBandFilter *SB) : Ortho_3D_WT(*SB)
{
}

DOWT::~DOWT()
{
	
}

void DOWT::init(int nx, int ny, int nz, int _NbrScale)
{
	DataNx = nx;
	DataNy = ny;
	DataNz = nz;
	NbrScale = _NbrScale;
}

// ********************************************************************
// Order of the bands : high frequencies in ...
// 0 : xy
// 1 : xz
// 2 : xyz
// 3 : x
// 4 : yz
// 5 : z
// 6 : zy
// High frequency :
// -> x : <4
// -> y : !%2
// -> z : %3

// ********************************************************************
// Order of the bands : same as the undecimated case -> high frequencies in ...
// 0 : xyz
// 1 : xy
// 2 : xz
// 3 : x
// 4 : yz
// 5 : y
// 6 : z
// High frequency :
// -> x : <4
// -> y : /2!%2
// -> z : !%2

int DOWT::size_band_nx(int S, int b)
{
	int nxfine, xfine = DataNx;
	for(int s=0;s<S+1;s++)
	{
		nxfine = xfine;
		xfine = (int) floor((xfine+1)/2);
		nxfine -= xfine;
	}
	if(S==NbrScale-1)
		return nxfine+xfine;
	else
		return b<4 ? nxfine : xfine ;
}
int DOWT::start_band_nx(int S, int b)
{
	int nxfine, xfine = DataNx;
	for(int s=0;s<S+1;s++)
	{
		nxfine = xfine;
		xfine = (int) floor((xfine+1)/2);
		nxfine -= xfine;
	}
	if(S==NbrScale-1)
		return 0;
	else
		return b<4 ? xfine  : 0 ;
}
int DOWT::size_band_ny(int S, int b)
{
	int nyfine, yfine = DataNy;
	for(int s=0;s<S+1;s++)
	{
		nyfine = yfine;
		yfine = (int) floor((yfine+1)/2);
		nyfine -= yfine;
	}
	if(S==NbrScale-1)
		return nyfine+yfine;
	else
		return (b/2)%2	? yfine  : nyfine ;
}
int DOWT::start_band_ny(int S, int b)
{
	int nyfine, yfine = DataNy;
	for(int s=0;s<S+1;s++)
	{
		nyfine = yfine;
		yfine = (int) floor((yfine+1)/2);
		nyfine -= yfine;
	}
	if(S==NbrScale-1)
		return 0;
	else
		return (b/2)%2 ? 0 	 : yfine ;
}
int DOWT::size_band_nz(int S, int b)
{
	int nzfine,zfine = DataNz;
	for(int s=0;s<S+1;s++)
	{
		nzfine = zfine;
		zfine = (int) floor((zfine+1)/2);
		nzfine -= zfine;
	}
	if(S==NbrScale-1)
		return nzfine+zfine;
	else
		return b%2 ? zfine : nzfine ;
}
int DOWT::start_band_nz(int S, int b)
{
	int nzfine,zfine = DataNz;
	for(int s=0;s<S+1;s++)
	{
		nzfine = zfine;
		zfine = (int) floor((zfine+1)/2);
		nzfine -= zfine;
	}
	if(S==NbrScale-1)
		return 0;
	else
		return b%2 ? 0 : zfine ;
}

// ********************************************************************

void DOWT::threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"DOWT::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;
	
	float lvl = SigmaNoise * NSigma ;
	if(Verbose) cerr<<"Thresholding at lvl="<<lvl<<endl;
	int xcoarse = DataNx;
	int ycoarse = DataNy;
	int zcoarse = DataNz;

	for(int i=0;i<NbrScale-1;i++)
	{
		xcoarse =(int) floor((xcoarse+1)/2);
		ycoarse =(int) floor((ycoarse+1)/2);
		zcoarse =(int) floor((zcoarse+1)/2);
	}

	if(FilterType==FT_HARD)
	{
		for (int i=0; i < DataNx; i++)
			for (int j=0; j < DataNy; j++)
				for (int k=0; k < DataNz; k++)
					if(i>=xcoarse || j>=ycoarse || k>=zcoarse)
						if(abs(TabBand(i,j,k)) < lvl) TabBand(i,j,k) = 0;
	}
	if(FilterType==FT_SOFT)
		for (int i=0; i < DataNx; i++)
			for (int j=0; j < DataNy; j++)
				for (int k=0; k < DataNz; k++)
					if(i>=xcoarse || j>=ycoarse || k>=zcoarse)
						TabBand(i,j,k) = ( abs(TabBand(i,j,k)) < lvl ) ? 0 : TabBand(i,j,k) - (2*int(TabBand(i,j,k)>0)-1)*lvl ;
	
	if(Verbose) cerr<<"...End DOWT::threshold"<<endl;
}

/*********************************************************************/

void DOWT::wiener(fltarray &TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"DOWT::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
	
	float val;
	float noise2 = noise_lvl*noise_lvl;
	
	int xfine = DataNx;
	int yfine = DataNy;
	int zfine = DataNz;
	
	for(int s=0;s<NbrScale-1;s++)
	{
		int nxfine = xfine;
		int nyfine = yfine;
		int nzfine = zfine;
		xfine = (int) floor((xfine+1)/2);
		yfine = (int) floor((yfine+1)/2);
		zfine = (int) floor((zfine+1)/2);
		nxfine -= xfine;
		nyfine -= yfine;
		nzfine -= zfine;
		
		for(int b=0;b<7;b++)
		{
			// discriminate the 7 directions with <4 /2!%2 !%2 operators
			int nx = 	b<4 ? nxfine : xfine ;
			int xband = b<4 ? xfine  : 0 ;
			int ny = 	b%2	? yfine  : nyfine ;
			int yband = b%2 ? 0 	 : yfine ;
			int nz = 	b%3 ? nzfine : zfine ;
			int zband = b%3 ? zfine  : 0 ;
			// cerr<<"band ("<<s<<","<<b<<") x:x+nx "<<xband<<":"<<xband+nx-1<<" "<<yband<<":"<<yband+ny-1<<" "<<zband<<":"<<zband+nz-1<<endl;

			int Nx = ceil(float(nx)/float(LocalBS));
			int Ny = ceil(float(ny)/float(LocalBS));
			int Nz = ceil(float(nz)/float(LocalBS));

			fltarray coef_wiener(nx,ny,nz);
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
									int x = (bx + kx*LocalBS) % nx;
									int y = (by + ky*LocalBS) % ny;
									int z = (bz + kz*LocalBS) % nz;
									val = (TabBand)(xband+x,yband+y,zband+z);
									sigma2+=pow(val,2);
								}

						float sig = max( 0.0, sigma2/cnt - noise2 );
						float norm = sig / (sig+noise2);


					// Store the coef in the table
						for(int bx = 0; bx < LocalBS; bx++)
							for(int by = 0; by < LocalBS; by++)
								for(int bz = 0; bz < LocalBS; bz++)
								{
									int x = (bx + kx*LocalBS) % nx;
									int y = (by + ky*LocalBS) % ny;
									int z = (bz + kz*LocalBS) % nz;
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
									int x = (bx + kx*LocalBS) % nx;
									int y = (by + ky*LocalBS) % ny;
									int z = (bz + kz*LocalBS) % nz;
									(TabBand)(xband+x,yband+y,zband+z) *= coef_wiener(x,y,z);
								}
					}
		}// end band
	}// end scale

	if(Verbose) cerr<<"...End DOWT::wiener"<<endl;
}

/*********************************************************************/

void DOWT::fdr(fltarray &TabBand, float Alpha, float SigmaNoise)
{
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
//	bool LocVerbose = true & Verbose;
	

	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

/****************************************************************************/

void DOWT::stein_block_threshold(fltarray &TabBand, float SigmaNoise)
{
	if(Verbose) cerr<<"DOWT::stein_block_threshold(.)"<<endl;
	//bool LocVerbose = true & Verbose;

	if(Verbose) cerr<<"...End DOWT::stein_block_threshold"<<endl;
}

/*********************************************************************/

void DOWT::write(char *Name, fltarray &TabBand, bool Normalize)
{
	if(Verbose) cerr<<"DOWT::write("<<Name<<",.,"<<Normalize<<")"<<endl;

	if(Verbose) cerr<<"...end DOWT::write"<<endl;
}

// ********************************************************************

void DOWT::read(char *Name, fltarray &TabBand, bool *NormalizeInv)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;


	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}




// ********************************************************************
// ********************************************************************
// GLOBAL DOWT functions
// ********************************************************************
// ********************************************************************

void coef_import(DOWT *dwt, vector< vector< fltarray* > > &vTabBand, fltarray &TabBand)
{
	int _NbrScale = vTabBand.size();
	for(int s=0;s<_NbrScale;s++)
	{
		int nb = (s==_NbrScale-1) ? 1:7;
		for(int b=0;b<nb;b++)
		{
			int nx = dwt->size_band_nx(s,b);
			int ny = dwt->size_band_ny(s,b);
			int nz = dwt->size_band_nz(s,b);
			int x0 = dwt->start_band_nx(s,b);
			int y0 = dwt->start_band_ny(s,b);
			int z0 = dwt->start_band_nz(s,b);
			for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++)
				TabBand(i+x0,j+y0,k+z0) = (*vTabBand[s][b])(i,j,k);
		}
	}
}

void coef_export(DOWT *dwt, fltarray &TabBand, vector< vector< fltarray* > > &vTabBand, int _NbrScale, bool alloc)
{ 
	if(alloc) vTabBand.resize(_NbrScale);
	for(int s=0;s<_NbrScale;s++)
	{
		int nb = (s==_NbrScale-1) ? 1:7;
		if(alloc) vTabBand[s].resize(nb);
		for(int b=0;b<nb;b++)
		{
			int nx = dwt->size_band_nx(s,b);
			int ny = dwt->size_band_ny(s,b);
			int nz = dwt->size_band_nz(s,b);
			int x0 = dwt->start_band_nx(s,b);
			int y0 = dwt->start_band_ny(s,b);
			int z0 = dwt->start_band_nz(s,b);
			if(alloc) vTabBand[s][b] = new fltarray(nx,ny,nz);
			for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++)
				(*vTabBand[s][b])(i,j,k) = TabBand(i+x0,j+y0,k+z0);
		}
	}
}

void dwt3d_clear(vector< vector< fltarray* > > &C)
{// coherent with coef_export
	for(int i=0;i<C.size();i++)
	{
		for(int j=0;j<C[i].size();j++)
			delete C[i][j];
		C[i].clear();
	}
	C.clear();
}

void dwt3d_transform(fltarray &Data, fltarray &TabBand, int NbrScale3D, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	dwt->init(Data.nx(), Data.ny(), Data.nz(), NbrScale3D);
	dwt->transform(Data,TabBand,NbrScale3D);
	
	delete dwt;
	delete SB1D;
}

void dwt3d_recons(fltarray &TabBand, fltarray &Data, int NbrScale3D, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	dwt->init(Data.nx(), Data.ny(), Data.nz(), NbrScale3D);
	dwt->recons(TabBand,Data,NbrScale3D);
	
	delete dwt;
	delete SB1D;
}

void dwt3d_transform(fltarray &Data, vector< vector< fltarray* > > &vTabBand, int _NbrScale, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
// dwt allocation
	fltarray TabBand;
	DOWT *dwt = new DOWT(SB1D);
	
// Calculus
	dwt->init(Data.nx(), Data.ny(), Data.nz(), _NbrScale);
	dwt->transform(Data, TabBand, _NbrScale);
	
// vTabBand allocation and recopy
	coef_export(dwt,TabBand,vTabBand,_NbrScale,true);
				
	delete dwt;
}

void dwt3d_recons(vector< vector< fltarray* > > &vTabBand, fltarray &Data, int wavelet_type)
{
	int _NbrScale = vTabBand.size();
	
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	int Nx=vTabBand[_NbrScale-1][0]->nx(),Ny=vTabBand[_NbrScale-1][0]->ny(),Nz=vTabBand[_NbrScale-1][0]->nz();
	for(int s=0;s<_NbrScale-1;s++)
	{
		Nx+=vTabBand[s][2]->nx();
		Ny+=vTabBand[s][2]->ny();
		Nz+=vTabBand[s][2]->nz();
	}
	dwt->init(Nx,Ny,Nz,_NbrScale);
	fltarray TabBand(Nx,Ny,Nz);
	
// Coefficients recopy
	coef_import(dwt,vTabBand,TabBand);
	
	dwt->recons(TabBand, Data, _NbrScale);
	
	delete dwt;
	delete SB1D;
}

void dwt3d_threshold(fltarray &TabBand, int NbrScale3D, float threshold, filter_type FilterType)
{
// Wavelet initialisation (necessary but unused)
	type_sb_filter wavelet_type = F_MALLAT_7_9;
	FilterAnaSynt SelectFilter(wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	dwt->init(TabBand.nx(), TabBand.ny(), TabBand.nz(),NbrScale3D);
	if(FilterType==FT_HARD || FilterType==FT_SOFT) dwt->threshold(TabBand, threshold, 1, FilterType, false);
	else if(FilterType==FT_WIENER) dwt->wiener(TabBand, threshold, 3);
	else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet"<<endl;
	
	delete dwt;
	delete SB1D;
}

void dwt3d_threshold(vector< vector< fltarray* > > &vTabBand, float threshold, filter_type FilterType)
{
	int _NbrScale = vTabBand.size();
	
// Wavelet initialisation (necessary but unused)
	type_sb_filter wavelet_type = F_MALLAT_7_9;
	FilterAnaSynt SelectFilter(wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	int Nx=vTabBand[_NbrScale-1][0]->nx(),Ny=vTabBand[_NbrScale-1][0]->ny(),Nz=vTabBand[_NbrScale-1][0]->nz();
	for(int s=0;s<_NbrScale-1;s++)
	{
		Nx+=vTabBand[s][2]->nx();
		Ny+=vTabBand[s][2]->ny();
		Nz+=vTabBand[s][2]->nz();
	}
	dwt->init(Nx,Ny,Nz,_NbrScale);
	fltarray TabBand(Nx,Ny,Nz);
	
// Coefficients import
	coef_import(dwt,vTabBand,TabBand);
	
// Thresholding
	if(FilterType==FT_HARD || FilterType==FT_SOFT) dwt->threshold(TabBand, threshold, 1, FilterType, false);
	else if(FilterType==FT_WIENER) dwt->wiener(TabBand, threshold, 3);
	else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet"<<endl;
	
// Coefficients export
	coef_export(dwt,TabBand,vTabBand,_NbrScale,false);
	
	delete dwt;
	delete SB1D;
}


void dwt3d_filter(fltarray &Data, fltarray &Recons, int NbrScale3D, float threshold, filter_type FilterType, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	SB1D->Border = I_PERIOD;
	
	DOWT *dwt = new DOWT(SB1D);
	fltarray TabBand;
	dwt->init(Data.nx(), Data.ny(), Data.nz(),NbrScale3D);
	dwt->transform(Data,TabBand,NbrScale3D);
	if(FilterType==FT_HARD || FilterType==FT_SOFT) dwt->threshold(TabBand, threshold, 1, FilterType, false);
	else if(FilterType==FT_WIENER) dwt->wiener(TabBand, threshold, 3);
	else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet"<<endl;
	dwt->recons(TabBand,Recons,NbrScale3D);
	
	delete dwt;
	delete SB1D;
}


