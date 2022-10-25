
#include "FCur_SFrame.h"
extern Bool Verbose;

FCur_SFrame::FCur_SFrame()
{
	TabCur=NULL;
	threshold_coarse=false;
}

FCur_SFrame::~FCur_SFrame()
{
	
}

// ********************************************************************



void FCur_SFrame::alloc_from_coarse(int _NbrScale, int Nx, int Ny, int Nz, int _NbrDir2d, Bool _ExtendWT, Bool _IsotropWT, Bool _RealCur)
{
	DataNx = Nx;
	DataNy = Ny;
	DataNt = Nz;
	NbrScale = _NbrScale;
	NbrDir2d = _NbrDir2d;

	TabCur = new FCUR;
	TabCur->Verbose = False;
	TabCur->alloc_from_coarse(NbrScale, Ny, Nx, NbrDir2d, False, False, True);

// CubeSigma initialisation
	TabCur->get_norm_coeff(3);
	
/*	char filename[128];
	sprintf(filename,"%s/common/F2DSC_TabSigma/TabSigma_n%02d_d%d.fits",getenv("CUR"),NbrScale,NbrDir2d);
	if(access(filename, F_OK)==0)
		fits_read_fltarr(filename, TabSigma);
	else
	{
		cerr<<" Error : cannot load normalisation file "<<endl<<"  "<<filename<<endl;
		//exit(-7);
	}
*/}

// ********************************************************************

void FCur_SFrame::get_norm_coeff(float NSigma)
{
//	for(int i=0;i<DataNt;i++)
//		TabCur[i].get_norm_coeff(NSigma);
}

// ********************************************************************

void FCur_SFrame::transform(fltarray &Data, Ifloat **&TabBand, bool allocTB)
{
	if(Verbose) cerr<<"FCur_SFrame::transform(.,.,"<<allocTB<<")"<<endl;
	
	float* Ptr = Data.buffer();
	Ifloat Frame(DataNy,DataNx);
	
	Ifloat NewFrame;
	for(int i=0;i<DataNt;i++)
	{
		NewFrame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		Frame += NewFrame;
	}

	for(int x=0;x<Frame.nl();x++)
	for(int y=0;y<Frame.nc();y++)
		Frame(x,y) /= float(sqrt(DataNt));
	
	TabCur->cur_trans(Frame);
	
	TabBand = new Ifloat *[NbrScale];
	for(int s=0;s<NbrScale;s++)
	{
		int NbrBands=nbr_bands(s);
		TabBand[s] = new Ifloat [NbrBands];

		for(int b=0;b<NbrBands;b++)
		{
			TabBand[s][b].alloc(TabCur->size_band_nl(s,b),TabCur->size_band_nc(s,b));
			TabCur->get_band(s, b, TabBand[s][b]);
		}
	}

	if(Verbose) cerr<<"...End FCur_SFrame::transform"<<endl;
}

// ********************************************************************

void FCur_SFrame::recons(Ifloat **&TabBand, fltarray &Recons)
{
	if(Verbose) cerr<<"FCur_SFrame::recons(.,.)"<<endl;
	
	Recons.alloc(DataNx,DataNy,DataNt);
	float* Ptr = Recons.buffer();
	Ifloat Frame(DataNy,DataNx);
	
	for(int s=0;s<NbrScale;s++)
		for(int b=0;b<nbr_bands(s);b++)
			TabCur->put_band(s, b, TabBand[s][b]);

	TabCur->cur_recons(Frame);
	
	Ifloat Recons_frame;
	for(int i=0;i<DataNt;i++)
	{
		Recons_frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		Recons_frame += Frame;
	}
	for(int i=0;i<DataNx;i++)
	for(int j=0;j<DataNy;j++)
	for(int t=0;t<DataNt;t++)
		Recons(i,j,t)/=float(sqrt(DataNt));

	if(Verbose) cerr<<"...End FCur_SFrame::recons"<<endl;
}

// ********************************************************************

void FCur_SFrame::noise_calib(Ifloat **&TabBand, char* Outname)
{
	if(Verbose) cerr<<"FCur_SFrame::noise_calib(.,"<<Outname<<")"<<endl;
	double sigma2;
	fltarray TabStat;
	char filename[64];

	for(int s=0;s<NbrScale;s++)
	{
		int NbrBands=nbr_bands(s);
		for(int b=0;b<NbrBands;b++)
		{
			int nx=TabCur[0].size_band_nc(s,b);
			int ny=TabCur[0].size_band_nl(s,b);
			fltarray Band3d(nx,ny,DataNt);

			if(s==0 & b==0)
			{
				TabStat.alloc(NbrScale, NbrBands);
				TabStat.init(0);
			}

			sigma2=0;
			for(int i=0;i<nx;i++)
				for(int j=0;j<ny;j++)
					sigma2+=pow((TabBand[s][b])(i,j),2);
			if(Verbose) cerr<<"Calib Scale("<<s<<","<<b<<")="<<sqrt(sigma2/(nx*ny))<<endl;
			TabStat(s,b) = sqrt(sigma2/(nx*ny));
		}
		sprintf(filename,"%s_nsig.fits",Outname);
		writefltarr(filename, TabStat);
	}
	if(Verbose) cerr<<"...End FCur_SFrame::noise_calib"<<endl;
}

// ********************************************************************

void FCur_SFrame::extract_stat(Ifloat **&TabBand, char* Outname)
{
	if(Verbose) cerr<<"FCur_SFrame::extract_stat(.,"<<Outname<<")"<<endl;
	double sigma2;
	fltarray TabStat;
	char filename[64];

	for(int s=0;s<NbrScale-1;s++)
	{
		int NbrBands=nbr_bands(s);

		for(int b=0;b<NbrBands;b++)
		{
			int nx=TabCur->size_band_nc(s,b);
			int ny=TabCur->size_band_nl(s,b);

			if(s==0 & b==0)
			{
				TabStat.alloc(NbrScale, NbrBands);
				TabStat.init(0);
			}

				sigma2=0;
				for(int i=0;i<nx;i++)
					for(int j=0;j<ny;j++)
						sigma2+=pow((TabBand[s][b])(i,j)/TabCur->norm_band(s,b),2);
				if(Verbose) cerr<<"Stats Scale("<<s<<","<<b<<")="<<sqrt(sigma2/(nx*ny))<<",coef="<<TabCur->norm_band(s,b)<<endl;
				TabStat(s,b) = sqrt(sigma2/(nx*ny));
		}
		sprintf(filename,"%s_stat.fits",Outname);
		writefltarr(filename, TabStat);
	}
	if(Verbose) cerr<<"...End FCur_SFrame::extract_stat"<<endl;
}

// ********************************************************************

void FCur_SFrame::threshold(Ifloat **&TabBand, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"FCur_SFrame::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<")"<<endl;

	int NS = NbrScale-1 + threshold_coarse;
	for(int s=0;s<NS;s++)
	{
		for(int b=0;b<nbr_bands(s);b++)
		{
			double lvl = NSigma*SigmaNoise;
			for(int i=0;i<TabCur->size_band_nl(s,b);i++)
				for(int j=0;j<TabCur->size_band_nc(s,b);j++)
					(TabBand[s][b])(i,j) *= ( abs((TabBand[s][b])(i,j)) > lvl );
		}
	}
	
	if(Verbose) cerr<<"...End FCur_SFrame::threshold"<<endl;
}



