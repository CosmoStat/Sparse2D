
#include "FCur_TFrame.h"
extern Bool Verbose;

FCur_TFrame::FCur_TFrame()
{
	TabCur=NULL;
	kill_time_fine=0;
}

FCur_TFrame::~FCur_TFrame()
{
	
}

// ********************************************************************

void FCur_TFrame::alloc_from_coarse(int _NbrScale, int _NbrScale1D, int Nx, int Ny, int Nz, int _NbrDir2d, Bool _ExtendWT, Bool _IsotropWT, Bool _RealCur)
{
	DataNx = Nx;
	DataNy = Ny;
	DataNt = Nz;
	NbrScale = _NbrScale;
	NbrDir2d = _NbrDir2d;
	ExtendWT = _ExtendWT;
	RealBand = _RealCur;

	if(_NbrScale1D != 0) NbrScale1D=min(_NbrScale1D,(int)floor(log2(DataNt)));
	else NbrScale1D=(int)floor(log2(DataNt));

	TabCur = new FCUR[DataNt];
	for(int i=0;i<DataNt;i++)
	{
		TabCur[i].Verbose = False;
		TabCur[i].alloc_from_coarse(NbrScale, Ny, Nx, NbrDir2d, False, False, RealBand);
	}

// CubeSigma initialisation
	char filename[128];
	sprintf(filename,"%s/common/F2DTC_TabSigma/TabSigma_n%02d_d%d.fits",getenv("CUR"),NbrScale,NbrDir2d);
	if(access(filename, F_OK)==0)
		fits_read_fltarr(filename, TabSigma);
	else
	{
		cerr<<" Error : cannot load normalisation file "<<endl<<"  "<<filename<<endl;
		exit(-7);
	}
}

// ********************************************************************

void FCur_TFrame::get_norm_coeff(float NSigma)
{
//	for(int i=0;i<DataNt;i++)
//		TabCur[i].get_norm_coeff(NSigma);
}

// ********************************************************************

void FCur_TFrame::transform(fltarray &Data, fltarray ***&TabBand, bool allocTB)
{
	if(Verbose) cerr<<"FCur_TFrame::transform(.,.,"<<allocTB<<")"<<endl;
	Bool LVerbose=Verbose;
	
	float* Ptr = Data.buffer();
	Ifloat Frame;
	fltarray time(DataNt);
	MEYER_WT1D Wave;
	fltarray *timeW;
	Bool alloc_tW=True;
	
	for(int i=0;i<DataNt;i++)
	{
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		TabCur[i].cur_trans(Frame);
	}
	
	Verbose=False;
	Wave.init(NbrScale1D, DataNt);
	Verbose=LVerbose;
	TabBand = new fltarray **[NbrScale];

	for(int s=0;s<NbrScale;s++)
	{
		int NbrBands=nbr_bands(s);
		TabBand[s] = new fltarray *[NbrBands];

		for(int b=0;b<NbrBands;b++)
		{
			TabBand[s][b] = new fltarray[NbrScale1D];

			int nx=TabCur[0].size_band_nc(s,b);
			int ny=TabCur[0].size_band_nl(s,b);
			fltarray Band3d(nx,ny,DataNt);
			float* ptr = Band3d.buffer();
			Ifloat Band2d;
			for(int i=0;i<DataNt;i++)
			{
				Band2d.alloc(ptr+i*nx*ny, ny, nx);
				TabCur[i].get_band(s, b, Band2d);
			}

			for(int i=0;i<nx;i++)
				for(int j=0;j<ny;j++)
				{
					for(int t=0;t<DataNt;t++)
						time(t) = Band3d(i,j,t);
					Verbose=False;
					Wave.transform(time, timeW, alloc_tW);
					Verbose=LVerbose;
					alloc_tW=False;
					for(int s1=0;s1<NbrScale1D;s1++)
					{
						if(i==0 && j==0) 
							TabBand[s][b][s1].alloc(nx,ny,timeW[s1].nx());
						for(int t=0;t<timeW[s1].nx();t++)
							(TabBand[s][b][s1])(i,j,t)=(timeW[s1])(t);
					}
				}
		}
	}

	if(kill_time_fine)
		for(int s=0;s<NbrScale;s++)
			for(int b=0;b<nbr_bands(s);b++)
				for(int s1=0;s1<kill_time_fine;s1++)
					TabBand[s][b][s1].init(0.0);

	if(Verbose) cerr<<"...End FCur_TFrame::transform"<<endl;
}

// ********************************************************************

void FCur_TFrame::recons(fltarray ***&TabBand, fltarray &Recons)
{
	if(Verbose) cerr<<"FCur_TFrame::recons(.,.)"<<endl;
	Bool LVerbose=Verbose;
	
	fltarray time(DataNt);
	MEYER_WT1D Wave;
	Ifloat Frame;
	float* Ptr = Recons.buffer();
	fltarray *timeW;
	timeW = new fltarray[NbrScale1D];
	for(int s1=0;s1<NbrScale1D;s1++)
		timeW[s1].alloc(DataNt/int(pow((float)2.,(float)s1)));
	
	Verbose=False;
	Wave.init(NbrScale1D, DataNt);
	Verbose=LVerbose;
	
	for(int s=0;s<NbrScale;s++)
	{
		int NbrBands=nbr_bands(s);

		for(int b=0;b<NbrBands;b++)
		{
			int nx=TabCur[0].size_band_nc(s,b);
			int ny=TabCur[0].size_band_nl(s,b);
			fltarray Band3d(nx,ny,DataNt);

			for(int i=0;i<nx;i++)
				for(int j=0;j<ny;j++)
				{
					for(int s1=0;s1<NbrScale1D;s1++)
					{
						for(int t=0;t<timeW[s1].nx();t++)
							(timeW[s1])(t) = (TabBand[s][b][s1])(i,j,t);
					}
					Verbose=False;
					Wave.recons(timeW, time, False);
					Verbose=LVerbose;
					for(int t=0;t<DataNt;t++)
						Band3d(i,j,t) = time(t);
				}
			float* ptr = Band3d.buffer();
			Ifloat Band2d;
			for(int i=0;i<DataNt;i++)
			{
				Band2d.alloc(ptr+i*nx*ny, ny, nx);
				TabCur[i].put_band(s, b, Band2d);
			}
		}
	}
	
	for(int i=0;i<DataNt;i++)
	{
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		TabCur[i].cur_recons(Frame);
	}
	if(Verbose) cerr<<"...End FCur_TFrame::recons"<<endl;
}

// ********************************************************************

void FCur_TFrame::noise_calib(fltarray ***&TabBand, char* Outname)
{
	if(Verbose) cerr<<"FCur_TFrame::noise_calib(.,"<<Outname<<")"<<endl;
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
				TabStat.alloc(NbrScale, NbrBands, NbrScale1D-1);
				TabStat.init(0);
			}

			for(int s1=0;s1<NbrScale1D-1;s1++)
			{
				int nt=(TabBand[s][b][s1]).nz();
				sigma2=0;
				for(int i=0;i<nx;i++)
					for(int j=0;j<ny;j++)
						for(int t=0;t<nt;t++)
							sigma2+=pow((TabBand[s][b][s1])(i,j,t),2);
				if(Verbose) cerr<<"Calib Scale("<<s<<","<<b<<","<<s1<<")="<<sqrt(sigma2/(nt*nx*ny))<<endl;
				TabStat(s,b,s1) = sqrt(sigma2/(nt*nx*ny));
			}
		}
		sprintf(filename,"%s_nsig.fits",Outname);
		writefltarr(filename, TabStat);
	}
	if(Verbose) cerr<<"...End FCur_TFrame::noise_calib"<<endl;
}	

// ********************************************************************

void FCur_TFrame::extract_stat(fltarray ***&TabBand, char* Outname)
{
	if(Verbose) cerr<<"FCur_TFrame::extract_stat(.,"<<Outname<<")"<<endl;
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
				TabStat.alloc(NbrScale, NbrBands, NbrScale1D-1);
				TabStat.init(0);
			}

			for(int s1=0;s1<NbrScale1D-1;s1++)
			{
				int nt=(TabBand[s][b][s1]).nz();
				sigma2=0;
				for(int i=0;i<nx;i++)
					for(int j=0;j<ny;j++)
						for(int t=0;t<nt;t++)
							sigma2+=pow((TabBand[s][b][s1])(i,j,t)/TabSigma(s,b,s1),2);
				if(Verbose) cerr<<"Stats Scale("<<s<<","<<b<<","<<s1<<")="<<sqrt(sigma2/(nt*nx*ny))<<",coef="<<TabSigma(s,b,s1)<<endl;
				TabStat(s,b,s1) = sqrt(sigma2/(nt*nx*ny));
			}
		}
		sprintf(filename,"%s_stat.fits",Outname);
		writefltarr(filename, TabStat);
	}
	if(Verbose) cerr<<"...End FCur_TFrame::extract_stat"<<endl;
}	

// ********************************************************************

void FCur_TFrame::threshold(fltarray ***&TabBand, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"FCur_TFrame::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<")"<<endl;
	for(int s=0;s<NbrScale;s++)
	{
		for(int b=0;b<nbr_bands(s);b++)
		{
			int nx=TabCur[0].size_band_nc(s,b);
			int ny=TabCur[0].size_band_nl(s,b);
			fltarray Band3d(nx,ny,DataNt);
			
			for(int s1=0;s1<NbrScale1D;s1++)
			{
				if(s!=NbrScale-1 || s1!=NbrScale1D-1)
				{
					int nt=(TabBand[s][b][s1]).nz();
					double lvl = NSigma*SigmaNoise*TabSigma(s,b,s1);
					for(int i=0;i<nx;i++)
						for(int j=0;j<ny;j++)
							for(int t=0;t<nt;t++)
								(TabBand[s][b][s1])(i,j,t) *= ( abs((TabBand[s][b][s1])(i,j,t)) > lvl );
				}
			}
		}
	}
	
	if(Verbose) cerr<<"...End FCur_TFrame::threshold"<<endl;
}



