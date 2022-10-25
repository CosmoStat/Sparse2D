
#include "FCur_Frame.h"
extern Bool Verbose;

FCur_Frame::FCur_Frame()
{
	TabCur=NULL;
	
}

FCur_Frame::~FCur_Frame()
{
	
}

void FCur_Frame::alloc_from_coarse(int _NbrScale, int Nx, int Ny, int Nz, int _NbrDir2d, Bool _ExtendWT, Bool _IsotropWT, Bool _RealCur)
{
	DataNx = Nx;
	DataNy = Ny;
	DataNt = Nz;
	NbrScale = _NbrScale;
	NbrDir2d = _NbrDir2d;
	ExtendWT = _ExtendWT;
	RealBand = _RealCur;
	
	TabCur = new FCUR[Nz];
	for(int i=0;i<DataNt;i++)
	{
		TabCur[i].Verbose = Verbose;
		TabCur[i].alloc_from_coarse(NbrScale, Ny, Nx, NbrDir2d, False, False, RealBand);
	}
}

void FCur_Frame::get_norm_coeff(float NSigma)
{
	for(int i=0;i<DataNt;i++)
		TabCur[i].get_norm_coeff(NSigma);
}

void FCur_Frame::transform(fltarray &Data, bool allocTB)
{
	float* Ptr = Data.buffer();
	Ifloat Frame;
	
	for(int i=0;i<DataNt;i++)
	{
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		TabCur[i].cur_trans(Frame);
	}
}

void FCur_Frame::recons(fltarray &Data)
{
	Ifloat Frame;
	Data.resize(DataNx,DataNy,DataNt);
	float* Ptr = Data.buffer();
	
	for(int i=0;i<DataNt;i++)
	{
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		TabCur[i].cur_recons(Frame);
	}
}


void FCur_Frame::threshold(float SigmaNoise, float NSigma, filter_type FilterType)
{
	for(int i=0;i<DataNt;i++)
	{
		if(FilterType==FT_HARD) TabCur[i].threshold(SigmaNoise, NSigma);
		else if(FilterType==FT_WIENER) TabCur[i].wiener_filter(SigmaNoise);
		else cerr<<"Filtering method '"<<FilterType<<"'not implemented yet";
	}
}



