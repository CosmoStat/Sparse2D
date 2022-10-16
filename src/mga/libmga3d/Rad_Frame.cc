
#include "Rad_Frame.h"
extern Bool Verbose;

Rad_Frame::Rad_Frame()
{
	TabRadon=NULL;
	
}

Rad_Frame::~Rad_Frame()
{
	
}

void Rad_Frame::alloc(int nx, int ny, int nz)
{
	if(Verbose) cerr<<"Rad_Frame::alloc()..."<<endl;
	DataNx = nx;
	DataNy = ny;
	DataNz = nz;
	
	TabRadon = new Radon[DataNz+1];
	
	TabRadon[DataNz].alloc(DataNy, DataNx, DEF_RADON);
	RadNl = TabRadon[DataNz].radon_nl();
	RadNc = TabRadon[DataNz].radon_nc();
	
	if(Verbose) cerr<<"RadNlNc="<<RadNl<<","<<RadNc<<endl;
	
	if(Verbose) cerr<<" end Rad_Frame::alloc()"<<endl;
}

void Rad_Frame::transform(fltarray &Data, fltarray &TabBand, bool allocTB)
{
	if(Verbose) cerr<<"Rad_Frame::transform()..."<<endl;
	Ifloat Frame, RFrame;

	if(allocTB)
		TabBand.alloc(RadNc, RadNl, DataNz);
	
	float* Ptr = Data.buffer();
	float* RPtr = TabBand.buffer();
	
	for(int i=0;i<DataNz;i++)
	{
	// Frame pointers allocation
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		RFrame.alloc(RPtr+i*RadNl*RadNc,RadNl,RadNc);
		
	// Radon init
		if(i!=0) TabRadon[i].alloc(DataNy, DataNx, DEF_RADON);
		TabRadon[i].SSR.Verbose=Verbose;
		
	// Forward radon transform	
		TabRadon[i].transform(Frame, RFrame);
	}
	if(Verbose) cerr<<" end Rad_Frame::transform()"<<endl;
}

void Rad_Frame::recons(fltarray &TabBand, fltarray &Data)
{
	if(Verbose) cerr<<"Rad_Frame::recons()..."<<endl;
	Ifloat Frame, RFrame;
	Data.resize(DataNx,DataNy,DataNz);
	float* Ptr = Data.buffer();
	float* RPtr = TabBand.buffer();
	
	
	for(int i=0;i<DataNz;i++)
	{
	// Frame pointers allocation
		Frame.alloc(Ptr+i*DataNy*DataNx,DataNy,DataNx);
		RFrame.alloc(RPtr+i*RadNl*RadNc,RadNl,RadNc);
		
	// Forward radon transform	
		TabRadon[i].recons(RFrame, Frame);
	}
	if(Verbose) cerr<<" end Rad_Frame::recons()"<<endl;
}

void Rad_Frame::threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType)
{
	float lvl = sqrt(RadNl);
	cerr<<lvl;
	for(int k=0;k<DataNz;k++)
		for(int j=0;j<RadNl;j++)
			for(int i=0;i<RadNc;i++)
				TabBand(i,j,k) *= ( TabBand(i,j,k) > NSigma*SigmaNoise*lvl);
}


