
#include "DCT_Frame.h"
extern Bool Verbose;

DCT_Frame::DCT_Frame()
{
	
}

DCT_Frame::~DCT_Frame()
{
	
}

void DCT_Frame::alloc(int nx, int ny, int nz, int _BlockSize, bool _BlockOverlap)
{
	if(Verbose) cerr<<"DCT_Frame::alloc()..."<<endl;
	DataNx = nx;
	DataNy = ny;
	DataNz = nz;
	BlockSize = _BlockSize;
	BlockOverlap = _BlockOverlap;
	skip_order = -1;
	
	B2D.BlockOverlap = (Bool) BlockOverlap;
	B2D.WeightFirst = False;
	B2D.alloc(DataNy, DataNx, BlockSize);
	
	TransNx = B2D.nbr_block_nc()*BlockSize;
	TransNy = B2D.nbr_block_nl()*BlockSize;
	
	B2DTrans.BlockOverlap = False;
	B2DTrans.WeightFirst = False;
	B2DTrans.alloc(TransNy, TransNx, BlockSize);
	
	if(Verbose) cerr<<" end DCT_Frame::alloc()"<<endl;
}

void dctw2d(Ifloat &In, Ifloat &Out, int nx, int ny, bool backward)
{
//	if(Verbose) cerr<<"dctw2d..."<<endl;
	double norm;
	double* in = new double[nx];
	double* out = new double[nx];
	fftw_plan plan;

// dct on x axis
	norm = (double)sqrt(double(nx)*2.);
	for(int j=0;j<ny;j++)
	{
		for (int i = 0; i < nx; i++ )
			in[i] = In(j,i);

		// Transform
		if(backward)
			plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT01, FFTW_ESTIMATE);
		else
			plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
		fftw_execute ( plan );

		// Normalization
		for (int i = 0; i < nx; i++ )
			Out(j,i) = out[i]/norm;

	}

// dct on y axis
	norm = (double)sqrt(double(ny)*2.);
	for(int j=0;j<ny;j++)
	{
		for (int i = 0; i < nx; i++ )
			in[i] = Out(i,j);

		// Transform
		if(backward)
			plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT01, FFTW_ESTIMATE);
		else
			plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
		fftw_execute ( plan );

		// Normalization
		for (int i = 0; i < nx; i++ )
			Out(i,j) = out[i]/norm;

	}
	
	fftw_destroy_plan ( plan );
	fftw_free ( out );

//	if(Verbose) cerr<<" end dctw2d"<<endl;
}

void DCT_Frame::transform(fltarray &Data, fltarray &TabBand, bool allocTB)
{
	if(Verbose) cerr<<"2D_DCT_Frame::transform()..."<<endl;
	Ifloat Frame, TFrame;

	if(allocTB)
		TabBand.alloc(TransNx, TransNy, DataNz);
	TabBand.init(0.F);
	
	float* Ptr = Data.buffer();
	float* RPtr = TabBand.buffer();
	
	Ifloat ImaBlock(BlockSize, BlockSize);
	Ifloat TransBlock(BlockSize, BlockSize);
	
	for(int k=0;k<DataNz;k++)
	{
	// Frame pointers allocation
		Frame.alloc(Ptr+k*DataNy*DataNx,DataNx,DataNy);
		TFrame.alloc(RPtr+k*TransNy*TransNx,TransNx,TransNy);
		
		for (int Bi = 0; Bi < B2D.nbr_block_nc(); Bi++)
		for (int Bj = 0; Bj < B2D.nbr_block_nl(); Bj++)
		{
			B2D.get_block_ima(Bj, Bi, Frame, ImaBlock);
			
		// Forward DCT transform
//			cerr<<ImaBlock.nl()<<","<<ImaBlock.nc()<<","<<TransBlock.nl()<<","<<TransBlock.nc()<<","<<endl;
			dctw2d(ImaBlock, TransBlock, BlockSize, BlockSize, false);
			
			B2DTrans.put_block_ima(Bj, Bi, TFrame, TransBlock);
		}
	}
	if(Verbose) cerr<<" end 2D_DCT_Frame::transform()"<<endl;
}

void DCT_Frame::recons(fltarray &TabBand, fltarray &Data)
{
	if(Verbose) cerr<<"DCT_Frame::recons()..."<<endl;
	Ifloat Frame, TFrame;
	
	Data.resize(DataNx,DataNy,DataNz);
	Data.init(0.F);
	
	float* Ptr = Data.buffer();
	float* RPtr = TabBand.buffer();
	
	Ifloat ImaBlock(BlockSize, BlockSize);
	Ifloat TransBlock(BlockSize, BlockSize);

	for(int k=0;k<DataNz;k++)
	{
	// Frame pointers allocation
		Frame.alloc(Ptr+k*DataNy*DataNx,DataNy,DataNx);
		TFrame.alloc(RPtr+k*TransNy*TransNx,TransNy,TransNx);
		
		for (int Bi = 0; Bi < B2D.nbr_block_nc(); Bi++)
		for (int Bj = 0; Bj < B2D.nbr_block_nl(); Bj++)
		{
			B2DTrans.get_block_ima(Bj, Bi, TFrame, TransBlock);
ImaBlock.init();
		// Backward DCT transform
			dctw2d(TransBlock, ImaBlock, BlockSize, BlockSize, true);
			
			B2D.add_block_ima(Bj, Bi, Frame, ImaBlock);
		}
	}
	if(Verbose) cerr<<" end DCT_Frame::recons()"<<endl;
}

void DCT_Frame::threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"DCT_Frame::threshold(.,"<<SigmaNoise<<","<<NSigma<<")..."<<endl;
	for(int k=0;k<DataNz;k++)
		for(int j=0;j<TransNy;j++)
			for(int i=0;i<TransNx;i++)
				TabBand(i,j,k) *= ( abs(TabBand(i,j,k)) > NSigma*SigmaNoise);
	if(skip_order>-1)
	{
		skip_order = min(skip_order,BlockSize-1);
		
		for(int k=0;k<DataNz;k++)
			for (int Bi = 0; Bi < B2D.nbr_block_nl(); Bi++)
			for (int Bj = 0; Bj < B2D.nbr_block_nc(); Bj++)
				for(int x=0;x<=skip_order;x++)
				for(int y=0;y<=skip_order-x;y++)
					TabBand(Bi*BlockSize+y, Bj*BlockSize+x, k)=0;
	}
	
	if(Verbose) cerr<<" end DCT_Frame::threshold()"<<endl;
}






































