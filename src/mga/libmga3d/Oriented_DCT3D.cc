/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  27/05/2009
**    
**    File:  Oriented_DCT3D.cc
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  see .h
**    -----------  
**
******************************************************************************/

#include "Oriented_DCT3D.h"

extern Bool Verbose;

Oriented_DCT3D::Oriented_DCT3D()
{
	BlockOverlap=False;
	skip_order=-1;
}

Oriented_DCT3D::~Oriented_DCT3D()
{
	
}

void Oriented_DCT3D::init(int Nx, int Ny, int Nz, int BS, Bool _BlockOverlap)
{
	if(Verbose) cerr<<"Oriented_DCT3D::init("<<Nx<<","<<Ny<<","<<Nz<<","<<BS<<","<<BlockOverlap<<")..."<<endl;
	
// General parameters
	DataNx=Nx;
	DataNy=Ny;
	DataNz=Nz;
	BlockOverlap = _BlockOverlap;

// BlockSize constraints
	if(BS<25) BlockSize = 17;
	else BlockSize = 33;
	NbAngle = 3*BlockSize*BlockSize;

// Block init
	B3D.BlockOverlap = BlockOverlap;
	B3D.alloc(Nx,Ny,Nz,BlockSize);

	if(Verbose) cerr<<"...End Oriented_DCT3D::init()"<<endl;
}

extern void cosft2(float y[], int n, int isign);

void Oriented_DCT3D::transform(fltarray &Cube, fltarray &Trans)
{
	if(Verbose) cerr<<"Oriented_DCT3D::transform(.,.)..."<<endl;
	
	char filename[64];
	Trans.alloc(BlockSize,NbAngle*nbr_block(),1);
	
// Radon transform on each block
	#pragma omp parallel for 
	for(int x=0;x<nbr_block_nx();x++)
	#pragma omp parallel for 
	for(int y=0;y<nbr_block_ny();y++)
	#pragma omp parallel for 
	for(int z=0;z<nbr_block_nz();z++)
	{
		fltarray Block(BlockSize,BlockSize,BlockSize);
		Ifloat Block_Trans;//(NbAngle,BlockSize);
		Radon3D *Rad3D_loc = new Radon3D;
		
		B3D.get_block_cube(x,y,z, Cube, Block);
		Rad3D_loc->transform(Block, Block_Trans);

		for(int j=0;j<NbAngle;j++)
		for(int i=0;i<BlockSize;i++)
			Trans(i,j+NbAngle*num_block(x,y,z),0)=Block_Trans(j,i);
		
		delete Rad3D_loc;
	}

	if(Verbose)
	{
		sprintf(filename,"%s_trans.fits","out");
		writefltarr(filename, Trans);
	}
/*
FFTN_1D F;
fltarray d(BlockSize,1,1);
fltarray dd(BlockSize,1,1);
for(int j=0;j<NbAngle*nbr_block();j++)
{
for(int i=0;i<BlockSize;i++)
	d(i,0,0)=Trans(i,j,0);
F.dct1d(d,dd,False);
for(int i=0;i<BlockSize;i++)
	Trans(i,j,0)=d(i,0,0);//(16./2.);
}
*/


// DCT1D on all radon lines
	double *d = new double[BlockSize];
//	float *d = new float[BlockSize];
	for(int j=0;j<NbAngle*nbr_block();j++)
	{
		for(int i=0;i<BlockSize;i++)
			d[i]=Trans(i,j,0);
//		if(j==0)
//			cerr<<d[0]<<","<<d[1]<<","<<d[2]<<endl;
//		if(j==1)
//			cerr<<d[0]<<","<<d[1]<<","<<d[2]<<endl;
		dct1d(d,false);
//		if(j==0)
//			cerr<<d[0]<<","<<d[1]<<","<<d[2]<<endl;
//		if(j==1)
//			cerr<<d[0]<<","<<d[1]<<","<<d[2]<<endl;
		//cosft2(d-1, BlockSize-1, 1);
		for(int i=0;i<BlockSize;i++)
			Trans(i,j,0)=d[i];
	}

	if(Verbose)
	{
		sprintf(filename,"%s_transT.fits","out");
		writefltarr(filename, Trans);
	}

	if(Verbose) cerr<<"...End Oriented_DCT3D::transform(.,.)"<<endl;
}

void Oriented_DCT3D::recons(fltarray &Trans, fltarray &Recons)
{
	if(Verbose) cerr<<"Oriented_DCT3D::recons(.,.)..."<<endl;
	
	Recons.alloc(DataNx, DataNy, DataNz);
	
	char filename[64];
// iDCT2D on all radon lines
	double *d = new double[BlockSize];
//	float *d = new float[BlockSize];
	for(int j=0;j<NbAngle*nbr_block();j++)
	{
		for(int i=0;i<BlockSize;i++)
			d[i]=Trans(i,j,0);
		dct1d(d,true);
		//cosft2(d-1, BlockSize-1, -1);
		for(int i=0;i<BlockSize;i++)
			Trans(i,j,0)=d[i];//(16./2.);
	}



/*
FFTN_1D F;
fltarray d(BlockSize,1,1);
fltarray dd(BlockSize,1,1);
for(int j=0;j<NbAngle*nbr_block();j++)
{
for(int i=0;i<BlockSize;i++)
	d(i,0,0)=Trans(i,j,0);
F.dct1d(d,dd,True);
for(int i=0;i<BlockSize;i++)
	Trans(i,j,0)=dd(i,0,0);//(16./2.);
}
*/




	if(Verbose)
	{
		sprintf(filename,"%s_transR.fits","out");
		writefltarr(filename, Trans);
	}
	
// Radon inverse transform on each block
	#pragma omp parallel for 
	for(int x=0;x<nbr_block_nx();x++)
	#pragma omp parallel for 
	for(int y=0;y<nbr_block_ny();y++)
	#pragma omp parallel for 
	for(int z=0;z<nbr_block_nz();z++)
	{
		fltarray Block(BlockSize,BlockSize,BlockSize);
		Ifloat Block_Trans(NbAngle,BlockSize);
		Radon3D *Rad3D_loc = new Radon3D;
		for(int j=0;j<NbAngle;j++)
		for(int i=0;i<BlockSize;i++)
			Block_Trans(j,i)=Trans(i,j+NbAngle*num_block(x,y,z),0);
		Rad3D_loc->recons(Block_Trans, Block);
		B3D.add_block_cube(x,y,z, Recons, Block, (Bool)BlockOverlap);
		delete Rad3D_loc;
	}

	if(Verbose)
	{
		sprintf(filename,"%s_recons.fits","out");
		writefltarr(filename, Recons);
	}

	if(Verbose) cerr<<"...End Oriented_DCT3D::recons(.,.)"<<endl;
}

void Oriented_DCT3D::threshold(fltarray &Trans, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"Oriented_DCT3D::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<")..."<<endl;
	
	char filename[64];
	float lvl;
	float cnt=0,cnt2=0;
	
	lvl = SigmaNoise * NSigma * BlockSize; // norm coef : BlockSize
	
// Thresholding
	for(int j=0;j<NbAngle*nbr_block();j++)
	for(int i=1;i<BlockSize;i++)
	{
		cnt+=1;
		if( abs(Trans(i,j,0)) < lvl )
		{// Hard and soft put these to 0
			cnt2+=1;
			Trans(i,j,0)=0;
		}// soft lowers the remaining by lvl
		else if(FilterType==FT_SOFT) Trans(i,j,0) -= (2*float( Trans(i,j,0)>0 )-1)*lvl;
	}
	if(Verbose) cerr<<"n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;

// Coarse scale thresholding
	if(threshold_coarse) 
		for(int j=0;j<NbAngle*nbr_block();j++)
		{
			cnt+=1;
			if( abs(Trans(0,j,0)) < lvl )
			{// Hard and soft put these to 0
				cnt2+=1;
				Trans(0,j,0)=0;
			}// soft lowers the remaining by lvl
			else if(FilterType==FT_SOFT) Trans(0,j,0) -= (2*float( Trans(0,j,0)>0 )-1)*lvl;
		}
		if(Verbose) cerr<<"n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;


// Skipping orders
	if(skip_order > -1)
		for(int j=0;j<NbAngle*nbr_block();j++)
			for(int i=0;i<=skip_order;i++)
				Trans(i,j,0) = 0.F;
			
	if(Verbose)
	{
		sprintf(filename,"%s_transTT.fits","out");
		writefltarr(filename, Trans);
	}
	
	if(Verbose) cerr<<"...End Oriented_DCT3D::threshold()"<<endl;
}

void Oriented_DCT3D::normalize_self(fltarray &TabBand, bool inverse)
{
	if(Verbose) cerr<<"Oriented_DCT3D::normalize_self()..."<<endl;
	// The dct is normalized, so we only apply the Radon3D normalization
	
	float norm = BlockSize;
	float s2 = sqrt(2.);
	if(!inverse)
	{
		for(int j=0;j<NbAngle*nbr_block();j++)
		for(int i=0;i<BlockSize;i++)
			TabBand(i,j,0) /= norm;
		for(int j=0;j<NbAngle*nbr_block();j++)
			TabBand(0,j,0) /= s2;
	}
	else
	{
		for(int j=0;j<NbAngle*nbr_block();j++)
		for(int i=0;i<BlockSize;i++)
			TabBand(i,j,0) *= norm;
		for(int j=0;j<NbAngle*nbr_block();j++)
			TabBand(0,j,0) *= s2;
	}
	
	if(Verbose) cerr<<"end Oriented_DCT3D::normalize_self"<<endl;
}

void Oriented_DCT3D::dct1d(double* in, bool backward)
{
	int i;
	int n=BlockSize;
	double norm = (double)sqrt(double(n)*2.);
	double* out = new double[n];
	fftw_plan plan;

	// Transform
	if(backward)
		plan = fftw_plan_r2r_1d(n, in, out, FFTW_REDFT01, FFTW_ESTIMATE);
	else
		plan = fftw_plan_r2r_1d(n, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
	fftw_execute ( plan );
	
	// Normalization
	for ( i = 0; i < n; i++ )
		in[i] = out[i]/norm;

	fftw_destroy_plan ( plan );
	fftw_free ( out );

	return;
}

void Oriented_DCT3D::write(char *Name, fltarray &Trans, bool Normalize)
{
	writefltarr(Name,Trans);
}

void Oriented_DCT3D::temp(fltarray &Trans)
{
	for(int j=0;j<NbAngle*nbr_block();j++)
	for(int i=0;i<BlockSize;i++)
		Trans(i,j,0) = 0;
	
	for(int i=0;i<BlockSize;i++)
		Trans(i,i,0)=1;	
	
}
