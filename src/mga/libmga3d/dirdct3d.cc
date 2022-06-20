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
**    File:  dirdct3d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  see .h
**    -----------  
**
******************************************************************************/

#include "dirdct3d.h"

extern Bool Verbose;

DirDCT3D::DirDCT3D()
{
	BlockOverlap=False;
}

DirDCT3D::~DirDCT3D()
{
	
}


void DirDCT3D::PRadAlloc(PartialRadon3D &PRad3D)
{
	PRad3D.setTakeOnlyCompletePlane(False);
	PRad3D.setDoNotCorrectSigma(True);   
	PRad3D.setParam(BlockSize,0);  	   // unused second parameter
	PRad3D.alloc();
}

void DirDCT3D::init(int Nx, int Ny, int Nz, int BS, Bool _BlockOverlap)
{
	if(Verbose) cerr<<"DirDCT3D::init("<<Nx<<","<<Ny<<","<<Nz<<","<<BS<<","<<BlockOverlap<<")..."<<endl;
	
// General parameters
	DataNx=Nx;
	DataNy=Ny;
	DataNz=Nz;
	BlockSize = BS;
	BlockOverlap = _BlockOverlap;
	
// BlockSize constraints
	if(BS<18) BlockSize = 17;
	else BlockSize = 17;
	
// Block init
	B3D.BlockOverlap = BlockOverlap;
	B3D.alloc(Nx,Ny,Nz,BlockSize);

// Radon init
	PRadAlloc(PRad3D);
	NbPlane = PRad3D.getNbPlane();
//	cerr<<"NP="<<NbPlane<<endl;
	
	if(Verbose) cerr<<"...End DirDCT3D::init()"<<endl;
}

extern void cosft2(float y[], int n, int isign);


//void DirDCT3D::transform(fltarray &Cube, fltarray &Trans)
void DirDCT3D::transform(fltarray &Cube, fltarray &Trans)
{
	if(Verbose) cerr<<"DirDCT3D::transform(.,.)..."<<endl;
	
	Trans.alloc(BlockSize,BlockSize,NbPlane*nbr_block());
	
// Radon transform on each block
#pragma omp parallel for 
	for(int x=0;x<nbr_block_nx();x++)
#pragma omp parallel for 
	for(int y=0;y<nbr_block_ny();y++)
#pragma omp parallel for 
	for(int z=0;z<nbr_block_nz();z++)
	{
		fltarray Block(BlockSize,BlockSize,BlockSize);
		fltarray Block_Trans(BlockSize,BlockSize,NbPlane);
		PartialRadon3D *PRad3D_loc = new PartialRadon3D;
		PRadAlloc(*PRad3D_loc);
//cerr<<"S";
//cerr<<BlockSize<<":"<<Cube.nz()<<","<<Block.nz();
		B3D.get_block_cube(x,y,z, Cube, Block);
//cerr<<"Block("<<Block.nx()<<","<<Block.ny()<<","<<Block.nz()<<")"<<endl;
//cerr<<"Block_Trans("<<Block_Trans.nx()<<","<<Block_Trans.ny()<<","<<Block_Trans.nz()<<")"<<endl;
//cerr<<"Q";
		PRad3D_loc->transform(Block, Block_Trans);
		
//cerr<<"D";
		for(int k=0;k<NbPlane;k++)
		for(int j=0;j<BlockSize;j++)
		for(int i=0;i<BlockSize;i++)
			Trans(i,j,k+NbPlane*num_block(x,y,z))=Block_Trans(i,j,k);
		
		delete PRad3D_loc;
	}

//cerr<<"X";

// no odd-size dct for now
BlockSize-=1;



char filename[64];
//sprintf(filename,"%s_trans.fits","out");
//writefltarr(filename, Trans);


// DCT2D on all radon planes
	float *d = new float[BlockSize];
	// Transform on the two dimensions
	for(int j=0;j<NbPlane*nbr_block();j++)
	{
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(i,k,j);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Trans(i,k,j)=d[i];
		}
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(k,i,j);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Trans(k,i,j)=d[i];
		}
	}

//sprintf(filename,"%s_transT.fits","out");
//writefltarr(filename, Trans);
	
// no odd-size dct for now
BlockSize++;

	if(Verbose) cerr<<"...End DirDCT3D::transform(.,.)"<<endl;
}

void DirDCT3D::recons(fltarray &Trans, fltarray &Recons)
{
	if(Verbose) cerr<<"DirDCT3D::recons(.,.)..."<<endl;
	
	float *d = new float[BlockSize];
	Recons.alloc(DataNx, DataNy, DataNz);
//cerr<<BlockSize<<","<<Recons.nx()<<","<<Recons.ny()<<","<<Recons.nz()<<endl;


char filename[64];

// no odd-size dct for now
BlockSize--;

// DCT2D on all radon planes
	// Transform on the two dimensions
	for(int j=0;j<NbPlane*nbr_block();j++)
	{
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(i,k,j);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				Trans(i,k,j)=d[i];
		}
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(k,i,j);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				Trans(k,i,j)=d[i];
		}
	}
// Normalisation
	for(int k=0;k<NbPlane*nbr_block();k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
				Trans(i,j,k) /= (BlockSize*BlockSize/4);

// no odd-size dct for now
BlockSize++;

//sprintf(filename,"%s_transR.fits","out");
//writefltarr(filename, Trans);
	
//cerr<<"invDCT2D done"<<endl;

	
// Radon transform on each block
	for(int x=0;x<nbr_block_nx();x++)
	for(int y=0;y<nbr_block_ny();y++)
	for(int z=0;z<nbr_block_nz();z++)
	{
		fltarray Block(BlockSize,BlockSize,BlockSize);
		fltarray Block_Trans(BlockSize,BlockSize,NbPlane);
		PartialRadon3D *PRad3D_loc = new PartialRadon3D;
		PRadAlloc(*PRad3D_loc);
		
		for(int k=0;k<NbPlane;k++)
		for(int j=0;j<BlockSize;j++)
		for(int i=0;i<BlockSize;i++)
			Block_Trans(i,j,k) += Trans(i,j,k+NbPlane*num_block(x,y,z));
		
		PRad3D_loc->recons(Block_Trans, Block);
//cerr<<x<<","<<y<<","<<z<<":"<<Block.nx()<<","<<Block.ny()<<","<<Block.nz()<<":"<<Recons.nx()<<","<<Recons.ny()<<","<<Recons.nz()<<endl;
		B3D.put_block_cube(x,y,z, Recons, Block);
		
		delete PRad3D_loc;
	}
	

//sprintf(filename,"%s_recons.fits","out");
//writefltarr(filename, Recons);

	if(Verbose) cerr<<"...End DirDCT3D::recons(.,.)"<<endl;
}

void DirDCT3D::threshold(fltarray &Trans, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"DirDCT3D::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<")..."<<endl;
	
//	34.6 & 4.1
	float lvl;
	float cnt=0,cnt2=0;
	
	lvl = SigmaNoise * NSigma * 34.6;
	
	for(int k=0;k<NbPlane*nbr_block();k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
			{
						cnt+=1;
						if( abs((Trans)(i,j,k)) < lvl )
						{// Hard and soft put these to 0
							cnt2+=1;
							Trans(i,j,k)=0;
						}// soft lowers the remaining by lvl
						else if(FilterType==FT_SOFT) (Trans)(i,j,k) -= (2*int( (Trans)(i,j,k)>0 )-1)*lvl;
			}
	if(Verbose) cerr<<"n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;
	
//char filename[64];
//sprintf(filename,"%s_transTT.fits","out");
//writefltarr(filename, Trans);
	if(Verbose) cerr<<"...End DirDCT3D::threshold()"<<endl;
}

void DirDCT3D::write(char *Name, fltarray &Trans, bool Normalize)
{
	writefltarr(Name,Trans);
}


















void DirDCT3D::temp(fltarray &Trans)
{
  int i;
  double *in;
  double *in2;
  int n = 17;
  int nc;
  fftw_plan plan_backward;
  fftw_plan plan_forward;

  in = (double*)fftw_malloc ( sizeof ( double ) * n );
  in2 = (double*)fftw_malloc ( sizeof ( double ) * n );
  double* out2 = new double[n];

  for ( i = 0; i < n; i++ )
    in[i] = Trans(0,0,i);

  for ( i = 0; i < n; i++ )
    printf ( "  %4d  %12f\n", i, in[i] );

  plan_forward = fftw_plan_r2r_1d( n, in, out2, FFTW_REDFT11, FFTW_ESTIMATE);
  fftw_execute ( plan_forward );
  for ( i = 0; i < n; i++ )
	printf ( "  %4d  %12f\n", i, out2[i] );


  plan_backward = fftw_plan_r2r_1d( n, out2, in2, FFTW_REDFT11, FFTW_ESTIMATE);

  fftw_execute ( plan_backward );

  for ( i = 0; i < n; i++ )
    printf ( "  %4d  %12f\n", i, in2[i] / ( double ) ( 2.*n ) );

  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );

  fftw_free ( in );
  fftw_free ( in2 );
  fftw_free ( out2 );

  return;
}































