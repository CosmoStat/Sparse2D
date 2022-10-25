/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  08/09/2010
**    
**    File:  IM3D_Block.cc
**
**    Modification history:
**          A. Woiselle : Entierly rewritten 
**                        with floating overlapping factor between 0 and 0.5
**
*******************************************************************************
**
**    DESCRIPTION 3D block manadgement with partial overlapping
**    -----------  
**                 
******************************************************************************/

#include "IM3D_Block.h"

/***********************************************************************/
void Block3D::reset_param()
{
	Nx=Ny=Nz=0;
	BlockCubeSize=BlockSize=NbrBlock=0;
	Overlap=0;
	Nxb=Nyb=Nzb=0;
	BlockOverlap=False;
	Verbose=False;
}

void Block3D::alloc(int Nxc, int Nyc, int Nzc, int ParamBlockSize, float Overlapping)
{
	Overlap = Overlapping;
	if(Overlap>0)
		BlockOverlap=True;
	alloc(Nxc, Nyc, Nzc, ParamBlockSize);
}

void Block3D::alloc(int Nxc, int Nyc, int Nzc, int ParamBlockSize)
{
	int MaxXYZ;
	Nx=Nxc;
	Ny=Nyc;
	Nz=Nzc;
	BlockSize = ParamBlockSize;
	
	if(Overlap>0)
	{
		BlockOverlap=True;
		if(Overlap>0.5)
			Overlap=0.5;
	}
	else Overlap=0;
	if(BlockOverlap && (Overlap==0))
		Overlap=0.5;

	MaxXYZ = MAX(Nx,Ny);
	MaxXYZ = MAX(MaxXYZ,Nz);

	if (BlockSize==0 || BlockSize==MaxXYZ) 
	{
		Nxb = Nyb = Nzb = 1;
		BlockSize = MaxXYZ;
		BlockCubeSize = BlockSize;
		BlockOverlap = False;
		Overlap = 0;
		OverlapSize = 0;
	}
	else
	{
		BlockCubeSize = ceil(BlockSize*(1-Overlap));
		OverlapSize = BlockSize-BlockCubeSize;

		int ix=Nx-((Nx-OverlapSize)/BlockCubeSize)*BlockCubeSize;
		Nxb = ((Nx-OverlapSize)/BlockCubeSize) + (ix!=OverlapSize);
		int iy=Ny-((Ny-OverlapSize)/BlockCubeSize)*BlockCubeSize;
		Nyb = ((Ny-OverlapSize)/BlockCubeSize) + (iy!=OverlapSize);
		int iz=Nz-((Nz-OverlapSize)/BlockCubeSize)*BlockCubeSize;
		Nzb = ((Nz-OverlapSize)/BlockCubeSize) + (iz!=OverlapSize);
	}
	NbrBlock = Nxb*Nyb*Nzb;
	
	if (Verbose == True)
	{
		printf( "\nCube size = (%2d,%2d,%2d), BlockSize = %2d, BlockCubeSize = %2d, OverlapSize = %2d\n", Nx,Ny,Nz,BlockSize, BlockCubeSize,OverlapSize);
		printf( "Nbr Blocks = (%2d,%2d,%2d)\n", Nxb,Nyb,Nzb);
		printf( "NbrBlocks*BlockCubeSize + OverlapSize = (%2d,%2d,%2d)\n", 
			Nxb*BlockCubeSize+OverlapSize, Nyb*BlockCubeSize+OverlapSize, Nzb*BlockCubeSize+OverlapSize);
	}
}

/***********************************************************************/

float Block3D::get_weight(int PosPix, int CurrentBlock, int NbrBlock) 
{// not normalized when OverlapSize>0.5*BlockSize
	float ValReturn;
	float ValL=1.,ValR=1.;

	if((PosPix < OverlapSize) && (CurrentBlock > 0))
		ValL = (float) pow(sin((double) (PosPix+1) 
				/ (double) OverlapSize * PI/2.), 2.);

	if((PosPix > BlockCubeSize-1) && (CurrentBlock < NbrBlock-1))
		ValR = (float) pow(cos((double) (PosPix-(BlockCubeSize-1)) 
				/ (double) OverlapSize * PI/2.), 2.);
	
	ValReturn = ValL*ValR;
	
	return ValReturn;
}

/*************************************************************************/

void Block3D::get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight)
// Extract the block (Bi,Bj,Bk) from cube and put it in CubeBlock
{
	int k,l,m;
	int Depi,Depj,Depk;
	
	Depi = BlockCubeSize*Bi;
	Depj = BlockCubeSize*Bj;
	Depk = BlockCubeSize*Bk;

	if (Weight == False)
	{
		for (k = 0; k < BlockSize; k++)
		for (l = 0; l < BlockSize; l++)
		for (m = 0; m < BlockSize; m++)
		{
			int Indk = test_index_mirror(Depi+k,Nx);
			int Indl = test_index_mirror(Depj+l,Ny);
			int Indm = test_index_mirror(Depk+m,Nz);
			CubeBlock(k,l,m) = Cube(Indk,Indl,Indm);  
		}
	}
	else
	{
		for (k = 0; k < BlockSize; k++)
		for (l = 0; l < BlockSize; l++)
		for (m = 0; m < BlockSize; m++)
		{
			int Indk = test_index_mirror(Depi+k,Nx);
			int Indl = test_index_mirror(Depj+l,Ny);
			int Indm = test_index_mirror(Depk+m,Nz);
			CubeBlock(k,l,m) = Cube(Indk,Indl,Indm) 
				* get_weight(k, Bi, Nxb) * get_weight(l, Bj, Nyb) * get_weight(m, Bk, Nzb);  
		}
	}
}

/*************************************************************************/

void Block3D::put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
// Extract the block (Bi,Bj,Bk) from cube and put it in CubeBlock
{
	int k,l,m;
	int Depi,Depj,Depk;
	
	Depi = BlockCubeSize*Bi;
	Depj = BlockCubeSize*Bj;
	Depk = BlockCubeSize*Bk;

	for (k = 0; k < BlockSize; k++)
	for (l = 0; l < BlockSize; l++)
	for (m = 0; m < BlockSize; m++)
		if ((Depi+k >= 0) && (Depi+k < Nx) 
			&& (Depj+l >= 0) && (Depj+l < Ny)
			&& (Depk+m >= 0) && (Depk+m < Nz))
			Cube(Depi+k,Depj+l,Depk+m) = CubeBlock(k,l,m);  
}

/*************************************************************************/

void Block3D::add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight)
// Add the block (Bi,Bj,Bk) CubeBlock in Cube  with weighted values
{
	int k,l,m;
	int Depi,Depj,Depk;
	
	Depi = BlockCubeSize*Bi;
	Depj = BlockCubeSize*Bj;
	Depk = BlockCubeSize*Bk;

	if (Weight == False)
	{
		for (k = 0; k < BlockSize; k++)
		for (l = 0; l < BlockSize; l++)
		for (m = 0; m < BlockSize; m++)
		if ((Depi+k >= 0) && (Depi+k < Nx) 
			&& (Depj+l >= 0) && (Depj+l < Ny)
			&& (Depk+m >= 0) && (Depk+m < Nz))
			Cube(Depi+k,Depj+l,Depk+m) += CubeBlock(k,l,m);  
	}
	else
	{
		for (k = 0; k < BlockSize; k++)
		for (l = 0; l < BlockSize; l++)
		for (m = 0; m < BlockSize; m++)
		if ((Depi+k >= 0) && (Depi+k < Nx) 
			&& (Depj+l >= 0) && (Depj+l < Ny)
			&& (Depk+m >= 0) && (Depk+m < Nz))
			Cube(Depi+k,Depj+l,Depk+m) += CubeBlock(k,l,m) * 
				get_weight(k, Bi, Nxb) * get_weight(l, Bj, Nyb) * get_weight(m, Bk, Nzb);  
	}
}

/****************************************************************************/

// ################################################
//             Folding (for DCT)
// ################################################

float fff(int k, float w)
{
	return sin(M_PI/4.*(1.+k/w))*sin(M_PI/4.*(1.+k/w));
}

// takes [-1,1] values, returns decreasing function of x from 1 to 0
float fw(float y)
{
	float x = (y+1)/2.;
	float r,l;
	if(x<=0)
		r = 1;
	else if(x>=1)
		r = 0;
	else
	{
		l = exp(1-1/(1-exp(1-1/(1-x))));
		l *= l;
		r = exp(1-1/(1-exp(1-1/x)));
		r *= r;
		float norm = l+r;
		r /= norm;
	}
	return
//			y<0;
//			r;
//			(1.0-x);
			cos(M_PI/2.*(x))*cos(M_PI/2.*(x));
}

float coef(int k, float w, int BS)
{
//	float dk = float(k)/float(BS);
//	float dw = w/float(BS);
	
	if(k<-w) return 0.;
	else if(k<w) return fw(-float(2*k+1)/(2*float(w)+1));
	else if(k<BS-1-w) return 1.;
	else if(k<BS+w) return fw((2*(float(k)-BS)+1)/(2*float(w)+1)); //(2*(float(k)-(BS-1))-1)/(2*float(w)+1)
	else return 0.;
}

void Block3D::fold(fltarray &Data)
{
	int BS = BlockSize;
	fltarray temp(Nx);

	int i,j,k,it,jt,kt,i0,j0,k0;
	float w=4.;
	
	for (jt = 0; jt < Ny; jt++)
	for (kt = 0; kt < Nz; kt++)
	{
		j = jt%BS;
		j0 = jt-j;
		k = kt%BS;
		k0 = kt-k;
		for (it = 0; it < BS/2; it++)
		{
			i = it%BS;
			i0 = it-i;
			temp(it) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (it = BS/2; it < Nx-BS/2; it++)
		{
			i = it%BS;
			i0 = it-i;
			temp(it) = coef(i,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
					 + coef(2*BS-1-i,w,BS) * Data(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
					 + coef(-1-i,w,BS) * Data(test_index_mirror(i0-i-1,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (it = Nx-BS/2; it < Nx; it++)
		{
			i = it%BS;
			i0 = it-i;
			temp(it) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (it = 0; it < Nx; it++)
			Data(it,jt,kt)=temp(it);
	}
	
	for (it = 0; it < Nx; it++)
	for (kt = 0; kt < Nz; kt++)
	{
		i = it%BS;
		i0 = it-i;
		k = kt%BS;
		k0 = kt-k;
		for (jt = 0; jt < BS/2; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			temp(jt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (jt = BS/2; jt < Ny-BS/2; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			temp(jt) = coef(j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
					 + coef(2*BS-1-j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+k,Nz))
					 + coef(-1-j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0-j-1,Ny), test_index_mirror(k0+k,Nz));
		}
		for (jt = Ny-BS/2; jt < Ny; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			temp(jt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (jt = 0; jt < Nx; jt++)
			Data(it,jt,kt)=temp(jt);
	}

	for (jt = 0; jt < Ny; jt++)
	for (it = 0; it < Nx; it++)
	{
		i = it%BS;
		i0 = it-i;
		j = jt%BS;
		j0 = jt-j;
		for (kt = 0; kt < BS/2; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			temp(kt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (kt = BS/2; kt < Nz-BS/2; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			temp(kt) = coef(k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
					 + coef(2*BS-1-k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
					 + coef(-1-k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0-k-1,Nz));
		}
		for (kt = Nz-BS/2; kt < Nz; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			temp(kt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (kt = 0; kt < Nz; kt++)
			Data(it,jt,kt)=temp(kt);
	}
}

void Block3D::unfold(fltarray &Data)
{
	int BS = BlockSize;
	fltarray temp(Nx);

	int i,j,k,it,jt,kt,i0,j0,k0;
	float w=4.;
	
	for (jt = 0; jt < Ny; jt++)
	for (kt = 0; kt < Nz; kt++)
	{
		j = jt%BS;
		j0 = jt-j;
		k = kt%BS;
		k0 = kt-k;
		for (it = 0; it < BS/2; it++)
		{
			i = it%BS;
			i0 = it-i;
			temp(it) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (it = BS/2; it < Nx-BS/2; it++)
		{
			i = it%BS;
			i0 = it-i;
			if(i>BS/2) 
				temp(it) = (coef(i,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(2*BS-1-i,w,BS) * Data(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz)) )
						 / (coef(i,w,BS)*coef(i,w,BS) - coef(2*BS-1-i,w,BS)*coef(2*BS-1-i,w,BS));
			else 
				temp(it) = (coef(i,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(-1-i,w,BS) * Data(test_index_mirror(i0-i-1,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz)) )
						 / (coef(i,w,BS)*coef(i,w,BS) - coef(-1-i,w,BS)*coef(-1-i,w,BS));
		}
		for (it = Nx-BS/2; it < Nx; it++)
		{
			i = it%BS;
			i0 = it-i;
			temp(it) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (it = 0; it < Nx; it++)
			Data(it,jt,kt)=temp(it);
	}

	for (it = 0; it < Nx; it++)
	for (kt = 0; kt < Nz; kt++)
	{
		i = it%BS;
		i0 = it-i;
		k = kt%BS;
		k0 = kt-k;
		for (jt = 0; jt < BS/2; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			temp(jt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (jt = BS/2; jt < Ny; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			if(j>BS/2) 
				temp(jt) = (coef(j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(2*BS-1-j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+k,Nz)) )
						 / (coef(j,w,BS)*coef(j,w,BS) - coef(2*BS-1-j,w,BS)*coef(2*BS-1-j,w,BS));
			else 
				temp(jt) = (coef(j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(-1-j,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0-j-1,Ny), test_index_mirror(k0+k,Nz)) )
						 / (coef(j,w,BS)*coef(j,w,BS) - coef(-1-j,w,BS)*coef(-1-j,w,BS));
		}
		for (jt = Ny-BS/2; jt < Ny; jt++)
		{
			j = jt%BS;
			j0 = jt-j;
			temp(jt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (jt = 0; jt < Nx; jt++)
			Data(it,jt,kt)=temp(jt);
	}

	for (jt = 0; jt < Ny; jt++)
	for (it = 0; it < Nx; it++)
	{
		i = it%BS;
		i0 = it-i;
		j = jt%BS;
		j0 = jt-j;
		for (kt = 0; kt < BS/2; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			temp(kt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (kt = BS/2; kt < Nz; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			if(k>BS/2) 
				temp(kt) = (coef(k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(2*BS-1-k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+2*BS-1-k,Nz)) )
						 / (coef(k,w,BS)*coef(k,w,BS) - coef(2*BS-1-k,w,BS)*coef(2*BS-1-k,w,BS));
			else 
				temp(kt) = (coef(k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 - coef(-1-k,w,BS) * Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0-k-1,Nz)) )
						/ (coef(k,w,BS)*coef(k,w,BS) - coef(-1-k,w,BS)*coef(-1-k,w,BS));
		}
		for (kt = Nz-BS/2; kt < Nz; kt++)
		{
			k = kt%BS;
			k0 = kt-k;
			temp(kt) = Data(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz));
		}
		for (kt = 0; kt < Nz; kt++)
			Data(it,jt,kt)=temp(kt);
	}
}

/*

// folding in 3D (not separated in 1d foldings)
void Block3D::get_bc(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight)
// Extract the block (Bi,Bj,Bk) from cube and put it in CubeBlock
{
	int BS = BlockSize;
	int i,j,k;
	int i0 = (BS-BlockCubeSize)*Bi;
	int j0 = (BS-BlockCubeSize)*Bj;
	int k0 = (BS-BlockCubeSize)*Bk;

	float w=4.;
	for (i = 0; i < BS; i++)
	for (j = 0; j < BS; j++)
	for (k = 0; k < BS; k++)
	{
		float norm = coef(i,w,BS) * coef(j,w,BS) * coef(k,w,BS) 
				+ coef(2*BS-1-i,w,BS) * coef(j,w,BS) * coef(k,w,BS)
				+ coef(i,w,BS) * coef(2*BS-1-j,w,BS) * coef(k,w,BS)
				+ coef(i,w,BS) * coef(j,w,BS) * coef(2*BS-1-k,w,BS)
				+ coef(-i-1,w,BS) * coef(j,w,BS) * coef(k,w,BS)
				+ coef(i,w,BS) * coef(-j-1,w,BS) * coef(k,w,BS)
				+ coef(i,w,BS) * coef(j,w,BS) * coef(-k-1,w,BS);
		norm=1.;
		CubeBlock(i,j,k) = 
					// center
						coef(i,w,BS) * coef(j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						
					// 6 Faces : x,y,z, -x,-y,-z
						 + coef(2*BS-1-i,w,BS) * coef(j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(i,w,BS) * coef(2*BS-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(i,w,BS) * coef(j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0-i-1,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(i,w,BS) * coef(-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0-j-1,Ny), test_index_mirror(k0+k,Nz))
						 + coef(i,w,BS) * coef(j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0-k-1,Nz))
						
					// 12 Edges : xy,xz,yz, -xy,-xz,-yz, x-y,x-z,y-z, -x-y,-x-z,-y-z
						 + coef(2*BS-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(2*BS-1-i,w,BS) * coef(j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						 + coef(i,w,BS) * coef(2*BS-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(-1-i,w,BS) * coef(j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						 + coef(i,w,BS) * coef(-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						
						 + coef(2*BS-1-i,w,BS) * coef(-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(2*BS-1-i,w,BS) * coef(j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0-1-k,Nz))
						 + coef(i,w,BS) * coef(2*BS-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(-1-j,w,BS) * coef(k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0+k,Nz))
						 + coef(-1-i,w,BS) * coef(j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0+j,Ny), test_index_mirror(k0-1-k,Nz))
						 + coef(i,w,BS) * coef(-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						
					// 8 Corners : xyz, -xyz,x-yz,xy-z, -x-yz,-xy-z,x-y-z, -x-y-z
						 + coef(2*BS-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						 + coef(2*BS-1-i,w,BS) * coef(-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						 + coef(2*BS-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(-1-j,w,BS) * coef(2*BS-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0+2*BS-1-k,Nz))
						 + coef(-1-i,w,BS) * coef(2*BS-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0+2*BS-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						 + coef(2*BS-1-i,w,BS) * coef(-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0+2*BS-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						
						 + coef(-1-i,w,BS) * coef(-1-j,w,BS) * coef(-1-k,w,BS) *
							 Cube(test_index_mirror(i0-1-i,Nx), test_index_mirror(j0-1-j,Ny), test_index_mirror(k0-1-k,Nz))
						
						/norm;
	}
}
*/

