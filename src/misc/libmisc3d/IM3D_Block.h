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
**    Date:  28/10/2002 
**    
**    File:  Block3D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition for 3D block manadgement
**    ----------- 
******************************************************************************/

#ifndef _BLOCK3D_H_
#define _BLOCK3D_H_

#include "IM_Obj.h"
  
/***********************************************************************/

class Block3D {   
	private:
		int Nx, Ny, Nz;
		int BlockCubeSize;     // Cube block size  without taking account
                		  // the overlapping aera
		int BlockSize;         // Cube block size 	taking account
                		  // the overlapping aera		  
		int Nxb,Nyb,Nzb;       // Number of blocks in the x,y,z directions
		int NbrBlock;   
		void reset_param();
		
	public:
		Bool Verbose;
		Bool BlockOverlap;
		float Overlap;
		int OverlapSize;
		inline int nx() { return Nx;}
		inline int ny() { return Ny;}
		inline int nz() { return Nz;}
		inline int nbr_block_nx() { return Nxb;}
		inline int nbr_block_ny() { return Nyb;}
		inline int nbr_block_nz() { return Nzb;}
		inline int nbr_block()    { return NbrBlock;}
		inline int block_size()   { return BlockSize;}
		Block3D() {reset_param();} 
		void alloc(int Nxc, int Nyc, int Nzc, int ParamBlockSize);
		void alloc(int Nxc, int Nyc, int Nzc, int ParamBlockSize, float Overlaping);
		void get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False);
		void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock);
		void fold(fltarray &Data);
		void unfold(fltarray &Data);
		float get_weight(int PosPix, int CurrentBlock, int MaxBlockNumber); 
		void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False);
		~Block3D() {} 
};

/***********************************************************************/

#endif
