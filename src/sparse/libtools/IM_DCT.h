/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/02/00 
**    
**    File:  IM_DCT.h
**
**    Modification history :
**
******************************************************************************/

#ifndef _DCT_H_
#define _DCT_H_

#include"IM_Obj.h"
#include"IM_Block2D.h"

void ddct16x16s(int isgn, double **a);
void ddct8x8s(int isgn, double **a);

void ddst2d(int n1, int n2, int isgn, double **a, double *t,int *ip, double *w);
// ddst2d: Discrete Sine Transform
void ddct2d(int n1, int n2, int isgn, double **a, double *t,int *ip, double *w);
// ddct2d: Discrete Cosine Transform
void rdft2d(int n1, int n2, int isgn, double **a, double *t,int *ip, double *w);
// rdft2d: Real Discrete Fourier Transform
void cdft2d(int n1, int n2, int isgn, double **a, double *t,int *ip, double *w);
// cdft2d: Complex Discrete Fourier Transform


void im_dct(Ifloat &Ima, Ifloat &Trans, Bool Reverse=False);
void im_dst(Ifloat &Ima, Ifloat &Trans, Bool Reverse=False);
void im_blockdct(Ifloat &Ima, Ifloat &Trans, int BlockSize, Bool Overlap=False, Bool WeightF=False);
void im_blockdct(Ifloat &Ima, Ifloat &Trans, Block2D & B2D, Block2D & B2DTrans);
void im_invblockdct(Ifloat &Trans, Ifloat &Ima, int BlockSize, int Nl, int Nc, Bool Overlap=False, Bool WeightF=False);
void im_invblockdct(Ifloat &Trans, Ifloat &Ima, int BlockSize, Bool Overlap=False, Bool WeightF=False);
void im_invblockdct(Ifloat &Trans, Ifloat &Ima,  Block2D & B2D, Block2D & B2DTrans);

class LOCAL_DCT2D {
      Bool Overlap;
      int BlockSize;
      int Nl, Nc;
      int Nlt, Nct;
      void reset();
      public:
         Block2D B2D, B2DTrans;
	 Ifloat DCTIma;
         LOCAL_DCT2D() {reset();}
        ~LOCAL_DCT2D() {};
	void alloc(int Nl, int Nc, int Block_Size, Bool Overlapping=False, Bool WeightF=False);
        void alloc_from_trans(int Nli, int Nci, int BS, Bool Overlapping=False, Bool WeightF=False);
	int nl_ima() {return Nl;}
	int nc_ima() {return Nc;}
	int nl() {return Nlt;}
	int nc() {return Nct;}
	void set_dct(Ifloat &DCTData);
        inline float & operator() (int i, int j) {return DCTIma(i,j);}
        Ifloat & trans_ima()   { return DCTIma;}
	float * buffer() const { return DCTIma.buffer();}
	int block_size()  {return BlockSize;}
	void transform (Ifloat &Ima);
	void recons (Ifloat &Ima);
	float max() { return DCTIma.max();}
	float max_without_zerofreq();
	float maxfabs_without_zerofreq();
	inline float norm() { return 1.;}
};

#endif
