/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  11/06/03 
**    
**    File:  IM1D_Dct.h
**
**    Modification history :
**
******************************************************************************/

#ifndef _1DDCT_H_
#define _1DDCT_H_

//#include"IM_Obj.h"
#include "IM1D_Block.h"

void ddct16x16s (int isgn, double **a);
void ddct8x8s (int isgn, double **a);

void ddst1d(int n1, int isgn, double **a, double *t,int *ip, double *w);
// ddst1d: Discrete Sine Transform
void ddct1d(int n1, int isgn, double **a, double *t,int *ip, double *w);
// ddct1d: Discrete Cosine Transform
void rdft1d(int n1, int isgn, double **a, double *t,int *ip, double *w);
// rdft1d: Real Discrete Fourier Transform
void cdft1d(int n1, int isgn, double **a, double *t,int *ip, double *w);
// cdft1d: Complex Discrete Fourier Transform


void four1(float data[], unsigned long nn, int isign);
void realft(float data[], unsigned long n, int isign);
void cosft2(float y[], int n, int isign);


void sig_dct (fltarray& Sig, fltarray& Trans, Bool Reverse=False);
void sig_blockdct (fltarray& Sig, fltarray& Trans, int BlockSize, 
                   Bool Overlap, Bool WeightF=False);
void sig_blockdct (fltarray& Sig, fltarray& Trans, Block1D& B1D, 
                   Block1D & B1DTrans);
void sig_invblockdct (fltarray& Trans, fltarray& Sig, 
                      int BlockSize, int Nx, Bool Overlap, 
                      Bool WeightF=False);
void sig_invblockdct (fltarray& Trans, fltarray& Sig, int BlockSize, 
                      Bool Overlap, Bool WeightF=False);
void sig_invblockdct (fltarray& Trans, fltarray& Sig,  Block1D & B1D, 
                      Block1D & B1DTrans);

class LOCAL_DCT1D {
private:
   Bool _Overlap;
   int  _BlockSize;
   int  _Nx;
   int  _Nxt;
   void reset();
public:
   Block1D _B1D, _B1DTrans;
   fltarray _DCTSig;
   
   Bool SuppressDctCompCont;
   float COSSupportMin;
   float SigmaNoise;
   bool Verbose;
   bool Write;
   bool UseNormL1;
   float NSigmaSoft;
   float COS_Sensibility;

   LOCAL_DCT1D() {reset();}
   ~LOCAL_DCT1D() {};
   void alloc (int Nx, int Block_Size, Bool Overlapping=False, 
               Bool WeightF=False);
   void alloc_from_trans (int Nxi, int BS, Bool OverlappingFalse, 
                          Bool WeightF=False);
                          
   int nx_Sig() {return  _Nx;}
   int nx() {return  _Nxt;}
   void set_dct (fltarray& DCTData);
   inline float & operator() (int i) {return _DCTSig(i);}
   fltarray&  trans_Sig() {return _DCTSig;}
   float * buffer() const { return _DCTSig.buffer();}
   int block_size()  {return _BlockSize;}
   void threshold (float Lamba_Sigma);
   float update (float CoefSol, float Threshold, float SoftLevel);
   void transform (fltarray& Sig);
   void recons (fltarray& Sig);
   float getAbsMaxTransf (fltarray& Sig);
   inline float norm() { return 1.;}
};

#endif
