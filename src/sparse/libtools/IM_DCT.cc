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
**    File:  IM_DCT.cc
**
**     
**
******************************************************************************/


#include"IM_DCT.h"

extern void ddct2d (int, int, int, double **, double *, int *, double *);
extern void ddst2d (int, int, int, double **, double *, int *, double *);
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

void cosft2(float y[], int n, int isign)
{
	void realft(float data[], unsigned long n, int isign);
	int i;
	float sum,sum1,y1,y2,ytemp;
	double theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;

	theta=0.5*PI/n;
	wr1=cos(theta);
	wi1=sin(theta);
	wpr = -2.0*wi1*wi1;
	wpi=sin(2.0*theta);
	if (isign == 1) {
		for (i=1;i<=n/2;i++) {
			y1=0.5*(y[i]+y[n-i+1]);
			y2=wi1*(y[i]-y[n-i+1]);
			y[i]=y1+y2;
			y[n-i+1]=y1-y2;
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
		realft(y,n,1);
		for (i=3;i<=n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr-y[i+1]*wi;
			y2=y[i+1]*wr+y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		sum=0.5*y[2];
		for (i=n;i>=2;i-=2) {
			sum1=sum;
			sum += y[i];
			y[i]=sum1;
		}
	} else if (isign == -1) {
		ytemp=y[n];
		for (i=n;i>=4;i-=2) y[i]=y[i-2]-y[i];
		y[2]=2.0*ytemp;
		for (i=3;i<=n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr+y[i+1]*wi;
			y2=y[i+1]*wr-y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		realft(y,n,-1);
		for (i=1;i<=n/2;i++) {
			y1=y[i]+y[n-i+1];
			y2=(0.5/wi1)*(y[i]-y[n-i+1]);
			y[i]=0.5*(y1+y2);
			y[n-i+1]=0.5*(y1-y2);
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
	}
}

/******************************************************************************/

void im_dct(Ifloat &Ima, Ifloat &Trans, Bool Reverse)
{
    int Nl = Ima.nl();
    int Nc = Ima.nc();
    int i,j;
    float Norm = 2./Nc;
    int isign = (Reverse == False) ? 1: -1;
    float *Ptr = Trans.buffer() - 1;
    fltarray Tab;
    Trans = Ima;
    
    if (Reverse == False)
    {
       for (i = 0; i < Nl; i++) 
       {
          cosft2(Ptr, Nc, isign);
	  Ptr += Nc;
       }
       Tab.alloc(Nl);
       for (j=0; j < Nc; j++)
       {
          for (i = 0; i < Nl; i++) Tab(i) = Trans(i,j);
          Ptr = Tab.buffer() - 1;
          cosft2(Ptr, Nl, isign);
	  for (i = 0; i < Nl; i++) Trans(i,j) = Tab(i);
       }
    }
    else
    {
       Tab.alloc(Nl);
       for (j=0; j < Nc; j++)
       {
          for (i = 0; i < Nl; i++) Tab(i) = Ima(i,j);
          Ptr = Tab.buffer() - 1;
          cosft2(Ptr, Nl, isign);
	  for (i = 0; i < Nl; i++) Trans(i,j) = Tab(i);
       }
       Ptr = Trans.buffer() - 1;
       for (i = 0; i < Nl; i++) 
       {
          cosft2(Ptr, Nc, isign);
	  Ptr += Nc;
       }
    }
    for (i = 0; i < Nl; i++) 
    for (j=0; j < Nc; j++) Trans(i,j) *= Norm; 
}

/******************************************************************************/

void im_dct_old(Ifloat &Ima, Ifloat &Trans, Bool Reverse)
{
     int Nl = Ima.nl();
     int Nc = Ima.nc();
     int Dir = (Reverse == False) ? 1: -1;
     double **a, *t, *w, *dd;
     int i,j;
     // double Norm = 4. / (double)(Nl*Nc);
     double Norm = 1.;
     
     a = new double * [Nl];
     dd = new double [Nl*Nc];
     t = new double [2 * Nl];
     int n = MAX(Nl, Nc/2);
     int *ip = new int [2 + (int) sqrt(n + 0.5)];
     n = MAX(Nl * 5 / 4, Nc * 5 / 4) + Nc / 4;
     w = new double [n];
     ip[0] = 0;
 
     for (i=0; i < Nl; i++) a[i] =  &dd[i*Nc];
     for (i=0; i < Nl; i++)
     for (j=0; j < Nc; j++) a[i][j] = (double) Ima(i,j);
 
      ddct2d (Nl, Nc, Dir, a, t, ip, w);
  
     if (Reverse == True)
     {
         for (i=0; i < Nl; i++)
         for (j=0; j < Nc; j++) Trans(i,j) = (float) (a[i][j]*Norm);
     }
     else
     for (i=0; i < Nl; i++)
     for (j=0; j < Nc; j++) Trans(i,j) = (float) a[i][j];
  
     
     delete [] dd;
     delete [] a;
     delete [] w;
     delete [] t;
     delete [] ip;
 }

/******************************************************************************/

void im_dst(Ifloat &Ima, Ifloat &Trans, Bool Reverse)
{
    int Nl = Ima.nl();
    int Nc = Ima.nc();
    int Dir = (Reverse == False) ? 1: -1;
    double **a, *t, *w, *dd;
    int i,j;
    double Norm = 4. / (double)(Nl*Nc);

    a = new double * [Nl];
    dd = new double [Nl*Nc];
    t = new double [2 * Nl];
    int n = MAX(Nl, Nc/2);
    int *ip = new int [2 + (int) sqrt(n + 0.5)];
    n = MAX(Nl * 5 / 4, Nc * 5 / 4) + Nc / 4;
    w = new double [n];
    ip[0] = 0;

    for (i=0; i < Nl; i++) a[i] =  &dd[i*Nc]; 
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++) a[i][j] = (double) Ima(i,j);
    ddst2d (Nl, Nc, Dir, a, t, ip, w);

    if (Reverse == True)
    {
        for (i=0; i < Nl; i++)
        for (j=0; j < Nc; j++) Trans(i,j) = (float) (a[i][j]*Norm);
    }
    else
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++) Trans(i,j) = (float) a[i][j];


    delete [] dd;
    delete [] a;
    delete [] w;
    delete [] t;
    delete [] ip;
}

/******************************************************************************/

void im_blockdct(Ifloat &Ima, Ifloat &Trans, int BlockSize, Bool Overlap, Bool WeightFirst)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int i,j,Nlt, Nct, BS;
   Block2D B2D,B2DTrans;
   Ifloat ImaBlock, Block, TransBlock;
   // B2D.Verbose = True;
   B2D.BlockOverlap = Overlap;
   B2D.WeightFirst = WeightFirst;
   // B2DTrans.Verbose = True;
   
   B2D.alloc(Nl, Nc, BlockSize);
   BS = B2D.block_size();
   Nlt = B2D.nbr_block_nl() * BS;
   Nct = B2D.nbr_block_nc() * BS;
   
   if ((Trans.nl() != Nlt) || (Trans.nc() != Nct)) Trans.resize(Nlt,Nct);
   ImaBlock.alloc(BS,BS,"ImaBlock");
   TransBlock.alloc(BS,BS,"ImaTransBlock");
   B2DTrans.alloc(Nlt, Nct, BlockSize);
   
   for (i = 0; i < B2D.nbr_block_nl(); i++)
   for (j = 0; j < B2D.nbr_block_nc(); j++)
   {
      if (WeightFirst == False) B2D.get_block_ima(i, j, Ima, ImaBlock);
      else B2D.get_block_ima(i, j, Ima, ImaBlock);
      im_dct(ImaBlock,TransBlock);
      B2DTrans.put_block_ima(i, j, Trans, TransBlock);
   }
}

/******************************************************************************/

void im_blockdct(Ifloat &Ima, Ifloat &Trans, Block2D & B2D, Block2D & B2DTrans)
{
   int i,j,Nlt, Nct, BS;
   Ifloat ImaBlock, Block, TransBlock;
   BS = B2D.block_size();
   Nlt = B2D.l_nl();
   Nct = B2D.l_nc();
   
   if ((Trans.nl() != Nlt) || (Trans.nc() != Nct)) Trans.resize(Nlt,Nct);
   ImaBlock.alloc(BS,BS,"ImaBlock");
   TransBlock.alloc(BS,BS,"ImaTransBlock");   
   for (i = 0; i < B2D.nbr_block_nl(); i++)
   for (j = 0; j < B2D.nbr_block_nc(); j++)
   {
      B2D.get_block_ima(i, j, Ima, ImaBlock);
      im_dct(ImaBlock,TransBlock);
      B2DTrans.put_block_ima(i, j, Trans, TransBlock);
   }
}

/******************************************************************************/
 
void im_invblockdct(Ifloat &Trans, Ifloat &Ima, int BlockSize, Bool Overlap, Bool WeightFirst)
{
   int Nl = Trans.nl();
   int Nc = Trans.nc();
   if (Overlap == True)
   {
      Nl /=2;
      Nc /= 2;
   }
   im_invblockdct(Trans,Ima, BlockSize, Nl, Nc, Overlap, WeightFirst);
}

/******************************************************************************/

void im_invblockdct(Ifloat &Trans, Ifloat &Ima,  Block2D & B2D, Block2D & B2DTrans)
{
   int i,j, Nlt, Nct, BS;
   Ifloat ImaBlock, Block, TransBlock;
   int Nl = B2D.nl();
   int Nc = B2D.nc();
   BS = B2D.block_size();
   Nlt = B2D.nbr_block_nl() * BS;
   Nct = B2D.nbr_block_nc() * BS;
   if ((Trans.nl() != Nlt) || (Trans.nc() != Nct))
   {
       cout << "Error in im_invblockdct: Nlt = " <<  Nlt << " Nct = " <<  Nct  << endl;
       cout << "   Expected values are Nlt = " <<  Trans.nl() << " Nct = " <<  Trans.nc()  << endl;
       exit(-1);
   }
   if ((Ima.nl() != Nl) || (Ima.nc() != Nc)) Ima.resize(Nl,Nc);
   ImaBlock.alloc(BS,BS,"ImaBlock");
   TransBlock.alloc(BS,BS,"ImaTransBlock");
   Ima.init();
   for (i = 0; i < B2D.nbr_block_nl(); i++)
   for (j = 0; j < B2D.nbr_block_nc(); j++)
   {
      B2DTrans.get_block_ima(i, j, Trans, TransBlock);
      im_dct(TransBlock,ImaBlock,True);
      B2D.add_block_ima(i, j, Ima, ImaBlock);
   }
}

/******************************************************************************/

void im_invblockdct(Ifloat &Trans, Ifloat &Ima, int BlockSize, int Nl, int Nc, Bool Overlap, Bool WeightFirst)
{
   int i,j, Nlt, Nct, BS;
   Block2D B2D,B2DTrans;
   Ifloat ImaBlock, Block, TransBlock;
   
   // B2D.Verbose = True;
   B2D.BlockOverlap = Overlap;
   B2D.WeightFirst = WeightFirst;
   B2D.alloc(Nl, Nc, BlockSize);
   BS = B2D.block_size();
   Nlt = B2D.nbr_block_nl() * BS;
   Nct = B2D.nbr_block_nc() * BS;
   if ((Trans.nl() != Nlt) || (Trans.nc() != Nct))
   {
       cout << "Error in im_invblockdct: Nlt = " <<  Nlt << " Nct = " <<  Nct  << endl;
       cout << "   Expected values are Nlt = " <<  Trans.nl() << " Nct = " <<  Trans.nc()  << endl;
       exit(-1);
   }
   if ((Ima.nl() != Nl) || (Ima.nc() != Nc)) Ima.resize(Nl,Nc);
   ImaBlock.alloc(BS,BS,"ImaBlock");
   TransBlock.alloc(BS,BS,"ImaTransBlock");
   B2DTrans.alloc(Nlt, Nct, BlockSize);
   
   Ima.init();
   for (i = 0; i < B2D.nbr_block_nl(); i++)
   for (j = 0; j < B2D.nbr_block_nc(); j++)
   {
      B2DTrans.get_block_ima(i, j, Trans, TransBlock);
      im_dct(TransBlock,ImaBlock,True);
      B2D.add_block_ima(i, j, Ima, ImaBlock);
   }
}

/******************************************************************************/
 
void LOCAL_DCT2D::reset()
{
   Nl = Nc = Nlt =  Nlt = 0;
   BlockSize= 0;
   Overlap= False;
}

/******************************************************************************/

void LOCAL_DCT2D::alloc_from_trans(int Nli, int Nci, int BS, Bool Overlapping, Bool WeightF)
{
    int Nl1 = Nli;
    int Nc1 = Nci;
    if (Overlapping == True)
    {
       Nl1 /= 2;
       Nc1 /= 2;
    }
    alloc(Nl1, Nc1, BS, Overlapping, WeightF);
}

/******************************************************************************/

void LOCAL_DCT2D::alloc(int Nli, int Nci, int BS, Bool Overlapping, Bool WeightF)
{
   Nl= Nli;
   Nc = Nci;
   BlockSize = BS;
   Overlap = Overlapping;
   
   if ((Nl == 0) || (Nc == 0))
   {
      cout << "Error: LOCAL_DCT2D class is allocated for an image of zero dimension ... " << endl;
      cout << "       Image size is: Nl = " << Nl << " Nc = " << Nc <<  endl;
      exit(-1);
   }
   
   
   if ((BlockSize > MIN(Nl,Nc)) || (BlockSize <= 0)) BlockSize = MIN(Nl,Nc);
   if (is_power_of_2(BlockSize) == False)
   {
      cout << "Error: block size must be power of two ... " << endl;
      cout << "       BlockSize = " << BlockSize << endl;
      exit(-1);
   }
   B2D.BlockOverlap = Overlapping;
   B2D.WeightFirst = WeightF;
   B2D.alloc(Nl, Nc, BlockSize);
   BlockSize = B2D.block_size();
   Nlt = B2D.nbr_block_nl() * BlockSize;
   Nct = B2D.nbr_block_nc() * BlockSize;
   B2DTrans.alloc(Nlt, Nct, BlockSize);
   DCTIma.alloc(Nlt, Nct, "DCTIma");
}
	
/******************************************************************************/

void LOCAL_DCT2D::transform(Ifloat &Ima)
{
   if ((Nl == 0) || (Nc == 0))
   {
      cout << "Error: LOCAL_DCT2D class is not allocated ... " << endl;
      exit(-1);
   }
   
   if ((Ima.nl() != Nl) || (Ima.nc() != Nc)) 
   {
      cout << "Error: Bad image size in LOCAL_DCT2D::transform ... " << endl;
      cout << "       Image size is: Nl = " << Ima.nl() << " Nc = " << Ima.nc() <<  endl;
      cout << "       Expected size is: Nl = " << Nl << " Nc = " << Nc <<  endl;
      exit(-1);
   }
   im_blockdct(Ima, DCTIma, B2D, B2DTrans);
}

/******************************************************************************/

void LOCAL_DCT2D::recons(Ifloat &Ima)
{
   if ((Nl == 0) || (Nc == 0))
   {
      cout << "Error: LOCAL_DCT2D class is not allocated ... " << endl;
      exit(-1);
   }
   
   if ((Ima.nl() != Nl) || (Ima.nc() != Nc))  Ima.resize(Nl,Nc);
   im_invblockdct(DCTIma, Ima, B2D, B2DTrans);
}

/******************************************************************************/

void LOCAL_DCT2D::set_dct(Ifloat &DCTData)
{
   if ((DCTData.nl() != Nlt) || (DCTData.nc() != Nct)) 
   {
      cout << "Error: Bad image size in LOCAL_DCT2D::set_dct ... " << endl;
      cout << "       Image size is: Nl = " << DCTData.nl() << " Nc = " << DCTData.nc() <<  endl;
      cout << "       Expected size is: Nl = " << Nlt << " Nc = " << Nct <<  endl;
      exit(-1);
   }
   DCTIma = DCTData;
}

/******************************************************************************/

float LOCAL_DCT2D::max_without_zerofreq()
{
   int bi,bj;
   int BS = B2DTrans.block_size();
   Ifloat BlockCosIma(BS, BS, "blockima");
   float Max=0;
   float MaxBlock;
   for (bi = 0; bi < B2DTrans.nbr_block_nl(); bi++)
   for (bj = 0; bj < B2DTrans.nbr_block_nc(); bj++)
   {
       B2DTrans.get_block_ima(bi,bj, DCTIma, BlockCosIma);
       BlockCosIma(0,0) = 0.;
       MaxBlock = BlockCosIma.max();
       if (Max < MaxBlock) Max = MaxBlock;
   }      
   return Max;
}

/******************************************************************************/

float LOCAL_DCT2D::maxfabs_without_zerofreq()
{
   int bi,bj;
   int BS = B2DTrans.block_size();
   Ifloat BlockCosIma(BS, BS, "blockima");
   float Max=0;
   float MaxBlock;
   for (bi = 0; bi < B2DTrans.nbr_block_nl(); bi++)
   for (bj = 0; bj < B2DTrans.nbr_block_nc(); bj++)
   {
       B2DTrans.get_block_ima(bi,bj, DCTIma, BlockCosIma);
       BlockCosIma(0,0) = 0.;
       MaxBlock = BlockCosIma.maxfabs();
       if (ABS(Max) < ABS(MaxBlock)) Max = MaxBlock;
   }      
   return Max;
}
