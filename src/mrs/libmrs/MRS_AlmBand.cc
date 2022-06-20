/******************************************************************************
**                   Copyright (C) 2009 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  25/10/08
**    
**    File:  MRS_AlmBand.cc
**
*******************************************************************************
**
**    DESCRIPTION: Creation of band filter in the spherical harmonic domain  
**    ------------
**
******************************************************************************/


#include"MRS_Sparse.h"


// Band for the planck power spectrum binning
#define PS_NB 8
int  PS_From[PS_NB] = {0, 2,  11,  31, 151, 421,  1201, 2501};
int  PS_To[PS_NB] = {1, 10, 30, 150, 420, 1200, 2500, 3000};
int  PS_Step[PS_NB] = {2, 1,   2,   5,  10,  20,  50,  100};

// Band to limit the number of bans, but still whitening the CMB
int TNbrWin[9] = {22,22,11,8,17,11,5,3,4};
double TabWin2048[22] ={64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560,2688,2816,3000}; 
double TabWin1024[22] ={64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560,2688,2816,3000}; 
double TabWin512[11] = {64,128,192,256,320,384,512,768,1024,1280,1536};
double TabWin256[8] = {64,128,192,256,320,384,512,768};
double TabWin128[17] = {64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384};
double TabWin64[11] = {16,32,64,80,96, 112,128,144,160, 176, 192};
double TabWin32[5] = {8,16,32,64,92};
double TabWin16[3] = {8,16,48};
double TabWin8[4] = {4,8,16,24};


/*********************************************************************/

static double meyer_lowpass_fct(double x)
{
    double l,r;
    l = exp( 1-1 / (1-exp(1-1/(1-x))));
    r = exp(1-1/(1-exp(1-1/x)));
    return (l /= sqrt(l*l+r*r));
}

/***************************************************************************/

static void make_meyer_filter_lowpass(fltarray &H,  int PixStart2Zero, int PixEnd2Zero)
{
	int i;
	if (PixEnd2Zero <= PixStart2Zero) PixEnd2Zero = PixStart2Zero + 1;
	int Npix2zero = PixEnd2Zero - PixStart2Zero + 1;
    H.init();
    for (i=0; i < PixStart2Zero; i++)  H(i) = 1.;    
 	for (i=0; i < MIN(H.nx(),PixEnd2Zero); i++)
	{
        double x = i;
        double r =  x  / ( (double) Npix2zero -1);
        if ((PixStart2Zero+i) < H.nx())
        {
		    if (r <=0) H(PixStart2Zero+i) = 1.;
 		    else if (r >= 1)  H(PixStart2Zero+i) = 0.;
		    else H(PixStart2Zero+i) = meyer_lowpass_fct(1.-r);
        }
        // cout <<  i+1 << " nx = " << H.nx() << " r= " << r  <<  "  H = " << H(i) << endl;
 	}
 	for (i=PixEnd2Zero; i < H.nx(); i++)  H(i) = 0;
}

/***************************************************************************/

void get_planck_wp_meyer_filter(fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int Lmax, int LmaxT)
{
    int Nw=0;
    for (int s=0; s < PS_NB; s++)
        for (int f=PS_From[s] ; f < PS_To[s] ; f+=PS_Step[s])
        {
            if (f <= Lmax) Nw++;
        }
    //  cout << "Nbr Win = " << Nw << endl;
    TabH.resize(LmaxT+1, Nw);
    Win.resize(Nw,2);
    fltarray H(TabH.nx());
    
    int w=0;
    for (int s=0; s < PS_NB; s++)
        for (int f=PS_From[s] ; f < PS_To[s] ; f+=PS_Step[s])
        {
            int Beg = MAX(0,f - PS_Step[s]/2);
            int End = f + PS_Step[s]/2;
            if (f <= Lmax)
            {
                //  cout << w+1 << " " << Beg << "--> " << End << endl;
                make_meyer_filter_lowpass(H, f,  End);
                for (int l=0; l < H.nx(); l++)  TabH(l,w) = H(l);
                Win(w,0) = Beg;
                Win(w,1) = End;
                w++;
            }
        }
    
    WP_WPFilter.alloc(LmaxT+1, Nw+1);
    for (int b=0; b <= Nw; b++)
    {
        if (b == Nw)
            for (int l=0; l < H.nx(); l++)  WP_WPFilter(l,Nw) = TabH(l,0);
        else if (b == 0) for (int l=0; l < H.nx(); l++)  WP_WPFilter(l,0) = 1. - TabH(l,Nw-1);
        else  for (int l=0; l < H.nx(); l++)  WP_WPFilter(l,b) = TabH(l,Nw-b) - TabH(l,Nw-b-1);
    }
    
    //  fits_write_fltarr("xxz.fits", WP_WPFilter);
}


// Compute the filters for the Wavelet Packet decompositions.


void hgmey(int npix1,fltarray & H)
{
    float scale=0.5;
    int npix = 2 * npix1/2 + 1;
    H.resize(npix);
    int n1 = (int) (npix / POW(2.,scale+1) + 0.5);
    
    H.init(0);
    int N = n1/2;
    for (int i=0; i < n1; i++)
    {
        double x = (double) i;  
        double r = (ABS(x) - N) / (double) N;
        if (r <= 0) H(i) = 1;
        else if (r >= 1) H(i) = 0;
        else 
        {
            double xw = 1. - r;
            double lw = exp(1. -1./(1-exp(1.-1./(1-xw))));
            double rw = exp(1. -1./(1-exp(1.-1./xw)));
            double norm = double(sqrt(lw*lw+rw*rw));
            lw = lw/  norm;
            H(i) = lw;
        }
    }
}


void get_wp_meyer_filter(int nside, fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int  Lmax)
{
	int i,p,NbrBand=0;
	double *PtrTab=NULL;
	fltarray TabCl;
    // cout << "get_wp_meyer_filter: " <<nside << endl;
	
	if (nside == 2048)
	{
		NbrBand= TNbrWin[0];
		PtrTab = TabWin2048;
	}
	else if (nside == 1024)
	{
		NbrBand= TNbrWin[1];
		PtrTab = TabWin1024;
	}
	else if (nside == 512)
	{
		NbrBand= TNbrWin[2];
		PtrTab = TabWin512;
	}
	else if (nside == 256)
	{
		NbrBand= TNbrWin[3];
		PtrTab = TabWin256;
	}
	else if (nside == 128)
	{
		NbrBand= TNbrWin[4];
		PtrTab = TabWin128;
	}
	else if (nside == 64)
	{
		NbrBand= TNbrWin[5];
		PtrTab = TabWin64;
	}
	else if (nside == 32)
	{
		NbrBand= TNbrWin[6];
		PtrTab = TabWin32;
	}
	else if (nside == 16)
	{
		NbrBand= TNbrWin[7];
		PtrTab = TabWin16;
	}
	else if (nside == 8)
	{
		NbrBand= TNbrWin[8];
		PtrTab = TabWin8;
	}
	intarray TabB;
	
	int l=0;
	while ((l < NbrBand) && (PtrTab[l] <= Lmax)) l++;
	if (l < NbrBand) NbrBand = l;
	
	
 	if (PtrTab[NbrBand-1] == Lmax) TabB.alloc(NbrBand);
 	else TabB.alloc(NbrBand+1);
 	
	for (int i=0; i < NbrBand; i++)
	{
		TabB(i) = PtrTab[i];
        // cout << PtrTab[i] << endl;
	}
	if (PtrTab[NbrBand-1] != Lmax) 
	{
		TabB(NbrBand) = Lmax;
		NbrBand++;
	}
  	int LM = Lmax+1;
  	
 	fltarray H,Res;
 	Res.alloc(LM, NbrBand+1);
	Win.alloc(NbrBand,2);
    // cout << "get_wp_meyer_filter: " <<NbrBand << endl;
	
	H.resize(Lmax);
 	int Np = TabB(0);
	hgmey(Np, H);
    
    Res.init();
	for (p=0; p < H.nx(); p++) Res(p,0) = H(p);
	Win(0,0) = 0;
	Win(0,1) = TabB(0);
	Win(1,0) = TabB(0);
	Win(1,1) = TabB(1);
    
    // cout << " NbrBand = " << NbrBand << " " << H.nx() << endl;
	for (i=1; i < NbrBand; i++)
	{
	    //  cout << "   N = " << TabB(i-1) << " -> " << TabB(i) << endl;
        hgmey(TabB(i)-TabB(i-1), H);
   	    for (p=0; p < TabB(i-1); p++) Res(p,i) = 1.;
 	    for (p=0; p < H.nx(); p++) 
	        if (p+TabB(i-1) < Res.nx())  Res(p+TabB(i-1),i) = H(p);
	}
 	for (p=0; p < Lmax; p++) Res(p, NbrBand) = 1.;
  	
    
    // fltarray TabG;
    TabH.alloc(LM, NbrBand);
    TabG.alloc(LM, NbrBand);
    
    for (i=0; i < NbrBand; i++)
        for (p=0; p < LM; p++) 
        {
            TabH(p,i) = Res(p, i);
            if (i != NbrBand-1) TabG(p,i) = Res(p, i+1) - Res(p, i);
            else TabG(p,i) = 1. - Res(p, i);
        }
    // cout << " TabH = " << Lmax << " " << NbrBand << endl;
	
	for (i=2; i < NbrBand; i++)
	{
		Win(i,0) = TabB(i-2);
		Win(i,1) = TabB(i);	
	}
    
    WP_WPFilter.alloc(LM, NbrBand+1);
    for (int b=0; b <= NbrBand; b++)
    {
        if (b == NbrBand)
            for (int l=0; l < LM; l++)  WP_WPFilter(l,NbrBand) = TabH(l,0);
        else if (b == 0) for (int l=0; l < LM; l++)  WP_WPFilter(l,0) = 1. - TabH(l,NbrBand-1);
        else  for (int l=0; l < LM; l++)  WP_WPFilter(l,b) = TabH(l,NbrBand-b) - TabH(l,NbrBand-b-1);
    }
    
    // cout << " Win = " << Lmax << " " << NbrBand << endl;
    // fits_write_fltarr("xxh1.fits", TabH);
    // fits_write_fltarr("xxg1.fits", TabG);
    // fits_write_fltarr("xxw1.fits", Win);
    //  exit(-1);
}



 
