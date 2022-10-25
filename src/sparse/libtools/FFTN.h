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
**    Date:  28/07/00 
**    
**    File:  FFTN.h
**
*******************************************************************************
**
**    DESCRIPTION: FFT Class for 1D,2D,3D FFT with array of any size.  
**    ----------- 
**                 
**
******************************************************************************/

#ifndef _DEF_FFTN_H_
#define _DEF_FFTN_H_

#include "GlobalInc.h"
#ifdef USE_FFTW
#include <fftw3.h>
#endif

// void fft_free (void);

/* double precision routine */
// int fftn (int ndim, const int dims[], double Re[], double Im[],  int isign, double scaling);

/* float precision routine */
// int fftnf (int ndim, const int dims[], float Re[], float Im[],  int isign, double scaling);

#define NFACTOR	11
#define FFT_FLOAT 1
extern "C" {
    extern void fft_free (void);
    
    /* double precision routine */
    extern int fftn (int /* ndim */,
                     const int /* dims */[],
                     double /* Re */[],
                     double /* Im */[],
                     int /* isign */,
                     double /* scaling */);
    
    /* float precision routine */
    extern int fftnf (int /* ndim */,
                      const int /* dims */[],
                      float /* Re */[],
                      float /* Im */[],
                      int   /* isign */,
                      float /* scaling */);
}

#define  FORWARD_FFT -2
#define  REVERSE_FFT 2


class FFTN
{
		void fft_error() 
		{
		   cout << "Error in FFT computation ... " << endl;
		   exit(-1);
		}
		void scale_reverse(float *Buff, int N)
		{
		   double C = 1. / (double) N;
		   for (int i=0; i < 2*N; i++) Buff[i] *= C;
		}
		void scale_reverse(double *Buff, int N)
		{
		   double C = 1. / (double) N;
		   for (int i=0; i < 2*N; i++) Buff[i] *= C;
		}
	public:
		
	// Utilities
		int factor [NFACTOR];
		/* temp space, (void *) since both float and double routines use it */
		void *Tmp0;  /* temp space for real part */
		void *Tmp1;  /* temp space for imaginary part */
		void *Tmp2;  /* temp space for Cosine values */
		void *Tmp3;  /* temp space for Sine values */
		int  *Perm;  /* Permutation vector */

		// double precision routines
		/* 
        int fftn (int ndim, const int dims[], double Re[], double Im[], int isign, double scaling);
		int fftradix (double Re[], double Im[], size_t nTotal, size_t nPass, size_t nSpan, int isign, int max_factors, int max_perm);
		// float precision routines
		int fftnf (int ndim, const int dims[], float Re[], float Im[], int isign, double scaling);
		int fftradixf (float Re[], float Im[], size_t nTotal, size_t nPass, size_t nSpan, int isign, int max_factors, int max_perm);
        */
		// parameters for memory management 
		size_t SpaceAlloced;
		size_t MaxPermAlloced;

	// Transform functions
		void transform1d(float *Buff, int N, Bool Reverse=False, bool normalize=false)
		// Complex 1D FFT transform
		// N: IN = number of elements
		// CFBuff: IN-OUT = buffer [0..2N-1]
		//                   CFBuff[2i] = real part
		//                   CFBuff[2i+1] = imaginary part
		// Reverse: IN = True for inverse FFT transform
             {
                float *Pr = Buff;
                float *Pi = Buff + 1;
                int  Ndim = 1;
                int dims[1];
                double Scaling = 0.;
                int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT; 
                dims[0] = N;
                int ret = fftnf (Ndim, dims, Pr, Pi, isign, Scaling);
                if (Reverse == True) scale_reverse(Buff, N);
                if (ret != 0) fft_error();		 
                if(normalize)
                {
                    float val = (Reverse==True) ? sqrt(N) : 1./sqrt(N);
                    for(int i = 0;i<2*N;i++)
                        Buff[i] *= val;
                }
            };
            

    // Ditto with double
	void  transform1d(double *Buff, int N, Bool Reverse=False, bool normalize=false)
    {
        double *Pr = Buff;
        double *Pi = Buff + 1;
        int  Ndim = 1;
        int dims[1];
        double Scaling = 0.;
        int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT; 
        dims[0] = N;
        int ret = fftn (Ndim, dims, Pr, Pi, isign, Scaling);		 
        if (Reverse == True) scale_reverse(Buff, N);
        if (ret != 0) fft_error();		 
        if(normalize)
        {
            double val = (Reverse==True) ? sqrt(N) : 1./sqrt(N);
            for(int i = 0;i<2*N;i++)
                Buff[i] *= val;
        }
    };
    

    // Complex 2D FFT transform
    // Nx: IN = number of columns
    // Ny: IN = number of lines
    // CFBuff: IN-OUT = buffer [0..N-1] with N = 2NxNy 
    //                   CFBuff[2i] = real part
    //                   CFBuff[2i+1] = imaginary part
    // Reverse: IN = True for inverse FFT transform    
    void transform2d(float *Buff, int Nx,  int Ny, Bool Reverse=False, bool normalize=false)
    {
        float *Pr = Buff;
        float *Pi = Buff + 1;
        int  Ndim = 2;
        int dims[2];
        double Scaling = 0.;
        int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT;
        dims[0] = Nx;
        dims[1] = Ny;
        int ret = fftnf (Ndim, dims, Pr, Pi, isign, Scaling);
        if (Reverse == True) scale_reverse(Buff, Nx*Ny);
        if (ret != 0) fft_error();		 
        if(normalize)
        {
            float val = (Reverse==True) ? sqrt(Nx*Ny) : 1./sqrt(Nx*Ny);
            for(int i = 0;i<2*Nx*Ny;i++)
                Buff[i] *= val;
        }
    };

 		// Ditto with double
    void  transform2d(double *Buff, int Nx,  int Ny, Bool Reverse=False, bool normalize=false)
    {
        double *Pr = Buff;
        double *Pi = Buff + 1;
        int  Ndim = 2;
        int dims[2];
        double Scaling = 0.;
        int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT; 
        dims[0] = Nx;
        dims[1] = Ny;
        int ret = fftn (Ndim, dims, Pr, Pi, isign, Scaling);		 
        if (Reverse == True) scale_reverse(Buff, Nx*Ny);
        if (ret != 0) fft_error();		 
        if(normalize)
        {
            double val = (Reverse==True) ? sqrt(Nx*Ny) : 1./sqrt(Nx*Ny);
            for(int i = 0;i<2*Nx*Ny;i++)
                Buff[i] *= val;
        }
    };

 		// Complex 2D FFT transform
		// Nx,Ny,Nz: IN = array dimensions
		// CFBuff: IN-OUT = buffer [0..N-1] with N = 2NxNyNz
		//                   CFBuff[2i] = real part
		//                   CFBuff[2i+1] = imaginary part
		// Reverse: IN = True for inverse FFT transform
    void transform3d(float *Buff, int Nx,  int Ny, int Nz, Bool Reverse=False, bool normalize=false)
    {
        float *Pr = Buff;
        float *Pi = Buff + 1;
        int  Ndim = 3;
        int dims[3];
        double Scaling = 0.;
        int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT; 
        dims[0] = Nx;
        dims[1] = Ny;
        dims[2] = Nz;
        int ret = fftnf (Ndim, dims, Pr, Pi, isign, Scaling);
        if (Reverse == True) scale_reverse(Buff, Nx*Ny*Nz);
        if (ret != 0) fft_error();		 
        if(normalize)
        {
            float val = (Reverse==True) ? sqrt(Nx*Ny*Nz) : 1./sqrt(Nx*Ny*Nz);
            for(int i = 0;i<2*Nx*Ny*Nz;i++)
                Buff[i] *= val;
        }
    };
 		// Ditto with double
    void  transform3d(double *Buff, int Nx,  int Ny, int Nz, Bool Reverse=False, bool normalize=false)
    {
        double *Pr = Buff;
        double *Pi = Buff + 1;
        int  Ndim = 3;
        int dims[3];
        double Scaling = 0;
        int isign  = (Reverse == False) ? FORWARD_FFT: REVERSE_FFT; 
        dims[0] = Nx;
        dims[1] = Ny;
        dims[2] = Nz;
        int ret = fftn (Ndim, dims, Pr, Pi, isign, Scaling);		 
        if (Reverse == True) scale_reverse(Buff, Nx*Ny*Nz);
        if (ret != 0) fft_error();		 
        if(normalize)
        {
            double val = (Reverse==True) ? sqrt(Nx*Ny*Nz) : 1./sqrt(Nx*Ny*Nz);
            for(int i = 0;i<2*Nx*Ny*Nz;i++)
                Buff[i] *= val;
        }
    };
    FFTN()
    {
        SpaceAlloced = 0;
        MaxPermAlloced = 0;
        Tmp0 = NULL;
        Tmp1 = NULL;
        Tmp2 = NULL;
        Tmp3 = NULL;
        Perm = NULL;
    };
    
    ~FFTN()
    {
        SpaceAlloced = MaxPermAlloced = 0;
        if (Tmp0 != NULL)	{ free (Tmp0);	Tmp0 = NULL; }
        if (Tmp1 != NULL)	{ free (Tmp1);	Tmp1 = NULL; }
        if (Tmp2 != NULL)	{ free (Tmp2);	Tmp2 = NULL; }
        if (Tmp3 != NULL)	{ free (Tmp3);	Tmp3 = NULL; }
        if (Perm != NULL)	{ free (Perm);	Perm = NULL; }
    };
};


#endif
