
extern "C"  {                   

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define Re(arr,j) arr[(j<<1)]
#define Im(arr,j) arr[(j<<1)+1]

#define TWOPI 6.28318530717959


//SS  FFTmid0  ifftcol  fftshift

#define YR(i,j) yr[(i) + M*(j)]
#define YI(i,j) yi[(i) + M*(j)] 
#define ZR(i,j) zr[(i) + M*(j)] 
#define ZI(i,j) zi[(i) + M*(j)] 
#define XR(i,j) xr[(i) + N*(j)]
#define XI(i,j) xi[(i) + N*(j)] 

//Adj_PPFFT

#define X0R(i,j) x0r[(i) + M*(j)]
#define X0I(i,j) x0i[(i) + M*(j)] 
#define X1R(i,j) x1r[(i) + M*(j)]
#define X1I(i,j) x1i[(i) + M*(j)] 
#define Y0R(i,j) yr[(i) + M*(j)]
#define Y0I(i,j) yi[(i) + M*(j)] 
#define Y1R(i,j) yr[(i+N) + M*(j)]
#define Y1I(i,j) yi[(i+N) + M*(j)] 
#define Z0R(i,j) z0r[(i+HalfN) + M*(j)]
#define Z0I(i,j) z0i[(i+HalfN) + M*(j)] 
#define Z1R(i,j) z1r[(i) + N*(j+HalfN)]
#define Z1I(i,j) z1i[(i) + N*(j+HalfN)] 


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr



#define BYROW  1
#define BYCOL  0
#define REAL double


void SS( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N);

void PPFFT(REAL zr[],REAL zi[],REAL xr[], REAL xi[],
						 int N,int preval);
void IFFTmid0( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N,
                   int M,
				   int rowcol);                         
                         

void copy_sub( REAL yr[],REAL yi[],
               REAL xr[],REAL xi[], 
               int N);
                   


void FFFT( REAL zr[],REAL zi[],
                        REAL xr[],REAL xi[],
                        int N,
                        const REAL alpha);

void FFTmid0( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
				   int N,
				   int rowcol);

void ifftcol( REAL yr[], REAL yi[], REAL tmp[], 
             int M,
             int N);
void fftshift(REAL zr[],REAL zi[],
              int M,
			  int N); 

void four1(REAL data[], int nn, int isign);

void fftcol( REAL yr[], REAL yi[], REAL tmp[], 
             int M,
             int N);

void Adj_SS( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N);
void Inv_SS( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N, int Verb);

void Adj_PPFFT(REAL zr[],REAL zi[],REAL xr[], REAL xi[],
						 int N,int preval);

void FFTmid0NP( REAL yr[],REAL yi[],
                   int M,
				   int N);                         
                                            
                         
}
