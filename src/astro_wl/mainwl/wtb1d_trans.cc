/******************************************************************************
**                   Copyright (C) 1999 by CEA
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
**    File:  dct.cc
**
*******************************************************************************
**
**    DESCRIPTION  dct program
**    ----------- 
**                 
******************************************************************************/

#include "GlobalInc.h"
#include "IM_IO.h"
 

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool Reverse=False;
Bool Linear=False;
int N;

#define pi 3.1415926535898
#define maxf 7  /* le filtre est contenu dans [-maxf,+maxf]  */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


void tfd1D_br_free_V0(float *c,float *d,float *h1,float *g1)
{
    int j,i,l;
    float *k;
    k=(float*) malloc(sizeof(float)*(N+1));
    
    for (i=0;i<=N;i++) {
        k[i]=c[i];
    }
    for (j=N/2;j>1;j=j/2) {
        for (i=2;i<(j-1);i++) {
            d[i]=0.;
            for (l=-maxf;l<=maxf;l++) d[i]+=h1[l+maxf]*k[(2*i+l+maxf*(2*j+1))%(2*j+1)];
        }
        for (i=1;i<(j-1);i++) {
            d[j+1+i]=0.;
            for (l=-maxf;l<=maxf;l++) d[j+1+i]+=g1[l+maxf]*k[(2*i+l+maxf*(2*j+1))%(2*j+1)];
        }
        // exception pour les bords ds V0 !
        d[j+1]=-0.25*k[0]+0.5*k[1]-0.25*k[2];
        d[0]=k[0]+0.75*d[j+1];
        d[1]=0.25*k[1]+0.75*k[2]+0.25*k[3]-1./8.*k[4]-1./8.*d[0]-1./32.*d[j+1];
        
        d[2*j]=-0.25*k[2*j]+0.5*k[2*j-1]-0.25*k[2*j-2];
        d[j]=k[2*j]+0.75*d[2*j];
        d[j-1]=0.25*k[2*j-1]+0.75*k[2*j-2]+0.25*k[2*j-3]-1./8.*k[2*j-4]-1./8.*d[j]-1./32.*d[2*j];
        if (j==2) d[1]=0.25*k[2*j-1]+0.75*k[2*j-2]+0.25*k[2*j-3]-1./8.*d[0]-1./32.*d[j+1]-1./8.*d[j]-1./32.*d[2*j];
        
        for (i=0;i<=j;i++) k[i]=d[i];
    }
    d[0]=1./8.*k[0]+0.75*k[1]+1./8.*k[2];
    d[1]=7./12.*k[0]-0.5*k[1]-1./12.*k[2];
    d[2]=-1./12.*k[0]-0.5*k[1]+7./12.*k[2];
    free(k);
}

void tfd1Dinv_br_free_V0(float *c,float *d,float *h1,float *g1)
{
    int j,i,l;
    float *cw;
    cw=(float*) malloc(sizeof(float)*(N+1));
    
    for (i=0;i<=N;i++) {
        cw[i]=d[i];
        c[i]=d[i];
    }
    c[0]=d[0]+1.5*d[1];
    c[1]=d[0]-0.25*d[1]-0.25*d[2];
    c[2]=d[0]+1.5*d[2];
    for (j=2;j<N;j=2*j) {
        for (i=0;i<=j;i++) cw[i]=c[i];
        for (i=0;i<=2*j;i++) c[i]=0.;
        for (i=1;i<j;i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+l+maxf*(2*j+1))%(2*j+1)]+=h1[l+maxf]*cw[i];
            }
        }
        for (i=1;i<(j-1);i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+l+maxf*(2*j+1))%(2*j+1)]+=g1[l+maxf]*d[j+1+i];
            }
        }
        c[0]+=-0.75*d[j+1]+cw[0];
        c[1]+=23./16.*d[j+1]+0.5*cw[0];
        c[2]+=-3./8.*d[j+1];
        c[3]+=-3./16.*d[j+1];
        c[2*j]+=-0.75*d[2*j]+cw[j];
        c[2*j-1]+=23./16.*d[2*j]+0.5*cw[j];
        c[2*j-2]+=-3./8.*d[2*j];
        c[2*j-3]+=-3./16.*d[2*j];
    }
    free(cw);
}

/***********************   Ondelettes dans V1   **********************************************/


void tfd1D_br_V1(float *c,float *d,float *h1,float *g1)
{
    int j,i,l;
    float *k;
    k=(float*) malloc(sizeof(float)*(N+2));
    
    for (i=0;i<=(N+1);i++) {
        k[i]=c[i];
    }
    for (j=N/2;j>1;j=j/2) {
        for (i=2;i<j;i++) {
            d[i]=0.;
            for (l=-maxf;l<=maxf;l++) d[i]+=h1[l+maxf]*k[(2*i-1+l+maxf*(2*j+2))%(2*j+2)];
        }
        for (i=1;i<(j-1);i++) {
            d[j+2+i]=0.;
            for (l=-maxf;l<=maxf;l++) d[j+2+i]+=g1[l+maxf]*k[(2*i+1+l+maxf*(2*j+2))%(2*j+2)];
        }
        // exception pour les bords ds V1 !
        d[j+2]=-0.25*k[0]+0.5*k[1]-3./8.*k[2]+1./8.*k[3];
        d[0]=k[0]-d[j+2];
        d[1]=1.5*k[2]-0.5*k[3]+1.5*d[j+2];
        
        d[2*j+1]=-0.25*k[2*j+1]+0.5*k[2*j]-3./8.*k[2*j-1]+1./8.*k[2*j-2];
        d[j+1]=k[2*j+1]-d[2*j+1];
        d[j]=1.5*k[2*j-1]-0.5*k[2*j-2]+1.5*d[2*j+1];
        
        for (i=0;i<=(j+1);i++) k[i]=d[i];
    }
    d[0]=1./6.*k[0]+1./3.*k[1]+1./3.*k[2]+1./6.*k[3];
    d[1]=-1./8.*k[0]-0.25*k[1]+0.25*k[2]+1./8.*k[3];
    d[2]=7./12.*k[0]-5./6.*k[1]+1./6.*k[2]+1./12.*k[3];
    d[3]=1./12.*k[0]+1./6.*k[1]-5./6.*k[2]+7./12.*k[3];
    free(k);
}

void tfd1Dinv_br_V1(float *c,float *d,float *h1,float *g1)
{
    int j,i,l;
    float *cw;
    cw=(float*) malloc(sizeof(float)*(N+2));
    
    for (i=0;i<=N;i++) {
        cw[i]=d[i];
        c[i]=d[i];
    }
    c[0]=d[0]-2.*d[1]+d[2];
    c[1]=d[0]-d[1]-0.5*d[2];
    c[2]=d[0]+d[1]-0.5*d[3];
    c[3]=d[0]+2.*d[1]+d[3];
    for (j=2;j<N;j=2*j) {
        for (i=0;i<=(j+1);i++) cw[i]=c[i];
        for (i=0;i<=(2*j+1);i++) c[i]=0.;
        for (i=2;i<j;i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i-1+l+maxf*(2*j+2))%(2*j+2)]+=h1[l+maxf]*cw[i];
            }
        }
        for (i=1;i<(j-1);i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+1+l+maxf*(2*j+2))%(2*j+2)]+=g1[l+maxf]*d[j+2+i];
            }
        }
        c[0]+=d[j+2]+cw[0];
        c[1]+=7./4.*d[j+2]+0.5*cw[1]+0.5*cw[0];
        c[2]+=-9./8.*d[j+2]+0.75*cw[1];
        c[3]+=-3./8.*d[j+2]+0.25*cw[1];
        c[2*j+1]+=d[2*j+1]+cw[j+1];
        c[2*j]+=7./4.*d[2*j+1]+0.5*cw[j]+0.5*cw[j+1];
        c[2*j-1]+=-9./8.*d[2*j+1]+0.75*cw[j];
        c[2*j-2]+=-3./8.*d[2*j+1]+0.25*cw[j];
    }
    free(cw);
}


/***************************************/

static void usage(char *argv[]
)
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_data \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
   fprintf(OUTMAN, "         [-r]\n");
   fprintf(OUTMAN, "             Inverse transform.\n");
   manline();
    
    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "             Linear wavelets. Default is quadratic wavelets.\n");
    fprintf(OUTMAN, "               Linear wavelets: input signal must be of size N+1, where N is a power of 2.\n");
    fprintf(OUTMAN, "               Quadratic wavelets: input signal must be of size N+2, where N is a power of 2.\n");
    manline();
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();

    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv,"LrvzZ")) != -1) 
    {
	switch (c) 
        {
            case 'L': Linear = True; break;
            case 'r': Reverse = True; break;
            case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
 
       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

       if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}


 
 
/*********************************************************************/
 
int main(int argc, char *argv[])
{
    int i,j,k, ind=0;
      
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl; 
        if (Linear == True)  cout << "Linear wavelets " << endl; 
        else  cout << "Quadratic wavelets " << endl; 
    }
 
    fltarray Data;
    
    fits_read_fltarr(Name_Imag_In, Data);
    fltarray Result;

    float Min,Max;
    int Nx = Data.nx();
    N = Nx;
    if (Linear == False) N = Nx-2;
    else N=Nx-1;
    Result.alloc(Nx);
    if (Verbose == True) cout << "Nx = " << Nx << " N = " << N <<  endl;
    if (Verbose == True) cout << "Min = " << Data.min()  << "  Max = " << Data.max() <<  " Sigma = " << Data.sigma() << endl;
    
    
    float g1d[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
    float h1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
    float h1d[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
    float g1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};
    
    float g2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
    float h2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
    float h2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
    float g2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};
    float *h,*g;
    double NormH=0;
    double NormG=0;
    if (Linear == False)
    {
       h = h2d;
       g = g2d;
    }
    else 
    {
       h = h1d;    
       g = g1d; 
    }
    for (int i=0; i < 2*maxf+1; i++) NormH += h[i]*h[i];
    NormH = 1. / sqrt(NormH);
    for (int i=0; i < 2*maxf+1; i++) NormG += g[i]*g[i];
    NormG = 1. / sqrt(NormG);
    
    fltarray T;
    T.alloc(Nx);
    T = Data;
    float *c = Data.buffer();  // pointer to the data
    float *d = Result.buffer();
    
    int Np = Nx;
    int s=1;
    int Nw = N/2;
     
    for (int j=N/2;j>1;j=j/2) 
    {
       cout << "scale " << s << " " << Nw << " " << NormH << " " << NormG << endl;
       Nw /= 2;
       s++;
    }

    
    int Debug=0;
    if (Reverse == False)
    {
       if (Linear == False) tfd1D_br_V1(T.buffer(),d, h2d, g2d);
       else tfd1D_br_free_V0(T.buffer(),d, h1d, g1d);
       if (Debug == 1)
       {
         Result.info("RR");
          T.init();
          if (Linear == False)  tfd1Dinv_br_V1(T.buffer(),d, h2i, g2i);
          else tfd1Dinv_br_free_V0(T.buffer(),d, h1i, g1i);
          fits_write_fltarr("xx_wtb.fits",  Result);
          Data -= T;
          Data.info("Resi");
          fits_write_fltarr(Name_Imag_Out,  Data);
         exit(0);
       }
       
       // Normalization
       double ValN = 2.*NormG;
       int ind = Nx;
       int Nw = N/2;
       s=1;
       for (int j=N/2;j>1;j=j/2) 
       {
          cout << "scale " << s << " " << Nw << " " << NormH << " " << NormG << endl;
          for (int p=0; p < Nw; p++) Result(ind--) *= ValN;
          ValN *= NormH;
          Nw /= 2;
          s++;
        }
    }
    else 
    {
        if (Linear == False) tfd1Dinv_br_V1(d,c, h2i, g2i);
        else tfd1Dinv_br_free_V0(d,c, h1i, g1i);
    }
    
    
          
                
    fits_write_fltarr(Name_Imag_Out,  Result);
    exit(0);
}

/*********************************************************************/


