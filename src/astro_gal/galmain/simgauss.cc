/* Generates an initial density cube for the genus test.	
 * Second version -- smoothing of P(k)\sim k^n realization
 * by convolution.
 * ES, June 2003	*/

#include "Array.h"
#include "IM_IO.h"
#include "NR.h"

int N=32;                  /* grid size, power of 2*/
int N1;                 /* N+1                  */
int NSQ;                /* N*N                  */
int NP;                 /* number of points     */
int N2;                 /* N/2                  */
int N21;                /* N/2+1                */
int NCB;                /* N^3                  */
float TPINSQ;           /* (2\pi/N)^2           */
float POWA=1;				/* spectral amplitude	*/
float pown=-1;				/* spectral index		*/
float R=0;				/* filter radius		*/
float ***dens;          /* final density 		*/
float *dens1;           /* density for fourn    */
unsigned long *nn;               /* size array for fourn */
float *filter;			/* fourn-packed filter	*/
float meand;            /* mean density         */
float maxl;             /* right cube edge      */
float Norm;             /* FFT norm, 2/N^3      */
FILE *out;				/* density file			*/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

float spect();
void gauss(float disp, float *x, float *y);
float spect(float k);

char Name_Imag_Out[256]; /* output file name */
extern int  OptInd;
extern char *OptArg;

Bool Verbose = False;

/***************************************************************************/

void fourn(float *data, unsigned long *nn, int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP

/***************************************************************************/

/***************************************************************************/
 
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options out_result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-G SpectralIndex]\n");
    fprintf(OUTMAN, "            Power index in the simulated Gaussian field distribution.\n");
    fprintf(OUTMAN, "            Default is %5.2f.\n", pown);
    manline();    
    fprintf(OUTMAN, "         [-D Dimension]\n");
    fprintf(OUTMAN, "            Cube size dimension in the simulation.\n");
    fprintf(OUTMAN, "            Default is %2d.\n", N);
    manline();
    fprintf(OUTMAN, "         [-A SpectralAmplitude]\n");
    fprintf(OUTMAN, "             Spectral amplitude.\n");
    fprintf(OUTMAN, "             Default is %5.2f.\n", POWA);
    manline();        
    fprintf(OUTMAN, "         [-c Std]\n");
    fprintf(OUTMAN, "             Convolve the input data with a Gaussian width sigma=Std.\n");

    manline();
    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    manline();
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/

// void init(int argc, char **argv)
// {
//     long seed;
// 
//     if(argc<6) {
//       fprintf(stderr,"Use: %s N n A R out_file [seed]\n",argv[0]);
//       fprintf(stderr,"N: grid size;\n");
//       fprintf(stderr,"n: spectral index;\n");
//       fprintf(stderr,"A: spectral amplitude;\n");
//       fprintf(stderr,"R: filter radius;\n");
//       fprintf(stderr,"out_file: output density file;\n");
//       fprintf(stderr,"seed: (long int) optional random seed\n");
//             exit(1); 
//     }
//     N=atoi(argv[1]); 
//     pown=atof(argv[2]); POWA=atof(argv[3]); R=atof(argv[4]);
//     if((out=fopen(argv[5],"w"))==NULL) {
//        fprintf(stderr,"Cannot open output file %s\n",argv[5]);
//        exit(1);
//     }
//     if(argc==7) {
//         seed=atol(argv[6]); srand48(seed); 
//     }
//     N1=N+1; N2=N/2; 
//     N21=N2+1; 
//     NSQ=N*N; 
//     NP=NCB=N*N*N; 
//     Norm=2.0/NCB; 
//     TPINSQ=4*M_PI*M_PI/NSQ;
//     dens=f3tensor(1,N,1,N,1,N);
// }

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void siminit(int argc, char *argv[])
{
    int c; 
    float Val; 
    int seed;

    while ((c = GetOpt(argc,argv,"D:G:c:A:I:vzZ")) != -1) 
    {
	switch (c) 
        { 
          case 'D': 
	        if (sscanf(OptArg,"%d",&N) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad dimension parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'G': 
	        if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad spectral index parameter: %s\n", OptArg);
                    exit(-1);
                }
		pown = Val;
                break;
            case 'c': 
                if (sscanf(OptArg,"%f",&R) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Radius parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	    case 'A': 
                if (sscanf(OptArg,"%f",&POWA) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad spectral amplitude parameter: %s\n", OptArg);
                    exit(-1);
                }
                 break;
           case 'I': 
                if (sscanf(OptArg,"%d",&seed) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad seed parameter: %s\n", OptArg);
                    exit(-1);
                }
		srand48((long) seed);
                 break;
	   case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
       

    if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
    else usage(argv);

	/* make sure there are not too many parameters */
    if (OptInd < argc)
    {
       fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
       exit(-1);
    }
 
    N1=N+1; 
    N2=N/2; 
    N21=N2+1; 
    NSQ=N*N; 
    NP=NCB=N*N*N; 
    Norm=2.0/NCB; 
    TPINSQ=4*M_PI*M_PI/NSQ;
    dens=f3tensor(1,N,1,N,1,N);
}

/*********************************************************************/

void gen_filter()
{
    int i,j,k,ik,jk,kk,ir,jr,kr;
    int ii,iir;
    // int ii1,iir1;
    float iksq,jksq,ksq;
    // float kphys,p;
    float RSQ;
    float a,b,tt;
    int izero,jzero,kzero;

    filter= (float *) fvector(1,2*NCB);
	RSQ=R*R;
    for(i=ii=1;i<=N;i++) {
        ik=i-1; if(ik>N2) ik-=N; iksq=ik*ik;
        ir=(N1-i)%N+1;
        if(ik>=N2||ik<=-N2||(ik>0&&ik<1)||(ik<0&&ik>-1)) izero=1; else izero=0;
        for(j=1;j<=N;j++) {
            jk=j-1; if(jk>N2) jk-=N; jksq=iksq+jk*jk;
            jr=(N1-j)%N+1;
            if(jk>=N2||jk<=-N2||(jk>0&&jk<1)||(jk<0&&jk>-1)) jzero=1; 
	    	else jzero=0;
            for(k=1;k<=N21;k++) {
                ii=k+(j-1)*N+(i-1)*NSQ; ii=ii+ii-1;
                if(ii==1) { filter[1]=filter[2]=0.0; 
                      continue; }
                kk=k-1; kr=(N1-k)%N+1;
                iir=kr+(jr-1)*N+(ir-1)*NSQ; iir=iir+iir-1;
                if(kk>=N2||(kk>0&&kk<1)) kzero=1; else kzero=0;
                if(izero || jzero || kzero) a=b=0.0; 
                else { 
                    ksq=TPINSQ*(jksq+kk*kk);
// Gaussian is a real-valued filter
					a=exp(-RSQ*ksq/2.);
// We keep b for other possible cases 
					b=0.;
                    } 
                filter[ii++]=filter[iir++]=a;
                if(ii==iir) tt=0.; else tt=b;
// This line supposes a real filter in real space
                filter[ii]=tt; filter[iir]=-tt;
                }
            }
        }
}

void gen_dens()
{
    int i,j,k,ik,jk,kk,ir,jr,kr;
    int ii,iir;
    float iksq,jksq,ksq,kphys;
    float a,b,p,tt;
    int izero,jzero,kzero;
    int i1;
    float dre,dim,fre,fim;

    nn= (unsigned long *) fvector(1,3); nn[1]=nn[2]=nn[3]=N;
    dens1=(float*) fvector(1,2*NCB);
    for(i=ii=1;i<=N;i++) 
    {
        ik=i-1; 
	if(ik>N2) ik-=N; iksq=ik*ik;
        ir=(N1-i)%N+1;
        if(ik>=N2||ik<=-N2||(ik>0&&ik<1)||(ik<0&&ik>-1)) izero=1; 
	else izero=0;
        for(j=1;j<=N;j++) 
	{
            jk=j-1; if(jk>N2) jk-=N; jksq=iksq+jk*jk;
            jr=(N1-j)%N+1;
            if(jk>=N2||jk<=-N2||(jk>0&&jk<1)||(jk<0&&jk>-1)) jzero=1; 
	    else jzero=0;
	    
            for(k=1;k<=N21;k++) 
	    {
                ii=k+(j-1)*N+(i-1)*NSQ; ii=ii+ii-1;
                if(ii==1) 
		{ 
		    dens1[1]=dens1[2]=0.0; 
                    continue; 
		}
                kk=k-1; kr=(N1-k)%N+1;
                iir=kr+(jr-1)*N+(ir-1)*NSQ; 
		iir=iir+iir-1;
                if(kk>=N2||(kk>0&&kk<1)) kzero=1; 
		else kzero=0;
                if(izero || jzero || kzero) a=b=0.0; 
                else 
		{ 
                    ksq=jksq+kk*kk;
                    kphys=2*M_PI/N*sqrt(ksq);
                    p=spect(kphys);
                    gauss(0.5*p,&a,&b);
                } 
                dens1[ii++]=dens1[iir++]=a;
                if(ii==iir) tt=0.; 
		else tt=b;
                dens1[ii]=tt; 
		dens1[iir]=-tt;
            } // end for(k=1;k<=N21;k++)
        } // end for(j=1;j<=N;j++) 
    } // end for(i=ii=1;i<=N;i++)

   // convolution, valid for general complex filters 
   if (R > 0)
   {
      for(i=1;i<=2*NCB;i+=2) 
      {
         i1=i+1;
         dre=dens1[i]; dim=dens1[i1];	
         fre=filter[i]; fim=filter[i1];
         dens1[i]=dre*fre-dim*fim;
         dens1[i1]=dre*fim+fre*dim;
      }
   }

   // Fourier transform
    fourn(dens1,nn,3,-1);
    for(i=ii=1;i<=N;i++)  
    for(j=1;j<=N;j++) 
    for(k=1;k<=N;k++,ii+=2) 
        dens[i][j][k]=dens1[ii]; 
//    free_vector(dens1,1,2*NCB);
}

void gauss(float disp, float *x, float *y)    /* Recipes */
{
    float v1,v2,r,fac;
    // double drand48();

    do {
        v1=2.0*drand48()-1; 
        v2=2.0*drand48()-1; 
        r=v1*v1+v2*v2;
        }
    while(r>=1.0);
    fac=sqrt(-2*disp*log((double) r)/r);
    *x=v1*fac;
    *y=v2*fac;
}
    /* Power spectrum */
float spect(float k)
{
    float P;	

    P=POWA*pow((double) k,(double) pown);
    return(P);
}

// write the output array
// void output()
// {
// 	int i,j,k;
// 	float *dp;
// 	float d,mean,rms;
// 
// 	mean=rms=0.;
// 	for(i=1,dp=dens1+1;i<=N;i++)
// 	for(j=1;j<=N;j++) {
// 		for(k=1;k<=N;k++) {
// 			d=dens[i][j][k];
// 			mean+=d; rms+=d*d;
// //		*dp++=d;
// 			fprintf(out,"%g ",d);
// 			}
// 		fprintf(out,"\n",d);
// 		}
// 	mean/=NCB; rms=sqrt(rms/NCB-mean*mean);
// 	fprintf(stderr,"mean: %g  rms: %g\n",mean,rms);
// //	if(!fwrite(dens1+1,NCB*sizeof(float),1,out)) {
// //		fprintf(stderr,"Error while writing the output file\n");
// //		exit(1);
// //	}
// }

int main(int argc, char ** argv)
{
    int i,j,k;
    // init(argc,argv);
    siminit(argc,argv);
    if (Verbose == True)
    { 
        cout << endl << endl << "# PARAMETERS: " << endl << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
	cout << "# Spectral Index = " << pown   << endl;  
	cout << "# Spectral Amplitude = " << POWA  << endl;
	cout << "# Dimension = " << N  << endl;
	cout << "# Radius R = " << R << endl;
    }
    if (R > 0) gen_filter();
    gen_dens();
    fltarray Tab(N, N, N);	
    for(i=1;i<=N;i++)
    for(j=1;j<=N;j++) 
    for(k=1;k<=N;k++) Tab(i-1,j-1,k-1) = dens[i][j][k];
    if (Verbose == True)
    {
       cout << " Mean = " << Tab.mean() << " Sigma = " << Tab.sigma() << endl;
    }
    fits_write_fltarr(Name_Imag_Out, Tab);
    // output();
    exit(0);
}
