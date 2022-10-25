 
#include "FSS.h" 

extern "C"  {

int NUM_OF_ITER = 5; // For the Iterative Inversion
                   
/* Main functions are 
   SS - Forward 2D Radon Transform, gets an N*N image and computes a 
        2N*2N array of Radon coefficients
        xr, xi - arrays of the real and imaginary parts of the 2D image (N^2 long each)
		yr, yi - arrays to which the Radon Coefficients will be written (4N^2 long each)
		N - Image side length

  Adj_SS - Adjoint 2D Radon Transform, gets an 2N*2N array of Radon coefficients
        and computes an N*N array of the Adjoint Radon coefficients
        xr, xi - arrays of the real and imaginary parts of the 2D radon coefficients (4N^2 long each).
		yr, yi - arrays to which the Adjoint Radon Coefficients will be written (N^2 long each).
		N - Image side length 
		
  Inv_SS - Inverse 2D Radon Transform, gets an 2N*2N array of Radon coefficients
        and computes the inverse N*N image
        xr, xi - arrays of the real and imaginary parts of the 2D radon coefficients (4N^2 long each).
		yr, yi - arrays to which the Inverse Radon Coefficients will be written (N^2 long each).
		N - Image side length 

		The Inverse is computed by CG iterations, default number of iteration is 5
		it can be changed in the define line above
*/




/*  
    FastSlantStack fin_name fout_name n mode
	fin_name - input file name
	fout_name - output file name
	n - image array size 
	mode - 1 for forward transform 2 for adjoint transform  3 for the inverse transform */


/*
main (int argc, char *argv[])
{
	int n, mode;
	REAL *xR,*xI,*yR,*yI;
	FILE *fopen(), *fin,*fout;

	if (argc!=5){
		printf("Error! inproper numer of arguments. \n");
		exit(-1);
	}
	n=atoi(argv[3]);
	mode=atoi(argv[4]);
	if (mode == 1){
		xR=(REAL*)malloc(sizeof(REAL)*n*n);
		xI=(REAL*)calloc(n*n,sizeof(REAL));
		yR=(REAL*)malloc(sizeof(REAL)*n*n*4);
		yI=(REAL*)malloc(sizeof(REAL)*n*n*4);
		if ((fin=fopen(argv[1],"rb"))==NULL) {
			printf("Cannot open file. \n");
			exit(1);
			}
	
		fread(xR,sizeof(REAL),n*n,fin);
		fread(xI,sizeof(REAL),n*n,fin);
		fclose(fin);
    	SS(yR,yI,xR,xI,n); 
		fout=fopen(argv[2],"wb");
		fwrite(yR,sizeof(REAL),4*n*n,fout);
		fwrite(yI,sizeof(REAL),4*n*n,fout);
		fclose(fout);
	}

	if (mode == 2){
		xR=(REAL*)malloc(sizeof(REAL)*4*n*n);
		xI=(REAL*)calloc(4*n*n,sizeof(REAL));
		yR=(REAL*)malloc(sizeof(REAL)*n*n);
		yI=(REAL*)malloc(sizeof(REAL)*n*n);
		if ((fin=fopen(argv[1],"rb"))==NULL) {
			printf("Cannot open file. \n");
			exit(1);
			}
	
		fread(xR,sizeof(REAL),4*n*n,fin);
		fread(xI,sizeof(REAL),4*n*n,fin);
		fclose(fin);
    	Adj_SS(yR,yI,xR,xI,n); 
		fout=fopen(argv[2],"wb");
		fwrite(yR,sizeof(REAL),n*n,fout);
		fwrite(yI,sizeof(REAL),n*n,fout);
		fclose(fout);
	}
		
if (mode == 3){
		xR=(REAL*)malloc(sizeof(REAL)*4*n*n);
		xI=(REAL*)calloc(4*n*n,sizeof(REAL));
		yR=(REAL*)malloc(sizeof(REAL)*n*n);
		yI=(REAL*)malloc(sizeof(REAL)*n*n);
		if ((fin=fopen(argv[1],"rb"))==NULL) {
			printf("Cannot open file. \n");
			exit(1);
			}
	
		fread(xR,sizeof(REAL),4*n*n,fin);
		fread(xI,sizeof(REAL),4*n*n,fin);
		fclose(fin);
    	Inv_SS(yR,yI,xR,xI,n); 
		fout=fopen(argv[2],"wb");
		fwrite(yR,sizeof(REAL),n*n,fout);
		fwrite(yI,sizeof(REAL),n*n,fout);
		fclose(fout);
	}
free (xR);
free (xI);
free (yR);
free (yI);
}
*/
void SS( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N)
{  
    int i,j;
	REAL *zr,*zi;
    
    PPFFT(yr,yi,xr,xi,N,0);
    zr=(REAL*)malloc(sizeof(REAL)*N*N*4);
	zi=(REAL*)malloc(sizeof(REAL)*N*N*4); 
    
    copy_sub(zr,zi,yr,yi,N*N*4);
    for( j= 0; j < 2*N; j++)
				for( i=0; i< 2*N; i++){

					zr[(j) + 2*N*(i)] = yr[(i) + 2*N*(j)];
					zi[(j) + 2*N*(i)] = yi[(i) + 2*N*(j)];
				}
                
    IFFTmid0(yr,yi,zr,zi,2*N,2*N,BYCOL);
    
    free(zr);
    free(zi);
    }


void Adj_SS( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N)
{  
 
	REAL *zr, *zi;
	int i,j;

	FFTmid0NP(xr,xi,2*N,2*N);
    zr=(REAL*)malloc(sizeof(REAL)*N*N*4);
	zi=(REAL*)malloc(sizeof(REAL)*N*N*4); 
    
	copy_sub(zr,zi,xr,xi,N*N*4);
    for( j= 0; j < 2*N; j++)
				for( i=0; i< 2*N; i++){
					zr[(j) + 2*N*(i)] = xr[(i) + 2*N*(j)];
					zi[(j) + 2*N*(i)] = xi[(i) + 2*N*(j)];
				}
	
	Adj_PPFFT(yr,yi,zr,zi,2*N,0);

free(zr);
free(zi);    
}


void Inv_SS( REAL xr[],REAL xi[],
                   REAL yr[],REAL yi[], 
                   int N, int Verb)
{  
 
	REAL *zr,*zi,*rr, *ri, *pr,*pi,*tr,*ti,*Apr,*Api, weight,rn,rn0,npAp,an,bn;
	int i,j,M,It;
	M=2*N;
	tr=(REAL*)malloc(sizeof(REAL)*4*N*N);
	ti=(REAL*)malloc(sizeof(REAL)*4*N*N); 
	for( j= 0; j < M*M; j++){
		tr[j]=yr[j];
		ti[j]=yi[j];
	}

	FFTmid0NP(tr,ti,M,M);
	for( j= 0; j < M; j++)
		for( i=0; i< M; i++){
			if (i!=N)
				weight=abs(i-N);
            else
				weight=0.5;
			tr[(i) + M*(j)]*=weight;
			ti[(i) + M*(j)]*=weight;
		}
		// tr[0]*=0.5;
		// ti[0]*=0.5;
		//tr[M*(N)]*=0.5;
		// ti[M*(N)]*=0.5;
		
		zr=(REAL*)malloc(sizeof(REAL)*4*N*N);
	    zi=(REAL*)malloc(sizeof(REAL)*4*N*N); 
		for( j= 0; j < 2*N; j++)
				for( i=0; i< 2*N; i++){
					zr[(j) + 2*N*(i)] = tr[(i) + 2*N*(j)];
					zi[(j) + 2*N*(i)] = ti[(i) + 2*N*(j)];
				}
			

    rr=(REAL*)malloc(sizeof(REAL)*N*N);
	ri=(REAL*)malloc(sizeof(REAL)*N*N); 
	pr=(REAL*)malloc(sizeof(REAL)*N*N);
	pi=(REAL*)malloc(sizeof(REAL)*N*N); 
	Apr=(REAL*)malloc(sizeof(REAL)*N*N);
	Api=(REAL*)malloc(sizeof(REAL)*N*N); 	

	Adj_PPFFT(rr,ri,zr,zi,M,0);
	
	free(zr);
	free(zi);


	for (j=0;j<N*N;j++){
		pr[j]=rr[j];
		pi[j]=ri[j];
		xr[j]=0;
		xi[j]=0;}

// printf( "NUM_OF_ITER = %d\n",NUM_OF_ITER);
	for (It=0;It<NUM_OF_ITER;It++){
		PPFFT(tr,ti,pr,pi,N,0);
		for( j= 0; j < M; j++)
			for( i=0; i< M; i++){
				if (j!=N)
					weight=abs(j-N);
				else
					weight=0.5;
				tr[(i) + M*(j)]*=weight;
				ti[(i) + M*(j)]*=weight;
			}
		//tr[0]*=0.5;
		//ti[0]*=0.5;
		// tr[N]*=0.5;
		// ti[N]*=0.5;
		Adj_PPFFT(Apr,Api,tr,ti,M,0);
		rn0=0;
		npAp=0;
		for (j=0;j<N*N;j++){
			rn0+=(rr[j]*rr[j]+ri[j]*ri[j]);
			npAp+=(pr[j]*Apr[j]+pi[j]*Api[j]);
		}
		if (npAp != 0) an=rn0/npAp;
        else an=0;
		if (Verb) printf("It=%d,error=%g\n",It+1,rn0);
		for (j=0;j<N*N;j++){
			xr[j]+=an*pr[j];
			xi[j]+=an*pi[j];
			// if (xr[j] < 0) xr[j] = 0;
			xi[j]=0;
			rr[j]-=an*Apr[j];
			ri[j]-=an*Api[j];
		}
		rn=0;
		for (j=0;j<N*N;j++)
			rn+=(rr[j]*rr[j]+ri[j]*ri[j]);
        if (rn0 != 0) bn=rn/rn0;
        else bn =0;
		for (j=0;j<N*N;j++){
			pr[j]=rr[j]+bn*pr[j];
			pi[j]=ri[j]+bn*pi[j];
		}
		}
free(rr);
free(ri);
free(pr);
free(pi);
free(tr);
free(ti);
free(Apr);
free(Api);

}

void PPFFT( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
                   int N,
				   int preval)
{  
  /* variable declaration */
  // int r,j,HalfN,M;
  int r,j,M;

  REAL alpha,mult;
  REAL sr,si, *tr,*ti,*ur,*ui,*zr,*zi;
  /*REAL mularray[128];*/
  
	// HalfN = N>>1;
	M     = N<<1;
	


	tr  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ti  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ur  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ui  = (REAL *) calloc((int) N*4, sizeof(REAL));
	zr  = (REAL *) calloc((int) M*N, sizeof(REAL));
	zi  = (REAL *) calloc((int) M*N, sizeof(REAL));  
	
	if (tr == NULL || ti == NULL ||
		ur == NULL || ui == NULL ||
		zr == NULL || zi == NULL) {
		printf("PPFFT.c: could not allocate memory for tmp. \n");
		exit(1);
	}

/*	mularray[0] =  (REAL) preval;
	PrintArray(mularray, (REAL *) NULL, 1, 1);*/
/*	PrintArray(xr,xi,N,N); */
	

	
	
	FFTmid0(zr,zi,xr,xi,N,BYCOL);  
/*	PrintArray(zr,zi,M,N); */

	for(r=0; r < N; r++){
	
		alpha=((REAL)(N-r))/N;
		mult = preval ? sqrt(alpha) : 1;
		/*mularray[r] = mult; */
		
		/* t = mult .* FractionalFFT_mid0(z(r,:),alpha); */
		for( j=0; j<N; j++){
			ur[j] = ZR(r,j);
			ui[j] = ZI(r,j);
		}
		FFFT(tr,ti,ur,ui,N,alpha);
		
		/* 	Y(1:n,r) = t.';  */
		for( j=0; j<N; j++){
			YR(j,r) = tr[j]*mult;
			YI(j,r) = ti[j]*mult;
		}
		/* PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */
		
		alpha=((-1) * ((REAL) r))/((REAL) N);

		switch(preval){
			case 0:
				mult =1;
				break;
			case 1:
				mult = r!=0 ? sqrt(fabs(alpha)) : 1/(sqrt((REAL) 2*N));
				break;
			case 2:
				mult = r!=0 ? sqrt(fabs(alpha)) : 1/(2* sqrt((REAL) N));
				break;
		}
		/*mularray[N+r] = mult; */
		/* t = mult .* FractionalFFT_mid0(z(r+n,:),alpha); */
		for( j=0; j<N; j++){
			ur[j] = ZR(r+N,j);
			ui[j] = ZI(r+N,j);
		}
		FFFT(tr,ti,ur,ui,N,alpha);
	/*	PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */
		
		/* Y(1:n,r+n) = t.'; */
		for( j=0; j<N; j++){
			YR(j,r+N) = tr[j]*mult;
			YI(j,r+N) = ti[j]*mult;
		}
	}
	/*PrintArray(mularray,(REAL *) NULL,1,M);*/
 
/*  z = (fft_mid0([zeros(n,n2),X,zeros(n,n2)].'));    % FFT of padded rows */


FFTmid0(zr,zi,xr,xi,N,BYROW);    

/* z=[z(2:m,:);z(1,:)]; %Move left edge to right edge */
for(j=0; j < N; j++){
	sr = ZR(0,j); si = ZI(0,j);
	for( r=1; r<M; r++){
		ZR(r-1,j) = ZR(r,j);
		ZI(r-1,j) = ZI(r,j);
	}
	ZR(M-1,j) = sr;
	ZI(M-1,j) = si;
}
/* PrintArray(zr,zi,M,N); */
for(r=0; r < N; r++){

	alpha=((REAL)(N-r))/N;
	mult = preval ? sqrt(alpha) : 1;

	/* t = mult .* FractionalFFT_mid0(z(m+1-r,:),alpha); */
    for( j=0; j<N; j++){
		ur[j] = ZR(M-1-r,j);
		ui[j] = ZI(M-1-r,j);
	}
	FFFT(tr,ti,ur,ui,N,alpha);

	/* 		 Y(n+1:m,r) = t.';  */
    for( j=0; j<N; j++){
		YR(N+j,r) = tr[j]*mult;
		YI(N+j,r) = ti[j]*mult;
	}

	/*	PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */

	alpha=((-1) * (REAL) r)/N;

	switch(preval){
		case 0:
			mult = 1;
			break;
		case 1:
			mult = r!=0 ? sqrt(fabs(alpha)) : 1/(sqrt((REAL) 2*N));
			break;
		case 2:
			mult = r!=0 ? sqrt(fabs(alpha)) : 1/(2* sqrt((REAL) N));
			break;
	}
		 
	/* t = mult .* FractionalFFT_mid0(z(n+1-r,:),alpha); */
    for( j=0; j<N; j++){
		ur[j] = ZR(N-1-r,j);
		ui[j] = ZI(N-1-r,j);
	}
	FFFT(tr,ti,ur,ui,N,alpha);

	/*	PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */
	
	/* Y(n+1:m,n+r)  = t.' */
    for( j=0; j<N; j++){
		YR(N+j,r+N) = tr[j]*mult;
		YI(N+j,r+N) = ti[j]*mult;
	}
 }
 
 /* Y = Y ./ sqrt(2) n */
 
/*mult = 1.0 / (sqrt((REAL) 2)*N ); 
for(j=0; j < M*M; j++){
		yr[j] = yr[j]*mult;
		yi[j] = yi[j]*mult;
}
*/

free(tr); 
free(ti);
free(ur);
free(ui);
free(zr);
free(zi);


	
}

void IFFTmid0( REAL xr[],REAL xi[],
                    REAL yr[],REAL yi[],
                   int M,
                   int N,
				   int rowcol)
{  
  /* variable declaration */
  // int i,j,HalfN;
  int i,j;
  REAL *tmp,*zr,*zi;
  
  tmp = (REAL *) calloc(N*4, sizeof(REAL));
  zr  = (REAL *) calloc(N*M, sizeof(REAL));
  zi  = (REAL *) calloc(N*M, sizeof(REAL));
  if (tmp == NULL || zr == NULL || zi == NULL) {
	printf("IFFTmid0.c: could not allocate memory for tmp. \n");
		exit(1);
		}
 

	// HalfN = N>>1;
	copy_sub(zr,zi,yr,yi,N*M);
	fftshift(zr,zi,M,N);
	ifftcol(zr,zi,tmp,M,N);
	fftshift(zr,zi,M,N);  

	/* PrintArray(zr,zi,M,N); */
	
	switch(rowcol){
		case BYCOL:
			/* build X = Z matrix */
			for( j= 0; j < N; j++){
				for( i=0; i< M; i++){
					xr[(i) + M*(j)] = zr[(i) + M*(j)];
					xi[(i) + M*(j)] = zi[(i) + M*(j)];
				}
			}
			/* PrintArray(xr,xi,M,N); */
			break;
		case BYROW:
			/* build Y = X.' matrix */
			for( j= 0; j < N; j++){
				for( i=0; i< M; i++){
					xr[(j) + N*(i)] = zr[(i) + M*(j)];
					xi[(j) + N*(i)] = zi[(i) + M*(j)];
				}
			}
			/* PrintArray(xr,xi,N,M); */
	}
	free(tmp);
	free(zr);
	free(zi);
	

}


void copy_sub( REAL yr[],REAL yi[],
               REAL xr[],REAL xi[], 
               int N)
{  
  int j;
  

    for(j=0;j<N;j++){
      yr[j]  = xr[j]; 
      yi[j]  = xi[j]; 
    } 
}


void FFFT( REAL zr[],REAL zi[],
                        REAL xr[],REAL xi[],
                        int N,
                        const REAL alpha)
{  
  /* variable declaration */
  int j,TwoN;
  REAL phase,temp;
  REAL *seq1, *seq2;
  REAL E,re,im,Ninv;
  
  seq1 = (REAL *) calloc(N*4, sizeof(REAL));
  seq2 = (REAL *) calloc(N*4, sizeof(REAL));
  if (seq1 == NULL) {
    
	  printf("FFFT.c: could not allocate memory for seq1.. \n");
		exit(1);
		}
  if (seq2 == NULL){
		printf("FFFT.c: could not allocate memory for seq2.. \n");
		exit(1);
		}
    
	E = TWOPI/2 * alpha;
	Ninv = 1.0/((REAL) N);
	TwoN = N<<1;
     /* 1. point multiplication */

    for(j=0;j<N;j++){
      phase  = E * (j - j*j*Ninv);
      Re(seq1,j)   = xr[j]*cos(phase) - xi[j]*sin(phase);
      Im(seq1,j)   = xr[j]*sin(phase) + xi[j]*cos(phase); 
      
      Re(seq1, (N+j)) = 0;
      
      Im(seq1, (N+j)) = 0;
    }
    /*  (b)   second sequence  */
    Re(seq2,0) = 1; 
    Im(seq2,0) = 0;
    Re(seq2,N) = 0;
    Im(seq2,N) = 0;
    for(j=1;j<N;j++){
      phase = E/N * j*j;
      Re(seq2,j) = cos(phase); 
      Im(seq2,j) = sin(phase); 
      Re(seq2,(TwoN-j)) = cos(phase); 
      Im(seq2,(TwoN-j)) = sin(phase); 
    }
    /*  2.   convolve via FFT */
    four1(seq1-1,TwoN,-1);
    four1(seq2-1,TwoN,-1);
    for(j=0;j<TwoN;j++){
      temp       = Re(seq1,j)*Re(seq2,j) - Im(seq1,j)*Im(seq2,j); 
      Im(seq1,j) = Re(seq1,j)*Im(seq2,j) + Re(seq2,j)*Im(seq1,j); 
      Re(seq1,j) = temp; 
    }
    four1(seq1-1,TwoN,1);

    /* 3. point multiplication */
    for(j=0;j<N;j++){
      phase  = E*(j - ((REAL) N)/2 - j*(j*Ninv));
      re     = Re(seq1,j)/TwoN;
      im     = Im(seq1,j)/TwoN; 
      zr[j]  = re*cos(phase) - im*sin(phase); 
      zi[j]  = re*sin(phase) + im*cos(phase); 
    }
	free(seq1);
	free(seq2);
}

void FFTmid0( REAL yr[],REAL yi[],
                   REAL xr[],REAL xi[], 
				   int N,
				   int rowcol)
{  
  /* variable declaration */
  int i,j,HalfN,M;
  REAL *tmp;
  
  tmp = (REAL *) calloc(N*4, sizeof(REAL));

  
  if (tmp == NULL) {
	  printf("FFFT.c: could not allocate memory for tmp. \n");
		exit(1);
		}
    
	  
	HalfN = N>>1;
	M     = N<<1;
	
	switch(rowcol){
		case BYCOL:
			/* build Y = [0 ; X ; 0 ] matrix */
			for( j= 0; j < N; j++){
				for( i=0; i< HalfN; i++)
					YR(i,j) = YI(i,j) = 0;
				for( i=0; i< N; i++){
					YR(i+HalfN,j) = XR(i,j);
					YI(i+HalfN,j) = XI(i,j);
				}
				for( i=0; i< HalfN; i++)
					YR(i+N+HalfN,j) = YI(i+N+HalfN,j) = 0;
			}
			break;
		case BYROW:
			/* build Y = [0   X   0 ].' matrix */
			for( j= 0; j < N; j++){
				for( i=0; i< HalfN; i++)
					YR(i,j) = YI(i,j) = 0;
				for( i=0; i< N; i++){
					YR(i+HalfN,j) = XR(j,i);
					YI(i+HalfN,j) = XI(j,i);
				}
				for( i=0; i< HalfN; i++)
					YR(i+N+HalfN,j) = YI(i+N+HalfN,j) = 0;
			}
	}
	fftshift(yr,yi,M,N);
	fftcol(yr,yi,tmp,M,N);
	fftshift(yr,yi,M,N);
	free(tmp);  
}


void fftshift(REAL zr[],REAL zi[],
              int M,
			  int N)
{  
  /* variable declaration */
  REAL tempr; 
  int i,j,HalfN,HalfM; 

	HalfN = N>>1;
	HalfM = M>>1;
  /* fftshift columns */
  for(j=0;j<HalfN;j++){
   for(i=0;i<M;i++){
      SWAP(ZR(i,j), ZR(i,j+HalfN));
      SWAP(ZI(i,j), ZI(i,j+HalfN));
    }
  }
  /* fftshift rows */
  for(i=0;i<HalfM;i++){
   for(j=0;j<N;j++){
      SWAP(ZR(i,j), ZR(i+HalfM,j));
      SWAP(ZI(i,j), ZI(i+HalfM,j));
    }
  }
}


void four1(REAL data[], int nn, int isign)
{
        int n,mmax,m,j,istep,i;
        REAL wtemp,wr,wpr,wpi,wi,theta;
        REAL tempr,tempi;

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

void ifftcol( REAL yr[], REAL yi[], REAL tmp[], 
             int M,
             int N)
{  
  /* variable declaration */
  int i,j;

	for( j=0; j<N; j++){

		 /* 1. Interleave */	
		for(i=0;i<M;i++){
		  Re(tmp,i)   = YR(i,j);
		  Im(tmp,i)   = YI(i,j); 
		}
	   
		/*  2. Call four1 */
		four1(tmp-1,M,1);
		
		/*  3. De-Interleave */
		for(i=0;i<M;i++){
		  YR(i,j)  = Re(tmp,i)/M; 
		  YI(i,j)  = Im(tmp,i)/M; 
		}
	}
}

void fftcol( REAL yr[], REAL yi[], REAL tmp[], 
             int M,
             int N)
{  
  /* variable declaration */
  int i,j;

	for( j=0; j<N; j++){

		 /* 1. Interleave */	
		for(i=0;i<M;i++){
		  Re(tmp,i)   = YR(i,j);
		  Im(tmp,i)   = YI(i,j); 
		}
	   
		/*  2. Call four1 */
		four1(tmp-1,M,-1);
		
		/*  3. De-Interleave */
		for(i=0;i<M;i++){
		  YR(i,j)  = Re(tmp,i); 
		  YI(i,j)  = Im(tmp,i); 
		}
	}
}

void Adj_PPFFT( REAL xr[],REAL xi[],
                   REAL yr[],REAL yi[], 
                   int M,
				   int preval)
{  
  /* variable declaration */
  int r,j,HalfN,N;
  REAL alpha,mult,renorm;
  REAL sr,si, *tr,*ti,*ur,*ui;
  REAL *x0r,*x0i,*x1r,*x1i;
  REAL *z0r,*z0i,*z1r,*z1i;
  /*REAL mularray[128];*/
	N     = M>>1;  
	HalfN = N>>1;
	
	tr  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ti  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ur  = (REAL *) calloc((int) N*4, sizeof(REAL));
	ui  = (REAL *) calloc((int) N*4, sizeof(REAL));
	x0r  = (REAL *) calloc((int) M*N, sizeof(REAL));
	x0i  = (REAL *) calloc((int) M*N, sizeof(REAL));  
	x1r  = (REAL *) calloc((int) M*N, sizeof(REAL));
	x1i  = (REAL *) calloc((int) M*N, sizeof(REAL));  
	z0r  = (REAL *) calloc((int) M*N, sizeof(REAL));
	z0i  = (REAL *) calloc((int) M*N, sizeof(REAL));  
	z1r  = (REAL *) calloc((int) M*N, sizeof(REAL));
	z1i  = (REAL *) calloc((int) M*N, sizeof(REAL));  
	if (tr == NULL || ti == NULL ||
		ur == NULL || ui == NULL ||
		z1r== NULL || z1i == NULL||
		z0r== NULL || z0i == NULL||
		x0r== NULL || x0i == NULL||
		x1r== NULL || x1i == NULL) {
		printf("Adj_PPFFT.c: could not allocate memory for tmp. \n");
		exit(1);
	}
	

 /* Y = Y ./ sqrt(2) n */
 
renorm = (sqrt((REAL) 2)*N ); 

	/*	mularray[0] =  (REAL) preval;
	PrintArray(mularray, (REAL *) NULL, 1, 1);*/

	for(r=0; r < N; r++){
	
		alpha=((REAL)(N-r))/N;
		mult = preval ? sqrt(alpha) : 1;
		/*mularray[r] = mult; */
		
		/* t = mult .* FractionalFFT_mid0(Y0(:,r).',-alpha); */
		for( j=0; j<N; j++){
			ur[j] = Y0R(j,r)*renorm;
			ui[j] = Y0I(j,r)*renorm;
		}
		FFFT(tr,ti,ur,ui,N,-alpha);
		
		/* 	X0(r,:) = t.';  */
		for( j=0; j<N; j++){
			X0R(r,j) = tr[j]*mult;
			X0I(r,j) = ti[j]*mult;
		}
		/* PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */ 
		
		alpha=((-1) * ((REAL) r))/((REAL) N);

		switch(preval){
			case 0:
				mult =1;
				break;
			case 1:
				mult = r!=0 ? sqrt(fabs(alpha)) : 1/(sqrt((REAL) 2*N));
				break;
			case 2:
				mult = r!=0 ? sqrt(fabs(alpha)) : 1/(2* sqrt((REAL) N));
				break;
		}
		/*mularray[N+r] = mult; */
		/* t = mult .* FractionalFFT_mid0(Y0(:,r+n).',-alpha); */
		for( j=0; j<N; j++){
			ur[j] = Y0R(j,r+N)*renorm;
			ui[j] = Y0I(j,r+N)*renorm;
		}
		FFFT(tr,ti,ur,ui,N,-alpha);
		/* PrintArray(ur,ui,1,N);
		PrintArray(tr,ti,1,N); */ 
		
		/* X0(r+n,:) = t.'; */
		for( j=0; j<N; j++){
			X0R(r+N,j) = tr[j]*mult;
			X0I(r+N,j) = ti[j]*mult;
		}
	}
	/* PrintArray(x0r,x0i,M,N); */


/* X0 = ifft_mid0(X0); */
IFFTmid0(z0r,z0i,x0r,x0i,M,N,BYCOL);    
/* PrintArray(z0r,z0i,M,N); */

for(r=0; r < N; r++){

	alpha=((REAL)(N-r))/N;
	mult = preval ? sqrt(alpha) : 1;

	/* t = mult .* FractionalFFT_mid0(Y1(:,r).',-alpha); */
    for( j=0; j<N; j++){
		ur[j] = Y1R(j,r)*renorm;
		ui[j] = Y1I(j,r)*renorm;
	}
	FFFT(tr,ti,ur,ui,N,-alpha);

	/* 		 X1(m+1-r,:) = t;  */
    for( j=0; j<N; j++){
		X1R(M-1-r,j) = tr[j]*mult;
		X1I(M-1-r,j) = ti[j]*mult;
	}

	/* PrintArray(ur,ui,1,N);
	PrintArray(tr,ti,1,N); */ 

	alpha=((-1) * (REAL) r)/N;

	switch(preval){
		case 0:
			mult = 1;
			break;
		case 1:
			mult = r!=0 ? sqrt(fabs(alpha)) : 1/(sqrt((REAL) 2*N));
			break;
		case 2:
			mult = r!=0 ? sqrt(fabs(alpha)) : 1/(2* sqrt((REAL) N));
			break;
	}
		 
	/* t = mult .* FractionalFFT_mid0(Y1(:,n+r).',:),-alpha); */
    for( j=0; j<N; j++){
		ur[j] = Y1R(j,r+N)*renorm;
		ui[j] = Y1I(j,r+N)*renorm;
	}
	FFFT(tr,ti,ur,ui,N,-alpha);

	/* PrintArray(ur,ui,1,N);
	PrintArray(tr,ti,1,N);  */
	
	/* X1(n+1-r,:)  = t */
    for( j=0; j<N; j++){
		X1R(N-1-r,j) = tr[j]*mult;
		X1I(N-1-r,j) = ti[j]*mult;
	}
 }

/* X1=[X1(m,:) X1(1:m-1,:)]; %Move left edge to right edge */
for(j=0; j < N; j++){
	sr = X1R(M-1,j); si = X1I(M-1,j);
	for( r=M-1; r>0; --r){
		X1R(r,j) = X1R(r-1,j);
		X1I(r,j) = X1I(r-1,j);
	}
	X1R(0,j) = sr;
	X1I(0,j) = si;
}
/* PrintArray(x1r,x1i,M,N); */
/* PrintArray(zr,zi,M,N); */
/*  X1 =(ifft_mid0(X1)).';    % FFT of padded rows */
IFFTmid0(z1r,z1i,x1r,x1i,M,N,BYROW);    
/* PrintArray(z1r,z1i,N,M); */

/* Combine */

for( j=0; j < N; j++){
	for( r=0; r < N; r++){
		XR(r,j) = (Z0R(r,j)+Z1R(r,j))/(sqrt(2.0)*N); 
		XI(r,j) = (Z0I(r,j)+Z1I(r,j))/(sqrt(2.0)*N);
	}
}

free(tr); 
free(ti);
free(ur);
free(ui);
free(x0r);
free(x0i);
free(x1r);
free(x1i);
free(z0r);
free(z0i);
free(z1r);
free(z1i);
} 

void FFTmid0NP(REAL xr[],REAL xi[],
                   int M,
                   int N)
{  
  /* variable declaration */
  REAL *tmp;
  tmp = (REAL *) calloc(N*4, sizeof(REAL));
 
  if (tmp == NULL) {
	printf("IFFTmid0.c: could not allocate memory for tmp. \n");
		exit(1);
		}
 

	fftshift(xr,xi,M,N);
	fftcol(xr,xi,tmp,M,N);
	fftshift(xr,xi,M,N);  
    free(tmp);
}


                
}
