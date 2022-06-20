/*******************************************************
Program: 'lognormal.cc'
Author: P. Arnalte-Mur (OAUV)
Creation date: 04/06/2009
Last modification date: 09/09/2009

Program to generate a random log-Gaussian field inside a 
cubic box.

This is a modification of 'simulgausbox.cc' (see comments 
below), just making the last transformation from Gaussian
to log-Gaussian field.

The input power spectrum  **corresponds to the 'initial'
Gaussian field** (should be calculated somewhere else 
from the proper one). We also need as input the alpha^2
parameter (for normalization of the log-gaussian field),
and the sigma_0/alpha ratio (for relative normalization
of both fields, depends on how you calculated the P(k)
for the gaussian field).

Also, substitute routines and types from 'genus' by
my own to make it more portable.

Version history:

  V. 0.1 (04/06/2009): Initial version.

  V. 0.2 (04/06/2009): Use linear interpolation instead of
         spline: this will avoid problems with possible
         negative values of P(k).

  V. 0.3 (08/06/2009): Implement linear interpolation directly,
         to avoid need for external GSL routines.

  V. 0.4 (09/06/2009): Add initialization from computer time
         if there is not random seed provided!

  V. 0.5 (08/09/2009): Change function gen_dens() (coming from
         'genus' software) to know step by step what I'm doing.
         I will use FFT from library FFTW-3 (external).
         Use FFTW routines specific for real fields.
********************************************************/

/****************************************************
Comments for 'lognormal.cc' program:
-----------------------------------------------------
Program to generate a random Gaussian field inside a
cubic box, with a power spectrum given from an external
file.
 
This program was created from 'simgauss.cc' ( from the 'genus'
software package) written by Enn Saar and Jean-Luc Starck. 
Changes made to it:
allow for user-defined input powspect, do not smooth things.

P. Arnalte-Mur (OAUV)
Date created: 23/04/2009
Date last modified: 12/05/2009

Version 0.1 (23/04/2009):  
   Things needed to improve: correct scaling/units: for
      now, we assume that spacing between points in the
      grid equals 1 unit of distance. So, if using 
      'standard' input power spectrum, grid points are
      separated a distance of 1 Mpc/h

Version 0.2 (24/4/2009):
    Should now have the correct spacing (changed in the
      calculation of 'kphys' inside 'gens_dens' function).

Version 0.3 (27/04/2009):
    Scaling was incorrect in v. 0.2, need to also add a 
     re-scaling of the Fourier-space density, following
     eq. (12.1.8) in Numerical Recipes.

    Add correct normalization for the (inverse) Fourier
     transform: there was a N^-3 factor missing.

    We add these corrections outside other routines, so 
     they are not specially efficient (could improve this).

Version 0.4 (12/05/2009):
    Change scaling again. Now it should be consistent thorough
     with the following conventions:
    
     * For continuous FT:

       \hat{\delta}(\vec{k}) = 
           \int d^3x e^{i \vec{k}\dot\vec{x}} \delta(\vec{x})

       \delta(\vec{x}) = (2\pi)^{-3} 
           \int d^3k \hat{\delta}(\vect{k}) e^{-i \vec{k}\dot\vec{x}}

     * For discrete FT:

       \delta_{\vec{k}} = \sum_{\vec{j}} e^{i \vec{k}\dot\vec{x}_{\vec{j}}} 
                          \delta(\vec{x}_{\vec{j}})

       \delta(\vec{x}_{\vec{j}}) = N^{-3} \sum_{\vec{k}} 
           e^{-i \vec{k}\dot\vec{x}_{\vec{j}}} \delta_{\vec{k}}


Version 0.5 (12/05/2009):
    Change output to be written in ASCII format (had problems
     with binary FITS in 64bit machines).

Antoine Labatie
 Version 0.6 (05/02/2010):
 Change again output to write .fits, linking external .h to do this
 
 Version 0.7 (18/02/2010):
 Changed normalization of the lognormal field by adding each individual
 variance of independent modes.  
 Possible improvement for computing the mode variance (integrate P(k) on each
 cell volume [kx-dk/2,kx+dk/2]*[ky-dk/2,ky+dk/2]*[kz-dk/2,kz+dk/2] ?? Because the
 by taking the value in the middle of the cell the total power is underestimated
 because the function is concave.
 
*******************************************************/

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include "Array.h"
#include "IM_IO.h"
#include "lognormal.h"
#include <omp.h>

//For FFTW library
#include <fftw3.h>

int N=32;                  /* grid size, power of 2*/
int N2;                 /* N/2                  */
int N21;
int NCB;                /* N^3                  */

double *dens;              /* final density   */
FILE *out;				/* density file			*/
//FILE *outk;
//FILE *outg;

void gauss(double disp, double *x, double *y);
double spect(double k);


double Delta;        //Grid step (in Mpc/h): spacing between points in the grid
double DCB;          //Delta^3


//Added for log-normal field and different normalizations
double alpha2 = 1.;    //Normalization of log-gauss field
double saratio = 1.;   //Relative normalization of gaus and log-gaus fields


bool Verbose = false;
long seed;

char Name_Out_Cat_Suffix[256]; /* galaxy catalogue output file name */
char Name_Out_Cat_Prefix[256];
char Name_Out_Cat[256];

bool Catalogue_Random=false;


//Power spectrum tab
FILE *inpk;
char Name_Pk_In[256];
double *ktab, *pktab;
long npk=0;

//nbar tab
FILE *indens;
double *rtab, *denstab;
long ndens=0;


double globalvariance=0.0;



//maximum number of procs used for the loops
int Nproc_max=40;

//parameters to be set in lognormal.param
bool Use_Density=false;
char Name_Density[256];
bool Use_Mask=false;
char Name_Mask[256];
double L;       //Lenght of box side, in Mpc/h.
double x2min, y2min, z2min;
double nbar;


/***************************************************************************/
 
static void usage(char *argv[])
{

  fprintf(stderr, "Usage: %s options in_Pk_file out_Catalogue_file_suffix\n\n", argv[0]);
  fprintf(stderr, "   where options = \n");
	
  fprintf(stderr, "        [-d Dimension]\n");
  fprintf(stderr, "            Cube size dimension in the simulation.\n");
  fprintf(stderr, "            Default is %2d.\n\n", N);

  fprintf(stderr, "         [-s SigAlpRatio]\n");
  fprintf(stderr, "             Ratio of normalizations of gaussian to log-gaussian fields: sigma_0/alpha.\n");
  fprintf(stderr, "             Default is %.1f\n\n", saratio);

	
  fprintf(stderr, "         [-I InitRandomVal]\n");
  fprintf(stderr, "             Value used for random value generator initialization.\n\n");

  fprintf(stderr, "         [-r]\n");
  fprintf(stderr, "             Generate random catalogue (i.e. a catalogue with no fluctuations).\n\n");
	
  fprintf(stderr, "         [-v]\n");
  fprintf(stderr, "             Verbose.\n\n");
  

  fprintf(stderr, "\n");
  exit(-1);
}


void get_args(int argc, char *argv[])
{

  
  /* Require arguments (need at least one) !! */
  if(argc == 1){
    usage(argv);
  }
  /* Start at i = 1 to skip the command name. */
  int i=1;
  
    /* Check for a switch (leading "-"). */

  while(argv[i][0] == '-') {

      /* Use the next character to decide what to do. */

    switch (argv[i][1]) {

    case 'd': N = atoi(argv[++i]);
      break;

    case 's': saratio = atof(argv[++i]);
      break;

    case 'I': 
      seed = atol(argv[++i]);
      break;

	case 'r': Catalogue_Random = true;
		break;
			
    case 'v': Verbose = true;
      break;
			
    case '?': usage(argv);
      break;
      
    default:  usage(argv);
      break;
    }
    i++;
	if(i==argc) usage(argv);
  }



  if(i<(argc)-1){
      strcpy(Name_Pk_In, argv[i++]);
      strcpy(Name_Out_Cat_Suffix, argv[i++]);
  }
  else usage(argv);

  if(i < argc){
    fprintf(stderr, "Too many parameters: %s ...\n", argv[i]);
    usage(argv);
  }

}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void siminit(int argc, char *argv[])
{
    seed = time(NULL);

    get_args(argc,argv);

    srand48((long) seed);

    N2=N/2; 
    N21 = N2+1;
    NCB=N*N*N; 

    //Dimensions N^3
    dens = (double *) malloc(NCB*sizeof(double)); assert(dens);

    Delta = L/((double)N);
    DCB=Delta*Delta*Delta;

	//Read pk
    int i;
	int ret;
    double temp_k, temp_pk;
    inpk = fopen(Name_Pk_In, "r"); assert(inpk);

	while ( fscanf(inpk, "%lf %lf\n", &temp_k,&temp_pk)!=EOF ) 	npk++;
	
	inpk = fopen(Name_Pk_In, "r"); 
	ktab = (double *)malloc(npk*sizeof(double)); assert(ktab);
    pktab = (double *)malloc(npk*sizeof(double)); assert(pktab);

    for(i=0;i<npk;i++){
      ret=fscanf(inpk, "%lf %lf\n", &temp_k,&temp_pk);
      ktab[i]=temp_k; pktab[i]=temp_pk;
    }
   fclose(inpk);
	
	//Read nbar tab if external ascii file provided (i.e. if Use_Density==true in lognormal.param)
	if(Use_Density)
	{
		double temp_r, temp_dens;
		indens = fopen(Name_Density, "r"); assert(indens);
	
		while ( fscanf(indens, "%lf %lf\n", &temp_r,&temp_dens)!=EOF ) 	ndens++;	
	
		indens = fopen(Name_Density, "r"); 
		rtab = (double *)malloc(ndens*sizeof(double)); assert(rtab);
		denstab = (double *)malloc(ndens*sizeof(double)); assert(denstab);
	
		for(i=0;i<ndens;i++){
			ret=fscanf(indens, "%lf %lf\n", &temp_r,&temp_dens);
			rtab[i]=temp_r; denstab[i]=temp_dens;
		}
		fclose(indens);
	}
}

/*********************************************************************/

/* GET PARAMETERS */

void get_param()
{
	char Name_Param_File[256];
	sprintf(Name_Param_File, "../param/lognormal.param");	
	FILE *File=fopen(Name_Param_File,"r");
	
    if (File == NULL)
    {
		cerr << "Error: cannot open file "  <<  Name_Param_File << endl;
		exit(-1);
    }
	char Temp[256];
	int temp_bool;
	int ret;
	ret=fscanf(File, "%s\t%i\n", Temp, &temp_bool);	//Use external ascii file for density as a function of distance (0 means false otherwise means true)
	Use_Density = (temp_bool !=0);
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Density);	//Name of ascii file for density
	ret=fscanf(File, "%s\t%i\n", Temp, &temp_bool);		//Use external fits file for mask (0 means false otherwise means true)
	Use_Mask = (temp_bool !=0);
	ret=fscanf(File, "%s\t%s\n\n", Temp, Name_Mask);	//Name of fits file for mask

	//minimum value of x2,y2 and z2 axis. (x2,y2,z2) basis is adapted to geometry of SDSS.
	//when providing an external (lambda,eta) mask, these coordinates will be converted to regular (x,y,z) coordinates by the transform:
	//x = (x2-y2)/sqrt(2.0); 		y = (x2+y2)/sqrt(2.0);   			z = z2;
	//(x,y,z) are linked to (eta,lambda) by the transform:
	//  x = -D*sin(lambda) ;	y = D*cos(lambda)*cos(eta)  ;  z = D*cos(lambda)*sin(eta)
	//when no mask is provided, the output coordinates stay in (x2,y2,z2) basis
	ret=fscanf(File, "%s\t%lf\n", Temp, &x2min);		
	ret=fscanf(File, "%s\t%lf\n", Temp, &y2min);
	ret=fscanf(File, "%s\t%lf\n", Temp, &z2min);	
	
	ret=fscanf(File, "%s\t%lf\n", Temp, &L);	//size of cubic box in Mpc/h
	ret=fscanf(File, "%s\t%lf\n\n", Temp, &nbar); //constant mean density of points in the volume (use only when Use_Density==false in lognormal.param)

	ret=fscanf(File, "%s\t%s\n", Temp, Name_Out_Cat_Prefix);
	
	fclose(File);
	
}


/*********************************************************************/


void gen_dens()
{
  int ix,iy,iz, i;
  int kxind, kyind, kzind;
  double kmod, variance, a, b;
  double Deltak;

  int KDIM = N*N*N21;

  //To avoid errors due to int division
  double NCBf = NCB;
	
  fftw_complex *densk;
  double *densr;
  fftw_plan plan;

  densr = (double *) malloc(sizeof(double) * NCB);
  densk = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * KDIM);
	
  /*******************************/
  printf("Size of a double = %ld\n", sizeof(double));
  printf("Size of a fftw_complex = %ld\n", sizeof(fftw_complex));
  printf("Size of 'densk' = %ld\n", sizeof(densk));
  printf("Size of 'densr' = %ld\n", sizeof(densr));
  /********************************/

	
  //Use FFTW's 'complex2real' routine, but now the sign in the exponent is reversed wrt
  //our convention!!!!!
  //plan = fftw_plan_dft_c2r_3d(N, N, N, densk, densr, FFTW_ESTIMATE);
  plan = fftw_plan_dft_c2r_3d(N, N, N, densk, densr, FFTW_MEASURE);


  Deltak = 2.*M_PI/L;                   //Spacing between nodes in k grid
  fprintf(stderr, "\nDelta_k = %f\n\n",Deltak);

	
  //Generate k-modes as Gaussian distributed real and imaginary parts.
  //Will have, as independent k-modes 0<=ix<N; 0<=iy<N; 0<=iz<N21
	#ifdef _OPENMP
		int Nproc=omp_get_num_procs();
		if(Nproc>Nproc_max) Nproc=Nproc_max;
		printf("\n %u processors used for the loop \n\n",Nproc);
	#endif	
	
	#pragma omp parallel default(none)  shared(N,N2,N21,NCB,NCBf,DCB,densk,Deltak,stderr,globalvariance) \
	private(ix,iy,iz,kxind,kyind,kzind,kmod,variance,a,b) num_threads(Nproc)
	{
		#pragma omp for schedule(static)
		for(ix=0;ix<N;ix++)
		{
			if(ix<=N2) {kxind = ix;}   //Positive kx
			else {kxind = ix - N;}     //Negative kx
			for(iy=0;iy<N;iy++)
			{
				if(iy<=N2) {kyind = iy;}   //Positive ky
				else {kyind = iy - N;}     //Negative ky
				for(iz=0;iz<N21;iz++) 
				{
					kzind = iz;               //Only positive kz
					kmod = Deltak*sqrt( (kxind*kxind) + (kyind*kyind) + (kzind*kzind) );

					/*
					 Not all the modes generated in densk are independent normally. 
					 Indeed the modes (kxind,kyind,kzind) and (N-kxind,N-kyind,N-kzind) are always complex conjugates.
					 So the part kz>=N21 can always be deduced from the part kz<N21 and the missing part of the array is added following this.
					 But other parts are also redundant for example the modes (0,kyind,kzind) and (0,N-kyind,N-kzind).
					 fftw does a normal inverse FFT and finally keeps the real part only so we can still generate independent modes even if 
					 they are not (we just have to adjust the variance to have the correct amount of power for each mode).
				 
					 Each mode k verifies <|delta_k|^2|>=<a^2+b^2>=P(k)(dk)^3/(2*Pi)^3=P(k)/L^3
					 Because of the further normalization by N^3 we must multiply this by N^6 to get <a^2>=<b^2>=0.5*P(k)*N^6/L^3=0.5*P(k)*NCB/DCB
					 For the modes that we generate independently but which are not in reality independent (there is a pair which should be conjugate) 
					 the result will result in adding 2 indepdent variables X+Y instead of having 2*X so the variance will be underestimated by a factor 2,
					 that's why we take for them P(k)*NCB/DCB.
					 Finally for the modes that are always real like P(0) and P(N/2) if N is even then <|delta_k|^2|>=<a^2> and the variance of a must 
					 be also taken =P(k)*NCB/DCB
					 */
				
					if(iz==0 || (N/2==N/2.0 && iz ==N/2)) 
					{
						variance = spect(kmod)*NCB/DCB;
						# pragma omp critical
							globalvariance+= variance;
					}
					else 
					{
						variance = 0.5*spect(kmod)*NCB/DCB;
						# pragma omp critical
							globalvariance+= 4.0*variance;
					}
					gauss(variance, &a, &b); 
					densk[(N*N21*ix)+(N21*iy)+iz][0] = a;     //Real part
					densk[(N*N21*ix)+(N21*iy)+iz][1] = b;     //Imag. part
			
				}
			}
		}
	}

  
	//Make FFT transform
	fftw_execute(plan);
  
	//Now, values of \delta(x) (unnormalized) are stored in densr, we normalize it 
	//(same for 3 dimensions, real and imag. parts):
	for(i=0;i<NCB;i++)
		densr[i] = densr[i]/NCBf;

  globalvariance/=(NCBf*NCBf);
	
  //Calculate value of \alpha2 from the global variance
  alpha2=globalvariance;
  
  printf("\nGaussian field:\n");
  printf("Theoretical global variance = %f\n", globalvariance);

	
  /*************************************************/

  //Transformation to log-gaussian field, just local transform
    //(could do it with less loops, I know...)
  for(ix=0; ix<N; ix++){
    for(iy=0;iy<N;iy++){
      for(iz=0;iz<N;iz++){
	dens[(N*N*ix) + (N*iy) + iz] = exp( (saratio*densr[(N*N*ix) + (N*iy) + iz]) - (alpha2/2.) );
      }
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(densk);
  free(densr);
  
}



/*************************************************/

void gauss(double disp, double *x, double *y)    /* Recipes */
{
	//Create a 2D gaussian independent in x,y  => f(x,y)=1/(2*pi*sigma) e^-(x^2+y^2)/(2*sigma^2)
    double v1,v2,r,fac;

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
  

double spect(double k)
{

  //Change: assign P(0) = 0 -- Is this correct????
  if(k==0){
    return(0.);
  }
  else{

    //Implement here the linear interpolation:
    int i=0;
	  while(k > ktab[i] && i<npk-1)  i++;
    if(i==0 || i>=(npk-1)){
	  printf("ktab[0]: %f, k: %f, i: %u, npk: %lu \n",ktab[0],k,i,npk);
      fprintf(stderr, "Error: k-value outside of tabulated values!\n");
      exit(1);
    }
    
    double P;
    P = pktab[i-1] + ( (pktab[i] - pktab[i-1])*(k-ktab[i-1])/(ktab[i] - ktab[i-1]));

    return(P);
  }

}

int main(int argc, char ** argv)
{
    int i,j,k;
	
	get_param();
    siminit(argc,argv);
	strcpy(Name_Out_Cat, Name_Out_Cat_Prefix); strcat(Name_Out_Cat, Name_Out_Cat_Suffix);
	
    if (Verbose)
    { 
		printf("\n\n# PARAMETERS: \n\n");
		printf("# Input P(k) File = %s\n", Name_Pk_In);
		if(Catalogue_Random) printf("# Generate random catalogue with no fluctuations\n");
		printf("# Write catalogue in ascii = %s\n", Name_Out_Cat);
		printf("# Dimension = %d\n", N);
		printf("# Box Length = %f Mpc/h\n", L);    
		if(Use_Density) 
			printf("Use nbar provided in %s\n",Name_Density);
		else 
			printf("# nbar = %f (Mpc/h)^{-1}\n", nbar);
		if(Use_Mask) printf("Use mask provided in %s\n",Name_Mask);
		printf("# (x2min,y2min,z2min) = %f\t%f\t%f\n\n", x2min,y2min,z2min);
    }
	
	fltarray dens_array(N,N,N);
	double mean=0;
	double sumsq=0;
	double sigma=0;
	
	if(!Catalogue_Random)
	{
		gen_dens();
		
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				for(k=0;k<N;k++){
					mean = mean + dens[N*N*i+N*j+k];
					sumsq = sumsq + (dens[N*N*i+N*j+k]*dens[N*N*i+N*j+k]);
					dens_array(i,j,k)=dens[N*N*i+N*j+k];
				}
			}
		}
	}
	else //Generate random catalogue (no fluctuations)
	{
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				for(k=0;k<N;k++){
					dens_array(i,j,k)=1.0;
				}
			}
		}
	}
	
    if(Verbose){
      mean = mean/NCB;
      sigma = (sumsq/NCB) - (mean*mean); sigma=sqrt(sigma);
      if(!Catalogue_Random)
	  {
		  printf("\nLognormal field:\n");
		  printf("Sample Mean = %f , Sigma = %f\n ", mean, sigma);
	  }
    }
	
	
	//Do Poisson sampling on the density field to get a galaxy catalogue with a given density
	if(Verbose) printf("\nSample continuous and write galaxy catalogue in %s\n",Name_Out_Cat);
	long int Ngal=0;
	int number=-1;
	int n=0;
	
	//Open Out_Catalogue_File
	FILE *outcatalogue;
	outcatalogue = fopen(Name_Out_Cat, "w"); assert(outcatalogue);

	double shift_x,shift_y,shift_z;
	double x2,y2,z2;
	double x,y,z;
	
	double dmin,dmax,deltad;
	int indd;
	
	double d,lambda,eta,dlambda,deta;
	int indlambda,indeta;
	fltarray Mask;
	
	double r2d=180/M_PI;
		
	fprintf(outcatalogue, "x  y  z \n");

	if(Use_Density) 
	{
		deltad=(rtab[ndens-1]-rtab[0])/double(ndens-1.0);
		dmin=rtab[0]-deltad/2.0;
		dmax=rtab[ndens-1]+deltad/2.0;
	}
	if(Use_Mask)
	{
		fits_read_fltarr(Name_Mask,Mask);
		dlambda=360.0/double(Mask.nx());
		deta=180.0/double(Mask.ny());
	}
	

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				x2=i*Delta+x2min;  			y2=j*Delta+y2min;  				z2=k*Delta+z2min;
				
				if(Use_Mask)
				{
					x=(x2-y2)/sqrt(2.0); 		y=(x2+y2)/sqrt(2.0);   			z=z2;
					d=sqrt(x*x+y*y+z*z);
					if(y>0) lambda=-asin(x/d); 
					if(y<=0) lambda=M_PI+asin(x/d);
					eta=asin(z/(d*cos(lambda)));
					lambda*=r2d; 				eta*=r2d;
					
					if(lambda>180.0) 
						lambda-=360.0;
					
					indlambda=floor((lambda+180.0)/dlambda);
					indeta=floor((eta+90.0)/deta);
				
					if( (Mask(indlambda,indeta)==1) && (d>=dmin) && (d<=dmax) )
					{
						if(Use_Density)
						{
							indd=floor((d-dmin)/deltad);
							if(indd<0 || indd>=ndens) 
								nbar=0;
							else
								nbar=denstab[indd];
						}
						n=(int) poidev(nbar*pow(Delta,3.)*dens_array(i,j,k),&number);
						while(n>0)
						{
							shift_x=drand48(); 	shift_y=drand48();  shift_z=drand48();
						
							x2=(i+shift_x)*Delta+x2min; y2=(j+shift_y)*Delta+y2min; z2=(k+shift_z)*Delta+z2min;
							x=(x2-y2)/sqrt(2.0); y=(x2+y2)/sqrt(2.0); z=z2;

		
							
							fprintf(outcatalogue, "%f  %f  %f \n",x,y,z);
							Ngal++;
							n--;
						}
					}
				}
				else
				{
					d=sqrt(x2*x2+y2*y2+z2*z2);
					if(Use_Density)
					{
						indd=floor((d-dmin)/deltad);
						if(indd<0 || indd>=ndens) 
							nbar=0;
						else
							nbar=denstab[indd];
					}
					n=(int) poidev(nbar*pow(Delta,3.)*dens_array(i,j,k),&number);
					while(n>0)
					{
						shift_x=drand48(); 	shift_y=drand48();  shift_z=drand48();
						
						x2=(i+shift_x)*Delta+x2min; y2=(j+shift_y)*Delta+y2min; z2=(k+shift_z)*Delta+z2min;
						
						fprintf(outcatalogue, "%f  %f  %f \n",x2,y2,z2);
						Ngal++;
						n--;
					}
				}
			}
		}	
	}
	
	printf("N_galaxies= %i\n \n \n",Ngal);
	fclose(outcatalogue);
	
	
	free(dens);
    exit(0);
}
