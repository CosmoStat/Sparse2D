/****************************************************
   Program: 'ps_transform.c'
   Author: Pablo Arnalte-Mur (OAUV)
   Creation date: 01/06/2009
   Last modification date: 03/06/2009

   Program to make transformations on power-spectrum as
   first step of the simulation of log-normal fields. It
   takes the power spectrum we want for the final log-normal 
   field and transforms it to the power spectrum of the 
   corresponding gaussian field to simulate. 

   This is a implementation in 'C' of the program 
   'powesp_transform.py' written in Python. Here, we will 
   use the FFTLog Fortran library by A. Hamilton.


   Input: 
      * File (typically from CAMB) containing the P(k).
      * Value (in Mpc/h) for the scale of the filter to apply
        (we will use, for the moment, a top-hat filter).
      * 'q parameter' (bias exponent) for making the transforms

   Version history:
      V. 0.1 (03/06/2009): Initial version.

      V. 0.1.3 (03/06/2009): For the output P(k)'s, calculate
               3-point averages, and set negative values to zero.
               This should be used with low-resolution spectra
               (or this is my guess). 

      V. 0.2 (03/06/2009): Calculate also the value of sigma (alpha), 
             using the trapezoidal rule, this doesn't need to be too 
             precise. The value of alpha2 will be the output of the
             program (any other comments to stderr).
             Also, write some useful data to a log-file.
 
	  Antoine Labatie:
      V. 0.3 (18/02/2010) Removed all filterings, Kept only constraint
			 on positivity for output P(k). Check the precision by doing
			 the reverse of every step. Precision is checked by comparing Xi[0] 
			 with the integrated P(k).
			 Added a function to calculate power contained between two values kmin and kmax
			 to compare with variance of the gaussian field generated on a grid where kmin 
			 and kmax are known.
			 Added a zero-padding option which prevents the errors in FFTLog 
			 (because of the periodization). Little error remain on small k 
			 (normaly scales larger than usual samples), for which I don't know the reason, 
             but the rest of the curve perfectly fits the input P(k).
	
 
***************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>



#define NMAX  18000  /* Maximum size for k/r/xi/pk arrays */


// Default parameters
double RT = 0.;           /* Radius of filter */
double Q  = 0.;            /* Bias exponent */
double KR = 1;            /* Ratio k_c r_c for transforming the k/r arrays */
double saratio = 1.;  /* Ratio sigma_0 to alpha, for normalization of gaussian field.
                        It shouldn't be relevant for the final process. */
FILE *inpkfile;           /* Input file containing P(k) */
char nameinpkfile[80];    /* Name for input file */
char outprefix[80];		  /* Prefix for output files */
char namexifile[80], namexigausfile[80],  \
  namepkgausfile[80], namepkrecfile[80];   /* Names of different output files */



int ZeroPadding=-1; // -1 means padding on each side with as many points as input array
				   // Having more points than input array seems to cause problems ('NaN' points)
double bias=1.0;

//bool Verbose = False;
int Verbose = 0;

//Declare functions (extern are from fftlog.f)
extern void fhti_(int *n,double *mu,double *q,double *dlnr,double *kr,int *kropt,double *wsave,int *k);
extern void fht_(int *n,double *a,int *dir,double *wsave);

/*********************************************************************/

static void usage(char *argv[])
{

  fprintf(stderr, "Usage: %s options  in_Pk_file\n\n", argv[0]);
  fprintf(stderr, "   where options = \n");

  //Filter radius
  fprintf(stderr, "     [-r FilterRadius]\n");
  fprintf(stderr, "         Radius of the Top-Hat filter to use (in Mpc/h).\n");
  fprintf(stderr, "         Default is %g\n", RT);

  //Bias exponent
  fprintf(stderr, "     [-q BiasExponent]\n");
  fprintf(stderr, "         Bias exponent (q) for making the transforms in FFTLog.\n");
  fprintf(stderr, "         Default is %g\n", Q);

  //Zero padding
  fprintf(stderr, "     [-z ZeroPadding Number of points\n");
  fprintf(stderr, "         Change ZeroPadding number of points.\n");
  fprintf(stderr, "         Default is Same Number of points as input above and below \n");
  fprintf(stderr, "         0 would mean no ZeroPadding \n");

  //Bias
  fprintf(stderr, "     [-b BiasValue\n");
  fprintf(stderr, "         Default is 1.0.\n");
	

  //Verbosity
  fprintf(stderr, "     [-v]\n");
  fprintf(stderr, "         Verbose.\n");

  fprintf(stderr, "\n");
  exit(-1);
}

/*********************************************************************/

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

      //Filter radius
    case 'r':	RT = atof(argv[++i]);
      break;
      
    case 'q':	Q = atof(argv[++i]);
      break;
			
	case 'z':  
			ZeroPadding = atoi(argv[++i]);
	  break;
			
    case 'v': Verbose = 1;
      break;
			
	case 'b': bias = atof(argv[++i]);
	  break;
			
    case '?': usage(argv);
      break;
      
    default:	usage(argv);
      break;
    }
    i++;
  }

  if(i<argc){
    strcpy(nameinpkfile, argv[i]);
  }
  else usage(argv);

}



/*********************************************************************/

/* GET PARAMETERS */

void get_param()
{
	char Name_Param_File[256];
	sprintf(Name_Param_File, "../param/ps_transform.param");	
	FILE *File=fopen(Name_Param_File,"r");
	
    if (File == NULL)
    {
		fprintf(stderr, "Error: cannot open file %s \n\n",Name_Param_File);
		exit(-1);
    }
	char Temp[256];
	int ret;
	ret=fscanf(File, "%s\t%s\n", Temp, outprefix);	//Name_Out prefix
	fclose(File);
	
}

/*********************************************************************/

int write2file(double *xarray, double *yarray, char *filename, int Np)
{
  FILE *fout = fopen(filename, "w"); assert(fout);
	int i=0;
  
	for(i = 0; i<Np; i++){
	  fprintf(fout, "%0.10f  %0.10f\n", xarray[i], yarray[i]);
  }
  
  fclose(fout);
  return 0;
}

/*********************************************************************/

void xifrompk(double *xiout, double *rout, double *pkin, double *kin, double *wsave, int N)
{

  double temparr[NMAX];
  int i;
  //double temptemp;
  for(i=0; i<N;i++){
    //    temparr[i] = pkin[i]*sqrt(kin[i]);
    temparr[i] = pkin[i] * pow(kin[i],1.5);
	temparr[i] = temparr[i] * pow(kin[i],-Q);
    //temptemp = temparr[i];
    //  printf("%g\n", temptemp);
  }

  int dir=1;

  fht_(&N, (double *) temparr, &dir,  (double *) wsave);

  for(i=0; i < N; i++){
    //xiout[i] = pow(2*M_PI,-1.5) * pow(rout[i],-0.5)*temparr[i];
	xiout[i] = temparr[i] * pow(rout[i], -Q);
    xiout[i] = xiout[i] * pow(2*M_PI, -1.5) * pow(rout[i], -1.5);
    //temptemp = xiout[i];
    //printf("%g\n", temptemp);
  }

}

/*********************************************************************/

void pkfromxi(double *pkout, double *kout, double *xiin, double *rin, double *wsave, int N)
{
  double temparr[NMAX];
  int i, dir=-1;
  for(i=0; i < N; i++){
    //    temparr[i] = xiin[i]*sqrt(rin[i]);
    temparr[i] = xiin[i]*pow(2*M_PI, 1.5)*pow(rin[i], 1.5);
	temparr[i] = temparr[i]*pow(rin[i],Q);
  }

  fht_(&N, (double *) temparr, &dir, (double *) wsave);

  for(i=0; i < N; i++){
    //    pkout[i] = pow(2*M_PI, 1.5) * pow(kout[i], -0.5) * temparr[i];
	pkout[i] = temparr[i] * pow(kout[i],Q);
    pkout[i] = pkout[i] * pow(kout[i],-1.5);
  }


}

/*********************************************************************/

double getsigma2(double *pk, double *karr, int N)
{
  int i;
  double sigma2=0;
  //  double integrand0, integrand1;

  for(i=1; i < N; i++)
    {
      sigma2 = sigma2 + (0.5*( (pk[i]*karr[i]*karr[i]) + (pk[i-1]*karr[i-1]*karr[i-1]))*(karr[i] - karr[i-1]));
    }
  sigma2 = sigma2/(2.*M_PI*M_PI);

  return sigma2;
}

/*********************************************************************/

double getsigma2k(double *pk, double *karr, int N,double kmin, double kmax)
{
	int i;
	double sigma2=0;
	//  double integrand0, integrand1;
	
	for(i=1; i < N; i++)
    {
		if(karr[i]>kmin && karr[i]<kmax)
			sigma2 = sigma2 + (0.5*( (pk[i]*karr[i]*karr[i]) + (pk[i-1]*karr[i-1]*karr[i-1]))*(karr[i] - karr[i-1]));
    }
	sigma2 = sigma2/(2.*M_PI*M_PI);
	
	return sigma2;
}

/*********************************************************************/

void filter_pk(double *pkin, double *pkout, double *karr, double Rt, int Nk)
{
	
	int i;
	double x, filt;
	
	for(i=0; i < Nk; i++){
		
		x = karr[i]*Rt;
		//filt = (3./pow(x,3))*(sin(x) - (x*cos(x)));  //Top-Hat filter
		filt = exp(-x*x/2.0);
		pkout[i] = filt*filt*pkin[i];
	}
}

/*********************************************************************/

//Main program!
int main(int argc, char *argv[])
{

  int i=0;
  int Nkin,Nkino, kropt;
  int ok;
  double inpko[NMAX], inkarro[NMAX];  //Original Array Points
  double inpk[NMAX], inkarr[NMAX];    //Zero padded Array Points
  double rarr[NMAX], xilog[NMAX], xigaus[NMAX], pkgaus[NMAX], pkrec[NMAX];
  double pkgausav[NMAX], inkarrsav[NMAX];
  double sigma2, alpha2;
  double logkmax,logkmin;
  double logkc, dlogr,dlnr, logrc;
  double nc;
  int Nwarr;
  double *warray;

  /* Initialization: */

  //Read CL arguments
  get_args(argc,argv);

  //get parameters
  get_param();

  //Open files
  inpkfile = fopen(nameinpkfile, "r");   //Input
  assert(inpkfile);


  //Output: corr.func. of loggauss field
  strcpy(namexifile,outprefix); strcat(namexifile,"_xi.dat");   
  
  //Output: corr.func. of gauss field
  strcpy(namexigausfile, outprefix); strcat(namexigausfile, "_xigaus.dat");
 
  //Output: P(k) of gauss field (main output)
  strcpy(namepkgausfile, outprefix); strcat(namepkgausfile, "_pkgaus.dat");
  
  //Output: for checking purposes, should equal outpkfilt
  strcpy(namepkrecfile, outprefix); strcat(namepkrecfile, "_pkrec.dat");



  //Read input P(k) from file 
  while(!feof(inpkfile))
  {
	  fscanf(inpkfile,"%lf %lf\n", &inkarro[i], &inpko[i]);
	  i++;
  }
  Nkino=i;
	
  //Number of points in the input arrays
  if(ZeroPadding==-1) 
  {
	 ZeroPadding=i;
	 Nkin=3*i;
  }
  else if(ZeroPadding<-1) 
  {
	  ZeroPadding=0;
	  Nkin=i;
  }
  else Nkin = i+2*ZeroPadding;

  if(Verbose) fprintf(stderr, "Input file read, number of points with Zero Padding = %d\n\n", Nkin);


  //Parameters for spacing of points:
  //Follow definitions made in 'fftlogtest.f'
  //First do as if no zero padding
  logkmin = log10(inkarro[0]);
  logkmax = log10(inkarro[Nkino - 1]);
  nc = (Nkin + 1.)/2.;

  logkc = (logkmin + logkmax)/2.;
  dlogr = (logkmax - logkmin)/(Nkino - 1.);
  dlnr = dlogr*log(10.);
	
  //Update considering zero padding
  logkmin = logkmin-ZeroPadding*dlogr;
  logkmax = logkmax+ZeroPadding*dlogr;
	
  for(i=0;i<Nkino;i++)
  {
	  inkarr[ZeroPadding+i]=inkarro[i];
	  inpk[ZeroPadding+i]=bias*bias*inpko[i]; 
  }
  for(i=0;i<ZeroPadding;i++)
  {
	  inkarr[i]=pow(10.0,logkmin+i*dlogr);
	  inpk[i]=0;
	  inkarr[ZeroPadding+Nkino+i]=pow(10,logkmax-(ZeroPadding-(i+1))*dlogr);
	  inpk[ZeroPadding+Nkino+i]=0;
  }
	
  //filtering
  if(RT>0) filter_pk(inpk, inpk, inkarr, RT , Nkin);

  if(Verbose){
    fprintf(stderr, "Spacing parameters: \n");
    fprintf(stderr, "    logkc = %g\n", logkc);
    fprintf(stderr, "    dlogr = %g\n\n", dlogr);
  }

 

  //Initialize FFTLog functions
  double mu=0.5;
  if(Q==0.) Nwarr = (2*Nkin) + (2*(Nkin/2)) + 18;
  else      Nwarr = (2*Nkin) + (3*(Nkin/2)) + 19;
  warray = (double *) malloc(Nwarr*sizeof(double)); assert(warray);
  kropt=0;
  fhti_(&Nkin,&mu,&Q,&dlnr,&KR,&kropt,(double *) warray,&ok);

  if(Verbose)    fprintf(stderr, "kr changed to %.8f\n\n", KR);
  if(Verbose*ok) fprintf(stderr, "Initialization of FFTLog routines went OK.\n\n");
  
  if(!ok){
    fprintf(stderr, "Something went wrong with FFTLog initialization, fhti()! Exiting now!\n");
    exit(-1);
  }

	
  //Define r-array
  logrc = log10(KR) - logkc;
  for(i=0; i < Nkin; i++) rarr[i] = pow(10., logrc + ((i+1-nc)*dlogr));


	
  //Calculate Xi(r)-loggaus
  if(Verbose) fprintf(stderr, "Now calculating Xi(r) of log-gaussian (filtered) field.\n \n");
  xifrompk(xilog, rarr, inpk, inkarr, warray, Nkin);
	
  //Saving Xi(r)-loggaus
  write2file(rarr, xilog, namexifile, Nkin);
  if(Verbose) fprintf(stderr, "Xi(r) for log-gaussian field saved in %s\n\n", namexifile);
  		

  //See precision of FFTlog on log gauss field
  sigma2 = getsigma2(inpk, inkarr, Nkin);
  alpha2 = log(1. + sigma2);
	
  printf("See Precision of FFTLog on the logGauss field (Sigma^2 PS = Xilog[0] ?) \n");
  printf("Sigma^2 = %f \n",sigma2);
  printf("Alpha^2 =ln(1+sigma^2) = %f \n",alpha2);
  printf("Xilog[0]:%f \n \n",xilog[0]);
	
	
  //Transform to Xi(r)-gaussian
  for(i=0;i<Nkin;i++)  xigaus[i] = saratio*saratio*log(1. + xilog[i]);

  //Saving Xi(r)-gaussian
  write2file(rarr,xigaus, namexigausfile, Nkin);  
  if(Verbose)  fprintf(stderr, "Xi(r) for Gaussian field saved in %s\n\n", namexigausfile);
  
	
  //Calculate P(k)-gaussian
  if(Verbose) fprintf(stderr, "Now calculating P(k) of Gaussian field.\n\n");
  pkfromxi(pkgaus, inkarr, xigaus, rarr, warray, Nkin);


  //See precision of FFTlog on Gaussian field
  sigma2 = getsigma2(pkgaus, inkarr,Nkin);
  printf("See Precision of FFTLog on the Gaussian field (Sigma^2 PS = Xilog[0] ?) \n");
  printf("Sigma^2  = %f \n",sigma2);
  printf("Xigauss[0]:%f \n \n",xigaus[0]);
	
  //Processing before output
  //For output file, keep only part of the spectrum that was in the original range (undo zero padding) 
  for(i=0; i < Nkino ;i++) 
  {
	  inkarrsav[i]=inkarr[ZeroPadding+i];
	  pkgausav[i]=pkgaus[ZeroPadding+i];
  }
	
  //Constraint spectrum to be positive
  for(i=0; i < Nkino ;i++){
		if(pkgausav[i] < 0.) pkgausav[i] = 0;
  }
	
  //Saving processed Gaussian P(k)
  write2file(inkarrsav, pkgausav, namepkgausfile, Nkino);
  if(Verbose) fprintf(stderr, "Processed P(k) for Gaussian field saved in %s\n", namepkgausfile);
	
  //See global variance
  sigma2 = getsigma2(pkgausav, inkarrsav,Nkino);
  printf("\n Global variance contained in the output Gaussian Power Spectrum\n");
  printf("Sigma^2  = %f \n",sigma2);
  //sigma2 = getsigma2k(pkgausav, inkarrsav,Nkino,2*M_PI/200.0,2.5*2*M_PI);
  //printf("Sigma^2 in [2pi/200, 2pi*2.5] = %f \n \n",sigma2);

	
  //Reverse everything for checking precision & especially effect of positivity constraint
  //Enforce positivity constraint before reversing to see the effect present in output P(k)-gaussian
  for(i=0; i < Nkin ;i++){
		if(pkgaus[i] < 0.) pkgaus[i] = 0;
	}
	
  if(Verbose) fprintf(stderr, "\nNow inverse everything to see global precision.\n");
  xifrompk(xigaus, rarr, pkgaus, inkarr, warray, Nkin);	
	
  for(i=0;i<Nkin;i++) xilog[i]=exp(xigaus[i]/(saratio*saratio)) -1;
  pkfromxi(pkrec, inkarr, xilog,rarr,warray,Nkin);

  //Saving recovered P(k)-loggauss
  write2file(inkarr, pkrec, namepkrecfile, Nkin);	
  if(Verbose) fprintf(stderr, "Recovered P(k) for log-Gaussian field saved in %s\n\n", namepkrecfile);
  


  if(Verbose) fprintf(stderr, "Program finished!\n");

  return 0;
}


