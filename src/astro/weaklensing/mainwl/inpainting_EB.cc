/******************************************************************************
**                   Copyright (C) 2013 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Sandrine PIRES
**
**    Date: 15/11/13
**    
**    File:  inpainting_EB.cc
**
*******************************************************************************
**
**    DESCRIPTION: Inpainting of Incomplete shear maps using DCT
**       Options :
**          -i# number of iterations (default is 100)
**          -B# Cosinus Transform BlockSize
**          -n# number of wavelet scales (used in sigma bounded)
**          -S# First Threshold level
**          -s# Last Threshold Level
**          -L# 1- Linear decrease
**              2- Exponential decrease
**          -D SigmaBounded
**          -b Bmodes constraint 
**          -f Field in arcmin
**		    
******************************************************************************/
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "NR.h"
#include "IM_DCT.h"
#include <ctime>

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
int Nbr_Plan = -1;  /* number of scales */
float Field = 0.;
int Nbr_Iter = 100;
int DCT_BlockSize = -1;
double FirstThreshold = -1.; // DEF_MCA_FIRST_SOFT_THRESHOLD;
float LastThreshold = 0;
Bool SigmaBounded = False;
Bool Bmodes = False;
Bool Linear = False;
Bool Verbose = False;
type_transform Transform = DEFAULT_MR_TRANS; 

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

/****************************************************************************/
/****************************************************************************/
void keb2g(fltarray &Kappa, double &Field1, double &Field2, fltarray &Gamma)
{
  //time_t tbegin1,tend1;
  //double texec1=0.;
  //tbegin1=time(NULL); 
  int i, j;
  double ii, jj, n1, n2, l1, l2;
  double DPI = 3.1415926536;
  n1 = Kappa.nx();
  n2 = Kappa.ny();
  Ifloat KappaE(n1, n2, "KappaE");
  Ifloat KappaB(n1, n2, "KappaB");
  Icomplex_f KappaE_tilde(n1, n2, "kappa_fft");
  Icomplex_f KappaB_tilde(n1, n2, "kappa_fft");
  Icomplex_f PsiE_tilde(n1, n2, "psi");
  Icomplex_f PsiB_tilde(n1, n2, "psi");
  Icomplex_f GammaE1_tilde(n1, n2, "gammaE1_fft");
  Icomplex_f GammaE2_tilde(n1, n2, "gammaE2_fft");
  Icomplex_f GammaB1_tilde(n1, n2, "gammaB1_fft");
  Icomplex_f GammaB2_tilde(n1, n2, "gammaB2_fft");
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
	KappaE(j,i) = Kappa(i,j,0);
	KappaB(j,i) = Kappa(i,j,1);
      }
    }

  FFTN_2D FFT;
  Bool Reverse = False;
  FFT.fftn2d(KappaE, KappaE_tilde);
  FFT.fftn2d(KappaB, KappaB_tilde);

  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
      ii = double(i);
      jj = double(j);
      l1 = double(2.*DPI/Field1*(ii-(n1/2.)));
      l2 = double(2.*DPI/Field2*(jj-(n2/2.)));  
      if (l1*l1+l2*l2 == 0)
        {
	  PsiE_tilde(j, i) = complex_f(0., 0.);
 	  PsiB_tilde(j, i) = complex_f(0., 0.);
	}
      else
	{
	  PsiE_tilde(j,i) = complex_f( -(2./double(l1*l1+l2*l2))*KappaE_tilde(j,i).real(),
                                  -(2./double(l1*l1+l2*l2))*KappaE_tilde(j,i).imag());
	  PsiB_tilde(j,i)= complex_f( -(2./double(l1*l1+l2*l2))*KappaB_tilde(j,i).real(),
                                 -(2./double(l1*l1+l2*l2))*KappaB_tilde(j,i).imag());
        }
      }
    }
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
        ii = double(i);
        jj = double(j);
        l1 = double(2.*DPI/Field1*(ii-n1/2));
        l2 = double(2.*DPI/Field2*(jj-n2/2));

        GammaE1_tilde(j,i) = complex_f( -double(.5*(l1*l1-l2*l2))*PsiE_tilde(j,i).real(),
                                    -double(.5*(l1*l1-l2*l2))*PsiE_tilde(j,i).imag());
        GammaE2_tilde(j,i) = complex_f( -double(l1*l2)*PsiE_tilde(j,i).real(),
                                       -double(l1*l2)*PsiE_tilde(j,i).imag());
        
        GammaB1_tilde(j,i)= complex_f( -double(.5*(l1*l1-l2*l2))*PsiB_tilde(j,i).real(),
                                      -double(.5*(l1*l1-l2*l2))*PsiB_tilde(j,i).imag());
        GammaB2_tilde(j,i)= complex_f( -double(l1*l2)*PsiB_tilde(j,i).real(),
                                      -double(l1*l2)*PsiB_tilde(j,i).imag());
      }
    }
  Reverse = True;
  FFT.fftn2d(GammaE1_tilde, Reverse);
  FFT.fftn2d(GammaE2_tilde, Reverse);
  FFT.fftn2d(GammaB1_tilde, Reverse);
  FFT.fftn2d(GammaB2_tilde, Reverse);
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
      Gamma(i, j, 0) = GammaE1_tilde(j,i).real()+GammaB2_tilde(j,i).real();
      Gamma(i, j, 1) = GammaE2_tilde(j,i).real()-GammaB1_tilde(j,i).real();
      }
    }

  //tend1=time(NULL); 
  //texec1=difftime(tend1,tbegin1); 
  //cout << "====================================" << endl;
  //cout << "Execution Time keb2g = "  << texec1 << endl;
  //cout << "====================================" << endl;
}

/****************************************************************************/
/****************************************************************************/
void g2keb(fltarray &Gamma, double &Field1, double &Field2, fltarray &Kappa)
{
  //time_t tbegin2,tend2;
  //double texec2=0.;
  //tbegin2=time(NULL); 
  int i, j;
  double ii, jj, n1, n2, l1, l2;
  n1 = Gamma.nx();
  n2 = Gamma.ny();
  Ifloat Gamma1E(n1, n2, "Gamma1E");
  Ifloat Gamma2E(n1, n2, "Gamma2E");
  Ifloat Gamma1B(n1, n2, "Gamma1B");
  Ifloat Gamma2B(n1, n2, "Gamma2B");
  Icomplex_f KappaE_tilde(n1, n2, "kappaE_fft");
  Icomplex_f KappaB_tilde(n1, n2, "kappaB_fft");
  Icomplex_f Gamma1E_tilde(n1, n2, "gamma1E_fft");
  Icomplex_f Gamma2E_tilde(n1, n2, "gamma2E_fft");
  Icomplex_f Gamma1B_tilde(n1, n2, "gamma1B_fft");
  Icomplex_f Gamma2B_tilde(n1, n2, "gamma2B_fft");
  FFTN_2D FFT;
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
	Gamma1E(j,i) = Gamma(i,j,0);
	Gamma2E(j,i) = Gamma(i,j,1);
	Gamma1B(j,i) = -Gamma(i,j,1);
	Gamma2B(j,i) = Gamma(i,j,0);
      }
    }
  Bool Reverse = False;
  FFT.fftn2d(Gamma1E, Gamma1E_tilde);
  FFT.fftn2d(Gamma2E, Gamma2E_tilde);
  FFT.fftn2d(Gamma1B, Gamma1B_tilde);
  FFT.fftn2d(Gamma2B, Gamma2B_tilde);
  double DPI = 3.1415926536;
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
      ii = double(i);
      jj = double(j);
      l1 = double(2.*DPI/Field1*(ii-(n1/2.)));
      l2 = double(2.*DPI/Field2*(jj-(n2/2.)));  
      if (l1*l1+l2*l2 == 0)
        {
	  KappaE_tilde(j, i) = complex_f(0., 0.);
	  KappaB_tilde(j, i) = complex_f(0., 0.);
	}
      else
	{
          KappaE_tilde(j,i) = complex_f(  (double(l1*l1-l2*l2)*Gamma1E_tilde(j,i).real()+double(2.*l1*l2)*Gamma2E_tilde(j,i).real())/double(l1*l1+l2*l2),
              (double(l1*l1-l2*l2)*Gamma1E_tilde(j,i).imag()+double(2.*l1*l2)*Gamma2E_tilde(j,i).imag())/double(l1*l1+l2*l2));
          KappaB_tilde(j,i) = complex_f( (double(l1*l1-l2*l2)*Gamma1B_tilde(j,i).real()+double(2.*l1*l2)*Gamma2B_tilde(j,i).real())/double(l1*l1+l2*l2),
              (double(l1*l1-l2*l2)*Gamma1B_tilde(j,i).imag()+double(2.*l1*l2)*Gamma2B_tilde(j,i).imag())/double(l1*l1+l2*l2));
        }
      }
    }
  
  Reverse = True;
  FFT.fftn2d(KappaE_tilde, Reverse);
  FFT.fftn2d(KappaB_tilde, Reverse);
    
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
	Kappa(i, j, 0) = KappaE_tilde(j,i).real();
	Kappa(i, j, 1) = KappaB_tilde(j,i).real();
      }
    }
  //tend2=time(NULL); 
  //texec2=difftime(tend2,tbegin2); 
  //cout << "====================================" << endl;
  //cout << "Execution Time g2keb = "  << texec2 << endl;
  //cout << "====================================" << endl;
}

/****************************************************************************/
/****************************************************************************/

void residual_wl(fltarray &Gamma, double &Field1, double &Field2, fltarray &Kappa, Ifloat &Mask, fltarray &Residual, int &Nbr_Plan, Bool &SigmaBounded)
{
  int i, j, l, n1, n2;
  float Aint_1, Aint_2, Bint_1, Bint_2, Cint;
  float Aout_1, Aout_2, Bout_1, Bout_2, Cout;
  double sigint2_1, sigint2_2, sigout2_1, sigout2_2, sigint, sigout;
  n1 = Gamma.nx();
  n2 = Gamma.ny();
  Ifloat Gamma_i1(n1, n2, "Gamma1E");
  Ifloat Gamma_i2(n1, n2, "Gamma2E");
  MultiResol MR_Gamma1(n1, n2, Nbr_Plan, Transform, "MRModel1");
  MultiResol MR_Gamma2(n1, n2, Nbr_Plan, Transform, "MRModel2");
  fltarray Gamma_Rec(n1, n2, 2, "gamma_rec");
  fltarray Gamma_i(n1, n2, 2, "gamma_rec");


  keb2g(Kappa, Field1, Field2, Gamma_i);
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
      Gamma_i1(j,i) = Gamma_i(i,j,0);
      Gamma_i2(j,i) = Gamma_i(i,j,1);
      }
    }

  if (SigmaBounded == True) 
    {
    MR_Gamma1.transform(Gamma_i1);
    MR_Gamma2.transform(Gamma_i2);
   for (l=0; l<Nbr_Plan; l++)
      {
      for (i=0; i<n1; i++)
        {
        for (j=0; j<n2; j++)
          {
	    if (Mask(j, i) == 0)
            {
	    Aint_1 = Aint_1 + MR_Gamma1(l, j, i);
            Aint_2 = Aint_2 + MR_Gamma2(l, j, i);
            Bint_1 = Bint_1 + MR_Gamma1(l, j, i)*MR_Gamma1(l, j, i);
            Bint_2 = Bint_2 + MR_Gamma2(l, j, i)*MR_Gamma2(l, j, i);
            Cint = Cint + 1.;
            }
          else
            {
	    Aout_1 = Aout_1 + MR_Gamma1(l, j, i);
            Aout_2 = Aout_2 + MR_Gamma2(l, j, i);
            Bout_1 = Bout_1 + MR_Gamma1(l, j, i)*MR_Gamma1(l, j, i);
            Bout_2 = Bout_2 + MR_Gamma2(l, j, i)*MR_Gamma2(l, j, i);
            Cout = Cout + 1.;
	    }
	  }
	}
      sigint2_1 = Bint_1/Cint - (Aint_1/Cint)*(Aint_1/Cint);
      sigint2_2 = Bint_2/Cint - (Aint_2/Cint)*(Aint_2/Cint);
      sigint = sqrt(sigint2_1 + sigint2_2);
      sigout2_1 = Bout_1/Cout - (Aout_1/Cout)*(Aout_1/Cout);
      sigout2_2 = Bout_2/Cout - (Aout_2/Cout)*(Aout_2/Cout);
      sigout = sqrt(sigout2_1 + sigout2_2);
      if (sigint > sigout)
        {
        for (i=0; i<n1; i++)
          {
          for (j=0; j<n2; j++)
            {
	      if (Mask(j,i) == 0)
              {
	      MR_Gamma1(l, j, i) = MR_Gamma1(l, j, i)*(sigout/sigint);
	      MR_Gamma2(l, j, i) = MR_Gamma2(l, j, i)*(sigout/sigint);
              }
	    }
       	  }
        }
      }
    MR_Gamma1.recons(Gamma_i1);  
    MR_Gamma2.recons(Gamma_i2);  
    }

  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
	Gamma_Rec(i,j, 0) = (Gamma(i,j,0) - Gamma_i1(j, i))*Mask(j, i);
	Gamma_Rec(i,j, 1) = (Gamma(i,j,1) - Gamma_i2(j, i))*Mask(j, i);
      }
    }
   g2keb(Gamma_Rec, Field1, Field2, Residual);

 }

/****************************************************************************/
/****************************************************************************/
static void usage(char *argv[])
{    
    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-f Field]\n");
    fprintf(OUTMAN, "             Field in arcmin. Default is %f.\n", Field);    
    manline();

    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration. Default is %d.\n", Nbr_Iter);    
    manline();
 
    fprintf(OUTMAN, "         [-B DCT_ockockSize]\n");
    fprintf(OUTMAN, "             Local DCT block size.\n");
    fprintf(OUTMAN, "             By default, a local DCT (block size = 32) is used. \n");    
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the WT.\n");
    fprintf(OUTMAN, "             Default is derived from the data.\n");    
    manline();

    fprintf(OUTMAN, "         [-S FirstThresholdLevel]\n");
    fprintf(OUTMAN, "            First thresholding value.\n");
    fprintf(OUTMAN, "            Default is derived from the data.\n");    
    manline();

    fprintf(OUTMAN, "         [-s LastThresholdLevel]\n");
    fprintf(OUTMAN, "            Last thresholding value..\n");
    fprintf(OUTMAN, "            default is %f.\n", LastThreshold);    
    manline();

    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "             Linear descent\n");
    fprintf(OUTMAN, "             Default is exponential descent.\n");  
    manline();

    fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "             Sigma Bounded\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();   

    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Bmodes constraint\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();       
           
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n"); 
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
  {
  int c;
  /* get options */
  while ((c = GetOpt(argc,argv,"f:i:B:n:S:s:LDbv")) != -1) 
    {
    switch (c) 
      {
      case 'f':
      if (sscanf(OptArg,"%f",&Field) != 1) 
        {
        fprintf(OUTMAN, "bad Field size: %s\n", OptArg);
        exit(-1);
        }
      break;	 
      case 'i':
      if (sscanf(OptArg,"%d",&Nbr_Iter) != 1) 
        {
        fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
        exit(-1);
        }
      break;	 
      case 'B':
      if (sscanf(OptArg,"%d",&DCT_BlockSize) != 1) 
        {
        fprintf(OUTMAN, "bad cosinus transform block size: %s\n", OptArg);
        exit(-1);
        }
      break;
      case 'n':
      if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
	{
        fprintf(stderr, "Error: bad number of scales: %s\n", OptArg);
        exit(-1);
        }
      break; 
      case 'S':               
      if (sscanf(OptArg,"%lf",&FirstThreshold) != 1) 
	{
        fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
        exit(-1);
	}
      if (FirstThreshold < 0)
	{
        fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
        exit(-1);
	}
      break;
      case 's':               
      if (sscanf(OptArg,"%f",&LastThreshold) != 1) 
	{
        fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
        exit(-1);
	}
      break;
      case 'L': Linear = True; break;
      case 'D': SigmaBounded = True; break;	
      case 'b': Bmodes = True; break;
      case 'v': Verbose = True; break;
      case '?': usage(argv); break;
      default : usage(argv); break;
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
     fprintf(stderr, "Too many parameters: %s ...\n", argv[OptInd]);
     exit(-1);
     }
  }

/*********************************************************************/
int main(int argc, char *argv[])
  {
  time_t tbegin1,tend1;
  double texec1=0.;
  tbegin1=time(NULL); 
  time_t tbegin2,tend2;
  double texec2=0.;

  long i, j, k, n1, n2;
  int Ndct1, Ndct2, nn1, nn2, i0, j0, iter;
  float Threshold, lambda;
  //Ifloat Imag;    
  fltarray Gamma;
  /* support image creation */
  fitsstruct Header;
  char Cmd[256];
  extern softinfo Soft;
  Soft.mr2();
  Cmd[0] = '\0';
  for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
  lm_check(LIC_MR2);
  /* Get command line arguments, open input file(s) if necessary */
  filtinit(argc, argv);     
  
    // Read the input image
  io_3d_read_data(Name_Imag_In, Gamma, &Header);

  Header.origin = Cmd;

  cout << "Original Image size : n1 = " << Gamma.nx() << endl;
  
  n1 = Gamma.nx();
  n2 = Gamma.ny();
  
  if (( fmod(log(n1),log(2.)) > 1e-8) ||( fmod(log(n2),log(2.)) > 1e-8))
    {
    nn1 = pow(2.,  int(log(n1)/log(2.))+1);
    nn2 = pow(2.,  int(log(n2)/log(2.))+1);
    i0 = int((nn1 - n1)/2);
    j0 = int((nn2 - n2)/2);
    }
  else
    {
    nn1 = n1;
    nn2 = n2;
    i0 = 0;
    j0 = 0;
    }

  if (Bmodes == True)
    {
    n1 = pow(2., int(log(nn1)/log(2.))+1.);
    n2 = pow(2., int(log(nn2)/log(2.))+1.);
    i0 = i0 + int((n1 - nn1)/2);
    j0 = j0 + int((n2 - nn2)/2);
    Field = Field *2.;
    }
  else
    {
    n1 = nn1;
    n2 = nn2;
    }

  if (n1 != n2)
    {
    cout << "Image is not square" << endl;
    }

  Ifloat Mask_Total(n1, n2,"Mask_Total");
  Ifloat Mask_Border(n1, n2,"Mask_Border");
  fltarray Gamma_2(n1, n2, 2, "Gamma_2");

  for (i=0; i<Gamma.nx(); i++)
    {
      for (j=0; j<Gamma.ny(); j++)
      {
      Gamma_2(i+i0, j+j0, 0) = Gamma(i, j, 0);
      Gamma_2(i+i0, j+j0, 1) = Gamma(i, j, 1);
      Mask_Border(j+j0, i+i0) = 1;
      if (Gamma(i, j, 0) != 0) 
	{
	  Mask_Total(j+j0,i+i0) = 1.;
	}
      }
    }
  cout << "Processing Image size : n1 = " << n1 << endl;
 
  double Field1 = Field/60.;
  double Field2 = Field/60.;
  Ifloat KappaE(n1, n2,"KappaE");
  Ifloat KappaB(n1, n2,"KappaE");
  Ifloat ResultE(n1, n2,"Result");
  Ifloat ABSResultE(n1, n2,"Result");
  Ifloat ResultB(n1, n2,"Result");
  Ifloat ResultE_tempo(n1, n2,"Result");
  Ifloat ResultB_tempo(n1, n2,"Result");
  fltarray Kappa(n1, n2, 2, "Kappa");
  fltarray Kappa_out(Gamma.nx(), Gamma.ny(), 2, "Kappa_out");
  fltarray out_data(n1, n2, 2, "out_data");
  fltarray Residual(n1, n2, 2, "Residual");
  Ifloat tempo_dataE(n1, n2,"tempo_dataE");
  Ifloat tempo_dataB(n1, n2,"tempo_dataB");
  Bool BlockOverlap = True;
  Bool Weight = False;
  LOCAL_DCT2D LDCT;
  LOCAL_DCT2D LDCTE;
  LOCAL_DCT2D LDCTB;
  LDCT.alloc(n1, n2, DCT_BlockSize, BlockOverlap, Weight);
  LDCTE.alloc(n1, n2, DCT_BlockSize, BlockOverlap, Weight);
  LDCTB.alloc(n1, n2, DCT_BlockSize, BlockOverlap, Weight);
  Bool Reverse;

  //Mass Inversion
  g2keb(Gamma_2, Field1, Field2, Kappa);
  for (i=0; i<n1; i++)
    {
    for (j=0; j<n2; j++)
      {
	KappaE(j, i) = Kappa(i, j, 0); 
 	KappaB(j, i) = Kappa(i, j, 1); 
      }
    }
  if (Nbr_Plan == -1)
    {
      Nbr_Plan = int(log(n1)/log(2.))-3.-2.;
    }
  if (DCT_BlockSize == -1)
    {
      DCT_BlockSize = (n1 < n2) ? n1 : n2 ;
    }

  if (FirstThreshold == -1)
    {
    if ((DCT_BlockSize == n1) && (DCT_BlockSize == n2))
      {
      Reverse = False;
      im_dct(KappaE, ResultE, Reverse);
      ResultE(0, 0) = 0;
      }
    else
      {
      LDCT.transform(KappaE);
      ResultE = LDCT.trans_ima();
      cout << "ndct= "<< ResultE.nx() << endl;
      Ndct1 = int(ResultE.nx()/DCT_BlockSize);
      Ndct2 = int(ResultE.ny()/DCT_BlockSize);
      for (i=0; i < Ndct1; i++)
        {
	for (j=0; j < Ndct2; j++)
          {
	  ResultE(j*DCT_BlockSize, i*DCT_BlockSize) = 0;
	  }
        }
      }
    FirstThreshold = abs(ResultE.max()) > abs(ResultE.min()) ? abs(ResultE.max()) : abs(ResultE.min());
    }

if (Verbose == True)
    {  
    cout << "====================================" << endl;
    cout << "PARAMETERS: INPAINTING WITH DCT" << endl;
    cout << "====================================" << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Number of iterations = " << Nbr_Iter << endl;
    cout << "DCT BlockSize = " << DCT_BlockSize << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "First Threshold = " << FirstThreshold << endl;
    cout << "Last Threshold = " << LastThreshold << endl;
    cout << "Field in degrees = " << Field/60. << endl;
    }

  //Initialization
  FirstThreshold = double(FirstThreshold);
  lambda = FirstThreshold; 
  residual_wl(Gamma_2, Field1, Field2, out_data, Mask_Total, Residual, Nbr_Plan, SigmaBounded);

  for (iter = 0; iter < Nbr_Iter; iter++)
    {
    tbegin2=time(NULL); 
    cout << "====================================" << endl;
    cout << "Iteration n = " << iter << ", lambda = " << lambda << endl;
    for (i=0; i<n1; i++)
      {
      for (j=0; j<n2; j++)
        {
	  tempo_dataE(j, i) = out_data(i, j, 0) + Residual(i, j, 0);
          tempo_dataB(j, i) = out_data(i, j, 1) + Residual(i, j, 1);
	}
      }
    if (lambda > 0) 
      {
        Threshold = lambda;
	//THRESHOLD
      if ((DCT_BlockSize == n1) && (DCT_BlockSize == n2))
        {
        Reverse = False;
        im_dct(tempo_dataE, ResultE, Reverse);
        im_dct(tempo_dataB, ResultB, Reverse);
        float A = ResultE(0, 0);
        float B = ResultB(0, 0);
        for (i=0; i<ResultE.nx(); i++)
          {
          for (j=0; j<ResultE.ny(); j++)
            {
	    if (abs(ResultE(j, i)) < Threshold) 
	      { 
		ResultE(j, i) = 0;
	      }
	    if (abs(ResultB(j, i)) < Threshold) 
	      { 
		ResultB(j, i) = 0;
	      }
	    }
	  }
        ResultE(0,0) = A;
        ResultB(0,0) = B;
        Reverse = True;
        im_dct(ResultE, tempo_dataE, Reverse);
        im_dct(ResultB, tempo_dataB, Reverse);
        }
      else
	{
        LDCTE.transform(tempo_dataE);
        LDCTB.transform(tempo_dataB);
        ResultE = LDCTE.trans_ima();
        ResultB = LDCTB.trans_ima();
        for (i=0; i<ResultE.nx(); i++)
          {
          for (j=0; j<ResultE.ny(); j++)
            {
	    ResultE_tempo(j, i) = ResultE(j,i);
	    ResultB_tempo(j, i) = ResultB(j,i);
	    if (abs(ResultE(j, i)) < Threshold) 
	      { 
		ResultE(j, i) = 0;
	      }
	    if (abs(ResultB(j, i)) < Threshold) 
	      { 
		ResultB(j, i) = 0;
	      }
	    }
	  }
        for (i=0; i < Ndct1; i++)
          {
	  for (j=0; j < Ndct2; j++)
            {
	    ResultE(j*DCT_BlockSize, i*DCT_BlockSize) = ResultE_tempo(j*DCT_BlockSize, i*DCT_BlockSize);
	    ResultB(j*DCT_BlockSize, i*DCT_BlockSize) = ResultB_tempo(j*DCT_BlockSize, i*DCT_BlockSize);
	    }
          }
        LDCTE.set_dct(ResultE);
        LDCTB.set_dct(ResultB);
        LDCTE.recons(tempo_dataE);
        LDCTB.recons(tempo_dataB);
	}
      for (i=0; i<n1; i++)
        {
        for (j=0; j<n2; j++)
          {
	    out_data(i, j, 0) = tempo_dataE(j, i);
            out_data(i, j, 1) = tempo_dataB(j, i);
	  }
	}

      if (Bmodes == True)
	{
        for (i=0; i < n1; i++)
          {
	  for (j=0; j < n2; j++)
            {
	    if (Mask_Border(j, i) == 0)
	      {
		out_data(i, j, 1) = 0;
	      }
	    }
	  }
	 }
      residual_wl(Gamma_2, Field1, Field2, out_data, Mask_Total, Residual, Nbr_Plan, SigmaBounded);
      }
    //lambda evolution
    if (Linear == True)
      {
      lambda = lambda - (FirstThreshold - LastThreshold)/Nbr_Iter;
      }
    else
      {
      lambda = LastThreshold + (FirstThreshold - LastThreshold) *(1. -erff(2.8*float(iter+1)/Nbr_Iter));
      }
    if ((lambda < LastThreshold) || (iter == Nbr_Iter-2)) 
      {
      lambda = LastThreshold;
      }
  tend2=time(NULL); 
  texec2=difftime(tend2,tbegin2); 
  cout << "=> Loop Time = "  << texec2 << endl;

    }

  for (i=0; i<Gamma.nx(); i++)
    {
      for (j=0; j<Gamma.ny(); j++)
      {
      Kappa_out(i, j, 0) =  out_data(i+i0, j+j0, 0);
      Kappa_out(i, j, 1) =  out_data(i+i0, j+j0, 1);
      }
    }


  io_3d_write_data(Name_Imag_Out, Kappa_out);

  tend1=time(NULL); 
  texec1=difftime(tend1,tbegin1); 
  cout << "====================================" << endl;
  cout << "Total Execution Time = "  << texec1 << endl;
  cout << "====================================" << endl;
  exit(0);
} 
