/******************************************************************************
**                   Copyright (C) 2003 by CEA + Valencia observatory
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Enn Saar, Jean-Luc Starck and Vicent Martinez
**
**    Date:  27/02/03
**    
**    File:  genus.cc
**
*******************************************************************************
**
**    DESCRIPTION  genus program
**    ----------- 
**                 
******************************************************************************/

//ES are comments by ES

#include "Array.h"
#include "IM_IO.h"
#include "NR.h"
#include "DefPoint.h"
#include "CUBE_Segment.h"
#include "DefFunc.h"
#include "Genus3D.h"
#include "FFTN_3D.h"


#define WAVELET 1

#ifdef WAVELET    
#include "Atrou3D.h"
#endif

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Dens_Out[256]; /* output density file name */
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float Sigma=0;
//ES changed Vfstep to Mfstep
//float VfStep = 0.02;
float Mfstep = 0.01;
//ES added Periodic
Bool Periodic=True;
Bool SimuPois=False;
Bool SimuGauss=False;

#define SIMU_DIM 32
#define SIMU_GAUSS_INDEX  (double) -1.

double NuIndex = SIMU_GAUSS_INDEX;
float POWA=1;				/* spectral amplitude	*/
int SimuDim = SIMU_DIM;
float Step=1.;
Bool Normalize = False;
Bool Segmentation = False;
Bool WriteAll = False;
Bool OnlyPos = False;
#ifdef WAVELET    
int Nbr_Plan = 4;
#endif     

/***************************************/

extern void convolve(fltarray & Data, fltarray & Gauss); 

   
/***************************************************************************/
 
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
//ES changed nustep to vfstep
//ES we use Mfstep here
#ifdef WAVELET    
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Default is 4.\n");
#endif     
    
    fprintf(OUTMAN, "         [-M Mfstep]\n");
    fprintf(OUTMAN, "             Mfstep.Default is 0.01.\n");
//ES added the periodicity option
//    fprintf(OUTMAN, "         [-P]\n");
//    fprintf(OUTMAN, "             Periodic data. Default is nonperiodic.\n");
//    manline();
    fprintf(OUTMAN, "         [-c Std]\n");
    fprintf(OUTMAN, "             Convolve the input data with a Gaussian width sigma=Std.\n");
    fprintf(OUTMAN, "         [-l Lambda]\n");
    fprintf(OUTMAN, "             Convolve the input data with a Gaussian width sigma=sqrt(2)Lambda.\n");
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "            Poisson distribution simulation.\n");
    fprintf(OUTMAN, "         [-g]\n");
    fprintf(OUTMAN, "            Gaussian field distribution simulation.\n");
    fprintf(OUTMAN, "         [-G IndexVal]\n");
    fprintf(OUTMAN, "            Power index in the simulated Gaussian field distribution.\n");
    fprintf(OUTMAN, "            Default is -1.\n");
    fprintf(OUTMAN, "         [-D Dimension]\n");
    fprintf(OUTMAN, "            Cube size dimension in the simulation.\n");
    fprintf(OUTMAN, "            Default is 32.\n");
    fprintf(OUTMAN, "         [-s Step]\n");
    fprintf(OUTMAN, "             Bin size. \n");
    fprintf(OUTMAN, "             default is %f. \n", Step);    
    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "         [-w FileName]\n");
    fprintf(OUTMAN, "             Write the data convolved with the Gaussian kernel.\n");
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    fprintf(OUTMAN, " Output file contains the four minkovski functionals (MF) for \n");
    fprintf(OUTMAN, "different theshold levels:\n");
    fprintf(OUTMAN, " OUT[*,0] = MF(0) = Area of the surface.\n");
    fprintf(OUTMAN, " OUT[*,1] = MF(1) = Volume enclosed by the surface.\n");
    fprintf(OUTMAN, " OUT[*,2] = MF(2) = Integrated mean curvature of the surface.\n");
    fprintf(OUTMAN, " OUT[*,3] = MF(3) = Euler characteristic.\n");
    fprintf(OUTMAN, " OUT[*,4] = VF = Percent. of Volume larger than the level, \n");
    fprintf(OUTMAN, " OUT[*,5] = Threshold Level, \n");
    fprintf(OUTMAN, " OUT[*,6] = Nu = erf^{-1}(VF)  \n");
    fprintf(OUTMAN, "   Genus  = 1 - OUT[*,3] / 2.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void genusinit(int argc, char *argv[])
{
    int c; 
    float Val; 
    int seed;

    /* get options */
#ifndef WAVELET    
    while ((c = GetOpt(argc,argv,"MTw:A:I:SNl:s:D:G:c:pgvzZ")) != -1) 
#else
    while ((c = GetOpt(argc,argv,"MTw:A:I:SNl:s:D:G:c:pgn:vzZ")) != -1) 
#endif
    {
	switch (c) 
        { 
	  case 'T': OnlyPos = (OnlyPos == True) ? False: True; break;
	  case 'A': 
                if (sscanf(OptArg,"%f",&POWA) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad spectral amplitude parameter: %s\n", OptArg);
                    exit(-1);
                }
                 break;
#ifdef WAVELET
          case 'n': 
                if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad seed parameter: %s\n", OptArg);
                    exit(-1);
                }
		init_random((unsigned int) seed);
                break;
#endif
          case 'I': 
                if (sscanf(OptArg,"%d",&seed) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad seed parameter: %s\n", OptArg);
                    exit(-1);
                }
		init_random((unsigned int) seed);
                break;
          case 'S': Segmentation = True; break;
  	  case 'N': Normalize = True; break;
          case 's': 
	        if (sscanf(OptArg,"%f",&Step) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad step value: %s\n", OptArg);
                    exit(-1);
                }
                break;
          case 'D': 
	        if (sscanf(OptArg,"%d",&SimuDim) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad dimension parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'G': 
	        if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad index parameter: %s\n", OptArg);
                    exit(-1);
                }
		NuIndex = Val;
                break;
//ES added that
	   case 'P': Periodic = True; break;
	   case 'p': SimuPois = (SimuPois == True) ? False: True; break;
           case 'g': SimuGauss = (SimuGauss == True) ? False: True; break;
//ES replaced Vfstep by Mfstep
	   case 'M': 
                if (sscanf(OptArg,"%f",&Mfstep) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Mfstep parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
            case 'c': 
                if (sscanf(OptArg,"%f",&Sigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Std parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	    case 'l': 
                if (sscanf(OptArg,"%f",&Sigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Std parameter: %s\n", OptArg);
                    exit(-1);
                }
		Sigma *= sqrt(2.);
                break;
            case 'w':   
	           if (sscanf(OptArg,"%s", Name_Dens_Out) != 1) 
                   {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                   }
		   WriteAll = True; 
	           break;
	    case 'v': Verbose = True; break;
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
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    genusinit(argc, argv);
      
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
//ES replaced Vfstep by Mfstep
	cout << "# MFSTEP = " << Mfstep    << endl;  
	if (Sigma > 0) cout << "# Convolve the data with a Gaussian: Sigma = " << Sigma    << endl;  
        if ((SimuPois == True) || (SimuGauss == True))
	{
	    cout << "# Simulation mode: " << endl;
	    if (SimuGauss == True) cout << "#    Gaussian field with index power = " << NuIndex << endl;
	    if (SimuGauss == True) cout << "#    Gaussian field with spectral amplitude = " << POWA << endl;
	    else cout << "#    Poisson distribution " << endl;
	    cout << "#    Cube size for the simulated data = " << SimuDim << endl;
        }
    }
    
   
    fltarray Data;
    fltarray Result;
    fltarray Gauss;
    
    // Read the data oe simulate data
    if ((SimuPois == False) && (SimuGauss == False))
    {
        read_data(Data,Name_Imag_In,Step);
        // fits_read_fltarr(Name_Imag_In, Data, &Header);
    }
    else 
    {
        Data.alloc(SimuDim, SimuDim, SimuDim);
        if (SimuGauss == True)  gengauss(Data, NuIndex, POWA);
        else  genpois(Data);
	// if (WriteAll == True) fits_write_fltarr("xx_data.fits", Data);
    }
    
    // Convolution with a Gaussian
    if (Sigma > 0)
    {
        // Make the 3D Gaussian function
        Gauss.alloc(Data.nx(), Data.ny(), Data.nz());
        simu_gaussian_noise3d(Gauss, Sigma); 
	// if (WriteAll == True) fits_write_fltarr("xx_gauss.fits", Gauss);
	
	// Convolve the data with the Gaussian function
	convolve3d(Data, Gauss);
    }
    initfield(&Header);
    Header.origin = Cmd;
    Header.bitpix = -32;   
    Header.naxis = Data.naxis();
    Header.npix = Data.n_elem();   
    Header.width =  Data.nx();
    Header.height = Data.ny();
    for (int i=0; i < Header.naxis; i++) Header.TabAxis[i] = Data.axis(i+1);

    if (WriteAll == True) fits_write_fltarr(Name_Dens_Out, Data, &Header);
   // if (WriteAll == True) fits_write_fltarr(Name_Dens_Out, Data);

    if (OnlyPos == True)
    {
       for (int i=0; i < Data.nx(); i++)
       for (int j=0; j < Data.ny(); j++)
          if (Data(i,j,k) < 0) Data(i,j,k) = 0;
    }
    // Normalization
    if (Normalize == True)
    {
       float SigmaData = Data.sigma();
       if (Verbose == True) cout << " Sigma [convolved] Data = " << SigmaData << endl;
       if (SigmaData == 0) SigmaData = 1.;
       for (int i=0; i < Data.nx(); i++)
       for (int j=0; j < Data.ny(); j++)
       for (int k=0; k < Data.nz(); k++) Data(i,j,k) /= SigmaData;
    }
    
    // Calculate the curves
#ifndef WAVELET    
    mf_curves(Data, Result, Mfstep, Periodic,  Segmentation, Step, Verbose);
#else
    ATROUS_3D_WT AWT;
    fltarray *TabBand;
    AWT.alloc(TabBand, Data.nx(),Data.ny(),Data.nz(),  Nbr_Plan);
    
    if (Verbose == True) cout << "Wavelet transform: Number of scales = "  << Nbr_Plan << endl;
    AWT.transform(Data, TabBand, Nbr_Plan);
    
    for (int s=0;s < Nbr_Plan; s++)
    {
       fltarray ResScale;
       if (Verbose == True) cout << "Scale " << s+1 << endl;
       mf_curves(TabBand[s], ResScale, Mfstep, Periodic, Segmentation, Step, Verbose);
       if (Verbose == True) cout << "   Result size = " << ResScale.nx() <<  " " << ResScale.ny() << endl;
       if (s==0) Result.alloc(ResScale.nx(), ResScale.ny(), Nbr_Plan);
       for (int i=0; i < ResScale.nx(); i++)
       for (int j=0; j < ResScale.ny(); j++) 
               if (i < Result.nx()) Result(i,j,s) = ResScale(i,j);
    }
    AWT.free(TabBand,Nbr_Plan);
#endif

    // Save the finale result  
    HD1.origin = Cmd;
    HD1.bitpix = -32; 
    HD1.naxis = Result.naxis();
    HD1.npix = Result.n_elem();     
    HD1.width =  Result.nx();
    HD1.height = Result.ny();
    for (int i=0; i < HD1.naxis; i++) HD1.TabAxis[i] = Result.axis(i+1);
    fits_write_fltarr(Name_Imag_Out, Result, &HD1);
    
    // Calculate the theoritical MinFun curve
//     if (SimuGauss == True)
//     {
//         float Sig=Data.sigma();
//         fltarray TheoMinFun(Result.nx(), 5);
// 	// Take the same nu values as in the data
// 	for (int i=0; i < Result.nx(); i++) TheoMinFun(i,4) = Result(i,6);
// 	// theoritical MinFun curve
// 	make_theo_genus(TheoMinFun, NuIndex, Sigma);
// 	// Save the result in fits format
//  	if (WriteAll == True) fits_write_fltarr("xx_theo_genus.fits", TheoMinFun);
//     }
    exit(0);
}

