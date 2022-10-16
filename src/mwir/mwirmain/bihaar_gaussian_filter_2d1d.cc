/**************************************************
 * PROJECT : MultiWavelenImageResto (MWIR) - version II uniquely for 2D+1D image and gaussian case
 * CEA
 * Filename : maing2.cc
 * This is the main file of MWIR
 **************************************************/

#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "cdflib.h"
#include "GlobalInc.h"
#include "IM_IO.h"
#include "Array.h"
#include "MSVST_Filter.h"
#include "Wavelet1.h"
#include "randomc.h"             
          
using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];      // input image file
char   Name_Imag_Out[256]; // output image file
enum   MODE1 { GAUSS = 0, GAUSSFDR = 1 };  // denoising modes
enum   MODE2 { DIRECT = 0, HSD = 1 }; // iteration modes
MODE1  PROGMODE = GAUSS;    // denoising mode
MODE2  ITERMODE = DIRECT;   // iteration mode
int    NITER = 0;           // max. number of iterations
double PROBA = 0.0;         // individual proba.
double FDR = 0.01;          // individual FDR
int    NSCALEXY = 2;        // max. number of scales along X and Y
int    NSCALEZ  = 2;        // max. number of scales along Z
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, true}; // X, Y - undecimated; Z - decimated
bool   AUTOFDR = true;      // automatically choose the FDR mode
bool   FDRINDEP = false;    // independence parameter in FDR
double Sigma = -1.;         // noise standard deviation - decide automatically
type_border T_BORDER = I_MIRROR; // border type
bool BONF = false;          // Bonferroni's multi-test schema
int  WAVE = 0;              // Wavelet used in the denoising [0, 1]
int  MS = 0;                // Wavelet used in the iterative restoration [0, 1]
int    FSCALEXY = 1;          // First detection XY scale
int    FSCALEZ = 1;          // First detection Z scale
double GLEVEL = 0;          // Gauss-detection level (k-sigma)

// NOTICE : Bonferroni and FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name" << endl;
  cout << "NOTICE : This program is only for pseudo-3D (2D+1D) image" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-T value] wavelet used in the denoising" << endl;
  cout << "         T = 0 : Haar (default)" << endl;
  cout << "         T = 1 : Daubechies - 4" << endl;     
  cout << "    [-S value]" << endl;
  cout << "         S = 0 : use the same wavelet as in the denoising to the iterative restoration (default)" << endl;
  cout << "         S = 1 : use ODEGARD 9/7 wavelet in the iterative restoration" << endl;
  cout << "    [-g value] standard deviation used in denoising (default = decide automatically) " << endl;
  cout << "    [-D] to decimate X and Y directions in the wavelet transform (not default)" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : denoise (default)" << endl;
  cout << "         M = 1 : denoise with FDR" << endl;
  cout << "    [-E value] two sided Pb or FDR (default = decide automatically)" << endl;
  cout << "    [-s value] Gauss-detection level (default = automatically decided)" << endl;
  cout << "    [-c value] manual FDR setting, c = 0/1 means independence/dependence of coef. (default = decide automatically)" << endl;
  cout << "    [-n value] max. number of scales along X and Y (default = 2)" << endl;
  cout << "    [-N value] max. number of scales along Z (default = 2)" << endl;
  cout << "    [-F value] first detection X-Y scale (default = 1)" << endl;
  cout << "    [-f value] first detection Z scale (default = 1)" << endl;
  cout << "    [-I value] " << endl;
  cout << "         I = 0 : Direct iterative mode (default)" << endl;
  cout << "         I = 1 : Hybrid steepest descent iterative mode" << endl;
  cout << "    [-i value] max. of number of iteration (default = 0)" << endl; 
  cout << "    [-B value] border type" << endl;
  cout << "         B = 0 : symmetric border (default)" << endl;
  cout << "         B = 1 : periodic border" << endl;
  cout << "    [-R] Bonferroni's multi-test schema (not default and no effect in FDR)" << endl;
  cout << "    [-v] verbose mode" << endl;
  cout << endl;
}
 
//*********************************************************************

// GET COMMAND LINE ARGUMENTS
static void filtinit(int argc, char *argv[])
{
    int c;  

    // get options 
    while ((c = GetOpt(argc,argv,"T:S:g:DM:E:s:c:n:N:F:f:I:i:B:Rv")) != -1) 
    {
      switch (c) 
	{
	case 'S':
	  if (sscanf(OptArg, "%d", &MS) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  

	case 'T':
	  if (sscanf(OptArg, "%d", &WAVE) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;

	case 'g':
	  if (sscanf(OptArg, "%lf", &Sigma) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  

	case 'D':
	  DEC[0] = true; DEC[1] = true;
	  break;

	case 'M':
	  int pmd;
	  if (sscanf(OptArg, "%d", &pmd) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (pmd == 0)
	    PROGMODE = GAUSS;
	  else
	    PROGMODE = GAUSSFDR;
	  break;
	
	case 'E':
	  if (sscanf(OptArg, "%lf", &PROBA) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  FDR = PROBA;
	  break;
	  
	case 's':
	  if (sscanf(OptArg, "%lf", &GLEVEL) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
      PROBA = 2. * (1 - Utils<double>::cumNormal(GLEVEL));
	  FDR = PROBA;
	  break;

	case 'c':
	  int dep;
	  AUTOFDR = false;
	  if (sscanf(OptArg, "%d", &dep) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (dep != 0) FDRINDEP = false;
	  else FDRINDEP = true;
	  break;
	  
	case 'n':
	  if (sscanf(OptArg, "%d", &NSCALEXY) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;
	
	case 'N':
	  if (sscanf(OptArg, "%d", &NSCALEZ) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;

	case 'F':
	  if (sscanf(OptArg, "%d", &FSCALEXY) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;

	case 'f':
	  if (sscanf(OptArg, "%d", &FSCALEZ) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;

	case 'I':
	  int imd;
	  if (sscanf(OptArg, "%d", &imd) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (imd == 0)
	    ITERMODE = DIRECT;
	  else ITERMODE = HSD;
	  break;
	
	case 'i':
	  if (sscanf(OptArg, "%d", &NITER) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;
	
	case 'B':
	  int tbrd;
	  if (sscanf(OptArg, "%d", &tbrd) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (tbrd == 0)
	    T_BORDER = I_MIRROR;
	  else
	    T_BORDER = I_PERIOD;
	  break;

	case 'R':
	  BONF = true;
	  break;

	case 'v': 
	  VERBOSE = true;
	  break;
	
	case '?':
	default: 
	  usage(argv); 
	  exit(-1);
	}
    } 
       

    // get optional input file names from trailing parameters and open files 
    if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
    else { usage(argv); exit(-1); }

    if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
    else { usage(argv); exit(-1); }

    // make sure there are not too many parameters 
    if (OptInd < argc)
      {
	fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
	exit(-1);
      }
}

// max. scale of the data
int scaleOfData (int dim, int nx, int ny, int nz)
{
	int s = 0;
	
	if (dim == 1)
	{
		s = iilog2(nx);
	}
	else if (dim == 2)
	{
		s = MIN(iilog2(nx), iilog2(ny));
	}
	else // dim == 3
	{
		s = MIN(iilog2(nx), iilog2(ny));
		s = MIN(s, iilog2(nz));
	}
	return s;
}

// display the content of a data array (useful in debug mode)
template <typename DATATYPE>
void display (to_array<DATATYPE, true> &data)
{
  if (VERBOSE)
    {
      for (int z=0; z<data.nz(); z++)
	{
	  for (int y=0; y<data.ny(); y++)
	    {
	      for (int x=0; x<data.nx(); x++)
		cout << "(" << x << "," << y << "," << z << ")=" << data(x, y, z) << " ";
	      cout << endl << endl;
	    }
	  cout << endl << endl;
	}
    }
  cout << " data information : " << endl;
  data.info();
  cout << endl << endl;
}

// extend the data with specified border condition; the extension size at each border is returned;
// maxlen is the max. length of the filters
// this function garantees the even number of samples along each direction when decimated transform is used.
template <typename DATATYPE>
void extData(to_array<DATATYPE, true> &data, int maxlen[], type_border BORDERTYPE, int scale, bool dec[], int ext[6])
{
  int nx = data.nx(), ny = data.ny(), nz = data.nz();
  int dim = 3;
  for (int k=0; k<6; k++) ext[k] = 0;

  if ((BORDERTYPE != I_PERIOD) && ( (NITER == 0) || ((NITER != 0) && (MS == 0)) ))
  {  
  for (int k=0; k<dim; k++)
  {
    if (dec[k])
      {
	ext[2*k] = (maxlen[k] + 1) / 2;
	ext[2*k+1] = ext[2*k];
      }
    else
      {
	ext[2*k] = ((maxlen[k] + 1) / 2) * POW2(scale-1);
	ext[2*k+1] = ext[2*k];
      }
  }
  }

  int lext=ext[0], rext=ext[1], uext=ext[2], dext=ext[3], bext=ext[4], fext=ext[5];

  int newlen1 = nx + lext + rext, newlen2 = ny + uext + dext, newlen3 = nz + bext + fext;
  if ((newlen1 % 2 != 0) && dec[0]) { newlen1++; ext[1]++; rext++; }
  if ((newlen2 % 2 != 0) && dec[1]) { newlen2++; ext[3]++; dext++; }
  if ((newlen3 % 2 != 0) && dec[2]) { newlen3++; ext[5]++; fext++; }
      
  to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen1, newlen2, newlen3);
  for (int z=0; z<newlen3; z++)
    for (int y=0; y<newlen2; y++)
      for (int x=0; x<newlen1; x++)
	(*temp)(x, y, z) = data(x-lext, y-uext, z-bext, BORDERTYPE);
  data = *temp;
  delete temp; temp = NULL;
}

// extract the data from the extension version
template <typename DATATYPE>
void extractData (to_array<DATATYPE, true> &data, int ext[6])
{
  int lext = ext[0], rext = ext[1], uext = ext[2];
  int dext = ext[3], bext = ext[4], fext = ext[5];
  
  if ((lext != 0) || (rext != 0) || (uext != 0) || (dext != 0) || (bext != 0) || (fext != 0))
  {
      int finallen1 = data.nx() - lext - rext;
      int finallen2 = data.ny() - uext - dext;
      int finallen3 = data.nz() - bext - fext;
      to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen1, finallen2, finallen3);
      for (int x=0; x<finallen1; x++) 
	for (int y=0; y<finallen2; y++)
	  for (int z=0; z<finallen3; z++)	
	    (*tdata)(x, y, z) = data(lext+x, uext+y, bext+z);
      data = *tdata;
      delete tdata; tdata = NULL;
  }	  
}

template <typename DATATYPE>
void setConst (to_array<DATATYPE, true> &dest, to_array<DATATYPE, true> &ref, double cst)
{
	int dim = ref.naxis();
	int nx = ref.nx(), ny = ref.ny(), nz = ref.nz();
	if (dim == 1)
	{
		dest.resize(nx);
		for (int x=0; x<nx; x++) dest(x) = (DATATYPE)cst;
	}
	else if (dim == 2)
	{
		dest.resize(nx, ny);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
			dest(x, y) = (DATATYPE)cst;
	}
	else // dim == 3
	{
		dest.resize(nx, ny, nz);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
			dest(x, y, z) = (DATATYPE)cst;		
	}
}

// standard deviation estimation
double noiseEstimation (fltarray &data)
{
  int dim = data.naxis();
  double estimSigma;
  WaveletShrinkage<float, float> ws;

  if (dim == 1)
    {
      fltarray* ch = new fltarray;
      fltarray* cg = new fltarray;
      OrthogonalWaveletTransform<float>::dwt1D(data, OrthogonalWaveletTransform<float>::DAUB4, *ch, *cg);
      estimSigma = ws.gaussianStdDevEstim(*cg);
      delete ch; delete cg;
      ch = NULL; cg = NULL;
    }
  else if (dim == 2)
    {
      fltarray* chh = new fltarray;
      fltarray* chg = new fltarray;
      fltarray* cgh = new fltarray;
      fltarray* cgg = new fltarray;
      OrthogonalWaveletTransform<float>::dwt2D(data, OrthogonalWaveletTransform<float>::DAUB4, *chh, *chg, *cgh, *cgg);
      estimSigma = ws.gaussianStdDevEstim(*cgg);
      delete chh; delete chg; delete cgh; delete cgg;
      chh = NULL; chg = NULL; cgh = NULL; cgg = NULL;
    }
  else // dim == 3
    {
      fltarray* chhh = new fltarray;
      fltarray* chhg = new fltarray;
      fltarray* chgh = new fltarray;
      fltarray* chgg = new fltarray;
      fltarray* cghh = new fltarray;
      fltarray* cghg = new fltarray;
      fltarray* cggh = new fltarray;
      fltarray* cggg = new fltarray;
      OrthogonalWaveletTransform<float>::dwt3D(data, OrthogonalWaveletTransform<float>::DAUB4, *chhh, *chhg, *chgh, *chgg, *cghh, *cghg, *cggh, *cggg);
      estimSigma = ws.gaussianStdDevEstim(*cggg);
      delete chhh; delete chhg; delete chgh; delete chgg;
      delete cghh; delete cghg; delete cggh; delete cggg;
      chhh = NULL; chhg = NULL; chgh = NULL; chgg = NULL;
      cghh = NULL; cghg = NULL; cggh = NULL; cggg = NULL;
    }
  
  return estimSigma;
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void waveletDenoise (fltarray &data, SubBandFilter sbf[], bool dec[], to_array<SUPTYPE, true> *multiSup = NULL)
{
	int extxy[NSCALEXY][6];
	int extz[NSCALEZ][6];
	double fdrp;
	int N = data.n_elem();
	bool DEFAULT = (PROBA == 0);

	type_border TB = sbf[0].Border;
	int h0lenx = 0, h1lenx = 0, g0lenx = 0, g1lenx = 0;
	int h0leny = 0, h1leny = 0, g0leny = 0, g1leny = 0;
	int h0lenz = 0, h1lenz = 0, g0lenz = 0, g1lenz = 0;
	sbf[0].getFilterSize(h0lenx, h1lenx, g0lenx, g1lenx);
	sbf[1].getFilterSize(h0leny, h1leny, g0leny, g1leny);
	sbf[2].getFilterSize(h0lenz, h1lenz, g0lenz, g1lenz);
	int maxlenxy[3] = {MAX(h0lenx, MAX(h1lenx, MAX(g0lenx, g1lenx))), \
			   MAX(h0leny, MAX(h1leny, MAX(g0leny, g1leny))), 0};
	int maxlenz [3] = {0, 0, MAX(h0lenz, MAX(h1lenz, MAX(g0lenz, g1lenz)))};
	
	// wavelet filter configuration
	MWIRWaveletTransform *mwirWT = new MWIRWaveletTransform(sbf[0], sbf[1], sbf[2]);
	WaveletShrinkage<float, SUPTYPE> ws;


	// denoising
	if (VERBOSE) cerr << "Wavelet transform 2D plane-by-plane ... " << endl;

	fltarray *chh = new fltarray;
	fltarray *chg = new fltarray[NSCALEXY];
	fltarray *cgh = new fltarray[NSCALEXY];
	fltarray *cgg = new fltarray[NSCALEXY];
	for (int sxy=1; sxy<=NSCALEXY; sxy++)
	  {
	    extData(data, maxlenxy, TB, sxy, dec, extxy[sxy-1]); 

	    mwirWT->transform3DXY(data, *chh, chg[sxy-1], cgh[sxy-1], cgg[sxy-1], sxy, dec);
	    data = *chh;
	  }

	if (VERBOSE) cerr << "Transform 2D plane-by-plane complete." \
			  << endl << "Transform 1D along Z and thresholding ... " << endl;

	int msi = 0;
	fltarray *ch1 = new fltarray; fltarray *ch2 = new fltarray;
	fltarray *ch3 = new fltarray; fltarray *ch4 = new fltarray;
	fltarray *cg1 = new fltarray[NSCALEZ];
	fltarray *cg2 = new fltarray[NSCALEZ];
	fltarray *cg3 = new fltarray[NSCALEZ];
	fltarray *cg4 = new fltarray[NSCALEZ];

	for (int sxy=NSCALEXY; sxy>=1; sxy--)
	  {
	    for (int sz=1; sz<=NSCALEZ; sz++)
	      {
		extData(*chh, maxlenz, TB, sz, dec, extz[sz-1]); 
		extData(chg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		extData(cgh[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		extData(cgg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		
		mwirWT->transform3DZ(*chh, *ch1, cg1[sz-1], sz, dec[2]);
		mwirWT->transform3DZ(chg[sxy-1], *ch2, cg2[sz-1], sz, dec[2]);
		mwirWT->transform3DZ(cgh[sxy-1], *ch3, cg3[sz-1], sz, dec[2]);
		mwirWT->transform3DZ(cgg[sxy-1], *ch4, cg4[sz-1], sz, dec[2]);

		*chh = *ch1;
		chg[sxy-1] = *ch2;
		cgh[sxy-1] = *ch3;
		cgg[sxy-1] = *ch4;
	      }
	    
	    for (int sz=NSCALEZ; sz>=1; sz--)
	      {
		double dlen = ch1->n_elem();
		double pr = (BONF ? PROBA/dlen : PROBA);
		
        if ((sxy < FSCALEXY) && (sz < FSCALEZ))
        {
			if (sxy == NSCALEXY)
			{
               setConst(cg1[sz-1], cg1[sz-1], 0);
               if (multiSup != NULL) 
                  setConst(multiSup[msi++], cg1[sz-1], -Sigma);  
            }
            setConst(cg2[sz-1], cg2[sz-1], 0);
            setConst(cg3[sz-1], cg3[sz-1], 0);
            setConst(cg4[sz-1], cg4[sz-1], 0);
            if (multiSup != NULL)
            {
              setConst(multiSup[msi++], cg2[sz-1], -Sigma);
              setConst(multiSup[msi++], cg3[sz-1], -Sigma);
              setConst(multiSup[msi++], cg4[sz-1], -Sigma);
            }
			if (sz == NSCALEZ)
			  {
                setConst(*ch2, *ch2, 0);
                setConst(*ch3, *ch3, 0);
                setConst(*ch4, *ch4, 0);
                if (multiSup != NULL)
                {
                   setConst(multiSup[msi++], *ch2, -Sigma);
                   setConst(multiSup[msi++], *ch3, -Sigma);
                   setConst(multiSup[msi++], *ch4, -Sigma);
                }
              }
        }
		else if (PROGMODE == GAUSS)
		  {
    		int sc[3] = {NSCALEXY, NSCALEXY, NSCALEZ};
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;

		    if (multiSup != NULL)
		      {
            if (sxy == NSCALEXY)
			  pr = ws.gaussHardThreshold(cg1[sz-1], sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			pr = ws.gaussHardThreshold(cg2[sz-1], sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			pr = ws.gaussHardThreshold(cg3[sz-1], sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			pr = ws.gaussHardThreshold(cg4[sz-1], sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.gaussHardThreshold(*ch2, sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch3, sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch4, sc, pr, Sigma, &multiSup[msi++], N, DEFAULT);
			  }
		      }
		    else
		      {
            if (sxy == NSCALEXY)
		        pr = ws.gaussHardThreshold(cg1[sz-1], sc, pr, Sigma, NULL, N, DEFAULT);
			pr = ws.gaussHardThreshold(cg2[sz-1], sc, pr, Sigma, NULL, N, DEFAULT);
			pr = ws.gaussHardThreshold(cg3[sz-1], sc, pr, Sigma, NULL, N, DEFAULT);
			pr = ws.gaussHardThreshold(cg4[sz-1], sc, pr, Sigma, NULL, N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.gaussHardThreshold(*ch2, sc, pr, Sigma, NULL, N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch3, sc, pr, Sigma, NULL, N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch4, sc, pr, Sigma, NULL, N, DEFAULT);
			  }
		      }
		  } // if PROGMODE == GAUSS
		else if (PROGMODE == GAUSSFDR)
		  {
		    if (multiSup != NULL)
		      {
            if (sxy == NSCALEXY)
            {
			  fdrp = ws.gaussFDRThreshold(cg1[sz-1], Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			  if (VERBOSE)
			    cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
            }
			fdrp = ws.gaussFDRThreshold(cg2[sz-1], Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.gaussFDRThreshold(cg3[sz-1], Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.gaussFDRThreshold(cg4[sz-1], Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			if (sz == NSCALEZ)
			  {
			    fdrp = ws.gaussFDRThreshold(*ch2, Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.gaussFDRThreshold(*ch3, Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.gaussFDRThreshold(*ch4, Sigma, FDR, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			  }
		      }
		    else			  
		      {
            if (sxy == NSCALEXY)
            {
			fdrp = ws.gaussFDRThreshold(cg1[sz-1], Sigma, FDR, FDRINDEP);
			if (VERBOSE)
			  cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
            }
			fdrp = ws.gaussFDRThreshold(cg2[sz-1], Sigma, FDR, FDRINDEP);
			if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.gaussFDRThreshold(cg3[sz-1], Sigma, FDR, FDRINDEP);
			if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.gaussFDRThreshold(cg4[sz-1], Sigma, FDR, FDRINDEP);
			if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			if (sz == NSCALEZ)
			  {
			    fdrp = ws.gaussFDRThreshold(*ch2, Sigma, FDR, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.gaussFDRThreshold(*ch3, Sigma, FDR, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.gaussFDRThreshold(*ch4, Sigma, FDR, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			  }
		      } 
		  } // if PROGMODE == GAUSSFDR

		*chh = *ch1;
		chg[sxy-1] = *ch2;
		cgh[sxy-1] = *ch3;
		cgg[sxy-1] = *ch4;
		mwirWT->reconstruction3DZ(*chh, cg1[sz-1], *ch1, sz, dec[2]);
		mwirWT->reconstruction3DZ(chg[sxy-1], cg2[sz-1], *ch2, sz, dec[2]);
		mwirWT->reconstruction3DZ(cgh[sxy-1], cg3[sz-1], *ch3, sz, dec[2]);
		mwirWT->reconstruction3DZ(cgg[sxy-1], cg4[sz-1], *ch4, sz, dec[2]);
		extractData(*ch1, extz[sz-1]);
		extractData(*ch2, extz[sz-1]);
		extractData(*ch3, extz[sz-1]);
		extractData(*ch4, extz[sz-1]);
	      } // for sz
	    mwirWT->reconstruction3DXY(*ch1, *ch2, *ch3, *ch4, *chh, sxy, dec);
	    extractData(*chh, extxy[sxy-1]);
	  } // for sxy
	data = *chh;

	if (VERBOSE) cerr << "Thresholding complete." << endl;
	
	delete chh; chh = NULL; delete [] chg; chg = NULL; 
	delete [] cgh; cgh = NULL; delete [] cgg; cgg = NULL;
	delete ch1; delete ch2; delete ch3; delete ch4;
	ch1 = ch2 = ch3 = ch4 = NULL;
	delete [] cg1; cg1 = NULL; delete [] cg2; cg2 = NULL;
	delete [] cg3; cg3 = NULL; delete [] cg4; cg4 = NULL;
	delete mwirWT; mwirWT = NULL;
}

template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef, double lambda)
{
  int nx = data.nx(), ny = data.ny(), nz = data.nz();

  for (int x=0; x<nx; x++)
    for (int y=0; y<ny; y++)
      for (int z=0; z<nz; z++)
	{
	  data(x, y, z) = (coef(x, y, z) >= 0) ? origdata(x, y, z) : data(x, y, z);
	  if (ITERMODE == HSD)
	    { 
	      data(x, y, z) = soft_threshold(data(x, y, z), lambda);
	    }
	}
}

// put to zeros the negative values
void posProject (fltarray &data)
{	
	float value;
	int len = data.n_elem();
	
	for (int i=0; i<len; i++)
	{
		value = data(i);
		data(i) = MAX(value, 0);
	}
}

// solve with the iterative algo. in using the multiresolution support
void multiSupIter (fltarray &origdata, fltarray &solution, SubBandFilter sbf[], bool dec[], fltarray *multiSup, int niter)
{
  int extxy[NSCALEXY][6];
  int extz[NSCALEZ][6];

  type_border TB = sbf[0].Border;
  int h0lenx = 0, h1lenx = 0, g0lenx = 0, g1lenx = 0;
  int h0leny = 0, h1leny = 0, g0leny = 0, g1leny = 0;
  int h0lenz = 0, h1lenz = 0, g0lenz = 0, g1lenz = 0;
  sbf[0].getFilterSize(h0lenx, h1lenx, g0lenx, g1lenx);
  sbf[1].getFilterSize(h0leny, h1leny, g0leny, g1leny);
  sbf[2].getFilterSize(h0lenz, h1lenz, g0lenz, g1lenz);
  int maxlenxy[3] = {MAX(h0lenx, MAX(h1lenx, MAX(g0lenx, g1lenx))), \
		     MAX(h0leny, MAX(h1leny, MAX(g0leny, g1leny))), 0};
  int maxlenz [3] = {0, 0, MAX(h0lenz, MAX(h1lenz, MAX(g0lenz, g1lenz)))};
  
  // wavelet filter configuration
  MWIRWaveletTransform *mwirWT = new MWIRWaveletTransform(sbf[0], sbf[1], sbf[2]);
 
  double lambda = 0, delta = 0;
  if (ITERMODE == HSD)
    {
      if (niter <= 1) return;
      lambda = 1.;
      delta = lambda / (niter-1);
    }

  fltarray *tempd = new fltarray;
  fltarray *chh = new fltarray,           *chhd = new fltarray[NSCALEXY];
  fltarray *chg = new fltarray[NSCALEXY], *chgd = new fltarray[NSCALEXY];
  fltarray *cgh = new fltarray[NSCALEXY], *cghd = new fltarray[NSCALEXY];
  fltarray *cgg = new fltarray[NSCALEXY], *cggd = new fltarray[NSCALEXY];
  fltarray *ch1 = new fltarray, *ch1d = new fltarray, *ch2 = new fltarray, *ch2d = new fltarray;
  fltarray *ch3 = new fltarray, *ch3d = new fltarray, *ch4 = new fltarray, *ch4d = new fltarray;
  fltarray *cg1 = new fltarray[NSCALEZ], *cg1d = new fltarray[NSCALEZ];
  fltarray *cg2 = new fltarray[NSCALEZ], *cg2d = new fltarray[NSCALEZ];
  fltarray *cg3 = new fltarray[NSCALEZ], *cg3d = new fltarray[NSCALEZ];
  fltarray *cg4 = new fltarray[NSCALEZ], *cg4d = new fltarray[NSCALEZ];
  for (int i=0; i<niter; i++)
  {
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;

      *tempd = origdata;
      for (int sxy=1; sxy<=NSCALEXY; sxy++)
	{
	  extData(*tempd, maxlenxy, TB, sxy, dec, extxy[sxy-1]); 
	  extData(solution, maxlenxy, TB, sxy, dec, extxy[sxy-1]); 
	  
	  mwirWT->transform3DXY(solution, *chh, chg[sxy-1], cgh[sxy-1], cgg[sxy-1], sxy, dec);
	  mwirWT->transform3DXY(*tempd, chhd[sxy-1], chgd[sxy-1], cghd[sxy-1], cggd[sxy-1], sxy, dec);
	  solution = *chh;
	  *tempd = chhd[sxy-1];
	}

      int msi = 0;      
      for (int sxy=NSCALEXY; sxy>=1; sxy--)
	{
	  for (int sz=1; sz<=NSCALEZ; sz++)
	    {
	      extData(*chh, maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(chhd[sxy-1],maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(chg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(chgd[sxy-1],maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(cgh[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(cghd[sxy-1],maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(cgg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
	      extData(cggd[sxy-1],maxlenz, TB, sz, dec, extz[sz-1]);

	      mwirWT->transform3DZ(*chh, *ch1, cg1[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(chg[sxy-1], *ch2, cg2[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(cgh[sxy-1], *ch3, cg3[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(cgg[sxy-1], *ch4, cg4[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(chhd[sxy-1], *ch1d, cg1d[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(chgd[sxy-1], *ch2d, cg2d[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(cghd[sxy-1], *ch3d, cg3d[sz-1], sz, dec[2]);
	      mwirWT->transform3DZ(cggd[sxy-1], *ch4d, cg4d[sz-1], sz, dec[2]);

	      *chh = *ch1;
	      chg[sxy-1] = *ch2;
	      cgh[sxy-1] = *ch3;
	      cgg[sxy-1] = *ch4;
	      chhd[sxy-1] = *ch1d;
	      chgd[sxy-1] = *ch2d;
	      cghd[sxy-1] = *ch3d;
	      cggd[sxy-1] = *ch4d;
	    }

	  if (sxy == NSCALEXY) *ch1 = *ch1d;
	  for (int sz=NSCALEZ; sz>=1; sz--)
	    {
          if (sxy == NSCALEXY)
	        procArr<float>(cg1d[sz-1], cg1[sz-1], multiSup[msi++], lambda);
	      procArr<float>(cg2d[sz-1], cg2[sz-1], multiSup[msi++], lambda);
	      procArr<float>(cg3d[sz-1], cg3[sz-1], multiSup[msi++], lambda);
	      procArr<float>(cg4d[sz-1], cg4[sz-1], multiSup[msi++], lambda);
	      if (sz == NSCALEZ)
		{
		  procArr<float>(*ch2d, *ch2, multiSup[msi++], lambda);
		  procArr<float>(*ch3d, *ch3, multiSup[msi++], lambda);
		  procArr<float>(*ch4d, *ch4, multiSup[msi++], lambda);
		}
	  
	      *chh = *ch1;
	      chg[sxy-1] = *ch2;
	      cgh[sxy-1] = *ch3;
	      cgg[sxy-1] = *ch4;
	      mwirWT->reconstruction3DZ(*chh, cg1[sz-1], *ch1, sz, dec[2]);
	      mwirWT->reconstruction3DZ(chg[sxy-1], cg2[sz-1], *ch2, sz, dec[2]);
	      mwirWT->reconstruction3DZ(cgh[sxy-1], cg3[sz-1], *ch3, sz, dec[2]);
	      mwirWT->reconstruction3DZ(cgg[sxy-1], cg4[sz-1], *ch4, sz, dec[2]);
	      extractData(*ch1, extz[sz-1]);
	      extractData(*ch2, extz[sz-1]);
	      extractData(*ch3, extz[sz-1]);
	      extractData(*ch4, extz[sz-1]);
	    } // for sz
	  mwirWT->reconstruction3DXY(*ch1, *ch2, *ch3, *ch4, *chh, sxy, dec);
	  extractData(*chh, extxy[sxy-1]);
	} // for sxy
      
      solution = *chh;
      posProject(solution);
      lambda -= delta;
      
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " complete. " << endl;
  }
  
  delete tempd; tempd = NULL;
  delete chh; chh = NULL; delete [] chhd; chhd = NULL;
  delete [] chg; chg = NULL; delete [] chgd; chgd = NULL;
  delete [] cgh; cgh = NULL; delete [] cghd; cghd = NULL;
  delete [] cgg; cgg = NULL; delete [] cggd; cggd = NULL;
  delete  ch1; ch1 = NULL; delete  ch1d; ch1d = NULL;
  delete  ch2; ch2 = NULL; delete  ch2d; ch2d = NULL;
  delete  ch3; ch3 = NULL; delete  ch3d; ch3d = NULL;
  delete  ch4; ch4 = NULL; delete  ch4d; ch4d = NULL;
  delete [] cg1; cg1 = NULL; delete [] cg1d; cg1d = NULL;
  delete [] cg2; cg2 = NULL; delete [] cg2d; cg2d = NULL;
  delete [] cg3; cg3 = NULL; delete [] cg3d; cg3d = NULL;
  delete [] cg4; cg4 = NULL; delete [] cg4d; cg4d = NULL;
  delete mwirWT; mwirWT = NULL;
}

// wavelet direct denoising
void waveletGaussDenoise (fltarray &data)
{	
  // backup the original data
  fltarray *origData = new fltarray;
  *origData = data;
	
  // wavelet filter configuration
  SubBandFilter sbh = SubBandFilter(F_HAAR, NORM_L2);
  SubBandFilter sbd = SubBandFilter(F_DAUBE_8, NORM_L2);
  sbh.Border = sbd.Border = T_BORDER;

  SubBandFilter sbf[3] = { sbh, sbh, sbh };
  SubBandFilter sbf1[3] = { sbd, sbd, sbd };  
  
  fltarray *multiSup = NULL;
  if (NITER > 1)
    {
      multiSup = new fltarray[(NSCALEZ-1)*4+7 + ((NSCALEZ-1)*3+6)*(NSCALEXY-1)];
    }
  if (VERBOSE)
    cerr << "Initial denoising ... " << endl;
  if (Sigma < 0)
    {
      Sigma = noiseEstimation(data);
      if (VERBOSE) cerr << "Estimated sigma = " << Sigma << endl;
    }
  if (WAVE == 0)
    waveletDenoise<float>(data, sbf, DEC, multiSup);
  else
    waveletDenoise<float>(data, sbf1, DEC, multiSup);
  posProject(data);
  if (VERBOSE)
    cerr << "Entering into the iterative denoising ..." << endl;
  if (NITER > 1)
    {
      SubBandFilter sbo = SubBandFilter(F_ODEGARD_9_7);
      sbo.Border = T_BORDER;
      SubBandFilter mssbf[3] = { sbo, sbo, sbo };
      if (MS != 0)
	multiSupIter (*origData, data, mssbf, DEC, multiSup, NITER-1);
      else if (WAVE == 0)
	multiSupIter (*origData, data, sbf, DEC, multiSup, NITER-1);
      else
	multiSupIter (*origData, data, sbf1, DEC, multiSup, NITER-1);
    }
  if (VERBOSE)
    cerr << "Iteration complete." << endl;
  
  delete origData; origData = NULL;
  if (multiSup != NULL) { delete [] multiSup; multiSup = NULL; }
}

//*********************************************************************
int main(int argc, char *argv[])
{
   char cmd[512];
   cmd[0] = '\0';
   for (int k=0; k<argc; k++) 
   	sprintf(cmd, "%s %s", cmd, argv[k]);
    
   fitsstruct header;
   fltarray *data = new fltarray;
   
   // Get command line arguments, open input file(s) if necessary
   filtinit(argc, argv);
   fits_read_fltarr(Name_Imag_In, *data, &header);

   header.bitpix = BP_FLOAT;
   header.origin = cmd;

   if (data->naxis() != 3)
     {
       cerr << "Image is not 3D" << endl;
       delete data; data = NULL;
       exit(-1);
     }
   
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALEXY = MIN(NSCALEXY, scaleOfData(2, nx, ny, nz));
   NSCALEZ = MIN(NSCALEZ, scaleOfData(1, nz, nz, nz));
   if (AUTOFDR)
     {
       if (DEC[0] && DEC[1] && DEC[2]) FDRINDEP = true;
       else FDRINDEP = false;
     }
   if (PROGMODE == GAUSSFDR) BONF = false;

   if (VERBOSE)
   {
   	cout << "Input image : " << Name_Imag_In << endl;
   	cout << "Output image : " << Name_Imag_Out << endl;

	if (WAVE == 1)
	  cout << "Denoising wavelet : Daubechies - 4 " << endl;
	else
	  cout << "Denoising wavelet : Haar" << endl;

	if (DEC[0]) cout << " X - decimated ";
	else cout << " X - undecimated ";
	if (DEC[1]) cout << " Y - decimated ";
	else cout << " Y - undecimated ";
	if (DEC[2]) cout << " Z - decimated ";
	else cout << " Z - undecimated ";
	cout << endl;
	
	if (PROGMODE == GAUSS)
	  cout << "Mode : " << "k-sigma thresholding" << endl;
	else if (PROGMODE == GAUSSFDR)
	  {
	    cout << "Mode : " << "FDR thresholding" << endl;
	    cout << "FDR mode : ";
	    if (FDRINDEP) cout << "independence mode";
	    else cout << "dependence mode";
	    cout << endl;	    
	  }
	
	if (Sigma < 0) cout << "sigma will be estimated" << endl;
	else cout << "sigma = " << Sigma << endl;
	
	if (ITERMODE == DIRECT)
	  cout << "Iterative Mode : Direct" << endl;
	else
	  cout << "Iterative Mode : HSD" << endl;
	cout << "Max. Iteration : " << NITER << endl;
   	
	cout << "Max scale(s) X and Y: " << NSCALEXY << endl;
	cout << "Max scale(s) Z: " << NSCALEZ << endl;
    cout << "First detection X-Y scale : " << FSCALEXY << endl;
    cout << "First detection Z scale : " << FSCALEZ << endl;
	if (PROGMODE != GAUSSFDR)
	{
	  if (PROBA != 0)
	    cout << "Proba. : " << PROBA << endl;
	  else
	    cout << "Proba. : Universal deciding" << endl;
	}
	else
	{
	  cout << "FDR : " << FDR << endl;
	}
	
	cout << "Bonferroni's schema : " << (BONF ? "yes" : "no") << endl;
	if (MS == 0)
	  cout << "Iterative reconstruction wavelet : same as denoising case" << endl;
	else
	  cout << "Iterative reconstruction wavelet : Odegard 9/7" << endl;
	cout << "Border : " << (T_BORDER == I_MIRROR ? "symmetric" : "periodic") << endl;
   }
   
   // Denoising
   try{
     if (VERBOSE)
	 cerr << "Denoising ..." << endl;
     
     waveletGaussDenoise(*data);
     if (VERBOSE)
       cerr << "Writing denoising result file ... " << endl;
     fits_write_fltarr(Name_Imag_Out, *data, &header);
     if (VERBOSE)
       cerr << "Writing complete." << endl;
       
     if (VERBOSE)
	 cerr << "Denoising complete." << endl;
     
     delete data; data = NULL;
   }
   catch (MWIRException mwirExcept)
   {
   	cerr << "Exception : " << endl;
   	cerr << mwirExcept.getReason() << endl;
   	exit (-1);
   }

    return 0;
}
