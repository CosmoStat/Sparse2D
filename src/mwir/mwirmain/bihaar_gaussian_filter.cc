/**************************************************
 * PROJECT : MultiWavelenImageResto (MWIR) - version I for Gaussian case
 * CEA
 * Filename : maing1.cc
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
int    NSCALE = 3;          // max. number of scales
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, true}; // X, Y - undecimated; Z - decimated
bool   AUTOFDR = true;      // automatically choose the FDR mode
bool   FDRINDEP = false;    // independence parameter in FDR
bool   GENER = false;       // generating an image with additive Gaussian noise
double Sigma = -1.;         // noise standard deviation - decide automatically or 1 in generation
type_border T_BORDER = I_MIRROR; // border type
long int SEED = time(NULL); // random seed
bool BONF = false;          // Bonferroni's multi-test schema
int  WAVE = 0;              // Wavelet used in the denoising [0, 1]
int  MS = 0;                // Wavelet used in the iterative restoration [0, 1]
int    FSCALE = 1;          // First detection scale
double GLEVEL = 0;          // Gauss-detection level (k-sigma)

// NOTICE : Bonferroni and FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-T value] wavelet used in the denoising" << endl;
  cout << "         T = 0 : Haar (default)" << endl;
  cout << "         T = 1 : Daubechies - 4" << endl;     
  cout << "    [-S value]" << endl;
  cout << "         S = 0 : use the same wavelet as in the denoising to the iterative restoration (default)" << endl;
  cout << "         S = 1 : use ODEGARD 9/7 wavelet in the iterative restoration" << endl;
  cout << "    [-G] generate an image from input image with additive Gaussian noise" << endl;
  cout << "    [-g value] standard deviation used in generation (default = 1) or in denoising (default = decide automatically) " << endl;
  cout << "    [-e value] generating integer seed (default = current time)" << endl;
  cout << "    [-D] to decimate X and Y directions in the wavelet transform (not default)" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : denoise (default)" << endl;
  cout << "         M = 1 : denoise with FDR" << endl;
  cout << "    [-E value] two sided Pb or FDR (default = automatically estimated)" << endl;
  cout << "    [-s value] Gauss-detection level (default = automatically decided)" << endl;
  cout << "    [-c value] manual FDR setting, c = 0/1 means independence/dependence of coef. (default = decide automatically)" << endl;
  cout << "    [-n value] max. number of scales (default = 3)" << endl;
  cout << "    [-F value] first detection scale (default = 1)" << endl;
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
    while ((c = GetOpt(argc,argv,"T:S:Gg:e:DM:E:s:c:n:F:I:i:B:Rv")) != -1) 
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

	case 'G':
	  GENER = true;
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
	
	case 'e':
	  if (sscanf(OptArg, "%ld", &SEED) != 1)
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
	  if (sscanf(OptArg, "%d", &NSCALE) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;
	
	case 'F':
	  if (sscanf(OptArg, "%d", &FSCALE) != 1)
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
	else
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
	int dim = data.naxis();
	if (VERBOSE)
	  {
	    if (dim == 1)
	      {
		for (int x=0; x<data.nx(); x++)
		  cout << "(" << x << ")=" << data(x) << " ";
		cout << endl << endl;
	      }
	    else if (dim == 2)
	      {
		for (int y=0; y<data.ny(); y++)
		  {
		    for (int x=0; x<data.nx(); x++)
		      cout << "(" << x << "," << y << ")=" << data(x, y) << " ";
		    cout << endl << endl;
		  }
	      }
	    else
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
	  }
	cout << " data information : " << endl;
	data.info();
	cout << endl << endl;
}

// extend the data with specified border condition; the extension size at each bord is returned;
// maxlen is the max. length of the filters
// this function garantees the even number of samples along each direction when decimated transform is used.
template <typename DATATYPE>
void extData(to_array<DATATYPE, true> &data, int maxlen[], type_border BORDERTYPE, int scale, bool dec[], int ext[6])
{
  int nx = data.nx(), ny = data.ny(), nz = data.nz();
  int dim = data.naxis();
  for (int k=0; k<6; k++) ext[k] = 0;

  // WARNING : 
  // no filter-length compensated extension when different wavelet 
  // from that of denoising is used in the iterative restoration
  // So in this case, the perfect-border-reconstruction property 
  // may not be satisfied

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
  if (dim == 1)
    {
      int newlen = nx + lext + rext;
      if ((newlen % 2 != 0) && dec[0]) { newlen++; ext[1]++; rext++; }
      to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen);

      for (int x=0; x<newlen; x++)
	(*temp)(x) = data(x-lext, BORDERTYPE);
      data = *temp;
      delete temp; temp = NULL;
    }
  else if (dim == 2)
    {
      int newlen1 = nx + lext + rext, newlen2 = ny + uext + dext;
      if ((newlen1 % 2 != 0) && dec[0]) { newlen1++; ext[1]++; rext++; }
      if ((newlen2 % 2 != 0) && dec[1]) { newlen2++; ext[3]++; dext++; }

      to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen1, newlen2);
      for (int y=0; y<newlen2; y++)
	for (int x=0; x<newlen1; x++)
	  (*temp)(x, y) = data(x-lext, y-uext, BORDERTYPE);
      data = *temp;
      delete temp; temp = NULL;
    }
  else // dim == 3
    {
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
}

// extract the data from the extension version
template <typename DATATYPE>
void extractData (to_array<DATATYPE, true> &data, int ext[6])
{
  int dim = data.naxis();
  int lext = ext[0], rext = ext[1], uext = ext[2];
  int dext = ext[3], bext = ext[4], fext = ext[5];
  
  if ((dim == 1) && ((lext != 0) || (rext != 0)))
    {
      int finallen = data.nx() - lext - rext;
      to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen);
      for (int x=0; x<finallen; x++) (*tdata)(x) = data(lext+x);
      data = (*tdata);
      delete tdata; tdata = NULL;
    }
  else if ((dim == 2) && \
	   ((lext != 0) || (rext != 0) || (uext != 0) || (dext != 0)))
    {
      int finallen1 = data.nx() - lext - rext;
      int finallen2 = data.ny() - uext - dext;
      to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen1, finallen2);
      for (int x=0; x<finallen1; x++) 
	for (int y=0; y<finallen2; y++)
	  (*tdata)(x, y) = data(lext+x, uext+y);
      data = *tdata;
      delete tdata; tdata = NULL;
    }
  else if ((dim == 3) && ((lext != 0) || (rext != 0) || \
			  (uext != 0) || (dext != 0) || \
			  (bext != 0) || (fext != 0)))
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
	int dim = data.naxis();
	int ext[NSCALE][6];
	double fdrp;
	int N = data.n_elem();
	bool DEFAULT = (PROBA == 0);

	type_border TB = sbf[0].Border;
	int h0len = 0, h1len = 0, g0len = 0, g1len = 0, maxlen[3] = {0, 0, 0};
	for (int k=0; k<dim; k++)
	{
	    sbf[k].getFilterSize(h0len, h1len, g0len, g1len);
	    maxlen[k] = MAX(h0len, MAX(h1len, MAX(g0len, g1len)));
	}

	// wavelet filter configuration
	MWIRWaveletTransform *mwirWT = new MWIRWaveletTransform(sbf[0], sbf[1], sbf[2]);
	WaveletShrinkage<float, SUPTYPE> ws;

	if (dim == 1)
	{
	        if (VERBOSE) cerr << "Wavelet transform 1D ... " << endl;
		fltarray *ch = new fltarray;
		fltarray *cg = new fltarray[NSCALE];
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform1D(data, *ch, cg[s-1], s, dec[0]);
			data = *ch;
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;
		
		for (int s=NSCALE; s>=1; s--)
		{
			double dlen = cg[s-1].n_elem();
			double pr = (BONF ? PROBA/dlen : PROBA);

            if (s < FSCALE)
            {
                setConst(cg[s-1], cg[s-1], 0);
                if (multiSup != NULL)
                    setConst(multiSup[s-1], cg[s-1], -Sigma);  
            }
			else if (PROGMODE == GAUSS)
			{
  	          int sc[1] = {NSCALE};
			  if (multiSup != NULL)
				pr = ws.gaussHardThreshold(cg[s-1], sc, pr, Sigma, &multiSup[s-1], N, DEFAULT);
			  else 	pr = ws.gaussHardThreshold(cg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == GAUSSFDR)
			{
			  if (multiSup != NULL)
			    {
			        fdrp = ws.gaussFDRThreshold(cg[s-1], Sigma, FDR, FDRINDEP, &multiSup[s-1]);
				if (VERBOSE)
				  cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else 	
			    {
			      fdrp = ws.gaussFDRThreshold(cg[s-1], Sigma, FDR, FDRINDEP);
			      if (VERBOSE)
				cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			}

			mwirWT->reconstruction1D(*ch, cg[s-1], data, s, dec[0]);
			extractData(data, ext[s-1]);
			*ch = data;
		}
		delete ch; delete [] cg; 
		ch = NULL; cg = NULL;
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	}
	else if (dim == 2)
	{
	        if (VERBOSE) cerr << "Wavelet transform 2D ... " << endl;
		fltarray *chh = new fltarray;
		fltarray *chg = new fltarray[NSCALE];
		fltarray *cgh = new fltarray[NSCALE];
		fltarray *cgg = new fltarray[NSCALE];
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform2D(data, *chh, chg[s-1], cgh[s-1], cgg[s-1], s, dec);
			data = *chh;
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;

		for (int s=NSCALE; s>=1; s--)
		{
		        double dlen = chg[s-1].n_elem();
			double pr = (BONF ? PROBA/dlen : PROBA);

            if (s < FSCALE)
            {
                setConst(chg[s-1], chg[s-1], 0);
                setConst(cgh[s-1], cgh[s-1], 0);
                setConst(cgg[s-1], cgg[s-1], 0);
                if (multiSup != NULL)
                {
                 setConst(multiSup[3*(s-1)], chg[s-1], -Sigma);  
                 setConst(multiSup[3*(s-1)+1], cgh[s-1], -Sigma);  
                 setConst(multiSup[3*(s-1)+2], cgg[s-1], -Sigma);  
                }
            }
			else if (PROGMODE == GAUSS)			
			{
              int sc[2] = {NSCALE, NSCALE};
			  if (multiSup != NULL)
			    {
				pr = ws.gaussHardThreshold(chg[s-1], sc, pr, Sigma, &multiSup[3*(s-1)], N, DEFAULT);
				pr = ws.gaussHardThreshold(cgh[s-1], sc, pr, Sigma, &multiSup[3*(s-1)+1], N, DEFAULT);
				pr = ws.gaussHardThreshold(cgg[s-1], sc, pr, Sigma, &multiSup[3*(s-1)+2], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.gaussHardThreshold(chg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cgh[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cgg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == GAUSSFDR)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.gaussFDRThreshold(chg[s-1], Sigma, FDR, FDRINDEP, &multiSup[3*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgh[s-1], Sigma, FDR, FDRINDEP, &multiSup[3*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgg[s-1], Sigma, FDR, FDRINDEP, &multiSup[3*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.gaussFDRThreshold(chg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgh[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			}

			mwirWT->reconstruction2D(*chh, chg[s-1], cgh[s-1], cgg[s-1], data, s, dec);
			extractData(data, ext[s-1]);
			*chh = data;
		}
		delete chh; delete [] chg; delete [] cgh; delete [] cgg;
		chh = NULL; chg = NULL; cgh = NULL; cgg = NULL;
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	}
	else // dim == 3
	{
	        if (VERBOSE) cerr << "Wavelet transform 3D ... " << endl;
		fltarray *chhh = new fltarray;
		fltarray *chhg = new fltarray[NSCALE];
		fltarray *chgh = new fltarray[NSCALE];
		fltarray *chgg = new fltarray[NSCALE];
		fltarray *cghh = new fltarray[NSCALE];
		fltarray *cghg = new fltarray[NSCALE];
		fltarray *cggh = new fltarray[NSCALE];
		fltarray *cggg = new fltarray[NSCALE];
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform3D(data, *chhh, chhg[s-1], chgh[s-1], chgg[s-1], \
					    cghh[s-1], cghg[s-1], cggh[s-1], cggg[s-1], s, dec);
			data = *chhh;
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;

		for (int s=NSCALE; s>=1; s--)
		{
			double dlen = chhg[s-1].n_elem();
			double pr = (BONF ? PROBA/dlen : PROBA);

            if (s < FSCALE)
            {
                setConst(chhg[s-1], chhg[s-1], 0);
                setConst(chgh[s-1], chgh[s-1], 0);
                setConst(chgg[s-1], chgg[s-1], 0);
                setConst(cghh[s-1], cghh[s-1], 0);
                setConst(cghg[s-1], cghg[s-1], 0);
                setConst(cggh[s-1], cggh[s-1], 0);
                setConst(cggg[s-1], cggg[s-1], 0);
                if (multiSup != NULL)
                {                
                   setConst(multiSup[7*(s-1)], chhg[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+1], chgh[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+2], chgg[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+3], cghh[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+4], cghg[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+5], cggh[s-1], -Sigma);  
                   setConst(multiSup[7*(s-1)+6], cggg[s-1], -Sigma);  
                }
            }
			else if (PROGMODE == GAUSS)
			{
		      int sc[3] = {NSCALE, NSCALE, NSCALE};                         
			  if (multiSup != NULL)
			    {
				pr = ws.gaussHardThreshold(chhg[s-1], sc, pr, Sigma, &multiSup[7*(s-1)], N, DEFAULT);
				pr = ws.gaussHardThreshold(chgh[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+1], N, DEFAULT);
				pr = ws.gaussHardThreshold(chgg[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+2], N, DEFAULT);
				pr = ws.gaussHardThreshold(cghh[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+3], N, DEFAULT);
				pr = ws.gaussHardThreshold(cghg[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+4], N, DEFAULT);
				pr = ws.gaussHardThreshold(cggh[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+5], N, DEFAULT);
				pr = ws.gaussHardThreshold(cggg[s-1], sc, pr, Sigma, &multiSup[7*(s-1)+6], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.gaussHardThreshold(chhg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(chgh[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(chgg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cghh[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cghg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cggh[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
				pr = ws.gaussHardThreshold(cggg[s-1], sc, pr, Sigma, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == GAUSSFDR)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.gaussFDRThreshold(chhg[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(chgh[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(chgg[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cghh[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+3]);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cghg[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+4]);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cggh[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+5]);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cggg[s-1], Sigma, FDR, FDRINDEP, &multiSup[7*(s-1)+6]);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.gaussFDRThreshold(chhg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(chgh[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(chgg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cghh[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cghg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cggh[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cggg[s-1], Sigma, FDR, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			}

			mwirWT->reconstruction3D(*chhh, chhg[s-1], chgh[s-1], chgg[s-1], \
						 cghh[s-1], cghg[s-1], cggh[s-1], cggg[s-1], \
						 data, s, dec);
			extractData(data, ext[s-1]);
			*chhh = data;
		} // end for
		delete chhh; chhh= NULL; delete [] chhg; chhg = NULL; 
		delete [] chgh; chgh = NULL; delete [] chgg; chgg = NULL;
		delete [] cghh; cghh = NULL; delete [] cghg; cghg = NULL;
		delete [] cggh; cggh = NULL; delete [] cggg; cggg = NULL;		
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	} // end if
	delete mwirWT; mwirWT = NULL;
}

template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef, double lambda)
{
  int dim = data.naxis();
  int nx = data.nx(), ny = data.ny(), nz = data.nz();

  if (dim == 1)
    {
      for (int x=0; x<nx; x++)
      {
	data(x) = (coef(x) >= 0) ? origdata(x) : data(x);
	if (ITERMODE == HSD)
	  { 
	    data(x) = soft_threshold(data(x), lambda);
          }
      }
    }
  else if (dim == 2)
    {
      for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
      {
	data(x, y) = (coef(x, y) >= 0) ? origdata(x, y) : data(x, y);
	if (ITERMODE == HSD)
	  { 
	    data(x, y) = soft_threshold(data(x, y), lambda);
          }
      }
    }
  else // dim == 3
    {
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
  int dim = origdata.naxis();
  int ext[NSCALE][6];

  type_border TB = sbf[0].Border;
  int h0len = 0, h1len = 0, g0len = 0, g1len = 0, maxlen[3] = {0, 0, 0};
  for (int k=0; k<3; k++)
    {
      sbf[k].getFilterSize(h0len, h1len, g0len, g1len);
      maxlen[k] = MAX(h0len, MAX(h1len, MAX(g0len, g1len)));
    }

  MWIRWaveletTransform *mwirWT = new MWIRWaveletTransform(sbf[0], sbf[1], sbf[2]);
 
  double lambda = 0, delta = 0;
  if (ITERMODE == HSD)
    {
      if (niter <= 1) return;
      lambda = 1.;
      delta = lambda / (niter-1);
    }

  fltarray *chs = new fltarray;
  fltarray *cgs = new fltarray[NSCALE];
  fltarray *chhs = new fltarray;
  fltarray *chgs = new fltarray[NSCALE];
  fltarray *cghs = new fltarray[NSCALE];
  fltarray *cggs = new fltarray[NSCALE];
  fltarray *chhhs = new fltarray;
  fltarray *chhgs = new fltarray[NSCALE];
  fltarray *chghs = new fltarray[NSCALE];
  fltarray *chggs = new fltarray[NSCALE];
  fltarray *cghhs = new fltarray[NSCALE];
  fltarray *cghgs = new fltarray[NSCALE];
  fltarray *cgghs = new fltarray[NSCALE];
  fltarray *cgggs = new fltarray[NSCALE];
  fltarray *tempd = new fltarray;
  fltarray *ch = new fltarray;
  fltarray *cg = new fltarray;
  fltarray *chh = new fltarray;
  fltarray *chg = new fltarray;
  fltarray *cgh = new fltarray;
  fltarray *cgg = new fltarray;
  fltarray *chhh = new fltarray;
  fltarray *chhg = new fltarray;
  fltarray *chgh = new fltarray;
  fltarray *chgg = new fltarray;
  fltarray *cghh = new fltarray;
  fltarray *cghg = new fltarray;
  fltarray *cggh = new fltarray;
  fltarray *cggg = new fltarray;

  for (int i=0; i<niter; i++)
  {
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;
      *tempd = origdata;
      if (dim == 1)
	  {
	        if (VERBOSE) cerr << "Wavelet transform 1D and processing ... " << endl;
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(*tempd, maxlen, TB, s, dec, ext[s-1]);
			extData(solution, maxlen, TB, s, dec, ext[s-1]);
		        mwirWT->transform1D(*tempd, *ch, *cg, s, dec[0]);
			mwirWT->transform1D(solution, *chs, cgs[s-1], s, dec[0]);
			procArr<float>(*cg, cgs[s-1], multiSup[s-1], lambda);
			solution = *chs;
			*tempd = *ch;
		}
		*chs = *ch; // same approximation band
		if (VERBOSE) cerr << "Wavelet transform 1D complete. " << endl << "Wavelet reconstruction ... " << endl;

		for (int s=NSCALE; s>=1; s--)
		{
		    mwirWT->reconstruction1D(*chs, cgs[s-1], solution, s, dec[0]);
		    extractData(solution, ext[s-1]);
		    *chs = solution;
		}
	        if (VERBOSE) cerr << "Wavelet reconstruction complete. " << endl;
	  }
      else if (dim == 2)
	  {
	        if (VERBOSE) cerr << "Wavelet transform 2D and processing ... " << endl;
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(*tempd, maxlen, TB, s, dec, ext[s-1]);
			extData(solution, maxlen, TB, s, dec, ext[s-1]);	
		        mwirWT->transform2D(*tempd, *chh, *chg, *cgh, *cgg, s, dec);
			mwirWT->transform2D(solution, *chhs, chgs[s-1], cghs[s-1], cggs[s-1], s, dec);
			procArr<float>(*chg, chgs[s-1], multiSup[3*(s-1)], lambda);
			procArr<float>(*cgh, cghs[s-1], multiSup[3*(s-1)+1], lambda);
			procArr<float>(*cgg, cggs[s-1], multiSup[3*(s-1)+2], lambda);
			solution = *chhs;
			*tempd = *chh;
		}
		*chhs = *chh; // same approximation band
		if (VERBOSE) cerr << "Wavelet transform 2D complete. " << endl << "Wavelet reconstruction ... " << endl;

		for (int s=NSCALE; s>=1; s--)
		{
		    mwirWT->reconstruction2D(*chhs, chgs[s-1], cghs[s-1], cggs[s-1], solution, s, dec);
		    extractData(solution, ext[s-1]); 
		    *chhs = solution;
		}
		if (VERBOSE) cerr << "Wavelet reconstruction complete. " << endl;
	  }
      else // dim == 3
	  {
	        if (VERBOSE) cerr << "Wavelet transform 3D and processing ... " << endl;
		for (int s=1; s<=NSCALE; s++)
		{
		        extData(*tempd, maxlen, TB, s, dec, ext[s-1]);
			extData(solution, maxlen, TB, s, dec, ext[s-1]);	
		        mwirWT->transform3D(*tempd, *chhh, *chhg, *chgh, *chgg, *cghh, *cghg, *cggh, *cggg, s, dec);
			mwirWT->transform3D(solution, *chhhs, chhgs[s-1], chghs[s-1], chggs[s-1], cghhs[s-1], cghgs[s-1], cgghs[s-1], cgggs[s-1], s, dec);
			procArr<float>(*chhg, chhgs[s-1], multiSup[7*(s-1)], lambda);
			procArr<float>(*chgh, chghs[s-1], multiSup[7*(s-1)+1], lambda);
			procArr<float>(*chgg, chggs[s-1], multiSup[7*(s-1)+2], lambda);
			procArr<float>(*cghh, cghhs[s-1], multiSup[7*(s-1)+3], lambda);
			procArr<float>(*cghg, cghgs[s-1], multiSup[7*(s-1)+4], lambda);
			procArr<float>(*cggh, cgghs[s-1], multiSup[7*(s-1)+5], lambda);
			procArr<float>(*cggg, cgggs[s-1], multiSup[7*(s-1)+6], lambda);
			solution = *chhhs;
			*tempd = *chhh;
		}
		*chhhs = *chhh; // same approximation band
		if (VERBOSE) cerr << "Wavelet transform 3D complete. " << endl << "Wavelet reconstruction ... " << endl;
	
		for (int s=NSCALE; s>=1; s--)
		{
			mwirWT->reconstruction3D(*chhhs, chhgs[s-1], chghs[s-1], chggs[s-1], cghhs[s-1], cghgs[s-1], cgghs[s-1], cgggs[s-1], solution, s, dec);
			extractData(solution, ext[s-1]);
			*chhhs = solution;
		}
		if (VERBOSE) cerr << "Wavelet reconstruction complete. " << endl;
	  }

      posProject(solution);
      lambda -= delta;
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " complete. " << endl;
  }
  
  delete chs; delete [] cgs; chs = NULL; cgs = NULL;
  delete chhs; delete [] chgs; delete [] cghs; delete [] cggs;
  chhs = NULL; chgs = NULL; cghs = NULL; cggs = NULL;
  delete chhhs; chhhs= NULL; delete [] chhgs; chhgs = NULL; 
  delete [] chghs; chghs = NULL; delete [] chggs; chggs = NULL;
  delete [] cghhs; cghhs = NULL; delete [] cghgs; cghgs = NULL;
  delete [] cgghs; cgghs = NULL; delete [] cgggs; cgggs = NULL;
  delete mwirWT; mwirWT = NULL;
  delete tempd; tempd = NULL;
  delete ch; delete cg; ch = NULL; cg = NULL;
  delete chh; delete chg; delete cgh; delete cgg;
  chh = NULL; chg = NULL; cgh = NULL; cgg = NULL;
  delete chhh; delete chhg; delete chgh; delete chgg;
  delete cghh; delete cghg; delete cggh; delete cggg;
  chhh = NULL;  chhg = NULL;  chgh = NULL;  chgg = NULL;
  cghh = NULL;  cghg = NULL;  cggh = NULL;  cggg = NULL;
}

// wavelet direct denoising
void waveletGaussDenoise (fltarray &data)
{	
        int dim = data.naxis();
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
	    if (dim == 1)
	      multiSup = new fltarray[NSCALE];
	    else if (dim == 2)
	      multiSup = new fltarray[3*NSCALE];
	    else // dim == 3
	      multiSup = new fltarray[7*NSCALE];
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

   int dim = data->naxis();
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALE = MIN(NSCALE, scaleOfData(dim, nx, ny, nz));
   if (AUTOFDR)
     {
       if (DEC[0] && DEC[1] && DEC[2]) FDRINDEP = true;
       else FDRINDEP = false;
     }
   if (PROGMODE == GAUSSFDR) BONF = false;
   if (GENER)
     Sigma = ABS(Sigma);

   if (VERBOSE)
   {
   	cout << "Input image : " << Name_Imag_In << endl;
   	cout << "Output image : " << Name_Imag_Out << endl;
	if (!GENER)
	{
	  if (WAVE == 1)
	    cout << "Denoising wavelet : Daubechies - 4 " << endl;
	  else
	    cout << "Denoising wavelet : Haar" << endl;
	  if (DEC[0]) cout << "X - decimated ";
	  else cout << "X - undecimated ";
	  if (DEC[1]) cout << " Y - decimated ";
	  else cout << "Y - undecimated ";
	  if (DEC[2]) cout << " Z - decimated ";
	  else cout << "Z - undecimated ";
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
   	
	  cout << "Max scale(s) : " << NSCALE << endl;
	  cout << "First detection scale : " << FSCALE << endl;
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
	else
	{
	  cout << "Generating Mode" << endl;
	  cout << "Seed = " << SEED << " sigma = " << Sigma << endl;
	}
   }
   
   // Denoising
   try{
     if (VERBOSE)
     {
       if (!GENER)
	 cerr << "Denoising ..." << endl;
       else 
	 cerr << "Generating ..." << endl;
     }

     if (!GENER)
     {
       waveletGaussDenoise(*data);
       if (VERBOSE)
	 cerr << "Writing denoising result file ... " << endl;
       fits_write_fltarr(Name_Imag_Out, *data, &header);
       if (VERBOSE)
	 cerr << "Writing complete." << endl;
     }
     else
     {
       StochasticLib2 sto(SEED);
       char prename[256];
       int datalen = data->n_elem();

       for (int k=0; k<datalen; k++)
       {
	 (*data)(k) += (float) sto.Normal(0, Sigma);
       }

       strcpy(prename, Name_Imag_Out);
       if (VERBOSE)
	 cerr << "Writing generated result ... " << endl;
       fits_write_fltarr(strcat(prename, "_Gen"), *data, &header);
       if (VERBOSE)
	 cerr << "Writing complete." << endl;
     }

     if (VERBOSE)
     {
       if (!GENER)
	 cerr << "Denoising complete." << endl;
       else
	 cerr << "Generation complete." << endl;
     }

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
