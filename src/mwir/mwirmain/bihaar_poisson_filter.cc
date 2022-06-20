/**************************************************
 * PROJECT : MultiWavelenImageResto (MWIR) - version I
 * CEA
 * Filename : main1.cc
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
// #include "utils/mersenne.cpp"          
#include "stocc.h"               
// #include "utils/stoc1.cpp"             
// #include "utils/stoc2.cpp"  
// #include "utils/userintf.cpp"           
using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];       // input image file
char   Name_Imag_Out[256];  // output image file
char   Name_Model_In[256] = ""; // input model file image
const int METHOD_LEN = 8;       // total number of methods available
enum   MODEW { HAAR = 0, BIHAAR = 1, BSPLINE = 2, BSPLINE2 = 3, ODEGARD97 = 4 }; // wavelet modes
enum   MODE1 { KOLA = 0, MKOLA = 1, BIJA = 2, FDR = 3, ANSCOM = 4, FISZ = 5, MODEL = 6, MODELFDR = 7 };  // denoising modes
enum   MODE2 { DIRECT = 0, L1REG = 1 }; // iteration modes
MODEW  DENWAVE = BIHAAR;    // Wavelet used in the denoising
MODEW  ITRWAVE = BIHAAR;    // Wavelet used in the iterative restoration
MODE1  PROGMODE = KOLA;     // denoising mode
MODE2  ITERMODE = L1REG;    // iteration mode
int    CTRANS = 0;          // max. cyclic transitions
int    NITER = 0;           // max. number of iterations
double PROBA[METHOD_LEN] = {0., 0., 1e-3, 1e-2, 0., 0., 0., 1e-2}; // proba or FDR.
int    NSCALE = 3;          // max. number of scales
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, true}; // X, Y - undecimated; Z - decimated
bool   AUTOFDR = true;      // automatically choose the FDR dependence parameter
bool   FDRINDEP = false;    // dependence parameter in FDR
double BACKINTENS = -1.;    // a priori background intensity (unknown by default) 
bool   BONF = false;        // Bonferroni's multi-test schema
int    FSCALE = 1;          // First detection scale
double GLEVEL = 0;          // Gauss-detection level (k-sigma)
bool FDRTHRESH = false;
bool   GENER = false;       // generating an image with Poisson noise
double SCALF  = 1.;         // scaling factor 
long int SEED = time(NULL); // random seed

type_border T_BORDER = I_MIRROR; // border type

// NOTICE : Bonferroni and FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name [input_model_name]" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    =========== GENERATION OPTIONS ==========" << endl;
  cout << "    [-G] generate an Poisson noisy image from the input; background can be added by -l (not default)" << endl;
  cout << "    [-e value] generating integer seed (default = current time)" << endl;
  cout << "    [-a value] scaling factor of the input intensity image while generating (default = 1)" << endl;
  cout << endl;
  cout << "    =========== DENOISING OPTIONS ==========" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : Kolaczyk's threshold (default)" << endl;
  cout << "         M = 1 : Modified Kolaczyk's threshold" << endl;
  cout << "         M = 2 : Bijaoui-Jammal's threshold" << endl;
  cout << "         M = 3 : False Discovery Rate's threshold" << endl;
  cout << "         M = 4 : Anscombe denoising" << endl;  
  cout << "         M = 5 : Fisz denoising" << endl;
  cout << "         M = 6 : Model based p-value threshold" << endl;
  cout << "         M = 7 : Model based FDR threshold" << endl;
  cout << "    [-m] Use FDR in Anscombe or Fisz" <<endl;
  cout << "    [-T value] wavelet used for modes M = 0,1,2,3" << endl;
  cout << "         T = 0 : Haar" << endl;
  cout << "         T = 1 : BiHaar 2/6 (default)" << endl;     
  cout << "         T = 2 : Haar + B-Spline along X and Y; BiHaar along Z" << endl;
  cout << "         T = 3 : Haar + B-Spline-2 along X and Y; BiHaar along Z" << endl;
  cout << "    [-S value] wavelet in the iterative restoration" << endl;
  cout << "         S = 0,1,2,3 : see -T option (default = 1)" << endl;
  cout << "         S = 4 : ODEGARD 9/7" << endl;
  cout << "    [-D] to decimate X and Y directions; if W = 2,3 this option has no effect; (not default)" << endl;
  cout << "    [-l value] prior background intensity (default = automatically estimated)" << endl;
  cout << "    [-t value] max. cycle translation for M = 5 (default = 0)" << endl;
  cout << "    [-E value] two sided Pr or FDR > 0 (default = automatically decided)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = automatically decided)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = automatically decided)" << endl;
  cout << "    [-n value] number of scales (default = 3)" << endl;
  cout << "    [-F value] first detection scale (default = 1)" << endl;
  cout << "    [-I value] iteration modes for M = 0,1,2,3,6,7" << endl;
  cout << "                      I = 0 : Direct iterative mode" << endl;
  cout << "                      I = 1 : L1-regularized iterative mode (default)" << endl;
  cout << "    [-i value] number of iteration (default = 0)" << endl; 
  cout << "    [-B value] border type" << endl;
  cout << "                      B = 0 : symmetric border (default)" << endl;
  cout << "                      B = 1 : periodic border" << endl;
  cout << "    [-R] Bonferroni correction (not default and no effect in FDR)" << endl;
  cout << "    [-v] verbose mode" << endl;
  cout << endl;
}
 
//*********************************************************************

// GET COMMAND LINE ARGUMENTS
static void filtinit(int argc, char *argv[])
{
    int c;  

    // get options 
    while ((c = GetOpt(argc,argv,"GmDvRS:T:M:e:a:l:t:E:s:c:n:F:I:i:B:")) != -1) 
    {
      switch (c) 
	{	
	case 'S':
	  int iwave;
	  if (sscanf(OptArg, "%d", &iwave) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (iwave == 0)
	    ITRWAVE = HAAR;
	  else if (iwave == 1)
	    ITRWAVE = BIHAAR;
	  else if (iwave == 2)
	    ITRWAVE = BSPLINE;
	  else if (iwave == 3)
	    ITRWAVE = BSPLINE2;
	  else
	    ITRWAVE = ODEGARD97;
	  break;  

	case 'T':
	  int dwave;
	  if (sscanf(OptArg, "%d", &dwave) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (dwave == 0)
	    DENWAVE = HAAR;
	  else if (dwave == 1)
	    DENWAVE = BIHAAR;
	  else if (dwave == 2)
	    DENWAVE = BSPLINE;
	  else if (dwave == 3)
	    DENWAVE = BSPLINE2;
	  else
	    DENWAVE = ODEGARD97;
	  break;

	case 'G':
	  GENER = true;
	  break;

	case 'm':
	  FDRTHRESH = true;
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
	    PROGMODE = KOLA;
	  else if (pmd == 1) 
	    PROGMODE = MKOLA;
	  else if (pmd == 2)
	    PROGMODE = BIJA;
	  else if (pmd == 3)
	    PROGMODE = FDR;
	  else if (pmd == 4)
	    PROGMODE = ANSCOM;
	  else if (pmd == 5)
	    PROGMODE = FISZ;
      else if (pmd == 6)
        PROGMODE = MODEL;
      else
        PROGMODE = MODELFDR;
	  break;
	
	case 'e':
	  if (sscanf(OptArg, "%ld", &SEED) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  

	case 'a':
	  if (sscanf(OptArg, "%lf", &SCALF) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  

	case 'l':
	  if (sscanf(OptArg, "%lf", &BACKINTENS) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  
	
	case 't':
	  if (sscanf(OptArg, "%d", &CTRANS) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  break;  
	
	case 'E':
	  double proba;
	  if (sscanf(OptArg, "%lf", &proba) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  for (int i=0; i<METHOD_LEN; i++) PROBA[i] = proba;
	  break;
	  
	case 's':
	  if (sscanf(OptArg, "%lf", &GLEVEL) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
      for (int i=0; i<METHOD_LEN; i++) PROBA[i] = 2. * (1 - Utils<double>::cumNormal(GLEVEL));
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
	  else 
	    ITERMODE = L1REG;
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
    else { cerr << "Error : Input image required." << endl; usage(argv); exit(-1); }

    if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
    else { cerr << "Error : Output image required." << endl; usage(argv); exit(-1); }
    
    if (OptInd < argc) strcpy(Name_Model_In, argv[OptInd++]); 

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

  if ((BORDERTYPE != I_PERIOD) && ( (NITER == 0) || ((NITER != 0) && (ITRWAVE == DENWAVE)) ))
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

// cycle translation of an image
void cycleTrans (fltarray &data, int dx, int dy, int dz)
{
        fltarray *temp = new fltarray;
	*temp = data;
	int dim = data.naxis();
	int nx = data.nx(), ny = data.ny(), nz = data.nz();
	int mx, my, mz;

	if (dim == 1)
	{
		if (dx == 0) return;
		for (int x=0; x<nx; x++)
		{
		    mx = (x-dx) % nx;
		    if (mx < 0) mx += nx;
		    data(x) = (*temp)(mx);
	        }
	}
	else if (dim == 2)
	{
		if ((dx == 0) && (dy == 0)) return;

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
		    mx = (x-dx) % nx; 
		    my = (y-dy) % ny;
		    if (mx < 0) mx += nx;
		    if (my < 0) my += ny;
		    data(x, y) = (*temp)(mx, my);
		}
	}
	else // dim == 3
	{
		if ((dx == 0) && (dy == 0) &&(dz == 0)) return;

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
		    mx = (x-dx) % nx; 
		    my = (y-dy) % ny; 
		    mz = (z-dz) % nz;
		    if (mx < 0) mx += nx;
		    if (my < 0) my += ny;
		    if (mz < 0) mz += nz;
		    data(x, y, z) = (*temp)(mx, my, mz);
		}
	}
	delete temp; temp = NULL;
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void waveletDenoise (fltarray &data, SubBandFilter sbf[], bool dec[], \
		             to_array<SUPTYPE, true> *multiSup = NULL, fltarray *model = NULL)
{
	int dim = data.naxis();
	int ext[NSCALE][6];
	double fdrp;
	int N = data.n_elem();
	bool DEFAULT = (PROBA[PROGMODE] == 0);

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
		fltarray *mch = new fltarray, *ch = new fltarray;
		fltarray *mcg = new fltarray[NSCALE], *cg = new fltarray[NSCALE];

		for (int s=1; s<=NSCALE; s++)
		{
		    extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform1D(data, *ch, cg[s-1], s, dec[0]);
			data = *ch;
			if (model != NULL)
			{
		        extData(*model, maxlen, TB, s, dec, ext[s-1]); 
			    mwirWT->transform1D(*model, *mch, mcg[s-1], s, dec[0]);
			    *model = *mch;                
            }
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;
		
		fltarray *backprei = new fltarray, *backptr = NULL;
		for (int s=NSCALE; s>=1; s--)
		{
		    int sc[1] = {s};
			if ((BACKINTENS >= 0) && (model == NULL))
			{
				setConst(*backprei, *ch, BACKINTENS);
				backptr = backprei;
			}
			else
				backptr = ch;
			
			double dlen = cg[s-1].n_elem();
			double pr = (BONF ? PROBA[PROGMODE]/dlen : PROBA[PROGMODE]);
			
            if (s < FSCALE)
            {
                if (model == NULL)
                    setConst(cg[s-1], cg[s-1], 0);
                else
                    cg[s-1] = mcg[s-1];
                    
                if (multiSup != NULL)
                    setConst(multiSup[s-1], cg[s-1], -1);  
            }
			else if (PROGMODE == KOLA)
			{
			  if (multiSup != NULL)
				pr = ws.haarKolaThreshold(*backptr, cg[s-1], sc, pr, &multiSup[s-1], N, DEFAULT);
			  else 	pr = ws.haarKolaThreshold(*backptr, cg[s-1], sc, pr, NULL, N, DEFAULT);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == MKOLA)
			{
			  if (multiSup != NULL)
				pr = ws.haarMKolaThreshold(*backptr, cg[s-1], sc, pr, &multiSup[s-1], N, DEFAULT);
			  else 	pr = ws.haarMKolaThreshold(*backptr, cg[s-1], sc, pr, NULL, N, DEFAULT);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == BIJA)
			{
			  if (multiSup != NULL)
				pr = ws.haarBJThreshold(*backptr, cg[s-1], sc, pr, &multiSup[s-1]);
			  else 	pr = ws.haarBJThreshold(*backptr, cg[s-1], sc, pr);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
		        }
			else if (PROGMODE == FDR)
			{
			  if (multiSup != NULL)
			      fdrp = ws.haarFDRThreshold(*backptr, cg[s-1], sc, pr, FDRINDEP, &multiSup[s-1]);
			  else 	
			      fdrp = ws.haarFDRThreshold(*backptr, cg[s-1], sc, pr, FDRINDEP);
			  if (VERBOSE)
				cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			}
			else if ((PROGMODE == FISZ) || (PROGMODE == ANSCOM))
			{
              int sc[1] = {NSCALE};
			  // use the asymptotic estimation : sigma = 1
			  if (multiSup != NULL)
			  {
			    if (FDRTHRESH)
			    	fdrp = ws.gaussFDRThreshold (cg[s-1], 1., pr, FDRINDEP, &multiSup[s-1]);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], sc, pr, 1., &multiSup[s-1], N, DEFAULT);
			  }
			  else 
			  {
			    if (FDRTHRESH)
		  	        fdrp = ws.gaussFDRThreshold (cg[s-1], 1., pr, FDRINDEP);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], sc, pr, 1., NULL, N, DEFAULT);
			  }	  		  

			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << (FDRTHRESH ? fdrp : pr) << endl;
			}
			else if (PROGMODE == MODEL)
			{
			  if (multiSup != NULL)
			    pr = ws.haarModelThreshold(*mch, mcg[s-1], cg[s-1], sc, pr, &multiSup[s-1], N, DEFAULT);
			  else 
                pr = ws.haarModelThreshold(*mch, mcg[s-1], cg[s-1], sc, pr, NULL, N, DEFAULT);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;                 
            }
			else if (PROGMODE == MODELFDR)
			{
			  if (multiSup != NULL)
			    fdrp = ws.haarModelFDRThreshold(*mch, mcg[s-1], cg[s-1], sc, pr, FDRINDEP, &multiSup[s-1]);
			  else 
                fdrp = ws.haarModelFDRThreshold(*mch, mcg[s-1], cg[s-1], sc, pr, FDRINDEP);
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;                 
            }
			mwirWT->reconstruction1D(*ch, cg[s-1], data, s, dec[0]);
			extractData(data, ext[s-1]);
			*ch = data;
			if (model != NULL)
			{
   			    mwirWT->reconstruction1D(*mch, mcg[s-1], *model, s, dec[0]);
			    extractData(*model, ext[s-1]);
			    *mch = *model;
            }
		}
		delete ch; delete [] cg; delete backprei;
		ch = NULL; cg = NULL; backprei = NULL; backptr = NULL;
        delete mch; delete [] mcg; 
        mch = NULL; mcg = NULL;
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	}
	else if (dim == 2)
	{
	    if (VERBOSE) cerr << "Wavelet transform 2D ... " << endl;
		fltarray *mchh = new fltarray, *chh = new fltarray;
		fltarray *mchg = new fltarray[NSCALE], *chg = new fltarray[NSCALE];
		fltarray *mcgh = new fltarray[NSCALE], *cgh = new fltarray[NSCALE];
		fltarray *mcgg = new fltarray[NSCALE], *cgg = new fltarray[NSCALE];
		for (int s=1; s<=NSCALE; s++)
		{
		    extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform2D(data, *chh, chg[s-1], cgh[s-1], cgg[s-1], s, dec);
			data = *chh;
			if (model != NULL)
			{
		        extData(*model, maxlen, TB, s, dec, ext[s-1]); 
			    mwirWT->transform2D(*model, *mchh, mchg[s-1], mcgh[s-1], mcgg[s-1], s, dec);
			    *model = *mchh;                
            }
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;

		fltarray *backprei = new fltarray, *backptr = NULL;
		for (int s=NSCALE; s>=1; s--)
		{
		    int sc[2] = {s, s};
			if ((BACKINTENS >= 0) && (model == NULL))
			{
				setConst(*backprei, *chh, BACKINTENS);
				backptr = backprei;
			}
			else
				backptr = chh;

			double dlen = chg[s-1].n_elem();
			double pr = (BONF ? PROBA[PROGMODE]/dlen : PROBA[PROGMODE]);

            if (s < FSCALE)
            {
                if (model == NULL)
                {
                    setConst(chg[s-1], chg[s-1], 0);
                    setConst(cgh[s-1], cgh[s-1], 0);
                    setConst(cgg[s-1], cgg[s-1], 0);
                }
                else
                {
                    chg[s-1] = mchg[s-1];
                    cgh[s-1] = mcgh[s-1];
                    cgg[s-1] = mcgg[s-1];
                }
                if (multiSup != NULL)
                {
                 setConst(multiSup[3*(s-1)], chg[s-1], -1);  
                 setConst(multiSup[3*(s-1)+1], cgh[s-1], -1);  
                 setConst(multiSup[3*(s-1)+2], cgg[s-1], -1);  
                }
            }
			else if (PROGMODE == KOLA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarKolaThreshold(*backptr, chg[s-1], sc, pr, &multiSup[3*(s-1)], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cgh[s-1], sc, pr, &multiSup[3*(s-1)+1], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cgg[s-1], sc, pr, &multiSup[3*(s-1)+2], N, DEFAULT);
			    }
			  else
			    {
			      pr = ws.haarKolaThreshold(*backptr, chg[s-1], sc, pr, NULL, N, DEFAULT);
			      pr = ws.haarKolaThreshold(*backptr, cgh[s-1], sc, pr, NULL, N, DEFAULT);
			      pr = ws.haarKolaThreshold(*backptr, cgg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == MKOLA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarMKolaThreshold(*backptr, chg[s-1], sc, pr, &multiSup[3*(s-1)], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cgh[s-1], sc, pr, &multiSup[3*(s-1)+1], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cgg[s-1], sc, pr, &multiSup[3*(s-1)+2], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.haarMKolaThreshold(*backptr, chg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cgh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cgg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == BIJA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarBJThreshold(*backptr, chg[s-1], sc, pr, &multiSup[3*(s-1)]);
				pr = ws.haarBJThreshold(*backptr, cgh[s-1], sc, pr, &multiSup[3*(s-1)+1]);
				pr = ws.haarBJThreshold(*backptr, cgg[s-1], sc, pr, &multiSup[3*(s-1)+2]);
			    }
			  else
			    {
				pr = ws.haarBJThreshold(*backptr, chg[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cgh[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cgg[s-1], sc, pr);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == FDR)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.haarFDRThreshold(*backptr, chg[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cgh[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cgg[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.haarFDRThreshold(*backptr, chg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cgh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cgg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			}
			else if ((PROGMODE == FISZ) || (PROGMODE == ANSCOM))
			{
   		      int sc[2] = {NSCALE, NSCALE};
			  // use the asymptotic estimation : sigma = 1
			  if (multiSup != NULL)
			    {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chg[s-1], 1., pr, FDRINDEP, &multiSup[3*(s-1)]);
	  			    if (VERBOSE)
			   		cerr << "scale(hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgh[s-1], 1., pr, FDRINDEP, &multiSup[3*(s-1)+1]);
	  			    if (VERBOSE)
			   		cerr << "scale(gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgg[s-1], 1., pr, FDRINDEP, &multiSup[3*(s-1)+2]);
	  			    if (VERBOSE)
			   		cerr << "scale(gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	}
			        else
			        {	
				    pr = ws.gaussHardThreshold(chg[s-1], sc, pr, 1., &multiSup[3*(s-1)], N, DEFAULT);
				    pr = ws.gaussHardThreshold(cgh[s-1], sc, pr, 1., &multiSup[3*(s-1)+1], N, DEFAULT);
				    pr = ws.gaussHardThreshold(cgg[s-1], sc, pr, 1., &multiSup[3*(s-1)+2], N, DEFAULT);				
		  		    if (VERBOSE)
				    	cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
				}
			    }
			  else
			    {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chg[s-1], 1., pr, FDRINDEP);
	  			    if (VERBOSE)
			   		cerr << "scale(hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgh[s-1], 1., pr, FDRINDEP);
	  			    if (VERBOSE)
			   		cerr << "scale(gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgg[s-1], 1., pr, FDRINDEP);
	  			    if (VERBOSE)
			   		cerr << "scale(gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	}
			        else
			        {	
				    pr = ws.gaussHardThreshold(chg[s-1], sc, pr, 1., NULL, N, DEFAULT);
				    pr = ws.gaussHardThreshold(cgh[s-1], sc, pr, 1., NULL, N, DEFAULT);
				    pr = ws.gaussHardThreshold(cgg[s-1], sc, pr, 1., NULL, N, DEFAULT);				
				    if (VERBOSE)
				    	cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
				}
			    }
			}
			else if (PROGMODE == MODEL)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarModelThreshold (*mchh, mchg[s-1], chg[s-1], sc, pr, &multiSup[3*(s-1)], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchh, mcgh[s-1], cgh[s-1], sc, pr, &multiSup[3*(s-1)+1], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchh, mcgg[s-1], cgg[s-1], sc, pr, &multiSup[3*(s-1)+2], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.haarModelThreshold (*mchh, mchg[s-1], chg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchh, mcgh[s-1], cgh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchh, mcgg[s-1], cgg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
            }
			else if (PROGMODE == MODELFDR)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.haarModelFDRThreshold(*mchh, mchg[s-1], chg[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchh, mcgh[s-1], cgh[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchh, mcgg[s-1], cgg[s-1], sc, pr, FDRINDEP, &multiSup[3*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.haarModelFDRThreshold(*mchh, mchg[s-1], chg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchh, mcgh[s-1], cgh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchh, mcgg[s-1], cgg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }                 
            }
			mwirWT->reconstruction2D(*chh, chg[s-1], cgh[s-1], cgg[s-1], data, s, dec);
			extractData(data, ext[s-1]);
			*chh = data;
			if (model != NULL)
			{
			mwirWT->reconstruction2D(*mchh, mchg[s-1], mcgh[s-1], mcgg[s-1], *model, s, dec);
			extractData(*model, ext[s-1]);
			*mchh = *model;                      
            }
		}
		delete chh; delete [] chg; delete [] cgh; delete [] cgg; delete backprei;
		chh = NULL; chg = NULL; cgh = NULL; cgg = NULL; backprei = NULL; backptr = NULL;
		delete mchh; delete [] mchg; delete [] mcgh; delete [] mcgg; 
		mchh = NULL; mchg = NULL; mcgh = NULL; mcgg = NULL;
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	}
	else // dim == 3
	{
	    if (VERBOSE) cerr << "Wavelet transform 3D ... " << endl;
		fltarray *mchhh = new fltarray, *chhh = new fltarray;
		fltarray *mchhg = new fltarray[NSCALE], *chhg = new fltarray[NSCALE];
		fltarray *mchgh = new fltarray[NSCALE], *chgh = new fltarray[NSCALE];
		fltarray *mchgg = new fltarray[NSCALE], *chgg = new fltarray[NSCALE];
		fltarray *mcghh = new fltarray[NSCALE], *cghh = new fltarray[NSCALE];
		fltarray *mcghg = new fltarray[NSCALE], *cghg = new fltarray[NSCALE];
		fltarray *mcggh = new fltarray[NSCALE], *cggh = new fltarray[NSCALE];
		fltarray *mcggg = new fltarray[NSCALE], *cggg = new fltarray[NSCALE];
		for (int s=1; s<=NSCALE; s++)
		{
		    extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform3D(data, *chhh, chhg[s-1], chgh[s-1], chgg[s-1], \
				               cghh[s-1], cghg[s-1], cggh[s-1], cggg[s-1], s, dec);
			data = *chhh;
			
			if (model != NULL)
			{
		    extData(*model, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform3D(*model, *mchhh, mchhg[s-1], mchgh[s-1], mchgg[s-1], \
				               mcghh[s-1], mcghg[s-1], mcggh[s-1], mcggg[s-1], s, dec);
			*model = *mchhh;
            }
		}
		if (VERBOSE) cerr << "Transform complete." << endl << "Wavelet thresholding ... " << endl;

		fltarray *backprei = new fltarray, *backptr = NULL;
		for (int s=NSCALE; s>=1; s--)
		{
		        int sc[3] = {s, s, s};
			if ((BACKINTENS >= 0) && (model == NULL))
			{
				setConst(*backprei, *chhh, BACKINTENS);
				backptr = backprei;
			}
			else
				backptr = chhh;

			double dlen = chhg[s-1].n_elem();
			double pr = (BONF ? PROBA[PROGMODE]/dlen : PROBA[PROGMODE]);

            if (s < FSCALE)
            {
                if (model == NULL)
                {
                    setConst(chhg[s-1], chhg[s-1], 0);
                    setConst(chgh[s-1], chgh[s-1], 0);
                    setConst(chgg[s-1], chgg[s-1], 0);
                    setConst(cghh[s-1], cghh[s-1], 0);
                    setConst(cghg[s-1], cghg[s-1], 0);
                    setConst(cggh[s-1], cggh[s-1], 0);
                    setConst(cggg[s-1], cggg[s-1], 0);
                }
                else
                {
                    chhg[s-1] = mchhg[s-1];
                    chgh[s-1] = mchgh[s-1];
                    chgg[s-1] = mchgg[s-1];
                    cghh[s-1] = mcghh[s-1];
                    cghg[s-1] = mcghg[s-1];
                    cggh[s-1] = mcggh[s-1];
                    cggg[s-1] = mcggg[s-1];
                }
                if (multiSup != NULL)
                {                
                   setConst(multiSup[7*(s-1)], chhg[s-1], -1);  
                   setConst(multiSup[7*(s-1)+1], chgh[s-1], -1);  
                   setConst(multiSup[7*(s-1)+2], chgg[s-1], -1);  
                   setConst(multiSup[7*(s-1)+3], cghh[s-1], -1);  
                   setConst(multiSup[7*(s-1)+4], cghg[s-1], -1);  
                   setConst(multiSup[7*(s-1)+5], cggh[s-1], -1);  
                   setConst(multiSup[7*(s-1)+6], cggg[s-1], -1);  
                }
            }
			else if (PROGMODE == KOLA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarKolaThreshold(*backptr, chhg[s-1], sc, pr, &multiSup[7*(s-1)], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, chgh[s-1], sc, pr, &multiSup[7*(s-1)+1], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, chgg[s-1], sc, pr, &multiSup[7*(s-1)+2], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cghh[s-1], sc, pr, &multiSup[7*(s-1)+3], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cghg[s-1], sc, pr, &multiSup[7*(s-1)+4], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cggh[s-1], sc, pr, &multiSup[7*(s-1)+5], N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cggg[s-1], sc, pr, &multiSup[7*(s-1)+6], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.haarKolaThreshold(*backptr, chhg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, chgh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, chgg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cghh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cghg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cggh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarKolaThreshold(*backptr, cggg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == MKOLA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarMKolaThreshold(*backptr, chhg[s-1], sc, pr, &multiSup[7*(s-1)], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, chgh[s-1], sc, pr, &multiSup[7*(s-1)+1], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, chgg[s-1], sc, pr, &multiSup[7*(s-1)+2], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cghh[s-1], sc, pr, &multiSup[7*(s-1)+3], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cghg[s-1], sc, pr, &multiSup[7*(s-1)+4], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cggh[s-1], sc, pr, &multiSup[7*(s-1)+5], N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cggg[s-1], sc, pr, &multiSup[7*(s-1)+6], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.haarMKolaThreshold(*backptr, chhg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, chgh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, chgg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cghh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cghg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cggh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarMKolaThreshold(*backptr, cggg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == BIJA)
			{
			  if (multiSup != NULL)
			    {
				pr = ws.haarBJThreshold(*backptr, chhg[s-1], sc, pr, &multiSup[7*(s-1)]);
				pr = ws.haarBJThreshold(*backptr, chgh[s-1], sc, pr, &multiSup[7*(s-1)+1]);
				pr = ws.haarBJThreshold(*backptr, chgg[s-1], sc, pr, &multiSup[7*(s-1)+2]);
				pr = ws.haarBJThreshold(*backptr, cghh[s-1], sc, pr, &multiSup[7*(s-1)+3]);
				pr = ws.haarBJThreshold(*backptr, cghg[s-1], sc, pr, &multiSup[7*(s-1)+4]);
				pr = ws.haarBJThreshold(*backptr, cggh[s-1], sc, pr, &multiSup[7*(s-1)+5]);
				pr = ws.haarBJThreshold(*backptr, cggg[s-1], sc, pr, &multiSup[7*(s-1)+6]);
			    }
			  else
			    {
				pr = ws.haarBJThreshold(*backptr, chhg[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, chgh[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, chgg[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cghh[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cghg[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cggh[s-1], sc, pr);
				pr = ws.haarBJThreshold(*backptr, cggg[s-1], sc, pr);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			}
			else if (PROGMODE == FDR)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.haarFDRThreshold(*backptr, chhg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, chgh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, chgg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cghh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+3]);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cghg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+4]);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cggh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+5]);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cggg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+6]);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.haarFDRThreshold(*backptr, chhg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, chgh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, chgg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cghh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cghg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cggh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarFDRThreshold(*backptr, cggg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			}
			else if ((PROGMODE == FISZ) || (PROGMODE == ANSCOM))
			{
   		     int sc[3] = {NSCALE, NSCALE, NSCALE};
			  // use the asymptotic estimation : sigma = 1
			  if (multiSup != NULL)
			    {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chhg[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)]);
				    if (VERBOSE)
				        cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (chgh[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+1]);
				    if (VERBOSE)
				        cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (chgg[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+2]);
				    if (VERBOSE)
				        cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cghh[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+3]);
				    if (VERBOSE)
				        cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cghg[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+4]);
				    if (VERBOSE)
				        cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cggh[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+5]);
				    if (VERBOSE)
				        cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cggg[s-1], 1., pr, FDRINDEP, &multiSup[7*(s-1)+6]);
				    if (VERBOSE)
				        cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			        }
			        else
			        {
					pr = ws.gaussHardThreshold(chhg[s-1], sc, pr, 1., &multiSup[7*(s-1)], N, DEFAULT);
					pr = ws.gaussHardThreshold(chgh[s-1], sc, pr, 1., &multiSup[7*(s-1)+1], N, DEFAULT);
					pr = ws.gaussHardThreshold(chgg[s-1], sc, pr, 1., &multiSup[7*(s-1)+2], N, DEFAULT);	       	
					pr = ws.gaussHardThreshold(cghh[s-1], sc, pr, 1., &multiSup[7*(s-1)+3], N, DEFAULT);
					pr = ws.gaussHardThreshold(cghg[s-1], sc, pr, 1., &multiSup[7*(s-1)+4], N, DEFAULT);
					pr = ws.gaussHardThreshold(cggh[s-1], sc, pr, 1., &multiSup[7*(s-1)+5], N, DEFAULT);	     	
					pr = ws.gaussHardThreshold(cggg[s-1], sc, pr, 1., &multiSup[7*(s-1)+6], N, DEFAULT);	   	
				        if (VERBOSE)
					    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
			        }
			    }
			  else
			    {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chhg[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (chgh[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (chgg[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cghh[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cghg[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cggh[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cggg[s-1], 1., pr, FDRINDEP);
				    if (VERBOSE)
				        cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			        }
			        else
			        {
					pr = ws.gaussHardThreshold(chhg[s-1], sc, pr, 1., NULL, N, DEFAULT);
					pr = ws.gaussHardThreshold(chgh[s-1], sc, pr, 1., NULL, N, DEFAULT);
					pr = ws.gaussHardThreshold(chgg[s-1], sc, pr, 1., NULL, N, DEFAULT);				
					pr = ws.gaussHardThreshold(cghh[s-1], sc, pr, 1., NULL, N, DEFAULT);
					pr = ws.gaussHardThreshold(cghg[s-1], sc, pr, 1., NULL, N, DEFAULT);
					pr = ws.gaussHardThreshold(cggh[s-1], sc, pr, 1., NULL, N, DEFAULT);				
					pr = ws.gaussHardThreshold(cggg[s-1], sc, pr, 1., NULL, N, DEFAULT);				
				        if (VERBOSE)
					    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
				}
			    }
			}
			else if (PROGMODE == MODEL)
			{
 			  if (multiSup != NULL)
			    {
				pr = ws.haarModelThreshold (*mchhh, mchhg[s-1], chhg[s-1], sc, pr, &multiSup[7*(s-1)], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mchgh[s-1], chgh[s-1], sc, pr, &multiSup[7*(s-1)+1], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mchgg[s-1], chgg[s-1], sc, pr, &multiSup[7*(s-1)+2], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcghh[s-1], cghh[s-1], sc, pr, &multiSup[7*(s-1)+3], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcghg[s-1], cghg[s-1], sc, pr, &multiSup[7*(s-1)+4], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcggh[s-1], cggh[s-1], sc, pr, &multiSup[7*(s-1)+5], N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcggg[s-1], cggg[s-1], sc, pr, &multiSup[7*(s-1)+6], N, DEFAULT);
			    }
			  else
			    {
				pr = ws.haarModelThreshold (*mchhh, mchhg[s-1], chhg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mchgh[s-1], chgh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mchgg[s-1], chgg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcghh[s-1], cghh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcghg[s-1], cghg[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcggh[s-1], cggh[s-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mchhh, mcggg[s-1], cggg[s-1], sc, pr, NULL, N, DEFAULT);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(12) << pr << endl;
            }
            else if (PROGMODE == MODELFDR)
            {
			  if (multiSup != NULL)
			    {
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchhg[s-1], chhg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchgh[s-1], chgh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchgg[s-1], chgg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcghh[s-1], cghh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+3]);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcghg[s-1], cghg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+4]);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcggh[s-1], cggh[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+5]);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcggg[s-1], cggg[s-1], sc, pr, FDRINDEP, &multiSup[7*(s-1)+6]);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }
			  else			  
			    {
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchhg[s-1], chhg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hhg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchgh[s-1], chgh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mchgg[s-1], chgg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (hgg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcghh[s-1], cghh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcghg[s-1], cghg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ghg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcggh[s-1], cggh[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggh) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mchhh, mcggg[s-1], cggg[s-1], sc, pr, FDRINDEP);
				if (VERBOSE)
				  cerr << "scale (ggg) = " << s << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    }                 
            }
			mwirWT->reconstruction3D(*chhh, chhg[s-1], chgh[s-1], chgg[s-1], \
						 cghh[s-1], cghg[s-1], cggh[s-1], cggg[s-1], \
						 data, s, dec);
			extractData(data, ext[s-1]);
			*chhh = data;
			if (model != NULL)
			{
			mwirWT->reconstruction3D(*mchhh, mchhg[s-1], mchgh[s-1], mchgg[s-1], \
						 mcghh[s-1], mcghg[s-1], mcggh[s-1], mcggg[s-1], \
						 *model, s, dec);
			extractData(*model, ext[s-1]);
			*mchhh = *model;
            }
		}
		delete chhh; chhh= NULL; delete [] chhg; chhg = NULL; 
		delete [] chgh; chgh = NULL; delete [] chgg; chgg = NULL;
		delete [] cghh; cghh = NULL; delete [] cghg; cghg = NULL;
		delete [] cggh; cggh = NULL; delete [] cggg; cggg = NULL;
		delete backprei; backprei = NULL; backptr = NULL;
		delete mchhh; mchhh= NULL; delete [] mchhg; mchhg = NULL; 
		delete [] mchgh; mchgh = NULL; delete [] mchgg; mchgg = NULL;
		delete [] mcghh; mcghh = NULL; delete [] mcghg; mcghg = NULL;
		delete [] mcggh; mcggh = NULL; delete [] mcggg; mcggg = NULL;
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	} // end if
	delete mwirWT; mwirWT = NULL;
}

// for different modes of iteration
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
	if (ITERMODE == L1REG)
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
	if (ITERMODE == L1REG)
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
	    if (ITERMODE == L1REG)
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

// Anscombe denoising
void anscombeDenoise (fltarray &data)
{
  // wavelet filter configuration
  SubBandFilter sbf[3] = { SubBandFilter(F_DAUBE_8, NORM_L2), \
			   SubBandFilter(F_DAUBE_8, NORM_L2), \
			   SubBandFilter(F_DAUBE_8, NORM_L2) };
  sbf[0].Border = sbf[1].Border = sbf[2].Border = T_BORDER;
  
  // Anscombe transform
  if (VERBOSE)
    cerr << "Anscombe transform ... " << endl;
  Utils<float>::anscombeTransform(data);	
  if (VERBOSE)
    cerr << "Anscombe transform complete. " << endl;
	
  // denoising and restoration
  waveletDenoise<float>(data, sbf, DEC);
  Utils<float>::invAnscombeTransform(data);
  posProject(data);
}

// Fisz denosing
void fiszDenoise (fltarray &data, int DX, int DY, int DZ)
{
	int dim = data.naxis();
	int dlen = data.n_elem();
	
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;

	// a sum of data for calculate the mean
	fltarray *sumData = new fltarray;

	int c = 0;
	int dx = (int)floor(DX/2.), dy = (int)floor(DY/2.), dz = (int)floor(DZ/2.);
	for (int tz=-dz; tz<=dz; tz++)
	for (int ty=-dy; ty<=dy; ty++)
	for (int tx=-dx; tx<=dx; tx++)
	{
	  if (VERBOSE)
	    cerr << "Cycle transition : tx = " << tx << " ty = " << ty << " tz = " << tz << endl;
		cycleTrans (data, tx, ty, tz);
	
		// data extension to do the fisz transform
		FiszTransform<float> fisztr;	
		int ext[6];
		fisztr.dataExtension(data, ext);
    
		// wavelet filter configuration
		SubBandFilter sbf[3] = { SubBandFilter(F_DAUBE_8, NORM_L2), \
					 SubBandFilter(F_DAUBE_8, NORM_L2), \
					 SubBandFilter(F_DAUBE_8, NORM_L2) };
   		sbf[0].Border = sbf[1].Border = sbf[2].Border = T_BORDER;

		if (VERBOSE)
		  cerr << "Fisz transform ... " << endl;
		// Fisz transform
		if (dim == 1)
			fisztr.fisz1D(data);
		else if (dim == 2)
			fisztr.fisz2D(data);
		else // dim == 3
			fisztr.fisz3D(data);
		if (VERBOSE)
		  cerr << "Fisz transform complete. " << endl;
	
		// denoising and restoration
		waveletDenoise<float>(data, sbf, DEC);
		if (dim == 1)
			fisztr.ifisz1D(data);
		else if (dim == 2)
			fisztr.ifisz2D(data);
		else // dim == 3
			fisztr.ifisz3D(data);
		fisztr.dataExtraction(data, ext);
		posProject(data);
		
		cycleTrans (data, -tx, -ty, -tz);
		
		if (c == 0)
		  *sumData = data;
		else
		  *sumData += data;
		
		c++;
		data = *origData;
	}
	
	if (VERBOSE)
	  cerr << " Calculate mean of all the cycle spins ... " << endl;
	if (c != 0)
	{
		for (int j=0; j<dlen; j++)
           data(j) = (*sumData)(j) / (float)c;
	}
	if (VERBOSE)
	  cerr << " Calculation complete. " << endl;

	// free all the variables	
	delete origData;  delete sumData;
	origData = NULL;  sumData = NULL;
}

// solve with the iterative algo. in using the multiresolution support
void multiSupIter (fltarray &origdata, fltarray &solution, SubBandFilter sbf[], bool dec[], fltarray *multiSup, int niter, fltarray *model = NULL)
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
  if (ITERMODE == L1REG)
    {
      if (niter <= 1) return;
      lambda = 1.;
      delta = lambda / (niter-1);
    }
    
  if (model != NULL)
  {
      solution = *model; // use the prior model as initial solution
      solution -= *model;
      origdata -= *model;
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
    
      if (model != NULL)
      {
          solution += *model;
          posProject(solution);
          solution -= *model;
      }
      else
         posProject(solution);
      lambda -= delta;
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " complete. " << endl;
  }
  
  if (model != NULL)
  {
     solution += *model;
     origdata += *model;
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

// Haar denoising
void HaarDenoise (fltarray &data, fltarray *model = NULL)
{	
    int dim = data.naxis();
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;
	// wavelet filter configuration
	SubBandFilter sbb = SubBandFilter(F_BI2HAAR, NORM_L1);
	SubBandFilter sbh = SubBandFilter(F_HAAR, NORM_L1);
	SubBandFilter sbu = SubBandFilter(U_HAAR_B3S, NORM_L1);
	SubBandFilter sbu2 = SubBandFilter(U_HAAR_B3S2, NORM_L1);
	sbb.Border = sbh.Border = sbu.Border = sbu2.Border = T_BORDER;
	
	SubBandFilter sbf[3]  = { sbb, sbb, sbb };
	SubBandFilter sbf0[3] = { sbh, sbh, sbh };
	SubBandFilter sbf2[3] = { sbu, sbu, sbb };
	SubBandFilter sbf3[3] = { sbu2, sbu2, sbb };

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
	if (DENWAVE == HAAR)
	  waveletDenoise<float>(data, sbf0, DEC, multiSup, model);
	else if (DENWAVE == BIHAAR)
	  waveletDenoise<float>(data, sbf,  DEC, multiSup, model);
	else if (DENWAVE == BSPLINE)
	  waveletDenoise<float>(data, sbf2, DEC, multiSup, model);
	else
	  waveletDenoise<float>(data, sbf3, DEC, multiSup, model);
	posProject(data);
	if (VERBOSE)
	  cerr << "Entering into the iterative denoising ..." << endl;
	if (NITER > 1)
	  {
        SubBandFilter sbo = SubBandFilter(F_ODEGARD_9_7);
        sbo.Border = T_BORDER;
        SubBandFilter mssbf[3] = { sbo, sbo, sbo };
        
	    if (ITRWAVE == HAAR)
	      multiSupIter (*origData, data, sbf0, DEC, multiSup, NITER-1, model);
	    else if (ITRWAVE == BIHAAR)
	      multiSupIter (*origData, data, sbf, DEC, multiSup, NITER-1, model);
	    else if (ITRWAVE == BSPLINE)
	      multiSupIter (*origData, data, sbf2, DEC, multiSup, NITER-1, model);
	    else if (ITRWAVE == BSPLINE2)
	      multiSupIter (*origData, data, sbf3, DEC, multiSup, NITER-1, model);
	    else
	      multiSupIter (*origData, data, mssbf, DEC, multiSup, NITER-1, model);
	  }
	if (VERBOSE)
	  cerr << "Iteration complete." << endl;
    if (model != NULL)
     {
        data -= *model;
     }	
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

   fitsstruct header, mheader;
   fltarray *data = new fltarray;
   fltarray *model = NULL;
       
   // Get command line arguments, open input file(s) if necessary
   filtinit(argc, argv);
   fits_read_fltarr(Name_Imag_In, *data, &header);

   header.bitpix = BP_FLOAT;
   header.origin = cmd;

   if ((PROGMODE == MODEL) || (PROGMODE == MODELFDR))
   {
      if (Name_Model_In == "")
      {
          cerr << "Error : Model modes require a model file." << endl;
          delete data; data = NULL;
          exit (-1);
      }
      else
      {
          model = new fltarray;
          fits_read_fltarr(Name_Model_In, *model, &mheader);
      }
   }
   else if (Name_Model_In[0] != '\0')
   {
        cerr << "Model file = " << Name_Model_In << endl;
        cerr << "Error : only M = 6 or 7 can be used with model." << endl;
        exit(-1);
   }
     
   int dim = data->naxis();
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALE = MIN(NSCALE, scaleOfData(dim, nx, ny, nz));
   if (((DENWAVE == BSPLINE) || (DENWAVE == BSPLINE2)) && (PROGMODE != FISZ) && (PROGMODE != ANSCOM)) DEC[0] = DEC[1] = false;
   if (AUTOFDR)
     {
       if (DEC[0] && DEC[1] && DEC[2] && \
           (DENWAVE == HAAR) && \
           ( (BACKINTENS >= 0) || (PROGMODE == MODEL) || (PROGMODE == MODELFDR) ) 
          ) 
          FDRINDEP = true;
       else FDRINDEP = false;
     }
   if ((PROGMODE == FDR) || (PROGMODE == MODELFDR)) BONF = false;
   SCALF = MAX(SCALF, 0);

   if (VERBOSE)
   {
   	cout << "Input image : " << Name_Imag_In << endl;
   	cout << "Output image : " << Name_Imag_Out << endl;
   	cout << "Input model image : " << Name_Model_In << endl;

	if (!GENER)
	{
	  if ((PROGMODE != FISZ) && (PROGMODE != ANSCOM))
	    {
	      if (DENWAVE == HAAR)
		cout << " Denoising wavelet : Haar " << endl;
	      else if (DENWAVE == BSPLINE)
		cout << " Denoising wavelet : X and Y is Haar + B-Spline; Z is BiHaar" << endl;
	      else if (DENWAVE == BSPLINE2)
		cout << " Denoising wavelet : X and Y is Haar + B-Spline (2); Z is BiHaar" << endl;
	      else
		cout << " Denoising wavelet : BiHaar" << endl;
	    }
	  if (DEC[0]) cout << " X - decimated ";
	  else cout << " X - undecimated ";
	  if (DEC[1]) cout << " Y - decimated ";
	  else cout << " Y - undecimated ";
	  if (DEC[2]) cout << " Z - decimated ";
	  else cout << " Z - undecimated ";
	  cout << endl;

	  if ((PROGMODE != FISZ) && (PROGMODE != ANSCOM))
	  {
          cout << "Mode : ";
	      if (PROGMODE == KOLA)
		cout << "(Bi)orthogonal" << " Haar De-noising - Kolaczyk's threshold" << endl;
	      else if (PROGMODE == MKOLA)
		cout << "(Bi)orthogonal" << " Haar De-noising - Modified Kolaczyk's threshold" << endl;
	      else if (PROGMODE == BIJA)
		cout << "(Bi)orthogonal" << " Haar De-noising - Bijaoui-Jammal's threshold" << endl;
	      else if (PROGMODE == FDR)
		cout << "(Bi)orthogonal" << " Haar De-noising - FDR's threshold" << endl;
	      else if (PROGMODE == MODEL)
		cout << "(Bi)orthogonal" << " Model based Haar De-noising" << endl;
	      else if (PROGMODE == MODELFDR)
		cout << "(Bi)orthogonal" << " Model based FDR Haar De-noising" << endl;
        
        if ((PROGMODE == FDR) || (PROGMODE == MODELFDR))
        {
			cout << "FDR mode : ";
		    if (FDRINDEP) cout << "independence mode";
      		else cout << "dependence mode";
	    	cout << endl;	    
        }

	    if ((BACKINTENS >= 0) && (PROGMODE != MODEL) && (PROGMODE != MODELFDR))
	  	cout << "Background prior intensity : " << BACKINTENS << endl;
	    else if ((PROGMODE == MODEL) || (PROGMODE == MODELFDR))
        cout << "Background prior given by model " << endl;
        else
	  	cout << "Background prior intensity is estimated automatically " << endl;
	      
	      if (ITERMODE == DIRECT)
		cout << "Iterative Mode : Direct" << endl;
	      else
		cout << "Iterative Mode : L1 - Regularized" << endl;
	      cout << "Max. Iteration : " << NITER << endl;
	  }
	  else if (PROGMODE == FISZ)
	  {
	    cout << "Mode : " << "Fisz --> Daubechies (order 4) De-noising" << endl;
	    cout << "Cycle translation : " << CTRANS << endl;
      }
	  else if (PROGMODE == ANSCOM)
	    cout << "Mode : " << "Anscombe --> Daubechies (order 4) De-noising" << endl;
	    
   	
	  cout << "Max scale(s) : " << NSCALE << endl;
	  cout << "First detection scale : " << FSCALE << endl;
	  
	  if ((PROGMODE != FDR) && (PROGMODE != MODELFDR))
	  {
	    if (PROBA[PROGMODE] != 0)
	      cout << "Proba. : " << PROBA[PROGMODE] << endl;
	    else
	      cout << "Proba. : auto-deciding" << endl;
	  }
	  else
	    cout << "FDR : " << PROBA[PROGMODE] << endl;

	  cout << "Bonferroni's schema : " << (BONF ? "yes" : "no") << endl;
	  cout << "Iterative reconstruction wavelet : ";
	  if (ITRWAVE == 0)
	    cout << "Haar" << endl;
	  else if (ITRWAVE == 1)
	    cout << "BiHaar" << endl;
	  else if (ITRWAVE == 1)
	    cout << "Bspline" << endl;
	  else if (ITRWAVE == 1)
	    cout << "Bspline (2)" << endl;
      else	  
	    cout << "Odegard 9/7" << endl;

	  cout << "Border : " << (T_BORDER == I_MIRROR ? "symmetric" : "periodic") << endl;
	}
	else
	{
	  cout << "Generating Mode" << endl;
	  cout << "Seed = " << SEED << " scaling = " << SCALF << endl;
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
       if (PROGMODE == FISZ)
	 fiszDenoise(*data, MAX(0, MIN(CTRANS, nx-1)), MAX(0, MIN(CTRANS, ny-1)), MAX(0, MIN(CTRANS, nz-1)));
       else if (PROGMODE == ANSCOM)
	 anscombeDenoise(*data);
       else
	 HaarDenoise(*data, model);

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
       double dval;
       for (int k=0; k<datalen; k++)
       {
	 dval = (*data)(k);
	 (*data)(k) = dval * SCALF + MAX(0, BACKINTENS);
       }
       
       strcpy(prename, Name_Imag_Out);
       if (VERBOSE)
	 cerr << "Writing intensity image ... " << endl;
       fits_write_fltarr(strcat(prename, "_Src"), *data, &header);
       if (VERBOSE)
	 cerr << "Writing complete." << endl;

       if (VERBOSE)
	 cerr << "Generating image with noise ... " << endl;
       for (int k=0; k<datalen; k++)
       {
	 dval = (*data)(k);
	 (*data)(k) = (float) sto.Poisson( (double) dval );
       }
       if (VERBOSE)
	 cerr << "Generation complete." << endl;

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
     if ((PROGMODE == MODEL) || (PROGMODE == MODELFDR))
     {
        delete model; model = NULL;
     }
   }
   catch (MWIRException mwirExcept)
   {
   	cerr << "Exception : " << endl;
   	cerr << mwirExcept.getReason() << endl;
   	exit (-1);
   }

    return 0;
}
