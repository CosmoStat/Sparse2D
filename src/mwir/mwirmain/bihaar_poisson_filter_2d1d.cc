/**************************************************
 * PROJECT : MultiWavelenImageResto (MWIR) - version II uniquely for 2D+1D image
 * CEA
 * Filename : main2.cc
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
char   Name_Model_In[256] = ""; // input model file image
const int  METHOD_LEN = 8;          // total number of methods available
enum   MODEW { HAAR = 0, BIHAAR = 1, BSPLINE = 2, BSPLINE2 = 3, ODEGARD97 = 4 }; // wavelet modes
enum   MODE1 { KOLA = 0, MKOLA = 1, BIJA = 2, FDR = 3, ANSCOM = 4, FISZ = 5, MODEL = 6, MODELFDR = 7 };  // denoising modes
enum   MODE2 { DIRECT = 0, L1REG = 1 }; // iteration modes
MODEW  DENWAVE = BIHAAR;    // Wavelet used in the denoising
MODEW  ITRWAVE = BIHAAR;    // Wavelet used in the iterative restoration
MODE1  PROGMODE = KOLA;     // denoising mode
MODE2  ITERMODE = L1REG;    // iteration mode
int    CTRANS = 0;          // max. cyclic transitions
int    NITER = 0;           // max. number of iterations
double PROBA[METHOD_LEN] = {0., 0., 1e-3, 1e-2, 0., 0., 0., 1e-2}; // individual proba.
int    NSCALEXY = 2;        // max. number of scales along X and Y
int    NSCALEZ  = 2;        // max. number of scales along Z
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, true}; // X, Y - undecimated; Z - decimated
bool   AUTOFDR = true;      // automatically choose the FDR mode
bool   FDRINDEP = false;    // independence parameter in FDR
double BACKINTENS = -1.;    // a priori background intensity (unknown by default) 
bool   BONF = false;        // Bonferroni's multi-test schema
int    FSCALEXY = 1;          // First detection X-Y scale
int    FSCALEZ = 1;          // First detection Z scale
double GLEVEL = 0;          // Gauss-detection level (k-sigma)

type_border T_BORDER = I_MIRROR; // border type

// NOTICE : Bonferroni and FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name [input_model_name]" << endl;
  cout << "NOTICE : This program is only for pseudo-3D (2D+1D) image" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : Kolaczyk's threshold (default)" << endl;
  cout << "         M = 1 : Modified Kolaczyk's threshold" << endl;
  cout << "         M = 2 : Bijaoui-Jammal's threshold" << endl;
  cout << "         M = 3 : False Discovery Rate's threshold" << endl;
  cout << "         M = 4 : Anscombe denoising" << endl;  
  cout << "         M = 5 : Fisz denoising" << endl;
  cout << "         M = 6 : Model based p-value threshold" << endl;
  cout << "         M = 7 : Model based FDR threshold" << endl;
  cout << "    [-T value] wavelet used for modes M = 0,1,2,3" << endl;
  cout << "         T = 0 : Haar" << endl;
  cout << "         T = 1 : BiHaar 2/6 (default)" << endl;     
  cout << "         T = 2 : Haar + B-Spline along X and Y; BiHaar along Z" << endl;
  cout << "         T = 3 : Haar + B-Spline-2 along X and Y; BiHaar along Z" << endl;
  cout << "    [-S value] wavelet in the iterative restoration" << endl;
  cout << "         S = 0,1,2,3 : see -W option (default = 1)" << endl;
  cout << "         S = 4 : ODEGARD 9/7" << endl;
  cout << "    [-D] to decimate X and Y directions; if W = 2,3 this option has no effect; (not default)" << endl;
  cout << "    [-l value] prior background intensity (default = automatically estimated)" << endl;
  cout << "    [-t value] max. cycle translation for M = 5 (default = 0)" << endl;
  cout << "    [-E value] two sided Pr or FDR > 0 (default = automatically decided)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = automatically decided)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = automatically decided)" << endl;
  cout << "    [-n value] max. number of scales along X and Y (default = 2)" << endl;
  cout << "    [-N value] max. number of scales along Z (default = 2)" << endl;
  cout << "    [-F value] first detection X-Y scale (default = 1)" << endl;
  cout << "    [-f value] first detection Z scale (default = 1)" << endl;
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
    while ((c = GetOpt(argc,argv,"T:S:DM:l:t:E:s:c:n:N:F:f:I:i:B:Rv")) != -1) 
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

template <typename DATATYPE>
void setConst (to_array<DATATYPE, true> &dest, to_array<DATATYPE, true> &ref, double cst)
{
	int nx = ref.nx(), ny = ref.ny(), nz = ref.nz();
	dest.resize(nx, ny, nz);
	
	for (int x=0; x<nx; x++)
	  for (int y=0; y<ny; y++)
	    for (int z=0; z<nz; z++)
	      dest(x, y, z) = (DATATYPE)cst;		
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

// cycle translation of an image
void cycleTrans (fltarray &data, int dx, int dy, int dz)
{
        fltarray *temp = new fltarray;
	*temp = data;
	int nx = data.nx(), ny = data.ny(), nz = data.nz();
	int mx, my, mz;

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
	delete temp; temp = NULL;
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void waveletDenoise (fltarray &data, SubBandFilter sbf[], bool dec[], \
                     to_array<SUPTYPE, true> *multiSup = NULL, fltarray *model = NULL)
{
	int extxy[NSCALEXY][6];
	int extz[NSCALEZ][6];
	double fdrp;
	int N = data.n_elem();
	bool DEFAULT = (PROBA[PROGMODE] == 0);

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

	fltarray *mchh = new fltarray, *chh = new fltarray;
	fltarray *mchg = new fltarray[NSCALEXY], *chg = new fltarray[NSCALEXY];
	fltarray *mcgh = new fltarray[NSCALEXY], *cgh = new fltarray[NSCALEXY];
	fltarray *mcgg = new fltarray[NSCALEXY], *cgg = new fltarray[NSCALEXY];
	for (int sxy=1; sxy<=NSCALEXY; sxy++)
	  {
	    extData(data, maxlenxy, TB, sxy, dec, extxy[sxy-1]); 
	    mwirWT->transform3DXY(data, *chh, chg[sxy-1], cgh[sxy-1], cgg[sxy-1], sxy, dec);
	    data = *chh;
	    
	    if (model != NULL)
	    {
	    extData(*model, maxlenxy, TB, sxy, dec, extxy[sxy-1]); 
	    mwirWT->transform3DXY(*model, *mchh, mchg[sxy-1], mcgh[sxy-1], mcgg[sxy-1], sxy, dec);
	    *model = *mchh;
        }
	  }

	if (VERBOSE) cerr << "Transform 2D plane-by-plane complete." \
			  << endl << "Transform 1D along Z and thresholding ... " << endl;

	int msi = 0;
	fltarray *mch1 = new fltarray, *ch1 = new fltarray; 
    fltarray *mch2 = new fltarray, *ch2 = new fltarray;
	fltarray *mch3 = new fltarray, *ch3 = new fltarray; 
    fltarray *mch4 = new fltarray, *ch4 = new fltarray;
	fltarray *mcg1 = new fltarray[NSCALEZ], *cg1 = new fltarray[NSCALEZ];
	fltarray *mcg2 = new fltarray[NSCALEZ], *cg2 = new fltarray[NSCALEZ];
	fltarray *mcg3 = new fltarray[NSCALEZ], *cg3 = new fltarray[NSCALEZ];
	fltarray *mcg4 = new fltarray[NSCALEZ], *cg4 = new fltarray[NSCALEZ];
	fltarray *backprei = new fltarray, *backptr = NULL;
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
		
		  if (model != NULL)
		  {
		    extData(*mchh, maxlenz, TB, sz, dec, extz[sz-1]); 
		    extData(mchg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		    extData(mcgh[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		    extData(mcgg[sxy-1], maxlenz, TB, sz, dec, extz[sz-1]);
		
		    mwirWT->transform3DZ(*mchh, *mch1, mcg1[sz-1], sz, dec[2]);
		    mwirWT->transform3DZ(mchg[sxy-1], *mch2, mcg2[sz-1], sz, dec[2]);
		    mwirWT->transform3DZ(mcgh[sxy-1], *mch3, mcg3[sz-1], sz, dec[2]);
		    mwirWT->transform3DZ(mcgg[sxy-1], *mch4, mcg4[sz-1], sz, dec[2]);

		    *mchh = *mch1;
		    mchg[sxy-1] = *mch2;
		    mcgh[sxy-1] = *mch3;
		    mcgg[sxy-1] = *mch4;  
          }
	    } // end of for sz
	    
	    for (int sz=NSCALEZ; sz>=1; sz--)
	      {
		int sc[3] = {sxy, sxy, sz};
		if ((BACKINTENS >= 0) && (model == NULL)) 
		  {
		    setConst(*backprei, *ch1, BACKINTENS);
		    backptr = backprei;
		  }
		else
		  backptr = ch1;
	
		double dlen = ch1->n_elem();
		double pr = (BONF ? PROBA[PROGMODE]/dlen : PROBA[PROGMODE]);
		
        if ((sxy < FSCALEXY) && (sz < FSCALEZ))
        {
			if (sxy == NSCALEXY)
			{
               if (model == NULL)
                   setConst(cg1[sz-1], cg1[sz-1], 0);
               else
                   cg1[sz-1] = mcg1[sz-1];
                   
               if (multiSup != NULL) 
                 setConst(multiSup[msi++], cg1[sz-1], -1);  
            }
            if (model == NULL)
            {
                setConst(cg2[sz-1], cg2[sz-1], 0);
                setConst(cg3[sz-1], cg3[sz-1], 0);
                setConst(cg4[sz-1], cg4[sz-1], 0);
            }
            else
            {
                cg2[sz-1] = mcg2[sz-1];
                cg3[sz-1] = mcg3[sz-1];
                cg4[sz-1] = mcg4[sz-1];
            }
            if (multiSup != NULL)
            {
              setConst(multiSup[msi++], cg2[sz-1], -1);
              setConst(multiSup[msi++], cg3[sz-1], -1);
              setConst(multiSup[msi++], cg4[sz-1], -1);
            }
			if (sz == NSCALEZ)
			{
                if (model == NULL)
                {
                    setConst(*ch2, *ch2, 0);
                    setConst(*ch3, *ch3, 0);
                    setConst(*ch4, *ch4, 0);
                }
                else
                {
                    *ch2 = *mch2;
                    *ch3 = *mch3;
                    *ch4 = *mch4;
                }
                if (multiSup != NULL)
                {
                   setConst(multiSup[msi++], *ch2, -1);
                   setConst(multiSup[msi++], *ch3, -1);
                   setConst(multiSup[msi++], *ch4, -1);
                }
            }
        }
        else if (PROGMODE == KOLA)
		  {
		    if (multiSup != NULL)
		    {
			if (sxy == NSCALEXY)
		      pr = ws.haarKolaThreshold(*backptr, cg1[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg2[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg3[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg4[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarKolaThreshold(*backptr, *ch2, sc, pr, &multiSup[msi++], N, DEFAULT);
			    pr = ws.haarKolaThreshold(*backptr, *ch3, sc, pr, &multiSup[msi++], N, DEFAULT);
			    pr = ws.haarKolaThreshold(*backptr, *ch4, sc, pr, &multiSup[msi++], N, DEFAULT);
			  }
		    }
		    else
		    {
			if (sxy == NSCALEXY)
	          pr = ws.haarKolaThreshold(*backptr, cg1[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg2[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg3[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarKolaThreshold(*backptr, cg4[sz-1], sc, pr, NULL, N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarKolaThreshold(*backptr, *ch2, sc, pr, NULL, N, DEFAULT);
			    pr = ws.haarKolaThreshold(*backptr, *ch3, sc, pr, NULL, N, DEFAULT);
			    pr = ws.haarKolaThreshold(*backptr, *ch4, sc, pr, NULL, N, DEFAULT);
			  }
		    }
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;
		  } // if PROGMODE == KOLA
		else if (PROGMODE == MKOLA)
		  {
		    if (multiSup != NULL)
		      {
			if (sxy == NSCALEXY)
  		      pr = ws.haarMKolaThreshold(*backptr, cg1[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg2[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg3[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg4[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarMKolaThreshold(*backptr, *ch2, sc, pr, &multiSup[msi++], N, DEFAULT);
			    pr = ws.haarMKolaThreshold(*backptr, *ch3, sc, pr, &multiSup[msi++], N, DEFAULT);
			    pr = ws.haarMKolaThreshold(*backptr, *ch4, sc, pr, &multiSup[msi++], N, DEFAULT);
			  }
		      }
		    else
		      {
			if (sxy == NSCALEXY)
		      pr = ws.haarMKolaThreshold(*backptr, cg1[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg2[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg3[sz-1], sc, pr, NULL, N, DEFAULT);
			pr = ws.haarMKolaThreshold(*backptr, cg4[sz-1], sc, pr, NULL, N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarMKolaThreshold(*backptr, *ch2, sc, pr, NULL, N, DEFAULT);
			    pr = ws.haarMKolaThreshold(*backptr, *ch3, sc, pr, NULL, N, DEFAULT);
			    pr = ws.haarMKolaThreshold(*backptr, *ch4, sc, pr, NULL, N, DEFAULT);
			  }
		      }
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;
		  } // if PROGMODE == MKOLA
		else if (PROGMODE == BIJA)
		  {
		    if (multiSup != NULL)
		      {
			if (sxy == NSCALEXY)
  			  pr = ws.haarBJThreshold(*backptr, cg1[sz-1], sc, pr, &multiSup[msi++]);
			pr = ws.haarBJThreshold(*backptr, cg2[sz-1], sc, pr, &multiSup[msi++]);
			pr = ws.haarBJThreshold(*backptr, cg3[sz-1], sc, pr, &multiSup[msi++]);
			pr = ws.haarBJThreshold(*backptr, cg4[sz-1], sc, pr, &multiSup[msi++]);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarBJThreshold(*backptr, *ch2, sc, pr, &multiSup[msi++]);
			    pr = ws.haarBJThreshold(*backptr, *ch3, sc, pr, &multiSup[msi++]);
			    pr = ws.haarBJThreshold(*backptr, *ch4, sc, pr, &multiSup[msi++]);
			  }
		      }
		    else
		      {
			if (sxy == NSCALEXY)
			  pr = ws.haarBJThreshold(*backptr, cg1[sz-1], sc, pr);
			pr = ws.haarBJThreshold(*backptr, cg2[sz-1], sc, pr);
			pr = ws.haarBJThreshold(*backptr, cg3[sz-1], sc, pr);
			pr = ws.haarBJThreshold(*backptr, cg4[sz-1], sc, pr);
			if (sz == NSCALEZ)
			  {
			    pr = ws.haarBJThreshold(*backptr, *ch2, sc, pr);
			    pr = ws.haarBJThreshold(*backptr, *ch3, sc, pr);
			    pr = ws.haarBJThreshold(*backptr, *ch4, sc, pr);
			  }
		      }
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;
		  } // if PROGMODE == BIJA
		else if (PROGMODE == FDR)
		  {
		    if (multiSup != NULL)
		    {
			if (sxy == NSCALEXY)
			{
			  fdrp = ws.haarFDRThreshold(*backptr, cg1[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
			  if (VERBOSE)
			    cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			         << " cut-off p-value = " << setprecision(12) << fdrp << endl;
            }
			fdrp = ws.haarFDRThreshold(*backptr, cg2[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.haarFDRThreshold(*backptr, cg3[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.haarFDRThreshold(*backptr, cg4[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
			if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			if (sz == NSCALEZ)
			  {
			    fdrp = ws.haarFDRThreshold(*backptr, *ch2, sc, pr, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.haarFDRThreshold(*backptr, *ch3, sc, pr, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.haarFDRThreshold(*backptr, *ch4, sc, pr, FDRINDEP, &multiSup[msi++]);
			    if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			  }
		      }
		    else			  
		      {
			if (sxy == NSCALEXY)
			{
			fdrp = ws.haarFDRThreshold(*backptr, cg1[sz-1], sc, pr, FDRINDEP);
			if (VERBOSE)
			  cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
            }
			fdrp = ws.haarFDRThreshold(*backptr, cg2[sz-1], sc, pr, FDRINDEP);
			if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.haarFDRThreshold(*backptr, cg3[sz-1], sc, pr, FDRINDEP);
			if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			fdrp = ws.haarFDRThreshold(*backptr, cg4[sz-1], sc, pr, FDRINDEP);
			if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			if (sz == NSCALEZ)
			  {
			    fdrp = ws.haarFDRThreshold(*backptr, *ch2, sc, pr, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.haarFDRThreshold(*backptr, *ch3, sc, pr, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			    fdrp = ws.haarFDRThreshold(*backptr, *ch4, sc, pr, FDRINDEP);
			    if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			  }
		      } 
		  } // if PROGMODE == FDR
		else if ((PROGMODE == FISZ) || (PROGMODE == ANSCOM))
		  {
    		int sc[3] = {NSCALEXY, NSCALEXY, NSCALEZ};
		    // use the asymptotic estimation : sigma = 1
		    if (multiSup != NULL)
		    {
			if (sxy == NSCALEXY)
			  pr = ws.gaussHardThreshold(cg1[sz-1], sc, pr, 1., &multiSup[msi++], N, DEFAULT);
			pr = ws.gaussHardThreshold(cg2[sz-1], sc, pr, 1., &multiSup[msi++], N, DEFAULT);
			pr = ws.gaussHardThreshold(cg3[sz-1], sc, pr, 1., &multiSup[msi++], N, DEFAULT);	       	
			pr = ws.gaussHardThreshold(cg4[sz-1], sc, pr, 1., &multiSup[msi++], N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.gaussHardThreshold(*ch2, sc, pr, 1., &multiSup[msi++], N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch3, sc, pr, 1., &multiSup[msi++], N, DEFAULT);	     	
			    pr = ws.gaussHardThreshold(*ch4, sc, pr, 1., &multiSup[msi++], N, DEFAULT);	   	
			  }
		    }
		    else
		    {
			if (sxy == NSCALEXY)
			  pr = ws.gaussHardThreshold(cg1[sz-1], sc, pr, 1., NULL, N, DEFAULT);
			pr = ws.gaussHardThreshold(cg2[sz-1], sc, pr, 1., NULL, N, DEFAULT);
			pr = ws.gaussHardThreshold(cg3[sz-1], sc, pr, 1., NULL, N, DEFAULT);	       	
			pr = ws.gaussHardThreshold(cg4[sz-1], sc, pr, 1., NULL, N, DEFAULT);
			if (sz == NSCALEZ)
			  {
			    pr = ws.gaussHardThreshold(*ch2, sc, pr, 1., NULL, N, DEFAULT);
			    pr = ws.gaussHardThreshold(*ch3, sc, pr, 1., NULL, N, DEFAULT);	     	
			    pr = ws.gaussHardThreshold(*ch4, sc, pr, 1., NULL, N, DEFAULT);	   	
			  }
		    }
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;
		  } // if PROGMODE == FISZ or ANSCOM
        else if (PROGMODE == MODEL)
        {
 			  if (multiSup != NULL)
			  {
       			if (sxy == NSCALEXY)
				  pr = ws.haarModelThreshold (*mch1, mcg1[sz-1], cg1[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg2[sz-1], cg2[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg3[sz-1], cg3[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg4[sz-1], cg4[sz-1], sc, pr, &multiSup[msi++], N, DEFAULT);
			    if (sz == NSCALEZ)
			    {				
				pr = ws.haarModelThreshold (*mch1, *mch2, *ch2, sc, pr, &multiSup[msi++], N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, *mch3, *ch3, sc, pr, &multiSup[msi++], N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, *mch4, *ch4, sc, pr, &multiSup[msi++], N, DEFAULT);
                }
			  }
			  else
			  {
       			if (sxy == NSCALEXY)
				  pr = ws.haarModelThreshold (*mch1, mcg1[sz-1], cg1[sz-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg2[sz-1], cg2[sz-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg3[sz-1], cg3[sz-1], sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, mcg4[sz-1], cg4[sz-1], sc, pr, NULL, N, DEFAULT);
			    if (sz == NSCALEZ)
			    {				
				pr = ws.haarModelThreshold (*mch1, *mch2, *ch2, sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, *mch3, *ch3, sc, pr, NULL, N, DEFAULT);
				pr = ws.haarModelThreshold (*mch1, *mch4, *ch4, sc, pr, NULL, N, DEFAULT);
                }
			  }
		    if (VERBOSE)
		      cerr << "scalexy = " << sxy << " scalez = " << sz \
			   << " cut-off p-value = " << setprecision(12) << pr << endl;
        }
        else if (PROGMODE == MODELFDR)
        {
			  if (multiSup != NULL)
			  {
       			if (sxy == NSCALEXY)
       			{
				  fdrp = ws.haarModelFDRThreshold(*mch1, mcg1[sz-1], cg1[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
				  if (VERBOSE)
			  cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
                }
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg2[sz-1], cg2[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
				if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg3[sz-1], cg3[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
				if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg4[sz-1], cg4[sz-1], sc, pr, FDRINDEP, &multiSup[msi++]);
            	if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			
			    if (sz == NSCALEZ)
			    {				
				  fdrp = ws.haarModelFDRThreshold(*mch1, *mch2, *ch2, sc, pr, FDRINDEP, &multiSup[msi++]);
				  if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				  fdrp = ws.haarModelFDRThreshold(*mch1, *mch3, *ch3, sc, pr, FDRINDEP, &multiSup[msi++]);
			  	  if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			  	  fdrp = ws.haarModelFDRThreshold(*mch1, *mch4, *ch4, sc, pr, FDRINDEP, &multiSup[msi++]);
				  if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
                }
			  }
			  else			  
			  {
       			if (sxy == NSCALEXY)
       			{
				  fdrp = ws.haarModelFDRThreshold(*mch1, mcg1[sz-1], cg1[sz-1], sc, pr, FDRINDEP);
				  if (VERBOSE)
			  cerr << "g1 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
                }
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg2[sz-1], cg2[sz-1], sc, pr, FDRINDEP);
				if (VERBOSE)
			  cerr << "g2 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg3[sz-1], cg3[sz-1], sc, pr, FDRINDEP);
				if (VERBOSE)
			  cerr << "g3 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				fdrp = ws.haarModelFDRThreshold(*mch1, mcg4[sz-1], cg4[sz-1], sc, pr, FDRINDEP);
            	if (VERBOSE)
			  cerr << "g4 : scalexy = " << sxy << " scalez = " << sz \
			       << " cut-off p-value = " << setprecision(12) << fdrp << endl;
			
			    if (sz == NSCALEZ)
			    {				
				  fdrp = ws.haarModelFDRThreshold(*mch1, *mch2, *ch2, sc, pr, FDRINDEP);
			      if (VERBOSE)
			      cerr << "h2 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				  fdrp = ws.haarModelFDRThreshold(*mch1, *mch3, *ch3, sc, pr, FDRINDEP);
				  if (VERBOSE)
			      cerr << "h3 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
				  fdrp = ws.haarModelFDRThreshold(*mch1, *mch4, *ch4, sc, pr, FDRINDEP);
				  if (VERBOSE)
			      cerr << "h4 : scalexy = " << sxy << " scalez = " << sz \
				   << " cut-off p-value = " << setprecision(12) << fdrp << endl;
               }
			 }                 
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

        if (model != NULL)
        {
		*mchh = *mch1;
		mchg[sxy-1] = *mch2;
		mcgh[sxy-1] = *mch3;
		mcgg[sxy-1] = *mch4;
		mwirWT->reconstruction3DZ(*mchh, mcg1[sz-1], *mch1, sz, dec[2]);
		mwirWT->reconstruction3DZ(mchg[sxy-1], mcg2[sz-1], *mch2, sz, dec[2]);
		mwirWT->reconstruction3DZ(mcgh[sxy-1], mcg3[sz-1], *mch3, sz, dec[2]);
		mwirWT->reconstruction3DZ(mcgg[sxy-1], mcg4[sz-1], *mch4, sz, dec[2]);
		extractData(*mch1, extz[sz-1]);
		extractData(*mch2, extz[sz-1]);
		extractData(*mch3, extz[sz-1]);
		extractData(*mch4, extz[sz-1]);
        }
	      } // for sz
	    mwirWT->reconstruction3DXY(*ch1, *ch2, *ch3, *ch4, *chh, sxy, dec);
	    extractData(*chh, extxy[sxy-1]);
        
        if (model != NULL)
        {
	    mwirWT->reconstruction3DXY(*mch1, *mch2, *mch3, *mch4, *mchh, sxy, dec);
	    extractData(*mchh, extxy[sxy-1]);
        }
	  } // for sxy
	data = *chh;
    if (model != NULL)
       *model = *mchh;
	if (VERBOSE) cerr << "Thresholding complete." << endl;
	
	delete chh; chh = NULL; delete [] chg; chg = NULL; 
	delete [] cgh; cgh = NULL; delete [] cgg; cgg = NULL;
	delete ch1; delete ch2; delete ch3; delete ch4;
	ch1 = ch2 = ch3 = ch4 = NULL;
	delete [] cg1; cg1 = NULL; delete [] cg2; cg2 = NULL;
	delete [] cg3; cg3 = NULL; delete [] cg4; cg4 = NULL;

	delete mchh; mchh = NULL; delete [] mchg; mchg = NULL; 
	delete [] mcgh; mcgh = NULL; delete [] mcgg; mcgg = NULL;
	delete mch1; delete mch2; delete mch3; delete mch4;
	mch1 = mch2 = mch3 = mch4 = NULL;
	delete [] mcg1; mcg1 = NULL; delete [] mcg2; mcg2 = NULL;
	delete [] mcg3; mcg3 = NULL; delete [] mcg4; mcg4 = NULL;

	delete backprei; backprei = NULL; backptr = NULL;
	delete mwirWT; mwirWT = NULL;
}

// for different modes of iteration in BHAAR case
template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef, double lambda)
{
  int nx = data.nx(), ny = data.ny(), nz = data.nz();

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
void multiSupIter (fltarray &origdata, fltarray &solution, SubBandFilter sbf[], bool dec[], fltarray *multiSup, int niter, fltarray *model = NULL)
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
	  fisztr.fisz3D(data);
	  if (VERBOSE)
	    cerr << "Fisz transform complete. " << endl;
	  
	  // denoising and restoration
	  waveletDenoise<float>(data, sbf, DEC);
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
  delete origData; delete sumData;
  origData = NULL; sumData = NULL;
}

// Haar direct denoising
void HaarDenoise (fltarray &data, fltarray *model = NULL)
{	
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
      multiSup = new fltarray[(NSCALEZ-1)*4+7 + ((NSCALEZ-1)*3+6)*(NSCALEXY-1)];
    }
  if (VERBOSE)
    cerr << "Initial denoising ... " << endl;
  if (DENWAVE == HAAR)
    waveletDenoise<float>(data, sbf0, DEC, multiSup, model);
  else if (DENWAVE == BIHAAR)
    waveletDenoise<float>(data, sbf, DEC, multiSup, model);
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
          if (model->naxis() != 3)
          {
           cerr << "Error : Model is not 3D" << endl;
	       delete data; data = NULL;
           delete model; model = NULL;
           exit(-1);
          }      
       }
   }
   else if (Name_Model_In[0] != '\0')
   {
        cerr << "Model file = " << Name_Model_In << endl;
        cerr << "Error : only M = 6 or 7 can be used with model." << endl;
        exit(-1);
   }

   if (data->naxis() != 3)
     {
       cerr << "Image is not 3D" << endl;
       delete data; data = NULL;
       exit(-1);
     }
   
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALEXY = MIN(NSCALEXY, scaleOfData(2, nx, ny, nz));
   NSCALEZ = MIN(NSCALEZ, scaleOfData(1, nz, nz, nz));
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

   if (VERBOSE)
   {
   	cout << "Input image : " << Name_Imag_In << endl;
   	cout << "Output image : " << Name_Imag_Out << endl;
   	cout << "Input model image : " << Name_Model_In << endl;

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
   	
	cout << "Max scale(s) X and Y: " << NSCALEXY << endl;
	cout << "Max scale(s) Z: " << NSCALEZ << endl;
    cout << "First detection X-Y scale : " << FSCALEXY << endl;
    cout << "First detection Z scale : " << FSCALEZ << endl;

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
   
   // Denoising
   try{
     if (VERBOSE)
       {
	 cerr << "Denoising ..." << endl;
       }
     
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
     
     if (VERBOSE)
     {
       cerr << "Denoising complete." << endl;
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
