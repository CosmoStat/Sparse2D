/**************************************************
 * PROJECT : Ansc and Fisz
 * CEA
 * Filename : AtrousPoisDenoise.cc
 * This is the main file of PoisMSVST
 **************************************************/
 
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "cdflib.h"
#include "GlobalInc.h"
#include "IM_IO.h"
#include "IM_Obj.h"
#include "PoisMSVSTException.h"
#include "Atrous.h"
#include "B3VSTAtrous.h"
#include "Wavelet.h"
#include "ImLib.h"
#include "Border.h"
#include "Fisz.h"

using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];   // input image file
char   Name_Imag_Out[256];  // output image file
const int METHOD_LEN = 2;   // 2 thresholding methods available (p-value and FDR)
enum   MODE1 { ANSC = 0, FISZ = 1 };  // denoising modes
enum   MODE2 { DIRECT = 0, L1REG = 1 }; // iteration modes
MODE1  PROGMODE = ANSC;     // denoising mode
MODE2  ITERMODE = L1REG;    // iteration mode
int    CTRANS = 0;          // max. cyclic transitions
int    NITER = 0;          // max. number of iterations
double PROBA[METHOD_LEN] = {0.000465, 0.1}; // default cut p-value at 3.5*sigma and FDR = 0.1
int    NSCALE = 5;          // max. number of scales
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, false}; // undecimated
bool   FDRINDEP = false;    // dependence parameter in FDR
bool   KILLLAST = false;    // ignore the last approximation scale
bool   DETPOS   = false;    // detect only positive coefficients
int    FSCALE = 1;          // First detection scale
double GLEVEL = 3.5;          // Gauss-detection level (k-sigma)
bool   TRANSF2 = false;     // ver2 in the transform
bool   FDRTHRESH = false;   // FDR threshold 
bool   POSPROJ = true;      // positivity projection
bool   CBIAS   = true;     // correct the Anscombe's bais
// NOTICE : FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Anscombe and Fisz estimations using B3-Spline" << endl;
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-T] use g=Id-h*h as iteration band-pass filter" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : Anscombe (default)" << endl;
  cout << "         M = 1 : Fisz" << endl;
  cout << "    [-b] do NOT correct the Anscombe's bias" <<endl;  
  cout << "    [-m] Use FDR in Anscombe or Fisz" <<endl;
  cout << "    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)" << endl;
  cout << "    [-n value] number of scales (default = 5)" << endl;
  cout << "    [-F value] first detection scale (default = 1)" << endl;
  cout << "    [-K] ignore the last approximation band (used with iteration)" << endl;
  cout << "    [-p] detect only the positive coefficients (used with iteration)" << endl;
  cout << "    [-N] do NOT do positivity projection" << endl;
  cout << "    [-t value] max. cycle translation for Fisz (default = 0)" << endl;
  cout << "    [-I value] iteration modes for M = 0,1" << endl;
  cout << "                      I = 0 : Direct iterative mode" << endl;
  cout << "                      I = 1 : L1-regularized iterative mode (default)" << endl;
  cout << "    [-i value] number of iteration (default = 0)" << endl; 
  cout << "    [-v] verbose mode" << endl;
  cout << endl;
}
 
//*********************************************************************

// GET COMMAND LINE ARGUMENTS
static void filtinit(int argc, char *argv[])
{
    int c;  

    // get options 
    while ((c = GetOpt(argc,argv,"M:E:t:s:c:n:F:I:i:TbKpNmv")) != -1) 
    {
       switch (c) 
	   {	
    case 'M':
	  int pmd;
	  if (sscanf(OptArg, "%d", &pmd) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  if (pmd == 0)
	    PROGMODE = ANSC;
	  else 
        PROGMODE = FISZ;
	  break;
	
	case 'm':
	  FDRTHRESH = true;
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
	  
	case 't':
	  if (sscanf(OptArg, "%d", &CTRANS) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
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
	
	case 'K': 
	  KILLLAST = true;
	  break;

	case 'b': 
	  CBIAS = false;
	  break;

	case 'N': 
	  POSPROJ = false;
	  break;

	case 'T': 
	  TRANSF2 = true;
	  break;

	case 'p': 
	  DETPOS = true;
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
void display (convert_to_array<DATATYPE, true> &data)
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

template <typename DATATYPE>
void setConst (convert_to_array<DATATYPE, true> &dest, convert_to_array<DATATYPE, true> &ref, double cst)
{
	int nx = ref.nx(), ny = ref.ny(), nz = ref.nz();
	int len = ref.n_elem();

	dest.resize(nx, ny, nz);
    for (int i=0; i<len; i++)
		dest(i) = (DATATYPE)cst;		
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void b3SplineDenoise (fltarray &data, convert_to_array<SUPTYPE, true> *multiSup)
{
	int dim = data.naxis();
	double pr, fdrp; // threshold p-value and FDR-threshold p-value
	int N = data.n_elem();
    bool ms = (multiSup != NULL);
    
	// wavelet filter configuration
	SplineAtrous atrous(3);
	WaveletShrinkage<float, SUPTYPE> ws;

	fltarray *ch = new fltarray;
	fltarray *cg = new fltarray[NSCALE];
    
	for (int s=1; s<=NSCALE; s++)
	{
		atrous.transform(data, *ch, cg[s-1], s);
		data = *ch;
		    
        if (s < FSCALE)
        {
            setConst(cg[s-1], data, 0);
            if (ms)
            {
                multiSup[s-1].resize(cg[s-1].nx(), cg[s-1].ny(), cg[s-1].nz());
                for (int i=0; i<N; i++) (multiSup[s-1])(i) = -1;
            }
        }		
	    else if ((PROGMODE == FISZ) || (PROGMODE == ANSC))
		{
			  // use the asymptotic estimation : sigma = 1 at the initial scale = 0
			  double sig = sqrt(Utils<double>::b3CoefVar (dim, s));
			  if (VERBOSE)
			     cerr << "sig(" << s << ") = " << sig << endl;
			     
			  if (multiSup != NULL)
			  {
			    if (FDRTHRESH)
			    	fdrp = ws.gaussFDRThreshold (cg[s-1], sig, PROBA[1], FDRINDEP, &multiSup[s-1]);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], PROBA[0], sig, &multiSup[s-1]);
			  }
			  else 
			  {
			    if (FDRTHRESH)
		  	        fdrp = ws.gaussFDRThreshold (cg[s-1], sig, PROBA[1], FDRINDEP);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], PROBA[0], sig, NULL);
			  }	  		  

			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << (FDRTHRESH ? fdrp : pr) << endl;
		}
	}
	for (int s=NSCALE; s>=1; s--)
	{
        atrous.recons(*ch, cg[s-1], data, s);
        *ch = data;
    }
    
	delete ch; delete [] cg; ch = NULL; cg = NULL; 
}

// for different modes of iteration
template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, convert_to_array<SUPTYPE, true> &coef, double lambda)
{
  int n = data.n_elem();

  for (int x=0; x<n; x++)
  {
      if (DETPOS)
          data(x) = ((coef(x) >= 0) && (origdata(x) >= 0)) ? origdata(x) : data(x);
      else
          data(x) = (coef(x) >= 0) ? origdata(x) : data(x);
      if (ITERMODE == L1REG)
           data(x) = soft_threshold(data(x), lambda);
  }
}

// put to zeros the negative values
void posProject (fltarray &data)
{	
    if (POSPROJ)
    {
    	float value;
	    int len = data.n_elem();
	
    	for (int i=0; i<len; i++)
	    {
	     	value = data(i);
	   	    data(i) = MAX(value, 0);
	    }
    }
}

// solve with the iterative algo. in using the multiresolution support
void multiSupIter (fltarray &origdata, fltarray &solution, fltarray *multiSup, int niter)
{
  int dim = origdata.naxis();
  SplineAtrous atrous = SplineAtrous(3, 1., I_MIRROR, TRANSF2);
  double lambda = 0, delta = 0;
  
  if (ITERMODE == L1REG)
  {
      if (niter <= 1) return;
      lambda = 1.;
      delta = lambda / (niter-1);
  }
    
  fltarray *chs = new fltarray;
  fltarray *cgs = new fltarray[NSCALE];
  fltarray *tempd = new fltarray;
  fltarray *ch = new fltarray;
  fltarray *cg = new fltarray;

  for (int i=0; i<niter; i++)
  {
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;
      *tempd = origdata;
      
	  for (int s=1; s<=NSCALE; s++)
	  {
	    atrous.transform(*tempd, *ch, *cg, s);
		atrous.transform(solution, *chs, cgs[s-1], s);
		procArr<float>(*cg, cgs[s-1], multiSup[s-1], lambda);
		solution = *chs;
		*tempd = *ch;
	  }
	  if ((KILLLAST) && (i==niter-1))
        setConst(*chs, *ch, 0.);    
	  else
        *chs = *ch; // same approximation band

	  for (int s=NSCALE; s>=1; s--)
	  {
	    atrous.recons(*chs, cgs[s-1], solution, s);
	    *chs = solution;
      }
      
      posProject(solution);
      lambda -= delta;
  }
  
  delete chs; delete [] cgs; chs = NULL; cgs = NULL;
  delete tempd; tempd = NULL;
  delete ch; delete cg; ch = NULL; cg = NULL;
}

// Anscombe denoising
template <typename SUPTYPE>
void anscombeDenoise (fltarray &data, convert_to_array<SUPTYPE, true> *multiSup=NULL)
{
  // Anscombe transform
  if (VERBOSE)
    cerr << "Anscombe transform ... " << endl;
  Utils<float>::anscombeTransform(data);	
  if (VERBOSE)
    cerr << "Anscombe transform complete. " << endl;
	
  // denoising and restoration
  b3SplineDenoise<float>(data, multiSup);
  Utils<float>::invAnscombeTransform(data, CBIAS);
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
		b3SplineDenoise<float>(data, NULL);
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
           data(j) = ((*sumData)(j)) / ((float)c);
	}
	if (VERBOSE)
	  cerr << " Calculation complete. " << endl;

	// free all the variables	
	delete origData;  delete sumData;
	origData = NULL;  sumData = NULL;
}

void denoise (fltarray &data)
{	
    int dim = data.naxis();
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;
    
    fltarray *multiSup = NULL;
    if ((NITER > 1) && (PROGMODE == ANSC))
    {
 	    multiSup = new fltarray[NSCALE];
    }

	if (VERBOSE)
	  cerr << "Initial denoising ... " << endl;
	  
	if (PROGMODE == ANSC)
    { 
        anscombeDenoise<float>(data, multiSup);
    }
    else if (PROGMODE == FISZ)
    {
        fiszDenoise(data, MAX(0, MIN(CTRANS, nx-1)), MAX(0, MIN(CTRANS, ny-1)), MAX(0, MIN(CTRANS, nz-1)));
    } 
    	
    if (multiSup != NULL)
    {
    	if (VERBOSE)
	       cerr << "Entering into the iterative denoising ..." << endl;
        multiSupIter (*origData, data, multiSup, NITER);
	    if (VERBOSE)
	       cerr << "Iteration complete." << endl;
    }
    
	delete origData; origData = NULL; 
	if (multiSup != NULL)
	{
        delete [] multiSup; multiSup = NULL;
    }
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
   if (VERBOSE)
   {
       	cout << "Input image : " << Name_Imag_In << endl;
       	cout << "Output image : " << Name_Imag_Out << endl;
   }
   fits_read_fltarr(Name_Imag_In, *data, &header);

   header.bitpix = BP_FLOAT;
   header.origin = cmd;

   int dim = data->naxis();
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALE = MIN(NSCALE, scaleOfData(dim, nx, ny, nz));
   if ((NITER > 0) && (NITER <= 9))
        NITER = 10;
   if ((KILLLAST) && (NITER <= 9))
        NITER = 10;
   if ((DETPOS) && (NITER <= 9))
        NITER = 10;
   
   if (VERBOSE)
   {
        cout << "Iteration Filter g = " << ((TRANSF2) ? "Id-h*h" : "Id-h") << endl;
        cout << "Mode : ";
        if (PROGMODE == ANSC)
    		cout << "Anscombe" << endl;
        else if (PROGMODE == FISZ)
    		cout << "Fisz" << endl;
            
        if (FDRTHRESH)
        {
    	    if (FDRINDEP) cout << "independence mode" << endl;
       		else cout << "dependence mode" << endl;
        }
    
    	if (ITERMODE == DIRECT)
    		cout << "Iterative Mode : Direct" << endl;
    	else
    		cout << "Iterative Mode : L1 - Regularized" << endl;
    	cout << "Max. Iteration : " << NITER << endl;
    
        cout << "Max scale(s) : " << NSCALE << endl;
	    cout << "First detection scale : " << FSCALE << endl;
	    cout << "Ignore the last approx. band : " << (KILLLAST ? "true" : "false") << endl;
	    cout << "Detect only positive coefficients : " << (DETPOS ? "true" : "false") << endl;
   }
   
   // Denoising
   try {
	   denoise(*data);

       if (VERBOSE)
    	   cerr << "Writing denoising result file ... " << endl;
       fits_write_fltarr(Name_Imag_Out, *data, &header);
       if (VERBOSE)
    	 cerr << "Writing complete." << endl;

       delete data; data = NULL; 
   }
   catch (PoisMSVSTException msvstExcept)
   {
   	cerr << "Exception : " << endl;
   	cerr << msvstExcept.getReason() << endl;
   	exit (-1);
   }

   return 0;
}
