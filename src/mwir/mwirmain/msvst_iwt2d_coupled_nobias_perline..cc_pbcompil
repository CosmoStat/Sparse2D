/**************************************************
 * PROJECT : PoisMSVST - Corrected Coupled detection and estimation of each line of a 2D image
 * CEA
 * Filename : mainlinecc.cc
 * This is the main file of PoisMSVST
 **************************************************/
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "cdflib.h"
#include "GlobalInc.h"
#include "IM_IO.h"
#include "Array.h"
#include "PoisMSVSTException.h"
#include "Atrous.h"
#include "B3VSTAtrous.h"
#include "Wavelet.h"
#include "ImLib.h"
#include "Border.h"

using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];   // input image file
char   Name_Imag_Out[256];  // output image file
const int METHOD_LEN = 2;   // total number of methods available
enum   MODE1 { PV = 0, FDR = 1 };  // denoising modes
enum   MODE2 { DIRECT = 0, L1REG = 1 }; // iteration modes
MODE1  PROGMODE = PV;     // denoising mode
MODE2  ITERMODE = L1REG;    // iteration mode
int    NITER = 10;          // max. number of iterations
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
bool   MODELIM = false;     // model image
bool   POSPROJ = true;
bool   INPAINT = false;     // inpainting mode
bool   CBIAS   = true;     // correct the VST's bais
bool   INTENSITY = true; // intensity or scale mode
// NOTICE : FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "MS-VST : Coupled detection and estimation of each line of a 2D image" << endl;
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-T] use g=Id-h*h as iteration band-pass filter" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : p-value threshold (default)" << endl;
  cout << "         M = 1 : FDR threshold" << endl;
  cout << "    [-S] scale mode instead of intensity mode" <<endl;  
  cout << "    [-b] do NOT correct the VST's bias" <<endl;  
  cout << "    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)" << endl;
  cout << "    [-n value] number of scales (default = 5)" << endl;
  cout << "    [-F value] first detection scale (default = 1)" << endl;
  cout << "    [-K] ignore the last approximation band (used with iteration)" << endl;
  cout << "    [-p] detect only the positive coefficients (used with iteration)" << endl;
  cout << "    [-N] do NOT do positivity projection" << endl;
//  cout << "    [-I value] iteration modes for M = 0,1" << endl;
//  cout << "                      I = 0 : Direct iterative mode" << endl;
//  cout << "                      I = 1 : L1-regularized iterative mode (default)" << endl;
  cout << "    [-i value] number of iteration (default = 10)" << endl; 
  cout << "    [-v] verbose mode" << endl;
  cout << endl;
}
 
//*********************************************************************

// GET COMMAND LINE ARGUMENTS
static void filtinit(int argc, char *argv[])
{
    int c;  

    // get options 
    while ((c = GetOpt(argc,argv,"M:E:s:c:n:F:i:TKbpNSv")) != -1) 
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
	    PROGMODE = PV;
	  else 
        PROGMODE = FDR;
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

//	case 'I':
//	  int imd;
//	  if (sscanf(OptArg, "%d", &imd) != 1)
//	    {
//	      cerr << "Bad or missing parameter " << OptArg << endl;
//	      exit (-1);
//	    }
//	  if (imd == 0)
//	    ITERMODE = DIRECT;
//	  else 
//	    ITERMODE = L1REG;
//	  break;
	
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

	case 'T': 
	  TRANSF2 = true;
	  break;

	case 'b': 
	  CBIAS = false;
	  break;

	case 'p': 
	  DETPOS = true;
	  break;

	case 'N': 
	  POSPROJ = false;
	  break;

	case 'S': 
	  INTENSITY = false;
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
	int nx = ref.nx(), ny = ref.ny(), nz = ref.nz();
	int len = ref.n_elem();

	dest.resize(nx, ny, nz);
    for (int i=0; i<len; i++)
		dest(i) = (DATATYPE)cst;		
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void b3SplineDenoise (fltarray &data, to_array<SUPTYPE, true> *multiSup)
{
	int dim = data.naxis();
	double pr, fdrp; // threshold p-value and FDR-threshold p-value
	int N = data.n_elem();
    bool ms = (multiSup != NULL);
    
	// wavelet filter configuration
	B3VSTAtrous atrous;
	WaveletShrinkage<float, SUPTYPE> ws;

	fltarray *ch = new fltarray, *chvst = new fltarray, *inttemp = new fltarray;;
	fltarray *cg = new fltarray[NSCALE];
    
    if (INTENSITY)
    {
	    for (int s=1; s<=NSCALE; s++)
	    {
		    atrous.transform(data, *ch, *chvst, cg[s-1], s);
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
    		else if (PROGMODE == PV)
    		{
    		  if (VERBOSE)
                  cerr << "Var = " << Utils<double>::b3VSTCoefVar (dim, s) << endl;
              pr = ws.gaussHardThreshold(cg[s-1], PROBA[PROGMODE], 
                                         sqrt(Utils<double>::b3VSTCoefVar (dim, s)), 
                                         ms ? &multiSup[s-1] : NULL);
    		  if (VERBOSE)
    		    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << pr << endl;                 
            }
    		else if (PROGMODE == FDR)
    		{
    		  if (VERBOSE)
                  cerr << "Var = " << Utils<double>::b3VSTCoefVar (dim, s) << endl;
              fdrp = ws.gaussFDRThreshold(cg[s-1], sqrt(Utils<double>::b3VSTCoefVar (dim, s)), 
                                          PROBA[PROGMODE], FDRINDEP, ms ? &multiSup[s-1] : NULL);
    		  if (VERBOSE)
    		    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;                 
            }
    	}
    	for (int s=NSCALE; s>=1; s--)
    	{
            atrous.recons(*chvst, cg[s-1], *ch, s);
            *chvst = *ch;
        }
        atrous.InvMSVST(*chvst, data, 0, CBIAS);
    }
    else
    {
    	for (int s=1; s<=NSCALE; s++)
    	{
    		atrous.transform(data, *ch, *chvst, cg[s-1], s);
    		data = *ch;
        }
    	for (int s=NSCALE; s>=1; s--)
        {	    
            if (s < FSCALE)
            {
                setConst(cg[s-1], data, 0);
                if (ms)
                {
                    multiSup[s-1].resize(cg[s-1].nx(), cg[s-1].ny(), cg[s-1].nz());
                    for (int i=0; i<N; i++) (multiSup[s-1])(i) = -1;
                }
            }		
    		else if (PROGMODE == PV)
    		{
              pr = ws.ccHardThreshold(cg[s-1], *ch, s, PROBA[PROGMODE],                                       
                                         ms ? &multiSup[s-1] : NULL);
    		  if (VERBOSE)
    		    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << pr << endl;                 
            }
    		else if (PROGMODE == FDR)
    		{
              fdrp = ws.ccFDRThreshold(cg[s-1], *ch, s, PROBA[PROGMODE], FDRINDEP, ms ? &multiSup[s-1] : NULL);
    		  if (VERBOSE)
    		    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;                 
            }
    	
    	    *inttemp = *ch;
            atrous.recons(*chvst, cg[s-1], *ch, s);
            *chvst = *ch;
            atrous.InvMSVST(*chvst, *ch, s-1, CBIAS, inttemp);
        }
        
        data = *ch;
    }
	delete ch; delete chvst; delete inttemp; delete [] cg; 
    ch = NULL; inttemp = NULL; chvst = NULL; cg = NULL;
}

// for different modes of iteration
template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef, double lambda)
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

template <typename SUPTYPE>
void procArr (fltarray &data, double lambda)
{
  int n = data.n_elem();

  for (int x=0; x<n; x++)
  {
      if (ITERMODE == L1REG)
           data(x) = soft_threshold(data(x), lambda);
  }
}

template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef)
{
  int n = data.n_elem();

  for (int x=0; x<n; x++)
  {
      if (DETPOS)
          data(x) = ((coef(x) >= 0) && (origdata(x) >= 0)) ? data(x) : 0.;
      else
          data(x) = (coef(x) >= 0) ? data(x) : 0.;
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
void multiSupIter (fltarray &origdata, fltarray &solution, fltarray *model, fltarray *mask, fltarray *multiSup, int niter)
{
  int dim = origdata.naxis();
  SplineAtrous atrous = SplineAtrous(3, 1., I_MIRROR, TRANSF2);
  double lambda = 0, delta = 0;
  
  if (MODELIM)
  {
        solution -= *model;
  }
  
  if (!INPAINT)
  {  
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
          if (MODELIM) *tempd -= *model;
          
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
          
          if (MODELIM)
          {
                solution += *model;
          }
          
          posProject(solution);
          
          if (MODELIM)
          {
                solution -= *model;
          }    
          lambda -= delta;
      }
      
      delete chs; delete [] cgs; chs = NULL; cgs = NULL;
      delete tempd; tempd = NULL;
      delete ch; delete cg; ch = NULL; cg = NULL;
  }
  else
  {
      if (niter <= 1) return;

      fltarray *fsolution = new fltarray;
      fltarray *residual = new fltarray;
      fltarray *chr = new fltarray;
      fltarray *cgr = new fltarray[NSCALE];
      fltarray *ch = new fltarray;
      fltarray *cg = new fltarray;
      fltarray *tempd = new fltarray;

      for (int i=0; i<niter; i++)
      {
          if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;
          
          // reinforce the detections
          *residual = origdata;
          *fsolution = solution;
          if (MODELIM) *fsolution += *model;
          *fsolution *= *mask;
          *residual -= *fsolution;
          
          *tempd = origdata;
       	  for (int s=1; s<=NSCALE; s++)
    	  {
    	    atrous.transform(*tempd, *ch, *cg, s);
    		atrous.transform(*residual, *chr, cgr[s-1], s);
    		procArr<float>(*cg, cgr[s-1], multiSup[s-1]);
    		*residual = *chr;
    		*tempd = *ch;
    	  }
    	  for (int s=NSCALE; s>=1; s--)
    	  {
    	    atrous.recons(*chr, cgr[s-1], *residual, s);
    	    *chr = *residual;
          }
    	  solution += *residual;
 
          // soft thresholding
          *tempd = origdata;
          if (MODELIM) *tempd -= *model;
    	  for (int s=1; s<=NSCALE; s++)
    	  {
    	    atrous.transform(*tempd, *ch, *cg, s);
    		atrous.transform(solution, *chr, cgr[s-1], s);
    		if (i==0)
    		{
                lambda = MAX(lambda, ABS(cgr[s-1].max())+0.1); // in case of numerical errors
                lambda = MAX(lambda, ABS(cgr[s-1].min())+0.1); // in case of numerical errors
                delta = lambda / (niter-1);
            }
    		procArr<float>(cgr[s-1], lambda);
    		solution = *chr;
    		*tempd = *ch;
    	  }
	      if ((KILLLAST) && (i==niter-1))
            setConst(*chr, *ch, 0.);    
    	  else
           *chr = *ch; // same approximation band    
    	  for (int s=NSCALE; s>=1; s--)
    	  {
    	    atrous.recons(*chr, cgr[s-1], solution, s);
    	    *chr = solution;
          }
          
          if (MODELIM)
          {
                solution += *model;
          }
          
          posProject(solution);
          
          if (MODELIM)
          {
                solution -= *model;
          }    
          lambda -= delta;
      }
 
      delete fsolution; fsolution = NULL;
      delete residual; residual = NULL;
      delete chr; chr = NULL;
      delete [] cgr; cgr = NULL;
      delete ch; ch = NULL;
      delete cg; cg = NULL;
      delete tempd; tempd = NULL;
  }
}

void denoise (fltarray &data)
{	
    int dim = data.naxis();
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;
    
    fltarray *multiSup = NULL;
    if (NITER > 1)
    {
 	    multiSup = new fltarray[NSCALE];
    }

	if (VERBOSE)
	  cerr << "Initial denoising ... " << endl;
	  
    b3SplineDenoise<float>(data, multiSup);
    posProject(data);
    	
    if (NITER > 1)
    {
    	if (VERBOSE)
	       cerr << "Entering into the iterative denoising ..." << endl;
        multiSupIter (*origData, data, NULL, NULL, multiSup, NITER);
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

   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALE = MIN(NSCALE, scaleOfData(1, nx, ny, nz)); // one dimension
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
        if (PROGMODE == PV)
    		cout << "Individual tests" << endl;
        else if (PROGMODE == FDR)
    		cout << "FDR multiple tests" << endl;
            
        if (PROGMODE == FDR)
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
        fltarray dataline(nx);
        for (int y=0; y<ny; y++)
        {
            for (int x=0; x<nx; x++)
            {
                dataline(x) = (*data)(x, y);
            }
       	    denoise(dataline);
            for (int x=0; x<nx; x++) (*data)(x, y) = dataline(x);
        }

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
