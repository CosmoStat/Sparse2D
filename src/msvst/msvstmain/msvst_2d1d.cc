/**************************************************
 * PROJECT : PoisMSVST - MS-VST 2D+1D
 * CEA
 * Filename : main2.cc
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
#include "Atrous2D1D.h"
#include "B3VSTAtrous2D1D.h"
#include "Wavelet.h"
#include "ImLib.h"
#include "Border.h"

using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];   // input image file
char   Name_Imag_Out[256];  // output image file
char   Name_Imag_Model[256];  // output model file
char   Name_Imag_SNR[256];  // output model file
const int METHOD_LEN = 2;   // total number of methods available
enum   MODE1 { PV = 0, FDR = 1 };  // denoising modes
enum   MODE2 { DIRECT = 0, L1REG = 1 }; // iteration modes
MODE1  PROGMODE = PV;     // denoising mode
MODE2  ITERMODE = L1REG;    // iteration mode
int    NITER = 10;          // max. number of iterations
double PROBA[METHOD_LEN] = {0.000465, 0.1}; // default cut p-value at 3.5*sigma and FDR = 0.1
int    NSCALEXY = 3;          // max. number of scalexy
int    NSCALEZ = 5;           // max. number of scalez
bool   VERBOSE = false;     // verbose mode
bool   DEC[3] = {false, false, false}; // undecimated
bool   FDRINDEP = false;    // dependence parameter in FDR
bool   KILLLAST = false;    // ignore the last approximation scale
bool   DETPOS   = false;    // detect only positive coefficients
int    FSCALEXY = 1;        // First detection scalexy
int    FSCALEZ = 1;         // First detection scalez
double GLEVEL = 3.5;        // Gauss-detection level (k-sigma)
bool   TRANSF2 = false;     // ver2 in the transform
bool   MODELIM = false;     // model image
bool   WRITESNR = false;

// NOTICE : FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name model_image_name" << endl;
  cout << "MS-VST 2D+1D OPTIONS are" << endl;  
  cout << "    [-T] use g=Id-h*h as iteration band-pass filter" << endl;
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : p-value threshold (default)" << endl;
  cout << "         M = 1 : FDR threshold" << endl;
  cout << "    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)" << endl;
  cout << "    [-n value] number of scalexy (default = 3)" << endl;
  cout << "    [-N value] number of scalez (default = 5)" << endl;
  cout << "    [-F value] first detection scalexy (default = 1)" << endl;
  cout << "    [-f value] first detection scalez (default = 1)" << endl;
  cout << "    [-K] ignore the last approximation band (used with iteration)" << endl;
  cout << "    [-p] detect only the positive coefficients (used with iteration)" << endl;
  cout << "    [-I value] iteration modes for M = 0,1" << endl;
  cout << "                      I = 0 : Direct iterative mode" << endl;
  cout << "                      I = 1 : L1-regularized iterative mode (default)" << endl;
  cout << "    [-i value] number of iteration (default = 10)" << endl; 
  cout << "    [-Q value] write SNR file for every band" << endl;
  cout << "    [-v] verbose mode" << endl;
  cout << endl;
}
 
//*********************************************************************

// GET COMMAND LINE ARGUMENTS
static void filtinit(int argc, char *argv[])
{
    int c;  

    // get options 
    while ((c = GetOpt(argc,argv,"M:E:s:c:n:N:F:f:I:i:Q:TKpv")) != -1) 
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
	
	case 'Q': 
	  if (sscanf(OptArg, "%s", Name_Imag_SNR) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }	  
	  WRITESNR = true;  
	  break;

	case 'K': 
	  KILLLAST = true;
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
    
    if (OptInd < argc) 
    {
        strcpy(Name_Imag_Model, argv[OptInd++]);
        MODELIM = true;
    }

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

// denoise a band
template <typename SUPTYPE>
void bandDenoise(fltarray &band, double sigma, to_array<SUPTYPE, true> *ms, int sxy, int sz)
{
	WaveletShrinkage<float, SUPTYPE> ws;
	double pr, fdrp; // threshold p-value and FDR-threshold p-value
    char snr_filename[256];

    if (PROGMODE == PV)
    {
        pr = ws.gaussHardThreshold(band, PROBA[PROGMODE], sigma, ms);
        if (VERBOSE)
	      cerr << "(scalexy, scalez) = (" << sxy << ", " << sz << ") cut-off p-value = " << setprecision(6) << pr << endl;                 
        if (WRITESNR)
        {
              sprintf(snr_filename, "%s_XY%d_Z%d.fits", Name_Imag_SNR, sxy, sz);
              fltarray* temp = new fltarray;
              *temp = band;
              int temp_len = temp->n_elem();
              for (int idx=0; idx<temp_len; idx++)
                 (*temp)(idx) = ABS((*temp)(idx)) / sigma;
              fits_write_fltarr(snr_filename, *temp);
              delete temp;
        }  
    }
    else if (PROGMODE == FDR)
    {
        fdrp = ws.gaussFDRThreshold(band, sigma, PROBA[PROGMODE], FDRINDEP, ms);
        if (VERBOSE)
          cerr << "(scalexy, scalez) = (" << sxy << ", " << sz << ") cut-off p-value = " << setprecision(6) << fdrp << endl;                 
        if (WRITESNR)
        {
              sprintf(snr_filename, "%s_XY%d_Z%d.fits", Name_Imag_SNR, sxy, sz);
              fltarray* temp = new fltarray;
              *temp = band;
              int temp_len = temp->n_elem();
              for (int idx=0; idx<temp_len; idx++)
                 (*temp)(idx) = ABS((*temp)(idx)) / sigma;
              fits_write_fltarr(snr_filename, *temp);
              delete temp;
        }  
    }
}

// B3 Wavelet denoising - general process
template <typename SUPTYPE>
void b3SplineDenoise (fltarray &data, to_array<SUPTYPE, true> *multiSup)
{
    double sigma;
    
	// wavelet filter configuration
	B3VSTAtrous2D1D atrous;

	fltarray *cxy = new fltarray, *dxy = new fltarray;
	fltarray *ccxyz = new fltarray[NSCALEZ+1];
	fltarray *ccpxyz = new fltarray[NSCALEZ+1];
	fltarray *cdxyz = new fltarray; 
    fltarray *dcxyz = new fltarray;	
    fltarray *ddxyz = new fltarray;
    fltarray *temp = new fltarray;
    
    if (VERBOSE) cerr << "Atrous transform ... " << endl;

    // create the ccpxyz of sxy = 0
    *temp = data;
    ccpxyz[0] = *temp;
    for (int sz=1; sz<=NSCALEZ; sz++)
    {
        atrous.transformZ_Approx(*temp, ccpxyz[sz], *cdxyz, 0, sz);
        *temp = ccpxyz[sz];
    }
    
    int s = 0;
    // main loop
	for (int sxy=1; sxy<=NSCALEXY; sxy++)
	{
        // XY transform
		atrous.transformXY(data, *cxy, *dxy, sxy, 0);
		data = *cxy;
		
        // create ccxyz
        *temp = *cxy;
        ccxyz[0] = *temp;
        for (int sz=1; sz<=NSCALEZ; sz++)
        {
            atrous.transformZ_Approx(*temp, ccxyz[sz], *cdxyz, sxy, sz);
            *temp = ccxyz[sz];
    		// the cdxyz may require denoising if last scalexy
	       	if ((sxy == NSCALEXY) && ((sxy >= FSCALEXY) || (sz >= FSCALEZ)))
	       	     
	       	{
                sigma = sqrt(Utils<double>::b32D1DVSTCoefZVar (3, sxy, sz));
                bandDenoise<SUPTYPE>(*cdxyz, sigma, &multiSup[s], sxy, sz);   
               if ((sxy == NSCALEXY) && (KILLLAST))  setConst(*cdxyz, *cdxyz, 0.);    
                s++;
            }
        }
		
		// denoise of dc and dd
		for (int sz=1; sz<=NSCALEZ; sz++)
		{
    		atrous.transformZ_Detail(ccpxyz[sz-1], ccpxyz[sz], ccxyz[sz-1], ccxyz[sz], 
                                        *dcxyz, *ddxyz, sxy, sz);
            if ((sxy >= FSCALEXY) && (sz >= FSCALEZ))
            {
                // denoise dcxyz if last scalez
                if (sz == NSCALEZ)
                {
                    sigma = sqrt(Utils<double>::b32D1DVSTCoefXYVar (3, sxy, sz));
                    bandDenoise<SUPTYPE>(*dcxyz, sigma, &multiSup[s], sxy, sz);      
                    s++;
                }
                // denoise ddxyz
                sigma = sqrt(Utils<double>::b32D1DVSTCoefXYVar (3, sxy, sz-1));
                // sqrt(tau_2(g_1D)) = 0.7235; ignore correlation
                bandDenoise<SUPTYPE>(*ddxyz, 0.7235*sigma, &multiSup[s], sxy, sz);      
                s++;                
            }
        }
        
        // update ccpxyz
        for (int i=0; i<=NSCALEZ; i++) ccpxyz[i] = ccxyz[i];        
	}
	
	delete cxy; delete dxy; cxy = NULL; dxy = NULL;
	delete [] ccxyz; ccxyz = NULL;
	delete [] ccpxyz; ccpxyz = NULL;
	delete cdxyz; cdxyz = NULL;
    delete dcxyz; dcxyz = NULL;	
    delete ddxyz; ddxyz = NULL;
    delete temp; temp = NULL;
}

// for different modes of iteration
template <typename SUPTYPE>
void procArr (fltarray &origdata, fltarray &data, to_array<SUPTYPE, true> &coef, double lambda, int itr)
{
  int n = data.n_elem();

  for (int x=0; x<n; x++)
  {
      if (itr == 0)
      {
        data(x) = (coef(x) < 0) ? 0. : data(x);
        if (DETPOS)
            data(x) = (origdata(x) < 0) ? 0. : data(x);
      }
      if (DETPOS)
          data(x) = ((coef(x) >= 0) && (origdata(x) >= 0)) ? origdata(x) : data(x);
      else
          data(x) = (coef(x) >= 0) ? origdata(x) : data(x);
      if (ITERMODE == L1REG)
          data(x) = soft_threshold(data(x), lambda);
  }
}
template <typename SUPTYPE>
void procArr (fltarray &data, double lambda, int itr)
{
  int n = data.n_elem();

  for (int x=0; x<n; x++)
  {
     if (itr == 0)
     {
        data(x) = 0.;
     }
     if (ITERMODE == L1REG)
        data(x) = soft_threshold(data(x), lambda);
  }
}

// put to zeros the negative values
void posProject (fltarray &data)
{	
	int len = data.n_elem();
	
	for (int i=0; i<len; i++)
	{
		data(i) = MAX(data(i), 0);
	}
}

// solve with the iterative algo. in using the multiresolution support
void multiSupIter (fltarray &origdata, fltarray &solution, fltarray *model, fltarray *multiSup, int niter)
{
  Atrous2D1D atrous(I_MIRROR, TRANSF2);
  double lambda = 0, delta = 0;
  
  if (ITERMODE == L1REG)
  {
      if (niter <= 1) return;
      lambda = 1.;
      delta = lambda / (niter-1);
  }
    
  if (MODELIM)
  {
        solution -= *model;
  }

  setConst(solution, solution, 0.);
  
  fltarray *chxys = new fltarray;
  fltarray *cgxys = new fltarray[NSCALEXY];
  fltarray *chzs = new fltarray;
  fltarray *cgzs = new fltarray[NSCALEZ];
  
  fltarray *tempd = new fltarray;
  fltarray *chxy = new fltarray;
  fltarray *cgxy = new fltarray;
  fltarray *chz = new fltarray;
  fltarray *cgz = new fltarray;

  for (int i=0; i<niter; i++)
  {
      int s = 0;
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;
      *tempd = origdata;
      if (MODELIM) *tempd -= *model;
      
      // decompose XY
	  for (int sxy=1; sxy<=NSCALEXY; sxy++)
	  {
	    atrous.transformXY(*tempd, *chxy, *cgxy, sxy, 0);
		atrous.transformXY(solution, *chxys, cgxys[sxy-1], sxy, 0);
		// process cxy_dz if last scalexy
		if (sxy == NSCALEXY)
		{
            // decompose Z of cxy (last scalexy)
    		for (int sz=1; sz<=NSCALEZ; sz++)
	       	{
	       		
    	        atrous.transformZ(*chxy, *chz, *cgz, sxy, sz);
    	       	atrous.transformZ(*chxys, *chzs, cgzs[sz-1], sxy, sz);
    	       	
    	       	if (!KILLLAST)
    	       	{
    	       	if ((sxy >= FSCALEXY) && (sz >= FSCALEZ))    	       	
    	       	{
            		procArr<float>(*cgz, cgzs[sz-1], multiSup[s], lambda, i);
        	       	s++;
                }
            	else
            		procArr<float>(cgzs[sz-1], lambda, i);
    	       	}
        		*chxy = *chz;
        		*chxys = *chzs;
            }
            // if ((KILLLAST) && (i==0))  setConst(*chzs, *chz, 0.);
            // if (! KILLLAST)  *chzs = *chz; // same approximation band
            
            if (!KILLLAST)  *chzs = *chz;
            
            //if ((KILLLAST) && (i==niter-1))
            //    setConst(*chzs, *chz, 0.);    
      	    // else
            //      *chzs = *chz; // same approximation band
            
                 
            // reconstruct Z of cxy       
            for (int sz=NSCALEZ; sz>=1; sz--)
            {
	          	atrous.reconsZ(*chzs, cgzs[sz-1], *chxys, sxy, sz);
    	       	*chzs = *chxys;
            }            
        }
		
		// decompose Z of dxy
		for (int sz=1; sz<=NSCALEZ; sz++)
		{
    	    atrous.transformZ(*cgxy, *chz, *cgz, sxy, sz);
	       	atrous.transformZ(cgxys[sxy-1], *chzs, cgzs[sz-1], sxy, sz);
	       	// process dxy_cz if last scalez
	       	if ((sxy != NSCALEXY) || (!KILLLAST))
	       	{
	       	if (sz == NSCALEZ)
	       	{
                if ((sxy >= FSCALEXY) && (sz >= FSCALEZ))    	  
                {
        		    procArr<float>(*chz, *chzs, multiSup[s], lambda, i);
        		    s++;
                }
                else
                    procArr<float>(*chzs, lambda, i);
            }
            if ((sxy >= FSCALEXY) && (sz >= FSCALEZ))    	  
            {
        		procArr<float>(*cgz, cgzs[sz-1], multiSup[s], lambda, i);
        		s++;
            }
            else
                procArr<float>(cgzs[sz-1], lambda, i);
	       	}
    		*cgxy = *chz;
    		cgxys[sxy-1] = *chzs;
        }
        // reconstruct Z of dxy
        for (int sz=NSCALEZ; sz>=1; sz--)
        {
	       	atrous.reconsZ(*chzs, cgzs[sz-1], cgxys[sxy-1], sxy, sz);
    		*chzs = cgxys[sxy-1];
        }
        
		*tempd = *chxy;
		solution = *chxys;
	  }
	  // reconstruct XY
	  for (int sxy=NSCALEXY; sxy>=1; sxy--)
	  {
	    atrous.reconsXY(*chxys, cgxys[sxy-1], solution, sxy, 0);
	    *chxys = solution;
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
  
  delete chxys; delete [] cgxys; chxys = NULL; cgxys = NULL;
  delete chzs; delete [] cgzs; chzs = NULL; cgzs = NULL;
  delete tempd; tempd = NULL;
  delete chxy; delete cgxy; chxy = NULL; cgxy = NULL;
  delete chz; delete cgz; chz = NULL; cgz = NULL;
}

void denoise (fltarray &data, fltarray *model)
{	
    int dim = data.naxis();
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;

 	fltarray *multiSup = new fltarray[NSCALEXY*(NSCALEZ+1) + NSCALEZ];

	if (VERBOSE)
	  cerr << "Initial denoising ... " << endl;
	  
    b3SplineDenoise<float>(data, multiSup);
	data = *origData;
	
	if (VERBOSE)
	  cerr << "Entering into the iterative denoising ..." << endl;
    multiSupIter (*origData, data, model, multiSup, NITER);
	if (VERBOSE)
	  cerr << "Iteration complete." << endl;

	delete origData; origData = NULL; 
    delete [] multiSup; multiSup = NULL;
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
       	if (MODELIM) cout << "Model image : " << Name_Imag_Model << endl;
   }
   fits_read_fltarr(Name_Imag_In, *data, &header);
   if (MODELIM) 
   {
        model = new fltarray;
        fits_read_fltarr(Name_Imag_Model, *model, &header);
   }

   header.bitpix = BP_FLOAT;
   header.origin = cmd;

   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALEXY = MIN(NSCALEXY, scaleOfData(2, nx, ny, nz));
   NSCALEZ  = MIN(NSCALEZ, scaleOfData(1, nz, nz, nz));
   // if (NITER < 10) NITER = 10;
   
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
    
        cout << "Max scalexy(s) : " << NSCALEXY << endl;
        cout << "Max scalez(s) : " << NSCALEZ << endl;
	    cout << "First detection scalexy : " << FSCALEXY << endl;
	    cout << "First detection scalez : " << FSCALEZ << endl;
	    cout << "Ignore the last approx. band : " << (KILLLAST ? "true" : "false") << endl;
	    cout << "Detect only positive coefficients : " << (DETPOS ? "true" : "false") << endl;
  }
   
   // Denoising
   try {
	   denoise(*data, model);

       if (VERBOSE)
    	   cerr << "Writing denoising result file ... " << endl;
       fits_write_fltarr(Name_Imag_Out, *data, &header);
       if (VERBOSE)
    	 cerr << "Writing complete." << endl;

       if (MODELIM)
       {
            delete model; model = NULL;
       }
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
