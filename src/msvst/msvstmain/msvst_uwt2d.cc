/**************************************************
 * PROJECT : PoisMSVST - version sw (1D/2D Separated wavelet Poisson detection and estiamtion)
 * CEA
 * Filename : mainsw.cc
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
#include "MSVST_Filter.h"
#include "Wavelet.h"
#include "Fisz.h"
#include "ImLib.h"
#include "Border.h"

using namespace std;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char   Name_Imag_In[256];   // input image file
char   Name_Imag_Out[256];  // output image file
const int METHOD_LEN = 3;   // total number of methods available
enum   MODE1 { MSVST = 0, ANSC = 1, FISZ = 2}; // denoising methods
bool   FDRTHRESH = false; // PV or FDR
MODE1  PROGMODE = MSVST;     // denoising mode
int    NITER = 10;          // max. number of iterations
int    CTRANS = 0;          // max. cyclic transitions
double PROBA[2] = {0.000465, 0.1}; // default cut p-value at 3.5*sigma and FDR = 0.1
int    NSCALE = 5;          // max. number of scales
bool   VERBOSE = false;     // verbose mode
bool   DEC[2] = {false, false}; // undecimated
bool   FDRINDEP = false;    // dependence parameter in FDR
int    FSCALE = 1;          // First detection scale
bool   KILLLAST = false;    // ignore the last approximation scale
bool   DETPOS   = false;    // detect only positive coefficients
double GLEVEL = 3.5;        // Gauss-detection level (k-sigma)
type_border T_BORDER = I_MIRROR; // border type

// NOTICE : FDR are all band-by-band tested

//***************************************
static void usage(char *argv[])
{
  cout << "MS-VST : 1D/2D Separated detection and estiamtion" << endl;
  cout << "Usage: " << argv[0] << " OPTIONS input_image_name output_image_name" << endl;
  cout << "OPTIONS are" << endl;  
  cout << "    [-M value]" << endl;
  cout << "         M = 0 : MS-VST (default)" << endl;
  cout << "         M = 1 : Anscombe" << endl;
  cout << "         M = 2 : Haar-Fisz" << endl;
  cout << "    [-m] Use FDR thresholding" <<endl;
  cout << "    [-t value] max. cycle translation for Haar-Fisz (default = 0)" << endl;
  cout << "    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)" << endl;
  cout << "    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)" << endl;
  cout << "    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)" << endl;
  cout << "    [-n value] number of scales (default = 5)" << endl;
  cout << "    [-F value] first detection scale (default = 1)" << endl;
  cout << "    [-K] ignore the last approximation band (used with iteration)" << endl;
  cout << "    [-p] detect only the positive coefficients (used with iteration)" << endl;
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
    while ((c = GetOpt(argc,argv,"M:t:E:s:c:n:F:i:Kmpv")) != -1) 
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
	    PROGMODE = MSVST;
	  else if (pmd == 1)
        PROGMODE = ANSC;
      else
        PROGMODE = FISZ;
	  break;
	
	case 'E':
	  double proba;
	  if (sscanf(OptArg, "%lf", &proba) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	  for (int i=0; i<2; i++) PROBA[i] = proba;
	  break;
	  
	case 's':
	  if (sscanf(OptArg, "%lf", &GLEVEL) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
      for (int i=0; i<2; i++) PROBA[i] = 2. * (1 - Utils<double>::cumNormal(GLEVEL));
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

	case 'i':
	  if (sscanf(OptArg, "%d", &NITER) != 1)
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

	case 'K': 
	  KILLLAST = true;
	  break;

	case 'm': 
	  FDRTHRESH = true;
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

  if (BORDERTYPE != I_PERIOD)
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

template <typename DATATYPE>
void msvst (to_array<DATATYPE, true> &data, to_array<DATATYPE, true> &vstData, int scale)
{
     double b, c;
     int dim = data.naxis();
     Utils<double>::antonini79VSTCoefs(dim, scale, b, c);
     
 	 int len = data.n_elem();
 	 double val;
 	 vstData.resize(data.nx(), data.ny(), data.nz());
     for (int i=0; i<len; i++)
     {
         val = data(i) + c;
		 vstData(i) = (DATATYPE) (sign(val) * b * sqrt(ABS(val)));		
     }
}

// Wavelet denoising - general process
template <typename SUPTYPE>
void waveletDenoise (fltarray &data, SubBandFilter sbf[], bool dec[], \
		             to_array<SUPTYPE, true> *multiSup = NULL)
{
	int dim = data.naxis();
	int ext[NSCALE][6];
	double fdrp;
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
	double stdg, stdhg, stdgh, stdgg;

    if ((PROGMODE==ANSC || PROGMODE == FISZ) && dim == 1)
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
            stdg = Utils<double>::antonini79SepCoefStd(s);
			double pr = FDRTHRESH ? PROBA[1] : PROBA[0];
			
            if (s < FSCALE)
            {
                setConst(cg[s-1], cg[s-1], 0);
                if (multiSup != NULL)
                    setConst(multiSup[s-1], cg[s-1], -1);  
            }
			else
			{
			  if (multiSup != NULL)
			  {
			    if (FDRTHRESH)
			    	fdrp = ws.gaussFDRThreshold (cg[s-1], stdg, pr, FDRINDEP, &multiSup[s-1]);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], pr, stdg, &multiSup[s-1]);
			  }
			  else 
			  {
			    if (FDRTHRESH)
		  	        fdrp = ws.gaussFDRThreshold (cg[s-1], stdg, pr, FDRINDEP);
			    else
			        pr = ws.gaussHardThreshold(cg[s-1], pr, stdg);
			  }	  		  

			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << (FDRTHRESH ? fdrp : pr) << endl;
			}
			mwirWT->reconstruction1D(*ch, cg[s-1], data, s, dec[0]);
			extractData(data, ext[s-1]);
			*ch = data;
		}
		delete ch; delete [] cg; 
		ch = NULL; cg = NULL; 
		if (VERBOSE) cerr << "Thresholding complete." << endl;
	}
	else if ((PROGMODE==ANSC || PROGMODE == FISZ) && dim == 2)
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
			double pr = FDRTHRESH ? PROBA[1] : PROBA[0];
            Utils<double>::antonini79SepCoefStd(s, stdhg, stdgh, stdgg);
            if (s < FSCALE)
            {
                setConst(chg[s-1], chg[s-1], 0);
                setConst(cgh[s-1], cgh[s-1], 0);
                setConst(cgg[s-1], cgg[s-1], 0);
                if (multiSup != NULL)
                {
                 setConst(multiSup[3*(s-1)], chg[s-1], -1);  
                 setConst(multiSup[3*(s-1)+1], cgh[s-1], -1);  
                 setConst(multiSup[3*(s-1)+2], cgg[s-1], -1);  
                }
            }
			else
			{
			  // use the asymptotic estimation : sigma = 1
			  if (multiSup != NULL)
			  {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chg[s-1], stdhg, pr, FDRINDEP, &multiSup[3*(s-1)]);
	  			        if (VERBOSE)
			   		       cerr << "scale(hg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgh[s-1], stdgh, pr, FDRINDEP, &multiSup[3*(s-1)+1]);
	  			        if (VERBOSE)
			   		       cerr << "scale(gh) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgg[s-1], stdgg, pr, FDRINDEP, &multiSup[3*(s-1)+2]);
	  			        if (VERBOSE)
			   		       cerr << "scale(gg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	}
			        else
			        {	
				      pr = ws.gaussHardThreshold(chg[s-1], pr, stdhg, &multiSup[3*(s-1)]);
				      pr = ws.gaussHardThreshold(cgh[s-1], pr, stdgh, &multiSup[3*(s-1)+1]);
				      pr = ws.gaussHardThreshold(cgg[s-1], pr, stdgg, &multiSup[3*(s-1)+2]);				
		  		      if (VERBOSE)
				    	cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << pr << endl;
				    }
			  }
			  else
			  {
			        if (FDRTHRESH)
			        {
			    	    fdrp = ws.gaussFDRThreshold (chg[s-1], stdhg, pr, FDRINDEP);
	  			        if (VERBOSE)
			   		       cerr << "scale(hg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgh[s-1], stdgh, pr, FDRINDEP);
    	  			    if (VERBOSE)
			   		       cerr << "scale(gh) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	    fdrp = ws.gaussFDRThreshold (cgg[s-1], stdgg, pr, FDRINDEP);
	  			        if (VERBOSE)
			   		       cerr << "scale(gg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    	}
			        else
			        {	
				      pr = ws.gaussHardThreshold(chg[s-1], pr, stdhg);
				      pr = ws.gaussHardThreshold(cgh[s-1], pr, stdgh);
				      pr = ws.gaussHardThreshold(cgg[s-1], pr, stdgg);				
				      if (VERBOSE)
				    	cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << pr << endl;
				    }
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
    else if ((PROGMODE == MSVST) && (dim == 1))
    {
	    if (VERBOSE) cerr << "Wavelet transform 1D ... " << endl;	    
	    fltarray *vstData = new fltarray;
		fltarray *chvst = new fltarray, *ch = new fltarray;
		fltarray *cg = new fltarray[NSCALE];
                 
		for (int s=1; s<=NSCALE; s++)
		{
            msvst(data, *vstData, s-1);
		    extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform1D(data, *ch, cg[s-1], s, dec[0]);
			data = *ch;
			extData(*vstData, maxlen, TB, s, dec, ext[s-1]);
			mwirWT->transform1D(*vstData, *chvst, cg[s-1], s, dec[0]);

            stdg = Utils<double>::antonini79VSTSepCoefStd(s);
            if (s < FSCALE)
            {
                if (multiSup != NULL)
                    setConst(multiSup[s-1], cg[s-1], -1);  
            }
			else if (!FDRTHRESH)
			{
			  if (multiSup != NULL)
			    {
				   ws.gaussHardThreshold(cg[s-1], PROBA[0], stdg, &multiSup[s-1]);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << PROBA[0] << endl;
			}
			else if (FDRTHRESH)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.gaussFDRThreshold(cg[s-1], stdg, PROBA[1], FDRINDEP, &multiSup[s-1]);
				if (VERBOSE)
				  cerr << "scale (g) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    }
            }
		}
		delete ch; delete [] cg; delete chvst; delete vstData;
		ch = NULL; cg = NULL; chvst = NULL; vstData = NULL;
    }    
	else if ((PROGMODE == MSVST) && (dim == 2))
	{
	    if (VERBOSE) cerr << "Wavelet transform 2D ... " << endl;
	    fltarray *vstData = new fltarray;
		fltarray *chhvst = new fltarray, *chh = new fltarray;
		fltarray *chg = new fltarray[NSCALE], \
                 *cgh = new fltarray[NSCALE], *cgg = new fltarray[NSCALE];
                 
		for (int s=1; s<=NSCALE; s++)
		{
            msvst(data, *vstData, s-1);
		    extData(data, maxlen, TB, s, dec, ext[s-1]); 
			mwirWT->transform2D(data, *chh, chg[s-1], cgh[s-1], cgg[s-1], s, dec);
			data = *chh;
			extData(*vstData, maxlen, TB, s, dec, ext[s-1]);
			mwirWT->transform2D(*vstData, *chhvst, chg[s-1], cgh[s-1], cgg[s-1], s, dec);

            Utils<double>::antonini79VSTSepCoefStd(s, stdhg, stdgh, stdgg);
            if (s < FSCALE)
            {
                if (multiSup != NULL)
                {
                 setConst(multiSup[3*(s-1)], chg[s-1], -1);  
                 setConst(multiSup[3*(s-1)+1], cgh[s-1], -1);  
                 setConst(multiSup[3*(s-1)+2], cgg[s-1], -1);  
                }
            }
			else if (!FDRTHRESH)
			{
			  if (multiSup != NULL)
			    {
				ws.gaussHardThreshold(chg[s-1], PROBA[0], stdhg, &multiSup[3*(s-1)]);
				ws.gaussHardThreshold(cgh[s-1], PROBA[0], stdgh, &multiSup[3*(s-1)+1]);
				ws.gaussHardThreshold(cgg[s-1], PROBA[0], stdgg, &multiSup[3*(s-1)+2]);
			    }
			  if (VERBOSE)
			    cerr << "scale = " << s << " cut-off p-value = " << setprecision(6) << PROBA[0] << endl;
			}
			else if (FDRTHRESH)
			{
			  if (multiSup != NULL)
			    {
				fdrp = ws.gaussFDRThreshold(chg[s-1], stdhg, PROBA[1], FDRINDEP, &multiSup[3*(s-1)]);
				if (VERBOSE)
				  cerr << "scale (hg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgh[s-1], stdgh, PROBA[1], FDRINDEP, &multiSup[3*(s-1)+1]);
				if (VERBOSE)
				  cerr << "scale (gh) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
				fdrp = ws.gaussFDRThreshold(cgg[s-1], stdgg, PROBA[1], FDRINDEP, &multiSup[3*(s-1)+2]);
				if (VERBOSE)
				  cerr << "scale (gg) = " << s << " cut-off p-value = " << setprecision(6) << fdrp << endl;
			    }
            }
		}
		delete chh; delete [] chg; delete [] cgh; delete [] cgg; delete chhvst; delete vstData;
		chh = NULL; chg = NULL; cgh = NULL; cgg = NULL; chhvst = NULL; vstData = NULL;
	}
	if (VERBOSE) cerr << "Thresholding complete." << endl;
	delete mwirWT; mwirWT = NULL;
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
  SubBandFilter sbf[3] = { SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2) };
  sbf[0].Border = sbf[1].Border = sbf[2].Border = T_BORDER;
  
  // Anscombe transform
  if (VERBOSE)
    cerr << "Anscombe transform ... " << endl;
  Utils<float>::anscombeTransform(data);	
  if (VERBOSE)
    cerr << "Anscombe transform complete. " << endl;
	
  // denoising and restoration
  waveletDenoise<float>(data, sbf, DEC);
  Utils<float>::invAnscombeTransform(data, true);
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
        SubBandFilter sbf[3] = { SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2) };
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
      data(x) = soft_threshold(data(x), lambda);
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
  if (niter <= 1) return;
  lambda = 1.;
  delta = lambda / (niter-1);
    
  fltarray *chs = new fltarray;  
  fltarray *cgs = new fltarray[NSCALE];
  fltarray *ch = new fltarray;
  fltarray *cg = new fltarray;
  
  fltarray *chhs = new fltarray;
  fltarray *chgs = new fltarray[NSCALE];
  fltarray *cghs = new fltarray[NSCALE];
  fltarray *cggs = new fltarray[NSCALE];
  fltarray *tempd = new fltarray;
  fltarray *chh = new fltarray;
  fltarray *chg = new fltarray;
  fltarray *cgh = new fltarray;
  fltarray *cgg = new fltarray;

  for (int i=0; i<niter; i++)
  {
      if (VERBOSE) cerr << "Iteration = " << (i+1) << " ... " << endl;
      *tempd = origdata;
      
      if (dim == 1)
      {
          for (int s=1; s<=NSCALE; s++)
		  {
		      extData(*tempd, maxlen, TB, s, dec, ext[s-1]);
			  extData(solution, maxlen, TB, s, dec, ext[s-1]);	
		      mwirWT->transform1D(*tempd, *ch, *cg, s, dec[0]);
			  mwirWT->transform1D(solution, *chs, cgs[s-1], s, dec[0]);
			  procArr<float>(*cg, cgs[s-1], multiSup[s-1], lambda, i);
			  solution = *chs;
			  *tempd = *ch;
		  }
	      if ((KILLLAST) && (i==niter-1))
              setConst(*chs, *ch, 0.);    
	      else
              *chs = *ch; // same approximation band
		  if (VERBOSE) cerr << "Wavelet transform 1D complete. " << endl << "Wavelet reconstruction ... " << endl;

		  for (int s=NSCALE; s>=1; s--)
		  {
		      mwirWT->reconstruction1D(*chs, cgs[s-1], solution, s, dec[0]);
		      extractData(solution, ext[s-1]); 
		      *chs = solution;
		  }
		  if (VERBOSE) cerr << "Wavelet reconstruction complete. " << endl;
      
          posProject(solution);
          lambda -= delta;
      }
      else if (dim == 2)
      {
          for (int s=1; s<=NSCALE; s++)
		  {
		      extData(*tempd, maxlen, TB, s, dec, ext[s-1]);
			  extData(solution, maxlen, TB, s, dec, ext[s-1]);	
		      mwirWT->transform2D(*tempd, *chh, *chg, *cgh, *cgg, s, dec);
			  mwirWT->transform2D(solution, *chhs, chgs[s-1], cghs[s-1], cggs[s-1], s, dec);
			  procArr<float>(*chg, chgs[s-1], multiSup[3*(s-1)], lambda, i);
			  procArr<float>(*cgh, cghs[s-1], multiSup[3*(s-1)+1], lambda, i);
			  procArr<float>(*cgg, cggs[s-1], multiSup[3*(s-1)+2], lambda, i);
			  solution = *chhs;
			  *tempd = *chh;
		  }
	      if ((KILLLAST) && (i==niter-1))
              setConst(*chhs, *chh, 0.);    
	      else
              *chhs = *chh; // same approximation band
		  if (VERBOSE) cerr << "Wavelet transform 2D complete. " << endl << "Wavelet reconstruction ... " << endl;

		  for (int s=NSCALE; s>=1; s--)
		  {
		      mwirWT->reconstruction2D(*chhs, chgs[s-1], cghs[s-1], cggs[s-1], solution, s, dec);
		      extractData(solution, ext[s-1]); 
		      *chhs = solution;
		  }
		  if (VERBOSE) cerr << "Wavelet reconstruction complete. " << endl;
      
          posProject(solution);
          lambda -= delta;
      }
  }

  delete ch; delete cg; delete chs; delete [] cgs;
  ch = NULL; cg = NULL; chs = NULL; cgs = NULL;  
  delete chhs; delete [] chgs; delete [] cghs; delete [] cggs;
  chhs = NULL; chgs = NULL; cghs = NULL; cggs = NULL;
  delete tempd; tempd = NULL;
  delete chh; delete chg; delete cgh; delete cgg;
  chh = NULL; chg = NULL; cgh = NULL; cgg = NULL;
}

void denoise (fltarray &data)
{	
    int dim = data.naxis();
	// backup the original data
	fltarray *origData = new fltarray;
	*origData = data;
    
    // design the fiters
    SubBandFilter sbf[3] = { SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2), \
			   SubBandFilter(F_MALLAT_7_9, NORM_L2) };
    sbf[0].Border = sbf[1].Border = sbf[2].Border = T_BORDER;
	
    fltarray *multiSup = NULL;
    if (NITER > 1)
    {
        if (dim == 1)
           multiSup = new fltarray[NSCALE];
        else 
 	       multiSup = new fltarray[3*NSCALE];
    }

	if (VERBOSE)
	  cerr << "Initial denoising ... " << endl;
	  
    waveletDenoise<float>(data, sbf, DEC, multiSup);
    data = *origData;
    	
    if (NITER > 1)
    {
    	if (VERBOSE)
	       cerr << "Entering into the iterative denoising ..." << endl;
        multiSupIter (*origData, data, sbf, DEC, multiSup, NITER);
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

   int dim = data->naxis();
   int nx = data->nx(), ny = data->ny(), nz = data->nz();
   NSCALE = MIN(NSCALE, scaleOfData(dim, nx, ny, nz));
   if ((PROGMODE == MSVST) && (NITER <= 9)) NITER = 10;

   if (VERBOSE)
   {
        cout << "Mode : ";
        if (PROGMODE == MSVST)
    		cout << "MSVST Denoising" << endl;
        else if (PROGMODE == FISZ)
	        cout << "Fisz Denoising, Cycle translation : " << CTRANS << endl;
   	    else if (PROGMODE == ANSC)
	        cout << "Anscombe Denoising" << endl;
            
        if (FDRTHRESH)
        {
    	    if (FDRINDEP) cout << "independence mode" << endl;
       		else cout << "dependence mode" << endl;
        }
    
   		cout << "Iterative Mode : L1 - Regularized" << endl;
	    if (PROGMODE == MSVST)
	    	cout << "Max. Iteration : " << NITER << endl;
    
        cout << "Max scale(s) : " << NSCALE << endl;
	    cout << "First detection scale : " << FSCALE << endl;
	    cout << "Ignore the last approx. band : " << (KILLLAST ? "true" : "false") << endl;
	    cout << "Detect only positive coefficients : " << (DETPOS ? "true" : "false") << endl;
   }
   
   // Denoising
   try {
       if (PROGMODE == FISZ)
     	 fiszDenoise(*data, MAX(0, MIN(CTRANS, nx-1)), MAX(0, MIN(CTRANS, ny-1)), MAX(0, MIN(CTRANS, nz-1)));
       else if (PROGMODE == ANSC)
         anscombeDenoise(*data);
       else
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
