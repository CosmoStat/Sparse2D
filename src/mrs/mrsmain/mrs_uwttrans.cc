/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  / Florent Sureau
**
**    Date:  25/10/08 and 27/09/11
**    
**    File:  fcs_mrs_uwttrans.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Isotropic Undecimated Wavelet Transform on the sphere (Through Spherical Harmonics)
**    -------------------
**
******************************************************************************
**
**    Modifications from mrs_uwttrans.c : Original version only for testing that
**		forward and inverse transform were matching. 
**
**      27/09/11 : FCS : added Lmax input, PlanckFlag and InverseFlag to allow 
**		separate forward and backward transform , BandLimit to enforce or not 
**  	bandlimited signals, Square root filters for tight frame decomposition,
**		and choice of number of wavelet scales (not the full number).
**
******************************************************************************/

#include"MRS_Sparse.h"
#include <fstream>
using namespace std;

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
  
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

bool Verbose=false;
bool InverseFlag=false;
bool PlanckFlag=false;
bool BandLimit=false;
bool SqrtFilters=false;
bool TightFrame=false;
bool UseMeyer=false;
float SigmaNoise=0.;
 
int NbrScale=0;
int Lmax=0;
int ALM_iter=0;

/***************************************************************************/

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  Default is automatically estimated\n");
    fprintf(OUTMAN, "         [-M\n");
    fprintf(OUTMAN, "             Use Meyer wavelets instead of Spline wavevelets.  \n");
    fprintf(OUTMAN, "             Number of scales is automatically estimated according to Nside parameter.  \n");
    fprintf(OUTMAN, "         [-T]\n");
    fprintf(OUTMAN, "             Use square root filters which ensure having a tight frame (better for inverse problems).  \n");
    fprintf(OUTMAN, "             Default is taking instead full filters, i.e. only summing wavelet coefficients for reconstruction.  \n");
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Enforce band limiting of the signal to Lmax (finest scale).  \n");
    fprintf(OUTMAN, "             Do not consider signal at multipoles > Maximal_multipole (option -l).  \n");
    fprintf(OUTMAN, "         [-a N_ITER_ALM]\n");
    fprintf(OUTMAN, "             Number of iterations N_ITER_ALM for iterative spherical harmonics transform.  \n");
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Maximal multipole for spherical harmonic transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated according to the input map nside: min(3*Nside,%d). \n", ALM_MAX_L);
    fprintf(OUTMAN, "         [-r\n");
    fprintf(OUTMAN, "             Perform inverse wavelet transform instead of forward transform.  \n");
    fprintf(OUTMAN, "             In that case, option -n is not taken into account. Other options (i.e. -T,-M,-l,-p) used for the transformation should be set as well for the reconstruction. \n");
    fprintf(OUTMAN, "         [-p\n");
    fprintf(OUTMAN, "             Use wavelet filters defined for Planck.  \n");
    fprintf(OUTMAN, "   Note that the output map containing the wavelet decomposition is in RING format.\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "Mn:brTpa:l:vzZ")) != -1)
    {
	switch (c)  { 
         case 'M': UseMeyer = true; break; 
         case 'a':
 	     	  if (sscanf(OptArg,"%d",&ALM_iter) != 1)  
               {
		           fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
		            exit(-1); 
	            	}
 		          break;
        
        case 'n':
 	     	  if (sscanf(OptArg,"%d",&NbrScale) != 1)  
                  {
		            fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
		            exit(-1); 
	            	}
 		          break;
	    case 'l':
 	     	  if (sscanf(OptArg,"%d",&Lmax) != 1)  {
		            fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
		            exit(-1); 
	          }
 		  break; 
	    case 'b':
 	      BandLimit=true;
 		  break;
 		case 'T':
            SqrtFilters=true;TightFrame=true;
 		  break;
 		case 'p':
 	     	  PlanckFlag=true;
 		  break;
	    case 'r':
 	     	 InverseFlag=true;
 		  break;
     	    case 'v': Verbose = true; break;
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
/*********************************************************************/


int main(int argc, char *argv[])
{
    int k;
    double Min, Max;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (NbrScale > 0)  cout << "# NbrScale = " <<  NbrScale << endl;
        if (TightFrame==true) cout << "Use Sqrt Filters (Tight Frame)" << endl;
     }

   int Nside;
   
   C_UWT WT;
   WT.All_WP_Band=PlanckFlag; //Use Planck Wavelet Packets or not (see MRS_AlmBand.cc)
   WT.set_alm_iter(ALM_iter); 
   if(InverseFlag ==false)
   {
		Hdmap Map;   
		Map.read(Name_Imag_In);   
   		Nside = Map.Nside();
    	// if(Lmax == 0) Lmax=min(2*Nside,ALM_MAX_L);//default value for Lmax
        int LMAX = mrs_get_lmax (Lmax,  Nside, 0.0);
        Lmax =LMAX;
        WT.Verbose = Verbose;
  		if (UseMeyer == true)  WT.wp_alloc(Nside, Lmax, DEF_MRS_ORDERING);//Ring Format (quicker spherical harmonic transform)
        else WT.wt_alloc(Nside, NbrScale, Lmax, DEF_MRS_ORDERING, TightFrame);
		if (Verbose == true){
		   cout << "Forward Transform" << endl;
   		   Map.minmax(Min,Max);
   		   // Map.info((char*) "Input MAP");
   		   if (Verbose == true) 
           {
               cout << "Lmax=" << Lmax <<endl;
               cout << "NbrScale = " <<  WT.nscale() << endl;
               if (TightFrame==true) cout << "Use Sqrt Filters (Tight Frame)" << endl;
               cout << "BandLimit = " <<  BandLimit << endl;
           }
   		} 
 		WT.transform(Map,BandLimit,SqrtFilters,NbrScale); 
 		if(NbrScale <=0) NbrScale=WT.nscale();
		fits_write_dblarr(Name_Imag_Out,WT.WTTrans);
       // fits_write_fltarr("xx_wpfilter.fits", WT.WP_WPFilter);
	}
   else
   {
		Hdmap MapOut;
		dblarray InData;
        fits_read_dblarr(fitsname(Name_Imag_In),InData);
        cout << "Input: " <<InData.nx() << " " << InData.ny() << endl;

		Nside=sqrt(InData.nx()/12l);
		NbrScale=InData.ny();
		if(Nside <= 0) {
			printf("Problem with input wavelet transform coefficients: Nside<=0\n");
			exit(1);	
		}
        Lmax = mrs_get_lmax (Lmax,  Nside, 0.0);
		// if(Lmax == 0) Lmax=min(2*Nside,ALM_MAX_L);//default value for Lmax

       if (UseMeyer == true)  WT.wp_alloc(Nside, Lmax, DEF_MRS_ORDERING);//Ring Format (quicker spherical harmonic transform)
       else WT.wt_alloc(Nside, NbrScale, Lmax, DEF_MRS_ORDERING, TightFrame);
       
		WT.WTTrans=InData;
		if (Verbose == true) { 
 			cout << "Inverse Transform" << endl;
			cout << "NScales=" << NbrScale << endl;
   		    cout << "Lmax=" << Lmax <<endl;
            cout << "BandLimit = " <<  BandLimit << endl;
            if (TightFrame==true) cout << "Use Sqrt Filters (Tight Frame)" << endl;
		}
		
   		MapOut.alloc(Nside, DEF_MRS_ORDERING);
   		MapOut.fill(0.);
       
   		if ((UseMeyer == false) || (SqrtFilters==false)) WT.recons(MapOut,BandLimit);
   		else  WT.recons(MapOut,BandLimit,SqrtFilters,NbrScale-1);
   		MapOut.write(Name_Imag_Out);
	}
  
   exit(0);
}

