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

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
  
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

bool Verbose=false;
bool InverseFlag=false;
bool PlanckFlag=false;
bool BandLimit=true;
bool SqrtFilters=false;
float SigmaNoise=0.;
int ALM_iter=0;

int NbrScale=0;
int Lmax=0;

/***************************************************************************/

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated according to Nside parameter.  \n");
    fprintf(OUTMAN, "         [-s Square_root_filter]\n");
    fprintf(OUTMAN, "             Use square root filters which ensure having a tight frame (better for inverse problems).  \n");
    fprintf(OUTMAN, "             Default is taking instead full filters, i.e. only summing wavelet coefficients for reconstruction (advantage: 1 spherical harmonic transform against 2 for sqrt filters.  \n");
    fprintf(OUTMAN, "         [-b BandLimit_Mode]\n");
    fprintf(OUTMAN, "             Do not enforce band limiting of the signal to Lmax (finest scale).  \n");
    fprintf(OUTMAN, "             Default is taking no signal at multipoles > Maximal_multipole (option -l).  \n");
    fprintf(OUTMAN, "         [-a N_ITER_ALM]\n");
    fprintf(OUTMAN, "             Number of iterations N_ITER_ALM for iterative spherical harmonics transform.  \n");
    fprintf(OUTMAN, "         [-l Maximal_multipole]\n");
    fprintf(OUTMAN, "             Maximal multipole for spherical harmonic transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated according to the input map nside: min(2*Nside,3000).  \n");
    fprintf(OUTMAN, "         [-i\n");
    fprintf(OUTMAN, "             Perform inverse wavelet transform instead of forward transform.  \n");
    fprintf(OUTMAN, "             In that case, option -n is not taken into account.  \n");
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
    while ((c = GetOpt(argc,argv,(char *) "a:n:bispl:vzZ")) != -1) 
    {
	switch (c)  {
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
 	      BandLimit=false;
 		  break;
 		case 's':
 	      SqrtFilters=true;
 		  break;
 		case 'p':
 	     	  PlanckFlag=true;
 		  break;
	    case 'i':
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
        if(SqrtFilters==true) cout << "Use Sqrt Filters" << endl;
     }

   int Nside;
   
   C_PWT WT;
   WT.All_WP_Band=PlanckFlag; //Use Planck Wavelet Packets or not (see MRS_AlmBand.cc)
   WT.set_alm_iter(ALM_iter); 
   if(InverseFlag ==false)  {
		Hdmap Map;   
		Map.read(Name_Imag_In);   
   		Nside = Map.Nside();
    	if(Lmax == 0) Lmax=min(2*Nside,ALM_MAX_L);//default value for Lmax
  		WT.wp_alloc(Nside, Lmax, DEF_MRS_ORDERING);//Ring Format (quicker spherical harmonic transform)
		if (Verbose == true){
		   cout << "Forward Transform" << endl;
   		   Map.minmax(Min,Max);
   		   Map.info((char*) "Input MAP");
   		   cout << "Lmax=" << Lmax <<endl;
   		   if (NbrScale <= 0)  cout << "NbrScale = " <<  WT.nscale() << endl;
   		} 
 		WT.transform(Map,BandLimit,SqrtFilters,NbrScale); 
 		if(NbrScale <=0) NbrScale=WT.nscale();
		if(NbrScale<=0) {
			printf("NbrScale should be >0\n");
			exit(-1);
		} 
  		WT.write_scales(Name_Imag_Out,NbrScale);
	} else {
		Hdmap MapOut;
		WT.read_scales(Name_Imag_In,Lmax,NbrScale,DEF_MRS_ORDERING);
		if (Verbose == true) { 
 			cout << "Inverse Transform" << endl;
			cout << "NScales=" << NbrScale << endl;
   		    cout << "Lmax=" << Lmax <<endl;
		}
   		MapOut.alloc(WT.nside(), DEF_MRS_ORDERING);
   		MapOut.fill(0.);
   		WT.recons(MapOut,BandLimit,SqrtFilters,NbrScale);
   		MapOut.write(Name_Imag_Out);
	}
  
   exit(0);
}

