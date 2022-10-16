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
char  Name_Mask_Out[512];

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

bool Verbose=false;
bool InverseFlag=false;
bool PlanckFlag=false;
bool BandLimit=false;
bool SqrtFilters=false;
bool UseMeyer=false;
float SigmaNoise=0.;
 
int NbrScale=0;
int Lmax=0;
int ALM_iter=0;

float NSigma=5.;
bool UseMad=true;
bool CMBDel = true;
bool KillLastScale = true;
bool UseMask=false;

/***************************************************************************/

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map [out_mask]  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Maximal multipole for spherical harmonic transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated according to the input map nside: min(3*Nside,%d). \n", ALM_MAX_L);
    
    fprintf(OUTMAN, "         [-s Nsigma]\n");
    fprintf(OUTMAN, "             Nsigma Detection level. Default is 5.\n");
    

    fprintf(OUTMAN, "   Note that the output map is in RING format.\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "n:s:l:vzZ")) != -1) 
    {
	switch (c)  { 
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
  		case 'g':
            if (sscanf(OptArg,"%f",&SigmaNoise) != 1)  
            {
                fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                exit(-1); 
            }
            break;
   		case 's':
            if (sscanf(OptArg,"%f",&NSigma) != 1)  
            {
                fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                exit(-1); 
            }
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

    if (OptInd < argc) 
    {
       strcpy(Name_Mask_Out, argv[OptInd++]);
       UseMask = true;
    }

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
    }

   int Nside;
   
   C_UWT WT;
   WT.set_alm_iter(ALM_iter); 
    
   Hdmap Map, MapOut;   
   Map.read(Name_Imag_In);   
   Nside = Map.Nside();
   int Lmax = mrs_get_lmax (Lmax,  Nside, 0.0);
   if (NbrScale < 2)  NbrScale = log((double) Nside) / log((double) 2.);
     
   WT.wt_alloc(Nside, NbrScale, Lmax, DEF_MRS_ORDERING);
   
   if (Verbose == true)
    {
        Map.minmax(Min,Max);
        if (Verbose == true) 
        {
            Map.info((char*) "Input MAP");
            cout << "Lmax=" << Lmax <<endl;
            cout << "NbrScale = " <<  WT.nscale() << endl;
        }
    } 
    
    MapOut.alloc(Nside, DEF_MRS_ORDERING);
    MapOut = Map;
    WT.hard_thresholding( MapOut, NSigma, SigmaNoise, UseMad, KillLastScale);
    for (int p=0; p < MapOut.Npix(); p++) Map[p]-= MapOut[p];		
    Map.write(Name_Imag_Out);
    
    if (UseMask == true)
    {
        Hdmap Mask;   
        Mask.alloc(Nside, DEF_MRS_ORDERING);
        for (int p=0; p < MapOut.Npix(); p++) 
        {
           if (ABS( MapOut[p] ) > 0) Mask[p]=0;
           else Mask[p]=1;
        }
        Mask.write(Name_Mask_Out);
    }
   exit(0);
}

