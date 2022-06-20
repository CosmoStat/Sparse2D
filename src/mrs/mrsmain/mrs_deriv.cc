/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  28/11/08
**    
**    File:  mrs_deriv.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Compute the derivatives of a sperical image in theta and phi
**    -----------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out1[512]; /* output file name */
char Name_Imag_Out2[512]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool OptFilter = False;
float Fwhm=0.;
Bool NormALM = False;
float  ZeroPadding=0.;
int Lmax = 0.;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map out_map_der_theta  out_map_der_phi\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-F Filter_FileName]\n");
    fprintf(OUTMAN, "             Filter file name . \n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Full width at half maximum (in arc minute). Default is 60 arc minutes.\n");
	
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");

    exit(-1);
}
 
/*********************************************************************/
 
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "f:l:p:vzZ")) != -1) 
    {
	switch (c) 
        { 
          case 'f':
 	     	      if (sscanf(OptArg,"%f",&Fwhm) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad FWHM  value: %s\n", OptArg);
		            exit(-1);
	            	}
                    OptFilter  = True;
 		          break; 	
            case 'p':
 	     	      if (sscanf(OptArg,"%f",&ZeroPadding) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad zero padding  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
	           case 'l':
 	     	      if (sscanf(OptArg,"%d",&Lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
   	    case 'v': Verbose = True; break;
        case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
       

       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out1, argv[OptInd++]);
         else usage(argv);

    if (OptInd < argc) strcpy(Name_Imag_Out2, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
  
 
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
        cout << "# File Name Out Der Theta = " << Name_Imag_Out1 << endl;
        cout << "# File Name Out Der Phi = " << Name_Imag_Out2 << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
		if (OptFilter == True)  cout << "# Fwhm = " <<  Fwhm << endl;
    }

   Hdmap Map;
   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP");
   }

   int Nside = Map.Nside();
   if (Lmax <= 0) Lmax = 3 * Nside;
   int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;

   Hdmap mapdth;
   Hdmap mapdph;
   CAlmR  ALM;
   // ALM.Verbose = Verbose;
   ALM.Norm = False; 
   ALM.UseBeamEff = (ZeroPadding  > 0) ? True: False;
   ALM.alloc(Nside, Lmax_Tot);
   if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);

   // fits_write_fltarr("beameff.fits", ALM.BeamEff);
   
   ALM.alm_trans(Map);

   ALM.deriv(Map,  mapdth,   mapdph, Fwhm);
   // mrs_write_3maps("xxder", Map,  mapdth,   mapdph);
   mapdth.write(Name_Imag_Out1);
   mapdph.write(Name_Imag_Out2);
   exit(0);
}

