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
**    File:  mrs_getcmb.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Get CMB realizations using a given cosmological model    
**    -----------
**
******************************************************************************/

#include"HealpixClass.h"
#include"omp.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool NormALM = False;
float  ZeroPadding=0.;
int Lmax = 10.;
 
int Nside=256;
int NbrCMBMap=1;

float Fwhm_arcmin=5.;
Bool Old = False;
int NsideAngle=32;
Bool ALMin=False;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "             Input file is an Alm file, and not a map.\n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is %d.\n", Lmax);
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");
    
    fprintf(OUTMAN, "         [-n Nside]\n");
    fprintf(OUTMAN, "             Nside resolution parameter. Default is %d.\n", NsideAngle);
    
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "oaI:N:n:C:f:l:p:vzZ")) != -1) 
    {
	switch (c) 
        { 
            case 'a': ALMin = True; break;
            case 'n':
                if (sscanf(OptArg,"%d",&NsideAngle) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
                }      
                break;
             case 'o': Old = True; break;
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
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (NsideAngle > 0)  cout << "# NsideAngle = " <<  NsideAngle << endl;
    }
   omp_set_num_threads(1);
   
   int Nside=0, Lmax_Tot;
   CAlmR  ALM;
   if (ALMin == False)
   {
      Hdmap Map;
      Map.read(Name_Imag_In);
      Nside = Map.Nside();
      if (Lmax <= 0) Lmax = 3 * Nside;
      Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
       // ALM.Verbose = Verbose;
       ALM.Norm = False; 
       ALM.UseBeamEff = (ZeroPadding  > 0) ? True: False;
       ALM.alloc(Nside, Lmax_Tot);
       if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);
       ALM.alm_trans(Map);       
       // if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
       if (Verbose == True) cout << "Read map file " << Name_Imag_In << ", Nside = " <<  Nside << endl;
   } 
   else 
   {
       Nside = NsideAngle;
       ALM.read(Name_Imag_In, NsideAngle);
       Lmax_Tot = ALM.Lmax();
       if (Verbose == True) cout << "Read Alm file " << Name_Imag_In << ", Lmax = " <<  ALM.Lmax() << endl;
   }
    

    
    CAlmR  RotALM;
    RotALM.Norm = ALM.Norm; 
    RotALM.UseBeamEff = ALM.UseBeamEff;
    RotALM.alloc(Nside, Lmax_Tot);
   
   NbrCMBMap = NsideAngle * NsideAngle * 12;
   pointing Angle;
   fltarray TabThetaPhi;
   TabThetaPhi.alloc(NbrCMBMap, 5);
   fltarray Are, Aim;
   Are.alloc(Lmax_Tot+1, Lmax_Tot+1, NbrCMBMap);
   Aim.alloc(Lmax_Tot+1, Lmax_Tot+1, NbrCMBMap);
  
   Hdmap MapResol;
   MapResol.alloc (NsideAngle);
   ; MapResol.info();
   vec3 Vect;
   for (int i=0; i < NbrCMBMap; i++)
   {
      double psi=1.;
      Angle = MapResol.pix2ang(i);
      Vect =  MapResol.pix2vec(i);
      TabThetaPhi(i,0) = Angle.theta;
      TabThetaPhi(i,1) = Angle.phi;
      TabThetaPhi(i,2) = Vect.x;  
      TabThetaPhi(i,3) = Vect.y;
      TabThetaPhi(i,4) = Vect.z;

      for (int l=0; l <= ALM.Lmax(); l++)
      for (int m=0; m <= l; m++) 
      {
         RotALM(l,m) = ALM(l,m);
       }
      
      if (Old == True) rotate_alm(RotALM, psi, (double) Angle.theta, (double)  Angle.phi);
      else rotate_alm(RotALM, (double) Angle.phi, (double) Angle.theta, (double)  psi);
    
      for (int l=0; l <= ALM.Lmax(); l++)
      for (int m=0; m <= l; m++) 
      {
               Are(l,m,i) = real(RotALM(l,m));
               Aim(l,m,i) = imag(RotALM(l,m));
               // if ((l == 192) && (m == 10)) cout <<  A(l,m,0) << " " <<  A(l,m,1) << endl;
      }
   }
   char FN[512]; /* input file image */
   sprintf(FN, "%s_almre.fits", Name_Imag_Out);
   fits_write_fltarr(FN, Are);
   sprintf(FN, "%s_almim.fits", Name_Imag_Out);
   fits_write_fltarr(FN, Aim);
   sprintf(FN, "%s_thetaphi.fits", Name_Imag_Out);
   fits_write_fltarr(FN, TabThetaPhi);
   exit(0);
}

