/******************************************************************************
**                   Copyright (C) 2006 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  im_fdr.cc
** 
**
******************************************************************************/
 

#include <stdio.h>
#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Math.h"
#include "NR.h"
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);
char Name_Imag_In[100];
char Name_Imag_Out[100];
Bool Correl = False;
float Alpha=0.05;
Bool OutFile = False;
Bool Verbose = False;
float Noise_Ima=0.;

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Image [FileOut] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    fprintf(OUTMAN, "         [-c]\n");
    fprintf(OUTMAN, "             Data are correlated.\n");
    fprintf(OUTMAN, "             Default is no.\n");
    fprintf(OUTMAN, "         [-g Sigma]\n");
    fprintf(OUTMAN, "             Gaussian noise standard deviation.\n");    
    fprintf(OUTMAN, "         [-a Alpha Value]\n");
    fprintf(OUTMAN, "             Default is %f.\n", Alpha);
    vm_usage();
    verbose_usage();
    manline();
    exit(-1);
}

/*********************************************************************/

static void infinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
    //Bool Verbose = False;

    /* get options */
    while ((c = GetOpt(argc,argv,"g:ca:vzZ:")) != -1) 
    {
        switch (c) 
        {
            case 'g': /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
		{
                    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
                   exit(-1);
                }
               break;
	   case 'v': Verbose = True; break;
	   case 'c': Correl = True;  break;
           case 'a': if (sscanf(OptArg,"%f",&Alpha )  < 1) 
		     {
                       fprintf(OUTMAN, "Error: bad alpha value ...\n");
                       exit(-1);
                     } 
                     break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
            case '?': usage(argv); break;
	    default: usage(argv); break;
       }
    }
    if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
    else usage(argv);

    if (OptInd < argc)
    {
       strcpy(Name_Imag_Out, argv[OptInd++]);
       OutFile=True;
       if (OptInd < argc) usage(argv);
    }
    
    
    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
        fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
        usage(argv);
    }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

double get_fdr_pvalue(double *TabPVal, int N, double Alpha, Bool Correl)
{
   int i;
   dblarray PVal(N);
   double p_cutoff=0;
   
   for (i = 0; i < N; i++) PVal(i) = TabPVal[i];
   double *P_val_l_tri = PVal.buffer()-1;
   sort(N,P_val_l_tri);
   
   for(i = 0; i < N; i++)
   {
      double  j_alpha = Alpha * (double) i / (double) N;
      if (PVal(i) <= j_alpha) p_cutoff = PVal(i);
   }
   //cout << "p_cutoff =" << p_cutoff << endl;
   for(i = 0; i < N; i++) TabPVal[i] = (  TabPVal[i] < p_cutoff ) ? 1 : 0;
   return p_cutoff;
}


/*********************************************************************/

int main(int argc, char *argv[])
{
   lm_check(LIC_MR1);
   infinit(argc, argv); 
   double PDet=0.;
   double cn = 0.;
   if (Verbose == True)
   { 
     cout << endl << endl << "PARAMETERS: " << endl << endl;
     cout << "File Name in = " << Name_Imag_In << endl;
     cout << "Alpha = " <<    Alpha << endl;  
     if (Correl == True) cout << "Correlated data " << endl; 
     cout << endl; 
   }  
    if (io_detect_format(Name_Imag_In) == F_FITS)
    {
       //dblarray Dat;
       //fits_read_dblarr(Name_Imag_In, Dat);
       dblarray Dat;
       fits_read_dblarr(Name_Imag_In, Dat);
       if (Noise_Ima > 0)
       {
          for (int i=0; i < Dat.n_elem() ; i++)   Dat(i) = 1.-erf(abs( (double) Dat(i))/(sqrt(2.)* (double) Noise_Ima));
       }
        if (Correl == True)
         { 
	 for (int j=1; j <= Dat.n_elem() ; j++)  cn = cn + 1./((double) j); 
 	 }
       else cn = 1;
       if (Verbose == True) cout << "cn = " << cn << endl; 
       PDet=get_fdr_pvalue(Dat.buffer(), Dat.n_elem(), Alpha/cn, Correl);
       if (Verbose == True) printf("FDR   PDet = %1.10f\n",   PDet);
       if (Verbose == True)  cout << "Nbr detect = " << Dat.total()  << endl;
       
       if (OutFile == True) 
       {
        //dblarray R(1);
       //R(0) = PDet;
       //fits_write_dblarr(Name_Imag_Out, R);
       fits_write_dblarr(Name_Imag_Out, Dat);
       }
       
    }    
    else
    { 
       Ifloat Dat;
       io_read_ima_float(Name_Imag_In, Dat);
       //dblarray Dat1(Dat.n_elem());
       dblarray Dat1(Dat.n_elem());
       
       if (Noise_Ima > 0)
       {
          for (int i=0; i < Dat.n_elem() ; i++)   Dat1(i) = 1.-erf(abs( (double) Dat(i))/(sqrt(2.)* (double) Noise_Ima));
       }
       else  for (int i=0; i < Dat.n_elem() ; i++) Dat1(i) = Dat(i);
       
       if (Correl == True)
         { 
 	 for (int j=1; j <= Dat.n_elem() ; j++)  cn = cn + 1./((double) j); 
 	 }
       else cn = 1;
       if (Verbose == True)     cout << "cn = " << cn << endl;
        PDet=get_fdr_pvalue(Dat1.buffer(), Dat1.n_elem(), Alpha/cn, Correl);
	for (int i=0; i < Dat.n_elem() ; i++)   Dat(i) = (float) Dat1(i);
    
        if (Verbose == True) printf("FDR   PDet = %1.10f\n",   PDet);
        if (Verbose == True)  cout << "Nbr detect = " << Dat.total()  << endl;
        if (OutFile == True) 
        {
          //dblarray R(1);
          //R(0) = PDet;
          //fits_write_dblarr(Name_Imag_Out, R);
          io_write_ima_float(Name_Imag_Out, Dat);
	}
    }
    exit(0);
} 

