/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  im_info.cc
**
*******************************************************************************
**
**    DESCRIPTION  print informations about an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/
 

#include <stdio.h>
#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Edge.h"
#include "IM_Math.h"
#include "DefFunc.h"

#include "fitsio2.h"

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);
char Name_Imag_In[100];

float BinHisto = 1.;
Bool Entrop = False;
Bool MSup = False;
Bool NewStat = False;
Bool Survival = False;

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    fprintf(OUTMAN, "         [-e]\n");
    fprintf(OUTMAN, "             Entropy and Fisher Information calculation.\n");
    fprintf(OUTMAN, "             Histogram bin size value is 1.\n");
    manline();
    fprintf(OUTMAN, "         [-E HistoBinSize]\n");
    fprintf(OUTMAN, "             Entropy and Fisher Information calculation.\n");
    manline();
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Skewness, Kurtosis, HC and HC^+  calculation.\n");
    manline();
    vm_usage();
    manline();
    verbose_usage();
    manline();
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
    Bool Verbose = False;

    /* get options */
    while ((c = GetOpt(argc,argv,"SMeE:vzZ:")) != -1) 
    {
        switch (c) 
        {
	   case 'S': Survival = True;  break;
	   case 'M': MSup = True;  NewStat = True; break;
           case 'E': Entrop = True;
  	             if (sscanf(OptArg,"%f",&BinHisto )  < 1) 
		     {
                       fprintf(OUTMAN, "Error: bad histogram bin size ...\n");
                       exit(-1);
                     } 
                     break;
           case 'e': Entrop = True;
                     break;
	   case 'v': Verbose = True;
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

void fisher_entrop_information(Ifloat &Ima, double & FishVal, 
                               double & EntropFried, float FluxIma)
{
    int i,j,Nl = Ima.nl();
    int Nc = Ima.nc();
    type_edge TypeEdge = ED_SEPPIXDIFF;
    EDGE CED(TypeEdge);
    float Flux = (FluxIma < 0) ? flux(Ima): FluxIma;
    Ifloat  Grad(Nl,Nc,"grad");
    double Pix;
    int P=0;

    if (Flux < FLOAT_EPSILON)
    {
        cout << "Error: total of image must be larger than zero ... " << endl;
        exit(-1);
    }
    EntropFried  = 0.;
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++) 
    {
       if (Ima(i,j) < 0) 
       {
          if (P ==0)
          {
            cout << "Warning: input image pixel contains negative values ... " << endl;
            cout << "         ==> negative pixels are set to zero ... " << endl;
            P = 1;
          }
          Ima(i,j) = 0;
       }
       else if (Ima(i,j) > 0)
       {
          Pix = Ima(i,j) / Flux;
          EntropFried += Pix*log(Pix);
          Ima(i,j) = sqrt(Pix);
       }
    }

    CED.Bord = I_MIRROR;
    CED.detection(Ima, Grad);
    // INFO(Grad, "grad");
    FishVal = 0;
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++) FishVal += Grad(i,j)*Grad(i,j);
    FishVal *=  4.;  
}

 
/*********************************************************************/


int main(int argc, char *argv[])
{

   lm_check(LIC_MR1);
   infinit(argc, argv);
     
/* temporarily add test for correct byteswapping.  Remove this after
   SADAM is more fully tested on all machines.
*/
   union u_tag {
     short ival;
     char cval[2];
   } u;

   u.ival = 1;
   if  ((BYTESWAPPED == True && u.cval[0] != 1) ||
        (BYTESWAPPED == False && u.cval[1] != 1) )
   {
   printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
   printf(" Byteswapping is not being done correctly on this system.\n");
   printf(" Check the MACHINE and BYTESWAPPED definitions in fitsio2.h\n");
   printf(" Please report this problem to the author at\n");
   printf("     pence@tetra.gsfc.nasa.gov\n");
   printf(  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
   }

/* end of temporary byteswap checking code */


    cout << "File name = " << Name_Imag_In   << endl << endl;
    if (io_detect_format(Name_Imag_In) == F_FITS)
    {
       fltarray Dat;
       double Mean,Sigma,Skew,Curt;
       float Min,Max;
       fits_read_fltarr(Name_Imag_In, Dat);
       moment4(Dat.buffer(),Dat.n_elem(), Mean,  Sigma, Skew, Curt, Min,Max);
      
       if (Dat.naxis() == 1)
          cout << "  Np = " << Dat.nx() << endl;
       else if (Dat.naxis() == 2)
         cout << "  Nl = " << Dat.ny()  << "  Nc = " << Dat.nx() << endl;
       else if (Dat.naxis() == 3)
         cout << "  Nx = " << Dat.nx() << 
	         "  Ny = " << Dat.ny() << 
		 "  Nz = " << Dat.nz() << endl;
       float FluxIma = Mean*Dat.n_elem();
       float Energy = (float) (Sigma*Sigma+Mean*Mean)*Dat.n_elem();
       cout << "  Min = " << Min << "  Max = " << Max << endl;
       cout << "  Mean = " << Mean << "  sigma = " << Sigma << endl;
       cout << "  Flux = " <<  FluxIma <<  "  Energy = " << Energy << endl;
       if (MSup == True)
       {
	 cout << "  Skew = " << Skew << "  Kurtosis = " << Curt << endl;
       }
       if (NewStat == True)
       {
         float HC1,HC2;
         hc_test(Dat.buffer(),Dat.n_elem(), HC1,HC2);
	 cout << "  HC1 = " << HC1 << "  HC2 = " << HC2 << endl;
	 //gausstest(Dat.buffer(),Dat.n_elem(),HC1,HC2);
	 //cout << "  NT1 = " << HC1 << "  NT2 = " << HC2 << endl;
       }
       
       // Attention fisher_entrop_information change data values !!!!
       if (Entrop == True)
       {
         cout << "   Shannon Entropy = " <<  entropy (Dat.buffer(), Dat.n_elem(), BinHisto) << endl;
         if (Dat.naxis() == 2)
         {
            Ifloat Frame;
            Frame.alloc (Dat.buffer(),  Dat.ny(),  Dat.nx());
            double FishInfo, FriedInfo;
            fisher_entrop_information (Frame,FishInfo, FriedInfo, FluxIma);
            cout << "  Fisher information = " << FishInfo << endl;
            cout << "  Frieden entropy = " << FriedInfo  << endl;
         }
       }
       
       if (Survival == True)
       {
          fltarray Survival;
          survival (Dat, Survival);
	  io_1d_write_data ("xx_survival", Survival);  
       }
      }
    else
    {
       Ifloat Dat;
       float FluxIma;
       double Mean,Sigma,Skew,Curt;
       float Min,Max;
       io_read_ima_float(Name_Imag_In, Dat);
       moment4(Dat.buffer(),Dat.n_elem(), Mean,  Sigma, Skew, Curt, Min,Max);
       FluxIma =Mean*Dat.n_elem();
       float Energy = (float) (Sigma*Sigma+Mean*Mean)*Dat.n_elem();
       cout << "  Nl = " << Dat.nl()  << "  Nc = " << Dat.nc() << endl;
       cout << "  Min = " << Min << "  Max = " << Max << endl;
       cout << "  Mean = " << Mean << "  Sigma = " << Sigma << endl;
       cout << "  Flux = " << FluxIma <<  "  Energy = " << Energy << endl;
       if (MSup == True)
  	 cout << "  Skew = " << Skew << "  Kurtosis = " << Curt << endl;
       if (NewStat == True)
       {
         float HC1,HC2;
         im_hc_test(Dat,HC1,HC2);
	 cout << "  HC1 = " << HC1 << "  HC2 = " << HC2 << endl;
	 //im_gaussianity_test(Dat,HC1,HC2);
	 //cout << "  NT1 = " << HC1 << "  NT2 = " << HC2 << endl;
       }       
       // Attention fisher_entrop_information change data values !!!!
       if (Entrop == True)
       {
            double FishInfo, FriedInfo;
            cout << "  Shannon Entropy = " <<  entropy (Dat.buffer(), Dat.n_elem(), BinHisto) << endl;
            fisher_entrop_information(Dat,FishInfo, FriedInfo, FluxIma);
            cout << "  Fisher information = " << FishInfo << endl;
            cout << "  Frieden entropy = " << FriedInfo  << endl;
       }
    }
 
    exit(0);
} 

