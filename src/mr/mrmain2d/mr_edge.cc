/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12
**    
**    File:  mr_edge.cc
**
*******************************************************************************
**
**    DESCRIPTION  detect the multiscale edges in an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_edge.cc 1.0 98/01/12 CEA 1998 @(#)";
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Edge.h"
 
char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Write_Sup[256];          /* output support file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
type_transform Transform = TO_DIADIC_MALLAT; /* type of transform */
  
char Name_RMSMap[256]; 
Bool UseRMSMap=False;
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose=False;
 
/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_mr_file\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    nbr_scale_usage(Nbr_Plan);
    manline();
    vm_usage();
    manline();
    verbose_usage();
    manline();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    Bool NscaleOpt=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
     
    /* get options */
    while ((c = GetOpt(argc,argv,"n:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break; 
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		NscaleOpt = True;
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
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}


/****************************************************************************/

int main(int argc, char *argv[])
{
    int b,j,i,k;
    Ifloat Data;
     
     /* support image creation */
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    int Nl = Data.nl();
    int Nc = Data.nc();
 
    // noise model class initialization
    MultiResol MR_Data (Nl,Nc, Nbr_Plan, Transform, "MR_Transform");
    MR_Data.Border = I_MIRROR;
    MR_Data.transform(Data);
    
    // Compute the thresholded wavelet coefficents  
  
     float ValB1,ValB2,Mod,Phase;
     EDGE CED(ED_CANNY);
              
     for (b=0; b < MR_Data.nbr_band()-1; b += 2)
     {
          // compute the modulus and phase
          for (i=0; i< MR_Data.size_band_nl(b); i++) 
          for (j=0; j< MR_Data.size_band_nc(b); j++) 
          {
             ValB1 = MR_Data(b,i,j);
	     ValB2 = MR_Data(b+1,i,j);
             Mod = sqrt( ValB1*ValB1 + ValB2*ValB2);
 	     ARG (ValB1, ValB2, Phase);
	     // if (ValB1 < 0) Phase = PI - Phase;
	     MR_Data(b,i,j) = Mod;
	     MR_Data(b+1,i,j) = Phase;
	      
          }
          CED.edge_from_1deriv(MR_Data.band(b), MR_Data.band(b+1)); 
      }
      MR_Data.write(Name_Imag_Out);
      exit(0);
}
