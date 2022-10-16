/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  J.L. Starck
**
**    Date:  98/02/03
**    
**    File:  im_fusion.cc
**
*******************************************************************************
**
**    DESCRIPTION  large image test 
**    ----------- 
**
**
**    PARAMETRES    
**    ----------    
** 
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "FFTN_2D.h"
#include "IM_Rot.h"

char Name_In[256];  /*  list of input image names  */
char Name_Out[256];  /*  list of input image names  */
 
Bool Verbose = False;
// int Nbr_Plan = 3;
// int Si=0;
// int Sj=0;
float Resol = 1.;
float Dens = 1.;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

/*************************************************************************/
 
/****************************************************************************/

static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input_image_list\n", argv[0]);

    fprintf(OUTMAN, "        [-r Resol]\n");
    fprintf(OUTMAN, "            Resolution for the Power Spectrum calculation.\n");
    fprintf(OUTMAN, "            Default is 1.\n"); 
    manline();

    fprintf(OUTMAN, "        [-d Density]\n");
    fprintf(OUTMAN, "            Density parameter.\n");
    fprintf(OUTMAN, "            Default is 1.\n"); 
    manline();
        
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose. \n");
    fprintf(OUTMAN, "              Default is no\n");     
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void imfusinit(int argc, char *argv[])
{
    int c;
 
    /* get options */
    while ((c = GetOpt(argc,argv,"vr:d:")) != -1) 
    {
	switch (c) 
        {
	    case 'd': 
	        if (sscanf(OptArg,"%f",&Dens) != 1) 
                {
                    fprintf(OUTMAN, "bad density parameter: %s\n", OptArg);
                    exit(-1);
                }
	        break;
            case 'r':
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&Resol) != 1) 
                {
                    fprintf(OUTMAN, "bad resolution paramter: %s\n", OptArg);
                    exit(-1);
                }
	        break;
	case 'v':
		/* verbose flag -v */
		Verbose = True;
		break;
	   case '?':
		usage(argv);
	}
    } 

    /* get optional input file names from trailing 
       parameters and open files */
    
    if (OptInd < argc)  strcpy(Name_In, argv[OptInd++]);
    else usage(argv);
    
   if (OptInd < argc)   strcpy(Name_Out, argv[OptInd++]);
    else usage(argv);

    
    /* make sure there are not too many parameters */ 
    if (OptInd < argc)
    {
	fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
    }
}

/*********************************************************************/

int main(int argc, char *argv[]) 
{
   char Cmd[512];  
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   imfusinit(argc, argv);
  
   if (Verbose == True)
   {       
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Name in = " << Name_In << endl;
      cout << "File Name Out = " << Name_Out << endl;  
   }
  
  Ifloat Data;
  io_read_ima_float(Name_In, Data);
  fltarray Spectrum;
  get_isotropic_spectrum(Data, Spectrum, Resol,Dens);
  io_1d_write_data(Name_Out, Spectrum);
  exit(0);	      
}

/****************************************************************************/
