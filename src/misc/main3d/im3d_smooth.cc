/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  23/11/99
**    
**    File:  im3d_smooth.cc
**
*******************************************************************************
**
**    DESCRIPTION  Smmoth a cube by a Gaussian
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
 

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "FFTN_3D.h"
#include "DefFunc.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
float Sigma=3.;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options cube_in cube_out\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-c Std]\n");
    fprintf(OUTMAN, "             Convolve the input data with a Gaussian width sigma=Std.\n");
    fprintf(OUTMAN, "             Default is %f\n", Sigma);
    manline();        
    
    vm_usage();
    manline();
    verbose_usage();
    fprintf(OUTMAN, "               \n\n");
    exit(0);
}

/****************************************************************************/

static void init(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
 
    /* get options */
    while ((c = GetOpt(argc,argv,"c:vzZ:")) != -1) 
    {
        switch (c) 
        {
            case 'v': Verbose = True;break;
	    case 'c': 
                if (sscanf(OptArg,"%f",&Sigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Std parameter: %s\n", OptArg);
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
	usage(argv);
    }
    
    #ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    fitsstruct Header;
    fltarray Dat,Gauss;
    int i,j,k;
    Bool NewIma = False;
    char Cmd[512];
    extern softinfo Soft;

    Soft.mr3();
	    
    lm_check(LIC_MR3);
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    init(argc, argv);

    io_3d_read_data(Name_Imag_In, Dat, &Header);
    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();
    
    if (Verbose == True )
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << " Nx = " << Dat.nx();
       cout << " Ny = " << Dat.ny();
       cout << " Nz = " << Dat.nz()  << endl;
       cout << " Sigma = " << Sigma  << endl;
    }
    FFTN_3D FFT;
    Gauss.alloc(Nx,Ny,Nz);
    make_gaussian3d(Gauss,Sigma);
    FFT.convolve(Dat, Gauss);
    io_3d_write_data (Name_Imag_Out, Dat, &Header);
    exit(0);
}

