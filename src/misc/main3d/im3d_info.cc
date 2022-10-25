/******************************************************************************
**                   Copyright (C) 1999 by CEA
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
**    File:  im3d_info.cc
**
*******************************************************************************
**
**    DESCRIPTION  print information about a multidimensional image
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

Bool NormByt = False;
int MaxVal=1;
int MinVal=0;

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
Bool InfoFrame=False;
 
/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options input_file \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "             Print also information about all individual frames.\n");
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
    while ((c = GetOpt(argc,argv,"avzZ:")) != -1) 
    {
        switch (c) 
        {
	  case 'a': InfoFrame=True;break;
  	  case 'v': Verbose = True; break;
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
    fltarray Dat;
    int k;
    Ifloat Frame;
    
    lm_check(LIC_MR3);
    init(argc, argv);
 
    io_3d_read_data(Name_Imag_In, Dat);
    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();

    cout << "Name File in = " << Name_Imag_In << endl ;
    cout << "   Nx = " << Dat.nx();
    cout << " Ny = " << Dat.ny();
    cout << " Nz = " << Dat.nz()  << endl;    
    cout << "   Min = " << Dat.min() << " Max = " << Dat.max() << endl;
    cout << "   Mean = " << Dat.mean() << "  Sigma = " << Dat.sigma() << endl;
    cout << "   Flux = " << Dat.total() <<  " Energy = " << (Dat*Dat).total() << endl << endl;
    // float *Ptr = (float *) Dat.buffer();
    // printf("P(0) = %f\n", (float) (*Ptr));

    if (InfoFrame == True)
    {
       float *Ptr= Dat.buffer();
       for (k=0; k < Nz; k++)
       {
          //cout << "P(0)" << Ptr[0] << endl;
          Frame.alloc(Ptr, Ny, Nx);
          cout << "Frame " << k+1 << "  Nl = " << Frame.nl() << " Nc = " << Frame.nc() << endl;
          cout << "   Min = " << Frame.min() << " Max = " << Frame.max() << endl;
          cout << "   Mean = " << Frame.mean() << " Sigma = " <<  Frame.sigma() << endl;
          cout << "   Flux = " << Frame.total() <<  " Energy = " << energy (Frame) << endl;
          //cout << "P(0)" << Ptr[0] << endl;
	  Ptr += Nx*Ny;
       }
    }
  
    exit(0);
} 

