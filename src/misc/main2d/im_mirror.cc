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
**    Date:  20/12/98
**    
**    File:  im_mirror.cc
**
*******************************************************************************
**
**    DESCRIPTION  image mirror extension 
**    ----------- 
**                 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Deconv.h"

char Name_Imag_In[256];    // input file image 
char Name_Imag_Out[256];   // output file name 

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose=False;
Bool Neg=False;

int Nlz=0;
int Ncz=0;

/*********************************************************************/


static void usage(char *argv[])
{
    
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
         
    fprintf(OUTMAN, "         [-x Number_of_Columns]\n");
    fprintf(OUTMAN, "             Number of columns in the output image.\n");
    fprintf(OUTMAN, "             By default, the next power of 2.\n");
    manline();
    fprintf(OUTMAN, "         [-y Number_of_Lines]\n");
    fprintf(OUTMAN, "             Number of lines in the output image.\n");
    fprintf(OUTMAN, "             By default, the next power of 2.\n");
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
static void decinit(int argc, char *argv[])
{
    int c;    
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif  
   
    /* get options */
    while ((c = GetOpt(argc,argv,"nx:vy:zZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'n': Neg = True; break;
	   case 'v': Verbose = True; break;
	   case 'x':
		/* -d <type> type of deconvolution */
		if (sscanf(OptArg,"%d",&Ncz ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of columns: %s\n", OptArg);
	            exit(-1);
 		}
		break;
 	  case 'y':
		/* -d <type> type of deconvolution */
		if (sscanf(OptArg,"%d",&Nlz ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of lines: %s\n", OptArg);
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
		fprintf(OUTMAN, "Error: Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif	    
}

/**************************************************************************/

static int test_ind (int ind, int N)
{
    int Val;

    if (ind < 0) Val = - ind;
    else
    {
        if (ind >= N) Val = 2 * (N - 1) - ind;
        else Val = ind;
    }
    if ((Val >= N) || (Val < 0)) Val = -1;
    return (Val);
}

/****************************************************************************/

void im_neg_mirror (const Ifloat &Imag, Ifloat &Imag_Out)
{
    int Nl0 = Imag.nl();
    int Nc0 = Imag.nc();
    int Nl1 = Imag_Out.nl();
    int Nc1 = Imag_Out.nc();
    int i1,j1,i0,j0, Depi, Depj;

    Depi = (Nl1 - Nl0) / 2;
    Depj = (Nc1 - Nc0) / 2;

    for (i1 = 0; i1 < Nl1; i1++)
    for (j1 = 0; j1 < Nc1; j1++)
    {
       i0 = i1 - Depi;
       j0 = j1 - Depj;
       if ((i0 < 0) || (j0 < 0) || (i0 >= Nl0) || (j0  >= Nc0)) 
       {
          int ii0 = test_ind (i1 - Depi, Nl0);
          int ij0 = test_ind (j1 - Depj, Nc0);
          if ((i0 < 0) || (j0 < 0)) Imag_Out (i1,j1) = 0.;
          else Imag_Out (i1,j1) = - Imag (ii0,ij0);
       }
       else Imag_Out (i1,j1) = Imag(i0,j0);
    }
}

/**************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    Ifloat Imag, Imag_Out;
    char Cmd[256];

    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR1);
    
     /* Get command line arguments, open input file(s) if necessary */
    decinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
     }

    /* read input image */
    io_read_ima_float(Name_Imag_In, Imag);
    
    if ((Nlz == 0) && (Ncz == 0))
    {
       dec_line_column (Imag.nl(), Nlz);
       dec_line_column (Imag.nc(), Ncz);
       if (Nlz < Ncz) Nlz = Ncz;
    }
    else 
    {
       if (Nlz == 0) Nlz = Imag.nl();
       if (Ncz == 0) Ncz = Imag.nc();
    }
     
    if ((Nlz < Imag.nl()) || (Ncz < Imag.nc()))
    {
       cerr << "Error: output image dimensions must be greater than input image dimensions ... " << endl;
       exit(-1);
    }
    Imag_Out.alloc(Nlz, Ncz, "zoom");
    if (Neg == True) im_neg_mirror ( Imag,  Imag_Out);
    im_extend (Imag, Imag_Out);
    io_write_ima_float(Name_Imag_Out, Imag_Out);
    exit(0);
}

