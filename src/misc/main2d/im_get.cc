/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/06/98
**    
**    File:  im_get.cc
**
*******************************************************************************
**
**    DESCRIPTION  extract a sub-image
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

 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
 
int I1=-1;
int I2=-1;
int J1=-1;
int J2=-1;
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options in_image out_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-x Xstart:XEnd]\n");
    fprintf(OUTMAN, "             First and last columns to extract.\n");
    manline();
    fprintf(OUTMAN, "         [-y Ystart:YEnd]\n");
    fprintf(OUTMAN, "             First and last lines to extract.\n");
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
static void infinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
   
    /* get options */
    while ((c = GetOpt(argc,argv,"x:y:vzZ:")) != -1) 
    {
	switch (c) 
        {
 	   case 'v': Verbose = True; break;         
            case 'x':
	        if (sscanf(OptArg,"%d:%d", &J1,&J2) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break; 
            case 'y':
                if (sscanf(OptArg,"%d:%d", &I1,&I2) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameters: %s\n", OptArg);
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

/*********************************************************************/


int main(int argc, char *argv[])
{
    Ifloat Dat1,Dat2;
    int i,j,Nl,Nc;
    fitsstruct Header;
  
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);

    io_read_ima_float(Name_Imag_In, Dat1, &Header);

    
    if (I1 < 0) I1 = 0;
    
    if ((I2 < 0) || (I2 >= Dat1.nl())) I2 = Dat1.nl() - 1;
      
    if (J1 < 0) J1 = 0;
    if ((J2 < 0) ||(J2 >= Dat1.nc())) J2 = Dat1.nc() - 1;
    Nl = I2 - I1 + 1;
    Nc = J2 - J1 + 1;
    Dat2.alloc(Nl,Nc,"get");
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "  Nl = " << Dat1.nl() << endl;
       cout << "  Nc = " << Dat1.nc() << endl;    
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Pos X = " << J1 << " Pos Y = " << J1 << endl;
       cout << "  Nl = " << Nl << endl;
       cout << "  Nc = " << Nc << endl;
       cout << endl;
    }       
    for (i=0;i < Nl; i++)
    for (j=0;j < Nc; j++) Dat2(i,j) = Dat1(i+I1,j+J1);
    
    io_write_ima_float (Name_Imag_Out, Dat2);
    exit(0);
} 

