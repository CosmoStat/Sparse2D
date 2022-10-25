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
**    File:  im_put.cc
**
*******************************************************************************
**
**    DESCRIPTION  insert a sub-image
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
int J1=-1;
Bool Verbose = False;
 
/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options in_image inout_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-x Xstart]\n");
    fprintf(OUTMAN, "             x-coordinate of the image to insert.\n");
    manline();
    fprintf(OUTMAN, "         [-y Ystart]\n");
    fprintf(OUTMAN, "             y-coordinate of the image to insert.\n");
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
	        if (sscanf(OptArg,"%d", &J1) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameter: %s\n", OptArg);
		    exit(-1);
		}
 		break; 
            case 'y':
                if (sscanf(OptArg,"%d", &I1) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameter: %s\n", OptArg);
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
    int k;
    char Cmd[512];
	
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);

    io_read_ima_float(Name_Imag_In, Dat1);
    io_read_ima_float(Name_Imag_Out, Dat2, &Header);
    Header.origin = Cmd;

    if (I1 < 0) I1 = 0;  
    if (J1 < 0) J1 = 0;
    if (I1 >= Dat2.nl()) I1 = Dat2.nl() - 1;
    if (J1 >= Dat2.nc()) J1 = Dat2.nc() - 1;
    Nl = Dat1.nl();
    Nc = Dat1.nc();
    
    if ((I1+Nl > Dat2.nl()) || (J1+Nc > Dat2.nc()))  
    {
       cerr << "Error: first image cannot be inserted in the second one ... " << endl;
       cerr << "       First image: Nl = "  << Nl << " Nc = " << Nc << endl;
       cerr << "       Second image: Nl = " << Dat2.nl()  << " Nc = " << Dat2.nc()  << endl;
       cerr << "       Insertion position : y = " << J1  << "  x = " << I1  << endl;
       exit(-1);
    }
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "  Nl = " << Dat1.nl() << endl;
       cout << "  Nc = " << Dat1.nc() << endl;    
       cout << "Name File in-out 2  = " << Name_Imag_Out << endl ;
       cout << "  Pos X = " << J1 << " Pos Y = " << J1 << endl;
       cout << "  Nl = " << Nl << endl;
       cout << "  Nc = " << Nc << endl;
       cout << endl;
    }              
    for (i=0;i < Nl; i++)
    for (j=0;j < Nc; j++) Dat2(i+I1,j+J1) = Dat1(i,j);
    
    io_write_ima_float (Name_Imag_Out, Dat2,&Header);
    exit(0);
} 

