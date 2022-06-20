/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: PS
**
**    Date:  97/08/22
**    
**    File:  im_test.cc
**
*******************************************************************************
**
**    DESCRIPTION  test
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
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    imopen_usage();
    manline();
    
    /*immorpho_usage(Elstr_Size);
    manline();*/
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
    while ((c = GetOpt(argc,argv,"vzZ:")) != -1) 
    {
	switch (c) 
        {
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
    Ifloat Dat,Dat_out;
    int i,j,k;
    float min, max;
    char Cmd[512];
    fitsstruct Header;

    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

    lm_check(LIC_MR1);
    infinit(argc, argv);

    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Header.origin = Cmd;
    
    if (Verbose == True )
    {
       cout << "Name File = " << Name_Imag_In << endl ;
       cout << "  Nl = " << Dat.nl() << endl;
       cout << "  Nc = " << Dat.nc() << endl << endl;
    }
    
    Dat_out.alloc(Dat.nl(),Dat.nc(), "sortie");
    min= Dat_out.min();
    max= Dat_out.max();
    cout << min << "  " << max << endl;
    cout << "  min= " << Dat_out.min() << endl;
    cout << "  max= " << Dat_out.max() << endl;
    cout << "  Nl = " << Dat.nl() << endl;
    cout << "  Nc = " << Dat.nc() << endl << endl;
    for (i=0; i < Dat.nl(); i++)
    {
    	for (j=0; j < Dat.nc(); j++)
    	{
       	Dat_out(i,j)= Dat(i,Dat.nc()-1-j);
       	}
    }

   io_write_ima_float (Name_Imag_Out, Dat_out, &Header);
    exit(0);
} 

