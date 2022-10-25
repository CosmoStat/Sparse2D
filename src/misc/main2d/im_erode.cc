/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  97/08/22
**    
**    File:  im_erode.cc
**
*******************************************************************************
**
**    DESCRIPTION  image segementation
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

int Order = 1;
 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
int Elstr_Size = 3;
int Elstr_Shape = 1;
int TypeElem = 3;
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    imerode_usage();

    immorpho_usage(Elstr_Size);
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
    while ((c = GetOpt(argc,argv,"n:d:s:vzZ:")) != -1) 
    {
	switch (c) 
        {
   	   case 'v': Verbose = True; break;       
            case 's':
                if (sscanf(OptArg,"%d",&TypeElem) != 1) 
                {
		    fprintf(OUTMAN, "bad structural element type: %s\n", OptArg);
	            usage(argv);
 		}
                if ((Order <= 0) && (Order > 3))   
                {
		    fprintf(OUTMAN, "bad structural element type: %s\n", OptArg);
	            usage(argv);
 		}
 		break; 
            case 'n':
                if (sscanf(OptArg,"%d",&Order) != 1) 
                {
		    fprintf(OUTMAN, "bad erosion number: %s\n", OptArg);
	            usage(argv);
 		}
                if ((Order <= 0) && (Order > 50))   
                {
		    fprintf(OUTMAN, "bad erosion number: %s\n", OptArg);
	            usage(argv);
 		}
 		break;           
 	     case 'd':
		/* -d <Elstr_Size> */
		if (sscanf(OptArg,"%d",&Elstr_Size) != 1) 
                {
		    fprintf(OUTMAN, "bad Elstr_Size: %s\n", OptArg);
		    usage(argv);
		}
                if ((Elstr_Size < 3.) || (Elstr_Size > 100.)) 
                                        Elstr_Size = 5;
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
    Ifloat Dat,Erode;
    int k,i;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
 
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);
    
    io_read_ima_float(Name_Imag_In, Dat, &Header);
     Header.origin = Cmd;
   
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Nl = " << Dat.nl() << endl;
       cout << "  Nc = " << Dat.nc() << endl << endl;    
       cout << "  Erosion Number = " << Order << endl << endl;
       if (TypeElem == 1) 
           cout << "  Structural element = square" << endl;
       else if (TypeElem == 2)
            cout << "  Structural element =  cross" << endl;
       else cout << "  Structural element =  circle" << endl;
       cout << endl;
    }
     
    Erode.alloc(Dat.nl(),Dat.nc(), "Erosion");
    for (i=0; i < Order; i++)
    {
       switch(TypeElem)
       {
       case 1: morpho_erosion (Dat, Erode, Elstr_Size);break;
       case 2: morpho4_erosion (Dat, Erode); break;
       case 3: morpho_cercle_erosion (Dat, Erode, Elstr_Size); break;
       default: break;
       }
       if (i != Order - 1) Dat = Erode;
    }
    
    io_write_ima_float (Name_Imag_Out, Erode, &Header);
    exit(0);
} 

