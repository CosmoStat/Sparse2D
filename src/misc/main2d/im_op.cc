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
char Name_Oper[100];
char Name_Imag_In2[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options in_im1 [+-*/] in_im2  out_image \n\n", argv[0]);
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
    if (OptInd < argc) strcpy(Name_Oper, argv[OptInd++]);
    else usage(argv);

    if (OptInd < argc) strcpy(Name_Imag_In2, argv[OptInd++]);
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
    fitsstruct Header;
    int k;
    char Cmd[512];
	
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
      
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);

    if (strlen(Name_Oper) != 1)
    {
       cerr << "Error: bad operator parameter ... " << endl;
       cerr << "       operators are +,-,*,/ " << endl;
       exit(-1);
    }
    io_read_ima_float(Name_Imag_In, Dat1, &Header);
    io_read_ima_float(Name_Imag_In2, Dat2);
    
    if ((Dat1.nl() != Dat2.nl()) || (Dat1.nc() != Dat2.nc()))
    {
       cerr << "Error: images must have the same size ... " << endl;
       exit(-1);
    }
    Header.origin = Cmd;
   
   switch(Name_Oper[0])
   {
      case '+': Dat1 += Dat2; break;
      case '-': Dat1 -= Dat2; break;
      case '*': Dat1 *= Dat2; break;
      case '/': Dat1 = Dat1 / Dat2; break;
      default:
         cerr << "Error: bad operator parameter ... " << endl;
         cerr << "       operators are +,-,*,/ " << endl;
         exit(-1);
      break;
   }    
   if (Verbose == True )
   {
      cout << "Name File in 1 = " << Name_Imag_In << endl ;
      cout << "Name File in 2 = " << Name_Imag_In2 << endl ;
      cout << "Name File out = " << Name_Imag_Out << endl ;
      cout << "Operation: " << Name_Oper[0] << endl;
   }
    
    io_write_ima_float (Name_Imag_Out, Dat1, &Header);
    exit(0);
} 

