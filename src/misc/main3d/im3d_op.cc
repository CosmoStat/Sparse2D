/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe
**
**    Date:  20/02/01
**    
**    File:  im3d_op.cc
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
#include "IM3D_IO.h"

 
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
    fprintf(OUTMAN, "\n\nUsage: %s options in_cub1 [+-*/] in_cub2  out_cube \n\n", argv[0]);
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


int main(int argc, char *argv[]) {

   fitsstruct Header;
   int k;
   char Cmd[512];
   extern softinfo Soft;

   Soft.mr3();
	
   Cmd[0] = '\0';
   for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
      
   /* Get command line arguments, open input file(s) if necessary */
   lm_check(LIC_MR3);
   infinit(argc, argv);

   if (strlen(Name_Oper) != 1) {
      cerr << "Error: bad operator parameter ... " << endl;
      cerr << "       operators are +,-,*,/ " << endl;
      exit(-1);
   }
    
   fltarray Dat1, Dat2;
   io_3d_read_data(Name_Imag_In, Dat1, &Header);
   io_3d_read_data(Name_Imag_In2, Dat2);
    
   if (     ((Dat1.nx() ==0) || (Dat1.nx() != Dat2.nx())) 
         || ((Dat1.ny() ==0) || (Dat1.ny() != Dat2.ny()))
         || ((Dat1.nz() ==0) || (Dat1.nz() != Dat2.nz()))) {
      cerr << "Error: cubes must have the same size ... " << endl;
      exit(-1);
   }
   Header.origin = Cmd;
   
   switch(Name_Oper[0]) {
      
      case '+': Dat1 += Dat2; break;
      case '-': Dat1 -= Dat2; break;
      case '*': Dat1 *= Dat2; break;
      case '/': Dat1 /= Dat2; break;
      default:
         cerr << "Error: bad operator parameter ... " << endl;
         cerr << "       operators are +,-,*,/ " << endl;
         exit(-1);
      break;
   }    
   
   if (Verbose == True ) {
      cout << "Name File in 1 = " << Name_Imag_In << endl ;
      cout << "Name File in 2 = " << Name_Imag_In2 << endl ;
      cout << "Name File out = " << Name_Imag_Out << endl ;
      cout << "Operation: " << Name_Oper[0] << endl;
   }
    
   io_3d_write_data (Name_Imag_Out, Dat1, &Header);
   exit(0);
} 

