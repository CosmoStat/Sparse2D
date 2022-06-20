/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/00
**    
**    File:  cur_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet filtering program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Usage.h"
#include "PCur.h"
#include "FCur.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
Bool Simu=False;
int NbrDirection=DEF_CUR_BLOCK_SIZE;
  

Bool ExtendWT=False;
Bool IsotropWT=False;
Bool Reverse=False;
Bool Real=False; 
/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options input output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
 

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Use an isotropic 2D wavelet transform.\n");
    fprintf(OUTMAN, "         [-e]\n");
    fprintf(OUTMAN, "             Extended wavelet transform (redundancy of 7 instead of 4).\n");
    fprintf(OUTMAN, "         [-d NbrDirection]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrDirection);    
    manline();
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "             Real curvelet.\n");
    manline();
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Reconstruction.\n");
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
static void filtinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"Rd:rein:vzZ")) != -1) 
    {
	switch (c) 
        { 
	   case 'R':  Real = (Real== True) ? False: True; break;
	   case 'r':  Reverse = (Reverse== True) ? False: True; break;
	   case 'e':  ExtendWT = (ExtendWT== True) ? False: True; break;
	   case 'i':  IsotropWT = (IsotropWT== True) ? False: True; break;
 	   case 'd':
                if (sscanf(OptArg,"%d",&NbrDirection) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'v': Verbose = True;break;
	   case 'n':
                /* -n <NbrScale2D> */
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
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
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/*************************************************************************/


int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        if (NbrDirection > 0)  cout << "NbrDirection = " << NbrDirection <<endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
     }

     FCUR  Cur;
     Cur.Verbose = Verbose;   
      if (Reverse == True)
      {
         Ifloat Data;
         Cur.read(Name_Imag_In);
         Cur.cur_recons(Data);
         io_write_ima_float(Name_Imag_Out, Data);
      }
      else
      {
         Ifloat Data;
	 
         io_read_ima_float(Name_Imag_In, Data, &Header);
         Header.origin = Cmd;
         Cur.alloc_from_fine(NbrScale2D, Data.nl(), Data.nc(), NbrDirection, ExtendWT,IsotropWT, Real);
         Cur.cur_trans(Data);
	 if (Verbose == True) cout << "Write to " << Name_Imag_Out << endl;  
         Cur.write(Name_Imag_Out);
	 
// 	 if (Verbose == True) cout << "read " << endl;  
//  	 Ifloat Rec;
// 	 Rec=Data;
//  	 Cur.read(Name_Imag_Out);
// 	 if (Verbose == True) cout << "rec " << endl; 
//          Cur.cur_recons(Rec);
//  	 io_write_ima_float("tt.fits", Rec);
// 	 Data -= Rec;
// 	 INFO(Data, "residu");
      }
      exit(0);
}

