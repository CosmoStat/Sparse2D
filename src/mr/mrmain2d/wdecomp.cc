/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  16/03/2000 
**    
**    File:  wdecomp.cc
**
*******************************************************************************
**
**    DECRIPTION  decompression program
**    ---------- 
**
******************************************************************************/
 
   
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Noise.h"
#include "IM_Comp.h"
#include "MR_Noise.h"
#include "MR_Comp.h"


int Resol=-1;
char File_Name_Imag[80], File_Name_Transform[80];
Bool Verbose=False;
 
extern void initfield(fitsstruct *Header); /* initialize the structure FITSSTRUCT */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

extern Bool InputFromStdin;
char OutputType='u';
Bool SimuNoise = False;
int IterRec=0;
int NumBlock=-1;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options CompressedFileName OutputFileName \n\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    fprintf(OUTMAN, "        [-B BlockNbr]\n");
    fprintf(OUTMAN, "              Decompress only one block. \n");
    fprintf(OUTMAN, "              BlockNbr is the block number to decompress.\n");
    fprintf(OUTMAN, "              Default is no.\n");
    manline();    
  
    fprintf(OUTMAN, "        [-r resolution]\n");
    fprintf(OUTMAN, "          resol = 0..nbr_of_scale-1 \n");
    fprintf(OUTMAN, "          resol = 0 for full resolution (default) \n");
    fprintf(OUTMAN, "          resol = nbr_of_scale-1 for the worse resol.\n");
    manline();    

//     fprintf(OUTMAN, "        [-t] output type\n");
//     fprintf(OUTMAN, "              if the input image was a fits image, \n");
//     fprintf(OUTMAN, "              the image output type can be fixed by the user \n");
//     fprintf(OUTMAN, "              to 'f' for float, to 'i' for integer, or 's' \n");
//     fprintf(OUTMAN, "              for short. By default, the output type \n");
//     fprintf(OUTMAN, "              is the same as the type of the original image. \n");
//     manline();

    vm_usage();
    manline();   
    verbose_usage();
    manline();
    manline();
    manline(); 
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void hcinit(int argc, char *argv[])
{
    int c, L;
    char *Ptr;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif

    strcpy(File_Name_Transform, "stdin");
    InputFromStdin = True;

    /* get options */
    while ((c = GetOpt(argc,argv,"vr:B:zZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v':
		/* Verbose flag -v */
		Verbose = True;
		break;
            case 'B': 
                if (sscanf(OptArg,"%d",&NumBlock) != 1)
                {
		   fprintf(OUTMAN, "Error: bad block parameter: %s\n", OptArg);
		   exit(-1);
		}
		break;
             case 'r': /* -r < resolution> */
		if (sscanf(OptArg,"%d",&Resol) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad Resolution parameter: %s\n", OptArg);
		    exit(-1);
		}
                if ((Resol < 0) || (Resol > 20))
                {
		    fprintf(OUTMAN, "Error: bad Resolution parameter: %s\n", OptArg);
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

    /*  MRC file name */
    if (OptInd < argc) strcpy(File_Name_Transform, argv[OptInd++]);
    else usage(argv);
    L = strlen (File_Name_Transform);
    if ((L == 1) && (File_Name_Transform[0] == '-')) InputFromStdin = True;
    else 
    {
       InputFromStdin = False;
       Ptr = File_Name_Transform;
       if (L < 5) strcat (File_Name_Transform, ".MRC");
       else if ((Ptr[L-1] != 'C') || (Ptr[L-2] != 'R') || (Ptr[L-3] != 'M')) 
                         strcat (File_Name_Transform, ".MRC");
    }		
	
    /* output file name */		
    if (OptInd < argc) strcpy(File_Name_Imag, argv[OptInd++]);
    else usage(argv);

    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
	fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
    }
    
//     if (Opt_s_or_o == False)
//     {
// 	fprintf(OUTMAN, "Error: -s or -o option  has to be set ...\n");
// 	usage(argv);
//     }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***************************************************/

int main (int argc, char *argv[])
{ 
    int k;
    char Cmd[512];
	
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments*/

    // lm_check(LIC_MR3); 
    hcinit(argc, argv);

  MR_CompData DecIma;
    DecIma.Cmd = Cmd;
    DecIma.File_Name_Imag=File_Name_Imag;
    DecIma.File_Name_Transform=File_Name_Transform;
    DecIma.Resol=Resol;
    DecIma.OutputType=OutputType;
    DecIma.UseLiftingInsteadOfFilter = True;
    if (SimuNoise == True) DecIma.SimuNoise=True;
    else DecIma.SimuNoise=False;
    if (Verbose==True) DecIma.Verbose = True;

    // apply the decompression
    if (NumBlock > 0) DecIma.decompress(NumBlock);
    else DecIma.decompress();
    
    // save the result (field Dat) in DecIma.File_Name_Imag
    // in case where the full image is decompressed with block compression
    // options, the result is already in the output file.
    if ((NumBlock > 0) || (DecIma.NbrBlock <= 1)) DecIma.write();
    exit(0);
}



 


