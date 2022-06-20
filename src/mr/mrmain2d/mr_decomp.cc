/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  mr_decomp.cc
**
*******************************************************************************
**
**    DECRIPTION  decompression program
**    ---------- 
**
*******************************************************************************
**
**    PARAMETRES
**    ----------
**
**        Input_File: image file name to compress
**        [Output_File]: compressed image file
**                       by default, "Input_File.P"
**           [-v] : Verbose
**
**		"Usage: %s [-v] [-r resolution] input output\n",
**        [-v] : Verbose
**        [-r resolution] : we are always interested to decompress
**                          the image at full resolution. This 
**                          parameter allows to extract from the
**                          compressed file, the image at worse 
**                          resolution, and produces by this way 
**                          a smaller image than the original
**                          Even is the compressed file contained
**                          the noise of the input image, the noise
**                          will not be decompressed.
**                          resolution must be >= 0 and < number of scales
**                          of the transform
**                          be default, the image is reconstructed at
**                          full resolution with its noise if it exists
**                          in the compressed file.
**         input : compressed file name
**        output : decompressed file name
**      
**        [-s] : the decompressed image is send to the 
**                standard output (screen). 
**
**        [-i] : if this option is not set, input 
**                data are read from the standard input. 
**
**        [-o] : if this option is set, the decompressed
**                image is written to a file. If not, the
**                -s option has to be set\n");
**
**         [-t] output type\n");
**                if the input image was a fits image, 
**                the image output type can be fixed to 
**                by the user to 'i' for integer, or 's' 
**                for short. By default, the output type 
**                is the same as the type of the original image
**
**         [-g] \n");
**                add a simulated noise to the decompressed
**                image with the same properties as in the
**                original image. So they look very similar.
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_decomp.cc 3.3 96/05/07 CEA 1995 @(#)";

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Noise.h"
#include "IM_Comp.h"
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

    mrdecomp_option_usage();
    manline();
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
    while ((c = GetOpt(argc,argv,"vr:t:gI:B:zZ:")) != -1) 
    {
	switch (c) 
        {
           case 'I': if (sscanf(OptArg,"%d",&IterRec) != 1) 
                     {
		    fprintf(OUTMAN, "Error: bad iteration number parameter: %s\n", OptArg);
		    exit(-1);
		     } 
		     break;
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
            case 't':
               /* output can be interger or short */
		if (sscanf(OptArg,"%c",&OutputType) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad output type: %s\n", OptArg);
		   exit(-1);
		}
		if ((OutputType != 'i') && (OutputType != 's') && (OutputType != 'f'))
		{
		    fprintf(OUTMAN, "Error: bad output type: %s\n", OptArg);
		    fprintf(OUTMAN, "       output type has to equal to 'i','s' or 'f'\n");
		    exit(-1);
		}
		break;
            case 'g':
		/*simulates the noise*/
		SimuNoise=True;
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

    hcinit(argc, argv);

 MR_CompData DecIma;
    DecIma.IterRec = IterRec;
    DecIma.Cmd = Cmd;
    DecIma.File_Name_Imag=File_Name_Imag;
    DecIma.File_Name_Transform=File_Name_Transform;
    DecIma.Resol=Resol;
    DecIma.OutputType=OutputType;
    if (SimuNoise == True) DecIma.SimuNoise=True;
    else DecIma.SimuNoise=False;
    if (Verbose==True) DecIma.Verbose = True;

    // apply the decompression
    if (NumBlock > 0) DecIma.decompress(NumBlock);
    else DecIma.decompress();
    
    // save the result (field Dat) in DecIma.File_Name_Imag
    // in case where the full image is decompressed with block compression
    // options, the result is already in the output file.
    if ((NumBlock > 0) || (DecIma.UseBlock == False)) DecIma.write();
    exit(0);
}



 


