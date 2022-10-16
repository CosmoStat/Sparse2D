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
**           [-v] : verbose
**
**		"Usage: %s [-v] [-r resolution] input output\n",
**        [-v] : verbose
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

int Resol=0;
char File_Name_Imag[80], File_Name_Transform[80];

int Verbose;
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

extern Bool InputFromStdin;
int IterRec=0;
int NumBlock=-1;
Bool GetDebug=False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options [Input MR File] [Output coeff file]\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    mrdecomptool_option_usage();
    manline();
    fprintf(OUTMAN, "        [-r resolution]\n");
    fprintf(OUTMAN, "          resolution  = -1,0..nbr_of_scale] \n");
    fprintf(OUTMAN, "          resolution  = 0 for full resolution (default)  \n");
    fprintf(OUTMAN, "          resolution = -1 for residual map.\n");
    fprintf(OUTMAN, "          resolution = nbr_of_scale-1 for the worse resol.\n");
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
 
    Verbose = 0;
    strcpy(File_Name_Transform, "stdin");
    InputFromStdin = True;

    /* get options */
    while ((c = GetOpt(argc,argv,"xvr:B:")) != -1) 
    {
	switch (c) 
        {
	   case 'x': GetDebug=True; break;
           case 'I': if (sscanf(OptArg,"%d",&IterRec) != 1) 
                     {
		    fprintf(OUTMAN, "bad iteration number parameter: %s\n", OptArg);
		    usage(argv);
		     } 
		     break;
	   case 'v':
		/* verbose flag -v */
		Verbose = 1;
		break;
            case 'B': 
                if (sscanf(OptArg,"%d",&NumBlock) != 1)
                {
		   fprintf(OUTMAN, "bad block parameter: %s\n", OptArg);
		   usage(argv);
		}
		break;
             case 'r': /* -r < resolution> */
		if (sscanf(OptArg,"%d",&Resol) != 1) 
                {
		    fprintf(OUTMAN, "bad Resolution parameter: %s\n", OptArg);
		    usage(argv);
		}
                if ((Resol < -1) || (Resol > 20))
                {
		    fprintf(OUTMAN, "bad Resolution parameter: %s\n", OptArg);
		    usage(argv);
		}
               break;
	   case '?':
		usage(argv);
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

    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
	fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
    }
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
    if (Verbose) DecIma.Verbose = True;
      
     if (NumBlock < 1) NumBlock=1;
     
     char *Buff;
     int BuffSize;
     DecIma.get_resol(Resol, BuffSize, Buff, NumBlock);
     
     if (Verbose) cerr << "BuffSize = " << BuffSize << endl;
     
     FILE *Output;
     if (std_inout(File_Name_Imag) == True) Output = stdout;
     else 
     {
        if (!(Output = fopen(File_Name_Imag, "wb")))
        {
           fprintf(stderr,"Error writing on file %s\n", File_Name_Imag);
           exit (-1);
        }
    }
 
    writeint(Output, BuffSize);
    int Nbr = fwrite ((void *) Buff, sizeof (char),  BuffSize, Output);
    if (Nbr <= 0)
    {
        fprintf(stderr,"Error writting in file %s\n", File_Name_Imag);
        exit (-1);
    }

    if (Output != stdout) fclose(Output);
    
    if (GetDebug == True)
    {
       int *Tab;
       int Nl, Nc, ind=0;
       float scale;
       Output = fopen(File_Name_Imag, "rb");
       decode(Output, &Tab, &Nc,&Nl, &scale);
       cout << "Nl = " << Nl << " Nc = " << Nc << " scale = " << scale << endl;
       Ifloat Band(Nl,Nc,"band");
       if (scale > FLOAT_EPSILON)
       {
          for (int i=0; i < Nl; i++)
          for (int j=0; j < Nc; j++) Band(i,j) = Tab[ind++]*scale;
       }
       else
       {
          for (int i=0; i < Nl; i++)
          for (int j=0; j < Nc; j++) Band(i,j) = Tab[0];
       }
       io_write_ima_float("xx_band.fits", Band);
    }
    exit(0);
}



 


