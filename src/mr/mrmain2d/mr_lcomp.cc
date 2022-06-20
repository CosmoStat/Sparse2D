/******************************************************************************
**                   Copyright (C) 1998 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  mr_lcomp.cc
**
*******************************************************************************
**
**    DECRIPTION  compression program
**    ---------- 
**
*******************************************************************************
**
**    PARAMETRES
**    ----------
**
**        Input_File: image file name to compress
**        [Output_File]: compressed image file
**                       by default, "Input_File.MRC"
**           [-v] : vebose
**
**
**           [-n number of scales] : number of scales used by the 
**                                 multiresolution transform
**                                 by default number_of_scales = 4
**      
******************************************************************************/

 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"

char File_Name_Imag[256];  /* input file image */
char File_Name_Transform[256]; /* output file name */
int Nbr_Plan=DEFAULT_CMP_NBR_SCALE;  /* number of scales */
int KeepFitsHeader = 0;
Bool NoBscale=False; // for integer coding. Fits image are not divided
                     // by the BScale fits keyword
Bool Verbose=False;
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

#define WRITE_PARAM 0
 
type_comp Comp_Method = COMP_INT_LIFTING;
 
Bool UseBudget=False;
Bool NoiseInData = False;
int BlockSize=0;
Bool UseBlock=False;
type_lift LiftingTrans = DEF_LIFT;

#define NBR_OK_METHOD 5
static type_lift TabMethod[NBR_OK_METHOD] = {TL_MEDIAN, TL_INT_HAAR, TL_INT_CDF, 
                                  TL_INT_4_2, TL_INT_F79};

/*********************************************************************/

static void usage(char *argv[])
{
   int i;
   
    fprintf(OUTMAN, "Usage: %s options in_image [out_file]\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-m Compression_Method]\n");
    for (i = 0; i < NBR_OK_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringLSTransform(TabMethod[i]));
    fprintf(OUTMAN, "              default is %s\n",StringLSTransform (LiftingTrans));
    manline();    

    nbr_scale_usage(Nbr_Plan);
    manline();
    fprintf(OUTMAN, "         [-f ]\n");
    fprintf(OUTMAN, "              Keep all the fits header. Default is no.\n");
    manline();    

    fprintf(OUTMAN, "         [-B]\n");
    fprintf(OUTMAN, "              do not apply BSCALE operation\n");
    fprintf(OUTMAN, "              in case of fits images.\n");
    manline();    
    fprintf(OUTMAN, "         [-C BlockSize]\n");
    fprintf(OUTMAN, "              Compress by block. \n");
    fprintf(OUTMAN, "              BlockSize = size of each block.\n");
    fprintf(OUTMAN, "              Default is No.\n");
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
    /* get options */
    while ((c = GetOpt(argc,argv,"m:vn:fBC:zZ:")) != -1) 
    {
	switch (c) 
        {
           case 'C': 
              if (sscanf(OptArg,"%d",&BlockSize ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
	            exit(-1);
                    
		}
		UseBlock = True;
                break;
	   case 'm':
		/* -f <type> type of compression */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of compression: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <=  NBR_OK_METHOD)) 
		    LiftingTrans = (type_lift) (TabMethod[c-1]);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of compression: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	    case 'v':
		/* Verbose flag -v */
		Verbose = True;
            case 'f':
                /* keep all keywords in fits header */
                KeepFitsHeader = 1;
               break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
		    exit(-1);
		}
 		break;
            case 'B':
                NoBscale = True;
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
	if (OptInd < argc) 
             strcpy(File_Name_Imag, argv[OptInd++]);
         else usage(argv);

        strcpy (File_Name_Transform, File_Name_Imag);

	if (OptInd < argc) 
             strcpy(File_Name_Transform, argv[OptInd++]);

	L = strlen (File_Name_Transform);
	Ptr = File_Name_Transform;
	if (L < 5) strcat (File_Name_Transform, ".MRC");
	else if ((Ptr[L-1] != 'C') || (Ptr[L-2] != 'R') || (Ptr[L-3] != 'M')) 
               strcat (File_Name_Transform, ".MRC");

  
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***************************************************/

int main(int argc, char *argv[])
{
    int k;
    char Cmd[256];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, 
       open input file(s) if necessary */
    lm_check(LIC_MR1);
    hcinit(argc, argv);
    MR_CompData CompIma;
    
    if (UseBlock == True)
    {
         CompIma.UseBlock = True;
         CompIma.BlockSize = BlockSize;
    }
    CompIma.NoiseInData = NoiseInData;
    CompIma.NoBscale=NoBscale;
    CompIma.Cmd=Cmd;
    CompIma.File_Name_Imag=File_Name_Imag;
    CompIma.File_Name_Transform=File_Name_Transform;
    CompIma.Nbr_Plan=Nbr_Plan;
    CompIma.KeepFitsHeader = KeepFitsHeader;
    CompIma.Comp_Method = Comp_Method;
    CompIma.LiftingTrans = LiftingTrans;
    if (Verbose == True) CompIma.Verbose = True;
    CompIma.compress();

    exit(0);
}


/***************************************************/

 

