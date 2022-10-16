/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  04/09/99
**    
**    File:  mr3d_insert.cc
**
*******************************************************************************
**
**    DESCRIPTION  insert a band in a multiresolution transform  
**    ----------- 
**                 
**    Usage: mr3d_insert options cube output
**
******************************************************************************/


#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"

/****************************************************************************/

char Name_Cube_In[256];
char Name_Cube_Out[256];
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose=False;
int BandNumber=1; // band number of extract

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_out_mr_file in_cube\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "              Band number to insert.\n");
    manline();    
    verbose_usage();    
    manline();   

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
    while ((c = GetOpt(argc,argv,"b:vzZ:")) != -1) 
    {
	switch (c) 
        {
	    case 'b':
 		if (sscanf(OptArg,"%d",&BandNumber) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad band number: %s\n", OptArg);
		   exit(-1);
		}
                if (BandNumber <= 0) BandNumber = 1;
 		break;		
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
 	   case '?':
			usage(argv);
		}
	}
     
	BandNumber --;
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Cube_Out, argv[OptInd++]);
        else usage(argv);

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

/*********************************************************************/

int main(int argc, char *argv[])
{
   fitsstruct Header;
   char Cmd[512];
   extern softinfo Soft;
   Soft.mr3();
   int k;
   
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    transinit(argc, argv);
 
    if (Verbose == True)
    {
       cout << "MR Filename = " << Name_Cube_In << endl;
       cout << "Filename in = " << Name_Cube_Out  << endl;
       cout << "Band number to extract = " << BandNumber+1 << endl;
    }    
    fltarray DatBand;
    MR_3D MR_Data;

    if (Verbose == True) cout << "Read ... " << endl;
    MR_Data.read(Name_Cube_In);
    int Nx = MR_Data.size_cube_nx();
    int Ny = MR_Data.size_cube_ny();
    int Nz = MR_Data.size_cube_nz();
    if (Verbose == True) 
    {
       cout << "Number of bands = " << MR_Data.nbr_band() << endl;
       cout << "Nx = " << Nx << " Ny = " << Ny << " Nz = " << Nz << endl;
    }
    if ((BandNumber < 0) || (BandNumber >= MR_Data.nbr_band()))
    {
       cerr << "Error: bad band number ... " << endl;
       exit(-1);
    }
    if (Verbose == True) 
    {
       cout << "Band size = " << MR_Data.size_band_nx(BandNumber);
       cout << ", " << MR_Data.size_band_ny(BandNumber);
       cout << ", " << MR_Data.size_band_nz(BandNumber) << endl;
    }    
    
    io_3d_read_data(Name_Cube_Out, DatBand);
    if (Verbose == True) cout << "insert band ... " << endl;

    MR_Data.insert_band(BandNumber, DatBand);
    if (Verbose == True) cout << "Write ... " << endl;
    MR_Data.write(Name_Cube_In);
    exit(0);
}

