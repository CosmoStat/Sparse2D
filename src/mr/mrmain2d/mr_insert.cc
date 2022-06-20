/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  mr_insert.cc
**
*******************************************************************************
**
**    DESCRIPTION  insert an image in a multiresoluiton transform
**    ----------- 
**                 
**    PARAMETRES    
**    ----------  
**
**    USAGE: mr_insert option multiresolution_file input_image
**        multiresolution_file = file (.mr) which contains the multiresolution
**                               transform.
**        input_image = input scale file name
**  
**        where options = 
**   
**           [-s scale_number]
**                scale number to insert. default is 1.
** 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"

char Name_MR_In[256];
char Name_Imag_In[256];
int ScaleNumber = 1;
int BandNumber = 1;
Bool InsertBand = False;
Bool InsertScale = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_out_mr_file in_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options are =  \n");
    manline();

    read_band_usage();
    manline();    
    scale_number_insert_usage();
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
static void insinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif        
    /* get options */
    while ((c = GetOpt(argc,argv,"s:b:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
  	   case 's':
		/* -n <scale number> */
		if (sscanf(OptArg,"%d",&ScaleNumber) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
		    exit(-1);
		}
                if (ScaleNumber <= 0)   
                 {
		    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
  		    exit(-1);
		}
                InsertScale = True;
                if (InsertBand  == True)
                {
                     fprintf(OUTMAN, "Error: options -b and -s are exclusive ... \n");
                     exit(-1);
                }
		break;
	    case 'b': 
		/* -n <scale number> */
		if (sscanf(OptArg,"%d",&BandNumber) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
		    exit(-1);
		}
                if (BandNumber <= 0)   
                 {
		    fprintf(OUTMAN, "Error: bad band number: %s\n", OptArg);
  		    exit(-1);
		}
                InsertBand = True; 
                if (InsertScale == True)
                {
                     fprintf(OUTMAN, "Error: options -b and -s are exclusive ... \n");
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
	if (OptInd < argc) strcpy(Name_MR_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
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
    MultiResol MR_Data;
    Ifloat Dat;
    int s;

    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    insinit(argc, argv);

    if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Name in/out = " << Name_MR_In  << endl;
    }

    MR_Data.Verbose = Verbose;
    MR_Data.read (Name_MR_In);
    if (InsertScale == False) InsertBand = True;
    if (Verbose == True) MR_Data.print_info();

    if (InsertBand == False)
    {
       s = ScaleNumber;
       if ((s < 1) || (s > MR_Data.nbr_scale()))
       {
           cerr << "Error: illegal scale number ..." << endl;
           cerr << " 0 < scale_number <= " << MR_Data.nbr_scale() << endl;
           exit (0);
       }
       MR_Data.read(Name_Imag_In, s-1);
    }
    else
    {
       s = BandNumber;
       if ((s < 1) || (s > MR_Data.nbr_band()))
       {
           cerr << "Error: illegal band number ..." << endl;
           cerr << " 0 < band_number <= " << MR_Data.nbr_band() << endl;
           exit (0);
       }
       MR_Data.read_band(Name_Imag_In, s-1);
    }    
    MR_Data.write (Name_MR_In);
    exit(0);
} 

