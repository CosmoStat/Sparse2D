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
**    File:  mr_extract.cc
**
**    History : modification by R Gastaud , 18 Feb 1997, add write
**
*******************************************************************************
**
**    DESCRIPTION  extract a scale from a multiresolution transform
**    ----------- 
**
**    PARAMETRES    
**    ---------- 
**
**    USAGE: mr_extract option multiresolution_file  output_image
**        multiresolution_file = file (.mr) which contains the multiresolution
**                               transform.
**        output_image = output scale file name
**  
**        where options = 
**   
**           [-s scale_number]
**                scale number to extract. default is 1.
**
**           [-a]
**                extract all the scales
**
**           [-b] 
**                interpolate by block the scale, in order to have the same
**                size as the original image
**                This option is valid only if the choosen multiresolution 
**                transform is pyramidal (6,7,8,9,10,11,12)
**
**           [-i] 
**                interpolate by B3-spline the scale ,in order to have the same
**                size as the original image.
**                This option is valid only if the choosen multiresolution 
**                transform is pyramidal (6,7,8,9,10,11,12)
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include <iostream>

char Name_Imag_Out[256]; /* output file name */
char Name_Imag_In[256]; /* input file image */
Bool WriteScaleb = False;          /* interpolate with b option */
Bool WriteScalei = False;          /* interpolate with i option */
Bool WriteAll = False;
int ScaleNumber = 1;
int BandNumber = 1;
Bool ExtractBand = False;
Bool ExtractScale = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_mr_file out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

    write_band_usage();
    manline();
    scale_number_usage();
    manline();
    write_scales_x_band_usage();
    manline();
    fprintf(OUTMAN, "         [-B]\n");
    fprintf(OUTMAN, "             Interpolate the extracted band by block to the original image size.\n");
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
static void extrinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
    
    /* get options */
    while ((c = GetOpt(argc,argv,"xs:b:BvzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
  	   case 's':
		/* -n <scale number> */
		if (sscanf(OptArg,"%d",&ScaleNumber) != 1) 
                {
		    fprintf(OUTMAN, "bad scale number: %s\n", OptArg);
		    usage(argv);
		}
                if (ScaleNumber <= 0) 
                 {
		    fprintf(OUTMAN, "bad scale or band number: %s\n", OptArg);
  		    usage(argv);
		}
                ExtractScale = True;
		break;
	    case 'b':
                if (sscanf(OptArg,"%d",&BandNumber) != 1) 
                {
		    fprintf(OUTMAN, "bad band number: %s\n", OptArg);
		    usage(argv);
		}
                if (BandNumber  <= 0) 
                 {
		    fprintf(OUTMAN, "bad   band number: %s\n", OptArg);
  		    usage(argv);
		}
		ExtractBand = True;
  		break;
            case 'B': WriteScaleb = True; break;
            case 'x': WriteAll = True; ExtractBand = True; break;
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
        if ((ExtractScale == True) && 
            ((ExtractBand == True) || (WriteAll  == True) || (WriteScaleb == True)))
        {
           fprintf(OUTMAN, "Error: -s option is not compatible b,B,x options ... \n");
           exit(-1);
        }

       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) 
             strcpy(Name_Imag_Out, argv[OptInd++]);
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
    int s,Nl,Nc,Nbr_Plan,NbrBand;
    type_transform Transform;
    set_transform Set_Transform;

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    extrinit(argc, argv);
    if (ExtractScale == False) ExtractBand = True;
 
    if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Name in = " << Name_Imag_In << endl;
      cout << "File Name Out = " << Name_Imag_Out << endl;
    }
    MR_Data.Verbose = Verbose;
    MR_Data.read (Name_Imag_In);
    Nl = MR_Data.size_ima_nl();
    Nc = MR_Data.size_ima_nc();
    Nbr_Plan = MR_Data.nbr_scale();
    NbrBand = MR_Data.nbr_band();
    Transform = MR_Data.Type_Transform;
    Set_Transform = MR_Data.Set_Transform;

    if (Verbose == True) MR_Data.print_info();
 
    int StartScale = (ExtractScale == True) ? ScaleNumber-1: BandNumber -1;
    int EndScale =  (ExtractScale == True) ? ScaleNumber-1: BandNumber -1;
    if (WriteAll == True)
    {
       StartScale = 0;
       EndScale = MR_Data.nbr_band()-1;
    }
           
   /*********   check the arguments  *************************/

    if ((MR_Data.FormatInputImag != F_UNKNOWN) &&
        (io_which_format(Name_Imag_Out) == F_UNKNOWN))
                             io_set_format(MR_Data.FormatInputImag);
  
     if (WriteAll == False)
     {
        if ((ExtractBand == True) && ((BandNumber < 1) || (BandNumber > NbrBand)))
        {
           cerr << "Error: illegal band number " << endl;
           cerr << " 0 < band_number <= " << Nbr_Plan << endl;
           exit (0);
        }
        if ((ExtractScale == True) && ((ScaleNumber < 1) || (ScaleNumber > Nbr_Plan)))
        {
           cerr << "Error: illegal scale number " << endl;
           cerr << " 0 < scale_number <= " <<  Nbr_Plan << endl;
           exit (0);
        }
     }
 
     /*********   write the output data *************************/

    if (ExtractScale == True)
    {
        char NameImag[80];
        strcpy(NameImag, Name_Imag_Out);
	io_set_format(Name_Imag_Out);
        MR_Data.write(NameImag,ScaleNumber-1);
    }
    else for (s = StartScale; s <= EndScale; s++)
    {
       char NameImag[80];
       if (StartScale == EndScale) strcpy(NameImag, Name_Imag_Out);
       else // create a filename from the band number
       {
          io_strcpy_prefix(NameImag,Name_Imag_Out);
          sprintf(NameImag, "%s_band_%d", NameImag, s+1);
	  io_set_format(Name_Imag_Out);
        }
        if (Verbose == True) cout << "Create " << NameImag << endl;
        if (WriteScaleb==True) 
        {
          Ifloat Dats (Nl, Nc, "Interp");
          im_block_extend(MR_Data.band(s), Dats);
          io_write_ima_float(NameImag, Dats);
        }
        else  MR_Data.write_band(NameImag,s);
    }
    //exit(0);
} 

