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
**    File:  mr3d_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube
**    ----------- 
**                 
**    Usage: mr3d_trans options cube output
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "MR1D_Obj.h"
#include "MR_Obj.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
#include "DefFunc.h"
#include "MR2D1D.h"

/****************************************************************************/
#define TESTIND 1

 
/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/

/****************************************************************************/

char Name_Cube_In[256];
char Name_Out[256];
 
int NbrScale2d = 5;
int Nbr_Plan=4;
type_transform  Transform=TO_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

// sb_type_norm Norm = NORM_L2;
// type_sb_filter SB_Filter = F_HAAR; // F_MALLAT_7_9;
// type_lift LiftingTrans = DEF_LIFT;
Bool Verbose=False;
Bool GetMax = False;
float N_Sigma=3;                // number of sigma (for the noise)  
Bool Normalize=False;       // normalize data in
 
type_border Bord = I_MIRROR;
Bool Reverse = False;

/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options cube output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    all_transform_usage(Transform);
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales_2D]\n");
    fprintf(OUTMAN, "              number of scales used in the 2D wavelet transform\n");
    manline();
    fprintf(OUTMAN, "         [-N number_of_scales_1D]\n");
    fprintf(OUTMAN, "              number of scales used in the 1D wavelet transform\n");
    manline();
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Normalize the data. Default is no. \n");    
    manline();
    manline();
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Reconstruction. Default is no. \n");    
    manline();       
    verbose_usage();    
    manline();
    vm_usage();
    manline();    
    manline();
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
    while ((c = GetOpt(argc,argv,"rN:t:n:MvzZ:")) != -1) 
    {
	switch (c) 
        {
           case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM+1)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}                
 		break;
	   case 'r': Reverse = True; break;
       case 'M': Normalize = (Normalize == True) ? False: True; break;
 	   case 'v': Verbose = True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrScale2d) != 1) 
                {
		    fprintf(OUTMAN, "bad number of 2D scales : %s\n", OptArg);
		    exit(-1);
		}
                if ((NbrScale2d <= 1) || (NbrScale2d > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
  		    exit(-1);
		}
		break;
	case 'N':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of 1D scales : %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 0) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
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
 	   case '?':
			usage(argv);
		}
	}
      
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Out, argv[OptInd++]);
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
    fltarray Dat;
    /* Get command line arguments, open input file(s) if necessary */
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   
    lm_check(LIC_MR3);
    transinit(argc, argv);
 

    
    if (Reverse == False)
    {
        if (Verbose == True)
        {
           cout << "Filename in = " << Name_Cube_In << endl;
           cout << "Filename out = " << Name_Out  << endl;
           cout << "Transform = " << StringTransform((type_transform) Transform) << endl;
           cout << "NbrScale2d = " << NbrScale2d<< endl;
           cout << "NbrScale1d = " << Nbr_Plan<< endl;
            if (Normalize == True) cout << "Normalize the input data: " << "yes" << endl;
            else cout << "Normalize the input data: " << "no" << endl;
        }
    
       io_3d_read_data(Name_Cube_In, Dat, &Header);
 
       int Nx = Dat.nx();
       int Ny = Dat.ny();
       int Nz = Dat.nz();
       if (Verbose == True) cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
     
       if (Normalize == True)
       {
         double Mean = Dat.mean();
         double Sigma = Dat.sigma();
         for (int i=0;i<Nx;i++)
         for (int j=0;j<Ny;j++)
         for (int k=0;k<Nz;k++) Dat(i,j,k) = (Dat(i,j,k)-Mean)/Sigma;
       }    
    
       MR2D1D WT;
       if (Verbose == True) cout << "Alloc ...  " << endl;
       WT.alloc(Nx, Ny, Nz, Transform, NbrScale2d, Nbr_Plan);

       if (Verbose == True) cout << "Transform ...  " << endl;
       WT.transform (Dat);

       if (Verbose == True)cout << "Write result ...  " << endl;
       WT.write(Name_Out);
              
       if (Verbose == True)
       {
          for (int s2 = 0; s2 < WT.nbr_band_2d (); s2++)
	      for (int s1 = 0; s1 < WT.nbr_band_1d (); s1++)
	      {
	          cout << "  Band " << s2 << ", " << s1 << ": " << " Nx = " << WT.size_band_nx(s2,s1) << ", Ny = " << WT.size_band_ny(s2,s1) <<  ", Nz = " << WT.size_band_nz(s2,s1) << endl;
	          fltarray Band;
	          Band = WT.get_band(s2, s1);
	          cout << "  Sigma = " << Band.sigma() << " Min = " << Band.min() << " Max = " << Band.max() << endl;
        // char Name[512];
        // sprintf(Name, "Band_%d_%d.fits", s2+1, s1+1);
        //  fits_write_fltarr(Name, Band);
	      }
       }
//        cout << endl << "READ " << endl;
//        MR2D1D WT1;
//        WT1.read(Name_Out);
//        fltarray Dat1;
//        WT1.recons (Dat1);
//        fits_write_fltarr (Name_Out, Dat1);
    }
    else
    {
       MR2D1D WT;
       WT.read(Name_Cube_In);
       WT.recons (Dat);
 
       if (Verbose == True) cout << "Write result ...  " << endl;
       fits_write_fltarr(Name_Out, Dat);
    }
    // Header.origin = Cmd;	 
    // fits_write_fltarr (Name_Out, StaInfo, &Header);
     
    exit(0);
}
