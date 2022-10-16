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
**    File:  mr1d_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform  an a signal
**    ----------- 
**                 
**    Usage: mr1d_trans options image output
**
**   where options =
**
**         [-t type_of_multiresolution_transform]
**              1: linear wavelet transform: a trous algorithm
**              2: b1spline wavelet transform: a trous algorithm
**              3: b3spline wavelet transform: a trous algorithm
**              4: Morlet's wavelet transform
**              5: morphological median transform
**              6: mexican hat wavelet transform
**              7: french hat wavelet transform
**              8: pyramidal linear wavelet transform
**              9: pyramidal b3spline wavelet transform
**              10: pyramidal median transform
**         [-n number_of_scales]
**              number of scales used in the multiresolution transform
**        [-r]
**              rebin all scales to the input size
**              (for pyramidal median transform only)
**
**
******************************************************************************/

/*
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR1D_Obj.h"
#include "IM_Errors.h"
#include "MR_NoiseModel.h"
#include "MR_Sigma.h"
*/

#include "GlobalInc.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"

#define DEF_ITER 10

char Name_Imag_In[256];
char Name_Imag_Out[256];
char NameInfo[256];

int Nbr_Plan=0;
type_trans_1d Transform=TO1_PAVE_B3SPLINE;
Bool Rebin=False;
Bool KillBord=False;
Bool WriteInfo = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_lift LiftingTrans = DEF_LIFT;
Bool Verbose = False;
type_undec_filter U_Filter = DEF_UNDER_FILTER;

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
    int i;

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_TRANS_1D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf1D((type_trans_1d)i));

    manline();
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();
    
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "              rebin all scales to the input size\n");
    fprintf(OUTMAN, "              (for pyramidal transforms only)\n");
    manline();
    
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "              set to 0 pixels contaminated by the border problem.\n");
    manline();
    
    sb_usage(SB_Filter);
    manline();
    usb_usage(U_Filter);
    manline();
    
    lifting_usage(LiftingTrans);
    manline();
    
    fprintf(OUTMAN, "         [-w InfoFileName]\n");
    fprintf(OUTMAN, "             write in a file the size and the starting index of each band.\n");
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
    Bool OptL = False, Optf = False;

    /* get options */
    while ((c = GetOpt(argc,argv,(char*)"U:vt:n:rkw:Ll:T:")) != -1) 
    {
	switch (c) 
        {
	   case 'U': 
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "ErrorL bad type of filters: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_UNDEC_FILTER)) 
                                        U_Filter = (type_undec_filter) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of filters: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	   case 'v': Verbose = True; break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
	   case 'l':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_LIFT)) LiftingTrans = (type_lift) (c);
                else  
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	   case 'T':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad filter: %s\n", OptArg);
	            exit(-1);
                    
		}
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
 		break;
	   case 'w': 
	      if (sscanf(OptArg,"%s",NameInfo ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	            exit(-1);
                    
		}
		WriteInfo = True; break;
           case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANS_1D)) 
                                        Transform = (type_trans_1d) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    exit(-1);
		}
		break;
 	   case 'r':
                /* rebin */
                Rebin=True;
                break;
	   case 'k':
	        KillBord=True;
                break;
	   case '?':
			usage(argv);
		}
	}

//        if ( (WriteInfo == True) && 
//             (which_set_is_trans1d(Transform) != TRANS1_MALLAT)
// 	  && (which_set_is_trans1d(Transform) != WP1_MALLAT))
//        {
//           fprintf(OUTMAN, "Error:  .\n");
//           exit(-1);;
//        }       
       
       if ( (Rebin== True) && (which_set_is_trans1d(Transform) != TRANS1_PYR))
       {
          fprintf(OUTMAN, "Error: rebin option must be used with pyramidal transforms.\n");
          exit(-1);;
       }
        if (((Transform != TO1_MALLAT) && (Transform != TU1_MALLAT) && (Transform != WP1_MALLAT)) 
	    && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with transforms 7, 15 and 17  ... \n");
           exit(0);
	}
 	if (((Transform != TO1_LIFTING) && (Transform != WP1_LIFTING)) 
	    && ( LiftingTrans != DEF_LIFT))
	{
	   fprintf(OUTMAN, "Error: option -l is only valid with transforms 16 and 18  ... \n");
           exit(0);
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
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    fltarray Dat;
    fitsstruct FitsHeader;
    extern type_1d_format IO_1D_Format;

    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_M1D);
    transinit(argc, argv);
    char Cmd[256];
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
  
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "Transform = " << StringTransf1D(Transform) << endl;
       if ((Transform == TO1_MALLAT) || (Transform == TU1_MALLAT) ||(Transform ==  WP1_MALLAT))
       {
          cout << StringSBFilter(SB_Filter) << endl;
          if (Norm == NORM_L2) cout << "L2 normalization" << endl;
       }
       if ((Transform == TO1_LIFTING) || (Transform == WP1_LIFTING)) 
          cout << "Lifting method = " <<  StringLSTransform(LiftingTrans) << endl;
       if (Transform == TU1_UNDECIMATED_NON_ORTHO)  cout << "Undec. Filter Bank = " << StringUndecFilter(U_Filter) << endl;

    }

    // fits_read_fltarr(Name_Imag_In, Dat);
    io_1d_read_data(Name_Imag_In, Dat, &FitsHeader);
    FitsHeader.origin = Cmd;
    reform_to_1d(Dat);
    int Np = Dat.nx();
       
    //if (Nbr_Plan ==0) Nbr_Plan = max_scale_number(Dat.n_elem());
    //cout << "Number of scales = " << Nbr_Plan << endl;

    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TU1_MALLAT) || (Transform == TO1_MALLAT) || (Transform == WP1_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_1D MR_Data;
    MR_Data.U_Filter = U_Filter;
    MR_Data.alloc (Dat.n_elem(), Transform, Nbr_Plan , PtrFAS, Norm, Rebin);
    // MR_1D MR_Data (Dat.n_elem(), Transform, "MR_Transform", Nbr_Plan, Rebin);
    if (KillBord == True) MR_Data.KillBorder = True;
    if ((Transform == TO1_LIFTING)  || (Transform == WP1_LIFTING))
             MR_Data.LiftingTrans = LiftingTrans;

//     if ((Transform == TO1_MALLAT) || (Transform == WP1_MALLAT))
//     {
//         MR_Data.SB_Filter = SB_Filter;
// 	MR_Data.Norm = Norm;
//     }
 
    if (Verbose == True)
        cout << "Number of scales = " << MR_Data.nbr_scale() << endl;
    
    MR_Data.transform (Dat);

//     for (int s =0; s <   MR_Data.nbr_band(); s++)
//     {
//        cout << "Band " << s << " Border size = " << MR_Data.border_size(s);		   
//        MR_Data.scale(Dat, s);
//        cout << " sigma = " << Dat.sigma() << endl;
//     }
    
    // fits_write_fltarr(Name_Imag_Out, MR_Data.image());
    Dat = MR_Data.image();

    if (IO_1D_Format == F1D_FITS)
            fits_write_fltarr(Name_Imag_Out, Dat);
    else io_write2d_ascii(Name_Imag_Out, Dat);
       
       
           int SizeUserName = 0;
    if (SB_Filter == F_USER)
    {
       if (UserFilterFileName == NULL)
       {
          cout << "Error: user filter name is not defined ... " << endl;
          exit(-1);
       }
       SizeUserName = strlen(UserFilterFileName);
    }
    if (WriteInfo == True)
    {
       fltarray Info(2, MR_Data.nbr_band()+4+SizeUserName);
       Info(0,0) = (int) Transform + 1;
       Info(1,0) =  MR_Data.nbr_scale();
       Info(0,1) = (int) MR_Data.Set_Transform + 1;
       Info(1,1) =  MR_Data.nbr_band();
       Info(0,2) =  Np;
       Info(1,2) = (int) MR_Data.LiftingTrans;
       if (Transform != TU1_UNDECIMATED_NON_ORTHO) Info(0,3) = (int) MR_Data.SB_Filter;
       else  Info(0,3) = (int) MR_Data.U_Filter;
       Info(1,3) = (int) MR_Data.Norm;
       
       for (int s =1; s <= MR_Data.nbr_band(); s++)
       {
          Info(0,s+3) = MR_Data.size_scale_np(s-1);
	  Info(1,s+3) = MR_Data.pos_band(s-1);
       }
       int FirstChar = MR_Data.nbr_band() + 4;
       if (SB_Filter == F_USER)
       {
          Info(1,FirstChar) = SizeUserName;
          // cout << " SizeUserName = " << SizeUserName << endl;
          // cout << " FirstChar  = " <<  FirstChar << endl;

          for (int c =0; c < SizeUserName; c++)
          {
             Info(0,c+FirstChar) = (int) UserFilterFileName[c];
	     if (c > 0) Info(1,c+FirstChar) = 0.;
          }
       }
       if (IO_1D_Format == F1D_FITS)
            fits_write_fltarr(NameInfo, Info);
       else io_write2d_ascii(NameInfo, Info);
    }   
    exit(0);
} 

