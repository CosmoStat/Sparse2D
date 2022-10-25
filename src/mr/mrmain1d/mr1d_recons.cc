/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  23/08/98
**    
**    File:  mr1d_recons.cc
**
*******************************************************************************
**
**    DESCRIPTION  reconstruction of a signal from its transform
**    ----------- 
**                 
**    Usage: mr1d_recons transform_file info_file output
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
// #include "MR_Sigma.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char NameInfo[256];

int Nbr_Plan=0;
type_trans_1d Transform=TO1_PAVE_B3SPLINE;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose = False;
type_undec_filter U_Filter = DEF_UNDER_FILTER;

/*********************************************************************/

static void usage(char *argv[])
{

    fprintf(OUTMAN, "Usage: %s  TransformFileName InfoFileName OutputFileName\n\n", argv[0]);
    // fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"v")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
	   case '?': usage(argv); break;
	   default: usage(argv);break;
	}
   }
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(NameInfo, argv[OptInd++]);
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
    fltarray Dat,Info;
    sb_type_norm Norm = NORM_L1;
    type_sb_filter SB_Filter = F_MALLAT_7_9;
    type_lift LiftingTrans = DEF_LIFT;

    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_M1D);
    transinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name out = " <<  Name_Imag_Out << endl;
       cout << "Info file = " <<  NameInfo  << endl;
   }

    //Dat.fits_read(Name_Imag_In);
    //Info.fits_read(NameInfo);
    type_1d_format  DatFormat = io_detect_1dformat(Name_Imag_In);
    if (DatFormat == F1D_FITS)
    {        
       fits_read_fltarr(Name_Imag_In, Dat);
       fits_read_fltarr(NameInfo, Info);
    }
    else
    {
       io_read2d_ascii(Name_Imag_In, Dat);
       io_read2d_ascii(NameInfo, Info);
    }
    
     //  fits_read_fltarr(Name_Imag_In, Dat);
    Transform = (type_trans_1d) (Info(0,0) - 1);
    int Np = (int) Info(0,2);
    Nbr_Plan = (int) Info(1,0);
    int NbrBand = (int) Info(1,1);
    Norm =  (sb_type_norm)  Info(1,3);
    if (Transform != TU1_UNDECIMATED_NON_ORTHO)  SB_Filter = (type_sb_filter) Info(0,3);
    else U_Filter = (type_undec_filter) Info(0,3);
    
    LiftingTrans =  (type_lift) Info(1,2);
    int FirstChar = NbrBand + 4;
    if (SB_Filter == F_USER)
    {
        int SizeUserName = (int) Info(1,FirstChar);
        // cout << " SizeUserName = " << SizeUserName << endl;
        // cout << " FirstChar  = " <<  FirstChar << endl;
        UserFilterFileName = new char[SizeUserName+1];
        for (int c =0; c < SizeUserName; c++)
        {
            UserFilterFileName[c] = (char) Info(0,c+FirstChar);
        } 
        UserFilterFileName[SizeUserName] = '\0';
        // cout << "UserFilterFileName = " << UserFilterFileName << endl;
    }    

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
    MR_Data.alloc (Np, Transform, Nbr_Plan, PtrFAS, Norm, False);
    MR_Data.LiftingTrans = LiftingTrans;

    if (Verbose == True)
    {
        cout << "Transform = " << StringTransf1D(Transform) << endl;
        cout << "Number of scales = " << Nbr_Plan  << endl;

       if ((Transform == TO1_MALLAT) || (Transform == TU1_MALLAT) ||(Transform ==  WP1_MALLAT))
       {
          cout << StringSBFilter(SB_Filter) << endl;
          if (Norm == NORM_L2) cout << "L2 normalization" << endl;
       }
       if ((Transform == TO1_LIFTING) || (Transform == WP1_LIFTING)) 
          cout << "Lifting method = " <<  StringLSTransform(LiftingTrans) << endl;
       if (Transform == TU1_UNDECIMATED_NON_ORTHO)  cout << "Undec. Filter Bank = " << StringUndecFilter(U_Filter) << endl;
    }

    for (int s =0; s < MR_Data.nbr_band(); s++)
    {
       int Pos = MR_Data.pos_band(s);
       for (int i =0; i < MR_Data.size_scale_np(s); i++)
       {
          if ((MR_Data.Set_Transform == TRANS1_MALLAT) ||
             (MR_Data.Set_Transform == TRANS1_WP_MALLAT))
          {
	     MR_Data(s,i) = Dat(i+Pos);
	  }
	  else MR_Data(s,i) = Dat(i,s);
       }
    }
    fltarray Rec(Np);
    MR_Data.recons (Rec);
    // fits_write_fltarr (Name_Imag_Out, Rec);
    io_1d_write_data(Name_Imag_Out, Rec);
   exit(0);
} 

