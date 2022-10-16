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
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
#include "DefFunc.h"

/****************************************************************************/

char Name_Cube_In[256];
char Name_Out[256];
char NameSurv[256];

int Nbr_Plan=9;
type_trans_3d Transform=TO3_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

sb_type_norm Norm = NORM_L2;
type_sb_filter SB_Filter = F_HAAR; // F_MALLAT_7_9;
type_lift LiftingTrans = DEF_LIFT;
Bool Verbose=False;
Bool InfoBand=False;
char FileSimu[180];
Bool WithFileSimu=False;
int NbrCoeffMax = 10;
Bool GetMax = False;
float N_Sigma=3;                // number of sigma (for the noise)  
Bool Normalize=False;       // normalize data in
Bool Survival=False;        


type_border Bord = I_MIRROR;

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

//     fprintf(OUTMAN, "\n");
//     fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
//     for (i = 0; i < NBR_TRANS_3D; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i+1,
//                                             StringTransf3D((type_trans_3d)i));
//     fprintf(OUTMAN, "              default is %s\n", StringTransf3D((type_trans_3d) Transform));
    
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();

    sb_usage(SB_Filter);
    manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Normalize the data. Default is no. \n");    
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
    Bool OptL = False, Optf = False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
    while ((c = GetOpt(argc,argv,"VNs:m:n:S:Ll:T:vzZ:")) != -1) 
    {
	switch (c) 
        {
	  case 'V': Survival = (Survival == True) ? False: True; break;
          case 'N': Normalize = (Normalize == True) ? False: True; break;
 	   case 'v': Verbose = True; break;
	   case 'L': Norm = NORM_L1; OptL = True; break;
	   case 'T':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad filter: %s\n", OptArg);
	            exit(-1);
                    
		}
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		// SB_Filter = (type_sb_filter) (c);
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
	   case 'S':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%s",FileSimu) != 1) 
                {
		    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		    exit(-1);
		}
		WithFileSimu = True;
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
    int N,b;
    /* Get command line arguments, open input file(s) if necessary */
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   
    lm_check(LIC_MR3);
    transinit(argc, argv);
 
    if (Verbose == True)
    {
       cout << "Filename in = " << Name_Cube_In << endl;
       cout << "Filename out = " << Name_Out  << endl;
       // cout << "Transform = " << StringTransf3D((type_trans_3d) Transform) << endl;
    }
    io_3d_read_data(Name_Cube_In, Dat);
 
    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();
    if (Verbose == True)
    {
        cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
    }
    if (Normalize == True)
    {
      double Mean = Dat.mean();
      double Sigma = Dat.sigma();
      for (int i=0;i<Nx;i++)
      for (int j=0;j<Ny;j++)
      for (int k=0;k<Nz;k++) Dat(i,j,k) = (Dat(i,j,k)-Mean)/Sigma;
    }    

   // NAXIS   =                    3 /                                                
   // NAXIS1  =                  255 /                                                
   // NAXIS2  =                  134 /                                                
   // NAXIS3  =                   32 /       
    int NbrScale2d = 5;
    fltarray *TabWT;
    TabWT= new fltarray [NbrScale2d];
    
    if (Verbose == True) cout << "Alloc 2D ... " << Nx << endl;
    int NbrFrame = Nx;
    int Nl = Nz;
    int Nc = Ny;
    for (b=0; b < NbrScale2d; b++) (TabWT[b]).alloc(Nc, Nl, NbrFrame);   //  = alloc(Ny, Nz, Nx)
    
    if (Verbose == True) cout << "Alloc Class AT ... " << Nx << endl;
    ATROUS_2D_WT AT;
    Ifloat *TabTrans;
    AT.alloc(TabTrans, Nl,  Nc, NbrScale2d);
    AT.Bord = I_ZERO;
    
   
    if (Verbose == True) cout << "Transform 2D ... " << NbrScale2d << endl;
    Ifloat Frame(Nl, Nc);
    int i,j,z;
    for (z=0; z < NbrFrame; z++)
    {
       for (i =0; i < Nl; i++)
       for (j =0; j < Nc; j++) Frame(i,j) = Dat(z,j,i);
       AT.transform(Frame, TabTrans, NbrScale2d);
       for (b=0; b < NbrScale2d; b++) 
       for (i =0; i < Nl; i++)
       for (j =0; j < Nc; j++) (TabWT[b])(j,i,z) = (TabTrans[b])(i,j);
    }
    // io_3d_write_data("xx_2d.fits", Trans3D);   
    
    if (Verbose == True) cout << " Transform 1D ... " << Nbr_Plan << endl;
    FilterAnaSynt FAS;
    // FAS.Verbose = Verbose;
    FAS.alloc(SB_Filter);
    SubBandFilter SBF(FAS,  Norm);
    SBF.Border = Bord;
    if (Verbose == True) cout << StringSBFilter(SB_Filter) << endl;
    MR_1D MR_Data;
    Bool Rebin=False;
    MR_Data.alloc (NbrFrame, TO1_MALLAT, Nbr_Plan , &FAS, Norm, Rebin);
    intarray TabFrom(Nbr_Plan);
    intarray TabNpix(Nbr_Plan);
    intarray TabCoefPerBand(NbrScale2d,Nbr_Plan);
    for (b=0; b < MR_Data.nbr_band(); b++)
    {
        int NpixBand = MR_Data.size_scale_np(b);
	if (b == 0) TabFrom(b) = 0;
	else TabFrom(b) = TabFrom(b-1) + TabNpix(b-1);
	TabNpix(b) = NpixBand;
	if (Verbose == True) cout << " Band " << b+1 << " From = " << TabFrom(b) << " Npix = " << TabNpix(b) << endl;
    }
    
    if (Verbose == True) cout << "Transform 1D: " << Nbr_Plan << " NbrFrame = " << NbrFrame << endl;
    fltarray  Vect(NbrFrame);
    fltarray  TransVect(NbrFrame);
    char Name[256];
    for (int s2d=0; s2d < NbrScale2d; s2d++)
    {
      for (i =0; i < Nl; i++)
      for (j =0; j < Nc; j++) 
      {
         for (z=0; z < NbrFrame; z++) Vect(z) = (TabWT[s2d])(j,i,z);
         MR_Data.transform(Vect);
         for (b=0; b < MR_Data.nbr_band(); b++)
         {
           for (z=0;z < MR_Data.size_scale_np(b); z++) 
	   {
	      (TabWT[s2d])(j,i, TabFrom(b)+z) = MR_Data(b,z);
	       TabCoefPerBand(s2d,b) ++;
	   }
         }
      }
      sprintf(Name, "wscale%d.fits", s2d+1);
      io_3d_write_data(Name, TabWT[s2d]);
    }
    
    if (Verbose == True) cout << "Compute the statistics ... " << endl;
    int BordSize=0;
    double Mean, Sigma, Skew, Curt;
    float  Min, Max; 
    fltarray ResSurvival;
    fltarray TabSurvival;
    fltarray TabSurvivalX;	
    fltarray StaInfo;
    int NbrStat = 6;
    StaInfo.alloc(NbrScale2d, MR_Data.nbr_band(), NbrStat);
    Bool Survival=True;
    
    for (int s2d=0; s2d < NbrScale2d; s2d++)
    {
       for (b=0; b < MR_Data.nbr_band(); b++)
       {
          int From = TabFrom(b);
	  int Np = TabNpix(b);
	  if (Verbose == True) cout << "Band 2D " << s2d+1 << " Band 1D " << b+1 << " NbrCoef = " <<  TabCoefPerBand(s2d,b) << endl;
	  fltarray Band(TabCoefPerBand(s2d,b));
	  int ind=0;
	  for (i = BordSize; i < Nl-BordSize; i++)
          for (j = BordSize; j < Nc-BordSize; j++) 
	  for (z=0;z < Np; z++)   Band(ind++) = (TabWT[s2d])(j,i, From+z);
 	  
          int N =  ind;
          if (Survival == True)
	  {
            float NSig=50.;
	    int Np=4001;
	    survival (Band, ResSurvival, Np, NSig);
	    if ((s2d==0) && (b == 0)) TabSurvival.alloc(ResSurvival.nx(), MR_Data.nbr_band(), NbrScale2d);
	    if ((s2d==0) && (b == 0)) TabSurvivalX.alloc(ResSurvival.nx(), MR_Data.nbr_band(), NbrScale2d);
	    for (int s=0;s< ResSurvival.nx();s++)
	    {
	       TabSurvivalX(s, b, s2d) = ResSurvival(s,0);
	       TabSurvival(s, b, s2d) = ResSurvival(s,1);
	    }
	  }
	   if (Verbose == True)cout << "moment " << endl;
          moment4(Band.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
	  
          StaInfo(s2d, b,0) = (float) Mean;
	  StaInfo(s2d, b,1) = (float) Sigma;
	  StaInfo(s2d, b,2) = (float) Skew;
	  StaInfo(s2d, b,3) = (float) Curt;
	  StaInfo(s2d, b,4) = Min;
	  StaInfo(s2d, b,5) = Max;
	 
          if (Verbose == True)
	  {
            printf("  Band %d %d:  Min = %5.3f, Max = %5.3f\n",
	           s2d+1, b+1,  StaInfo(s2d, b, 4), StaInfo(s2d, b, 5));
            printf("           Mean = %5.3f, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f \n",
	            StaInfo(s2d, b, 0), StaInfo(s2d, b,1),  StaInfo(s2d, b, 2),StaInfo(s2d, b,3));
          }
       }
       BordSize=BordSize*2;
    }

     if (Verbose == True)cout << "Write result ...  " << endl;
    Header.hd_fltarray(StaInfo);
    Header.origin = Cmd;	 
    fits_write_fltarr (Name_Out, StaInfo, &Header);
    sprintf(NameSurv, "%s_surv.fits", Name_Out);
    if (Survival == True) fits_write_fltarr (NameSurv, TabSurvival);
    sprintf(NameSurv, "%s_survX.fits", Name_Out);
    if (Survival == True) fits_write_fltarr (NameSurv, TabSurvivalX);
    
    exit(0);
}
