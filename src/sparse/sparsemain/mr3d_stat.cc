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
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
#include "DefFunc.h"

/****************************************************************************/

char Name_Cube_In[256];
char Name_Out[256];
char NameSurv[256];

int Nbr_Plan=4;
type_trans_3d Transform=TO3_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
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
    for (i = 0; i < NBR_TRANS_3D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf3D((type_trans_3d)i));
    fprintf(OUTMAN, "              default is %s\n", StringTransf3D((type_trans_3d) Transform));
    
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();

    sb_usage(SB_Filter);
    manline();
    
    lifting_usage(LiftingTrans);
    manline();
    
    fprintf(OUTMAN, "         [-S Simu File Name]\n");
    fprintf(OUTMAN, "              Sigma of no strutured simulation for all scale. \n");    
    manline();

    nsigma_usage(N_Sigma);
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
    while ((c = GetOpt(argc,argv,"VNs:m:t:n:S:Ll:T:vzZ:")) != -1) 
    {
	switch (c) 
        {
	  case 'V': Survival = (Survival == True) ? False: True; break;
          case 'N': Normalize = (Normalize == True) ? False: True; break;
          case 's': 
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)  N_Sigma = 3;
		break;
          case 'm': 
 		if (sscanf(OptArg,"%d",&NbrCoeffMax ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of coeff: %s\n", OptArg);
	            exit(-1);
 		}
		if (c <= 0)
                {
		    fprintf(OUTMAN, "Error: bad number of coeff: %s\n", OptArg);
	            exit(-1);
 		}	              
                GetMax = True;
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
		// SB_Filter = (type_sb_filter) (c);
		break;
           case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANS_3D)) 
                                        Transform = (type_trans_3d) (c-1);
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
      
 
        if ((Transform != TO3_MALLAT) && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -f and -L are only valid with transforms 14 and 16  ... \n");
           exit(0);
	}
 	if ((Transform != TO3_LIFTING) && ( LiftingTrans != DEF_LIFT))
	{
	   fprintf(OUTMAN, "Error: option -f and -L are only valid with transforms 14 and 16  ... \n");
           exit(0);
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
       cout << "Transform = " << StringTransf3D((type_trans_3d) Transform) << endl;
       if (Transform == TO3_MALLAT)  
       {
         cout << "Filter = " <<  StringSBFilter(SB_Filter) << endl;
         if (Norm == NORM_L1) cout << "Norm L1 " << endl;
         else cout << "Norm L2 " << endl;
       }
       if (Transform == TO3_LIFTING)  
         cout << "Lifting method = " <<  StringLSTransform(LiftingTrans) << endl;
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
    
    MR_3D MR_Data;
    MR_Data.Verbose=Verbose;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if (Transform == TO3_MALLAT) 
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_Data.alloc (Nx, Ny, Nz, Transform, Nbr_Plan, PtrFAS, Norm);
    int NbrBand = MR_Data.nbr_band();
    
    if (Transform == TO3_LIFTING)  MR_Data.LiftingTrans = LiftingTrans;
        
    if (Verbose == True)
    {
        cout << "Number of scales = " << MR_Data.nbr_scale() << endl;
	// cout << "Number of bands = " << NbrBand << endl;
    }
    if (Verbose == True) MR_Data.info_pos_band();
    
    if (Verbose == True) cout << "3D transform ... " << endl << endl;
    MR_Data.transform (Dat);
    // if (Verbose == True) cout << "Write ... " << endl;
    // io_3d_write_data(Name_Cube_Out, MR_Data.cube());
    // MR_Data.write(Name_Cube_Out);
     	
    fltarray Simu;
    fltarray StaInfo;
    int NbrStat = 8;
    if (WithFileSimu) 
    {
        NbrStat++;
	StaInfo.alloc(NbrBand,NbrStat);
	fits_read_fltarr(FileSimu, Simu);
	if (Simu.nx() !=  NbrBand)
	{
	    cout << "Error: simulation file has not the correct size ... " << endl;
	    cout << "       it should have " << NbrBand << " entries " << endl;
	    cout << "       Nx = " << Simu.nx() << endl;
	    exit(-1);
	}
    }
    else StaInfo.alloc(NbrBand,NbrStat);
    fltarray MaxVal(NbrBand, NbrCoeffMax);
	
    fltarray Band;
    fltarray TabSurvival;
    for (b=0;b<NbrBand;b++) 
    {
	double Mean, Sigma, Skew, Curt;
        float  Min, Max;
	fltarray ResSurvival;
	
	MR_Data.get_band(b, Band);
	N = Band.n_elem();
        if (Survival == True)
	{
	    // void survival (fltarray & Data, fltarray &Survival, int NpixSurv=1001, 
	    //               float NSigma=10., Bool Norm=True, float SigmaNorm=0., float MeanNorm=-1.);
            float NSig=50.;
	    int Np=4001;
	    survival (Band, ResSurvival, Np, NSig);
	    if (b == 0) TabSurvival.alloc(ResSurvival.nx(), 2, NbrBand);
	    for (int s=0;s< ResSurvival.nx();s++)
	    {
	       TabSurvival(s, 0, b) = ResSurvival(s,0);
	       TabSurvival(s, 1, b) = ResSurvival(s,1);
	    }
	}
        // if (Verbose == True) cout << "BAND " << b+1 << " Nx = " << Band.nx() << " Ny = " << Band.ny() <<  " Nz = " << Band.nz() << endl;
	float SigmaSimu=0.;
	if (WithFileSimu == True) 
	{
	      SigmaSimu = Simu(b);
	      if (Verbose == True) cout << "SigmaSimu = " << SigmaSimu << endl;
	}
        moment4(Band.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);

        StaInfo(b,0) = (float) Mean;
	StaInfo(b,1) = (float) Sigma;
	StaInfo(b,2) = (float) Skew;
	StaInfo(b,3) = (float) Curt;
	StaInfo(b,4) = Min;
	StaInfo(b,5) = Max;
        StaInfo(b,6) = Band.energy();
        StaInfo(b,7)=0;
	
	for (int i=0;i< Band.n_elem();i++) 
	{
	   float Coef = Band(i);
	   if (fabs(Coef- StaInfo(b,0)) > N_Sigma*StaInfo(b,1))
	        StaInfo(b,7) += Coef*Coef;
        }

        if (WithFileSimu) 
	{
	    StaInfo(b,8)=0;
	    for (int i=0;i< Band.n_elem();i++) 
	    {
	       float Coef = Band(i);
	       if (fabs(Coef - StaInfo(b,0)) > N_Sigma*SigmaSimu)
	            StaInfo(b,8) += Coef*Coef;
	    }
        }	
	   
        //!!!!!!!!!!!!!!! modify the Band values !!!!!!!!!!!!!!!!!!!!!!!
	if (GetMax == True)
	{
	   for (int i=0;i< NbrCoeffMax;i++) 
	   {
	       int IndMax=0;
	       MaxVal(b,i) = Band.maxfabs(IndMax);
	       Band(IndMax)=0.0;
	   }
        }	 
        if (Verbose == True)
	{
            printf("  Band %d: Nx = %d, Ny = %d, Nz = %d,  Min = %5.3f, Max = %5.3f\n",
	           b+1, Band.nx(), Band.ny(), Band.nz(),  StaInfo(b, 4), StaInfo(b, 5));
            printf("           Mean = %5.3f, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Ener3sig = %5.3f\n",
	            StaInfo(b, 0), StaInfo(b,1),  StaInfo(b, 2),StaInfo(b,3), StaInfo(b,7));
	    if (WithFileSimu == True) printf("           Ener3sigNoise = %5.3f\n",  StaInfo(b,8));
        }
    }
    Header.hd_fltarray(StaInfo);
    Header.origin = Cmd;	 
    fits_write_fltarr (Name_Out, StaInfo, &Header);
    sprintf(NameSurv, "%s_surv.fits", Name_Out);
    if (Survival == True) fits_write_fltarr (NameSurv, TabSurvival);
    
    if (GetMax == True) fits_write_fltarr ("MaxWavelet.fits", MaxVal);
    exit(0);
}
