/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/03
**    
**    File:  mbase.cc
**
*******************************************************************************
**
**    DESCRIPTION  Decomposition of an image on multiple bases
**    ----------- 
**                 
******************************************************************************/


#include "CB_MCA.h" 
#include "MR_SoftRegul.h" 

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale1D = -1;
int NbrScale2D = 5;
float Noise_Ima=0.5;
Bool PoissonNoise=False;

Bool PositivSol=True;
int BlockSize=0;
int CosBlockSize=0;
int Nbr_Iter = DEF_MCA_NBR_ITER;
Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
int FirstDetectScale=0;
Bool KillLastScale = False;
 
Bool TabBase[NBR_TOT_MCA];
type_mcabase TabSelect[NBR_TOT_MCA];
int NbrBase = 0;

int WT_NbrUndecimatedScale = -1;
int WP_NbrUndecimatedScale = 2;
Bool WriteAllRec=True;
float FirstSoftThreshold = -1.; // DEF_MCA_FIRST_SOFT_THRESHOLD;
float LastSoftThreshold = DEF_MCA_LAST_SOFT_THRESHOLD;
Bool RidKillLastScale = False;
Bool TV = False;
Bool UsePsf = False;
float LambdaTV= 0.;
Bool  UseNormL1 = False;
float SoftLevel = 1.5;
int DWT_NbrAngle = 10;
Bool TestMode = False;
float COS_Sensibility = 1.;
float COS_Min  = 0.;
Bool RemoveLastScale = False;
Bool COS_WeightFirst  =False;
Bool CosIsotrop =False;
Bool UseMad =False;
Bool UseMask =False;
Bool UseZoom=False;
Bool WriteAll=False;
Bool Bounded=False;
Bool SigmaBounded=False;
int TV_NbrScale2D=2;
Bool Linear = True;
float CorrelConstraint=0.;

int NDir=32;
int ZoomFactor=4;
type_threshold_decrease T_Decreasing = MCA_THRES_MOM;  // MCA_THRES_LINEAR;

/***********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t TransformSelection]\n");
    for (int i = 0; i < NBR_MCA; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMCABase ((type_mcabase) i));
    manline();
    
		fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "             Min-Max Bounded\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();
    
		fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "             Sigma Bounded\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();
		
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the WT, the a trous, the PMT and the curvelet transform.\n");
    fprintf(OUTMAN, "             default is 5.\n");    
    manline();
 
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size in the ridgelet transform.\n");
    fprintf(OUTMAN, "             Default is image size. \n");    
    manline();

    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration. Default is %d.\n", Nbr_Iter);    
    manline();
 
    fprintf(OUTMAN, "         [-B DCT_BlockSize]\n");
    fprintf(OUTMAN, "             Local DCT block size.\n");
    fprintf(OUTMAN, "             By default, a local DCT (block size = 32) is used. \n");    
    manline();

    fprintf(OUTMAN, "         [-S FirstThresholdLevel]\n");
    fprintf(OUTMAN, "            First thresholding value.\n");
    fprintf(OUTMAN, "            Default is derived from the data.\n");    
    manline();

    fprintf(OUTMAN, "         [-s LastThresholdLevel]\n");
    fprintf(OUTMAN, "            Last thresholding value..\n");
    fprintf(OUTMAN, "            default is %f.\n", LastSoftThreshold);    
    manline();

//     fprintf(OUTMAN, "         [-u]\n");
//     fprintf(OUTMAN, "             Number of undecimated scales in the WT.\n");
//     fprintf(OUTMAN, "             default is all scales. \n");    
//     manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Minimize the L1 norm.\n");
    fprintf(OUTMAN, "             Default is L0 norm.\n");    
    manline();

    fprintf(OUTMAN, "         [-L Decrease]\n");
    fprintf(OUTMAN, "             1- Linear descent\n");
    fprintf(OUTMAN, "             2- Exponential descent\n");
    fprintf(OUTMAN, "             3- MOM descent\n");  
    fprintf(OUTMAN, "             Default is linear descent.\n");  
    manline();
    
    fprintf(OUTMAN, "         [-l]\n");
    fprintf(OUTMAN, "             Remove last scale. Default is no. \n");    
    manline();
     
    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation. \n");
    fprintf(OUTMAN, "             Default is 0.5 (quantification noise).\n");
    fprintf(OUTMAN, "             if sigma is set to 0, noise standard deviation\n");
    fprintf(OUTMAN, "             is automatically estimated.\n");
    manline();  

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Supress the block overlapping. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Supress the positivity constraint. Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-G RegulVal[,NbrScale]]\n");
    fprintf(OUTMAN, "             Total Variation regularization term. Default is 0.\n");    
    fprintf(OUTMAN, "             NbrScale = number of scales used in Haar TV regularizarion. \n");
    fprintf(OUTMAN, "                        default is 2. \n");   
    manline();
    
    fprintf(OUTMAN, "         [-H]\n");
    fprintf(OUTMAN, "             Data contained masked area (must have a zero value). Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-I]\n");
    fprintf(OUTMAN, "             Interpolate the data (super-resolution). Default is no. \n");    
    manline();
    vm_usage();
    manline();    
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}

/***********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
   // float Val;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"c:ldDHMNIWL:C:B:Xa:PG:s:S:wu:U:t:KOi:b:g:n:vzZ")) != -1) 
    {
	switch (c) 
        {  
   case 'd': Bounded = True; break;	   
	 case 'D': SigmaBounded = True; break;			  	 
	 case 'w': WriteAll = True; break;
	 case 'L': if (sscanf(OptArg,"%d",&c) != 1) 
                {
                fprintf(OUTMAN, "Error: bad type of decrease %s\n", OptArg);
                exit(-1);
                }
                if ((c > 0) && (c <= NBR_THRES_DEC))
                  {
                   T_Decreasing = (type_threshold_decrease)(c-1);
                  }
                else  
                  {
                   fprintf(OUTMAN, "Error: bad type of decrease: %s\n", OptArg);
                   exit(-1);
                  }
                break;
	 case 'M': UseMad = True; break;
	 case 'H': UseMask = True; break;
	 case 'I': UseZoom = True;
	           // CosIsotrop = True; 
	           break;
	 case 'l': RemoveLastScale = (RemoveLastScale == True) ? False:True; 
	           PositivSol = False; break;
	 case 'W': COS_WeightFirst   = (COS_WeightFirst   == True) ? False: True; break;
	 case 'C':
                if (sscanf(OptArg,"%f",&COS_Min ) != 1) 
                {
                    fprintf(OUTMAN, "bad cosinus sensibility: %s\n", OptArg);
                    exit(-1);
                }
                break;
	 case 'c':
                if (sscanf(OptArg,"%f",&CorrelConstraint ) != 1) 
                {
                    fprintf(OUTMAN, "bad correlation constraint parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	 case 'X':TestMode = True; break;
         case 'a':
                if (sscanf(OptArg,"%d",&DWT_NbrAngle) != 1) 
                {
                    fprintf(OUTMAN, "bad number of angles: %s\n", OptArg);
                    exit(-1);
                }
                break;
          case 'N': 
	        UseNormL1 = (UseNormL1  == True) ? False: True;
		break;
  	  case 'G':   
	       {
	           float LT=0;
		   int NS=0;
	           int N = sscanf(OptArg,"%f,%d",&LT, &NS);
                   if (N < 1)
                   {
                      fprintf(OUTMAN, "bad regul param: %s\n",  OptArg);
                      exit(-1);
                   }
		   else
		   {
		      LambdaTV = LT;
		      if (LambdaTV < 0)
	              {
			    fprintf(OUTMAN, "Error: bad TV regul param: %s\n", OptArg);
			    exit(-1);
		      }
 	              TV = True; 
		      if (N > 1) 
                      {
		         TV_NbrScale2D = NS;
			 if (TV_NbrScale2D < 2)
			 {
			    fprintf(OUTMAN, "Error: bad number of TV scales: %s\n", OptArg);
			    exit(-1);
			 }
		      }
		   }
		}
		break;
          case 'k': RidKillLastScale = True; break;
          case 'S':               
                if (sscanf(OptArg,"%f",&FirstSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0)
                {
                   fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                   exit(-1);
                }
                break;
          case 's':               
                if (sscanf(OptArg,"%f",&LastSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
                    exit(-1);
                }
//                 if (FirstSoftThreshold < 0)
//                 {
//                    fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
//                    exit(-1);
//                 } 
                break;
          case 'P': PositivSol = False; break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_TOT_MCA))
                {
                   if (TabBase[c-1] == True)
                   {
                       fprintf(OUTMAN, "Error: transform already selected: %s\n", OptArg);
                       exit(-1);
                   }
                   else
                   {
                      TabBase[c-1] = True;
                      TabSelect[NbrBase++] = (type_mcabase) (c-1);
                   }
                }
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'K': KillLastScale = True; break;
           case 'u': 
                if (sscanf(OptArg,"%d",&WT_NbrUndecimatedScale) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                    exit(-1);
                }
                if (WT_NbrUndecimatedScale < 0)  
                {
                   fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                   exit(-1);
                }
                break;
            case 'U': 
                if (sscanf(OptArg,"%d",&WP_NbrUndecimatedScale) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                    exit(-1);
                }
                if (WP_NbrUndecimatedScale < 0)  
                {
                   fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                   exit(-1);
                }
                break;
	   case 'F': 
                if (sscanf(OptArg,"%d",&FirstDetectScale ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first detection scale: %s\n", OptArg);
                    exit(-1);
                }
                FirstDetectScale --;
                if (FirstDetectScale < 0)  
                {
                   fprintf(OUTMAN, "Error: bad first detection scale: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
           case 'i':
                if (sscanf(OptArg,"%d",&Nbr_Iter) != 1) 
                {
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
//            case 'p': PoissonNoise=True;
//                      Noise_Ima = 1.;
//                      break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'B':
                if (sscanf(OptArg,"%d",&CosBlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad cosinus transform block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
     case 'v': Verbose = True;break;
	   case 'n':
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D  > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
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
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/***************************************************************************/

static void abs_morpho_dilation (Ifloat &Imag1, Ifloat &Imag2, int Window_Size)
{
    int i,j,k,l;
    float Max;
    int Nl = Imag1.nl();
    int Nc = Imag1.nc();
    int Window2 = (Window_Size - 1) / 2;

    for (i = 0; i < Nl; i++) 
    for (j = 0; j < Nc; j++) 
    { 
        Max =  Imag1(i,j);
        for (k = i - Window2; k <= i + Window2; k++)
        for (l = j - Window2; l <= j + Window2; l++)
            if (ABS(Imag1(k,l, I_CONT)) > ABS(Max)) Max = Imag1(k,l, I_CONT);
        Imag2(i,j) = Max;
    }
}


void first_guess_interpol(Ifloat &Data, int NbrMissing)
{
   Ifloat Buff_Dilat(Data.nl(), Data.nc());
   int Nl = Data.nl();
   int Nc = Data.nc();
   int Nm = NbrMissing;
   while (Nm > 0)
   {
       abs_morpho_dilation (Data, Buff_Dilat, 3);
       Nm = 0;
       for (int i=0; i < Nl; i++)
       for (int j=0; j < Nc; j++) 
       {
          if (Data(i,j) == 0) Data(i,j) = Buff_Dilat(i,j);
	  if (Data(i,j) == 0) Nm ++;
       }
       // cout << Nm << endl;
       io_write_ima_float("xx_temp.fits", Data);
   }
}


/***************************************************************************/

int main(int argc, char *argv[])
{
    int i,j,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR3);
    for (i = 0; i < NBR_TOT_MCA; i++) TabBase[i] = False;
    filtinit(argc, argv);
    if (TestMode == True)
    {
     Ifloat Data, DataTrans, Rec;
     type_sb_filter SB_Filter = F_MALLAT_7_9;
     FilterAnaSynt *WT_SelectFilter; // Selected filter bank
     SubBandFilter *WT_SB1D;  // Pointer to the fiter bank decomposition class
     WT_SelectFilter = new FilterAnaSynt(SB_Filter);
     WT_SB1D = new SubBandFilter(*WT_SelectFilter, NORM_L2);
     DirectionalLineCol LC(*WT_SB1D);
     // float AngleRot = Noise_Ima;
     //Bool RadianUnitAngle=False;
     int NbrScaleLine = 0;
     int NbrScaleCol = 5;
     LC.NbrUndecimatedScaleCol = 2;
     io_read_ima_float(Name_Imag_In, Data, &Header);
     cout << "trans " << NbrScaleLine << " " <<  NbrScaleCol  << endl;        
     MDCT M_DCT;
     int FirstBlockSize=8;
     M_DCT.Verbose = True;
     M_DCT.BlockOverlap = False;
     M_DCT.alloc(Data.nl(), Data.nc(), NbrScaleCol, FirstBlockSize);
     M_DCT.transform(Data);
     M_DCT.threshold(Noise_Ima, 3.);
     Data.init();  
     M_DCT.recons(Data);
//      LC.transform(Data, DataTrans, NbrScaleLine, NbrScaleCol);
//      // LC.directional_transform(Data, DataTrans, AngleRot, NbrScaleLine, NbrScaleCol, RadianUnitAngle);
//      cout << "end trans" << endl;
//      io_write_ima_float("xx_trans", DataTrans);
//      Data.init();
//      cout << "rec" << NbrScaleLine << " " <<  NbrScaleCol  << endl;
//      // LC.directional_recons(DataTrans, Data, AngleRot,NbrScaleLine, NbrScaleCol, RadianUnitAngle);
//      LC.recons(DataTrans, Data, NbrScaleLine, NbrScaleCol);

     cout << "end rec" << endl;
     io_write_ima_float("xx_rec", Data);
     exit(-1);
    }     
    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
	if (CosBlockSize > 0)  cout << "CosBlockSize = " << CosBlockSize <<endl;
        cout << "NbrScale2D = " << NbrScale2D << endl;
        if (FirstDetectScale > 0) cout << "FirstDetectScale = " << FirstDetectScale << endl;
        if (FirstSoftThreshold > 0)
        cout << "FirstThreshold = " << FirstSoftThreshold<< endl;
       cout << "LastThreshold = " << LastSoftThreshold<< endl;
       if (Noise_Ima > FLOAT_EPSILON) cout << " Noise = " << Noise_Ima << endl;
       if (RemoveLastScale == True) cout << "Remove last scale " << endl;    
       cout << "Number of iteration = " << Nbr_Iter << endl; 
       if (T_Decreasing == MCA_THRES_LINEAR)  cout << "Linear step descent" << endl;
       else if (T_Decreasing == MCA_THRES_EXP)  cout << "Exponential step descent" << endl;
       else if (T_Decreasing == MCA_THRES_MOM) cout << "MOM step descent" << endl;
       else cout << "Exponential step descent" << endl;  
       if (UseNormL1)   cout << "Used L1 norm in threshold" << endl;
       if (TV == True) cout << "TV Regularization: LambdaTV = " << LambdaTV << " NbrHaarScale = " << TV_NbrScale2D << endl;
       if (UseZoom == True) cout << "Zoom the image " << endl;
			 if (Bounded == True) cout << "Bounded " << endl;
			 if (SigmaBounded == True) cout << "SigmaBounded " << endl;
    }
    if (NbrBase == 0)
    {
       NbrBase = 2;
       TabSelect[1] = MCA_CUR;
       TabSelect[0] = MCA_ATROU;
    }
    Ifloat Tab[NbrBase];
    MCA MB;  // Multiple base decomposition Class
    Ifloat Data;   
    io_read_ima_float(Name_Imag_In, Data, &Header);       
    if (UseZoom == True)
    {
        Ifloat Data1( Data.nl()* ZoomFactor, Data.nc()*ZoomFactor, "ZZ");
        MB.UnZoomData.resize(Data.nl(), Data.nc());
	      MB.UnZoomData = Data;
	      im_zoom(Data, Data1);
        Data.resize(Data1.nl(), Data1.nc());
	      Data = Data1;
    }
    Header.origin = Cmd;
    int Nl = Data.nl();
    int Nc = Data.nc();
   
    if (Noise_Ima < FLOAT_EPSILON) 
    {
       Noise_Ima = detect_noise_from_med (Data);
       if (Verbose == True) cout << " Estimated Noise = " << Noise_Ima << endl;
    }
//     if (FirstSoftThreshold < LastSoftThreshold)
//     {
//        FirstSoftThreshold = 0.8*(Data.max() - Data.min()) / Noise_Ima;
//        if (Verbose == True) cout << " FirstThreshold = " << FirstSoftThreshold<< endl;
//     }
    for (i=0; i < NbrBase; i++) Tab[i].alloc(Nl, Nc, "alloc");

    // for (i = 0; i < NBR_MCA; i++) MB.TabBase[i] = TabBase[i];
    MB.NbrBase = NbrBase;
    for (i = 0; i < NbrBase; i++) MB.TabSelect[i] = TabSelect[i];
    MB.DataSigmaNoise = Noise_Ima;
    // MB.PoissonNoise = PoissonNoise;
    MB.Verbose = Verbose;
		MB.Bounded = Bounded;
		MB.SigmaBounded = SigmaBounded;
    MB.Nbr_Iter = Nbr_Iter;
    MB.PositivRecIma = PositivSol;
    MB.WT_PositivRecIma= PositivSol;
    MB.CUR_PositivRecIma= PositivSol;
    MB.RID_PositivRecIma = PositivSol;
    MB.AT_PositivRecIma = PositivSol;
    MB.FCUR_PositivRecIma = PositivSol;
    MB.PCUR_PositivRecIma = PositivSol;
    
    MB.FirstSoftThreshold = FirstSoftThreshold;
    MB.LastSoftThreshold = LastSoftThreshold;
    MB.UseNormL1 = UseNormL1;
    MB.NSigmaSoft = SoftLevel;
    
    MB.RID.RidTrans = RidTrans;
    MB.RID.BlockOverlap = BlockOverlap;
    MB.RID_FirstDetectScale = FirstDetectScale;
    // MB.RID.KillLastScale = RidKillLastScale;
    if (BlockSize > 0) MB.RID_BlockSize = BlockSize;

    MB.AT_NbrScale2D = NbrScale2D;
    MB.AT_KillLastScale  = KillLastScale;

    MB.WT_NbrScale2D = NbrScale2D;
    MB.WT_NbrUndecimatedScale = WT_NbrUndecimatedScale;
 
    MB.CUR_NbrScale2D = NbrScale2D;
    // MB.CUR_BlockSize = BlockSize;
    MB.CUR_BlockOverlap = BlockOverlap;
    MB.CUR_KillLastScale  = KillLastScale;
    
    MB.PCUR_NbrScale2D = NbrScale2D;
    MB.PCUR_BlockOverlap = BlockOverlap;
    MB.PCUR_KillLastScale  = KillLastScale;

    MB.FCUR_NbrScale2D = NbrScale2D;
    // MB.FCUR_NDIR = BlockSize;
    MB.FCUR_KillLastScale  = KillLastScale;
    
    MB.TotalVariation = TV;
    MB.LambdaTV = LambdaTV;
    
    MB.COS_PositivRecIma = False;
    if (CosBlockSize > 0) MB.COS_BlockSize = CosBlockSize;
    MB.COS_Overlapping = BlockOverlap;
    MB.COS_Sensibility = COS_Sensibility;
    MB.COS_Min = COS_Min ;
    MB.COS_WeightFirst = COS_WeightFirst;
    MB.COS_Isotrop = CosIsotrop;
    
    MB.WP_NbrUndecimatedScale = WP_NbrUndecimatedScale;
    MB.WP_NbrScale2D = NbrScale2D;
    MB.WP_PositivRecIma = False;
    MB.WP_Filter = F_MALLAT_7_9;
    if (COS_Min > 0) MB.WP_KillLastScale = True;
    else MB.WP_KillLastScale = False;
    
    MB.MCOS_Overlapping = BlockOverlap;
    MB.MCOS_PositivRecIma = False;
    if (CosBlockSize > 0) MB.MCOS_FirstBlockSize = CosBlockSize;
    MB.MCOS_NbrScale2D = NbrScale2D;
    
    MB.DWT_NbrScale2Di = NbrScale2D;
    MB.DWT_NbrScale2Dj = 2;
    MB.DWT_PositivRecIma = PositivSol;
    MB.DWT_NbrAngle = DWT_NbrAngle;

    MB.RemoveLastScale = RemoveLastScale;
    MB.UseMad = UseMad;
    MB.UseMask = UseMask;
    if (CorrelConstraint > 0)
    {
    MB.CorrelConstraint = CorrelConstraint;
	  MB.UseCorrelConstraint = True;
    }
    MB.Verbose = Verbose;
    MB.WriteAll = WriteAll;
    MB.UseZoom = UseZoom;
    MB.ZoomFactor = ZoomFactor;
    MB.TV_NbrScale2D = TV_NbrScale2D;
    // MB.Linear = Linear;
    MB.T_Decreasing = T_Decreasing;
		MB.Bounded = Bounded;
		MB.SigmaBounded = SigmaBounded;
    if ((NbrBase == 1) && (MB.T_Decreasing == MCA_THRES_MOM)) MB.T_Decreasing = MCA_THRES_LINEAR;
    MB.alloc(Nl, Nc, Tab);
    if (UseMask == True)
    {
       MB.MaskedData.alloc(Nl, Nc, "Mask");
       if (UseZoom == True)
         {
          for (i=0; i <  Data.nl()/2; i++)   
	        for (j=0; j <  Data.nc(); j++) MB.MaskedData(2*i,j) = 1;
	        for (i=0; i <  Data.nl(); i++)   
	        for (j=0; j <  Data.nc()/2; j++) MB.MaskedData(i,2*j) = 1;
         }
       else
         {
          float Mean=0.;
	        int Npix=0;
          for (i=1; i < Nl-1; i++)
          for (j=1; j < Nc-1; j++)
          {
       if (Data(i,j) == 0) MB.MaskedData(i,j) = 0;
	     else 
	       {
	        MB.MaskedData(i,j) = 1;
		      Mean += Data(i,j);
		      Npix ++;
	        for (int l=-1; l <= 1; l++)
	        for (int m=-1; m <= 1; m++)
          if ((l!=0) && (m!=0) && (Data(i+l,j+m, I_MIRROR) == 0))  MB.MaskedData(i,j) = 0.5;
	       }
	    }
 	    Mean /= (float) Npix;
      int NbrMissing =  Nl*Nc - MB.MaskedData.total();
      first_guess_interpol(Data, NbrMissing);
 	  
      if (Verbose == True)
      cout << "Number of missing data = " << Nl*Nc - MB.MaskedData.total() << endl;
//           Ifloat Background(Nl,Nc,"bgr");
//           int NPix = 16;
//           int Ns = MIN(Nl,Nc);
//           int BGR_Nbr_Plan = 1;
//           while (NPix < Ns) 
//           {
//            BGR_Nbr_Plan ++;
//            Ns /= 2;
//          }
//          if (BGR_Nbr_Plan < 3) BGR_Nbr_Plan = 3;
//          if (Verbose == True) cout << "BGR estimation: Number of scales = " << BGR_Nbr_Plan << endl;
//          MultiResol BGR_MR_Data (Nl, Nc, BGR_Nbr_Plan, TM_TO_PYR, "MR_Transform");   // TM_PYR_MEDIAN
//          BGR_MR_Data.Border=I_MIRROR;
//          BGR_MR_Data.transform (Data);
//          for (int b=0; b < BGR_Nbr_Plan-1; b++) BGR_MR_Data.band(b).init();
//          BGR_MR_Data.recons(Background);
//          // io_write_ima_float("xx_bgr.fits", Background, &Header);
//          // INFO(Background, "BGR");
//          for (i=0; i < Nl; i++)
//          for (j=0; j < Nc; j++) if (Data(i,j) == 0) Data(i,j) = Mean; // Background(i,j);
         io_write_ima_float("xx_interp.fits", Data, &Header);
       }
    }
    MB.decomposition(Data);
    MB.write_allima(Name_Imag_Out, Header);
    Data.init();
    MB.reconstruction(Data);
    if (MB.RemoveLastScale == True) Data += MB.LastScale;
    io_write_ima_float(Name_Imag_Out, Data, &Header);
    exit(0);
}
