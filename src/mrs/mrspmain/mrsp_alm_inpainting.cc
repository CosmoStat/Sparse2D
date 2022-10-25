/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  25/10/08
**    
**    File:  mrs_inpainting.cc 
**
*******************************************************************************
**
**    DESCRIPTION:  Inpainting + Iterative wiener filtering
**    -------------------
**
******************************************************************************/

#include"HealpixClass.h"
#include"MRS_Sparse.h"
#include"MRSP_Inp.h"
// #include"MRSP_Inp.cc"

char Name_Imag_In_T[512]; /* input file image */
char Name_Imag_In_Q[512]; /* input file image */
char Name_Imag_In_U[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Name_Mask_T[512]; 
char Name_Mask_Q[512]; 
char Name_Mask_U[512]; 
char Name_EB_Out[512];
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0;
#define  DEF_IWIENER_NBR_ITER  40
#define  DEF_ITER_GAUSSIANITY 10
 
int Max_Iter=DEF_IWIENER_NBR_ITER;
Bool UseRMSMap = False;   // Apply a denoising using an iterative Wiener filtering using a RMS map.
Bool UseMask = True;          // Apply an inpainting
Bool ThresholdResi = False;   // Threshold the Alm of the residual instead of the Alm of the solution
Bool WienerOK = False;        // Apply a denoising using an iterative Wiener filtering
Bool HardThreshold = False;  // Apply an Alm iterative hard threshold for the inpainting
Bool SoftThreshold = True;
Bool Use_ZeroMeanConstraint = True; // Apply a zero mean constraint on the solution
Bool Use_WP_Constraint = True;   // Apply a variance constraint in the wavelet packet of the solution
Bool EstimPowSpec = False;
Bool All_WP_Band = False;
Bool IsotropyCst = False;
Bool ProjectConstraint = False;
Bool EqualInMask = False;

float  ZeroPadding=0.0;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;
Bool UpdatePowSpec = True;

float Eps=1.;
Bool GaussianityConstraint = False;
float GaussCst_ProjectCst_Lambda = 5.;
Bool Acceleration=False;
Bool Pos=False;

Bool Analysis = False;
Bool OptMat = False;
Bool OldAnalysis = False;
Bool OptRevMat = False;
Bool EBWrite=False;

type_pol_inpainting InpaintMethod = DEF_POL_INP_METHOD;
Bool Denoising = False;
int AccelationLevel=1;
Bool UseBeam=False;
Bool EBMode=True;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map_T in_map_Q in_map_U in_mask_T in_mask_Q in_mask_U out_map [out_eb] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

//    fprintf(OUTMAN, "         [-m Inpainting method]\n");
//    for (int m=0; m < NBR_POLA_INP_METHOD; m ++)
//        fprintf(OUTMAN, "             %d: %s.\n", m+1, (char *) StringPolInpaintingMethod( (type_pol_inpainting) m)); 
//    fprintf(OUTMAN, "             default is %s\n", StringPolInpaintingMethod (InpaintMethod));

    
    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
        
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    
    fprintf(OUTMAN, "         [-A AccelerationLevel]\n");
    fprintf(OUTMAN, "             If set, the WP constraint is applied only once in N iterations, with N=AccelerationLevel. \n");
    fprintf(OUTMAN, "             When -A 0 is set, no WP constraint is used.  \n");   // Only valid with -m1 method. 

    fprintf(OUTMAN, "         [-Q]\n");
    fprintf(OUTMAN, "             Apply the inpainting on the Alm coef of Q and U maps instead of the Alm Coef of E and B modes. \n");
    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "             Force the output to be equal to the input in the mask. \n");
    fprintf(OUTMAN, "         [-W]\n");
    fprintf(OUTMAN, "             Don't appy the wavelet packet constraint. \n");

    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptI = False;
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "m:QECWA:l:pi:vzZ")) != -1) 
    {
	switch (c) 
        {   
	        case 'A':  if (sscanf(OptArg,"%d",&AccelationLevel) != 1) 
                       {
                          fprintf(OUTMAN, "Error: bad AccelationLevel: %s\n", OptArg);
                          exit(-1);
                       } 
                      Acceleration = True;  
                      break;
             case 'm':
				if (sscanf(OptArg,"%d",&c ) != 1) 
				{
					fprintf(OUTMAN, "Error: bad type of inpainting method: %s\n", OptArg);
					exit(-1);
 				}
                if ((c > 0) && (c <= NBR_POLA_INP_METHOD))  InpaintMethod = (type_pol_inpainting) (c-1);
                else  
                {
					fprintf(OUTMAN, "Error: bad type of inpainting method: %s\n", OptArg);
					exit(-1);
				}
                break;
            case 'E': EqualInMask = (EqualInMask == True) ? False: True;  break;
            case 'Q': EBMode = (EBMode == True) ? False: True;  break;
            case 'C': Use_ZeroMeanConstraint = (Use_ZeroMeanConstraint == True) ? False: True;  break;
         	case 'W': Use_WP_Constraint = (Use_WP_Constraint == True) ? False: True;  
        	           break;
            case 'a': All_WP_Band = (All_WP_Band == True) ? False: True;  
        	           break;
            case 'p': Pos=True; Use_ZeroMeanConstraint=False;break;
	        case 'l':
 	     	      if (sscanf(OptArg,"%d",&Lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
             case 'i':  
                   if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                   {
		            fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
		    		exit(-1);
		           } 
                   OptI = True;
                   if (Max_Iter <= 0)   Max_Iter =  DEF_IWIENER_NBR_ITER;
 	        	break;
    	    case 'v': Verbose = True; break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 

  
     if (OptI == False)  Max_Iter =  DEF_ITER_GAUSSIANITY;

      
        /* get optional input file names from trailing 
          parameters and open files */
    if (OptInd < argc) strcpy(Name_Imag_In_T, argv[OptInd++]);
    else usage(argv);          
    if (OptInd < argc) strcpy(Name_Imag_In_Q, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Imag_In_U, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Mask_T, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Mask_Q, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Mask_U, argv[OptInd++]);
    else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
    else usage(argv);

 
    if (OptInd < argc) 
	{
        strcpy(Name_EB_Out, argv[OptInd++]);
        EBWrite = True;
	}
    

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
  
/*********************************************************************/

class C_VectorField_Inpainting      
{
    CAlmR ALM_Rec, ALM_Rec1, ALM_Rec2;
    Hdmap Resi1, Resi2;
    fltarray WP_H, WP_G, WP_W, WP_WPFilter, TabPsi2;
    int NbrWP_Band;
    CAlmR  ALM_Band;
    Hdmap Band;
	arr<double> weight_T;

    double SigData1, SigData2;
    double MeanData1, MeanData2;
    double MinData1, MinData2;
    double MaxData1, MaxData2;
    int NbrTotMask1, NbrTotMask2; 
    double FSky1, FSky2;
    int Nside;
    int Lmax_Tot;
    double MinSigmaNoise;
    void init_inpainting();
    void exit_inpainting(CAlmR & ALM, Hdmap &Map, Hdmap &MaskMap, Hdmap &Result, double MinData, double MeanData);
    double get_residual();
    void alm_denoising(CAlmR & A, Bool Wiener=False);
    void init_wavelet_packet();
    void wp_band_constraint(Hdmap &Band, Hdmap & MaskMap,  int NbrTotMask, int b);
    void wp_constraint(CAlmR & ALM, Hdmap &Resi, Hdmap &Result, Hdmap & MaskMap, int NbrTotMask, int ALM_FIRST_L);

public:
     
    type_pol_inpainting InpaintMethod;
    
    Bool ThresholdResi;   // Threshold the Alm of the residual instead of the Alm of the solution
    Bool Verbose;         // Verbose mode
    int Max_Iter;         // Number of iteration in the inpainting algorithm
    Bool HardThreshold;   //  Apply an Alm iterative hard threshold for the inpainting
    Bool SoftThreshold;   // Apply an Alm iterative soft threshold for the inpainting
    float SigmaNoise;     // Noise standard deviation
     
    // Zero mean constrain on the solution 
    Bool Use_ZeroMeanConstraint; // Apply a zero mean constraint on the solution
    
     
    // Parameters for WP decomposition constraints
    Bool Use_WP_Constraint;      // Apply a variance constraint in the wavelet packet of the solution
    Bool All_WP_Band;
    float  ZeroPadding;   
    int Lmax;
    float Eps;
    
    bool Acceleration;
    int AccelationLevel;
    bool Pos;
    bool InputPowSpec; 
    Bool EqualInMask;
    
    Bool MCA;
    Hdmap *MCA_TabResult;  // MCA component maps
    CAlmR  MCA_ALM;  
    C_UWT  MCA_WT;  
    int MCA_WT_NbrScale;
    int MCA_Alm_Lmax;
    int NbrMCA;
    Bool MCA_InitCleanCMB;

    Hdmap MapT;     // Map to inpainted
    Hdmap Map1;     // Map to inpainted
    Hdmap Map2;     // Map to inpainted
    CAlmR  ALMT, ALM1, ALM2;  
    Hdmap MaskMap1; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    Hdmap MaskMap2; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    Hdmap MaskMapT;
    Hdmap ResultT; 
    Hdmap Result1;  // inpainted map
    Hdmap Result2;  // inpainted map
    Bool EBTrans;   // If true, the sparsity is applied in the EB domain.
    
    C_VectorField_Inpainting() {Max_Iter=DEF_IWIENER_NBR_ITER; ThresholdResi = False; HardThreshold = True; SoftThreshold = False; InpaintMethod = DEF_POL_INP_METHOD;
        Use_ZeroMeanConstraint = True; Use_WP_Constraint = False; All_WP_Band = False; MCA=False; NbrWP_Band=0;Nside=0;Lmax_Tot=0; NbrMCA=0;
        ZeroPadding=0.0; Lmax = 0.; Eps=1.;  EBTrans=True;
        Pos=False;Verbose=False;  SigmaNoise=0;MinSigmaNoise=0;EqualInMask=False;
        AccelationLevel=1; Acceleration=False; 
        MCA_InitCleanCMB=False;}
        
    void analyse_alm_inpainting();
    void forward_transform(Hdmap & D1, Hdmap & D2, CAlmR & A1, CAlmR & A2, Bool ToEB);
    void backward_transform(CAlmR & A1, CAlmR & A2, Hdmap & D1, Hdmap & D2, Bool FromEB);

    
    //   void analyse_alm_inpainting();
    //   void synthesis_alm_inpainting();
    //   void synthesis_alm_inpainting_with_master_decconv();
    
    ~C_VectorField_Inpainting(){};
};


/*********************************************************************/

/*********************************************************************/


void C_VectorField_Inpainting::init_inpainting()
{
    NbrTotMask1 = NbrTotMask2 = 0; 
    SigData1 = SigData2 = 0.;
    MeanData1 = MeanData2 = 0.;
    FSky1 = FSky2  = 1.;
    
    Map1.minmax(MinData1,MaxData1);
    Map1.minmax(MinData2,MaxData2);

    if (Verbose == True) Map1.info((char*) "Input MA1P");
    // Get the minimum, the stdev and the mean of the data, but in the mask
    for (int p=0; p<  Map1.Npix(); p++)  
    {
        if (MaskMap1[p] > 0)  
        {
            NbrTotMask1 ++;
            MeanData1 += Map1[p];
            if (MinData1 > Map1[p]) MinData1 = Map1[p];
         }
        if (MaskMap2[p] > 0)  
        {
            NbrTotMask2 ++;
            MeanData2 += Map2[p];
            if (MinData2 > Map2[p]) MinData2 = Map2[p];
        }
    }
    MeanData1 /= (double) NbrTotMask1;
    MeanData2 /= (double) NbrTotMask2;

    for (int p=0; p<  Map1.Npix(); p++) 
    { 
        if (MaskMap1[p] > 0)   SigData1 += (Map1[p]-MeanData1)*(Map1[p]-MeanData1);
        if (MaskMap2[p] > 0)   SigData2 += (Map2[p]-MeanData2)*(Map2[p]-MeanData2);
    }
    SigData1 = sqrt(SigData1 / (double) NbrTotMask1);
    SigData2 = sqrt(SigData2 / (double) NbrTotMask2);
    FSky1 = (double) Map1.Npix() / NbrTotMask1;
    FSky2 = (double) Map2.Npix() / NbrTotMask2;
    if (Verbose == True) cout << "Percentage of good pixels in Map 1: " << 1. / FSky1 * 100. <<  " MinData1 = " << MinData1 << endl;
    if (Verbose == True) cout << "Percentage of good pixels in Map 2: " << 1. / FSky2 * 100. <<  " MinData2 = " << MinData2 << endl;

    if (MCA == True)
    {
     	Use_ZeroMeanConstraint = False;
    }
    
    // We start with a zero mean image
    if (Use_ZeroMeanConstraint == True)
    {
        for (int p=0; p<  Map1.Npix(); p++) if (MaskMap1[p] > 0) Map1[p] -= MeanData1;
        for (int p=0; p<  Map2.Npix(); p++) if (MaskMap2[p] > 0) Map2[p] -= MeanData2;
    }  
    /*  if (Pos == True) 
     {
     for (int p=0; p<  Map.Npix(); p++) 
     if (MaskMap[p] > 0) Map[p] -= MinData;  
     }  */   
    Nside = Map1.Nside();
    Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT = " <<  Lmax_Tot << endl;
    
    ALM1.alloc(Nside, Lmax_Tot);
    // ALM.Norm = (IsotropyCst == True) ? True: False;
    ALM1.Norm = True;
    // ALM.Norm = False;
    ALM1.set_beam_eff(Lmax, Lmax_Tot);
    ALM1.UseBeamEff = False;
    ALM_Rec1.alloc(Nside, Lmax_Tot);
    ALM_Rec1.Norm = ALM1.Norm;
    ALM_Rec1.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Rec1.UseBeamEff = ALM1.UseBeamEff;   
    ALM2.alloc(Nside, Lmax_Tot);
    ALM2.Norm = ALM1.Norm;
    ALM2.set_beam_eff(Lmax, Lmax_Tot);
    ALM2.UseBeamEff = ALM1.UseBeamEff;
    ALM_Rec2.alloc(Nside, Lmax_Tot);
    ALM_Rec2.Norm = ALM1.Norm;
    ALM_Rec2.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Rec2.UseBeamEff = ALM1.UseBeamEff; 
    Result1.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    Result1.fill(0.);
    Resi1.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    Result2.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    Result2.fill(0.);
    Resi2.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    
/*    if (MCA == True)
    {
    	NbrMCA = 2;
    	MCA_ALM.alloc(Nside, Lmax_Tot);
    	MCA_ALM.Norm = ALM.Norm;
        if (MCA_Alm_Lmax ==0) MCA_Alm_Lmax = Lmax_Tot;
        MCA_TabResult = new Hdmap[2];
        for (int i=0; i < NbrMCA; i++)
        {
            MCA_TabResult[i].SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
            MCA_TabResult[i].fill(0.);
        }
        // WT init
        if ( MCA_WT_NbrScale <= 1)  MCA_WT_NbrScale=4;
        if (DEF_MRS_ORDERING != RING)  MCA_WT.wt_alloc(Nside, MCA_WT_NbrScale, Lmax_Tot, true);
        else  MCA_WT.wt_alloc(Nside, MCA_WT_NbrScale, Lmax_Tot);
        // wt_alloc(Nside, Lmax, DEF_MRS_ORDERING)
        if (MCA_WT_NbrScale <=0) MCA_WT_NbrScale = MCA_WT.nscale();
        if (Verbose ==True)  cout << "MCA WT NbrScale = " << MCA_WT.nscale() << endl;
        
        // Clean CMB first
        // Bool MCA_InitCleabCMB=False;
        if (MCA_InitCleabCMB == True)
        {
            C_UWT WT;
            
            int Ns = log((double) Nside) / log((double) 2.);
            MCA_TabResult[1] = Map;
            WT.wt_alloc(Nside, Ns, Lmax_Tot, DEF_MRS_ORDERING);
            float Sig=1.;
            WT.hard_thresholding(MCA_TabResult[1], (float) 5., Sig, true, true);
            (MCA_TabResult[1]).write("xx_wt_0.fits");
            for (int p=0; p < Map.Npix(); p++) 
            {
                if (ABS( (MCA_TabResult[1])[p] ) > 0) MaskMap[p]=0;
            }             
        }
    } */
    
    if (SigmaNoise == 0) SigmaNoise = 1.;
    MinSigmaNoise = SigmaNoise;
    // else if (ALM.Norm == false) SigmaNoise /= ALM.NormVal;
    
    // Find the minimum value in RMSMap in the valid part of the data
    
    
  // mrs_alloc_powspec(PowSpecSignal, Lmax_Tot);
        
    weight_T.alloc( 2*Nside );
    weight_T.fill(1);

#ifdef MRSDEBUG                      
    cout << "ALM.Norm = " << ALM.Norm << endl;
    cout << "ALM.NormVal = " << ALM.NormVal << endl;
    cout << "MinSigmaNoise = " << MinSigmaNoise << endl;
    cout << "Lmax_Tot  = " << Lmax_Tot << endl;
    cout << "MeanData  = " << MeanData1 << endl;
    cout << "MinData  = " << MinData1 << endl;
#endif
    
}


/*********************************************************************/
 
void C_VectorField_Inpainting::init_wavelet_packet()
{
    ALM_Band.alloc(Nside, Lmax_Tot);
    ALM_Band.Norm = ALM1.Norm;
    ALM_Band.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Band.UseBeamEff =  ALM1.UseBeamEff ;
    
    if (All_WP_Band == False) get_wp_meyer_filter(Nside, WP_H, WP_G, WP_W, WP_WPFilter, Lmax_Tot);
    else get_planck_wp_meyer_filter(WP_H, WP_G, WP_W, WP_WPFilter, Lmax, Lmax_Tot);
    NbrWP_Band = WP_H.ny();
    if (Verbose == True) cout << " NbrWP_Band = " << NbrWP_Band << endl;
    Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    
    if (WP_H.nx() != Lmax_Tot+1)
    {
        cout << "Error: bad dimension of the H filter " << endl;
        exit(-1);
    }
}

/*********************************************************************/

void C_VectorField_Inpainting::forward_transform(Hdmap & D1, Hdmap & D2, CAlmR & A1, CAlmR & A2, Bool ToEB)
{
   if (ToEB == False)
   {
     A1.alm_trans(D1);
     A2.alm_trans(D2);
   }
   else
   {
       // map2alm_pol_iter_QU( D1, D2, A1, A2, 0, weight_T);
       map2alm_spin(D1,D2,A1,A2,2,weight_T,false);
   }
}
/*********************************************************************/

void C_VectorField_Inpainting::backward_transform(CAlmR & A1, CAlmR & A2, Hdmap & D1, Hdmap & D2, Bool FromEB)
{
    if (FromEB == False)
    {
        A1.alm_rec(D1, True);
        A2.alm_rec(D2, True);
    }
    else
    {
        // alm2map_pol_QU( A1, A2, D1, D2 );
        alm2map_spin(A1,A2,D1,D2,2);
    }
}
  
/*********************************************************************/

double C_VectorField_Inpainting::get_residual()
{
    double Err = 0.;
    
    for (int p=0; p < Result1.Npix(); p++) 
    {
      /*  if (MCA == True)
        {
            Result[p] = 0;
            for (int i=0; i < NbrMCA; i++) Result[p] += (MCA_TabResult[i])[p];
        }*/
        Resi1[p] =  Map1[p] - Result1[p];
        if (MCA != True) Resi1[p] *= MaskMap1[p];
        
        Err += Resi1[p]*Resi1[p];
    }
    for (int p=0; p < Result2.Npix(); p++) 
    {
        /*  if (MCA == True)
         {
         Result[p] = 0;
         for (int i=0; i < NbrMCA; i++) Result[p] += (MCA_TabResult[i])[p];
         }*/
        Resi2[p] =  Map2[p] - Result2[p];
        if (MCA != True) Resi2[p] *= MaskMap2[p];
        
        Err += Resi2[p]*Resi2[p];
    }
    Err = sqrt(Err / (double) (NbrTotMask1 + NbrTotMask2)) / MinSigmaNoise;
    return Err;
}

/*********************************************************************/

void C_VectorField_Inpainting::wp_band_constraint(Hdmap &Band, Hdmap & MaskMap, int NbrTotMask, int b)
{
    double EnergyInMap = 0.;
    for (int p=0; p < MaskMap.Npix(); p++)  if (MaskMap[p] != 0) EnergyInMap += Band[p]*Band[p];
    double Sig_InMap  = sqrt(EnergyInMap  /  (double) (NbrTotMask));
    
    
    // Compute the expected power in the band
    double EnergyBand=0.;
    double Np = (double) MaskMap.Npix();
    // if (InputPowSpec == True)  EnergyBand = Np * TabPsi2(b);
    // else 
    EnergyBand = EnergyInMap * Np / (double) NbrTotMask;
    double SigBand = (EnergyBand > 0) ? sqrt(EnergyBand / Np): 0;
    
    // double EnergyBand = Map.Npix() * SigBand*SigBand;
    
    double EnergyInMask = 0.;
    for (int p=0; p < MaskMap.Npix(); p++) 
    {
        // float Level = GaussCst_ProjectCst_Lambda * Sig_InMap;
        if (MaskMap[p] == 0) 
        {
           EnergyInMask += Band[p]*Band[p];
        } 
    }
    double Sig_Mask = sqrt( EnergyInMask /  (double) (Np - NbrTotMask));
    
    // Band renormalization
    float Coeff =  1.;
    double ExpectedEnergyInMask = EnergyBand - EnergyInMap;
    if (ExpectedEnergyInMask <= 0) ExpectedEnergyInMask = EnergyBand * (Np -  NbrTotMask) / Np;
    double E = EnergyInMap  * (Np -  NbrTotMask) / (double) NbrTotMask;
    // cout << "      TabPsi2 = " << sqrt(TabPsi2(b)) << ", ExpectedSigmaInMask = " << sqrt(ExpectedEnergyInMask / (Np -  NbrTotMask)) << " " << sqrt(E /  (Np -  NbrTotMask)) <<  " " << Sig_Mask << endl;
    if (EnergyInMask > 0) 
        Coeff =  sqrt(ExpectedEnergyInMask  / EnergyInMask);
    if (Coeff > 5) Coeff = 5.;
    
    // if ((Sig_Mask > 0) && (Sig_InMap > 0)) Coeff = Sig_InMap / Sig_Mask;
    //  if (Verbose == True) cout << "        WP:Band " << b+1 << ",   Coeff " <<  Coeff << " Sig_Mask = " << Sig_Mask << ",  Sig_InMap = " << Sig_InMap <<  "   Npix =  " <<  Map.Npix() << endl;
    
    for (int p=0; p < MaskMap.Npix(); p++)  if (MaskMap[p] == 0)  Band[p] *= Coeff;
}

/*********************************************************************/


void C_VectorField_Inpainting::wp_constraint(CAlmR & ALM, Hdmap &Resi, Hdmap &Result, Hdmap & MaskMap, int NbrTotMask, int ALM_FIRST_L)
{
   Resi.init();
   for (int b=0; b <= NbrWP_Band; b++)
   {  
       int LMin = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,0): WP_W(0,0);
       int LMax = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,1): WP_W(0,1);
    
       //  if (Verbose == True) cout << "        WP:Band " << b+1 << ", Lmin = " << LMin << ", Lmax = " << LMax << endl; 
           // If the ALM_FIRST_L < to the first non zero coef of the filter then band is zero.
       if  (LMin > ALM_FIRST_L) Band.init();
       else 
       {
            if (b == NbrWP_Band) Band = Result;
            else   
            {
               ALM_Band.Set(LMax,  LMax);
               ALM_Band.SetToZero();

               // Compute the solution at a given resolution (convolution with H)
               int FL = MIN(ALM_Band.Lmax(),WP_H.nx()-1);
               if (FL > ALM.Lmax()) FL = ALM.Lmax();
               for (int l=0; l <= FL; l++)
               for (int m=0; m <= l; m++)
               {
                    //  ALM_Band(l,m).re = ALM(l,m).re * WP_H(l, NbrWP_Band-1-b);
                    // ALM_Band(l,m).im = ALM(l,m).im * WP_H(l, NbrWP_Band-1-b);
                    ALM_Band(l,m) = ALM(l,m) * xcomplex<REAL>(WP_H(l, NbrWP_Band-1-b), 0.);
                }
                ALM_Band.alm_rec(Band);
                    
                    /* {
                     char FN[256];
                     sprintf(FN, "xx_r%d.fits", b+1);
                     Result.write(FN);
                     sprintf(FN, "xx_l%d.fits", b+1);
                     Band.write(FN);
                     } */
                    
                    
                // Compute the coefficients and standard deviations in and out the mask
                for (int p=0; p < Result.Npix(); p++) 
                {
                    double Val = Band[p];
                    Band[p] = Result[p] - Val;
                    Result[p] = Val;
                }
            } // endelse if (b == NbrWP_Band)
            
            /*{
             char FN[256];
             sprintf(FN, "xx_b%d.fits", b+1);
             Band.write(FN);
             }*/
            
            // Wavelet packet constraint on the band b
           wp_band_constraint(Band, MaskMap, NbrTotMask, b);

            // coadd the band to new new solution
            Resi += Band;
       } // else if  (LMin > ALM_FIRST_L) 
   } // endfor b
   Result = Resi;
}


/*********************************************************************/

void C_VectorField_Inpainting::exit_inpainting(CAlmR & ALM, Hdmap &Map, Hdmap &MaskMap, Hdmap &Result, double MinData, double MeanData)
{
    double ValAdd = 0.;
    if (Pos == True) ValAdd += MinData1;
    if (Use_ZeroMeanConstraint == True) ValAdd += MeanData;
    
    if (ValAdd != 0) 
    {
        for (int p=0; p<  Map.Npix(); p++)
        {
            if (MaskMap[p] > 0) Map[p] += ValAdd;
            Result[p] += ValAdd;
        }
    }
    
    // We reinsert the correct pixel values 
    if (EqualInMask == True)
       for (int p=0; p < Map.Npix(); p++) if (MaskMap[p] > 0) Result[p] = Map[p];
    
    if (ALM.Norm == True)
    {
        double CoefN = ALM.NormVal;
        for (int l=0; l <= ALM.Lmax(); l++)
            for (int m=0; m <= l; m++)
            {
                ALM(l,m) /= CoefN;
            }
    }
}


/*********************************************************************/


void C_VectorField_Inpainting::analyse_alm_inpainting()
{
    double Sig_InMap, Sig_Mask;
    Bool Debug = False;
#ifdef MRSDEBUG                      
    Debug = True;
#endif
    
    init_inpainting();
    if (Use_WP_Constraint == True) init_wavelet_packet();

  
    // if (Verbose == True)  Map.info((char *) "INIT");
    // Map.write("xx_tt.fits");
    
    if (Debug == True) Map1.info("Map1 ");
    if (Debug == True) Map2.info("Map2 ");
    
    forward_transform(Map1, Map2, ALM1, ALM2, EBTrans);
    // ALM.write((char *) "xx_init_alm.fits", True);
    
    double MaxAlm1 =ALM1.max_absalm();
    double MaxAlm2 =ALM2.max_absalm();
    double MaxAlm = MAX(MaxAlm1,MaxAlm2); 
    // IDL: NormVal =  sqrt(nel / (4.*!DPI))
    MaxAlm *= 0.9999999;
    if (ALM1.Norm == False) SigmaNoise /= ALM1.NormVal;  // Normally ALM.Norm should be set to True,
    
    //  if (Debug == True)
    {
        if (ALM1.Norm == True) cout << " Alm norm: max =  " << MaxAlm <<  endl;
        else cout << " Alm nonorm: max =  " << MaxAlm << endl;
        printf(" MaxAlm = %5.4f, SigNoise = %f\n", (float) MaxAlm / SigmaNoise, SigmaNoise);
    }

      
    double FirstSoftThreshold = MaxAlm / SigmaNoise * sqrt(2.);  // sqrt(2.) is here, to be equivalent to the IDL code.
    // this has no effect since we redivide latter by the same value.
    double LastSoftThreshold = 0.;
    double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    // double StepL = DeltaThreshold / (float) (Max_Iter -1);
    double Lambda = FirstSoftThreshold;
    if (Debug == True)  cout << "==> Max Alm = " << MaxAlm << " " << (MaxAlm * SigmaNoise)  << ", N = " << ALM1.NormVal <<  ", First Threshold = " << Lambda << endl;
    
    double Err, Err1=Map1.rms();
    double Err2=Map2.rms();
    int MaxNonZeroL1=0; 
    int MaxNonZeroL2=0; 
    int ALM_FIRST_L = (Acceleration == True) ? 100: Lmax_Tot;
    if (IsotropyCst == True) ALM_FIRST_L = Lmax_Tot;
    // int Ni = Max_Iter;
    // int Fi = 0;
    
    // Use for isotropic constrain only
    float GaussCst_LambdaFirst = 0.5;
    float GaussCst_StepLambda = (GaussCst_ProjectCst_Lambda - GaussCst_LambdaFirst)  / float(Max_Iter);
    float GaussCst_Lambda = GaussCst_LambdaFirst;
    
    // fits_write_fltarr((char *)"xxpsi.fits", TabPsi2);
    // wp_trans(Map,  Band,  ALM,  ALM_Band,  WP_W,  WP_H,  NbrWP_Band,  ALM_FIRST_L, TabPsi2);
    // exit(0);
    
    // MAIN ITER 
    if (Debug == True)
    {
        char NN[256];
        sprintf(NN, "xx_input1.fits");
        Map1.write(NN);
    }
    Lambda = FirstSoftThreshold;
    
   
    for (int i=0; i < Max_Iter; i++)
    { 
        Err = get_residual();
        if (ThresholdResi == False) 
        {
            Resi1 += Result1;
            Resi2 += Result2;
        }
        
        if (Debug == True)
        {
            char NN[256];
            sprintf(NN, "xx_resi1_it_%d.fits", i+1);
            if (i >= 0) Resi1.write(NN);
        }
        int CptAlm=0;
        
        if (Acceleration == True) 
        {
           ALM1.Set(ALM_FIRST_L,ALM_FIRST_L);
           ALM2.Set(ALM_FIRST_L,ALM_FIRST_L);
        }
        forward_transform(Resi1, Resi2, ALM1, ALM2, EBTrans);

        // MaxAlm =ALM.max_absalm() / SigmaNoise;
        // cout << " Max ABS(ALM) Resi  " << MaxAlm  << endl;
        if (i > 0) Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * ((float) i / (float) Max_Iter)));  // 4 versus 2.8  in IDL code
        if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;	      
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        if (Debug == True) printf(" Lambda = %5.4f\n", (float) Lambda);
        if (Verbose == True) cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
        if (Debug == True) Resi1.info("RESI1 ");
        if (Debug == True) Resi1.info("RESI2 ");

        // if (Denoising == True) alm_denoising(ALM, WienerFilter);
        
        if ((HardThreshold == True)  || (SoftThreshold == True))  // Alm L1 minimization
        { 
            float ThresholdLevel = SigmaNoise*Lambda / sqrt(2.);
            if (Lambda > 0)
            {
                if (HardThreshold == True) 
                {
                    // cout << "Max AlmResi = " << ALM.max_absalm() << endl;
                    // ALM.write("xxalm.fits", true);
                    CptAlm = ALM1.hard_threshold(ThresholdLevel, MaxNonZeroL1);
                    CptAlm += ALM2.hard_threshold(ThresholdLevel, MaxNonZeroL2);
                }
                else 
                {
                   if (SoftThreshold == True)   CptAlm = ALM1.soft_threshold(ThresholdLevel/sqrt(2.), MaxNonZeroL2);
                   if (SoftThreshold == True)   CptAlm += ALM2.soft_threshold(ThresholdLevel/sqrt(2.), MaxNonZeroL2);
                }
                if (Debug == True) cout << "  Number of significant Alm : = " << CptAlm <<  " MaxL = " << MaxNonZeroL1 << " " << MaxNonZeroL2 << endl; 
            }
            else  MaxNonZeroL1 = MaxNonZeroL2 = Lmax_Tot;
            // ALM.info((char *) "ALM INFO  1");
            
            
            // Update the solution
            if (ThresholdResi == True) 
            {
                backward_transform(ALM1, ALM2, Resi1, Resi2, EBTrans);
                for (int p=0; p<  Result1.Npix(); p++)   Result1[p]  +=  Eps*Resi1[p];
                for (int p=0; p<  Result2.Npix(); p++)   Result2[p]  +=  Eps*Resi2[p];
                forward_transform(Result1, Result2, ALM1, ALM2, EBTrans);

            }   
            else backward_transform(ALM1, ALM2, Result1, Result2, EBTrans);
        
            // To accelerate, we compute the Alm only to MaxNonZeroL*2
            // if (Acceleration == True)  
            {
                int MaxNonZeroL = MAX(MaxNonZeroL1, MaxNonZeroL2);
                if (2*MaxNonZeroL >  ALM_FIRST_L)  ALM_FIRST_L = 2*MaxNonZeroL;
                ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
            }
        } // end l1 minimization i.e.  ((HardThreshold == True)  || (SoftThreshold == True))
	    
                
        //             ALM1.info((char *) "ALM INFO  2");
        //             ALM2.info((char *) "ALM INFO 3");
        
        // Apply constraint per wavelet packet band. We start these constraint only at half the total number of iterations.
        // Use_WP_Constraint = False;
        
        
        if ((Use_WP_Constraint == True) && (i > 0) && (i % AccelationLevel == 0))
        {
            if (EBTrans == True) forward_transform(Result1, Result2, ALM1, ALM2, False);
            if (Debug == True) cout << "        Wavelet Packets constraint ... " << endl;
            // ALM.info((char *) "ALM INFO 4");
            
            // We compute the WP transform:
            
            //   ITER: Result contains the high resolution   (c_j)
            //   We compute Band, which contains the low resolution (c_{j+1}
            //   Coef =  Result - Band
            //   Result = Band
            wp_constraint(ALM1, Resi1, Result1, MaskMap1, NbrTotMask1, ALM_FIRST_L);
            wp_constraint(ALM2, Resi2, Result2, MaskMap2, NbrTotMask2, ALM_FIRST_L);
            
            // Here Result1 and 2 are the E and B maps and we want to reconstruct the Q and U map
            // if (EBTrans == True)  forward_transform(Result1, Result2, ALM1, ALM2, EBTrans);
        }   // endif Use_WP_Constraint
         
        
        // We want a zero mean solution
        if (Use_ZeroMeanConstraint == True)
        {
            if (Debug == True)  cout << "ZERO" << endl;
            double MeanRes = Result1.average();
            for (int p=0; p < Map1.Npix(); p++)  Result1[p] -= MeanRes;
            MeanRes = Result2.average();
            for (int p=0; p < Map2.Npix(); p++)  Result2[p] -= MeanRes;
        }
        
        if (Pos == True)
        {
            if (Debug == True)  cout << "POS" << endl;
            for (int p=0; p < Map1.Npix(); p++) if (Result1[p] < 0) Result1[p] = 0.;
            for (int p=0; p < Map2.Npix(); p++) if (Result2[p] < 0) Result2[p] = 0.;
        }
                
        if (Debug == True)
        {
            char NN[256];
            sprintf(NN, "xx_res1_it_%d.fits", i+1);
            if (i >= 0) Result1.write(NN);
            sprintf(NN, "xx_res2_it_%d.fits", i+1);
            if (i >= 0) Result2.write(NN);
        }
        if (Verbose == True) Result1.info((char *) "  REC");
        if (Verbose == True) Result2.info((char *) "  REC");
    }
 
    exit_inpainting(ALM1, Map1, MaskMap1, Result1, MinData1, MeanData1);
    exit_inpainting(ALM2, Map2, MaskMap2, Result2, MinData2, MeanData2);

}

/*********************************************************************/


int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);

    
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in Q = " << Name_Imag_In_Q << endl;
        cout << "# File Name in U = " << Name_Imag_In_U << endl;
        cout << "# Prefixe File Name Out = " << Name_Imag_Out << endl;   
        cout << "# Number of iterations = " <<  Max_Iter << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;    
        if (EBMode== False)  cout << "# EBmode is off " << endl;
        if (Use_WP_Constraint == True) cout << "# Use WP constraint " << endl;
    }
           
#ifdef MRSDEBUG                      
        cout << "DEBUG MODE " << endl;
#endif
     
    C_VectorField_Inpainting INP;
    
    INP.MapT.read(Name_Imag_In_T);
    INP.Map1.read(Name_Imag_In_Q);
    INP.Map2.read(Name_Imag_In_U);
    INP.MaskMapT.read(Name_Mask_T);
    INP.MaskMap1.read(Name_Mask_Q);
    INP.MaskMap2.read(Name_Mask_U);
    
	if( INP.Map2.nside() != INP.Map1.nside() )
	{
		fprintf( OUTMAN, "Error: Maps have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
    if( INP.MaskMap1.nside() != INP.Map1.nside() )
	{
		fprintf( OUTMAN, "Error: Q Map and Mask have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
    if( INP.MaskMap2.nside() != INP.Map2.nside() )
	{
		fprintf( OUTMAN, "Error: U Map and Mask have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
	if( INP.MaskMapT.nside() != INP.MapT.nside() )
	{
		fprintf( OUTMAN, "Error: T Map and Mask have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
    for (int p=0; p < INP.MaskMapT.Npix(); p++)  
   {
      if (INP.MaskMapT[p] != 1) 
      {
         INP.MaskMapT[p] = 0;
         INP.MapT[p] = 0;
      }
   }

   for (int p=0; p < INP.MaskMap1.Npix(); p++)  
   {
      if (INP.MaskMap1[p] != 1) 
      {
         INP.MaskMap1[p] = 0;
         INP.Map1[p] = 0;
      }
   }
   for (int p=0; p < INP.MaskMap2.Npix(); p++)  
   {
      if (INP.MaskMap2[p] != 1) 
      {
        INP.MaskMap2[p] = 0;
        INP.Map2[p] = 0;
      }
   }
  
   int Nside = INP.Map1.Nside();
   int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
   if (Verbose == True) cout << "# Used Lmax = " <<  Lmax_Tot << endl;
   // mrsp_alloc_powspec(INP.PolPowSpecData, Lmax_Tot);
      
   INP.Verbose = Verbose;
   INP.HardThreshold = HardThreshold;
   INP.SoftThreshold = SoftThreshold;
   INP.ThresholdResi = ThresholdResi;
   INP.SigmaNoise = SigmaNoise;
   INP.Use_ZeroMeanConstraint = Use_ZeroMeanConstraint;
   INP.Use_WP_Constraint = Use_WP_Constraint;
   INP.All_WP_Band = All_WP_Band;
   INP.ZeroPadding = ZeroPadding;
   INP.Lmax = Lmax;
   INP.Eps=Eps;
   INP.Acceleration=Acceleration;
   INP.AccelationLevel = AccelationLevel;
   INP.Pos=Pos;
   INP.Max_Iter = Max_Iter;
   INP.EBTrans  = EBMode;
   INP.EqualInMask = EqualInMask;
   
   switch (InpaintMethod)
   {
       case POL_INP_L1_ALM_ANALYSIS_CST_VAR:
            Analysis = True;
            if (AccelationLevel == 0) INP.Use_WP_Constraint = False;
            INP.analyse_alm_inpainting();
            break;
        default: cout << "Error: method not implemented ... " << endl;
            exit(-1);
            break;
    }

    char FN[512];
    
    sprintf(FN, "%s_mapQ.fits", Name_Imag_Out);
    INP.Result1.write(FN);
    sprintf(FN, "%s_mapU.fits", Name_Imag_Out);
    INP.Result2.write(FN);

    if (EBWrite == True)  
    {
        // cout << "EstimPowSpec " << endl;  
        INP.forward_transform(INP.Result1, INP.Result2, INP.ALM1, INP.ALM2, True);
        PowSpec SigPowSpec1(1, Lmax_Tot);
        INP.ALM1.alm2powspec(SigPowSpec1);
        sprintf(FN, "%s_psE.fits", Name_EB_Out);
        mrs_write_powspec( FN,  SigPowSpec1);    

        PowSpec SigPowSpec2(1, Lmax_Tot);
        INP.ALM2.alm2powspec(SigPowSpec2);
        sprintf(FN, "%s_psB.fits", Name_EB_Out);
        mrs_write_powspec( FN,  SigPowSpec2);     
       
        sprintf(FN, "%s_almE.fits", Name_EB_Out);
        INP.ALM1.write(FN);
        
        sprintf(FN, "%s_almB.fits", Name_EB_Out);
        INP.ALM2.write(FN);
       
        INP.ALM1.alm_rec(INP.Result1, True);
        sprintf(FN, "%s_mapE.fits", Name_EB_Out);
        INP.Result1.write(FN);
        
        INP.ALM2.alm_rec(INP.Result2, True);
        sprintf(FN, "%s_mapB.fits", Name_EB_Out);
        INP.Result2.write(FN);
   }
   
   exit(0);
}

