/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck with addition by Florent Sureau
**
**    Date:  25/10/08 // 05/2017
**    
**    File:  mrsp_alm_inpainting_spectra.cc 
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
#include "PolaHealpixClass.h"
// #include"MRSP_Inp.cc"
#include "mrsp_lib_matmask.h"

char Name_Imag_In_TQU[1024]; /* input file image */
char Name_Imag_Out[1024]; /* output file name */
char Name_Mask_TQU[1024]; /* Name of TQU mask */
char Name_TEB_Out[1024]; /* Name for optional TEB map out */
char Name_InSpectra[1024];/* Name for optional input spectra*/
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0;
#define DEF_IWIENER_NBR_ITER  40
#define DEF_ITER_GAUSSIANITY 10
//#define MRS_SPEC_DBG
 
int Max_Iter=DEF_IWIENER_NBR_ITER;
// Apply a denoising using an iterative Wiener filtering using a RMS map.
Bool UseRMSMap = False;
// Apply an inpainting
Bool UseMask = True;
// Threshold the Alm of the residual instead of the Alm of the solution
Bool ThresholdResi = False;
// Apply a denoising using an iterative Wiener filtering
Bool WienerOK = False;
// Apply an Alm iterative hard threshold for the inpainting
Bool HardThreshold = False;
// Apply an Alm iterative soft threshold for the inpainting
Bool SoftThreshold = True;
// Apply a zero mean constraint on the solution
Bool Use_ZeroMeanConstraint = False;
Bool ForceMonopDip = False;
// Apply a variance constraint in the wavelet packet of the solution
Bool Use_WP_Constraint = True;
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
Bool TEBWrite=False;

type_pol_inpainting InpaintMethod = DEF_POL_INP_METHOD;
Bool Denoising = False;
int AccelationLevel=1;
Bool UseBeam=False;
Bool EBMode=True;
Bool AlmFast=False;
Bool SpectraSVD=False;
Bool RestrictSpectraSVD=False;
Bool ThresholdPerChannel=False;
Bool FixedSpectra=False;

/******************************************************************************/
static void usage(char *argv[])
{
    fprintf(OUTMAN, 
      "Usage: %s options in_map_TQU in_mask_TQU out_map [out_eb]\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN, "         [-a Use all wavelet packet bands]\n");
    fprintf(OUTMAN, "             Default is reduced number of bands.\n");
    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    fprintf(OUTMAN, "         [-F fast]\n");
    fprintf(OUTMAN, "             Fast Alm Transform \n");
    fprintf(OUTMAN, "         [-A AccelerationLevel]\n");// Only with -m1.
    fprintf(OUTMAN, "             If set, the WP constraint is applied only\n");
    fprintf(OUTMAN, "             every N iter, with N=AccelerationLevel.");
    fprintf(OUTMAN, "             When -A 0 set, no WP constraint is used.\n");
    fprintf(OUTMAN, "         [-Q]\n");
    fprintf(OUTMAN, "             Apply the inpainting on Alm coefs Q/U.\n");
    fprintf(OUTMAN, "             instead of the Alm Coef of E/B. \n");
    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "             Force output to input in the mask. \n");
    fprintf(OUTMAN, "         [-W]\n");
    fprintf(OUTMAN, "             No wavelet packet constraint. \n");
    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "             Apply spectral constraints on SVD space. \n");
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "             Set TB,EB spectra to 0 \n");
    fprintf(OUTMAN, "         [-T]\n");
    fprintf(OUTMAN, "             Threshold Pola and Temperature differently\n");
    fprintf(OUTMAN, "         [-s FSpectra]\n");
    fprintf(OUTMAN, "             Used fix set of spectra given in FSpectra\n");

    exit(-1);
}
 
/******************************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptI = False;
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "am:QEFCDWSRTA:l:pi:s:vzZ")) != -1) {
        switch (c) {
            case 'A':  
                if (sscanf(OptArg,"%d",&AccelationLevel) != 1) {
                    fprintf(OUTMAN, "Error: bad AccelationLevel: %s\n", OptArg);
                    exit(-1);
                }
                Acceleration = True;  
                break;
            case 'm':
                if (sscanf(OptArg,"%d",&c ) != 1) {
                    fprintf(OUTMAN,"Error: bad type of inpainting method: %s\n",
                                                                        OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_POLA_INP_METHOD))  
                                    InpaintMethod = (type_pol_inpainting) (c-1);
                else {
                    fprintf(OUTMAN,"Error: bad type of inpainting method: %s\n",
                                                                        OptArg);
                    exit(-1);
                }
                break;
            case 'E': EqualInMask = (EqualInMask == True) ? False: True;  break;
            case 'F': AlmFast = True;break;
            case 'R': RestrictSpectraSVD = True;break;
            case 'S': SpectraSVD = True;break;
            case 'Q': EBMode = False; break;
            case 'C': Use_ZeroMeanConstraint = True;
                break;
            case 'D': ForceMonopDip = True;
                break;
            case 'W': Use_WP_Constraint = False; break;
            
            case 'T': ThresholdPerChannel=True;break;
            case 'a': All_WP_Band = True; break;
            case 'p': Pos=True; Use_ZeroMeanConstraint=False;break;
            case 'l':
                if (sscanf(OptArg,"%d",&Lmax) != 1) {
                    fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
                    exit(-1);
                }
                break;
            case 's':
                if (sscanf(OptArg,"%s",Name_InSpectra) != 1) {
                    fprintf(OUTMAN, "Error: bad s value: %s\n", OptArg);
                    exit(-1);
                }
                FixedSpectra=True;
                break;
            case 'i':  
                if (sscanf(OptArg,"%d",&Max_Iter) != 1) {
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
    /* get optional input file names from trailing parameters and open files */
    if (OptInd < argc) strcpy(Name_Imag_In_TQU, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Mask_TQU, argv[OptInd++]);
    else usage(argv);

    if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
    else usage(argv);
    
    if (OptInd < argc) {
        strcpy(Name_TEB_Out, argv[OptInd++]);
        TEBWrite = True;
    }
    /* make sure there are not too many parameters */
    if (OptInd < argc){
        fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
        exit(-1);
    }
}

/******************************************************************************/
class Pola_VectorField_Inpainting {
    //CAlmR ALM_Rec, ALM_Rec1, ALM_Rec2;
    PolaAlmR ALM_Rec;
    PolaHdmap PolaResi;
    //Hdmap Resi1, Resi2;
    fltarray WP_H, WP_G, WP_W, WP_WPFilter, TabPsi2;
    int NbrWP_Band;
    CAlmR ALM_Band;
    Hdmap Band;
    arr<double> weight_T;
    dblarray* ProjAlm;
    //double SigData1, SigData2;
    double SigDataPola[3];
    //double MeanData1, MeanData2;
    double MeanDataPola[3];
    //double MinData1, MinData2;
    double MinDataPola[3];
    //double MaxData1, MaxData2;
    double MaxDataPola[3];
    //int NbrTotMask1, NbrTotMask2; 
    unsigned long NbrTotMask[3];
    unsigned long NbrOverallMask;
    //double FSky1, FSky2;
    double FSky[3];
    int Nside;
    int Lmax_Tot;
    double MinSigmaNoise;
    void init_inpainting();
    void exit_inpainting(CAlmR &ALM,Hdmap &Map,Hdmap &MaskMap,Hdmap &Result,
                                                double MinData,double MeanData);
    double get_residual();
    void alm_denoising(CAlmR & A, Bool Wiener=False);
    void init_wavelet_packet();
    void wp_band_constraint(Hdmap &Band,Hdmap & MaskMap,int NbrTotMask,int b);
    void wp_constraint(CAlmR &ALM,Hdmap &Resi,Hdmap &Result,Hdmap &MaskMap,
                                                int NbrTotMask,int ALM_FIRST_L);
    void getSpectraEigenVectors(PolaAlmR &A1);
    void projectToEigenVectors(PolaAlmR &A1,bool back);

public:
     
    type_pol_inpainting InpaintMethod;
    
    Bool ThresholdResi;   // Threshold Alm residual instead of solution
    Bool Verbose;         // Verbose mode
    int Max_Iter;         // Number of iteration in the inpainting algorithm
    Bool HardThreshold;   // Alm iterative hard threshold for the inpainting
    Bool SoftThreshold;   // Alm iterative soft threshold for the inpainting
    //float SigmaNoise;   // Noise standard deviation
    float SigmaNoise[3];  // Noise standard deviation
     
    Bool Use_ZeroMeanConstraint; // Apply a zero mean constraint on the solution
    
     
    // Parameters for WP decomposition constraints
    Bool Use_WP_Constraint; //Variance constraint via wavelet packet on solution
    Bool All_WP_Band;
    float ZeroPadding;
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
   
    //Hdmap Map1;     // Map to inpainted
    //Hdmap Map2;     // Map to inpainted
    PolaHdmap Map_TQU;
    PolaAlmR  ALMPola;  
    Healpix_PowSpec FixedSpectraArray;
    //CAlmR  ALM1, ALM2;  
    PolaHdmap Mask_TQU;
    //Hdmap MaskMap1; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    //Hdmap MaskMap2; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    //Hdmap Result1;  // inpainted map
    //Hdmap Result2;  // inpainted map
    PolaHdmap Result_TQU;
    Bool EBTrans;   // If true, the sparsity is applied in the EB domain.
    
    Pola_VectorField_Inpainting() {
        Max_Iter=DEF_IWIENER_NBR_ITER; ThresholdResi = False;
        HardThreshold = True; SoftThreshold = False; 
        InpaintMethod = DEF_POL_INP_METHOD;
        Use_ZeroMeanConstraint = True; Use_WP_Constraint = False; 
        All_WP_Band = False; MCA=False; NbrWP_Band=0;Nside=0;Lmax_Tot=0; 
        NbrMCA=0;ZeroPadding=0.0; Lmax = 0.; Eps=1.;  EBTrans=True;
        Pos=False;Verbose=False;EqualInMask=False;
        SigmaNoise[0]=SigmaNoise[1]=SigmaNoise[2]=0;MinSigmaNoise=0;
        AccelationLevel=1; Acceleration=False; ProjAlm=NULL;
        MCA_InitCleanCMB=False;}
    void analyse_alm_inpainting();
    void forward_transform(PolaHdmap &D1, PolaAlmR &A1);
    void backward_transform(PolaAlmR &A1,PolaHdmap &D1);
    void read_spectra(char *Name, int lLmax);

    
    //   void analyse_alm_inpainting();
    //   void synthesis_alm_inpainting();
    //   void synthesis_alm_inpainting_with_master_decconv();
    
    ~ Pola_VectorField_Inpainting(){
        if(ProjAlm!=NULL) delete ProjAlm;
    };
};


/*********************************************************************/
void Pola_VectorField_Inpainting::init_inpainting()
{
    if(Mask_TQU.Npix()!=Map_TQU.Npix()){
        printf("Inconsistency: Npix Map=%ld, while Npix Mask=%ld\n",
                                               Map_TQU.Npix(), Mask_TQU.Npix());
        exit(EXIT_FAILURE);
    }
    for(int kmap=0;kmap<3;++kmap){
        NbrTotMask[kmap]=0;
        SigDataPola[kmap]= MeanDataPola[kmap]=0.0;
        FSky[kmap]=1.0;
        Map_TQU.get_map(kmap)->minmax(MaxDataPola[kmap],MinDataPola[kmap]);
    }
    if (Verbose == True) Map_TQU.get_map(0)->info((char*) "Input MAP 1");
    // Get the minimum, the stdev and the mean of the data, but in the mask
    NbrOverallMask=0;
    for (int kmap=0; kmap<3;++kmap) {
        Hdmap* CurMask=Mask_TQU.get_map(kmap);
        Hdmap* CurMap=Map_TQU.get_map(kmap);
        for (int p=0; p<Map_TQU.Npix();p++) {
            if ((*CurMask)[p] > 0){
                NbrTotMask[kmap]++;
                MeanDataPola[kmap] += (*CurMap)[p];
                if (MinDataPola[kmap] > (*CurMap)[p]) 
                                                MinDataPola[kmap]=(*CurMap)[p];
                if (MaxDataPola[kmap] < (*CurMap)[p])
                                                MaxDataPola[kmap]=(*CurMap)[p];
            }
        }
        MeanDataPola[kmap]/=(double)NbrTotMask[kmap];
        for (int p=0; p<Map_TQU.Npix();p++) {
            if((*CurMask)[p]>0) 
                    SigDataPola[kmap]+=((*CurMap)[p]-MeanDataPola[kmap])*
                                              ((*CurMap)[p]-MeanDataPola[kmap]);
        }
        SigDataPola[kmap]=sqrt(SigDataPola[kmap]/(double) NbrTotMask[kmap]);
        FSky[kmap]= (double) CurMap->Npix()/NbrTotMask[kmap];
        NbrOverallMask+= NbrTotMask[kmap];
        if (Verbose == True) 
            std::cout << "Percentage of good pixels in Map "<<kmap<<" : "<<
                1.0/FSky[kmap]*100<<" MinData = " << MinDataPola[kmap] <<
                                " MaxData = " << MaxDataPola[kmap] <<std::endl;
    }
    if((MCA == True)||(ForceMonopDip)){
        Use_ZeroMeanConstraint = False;
    }
    
    // We start with a zero mean image
    if (Use_ZeroMeanConstraint == True){
        for (int kmap=0; kmap<3;++kmap) {
            Hdmap* CurMask=Mask_TQU.get_map(kmap);
            Hdmap* CurMap=Map_TQU.get_map(kmap);
            for (int p=0; p<  Map_TQU.Npix(); p++) 
                if ((*CurMask)[p] > 0) (*CurMap)[p] -= MeanDataPola[kmap];
        }
    }
    /*if (Pos == True) {
        for (int kmap=0; kmap<3;++kmap) {
            Hdmap CurMask=Mask_TQU.get_map(kmap);
            Hdmap CurMap=Map_TQU.get_map(kmap);
            for (int p=0; p<  Map_TQU.Npix(); p++) 
                if (CurMask[p] > 0) CurMap[p] -= MinDataPola[kmap];
        }
    }*/
    Nside = Map_TQU.get_nside();
    Lmax_Tot = mrs_get_lmax (Lmax,Nside, ZeroPadding);
    if (Verbose == True) 
        std::cout<<"==> Use Lmax="<<Lmax<<", Lmax TOT="<<Lmax_Tot<<std::endl;
    
    #ifdef MRS_SPEC_DBG
        std::cout<<"Allocate Data"<<std::endl;
    #endif

    ALMPola.alloc(Nside, Lmax_Tot,AlmFast);
    ALMPola.set_flag_NormALM(True);
    //ALM.Norm = (IsotropyCst == True) ? True: False;
    //ALM1.set_beam_eff(Lmax, Lmax_Tot);
    //ALM1.UseBeamEff = False;
    ALM_Rec.alloc(Nside, Lmax_Tot,AlmFast);
    ALMPola.set_flag_NormALM(ALMPola.flag_NormALM());
    Result_TQU.pola_set_nside((int) Nside, False, DEF_MRS_ORDERING);
    Result_TQU.pola_fill(0.);
    PolaResi.pola_set_nside((int) Nside, False, DEF_MRS_ORDERING);
    PolaResi.pola_fill(0.);

    
/*    if (MCA == True){
        NbrMCA = 2;
        MCA_ALM.alloc(Nside, Lmax_Tot);
        MCA_ALM.Norm = ALM.Norm;
        if (MCA_Alm_Lmax ==0) MCA_Alm_Lmax = Lmax_Tot;
        MCA_TabResult = new Hdmap[3];
        for (int i=0; i < NbrMCA; i++)
        {
            MCA_TabResult[i].SetNside ((int) Nside,
                                    (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
            MCA_TabResult[i].fill(0.);
        }
        // WT init
        if ( MCA_WT_NbrScale <= 1)  MCA_WT_NbrScale=4;
        if (DEF_MRS_ORDERING != RING) 
                        MCA_WT.wt_alloc(Nside, MCA_WT_NbrScale, Lmax_Tot, true);
        else  MCA_WT.wt_alloc(Nside, MCA_WT_NbrScale, Lmax_Tot);
        if (MCA_WT_NbrScale <=0) MCA_WT_NbrScale = MCA_WT.nscale();
        if (Verbose ==True) std::cout<<"MCA WT NbrScale="<<MCA_WT.nscale()
                                                                    <<std::endl;
        // Clean CMB first
        // Bool MCA_InitCleabCMB=False;
        if (MCA_InitCleanCMB == True){
            C_UWT WT;
            int Ns = log((double) Nside) / log((double) 2.);
            MCA_TabResult[1] = Map;
            WT.wt_alloc(Nside, Ns, Lmax_Tot, DEF_MRS_ORDERING);
            float Sig=1.;
            WT.hard_thresholding(MCA_TabResult[1],(float) 5., Sig, true, true);
            (MCA_TabResult[1]).write("xx_wt_0.fits");
            for (int p=0; p < Map.Npix(); p++) 
            {
                if (ABS( (MCA_TabResult[2])[p] ) > 0) MaskMap[p]=0;
            }             
        }
    } */
    for(int kmap=0;kmap<3;++kmap)
        if (SigmaNoise[kmap] == 0) SigmaNoise[kmap] = 1.;
    MinSigmaNoise = SigmaNoise[0];
    MinSigmaNoise=MIN(MinSigmaNoise,MIN(SigmaNoise[1], SigmaNoise[2]));
    // else if (ALM.Norm == false) SigmaNoise /= ALM.NormVal;
    
    // Find the minimum value in RMSMap in the valid part of the data
    
    
  // mrs_alloc_powspec(PowSpecSignal, Lmax_Tot);
    weight_T.alloc( 2*Nside );
    if( AlmFast == false ){
        #ifdef MRS_SPEC_DBG
            std:cout<<"SET WEIGHTS"<<std::endl;
        #endif
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[1024];
        sprintf(FN, "%s/data/", HealpixFN);
        read_weight_ring(FN, Nside, weight_T );
        for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
    } else weight_T.fill(1);

#ifdef MRS_SPEC_DBG
    std::cout << "ALMPola.Norm = " << ALMPola.flag_NormALM()<< std::endl;
    std::cout << "ALMPola.NormVal = " << ALMPola.get_normval() << std::endl;
    std::cout << "MinSigmaNoise = " << MinSigmaNoise << std::endl;
    std::cout << "Lmax_Tot  = " << Lmax_Tot << std::endl;
    for(int kmap=0;kmap<3;kmap++){
        std::cout << "ALM for map " << kmap << std::endl;
        std::cout << "MeanData  = " << MeanDataPola[kmap] << std::endl;
        std::cout << "MinData  = " << MinDataPola[kmap] << std::endl;
        std::cout << "MaxData  = " << MaxDataPola[kmap] << std::endl;
        std::cout<<std::endl;
    }
#endif
    
}

/*********************************************************************/
void Pola_VectorField_Inpainting::init_wavelet_packet(){
    ALM_Band.alloc(Nside, Lmax_Tot);
    ALM_Band.Norm = ALMPola.flag_NormALM();
    ALM_Band.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Band.UseBeamEff = UseBeam;
    
    if (All_WP_Band == False) {
        get_wp_meyer_filter(Nside,WP_H,WP_G,WP_W,WP_WPFilter,Lmax_Tot);
        #ifdef MRS_SPEC_DBG
            fits_write_fltarr("Meyer_WP_bands.fits", WP_WPFilter);
        #endif
    } else {
        get_planck_wp_meyer_filter(WP_H,WP_G,WP_W,WP_WPFilter,Lmax,Lmax_Tot);
        #ifdef MRS_SPEC_DBG
            fits_write_fltarr("Planck_Meyer_WP_bands.fits", WP_WPFilter);
        #endif
    }
    NbrWP_Band = WP_H.ny();
    if (Verbose == True) cout << " NbrWP_Band = " << NbrWP_Band << endl;
    Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    
    if (WP_H.nx() != Lmax_Tot+1){
        cout << "Error: bad dimension of the H filter " << endl;
        exit(-1);
    }
}

/*********************************************************************/
void Pola_VectorField_Inpainting::read_spectra(char *Name,int lLmax){
    char *NameFits=fitsname(Name);
    std::string infile=std::string(NameFits);
    free(NameFits);
    PowSpec powspec;
    fitshandle out;
    out.open(infile);
    out.goto_hdu(2);
    int ncols_file=out.ncols();
    long long nrows_file=out.nrows();
    out.close();
    std::cout<<"USE "<<ncols_file<<" FIXED INPUT SPECTRA"<<std::endl;
    std::cout<<"LMAX= "<<nrows_file-1 <<std::endl;

    read_powspec_from_fits(infile,powspec,ncols_file, lLmax);

    std::cout<<"SET PSPEC: NSPEC="<<powspec.Num_specs()<<std::endl;
    std::cout<<"SET PSPEC: Lmax="<<powspec.Lmax()<<std::endl;

    FixedSpectraArray = Healpix_PowSpec(ncols_file, lLmax);
    arr<double> arr_temp;
    if(ncols_file>=1) {
        arr_temp=powspec.tt();
        std::cout<<"TT Size="<<arr_temp.size()<<" vs "<<powspec.tt().size()<<std::endl;
        int test=FixedSpectraArray.Set_TT(arr_temp);
        std::cout<<" test="<<test<<std::endl;
    }
    if(ncols_file>=4) {
        arr_temp=powspec.gg();
        FixedSpectraArray.Set_EE(arr_temp);
        arr_temp=powspec.cc();
        FixedSpectraArray.Set_BB(arr_temp);
        arr_temp=powspec.tg();
        FixedSpectraArray.Set_TE(arr_temp);
    }
    if(ncols_file>=6){
        arr_temp=powspec.tc();
        FixedSpectraArray.Set_TB(arr_temp);
        arr_temp= powspec.gc();
        FixedSpectraArray.Set_EB(arr_temp);
    }
    #ifdef MRS_SPEC_DBG
        FixedSpectraArray.write("test_spectra_read.fits");
    #endif
}

/*********************************************************************/
void Pola_VectorField_Inpainting::getSpectraEigenVectors(PolaAlmR &A1){
    Healpix_PowSpec* spec;
    if(ProjAlm!=NULL) delete ProjAlm;
    ProjAlm = new dblarray(3,3,Lmax_Tot+1);
    ProjAlm->init(0.0);
    MatOper* matOp=new MatOper();
    if(RestrictSpectraSVD){
        if(FixedSpectra) spec= &FixedSpectraArray;
        else {
            std::cout<<"ESTIMATE SPECTRA"<<std::endl;
            spec=new Healpix_PowSpec();
            A1.pola_alm2powspec(*spec);//TT, EE, BB and TE spectra
        }
        dblarray TempTE(2,2),U,S,Vt;//
        for(int l=0; l <= Lmax_Tot; l++){
            TempTE(0,0)= spec->tt(l);
            TempTE(1,1)= spec->ee(l);
            TempTE(0,1)= spec->te(l);
            TempTE(1,0)= spec->te(l);
            matOp->svd(TempTE,U,S,Vt);
            for(int kx=0;kx<2;++kx)
                for(int ky=0;ky<2;++ky) (*ProjAlm)(kx,ky,l)=U(kx,ky);
            (*ProjAlm)(2,2,l)=1.0;
        }
    } else {
        if(FixedSpectra) spec= &FixedSpectraArray;
        else {
            std::cout<<"ESTIMATE FULL SPECTRA"<<std::endl;
            spec=new Healpix_PowSpec();
            A1.pola_alm2powspec_all(*spec);
        }
        dblarray TempTEB(3,3),U,S,Vt,Ut;//
        for(int l=0; l <= Lmax_Tot; l++){
            TempTEB(0,0)= spec->tt(l);
            TempTEB(1,1)= spec->ee(l);
            TempTEB(0,1)= spec->te(l);
            TempTEB(1,0)= spec->te(l);
            TempTEB(2,2)= spec->bb(l);
            TempTEB(2,0)= spec->tb(l);
            TempTEB(0,2)= spec->tb(l);
            TempTEB(2,1)= spec->eb(l);
            TempTEB(1,2)= spec->eb(l);
            matOp->svd(TempTEB,U,S,Vt);
            matOp->transpose(U,Ut);
            for(int kx=0;kx<3;++kx)
                for(int ky=0;ky<3;++ky) (*ProjAlm)(kx,ky,l)=Ut(kx,ky);
        }
    }
    if(!FixedSpectra) delete spec;
    #ifdef MRS_SPEC_DBG
        char TempName[1024];
        sprintf(TempName,"ProjAlmMatrix.fits");
        fits_write_dblarr(TempName,*ProjAlm);
    #endif
    delete matOp;
}

/*********************************************************************/
void Pola_VectorField_Inpainting::projectToEigenVectors(PolaAlmR &A1,bool back){
    Healpix_PowSpec spec;
    //Reference to Stokes ALM
    CAlmR* almT=A1.get_alm(0);
    CAlmR* almP1=A1.get_alm(1);
    CAlmR* almP2=A1.get_alm(2);
    //Structure to store T/P alms for a given l
    dblarray EigenLArr(3,3);
    double *TempLArrBuf, *EigenBuf, *EigenProjBuf,*EigenLArrBuff;
    double Temp=0.;
    MatOper* matOp=new MatOper();
    EigenBuf =ProjAlm->buffer();
    #ifdef MRS_SPEC_DBG
        char TempName[1024];
        sprintf(TempName,"AlmArrayToProj.fits");
        A1.write_array(TempName);
    #endif
    for(int l=0; l <= Lmax_Tot;++l){
        EigenLArrBuff= EigenLArr.buffer();
        for(int kk=0;kk<9;++kk,++EigenBuf,++EigenLArrBuff) 
                                                       *EigenLArrBuff=*EigenBuf;
        if(back){
            double TempValue;
            for(int kx=0;kx<3;++kx){
                for(int ky=kx+1;ky<3;++ky){
                    TempValue =EigenLArr(kx,ky);
                    EigenLArr(kx,ky)=EigenLArr(ky,kx);
                    EigenLArr(ky,kx)= TempValue;
                }
            }
        }
        #ifdef MRS_SPEC_DBG
            if(l==2) {
                std::cout<<"PROJMAT="<<std::endl<<std::flush;
                EigenLArr.display(9);
            }
        #endif
        //Real part
        dblarray TempLArr(l+1,3);
        dblarray EigenProj(l+1,3);
        TempLArrBuf= TempLArr.buffer();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almT)(l,m).real();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almP1)(l,m).real();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almP2)(l,m).real();
        #ifdef MRS_SPEC_DBG
            if(l==2){
                std::cout<<"TEMPLARRBUFF.RE="<<std::endl;
                TempLArrBuf= TempLArr.buffer();
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
            }
        #endif
        matOp->mat_mult(EigenLArr, TempLArr, EigenProj);
        #ifdef MRS_SPEC_DBG
            if(l==2) {
                std::cout<<"PROJMAT2="<<std::endl;
                EigenLArr.display(9);
                std::cout<<"TempLArr.Re="<<std::endl;
                TempLArr.display(9);
                std::cout<<"EigenProj.Re="<<std::endl;
                EigenProj.display(9);
            }
        #endif
        EigenProjBuf= EigenProj.buffer();
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almT)(l,m) =  xcomplex<REAL>(*EigenProjBuf,(*almT)(l,m).imag());
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almP1)(l,m) = xcomplex<REAL>(*EigenProjBuf,(*almP1)(l,m).imag());
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almP2)(l,m) = xcomplex<REAL>(*EigenProjBuf,(*almP2)(l,m).imag());
        #ifdef MRS_SPEC_DBG
            if(l==2){
                std::cout<<"ALM.RE="<<std::endl;
                EigenProjBuf= EigenProj.buffer();
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
            }
        #endif
        //Imaginary part
        TempLArrBuf= TempLArr.buffer();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almT)(l,m).imag();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almP1)(l,m).imag();
        for(int m=0; m <= l; ++m,++TempLArrBuf) *TempLArrBuf= (*almP2)(l,m).imag();
        #ifdef MRS_SPEC_DBG
            if(l==2){
                std::cout<<"TEMPLARRBUFF.IM="<<std::endl;
                TempLArrBuf= TempLArr.buffer();
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++TempLArrBuf) 
                                                std::cout<<*TempLArrBuf<<" ";
                std::cout<<std::endl;
            }
        #endif
        matOp->mat_mult(EigenLArr, TempLArr, EigenProj);
        EigenProjBuf= EigenProj.buffer();
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almT)(l,m) = xcomplex<REAL>((*almT)(l,m).real(), *EigenProjBuf);
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almP1)(l,m) = xcomplex<REAL>((*almP1)(l,m).real(), *EigenProjBuf);
        for(int m=0;m<=l;++m,++EigenProjBuf) (*almP2)(l,m) = xcomplex<REAL>((*almP2)(l,m).real(), *EigenProjBuf);
        #ifdef MRS_SPEC_DBG
            if(l==2){
                std::cout<<"ALM.IM="<<std::endl;
                EigenProjBuf= EigenProj.buffer();
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
                for(int m=0; m <= l; ++m,++EigenProjBuf) 
                                                std::cout<<*EigenProjBuf<<" ";
                std::cout<<std::endl;
            }
        #endif
    }
    if(ForceMonopDip){
        (*almT)(0,0) = xcomplex<REAL>(0.,0.);
        (*almT)(1,0) = xcomplex<REAL>(0.,0.);
        (*almT)(1,1) = xcomplex<REAL>(0.,0.);
        (*almT)(1,1) = xcomplex<REAL>(0.,0.);
        // JLS:  HEALPIX 3.60 modif
        //        (*almT)(0,0).re=0.;
        //        (*almT)(1,0).re=0.;
        //        (*almT)(1,1).re=0.;
        //        (*almT)(1,1).im=0.;
    }
    delete matOp;
    #ifdef MRS_SPEC_DBG
        sprintf(TempName,"AlmArrayProjected.fits");
        A1.write_array(TempName);
    #endif

}

/*********************************************************************/
void Pola_VectorField_Inpainting::forward_transform(PolaHdmap &D1,PolaAlmR &A1){
    if(EBTrans) A1.pola_alm_trans(D1);
    else {
        //All Maps are assumed in ring format
        //Map is in TQU, ALM in TQU	
        int maxL=A1.get_lmax();
        {
            Hdmap TMap= D1.get_map_T();
            CAlmR* AlmT= A1.get_alm(0);
            double avg = TMap.average();
            TMap.Add(-avg);
            #ifdef MRS_SPEC_DBG
                std::cout<<"Process T map"<<std::endl;
                std::cout<<"avg="<<avg<<std::endl;
                std::cout<<"maxL ="<<maxL <<std::endl;
                std::cout<<"niter ="<<A1.get_niter() <<std::endl;
                std::cout<<"Norm ="<<A1. get_normval() <<std::endl;
            #endif
            map2alm_iter(TMap, *AlmT, A1.get_niter(), weight_T);
            (*AlmT)(0,0) += avg*sqrt(fourpi);
            if(A1.get_normval() == true ){
                for(int l=0; l <= maxL; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmT)(l,m) *= A1.get_normval();
                    }
                }
            }
        }
        {
            Hdmap QMap= D1.get_map_Q();
            CAlmR* AlmQ= A1.get_alm(1);
            map2alm_iter(QMap, *AlmQ, A1.get_niter(), weight_T);
            if(  A1.get_normval()  == true ){
                for(int l=0; l <= maxL; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmQ)(l,m) *= A1.get_normval();
                    }
                }
            }
        }
        {
            Hdmap UMap= D1.get_map_U();
            CAlmR* AlmU= A1.get_alm(2);
            map2alm_iter(UMap, *AlmU, A1.get_niter(), weight_T);
            if( A1.get_normval() == true ){
                for(int l=0; l <= maxL; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmU)(l,m) *= A1.get_normval();
                    }
                }
            }
        }
    }
    if(ForceMonopDip){
        CAlmR* almT= A1.get_alm(0);
        (*almT)(0,0) = xcomplex<REAL>(0.,0.);
         (*almT)(1,0) = xcomplex<REAL>(0.,0.);
         (*almT)(1,1) = xcomplex<REAL>(0.,0.);
         (*almT)(1,1) = xcomplex<REAL>(0.,0.);
// JLS:  HEALPIX 3.60 modif
//        (*almT)(0,0).re=0.;
//        (*almT)(1,0).re=0.;
//        (*almT)(1,1).re=0.;
//        (*almT)(1,1).im=0.;
    }
}

/*********************************************************************/
void Pola_VectorField_Inpainting::backward_transform(PolaAlmR &A1,PolaHdmap &D1)
{
    if(EBTrans) A1.pola_alm_rec(D1);
    else {
        //All Maps are assumed in ring format
        //Map is in TQU, ALM in TQU	
        {
            Hdmap* TMap= D1.get_map(0);
            CAlmR* AlmT= A1.get_alm(0);
            if( A1.get_normval() == true ){
                for(int l=0; l <= Lmax_Tot; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmT)(l,m) /= A1.get_normval();
                    }
                }
            }
            double avg = (*AlmT)(0,0).real()/sqrt(fourpi);
            (*AlmT)(0,0) = 0;
            alm2map(*AlmT,*TMap);
            TMap->Add(avg);
        }
        {
            Hdmap* QMap= D1.get_map(1);
            CAlmR* AlmQ= A1.get_alm(1);
            if( A1.get_normval() == true ){
                for(int l=0; l <= Lmax_Tot; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmQ)(l,m) /= A1.get_normval();
                    }
                }
            }
            alm2map(*AlmQ,*QMap);
        }
        {
            Hdmap* UMap= D1.get_map(2);
            CAlmR* AlmU= A1.get_alm(2);
            if( A1.get_normval() == true ){
                for(int l=0; l <= Lmax_Tot; l++){
                    for(int m=0; m <= l; m++) {
                        (*AlmU)(l,m) /= A1.get_normval();
                    }
                }
            }
            alm2map(*AlmU,*UMap);
        }
    }
}

/*********************************************************************/
double Pola_VectorField_Inpainting::get_residual()
{
    double Err = 0.;
    for (int kmap=0; kmap< 3; ++ kmap) {
        Hdmap* CurResi=PolaResi.get_map(kmap);
        Hdmap* CurMap=Map_TQU.get_map(kmap);
        Hdmap* CurMask=Mask_TQU.get_map(kmap);
        Hdmap* CurRes= Result_TQU.get_map(kmap);
        for (int p=0; p < PolaResi.Npix(); p++) {
            (*CurResi)[p]= (*CurMap)[p]-(*CurRes)[p];
            if (MCA != True) (*CurResi)[p] *= (*CurMask)[p];
            /*  else{
                Result[p] = 0;
                for (int i=0; i < NbrMCA; i++) Result[p]+=(MCA_TabResult[i])[p];
            }*/
            Err += (*CurResi)[p]* (*CurResi)[p]/SigmaNoise[kmap];
        }
    }
    Err = sqrt(Err /(double) NbrOverallMask);
    return Err;
}

/*********************************************************************/
void Pola_VectorField_Inpainting::wp_band_constraint(Hdmap &Band, 
                                        Hdmap & MaskMap, int NbrTotMask, int b){
    double EnergyInMap = 0.;
    for (int p=0; p < MaskMap.Npix(); p++)  
        if (MaskMap[p] != 0) EnergyInMap += Band[p]*Band[p];
    double Sig_InMap  = sqrt(EnergyInMap/(double) (NbrTotMask));
    
    // Compute the expected power in the band
    double EnergyBand=0.;
    double Np = (double) MaskMap.Npix();
    // if (InputPowSpec == True)  EnergyBand = Np * TabPsi2(b);
    // else 
    EnergyBand = EnergyInMap * Np / (double) NbrTotMask;
    double SigBand = (EnergyBand > 0) ? sqrt(EnergyBand / Np): 0;
    
    // double EnergyBand = Map.Npix() * SigBand*SigBand;
    
    double EnergyInMask = 0.;
    for (int p=0; p < MaskMap.Npix(); p++) {
        // float Level = GaussCst_ProjectCst_Lambda * Sig_InMap;
        if (MaskMap[p] == 0) EnergyInMask += Band[p]*Band[p];
    }
    double Sig_Mask = sqrt( EnergyInMask / (double) (Np - NbrTotMask));
    
    // Band renormalization
    float Coeff =  1.;
    double ExpectedEnergyInMask = EnergyBand - EnergyInMap;
    if (ExpectedEnergyInMask <= 0) //Can never occur
        ExpectedEnergyInMask = EnergyBand * (Np -  NbrTotMask) / Np;
    double E = EnergyInMap  * (Np -  NbrTotMask) / (double) NbrTotMask;
    //std::cout<<"      TabPsi2="<<sqrt(TabPsi2(b))<<", ExpectedSigmaInMask="
    //     << sqrt(ExpectedEnergyInMask/(Np-NbrTotMask))<<" "<<
    //        sqrt(E/(Np-NbrTotMask))<<" "<<Sig_Mask<<std::endl;
    if (EnergyInMask > 0) 
        Coeff =  sqrt(ExpectedEnergyInMask  / EnergyInMask);
    if (Coeff > 5) Coeff = 5.;
    
    //if((Sig_Mask>0)&&(Sig_InMap>0)) Coeff=Sig_InMap/Sig_Mask;
    //if (Verbose == True) std::cout<<"        WP:Band "<<b+1<<",   Coeff "<<
    //          Coeff <<" Sig_Mask="<<Sig_Mask<<", Sig_InMap="<<Sig_InMap<<
    //              "   Npix="<<Map.Npix()<<std::endl;
    
    for (int p=0; p < MaskMap.Npix(); p++) 
                                        if (MaskMap[p] == 0)  Band[p] *= Coeff;
}

/*********************************************************************/
void Pola_VectorField_Inpainting::wp_constraint(CAlmR &ALM,Hdmap &Resi,
                   Hdmap &Result,Hdmap &MaskMap,int NbrTotMask,int ALM_FIRST_L){
    Resi.init();
    for (int b=0; b <= NbrWP_Band; b++){
        int LMin=(b != NbrWP_Band)?WP_W(NbrWP_Band-1-b,0):WP_W(0,0);
        int LMax=(b != NbrWP_Band)?WP_W(NbrWP_Band-1-b,1):WP_W(0,1);
        //if (Verbose == True) std::cout<<"        WP:Band "<<b+1<<", Lmin="
        //                                  <<LMin<<", Lmax="<<LMax<<std::endl;
        // If ALM_FIRST_L < first non zero coef of the filter then band is zero.
        if(LMin > ALM_FIRST_L) Band.init();
        else {
            if (b == NbrWP_Band) Band = Result;
            else {
                ALM_Band.Set(LMax,LMax);
                ALM_Band.SetToZero();
                // Compute solution at given resolution (convolution with H)
                int FL = MIN(ALM_Band.Lmax(),WP_H.nx()-1);
                if(FL > ALM.Lmax()) FL = ALM.Lmax();
                for (int l=0;l<=FL;++l)
                    for (int m=0;m<=l;++m){
                        ALM_Band(l,m) = ALM(l,m) * (REAL) WP_H(l,NbrWP_Band-1-b);
                    }
                ALM_Band.alm_rec(Band);
                #ifdef MRS_SPEC_DBG
                    char FN[512];
                    sprintf(FN, "xx_r%d.fits", b+1);
                    Result.write(FN);
                    sprintf(FN, "xx_l%d.fits", b+1);
                    Band.write(FN);
                #endif
                // Compute coeffs/standard deviations in and out the mask
                for (int p=0; p < Result.Npix(); p++) {
                    double Val = Band[p];
                    Band[p] = Result[p] - Val;
                    Result[p] = Val;
                }
            }// endelse if (b == NbrWP_Band)
            #ifdef MRS_SPEC_DBG
                char FN[512];
                sprintf(FN, "xx_b%d.fits", b+1);
                Band.write(FN);
            #endif
            // Wavelet packet constraint on the band b
            wp_band_constraint(Band,MaskMap,NbrTotMask,b);
            // coadd the band to new new solution
            Resi+=Band;
        } // else if  (LMin > ALM_FIRST_L) 
    } // endfor b
    Result=Resi;
}


/*********************************************************************/
void Pola_VectorField_Inpainting::exit_inpainting(CAlmR &ALM,Hdmap &Map,
                   Hdmap &MaskMap,Hdmap &Result,double MinData,double MeanData){
    double ValAdd=0.;
    if(Pos == True) ValAdd+=MinData;
    if(Use_ZeroMeanConstraint==True) ValAdd+=MeanData;
    
    if (ValAdd != 0){
        for (int p=0;p<Map.Npix();++p){
            if(MaskMap[p]>0) Map[p]+=ValAdd;
            Result[p]+=ValAdd;
        }
    }
    // We reinsert the correct pixel values 
    if (EqualInMask == True)
       for (int p=0; p < Map.Npix(); p++) if(MaskMap[p]>0) Result[p]=Map[p];
    
    if (ALM.Norm == True){
        double CoefN = ALM.NormVal;
        for (int l=0; l <= ALM.Lmax(); l++)
            for (int m=0; m <= l; m++){
                ALM(l,m) /= CoefN;
            }
    }
}


/*********************************************************************/
void Pola_VectorField_Inpainting::analyse_alm_inpainting(){
    double Sig_InMap, Sig_Mask;
    Bool Debug = False;

    #ifdef MRS_SPEC_DBG
        std::cout<<"Init inpainting"<<std::endl;
    #endif
    init_inpainting();
    if (Use_WP_Constraint==True){
        #ifdef MRS_SPEC_DBG
            std::cout<<"Init Wavelet Packets"<<std::endl;
        #endif
         init_wavelet_packet();
    }
    // if (Verbose == True)  Map.info((char *) "INIT");
    // Map.write("xx_tt.fits");
    
    #ifdef MRS_SPEC_DBG
        for(int kmap=0;kmap<3;++kmap) {
            char MapName[64];
            sprintf(MapName,"Map%d",kmap);
            Map_TQU.get_map(kmap)->info(MapName);
        }
    #endif
    forward_transform(Map_TQU, ALMPola);
    #ifdef MRS_SPEC_DBG
        ALMPola.write_array((char *) "xx_init_alm.fits");
    #endif
    
    double MaxAlm = ALMPola.max_absalm_TEB();
    // IDL: NormVal =  sqrt(nel / (4.*!DPI))
    MaxAlm *= 0.9999999;
    // Normally ALM.Norm should be set to True
    if (ALMPola.flag_NormALM()==False) 
         for(int kmap=1;kmap<3;++kmap) SigmaNoise[kmap]/=ALMPola.get_normval();
    
    if (ALMPola.flag_NormALM() == True) 
                            std::cout<<" Alm norm: max =  "<<MaxAlm<<std::endl;
    else std::cout<<" Alm nonorm: max =  " <<MaxAlm<<std::endl;
      
    double FirstSoftThreshold[3],LastSoftThreshold[3],DeltaThreshold[3],
                                                                    Lambda[3];
    if(!ThresholdPerChannel){
        for(int kmap=0;kmap<3;++kmap){
            FirstSoftThreshold[kmap]=MaxAlm / SigmaNoise[kmap] * sqrt(2.);
                          // sqrt(2.) is here, to be equivalent to the IDL code.
            LastSoftThreshold[kmap]=0.;
            DeltaThreshold[kmap]=
                               FirstSoftThreshold[kmap]-LastSoftThreshold[kmap];
            Lambda[kmap]=FirstSoftThreshold[kmap];
            // this has no effect since we redivide latter by the same value.
            #ifdef MRS_SPEC_DBG
                std::cout<<"==>Max Alm="<<MaxAlm<<" "<<(MaxAlm*SigmaNoise[kmap])
                    <<", N="<<ALMPola.get_normval()<<", First Threshold="<<
                                            FirstSoftThreshold[kmap]<<std::endl;
            #endif
        }
    } else {
        for(int kmap=0;kmap<3;++kmap){
            FirstSoftThreshold[kmap]=(ALMPola.get_alm(kmap))->max_absalm()
                                                  / SigmaNoise[kmap] * sqrt(2.);
                          // sqrt(2.) is here, to be equivalent to the IDL code.
            LastSoftThreshold[kmap]=0.;
            DeltaThreshold[kmap]=
                               FirstSoftThreshold[kmap]-LastSoftThreshold[kmap];
            Lambda[kmap]=FirstSoftThreshold[kmap];
            // this has no effect since we redivide latter by the same value.
            #ifdef MRS_SPEC_DBG
                 std::cout<<"==>Max Alm="<<(ALMPola.get_alm(kmap))->max_absalm()
                   <<" "<<MaxAlm*SigmaNoise[kmap]<<", N="<<ALMPola.get_normval()
                   <<", First Threshold="<<FirstSoftThreshold[kmap]<<std::endl;
                 std::cout<<"Lmax="<<(ALMPola.get_alm(kmap))->Lmax()<<std::endl;
            #endif
        }
    }

    double Err;
    int MaxNonZeroL[3];
    for(int kx=0;kx<3;++kx) MaxNonZeroL[kx]=0; 
    int ALM_FIRST_L = (Acceleration == True) ? 100: Lmax_Tot;
    if (IsotropyCst == True) ALM_FIRST_L = Lmax_Tot;
    
    // Use for isotropic constrain only
    float GaussCst_LambdaFirst=0.5;
    float GaussCst_StepLambda=(GaussCst_ProjectCst_Lambda-GaussCst_LambdaFirst)
                                                               /float(Max_Iter);
    float GaussCst_Lambda = GaussCst_LambdaFirst;
    
    // fits_write_fltarr((char *)"xxpsi.fits",TabPsi2);
    // wp_trans(Map,Band,ALM,ALM_Band,WP_W,WP_H,NbrWP_Band,ALM_FIRST_L,TabPsi2);
    // exit(0);

    #ifdef MRS_SPEC_DBG
        std::cout<<"Lmax_Tot ="<<Lmax_Tot<<std::endl;
        std::cout<<"ALM_FIRST_L="<<ALM_FIRST_L<<std::endl;
    #endif
    // MAIN ITER 
    #ifdef MRS_SPEC_DBG
        char NN[256];
        sprintf(NN, "xx_input.fits");
        Map_TQU.write(NN);
        sprintf(NN, "xx_res_init.fits");
        Result_TQU.write(NN);
    #endif
    
    #ifdef MRS_SPEC_DBG
        std::cout<<"Start Main Loop"<<std::endl;
        std::cout<<"Threshold Residual ? "<<ThresholdResi<<std::endl;
    #endif
    //MAIN LOOP
    for (int i=0; i < Max_Iter; i++){
        #ifdef MRS_SPEC_DBG
            std::cout<<"Start Iteration "<<i<<std::endl;
            std::cout<<"Compute Residual "<<std::endl;
        #endif
        Err = get_residual();
        
        if (ThresholdResi == False){
            for(int kmap=0;kmap<3;++kmap){
                Hdmap* CurResi=PolaResi.get_map(kmap);
                Hdmap* CurRes= Result_TQU.get_map(kmap);
                *CurResi+= *CurRes;//Add residual to current result
            }
            #ifdef MRS_SPEC_DBG
                char NN[256];
                sprintf(NN, "xx_resi_it_%d.fits", i+1);
                if (i >= 0) PolaResi.write(NN);
            #endif
        }
        int CptAlm=0;
        
        if (Acceleration == True) ALMPola.alloc(ALMPola.get_nside(),ALM_FIRST_L,
                                                                       AlmFast);
        if(Verbose) std::cout<<"Forward Transform Data "<<std::endl;
        forward_transform(PolaResi, ALMPola);//Get current Alms
        //Force Monopole and Dipole to 0
        if(ForceMonopDip){
            CAlmR* almT=ALMPola.get_alm(0);
            (*almT)(0,0) = xcomplex<REAL>(0.,0.);
             (*almT)(1,0) = xcomplex<REAL>(0.,0.);
             (*almT)(1,1) = xcomplex<REAL>(0.,0.);
             (*almT)(1,1) = xcomplex<REAL>(0.,0.);
            // JLS modif HEALPIX 3.60
            // (*almT)(0,0).re=0.;
            // (*almT)(1,0).re=0.;
            // (*almT)(1,1).re=0.;
            // (*almT)(1,1).im=0.;
        }
        // MaxAlm =ALM.max_absalm() / SigmaNoise;
        // cout << " Max ABS(ALM) Resi  " << MaxAlm  << endl;
        #ifdef MRS_SPEC_DBG
            std::cout<<"Compute Threshold "<<std::endl;
        #endif
        if(i>0){
            for(int kmap=0;kmap<3;++kmap){
                Lambda[kmap]=LastSoftThreshold[kmap]+DeltaThreshold[kmap]*
                                       (1.-erf(2.8*((float)i/(float)Max_Iter)));
                                       // 4 versus 2.8  in IDL code
                if((Lambda[kmap]<LastSoftThreshold[kmap])||(i==Max_Iter-1))
                                        Lambda[kmap] = LastSoftThreshold[kmap];
                #ifdef MRS_SPEC_DBG
                    printf(" Lambda[%d] = %5.4f\n", kmap,(float) Lambda[kmap]);
                    Hdmap *CurResi= PolaResi.get_map(kmap);
                    char NN[64];
                    sprintf(NN, "Resi%d",i);
                    CurResi->info(NN);
                #endif
            }
        } 
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        if (Verbose == True) std::cout<<"Iter "<<i+1<<", Lambda="<<Lambda[0]<<
                          ","<<Lambda[1]<<","<<Lambda[2]<<" Err = "<<Err<<endl;

        if(SpectraSVD){
            #ifdef MRS_SPEC_DBG
                std::cout<<"Use constraints in SVD space "<<std::endl;
                char NN[256];
                sprintf(NN, "xx_almpola_it_%d.fits", i+1);
                ALMPola.write(NN);
            #endif
            getSpectraEigenVectors(ALMPola);
            projectToEigenVectors(ALMPola,false);
        }

        // if (Denoising == True) alm_denoising(ALM, WienerFilter);
        #ifdef MRS_SPEC_DBG
            std::cout<<"Threshold "<<std::endl;
        #endif
        if((HardThreshold==True)||(SoftThreshold==True)){
            // Alm L1 minimization
            for(int kmap=0;kmap<3;++kmap){
                float ThresholdLevel = SigmaNoise[kmap]*Lambda[kmap] / sqrt(2.);
                CAlmR* CurAlm=ALMPola.get_alm(kmap);
                if (Lambda[kmap]>0){
                    if (HardThreshold == True) CptAlm= 
                     (*CurAlm).hard_threshold(ThresholdLevel,MaxNonZeroL[kmap]);
                    else if (SoftThreshold == True) CptAlm=
                     (*CurAlm).soft_threshold(ThresholdLevel/sqrt(2.),
                                                            MaxNonZeroL[kmap]);
                } else  MaxNonZeroL[kmap]= Lmax_Tot;
                #ifdef MRS_SPEC_DBG
                    std::cout<<"Map "<<kmap<<std::endl;
                    std::cout<<"  Number of significant Alm="<<CptAlm <<" MaxL="
                                            <<MaxNonZeroL[kmap] <<std::endl; 
                #endif
            }
            // Update the solution
            if (ThresholdResi == True){
                backward_transform(ALMPola,PolaResi);
                for(int kmap=0;kmap<3;++kmap){
                    Hdmap* CurResi=PolaResi.get_map(kmap);
                    Hdmap* CurRes=Result_TQU.get_map(kmap);
                    for (int p=0; p<  Result_TQU.Npix(); p++) 
                                               (*CurRes)[p]+=Eps* (*CurResi)[p];
                }
                forward_transform(Result_TQU,ALMPola);
            } else backward_transform(ALMPola, Result_TQU);
            int OverallMaxNonZeroL = MAX(MaxNonZeroL[0], MAX(MaxNonZeroL[1],
                                                               MaxNonZeroL[2]));
            if (2*OverallMaxNonZeroL >  ALM_FIRST_L)  
                                            ALM_FIRST_L = 2*OverallMaxNonZeroL;
            ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
            #ifdef MRS_SPEC_DBG
                char NN[256];
                if(HardThreshold==True)
                    sprintf(NN, "xx_hardthr_it_%d.fits", i+1);
                else 
                    sprintf(NN, "xx_softthr_it_%d.fits", i+1);
                Result_TQU.write(NN);
            #endif
        } // end l1 minimization i.e.  ((HardThreshold == True)  || (SoftThreshold == True))

        // Apply constraint per wavelet packet band. We start these constraint only at half the total number of iterations.
        // Use_WP_Constraint = False;
        if((Use_WP_Constraint==True)&&(i>0)&&(i%AccelationLevel==0)){
            if(Verbose==True) 
                std::cout<<"        Wavelet Packets constraint ... "<<std::endl;
            // We compute the WP transform:
            //   ITER: Result contains the high resolution   (c_j)
            //   We compute Band, which contains the low resolution (c_{j+1}
            //   Coef =  Result - Band
            //   Result = Band
            #ifdef MRS_SPEC_DBG
                PolaHdmap Result_TQU_temp;
                backward_transform(ALMPola, Result_TQU_temp);
                char NN[256];
                sprintf(NN, "xx_proj_svd_it_%d.fits", i+1);
                Result_TQU_temp.write(NN);
            #endif
            for(int kmap=0;kmap<3;++kmap){
                Hdmap* CurResi=PolaResi.get_map(kmap);
                Hdmap* CurMask=Mask_TQU.get_map(kmap);
                Hdmap* CurRes= Result_TQU.get_map(kmap);
                CAlmR* CurAlm=ALMPola.get_alm(kmap);
                wp_constraint(*CurAlm,*CurResi,*CurRes,*CurMask,
                                            NbrTotMask[kmap], ALM_FIRST_L);
            }
            #ifdef MRS_SPEC_DBG
                sprintf(NN, "xx_projband_svd_it_%d.fits", i+1);
                Result_TQU.write(NN);
            #endif
        }// endif Use_WP_Constraint
        if(SpectraSVD){
            forward_transform(Result_TQU,ALMPola);
            projectToEigenVectors(ALMPola,true);
            backward_transform(ALMPola,Result_TQU);
        }
        
        // We want a zero mean solution
        if (Use_ZeroMeanConstraint == True){
            if(Verbose==True) std::cout << "ZERO MEAN CONSTRAINT" <<std::endl;
            for(int kmap=0;kmap<3;++kmap){
                Hdmap* CurRes= Result_TQU.get_map(kmap);
                double MeanRes = CurRes->average();
                for (int p=0;p< CurRes->Npix(); p++) (*CurRes)[p]-=MeanRes;
            }
        }
        // We want a positive solution
        if (Pos == True){
            if(Verbose==True) std::cout <<"POSITIVITY CONSTRAINT"<<std::endl;
            for(int kmap=0;kmap<3;++kmap){
                Hdmap* CurRes= Result_TQU.get_map(kmap);
                for (int p=0; p < CurRes->Npix(); p++) 
                                            if ((*CurRes)[p]<0) (*CurRes)[p]=0.;
            }
        }
        #ifdef MRS_SPEC_DBG
            char NN[256];
            sprintf(NN, "xx_res_it_%d.fits", i+1);
            Result_TQU.write(NN);
        #endif
        
        #ifdef MRS_SPEC_DBG
            for(int kmap=0;kmap<3;++kmap){
                Hdmap* CurRes= Result_TQU.get_map(kmap);
                std::cout<<"MAP "<<kmap<<std::endl;
                CurRes->info((char *) "  REC");
            }
        #endif
    }//END MAIN LOOP
    for(int kmap=0;kmap<3;++kmap){
        Hdmap* CurMask=Mask_TQU.get_map(kmap);
        Hdmap* CurRes= Result_TQU.get_map(kmap);
        Hdmap* CurMap=Map_TQU.get_map(kmap);
        CAlmR* CurAlm=ALMPola.get_alm(kmap);
        exit_inpainting(*CurAlm,*CurMap,*CurMask,*CurRes, MinDataPola[kmap],
                                                           MeanDataPola[kmap]);
    }
}

/*********************************************************************/
int main(int argc, char *argv[]){
    int k;
    fitsstruct Header, HD1;
    char Cmd[2048];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    
    Pola_VectorField_Inpainting INP;
    if(FixedSpectra){
        std::cout<<"READ SPECTRA IN "<<Name_InSpectra<<std::endl;
        INP.read_spectra(Name_InSpectra,Lmax);
        std::cout<<"Nspecs= "<<INP.FixedSpectraArray.Num_specs()<<std::endl;
        Lmax=MIN(Lmax,INP.FixedSpectraArray.Lmax());
    }
    if(Verbose == True){ 
        std::cout << "# PARAMETERS: " <<std::endl;
        std::cout << "# File Name TQU = "<<Name_Imag_In_TQU << std::endl;
        std::cout << "# File Name Mask TQU = "<<Name_Mask_TQU << std::endl;
        std::cout << "# File Out Prefix = "<<Name_Imag_Out << std::endl;   
        std::cout << "# Number of iterations = "<<Max_Iter << std::endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << std::endl;
        if (EBMode== False) std::cout << "# EBmode is off " <<std::endl;
        if (Use_WP_Constraint == True) 
                                   std::cout<<"# Use WP constraint "<<std::endl;
        if (SpectraSVD == True) 
                                   std::cout<<"# Use SpectraSVD "<<std::endl;
        if (RestrictSpectraSVD == True) 
                                   std::cout<<"# Only TT/EE/BB/TE "<<std::endl;
                                  
    }

    #ifdef MRS_SPEC_DBG
        std::cout << "DEBUG MODE " <<std::endl;
    #endif
    if(Verbose) std::cout<<"READ MAP"<<std::endl;
    INP.Map_TQU.read(Name_Imag_In_TQU,false,true);
    if(Verbose) std::cout<<"READ MASK"<<std::endl;
    INP.Mask_TQU.read(Name_Mask_TQU,false,true);
    
    if(INP.Mask_TQU.get_nside()!=INP.Map_TQU.get_nside()){
        fprintf(OUTMAN,
        "Error: Maps and Masks have not the same Nside in: %s ...\n", argv[0] );
        exit(-1);
    }
    for(int kmap=0;kmap<3;++kmap){
        Hdmap* CurMask=INP.Mask_TQU.get_map(kmap);
        Hdmap* CurMap=INP.Map_TQU.get_map(kmap);
        for (int p=0; p < INP.Map_TQU.Npix(); p++){
            if ((*CurMask)[p] != 1) {
                (*CurMask)[p] = 0;
                (*CurMap)[p] = 0;
            }
        }
    }
    
    int Nside = INP.Map_TQU.get_nside();
    int Lmax_Tot = mrs_get_lmax (Lmax, Nside, ZeroPadding);
    if (Verbose == True) cout << "# Used Lmax = " <<  Lmax_Tot << endl;
    // mrsp_alloc_powspec(INP.PolPowSpecData, Lmax_Tot);
    INP.All_WP_Band=All_WP_Band;
    INP.Verbose = Verbose;
    INP.HardThreshold=HardThreshold;
    INP.SoftThreshold=SoftThreshold;
    INP.ThresholdResi=ThresholdResi;
    for(int kmap=0;kmap<3;++kmap) INP.SigmaNoise[kmap]=SigmaNoise;
    INP.Use_ZeroMeanConstraint=Use_ZeroMeanConstraint;
    INP.Use_WP_Constraint=Use_WP_Constraint;
    INP.All_WP_Band=All_WP_Band;
    INP.ZeroPadding=ZeroPadding;
    INP.Lmax = Lmax;
    INP.Eps=Eps;
    INP.Acceleration=Acceleration;
    INP.AccelationLevel = AccelationLevel;
    INP.Pos=Pos;
    INP.Max_Iter = Max_Iter;
    INP.EBTrans = EBMode;
    INP.EqualInMask = EqualInMask;
    switch (InpaintMethod){
        case POL_INP_L1_ALM_ANALYSIS_CST_VAR:
            Analysis = True;
            if (AccelationLevel == 0) INP.Use_WP_Constraint = False;
            INP.analyse_alm_inpainting();
            break;
        default: cout << "Error: method not implemented ... " << endl;
            exit(-1);
            break;
    }
    char FN[2048];
    sprintf(FN, "%s_mapTQU.fits", Name_Imag_Out);
    INP.Result_TQU.write(FN);
    if(TEBWrite == True){
        INP.EBTrans=True;
        INP.forward_transform(INP.Result_TQU, INP.ALMPola);
        Healpix_PowSpec pspec;
        INP.ALMPola.pola_alm2powspec(pspec);
        sprintf(FN, "%s_TEBspec.fits", Name_Imag_Out);
        pspec.write(FN);
        sprintf(FN, "%s_alms.fits", Name_Imag_Out);
        INP.ALMPola.write_array(FN);
        INP.EBTrans=False;
        INP.backward_transform(INP.ALMPola,INP.Result_TQU);
        sprintf(FN, "%s_mapTEB.fits", Name_Imag_Out);
        INP.Result_TQU.write(FN);
    }
    exit(0);
}

