/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date: 19/09/08
**    
**    File:  HealpixClass.h
**
*******************************************************************************
**
**    DESCRIPTION  Class definition for inpainting   
**    ----------- 
**                 
******************************************************************************/

#ifndef _INP_H_
#define _INP_H_


#include "SB_Filter.h"
// #include "GMCA.h"
#include "Array.h"
#include "NR.h"
#include "IM_IO.h"
#include <cmath>
#include "MatMask.h"

#define NBR_INP_METHOD 4
enum type_inpainting {INP_L1_ALM_ANALYSIS_CST_VAR, INP_L1_ALM_SYNTHESIS, INP_L1_ALM_SYNT_MASTER, INP_ALM_ISOTROP, INP_L2_ALM};
#define DEF_INP_METHOD INP_L1_ALM_SYNT_MASTER
const char * StringInpaintingMethod (type_inpainting TypeMethod);


void isotropy_constraint(CAlmR &ALM, PowSpec & P, float Threshold, Bool UpdatePowSpec);

#define  DEF_IWIENER_NBR_ITER  40
#define  DEF_ITER_GAUSSIANITY 10

class CInpainting {
    CAlmR ALM_Rec;
    Hdmap Resi;
    PowSpec PowSpecNoise;
     // Wavelet packet variables
    fltarray WP_H, WP_G, WP_W, WP_WPFilter, TabPsi2;
    int NbrWP_Band;
    CAlmR  ALM_Band;
    Hdmap Band;
    void init_wavelet_packet(); // initialization for the wavelet packet decomposition
    void wp_band_constraint(Hdmap &Band, int b);
    void forward_transform(Hdmap & Map1, Hdmap & Map2, CAlmR & ALM1, CAlmR & ALM2, Bool ToEB);

    void init_inpainting();
    void exit_inpainting();
    double get_residual();
    void alm_denoising(CAlmR & A, Bool Wiener=False);
    void alm_powspec_constraint(CAlmR & A);
    void alm_l2_constraint(CAlmR & A, double Lambda);

    // Alm denoising using the Wiener or the Cole filter.
    //          Den = H2 * (PS_Signal  + PS_Noise )
    //   WienerFilter = PS_Signal  / Den
    //   ColeFilter = sqrt(PS_Signal  / Den)
    double Min, Max;
    int NbrTotMask; 
    double SigData;
    double MeanData;
    double FSky;
    double MinData, MaxData;
    int Nside;
    int Lmax_Tot;
    
    public:
    Bool MCA;
    
    type_inpainting InpaintMethod;

    Bool ThresholdResi;   // Threshold the Alm of the residual instead of the Alm of the solution
    Bool Verbose;         // Verbose mode
    int Max_Iter;         // Number of iteration in the inpainting algorithm
    Bool HardThreshold;   //  Apply an Alm iterative hard threshold for the inpainting
    Bool SoftThreshold;   // Apply an Alm iterative soft threshold for the inpainting
    float SigmaNoise;     // Noise standard deviation
    dblarray MatMask;     // Matrix related to the Mask
    dblarray RevMatMask;  // Reverse  Matrix related to the Mask

    // Zero mean constrain on the solution 
    Bool Use_ZeroMeanConstraint; // Apply a zero mean constraint on the solution

    // Isotropy constraint
    Bool IsotropyCst;   // if true, inpainting using an isotropy constraint. Otherwise, inpainting using l1 minimization.
    Bool Use_Realization_Powspec;
    
    // Parameters for WP decomposition constraints
    Bool Use_WP_Constraint;      // Apply a variance constraint in the wavelet packet of the solution
    Bool All_WP_Band;
    Bool ProjectConstraint;      // Add a Gaussianity projection constraint. Apply a wavelet packet decomposition, and set to 
                                 // to zero all coeff larger than  Gauss_Cst_KSigma*SigBand
    Bool GaussianityConstraint;  // Apply a Gaussianity constraint, by soft-thresholding all coef in WP decomposition at each iter.
    Bool L2Constraint;   
    Bool Denoising;              // If this term is true, a denoising a also performed.
    float GaussCst_ProjectCst_Lambda; // using a thresholding equal to GaussCst_ProjectCst_Lambda * Sig_InMap
    double Gauss_Cst_KSigma;

    float  ZeroPadding; // 
    int Lmax;
    float Eps;
                              
    bool Acceleration;
    int AccelationLevel;
    bool Pos;
    bool InputPowSpec; 

    Hdmap Map;     // Map to inpainted
    CAlmR  ALM;  
    Hdmap MaskMap; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    Hdmap Result;  // inpainted map
    Hdmap *MCA_TabResult;  // MCA component maps
    CAlmR  MCA_ALM;  
    C_UWT  MCA_WT;  
    int MCA_WT_NbrScale;
    int MCA_Alm_Lmax;
    int NbrMCA;
    Hdmap RMSMap;  // Root mean square map, in case of a denoising.
    double MinSigmaNoise;
    Bool UseRMSMap; 
    PowSpec PowSpecData;
    PowSpec PowSpecSignal;
    PowSpec PowSpecRealization;
    fltarray Beam;
    Bool DeconvBeam;
    Bool LogNormal;
    int NbrMasterIter;
    Bool MCA_InitCleabCMB;
    
    Bool WienerFilter; // If denoising, use the Wiener or the Cole Filter
    CInpainting() {Max_Iter=DEF_IWIENER_NBR_ITER; ThresholdResi = True;HardThreshold = True; SoftThreshold = False; Gauss_Cst_KSigma=5.;InpaintMethod = DEF_INP_METHOD;
                  Use_ZeroMeanConstraint = True; Use_WP_Constraint = False; All_WP_Band = False; MCA=False; NbrWP_Band=0;Nside=0;Lmax_Tot=0;L2Constraint=False; NbrMCA=0;
                  IsotropyCst = False; ProjectConstraint = False;ZeroPadding=0.0; Lmax = 0.; Eps=1.; GaussianityConstraint = False;NbrMasterIter=DEF_NBR_MASTER_ITER;
                  GaussCst_ProjectCst_Lambda = 5.; Acceleration=True; Pos=False;Verbose=False; SigmaNoise=0;InputPowSpec=False;Use_Realization_Powspec= False;
                  AccelationLevel=1; LogNormal=False; Acceleration=false;MCA_InitCleabCMB=False; Denoising=False;UseRMSMap=False;MinSigmaNoise=0;WienerFilter=False;DeconvBeam=False;}
    void analyse_alm_inpainting();
    void synthesis_alm_inpainting();
    void synthesis_alm_inpainting_with_master_decconv();

    void old_analyse_alm_inpainting();
    ~CInpainting(){};
    };


#endif
