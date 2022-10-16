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
***************************************************************************** */
  
#ifndef _POLINP_H_
#define _POLINP_H_

#include "MRS_Inp.h"
#include "PolaHealpixClass.h"

#define NBR_POLA_INP_METHOD 1
enum type_pol_inpainting {POL_INP_L1_ALM_ANALYSIS_CST_VAR, POL_INP_L1_ALM_SYNTHESIS, POL_INP_L1_ALM_SYNT_MASTER, POL_INP_ALM_ISOTROP};
#define DEF_POL_INP_METHOD  POL_INP_L1_ALM_ANALYSIS_CST_VAR

const char * StringPolInpaintingMethod (type_pol_inpainting TypeMethod);


class CPInpainting : public CInpainting  
{
 //   CAlmR ALM_Rec;
 //   Hdmap Resi;
     
//    void init_inpainting();
//    void exit_inpainting();
//    double get_residual();
//    void alm_denoising(CAlmR & A, Bool Wiener=False);
    // Alm denoising using the Wiener or the Cole filter.
    //          Den = H2 * (PS_Signal  + PS_Noise )
    //   WienerFilter = PS_Signal  / Den
    //   ColeFilter = sqrt(PS_Signal  / Den)
    double SigDataT, SigDataQ, SigDataU;
    double MeanDataT , MeanDataQ, MeanDataU;
    double MinDataT, MinDataQ, MinDataU;
    int NbrTotMaskT, NbrTotMaskQ, NbrTotMaskU; 

    public:
    PolaHdmap PolMap;     // Map to inpainted
    PolaAlmR  PolALM;  
    PolaHdmap PolMaskMap; // Mask (MaskMap[k] = 0 if not data, and 1 otherwise)
    PolaHdmap PolResult;  // inpainted map
    PolaHdmap PolResi;
    PolaHdmap PolRMSMap;  // Root mean square map, in case of a denoising.
 
    Healpix_PowSpec PolPowSpecData;
    Healpix_PowSpec PolPowSpecSignal;

      CPInpainting(): CInpainting() {}
 //   void analyse_alm_inpainting();
 //   void synthesis_alm_inpainting();
 //   void synthesis_alm_inpainting_with_master_decconv();

     ~CPInpainting(){};
    };




#endif
