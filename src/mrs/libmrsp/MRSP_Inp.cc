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
**    File:  MRS_Inp.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Code for inpainting
**    -------------------
**
******************************************************************************/


#include"MRS_Sparse.h"
#include"MRSP_Inp.h"
#include"MatMask.h"

/*********************************************************************/

const char * StringPolInpaintingMethod (type_pol_inpainting TypeMethod)
{
    switch (TypeMethod)
    {
        case POL_INP_L1_ALM_ANALYSIS_CST_VAR: 
			return ("ALM l1 minimization + wavelet packet variance regularization.");break;
        case POL_INP_L1_ALM_SYNTHESIS: 
			return ("ALM l1 minimization without constraint.");break;
        case POL_INP_L1_ALM_SYNT_MASTER: 
			return ("ALM l1 minimization + Power spectrum deconvolution");break;
        case POL_INP_ALM_ISOTROP: 
			return ("Isotropic Constraint");break;
		default:
			return ("Undefined decomposition");
			break;
    }
}

/*********************************************************************/

#ifdef NOCPP

void CInpainting::pol_init_inpainting()
{
    Min=Max=0;
    NbrTotMask = 0; 
    SigData = SigDataQ = SigDataU = 0.;
    MeanData = MeanDataQ = MeanDataU = 0.;
    FSky = 1.;
    
    if (Verbose == True) PolMap.info((char*) "Input MAP");
    
    // Get the minimum, the stdev and the mean of the data, but in the mask
    MinData=  PolMap.map_T[0];
    MinDataQ =  PolMap.map_Q[0];
    MinDataU =  PolMap.map_U[0];
    for (int p=0; p<  Map.Npix(); p++)  
    {
        if (PolMaskMap.map_T[p] > 0)  
        {
            NbrTotMaskT ++;
            MeanDataT += PolMap.map_T[p];
            if (MinDataT > Map[p]) MinData = PolMap.map_T[p];
            if ((Denoising == True) && (UseRMSMap == True) && (PolRMSMap.map_T[p] == 0))
            {
               cout << "Error: RMS map is equal to zero where the mask is equal to 1 at pixel position in T : " << p << endl;
               exit(-1);
            }
        }
        if (PolMaskMap.map_Q[p] > 0)  
        {
            NbrTotMaskQ ++;
            MeanDataQ += PolMap.map_Q[p];
             if (MinDataQ > PolMap.map_Q[p]) MinDataQ = PolMap.map_Q[p];
 
            if ((Denoising == True) && (UseRMSMap == True) && (PolRMSMap.map_Q[p] == 0))
            {
               cout << "Error: RMS map is equal to zero where the mask is equal to 1 at pixel position in Q: " << p << endl;
               exit(-1);
            }
        }
  
        if (PolMaskMap.map_U[p] > 0)  
        {
            NbrTotMaskU ++;
            MeanDataU += PolMap.map_U[p];
            if (MinDataU > PolMap.map_U[p]) MinDataU = PolMap.map_U[p];
            if ((Denoising == True) && (UseRMSMap == True) && (PolRMSMap.map_U[p] == 0))
            {
               cout << "Error: RMS map is equal to zero where the mask is equal to 1 at pixel position in U : " << p << endl;
               exit(-1);
            }
        }
    }
    MeanDataT /= (double) NbrTotMask;
    MeanDataQ /= (double) NbrTotMask;
    MeanDataU /= (double) NbrTotMask;

    for (int p=0; p<  PolMap.map_T.Npix(); p++)  
         if (PolMaskMap.map_T[p] > 0)  SigDataT += (PolMap.map_T[p]-MeanDataT)*(PolMap.map_T[p]-MeanDataT);
        
    for (int p=0; p<  PolMap.map_Q.Npix(); p++)  
          if (PolMaskMap.map_Q[p] > 0)   SigDataQ += (PolMap.map_Q[p]-MeanDataQ)*(PolMap.map_Q[p]-MeanDataQ);
          
    for (int p=0; p<  PolMap.map_U.Npix(); p++)  
         if (PolMaskMap.map_U[p] > 0)    SigDataU += (PolMap.map_U[p]-MeanDataU)*(PolMap.map_U[p]-MeanDataU);
    }
    SigDataT = sqrt(SigDataT / (double) NbrTotMaskT);
    SigDataQ = sqrt(SigDataQ / (double) NbrTotMaskQ);
    SigDataU = sqrt(SigDataU / (double) NbrTotMaskU);

    FSky = (double) (3*PolMap.map_T.Npix())  / (NbrTotMaskT+ NbrTotMaskQ+ NbrTotMaskU);
    if (Verbose == True) cout << "Percentage of good pixels: " << 1. / FSky * 100. << endl;
    
    // We start with a zero mean image
    if (Use_ZeroMeanConstraint == True)
    {
        for (int p=0; p<  PolMaskMap.map_T.Npix(); p++) if (PolMaskMap.map_T[p] > 0) PolMap.map_T[p] -= MeanDataT;
        for (int p=0; p<  PolMaskMap.map_Q.Npix(); p++) if (PolMaskMap.map_Q[p] > 0) PolMap.map_Q[p] -= MeanDataQ;
        for (int p=0; p<  PolMaskMap.map_U.Npix(); p++) if (PolMaskMap.map_U[p] > 0) PolMap.map_U[p] -= MeanDataU;
    }    
  /*  if (Pos == True) 
    {
        for (int p=0; p<  Map.Npix(); p++) 
            if (MaskMap[p] > 0) Map[p] -= MinData;  
    }  */   
    Nside = PolMap.map_T.Nside();
    Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT = " <<  Lmax_Tot << endl;
    
    PolALM.alloc(Nside, Lmax_Tot, True);
    // ALM.Norm = (IsotropyCst == True) ? True: False;
    // PolALM.Norm = True;
    // ALM.Norm = False;
    // PolALM.set_beam_eff(Lmax, Lmax_Tot);
    // ALM.UseBeamEff = False;
    PolALM_Rec.alloc(Nside, Lmax_Tot, True);
    // PolALM_Rec.Norm = PolALMALM.Norm;
    // PolALMALM_Rec.set_beam_eff(Lmax, Lmax_Tot);
    // PolALMALM_Rec.UseBeamEff = ALM.UseBeamEff;   
    
    PolResult.alloc((int) Nside,   false,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    PolResi.alloc((int) Nside,   false,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
     
    if (SigmaNoise == 0) SigmaNoise = 1.;
    // else if (ALM.Norm == false) SigmaNoise /= ALM.NormVal;
    
    // Find the minimum value in RMSMap in the valid part of the data
    
    if ((Denoising == True) && (UseRMSMap == True))
    {
      double MinSigmaNoise= 0.;
       SigmaNoiseT =  SigmaNoiseQ  = SigmaNoiseU = 0.;
       MinSigmaNoiseT = MinSigmaNoiseQ = MinSigmaNoiseU = 0.;
       for (int p=0; p<  PolRMSMap.map_T.Npix(); p++) 
       { 
          if ((PolMaskMap.map_T[p] > 0) && ((MinSigmaNoiseT > RMSMap[p]) || (MinSigmaNoiseT == 0))) MinSigmaNoise = PolRMSMap.map_T[p];
       }
       SigmaNoiseT = MinSigmaNoise;
       
    }
    else  MinSigmaNoise = SigmaNoise;
    
    if (InputPowSpec == false) mrs_alloc_powspec(PowSpecSignal, Lmax_Tot);
    else 
    {
       double CoefN = ALM.NormVal*ALM.NormVal;
       //cout << "RENORM " << CoefN << endl << endl;
       if (ALM.Norm == True) for (int l=0; l <= ALM.Lmax(); l++) PowSpecSignal.tt(l) *= CoefN;
    }

    mrs_alloc_powspec(PowSpecNoise, Lmax_Tot);
    mrs_alloc_powspec(PowSpecData, Lmax_Tot);
    
#ifdef MRSDEBUG                      
    cout << "ALM.Norm = " << ALM.Norm << endl;
    cout << "ALM.NormVal = " << ALM.NormVal << endl;
    cout << "MinSigmaNoise = " << MinSigmaNoise << endl;
    cout << "Lmax_Tot  = " << Lmax_Tot << endl;
    cout << "MeanData  = " << MeanData << endl;
    cout << "MinData  = " << MinData << endl;
    if (WienerFilter == True) cout << "Wiener filter = " <<  endl;
    else cout << "Cole filter = " <<  endl;
#endif

}


/*********************************************************************/

void CInpainting::exit_inpainting()
{
   double ValAdd = 0.;
   if (Pos == True) ValAdd += MinData;
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
    if (Denoising == False) 
         for (int p=0; p<  Map.Npix(); p++) if (MaskMap[p] > 0) Result[p] = Map[p];

    if (ALM.Norm == True)
    {
        double CoefN = ALM.NormVal;
        for (int l=0; l <= ALM.Lmax(); l++)
	for (int m=0; m <= l; m++)
	{
	    ALM(l,m).re /= CoefN;
            ALM(l,m).im /= CoefN;
	}
    }
}


/*********************************************************************/

void CInpainting::init_wavelet_packet()
{
    ALM_Band.alloc(Nside, Lmax_Tot);
    ALM_Band.Norm = ALM.Norm;
    ALM_Band.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Band.UseBeamEff =  ALM.UseBeamEff ;
    
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
    
    // fits_write_fltarr((char *)"xxw.fits", WP_W);
    //   fits_write_fltarr((char *)"xxh.fits", WP_H);
    // exit(0);
    
    if ((InputPowSpec == True) && (PowSpecSignal.Lmax() > 0))
    {
        if (Verbose == True) cout << " Powspec lmax = " << PowSpecSignal.Lmax() << endl;
        // for (l=0; l <= PowSpecData.Lmax(); l++)  PowSpecData.tt(l) *= ALM.NormVal;
        if (Use_WP_Constraint == True)
        {
	    double CoefN = ALM.NormVal*ALM.NormVal;
	    if (ALM.Norm == False) CoefN = 1.;
            double EnerBand;
            TabPsi2.alloc(NbrWP_Band+1);
            for (int b=0; b <= NbrWP_Band; b++)
            {
                EnerBand = 0.;
                for (int l=0; l <= MIN(PowSpecData.Lmax(), WP_H.nx()-1); l++) 
		{
		   double P = (PowSpecSignal.tt)(l) / CoefN; 
                    EnerBand += (2.*l+1)* P * WP_WPFilter(l, b)*WP_WPFilter(l, b);
                }
		TabPsi2(b)  = EnerBand  /  (4. * PI);
		// cout << b+1 << " " << TabPsi2(b) << " " << (PowSpecSignal.tt)(10) << endl;
            }
        }    
    }
}
/*********************************************************************/


void CInpainting::wp_band_constraint(Hdmap &Band, int b)
{
    double EnergyInMap = 0.;
    for (int p=0; p < Map.Npix(); p++)  if (MaskMap[p] != 0) EnergyInMap += Band[p]*Band[p];
    double Sig_InMap  = sqrt(EnergyInMap  /  (double) (NbrTotMask));
    
    
    // Compute the expected power in the band
    double EnergyBand=0.;
    double Np = (double) Map.Npix();
    // if (InputPowSpec == True)  EnergyBand = Np * TabPsi2(b);
    // else 
    EnergyBand = EnergyInMap * Np / (double) NbrTotMask;
    double SigBand = (EnergyBand > 0) ? sqrt(EnergyBand / Np): 0;
    
    // double EnergyBand = Map.Npix() * SigBand*SigBand;
    
    // projection. All  coeff > Gauss_Cst_KSigma*SigBand are set to zero
    if (ProjectConstraint ==  True)
        for (int p=0; p < Map.Npix(); p++) 
            if ((MaskMap[p] == 0) &&  (ABS(Band[p]) > Gauss_Cst_KSigma * SigBand))  Band[p] = 0.;
    
    double EnergyInMask = 0.;
    for (int p=0; p < Map.Npix(); p++) 
    {
        // float Level = GaussCst_ProjectCst_Lambda * Sig_InMap;
        if (MaskMap[p] == 0) 
        {
            if  ((GaussianityConstraint ==  True) && (b != NbrWP_Band))
            {
                // if (ABS(Band[p]) > Level) Band[p] = 0.;
                // float SoftT = ABS(Band[p]) - GaussCst_Lambda *  Sig_InMap;
                float SoftT = ABS(Band[p]) - GaussCst_ProjectCst_Lambda *  Sig_InMap;
                if (SoftT < 0) SoftT = 0.;
                Band[p] = soft_threshold(Band[p], SoftT);
            }
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
    
    for (int p=0; p < Map.Npix(); p++)  if (MaskMap[p] == 0)  Band[p] *= Coeff;
}

/*********************************************************************/
 

void CInpainting::alm_denoising(CAlmR & A, Bool Wiener)
{

    double Coef;
    double H2=1.; // in case of PSF, H is the square of the beam
 
    A.alm2powspec(PowSpecData);
    if (InputPowSpec == false)
   {
      for (int l=0; l <= A.Lmax(); l++) 
      {
         PowSpecNoise.tt(l) = MinSigmaNoise*MinSigmaNoise;
         PowSpecSignal.tt(l) = MAX(0, PowSpecData.tt(l) - PowSpecNoise.tt(l));
      }
   }
   else
   {
       float NormVal = A.NormVal*A.NormVal;
	
       // cout << "DENOISING:  " <<  PowSpecData.tt(10) / NormVal << "  " <<  PowSpecSignal.tt(10) / NormVal<< endl;
       for (int l=0; l <= A.Lmax(); l++) PowSpecNoise.tt(l) = MAX(0, PowSpecData.tt(l) - PowSpecSignal.tt(l));
       // for (int l=0; l <= A.Lmax(); l++) PowSpecNoise.tt(l) = 0.;
   }
   
    for (int l=0; l <= A.Lmax(); l++)
    {
        double Den = H2*PowSpecData.tt(l);  // (PowSpecSignal.tt(l) + PowSpecNoise.tt(l));
        if (Wiener == True) Coef = (Den > 0) ? PowSpecSignal.tt(l) / Den: 0.; 
        else Coef = (Den > 0) ?  sqrt( PowSpecSignal.tt(l) / Den): 0.;    // Cole Filter
 	if (l == 0) Coef=1.;
	for (int m=0; m <= l; m++)
        {
            A(l,m).re = A(l,m).re * Coef;
            A(l,m).im = A(l,m).im * Coef;
        }
    }
}


/*********************************************************************/

// Compute the residual
double CInpainting::get_residual()
{
    double Err = 0.;
    
    for (int p=0; p < Result.Npix(); p++) 
    {
        Resi[p] =  Map[p] - Result[p];
        Resi[p] *= MaskMap[p];
        if ((Denoising == True) && (MaskMap[p] > 0) &&  (UseRMSMap == True) && (RMSMap[p] != 0))
        {
            Resi[p] *= MinSigmaNoise / RMSMap[p];
        }
        Err += Resi[p]*Resi[p];
    }
    Err = sqrt(Err / (double) (NbrTotMask)) / MinSigmaNoise;
    return Err;
}

/*********************************************************************/


void CInpainting::analyse_alm_inpainting()
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

    if (Debug == True) Map.info("Map ");
    ALM.alm_trans(Map);
     
    // ALM.write((char *) "xx_init_alm.fits", True);
    
    double MaxAlm =ALM.max_absalm();
    // IDL: NormVal =  sqrt(nel / (4.*!DPI))
    MaxAlm *= 0.9999999;
    if (ALM.Norm == False) SigmaNoise /= ALM.NormVal;  // Normally ALM.Norm should be set to True,
    
   //  if (Debug == True)
    {
       if (ALM.Norm == True) cout << " Alm norm: max =  " << MaxAlm <<  endl;
       else cout << " Alm nonorm: max =  " << MaxAlm << endl;
       printf(" MaxAlm = %5.4f, SigNoise = %f\n", (float) MaxAlm / SigmaNoise, SigmaNoise);
    }
    
    double FirstSoftThreshold = MaxAlm / SigmaNoise * sqrt(2.);  // sqrt(2.) is here, to be equivalent to the IDL code.
                                                                 // this has no effect since we redivide latter by the same value.
    double LastSoftThreshold = 0.;
    double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    // double StepL = DeltaThreshold / (float) (Max_Iter -1);
    double Lambda = FirstSoftThreshold;
    if (Debug == True)  cout << "==> Max Alm = " << MaxAlm << " " << (MaxAlm * SigmaNoise)  << ", N = " << ALM.NormVal <<  ", First Threshold = " << Lambda << endl;

    double Err=Map.rms();
    int MaxNonZeroL=0; 
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
        sprintf(NN, "xx_input.fits");
        Map.write(NN);
    }
    Lambda = FirstSoftThreshold;
    for (int i=0; i < Max_Iter; i++)
    { 
       Err = get_residual();
       if ((ThresholdResi == False) || (IsotropyCst == True)) Resi += Result;

       if (Debug == True)
       {
            char NN[256];
            sprintf(NN, "xx_resi_it_%d.fits", i+1);
            if (i >= 0) Resi.write(NN);
       }
        int CptAlm=0;
        
        if (Acceleration == True) ALM.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM.alm_trans(Resi);
        // MaxAlm =ALM.max_absalm() / SigmaNoise;
        // cout << " Max ABS(ALM) Resi  " << MaxAlm  << endl;
        if (i > 0) Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * ((float) i / (float) Max_Iter)));  // 4 versus 2.8  in IDL code
        if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;	      
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        if (Debug == True) printf(" Lambda = %5.4f\n", (float) Lambda);
        if ((Verbose == True)  && (IsotropyCst == False)) cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
        if ((Verbose == True)  && (IsotropyCst == True)) cout << "Iter  " << i+1 << ", Lambda = " << GaussCst_Lambda  <<  " Err = " << Err << endl;
        if (Debug == True) Resi.info("RESI ");
        
        if (Denoising == True) alm_denoising(ALM, WienerFilter);

        if ((HardThreshold == True)  || (SoftThreshold == True))  // Alm L1 minimization
        { 
            float ThresholdLevel = SigmaNoise*Lambda / sqrt(2.);
            if (Lambda > 0)
            {
                if (HardThreshold == True) 
                {
                 // cout << "Max AlmResi = " << ALM.max_absalm() << endl;
                // ALM.write("xxalm.fits", true);
                   CptAlm = ALM.hard_threshold(ThresholdLevel, MaxNonZeroL);
                }
                else if (SoftThreshold == True)   CptAlm = ALM.soft_threshold(ThresholdLevel/sqrt(2.), MaxNonZeroL);
                if (Debug == True) cout << "  Number of significant Alm : = " << CptAlm <<  " MaxL = " << MaxNonZeroL << endl; 
                
                

                
             /*   if (MaxNonZeroL >= 1) 
	            {
                    if (MaxNonZeroL >= ALM_FIRST_L) ALM_Rec.alm_rec(Resi);
                    else
                    {
                        ALM_Rec.Set(MaxNonZeroL,  MaxNonZeroL);
                        for (int l=0; l <= MaxNonZeroL; l++)
                            for (int m=0; m <= l; m++) ALM_Rec(l,m) = ALM(l,m);
                        ALM_Rec.alm_rec(Resi);
                    }
                }
	            else Resi.init();
             */   
            }
            else  MaxNonZeroL = Lmax_Tot;
             // ALM.info((char *) "ALM INFO  1");
  

            // Update the solution
            if (ThresholdResi == True) 
            {
                ALM.alm_rec(Resi);
                for (int p=0; p<  Result.Npix(); p++)   Result[p]  +=  Eps*Resi[p];
                ALM.alm_trans(Result);
            }   
	    else ALM.alm_rec(Result, True);
             
            // To accelerate, we compute the Alm only to MaxNonZeroL*2
            // if (Acceleration == True)  
            {
                if (2*MaxNonZeroL >  ALM_FIRST_L)  ALM_FIRST_L = 2*MaxNonZeroL;
                ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
            }
        } // end l1 minimization i.e.  ((HardThreshold == True)  || (SoftThreshold == True))
	    
        if ((IsotropyCst == True) && (InputPowSpec == True))
        {
            isotropy_constraint(ALM, PowSpecSignal, GaussCst_Lambda, False);
        }
	
//             ALM.info((char *) "ALM INFO  2");
        
//                ALM.info((char *) "ALM INFO 3");

        // Apply constraint per wavelet packet band. We start these constraint only at half the total number of iterations.
        // Use_WP_Constraint = False;
        if ((Use_WP_Constraint == True) && (i > 0) && (i % AccelationLevel == 0))
        {
            if (Debug == True) cout << "        Wavelet Packets constraint ... " << endl;
            // ALM.info((char *) "ALM INFO 4");

            // We compute the WP transform:
            
            //   ITER: Result contains the high resolution   (c_j)
            //   We compute Band, which contains the low resolution (c_{j+1}
            //   Coef =  Result - Band
            //   Result = Band
            
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
                             ALM_Band(l,m).re = ALM(l,m).re * WP_H(l, NbrWP_Band-1-b);
                             ALM_Band(l,m).im = ALM(l,m).im * WP_H(l, NbrWP_Band-1-b);
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
                        for (int p=0; p < Map.Npix(); p++) 
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
                    wp_band_constraint(Band,  b);
                                  
                    // coadd the band to new new solution
                    Resi += Band;
                } // else if  (LMin > ALM_FIRST_L) 
            } // endfor b
            Result = Resi;
            
        }   // endif Use_WP_Constraint
	else if (IsotropyCst == True)  ALM.alm_rec(Result);
        
        
        // We want a zero mean solution
        if (Use_ZeroMeanConstraint == True)
        {
            if (Debug == True)  cout << "ZERO" << endl;
        
            double MeanRes = Result.average();
            for (int p=0; p < Map.Npix(); p++)  Result[p] -= MeanRes;
        }
       
         if (Pos == True)
        {
            if (Debug == True)  cout << "POS" << endl;
         for (int p=0; p < Map.Npix(); p++) if (Result[p] < 0) Result[p] = 0.;
        }
 /*     {
              if (Debug == True)  cout << "POS" << endl;

            double ExpMean = MeanData - MinData;
            double MeanMask=0.;
            for (int p=0; p < Map.Npix(); p++)  
            {
                if (Result[p] < 0) Result[p] = 0.;
                if (MaskMap[p] == 0) MeanMask += Result[p];
            }
            MeanMask /= (double)  (Map.Npix() -  NbrTotMask);
            for (int p=0; p < Map.Npix(); p++)  
                if (MaskMap[p] == 0)   Result[p] += ExpMean - MeanMask;
        } */
        
        if (Debug == True)
        {
         char NN[256];
         sprintf(NN, "xx_res_it_%d.fits", i+1);
         if (i >= 0) Result.write(NN);
         }
        if (Verbose == True) Result.info((char *) "  REC");
    }
   
    exit_inpainting();
    
}

/*********************************************************************/

void CInpainting::old_analyse_alm_inpainting()
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

    if (Debug == True) Map.info("Map ");
    ALM.alm_trans(Map);
     
    // ALM.write((char *) "xx_init_alm.fits", True);
    
    double MaxAlm =ALM.max_absalm();
    // IDL: NormVal =  sqrt(nel / (4.*!DPI))
    MaxAlm *= 0.9999999;
    if (ALM.Norm == False) SigmaNoise /= ALM.NormVal;  // Normally ALM.Norm should be set to True,
    
    if (Debug == True)
    {
       if (ALM.Norm == True) cout << " Alm norm: max =  " << MaxAlm <<  endl;
       else cout << " Alm nonorm: max =  " << MaxAlm << endl;
       printf(" MaxAlm = %5.4f, SigNoise = %f\n", (float) MaxAlm / SigmaNoise, SigmaNoise);
    }
    
    double FirstSoftThreshold = MaxAlm / SigmaNoise * sqrt(2.);  // sqrt(2.) is here, to be equivalent to the IDL code.
                                                                 // this has no effect since we redivide latter by the same value.
    double LastSoftThreshold = 0.;
    double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    // double StepL = DeltaThreshold / (float) (Max_Iter -1);
    double Lambda = FirstSoftThreshold;
    if (Debug == True)  cout << "==> Max Alm = " << MaxAlm << " " << (MaxAlm * SigmaNoise)  << ", N = " << ALM.NormVal <<  ", First Threshold = " << Lambda << endl;

    double Err=Map.rms();
    int MaxNonZeroL=0; 
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
        sprintf(NN, "xx_input.fits");
        Map.write(NN);
    }
    Lambda = FirstSoftThreshold;
    for (int i=0; i < Max_Iter; i++)
    { 
       Err = get_residual();
       if ((ThresholdResi == False) || (IsotropyCst == True)) Resi += Result;

       if (Debug == True)
       {
            char NN[256];
            sprintf(NN, "xx_resi_it_%d.fits", i+1);
            if (i >= 0) Resi.write(NN);
       }
        int CptAlm=0;
        
        if (Acceleration == True) ALM.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM.alm_trans(Resi);
        // MaxAlm =ALM.max_absalm() / SigmaNoise;
        // cout << " Max ABS(ALM) Resi  " << MaxAlm  << endl;
        if (i > 0) Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * ((float) i / (float) Max_Iter)));  // 4 versus 2.8  in IDL code
        if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;	      
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        if (Debug == True) printf(" Lambda = %5.4f\n", (float) Lambda);
        if ((Verbose == True)  && (IsotropyCst == False)) cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
        if ((Verbose == True)  && (IsotropyCst == True)) cout << "Iter  " << i+1 << ", Lambda = " << GaussCst_Lambda  <<  " Err = " << Err << endl;
        if (Debug == True) Resi.info("RESI ");
        
        if (IsotropyCst == True)
        {
            // ALM.info((char *) "ALM bef");
            for (int l=0; l <=  ALM.Lmax(); l++)
            for (int m=0; m <= l; m++) ALM_Rec(l,m) = ALM(l,m);
            Bool UpdatePowSpec =(InputPowSpec == true) ? False: True;
            isotropy_constraint(ALM, PowSpecData, GaussCst_Lambda, UpdatePowSpec);
            
            // for (int l=0; l <=  ALM.Lmax(); l++)
            // for (int m=0; m <= l; m++) ALM_Rec(l,m) -= ALM(l,m);
            // ALM_Rec.info((char *) "ALM Diff");
        }
        else   // Alm L1 minimization
        { 
            float ThresholdLevel = SigmaNoise*Lambda / sqrt(2.);
            if (Lambda > 0)
            {
                if (HardThreshold == True) 
                {
                 // cout << "Max AlmResi = " << ALM.max_absalm() << endl;
                // ALM.write("xxalm.fits", true);
                   CptAlm = ALM.hard_threshold(ThresholdLevel, MaxNonZeroL);
                }
                else if (SoftThreshold == True)   CptAlm = ALM.soft_threshold(ThresholdLevel/sqrt(2.), MaxNonZeroL);
                if (Debug == True) cout << "  Number of significant Alm : = " << CptAlm <<  " MaxL = " << MaxNonZeroL << endl; 
                
                

                
             /*   if (MaxNonZeroL >= 1) 
	            {
                    if (MaxNonZeroL >= ALM_FIRST_L) ALM_Rec.alm_rec(Resi);
                    else
                    {
                        ALM_Rec.Set(MaxNonZeroL,  MaxNonZeroL);
                        for (int l=0; l <= MaxNonZeroL; l++)
                            for (int m=0; m <= l; m++) ALM_Rec(l,m) = ALM(l,m);
                        ALM_Rec.alm_rec(Resi);
                    }
                }
	            else Resi.init();
             */   
            }
            else  MaxNonZeroL = Lmax_Tot;
            //  ALM.info((char *) "ALM INFO  1");
           
            // if (Denoising == True) cout << "DENOISING ! " << endl;
            if (Denoising == True) alm_denoising(ALM, WienerFilter);
  
            // MaxNonZeroL = Lmax_Tot;
            ALM.alm_rec(Resi, True);
            // ALM.info((char *) "ALM INFO 2");
            
            // Update the solution
            // if (MaxNonZeroL > 1)
            {
                if (ThresholdResi == False) 
                {
                    for (int p=0; p<  Result.Npix(); p++)   Result[p] = (1-Eps)*Result[p] + Eps*Resi[p];
                }
                else for (int p=0; p<  Result.Npix(); p++)   Result[p]  +=  Eps*Resi[p];
            }   
            
            
            // To accelerate, we compute the Alm only to MaxNonZeroL*2
            // if (Acceleration == True)  
            {
                if (2*MaxNonZeroL >  ALM_FIRST_L)  ALM_FIRST_L = 2*MaxNonZeroL;
                ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
            }
        } // end l1 minimization
         //     ALM.info((char *) "ALM INFO 3");
	    
	    
        // Apply constraint per wavelet packet band. We start these constraint only at half the total number of iterations.
        // Use_WP_Constraint = False;
        if ((Use_WP_Constraint == True) && (i > 0) && (i % AccelationLevel == 0))
        {
            if (Debug == True) cout << "        Wavelet Packets constraint ... " << endl;
            // ALM.info((char *) "ALM INFO");
            // We compute the WP transform:
            
            //   ITER: Result contains the high resolution   (c_j)
            //   We compute Band, which contains the low resolution (c_{j+1}
            //   Coef =  Result - Band
            //   Result = Band
            
            if (IsotropyCst == True)  
            {
                ALM.alm_rec(Result);
                for (int p=0; p<  Map.Npix(); p++) if (MaskMap[p] > 0) Result[p] = Map[p];
                ALM.alm_trans(Result);
            }
            Resi.init();
            for (int b=0; b <= NbrWP_Band; b++)
            {  
                int LMin = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,0): WP_W(0,0);
                int LMax = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,1): WP_W(0,1);
                
                // if (Verbose == True) cout << "        WP:Band " << b+1 << ", Lmin = " << LMin << ", Lmax = " << LMax << endl; 
                
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
                             ALM_Band(l,m).re = ALM(l,m).re * WP_H(l, NbrWP_Band-1-b);
                             ALM_Band(l,m).im = ALM(l,m).im * WP_H(l, NbrWP_Band-1-b);
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
                        for (int p=0; p < Map.Npix(); p++) 
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
                    wp_band_constraint(Band,  b);
                                  
                    // coadd the band to new new solution
                    Resi += Band;
                } // else if  (LMin > ALM_FIRST_L) 
            } // endfor b
            Result = Resi;
            
        }   // endif Use_WP_Constraint
	    else if (IsotropyCst == True)  ALM.alm_rec(Result);
        
        
        // We want a zero mean solution
        if (Use_ZeroMeanConstraint == True)
        {
            if (Debug == True)  cout << "ZERO" << endl;
        
            double MeanRes = Result.average();
            for (int p=0; p < Map.Npix(); p++)  Result[p] -= MeanRes;
        }
        
        if (Pos == True)
        {
              if (Debug == True)  cout << "POS" << endl;

            double ExpMean = MeanData - MinData;
            double MeanMask=0.;
            for (int p=0; p < Map.Npix(); p++)  
            {
                if (Result[p] < 0) Result[p] = 0.;
                if (MaskMap[p] == 0) MeanMask += Result[p];
            }
            MeanMask /= (double)  (Map.Npix() -  NbrTotMask);
            for (int p=0; p < Map.Npix(); p++)  
                if (MaskMap[p] == 0)   Result[p] += ExpMean - MeanMask;
        }
        
        if (Debug == True)
        {
         char NN[256];
         sprintf(NN, "xx_res_it_%d.fits", i+1);
         if (i >= 0) Result.write(NN);
         }
        if (Verbose == True) Result.info((char *) "  REC");
    }
    
    exit_inpainting();
     
    
}


/*********************************************************************/
// ../cxx/mrs_alm_inpainting -v -R -s -i100 -H cmb512mask.fits mask512.fits tt ttp
//  mrs_alm_inpainting -v -H -R -s -i100 -M union_matmask.fits  lensed_map-11187_Lobb_hfi6chan-29532.fits  union_mask.fits  tts1 ttps1
// no threshold resi, synthesis, 100iter
//   ../cxx/mrs_alm_inpainting -v -H -r -s -i100 -M matmask512.fits -R revmatmask512.fits cmb512mask.fits mask512.fits tt ttp
void CInpainting::synthesis_alm_inpainting()
{
    double Sig_InMap, Sig_Mask;
    if (Verbose == True)  cout << " Synthesis ALM inpainting: Niter = " << Max_Iter << endl;
     
    init_inpainting();
    
    if (Verbose == True)  Map.info((char *) "S INIT");
    // ALM.Norm = False;
    ALM.alm_trans(Map);
    
    double MaxAlm =ALM.max_absalm();
    MaxAlm *= 0.9999999;
    if (Verbose == True)  cout << "==> Max Alm = " << MaxAlm << " N = " << ALM.NormVal << endl;
    
    double FirstSoftThreshold = MaxAlm / SigmaNoise;
    if (SoftThreshold == True) FirstSoftThreshold /= sqrt(2.);
    double LastSoftThreshold = 0.;
    double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    // double StepL = DeltaThreshold / (float) (Max_Iter -1);
    double Lambda = FirstSoftThreshold;
    
    double Err=Map.rms();
    int MaxNonZeroL=0; 
    Acceleration = False;
    int ALM_FIRST_L = (Acceleration == True) ? 100: Lmax_Tot;
    if (IsotropyCst == True) ALM_FIRST_L = Lmax_Tot;
    // int Ni = Max_Iter;
    // int Fi = 0;
    
    // Use for isotropic constrain only
    float GaussCst_LambdaFirst = 0.;
    float GaussCst_StepLambda = (GaussCst_ProjectCst_Lambda - GaussCst_LambdaFirst)  / float(Max_Iter);
    float GaussCst_Lambda = GaussCst_LambdaFirst;
    
    
    Result.init();
    
    CAlmR ALM_Sol;
    ALM_Sol.alloc(Nside, Lmax_Tot);
    ALM_Sol.Norm = ALM.Norm;
    ALM_Sol.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Sol.UseBeamEff =  ALM.UseBeamEff ;
    
    long Nalm = (Lmax_Tot+1)*(Lmax_Tot+1) / 2;
    for (int i=0; i < Max_Iter; i++)
    { 
        int CptAlm=0;
        // cout << i+1 << endl;
        ALM_Rec.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM_Rec.Norm = ALM.Norm;
	
	int Test=1;
 	if (Test == 1)
	{
        Err = get_residual();
        ALM_Rec.alm_trans(Resi);
        for (int l=0; l <=  ALM_Rec.Lmax(); l++)
        for (int m=0; m <= l; m++) 
        {
            ALM_Sol(l,m).re += ALM_Rec(l,m).re;
            ALM_Sol(l,m).im += ALM_Rec(l,m).im;
        }
	}
	else
	{
        if (i > 0)
        {
            for (int l=0; l <=  ALM_Rec.Lmax(); l++)
            for (int m=0; m <= l; m++) ALM_Rec(l,m) = ALM_Sol(l,m);
            ALM_Rec.alm_rec(Result);
            // {
            //    char NN[256];
            //    sprintf(NN, "xx_res_it_%d.fits", i+1);
            //    if (i >= 0) Result.write(NN);
            // }
            
            if (Verbose == True) Result.info((char *) "  RES");
            for (int p=0; p < Result.Npix(); p++)  Result[p] =  Result[p]*MaskMap[p];
            ALM_Rec.alm_trans(Result);
        }
        
        for (int l=0; l <=  ALM_Rec.Lmax(); l++)
            for (int m=0; m <= l; m++) 
            {
                ALM_Rec(l,m).re = ALM(l,m).re - ALM_Rec(l,m).re;
                ALM_Rec(l,m).im = ALM(l,m).im - ALM_Rec(l,m).im;
                Err += norm( ALM_Rec(l,m) );
                ALM_Sol(l,m).re += ALM_Rec(l,m).re;
                ALM_Sol(l,m).im += ALM_Rec(l,m).im;
            }
        Err /= (float) Nalm;
        Err = sqrt(Err);
	}
	
	
        
        // cout << FirstSoftThreshold << " " << LastSoftThreshold  << endl;
        Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * ((float) i / (float) Max_Iter)));  // 4 versus 2.8  in IDL code
        if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;	  
        
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        
        if ((Verbose == True)  && (IsotropyCst == False)) cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
        if ((Verbose == True)  && (IsotropyCst == True)) cout << "Iter  " << i+1 << ", Lambda = " << GaussCst_Lambda  <<  " Err = " << Err << endl;
        
        if (Denoising == True) alm_denoising(ALM_Sol, WienerFilter);

        float ThresholdLevel = SigmaNoise*Lambda;
        if (Lambda > 0)
        {
            if (HardThreshold == True) CptAlm = ALM_Sol.hard_threshold(ThresholdLevel, MaxNonZeroL);
            else if (SoftThreshold == True)   CptAlm = ALM_Sol.soft_threshold(ThresholdLevel*sqrt(2.), MaxNonZeroL);
            // if (Verbose == True) cout << "  Number of significant Alm : = " << CptAlm <<  " MaxL = " << MaxNonZeroL << endl; 
        }
        else  MaxNonZeroL = Lmax_Tot;
        
        // To accelerate, we compute the Alm only to MaxNonZeroL*2
        if (Acceleration == True)  
        {
            if (2*MaxNonZeroL >  ALM_FIRST_L)  ALM_FIRST_L = 2*MaxNonZeroL;
            ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
        }
     
       for (int l=0; l <=  ALM.Lmax(); l++)
       for (int m=0; m <= l; m++) ALM(l,m) = ALM_Sol(l,m);
       ALM.alm_rec(Result, True);
    }
    
    // Apply constraint per wavelet packet band. We start these constraint only at half the total number of iterations.
    // Use_WP_Constraint = False;
    if (Use_WP_Constraint == True)  
    {
        init_wavelet_packet();
        if (Verbose == True) cout << "        Wavelet Packets constraint ... " << endl;
        
        Resi.init();
        for (int b=0; b <= NbrWP_Band; b++)
        {  
            int LMin = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,0): WP_W(0,0);
            int LMax = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,1): WP_W(0,1);
            
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
                    ALM_Band(l,m).re = ALM(l,m).re * WP_H(l, NbrWP_Band-1-b);
                    ALM_Band(l,m).im = ALM(l,m).im * WP_H(l, NbrWP_Band-1-b);
                }
                ALM_Band.alm_rec(Band);                                    
                for (int p=0; p < Map.Npix(); p++) 
                {
                        double Val = Band[p];
                        Band[p] = Result[p] - Val;
                        Result[p] = Val;
                }
            } // endelse if (b == NbrWP_Band)
            // Wavelet packet constraint on the band b
            wp_band_constraint(Band,  b);
                
            // coadd the band to new new solution
            Resi += Band;
        } // endfor b
        Result = Resi;
    }   // endif Use_WP_Constraint

    
    if (Pos == True)
    {
         for (int p=0; p < Result.Npix(); p++) if (Result[p] < 0) Result[p] = 0.;
    }
    
    exit_inpainting();
}


/*********************************************************************/

/*********************************************************************/

/*********************************************************************/
// ../cxx/mrs_alm_inpainting -v -R -s -i100 -H cmb512mask.fits mask512.fits tt ttp
//  mrs_alm_inpainting -v -H -R -s -i100 -M union_matmask.fits  lensed_map-11187_Lobb_hfi6chan-29532.fits  union_mask.fits  tts1 ttps1
// no threshold resi, synthesis, 100iter
//   ../cxx/mrs_alm_inpainting -v -H -r -s -i100 -M matmask512.fits -R revmatmask512.fits cmb512mask.fits mask512.fits tt ttp
void CInpainting::synthesis_alm_inpainting_with_master_decconv()
{
    double Sig_InMap, Sig_Mask;
    if (Verbose == True)  cout << " Synthesis ALM inpainting: Niter = " << Max_Iter << endl;
    Verbose = True;
    
    init_inpainting();
    if (Use_WP_Constraint == True) init_wavelet_packet();
    
    if (Verbose == True)  Map.info((char *) "S INIT");
    ALM.Norm = False;
     
    ALM.alm_trans(Map);
    
    // Compute the matrix related to the mask
    if (MatMask.naxis() == 2) 
    {
       if ((MatMask.nx() != Lmax_Tot+1) || (MatMask.ny() != Lmax_Tot+1))
       {
          cout << "Error: the input mask matrix has not the correct dimension: (" << MatMask.nx() << "," << MatMask.ny() << "), instead of " <<  Lmax_Tot << endl;
          exit(-1);
       }
    }
    else 
    {
        if (Verbose == True)  cout << " Compute the mask matrix: Lmax = "  << Lmax_Tot << endl;
        
        MatMask = mrs_matmask( MaskMap, Lmax_Tot);
        if (Verbose == True)  cout << " Computed Matrix: Nx = " << MatMask.nx() << ", Ny = " << MatMask.ny() << endl;
        fits_write_dblarr("xx_matmask.fits", MatMask );
    }
    if (RevMatMask.naxis() == 2) 
    {
        if ((RevMatMask.nx() != Lmax_Tot+1) || (RevMatMask.ny() != Lmax_Tot+1))
        {
            cout << "Error: the input reverse mask matrix has not the correct dimension: (" << RevMatMask.nx() << "," << RevMatMask.ny() << "), instead of " <<  Lmax_Tot << endl;
            exit(-1);
        }
    }
    else 
    {
        if (Verbose == True)  cout << " Compute the mask matrix: Lmax = "  << Lmax_Tot << endl;
        for (int p=0; p < MaskMap.Npix(); p++)  Result[p] = 1. - MaskMap[p];
        RevMatMask = mrs_matmask( Result, Lmax_Tot);
        if (Verbose == True)  cout << " Computed Reverse Mask Matrix: Nx = " << RevMatMask.nx() << ", Ny = " << RevMatMask.ny() << endl;
        fits_write_dblarr("xx_revmatmask.fits", RevMatMask );
    }
    
    PowSpec ClData, ClResi, ClResi1, ClSol, ClUpdateSol, ClConv, ClM0_Exp, ClM0, ClM1, ClM0M1;
    mrs_alloc_powspec(ClSol, Lmax_Tot);
    mrs_alloc_powspec(ClConv, Lmax_Tot);
    mrs_alloc_powspec(ClResi, Lmax_Tot);
    mrs_alloc_powspec(ClData, Lmax_Tot);
    mrs_alloc_powspec(ClResi1, Lmax_Tot);
    mrs_alloc_powspec(ClM0_Exp, Lmax_Tot);
    mrs_alloc_powspec(ClM0, Lmax_Tot);
    mrs_alloc_powspec(ClM1, Lmax_Tot);
    mrs_alloc_powspec(ClM0M1, Lmax_Tot);
    mrs_alloc_powspec(ClUpdateSol, Lmax_Tot);
    ALM.alm2powspec(ClData);

    // ALM.write((char *) "xx_init_alm.fits", True);
    
    double MaxAlm =ALM.max_absalm();
    MaxAlm -= 0.00001 * MaxAlm;
    if (Verbose == True)  cout << "==> Max Alm = " << MaxAlm << " N = " << ALM.NormVal << endl;
    
    double FirstSoftThreshold = MaxAlm / SigmaNoise;
    if (SoftThreshold == True) FirstSoftThreshold /= sqrt(2.);
    double LastSoftThreshold = 0.;
    double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    // double StepL = DeltaThreshold / (float) (Max_Iter -1);
    double Lambda = FirstSoftThreshold;
    
    double Err=Map.rms();
    int MaxNonZeroL=0; 
    Acceleration = False;
    int ALM_FIRST_L = (Acceleration == True) ? 100: Lmax_Tot;
    if (IsotropyCst == True) ALM_FIRST_L = Lmax_Tot;
    // int Ni = Max_Iter;
    // int Fi = 0;
    
    // Use for isotropic constrain only
    float GaussCst_LambdaFirst = 0.;
    float GaussCst_StepLambda = (GaussCst_ProjectCst_Lambda - GaussCst_LambdaFirst)  / float(Max_Iter);
    float GaussCst_Lambda = GaussCst_LambdaFirst;
    
    // fits_write_fltarr((char *)"xxpsi.fits", TabPsi2);
    // wp_trans(Map,  Band,  ALM,  ALM_Band,  WP_W,  WP_H,  NbrWP_Band,  ALM_FIRST_L, TabPsi2);
    // exit(0);
    
    // MAIN ITER 
 
    Result.init();
    
    CAlmR  ALM_M0, ALM_M1, ALM_Sol;
    ALM_M0.alloc(Nside, Lmax_Tot);
    ALM_M0.Norm = ALM.Norm;
    ALM_M0.set_beam_eff(Lmax, Lmax_Tot);
    ALM_M0.UseBeamEff =  ALM.UseBeamEff ;
    ALM_M1.alloc(Nside, Lmax_Tot);
    ALM_M1.Norm = ALM.Norm;
    ALM_M1.set_beam_eff(Lmax, Lmax_Tot);
    ALM_M1.UseBeamEff =  ALM.UseBeamEff ;
    ALM_Sol.alloc(Nside, Lmax_Tot);
    ALM_Sol.Norm = ALM.Norm;
    ALM_Sol.set_beam_eff(Lmax, Lmax_Tot);
    ALM_Sol.UseBeamEff =  ALM.UseBeamEff ;
    
    long Nalm = (Lmax_Tot+1)*(Lmax_Tot+1) / 2;
    for (int i=0; i < Max_Iter; i++)
    { 
        int CptAlm=0;
    // cout << i+1 << endl;
        Err = 0.;
        ALM_Rec.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM_Rec.Norm = ALM.Norm;
        ALM_M1.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM_M1.Norm = ALM.Norm;
        ALM_M0.Set(ALM_FIRST_L,ALM_FIRST_L);
        ALM_M0.Norm = ALM.Norm;
        
        if (i > 0)
        {
            for (int l=0; l <=  ALM_Rec.Lmax(); l++)
            for (int m=0; m <= l; m++) ALM_Rec(l,m) = ALM_Sol(l,m);
            ALM_Rec.alm_rec(Result);
            //{
            //     char NN[256];
            //     sprintf(NN, "xx_res_it_%d.fits", i+1);
            //     if (i >= 0) Result.write(NN);
           // }
            
            if (Verbose == True) Result.info((char *) "  RES");
            for (int p=0; p < Result.Npix(); p++)  Result[p] =  Result[p]*MaskMap[p];
            ALM_M1.alm_trans(Result);
        }
        
        for (int l=0; l <=  ALM_Rec.Lmax(); l++)
        for (int m=0; m <= l; m++) 
        {
            ALM_Rec(l,m).re = ALM(l,m).re - ALM_M1(l,m).re;
            ALM_Rec(l,m).im = ALM(l,m).im - ALM_M1(l,m).im;
            Err += norm( ALM_Rec(l,m) );
            ALM_M0(l,m).re = ALM_Sol(l,m).re - ALM_M1(l,m).re;
            ALM_M0(l,m).im = ALM_Sol(l,m).im - ALM_M1(l,m).im;
            
            ALM_Sol(l,m).re = ALM(l,m).re + ALM_M0(l,m).re;
            ALM_Sol(l,m).im = ALM(l,m).im + ALM_M0(l,m).im;
        }
        ALM_Sol.alm2powspec(ClUpdateSol);
        Err /= (float) Nalm;
        Err = sqrt(Err);
          
        // Power spectrum deconvolution
        matmask_mult_cl(ClSol, MatMask, ClConv);
        // cout << "ok matmask_mult_cl " << endl;
        double ErrSpec=0.;
        
        fltarray MatErr;
        MatErr.alloc(ALM.Lmax()+1);
        for (int l=0; l <=  ALM.Lmax(); l++) 
        {
            ClResi.tt(l) = ClData.tt(l)   - ClConv.tt(l) ;
            MatErr(l) = ClResi.tt(l)*l*(l+1);
            ErrSpec += (ClResi.tt(l)*ClResi.tt(l))*l*(l+1);
        }
        if (Verbose == True) cout << "   ResiSpec = " << sqrt( ErrSpec / (float) ALM.Lmax() )  << endl;
        MatErr.info("SpecData-M#SpecSol");
        matmask_mult_cl(ClResi, MatMask, ClResi1, true);
        for (int l=0; l <=  ALM.Lmax(); l++) ClSol.tt(l) = MAX(0, ClSol.tt(l) + ClResi1.tt(l));

        extract_crosspowspec(ALM_M0, ALM_M1, ClM0M1);
        ALM_M0.alm2powspec(ClM0);
        ALM_M1.alm2powspec(ClM1);
        
        double Err1 = 0.;
        double Err2 = 0.;
        for (int l=0; l <=  ALM.Lmax(); l++) Err1 +=  (ClUpdateSol.tt(l) - ClSol.tt(l))* (ClUpdateSol.tt(l) - ClSol.tt(l));
        
          
        matmask_mult_cl(ClSol, RevMatMask, ClM0_Exp);
       
        for (int l=1; l <=  ALM_Rec.Lmax(); l++)
        {
            double Coef =(ClM0.tt(l) > 0) ? sqrt( ClM0_Exp.tt(l) / ClM0.tt(l)): 1.;
            if (Coef > 2) Coef = 2;
            for (int m=0; m <= l; m++) 
            {
                ALM_M0(l,m).re  *= Coef;
                ALM_M0(l,m).im  *= Coef;
                
            //    ALM_Sol(l,m).re = ALM(l,m).re + ALM_M0(l,m).re;
            //    ALM_Sol(l,m).im = ALM(l,m).im + ALM_M0(l,m).im;
            } 
        }
         
        
        ALM_Sol.alm2powspec(ClUpdateSol);
        for (int l=0; l <=  ALM.Lmax(); l++) Err2 +=  (ClUpdateSol.tt(l) - ClSol.tt(l))* (ClUpdateSol.tt(l) - ClSol.tt(l));
        ALM_Sol.alm2powspec(ClSol);
        // cout << "  Update Err1 = " << Err1 << ", Err2 = " << Err2 << endl;
        
        // ALM_Rec.info((char *) "ALM_Rec bef");

        // cout << FirstSoftThreshold << " " << LastSoftThreshold  << endl;
        Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * ((float) i / (float) Max_Iter)));  // 4 versus 2.8  in IDL code
        if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;	  
            
        if (i != 0) GaussCst_Lambda += GaussCst_StepLambda;
        
        if ((Verbose == True)  && (IsotropyCst == False)) cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
        if ((Verbose == True)  && (IsotropyCst == True)) cout << "Iter  " << i+1 << ", Lambda = " << GaussCst_Lambda  <<  " Err = " << Err << endl;
        
       
       // isotropy_constraint(ALM, PowSpecData, GaussCst_Lambda, UpdatePowSpec);

        if (IsotropyCst == True) 
        {
           if (InputPowSpec == true) isotropy_constraint(ALM_Sol, PowSpecData, GaussCst_Lambda, False);
           else isotropy_constraint(ALM_Sol, ClSol, GaussCst_Lambda, False);
        }
        else   // Alm L1 minimization
        { 
            float ThresholdLevel = SigmaNoise*Lambda;
            if (Lambda > 0)
            {
                if (HardThreshold == True) CptAlm = ALM_Sol.hard_threshold(ThresholdLevel, MaxNonZeroL);
                else if (SoftThreshold == True)   CptAlm = ALM_Sol.soft_threshold(ThresholdLevel*sqrt(2.), MaxNonZeroL);
                // if (Verbose == True) cout << "  Number of significant Alm : = " << CptAlm <<  " MaxL = " << MaxNonZeroL << endl; 
            }
            else  MaxNonZeroL = Lmax_Tot;
                        
             // To accelerate, we compute the Alm only to MaxNonZeroL*2
            if (Acceleration == True)  
            {
                if (2*MaxNonZeroL >  ALM_FIRST_L)  ALM_FIRST_L = 2*MaxNonZeroL;
                ALM_FIRST_L = MIN(ALM_FIRST_L, Lmax_Tot);
            }
        } // end l1 minimization
	    
         ALM_Sol.alm2powspec(ClSol);
     }
    
    for (int l=0; l <=  ALM.Lmax(); l++)
    for (int m=0; m <= l; m++) ALM(l,m) = ALM_Sol(l,m);

    ALM_Sol.alm_rec(Result);
    exit_inpainting();
}

#endif
/*********************************************************************/


/*********************************************************************/
