/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  / Florent Sureau Modified
**
**    Date:  25/10/08 / 04/10/11
**    
**    File:  MRS_Sparse.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Code for sparse representation
**    -------------------
**
******************************************************************************/


#include"MRS_Sparse.h"
#include"WT1D_FFT.h"

//FCS MODIFIED: adding C_UWT2D transforms for lGMCA inversion
void C_UWT2D_ATROUS::alloc(int NsideIn, int NScale, bool nested, type_border TBorder, int lnb_thr)
{
    nest = nested;
    Nside = NsideIn;
    if (NScale >= 2) NbrScale=NScale;
    else  NbrScale = (int) (log((double) NsideIn) / log(2.)+1);

    UWT2D->Bord=TBorder;
    WTTrans=new Hmap<REAL> [NbrScale];
    for (int p=0;p<NbrScale;p++) (WTTrans[p]).alloc(Nside, nest);
    AllocMem=1;
    nb_thr=lnb_thr;
}

void C_UWT2D_ATROUS::transform(Hmap<REAL> & DataIn)
{
	int f,p;
	fltarray Face,FaceW;
	Ifloat IFace, *TabTrans;
	TabTrans=new Ifloat[NbrScale];
	
	Face.alloc(Nside,Nside);
	IFace.alloc(Face.buffer(), Face.ny(), Face.nx());
	for (p=0;p<NbrScale;p++) TabTrans[p].alloc(Face.ny(), Face.nx());

	//#pragma omp parallel for default(none) private(f,p,Face,IFace,TabTrans,FaceW) 
    for (f=0; f < 12; f++){	
    	DataIn.get_one_face(Face, f);
		UWT2D->transform(IFace,TabTrans,NbrScale,nb_thr);
        for (p=0;p<NbrScale;p++) {
        	FaceW.alloc(TabTrans[p].buffer(),TabTrans[p].ny(),TabTrans[p].nx());
        	(WTTrans[p]).put_one_face(FaceW, f);
        }
    }
    delete [] TabTrans;
}

void C_UWT2D_ATROUS::recons(Hmap<REAL> & DataOut)
{
    fltarray *Face,FaceR;
	Ifloat *IFace, TabTrans;
	int p,f;
    
    Face= new fltarray[NbrScale];
    for (p=0;p<NbrScale;p++) (Face[p]).alloc(nside(),nside());
    IFace= new Ifloat[NbrScale];
    for (p=0;p<NbrScale;p++) (IFace[p]).alloc((Face[p]).buffer(),(Face[p]).ny(),(Face[p]).nx()); //map same region as Face, but Ifloat vs fltarray
	FaceR.alloc(TabTrans.buffer(),TabTrans.ny(),TabTrans.nx());


 	//#pragma omp parallel for default(none) private(f,p) shared (NbrScale)
	 for (int f=0; f < 12; f++) {
    	for (p=0;p<NbrScale;p++) (WTTrans[p]).get_one_face(Face[p], f);    	//For each wavelet scale, get the face
        UWT2D->recons(IFace,TabTrans,NbrScale,nb_thr); //perform reconstruction
        DataOut.put_one_face(FaceR, f); //put the result in DataOut
    }
    delete [] Face;
    delete [] IFace;
}
//End FCSMODIFIED C_UWT2D

void C_OWT::alloc(int NsideIn, int NScale, bool nested)
{
    int x=0;
    nest = nested;
    Nside = NsideIn;
    WTTrans.alloc(NsideIn, nested);
    T_FilterBank = F_MALLAT_7_9; 
    FAS = new FilterAnaSynt;
    FAS->alloc(T_FilterBank);
    SBF = new SubBandFilter(*FAS, NORM_L2);
    SBF->setBorder(I_MIRROR);
    
    if (NScale >= 2) NbrScale=NScale;
    else  NbrScale = (int) (log((double) NsideIn) / log(2.)+1);
    OWT2D = new Ortho_2D_WT(*SBF);
}

void C_OWT::transform(Hmap<REAL> & DataIn)
{
   fltarray Face;
   Ifloat IF;
   Face.alloc(Nside,Nside);
   IF.alloc(Face.buffer(), Face.ny(), Face.nx());
   //cout << "Transform " << endl;
    for (int f=0; f < 12; f++)
    {
       DataIn.get_one_face(Face, f);
       //cout << "  get " << NbrScale << endl;
       OWT2D->transform(IF,NbrScale);
       //cout << "  put " << endl;
        WTTrans.put_one_face(Face, f);
       //cout << "  ok " << endl;
    }
    //cout << "transform ok " << endl;
}

void C_OWT::recons(Hmap<REAL> & DataOut)
{
    fltarray Face;
    Ifloat IF;
    Face.alloc(nside(),nside());
    IF.alloc(Face.buffer(), Face.ny(), Face.nx());
    
    for (int f=0; f < 12; f++)
    {
        WTTrans.get_one_face(Face, f);
        OWT2D->recons(IF,NbrScale);
        DataOut.put_one_face(Face, f);
    }
}




//FCS ADDED: Pyramidal wavelet transform
void wp_trans(Hdmap & Map, dblarray **TabCoef,intarray &NsideBand, fltarray & WP_W, fltarray & WP_H,  int NbrWP_Band, int Lmax,bool BandLimit, bool SqrtFilter, int NScale, int ALM_iter) //FCS Added
{
    int Nside = Map.Nside();
    long Npix;
    Hdmap Band,Result;
    
    if(BandLimit==false) Result=Map; //FCS Added First wavelet scale: difference between input image and first approximation scale
    CAlmR  ALM;
    ALM.alloc(Nside, Lmax);
    ALM.Norm = False;
    ALM.UseBeamEff = False;
    ALM.set_beam_eff(Lmax, Lmax);
    if((NScale==0)||(NScale>NbrWP_Band))  NScale=NbrWP_Band;
    *TabCoef=new dblarray[NScale+1];
    ALM.Niter=ALM_iter;//Nb iterations in ALM transform
    ALM.alm_trans(Map);
    for (int b=0; b <= NScale; b++) {
        printf("Process band %d\n",b);
        // int LMin = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,0): WP_W(0,0);
        CAlmR ALM_Band;
        int LMax;
        if (b == NScale) LMax=(int)WP_W(NbrWP_Band-b,1);
        else if (b ==0) LMax=Lmax;
        else LMax=(int)WP_W(NbrWP_Band-b,1);
        Npix=(unsigned long) NsideBand(b) * NsideBand(b) * 12ul;
        (*TabCoef)[b].alloc(Npix);
        ALM_Band.alloc(NsideBand(b), LMax);//Beware, using ring weighting implies that the ALM are correctly initialized
        ALM_Band.SetToZero();
        ALM_Band.Norm = False;
        ALM_Band.UseBeamEff = False;
        ALM_Band.set_beam_eff(LMax, LMax);
        // Compute the solution at a given resolution (convolution with H)
        int FL = MIN(ALM_Band.Lmax(),WP_H.nx()-1);
        if (FL > ALM.Lmax()) FL = ALM.Lmax();
        if (SqrtFilter == false)
        {
            if((b==0)&&(BandLimit==false))
            {
                  Result.SetNside ((int) NsideBand(b),  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
                ALM.alm_rec(Result,true,NsideBand(b));
            }
            //Modified such that we can use pyramidal transform : substraction in ALM SPACE
            //Note that using filter banks results in lower precision for coarsest (and more energetic) scales
            if(b == NScale)
            {
                for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++)
                           ALM_Band(l,m) = ALM(l,m) * (REAL) WP_H(l, NbrWP_Band-NScale);
             } else if (b==0) {
                   for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++)
                           ALM_Band(l,m) = ALM(l,m) * (REAL) (1L - (double) WP_H(l, NbrWP_Band-1));
             } else {
                 for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++)
                           ALM_Band(l,m) = ALM(l,m) * (REAL) (WP_H(l, NbrWP_Band-b)-WP_H(l, NbrWP_Band-1-b));
             }
             ALM_Band.alm_rec(Band,false,NsideBand(b));
             if((b==0)&&(BandLimit==false)) for (int p=0; p < Band.Npix(); p++) Band[p]+=Map[p]-Result[p];//FCS Modified: First band is now difference between band-limited Image at Lmax and Lowpass
         } else { //FCS Added: Sqrt filters
             if(b == NScale) {
                 for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL) sqrt(WP_H(l,NbrWP_Band-NScale));
             } else if (b == 0) {
                    for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL)  sqrt(1.-WP_H(l, NbrWP_Band-1));
             } else {
                   for (int l=0; l <= FL; l++)
                       for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL)  sqrt(WP_H(l, NbrWP_Band-b)-WP_H(l, NbrWP_Band-1-b));
             }
             ALM_Band.alm_rec(Band,false,NsideBand(b));
             if((b==0)&&(BandLimit==false)) {//FCS Modified: add Info at multipoles > lmax
                 Result.SetNside ((int) NsideBand(b),  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
                   ALM.alm_rec(Result,true);
                for (int p=0; p < Band.Npix(); p++){
                     Result[p]=Map[p]-Result[p];
                     Band[p]+= Result[p];
                }
             }
        } // endelse if (b == NbrWP_Band)
        for (int p=0; p < Band.Npix(); p++) ((*TabCoef)[b])(p) = Band[p];
    }
    // fits_write_fltarr("xxcoef.fits", TabCoef);

}



/*ooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooo*/
//Undecimated transforms
/****************************************************************************/

//General functions

void wp_trans(Hdmap & Map, dblarray & TabCoef, fltarray & WP_W, fltarray & WP_H, int NbrWP_Band, int Lmax,bool BandLimit, bool SqrtFilter, int NScale, int ALM_iter) //FCS Added 
{
    int Nside = Map.Nside();
    Hdmap Band,Result;
    
    Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    Result.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    if(BandLimit==false) Result=Map; //FCS Added First wavelet scale: difference between input image and first approximation scale
    CAlmR  ALM,ALM_Band;  
    ALM.alloc(Nside, Lmax);
    ALM_Band.alloc(Nside, Lmax);
    ALM.Norm = ALM_Band.Norm = False;
    ALM.UseBeamEff = ALM_Band.UseBeamEff = False;
    ALM.set_beam_eff(Lmax, Lmax);
    ALM_Band.set_beam_eff(Lmax, Lmax);
         
    if((NScale==0)||(NScale>NbrWP_Band))  NScale=NbrWP_Band;
    TabCoef.resize(Map.Npix(),NScale+1);
    ALM.Niter=ALM_iter;//Nb iterations in ALM transform
    ALM.alm_trans(Map);
    
    
    for (int b=0; b <= NScale; b++)
    {  
        // int LMin = (b != NbrWP_Band) ? WP_W(NbrWP_Band-1-b,0): WP_W(0,0);
        int LMax;
        if (b == NScale) LMax=(int)WP_W(NbrWP_Band-b,1);
        else if (b ==0) LMax=Lmax;
        else LMax=(int)WP_W(NbrWP_Band-b,1);
        
       // if (Verbose == True) 
         cout << "        WP:Band " << b+1 << ",  Lmax = " << LMax << endl;
        
        if ((b == NScale) && (SqrtFilter == false)) Band = Result;
        else 
        {
            ALM_Band.Set(LMax,  LMax);
            ALM_Band.SetToZero();
                
            // Compute the solution at a given resolution (convolution with H)
            int FL = MIN(ALM_Band.Lmax(),WP_H.nx()-1);
            if (FL > ALM.Lmax()) FL = ALM.Lmax();
            if (SqrtFilter == false) 
            {
                for (int l=0; l <= FL; l++) 
                for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL) WP_H(l, NbrWP_Band-1-b);
                ALM_Band.alm_rec(Band);
                if((b==0)&&(BandLimit==true)) ALM.alm_rec(Result);//FCS Modified: First band is now difference between band-limited Image at Lmax and Lowpass
                
                // Compute the coefficients and standard deviations in and out the mask
                for (int p=0; p < Band.Npix(); p++) 
                {
                   	double Val = Band[p];
                   	Band[p] = Result[p] - Val;
                   	Result[p] = Val;
                }
             } 
             else 
             { //FCS Added: Sqrt filters 
                if(b == NScale) 
                {
                	for (int l=0; l <= FL; l++) 
                   		for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL) sqrt(WP_H(l,NbrWP_Band-NScale));
                } 
                else if (b == 0) 
                {
                	 for (int l=0; l <= FL; l++) 
                   		for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL) sqrt(1.-WP_H(l, NbrWP_Band-1));
                } 
                else 
                {
                	for (int l=0; l <= FL; l++) 
                   		for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL)  sqrt(WP_H(l, NbrWP_Band-b)-WP_H(l, NbrWP_Band-1-b));
                }                	
                ALM_Band.alm_rec(Band);
                if((b==0)&&(BandLimit==false)) {//FCS Modified: add Info at multipoles > lmax
                	ALM.alm_rec(Result,true);
                	for (int p=0; p < Band.Npix(); p++){
                		 Result[p]=Map[p]-Result[p];
                		 Band[p]+= Result[p];
                	}
                }
             }
        } // endelse if (b == NbrWP_Band)
        for (int p=0; p < Band.Npix(); p++) TabCoef(p, b) = Band[p];
    }
    // fits_write_fltarr("xxcoef.fits", TabCoef);
}

inline void wp_trans(Hdmap & Map, dblarray & TabCoef, fltarray & WP_W, fltarray & WP_H, int NbrWP_Band, int Lmax){wp_trans(Map,TabCoef,WP_W,WP_H,NbrWP_Band,Lmax,false,false,0,0);} //FCS Added
inline void wp_trans(Hdmap & Map, dblarray **TabCoef,intarray &NsideBand, fltarray & WP_W, fltarray & WP_H, int NbrWP_Band, int Lmax){wp_trans(Map,TabCoef,NsideBand,WP_W,WP_H,NbrWP_Band,Lmax,false,false,0,0);} //FCS Added


/****************************************************************************/


void get_wt_bspline_filter(int nside, fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int  Lmax, int NbrBand, bool TightFrame)
{
   WT1D_FFT W;
   WT1D_FFT WTight(WT_PYR_FFT_DIFF_SQUARE);
   //  enum type_wt_filter {WT_PHI,WT_PSI,WT_H,WT_H_TILDE,WT_G,WT_G_TILDE};

   type_wt_filter Filter;
   int LM = Lmax+1;
   WP_WPFilter.alloc(LM, NbrBand);
   int Nstep = NbrBand - 1;
   fltarray TabHH;
   
    TabH.alloc(LM, Nstep);
    TabHH.alloc(LM, Nstep);
    TabG.alloc(LM, Nstep);
    Win.alloc(Nstep,2);

    int b,p,Np = 2*LM;
    for (b=0; b < Nstep; b++)
    {       
        // cout << Np << endl;
        for (p=0; p < LM; p++)  
        {
            if (TightFrame)
            {
                TabH(p, Nstep-1-b) = WTight.filter(WT_H, p, Np);
                TabG(p, Nstep-1-b) = WTight.filter(WT_G, p, Np);
            }
            else
            {
                TabH(p, Nstep-1-b) = W.filter(WT_H, p, Np);
                TabG(p, Nstep-1-b) = W.filter(WT_G, p, Np);
            }
           if (b > 0) TabHH(p, b) =  TabH(p, Nstep-1-b) * TabHH(p,b-1);
           else TabHH(p, 0) =  TabH(p, Nstep-1-b) ;
        }
        Np /= 2;
    }
 
    Np = LM;
    for (b=0; b < Nstep; b++) 
    {
        Win(Nstep-1-b,0) = 0;
        Win(Nstep-1-b,1) = Np;
        Np /= 2;
    }

    if (TightFrame)
    {
        // cout << " TightFrame " << endl;
        for (int l=0; l < LM; l++) WP_WPFilter(l,0) = TabG(l,NbrBand-2);
    }
    else
        for (int l=0; l < LM; l++) WP_WPFilter(l,0) = 1. - TabH(l,NbrBand-2);
    
    for (int b=1; b < Nstep; b++)
    for (int l=0; l < LM; l++)  WP_WPFilter(l,b) = TabHH(l,b-1) * TabG(l,Nstep-b-1);
    b = Nstep;
    for (int l=0; l < LM; l++)  WP_WPFilter(l,b) = TabHH(l,b-1);

/*
    fits_write_fltarr("xx_hh.fits", TabHH);
    cout << " Win = " << Lmax << " " << NbrBand << endl;
    fits_write_fltarr("xx_h.fits", TabH);
    fits_write_fltarr("xx_g.fits", TabG);
    fits_write_fltarr("xx_w2.fits", Win);
    fits_write_fltarr("xx_filters.fits", WP_WPFilter);
*/
     // exit(-1);
}

/****************************************************************************/

void mrs_wt_trans(Hdmap & Map, dblarray & TabCoef, fltarray & TabFilter,  int Lmax, int NScale, int ALM_iter)
{
    int Nside = Map.Nside();
    Hdmap Band,SmoothMap;
    
    // cout << "mrs_wt_transL Lmax =  " << Lmax << ", Nside = "  << Nside << ", Nscale = " << NScale << endl;
    Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    SmoothMap.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    SmoothMap=Map;  
    
    CAlmR  ALM;  
    ALM.alloc(Nside, Lmax);
    ALM.Norm = False;
    ALM.UseBeamEff =  False;
    ALM.set_beam_eff(Lmax, Lmax);
     
    TabCoef.resize(Map.Npix(),NScale);
    ALM.Niter=ALM_iter;//Nb iterations in ALM transform

    int Nstep = NScale - 1;
    for (int p=0; p < Band.Npix(); p++)  Band[p] =  Map[p];
    // cout << "ALM =  "<< endl;
    ALM.alm_trans(Map);
    // Map.info((char*) "Input MAP");
    // fits_write_fltarr("xxf3.fits", TabFilter);

    for (int b=0; b < Nstep; b++)
    {  
        CAlmR ALM_Band;
        // cout << "AllocBand =  "<< b+1 << endl;

        ALM_Band.alloc(Nside, Lmax);
        // cout << "BandAlloc ok  "<< b+1 << endl;
        ALM_Band.SetToZero();
        // cout << "   SetToZero =  "<< b+1 << endl;

        ALM_Band.Norm = False;
        ALM_Band.UseBeamEff = False;
    	ALM_Band.set_beam_eff(Lmax, Lmax);  
         
        // cout << "        WP:Band " << b+1 << " " << NScale-2-b << ",  Lmax = " << Lmax << " TabFilter.nx = " << TabFilter.nx() << endl;
        
        // Compute the solution at a given resolution (convolution with H)
     
        for (int l=0; l <= Lmax; l++) 
        for (int m=0; m <= l; m++) ALM_Band(l,m) = ALM(l,m) * (REAL) TabFilter(l, Nstep-1-b);
        ALM_Band.alm_rec(SmoothMap);
        // SmoothMap.info((char*) "SmoothMap ");
        for (int p=0; p < Band.Npix(); p++) Band[p] = Band[p] - SmoothMap[p];
        // Band.info((char*) "Band ");

        for (int p=0; p < Band.Npix(); p++)  
        {
           TabCoef(p, b) =  Band[p];
           Band[p] = SmoothMap[p];
        }
     }
     int b = Nstep;
     for (int p=0; p < Band.Npix(); p++)  TabCoef(p, b) =  SmoothMap[p];

     
    // fits_write_fltarr("xxcoef.fits", TabCoef);
}
/****************************************************************************/


    // fits_write_fltarr("xxcoef.fits", TabCoef);
/****************************************************************************/
//Undecimated class
void C_UWT::wt_alloc(int NsideIn, int NScale, int LM, bool nested, bool TFrame)
{
    nest = nested;
    Nside = NsideIn;
    Lmax=LM;
    TightFrame=TFrame;
    
    T_Fil = F_ALM_SPLINE;
    MeyerWT = false;
    if (T_Fil == F_ALM_MEYER)
    {
        wp_alloc(NsideIn, LM, nested);
    }
    else // SPLINE Filter
    {
        if (NScale >= 2) NbrScale=NScale;
        else  NbrScale = (int) (log((double) NsideIn) / log(2.)-1);
        NpixPerBand = NsideIn*NsideIn*12;

        // WTTrans = new Hmap<REAL> [NbrScale];
         // for (int b=0; b < NbrScale; b++) (WTTrans[b].alloc)(NsideIn, nested);
        WTTrans.alloc(NpixPerBand, NbrScale);
        
        get_wt_bspline_filter(Nside, WP_H, WP_G, WP_W, WP_WPFilter, Lmax, NbrScale,TightFrame);
        
        double EnerBand;
        TabNorm.alloc(NbrScale);
        TabMad.alloc(NbrScale);
        for (int b=0; b < NbrScale; b++)
        {
            EnerBand = 0.;
            for (int l=0; l <= MIN(Lmax, WP_WPFilter.nx()-1); l++) 
            {
                EnerBand += (2.*l+1)*  WP_WPFilter(l, b)* WP_WPFilter(l, b);
            }
            TabNorm(b)  = sqrt( EnerBand  /  (double) NpixPerBand);
            if (Verbose) cout << "Band " << b+1 << ", Norm = " << TabNorm(b) << endl;
        }
    }
}

/****************************************************************************/

void C_UWT::wp_alloc(int NsideIn, int LM, bool nested)
{
    nest = nested;
    Nside = NsideIn;
    Lmax=LM;
	
    MeyerWT = true;
    if (All_WP_Band == False) get_wp_meyer_filter(Nside, WP_H, WP_G, WP_W, WP_WPFilter, Lmax);
    else get_planck_wp_meyer_filter(WP_H, WP_G, WP_W, WP_WPFilter, Lmax, Lmax);
    
        fits_write_fltarr("xx_h.fits", WP_H);
        fits_write_fltarr("xx_g.fits", WP_G);
        fits_write_fltarr("xx_filters.fits", WP_WPFilter);
    
    NbrScale = WP_H.ny();
    // cout << " NN = " << NbrScale << " Lmax = " << Lmax << ", WP_WPFilter:  " << WP_WPFilter.nx() << " " <<  WP_WPFilter.ny() << endl;
    NpixPerBand = NsideIn*NsideIn*12;
    WTTrans.alloc(NpixPerBand, NbrScale);
    // WTTrans = new Hmap<REAL> [NbrScale];
    // for (int b=0; b < NbrScale; b++) (WTTrans[b].alloc)(NsideIn, nested);
    
    // double CoefN = 1.;
    double EnerBand;
    TabNorm.alloc(NbrScale);
    TabMad.alloc(NbrScale);
    for (int b=0; b < NbrScale; b++)
    {
        EnerBand = 0.;
        for (int l=0; l <= MIN(Lmax, WP_WPFilter.nx()-1); l++) 
		{
            EnerBand += (2.*l+1)*  WP_WPFilter(l, b)* WP_WPFilter(l, b);
        }
		TabNorm(b)  = sqrt( EnerBand  /  (double) NpixPerBand);
       // cout << "Band " << b+1 << ", Norm = " << TabNorm(b) << endl;
    }
}  
 
/****************************************************************************/

void C_UWT::transform(Hmap<REAL> & DataIn, bool BandLimit)
{
    if ((MeyerWT == false) && (TightFrame == false))
            mrs_wt_trans(DataIn, WTTrans, WP_H, Lmax, NbrScale, ALM_iter);
    else transform(DataIn, BandLimit, TightFrame, NbrScale);
        // wp_trans(DataIn, WTTrans, WP_W, WP_H,  NbrScale, Lmax);
}


/****************************************************************************/

void C_UWT::transform(Hmap<REAL> & DataIn, bool BandLimit, bool SqrtFilter, int NScale)
{
    if (MeyerWT == false)
    {
        if (TightFrame == false) mrs_wt_trans(DataIn, WTTrans, WP_H, Lmax, NbrScale, ALM_iter);
        else
        {
            // cout << "TRANS TIGHT: " << WP_W.nx() << " " << WP_H.nx() << " Lmax = " << Lmax << endl;
            Hdmap Band,Result;
            Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
            Result.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
            if(BandLimit==false) Result=DataIn;
            CAlmR  ALM,ALM_Band_H, ALM_Band_G;
            ALM.alloc(Nside, Lmax);
            ALM_Band_H.alloc(Nside, Lmax);
            ALM_Band_G.alloc(Nside, Lmax);
            ALM.Norm = ALM_Band_H.Norm = ALM_Band_G.Norm= False;
            ALM.UseBeamEff = ALM_Band_H.UseBeamEff = ALM_Band_G.UseBeamEff =False;
            ALM.set_beam_eff(Lmax, Lmax);
            ALM_Band_H.set_beam_eff(Lmax, Lmax);
            ALM_Band_G.set_beam_eff(Lmax, Lmax);

            WTTrans.resize(DataIn.Npix(),NbrScale);
            ALM.Niter= ALM_Band_H.Niter = ALM_Band_G.Niter= ALM_iter; //Nb iterations in ALM transform
            ALM.alm_trans(DataIn);
            
            int LMax=Lmax;
            int Nstep=NbrScale-1;
            
            if (BandLimit==false)
            {// Keep Info at multipoles > lmax
                ALM.alm_rec(Result,true);
                for (int p=0; p < Band.Npix(); p++)
                      Result[p]=DataIn[p]-Result[p];
            }
            
            for (int b=0; b < Nstep; b++)
            {
                 if (Verbose == True)
                   cout << "        WT:Band " << b+1 << ",  Lmax = " << LMax << endl;
                
                 ALM_Band_H.Set(LMax,  LMax);
                 ALM_Band_H.SetToZero();
                 ALM_Band_G.Set(LMax,  LMax);
                 ALM_Band_G.SetToZero();
                // Compute the solution at a given resolution (convolution with H)
                int FL = MIN(ALM_Band_H.Lmax(),WP_H.nx()-1);
                if (FL > ALM.Lmax()) FL = ALM.Lmax();
                for (int l=0; l <= FL; l++)
                for (int m=0; m <= l; m++)
                {
                    ALM_Band_G(l,m) = ALM(l,m) * (REAL) WP_G(l,Nstep-1-b);
                    ALM_Band_H(l,m) = ALM(l,m) * (REAL) WP_H(l,Nstep-1-b);
                    ALM(l,m) = ALM_Band_H(l,m);
                }
                ALM_Band_G.alm_rec(Band);
                for (int p=0; p < Band.Npix(); p++) WTTrans(p, b) = Band[p];
            }
            ALM_Band_H.alm_rec(Band);
            for (int p=0; p < Band.Npix(); p++) WTTrans(p, NbrScale-1) = Band[p];
            
            if (BandLimit==false)
            {// add Info at multipoles > lmax
                 for (int p=0; p < Band.Npix(); p++)
                      WTTrans(p, 0) += Result[p];
            }
        }
    }
    else wp_trans(DataIn, WTTrans, WP_W, WP_H,  NbrScale, Lmax,BandLimit,SqrtFilter,NScale,ALM_iter);
}

/****************************************************************************/

void C_UWT::recons(Hmap<REAL> & DataOut, bool BandLimit) {
	// for (int b=0; b <  NbrScale; b++)
	
	NbrScale=WTTrans.ny();
    
    if (TightFrame == false)
    {
        DataOut.fill(0.);
        for (int b=0; b <  NbrScale; b++)
            for (int p=0; p < NpixPerBand; p++) DataOut[p] += WTTrans(p, b);
    }
    else
    {
        // cout << "Tight Frame WT reconstruction" << endl;
        int NScale=WTTrans.ny()-1;
        int NS = sqrt(WTTrans.nx()/12l);
        if (NS != Nside)
        {
            cout << "Error: nside = " << NS << ", expected nside = " <<Nside << endl;
        }
        unsigned long Npix=WTTrans.nx();
        Hdmap Band,Rec;
        Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
        DataOut.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
        Rec.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
        DataOut.fill(0.);
        CAlmR  ALM,ALM_Band;
        ALM.alloc(Nside, Lmax);
        ALM_Band.alloc(Nside, Lmax);
        ALM.Norm = ALM_Band.Norm = False;
        ALM.UseBeamEff = ALM_Band.UseBeamEff = False;
        ALM.SetToZero();
        ALM_Band.set_beam_eff(Lmax, Lmax);
        ALM_Band.Set(Lmax,  Lmax);
        ALM_Band.SetToZero();
        
        // Coarsest scale
        int LM = WP_W(0,1);
        // cout << "Coarse : Tight Frame WT reconstruction: lmax = " << LM << endl;
        ALM_Band.alloc(Nside, LM);
        ALM_Band.set_beam_eff(LM, LM);
        ALM_Band.SetToZero();
        ALM_Band.Niter=ALM_iter;//Nb iterations in ALM transform
        int b=NScale;
        for(unsigned long p=0;p<Npix;p++) Band[p]=WTTrans(p,b);
        // cout << "ALMtrans"<< endl;
        ALM_Band.alm_trans(Band);
        for (int l=0; l <= LM; l++)
        for (int m=0; m <= l; m++) ALM (l,m) = ALM_Band(l,m);
        // cout << "loop"<< WP_W.nx() << " " <<   WP_W.ny() << endl;
        // fits_write_fltarr("xxw.fits",WP_W);
        for (int i=0; i < NScale; i++)
        {
            b=NScale-i-1;
            LM = WP_W(i,1);
            // cout << i << ": Band " << b+1 << ": Tight Frame WT reconstruction: lmax = " << LM << endl;
            ALM_Band.alloc(Nside, LM);
            ALM_Band.set_beam_eff(LM, LM);
            ALM_Band.SetToZero();
            ALM_Band.Niter=ALM_iter;//Nb iterations in ALM transform
            for(unsigned long p=0;p<Npix;p++) Band[p]=WTTrans(p,b);
            ALM_Band.alm_trans(Band);
            for (int l=0; l <= LM; l++)
            for (int m=0; m <= l; m++)
                   ALM(l,m) =  (ALM(l,m) * (REAL) WP_H(l,NbrScale-NScale)) + (ALM_Band(l,m) * (REAL) WP_G(l,NbrScale-NScale));
        }
        ALM.alm_rec(DataOut);
        if(BandLimit==false)
        {
            ALM_Band.alm_rec(Rec);
            for (int p=0; p < Band.Npix(); p++) Rec[p]= Band[p] - Rec[p];
            for (int p=0; p < DataOut.Npix(); p++) DataOut[p]+=Rec[p];
        }
    }

        
//        for (int b=0; b <= NScale; b++)
//        {//NbrWP_Band
//            cout << "band " << b << endl;
//
//            for(unsigned long p=0;p<Npix;p++) Band[p]=WTTrans(p,b);
//            int LMax;
//            if(b==NScale) LMax=(int) WP_W(NbrScale-NScale,1);
//            else if (b==0) LMax=Lmax;
//            else LMax=(int) WP_W(NbrScale-b,1);
//            ALM.alloc(Nside, LMax);
//            ALM.set_beam_eff(LMax, LMax);
//            ALM.SetToZero();
//            ALM.Niter=ALM_iter;//Nb iterations in ALM transform
//            ALM.alm_trans(Band);
//            int FL = MIN(ALM.Lmax(),WP_H.nx()-1);
//            if (FL > ALM.Lmax()) FL = ALM.Lmax();
//            if(b==NScale)
//            {//Low pass filter
//                for (int l=0; l <= FL; l++)
//                    for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) WP_H(l,NbrScale-NScale));
//            } else {
//                for (int l=0; l <= FL; l++)
//                    for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) WP_WPFilter(l,b));
//            }
//
//            if((b==0)&&(BandLimit==false)) {//FCS Modified: add Info at multipoles > lmax
//                ALM.alm_rec(Rec);
//                for (int p=0; p < Band.Npix(); p++) Rec[p]= Band[p] - Rec[p];
//            }
//        }
//        ALM_Band.alm_rec(DataOut);
//        if(BandLimit==false) for (int p=0; p < DataOut.Npix(); p++) DataOut[p]+=Rec[p];
//    }
}

/****************************************************************************/

void C_UWT::recons(Hmap<REAL> & DataOut, bool BandLimit, bool SqrtFilter, int NScale) //FCS added
{
    DataOut.fill(0.);

	if(SqrtFilter==false) recons(DataOut);
	else {
		if(NScale==0) NScale=WTTrans.ny()-1;
		int Nside = sqrt(WTTrans.nx()/12l);
		unsigned long Npix=WTTrans.nx();
    	Hdmap Band,Rec;    	
    	Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
		DataOut.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
		DataOut.fill(0.);
    	CAlmR  ALM,ALM_Band;  
    	ALM_Band.alloc(Nside, Lmax);
    	ALM.Norm = ALM_Band.Norm = False;
    	ALM.UseBeamEff = ALM_Band.UseBeamEff = False;
    	ALM_Band.set_beam_eff(Lmax, Lmax); 
        ALM_Band.Set(Lmax,  Lmax);
        ALM_Band.SetToZero();
    	for (int b=0; b <= NScale; b++) {//NbrWP_Band
    		for(unsigned long p=0;p<Npix;p++) Band[p]=WTTrans(p,b);
    		int LMax;
    		if(b==NScale) LMax=(int) WP_W(NbrScale-NScale,1);
    		else if (b==0) LMax=Lmax;
    		else LMax=(int) WP_W(NbrScale-b,1);
    		ALM.alloc(Nside, LMax);
    		ALM.set_beam_eff(LMax, LMax);
    		ALM.SetToZero();
        	ALM.Niter=ALM_iter;//Nb iterations in ALM transform
    		ALM.alm_trans(Band);
            int FL = MIN(ALM.Lmax(),WP_H.nx()-1);
            if (FL > ALM.Lmax()) FL = ALM.Lmax();
            if(b==NScale) {//Low pass filter
            	for (int l=0; l <= FL; l++) 
                	for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) sqrt(WP_H(l,NbrScale-NScale)));
            } else {
            	for (int l=0; l <= FL; l++) 
                	for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) sqrt(WP_WPFilter(l,b)));
            }
            if((b==0)&&(BandLimit==false)) {//FCS Modified: add Info at multipoles > lmax
            	ALM.alm_rec(Rec);
            	for (int p=0; p < Band.Npix(); p++) Rec[p]= Band[p] - Rec[p];
        	}    
    	}
    	ALM_Band.alm_rec(DataOut);
        if(BandLimit==false) for (int p=0; p < DataOut.Npix(); p++) DataOut[p]+=Rec[p];    	
    }
    // fits_write_fltarr("xxcoef.fits", TabCoef);
} 

/****************************************************************************/

void C_UWT::hard_thresholding(int b, float NSigma, float & SigmaNoise, bool UseMad)
{
    float Level = SigmaNoise * NSigma;
    if (UseMad == true)
    {
       fltarray Tab;
       Tab.alloc(WTTrans.nx());
       for (int i=0; i < Tab.nx(); i++) Tab(i) = WTTrans(i,b);
       SigmaNoise = get_sigma_mad(Tab.buffer(), Tab.n_elem() );
       Level = SigmaNoise * NSigma;
       TabMad(b) = SigmaNoise;
    }
    else Level = SigmaNoise * NSigma * TabNorm(b);
     
    for (int i=0; i < WTTrans.nx(); i++)
    {
       if (ABS(WTTrans(i,b)) < Level) WTTrans(i,b) = 0;
    }
}

/****************************************************************************/

void C_UWT::set_band(int b, float Value)
{
    for (int i=0; i < WTTrans.nx(); i++) WTTrans(i,b) = Value;
}

/****************************************************************************/

void C_UWT::hard_thresholding(Hmap<REAL> & DataIn, float NSigma, float & SigmaNoise, bool UseMad, bool KillLastScale, int FirstDetectScale)
{
    transform(DataIn);
    for (int b=0; b < NbrScale-1; b++) hard_thresholding(b,  NSigma,  SigmaNoise,  UseMad);
    if (KillLastScale == true) set_band(NbrScale-1, 0.);
    for (int b=0; b < FirstDetectScale; b++) set_band(b, 0.);
    // fits_write_dblarr("xx_wt.fits", WTTrans );
    recons(DataIn);
}


/*ooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooo*/
//Pyramidal transforms
/****************************************************************************/

void C_PWT::wp_alloc(int NsideIn, int LM, bool nested)
{
    int LMax,b;
    nest = nested;
    Nside = NsideIn;
    Lmax=LM;

    if (All_WP_Band == False) get_wp_meyer_filter(Nside, WP_H, WP_G, WP_W, WP_WPFilter, Lmax);
    else get_planck_wp_meyer_filter(WP_H, WP_G, WP_W, WP_WPFilter, Lmax, Lmax);
    NbrScale = WP_H.ny();
    NsidePerBand.alloc(NbrScale+1);
    NpixPerBand.alloc(NbrScale+1);
    for(b=0;b<=NbrScale;b++) {
        if (b == NbrScale) LMax=(int)WP_W(0,1);
        else if (b ==0) LMax=Lmax;
        else LMax=(int)WP_W(NbrScale-b,1);
        NsidePerBand(b)=max((int) pow(2.,ceil(log( (double) LMax)/log(2.))-1.),64);//Compute the Nside parameter (power of two such that 2*Nside > Lmax)
        NpixPerBand(b)=(unsigned long) NsidePerBand(b) * NsidePerBand(b) * 12ul;
    }
}  


/****************************************************************************/
void C_PWT::transform(Hmap<REAL> & DataIn, bool BandLimit, bool SqrtFilter, int NScale)
{
   wp_trans(DataIn, &PWTTrans,NsidePerBand, WP_W, WP_H,  NbrScale, Lmax,BandLimit,SqrtFilter,NScale,ALM_iter);
   WTalloc=true;
}

/****************************************************************************/
void C_PWT::recons(Hmap<REAL> & DataOut, bool BandLimit, bool SqrtFilter, int NScale) //FCS added
{
   	if(NScale==0) {
   		printf("Should specify the number of scales for Pyramidal wavelet transform\n");
   		exit(-1);	
   	}
   	Hdmap Band,Rec;    	
	DataOut.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
	DataOut.fill(0.);
   	CAlmR  ALM,ALM_Band;  
   	ALM_Band.alloc(Nside, Lmax);//Beware of ring weighting (Nside and not NsidePerBand(b))
   	ALM.Norm = ALM_Band.Norm = False;
   	ALM.UseBeamEff = ALM_Band.UseBeamEff = False;
   	ALM_Band.set_beam_eff(Lmax, Lmax); 
    ALM_Band.Set(Lmax,  Lmax);
    ALM_Band.SetToZero();
   	for (int b=0; b <= NScale; b++)
    {//NbrWP_Band
   	    Band.SetNside ((int) NsidePerBand(b),  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
   		for(long p=0;p<Band.Npix();p++) Band[p]=(PWTTrans[b])(p);
   		int LMax;
   		if(b==NScale) LMax=(int) WP_W(NbrScale-NScale,1);
   		else if (b==0) LMax=Lmax;
   		else LMax=(int) WP_W(NbrScale-b,1);
   		ALM.alloc( NsidePerBand(b), LMax);//Beware of ring weighting (NsidePerBand(b) and not Nside)
   		ALM.set_beam_eff(LMax, LMax);
   		ALM.SetToZero();
   		ALM.Niter=ALM_iter;//Nb iterations in ALM transform
   		ALM.alm_trans(Band);
           int FL = MIN(ALM.Lmax(),WP_H.nx()-1);
           if (FL > ALM.Lmax()) FL = ALM.Lmax();
           if(SqrtFilter==false) {
           	for (int l=0; l <= FL; l++) 
               	for (int m=0; m <= l; m++) ALM_Band(l,m) +=  ALM(l,m) ;
   		} else {
           	if(b==NScale) {//Low pass filter
           		for (int l=0; l <= FL; l++) 
               		for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) sqrt(WP_H(l,NbrScale-NScale)));
           	} else {
           		for (int l=0; l <= FL; l++) 
               		for (int m=0; m <= l; m++) ALM_Band(l,m) +=  (ALM(l,m) * (REAL) sqrt(WP_WPFilter(l,b)));
           	}
   		}
           if((b==0)&&(BandLimit==false)) {//FCS Modified: add Info at multipoles > lmax
           	printf("NO BAND LIMIT\n");
           	ALM.alm_rec(Rec);
           	for (int p=0; p < Band.Npix(); p++) Rec[p]= Band[p] - Rec[p];
       	}    
   	}
   	ALM_Band.alm_rec(DataOut);
    if(BandLimit==false) for (int p=0; p < DataOut.Npix(); p++) DataOut[p]+=Rec[p];   	
} 


/****************************************************************************/
void C_PWT::read_scales(char *Name_Imag_In, int &LM,int &NScale, bool order) {
	
	char prefix_inv[512],NameIn[512],suffix_inv[512],TempInv[512], * pos_substr, * pos_substr2;
	int len_substr,offset_str;
	if (strstr(Name_Imag_In, ".fit") != NULL)  strcpy(NameIn, Name_Imag_In);
   	else sprintf(NameIn, "%s.%s", Name_Imag_In, "fits");
   	if (strstr(Name_Imag_In, "_scale0_") == NULL) {
   		printf("Wavelet scales do not follow proper string formatting (should contain _scale0_)\n");
   		exit(-1);	
   	}
	pos_substr=strstr(NameIn,"_scale");
	len_substr=pos_substr-NameIn;
	strncpy(prefix_inv,NameIn,len_substr);
	prefix_inv[len_substr]='\0';
	if ((pos_substr=strstr(NameIn,"_over")) == NULL) {
   		printf("Wavelet scales do not follow proper string formatting (should contain _over)\n");
   		exit(-1);	
   	}
	pos_substr2=strstr(NameIn,".fits");
	offset_str=(pos_substr-NameIn);
	len_substr=pos_substr2-pos_substr;
	memcpy(suffix_inv,&(NameIn[offset_str]),len_substr);
	suffix_inv[len_substr]='\0';
	
	if(order == RING) nest=false;
	else nest=true;
	
	//How many scales? Check according to the name
	ifstream infile;
	int kb=0;
	sprintf(TempInv,"%s_scale%d%s.fits",prefix_inv,kb,suffix_inv);		
	infile.open(TempInv);
	while(infile.is_open()) {
		infile.close();
		kb++;
		sprintf(TempInv,"%s_scale%d%s.fits",prefix_inv,kb,suffix_inv);		
		infile.open(TempInv);
	}
	NScale=kb-1;
	printf("NScale=%d\n",NScale);
		
	PWTTrans = new dblarray[NScale+1];
	WTalloc=true;
	for(int kb=0;kb<=NScale;kb++) {
		sprintf(TempInv,"%s_scale%d%s.fits",prefix_inv,kb,suffix_inv);		
		printf("Process file %s\n",TempInv);
		fits_read_dblarr(TempInv,(PWTTrans)[kb]);
	}
	Nside=sqrt((PWTTrans[0]).nx()/ 12ul);
	if(LM == 0) LM=min(2*Nside,ALM_MAX_L);//default value for Lmax
	wp_alloc(Nside, LM, DEF_MRS_ORDERING);//Ring Format (quicker spherical harmonic transform)

}

/****************************************************************************/
void C_PWT::write_scales(char *Name_Imag_Out, int NScale) {
	char prefix_inv[512],NameOut[512],TempInv[512], * pos_substr;
	int len_substr;
	if (strstr(Name_Imag_Out, ".fit") != NULL)  strcpy(NameOut, Name_Imag_Out);
   	else sprintf(NameOut, "%s.%s", Name_Imag_Out, "fits");
	pos_substr=strstr(NameOut,".fits");
	len_substr=pos_substr-NameOut;
	strncpy(prefix_inv,NameOut,len_substr);
	prefix_inv[len_substr]='\0';
	for(int kb=0;kb<=NScale;kb++) {
		sprintf(TempInv,"%s_scale%d_over%d.fits",prefix_inv,kb,NScale);		
		fits_write_dblarr(TempInv,(PWTTrans)[kb]);
	}
}

/****************************************************************************/
/*
void MRS_GMCA::transform_sources(fltarray &Data, fltarray &DataOut, bool reverse)
{
   for (int c=0; c < Data.ny(); c++)
   {
      for (int f=0; f < 12; f++)
      {
         fltarray Face;
         Ifloat IF;
         Face.alloc(nside(),nside());
         IF.alloc(Face.buffer(), Face.ny(), Face.nx());
               
         for (int x=0; x < nside(); x++)
         for (int y=0; y < nside(); y++) 
                  Face(y,x) = Data( WTTrans.xyf2ind(x,y,f), c);
         if (reverse == false) OWT2D->transform(IF,nscale());
         else OWT2D->recons(IF,nscale());
         for (int x=0; x < nside(); x++)
         for (int y=0; y < nside(); y++) 
                  DataOut( WTTrans.xyf2ind(x,y,f), c) = Face(y,x);
         }
    }
}
 */
/****************************************************************************/
 
 
/*********************************************************************/


 
