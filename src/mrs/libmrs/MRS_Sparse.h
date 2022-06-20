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
**    DESCRIPTION  Class interface to the Healpix package   
**    ----------- 
**                 
******************************************************************************/

#ifndef _SPARSE_H_
#define _SPARSE_H_

#include "HealpixClass.h"
#include "SB_Filter.h"
// #include <sparse2d/GMCA.h>
#include <fstream>

class SparseTrans {
    REAL Buf1;
public:
    Bool Verbose;
    SparseTrans(){}
    ~SparseTrans(){}
    virtual int n_elements() = 0;
    virtual REAL  & operator() (int Num) = 0; 
    virtual void transform(Hmap<REAL> & DataIn) = 0;
    virtual void recons(Hmap<REAL> & DataOut)=0;
};


// class C_OWT: public SparseTrans {
class C_OWT {
    bool nest;
    int Nside;
    type_sb_filter T_FilterBank;
    int NbrScale;
    FilterAnaSynt *FAS;
    SubBandFilter *SBF;
public:
    int nscale() { return NbrScale; }
    int nside() { return Nside; }
    Ortho_2D_WT *OWT2D;
    C_OWT() {NbrScale=0;}
    ~C_OWT(){};
    Hmap<REAL> WTTrans;
    int n_elements(){return WTTrans.Npix();}
    void alloc(int NsideIn, int NScale, bool nested=false);
    REAL * buffer () { return WTTrans.buffer();}
    REAL & operator()(int i) { return WTTrans[i];} 
    void read(char * Name, bool set_ring=true) {}   // read a Healpix map
    void write(char *Name){}  // write a Healpix map 
    void transform(Hmap<REAL> & DataIn);
    void recons(Hmap<REAL> & DataOut);
};


class C_UWT2D_ATROUS { //FCS ADDED for lGMCA inversion
    bool nest;
    int Nside;
    int NbrScale;
    int nb_thr;
public:
   	Hmap<REAL> *WTTrans;
    int AllocMem;
    int nscale() { return NbrScale; }
    int nside() { return Nside; }
    ATROUS_2D_WT *UWT2D;
    C_UWT2D_ATROUS() {NbrScale=0;AllocMem=0;Nside=0;nb_thr=0;nest=true; UWT2D = new ATROUS_2D_WT(); 
}
    int n_elements(int Plane=0){
    	int nelem;
    	if(Plane < NbrScale) nelem=(WTTrans[Plane]).Npix();
    	else nelem=0;
    	return nelem;
    }
    void alloc(int NsideIn, int NScale, bool nested=true, type_border TBorder=I_CONT, int lnb_thr=0);
    inline void alloc(int NsideIn, int NScale, bool nested=true, int lnb_thr=0) { alloc(NsideIn,NScale,nested,I_CONT,lnb_thr);};
    REAL * buffer (int p=0) { return (WTTrans[p]).buffer();}
    REAL & operator()(int i, int p=0) { return (WTTrans[p])[i];} 
    void transform(Hmap<REAL> & DataIn);
    void recons(Hmap<REAL> & DataOut);
    ~C_UWT2D_ATROUS(){delete [] UWT2D;if(AllocMem == 1) delete [] WTTrans;};
};





enum type_alm_filter {F_ALM_SPLINE, F_ALM_MEYER};


//Pyramidal WT class
class C_PWT {
    bool nest;
    int Nside;
    int NbrScale;//total (maximal) number of scales for the given Nside
    fltarray WP_W; // WP_W(0:NbrScale-1,0:1) = lmin and lmax of each band (corser resol = WP_H(*,0)
    fltarray WP_H; // WP_H(0:lmax, 0:NbrScale-1) = filter h (corser resol = WP_H(*,0) 
    fltarray WP_G; // WP_G(0:lmax, 0:NbrScale-1) = filter g (corser resol = WP_G(*,0) 
    fltarray TabPsi;
    fltarray WP_WPFilter;
    intarray NsidePerBand;
    intarray NpixPerBand;
    int Lmax;
    bool WTalloc;
    int ALM_iter;
  public:
	type_alm_filter T_Fil;
    bool All_WP_Band;   // if true use the planck binning powspec band
    int nscale() { return NbrScale; }
    int nside() { return Nside; }
    void set_alm_iter(int ALM_IT) { ALM_iter=ALM_IT; }
    C_PWT() {NbrScale=0;PWTTrans=NULL;T_Fil=F_ALM_MEYER;Nside=0;nest=false;All_WP_Band=false;
              Lmax=0;WTalloc=false;ALM_iter=0;}
    ~C_PWT(){if(WTalloc==true) delete [] PWTTrans;}
    dblarray *PWTTrans;
    void wp_alloc(int NsideIn, int LM, bool nested=false);
    REAL * buffer () { return (PWTTrans[0]).buffer();}
    REAL & operator()(unsigned long i, int b) { return (PWTTrans[b])(i);} 
    void read_scales(char * Name_Imag_In, int &LM,int &NScale, bool order=DEF_MRS_ORDERING);   // read a Healpix map
    void write_scales(char *Name_Imag_Out, int NScale);  // write a Healpix map 
    void transform(Hmap<REAL> & DataIn, bool BandLimit, bool SqrtFilter, int NScale); //FCS Added
    void recons(Hmap<REAL> & DataOut, bool BandLimit, bool SqrtFilter, int NScale);//FCS Added
};

//Undecimated WT class
void mrs_wt_trans(Hdmap & Map, dblarray & TabCoef, fltarray & WP_H,  int Lmax, int NScale, int ALM_iter=0);
   

class C_UWT {
    bool nest;
    bool TightFrame;
    int Nside;
    int NbrScale;
    int NpixPerBand;
public:
    fltarray WP_W; // WP_W(0:NbrScale-1,0:1) = lmin and lmax of each band (corser resol = WP_H(*,0)
    fltarray WP_H; // WP_H(0:lmax, 0:NbrScale-1) = filter h (corser resol = WP_H(*,0) 
    fltarray WP_G; // WP_G(0:lmax, 0:NbrScale-1) = filter g (corser resol = WP_G(*,0) 
    fltarray TabPsi;
    fltarray WP_WPFilter;
    int Lmax;
    int ALM_iter;
    bool Verbose;

    dblarray TabNorm;
    fltarray TabMad;
    type_alm_filter T_Fil;
    bool MeyerWT;
    bool All_WP_Band;   // if true use the planck binning powspec band
    int nscale() { return NbrScale; }
    int nside() { return Nside; }
    void set_alm_iter(int ALM_IT) { ALM_iter=ALM_IT; }
     C_UWT() {NbrScale=0;T_Fil=F_ALM_MEYER;Nside=0;nest=false;All_WP_Band=false;
         Lmax=0;ALM_iter=0;Verbose=false;}
    ~C_UWT(){}
    dblarray WTTrans;
    int n_elements(){return NpixPerBand*NbrScale;}
    void wt_alloc(int NsideIn, int NScale, int LM, bool nested=false, bool Tight_Frame=false);
    void wp_alloc(int NsideIn, int LM, bool nested=false);
    REAL * buffer () { return WTTrans.buffer();}
    REAL & operator()(int i) { return  WTTrans(i % NpixPerBand, i / NpixPerBand);} 
    REAL & operator()(int b, int i) { return WTTrans(i,b);} 
    void read(char * Name, bool set_ring=true) {}   // read a Healpix map
    void write(char *Name){}  // write a Healpix map 
    void transform(Hmap<REAL> & DataIn, bool BandLimit=false);
    void transform(Hmap<REAL> & DataIn, bool BandLimit, bool SqrtFilter, int NScale); //FCS Added
    void recons(Hmap<REAL> & DataOut, bool BandLimit, bool SqrtFilter, int NScale);//FCS Added
    void recons(Hmap<REAL> & DataOut, bool BandLimit=false);
    void set_band(int b, float Value=0.);
    void hard_thresholding(int b, float NSigma, float & SimgaNoise, bool UseMad=false);
    void hard_thresholding(Hmap<REAL> & DataIn, float NSigma, float & SigmaNoise, bool UseMad=false, bool KillLastScale=false, int FirstDetectScale=0);
};


//
void wp_trans(Hdmap & Map, fltarray & TabCoef, fltarray & WP_W, fltarray & WP_H, int NbrWP_Band, int Lmax, bool BandLimit, bool SqrtFilter, int ALM_iter=0);//FCS added BandLimit Option
void wp_trans(Hdmap & Map, dblarray **TabCoef,intarray &NsideBand, fltarray & WP_W, fltarray & WP_H, int NbrWP_Band, int Lmax,bool BandLimit, bool SqrtFilter, int NScale, int ALM_iter=0);//FCS ADDED
// Compute the wavelet packet decomposition using a set of predefined filters

void get_wp_meyer_filter(int nside, fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int  Lmax);
// Compute the large set of filter bands, corresponding to the PLANCK binning of the powspec

void get_planck_wp_meyer_filter(fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int Lmax, int LmaxT);
// Compute the filter bands for the wavelet packet decomposition

void get_wt_bspline_filter(int nside, fltarray &TabH, fltarray &TabG, fltarray &Win, fltarray & WP_WPFilter, int  Lmax, int NbrScale, bool TightFrame=false);

/*
class MRS_GMCA: public C_OWT, public GMCA 
{
public:
    bool Nested;
    MRS_GMCA():GMCA() {Nested=false;};
    void run(fltarray &TabCannels, fltarray & TabSource, fltarray & InvMixingMat);
    void transform_sources(fltarray &DataIn, fltarray &DataOut, bool reverse);
    void transrecons_sources(fltarray &TabVect,fltarray &Recdata) {};
    ~MRS_GMCA() {} ;
};
*/

/*********************************************************************/
/*
template<typename T> class MRS_PWT: public  SparseTrans<T>
{
public:
    inline int nside() {  return Healpix_Map<T>::Nside(); }
    T * buffer ();
    void alloc( int NsideIn, bool nested );
    
    void read(char * Name, bool set_ring=true);   // read a Healpix map
    void write(char *Name);  // write a Healpix map 
    void transform(Hmap<T> & DataIn, Hmap<T> & WTOut, int NbrScale)
    void recons(Hmap<T> & WTin, Hmap<T> & DataOut, int NbrScale)
};

template<typename T> class MRS_UWT: public SparseTrans<T> 
{
public:
    void alloc( int NsideIn, bool nested );
    
    void read(char * Name, bool set_ring=true);   // read a Healpix map
    void write(char *Name);  // write a Healpix map 
    void transform(Hmap<T> & DataIn, Hmap<T> & WTOut, int NbrScale)
    void recons(Hmap<T> & WTin, Hmap<T> & DataOut, int NbrScale)
};

template<typename T> class MRS_CUR: public  SparseTrans<T> 
{
public:
     T * buffer ();
    void alloc( int NsideIn, bool nested );
    
    void read(char * Name, bool set_ring=true);   // read a Healpix map
    void write(char *Name);  // write a Healpix map 
    void transform(Hmap<T> & DataIn, Hmap<T> & WTOut, int NbrScale)
    void recons(Hmap<T> & WTin, Hmap<T> & DataOut, int NbrScale)
};
*/
/*********************************************************************/

#endif
