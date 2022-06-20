/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.-L.  Starck
**
**    Date: 24/06/2020
**
**    File:  WLS_MassMapping.h
**
*******************************************************************************
**
**    DESCRIPTION  Class interface to the Healpix package
**    -----------
**
******************************************************************************/

#ifndef _WLS_MassMapping_
#define _WLS_MassMapping_

// #include "Healpix_PowSpec.h"
#include "HealpixClass.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "compatibility_healpix.h"
#include"MRS_Sparse.h"

#define WL_INFINITE_COV_VALUE 1e9
// ALM_DEF_NITER = 0 in Healpix Class. We may want to
// have a different setting gor WL applications.
#define WL_DEF_ALM_NITER 0

/*********************************************************************/

/* A HEALPix   shear class (2 maps) of a given datatype */

template<typename T> class SP_WL_Field
{
    private:
        bool nested;
        int Nside;
        bool FlagEB;

    public:
        bool Verbose;
        Hmap<T> map_G1;
        Hmap<T> map_G2;
        
        /*! Constructs an unallocated map. */
        SP_WL_Field () {Nside=0;FlagEB=false;nested=false;Verbose=True;}
                
        int get_nside() { return Nside; }
        inline bool flag_eb() { return FlagEB; }
        inline bool flag_nested() { return nested; }
        long Npix() {return (long) map_G1.Npix(); }

        void alloc( int NsideIn, bool flag_eb=false, bool flag_nested=false);
        void alloc(Hmap<T> &MapG1, Hmap<T> &MapG2, bool flag_eb=false, bool flag_nested=false);
        void read( char * NameG1, char * NameG2, bool flag_eb=false, bool set_ring=false );   // read a Healpix map
        void write( char *NameG1, char * NameG2);  // write a Healpix map
        inline Hmap<T> get_map_G1() { return map_G1; }
        inline Hmap<T> get_map_G2() { return map_G2; }
        Hmap<T>* get_map(int index) {
            switch (index){
                case 1: return &map_G1;break;
                case 2: return &map_G2;break;
                default:
                    printf("ONLY 2 MAPS AVAILABLE (G1,G2)\n");
                    exit(EXIT_FAILURE);
            }
        }
        
        REAL & operator()(int PolType, int PixNum)
               { if(PolType == 1) {return map_G1[PixNum];}
             else if (PolType == 2) return map_G2[PixNum];
               }

        int set_map_G1( Hmap<T> & map );
        int set_map_G2( Hmap<T> & map );

        const SP_WL_Field<T>  & operator = (SP_WL_Field<T> & S_Map);
         const SP_WL_Field<T> & operator += (const SP_WL_Field<T> & S_Map);
         const SP_WL_Field<T> & operator *= (const SP_WL_Field<T> & S_Map);
         const SP_WL_Field<T> & operator -= (const SP_WL_Field<T> & S_Map);
         const SP_WL_Field<T> & operator /= (const SP_WL_Field<T> & S_Map);
         
         void swap_gamma_eb( bool fast );
         void swap_nested_ring();
         
         // Deletes the old map and creates a new map  with a given nside and the ordering scheme.
        void shear_set_nside( int nside, bool flag_eb=false, bool flag_nested=false);
        //Set to val all pixels
        void shear_fill(double value );
        
        /*! Imports the map \a orig into the current map, adjusting the ordering scheme and the map resolution if necessary.
        When downgrading, \a pessimistic determines whether or not pixels are set to \a Healpix_undef when not all of the corresponding high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
        void import( SP_WL_Field<T> &S_Map, bool pessimistic=false );
        int import_via_alm( SP_WL_Field<T> &S_Map_in, bool fast=true );

         void info();  // print statistical information
};

#define WLS_Field SP_WL_Field<REAL>      // double healpix map

// ===========================
 
// ===========================

class ShearAlm
{
    private:
        CAlmR ALM_E;        //AlmR ALM_E; //FCS MOD
        CAlmR ALM_B;        //AlmR ALM_B; //FCS MOD
        
        int ShearNside;
        int Shear_Lmax;
        bool ShearFastALM;
        int niter;
        
        bool Norm;  // if True, normalize the Alm coefficients such a Gaussian noise noise with variance S^2, produce Alm with variance 1/2 S^2
        double NormVal; // Alm Nomalization value
                
        arr<double> weight_T; // internal variable
    public:
        dblarray WienerFilter_E;
        dblarray WienerFilter_B;
        ShearAlm(){}
        inline int get_nside(){ return ShearNside; }
        inline int get_lmax(){ return Shear_Lmax; }
        inline bool flag_FastALM(){ return ShearFastALM; }
        inline int get_niter(){ return niter; }
        inline double get_normval(){ return NormVal; }
        inline bool flag_NormALM(){ return Norm; }
        
        inline CAlmR* get_alm(int index) {
            switch (index){
                case 1: return &ALM_E;break;
                case 2: return &ALM_B;break;
                default:
                    printf("ONLY 2 ALM COMPONENTS AVAILABLE (FlagEB)\n");
                    exit(EXIT_FAILURE);
            }
        }
        dblarray get_WienerFilter_E(){ return WienerFilter_E; }
        dblarray get_WienerFilter_B(){ return WienerFilter_B; }
        
        void set_niter( int nb_iter ){ niter = nb_iter; }
        void set_flag_NormALM( bool flag_norm );
        
        int set_Alm_E( AlmR & ALM );
        int set_Alm_B( AlmR & ALM );
        
        void alloc(int Nside, int u_lmax, bool Fast);
    
        void read( char *Name, int Nside, int u_lmax=0, bool Fast=DEF_ALM_FAST  ); // read the Alm from a fits file
        void read_array( char *Name );
        void write( char *Name );  // write the Alm to fits file
        void write_array( char *Name );    // write Alm "array" format in files "ALM_T_Name.xxx" "ALM_E_Name.xxx" "ALM_B_Name.xxx"
        void write_array2( char *Name );
        
        void shear_alm_trans( WLS_Field & SMap, bool UseSpin2=true); // Alm transformation of a Healpix map
        void shear_alm_rec( WLS_Field & SMap, bool UseSpin2=true, bool KeepAlm=false, int RecNside=0); // Alm inverse transform
        
        void shear_alm2powspec( PowSpec & specE, PowSpec & specB);  // compute the power spectrum from the Alm
        
        void convol(fltarray &Filter_E, fltarray &Filter_B ); // convol the Alm's with a given filter for each component
        void convol( fltarray &Filter ); // convol the Alm's with the same filter for all components
        void convol( float Fwhm );
        
        void wiener( PowSpec & ps_noise, PowSpec & ps_signal ); // Apply a wiener filtering to the Alm, knowing the signal and noise power spectrum
        void set_wiener_filter(PowSpec & ps_noise, PowSpec & ps_signal);
        void alm_mult_wiener();
        void mult_wiener(WLS_Field & SMap, bool UseSpin2Trans=true, bool UseSpin2Rec=true);
        void wiener(WLS_Field & SMap, PowSpec & ps_noise, PowSpec & ps_signal, bool UseSpin2=true);

        double max_absalm_E();
        double max_absalm_B();
        double max_absalm_FlagEB();

        int hard_threshold(float lambda_e, float lambda_b, int & MaxNonZeroL_E, int & MaxNonZeroL_B );//threshold coeff lower thant T and return the number of non zero coeff after thresholding
        int soft_threshold(float lambda_e, float lambda_b, int & MaxNonZeroL_E, int & MaxNonZeroL_B );
        
        void info();
};


// ===========================

enum wl_type_weight {NO_WEIGHT, SIGMA_WEIGHT, VAR_WEIGH};
                    
 class WLS_MassMapping
 {
  private:
    int Nside;
    int Npix;
    int NbrScale;
    Hdmap WeightGamma_Wiener;
    Hdmap WeightGamma_Sparse;
    Hdmap Mask;
    C_UWT WT;             // Class for wavelet decomposition
    Hdmap CovMat;    // CovMat is the diagonal cov matrix on each indivual shear component G1 and G2,
                     // and NOT the diagonal cov matrix on the complex shear field.
                     // covmat(complex shear field) = 2 * covmat(G1) =  2 * covmat(G2)
    
     void get_residual_gamma(WLS_Field & Kappa, WLS_Field & Resi_Gamma, wl_type_weight TypeWeight=NO_WEIGHT);
    void get_residual_eb(WLS_Field & Kappa, WLS_Field & Resi_EB, wl_type_weight TypeWeight=NO_WEIGHT);
    WLS_Field FieldTemp;
    double MinCovMat;
     void mult_sparse(WLS_Field & KappaSparse, bool NoSparseBMode);
     void get_active_coef(WLS_Field & Kappa, float Nsigma, bool KillLastScale=false, bool OnlyPos=true, bool NoSparseMode=true, int FirstDetectScale=0);
  
 public:
     bool WT_BandLimit;
     ShearAlm CAlm;
     dblarray TabActivCoefE;
     dblarray TabActivCoefB;
     bool Verbose;
     WLS_Field GammaData; // Shear field
     WLS_MassMapping(){WT_BandLimit=true;Nside=0;NbrScale=0;Npix=0;Verbose=false;MinCovMat = WL_INFINITE_COV_VALUE;}
     void alloc(Hdmap &G1, Hdmap &G2, Hdmap &CovMatrix, Hdmap &MaskData, int NbrScale=0, int lmax=0);
     void alloc(Hdmap &G1, Hdmap &G2, Hdmap &CovMatrix, int NbrScale=0, int lmax=0);
     void alloc(Hdmap &G1, Hdmap &G2, float SigmaNoise, Hdmap &MaskData, int NbrScale=0, int lmax=0);
     void alloc(Hdmap &G1, Hdmap &G2, float SigmaNoise, int NbrScale=0, int lmax=0);
     void set_niter( int nb_iter ){ CAlm.set_niter(nb_iter); }

     void gamma2eb(WLS_Field & Gamma, WLS_Field & Kappa, bool KeepAlm=false);
     void eb2gamma(WLS_Field & Kappa, WLS_Field & Gamma, bool KeepAlm=false);
     void sp_kaiser_squires(WLS_Field & Gamma, WLS_Field & Kappa, float Fwhm=0);
     void sp_kaiser_squires(WLS_Field & Kappa, float Fwhm=0)
     {
         sp_kaiser_squires(GammaData, Kappa, Fwhm);
     }
     void wiener(WLS_Field & SMap, PowSpec & ps_noise, PowSpec & ps_signal, bool UseSpin2Trans=true, bool UseSpin2Rec=true);

     void iter_wiener(WLS_Field & Gamma, PowSpec & ps_signal, WLS_Field & Kappa, int NiterWiener=10);
     
     void sparse_reconstruction(WLS_Field & Gamma, WLS_Field & Kappa, float NSigma, int NiterSparse=10, int FirstDetectScale=0);

     void mcalens(WLS_Field & Gamma, PowSpec & ps_signal, WLS_Field & Kappa, WLS_Field & KappaSparse, float NSigma, int NiterSparse=10, bool SparsePositivityConstraint=true, int FirstDetectScale=0);
};

// ===========================

#endif

