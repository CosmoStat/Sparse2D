#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>
#include "powspec_fitsio.h"
#include "PolaHealpixClass.h"
#include "MatrixOper.h"
#include <fstream> 

#define DEF_NITER_ITER_INV 100

template< class PARAM> 
struct is_double { static const bool value = false;}; 
template<> 
struct is_double< double> { static const bool value = true;}; 

enum PSPEC_TYPE{PSPEC_ALL,PSPEC_TT,PSPEC_EE,PSPEC_BB,PSPEC_TE,PSPEC_TB,PSPEC_EB};
enum MLL_TYPE{MLL_TTTT,MLL_EEEE,MLL_EEBB,MLL_TETE,MLL_EBEB,MLL_EB_4BLOCKS};


template<typename T> class MasterPola: MatOper
{

    private:
        PolaHmap<T> Map_TQU;
        PolaHmap<T> NoiseMap_TQU;
        Hmap<T> MaskT;
        Hmap<T> MaskP;
        int Nside_map;    
        Healpix_PowSpec MskMapPowSpec;
        Healpix_PowSpec MskNoisePowSpec;
        arr<double> MaskTT_Powspec; 
        arr<double> MaskPP_Powspec;
        arr<double> MaskTP_Powspec;
        dblarray  MAT_TT_TT, MAT_EE_EE, MAT_EE_BB, MAT_TE_TE, MAT_EB_EB; // Master Coupling matrices 
        dblarray  inv_MAT_TT_TT, inv_MAT_EE_EE, inv_MAT_EE_BB, inv_MAT_TE_TE, inv_MAT_EB_EB, inv_MAT_EB_4Blocks;//inv matrices if SVD decomposition
        double gamma_TT_TT, gamma_EE_EE, gamma_EE_BB,gamma_TE_TE, gamma_EB_EB, gamma_EB_4BLOCKS; //spectral radii useful for iterative inversion
        arr<double> Master_TT_Powspec; 
          arr<double> Master_EE_Powspec;
          arr<double> Master_BB_Powspec;
          arr<double> Master_TE_Powspec;
          arr<double> Master_TB_Powspec;
          arr<double> Master_EB_Powspec;
        int Nmaps;
        int lLmax;
        bool MaskTFlag;
        bool MaskPFlag;
        bool PMaskFlag;
        bool PMapFlag;
        bool PNoiseMapFlag;
        bool Verbose;
        bool ZeroFirst2L;
        bool Flag_Work_EB;
        bool FlagUncorrNoise;
        bool Timer;
    public:
        MasterPola():  Nside_map(0),
            gamma_TT_TT(0), gamma_EE_EE(0), gamma_EE_BB(0), gamma_TE_TE(0), gamma_EB_EB(0), gamma_EB_4BLOCKS(0), Nmaps(0),lLmax(0), MaskTFlag(false), MaskPFlag(false), PMaskFlag(false), PMapFlag(false), PNoiseMapFlag(false), Verbose(false),  ZeroFirst2L(false), Flag_Work_EB(false), FlagUncorrNoise(false),Timer(false) { }

        //Setting maps or useful constants directly
        void set_Map_TQU(PolaHmap<T> cMap_TQU) { 
            Map_TQU= cMap_TQU;
            Nside_map = Map_TQU.get_nside();
            if(Verbose) printf("Nside Map:%ld\n", Nside_map);
            Nmaps=3;
        }
        void set_Lmax(const int clLmax) {lLmax= clLmax;}
        void set_First2L_tozero( bool Flag_Zero) {ZeroFirst2L = Flag_Zero;}
        void set_WORK_EB(bool cEBWork) {
            if(Verbose) if(cEBWork) printf("MASKING TEB MAPS\n");
            Flag_Work_EB = cEBWork;
        }
        void set_Uncorr_Noise(bool cNoiseUncorr) {FlagUncorrNoise = cNoiseUncorr;}
        void set_Verbose(bool cVerbose) {Verbose = cVerbose;}// MatOper::Verbose=(Verbose_ ? True: False); }
        void set_Timer(bool cTimer) {Timer = cTimer;}
        void set_Map_T(const Hmap<T> cMap_T) {Map_TQU.map_T= cMap_T;Nside_map = cMap_T.Nside();if(Verbose) printf("Nside Map:%d\n", Nside_map);Nmaps=1;}
        void set_MaskT(const Hmap<T>  cMask_T) {MaskT = cMask_T; MaskTFlag=true;}
        void set_MaskP(const Hmap<T>  cMask_P) {MaskP = cMask_P; MaskPFlag=true;}
        void set_SpecRad(const double & gamma, const MLL_TYPE selspec) {
            switch (selspec){
                case MLL_TTTT: gamma_TT_TT=gamma;break;
                case MLL_TETE: gamma_TE_TE=gamma;break;
                case MLL_EEEE: gamma_EE_EE=gamma;break;
                case MLL_EEBB: gamma_EE_BB=gamma;break;
                case MLL_EB_4BLOCKS: gamma_EB_4BLOCKS =gamma;break;//same as gamma_EE_EE
                case MLL_EBEB: gamma_EB_EB=gamma;break;//same as gamma_EE_EE
                default:
                    printf("Can only set individual spectral radii from TTTT to EBEB\n");
            }
        } 
        
        //Setting maps from fits files
        void read_PolaHmap(char * Name, bool flag_teb, bool set_ring=false, bool PolaFastALM=true ) {Map_TQU.read(Name, flag_teb,set_ring);
            Nside_map = Map_TQU.map_T.Nside();if(Verbose) printf("Nside Map:%d\n", Nside_map);Nmaps=3;
            if((Flag_Work_EB) && (flag_teb == false))  Map_TQU.swap_tqu_teb(PolaFastALM);
            if((! Flag_Work_EB) && (flag_teb == true))  Map_TQU.swap_tqu_teb(PolaFastALM);
        }
        void read_NoisePolaHmap(char * Name, bool flag_teb, bool set_ring=false, bool PolaFastALM=true) {
            NoiseMap_TQU.read(Name, flag_teb,set_ring);
            if(NoiseMap_TQU.map_T.Nside() != Nside_map) printf("BEWARE MISMATCH IN NSIDE FROM MAPS (%d) AND NOISE (%d)\n",Nside_map,NoiseMap_TQU.map_T.Nside());
            if((Flag_Work_EB) && (flag_teb == false))  NoiseMap_TQU.swap_tqu_teb(PolaFastALM);
            if((! Flag_Work_EB) && (flag_teb == true))  NoiseMap_TQU.swap_tqu_teb(PolaFastALM);
        }
        void read_HmapT(char * Name, bool set_ring=false) {Map_TQU.map_T.read(Name, set_ring);
            Nside_map = Map_TQU.map_T.Nside();if(Verbose) printf("Nside Map:%d\n", Nside_map);Nmaps=1;}
        void read_NoiseHmapT(char * Name, bool set_ring=false) {NoiseMap_TQU.map_T.read(Name, set_ring);
            if(NoiseMap_TQU.map_T.Nside() != Nside_map) printf("BEWARE MISMATCH IN NSIDE FROM MAPS (%d) AND NOISE (%d)\n",Nside_map,NoiseMap_TQU.map_T.Nside());
        }
        void read_MaskT(char * Name, bool set_ring=false) {MaskT.read(Name, set_ring); MaskTFlag=true;}
        void read_MaskP(char * Name, bool set_ring=false) {MaskP.read(Name, set_ring); MaskPFlag=true;}
                
        //Computing TEB Pspec and XPspec
        void get_TEB_Pspec_from_map(bool NormALM=false, bool PolaFastALM=true);
        void get_TEB_Pspec_from_Noisemap(bool NormALM=false, bool PolaFastALM=true);
        void get_Pspec_from_masks(bool NormALM=false, bool PolaFastALM=true);

        //Master Routines
        void get_MASTER_pspec(unsigned short xsubi[3],PSPEC_TYPE selspec=PSPEC_ALL, bool iter=true,int Niter= DEF_NITER_ITER_INV,bool Positivity = true, bool Fast_Radius = true);
        double Check_Matrix_2Norm(dblarray &Matrix,unsigned short xsubi[3], int nit_max=100);
        void make_mll_blocks_c(); //C code adapted from M. Tristram //Computing Coupling Mask matrices
        void make_mll_blocks_c_fast();
        void Mult_cl_matmask_1block(const arr<double> & ClIn,const dblarray & MatMask,arr<double> & ClOut,bool transpose=false);
        void SingleBlock_deconv_iter(const arr<double> & ClData, arr<double> & ClSol, dblarray & MatMask, const int Niter, const double mu=1., bool Positivity=0,bool ZeroFirst2Flag=false);
        void FourSymBlock_deconv_iter(const arr<double> & ClData1,const arr<double> & ClData2, arr<double> & ClSol1, arr<double> & ClSol2, dblarray & MatMask1, dblarray & MatMask2, const int Niter, const double mu=1., bool Positivity=0);

        //Read Intermediary products
        void read_TEB_Pspec_from_file(char *Name ) {
            MskMapPowSpec.read(Name);
            Nmaps=(MskMapPowSpec.num_specs == 6) ? 3:1;
            printf("Nmaps=%d, %d\n",Nmaps,MskMapPowSpec.num_specs);
            if(lLmax >0) lLmax=MIN(MskMapPowSpec.Lmax(),lLmax);
            else lLmax=MskMapPowSpec.Lmax();
            PMapFlag=true;
        }
        void read_NoiseTEB_Pspec_from_file(char *Name ) {
            MskNoisePowSpec.read(Name);
            int cNmaps=(MskNoisePowSpec.num_specs == 6) ? 3:1;
            printf("Nmaps=%d, %d\n", cNmaps, MskNoisePowSpec.num_specs);
            if(lLmax >0) lLmax=MIN(MskNoisePowSpec.Lmax(),lLmax);
            else lLmax=MskNoisePowSpec.Lmax();
            PNoiseMapFlag=true;
        }
        void read_Mask_Pspec_from_file(char *Name ) {
            dblarray spec;
            char* FitsName= fitsname(Name);
            fits_read_dblarr(FitsName, spec );
            free(FitsName);
            int size = spec.nx();
            int num_specs_mask = spec.ny();
            if(lLmax+1 > size) {
                printf("Change lLmax from %d to %d (powerspec and maps not consistent\n",lLmax,size-1);
                lLmax=size-1;
            } else printf("lMax used for Mask PS= %d\n",lLmax);
            MaskTT_Powspec.alloc(lLmax+1);
            for(int l=0; l <= lLmax; l++) MaskTT_Powspec[l]=spec(l,0);
            if(num_specs_mask==3) {
                MaskTP_Powspec.alloc(lLmax+1);
                MaskPP_Powspec.alloc(lLmax+1);
                for(int l=0; l <= lLmax; l++) MaskTP_Powspec[l]=spec(l,1);
                for(int l=0; l <= lLmax; l++) MaskPP_Powspec[l]=spec(l,2);
            } else {
                MaskTP_Powspec=MaskTT_Powspec;
                MaskPP_Powspec= MaskTT_Powspec;
            }
            PMaskFlag=true;
        }
        void read_coupling_mat(char * Name);//Reading Coupling matrices
        void read_spec_radii(char * Name);//Reading spectral radii (double array)
        void read_inv_mat(char * Name);//Reading Inverse matrices

        //Write Intermediary products
        void write_PolaHmap_EB(char * Name, bool PolaFastALM=true) {
            if(Nmaps >0) { 
                PolaHmap<T> Map_Save;
                Map_Save.alloc(Map_TQU.get_nside(), Map_TQU.flag_teb(), Map_TQU.flag_nested());
                Map_Save.set_map_T(Map_TQU.map_T);
                Map_Save.set_map_Q(Map_TQU.map_Q);
                Map_Save.set_map_U(Map_TQU.map_U);
                if(! Map_TQU.flag_teb()) Map_Save.swap_tqu_teb(PolaFastALM);
                Map_Save.write(Name); 
            } else printf("No TQU maps to write\n");} //writing TQU Maps
        void write_HmapT(char * Name) {if(Nmaps >0) Map_TQU.map_T.write(Name); else printf("No T map to write\n");} //write T Maps
        void write_MaskT(char * Name) {if(MaskTFlag) MaskT.write(Name); else printf("No T mask to write\n");} //write Mask T
        void write_MaskP(char * Name) {if(MaskPFlag) MaskP.write(Name); else printf("No P mask to write\n");} //write Mask P
        void write_MaskTEB_Pspec(char * Name) {
            if(PMapFlag) {
                if(Verbose) printf("MasterPola: Masked Map PS [%d]\n",(int) MskMapPowSpec.tt_.size()); 
                if(Verbose) printf("MasterPola: Masked Map PS dim [%d]\n",(int) MskMapPowSpec.num_specs); 
                MskMapPowSpec.write(Name);
            } else printf("No TEB PSpec to write\n");
        } //writing TEB Pspec
        void write_MaskNoiseTEB_Pspec(char * Name) {
            if(PNoiseMapFlag) {
                if(Verbose) printf("MasterPola: Masked Noise Map PS [%d]\n",(int) MskNoisePowSpec.tt_.size()); 
                if(Verbose) printf("MasterPola: Masked Noise Map PS dim [%d]\n",(int) MskNoisePowSpec.num_specs); 
                MskNoisePowSpec.write(Name);
            } else printf("No TEB NOISE PSpec to write\n");
        } //writing TEB Pspec
        void write_Mask_Pspec(char * Name) { //Write Masks T/P Pspec 
            if(PMaskFlag) {
                int size = MaskTT_Powspec.size();
                dblarray spec( size, 3);
                for(int l=0; l < size; l++) {
                    spec( l, 0 ) = MaskTT_Powspec[l];
                    spec( l, 1 ) = MaskTP_Powspec[l];
                    spec( l, 2 ) = MaskPP_Powspec[l];
                }    
                 char* FitsName= fitsname(Name);
                   fits_write_dblarr(FitsName, spec );
                free(FitsName);
            } else printf("No Mask PSpec to write\n");
        }
        void write_coupling_mat(char * Name); //writing coupling matrices (double array)
        void write_spec_radii(char * Name); //writing spectral radii (double array)
        void write_inv_mat(char * Name); //writing inverse matrices (double array)
        void write_master_PSPEC(char *Name); //writing MASTER Pspec (double array), one per hdu
        void write_master_allPSPEC(char *Name); //writing MASTER Pspec (double array), all in one hdu

        //IDL/C++ Wrappers
        //Start with all routines to set data
        void idl_set_map_TQU(T *cMap, long Npix, bool cflag_tqu,bool onlyT, bool cflag_nested=true, bool PolaFastALM=true );
        void idl_set_Noisemap_TQU(T *cNoiseMap, long Npix, bool cflag_tqu,bool onlyT, bool cflag_nested=true, bool PolaFastALM=true );
        void idl_set_MaskT(T *MaskT, long Npix, bool cflag_nested=true, bool set_ring=false);
        void idl_set_MaskP(T *MaskP, long Npix, bool cflag_nested=true, bool set_ring=false);
        void idl_set_spec_TEB(double* PSPEC_MSKMAP, bool onlyT);
        void idl_set_spec_NoiseTEB(double* PSPEC_MSKNOISE, bool onlyT);
        void idl_set_spec_mask(double* PSPEC_MASKS, bool onlyT);
        void idl_set_coupling_mat(double* KMat, bool onlyT);
        void idl_set_spec_radii(double* SpecRad);
        void idl_set_inv_coupling_mat(double* InvMat);
        //routines to get the intermediary and final results
        void idl_get_TEBmap(T *EBMap, bool PolaFastALM=true); 
        void idl_get_spec_masks(T *PSPEC_MASKS, bool onlyT); 
        void idl_get_spec_mskmap(T *PSPEC_MSKMAP);
        void idl_get_spec_msknoisemap(T *PSPEC_MSKNOISE); 
        void idl_get_coupling_mat(double* KMat);
        void idl_get_spec_radii(double* SpecRad);
        void idl_get_inv_coupling_mat(double* invMat);
        void idl_get_master_pspec(double* MasterPspec);
};
 
