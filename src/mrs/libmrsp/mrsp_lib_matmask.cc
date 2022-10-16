#include "mrsp_lib_matmask.h"
extern "C"  {
//=========================================================================
//=========================================================================
//=========================================================================

//=========================================================================
// Routines for Mll computation
//
// input:
//    spectra from the masks applied in temperature and polarization from l=0 to lmax
//    well = [well_TT,well_TP,well_PP]
//
// output:
//    mll matrices in RowMajorOrder, 5 matrices of size (lmax+1,lmax+1)
//    mll = [mll_TT_TT, mll_EE_EE, mll_EE_BB, mll_TE_TE, mll_EB_EB
//
// dependancy:
//    Need FORTRAN routine wig3j_f.f to compute 3j-wigner
//
// comment:
//    loop on l3 should be 2times the lmax of the Mll matrices
//    in practice if you want to compute Mlls up to 3*nside then you can NOT
//    get Well up to 2*(3*nside) due to HEALPix pixel sizes
//    At first order, you can neglect the points above 3*nside in the l3 loop as
//    the wigner for large l3 get smaller and smaller...
//
// M. Tristram    - sep 2004 -  Polarized version of make_mll
//                - dec 2004 -  Add different mask for temperature and polarization
//                              clean memory in wig3j
//=========================================================================

#ifdef _OPENMP
     extern int inner_loop_threads;
      extern int outer_loop_threads;
#endif

extern void wig3j_( double *L2, double *L3, double *M2, double *M3,
         double *L1MIN, double *L1MAX, double *THRCOF, int *NDIM, int *IER);
extern void wig3j_c( long l2, long l3, long m2, long m3, double *wigner);
//void make_mll_pol( long lmax, double *well, double *mll);
}

//****************************************************************************************************************************//
template <class T>
double  MasterPola <T>::Check_Matrix_2Norm(dblarray &Matrix,unsigned short xsubi[3],int nit_max) {
//Power Method
//see Estimating the matrix p-norm, NJ Higham, Numer. Math. 62:539-555, 1992

    long Nx=Matrix.axis(1), Ny=Matrix.axis(2);
    double init_vector[Nx],vec_fwd[Ny],vec_bwd[Nx],vec_bwd2[Nx];
    long kl,kl1,kl2, offset;
    double gamma=-1.;
    double Nrm2In=0.,Nrm2fwd,Nrm2bwd,Totalbwd;
    double *buffer_mat=Matrix.buffer();

    for(kl=0;kl<Nx;kl++) {
        init_vector[kl]=erand48(xsubi);
        Nrm2In+=init_vector[kl]*init_vector[kl];
    }

    Nrm2In =sqrt(Nrm2In);
    for(kl=0;kl<Nx;kl++) init_vector[kl] /= Nrm2In;

    for(int kit=0;kit < nit_max;kit++) {
        Nrm2fwd=0.;
        #ifdef _OPENMP
            #pragma omp parallel for default(none) reduction(+: Nrm2fwd)  shared(vec_fwd,init_vector, buffer_mat,Nx,Ny) private(kl1,kl2,offset)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
        #endif
        for(kl1=0;kl1<Ny;kl1++) {
            vec_fwd[kl1]=0.;
            offset=kl1* Nx;
            for(kl2=0;kl2<Nx;kl2++) vec_fwd[kl1] += buffer_mat[offset+kl2] * init_vector[kl2];
            Nrm2fwd+= vec_fwd[kl1]* vec_fwd[kl1];
        }
        Nrm2fwd =sqrt(Nrm2fwd);
        for(kl1=0;kl1<Ny;kl1++) vec_fwd[kl1] /= Nrm2fwd;

        for(kl1=0;kl1<Nx;kl1++) vec_bwd[kl1]=0.;

        #ifdef _OPENMP
            #pragma omp parallel for default(none) shared(vec_bwd, vec_fwd, buffer_mat,Nx,Ny) private(kl1,kl2)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
        #endif
        for(kl1=0;kl1<Nx;kl1++) {
            for(kl2=0;kl2<Ny;kl2++) vec_bwd[kl1]+= buffer_mat[kl2* Nx+kl1]* vec_fwd[kl2];
        }

        Nrm2bwd =0.;
        for(kl1=0;kl1<Nx;kl1++) Nrm2bwd += vec_bwd[kl1]* vec_bwd[kl1];

        Nrm2bwd =sqrt(Nrm2bwd);
        Totalbwd=0.;
        for(kl1=0;kl1<Nx;kl1++) {
            vec_bwd2[kl1] = vec_bwd[kl1]* init_vector[kl1];
            Totalbwd+=vec_bwd2[kl1];
        }
        if(Nrm2bwd < Totalbwd)     break;
        for(kl1=0;kl1<Nx;kl1++) init_vector[kl1] = vec_bwd[kl1] / Nrm2bwd;
        if(this->Verbose) printf("Spectral Radius Estimate it[%d]=%3.10g\n",kit, Nrm2fwd);
    }
    if(gamma < 0.) if(this->Verbose) printf("Spectral Radius computing as not converged in %d iterations\n",nit_max);
    gamma = Nrm2fwd;
    return gamma;
}

//****************************************************************************************************************************//
//****************************************************************************************************************************//
// PSPEC POLA ROUTINES
//****************************************************************************************************************************//
//****************************************************************************************************************************//

template <class T>
void MasterPola <T>::get_Pspec_from_masks(bool NormALM, bool PolaFastALM) {

    struct timeval start,end,diff;
    if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
    CAlmR AlmMskT, AlmMskP;
    if( this->lLmax <= 0 ) this->lLmax = 3 * Nside_map;
    if(this->Verbose) printf("MasterPola: set Lmax\n");
    if((!this->MaskTFlag)&&(!this->MaskPFlag)) {
        printf("Should first specify either MaskT or MaskP\n");
        exit(EXIT_FAILURE);
    }

    if(this->Verbose) printf("MasterPola: Compute Mask Power Spectra\n");
    this->MaskTT_Powspec.alloc(this->lLmax+1);
    this->MaskPP_Powspec.alloc(this->lLmax+1);
    this->MaskTP_Powspec.alloc(this->lLmax+1);

    if(this->MaskTFlag) {
          long NpixMask=this->MaskT.Npix();
        int NsideMask=this->MaskT.Nside();
        AlmMskT.alloc(NsideMask,this->lLmax, PolaFastALM);
        Hdmap mask_dbl;
        int limit;
        if(this->Verbose) printf("MasterPola: Copy MaskT Array %ld\n", NpixMask);
        mask_dbl.alloc(NsideMask,(this->MaskT.Scheme()==NEST) ? true : false);
        for(long kl=0;kl<NpixMask;kl++) mask_dbl[kl]=(double) this->MaskT[kl];
        if(this->Verbose) printf("MasterPola: AlmX for MaskT\n");
        AlmMskT.alm_trans(mask_dbl);
        if(this->Verbose) printf("MasterPola: PS for MaskT\n");
        for(int l=0; l <= this->lLmax; ++l) {
            this->MaskTT_Powspec[l] = norm( AlmMskT(l,0));
            limit = min( l,AlmMskT.Mmax() );
               for(int m=1; m <= limit; ++m) this->MaskTT_Powspec[l] += 2*norm( AlmMskT(l,m) );
               this->MaskTT_Powspec[l]/=(2.*l+1.);
        }
    }
    if(this->MaskPFlag) {
          long NpixMask=this->MaskT.Npix();
        int NsideMask=this->MaskT.Nside();
        AlmMskP.alloc(this->MaskP.Nside(),this->lLmax, PolaFastALM);
        Hdmap mask_dbl;
        int limit;
        if(this->Verbose) printf("MasterPola: Copy MaskP Array\n");
        mask_dbl.alloc(NsideMask,(this->MaskP.Scheme()==NEST) ? true : false);
        for(long kl=0;kl<NpixMask;kl++) mask_dbl[kl]=(double) this->MaskP[kl];
        if(this->Verbose) printf("MasterPola: AlmX for MaskP\n");
        AlmMskP.alm_trans(mask_dbl);
        if(this->Verbose) printf("MasterPola: PS for MaskP\n");
        for(int l=0; l <= this->lLmax; ++l) {
            this->MaskPP_Powspec[l] = norm( AlmMskP(l,0));
            limit = min( l,AlmMskP.Mmax() );
               for(int m=1; m <= limit; ++m) this->MaskPP_Powspec[l] += 2*norm( AlmMskP(l,m) );
               this->MaskPP_Powspec[l]/=(2.*l+1.);
        }
    }
    if((this->MaskTFlag)&&(this->MaskPFlag))    {
        if(this->Verbose) printf("MasterPola: Check MaskP/MaskT compatibility\n");
        planck_assert(AlmMskT.conformable(AlmMskP), "get_Pspec_from_masks: a_lms are not conformable" );
        int limit;
        for(int l=0; l <= this->lLmax; ++l) {
            this->MaskTP_Powspec[l] = ( AlmMskT(l,0)*conj( AlmMskP(l,0) ) ).real();
            limit = min( l,AlmMskT.Mmax() );
               for(int m=1; m <= limit; ++m) this->MaskTP_Powspec[l] += 2*( AlmMskT(l,m)*conj( AlmMskP(l,m) ) ).real();
               this->MaskTP_Powspec[l]/=(2.*l+1.);
        }
    } else {
        if(this->MaskTFlag) {
            this->MaskTP_Powspec=this->MaskTT_Powspec;
            this->MaskPP_Powspec=this->MaskTT_Powspec;
        }
        else {
            this->MaskTP_Powspec=this->MaskPP_Powspec;
            this->MaskTT_Powspec=this->MaskPP_Powspec;
        }
    }
    PMaskFlag=true;

    if(this->Timer==True) {
        gettimeofday(&end, (struct timezone *) NULL);
        timersub(&end,&start,&diff);
           printf("Time to mask power spectra =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
    }


}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::get_TEB_Pspec_from_map(bool NormALM, bool PolaFastALM) {
    if( this->lLmax <= 0 ) this->lLmax = 3 * Nside_map;


    if(Nmaps ==3) {
          long NpixMap= this->Map_TQU.Npix();
        int NsideMap=this->Map_TQU.get_nside();
         struct timeval start,end,diff;
         if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);

        if(this->Verbose) printf("PS for all 3 maps\n");
        if((! this->MaskTFlag)&&(! this->MaskPFlag)) {
            printf("Should specify either MaskT or MaskP\n");
            exit(EXIT_FAILURE);
        }
        if(( this->MaskTFlag)&&( this->MaskT.Nside()!= this->Nside_map)) {
            printf("MaskT should have same size as Maps (Ns=%d vs %d)\n", this->MaskT.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
        if(( this->MaskPFlag)&&( this->MaskP.Nside()!= this->Nside_map)) {
            printf("MaskP should have same size as Maps (Ns=%d vs %d)\n", this->MaskP.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
        this->MskMapPowSpec.num_specs=6;
        PolaHmap<double> Map_TQU_msk;
        Map_TQU_msk.alloc(NsideMap ,this->Map_TQU.flag_teb(), this->Map_TQU.flag_nested());
        if( this->MaskTFlag) for(int kl=0;kl< NpixMap;kl++) Map_TQU_msk.map_T[kl]=(double)  this->Map_TQU.map_T[kl]* this->MaskT[kl];
        else for(long kl=0;kl<NpixMap;kl++) Map_TQU_msk.map_T[kl]=(double) this->Map_TQU.map_T[kl]*this->MaskP[kl];
        if( this->MaskPFlag) {
            for(long kl=0;kl<NpixMap;kl++) Map_TQU_msk.map_Q[kl]=(double) this->Map_TQU.map_Q[kl]* this->MaskP[kl];
            for(long kl=0;kl<NpixMap;kl++) Map_TQU_msk.map_U[kl]=(double) this->Map_TQU.map_U[kl]* this->MaskP[kl];
        } else {
            for(long kl=0;kl<NpixMap;kl++) Map_TQU_msk.map_Q[kl]=(double) this->Map_TQU.map_Q[kl]* this->MaskT[kl];
            for(long kl=0;kl<NpixMap;kl++) Map_TQU_msk.map_U[kl]=(double) this->Map_TQU.map_U[kl]* this->MaskT[kl];
        }
        PolaAlmR AlmX;
        AlmX.alloc( NsideMap,this->lLmax, PolaFastALM );
        AlmX.pola_alm_trans( Map_TQU_msk);
          AlmX.pola_alm2powspec_all(  this->MskMapPowSpec );

        if(this->Timer==True) {
            gettimeofday(&end, (struct timezone *) NULL);
            timersub(&end,&start,&diff);
               printf("Time to get masked map power spectra  =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
        }

    } else { //only TT
          long NpixMap=this->Map_TQU.map_T.Npix();
        int NsideMap=this->Map_TQU.map_T.Nside();

         struct timeval start,end,diff;
         if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);

        if(this->Verbose) printf("PS for map T [lLmax=%d, NpixMap=%ld, NsideMap=%d]\n",this->lLmax,NpixMap,NsideMap);
        if(! this->MaskTFlag) {
            printf("Should specify MaskT\n");
            exit(EXIT_FAILURE);
        }
        if(( this->MaskTFlag)&&( this->MaskT.Nside()!= this->Nside_map)) {
            printf("MaskT should have same size as Maps (Ns=%d vs %d)\n", this->MaskT.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
            this->MskMapPowSpec.num_specs=1;
        PowSpec powspec_TT;
         CAlmR AlmX;
        AlmX.alloc( NsideMap,this->lLmax, PolaFastALM);
        Hdmap map_T_dbl;
        map_T_dbl.alloc(NsideMap,(this->Map_TQU.map_T.Scheme()==NEST) ? true : false);
        if(this->Verbose) printf("MasterPola: Copy Map T Array [%ld]\n", NpixMap);
        for(long kl=0;kl< NpixMap;kl++) map_T_dbl[kl]=(double) this->Map_TQU.map_T[kl] * this->MaskT[kl];
        if(this->Verbose) printf("MasterPola: AlmX for Map T\n");
        AlmX.alm_trans(map_T_dbl);
        if(this->Verbose) printf("AlmX.lmax=%d\n",AlmX.Lmax());
        AlmX.alm2powspec(powspec_TT);
        if(this->Verbose) printf("MasterPola: PS for Map T [%d]\n",powspec_TT.Lmax()+1);
        this->MskMapPowSpec.tt_.alloc(powspec_TT.Lmax()+1);
        for(long kl=0;kl<=powspec_TT.Lmax();kl++) this->MskMapPowSpec.tt_[kl]=powspec_TT.tt(kl);

        if(this->Timer==True) {
            gettimeofday(&end, (struct timezone *) NULL);
            timersub(&end,&start,&diff);
               printf("Time to get masked map TT power spectrum  =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
        }
    }
    PMapFlag=true;
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::get_TEB_Pspec_from_Noisemap(bool NormALM, bool PolaFastALM) {
    if( this->lLmax <= 0 ) this->lLmax = 3 * Nside_map;

    if(Nmaps ==3) {
          long NpixNoiseMap=this->NoiseMap_TQU.Npix();
        int NsideNoiseMap=this->NoiseMap_TQU.get_nside();
         struct timeval start,end,diff;
         if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);

        if(this->Verbose) printf("PS for all 3 maps\n");
        if((! this->MaskTFlag)&&(! this->MaskPFlag)) {
            printf("Should specify either MaskT or MaskP\n");
            exit(EXIT_FAILURE);
        }
        if(( this->MaskTFlag)&&( this->MaskT.Nside()!= this->Nside_map)) {
            printf("MaskT should have same size as Maps (Ns=%d vs %d)\n", this->MaskT.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
        if(( this->MaskPFlag)&&( this->MaskP.Nside()!= this->Nside_map)) {
            printf("MaskP should have same size as Maps (Ns=%d vs %d)\n", this->MaskP.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
        this->MskNoisePowSpec.num_specs=6;
        PolaHmap<double> NoiseMap_TQU_msk;
        NoiseMap_TQU_msk.alloc(NsideNoiseMap ,this->NoiseMap_TQU.flag_teb(), this->NoiseMap_TQU.flag_nested());
        if( this->MaskPFlag) for(int kl=0;kl< NpixNoiseMap;kl++) NoiseMap_TQU_msk.map_T[kl]=(double)  this->NoiseMap_TQU.map_T[kl]* this->MaskP[kl];
        else for(long kl=0;kl<NpixNoiseMap;kl++)  NoiseMap_TQU_msk.map_T[kl]=(double) this->NoiseMap_TQU.map_T[kl]*this->MaskT[kl];
        if( this->MaskPFlag) {
            for(long kl=0;kl<NpixNoiseMap;kl++) NoiseMap_TQU_msk.map_Q[kl]=(double) this->NoiseMap_TQU.map_Q[kl]* this->MaskP[kl];
            for(long kl=0;kl<NpixNoiseMap;kl++) NoiseMap_TQU_msk.map_U[kl]=(double) this->NoiseMap_TQU.map_U[kl]* this->MaskP[kl];
        } else {
            for(long kl=0;kl<NpixNoiseMap;kl++) NoiseMap_TQU_msk.map_Q[kl]=(double) this->NoiseMap_TQU.map_Q[kl]* this->MaskT[kl];
            for(long kl=0;kl<NpixNoiseMap;kl++) NoiseMap_TQU_msk.map_U[kl]=(double) this->NoiseMap_TQU.map_U[kl]* this->MaskT[kl];
        }
        PolaAlmR AlmX;
        AlmX.alloc(NsideNoiseMap,this->lLmax, PolaFastALM );
        AlmX.pola_alm_trans( NoiseMap_TQU_msk);
          AlmX.pola_alm2powspec_all(  this->MskNoisePowSpec );

        if(this->Timer==True) {
            gettimeofday(&end, (struct timezone *) NULL);
            timersub(&end,&start,&diff);
               printf("Time to get masked noise map power spectra  =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
        }

    } else { //only TT
        struct timeval start,end,diff;
        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
          long NpixNoiseMap=this->NoiseMap_TQU.map_T.Npix();
        int NsideNoiseMap=this->NoiseMap_TQU.map_T.Nside();

        if(this->Verbose) printf("PS for map T [lLmax=%d]\n",this->lLmax);
        if(! this->MaskTFlag) {
            printf("Should specify MaskT\n");
            exit(EXIT_FAILURE);
        }
        if(( this->MaskTFlag)&&( this->MaskT.Nside()!= this->Nside_map)) {
            printf("MaskT should have same size as Maps (Ns=%d vs %d)\n", this->MaskT.Nside(), this->Nside_map);
            exit(EXIT_FAILURE);
        }
        this->MskMapPowSpec.num_specs=1;
        PowSpec powspec_TT;
         CAlmR AlmX;
        AlmX.alloc(NsideNoiseMap,this->lLmax, PolaFastALM);
        Hdmap map_T_dbl;
        map_T_dbl.alloc(NsideNoiseMap,(this->NoiseMap_TQU.map_T.Scheme()==NEST) ? true : false);
        if(this->Verbose) printf("MasterPola: Copy Noise Map T Array [%ld]\n", NpixNoiseMap);
        for(long kl=0;kl<NpixNoiseMap;kl++) map_T_dbl[kl]=(double) this->NoiseMap_TQU.map_T[kl] * this->MaskT[kl];
        if(this->Verbose) printf("MasterPola: AlmX for Map T\n");
        AlmX.alm_trans(map_T_dbl);
        if(this->Verbose) printf("AlmX.lmax=%d\n",AlmX.Lmax());
        AlmX.alm2powspec(powspec_TT);
        if(this->Verbose) printf("MasterPola: PS for Map T [%d]\n",powspec_TT.Lmax()+1);
        this->MskNoisePowSpec.tt_.alloc(powspec_TT.Lmax()+1);
        for(long kl=0;kl<=powspec_TT.Lmax();kl++) this->MskNoisePowSpec.tt_[kl]=powspec_TT.tt(kl);

        if(this->Timer==True) {
            gettimeofday(&end, (struct timezone *) NULL);
            timersub(&end,&start,&diff);
               printf("Time to get masked noise map TT power spectrum  =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
        }
    }
    PNoiseMapFlag=true;
}

//****************************************************************************************************************************//
//****************************************************************************************************************************//
// MASTER POLA ROUTINES
//****************************************************************************************************************************//
//****************************************************************************************************************************//

template <class T>
void MasterPola <T>::make_mll_blocks_c( ) {
     double sum_TT, sum_TE, sum_EE_EE, sum_EE_BB, sum_EB;
     double *wigner0, *wigner2;
     long ndim, maxl3;

     struct timeval start,end,diff;
     if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);

     if(this->Nmaps < 1) {
        printf("There should be at least one map at this stage (%d)\n",this->Nmaps);
        exit(EXIT_FAILURE);
     }

    this->MAT_TT_TT.alloc(this->lLmax+1,this->lLmax+1);
    if(this->Nmaps == 3) {
        this->MAT_EE_EE.alloc(this->lLmax+1,this->lLmax+1);
        this->MAT_EE_BB.alloc(this->lLmax+1,this->lLmax+1);
        this->MAT_TE_TE.alloc(this->lLmax+1,this->lLmax+1);
        this->MAT_EB_EB.alloc(this->lLmax+1,this->lLmax+1);
    }

     /* Allocation set out of the loop */
    wigner0 = (double *) malloc( (2*this->lLmax+1) * sizeof(double));
    if(this->Verbose) printf("wigner0 allocated\n");
    if(this->Nmaps == 3) {
        if(this->Verbose) printf("Wigner 2 allocated\n");
        wigner2 = (double *) malloc( (2*this->lLmax+1) * sizeof(double));
    }

     /* loop over the matrice elements */
    for( long l1=0; l1<=(long) this->lLmax; l1++) {
        //if (this->Verbose) printf("MasterPola: computing MLL[%ld] (out of %d) \n",l1,this->lLmax);
        for( long l2=0; l2<= (long) this->lLmax; l2++) {
            /* alloc wigners */
            ndim = l2+l1-abs(l1-l2)+1;
            if(ndim > (2*this->lLmax+1)) printf("BEWARE OF NDIM !!!! [%ld]\n",ndim);

            /* compute wigners */
            wig3j_c( l1, l2, (long) 0,  (long) 0, wigner0);
            if(this->Nmaps == 3) wig3j_c( l1, l2,  (long)-2,  (long)2, wigner2);

            /* loop on l3 */
            maxl3 = (l1+l2 < this->lLmax) ? l1+l2 : this->lLmax;
            if(this->Nmaps !=3) {
                /* initialization */
                sum_TT    = 0.;
                for( long l3=abs(l1-l2); l3<=maxl3; l3++)
                if( (l1+l2+l3)%2 == 0) sum_TT    += MaskTT_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner0[l3];
                this->MAT_TT_TT(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_TT;
            } else {
                /* initialization */
                sum_TT    = 0.;
                sum_TE    = 0.;
                sum_EE_EE = 0.;
                sum_EE_BB = 0.;
                sum_EB    = 0.;

                for( long l3=abs(l1-l2); l3<=maxl3; l3++) {
                    if( (l1+l2+l3)%2 == 0) sum_TT    += MaskTT_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner0[l3];
                    if( (l1+l2+l3)%2 == 0) sum_TE    += MaskTP_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner2[l3];
                    if( (l1+l2+l3)%2 == 0) sum_EE_EE += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                    if( (l1+l2+l3)%2 != 0) sum_EE_BB += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                    sum_EB += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                }
                this->MAT_TT_TT(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_TT; //BEWARE: M_ll' = M(l',l) (first index=column)
                this->MAT_EE_EE(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_EE_EE;
                this->MAT_EE_BB(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_EE_BB;
                this->MAT_TE_TE(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_TE;
                this->MAT_EB_EB(l2,l1) = (2.*double(l2)+1.)/(4.*M_PI) * sum_EB;
            }
            for( long l3=0; l3< ndim; l3++) wigner0[l3]=0.;
            if(this->Nmaps == 3) for( long l3=0; l3< ndim; l3++) wigner2[l3]=0.;
        } //end loop l2
    } //end loop l1
    free(wigner0);
    if(this->Nmaps == 3) free(wigner2);

    if(this->Timer==True) {
        gettimeofday(&end, (struct timezone *) NULL);
        timersub(&end,&start,&diff);
           printf("Time to get coupling matrices =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
    }

}

template <class T>
void MasterPola <T>::make_mll_blocks_c_fast( ) {
     double sum_TT, sum_TE, sum_EE_EE, sum_EE_BB, sum_EB;
     double *wigner0, *wigner2,C3, *normfact;
     int kthr;
     long ndim, maxl3,minl3_odd,minl3_even,offset_thr,l1,l2;
     long Nls=this->lLmax+1;
     struct timeval start,end,diff;
     if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);

     if(this->Nmaps < 1) {
        printf("There should be at least one map at this stage (%d)\n",this->Nmaps);
        exit(EXIT_FAILURE);
     }

    this->MAT_TT_TT.alloc(Nls, Nls);
    if(this->Nmaps == 3) {
        this->MAT_EE_EE.alloc(Nls, Nls);
        this->MAT_EE_BB.alloc(Nls, Nls);
        this->MAT_TE_TE.alloc(Nls, Nls);
        this->MAT_EB_EB.alloc(Nls, Nls);
    }

     /* Allocation set out of the loop */
    #ifdef _OPENMP
        wigner0 = (double *) malloc( (2*this->lLmax+1)*outer_loop_threads*inner_loop_threads * sizeof(double));
        if(this->Verbose) printf("wigner0 allocated\n");
        if(this->Nmaps == 3) {
            if(this->Verbose) printf("Wigner 2 allocated\n");
            wigner2 = (double *) malloc( (2*this->lLmax+1)*outer_loop_threads*inner_loop_threads * sizeof(double));
        }
    #else
        wigner0 = (double *) malloc( (2*this->lLmax+1) * sizeof(double));
        if(this->Verbose) printf("wigner0 allocated\n");
        if(this->Nmaps == 3) {
            if(this->Verbose) printf("Wigner 2 allocated\n");
            wigner2 = (double *) malloc( (2*this->lLmax+1) * sizeof(double));
        }
    #endif
    normfact = (double *) malloc((this->lLmax+1) * sizeof(double));
    for( long l1=0; l1<Nls; l1++) normfact[l1]= (2.*double(l1)+1.)/(4.*M_PI);

    long offset[Nls];
    for(l1=0;l1<Nls;l1++) offset[l1]=l1*Nls;

    //Decompose into two steps: build half matrices and use symmetries, then apply normalization
     /* loop over the matrice elements */
    double *MAT_TT_TT_Buffer=this->MAT_TT_TT.buffer(),*MAT_EE_EE_Buffer= this->MAT_EE_EE.buffer(),*MAT_EE_BB_Buffer= this->MAT_EE_BB.buffer(),*MAT_TE_TE_Buffer= this->MAT_TE_TE.buffer(),*MAT_EB_EB_Buffer= this->MAT_EB_EB.buffer() ;
    #ifdef _OPENMP
        #pragma omp parallel for default(none) private(l1,l2,ndim,maxl3,sum_TT,sum_TE, sum_EE_EE, sum_EE_BB, sum_EB, minl3_even, minl3_odd,C3,kthr,offset_thr) shared(normfact, wigner0, wigner2, Nls, offset, MAT_TT_TT_Buffer, MAT_EE_EE_Buffer, MAT_EE_BB_Buffer, MAT_TE_TE_Buffer, MAT_EB_EB_Buffer) num_threads(outer_loop_threads*inner_loop_threads) schedule(dynamic)
    #endif
    for(l1=0; l1<Nls; l1++) {
        #ifdef _OPENMP
            kthr=omp_get_thread_num();
        #else
            kthr=0;
        #endif
        offset_thr= kthr*(2*this->lLmax+1);
        for(l2=l1; l2<Nls; l2++) {
            /* alloc wigners */
            ndim = l2+l1-abs(l1-l2)+1;
            if(ndim > (2*this->lLmax+1)) printf("BEWARE OF NDIM !!!! [%ld]\n",ndim);

            /* compute wigners */
            wig3j_c( l1, l2, (long) 0,  (long) 0, &wigner0[offset_thr]);
            if(this->Nmaps == 3) wig3j_c( l1, l2,  (long)-2,  (long)2, &wigner2[offset_thr]);

            /* loop on l3 */
            maxl3 = (l1+l2 < this->lLmax) ? l1+l2 : this->lLmax;
            if(this->Nmaps !=3) {
                /* initialization */
                sum_TT    = 0.;
                if(((l1+l2+(long) abs(l1-l2)) %2)==0) minl3_even =abs(l1-l2);
                else minl3_even =(long) abs(l1-l2)+1;
                for( long l3= minl3_even; l3<=maxl3; l3+=2) {
                    sum_TT += MaskTT_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3+offset_thr] * wigner0[l3+offset_thr];
                }
                MAT_TT_TT_Buffer[l2+ offset[l1]] = normfact[l2] * sum_TT;
                MAT_TT_TT_Buffer[l1+ offset[l2]] = normfact[l1] * sum_TT;    //Symmetry of wigner0 with respect to l1 and l2
            } else {
                /* initialization */
                sum_TT    = 0.;
                sum_TE    = 0.;
                sum_EE_EE = 0.;
                sum_EE_BB = 0.;
                sum_EB    = 0.;

                /*for( long l3=abs(l1-l2); l3<=maxl3; l3++) {
                    if( (l1+l2+l3)%2 == 0) sum_TT    += MaskTT_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner0[l3];
                    if( (l1+l2+l3)%2 == 0) sum_TE    += MaskTP_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner2[l3];
                    if( (l1+l2+l3)%2 == 0) sum_EE_EE += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                    if( (l1+l2+l3)%2 != 0) sum_EE_BB += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                    sum_EB += MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
                }*/

                if(((l1+l2+(long) abs(l1-l2)) %2)==0) {
                    minl3_even =(long) abs(l1-l2);
                    minl3_odd =(long) abs(l1-l2)+1;
                } else {
                    minl3_even =(long) abs(l1-l2)+1;
                    minl3_odd =(long)  abs(l1-l2);
                }
                for( long l3= minl3_even; l3<=maxl3; l3+=2) {
                    sum_TT    += MaskTT_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3+offset_thr] * wigner0[l3+offset_thr];
                    sum_TE    += MaskTP_Powspec[l3] * (double)(2.*l3+1.) * wigner0[l3+offset_thr] * wigner2[l3+offset_thr];
                    C3=MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3+offset_thr] * wigner2[l3+offset_thr];
                    sum_EE_EE +=C3 ;
                    sum_EB += C3;
                }
                for( long l3= minl3_odd; l3<=maxl3; l3+=2) {
                    C3= MaskPP_Powspec[l3] * (double)(2.*l3+1.) * wigner2[l3+offset_thr] * wigner2[l3+offset_thr];
                    sum_EE_BB +=C3;
                    sum_EB +=C3;
                }
                MAT_TT_TT_Buffer[l2+offset[l1]] = normfact[l2]* sum_TT; //BEWARE: M_ll' = M(l',l) (first index=column)
                MAT_EE_EE_Buffer[l2+offset[l1]]  = normfact[l2] * sum_EE_EE;
                MAT_EE_BB_Buffer[l2+offset[l1]]  = normfact[l2] * sum_EE_BB;
                MAT_TE_TE_Buffer[l2+offset[l1]] = normfact[l2] * sum_TE;
                MAT_EB_EB_Buffer[l2+offset[l1]]  = normfact[l2] * sum_EB;
                if(l1 != l2) {//wigner2 and wigner0 symmetric with respect to l1 and l2
                    MAT_TT_TT_Buffer[l1+offset[l2]]  = normfact[l1]* sum_TT;
                    MAT_EE_EE_Buffer[l1+offset[l2]]  = normfact[l1] * sum_EE_EE;
                    MAT_EE_BB_Buffer[l1+offset[l2]] = normfact[l1] * sum_EE_BB;
                    MAT_TE_TE_Buffer[l1+offset[l2]]  = normfact[l1] * sum_TE;
                    MAT_EB_EB_Buffer[l1+offset[l2]] = normfact[l1] * sum_EB;
                }
            }
            for( long l3=0; l3< ndim; l3++) wigner0[l3+offset_thr]=0.;
            if(this->Nmaps == 3) for( long l3=0; l3< ndim; l3++) wigner2[l3+offset_thr]=0.;
        } //end loop l2
    } //end loop l1
    free(wigner0);
    free(normfact);
    if(this->Nmaps == 3) free(wigner2);

    if(this->Timer==True) {
        gettimeofday(&end, (struct timezone *) NULL);
        timersub(&end,&start,&diff);
           printf("Time to get coupling matrices =%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
    }

}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::Mult_cl_matmask_1block(const arr<double> & ClIn,const dblarray & MatMask,arr<double> & ClOut,bool transpose) {

    int Nls = ClIn.size();
    long offset,kl1,kl2;
    double* MatMask_buffer= MatMask.buffer(), *ClInBuffer=(double *) ClIn.begin(), *ClOutBuffer=(double *) ClOut.begin();
    if((ClOut.size()!= (unsigned int) Nls) || (MatMask.nx() != Nls) || (MatMask.ny() != Nls)) {
        printf("Size of Matrices [%d,%d] and vectors [%d] or [%d] do not agree\n",MatMask.nx(),MatMask.ny(),Nls,(int) ClOut.size());
        exit(EXIT_FAILURE);
    }

    if(transpose) {

        for(kl1=0;kl1<Nls;kl1++) ClOutBuffer[kl1]=0.;
        #ifdef _OPENMP
            #pragma omp parallel for default(none) private(kl1,kl2) shared(Nls, MatMask_buffer, ClInBuffer, ClOutBuffer)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
        #endif
        for(kl1=0;kl1<Nls;kl1++) {
            for(kl2=0;kl2<Nls;kl2++)  ClOutBuffer[kl1]+= MatMask_buffer[kl1+kl2*Nls]* ClInBuffer[kl2];  //BEWARE: M_ll' = M(l',l) (first index=column)
        }
    } else {
        #ifdef _OPENMP
            #pragma omp parallel for default(none)  private(offset,kl1,kl2) shared(Nls, MatMask_buffer, ClInBuffer, ClOutBuffer) num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
        #endif
        for(kl1=0;kl1<Nls;kl1++) {
            ClOutBuffer[kl1]=0.;
            offset=kl1*Nls;
            for(kl2=0;kl2<Nls;kl2++) ClOutBuffer[kl1]+= MatMask_buffer[kl2+ offset]* ClInBuffer[kl2];
        }
    }
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::get_MASTER_pspec(unsigned short xsubi[3], PSPEC_TYPE selspec , bool iter, int Niter,bool Positivity, bool Fast_Radius) {
        int Nls=this->lLmax+1.;
        long kl2;

    //Note: the fact that entire rows/columns of the polarized matrices are equal to 0 is not a problem with the iterative approach (monopole/dipole will be zero), but can generate instabilities for SVD inversion due to a non-trivial null-space
        if((!PMapFlag)) {
            printf("Power spectra of the masked maps should first be computed\n");
            exit(EXIT_FAILURE);
        }
        if(Niter <= 0) Niter=DEF_NITER_ITER_INV;

        switch (selspec)
        {
            case PSPEC_TT: {
                    if(this->Verbose)  printf("PSPEC_TT processing\n");
                    arr<double> PS_TT_2process(Nls);
                    if(PNoiseMapFlag) {
                        if(this->Verbose)  printf("First Subtract Masked Noise Power Spectra\n");
                        for(int kl=0;kl<Nls;kl++)  PS_TT_2process[kl] =this->MskMapPowSpec.tt_[kl]-this->MskNoisePowSpec.tt_[kl];
                    } else {
                        if(this->Verbose)  printf("No Noise Subtraction\n");
                        PS_TT_2process=this->MskMapPowSpec.tt_;
                    }
                    if ( iter ) {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if(this->gamma_TT_TT==0) this->gamma_TT_TT=Check_Matrix_2Norm(this->MAT_TT_TT, xsubi);
                        double mu=2. / (1.1 * gamma_TT_TT* gamma_TT_TT);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TT SpecRad (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                        if(this->Verbose) printf("MasterPola:Spectral Radius for TTTT matrix: %3.5g, relaxation parameter: %3.5g [%d]\n",gamma_TT_TT, mu, Niter);
                        SingleBlock_deconv_iter(PS_TT_2process,this-> Master_TT_Powspec,this->MAT_TT_TT, Niter,mu, Positivity,this->ZeroFirst2L);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TT Powspec (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    } else {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if((this->inv_MAT_TT_TT.nx() != Nls)||(this->inv_MAT_TT_TT.ny()!= Nls)) {
                            if(this->ZeroFirst2L) {
                                dblarray TempMat(Nls-2,Nls-2), invTempMat; //Note: TempMat(Nls-2,Nls) should theoretically be used, but because of large cosmic variance - even for a full sky field with zero monopole and dipole, the masked sky has high probability of having quite large non-zero monopole or dipoles, we do not want even to take into account contribution from l>=2 multipoles to the monopole/dipole obtained after masking.
                                double *TempBuffer= TempMat.buffer(), *MAT_TT_TT_Buffer=this->MAT_TT_TT.buffer();
                                long offset_Nlsm2[Nls-2], offset_Nls[Nls],offset1,offset2;
                                for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                                for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, MAT_TT_TT_Buffer,TempBuffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1= offset_Nlsm2[kl2];
                                    offset2=2+offset_Nls[kl2+2];
                                    for(long kl1=0;kl1<Nls-2;kl1++) TempBuffer[kl1+ offset1]= MAT_TT_TT_Buffer[kl1+offset2];
                                }
                                inv_mat_svd(TempMat, invTempMat);
                                TempMat.free();
                                this->inv_MAT_TT_TT.alloc(Nls,Nls);
                                double *invMatBuffer= invTempMat.buffer(), *invMAT_TT_TT_Buffer=this->inv_MAT_TT_TT.buffer();
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, invMAT_TT_TT_Buffer, invMatBuffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=offset_Nlsm2[kl2];
                                    offset2=2+offset_Nls[kl2+2];
                                    for(long kl1=0;kl1<Nls-2;kl1++) {
                                        invMAT_TT_TT_Buffer[kl1+offset2]= invMatBuffer[kl1+ offset1];
                                    }
                                }
                            } else inv_mat_svd(this->MAT_TT_TT, this->inv_MAT_TT_TT);
                        }
                        if(this->Verbose) printf("MasterPola: inv_MAT_TT_TT size: [%d, %d]\n", this->inv_MAT_TT_TT.nx(), this->inv_MAT_TT_TT.ny());
                        this->Master_TT_Powspec.alloc(Nls);
                        this->Master_TT_Powspec.fill(0.);
                        Mult_cl_matmask_1block(PS_TT_2process, this->inv_MAT_TT_TT, this->Master_TT_Powspec);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TT Powspec (SVD)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    }
                }
                break;
            case PSPEC_EE:
            case PSPEC_BB: {
                    if(this->Verbose) printf("PSPEC_EE AND BB processing\n");
                    arr<double> PS_EE_2process(Nls);
                    arr<double> PS_BB_2process(Nls);
                    if(PNoiseMapFlag) {
                        if(this->Verbose)  printf("First Subtract Masked Noise Power Spectra\n");
                        for(int kl=0;kl<Nls;kl++)  PS_EE_2process[kl] =this->MskMapPowSpec.ee_[kl]-this->MskNoisePowSpec.ee_[kl];
                        for(int kl=0;kl<Nls;kl++)  PS_BB_2process[kl] =this->MskMapPowSpec.bb_[kl]-this->MskNoisePowSpec.bb_[kl];
                    } else {
                        if(this->Verbose)  printf("No Noise Subtraction\n");
                        PS_EE_2process =this->MskMapPowSpec.ee_;
                        PS_BB_2process =this->MskMapPowSpec.bb_;
                    }
                    if ( iter ) {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        double mu;
                        if(Fast_Radius) { //Compute upper bound for spectral radius from sub_matrices
                            if(this->gamma_EE_EE==0) this->gamma_EE_EE=Check_Matrix_2Norm(this->MAT_EE_EE, xsubi);
                            if(this->gamma_EE_BB==0) this->gamma_EE_BB=Check_Matrix_2Norm(this->MAT_EE_BB, xsubi);
                            mu=2. / (1.1 * (this->gamma_EE_EE+ this->gamma_EE_BB)* (this->gamma_EE_EE+ this->gamma_EE_BB));
                            if(this->Verbose) printf("MasterPola: FAST Spectral Radius for EEEE matrix: %3.5g, relaxation parameter: %3.5g\n",this->gamma_EE_EE, mu);
                            if(this->Verbose) printf("MasterPola: FAST Spectral Radius for EEBB matrix: %3.5g, relaxation parameter: %3.5g\n",this->gamma_EE_BB, mu);
                        } else {//Compute upper bound for spectral radius from the total matrix
                            if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                            if(this->gamma_EB_4BLOCKS == 0) {
                                dblarray MAT_EB_4Blocks(2* Nls,2* Nls);
                                double *Mat_EB_4Buffer= MAT_EB_4Blocks.buffer(), *MAT_EE_EE_Buffer=this->MAT_EE_EE.buffer(), *MAT_EE_BB_Buffer =this->MAT_EE_BB.buffer();
                                long offset_2Nls[2*Nls], offset_Nls[Nls];
                                for(kl2=0;kl2<2*Nls;kl2++) offset_2Nls[kl2]=kl2*2*Nls;
                                for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                                long offset1,offset2,offset3,offset4,offset5;

                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,offset3,offset4,offset5,kl2) shared(Nls, MAT_EE_EE_Buffer, MAT_EE_BB_Buffer, Mat_EB_4Buffer,offset_2Nls, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls;kl2++) {
                                    offset1=offset_Nls[kl2];
                                    offset2=offset_2Nls[kl2];
                                    offset3=offset_2Nls[kl2+Nls];
                                    offset4=Nls+offset2;
                                    offset5=Nls+offset3;
                                    for(long kl1=0;kl1<Nls;kl1++) Mat_EB_4Buffer[kl1+offset2]= MAT_EE_EE_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls;kl1++) Mat_EB_4Buffer[kl1+offset4]= MAT_EE_BB_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls;kl1++) Mat_EB_4Buffer[kl1+offset3]= MAT_EE_BB_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls;kl1++) Mat_EB_4Buffer[kl1+offset5]= MAT_EE_EE_Buffer[kl1+offset1];
                                }
/*                                for(int kl2=0;kl2<Nls;kl2++) {
                                    for(int kl1=0;kl1<Nls;kl1++) {
                                        MAT_EB_4Blocks(kl1,kl2)=this->MAT_EE_EE(kl1,kl2);
                                        MAT_EB_4Blocks(kl1+Nls,kl2)=this->MAT_EE_BB(kl1,kl2);
                                        MAT_EB_4Blocks(kl1,kl2+Nls)=this->MAT_EE_BB(kl1,kl2);
                                        MAT_EB_4Blocks(kl1+Nls,kl2+Nls)=this->MAT_EE_EE(kl1,kl2);
                                    }
                                }*/
                                this->gamma_EB_4BLOCKS =Check_Matrix_2Norm(MAT_EB_4Blocks, xsubi);
                            }
                            mu=2. / (1.1 * gamma_EB_4BLOCKS* gamma_EB_4BLOCKS);
                            if(this->Verbose) printf("MasterPola:Spectral Radius for EB_4Blocks matrix: %3.5g, relaxation parameter: %3.5g\n", gamma_EB_4BLOCKS, mu);
                        }
                        if(this->Verbose) printf("Done with radius\n");
                        this->Master_EE_Powspec.alloc(Nls);
                        this->Master_EE_Powspec.fill(0.);
                        this->Master_BB_Powspec.alloc(Nls);
                        this->Master_BB_Powspec.fill(0.);
                        FourSymBlock_deconv_iter(PS_EE_2process, PS_BB_2process, this->Master_EE_Powspec, this->Master_BB_Powspec, this->MAT_EE_EE, this->MAT_EE_BB, Niter,mu,Positivity);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master EE and BB Powspec (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    } else {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if((this-> inv_MAT_EB_4Blocks.nx() != 2*Nls)||(this-> inv_MAT_EB_4Blocks.ny()!= 2*Nls)) {
                            if(Fast_Radius) { //Compute inverse matrices from sub_matrices: ((A,B),(B,A)) where A= MAT_EE_EE, B= MAT_EE_BB
                                dblarray Mat_EE_EE_r(Nls-2,Nls-2),Mat_EE_BB_r(Nls-2,Nls-2),inv_Mat_EE_EE_r;
                                double *Mat_EE_EE_r_Buffer= Mat_EE_EE_r.buffer(),*Mat_EE_BB_r_Buffer= Mat_EE_BB_r.buffer(),*Mat_EE_EE_Buffer= this->MAT_EE_EE.buffer(),*Mat_EE_BB_Buffer=this->MAT_EE_BB.buffer();
                                long offset1,offset2;
                                long offset_Nls[Nls],offset_Nlsm2[Nls-2];
                                for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                                for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                            #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, Mat_EE_EE_Buffer, Mat_EE_EE_r_Buffer, Mat_EE_BB_r_Buffer, Mat_EE_BB_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1= offset_Nlsm2[kl2];
                                    offset2=2+offset_Nls[kl2+2];
                                    for(long kl1=0;kl1<Nls-2;kl1++)  {
                                        Mat_EE_EE_r_Buffer[kl1+offset1]=Mat_EE_EE_Buffer[kl1+offset2];
                                        Mat_EE_BB_r_Buffer[kl1+offset1]=Mat_EE_BB_Buffer[kl1+offset2];
                                    }
                                }
                                inv_mat_svd(Mat_EE_EE_r, inv_Mat_EE_EE_r);//A^-1
                                (this->inv_MAT_EE_EE).alloc(Nls,Nls);
                                double *inv_Mat_EE_EE_r_Buffer= inv_Mat_EE_EE_r.buffer(),*inv_Mat_EE_EE_Buffer= this->inv_MAT_EE_EE.buffer();

                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, inv_Mat_EE_EE_Buffer, inv_Mat_EE_EE_r_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=offset_Nlsm2[kl2];
                                    offset2=2+offset_Nls[kl2+2];
                                    for(long kl1=0;kl1<Nls-2;kl1++)  inv_Mat_EE_EE_Buffer[kl1+offset2]= inv_Mat_EE_EE_r_Buffer[kl1+offset1];
                                }
                                dblarray tempMat1(Nls-2,Nls-2);
                                double *TempMat1Buffer= tempMat1.buffer();

                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, TempMat1Buffer, Mat_EE_BB_r_Buffer, inv_Mat_EE_EE_r_Buffer, offset_Nlsm2)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=offset_Nlsm2[kl2];
                                    for(long kl3=0;kl3<Nls-2;kl3++) {
                                        offset2=offset_Nlsm2[kl3];
                                        for(long kl1=0;kl1<Nls-2;kl1++) TempMat1Buffer[kl1+offset1] +=Mat_EE_BB_r_Buffer[kl3+offset1]* inv_Mat_EE_EE_r_Buffer[kl1+offset2];
                                    }
                                } //mat_mult(Mat_EE_BB_r, inv_Mat_EE_EE_r, tempMat1); //BA^-1


                                dblarray tempMat2(Nls-2,Nls-2);
                                double *TempMat2Buffer= tempMat2.buffer();
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, TempMat2Buffer, TempMat1Buffer, Mat_EE_BB_r_Buffer, Mat_EE_EE_r_Buffer, offset_Nlsm2)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=offset_Nlsm2[kl2];
                                    for(long kl3=0;kl3<Nls-2;kl3++) {
                                        offset2=offset_Nlsm2[kl3];
                                        for(long kl1=0;kl1<Nls-2;kl1++) TempMat2Buffer[kl1+offset1] += TempMat1Buffer[kl3+offset1]* Mat_EE_BB_r_Buffer[kl1+offset2];
                                    }
                                    for(long kl1=0;kl1<Nls-2;kl1++) TempMat2Buffer[kl1+offset1]= Mat_EE_EE_r_Buffer[kl1+offset1]- TempMat2Buffer[kl1+offset1];
                                }
                                //mat_mult(tempMat1, Mat_EE_BB_r, tempMat2); //B A^-1 B
                                //for(int kl1=0;kl1<Nls-2;kl1++) for(int kl2=0;kl2<Nls-2;kl2++) tempMat2(kl1,kl2)= Mat_EE_EE_r(kl1,kl2)-tempMat2(kl1,kl2);//A - BA^-1B

                                dblarray BlockDiag;
                                inv_mat_svd(tempMat2, BlockDiag); //(A - BA^-1B)^-1

                                dblarray BlockOffDiag(Nls-2,Nls-2);
                                double *BlockOffDiagBuffer= BlockOffDiag.buffer(), *BlockDiagBuffer= BlockDiag.buffer();
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, BlockOffDiagBuffer, BlockDiagBuffer, TempMat1Buffer, offset_Nlsm2)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=offset_Nlsm2[kl2];
                                    for(long kl3=0;kl3<Nls-2;kl3++) {
                                        offset2=offset_Nlsm2[kl3];
                                        for(long kl1=0;kl1<Nls-2;kl1++) BlockOffDiagBuffer[kl1+offset1] += BlockDiagBuffer[kl3+offset1]* TempMat1Buffer[kl1+offset2];
                                    }
                                } //mat_mult(BlockDiag, tempMat1, BlockOffDiag); //(A-BA^-1B)^-1 B A^-1

                                long offset_2Nls[2*Nls];
                                for(kl2=0;kl2<2*Nls;kl2++) offset_2Nls[kl2]=kl2*2*Nls;

                                long offset3, offset4,offset5;
                                (this->inv_MAT_EB_4Blocks).alloc(2*Nls,2*Nls);
                                double *inv_MAT_EB_4BlocksBuffer= inv_MAT_EB_4Blocks.buffer();
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,offset3,offset4,offset5,kl2) shared(Nls, BlockDiagBuffer, inv_MAT_EB_4BlocksBuffer, BlockOffDiagBuffer, offset_Nlsm2, offset_2Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1= offset_Nlsm2[kl2];
                                    offset2=2+ offset_2Nls[kl2+2];
                                    offset3=offset2+Nls;
                                    offset4=2+ offset_2Nls[kl2+2+Nls];
                                    offset5= offset4+Nls;
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+offset2]= BlockDiagBuffer[kl1+ offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+offset3]=-BlockOffDiagBuffer[kl1+ offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+offset4]=-BlockOffDiagBuffer[kl1+ offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+offset5]= BlockDiagBuffer[kl1+ offset1];
                                }
                                /*for(int kl1=0;kl1<Nls-2;kl1++) {
                                    for(int kl2=0;kl2<Nls-2;kl2++) {
                                        this->inv_MAT_EB_4Blocks(kl1+2,kl2+2)= BlockDiag(kl1,kl2);
                                        this->inv_MAT_EB_4Blocks(kl1+Nls+2,kl2+2)=-BlockOffDiag(kl1,kl2);
                                        this->inv_MAT_EB_4Blocks(kl1+2,kl2+Nls+2)=-BlockOffDiag(kl1,kl2);
                                        this->inv_MAT_EB_4Blocks(kl1+Nls+2,kl2+Nls+2)= BlockDiag(kl1,kl2);
                                    }
                                }*/
                                if(this->Verbose) printf("MasterPola: FAST Inversion computed for EB block matrix\n");

                            } else {
                                dblarray MAT_EB_4Blockst_nomonop_nodip(2* (Nls-2),2* (Nls-2));
                                double *MAT_EB_4Blockst_nomonop_nodipBuffer= MAT_EB_4Blockst_nomonop_nodip.buffer();
                                double *Mat_EE_EE_Buffer= this->MAT_EE_EE.buffer(),*Mat_EE_BB_Buffer=this->MAT_EE_BB.buffer();

                                long offset_Nls[Nls],offset_2Nls[2*Nls],offset_Nlsm2[Nls-2],offset_2Nlsm4[2* (Nls-2)];
                                for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                                for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;
                                for(kl2=0;kl2<2*Nls;kl2++) offset_2Nls[kl2]=kl2*2*Nls;
                                for(kl2=0;kl2<2* (Nls-2);kl2++) offset_2Nlsm4[kl2]=kl2*2* (Nls-2);

                                long offset1,offset2,offset3, offset4,offset5;
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,offset3,offset4,offset5,kl2) shared(Nls, MAT_EB_4Blockst_nomonop_nodipBuffer, Mat_EE_EE_Buffer, Mat_EE_BB_Buffer, offset_2Nlsm4, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1=2+ offset_Nls[kl2+2];
                                    offset2=offset_2Nlsm4[kl2];
                                    offset3=offset_2Nlsm4[kl2+Nls-2];
                                    offset4=offset2+Nls-2;
                                    offset5=offset3+Nls-2;
                                    for(long kl1=0;kl1<Nls-2;kl1++) MAT_EB_4Blockst_nomonop_nodipBuffer[kl1+offset2]= Mat_EE_EE_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) MAT_EB_4Blockst_nomonop_nodipBuffer[kl1+offset4]= Mat_EE_BB_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) MAT_EB_4Blockst_nomonop_nodipBuffer[kl1+offset3]= Mat_EE_BB_Buffer[kl1+offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) MAT_EB_4Blockst_nomonop_nodipBuffer[kl1+offset5]= Mat_EE_EE_Buffer[kl1+offset1];
                                }

                                dblarray invMAT_EB_4Blockst_nomonop_nodip(2* (Nls-2),2* (Nls-2));
                                inv_mat_svd(MAT_EB_4Blockst_nomonop_nodip, invMAT_EB_4Blockst_nomonop_nodip);
                                (this->inv_MAT_EB_4Blocks).alloc(2*Nls,2*Nls);
                                double *invMAT_EB_4Blockst_nomonop_nodipBuffer= invMAT_EB_4Blockst_nomonop_nodip.buffer();
                                double *inv_MAT_EB_4BlocksBuffer= this->inv_MAT_EB_4Blocks.buffer();
                                long offset6,offset7,offset8;
                                #ifdef _OPENMP
                                    #pragma omp parallel for default(none) private(offset1,offset2,offset3,offset4,offset5,offset6,offset7,offset8,kl2) shared(Nls, inv_MAT_EB_4BlocksBuffer, invMAT_EB_4Blockst_nomonop_nodipBuffer, offset_2Nlsm4, offset_2Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                                #endif
                                for(kl2=0;kl2<Nls-2;kl2++) {
                                    offset1= offset_2Nlsm4[kl2];
                                    offset2= offset_2Nlsm4[kl2+Nls-2];
                                    offset3=offset1+Nls-2;
                                    offset4=offset2+Nls-2;
                                    offset5=2+ offset_2Nls[kl2+2];
                                    offset6=2+ offset_2Nls[kl2+2+Nls];
                                    offset7=offset5+Nls;
                                    offset8=offset6+Nls;
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+ offset5]= invMAT_EB_4Blockst_nomonop_nodipBuffer[kl1+ offset1];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+ offset7]= invMAT_EB_4Blockst_nomonop_nodipBuffer[kl1+ offset3];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+ offset6]= invMAT_EB_4Blockst_nomonop_nodipBuffer[kl1+ offset2];
                                    for(long kl1=0;kl1<Nls-2;kl1++) inv_MAT_EB_4BlocksBuffer[kl1+ offset8]= invMAT_EB_4Blockst_nomonop_nodipBuffer[kl1+ offset4];
                                }
                                if(this->Verbose) printf("MasterPola: Inversion computed for EB block matrix\n");
                            }
                        }
                        if(this->Verbose) printf("MasterPola: inv_MAT_EB_4Blocks size: [%d, %d]\n", this->inv_MAT_EB_4Blocks.nx(), this->inv_MAT_EB_4Blocks.ny());
                        arr<double> ClsIn(2*Nls), ClsOut(2*Nls);
                        for(long kl=0;kl<Nls;kl++) {
                            ClsIn[kl]= PS_EE_2process[kl];
                            ClsIn[kl+Nls]= PS_BB_2process[kl];
                        }
                        Mult_cl_matmask_1block(ClsIn,this->inv_MAT_EB_4Blocks,ClsOut);
                        this->Master_EE_Powspec.alloc(Nls);
                        this->Master_BB_Powspec.alloc(Nls);
                        for(long kl=0;kl<Nls;kl++) {
                            this->Master_EE_Powspec[kl]= ClsOut[kl];
                            this->Master_BB_Powspec[kl]= ClsOut[kl+Nls];
                        }
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master EE and BB Powspec (SVD)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    }
                }
                break;
            case PSPEC_TE: {
                    if(this->Verbose) printf("PSPEC_TE processing\n");
                    arr<double> PS_TE_2process(Nls);
                    if((PNoiseMapFlag)&&(!FlagUncorrNoise)) {
                        if(this->Verbose)  printf("First Subtract Masked Noise Power Spectra\n");
                        for(int kl=0;kl<Nls;kl++)  PS_TE_2process[kl] =this->MskMapPowSpec.te_[kl]-this->MskNoisePowSpec.te_[kl];
                    } else {
                        if(this->Verbose)  printf("No Noise Subtraction\n");
                        PS_TE_2process=this->MskMapPowSpec.te_;
                    }
                    if ( iter ) {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if(this->gamma_TE_TE==0) this->gamma_TE_TE=Check_Matrix_2Norm(this->MAT_TE_TE, xsubi);
                        double mu=2. / (1.1 * this->gamma_TE_TE* this->gamma_TE_TE);
                        if(this->Verbose) printf("MasterPola:Spectral Radius for TETE matrix: %3.5g, relaxation parameter: %3.5g\n",this->gamma_TE_TE, mu);
                        SingleBlock_deconv_iter(PS_TE_2process, this->Master_TE_Powspec, this->MAT_TE_TE, Niter,mu, Positivity,this-> ZeroFirst2L);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TE Powspec (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    } else {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if((this->inv_MAT_TE_TE.nx() != Nls)||(this->inv_MAT_TE_TE.ny()!= Nls)) {
                            dblarray Mat_TE_TE_r(Nls-2,Nls-2),invMat_TE_TE_r(Nls-2,Nls-2);
                            (this->inv_MAT_TE_TE).alloc(Nls,Nls);
                            double *Mat_TE_TE_r_Buffer= Mat_TE_TE_r.buffer(),*Mat_TE_TE_Buffer=this->MAT_TE_TE.buffer();
                            long offset1,offset2;
                            long offset_Nls[Nls],offset_Nlsm2[Nls-2];
                            for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                            for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, Mat_TE_TE_r_Buffer, Mat_TE_TE_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) Mat_TE_TE_r_Buffer[kl1+offset1]= Mat_TE_TE_Buffer[kl1+ offset2];
                            }
                            inv_mat_svd(Mat_TE_TE_r, invMat_TE_TE_r);

                            double *invMat_TE_TE_r_Buffer=invMat_TE_TE_r.buffer(),*invMat_TE_TE_Buffer=this->inv_MAT_TE_TE.buffer();
                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, invMat_TE_TE_r_Buffer, invMat_TE_TE_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) invMat_TE_TE_Buffer[kl1+offset2]= invMat_TE_TE_r_Buffer[kl1+offset1];
                            }
                        }
                        if(this->Verbose) printf("MasterPola: inv_MAT_TE_TE size: [%d, %d]\n", this->inv_MAT_TE_TE.nx(), this->inv_MAT_TE_TE.ny());
                        this->Master_TE_Powspec.alloc(Nls);
                        Mult_cl_matmask_1block(PS_TE_2process, this->inv_MAT_TE_TE, this->Master_TE_Powspec);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TE Powspec (SVD)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    }
                }
                break;
            case PSPEC_TB: {
                    if(this->Verbose) printf("PSPEC_TB processing\n");
                    arr<double> PS_TB_2process(Nls);
                    if((PNoiseMapFlag)&&(!FlagUncorrNoise)) {
                        if(this->Verbose)  printf("First Subtract Masked Noise Power Spectra\n");
                        for(int kl=0;kl<Nls;kl++)  PS_TB_2process[kl] =this->MskMapPowSpec.tb_[kl]-this->MskNoisePowSpec.tb_[kl];
                    } else {
                        if(this->Verbose)  printf("No Noise Subtraction\n");
                        PS_TB_2process=this->MskMapPowSpec.tb_;
                    }
                    if ( iter ) {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if(this->gamma_TE_TE==0) this->gamma_TE_TE=Check_Matrix_2Norm(this->MAT_TE_TE, xsubi);
                        double mu=2. / (1.1 * this->gamma_TE_TE* this->gamma_TE_TE);
                        if(this->Verbose) printf("MasterPola:Spectral Radius for TETE matrix: %3.5g, relaxation parameter: %3.5g\n",this->gamma_TE_TE, mu);
                        SingleBlock_deconv_iter(PS_TB_2process, this->Master_TB_Powspec, this->MAT_TE_TE, Niter,mu, Positivity,this-> ZeroFirst2L);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TB Powspec (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    } else {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if((this->inv_MAT_TE_TE.nx() != Nls)||(this->inv_MAT_TE_TE.ny()!= Nls)) {

                            dblarray Mat_TE_TE_r(Nls-2,Nls-2),invMat_TE_TE_r(Nls-2,Nls-2);
                            (this->inv_MAT_TE_TE).alloc(Nls,Nls);
                            double *Mat_TE_TE_r_Buffer= Mat_TE_TE_r.buffer(),*Mat_TE_TE_Buffer=this->MAT_TE_TE.buffer();
                            long offset1,offset2;
                            long offset_Nls[Nls],offset_Nlsm2[Nls-2];
                            for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                            for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, Mat_TE_TE_r_Buffer, Mat_TE_TE_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) Mat_TE_TE_r_Buffer[kl1+offset1]= Mat_TE_TE_Buffer[kl1+ offset2];
                            }
                            inv_mat_svd(Mat_TE_TE_r, invMat_TE_TE_r);

                            double *invMat_TE_TE_r_Buffer=invMat_TE_TE_r.buffer(),*invMat_TE_TE_Buffer=this->inv_MAT_TE_TE.buffer();
                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, invMat_TE_TE_r_Buffer, invMat_TE_TE_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) invMat_TE_TE_Buffer[kl1+offset2]= invMat_TE_TE_r_Buffer[kl1+offset1];
                            }
                        }
                        if(this->Verbose) printf("MasterPola: inv_MAT_TE_TE size: [%d, %d]\n", this->inv_MAT_TE_TE.nx(), this->inv_MAT_TE_TE.ny());
                        this->Master_TB_Powspec.alloc(Nls);
                        Mult_cl_matmask_1block(PS_TB_2process, this->inv_MAT_TE_TE, this->Master_TB_Powspec);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master TB Powspec (SVD)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    }
                }
                break;
            case PSPEC_EB: {
                    if(this->Verbose) printf("PSPEC_EB processing\n");
                    arr<double> PS_EB_2process(Nls);
                    if((PNoiseMapFlag)&&(!FlagUncorrNoise)) {
                        if(this->Verbose)  printf("First Subtract Masked Noise Power Spectra\n");
                        for(int kl=0;kl<Nls;kl++)  PS_EB_2process[kl] =this->MskMapPowSpec.eb_[kl]-this->MskNoisePowSpec.eb_[kl];
                    } else {
                        if(this->Verbose)  printf("No Noise Subtraction\n");
                        PS_EB_2process =this->MskMapPowSpec.eb_;
                    }
                    if ( iter ) {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if(this->gamma_EB_EB==0) this->gamma_EB_EB=Check_Matrix_2Norm(this->MAT_EB_EB, xsubi);
                        double mu=2. / (1.1 * this-> gamma_EB_EB* this-> gamma_EB_EB);
                        if(this->Verbose) printf("MasterPola:Spectral Radius for TETE matrix: %3.5g, relaxation parameter: %3.5g\n",this->gamma_EB_EB, mu);
                        SingleBlock_deconv_iter(PS_EB_2process, this->Master_EB_Powspec, this->MAT_EB_EB, Niter,mu, Positivity,this-> ZeroFirst2L);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master EB Powspec (ITER)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    } else {
                        struct timeval start,end,diff;
                        if(this->Timer==True) gettimeofday(&start, (struct timezone *) NULL);
                        if((this->inv_MAT_EB_EB.nx() != Nls)||(this->inv_MAT_EB_EB.ny()!= Nls)) {

                            dblarray Mat_EB_EB_r(Nls-2,Nls-2), invMat_EB_EB_r(Nls-2,Nls-2);
                            (this-> inv_MAT_EB_EB).alloc(Nls,Nls);
                            double *Mat_EB_EB_r_Buffer= Mat_EB_EB_r.buffer(),*Mat_EB_EB_Buffer=this->MAT_EB_EB.buffer();
                            long offset1,offset2;
                            long offset_Nls[Nls],offset_Nlsm2[Nls-2];
                            for(kl2=0;kl2<Nls-2;kl2++) offset_Nlsm2[kl2]=kl2*(Nls-2);
                            for(kl2=0;kl2<Nls;kl2++) offset_Nls[kl2]=kl2*Nls;

                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, Mat_EB_EB_r_Buffer, Mat_EB_EB_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) Mat_EB_EB_r_Buffer[kl1+offset1]= Mat_EB_EB_Buffer[kl1+ offset2];
                            }
                            inv_mat_svd(Mat_EB_EB_r, invMat_EB_EB_r);

                            double *invMat_EB_EB_r_Buffer=invMat_EB_EB_r.buffer(),*invMat_EB_EB_Buffer=this->inv_MAT_EB_EB.buffer();
                            #ifdef _OPENMP
                                #pragma omp parallel for default(none) private(offset1,offset2,kl2) shared(Nls, invMat_EB_EB_r_Buffer, invMat_EB_EB_Buffer, offset_Nlsm2, offset_Nls)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
                            #endif
                            for(kl2=0;kl2<Nls-2;kl2++) {
                                offset1= offset_Nlsm2[kl2];
                                offset2=2+ offset_Nls[kl2+2];
                                for(long kl1=0;kl1<Nls-2;kl1++) invMat_EB_EB_Buffer[kl1+offset2]= invMat_EB_EB_r_Buffer[kl1+offset1];
                            }
                        }
                        if(this->Verbose) printf("MasterPola: inv_MAT_EB_EB size: [%d, %d]\n", this->inv_MAT_EB_EB.nx(), this->inv_MAT_EB_EB.ny());
                        this->Master_EB_Powspec.alloc(Nls);
                        Mult_cl_matmask_1block(PS_EB_2process, this->inv_MAT_EB_EB, this->Master_EB_Powspec);
                        if(this->Timer==True) {
                            gettimeofday(&end, (struct timezone *) NULL);
                            timersub(&end,&start,&diff);
                               printf("Time to get Master EB Powspec (SVD)=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
                        }
                    }
                }
                break;
          default:
                 this->get_MASTER_pspec(xsubi, PSPEC_TT , iter, Niter, Positivity, Fast_Radius);
                 this->get_MASTER_pspec(xsubi, PSPEC_EE , iter, Niter, Positivity, Fast_Radius);
                 this->get_MASTER_pspec(xsubi, PSPEC_TE , iter, Niter, Positivity, Fast_Radius);
                 this->get_MASTER_pspec(xsubi, PSPEC_TB , iter, Niter, Positivity, Fast_Radius);
                 this->get_MASTER_pspec(xsubi, PSPEC_EB , iter, Niter, Positivity, Fast_Radius);
    }

}


//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::SingleBlock_deconv_iter(const arr<double> & ClData, arr<double> & ClSol, dblarray & MatMask,const int Niter,const double mu, bool Positivity,bool ZeroFirst2Flag) {

    int kl, Nls=ClData.size();
    arr<double> Fwd(Nls);
    arr<double> Resi(Nls);
    double total_energy=0.;
    arr<double> Prop(Nls);
    Prop.fill(0.);

    if((MatMask.nx() != Nls) || (MatMask.ny() !=Nls)) {
        printf("MasterPola: MatMask size [%d,%d] and Power spectrum size %d do not agree\n",MatMask.nx(),MatMask.ny(),Nls);
        throw 1;
    }
    ClSol.alloc(Nls);
    ClSol.fill(0.);
    for(kl=2;kl<Nls;kl++) if(fabs(MatMask(kl,kl)) > 1e-12 ) ClSol[kl]=ClData[kl]/MatMask(kl,kl); else  ClSol[kl]= ClData[kl];
    //const double *ClDataBuffer=ClData.begin(),*ClSolBuffer=ClSol.begin(), *ResiBuffer= Resi.begin(), *FwdBuffer= Fwd.begin(), *PropBuffer=Prop.begin();
    double *ClDataBuffer=(double *)ClData.begin();
    double *ClSolBuffer=(double *) ClSol.begin(), *ResiBuffer=(double *)Resi.begin(), *FwdBuffer=(double *) Fwd.begin(), *PropBuffer=(double *) Prop.begin();

    if(Positivity) {
        if(this->Verbose) {
            for(int kit=0;kit<Niter;kit++) {
                Mult_cl_matmask_1block(ClSol, MatMask, Fwd);
                for(kl=0;kl<Nls;kl++) ResiBuffer[kl]= ClDataBuffer[kl]-FwdBuffer[kl];
                Mult_cl_matmask_1block(Resi, MatMask, Prop,true);
                for(kl=0;kl<Nls;kl++) ClSolBuffer[kl]=MAX(ClSolBuffer[kl]+mu* PropBuffer[kl],0.);
                if( ZeroFirst2Flag)  {
                    ClSolBuffer[0]=0;
                    ClSolBuffer[1]=0;
                }
                total_energy=0.;
                for(kl=0;kl<Nls;kl++) total_energy += pow(ResiBuffer[kl],2.);
                printf("[%d] Residual Energy=%3.10g\n",kit, total_energy);
            }
        } else {
            for(int kit=0;kit<Niter;kit++) {
                Mult_cl_matmask_1block(ClSol, MatMask, Fwd);
                for(kl=0;kl<Nls;kl++) ResiBuffer[kl]= ClDataBuffer[kl]-FwdBuffer[kl];
                Mult_cl_matmask_1block(Resi, MatMask, Prop,true);
                for(kl=0;kl<Nls;kl++) ClSolBuffer[kl]=MAX(ClSolBuffer[kl]+mu* PropBuffer[kl],0.);
                if( ZeroFirst2Flag)  {
                    ClSolBuffer[0]=0;
                    ClSolBuffer[1]=0;
                }
            }
        }
    } else{
        if(this->Verbose) {
            for(int kit=0;kit<Niter;kit++) {
                Mult_cl_matmask_1block(ClSol, MatMask, Fwd);
                for(kl=0;kl<Nls;kl++) ResiBuffer[kl]= ClDataBuffer[kl]-FwdBuffer[kl];
                Mult_cl_matmask_1block(Resi, MatMask, Prop,true);
                for(kl=0;kl<Nls;kl++) ClSolBuffer[kl]+=mu* PropBuffer[kl];
                if( ZeroFirst2Flag)  {
                    ClSolBuffer[0]=0;
                    ClSolBuffer[1]=0;
                }
                total_energy=0.;
                for(kl=0;kl<Nls;kl++) total_energy += pow(ResiBuffer[kl],2.);
                printf("[%d] Residual Energy=%3.10g\n",kit, total_energy);
            }
        } else {
            for(int kit=0;kit<Niter;kit++) {
                Mult_cl_matmask_1block(ClSol, MatMask, Fwd);
                for(kl=0;kl<Nls;kl++) ResiBuffer[kl]= ClDataBuffer[kl]-FwdBuffer[kl];
                Mult_cl_matmask_1block(Resi, MatMask, Prop,true);
                for(kl=0;kl<Nls;kl++) ClSolBuffer[kl]+=mu* PropBuffer[kl];
                if( ZeroFirst2Flag)  {
                    ClSolBuffer[0]=0;
                    ClSolBuffer[1]=0;
                }
            }
        }
    }
}


//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: FourSymBlock_deconv_iter(const arr<double> & ClData1,const arr<double> & ClData2, arr<double> & ClSol1, arr<double> & ClSol2, dblarray & MatMask1, dblarray & MatMask2, const int Niter, const double mu, bool Positivity) {

    long Nls=ClData1.size(),kl2;
    arr<double> CLIn(2*Nls),ClOut(2*Nls);
    CLIn.fill(0.);
    ClOut.fill(0.);
    if((MatMask1.nx() != (int) Nls) || (MatMask1.ny() !=(int) Nls)) {
        printf("MatMask1 size [%d,%d] and Power spectrum size %d do not agree\n",MatMask1.nx(),MatMask1.ny(),(int) Nls);
        throw 1;
    }
    if((MatMask2.nx() != (int) Nls) || (MatMask2.ny() !=(int) Nls)) {
        printf("MatMask2 size [%d,%d]  and Power spectrum size %d do not agree\n",MatMask2.nx(),MatMask2.ny(),(int) Nls);
        throw 1;
    }
    if((ClData2.size() != (unsigned int) Nls)) {
        printf("Power spectrum size %ld and %d do not agree\n",ClData2.size(),(int) Nls);
        throw 1;
    }

    dblarray MAT_4Blocks(2* Nls,2* Nls);
    double *MatMask1Buffer= MatMask1.buffer(),*MatMask2Buffer= MatMask2.buffer(), *MAT_4BlocksBuffer= MAT_4Blocks.buffer();
    double *ClData1Buffer= (double *) ClData1.begin(),*ClSol1Buffer=(double *) ClSol1.begin(),*ClData2Buffer=(double *) ClData2.begin(),*ClSol2Buffer=(double *) ClSol2.begin(), *CLInBuffer=(double *) CLIn.begin(),*CLOutBuffer =(double *) ClOut.begin();
    long offset1,offset2,offset3,offset4,offset5;

    #ifdef _OPENMP
        #pragma omp parallel for default(none) private(offset1,offset2,offset3,offset4,offset5,kl2) shared(Nls, MAT_4BlocksBuffer, MatMask1Buffer, MatMask2Buffer)  num_threads(outer_loop_threads*inner_loop_threads) schedule(static)
    #endif
    for(kl2=0;kl2<Nls;kl2++) {
        offset1 =kl2*Nls;
        offset2 =kl2*2*Nls;
        offset3 =(kl2+Nls)*2*Nls;
        offset4 = offset2+Nls;
        offset5 = offset3+Nls;
        for(long kl1=0;kl1<Nls;kl1++) MAT_4BlocksBuffer[kl1+ offset2]= MatMask1Buffer[kl1+offset1];
        for(long kl1=0;kl1<Nls;kl1++) MAT_4BlocksBuffer[kl1+ offset4]= MatMask2Buffer[kl1+offset1];
        for(long kl1=0;kl1<Nls;kl1++) MAT_4BlocksBuffer[kl1+ offset3]= MatMask2Buffer[kl1+offset1];
        for(long kl1=0;kl1<Nls;kl1++) MAT_4BlocksBuffer[kl1+ offset5]= MatMask1Buffer[kl1+offset1];
    }

    for(long kl1=0;kl1<Nls;kl1++) {
        CLInBuffer[kl1]= ClData1Buffer[kl1];
        CLInBuffer[kl1+Nls]= ClData2Buffer[kl1];
    }
    SingleBlock_deconv_iter(CLIn, ClOut, MAT_4Blocks, Niter,mu, Positivity);
    ClSol1.alloc(Nls);
    ClSol2.alloc(Nls);
    for(long kl1=0;kl1<Nls;kl1++) {
        ClSol1Buffer[kl1]= CLOutBuffer[kl1];
        ClSol2Buffer[kl1]= CLOutBuffer[kl1+Nls];
    }

}

//****************************************************************************************************************************//
//****************************************************************************************************************************//
// READING ROUTINES
//****************************************************************************************************************************//
//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: read_coupling_mat(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
    int  nfound, anynull,nhdu;
    double nullval;
      long naxes[2],fpixel=1;
       int status = 0,hdutype=0;
      char FitsName[MAXCHAR];
    char *NameFits=fitsname(Name);
     strcpy(FitsName, NameFits);
    free(NameFits);
       if ( fits_open_file(&fptr, FitsName, (int) READONLY, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
    if (fits_get_num_hdus(fptr, &nhdu, &status) ) {
        printf("Error: cannot get number of hdus for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
    if(nhdu != 5) {
        printf("%s is not a proper filename for coupling matrices\n",Name);
        exit(EXIT_FAILURE);
    }

    //READ TTTT inv Matrix
       naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for TTTT inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for TTTT=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->MAT_TT_TT.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->MAT_TT_TT.n_elem(),  (void *) &nullval, this->MAT_TT_TT.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s  for TTTT Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ EEEE inv Matrix
    if(fits_movabs_hdu(fptr, 2, &hdutype,&status)) {
        printf("Error: cannot move to hdu 2 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EEEE Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for EEEE=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->MAT_EE_EE.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->MAT_EE_EE.n_elem(),  (void *) &nullval, this->MAT_EE_EE.buffer(), &anynull, &status) ){
                   printf("\n error in fits_read_img %s for EEEE Matrix", FitsName);
                   exit(EXIT_FAILURE);
            }
    }

    //READ EEBB inv Matrix
    if(fits_movabs_hdu(fptr, 3, &hdutype,&status)) {
        printf("Error: cannot move to hdu 3 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
        if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
            printf(" Error: cannot read NAXIS keyword for EEBB Matrix");
        fits_close_file(fptr, &status);
            exit(EXIT_FAILURE);
    }
    printf("Naxes for EEBB=%ld %ld \n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->MAT_EE_BB.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->MAT_EE_BB.n_elem(),  (void *) &nullval, this->MAT_EE_BB.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EEBB Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ TETE inv Matrix
    if(fits_movabs_hdu(fptr, 4, &hdutype,&status)) {
        printf("Error: cannot move to hdu 2 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword  for TETE Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for TETE=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->MAT_TE_TE.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->MAT_TE_TE.n_elem(),  (void *) &nullval, this->MAT_TE_TE.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s  for TETE Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ EBEB inv Matrix
    if(fits_movabs_hdu(fptr, 5, &hdutype,&status)) {
        printf("Error: cannot move to hdu 5 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EBEB Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for EBEB=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->MAT_EB_EB.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->MAT_EB_EB.n_elem(),  (void *) &nullval, this->MAT_EB_EB.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EBEB Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }

}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::read_spec_radii(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
    int  nfound, anynull;
    double nullval;
      long naxes[1],fpixel=1;
       int status = 0;
      char FitsName[MAXCHAR];
     char *NameFits=fitsname(Name);
    strcpy(FitsName, NameFits);
    free(NameFits);
       if ( fits_open_file(&fptr, FitsName, (int) READONLY, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(-1);
      }
    //SAVE TTTT inv Matrix
       naxes[0] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 1, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    if(naxes[0] != 6) {
        printf("Should have 6 values in file %s instead of %ld\n", Name,naxes[0]);
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    double spec_radii[6];
    nullval= 0;
    if(fits_read_img(fptr, TDOUBLE, fpixel, 6,  (void *) &nullval, spec_radii, &anynull, &status) ){
       printf("\n error in fits_read_img %s", FitsName);
       exit(EXIT_FAILURE);
    }

    this->gamma_TT_TT= spec_radii[0];
    this->gamma_EE_EE= spec_radii[1];
    this->gamma_EE_BB= spec_radii[2];
    this->gamma_TE_TE= spec_radii[3];
    this->gamma_EB_EB= spec_radii[4];
    this->gamma_EB_4BLOCKS= spec_radii[5];

    if(this->Verbose) printf("SpecRadii= TTTT %f, EEEE %f, EEBB %f, TETE %f,EBEB %f, EB_4B %f\n",this->gamma_TT_TT,this->gamma_EE_EE,this->gamma_EE_BB, this->gamma_TE_TE,this->gamma_EB_EB,this->gamma_EB_4BLOCKS);

    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: read_inv_mat(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
    int nfound, anynull,nhdu;
    double nullval;
      long naxes[2],fpixel=1;
       int status = 0,hdutype=0;
      char FitsName[MAXCHAR];
      char *NameFits=fitsname(Name);
    strcpy(FitsName, NameFits);
    free(NameFits);
       if ( fits_open_file(&fptr, FitsName, (int) READONLY, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
    if (fits_get_num_hdus(fptr, &nhdu, &status) ) {
        printf("Error: cannot get number of hdus for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
    if(nhdu != 6) {
        printf("%s is not a proper filename for invert matrices\n",Name);
        exit(EXIT_FAILURE);
    }

    //READ TTTT inv Matrix
       naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for TTTT inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invTTTT=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_TT_TT.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_TT_TT.n_elem(),  (void *) &nullval, this->inv_MAT_TT_TT.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s  for TTTT inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }


    //READ EEEE inv Matrix
    if(fits_movabs_hdu(fptr, 2, &hdutype,&status)) {
        printf("Error: cannot move to hdu 3 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EEEE inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invEEEE=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_EE_EE.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_EE_EE.n_elem(),  (void *) &nullval, this->inv_MAT_EE_EE.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EEEE inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ EEBB inv Matrix
    if(fits_movabs_hdu(fptr, 3, &hdutype,&status)) {
        printf("Error: cannot move to hdu 4 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EEBB inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invEEBB=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_EE_BB.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_EE_BB.n_elem(),  (void *) &nullval, this->inv_MAT_EE_BB.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EEBB inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ TETE inv Matrix
    if(fits_movabs_hdu(fptr, 4, &hdutype,&status)) {
        printf("Error: cannot move to hdu 2 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword  for TETE inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invTETE=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_TE_TE.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_TE_TE.n_elem(),  (void *) &nullval, this->inv_MAT_TE_TE.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s  for TETE inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }
    //READ EBEB inv Matrix
    if(fits_movabs_hdu(fptr, 5, &hdutype,&status)) {
        printf("Error: cannot move to hdu 5 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EBEB inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invEBEB=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_EB_EB.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_EB_EB.n_elem(),  (void *) &nullval, this->inv_MAT_EB_EB.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EBEB inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    //READ EB 4 Blocks inv Matrix
    if(fits_movabs_hdu(fptr, 6, &hdutype,&status)) {
        printf("Error: cannot move to hdu 6 for file %s %d ", FitsName, status);
        exit(EXIT_FAILURE);
      }
      naxes[0] = 0;  // int converted to long
       naxes[1] = 0;  // int converted to long
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {
        printf(" Error: cannot read NAXIS keyword for EB 4 Blocks inv Matrix");
        fits_close_file(fptr, &status);
        exit(EXIT_FAILURE);
    }
    printf("Naxes for invEB 4 Blocks=%ld %ld\n",naxes[0],naxes[1]);
    if(naxes[0]*naxes[1] !=0) {
        this->inv_MAT_EB_4Blocks.alloc(naxes[0],naxes[1]);
        if(fits_read_img(fptr, TDOUBLE, fpixel, this->inv_MAT_EB_4Blocks.n_elem(),  (void *) &nullval, this->inv_MAT_EB_4Blocks.buffer(), &anynull, &status) ){
               printf("\n error in fits_read_img %s for EB 4 Blocks inv Matrix", FitsName);
               exit(EXIT_FAILURE);
        }
    }

    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }

}


//****************************************************************************************************************************//
//****************************************************************************************************************************//
// WRITING ROUTINES
//****************************************************************************************************************************//
//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::write_coupling_mat(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
      long naxes[2], firstpixel=1, npixels;
       int status = 0;
      int bitpix = -64;
      char FitsName[MAXCHAR];
       char *NameFits=fitsname(Name);
    strcpy(FitsName, NameFits);
    free(NameFits);

       remove(FitsName);
       if ( fits_create_file(&fptr, FitsName, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(-1);
      }

    //SAVE TTTT Matrix
       naxes[0] = (long) this->MAT_TT_TT.axis(1);  // int converted to long
       naxes[1] = (long) this->MAT_TT_TT.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->MAT_TT_TT.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->MAT_TT_TT.buffer(),&status) ) {
               printf("Error: cannot write TTTT matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    //SAVE EEEE Matrix
       naxes[0] = (long) this->MAT_EE_EE.axis(1);  // int converted to long
       naxes[1] = (long) this->MAT_EE_EE.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->MAT_EE_EE.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->MAT_EE_EE.buffer(),&status) ) {
               printf("Error: cannot write EEEE matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
       //SAVE EEBB Matrix
       naxes[0] = (long) this->MAT_EE_BB.axis(1);  // int converted to long
       naxes[1] = (long) this->MAT_EE_BB.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->MAT_EE_BB.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->MAT_EE_BB.buffer(),&status) ) {
               printf("Error: cannot write EEBB matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    //SAVE TETE Matrix
       naxes[0] = (long) this->MAT_TE_TE.axis(1);  // int converted to long
       naxes[1] = (long) this->MAT_TE_TE.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->MAT_TE_TE.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->MAT_TE_TE.buffer(),&status) ) {
               printf("Error: cannot write TETE matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    //SAVE EBEB Matrix
       naxes[0] = (long) this->MAT_EB_EB.axis(1);  // int converted to long
       naxes[1] = (long) this->MAT_EB_EB.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->MAT_EB_EB.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->MAT_EB_EB.buffer(),&status) ) {
               printf("Error: cannot write EBEB matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::write_spec_radii(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
      long naxes[1];
       int status = 0;
      int bitpix = -64;
      char FitsName[MAXCHAR];
    char *NameFits=fitsname(Name);
    strcpy(FitsName, NameFits);
    free(NameFits);
       remove(FitsName);
       if ( fits_create_file(&fptr, FitsName, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(-1);
      }
    //SAVE TTTT inv Matrix
       naxes[0] = 6;  // int converted to long

    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
    double spec_radii[6];
    spec_radii[0]= this->gamma_TT_TT;
    spec_radii[1]= this->gamma_EE_EE;
    spec_radii[2]= this->gamma_EE_BB;
    spec_radii[3]= this->gamma_TE_TE;
    spec_radii[4]= this->gamma_EB_EB;
    spec_radii[5]= this->gamma_EB_4BLOCKS;

    if ( fits_write_img( fptr, TDOUBLE, 1, 6, spec_radii,&status) ) {
      printf("Error: cannot write matrix spec. radii %s %d ", FitsName, status);
      exit(-1);
       }
    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::write_inv_mat(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
      long naxes[2], firstpixel=1, npixels;
       int status = 0, hdunum,nhdu,res;
      int bitpix = -64;
      char FitsName[MAXCHAR];
    char *NameFits=fitsname(Name);
     strcpy(FitsName, NameFits);
    free(NameFits);
       remove(FitsName);
       if ( fits_create_file(&fptr, FitsName, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(-1);
      }
    //SAVE TTTT inv Matrix
       naxes[0] = (long) this->inv_MAT_TT_TT.axis(1);  // int converted to long
       naxes[1] = (long) this->inv_MAT_TT_TT.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_TT_TT.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_TT_TT.buffer(),&status) ) {
               printf("Error: cannot write inv TTTT matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("TTTT: %d out of %d\n", hdunum, nhdu);
    }

    //SAVE EEEE inv Matrix
       naxes[0] =(long) this->inv_MAT_EE_EE.axis(1);  // int converted to long
       naxes[1] =(long) this->inv_MAT_EE_EE.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_EE_EE.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_EE_EE.buffer(),&status) ) {
               printf("Error: cannot write inv EEEE matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("EEEE: %d out of %d\n", hdunum, nhdu);
    }

       //SAVE EEBB inv Matrix
       naxes[0] = (long) this->inv_MAT_EE_BB.axis(1);  // int converted to long
       naxes[1] = (long) this->inv_MAT_EE_BB.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_EE_BB.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_EE_BB.buffer(),&status) ) {
               printf("Error: cannot write inv EEBB matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("EEBB: %d out of %d\n", hdunum, nhdu);
    }


    //SAVE TETE inv Matrix
       naxes[0] = (long) this->inv_MAT_TE_TE.axis(1);  // int converted to long
       naxes[1] = (long) this->inv_MAT_TE_TE.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_TE_TE.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_TE_TE.buffer(),&status) ) {
               printf("Error: cannot write inv TETE matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("TETE: %d out of %d\n", hdunum, nhdu);
    }


    //SAVE EBEB inv Matrix
       naxes[0] =  (long) this->inv_MAT_EB_EB.axis(1);  // int converted to long
       naxes[1] =  (long) this->inv_MAT_EB_EB.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_EB_EB.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_EB_EB.buffer(),&status) ) {
               printf("Error: cannot write inv EBEB matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("EBEB: %d out of %d\n", hdunum, nhdu);
    }

    //SAVE EB 4 blocks inv Matrix
       naxes[0] = (long) this->inv_MAT_EB_4Blocks.axis(1);  // int converted to long
       naxes[1] = (long) this->inv_MAT_EB_4Blocks.axis(2);
    if (fits_create_img(fptr,bitpix, 2, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->inv_MAT_EB_4Blocks.n_elem();
    if(npixels >=0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, this->inv_MAT_EB_4Blocks.buffer(),&status) ) {
               printf("Error: cannot write inv  EB 4Blocks matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if(this->Verbose) {
        res =fits_get_hdu_num(fptr,&hdunum);
        res =fits_get_num_hdus(fptr,&nhdu,&status);
        printf("EB 4 blocks: %d out of %d\n", hdunum, nhdu);
    }

    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: write_master_allPSPEC(char * Name) {
    char *NameFits=fitsname(Name);
    std::string infile=std::string(NameFits);

    //Get number of spectra
    long nspec = 0;
    long szspec[6];
    szspec[0] = (long) Master_TT_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;
    szspec[1] = (long) Master_EE_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;
    szspec[2] = (long) Master_BB_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;
    szspec[3] = (long) Master_TE_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;
    szspec[4] = (long) Master_TB_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;
    szspec[5] = (long) Master_EB_Powspec.size();  // int converted to long
    if(szspec>(void *)0) ++nspec;

    PowSpec HPspec= PowSpec(nspec, lLmax);
    arr<double> arr_TT_temp=Master_TT_Powspec;
    if(nspec==1) HPspec.Set(arr_TT_temp);
    else if(nspec>1){
        arr<double> arr_EE_temp=Master_EE_Powspec;
        arr<double> arr_BB_temp=Master_BB_Powspec;
        arr<double> arr_TE_temp=Master_TE_Powspec;
        if(nspec==4){
            HPspec.Set(arr_TT_temp, arr_EE_temp, arr_BB_temp, arr_TE_temp);
        } else {
            arr<double> arr_TB_temp=Master_TB_Powspec;
            arr<double> arr_EB_temp=Master_EB_Powspec;
            HPspec.Set(arr_TT_temp, arr_EE_temp, arr_BB_temp, arr_TE_temp,
                                                    arr_TB_temp, arr_EB_temp);
        }
    }
    fitshandle out;
    std::ifstream intest(NameFits);
    if (intest.good()) out.delete_file(infile);
    write_powspec_to_fits (infile,HPspec,nspec);
    free(NameFits);
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: write_master_PSPEC(char * Name) {
    fitsfile *fptr;   // pointer to the FITS file
      long naxes[1], firstpixel=1, npixels;
       int status = 0;
      int bitpix = -64;
      char FitsName[MAXCHAR];
    char *NameFits=fitsname(Name);
    strcpy(FitsName, NameFits);
    free(NameFits);
    if ( fits_create_file(&fptr, FitsName, &status) ) {
        printf("Error: cannot open file %s %d ", FitsName, status);
        exit(-1);
     }
     remove(FitsName);

    //SAVE Master_TT_Powspec
       naxes[0] = (long) this->Master_TT_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_TT_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_TT_Powspec.begin(),&status) ) {
               printf("Error: cannot write Master_TT_Powspec %s %d ", FitsName, status);
               exit(-1);
           }
    }

    //SAVE Master_EE_Powspec
       naxes[0] = (long) this->Master_EE_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_EE_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_EE_Powspec.begin(),&status) ) {
               printf("Error: cannot write Master_EE_Powspec %s %d ", FitsName, status);
               exit(-1);
           }
    }
       //SAVE Master_BB_Powspec
       naxes[0] = (long) this->Master_BB_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_BB_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_BB_Powspec.begin(),&status) ) {
               printf("Error: cannot write Master_BB_Powspec %s %d ", FitsName, status);
               exit(-1);
           }
    }

    //SAVE Master_TE_Powspec
       naxes[0] = (long) this->Master_TE_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_TE_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_TE_Powspec.begin(),&status) ) {
               printf("Error: cannot write Master_TE_Powspec %s %d ", FitsName, status);
               exit(-1);
           }
    }
    //SAVE Master_TB_Powspec
       naxes[0] = (long) this->Master_TB_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_TB_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_TB_Powspec.begin(),&status) ) {
               printf("Error: cannot write Master_TE_Powspec %s %d ", FitsName, status);
               exit(-1);
           }
    }

    //SAVE Master_EB_Powspec
       naxes[0] = (long) this->Master_EB_Powspec.size();  // int converted to long
    if (fits_create_img(fptr,bitpix, 1, naxes, &status) ) {
        printf("Error: cannot write header %s %d", FitsName, status);
        exit(-1);
       }
       npixels = (long) this->Master_EB_Powspec.size();
    if(npixels >0) {
        if ( fits_write_img( fptr, TDOUBLE, firstpixel, npixels, Master_EB_Powspec.begin(),&status) ) {
               printf("Error: cannot write inv EBEB matrix %s %d ", FitsName, status);
               exit(-1);
           }
    }
    if ( fits_close_file(fptr, &status) ) {
       printf("Error: cannot close %s %d ", FitsName, status);
       exit(-1);
    }
}
