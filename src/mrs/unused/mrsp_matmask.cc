/******************************************************************************
**                   Copyright (C) 2013 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Florent Sureau
**
**    Date:  Oct 13
**
**    File:  mrsp_matmask.cc
**
*******************************************************************************
**
**    DESCRIPTION Polarized MASTER deconvolution to obtain power spectra for
**		polarized masked maps on the sphere
**    -----------
**
**    Usage: mrsp_matmask options TQU/TEB PolaMap_File_input MaskT_input spec_file_output
**
******************************************************************************
**
**    HISTORY
**    -----------
**
**    Core is computation of standart coupling matrices
**	  from M. Tristram routines (i.e. wrapper to wigner-3j fortran routines,
**	and calling routine to compute the matrices commented below)
**    This code was rewritten for better runtime performance and
**	  added input/outputs routines, own inversion, open-mp parallelization
**
******************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// #include "Array.h"
#include "Array.h"
#include "HealpixClass.h"
#include "PolaHealpixClass.h"
#include "mrsp_lib_matmask.cc"
#include "fs_omp_routines.h"

extern"C"  {

void wig3j_( double *L2, double *L3, double *M2, double *M3,
	     double *L1MIN, double *L1MAX, double *THRCOF, int *NDIM, int *IER);
void wig3j_c( long l2, long l3, long m2, long m3, double *wigner);
void make_mll_pol( long lmax, double *well, double *mll);


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


    //Wrapper to FORTRAN routine wig3j_f.f
    //Be careful, for vanishing wigner coeff, FORTRAN routine does not return 0.
    void wig3j_c( long l2, long l3, long m2, long m3, double *wigner)
    {
        double *THRCOF;
        double l1min, l1max;
        double dl2, dl3, dm2, dm3;
        int ndim, ier=0, l;

        l1min = MAX( abs(l2-l3), abs(m2+m3));
        l1max = l2 + l3;
        ndim = (int)(l1max-l1min+1);

        dl2 = (double)l2;
        dl3 = (double)l3;
        dm2 = (double)m2;
        dm3 = (double)m3;

        // wigner=0 for l<l1min
        for( l=0; l<l1min; l++) wigner[l] = 0.;

        THRCOF = wigner + (int)l1min;    // wig3j start at max( abs(l2-l3), abs(m2+m3))
        if( l2<abs(m2) || l3<abs(m3) || (l1max-l1min)<0 ) for( l=0; l<ndim; l++) THRCOF[l] = 0.;
        else wig3j_( &dl2, &dl3, &dm2, &dm3, &l1min, &l1max, THRCOF, &ndim, &ier);

        if( ier) {
            for( l=0; l<ndim; l++) THRCOF[l] = 0.;
            printf( "err=%d  l2=%d, l3=%d, m2=%d, m3=%d : ", ier, l2, l3, m2, m3);
            switch( ier)
            {
                case 1 : printf( "Either L2.LT.ABS(M2) or L3.LT.ABS(M3)\n"); break;
                case 2 : printf( "Either L2+ABS(M2) or L3+ABS(M3) non-integer\n"); break;
                case 3 : printf( "L1MAX-L1MIN not an integer (l1min=%d, l1max=%d)\n", (int)l1min, (int)l1max); break;
                case 4 : printf( "L1MAX less than L1MIN (l1min=%d, l1max=%d)\n", (int)l1min, (int)l1max); break;
                case 5 : printf( "NDIM less than L1MAX-L1MIN+1 (ndim=%d)\n", ndim); break;
            }
            fflush(stdout);
        }

    }



/*void make_mll_pol( long lmax, double *well, double *mll)
{
  double *mll_TT_TT, *mll_EE_EE, *mll_EE_BB, *mll_TE_TE, *mll_EB_EB;
  double *well_TT, *well_TP, *well_PP;
  double sum_TT, sum_TE, sum_EE_EE, sum_EE_BB, sum_EB;
  double *wigner0, *wigner2;
  long ndim, maxl3;
  long n=lmax+1;
  long well_lmax = lmax;   // well supposed to be lmax+1 size long also

  // split Mll into components
  mll_TT_TT = &mll[ 0*n*n];
  mll_EE_EE = &mll[ 1*n*n];
  mll_EE_BB = &mll[ 2*n*n];
  mll_TE_TE = &mll[ 3*n*n];
  mll_EB_EB = &mll[ 4*n*n];

  // split well into components
  well_TT = &well[ 0*(well_lmax+1)];
  well_TP = &well[ 1*(well_lmax+1)];
  well_PP = &well[ 2*(well_lmax+1)];

  // loop over the matrice elements
  for( long l1=0; l1<=lmax; l1++) {
    for( long l2=0; l2<=lmax; l2++) {

      // initialization
      sum_TT    = 0.;
      sum_TE    = 0.;
      sum_EE_EE = 0.;
      sum_EE_BB = 0.;
      sum_EB    = 0.;

      // alloc wigners
      ndim = l2+l1-abs(l1-l2)+1;
      wigner0 = (double *) malloc( ndim * sizeof(double));
      wigner2 = (double *) malloc( ndim * sizeof(double));

      // compute wigners
      wig3j_c( l1, l2, (long) 0,  (long) 0, wigner0);
      wig3j_c( l1, l2,  (long)-2,  (long)2, wigner2);

      // loop on l3
      maxl3 = (l1+l2 < well_lmax) ? l1+l2 : well_lmax;
      for( long l3=abs(l1-l2); l3<=maxl3; l3++) {

	if( (l1+l2+l3)%2 == 0) sum_TT    += well_TT[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner0[l3];
	if( (l1+l2+l3)%2 == 0) sum_TE    += well_TP[l3] * (double)(2.*l3+1.) * wigner0[l3] * wigner2[l3];
	if( (l1+l2+l3)%2 == 0) sum_EE_EE += well_PP[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
	if( (l1+l2+l3)%2 != 0) sum_EE_BB += well_PP[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
	sum_EB += well_PP[l3] * (double)(2.*l3+1.) * wigner2[l3] * wigner2[l3];
      }

      mll_TT_TT[ l1*n+l2] = (2.*double(l2)+1.)/(4.*M_PI) * sum_TT;
      mll_EE_EE[ l1*n+l2] = (2.*double(l2)+1.)/(4.*M_PI) * sum_EE_EE;
      mll_EE_BB[ l1*n+l2] = (2.*double(l2)+1.)/(4.*M_PI) * sum_EE_BB;
      mll_TE_TE[ l1*n+l2] = (2.*double(l2)+1.)/(4.*M_PI) * sum_TE;
      mll_EB_EB[ l1*n+l2] = (2.*double(l2)+1.)/(4.*M_PI) * sum_EB;

      free( wigner0);
      free( wigner2);

    } //end loop l2
  } //end loop l1
}*/


}



/*********************************************************************/
char flag_tqu_teb[4];
char Name_file_Map[1024];
char Name_file_EB[1024];
char Name_file_MaskT[1024];
char Name_file_MaskP[1024];
char Name_file_masked_map_spec[1024];
char Name_file_mask_spec[1024];
char Name_file_coupling_mat[1024];
char Name_file_inv_coupling_mat[1024];
char Name_file_spec_radii[1024];
char Name_file_MASTER_spec[1024];
char Name_file_noise[1024];
char Name_file_spec_noise[1024];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool NormALM = false;
bool PolaFastALM = true;

bool Flag_MaskP =false;
bool Flag_Save_EB =false;
bool Flag_Work_EB =false;
bool Flag_masked_map_spec = false;
bool Flag_mask_spec = false;
bool Flag_coupling_mat = false;
bool Flag_fastinv = false;
bool Flag_SVD = false;
bool Flag_inv_coupling_mat = false;
bool Flag_spec_radii = false;
bool Flag_onlyT=false;
bool Flag_positivity=false;
bool Verbose=false;
bool Flag_timer=false;
bool Zero_first2l=false;
bool FlagNoise=false;
bool FlagSpecNoise=false;
bool FlagUncorrNoise=false;
bool WriteIntFile = false;
bool set_threads=false;

float ZeroPadding=0;
int Lmax = 0;
int Nit_inv = 0;

#ifdef _OPENMP
 	extern int inner_loop_threads;
  	extern int outer_loop_threads;
#endif

/***************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options TQU/TEB PolaMap_File_input MaskT_input spec_file_output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-a Nout:Ninn\n");
    fprintf(OUTMAN, "             Specify number of threads used (outer loops and inner loops).\n");

    fprintf(OUTMAN, "         [-e] \n");
    fprintf(OUTMAN, "             Maps are masked in TEB. Default is maps are masked in TQU.\n");

    fprintf(OUTMAN, "         [-f] \n");
    fprintf(OUTMAN, "             Fast inversion. Default is false (faster upper bound for iterative approach, or fast block inverse\n");

    fprintf(OUTMAN, "         [-i] \n");
    fprintf(OUTMAN, "             Invert Coupling Matrix with SVD. Default is no SVD (iterative approach).\n");

    fprintf(OUTMAN, "         [-k Nits]\n");
    fprintf(OUTMAN, "             Number of iterations for iterative approach (default 40).\n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is 3*nside \n");

    fprintf(OUTMAN, "         [-m Name_File_Polarized Mask.\n");
    fprintf(OUTMAN, "             Mask for Polarization (default : same as temperature).\n");

    fprintf(OUTMAN, "         [-n]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");

    fprintf(OUTMAN, "         [-o]\n");
    fprintf(OUTMAN, "             Timer on.\n");

    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Positivity. Default is no.\n");

    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");

    fprintf(OUTMAN, "         [-t]\n");
    fprintf(OUTMAN, "             Only Temperature. Default is all.\n");

    fprintf(OUTMAN, "         [-u]\n");
    fprintf(OUTMAN, "             Uncorrelated Noise in Pixel space (i.e. noise only impacts TT, EE and BB).\n");

    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbosity.\n");

    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Write Intermediary Files.\n");

    fprintf(OUTMAN, "         [-z]\n");
    fprintf(OUTMAN, "             zero monopole and dipole in solution (temperature).\n");

    fprintf(OUTMAN, "         [-E Name_File_EB]\n");
    fprintf(OUTMAN, "         	  FileName to save EB map (if input map is TQU).\n");

    fprintf(OUTMAN, "         [-I Name_file_inv_coupling_matrix]\n");
    fprintf(OUTMAN, "             FileName where the inverse coupling matrix is stored.\n");

    fprintf(OUTMAN, "         [-K Name_file_coupling_matrix]\n");
    fprintf(OUTMAN, "             FileName where the coupling matrix is stored.\n");

    fprintf(OUTMAN, "         [-M Name_file_mask_spectra]\n");
    fprintf(OUTMAN, "             FileName where the power spectra of the masks are stored.\n");

    fprintf(OUTMAN, "         [-N Name_file_noise]\n");
    fprintf(OUTMAN, "             FileName where the noise maps are stored.\n");

    fprintf(OUTMAN, "         [-O Name_file_mask_spectra]\n");
    fprintf(OUTMAN, "             FileName where the power spectra of the noise maps are stored.\n");

    fprintf(OUTMAN, "         [-P Name_file_masked_map_spectra]\n");
    fprintf(OUTMAN, "             FileName where the power spectra of the masked maps are saved.\n");

    fprintf(OUTMAN, "         [-R Name_file_spec_radii]\n");
    fprintf(OUTMAN, "             Name of file where spectral radii of the submatrices are saved (iterative approach).\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    /* get options */
    while( (c = GetOpt( argc,argv,(char *) "a:efik:l:m:nopstuvwzE:I:K:M:N:O:P:R:" )) != -1 ) {
	switch (c) {
		case 'a':
			#ifdef _OPENMP
	        	if (sscanf(OptArg,"%d:%d",&outer_loop_threads,&inner_loop_threads) < 1)
				{
		   			fprintf(OUTMAN, "Error: syntaxe is outer_loop_threads:inner_loop_threads ... \n");
		   			exit(-1);
				}
	        	set_threads=true;
			#endif
	            break;
		case 'e':
			Flag_Work_EB = true;
			break;
		case 'f':
			Flag_fastinv = true;
			break;
		case 'i':
			Flag_SVD = true;
			 break;
		case 'k':
			if( sscanf(OptArg,"%d",&Nit_inv) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Nit_inv  value: %s\n", OptArg);
		            exit(-1);
			}
			break;
		case 'l':
			if( sscanf(OptArg,"%d",&Lmax) != 1 ) {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
			}
			break;
		case 'm' :
			if( sscanf(OptArg,"%s",Name_file_MaskP) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for Mask for Polarization: %s\n", OptArg);
		            exit(-1);
			}
			Flag_MaskP=true;
			break;
		case 'n':
           		NormALM = true;
 		        break;
		case 'o':
           		Flag_timer = true;
 		        break;
		case 'p':
           		Flag_positivity = true;
 		        break;
 		case 's':
	   		PolaFastALM = false;
 		        break;
 		case 't':
	   		Flag_onlyT = true;
 		        break;
 		case 'u':
	   		FlagUncorrNoise = true;
 		        break;
 		case 'v':
	   		Verbose = true;
 		        break;
 		case 'w':
	   		WriteIntFile = true;
 		     break;

 		case 'z':
	   		Zero_first2l = true;
 		    break;
		case 'E':
			if( sscanf(OptArg,"%s",Name_file_EB) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for EBmap: %s\n", OptArg);
		            exit(-1);
			}
			Flag_Save_EB=true;
			break;
		case 'I':
  	              if( sscanf(OptArg,"%s",Name_file_inv_coupling_mat) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for inverse Coupling Matrix: %s\n", OptArg);
		            exit(-1);
			}
			Flag_inv_coupling_mat=true;
			break;
		case 'K':
  	              if( sscanf(OptArg,"%s",Name_file_coupling_mat) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for Coupling Matrix: %s\n", OptArg);
		            exit(-1);
			}
			Flag_coupling_mat=true;
			break;
		case 'M':
  	              if( sscanf(OptArg,"%s",Name_file_mask_spec) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for Mask for Polarization: %s\n", OptArg);
		            exit(-1);
			}
			Flag_mask_spec=true;
			break;
		case 'N':
  	              if( sscanf(OptArg,"%s",Name_file_noise) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for Mask for Polarization: %s\n", OptArg);
		            exit(-1);
			}
			FlagNoise =true;
			break;
		case 'O':
  	              if( sscanf(OptArg,"%s",Name_file_spec_noise) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for Mask for Polarization: %s\n", OptArg);
		            exit(-1);
			}
			FlagSpecNoise =true;
			break;
		case 'P':
			if( sscanf(OptArg,"%s",Name_file_masked_map_spec) != 1 ) {
		            fprintf(OUTMAN, "Error: bad Filename for EBmap: %s\n", OptArg);
		            exit(-1);
			}
			Flag_masked_map_spec=true;
			break;
 	     	default:
 	     		usage(argv);
 	     		break;
 		case 'R':
  	        if( sscanf(OptArg,"%s",Name_file_spec_radii) != 1 ) {
		    	fprintf(OUTMAN, "Error: bad Filename for Spectral Radii: %s\n", OptArg);
		    	exit(-1);
			}
			Flag_spec_radii=true;
			break;
		}
	}


	// get optional input file names from trailing parameters and open files
	if( OptInd < argc ) {
		strcpy( flag_tqu_teb, argv[OptInd] );
		OptInd++;
	} else usage(argv);

	if( OptInd < argc ) {
	    	strcpy( Name_file_Map, argv[OptInd] );
    		OptInd++;
	} else usage(argv);

	if( OptInd < argc ) {
	    	strcpy( Name_file_MaskT, argv[OptInd] );
    		OptInd++;
	} else usage(argv);

	if( OptInd < argc ) {
    		strcpy( Name_file_MASTER_spec, argv[OptInd] );
    		OptInd++;
	} else usage(argv);

	/* make sure there are not too many parameters */
	if( OptInd < argc ) {
		fprintf( OUTMAN, "Too many parameters: %s ...\n", argv[OptInd] );
		exit(-1);
	}
}


int main( int argc, char *argv[], char** envp )
{
extern int GetOpt(int argc, char **argv, char *opts);

	//char** env;
 	/*for (env = envp; *env != 0; env++) {
    	char* thisEnv = *env;
    	printf("%s\n", thisEnv);
  	}*/
	printf("Inside Main\n");
	sinit(argc, argv);

	#ifdef _OPENMP
		if(! set_threads) auto_set_inner_outer_threads(0,0);
		else {
			if((outer_loop_threads<=0) || (inner_loop_threads<=0)) {
				printf("Not Valid Thread numbers: %d, %d\n Automatic Initialization instead\n", outer_loop_threads, inner_loop_threads);
				auto_set_inner_outer_threads(0,0);
			} else  set_nb_threads(outer_loop_threads,inner_loop_threads);
		}
		openmp_params();
	#endif

	MasterPola<double> MasterP;
	MasterP.set_Verbose(Verbose);
	MasterP.set_Timer(Flag_timer);
	MasterP.set_WORK_EB(Flag_Work_EB);

	MasterP.set_First2L_tozero(Zero_first2l);
	MasterP.set_Uncorr_Noise(FlagUncorrNoise);
	if( (strcmp( flag_tqu_teb, "TQU" ) != 0) && (strcmp( flag_tqu_teb, "TEB" ) != 0) ){
		fprintf(OUTMAN, "Wrong flag parameter in %s call \n\n", argv[0]);
		fprintf(OUTMAN, "Usage: %s TQU PolaMap_File_input spec_file_output \n\n OR \n\n %s TEB PolaMap_File_input spec_file_output \n\n", argv[0], argv[0]);
		exit(-1);
	}

	//get root name for intermediary files
	int len_substr;
	char *pos_substr, prefix[1024],tempname[1024];
	pos_substr=strrchr(Name_file_MASTER_spec,'.');
	if(pos_substr==NULL) pos_substr=strrchr(Name_file_MASTER_spec,'\0');
	len_substr=pos_substr-Name_file_MASTER_spec;
	strncpy(prefix,Name_file_MASTER_spec,len_substr);
	prefix[len_substr]='\0';
	if(Verbose) printf("prefix=%s\n",prefix);


	if( Lmax > 0 ) MasterP.set_Lmax(Lmax);
	if(! Flag_masked_map_spec) { //Need to load maps
		if(Flag_onlyT) MasterP.read_HmapT(Name_file_Map);
		else {
			if( strcmp( flag_tqu_teb, "TQU" ) == 0 ) {
				if(Verbose) printf("Load TQU maps\n");
				MasterP.read_PolaHmap(Name_file_Map,false,false,false);
				if(Flag_Save_EB) MasterP.write_PolaHmap_EB(Name_file_EB,PolaFastALM);//Save EB map if QU input
			} else {
				if(Verbose) printf("Load TEB maps\n");
				MasterP.read_PolaHmap(Name_file_Map,true,false,false);
			}
		}
		if(Verbose) printf("Loaded Input Map\n");
	 } else {
		MasterP.read_TEB_Pspec_from_file(Name_file_masked_map_spec);
		if(Verbose) printf("Loaded Input Map Power Spectra\n");
	}

	if((!FlagSpecNoise)&&(FlagNoise)) { //Need to load Noise maps
		if(Flag_onlyT) MasterP.read_NoiseHmapT(Name_file_noise);
		else {
			if( strcmp( flag_tqu_teb, "TQU" ) == 0 ) {
				if(Verbose) printf("Load TQU Noise maps\n");
				MasterP.read_NoisePolaHmap(Name_file_noise,false,false,false);
			} else {
				if(Verbose) printf("Load TEB Noise maps\n");
				MasterP.read_NoisePolaHmap(Name_file_noise,true,false,false);
			}
		}
		if(Verbose) printf("Loaded Noise Map\n");
	 } else if (FlagSpecNoise) {
		MasterP.read_NoiseTEB_Pspec_from_file(Name_file_spec_noise);
		if(Verbose) printf("Loaded Noise Map Power Spectra\n");
	}

	if((!Flag_mask_spec)||(! Flag_masked_map_spec)||((!FlagSpecNoise)&&(FlagNoise))) {//Need to load Masks
		if(Verbose) printf("SET MASK DATA\n");
		MasterP.read_MaskT(Name_file_MaskT);
		if(Verbose) printf("Loaded MaskT\n");
		if(Flag_MaskP) {
			MasterP.read_MaskP(Name_file_MaskP);
			if(Verbose) printf("Loaded MaskP\n");
		}
	}

	if((!Flag_coupling_mat)&&(!Flag_inv_coupling_mat)) {
		if(! Flag_mask_spec) { //load mask(s) and compute the mask power spectra
			if(Verbose) printf("Masks power spectra needs computing\n");
			MasterP.get_Pspec_from_masks(NormALM,PolaFastALM);
			if(Verbose) printf("Masks power spectra computed\n");
			sprintf(tempname,"%s_mask_pspec.fits",prefix);
			if(WriteIntFile) MasterP.write_Mask_Pspec(tempname);
		} else {//load mask power spectra
			MasterP.read_Mask_Pspec_from_file(Name_file_mask_spec);
			if(Verbose) printf("Loaded Input Mask Power Spectra\n");
		}
	}
	if(! Flag_masked_map_spec) {//Compute also the masked map power spectra
		if(Verbose) printf("Computing Masked Map Power Spectra\n");
		MasterP.get_TEB_Pspec_from_map(NormALM,PolaFastALM);
		if(Verbose) printf("Masked Map Power Spectra Computed\n");
		sprintf(tempname,"%s_maskedmap_pspec.fits",prefix);
		if(WriteIntFile) MasterP.write_MaskTEB_Pspec(tempname);
	}
	if((!FlagSpecNoise)&&(FlagNoise)) { //Need to Compute Masked Noise maps
		if(Verbose) printf("Computing Masked Noise Map Power Spectra\n");
		MasterP.get_TEB_Pspec_from_Noisemap(NormALM,PolaFastALM);
		if(Verbose) printf("Masked Noise Map Power Spectra Computed\n");
		sprintf(tempname,"%s_maskednoisemap_pspec.fits",prefix);
		if(WriteIntFile) MasterP.write_MaskNoiseTEB_Pspec(tempname);
	}

	if(! Flag_inv_coupling_mat) {//Need to compute inverse coupling matrix and apply to the data
		if( !Flag_coupling_mat) {//First compute coupling matrix
			if(Verbose) printf("Computing Coupling Matrix\n");
			//MasterP.make_mll_blocks_c();
			MasterP.make_mll_blocks_c_fast();
			if(Verbose) printf("Coupling Matrix Computed\n");
			sprintf(tempname,"%s_coupling_matrices.fits",prefix);
			if(WriteIntFile) MasterP.write_coupling_mat(tempname);
		} else {
			MasterP.read_coupling_mat(Name_file_coupling_mat);
			if(Verbose) printf("Loaded Coupling Matrix\n");
		}

		unsigned short xsubi[3]={1980,1947,1978};
		bool Flag_iter=!Flag_SVD;

		if((Flag_iter)&&(Flag_spec_radii)) {
			MasterP.read_spec_radii(Name_file_spec_radii);
			if(Verbose) printf("Loaded Spectral Radii\n");
		}
		if((Verbose)&&(Flag_iter)) printf("ITERATIVE MASTER of PowerSpectra, iter=%d \n", Nit_inv);
		if((Verbose)&&(!Flag_iter)) printf("SVD MASTER of PowerSpectra\n");

		if(Flag_onlyT) MasterP.get_MASTER_pspec(xsubi,PSPEC_TT,Flag_iter,Nit_inv,Flag_positivity,Flag_fastinv);
		else MasterP.get_MASTER_pspec(xsubi,PSPEC_ALL,Flag_iter,Nit_inv,Flag_positivity,Flag_fastinv);

		if(Flag_iter) {
			if(Flag_fastinv) sprintf(tempname,"%s_mask_spec_radii_fast.fits",prefix);
			else sprintf(tempname,"%s_mask_spec_radii.fits",prefix);
			if(WriteIntFile) MasterP.write_spec_radii(tempname);
		} else {
			if(Flag_fastinv) sprintf(tempname,"%s_mask_invmat_fast.fits",prefix);
			else sprintf(tempname,"%s_mask_invmat.fits",prefix);
			if(WriteIntFile) MasterP.write_inv_mat(tempname);
		}
	} else {
		unsigned short xsubi[3]={1980,1947,1978};
		MasterP.read_inv_mat(Name_file_inv_coupling_mat);
		if(Flag_onlyT) MasterP.get_MASTER_pspec(xsubi,PSPEC_TT,false,Nit_inv,Flag_positivity,Flag_fastinv);
		else MasterP.get_MASTER_pspec(xsubi,PSPEC_ALL,false,Nit_inv,Flag_positivity,Flag_fastinv);
	}

	//Finally Save Spectra
	//MasterP.write_master_PSPEC(Name_file_MASTER_spec);
	MasterP.write_master_allPSPEC(Name_file_MASTER_spec);

	return 0;
}



/*********************************************************************/
