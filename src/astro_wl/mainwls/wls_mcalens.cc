/******************************************************************************************
**                   Copyright (C) 2009 by CEA  
*******************************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Olivier Fourt
**
**    Date:  26/06/09
**    
**    File:  mrsp_qu_to_eb.cc
**
********************************************************************************************
**
**    DESCRIPTION  convert Healpix polarized map from Q and U fields to E and B fields
**    ----------- 
**                 
*******************************************************************************************/

#include "HealpixClass.h"
#include "WLS_MassMapping.h"
#include "MRS_Sparse.h"
#include "Healpix_PowSpec.h"

// #include "map_alm_qu_eb.h"

char Name_file_Map1[512];
char Name_file_Map2[512];
char Name_file_Map3[512];
char Name_file_Map4[512];
char SignalPowSpec_FileName [512];
char Name_file_Noise_spec[512];
char Name_file_NoiseCovMat[512];
char Name_file_Mask[512];

Bool OptS_PS = False;
Bool OptN_CovMat = False;
Bool OptMask = False;
Bool MethodOpt=False;

bool input_spectrum = false;
bool NormALM = false;

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool PolaFastALM = true;
Bool Verbose=False;

float  ZeroPadding=0.0;

float SigmaNoise=0.;
int NbrScale=0;
int Lmax=0;
int ALM_iter=0;
int Niter=50;
float NSigma=3.;
int FirstDetectScale=0;
bool SparsePositivityConstraint=true;

#define NBR_MASS_MAPPING_METHOD 4
enum type_massmapping {KAISER_SQUIRES, WIENER, SPARSE, MCAlens};
type_massmapping MassMethod = KAISER_SQUIRES;

const char * StringMassMethod (type_massmapping type)
{
    switch (type)
    {
        case  KAISER_SQUIRES:
            return ("Spherical Kaiser-Squires");break;
        case  WIENER:
            return ("Proximal Wiener Filtering");break;
        case  SPARSE:
            return ("Sparse Recovery");break;
        case MCAlens:
            return ("Sparse-Wiener Reconstruction (MCAlens)");break;
        default:
            return ("Undefined method");
            break;
    }
}

inline void all_massmapping_usage(type_massmapping MassMethod)
{
     fprintf(OUTMAN, "        [-m mass_mapping_method]\n");
     for (int i = 0; i < NBR_MASS_MAPPING_METHOD; i++)
     fprintf(OUTMAN, "              %d: %s \n",i+1,
                                             StringMassMethod((type_massmapping)i));
     fprintf(OUTMAN, "             default is %s\n",  StringMassMethod((type_massmapping) MassMethod));
}

/***********************************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Map_G1_File_input Map_G2_File_input Suffixe_File_output  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    all_massmapping_usage(MassMethod);
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");
    fprintf(OUTMAN, "         [-i Number_of_Iterations]\n");
    fprintf(OUTMAN, "             Number of iterations (for iter algorithms). Default is %d.\n", Niter);
    fprintf(OUTMAN, "         [-a N_ITER_ALM]\n");
    fprintf(OUTMAN, "             Number of iterations N_ITER_ALM for iterative spherical harmonics transform.  \n");
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    fprintf(OUTMAN, "         [-g GaussianNoise]\n");
    fprintf(OUTMAN, "             Noise Standard deviation.  \n");
    fprintf(OUTMAN, "         [-s Nsigma]\n");
    fprintf(OUTMAN, "             Nsigma Detection level.  \n");
    fprintf(OUTMAN, "         [-F FirstDetectionScale]\n");
    fprintf(OUTMAN, "             First detection scale.  \n");
    fprintf(OUTMAN, "         [-S Signal_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Input Signal Power Spectrum file name (for Wiener filtering). Default is automatically estimated. \n");
    fprintf(OUTMAN, "         [-N Noise_Covariance_Matrix_FileName]\n");
    fprintf(OUTMAN, "             Input Cov. Mat file name (for Wiener filtering).\n");
    fprintf(OUTMAN, "         [-M Mask_FileName]\n");
    fprintf(OUTMAN, "             Input Mask file name.\n");
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Remove the positivy contraint in the MCAlens sparse component.\n");


  //  fprintf(OUTMAN, "         [-P]\n");
  //  fprintf(OUTMAN, "             Input signal spectrum. Default is no, spectrum estimated from image in spectrum and noise spectrum.\n");
	fprintf(OUTMAN, "         [-v Verbose]\n");
    exit(-1);
}

/************************************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c,m;
    // get options
     
    while( (c = GetOpt( argc,argv,(char *) "pF:i:m:g:M:N:a:S:n:s:i:l:fv" )) != -1 )
    {
		switch (c) 
        {
            case 'p':
                SparsePositivityConstraint=false;
                break;
            case 'F':
                if (sscanf(OptArg,"%d",&FirstDetectScale) != 1)
                {
                    fprintf(OUTMAN, "Error: bad first detection scale: %s\n", OptArg);
                    exit(-1);
                }
                FirstDetectScale --;
                if (FirstDetectScale < 0) FirstDetectScale=0;
                break;
            case 'm':
                 if (sscanf(OptArg,"%d",&m) != 1)
                 {
                     fprintf(OUTMAN, "Error, bad type of mass mapping method: %s\n", OptArg);
                     exit(-1);
                 }
                if ((m > 0) && (m <= NBR_MASS_MAPPING_METHOD))
                     MassMethod = (type_massmapping) (m-1);
                 else
                 {
                     fprintf(OUTMAN, "Error: bad type of mass mapping method: %s\n", OptArg);
                     exit(-1);
                 } 
                 MethodOpt = True;
                 break;
            case 'i':
                if (sscanf(OptArg,"%d",&Niter) != 1)
                {
                    fprintf(OUTMAN, "Error: bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
            case 'g':
                    if (sscanf(OptArg,"%f",&SigmaNoise) != 1)
                    {
                        fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                        exit(-1);
                    }
                    break;
            case  'M':
                     if (sscanf(OptArg,"%s", Name_file_Mask) != 1)
                     {
                       fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                      exit(-1);
                      }
                   OptMask = True;
                   break;
            case  'N':
                     if (sscanf(OptArg,"%s", Name_file_NoiseCovMat) != 1)
                     {
                       fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                      exit(-1);
                      }
                   OptN_CovMat = True;
                   break;
            case 'a':
                     if (sscanf(OptArg,"%d",&ALM_iter) != 1)
                      {
                          fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                           exit(-1);
                           }
                         break;
            case  'S':
                     if (sscanf(OptArg,"%s", SignalPowSpec_FileName) != 1)
                     {
                       fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                      exit(-1);
                      }
                   input_spectrum = true;
                   OptS_PS = True;
                   break;
             case 'n':
                  if (sscanf(OptArg,"%d",&NbrScale) != 1)
                      {
                        fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                        exit(-1);
                        }
                      break;
           case 's':
               if (sscanf(OptArg,"%f",&NSigma) != 1)
               {
                   fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                   exit(-1);
               }
               break;
 		    case 'f':
				OptInd--;
           		PolaFastALM = false;
 		        break;
	        case 'l':
                if (sscanf(OptArg,"%d",&Lmax) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
     	    case 'v': Verbose = True; break;
 	     	default:
 	     		usage(argv);
 	     		break;
 		}
	} 

	//get optional input file names from trailing parameters and open files
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map1, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map2, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map3, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
	// make sure there are not too many parameters
	if( OptInd < argc )
    {
		fprintf( OUTMAN, "Too many parameters: %s ...\n", argv[OptInd] );
		exit(-1);
	}
}
  
/*********************************************************************/

int main( int argc, char *argv[] )
{
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

	sinit(argc, argv);

    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        if (NbrScale > 0)  cout << "# NbrScale = " <<  NbrScale << endl;
        cout << "# File Name in 1 = " << Name_file_Map1 << endl;
        cout << "# File Name in 2 = " << Name_file_Map2 << endl;
        cout << "# Prefix File Name Out = " << Name_file_Map3 << endl;
        if (OptN_CovMat == True)
            cout << "# File Name Cov Mat = " << Name_file_NoiseCovMat << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  SignalPowSpec_FileName << endl;
        cout << "# MassMapping Method = " << StringMassMethod(MassMethod) << endl;
        cout << "# Niter = " << Niter << endl;
        if (SigmaNoise > 0) cout << "# SigmaNoise = " << SigmaNoise << endl;
        if (Lmax > 0) cout << "# Lmax = " << Lmax << endl;
        if (FirstDetectScale > 0) cout << "# FirstDetectScale = " << FirstDetectScale+1 << endl;
        cout << endl;
    }
            
	Hdmap Map_G1, Map_G2, Map_E, Map_B;
    
	Map_G1.read( Name_file_Map1 );
	Map_G2.read( Name_file_Map2 );
	
	if( Map_G1.Scheme() != Map_G2.Scheme() )
	{
		fprintf( OUTMAN, "Maps are not in the same scheme in: %s ...\n", argv[0] );
		exit(-1);
	}
    if (Verbose == True) Map_G1.info((char*) "Input G1");

	if( Map_G1.nside() != Map_G2.nside() )
	{
		fprintf( OUTMAN, "Maps have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
    if (Verbose == True) Map_G2.info((char*) "Input G2");

	int Nside = Map_G1.nside();
    int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    // cout << " Lmax_Tot = " << Lmax_Tot << endl;
    
    Map_E.alloc(Nside);
    Map_B.alloc(Nside);
 	
    PowSpec SigPowSpec(1, Lmax_Tot);
    if (OptS_PS == True)
    {
          mrs_read_powspec(SignalPowSpec_FileName, SigPowSpec);
    }
    else
    {
        if  ((MassMethod == MCAlens) || (MassMethod == WIENER))
        {
           cout << " Error: a theoretical signal power spectrum is required. " << endl;
           exit(-1);
        }
    }
    Hdmap CovMat;
    if (OptN_CovMat == True)
    {
        CovMat.read( Name_file_NoiseCovMat );
    }
    Hdmap Mask;
    if (OptMask == True)
    {
        Mask.read( Name_file_Mask );
    }
    
    /*
    if( input_spectrum == true )
    {
        signal_spec.read( Name_file_Signal_spec );
    }
    noise_spec.read( Name_file_Noise_spec );
    */
    WLS_MassMapping MM;
    WLS_Field Kappa, RecG,KappaSparse;
    MM.Verbose = Verbose;
    if (OptN_CovMat == False)
    {
        if (OptMask == True)
            MM.alloc(Map_G1, Map_G2, SigmaNoise, Mask, NbrScale, Lmax_Tot);
        else
            MM.alloc(Map_G1, Map_G2, SigmaNoise, NbrScale, Lmax_Tot);
    }
    else
    {
        if (OptMask == True)
            MM.alloc(Map_G1, Map_G2, CovMat, Mask,NbrScale, Lmax_Tot);
        else
            MM.alloc(Map_G1, Map_G2, CovMat, NbrScale, Lmax_Tot);
    }
    if (ALM_iter > 0) MM.set_niter(ALM_iter);
    
    int NiterWiener=Niter;
    int NiterSparse=Niter;
    int NiterMCAlens=Niter;

    switch(MassMethod)
    {
        case KAISER_SQUIRES:
             MM.sp_kaiser_squires(MM.GammaData, Kappa);
            break;
        case WIENER:
            MM.iter_wiener(MM.GammaData, SigPowSpec, Kappa, NiterWiener);
            break;
        case SPARSE:
            MM.sparse_reconstruction(MM.GammaData, Kappa, NSigma, NiterSparse, FirstDetectScale);
            break;
        case MCAlens:
            MM.mcalens(MM.GammaData, SigPowSpec, Kappa,  KappaSparse, NSigma, NiterMCAlens, SparsePositivityConstraint, FirstDetectScale);
            break;
        default:
            break;
    }
    
     char FN[512];
     sprintf(FN, "%s_e.fits", Name_file_Map3);
     // Kappa.map_G1.info("Kappa");
     Kappa.map_G1.write(FN);// Write map E
     sprintf(FN, "%s_b.fits", Name_file_Map3);
     Kappa.map_G2.write(FN);// Write map E
    if (MassMethod == MCAlens)
    {
        sprintf(FN, "%s_sparse_e.fits", Name_file_Map3);
        KappaSparse.map_G1.write(FN);// Write map E
        sprintf(FN, "%s_sparse_b.fits", Name_file_Map3);
        KappaSparse.map_G2.write(FN);// Write map B
        sprintf(FN, "%s_active_coef_e.fits", Name_file_Map3);
        fits_write_dblarr(FN, MM.TabActivCoefE);
        sprintf(FN, "%s_active_coef_b.fits", Name_file_Map3);
        fits_write_dblarr(FN, MM.TabActivCoefB);
    }
    if ((MassMethod == MCAlens) || (MassMethod == WIENER))
    {
        sprintf(FN, "%s_filter_wiener_e.fits", Name_file_Map3);
        fits_write_dblarr( FN, MM.CAlm.WienerFilter_E );
        sprintf(FN, "%s_filter_wiener_b.fits", Name_file_Map3);
        fits_write_dblarr( FN, MM.CAlm.WienerFilter_B );
    }

 /*
  int Mmax = Lmax_Tot;

    Hdmap Map;
    //  Map.read(Name_Imag_In);
    // Nside = Map.Nside();
    // int Lmax = mrs_get_lmax (Lmax,  Nside, 0.0);
    // WT.wt_alloc(Nside, NbrScale, Lmax, DEF_MRS_ORDERING);
    // WT.hard_thresholding( Map, NSigma, SigmaNoise, UseMad, KillLastScale);

	Alm<xcomplex<double> > ALM_E(Lmax_Tot,Mmax), ALM_B(Lmax_Tot,Mmax);
	
	// map2alm_pol_iter_QU( Map_G1, Map_G2, ALM_E, ALM_B, 0, weight_T );
    Map_G1.info("Map_G1:");
    Map_G2.info("Map_G2:");
    map2alm_spin(Map_G1,Map_G2,ALM_E,ALM_B,2,weight_T,false);
    alm2map( ALM_E, Map_E );
    alm2map( ALM_B, Map_B );

    // WLS_MassMapping G;
    // ShearAlm AlmG;
    // G.alloc(Map_G1,Map_G2,SigmaNoise);
    
	Map_E.write( Name_file_Map3 );// Write map E
	Map_B.write( Name_file_Map4 );// Write map B
	
    Map_E.info("Map_True");

    Map_E -= Kappa.map_G1;
    Map_E.info("Resi mode E:");
*/
	exit(0);
}

