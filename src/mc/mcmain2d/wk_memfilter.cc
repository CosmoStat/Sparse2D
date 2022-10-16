/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.4
**
**    Author: 07/01/99
**
**    Date:  99/07/05
**    
**    File:  mc_mwfilter.cc
**
*******************************************************************************/

#include "wk_memfilter.h"
#include "MW_mcFilter.h"
#include "MR_Filter.h"
#include "MR_Edge.h"

#define MAX_ITER 100
#define DEF_SIGMA_PCA -1
#define DEF_SIGMA_MR2D 3

#define MAX_IMAG_VALUE 255

#define DEF_CORREL_COMP E_LOCAL
#define DEF_CORREL_NOISE E_THRESHOLD
#define DEF_CONV_PARAM DEF_MEM_FILT_CONGER
#define DEF_REGUL_PARAM 1
#define MAX_NBR_TYPE_OPT 3
#define DEF_SIZE_BLK DEFAULT_SIZE_BLOCK_SIG
#define DEF_ITER_SIGCLP 1
#define DEF_ITER_ENTROP DEF_MEM_FILT_MAX_ITER

extern int  OptInd;
extern char *OptArg;

/***********************************************************************************/
char* getTypeOpt (int ind) {
   switch (ind) {
      case 1 : return ("fixed user Alpha value");break;
      case 2 : return ("Estimate the optimal Alpha");break;
      case 3 : return ("Estimate one  Alpha value per band");break;
   default : return ("Error: bad optim type");break;
   }
}


/***********************************************************************************/
void mr_getmodel_from_edge(MultiResol & MR_Data, MultiResol &MR_Edge,
                           MRNoiseModel & ModelData) {
/* Find a multiresolution model for edges:
    -- create an SNR edge image: ImaEdge
    -- threshold ImaEdge if ImaEdge < NSigma
    -- kill isolated pixel in ImaEdge
    -- in each band, 
        . threshold wavelet coefficient if there is no edge
        .average the wavelet coefficient
         with the two other values in the edge direction
   !! This routine works only if the a-trou algorithm is choosen
*/

   int Nbr_Band =  MR_Data.nbr_band()-1;
   int i,j,Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   Ifloat ImaEdge(Nl,Nc,"ImaEdge");
   Iint ImaAngle(Nl,Nc,"ImaAngle");

   mr_get_edge( MR_Data,  MR_Edge,  ImaEdge, ImaAngle,  ModelData);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
   {
      ImaEdge(i,j) /= ModelData.NSigma[1];
      if (ImaEdge(i,j) > 1) ImaEdge(i,j) = 1.;
      else ImaEdge(i,j) = 0.;
   }
   for (i=0;i < Nl; i++)
   for (j=0;j < Nc; j++) 
       if (isolated_pixel(ImaEdge,i,j,I_MIRROR) == True) ImaEdge(i,j) = 0.;  

   for (int b = 0; b < Nbr_Band; b++) 
   {
      int Nlb = MR_Data.size_band_nl(b);
      int Ncb = MR_Data.size_band_nc(b);
      int Step = (MR_Data.band_to_scale(b) == b) ? b : 0;
  
      for (i = 0; i < Nlb; i++)
      for (j = 0; j < Ncb; j++)
      {
         if (ImaEdge(i,j) > FLOAT_EPSILON)
               MR_Edge(b,i,j) = val_contour_min(MR_Data.band(b), ImaEdge, ImaAngle,
                                              i,j, FLOAT_EPSILON, Step);
         else  MR_Edge(b,i,j) = 0.;
      }
   }
}

int main(int argc, char *argv[]) 
{
   extern softinfo Soft;
   Soft.mr3();
   lm_check(LIC_MR3);

  // licence manager
  //lm2();

  // init local var
  mw_Param o_Param;

  // Read param IN
  o_Param (argc, argv);
  
  // init Result and Iter classes
  mw_Result o_Result (&o_Param);   
  mw_Iter o_Iter (&o_Result); 
     
  for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {} 
      
  // write result
  o_Result.res_Write ();
        
  exit(0);
}





/******************************************************************************
classe mw_Param
******************************************************************************/
 
mw_Param::mw_Param () {

  //default init of read var
  e_Transform =           DEFAULT_TRANSFORM;    
  e_Norm =                NORM_L1; 
  e_SBFilter =            F_MALLAT_7_9;
  po_FAS =                (FilterAnaSynt*)NULL;     
  i_NbPlan =              DEFAULT_NBR_SCALE;
  i_NbIter =              1;
  e_WithNoiseModel =      True;         // always use noise model, def gaussian
  e_UseReconsAdjoint =    False;
  e_TypeNoise =           DEFAULT_STAT_NOISE;
  e_PositivImag =         False;
  e_MaxImag =             False;  
  e_NormCorrelMatrix =    True;
  e_CorrelComp =          DEF_CORREL_COMP;  
  po_ImportCorMat =       NULL;
  e_CorrelNoise =         DEF_CORREL_NOISE;  
  f_NoiseIma =            0.;
  f_Gain =                0.;
  f_SigmaGauss =          0.;
  f_MeanGauss =           0.;
  f_NSigmaMr2d =          DEF_SIGMA_MR2D; 
  f_NSigmaPCA =           DEF_SIGMA_PCA;
  e_Verbose =             False;
  e_WriteTransf =         False;
  e_WriteEigenOnDisk =    False;
  e_WriteSupportMr2d =    False;
  e_WriteSupportPCA =     False;
  e_TraceParamIn =        False;
   
  f_CvgParam =            DEF_CONV_PARAM;
  f_RegulVal =            1.;
  i_TypeOpt =             DEF_MEM_ALPHA_CST;  
  e_DataEdge =            False;
  e_DataSNR =             False;
  e_UseRMSMap =           False;
  i_SizeBlock =           DEF_SIZE_BLK;
  i_NiterClip =           DEF_ITER_SIGCLP;  
  i_NbIterEntrop =        DEF_ITER_ENTROP; 
  e_NscaleOpt =           False;
  e_TransfOpt =           False;
#ifdef LARGE_BUFF
  e_OptZ =                False;
  i_VMSSize =             -1; 
#endif                
  
  for (int b=0;b<MAX_NB_BAND;b++) {
    ti_NbrUseEigen[b] = 0;
    for (int i=0;i<MAX_NBR_PCA_IMA;i++)
      tti_TabKill[i][b] = 0;
  }
}

void mw_Param::operator () (int argc, char *argv[]) {

  int c,i,b; 
  Bool readfile = False, Optf = False, OptL = False;;
  while ((c = GetOpt(argc,argv,"t:m:g:c:n:s:S:x:O:y:WwrkRapvLbPF:K:C:G:U:T:M:B:N:i:ADzZ:")) != -1) {
	
    switch (c) {

    case 't': /* -d <type> type of transform */
      readint (i, "bad type of multiresolution transform: %s\n");
      if (testborn (1 ,NBR_TRANSFORM, i, "bad type of transform: %s\n"))
	e_Transform = (type_transform) (i-1);
      e_TransfOpt = True;
      break;
    
    case 'T':
      Optf = True;
      e_SBFilter = get_filter_bank(OptArg);
      break;     
    
    case 'L':
      e_Norm = NORM_L2;OptL = True;break; 
                    
    case 'm':
      readint (i, "Error: bad type of noise: %s\n");
      if (testborn(1, NBR_NOISE, i, "Error: bad type of noise: %s\n"))
	e_TypeNoise = (type_noise) (i-1);
      e_WithNoiseModel = True;	
      break;

    case 'n': /* -n <NbPlan> */
      readint (i_NbPlan, "bad number of scales: %s\n");
      testborn (2, MAX_SCALE, i_NbPlan, "bad number of scales: %s\n");
      e_NscaleOpt = True;
      break;
      
    //case 'i': /* -i <NbIter> */
    //  readint (i_NbIter, "bad number of iterarions: %s\n");
    //  testborn (1, MAX_ITER, i_NbIter, "bad number of iterations: %s\n");
    //  break;

    case 'g': /* -g <sigma_noise> */
      readfloat (f_NoiseIma,  "Error: bad sigma noise: %s\n");
      e_TypeNoise = NOISE_GAUSSIAN;
      e_WithNoiseModel = True;
      break;

    case 'c': /* -c <gain sigma mean> */
      readfloat (f_Gain,  "Error: bad sigma noise: %s\n");
      readfloat (f_SigmaGauss,  "Error: bad sigma noise: %s\n");
      readfloat (f_MeanGauss,  "Error: bad sigma noise: %s\n");
      e_TypeNoise = NOISE_POISSON;
      e_WithNoiseModel = True;
      f_NoiseIma = 1.;
      break;

    case 's': /* -s <nsigma> */
      readfloat (f_NSigmaMr2d, "Error: bad NSigma: %s\n");
      break;
      
    case 'S': /* -S <nsigma> */
      readfloat (f_NSigmaPCA, "Error: bad NSigma: %s\n");
      break;
       
    case 'F':
      readint (i, "Error: bad number of eigen vectors: %s\n");
      for (b=0;b<MAX_NB_BAND;b++) ti_NbrUseEigen[b] = i;
      break;

    case 'K':
      readint (i, "Error: bad eigen vector number: %s\n");
      testborn (1, MAX_NBR_PCA_IMA, i, "Error: bad eigen vector number: %s\n");
      for (b=0;b<MAX_NB_BAND;b++) tti_TabKill[i-1][b] = 1;
      break;
          
    case 'C':
      readfloat (f_CvgParam, "Error: bad Cvg Parameter: %s\n");                
      if (f_CvgParam <= 0.) f_CvgParam = 0.1;
      break;
      
    case 'G':
      readfloat (f_RegulVal, "Error: bad Regularization Parameter: %s\n");                
      if (f_RegulVal < 0.) f_RegulVal = 1.;
      break;
      
    case 'U':
      readint (i_TypeOpt, "bad type of regularization %s\n");
      testborn (1, MAX_NBR_TYPE_OPT, i_TypeOpt, "Error: regularization type must be in [1,3]");
      break;
      
    case 'M':
      if (sscanf(OptArg,"%s", tc_NameRMSMap) != 1) {
         fprintf(stderr, "Error: bad file name: %s\n", OptArg); exit(-1);
      } 
      e_UseRMSMap = True;
      break;
      
    case 'B':
      readint (i_SizeBlock, "bad block size (sigma clipping): %s\n");
      break;
    
    case 'N':
      readint (i_NiterClip, "Error: ad it. number for the 3sigma clipping: %s\n");
      break;
      
    case 'i':
      readint (i_NbIterEntrop, "Error: bad number of iterations: %s\n");
      break;
      
    case 'x': /* -x  type of compute correlation */
      readint (i, "bad type of correlation compute : %s\n");
      if (testborn (0 ,MAX_CORREL_COMP-1, i, "bad type of correlation compute: %s\n"))
	e_CorrelComp = (CorrelComputeType) (i);
      break;
      
    case 'O': /* -C correaltion matrix name file */
      if (sscanf(OptArg,"%s", tc_ImportCorMatName) != 1) {
         fprintf(stderr, "Error: bad correlation matrix file name: %s\n", OptArg); exit(-1);
      }       
      readfile = True;
      break;
         
    case 'y': /* -x  type of noise correlation */
      readint (i, "bad type of correlation noise : %s\n");
      if (testborn (0 ,MAX_CORREL_NOISE-1, i, "bad type of correlation noise: %s\n"))
	e_CorrelNoise = (CorrelNoiseType) (i);
      break;
                 
#ifdef LARGE_BUFF
    case 'z':
      if (e_OptZ == True) {
        fprintf(OUTMAN, "Error: Z option already set...\n"); exit(-1);
      }
      e_OptZ = True; break;  
    case 'Z':
      if (sscanf(OptArg,"%d:%s",&i_VMSSize, tc_VMSName) < 1) {
        fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n"); exit(-1);
      }
      if (e_OptZ == True) {
        fprintf(OUTMAN, "Error: z option already set...\n"); exit(-1);
      }
      e_OptZ = True;
      break; 
#endif   
    
    case 'p':
      e_NormCorrelMatrix = False;break;                    
    case 'P':
      e_PositivImag = True;break;    
    case 'b':
      e_MaxImag = True;break;
    case 'a':
      e_UseReconsAdjoint = True;break;
    case 'W':
      e_WriteTransf = True;break;
    case 'w':
      e_WriteEigenOnDisk = True;break;
    case 'r':
      e_WriteSupportMr2d = True;break;
    case 'R':
      e_WriteSupportPCA = True;break;           
    case 'A' : 
      e_DataEdge = True; break;
    case 'D' : 
      e_DataSNR = True;  break;
    case 'k':
      e_TraceParamIn = True;break;
    case 'v':
      e_Verbose = True;break;
    case '?':
      usage(argv);
    }
  }
  
  if (      (e_Transform != TO_UNDECIMATED_MALLAT) 
        &&  (e_Transform != TO_MALLAT) 
        && ((OptL == True) || (Optf == True))) {
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with Mallat transform ... \n");
           exit(0);
  }
  
  if (   (e_Transform == TO_UNDECIMATED_MALLAT)
      || (e_Transform == TO_MALLAT) ) {
    po_FAS = new FilterAnaSynt(); /* mem loss, never free ... */
    po_FAS->alloc(e_SBFilter);
    po_FAS->Verbose = e_Verbose;
  } 
  
  // get optional input file names from trailing parameters and open files 
  // read input image names */
  if (OptInd < argc) strcpy(tc_NameIn, argv[OptInd++]);
  else usage(argv);

  if (OptInd < argc) strcpy(tc_NameOut, argv[OptInd++]);
  else usage(argv);
    
  // make sure there are not too many parameters  
  if (OptInd < argc) {
    fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  } 	
  
  
  // Test inconsistencies in the option call 
  if (e_TypeNoise == NOISE_EVENT_POISSON) {
    if ((e_TransfOpt == True) && (e_Transform != TO_PAVE_BSPLINE)) {
      cerr << "WARNING: with this noise model, only the BSPLINE A TROUS can be used ... " << endl;
      cerr << "        Type transform is set to: BSPLINE A TROUS ALGORITHM " << endl;
    }
    e_Transform = TO_PAVE_BSPLINE;
    if (e_NscaleOpt != True) i_NbPlan = DEF_N_SCALE;
  }
        
  if ((e_TypeNoise == NOISE_CORREL) && (e_UseRMSMap != True)) {
    cerr << endl << endl;
    cerr << "  Error: this noise model need a noise map (-R option) " << endl;
    exit(-1);
  }  
  
  if (e_UseRMSMap == True) {
    if ((e_TypeNoise != NOISE_NON_UNI_ADD) && (e_TypeNoise !=  NOISE_CORREL)) {
      cerr << "Error: this noise model is not correct when RMS map option is set." << endl;
      cerr << "       Valid models are: " << endl;
      cerr << "        " << StringNoise(NOISE_NON_UNI_ADD) << endl;
      cerr << "        " << StringNoise(NOISE_CORREL) << endl;
      exit(-1);
    }
    // Stat_Noise = NOISE_NON_UNI_ADD;
  }
        
  if ((isotrop(e_Transform) == False)
      && ((e_TypeNoise == NOISE_NON_UNI_ADD) || (e_TypeNoise  == NOISE_NON_UNI_MULT))) {
    cerr << endl << endl;
    cerr << "  Error: with this transform, non stationary noise models are not valid : " << StringFilter(FILTER_THRESHOLD) << endl;
    exit(-1);
  }
  
#ifdef LARGE_BUFF
    if (e_OptZ == True) vms_init(i_VMSSize, tc_VMSName, e_Verbose);
#endif 	
	
  read_image (argv);  
  
  // matrix file name
  if (e_CorrelComp == E_IMPORT) {
  
     if (readfile != True) {
        cout << " -O Matrix_File_Name not initialized" << endl;
	exit (-1);
     }
     
     fitsstruct Header;
     po_ImportCorMat = new Ifloat(i_NbImage, i_NbImage, "");
     io_read_ima_float(tc_ImportCorMatName, *po_ImportCorMat, &Header);
     int Nl = po_ImportCorMat->nl();
     int Nc = po_ImportCorMat->nc();     
     if (Nl != Nc || Nl != i_NbImage) {
        cout << " bad size of Matrix File Name" << endl;
	exit (-1);     
     }
  } 
  
  if (e_TraceParamIn) trace ();
  
}

void mw_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options input_image output_image\n", argv[0]);
  manline(); 

  transform_usage (e_Transform);      manline();  
  wk_filter (F_MALLAT_7_9);           manline();
  //wk_noise_usage ();                     manline();
  nbr_scale_usage (i_NbPlan);         manline();
  gauss_usage ();                     manline();
  //ccd_usage ();                       manline();
  // normalize_correl_matrix ();         manline();
  correl_compute (DEF_CORREL_COMP);   manline(); 
  import_correl_file();               manline();
  correl_noise (DEF_CORREL_NOISE);    manline();
  nsigma_mr2d (DEF_SIGMA_MR2D);       manline();
  //nsigma_pca (DEF_SIGMA_PCA);         manline();
  positiv_imag_constraint();          manline();
  max_imag_constraint();              manline();
  write_eigen_on_disk ();             manline();
  nb_eigen_used ();                   manline();
  vector_not_used ();                 manline();
   
  wk_conv_param (DEF_CONV_PARAM);        manline();
  wk_regul_param (DEF_REGUL_PARAM);      manline();
  wk_type_opt ();                        manline();
  // wk_model_edge ();                      manline();
  wk_data_snr ();                        manline();
  // wk_rms_map () ;                        manline();
  // wk_size_block_sig_clip (DEF_SIZE_BLK); manline();
  // wk_nb_iter_sig_clip (DEF_ITER_SIGCLP); manline();
  wk_nb_iter_entrop (DEF_ITER_ENTROP);   manline();
  
  //write_support_mr2d ();              manline();
  //write_support_pca ();               manline();
  verbose ();                         manline();
  exit(-1);
}

void mw_Param::read_image (char *argv[]) {

  if (e_Verbose) cout << "Name file IN : " << tc_NameIn << endl;
  
  // read the data (-> fltarray)
  fits_read_fltarr (tc_NameIn, ao_3dDataIn);
  i_NbCol = ao_3dDataIn.nx(); //col
  i_NbLin = ao_3dDataIn.ny(); //lin
  i_NbImage = ao_3dDataIn.nz(); //image

  // convert data -> Ifloat
  po_TabIma = new Ifloat [i_NbImage];
  for (int k=0;k<i_NbImage;k++) {
    po_TabIma[k].alloc (i_NbLin, i_NbCol, "");
    for (int i=0;i<i_NbLin;i++)
      for (int j=0;j<i_NbCol;j++)
	po_TabIma[k](i,j) = ao_3dDataIn(j,i,k);
  }

  // verify that all Mr2d have same size
  if (!control_image (i_NbImage, po_TabIma)) exit(-1);
  
  // set by default the number of eigen vectors used for the reconstruction
  // to the number of images
  for (int b=0;b<MAX_NB_BAND;b++)
    if (ti_NbrUseEigen[b] < 1 || ti_NbrUseEigen[b] >= i_NbImage) 
      ti_NbrUseEigen[b] = i_NbImage;
}


void mw_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransform(e_Transform) << endl;
   if (e_WithNoiseModel) cout << " -- [-m:] : with noise model : " << StringNoise(e_TypeNoise) << endl;
   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   //cout << " -- [-i:] : number of iteration : " << i_NbIter << endl;
   cout << "          : number of iteration : " << i_NbIter << endl;
   if (f_NoiseIma != 0.) cout << " -- [-g:] : gauss noise : " << f_NoiseIma << endl;
   if (f_Gain != 0.) cout << " -- [-c:] : gain : " << f_Gain << ", sigma : " << f_SigmaGauss << ", mean : " << f_MeanGauss << endl;
   if (e_NormCorrelMatrix) cout << " -- [-p]  : normalize correl matrix" << endl;
   if (e_WithNoiseModel) {
      cout << " -- [-x:] : correlation compute : " << CorCompTransform(e_CorrelComp) << endl;
      if (e_CorrelComp == E_IMPORT)
         cout << " -- [-O:] : import correlation file : " << tc_ImportCorMatName << endl;cout << " -- [-y:] : correlation noise : " << CorNoiseTransform(e_CorrelNoise) << endl;         
      cout << " -- [-s:] : multiresol nsigma : " << f_NSigmaMr2d << endl;
      cout << " -- [-S:] : pca nsigma : " << f_NSigmaPCA << endl;
   }  
   cout << " -- [-F:] : number of eigen vector used : " << ti_NbrUseEigen[0] << endl; 
   
   cout << " -- [-C:] : convergence parametre : " << f_CvgParam << endl;
   cout << " -- [-G:] : regulation parametre : " << f_RegulVal << endl;
   cout << " -- [-T:] : type optimisation : " << getTypeOpt (i_TypeOpt) << endl;
   if (e_DataEdge) cout << " -- [-A]  : use model edge" << endl;
   if (e_DataSNR)  cout << " -- [-D]  : use data SNR" << endl;
   if (e_UseRMSMap) cout << "-- [-M]  : RMS Map file : " << tc_NameRMSMap << endl;
   cout << " -- [-B:] : block size sigma clipping : " << i_SizeBlock << endl;
   cout << " -- [-N:] : number of iteration sigma clipping : " << i_NiterClip << endl;
   cout << " -- [-i:] : number of iteration for entrop programm : " << i_NbIterEntrop << endl;
   
   if (e_UseReconsAdjoint) cout << " -- [-a]  : use rec_adjoint in recons process" << endl;
   if (e_PositivImag) cout << " -- [-P]  : positivity constraint" << endl;
   if (e_MaxImag) cout << " -- [-b]  : max image constraint (255)" << endl;   
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect images" << endl;
   if (e_WriteEigenOnDisk) cout << " -- [-w]  : write eigen vector images" << endl;  
   if (e_WriteSupportMr2d) cout << " -- [-r]  : write Multiresol support" << endl;
   if (e_WriteSupportPCA) cout << " -- [-R]  : write PCA support" << endl;
   if (e_Verbose) cout << " -- [-v]  : verbose" << endl;
}



/******************************************************************************
classe mw_Result
******************************************************************************/

mw_Result::mw_Result (mw_Param* ppo_Param) {

   po_Param = ppo_Param;
   // allocate fltarray for write result
   o_3dDataOut.alloc (po_Param->i_NbCol, po_Param->i_NbLin, po_Param->i_NbImage);
   // allocate the memory for all multiresol
   po_TabMr2d = new MultiResol [po_Param->i_NbImage];
   for (int k=0; k<po_Param->i_NbImage; k++)  
      po_TabMr2d[k].alloc (po_Param->i_NbLin, po_Param->i_NbCol, po_Param->i_NbPlan, 
                           po_Param->e_Transform, po_Param->po_FAS, po_Param->e_Norm);
  
   // init number of band
   po_Param->i_NbBand = po_TabMr2d[0].nbr_band();
 
   
   // class for image correlation analysis
   o_Mr2dPca.alloc (po_Param->i_NbImage, po_Param->i_NbBand);
   for (int b=0; b<po_Param->i_NbBand-1; b++) {
      o_Mr2dPca.MatCor[b].Verbose = po_Param->e_Verbose;
      o_Mr2dPca.MatCor[b].NormAna = po_Param->e_NormEigen;
   }
   			   
			   
   // allocate memory for Noise Model classes
   if (po_Param->e_WithNoiseModel) 
      po_TabNoiseModel = new MRNoiseModel [po_Param->i_NbImage];
   else 
      po_TabNoiseModel = (MRNoiseModel*)NULL;
}


mw_Result::~mw_Result () {
  delete [] po_TabMr2d;
  if (po_Param->e_WithNoiseModel) delete [] po_TabNoiseModel;
}
 
 
void mw_Result::res_FirstCompute () {

  // create all multiresol
  if (po_Param->e_Verbose) cout << "Multiresol transform of data ..." << endl;
  for (int k=0; k<po_Param->i_NbImage; k++)
    po_TabMr2d[k].transform (po_Param->po_TabIma[k]);
    
  if (po_Param->e_WriteTransf) {
     char Mr2dRecons[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(Mr2dRecons, "mr2d_Transf_%d",k+1);
        po_TabMr2d[k].write (Mr2dRecons);
     } 
  }

  // create all noise model
  if (po_Param->e_WithNoiseModel) 
     res_InitNoiseModelMr2d();
  
  if (po_Param->e_WriteTransf) {
     char Mr2dRecons[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(Mr2dRecons, "mr2d_Transf_Thresh_%d",k+1);
        po_TabMr2d[k].write (Mr2dRecons);
     } 
  } 
  
  // calculate the correlation matrix and its eigen values
  if (po_Param->e_Verbose) cout << "Principal component analysis ..." << endl;
  if (po_Param->e_WithNoiseModel) 
     o_Mr2dPca.compute (po_TabMr2d, po_TabNoiseModel,
                        po_Param->e_NormCorrelMatrix,
                        po_Param->e_CorrelComp, 
                        po_Param->e_CorrelNoise, 
			po_Param->po_ImportCorMat);
  else
     o_Mr2dPca.compute(po_TabMr2d, (MRNoiseModel*)NULL);
 
  // print the result to the standard output
  if (po_Param->e_Verbose) o_Mr2dPca.print();
  
  // calculate the eigen vector images 
  if (po_Param->e_Verbose) cout << "Compute eigen vector images ..." << endl;
  o_Mr2dPca.pca_TransfSignal (po_TabMr2d, po_TabMr2d);
  
        if (po_Param->e_WriteTransf) {
          char NameEigen[256];
             for (int k=0; k < po_Param->i_NbImage; k++) {
                sprintf(NameEigen, "pca_%dEigenVectors_Aft",k+1);
                po_TabMr2d[k].write (NameEigen);
             }
          }
    
  if (po_Param->e_WithNoiseModel) {
     res_InitNoiseModelPca ();
     //o_Mr2dPca.pca_Thresold (po_TabMr2d, po_TabNoiseModel);
     
     // entrop filter
     // -------------
     if (po_Param->e_Verbose) cout << "Multiscale Entropy Filtering ..." << endl;
     for (int k=0; k<po_Param->i_NbImage; k++) {
     
        fltarray o_TabAlpha(po_Param->i_NbBand-1);
	for (int b=0; b < po_Param->i_NbBand-1; b++) o_TabAlpha(b) = po_Param->f_RegulVal;
	
	MultiResol  *o_MR_Model = NULL;
        if ((po_Param->e_DataEdge == True) && (po_Param->e_Transform == TO_PAVE_BSPLINE)) {
	   o_MR_Model = new MultiResol(po_Param->i_NbLin, po_Param->i_NbCol,
	                po_Param->i_NbPlan, po_Param->e_Transform, "Model");
           mr_getmodel_from_edge(po_TabMr2d[k], *o_MR_Model, po_TabNoiseModel[k]);
	   (*o_MR_Model).write("xx_model.mr");	   
	}
	
	for (int b=0; b<po_Param->i_NbBand-1; b++) {
	   
           mw_mcfilter (b, po_TabMr2d[k], po_TabNoiseModel[k], po_Param->i_TypeOpt, 
	                o_TabAlpha, po_Param->e_DataSNR, po_Param->e_DataEdge,
			o_MR_Model, po_Param->f_CvgParam, po_Param->i_NbIterEntrop, 
			po_Param->e_PositivImag, po_Param->e_Verbose);
			
	}
     }
  }
  
  if (po_Param->e_WriteEigenOnDisk) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "wk_ev_%d",k+1);
        po_TabMr2d[k].write (NameEigen);
     }
  }
}



void mw_Result::res_CurrentCompute () {

   // not used
   cout << "!!!!! NOT USED !!!!!" << endl;
}


void mw_Result::res_Recons () {

   if (po_Param->e_Verbose) cout << "Reconstruction process ..." << endl;
   int ai_NbEigen=0;
   for (int i=0; i<po_Param->ti_NbrUseEigen[0]; i++) //same value for all bands
      if (po_Param->tti_TabKill[i][0] == 0) ai_NbEigen++;
   if (po_Param->e_Verbose) cout << "Number of eigen vector used for the reconstruction : " ;
   if (po_Param->e_Verbose) cout << "Ne = " << ai_NbEigen << " (in all band)" << endl;
  
   // apply the inverse reconstruction from a subset of eigen vectors
   o_Mr2dPca.invsubtransform (po_TabMr2d, po_TabMr2d, po_Param->ti_NbrUseEigen,
                             po_Param->tti_TabKill);
			     		      
   // inv transform Mr2d
   for (int k=0; k<po_Param->i_NbImage; k++)
      if (po_Param->e_UseReconsAdjoint) po_TabMr2d[k].rec_adjoint(po_Param->po_TabIma[k]);
      else po_TabMr2d[k].recons(po_Param->po_TabIma[k]);
      
   // Write Mr2d recons
   if (po_Param->e_WriteTransf == True) {
      char Mr2dRecons[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
         sprintf(Mr2dRecons, "mr2d_pca_InvTransf_%d",k+1);
         po_TabMr2d[k].write (Mr2dRecons);
         sprintf(Mr2dRecons, "mr2d_pca_InvImag_%d", k+1);
         io_write_ima_float(Mr2dRecons, po_Param->po_TabIma[k]);
      }
   }    
}

void mw_Result::res_Write () {

  if (po_Param->e_Verbose) cout << "Write result ..." << endl;
  for (int k=0;k<po_Param->i_NbImage;k++)
    for (int i=0;i<po_Param->i_NbLin;i++)
      for (int j=0;j<po_Param->i_NbCol;j++)
	o_3dDataOut(j,i,k) = po_Param->po_TabIma[k](i,j);

  // write the data
  fits_write_fltarr (po_Param->tc_NameOut, o_3dDataOut);
} 

void mw_Result::res_InitNoiseModelMr2d () {

   // noise model class initialization
   if (po_Param->e_Verbose) cout << "Create mr2d noise model ..." << endl;
   for (int k=0; k<po_Param->i_NbImage; k++) {

      po_TabNoiseModel[k].alloc (po_Param->e_TypeNoise, po_Param->i_NbLin, 
                                 po_Param->i_NbCol, po_Param->i_NbPlan,
                                 po_Param->e_Transform);

      if (po_Param->f_NoiseIma > FLOAT_EPSILON) 
         po_TabNoiseModel[k].SigmaNoise = po_Param->f_NoiseIma;
      po_TabNoiseModel[k].CCD_Gain = po_Param->f_Gain;
      po_TabNoiseModel[k].CCD_ReadOutSigma = po_Param->f_SigmaGauss;
      po_TabNoiseModel[k].CCD_ReadOutMean = po_Param->f_MeanGauss;
      
  // entrop
      po_TabNoiseModel[k].NiterSigmaClip = po_Param->i_NiterClip;
      po_TabNoiseModel[k].SizeBlockSigmaNoise = po_Param->i_SizeBlock;
  
      if (po_Param->e_UseRMSMap == True) {
        po_TabNoiseModel[k].UseRmsMap = True;
        io_read_ima_float(po_Param->tc_NameRMSMap, po_TabNoiseModel[k].RmsMap);
      }
      
      if (po_Param->e_TypeNoise == NOISE_SPECKLE) 
         po_TabNoiseModel[k].SigmaApprox = True;
  // end entrop
  
      if (po_Param->f_NSigmaMr2d >= 0)
         for (int i=0; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaMr2d;

      //po_TabNoiseModel[k].SupIsol = po_Param->e_SupIsolPixel;
      
      po_TabNoiseModel[k].model (po_Param->po_TabIma[k]);
      
      if (po_Param->e_WriteSupportMr2d) {
         char NameSupport[256];
         sprintf(NameSupport, "mr2d_Support_%d",k+1);
         po_TabNoiseModel[k].write_support_mr (NameSupport);
      }
   }
	 
   //if (po_Param->e_ModNoiseModelMr2d) res_ModSupport (po_Param->f_ModNsigmaMr2d); 
}


void mw_Result::res_InitNoiseModelPca () {

   // noise model class initialization
   if (po_Param->e_Verbose) cout << "Create pca noise model ..." << endl;
   int k;
   o_Mr2dPca.pca_ComputeNoise (po_TabMr2d, po_TabNoiseModel, po_TabNoiseModel);
   
   for (k=0; k<po_Param->i_NbImage; k++) {
   
      //if (po_Param->e_DilateSupportPCA) po_TabNoiseModel[k].DilateSupport = True;
      if (po_Param->f_NSigmaPCA >= 0)
         for (int i=0; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaPCA;
      else {
         po_TabNoiseModel[k].NSigma[0] = 4;
         for (int i=1; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = 3;
      }
                   	 
      po_TabNoiseModel[k].set_support (po_TabMr2d[k]);
     
      
      //if (po_Param->e_DestroyRings) res_DestroyRingsInPcaSup ();
   }
   //if (po_Param->e_ModNoiseModelEntrop) res_ModSupport (po_Param->f_ModNsigmaEntrop);
    
   for (k=0; k<po_Param->i_NbImage; k++)
      if (po_Param->e_WriteSupportPCA) {
      char NameSupport[256];
      sprintf(NameSupport, "pca_Support_%d",k+1);
      po_TabNoiseModel[k].write_support_mr (NameSupport);
   }
}


to_Param* mw_Result::res_GetpParam () {return po_Param;}




/******************************************************************************
classe pca_Iter
******************************************************************************/
Bool mw_Iter::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > po_Result->po_Param->i_NbIter) { 
   
      for (int k=0;k<po_Result->po_Param->i_NbImage;k++) 
         for (int i=0;i<po_Result->po_Param->i_NbLin;i++)
            for (int j=0;j<po_Result->po_Param->i_NbCol;j++) {
               po_Result->po_Param->po_TabIma[k](i,j) = o_3dResidu(i,j,k); 
	       if (   po_Result->po_Param->e_PositivImag 
	           && po_Result->po_Param->po_TabIma[k](i,j)<0.) 
                  po_Result->po_Param->po_TabIma[k](i,j) = 0.0;
	       if (   po_Result->po_Param->e_MaxImag 
	           && po_Result->po_Param->po_TabIma[k](i,j)>MAX_IMAG_VALUE) 
                  po_Result->po_Param->po_TabIma[k](i,j) = MAX_IMAG_VALUE;		  
		  
	    }
	            
      return False;
   } else return True;
}
 
void mw_Iter::iter_Residu () {

   for (int k=0;k<po_Result->po_Param->i_NbImage;k++)
      for (int i=0;i<po_Result->po_Param->i_NbLin;i++)
         for (int j=0;j<po_Result->po_Param->i_NbCol;j++) {
	 
 	    // positiv constraint
	    if (    po_Result->po_Param->e_PositivImag 
	         && po_Result->po_Param->po_TabIma[k](i,j)<0.) 
               po_Result->po_Param->po_TabIma[k](i,j) = 0.0; 
	    if (    po_Result->po_Param->e_MaxImag 
	         && po_Result->po_Param->po_TabIma[k](i,j)>MAX_IMAG_VALUE) 
               po_Result->po_Param->po_TabIma[k](i,j) = MAX_IMAG_VALUE;
	       	          
            // I(i) = I(i) + E(i)
            o_3dResidu(i,j,k) += po_Result->po_Param->po_TabIma[k](i,j);
	                
	    // E(i) = Einit - I(i)
	    po_Result->po_Param->po_TabIma[k](i,j) = o_3dInit(i,j,k) - o_3dResidu(i,j,k);
	    
   }    
}

to_Result* mw_Iter::iter_getpResult() {return po_Result;}




