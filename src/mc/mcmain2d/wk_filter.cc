/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.6
**
**    Author: 02/29/00
**
**    Date:  00/03/30
**    
**    File:  mc_pcafilter.cc
**
*******************************************************************************/

#include "wk_filter.h"

#define MAX_ITER 100
#define DEF_SIGMA_PCA -1
#define DEF_SIGMA_MR2D 3

#define MAX_IMAG_VALUE 255

#define DEF_CORREL_COMP (CorrelComputeType) E_LOCAL
#define DEF_CORREL_NOISE (CorrelNoiseType) E_THRESHOLD

extern char *OptArg;
extern int   OptInd;

int main(int argc, char *argv[]) {

   extern softinfo Soft;
   Soft.mr3();
   lm_check(LIC_MR3);
   
  // licence manager
  //lm2();

  // init local var
  pca_Param o_Param;

  // Read param IN
  o_Param (argc, argv);

  {  
  // CPUTime t;  
   
  // init Result and Iter classes
  pca_Result o_Result (&o_Param);   
  pca_Iter o_Iter (&o_Result); 
       
  for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {} 

  // write result
  o_Result.res_Write ();
  }     
  exit(0);
}





/******************************************************************************
classe pca_Param
******************************************************************************/
 
pca_Param::pca_Param () {

  //default init of read var
  e_Transform =           DEFAULT_TRANSFORM;   
  e_Norm =                NORM_L1; 
  e_SBFilter =            F_MALLAT_7_9;
  po_FAS =                (FilterAnaSynt*)NULL;  
  i_NbPlan =              DEFAULT_NBR_SCALE;
  i_NbIter =              1;
  e_WithNoiseModel =      True;         // always use noise model, def gaussian
  e_ModNoiseModelPCA =    False;  
  e_ModNoiseModelMr2d =   False;
  e_Mr2dLinearThreshold = False;
  e_PCALinearThreshold =  False;
  e_UseReconsAdjoint =    False;
  e_TypeNoise =           DEFAULT_STAT_NOISE;
  e_PositivImag =         False;
  e_MaxImag =             False;
  e_CorrelComp =          (CorrelComputeType)DEF_CORREL_COMP;
  po_ImportCorMat =       NULL;
  e_CorrelNoise =         (CorrelNoiseType)DEF_CORREL_NOISE;
  e_NormCorrelMatrix =    True;
  f_NoiseIma =            0.;
  f_Gain =                0.;
  f_SigmaGauss =          0.;
  f_MeanGauss =           0.;
  f_NSigmaMr2d =          0.; 
  f_NSigmaPCA =           DEF_SIGMA_PCA;
  f_ModNsigmaMr2d =       100;
  f_ModNsigmaPCA =        100;
  e_Verbose =             False;
  e_WriteTransf =         False;
  e_WriteEigenOnDisk =    False;
  e_NormEigen =           False;
  e_WriteSupportMr2d =    False;
  e_WriteSupportPCA =     False;
  e_DilateSupportPCA =    False;
  e_SupIsolPixel =        False;
  e_DestroyRings =        False;
  e_TraceParamIn =        False;
  e_WriteCorrelMatrix =   False;
  for (int b=0;b<MAX_NB_BAND;b++) {
    ti_NbrUseEigen[b] = 0;
    for (int i=0;i<MAX_NBR_PCA_IMA;i++)
      tti_TabKill[i][b] = 0;
  }
}

void pca_Param::operator () (int argc, char *argv[]) {

  int c,i,b;
  Bool readfile = False, Optf = False, OptL = False;
  while ((c = GetOpt(argc,argv,"t:m:g:c:n:s:S:M:N:x:T:y:O:WwrdkRpDaevulLJPbCF:K:")) != -1) {
	
    switch (c) {

    case 't': /* -d <type> type of transform */
      readint (i, "bad type of multiresolution transform: %s\n");
      if (testborn (1 ,NBR_TRANSFORM, i, "bad type of transform: %s\n"))
	e_Transform = (type_transform) (i-1);
      break;
            
    case 'm':
      readint (i, "Error: bad type of noise: %s\n");
      if (testborn(1, NBR_NOISE, i, "Error: bad type of noise: %s\n"))
	e_TypeNoise = (type_noise) (i-1);
      e_WithNoiseModel = True;	
      break;

    case 'n': /* -n <NbPlan> */
      readint (i_NbPlan, "bad number of scales: %s\n");
      testborn (2, MAX_SCALE, i_NbPlan, "bad number of scales: %s\n");
      break;
      
    case 'i': /* -i <NbIter> */
      readint (i_NbIter, "bad number of iterarions: %s\n");
      testborn (1, MAX_ITER, i_NbIter, "bad number of iterations: %s\n");
      break;

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
      
    case 'M': /* -M <nsigma> */
      e_ModNoiseModelPCA = True;
      readfloat (f_ModNsigmaPCA, "Error: bad PCA ModNSigma: %s\n");
      break;   
          
    case 'N': /* -N <nsigma> */
      e_ModNoiseModelMr2d = True;
      readfloat (f_ModNsigmaMr2d, "Error: bad Mr2d ModNSigma: %s\n");
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
      
   case 'x': /* -x  type of compute correlation */
      readint (i, "bad type of correlation compute : %s\n");
      if (testborn (0 ,MAX_CORREL_COMP-1, i, "bad type of correlation compute: %s\n"))
	e_CorrelComp = (CorrelComputeType) (i);
      break;
      
   case 'O': /* -C correaltion matrix name file */
      if (sscanf(OptArg,"%s", tc_ImportCorMatName) != 1) {
         fprintf(stderr, "Error: bad correlation file name: %s\n", OptArg); exit(-1);
      } 
      readfile = True;
      break; 
        
   case 'y': /* -x  type of noise correlation */
      readint (i, "bad type of correlation noise : %s\n");
      if (testborn (0 ,MAX_CORREL_NOISE-1, i, "bad type of correlation noise: %s\n"))
	e_CorrelNoise = (CorrelNoiseType) (i);
      break;
 	   
   case 'T':
      Optf = True;
      e_SBFilter = get_filter_bank(OptArg);
      break;     
    
    case 'L':
      e_Norm = NORM_L2;OptL = True;break;         
    case 'p':
      e_NormCorrelMatrix = True;break;                    
    case 'C':
      e_WriteCorrelMatrix = True;break;
    case 'P':
      e_PositivImag = True;break;
    case 'b':
      e_MaxImag = True;break;
    case 'l':
      e_Mr2dLinearThreshold = True;break;
    case 'J':
      e_PCALinearThreshold = True;break;
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
    case 'u':
      e_SupIsolPixel = True;break;
    case 'd':
      e_DilateSupportPCA = True;break;
    case 'D':
      e_DestroyRings = True;break;
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

void pca_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options input_image output_image\n", argv[0]);
  manline(); 

  transform_usage (e_Transform);      manline();
  wk_filter (F_MALLAT_7_9);           manline();
  //wk_noise_usage ();                  manline();
  nbr_scale_usage (i_NbPlan);         manline();
  gauss_usage ();                     manline();
  //ccd_usage ();                       manline();
  correl_compute (DEF_CORREL_COMP);   manline();
  // normalize_correl_matrix ();         manline();
  import_correl_file();               manline();
  correl_noise (DEF_CORREL_NOISE);    manline();
  nsigma_mr2d (DEF_SIGMA_MR2D);       manline();
  nsigma_pca (3.);         manline();
  positiv_imag_constraint();          manline();
  max_imag_constraint();              manline();  
  //dilate_support_pca ();              manline();
  //destroy_rings_pca();                manline();
  //destroy_isol_pix();                 manline();
  write_eigen_on_disk ();             manline();
  //noise_select ();                    manline();
  nb_eigen_used ();                   manline();
  //write_support_mr2d ();              manline();
  //write_support_pca ();               manline();
  vector_not_used ();                 manline();
  write_correl_matrix ();             manline();
  verbose ();                         manline();
  exit(-1);
}

void pca_Param::read_image (char *argv[]) {

  if (e_Verbose) cout << "Name file IN : " << tc_NameIn << endl;
  
  // read the data (-> fltarray)
  //fits_read_fltarr (tc_NameIn, ao_3dDataIn);
  io_3d_read_data(tc_NameIn, ao_3dDataIn);
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


void pca_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransform(e_Transform) << endl;
   if ((e_Transform == TO_MALLAT) || (e_Transform == TO_UNDECIMATED_MALLAT)) {
      cout << " -- [-L]  : Normalisation L2" << endl;
      cout << " -- [-f]  : " << StringSBFilter(e_SBFilter) << endl;
   }
   if (e_WithNoiseModel) cout << " -- [-m:] : with noise model : " << StringNoise(e_TypeNoise) << endl;
   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   cout << " -- [-i:] : number of iteration : " << i_NbIter << endl;
   if (f_NoiseIma != 0.) cout << " -- [-g:] : gauss noise : " << f_NoiseIma << endl;
   if (f_Gain != 0.) cout << " -- [-c:] : gain : " << f_Gain << ", sigma : " << f_SigmaGauss << ", mean : " << f_MeanGauss << endl;
   if (e_NormCorrelMatrix) cout << " -- [-p]  : normalize correl matrix" << endl;
   if (e_WithNoiseModel) {
      cout << " -- [-x:] : correlation compute : " << CorCompTransform(e_CorrelComp) << endl;
      if (e_CorrelComp == E_IMPORT)
         cout << " -- [-O:] : import correlation file : " << tc_ImportCorMatName << endl;
      cout << " -- [-y:] : correlation noise : " << CorNoiseTransform(e_CorrelNoise) << endl;  
      if (f_NSigmaMr2d < 1e-03 )
         cout << " -- [-s:] : multiresol nsigma : 3 " << endl;
      else
         cout << " -- [-s:] : multiresol nsigma : " << f_NSigmaMr2d << endl;
      cout << " -- [-S:] : pca nsigma : " << f_NSigmaPCA << endl;
   }  
   if (e_ModNoiseModelMr2d) 
      cout << " -- [-N:] : mr2d mod nsigma : " << f_ModNsigmaMr2d << endl;
   if (e_ModNoiseModelPCA) 
      cout << " -- [-M:] : pca mod nsigma : " << f_ModNsigmaPCA << endl;
   cout << " -- [-F:] : number of eigen vector used : " << ti_NbrUseEigen[0] << endl; 
   if (e_Mr2dLinearThreshold) cout << " -- [-l]  : Mr2d linear threshold" << endl;
   if (e_PCALinearThreshold) cout << " -- [-J]  : PCA linear threshold" << endl;  
   if (e_UseReconsAdjoint) cout << " -- [-a]  : use rec_adjoint in recons process" << endl;
   if (e_PositivImag) cout << " -- [-P]  : positivity constraint" << endl;
   if (e_MaxImag) cout << " -- [-b]  : max image constraint (255)" << endl;
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect images" << endl;
   if (e_WriteEigenOnDisk) cout << " -- [-w]  : write eigen vector images" << endl;  
   if (e_WriteSupportMr2d) cout << " -- [-r]  : write Multiresol support" << endl;
   if (e_WriteSupportPCA) cout << " -- [-R]  : write PCA support" << endl;
   if (e_DilateSupportPCA) cout << " -- [-d]  : dilate PCA support" << endl;
   if (e_SupIsolPixel) cout << " -- [-u]  : kill isol pixel in Mr2d support" << endl;
   if (e_DestroyRings) cout << " -- [-D]  : destroy rings in PCA support" << endl;
   if (e_WriteCorrelMatrix) cout << " -- [-C]  : write correlation matrix" << endl;
   if (e_Verbose) cout << " -- [-v]  : verbose" << endl;

}



/******************************************************************************
classe pca_Result
******************************************************************************/

pca_Result::pca_Result (pca_Param* ppo_Param) {

   po_Param = ppo_Param;
   // allocate fltarray for write result
   o_3dDataOut.alloc (po_Param->i_NbCol, po_Param->i_NbLin, po_Param->i_NbImage);
   // allocate the memory for all multiresol
   po_TabMr2d = new MultiResol [po_Param->i_NbImage];
   for (int k=0; k<po_Param->i_NbImage; k++) {
      po_TabMr2d[k].alloc (po_Param->i_NbLin, po_Param->i_NbCol, po_Param->i_NbPlan, 
                           po_Param->e_Transform, po_Param->po_FAS, po_Param->e_Norm);
      //po_TabMr2d[k].TypeNorm = po_Param->e_Norm;
      //po_TabMr2d[k].SBFilter = po_Param->e_SBFilter;  		   
   }
   
   // init number of band
   po_Param->i_NbBand = po_TabMr2d[0].nbr_band();
  			   
   // class for image correlation analysis
   o_Mr2dPca.alloc (po_Param->i_NbImage, po_Param->i_NbBand,
                    po_Param->e_WriteCorrelMatrix);
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


pca_Result::~pca_Result () {
  delete [] po_TabMr2d;
  if (po_Param->e_WithNoiseModel) delete [] po_TabNoiseModel;
}
 
 
void pca_Result::res_FirstCompute () {

  // create all multiresol
  if (po_Param->e_TraceParamIn)
     cout << "Multiresol transform of data ..." << endl;
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
  
  if (po_Param->e_WriteEigenOnDisk) {
     char Mr2dRecons[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(Mr2dRecons, "mr2d_Transf_Thresh_%d",k+1);
        po_TabMr2d[k].write (Mr2dRecons);
     } 
  } 
  
  // calculate the correlation matrix and its eigen values
  if (po_Param->e_TraceParamIn)
     cout << "Principal component analysis ..." << endl;
  if (po_Param->e_WithNoiseModel) 
     o_Mr2dPca.compute(po_TabMr2d, po_TabNoiseModel, 
                       po_Param->e_NormCorrelMatrix,
		       po_Param->e_CorrelComp,
                       po_Param->e_CorrelNoise, 
		       po_Param->po_ImportCorMat);
  else
     o_Mr2dPca.compute(po_TabMr2d, (MRNoiseModel*)NULL);
 
  // print the result to the standard output
  if (po_Param->e_Verbose) o_Mr2dPca.print();
  
  // calculate the eigen vector images 
  if (po_Param->e_TraceParamIn)
     cout << "Compute eigen vector images ..." << endl;
  o_Mr2dPca.pca_TransfSignal (po_TabMr2d, po_TabMr2d);
  
  if (po_Param->e_WriteTransf) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "wk_ev_%d",k+1);
	po_TabMr2d[k].write (NameEigen);
     }
  }  
  
  if (po_Param->e_WithNoiseModel) {
     res_InitNoiseModelPca ();
     o_Mr2dPca.pca_Thresold (po_TabMr2d, po_TabNoiseModel);
  }
  
  if (po_Param->e_WriteTransf) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "pca_%dEigenVectors",k+1);
        po_TabMr2d[k].write (NameEigen);
     }
  }
}


void pca_Result::res_CurrentCompute () {

   // create all multiresol
   if (po_Param->e_TraceParamIn)
      cout << "Multiresol transform of data ..." << endl;
   for (int k=0; k<po_Param->i_NbImage; k++) 
      po_TabMr2d[k].transform (po_Param->po_TabIma[k]);
      
   if (po_Param->e_WriteTransf) {
      char Mr2dRecons[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
         sprintf(Mr2dRecons, "mr2d_Transf_%d",k+1);
         po_TabMr2d[k].write (Mr2dRecons);
      } 
   }
   
   // calculate the eigen vector images and save them on the disk 
   if (po_Param->e_TraceParamIn)
      cout << "Compute eigen vector images ..." << endl;
   o_Mr2dPca.transform(po_TabMr2d, po_TabMr2d);
   if (po_Param->e_WithNoiseModel == True) 
      o_Mr2dPca.pca_Thresold (po_TabMr2d, po_TabNoiseModel); 	
  
   if (po_Param->e_WriteEigenOnDisk) {
      char NameEigen[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "pca_%dEigenVectors",k+1);
        po_TabMr2d[k].write (NameEigen);
      }
   } 
}


void pca_Result::res_Recons () {

   if (po_Param->e_TraceParamIn)
      cout << "Reconstruction process ..." << endl;
   int ai_NbEigen=0;
   for (int i=0; i<po_Param->ti_NbrUseEigen[0]; i++) //same value for all bands
      if (po_Param->tti_TabKill[i][0] == 0) ai_NbEigen++;
   if (po_Param->e_TraceParamIn)
      cout << "Number of eigen vector used for the reconstruction : " ;
   if (po_Param->e_TraceParamIn)
      cout << "Ne = " << ai_NbEigen << " (in all band)" << endl;
  
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


void pca_Result::res_Write () {

  if (po_Param->e_TraceParamIn) cout << "Write result ..." << endl;
  for (int k=0;k<po_Param->i_NbImage;k++)
    for (int i=0;i<po_Param->i_NbLin;i++)
      for (int j=0;j<po_Param->i_NbCol;j++)
	o_3dDataOut(j,i,k) = po_Param->po_TabIma[k](i,j);

  // write the data
  io_3d_write_data(po_Param->tc_NameOut, o_3dDataOut);
  // fits_write_fltarr (po_Param->tc_NameOut, o_3dDataOut);

  // prov write 
  /*char NameOut[512];
  for (int l=0; l<po_Param->i_NbImage; l++) {
    sprintf(NameOut, "%s_%d", po_Param->tc_NameOut, l+1);
    io_write_ima_float(NameOut, po_Param->po_TabIma[l]);
  }*/
} 


void pca_Result::res_InitNoiseModelMr2d () {

   // noise model class initialization
   if (po_Param->e_TraceParamIn) cout << "Create mr2d noise model ..." << endl;
   for (int k=0; k<po_Param->i_NbImage; k++) {

      po_TabNoiseModel[k].alloc (po_Param->e_TypeNoise, po_Param->i_NbLin, 
                                 po_Param->i_NbCol, po_Param->i_NbPlan,
                                 po_Param->e_Transform);

      if (po_Param->f_NoiseIma > FLOAT_EPSILON) 
         po_TabNoiseModel[k].SigmaNoise = po_Param->f_NoiseIma;
      po_TabNoiseModel[k].CCD_Gain = po_Param->f_Gain;
      po_TabNoiseModel[k].CCD_ReadOutSigma = po_Param->f_SigmaGauss;
      po_TabNoiseModel[k].CCD_ReadOutMean = po_Param->f_MeanGauss;

      if (po_Param->f_NSigmaMr2d > 1e-03)
         for (int i=0; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaMr2d;
      else {
         po_TabNoiseModel[k].NSigma[0] = 4.0;
         for (int i=1; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = 3.0;     
      }
      
      po_TabNoiseModel[k].SupIsol = po_Param->e_SupIsolPixel;
      
      po_TabNoiseModel[k].model (po_Param->po_TabIma[k],
                                 po_TabMr2d[k]);
      
      if (po_Param->e_WriteSupportMr2d) {
         char NameSupport[256];
         sprintf(NameSupport, "mr2d_Support_%d",k+1);
         po_TabNoiseModel[k].write_support_mr (NameSupport);
      }
   }
	 
   if (po_Param->e_ModNoiseModelMr2d) res_ModSupport (po_Param->f_ModNsigmaMr2d);
   
//    for (int k=0; k<po_Param->i_NbImage; k++) {
//       cout << "Sigma Mr2d :";
//       for (int i=0; i<po_Param->i_NbPlan; i++) {
//          cout << " p" << i << ":" << po_TabNoiseModel[k].sigma(i,1,1) << ", ";
//       }
//       cout << endl;
//    }
}

void pca_Result::res_InitNoiseModelPca () {

   // noise model class initialization
   if (po_Param->e_TraceParamIn) cout << "Create pca noise model ..." << endl;
   int k;
   o_Mr2dPca.pca_ComputeNoise (po_TabMr2d, po_TabNoiseModel, po_TabNoiseModel);
   
   for (k=0; k<po_Param->i_NbImage; k++) {
   
      if (po_Param->e_DilateSupportPCA) po_TabNoiseModel[k].DilateSupport = True;
      if (po_Param->f_NSigmaPCA >= 0)
         for (int i=0; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaPCA;
      else {
         po_TabNoiseModel[k].NSigma[0] = 4;
         for (int i=1; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = 3;            
      }
      
      po_TabNoiseModel[k].SupIsol = False;
      
      po_TabNoiseModel[k].set_support (po_TabMr2d[k]);
     
      
      if (po_Param->e_DestroyRings) res_DestroyRingsInPcaSup ();
   }
   if (po_Param->e_ModNoiseModelPCA) res_ModSupport (po_Param->f_ModNsigmaPCA);
    
   for (k=0; k<po_Param->i_NbImage; k++)
      if (po_Param->e_WriteSupportPCA) {
      char NameSupport[256];
      sprintf(NameSupport, "pca_Support_%d",k+1);
      po_TabNoiseModel[k].write_support_mr (NameSupport);
   }   
   
//    for (int k=0; k<po_Param->i_NbImage; k++) {
//       cout << "Sigma Pca :";
//       for (int i=0; i<po_Param->i_NbPlan; i++) {
//          cout << " p" << i << ":" << po_TabNoiseModel[k].sigma(i,1,1) << ", ";
//       }
//       cout << endl;
//    }
}


void pca_Result::res_ModSupport (float f_ModNSigma) {

   if (po_Param->e_TraceParamIn) cout << "Modify noise model..." << endl;
   int s;
   fltarray ModNsigma(po_Param->i_NbPlan);
   for (int i=0; i<po_Param->i_NbPlan; i++) 
      ModNsigma(i) = f_ModNSigma;
	  
   for (s=0; s<po_Param->i_NbPlan-1; s++)
      for (int k=0; k<po_Param->i_NbImage; k++) {
      
      	 int Nls = po_TabMr2d[k].size_band_nl(s);
         int Ncs = po_TabMr2d[k].size_band_nc(s);
         for (int i=0; i<Nls; i++)
            for (int j=0; j<Ncs; j++) {
	   
	       if (po_TabNoiseModel[k].support(s,i,j)==VAL_SupOK)
	       
	          //if (s==1) cout << "(" << k << "," << i << "," << j << ") ,"; 
	      
	          for (int l=0; l<po_Param->i_NbImage; l++) {
	      
	             if (    l!=k
		          && po_TabNoiseModel[l].signif((po_TabMr2d[l])(s,i,j), s,i,j, ModNsigma) == True
		          && po_TabNoiseModel[l].support(s,i,j)!=VAL_SupOK) {
		        //cout <<"==> Modif Support False -> True [IN:"<<k<<" [im:"<<l<<",b:"<<s<<",("<<i<<","<<j<<")]"<<endl;
                        po_TabNoiseModel[l].support(s,i,j)=VAL_SupOK*100;
		
	             }
		  }
	    }	          
      }
      
   for (s=0; s<po_Param->i_NbPlan-1; s++)
    for (int k=0; k<po_Param->i_NbImage; k++) {
      
      	 int Nls = po_TabMr2d[k].size_band_nl(s);
         int Ncs = po_TabMr2d[k].size_band_nc(s);
         for (int i=0; i<Nls; i++)
          for (int j=0; j<Ncs; j++) 
	   if (po_TabNoiseModel[k].support(s,i,j)==VAL_SupOK*100)
            po_TabNoiseModel[k].support(s,i,j)=VAL_SupOK;

   }

}



void pca_Result::res_DestroyRingsInPcaSup () {


   for (int k=0; k<po_Param->i_NbImage; k++) 
   for (int s=0; s<po_Param->i_NbPlan; s++) {
      
         Bool modif=True;
      	 int Nls = po_TabMr2d[k].size_band_nl(s);
         int Ncs = po_TabMr2d[k].size_band_nc(s);
	 
	 while (modif) {
	 
	    modif = False;
            for (int i=1; i<Nls-1; i++)
            for (int j=1; j<Ncs-1; j++) {

                  int compt=0;
                  for (int l=-1;l<2;l++)
	          for (int n=-1;n<2;n++) {
	      
	             if (po_TabNoiseModel[k].support(s,i+l,j+n)==VAL_SupOK)
		        compt++;
                  }
	       
	          if (compt >= 5) {
	             if (po_TabNoiseModel[k].support(s,i,j)!=VAL_SupOK) modif = True;
	             po_TabNoiseModel[k].support(s,i,j)=VAL_SupOK;
                  }
	    }
	 }
   }
}      


to_Param* pca_Result::res_GetpParam () {return po_Param;}


/******************************************************************************
classe pca_Iter
******************************************************************************/
Bool pca_Iter::iter_End () {

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

void pca_Iter::iter_Residu () {

   for (int k=0;k<po_Result->po_Param->i_NbImage;k++)
      for (int i=0;i<po_Result->po_Param->i_NbLin;i++)
         for (int j=0;j<po_Result->po_Param->i_NbCol;j++) {
	 
 	    // positiv constraint
	    if (    po_Result->po_Param->e_PositivImag 
	         && po_Result->po_Param->po_TabIma[k](i,j)<0.) 
               po_Result->po_Param->po_TabIma[k](i,j) = 0.0; 
	       
	    if (   po_Result->po_Param->e_MaxImag 
	        && po_Result->po_Param->po_TabIma[k](i,j)>MAX_IMAG_VALUE) 
               po_Result->po_Param->po_TabIma[k](i,j) = MAX_IMAG_VALUE;
		  		          
            // I(i) = I(i) + E(i)
            o_3dResidu(i,j,k) += po_Result->po_Param->po_TabIma[k](i,j);
	                
	    // E(i) = Einit - I(i)
	    po_Result->po_Param->po_TabIma[k](i,j) = o_3dInit(i,j,k) - o_3dResidu(i,j,k);
	    
   }    
}

to_Result* pca_Iter::iter_getpResult() {return po_Result;}







