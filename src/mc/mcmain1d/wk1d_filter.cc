/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.5
**
**    Author: Philippe
**
**    Date:  01/02/09
**    
**    File:  mc1d_pcafilter.cc
**
*******************************************************************************
**
******************************************************************************/

  
 

#include "wk1d_filter.h"

#define MAX_ITER 100
#define DEF_SIGMA_PCA -1
#define DEF_SIGMA_MR1D -1
#define DEFAULT_NBR_SCALE 4
#define DEF_TRANSFORM TO1_PAVE_B3SPLINE
#define DEFAULT_FIRST_SCALE 0
#define MAX_SCALE 10
#define MAX_SPEC_VALUE 255

#define DEF_CORREL_COMP (CorrelComputeType1d) E_LOCAL1D
#define DEF_CORREL_NOISE (CorrelNoiseType1d) E_THRESHOLD1D


extern int  OptInd;
extern char *OptArg;
//extern int  getopt(int argc, char *const*argv, char *opts);
extern Bool one_level_per_pos(type_noise TNoise);


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
  
cout << "!!!!!!!!!!!! " << argv[argc-2] << endl;
      
  // init Result and Iter classes
  pca_Result o_Result (&o_Param);   
     
  // loop ...
  if (!o_Param.e_OptimSoftThreshold) {
     
     pca_Iter o_Iter (&o_Result);
     for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {} 
     
  } else {
     
     pca_SoftIter o_Iter (&o_Result);
     for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {}  
         
  }
   
  // write result
  o_Result.res_Write ();
        
  exit(0);
}


 

/******************************************************************************
classe to_Param
******************************************************************************/
 
pca_Param::pca_Param () {

  //default init of read var
  e_Transform =           DEF_TRANSFORM;     
  e_Norm =                NORM_L1; 
  e_SBFilter =            F_MALLAT_7_9;
  po_FAS =                (FilterAnaSynt*)NULL; 
  e_WithNoiseModel =      True;         // always use noise model, def gaussian 
  e_TypeNoise =           DEFAULT_STAT_NOISE;
  e_PositivSpectre =      True;
  e_MaxSpectre =          False;  
  //e_Border
  i_NbPlan =              DEFAULT_NBR_SCALE;
  i_NbIter =              1;
  i_FirstScale = 	  DEFAULT_FIRST_SCALE;
  //i_NbrUndec =
  //i_NbBand
  f_NoiseSpectre =        0.;
  e_UseReconsAdjoint =    False;
  e_CorrelComp =          (CorrelComputeType1d)DEF_CORREL_COMP;
  po_ImportCorMat =       NULL;
  e_CorrelNoise =         (CorrelNoiseType1d)DEF_CORREL_NOISE;
  f_NSigmaMr1d =          DEF_SIGMA_MR1D; 
  f_NSigmaPCA =           DEF_SIGMA_PCA;
  e_Verbose =             False;
  e_NormCorrelMatrix =    True;
  e_WriteTransf =         False;
  e_NoDisplay =           True;
  e_WriteEigenOnDisk =    False;
  e_WriteSupportMr1d =    False;
  e_WriteSupportPCA =     False;
  //e_DilateSupportPCA =    False;
  e_WriteCorrelMatrix =   False;
  e_SupIsolPixel =        False;
  e_OnlyPositivDetect =   False;
  e_RemoveLastScale =     False;
  e_OptimSoftThreshold =  False;
  e_TresholdInMrSpace =   False;
  e_TraceParamIn =        False;
  for (int b=0;b<MAX_NB_BAND;b++) {
    ti_NbrUseEigen[b] = 0;
    for (int i=0;i<MAX_NBR_PCA_SPEC;i++)
      tti_TabKill[i][b] = 0;
  }
}
 
void pca_Param::operator () (int argc, char *argv[]) {

  int c,i,b;  
  Bool readfile = False, Optf = False;
  
  tc_Cmd[0] = '\0';
  for (int k =0; k < argc; k++) 
     sprintf(tc_Cmd, "%s %s", tc_Cmd, argv[k]);
  
  while ((c = GetOpt(argc,argv,"c:d:F:g:i:K:t:m:M:n:N:O:s:S:T:x:y:aAbCDefJklLopPrRuvWl")) != -1) {
	
    switch (c) {

    case 't': /* -d <type> type of transform */
      readint (i, "bad type of multiresolution transform: %s\n");
      if (testborn (1 ,NBR_TRANSFORM, i, "bad type of transform: %s\n"))
	e_Transform = (type_trans_1d) (i-1);
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
    
   case 'g': /* -g <sigma_noise> */
      readfloat (f_NoiseSpectre,  "Error: bad sigma noise: %s\n");
      e_TypeNoise = NOISE_GAUSSIAN;
      e_WithNoiseModel = True;
      break;  
    
   case 's': /* -s <nsigma> */
      readfloat (f_NSigmaMr1d, "Error: bad NSigma: %s\n");
      break;
      
   case 'S': /* -S <nsigma> */
      readfloat (f_NSigmaPCA, "Error: bad NSigma: %s\n");
      break;  
           
   case 'F':
      readint (i, "Error: bad number of eigen vectors: %s\n");
      for (b=0;b<MAX_NB_BAND;b++)
         ti_NbrUseEigen[b] = i;
      break;

   case 'K':
      readint (i, "Error: bad eigen vector number: %s\n");
      testborn (1, MAX_NBR_PCA_IMA, i, "Error: bad eigen vector number: %s\n");
      for (b=0;b<MAX_NB_BAND;b++) tti_TabKill[i-1][b] = 1;
      break;
            
   case 'x': /* -x  type of compute correlation */
      readint (i, "bad type of correlation compute : %s\n");
      if (testborn (0 ,MAX_CORREL_COMP-1, i, "bad type of correlation compute: %s\n"))
	e_CorrelComp = (CorrelComputeType1d) (i);
      break;
      
   case 'O': /* -C correaltion matrix name file */
      if (sscanf(OptArg,"%s", tc_ImportCorMatName) != 1) {
         fprintf(stderr, "Error: bad correlation file name: %s\n", OptArg); exit(-1);
      } 
      readfile = True;
      break; 
   
    case 'i': /* -i <NbIter> */
      readint (i_NbIter, "bad number of iterarions: %s\n");
      testborn (1, MAX_ITER, i_NbIter, "bad number of iterations: %s\n");
      break;
     
  // case 'u': /* -u number of undecimated scale */
   //   readint (i_NbrUndec, "Error: bad number of scales: %s\n");
   //   if (testborn (0 ,MAX_SCALE-1, i_NbrUndec, "bad number of scales: %s\n"))	
   //   break;
      
    case 'T': /* -T get filter... */
      Optf = True;
      e_SBFilter = get_filter_bank(OptArg);
      break;     
   
    case 'd': /* begin the detection */		
      readint (i_FirstScale, "Error: bad first scale detection %d\n");
      testborn (1, 10, i_FirstScale, "Error: bad filter type: %d\n");
      i_FirstScale--;
      break;
      	
    case 'o' :
      e_OptimSoftThreshold=True;break;
    case 'A':
      e_TresholdInMrSpace=True;break;
    case 'f':
      e_RemoveLastScale=True;break;           
    case 'P':           
      e_PositivSpectre=False;break;  
    case 'a':
      e_UseReconsAdjoint=True;break;
    case 'L':
      e_Norm = NORM_L2;break;     
    case 'p':
      e_NormCorrelMatrix = False;break;              
    case 'r':
      e_WriteSupportMr1d = True;break;
    case 'R':
      e_WriteSupportPCA = True;break;
    case 'u':
      e_SupIsolPixel = True;break;
    case 'b':
      e_OnlyPositivDetect=True;break;
    case 'C':
      e_WriteCorrelMatrix = True;break;
    case 'w':
      e_WriteEigenOnDisk = True;break;
    case 'W':
      e_WriteTransf = True;break;
    case 'k':
      e_TraceParamIn = True;break;
    case 'v':
      e_NoDisplay = False; 
      e_Verbose = True;break;
    case '?':
      usage(argv);
    }
  }  
  
  /*if (      (e_Transform != TO_UNDECIMATED_MALLAT) 
        &&  (e_Transform != TO_MALLAT) 
        && ((OptL == True) || (Optf == True))) {
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with Mallat transform ... \n");
           exit(0);
  }
  
  if (   (e_Transform == TO_UNDECIMATED_MALLAT)
      || (e_Transform == TO_MALLAT) ) {
    po_FAS = new FilterAnaSynt(); // mem loss, never free ... 
    po_FAS->alloc(e_SBFilter);
    po_FAS->Verbose = e_Verbose;
  }*/
   
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

  read_spectre (argv);

  if (e_CorrelComp == E_IMPORT1D) {
   
     if (readfile != True) {
        cout << " -O Matrix_File_Name not initialized" << endl;
	exit (-1);
     }
     
     fitsstruct Header;
     po_ImportCorMat = new Ifloat(i_NbSpectre, i_NbSpectre, "");
     io_read_ima_float(tc_ImportCorMatName, *po_ImportCorMat, &Header);
     int Nl = po_ImportCorMat->nl();
     int Nc = po_ImportCorMat->nc();
     if (Nl != Nc || Nl != i_NbSpectre) {
        cout << " bad size of Matrix File Name" << endl;
	exit (-1);     
     }
  }          
  
  if (e_TraceParamIn) trace ();
}

void pca_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options input_spectra output_spectra\n", argv[0]);
  manline(); 

  //transform_usage (e_Transform);      manline();
  //wk_noise_usage ();                     manline();
  nbr_scale_usage (i_NbPlan);         manline();
  gauss_usage ();                     manline();
  // write_eigen_on_disk ();             manline();
  write_correl_matrix ();             manline();
  //normalize_correl_matrix ();         manline();
  correl_compute1d (DEF_CORREL_COMP); manline();
  import_correl_file();               manline(); 
  correl_noise1d (DEF_CORREL_NOISE);  manline();
  nsigma_mr1d (DEF_SIGMA_MR1D);       manline();
  nsigma_pca (DEF_SIGMA_PCA);         manline();
  positiv_imag_constraint();          manline();
  // max_imag_constraint();              manline();    
  nb_eigen_used ();                   manline();
  vector_not_used ();                 manline();   
  
  
  verbose ();                         manline();
  exit(-1);
} 

void pca_Param::read_spectre (char *argv[]) {

  if (e_Verbose) cout << "Name file IN : " << tc_NameIn << endl;
  
  // read the data (-> fltarray)
  fits_read_fltarr (tc_NameIn, ao_2dDataIn, &ts_Header);
  i_NbPts =     ao_2dDataIn.nx(); //point
  i_NbSpectre = ao_2dDataIn.ny(); //spectre

  // convert data -> Ifloat
  po_TabSpectre = new fltarray [i_NbSpectre];
  for (int k=0;k<i_NbSpectre;k++) {
    (po_TabSpectre[k]).alloc (i_NbPts);
    for (int i=0;i<i_NbPts;i++)
      (po_TabSpectre[k])(i) = ao_2dDataIn(i,k);
  }

  // verify that all Mr2d have same size
  if (!control_spectre (i_NbSpectre, po_TabSpectre)) exit(-1);

  // set by default the number of eigen vectors used for the reconstruction
  // to the number of spectres
  for (int b=0;b<MAX_NB_BAND;b++)
    if (ti_NbrUseEigen[b] < 1 || ti_NbrUseEigen[b] >= i_NbSpectre) 
      ti_NbrUseEigen[b] = i_NbSpectre;
}


void pca_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransf1D(e_Transform) << endl;
   if (e_WithNoiseModel) cout << " -- [-m:] : with noise model : " << StringNoise(e_TypeNoise) << endl;
   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   cout << " -- [-i:] : number of iteration : " << i_NbIter << endl;
   if (e_WithNoiseModel) {
      if (f_NSigmaMr1d >= 0)
         cout << " -- [-s:] : multiresol nsigma : " << f_NSigmaMr1d << endl;
      else cout << " -- [-s:] : multiresol nsigma : 4 3 3 3 3 ..." << endl;
      if (f_NSigmaPCA >= 0)
      cout << " -- [-S:] : pca nsigma : " << f_NSigmaPCA << endl;
      else cout << " -- [-S:] : pca nsigma : 4 3 3 3 3 ..." << endl;
   }  
   if (e_WriteCorrelMatrix) cout << " -- [-C]  : write correlation matrix" << endl;
   cout << " -- [-x:] : correlation compute : " << CorCompTransform1d(e_CorrelComp) << endl;
   if (e_CorrelComp == E_IMPORT1D)
      cout << " -- [-O:] : import correlation file : " << tc_ImportCorMatName << endl;
   if (e_NormCorrelMatrix) cout << " -- [  ]  : normalize correl matrix" << endl;
   else                    cout << " -- [-p]  : do not normalize correl matrix" << endl;
   cout << " -- [-y:] : correlation noise : " << CorNoiseTransform1d(e_CorrelNoise) << endl;
   cout << " -- [-F:] : number of eigen vector used : " << ti_NbrUseEigen[0] << endl; 
   if (e_UseReconsAdjoint) cout << " -- [-a]  : use rec_adjoint in recons process" << endl;
   cout << " -- [-d]  : first scale : " << i_FirstScale << endl;
   if (e_PositivSpectre) cout << " -- [  ]  : positivity constraint" << endl;
   else                  cout << " -- [-P]  : no positivity constraint" << endl;
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect spectres" << endl;
   if (e_WriteSupportMr1d) cout << " -- [-r]  : write Multiresol support" << endl;
   if (e_WriteSupportPCA) cout << " -- [-R]  : write PCA support" << endl;
   if (e_SupIsolPixel) cout << " -- [-u]  : suppres isolated pixel" << endl;
   if (e_OnlyPositivDetect) cout << " -- [-b]  : only detect positiv pixels" << endl;
   if (e_RemoveLastScale) cout << " -- [-f]  : remove last scale" << endl;
   if (e_OptimSoftThreshold) cout << " -- [-o]  : optim : soft threshold..." << endl;
   if (e_TresholdInMrSpace) cout << " -- [-A]  : threshold wavelet coef in mr space" << endl;
   if (e_Verbose) cout << " -- [-v]  : verbose" << endl;

}


/******************************************************************************
classe to_Result
******************************************************************************/

pca_Result::pca_Result (pca_Param* ppo_Param) {

   po_Param = ppo_Param;
   LoopNumber = 0;
   // allocate fltarray for write result
   o_2dDataOut.alloc (po_Param->i_NbPts, po_Param->i_NbSpectre);
   // allocate the memory for all multiresol
   po_TabMr1d = new MR_1D [po_Param->i_NbSpectre];
   for (int k=0; k<po_Param->i_NbSpectre; k++) {
      //po_TabMr1d[k].alloc (po_Param->i_NbPts, po_Param->e_Transform, 
      //                     "PCA", po_Param->i_NbPlan);
      po_TabMr1d[k].alloc (po_Param->i_NbPts, po_Param->e_Transform,
                           po_Param->i_NbPlan, po_Param->po_FAS, 
			   po_Param->e_Norm);      
   }
   
   // init number of band
   po_Param->i_NbBand = po_TabMr1d[0].nbr_band();			   
			   
   // class for spectral correlation analysis
   o_Mr1dPca.alloc (po_Param->i_NbSpectre, po_Param->i_NbBand,
                    po_Param->e_WriteCorrelMatrix);
    
   for (int b=0; b<po_Param->i_NbBand; b++) {
      o_Mr1dPca.MatCor[b].Verbose = po_Param->e_Verbose;
   }  		   
   
   // allocate memory for Noise Model classes
   if (po_Param->e_WithNoiseModel) 
      po_TabNoiseModel = new MR1DNoiseModel [po_Param->i_NbSpectre];
   else 
      po_TabNoiseModel = (MR1DNoiseModel*)NULL;
      
}


pca_Result::~pca_Result () {
  delete [] po_TabMr1d;
  if (po_Param->e_WithNoiseModel) delete [] po_TabNoiseModel;
}


void pca_Result::res_FirstCompute () {

  // call number
  LoopNumber = ++LoopNumber;

  // create all multiresol
  if (po_Param->e_TraceParamIn)
     cout << "Multiresol transform of data ..." << endl;
  for (int k=0; k<po_Param->i_NbSpectre; k++)
    po_TabMr1d[k].transform (po_Param->po_TabSpectre[k]);
    
  if (po_Param->e_WriteTransf) {
     char Mr1dRecons[256];
     for (int k=0; k < po_Param->i_NbSpectre; k++) {
        sprintf(Mr1dRecons, "mr1d_Transf_loop_%d_spec_%d", LoopNumber, k+1);
        fits_write_fltarr (Mr1dRecons, po_TabMr1d[k].image());
     } 
  }

  // create all noise model
  if (po_Param->e_WithNoiseModel) {
     res_InitNoiseModelMr1d();
     //res_ComputeCommonSupport();
  }
 
  // calculate the correlation matrix and its eigen values
  if (po_Param->e_TraceParamIn)
      cout << "Principal component analysis ..." << endl;  
  if (po_Param->e_WithNoiseModel) 
     o_Mr1dPca.compute (po_TabMr1d, (MR1DNoiseModel*)NULL,
                        po_Param->e_NormCorrelMatrix,
                        po_Param->e_CorrelComp,
                        DEF_CORREL_NOISE, 
		        po_Param->po_ImportCorMat); 
  else
     o_Mr1dPca.compute(po_TabMr1d, (MR1DNoiseModel*)NULL);			
 
  // print the result to the standard output
  if (po_Param->e_Verbose) o_Mr1dPca.print();
  
  // Write signal in
  if (po_Param->e_WriteTransf == True) {
      char SigIn[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(SigIn, "mr1d_SignalIn_loop_%d_spec_%d", LoopNumber, k+1);
	 fits_write_fltarr (SigIn, po_Param->po_TabSpectre[k]);
      }
   }   
       
   // create all multiresol
   //if (po_Param->e_TraceParamIn)
      //cout << "Multiresol transform of data ..." << endl
   //   for (int k=0; k<po_Param->i_NbSpectre; k++) 
   //      po_TabMr1d[k].transform (po_Param->po_TabSpectre[k]);  
  
  // calculate the eigen vector spectres 
  if (po_Param->e_TraceParamIn)
     cout << "Pca transform of data ..." << endl; 
  o_Mr1dPca.pca_TransfSignal (po_TabMr1d, po_TabMr1d);
  
  // trace.... 
  if (po_Param->e_WriteTransf) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbSpectre; k++) {
        sprintf(NameEigen, "pca_Transf_loop_%d_spec_%d", LoopNumber, k+1);
        fits_write_fltarr (NameEigen, po_TabMr1d[k].image());
     }
  }  
  
  // threshold in pca space
  if (po_Param->e_WithNoiseModel) {
     res_InitNoiseModelPca ();
     if (po_Param->e_OptimSoftThreshold) {
        res_SoftInitAlloc ();
        res_SoftThreshold();
     } else
        o_Mr1dPca.pca_Thresold (po_TabMr1d, po_TabNoiseModel);
  }
   
  if (po_Param->e_WriteTransf) {
      char Pca_ThreshTransf[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(Pca_ThreshTransf, "pca_Threshold_Transf_loop_%d_spec_%d", LoopNumber, k+1);
         fits_write_fltarr (Pca_ThreshTransf, po_TabMr1d[k].image());
      }
   }      
}
    

void pca_Result::res_CurrentCompute () {
   
   
   // call number
   LoopNumber = ++LoopNumber;   
   
   // Write signal in
   if (po_Param->e_WriteTransf == True) {
      char SigIn[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(SigIn, "mr1d_SignalIn_loop_%d_spec_%d", LoopNumber, k+1);
	 fits_write_fltarr (SigIn, po_Param->po_TabSpectre[k]);
      }
   }   
       
   // create all multiresol
   if (po_Param->e_TraceParamIn)
      cout << "Multiresol transform of data ..." << endl;
   for (int k=0; k<po_Param->i_NbSpectre; k++) 
      po_TabMr1d[k].transform (po_Param->po_TabSpectre[k]);

   if (po_Param->e_WriteTransf) {
      char Mr1dRecons[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(Mr1dRecons, "mr1d_Transf_loop_%d_spec_%d", LoopNumber, k+1);
         fits_write_fltarr (Mr1dRecons, po_TabMr1d[k].image());
      } 
   }   
   
   // calculate the eigen vector spectres 
   if (po_Param->e_TraceParamIn)
      cout << "Pca transform of data ..." << endl; 
   o_Mr1dPca.pca_TransfSignal (po_TabMr1d, po_TabMr1d);   
  
   if (po_Param->e_WriteTransf) {
      char NameEigen[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(NameEigen, "pca_Transf_loop_%d_spec_%d", LoopNumber, k+1);
         fits_write_fltarr (NameEigen, po_TabMr1d[k].image());
      }
   }  
  
   if (po_Param->e_WithNoiseModel) {
      if (po_Param->e_OptimSoftThreshold)
         res_SoftThreshold();
      else
         o_Mr1dPca.pca_Thresold (po_TabMr1d, po_TabNoiseModel);
   }
   
   if (po_Param->e_WriteTransf) {
      char Pca_ThreshTransf[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(Pca_ThreshTransf, "pca_Threshold_Transf_loop_%d_spec_%d", LoopNumber, k+1);
         fits_write_fltarr (Pca_ThreshTransf, po_TabMr1d[k].image());
      }
   }    
}


void pca_Result::res_Recons () {

   if (po_Param->e_Verbose) cout << "Reconstruction process ..." << endl;
   int ai_NbEigen=0;
   for (int i=0; i<po_Param->ti_NbrUseEigen[0]; i++) //same value for all bands
      if (po_Param->tti_TabKill[0][i] == 0) ai_NbEigen++;
   if (ai_NbEigen != po_Param->i_NbSpectre) {
      if (po_Param->e_Verbose) cout << "  Number of eigen vector used for the reconstruction : " ;
      if (po_Param->e_Verbose) cout << "  Ne = " << ai_NbEigen << " (in all band)" << endl;
   }
    
   // apply the inverse reconstruction from a subset of eigen vectors
   if (po_Param->e_Verbose) cout << "  Inverse PCA Transform ..." << endl;
   o_Mr1dPca.invsubtransform (po_TabMr1d, po_TabMr1d,
                              po_Param->ti_NbrUseEigen,
                              po_Param->tti_TabKill);
			      
   // accept wavelet coef only if sigma is .....
   if (po_Param->e_WithNoiseModel && po_Param->e_TresholdInMrSpace) {
      int Ind;
      for (int k=0; k<po_Param->i_NbSpectre; k++)
         for (int s=0; s<po_Param->i_NbBand-1;s++)
            for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++) {
	       if (one_level_per_pos(po_Param->e_TypeNoise) == True) 
	            Ind = po_TabNoiseModel[k].get_index(s, i);
               else Ind = s;
	       if ((po_TabMr1d[k])(s,i) < MrSigmaLevel(k,Ind)) 
	          (po_TabMr1d[k])(s,i) = 0.0;  
	    }
   }
   
   // smoth constraint
   if (po_Param->e_WithNoiseModel && po_Param->e_OptimSoftThreshold) {
      for (int k=0; k<po_Param->i_NbSpectre; k++)
         for (int s=0; s<po_Param->i_NbBand-1;s++)
            for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++) {
	       //s serach for sigma at sacle s and position i...
	       // => (MrSigmaLevel(k,Ind))   
	       int Ind;    
	       if (one_level_per_pos(po_Param->e_TypeNoise) == True) 
	            Ind = po_TabNoiseModel[k].get_index(s, i);
               else Ind = s;
	       // soft threshold
	       float SoftThres =   (po_Param->i_NbIter+1 - LoopNumber)
	                         / po_Param->i_NbIter * MrSigmaLevel(k,Ind); 
	       (po_TabMr1d[k])(s,i)  = soft_threshold ((po_TabMr1d[k])(s,i),
		                                       SoftThres);	    
	    }   
   }
   
   //remove negative detection in recons mr (by inverse PCA)
   if (po_Param->e_WithNoiseModel && po_Param->e_OnlyPositivDetect) {
      for (int k=0; k<po_Param->i_NbSpectre; k++)
         for (int s=0; s<po_Param->i_NbBand-1;s++)
            for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++)
	       if ((po_TabMr1d[k])(s,i) < 0.0) (po_TabMr1d[k])(s,i)=0.0;
   }
   
   if (po_Param->e_WriteTransf) {
      char Mr1dRecons[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(Mr1dRecons, "mr1d_pca_InvTransf_loop_%d_spec_%d", LoopNumber, k+1);
         fits_write_fltarr (Mr1dRecons, po_TabMr1d[k].image());
      } 
   }    
   
    		     		      
   // inv transform Mr1d
   if ((po_Param->e_Verbose) && (po_Param->e_TraceParamIn))
      cout << "  Inverse Multiresol Transform ..." << endl;
   for (int k=0; k<po_Param->i_NbSpectre; k++) {
      if (po_Param->e_UseReconsAdjoint) {
         if (po_Param->e_RemoveLastScale) 
            po_TabMr1d[k].rec_adjoint((po_Param->po_TabSpectre)[k], False);
	 else
	    po_TabMr1d[k].rec_adjoint((po_Param->po_TabSpectre)[k]);
      } else {
         if (po_Param->e_RemoveLastScale)
            for (int i=0;i<po_TabMr1d[k].size_scale_np(po_TabMr1d[k].nbr_band()-1);i++) 
	       (po_TabMr1d[k])(po_TabMr1d[k].nbr_band()-1,i)=0.0;
         po_TabMr1d[k].recons((po_Param->po_TabSpectre)[k]); 
      }
   }   
      
   // Write Mr1d recons
   if (po_Param->e_WriteTransf == True) {
      char Mr1dRecons[256];
      for (int k=0; k < po_Param->i_NbSpectre; k++) {
         sprintf(Mr1dRecons, "mr1d_InvTransf_loop_%d_spec_%d", LoopNumber, k+1);
	 fits_write_fltarr (Mr1dRecons, po_Param->po_TabSpectre[k]);
      }
   }    
}



void pca_Result::res_Write () {

  if (po_Param->e_TraceParamIn)
     cout << "Write result ..." << endl;
  for (int k=0;k<po_Param->i_NbSpectre;k++)
    for (int i=0;i<po_Param->i_NbPts;i++)
      o_2dDataOut(i,k) = po_Param->po_TabSpectre[k](i);

  // write the data
  po_Param->ts_Header.origin = po_Param->tc_Cmd;
  fits_write_fltarr (po_Param->tc_NameOut, o_2dDataOut,
                     &po_Param->ts_Header);
  
} 




void pca_Result::res_InitNoiseModelMr1d () {

   // noise model class initialization
   if (po_Param->e_TraceParamIn)
      cout << "Create mr1d noise model ..." << endl;
   for (int k=0; k<po_Param->i_NbSpectre; k++) {

      po_TabNoiseModel[k].alloc (po_Param->e_TypeNoise, po_Param->i_NbPts, 
                                 po_Param->i_NbPlan, po_Param->e_Transform);

      if (po_Param->f_NoiseSpectre > FLOAT_EPSILON) 
         po_TabNoiseModel[k].SigmaNoise = po_Param->f_NoiseSpectre;

      if (po_Param->f_NSigmaMr1d >= 0) {
         for (int i=0; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaMr1d;
      } else {
         po_TabNoiseModel[k].NSigma[0] = 4.0;
         for (int i=1; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = 3.0;     
      }	    
      
      po_TabNoiseModel[k].FirstDectectScale = po_Param->i_FirstScale;
      //po_TabNoiseModel[k].OnlyPositivDetect = po_Param->e_OnlyPositivDetect;
      po_TabNoiseModel[k].SupIsol = po_Param->e_SupIsolPixel;
      
      po_TabNoiseModel[k].model (po_Param->po_TabSpectre[k], 
                                 po_TabMr1d[k] );
      
      if (po_Param->e_WriteSupportMr1d) {
         char NameSupport[256];
         sprintf(NameSupport, "mr1d_Support_%d",k+1);     
         if (po_TabNoiseModel[k].set_trans() == TRANS1_PAVE) {
            fltarray Support(po_Param->i_NbPts);
	    Support.init (0.);
            for (int s=0; s<po_Param->i_NbBand-1;s++)
               for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++) {
	          if (po_TabNoiseModel[k].signif((po_TabMr1d[k])(s,i),s,i) != 0){ 
	             Support(i) += 2^(po_Param->i_NbBand-s);
                  }
	       }
	    fits_write_fltarr (NameSupport, Support);       
         } else {
            po_TabNoiseModel[k].write_support_mr (NameSupport);
	 }
      }
   }
   
   // alloc MrSigmaLevel
   if (po_Param->e_WithNoiseModel) {
      if (one_level_per_pos(po_Param->e_TypeNoise) == True) 
         MrSigmaLevel.alloc (po_Param->i_NbSpectre, po_TabNoiseModel->get_size());
      else MrSigmaLevel.alloc (po_Param->i_NbSpectre, po_Param->i_NbBand);
   }      
      
   for (int k=0; k<po_Param->i_NbSpectre; k++) {
   // memorize Tab_level (sigma level) for future use in recons.....
      if (one_level_per_pos(po_Param->e_TypeNoise) == True) 
         for (int i=0;i<po_TabNoiseModel[k].get_size();i++) 
	    MrSigmaLevel(k,i) = po_TabNoiseModel[k].sigma(i);
      else
         for (int s=0;s<po_Param->i_NbBand;s++) 
	    MrSigmaLevel(k,s) = po_TabNoiseModel[k].sigma(s,0);
   }
}

void pca_Result::res_InitNoiseModelPca () {

   // noise model class initialization
   if (po_Param->e_TraceParamIn)
      cout << "Create pca noise model ..." << endl;
   int k;
   
//  for (k=0; k<po_Param->i_NbSpectre; k++) {
//    cout << "Before Spectre " << k << endl;
//    cout << "  ";
//    for (int i=0;i<po_Param->i_NbPlan-1; i++) 
//       cout << po_TabNoiseModel[k].sigma(i,0) << ",";
//    cout << endl;
// } 
   
   o_Mr1dPca.pca_ComputeNoise (po_TabMr1d, po_TabNoiseModel, po_TabNoiseModel);
   
   for (k=0; k<po_Param->i_NbSpectre; k++) {
   
      //if (o_Param.e_DilateSupportPCA) po_TabNoiseModel[k].DilateSupport = True;
      if (po_Param->f_NSigmaPCA >= 0)
         for (int i=0; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaPCA;
      else {
         po_TabNoiseModel[k].NSigma[0] = 4;
         for (int i=1; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = 3;         
      }
     
      po_TabNoiseModel[k].FirstDectectScale = po_Param->i_FirstScale;
      //po_TabNoiseModel[k].OnlyPositivDetect = po_Param->e_OnlyPositivDetect;
      po_TabNoiseModel[k].SupIsol = po_Param->e_SupIsolPixel;
      
      po_TabNoiseModel[k].set_support (po_TabMr1d[k]);
     
      
      //if (o_Param.e_DestroyRings) res_DestroyRingsInPcaSup ();
   }
    
   for (k=0; k<po_Param->i_NbSpectre; k++) {
      // cout << "After Spectre " << k << endl;
      // cout << "  ";
      // for (int i=0;i<po_Param->i_NbPlan-1; i++) 
      //    cout << po_TabNoiseModel[k].sigma(i,0) << ",";
      // cout << endl;
      if (po_Param->e_WriteSupportPCA) {
         char NameSupport[256];
         sprintf(NameSupport, "pca1d_Support_%d",k+1);     
         if (po_TabNoiseModel[k].set_trans() == TRANS1_PAVE) {
            fltarray Support(po_Param->i_NbPts);
	    Support.init (0.);
            for (int s=0; s<po_Param->i_NbBand-1;s++)
               for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++) {
	          if (po_TabNoiseModel[k].signif((po_TabMr1d[k])(s,i),s,i) != 0) 
	             Support(i) += 2^(po_Param->i_NbBand-s);
	       }
	    fits_write_fltarr (NameSupport, Support);       
         } else {
            po_TabNoiseModel[k].write_support_mr (NameSupport);
	 }
      }
   }
}

to_Param1d* pca_Result::res_GetpParam () {return po_Param;}


void pca_Result::res_ComputeCommonSupport () {
      
   fltarray ComSupport(po_Param->i_NbBand, po_Param->i_NbPts);
   ComSupport.init(0.);
   for (int k=0; k<po_Param->i_NbSpectre; k++) {
   
      if (po_TabNoiseModel[k].set_trans() == TRANS1_PAVE) {
         
	 for (int s=0; s<po_Param->i_NbBand-1;s++)
            for (int i=0;i<po_TabMr1d[k].size_scale_np(s);i++) {
	       if (po_TabNoiseModel[k].signif((po_TabMr1d[k])(s,i),s,i)==1) {
	          ComSupport(s,i) = 1;
	       }  
            }      
      } else {
         cout << " Not yet implemented ... " << endl;
         exit(-1);
      }         
   }
   
   if (po_Param->e_WriteSupportMr1d) {
      fits_write_fltarr ("mr1d_ComSupport", ComSupport);
   }
   
   //reduce support
   for (int s=0; s<po_Param->i_NbBand-1;s++) {
      
      int CurrentIndice=0;
      // throught the band s
      while (CurrentIndice < po_Param->i_NbPts) {
      
         // is it support ?
	 if (ComSupport(s,CurrentIndice) == 1) {
	 
	    //compute length of support
	    int BegIndSup = CurrentIndice;
	    int Compt = 1;
	    while (ComSupport(s,CurrentIndice++) == 1) {
	       Compt++;
	    }
	    
	    //if Length too short => no support
	    if (Compt < 20)
	       for (int i=BegIndSup; i<BegIndSup+Compt; i++)
	          ComSupport(s,i) = 0;
	    
	 // else increment CurrentIndice
	 } else {
	    CurrentIndice++;
	 }
      }
   }
       
   if (po_Param->e_WriteSupportMr1d) {
      fits_write_fltarr ("mr1d_RedComSupport", ComSupport);
   } 
   
   //compute segment number
   fltarray NbSegment(po_Param->i_NbBand);
   fltarray BegSegment(po_Param->i_NbPts,po_Param->i_NbBand);
   fltarray EndSegment(po_Param->i_NbPts,po_Param->i_NbBand);
   
   
   // throught the band s
   for (int s=0; s<po_Param->i_NbBand-1;s++) {
   
      int CurrentIndice=0;
      BegSegment(0,s)=0;
      NbSegment(s)=-1;
      // throught the band s
      while (CurrentIndice < po_Param->i_NbPts-1) {
      
         // new segment and begin
         NbSegment(s) =  NbSegment(s) + 1;
         BegSegment(NbSegment(s),s)=CurrentIndice;
	 
         while (    ComSupport(s,CurrentIndice+1) 
	         == ComSupport(s,CurrentIndice)) {
	    CurrentIndice++;
	 }
	 
	 // end segment
	 EndSegment(NbSegment(s),s) = CurrentIndice;    
	 CurrentIndice++; 
      }
   }
   
   if (po_Param->e_WriteSupportMr1d) {
      fits_write_fltarr ("mr1d_segment_nb.fits", NbSegment);
      fits_write_fltarr ("mr1d_beg_segment.fits", BegSegment);
      fits_write_fltarr ("mr1d_end_segment.fits", EndSegment);
   }   
}


/******************************************************************************
classe SoftIter
******************************************************************************/



void pca_Result::res_SoftInitAlloc () {
  
   if (po_Param->e_Transform == TO1_PAVE_MORLET) {
      
      cout << "Not Yet impemented..." << endl;
      exit(-1);  
      
      /*FirstPcaCoef.alloc (po_Param->i_NbSpectre, 
                          po_Param->i_NbPts, 
			  2*po_Param->i_NbPlan);
      ResidualPcaCoef.alloc (po_Param->i_NbSpectre, 
                             po_Param->i_NbPts, 
		             2*po_Param->i_NbPlan);
      for (int k=0; k<po_Param->i_NbSpectre; k++)
         for (int s=0; s<2*po_Param->i_NbPlan; s++)
	    for (int i=0; i<po_Param->i_NbPts; i++) {
	       FirstPcaCoef(k,i,s) = (po_TabMr1d[k])(s,i);
	       (po_TabMr1d[k])(s,i) = 0.;
	    }
      */
   } else 
   if (   (which_set_is_trans1d(po_Param->e_Transform) == TRANS1_MALLAT)
        ||(which_set_is_trans1d(po_Param->e_Transform) == TRANS1_WP_MALLAT)) {
       
      cout << "Not Yet impemented..." << endl;
      exit(-1); 
         
      /*FirstPcaCoef.alloc (po_Param->i_NbSpectre, po_Param->i_NbPts); 
      ResidualPcaCoef.alloc (po_Param->i_NbSpectre,po_Param->i_NbPts); 		  
      for (int k=0; k<po_Param->i_NbSpectre; k++) {
         int j=0;
         for (int s=0; s<po_Param->i_NbPlan; s++)
	    for (int i=0; i<po_TabMr1d[k].size_scale_np(s); i++) {
	       FirstPcaCoef(k,j++) = (po_TabMr1d[k])(s,i);
	       (po_TabMr1d[k])(s,i) = 0.;
	    }
      }
      */          
   } else {
   
      FirstPcaCoef.alloc (po_Param->i_NbSpectre,
                          po_Param->i_NbPts, 
                          po_Param->i_NbBand);
      ResidualPcaCoef.alloc (po_Param->i_NbSpectre,
                             po_Param->i_NbPts, 
                             po_Param->i_NbBand);	   
      for (int k=0; k<po_Param->i_NbSpectre; k++)
         for (int s=0; s<po_Param->i_NbBand; s++)
	    for (int i=0; i<po_Param->i_NbPts; i++) {
	       FirstPcaCoef(k,i,s) = (po_TabMr1d[k])(s,i);
	       (po_TabMr1d[k])(s,i) = 0.;
	    }
   } 
}

void pca_Result::res_SoftThreshold () {

   // compute coef for lamba at current iter 
   // Lambda = Lambda * (1-currrent iter)/NbIter
   //-------------------------------------------
   float CoefLambda =   (1. - (LoopNumber-1.)/po_Param->i_NbIter);
   if (po_Param->e_TraceParamIn) cout << "Coef Lambda : " << CoefLambda << endl;
   
   // threshold
   //----------
   if (po_Param->e_Transform == TO1_PAVE_MORLET) {
   
      cout << "Not Yet impemented..." << endl;
      exit(-1);
    
      /*for (int k=0; k<po_Param->i_NbSpectre; k++) {
         CoefLambda = CoefLambda * po_TabNoiseModel[0].NSigma[k];
         for (int s=0; s<2*po_Param->i_NbPlan-1; s++)
	    for (int i=0; i<po_Param->i_NbPts; i++) {
	    
	       float LocalSigma = po_TabNoiseModel[k].sigma(s,i);
              
	       // compute residual on PCA coef
               //----------------------------- 	    
	       ResidualPcaCoef(k,i,s) = 
	          FirstPcaCoef(k,i,s) - (po_TabMr1d[k])(s,i);
	    
	       // first rule, if FirstPcaCoef is significatif and
               // abs(ResidualPcaCoef) > eps sigma then coef = FirstPcaCoef
               //----------------------------------------------------------
	       if (    (po_TabNoiseModel[k].signif (FirstPcaCoef(k,i,s),s,i))             
	            && (ABS(ResidualPcaCoef(k,i,s)) > LocalSigma/2.)) {          
		  (po_TabMr1d[k])(s,i) = FirstPcaCoef(k,i,s);
	       } else {
	          (po_TabMr1d[k])(s,i) = 0.0;
	       }
	       
	       // second rule, soft threshold at sigma*CoefLambda
	       (po_TabMr1d[k])(s,i) = soft_threshold ((po_TabMr1d[k])(s,i),
	                                              CoefLambda*LocalSigma);
	    }
      } */    	
   } else
   if (   (which_set_is_trans1d(po_Param->e_Transform) == TRANS1_MALLAT)
        ||(which_set_is_trans1d(po_Param->e_Transform) == TRANS1_WP_MALLAT)) {
      
      cout << "Not Yet impemented..." << endl;
      exit(-1);
           
      /*for (int k=0; k<po_Param->i_NbSpectre; k++) {
         int j=0;
	 CoefLambda = CoefLambda * po_TabNoiseModel[0].NSigma[k];
         for (int s=0; s<po_Param->i_NbPlan; s++)
	    for (int i=0; i<po_TabMr1d[k].size_scale_np(s); i++) {
	       
	       float LocalSigma = po_TabNoiseModel[k].sigma(s,i);
	       
	       // compute residual on PCA coef
               //----------------------------- 	    
	       ResidualPcaCoef(k,j) = 
	          FirstPcaCoef(k,j) - (po_TabMr1d[k])(s,i); 
	     
	       // first rule, if FirstPcaCoef is significatif and
               // abs(ResidualPcaCoef) > eps sigma then coef = FirstPcaCoef
               //----------------------------------------------------------
	       if (po_TabNoiseModel[k].signif (FirstPcaCoef(k,i,s),s,i)) {
	                    
	          // Pix is signif, if abs(res) > sig/2 => sol = init 
		  if ((ABS(ResidualPcaCoef(k,j)) > LocalSigma/2.)) {          
		     (po_TabMr1d[k])(s,i) = FirstPcaCoef(k,j);
		  }
		  
	       } else {
	       // pix is not signif, abs(sol) < k*sigma
	          float kSigma =   po_TabMr1d[k].Nsigma[s] * LocalSigma; 
		  
	       }
	       
	       // second rule, soft threshold at sigma*CoefLambda
	       (po_TabMr1d[k])(s,i) = soft_threshold ((po_TabMr1d[k])(s,i),
	                                              CoefLambda*LocalSigma);	     
	       j++;
	    }   
      }*/
                
   } else {
      
      for (int k=0; k<po_Param->i_NbSpectre; k++) {
         for (int s=0; s<po_Param->i_NbPlan; s++) {
	    CoefLambda = CoefLambda * po_TabNoiseModel[k].NSigma[s];
	    for (int i=0; i<po_Param->i_NbPts; i++) {
	        
	       float LocalSigma = po_TabNoiseModel[k].sigma(s,i);
	       
	       // compute residual on PCA coef
               //----------------------------- 	    
	       ResidualPcaCoef(k,i,s) =           
	          FirstPcaCoef(k,i,s) - (po_TabMr1d[k])(s,i);
	       
	       // first rule, if FirstPcaCoef is significatif and
               // abs(ResidualPcaCoef) > eps sigma then coef = FirstPcaCoef
               //----------------------------------------------------------
	       if (po_TabNoiseModel[k].signif (FirstPcaCoef(k,i,s),s,i)) {
	              
		  // Pix is signif, if abs(res) > sig/2 => sol = init      
	          if ((ABS(ResidualPcaCoef(k,i,s)) > LocalSigma/2.)) { 
		     (po_TabMr1d[k])(s,i) = FirstPcaCoef(k,i,s);
                  }
		  
	       } else {
	          // pix is not signif, abs(sol) < k*sigma
		  /*float kSigma =   po_TabNoiseModel[k].NSigma[s] * LocalSigma;
	          if (    ((po_TabMr1d[k])(s,i) >= 0) 
		       && ((po_TabMr1d[k])(s,i) > kSigma)) 
		     (po_TabMr1d[k])(s,i) = kSigma;
		  if (    ((po_TabMr1d[k])(s,i) < 0)
		       && ((po_TabMr1d[k])(s,i) < kSigma)) 		  
		     (po_TabMr1d[k])(s,i) = -kSigma;*/
	       }
	       
	       // second rule, soft threshold at sigma*CoefLambda
	       // (po_TabMr1d[k])(s,i) = soft_threshold ((po_TabMr1d[k])(s,i),
	       //                                       CoefLambda*LocalSigma);		  
	    }
	 }
      }
   }   
}
   

/******************************************************************************
classe to_Iter
******************************************************************************/

Bool pca_Iter::iter_End () {
   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > po_Result->po_Param->i_NbIter) { 
  
   
      for (int k=0;k<po_Result->po_Param->i_NbSpectre;k++) 
         for (int i=0;i<po_Result->po_Param->i_NbPts;i++) {
               po_Result->po_Param->po_TabSpectre[k](i) = o_2dResidu(i,k); 
	       if (   po_Result->po_Param->e_PositivSpectre 
	           && po_Result->po_Param->po_TabSpectre[k](i)<0.) 
                  po_Result->po_Param->po_TabSpectre[k](i) = 0.0;
	       if (   po_Result->po_Param->e_MaxSpectre 
	           && po_Result->po_Param->po_TabSpectre[k](i)>MAX_SPEC_VALUE) 
                  po_Result->po_Param->po_TabSpectre[k](i) = MAX_SPEC_VALUE;		  
		  
	    } 
	            
      return False;
   } else return True;
}

void pca_Iter::iter_Residu () {

   for (int k=0;k<po_Result->po_Param->i_NbSpectre;k++) {
      for (int i=0;i<po_Result->po_Param->i_NbPts;i++) {
	 
	    if (   po_Result->po_Param->e_MaxSpectre 
	        && po_Result->po_Param->po_TabSpectre[k](i)>MAX_SPEC_VALUE) 
               po_Result->po_Param->po_TabSpectre[k](i) = MAX_SPEC_VALUE;
	     
		  		          
            // I(i) = I(i) + E(i)
            o_2dResidu(i,k) += po_Result->po_Param->po_TabSpectre[k](i);
	    
	    // positiv constraint
	    if (    po_Result->po_Param->e_PositivSpectre 
	         && o_2dResidu(i,k) < 0.) 
               o_2dResidu(i,k) = 0.0; 
	                
	    // E(i) = Einit - I(i)
	    po_Result->po_Param->po_TabSpectre[k](i) = o_2dInit(i,k) - o_2dResidu(i,k);
	    
      }
      
      // cout << "energy[" << k << "]=" << po_Result->po_Param->po_TabSpectre[k].energy() << endl;
      // write  sol at current iteration	
      if (po_Result->po_Param->e_WriteTransf == True) {
          char Mr1dRecons[256];
          sprintf(Mr1dRecons, "current sol_loop_%d_spec_%d",NbIter , k+1);
	  fltarray prov(po_Result->po_Param->i_NbPts);
	  for (int l=0;l<po_Result->po_Param->i_NbPts;l++)
	     prov(l) = o_2dResidu(l,k);
	  fits_write_fltarr (Mr1dRecons, prov);
      }  
               
   }      
}

to_Result1d* pca_Iter::iter_getpResult() {return po_Result;}


/******************************************************************************
classe to_Iter
******************************************************************************/

Bool pca_SoftIter::iter_End () {   NbIter++;
      
   if (NbIter > po_Result->po_Param->i_NbIter) { 
  
      for (int k=0;k<po_Result->po_Param->i_NbSpectre;k++) 
         for (int i=0;i<po_Result->po_Param->i_NbPts;i++) {
						      	 
	       if (   po_Result->po_Param->e_PositivSpectre 
	           && po_Result->po_Param->po_TabSpectre[k](i)<0.) 
                  po_Result->po_Param->po_TabSpectre[k](i) = 0.0;
	       
	       if (   po_Result->po_Param->e_MaxSpectre 
	           && po_Result->po_Param->po_TabSpectre[k](i)>MAX_SPEC_VALUE) 
                  po_Result->po_Param->po_TabSpectre[k](i) = MAX_SPEC_VALUE;		  
		  
	    } 
	            
      return False;
   } else return True;
}
   
void pca_SoftIter::iter_ActionOnSol() {   

   for (int k=0;k<po_Result->po_Param->i_NbSpectre;k++) {
      for (int i=0;i<po_Result->po_Param->i_NbPts;i++) {
		  	 
	 if (   po_Result->po_Param->e_PositivSpectre 
	     && po_Result->po_Param->po_TabSpectre[k](i)<0.) 
            po_Result->po_Param->po_TabSpectre[k](i) = 0.0;
  
         if (   po_Result->po_Param->e_MaxSpectre 
	     && po_Result->po_Param->po_TabSpectre[k](i)>MAX_SPEC_VALUE) 
            po_Result->po_Param->po_TabSpectre[k](i) = MAX_SPEC_VALUE;
      }    
      
      cout << "energy[" << k << "]=" << po_Result->po_Param->po_TabSpectre[k].energy() << endl;
      // write  sol at current iteration	
      if (po_Result->po_Param->e_WriteTransf == True) {
          char Mr1dRecons[256];
          sprintf(Mr1dRecons, "current sol_loop_%d_spec_%d",NbIter , k+1);
          fits_write_fltarr (Mr1dRecons, po_Result->po_Param->po_TabSpectre[k]);
      }
   }                
}    

to_Result1d* pca_SoftIter::iter_getpResult() {return po_Result;}






















