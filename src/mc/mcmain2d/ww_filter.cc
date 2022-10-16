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
**    File:  mc_wthfilter.cc
**
*******************************************************************************/


#include "ww_filter.h"
//#include "mc_usage.h"

#define DEF_SIGMA_WTH 3
#define DEF_SIGMA_MR2D 3

#define MAX_ITER 100
#define NBR_FILTER 4

extern char *OptArg;
extern int   OptInd;

int main(int argc, char *argv[]) {

   extern softinfo Soft;
   Soft.mr3();
   lm_check(LIC_MR3);
   
  // licence manager
  //lm2();

  // init local var
  wth_Param o_Param;

  // Read param IN
  o_Param (argc, argv);
      
  // init Result and Iter classes
  wth_Result o_Result (&o_Param);   
  wth_Iter o_Iter (&o_Result); 
     
  for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {} 

  // write result
  o_Result.res_Write ();
        
  exit(0);
}





/******************************************************************************
classe wth_Param
******************************************************************************/
 
wth_Param::wth_Param () {

  //default init of read var
  e_Transform =           DEFAULT_TRANSFORM;     
  e_Norm =                NORM_L1; 
  e_SBFilter =            F_MALLAT_7_9;
  po_FAS =                (FilterAnaSynt*)NULL;  
  i_NbPlan =              DEFAULT_NBR_SCALE;
  i_NbIter =              1;
  e_UseReconsAdjoint =    False;
  e_Verbose =             False;
  e_WriteTransf =         False;
  e_TraceParamIn =        False;
  
  // noise model
  e_WithNoiseModel =      True;       // always use noise model, def Gaussian
  e_TypeNoise =           DEFAULT_STAT_NOISE;
  f_NoiseIma =            0.;
  f_Gain =                0.;
  f_SigmaGauss =          0.;
  f_MeanGauss =           0.;  
  f_NSigmaMr2d  =         DEF_SIGMA_MR2D; 
  f_NSigmaWth =           DEF_SIGMA_WTH; 
  e_WriteSupportMr2d =    False;
  e_WriteSupportWth =     False;  
  e_LiftingTrans =        DEF_LIFT;
  e_Norm  =               NORM_L1;
  e_SBFilter =            F_MALLAT_7_9;  
  KillLastScale = False;
  e_PositivImag =         False;
}

void wth_Param::operator () (int argc, char *argv[]) {
    
  int c,i;
  Bool OptL = False, Optf = False;
  while ((c = GetOpt(argc,argv,"Kt:m:n:g:c:s:S:l:f:C:T:LrRWkPav")) != -1) {
	
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
      readfloat (f_NSigmaWth, "Error: bad NSigma: %s\n");
      break;	  
      
    case 'l':
      readint (i, "Error: bad lifting method: %s\n");
      if (testborn(1, NBR_LIFT, i, "Error: bad type of noise: %s\n"))
	e_LiftingTrans = (type_lift) (i-1);
      break;   
      
    case 'f':
      readint (i, "Error: bad filter: %s\n");
      if (testborn(1, NBR_FILTER, i, "Error: bad type of filter: %s\n")) {
	e_SBFilter = (type_sb_filter) (i);		
	Optf = True;
      }
      break;  
   case 'T':
      Optf = True;
      e_SBFilter = get_filter_bank(OptArg);
      break;     
    case 'P':
      e_PositivImag = True;break;      
    case 'K':
      KillLastScale = True;break;      
    case 'L':
      e_Norm = NORM_L2; OptL = True; break;
    case 'r':
      e_WriteSupportMr2d = True;break;
    case 'R':
      e_WriteSupportWth = True;break;                  
    case 'a':
      e_UseReconsAdjoint = True;break;
    case 'W':
      e_WriteTransf = True;break;
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
  if (e_TraceParamIn) trace ();
}

void wth_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options input_image output_image\n", argv[0]);
  manline(); 

  transform_usage (e_Transform);      manline();
  //lifting_usage (e_LiftingTrans);     manline();
  sb_usage (e_SBFilter);              manline();  
  nbr_scale_usage (i_NbPlan);         manline();
  //nbr_iteration (i_NbIter);           manline();
  //noise
  //wk_noise_usage ();                  manline();
  gauss_usage ();                     manline();
  //ccd_usage ();                       manline();
  nsigma_mr2d (DEF_SIGMA_MR2D);       manline();
  nsigma_wth (DEF_SIGMA_WTH);         manline(); 
  //write_support_mr2d ();              manline();
  //write_support_wth ();               manline();
  positiv_imag_constraint();          manline();
     
  verbose ();                         manline();
  exit(-1);
}

void wth_Param::read_image (char *argv[]) {

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
}


void wth_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransform(e_Transform) << endl;
   if (e_WithNoiseModel) cout << " -- [-m:] : with noise model : " << StringNoise(e_TypeNoise) << endl;   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   cout << " -- [-i:] : number of iteration : " << i_NbIter << endl;
   if (f_NoiseIma != 0.) cout << " -- [-g:] : gauss noise : " << f_NoiseIma << endl;
   if (f_Gain != 0.) cout << " -- [-c:] : gain : " << f_Gain << ", sigma : " << f_SigmaGauss << ", mean : " << f_MeanGauss << endl;
   if (e_WithNoiseModel) {
      cout << " -- [-s:] : multiresol nsigma : " << f_NSigmaMr2d << endl;
      cout << " -- [-S:] : wth nsigma : " << f_NSigmaWth << endl;
   } 
   if (e_Transform == 14) cout << " -- [-f:] : filter : " << (int)e_SBFilter << endl;
   if (e_Transform == 24) cout << " -- [-l:] : lifting : " << e_LiftingTrans-1 << endl;
   if (e_UseReconsAdjoint) cout << " -- [-a]  : use rec_adjoint in recons process" << endl;
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect images" << endl;
   if (e_WriteSupportMr2d) cout << " -- [-r]  : write Multiresol support" << endl;
   if (e_WriteSupportWth) cout << " -- [-R]  : write WTH support" << endl;
   if (e_Norm) cout << " -- [-L]  : Norm = L2" << endl;   
   if (e_PositivImag) cout << " -- [-P]  : positivity constraint" << endl;
   if (e_Verbose) cout << " -- [-v]  : verbose" << endl;
}



/******************************************************************************
classe wth_Result
******************************************************************************/

wth_Result::wth_Result (wth_Param* ppo_Param) {

   po_Param = ppo_Param;
   // allocate fltarray for write result
   o_3dDataOut.alloc (po_Param->i_NbCol, po_Param->i_NbLin, po_Param->i_NbImage);
   // allocate the memory for all multiresol
   po_TabMr2d = new MultiResol [po_Param->i_NbImage];
   for (int k=0; k<po_Param->i_NbImage; k++) { 
      po_TabMr2d[k].alloc (po_Param->i_NbLin, po_Param->i_NbCol, po_Param->i_NbPlan, 
                           po_Param->e_Transform, po_Param->po_FAS, po_Param->e_Norm);
      //if (po_Param->e_Transform == TO_LIFTING) 
      //	  po_TabMr2d[k].LiftingTrans = po_Param->e_LiftingTrans;
      //if (po_Param->e_Transform == TO_MALLAT) {
      //   po_TabMr2d[k].SBFilter = po_Param->e_SBFilter;
      //   po_TabMr2d[k].TypeNorm = po_Param->e_Norm;
      //}
   }
   
   // init number of band
   po_Param->i_NbBand = po_TabMr2d[0].nbr_band();
			   
   // class for WTH analysis
   o_Mc2dWth.filt_Alloc (po_Param->i_NbImage);   
   
   // allocate memory for Noise Model classes
   if (po_Param->e_WithNoiseModel) 
      po_TabNoiseModel = new MRNoiseModel [po_Param->i_NbImage];
   else 
      po_TabNoiseModel = (MRNoiseModel*)NULL;
}


wth_Result::~wth_Result () {
  delete [] po_TabMr2d;
}
 
 
void wth_Result::res_FirstCompute () {

  // create all multiresol
  //cout << "Multiresol transform of data ..." << endl;
  for (int k=0; k<po_Param->i_NbImage; k++)
  {
    po_TabMr2d[k].Border = I_MIRROR;
    po_TabMr2d[k].transform (po_Param->po_TabIma[k]);
  }
  if (po_Param->e_WriteTransf) {
     char Mr2dRecons[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(Mr2dRecons, "mr2d_Transf_%d",k+1);
        po_TabMr2d[k].write (Mr2dRecons);
     } 
  }
    
  // create mr2d noise model
  if (po_Param->e_WithNoiseModel) res_InitNoiseModelMr2d();
  
  // calculate the wth transform
  //cout << "WTH analysis ..." << endl;
  o_Mc2dWth.filt_TransfSignal (po_TabMr2d, po_TabMr2d);
 
  // print the result to the standard output
  if (po_Param->e_Verbose) o_Mc2dWth.filt_Print();
  
  // create wth noise model, threshold  
  if (po_Param->e_WithNoiseModel) {
     res_InitNoiseModelWth ();
     o_Mc2dWth.filt_Threshold (po_TabMr2d, po_TabNoiseModel);
  } 
   
  if (po_Param->e_WriteTransf) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "wth_Transf_%d",k+1);
        po_TabMr2d[k].write (NameEigen);
     }
  }
}


void wth_Result::res_CurrentCompute () {

   // create all multiresol
   //cout << "Multiresol transform of data ..." << endl;
   for (int k=0; k<po_Param->i_NbImage; k++) 
      po_TabMr2d[k].transform (po_Param->po_TabIma[k]);
      
   if (po_Param->e_WriteTransf) {
      char Mr2dRecons[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
         sprintf(Mr2dRecons, "mr2d_Transf_%d",k+1);
         po_TabMr2d[k].write (Mr2dRecons);
      } 
   }
   
   // calculate WTH transform 
   //cout << "WTH Analysis..." << endl;
   o_Mc2dWth.filt_TransfSignal (po_TabMr2d, po_TabMr2d);
   
   if (po_Param->e_WithNoiseModel == True) 
      o_Mc2dWth.filt_Threshold (po_TabMr2d, po_TabNoiseModel); 
      	  
   if (po_Param->e_WriteTransf) {
      char NameEigen[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
        sprintf(NameEigen, "wth_Transf_%d",k+1);
        po_TabMr2d[k].write (NameEigen);
      }
   } 
}


void wth_Result::res_Recons () {

   //cout << "Reconstruction process ..." << endl;

   // apply the inverse WTH
   o_Mc2dWth.filt_InvTransform (po_TabMr2d, po_TabMr2d, po_Param->KillLastScale);
   			     		      
   // inv transform Mr2d
   for (int k=0; k<po_Param->i_NbImage; k++)
   {
      if (po_Param->KillLastScale == True) po_TabMr2d[k].band( po_TabMr2d[k].nbr_band()-1).init();
      if (po_Param->e_UseReconsAdjoint) po_TabMr2d[k].rec_adjoint(po_Param->po_TabIma[k]);
      else po_TabMr2d[k].recons(po_Param->po_TabIma[k]);
   }   
   // Write Mr2d recons
   if (po_Param->e_WriteTransf == True) {
      char Mr2dRecons[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
         sprintf(Mr2dRecons, "mr2d_wth_InvTransf_%d",k+1);
         po_TabMr2d[k].write (Mr2dRecons);
         sprintf(Mr2dRecons, "mr2d_wth_InvImag_%d", k+1);
         io_write_ima_float(Mr2dRecons, po_Param->po_TabIma[k]);
      }
   }    
}


void wth_Result::res_Write () {

  //cout << "Write result ..." << endl;
  for (int k=0;k<po_Param->i_NbImage;k++)
    for (int i=0;i<po_Param->i_NbLin;i++)
      for (int j=0;j<po_Param->i_NbCol;j++)
	o_3dDataOut(j,i,k) = po_Param->po_TabIma[k](i,j);

  // write the data
  fits_write_fltarr (po_Param->tc_NameOut, o_3dDataOut);
} 


void wth_Result::res_InitNoiseModelMr2d () {

   // noise model class initialization
   //cout << "Create mr2d noise model ..." << endl;
   for (int k=0; k<po_Param->i_NbImage; k++) {

      po_TabNoiseModel[k].alloc (po_Param->e_TypeNoise, po_Param->i_NbLin, 
                                 po_Param->i_NbCol, po_Param->i_NbPlan,
                                 po_Param->e_Transform);

      if (po_Param->f_NoiseIma > FLOAT_EPSILON) 
         po_TabNoiseModel[k].SigmaNoise = po_Param->f_NoiseIma;
      po_TabNoiseModel[k].CCD_Gain = po_Param->f_Gain;
      po_TabNoiseModel[k].CCD_ReadOutSigma = po_Param->f_SigmaGauss;
      po_TabNoiseModel[k].CCD_ReadOutMean = po_Param->f_MeanGauss;

      if (po_Param->f_NSigmaMr2d >= 0)
         for (int i=0; i<po_Param->i_NbPlan; i++) 
	    po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaMr2d;

      po_TabNoiseModel[k].model (po_Param->po_TabIma[k], po_TabMr2d[k]);
      
      if (po_Param->e_WriteSupportMr2d) {
         char NameSupport[256];
         sprintf(NameSupport, "ww_sup_%d",k+1);
         po_TabNoiseModel[k].write_support_mr (NameSupport);
      }
   }
}


void wth_Result::res_InitNoiseModelWth () {

   // noise model class initialization
   //cout << "Create wth noise model ..." << endl;
   int k;
   o_Mc2dWth.filt_ComputeNoise (po_TabMr2d, po_TabNoiseModel, po_TabNoiseModel);
   
   for (k=0; k<po_Param->i_NbImage; k++) {
   
      if (po_Param->f_NSigmaWth >= 0) {
         for (int i=0; i<po_Param->i_NbPlan; i++)
            po_TabNoiseModel[k].NSigma[i] = po_Param->f_NSigmaWth;
	 po_TabNoiseModel[k].NSigma[0] = 3;
      }
	 
      po_TabNoiseModel[k].set_support (po_TabMr2d[k]);
   }
    
   for (k=0; k<po_Param->i_NbImage; k++)
      if (po_Param->e_WriteSupportWth) {
      char NameSupport[256];
      sprintf(NameSupport, "wth_Support_%d",k+1);
      po_TabNoiseModel[k].write_support_mr (NameSupport);
   }
}


to_Param* wth_Result::res_GetpParam () {return po_Param;}


/******************************************************************************
classe pca_Iter
******************************************************************************/
Bool wth_Iter::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > po_Result->po_Param->i_NbIter) { 
   
      for (int k=0;k<po_Result->po_Param->i_NbImage;k++) 
         for (int i=0;i<po_Result->po_Param->i_NbLin;i++)
            for (int j=0;j<po_Result->po_Param->i_NbCol;j++) {
               po_Result->po_Param->po_TabIma[k](i,j) = o_3dResidu(i,j,k); 
	       //if (   po_Result->po_Param->e_PositivImag 
	       //    && po_Result->po_Param->po_TabIma[k](i,j)<0.) 
               //   po_Result->po_Param->po_TabIma[k](i,j) = 0.0;
	    }
	            
      return False;
   } else return True;
}
 
void wth_Iter::iter_Residu () {

   for (int k=0;k<po_Result->po_Param->i_NbImage;k++)
      for (int i=0;i<po_Result->po_Param->i_NbLin;i++)
         for (int j=0;j<po_Result->po_Param->i_NbCol;j++) {
	 
 	    // positiv constraint
	    if (    po_Result->po_Param->e_PositivImag 
	         && po_Result->po_Param->po_TabIma[k](i,j)<0.) 
               po_Result->po_Param->po_TabIma[k](i,j) = 0.0; 
	          
            // I(i) = I(i) + E(i)
            o_3dResidu(i,j,k) += po_Result->po_Param->po_TabIma[k](i,j);
	                
	    // E(i) = Einit - I(i)
	    po_Result->po_Param->po_TabIma[k](i,j) = o_3dInit(i,j,k) - o_3dResidu(i,j,k);
	    
   }    
}

to_Result* wth_Iter::iter_getpResult() {return po_Result;}







