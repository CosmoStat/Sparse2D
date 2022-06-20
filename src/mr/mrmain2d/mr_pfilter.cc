/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.6
**
**    Author: Philippe Querre
**
**    Date:  00/03/30
**    
**    File:  mr_newpfilter.cc
**
*******************************************************************************/

#include "mr_pfilter.h"
#include "MR_Filter.h"

#define MAX_ITER 10
#define MAX_Nb_EVENTS 100

#define DEFAULT_NUMBER_SCALE 6
#define DEFAULT_NUMBER_MIN_EVENT 4
#define DEFAULT_FIRST_SCALE 0
#define DEFAULT_CONV_PARAM 1e-5
#define DEFAULT_FILTER_TYPE (type_filter)FILTER_ITER_THRESHOLD

#define DEF_TRANS (type_transform)TO_PAVE_BSPLINE
#define DEF_NOISE (type_noise)NOISE_EVENT_POISSON
#define DEF_BORDER_ (type_border)I_CONT



extern char *OptArg;
extern int   OptInd;

extern void mr2d_psupport(MultiResol&     Mr2d_Data, 
                   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border,
		   Bool            WriteAllInfo);
		   
extern void mr2d_psupport(Iint&    EventSignal, 
                   MultiResol&     Mr2d_Data, 
                   Ifloat&         Abaque, 
		   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border, 
		   Bool            WriteAllInfo);		   
		   

int main(int argc, char *argv[]) {


  // licence manager
  //lm2();
  lm_check(LIC_MR1);

  // init local var
  pfilt_Param o_Param;

  // Read param IN
  o_Param (argc, argv);
   
  // init Result and Iter classes
  pfilt_Result o_Result (&o_Param);   

  if (!o_Param.e_OptimSoftThreshold) {
     
     pfilt_Iter o_Iter (&o_Result);
     for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {} 
     
  } else {
     
     pfilt_SoftIter o_Iter (&o_Result);
     for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {}  
         
  }  
  
  if (o_Param.e_TwoStep) {
  
     o_Param.e_CorImag = True;
     if (o_Param.i_NbIter == 0) o_Param.i_NbIter = 10;
        
     if (!o_Param.e_OptimSoftThreshold) {
     
        pfilt_Iter o_Iter (&o_Result);
        for (o_Iter.iter_Begin(); 
             o_Iter.iter_End(); o_Iter++) {} 
     
     } else {
     
        pfilt_SoftIter o_Iter (&o_Result);
        for (o_Iter.iter_Begin(); o_Iter.iter_End(); o_Iter++) {}  
         
     }    
  }

  // write result
  o_Result.res_Write ();
        
  exit(0);
}





/******************************************************************************/
/** classe pca_Param                                                          */
/******************************************************************************/

pfilt_Param::pfilt_Param () {

   //default init of read var
   e_Ascii =               False;
   e_Imag =                False;
   i_NbScale =             DEFAULT_NUMBER_SCALE;     
   i_MinNumberEvent =      DEFAULT_NUMBER_MIN_EVENT;
   i_FirstScale = 	   DEFAULT_FIRST_SCALE;
   i_NbIter =              MAX_ITER;
   f_ConvParam =  	   DEFAULT_CONV_PARAM;
   e_FilterType =          DEFAULT_FILTER_TYPE;
   e_UseAdjoint =          True;
   e_DetectOnlyPos =	   False;
   e_PositivImag =         True;
   e_DilSupport =	   False;
   e_DelIsolPix =	   False;
   e_KillLastScale =	   False;
   e_OptimSoftThreshold =  False;   
   e_CorImag =             False;
   e_TwoStep =             False;
   e_WriteAll = 	   False;
   e_OptZ =		   False;
   e_Border =              DEF_BORDER_;
   i_VMSSize =		   0;
   e_Trace =               False;
   e_Verbose =             False;
   tc_Cmd[0] =                '\0';
}

void pfilt_Param::operator () (int argc, char *argv[]) {

  strcpy(tc_NameAbaque, Name_Abaque_Default);
  for (int k=0; k<argc; k++) sprintf(tc_Cmd, "%s %s", tc_Cmd, argv[k]);
  
  int c,i;
  while ((c = GetOpt(argc,argv,"a:I:f:F:n:e:E:q:i:c:b:TAlkKpwozZ:v")) != -1) {
	
    switch (c) {

    case 'a': /* ascii file */
      if (sscanf(OptArg,"%s", tc_NameIn) != 1) {
         fprintf(stderr, "Error: bad ascii file name: %s\n", OptArg); exit(-1);
      } 
      e_Ascii = True;
      break; 
    case 'I': /* image_filename */
      if (sscanf(OptArg,"%s", tc_NameIn) != 1) {
         fprintf(stderr, "Error: bad image file name: %s\n", OptArg); exit(-1);
      } 
      e_Imag = True;
      break;  		
    case 'f': /* -f <type> type of filtering */
      readint (i, "Error: bad filter type: %i\n");
      testborn (1, 3, i, "Error: bad filter type: %s\n");
      e_FilterType = type_filter(i-1);
      break;    
    case 'F': /* begin the detection */		
      readint (i_FirstScale, "Error: bad first scale detection %d\n");
      testborn (1, 10, i_FirstScale, "Error: bad filter type: %d\n");
      i_FirstScale--;
      break;	
    case 'n': /* -n <NbPlan> */
      readint (i_NbScale, "Error: bad number of scales: %d\n");
      testborn (2, MAX_SCALE, i_NbScale, "bad number of scales: %d\n");
      break;      	
    case 'e': /* -n <min numbre of event> */
      readint (i_MinNumberEvent, "Error: bad min number of event: %d\n");
      testborn (2, MAX_Nb_EVENTS, i_MinNumberEvent, "Error: bad min number of event: %d\n");      
      break;
    case 'E': /* -e < Convergence parameter> */
      readfloat (f_ConvParam, "Error: bad convergence parameter %f\n");
      break;
    case 'q': /* abaque file */
      if (sscanf(OptArg,"%s", tc_NameAbaque) != 1) {
         fprintf(stderr, "Error: bad abaque file name: %s\n", OptArg); exit(-1);
      } 
      break;  	 
    case 'i': /* -i <NbIter> */
      readint (i_NbIter, "bad number of iterarions: %s\n");
      testborn (1, MAX_ITER, i_NbIter, "bad number of iterations: %s\n");
      break;	 
    case 'c': /* corected image file */
      if (sscanf(OptArg,"%s", tc_CorectedImag) != 1) {
         fprintf(stderr, "Error: bad corected imag file name: %s\n", OptArg); exit(-1);
      } 
      e_TwoStep = True;
      break;  
    case 'b':/* -b <border> */
      int i;
      if ((sscanf(OptArg,"%d",&i) != 1)||(i<=0)||(i>3)){
         fprintf(OUTMAN, "bad number for border: %s\n", OptArg);exit(-1);
      }
      switch(i){
         case 1:e_Border=I_CONT;break;
         case 2:e_Border=I_MIRROR;break;
         case 3:e_Border=I_ZERO;break;
      }
    break;
    	            
    case 'o' :
      e_OptimSoftThreshold=True;break;
    case 'T':
      e_Trace = True; break;
    case 'A':
       e_UseAdjoint = False; break;
    case 'p':
      e_DetectOnlyPos = True; break;  
    case 'l':
      e_DilSupport = True; break; 
    case 'k':
      e_DelIsolPix = True; break;
    case 'K':
      e_KillLastScale = True; break;
    case 'w':
      e_WriteAll = True; break; 
#ifdef LARGE_BUFF      
    case 'z':      
      e_OptZ = True;	        
    case 'Z':     
      if (sscanf(OptArg,"%d:%s",&i_VMSSize, tc_VMSName) < 1) {
         fprintf(stderr, "Error: syntaxe is Size:Directory ... \n"); exit(-1);
      }
      e_OptZ = True;
      break;
#endif
    case 'v':
      e_Verbose = True;break;
    case '?': 
      usage(argv); break;
    }
  }

  // get optional input file names from trailing parameters and open files 
  // read output image names */
  if (OptInd < argc) strcpy(tc_NameOut, argv[OptInd++]);
  else usage(argv);
    
  // make sure there are not too many parameters  
  if (OptInd < argc) {
    fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  }
  
  if (e_Ascii == False && e_Imag==False) {
      fprintf(OUTMAN, "\n\nError: option -a or -i has to be set\n");exit(-1);
  }
  
  read_image_and_abaque ();
  i_NbCol = o_Imag.nc();
  i_NbLin = o_Imag.nl();
//cout << "Line:" << i_NbLin << ",Col:" << i_NbCol << endl;
  check_scale(i_NbLin, i_NbCol, i_NbScale);
  
#ifdef LARGE_BUFF
  if (e_OptZ == True) vms_init(i_VMSSize, tc_VMSName, e_Verbose);
#endif
  
  if (e_Trace) trace ();
}


void pfilt_Param::read_image_and_abaque () {

   if (e_Verbose) {
      if (!e_TwoStep)
         cout << "Imag file IN         : " << tc_NameIn << endl;
      else {
         cout << "compute support on   : " << tc_NameIn << endl;
	 cout << "Imag file IN         : " << tc_CorectedImag << endl;
      }
      cout << "Abaque file IN       : " << tc_NameAbaque << endl;
   }
     
   // read image
   if (e_Ascii) {
      building_imag_ascii(tc_NameIn, o_EventImage, o_Imag);
      //building_imag_ascii(tc_CorectedImag, o_EventCorImag, o_CorImag);
   } else if (e_Imag) {
      io_read_ima_float(tc_NameIn, o_Imag, &s_Header);
 
      building_imag_imag(o_Imag, o_EventImage);       
      if (e_Trace) 
         io_write_ima_float ("imag_inv0.fits", o_Imag, &s_Header); 
      
      //io_write_ima_float("test_proj.fits", o_Imag);
      //io_read_ima_float(tc_NameIn, o_CorImag, &s_Header); 
      //building_imag_imag(o_Imag, o_EventCorImag); 
//io_write_ima_float("totest.fits", o_Imag);         
   } 
   
   if (e_TwoStep) {
      io_read_ima_float(tc_CorectedImag, o_CorImag, &s_Header);
      if (  o_CorImag.nx() != o_Imag.nx() 
         || o_CorImag.ny() != o_Imag.ny()) {
	 cout << "Images must have same size" << endl; exit(-1);
      }
      Iint o_Null;
      building_imag_imag(o_CorImag, o_Null);  
      if (e_Trace) 
         io_write_ima_float ("corected_inv0.fits", o_CorImag, &s_Header);
   }
   
   // update Header
   //--------------
   s_Header.origin = tc_Cmd;  
   
   // read abaque
   io_read_ima_float(tc_NameAbaque, o_Abaque);
   io_set_format(F_UNKNOWN);   
   
   if (e_Trace) {
      io_write_ima_float ("Imag_In.fits", o_Imag, &s_Header);
      io_write_ima_float ("Abaque_In.fits", o_Abaque);
      io_write_ima_int ("Event_In.fits", o_EventImage); 
   } 
}




void pfilt_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options out_signal\n", argv[0]);
  manline(); 
  poisson_usage();                    manline();
  corected_imag();                    manline();
  filterp_usage(e_FilterType);        manline();  
  first_detect_scale_usage();         manline();  
  nbr_scale_usage(i_NbScale);         manline();  
  min_event_usage(i_MinNumberEvent);  manline();  
  convergp_param_usage(f_ConvParam);  manline();
  abaque_file_usage();                manline(); 
  border_usage(e_Border);             manline(); 
  max_iter_usage(MAX_ITER);           manline();  
  // use_adjoint();                   manline();
  detect_pos_usage();                 manline(); 
  dilate_sup_usage();                 manline(); 
  kill_isol_pix_usage();              manline();   
  kill_last_scale_usage();            manline(); 
  vm_usage();                         manline();      
  verbose ();                         manline();
  exit(-1);
}



void pfilt_Param::trace () {

   cout << "Param IN :" << endl;
   if (e_Ascii) {
      cout << " -- [-a:] : ascii in : " << tc_NameIn << endl;
   }
   if (e_Imag)  {
      cout << " -- [-I:] : imag in : " << tc_NameIn << endl;
   }
   if (e_CorImag) 
      cout << " -- [-c:] : corrected image : " << tc_CorectedImag << endl;
   cout << " -- [-f]  : filter type : " << StringFilter(e_FilterType) << endl;
   cout << " -- [-F]  : first scale : " << i_FirstScale << endl;
   cout << " -- [-n:] : number of plan : " << i_NbScale << endl;
   cout << " -- [-e:] : min number of event : " << i_MinNumberEvent << endl;
   cout << " -- [-E:] : convergence parameter : " << f_ConvParam << endl;  
   cout << " -- [-q:] : abaque file name : " <<  tc_NameAbaque << endl;
   cout << " -- [-i:] : number of iteration : " << i_NbIter << endl;
   cout << " -- [-b:] : border              : " << e_Border << endl;
   if (e_OptimSoftThreshold) cout << " -- [-o]  : optim : soft threshold..." << endl;
   if (e_UseAdjoint)
      cout << " --       : reconstruction with adjoint operator" << endl;
   else
      cout << " -- [-A]  : reconstruction without adjoint operator" << endl; 
   if (e_DetectOnlyPos) 
      cout << " -- [-p]  : detect only positiv structure" << endl;  
   if (e_DilSupport)
      cout << " -- [-l]  : dilate support" << endl; 
   if (e_DelIsolPix)
      cout << " -- [-k]  : remove isolated pixels" << endl;   
   if (e_KillLastScale)
      cout << " -- [-K]  : kill last scale" << endl;         
   if (e_WriteAll)
      cout << " -- [-w]  : write MR transform and Threshold levels" << endl;            
   if (e_Verbose) 
      cout << " -- [-v]  : verbose" << endl;
}




/******************************************************************************
classe pca_Result
******************************************************************************/

pfilt_Result::pfilt_Result (pfilt_Param* ppo_Param) {

   po_Param = ppo_Param; 
   LoopNumber = 0;
   
   // allocate fltarray for write result
   o_DataOut.init(po_Param->o_Imag);
    
   // allocate the memory for all multiresol
   po_Mr2d = new MultiResol ();
   po_Mr2d->alloc (po_Param->i_NbLin, po_Param->i_NbCol, 
                   po_Param->i_NbScale, DEF_TRANS, "pfilter");
   po_Mr2d->TypeNorm = NORM_L1;
   po_Mr2d->SBFilter = F_MALLAT_7_9;
		  
   //init noise model
   o_NoiseModel.NoCompDistrib = True;	  
   o_NoiseModel.alloc (DEF_NOISE, po_Param->i_NbLin, 
                       po_Param->i_NbCol, po_Param->i_NbScale,
                       DEF_TRANS);
			
   o_NoiseModel.MinEventNumber = po_Param->i_MinNumberEvent;
   o_NoiseModel.MinEvent = True;
   o_NoiseModel.TransImag = False;
   o_NoiseModel.SigmaNoise = 1.;
   o_NoiseModel.OnlyPositivDetect = po_Param->e_DetectOnlyPos;
   o_NoiseModel.FirstDectectScale = po_Param->i_FirstScale;
		  		   
   //init sigma noise
   SigmaNoise = o_NoiseModel.SigmaNoise;
   CurSigma   = energy(po_Param->o_Imag);
   OldSigma   = CurSigma;
}


void pfilt_Result::res_FirstCompute () {  

  // select image to process...
  if (po_Param->e_CorImag)  
     po_Param->o_Imag = po_Param->o_CorImag;
       
  //init attribut.....
  LoopNumber = 0;
  SigmaNoise = o_NoiseModel.SigmaNoise;
  CurSigma   = energy(po_Param->o_Imag);
  OldSigma   = CurSigma;
  
  // call number
  LoopNumber = ++LoopNumber;   
  
  if (!po_Param->e_CorImag) {
  
     // boder type...
     type_border NM_border = po_Param->e_Border;
  
     // create multiresol
     if (po_Param->e_Trace)
        cout << "Multiresol transform of data (mr_support)..." << endl;
     po_Mr2d->transform (po_Param->o_Imag, NM_border);
    
     if (po_Param->e_Trace) {
        char Mr2dRecons[256];
        sprintf(Mr2dRecons, "mr2d_FirstTransf");
        po_Mr2d->write (Mr2dRecons);
     } 

     //compute mr_support
     mr2d_psupport(po_Param->o_EventImage, *po_Mr2d, 
                   po_Param->o_Abaque, o_NoiseModel, 
                   NM_border, po_Param->e_WriteAll);

     if (po_Param->e_WriteAll) 
        o_NoiseModel.write_support_ima("xx_Support.fits"); 
  } 
     
  // set OnlyPositivDetect to true (case poisson noise...)
  //po_Param->e_DetectOnlyPos = True;
  //o_NoiseModel.OnlyPositivDetect = True;
  
  // create all multiresol    
  if (po_Param->e_Trace)
     cout << "Multiresol transform of data ..." << endl;
  po_Mr2d->transform (po_Param->o_Imag, po_Param->e_Border);
     
   
  if (po_Param->e_Trace) {
      char Mr2dRecons[256];
      sprintf(Mr2dRecons, "mr2d_Transf_%d", LoopNumber);
      po_Mr2d->write (Mr2dRecons);
  } 
  
  // theshold 
  if (po_Param->e_OptimSoftThreshold) {  
     res_SoftInitAlloc ();
     res_SoftThreshold();
  } else {
     if (po_Param->e_CorImag)
        o_NoiseModel.threshold (*po_Mr2d, False);
     else 
	o_NoiseModel.threshold (*po_Mr2d, True);
  }
  
  if (po_Param->e_KillLastScale)
     po_Mr2d->scale(po_Mr2d->nbr_scale()-1).init();  
     
  if (po_Param->e_Trace) {
      char Mr2dRecons[256];
      sprintf(Mr2dRecons, "mr2d_ThreshTransf_%d", LoopNumber);
      po_Mr2d->write (Mr2dRecons);
  } 
}


void pfilt_Result::res_CurrentCompute () {   
  
   // call number
   LoopNumber = ++LoopNumber;
    
   // create multiresol
   if (po_Param->e_Trace)
      cout << "Multiresol transform of data ..." << endl;
   po_Mr2d->transform (po_Param->o_Imag, po_Param->e_Border);
      
   if (po_Param->e_Trace) {
      char Mr2dRecons[256];
      sprintf(Mr2dRecons, "mr2d_Transf_%d", LoopNumber);
      po_Mr2d->write (Mr2dRecons);
   }
   
   // theshold   
   if (po_Param->e_OptimSoftThreshold) { 
      res_SoftThreshold();
   } else {
      if (po_Param->e_CorImag) 
         o_NoiseModel.threshold (*po_Mr2d, False);
      else 
         o_NoiseModel.threshold (*po_Mr2d, True);
   }
   
   if (po_Param->e_KillLastScale)
      po_Mr2d->scale(po_Mr2d->nbr_scale()-1).init();   
      
   if (po_Param->e_Trace) {
      char Mr2dRecons[256];
      sprintf(Mr2dRecons, "mr2d_ThreshTransf_%d", LoopNumber);
      po_Mr2d->write (Mr2dRecons);
  }   
}


void pfilt_Result::res_Recons () {

   if (po_Param->e_Trace)
      cout << "Reconstruction process ..." << endl;    
  
   // inv transform Mr2d
   if (!po_Param->e_UseAdjoint) po_Mr2d->recons(po_Param->o_Imag);
   else po_Mr2d->rec_adjoint(po_Param->o_Imag);
            
   // Write Mr2d recons
   if (po_Param->e_Trace == True) {
      char Mr2dRecons[256];
      sprintf(Mr2dRecons, "RecImage_%d", LoopNumber);
      io_write_ima_float (Mr2dRecons, po_Param->o_Imag);
   } 
}


void pfilt_Result::res_Write () {

  if (po_Param->e_Trace) cout << "Write result ..." << endl;

  o_DataOut = po_Param->o_Imag;

  // write the data  
  po_Param->s_Header.bitpix = BP_FLOAT;
  io_write_ima_float (po_Param->tc_NameOut, o_DataOut, &po_Param->s_Header);
} 


pfilt_Param* pfilt_Result::res_GetpParam () {return po_Param;}



/******************************************************************************
classe SoftIter
******************************************************************************/

void pfilt_Result::res_SoftInitAlloc () {
   
   if (po_Param->e_Trace)
      cout << "Soft threshold alloc ..." << endl;
      
   FirstCoef.alloc (po_Mr2d->nbr_band(),
                    po_Param->i_NbLin,po_Param->i_NbCol);
   ResidualCoef.alloc (po_Mr2d->nbr_band(),
                       po_Param->i_NbLin,po_Param->i_NbCol);	   
   for (int s=0; s<po_Mr2d->nbr_band(); s++) {
      for (int i=0; i<po_Param->i_NbLin; i++) {
         for (int j=0; j<po_Param->i_NbCol; j++) {
      
            FirstCoef(s,i,j) = (*po_Mr2d)(s,i,j);
	    (*po_Mr2d)(s,i,j) = 0.;
         }
      }
   } 
}

void pfilt_Result::res_SoftThreshold () {

   if (po_Param->e_Trace)
      cout << "Soft threshold process ..." << endl;
         
   // compute coef for lamba at current iter 
   // Lambda = Lambda * (1-currrent iter)/NbIter
   //-------------------------------------------
   float CoefLambda=1.;
   if (po_Param->i_NbIter != 1)
      CoefLambda =   (1. - (LoopNumber-1.)/(po_Param->i_NbIter-1));
   if (po_Param->e_Trace) cout << "Coef Lambda : " << CoefLambda << endl;
   
   // threshold
   //----------
   for (int s=0; s<po_Mr2d->nbr_band()-1; s++) {
	    
      for (int i=0; i<po_Param->i_NbLin; i++) {
         for (int j=0; j<po_Param->i_NbCol; j++) {
	         
            float LocalSigma = o_NoiseModel.sigma(s,i,j);
	       
	    // compute residual on PCA coef
            //----------------------------- 	    
	    ResidualCoef(s,i,j) =           
	       FirstCoef(s,i,j) - (*po_Mr2d)(s,i,j);
	       
	    // first rule, if FirstCoef is significatif and
            // abs(ResidualCoef) > eps sigma then coef = FirstCoef
            //----------------------------------------------------------	 
	    if (o_NoiseModel(s,i,j)) {
	              
	       // Pix is signif, if abs(res) > sig/2 => sol = init      
	       if ((ABS(ResidualCoef(s,i,j)) > LocalSigma/2.)) { 
	           (*po_Mr2d)(s,i,j) = FirstCoef(s,i,j);
               } 
	    }
	    
	    //soft threshold
	    //--------------	    
	    float SoftThres = CoefLambda * o_NoiseModel.sigma(s,i,j); 
	    (*po_Mr2d)(s,i,j)  = soft_threshold ((*po_Mr2d)(s,i,j),SoftThres);
	 }
      }
   }
   for (int i=0; i<po_Param->i_NbLin; i++) {
      for (int j=0; j<po_Param->i_NbCol; j++) {
         (*po_Mr2d)(po_Mr2d->nbr_band()-1,i,j) = FirstCoef(po_Mr2d->nbr_band()-1,i,j);
      }
   }
}
   
/******************************************************************************
classe pfilt_Iter
******************************************************************************/
void pfilt_Iter::iter_Begin () {
  
  if (!po_Result->po_Param->e_CorImag)         
     o_Init = po_Result->po_Param->o_Imag; 
  else                                        
     o_Init = po_Result->po_Param->o_CorImag; 

  o_FilteredImag.init(iter_getpResult()->res_GetpParam()->o_Imag);
  //o_Init.init(iter_getpResult()->res_GetpParam()->o_Imag); 
  
  NbIter = 0;

  //o_Init = iter_getpResult()->res_GetpParam()->o_Imag;
  o_FilteredImag.init(0.0);  
                                         // I(0) = 0
//char File[80];
//sprintf(File, "ImagIn_%d.fits", NbIter+1);
//io_write_ima_float (File, iter_getpResult()->res_GetpParam()->o_Imag);         
  
  NbIter++;
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_Residu();

}

Bool pfilt_Iter::iter_End () {

//   cout << " ==> END" << endl;
   
   po_Result->OldSigma = po_Result->CurSigma;
   po_Result->CurSigma = energy (po_Result->po_Param->o_Imag);
   float Nl = po_Result->po_Param->o_Imag.nl();
   float Nc = po_Result->po_Param->o_Imag.nc();
   float Delta =   (po_Result->OldSigma - po_Result->CurSigma) 
                 / (Nl*Nc*po_Result->SigmaNoise*po_Result->SigmaNoise);
   float Xi    =   po_Result->CurSigma
                 / (Nl*Nc*po_Result->SigmaNoise*po_Result->SigmaNoise);
   cout << "Iter:" << NbIter << ", Xi=" << Xi << " ==> Delta=" << 
                                           Delta << endl;
   
   NbIter++;
   
//cout << "end loop: " << po_Result->SigmaNoise << ", " << po_Result->CurSigma << ", " <<  po_Result->OldSigma << endl; 
//cout << "          nrj(In) :" << energy(po_Result->po_Param->o_Imag) << endl;   
//cout << "          nrj(Out):" << energy(o_FilteredImag) << endl;  
      
   if (    NbIter >=po_Result->po_Param->i_NbIter+1
        || (   !po_Result->po_Param->e_CorImag 
	     && Delta  < po_Result->po_Param->f_ConvParam)) { 
   
      po_Result->po_Param->o_Imag = o_FilteredImag;
      return False;  
   } else {
      return True;
   }
}

void pfilt_Iter::iter_Residu () {
		  		          
   // add new (current) filtered image to old (previous iter) filtered image 
   if (po_Result->po_Param->e_Trace) {
      char IterName[256];
      sprintf(IterName, "FiltRes_%d.fits", NbIter);
      io_write_ima_float (IterName, po_Result->po_Param->o_Imag);
   }   
   o_FilteredImag += po_Result->po_Param->o_Imag;
   if (po_Result->po_Param->e_Trace) {
      char IterName[256];
      sprintf(IterName, "FiltRecons_%d.fits", NbIter);
      io_write_ima_float (IterName, o_FilteredImag);
   }
   
   // positiv constraint on filtered image
   if (po_Result->po_Param->e_PositivImag) {
      threshold(o_FilteredImag);
      if (po_Result->po_Param->e_Trace) {
         char IterName[256];
         sprintf(IterName, "FiltReconsThresh_%d.fits", NbIter);
         io_write_ima_float (IterName, o_FilteredImag);
      }      
   }
         
   // compute residu (init - filtered image)
   po_Result->po_Param->o_Imag = o_Init - o_FilteredImag;
   if (po_Result->po_Param->e_Trace) {
      char IterName[256];
      sprintf(IterName, "Residu_%d.fits", NbIter);
      io_write_ima_float (IterName, po_Result->po_Param->o_Imag);
   }   
}

void pfilt_Iter::operator++ () {(*this)++;}
void pfilt_Iter::operator++ (int) {
      
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_Residu ();
}

pfilt_Result* pfilt_Iter::iter_getpResult() {return po_Result;}

/******************************************************************************
classe to_SoftIter
******************************************************************************/
void pfilt_SoftIter::iter_Begin () {
  //cout << " ==> BEGIN" << endl;
  NbIter = 0;                                        
   
  NbIter++;
  
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_ActionOnSol();
}

Bool pfilt_SoftIter::iter_End () {   NbIter++;
      
   if (NbIter > po_Result->po_Param->i_NbIter) {             
      return False;
   } else return True;
}

void pfilt_SoftIter::operator++ () {(*this)++;}
void pfilt_SoftIter::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_ActionOnSol ();
}
  
void pfilt_SoftIter::iter_ActionOnSol() {   

   // positiv if  necessary...
   if (po_Result->po_Param->e_DetectOnlyPos) 
      threshold(po_Result->po_Param->o_Imag);   
      
   cout << "energy=" << po_Result->po_Param->o_Imag.energy() << endl;
      
   // write  sol at current iteration	
   if (po_Result->po_Param->e_WriteAll== True) {
       char Mr1dRecons[256];
       sprintf(Mr1dRecons, "current sol_loop_%d",NbIter);
       io_write_ima_float (Mr1dRecons, po_Result->po_Param->o_Imag);
   }
}    

pfilt_Result* pfilt_SoftIter::iter_getpResult() {return po_Result;}
























