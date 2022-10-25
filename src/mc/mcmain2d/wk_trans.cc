/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.8
**
**    Author: 01/23/01
**
**    Date:  01/02/09
**    
**    File:  mc_pca.cc
**
*******************************************************************************/

#include "wk_trans.h"

#define MAX_ITER 100
#define DEF_SIGMA_PCA -1
#define DEF_SIGMA_MR2D 3

#define DEF_CORREL_COMP (CorrelComputeType) E_LOCAL
#define DEF_CORREL_NOISE (CorrelNoiseType) E_WITHOUT

extern char *OptArg;
extern int   OptInd; 

int main(int argc, char *argv[]) 
{
   extern softinfo Soft;
   Soft.mr3();
   lm_check(LIC_MR3);
   
  // licence manager
  //lm2();

  // init local var
  pca_Param o_Param;

  // Read param IN
  o_Param (argc, argv);
      
  // init Result and Iter classes
  pca_Result o_Result (&o_Param);   

  // compute
  o_Result.res_FirstCompute ();

  // write result
  o_Result.res_Write ();
        
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
  e_Bord =                I_CONT;
  po_FAS =                (FilterAnaSynt*)NULL;
  i_NbrUndec =            -1;
  i_NbPlan =              DEFAULT_NBR_SCALE;
  e_NormCorrelMatrix =    True;
  e_Verbose =             False;
  e_NoDisplay =           True;
  e_WriteEigenOnDisk =    True;
  e_WriteTransf =         False;
  e_TraceParamIn =        False;
  e_WriteCorrelMatrix =   False;  
  e_CorrelComp =          (CorrelComputeType)DEF_CORREL_COMP;
  po_ImportCorMat =       NULL;
}

void pca_Param::operator () (int argc, char *argv[]) {

  int c,i;  
  Bool readfile = False, Optf = False;
  while ((c = GetOpt(argc,argv,"t:n:WkvpLCx:O:u:T:")) != -1) {
	
    switch (c) {

    case 't': /* -d <type> type of transform */
      readint (i, "bad type of multiresolution transform: %s\n");
      if (testborn (1 ,NBR_TRANSFORM, i, "bad type of transform: %s\n"))
	e_Transform = (type_transform) (i-1);
      break;
           
    case 'n': /* -n <NbPlan> */
      readint (i_NbPlan, "bad number of scales: %s\n");
      testborn (2, MAX_SCALE, i_NbPlan, "bad number of scales: %s\n");
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
      
   case 'u': /* -u number of undecimated scale */
      readint (i_NbrUndec, "Error: bad number of scales: %s\n");
      if (testborn (0 ,MAX_SCALE-1, i_NbrUndec, "bad number of scales: %s\n"))	
      break;
      
   case 'T': /* -T get filter... */
      Optf = True;
      e_SBFilter = get_filter_bank(OptArg);
      break;     
            
    case 'L':
      e_Norm = NORM_L2;break;     
    case 'p':
      e_NormCorrelMatrix = False;break;              
    case 'C':
      e_WriteCorrelMatrix = True;break;
    case 'w':
      e_WriteEigenOnDisk = True;break;
    case 'W':
      e_WriteTransf = True;break;
    case 'k':
      e_TraceParamIn = True;break;
    case 'v':
      e_Verbose = True;
      e_NoDisplay = False;
      break;
    case '?':
      usage(argv);
    }
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
  
  if (   (e_Transform != TO_UNDECIMATED_MALLAT)
      && (e_Transform != TO_MALLAT)
      && (Optf == True) ) {
       cout << "Error: option -T is only valid with Mallat transform ... \n" << endl;
       exit (-1);
    } 
    
  if (   (e_Transform == TO_UNDECIMATED_MALLAT)
      || (e_Transform == TO_MALLAT) ) {
      
    po_FAS = new FilterAnaSynt(); /* mem loss, never free ... */
    po_FAS->alloc(e_SBFilter);
    po_FAS->Verbose = e_Verbose;
  }
  
  if (e_TraceParamIn) trace ();
}

void pca_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options input_cube output_prefix\n", argv[0]);
  manline(); 

  transform_usage (e_Transform);      manline();
  wk_filter (F_MALLAT_7_9);           manline();
  //nbr_nbr_undec_usage(i_NbrUndec);    manline();
  nbr_scale_usage (i_NbPlan);         manline();
  // write_eigen_on_disk ();             manline();
  write_correl_matrix ();             manline();
  // normalize_correl_matrix ();         manline();
  verbose ();                         manline();  
  correl_compute (DEF_CORREL_COMP);   manline();
  import_correl_file();               manline();
  exit(-1);
}

void pca_Param::read_image (char *argv[]) {

  if (e_Verbose) cout << "Name file IN : " << tc_NameIn << endl;
  if (e_Verbose) cout << "Name file OUT : " << tc_NameOut << endl;

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
}


void pca_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransform(e_Transform) << endl;
   if ((e_Transform == TO_MALLAT) || (e_Transform == TO_UNDECIMATED_MALLAT)) {
      cout << " -- [-L]  : Normalisation L2" << endl;
      cout << " -- [-f]  : " << StringSBFilter(e_SBFilter) << endl;
      cout << " -- [-u:] : number of undecimated scales : " << i_NbrUndec << endl;
   }  
   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect images" << endl;
   if (e_WriteEigenOnDisk) cout << " -- [-w]  : write eigen vector images" << endl;  
   if (e_WriteCorrelMatrix) cout << " -- [-C]  : write correlation matrix" << endl;
   cout << " -- [-x:] : correlation compute : " << CorCompTransform(e_CorrelComp) << endl;
   if (e_CorrelComp == E_IMPORT)
      cout << " -- [-O:] : import correlation file : " << tc_ImportCorMatName << endl;
   if (e_NormCorrelMatrix) cout << " -- [-p]  : normalize correl matrix" << endl;
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
      po_TabMr2d[k].alloc (po_Param->i_NbLin, 
                           po_Param->i_NbCol, 
			   po_Param->i_NbPlan, 
                           po_Param->e_Transform, 
			   po_Param->po_FAS, 
			   po_Param->e_Norm,
			   po_Param->i_NbrUndec);
   }
			   			   
   // init number of band
   po_Param->i_NbBand = po_TabMr2d[0].nbr_band();
   
   // class for image correlation analysis
   o_Mr2dPca.alloc (po_Param->i_NbImage, po_Param->i_NbBand, 
                    po_Param->e_WriteCorrelMatrix);
   for (int b=0; b<po_Param->i_NbBand; b++) {
      o_Mr2dPca.MatCor[b].Verbose = po_Param->e_Verbose;
   }
}


pca_Result::~pca_Result () {
  delete [] po_TabMr2d;
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
  
   // calculate the correlation matrix and its eigen values
   if (po_Param->e_TraceParamIn)
      cout << "Principal component analysis ..." << endl;
   o_Mr2dPca.compute (po_TabMr2d, (MRNoiseModel*)NULL, 
                      po_Param->e_NormCorrelMatrix,
                      po_Param->e_CorrelComp,
                      DEF_CORREL_NOISE, 
		      po_Param->po_ImportCorMat);
 
   // print the result to the standard output
   if (!po_Param->e_NoDisplay) o_Mr2dPca.print();
   
   // write eigen vectors on disk
   if (po_Param->e_WriteEigenOnDisk) {
  
      // calculate the eigen vector images 
      if (po_Param->e_TraceParamIn)
         cout << "Compute eigen vector images ..." << endl;
      o_Mr2dPca.pca_TransfSignal (po_TabMr2d, po_TabMr2d);
   
      if (po_Param->e_WriteTransf) {
         char Mr2dRecons[256];
         for (int k=0; k < po_Param->i_NbImage; k++) {
            sprintf(Mr2dRecons, "mr2d_PcaTransf_%d",k+1);
            po_TabMr2d[k].write (Mr2dRecons);
         } 
      }  
   
      char NameEigen[256];
      for (int k=0; k < po_Param->i_NbImage; k++) {
         sprintf(NameEigen, "%s_ev_%d",po_Param->tc_NameOut,k+1);
         po_TabMr2d[k].write (NameEigen);
      }
   }
}


to_Param* pca_Result::res_GetpParam () {return po_Param;}

