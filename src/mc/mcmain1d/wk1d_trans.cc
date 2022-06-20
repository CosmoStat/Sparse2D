/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.5
**
**    Author: 02/29/00
**
**    Date:  01/02/09
**    
**    File:  mc1d_pcafilter.cc
**
*******************************************************************************
**
******************************************************************************/

#include "wk1d_trans.h"


#define MAX_ITER 100
#define DEF_SIGMA_PCA -1
#define DEF_SIGMA_MR1D 3
#define DEFAULT_NBR_SCALE 4
#define DEF_TRANSFORM TO1_PAVE_B3SPLINE
#define MAX_SCALE 10
#define DEF_CORREL_COMP (CorrelComputeType1d) E_LOCAL1D
#define DEF_CORREL_NOISE (CorrelNoiseType1d) E_WITHOUT1D


extern int  OptInd;
extern char *OptArg;
//extern int  getopt(int argc, char *const*argv, char *opts); 


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
    
   // init Result and Iter classes
  pca_Result o_Result (&o_Param);   

  // compute
  o_Result.res_FirstCompute ();

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
  e_CorrelComp =          (CorrelComputeType1d)DEF_CORREL_COMP;
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
	e_Transform = (type_trans_1d) (i-1);
      break;
           
    case 'n': /* -n <NbPlan> */
      readint (i_NbPlan, "bad number of scales: %s\n");
      testborn (2, MAX_SCALE, i_NbPlan, "bad number of scales: %s\n");
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
      e_NoDisplay = False;
      e_Verbose = True;break;
    case '?':
      usage(argv);
    }
  }  

  // get optional input file names from trailing parameters and open files 
  // read input spetral names */
  if (OptInd < argc) strcpy(tc_NameIn, argv[OptInd++]);
  else usage(argv);
  if (OptInd < argc) strcpy(tc_NameOut, argv[OptInd++]);
  else usage(argv);
  
  // make sure there are not too many parameters  
  if (OptInd < argc) {
    fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  }

  read_spectre (argv);  // matrix file name
  
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

  fprintf(OUTMAN, "Usage: %s options in_spectra output_prefix\n", argv[0]);
  manline(); 

  //transform_usage (e_Transform);      manline();
  //wk_noise_usage ();                     manline();
  nbr_scale_usage (i_NbPlan);         manline();
  // write_eigen_on_disk ();             manline();
  write_correl_matrix ();             manline();
  //normalize_correl_matrix ();         manline();
  correl_compute1d (DEF_CORREL_COMP);   manline();
  import_correl_file();               manline(); 
  verbose ();                         manline();
  exit(-1);
}

void pca_Param::read_spectre (char *argv[]) {

  if (e_Verbose) cout << "Name file IN : " << tc_NameIn << endl;
  if (e_Verbose) cout << "Name file OUT : " << tc_NameOut << endl;
 
  // read the data (-> fltarray)
  fits_read_fltarr (tc_NameIn, ao_2dDataIn);
  i_NbPts     = ao_2dDataIn.nx(); //point
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
}


void pca_Param::trace () {
   
   cout << "Param IN :" << endl;
   cout << " -- [-t:] : type of transform : " << StringTransf1D(e_Transform) << endl;
   if ((e_Transform == TO1_MALLAT)) {
      cout << " -- [-L]  : Normalisation L2" << endl;
      cout << " -- [-f]  : " << StringSBFilter(e_SBFilter) << endl;
      cout << " -- [-u:] : number of undecimated scales : " << i_NbrUndec << endl;
   }  
   cout << " -- [-n:] : number of plan : " << i_NbPlan << endl;
   if (e_WriteTransf) cout << " -- [-W]  : write transf vect images" << endl;
   if (e_WriteEigenOnDisk) cout << " -- [-w]  : write eigen vector images" << endl;  
   if (e_WriteCorrelMatrix) cout << " -- [-C]  : write correlation matrix" << endl;
   cout << " -- [-x:] : correlation compute : " << CorCompTransform1d(e_CorrelComp) << endl;
   if (e_CorrelComp == E_IMPORT1D)
      cout << " -- [-O:] : import correlation file : " << tc_ImportCorMatName << endl;
   if (e_NormCorrelMatrix) cout << " -- [-p]  : normalize correl matrix" << endl;
   if (e_Verbose) cout << " -- [-v]  : verbose" << endl;
}



/******************************************************************************
classe to_Result
******************************************************************************/

pca_Result::pca_Result (pca_Param* ppo_Param) {

   po_Param = ppo_Param;
   // allocate fltarray for write result
   o_2dDataOut.alloc (po_Param->i_NbPts, po_Param->i_NbSpectre);
   // allocate the memory for all multiresol
   po_TabMr1d = new MR_1D [po_Param->i_NbSpectre];
   for (int k=0; k<po_Param->i_NbSpectre; k++) {
      po_TabMr1d[k].alloc(po_Param->i_NbPts, po_Param->e_Transform, po_Param->i_NbPlan);
   }
   po_Param->i_NbBand = po_TabMr1d[0].nbr_band();			   
		   
   // class for spectral correlation analysis
   o_Mr1dPca.alloc (po_Param->i_NbSpectre, po_Param->i_NbPlan,
                    po_Param->e_WriteCorrelMatrix);
    
   for (int b=0; b<po_Param->i_NbBand; b++) {
      o_Mr1dPca.MatCor[b].Verbose = po_Param->e_Verbose;
   }  		   
}


pca_Result::~pca_Result () {
  delete [] po_TabMr1d;
}


void pca_Result::res_FirstCompute () {

  // create all multiresol
  if (po_Param->e_TraceParamIn)
     cout << "Multiresol transform of data ..." << endl;
  for (int k=0; k<po_Param->i_NbSpectre; k++)
    po_TabMr1d[k].transform (po_Param->po_TabSpectre[k]);
    
  if (po_Param->e_WriteTransf) {
     char Mr1dRecons[256];
     for (int k=0; k <po_Param->i_NbSpectre; k++) {
        sprintf(Mr1dRecons, "mr1d_Transf_%d",k+1);
        fits_write_fltarr (Mr1dRecons, po_TabMr1d[k].image());
     } 
  }
  
  // calculate the correlation matrix and its eigen values
  if (po_Param->e_TraceParamIn)
      cout << "Principal component analysis ..." << endl;  
  o_Mr1dPca.compute (po_TabMr1d, (MR1DNoiseModel*)NULL,
                      po_Param->e_NormCorrelMatrix,
                      po_Param->e_CorrelComp,
                      DEF_CORREL_NOISE, 
		      po_Param->po_ImportCorMat);
 
  // print the result to the standard output
  if (!po_Param->e_NoDisplay) o_Mr1dPca.print();
  
  // calculate the eigen vector spectres 
  if (po_Param->e_TraceParamIn)
     cout << "Compute eigen vector ..." << endl; 
  o_Mr1dPca.pca_TransfSignal (po_TabMr1d, po_TabMr1d);
  
  if (po_Param->e_WriteEigenOnDisk) {
     char NameEigen[256];
     for (int k=0; k < po_Param->i_NbSpectre; k++) {
        sprintf(NameEigen, "%s_ev_%d",po_Param->tc_NameOut,k+1);
        fits_write_fltarr (NameEigen, po_TabMr1d[k].image());
     }
  }
}

to_Param1d* pca_Result::res_GetpParam () {return po_Param;}


