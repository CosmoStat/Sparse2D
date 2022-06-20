/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.3
**
**    Author: Phlippe Querre
**
**    Date:  00/08/11
**    
**    File:  mr_abaque.cc
**
******************************************************************************/

#include "IM_Obj.h"
#include "Mr_FewEvent.h"

class nfe_Param    {
public:   
   char tc_AbaqueName[256];      /* abaque name */
   float f_Epsilon;              /* precision number */
   int i_NbAutoConv;             /* number max of event (of autoconvolution) */
   Bool e_Trace;                 /* trace param */
   Bool e_WriteAllFiles;         /* write all inter files */
   Bool e_AllDefTrue;            /* put all def at true */
   Bool e_Verbose;               /* verbose exec */
public:
   nfe_Param ();
   virtual void operator () (int argc, char *argv[]); /* read parameters */
private:
   virtual void usage (char *argv[]);      /* trace usage */
   virtual void trace ();                  /* trace parameters */    
};



#define DEF_EPSILON  1e-3
#define DEF_ABAQUENAME "Abaque.fits"
#define DEF_NBR_AUTOCONV 25

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

#define MaxVal_in_Tab_Poisson 9
float Tab_Poisson[MaxVal_in_Tab_Poisson][2] ={
                             {1e-1, 1.2816},
                             {1e-2, 2.3263},
                             {1e-3, 3.0902},
                             {1e-4, 3.7190},
                             {1e-5, 4.2649},
                             {1e-6, 4.7534},
                             {1e-7, 5.1993},
                             {1e-8, 5.6120},
                             {1e-9, 5.9978}};
			     
int main(int argc, char *argv[]) {


   // licence ....
   lm_check(LIC_MR1);
   
   // Read param IN
   nfe_Param o_Param;
   o_Param (argc, argv);   

   // create Few Event instance
   FewEventPoisson o_FewEventPoisson (o_Param.i_NbAutoConv);
   
   // set format type = f(format type of abaque file)
   //type_format FormatData = io_detect_format(o_Param.tc_AbaqueName);
   //io_set_format(FormatData);
   
   // compute distribution
   o_FewEventPoisson.compute_distribution (o_Param.e_WriteAllFiles);
   
   // compute threshold
   o_FewEventPoisson.find_threshold (o_Param.f_Epsilon, o_Param.e_WriteAllFiles);
   
   // write result
   io_write_ima_float (o_Param.tc_AbaqueName, o_FewEventPoisson._Threshold);

   return (0);
}





			     
inline void readint (int& pi_IntRead, char* ppc_msg) {
  if (sscanf(OptArg,"%d",&pi_IntRead) != 1) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit(-1);
  }
}

inline void readfloat (float& pf_FloatRead, char* ppc_msg) {
  if (sscanf(OptArg,"%f",&pf_FloatRead) != 1) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit(-1);
  }
}

inline Bool testborn (int pi_Min, int pi_Max, int pi_Val, char* ppc_msg) {
  if ((pi_Val<pi_Min) && (pi_Val>pi_Max)) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit (-1);
  }
  return (True);
}			     

inline Bool testborn (float pf_Min, float pf_Max, float pf_Val, char* ppc_msg) {
  if ((pf_Val<pf_Min) && (pf_Val>pf_Max)) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit (-1);
  }
  return (True);
}			     
			     
nfe_Param::nfe_Param () {   

   strcpy (tc_AbaqueName, DEF_ABAQUENAME);      /* abaque name */
   f_Epsilon =            DEF_EPSILON;          /* precision number */
   i_NbAutoConv =         DEF_NBR_AUTOCONV;         /* number max of event */
   e_WriteAllFiles =      False;                /* write all inter files */
   e_AllDefTrue =         False;                /* put all def at true */
   e_Trace =              False;                /* trace param */
   e_Verbose =            False;                /* verbose exec */
}			     
			     
			     
			     
void nfe_Param::operator () (int argc, char *argv[]) {

  int c;
  while ((c = GetOpt(argc,argv,"e:n:wdkv?")) != -1) {
	
    switch (c) {
    
    case 'e': /* -n <Epsilon> */
      readfloat (f_Epsilon, "bad epsilon: %f\n");
      testborn (0., 1., f_Epsilon, "bad epsilon: %s\n");
      break;    
    
    case 'n': /* -n <NbAutoConv> */
      readint (i_NbAutoConv, "bad autoconv number: %s\n");
      testborn (1, DEF_NBR_AUTOCONV, i_NbAutoConv, "bad autoconv number: %s\n");
      break;
   
    case 'd':
      e_AllDefTrue = True; break;
    case 'w':
      e_WriteAllFiles = True; break;
    case 'k':
      e_Trace = True; break;
    case 'v':
      e_Verbose = True; break;
    case '?':
      usage(argv);
    }
  }
	
  // get output file name */
  if (OptInd < argc) strcpy(tc_AbaqueName, argv[OptInd++]);
  else if (!e_AllDefTrue) usage(argv);
  
  // make sure there are not too many parameters  
  if (OptInd < argc) {
    fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  }
   
  if (e_Trace) trace ();
}			     
			     
			     
			     
void nfe_Param::usage (char *argv[]) {

  fprintf(OUTMAN, "Usage: %s options [output]\n\n", argv[0]);
  fprintf(OUTMAN, "   where options =  \n");  
  manline(); 
  
  fprintf(OUTMAN, "         [-e Epsilon]\n");
  fprintf(OUTMAN, "             Epsilon = precision for computing thresholds\n");
  fprintf(OUTMAN, "             default is %e \n\n", DEF_EPSILON);
  fprintf(OUTMAN, "           Correspondance for a large number of photons (Gaussian noise) \n");
  fprintf(OUTMAN, "             Espilon <--> NSigma\n");
  for (int i=0; i < MaxVal_in_Tab_Poisson; i++)
     fprintf(OUTMAN, "             %e <--> %f\n", Tab_Poisson[i][0], 
                     Tab_Poisson[i][1]);
  
  abaque_option_usage(DEF_ABAQUENAME);  manline(); 
  exit(-1);
}
		     
			     
			     
void nfe_Param::trace () {

   cout << "Param IN :" << endl;
   cout << " -- [-e:] : epsilon = " << f_Epsilon << endl;
   cout << " -- [-n:] : autoconv number : " << i_NbAutoConv << endl;   
   if (e_AllDefTrue)     cout << " -- [-d]  : all def true" << endl;   
   if (e_WriteAllFiles)  cout << " -- [-w]  : write all files" << endl;        
   if (e_Trace)   cout << " -- [-k]  : trace param" << endl;   
   if (e_Verbose)        cout << " -- [-v]  : verbose" << endl;

}
