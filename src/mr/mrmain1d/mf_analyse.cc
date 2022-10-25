 
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
#include "fractal.h"

extern int OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose = False;
Bool Optimize = False;

float xf_MinQ=-10.00;
float xf_MaxQ=10.00;
float xf_MinAlpha=0.0;
float xf_MaxAlpha=1.00;
int xi_NbAlpha=20;
int xi_IndMinScale=-1;
int xi_IndMaxScale=-1;
char stc_FileNameIn[256];


static void usage(char *argv[])
{

    manline();
    fprintf(OUTMAN, "Usage: %s options FileNameIn\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

    fprintf(OUTMAN, "         [-d alpha_minq]\n");
    fprintf(OUTMAN, "              value min of vector alpha\n");
    fprintf(OUTMAN, "              default is 0.00 \n");
    manline();

    fprintf(OUTMAN, "         [-e alpha_max]\n");
    fprintf(OUTMAN, "              value max of vector alpha\n");
    fprintf(OUTMAN, "              default is 1.00 \n");
    manline();

    fprintf(OUTMAN, "         [-f number_of_alpha]\n");
    fprintf(OUTMAN, "              number values of vector alpha\n");
    fprintf(OUTMAN, "              default is 20 \n");
    manline();

    fprintf(OUTMAN, "         [-m ind_scale_min]\n");
    fprintf(OUTMAN, "              ind scale min for computing tau\n");
    fprintf(OUTMAN, "              default is 30 \n");
    manline();

    fprintf(OUTMAN, "         [-M ind_scale_max]\n");
    fprintf(OUTMAN, "              ind scale max for computing tau\n");
    fprintf(OUTMAN, "              default is scale_max-2 \n");
    manline();
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n");
    manline();
    
    fprintf(OUTMAN, "   Z(q,s) are read from Z_<FileNameIn> file\n");
    fprintf(OUTMAN, "   q are read from q_<FileNameIn> file\n");
    fprintf(OUTMAN, "   s are read from s_<FileNameIn> file\n");
    manline();

    exit(-1);
}



static void init(int argc, char *argv[]) {
    
  int c;
   
  /* get options */
  while ((c = GetOpt(argc,argv,"d:e:f:m:M:vo")) != -1) {
    
    switch (c) {
    case 'v': Verbose = True; 
              break;
    case 'd': 		/* -d min_alpha */
      if (sscanf(OptArg,"%f",&xf_MinAlpha ) != 1) {
        fprintf(OUTMAN, "bad value of min of q: %s\n", OptArg);
        exit(-1);        
      }
      //cout << xf_MinAlpha << endl;
      break;

    case 'e': 		/* -e max_alpha */
      if (sscanf(OptArg,"%f",&xf_MaxAlpha ) != 1) {
        fprintf(OUTMAN, "bad value of max of q: %s\n", OptArg);
        exit(-1);        
      }
      //cout << xf_MaxAlpha << endl;
      break;

    case 'f': 		/* -f number_of_alpha */
      if (sscanf(OptArg,"%d",&xi_NbAlpha ) != 1) {
        fprintf(OUTMAN, "bad number of values of q: %s\n", OptArg);
        exit(-1);        
      }
      //cout << xi_NbAlpha << endl;
      break;

    case 'm': 		/* -m ind scale min in compute tau */ 
      if (sscanf(OptArg,"%d",&xi_IndMinScale) != 1) {
        fprintf(OUTMAN, "bad Prob: %s\n", OptArg);
        exit(-1);        
      }
      //cout << xi_IndMinScale << endl;
      break;

    case 'M': 		/* -m ind scale max in compute tau */ 
      if (sscanf(OptArg,"%d",&xi_IndMaxScale) != 1) {
        fprintf(OUTMAN, "bad Prob: %s\n", OptArg);
        exit(-1);        
      }
      //cout << xi_IndMaxScale << endl;
      break;
    case 'o':
      Optimize=True; break;
    case '?':
      usage(argv); break;

    default: break;
    }
  }

  /* get optional input file names from trailing parameters and open files */
  if (OptInd < argc) strcpy(stc_FileNameIn, argv[OptInd++]);
  else usage(argv);

  /* make sure there are not too many parameters */
  if (OptInd < argc) {
    fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
    exit(-1);
  }
}
 



int main (int argc, char *argv[]) {

  lm_check(LIC_M1D);
  init (argc, argv);

  if (Verbose == True)
  {
     cout << "min q : " << xf_MinQ << ", max q : " << xf_MaxQ;
     cout << "min alpha : " << xf_MinAlpha << ", max alpha : " << xf_MaxAlpha;
     cout << ", number of point : " << xi_NbAlpha << "" << endl;
  }

  // Themodyn formalism
  // ------------------  
  fltarray ao_alpha(xi_NbAlpha);   
  for (int i=0; i<xi_NbAlpha; i++) 
    ao_alpha(i) = xf_MinAlpha + (xf_MaxAlpha-xf_MinAlpha)/xi_NbAlpha*i;
  //type_1d_format  DatFormat = io_detect_1dformat(stc_FileNameIn);

  to_ThermoDynAnalysis ao_ThermoDyn (stc_FileNameIn);

  if (xi_IndMinScale == -1) xi_IndMinScale = 30;
  if (xi_IndMaxScale == -1) xi_IndMaxScale = ao_ThermoDyn.i_NbScale-1;
  
  if (Verbose == True)
  {
     cout << "ind scale min : " << xi_IndMinScale << ", ind scale max : ";
     cout << xi_IndMaxScale << endl;
  }
  ao_ThermoDyn.compute_tau(xi_IndMinScale, xi_IndMaxScale);
  ao_ThermoDyn.display_tau();
  ao_ThermoDyn.compute_frac_spectrum(&ao_alpha);
  ao_ThermoDyn.display_frac_spectrum();
  if (Optimize==True) ao_ThermoDyn.alpha_for_q_equals_zero();
  exit(0);
}
  




  





