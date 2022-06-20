 
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
#include "fractal.h"

extern int OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool xe_WriteFile=False;
float xf_MinQ=-10.00;
float xf_MaxQ=10.00;
int xi_NbQ=20;
int xi_IndMinScale=-1;
int xi_IndMaxScale=-1;
char stc_File_Name_Max[256];
char stc_File_Name_SupportMax[256];
char stc_NameImagOut[256];

Bool Verbose=False;

static void usage(char *argv[])
{

    manline();
    fprintf(OUTMAN, "Usage: %s options WaveMax SupportWaveMax FileNameOut\n\n", 
            argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

//    fprintf(OUTMAN, "         [-w]\n");
//    fprintf(OUTMAN, "              write file Z[.,s] for all q\n");
//    fprintf(OUTMAN, "              default is False \n");
//    manline();

    fprintf(OUTMAN, "         [-a q_minq]\n");
    fprintf(OUTMAN, "              value min of vector q\n");
    fprintf(OUTMAN, "              default is -10.00 \n");
    manline();

    fprintf(OUTMAN, "         [-b q_max]\n");
    fprintf(OUTMAN, "              value max of vector q\n");
    fprintf(OUTMAN, "              default is +10.00 \n");
    manline();

    fprintf(OUTMAN, "         [-c number_of_q]\n");
    fprintf(OUTMAN, "              number values of vector q\n");
    fprintf(OUTMAN, "              default is 20 \n");
    manline();

//    fprintf(OUTMAN, "         [-m ind_scale_min]\n");
//    fprintf(OUTMAN, "              ind scale min \n");
//    fprintf(OUTMAN, "              used only if -w is set, default is 30 \n");
//    manline();

//    fprintf(OUTMAN, "         [-M ind_scale_max]\n");
//    fprintf(OUTMAN, "              ind scale max \n");
//    fprintf(OUTMAN, "              used only if -w is set, default is scale_max-2 \n");
//    manline();

    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n");
    manline();
	
    fprintf(OUTMAN, "   Z(q,s) file name is Z_<FileNameOut>\n");
    fprintf(OUTMAN, "   q file name is q_<FileNameOut>\n");
    fprintf(OUTMAN, "   s file name is s_<FileNameOut>\n");
    manline();

    exit(-1);
}



static void init(int argc, char *argv[]) {
    
  int c;
   
  /* get options */
  while ((c = GetOpt(argc,argv,"a:b:c:m:M:wv")) != -1) {
    
    switch (c) {
    case 'v': Verbose = True; 
              break;
    case 'a': 		/* -a min_q */
      if (sscanf(OptArg,"%f",&xf_MinQ ) != 1) {
	fprintf(OUTMAN, "bad value of min of q: %s\n", OptArg);
	exit(-1);        
      }
      //cout << xf_MinQ << endl;
      break;

    case 'b': 		/* -b max_q */
      if (sscanf(OptArg,"%f",&xf_MaxQ ) != 1) {
	fprintf(OUTMAN, "bad value of max of q: %s\n", OptArg);
	exit(-1);        
      }
      //cout << xf_MaxQ << endl;
      break;

    case 'c': 		/* -c number_of__q */
      if (sscanf(OptArg,"%d",&xi_NbQ ) != 1) {
	fprintf(OUTMAN, "bad number of values of q: %s\n", OptArg);
	exit(-1);        
      }
      //cout << xi_NbQ << endl;
      break;

    case 'm': 		/* -m ind scale min for writen file */ 
      if (sscanf(OptArg,"%d",&xi_IndMinScale) != 1) {
	    fprintf(OUTMAN, "bad Prob: %s\n", OptArg);
	    exit(-1);        
      }
      break;

    case 'M': 		/* -m ind scale max for writen file */ 
      if (sscanf(OptArg,"%d",&xi_IndMaxScale) != 1) {
	    fprintf(OUTMAN, "bad Prob: %s\n", OptArg);
	    exit(-1);        
      }
      break;

    case 'w':
      xe_WriteFile = True; break;

    case '?':
      usage(argv); break;

    default: break;
    }
  }
  if ((xe_WriteFile == False) &&
            ( (xi_IndMinScale > 0) || (xi_IndMaxScale > 0)))
  {
     fprintf(OUTMAN, "Error: -m and -M option only valid is -w is set ...\n");
     exit(-1);
  }
  
  /* get optional input file names from trailing parameters and open files */
  if (OptInd < argc) strcpy(stc_File_Name_Max, argv[OptInd++]);
  else usage(argv);

  if (OptInd < argc) strcpy(stc_File_Name_SupportMax, argv[OptInd++]);
  else usage(argv);

  if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
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
     cout << ", number of point : " << xi_NbQ << "" << endl;
  }
  
  fltarray ao_Max, ao_SupportMax;
  type_1d_format  DatFormat = io_detect_1dformat(stc_File_Name_Max);
  if (DatFormat == F1D_FITS) 
  {
     fits_read_fltarr (stc_File_Name_Max, ao_Max);
     fits_read_fltarr (stc_File_Name_SupportMax, ao_SupportMax);
  }
  else
  {
     io_read2d_ascii(stc_File_Name_Max, ao_Max);
     io_read2d_ascii(stc_File_Name_SupportMax, ao_SupportMax);
  }
  
  if (xe_WriteFile) 
  {
    if (xi_IndMinScale == -1) xi_IndMinScale = 30;
    if (xi_IndMaxScale == -1) xi_IndMaxScale = ao_Max.ny()-1;
    
    if (Verbose == True)
    {
       cout << "ind scale min : " << xi_IndMinScale << ", ind scale max : ";
       cout << xi_IndMaxScale << endl;
    }
  }

  // Themodyn formalism
  // ------------------  
  fltarray ao_q(xi_NbQ);   
  for (int i=0; i<xi_NbQ; i++) 
    ao_q(i) = xf_MinQ + (xf_MaxQ-xf_MinQ)/xi_NbQ*i;

  to_ThermoDynRepartFunction ao_ThermoDyn (ao_Max, ao_SupportMax);
  ao_ThermoDyn.compute_thermo_partition(&ao_q, stc_NameImagOut);

  if (xe_WriteFile)
    ao_ThermoDyn.write_thermo_partition(xi_IndMinScale, xi_IndMaxScale);
  exit(0);
}
  




  





