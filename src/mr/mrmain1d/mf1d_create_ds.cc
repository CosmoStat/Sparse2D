 
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
#include "fractal.h"


// extern var
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool xe_Verbose = False;
int xi_NbPoint=243;
int xi_NbInterval=3;
fltarray* xpo_Prob= new fltarray (xi_NbInterval);
int xi_NbPBorder=81;
char stc_NameImagOut[256];

static int si_Compt=0;


static void usage(char *argv[])
{

    manline();
    fprintf(OUTMAN, "Usage: %s options image out\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

    fprintf(OUTMAN, "         [-n number_of_points]\n");
    fprintf(OUTMAN, "              Number of points of DevilStairCase\n");
    manline();    

//     fprintf(OUTMAN, "         [-i number_of_intervals]\n");
//     fprintf(OUTMAN, "              number of intervals in DevilStairCase)\n");
//     manline();    
    
    fprintf(OUTMAN, "         [-p Prob1]\n");
    fprintf(OUTMAN, "             Probability P1. Default is 0.5.\n");
    manline();    
    
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "              Add a border. Default is no. \n");
    manline();    
    
    fprintf(OUTMAN, "         [-v] \n");
    fprintf(OUTMAN, "              verbose\n");
    fprintf(OUTMAN, "              default is False\n");
    manline();    

    exit(-1);
}

static void init(int argc, char *argv[]) {
    
  int c;
  Bool Optb = False;
  Bool Optp = False;
  float P1;
   
  /* get options */
  while ((c = GetOpt(argc,argv,"n:p:m:M:bv")) != -1) {
    
    switch (c) 
    {

    case 'n': 		/* -n number of points */
      if (sscanf(OptArg,"%d",&xi_NbPoint ) != 1) {
	fprintf(OUTMAN, "Error: bad number of Point: %s\n", OptArg);
	exit(-1);        
      }
       break;

//     case 'i': 		 
//       if (sscanf(OptArg,"%d",&xi_NbInterval ) != 1) {
// 	fprintf(OUTMAN, "bad number of Interval: %s\n", OptArg);
// 	exit(-1);        
//       }
//       xpo_Prob = new fltarray (xi_NbInterval);
//       break;

     case 'p': 	  
      if (sscanf(OptArg,"%f",&P1) != 1)
      {
          fprintf(OUTMAN, "Error: bad Prob: %s\n", OptArg);
          usage(argv);
      }
      if ( (P1 < FLOAT_EPSILON) || (P1 > 1. - FLOAT_EPSILON))
      {
          fprintf(OUTMAN, "Error: bad Probability ... \n");
	  fprintf(OUTMAN, "       P must verify:  0 < P < 1 \n");
          usage(argv);
      }
      (*xpo_Prob)(0) = P1;
      (*xpo_Prob)(1) = 0;
      (*xpo_Prob)(2) = 1. - P1;      
      si_Compt = 3;
      Optp = True;
      break;
    case 'b':
      Optb = True;
      break;
    case 'v': 		/* -v verbose */ 
      xe_Verbose = True;
      break;

    case '?': usage(argv); break;

    default: break;
    }
  }  
  if (Optb == True) xi_NbPBorder = xi_NbPoint / 3;
  else xi_NbPBorder = 0;
  
  if (Optp == False)
  {
     (*xpo_Prob)(0) = 0.5;
     (*xpo_Prob)(1) = 0.;
     (*xpo_Prob)(2) = 0.5;
     si_Compt = 3;
  }
  
  if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
  else usage(argv);

  /* make sure there are not too many parameters */
  if (OptInd < argc) {
    fprintf (OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  }
}





int main (int argc, char *argv[]) 
{
  lm_check(LIC_M1D);
  init (argc, argv);

  if (xe_Verbose == True)
  {
     cout << "number of points : " << xi_NbPoint+2*xi_NbPBorder << endl;
     cout << "number of intervals : " << xi_NbInterval << endl;
     for (int i=0; i<xi_NbInterval; i++)
     cout << " level of coef n°" << i+1 << " : " << (*xpo_Prob)(i) << endl;
     cout << "number of points of each border : " << xi_NbPBorder << endl;
     cout << "create DevilStairCase function..." << endl;
  }

  // data initialisation
  // -------------------

  to_GenDevilStairCase ao_class (xi_NbPoint, xi_NbPBorder, xpo_Prob);
  xi_NbPoint = xi_NbPoint + 2*xi_NbPBorder;
  fltarray ao_Level(xi_NbPoint);         // signal for multiresol analysis
  ao_Level = ao_class.get_level();       // get signal form Devil.. classe

  // if (xe_Verbose== True) 
  //  for (int i=0;i<xi_NbPoint;i++) cout << ao_Level(i) << " ";

  io_1d_write_data(stc_NameImagOut, ao_Level);
  exit(0);
}
