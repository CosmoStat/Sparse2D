 
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
int xi_TypeTransform = 0;
int xi_RemoveChainLength = 0;
float xf_RemoveDynRange = 0.0;
float xf_DistMinBetweenMax = 0.0;
char stc_NameImagOut[256];
char stc_NameImagIn[256];
Bool xe_DrawChain = False;
Bool xe_SortMax = False;
type_trans_1d xe_Transform = TO1_PAVE_DERIV_GAUSS;


static void usage(char *argv[])
{

    manline();
    fprintf(OUTMAN, "Usage: %s options signal_in\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

    fprintf(OUTMAN, "         [-t type of transform]\n");
    fprintf(OUTMAN, "              0 : TO1_PAVE_DERIV_GAUSS\n");
    fprintf(OUTMAN, "              1 : TO1_PAVE_MEX\n");
    fprintf(OUTMAN, "              default is TO1_PAVE_DERIV_GAUSS (0)\n");
    manline();    

//    fprintf(OUTMAN, "         [-p] \n");
//    fprintf(OUTMAN, "              sort maxima along line, take sup Max for\n");
//    fprintf(OUTMAN, "              all scale < current scale\n");
//    fprintf(OUTMAN, "              default is False\n");
//    manline();    

    fprintf(OUTMAN, "         [-r min_length_chain]\n");
    fprintf(OUTMAN, "              remove chains with length < min_leng..)\n");
    fprintf(OUTMAN, "              default is 0\n");
    manline();    

    fprintf(OUTMAN, "         [-s dynamique(%%)]\n");
    fprintf(OUTMAN, "              remove points of skelmap\n");
    fprintf(OUTMAN, "              if level < dynamique * max(current scale)\n");
    fprintf(OUTMAN, "              default is 0.\n");
    manline();    
 
//     fprintf(OUTMAN, "         [-c] \n");
//     fprintf(OUTMAN, "              draw chains of max\n");
//     fprintf(OUTMAN, "              default is False\n");
//     manline();    

    fprintf(OUTMAN, "         [-v] \n");
    fprintf(OUTMAN, "              verbose\n");
    fprintf(OUTMAN, "              default is False\n");
    manline();    

    fprintf(OUTMAN, "   Wavelet transform file name is WTMM_< signal_in>\n");
    fprintf(OUTMAN, "   Skeletton file name is WTMMskel_< signal_in>\n");

    exit(-1);
}

static void init(int argc, char *argv[]) {
    
  int c;
   
  /* get options */
  while ((c = GetOpt(argc,argv,"t:r:s:d:cpv")) != -1) {
    
    switch (c) {

    case 't': 		/* -t type of transform */ 
      if (sscanf(OptArg,"%d",&xi_TypeTransform) != 1) {
        fprintf(OUTMAN, "bad type of transform : %s\n", OptArg);
        exit(-1);        
      }
      if (xi_TypeTransform>1 || xi_TypeTransform<0) xi_TypeTransform=0;
      break;

    case 'r': 		/* -r min length of chain in support */ 
      if (sscanf(OptArg,"%d",&xi_RemoveChainLength) != 1) {
        fprintf(OUTMAN, "bad min chain length : %s\n", OptArg);
        exit(-1);        
      }
      //cout << xi_RemoveChainLength << endl;
      break;

    case 's': 		/* -s dynamique range */ 
      if (sscanf(OptArg,"%f",&xf_RemoveDynRange) != 1) {
        fprintf(OUTMAN, "bad min chain length : %s\n", OptArg);
        exit(-1);        
      }
      //cout << xf_RemoveDynRange << endl;
      break;
      
    case 'd': 		/* -d dist min between max */ 
      if (sscanf(OptArg,"%f",&xf_DistMinBetweenMax) != 1) {
        fprintf(OUTMAN, "bad dist min : %s\n", OptArg);
        exit(-1);        
      }
      //cout << xf_DistMinBetweenMax << endl;
      break;
      
    case 'p': 		/* -p sort Maxima line */ 
      xe_SortMax = True;
      break;

    case 'c': 		/* -c do not draw chain of maxima */ 
      xe_DrawChain = True;
      break;

    case 'v': 		/* -v verbose */ 
      xe_Verbose = True;
      break;

    case '?': usage(argv); break;

    default: break;
    }
  }  
    
  /* get optional input file names from trailing parameters and open files */
  if (OptInd < argc) strcpy(stc_NameImagIn, argv[OptInd++]);
  else usage(argv);

  //if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
  //else usage(argv);

  /* make sure there are not too many parameters */
  if (OptInd < argc) {
    fprintf (OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
    usage(argv);
  }
}





int main (int argc, char *argv[]) 
{
  // local var
  fltarray ao_Level;
  char atc_Wave[256];
  char atc_Skel[256];
  char *apc_ShortName, atc_Dir[256];
  extern type_1d_format IO_1D_Format;

  lm_check(LIC_M1D);
  init (argc, argv);
  if (xe_Verbose) 
  {
    if (!xe_DrawChain)
      cout << "compute support of maxima wavelet transform" << endl;
    else
      cout << "compute all chains lines in wavelet transform support" << endl;
  }

  // data read
  // ---------

  io_1d_read_data(stc_NameImagIn, ao_Level);
  reform_to_1d (ao_Level);
  int ai_Np = ao_Level.n_elem();
  if ((apc_ShortName = strrchr((char*)stc_NameImagIn, '/')) != (char*)NULL) {
     apc_ShortName = apc_ShortName + 1;
     for (int i=0; i<(int)strlen(stc_NameImagIn)-(int)strlen(apc_ShortName); i++)
        atc_Dir[i] = stc_NameImagIn[i];
     atc_Dir[(int)strlen(stc_NameImagIn)-(int)strlen(apc_ShortName)] = '\0';
     sprintf (atc_Wave, "%s/WTMM_%s", (char*)atc_Dir, apc_ShortName);
     sprintf (atc_Skel, "%s/WTMMskel_%s", (char*)atc_Dir, apc_ShortName);
  } else {
     sprintf (atc_Wave, "WTMM_%s", stc_NameImagIn);
     sprintf (atc_Skel, "WTMMskel_%s", stc_NameImagIn);
  }
  
  // init transform
  // --------------
  if (xi_TypeTransform == 0) xe_Transform=TO1_PAVE_DERIV_GAUSS;
  if (xi_TypeTransform == 1) xe_Transform=TO1_PAVE_MEX;


  // data transform + maxima
  // -----------------------

//TO1_PAVE_DERIV_GAUSS, TO1_PAVE_MEX
  if (xe_Verbose) 
     cout << "compute wavelet transform..." << endl;
  MR_1D ao_Mr1dData (ai_Np, xe_Transform, "");
  ao_Mr1dData.transform (ao_Level, I_MIRROR);// compute Mr1d transform
  ao_Mr1dData.loc_maxima();             // compute local maximum of Wave IN
  if (xe_Verbose) cout << "Number of scale : " << ao_Mr1dData.nbr_scale() 
                       << endl;
  fltarray ao_Temp;
  ao_Temp = ao_Mr1dData.image();
  if (IO_1D_Format == F1D_FITS) fits_write_fltarr(atc_Wave, ao_Temp);
  else io_write2d_ascii(atc_Wave, ao_Temp);

  // Compute chain of max
  // --------------------
  if (xe_Verbose) 
    cout << "compute support of maxima wavelet transform..." << endl;
  to_Skeleton ao_Skel (ao_Mr1dData, 1000);
  ao_Skel.Verbose = xe_Verbose;
  ao_Skel.remove_all_max_at_min_dist (xf_DistMinBetweenMax);
  ao_Skel.compute ();
  ao_Skel.remove_thin_skel (xi_RemoveChainLength);
  ao_Skel.remove_low_level (xf_RemoveDynRange);
  if (xe_SortMax) ao_Skel.take_max_along_inf_scale ();
  if (xe_Verbose) 
    cout << "Number of maxima line : " << ao_Skel.number_of_chain () << endl;
      
  if (!xe_DrawChain) {
    if (IO_1D_Format == F1D_FITS) 
         fits_write_fltarr(atc_Skel, ao_Skel.get_skel());
    else io_write2d_ascii(atc_Skel, ao_Skel.get_skel());
  } else {

    // write max chain
    // ---------------
    if (IO_1D_Format == F1D_FITS) 
         fits_write_fltarr(atc_Skel, ao_Skel.draw_lines (xe_Verbose));
    else io_write2d_ascii(atc_Skel, ao_Skel.draw_lines (xe_Verbose));
  }
  
  exit(0);
}
