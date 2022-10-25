/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.2
**
**    Author: 
**
**    Date:  98/07/30
**    
**    File:  mr1d_max
**
*******************************************************************************
**
******************************************************************************/

// static char sccsid[] = "@(#)mr1d_trans.cc 3.1 96/05/02 CEA 1994 @(#)";
 
#include "GlobalInc.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
 
char stc_NameImagIn[256];      /* input file image */
char stc_NameImagOut[256];     /* output file image */

int si_NbrPlan=0;              /* nb plan of decomp */
Bool se_SetPlan = False;  
Bool Verbose = False;     
type_trans_1d se_Transform=TO1_PAVE_B3SPLINE;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);



static void usage(char *argv[])
{

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();

    verbose_usage();
    manline();

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");

    exit(-1);
}




/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[]) {
    int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"n:t:v")) != -1) {
	
         switch (c) {
		case 't': /* -d <type> type of transform */
			if (sscanf(OptArg,"%d",&c ) != 1) {
		    		fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            		exit(-1);  
			}
                	if ((c > 0) && (c <= NBR_TRANS_1D)) se_Transform = (type_trans_1d) (c-1);
                	else {
		    		fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            		exit(-1);
 			}
			break;
		case 'n': /* -n <Nbr_Plan> */
			if (sscanf(OptArg,"%d",&si_NbrPlan) != 1) {
				fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
				usage(argv);
			}
			if ((si_NbrPlan <= 1) || (si_NbrPlan > MAX_SCALE_1D)) {
				fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
				fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
				usage(argv);
			}
			se_SetPlan=True;
			break;
	

		case 'v':
			Verbose = True;
			break;
		case '?':
			usage(argv);
		}
	}

	/* get optional input file names from trailing parameters and open files */
	if (OptInd < argc) strcpy(stc_NameImagIn, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
	else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}




/******************************************************************************/
static int max_scale_number (int pi_Nc) {

    int ai_Nmin = pi_Nc;
    int ai_ScaleMax;

	ai_ScaleMax = iround(log((float)ai_Nmin / (float)MAX_SIZE_LAST_SCALE) 
                                            / log(2.)+ 1.);
	return (ai_ScaleMax);
}

/******************************************************************************/



int main(int argc, char *argv[])
{
    extern type_1d_format IO_1D_Format;

  // local var
  fltarray  ao_Dat;

  /* Get command line arguments, open input file(s) if necessary */
  lm_check(LIC_M1D);
  transinit(argc, argv);
 
  // read the data
  // fits_read_fltarr (stc_NameImagIn, ao_Dat);
  io_1d_read_data(stc_NameImagIn, ao_Dat);
  reform_to_1d (ao_Dat);
  int ai_Np = ao_Dat.n_elem();

  // get nb Plan
  if (!se_SetPlan) si_NbrPlan = max_scale_number(ai_Np);

  fltarray ao_SignalIn(ai_Np);
  if (Verbose) {
    for (int i=0; i<ai_Np; i++) ao_SignalIn(i)=ao_Dat(i);
  }

  // transform
  MR_1D ao_Mr1dData (ai_Np, se_Transform, "Mr1d Transform", si_NbrPlan);
  ao_Mr1dData.Border = I_MIRROR;
  ao_Mr1dData.transform (ao_Dat);

  // trace
  if (Verbose) {
    cout << "Number of scales = " << ao_Mr1dData.nbr_scale() << endl;
    //fits_write_fltarr("WTransf", ao_Mr1dData.image());    
  }  
  
  // --------------------------------------------------- 
  // ---------------------------------------------------
  // compute local maximun of wavelet of IN signal
  ao_Mr1dData.loc_optima();
   if (IO_1D_Format == F1D_FITS)
        fits_write_fltarr(stc_NameImagOut, ao_Mr1dData.image());
   else io_write2d_ascii(stc_NameImagOut, ao_Mr1dData.image());

  exit(0);
}


