/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 
**
**    Author: P. Querre
**
**    Date:  8/8/1998
**    
**    File:  mr_rec_edge.cc
**
*******************************************************************************
**
**    DESCRIPTION  reconstitution from modulus maxima
**    ----------- 
**                 
**    Usage: mr2d_recons_max options image output
**
**   where options =
**
** 
**         [-k number of iteration]
**              number of iteration used in reconstruction
**
**         [-v]
**				Verbose 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Edge.h"

#define DEF_ITER 10
#define OK 0
#define NOT_OK 1
#define EPS 1e-07
#define MAX_ITER 100
#define MAX_SCALE 10

// static var
char stc_NameImagIn[256];
char stc_NameImagOut[256];
int  si_MaxIter=DEF_ITER;
Bool Verbose = False; 



// extern var
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);



/*********************************************************************/
static void usage(char *argv[])
{
    fprintf (OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf (OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-i number_of_iterations]\n");
    fprintf(OUTMAN, "              number of iteration\n");
    fprintf(OUTMAN, "              default is %d\n", DEF_ITER);
    manline();
    vm_usage();
    manline();
    verbose_usage();
    manline();

    manline();
    exit(-1);
}

/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
	int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
   
	/* get options */
	while ((c = GetOpt(argc,argv,"i:vzZ:")) != -1) {

		switch (c) {

		case 'i': /* -i < Number of iterations> */
			if (sscanf(OptArg,"%d",&si_MaxIter) != 1) {
				fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
				usage(argv);
			}
			if (si_MaxIter <= 0)  si_MaxIter = DEF_ITER;
			break;

		case 'v':
			Verbose = True;
			break;

#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
            case '?': usage(argv); break;
	    default: usage(argv); break;
  		}
	}

	/* get optional input file names from trailing parameters and open files */
	if (OptInd < argc) strcpy(stc_NameImagIn, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
	else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf (OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage (argv);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}





/******************************************************************************/

 

int main(int argc, char *argv[])
{

  // inter reconstruction file

  /* Get command line arguments, open input file(s) if necessary */
  lm_check(LIC_MR1);
  transinit(argc, argv);


  //---------------------------------------------------------------------------
  //--------- IMAGE IN --------------------------------------------------------
  //---------------------------------------------------------------------------
 
  // create multiresol struct of image IN (init)  
  // and read multiresolution IN (format Ifloat)
  // -------------------------------------------
  MultiResol ao_Mr2dData;
  ao_Mr2dData.read(stc_NameImagIn);
  int si_NbrPlan = ao_Mr2dData.nbr_scale();
  int ai_Nl = ao_Mr2dData.size_ima_nl();
  int ai_Nr = ao_Mr2dData.size_ima_nc();

  if (Verbose==True) {
	cout << "Number of band   = " << ao_Mr2dData.nbr_band() << endl;
    cout << "Number of scales = " << si_NbrPlan << endl;
    cout << "Number of line   = " << ai_Nl << endl;
    cout << "Number of row    = " << ai_Nr << endl; 
  } 


  // save Modulus and Phase
  // ----------------------
  Ifloat* apo_ModData = new Ifloat [si_NbrPlan-1];
  Ifloat* apo_PhaData = new Ifloat [si_NbrPlan-1];
  int s;
  for (s=0; s<si_NbrPlan-1; s++) {
      apo_ModData[s].alloc (ai_Nl, ai_Nr, "");
      apo_ModData[s] = ao_Mr2dData.extract_band(2*s);
      apo_PhaData[s].alloc (ai_Nl, ai_Nr, "");
      apo_PhaData[s] = ao_Mr2dData.extract_band(2*s+1);      
  }
  
  
  //---------------------------------------------------------------------------
  //--------- SEARCH IND MAX --------------------------------------------------
  //---------------------------------------------------------------------------
  // search coord (s,i,j) of all max max modulus 
  // vector ao_NumberMaxInLine    : number of max on line i at scale s
  // vector ao_NumberMaxInRow     : number of max on row j at scale s
  // intarray attpo_MaxModInLine  : coord of all max on line i at scale s
  // intarray attpo_MaxModInRow   : coord of all max on row j at scale s
  intarray ao_NumberMaxInLine (si_NbrPlan, ai_Nl);
  intarray ao_NumberMaxInRow (si_NbrPlan, ai_Nr);
  intarray* attpo_MaxModInLine [si_NbrPlan][ai_Nl];
  intarray* attpo_MaxModInRow [si_NbrPlan][ ai_Nr];  
  search_number_max (apo_ModData, ao_NumberMaxInLine, ao_NumberMaxInRow, si_NbrPlan);
  init_coord_max (apo_ModData, ai_Nl, ai_Nr, ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine,
                                             ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, si_NbrPlan);


  // compute level of Wavelet Transform
  // ----------------------------------
  for (s=0; s<si_NbrPlan-1; s++) {

    for (int i=0; i<ai_Nl; i++) {
       for (int j=0; j<ai_Nr; j++) {

         float x = ao_Mr2dData(2*s, i, j) * cos (ao_Mr2dData(2*s+1, i, j));
         float y = ao_Mr2dData(2*s, i, j) * sin (ao_Mr2dData(2*s+1, i, j));
         ao_Mr2dData(2*s, i, j) = x;
	 ao_Mr2dData(2*s+1, i, j) = y;
      }
    }
  }

  //ao_Mr2dData.write("WTransfRec");



   

  
 

  // Init Mutltiresol OUT
  // --------------------
  MultiResol ao_Mr2dRecData (ai_Nl, ai_Nr, si_NbrPlan, TO_DIADIC_MALLAT, "Mr2d Signal OUT");
  init_multiresol_out (ai_Nl, ai_Nr, si_NbrPlan, ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine,
                       ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, ao_Mr2dData, ao_Mr2dRecData);

  // init Multiresol IN
  // ------------------
  int b;
  for (b=0; b<ao_Mr2dData.nbr_band(); b++) 
    for (int i=0; i<ai_Nl; i++) 
      for (int j=0; j<ai_Nr; j++)
        ao_Mr2dData (b,i,j) = ao_Mr2dRecData (b,i,j);

   // ao_Mr2dData.write("trace.fits");






  // construct Ifloat[] modulus and phase struct 
  // -------------------------------------------
	// => OK deja fait
  
  //---------------------------------------------------------------------------
  //--------- SEARCH IND MAX --------------------------------------------------
  //---------------------------------------------------------------------------
  // search coord (s,i,j) of all max max modulus 
  // vector ao_NumberMaxInLine    : number of max on line i at scale s
  // vector ao_NumberMaxInRow     : number of max on row j at scale s
  // intarray attpo_MaxModInLine  : coord of all max on line i at scale s
  // intarray attpo_MaxModInRow   : coord of all max on row j at scale s
   	// => OK deja fait

  //---------------------------------------------------------------------------
  //--------- RECONSTRUCT IMAGE -----------------------------------------------
  //---------------------------------------------------------------------------

  // init recons signal struct (format Ifloat)
  Ifloat ao_RecData(ai_Nl, ai_Nr, "Recons");// used in recons of ao_Mr2dRecData

  // create multiresol struct of reconstruct image (init)
  	// => OK deja fait

  // raz of all coef of recons signal (except last scale)
  for (b=0; b<2*(si_NbrPlan-1); b++)
    for (int i=0; i<ai_Nl; i++) 
      for (int j=0; j<ai_Nr; j++)
  	ao_Mr2dRecData (b,i,j) = 0.0;

  // int last scale
  init_last_scale (si_NbrPlan, ai_Nl, ai_Nr, ao_Mr2dData, ao_Mr2dRecData);

  // --------------------------------------------------------------------------
  // --------- ITERATION FOR RECONSTRUCTION -----------------------------------
  // --------------------------------------------------------------------------
  // iteration for recons
  if (Verbose==True) cout << "number of iteration:" << si_MaxIter << endl;
  for (int iter=0; iter<si_MaxIter; iter++){

    if (Verbose==True) cout << "iteration:" << iter+1 << endl;;

    // --------------------------------------------------- 
    // ---------------------------------------------------
    interpolate (ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine, 
                 ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, 
                 ao_Mr2dData, ao_Mr2dRecData);

    // modify the modulus maximum of Mr2dRecDat : 
    // max (Mr2dRecData(s,i)) =  max (Mr2dData(s,i))
    init_max (si_NbrPlan, ai_Nl, ai_Nr, 
	      ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine,
	      ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, 
	      ao_Mr2dData, ao_Mr2dRecData);

    // int last scale
    init_last_scale (si_NbrPlan, ai_Nl, ai_Nr, ao_Mr2dData, ao_Mr2dRecData);

    // --------------------------------------------------- 
    // ---------------------------------------------------
    // invers and foward transform of Mr2dRecDat
    ao_Mr2dRecData.recons (ao_RecData);
//     if (Verbose==True) {
//       char atc_FileName[80];
//       sprintf (atc_FileName, "tfic_%d", iter);
//       io_write_ima_float (atc_FileName, ao_RecData);
//     }
    ao_Mr2dRecData.transform (ao_RecData);

  }
      
  // --------------------------------------------------- 
  // ---------------------------------------------------
  // modify the modulus maximum of Mr2dRecDat : 
  // max (Mr1dRecData(s,i,j)) =  max (Mr1dData(s,i,j))
  init_max (si_NbrPlan, ai_Nl, ai_Nr, 
	    ao_NumberMaxInLine, (intarray**)attpo_MaxModInLine,
	    ao_NumberMaxInRow,  (intarray**)attpo_MaxModInRow, 
        ao_Mr2dData, ao_Mr2dRecData);

  // int last scale
  init_last_scale (si_NbrPlan, ai_Nl, ai_Nr, ao_Mr2dData, ao_Mr2dRecData);

  // reconstruction IN omage
  ao_Mr2dRecData.recons (ao_RecData);
  io_write_ima_float (stc_NameImagOut, ao_RecData);


 
  exit(0);


} 

