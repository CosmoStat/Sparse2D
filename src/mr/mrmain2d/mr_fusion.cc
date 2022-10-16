/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: J.P. Djamdji & J.L. Starck
**
**    Date:  96/07/08
**    
**    File:  mr_fusion.cc
**
*******************************************************************************
**
**    DESCRIPTION   Image Registration
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_register option image_ref image_in image_out
**        where options = 
**
**           [-p]
**                Poisson noise
**                default is gaussian noise
**
**           [-g sigma]
**                Gaussian noise
**                  sigma = noise standard deviation 
**                by default, the noise is gaussian, and the standard 
**                devaition is automatically estimated. 
**
**           [-c gain,sigma,mean]
**                case of a CCD: noise = Poisson noise + read-out noise
**                  gain = CCD gain 
**                  sigma = standard deviation of the read-out noise
**                  mean = mean of the read-out noise
**                if this option is set, 
**                           Noise = Poisson + Gaussian read-out Noise
**                it is generally the case with the CCD.
**                Attention, these parameters must be separated by a comma 
**                without space. example: -c 0.133,1.733,0.
**                If mean, or sigma and mean are omitted, default values are 0.
**                gain can not be omitted. 
**
**           [-n number_of_scales]
**                number of scales used in the multiresolution transform
**                default is 4
**
**           [-s NSigma]
**                Thresolding at NSigma * SigmaNoise at each scale
**                default is 3
**
**           [-r res_min]
**                Miminum resolution for reconstruction
**                default is 1
**
**           [-d dist_max]
**                Maximum estimated distance between
**                two identical points in both images
**
**           [-l]
**                Registration choice:
**                    0:  Sub-scene registration 
**                    1:  Sub-scene and scene registration 
**                Default is sub-scene registration 
**
**           [-o]
**                Manual Options specifications:
**                    - Matching distance
**                    - Threshold level
**                    - Type of registration deformation model
**                    - Type of interpolation
**                    0:  None
**                    1:  Matching distance
**                    2:  Threshold level
**                    3:  Type of deformation model
**                    4:  Matching distance
**                        Threshold level
**                        Type of deformation model
**                    5:  Matching distance
**                        Threshold level
**                        Type of deformation model
**                        Type of interpolation
**                Default is None
**
**           [-w]
**                Write the following files
**                    - Deformation model parameters for each scale
**                    - control points for each scale
**                Default is None
**
**
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_fusion.cc 3.3 96/07/08 CEA 1996 @(#)";


#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "MR_Fusion.h"

char Name_Imag_In[80];  /* input file image to be warped */
char Name_Imag_Ref[80]; /* input reference file image  */
char Name_Imag_Out[80]; /* output warped file image */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=5;                 /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = TO_PAVE_BSPLINE; /* type of transform */

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
   
int ResMin = 1;
int ResMax = -1;
int DistMax = 0;
Bool WorkMinResol = False;
Bool WriteAll = False;

enum fusion_opt { DEF_OPT, MATCHING_LEVEL, THRESHOLD_LEVEL, REG_EQUA, MAN3_OPT, MAN4_OPT};
fusion_opt OptFusion = DEF_OPT;


#define LOW_RES 0
#define HIGH_RES 1
#define PARAM_CUBIC_INTERPOL -0.5
#define DEFAULT_INTERP T_REC_CUBIC_INTERPOL   
                          //  T_REC_CLOSER_PIXEL      ==> nearest neighbour
                          //  T_REC_BILINEAR_INTERPOL ==> bilinear
                          //  T_REC_CUBIC_INTERPOL    ==> cubic
                          
#define DEFAULT_DEFORM_MODEL T_EQUA_1   // T_EQUA_0: Pol first order type 1
 					// T_EQUA_1: Pol first order type 2
 					//  T_EQUA_2: Pol second order
 					//  T_EQUA_3 3: Pol third order

int OptImageType = LOW_RES;
int InterPolType = DEFAULT_INTERP;

int DeforModel = DEFAULT_DEFORM_MODEL;
Bool Verbose = False;

/*********************************************************************/


static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options ref_image in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    gauss_usage();
    manline();  
     
    poisson_noise_usage();
    manline();
    
    ccd_usage();
    manline();    

    nbr_scale_usage(Nbr_Plan);
    manline();
 
    nsigma_usage(N_Sigma);
    manline();
     
    fusion_option_usage(ResMin);
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
static void reginit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
    /* get options */
    while ((c = GetOpt(argc,argv,"wg:c:n:s:pr:D:l:o:i:d:vzZ:")) != -1) 
    {
	switch (c) 
        {
 	   case 'v': Verbose = True; break;
          // case 'p':
            //    WorkMinResol = True;
             //  break;
            case 'w':
                WriteAll = True;
               break;
            case 'p':
                /* Poisson noise */
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
               break;
            case 'i':
               if (sscanf(OptArg,"%d",&InterPolType) != 1) 
               {
		    fprintf(OUTMAN, "Error: bad interpolation type: %s\n", OptArg);
		    exit(-1);
	       }
	       if ( (InterPolType < 0) || (InterPolType > 2))
	       {
		    fprintf(OUTMAN, "Error: bad interpolation type: %s\n", OptArg);
                    fprintf(OUTMAN, "       it must be between 0 and 2.\n");
		    exit(-1);
	       }
	       break;
	    case 'd':
	       if (sscanf(OptArg,"%d",&DeforModel) != 1) 
               {
		    fprintf(OUTMAN, "Error: bad DeforModel: %s\n", OptArg);
		    exit(-1);
	       }
	       if ( ( DeforModel < 0) || ( DeforModel > 3))
	       {
		    fprintf(OUTMAN, "Error: bad deformation model: %s\n", OptArg);
                    fprintf(OUTMAN, "       it must be between 0 and 3.\n");
		    exit(-1);
	       }
	       break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		break;
             case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		break;
	    case 'r':
		/* -r <res_min> */
		if (sscanf(OptArg,"%d",&ResMin) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad minimum resolution: %s\n", OptArg);
		    exit(-1);
		}
		break;
            case 'D':
 		/* -d <dist_max> */
		if (sscanf(OptArg,"%d",&DistMax) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad maximum distance: %s\n", OptArg);
		    exit(-1);
		}
		break;
	    case 'o':
		/* -o <option> */
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad option: %s\n", OptArg);
		    exit(-1);
		}
		switch (c)
		{
		   case 0: OptFusion = DEF_OPT; break;
		   case 1: OptFusion = MATCHING_LEVEL; break;
		   case 2: OptFusion = THRESHOLD_LEVEL; break;
		   case 3: OptFusion = REG_EQUA; break;
		   case 4: OptFusion = MAN3_OPT; break;
		   case 5: OptFusion = MAN4_OPT; break;
		   default: fprintf(OUTMAN, "Error: bad option: %s\n", OptArg);
		            exit(-1);
		            break;
		}
		break;
	    case 'l': 
	        OptImageType = HIGH_RES; 
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

       /* get optional input file names from trailing 
          parameters and open files */

	if (OptInd < argc) strcpy(Name_Imag_Ref, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/
/*********************************************************************/
 
void get_param_polynomial(int res_min, int res_max, int *type_equa1)
{
     printf("\n");
     printf("Type of deformation model for the registration\n");
     printf("\n") ;
     printf("0    Polynomial of the first order and of type I \n") ;
     printf("1    Polynomial of the first order and of type II \n") ;
     printf("2    Polynomial of the second order\n") ;
     printf("3    Polynomial of the third order\n") ;
     printf("\n") ;

     for (int i = res_min ; i < res_max ; i++)
     {
        printf("resolution  :  %d            your choice  :  ",i+1) ;
        scanf("%d",&type_equa1[i]) ;
     }
}

/*********************************************************************/

void get_param_threshold(int res_min, int res_max, 
              float *niveau_seuil_ref, float *niveau_seuil_trv)
{
   int i;
   
   for (i = res_min; i < res_max ; i++) 
   {
      printf("\n");
      printf("Ref Image   res : %d     Threshold level[%d]  :  ", i+1, i);
      scanf("%g", &niveau_seuil_ref[i]);
   } 
   
   for (i = res_min; i < res_max ; i++) 
   {
       printf("\n");
       printf("In Image   res : %d     Threshold level[%d]  :  ", i+1, i);
       scanf("%g", &niveau_seuil_trv[i]);
   }
}

/*********************************************************************/
     
void get_param_distance(int res_min, int res_max, float *distance0)
{
   printf("Matching distance for the different wavelet scales  \n");
   for (int i = res_min; i < res_max ; i++)
   {
       printf("res : %d      distance[%d]  :  ", i+1, i);
       scanf("%f", &distance0[i]);
   }
}

/*********************************************************************/

void get_param_interpol(int &type_rec)
{
  printf("  \n");
  printf("Type of interpolation for the geometric registration\n");
  printf("  \n");
  printf("0    Zero order interpolation - nearest neighbour \n");
  printf("1    First order interpolation - bilineaire \n");
  printf("2    Second order interpolation - bicubique with parameter a \n");
  printf("  \n");
  printf("Your choice  :  ");
  scanf("%d", &type_rec);
}

/*********************************************************************/

void get_param_low_resol(int &depart_br_x, int &depart_br_y, 
                         float & pas_br_x, float &pas_br_y)
{
      printf("Low resolution image \n") ;
      printf("---------------------- \n") ;
      printf("\n") ;
      printf("starting coodinates in X    :   ") ;
      scanf("%d", &depart_br_x) ;
      printf("starting coodinates in Y    :   ") ;
      scanf("%d", &depart_br_y) ;
      printf("sampling step in X    :   ") ;
      scanf("%f", &pas_br_x) ;
      printf("sampling step in Y    :   ") ;
      scanf("%f", &pas_br_y) ;
}

/*********************************************************************/

void get_param_high_resol(char *name_entree_1, char *name_sortie_1, 
                     int & nb_lig1, int & nb_col1,
                     int &depart_x, int &depart_y, float &pas_x, float &pas_y)
                          
{

      printf("\n") ;
      printf("\n") ;
      printf("High resolution image \n") ;
      printf("---------------------- \n") ;
      printf("\n") ;
      printf("high resolution name file to register         :  ");
      scanf("%s",name_entree_1) ;
      sprintf(name_sortie_1,"%s_rec",name_entree_1) ;
      printf("\n") ;
      printf("\n") ;
      printf("number of lines             :   ") ;
      scanf("%d", &nb_lig1) ;
      printf("number of columns           :   ") ;
      scanf("%d", &nb_col1) ;
      printf("starting coordinates in X   :   ") ;
      scanf("%d", &depart_x) ;
      printf("starting coordinates in Y   :   ") ;
      scanf("%d", &depart_y) ;
      printf("sampling step in X   :   ") ;
      scanf("%f", &pas_x) ;
      printf("sampling step in Y   :   ") ;
      scanf("%f", &pas_y) ;  
} 

/*********************************************************************/

int main(int argc, char *argv[])
{
    int i,j,k, Nl, Nc;
    Ifloat Imag_Ref;
    Ifloat Imag_In;
    fitsstruct Header;
    char Cmd[256];

    Ifloat Imag_IN_HR;
    Ifloat Imag_Reg_HR;
    char Name_Imag_Reg_HR[256];
    char Name_Imag_IN_HR[256];


  float   *x1_e, *y1_e, *x2_e, *y2_e,
              *x1_s, *y1_s,
              *x00_s, *y00_s,
              *x10_s, *y10_s;

  float   *f1_e, *f2_e ;

  int         res_max, res,
              nmax, nmax1, nmax2, res_min, type_rec,
              depart_x, depart_y, depart_br_x=0, depart_br_y=0,
              nb_lig1, nb_col1;
  int         *type_equa1 ;
  float       *distance0 ;
  int	      dist_max;
  float       pas_x, pas_y, pas_br_x=0, pas_br_y=0,a;




  float       *niveau_seuil_ref ;
  float       *niveau_seuil_trv ;

  FILE        *file_param ;
  char	      name_mouch_param[1024] ;
  char        name_entree_1[256], name_sortie_1[256];

 parametre    par1, par1_cr ;

  char        name_incrus[31][256] ;



    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    reginit(argc, argv);

 if (Verbose == True)
 {
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Ref = " << Name_Imag_Ref << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "Noise Type = " << StringNoise(Stat_Noise) << endl;

    if (Stat_Noise == NOISE_GAUSS_POISSON)
    {
           cout << "Type of Noise = POISSON" << endl;
           cout << "  Gain = " << PasCodeur << endl;
           cout << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
           cout << "  Read-out Mean = " << MeanGauss << endl;
    }
    cout << "Sigma Noise = " << Noise_Ima << endl;
    cout << "N_Sigma = " << N_Sigma << endl;
 }



/*********************************************************************/
/***********  specifications of the option parameters   **************/
/*********************************************************************/

   res_max = Nbr_Plan - 1 ;
   res_min = ResMin -1;
   if (res_min < 0) res_min = 0;
   dist_max = DistMax ;
  
   if (dist_max != 0)
   {
     res_max = (int) (log((float)dist_max)/log(2.0) + 1.);
     Nbr_Plan = res_max + 1  ;
     cout << "Number of scales: " << Nbr_Plan << endl; 
   }
   if (res_min == res_max)
   {
      cerr << "Error: minimum resolution parameter must be smaller than the number of scales ... " << endl;
      exit(-1);
   }
       
   niveau_seuil_ref = new float[Nbr_Plan] ;
   niveau_seuil_trv = new float[Nbr_Plan] ;
   distance0    = new float[Nbr_Plan] ;
   type_equa1   = new int[Nbr_Plan] ;

   if (OptImageType != LOW_RES)
   {
      get_param_low_resol(depart_br_x, depart_br_y, pas_br_x, pas_br_y);
      get_param_high_resol(name_entree_1, name_sortie_1, 
                           nb_lig1,nb_col1,depart_x,depart_y, pas_x,pas_y);
   }

  // Default options
  for (i = 0; i <  Nbr_Plan; i++)
  {
      distance0[i] = (pow((float)2.0,(float)i+2)) + 1;
      if (i == 0) niveau_seuil_ref[i] = N_Sigma+1;
      else niveau_seuil_ref[i] = N_Sigma;
      niveau_seuil_trv[i] = niveau_seuil_ref[i];
      type_equa1[i] = DeforModel;
  }

  // By default, cubic iterpolation is taken
  type_rec = InterPolType;
  a =  PARAM_CUBIC_INTERPOL;

  if (OptFusion != DEF_OPT)
  {
     printf("\n");
     printf("PARAMETERS ENTRY:  \n");
     printf("\n");
  }
  
  if ((OptFusion == MATCHING_LEVEL)  
           || (OptFusion == MAN3_OPT) || (OptFusion == MAN4_OPT))
       get_param_distance(res_min, res_max, distance0);
 
  if ((OptFusion == THRESHOLD_LEVEL) 
           || (OptFusion == MAN3_OPT) || (OptFusion == MAN4_OPT))
      get_param_threshold(res_min, res_max, niveau_seuil_ref, niveau_seuil_trv);
 
  if ((OptFusion == REG_EQUA) 
           || (OptFusion == MAN3_OPT) || (OptFusion == MAN4_OPT))
      get_param_polynomial(res_min, res_max, type_equa1);
  
  if (OptFusion == MAN4_OPT)
      get_param_interpol(type_rec);
 
 if (WriteAll == True)
 {
     for (i = 0 ; i < res_max; i++)
         sprintf(name_incrus[i],"scale_%d_control_points.dat",i+1) ;
 }

/*********************************************************************/


    // read the reference image
    io_read_ima_float(Name_Imag_Ref, Imag_Ref, &Header);
    Header.origin = Cmd;
    Nl = Imag_Ref.nl();
    Nc = Imag_Ref.nc();

    // allocates the registered output image
    Ifloat Imag_Reg (Imag_Ref.nl(), Imag_Ref.nc(), "Result register");

    // read the input image
    io_read_ima_float(Name_Imag_In, Imag_In);
    check_scale(Imag_In.nl(), Imag_In.nc(), Nbr_Plan);



    /* --------------------------------------------------------- */
    /* --                Reference Image                      -- */
    /* --------------------------------------------------------- */
    // noise modelisation object for the reference image
    MRNoiseModel ModelDataRef(Stat_Noise, Imag_Ref.nl(), Imag_Ref.nc(), 
                              Nbr_Plan, Transform);
    ModelDataRef.CCD_Gain = PasCodeur;
    ModelDataRef.CCD_ReadOutSigma = SigmaGauss;
    ModelDataRef.CCD_ReadOutMean = MeanGauss;

    ModelDataRef.OnlyPositivDetect = True; // positive thresholding only
    for (i=0; i<res_max;i++)            // threshold specification for
      {                                    // noise detection at each scale
         ModelDataRef.NSigma[i] = niveau_seuil_ref[i] ;
      }

    // Multiresolution object declaration
    MultiResol MR_DataRef(Nl, Nc, Nbr_Plan, Transform, "Reference Image");

    ModelDataRef.model(Imag_Ref, MR_DataRef);
     
    if (Stat_Noise == NOISE_GAUSSIAN)
        cout << endl << "SigmaNoise of the reference image = " << ModelDataRef.SigmaNoise << endl;

    // Thresolholding of the wavelet transform
    ModelDataRef.threshold(MR_DataRef);


    /* --------------------------------------------------------- */
    /* --                  Input Image                        -- */
    /* --------------------------------------------------------- */
    // noise modelisation object for the reference image
    MRNoiseModel ModelDataIn(Stat_Noise, Imag_In.nl(), Imag_In.nc(), 
                             Nbr_Plan, Transform);
    ModelDataIn.CCD_Gain = PasCodeur;
    ModelDataIn.CCD_ReadOutSigma = SigmaGauss;
    ModelDataIn.CCD_ReadOutMean = MeanGauss;

    ModelDataIn.OnlyPositivDetect = True; // positive thresholding only
    for (i=0; i<res_max;i++)            // threshold specification for
      {                                    // noise detection at each scale
         ModelDataIn.NSigma[i] = niveau_seuil_trv[i] ;
      }

   // Multiresolution object declaration
   MultiResol MR_DataIn(Imag_In.nl(),Imag_In.nc(), Nbr_Plan, 
                        Transform, "Input Image");
    ModelDataIn.model(Imag_In, MR_DataIn);
    if (Stat_Noise == NOISE_GAUSSIAN)
        cout << "SigmaNoise of the second image = " << ModelDataIn.SigmaNoise << endl << endl;
 
   // Thresolholding of the wavelet transform
   ModelDataIn.threshold(MR_DataIn);


  /*======================================================================
    =====                                                            =====
    =====   Procedures call :                                        =====
    =====                      * Detection                           =====
    =====                      * Matching                            =====
    =====                      * Least squares                       =====
    =====                      * Registration                        =====
    =====                      * Real coordinates registration       =====
    =====                      * control points copy                 =====
    =====                                                            =====
    ======================================================================*/

    #if WRITE_PARAM
    printf("\n");
    printf("REGISTRATION PROCEDURE\n");
    #endif

  for (res = res_max-1; res >= res_min; res--) 
  {
    #if WRITE_PARAM
    printf("\n");
    printf("\n");
    printf("Working level scale: %d\n", res+1);

    printf("Detection procedure: Reference image \n");
    #endif
    detection(MR_DataRef.band(res), &x1_e, &y1_e, &f1_e, &nmax1);

    #if WRITE_PARAM
    printf("Detection procedure: Input image \n");
    #endif
    detection(MR_DataIn.band(res),  &x2_e, &y2_e, &f2_e, &nmax2);

    #if WRITE_PARAM
    cout <<   "Number of maxima in ImaRef = " << nmax1 << endl;
    cout <<   "Number of maxima in ImaIn = " << nmax2 << endl;
    #endif

    if (res != res_max-1) 
    {
      #if WRITE_PARAM
      printf("maxima coordinate transpose \n");
      #endif
      transpo_coor_max(x1_e, y1_e,&par1, nmax1,type_equa1[res+1], &x1_s, &y1_s);
    }
    else  
    {
      x1_s = new float[nmax1] ;
      y1_s = new float[nmax1] ;

      for (j = 0; j < nmax1; j++) 
      {
	x1_s[j] = x1_e[j];
	y1_s[j] = y1_e[j];
      }
    }
 

    if ( OptImageType == HIGH_RES && res == res_min)
     {
       #if WRITE_PARAM
       printf("matching procedure in real coodinates\n") ;
       printf("for the last scale\n") ;
       #endif
     }
    else
     #if WRITE_PARAM
     printf("matching procedure \n") ;
     #endif

    mise_en_correspondance_coor_reel(res,res_min, OptImageType,
                                     depart_br_x,depart_br_y,pas_br_x,pas_br_y,
                         	     x1_e, y1_e, x2_e, y2_e, x1_s, y1_s,
                          	     distance0[res],
			  	     nmax1, nmax2, type_equa1[res],
                          	     &x00_s, &y00_s, &x10_s, &y10_s,
                          	     &nmax,
				     &par1, 
				     &par1_cr);

    if (WriteAll == True)
   {

     #if WRITE_PARAM
     printf("Ground control points to file\n") ;
     #endif
     copie_pts_fichier(name_incrus[res],x00_s, y00_s, x10_s, y10_s, nmax) ;
    }
 
    if (nmax <= 2)
     {
       printf("Error: not enough detected points (less then 3 control points) ...\n") ;
       exit(-1);
     }


   if (res == res_min) 
   {
     #if WRITE_PARAM
     printf("\n") ;
     printf("\n") ;
     printf("\n") ;
     printf("\n") ;
     printf("Geometrical registration procedure\n");
     #endif
     recalage (&par1, type_equa1[res], type_rec, a,  Imag_In, Imag_Reg);

    // write registerd input file
    io_write_ima_float(Name_Imag_Out, Imag_Reg, &Header);


   if (WriteAll == True)
   {
    sprintf(name_mouch_param,"deform_model.txt") ;
    file_param = fopen(name_mouch_param, "w") ;
 
    if (type_equa1[res] ==  T_EQUA_0)
     {
       fprintf(file_param,"First order deformation model of type I\n\n");
       fprintf(file_param,"a = %g   b = %g\n",par1.ax,par1.bx) ;
       fprintf(file_param,"cx = %g  cy = %g\n",par1.cx, par1.cy) ;
     }

    if (type_equa1[res] ==  T_EQUA_1)
     {
       fprintf(file_param,"First order deformation model of type II\n\n");
       fprintf(file_param,"a = %g   b = %g  c  = %g\n",par1.ax,par1.bx,par1.cx); ;
       fprintf(file_param,"d = %g   e = %g  f  = %g\n",par1.ay,par1.by,par1.cy);
     }

    if (type_equa1[res] ==  T_EQUA_2)
     {
        fprintf(file_param,"Second order deformation model \n\n");
        fprintf(file_param,"a = %g   b = %g  c  = %g\n",par1.ax,par1.bx,par1.cx); ;
        fprintf(file_param,"d = %g   e = %g  f  = %g\n",par1.dx,par1.ex,par1.fx);
       
        fprintf(file_param,"g = %g   h = %g  i  = %g\n",par1.ay,par1.by,par1.cy); ;
        fprintf(file_param,"j = %g   k = %g  l  = %g\n",par1.dy,par1.ey,par1.fy);
     }

    if (type_equa1[res] ==  T_EQUA_3)
     {
       fprintf(file_param,"Third order deformation model \n\n");
       fprintf(file_param,"Coef of the first polynome are \n");
       fprintf(file_param,"%g\n",par1.ax) ;
       fprintf(file_param,"%g\n",par1.bx) ;
       fprintf(file_param,"%g\n",par1.cx) ;
       fprintf(file_param,"%g\n",par1.dx) ;
       fprintf(file_param,"%g\n",par1.ex) ;
       fprintf(file_param,"%g\n",par1.fx) ;
       fprintf(file_param,"%g\n",par1.gx) ;
       fprintf(file_param,"%g\n",par1.hx) ;
       fprintf(file_param,"%g\n",par1.ix) ;
       fprintf(file_param,"%g\n",par1.jx) ;

       fprintf(file_param,"Coef of the second polynome are \n");
       fprintf(file_param,"%g\n",par1.ay) ;
       fprintf(file_param,"%g\n",par1.by) ;
       fprintf(file_param,"%g\n",par1.cy) ;
       fprintf(file_param,"%g\n",par1.dy) ;
       fprintf(file_param,"%g\n",par1.ey) ;
       fprintf(file_param,"%g\n",par1.fy) ;
       fprintf(file_param,"%g\n",par1.gy) ;
       fprintf(file_param,"%g\n",par1.hy) ;
       fprintf(file_param,"%g\n",par1.iy) ;
       fprintf(file_param,"%g\n",par1.jy) ;
     } // endif

   fclose(file_param) ;
  } // end if writeall

     if (OptImageType == HIGH_RES)
      {
        // read the high resolution input image
        io_read_ima_float(Name_Imag_IN_HR, Imag_IN_HR, &Header);
        Header.origin = Cmd;
        nb_lig1 = Imag_IN_HR.nl();
        nb_col1 = Imag_IN_HR.nc();

        #if WRITE_PARAM
        printf("\n") ;
        printf("\n") ;
        printf("\n") ;
        printf("\n") ;
        printf("High resolution geometrical registration procedure\n");
        #endif
        recalage_coor_reelle(&par1_cr,
                             type_equa1[res], type_rec, a,
                             Imag_IN_HR, Imag_Reg_HR,
                             depart_x, depart_y, pas_x, pas_y);

        // write the high resolution registered input image 
        io_write_ima_float(Name_Imag_Reg_HR, Imag_Reg_HR, &Header);

      } // endif 

    } // end if res==resmin

    if (x1_e != NULL) delete [] x1_e;
    if (y1_e != NULL) delete [] y1_e;
    if (x2_e != NULL) delete [] x2_e;
    if (y2_e != NULL) delete [] y2_e;
    if (x1_s != NULL) delete [] x1_s;
    if (y1_s != NULL) delete [] y1_s;
    if (x00_s != NULL) delete [] x00_s;
    if (y00_s != NULL) delete [] y00_s;
    if (x10_s != NULL) delete [] x10_s;
    if (y10_s != NULL) delete [] y10_s;
    if (f1_e != NULL) delete [] f1_e;
    if (f2_e != NULL) delete [] f2_e;

   }

   if (WriteAll == True)
   {
     Ifloat grille(Imag_Reg.nl(),Imag_Reg.nc(),"grille in");
     Ifloat grille_out(Imag_Reg.nl(),Imag_Reg.nc(),"grille out");
     int ModVal =8;

     for (i=0; i < grille.nl(); i++)
     for (j=0; j < grille.nc(); j++)
     {
         grille(i,j) = 0;
         if ((i % ModVal == 0) || (j % ModVal == 0)) grille(i,j) = 1;

     }
     recalage (&par1, type_equa1[0], type_rec, a, grille, grille_out);
     io_write_ima_float("xx_grill_in.fits", grille);
     io_write_ima_float("xx_grill_out.fits", grille_out);
   }
   exit(0);
} 

