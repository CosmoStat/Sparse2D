/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Hubert Druesne, Philippe Querre
**
**    Date:  02/06/2003 
**    
**    File:  cb_mca1d.cc
**
*******************************************************************************
**
**    DESCRIPTION  Decomposition of a signal on multiple bases
**    ----------- 
**                 
******************************************************************************/


#include "GlobalInc.h"
 #include "Usage.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Sigma.h"
#include "MR1D_Obj.h"
#include "MR_SoftRegul.h" 
#include "IM1D_Dct.h"
#include "IM1D_ALDCT.h"
#include "MDCT1D.h"
#include "MCA1D.h" 
#include <vector>
// la classe vector est la classe de la STL, pour 
// l'utiliser tu dois avoir un #include <vector>
// dans le fichier qui l'utilise (float c'est pour indiquer 
// un vecteur de float)
// VectMax.push_back : ajoute simplement un  élément dans le vecteur VectMax,,,,
// AWT.getAbsMaxTransf (AT_Trans) : l'élément rajouté est la valeur abs du max de la
// tarnsformé de AT_Trans
// 
// si le pb persiste, tu peux remplacer
// vector<float> VectMax;
// par fltarray VectMax (nombre de transformées),
// puis
// VectMax.push_back( AWT.getAbsMaxTransf (AT_Trans));
// par VectMax [AWT] = AWT.getAbsMaxTransf (AT_Trans)
// de meme pour les autres transf


char nameSignalIn[MAXCHAR];  // input file image  
char nameSignalOut[MAXCHAR]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool RemoveSmoothPlane=False; 
int NbRemoveSmoothPlane=MAX_SMOOTH_REMOVE;
Bool Write=False; 
int NbrScale1D = 0;

float Noise_Sig=-1;
Bool PoissonNoise=False;

int CosBlockSize=256;
int Nbr_Iter = DEF_MCA1D_NBR_ITER;
Bool COS_Overlapping=False;
Bool MCOS_Overlapping=False;
Bool SuppressDctCompCont = False;
float COS_Sensibility = 1.;
float COSMin = 0.;
Bool Linear = True;
Bool CosWeightFirst=False;
int FirstDetectScale=0;

Bool KillLastScale = False;

Bool SuppressIsolatedPixel = False;
Bool OnlyPositivDetect = False;
Bool PositivSol=False;
Bool UseMask =False;
Bool TabBase[NBR_MCA1D];
type_mca1dbase TabSelect[NBR_MCA1D];
int NbrBase = 0;

int NbrUndecimatedScale = -1;
Bool WriteAllRec=True;
float FirstSoftThreshold= -1;
float LastSoftThreshold = DEF_MCA1D_LAST_SOFT_THRESHOLD;
Bool TV = False;
Bool UsePsf = False;
float LambdaTV= 0.;
Bool  UseNormL1 = False;

float SoftLevel = 1.7;

int Infocost = 0;
Bool AdaptDCT=False;
Bool NonOrthoFilterBank = False;

float CFAR_NsigmaMad=5.;
Bool CFAR=False;       // If true, the threshold at each iteration is derived from
                       // T= MAD(Coef)*CFAR_NsigmaMad
Bool CFDR=False;       // If true, the threshold at each iteration is derived from
                       // from the FDR. 
float CFDR_Qparam=0.25;

/*********************************************************************/

int maxscalenumber (int Nc) 
{
    int Nmin = Nc;
    int ScaleMax;
    
    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}

/***********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_signal result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t TransformSelection]\n");
    fprintf(OUTMAN, "             1: A Trous transform \n");    
    fprintf(OUTMAN, "             2: Half decimated WT \n");  
		fprintf(OUTMAN, "             3: Local Discrete Cosinus Transform\n");    
    fprintf(OUTMAN, "             4: Multiscale Local Discrete Cosinus Transform \n");  
    fprintf(OUTMAN, "             5: Adaptive Local Discrete Cosinus Transform \n");  		
    fprintf(OUTMAN, "             Default is 1. \n");    
    manline(); 

    fprintf(OUTMAN, "         [-H]\n");
    fprintf(OUTMAN, "             Data contained masked area (must have a zero value). Default is no. \n");    
    manline();
   
	  fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "             Non-Linear step descent. Default is Linear. \n");    
    manline();


    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "            Replace the local DCT by the adaptive local DCT.\n");
    manline();
    
    fprintf(OUTMAN, "         [-S FirstThresholdLevel]\n");
    fprintf(OUTMAN, "            First thresholding value.\n");
    fprintf(OUTMAN, "            Default is automatically found.\n");    
    manline();

    fprintf(OUTMAN, "         [-s LastThresholdLevel]\n");
    fprintf(OUTMAN, "            Last thresholding value..\n");
    fprintf(OUTMAN, "            Default is %f.\n", LastSoftThreshold);    
    manline();
   
    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration. Default is %d.\n", Nbr_Iter);    
    manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Minimize the L1 norm. Default is no.\n");    
    manline();
    gauss_usage();
    manline();  
    
    fprintf(OUTMAN, "         [-l]\n");
    fprintf(OUTMAN, "             Supression of the smooth plane \n");
    fprintf(OUTMAN, "             calculations. Default is no. \n");    
    manline();
 
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales in the Wavelet Transform.\n");
    fprintf(OUTMAN, "             Default is automatically estimated.\n");    
    manline();     

    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Positivity constraint on the wavelet component.\n");
    fprintf(OUTMAN, "             Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-F FirstDetectScale]\n");
    fprintf(OUTMAN, "             First used scale in the WT. Default is 1.\n");
    manline();
    
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "             Suppress isolated pixels in the WT.\n");
    fprintf(OUTMAN, "             Default is false.\n");    
    manline(); 
                
    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Block overlapping in the DCT.\n"); 
    fprintf(OUTMAN, "             Default is no. \n");    
    manline();
  
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Detect only positive wavelet coeff.\n");
    fprintf(OUTMAN, "             Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-B CosinusBlockSize]\n");
    fprintf(OUTMAN, "             Cosinus transform block size. \n");    
    fprintf(OUTMAN, "             Default is %d.\n", CosBlockSize);    
    manline(); 
  
    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "            Kill last scale in the atrous algorithm\n");
    fprintf(OUTMAN, "            Default is no.\n");    
    manline();  

    fprintf(OUTMAN, "         [-h]\n");
    fprintf(OUTMAN, "            Suppress the zero frequency in the DCT.\n");
    fprintf(OUTMAN, "            Default is false.\n");    
    manline();
                      
//     manline();    
//     fprintf(OUTMAN, "   HDWT settings\n");
//     fprintf(OUTMAN, "   -------------\n");
//     manline(); 
//      
//     fprintf(OUTMAN, "         [-u]\n");
//     fprintf(OUTMAN, "             Number of undecimated scales.\n");
//     fprintf(OUTMAN, "             default is all scales. \n");    
//     manline();
//     
      
    fprintf(OUTMAN, "         [-I TypeCriter]\n");
    fprintf(OUTMAN, "             Criterion for the adaptive local DCT: \n");
    fprintf(OUTMAN, "             1: Entropy\n");    
    fprintf(OUTMAN, "             2: Norm L1\n");    
    fprintf(OUTMAN, "             Default is Entropy. \n");    
    manline();
    
    vm_usage();
    manline();    
    verbose_usage();
    manline();          
    manline();    
    exit(-1);
}
   

/***********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[]) 
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c= GetOpt 
            (argc,argv,"HQ:q:cCUNAWLB:Ps:S:wu:t:KhF:Oi:g:n:vlR:zZkpD:I:")) != -1){
	switch (c){	 
   case 'H': UseMask = True; break; 
	 case 'C': CFDR=True; break;
	 case 'c': CFAR=True; break;
	 case 'Q':  if (sscanf(OptArg,"%f",&CFDR_Qparam) != 1)
	            {
                       fprintf(OUTMAN, "Error: bad FDR parameter: %s\n",  OptArg);
                       exit(-1);
                    }
		break;
	 case 'q': if (sscanf(OptArg,"%f",&CFAR_NsigmaMad) != 1)
	            {
                       fprintf(OUTMAN, "Error: bad CFAR parameter: %s\n",  OptArg);
                       exit(-1);
                    }
		break;
	 case 'U': NonOrthoFilterBank = True; break;
	 case 'N': UseNormL1 = True; break;
	 case 'A': AdaptDCT = (AdaptDCT == True) ? False: True; break;
	 case 'L': Linear = (Linear == True) ? False: True; break;
	 case 'W': CosWeightFirst = (CosWeightFirst == True) ? False: True; 
                   break;
// 	 case 'C':
//                 if (sscanf(OptArg,"%f",&COSMin) != 1){
//                     fprintf(OUTMAN, "bad cosinus sensibility: %s\n", OptArg);
//                     exit(-1);
//                 }
//                 break;
          case 'k': SuppressIsolatedPixel = True; break;
          case 'h': SuppressDctCompCont = True; break;
          case 'S':               
                if (sscanf(OptArg,"%f",&FirstSoftThreshold) != 1){
                    fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", 
                    OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0){
                   fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", 
                   OptArg);
                   exit(-1);
                }
                break;
          case 's':               
                if (sscanf(OptArg,"%f",&LastSoftThreshold) != 1){
                    fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", 
                    OptArg);
                    exit(-1);
                }
                if (LastSoftThreshold < 0){
                   fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", 
                   OptArg);
                   exit(-1);
                }
                break;
          case 'P': PositivSol = True; break;
          case 'p': OnlyPositivDetect = True; break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1){
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", 
                    OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_MCA1D)){
                   if (TabBase[c-1] == True){
                       fprintf(OUTMAN, "Error: transform already selected: %s\n", 
                       OptArg);
                       exit(-1);
                   }
                   else{
                      TabBase[c-1] = True;
                      TabSelect[NbrBase++] = (type_mca1dbase) (c-1);
                   }
                }
                else{
                   fprintf(OUTMAN, "Error: bad type of transform: %s\n", 
                   OptArg);
                   exit(-1);
                }
                break;
           case 'I': 
	        if (   (sscanf(OptArg,"%d",&Infocost) !=1) 
                    || (Infocost<1) || (Infocost>2))
		{
                   fprintf(OUTMAN, "Error: bad type of cost for ALDCT: %s\n", 
                   OptArg);
                   exit(-1);
		}
		Infocost --;
		break;	
           case 'K': KillLastScale = True; break;
           case 'u': 
                if (sscanf(OptArg,"%d",&NbrUndecimatedScale) != 1){
                    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrUndecimatedScale < 0){
                   fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'F': 
                if (sscanf(OptArg,"%d",&FirstDetectScale ) != 1){
                    fprintf(OUTMAN, "Error: bad first detection scale: %s\n", 
                    OptArg);
                    exit(-1);
                }
                break;
           case 'O': COS_Overlapping =( COS_Overlapping== True) ? False: True; 
                     break; 
           case 'i':
                if (sscanf(OptArg,"%d",&Nbr_Iter) != 1){
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'B':
                if (sscanf(OptArg,"%d",&CosBlockSize) != 1){
                    fprintf(OUTMAN, "bad cosinus transform block size: %s\n", 
                    OptArg);
                    exit(-1);
                }
                break;
           case 'D':
                if (sscanf(OptArg,"%f",&COS_Sensibility) != 1){
                    fprintf(OUTMAN, "bad cosinus sensibility: %s\n", OptArg);
                    exit(-1);
                }
                break;
	         case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Sig) != 1){
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
           case 'v': Verbose = True;break;
           case 'w': Write = True;break;
           case 'l': RemoveSmoothPlane= True; break;
           case 'R': 
                if (sscanf(OptArg,"%d",&NbRemoveSmoothPlane) != 1){
                    fprintf(OUTMAN, "bad number of number of smooth remove : %s\n", OptArg);
                    exit(-1);
                }
                break;
	         case 'n':
                if (sscanf(OptArg,"%d",&NbrScale1D) != 1){
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale1D  > MAX_MGA_SCALE_1D){
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_MGA_SCALE_1D);
                    exit(-1);
                }
                break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True){
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1){
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True){
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
       if (OptInd < argc) strcpy(nameSignalIn, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(nameSignalOut, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc){
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/***************************************************************************/

int main(int argc, char *argv[]) 
{
    int i,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    
    MCA1D MB;  // Multiple base decomposition Class
    
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR3);
    
    for (int i=0; i<NBR_MCA1D; i++) TabBase[i] = False;

    filtinit(argc, argv);
    
    if (NbrBase == 0) 
    {
       NbrBase = 2;
       TabSelect[1] = MCA1D_ATROU;
       if (AdaptDCT == False) TabSelect[0] = MCA1D_COS;
       else TabSelect[0] = MCA1D_ALDCT;
    }
    
    fltarray Tab[NbrBase];
    
    fltarray Data;
    io_1d_read_data(nameSignalIn, Data, &Header);
    reform_to_1d(Data);
    Header.origin = Cmd;
    int Nx = Data.nx();
    
    if (NbrScale1D == 0) NbrScale1D=maxscalenumber(Nx);
    if (NbrScale1D > MAX_MGA_SCALE_1D) NbrScale1D = MAX_MGA_SCALE_1D;
    
    fltarray Result(Nx,"result");
    MB.NbrBase = NbrBase;
    MB.Verbose = Verbose;
    MB.Nbr_Iter = Nbr_Iter;
    MB.PositivRecSig = PositivSol;
    MB.Nx=Nx;
    MB.Linear = Linear;
    MB.AT_PositivRecSig = PositivSol;
    MB.AT_NbrScale1D = NbrScale1D;
    MB.AT_FirstDetectScale=FirstDetectScale;
    MB.AT_KillLastScale=KillLastScale;
    MB.AT_OnlyPositivDetect = OnlyPositivDetect;
    MB.AT_SuppressIsolatedPixel = SuppressIsolatedPixel;
    
    MB.WT_FirstDetectScale=FirstDetectScale;
    MB.WT_NbrScale1D = NbrScale1D;
    MB.WT_PositivRecSig= PositivSol;
    MB.WT_OnlyPositivDetect = OnlyPositivDetect;
    MB.WT_SuppressIsolatedPixel = SuppressIsolatedPixel;
    MB.WT_NonOrthoFilterBank = NonOrthoFilterBank;
    // MB.WT_KillLastScale=KillLastScale;
    
    MB.UseNormL1 = UseNormL1;    
    for (i=0; i < NbrBase; i++) Tab[i].alloc(Nx, "alloc");
    for (i = 0; i < NbrBase; i++) MB.TabSelect[i] = TabSelect[i];
    
    MB.COS_Overlapping = COS_Overlapping;
    MB.COSMin = COSMin;
    MB.COS_Sensibility= COS_Sensibility;
    MB.COS_SuppressDctCompCont = SuppressDctCompCont;
    if (CosBlockSize > 0) MB.COS_BlockSize = CosBlockSize;
    
    MB.MCOS_FirstDetectScale=FirstDetectScale;
    MB.MCOS_Overlapping = COS_Overlapping;  //For the same option -O than LDCT
    MB.MCOS_PositivRecSig = False;
    if (CosBlockSize > 0) MB.MCOS_FirstBlockSize = CosBlockSize;
    MB.MCOS_NbrScale1D = NbrScale1D;
    
    MB.ALDCT_NbrScale1D = NbrScale1D;
    MB.ALDCT_Sensibility = COS_Sensibility;
    MB.ALDCT_InfoCost = Infocost;
        
    if (Noise_Sig == -1) 
       Noise_Sig = mr1d_detect_noise (Data);
 
       
    MB.DataSigmaNoise = Noise_Sig;
    MB.CFAR = CFAR;
    MB.CFDR = CFDR;
    MB.CFAR_NsigmaMad = CFAR_NsigmaMad;
    MB.CFDR_Qparam = CFDR_Qparam;
    
    // MB.alloc (Nx, Tab); 
    MB.alloc (Nx);
    
    if (UseMask == True)
      {
      MB.MaskedData.alloc(Nx, "Mask");
        for (i=0; i < Nx; i++)
          {
          if (Data(i) == 0) MB.MaskedData(i) = 0;
	        else 
	          {
	          MB.MaskedData(i) = 1;
	          }
	        }
      }
 
    if (RemoveSmoothPlane)
    {
       if (NbRemoveSmoothPlane<0) NbRemoveSmoothPlane=MAX_SMOOTH_REMOVE;
       if (Verbose == True)
          cout << "NbRemoveSmoothPlane == " << NbRemoveSmoothPlane <<endl;
       MB.remove_smooth_plane (Data,NbRemoveSmoothPlane);  
    }
    
    if (FirstSoftThreshold < 0) 
         MB.FirstSoftThreshold=MB.compute_first_lambda (Data);
    else MB.FirstSoftThreshold = FirstSoftThreshold;
      
    MB.LastSoftThreshold = LastSoftThreshold;  
    
    MB.infoVerbose (nameSignalIn, nameSignalOut);

    MB.decomposition (Data);
    MB.reconstruction (Result);
    
    if (RemoveSmoothPlane) Result += MB.Smooth;

    io_1d_write_data(nameSignalOut, Result, &Header);

    if (WriteAllRec == True) MB.write_allima (nameSignalOut, Header);
    exit(0);
}






 
