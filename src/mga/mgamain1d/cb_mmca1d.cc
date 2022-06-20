/******************************************************************************
**                   Copyright (C) 2004 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck && Yassir Moudden
**
**    Date:  02/06/2003 
**    
**    File:  cb_mca1d.cc
**
*******************************************************************************
**
**    DESCRIPTION  Decomposition of a multichannel signal on multiple bases
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
#include "MatrixOper.h" 
#include <vector>

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

const int NBR_MAX_COMPONENTS = 100;
Bool TabBase[NBR_MAX_COMPONENTS];
type_mca1dbase TabSelect[NBR_MAX_COMPONENTS];

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

int NbrSource=0;

intarray TabBasisFirst(NBR_MCA1D);
intarray TabBasisLast(NBR_MCA1D);

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

    fprintf(OUTMAN, "         [-t TransformSelection[,FirstComponent[:LastComponent]]\n");
    for (int i = 0; i < NBR_MCA1D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMCA1DBase ((type_mca1dbase) i));
    fprintf(OUTMAN, "              default : %s and %s", 
                     StringMCA1DBase ((type_mca1dbase)0), 
                     StringMCA1DBase ((type_mca1dbase)2));
    manline();    
    
    fprintf(OUTMAN, "         [-c Number_of_Components]\n");
    fprintf(OUTMAN, "            Number of Components. Default is the number of channels.\n");
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
    
//     fprintf(OUTMAN, "         [-l]\n");
//     fprintf(OUTMAN, "             Supression of the smooth plane \n");
//     fprintf(OUTMAN, "             calculations. Default is no. \n");    
//     manline();
 
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
    fprintf(OUTMAN, "         [-u]\n");
    fprintf(OUTMAN, "             Number of undecimated scales.\n");
    fprintf(OUTMAN, "             default is all scales. \n");    
    manline();
    
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
            (argc,argv,"NAWLB:Ps:S:wu:t:KhF:Oi:g:n:vlR:zZkpc:I:")) != -1){
	switch (c){ 
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
	       {
	          int First, Last;
	          int Nr = sscanf(OptArg,"%d,%d:%d",&c, &First, &Last);
		  
                  if (Nr < 1)  
	          {
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", 
                    OptArg);
                    exit(-1);
                  }
                  if ((c > 0) && (c <= NBR_MCA1D))
		  {
                     if (TabBase[c-1] == True)
		     {
                       fprintf(OUTMAN, "Error: transform already selected: %s\n", 
                       OptArg);
                       exit(-1);
                      }
                      else
		      {
                      TabBase[c-1] = True;
                      TabSelect[NbrBase++] = (type_mca1dbase) (c-1);
                      }
                   }
                   else
		   {
                      fprintf(OUTMAN, "Error: bad type of transform: %s\n", 
                      OptArg);
                      exit(-1);
                   }
		   if ((Nr > 1) && (First > 0)) TabBasisFirst(c-1) = First;
		   if ((Nr > 2) && (Last > 0)) TabBasisLast(c-1) = Last;
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
           case 'c':
                if (sscanf(OptArg,"%d",&NbrSource) != 1){
                    fprintf(OUTMAN, "bad number of sources: %s\n", OptArg);
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

class MMCA1D  {
    MatOper MatrixOper;
   public:
     MCA1D *TabMCA;
     int NbrChannel;
     int NbrSource;
     dblarray MixingMatrix;       // Mixing matrix
     dblarray InvMixingMatrix;    // Inverse of the mixing Matrix
     int NpixData;
     intarray NbrBasePerSource;  // number of bases for each source
     type_mca1dbase **TabSelectPerSource; // selected bases for each component
     fltarray **TabComponentPerSources; // reconstructed signal for each bases
     fltarray TabSources; // reconstructed signal for each bases
     int NbrUserBases;
     type_mca1dbase TabBasisSelect[NBR_MCA1D]; // selected bases
     
     fltarray MChannelData;  // Multichannel data
     fltarray RecSources;    // Reconstructed sources
     float FirstSoftThreshold;
     float LastSoftThreshold;
     int Nbr_Iter;
     
     void init(fltarray &Data, int Nbr_Source, int NbrUserBases, type_mca1dbase * TabUserSelect, intarray & TabBasisFirst, intarray & TabBasisLast);
     void recsol_to_data(fltarray &Rec, dblarray &A, fltarray &Image);
     // calculate the estimated data from the mixing matrix and the reconstructed sources
     void get_residual(fltarray &Rec, dblarray &A, fltarray &Resi);
     // calculate the residual
     
     void recons_all_sources(fltarray &Rec);
     // Reconstruct the 2D reconstructed sources from their components
     void recons_one_source(int s, fltarray &Rec);
     // Reconstruct a 1D source
     void decomposition (fltarray& Signal);
     void update_mixing_matrix(float Lambda);
     void write_all_sources(char *Name, int Iter);
     void write_all_sources(char *Name,  fitsstruct & Header);
};

/***************************************************************************/

void MMCA1D::update_mixing_matrix(float Lambda)
{
    // Yassir New MixingMatrix calculation
     
    // Inverse matrix calulation. Works only when the number of lines >= number of columns
    MatrixOper.inv_mat_svd(MixingMatrix,InvMixingMatrix);
    MatrixOper.mat_print(InvMixingMatrix, "inv matrix");
    
    dblarray MatTest;
    MatrixOper.mat_mult(InvMixingMatrix, MixingMatrix, MatTest);
    MatrixOper.mat_print(MatTest, "Test matrix");
}

/***************************************************************************/

void MMCA1D::recons_one_source(int s, fltarray &Rec)
{
   int b,p;
   if (Rec.nx() != NpixData) Rec.alloc(NpixData);
   else Rec.init();
   for (b =0; b < NbrBasePerSource(s); b++)
   for (p=0; p<NpixData; p++) Rec(p) += (TabMCA[s].TabSigRec[b])(p);
}

/***************************************************************************/

void MMCA1D::recons_all_sources(fltarray &Rec)
{
   int s,b,p;
   if ((Rec.nx() != NpixData) || (Rec.ny() != NbrSource)) Rec.alloc(NpixData, NbrSource);
   else Rec.init();
   for (s =0; s < NbrSource; s++)
   {
       for (b =0; b <  NbrBasePerSource(s); s++)
       for (p=0; p<NpixData; p++) Rec(p,s) += (TabMCA[s].TabSigRec[b])(p);
   }
}

/***************************************************************************/

void MMCA1D::recsol_to_data(fltarray &Rec, dblarray &A, fltarray &Image)
{
   int p,i;
   dblarray Vectin(NbrSource);  
   dblarray Vectout(NbrChannel);
   for (p=0; p < NpixData; p++)
   {
      for (i=0; i < NbrSource; i++) Vectin(i) = Rec(p,i);
      MatrixOper.mat_mult(A,Vectin,Vectout);
      for (i=0; i < NbrChannel; i++)  Image(p,i) = Vectout(i);
   }
}

/***************************************************************************/

void MMCA1D::get_residual(fltarray &Rec, dblarray &A, fltarray &Resi)
{
   int p,i;
   fltarray Image(NpixData, NbrChannel);
   recsol_to_data(Rec,A,Image);
   for (p=0; p < NpixData; p++)
   for (i=0; i < NbrChannel; i++) Resi(p,i) = MChannelData(p,i) - Image(p,i);
}

/***************************************************************************/

void MMCA1D::write_all_sources(char *Name, int Iter)
{
    char NameIter[256];
    fltarray Vect(NpixData);
    for (int s=0;s < NbrSource; s++)
    {
       sprintf(NameIter, "%s_source%d_iter%d_", Name, s+1, Iter+1);
       for (int p=0;s < NpixData; p++) Vect(p) = RecSources(p,s);
       io_1d_write_data(NameIter,Vect);  
       
       sprintf(Name, "xx_source%d",s+1),
       TabMCA[s].write_allima(Name, Iter);
    }
}

/***************************************************************************/

void MMCA1D::write_all_sources(char *Name,  fitsstruct & Header)
{
    char NameIter[256];
    fltarray Vect(NpixData);
    for (int s=0;s < NbrSource; s++)
    {
       sprintf(NameIter, "%s_%d", Name, s+1);
       for (int p=0;s < NpixData; p++) Vect(p) = RecSources(p,s);
       io_1d_write_data(NameIter,Vect,&Header);
       TabMCA[s].write_allima(Name, Header);  
    }
}

/***************************************************************************/

void MMCA1D::init(fltarray &Data, int Nbr_Source, int N_UserBases, type_mca1dbase * TabUserBasisSelect, intarray & TabBasisFirst, intarray & TabBasisLast)
{
    int c,c1,s,s1,i,j;
    
    NbrChannel = Data.ny();
    NpixData =  Data.nx();
    cout << "   Number of channels = " << NbrChannel  << endl;
    cout << "   Number of samples per channel = " << NpixData << endl;
    NbrSource = Nbr_Source;
    if ((NbrSource >  NpixData) || (NbrSource <= 0)) NbrSource = NbrChannel;
    cout << "   Number of sources in the solution = " << NbrSource << endl;
    
    // initialize the matrix
    MixingMatrix.alloc(NbrSource, NbrChannel);
    for (i=0; i < NbrChannel; i++) 
    for (j=0; j < NbrSource; j++) MixingMatrix(i,j) = 1. / sqrt((float) NbrChannel);
    
    NbrBasePerSource.alloc(NbrSource);
    NbrUserBases = N_UserBases;
    for (c=0; c < NbrUserBases; c++) TabBasisSelect[c] = TabUserBasisSelect[c];
    
    cout << "Nbr Selected bases = " << NbrUserBases << endl;
    for (c=0; c < NbrUserBases; c++) 
      cout << StringMCA1DBase(TabBasisSelect[c]) << " from " << TabBasisFirst(c) << " to " << TabBasisLast(c) << endl;

    // Count the number of components per sources
    for (c = 0; c < NbrUserBases; c++)
    {
       if ((TabBasisFirst(c)  <= 0) || (TabBasisFirst(c) > NbrSource)) TabBasisFirst(c) = 1;
       if ((TabBasisLast(c) <= 0) || (TabBasisLast(c) > NbrSource)) TabBasisLast(c) = NbrSource;
       for (c1 = TabBasisFirst(c)-1; c1 < TabBasisLast(c); c1++) NbrBasePerSource(c1) ++;  
       cout << " basis: " << StringMCA1DBase(TabBasisSelect[c]) << ": First = " <<  TabBasisFirst(c) << " Last = " << TabBasisLast(c) << endl;
    }
    
    // Allocate TabSelectPerSource
    TabSelectPerSource = new type_mca1dbase * [NbrSource]; 
    for (s = 0; s <  Nbr_Source; s++)
    {
       TabSelectPerSource[s] = new type_mca1dbase[ NbrBasePerSource(s)];
       NbrBasePerSource(s) = 0;      
    }
    
    // Set TabSelectPerSource
    for (c = 0; c < NbrUserBases; c++)
    {
        type_mca1dbase BaseChoice = (type_mca1dbase) c;
        for (s = TabBasisFirst(c)-1; s < TabBasisLast(c); s++) 
        {
           TabSelectPerSource[s][NbrBasePerSource(s)] = BaseChoice;
	   NbrBasePerSource(s) ++;
	}
    }
 
    for (s = 0; s < Nbr_Source; s++)
    {
        cout << "Source " << s+1 << " NbrComponents = " << NbrBasePerSource(s) << endl;
	for (s1 = 0; s1 < NbrBasePerSource(s); s1++)
	  cout << "     Base " << s1+1 << ": " << StringMCA1DBase(TabSelectPerSource[s][s1]) << endl;
    }
    
    MChannelData = Data;
    RecSources.alloc(NpixData, NbrSource);
    
    TabMCA = new MCA1D[NbrSource];
    
    MatrixOper.inv_mat_svd(MixingMatrix,InvMixingMatrix);
    MatrixOper.mat_print(InvMixingMatrix, "inv matrix");
    
    dblarray MatTest;
    MatrixOper.mat_mult(InvMixingMatrix, MixingMatrix, MatTest);
    MatrixOper.mat_print(MatTest, "Test matrix");
}
 
/***************************************************************************/

void MMCA1D::decomposition (fltarray& Signal) 
{ 
   int i,b,s;
   // RegulIma RIM;
   MR_Regul RIM;
   RIM.NbrScale = 5;
   RIM.ExpDecreasingLambda = False;
   int Nx = Signal.nx();
   int Ny = Signal.ny();
   fltarray ResiOneComponent(Nx);
   fltarray ResiComponent(Nx,NbrSource);
   fltarray ResiData(Nx,Ny);

   float Lambda = FirstSoftThreshold;
   if (Nbr_Iter == 1) Lambda = LastSoftThreshold;
   float StepL = (FirstSoftThreshold - LastSoftThreshold ) / (float) (Nbr_Iter-1);

   float DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
   // if (Linear == False) StepL = DeltaThreshold / (float) (Nbr_Iter-1);

   // int NbrScaleTV = 5;
   // RIM.GradOperType = OPER_ENERGY;
   // RIM.GradOperType = OPER_LAPLACIAN;
   fltarray TVSig;
   
   if (Signal.nx() != NpixData)  
   {  
       cout << "Error: MCA class not initialized with this signal size" 
            << Signal.nx() << endl;
       cout << "       Expected size is " << NpixData << endl;
       exit(-1);
   }
   // YASSIR:  Resi = Signal;  // residual calculation
          
   for (int Iter=0; Iter < Nbr_Iter; Iter++) 
   {
       if (Iter>0) 
       {
          if (Linear == True) Lambda -= StepL; 
          else Lambda =  LastSoftThreshold 
                          + DeltaThreshold  *(1.-erf(2.8*Iter/ Nbr_Iter));
          if (Lambda < LastSoftThreshold) Lambda = LastSoftThreshold;
       }
        
       if (Verbose == True)
           cout << endl << "Iteration " << Iter+1 << " Lambda = " << Lambda << endl;
     
       for (s=0; s < NbrSource; s++)
       {
         // YASSIR: must estimated the residual in this source ==> ResiOneComponent
	 // Noise in TabMCA[s].DataSigmaNoise
         for (b=0; b < NbrBasePerSource(s); b++) 
         {
             TVSig = TabMCA[s].TabSigRec[b];
	  
             TabMCA[s].TabSigRec[b] += ResiOneComponent;
             switch(TabSelect[b]) 
	     {
 	     case MCA1D_ATROU:
                if (Verbose == True) cout << "Atrou proj " << endl;
		TabMCA[s].atrou_proj(TabMCA[s].TabSigRec[b], Lambda, Iter);   
		break;
              case MCA1D_WT:
                if (Verbose == True) cout << "WT proj " << endl;
                TabMCA[s].wt_proj(TabMCA[s].TabSigRec[b],  Lambda, Iter);                    
                break;
 	     case MCA1D_COS:
                if (Verbose == True) cout << "Local DCT proj " << endl;
                TabMCA[s].cos_proj (TabMCA[s].TabSigRec[b], Lambda, Iter);   
                break;
 	     case MCA1D_MCOS:                 
                if (Verbose == True) cout << "Local MDCT proj " << endl;
                TabMCA[s].mcos_proj(TabMCA[s].TabSigRec[b], Lambda, Iter);                    
                break;
                break;
 	     case MCA1D_ALDCT:
                if (Verbose == True) cout << "ALDCT proj " << endl;
		TabMCA[s].aldct_proj(TabMCA[s].TabSigRec[b], Lambda, Iter);   
		break;                
 	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }
          
  	 // Positivity constraint
	 switch(TabMCA[s].TabSelect[b]) 
	 {
             case MCA1D_ATROU: 
	        if (TabMCA[s].AT_PositivRecSig == True) (TabMCA[s].TabSigRec[b]).inf_threshold(0.0);
                break;
             case MCA1D_WT:
                if (TabMCA[s].WT_PositivRecSig == True) (TabMCA[s].TabSigRec[b]).inf_threshold(0.0);
                break;
	     case MCA1D_MCOS:
	        break;
	     case MCA1D_COS: 
 	        break;
	     case MCA1D_ALDCT:
	        break;                 
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }     
	  
	  // YASSIR: Update the residual
	  //  make_residual (Signal, Resi);
 	  
        } // ENDFOR (b=0 ...)
     } // ENDFOR (s=0 ...)
     
     // Yassir, here update the mixing matrix
     update_mixing_matrix(Lambda );
      
     // write all components all 10 iterations
     if (Write && Iter % 10 == 0) 
     {
          char Name[256];
	  sprintf(Name, "xx_source%d_",s+1),
          TabMCA[s].write_allima(Name, Iter);
	  sprintf(Name, "xx_source%d_",s+1),
	  write_all_sources("xx_rec", Iter);
      }
   } // ENDFOR (Iter=...)
} 

 
/***************************************************************************/

int main(int argc, char *argv[]) 
{
    int i,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    
    MMCA1D MB;  // Multiple base decomposition Class
    
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR3);
    
    for (int i=0; i<NBR_MCA1D; i++) TabBase[i] = False;

    filtinit(argc, argv);
    
    fltarray Data;
    fits_read_fltarr(nameSignalIn, Data, &Header);
    Header.origin = Cmd;
    int Nx = Data.nx();
    int Ny = Data.ny();
    
    if (Data.naxis() != 2)
    {
        cout << "Error: input data must have dimension 2 (multichannel 1D data) " << endl;
	exit(-1);
    } 
      
    if (NbrBase == 0) 
    {
       NbrBase = 2;
       TabSelect[0] = MCA1D_ATROU;
       if (AdaptDCT == False) TabSelect[1] = MCA1D_COS;
       else TabSelect[1] = MCA1D_ALDCT;
       
    }
    MB.init(Data, NbrSource, NbrBase, TabSelect, TabBasisFirst, TabBasisLast);
    if (NbrScale1D == 0) NbrScale1D=maxscalenumber(Nx);
    if (NbrScale1D > MAX_MGA_SCALE_1D) NbrScale1D = MAX_MGA_SCALE_1D;
    
    // Initiaze MCA for each component
    for (int s=0; s < MB.NbrSource; s++)
    {
        MB.TabMCA[s].NbrBase = MB.NbrBasePerSource(s);
        MB.TabMCA[s].Verbose = Verbose;
        MB.TabMCA[s].Nbr_Iter = Nbr_Iter;
        MB.TabMCA[s].PositivRecSig = PositivSol;
        MB.TabMCA[s].Nx=Nx;
        MB.TabMCA[s].Linear = Linear;
        MB.TabMCA[s].AT_PositivRecSig = PositivSol;
        MB.TabMCA[s].AT_NbrScale1D = NbrScale1D;
        MB.TabMCA[s].AT_FirstDetectScale=FirstDetectScale;
        MB.TabMCA[s].AT_KillLastScale=KillLastScale;
        MB.TabMCA[s].AT_OnlyPositivDetect = OnlyPositivDetect;
        MB.TabMCA[s].AT_SuppressIsolatedPixel = SuppressIsolatedPixel;
    
        MB.TabMCA[s].WT_NbrScale1D = NbrScale1D;
        MB.TabMCA[s].WT_PositivRecSig= PositivSol;
        MB.TabMCA[s].WT_NbrUndecimatedScale=NbrUndecimatedScale;
        MB.TabMCA[s].WT_OnlyPositivDetect = OnlyPositivDetect;
        MB.TabMCA[s].WT_SuppressIsolatedPixel = SuppressIsolatedPixel;
        // MB.WT_KillLastScale=KillLastScale;
    
        MB.TabMCA[s].UseNormL1 = UseNormL1;    
        for (i = 0; i < MB.NbrBasePerSource(s); i++) MB.TabMCA[s].TabSelect[i] = MB.TabSelectPerSource[s][i];
    
        MB.TabMCA[s].COS_Overlapping = COS_Overlapping;
        MB.TabMCA[s].COSMin = COSMin;
        MB.TabMCA[s].COS_Sensibility= COS_Sensibility;
        MB.TabMCA[s].COS_SuppressDctCompCont = SuppressDctCompCont;
        if (CosBlockSize > 0) MB.TabMCA[s].COS_BlockSize = CosBlockSize;
    
        MB.TabMCA[s].MCOS_FirstDetectScale=FirstDetectScale;
        MB.TabMCA[s].MCOS_Overlapping = COS_Overlapping;  //For the same option -O than LDCT
        MB.TabMCA[s].MCOS_PositivRecSig = False;
        if (CosBlockSize > 0) MB.TabMCA[s].MCOS_FirstBlockSize = CosBlockSize;
        MB.TabMCA[s].MCOS_NbrScale1D = NbrScale1D;
    
        MB.TabMCA[s].ALDCT_NbrScale1D = NbrScale1D;
        MB.TabMCA[s].ALDCT_Sensibility = COS_Sensibility;
        MB.TabMCA[s].ALDCT_InfoCost = Infocost;
        
        if (Noise_Sig == -1)  Noise_Sig = 1.; // mr1d_detect_noise (Data);
 
        MB.TabMCA[s].DataSigmaNoise = Noise_Sig;
        MB.TabMCA[s].alloc (Nx); 
    
//     if (RemoveSmoothPlane)
//     {
//        if (NbRemoveSmoothPlane<0) NbRemoveSmoothPlane=MAX_SMOOTH_REMOVE;
//        if (Verbose == True)
//           cout << "NbRemoveSmoothPlane == " << NbRemoveSmoothPlane <<endl;
//        MB.remove_smooth_plane (Data,NbRemoveSmoothPlane);  
//     }
    
//       if (FirstSoftThreshold < 0) 
//              MB.TabMCA[s].FirstSoftThreshold=MB.compute_first_lambda (Data);
//       else   MB.TabMCA[s].FirstSoftThreshold = FirstSoftThreshold;
       if (FirstSoftThreshold < 0) MB.TabMCA[s].FirstSoftThreshold=1.;
       MB.TabMCA[s].LastSoftThreshold = LastSoftThreshold;  
    }
    MB.FirstSoftThreshold=1.;
    MB.LastSoftThreshold = LastSoftThreshold;  
    MB.Nbr_Iter = Nbr_Iter;
    MB.decomposition (Data);
    
    // MB.infoVerbose (nameSignalIn, nameSignalOut);
     // if (RemoveSmoothPlane) Result += MB.Smooth;
    fltarray Result;
    MB.recons_all_sources (Result);
    fits_write_fltarr(nameSignalOut, Result, &Header);
    if (WriteAllRec == True) MB.write_all_sources(nameSignalOut,  Header);
    exit(0);
}






 
