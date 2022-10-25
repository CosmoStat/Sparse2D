/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/00
**    
**    File:  cur_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet transform program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Curvelet.h"
#include "IM_Prob.h"
#include "IM_Simu.h"
#include "IM_Sigma.h"
#include "FFTN_2D.h"
#include "IM_Deconv.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool Reverse=False; // Reverse transform
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
int NbrScale1D = -1;
int BlockSize= DEF_CUR_BLOCK_SIZE;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool BlockOverlap = True;
Bool DebugTrans = False;

Bool UsePSF=False;
char PSFFile[256];
Bool Deconv = False;
float N_Sigma = DEFAULT_N_SIGMA;
unsigned int InitRnd = 100;
float Eps= 1.e-3;
Bool NoiseIma = False;
float Noise_Ima=0.;
char NoiseFile[256];  /*  list of input image names  */

/***************************************/

/***************************************/
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is %d. \n", BlockSize);    
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             No block overlapping. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "        [-N NoiseFileName]\n");
    manline();
    gauss_usage();
    manline(); 
    nsigma_usage(N_Sigma);
    manline();
    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "             default is 100. \n\n");    
    manline();
    fprintf(OUTMAN, "         [-e Epsilon]\n");
    fprintf(OUTMAN, "             Frequencies where the TF or the |PSF|^2  is.\n");
    fprintf(OUTMAN, "             smaller than Eps are set to 0. \n"); 
    fprintf(OUTMAN, "             Default is %f. \n\n", Eps);     
    manline();
//     fprintf(OUTMAN, "         [-f]\n");
//     fprintf(OUTMAN, "             Filter each scan of the Radon transform.\n");  
//     manline();
// 
//     fprintf(OUTMAN, "         [-w]\n");
//     fprintf(OUTMAN, "             Filter width.\n");  
//     fprintf(OUTMAN, "             Default is 100.\n");  
//     manline();  
// 
//     fprintf(OUTMAN, "         [-s]\n");
//     fprintf(OUTMAN, "             Sigma paameter for the filtering.\n");  
//     fprintf(OUTMAN, "             Default is 10.\n");  
//     manline();  

    vm_usage();
    manline();    
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}
  
/*********************************************************************/

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
    while ((c = GetOpt(argc,argv,"P:e:s:I:g:XOb:n:N:vzZ")) != -1) 
    {
	switch (c) 
        { 
	      case 'e': /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&Eps) != 1) {
                fprintf(OUTMAN, "Error: bad Eps: %s\n", OptArg);
                exit(-1);
                }
                if (Eps <= 0.)  Eps = 0.001;
                break;
              case 's': /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) {
                fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
                exit(-1);
                }
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
                break;
       	     case 'I':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		    exit(-1);
		}
                InitRnd = (unsigned int) c;
		break;
            case 'g': /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) {
                    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
                   exit(-1);
                }
               break;
	   case 'P':
	        if (sscanf(OptArg,"%s",PSFFile) != 1) 
                {
                    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
                    exit(-1);
                }
		Deconv = UsePSF = True;
	        break;	
	   case 'N':
	        if (sscanf(OptArg,"%s",NoiseFile) != 1) 
                {
                    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
                    exit(-1);
                }
		NoiseIma = True;
	        break;	
		
 	   case 'X': DebugTrans = True; break;
           case 'O': BlockOverlap= (BlockOverlap==True) ? False: True;break;
           case 'v': Verbose = True;break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'n':
                /* -n <NbrScale2D> */
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
                    exit(-1);
                }
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
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
	
   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/****************************************************************************/

void dec_inverse1 (Ifloat &Imag, Icomplex_f &Psf_cf, Ifloat &Result, float Eps, Bool Apodization)
{
   FFTN_2D FFT2D;
   int i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nl1 = Psf_cf.nl();
   int Nc1 = Psf_cf.nc();
   Icomplex_f O_cf (Nl1, Nc1, "Buffer conv");
   // float FluxPsf = (FluxNorm == False) ? 1.: Psf_cf(Nl1/2,Nc1/2).real();
   float FluxPsf = 1. / Psf_cf(Nl1/2,Nc1/2).real();
   
   if ((Nl != Nl1) || (Nc != Nc1)) im_extend (Imag, O_cf);
   else 
   {
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++) O_cf(i,j) = complex_f(Imag(i,j),0);
   }
   if ((Apodization == True) && (Nl1 == 2*Nl) && (Nc1 == 2*Nc))
   {
     for (i = 0; i < Nl1; i++) 
     for (j = 0; j < Nc1; j++) 
     {
      float Coefi=1.;
      float Coefj=1.;
      float Nl4 = Nl1 / 4.;
      float Nl8 = Nl4 / 2.;
      float Nc4 = Nc1 / 4.;
      float Nc8 = Nc4 / 2.;
      if ((i >= Nl4) && (i < Nl1-Nl4))  Coefi=1.;
      else if ((i < Nl8) || (i > Nl1-Nl8))   Coefi = 0.;
      else if (i < Nl4) Coefi = ABS(sin( (Nl4-i-Nl8) / Nl8 * PI/2.));
      else if (i > Nl1-Nl4) Coefi = ABS(sin( (i-Nl4-Nl8) / Nl8 * PI/2.));
      
      if ((j >= Nc4) && (j < Nc1-Nc4))  Coefj=1.;
      else if ((j < Nc8) || (j > Nc1-Nc8))   Coefj = 0.;
      else if (j < Nc4) Coefj = ABS(sin( (Nc4-j-Nc8) / Nc8 * PI/2.));
      else if (j > Nc1-Nc4) Coefj = ABS(sin( (j-Nc4-Nc8) / Nc8 * PI/2.));
      Coefi=Coefj=1.;
      O_cf(i,j) = complex_f(O_cf(i,j).real()*Coefi*Coefj, 0.);
     }
   }
   io_write_ima_complex_f("xx_cf", O_cf);
   FFT2D.fftn2d(O_cf);
   // fft2d (O_cf);
   float Den;

   for (i = 0; i < Nl1; i++) 
   for (j = 0; j < Nc1; j++) 
   {
       Den = norm (Psf_cf(i,j));
       O_cf(i,j) *= conj(Psf_cf(i,j));
       if (Den > Eps) O_cf(i,j) /= Den;
       else O_cf(i,j) = complex_f(0.,0.);
   }
   FFT2D.fftn2d(O_cf, True);
   // fft2d (O_cf, -1);
   if ((Nl != Nl1) || (Nc != Nc1)) im_extract (O_cf, Result, FluxPsf);
   else for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++)  Result(i,j) = O_cf(i,j).real() / FluxPsf;
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
        if (Noise_Ima  >  FLOAT_EPSILON) 
	  cout << "Sigma Noise = " <<   Noise_Ima << endl;
	if (Deconv == True)
        {     
           cout << "Eps = " <<    Eps << endl; 
	   cout << "PSF File Name = " << PSFFile << endl;
        }
    }

    if (NbrScale2D < 1) 
    {
        cout << "Error: number of scales must be larger than 1 ... " << endl;
        exit(-1);
    }
    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    Curvelet Cur(SB1D);
    Cur.RidTrans = RidTrans;
    Cur.NbrScale2D = NbrScale2D;
    if (NbrScale1D >= 0)
    {
        Cur.NbrScale1D = NbrScale1D;
        Cur.GetAutoNbScale = False;
        if (NbrScale1D < 2) Cur.BlockOverlap = False;
        else Cur.BlockOverlap = BlockOverlap;
    }
    else Cur.BlockOverlap = BlockOverlap;
    Cur.Verbose = False;
 
    Ifloat Data, Result, NoiseData, GaussNoise;
    fltarray Trans;
    
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    int Nl = Data.nl();
    int Nc = Data.nc();
  
 //        if ((is_power_of_2(Data.nl()) == False) || (is_power_of_2(Data.nc()) == False))
//        {
//            cout << "Error: image size must be power of two ... " << endl;
//            exit(-1);
//        }
  if (NoiseIma == True) 
  {
      io_read_ima_float(NoiseFile, NoiseData);
      if ((NoiseData.nl() != Nl) || (NoiseData.nl() != Nc))
      {
         cout << "Error: the noise map must have the same size";
         cout << "       as the input data ... " << endl;
      }
      if (Noise_Ima  >  FLOAT_EPSILON) 
      {
         GaussNoise.alloc(Nl, Nc, "gnoise");
	 im_noise_gaussian (GaussNoise, Noise_Ima, InitRnd);
	 NoiseData += GaussNoise;
      }
  }
  else
  {
     if (Noise_Ima < FLOAT_EPSILON)  
     {
         Noise_Ima = detect_noise_from_med (Data);	  
         if (Verbose == True) cout << "Estimated Noise Std = " <<   Noise_Ima << endl;
     }
     NoiseData.alloc(Nl, Nc, "noise");
     im_noise_gaussian (NoiseData, Noise_Ima, InitRnd);
  }
  int i,j, NumUndec = -1;
 
  if (Deconv == True)
  {
     cout << "Deconvolution ... " << endl;
     FFTN_2D FFT2D;
     int Nl1 = 2*Nl;
     int Nc1 = 2*Nc;     
     Ifloat OutPsf(Nl1, Nc1, "OPsf");   
     Ifloat ImaPSF;
     Icomplex_f Psf_cf(Nl1, Nc1, "OPsfcf");
     io_read_ima_float(PSFFile, ImaPSF);
     cout << "dec_center_psf" << endl;
     dec_center_psf (ImaPSF, OutPsf);
     FFT2D.fftn2d(OutPsf, Psf_cf);
     cout << "dec_inverse1 ... " << endl;
     dec_inverse1 (NoiseData, Psf_cf, NoiseData, Eps, False);  
  io_write_ima_float("xx_invn.fits", NoiseData);
     dec_inverse1 (Data, Psf_cf, Data, Eps, False);
  }
  io_write_ima_float("xx_inv.fits", Data);
    cout << "Curvelet transform ... " << endl;
    if (BlockSize > Data.nc())
    {
        cout << "Error: Block size must lower than the image size ... " << endl;
        exit(-1);
    } 
    Cur.alloc(Data.nl(), Data.nc(), BlockSize);
    int NbrBand = Cur.nbr_band();
    fltarray TabMin(NbrBand);
    fltarray TabMax(NbrBand);
    CImaProb CP;   
    Cur.transform(NoiseData,Trans);
    fltarray Band;  
    cout << "Noise transform ... " << endl; 
    for (int b=0; b <NbrBand; b++) 
    {
       int s2d, s1d;
       int Nlb, Ncb;       
       double LMin, LMax;
       float SigmaBand;
       
       Ifloat IBand;
       Cur.get_scale_number(b, s2d, s1d);
       Cur.get_band(Trans, b, Band);
       Nlb= Band.ny();
       Ncb = Band.nx();
       IBand.alloc(Band.buffer(),Nlb, Ncb);
       // cout << "Band " << b+1 << " Size = " << IBand.nl() << " " << IBand.nc() << endl;
       CP.set(IBand);
       if (s1d == 0) CP.find_gthreshold(N_Sigma+1, LMin, LMax); 
       else CP.find_gthreshold(N_Sigma, LMin, LMax);
//        SigmaBand = Band.sigma();
//        if (s1d == 0)  SigmaBand *= (N_Sigma+1);
//        else SigmaBand *= N_Sigma;
//        LMin = -SigmaBand;
//        LMax = SigmaBand;
       TabMin(b) = LMin;
       TabMax(b) = LMax;
    }
    cout << "Data transform ... " << endl; 
    Cur.transform(Data,Trans);
    for (int b=0; b <NbrBand-1; b++) 
    {
       int s2d, s1d;
       int Nlb, Ncb;
       Ifloat IBand;
       Cur.get_scale_number(b, s2d, s1d);
       Cur.get_band(Trans, b, Band);
       Nlb= Band.ny();
       Ncb = Band.nx();
       IBand.alloc(Band.buffer(),Nlb, Ncb);
       cout << "Band " << b+1 << " Size = " << IBand.nl() << " " << IBand.nc();
       cout <<  " LMin = " <<  TabMin(b) << "LMax = " <<  TabMax(b) << endl;
       for (i=0;i < Nlb; i++)
       for (j=0;j < Ncb; j++) 
          if ((IBand(i,j) < TabMax(b)) && (IBand(i,j) > TabMin(b))) IBand(i,j) = 0.;
       Cur.put_band(Trans, b, Band);
    }    
    Result.alloc(Data.nl(), Data.nc(), "REC");
    Cur.recons(Trans,Result, False);
    threshold(Result);
    io_write_ima_float(Name_Imag_Out, Result, &Header);
    exit(0);
}
