
#include "MR_Sat.h"
#include "IM_Sigma.h"
#include "IM_Math.h"
#include <time.h>
 #include <sys/types.h>
  
char Name_Imag_In[256]; 
char Name_Imag_Out[256]; 
int Nbr_Plan=DEFAULT_NBR_SCALE; 

const int MAX_DIRECTION = 36;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
Bool Dbg = False;
type_border Bord = I_CONT;
sb_type_norm Norm = NORM_L1;
int DirectionNumber = 4;
DirManager::DirectionType DirectionType = DirManager::STANDARD;
Bool Overlap = False;
Bool CompSigma = False;
Bool SuppressPos = False;
Bool SuppressLastScale = False;
Bool SuppressIsolPixel = False;
Bool WriteSupport = False;
Bool WriteNumSim = False;
Bool ReconsInFourier = False;
Bool IncreaseDirNumber = False; 
Bool CteFalseDetection = False;
Bool ComputeNumSigma = False;
int ImposedSNRMaxLevel = -1;
float ImposedSigmaImage = -1;
Bool WriteNormTransf = False;
Bool UsedRandomnImag = False;
Bool DetectOnlyPositive = False;
Bool UseNSigma = False;
float N_Sigma = 0;

Bool Filter = False;

#define TR_DEBUG 0 
 
/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    //nbr_scale_usage(Nbr_Plan);
    //manline();
    //vm_usage();
    //manline();
    //verbose_usage();    
    //manline();
    
    
    manline();

    fprintf(OUTMAN, "         [-n ScaleNumber] \n");     
    fprintf(OUTMAN, "             default is 4 scales \n");     
    manline();
   
    fprintf(OUTMAN, "         [-d DirectionNumber] \n");     
    fprintf(OUTMAN, "             number default is 4 directions \n");     
    manline();
    
    
    fprintf(OUTMAN, "         [-e] \n"); 
    fprintf(OUTMAN, "             reconstruction is done is the Fourier space, \n");
    fprintf(OUTMAN, "             default is each sacle is reconstruct by sommation of directions \n");
    fprintf(OUTMAN, "             in direct space : $s = \\sigma_d ((FFT^(-1)(mask_d*FFT(s)))$,\n");
    manline();
    
    fprintf(OUTMAN, "         [-o] \n");     
    fprintf(OUTMAN, "             overlap in direction, default is no \n");     
    fprintf(OUTMAN, "             if set, it is used in the transformation by \n");     
    fprintf(OUTMAN, "             $d (scale s) = FFT^{-1 } {FFT(s) * mask_d_with_overlap} $ \n");     
    fprintf(OUTMAN, "             but not in the reconstruction (at least for the moment...) \n");     
    manline();
 
    fprintf(OUTMAN, "         [-O Overlap Type] \n");     
    fprintf(OUTMAN, "             overlap type use for the direction \n");     
    fprintf(OUTMAN, "             default is standard (0) \n");     
    manline();
    
    
    fprintf(OUTMAN, "         [-h] \n");     
    fprintf(OUTMAN, "             if set increase the number of direction with scale \n");     
    fprintf(OUTMAN, "             the direction number is multiply by 2 each 2 sacles  \n");     
    fprintf(OUTMAN, "             for scale s we have $d(s) = 2^{s/2}$  \n");     
    manline();
    
    fprintf(OUTMAN, "         [-f] \n");     
    fprintf(OUTMAN, "             filtering...\n");     
    manline();
    
    fprintf(OUTMAN, "         [-b ImposedLevelToThreshold] \n");     
    fprintf(OUTMAN, "             only with -f\n");   
    fprintf(OUTMAN, "             used to compare result with a standartd wavelet filtering if -b0  \n"); 
    manline();

    fprintf(OUTMAN, "         [-m] \n");     
    fprintf(OUTMAN, "             compute sigma with simulations \n");   
    manline();
    fprintf(OUTMAN, "         [-G] \n");     
    fprintf(OUTMAN, "             used randomn image in transformation \n");   
    manline();
    fprintf(OUTMAN, "         [-y] \n");     
    fprintf(OUTMAN, "             used randomn image in transformation \n");   
    manline();
    fprintf(OUTMAN, "             write sulation transf ans sigma for all direction and all scale \n");
    fprintf(OUTMAN, "             in file numsim_transf_s<num of scale>.fits \n"); 
    fprintf(OUTMAN, "             in file numsim_sigma_s<num of scale>.fits \n"); 
    manline();
     
    fprintf(OUTMAN, "         [-x] \n");    
    fprintf(OUTMAN, "             debug mode, trace and write on files all results \n"); 
    fprintf(OUTMAN, "             ==> Traces: \n"); 
    fprintf(OUTMAN, "             TRANSFORM \n"); 
    fprintf(OUTMAN, "             info on each scale, and on each FFT of scale \n"); 
    fprintf(OUTMAN, "             for each scale, and each direction \n"); 
    fprintf(OUTMAN, "                a)info on FFT(s)*mask_direction \n"); 
    fprintf(OUTMAN, "                b)info on FFT^{-1}(FFT(s)*mask_direction)\n");
    fprintf(OUTMAN, "                c)info on direction = real (FFT^{-1}(FFT(s)*mask_direction)\n");
    fprintf(OUTMAN, "             RECONS \n");  
    fprintf(OUTMAN, "             for each scale, and each direction \n");
    fprintf(OUTMAN, "                A)info on direction \n");
    fprintf(OUTMAN, "             info on each scale (sum of direction) \n");
    fprintf(OUTMAN, "             if [-e] \n");
    fprintf(OUTMAN, "                RECONS \n");     
    fprintf(OUTMAN, "                for each scale, and each direction \n");  
    fprintf(OUTMAN, "                   A)info on direction d  \n");
    fprintf(OUTMAN, "                   B)info on FFT(direction)  \n");
    fprintf(OUTMAN, "                   C)info on FFT(direction)*mask_direction   \n"); 
    fprintf(OUTMAN, "                   D)info on FFT(scale) \n"); 
    fprintf(OUTMAN, "                for each scale  \n"); 
    fprintf(OUTMAN, "                   D)info on FFT(scale) \n"); 
    fprintf(OUTMAN, "                   E)info on scale s (inv FFT) \n"); 
    fprintf(OUTMAN, "             if [-f] \n"); 
    fprintf(OUTMAN, "                see in the Files section  \n");      
    fprintf(OUTMAN, "             ==> Files: \n"); 
    fprintf(OUTMAN, "             \"dir_support<direction_number>_<direction_number_in_scale>.fits\" support of direction \n"); 
    fprintf(OUTMAN, "             \"mr_fft<scale_number>_re(or im).fits\" fft of scale_number \n");  
    fprintf(OUTMAN, "             if [-f] \n");    
    fprintf(OUTMAN, "                all the information is write on \"info.dat\" file \n");  
    manline();
   
    fprintf(OUTMAN, "         [-W] \n");     
    fprintf(OUTMAN, "             write the normalized transf for all direction and all scale \n");
    fprintf(OUTMAN, "             in file norm_transf_s<num of scale>.fits \n"); 
    manline();
   
   
    fprintf(OUTMAN, "         [-g] \n");     
    fprintf(OUTMAN, "             sigma = noise standard deviation \n");
    fprintf(OUTMAN, "             default is automatically estimated \n"); 
    manline();
   
   
    fprintf(OUTMAN, "         [-a] \n");     
    fprintf(OUTMAN, "             used to compute the sigma value for the couple  (scale number n, direction number d) \n");
    fprintf(OUTMAN, "             compute only the transformation, result are in the file \"transf_s<scale_number>.fits\", \n");
    fprintf(OUTMAN, "             dimension are (x size, y size, 2*direction_number-1)  \n"); 
    fprintf(OUTMAN, "             in image may be a uniform random noise with sigma=1 \n"); 
    manline();
  
    fprintf(OUTMAN, "         [-c] \n");     
    fprintf(OUTMAN, "             write support for each scale (for each direction 2^(2*d-1) dir)\n");
    fprintf(OUTMAN, "             in file \"support_s<scale_number>.fits\" \n"); 
    manline();
    
    fprintf(OUTMAN, "         [-v] \n");   
    fprintf(OUTMAN, "             verbose mode, trace what the program is doing,  \n"); 
    fprintf(OUTMAN, "             trace the support of the directions  \n");
    fprintf(OUTMAN, "             if [-o] \n");    
    fprintf(OUTMAN, "                the support of overlap is olso traced [overlap]U[support]U[overlap]\n");
    fprintf(OUTMAN, "             write the support of the horizontal direction \n");
    manline();
    
    fprintf(OUTMAN, "         default : \n");   
    fprintf(OUTMAN, "             compute only the transformation with 4 directions,\n");
    fprintf(OUTMAN, "             write the reconstruct image \n");
    
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
    while ((c = GetOpt(argc,argv,"g:n:d:b:LoO:s:faerhcmGylvPkKWxzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'L': Norm = NORM_L2; break;
	   case 'f': Filter = True; break;
	   case 'v': Verbose = True; break;
	   case 'x': Dbg = True; break;
	   case 'o': Overlap = True; break;
	   case 'P': SuppressPos = True; break;
           case 'K': SuppressLastScale = True; break;
           case 'k': SuppressIsolPixel = True; break;
           case 'a': CompSigma = True; break;
           case 'c': WriteSupport = True; break;
           case 'e': ReconsInFourier = True; break;
           case 'h': IncreaseDirNumber = True; break;
	   case 'r': CteFalseDetection = True; break;
	   case 'W': WriteNormTransf = True; break;
           case 'm': ComputeNumSigma = True; break;
           case 'G': UsedRandomnImag = True; break;
           case 'y': WriteNumSim = True; break;
	   case 'l': DetectOnlyPositive = True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		break;
	   case 'd':
		/* -d <Nbr_Direction> */
		if (sscanf(OptArg,"%d",&DirectionNumber) != 1) 
                {
		    fprintf(OUTMAN, "bad number of directions: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_DIRECTION)) 
                 {
		    fprintf(OUTMAN, "bad number of directions: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < number of directions <= %d\n", MAX_DIRECTION);
 		    exit(-1);
		}
		break;
 	   case 'b':
		/* -d <Imposed SNR Max level> */
		if (sscanf(OptArg,"%d",&ImposedSNRMaxLevel) != 1) 
                {
		    fprintf(OUTMAN, "bad imposed SNR max level: %s\n", OptArg);
		    exit(-1);
		}
		break;               
 	   case 'g':
		/* -s <Sigma imposed> */
		if (sscanf(OptArg,"%f",&ImposedSigmaImage) != 1) 
                {
		    fprintf(OUTMAN, "bad imposed sigma: %s\n", OptArg);
		    exit(-1);
		}
		break;               
 	   case 'O':
		/* -O <Overlap type> */
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "bad overlap type: %s\n", OptArg);
		    exit(-1);
		}
		if( c >= 0 && c <= 1 ) DirectionType = (DirManager::DirectionType)c;
		break;               
	    case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                UseNSigma = True;
                if (N_Sigma <= 0.)  N_Sigma = 1;
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
		usage(argv);
	}


#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}

/*********************************************************************/

int main(int argc, char *argv[])
{

   transinit(argc, argv);
   Ifloat Dat;
   io_read_ima_float(Name_Imag_In, Dat);

   float sigmaImag=0;
   if (ImposedSigmaImage < 0) {
      sigmaImag = detect_noise_from_med (Dat);
      if (sigmaImag == 0) sigmaImag=1;
   } else 
      sigmaImag = ImposedSigmaImage;

   if (Verbose == True) {
      cout << endl << "  PARAMETERS  " << endl;
      cout         << "--------------" << endl;
      cout << "File Name in = " << Name_Imag_In << endl;
      cout << "File Name Out = " << Name_Imag_Out << endl;
      if (Norm == NORM_L2) cout << "L2 normalization" << endl;
      if (Overlap) cout << "overlap : ON" << endl;
      else cout << "overlap : NO" << endl;
      if( DirectionType == DirManager::STANDARD ) cout << "overlap type : STANDARD" << endl;
      if( DirectionType == DirManager::NEW ) cout << "overlap type : NEW" << endl;
      if (SuppressPos) cout << "Suppress positivity constraint" << endl;
      cout << "Number of scales = " << Nbr_Plan << endl;
      if (SuppressLastScale) cout << "Suppress last scale" << endl;
      if (SuppressIsolPixel) cout << "Suppress isolated pixel" << endl;
      if (ImposedSNRMaxLevel >= 0) 
         cout << "Imposed level of max SNR is " << ImposedSNRMaxLevel << endl;
      if (IncreaseDirNumber) cout << "Increase direction number in scale" << endl;
      if (CteFalseDetection) cout << "Work with constant false detection" <<endl;
      if (ReconsInFourier) cout << "Recons in Fourier" << std::endl;
      cout << "Border type = " << Bord << endl;
      cout << "Number of direction = " << DirectionNumber << endl;
      cout << "sigma of image = " << sigmaImag << endl;
      if( DetectOnlyPositive ) cout << "Detect only positive structure" << endl;
      cout << endl;
   }
   
   MR_Sat sat;
   sat.set_verbose(Verbose);
   sat.set_debug(Dbg);
   sat.set_suppress_pos(SuppressPos);
   sat.set_suppress_last_scale(SuppressLastScale);   
   sat.set_suppress_isol_pixel(SuppressIsolPixel);
   sat.set_imposed_level_of_max_snr(ImposedSNRMaxLevel);
   sat.set_sigma_image(sigmaImag);
   sat.set_write_support(WriteSupport);
   sat.set_write_numsim(WriteNumSim);
   sat.set_recons_in_fourier(ReconsInFourier);
   sat.set_increase_direction_number(IncreaseDirNumber);
   sat.set_cte_false_detect(CteFalseDetection);
   sat.set_write_norm_transf(WriteNormTransf);
   sat.alloc( Dat.nl(), Dat.nc(), Nbr_Plan, DirectionNumber, Overlap, DirectionType );
   sat.set_border(Bord);
   sat.set_norm(Norm);
   
   Ifloat randNorm; 
   //unsigned int seed=1234567l;
   if (ComputeNumSigma) {
      if (Verbose) {
         cout << "  NUMERICAL SIMULATION FOR SIGMA ESTIMATION " << endl;
         cout << "---------------------------------------------" << endl;  
      }
      sat.set_comp_sigma(True);
      randNorm.alloc(Dat.nx(), Dat.ny());
      time_t tloc;
      int seed = time(&tloc);
      srand(seed);
      im_noise_gaussian (randNorm, 1.,seed);
 std::cout << "sigma=" << randNorm.sigma() << std::endl;
      /*for (long i=0; i<randNorm.n_elem(); i++) {
         if( i/100*100 == i ) 
	    std::cout << i << "/" << randNorm.n_elem() << std::endl;
         randNorm(i)/=randNorm.sigma(); 
      }*/
      sat.transform(randNorm);
      sat.threshold(0, False);
      sat.set_simulated_sigma();
      sat.set_comp_sigma(False);
   }
   sat.set_comp_sigma(CompSigma);
   sat.set_detect_only_positive( DetectOnlyPositive );


   // TRANSFORM  
   //==========
   if (ComputeNumSigma && UsedRandomnImag) {
      if (Verbose) {
         cout << endl;      
         cout << "  TRANSFORM  " << endl;
         cout << "-------------" << endl;  
      }
      sat.set_sigma_image(randNorm.sigma());
      sat.transform(randNorm);
   } else {
      if (Verbose) {
         cout << endl;      
         cout << "  TRANSFORM  " << endl;
         cout << "-------------" << endl;  
      }
      sat.transform(Dat);
   }
     
   
   // FILTER
   //=======
   if (Filter || UsedRandomnImag) {
      if (Verbose) {
         cout << endl;
         cout << "  THRESHOLD  " << endl;
         cout << "-------------" << endl;
      }
      if( ! UseNSigma ) N_Sigma = 3.;
      sat.threshold( N_Sigma );
   }

   
   // RECONS
   //=======
   if (!CompSigma) {
      if (Verbose) {
         cout << endl;
         cout << "  RECONS  " << endl;
         cout << "----------" << endl;
      }
      Ifloat Recons(Dat.nl(), Dat.nc());
      sat.recons(Recons);
      io_write_ima_float(Name_Imag_Out, Recons);
   }

   exit(0); 
} 
 
