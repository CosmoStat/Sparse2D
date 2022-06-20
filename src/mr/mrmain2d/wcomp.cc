/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  16/03/2000 
**    
**    File:  wcomp.cc
**
*******************************************************************************
**
**    DECRIPTION  compression program
**    ---------- 
**     
******************************************************************************/


#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"
 
char File_Name_Imag[256];  /* input file image */
char File_Name_Transform[256]; /* output file name */
int Nbr_Plan=DEFAULT_CMP_NBR_SCALE;  /* number of scales */
float N_Sigma=1;   /* number of sigma (for the noise) */
int KeepResi = 0;
int KeepFitsHeader = 0;
Bool Verbose=False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

#define WRITE_PARAM 0

type_comp Comp_Method = COMP_MALLAT; // COMP_MRMEDIAN;
  
int BlockSize=0;
Bool UseBlock=False;
Bool Eval = False;
float SignalQuantif = DEFAULT_SIGNAL_QUANTIF;
Bool UnifQuant=False;

/*********************************************************************/

static void usage(char *argv[])
{
    
    fprintf(OUTMAN, "Usage: %s options in_image [out_file]\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-g QuantifParam]\n");
    fprintf(OUTMAN, "              Quantification parameter.\n");
    manline();    
             
    nbr_scale_usage(Nbr_Plan);
    manline();
    
    fprintf(OUTMAN, "         [-C BlockSize]\n");
    fprintf(OUTMAN, "              Compress by block. \n");
    fprintf(OUTMAN, "              BlockSize = size of each block.\n");
    fprintf(OUTMAN, "              Default is No.\n");
    manline();  

    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "              Evaluation mode: Decompression and PSNR estimation.\n");
    manline();   
    
    manline();    
    vm_usage();
    manline();           
    verbose_usage();
           
    manline();
    manline();
    manline(); 
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void hcinit(int argc, char *argv[])
{
    int c, L;
    char *Ptr; 
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif

    /* get options */
//    while ((c = GetOpt(argc,argv,"pm:vs:n:rc:g:R:NB:zZ:")) != -1) 
     while ((c = GetOpt(argc,argv,"us:EC:vn:g:zZ:")) != -1) 
   {
	switch (c) 
        {
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)  N_Sigma = 2;
 		break;
           case 'E': Eval = (Eval==True) ? False: True;break;
           case 'u': UnifQuant = (UnifQuant==True) ? False: True;break;
           case 'C': 
              if (sscanf(OptArg,"%d",&BlockSize ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
	            exit(-1);
                    
		}
		UseBlock = True;
                break;
	    case 'v':
		/* Verbose flag -v */
		Verbose = True;
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
	   case 'g':
		/* -c <SignalQuantif> */
		if (sscanf(OptArg,"%f",&SignalQuantif) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad SignalQuantif: %s\n", OptArg);
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
	if (OptInd < argc) 
             strcpy(File_Name_Imag, argv[OptInd++]);
         else usage(argv);
 
        strcpy (File_Name_Transform, File_Name_Imag);
 
	if (OptInd < argc) 
             strcpy(File_Name_Transform, argv[OptInd++]);

	L = strlen (File_Name_Transform);
	Ptr = File_Name_Transform;
	if (L < 5) strcat (File_Name_Transform, ".MRC");
	else if ((Ptr[L-1] != 'C') || (Ptr[L-2] != 'R') || (Ptr[L-3] != 'M')) 
               strcat (File_Name_Transform, ".MRC");
  
        if ((UseBlock == True) && (Eval == True))
	{
	   cout << "Error: -E and -C options cannot be used together ... " << endl;
	   exit(-1);
	}
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

/***************************************************/
/***************/

void print_cmp(Ifloat &ref, Ifloat &Imag)
{
    int i,j;
    int Nx = ref.nc();
    int Ny = ref.nl();
    long Nelem = Nx*Ny;
    double Snr,S1,S2;
    double Error;
    double Sum_X2= 0.,Sum_Y2= 0.,Sum_XY= 0.;
    double Moy = 0., Ecart = 0.;
    double  Moy_Err = 0., Ecart_Err = 0.;
    double Moy_Abs_Err = 0., Ecart_Abs_Err = 0.;
    double Snr_Db, Correl;
 
    for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
    {
        /* Calcul Correlation */
        Sum_X2 += ref (i,j) * ref (i,j);
        Sum_Y2 += Imag (i,j) * Imag (i,j);
        Sum_XY += ref (i,j) * Imag (i,j);

        /* Calcul of the standard deviation of ref */
        Moy += ref (i,j);
        Ecart += ref (i,j) * ref (i,j);

        /* Calcul of the standard deviation of the error */
        Error = ref(i,j)-Imag(i,j);
        Moy_Err += Error;
        Ecart_Err += Error * Error;

        Error = ABS(Error);
        Moy_Abs_Err += Error;
        Ecart_Abs_Err += Error*Error;
     }

    /* Correlation */
    Correl = Sum_XY / sqrt (Sum_X2 * Sum_Y2);

    /* variance of ref */
    Moy /= (double) Nelem;
    Ecart /= (double)  Nelem;
    S1 = Ecart - Moy*Moy;

    /* variance of the error */
    Moy_Err /= (double)  Nelem;
    Ecart_Err /= (double)  Nelem;
    S2 = Ecart_Err - Moy_Err*Moy_Err;

    /* variance of the absolute error */
    Moy_Abs_Err /=  (double) Nelem;
    Ecart_Abs_Err /=  (double) Nelem;
    Ecart_Abs_Err = Ecart_Abs_Err - Moy_Abs_Err*Moy_Abs_Err;

    /* signal to noise ration */
    Snr = S1 / S2;

    /* signal to noise ration (dB) */
    Snr_Db = 10. * log10 (Snr);

    ref -= Imag;
    float PSNR = sigma(ref);
    PSNR = 20. * log10 (255. / PSNR);

   printf ("     Sigma ( ABS(error) ) =  %f\n",  sqrt(Ecart_Abs_Err));
   printf ("     RMS =  %f\n", sqrt(S2));
   printf ("     SNR = variance(Ref) / variance(Error) = %f\n", Snr);
   printf ("     SNRdb = 10 log10(SNR) = %f dB\n", Snr_Db);
   printf ("     PSNR = 10 log10(SNR) = %f dB\n", PSNR);
}

/****************************************************************************/
int main(int argc, char *argv[])
{
    int k;
    char Cmd[256];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, 
       open input file(s) if necessary */
    lm_check(LIC_MR1);
    hcinit(argc, argv);
    MR_CompData CompIma;
    type_transform Transform  = TO_MALLAT; 
    if (UseBlock == True)
    {
         CompIma.UseBlock = True;
         CompIma.BlockSize = BlockSize;
    }
    CompIma.UseLiftingInsteadOfFilter = True;
    CompIma.UnifQuant = UnifQuant;
    CompIma.NoiseInData = False;
    CompIma.Cmd=Cmd;
    CompIma.File_Name_Imag=File_Name_Imag;
    CompIma.File_Name_Transform=File_Name_Transform;
    CompIma.Nbr_Plan=Nbr_Plan;
    CompIma.N_Sigma= N_Sigma;
    CompIma.Transform = Transform;
    CompIma.SignalQuantif = SignalQuantif;
    CompIma.KeepResi = KeepResi;
    CompIma.Comp_Method = Comp_Method;
    if (Verbose == True) CompIma.Verbose = True;
    CompIma.compress();


   if (Eval == True)
   {
       cout << endl << endl << "Decompression ... " << endl;
       MR_CompData DecIma;
       DecIma.Cmd=Cmd;
       DecIma.File_Name_Imag="xx_dec";
       DecIma.File_Name_Transform=File_Name_Transform;
       DecIma.Resol=0;
       DecIma.SimuNoise=False;
       DecIma.UseLiftingInsteadOfFilter = True;
       // if (Verbose==True) DecIma.Verbose = True;

      // apply the decompression
      DecIma.decompress();
    
      // save the result (field Dat) in DecIma.File_Name_Imag
      // in case where the full image is decompressed with block compression
      // options, the result is already in the output file.
      DecIma.write();
      cout << "Evaluation ... " << endl;
      Ifloat I1,I2;
      io_read_ima_float(File_Name_Imag, I1);
      io_read_ima_float(DecIma.File_Name_Imag, I2);
      print_cmp(I1, I2);
    }
    exit(0);
}


/***************************************************/

 

