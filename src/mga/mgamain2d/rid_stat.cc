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
#include "FFTN_2D.h"
#include "IM_Math.h"
#include "MR_Contrast.h"
#include "Curvelet.h"
#include "IM_Math.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
char Name_Dir_Out[256]; //  output file name  

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool Reverse=False; // Reverse transform
int NbrScale1D = -1;
int BlockSize= DEF_CUR_BLOCK_SIZE;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool BlockOverlap = False;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
Bool Clip = True;
Bool Sature = False;
Bool AngleStat = False;
float AngleVal = 0;
Bool DirectionAnalysis = False;
Bool BorderSize=False;
Bool WriteStat = False;
Bool MadAngleNorm= False;
Bool AllStat = False;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image [out_file]\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID_DECIMATED; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
    manline();
    
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is %d. \n", BlockSize);    
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Use overlapping block. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-A Angle]\n");
    fprintf(OUTMAN, "             Statistic for a given angle . Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-D DirectFileName]\n");
    fprintf(OUTMAN, "             Directional Analysis File Name. \n");    
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
static void filtinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"SMBD:A:t:Ob:n:vzZ")) != -1) 
    {
	switch (c) 
        {       
	   case 'S': AllStat = True; break;
	   case 'M': MadAngleNorm= True; break;
	   case 'B': BorderSize= True; break;
           case 'D': 	        
	        if (sscanf(OptArg,"%s",Name_Dir_Out) < 1)
		{
		   fprintf(OUTMAN, "Error: file name ... \n");
		   exit(-1);
		}
		DirectionAnalysis = True; 
		break;    
           case 'A': 
	        if (sscanf(OptArg,"%f",&AngleVal ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Angle: %s\n", OptArg);
                    exit(-1);
                }
	       AngleStat = True;
	       break;
           case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID_DECIMATED)) RidTrans = (type_ridgelet_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
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
                 if (sscanf(OptArg,"%d",&NbrScale1D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
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
	 
	if (OptInd < argc)
	{
	   strcpy(Name_Imag_Out, argv[OptInd++]);
	   WriteStat = True;
	}

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
 
/*************************************************************************/
    
int main(int argc, char *argv[])
{
    int b,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
 	if (NbrScale1D >= 0)
	 cout << "Nbr Scale 1D = " <<   NbrScale1D << endl;  
        if (AngleStat == True) cout << "Angle = " << AngleVal  << endl;
    }
    Ifloat Data;
    Ifloat Result;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    AngleVal = AngleVal / 180. * PI;
    
    if (RidTrans ==  RID_FINITE) 
    {    
        CPRIME_NUMBER CPN;

        if ((BlockSize > 0) && ( CPN.is_prime_number((unsigned long) BlockSize) == False))
        {
           cout << "Error: Block size must be a prime number ... " << endl;
           cout << "       Previous prime number = " << CPN.previous_prime_number((unsigned long) BlockSize) << endl;
           cout << "       Next prime number = " << CPN.next_prime_number((unsigned long) BlockSize) << endl;
           exit(-1);
        }
    }

    if (BlockSize > Data.nc())
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }
 
    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    Ridgelet RG(SB1D);
    RG.RidTrans = RidTrans;
    if (NbrScale1D  != DEF_RID_NBR_SCALE) 
    {
       RG.NbrScale = NbrScale1D;
       RG.GetAutoNbScale = False;
       if (NbrScale1D < 2)
       {
          RG.BlockOverlap = False;
       }
       else RG.BlockOverlap = BlockOverlap;
    }         
    else RG.BlockOverlap = BlockOverlap;
    RG.NbrScale = NbrScale1D;
    RG.Verbose = Verbose;
    RG.transform(Data,Result,BlockSize);
     
    // Calculate the statistic at each scale 
    int NbrBand = RG.NbrScale;
    int NbrStatPerBand = 5; // moment of order 2,3,4  
    if (AllStat == True) NbrStatPerBand = 9;
    
    fltarray TabStat(NbrBand-1, NbrStatPerBand);
    if (Verbose == True) cout << "NbrBand = " <<  NbrBand << endl; 
    
    Ifloat BandDir(Result.nl(), Result.nc(), "scale");
    
    fltarray DirectionStat;
    intarray TabNbrDir(NbrBand);
    if (DirectionAnalysis == True)
    {
       for (b=0; b < NbrBand-1; b++) 
       {        
 	  TabNbrDir(b) = RG.rid_one_block_nl();
       }
       DirectionStat.alloc(NbrBand-1, TabNbrDir.max(), NbrStatPerBand+1);
    }
    
    for (b=0; b < NbrBand-1; b++) 
    {
        int N,Nx,Ny;
	double Mean, Sigma, Skew, Curt;
	float  Min, Max;
	float HC1, HC2, TC1, TC2;
	if (MadAngleNorm == True) RG.mad_normalize_per_block(Result);
        // RG.mad_normalize_per_angle(Result);
	
        if (AngleStat == False) 
	{
	   Ifloat Band; 
 	   if (BorderSize == False) RG.get_scale(Result,Band,b);
	   else RG.get_scale_without_bord(Result, b, Band);
 
           N = Band.n_elem();
	   Nx = Band.nc();
	   Ny = Band.nl();
	   moment4(Band.buffer(), N, Mean, Sigma,  Skew, Curt, Min, Max);
	   if (AllStat == True) 
	   {
	      float *Buff = Band.buffer();
	      if (b != NbrBand-1)  hc_test(Buff, N,  HC1,  HC2, 0);
	      else hc_test(Buff, N,  HC1,  HC2);
              gausstest(Buff, N, TC1, TC2);
	   }
	}
	else 
	{
	   RG.get_scale_direction(Result, BandDir, b, AngleVal);
	   N = BandDir.nl()*BandDir.nc();
	   Nx = BandDir.nx();
	   Ny = BandDir.ny();
           float *Ptr = BandDir.buffer();
	   moment4(Ptr, N, Mean, Sigma, Skew, Curt, Min, Max);
	   if (AllStat == True) 
	   {
 	      hc_test(Ptr, N, HC1, HC2, 0.);
	      gausstest(Ptr, N, TC1, TC2);
	   }
        }
	
        TabStat(b, 0) = (float) Sigma;
	TabStat(b, 1) = (float) Skew;
	TabStat(b, 2) = (float) Curt;
        TabStat(b, 3) = (float) Min;
        TabStat(b, 4) = (float) Max; 
	if (Verbose == True)
            printf("  Band %d: Nl = %d, Nc = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Min = %5.3f, Max = %5.3f\n",
	           b+1, Ny, Nx, TabStat(b,0), 
		   TabStat(b, 1),TabStat(b,2),TabStat(b, 3), TabStat(b, 4));

        if (AllStat == True) 
	{
           TabStat(b, 5) = (float) HC1;
           TabStat(b, 6) = (float) HC2;
	   TabStat(b, 7) = (float) TC1;
           TabStat(b, 8) = (float) TC2;
	   if (Verbose == True)
            printf("   HC1 = %5.3f, HC2 = %5.3f, TC1 = %5.3f, TC2 = %5.3f\n",
 		   TabStat(b, 5), TabStat(b, 6), TabStat(b, 7), TabStat(b, 8));
	}
	
     	if (DirectionAnalysis == True)
	{	   
 	   // Ifloat  IBand;
 	   // RG.get_imablocksigma(Result, s1d, IBand);
	   int Ni = RG.rid_one_block_nl();
 	   for (int i = 0; i < Ni; i++)
	   {
	      fltarray VectPix;
	      DirectionStat(b, i, 5) = RG.get_angle_direction(i, True);
  	      RG.get_anglevect(Result, b, i, VectPix);
	      float *Ptr = VectPix.buffer();
 	      moment4(Ptr, VectPix.nx(), Mean, Sigma, Skew, Curt, Min, Max);
	      hc_test(Ptr, VectPix.nx(), HC1, HC2, 0.);
	      if (i == 0)  
	          printf("  0: N  = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Min = %5.3f, Max = %5.3f, HC1 = %5.3f, HC2 = %5.3f\n",
	             VectPix.nx(), (float) Sigma, (float) Skew,(float) Curt,(float) Min, (float) Max, (float) HC1, (float) HC2);
	      DirectionStat(b, i, 0) = (float) Sigma;
	      DirectionStat(b, i, 1) = (float) Skew;
	      DirectionStat(b, i, 2) = (float) Curt;
              DirectionStat(b, i, 3) = (float) Min;
              DirectionStat(b, i, 4) = (float) Max; 
	      if (AllStat == True) 
              {
		 hc_test(VectPix.buffer(), VectPix.nx(), HC1, HC2, 0.);
		 gausstest(VectPix.buffer(), VectPix.nx(), TC1, TC2);
 	         DirectionStat(b, i, 5) = (float) HC1;
                 DirectionStat(b, i, 6) = (float) HC2;
		 DirectionStat(b, i, 7) = (float) TC1;
                 DirectionStat(b, i, 8) = (float) TC2;
	      }
	   }
	}
    }
 
    if (WriteStat == True) fits_write_fltarr(Name_Imag_Out, TabStat); 
    if (DirectionAnalysis == True)  
        fits_write_fltarr(Name_Dir_Out, DirectionStat);
    exit(0);
}
