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
**    File:  ridgelet.cc
**
*******************************************************************************
**
**    DESCRIPTION  ridgelet transform program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Ridgelet3D.h"
#include "IM3D_Block.h"

char Name_Dat_In[256];  // input file image  
char Name_Dat_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool Reverse=False; // Reverse transform
int Nbr_Plan = DEF_RID3D_NBR_SCALE;
int BlockSize=0;
Bool BlockOverlap = False;
type_ridgelet3d_WTtrans RidTrans = DEF_RID3D_TRANS;
Bool ExtractBand = False;
char FileSimu[180];
Bool WithFileSimu=False;
int NbrCoeffMax = 10;
Bool GetMax = False;
float N_Sigma=3;                // number of sigma (for the noise)  
Bool WeightedBlock = False;
Bool Normalize=False;       // normalize data in

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID3D_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRid3DTransform((type_ridgelet3d_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRid3DTransform(RidTrans));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             Default is automatically calculated. \n");    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size. Default is cube size. \n");    
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Block Overlapping. Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Block apodisation. Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Normalize the data. Default is no. \n");    
    manline();
       
    nsigma_usage(N_Sigma);
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
    while ((c = GetOpt(argc,argv,"Nws:m:t:Ob:n:S:vzZ")) != -1) 
    {
	switch (c) 
        {
          case 'N': Normalize = (Normalize == True) ? False: True; break;
	  case 'w': WeightedBlock = (WeightedBlock == True) ? False: True; break;
          case 's': 
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)  N_Sigma = 3;
		break;
          case 'm': 
 		if (sscanf(OptArg,"%d",&NbrCoeffMax ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of coeff: %s\n", OptArg);
	            exit(-1);
 		}
		if (c <= 0)
                {
		    fprintf(OUTMAN, "Error: bad number of coeff: %s\n", OptArg);
	            exit(-1);
 		}	              
                GetMax = True;
                break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID3D_TRANS)) RidTrans = (type_ridgelet3d_WTtrans) (c);
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
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                break;	   
           case 'S':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%s",FileSimu) != 1) 
                {
		    fprintf(OUTMAN, "bad simu file: %s\n", OptArg);
		    exit(-1);
		}
		WithFileSimu = True;
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
       if (OptInd < argc) strcpy(Name_Dat_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Dat_Out, argv[OptInd++]);
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

/***************************************************************************/
 
int main(int argc, char *argv[])
{
    int i,j,k,s;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    extern softinfo Soft;
    Soft.mr4();
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);
    
    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Dat_In << endl;
        cout << "File Name Out = " << Name_Dat_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRid3DTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
	if (WeightedBlock == True)  cout << "Block apodisation  " << endl;
    }

    fltarray Data;
    Ifloat Result;
    io_3d_read_data(Name_Dat_In, Data);
    int Nx = Data.nx();
    int Ny = Data.ny(); 
    int Nz = Data.nz();
    
    if (Verbose == True)
    { 
        printf( "Input Cube size = (%2d,%2d,%2d)\n", Nx,Ny,Nz);
    }  
    if (Normalize == True)
    {
       double Mean = Data.mean();
       double Sigma = Data.sigma();
       for (int i=0;i<Nx;i++)
       for (int j=0;j<Ny;j++)
       for (int k=0;k<Nz;k++) Data(i,j,k) = (Data(i,j,k)-Mean)/Sigma;
    }       
    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    Ridgelet3D RG(SB1D);
    RG.RidTrans = RidTrans;
    if (Nbr_Plan != DEF_RID3D_NBR_SCALE) 
    {
       RG.NbrScale = Nbr_Plan;
       RG.GetAutoNbScale = False;
//        if (Nbr_Plan < 2)
//        {
//           RG.BlockOverlap = False;
//        }
//        else RG.BlockOverlap = BlockOverlap;
    }         
//     else RG.BlockOverlap = BlockOverlap;
    RG.NbrScale = Nbr_Plan;
    RG.Verbose = Verbose;
    RG.set_block_param(BlockSize,BlockSize,BlockSize,BlockSize);

       
    if (BlockSize > Data.nx())
    {
        cout << "Warning: Block size must lower than the cube size ... " << endl;
        cout << "         Block size is set to cube size " << endl;
        BlockSize = 0;
    }
    if (BlockSize == 0)
    {
       BlockSize = MIN(Nx,Ny);
       BlockSize = MIN(BlockSize,Nz);
    }    

    Block3D B3D;
    B3D.BlockOverlap = BlockOverlap;
    B3D.Verbose = Verbose;
    B3D.alloc(Nx,Ny,Nz,BlockSize);
    int Nbx = B3D.nbr_block_nx();
    int Nby = B3D.nbr_block_ny();
    int Nbz = B3D.nbr_block_nz();
    if (Nx % BlockSize != 0) Nbx--;
    if (Ny % BlockSize != 0) Nby--;
    if (Nz % BlockSize != 0) Nbz--;
    
    fltarray Block;
    Block.alloc(BlockSize,BlockSize,BlockSize);
    fltarray StaInfo;
    fltarray Simu;
    int NbrStat = 8;
    int NbrBand = RG.NbrScale;
    if (WithFileSimu) 
    {
        NbrStat++;
	StaInfo.alloc(NbrBand,NbrStat);
	fits_read_fltarr(FileSimu, Simu);	
    }
    else StaInfo.alloc(NbrBand,NbrStat);
    fltarray MaxVal(NbrBand, NbrCoeffMax);
    dblarray x1(NbrBand);
    dblarray x2(NbrBand);
    dblarray x3(NbrBand);
    dblarray x4(NbrBand);
    dblarray Mean(NbrBand);
    dblarray Sigma(NbrBand); 
    dblarray Skew(NbrBand);
    dblarray Curt(NbrBand);
    fltarray Min(NbrBand);
    fltarray Max(NbrBand);
    intarray Np(NbrBand); 
    dblarray Energy(NbrBand);
    dblarray EnergyNSig(NbrBand);
    dblarray EnergyNSigSimu(NbrBand);    
            
    Ifloat Band;  
    for (i=0; i < Nbx; i++)
    for (j=0; j < Nby; j++)
    for (k=0; k < Nbz; k++)
    {    
       float *Dat;
       int N;
       // if (Verbose == True) printf("Block: %d, %d, %d\n", i,j,k);
       B3D.get_block_cube(i, j, k, Data,Block,WeightedBlock);
       // Block.info("BBlock");
       RG.transform(Block,Result,BlockSize);
       RG.Verbose = False;
       B3D.Verbose = False;
       for (s=0; s< NbrBand;s++)
       { 
          RG.get_scale(Result,Band,s);
          Dat = Band.buffer();
          N = Band.n_elem();
          Np(s) += N;
	  float M1 = Band.mean();
	  float M2 = Band.sigma();
          for (int p=0; p < N; p++)
          {
	     float Coef = Dat[p];
	     double E2 = Coef*Coef;
             x1(s) +=  Coef;
             x2(s) +=  E2;
             x3(s) +=  E2*Coef;
             x4(s) +=  E2*E2;
	     if (p == 0) Min(s) = Max(s) = Coef;
	     else
	     {
	        if (Min(s) > Coef) Min(s) = Coef;
	        if (Max(s) < Coef) Max(s) = Coef;
 	     }
	     if (fabs(Coef-M1) > N_Sigma*M2)
	                            EnergyNSig(s) += E2;
	     if (WithFileSimu == True)
	          if (fabs(Coef-M1)> N_Sigma*Simu(s))
	                            EnergyNSigSimu(s) += E2;
	  }
       }
    }
    
    
    for (s=0; s< NbrBand;s++)
    { 
       Energy(s) = x2(s);
       x1(s) /= (double) Np(s);
       x2(s) /= (double) Np(s);
       x3(s) /= (double) Np(s);
       x4(s) /= (double) Np(s);
       Sigma(s) = x2(s) - x1(s)*x1(s);
       if (Sigma(s) > 0.)
       { 
          double x1_2 = x1(s)*x1(s);
          double x1_4 = x1_2*x1_2;
          Sigma(s) = sqrt(Sigma(s));
          Skew(s) = 1. / pow(Sigma(s),(double) 3.) * (x3(s) - 3.*x1(s)*x2(s)+2.*x1(s)*x1(s)*x1(s));
          Curt(s) = 1. / pow(Sigma(s),(double) 4.) * (x4(s) -4*x1(s)*x3(s) + 6.*x2(s)*x1_2 -3.*x1_4 ) - 3.;
       }
       else Curt(s) = 0.;
       StaInfo(s,0) = (float) x1(s);
       StaInfo(s,1) = (float) Sigma(s);
       StaInfo(s,2) = (float) Skew(s);
       StaInfo(s,3) = (float) Curt(s);
       StaInfo(s,4) = Min(s);
       StaInfo(s,5) = Max(s);
       StaInfo(s,6) = Energy(s);
       StaInfo(s,7)= EnergyNSig(s);
       if (WithFileSimu == True)
              StaInfo(s,8)=EnergyNSigSimu(s);
      
	//!!!!!!!!!!!!!!! modify the Band values !!!!!!!!!!!!!!!!!!!!!!!
//         if (GetMax == True)
// 	{
// 	    for (int i=0; i <  NbrCoeffMax;i++) 
// 	    {
// 	        int IndMax=0;
// 	        MaxVal(s,i) = Band.maxfabs(IndMax);
// 	        Band(IndMax)=0.0;
// 	    }   
// 	}
        if (Verbose == True)
	{
            printf("  Band %d: Nx = %d, Ny = %d,  Min = %5.3f, Max = %5.3f\n",
	           s+1, Band.nx(), Band.ny(),  StaInfo(s, 4), StaInfo(s, 5));
            printf("           Mean = %5.3f, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Ener3sig = %5.3f\n",
	            StaInfo(s, 0), StaInfo(s,1),  StaInfo(s, 2),StaInfo(s,3), StaInfo(s,7));
	    if (WithFileSimu == True) printf("           Ener3sigNoise = %5.3f\n",  StaInfo(s,8));
        }
   }			 
   fits_write_fltarr (Name_Dat_Out, StaInfo);
//   if (GetMax == True) fits_write_fltarr ("MaxRidgelet.fits", MaxVal);		
   exit(0);
}

