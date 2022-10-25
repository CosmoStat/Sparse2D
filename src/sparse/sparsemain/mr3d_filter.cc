/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  04/09/99
**    
**    File:  mr3d_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution filtering of cube
**    ----------- 
**                 
**    Usage: mr3d_filter options cube_in cube_out
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
 
/****************************************************************************/
#define MAX_SCALE3D 10
static double TabNorm3DPaveBspline[MAX_SCALE3D] = 
   {0.956, 0.12, 0.035, 0.012, 0.004, 0.001, 0.0005, 0.,0., 0.};
  
char Name_Cube_In[256];
char Name_Cube_Out[256];

int Nbr_Plan=4;
type_trans_3d Transform=TO3_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

sb_type_norm Norm = NORM_L2;
type_sb_filter SB_Filter = F_MALLAT_7_9;
Bool Verbose=False;
int Max_Iter = 0;               // Maximum number of iteration  
float N_Sigma=3;                // number of sigma (for the noise)  
float Noise_Ima=0.;             // noise standard deviation  
Bool MAD = False;

const int NBR_OK_TRANS3D_METHOD = 2;
static  type_trans_3d TabMethod[NBR_OK_TRANS3D_METHOD] = {TO3_MALLAT,TO3_ATROUS};

/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/

static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options in_cube out_cube\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_OK_TRANS3D_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf3D((type_trans_3d)TabMethod[i]));
    fprintf(OUTMAN, "              default is %s\n", StringTransf3D((type_trans_3d) Transform));
					    
 
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();
    
    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n\n", StringSBFilter ((type_sb_filter) SB_Filter));
    
    manline();
    gauss_usage();
    manline();
    nsigma_usage(N_Sigma);
    manline();
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "             Correlated noise.\n");
    manline();
    vm_usage();
    manline();
    //max_iter_usage(Max_Iter);
    //manline();
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
    while ((c = GetOpt(argc,argv,"t:Cs:i:g:n:T:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 't': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad transform : %s\n", OptArg);
                    exit(-1);
                    
                }        
                if ((c <= 0) || (c > NBR_OK_TRANS3D_METHOD))
                {
                   fprintf(OUTMAN, "Error: bad transform type : %s\n", OptArg);
                   exit(-1);
                }
		 Transform  = (type_trans_3d) (TabMethod[c-1]);           
	      break;
	   case 'v': Verbose = True; break;
           case 'C': MAD = True; break;
	   case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
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
                if (N_Sigma <= 0.)  N_Sigma = 3;
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0)   Max_Iter = 0;
 		break;		
	    case 'T':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad filter: %s\n", OptArg);
	            exit(-1);
                    
		}
		SB_Filter = get_filter_bank(OptArg);
 		// SB_Filter = (type_sb_filter) (c);
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
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
 	   case '?':
			usage(argv);
		}
	}
   
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Cube_Out, argv[OptInd++]);
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

// static void mr_regul_cube_rec(MR_3D &MR_Step, fltarray &Result, int MaxIter, float LevelScale, float *TabLevel)
// {
//    int Nx = Result.nx();
//    int Ny = Result.ny();
//    int Nz = Result.nz();
//    int Nw = Nx*Ny*Nz;
//    int b,i,j,k,w,Iter=0;
//    MR_3D MR_Iter (Nx, Ny, Nz, MR_Step.Type_Transform, MR_Step.nbr_scale());
//    fltarray Lap(Nx, Ny, Nz);
//    MR_Iter.SB_Filter = MR_Step.SB_Filter;
//    MR_Iter.Norm =  MR_Step.Norm;
//    MR_Iter.LiftingTrans = MR_Step.LiftingTrans;
//    // MR_Iter.NormHaar = MR_Step.NormHaar;
//    
//    unsigned char *TabSup = new unsigned char [Nw];
//    float F;
//    float F_old=0; 
//    float Delta_F=10;;
//    float Level = LevelScale;
//    
//    // type_border Bord=I_CONT;
//    float *TabSave = new float [Nw];
//    float Delta_w,LevelQ;
//      
//    w =0;
//    for (b=0; b < MR_Step.nbr_band()-1; b++)
//    for (i=0; i < MR_Step.size_band_nx(b)-1; i++)
//    for (j=0; j < MR_Step.size_band_ny(b)-1; j++)
//    for (k=0; k < MR_Step.size_band_nz(b)-1; k++)
//    {
//       // Save initial value of dequantized coefficients for range constraint
//       TabSave[w]=MR_Step(b,i,j,k);
//       
//       // compute the multiresolution support
//       TabSup[w] = (ABS(MR_Step(b,i,j,k)) > FLOAT_EPSILON) ? 1 : 0; 
//       w++;
//    }
// 	  
//    while (Iter<MaxIter)
//    // while (((Delta_F*Delta_F) > .001 )&&(Iter<MaxIter))
//    {
//        MR_Step.recons(Result);
//        
//        // compute the Laplacian
//        for (i=1;i<Nx-1;i++)
//        for (j=1;j<Ny-1;j++)
//        for (k=1;k<Nz-1;k++)
//        
// 	   Lap(i,j,k) =  6*Result(i,j,k) - (Result(i-1,j,k) +
// 					Result(i,j-1, k) + 
// 					Result(i,j+1, k) +
// 					Result(i+1,j, k) +
// 					Result(i,j, k+1) +
// 					Result(i,j, k-1)) ;
// 
//        // positivity constraint on reconstructed image
//        // threshold(Result);
//        
//        // compute the multiresolution transform of the laplacian 
//        MR_Iter.transform(Lap); 
//        F=0.;   
//        // compute the  new solution
//        w=0;
//        for (b=0; b < MR_Step.nbr_band()-1; b++)
//        for (i=0; i < MR_Step.size_band_nx(b)-1; i++)
//        for (j=0; j < MR_Step.size_band_ny(b)-1; j++)
//        for (k=0; k < MR_Step.size_band_nz(b)-1; k++)       
//        {
//           if (TabLevel != NULL) Level = TabLevel[w];
// 	  
// 	  //main constraint
// 	  MR_Step(b,i,j,k) -= 0.2 * MR_Iter(b,i,j,k);
// 	  F += MR_Iter(b,i,j,k);
// 
// 	  // range constraint on thresholded-restored coefficients 
// 	  if ((ABS(MR_Step(b,i,j,k))>Level) && (TabSup[w] == 0)) 
// 	             MR_Step(b,i,j,k) = (MR_Step(b,i,j,k)>0.0)? Level:-Level;
// 	   
// 	   // range constraint on quantized-restored coefficients 
// 	   LevelQ =  Level / 2.;
//   	   Delta_w =  MR_Step(b,i,j,k) - TabSave[w];
// 	   if ((ABS(Delta_w) > LevelQ/2.)&&(TabSup[w]==1)) 
// 	     MR_Step(b,i,j,k) = (Delta_w>0.0)? TabSave[w]+LevelQ:TabSave[w]-LevelQ;
//            w++;
//        }
//        Delta_F = (F-F_old)/ w;
//  
//        // cerr << "Iter = " << Iter+1<< " F = " << F << "  Delta_F =" <<  Delta_F*Delta_F <<  endl;
//        F_old=F;
//        Iter++;
//    }
//    MR_Step.recons(Result);
//    // positivity constraint on reconstructed image
//    // threshold(Result);
// 
//    delete TabSup;
//    delete TabSave;  
// }
//  
/*********************************************************************/

int threshold_3d_band(MR_3D & MR_Data, int b, float Level)
{
   int i,j,k;
   int Nxb = MR_Data.size_band_nx(b);
   int Nyb = MR_Data.size_band_ny(b);
   int Nzb = MR_Data.size_band_nz(b);
   int Cpt = 0;
   
   for (i=0; i < Nxb; i++)
   for (j=0; j < Nyb; j++)
   for (k=0; k < Nzb; k++)
   {
      if (ABS(MR_Data(b,i,j,k)) < Level) MR_Data(b,i,j,k) = 0.;
      else Cpt ++;
   }
   return Cpt;
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    fitsstruct Header;
    fltarray Dat;
    int b,k;
    char Cmd[512];
    /* Get command line arguments, open input file(s) if necessary */

    extern softinfo Soft;

    Soft.mr3();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
      
    lm_check(LIC_MR3);
    transinit(argc, argv);
     
    io_3d_read_data(Name_Cube_In, Dat, &Header);
    Header.origin = Cmd;

    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();
    if (Verbose == True)
    {
       cout << "Filename in = " << Name_Cube_In << endl;
       cout << "Filename out = " << Name_Cube_Out  << endl;
       cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
    }
    
    MR_3D MR_Data;
    MR_Data.Verbose=Verbose;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if (Transform == TO3_MALLAT) 
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_Data.alloc (Nx, Ny, Nz, Transform, Nbr_Plan, PtrFAS, Norm);
          
    if (Verbose == True)
    {
       cout << "Transform = " << StringTransf3D((type_trans_3d) Transform) << endl;
      if (Transform == TO3_MALLAT)  
      {
         cout << "Filter = " <<  StringSBFilter(MR_Data.SBFilter) << endl;
         if (MR_Data.TypeNorm == NORM_L1) cout << "Norm L1 " << endl;
         else cout << "Norm L2 " << endl;
      }
      cout << "Number of scales = " << MR_Data.nbr_scale() << endl;
      if ((MAD == False) && (Noise_Ima > FLOAT_EPSILON))
                cout << "Sigma Noise = " << Noise_Ima << endl;
   }
    
    if (Verbose == True) cout << "3D transform ... " << endl;
    MR_Data.transform (Dat);
    if (Verbose == True) cout << "Thresholding ... " << endl;
    
    if (Noise_Ima < FLOAT_EPSILON)
    {
       if (Verbose == True) cout << "Automatic noise estimation ... " << endl;
       fltarray DatBand;
       MR_Data.get_band(0, DatBand);
       Noise_Ima = DatBand.sigma_clip();
       if (Verbose == True) cout << "Sigma noise = " << Noise_Ima << endl;
    }
    int Ndetect = 0;
    for (b=0; b < MR_Data.nbr_band()-1; b++)
    {
       float Level = Noise_Ima*N_Sigma;
       if (MAD == True)
       {
          int i,j,k;
          int N1 = MR_Data.size_band_nx(b);
	  int N2 = MR_Data.size_band_ny(b);
	  int N3 = MR_Data.size_band_nz(b);
          fltarray Buff(N1*N2*N3);
          int Ind = 0;             
	  for (i=0;i<N1;i++)
          for (j=0;j<N2;j++) 
	  for (k=0;k<N3;k++) Buff(Ind++) = ABS(MR_Data(b,i,j,k));
          Level = N_Sigma*get_median(Buff.buffer(), N1*N2*N3) / 0.6745;
	}      
	else if (Transform == TO3_ATROUS)
	{
	  Level *= TabNorm3DPaveBspline[b] / TabNorm3DPaveBspline[0];
	}
	int Cpt = threshold_3d_band(MR_Data, b, Level);
       Ndetect += Cpt;
       if (Verbose == True)
         cout << "Band " << b+1 << " Number of significant coefficients = " << Cpt << endl;
    }
    if (Verbose == True)
       cout << "Total number of significant coefficients = " << Ndetect << endl;
 
    if (Verbose == True) cout << "Reconstruction ... " << endl;
    
    MR_Data.recons(Dat);
    Header.origin = Cmd;
    io_3d_write_data(Name_Cube_Out, Dat, &Header); 
       
    exit(0);
} 

