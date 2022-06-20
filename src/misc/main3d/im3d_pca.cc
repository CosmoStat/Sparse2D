/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Yves Bobichon && J.L. Starck
**
**    Date:  98/02/03
**    
**    File:  im_pca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image Principal Component Analysis
**    ----------- 
**                 Image Principal Component Analysis (see 'Analysis of rapid 
**                 variation in the spectra of alpha Col by Cross Correlation'
**                 A. Bijaoui and V. Doazan, A&A 70, pp 285-291,  1979).
**
**
**    PARAMETRES    
**    ----------    
** 
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "NR.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "IM3D_IO.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char RecName[256]; /* output file name */
int N_Image = 0 ;
Bool Verbose = False;
Bool WriteEigen = True;
Bool NormEigen = False;
Bool MeanSub = True;
Bool WriteCorrelMatrix=False;
Bool NoDisplay = True;
int NbrUseEigen = 1;
Bool Rec = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

 
/*************************************************************************/

static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input_cube output_cube\n", argv[0]);

    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Do not subtract the mean. \n");

    fprintf(OUTMAN, "         [-r ReconstructedFileName]\n");
    fprintf(OUTMAN, "             Reconstruct the data from a subset of eigen vector. \n"); 

    fprintf(OUTMAN, "         [-F NbrEigen]\n");
    fprintf(OUTMAN, "             Number of eigen vectors used in the reconstruction. \n");
    fprintf(OUTMAN, "             Default is 1.\n"); 

//     fprintf(OUTMAN, "         [-N]\n");
//     fprintf(OUTMAN, "             Normalize the eigen vectors. \n");
//     fprintf(OUTMAN, "             Default is no.\n");     
        
//     fprintf(OUTMAN, "         [-v]\n");
//     fprintf(OUTMAN, "              Verbose. \n");
//     fprintf(OUTMAN, "              Default is no\n");     
    vm_usage();
    manline();
    verbose_usage();    
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void imfusinit(int argc, char *argv[])
{
    int c,  i;
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"r:F:dvNMCzZ:")) != -1) 
    {
	switch (c) 
        {
           case 'F':
                if (sscanf(OptArg,"%d",&NbrUseEigen) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of eigen vectors: %s\n", OptArg);
		    usage(argv);
		}
                break;
           case 'C':
                 WriteCorrelMatrix = True;break;	    
	     case 'v':
		/* verbose flag -v */
		Verbose = True;
		break;
                  break;
            case 'N':
                NormEigen = True;
		break;    
	    case 'd':
               NoDisplay = True;break;
	    case 'M':
                MeanSub= False;
		break;
	    case 'r':
	         if (sscanf(OptArg,"%s", RecName) != 1)
		{
		   fprintf(OUTMAN, "Error: Bad file name ... \n");
		   exit(-1);
		}
		Rec = True;
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
    if ((Rec == True) && (NbrUseEigen < 1))
    {
       fprintf(OUTMAN, "Error: Nbr of eigen must be specified (-F option) ...\n");
       exit(-1); 
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
	fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
     }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}


/*********************************************************************/

int main(int argc, char *argv[]) 
{
  int i,j,k;
  int Nl=0, Nc=0;
  Ifloat   Dat_1, Dat_2;
  fitsstruct Header;
  extern void initfield(fitsstruct *Header);
  char Cmd[512];
 
  Cmd[0] = '\0';
  for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
  
  /* Get command line arguments*/
  lm_check(LIC_MR3);
  imfusinit(argc, argv);

  fltarray Dat,DatOut;
  io_3d_read_data(Name_Imag_In, Dat, &Header);
  int Nx = Dat.nx();
  int Ny = Dat.ny();
  int Nz = Dat.nz();
        
  // allocate the memory for all images
  Ifloat *TabIma = new Ifloat [Nz];
  N_Image = Nz;
  Nl = Ny;
  Nc = Nx;
  
   for(k=0; k < N_Image; k++) 
   {
       TabIma[k].alloc(Nl,Nc);
       for (i=0; i < Nl; i++)
       for (j=0; j < Nc; j++) TabIma[k](i,j) = Dat(j,i,k);
   }
   
   
   // class for image correlation analysis
   CorrelMatAna2D CMA(N_Image);
   CMA.MatCor.Verbose = Verbose;
   CMA.MatCor.NormAna = NormEigen;
   CMA.MatCor.SortEigen = True;

   // the subtract the mean of each image
   if (MeanSub == True) CMA.subtract_mean(TabIma);
   
   // calculate the correlation matrix and its eigen vectors
   CMA.compute(TabIma);
   
   // print the result to the standard output
   if (!NoDisplay) CMA.print();
   if (WriteCorrelMatrix) {
      fltarray CM = fltarray (N_Image,N_Image);
      for(int i=0; i<N_Image; i++) {
	 for(int j=0; j<N_Image; j++) {
	    CM(i,j) = CMA.CorrelMat(i,j);
	 }
      } 
      fits_write_fltarr ("Correl_Matrix.fits", CM);
   }
     
   // calculate the eigen vector images and save them on the disk 
   if (WriteEigen == True)
   {
      char NameEigen[256];
      CMA.transform(TabIma, TabIma);
      for(k=0; k < N_Image; k++)
      {
          for (int i=0; i < Nl; i++)
          for (int j=0; j < Nc; j++) Dat(j,i,k) = TabIma[k](i,j);
      }
   }
   io_3d_write_data (Name_Imag_Out, Dat, &Header);
   
   if ((Rec == True) && (NbrUseEigen > 0))
   { 
       CMA.invsubtransform(TabIma, TabIma, NbrUseEigen);
       if (MeanSub == True) CMA.add_mean(TabIma);
       for(k=0; k < N_Image; k++)
       {
          for (int i=0; i < Nl; i++)
          for (int j=0; j < Nc; j++) Dat(j,i,k) = TabIma[k](i,j);
       }
       io_3d_write_data (RecName, Dat, &Header);
   }
   
   exit(0);
}
