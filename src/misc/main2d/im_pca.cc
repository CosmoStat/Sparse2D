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

#define MAX_NBR_PCA_IMA 10
char Name_In[MAX_NBR_PCA_IMA ][256];  /*  list of input image names  */
     
int N_Image = 0 ;
Bool Verbose = False;
Bool WriteEigen = False;
Bool NormEigen = False;
Bool MeanSub = True;
Bool WriteCorrelMatrix=False;
Bool NoDisplay = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

 
/*************************************************************************/

static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input_image_list\n", argv[0]);

    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN ,"         [-w]\n");
    fprintf(OUTMAN, "             Write to the disk the eigen vector images.\n");
    fprintf(OUTMAN, "               file names are: pca_x \n");
    fprintf(OUTMAN, "               where x is the eigen vector number.\n"); 
    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Do not subtract the mean. \n");
 
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
    while ((c = GetOpt(argc,argv,"dvwNMCzZ:")) != -1) 
    {
	switch (c) 
        {
             case 'C':
                 WriteCorrelMatrix = True;break;	    
	     case 'v':
		/* verbose flag -v */
		Verbose = True;
		break;
            case 'w':
                WriteEigen = True;
                break;
            case 'N':
                NormEigen = True;
		break;    
	    case 'd':
               NoDisplay = True;break;
	    case 'M':
                MeanSub= False;
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
    
    /* read input image names */
    i=0;
    N_Image = 0;
    while (OptInd < argc)
    {
	strcpy(Name_In[i], argv[OptInd++]);
	i++;
	N_Image++;
    }

    if (N_Image == 0) usage(argv);
    
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
  int i,k;
  int Nl=0, Nc=0;
  Ifloat   Dat_1, Dat_2;
  fitsstruct Header;
  extern void initfield(fitsstruct *Header);
  char Cmd[512];
 
  Cmd[0] = '\0';
  for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
  
  /* Get command line arguments*/
  lm_check(LIC_MR1);
  imfusinit(argc, argv);

  // print input image list
  if (Verbose == True)
  {
      cout << "Input images : " ;
      for(i=0; i< N_Image; i++) cout << Name_In[i] << ", ";
      cout << endl;
  }

  // allocate the memory for all images
  Ifloat *TabIma = new Ifloat [N_Image];
   
   // read the data, and verify that all images have the same size
   for(k=0; k < N_Image; k++) 
   {
      if (k == 0)
      {
         io_read_ima_float(Name_In[k],  TabIma[k], &Header);
         Nl = TabIma[0].nl();
         Nc = TabIma[0].nc();
         Header.origin = Cmd;
      }
      else
      {
         io_read_ima_float(Name_In[k],  TabIma[k]);
         {
           if (( TabIma[k].nl()!= TabIma[0].nl()) || ( TabIma[k].nc()!= TabIma[0].nc()))
           {
	      cerr << "Error : image of different size ..." << endl ;
	      cerr << "   image 1: " << Name_In[0]  << " " <<  TabIma[0].nl() << "X"  <<  TabIma[0].nc() << endl ;
	      cerr << "   image " << k+1 << ": " << Name_In[k]  << " " <<  TabIma[k].nl() << "X"  <<  TabIma[k].nc() << endl ;
	      exit(-1);
           }
         }
      }
   }
   
   
   // class for image correlation analysis
   CorrelMatAna2D CMA(N_Image);
   CMA.MatCor.Verbose = Verbose;
   CMA.MatCor.NormAna = NormEigen;

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
         sprintf(NameEigen, "pca_%d",k+1);
         io_write_ima_float(NameEigen, TabIma[k]);
      }
   }
   exit(0);
}
