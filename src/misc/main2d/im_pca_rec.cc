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
**    Date:  98/02/04
**    
**    File:  im_pca_rec.cc
**
*******************************************************************************
**
**    DESCRIPTION  images reconstruction from their principal components
**    ----------- 
**
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
// #include "MR_Obj.h"
#include "IM_Noise.h"
#include "NR.h"
#include "CPca.h"
#include "IM_Pca.h"

int TabKill[MAX_NBR_PCA_IMA]; 
char Name_In[MAX_NBR_PCA_IMA ][256];  /*  list of input image names  */
char NameImagOut[256]; /* file name used for the reconstructed image names */

int N_Image = 0 ;
Bool Verbose = False;
Bool WriteEigen = False;
Bool NormEigen = False;
Bool MeanSub = True;
int NbrUseEigen = -1;
Bool NoDisplay = False;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

/*************************************************************************/

static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input_image_list output_prefix_file_name\n", argv[0]);

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
//     fprintf(OUTMAN, "\n");
    
    fprintf(OUTMAN, "         [-F NbrEigenVect]\n");
    fprintf(OUTMAN, "             Number of eigen vector used for the reconstruction. \n");
    fprintf(OUTMAN, "             Default is set to the number of images (==> output=input).\n");     
    fprintf(OUTMAN, "\n");
    
    fprintf(OUTMAN, "         [-K EigenVect_Number]\n");
    fprintf(OUTMAN, "             Eigen vector number which will not be used for the reconstruction. \n"); 
    fprintf(OUTMAN, "\n");
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
    while ((c = GetOpt(argc,argv,"dvwNF:K:MzZ:")) != -1) 
    {
	switch (c) 
        {
	    case 'v':
		/* verbose flag -v */
		Verbose = True;
		break;
	     case 'M':
                MeanSub= False;
		break;            
	     case 'w':
                WriteEigen = True;
                break;
            case 'N':
                NormEigen = True;
 	    case 'd':
               NoDisplay = True;break;
           case 'F':
                if (sscanf(OptArg,"%d",&NbrUseEigen) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of eigen vectors: %s\n", OptArg);
		    usage(argv);
		}
                break;
            case 'K':
                if (sscanf(OptArg,"%d",&i) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad eigen vector number: %s\n", OptArg);
		    usage(argv);
		}
		if ((i < 1) || (i > MAX_NBR_PCA_IMA))
		{
		    fprintf(OUTMAN, "Error: bad eigen vector number: %s\n", OptArg);
                    fprintf(OUTMAN, "           0 < VectNbr <= %d\n", MAX_NBR_PCA_IMA);
		    usage(argv);
		}
		TabKill[i-1] = 1;
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
	if (N_Image > MAX_NBR_PCA_IMA)
	{
	   cerr << "Error: too many images to analyze ... " << endl;
	   cerr << "       the maximum is " << MAX_NBR_PCA_IMA << endl;
	}
    }
    if (N_Image < 3) usage(argv);
    N_Image --;
    
    strcpy(NameImagOut, Name_In[N_Image]);

    // set by default the number of eigen vectors used for the reconstruction
    // to the number of images
    if ((NbrUseEigen < 1) || (NbrUseEigen >= N_Image)) NbrUseEigen = N_Image;
    
    
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
  for (k =0; k < MAX_NBR_PCA_IMA; k++) TabKill[k] = 0;
  
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
   // CMA.MatCor.Verbose = True;
   CMA.MatCor.NormAna = NormEigen;

   // the subtract the mean of each image
   if (MeanSub == True) CMA.subtract_mean(TabIma);
   
   // calculate the correlation matrix and its eigen vectors
   CMA.compute(TabIma);
   
   // print the result to the standard output
   if (!NoDisplay) CMA.print();
   
   // transform the data into their principal components
   CMA.transform(TabIma, TabIma);
   
   //  save  the eigen vector images on the disk 
   if (WriteEigen == True)
   {
      char NameEigen[256];
      for(k=0; k < N_Image; k++)
      {
         sprintf(NameEigen, "pca_%d",k+1);
         io_write_ima_float(NameEigen, TabIma[k], &Header);
      }
   }
   
   // Calculate the number of eigen vector used for the
   // reconstruction
   int Ne = 0;
   for(k=0; k <  NbrUseEigen; k++) if (TabKill[k] == 0) Ne++;
   if (!NoDisplay)
      cout << "Number of eigen vector used for the reconstruction: Ne = " << Ne  << endl;
   if (Ne == N_Image)
      cout << "Warning: Ne = Nbr input images ==> the ouput is equal to the input ... " << endl;

   // apply the inverse reconstruction from a subset of eigen vectors
   CMA.invsubtransform(TabIma, TabIma, NbrUseEigen, TabKill);
   
   // add the mean value
   if (MeanSub == True) CMA.add_mean(TabIma);
   
   // save the results
   char NameOut[512];
   for(k=0; k < N_Image; k++)
   {
       sprintf(NameOut, "%s_%d",NameImagOut,k+1);
       io_write_ima_float(NameOut, TabIma[k], &Header);
   }
   exit(0);  
}
