/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  23/11/99
**    
**    File:  im3d_convert.cc
**
*******************************************************************************
**
**    DESCRIPTION  convert an image at a given format to another format
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/
 

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"

Bool NormByt = False;
int MaxVal=1;
int MinVal=0;

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
Bool FlipX=False;
Bool FlipY=False;
Bool Rotate=False;
int Nx=0;
int Ny=0;
int Nz=0;
int NSkip = 0;
Bool Ascii=False;

/****************************************************************************/

Bool CV_RGB_to_LUV = False;
Bool CV_RGB_to_HSV = False;
Bool CV_RGB_to_YUV = False;
Bool CV_LUV_to_RGB = False;
Bool CV_HSV_to_RGB = False;
Bool CV_YUV_to_RGB = False;
    
/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options image_in image_out\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Normalize the data between 0 and 255.\n");
    manline();
    
    fprintf(OUTMAN, "         [-x]\n");
    fprintf(OUTMAN, "             Flip x-axis.\n");
    manline();
        
    fprintf(OUTMAN, "         [-y]\n");
    fprintf(OUTMAN, "             Flip y-axis.\n");
    manline();
        
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Rotate.\n");
    manline();

    fprintf(OUTMAN, "         [-a Nx,Ny,Nz,skip]\n");
    fprintf(OUTMAN, "             The input file is an ascii file with size Nx x Ny x Nz.\n");
    fprintf(OUTMAN, "             skip is the number of lines to skip in the input file.\n");
    fprintf(OUTMAN, "             ex: -a 256,256,256,5 for a 256x256x256 cube, where the first five lines are not considered.\n");
    manline();
                
    vm_usage();
    manline();
    verbose_usage();
    fprintf(OUTMAN, "               \n\n");
    exit(0);
}

/****************************************************************************/

static void init(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
 
    /* get options */
    while ((c = GetOpt(argc,argv,"a:lhuLHUxyrbvzZ:")) != -1) 
    {
        switch (c) 
        {
          case 'a': 
	          if (sscanf(OptArg,"%d,%d,%d,%d",&Nx,&Ny,&Nz,&NSkip )  < 3) 
		  {
		     fprintf(OUTMAN, "Error: bad ascii file format ...\n");
                     exit(-1);
		  }
	          Ascii=True;
	          break;	  
	  case 'L': CV_RGB_to_LUV=True;break;
	  case 'H': CV_RGB_to_HSV=True;break;
	  case 'U': CV_RGB_to_YUV=True;break;
	  case 'l': CV_LUV_to_RGB = True;break;
	  case 'h': CV_HSV_to_RGB = True;break;
	  case 'u': CV_YUV_to_RGB = True;break;
	  case 'x': FlipX=True;break;
	  case 'y': FlipY=True;break;
	  case 'r': Rotate=True;break;
 	  case 'v': Verbose = True; break;
	  case 'b': NormByt=True; MinVal=0;MaxVal=255;break;
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

/****************************************************************************/

void io_read_cube_ascii(char *Name_Dat_In,  fltarray &Dat, int Nx, int Ny, int Nz, int NbrSkip)
{
   FILE *input=NULL;
   float Val;
   int c,i,j,k;
    
   input = fopen(Name_Dat_In,"r");
   if (input == NULL) 
   {
     cout << "Error: cannot open  file " <<  Name_Dat_In << endl;
     exit(-1);
   }
   for (i=0; i < NbrSkip; i++) while ( (c=getc(input)) != '\n' );

   Dat.alloc(Nx,Ny,Nz);
   for (k=0; k < Nz; k++)
   for (i=0; i < Ny; i++)
   for (j=0; j < Nx; j++)
   
   {
      if (fscanf(input,"%f",&Val) == 1) Dat(j,i,k) = Val;
      else 
      {
        cout << "Error: cannot read  in " <<  Name_Dat_In  << endl;
        exit(-1);
      }
   }
   fclose(input);
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    fitsstruct Header;
    fltarray Dat,DatOut;
    int i,j,k;
    Bool NewIma = False;
    char Cmd[512];
    extern softinfo Soft;

    Soft.mr3();
	    
    lm_check(LIC_MR3);
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    init(argc, argv);

    
    if (Ascii == True) io_read_cube_ascii(Name_Imag_In, Dat, Nx,  Ny, Nz, NSkip);
    else  io_3d_read_data(Name_Imag_In, Dat, &Header);
    
    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();
    
    if ((FlipX==True) || (FlipY==True) || (Rotate==True)) NewIma = True;
    
    if (Rotate == True) DatOut.alloc(Ny,Nx,Nz);
    else if (NewIma == True) DatOut.alloc(Nx,Ny,Nz);
    
    if (NormByt==True)
    {
       float Min = Dat.min();
       float Scale = (Dat.max() - Min) / (MaxVal - MinVal);
       for (i=0; i < Ny; i++) 
       for (j=0; j < Nx; j++)
       for (k=0; k < Nz; k++) 
 	             Dat(j,i,k) = (Dat(j,i,k) - Min) / Scale + MinVal;  
    }
    
    if (NewIma == True)
    {
       for (i=0; i < Ny; i++) 
       for (j=0; j < Nx; j++)  
       for (k=0; k < Nz; k++)
       {
          int Indi = (FlipY == False) ? i: Ny - i - 1;
	  int Indj = (FlipX == False) ? j: Nx - j - 1;
	  if (Rotate == True)
	  {
	     int x = Indj;
	     Indj = Ny - Indi - 1;
	     Indi = x;
	  }
          DatOut(Indj,Indi,k) = Dat(j,i,k);  
       } 
    }
    
    if (Verbose == True )
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Nx = " << Dat.nx();
       cout << " Ny = " << Dat.ny();
       cout << " Nz = " << Dat.nz()  << endl;
       cout << "CONVERT " << String3DFormat(io_detect_3dformat(Name_Imag_In));
       cout << " to " << String3DFormat(io_detect_3dformat(Name_Imag_Out)) << endl;
       if (NormByt==True)
         cout << "Normalization between 0 and 255" << endl;
    }
    
    if (CV_RGB_to_LUV ==True) rgb_to_luv(Dat);
    else if (CV_RGB_to_HSV ==True) rgb_to_hsv(Dat);
    else if (CV_RGB_to_YUV ==True) rgb_to_yuv(Dat);
    else if (CV_LUV_to_RGB ==True) luv_to_rgb(Dat);
    else if (CV_HSV_to_RGB ==True) hsv_to_rgb(Dat);
    else if (CV_YUV_to_RGB ==True) yuv_to_rgb(Dat);
    
    io_3d_set_format(Name_Imag_Out);
    if (NewIma == False)     io_3d_write_data(Name_Imag_Out, Dat);
    else  io_3d_write_data (Name_Imag_Out, DatOut, &Header);
    exit(0);
} 

