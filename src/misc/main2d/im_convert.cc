/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  im_convert.cc
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
Bool Ascii=False;
int Nla=0;
int Nca=0;
int NSkip = 0;

extern int FITS_HDU_Number;

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
    
    fprintf(OUTMAN, "         [-h FitsHduNbr]\n");
    fprintf(OUTMAN, "             FiTS HDU number.\n");
    manline();

    fprintf(OUTMAN, "         [-a Nl,Nc,skip]\n");
    fprintf(OUTMAN, "             The input file is an ascii file with Nl lines and Nc columns.\n");
    fprintf(OUTMAN, "             skip is the number of lines to skip in the input file.\n");
    fprintf(OUTMAN, "             ex: -a 256,256,5 for a 256x256 image, where the first five lines are not considered.\n");
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
    while ((c = GetOpt(argc,argv,"a:h:xyrbvzZ:")) != -1) 
    {
        switch (c) 
        {
          case 'a': 
	          if (sscanf(OptArg,"%d,%d,%d",&Nla,&Nca,&NSkip )  < 2) 
		  {
		     fprintf(OUTMAN, "Error: bad ascii file format ...\n");
                     exit(-1);
		  }
	          Ascii=True;
	          break;
          case 'h':  
                 if (sscanf(OptArg,"%d",&FITS_HDU_Number )  < 1) 
                 {
                    fprintf(OUTMAN, "Error: bad HDU number ...\n");
                    exit(-1);
                 } 
		 if (FITS_HDU_Number < 0)
		 {
		    fprintf(OUTMAN, "Error: bad HDU number ...\n");
                    exit(-1);
 		 }
 	         break;
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

void io_read_ima_ascii(char *Name_Dat_In,  Ifloat &Dat, int Nl, int Nc, int NbrSkip)
{
   FILE *input=NULL;
   float Val;
   int c,i,j;
    
   input = fopen(Name_Dat_In,"r");
   if (input == NULL) 
   {
     cout << "Error: cannot open  file " <<  Name_Dat_In << endl;
     exit(-1);
   }
   for (i=0; i < NbrSkip; i++) while ( (c=getc(input)) != '\n' );

   Dat.alloc(Nl, Nc);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
      if (fscanf(input,"%f",&Val) == 1) Dat(i,j) = Val;
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
    Ifloat Dat,DatOut;
    int i,j;
    Bool NewIma = False;
    
    lm_check(LIC_MR1);
    init(argc, argv);
 
    if (Ascii == True) io_read_ima_ascii(Name_Imag_In, Dat, Nla,  Nca, NSkip);
    else  io_read_ima_float(Name_Imag_In, Dat);
    
    
    if ((FlipX==True) || (FlipY==True) || (Rotate==True)) NewIma = True;
    
    if (Rotate == True) DatOut.alloc(Dat.nc(),Dat.nl(),"CONV");
    else if (NewIma == True) DatOut.alloc(Dat.nl(),Dat.nc(),"CONV");
    
    if (NormByt==True)
    {
       float Min = min(Dat);
       float Scale = (max(Dat) - Min) / (MaxVal - MinVal);
       for (i=0; i < Dat.nl(); i++) 
       for (j=0; j < Dat.nc(); j++)  
	             Dat(i,j) = (Dat(i,j) - Min) / Scale + MinVal;  
    }
    
    if (NewIma == True)
    {
       for (i=0; i<  Dat.nl(); i++) 
       for (j=0; j<  Dat.nc(); j++)  
       {
          int Indi = (FlipY == False) ? i: Dat.nl() - i - 1;
	  int Indj = (FlipX == False) ? j: Dat.nc() - j - 1;
	  if (Rotate == True)
	  {
	     int x = Indj;
	     Indj = Dat.nl() - Indi - 1;
	     Indi = x;
	  }
          DatOut(Indi,Indj) = Dat(i,j);  
       } 
    }
           
    
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Nl = " << Dat.nl() << endl;
       cout << "  Nc = " << Dat.nc() << endl << endl;
       cout << "CONVERT " << StringFormat(io_detect_format(Name_Imag_In));
       cout << " to " << StringFormat(io_detect_format(Name_Imag_Out)) << endl;
       if (NormByt==True)
         cout << "Normalization between 0 and 255" << endl;
    }
    io_set_format(Name_Imag_Out);
    if (NewIma == False) io_write_ima_float (Name_Imag_Out, Dat);
    else  io_write_ima_float (Name_Imag_Out, DatOut);
    exit(0);
} 

