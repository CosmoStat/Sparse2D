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
**    Date:  03/02/00
**    
**    File:  radon.cc
**
*******************************************************************************
**
**    DESCRIPTION  Radon program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Radon.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
type_radon RadonMethod = DEF_RADON;
int OutputNl=-1;  // Number of project
int OutputNc=-1;    // Number of pixel per scan
Bool Filter=False;  // apply a smoothing
Bool Reverse=False; // Reverse transform
int FilterWidth=-1; // Width filter size
float FilterSigma=-1; // Sigma parameter for the filtering
Bool TestProg=False;
int NbrIter= DEF_RADON_FSS_NBR_ITER;
Bool UseLeviCode = True;
Bool OptIter = True;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    manline();
    fprintf(OUTMAN, "         [-m type_of_radon_method]\n");
    for (int i = 1; i <= NBR_ALL_RADON; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,StringRadon((type_radon)(i-1)));
    fprintf(OUTMAN, "              default is %s.\n", StringRadon(RadonMethod));
    manline();    

    fprintf(OUTMAN, "         [-y OutputLineNumber]\n");
    fprintf(OUTMAN, "             For the RADON transform, OutputLineNumber = number of projection,\n"); 
    fprintf(OUTMAN, "               and default is twice the input image lines number.\n"); 
//    fprintf(OUTMAN, "             For the inv. RADON transform, OutputLineNumber = number of lines,\n");
//    fprintf(OUTMAN, "               and default is the input image column number.\n"); 
    manline();    

    fprintf(OUTMAN, "         [-x OutputColumnNumber]\n");
    fprintf(OUTMAN, "             For the RADON transform, OutputLineNumber = number of pixels per projection,\n"); 
    fprintf(OUTMAN, "               and default is the input image column number.\n"); 
//    fprintf(OUTMAN, "             For the inv. RADON transform, OutputLineNumber = number of column,\n");
//    fprintf(OUTMAN, "               and default is the input image column number.\n"); 
    manline(); 

    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse RADON transform.\n");  
    manline();  

    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Filter each scan of the Radon transform.\n");  
    manline();

    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Filter width.\n");  
    fprintf(OUTMAN, "             Default is 100.\n");  
    manline();  

    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             Sigma parameter for the filtering.\n");  
    fprintf(OUTMAN, "             Default is 10.\n");  
    manline();  
    fprintf(OUTMAN, "         [-R NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration for the Slant Stack reconstruction.\n");  
    fprintf(OUTMAN, "             Default is %d iterations.\n", NbrIter);  
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
    while ((c = GetOpt(argc,argv,"otlR:w:s:x:y:rm:fvzZ")) != -1) 
    {
	switch (c) 
        { 
	    case 'o': OptIter = (OptIter == True) ? False: True; break;
	    case 'l': UseLeviCode = (UseLeviCode == True) ? False: True;break;
	    case 't': TestProg=True;break;
            case 'r': Reverse=True;break;
            case 'f': Filter=True;break;
            case 'w': 
                if (sscanf(OptArg,"%d",&FilterWidth) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad filter size: %s\n", OptArg);
                    exit(-1);
                    
                }
                Filter=True;
                break;
            case 's': 
                if (sscanf(OptArg,"%f",&FilterSigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad sigma value: %s\n", OptArg);
                    exit(-1);
                    
                }
                Filter=True;
                break;
            case 'R': 
                if (sscanf(OptArg,"%d",&NbrIter) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of iterations: %s\n", OptArg);
                    exit(-1);
                    
                }
                break;
	    case 'y': 
                if (sscanf(OptArg,"%d",&OutputNl) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of lines: %s\n", OptArg);
                    exit(-1);
                    
                }
                break;
            case 'x': 
                if (sscanf(OptArg,"%d",&OutputNc) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of columns: %s\n", OptArg);
                    exit(-1);
                    
                }
                break;
           case 'm': 
                if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of radon method: %s\n", OptArg);
                    exit(-1);
                    
                }
                if ((c > 0) && (c <= NBR_ALL_RADON)) RadonMethod = (type_radon) (c-1);
                else  
                {
                    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
                    exit(-1);
                }
                 break;
           case 'v': Verbose = True;break;
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
       if ((Reverse == False) && (Filter == True))
       {
          cout << "Error: f,w,s options are only valid for the inverse transform ... " << endl;
	  exit(-1);
       }
       if ((OutputNl > 0) || (OutputNl > 0))
       {
          if (OutputNl <= 0)
          {
              cout << "Error: number of lines is not specified  ..." << endl;
              exit(-1);
          }
          if (OutputNc <= 0)
          {
              cout << "Error: number of columns is not specified ..." << endl;
              exit(-1);
          }
          if ((RadonMethod != RADON_PROJECT_BACK) && (RadonMethod != RADON_PROJECT_FFT))
          {
              cout << "Error: column and line number are fixed with this method ..." << endl;
              exit(-1);
          }
       }
       if ((Filter==True) && (RadonMethod != RADON_PROJECT_FFT))
       {
           cout << "Error: The filtering is not available for this Radon method ..." << endl;
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
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/***************/

int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR1);
    filtinit(argc, argv);


    Ifloat Data;
    Ifloat Result;
    int Nl,Nc;
 
    if (Reverse == False)
    {
        if (Verbose == True)
        { 
           cout << endl << endl << "PARAMETERS: " << endl << endl;
           cout << "File Name in = " << Name_Imag_In << endl;
           cout << "File Name Out = " << Name_Imag_Out << endl;  
           if (OutputNl > 0) cout << "OutputNl  = " << OutputNl << endl; 
           if (OutputNc > 0) cout << "OutputNc  = " << OutputNc << endl; 
	   cout << "Radon method = " <<  StringRadon(RadonMethod) << endl;
        }
 
       io_read_ima_float(Name_Imag_In, Data, &Header);
       Radon RD;
       RD.Verbose = Verbose;
       Header.origin = Cmd;
       Nl = Data.nl();
       Nc = Data.nc();
       if ((OutputNl > 0) && (OutputNc > 0))
             RD.alloc(Nl, Nc, RadonMethod, OutputNl, OutputNc);
        else RD.alloc(Nl, Nc, RadonMethod);
        Result.alloc(RD.radon_nl(), RD.radon_nc(), "rad");
        RD.SSR.Verbose = Verbose;
	RD.SSR.UseLeviCode = UseLeviCode;
        RD.transform(Data,Result);
	RD.write(Name_Imag_Out, Result);
       // Header.width = Result.nc();
        // Header.height = Result.nl();
        // Header.npix = Header.width*Header.height;
        // io_write_ima_float(Name_Imag_Out, Result, &Header);	
	if (TestProg == True)
	{
	   cout << endl << "Reconstruction ... " << endl;
	   Ifloat Rec(Nl, Nc,"rec");
	   RD.recons(Result, Rec);
           INFO_X(Rec, "Recons");
	   io_write_ima_float("xx_rec.fits", Rec);
	   Data -= Rec;
	   INFO_X(Data, "Err");
	   io_write_ima_float("xx_err.fits", Rec);
	}
    }
    else
    {
        Radon RD;
        RD.Verbose = Verbose;
        RD.read(Name_Imag_In, Data);       
 	if (Verbose == True)
	{
            cout << endl << endl << "PARAMETERS: " << endl << endl;
            cout << "File Name in = " << Name_Imag_In << endl;
            cout << "File Name Out = " << Name_Imag_Out << endl;  
            cout << "Radon method = " <<  StringRadon(RD.radon_method()) << endl;
	    if (Reverse == True) cout << "Inverse Radon transform  " <<  endl;  
            if (FilterSigma > 0) cout << "FilterSigma = " << FilterSigma<< endl; 
            if (FilterWidth > 0) cout << "FilterWidth  = " << FilterWidth << endl;  
        }
	
        // io_read_ima_float(Name_Imag_In, Data, &Header);
        Nl = Data.nl();
        Nc = Data.nc();
        if (Filter == True)
        {
           if (FilterSigma > 0) RD.FilterSigma = FilterSigma;
           if (FilterWidth > 0) RD.FilterWidth = FilterWidth ;
           RD.filter(Data);
        }
        // if ((OutputNl > 0) && (OutputNc > 0))
//              RD.alloc(OutputNl, OutputNc, RadonMethod, Nl, Nc);
//         else RD.alloc(Nc, Nc, RadonMethod, Nl, Nc);
        Result.alloc(RD.imag_nl(), RD.imag_nc(), "inv");
	RD.NbrIter = NbrIter;
	RD.SSR.OptIter = OptIter;
	RD.SSR.Verbose = Verbose;
	RD.SSR.UseLeviCode = UseLeviCode;
        RD.recons(Data, Result);
        Header.origin = Cmd;
        Header.width = Result.nc();
        Header.height = Result.nl();
        Header.bitpix = BP_FLOAT;  
        Header.npix = Header.width*Header.height;
        type_format FormatData = io_which_format(Name_Imag_Out);
        if ((FormatData == F_UNKNOWN) && (RD.FormatInputImag != F_UNKNOWN))  
                                           io_set_format(RD.FormatInputImag);
        io_write_ima_float(Name_Imag_Out, Result);
    }
    exit(0);
}
