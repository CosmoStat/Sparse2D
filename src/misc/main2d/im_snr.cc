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
**    Date:  03/09/99
**    
**    File:  tfil.cc
**
*******************************************************************************
**
**    DESCRIPTION  filtering program
**    ----------- 
**                 
******************************************************************************/
   


#include "IM_Obj.h"
#include "IM_IO.h"
 
extern int  GetOpt(int argc, char **argv, char *opts);
extern int  OptInd;
extern char *OptArg;


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
Bool Verbose=False;
Bool SnrVal=False;
Bool Write=False;

#define NBR_ERROR_TYPE 5
enum error_type {LEAST_SQUARE_ERROR=1,
                 NORMALIZED_LEAST_SQUARE_ERROR=2,
                 PEAK_LEAST_SQUARE_ERROR=3,
		 PEAK255_LEAST_SQUARE_ERROR=4,
		 VAR_REF_ON_VAR_ERROR=5};
		 
error_type ErrorType = PEAK255_LEAST_SQUARE_ERROR;

inline char * StringErrorType (error_type ET) {
   switch (ET) {
      case LEAST_SQUARE_ERROR :            return ("Least Square Error (db) = -10 log10(Variance(Error))");break;
      case NORMALIZED_LEAST_SQUARE_ERROR : return ("Normalized LSE (db)     = -10 log10(Total(Error^2) / Total(Ref2))");break;
      case PEAK_LEAST_SQUARE_ERROR :       return ("Peak LSE (db)           = -10 log10(Variance(Error) / Max(Ref1^2) ");break;
      case VAR_REF_ON_VAR_ERROR :          return ("VRVE                    = Variance(Reference) / Variance(Error)");break;
      case PEAK255_LEAST_SQUARE_ERROR :    return ("PLSE255 (db)            = -10 log10(Variance(Error)/255^2)");break;
      default: return ("Undefined transform");break;  
   }
}
   
     
/***************************************/

static void usage(char *argv[]) {


   
    fprintf(OUTMAN, "Usage: %s options Reference Data\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    
    fprintf(OUTMAN, "         [-s SNR type]\n"); 
    for (int i=1;i<=NBR_ERROR_TYPE;i++)              
       fprintf(OUTMAN, "              %d: %s \n",i,StringErrorType((error_type)i));
    fprintf(OUTMAN, "             default is 4.\n");       
    verbose_usage();
 
    fprintf(OUTMAN, "         [-w]\n"); 
    fprintf(OUTMAN, "            write the result in file SNR.fits\n"); 
    
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void infinit(int argc, char *argv[]) {
    
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    int v;
    while ((c = GetOpt(argc,argv,"ws:kvzZ")) != -1) {
	
	switch (c) {  
        case 'w': Write = True;break;
	case 'v': Verbose = True;break;
	case 'k': SnrVal = True;break;
	case 's':
	        if (sscanf(OptArg,"%d",&v ) != 1) {
		   fprintf(OUTMAN, "bad type of error choice!!! : %s\n", OptArg);
	           exit(-1);
		}
		if ((v<1) || (v>NBR_ERROR_TYPE)) ErrorType = LEAST_SQUARE_ERROR;
		else ErrorType = (error_type)v;
		cout << StringErrorType(ErrorType) << endl;
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
       

    /* get optional input file names from trailing parameters and open files */
    if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
    else usage(argv);

    if (OptInd < argc)  strcpy(Name_Imag_Out, argv[OptInd++]);
    else usage(argv);
 
    /* make sure there are not too many parameters */
    if (OptInd < argc) {
       fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
       exit(-1);
    }

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/****************************************************************************/

int main(int argc, char *argv[]) {
 
     
    /* Get command line arguments, open input file(s) if necessary */
    // lm_check(LIC_MR1); 
    infinit(argc, argv);
    Ifloat IRef,IDat;

    fltarray Data;
    fltarray Ref;
    
    if (io_detect_format(Name_Imag_In) == F_FITS)
    {
       fits_read_fltarr(Name_Imag_In,  Ref);
       fits_read_fltarr(Name_Imag_Out, Data);
    }
    else
    {
       io_read_ima_float(Name_Imag_In, IRef);
       io_read_ima_float(Name_Imag_Out, IDat);
       Ref.alloc(IRef.buffer(), IRef.ny(), IRef.nx());
       Data.alloc(IDat.buffer(), IDat.ny(), IDat.nx());
    }
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    
    if ((Nx != Ref.nx()) || (Ny != Ref.ny()) || (Nz != Ref.nz())) {
       cerr << "Data must have same size!!" << endl;
       exit (-1);
    }
    
    int SizeOneD = (Nx > 0 ? Nx : 1) * (Ny > 0 ? Ny : 1) * (Nz > 0 ? Nz : 1);
    Data.reform(SizeOneD);
    Ref.reform(SizeOneD);
    
    fltarray Error(SizeOneD);
    fltarray Error2(SizeOneD);
    fltarray Ref2(SizeOneD);    

    for (int i=0;i<SizeOneD;i++) {
       Error(i)  = Ref(i)-Data(i);
       Error2(i) = Error(i) * Error(i); 
       Ref2(i)   = Ref(i) * Ref(i);
    }
    
    float MeanRef     = Ref.mean();
    float SigmaRef    = Ref.sigma();
    float MeanData    = Data.mean();
    float SigmaData   = Data.sigma();   
    //float MeanError   = Error.mean();
    float SigmaError  = Error.sigma();
    //float MeanError2  = Error2.mean();
    //float SigmaError2 = Error2.sigma();   

//cout << MeanRef << " | " << SigmaRef << " | " << MeanData << " | " << SigmaData << endl;
//cout << MeanError << "|" << SigmaError << "|" << MeanError2 << "|" << SigmaError2 << endl;
  
    float SNR;
    float M2 = Ref.max();
    M2 *= M2;
    float LSE  = -10.*log10(Error2.total());
    float NLSE = -10.*log10(Error2.total()/Ref2.total());
    float PLSE = -10.*log10(Error2.total()/M2);
    float P255LSE =  -20.*log10(SigmaError/255.);
    float VARERR = -10.*log(SigmaError*SigmaError/SigmaRef/SigmaRef)/log(10.);
    switch (ErrorType) 
    {
      case LEAST_SQUARE_ERROR : SNR = LSE; break;
      case NORMALIZED_LEAST_SQUARE_ERROR : SNR = NLSE;break;
      case PEAK_LEAST_SQUARE_ERROR: SNR = PLSE; break;
      case PEAK255_LEAST_SQUARE_ERROR: SNR = P255LSE;  break;
      case VAR_REF_ON_VAR_ERROR: SNR = VARERR; break;
      default : cerr << "bad error type..." << endl; exit(-1);break;
   }
   
   // Universal Image Quality Index (UIQI)
   float Sm = MeanRef;
   float VS = 0;
   float SRm = MeanData;
   float VSR = 0;
   float VSSR=0.;
   for (int i=0;i<SizeOneD;i++)  
   {
      VSSR  += (Ref(i)-Sm) * (Data(i) - SRm);
      VS += (Ref(i)-Sm)*(Ref(i)-Sm);
      VSR += (Data(i) - SRm)*(Data(i) - SRm);
   }
   VSSR /= (float) (SizeOneD-1);
   VS /= (float) (SizeOneD-1);
   VSR /= (float) (SizeOneD-1);
   float UIQI= 4*VSSR*Sm*SRm / ( (VS+VSR)*(Sm*Sm+SRm*SRm));
   
     
   if (Verbose == True) 
   {
       cout << "   Sigma(Err) " << Error.sigma() << " Min(Err) = " << Error.min() << " Max(Err) = " << Error.max() << endl;
      cout << "    LSE     = " << LSE << endl;
      cout << "    NLSE    = " << NLSE << endl;
      cout << "    PLSE    =  " << PLSE << endl;
      cout << "    P255LSE = " <<  P255LSE << endl;
      cout << "    VR/VE   = " <<  VARERR << endl;
      cout << "    UIQI   = " <<  UIQI << endl;
   }
   else cout << "SNR (" << StringErrorType(ErrorType) << ") = "  << SNR << " db" << endl; 
   fltarray SNRRes(1);
   SNRRes(0) = SNR;
   if (Write == True)
      fits_write_fltarr("SNR.fits", SNRRes);
    
   exit(0);
}
