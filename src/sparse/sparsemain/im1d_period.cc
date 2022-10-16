/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre & Jean-Luc Starck
**
**    Date:  19/03/03
**    
**    File:  im1d_period.cc
**
*******************************************************************************
**
**    DESCRIPTION  periodogram  calcualation
**    ----------- 
**                 
******************************************************************************/
   


#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "FFTN_1D.h"
 
extern int  GetOpt(int argc, char **argv, char *opts);
extern int  OptInd;
extern char *OptArg;


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
Bool Verbose=False;
Bool SubstractMean = False;
int DanielHalfSize=5;
int BarlettSegmentNumber=8;
int WelchShift=0;
int WelchNbEchant=0;
int NbPointOut=0;
float WinParam = 0.5;
type_std_win WinType = W_HANNING;

#define NBR_PERIODOGRAM_TYPE 4 
enum periodogram_type {STANDARD=1,
                       BARLETT=2,
                       WELCH=3,
		       DANIEL=4};
periodogram_type PeriodogramType=WELCH;		 

inline char * StringPeriodogramType (periodogram_type PT) 
{
   switch (PT) 
   {
      case STANDARD :   return ("Standard periodogram");break;
      case BARLETT  :   return ("Barlett periodogram");break;
      case WELCH    :    return ("Welch periodogram");break;
      case DANIEL   :    return ("Daniel periodogram");break;
      default: return ("Undefined periodogram type");
      exit(-1);break;  
   }
}
   
     
/***************************************/

static void usage(char *argv[]) {


   
    fprintf(OUTMAN, "Usage: %s options Reference Data\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    
    fprintf(OUTMAN, "         [-m Periodogram method]\n"); 
    for (int i=1;i<=NBR_PERIODOGRAM_TYPE;i++)              
       fprintf(OUTMAN, "             %d: %s \n",i,StringPeriodogramType((periodogram_type)i));
    fprintf(OUTMAN, "             default is 3.\n");   
    manline();
    
    fprintf(OUTMAN, "         [-h HalfSizeWindow]\n");         
    fprintf(OUTMAN, "             Window half size used in the Daniel periodogram.\n");         
    fprintf(OUTMAN, "             default is 5.\n");         
    manline();
        
    fprintf(OUTMAN, "         [-q SegmentNbr]\n");         
    fprintf(OUTMAN, "             Number of segment used in the Barlett periodogram.\n");         
    fprintf(OUTMAN, "             default is 8\n");         
    manline(); 
       
    fprintf(OUTMAN, "         [-r Shift]\n"); 
    fprintf(OUTMAN, "             Shift used in Welch periodogram.\n");        
    fprintf(OUTMAN, "             default is signal_point_number / 10\n");         
    manline(); 
    
    fprintf(OUTMAN, "         [-s SegmentNbr]\n");         
    fprintf(OUTMAN, "             Number of segment used in Welch periodogram.\n");
    fprintf(OUTMAN, "             default is signal_point_number / 2\n");         
    manline();
    
    fprintf(OUTMAN, "         [-t PeriodPointNbr]\n");
    fprintf(OUTMAN, "             Number of point of periodogram.\n");
    fprintf(OUTMAN, "             It must be greater or equal to signal_point_number\n");                
    fprintf(OUTMAN, "             default is signal_point_number.\n");         
    manline();  
    
    fprintf(OUTMAN, "         [-T type_of_window ]\n");
    for (int i = 0; i < NBR_STD_WIN ; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringSTDWin((type_std_win )i));
    fprintf(OUTMAN, "             default is %s.\n", StringSTDWin(WinType));
    manline();
   
    fprintf(OUTMAN, "         [-W window_param (time domain)]\n");
    fprintf(OUTMAN, "             Window parameter. Default is 0.5.\n");
    manline();
    
    fprintf(OUTMAN, "         [-M]");
    fprintf(OUTMAN, "             Substract mean (default is no)");
    manline(); 
       
    verbose_usage();
    
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
    while ((c = GetOpt(argc,argv,"m:h:q:r:s:t:T:W:MvzZ")) != -1) {
	
	switch (c) {  
	case 'v': Verbose = True; break;
	case 'M': SubstractMean = True; break;
	case 'm':
	        if (sscanf(OptArg,"%d",&v ) != 1) {
		   fprintf(OUTMAN, "bad periodogram type !!! : %s\n", OptArg);
	           exit(-1);
		}
		if ((v<1) || (v>NBR_PERIODOGRAM_TYPE)) 
		   PeriodogramType = WELCH;
		else PeriodogramType = (periodogram_type)v;
		cout << StringPeriodogramType(PeriodogramType) << endl;
		break;
        case 'h':
	        if (sscanf(OptArg,"%d",&DanielHalfSize) != 1) {
		   fprintf(OUTMAN, "bad Daniel periodogram half size value!!! : %s\n", 
		           OptArg);
	           exit(-1);
		}
		break;
        case 'q':
	        if (sscanf(OptArg,"%d",&BarlettSegmentNumber) != 1) {
		   fprintf(OUTMAN, "bad Barlett segment number!!! : %s\n", 
		           OptArg);
	           exit(-1);
		}
		break;
        case 'r':
	        if (sscanf(OptArg,"%d",&WelchShift) != 1) {
		   fprintf(OUTMAN, "bad Welch shift!!! : %s\n", 
		           OptArg);
	           exit(-1);
		}
		break;
        case 's':
	        if (sscanf(OptArg,"%d",&WelchNbEchant) != 1) {
		   fprintf(OUTMAN, "bad Welch point number!!! : %s\n", 
		           OptArg);
	           exit(-1);
		}
		break;
        case 't':
	        if (sscanf(OptArg,"%d",&NbPointOut) != 1) {
		   fprintf(OUTMAN, "bad out point number!!! : %s\n", 
		           OptArg);
	           exit(-1);
		}
		break;
        case 'T':
	        if (sscanf(OptArg,"%d",&v) != 1) {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
		    exit(-1);
		}
                if ((v <= 0) && (v >  NBR_STD_WIN)) {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
	            usage(argv);
 		}
		WinType = (type_std_win) (v-1);
                break;
        case 'W':
	        if (sscanf(OptArg,"%f",&WinParam) != 1) {
		    fprintf(OUTMAN, "bad window parameter: %s\n", OptArg);
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

int main(int argc, char *argv[]) 
{
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1); 
    infinit(argc, argv);
    Ifloat IDat;
    fltarray Data;
    
    io_1d_read_data(Name_Imag_In, Data);
    int Nx = Data.nx();
    if (SubstractMean) {
       float Mean = Data.mean();
       for (int i=0;i<Nx;i++) 
          Data(i) -= Mean;
    }
    
    if ((Data.ny() != 0) || (Data.nz() != 0)) 
    {
       cerr << "Error: input data must be one dimension data ..." << endl;
       exit (-1);
    }
    
    int NbEchant;
    int Decal;
    int NbPeriod;
    fltarray Window;
    float NrjWind=1;
    
    switch (PeriodogramType) 
    {  
       case STANDARD :
          NbEchant = Nx;
	  Decal    = Nx;
	  NbPeriod = 1; 
	  Window.alloc(NbEchant);
	  Window.init(1);  
          break;
       case BARLETT :
          if (Verbose) 
	     cout << "Segment number : " << BarlettSegmentNumber << endl;
          NbEchant = Nx/BarlettSegmentNumber;
	  Decal    = NbEchant;
	  NbPeriod = BarlettSegmentNumber;	  
	  Window.alloc(NbEchant);
	  Window.init(1);              
          break;
       case WELCH : 
          if (WelchNbEchant == 0) WelchNbEchant=Nx/2; 
          if (WelchShift == 0) WelchShift=Nx/10;
	  if (Verbose) 
	  {
	     cout << "Shift number : " << WelchShift << endl;
	     cout << "Point number : " << WelchNbEchant << endl;	     
	  }
	  NbEchant = WelchNbEchant;
	  Decal    = WelchShift;
	  NbPeriod = (Nx-WelchNbEchant)/WelchShift + 1;
	  
	  Window.alloc(NbEchant);
	  STD_Window STD;
	  STD.make_win (Window, WinType, WinParam);
	  NrjWind = 0.;
	  for (int i=0;i<NbEchant;i++) NrjWind += Window(i)*Window(i);
	  NrjWind = NrjWind/NbEchant;
	  if (Verbose) 
	  {
	     cout << "Type of window : " << StringSTDWin(WinType) << endl;
	     cout << "Param of window : " << WinParam << endl;
	     cout << "Energy of window : " << NrjWind << endl;
	  }
          break;
       case DANIEL : 
          if (Verbose)
	     cout << "Half size window used in Daniel periodogram" 
	          << DanielHalfSize << endl;         
	  NbEchant = Nx;
	  Decal    = Nx;
	  NbPeriod = 1; 
	  Window.alloc(NbEchant);
	  Window.init(1);  
          break;
       default: cout << "Undefined periodogram type";exit(-1);break;  
    }
    if (NbPointOut < Nx) NbPointOut=Nx;
    
    if (Verbose) {
       cout << "Number of signal point used in one periodogram : " << NbEchant << endl;
       cout << "shift between periodograms : " << Decal << endl;
       cout << "Number of periodogram : " << NbPeriod << endl;
       cout << "Number of point in periodogram : " << NbPointOut << endl;
    }    
   
 
    // compute presiodograms...
    //-------------------------
    fltarray TabDsp(NbPeriod, MAX(NbEchant,NbPointOut));
    fltarray LocSig(MAX(NbEchant,NbPointOut));
    complex_f* idsp = new complex_f [NbPointOut];
    fltarray Periodogram(NbPointOut);
    FFTN_1D CFFT1D;

    CFFT1D.CenterZeroFreq = True;

    for (int i=0;i<NbPeriod;i++) 
    {
       if (Verbose) cout << "periodogram number : " << i+1 << "/" << NbPeriod << endl;;
       for (int j=0;j<NbEchant;j++) 
          LocSig(j)=Data(i*Decal+j)*Window(j);
       if (NbPointOut > NbEchant)    // add sero (zero padding)
          for (int j=NbEchant;j<NbPointOut;j++) LocSig(j)=0.;
       CFFT1D.fftn1d(LocSig, idsp);
       for (int j=0;j<NbPointOut;j++) 
          TabDsp(i,j)=(pow((double) real(idsp[j]),(double)2.) + pow((double) imag(idsp[j]),(double)2.));
    }
    
    for (int i=0;i<NbPeriod;i++)
    for (int j=0;j<NbPointOut;j++)
          Periodogram(j) += TabDsp(i,j);
    for (int j=0;j<NbPointOut;j++)
       Periodogram(j) *= 1./NbPeriod/NrjWind/NbPointOut/NbPointOut;
    
    if (PeriodogramType == DANIEL) 
    {
       fltarray SmothPeriod(NbPointOut); 
       //for (int i=DanielHalfSize;i<NbPointOut-DanielHalfSize;i++)
       for (int i=0;i<NbPointOut;i++)
       for (int j=-DanielHalfSize;j<=DanielHalfSize;j++)
	     SmothPeriod(i) += Periodogram(i+j, I_CONT);
       for (int i=0;i<NbPointOut;i++)  
          Periodogram(i) = 1./(2*DanielHalfSize+1)*SmothPeriod(i);
    }
    
    // fits_write_fltarr(Name_Imag_Out, Periodogram);
   io_1d_write_data(Name_Imag_Out, Periodogram);

   exit(0);
}
