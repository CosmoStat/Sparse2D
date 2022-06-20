/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/03
**    
**    File:  im_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  filter an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: im_filter option image output
**        where options = 
**          [-f type_of_filtering]
**                 1: Median filtering 
**                 2: Average filtering 
**                 3: Bspline filtering 
**                 4: Anisotropic diffusion filtering 
**
**                 default is 2
**    
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
  
Bool PositivIma = True;
Bool MaxIma = False;
 
int SmoothWindowSize = 5;
type_border Border = DEFAULT_BORDER;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool PosOpt=False;
Bool WindowOpt=False;    
Bool Verbose=False;
    
#define NBR_IM_FILTERING 3
#define DEF_MAX_IN_IMA 255

enum type_im_filter {F_MEDIANE, F_AVERAGE, F_BSPLINE, F_ANISOTROPIC};
type_im_filter Filter = F_MEDIANE;  

static char * StringIMFilter (type_im_filter type)
{
    switch (type)
    {
        case F_MEDIANE: 
              return ("Median Filtering");break;
        case F_AVERAGE:
              return ("Moving Average Filtering");break;
        case F_BSPLINE:
              return ("Bspline Filtering");break;
        case F_ANISOTROPIC:
              return ("Anisotropic diffusion filtering");break;
    }
   return ("Error: bad type of filtering");
}

static void im_filter_usage(type_im_filter Filter)
{
    fprintf(OUTMAN, "         [-f type_of_filtering]\n");
    for (int i = 0; i < NBR_IM_FILTERING; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringIMFilter((type_im_filter)i));
    fprintf(OUTMAN, "              default is %s.\n", StringIMFilter(Filter));
}

static  void MaxThreshold (Ifloat &Image) {
    int i,j;
    for (i = 0; i < Image.nl(); i++)
    for (j = 0; j < Image.nc(); j++) if (Image(i,j) > DEF_MAX_IN_IMA) Image(i,j) =  DEF_MAX_IN_IMA;
}
 
/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    im_filter_usage(Filter);
    manline(); 
   
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
    manline();
    
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Add the maximum level constraint.\n");
    fprintf(OUTMAN, "             Max value is 255. Default is no.\n");
    manline();   

//    fprintf(OUTMAN, "         [-G RegulParam]\n");
//    fprintf(OUTMAN, "              Regularization parameter for the TV method.\n");
//    fprintf(OUTMAN, "              default is %f\n", RegulParam);
//    manline();
     
    window_size_usage(SmoothWindowSize);    
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
    while ((c = GetOpt(argc,argv,"pbm:f:t:PW:vzZ:G:")) != -1) 
    {
	switch (c) 
        {
//  	   case 'G':
// 		/* -s <nsigma> */
// 		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
//                 {
// 		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
// 		    exit(-1);
// 		}
//                 if (RegulParam  < 0.) RegulParam = 0.1;
// 		OptRegul = True;
//  		break;	   
	   case 'f':
		/* -f <type> type of filtering */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_IM_FILTERING)) Filter = (type_im_filter) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
 		}
 		break;
            case 'v': Verbose = True;break;
            case 'P':
                PositivIma=False;
                break;            
	    case 'b':
                MaxIma=True;
               break;
	    case 'W':
		if (sscanf(OptArg,"%d", &SmoothWindowSize) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad window size: %s\n", OptArg);
		   exit(-1);
		}
		WindowOpt = True;
		if (SmoothWindowSize < 3)
		{
		   fprintf(OUTMAN, "Error: bad window size parameter: %s \n", OptArg);
		   fprintf(OUTMAN, "           SmoothWindowSize  > 2\n");
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


/****************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    Ifloat Data;
     
     /* support image creation */
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
       cout << "Filter = " << StringIMFilter(Filter) << endl;
       if (PositivIma == True) cout << "Positivity constraint" << endl;
       else cout << "No positivity constraint" << endl;    
       if (MaxIma == True) cout << "Maximum level constraint" << endl;
       else cout << "No maximum level constraint" << endl;
    }
 
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    Ifloat Result (Data.nl(), Data.nc(), "Result Filtering");
    switch (Filter)
    {
        case F_MEDIANE:
 	     smooth_mediane (Data, Result, Border, 0, SmoothWindowSize);
            break;
       case F_AVERAGE:
              smooth_average (Data, Result,  Border, 0, SmoothWindowSize); 
             break;
        case F_BSPLINE:
            smooth_bspline (Data, Result,  Border, 0);
            break;
        default:
           cerr << "Error: this filtering method is not implemented in this routine ... " << endl;
           exit(-1);
           break;
    }
    
    if (PositivIma == True) threshold(Result);
    if (MaxIma == True) MaxThreshold(Result);  
    io_write_ima_float(Name_Imag_Out, Result, &Header);
    exit(0);
} 

