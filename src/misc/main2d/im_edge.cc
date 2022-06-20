/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  97/08/22
**    
**    File:  im_edge.cc
**
*******************************************************************************
**
**    DESCRIPTION  edge detection
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
#include "IM_Edge.h"

#include "IM_IO.h"
 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);


/*********************************************************************/
  
type_edge TypeEdge  =  DEF_EDGE_METHOD;
Bool KeepMaxima = False;
float Threshold = 0;
float ScaleGauss=-1;
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    edge_usage();
    
    manline();
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "             Select maxima or zero-crossing.\n");
    manline();

    fprintf(OUTMAN, "         [-t ThresholdValue]\n");
    fprintf(OUTMAN, "             Threshold the edge map.\n");
    manline();

    fprintf(OUTMAN, "         [-S Sigma]\n");
    fprintf(OUTMAN, "             Scale parameter. Used with DroG and LoG methods.\n");
    fprintf(OUTMAN, "             Default is 1/sqrt(3) \n");
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
static void infinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   
    /* get options */
    while ((c = GetOpt(argc,argv,"t:M:kS:vzZ:")) != -1) 
    {
	switch (c) 
        {
  	    case 'v': Verbose = True; break;        
            case 'M':
                if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "bad edge detection method: %s\n", OptArg);
	            usage(argv);
 		}
                if ((c <= 0) && (c >  NBR_EDGE))   
                {
		    fprintf(OUTMAN, "bad structural element type: %s\n", OptArg);
	            usage(argv);
 		}
		TypeEdge = (type_edge) (c-1);
 		break; 
 	    case 'k':
	        KeepMaxima = True;
		break;
	   case 't':
	        if (sscanf(OptArg,"%f",&Threshold) != 1) 
                {
		    fprintf(OUTMAN, "bad threshold value: %s\n", OptArg);
	            usage(argv);
 		}
		break;
	   case 'S':
	        if (sscanf(OptArg,"%f",&ScaleGauss) != 1) 
                {
		    fprintf(OUTMAN, "bad scale value: %s\n", OptArg);
	            usage(argv);
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
    if ((ScaleGauss > 0) &&
        (TypeEdge != ED_LOG ) && (TypeEdge != ED_CANNY))
    {
       fprintf(OUTMAN, "Error: -S option valid only with DroG and LoG method ... \n");
       exit(-1);
    }
       
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

/*********************************************************************/


int main(int argc, char *argv[])
{
    Ifloat Dat,Edge;
    fitsstruct Header;
    char Cmd[512];
    int k;
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);

    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Header.origin = Cmd;

    if (Verbose == True )
    {    
       cout << "Name File in  = " << Name_Imag_In << endl ;
       cout << "Name File out = " <<  Name_Imag_Out << endl ;
       cout << "  Nl = " << Dat.nl() << endl;
       cout << "  Nc = " << Dat.nc() << endl << endl;   
       cout << "  Edge detection method = " << StringEdge (TypeEdge) << endl;  
       if (Threshold > 0)
         cout << "  Threshold = " << Threshold << endl;
       if (ScaleGauss > 0)  
         cout << "  Scale parameter = " << ScaleGauss << endl;
       if (KeepMaxima == True)
          cout << "  Select local maxima in gradient direction "<< endl;
       cout << endl;
    }
 
    Edge.alloc(Dat.nl(),Dat.nc(), "Erosion");
    EDGE CED(TypeEdge);
    CED.SelectMaxima = KeepMaxima;
    if (ScaleGauss > 0)  CED.Sigma = ScaleGauss;
    CED.Bord = I_MIRROR;
    CED.detection(Dat, Edge);
    if (Threshold > 0) threshold(Edge, Threshold);
    io_write_ima_float (Name_Imag_Out, Edge, &Header);
    exit(0);
} 

