/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/06/98
**    
**    File:  im_random.cc
**
*******************************************************************************
**
**    DESCRIPTION  generates a noise
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
#include "IM_Math.h"
#include "IM_IO.h"

 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
 
int Nl = 256;
int Nc = 256;
int J1=-1;
int J2=-1;

int TypeNoise = 1;
float Sigma = 1.;
Bool Verbose = False;
unsigned int InitRnd = 100;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options out_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-x Nc]\n");
    fprintf(OUTMAN, "             Number of columns.\n");
    fprintf(OUTMAN, "             Default is %d.\n",Nc);
    manline();
    fprintf(OUTMAN, "         [-y Nl]\n");
    fprintf(OUTMAN, "             Number of lines.\n");
    fprintf(OUTMAN, "             Default is %d.\n",Nl);
    manline();
    fprintf(OUTMAN, "         [-t DistribLaw]\n");
    fprintf(OUTMAN, "             1: Gaussian.\n");
    fprintf(OUTMAN, "             2: Poisson.\n");
    fprintf(OUTMAN, "             3: Rayleigh.\n");
    fprintf(OUTMAN, "             4: Laplace.\n");
    fprintf(OUTMAN, "             Default is Gaussian.\n");
    manline();
    fprintf(OUTMAN, "         [-s DistribParam]\n");
    fprintf(OUTMAN, "             Distribution paramter.\n");
    fprintf(OUTMAN, "             Default is 1.\n");
    manline();
    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "             default is 100. \n\n");    manline();
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
    while ((c = GetOpt(argc,argv,"x:y:s:t:vzZ:I:")) != -1) 
    {
	switch (c) 
        {
	     case 'I':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		    exit(-1);
		}
                InitRnd = (unsigned int) c;
		break;
	    case 'v': Verbose = True; break;        
            case 'x':
	        if (sscanf(OptArg,"%d", &Nc) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameter: %s\n", OptArg);
		    exit(-1);
		}
 		break; 
            case 'y':
                if (sscanf(OptArg,"%d", &Nl) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break;           
 	   case 't':
                if (sscanf(OptArg,"%d", &TypeNoise) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad distribution law parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break;   
	   case 's':
	        if (sscanf(OptArg,"%f", &Sigma) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad distribution parameter: %s\n", OptArg);
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
	    default: usage(argv);
	}
    }

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
    Ifloat Dat;
    int i;
   
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    infinit(argc, argv);

    Dat.alloc(Nl,Nc,"noise");
    if (Verbose == True )
    {
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Nl = " << Dat.nl() << endl;
       cout << "  Nc = " << Dat.nc() << endl;     
       if (TypeNoise==1)
         cout << "  Gaussian noise " << endl;
       else if (TypeNoise==2)
         cout << "  Poisson noise " << endl;
       else if (TypeNoise==3)
         cout << "  Rayleigh  noise " << endl; 
       else  cout << "  Laplace  noise " << endl;
    }
    switch(TypeNoise)
    {
       case 1: im_noise_gaussian (Dat, Sigma, InitRnd); break;
       case 2: Dat.init(Sigma);
               im_noise_poisson (Dat, 1., InitRnd);
	       break;
       case 3: im_noise_rayleigh(Dat, 1, InitRnd);
               for (i=0;i<Nl*Nc;i++) Dat(i)*=Sigma;
	       break;
       case 4: im_noise_laplace(Dat, 1, InitRnd);
               for (i=0;i<Nl*Nc;i++) Dat(i)*=Sigma;
	       break;
       default: cerr <<  "Error: bad distribution law parameters ... " << endl;
                exit(-1);
		break;
    }
    
    io_write_ima_float(Name_Imag_Out, Dat);
    exit(0);
} 

