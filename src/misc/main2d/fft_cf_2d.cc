


#include "IM_Obj.h"
#include "IM_IO.h"
#include "FFTN_2D.h"

char Name_Imag_In[100];
char Name_Imag_Out[100];
int Dir = 1;
Bool Verbose = False;
Bool CenterZeroFreq = True;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    manline();
    manline();
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse Fourier transform\n");
    manline();
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             Zero frequency is by default in the middle of the image.\n");
    fprintf(OUTMAN, "             If -s is set, zero frequency is at the bottom left.\n");
    manline();
    vm_usage();
    manline();
    
    verbose_usage();
    manline();
    manline();
    exit(-1);
}

/*********************************************************************/

static void infinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif

    /* get options */
    while ((c = GetOpt(argc,argv,"srvzZ:")) != -1) 
    {
        switch (c) 
        {
	   case 's': CenterZeroFreq = (CenterZeroFreq == True) ? False: True ; break;
	   case 'r': Dir = -1; break;
	   case 'v': Verbose = True; break;
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
            case '?': usage(argv);
       }
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
 
int main (int argc,char *argv[])
{
 
    lm_check(LIC_MR1);
    infinit(argc, argv);
    
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "Name File out = " << Name_Imag_Out << endl ;
       if (Dir == 1)
         cout << "Direct Fourier transform" << endl;
       else cout << "Inverse Fourier transform" << endl;
    }
    cfarray cData;
    fits_read_cfarr2d(Name_Imag_In, cData);

    Icomplex_f Pict;
    Pict.alloc(cData.buffer(), cData.ny(),cData.nx());
    
    // Fourier transform 
    FFTN_2D CF;
    Bool Reverse = (Dir == 1) ? False : True;
    CF.CenterZeroFreq = CenterZeroFreq;
    CF.fftn2d(Pict, Reverse);
    fits_write_cfarr2d (Name_Imag_Out, cData);
   // io_write_ima_complex_f (Name_Imag_Out, Pict);
    exit(0);
}
