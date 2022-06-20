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
**    File:  im_rot.cc
**
*******************************************************************************
**
**    DESCRIPTION  rotation-translation program
**    ----------- 
**                 
******************************************************************************/
   
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Rot.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float AnglePar=0.;
float Dx=0.;
float Dy=0.;
Bool OptReversible=False;
Bool OptX=False;
Bool OptY=False;
Bool OptA=False;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    manline();
    fprintf(OUTMAN, "         [-r RotAngle]\n");
    fprintf(OUTMAN, "             Rotation angle in degrees\n");

    manline();
    fprintf(OUTMAN, "         [-x Dx] \n");
    fprintf(OUTMAN, "             x-axis translation offset\n");

    manline();
    fprintf(OUTMAN, "         [-y Dy] \n");
    fprintf(OUTMAN, "             y-axis translation offset\n");
    manline(); 

    fprintf(OUTMAN, "         [-f] \n");
    fprintf(OUTMAN, "             Use the Fourier transform in order to\n");
    fprintf(OUTMAN, "             have a reversible transformation.\n");
    // fprintf(OUTMAN, "             Image size must be a power of 2.\n");

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
    while ((c = GetOpt(argc,argv,"x:y:r:fvzZ")) != -1) 
    {
	switch (c) 
        { 
           case 'f':  OptReversible = True;break;
           case 'r': 
                if (sscanf(OptArg,"%f",&AnglePar) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad angle: %s\n", OptArg);
                    exit(-1);
                    
                }
                OptA=True;
                break;
           case 'x': 
                if (sscanf(OptArg,"%f",&Dx ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad x-axis offset: %s\n", OptArg);
                    exit(-1);
                    
                }
                OptX=True;
                 break;
           case 'y': 
                if (sscanf(OptArg,"%f",&Dy ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad y-axis offset: %s\n", OptArg);
                    exit(-1);
                    
                }
                OptY=True;
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


/***************************************************************************/

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

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        if (OptA == True) cout << "Angle = " << AnglePar  << endl; 
        if (OptReversible == True) cout << "Reversible operator " << endl; 
        if (OptX == True) cout << "Dx = " << Dx  << endl;  
        if (OptY == True) cout << "Dy = " << Dy  << endl;  
    }

    Rotation ROT;
    Ifloat Data, DataRot, DataRec;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;

    ROT.Verbose = Verbose;
    if (OptReversible == False)
                 ROT.im_move_rotate(Data,DataRot,AnglePar,Dx,Dy);
    else
    {
//        if ((is_power_of_2(Data.nl()) == False) || (is_power_of_2(Data.nc()) == False))
//        {
//            cout << "Error: image size must be power of two ... " << endl;
//            exit(-1);
//        }
       if (Verbose == True) cout << "Translation ... " << endl;
       if ((OptX == True) || (OptY == True)) translate2d(Data, Dx, Dy);
       if (ABS(AnglePar) > FLOAT_EPSILON) 
       {
          if (Verbose == True) cout << "Rotation ... " << endl;
          ROT.im_rotate(Data, DataRot,AnglePar);
       }
       else DataRot = Data;
    }
    io_write_ima_float (Name_Imag_Out, DataRot, &Header);
    exit(0);
}
