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
#include "IM3D_Radon.h"
#include "IM3D_IO.h"

char Name_Cube_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
type_radon3d RadonMethod = DEF_RADON3D;
Bool Reverse=False; // Reverse transform
 
/***************************************/
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    manline();
//     fprintf(OUTMAN, "         [-m type_of_radon_method]\n");
//     for (int i = 1; i <= NBR_RADON3D; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i,StringRadon3d((type_radon3d)(i-1)));
//     fprintf(OUTMAN, "              default is %s.\n", StringRadon3d(RadonMethod));
//     manline();    

    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse RADON transform.\n");  
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
    while ((c = GetOpt(argc,argv,"rvzZ")) != -1) 
    {
	switch (c) 
        { 
            case 'r': Reverse=True;break;
           case 'm': 
                if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of radon method: %s\n", OptArg);
                    exit(-1);
                    
                }
                if ((c > 0) && (c <= NBR_RADON3D)) RadonMethod = (type_radon3d) (c-1);
                else  
                {
                    fprintf(OUTMAN, "Error: bad type of radon method: %s\n", OptArg);
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
       

       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
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
    extern softinfo Soft;
    Soft.mr4();
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Cube_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        // cout << "Radon method = " <<  StringRadon3D(RadonMethod) << endl;  
        if (Reverse == True) cout << "Inverse Radon transform  " <<  endl;  
    }
 
    Radon3D RD;
    RD.Verbose = Verbose;
    if (Reverse == False)
    {
       fltarray Data;
       Ifloat Result;
       io_3d_read_data(Name_Cube_In, Data);
       int Nx = Data.nx();
       int Ny = Data.ny();
       int Nz = Data.nz();
       if ((Nx != Ny) || (Nx != Nz))
       {
           cout << "Error: the three dimensions must be identical ... " << endl;
           exit(-1);
       }

//         if ((is_power_of_2(Nx) == False) || (is_power_of_2(Ny) == False)|| (is_power_of_2(Nz) == False))
//         {
//            if (RadonMethod == RADON3D_FFT)
//            {
//               cout << "Error: image size must be power of two ... " << endl;
//               exit(-1);
//            }
//         }
        RD.alloc(Nx, Ny, Nz, RadonMethod);
        Result.alloc(RD.radon_nl(), RD.radon_nc(), "rad");
        RD.transform(Data,Result);
	io_write_ima_float(Name_Imag_Out, Result);
    }
    else
    {
        Ifloat Data;
	fltarray Result;
        io_read_ima_float(Name_Cube_In, Data);
	int Nx,Ny,Nz;
	Ny = Data.nc()/2;
	Nx = Nz = Ny;
        RD.alloc(Nx, Ny, Nz, RadonMethod);
        Result.alloc(Nx, Ny, Nz);
        RD.recons(Data, Result);
	io_3d_write_data(Name_Imag_Out, Result);
    }
    exit(0);
}
