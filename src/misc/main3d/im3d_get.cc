/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck, Philippe Querre
**
**    Date:  28/02/01
**    
**    File:  im3d_get.cc
**
*******************************************************************************
**
**    DESCRIPTION  extract a sub-image
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
#include "IM1D_IO.h"
#include "IM_IO.h"
#include "IM3D_IO.h"

 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
 
int I1=-1;
int I2=-1;
int J1=-1;
int J2=-1;
int K1=-1;
int K2=-1;
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[]) {
    fprintf(OUTMAN, "\n\nUsage: %s options in_cube out_cube \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");



    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-x Xstart:XEnd]\n");
    fprintf(OUTMAN, "             First and last columns to extract.\n");
    manline();
    fprintf(OUTMAN, "         [-y Ystart:YEnd]\n");
    fprintf(OUTMAN, "             First and last lines to extract.\n");
    manline();
    fprintf(OUTMAN, "         [-z Zstart:ZEnd]\n");
    fprintf(OUTMAN, "             First and last lines to extract.\n");
    manline();   
    
    /*vm_usage();
    manline();*/
    
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
    while ((c = GetOpt(argc,argv,"x:y:z:v")) != -1) 
    {
	switch (c) 
        {
 	    case 'v': Verbose = True; break;
	    case 'z':
	        if (sscanf(OptArg,"%d:%d", &K1,&K2) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break; 	             
            case 'x':
	        if (sscanf(OptArg,"%d:%d", &I1,&I2) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break; 
            case 'y':
                if (sscanf(OptArg,"%d:%d", &J1,&J2) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameters: %s\n", OptArg);
		    exit(-1);
		}
 		break;           
#ifdef LARGE_BUFF
/*	    case 'z':
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
                break;*/
#endif
 	    case '?': usage(argv); break;
	    default: usage(argv); break;
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


int main(int argc, char *argv[])
{
    fltarray Dat1,Dat2;
    fitsstruct Header;
    extern softinfo Soft;

    Soft.mr3();
  
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    infinit(argc, argv);

    io_3d_read_data(Name_Imag_In, Dat1, &Header);
    
    if ((Dat1.nx() ==0) || (Dat1.ny() ==0) ||(Dat1.nz() ==0)) {
      cerr << "Error: bad input dimension ... " << endl;
      cerr << "       only 3d data " << endl;
      exit(-1);   
    }
     
    if (I1 < 0) I1 = 0; 
    if ((I2 < 0) || (I2 >= Dat1.nx())) I2 = Dat1.nx() - 1;
      
    if (J1 < 0) J1 = 0;
    if ((J2 < 0) ||(J2 >= Dat1.ny())) J2 = Dat1.ny() - 1;
    
    if (K1 < 0) K1 = 0;
    if ((K2 < 0) ||(K2 >= Dat1.nz())) K2 = Dat1.nz() - 1;    
    
    int Nx = I2 - I1 + 1;
    int Ny = J2 - J1 + 1;
    int Nz = K2 - K1 + 1;
        
    Dat2.alloc(Nx,Ny,Nz);
    
    if (Verbose == True ) {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "  Nx = " << Dat1.nx() << endl;
       cout << "  Ny = " << Dat1.ny() << endl;
       cout << "  Nz = " << Dat1.nz() << endl;    
       cout << "Name File out = " << Name_Imag_Out << endl ;
       cout << "  Nx = " << Dat2.nx() << endl;
       cout << "  Ny = " << Dat2.ny() << endl;
       cout << "  Nz = " << Dat2.nz() << endl;        
       cout << endl;
    }       
    
    for (int i=0;i<Nx;i++)
       for (int j=0;j<Ny;j++)
          for (int k=0;k<Nz;k++)
	     Dat2(i,j,k) = Dat1(i+I1,j+J1,k+K1);
	     
    if (Nx==1)
       if      (Ny==1) {         Dat2.reform(Nz);    
                                 io_1d_write_data (Name_Imag_Out, Dat2);}
       else if (Nz==1) {         Dat2.reform(Ny);    
                                 io_1d_write_data (Name_Imag_Out, Dat2);}
       else            {         Dat2.reform(Ny,Nz);
                                 fits_write_fltarr(Name_Imag_Out, Dat2);} 
    else
       if ((Ny==1) && (Nz==1)) { Dat2.reform(Nx);    
                                 io_1d_write_data (Name_Imag_Out, Dat2);}
       else if (Ny==1)         { Dat2.reform(Nx,Nz);
                                 fits_write_fltarr(Name_Imag_Out, Dat2);}
       else if (Nz==1)         { Dat2.reform(Nx,Ny); 
                                 fits_write_fltarr(Name_Imag_Out, Dat2);}
       else                    { io_3d_write_data (Name_Imag_Out, Dat2, &Header);}
    
    exit(0);
} 

 
