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
**    File:  im3d_put.cc
**
*******************************************************************************
**
**    DESCRIPTION  insert a sub-image
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
#include "IM3D_IO.h"

 
char Name_Imag_In[100];
char Name_Imag_Out[100];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
 
int I1=-1;
int J1=-1;
int K1=-1;
Bool Verbose = False;
 
/*********************************************************************/

static void usage(char *argv[]) {
    
    fprintf(OUTMAN, "\n\nUsage: %s options in_cube inout_cube \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-x Xstart]\n");
    fprintf(OUTMAN, "             x-coordinate of the image to insert.\n");
    manline();
    
    fprintf(OUTMAN, "         [-y Ystart]\n");
    fprintf(OUTMAN, "             y-coordinate of the image to insert.\n");
    manline();
        
    fprintf(OUTMAN, "         [-z Zstart]\n");
    fprintf(OUTMAN, "             z-coordinate of the image to insert.\n");
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
            case 'x':
	        if (sscanf(OptArg,"%d", &I1) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad column parameter: %s\n", OptArg);
		    exit(-1);
		}
 		break; 
            case 'y':
                if (sscanf(OptArg,"%d", &J1) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameter: %s\n", OptArg);
		    exit(-1);
		}
 		break;  
            case 'z':
                if (sscanf(OptArg,"%d", &K1) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad line parameter: %s\n", OptArg);
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


int main(int argc, char *argv[]) {
    
    
    fltarray Dat1,Dat2;
    fitsstruct Header;
    char Cmd[512];
    extern softinfo Soft;

    Soft.mr3();	
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    infinit(argc, argv);
  
    io_3d_read_data(Name_Imag_In, Dat1);
    io_3d_read_data(Name_Imag_Out, Dat2, &Header);
    Header.origin = Cmd;
   
//     if (     ((Dat1.nx() ==0) || (Dat1.nx() != Dat2.nx())) 
//          || ((Dat1.ny() ==0) || (Dat1.ny() != Dat2.ny()))
//          || ((Dat1.nz() ==0) || (Dat1.nz() != Dat2.nz()))) {
//        cerr << "Error: cubes must have the same size ... " << endl;
//        exit(-1);
//     }
    
    if (I1 < 0) I1 = 0;  
    if (J1 < 0) J1 = 0;
    if (K1 < 0) K1 = 0;
    if (I1 >= Dat2.nx()) I1 = Dat2.nx() - 1;
    if (J1 >= Dat2.ny()) J1 = Dat2.ny() - 1;
    if (K1 >= Dat2.nz()) K1 = Dat2.nz() - 1;   
    int Nx = Dat1.nx();
    int Ny = Dat1.ny();
    int Nz = Dat1.nz();
    
    if ((I1+Nx > Dat2.nx()) || (J1+Ny > Dat2.ny()) || (K1+Nz > Dat2.nz()))  
    {
       cerr << "Error: first cube cannot be inserted in the second one ... " << endl;
       cerr << "       First image: Nx = " << Nx 
                               << " Ny = " << Ny 
			       << " Nz = " << Nz << endl;
       cerr << "       Second image: Nx = " << Dat2.nx()  
                                << " Ny = " << Dat2.ny()
				<< " Nz = " << Dat2.nz() << endl;
       cerr << "       Insertion position : x = " << I1  
                                       << " y = " << J1  
				       << " z = " << K1 << endl;
       exit(-1);
    }
    if (Verbose == True )
    {
       cout << "Name File in = " << Name_Imag_In << endl ;
       cout << "  Nx = " << Dat1.nx() << endl;      
       cout << "  Ny = " << Dat1.ny() << endl;
       cout << "  Nz = " << Dat1.nz() << endl;    
       cout << "Name File in-out 2  = " << Name_Imag_Out << endl ;
       cout << "  Pos x = " << I1 
             << " Pos y = " << J1 
	     << " Pos z = " << K1 << endl;
       cout << "  Nx = " << Nx << endl;
       cout << "  Ny = " << Ny << endl;      
       cout << "  Nz = " << Nz << endl;  
       cout << endl;
    }              
    for (int i=0;i<Nx;i++)
       for (int j=0;j<Ny;j++) 
          for (int k=0;k<Nz;k++) 
	     Dat2(i+I1,j+J1,k+K1) = Dat1(i,j,k);
    
    io_3d_write_data (Name_Imag_Out, Dat2,&Header);
    exit(0);
} 

