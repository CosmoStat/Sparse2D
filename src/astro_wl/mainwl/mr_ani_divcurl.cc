/******************************************************************************
**                   Copyright (C) 2007 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Erwan Deriaz & Jean-Luc Starck
**
**    Date:  08/09/08
**    
**    File:  mr_ani_divcurl.cc
**
*******************************************************************************
**
**    DESCRIPTION  anisotropic div and curl decomposition program  
**    ----------- 
**                 
******************************************************************************/

#include "MR_DivCurl.h"

char Name_Imag_In[256]; /* input data file   */
char Name_Imag_Out[256]; /* output data file */
 
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;   // Verbose mode
int NbrScale=0;       // Number of scales, default is automatically estimated

Bool CurlTrans = False;             // If true, a rotational-free transform is applied instead of the div-free transform
Bool Reverse = False;               // True if we apply a reconstruction. By default, it is a transformation
Bool WriteComp = False;             // if true, we write also the two wavelet transforms
Bool V0Proj = True;
int Debug = 0;
Bool WaveletBorder = False;

/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
        
    //fprintf(OUTMAN, "         [-n NbrScale]\n");
    //fprintf(OUTMAN, "             Number of wavelet scales. Default is %d\n", NbrScale);
    //manline();
	
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse transform.\n");
    fprintf(OUTMAN, "         [-c]\n");
    fprintf(OUTMAN, "             Curl transform.\n");
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Do not apply the V0 projection or the Point Value interpolation.\n");	
//	fprintf(OUTMAN, "         [-w]\n");
//    fprintf(OUTMAN, "             write the two reconstructed components:\n");
//    fprintf(OUTMAN, "             File names are  xx_div for the divergence free decomposition.\n ");
//    fprintf(OUTMAN, "                             xx_curl for curl free decomposition.\n ");
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv, (char *) "bdPn:wru:cx:y:vzZ")) != -1) 
    {
	switch (c) 
        {
        case 'b': WaveletBorder = (WaveletBorder == True) ? False: True; break;
        case 'd': Debug = True; break;
		case 'P': V0Proj = (V0Proj == True) ? False: True; break;
	    case 'w': WriteComp = True; break;
 	    case 'c': CurlTrans = (CurlTrans == True) ? False: True; break;
 	    case 'r': Reverse = True; break;
	    case 'n': 
                if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
 		break;
          case 'v': Verbose = True;break;
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
}
 
/***************************************/
 
int main(int argc, char *argv[])
{
    int k;
    
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
		if (CurlTrans == True)  cout << "CURL-FREE DECOMPOSITION" << endl;
 	    else cout << "DIV-FREE DECOMPOSITION" << endl; 
    }
 
    fltarray Data;
    fltarray Result;
    
    fits_read_fltarr(Name_Imag_In, Data);
 	
    // float Min,Max;
    int Nx = Data.nx();
    int Ny = Data.ny();
    // float *Ptr = Data.buffer();  // pointer to the data
    if (Verbose == True)
    { 
       cout << "Data size: Nx = " << Nx <<  " Ny = " << Ny  <<  endl; 
       Data.info("Input stat: ");
       // cout << "Min = " << Data.min()  << "  Max = " << Data.max() <<  " Sigma = " << Data.sigma() << endl;
    }
    fltarray TabRes;
    fltarray Tabd;
    
    fltarray D1,D2;
    D1.alloc(Data.nx(), Data.ny()-1);
    D2.alloc(Data.nx()-1, Data.ny());
    for (int i=0; i < D1.nx(); i++)
    for (int j=0; j < D1.ny(); j++) D1(i,j) = Data(i,j,0);
    for (int i=0; i < D2.nx(); i++)
    for (int j=0; j < D2.ny(); j++) D2(i,j) = Data(i,j,1);
    MR_ANI_DIVCURL MRD1;
    MRD1.alloc(Ny, Nx, NbrScale);
    
/*    fltarray R1, R2, T1, T2;
    cout << "D1: " << D1.nx() << " " << D1.ny() << endl;
    cout << "D2: " << D2.nx() << " " << D2.ny() << endl;

    MRD1.curl_bord_transform (D1, D2, T1, T2);
    cout << "T1: " << T1.nx() << " " << T1.ny() << endl;
    cout << "T2: " << T2.nx() << " " << T2.ny() << endl;

    MRD1.curl_bord_recons (R1, R2, T1, T2);
    
    cout << "R1: " << R1.nx() << " " << R1.ny() << endl;
    cout << "R2: " << R2.nx() << " " << R2.ny() << endl;
    
    D1.info("D1");
    D2.info("D2");
    
    D1 -= R1;
    D2 -= R2;
    D1.info("D1");
    D2.info("D2");
    exit(0);*/


    MR_ANI_DIVCURL MRD;
    MRD.Verbose = Verbose;
	// cout << " ALLOC  " << endl;	
    MRD.alloc(Ny, Nx, NbrScale);
	MRD.V0_Proj = V0Proj;
	MRD.PointValueRec = V0Proj;
    MRD.BorderWavelet = WaveletBorder;
    Bool DivTrans = (CurlTrans == True) ? False: True;
	Bool DebugTest = False;
    MRD.Debug=Debug;
	// cout << " RUN " << endl;
	
	if (Reverse == False) // Transformation
	{
		if (Verbose == True) cout << "Transform .... " << endl;
	    MRD.transform(Data, DivTrans);
 	    
	    // if (Verbose == True) MRD.info();
	    if (Verbose == True) cout << "Write the result .... " << endl;

        
        fits_write_fltarr(Name_Imag_Out, MRD.RotDiv_Trans);
	    
        fltarray Rec;
        Rec = Data;
        Rec.init();
        MRD.recons(Rec, DivTrans);
        Data -= Rec;
        Data.info("REC");
        fits_write_fltarr("xx_tt.fits", Rec);
        fits_write_fltarr("xx_resi.fits", Data);

 	    // if (WriteComp == True)
	    // {
	    //   Data.init();
	    //   MRD.component_recons(Data);
	    //   if (DivTrans == True) fits_write_fltarr((char *) "xx_div", Data);
	   //    else fits_write_fltarr((char *) "xx_curl", Data);
	   //  }
  	}
	else // here, it is the reconstruction
	{
	   if (Verbose == True) cout << "Reconstruction .... " << endl;
	   // MRD.RotDiv_Trans.alloc(Data.nx(), Data.ny(), Data.nz());
	   MRD.RotDiv_Trans = Data;
	   MRD.recons(Data, DivTrans);
	   if (Verbose == True) cout << "Write the result .... " << endl;
	   fits_write_fltarr(Name_Imag_Out, Data);
	}
 	// if (DebugTest == True) TabRes.info("Output stat: ");
   
	
    exit(0);
}
