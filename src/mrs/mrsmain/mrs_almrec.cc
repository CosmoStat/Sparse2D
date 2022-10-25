/******************************************************************************
**                   Copyright (C) 2011 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  02/11/2011
**    
**    File:  mrs_almrec.cc
**
*******************************************************************************
**
**    DESCRIPTION  genus program
**    ----------- 
**                 
******************************************************************************/


#include"HealpixClass.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;

Bool NormALM = False;
float  ZeroPadding=0.0;
int Lmax = 0.;
Bool OptS_PS = False;
Bool CFIMA= False;
int Nside=16;
int NbrIter=0;

/***************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-n Nside]\n");
    fprintf(OUTMAN, "             Default is 16.\n");
        
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");
    
	fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration for the reconstruction. Default is no iteraton\n");
    
    fprintf(OUTMAN, "         [-T]\n");
    fprintf(OUTMAN, "             Read the input as a complex image.\n");
    
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    float Val; 
    int seed;

    while ((c = GetOpt(argc,argv,(char *) "n:TNi:l:p:vzZ")) != -1) 
    {
        switch (c) 
        { 
		    case 'T':
                CFIMA=True;
                break;
            case 'N':
                NormALM=True;
                break;
            case 'p':
                if (sscanf(OptArg,"%f",&ZeroPadding) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad zero padding  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
            case 'i':
                if (sscanf(OptArg,"%d",&NbrIter) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad nside  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
            case 'n':
                if (sscanf(OptArg,"%d",&Nside) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad nside  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
            case 'l':
                if (sscanf(OptArg,"%d",&Lmax) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
            case 'v': Verbose = True; break;
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
  
/*********************************************************************/
/* void mrs_almtrans(Hdmap  map, CAlmR & A)
{
	map.info();
    double avg=map.average();
    map.add(-avg);
    if (map.Scheme()==NEST) map.swap_scheme();
    
    map2alm_iter(map, A, A.Niter, A.weight_T);
   
   // template<typename T> void map2alm_iter (const Healpix_Map<T> &map,
   // Alm<xcomplex<T> > &alm, int num_iter, const arr<double> &weight);
  
   int Niter = A.Niter;
    Alm<xcomplex<double> > A1(A.Lmax(),A.Mmax());
    map2alm_iter(map, A1, Niter, A.weight_T);
    A(0,0) += avg*sqrt(fourpi);
    map.info();	
}
*/
/*********************************************************************/


int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (NbrIter > 0) cout << "# Nbr of iter = " << NbrIter << endl;
    }
//   /Users/Shared/software/Healpix_2.01/src/cxx/Healpix_cxx    
//  alm_map_tools.cc et .h

    //  cout << "TEST " << Name_Imag_In  << endl;
   // Healpix_Map<float> map;
   Hdmap map;
   Hdmap Result;
   CAlmR  ALM;
   
    if (CFIMA == False)
    {
        ALM.read(Name_Imag_In, Nside);
        if (Verbose == True) cout << "Lmax = " << ALM.Lmax() << endl;
    }
    else
    {
        dblarray A;
        fits_read_dblarr(Name_Imag_In, A);
        Lmax = A.nx() - 1;
        int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);

        ALM.alloc(Nside, Lmax_Tot);
        for (int l=0; l <= ALM.Lmax(); l++)
        for (int m=0; m <= l; m++) 
        {
            ALM(l,m) = xcomplex<REAL>( A(l,m,0),A(l,m,1));
        }
   }
   if (NbrIter > 0) ALM.Niter = NbrIter;
   ALM.Norm = NormALM;
   ALM.alm_rec(Result);
    if (Verbose == True)
    {   
      double Min, Max;
      Result.minmax(Min,Max);
      cout << " Npix = " << Result.Npix() << endl;
      cout << " Min = " << Min << " Max = " << Max << endl;
   }
   Result.write(Name_Imag_Out);
   exit(0);
}

