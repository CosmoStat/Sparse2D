/******************************************************************************
 **                   Copyright (C) 2007 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Jean-Luc Starck - Simon Beckouche - Francois Lanusse
 **
 **    Date:  08 Nov 2013
 **
 **    File:  dl_omp.cc
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Orthogonal Matching Pursuit for Dictionary Learning
 **    -----------
 **
 **    Usage: Do not eat
 **
 ******************************************************************************/

#include "DefMath.h"
#include "MatrixOper.h"  // defined in $TOOLS
#include <cmath>
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"
#include "C_OMP.h"
/****************************************************************************/

char Name_Dico_In[256];  // 2D array
char Name_Vect_In[256];  // 2D array
char Name_Coef[256];      //
int defaultSparsityTarget = 5;
int SparsityTarget = defaultSparsityTarget;
bool UseSparsity = True;
double ErrorTarget = -1.;
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose=False;

/*********************************************************************/

static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options input_dico input_vect output_coef \n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-S TargetSparsity]\n");
	fprintf(OUTMAN, "         [-E]\n");
	fprintf(OUTMAN, "                 Use error reconstruction instead of sparsity target.\n");
	verbose_usage();
	vm_usage();
	manline();
	exit(-1);
}

/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
	int c;
#ifdef LARGE_BUFF
	int VMSSize=-1;
	Bool OptZ = False;
	char VMSName[1024] = "";
#endif

	/* get options */
	while ((c = GetOpt(argc,argv,"S:E:vzZ:")) != -1)
	{
		switch (c)
		{
		case 'S':
			/* -S <SparsityTarget> */
			if (sscanf(OptArg,"%d", &SparsityTarget) != 1)
			{
				fprintf(OUTMAN, "bad value: %s\n", OptArg);
				exit(-1);
			}
			break;
		case 'E':
			sscanf(OptArg,"%lf", &ErrorTarget);
			UseSparsity = false; break;
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
		case '?':
			usage(argv);
		}
	}



	/* get optional input file names from trailing
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Dico_In, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(Name_Vect_In, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(Name_Coef, argv[OptInd++]);
	else usage(argv);


	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
	if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

int main(int argc, char *argv[])
{
	/* Get command line arguments, open input file(s) if necessary */
	fitsstruct Header;
	char Cmd[512];
	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

	transinit(argc, argv);

	if (Verbose == True)
	{
		cout << "Dico Filename in = " << Name_Dico_In << endl;
		cout << "Data Filename in = " << Name_Vect_In  << endl;
		cout << "Coef Filename out = " << Name_Coef  << endl;
		if (UseSparsity  == True)
			if (SparsityTarget  != 0) cout << "Use Sparsity target  =  " << SparsityTarget << endl;
			else cout << " Use default sparsity target  =  " << defaultSparsityTarget << endl;
		else cout << "Use Error target " << endl;

	}

	if (Verbose == True) cout << "\nReading the data"<< endl;


	dblarray input_Dico;

	fits_read_dblarr (Name_Dico_In, input_Dico);
	if (Verbose == true)
		cout << "Dictionary: " << input_Dico.ny() << " atoms of length " << input_Dico.nx() << endl;

	dblarray Vect;
	fits_read_dblarr (Name_Vect_In, Vect);
	cout << "Training set: " << Vect.ny() << " samples of length " << Vect.nx() << endl;

	/*fltarray flVect;
	fitsstruct FitsHeader;
	io_1d_read_data(Name_Vect_In, flVect, &FitsHeader);
	reform_to_1d(flVect);
	for(int i=0; i < flVect.nx(); i++) Vect(i) = flVect(i);*/

	dblarray Coef;
	if (UseSparsity == false) SparsityTarget = input_Dico.nx();
	else ErrorTarget = -1;
	C_OMP coder(input_Dico);
	Coef = coder.omp(Vect,SparsityTarget,ErrorTarget,Verbose);
	/*cout << "Dico.nx = " << Dico.nx()<<" Dico.ny = "<< Dico.ny()<< endl;
	cout << "Vect.nx = " << Vect.nx()<<" Vect.ny = "<< Vect.ny()<< endl;
	cout << "Coef.nx = " << Coef.nx()<<" Coef.ny = "<< Coef.ny()<< endl;*/


	fits_write_dblarr (Name_Coef, Coef);
	exit(0);
}




//    C_OMP OMP;
//    cout << " Def Sparse Value = " << OMP.SparsityTarget << endl;
//    if (SparsityTarget  != 0)    OMP.SparsityTarget  = SparsityTarget;
//    cout << " Def Sparse Value = " << OMP.SparsityNumber << endl;
//    if (UseSparsity  == false) OMP.UseSparsityNumber = false;
//
//    fits_read_dblarr (Name_Dico_In, OMP.TabDico);
//    int Nx = OMP.TabDico.nx();
//    int Ny = OMP.TabDico.ny();
//    if (Verbose == True) cout << "Nx = " << Nx << " Ny = " << Ny << endl;
//
//    fits_read_dblarr (Name_Vect_In, OMP.TabVect);
//    Nx = OMP.TabVect.nx();
//    if (Verbose == True) cout << "Nx = " << Nx << endl;
//
//    OMP.alloc();
//
//    {
//        double Mean = OMP.TabDico.mean();
//        double Sigma = OMP.TabDico.sigma();
//        for (int i=0;i<Nx;i++)
//        for (int j=0;j<Ny;j++) OMP.TabDico(i,j) = (OMP.TabDico(i,j)-Mean)/Sigma;
//        OMP.TabDico.info();
//    }





/*
   // dblarray TabDico;
  //   dblarray TabVect;

       io_3d_read_data(Name_Cube_In, Dat, &Header);

       int Nx = Dat.nx();
       int Ny = Dat.ny();
       int Nz = Dat.nz();
       if (Verbose == True) cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;

       if (Normalize == True)
       {
         double Mean = Dat.mean();
         double Sigma = Dat.sigma();
         for (int i=0;i<Nx;i++)
         for (int j=0;j<Ny;j++)
         for (int k=0;k<Nz;k++) Dat(i,j,k) = (Dat(i,j,k)-Mean)/Sigma;
       }


       // MR2D1D WT;
       GMCA_2D WT;

       if (Verbose == True) cout << "Alloc ...  " << endl;
       WT.alloc(Nx, Ny, Nz, Transform, NbrScale2d, 1); // On ne regularise que la MixingMat

       if (Verbose == True)	cout << "2d1d_trans ... "<< endl;

       //WT.transform (Dat);  // Pas utile

       // WT.write(Name_Out);
       // Compute the 2D1D transform
       fltarray TabVect;
       WT.transform_to_vectarray(Dat, TabVect);
       // fits_write_fltarr ("xx_tabvect.fits", TabVect);

       // Initalize the class for GMCA


       if (UseMask == True)
       {
          fits_read_fltarr (Name_Mask, WT.Mask);
       }

        fits_write_fltarr(Name_Out, EstSources);

 */
