/******************************************************************************
 **                   Copyright (C) 2007 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Jean-Luc Starck - Jerome Bobin
 **
 **    Date:  11/30/07
 **
 **    File:  mr_gmca.cc
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Generalized Morphological Component Analysis
 **    -----------
 **
 **    Usage: mr_gmca options cube output
 **
 ******************************************************************************/
#include "DefMath.h"
#include "MatrixOper.h"  // defined in $TOOLS
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_Block.h"
#include "C_OMP.h"
#include "C_DL1D.h"
#include "C_DENOISE_IMAGE.h"
Bool Ov = True;
bool UseSparsity=true;
/****************************************************************************/
char Name_Dico[256];
char Name_Image_in[256];
char Name_Data_out[256];
double ErrorTarget = -1.;
extern int  OptInd;
extern char *OptArg;
dblarray Dico;
dblarray Image_out;
dblarray Image_in;
dblarray stacked_patches;
int defaultSparsityTarget = 5;
int SparsityTarget = defaultSparsityTarget;
int OverlapNumber = 1;
extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose=False;
int Datanx,Datany,Datanz,Diconx,Dicony,Diconz,indData,indDico;
dblarray oldDico,oldData,oldD;
/****************************************************************************/
static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options noisy_Image Dictionary  denoised_Image \n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-S TargetSparsity]\n");
	fprintf(OUTMAN, "         [-E ErrorTarget, TargetSparsity ignored]\n");
	fprintf(OUTMAN, "         [-Q Overlaping number, default 1]\n");
	verbose_usage();
	vm_usage();
	manline();
	exit(-1);
}
/****************************************************************************/
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
	while ((c = GetOpt(argc,argv,"OS:E:Q:vzZ:")) != -1)
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
			if (ErrorTarget > 0) UseSparsity = false;
			break;
		case 'Q':
			sscanf(OptArg,"%d", &OverlapNumber);
			break;
		case 'O': Ov = (Ov ==True) ? False: True; break;
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
	if (OptInd < argc) strcpy(Name_Dico, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Image_in, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Data_out, argv[OptInd++]);
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

/****************************************************************************/
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
		cout << "Dictionary Filename = " << Name_Dico << endl;
		cout << "Image Filename in = " << Name_Image_in  << endl;
		cout << "Image Filename out = " << Name_Data_out  << endl;
		if (UseSparsity  == True) cout << "Use Sparsity target = " << SparsityTarget << endl;
		else cout <<"Use Error target " << ErrorTarget << endl;
	}
	if (Verbose == True) cout << "Reading the data"<< endl;
	dblarray Dico, Image_in, TabPatch, SigBlock;
	fits_read_dblarr (Name_Dico,Dico);
	int Nx = Dico.nx();
	int Ny = Dico.ny();
	int PatchSize = sqrt( double( Nx) );
	if (Verbose == True)
		if (Dico.naxis() == 2)
			cout << "Dictionary: " << Nx << " atoms of length " << Ny << endl;
		else cout << "Dictionary: Nx = " << Nx << " Ny = " << Ny << " Nz = " << Dico.nz() << endl;
	fits_read_dblarr (Name_Image_in, Image_in);
	Nx = Image_in.nx();
	Ny = Image_in.ny();
	int Nz = Image_in.nz();
	if (Verbose == True)
		cout << "Noisy image: Nx = " << Image_in.nx() << " Ny = " << Image_in.ny() << endl;


/*	// Determining if dictionary is 1d or 2d. If 2d, converting it into 1d by stacking atoms column after column
	is2dDico = (Dico.naxis() > 2);
	if (is2dDico == true)
	{
		if (Verbose == true)
			cout << "Converting initial 2d dictionary  into 1d dictionary" << endl;
		Diconx = Dico.nx(); // number of lines in each atom
		Dicony = Dico.ny(); // number of column in each atom
		Diconz = Dico.nz(); // number of atoms
		oldDico = Dico;
		Dico.free();
		Dico.alloc(Diconx*Dicony,Diconz);
		for (int k=0;k<Diconz;k++)
		{
			indDico = 0;
			for (int j=0;j<Dicony;j++)
				for (int i=0;i<Diconx;i++)
				{
					Dico(indDico,k) = oldDico(i,j,k);
					indDico++;
				}
		}
	}
	*/
	if (UseSparsity == false) SparsityTarget = Dico.nx();
	else ErrorTarget = -1;
	/****************************************************************************/
	// Calling main function
	C_DENOISE_IMAGE denoiser(Image_in,Dico);
	stacked_patches = denoiser.extract_patches(Image_in,OverlapNumber);
	if (Verbose == true)
		cout << stacked_patches.ny() << " noisy patches extracted stacked with an overlap of " << OverlapNumber << endl;
	intarray full_image_size = denoiser.compute_full_image_size(Image_in,OverlapNumber);
	Image_out = denoiser.denoise_image(stacked_patches,Dico,OverlapNumber,SparsityTarget,ErrorTarget,full_image_size,Verbose);
	fits_write_dblarr (Name_Data_out, Image_out);
	exit(0);
}
/****************************************************************************/
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
//
///****************************************************************************/
//
//void C_OMP::alloc()
//{
//    Na = TabDico.ny();
//    Npix = TabDico.nx();
//    if (Npix != TabVect.nx())
//    {
//       cout << "Error: vector and dictonary dimensions do not match. " << endl;
//       cout << "   Dico.nx = " << Npix << ", Vect.nx = " <<TabVect.nx() << endl;
//    }
//    TabCoef.alloc(Na);
//}
//
///****************************************************************************/
//
//void C_OMP::omp()
//{
//    cout << "TO BE DONE" << endl;
//}
//class C_OMP
//{
//    public:
//      int Na;  // number of atoms
//      int Npix; // number of pixels per atom
//      bool UseSparsityTarget;
//      int SparsityTarget; // target sparsity
//      double RecError;    // target reconstruction error
//
//      C_OMP() {RecError=0; SparsityTarget=DEF_DL_SPARSITY_NUMBER; Npix=0; Na; UseSparsityTarget=true;}
//      dblarray TabDico;
//      dblarray TabVect;
//      dblarray TabCoef;
//      void alloc();
//      void omp();
//     ~C_OMP() {} ;
//};
//
//cout << " Def Sparse Value = " << OMP.SparsityTarget << endl;
//if (SparsityTarget  != 0)    OMP.SparsityTarget  = SparsityTarget;
//cout << " Def Sparse Value = " << OMP.SparsityTarget << endl;
//if (UseSparsity  == false) OMP.UseSparsityTarget = false;
/*
dblarray Dico, Image_in, TabPatch, SigBlock;
	fits_read_fltarr (Name_Dico,Dico);
	int Nx = Dico.nx();
	int Ny = Dico.ny();
	int PatchSize = sqrt( double( Nx) );
	if (Verbose == True) cout << "Dico: Nx = " << Nx << " Ny = " << Ny << endl;

	fits_read_fltarr (Name_Image_in, Image_in);
	int dim = Image_in.naxis();
	Nx = Image_in.nx();
	Block1D Block;
	if (dim == 1)
	{
		if (Verbose == True) cout << "Data:  Nx = " << Nx << endl;
		Block._WeightFirst = False;
		// Block.BLOCKSIZE_Power_of_2 = False;
		Block._Verbose = True;
		Block._BlockOverlap = Ov;
		Block.alloc (Nx, PatchSize, False);
		if (Verbose == True) cout << "Block:  NbrBlock = " << Block.nbr_block() << ", SizeBlock = " <<  Block.block_size() << endl;
		if (Block._BlockOverlap  == True) cout << " OVERLAP " << endl;

		TabPatch.alloc(Block.block_size(), Block.nbr_block());
		SigBlock.alloc(Block.block_size());

		// Extract the patches and put them in TabPatch
		for (int b=0; b < Block.nbr_block(); b++)
		{
			Block.get_block_sig (b, Image_in, SigBlock);
			for (int i=0; i < Block.block_size(); i++) TabPatch(i,b) = SigBlock(i);
		}

		// For testing, reconstruct the signal from the patches
		fltarray DataRec;
		DataRec.alloc(Image_in.nx());
		for (int b=0; b < Block.nbr_block(); b++)
		{
			for (int i=0; i < Block.block_size(); i++) SigBlock(i) = TabPatch(i,b);
			Block.add_block_sig (b, DataRec, SigBlock);
		}
		fits_write_fltarr("xx_rec", DataRec);
		Image_in -= DataRec;
		Image_in.info("RESI");
		fits_write_fltarr("xx_resi", Image_in);
	}
	else
	{
		if (Verbose == True) cout << "Data:  Nx = " << Nx << ", Ny = " << Image_in.ny() << endl;
	}*/

/*	// Convert float array into double array
	initial_D.alloc(Dico.nx(),Dico.ny());
	training_set.alloc(Image_in.nx(),Image_in.ny());
	for (int i=0;i<Dico.nx();i++)
		for (int j=0;j<Dico.ny();j++)
			initial_D(i,j) = Dico(i,j);
	for (int i=0;i<Image_in.nx();i++)
		for (int j=0;j<Image_in.ny();j++)
			training_set(i,j) = Image_in(i,j);
			// Calling main function
	C_DL1D dl1d(training_set);
	learned_D = dl1d.dl1d(training_set, initial_D, IterationNumber,SparsityTarget,ErrorTarget);
 */
/*io_read_ima_float(Name_Image_in, Image_in);*/
