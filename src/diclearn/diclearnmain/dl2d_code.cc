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
#include "IM_Block2D.h"

bool UseSparsity=true;
int SparsityNumber=0;
/****************************************************************************/

#define DEF_DL_SPARSITY_NUMBER 5

class C_OMP
{
    public:
      int Na;  // number of atoms
      int Npix; // number of pixels per atom
      bool UseSparsityNumber;
      int SparsityNumber; // target sparsity
      double RecError;    // target reconstruction error
      
      C_OMP() {RecError=0; SparsityNumber=DEF_DL_SPARSITY_NUMBER; Npix=0; Na; UseSparsityNumber=true;}
      dblarray TabDico;
      dblarray TabVect;
      dblarray TabCoef;
      void alloc();
      void omp();
     ~C_OMP() {} ;
};


class C_Patch
{
public:
      C_Patch() {}
      ~C_Patch() {} ;
};

class C_Learning
{
public:
    C_Learning() {}
    ~C_Learning() {} ;
};

/****************************************************************************/

void C_OMP::alloc()
{
    Na = TabDico.ny();
    Npix = TabDico.nx();
    if (Npix != TabVect.nx())
    {
       cout << "Error: vector and dictonary dimensions do not match. " << endl;
       cout << "   Dico.nx = " << Npix << ", Vect.nx = " <<TabVect.nx() << endl;
    }
    TabCoef.alloc(Na);
}

/****************************************************************************/

void C_OMP::omp()
{
    cout << "TO BE DONE" << endl;
}

/****************************************************************************/

char Name_Dico_In[256];  // 2D array
char Name_Vect_In[256];  // 2D array
char Name_Coef[256];      // 

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose=False;

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options inDico inSignal_or_TabSignal output_coef \n\n", argv[0]);
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
    while ((c = GetOpt(argc,argv,"S:EvzZ:")) != -1) 
    {
	switch (c) 
        {
	  case 'S':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d", &SparsityNumber) != 1) 
                {
		    fprintf(OUTMAN, "bad value: %s\n", OptArg);
		    exit(-1);
		}
  		break;
     case 'E': UseSparsity = false; break;
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
	       if (UseSparsity  == True) cout << "Use Sparsity target " << endl;
           else cout << "Use Error target " << endl;
	       if (SparsityNumber  != 0) cout << " Sparsity target  =  " << SparsityNumber << endl;
    }
	
    if (Verbose == True) cout << "\n Reading the data"<< endl;

    C_OMP OMP;
    cout << " Def Sparse Value = " << OMP.SparsityNumber << endl;
    if (SparsityNumber  != 0)    OMP.SparsityNumber  = SparsityNumber;
    cout << " Def Sparse Value = " << OMP.SparsityNumber << endl;
    if (UseSparsity  == false) OMP.UseSparsityNumber = false;
    
    fltarray Dico, TabPatch;
    Ifloat Data,  SigBlock;
    fits_read_fltarr (Name_Dico_In,Dico);
    int Nx = Dico.nx();
    int Natoms = Dico.ny();
    int PatchSize = 64;  // (int) sqrt((double) Nx);
    if (Verbose == True) cout << "Dico: Nx = " << Nx << " Natoms = " << Natoms << ", PatchSize = " << PatchSize <<  endl;
    
    io_read_ima_float (Name_Vect_In, Data, &Header);
    Header.origin = Cmd;
    int dim = Data.naxis();
    int Nl = Data.nl();
    int Nc = Data.nc();
    
    Block2D Block;
    Block.Verbose = True;
    if (dim == 2)
    {
        if (Verbose == True) cout << "Data:  Nl = " << Nl << ", Nc = " << Nc << endl;
       Block.WeightFirst = False;
       Block.BlockOverlap = True;
       Block.alloc (Nl, Nc, PatchSize);
       if (Verbose == True) cout << "Block:  NbrBlock = " << Block.nbr_block() << ", SizeBlock = " <<  Block.block_size() << endl;

       int NpixPatch = Block.block_size() * Block.block_size();
       TabPatch.alloc(NpixPatch, Block.nbr_block());
       SigBlock.alloc(Block.block_size(),Block.block_size());
       
       // Extract the patches and put them in TabPatch 
       int b=0;
       for (int bl=0; bl < Block.nbr_block_nl(); bl++)
       for (int bc=0; bc < Block.nbr_block_nc(); bc++)
       {
          Block.get_block_ima (bl, bc, Data, SigBlock);
          int Ind=0;
          for (int i=0; i < Block.block_size(); i++) 
          for (int j=0; j < Block.block_size(); j++) TabPatch(Ind++, b) = SigBlock(i,j);
          b++;
       }
       
       // For testing, reconstruct the signal from the patches
       Ifloat DataRec;
       DataRec.alloc(Data.nl(), Data.nc());
       b=0;
       for (int bl=0; bl < Block.nbr_block_nl(); bl++)
        for (int bc=0; bc < Block.nbr_block_nc(); bc++)
        {
            int Ind=0;
            for (int i=0; i < Block.block_size(); i++) 
            for (int j=0; j < Block.block_size(); j++) SigBlock(i,j) = TabPatch(Ind++, b);
            Block.add_block_ima (bl, bc, DataRec, SigBlock);
            b++;
        }
        io_write_ima_float("xx_rec", DataRec, &Header);
        Data -= DataRec;
        Data.info("RESI");
        io_write_ima_float("xx_resi", Data, &Header);
        fits_write_fltarr("xx_p", TabPatch);
    }
    else 
    {
        if (Verbose == True) cout << "Data:  Nx = " << Nx << ", Ny = " << Data.ny() << endl;
    }

    exit(0);
}


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