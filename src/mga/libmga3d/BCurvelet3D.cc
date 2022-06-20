/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  2008
**    
**    File:  BCurvelet3D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Curvelet transform and reconstruction
**    -----------  
**                 
******************************************************************************/
#include "BCurvelet3D.h"
#include "IM_IO.h"
#include "PrimeNumber.h"
#include "IM_Prob.h"
#include "DefFunc.h"
#include "Atrou3D.h"

/****************************************************************************/

extern Bool Verbose;
static inline void PrintError( int status)
{
    // ***************************************************** 
    // * Print out cfitsio error messages and exit program * 
    // ***************************************************** 

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        // get the error status description 
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  // get first message; null if stack is empty 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  // get remaining messages 
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );
}

/****************************************************************************/

static inline void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/
static inline double prob_noise(float Val)
{
	return (abs(Val) < FLOAT_EPSILON) ? 1.:erfc (abs(Val));
}

/****************************************************************************/
// Noise lvl in bands
// see estim_norm_cur1g.pro for details on the coefficients
	static float TabBCur_Meyer_17[4][7] = 
		{{	3.78376, 4.09430, 4.10529, 3.17828, 3.52315, 3.50892, 2.72378 },
		{	1.64544, 2.58945, 2.60306, 2.95049, 3.32995, 3.32210, 2.46549 },
		{	1.05113, 1.73694, 1.77048, 2.43036, 0.00000, 0.00000, 0.00000 },
		{	0.842366, 1.61889, 1.63248, 2.38188, 0.00000, 0.00000, 0.00000 }};

	static float TabBCur_Meyer_17_3sig[4][7] = 
		{{	3.92788, 4.20274, 4.21170, 3.25643, 3.80357, 3.78887, 3.55010 },
		{	1.75088, 2.68958, 2.69785, 3.06782, 3.58265, 3.57191, 3.22947 },
		{	1.14747, 1.84979, 1.87310, 2.55453, 0.00000, 0.00000, 0.00000 },
		{	0.915724, 1.70667, 1.71718, 2.51870, 0.00000, 0.00000, 0.00000 }};

	static float TabBCur_Meyer_33[4][10] = // old
		{{	5.47129, 5.71601, 5.72477, 4.59588, 4.45688, 4.44795, 2.06938, 2.23711, 2.22759, 2.52496 },
		{	2.64297, 3.87423, 3.88361, 4.41954, 4.32986, 4.31881, 1.92685, 2.03634, 2.03009, 2.26586 },
		{	1.78590, 2.69517, 2.71140, 2.84519, 3.28725, 3.28759, 2.68110, 0.00000, 0.00000, 0.00000 },
		{	1.61798, 2.61111, 2.62907, 3.02098, 3.40020, 3.36733, 2.13804, 0.00000, 0.00000, 0.00000 }};

	static float TabBCur_Meyer_33_3sig[4][10] = // old
		{{	5.62827, 5.81000, 5.81700, 4.66501, 4.80762, 4.79771, 2.35064, 2.73436, 2.72076, 3.45612 },
		{	2.76736, 3.96562, 3.97162, 4.50437, 4.64329, 4.63142, 2.23716, 2.50872, 2.49831, 3.11943 },
		{	1.89978, 2.80349, 2.81548, 2.98647, 3.57948, 3.57460, 3.46197, 0.00000, 0.00000, 0.00000 },
		{	1.73262, 2.71471, 2.72721, 3.11356, 3.62386, 3.60052, 2.91144, 0.00000, 0.00000, 0.00000 }};

	static float TabBCur_Atrou_17[4][7] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_Atrou_17_3sig[4][7] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_Atrou_33[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_Atrou_33_3sig[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_PMeyer_17[4][7] = 
		{{	3.8940813, 4.1954791, 4.2039885, 3.1209968, 3.3486187, 3.3348118, 2.5527966 },
		{	3.0904676, 4.6892262, 4.7207080, 4.9903796, 5.4292305, 5.4044487, 3.8978228 },
		{	1.9801558, 3.1622190, 3.2439873, 4.0282394, 0.0000000, 0.0000000, 0.0000000 },
		{	2.1948695, 3.4274286, 3.5534175, 3.9948936, 0.0000000, 0.0000000, 0.0000000 }};

	static float TabBCur_PMeyer_17_3sig[4][7] = // copy of TabBCur_PMeyer_17
		{{	3.8940813, 4.1954791, 4.2039885, 3.1209968, 3.3486187, 3.3348118, 2.5527966 },
		{	3.0904676, 4.6892262, 4.7207080, 4.9903796, 5.4292305, 5.4044487, 3.8978228 },
		{	1.9801558, 3.1622190, 3.2439873, 4.0282394, 0.0000000, 0.0000000, 0.0000000 },
		{	2.1948695, 3.4274286, 3.5534175, 3.9948936, 0.0000000, 0.0000000, 0.0000000 }};

	static float TabBCur_PMeyer_33[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_PMeyer_33_3sig[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_PMeyerC_17[4][7] = 
		{{	3.8940813, 4.1954791, 4.2039885, 3.1209968, 3.3486187, 3.3348118, 2.5527966 },
		{	3.0904676, 4.6892262, 4.7207080, 4.9903796, 5.4292305, 5.4044487, 3.8978228 },
		{	1.9801558, 3.1622190, 3.2439873, 4.0282394, 0.0000000, 0.0000000, 0.0000000 },
		{	2.1948695, 3.4274286, 3.5534175, 3.9948936, 0.0000000, 0.0000000, 0.0000000 }};

	static float TabBCur_PMeyerC_17_3sig[4][7] = // copy of TabBCur_PMeyer_17
		{{	3.8940813, 4.1954791, 4.2039885, 3.1209968, 3.3486187, 3.3348118, 2.5527966 },
		{	3.0904676, 4.6892262, 4.7207080, 4.9903796, 5.4292305, 5.4044487, 3.8978228 },
		{	1.9801558, 3.1622190, 3.2439873, 4.0282394, 0.0000000, 0.0000000, 0.0000000 },
		{	2.1948695, 3.4274286, 3.5534175, 3.9948936, 0.0000000, 0.0000000, 0.0000000 }};

	static float TabBCur_PMeyerC_33[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

	static float TabBCur_PMeyerC_33_3sig[4][10] = 
		{{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

/****************************************************************************/

BCurvelet3D::BCurvelet3D()
{
	Ptr_SB1D=NULL;
	SB1D=NULL;
	TabBlockSize=NULL;
	TabBeamlet=NULL;
	CubeSigma=NULL;
	InitClass=false;
	GetAutoNbScale=False;
	eval_stat=false;
	no_recons=false;
	TabStat=NULL;
	Use_3sigma=false;
	no_coarse=false;
	keep_energy=false;
	BlockOverlap=0;
	CurTrans=(type_linelet3d_WTtrans)0;
	Type_PMW3D=(type_pmeyer3D)DEF_TYPE_PMW3D;
	
	AllocClass=False;
}

BCurvelet3D::~BCurvelet3D()
{
	if(Verbose) cerr<<"BCurvelet3D::~BCurvelet3D()"<<endl;
	dealloc();
	deinit();
}

void BCurvelet3D::dealloc()
{
	if(TabBeamlet!=NULL) delete [] TabBeamlet;
	AllocClass=False;
}

void BCurvelet3D::deinit()
{
	if(CubeSigma!=NULL) delete [] CubeSigma;
	if(TabStat!=NULL) delete [] TabStat;
	if(TabBlockSize!=NULL) delete [] TabBlockSize;
	InitClass=false;
}

/****************************************************************************/

int BCurvelet3D::s2d_y0(int s3, int s)
{
	int NbPoint=TabBlockSize[s3];
	int ScaleNbPt,BegInd;
	
	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*get_nbscale2d(s3)-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return BegInd;
		else if(s%3==1) return BegInd;
		else return 0;//(s%3==2)
	}
	else // last scale
		return 0;
}

int BCurvelet3D::s2d_yn(int s3, int s)
{
	int NbPoint=TabBlockSize[s3];
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*get_nbscale2d(s3)-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return ScaleNbPt;
		else if(s%3==1) return ScaleNbPt;
		else return BegInd;//(s%3==2)
	}
	else // last scale
		return BegInd;
}

int BCurvelet3D::s2d_x0(int s3, int s)
{
	int NbPoint=TabBlockSize[s3];
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*get_nbscale2d(s3)-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return BegInd;
		else if(s%3==1) return 0;
		else return BegInd;//(s%3==2)
	}
	else // last scale
		return 0;
}

int BCurvelet3D::s2d_xn(int s3, int s)
{
	int NbPoint=TabBlockSize[s3];
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*get_nbscale2d(s3)-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return ScaleNbPt;
		else if(s%3==1) return BegInd;
		else return ScaleNbPt;//(s%3==2)
	}
	else // last scale
		return BegInd;
}

/****************************************************************************/

void BCurvelet3D::calib_noise_nsig(fltarray * TabBand, float N_Sigma, char* Name_Imag_Out)
{
	if(Verbose) cerr << "BCurvelet3D::calib_noise_nsig..." <<  endl;
//	bool LocVerbose = false & Verbose;
	double LMax;
	fltarray Band;
	
	TabSigma.init(0); // useless when there is no bug
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int b=0;b<get_nbband(s3);b++)
		{
			//float Nsig = (s3 == 0) ? N_Sigma+1: N_Sigma;
			float Nsig = N_Sigma;
			get_band(s3,b,TabBand, Band);
			
			CArrProb *TabCP;
			TabCP = new CArrProb;
			
			TabCP->set(Band,true);
			TabCP->find_gthreshold(Nsig, LMax);
			TabSigma(s3,b) = LMax / N_Sigma;
			delete TabCP;
			if(Verbose) cerr << "Band (" << s3<<","<<b << ") Lmax = " << TabSigma(s3,b)<< endl;
		}
		
	char filename[64];
	sprintf(filename,"%s_nsig.fits",Name_Imag_Out);
	writefltarr(filename, TabSigma);
	
	if(Verbose) cerr << "End BCurvelet3D::calib_noise_nsig" <<  endl;
}

/****************************************************************************/
 
void BCurvelet3D::tab_block_size_init()
{
	int s, B;
	
	TabBlockSize = new int[NbrScale3D-1];

// Block size must be a power of 2
	B = (int) next_power_of_2((int) BlockSize-1);

	switch (TypeBlock)
	{
		case BL3D_CONST: 
			for (s=0; s < NbrScale3D-1; s++) TabBlockSize[s] = B+1;
			break;
		case BL3D_UP:
			for (s=0; s < NbrScale3D-1; s++)
			{
				TabBlockSize[s] = B+1;
				B *= 2;
			}
			break;
		case BL3D_UP2: // 1 2 2 4 4 8 8
			for (s=0; s < NbrScale3D-1; s++)
			{
				TabBlockSize[s] = B+1;
				if (s % 2 == 0) B *= 2;
			}
			break;
		case BL3D_DOWN:
			for (s=0; s < NbrScale3D-1; s++)
			{
				TabBlockSize[s] = B+1;
				B /= 2;
				if (B < 8) B = 8;
			}
			break;
		case BL3D_DOWN2: // 8 4 4 2 2 1 1
			for (s=0; s < NbrScale3D-1; s++)
			{
				TabBlockSize[s] = B+1;
				if (s % 2 != 0) 
				{
					B /= 2;
					if (B < 8) B = 8;
				}
			}
			break;
	}
}

/****************************************************************************/

void BCurvelet3D::wavelet_transform(fltarray & Cube, fltarray * & TabWaveBand, bool alloc)
{
	if(TypeW3D==W3D_ATROU)
	{
		if(Verbose) cerr << " Atrou WT"<<endl;
		ATROUS_3D_WT atrou;
		atrou.alloc(TabWaveBand, NxCube, NyCube, NzCube, NbrScale3D);
		atrou.transform(Cube, TabWaveBand, NbrScale3D);
	}
	else if(TypeW3D==W3D_MEYER)
	{
		if(Verbose) cerr << " FFT based Meyer WT, size = "<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		MEYER_WT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.transform(Cube, TabWaveBand, (Bool)alloc);
		DataW.normalize_self(TabWaveBand);
	}
	else if(TypeW3D==W3D_PMEYER)
	{
		if(Verbose) cerr << " FFT based Poisson Meyer WT, size = "<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		POISSON_MWT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.Type_PMW3D=Type_PMW3D;
		DataW.transform(Cube, TabWaveBand, (Bool)alloc);
		DataW.normalize_self(TabWaveBand);
	}
	else
	{
		cerr<< "Undefined wavelet3D transform"<<endl;
		exit(-1);
	}
}

/****************************************************************************/

void BCurvelet3D::reset()
{
	NxCube=0;
	NyCube=0;
	NzCube=0;
	BlockSize=0;
	NbrScale3D=0;
	TypeW3D=(type_wavelet3D)0;
	SB1D=NULL;
	Border=I_MIRROR;
	TypeBlock=(type_3dcurvelet_block)0;
 //pas en entrée de init
	//BlockOverlap=0;
	//CurTrans=(type_linelet3d_WTtrans)0;
	//Type_PMW3D=(type_pmeyer3D)DEF_TYPE_PMW3D;
}

void BCurvelet3D::init(int Nx, int Ny, int Nz, int _NbrScale3D, int BS, type_wavelet3D _TypeW3D, SubBandFilter* _SB1D)
{
	if (Verbose) cout << "BCurvelet3D::init..." << endl;
	
	if(InitClass) deinit();
	InitClass=true;
	reset();
	
// Set parameters / Default parameters
	NxCube		= Nx;
	NyCube		= Ny;
	NzCube		= Nz;
	BlockSize	= BS > 0 ? BS : DEF_3DCUR_BLOCK_SIZE;
	NbrScale3D	= _NbrScale3D > 0 ? _NbrScale3D : DEF_3DCUR_NBR_SCALE;
	TypeW3D		= _TypeW3D;
	SB1D		= _SB1D;
	
	if(AllocClass) dealloc();

// Constrain the BlockSize to two values
	if(BlockSize>33) BlockSize = 33;
	if(BlockSize<33) BlockSize = 17;
	
	switch(TypeW3D)
	{
		case W3D_MEYER:
		case W3D_PMEYER:
//		case W3D_PMEYERBIS:
			TypeBlock = BL3D_DOWN2;
			break;
			
		case W3D_ATROU:
			TypeBlock = BL3D_UP2;
			break;
			
		default:
			TypeBlock = DEF_3DCUR_TYPE_BLOCK;
			break;
	}
	if(CurTrans<0 || CurTrans > NBR_LIN3D_TRANS) CurTrans = DEF_LIN3D_TRANS;
	if(TypeW3D<0 || TypeW3D > NBR_TYPE_W3D) TypeW3D = DEF_TYPE_W3D;

// Block Sizes calculus
	tab_block_size_init();
	
// TabSigma initialisation
	int NbrBand2D = 3*(iround((float)log((float) (TabBlockSize[0] / 4. * 3.) / log(2.)))-1)+1; // 3*NbScale2D+1
	TabSigma.resize(NbrScale3D-1,NbrBand2D);

	if(!Use_3sigma)
	{
		if (TypeW3D==W3D_ATROU)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Atrou_33[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Atrou_17[i][j];
		}
		else if(TypeW3D==W3D_PMEYER) 
		{
			if(Type_PMW3D==PMW3D_SUM)
			{
				if (BlockSize==33)
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyer_33[i][j];
				else //(BlockSize==17) and default
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyer_17[i][j];
			}
			else if(Type_PMW3D==PMW3D_COHERENT)
			{
				if (BlockSize==33)
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyerC_33[i][j];
				else //(BlockSize==17) and default
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyerC_17[i][j];
			}
		}
		else // W3D_MEYER and default
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Meyer_33[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Meyer_17[i][j];
		}
	}
	else // use 3sigma cumulative estimation
	{
		if (TypeW3D==W3D_ATROU)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Atrou_33_3sig[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Atrou_17_3sig[i][j];
		}
		else if(TypeW3D==W3D_PMEYER)
		{
			if(Type_PMW3D==PMW3D_SUM)
			{
				if (BlockSize==33)
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyer_33_3sig[i][j];
				else //(BlockSize==17) and default
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyer_17_3sig[i][j];
			}
			else if(Type_PMW3D==PMW3D_COHERENT)
			{
				if (BlockSize==33)
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyerC_33_3sig[i][j];
				else //(BlockSize==17) and default
					for(int i=0;i<NbrScale3D-1;i++)
						for(int j=0;j<NbrBand2D;j++)
							TabSigma(i,j)=TabBCur_PMeyerC_17_3sig[i][j];
			}
		}
		else // W3D_MEYER and default
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Meyer_33_3sig[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrBand2D;j++)
						TabSigma(i,j)=TabBCur_Meyer_17_3sig[i][j];
		}
	}

// CubeSigma initialisation
	char filename[128];
	CubeSigma = new fltarray[4];
	for(int s3=0;s3<4;s3++)
	{
		sprintf(filename,"%s/common/BC_cubesigma/BS%d/b%d_mean.fits",getenv("CUR"),BlockSize,s3);
			if(access(filename, F_OK)==0)
				fits_read_fltarr(filename, CubeSigma[s3]);
			else 
			{
				if(Verbose) cerr<<"Sigma map "<<filename<<" not found."<<endl;
				// allocate cubesigma with the correct size, set to 1
			}
	}
	
// TabStat initialisation
	if(eval_stat) TabStat = new dblarray [NbrScale3D+1];
	else TabStat=NULL;
	
	if (Verbose) cerr << "...End BCurvelet3D::init" << endl;
}

/****************************************************************************/

// Initialize the beamlet transform
void BCurvelet3D::beamlet_init(Linelet3D &Beam)
{
	Beam.set_BlockOverlap((BlockOverlap>0)?True:False);
	Beam.set_OverlapFactor(BlockOverlap);
	Beam.set_Verbose(Verbose);
	Beam.set_LinTransf(CurTrans, SB1D);
	Beam.set_keep_energy(keep_energy);
}

/****************************************************************************/

void BCurvelet3D::alloc()
{
	bool LocVerbose = false & Verbose;
	if (Verbose) cout << "BCurvelet3D::alloc..." << endl;
	
	AllocClass = True;

// Beamlet initialisation
	if (LocVerbose) cerr << " Beamlet init" << endl;
	TabBeamlet = new Linelet3D[NbrScale3D-1];
	for (int s = 0; s <  NbrScale3D-1; s++) 
		beamlet_init(TabBeamlet[s]);

	if (Verbose) cerr << "...End BCurvelet3D::alloc" << endl;
}


/****************************************************************************/

void BCurvelet3D::transform(fltarray &Cube, fltarray* & TabBand, bool TBalloc)
{
	if(Verbose) cerr << "BCurvelet3D::transform..." <<  endl;
	
	fltarray* TabWaveBand;
	
// Wavelet transform of the Cube
	wavelet_transform(Cube, TabWaveBand, true);
	
// Multiscale ridgelet transform
	transform(TabWaveBand, TabBand, TBalloc);
	
	delete[] TabWaveBand;
	
	if(Verbose) cerr << "...End BCurvelet3D::transform" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::transform(fltarray * & TabWaveBand, fltarray* & TabBand, bool TBalloc)
{
	bool LocVerbose=false & Verbose;
	if(Verbose) cerr << "BCurvelet3D::transform_band..." <<  endl;
	
// Dimentions of the coarsest scale
	TabSizeW.alloc(NbrScale3D,3,1);
	for(int i=0;i<NbrScale3D;i++)
	{
		TabSizeW(i,0,0)=TabWaveBand[i].nx();
		TabSizeW(i,1,0)=TabWaveBand[i].ny();
		TabSizeW(i,2,0)=TabWaveBand[i].nz();
	}
	int Nx = TabWaveBand[NbrScale3D-1].nx();
	int Ny = TabWaveBand[NbrScale3D-1].ny();
	int Nz = TabWaveBand[NbrScale3D-1].nz();
	
	alloc();

// Allocation TabBand
	if(TBalloc)
	{
		TabBand = new fltarray [NbrScale3D];
		TabBand[NbrScale3D-1].alloc(Nx,Ny,Nz);
	}
	
// Beamlet transformation on each scale except the coarsest	
	#if USE_OMP_BC
	omp_set_nested(1);
	#pragma omp parallel for 
	#endif
	for (int s = 0; s <  NbrScale3D-1; s++)
	{
		if (LocVerbose) cout << " Scale " << s << " Block size calcul = " << TabBlockSize[s] <<", 1stW_coef="<<TabWaveBand[s](0,0,0)<< endl;
		
		TabBeamlet[s].transform(TabWaveBand[s], TabBand[s], TabBlockSize[s]);
		
		if(LocVerbose)
		{
			char filename[64];
			sprintf(filename,"TWB_%d.fits",s);
			writefltarr(filename, TabWaveBand[s]);
			sprintf(filename,"TB_%d.fits",s);
			writefltarr(filename, TabBand[s]);
		}
	}
	
// keep the coarsest scale unchanged
	int s = NbrScale3D-1;
	TabBand[s] = TabWaveBand[s];
	if(LocVerbose)
	{
		char filename[64];
		sprintf(filename,"TWB_%d.fits",s);
		writefltarr(filename, TabWaveBand[s]);
	}
    
	if(Verbose) cerr << endl << "...End BCurvelet3D::transform_band" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::filter(fltarray & Cube, fltarray & Recons, float SigmaNoise, float NSigma, filter_type FilterType, 
		float Alpha, // FDR
		int LocalBS, // WIENER
		bool force4, 
		bool UseCubeSigma, 
		char* Outname)
{
	bool LocVerbose=false & Verbose;
	if(Verbose) cerr << "BCurvelet3D::filter..." <<  endl;
	
	fltarray* TabWaveBand;
	fltarray* TabWaveBandOut = new fltarray[NbrScale3D];
	
// Wavelet transform of the Cube
	wavelet_transform(Cube, TabWaveBand, true);
	
// Dimentions of the coarsest scale
	TabSizeW.alloc(NbrScale3D,3,1);
	for(int i=0;i<NbrScale3D;i++)
	{
		TabSizeW(i,0,0)=TabWaveBand[i].nx();
		TabSizeW(i,1,0)=TabWaveBand[i].ny();
		TabSizeW(i,2,0)=TabWaveBand[i].nz();
	}
	alloc();
	
// Beamlet transformation on each scale except the coarsest	
	for (int s3 = 0; s3 <  NbrScale3D-1; s3++)
	{
		if (LocVerbose) cout << " Scale " << s3 << " Block size calcul = " << TabBlockSize[s3] << endl;
		
		TabBeamlet[s3].set_no_recons(no_recons);
		TabBeamlet[s3].filter(TabWaveBand[s3], TabWaveBandOut[s3], TabBlockSize[s3], 
						s3, SigmaNoise, NSigma, FilterType, Alpha, LocalBS, force4, UseCubeSigma, 
						&TabSigma, CubeSigma, TabStat, NbrScale3D);
		
		if(LocVerbose)
		{
			char filename[64];
			sprintf(filename,"TWB_%d.fits",s3);
			writefltarr(filename, TabWaveBand[s3]);
			sprintf(filename,"TWBO_%d.fits",s3);
			writefltarr(filename, TabWaveBandOut[s3]);
		}
	}

// keep the coarsest scale unchanged
	int s = NbrScale3D-1;
	TabWaveBandOut[s] = TabWaveBand[s];
	if(no_coarse)
		TabWaveBandOut[s].init(0.0);
	
	if(LocVerbose)
	{
		char filename[64];
		sprintf(filename,"TWB_%d.fits",s);
		writefltarr(filename, TabWaveBand[s]);
	}
    
	if(eval_stat)
	{
	// Output stat file
		char Statname[250];

		strcpy(Statname, Outname);
		strcat(Statname, "_stat.dat");
		fstream cstat;
		cstat.open (Statname, fstream::out);
	
	// Coarse scale stats
		int s3=NbrScale3D-1;
		float Mini,Maxi;
		TabStat[s3].alloc(1,1,6);
		moment4(TabWaveBandOut[s3].buffer(),TabWaveBandOut[s3].n_elem(),TabStat[s3](0,0,0), TabStat[s3](0,0,1), TabStat[s3](0,0,2),
				TabStat[s3](0,0,3), Mini, Maxi);
		TabStat[s3](0,0,4)=Mini;
		TabStat[s3](0,0,5)=Maxi;
		if(LocVerbose) if(Outname!=NULL) cstat << s3 << TabStat[s3](0,0,0) <<"\t"<<  TabStat[s3](0,0,1) <<"\t"<< TabStat[s3](0,0,2)
				<<"\t"<< TabStat[s3](0,0,3) <<"\t"<< TabStat[s3](0,0,4) <<"\t"<< TabStat[s3](0,0,5) <<endl;
		
		// Global stats
		TabStat[NbrScale3D](s3,0,0)=TabStat[s3](0,0,0);
		TabStat[NbrScale3D](s3,0,1)=TabStat[s3](0,0,1);
		TabStat[NbrScale3D](s3,0,2)=TabStat[s3](0,0,2);
		TabStat[NbrScale3D](s3,0,3)=TabStat[s3](0,0,3);
		TabStat[NbrScale3D](s3,0,4)=Mini;
		TabStat[NbrScale3D](s3,0,5)=Maxi;
	
	// Centering the global stats
		for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s=0;s<get_nbband(s3);s++)
		{
			int N = nbr_block(s3);
			int Ntot = get_nbplane(s3) * s2d_xn(s3,s)*s2d_yn(s3,s) ;
			moment4_center(Ntot, TabStat[NbrScale3D](s3,s,0)/N, TabStat[NbrScale3D](s3,s,1)/N, TabStat[NbrScale3D](s3,s,2)/N, TabStat[NbrScale3D](s3,s,3)/N,
					TabStat[NbrScale3D](s3,s,1), TabStat[NbrScale3D](s3,s,2), TabStat[NbrScale3D](s3,s,3));
			TabStat[NbrScale3D](s3,s,0) /= N;
			if(LocVerbose)
				cerr << s3 <<"\t"<< s <<"\t"<< TabStat[NbrScale3D](s3,s,0) <<"\t"<< TabStat[NbrScale3D](s3,s,1) <<"\t"<< TabStat[NbrScale3D](s3,s,2) <<"\t"<<
						 TabStat[NbrScale3D](s3,s,3) <<"\t"<< TabStat[NbrScale3D](s3,s,4) << "\t" << TabStat[NbrScale3D](s3,s,5) << endl;
			if(Outname!=NULL)
				cstat << s3 <<"\t"<< s <<"\t"<< TabStat[NbrScale3D](s3,s,0) <<"\t"<< TabStat[NbrScale3D](s3,s,1) <<"\t"<< TabStat[NbrScale3D](s3,s,2) <<"\t"<<
						 TabStat[NbrScale3D](s3,s,3) <<"\t"<< TabStat[NbrScale3D](s3,s,4) << "\t" << TabStat[NbrScale3D](s3,s,5) << endl;
		}
		
	// Write the complete stat file
		if(Outname!=NULL) 
		{
			char filename[64];
			for(int s3=0;s3<NbrScale3D;s3++)
			{
				sprintf(filename,"%s_stat_full_s%d.fits",Outname,s3);
				writefltarr(filename, TabStat[s3]);
			}
			sprintf(filename,"%s_stat.fits",Outname);
			writefltarr(filename, TabStat[NbrScale3D]);
			cstat.close();
		}
	}

	if(!no_recons)
	{
		Recons.alloc(Cube.nx(),Cube.ny(),Cube.nz());
		// Wavelet inverse transform of the Cube
		wavelet_recons(TabWaveBandOut, Recons);
	}
	delete [] TabWaveBand;
	delete [] TabWaveBandOut;

	if(Verbose) cerr << endl << "...End BCurvelet3D::filter" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::wavelet_recons(fltarray * & TabWaveBand, fltarray & Cube)
{
	if(TypeW3D==W3D_ATROU)
	{
		if(Verbose) cerr << " Atrou WT"<<endl;
		ATROUS_3D_WT atrou;
		atrou.recons(TabWaveBand, Cube, NbrScale3D);
	}
	else if(TypeW3D==W3D_MEYER)
	{
		if(Verbose) cerr << " FFT based Meyer WT, size ="<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		MEYER_WT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.normalize_self(TabWaveBand, true);
		DataW.recons(TabWaveBand,Cube);
	}
	else if(TypeW3D==W3D_PMEYER)
	{
		if(Verbose) cerr << " FFT based Meyer WT, size ="<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		POISSON_MWT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.Type_PMW3D=Type_PMW3D;
		DataW.normalize_self(TabWaveBand, true);
		DataW.recons(TabWaveBand,Cube);
	}
	else
	{
		cerr<< "Undefined wavelet3D transform"<<endl;
		exit(-1);
	}
}

/****************************************************************************/

void BCurvelet3D::init_recons(fltarray * & TabWaveBand, fltarray* &TabBand)
{
	if(Verbose) cerr << "BCurvelet3D::init_recons..." <<  endl;

// Allocation TabWaveBand
//	TabWaveBand = new fltarray [NbrScale3D];
	for (int s = 0; s <  NbrScale3D-1; s++)
	{
//cerr<<"init tbeam size "<<TabBeamlet[s].get_SizeX()<<" "<<TabBeamlet[s].get_SizeY()<<" "<<TabBeamlet[s].get_SizeZ()<<endl;
		TabWaveBand[s].alloc(TabBeamlet[s].get_SizeX(),TabBeamlet[s].get_SizeY(),TabBeamlet[s].get_SizeZ());
//cerr<<" init twb size TWB("<<s<<") "<<TabWaveBand[s].nx()<<" "<<TabWaveBand[s].ny()<<" "<<TabWaveBand[s].nz()<<endl;
	}
	TabWaveBand[NbrScale3D-1].alloc(TabBand[NbrScale3D-1].nx(),TabBand[NbrScale3D-1].ny(),TabBand[NbrScale3D-1].nz());
	if(Verbose) cerr << "End BCurvelet3D::init_recons" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::recons(fltarray* &TabBand, fltarray* &TabWaveBand, bool TWBalloc)
{
	bool LocVerbose=false && Verbose;
	if(Verbose) cerr << "BCurvelet3D::recons_band..." <<  endl;
	
	
	if(TWBalloc)
	{
		TabWaveBand = new fltarray [NbrScale3D];
		for (int s = 0; s <  NbrScale3D-1; s++)
		{
	//cerr<<"init tbeam size "<<TabBeamlet[s].get_SizeX()<<" "<<TabBeamlet[s].get_SizeY()<<" "<<TabBeamlet[s].get_SizeZ()<<endl;
			TabWaveBand[s].alloc(TabBeamlet[s].get_SizeX(),TabBeamlet[s].get_SizeY(),TabBeamlet[s].get_SizeZ());
	//cerr<<" init twb size TWB("<<s<<") "<<TabWaveBand[s].nx()<<" "<<TabWaveBand[s].ny()<<" "<<TabWaveBand[s].nz()<<endl;
		}
		TabWaveBand[NbrScale3D-1].alloc(TabBand[NbrScale3D-1].nx(),TabBand[NbrScale3D-1].ny(),TabBand[NbrScale3D-1].nz());
//		init_recons(TabWaveBand, TabBand);
//		AllocRecon=True;
	}
	
// Beamlet inverse transformation on each scale except the coarsest	
#if USE_OMP_BC
		#pragma omp parallel for 
#endif
	for (int s = 0; s <  NbrScale3D-1; s++)
	{
		if (LocVerbose) cout << " Scale " << s << " Block size calcul = " << TabBlockSize[s] << endl;

		TabBeamlet[s].recons(TabBand[s], TabWaveBand[s]);
		if(LocVerbose)
		{
			char filename[64];
			sprintf(filename,"TWB_recons_%d.fits",s);
			writefltarr(filename, TabWaveBand[s]);
			sprintf(filename,"TB_read_%d.fits",s);
			writefltarr(filename, TabBand[s]);
		}
	}
	
// keep the coarsest scale unchanged
	int s = NbrScale3D-1;
	TabWaveBand[s] = TabBand[s];
	if(LocVerbose)
	{
		char filename[64];
		sprintf(filename,"TWB_recons_%d.fits",s);
		writefltarr(filename, TabWaveBand[s]);
	}
	
	if(Verbose) cerr << endl << "...End BCurvelet3D::recons_band" <<  endl;	
}

/****************************************************************************/

void BCurvelet3D::recons(fltarray* &TabBand, fltarray &Cube, bool Cubealloc)
{
	if(Verbose) cerr << "BCurvelet3D::recons..." <<  endl;
	
	fltarray* TabWaveBand;
	
	if(Cubealloc) Cube.alloc(NxCube,NyCube,NzCube);
	
// Multiscale beamlet inverse transform
	recons(TabBand, TabWaveBand);
	
// Wavelet inverse transform of the Cube
	wavelet_recons(TabWaveBand, Cube);
	
	delete [] TabWaveBand;
	
	if(Verbose) cerr << "...End BCurvelet3D::recons" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::get_band(int s3, int s, fltarray *TabBand, fltarray &band)
{
//	if(Verbose) cerr << "BCurvelet3D::get_band..." <<  endl;

	if(s3==NbrScale3D-1)
		band=TabBand[s3];
	else
	{
		int sx=s2d_x0(s3,s);
		int sy=s2d_y0(s3,s);
		int nx=s2d_xn(s3,s);
		int ny=s2d_yn(s3,s);
		
		band.resize(nx,ny,get_nbplane(s3));
		
		for(int x = sx; x < sx+nx; x++)
			for(int y = sy; y < sy+ny; y++)
				for(int z = 0; z < get_nbplane(s3); z++)
					band(x-sx,y-sy,z) = (TabBand[s3])(x,y,z);
	}
//	if(Verbose) cerr << "End BCurvelet3D::get_band" <<  endl;
}

/****************************************************************************/

void BCurvelet3D::noise_calibration(fltarray *TabBand, char* Outname)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"Noise_calibration..."<<endl;
		
	double ** TabStat = new double*[NbrScale3D-1];
	
// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_noise_calib.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	
	
	// fine scales : 0..N-2
	// Parallélisable : échelles indépendantes
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		if(LocVerbose) cerr<<"scale s3="<<s3<<" nb_scale="<<get_nbscale2d(s3)<<" nb_band="<<get_nbband(s3)<<endl;
		// Allocation
		TabStat[s3] = new double[get_nbband(s3)];// 3 dir each fine scale + 1 coarse scale
		int NbPlane=get_nbplane(s3);
		
		for(int s=0;s<get_nbband(s3);s++)
		{
			//int s=3*s2;
			double mean=0;	// mean
			double var=0;	// variance

			// mean
			for(int i = s2d_x0(s3,s); i < s2d_x0(s3,s)+s2d_xn(s3,s); i++)
				for(int j = s2d_y0(s3,s); j < s2d_y0(s3,s)+s2d_yn(s3,s); j++)
					for(int k = 0; k < NbPlane; k++)
					{
						mean += (TabBand[s3])(i,j,k);
						var += pow( (TabBand[s3])(i,j,k) , 2 );
					}
			mean /= (s2d_xn(s3,s)*s2d_yn(s3,s)*NbPlane) ;
			var  /= (s2d_xn(s3,s)*s2d_yn(s3,s)*NbPlane) ;

			// Standard deviation
			TabStat[s3][s] = sqrt(var);

			if(LocVerbose) cerr<<" ("<<s2d_x0(s3,s)<<":"<<s2d_x0(s3,s)+s2d_xn(s3,s)-1<<")x("<<s2d_y0(s3,s)<<":"<<s2d_y0(s3,s)+s2d_yn(s3,s)-1<<")";
			if(Verbose)	cerr<< "\tStat echelle ("<<s3<<","<<s<<") : (mean,std) = ("<<mean<<", \t"<<TabStat[s3][s]<<")"<<endl;
			cstat <<mean<<"\t"<<TabStat[s3][s]<<endl;
		}
	}
	
	// coarse scale : N-1
	// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats
	
	cstat.close();
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		delete [] TabStat[s3];
	delete [] TabStat;
	
	if(Verbose) cerr<<"...End Noise_calibration"<<endl;
}

/****************************************************************************/

void BCurvelet3D::extract_stat(fltarray *TabBand, char* Outname, bool UseCubeSigma, bool normalize)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"Extract_stat..."<<endl;
	float val;
	
// Output stat file
	char Statname[250];
	fstream cstat;
	if(Outname!=NULL)
	{
		strcpy(Statname, Outname);
		strcat(Statname, "_stat.dat");
		cstat.open (Statname, fstream::out);
	}
	
// Global Stats
	TabStat[NbrScale3D].alloc(NbrScale3D,get_nbband(0),6);

// fine scales : 0..N-2
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		if(LocVerbose) cerr<<"scale s3="<<s3<<" nb_scale="<<get_nbscale2d(s3)<<" nb_band="<<get_nbband(s3)<<endl;
		TabStat[s3].alloc(nbr_block(s3),get_nbband(s3),6);
		
		int depth = 3*TabBlockSize[s3]*TabBlockSize[s3];
		
		// note: parallélisable
		for(int n=0;n<nbr_block(s3);n++)
		for(int s=0;s<get_nbband(s3);s++)
		{
			fltarray temp(s2d_xn(s3,s),s2d_yn(s3,s),depth);
			float norm = 1;
			// Normalization and recopy
			for(int i = 0; i < s2d_xn(s3,s); i++)
				for(int j = 0; j < s2d_yn(s3,s); j++)
					for(int k = 0; k < depth; k++)
					{
						if(normalize)
						{
							if(UseCubeSigma) norm = (CubeSigma[s3])(i+s2d_x0(s3,s),j+s2d_y0(s3,s),k);
							else norm = TabSigma(s3,s);
						}
						temp(i,j,k)=(TabBand[s3])(i+s2d_x0(s3,s),j+s2d_y0(s3,s),k+n*depth)/norm;
					}
			double Mean, Sig, Skew, Curt;
			float Mini, Maxi;
			moment4(temp.buffer(), temp.n_elem(), Mean, Sig, Skew, Curt, Mini, Maxi, true); // not centered
			TabStat[s3](n,s,0)=Mean;
			TabStat[s3](n,s,4)=Mini;
			TabStat[s3](n,s,5)=Maxi;
			int Ntot = depth * s2d_xn(s3,s)*s2d_yn(s3,s) ;
			moment4_center(Ntot, Mean, Sig, Skew, Curt, TabStat[s3](n,s,1), TabStat[s3](n,s,2), TabStat[s3](n,s,3));
			if(LocVerbose)
			{
				cerr<<" ("<<s2d_x0(s3,s)<<":"<<s2d_x0(s3,s)+s2d_xn(s3,s)-1<<")x("<<s2d_y0(s3,s)<<":"<<s2d_y0(s3,s)+s2d_yn(s3,s)-1<<")";
				if(Verbose) cerr<< "\tStat echelle ("<<s3<<","<<s<<") : (mean,std) = ("<<Mean<<", \t"<<Sig<<")"<<endl;
				if(Verbose) cerr << " MaxCoef Value = "<<Mini<<","<<Maxi << endl;
				if(Outname!=NULL)
				cstat << s3 << "\t" << n << "\t" << s <<"\t"<< 
						TabStat[s3](n,s,0) <<"\t"<< TabStat[s3](n,s,1) <<"\t"<< TabStat[s3](n,s,2)<<"\t"<<TabStat[s3](n,s,3)<<"\t"<<
						TabStat[s3](n,s,4) <<"\t"<< TabStat[s3](n,s,5) <<"\t"<< endl;
			}
			
			// Global stats
			TabStat[NbrScale3D](s3,s,0)+=Mean;
			TabStat[NbrScale3D](s3,s,1)+=Sig;
			TabStat[NbrScale3D](s3,s,2)+=Skew;
			TabStat[NbrScale3D](s3,s,3)+=Curt;
			if(n==0)
			{
				TabStat[NbrScale3D](s3,s,4) = Mini;
				TabStat[NbrScale3D](s3,s,5) = Maxi;
			}
			else
			{
				TabStat[NbrScale3D](s3,s,4) = TabStat[NbrScale3D](s3,s,4) < Mini ? TabStat[NbrScale3D](s3,s,4) : Mini;
				TabStat[NbrScale3D](s3,s,5) = TabStat[NbrScale3D](s3,s,5) > Maxi ? TabStat[NbrScale3D](s3,s,5) : Maxi;
			}
		}
	}
// Coarse scale
	{
		int s3=NbrScale3D-1;
		float Mini,Maxi;
		TabStat[s3].alloc(1,1,6);
		moment4(TabBand[s3].buffer(),TabBand[s3].n_elem(),TabStat[s3](0,0,0), TabStat[s3](0,0,1), TabStat[s3](0,0,2),
				TabStat[s3](0,0,3), Mini, Maxi );
		TabStat[s3](0,0,4)=Mini;
		TabStat[s3](0,0,5)=Maxi;
		if(LocVerbose) if(Outname!=NULL) cstat << s3 << TabStat[s3](0,0,0) <<"\t"<<  TabStat[s3](0,0,1) <<"\t"<< TabStat[s3](0,0,2)
				<<"\t"<< TabStat[s3](0,0,3) <<"\t"<< TabStat[s3](0,0,4) <<"\t"<< TabStat[s3](0,0,5) <<endl;
		
		// Global stats
		TabStat[NbrScale3D](s3,0,0)=TabStat[s3](0,0,0);
		TabStat[NbrScale3D](s3,0,1)=TabStat[s3](0,0,1);
		TabStat[NbrScale3D](s3,0,2)=TabStat[s3](0,0,2);
		TabStat[NbrScale3D](s3,0,3)=TabStat[s3](0,0,3);
		TabStat[NbrScale3D](s3,0,4)=Mini;
		TabStat[NbrScale3D](s3,0,5)=Maxi;
	}
	
// Centering the global stats
	for(int s3=0;s3<NbrScale3D-1;s3++)
	for(int s=0;s<get_nbband(s3);s++)
	{
		int N = nbr_block(s3);
		int Ntot = get_nbplane(s3) * s2d_xn(s3,s)*s2d_yn(s3,s) ;
		moment4_center(Ntot, TabStat[NbrScale3D](s3,s,0)/N, TabStat[NbrScale3D](s3,s,1)/N, TabStat[NbrScale3D](s3,s,2)/N, TabStat[NbrScale3D](s3,s,3)/N,
				TabStat[NbrScale3D](s3,s,1), TabStat[NbrScale3D](s3,s,2), TabStat[NbrScale3D](s3,s,3));
		TabStat[NbrScale3D](s3,s,0) /= N;
		if(Verbose)
			cerr << s3 <<"\t"<< s <<"\t"<< TabStat[NbrScale3D](s3,s,0) <<"\t"<< TabStat[NbrScale3D](s3,s,1) <<"\t"<< TabStat[NbrScale3D](s3,s,2) <<"\t"<<
					 TabStat[NbrScale3D](s3,s,3) <<"\t"<< TabStat[NbrScale3D](s3,s,4) << "\t" << TabStat[NbrScale3D](s3,s,5) << endl;
		if(Outname!=NULL)
			cstat << s3 <<"\t"<< s <<"\t"<< TabStat[NbrScale3D](s3,s,0) <<"\t"<< TabStat[NbrScale3D](s3,s,1) <<"\t"<< TabStat[NbrScale3D](s3,s,2) <<"\t"<<
					 TabStat[NbrScale3D](s3,s,3) <<"\t"<< TabStat[NbrScale3D](s3,s,4) << "\t" << TabStat[NbrScale3D](s3,s,5) << endl;
	}
	
// Write the complete stat file
	if(Outname!=NULL) 
	{
		char filename[64];
		for(int s3=0;s3<NbrScale3D;s3++)
		{
			sprintf(filename,"%s_stat_full_s%d.fits",Outname,s3);
			writefltarr(filename, TabStat[s3]);
		}
		sprintf(filename,"%s_stat.fits",Outname);
		writefltarr(filename, TabStat[NbrScale3D]);
		cstat.close();
	}
	
	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}



/****************************************************************************/

void BCurvelet3D::normalize(fltarray *TabBand, fltarray *TabBandNorm, bool UseCubeSigma)
{
	if(Verbose) cerr<<"BCurvelet3D::normalize..."<<endl;
	float val;
	
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		int depth = 3*TabBlockSize[s3]*TabBlockSize[s3];
		for(int s=0;s<get_nbband(s3);s++)
			for(int i = s2d_x0(s3,s); i < s2d_x0(s3,s)+s2d_xn(s3,s); i++)
				for(int j = s2d_y0(s3,s); j < s2d_y0(s3,s)+s2d_yn(s3,s); j++)
					for(int k = 0; k < get_nbplane(s3); k++)
					{
						if(UseCubeSigma) val = (TabBand[s3])(i,j,k) / (CubeSigma[s3])(i,j,k%depth);
						else val = (TabBand[s3])(i,j,k) / TabSigma(s3,s);
						(TabBandNorm[s3])(i,j,k) = val;
					}
	}

	// coarse scale : N-1
		// No normalisation at coarse scale
	
	if(Verbose) cerr<<"...End BCurvelet3D::normalize"<<endl;
}

/****************************************************************************/

void BCurvelet3D::normalize_self(fltarray *TabBand, bool inverse, bool UseCubeSigma)
{
	if(Verbose) cerr<<"BCurvelet3D::normalize_self..."<<endl;
	float val;
	
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		int depth = 3*TabBlockSize[s3]*TabBlockSize[s3];
		for(int s=0;s<get_nbband(s3);s++)
			for(int i = s2d_x0(s3,s); i < s2d_x0(s3,s)+s2d_xn(s3,s); i++)
				for(int j = s2d_y0(s3,s); j < s2d_y0(s3,s)+s2d_yn(s3,s); j++)
					for(int k = 0; k < get_nbplane(s3); k++)
						if(inverse)
						{
							if(UseCubeSigma) val = (TabBand[s3])(i,j,k) * (CubeSigma[s3])(i,j,k%depth);
							else val = (TabBand[s3])(i,j,k) * TabSigma(s3,s);
							(TabBand[s3])(i,j,k) = val;
						}
						else
						{
							if(UseCubeSigma) val = (TabBand[s3])(i,j,k) / (CubeSigma[s3])(i,j,k%depth);
							else val = (TabBand[s3])(i,j,k) / TabSigma(s3,s);
							(TabBand[s3])(i,j,k) = val;
						}
	}

	// coarse scale : N-1
		// No normalisation at coarse scale
	
	if(Verbose) cerr<<"...End BCurvelet3D::normalize_self"<<endl;
}

/****************************************************************************/

void BCurvelet3D::threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4, bool UseCubeSigma)
{
	if(Verbose) cerr<<"BCurvelet3D::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<","<<UseCubeSigma<<")"<<endl;
	
	float lvl;
	
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		int depth = 3*TabBlockSize[s3]*TabBlockSize[s3];
		for(int s=0;s<get_nbband(s3);s++)
		{
			float cnt=0,cnt2=0;
			for(int i = s2d_x0(s3,s); i < s2d_x0(s3,s)+s2d_xn(s3,s); i++)
				for(int j = s2d_y0(s3,s); j < s2d_y0(s3,s)+s2d_yn(s3,s); j++)
					for(int k = 0; k < get_nbplane(s3); k++)
					{
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						if(UseCubeSigma) lvl = SigmaNoise * Nsig * (CubeSigma[s3])(i,j,k%depth);
						else lvl = SigmaNoise * Nsig * TabSigma(s3,s);

						cnt+=1;
						if( abs((TabBand[s3])(i,j,k)) < lvl )
						{// Hard and soft put these to 0
							cnt2+=1;
							(TabBand[s3])(i,j,k)=0;
						}// soft lowers the remaining by lvl
						else if(FilterType==FT_SOFT) (TabBand[s3])(i,j,k) -= (2*int( (TabBand[s3])(i,j,k)>0 )-1)*lvl;
					}
			if(Verbose) cerr<<"Band "<<s3<<","<<s<<") n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;
		}
	}
	
	// coarse scale : N-1
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);
	
	if(Verbose) cerr<<"...End BCurvelet3D::threshold"<<endl;
}

/****************************************************************************/

int BCurvelet3D::get_BC_pos(int kx,int ky,int bx,int by,int LocalBS,int N)
{
	int Px = ( bx + kx*LocalBS ) % (2*N) ;
	int Py = ( by + ky*LocalBS ) ;
	int s,f; // slow and fast components in terms of exterior and interior loop
	
	if(Py >= N)
	{
		if(Px < N)
		{
			f = Px ;
			s = (2*N - 1 - Py)*N ;
		}
		else // Px >= N
		{
			f = 2*N - 1 - Py ;
			s = (2*N - 1 - Px)*N ;
		}
	//cerr<<" kv kh by bx -> h v y = "<<kv<<","<<kh<<","<<by<<","<<bx<<" -> "<<h<<" "<<v<<" "<<y<<"   "<<endl;
	}
	else 
	{
		if(Px < N)
		{
			f = Py ;
			s = (Px + 2*N)*N ;
		}
		else // Px >= N
		{
			f = Py ;
			s = (3*N - Px)*N - 1 ;
		}
	}
	return s+f;
}
int xxx(int kx,int ky,int bx,int by,int LocalBS,int N)
{
	int Px = ( bx + kx*LocalBS ) % (2*N) ;
	int Py = ( by + ky*LocalBS ) ;
	int s,f; // slow and fast components in terms of exterior and interior loop
	cerr<<Px<<","<<Py;
	if(Py >= N)
	{
	cerr<<">=N"<<endl;
		if(Px < N)
		{cerr<<"Px<N"<<endl;
			f = Px ;
			s = (2*N - 1 - Py)*N ;
		}
		else // Px >= N
		{cerr<<"Px>N"<<endl;
			f = 2*N - 1 - Py ;
			s = (2*N - 1 - Px)*N ;
		}
	//cerr<<" kv kh by bx -> h v y = "<<kv<<","<<kh<<","<<by<<","<<bx<<" -> "<<h<<" "<<v<<" "<<y<<"   "<<endl;
	}
	else 
	{
	cerr<<"<N"<<endl;
		if(Px < N)
		{cerr<<"Px<N"<<endl;
			f = Py ;
			s = (Px + 2*N)*N ;
		}
		else // Px >= N
		{cerr<<"Px>N"<<endl;
			f = Py ;
			s = (3*N - Px)*N - 1 ;
		}
	}
	return s+f;
}

/****************************************************************************/

void BCurvelet3D::wiener(fltarray *TabBand, float noise_lvl, int LocalBS, bool UseCubeSigma)
{
	if(Verbose) cerr<<"BCurvelet3D::wiener("<<noise_lvl<<","<<LocalBS<<","<<UseCubeSigma<<")..."<<endl;
//	bool LocVerbose = true & Verbose;
	float val;
	float noise2 = noise_lvl*noise_lvl;
	
	// Fine wavelet scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
//cerr<<" 3d scale "<<s3<<endl;							
		// Allocation
		int NbS2D = get_nbband(s3);

		int N = TabBlockSize[s3]; // for height : number of angles 
		
		// Linelet scales
		for(int s=0;s<NbS2D;s++)
		{
			int Nax = ceil(float(2*N)/float(LocalBS)); // horizontal angle
			int Nay = ceil(float(N)/float(LocalBS)); // vertical angle
			int Npx = ceil(float(s2d_xn(s3,s))/float(LocalBS)); // horizontal space position
			int Npy = ceil(float(s2d_yn(s3,s))/float(LocalBS)); // vertical space position
			//cerr<<"s N nax Nay xn yn Npx Npy = ("<<s<<","<<N<<","<<Nax<<","<<Nay<<","<<s2d_xn(s3,s)<<","<<s2d_yn(s3,s)<<","<<Npx<<","<<Npy<<")"<<endl;

			// Spatial blocks B
//			#if USE_OMP_BC // Not parallelisable, because some blocks do overlap (but slight error only)
//				#pragma omp parallel for 
//			#endif
			for(int B=0;B<nbr_block(s3);B++)
			{
				fltarray coef_wiener(s2d_xn(s3,s),s2d_yn(s3,s),3*N*N); // 3*N² angles and s2d_xn(s3,s)*s2d_yn(s3,s) positions
				
				coef_wiener.init(-2);
				int z0 = B*3*N*N ;

				// Wiener blocks (angle^2 and space^2) for mainly horizontal Beamlets
				for(int ky=0 ; ky < Nay ; ky++)
					for(int kx=0 ; kx < Nax ; kx++)
						for(int kpy = 0; kpy < Npy; kpy++)
							for(int kpx = 0; kpx < Npx; kpx++)
							{
								// cerr<<endl<<"  B ky kx kpy kpx = "<<B<<" ("<<ky<<","<<kx<<","<<kpy<<","<<kpx<<") z ";
								double sigma2 = 0.0;
								float cnt = 0;
								
							// Sigma calculation
								// Pixels in a wiener block (angle)
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										// angular position
										int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);
										
										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												// if(B==0) cerr<<bx<<" "<<by<<" "<<bpx<<" "<<bpy<<endl;
												cnt+=1;
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if(UseCubeSigma)	val = (TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) 
																	  / (CubeSigma[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z);
												else 				val = (TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) / TabSigma(s3,s);

												sigma2+=pow(val,2);
											}
									}
								float sig2 = max( 0.0, sigma2/cnt - noise2 );
								float norm = sig2 / (sig2+noise2);
								
							// Store the coef in the table
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);

										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if( coef_wiener(x,y,z) < -1 )
													coef_wiener(x,y,z) = norm;
											}
									}
							}// end wiener blocks

				// Wiener blocks (angle^2 and space^2) for mainly vertical Beamlets
				for(int ky=0 ; ky < Nay ; ky++)
					for(int kx=0 ; kx < Nay ; kx++) // Nay becuse it's a square
						for(int kpy = 0; kpy < Npy; kpy++)
							for(int kpx = 0; kpx < Npx; kpx++)
							{
								double sigma2 = 0.0;
								float cnt = 0;
								
							// Sigma calculation
								// Pixels in a wiener block (angle)
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										// angular position
										int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));
										//z=0;
										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												cnt+=1;
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if(UseCubeSigma)	val = (TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) 
																	  / (CubeSigma[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z);
												else				val = (TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) / TabSigma(s3,s);

												sigma2+=pow(val,2);
											}
									}
								float sig = sqrt(max( 0.0, sigma2/cnt - pow(noise_lvl,2) ));
								float norm = sig / (sig+noise_lvl);
							
							// Store the coef in the table
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));

										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if( coef_wiener(x,y,z) < -1 )
													coef_wiener(x,y,z) = norm;
											}
									}
							}// end wiener blocks

//		char filename[64];
//		sprintf(filename,"coef_w_s%d_b%d.fits",s,B);
//		writefltarr(filename, coef_wiener);

			// Apply the wiener coefficient
				// Wiener blocks (int angle2 and space) for mainly horizontal beamlets
				for(int ky=0 ; ky < Nay ; ky++)
					for(int kx=0 ; kx < Nax ; kx++)
						for(int kpy = 0; kpy < Npy; kpy++)
							for(int kpx = 0; kpx < Npx; kpx++)
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										// angular position
										int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);
										
										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if( coef_wiener(x,y,z) > -1 )
												{
													(TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) *= coef_wiener(x,y,z);
													coef_wiener(x,y,z) = -2;
												}
											}
									}

				// Wiener blocks (int angle2 and space) for mainly vertical bebamlets
				for(int ky=0 ; ky < Nay ; ky++)
					for(int kx=0 ; kx < Nay ; kx++)
						for(int kpy = 0; kpy < Npy; kpy++)
							for(int kpx = 0; kpx < Npx; kpx++)
								for(int by = 0; by < LocalBS; by++)
									for(int bx = 0; bx < LocalBS; bx++)
									{
										// angular position
										int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));

										// pixels in a wiener block (space)
										for(int bpy = 0; bpy < LocalBS; bpy++)
											for(int bpx = 0; bpx < LocalBS; bpx++)
											{
												// spatial position = position in the block + block_position
												int x = ((bpx + kpx*LocalBS) % s2d_xn(s3,s)) ;
												int y = ((bpy + kpy*LocalBS) % s2d_yn(s3,s)) ;
												if( coef_wiener(x,y,z) > -1 )
												{
													(TabBand[s3])(x+s2d_x0(s3,s), y+s2d_y0(s3,s), z+z0) *= coef_wiener(x,y,z);
													coef_wiener(x,y,z) = -2;
												}
											}
									}
			}// end physical block
		}// end 2D scale
	}// end 3D scale
			
	// Coarse scale
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);

	if(Verbose) cerr<<"...End BCurvelet3D::wiener"<<endl;	
}

/****************************************************************************/

void BCurvelet3D::fdr(fltarray * &TabBand, float Alpha, float SigmaNoise)
{// Does not take --usecubsigma into account
	
	if(Verbose) cerr<<"BCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s=0;s<get_nbband(s3);s++)
		{
			
		//getband
			fltarray Band;
			get_band(s3, s, TabBand, Band);

		// Cumulative of the band
			CArrProb *TabCP;
			TabCP = new CArrProb;
			TabCP->set(Band,false);
			
			float lvl = TabSigma(s3,s)*SigmaNoise;
			
		// P-Values
			fltarray PVal(Band.nx(),Band.ny(),Band.nz());
			for (int i=0; i < Band.nx(); i++)
				for (int j=0; j < Band.ny(); j++)
					for (int k=0; k < Band.nz(); k++)
						PVal(i,j,k) = prob_noise(Band(i,j,k)/lvl);
	
			double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);
			
			cerr<<" band ("<<s3<<","<<s<<") : PDet="<<PDet<<endl;

		// threshold
			int i0=s2d_x0(s3,s);
			int j0=s2d_y0(s3,s);
			for(int i = 0; i < s2d_xn(s3,s); i++)
				for(int j = 0; j < s2d_yn(s3,s); j++)
					for(int k = 0; k < get_nbplane(s3); k++)
						if( PVal(i,j,k) > PDet) (TabBand[s3])(i+i0, j+j0, k) = 0;
			
			delete TabCP;
		}
	
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);

	if(Verbose) cerr<<"...End BCurvelet3D::fdr"<<endl;
}

/****************************************************************************/

// writes the values of the transform at positions givent in filename
void BCurvelet3D::values_at(fltarray *TabBand, char * filename, char* Outname)
{
	int S3,S2,x,y,z;
	float m;
	
	// Output stream
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_values_at.dat");
	fstream cval;
	cval.open (MaxCoefname, fstream::out);
	
	FILE* fic;
	fic=fopen(filename, "r");
	if(fic)
	{
		for(int s3=0;s3<NbrScale3D-1;s3++)
		{
		// the 3 fine ones at each 2d scale
			for(int s2=0; s2<get_nbscale2d(s3)-1; s2++)
			{
				if(!fscanf(fic,"%d %d %f %d %d %d",&S3,&S2,&m,&x,&y,&z)) cerr<<"Error while reading coordinates in "<<filename<<endl;
				cval<<s3<<"\t"<<3*s2<<"\t"<<abs((TabBand[s3])(x,y,z))<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
				if(!fscanf(fic,"%d %d %f %d %d %d",&S3,&S2,&m,&x,&y,&z)) cerr<<"Error while reading coordinates in "<<filename<<endl;
				cval<<s3<<"\t"<<3*s2+1<<"\t"<<abs((TabBand[s3])(x,y,z))<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
				if(!fscanf(fic,"%d %d %f %d %d %d",&S3,&S2,&m,&x,&y,&z)) cerr<<"Error while reading coordinates in "<<filename<<endl;
				cval<<s3<<"\t"<<3*s2+2<<"\t"<<abs((TabBand[s3])(x,y,z))<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
			}
		// the coarse one
			int s2=get_nbband(s3)-1; 
			if(!fscanf(fic,"%d %d %f %d %d %d",&S3,&S2,&m,&x,&y,&z)) cerr<<"Error while reading coordinates in "<<filename<<endl;
			cval<<s3<<"\t"<<s2<<"\t"<<abs((TabBand[s3])(x,y,z))<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
		}
		fclose(fic);
	}
	else cerr<<"Warnig: File "<<filename<<" not found, for use in BCurvelet3D::values_at"<<endl;
	cval.close();
}

/****************************************************************************/

void BCurvelet3D::write (char *Name, fltarray * & TabBand, bool Normalize)
{
	if(Verbose) cerr<<"BCurvelet3D::write("<<Name<<",.,"<<Normalize<<")"<<endl;
	
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1; 

	// .mr extention
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 
	
// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;

// write first header part (parameters)
	naxis=0;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D BeamCurvelet", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale3D, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of 3d blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Size of 3d blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"TypeW3D", (long) TypeW3D, (char*)"Type of wavelets used", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"CurTrans", (long) CurTrans, (char*)"Type of BCurvelet Transform(Beamlet type)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) NxCube, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) NyCube, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) NzCube, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	for (int s=0; s < NbrScale3D; s++)
	{
		naxis=3;
		naxes[0] = TabBand[s].nx();
		naxes[1] = TabBand[s].ny();
		naxes[2] = TabBand[s].nz();
//cerr<<"save : "<<s<<" "<<TabBand[s].nx()<<" "<<TabBand[s].ny()<<" "<<TabBand[s].nz()<<endl;
		if(ffcrhd(fptr,&status))
			PrintError( status );
		if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
			PrintError( status );
		if ( ffpkyj(fptr, (char*)"Nx", (long) naxes[0], (char*)"x size of the band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny", (long)  naxes[1], (char*)"y size of the band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz", (long)  naxes[2], (char*)"z size of the band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nxo", (long) TabSizeW(s,0,0), (char*)"x size of the wavelet scale", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nyo", (long) TabSizeW(s,1,0), (char*)"y size of the wavelet scale", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nzo", (long) TabSizeW(s,2,0), (char*)"z size of the wavelet scale", &status))
			PrintError( status );  
		
	// save the data
		if ( ffppre(fptr, group, 1, TabBand[s].nx()*TabBand[s].ny()*TabBand[s].nz(), (TabBand[s]).buffer(), &status) )
			PrintError( status );
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end BCurvelet3D::write"<<endl;
}

/*********************************************************************/

void BCurvelet3D::read(char *Name, fltarray * & TabBand, SubBandFilter * SB1D, bool *Normalize)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	char comment[FLEN_COMMENT];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );

// get number of Headers
	int nhead;
	fits_get_num_hdus(fptr, &nhead, &status);
	
// read primary header
	if ( ffmahd(fptr, 1, &hdutype, &status) ) PrintError( status );
	
	// Read params
	int _NbrScale,_BlockSize,_Nx,_Ny,_Nz;
	float _BlockOverlap;
	type_linelet3d_WTtrans _CurTrans;
	type_wavelet3D _TypeW3D;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*Normalize = (bool)mon_long;
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	_BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	_BlockOverlap = (float)mon_long;
	if (ffgkyj(fptr,(char*)"CurTrans", &mon_long, comment, &status)) PrintError( status );
	_CurTrans = (type_linelet3d_WTtrans)mon_long;
	if (ffgkyj(fptr,(char*)"TypeW3D", &mon_long, comment, &status)) PrintError( status );
	_TypeW3D = (type_wavelet3D)mon_long;
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
//	init the structure and the output vector
	CurTrans=_CurTrans;
	BlockOverlap=_BlockOverlap;
	init(_Nx, _Ny, _Nz, _NbrScale, _BlockSize, _TypeW3D, SB1D);
	alloc();
	TabBand = new fltarray[NbrScale3D];
//TESTcomment
//	for(int s=0; s < NbrScale3D-1; s++)
//		TabBand[s].alloc(TabBeamlet[s].get_BlockSize(), TabBeamlet[s].get_BlockSize(), get_nbplane(s));
	
//	TabWaveBand = new fltarray [NbrScale3D];
	
// read data
	for(int s=0;s<_NbrScale;s++)
	{
//cerr<<"save : "<<s<<" "<<TabBeamlet[s].get_BlockSize()<<" "<< TabBeamlet[s].get_BlockSize()<<" "<< get_nbplane(s3) <<endl;
		int NX,NY,NZ,NXo,NYo,NZo;
		if (fits_movabs_hdu(fptr, s+2, NULL, &status)) PrintError( status );
		
		if (ffgkyj(fptr,(char*)"Nx", &mon_long, comment, &status)) PrintError( status );
		NX = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Ny", &mon_long, comment, &status)) PrintError( status );
		NY = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nz", &mon_long, comment, &status)) PrintError( status );
		NZ = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nxo", &mon_long, comment, &status)) PrintError( status );
		NXo = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nyo", &mon_long, comment, &status)) PrintError( status );
		NYo = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nzo", &mon_long, comment, &status)) PrintError( status );
		NZo = (int)mon_long;
		
		if(Verbose) cerr<<" read TB("<<s<<") : "<<NX<<" "<<NY<<" "<<NZ;
		if(Verbose) if(s<_NbrScale-1) cerr<<" bs "<<TabBlockSize[s];
		if(Verbose) cerr<<endl;
		TabBand[s].alloc(NX,NY,NZ);
		if (ffgpve(fptr, 1, 1, NX*NY*NZ, nulval, (TabBand[s]).buffer(), &anynul, &status)) PrintError( status );
		if(Verbose) cerr<<" read TWB("<<s<<"): "<<NXo<<" "<<NYo<<" "<<NZo<<endl;

		if(s<_NbrScale-1)
		{
			TabBeamlet[s].set_attribut (NXo,NYo,NZo, TabBlockSize[s]);
			TabBeamlet[s].set_block_param();
		}
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}

/****************************************************************************/

void BCurvelet3D::temp(fltarray *TabBand)
{
/*	for(int s3=0;s3<NbrScale3D;s3++)
	{
		for(int i=0; i<TabBand[s3].nx(); i++)
			for(int j=0; j<TabBand[s3].ny(); j++)
				for(int k=0; k<TabBand[s3].nz(); k++)
					(TabBand[s3])(i,j,k)=0;
	}
*/
			
/*	cerr<<"nb planes/blocks="<<get_nbplane(1)<<"/"<<nbr_block(1)<<",nscale2d="<<get_nbscale2d(1)<<endl;
	cerr<<s2d_x0(1,2)<<"x"<<s2d_y0(1,2)<<endl;
	cerr<<s2d_x0(1,1)<<"x"<<s2d_y0(1,1)<<endl;
	cerr<<s2d_x0(1,0)<<"x"<<s2d_y0(1,0)<<endl;
//	(TabBand[1])(6,3,27444)=1; // 1 dir
	(TabBand[1])(7,7,18640)=1; // diag
*/
//	for(int x=0;x<NbrScale3D-1;x++) cerr<<get_nbscale2d(x)<<endl;
//	int s3,s; float p,q,a;
//	s3=0;s=6;p=0.5,q=0.5, a=0.62;
//	TabBand[s3](s2d_x0(s3,s)+int(s2d_xn(s3,s)*p), s2d_y0(s3,s)+int(s2d_yn(s3,s)*q), int(a*get_nbradplane(s3))+get_nbradplane(s3)*nbr_block(s3)/2)=1;//num_block(s3,nbr_block_nx(s3)/2,nbr_block_ny(s3)/2,nbr_block_nz(s3)/2))=1;
//	s3=0;s=6;p=0.5,q=0.5, a=0.5;
//	TabBand[s3](s2d_x0(s3,s)+int(s2d_xn(s3,s)*p), s2d_y0(s3,s)+int(s2d_yn(s3,s)*q), int(a*get_nbradplane(s3))+get_nbradplane(s3)*num_block(s3,nbr_block(s3)/2,nbr_block(s3)/2,nbr_block(s3)/2))=1;
//	s3=0;s=0;p=0.5,q=0.5, a=0.5;
//	TabBand[s3](s2d_x0(s3,s)+int(s2d_xn(s3,s)*p), s2d_y0(s3,s)+int(s2d_yn(s3,s)*q), int(a*get_nbradplane(s3))+get_nbradplane(s3)*num_block(s3,nbr_block_nx(s3)/2,nbr_block_ny(s3)/2,nbr_block_nz(s3)/2))=1;
}



void bcur_clear(vector< fltarray* > &C)
{// corresponds to BFCurvelet3D::transform and FCurvelet3D::cur_trans
	for(int i=0;i<(int)C.size();i++)
		delete C[i];
	C.clear();
}	



void bcur_transform(fltarray &Data, vector< fltarray* > &vTabBand, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D, type_pmeyer3D Type_PMW3D)
{
//	Verbose = True;
	if(Verbose) cerr << "bcur_transform..." <<  endl;
	
	fltarray *TabBand;
	
// BCurvelet initialisation	
	BCurvelet3D *DataC = new BCurvelet3D();
	DataC->CurTrans=DEF_LIN3D_TRANS;
	DataC->Type_PMW3D=Type_PMW3D;
	DataC->BlockOverlap=BlockOverlap;
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	DataC->init(Data.nx(), Data.ny(), Data.nz(), NbrScale3D, BlockSize, W3D, SB1D);
	
	DataC->transform(Data, TabBand);
	DataC->normalize_self(TabBand);
	
	vTabBand.resize(NbrScale3D);

	for(int s=0;s<NbrScale3D;s++)
	{
		vTabBand[s] = new fltarray;
		vTabBand[s]->alloc(TabBand[s].nx(),TabBand[s].ny(),TabBand[s].nz());
		float *buf=TabBand[s].buffer();
		float *vbuf=vTabBand[s]->buffer();
		for(int i=0;i<TabBand[s].n_elem();i++)
			vbuf[i] = buf[i];
	}
	
	delete [] TabBand;
	delete DataC;
	delete SB1D;
	if(Verbose) cerr << "...End bcur_transform" <<  endl;
}





void bcur_recons(vector< fltarray* > &vTabBand, fltarray &Recons, int Nx, int Ny, int Nz, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D, type_pmeyer3D Type_PMW3D)
{
//	Verbose = True;
	if(Verbose) cerr << "bcur_recons..." <<  endl;
	
// Load data
	fltarray *TabBand = new fltarray[NbrScale3D];
	for(int s=0;s<NbrScale3D;s++)
		TabBand[s].alloc(vTabBand[s]->buffer(),vTabBand[s]->nx(),vTabBand[s]->ny(),vTabBand[s]->nz());
	
// Initialize BCurvelet3D
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	BCurvelet3D *DataC = new BCurvelet3D();
	DataC->CurTrans=DEF_LIN3D_TRANS;
	DataC->Type_PMW3D=Type_PMW3D;
	DataC->BlockOverlap=BlockOverlap;
	DataC->init(Nx, Ny, Nz, NbrScale3D, BlockSize, W3D, SB1D);
	DataC->alloc();
	double DNX = Nx;
	double DNY = Ny;
	double DNZ = Nz;
	int NXo, NYo, NZo;

	for(int s=0;s<NbrScale3D;s++)
	{
		NXo = int( DNX + 0.5);
		NYo = int( DNY + 0.5);
		NZo = int( DNZ + 0.5);
		
		if(s<NbrScale3D-1)
		{
			DataC->TabBeamlet[s].set_attribut(NXo,NYo,NZo, DataC->TabBlockSize[s]);
			DataC->TabBeamlet[s].set_block_param();
		}

		DNX /= 2.;
		DNY /= 2.;
		DNZ /= 2.;
	}

// inverse normalization
	DataC->normalize_self(TabBand, true); 
	
// Reconstruction
	DataC->recons(TabBand, Recons, true);
	
	delete [] TabBand;
	delete DataC;
	delete SB1D;
	if(Verbose) cerr << "...End bcur_recons" <<  endl;
}



void bcur_threshold(vector< fltarray* > &vTabBand, float threshold, filter_type FilterType)
{
	for(int s=0;s<(int)vTabBand.size()-1;s++)
	{
		float* buf=vTabBand[s]->buffer();
		if(FilterType==FT_HARD)
		{
			for(int i=0;i<vTabBand[s]->n_elem();i++)
				if(abs(buf[i])<threshold) buf[i] = 0;
		}
		else if(FilterType==FT_SOFT)
		{
			for(int i=0;i<vTabBand[s]->n_elem();i++)
				buf[i] = ::soft_threshold(buf[i],threshold);
		}
		else cerr<<"Filtering method not implemented yet."<<endl;
	}
}


void bcur_filter(fltarray &Data, fltarray &recons, int NbrScale3D, int BlockSize, float BlockOverlap, float lvl, filter_type FilterType, type_wavelet3D W3D, type_pmeyer3D Type_PMW3D)
{
//	Verbose = false;
	if(Verbose) cerr << "bcur_filter..." <<  endl;
	
// BCurvelet initialisation	
	BCurvelet3D *DataC = new BCurvelet3D();
	DataC->CurTrans=DEF_LIN3D_TRANS;
	DataC->Type_PMW3D=Type_PMW3D;
	DataC->BlockOverlap=BlockOverlap;
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	DataC->init(Data.nx(), Data.ny(), Data.nz(), NbrScale3D, BlockSize, W3D, SB1D);
	
	DataC->BCurvelet3D::filter(Data, recons, lvl, 1, FilterType, 0,3,false, false, NULL);
	
	delete SB1D;
	delete DataC;
}



