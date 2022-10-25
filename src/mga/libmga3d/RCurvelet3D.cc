/******************************************************************************
**                   Copyright (C) 20XX by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: 
**
**    Date:  
**    
**    File:  
**
**    Modification history:
**
******************************************************************************
**
**    DESCRIPTION RidCurvelet transform and reconstruction
**    -----------  
**                 
******************************************************************************/

#include "RCurvelet3D.h"

extern Bool Verbose;

/****************************************************************************/

// Scale limits : for each blocksize, the maximum ridgelet scale where the 
// coefficients are significant
static inline double prob_noise(float Val)
{
	return (abs(Val) < FLOAT_EPSILON) ? 1.:erfc (abs(Val));
}
static int NScale1d_meyer_bs3d[2][4] = // BS = 17 33
	{{	4,4,3,2},
	{	5,5,4,2}};

// Noise lvl in bands
// see estim_norm_cur1g.pro for details on the coefficients
	static float TabRCur_Meyer_17[4][4] = 
		{{	13.866721, 1.2610560  , 0.57451532, 1.24029260	},
		{	16.756738, 1.7641779  , 0.77306028, 1.8684821	},
		{	8.7640957, 1.0544143  , 1.3563040 , 0.0000000	},
		{	9.5847841, 0.027607206, 0.000000  , 0.0000000	}};
	
	static float TabRCur_Meyer_33[4][5] = // wrong : divide by TabMeyer or recalculate
		{{	25.8605, 1.77663, 0.973708, 0.483488, 14.2408 },
		{	19.2258, 1.55780, 0.822433, 0.407946, 10.5549 },
		{	9.91600, 1.40789, 0.663052, 1.38914, 0.00000 },
		{	9.91979, 0.432196, 0.00000, 0.00000, 0.00000 }};

	static float TabRCur_Meyer_17_3sig[4][4] = // wrong : divide by TabMeyer or recalculate
		{{	13.7017, 1.31611, 0.575103, 1.05334 },
		{	9.81526, 1.12459, 0.474165, 0.586913 },
		{	5.17643, 1.16421, 1.00521, 0.00000 },
		{	4.92630, 0.0161400, 0.00000, 0.00000 }};

	static float TabRCur_Meyer_33_3sig[4][5] = // wrong : divide by TabMeyer or recalculate
		{{	26.2396, 1.82170, 0.988166, 0.484695, 22.8454 },
		{	19.4863, 1.59921, 0.817849, 0.394834, 16.9041 },
		{	10.0055, 1.49199, 0.697501, 0.876224, 0.00000 },
		{	9.93847, 0.439465, 0.00000, 0.00000, 0.00000 }};

	
	static float TabRCur_Atrou_17[4][4] = // on 64
		{{	17.7975, 2.12497, 0.789833,  1.09914}, 
		{	6.07268, 3.25467,  1.33222, 0.926225}, 
		{	1.74982, 2.41491,  1.87934,  1.33824}, 
		{	1.05984, 2.00479,  3.03320,  3.15215}};    

	static float TabRCur_Atrou_33[4][5] = // on 64
		{{	33.8855, 3.52858, 1.30910, 0.503704,  31.5298}, 
		{	12.0137, 6.42699, 2.27989, 0.738415, 0.488817}, 
		{	3.35609, 4.86843, 3.61662,  1.47286, 0.621251}, 
		{	1.61813, 4.70628, 8.33424,  6.72083,  3.71307}};

	static float TabRCur_Atrou_17_3sig[4][4] = // WRONG : copy TabRCur_Atrou_17
		{{	17.7975, 2.12497, 0.789833,  1.09914}, 
		{	6.07268, 3.25467,  1.33222, 0.926225}, 
		{	1.74982, 2.41491,  1.87934,  1.33824}, 
		{	1.05984, 2.00479,  3.03320,  3.15215}};    

	static float TabRCur_Atrou_33_3sig[4][5] = // WRONG copy TabRCur_Atrou_33
		{{	33.8855, 3.52858, 1.30910, 0.503704,  31.5298}, 
		{	12.0137, 6.42699, 2.27989, 0.738415, 0.488817}, 
		{	3.35609, 4.86843, 3.61662,  1.47286, 0.621251}, 
		{	1.61813, 4.70628, 8.33424,  6.72083,  3.71307}};

	
	static float TabRCur_PMeyer_17[4][4] = 
		{{	13.771378, 1.2791816, 0.59828704, 1.2583345},
		{	17.024107, 1.8568499, 0.84879462, 1.7062063},
		{	9.0080841, 1.0083679  , 1.1489091, 0.0000000},
		{	9.7039889, 0.028214873, 0.0000000, 0.0000000}};
	
	static float TabRCur_PMeyer_33[4][5] =
		{{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1}};

	static float TabRCur_PMeyer_17_3sig[4][4] = // copy from TabRCur_PMeyer_17[4][4]
		{{	13.771378, 1.2791816, 0.59828704, 1.2583345},
		{	17.024107, 1.8568499, 0.84879462, 1.7062063},
		{	9.0080841, 1.0083679  , 1.1489091, 0.0000000},
		{	9.7039889, 0.028214873, 0.0000000, 0.0000000}};

	static float TabRCur_PMeyer_33_3sig[4][5] = 
		{{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1},
		{	1, 1, 1, 1, 1}};


/****************************************************************************/

RCurvelet3D::RCurvelet3D()
{
	Ptr_SB1D=NULL;
	CurNl=NULL;
	CurNc=NULL;
	ScaleNx=NULL;
	ScaleNy=NULL;
	ScaleNz=NULL;
	TabBlockSize=NULL;
	TabRidgelet=NULL;
	CubeSigma=NULL;
	AllocClass=False;
	GetAutoNbScale=False;
	InitClass=false;
	no_coarse=false;
}

RCurvelet3D::~RCurvelet3D()
{
	dealloc();
	deinit();
}

void RCurvelet3D::dealloc()
{
	if(CurNl!=NULL)		delete [] CurNl;
	if(CurNc!=NULL)		delete [] CurNc;
	if(ScaleNx!=NULL)	delete [] ScaleNx;
	if(ScaleNy!=NULL)	delete [] ScaleNy;
	if(ScaleNz!=NULL)	delete [] ScaleNz;
	if(TabRidgelet!=NULL) delete [] TabRidgelet;
	AllocClass=False;
}

void RCurvelet3D::deinit()
{
	if(CubeSigma!=NULL) {for(int i=0;i<4;i++) delete [] CubeSigma[i]; delete [] CubeSigma;}
	if(TabBlockSize!=NULL) delete [] TabBlockSize;
	InitClass=false;
}

/****************************************************************************/

void RCurvelet3D::calib_noise_nsig(fltarray * TabBand, float N_Sigma, char* Name_Imag_Out)
{
	bool LocVerbose = false & Verbose;
	double LMax;
	fltarray Band;

	TabSigma.init(0); // useless when there is no bug

	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
		{
			// suppressed, because automatic NSigma+1 at lowest scale in thresholding, unless specified
			// float Nsig = (s3 == 0) ? N_Sigma+1: N_Sigma;
			float Nsig=N_Sigma;
			get_band(s3,s1,TabBand, Band);

			CArrProb *TabCP;
			TabCP = new CArrProb;

			TabCP->set(Band,true);
			TabCP->find_gthreshold(Nsig, LMax);
			TabSigma(s3,s1) = LMax / N_Sigma;
			if(LocVerbose) cerr << "Band (" << s3<<","<<s1 << ") Lmax = " << TabSigma(s3,s1)<< endl;

			delete TabCP;
		}

	char filename[64];
	sprintf(filename,"%s_nsig.fits",Name_Imag_Out);
	writefltarr(filename, TabSigma);
}

/****************************************************************************/

void RCurvelet3D::tab_block_size_init()
{
	int s, B;

	TabBlockSize = new int[NbrScale3D-1];

// Block size must be odd for good ridgelet reconstruction
	B = (int) next_power_of_2((int) BlockSize-1)+1;

	switch (TypeBlock)
	{
		case BL3D_CONST: 
			for (s=0; s < NbrScale3D-1; s++)
				TabBlockSize[s] = B;
			break;
		case BL3D_UP:
			for (s=0; s < NbrScale3D-1; s++)
			{
				if(B%2==0) B-=1;
				TabBlockSize[s] = B;
				B *= 2;
			}
			break;
		case BL3D_UP2:
			for (s=0; s < NbrScale3D-1; s++)
			{
				if(B%2==0) B-=1;
				TabBlockSize[s] = B;
				if (s % 2 == 0) B *= 2;
			}
			break;
		case BL3D_DOWN:
			for (s=0; s < NbrScale3D-1; s++)
			{
				if(B%2==0) B+=1;
				TabBlockSize[s] = B;
				B /= 2;
				if (B < 9) B = 9;
			}
			break;
		case BL3D_DOWN2:
			for (s=0; s < NbrScale3D-1; s++)
			{
				if(B%2==0) B+=1;
				TabBlockSize[s] = B;
				if (s % 2 != 0) 
				{
					B /= 2;
					if (B < 9) B = 9;
				}
			}
			break;
	}
}

/****************************************************************************/

void RCurvelet3D::wavelet_transform(fltarray & Cube, fltarray * & TabWaveBand, bool alloc)
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
		if(Verbose) cerr << " FFT based Meyer WT, size ="<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		MEYER_WT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.transform(Cube,TabWaveBand, (Bool)alloc);
		DataW.normalize_self(TabWaveBand);
	}
	else if(TypeW3D==W3D_PMEYER)
	{
		if(Verbose) cerr << " FFT based Poisson Meyer WT, size ="<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		POISSON_MWT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.transform(Cube,TabWaveBand, (Bool)alloc);
		DataW.normalize_self(TabWaveBand);
	}
	else
	{
		cerr<< "Undefined wavelet3D transform"<<endl;
		exit(-1);
	}
}

/****************************************************************************/

void RCurvelet3D::reset()
{
	NxCube=0;
	NyCube=0;
	NzCube=0;
	BlockSize=0;
	NbrScale1D=0;
	NbrScale3D=0;
	TypeW3D=(type_wavelet3D)0;
	Border=I_MIRROR;
	TypeBlock=(type_3dcurvelet_block)0;

 //pas en entrée de init
	//CurTrans=(type_ridgelet3d_WTtrans)0;
	//BlockOverlap=0;
}

/****************************************************************************/

void RCurvelet3D::init(int Nx, int Ny, int Nz, int _NbrScale3D, int _NbrScale1D, int BS, type_wavelet3D _TypeW3D, bool use3sig)
{
	if (Verbose) cout << "RCurvelet3D::init("<<Nx<<","<<Ny<<","<<Nz<<","<<_NbrScale3D<<","<<_NbrScale1D<<","<<BS<<","<<_TypeW3D<<","<<use3sig<<")" << endl;

	if(InitClass) dealloc();
	InitClass=true;
	reset();

// Set parameters / Default parameters
	NxCube		= Nx;
	NyCube		= Ny;
	NzCube		= Nz;
	BlockSize	= BS > 0 ? BS : DEF_3DCUR_BLOCK_SIZE;
// Constrain the BlockSize to two values
	if(BlockSize>25) BlockSize = 33;
	else BlockSize = 17;
	NbrScale1D	= _NbrScale1D > 0 ? _NbrScale1D : (int)floor(log2(BlockSize));
	NbrScale3D	= _NbrScale3D > 0 ? _NbrScale3D : DEF_3DCUR_NBR_SCALE;
	TypeW3D		= _TypeW3D;
	Use_3sigma	= use3sig;

	switch(TypeW3D)
	{
		case W3D_MEYER:
		case W3D_PMEYER:
			TypeBlock = BL3D_DOWN2;
			break;

		case W3D_ATROU:
			TypeBlock = BL3D_UP2;
			break;

		default:
			TypeBlock = DEF_3DCUR_TYPE_BLOCK;
			break;
	}
	if(CurTrans<=0 || CurTrans > NBR_RID3D_TRANS) CurTrans = DEF_RID3D_TRANS;
	if(TypeW3D<=0 || TypeW3D > NBR_TYPE_W3D) TypeW3D = DEF_TYPE_W3D;
	if(CurTrans == RID3D_OWT)
	{
		cerr<<endl<<"Ridgelet3D type : RID3D_CL_ORTHO not supported yet"<<endl;
		exit(-1);
	}

// Block Sizes calculus
	tab_block_size_init();

// TabSigma initialisation
	if(Verbose) cerr<<"TabSigma initialisation, Nscale1Dmax="<<NbrScale1D<<endl;
	TabSigma.resize(NbrScale3D-1,NbrScale1D);
	if(!use3sig)
	{
		if (TypeW3D==W3D_ATROU)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Atrou_33[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Atrou_17[i][j];
		}
		else if(TypeW3D==W3D_PMEYER)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_PMeyer_33[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_PMeyer_17[i][j];
		}
		else // W3D_MEYER and default
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Meyer_33[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Meyer_17[i][j];
		}
	}
	else // use 3sigma cumulative estimation
	{
		if (TypeW3D==W3D_ATROU)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Atrou_33_3sig[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Atrou_17_3sig[i][j];
		}
		else if(TypeW3D==W3D_PMEYER)
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_PMeyer_33_3sig[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_PMeyer_17_3sig[i][j];
		}
		else // W3D_MEYER and default
		{
			if (BlockSize==33)
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Meyer_33_3sig[i][j];
			else //(BlockSize==17) and default
				for(int i=0;i<NbrScale3D-1;i++)
					for(int j=0;j<NbrScale1D;j++)
						TabSigma(i,j)=TabRCur_Meyer_17_3sig[i][j];
		}
	}
	
// CubeSigma initialisation
	char filename[128];
	CubeSigma = new fltarray*[4];
	int bbs=BlockSize/16-1;
	for(int s3=0;s3<4;s3++)
	{
		int ns1d = min(NbrScale1D,NScale1d_meyer_bs3d[bbs][s3]);
		//cerr<<" scale "<<s3<<" : "<<NScale1d_meyer_bs3d[bbs][s3]<<endl;
		CubeSigma[s3] = new fltarray[ns1d];
		for(int s=0;s<ns1d;s++)
		{
			sprintf(filename,"%s/common/RC_cubesigma/BS%d/b%d%d_mean.fits",getenv("CUR"),BlockSize,s3,s);
			if(access(filename, F_OK)==0)
				fits_read_fltarr(filename, CubeSigma[s3][s]);
			else 
			{
				if(Verbose) cerr<<"Sigma map "<<filename<<" not found."<<endl;
				// allocate cubesigma with the correct size, set to 1
			}
		}
			
	}
	if (Verbose) cerr << "...End RCurvelet3D::init" << endl;
}

/****************************************************************************/

void RCurvelet3D::alloc(fltarray * & TabWaveBand)
{
	bool LocVerbose=true & Verbose;
	if (Verbose) cout << "RCurvelet3D::alloc..." << endl;

	if(AllocClass) dealloc();

	int s;

// Ridgelet initialisation
	if (LocVerbose) cerr << " Ridgelet init" << endl;
	TabRidgelet = new Ridgelet3D[NbrScale3D-1];
	for (s = 0; s <  NbrScale3D-1; s++)
	{
		ridgelet_init(TabRidgelet[s],TabBlockSize[s]);
		int n1dmax = 5;
		if(TypeW3D==W3D_MEYER || TypeW3D==W3D_PMEYER)
		{
			if(BlockSize==17) n1dmax=NScale1d_meyer_bs3d[0][s];
			else if(BlockSize==33) n1dmax=NScale1d_meyer_bs3d[1][s];
			else cerr<<"Wrong BlockSize"<<endl;
		}
		
		TabRidgelet[s].NbrScale = min( TabRidgelet[s].NbrScale, n1dmax);
		int Nmax = max(max(TabWaveBand[s].nx(),TabWaveBand[s].ny()),TabWaveBand[s].nz());
		if(TypeW3D==W3D_MEYER || TypeW3D==W3D_PMEYER)
			if(TabBlockSize[s] >= Nmax )// if one block at this scale, then two scales only, as more will have no energy due to the singularity of the transform (in Meyer modes)
			TabRidgelet[s].NbrScale=2;
		//cerr<<"after ridinit s"<<s<<", BS="<<BlockSize<<", Ns1d="<<TabRidgelet[s].NbrScale<<endl;
	}		
	AllocClass = True;

// Curvelet initialisation
	CurNl = new int [NbrScale3D-1];
	CurNc = new int [NbrScale3D-1];
	ScaleNx = new int [NbrScale3D];
	ScaleNy = new int [NbrScale3D];
	ScaleNz = new int [NbrScale3D];

	for (s=0; s < NbrScale3D-1; s++)
	{ 
		Block3D B3D;
		//B3D.BlockOverlap = BlockOverlap;
		B3D.BlockOverlap = (BlockOverlap>0)?True:False;
		B3D.Overlap = BlockOverlap;
		B3D.alloc(TabWaveBand[s].nx(),TabWaveBand[s].ny(),TabWaveBand[s].nz(),TabBlockSize[s]);
//		if (LocVerbose) cerr<<"  TabBlockSize["<<s<<"]="<<TabBlockSize[s]<<" B3D.BlkSize="<<B3D.block_size()<<endl;//NOTE
		TabBlockSize[s] = B3D.block_size();

		CurNl[s] = 3 * TabBlockSize[s] * TabBlockSize[s]; 
		CurNc[s] = TabBlockSize[s] * B3D.nbr_block();
		ScaleNx[s]=TabWaveBand[s].nx();
		ScaleNy[s]=TabWaveBand[s].ny();
		ScaleNz[s]=TabWaveBand[s].nz();
    	if (rid3d_class(TabRidgelet[s].RidTrans) == RID3D_CL_PYR) 
        	if ((GetAutoNbScale == True) || (NbrScale1D > 1))
        	   CurNc[s] *=2;

		if (LocVerbose)
			cerr<< " Scale "<<s<<" : CurNl CurNc "<<CurNl[s]<<" "<<CurNc[s]<<", nb blocs "<< B3D.nbr_block_nx()<<
					", blksize "<<B3D.block_size()<<", and Nscale1D="<<TabRidgelet[s].NbrScale<<endl;
	}

	{
		int s=NbrScale3D-1;
		ScaleNx[s]=TabWaveBand[s].nx();
		ScaleNy[s]=TabWaveBand[s].ny();
		ScaleNz[s]=TabWaveBand[s].nz();
	}
	if (Verbose) cerr << "...End RCurvelet3D::alloc" << endl;
}


/****************************************************************************/

// Initialize the ridgelet transform
void RCurvelet3D::ridgelet_init(Ridgelet3D &Rid, int BS)
{
	if (Ptr_SB1D != NULL) Rid.reset(Ptr_SB1D);
	Rid.BlockOverlap = (Bool)(BlockOverlap>0);
	Rid.OverlapFactor = BlockOverlap;
	Rid.RidTrans = CurTrans;
	Rid.Verbose = Verbose;
	Rid.NbrScale = min(NbrScale1D,(int)floor(log2(BS)));
	//cerr<<"? rid init ns="<<Rid.NbrScale<<",BS="<<BS<<endl;
	if(NbrScale1D>0) Rid.GetAutoNbScale = GetAutoNbScale;
	else Rid.GetAutoNbScale = True;
}

/****************************************************************************/

void RCurvelet3D::transform(fltarray &Cube, fltarray* & TabBand, bool TBalloc)
{
	if(Verbose) cerr << "RCurvelet3D::transform..." <<  endl;

// Wavelet transform of the Cube
	fltarray* TabWaveBand;
	wavelet_transform(Cube, TabWaveBand, true);

// Multiscale ridgelet transform
	transform(TabWaveBand, TabBand, TBalloc);

	delete[] TabWaveBand;

	if(Verbose) cerr << "...End RCurvelet3D::transform" <<  endl;
}

/****************************************************************************/

void RCurvelet3D::transform(fltarray * & TabWaveBand, fltarray* & TabBand, bool TBalloc)
{
	bool LocVerbose=False & Verbose;
	if(Verbose) cerr << "RCurvelet3D::transform_band..." <<  endl;

// Dimentions of the coarsest scale
	int Nx = TabWaveBand[NbrScale3D-1].nx();
	int Ny = TabWaveBand[NbrScale3D-1].ny();
	int Nz = TabWaveBand[NbrScale3D-1].nz();

	alloc(TabWaveBand);

// Allocation TabBand
	if(TBalloc)
	{
		TabBand = new fltarray [NbrScale3D];
		for (int s = 0; s <  NbrScale3D-1; s++)
			TabBand[s].alloc(CurNc[s],CurNl[s],1);
		TabBand[NbrScale3D-1].alloc(Nx,Ny,Nz);
	}

// Ridgelet transformation on each scale except the coarsest	
#if USE_OMP_RC
//	#pragma omp parallel for 
#endif
	for (int s = 0; s <  NbrScale3D-1; s++)
	{
//cerr<<"s"<<s<<" start..."<<endl;
		if (LocVerbose) cout << " Scale " << s << " Block size calcul = " << TabBlockSize[s] << endl;

		// Pointer to the 1st slice of the Cube, to make a referenced image
		float *Ptr = TabBand[s].buffer();
		Ifloat F;
		F.alloc(Ptr,TabBand[s].ny(),TabBand[s].nx());
		TabRidgelet[s].transform(TabWaveBand[s], F, TabBlockSize[s]);

		if(LocVerbose)
		{
			char filename[64];
			sprintf(filename,"TWB_%d.fits",s);
			writefltarr(filename, TabWaveBand[s]);
		}

//cerr<<"s"<<s<<" end"<<endl;
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

	if(Verbose) cerr << endl << "...End RCurvelet3D::transform_band" <<  endl;
}

/****************************************************************************/

void RCurvelet3D::wavelet_recons(fltarray * & TabWaveBand, fltarray & Cube)
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
		DataW.recons(TabWaveBand,Cube, False);
	}
	else if(TypeW3D==W3D_PMEYER)
	{
		if(Verbose) cerr << " FFT based Meyer WT, size ="<<Cube.nx()<<" "<<Cube.ny()<<" "<<Cube.nz()<<endl;
		POISSON_MWT3D DataW;
		DataW.init(NbrScale3D, Cube.nx(), Cube.ny(), Cube.nz() );
		DataW.normalize_self(TabWaveBand, true);
		DataW.recons(TabWaveBand,Cube, False);
	}
	else
	{
		cerr<< "Undefined wavelet3D transform"<<endl;
		exit(-1);
	}
}

/****************************************************************************/

void RCurvelet3D::recons(fltarray* &TabBand, fltarray* &TabWaveBand, bool TWBalloc) //fltarray &Cube)
{
	bool LocVerbose=false && Verbose;
	if(Verbose) cerr << "RCurvelet3D::recons_band..." <<  endl;

// Allocation TabWaveBand
	if(TWBalloc)
	{
		TabWaveBand = new fltarray [NbrScale3D];
		TabWaveBand[NbrScale3D-1].alloc(NxCube,NyCube,NzCube);
	}

// Ridgelet inverse transformation on each scale except the coarsest	
#if USE_OMP_RC
//	#pragma omp parallel for 
#endif
	for (int s = 0; s <  NbrScale3D-1; s++)
	{
		if (LocVerbose) cout << " Scale " << s << " Block size calcul = " << TabBlockSize[s] << endl;

		// Pointer to the 1st slice of the Cube, to make a referenced image
		float *Ptr = TabBand[s].buffer();
		Ifloat F;
		F.alloc(Ptr,TabBand[s].ny(),TabBand[s].nx());

		if(TWBalloc) TabWaveBand[s].resize(ScaleNx[s], ScaleNy[s], ScaleNz[s]);

		TabRidgelet[s].recons(F, TabWaveBand[s], TabBlockSize[s],TabWaveBand[s].nx(),TabWaveBand[s].ny(),TabWaveBand[s].nz());

		if(LocVerbose)
		{
			char filename[64];
			sprintf(filename,"TWB_recons_%d.fits",s);
			writefltarr(filename, TabWaveBand[s]);
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

	if(Verbose) cerr << endl << "...End RCurvelet3D::recons_band" <<  endl;	
}

/****************************************************************************/

void RCurvelet3D::recons(fltarray* &TabBand, fltarray &Cube, bool Cubealloc)
{
	if(Verbose) cerr << "RCurvelet3D::recons..." <<  endl;

	fltarray* TabWaveBand;

	if(Cubealloc) Cube.alloc(NxCube,NyCube,NzCube);

// Multiscale ridgelet inverse transform
	recons(TabBand, TabWaveBand, true);

// Wavelet inverse transform of the Cube
	wavelet_recons(TabWaveBand, Cube);

	delete [] TabWaveBand;

	if(Verbose) cerr << "...End RCurvelet3D::recons" <<  endl;
}

/****************************************************************************/

void RCurvelet3D::get_band(int s3, int s1, fltarray *TabBand, fltarray &band)
{
//	if(Verbose) cerr << "RCurvelet3D::get_band..." <<  endl;

	if(s3==NbrScale3D-1)
		band=TabBand[s3];
	else
	{
		band.resize(size_nc(s3,s1),size_nl(s3,s1),1);

		for(int x=jpos(s3,s1);x<jpos(s3,s1)+size_nc(s3,s1);x++)
			for(int y=ipos(s3,s1);y<ipos(s3,s1)+size_nl(s3,s1);y++)
				band(x-jpos(s3,s1),y-ipos(s3,s1),0) = (TabBand[s3])(x,y,0);
	}
//	if(Verbose) cerr << "End RCurvelet3D::get_band" <<  endl;
}

/****************************************************************************/

void RCurvelet3D::test_wavelet(fltarray &Cube, fltarray &Recons)
{
	bool LocVerbose=False & Verbose;
	if(Verbose) cerr << "RCurvelet3D::test_wavelet..." <<  endl;

// Wavelet transform of the Cube
	fltarray* TabWaveBand = new fltarray[NbrScale3D];
	wavelet_transform(Cube, TabWaveBand);

// Wavelet inverse transform of the Cube
	wavelet_recons(TabWaveBand, Recons);

	delete [] TabWaveBand;

	if(LocVerbose)
		writefltarr((char*)"TW_Recons.fits", Recons);

	if(Verbose) cerr << endl << "...End RCurvelet3D::test_wavelet" <<  endl;	
}

/****************************************************************************/

void RCurvelet3D::test_ridgelet(fltarray &Cube, fltarray &Recons)
{
	bool LocVerbose=False & Verbose;
	if(Verbose) cerr << "RCurvelet3D::test_wavelet..." <<  endl;

	Ifloat Rid;
	Ridgelet3D TabRidgelet;
//	Note
	TabRidgelet.BlockOverlap = (Bool)(BlockOverlap>0);
	TabRidgelet.set_OverlapFactor(BlockOverlap);
	TabRidgelet.RidTrans = CurTrans;
	TabRidgelet.Verbose = Verbose;
	TabRidgelet.NbrScale = NbrScale1D;

	int CurNl,CurNc;
	Block3D B3D;
//	B3D.BlockOverlap = BlockOverlap;
	B3D.BlockOverlap = (BlockOverlap>0)?True:False;
	B3D.Overlap = BlockOverlap;
	B3D.alloc(Cube.nx(),Cube.ny(),Cube.nz(),BlockSize);
	BlockSize = B3D.block_size();
	CurNl = 3 * BlockSize * BlockSize; 
	CurNc = BlockSize * B3D.nbr_block();
	CurNc *=2;

	Rid.alloc(CurNl, CurNc);

	if(LocVerbose)
		cerr<<"  Rid properties : "<<TabRidgelet.RidTrans<<" bs "<<BlockSize;
	
	TabRidgelet.transform(Cube, Rid, BlockSize);
	io_write_ima_float((char*)"TR_Ridtrans.fits",Rid);
	TabRidgelet.recons(Rid, Recons, BlockSize);

	if(LocVerbose)
		writefltarr((char*)"TR_Recons.fits", Recons);

	if(Verbose) cerr << endl << "...End RCurvelet3D::test_wavelet" <<  endl;	
}

/****************************************************************************/

// The same as extract_stat except that there is no normalisation 
void RCurvelet3D::noise_calibration(fltarray *TabBand, char* Outname)
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

	double var;
	double mean;
	int s3;

// fine scales : 0..N-2
	for(s3=0;s3<NbrScale3D-1;s3++)
	{
		// Allocation
		int NbS1D = get_nbscale1d(s3);
		TabStat[s3] = new double[NbS1D];

		for(int s1=0;s1<NbS1D;s1++)
		{
			mean=0;	// mean
			var=0;	// variance

			for(int i = ipos(s3,s1);		i < ipos(s3,s1) + size_nl(s3,s1); i++)
				for(int j = jpos(s3,s1);	j < jpos(s3,s1) + size_nc(s3,s1); j++)
				{
					mean += (TabBand[s3])(j,i,0);
					var += pow( (TabBand[s3])(j,i,0) , 2 );
				}
			mean /= (size_nl(s3,s1)*size_nc(s3,s1));
			var  /= (size_nl(s3,s1)*size_nc(s3,s1));
			TabStat[s3][s1] = sqrt(var);

		// Write output
			if(LocVerbose)	cerr<<" ("<<ipos(s3,s1)<<":"<<ipos(s3,s1) + size_nl(s3,s1)-1<<")x("<<
					jpos(s3,s1)<<":"<<jpos(s3,s1) + size_nc(s3,s1)-1<<") ";
			if(Verbose)	cerr<< "\tStat echelle ("<<s3<<","<<s1<<") : (mean,std) = ("<<
					mean<<",\t"<<TabStat[s3][s1]<<")"<<endl;
			cstat <<mean<<"\t"<<TabStat[s3][s1]<<endl;
		}
	}

// coarse scale : N-1
	// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats

	for(int s3=0;s3<NbrScale3D-1;s3++)
		delete [] TabStat[s3];
	delete [] TabStat;

	cstat.close();
	if(Verbose) cerr<<"...End Noise_calibration"<<endl;
}

/****************************************************************************/


/*void RCurvelet3D::extract_stat(fltarray *TabBand, char* Outname, bool UseCubeSigma)
{
	bool LocVerbose = True & Verbose;
	if(Verbose) cerr<<"Extract_stat..."<<endl;

	double ** TabStat = new double*[NbrScale3D];
	double ** MaxCoef = new double*[NbrScale3D-1];//we ommit the coarse scale
	int Maxi=0,Maxj=0;

// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_stat.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_maxcoef.dat");
	fstream ccoef;
	ccoef.open (MaxCoefname, fstream::out);

	double var;
	double mean;
	int s3;
	float val;
	
// fine scales : 0..N-2
	for(s3=0;s3<NbrScale3D-1;s3++)
	{
		// Allocation
		int NbS1D = get_nbscale1d(s3);
		TabStat[s3] = new double[NbS1D];
		MaxCoef[s3] = new double[NbS1D];

		for(int s1=0;s1<NbS1D;s1++)
		{
			mean=0;	// mean
			var=0;	// variance
			MaxCoef[s3][s1] = 0;

			for(int i = ipos(s3,s1);		i < ipos(s3,s1) + size_nl(s3,s1); i++)
				for(int j = jpos(s3,s1);	j < jpos(s3,s1) + size_nc(s3,s1); j++)
				{
				// Normalisatiton
					if(UseCubeSigma) val = (TabBand[s3])(j,i,0) / (CubeSigma[s3][s1])(i,(j-jpos(s3,s1))%size_nc(s3,s1,0));
					else val = (TabBand[s3])(j,i,0) / TabSigma(s3,s1);

		for(int s1=0;s1<NbS1D;s1++)
		{
			float cnt=0,cnt2=0;
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
					{
						cnt+=1.;
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						if(UseCubeSigma) lvl = SigmaNoise * Nsig * (CubeSigma[s3][s1])(i,j);
						else lvl = SigmaNoise * Nsig * TabSigma(s3,s1);





					
					mean += val;
					var  += pow( val , 2 );
					
					if( abs(val) > MaxCoef[s3][s1] )
					{
						MaxCoef[s3][s1] = abs(val);
						Maxi=i;
						Maxj=j;
					}
				}
			
			mean /= (size_nl(s3,s1)*size_nc(s3,s1));
			var /= (size_nl(s3,s1)*size_nc(s3,s1));
			TabStat[s3][s1] = sqrt(var);

			if(LocVerbose) cerr<<" ("<<ipos(s3,s1)<<":"<<ipos(s3,s1) + size_nl(s3,s1)-1<<")x("<<jpos(s3,s1)<<":"<<jpos(s3,s1) + size_nc(s3,s1)-1<<") ";
			if(Verbose) cerr<< "\tStat echelle ("<<s3<<","<<s1<<") : (mean,std) = ("<<mean<<", \t"<<TabStat[s3][s1]<<")"<<endl;
			if(Verbose)	cout << "MaxCoef Value = "<<MaxCoef[s3][s1] << endl;
			
			cstat <<mean<<"\t"<<TabStat[s3][s1]<<endl;
			ccoef << s3 << "\t" << s1 << "\t" << MaxCoef[s3][s1] << "\t" << Maxi << "\t" << Maxj << endl;
		}}
	}

// coarse scale : N-1
// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats
	s3=NbrScale3D-1;
	TabStat[s3] = new double[1];
	{
		mean=0;	// mean
		var=0;	// variance

		for(int i=0; i<TabBand[s3].nx(); i++)
			for(int j=0; j<TabBand[s3].ny(); j++)
				for(int k=0; k<TabBand[s3].nz(); k++)
				{
					mean += (TabBand[s3])(i,j,k);
					var += pow( (TabBand[s3])(i,j,k) , 2 );
				}
		mean /= (TabBand[s3].nx()*TabBand[s3].ny()*TabBand[s3].nz());
		var /= (TabBand[s3].nx()*TabBand[s3].ny()*TabBand[s3].nz());
		TabStat[s3][0] = sqrt(var);

		if(LocVerbose) cerr<<" ("<<TabBand[s3].nx()<<"x"<<TabBand[s3].ny()<<"x"<<TabBand[s3].nz()<<")";
		if(LocVerbose) cerr<< " Stat echelle ("<<s3<<") : (mean,std) = ("<<mean<<",\t"<<TabStat[s3][0]<<")"<<endl;
		
		cstat << mean<<"\t"<<TabStat[s3][0]<<endl;
	}

	cstat.close();
	ccoef.close();

	for(int s3=0;s3<NbrScale3D;s3++)
		delete [] TabStat[s3];
	delete [] TabStat;
	for(int s3=0;s3<NbrScale3D-1;s3++)
		delete [] MaxCoef[s3];
	delete [] MaxCoef;

	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}*/
void RCurvelet3D::extract_stat(fltarray *TabBand, char* Outname, bool UseCubeSigma, bool normalize)
{
	if(Verbose) cerr<<"Extract_stat..."<<endl;

	double ** MaxCoef = new double*[NbrScale3D-1];//we ommit the coarse scale

// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_stat.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	char MaxCoefname[250];
	strcpy(MaxCoefname, Outname);
	strcat(MaxCoefname, "_maxcoef.dat");
	fstream ccoef;
	ccoef.open (MaxCoefname, fstream::out);

	
	TabStat.alloc(NbrScale3D,get_nbscale1d(0),6);
// fine scales : 0..N-2
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		// Allocation
		int NbS1D = get_nbscale1d(s3);
		MaxCoef[s3] = new double[NbS1D];

		for(int s=0;s<NbS1D;s++)
		{
			double m1=0,m2=0,m3=0,m4=0;
			float val,norm=1;
			float Minni = (TabBand[s3])(jpos(s3,s),ipos(s3,s),0);
			float Maxxi = Minni;
			float Maxi=ipos(s3,s),Maxj=jpos(s3,s);
			MaxCoef[s3][s] = (TabBand[s3])(jpos(s3,s),ipos(s3,s),0);
//cerr<<"cubesig size="<<(CubeSigma[s3][s]).nl()<<","<<(CubeSigma[s3][s]).nc()<<endl;
//cerr<<"ijpos size="<<ipos(s3,s)<<"+"<<size_nl(s3,s)<<","<<jpos(s3,s)<<"+"<<size_nc(s3,s)<<endl;
			for(int i = 	ipos(s3,s);	i < ipos(s3,s) + size_nl(s3,s); i++)
				for(int j = jpos(s3,s);	j < jpos(s3,s) + size_nc(s3,s); j++)
				{
				// Normalisatiton
					if(normalize)
					{
						if(UseCubeSigma) norm = (CubeSigma[s3][s])( i-ipos(s3,s) , (j-jpos(s3,s)) % size_nc(s3,s,0) );
						else  norm = TabSigma(s3,s);
					}
					val = (TabBand[s3])(j,i,0) / norm;
					
				// Moments
					m1 += val;
					m2 += pow( val , 2 );
					m3 += pow( val , 3 );
					m4 += pow( val , 4 );
					
					if(val<Minni) Minni=val;
					if(val>Maxxi) Maxxi=val;
				}
			MaxCoef[s3][s] = max(abs(Maxxi),abs(Minni));
			
			int N = size_nl(s3,s)*size_nc(s3,s);
			
			TabStat(s3,s,0) = m1/N;
			moment4_center(N, TabStat(s3,s,0), m2/N, m3/N, m4/N, TabStat(s3,s,1), TabStat(s3,s,2), TabStat(s3,s,3));
			TabStat(s3,s,4) = Minni;
			TabStat(s3,s,5) = Maxxi;
			
			if(Verbose) cerr<<" ("<<ipos(s3,s)<<":"<<ipos(s3,s) + size_nl(s3,s)-1<<")x("<<jpos(s3,s)<<":"<<jpos(s3,s) + size_nc(s3,s)-1<<") ";
			if(Verbose)
				cerr << s3 <<"\t"<< s <<"("<<TabSigma(s3,s)<<")\t"<< TabStat(s3,s,0) <<"\t"<< TabStat(s3,s,1) <<"\t"<< TabStat(s3,s,2) <<"\t"<<
						 TabStat(s3,s,3) <<"\t"<< TabStat(s3,s,4) << "\t" << TabStat(s3,s,5) << endl;
			if(Outname!=NULL)
			{
				cstat << s3 <<"\t"<< s <<"\t"<< TabStat(s3,s,0) <<"\t"<< TabStat(s3,s,1) <<"\t"<< TabStat(s3,s,2) <<"\t"<<
						 TabStat(s3,s,3) <<"\t"<< TabStat(s3,s,4) << "\t" << TabStat(s3,s,5) << endl;
				ccoef << s3 << "\t" << s << "\t" << MaxCoef[s3][s] << "\t" << Maxi << "\t" << Maxj << endl;
			}
		}
	}

// coarse scale : N-1
// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats
	{
		int s3=NbrScale3D-1;
		int s=0;
		float Mini,Maxi;
		moment4(TabBand[s3].buffer(), TabBand[s3].n_elem(),TabStat(s3,s,0), TabStat(s3,s,1), TabStat(s3,s,2), TabStat(s3,s,3), 
				Mini,Maxi);
		TabStat(s3,s,4) = Mini;
		TabStat(s3,s,5) = Maxi;
		
		if(Verbose)
			cerr << s3 <<"\t"<< s <<"\t"<< TabStat(s3,s,0) <<"\t"<< TabStat(s3,s,1) <<"\t"<< TabStat(s3,s,2) <<"\t"<<
					 TabStat(s3,s,3) <<"\t"<< TabStat(s3,s,4) << "\t" << TabStat(s3,s,5) << endl;
		if(Outname!=NULL)
			cstat << s3 <<"\t"<< s <<"\t"<< TabStat(s3,s,0) <<"\t"<< TabStat(s3,s,1) <<"\t"<< TabStat(s3,s,2) <<"\t"<<
					 TabStat(s3,s,3) <<"\t"<< TabStat(s3,s,4) << "\t" << TabStat(s3,s,5) << endl;
	}

// Write the stat file
	if(Outname!=NULL) 
	{
		char filename[64];
		sprintf(filename,"%s_stat.fits",Outname);
		writefltarr(filename, TabStat);
	}
		
	cstat.close();
	ccoef.close();

	for(int s3=0;s3<NbrScale3D-1;s3++)
		delete [] MaxCoef[s3];
	delete [] MaxCoef;

	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}

/****************************************************************************/

void RCurvelet3D::normalize(fltarray *TabBand, fltarray *TabBandNorm, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::normalize..."<<endl;
	//bool LocVerbose = True & Verbose;
	float val;

	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
				for(int j = jpos(s3,s1); j < jpos(s3,s1) + size_nc(s3,s1); j++)
				{
					if(UseCubeSigma) val = (TabBand[s3])(j,i,0) / (CubeSigma[s3][s1])( i-ipos(s3,s1) , (j-jpos(s3,s1)) % (size_nl(s3,s1)/nbr_block(s3)) );
					else val = (TabBand[s3])(j,i,0) / TabSigma(s3,s1);
					(TabBandNorm[s3])(j,i,0) = val;
				}

	// Coarse scale
	// Do not normalize the coarsest scale
//	int s3=NbrScale3D-1;
//	for(int i=0; i<TabBand[s3].nx(); i++)
//		for(int j=0; j<TabBand[s3].ny(); j++)
//			for(int k=0; k<TabBand[s3].nz(); k++)
//				(TabBandNorm[s3])(i,j,k)=(TabBand[s3])(i,j,k);

	if(Verbose) cerr<<"...End RCurvelet3D::normalize"<<endl;
}

/****************************************************************************/

void RCurvelet3D::normalize_self(fltarray *TabBand,bool inverse, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::normalize_self..."<<endl;
	//bool LocVerbose = True & Verbose;
	
	float val;
	
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
				for(int j = jpos(s3,s1); j < jpos(s3,s1) + size_nc(s3,s1); j++)
					if(inverse)
					{
						if(UseCubeSigma) val = (TabBand[s3])(j,i,0) * (CubeSigma[s3][s1])( i-ipos(s3,s1) , (j-jpos(s3,s1)) % (size_nl(s3,s1)/nbr_block(s3)) );
						else val = (TabBand[s3])(j,i,0) * TabSigma(s3,s1);
						(TabBand[s3])(j,i,0) = val;
					}
					else
					{
						if(UseCubeSigma) val = (TabBand[s3])(j,i,0) / (CubeSigma[s3][s1])( i-ipos(s3,s1) , (j-jpos(s3,s1)) % (size_nl(s3,s1)/nbr_block(s3)) );
						else val = (TabBand[s3])(j,i,0) / TabSigma(s3,s1);
						(TabBand[s3])(j,i,0) = val;
					}

	// Coarse scale
	// Do not normalize the coarsest scale

	if(Verbose) cerr<<"...End RCurvelet3D::normalize_self"<<endl;
}

/*********************************************************************/

void RCurvelet3D::threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<","<<UseCubeSigma<<")"<<endl;

	float lvl;
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		// Allocation
		int NbS1D = get_nbscale1d(s3);
		for(int s1=0;s1<NbS1D;s1++)
		{
			float cnt=0,cnt2=0;
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
					{
						cnt+=1.;
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						if(UseCubeSigma) lvl = SigmaNoise * Nsig * (CubeSigma[s3][s1])(i,j);
						else lvl = SigmaNoise * Nsig * TabSigma(s3,s1);
						
						if( abs((TabBand[s3])(j+x0,i,0)) < lvl )
						{
							cnt2+=1.;
							(TabBand[s3])(j+x0,i,0)=0; // hard
						}
						else if(FilterType==FT_SOFT) (TabBand[s3])(j+x0,i,0) -= (2*int( (TabBand[s3])(j+x0,i,0)>0 )-1)*lvl;
						//(TabBand[s3])(j,i,0)+=4*s3+s1;
					}
			}
			if(Verbose) cerr<<"Band "<<s3<<","<<s1<<") n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;
		}
	}

	// Coarse scale
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);

	if(Verbose) cerr<<"...End RCurvelet3D::threshold"<<endl;
}

/****************************************************************************/

int RCurvelet3D::get_RC_pos(int kv,int kh,int by,int bx,int LocalBS,int N)
{
	int Ph = ( bx + kh*LocalBS ) % (2*N);
	int Pv = ( by + kv*LocalBS ) * (2*N);
	int h=-1, v=-1;
	//*
	if(Pv >= N*2*N)
	{
		int pv = Pv/(2*N)-(N-1);

		if(Ph <= N/2)
		{
			v = 2*N*N + (N/2 + Ph)*N ;
			h = N - pv ;
		}
		else if(Ph > N+N/2)
		{
			v = 2*N*N + (N/2 + 2*N - Ph)*N ;
			h = pv ;
		}
		else // N/2 < Ph <= N+N/2
		{
			v = 2*N*N + (N - pv)*N ;
			h = 3*N/2 - Ph ;
		}
//	int y = v + h ;
	//cerr<<" kv kh by bx -> h v y = "<<kv<<","<<kh<<","<<by<<","<<bx<<" -> "<<h<<" "<<v<<" "<<y<<"   "<<endl;
	}
	else 
	{
		v = Pv ;
		h = Ph ;
//	int y = v + h ;
	//cerr<<" kv kh by bx -> h v y = "<<kv<<","<<kh<<","<<by<<","<<bx<<" -> "<<h<<" "<<v<<" "<<y<<"   "<<endl;
	}
	return v + h ;
}


void RCurvelet3D::wiener(fltarray *TabBand, float noise_lvl, int LocalBS, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::wiener("<<noise_lvl<<","<<LocalBS<<","<<UseCubeSigma<<")..."<<endl;
//	bool LocVerbose = true & Verbose;
	float val;
	float noise2 = noise_lvl*noise_lvl;
	
	// Fine wavelet scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		// Allocation
		int NbS1D = get_nbscale1d(s3);

		int N = TabBlockSize[s3]; // for height : number of angles 

		// Ridgelet scales
		for(int s1=0;s1<NbS1D;s1++)
		{
			int Nw = (TabBlockSize[s3]-1)/(1<<s1) + 1; // for width : 3rd coordinate of wiener block // 33 17 9 5 3
			int Nv = ceil(float(N)/float(LocalBS));
			int Nh = ceil(float(2*N)/float(LocalBS));
			int Np = ceil(float(Nw)/float(LocalBS));
			//fltarray coef_wiener(Nv, Nh, Np);

			// Spatial blocks B
			#if USE_OMP_RC
				#pragma omp parallel for  
			#endif
			for(int B=0;B<nbr_block(s3);B++)
			//int B=0;
			{
				fltarray coef_wiener(Nw,3*N*N,1); // 3*N² angles and Nw positions
				
				coef_wiener.init(-2);

//				cerr<<"  N "<<N<<" kh<"<<ceil(float(2*N)/float(LocalBS))<<", kv<"<<ceil(float(N)/float(LocalBS))<<" = "<<
//					ceil(float(2*N)/float(LocalBS))*ceil(float(N)/float(LocalBS))<<
//					" LocalBS="<<LocalBS<<"x"<<LocalBS<<" => total="<<ceil(float(2*N)/float(LocalBS))*ceil(float(N)/float(LocalBS))*LocalBS*LocalBS<<endl;

				// Wiener blocks (int angle2 and space) for mainly vertical Ridgelets
				for(int kv=0 ; kv < Nv ; kv++)
					for(int kh=0 ; kh < Nh ; kh++)
						for(int kp = 0; kp < Np; kp++)
						{
							double sigma2 = 0.0;
							float cnt=0;

						// Sigma calculation
							// Pixels in a wiener block (angle)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int y = get_RC_pos(kv,kh,by,bx,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bp = 0; bp < LocalBS; bp++)
									{
										cnt+=1;
										// spatial position = position in the block + block_position
										int x = (bp + kp*LocalBS) % Nw;
										int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);

										//cerr<<" x,x0 y = "<<x<<","<<x0<<" "<<y<<endl;
										//(TabBand[s3])(x+x0,y,0) = B;
										//(TabBand[s3])(x+x0,y,0) *= sqrt(sig / (sig+TabSigma(s3,s1)*noise_lvl));
										
										if(UseCubeSigma) val = (TabBand[s3])(x+x0,y,0) / (CubeSigma[s3][s1])(x,y);
										else val = (TabBand[s3])(x+x0,y,0) / TabSigma(s3,s1);
										
										sigma2+=pow(val,2);
									}
								}

							float sig2 = max( 0.0, sigma2/cnt - noise2 );
							float norm = sig2 / (sig2+noise2);
							//if(LocVerbose) cerr<< " Vert norm("<<kh<<" "<<kv<<" "<<kp<<") = "<<norm;

						// Store the coef in the table
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int y = get_RC_pos(kv,kh,by,bx,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bp = 0; bp < LocalBS; bp++)
									{
										int x = (bp + kp*LocalBS) % Nw;
										if( coef_wiener(x,y,0) < -1 )
											coef_wiener(x,y,0) = norm;
									}
								}
						}// end wiener blocks




				// Wiener blocks (int angle2 and space) for mainly horizontal Ridgelets
				for(int kv=0 ; kv < Nv ; kv++)
					for(int kh=0 ; kh < Nv ; kh++)
						for(int kp = 0; kp < Np; kp++)
						{
							double sigma2 = 0.0;
							float cnt=0;

						// Sigma calculation
							// Pixels in a wiener block (angle)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int Ph = ( bx + kh*LocalBS );// NOTE %Nv no ?
									int Pv = 2*N*N + ( by + kv*LocalBS ) * N;
									int pv = by + kv*LocalBS;

									int y = Pv+Ph;
									//int y = get_RC_pos(kv,kh,by,bx,LocalBS,N);
									if(pv<N && Ph<N)
									{
										// pixels in a wiener block (space)
										for(int bp = 0; bp < LocalBS; bp++)
										{
											cnt+=1;
											// spatial position = position in the block + block_position
											int x = (bp + kp*LocalBS) % Nw;
											int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
											
											if(UseCubeSigma) val = (TabBand[s3])(x+x0,y,0) / (CubeSigma[s3][s1])(x,y);
											else val = (TabBand[s3])(x+x0,y,0) / TabSigma(s3,s1);
										
											sigma2+=pow(val,2);
										}
									}
								}

							float sig2 = max( 0.0, sigma2/cnt - pow(noise_lvl,2) );
							float norm = sig2 / (sig2+noise2);
							//if(LocVerbose) cerr<< " Horiz norm("<<kh<<" "<<kv<<" "<<kp<<") = "<<norm;

						// Store the coef in the table
							// Pixels in a wiener block (angle)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int Ph = ( bx + kh*LocalBS );
									int Pv = 2*N*N + ( by + kv*LocalBS ) * N;
									int pv = by + kv*LocalBS;
									int y = Pv+Ph;

									if(pv<N && Ph<N)
									{
										// pixels in a wiener block (space)
										for(int bp = 0; bp < LocalBS; bp++)
										{
											// spatial position = position in the block + block_position
											int x = (bp + kp*LocalBS) % Nw;
											if( coef_wiener(x,y,0) < -1 )
												coef_wiener(x,y,0) = norm;
										}
									}
								}
						}// end wiener blocks


			// Apply the wiener coefficient
				// Wiener blocks (int angle2 and space) for mainly vertical Ridgelets
				for(int kv=0 ; kv < Nv ; kv++)
					for(int kh=0 ; kh < Nh ; kh++)
						for(int kp = 0; kp < Np; kp++)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int y = get_RC_pos(kv,kh,by,bx,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bp = 0; bp < LocalBS; bp++)
									{
										// spatial position = position in the block + block_position
										int x = (bp + kp*LocalBS) % Nw;
										int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
										if( coef_wiener(x,y,0) > -1 )
										{
											(TabBand[s3])(x+x0,y,0) *= coef_wiener(x,y,0);
											coef_wiener(x,y,0) = -2;
										}

									}
								}

				for(int kv=0 ; kv < Nv ; kv++)
					for(int kh=0 ; kh < Nv ; kh++)
						for(int kp = 0; kp < Np; kp++)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int Ph = ( bx + kh*LocalBS );
									int Pv = 2*N*N + ( by + kv*LocalBS ) * N;
									int pv = by + kv*LocalBS;
									int y = Pv+Ph;
									if(pv<N && Ph<N)
									{
										// pixels in a wiener block (space)
										for(int bp = 0; bp < LocalBS; bp++)
										{
											// spatial position = position in the block + block_position
											int x = (bp + kp*LocalBS) % Nw;
											int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
											if( coef_wiener(x,y,0) > -1 )
											{
												(TabBand[s3])(x+x0,y,0) *= coef_wiener(x,y,0);
												coef_wiener(x,y,0) = -2;
											}
										}
									}
								}

			}// end physical block
		}// end 1D scale
	}// end 3D scale

	// Coarse scale
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);

	if(Verbose) cerr<<"...End RCurvelet3D::wiener"<<endl;
}

/****************************************************************************/

void RCurvelet3D::values_at(fltarray *TabBand, char * filename, char* Outname)
{
	int S1,S3,i,j;
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
			int NbS1D = get_nbscale1d(s3);
			for(int s1=0;s1<NbS1D;s1++)
			{
				if(!fscanf(fic,"%d %d %f %d %d",&S3,&S1,&m,&i,&j)) cerr<<"Error while reading coordinates in "<<filename<<endl;
				cval<<S3<<"\t"<<S1<<"\t"<<abs((TabBand[s3])(j,i,0))<<"\t"<<i<<"\t"<<j<<endl;
			}
		}
		fclose(fic);
	}
	else cerr<<"Warnig: File "<<filename<<" not found, for use in BCurvelet3D::values_at"<<endl;
	cval.close();
}

/****************************************************************************/

void RCurvelet3D::update_significant(fltarray *TabBand, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::update_significant(.,.,"<<SigmaNoise<<","<<NSigma<<","<<force4<<","<<UseCubeSigma<<")"<<endl;
	
	float lvl;
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
					{
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						//if(UseCubeSigma) lvl = threshold * Nsig * (CubeSigma[s3][s1])(i,j);
				// ?????//else 
							lvl = SigmaNoise * Nsig * TabSigma(s3,s1);
						
						if( abs((TabBandO[s3])(j+x0,i,0)) > lvl )
							if( abs((TabBand[s3])(j+x0,i,0)-(TabBandO[s3])(j+x0,i,0)) > SigmaNoise * TabSigma(s3,s1)/2)
								(TabBand[s3])(j+x0,i,0)=(TabBandO[s3])(j+x0,i,0);
					}
			}
	if(Verbose) cerr<<"...End RCurvelet3D::update_significant"<<endl;
}

/****************************************************************************/

void RCurvelet3D::select_significant(fltarray *TabBand, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::select_significant(.,.,"<<SigmaNoise<<","<<NSigma<<","<<force4<<","<<UseCubeSigma<<")"<<endl;
	
	float lvl;
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
					{
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						lvl = SigmaNoise * Nsig * TabSigma(s3,s1);
						
						if( abs((TabBandO[s3])(j+x0,i,0)) < lvl )
							(TabBand[s3])(j+x0,i,0)=0;
					}
			}
	if(Verbose) cerr<<"...End RCurvelet3D::select_significant"<<endl;
}

/****************************************************************************/

void RCurvelet3D::select_significant(fltarray *TabBand, intarray *support, fltarray *TabBandO)
{
	if(Verbose) cerr<<"RCurvelet3D::select_significant(.,.,.)"<<endl;
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
						(TabBand[s3])(j+x0,i,0) = (support[s3])(j+x0,i,0) ? (TabBandO[s3])(j+x0,i,0) : 0 ;
			}
	if(Verbose) cerr<<"...End RCurvelet3D::select_significant"<<endl;
}

/****************************************************************************/

void RCurvelet3D::support_significant(intarray *support, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4, bool UseCubeSigma)
{
	if(Verbose) cerr<<"RCurvelet3D::support_significant(.,.,"<<SigmaNoise<<","<<NSigma<<","<<force4<<","<<UseCubeSigma<<")"<<endl;
	
	float lvl;
	// Fine scales
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
			for(int B=0;B<nbr_block(s3);B++)
			{
				int x0 = jpos(s3,s1) + B*size_nc(s3,s1,B);
				for(int i = ipos(s3,s1); i < ipos(s3,s1) + size_nl(s3,s1); i++)
					for(int j = 0; j < size_nc(s3,s1,B); j++)
					{
						float Nsig=NSigma;
						if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
						lvl = SigmaNoise * Nsig * TabSigma(s3,s1);
						
						if( abs((TabBandO[s3])(j+x0,i,0)) > lvl )
							(support[s3])(j+x0,i,0)=1;
					}
			}
	if(Verbose) cerr<<"...End RCurvelet3D::support_significant"<<endl;
}

/****************************************************************************/

void RCurvelet3D::fdr(fltarray * &TabBand, float Alpha, float SigmaNoise)
{
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
	bool LocVerbose = true & Verbose;
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		for(int s1=0;s1<get_nbscale1d(s3);s1++)
		{
			
		//getband
			fltarray Band;
			get_band(s3, s1, TabBand, Band);

		// Cumulative of the band
			CArrProb *TabCP;
			TabCP = new CArrProb;
			TabCP->set(Band,false);
			
			float lvl = TabSigma(s3,s1)*SigmaNoise;
			
		// P-Values
			fltarray PVal(Band.nx(),Band.ny(),Band.nz());
			for (int i=0; i < Band.nx(); i++)
				for (int j=0; j < Band.ny(); j++)
					for (int k=0; k < Band.nz(); k++)
						PVal(i,j,k) = prob_noise(Band(i,j,k)/lvl);
	
			double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);
			
			if(LocVerbose) cerr<<" band ("<<s3<<","<<s1<<") : PDet="<<PDet<<endl;

		// threshold
			for (int i=0; i < Band.ny(); i++)
			{
				int j0 = jpos(s3,s1);
				for (int j=0; j < Band.nx(); j++)
				{
					if( PVal(j,i,0) > PDet) (TabBand[s3])(j+j0,i,0) = 0;
				}
			}
		}
	
	if(no_coarse)
		TabBand[NbrScale3D-1].init(0.0);

	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

/****************************************************************************/

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

void RCurvelet3D::write (char *Name, fltarray * & TabBand, bool Normalize)
{
	if(Verbose) cerr<<"RCurvelet3D::write("<<Name<<",.,"<<Normalize<<")"<<endl;

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
	if ( ffpkyj(fptr, (char*)"NScale3D", (long) NbrScale3D, (char*)"Number of 3D scales", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"NScale1D", (long) NbrScale1D, (char*)"Number of 1D scales", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of 3d blocks", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Overlaping blocks flag", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"TypeW3D", (long) TypeW3D, (char*)"Type of wavelets used", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"CurTrans", (long) CurTrans, (char*)"Type of RCurvelet Transform(Ridgelet type)", &status))
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
		if ( ffpkyj(fptr, (char*)"Nx", (long) naxes[0], (char*)"x size of the Band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny", (long) naxes[1], (char*)"y size of the Band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz", (long) naxes[2], (char*)"z size of the Band", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) ScaleNx[s], (char*)"x size of the scale", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) ScaleNy[s], (char*)"y size of the scale", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) ScaleNz[s], (char*)"z size of the scale", &status))
			PrintError( status );  

	// save the data
		if ( ffppre(fptr, group, 1, TabBand[s].nx()*TabBand[s].ny()*TabBand[s].nz(), (TabBand[s]).buffer(), &status) )
			PrintError( status );
	}

// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  

	if(Verbose) cerr<<"...end RCurvelet3D::write_multi"<<endl;
}

/*********************************************************************/

void RCurvelet3D::read(char *Name, fltarray * & TabBand, bool *Normalize)
{
	if(Verbose) cerr<<"RCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
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
	int _NbrScale3D,_NbrScale1D,_BlockSize,_Nx,_Ny,_Nz;
	float _BlockOverlap;
	type_ridgelet3d_WTtrans _CurTrans;
	type_wavelet3D _TypeW3D;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*Normalize = (bool)mon_long;
	if (ffgkyj(fptr,(char*)"NScale3D", &mon_long, comment, &status)) PrintError( status );
	_NbrScale3D = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NScale1D", &mon_long, comment, &status)) PrintError( status );
	_NbrScale1D = (int)mon_long;
	if(nhead!=_NbrScale3D+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale3D="<<_NbrScale3D<<endl; exit(0); }
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	_BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	_BlockOverlap = (float)mon_long;
	if (ffgkyj(fptr,(char*)"CurTrans", &mon_long, comment, &status)) PrintError( status );
	_CurTrans = (type_ridgelet3d_WTtrans)mon_long;
	if (ffgkyj(fptr,(char*)"TypeW3D", &mon_long, comment, &status)) PrintError( status );
	_TypeW3D = (type_wavelet3D)mon_long;
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;

//	init the structure 
	CurTrans=_CurTrans;
	BlockOverlap=_BlockOverlap;
	init(_Nx, _Ny, _Nz, _NbrScale3D, _NbrScale1D, _BlockSize, _TypeW3D);

	tab_block_size_init();

	TabRidgelet = new Ridgelet3D[NbrScale3D-1];
	for (int s = 0; s <  NbrScale3D-1; s++) 
		ridgelet_init(TabRidgelet[s],TabBlockSize[s]);
	AllocClass = True;
	CurNl = new int [NbrScale3D-1];
	CurNc = new int [NbrScale3D-1];
	ScaleNx = new int [NbrScale3D];
	ScaleNy = new int [NbrScale3D];
	ScaleNz = new int [NbrScale3D];



	TabBand = new fltarray[NbrScale3D];

// read data
	for(int s=0;s<_NbrScale3D;s++)
	{
//cerr<<"save : "<<s<<" "<<TabBeamlet[s].get_BlockSize()<<" "<< TabBeamlet[s].get_BlockSize()<<" "<< TabBeamlet[s].get_LinNbPlane() <<endl;
		int NX,NY,NZ;
		if (fits_movabs_hdu(fptr, s+2, NULL, &status)) PrintError( status );

		if (ffgkyj(fptr,(char*)"Nx", &mon_long, comment, &status)) PrintError( status );
		NX = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Ny", &mon_long, comment, &status)) PrintError( status );
		NY = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nz", &mon_long, comment, &status)) PrintError( status );
		NZ = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
		ScaleNx[s] = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
		ScaleNy[s] = (int)mon_long;
		if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
		ScaleNz[s] = (int)mon_long;


//cerr<<"load : "<<s<<" "<<NX<<" "<<NY<<" "<<NZ<<endl;
		TabBand[s].alloc(NX,NY,NZ);

		if(s<_NbrScale3D-1)
		{
			CurNl[s] = NY; 
			CurNc[s] = NX;
			TabRidgelet[s].set_block_param(ScaleNx[s],ScaleNy[s],ScaleNz[s],TabBlockSize[s]);
			//cerr<<"scale, tabsize, nbr_blk "<<s<<","<<TabRidgelet[s].rid_size(s)<<","<<TabRidgelet[s].rid_block_nbr()<<endl;
		}
		if (ffgpve(fptr, 1, 1, NX*NY*NZ, nulval, (TabBand[s]).buffer(), &anynul, &status)) PrintError( status );
	}

// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );

	if(Verbose) cerr<<"...end RCurvelet3D::read_multi"<<endl;
}

/****************************************************************************/

void RCurvelet3D::temp(fltarray *TabBand)
{

/*	
	for(int s3=0;s3<NbrScale3D-1;s3++)
	{
		cerr<<"s,nx,ny "<<s3<<","<<TabBand[s3].nx()<<","<<TabBand[s3].ny()<<endl;
		for(int i=0; i<TabBand[s3].nx(); i++)
			for(int j=0; j<TabBand[s3].ny(); j++)
				(TabBand[s3])(i,j,0)=0;
	}
	
	{
		int s3=NbrScale3D-1;
		for(int i=0; i<TabBand[s3].nx(); i++)
			for(int j=0; j<TabBand[s3].ny(); j++)
				for(int k=0; k<TabBand[s3].nz(); k++)
					(TabBand[s3])(i,j,k)=0;
	}

	// center object 
	//(TabBand[1])(1450,0,0)=1;
	
	
	//(TabBand[1])(1450,1,0)=1;
	
	for(int s3=0;s3<NbrScale3D-1;s3++)
		cerr<<get_nbscale1d(s3)<<",";
	
	int s3= 0;
	cerr<<"nblock="<<nbr_block_nx(s3)<<","<<nbr_block_ny(s3)<<","<<nbr_block_nz(s3)<<endl;
	int B = num_block(s3,7,7,7);
	int s = 0;
	float p=0.4, q=0.3;
//	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
//	s3=0; s=3; B=num_block(s3,7,7,7); p=0.3; q=0.6;
//	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
	s3=1; s=3; B=num_block(s3,3,3,3); p=0.3; q=0.6;
//	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
//	s3=2; s=2; B=num_block(s3,3,3,3); p=0.3; q=0.6;
//	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
//	s3=2; s=0; B=num_block(s3,3,3,3); p=0.3; q=0.6;
//	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
//	s3=2; s=1; B=num_block(s3,7,7,7); p=0.3; q=0.6;
	(TabBand[s3])(jpos(s3,s,B)+int(size_nc(s3,s,B)*p), ipos(s3,s)+int(size_nl(s3,s)*q), 0) = 1;
	
	//for(int s3=0;s3<NbrScale3D-1;s3++)
	//	for(int s
	//		cerr<<"("<<ipos(s3,s)<<":"<<ipos(s3,s)+size_nl(s3,s)<<")x("<<jpos(s3,s,B)<<":"<<jpos(s3,s,B)+size_nc(s3,s,B)<<")"<<endl;
	
*/	
}










