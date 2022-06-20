
#include "bfct.h"

extern bool Verbose;
extern bool ChangeTabDir;

static const double PI4 = M_PI/4.0;
static const double PI2 = M_PI/2.0;

// ********************************************************************
// ********************************************************************
//				PARAMETERS CLASS
// ********************************************************************
// ********************************************************************

BFCurvelet3D_params::BFCurvelet3D_params(fltarray a)
{
	reset_params();
	Nx=a.nx();
	Ny=a.ny();
	Nz=a.nz();
}

void BFCurvelet3D_params::reset_params()
{
// Transform properties
	Nx				=0;	
	Ny				=0;	
	Nz				=0;
	NbrScale		=3;
	NbrDir2d		=16;			
	BlockSize		=0;
	lapped			=false;
	BlockOverlap	=False;
	extract_type	=DEF_TYPE_EXTRACTION;
	RealData		=True;
	
// Memory properties
	lowmem			=false;
	reuse_input_data=false;
	
// Filtering properties
	no_norm			=false;
	no_fine			=false;
	no_coarse		=false;
	threshold_coarse=false;
	Positivity		=false;
	use_min			=false;
	use_max			=false;
	min_value		=0;
	max_value		=255;

	FilterType 		=FT_HARD;
	SigmaNoise		=-1;
	NSigma			=3;
	Alpha			=0.01;
	WienerBS		=3;
	Rho				=0.5;
	Kmin			=3.5;
	Kmean			=7;
	Lmax			=0.5;
	force4sigma		=false;
}

void BFCurvelet3D_params::import_params(BFCurvelet3D_params Q)
{
// Transform properties
	Nx				= Q.Nx;
	Ny				= Q.Ny;
	Nz				= Q.Nz;
	NbrScale		= Q.NbrScale;
	NbrDir2d		= Q.NbrDir2d;
	BlockSize		= Q.BlockSize;
	lapped			= Q.lapped;
	BlockOverlap	= Q.BlockOverlap;
	extract_type	= Q.extract_type;
	RealData		= Q.RealData;
	
// Memory properties
	lowmem			= Q.lowmem;
	reuse_input_data= Q.reuse_input_data;
	
// Filtering properties
	stat_import(Q);
}

void BFCurvelet3D_params::stat_import(BFCurvelet3D_params Q)
{
// Filtering properties
	no_norm			= Q.no_norm;
	no_fine			= Q.no_fine;
	no_coarse		= Q.no_coarse;
	threshold_coarse= Q.threshold_coarse;
	Positivity		= Q.Positivity;
	use_min			= Q.use_min;
	use_max			= Q.use_max;
	min_value		= Q.min_value;
	max_value		= Q.max_value;

	FilterType 		= Q.FilterType;
	SigmaNoise		= Q.SigmaNoise;
	NSigma			= Q.NSigma;
	Alpha			= Q.Alpha;
	WienerBS		= Q.WienerBS;
	Rho				= Q.Rho;
	Kmin			= Q.Kmin;
	Kmean			= Q.Kmean;
	Lmax			= Q.Lmax;
	force4sigma		= Q.force4sigma;	 
}

void BFCurvelet3D_params::print()
{
	cerr<<"Nx			    = "<<Nx 		 	 <<endl;
	cerr<<"Ny			    = "<<Ny 		 	 <<endl;
	cerr<<"Nz			    = "<<Nz 		 	 <<endl;
	cerr<<"NbrScale  	    = "<<NbrScale	 	 <<endl;
	cerr<<"NbrDir2d  	    = "<<NbrDir2d	 	 <<endl;
	cerr<<"BlockSize 	    = "<<BlockSize   	 <<endl;
	cerr<<"lapped		    = "<<lapped 	 	 <<endl;
	cerr<<"BlockOverlap     = "<<BlockOverlap	 <<endl;
	cerr<<"extract_type     = "<<extract_type	 <<endl;
	cerr<<"RealData  	    = "<<RealData	 	 <<endl;
	cerr<<" = "<<endl;
	cerr<<"lowmem		    = "<<lowmem 		 <<endl;
	cerr<<"reuse_input_data = "<<reuse_input_data<<endl;
	cerr<<" = "<<endl;
	cerr<<"no_norm  	    = "<<no_norm		 <<endl;
	cerr<<"no_fine  	    = "<<no_fine		 <<endl;
	cerr<<"no_coarse	    = "<<no_coarse  	 <<endl;
	cerr<<"threshold_coarse = "<<threshold_coarse<<endl;
	cerr<<"Positivity	    = "<<Positivity 	 <<endl;
	cerr<<"use_min  	    = "<<use_min		 <<endl;
	cerr<<"use_max  	    = "<<use_max		 <<endl;
	cerr<<"min_value	    = "<<min_value  	 <<endl;
	cerr<<"max_value	    = "<<max_value  	 <<endl;
	cerr<<" = "<<endl;
	cerr<<"FilterType	    = "<<FilterType 	 <<endl;
	cerr<<"SigmaNoise	    = "<<SigmaNoise 	 <<endl;
	cerr<<"NSigma		    = "<<NSigma 		 <<endl;
	cerr<<"Alpha		    = "<<Alpha  		 <<endl;
	cerr<<"WienerBS 	    = "<<WienerBS		 <<endl;
	cerr<<"Rho  		    = "<<Rho			 <<endl;
	cerr<<"Kmin 		    = "<<Kmin			 <<endl;
	cerr<<"Kmean		    = "<<Kmean  		 <<endl;
	cerr<<"Lmax 		    = "<<Lmax			 <<endl;
	cerr<<"force4sigma      = "<<force4sigma	 <<endl;
}

// ********************************************************************
// ********************************************************************
//				MAIN CLASS
// ********************************************************************
// ********************************************************************

BFCurvelet3D::BFCurvelet3D()
{
	reset_params();
	BSx=0;BSy=0;BSz=0;
	extern_norm=NULL;
	pointer=false;
	address=NULL;
}

BFCurvelet3D::~BFCurvelet3D()
{
	if(!lowmem) delete [] TabCur;
	else delete TabCur;
	if(TabStat != NULL) delete [] TabStat;
}

// ********************************************************************
void BFCurvelet3D::dealloc(fltarray*** &TabBand)
{
	int nb = lowmem ? 1 : nbr_block() ;
	for(int B=0;B<nb;B++)
	{
		for(int s=0;s<NbrScale;s++)
			delete [] TabBand[B][s];
		delete [] TabBand[B];
	}
	delete [] TabBand;
	TabBand = NULL;
}

// ********************************************************************

void BFCurvelet3D::set_tabsigmanoise(int s, float sig)
{
	if(!lowmem)
		for(int i=0;i<nbr_block();i++)
			TabCur[i].set_tabsigmanoise(s,sig);
	else
		TabCur->set_tabsigmanoise(s,sig);
}

// ********************************************************************

void BFCurvelet3D::alloc_from_coarse(fltarray Data, BFCurvelet3D_params &P)
{
	P.Nx = Data.nx();
	P.Ny = Data.ny();
	P.Nz = Data.nz();
	alloc_from_coarse(P);
// Automatic normalization
/*	if(!tabsigma)
	{
		estim_normalization(NULL);
		alloc_from_coarse(P);
	}
*/
}

// Nbr_Scale : nombre de direction total en 2D => NDir3DTotal = 6*(NbrDir/4)^2
void BFCurvelet3D::alloc_from_coarse(BFCurvelet3D_params &P)
{
	bool LocVerbose = true && Verbose ;
	if(Verbose) cerr<<"BFCurvelet3D::alloc_from_coarse()"<<endl;
	
// Read Parameters
	import_params(P);

// BlockSize check / B3D init
	if(BlockSize >= max(Nx,max(Ny,Nz)))
		BlockSize=0;
	if(BlockSize!=0)
	{
		if(BlockSize<25)
			BlockSize = 17;
		else if(BlockSize<45)
			BlockSize = 33;
		else BlockSize = 65;
	}
	if(BlockSize==0) { BSx=Nx; BSy=Ny; BSz=Nz;	}
	else { BSx=BlockSize; BSy=BlockSize; BSz=BlockSize;	}

	B3D.BlockOverlap = BlockOverlap;
	B3D.alloc(Nx,Ny,Nz,BlockSize);
	//cerr<<"B3D.size("<<nbr_block_nx()<<","<<nbr_block_ny()<<","<<nbr_block_nz()<<")"<<endl;
	
// Update parameters
//	Param.BlockSize = BlockSize;
//	P.BlockSize = BlockSize;
	
// TabCur allocation : classical way, memory usage
	if(!P.lowmem)
	{
		TabCur = new FCurvelet3D[nbr_block()];
		bool LVerbose = Verbose;
		for(int i=0;i<nbr_block();i++)
		{
			Verbose = LVerbose & (i==0);
			TabCur[i].set_positivity(P.Positivity);
			if(P.use_min) TabCur[i].set_min(P.min_value);
			if(P.use_max) TabCur[i].set_max(P.max_value);
			TabCur[i].set_no_coarse(P.no_coarse);
			TabCur[i].set_no_fine(P.no_fine);
			TabCur[i].set_threshold_coarse(P.threshold_coarse);
			if(no_norm) TabCur[i].set_no_norm();
			TabCur[i].alloc_from_coarse(P.NbrScale, BSx, BSy, BSz, P.NbrDir2d, False, False, P.RealData, P.extract_type);
		}
		Verbose = LVerbose;
	}
	// low memory way : only one curvelet saved
	else // lowmem
	{
		TabCur = new FCurvelet3D;
		TabCur[0].set_lowmem(P.BlockSize==0); // low memory on fct only if BS=0, else lowmem on blocks
		TabCur[0].set_positivity(P.Positivity);
		if(P.use_min) TabCur[0].set_min(P.min_value);
		if(P.use_max) TabCur[0].set_max(P.max_value);
		TabCur[0].set_no_coarse(P.no_coarse);
		TabCur[0].set_no_fine(P.no_fine);
		if(no_norm) TabCur[0].set_no_norm();
		TabCur[0].alloc_from_coarse(P.NbrScale,BSx, BSy, BSz, P.NbrDir2d, False, False, P.RealData, P.extract_type);
	}
	tabsigma = TabCur[0].isset_tabsigma();
	
// TabStat allocation
	TabStat = new dblarray[nbr_block()+1];
	for(int i=0;i<nbr_block()+1;i++)
		TabStat[i].alloc(P.NbrScale,nbr_band(0),6);
	
// Extern normalization initialisation
	if(extern_norm!=NULL)
		put_extern_norm();

	if(Verbose) cerr<<"...End BFCurvelet3D::alloc from coarse"<<endl;
}

// ********************************************************************

void BFCurvelet3D::transform(fltarray &Data, fltarray *** &TabBand, bool allocTB)
{
	//bool LocVerbose = false & Verbose;
	if (Verbose) cout << "BFCurvelet3D::transform... " << endl;
	
	if(!RealData)
	{
		cout<<"Compex data not taken into account"<<endl;
		exit(0);
	}
	if(allocTB)
	{
		// Creates the TabBand array
		TabBand = new fltarray**[nbr_block()];
		if(pointer)
		{
			// Allocates the vector, since the TB won't allocate anything
			address = new float[size_transform()];
//			cerr<<"Alloc pointer : "<<address<<":"<<&(address[size_transform()-1])<<", size="<<size_transform()<<endl;
		}
	}
	else if(pointer) // Sticks the address to the new TabBand to use.
		// The TB will be reallocated on address, thus must not change
		address = TabBand[0][0][0].buffer();
	
	bool LVerbose = Verbose;
	
if(lapped) B3D.fold(Data);

#if USE_OMP_BFC3D
omp_set_nested(1);
#pragma omp parallel for
#endif
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
//				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				fltarray CubeBlock(BSx,BSy,BSz);
				B3D.get_block_cube(Bi, Bj, Bk, Data, CubeBlock);
				int i=get_num_block(Bi,Bj,Bk);
//				Verbose = False;
				if(pointer)
				{
					int B = get_num_block(Bi,Bj,Bk);
					TabCur[i].set_pointer(get_pointer_of_block(B));
					//cerr<<"Pointer of block "<<B<<":"<<get_pointer_of_block(B)<<":"<<&((get_pointer_of_block(B)[TabCur[0].get_size_of_transform()-1]))<<endl;
				}
				TabCur[i].cur_trans(CubeBlock,TabBand[i],allocTB);
			}
	Verbose = LVerbose;
	
if(lapped) B3D.unfold(Data);
	
	Verbose = LVerbose;
	
	if (Verbose) cout << "...End BFCurvelet3D::transform" << endl;
}

void BFCurvelet3D::set_pointer(float* a)
{
	address=a;
	pointer=true;
//	cerr<<"Pointer:"<<a<<endl;
}
// It uses FCurvelet3D::get_single_wedge which sets the local TabBand correctly when allocTB is set
void BFCurvelet3D::transform(fltarray &Data, float * &out, bool alloc)
{
	if(Verbose) cerr<<"BFCurvelet3D::transform(f*)"<<endl;
	set_pointer(out);// activates pointer mode, and if !alloc makes loc_TB follow out/address
	if(!alloc) update_TabBand(out);
	transform(Data, local_TabBand, alloc);
	out = address;// updates out in case of alloc
	if(Verbose) cerr<<"BFCurvelet3D::transform(f*)"<<endl;
}
void BFCurvelet3D::recons(float * in, fltarray &out)
{
	update_TabBand(in);
	recons(local_TabBand, out);
}
void BFCurvelet3D::update_TabBand(float *in)
{
	set_pointer(in);
	for(int B=0;B<nbr_block();B++)
		TabCur[B].set_pointer(in);
	for(int B=0;B<nbr_block();B++)
		for(int s=0;s<NbrScale;s++)
			for(int b=0;b<nbr_band(s);b++)
				local_TabBand[B][s][b].alloc(get_pointer_of_band(B,s,b), size_band_nx(s,b), size_band_ny(s,b), size_band_nz(s,b));
}

// ********************************************************************

// Reconstruit dans Data, à partir de TabBand (curvelet), en passant par Tabcf_WT_Band, puis IFFT
void BFCurvelet3D::recons(fltarray *** TabBand, fltarray &Recons)
{
	//bool LocVerbose = true && Verbose;
	if (Verbose == True) cout << "BFCurvelet3D::recons()" << endl;
	
	Recons.alloc(Nx,Ny,Nz);
	Bool weight = BlockOverlap;
	
	bool LVerbose = Verbose;
	Verbose = False;

#if USE_OMP_BFC3D
omp_set_nested(1);
#pragma omp parallel for
#endif
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				fltarray CubeBlock(BSx,BSy,BSz);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].cur_recons(TabBand[i],CubeBlock);
				//add in case of overlapping
				B3D.add_block_cube(Bi, Bj, Bk, Recons, CubeBlock, weight);
			}
			
if(lapped) B3D.unfold(Recons);

	Verbose = LVerbose;
	
	if (Verbose == True) cout << "...End BFCurvelet3D::cur_recons" << endl;
}

// ********************************************************************

void BFCurvelet3D::get_band(int s, int b, fltarray ***TabBand, fltarray &band)
{
//	if(Verbose) cerr << "BFCurvelet3D::get_band("<<s<<","<<b<<".,.)..." <<  endl;
	//band=TabBand[s][b];
//	if(Verbose) cerr << "End BFCurvelet3D::get_band" <<  endl;
}

// ********************************************************************

void BFCurvelet3D::estim_normalization(char* Name_Imag_Out)
{
	if(Verbose) cerr << "BFCurvelet3D::estim_normalization("<<Name_Imag_Out<<")..."<< endl;
	
	// we must estimate the normalization coefficients on one block
	// at all scales, and half of one face of one block
	// the rest of the bands can be copied
	TabCur->estim_normalization(Name_Imag_Out);
	
// Output file
	char filename[256];
	int Ax, Ay, Az;
	if(BlockSize!=0)
	{
		Ax=BlockSize;
		Ay=BlockSize;
		Az=BlockSize;
	}
	else
	{
		Ax=Nx; Ay=Ny; Az=Nz;
	}
	if(Ax==Ay && Ax==Az)
		sprintf(filename,"%s/FastCur3D_norm/TabSigma_N%04d_n%02d_d%02d.fits", getenv("SAP_COM"), 
						Ax, NbrScale, NbrDir2d);
	else 
		sprintf(filename,"%s/FastCur3D_norm/TabSigma_N%04d_%04d_%04d_n%02d_d%02d.fits", getenv("SAP_COM"), 
						Ax, Ay, Az, NbrScale, NbrDir2d);
	Ifloat toto;
	toto.alloc(TabCur[0].TabSigma.buffer(), TabCur[0].TabSigma.ny(),TabCur[0].TabSigma.nx());
	io_write_ima_float(filename, toto);
	
	if(Verbose) cerr << "End BFCurvelet3D::estim_normalization"<<endl;
}

// ********************************************************************

void BFCurvelet3D::extract_stat(fltarray *** TabBand, char* Outname)
{
	//bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"BFCurvelet3D::Extract_stat..."<<endl;
	
// Output stat files
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_stat.dat");
	fstream cstat;
	cstat.open (Statname, fstream::out);
	int N=nbr_block();
	
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
	for(int Bj=0;Bj<nbr_block_ny();Bj++)
	for(int Bi=0;Bi<nbr_block_nx();Bi++)
	{
		int n=get_num_block(Bi,Bj,Bk);
		TabCur[n].extract_stat(TabBand[n],NULL, true); //no output, not centered
		TabStat[n]=TabCur[n].TabStat;
		TabStat[N]+=TabStat[n]; // Global statistics
		// Centering local stats and updating Min-Max
		for(int s=0;s<NbrScale;s++)
			for(int b=0;b<nbr_band(s);b++)
			{
				int Ntot = -1;
				moment4_center(Ntot, TabStat[n](s,b,0),TabStat[n](s,b,1),TabStat[n](s,b,2),TabStat[n](s,b,3),
						TabStat[n](s,b,1),TabStat[n](s,b,2),TabStat[n](s,b,3));
				if(2.*TabStat[n](s,b,4)<TabStat[N](s,b,4)) TabStat[N](s,b,4) = TabStat[n](s,b,4);
				else TabStat[N](s,b,4) -= TabStat[n](s,b,4);
				if(2.*TabStat[n](s,b,5)>TabStat[N](s,b,5)) TabStat[N](s,b,5) = TabStat[n](s,b,5);
				else TabStat[N](s,b,5) -= TabStat[n](s,b,5);
			}
	}
	// Centering global stats
	for(int s=0;s<NbrScale;s++)
		for(int b=0;b<nbr_band(s);b++)
		{
			int Ntot = -1 ;
			moment4_center(Ntot, TabStat[N](s,b,0)/N,TabStat[N](s,b,1)/N,TabStat[N](s,b,2)/N,TabStat[N](s,b,3)/N,
					TabStat[N](s,b,1),TabStat[N](s,b,2),TabStat[N](s,b,3));
			TabStat[N](s,b,0) /= N;
			if(Verbose) cerr << s <<"\t"<< b <<"\t"<< TabStat[N](s,b,0) <<"\t"<< TabStat[N](s,b,1) <<"\t"<< TabStat[N](s,b,2) <<"\t"<< TabStat[N](s,b,3) << endl;
			cstat << s <<"\t"<< b <<"\t"<< TabStat[N](s,b,0) <<"\t"<< TabStat[N](s,b,1) <<"\t"<< TabStat[N](s,b,2) <<"\t"<< TabStat[N](s,b,3) <<"\t"<<
					TabStat[N](s,b,4) << "\t" << TabStat[N](s,b,5) << endl;
		}
	
	cstat.close();
	if(Verbose) cerr<<"...End BFCurvelet3D::Extract_stat"<<endl;
}

// ********************************************************************
 
void BFCurvelet3D::threshold(fltarray ***TabBand, BFCurvelet3D_params &P)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"BFCurvelet3D::threshold(.,"<<P.SigmaNoise<<","<<P.NSigma<<","<<P.FilterType<<","<<P.force4sigma<<")"<<endl;
	
// Update filtering properties
	stat_import(P);
	
// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].threshold(TabBand[i], P.SigmaNoise, P.NSigma, P.FilterType, P.force4sigma);
			}
	Verbose = LVerbose;
	 
	if(Verbose) cerr<<"...End BFCurvelet3D::threshold"<<endl;
}
void BFCurvelet3D::soft_threshold(float* in, float lvl, bool th_coarse)
{
	if(Verbose) cerr<<"BFCurvelet3D::threshold(.,"<<lvl<<")"<<endl;
	
	FilterType = FT_SOFT;
	update_TabBand(in);

// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].soft_threshold(local_TabBand[i], lvl, th_coarse);
			}
	Verbose = LVerbose;
	
	if(Verbose) cerr<<"...End BFCurvelet3D::threshold"<<endl;
}
void BFCurvelet3D::hard_threshold(float* in, float lvl, bool th_coarse)
{
	if(Verbose) cerr<<"BFCurvelet3D::threshold(.,"<<lvl<<")"<<endl;
	
	FilterType = FT_SOFT;
	update_TabBand(in);

// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].hard_threshold(local_TabBand[i], lvl, th_coarse);
			}
	Verbose = LVerbose;
	
	if(Verbose) cerr<<"...End BFCurvelet3D::threshold"<<endl;
}

// ********************************************************************

void BFCurvelet3D::enhance(fltarray ***TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"BFCurvelet3D::enhance(.,"<<P.SigmaNoise<<","<<P.Rho<<","<<P.Kmin<<","<<P.Kmean<<","<<P.Lmax<<")"<<endl;
	
// Update filtering properties
	stat_import(P);
	
// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].enhance(TabBand[i], Rho, SigmaNoise*Kmin, SigmaNoise*Kmean, Lmax>0 ? Lmax : Lmax*SigmaNoise);
			}
	Verbose = LVerbose;
	
	if(Verbose) cerr<<"...End BFCurvelet3D::enhance"<<endl;
}

/*********************************************************************/

//void BFCurvelet3D::wiener(fltarray *** &TabBand, float noise_lvl, int WienerBS)
void BFCurvelet3D::wiener(fltarray *** &TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"BFCurvelet3D::wiener("<<P.SigmaNoise<<","<<P.WienerBS<<")..."<<endl;
	
// Update filtering properties
	stat_import(P);
	
// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].wiener(TabBand[i], P.SigmaNoise, P.WienerBS);
			}
	Verbose = LVerbose;

	if(Verbose) cerr<<"...End BFCurvelet3D::wiener"<<endl;
}

/*********************************************************************/

//void BFCurvelet3D::fdr(fltarray *** &TabBand, float Alpha, float SigmaNoise)
void BFCurvelet3D::fdr(fltarray *** &TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<P.Alpha<<","<<P.SigmaNoise<<")"<<endl;

// Update filtering properties
	stat_import(P);
	
// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].fdr(TabBand[i], P.Alpha, P.SigmaNoise);
			}
	Verbose = LVerbose;

	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

/****************************************************************************/

//void BFCurvelet3D::stein_block_threshold(fltarray *** &TabBand, float SigmaNoise)
void BFCurvelet3D::stein_block_threshold(fltarray *** &TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"BFCurvelet3D::stein_block_threshold(.)"<<endl;

// Update filtering properties
	stat_import(P);
	
// threshold
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].stein_block_threshold(TabBand[i], P.SigmaNoise);
			}
	Verbose = LVerbose;

	if(Verbose) cerr<<"...End BFCurvelet3D::stein_block_threshold"<<endl;
}

/*********************************************************************/

void BFCurvelet3D::put_extern_norm()
{
	for(int b=0;b<nbr_block();b++)
		TabCur[b].set_extern_norm(extern_norm[b]);
}

/*********************************************************************/

void BFCurvelet3D::filter(fltarray &Data, fltarray &Recons, BFCurvelet3D_params &P, char* Outname)
{
	if (Verbose) cout << "BFCurvelet3D::filter... " << endl;
	
// Update filtering properties
	stat_import(P);
	
// threshold
	char Statname[250];
	fstream cstat;
	if(Outname!=NULL && eval_stat)
	{
		strcpy(Statname, Outname);
		strcat(Statname, "_stat.dat");
		cstat.open (Statname, fstream::out);
	}
	
	bool LVerbose = Verbose;
	int N=nbr_block();
	
	if(reuse_input_data && !BlockOverlap) Recons.alloc(Data.buffer(),Data.nx(),Data.ny(),Data.nz());
	else Recons.alloc(Nx,Ny,Nz);
	
	Bool weight = BlockOverlap;
	
	Verbose = False;
	
// check wether to apply a curvelet-low-memory, or a block-low-memory
	if(nbr_block()==1)
	{
		fltarray CubeBlock(BSx,BSy,BSz);
		B3D.get_block_cube(0,0,0, Data, CubeBlock);
		TabCur[0].filter(CubeBlock, CubeBlock, SigmaNoise, NSigma, FilterType, Alpha, WienerBS, Rho, Kmin, Kmean, Lmax, force4sigma);
	}
#if USE_OMP_BFC3D
omp_set_nested(1);
#pragma omp parallel for
#endif
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
#if USE_OMP_BFC3D
#pragma omp parallel for
#endif
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
//cerr<<"Block("<<Bi<<","<<Bj<<","<<Bk<<")... lowmem ="<<TabCur[0].get_lowmem()<<endl;
				int n=get_num_block(Bi,Bj,Bk);
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				fltarray CubeBlock(BSx,BSy,BSz);
				B3D.get_block_cube(Bi, Bj, Bk, Data, CubeBlock);

				fltarray **TabBand=NULL;// = new fltarray*[1];
				FCurvelet3D* LocTabCur = new FCurvelet3D;
				LocTabCur[0].set_positivity(Positivity);
				if(use_min) LocTabCur[0].set_min(min_value);
				if(use_max) LocTabCur[0].set_max(max_value);
				LocTabCur[0].set_no_coarse(no_coarse);
				LocTabCur[0].set_no_fine(no_fine);
				LocTabCur[0].set_tabsigma(TabCur->TabSigma); // I think this is useless, but must be checked
				LocTabCur[0].alloc_from_coarse(NbrScale,BSx,BSy,BSz, NbrDir2d, False, False, RealData, extract_type);

				//LocTabCur[0].cur_trans(CubeBlock,TabBand,allocTB);
//cerr<<"Block("<<Bi<<","<<Bj<<","<<Bk<<") trans"<<endl;
				LocTabCur[0].cur_trans(CubeBlock,TabBand,true);
//cerr<<"Block("<<Bi<<","<<Bj<<","<<Bk<<") endtrans"<<endl;
//cerr<<CubeBlock(0,0,0)<<"->"<<TabBand[0][0](0,0,0);
				if(eval_stat)
				{
					LocTabCur[0].extract_stat(TabBand,NULL, true);// no output, not centered
					TabStat[n]=LocTabCur[0].TabStat;
					TabStat[N]+=TabStat[n]; // Global statistics
					// Centering local stats
					for(int s=0;s<NbrScale;s++)
						for(int b=0;b<nbr_band(s);b++)
						{
							moment4_center(-1,TabStat[n](s,b,0),TabStat[n](s,b,1),TabStat[n](s,b,2),TabStat[n](s,b,3),
											TabStat[n](s,b,1),TabStat[n](s,b,2),TabStat[n](s,b,3));
							if(2.*TabStat[n](s,b,4)<TabStat[N](s,b,4)) TabStat[N](s,b,4) = TabStat[n](s,b,4);
							else TabStat[N](s,b,4) -= TabStat[n](s,b,4);
							if(2.*TabStat[n](s,b,5)>TabStat[N](s,b,5)) TabStat[N](s,b,5) = TabStat[n](s,b,5);
							else TabStat[N](s,b,5) -= TabStat[n](s,b,5);
						}
				}
				if(FilterType==FT_HARD || FilterType==FT_SOFT) LocTabCur[0].threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma);
				else if(FilterType==FT_WIENER) LocTabCur[0].wiener(TabBand, SigmaNoise*NSigma/3., WienerBS);
				else if(FilterType==FT_FDR) LocTabCur[0].fdr(TabBand, Alpha, SigmaNoise);
				else if(FilterType==FT_SBT) LocTabCur[0].stein_block_threshold(TabBand, SigmaNoise);
				else if(FilterType==FT_SBT) LocTabCur[0].stein_block_threshold(TabBand, SigmaNoise);
				else if(FilterType==FT_CONTRAST) LocTabCur[0].enhance(TabBand, Rho, SigmaNoise*Kmin, SigmaNoise*Kmean, Lmax>0 ? Lmax : Lmax*SigmaNoise);

//cerr<<"Block("<<Bi<<","<<Bj<<","<<Bk<<") recons"<<endl;
//cerr<<"->"<<TabBand[0][0](0,0,0);
				LocTabCur[0].cur_recons(TabBand,CubeBlock);
//cerr<<"->"<<CubeBlock(0,0,0)<<" ";
//cerr<<"Block("<<Bi<<","<<Bj<<","<<Bk<<") end recons, del TB address="<<TabBand<<endl;
				delete LocTabCur;
				for(int i=0;i<NbrScale;i++)
					delete [] TabBand[i];
				delete [] TabBand;
				
				//add in case of overlapping
				if(BlockOverlap) B3D.add_block_cube(Bi, Bj, Bk, Recons, CubeBlock, weight);
				else B3D.put_block_cube(Bi, Bj, Bk, Recons, CubeBlock);
//cerr<<"  End Block "<<Bi<<","<<Bj<<","<<Bk<<endl;
			}
	// Centering global stats
	if(eval_stat)
	{
		for(int s=0;s<NbrScale;s++)
			for(int b=0;b<nbr_band(s);b++)
			{
				moment4_center(-1,TabStat[N](s,b,0)/N,TabStat[N](s,b,1)/N,TabStat[N](s,b,2)/N,TabStat[N](s,b,3)/N,
						TabStat[N](s,b,1),TabStat[N](s,b,2),TabStat[N](s,b,3));
				TabStat[N](s,b,0) /= N;
				if(Verbose)
					cerr << s <<"\t"<< b <<"\t"<< TabStat[N](s,b,0) <<"\t"<< TabStat[N](s,b,1) <<"\t"<< TabStat[N](s,b,2) <<"\t"<< TabStat[N](s,b,3) << endl;
				if(Outname!=NULL)
					cstat << s <<"\t"<< b <<"\t"<< TabStat[N](s,b,0) <<"\t"<< TabStat[N](s,b,1) <<"\t"<< TabStat[N](s,b,2) <<"\t"<< TabStat[N](s,b,3) <<"\t"<<
							TabStat[N](s,b,4) << "\t" << TabStat[N](s,b,5) << endl;
			}
	}
	
	Verbose = LVerbose;
	if(Outname!=NULL)
		cstat.close();
	
	if (Verbose) cout << "...End BFCurvelet3D::filter" << endl;
	
}

/*********************************************************************/

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
/*
void BFCurvelet3D::write(char *Name, fltarray *** TabBand)
{
	if(Verbose) cerr<<"BFCurvelet3D::write("<<Name<<",.,"<<")"<<endl;

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
//	long pcount   =   0;  // no group parameters 
//	long gcount   =   1;  // only a single image/group 
//	int  extend   =   False;

// write first header part (parameters)
	naxis=0;
	// Estimate the size of the transform
	int totalsize=0;
	for (int s=0; s < NbrScale; s++)
		for (int b=0; b < nbr_angle(s); b++)
		totalsize += TabBand[0][s][b].nx()*TabBand[0][s][b].ny()*TabBand[0][s][b].nz() ;
	totalsize *= nbr_block();
	
	naxis=1;
	naxes[0]=totalsize;
	naxes[1]=0;
	naxes[2]=0;
	
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D FastCurvelet", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrDir2d", (long) NbrDir2d, (char*)"Number of 2D directions (on a circle)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of the blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Overlaping blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"RealCur", (long) RealData, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) Nx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) Ny, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) Nz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	int cnt=1;
	for (int B=0;B<nbr_block();B++)
		for (int s=0; s < NbrScale; s++)
		{
			for (int b=0; b < nbr_angle(s); b++)
			{
				naxis=3;
				naxes[0] = TabBand[B][s][b].nx();
				naxes[1] = TabBand[B][s][b].ny();
				naxes[2] = TabBand[B][s][b].nz();

			// save the data
				if ( ffppre(fptr, group, cnt, naxes[0]*naxes[1]*naxes[2], (TabBand[B][s][b]).buffer(), &status) )
					PrintError( status );
				cnt += naxes[0]*naxes[1]*naxes[2];
			}
		}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	

	if(Verbose) cerr<<"...end BFCurvelet3D::write"<<endl;
}

// ********************************************************************

void BFCurvelet3D::read(char *Name, fltarray *** &TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
	bool LocVerbose = false & Verbose;
	
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
//	int _NbrScale, _Nx, _Ny, _Nz, _NbrDir2d, _BlockSize;
//	Bool _RealCur, _BlockOverlap;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	P.NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NbrDir2d", &mon_long, comment, &status)) PrintError( status );
	P.NbrDir2d = (int)mon_long;
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	P.BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	P.BlockOverlap = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"RealCur", &mon_long, comment, &status)) PrintError( status );
	P.RealData = (Bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Nz = (int)mon_long;
	
// Curvelet Alloc
	alloc_from_coarse(P);
		
// TabBand allocation
	TabBand = new fltarray**[nbr_block()];
	for(int B=0;B<nbr_block();B++)
	{
		TabBand[B] = new fltarray*[P.NbrScale];
		for(int s=0;s<P.NbrScale-1; s++)
		{
			TabBand[B][s] = new fltarray[nbr_angle(s)];
		}
		TabBand[B][P.NbrScale-1] = new fltarray[1];
	}

// Read data
	int cnt=1;
	for (int B=0;B<nbr_block();B++)
		for(int s=0;s<P.NbrScale;s++)
		{
			for (int b=0; b < nbr_angle(s); b++)
			{
				int NX,NY,NZ;
				NX = TabCur[B].size_band_nx(s,b);
				NY = TabCur[B].size_band_ny(s,b);
				NZ = TabCur[B].size_band_nz(s,b);
				if(LocVerbose) cerr<<" read TB("<<s<<","<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;

				TabBand[B][s][b].alloc(NX,NY,NZ);
				if (ffgpve(fptr, 1, cnt, NX*NY*NZ, nulval, (TabBand[B][s][b]).buffer(), &anynul, &status)) PrintError( status );
				cnt += NX*NY*NZ;
			}
		}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}

*/

void BFCurvelet3D::write(char *Name, fltarray *** TabBand)
{
	if(Verbose) cerr<<"BFCurvelet3D::write("<<Name<<",.,"<<")"<<endl;

	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3]; naxes[0]=0;naxes[1]=0;naxes[2]=0;
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
	naxis=1;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D FastCurvelet", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrDir2d", (long) NbrDir2d, (char*)"Number of 2D directions (on a circle)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of the blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Overlaping blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"RealCur", (long) RealData, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) Nx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) Ny, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) Nz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	for (int s=0; s < NbrScale; s++)
	{
		int cnt=1;
	// Size of the scale
		int totalsize=0;
		for (int b=0; b < nbr_angle(s); b++)
			totalsize += TabBand[0][s][b].nx()*TabBand[0][s][b].ny()*TabBand[0][s][b].nz() ;
		totalsize *= nbr_block();
		
		naxis=3;
		naxes[0] = totalsize;
		naxes[1] = 1;
		naxes[2] = 1;

		if(ffcrhd(fptr,&status))
			PrintError( status );
		if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
			PrintError( status );
		if ( ffpkyj(fptr, (char*)"NbrBands", (long) nbr_band(s), (char*)"Number of bands", &status))
			PrintError( status );  
		
		for (int B=0;B<nbr_block();B++)
			for (int b=0; b < nbr_angle(s); b++)
			{
				naxis=3;
				naxes[0] = TabBand[B][s][b].nx();
				naxes[1] = TabBand[B][s][b].ny();
				naxes[2] = TabBand[B][s][b].nz();

			// save the data
				if ( ffppre(fptr, group, cnt, naxes[0]*naxes[1]*naxes[2], (TabBand[B][s][b]).buffer(), &status) )
					PrintError( status );
				cnt += naxes[0]*naxes[1]*naxes[2];
			}
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  

	if(Verbose) cerr<<"...end BFCurvelet3D::write"<<endl;
}

// ********************************************************************

void BFCurvelet3D::read(char *Name, fltarray *** &TabBand, BFCurvelet3D_params &P)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
	bool LocVerbose = false & Verbose;
	
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
//	int _NbrScale, _Nx, _Ny, _Nz, _NbrDir2d, _BlockSize;
//	Bool _RealCur, _BlockOverlap;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	P.NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NbrDir2d", &mon_long, comment, &status)) PrintError( status );
	P.NbrDir2d = (int)mon_long;
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	P.BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	P.BlockOverlap = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"RealCur", &mon_long, comment, &status)) PrintError( status );
	P.RealData = (Bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	P.Nz = (int)mon_long;
	
// Curvelet Alloc
	alloc_from_coarse(P);
		
// TabBand allocation
	TabBand = new fltarray**[nbr_block()];
	for(int B=0;B<nbr_block();B++)
	{
		TabBand[B] = new fltarray*[P.NbrScale];
		for(int s=0;s<P.NbrScale-1; s++)
		{
			TabBand[B][s] = new fltarray[nbr_angle(s)];
		}
		TabBand[B][P.NbrScale-1] = new fltarray[1];
	}

// Read data
	for(int s=0;s<P.NbrScale;s++)
	{
		int cnt=1;
		if (fits_movabs_hdu(fptr, s+2, NULL, &status)) PrintError( status );
		for (int B=0;B<nbr_block();B++)
			for (int b=0; b < nbr_angle(s); b++)
			{
				int NX,NY,NZ;
				NX = TabCur[B].size_band_nx(s,b);
				NY = TabCur[B].size_band_ny(s,b);
				NZ = TabCur[B].size_band_nz(s,b);
				if(LocVerbose) cerr<<" read TB("<<s<<","<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;

				TabBand[B][s][b].alloc(NX,NY,NZ);
				if (ffgpve(fptr, 1, cnt, NX*NY*NZ, nulval, (TabBand[B][s][b]).buffer(), &anynul, &status)) PrintError( status );
				cnt += NX*NY*NZ;
			}
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}

// ********************************************************************

double BFCurvelet3D::get_max_coef(int s)
{
	double m=0;
	for(int b=0;b<nbr_band(s);b++)
	{
		m = max(m , max(abs(TabStat[nbr_block()](s,b,5)),abs(TabStat[nbr_block()](s,b,4))));
//		cerr<<abs(TabStat[nbr_block()](s,b,5))<<endl;
	}
	return m;
}

double BFCurvelet3D::get_max_coef()
{
	double m=0;
	for(int s=0;s<NbrScale-1;s++)
		m = max(m , get_max_coef(s));
	return m;
}

void BFCurvelet3D::temp(fltarray ***TabBand)
{
	for (int B=0;B<nbr_block();B++)
		TabCur->temp(TabBand[B]);
}

int BFCurvelet3D::size_transform()
{
	return TabCur[0].get_size_of_transform()*nbr_block();
}
float* BFCurvelet3D::get_pointer_of_band(int B, int scale, int band)
{
	return TabCur[B].get_pointer_of_band(scale, band);
}
float* BFCurvelet3D::get_pointer_of_block(int B)
{
	int a = TabCur[0].get_size_of_transform();
	return &(address[B*a]);
}

// ********************************************************************
// ********************************************************************
//					Global functions - for MEX
// ********************************************************************
// ********************************************************************

void fct3d_clear(vector< vector< vector< fltarray* > > > &C)
{// corresponds to BFCurvelet3D::transform and FCurvelet3D::cur_trans
	for(int i=0;i<C.size();i++)
	{
		for(int j=0;j<C[i].size();j++)
		{
			delete [] C[i][j][0];
			C[i][j].clear();
		}
		C[i].clear();
	}
	C.clear();
}	

int fct3d_transform(fltarray &Data, vector< vector< vector< fltarray* > > > &vTabBand, BFCurvelet3D_params &P)
{
//	P.print();
//fstream cx;
//cx.open ("log", fstream::out);
//cx<<"fdct3d_transform"<<endl;
	
// Parameters update
	P.Nx = Data.nx();
	P.Ny = Data.ny();
	P.Nz = Data.nz();
	
// FCurvelet initialisation
	BFCurvelet3D *DataC = new BFCurvelet3D;
	DataC->alloc_from_coarse(P);

// Noise calibration if not already done
	char filename[256];
	if(!DataC->isset_tabsigma())
	{
		fct3d_normalize(P, filename);
		// Update the coefficients for the current transform
		DataC->alloc_from_coarse(P);
	}
	
// Forward Transform
	fltarray*** TabBand;
	DataC->transform(Data, TabBand, true);

// vTabBand allocation
	vTabBand.resize(DataC->nbr_block());
	for(int B=0;B<DataC->nbr_block();B++)
	{
		vTabBand[B].resize(P.NbrScale);
		for(int s=0;s<P.NbrScale;s++)
		{
			vTabBand[B][s].resize(DataC->nbr_band(s));
			for(int b=0;b<DataC->nbr_band(s);b++)
			{
				vTabBand[B][s][b]=&TabBand[B][s][b];
				// We set the "e_GetBuffer" parameter of our arrays to prevent from deletion at 
				// the end of the function, as they will be covered by mxArrays
				// vTabBand[B][s][b].set_isbuffer(true);// doesn' seem to be enough
			}
		}
	}
	delete DataC;
	return 1;
}

int fct3d_recons(vector< vector< vector< fltarray* > > > &vTabBand, fltarray &Recons, BFCurvelet3D_params &P)
{
//	P.print();
//fstream cx;
//cx.open ("log", fstream::out);
//cx<<"fdct3d_recons"<<endl;

// Check params
	if(P.Nx==0 && P.Ny==0 && P.Nz==0)
		return -1; // "The size of the original data must be specified in the parameters.\nUse [C,P]=curvelet3d_transform(data,P) to update P."
	
// FCurvelet initialisation
	BFCurvelet3D *DataC = new BFCurvelet3D;
	DataC->alloc_from_coarse(P);

//cx<<"calib="<<!DataC->isset_tabsigma()<<endl;
	char filename[256];
	if(!DataC->isset_tabsigma())
		fct3d_normalize(P, filename);
	
// TabBand allocation
	fltarray*** TabBand;
	TabBand = new fltarray**[DataC->nbr_block()];
	for(int B=0;B<DataC->nbr_block();B++)
	{
		TabBand[B] = new fltarray*[P.NbrScale];
		for(int s=0;s<P.NbrScale;s++)
		{
			TabBand[B][s] = new fltarray[DataC->nbr_band(s)];
			for(int b=0;b<DataC->nbr_band(s);b++)
			{
				fltarray * A = vTabBand[B][s][b];
				TabBand[B][s][b].alloc(A->buffer(), A->nx(), A->ny(), A->nz());
			}
		}
	}
	
// Backward transform
	DataC->recons(TabBand, Recons);
	delete DataC;
	return 1;
}

int fct3d_filter(fltarray &Data, fltarray &Recons, BFCurvelet3D_params &P)
{
//	P.print();
// Parameters update
//	P.Nx = Data.nx();
//	P.Ny = Data.ny();
//	P.Nz = Data.nz();
	
// FCurvelet initialisation
	BFCurvelet3D *DataC = new BFCurvelet3D;
	DataC->alloc_from_coarse(P);

	char filename[256];
	if(!DataC->isset_tabsigma())
		fct3d_normalize(P, filename);
	
// FCurvelet filtering
	DataC->filter(Data,Recons,P, (char*)"");
	delete DataC;
	return 1;
}

int fct3d_threshold(vector< vector< vector< fltarray* > > > &vTabBand, BFCurvelet3D_params &P)
{
//	P.print();
// Check params
	if(P.Nx==0 && P.Ny==0 && P.Nz==0)
		return -1; // "The size of the original data must be specified in the parameters.\nUse [C,P]=curvelet3d_transform(data,P) to update P."
	
// FCurvelet initialisation
	BFCurvelet3D *DataC = new BFCurvelet3D;
	DataC->alloc_from_coarse(P);
	
// TabBand allocation
	fltarray*** TabBand;
	TabBand = new fltarray**[DataC->nbr_block()];
	for(int B=0;B<DataC->nbr_block();B++)
	{
		TabBand[B] = new fltarray*[P.NbrScale];
		for(int s=0;s<P.NbrScale;s++)
		{
			TabBand[B][s] = new fltarray[DataC->nbr_band(s)];
			for(int b=0;b<DataC->nbr_band(s);b++)
			{
				fltarray * A = vTabBand[B][s][b];
				TabBand[B][s][b].alloc(A->buffer(), A->nx(), A->ny(), A->nz());
			}
		}
	}
	
	if(P.FilterType==FT_HARD || P.FilterType==FT_SOFT)
		DataC->threshold(TabBand, P);
	else if(P.FilterType==FT_WIENER) DataC->wiener(TabBand, P);
	else if(P.FilterType==FT_FDR) DataC->fdr(TabBand, P);
	else if(P.FilterType==FT_SBT) DataC->stein_block_threshold(TabBand, P);
	else if(P.FilterType==FT_CONTRAST) DataC->enhance(TabBand, P);
	else cerr<<"Filtering method '"<<string_filter_type(P.FilterType)<<"'not implemented yet";
	
	delete DataC;
	return 1;
}

int fct3d_normalize(BFCurvelet3D_params &P, char filename[256])
{
//	P.print();
	P.Nx = (P.Nx/2)*2+1;
	P.Ny = (P.Ny/2)*2+1;
	P.Nz = (P.Nz/2)*2+1;
	fltarray Data(P.Nx,P.Ny,P.Nz);
	
// FCurvelet initialisation
	BFCurvelet3D *DataC = new BFCurvelet3D;
	DataC->alloc_from_coarse(P);
	
// Normalizing coefficients measurement
	DataC->estim_normalization(NULL);
	
	delete DataC;
	return 1;
}


















