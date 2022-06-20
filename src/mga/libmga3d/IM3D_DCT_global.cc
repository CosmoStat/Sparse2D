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
**    Date: 24th Oct. 2008
**    
**    File: IM3D_DCT.cc
**
**    Modification history:
**
******************************************************************************
**
**    DESCRIPTION 3D Block DCT transform and reconstruction
**    -----------  
**                 
******************************************************************************/

#include "IM3D_DCT.h"

extern Bool Verbose;

double TabSigma_BS_04 = 3.975;
double TabSigma_BS_08 = 9.531;
double TabSigma_BS_16 = 24.78;
double TabSigma_BS_32 = 67.03;
double TabSigma_BS_64 = 185.3;
double TabSigma_BS_128 = 517.2;

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

IM3D_DCT::IM3D_DCT()
{
	AllocClass=false;
	InitClass=false;
	reset();
}

IM3D_DCT::~IM3D_DCT()
{
	dealloc();
}

/****************************************************************************/

void IM3D_DCT::reset()
{
	if(AllocClass) dealloc();
	InitClass=false;
	NxCube=0;
	NyCube=0;
	NzCube=0;
	BlockSize=0;
	Border=I_MIRROR;
	Block1=NULL;
	Block3=NULL;
	dct_ip=NULL;
	dct_w=NULL;
	n1d=0;
	TabSigma = 0;
	lapped = false;
	skip_order=-1;

//pas en entrée de init
	//BlockOverlap=0;
}

/****************************************************************************/

void IM3D_DCT::alloc()
{
	if(AllocClass) dealloc();
	
// Temp Block allocation
	Block1 = new double[BlockSize*BlockSize*BlockSize];
	Block3 = new double**[BlockSize];
	for(int i=0;i<BlockSize;i++)
		Block3[i] = new double*[BlockSize];
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
			Block3[k][j] = &( Block1[ BlockSize*BlockSize*k + BlockSize*j ] );

// dct allocation
	n1d = 2 + (int) sqrt(BlockSize + 0.5);
	dct_ip = new int[n1d];
	dct_w = new double[BlockSize * 3 / 2];
}

/****************************************************************************/

void IM3D_DCT::dealloc()
{
	for(int k=0;k<BlockSize;k++)
		delete [] Block3[k];
	delete [] Block3;
	delete [] Block1;
	
	delete [] dct_ip;
	delete [] dct_w;
}

/****************************************************************************/

void IM3D_DCT::init(int Nx, int Ny, int Nz, int BS, Bool _BlockOverlap)
{
	if (Verbose) cout << "IM3D_DCT::init("<<Nx<<","<<Ny<<","<<Nz<<","<<BS<<","<<_BlockOverlap<<")" << endl;

	if(InitClass) reset();
	InitClass=true;

// Set parameters / Default parameters
	NxCube		= Nx;
	NyCube		= Ny;
	NzCube		= Nz;
	Border		= I_CONT;
	BlockOverlap= _BlockOverlap;
	if(BS==0) BlockSize = max(Nx,max(Ny,Nz));
	else BlockSize = BS > 0 ? (1<<(int)round(log2(BS))) : 8;

// Allocation
	alloc();

// Block init
	B3D.BlockOverlap = (Bool) BlockOverlap;
//cerr<<"Alloc B3D("<<Nx<<","<<Ny<<","<<Nz<<","<<BlockSize<<")"<<endl;
	B3D.BlockOverlap = (Bool) BlockOverlap;
	B3D.alloc(Nx,Ny,Nz,BlockSize);
//cerr<<"B3D.size("<<nbr_block_nx()<<","<<nbr_block_ny()<<","<<nbr_block_nz()<<")"<<endl;
	
//for(int i=0;i<BlockSize;i++)
//	cerr<<        "("<<2<<","<<i<<")="<<B3D.get_weight(i,2,15)<<
//			"\t("<<3<<","<<i<<")="<<B3D.get_weight(i,3,15)<<
//			"\t("<<4<<","<<i<<")="<<B3D.get_weight(i,4,15)<<
//			"\t("<<5<<","<<i<<")="<<B3D.get_weight(i,5,15)<<endl; 

// dct init
	dct_ip[0] = 0;

// TabSigma initialisation
	if(BlockSize==4) TabSigma = TabSigma_BS_04;
	else if(BlockSize==8) TabSigma = TabSigma_BS_08;
	else if(BlockSize==16) TabSigma = TabSigma_BS_16;
	else if(BlockSize==32) TabSigma = TabSigma_BS_32;
	else if(BlockSize==64) TabSigma = TabSigma_BS_64;
	else if(BlockSize==128) TabSigma = TabSigma_BS_128;
	else { TabSigma = 1; }
	
	if (Verbose) cerr << "...End IM3D_DCT::init" << endl;
}

/****************************************************************************/

void IM3D_DCT::transform_one_block(fltarray & Block, fltarray & BlockTrans)
{
	BlockTrans = Block;
	transform_one_block(BlockTrans);
}

/****************************************************************************/

extern void cosft2(float y[], int n, int isign);

void IM3D_DCT::transform_one_block(fltarray & Block)
{
	if(Verbose) cerr << "IM3D_DCT::transform_one_block..." <<  endl;

	float *d = new float[BlockSize];
	
// Transform on the three dimensions
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Block(i,j,k);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Block(i,j,k)=d[i];
		}

	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Block(j,i,k);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Block(j,i,k)=d[i];
		}

	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Block(k,j,i);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Block(k,j,i)=d[i];
		}
	
	delete [] d;
	
	if(Verbose) cerr << "...End IM3D_DCT::transform_one_block" <<  endl;
}

/****************************************************************************/

void IM3D_DCT::transform(fltarray &Cube, fltarray & CubeTrans)
{
	if(Verbose) cerr << "IM3D_DCT::transform..." <<  endl;
	
	CubeTrans.alloc(BlockSize*nbr_block_nx(),BlockSize*nbr_block_ny(),BlockSize*nbr_block_nz());
	
	fltarray Block(BlockSize,BlockSize,BlockSize);

	bool LVerbose = (bool) Verbose;
	Verbose = False;
	
	if(lapped)
		B3D.fold(Cube);

	//char filename[64];
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = (Bool) (LVerbose & (Bi==0) & (Bj==0) & (Bk==0));
				B3D.get_block_cube(Bi, Bj, Bk, Cube, Block);
//			sprintf(filename,"gblock%d%d%d.fits",Bi,Bj,Bk);
//			writefltarr(filename, Block);
				transform_one_block(Block);
				for(int i=0;i<BlockSize;i++)
				for(int j=0;j<BlockSize;j++)
				for(int k=0;k<BlockSize;k++)
					CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k) = Block(i,j,k);
//			sprintf(filename,"gtblock%d%d%d.fits",Bi,Bj,Bk);
//			writefltarr(filename, Block);
			}
	if(lapped)
	{
		B3D.unfold(Cube);
		writefltarr((char*)"unfold.fits", Cube);
	}
	
	Verbose = (Bool) LVerbose;
	if(Verbose) cerr << "...End IM3D_DCT::transform." <<  endl;
}

/****************************************************************************/

void IM3D_DCT::recons_one_block(fltarray & BlockTrans, fltarray & Block)
{
	Block = BlockTrans;
	recons_one_block(Block);
}

/****************************************************************************/

void IM3D_DCT::recons_one_block(fltarray &BlockTrans)
{
	if(Verbose) cerr << "IM3D_DCT::recons_one_block..." <<  endl;

	float *d = new float[BlockSize];

// Reconstruction on the three dimensions
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=BlockTrans(i,j,k);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				BlockTrans(i,j,k)=d[i];
		}

	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=BlockTrans(j,i,k);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				BlockTrans(j,i,k)=d[i];
		}

	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=BlockTrans(k,j,i);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				BlockTrans(k,j,i)=d[i];
		}

// Normalisation
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
				BlockTrans(i,j,k) /= (BlockSize*BlockSize*BlockSize/8);
	delete [] d;
	
	if(Verbose) cerr << "...End IM3D_DCT::recons_one_block" <<  endl;
}
/*
extern void ddct(int n, int isgn, double *a, int *ip, double *w);
extern void ddct2d (int, int, int, double **, double *, int *, double *);
void IM3D_DCT::recons_one_block(fltarray &BlockTrans)
{
	if(Verbose) cerr << "IM3D_DCT::recons_one_block..." <<  endl;

// Initialisation
	int cnt=0;
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
				Block1[cnt++] = BlockTrans(i,j,k);

// Reconstruction
	ddct3d(BlockSize,BlockSize,BlockSize,-1,Block3,NULL,dct_ip,dct_w);
	// 2d oomura : idem que cosft avec : ddct(BlockSize, -1, d, dct_ip,dct_w);
	
// Normalisation
	for(int i=0;i<BlockSize*BlockSize*BlockSize;i++)
		Block1[i]/=(BlockSize*BlockSize*BlockSize/8);
	for(int k=0;k<BlockSize;k++)
	{
		for(int i=0;i<BlockSize;i++)
			Block3[k][0][i] /= 2;
		for(int j=0;j<BlockSize;j++)
			Block3[k][j][0] /= 2;
	}
	for(int j=0;j<BlockSize;j++)
		for(int i=0;i<BlockSize;i++)
			Block3[0][j][i] /= 2;

//	// 2d oomura
//	for(int k=0;k<BlockSize;k++)
//		for(int j=0;j<BlockSize;j++)
//				BlockTrans(0,j,k) /= 2.;
//	for(int k=0;k<BlockSize;k++)
//		for(int j=0;j<BlockSize;j++)
//			for(int i=0;i<BlockSize;i++)
//				BlockTrans(i,j,k) /= BlockSize/2.;
	
// Recopy
	cnt=0;
	for(int k=0;k<BlockSize;k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
				BlockTrans(i,j,k) = Block1[cnt++];

	if(Verbose) cerr << "...End IM3D_DCT::recons_one_block" <<  endl;
}
*/
/****************************************************************************/

void IM3D_DCT::recons(fltarray &CubeTrans, fltarray &Recons)
{
	if(Verbose) cerr << "IM3D_DCT::recons..." <<  endl;

//	recons_one_block(TabBand,Cube);
	
	Recons.alloc(NxCube,NyCube,NzCube);
	
	fltarray Block(BlockSize,BlockSize,BlockSize);
//writefltarr((char*)"ctrans.fits", CubeTrans);
	bool LVerbose = (bool) Verbose;
	Verbose = False;
//char filename[64];
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = (Bool) (LVerbose & (Bi==0) & (Bj==0) & (Bk==0));
				for(int i=0;i<BlockSize;i++)
				for(int j=0;j<BlockSize;j++)
				for(int k=0;k<BlockSize;k++)
					Block(i,j,k) = CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k);
//			sprintf(filename,"tgtblock%d%d%d.fits",Bi,Bj,Bk);
//			writefltarr(filename, Block);
				recons_one_block(Block);
//			sprintf(filename,"tgblock%d%d%d.fits",Bi,Bj,Bk);
//			writefltarr(filename, Block);
				if(!BlockOverlap) B3D.put_block_cube(Bi, Bj, Bk, Recons, Block);
				else B3D.add_block_cube(Bi, Bj, Bk, Recons, Block, True);
			}
if(lapped)
B3D.unfold(Recons);
	Verbose = (Bool) LVerbose;

	if(Verbose) cerr << "...End IM3D_DCT::recons" <<  endl;
}

/****************************************************************************/

void IM3D_DCT::extract_stat(fltarray &TabBand, char* Outname)
{
//	bool LocVerbose = True & Verbose;
	if(Verbose) cerr<<"Extract_stat..."<<endl;

	if(Verbose) cerr<<"...End Extract_stat"<<endl;
}

/****************************************************************************/


void IM3D_DCT::normalize_self(fltarray TabBand, bool inverse)
{
	if(Verbose) cerr<<"IM3D_DCT::normalize_self(.,"<<inverse<<")"<<endl;
	for(int x=0;x<TabBand.nx();x++)
		for(int y=0;y<TabBand.ny();y++)
			for(int z=0;z<TabBand.nz();z++)
				if(inverse)
					TabBand(x,y,z) *= TabSigma;
				else
					TabBand(x,y,z) /= TabSigma;
	if(Verbose) cerr<<"end IM3D_DCT::normalize_self"<<endl;
}

/****************************************************************************/

void IM3D_DCT::threshold(fltarray &CubeTrans, float SigmaNoise, float NSigma, filter_type FilterType)
{
	if(Verbose) cerr<<"IM3D_DCT::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<")"<<endl;
	
	double lvl = SigmaNoise * NSigma * TabSigma;
	fltarray Block(BlockSize,BlockSize,BlockSize);
	float cnt=0,tot = BlockSize*BlockSize*BlockSize*nbr_block();
	
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				for(int k=0;k<BlockSize;k++)
					for(int j=0;j<BlockSize;j++)
						for(int i=0;i<BlockSize;i++)
							if(i!=0 | j!=0 | k!=0)
								if( abs(CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k)) < lvl )
								{
									cnt++;
									CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k)=0; // hard
								}
								else if(FilterType==FT_SOFT)
									CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k) -=
										(2*int(CubeTrans(Bi*BlockSize+i,Bj*BlockSize+j,Bk*BlockSize+k)>0)-1)*lvl;
			}
			
	if(skip_order>-1)
	{
		skip_order = min(skip_order,BlockSize-1);
		
		for(int Bk=0;Bk<nbr_block_nz();Bk++)
			for(int Bj=0;Bj<nbr_block_ny();Bj++)
				for(int Bi=0;Bi<nbr_block_nx();Bi++)
					for(int x=0;x<=skip_order;x++)
						for(int y=0;y<=skip_order-x;y++)
							for(int z=0;z<=skip_order-x-y;z++)
								CubeTrans(Bi*BlockSize+x,Bj*BlockSize+y,Bk*BlockSize+z)=0;
	}
	
	if(Verbose) cerr<<" n#proportion non seuillee ("<<lvl<<")="<<tot-cnt<<"#"<<(tot-cnt)/tot<<endl;
	
	if(Verbose) cerr<<"...End IM3D_DCT::threshold"<<endl;
}

/****************************************************************************/

void IM3D_DCT::wiener(fltarray &TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"IM3D_DCT::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
//	bool LocVerbose = true & Verbose;

	if(Verbose) cerr<<"...End IM3D_DCT::wiener"<<endl;
}

/****************************************************************************/

void IM3D_DCT::fdr(fltarray &TabBand, float Alpha, float SigmaNoise)
{
	if(Verbose) cerr<<"IM3D_DCT::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
	//bool LocVerbose = true & Verbose;
	
	if(Verbose) cerr<<"...End IM3D_DCT::fdr"<<endl;
}

/****************************************************************************/

void IM3D_DCT::write (char *Name, fltarray & TabBand, bool Normalize)
{
	if(Verbose) cerr<<"IM3D_DCT::write("<<Name<<",.)"<<endl;

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
	naxis=3;
	naxes[0] = TabBand.nx();
	naxes[1] = TabBand.ny();
	naxes[2] = TabBand.nz();
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
// write optional keyword to the header 
	if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D DCT", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) NxCube, (char*)"x size of the original cube", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) NyCube, (char*)"y size of the original cube", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) NzCube, (char*)"z size of the original cube", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of 3d blocks", &status))
		PrintError( status );  
	if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Overlaping blocks flag", &status))
		PrintError( status );  


// save the data
	if ( ffppre(fptr, group, 1, TabBand.nx()*TabBand.ny()*TabBand.nz(), (TabBand).buffer(), &status) )
		PrintError( status );

// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  

	if(Verbose) cerr<<"...end IM3D_DCT::write_multi"<<endl;
}

/*********************************************************************/

void IM3D_DCT::read(char *Name, fltarray & TabBand, bool *NormalizeInv)
{
	if(Verbose) cerr<<"IM3D_DCT::read("<<Name<<",.,.,.)"<<endl;
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
	int _BlockSize,_Nx,_Ny,_Nz;
	Bool _BlockOverlap;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	_BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	_BlockOverlap = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*NormalizeInv = (bool)mon_long;

//	init the structure 
	init(_Nx, _Ny, _Nz, _BlockSize, _BlockOverlap);

// tabband create
	TabBand.alloc(BlockSize*nbr_block_nx(),BlockSize*nbr_block_ny(),BlockSize*nbr_block_nz());

// read data
//	int NX,NY,NZ;

	if (ffgpve(fptr, 1, 1, BlockSize*nbr_block_nx()*BlockSize*nbr_block_ny()*BlockSize*nbr_block_nz(), 
					nulval, TabBand.buffer(), &anynul, &status)) PrintError( status );

// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );

	if(Verbose) cerr<<"...end IM3D_DCT::read_multi"<<endl;
}

/****************************************************************************/

void IM3D_DCT::temp(fltarray &TabBand)
{
	
}


