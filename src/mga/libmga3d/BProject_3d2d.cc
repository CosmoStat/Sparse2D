/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  27/05/2009
**    
**    File:  BProject_3d2d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  see .h
**    -----------  
**
******************************************************************************/

#include "BProject_3d2d.h"

extern bool Verbose;

BProject_3d2d::BProject_3d2d()
{
	no_coarse=false;
	no_fine=false;
	positivity=false;
	use_min=false; use_max=false;
	min_value=0; max_value=0;
	eval_stat=false;
	TabStat=NULL;
	threshold_coarse=false;
	Project=NULL;
}

BProject_3d2d::~BProject_3d2d()
{
}

// ********************************************************************

void BProject_3d2d::dealloc(Ifloat*** &TabBand)
{
/*	for(int B=0;B<nbr_block();B++)
	{
		for(int s=0;s<NbrScale;s++)
			delete [] TabBand[B][s];
		delete [] TabBand[B];
	}
*/	delete [] TabBand;
	TabBand = NULL;
}

// ********************************************************************

void BProject_3d2d::set_threshold_coarse(bool th)
{
//	for(int i=0;i<nbr_block();i++)
//		TabCur[i].set_threshold_coarse(th);
}

// ********************************************************************

// Nbr_Scale : nombre de direction total en 2D => NDir3DTotal = 6*(NbrDir/4)^2
void BProject_3d2d::alloc_from_coarse(int _NbrScale, int Nx, int Ny, int Nz, int _NbrDir2d, Bool _BlockOverlap)
{
	bool LocVerbose = true && Verbose ;
	if(Verbose) cerr<<"BProject_3d2d::alloc_from_coarse("<<_NbrScale<<","<<Nx<<","<<Ny<<","<<_NbrDir2d<<","<<_BlockOverlap<<")..."<<endl;
	
	NbrScale=_NbrScale;
	DataNx = Nx;
	DataNy = Ny;
	DataNz = Nz;
	NbrDir2d = max(_NbrDir2d,12);
	BlockOverlap=_BlockOverlap;
cerr<<Nx<<" "<<Ny<<" "<<Nz<<endl;
	if(Nz<33)
		BlockSize = 17;
	else if(Nz<65)
		BlockSize = 33;
	else 
	{
		cerr<<"too many frames"<<endl;
		exit(0);
	}
BlockSize = 17;
	BSx=BlockSize;
	BSy=BlockSize;
	BSz=BlockSize;

	B3D.BlockOverlap = (Bool) BlockOverlap;
	B3D.alloc(Nx,Ny,Nz,BlockSize);
cerr<<"BlockSize "<<BlockSize<<" B3D.size("<<nbr_block_nx()<<","<<nbr_block_ny()<<","<<nbr_block_nz()<<")="<<nbr_block()<<endl;
TabFrame = new Ifloat[nbr_block()];
	
// TabCur allocation : classical way, memory usage
	TabCur = new FCUR[nbr_block()];
	bool LVerbose = Verbose;
	for(int i=0;i<nbr_block();i++)
	{
		Verbose = LVerbose & (i==0);
//		TabCur[i].set_positivity(positivity);
//		if(use_min) TabCur[i].set_min(min_value);
//		if(use_max) TabCur[i].set_max(max_value);
//		TabCur[i].set_no_coarse(no_coarse);
//		TabCur[i].set_no_fine(no_fine);
		TabCur[i].alloc_from_coarse(NbrScale, BSx, BSy, NbrDir2d, False, False, True);
	}
	Verbose = LVerbose;
//	tabsigma = TabCur[0].isset_tabsigma();
	
// Project_3d2d allocation
	Project = new Project_3d2d;

// TabStat allocation
//	TabStat = new dblarray[nbr_block()+1];
//	for(int i=0;i<nbr_block()+1;i++)
//		TabStat[i].alloc(NbrScale,nbr_band(0),6);
	
// Extern normalization initialisation
//	if(extern_norm!=NULL)
//		put_extern_norm();

	if(Verbose) cerr<<"...End BProject_3d2d::alloc from coarse"<<endl;
}

// ********************************************************************

void BProject_3d2d::transform(fltarray &Data, Ifloat *** &TabBand, bool allocTB)
{
	//bool LocVerbose = false & Verbose;
	if (Verbose) cout << "BProject_3d2d::transform... " << endl;
cerr<<"Y";
	
	if(allocTB)
		TabBand = new Ifloat**[nbr_block()];
	
	bool LVerbose = Verbose;
	
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
//				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				fltarray CubeBlock(BSx,BSy,BSz);
cerr<<"Y";
				B3D.get_block_cube(Bi, Bj, Bk, Data, CubeBlock);
				int i=get_num_block(Bi,Bj,Bk);
				Verbose = False;

Project->setTakeOnlyCompletePlane(False);
Project->setDoNotCorrectSigma(True);// False ?
Project->setParam(BlockSize,NbrScale);
cerr<<"X";
Project->alloc();
cerr<<"B";
fltarray Partial(BlockSize, BlockSize, Project->getNbPlane());
Project->transform (CubeBlock, Partial);
cerr<<"A";
				Ifloat Frame(BlockSize, BlockSize);
				for(int i=0;i<BlockSize;i++)
				for(int j=0;j<BlockSize;j++)
					Frame(i,j) = Partial(i,j,0);
TabFrame[i] = Frame;
				TabCur[i].cur_trans(Frame);
cerr<<"C";
				
			// Saving to TabBand
				if(allocTB) TabBand[i] = new Ifloat *[NbrScale];
				for(int s=0;s<NbrScale;s++)
				{
					int NbrBands=nbr_bands(s);
					if(allocTB) TabBand[i][s] = new Ifloat [NbrBands];

					for(int b=0;b<NbrBands;b++)
					{
						TabBand[i][s][b].alloc(TabCur->size_band_nl(s,b),TabCur->size_band_nc(s,b));
						TabCur[i].get_band(s, b, TabBand[i][s][b]);
					}
				}
			}
	Verbose = LVerbose;
	
	if (Verbose) cout << "...End BProject_3d2d::transform" << endl;
}

// ********************************************************************

// Reconstruit dans Data, à partir de TabCF_Band (curvelet), en passant par Tabcf_WT_Band, puis IFFT
void BProject_3d2d::recons(Ifloat *** TabBand, fltarray &Recons)
{
	//bool LocVerbose = true && Verbose;
	if (Verbose == True) cout << "BProject_3d2d::recons()..." << endl;
cerr<<DataNz<<endl;
	Recons.alloc(DataNx,DataNy,DataNz);
	Bool weight = (Bool) BlockOverlap;
	
	bool LVerbose = Verbose;
	Verbose = False;

	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
cerr<<"EE";
//				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				fltarray CubeBlock(BSx,BSy,BSz);
				int i=get_num_block(Bi,Bj,Bk);
				
				for(int s=0;s<NbrScale;s++)
					for(int b=0;b<nbr_bands(s);b++)
						TabCur[i].put_band(s, b, TabBand[i][s][b]);

cerr<<"DD";
				Ifloat Frame;
				TabCur[i].cur_recons(Frame);

Frame = TabFrame[i];



Project->setTakeOnlyCompletePlane(False);
Project->setDoNotCorrectSigma(True);
Project->setParam(BlockSize,NbrScale);
Project->alloc();
fltarray Partial(BlockSize, BlockSize, Project->getNbPlane());
cerr<<"HH";
for(int i=0;i<BlockSize;i++)
for(int j=0;j<BlockSize;j++)
	Partial(i,j,0) = Frame(i,j);
cerr<<"GG";
Project->recons(Partial, CubeBlock);
//cerr<<"FF"<<endl;
				
			// = Project_3d2d->recons(CubeBlock)
			//add in case of overlapping
//cerr<<Recons.nx()<<Recons.ny()<<Recons.nz()<<","<<CubeBlock.nx()<<CubeBlock.ny()<<CubeBlock.nz()<<endl;
				B3D.add_block_cube(Bi, Bj, Bk, Recons, CubeBlock, weight);
			}
	Verbose = LVerbose;
	
	if (Verbose == True) cout << "...End BProject_3d2d::cur_recons" << endl;
}

// ********************************************************************

/*
void BProject_3d2d::get_band(int s, int b, Ifloat ***TabBand, fltarray &band)
{
//	if(Verbose) cerr << "BProject_3d2d::get_band("<<s<<","<<b<<".,.)..." <<  endl;
	//band=TabBand[s][b];
//	if(Verbose) cerr << "End BProject_3d2d::get_band" <<  endl;
}
*/
// ********************************************************************

/*void BProject_3d2d::get_norm_coeff(float N_Sigma, fltarray*** &TabBand, char* Name_Imag_Out, bool &allocTB)
{
	if(Verbose) cerr << "BProject_3d2d::get_norm_coeff("<<N_Sigma<<",.,"<<allocTB<<")..."<< endl;

	if(allocTB)
		TabBand = new fltarray**[nbr_block()];

	bool LVerbose = Verbose;
	Verbose = False;
	if(!lowmem) 
	{
		for(int Bk=0;Bk<nbr_block_nz();Bk++)
			for(int Bj=0;Bj<nbr_block_ny();Bj++)
				for(int Bi=0;Bi<nbr_block_nx();Bi++)
				{
					Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
					bool Lalloc = allocTB;
					int i=get_num_block(Bi,Bj,Bk);
					TabCur[i].get_norm_coeff(N_Sigma, TabBand[i], Name_Imag_Out, Lalloc);
				}
	}
	else
	{
		bool Lalloc = allocTB;
		TabCur[0].get_norm_coeff(N_Sigma, TabBand[0], Name_Imag_Out, Lalloc);
	}
	Verbose = LVerbose;
	allocTB = false;
	
	if(Verbose) cerr << "...End BProject_3d2d::get_norm_coeff" << endl;
}*/

// ********************************************************************

// The same as extract_stat except that there is no normalisation 
/*void BProject_3d2d::noise_calibration(fltarray *** TabBand, char* Outname)
{
	if(Verbose) cerr<<"BProject_3d2d::Noise_calibration(.,"<<Outname<<")..."<<endl;

	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].noise_calibration(TabBand[i], Outname);
			}
	
	if(Verbose) cerr<<"...End BProject_3d2d::Noise_calibration"<<endl;
}*/

// ********************************************************************
/*
void BProject_3d2d::extract_stat(fltarray *** TabBand, char* Outname)
{
	//bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"BProject_3d2d::Extract_stat..."<<endl;
	
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
	if(Verbose) cerr<<"...End BProject_3d2d::Extract_stat"<<endl;
}
*/
// ********************************************************************
/*
void BProject_3d2d::normalize_self(fltarray ***TabBand, bool inverse)
{
	if(Verbose) cerr<<"BProject_3d2d::normalize_self(."<<inverse<<")"<<endl;
	
	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].normalize_self(TabBand[i], inverse);
			}
	Verbose = LVerbose;

	if(Verbose) cerr<<"...End BProject_3d2d::normalize"<<endl;
}
*/
// ********************************************************************
/*
void BProject_3d2d::threshold(fltarray ***TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"BProject_3d2d::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;

	bool LVerbose = Verbose;
	Verbose = False;
	for(int Bk=0;Bk<nbr_block_nz();Bk++)
		for(int Bj=0;Bj<nbr_block_ny();Bj++)
			for(int Bi=0;Bi<nbr_block_nx();Bi++)
			{
				Verbose = LVerbose & (Bi==0) & (Bj==0) & (Bk==0);
				int i=get_num_block(Bi,Bj,Bk);
				TabCur[i].threshold(TabBand[i], SigmaNoise, NSigma, FilterType, force4);
			}
	Verbose = LVerbose;
	
	if(Verbose) cerr<<"...End BProject_3d2d::threshold"<<endl;
}
*/
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
void BProject_3d2d::write(char *Name, fltarray *** TabBand, bool Normalize)
{
	if(Verbose) cerr<<"BProject_3d2d::write("<<Name<<",.,"<<Normalize<<")"<<endl;

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
		if ( ffpkyj(fptr, (char*)"Normaliz", (long) Normalize, (char*)"1 if the transform is normalized, else 0", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrDir2d", (long) NbrDir2d, (char*)"Number of 2D directions (on a circle)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"BlkSize", (long) BlockSize, (char*)"Size of the blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Overlap", (long) BlockOverlap, (char*)"Overlaping blocks", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrBands", (long) nbr_tot_band(), (char*)"Total number of bands", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"RealCur", (long) RealBand, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) DataNx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) DataNy, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) DataNz, (char*)"z size of the original cube", &status))
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
	//	cerr<<"save ("<<s<<","<<b<<") "<<naxes[0]<<" "<<naxes[1]<<" "<<naxes[2]<<endl;

			// save the data
				if ( ffppre(fptr, group, cnt, naxes[0]*naxes[1]*naxes[2], (TabBand[B][s][b]).buffer(), &status) )
					PrintError( status );
				cnt += naxes[0]*naxes[1]*naxes[2];
	//			if ( fits_write_subset(fptr, group, ) )
	//				PrintError( status );

	//	cerr<<"saved"<<endl;
			}
		}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	

	if(Verbose) cerr<<"...end BProject_3d2d::write"<<endl;
}
*/
// ********************************************************************
/*
void BProject_3d2d::read(char *Name, fltarray *** &TabBand, bool *NormalizeInv)
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
	int _NbrScale, _Nx, _Ny, _Nz, _NbrDir2d, _BlockSize;
	Bool _RealCur, _BlockOverlap;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"Normaliz", &mon_long, comment, &status)) PrintError( status );
	*NormalizeInv = (bool)mon_long;
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NbrDir2d", &mon_long, comment, &status)) PrintError( status );
	_NbrDir2d = (int)mon_long;
	if (ffgkyj(fptr,(char*)"BlkSize", &mon_long, comment, &status)) PrintError( status );
	_BlockSize = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
	_BlockOverlap = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"RealCur", &mon_long, comment, &status)) PrintError( status );
	_RealCur = (Bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
//	init the structure
	alloc_from_coarse(_NbrScale,_Nx, _Ny, _Nz, _NbrDir2d, _BlockSize, _BlockOverlap, _RealCur);
	
// TabBand allocation
	TabBand = new fltarray**[nbr_block()];
	for(int B=0;B<nbr_block();B++)
	{
		TabBand[B] = new fltarray*[NbrScale];
		for(int s=0;s<NbrScale-1; s++)
		{
			TabBand[B][s] = new fltarray[nbr_angle(s)];
		}
		TabBand[B][NbrScale-1] = new fltarray[1];
	}

// Read data
	int cnt=1;
	for (int B=0;B<nbr_block();B++)
		for(int s=0;s<_NbrScale;s++)
		{
			for (int b=0; b < nbr_angle(s); b++)
			{
				int NX,NY,NZ;
				NX = TabCur[B].size_band_nx(s,b);
				NY = TabCur[B].size_band_ny(s,b);
				NZ = TabCur[B].size_band_nz(s,b);
				if(LocVerbose) cerr<<" read TB("<<s<<","<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;

				TabBand[B][s][b].alloc(NX,NY,NZ);
				TabCur[B].TabCF_Band[s][b].alloc(NX,NY,NZ);
				if (ffgpve(fptr, 1, cnt, NX*NY*NZ, nulval, (TabBand[B][s][b]).buffer(), &anynul, &status)) PrintError( status );
				cnt += NX*NY*NZ;
			}
		}
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}
*/
// ********************************************************************

