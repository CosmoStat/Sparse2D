/******************************************************************************
**                   Copyright (C) 2001 by CEA 
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe
**
**    Date:  25/08/2001 
**    
**    File:  Linelet3D3D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the 3D Linelet3D transform
**    ----------- 
******************************************************************************/

#ifndef _LINELET3D_H_
#define _LINELET3D_H_

#ifdef WINDOWS
	#define USE_OMP_L 0
#else
	#define USE_OMP_L 1
#endif

#if USE_OMP_L
	#include <omp.h>
#endif

#include "IM_Obj.h"
#include "IM3D_PartialRadon.h"
#include "IM3D_IO.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM3D_Block.h"
#include "MGA_Inc.h"

enum type_linelet3d_WTclass{LIN3D_CL_UNKNOW, 
							LIN3D_CL_PYR,
							LIN3D_CL_ORTHO};

#define NBR_LIN3D_TRANS 2
enum type_linelet3d_WTtrans{LIN3D_UNKNOWN,
							LIN3D_OWT,       // standard bi-orthogonal WT
							LIN3D_PYR_FFT};  // Pyramidal FFT transform
#define DEF_LIN3D_TRANS LIN3D_OWT

/*****************************************************************************/

const char * StringLin3DTransform (type_linelet3d_WTtrans  type);
type_linelet3d_WTclass lin3d_class(type_linelet3d_WTtrans Type);
void linelet3d_check_size_cube(fltarray &Cube, int & BlockSize);

/*****************************************************************************/

class Linelet3D {

/*****************************************************************************/
// private attributs
/*****************************************************************************/
private: 

	Bool AllocClass;
   // general
   Bool  _StatInfo;      // If true, statistical information of each band 
                         // is printed (default is no).
   fltarray* _pScale;  
   int   _NbPlaneUsedInStat; 
   int   _NbPlaneUsedInStatByBlock; 
   Bool  _TakeOnlyCompletePlane;	
   Bool  _DoNotCorrectSigma;
   Bool  _Verbose;       // Verbose mode  (default is no).
   Bool  _Control;       // control param in, attribut....
   Bool  _WriteFile;     // write all intermrdiary files
   type_linelet3d_WTtrans _LinTransf;
   SubBand1D* _ptrSB1D;  // pointer to SB1D for orthogonal transform
   
   // image in
   int   _Nx;            // size of cube in x direction
   int   _Ny;            // size of cube in y direction
   int   _Nz;            // size of cube in z direction
   
   // block info
   int   _BlockSize;     // sise of decomposition block
   Bool  _BlockOverlap;  // If True, Overlapped blocks are used.
   float _OverlapFactor; // Overlapping factor. Between 0 and 0.5
   int   _NbrBlock;      // number of block 
   int   _Nxb;           // Number of blocks in x direction
   int   _Nyb;           // Number of blocks in y direction
   int   _Nzb;           // Number of blocks in z direction  
    
   // linelet transf info
   int    NbrBand;       // Number of band
   int   _NbrScale;      // sacle number in WT
   int   *TabBandNx;     // Size of the bands on the x-axis
   int   *TabBandNy;     // Size of the bands on the y-axis
   int   *TabBandNz;     // Size of the bands on the z-axis
   int   *TabDepx;       // position of the bands on the x-axis
   int   *TabDepy;       // position of the bands on the y-axis
   int   *TabDepz;       // position of the bands on the z-axis 
     
   int   _RadNbPlane;    // plane number in radon transf = 3*_SizeCube^2
   int   _LinNbPlane;    // plane number in linelet transf 
   int   _LinPlaneSize;  // sise of each plane
   
   // struct used in transf one block and recons one block
   fltarray _PRadTransf;
   
   // partial radon transform
   PartialRadon3D PRad3D; 
   
   
	bool  no_recons;
	bool keep_energy;

   
/*****************************************************************************/
// private members
/*****************************************************************************/
private: 
  Block3D B3D;       // Class for 3D block management
		
  void reset();         // Reset all internal variables
  
  void set_tab_band();
                        // Initialize the size of each band.
			// TabBandNx, TabBandNy, TabBandNz
			  
  void set_recons_block_param();
                        // initialize the variables for the block transform 
                        // from the input cube size
                        // this routine is automatically called from transform
                        // recons, and filtering routines 
			
			  
  void block_transform (fltarray& Cube, fltarray& Transf);
                        // Linelet3D transform of a cube,  block per block  
  
  void sim_transform_one_block (fltarray& Cube, fltarray& Transf); 
  void transform_one_block (fltarray& Cube, fltarray& Transf);
                       // Linelet3D transform of a cube 		

  void  get_block_cube (int Bi, int Bj, int Bk, fltarray& Cube, 
                       fltarray& CurrentBlock);
                       // Extract a block from the cube			 
  
  void put_block_trans(int Bi, int Bj, int Bk,  fltarray& Transf,  
                       fltarray& CurrentBlock);
                       // Put a block from the ridgelet transformed image
                       // Transf and ImaTransBlock must be already allocated  
		    
  void block_recons (fltarray& Transf, fltarray& Cube);
                       // Cube reconstruction from its linelet transform
                       // Block per block	
		     	
  void get_block_trans (int Bi, int Bj, int Bk, fltarray& Transf, fltarray& Cube);
                       // Get a block from the linelet transformed cube 
                       // Transf and ImaTransBlock must be already allocated

  void sim_recons_one_block (fltarray& Transf, fltarray& Cube);		       
  void recons_one_block (fltarray& Transf, fltarray& Cube);
                       // Cube reconstruction from the linelet transform
		       
  void put_block_cube (int Bi, int Bj, int Bk, fltarray& Cube, 
                       fltarray& CubeBlock);  
  void add_block_cube (int Bi, int Bj, int Bk, fltarray& Cube, 
                       fltarray& CubeBlock);
                       // Put or Add (if OverLap) a block to a Cube
		       
  float get_weight(int Bi, int Bj, int Bk, int i, int j, int k);	
                       // return the weight value for pixel position i,j,k
		       // in the block Bi,Bj,Bk
  float lap_weight(int PosPix, int CurrentBlock, int MaxBlockNumber);
  
  void stat_alloc ();
  void stat_add (int Bi, int Bj, int Bk, fltarray& CubeBlock);
  void mr_io_fill_header(fitsfile *fptr);
  public:
  void mr_io_read (char *Name, fltarray& Cube);
  void mr_io_write (char* Name, fltarray& Cube);
	       	
		       	       		       			
/*****************************************************************************/
// public attributs
//*****************************************************************************/
public: 

/*****************************************************************************/
// public members
/*****************************************************************************/
public: 

	Linelet3D();
	~Linelet3D();

	void set_attribut(fltarray& Cube, int& BlockSize);
	void set_attribut(int nx,int ny,int nz, int & BlockSize);
                    	// set attribut before computing linelet...

	void set_block_param();
                    	// initialize the variables for the block transform 
                    	// from the input cube size
                    	// this routine is automatically called from transform
                    	// recons, and filtering routines




	void transform (fltarray& Cube, fltarray& Transf, int BlockSize);
                    	// Linelet3D transform of a cube whith NbrScale
                    	// aplied block per block
                    	// if BlockSize=0, BlockSize=cube size   

	void alloc(int NxCube, int NyCube, int NzCube, int BlockSize);
                	   // initialize the Class for a given cube siwe
			   // and a given block size
	void dealloc();

	void transform(fltarray &Cube, fltarray &Transf);
                	   // Linelet3D transform of a cube
			   // the class should have been first allocated
			   // with the alloc routine

	void recons (fltarray& Trans, fltarray& Cube);
                    	// Cube reconstruction from the ridgelet transform
                    	// applied block per block

	void write (char* Name, fltarray& Cube); 
                    	// write result on File name Name 

	void trace_param();  // trace in param 

	void info_stat (fltarray &Transf, char* FileSimu, Bool WithFileSimu=False);
                    	// print some statistical information, band per band
                    	// of the linelet transform   
	// Filtering methods
	void block_filter (fltarray& Cube, fltarray& Recons);
//	void filter(fltarray& Cube, fltarray& Recons, int BlockSize);
	void filter(fltarray &Cube, fltarray &Recons, int BlockSize, 
		int s3, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, float Alpha=0.01, int WienerBS=3, bool force4=false, bool UseCubeSigma=false, 
		fltarray *TabSigma=NULL, fltarray *CubeSigma=NULL, dblarray * TabStat=NULL, int NbrScale3D=0);
	void threshold_single(fltarray &Band, int BlockSize, 
		int s3, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4=false, bool UseCubeSigma=false, 
		fltarray *TabSigma=NULL, fltarray *CubeSigma=NULL);
	void wiener_single(fltarray &Band, int BlockSize, int s3, float noise_lvl, int LocalBS, 
		bool UseCubeSigma=false, fltarray *TabSigma=NULL, fltarray *CubeSigma=NULL);
	int get_BC_pos(int kx,int ky,int bx,int by,int LocalBS,int N);
//	void wiener_single(fltarray &Band, float noise_lvl, int LocalBS, bool UseCubeSigma=false);
//	void fdr_single(fltarray &Band, float Alpha, float SigmaNoise);

	int s2d_x0(int s);
	int s2d_y0(int s);
	int s2d_xn(int s);
	int s2d_yn(int s);
	// acessors
	void set_Verbose (Bool Verbose) {_Verbose=Verbose;}
	void set_WriteFile (Bool WriteFile) {_WriteFile=WriteFile;}
	void set_StatInfo (Bool StatInfo) {_StatInfo=StatInfo;}
	void set_TakeOnlyCompletePlane (Bool Flag) {_TakeOnlyCompletePlane=Flag;}
	void set_DoNotCorrectSigma (Bool Flag) {_DoNotCorrectSigma=Flag;}
	void set_BlockSize (int BlockSize) {_BlockSize=BlockSize;}
	void set_BlockOverlap (Bool BlockOverlap) {_BlockOverlap=BlockOverlap;}
	void set_OverlapFactor (float OverlapFactor) {_OverlapFactor=OverlapFactor;}
	void set_LinTransf (type_linelet3d_WTtrans LinTransf,SubBand1D* SB1D=NULL) {_LinTransf=LinTransf; _ptrSB1D=SB1D;}
	void set_removeplane(Bool RemovePlane) {PRad3D.setRemovePlane(RemovePlane);}
	inline int get_SizeX () {return _Nx;}
	inline int get_SizeY () {return _Ny;}
	inline int get_SizeZ () {return _Nz;}
	inline int get_BlockSize () {return _BlockSize;}
	int num_block(int Bx, int By, int  Bz) { return (Bz*(_Nxb*_Nyb)+By*_Nxb+Bx);}	
	inline int get_NbBlock () {return _NbrBlock;};
	inline int get_NbBlock_nx () {return _Nxb;};
	inline int get_NbBlock_ny () {return _Nyb;};
	inline int get_NbBlock_nz () {return _Nzb;};
	inline int get_NbScale () {return _NbrScale;}
	inline int get_RadNbPlane () {return _RadNbPlane;}
	inline int get_LinNbPlane () {return _LinNbPlane;}
	inline int nbr_band () {return NbrBand;}
	inline int size_band_nx (int b) {return TabBandNx[b];}
	inline int size_band_ny (int b) {return TabBandNy[b];}
	inline int size_band_nz (int b) {return TabBandNz[b];}
	inline int pos_band_nx (int b) {return TabDepx[b];}
	inline int pos_band_ny (int b) {return TabDepy[b];}
	inline int pos_band_nz (int b) {return TabDepz[b];}
	void get_band(fltarray& Transf, fltarray& Band, int NumBand);
	inline void set_no_recons(bool _no_recons) { no_recons=_no_recons; }
	inline void set_keep_energy(bool nc) { keep_energy=nc; }
};

/*****************************************************************************/

#endif
