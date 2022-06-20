/******************************************************************************
**                   Copyright (C) 2001 by CEA 
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  25/08/2001 
**    
**    File:  Ridgelet3D3D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the 3D Ridgelet3D transform
**    ----------- 
******************************************************************************/

#ifndef _RIDGELET3D_H_
#define _RIDGELET3D_H_

#ifdef WINDOWS
	#define USE_OMP_R 0
#else
	#define USE_OMP_R 1
#endif

#if USE_OMP_R
	#include <omp.h>
#endif

#include "IM_Obj.h"
#include "IM3D_Radon.h"
#include "IM3D_IO.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM3D_Block.h"

#ifndef CURALONE
#include "MR_Obj.h"
#include "IM_Sigma.h"
#else
#define MAX_SCALE 10
#define DEFAULT_N_SIGMA 3
#endif

#define DEF_RID3D_NBR_SCALE -1
#define DEF_RID3D_MIN_BLOCK_SIZE 7

enum type_ridgelet3d_WTclass {RID3D_CL_UNKNOW, RID3D_CL_PYR,RID3D_CL_ORTHO};

#define NBR_RID3D_TRANS 2
enum type_ridgelet3d_WTtrans {RID3D_UNKNOWN,
                              RID3D_OWT,       // standard bi-orthogonal WT
                              RID3D_PYR_FFT   // Pyramidal FFT transform
                             };
#define DEF_RID3D_TRANS  RID3D_PYR_FFT

const char * StringRid3DTransform (type_ridgelet3d_WTtrans  type);
type_ridgelet3d_WTclass rid3d_class(type_ridgelet3d_WTtrans Type);

#define RID3D_MAX_SCALE 10

/***********************************************************************/

inline void ridgelet3d_check_size_cube(int Nx, int Ny, int Nz, int & BlockSize)
{
    if ((is_power_of_2(Nx) == False) || (is_power_of_2(Ny) == False)
                                     || (is_power_of_2(Nz) == False))
    {
           cout << "Error: cube size must be power of two ... " << endl;
           exit(-1);
    }
    if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False))
    {
           cout << "Error: Block size must be a power of two ... " << endl;
           exit(-1);
    }
    if (BlockSize > Ny)
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }
}


/***********************************************************************/

class Ridgelet3D {
   SubBand1D *Ptr_SB1D; // pointer to the filter bank to use 
                        // in case of othogonal 1D WT

   void variance_stab(Ifloat &Transf);  
                        // variance stabilization in case of Poisson noise

   void inv_variance_stab(Ifloat &Transf); 
                       // inverse variance stabilization

   void transform_one_block(fltarray &Cube, Ifloat &Transf);
                       // Ridgelet3D transform of a cube whith NbrScale

   void block_transform(fltarray &Cube, Ifloat &Transf);
                      // Ridgelet3D transform of a cube,  block per block

   void block_recons(Ifloat &Transf, fltarray &Cube);
                  // Image reconstruction from its ridgelet transform
                     // Block per block

   void recons_one_block(Ifloat &Trans, fltarray &Cube);
                       // Cube reconstruction from the ridgelet transform

   // Parameter for ridgelet block analysis
   float get_weight(int Bi, int Bj, int Bk, int k, int l, int m);
       // return the weight value for pixel position k,l in the block Bi,Bj

   type_ridgelet3d_WTclass RidClass; // WT can be non redundant, pyramidal,
                                   // or pave
   Bool WTinFFT;      // True if the WT is calculated in the Fourier Domain

   // int Nx,Ny,Nz;         // Input image size 
   // int NxPow2, NyPow2, NzPow2;  // power 2 size which contains the image

   int UserSizeBlockParam; // User block size parameter
   // int BlockCubeSize;  // Cube block size whithout the overlapping size 
   // int Nxb,Nyb,Nzb;       // Number of blocks in both directions
   // int NbrBlock;
   // int BS;            // Cube block size including the overlapping size
   Block3D B3D;       // Class for 3D block management
   
   int NlRid,NcRid;      // Ridgelet3D transform image size
   int NlBlockTransSize; // number of lines of  a transformed block  
   int NcBlockTransSize; // number of columns of a transformed block
    
   int BlockCubeOverLap; // if BlockOverlap == True then BlockCubeOverLap
                        // is initialized to BlockImaSize/2
                        // else  BlockCubeOverLap == 0 

   int RidPixNbr;     // Number of pixels in the 1D WT transform
                      // for a given block 
                      // For an orthogonal 1DWT,   RidPixNbr=NcRadon
   intarray TabPos;   // starting pixel column index table
                      // inside a given block 
   intarray TabSize;  // size band table in column direction
                      // inside a given block 
   fltarray TabNorm;  // Normalization band table to apply for keeping
                      // a normalized noise level in the bands
   int set_tab_pos(int NbrScale, int N);
                      // initialize the previous table

   void  get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
                       {B3D.get_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
                      // Extract a block from the cube

  void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
                      {B3D.put_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}
                      // Put a block to the cube

   void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
                      {B3D.add_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
		     // Add a block to a cubeCube

   void get_block_trans(int Bi, int Bj, int Bk, Ifloat & Transf, Ifloat &  ImaTransBlock);
                    // Get a block from the ridgelet transformed image 
                    // Transf and ImaTransBlock must be already allocated
   void put_block_trans(int Bi, int Bj, int Bk, Ifloat & Transf, Ifloat & ImaTransBlock);
                    // Put a block from the ridgelet transformed image
                    // Transf and ImaTransBlock must be already allocated

   void filtering_one_block(fltarray &Imag, float NoiseIma, float N_Sigma, fltarray *TabN=NULL);
                      // Image filtering by using the ridgelet transform

   void ksig_threshold(Ifloat &Trans, float SigmaNoise, float NSigma);
                      // Ksigma Thresholding of the ridgelet transform
                      // Transf = one block of the ridgelet transform

   Bool SetTable;     // If False, table is not filled when filtering
                      // It implies that the table must have been filled
		      // before using the routine set_tab_norm_from_simu

   void mr_io_fill_header(fitsfile *fptr);

   public:           


   Radon3D RADON;        // Radon class used by the Ridgelet3D transform

//   float get_coef(Ifloat &Trans, int i, int j, int s, int Bi, int Bj, int Bk);
//                      // get a ridgelet coefficient in the block Bi,Bj,Bk
//                      // at scale s, and at position i,j

//   void put_coef(Ifloat &Trans, int i, int j, float Val, int s, int Bi, int Bj);
//                       // set a ridgelet coefficient in the block Bi,Bj
//                      // at scale s, and at position i,j to the value Val

   inline int nx() { return B3D.nx();}
   inline int ny() { return B3D.ny();}
   inline int nz() { return B3D.nz();}
   inline int nbr_block_nx() { return B3D.nbr_block_nx();}
   inline int nbr_block_ny() { return B3D.nbr_block_ny();}
   inline int nbr_block_nz() { return B3D.nbr_block_nz();}
   inline int nbr_block()    { return B3D.nbr_block();}
   inline int block_size()   { return B3D.block_size();}
   
   int size_scale_nl(int s) {return NlRid;} 
                        // return the number of lines of the scale s

   int size_scale_nc(int s) {return TabSize(s)*nbr_block();}
                        // return the number of columns of the scale s

   int ipos(int s) { return 0;};
                       // return the starting position (line number) 
                       // of the scale s (which contains all blocks)
                       // in the transform

   int jpos(int s) { return (rid_pos(s)*rid_block_nc());}
                       // return the starting position (column number) 
                       // of the scale s (which contains all blocks) 
                       // in the transform

   int ipos(int s, int Bi) { return 0;};
                       // return the starting position (line number) 
                       // of the scale s of the block Bi

   int jpos(int s, int Bi) { return (jpos(s) + Bi*rid_size(s));}
                       // return the starting position (column number) 
                       // of the scale s of the block Bi
   int num_block(int Bx, int By, int  Bz)
                           { return (Bz*(nbr_block_nx()*nbr_block_ny())+By*nbr_block_nx()+Bx);}
   void get_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s);
                        // extract the scale s from the ridgelet 
                        // transform ImaTrans. The return image contains
                        // the band s of all the blocks

   void put_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s);
                        // insert a scale into the ridgelet transform
                        // ImaScale must contains
                        // the band s of all the blocks

   int rid_nl () { return NlRid;}  // return the number of lines of the
                                   // ridgelet transform

   int rid_nc () { return NcRid;}  // return the number of columns 
                                   // of the ridgelet transform

   int rid_block_nl () { return 1;}  // return the number of blocks 
                                       // in the line direction

   int rid_block_nc () { return B3D.nbr_block();}  // return the number of blocks 
                                       // in the column direction

   int rid_block_nbr () { return B3D.nbr_block();}  // return the number of blocks

   int rid_np () { return  RidPixNbr;}  // return the number of pixels
                                        // in one line of the ridgelet 
                                        // transform for one block

   int rid_pos (int s) { return TabPos(s);}
                                       // Return the starting pixel column index 
                                       // of the scale s (for a given block)
                                       // (s=0,..,NbrScale-1)

   int rid_size (int s) { return TabSize (s);}
                                       // Return the number of pixels
                                       // of the scale s in a given line 
                                       // (s=0,..,NbrScale-1)
                                       // (for a given block)
 
   float rid_norm (int s) { return TabNorm (s);}
                                       // Normalization parameter
                                       // As the noise variance is not kept
                                       // rid_norm return the weight to 
                                       // apply to the noise standard deviation
                                       // in order to compare a coefficient
                                       // to its noise standard deviation

   Bool GetAutoNbScale; // If true the number of scale is automatically
                        // calculated by:
                        // NbrScale = fix( log( (3N/4) / log(2)));

   int NbrScale;       // Number of scales used by the 1D wavelet transform

   type_ridgelet3d_WTtrans RidTrans; // Ridgelet3D transform type

   Bool BlockOverlap;  // If True, Overlapped blocks are used.
   float OverlapFactor; // Overlapping factor. Between 0 and 0.5
   void set_OverlapFactor (float _OverlapFactor) {OverlapFactor=_OverlapFactor;}

   Bool Verbose;       // Verbose mode  (default is no).

   Bool StatInfo;      // If true, statistical information of each band 
                       // is printed (default is no).

   Bool VarianceStab;  // when True, the Anscombe variance stabilization
                       // is applied on the Radon coefficient 
                       // (Default is no).

   int FirstDetectScale; // First detection scale. All the finest scale
                         // are set to 0 when filtering
                         // By default, FirstDetectScale = 0
                         // ==> no scale is killed 

   Bool KillLastScale;  // Reset to 0 all coefficient of the last scale 
                        // when filtering. Default is no.

   Bool ColTrans;      // if RidTrans == RID_OWT then column can also
                       // be WT transformed. It is done when ColTrans=True

   Bool MadThreshold;  // when true, the threshold level is estimating from 
                       // the absolute median value instead of the ksigma
                       // thresholding
                       // (Default is no, and it produce poor results!)
                       // Thus parameter is used by the filtering routine

   void reset(SubBand1D *SB1D=NULL);
                      // Reset all internal variables

   Ridgelet3D(SubBand1D &SB1D)  {reset(&SB1D);}  
   Ridgelet3D()  {reset();}  

   void set_block_param(int Nx, int Ny, int Nz, int BlockSize);
                     // initialize the variables for the block transform from
                     // the input cube size
                     // this routine is automatically called from transform
                     // recons, and filtering routines

   void set_block_param_from_transform(int Nlt, int Nct, int BlockSize);
                     // initialize the variable for block reconstruction from
                     // the input ridgelet transform size
                     // this routine is automatically called from  
                     // recons routine 

   void transform(fltarray &Image, Ifloat &Transf, int BlockSize=0);
                      // Ridgelet3D transform of a cube whith NbrScale
                      // aplied block per block
                      // if BlockSize=0, BlockSize=image size

   void recons(Ifloat &Trans, fltarray &Cube, int BlockSize=0);
                       // Cube reconstruction from the ridgelet transform
                       // applied block per block

   void recons(Ifloat &Transf, fltarray &Cube, int BlockSize, 
                       int InputNx, int InputNy, int InputNz);
                       // Cube reconstruction from the ridgelet transform
                       // applied block per block
                       // InputN* are the dimensions of 
                       // the reconstructed cube

   void set_tab_norm(int NbrScale);
                      // initialize the normalization table

   void set_tab_norm_with_tab(int NbrScale, fltarray & Tab);
                     // initialize the band normalization table with Tab
                     // This routine is used by the curvelet transform
                     // because the noise inside a given band is colorated

   void set_tab_norm_from_simu(Ifloat &Simu, int NbrScale);
                     // initialize the band normalization table with
		     // simulated data. 

   void filtering(fltarray &Ima, float NoiseIma, float N_Sigma, int BlockSize=0, fltarray *TabN=NULL);
                      // cube filtering by using the ridgelet transform
                      // block per block
                      // BlockSize = 0 ==> only one block, i.e. the image
 
   void info_stat(Ifloat &Transf);
                     // print some statistical information, band per band
                     // of the ridgelet transform

   void thresholding(Ifloat &Transf, float NoiseIma, float N_Sigma);
                     // threshold a ridgelet transform 
                     // by a ksigma thresholding

   void mad_threshold(Ifloat &Transf, float N_Sigma);
                      // MAD thresholding of the ridgelet transform
                      // Level = Nsigma*get_abs_median(Band) / 0.6745;
                      // Transf = one block of the ridgelet transform
		      
   void init_scale(Ifloat &Transf, int s);
                    // set to zero the scale s of the transform

   void normalize_coef(Ifloat &Transf, Bool Inverse=False, float NoiseIma=1.);
                    // Normalize all coefficients so that noise variance is 
                    // equal to 1.
		    
    void read  (char *Name, Ifloat &Ima);
    void write (char *Name, Ifloat &Ima);

   ~Ridgelet3D(){reset();}
};

/***********************************************************************/

#endif
