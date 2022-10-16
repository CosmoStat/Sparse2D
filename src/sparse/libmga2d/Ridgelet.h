/******************************************************************************
**                   Copyright (C) 2000 by CEA + Santford University
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  Ridgelet.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the Ridgelet transform
**    ----------- 
******************************************************************************/

#ifndef _RIDGELET_H_
#define _RIDGELET_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "PrimeNumber.h"

#ifndef CURALONE
#include "MR_Obj.h"
#include "IM_Sigma.h"
#else
#define MAX_SCALE 10
#define DEFAULT_N_SIGMA 3
#endif

#define DEF_RID_NBR_SCALE -1
#define DEF_RID_MIN_BLOCK_SIZE 7

#define NBR_RID_TRANS 8
#define NBR_RID_DECIMATED 7
 
enum type_ridgelet_WTclass {RID_CL_UNKNOW, RID_CL_PYR,RID_CL_PAVE,RID_CL_ORTHO,RID_CL_FINITE};

enum type_ridgelet_WTtrans {RID_UNKNOWN,
                           RID_OWT,       // standard bi-orthogonal WT
                           RID_PYR_FFT,   // Pyramidal FFT transform
                           RID_PYR_DIRECT, // Pyramidal WT in direct space
                           RID_FINITE,     // using finite radon transform
                                          // (it is not really a ridget transform
                                          // and it produces poor results)
                            RID_FSS,       // Ridgelet transform using the FAST SLANT STACK transform + wt1d in Fourier space		  
                            RID_FSS_ORTHO, // Ridgelet transform using the FAST SLANT STACK transform + orthogonal wt
                           RID_FSS_PYR,   // Ridgelet transform using the FAST SLANT STACK transform + pyramidal WT in diract space
                           RID_UWT,       // undecimated starlet WT
                           RID_PAVE_FFT   // Undecimated FFT transform  (not implemented
                           };

#define DEF_RID_TRANS  RID_PYR_FFT
const char * StringRidTransform (type_ridgelet_WTtrans  type);
type_ridgelet_WTclass rid_class(type_ridgelet_WTtrans Type);

#define RID_MAX_SCALE 10

/***********************************************************************/

const Bool UseFSS (type_ridgelet_WTtrans);


inline void ridgelet_check_size_ima(int Nl, int Nc, int & BlockSize)
{
    if ((is_power_of_2(Nl) == False) || (is_power_of_2(Nc) == False))
    {
           cout << "Error: image size must be power of two ... " << endl;
           exit(-1);
    }
    if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False))
    {
           cout << "Error: Block size must be a power of two ... " << endl;
           exit(-1);
    }
    if (BlockSize > Nc)
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }
}

/***********************************************************************/

inline float enhance_contrast_function(float x, float Contrast_P_Param, 
                                 float Contrast_Q_Param,
				 float Contrast_M_Param, float Contrast_C_Param)
{
       float ValRet=0.;
       if (x > FLOAT_EPSILON)
       {
          if (ABS(x) < Contrast_C_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/Contrast_C_Param), (double) Contrast_P_Param);
          else if (ABS(x) < Contrast_M_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/ABS(x)), (double) Contrast_P_Param);
          else ValRet = 1.;
       }
       return ValRet;
}

/***********************************************************************/

class Ridgelet {
   SubBand1D *Ptr_SB1D; // pointer to the filter bank to use 
                        // in case of othogonal 1D WT
                        
   UndecSubBandFilter *Ptr_Undec_SB1D; // pointer to the filter bank to use in case of UWT
   
   void variance_stab(Ifloat &Transf);  
                        // variance stabilization in case of Poisson noise

   void inv_variance_stab(Ifloat &Transf); 
                       // inverse variance stabilization



   void block_transform(Ifloat &Image, Ifloat &Transf);
                      // Ridgelet transform of image block per block

   void block_recons(Ifloat &Trans, Ifloat &Image);
                     // Image reconstruction from its ridgelet transform
                     // Block per block



   // Parameter for ridgelet block analysis
   float get_weight(int Bi, int Bj, int k, int l);
       // return the weight value for pixel position k,l in the block Bi,Bj

   type_ridgelet_WTclass RidClass; // WT can be non redundant, pyramidal,
                                   // or pave
   Bool WTinFFT;      // True if the WT is calculated in the Fourier Domain

   int NlIma, NcIma;         // Input image size 
   int NlImaPow2, NcImaPow2;  // power 2 size which contains the image

   int BlockImaSize;  // Image block size whithout the overlapping size 
   int Nlb,Ncb;       // Number of blocks in both directions
   int BS;            // Image block size including the overlapping size

   int NlRid,NcRid;      // Ridgelet transform image size
   int NlBlockTransSize; // number of lines of  a transformed block  
   int NcBlockTransSize; // number of columns of a transformed block
   int BlockImaOverLap; // if BlockOverlap == True then BlockImaOverLap
                        // is initialized to BlockImaSize/2
                        // else  BlockImaOverLap == 0 

   int RidPixNbr;     // Number of pixels in the 1D WT transform
                      // for a given block 
                      // For an orthogonal 1DWT,   RidPixNbr=NcRadon
   dblarray TabAngle; // Angle relative to the different directions.
                      // TabAngle(i) = angle of the ith direction
		      //               i = 0..2*BS-1
   intarray TabPos;   // starting pixel column index table
                      // inside a given block 
   intarray TabSize;  // size band table in column direction
                      // inside a given block 
   fltarray TabNorm;  // Normalization band table to apply for keeping
                      // a normalized noise level in the bands
   int set_tab_pos(int NbrScale, int N);
                      // initialize the previous table


   void filtering_one_block(Ifloat &Imag, float NoiseIma, float N_Sigma, fltarray *TabN=NULL);
                      // Image filtering by using the ridgelet transform

   void undec_wt_filtering_one_block(Ifloat &Image, float NoiseIma, float N_Sigma);
                      // Image filtering by the ridgelet transform
                      // using an undecimated orthogonal wavelet transform
                      // if OrthoWT == True or an undecimated FFT based
                      // WT if  OrthoWT == False
                      // Results seem not better than with an OWT

   void radon_filtering_one_block(Ifloat &Imag, float NoiseIma, float N_Sigma);
                      // Image filtering by using the radon transform
                      // input image must be zero mean (i.e. wavelet scale)

   void ksig_threshold(Ifloat &Trans, float SigmaNoise, float NSigma, 
                       fltarray *TabN=NULL);
                      // Ksigma Thresholding of the ridgelet transform
                      // Transf = one block of the ridgelet transform
		      // TabN = normalization coefficient of the ridgelet
		      //        coefficients

   void iter_recons(Ifloat & RidIma, Ifloat &Data, float Lambda);
                     // Reconstruction with smoothness constraint

   Bool SetTable;     // If False, table is not filled when filtering
                      // It implies that the table must have been filled
		      // before using the routine set_tab_norm_from_simu
   
   void set_tab_angle(); // Initialize the TabAngle table.
   
   void mr_io_fill_header(fitsfile *fptr);

   void radon_transform_one_cst_block();
        // Apply the Radon transform to a constant block.
	// Routine used for the normalization of the Fast Slant Stack Radon Trans.
	// and initialize variable ImaCstRadonTransOneBlock
   Ifloat ImaCstRadonTransOneBlock; 
        // Radon transform of a cst block. Used with the  Fast Slant Stack Radon Trans.
   
   Bool AllocatedClass; // True if the class has been allocated
   public:       
    
   void transform_one_block(Ifloat &Image, Ifloat &Transf);
                       // Ridgelet transform of an image whith NbrScale  

   void recons_one_block(Ifloat &Trans, Ifloat &Image);
                       // Image reconstruction from the ridgelet transform
		       		         
   void get_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock);
                     // Extract a block from the image

   void put_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock);
                     // Put a block to the image

   void add_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock);
                     // Add a block to the image

   void get_block_trans(int Bi, int Bj, Ifloat & Transf, Ifloat &  ImaTransBlock);
                    // Get a block from the ridgelet transformed image 
                    // Transf and ImaTransBlock must be already allocated
   void put_block_trans(int Bi, int Bj, Ifloat & Transf, Ifloat & ImaTransBlock);
                    // Put a block from the ridgelet transformed image
                    // Transf and ImaTransBlock must be already allocated


   Radon RADON;        // Radon class used by the Ridgelet transform

   float get_coef(Ifloat &Trans, int i, int j, int s, int Bi, int Bj);
                      // get a ridgelet coefficient in the block Bi,Bj
                      // at scale s, and at position i,j

   void put_coef(Ifloat &Trans, int i, int j, float Val, int s, int Bi, int Bj);
                       // set a ridgelet coefficient in the block Bi,Bj
                      // at scale s, and at position i,j to the value Val
   int block_size() {return BS;} 
                        // return the block size
			
   int size_scale_nl(int s) {return NlRid;} 
                        // return the number of lines of the scale s

   int size_scale_nc(int s) {return TabSize(s)*Ncb;}
                        // return the number of columns of the scale s

   int ipos(int s) { return 0;};
                       // return the starting position (line number) 
                       // of the scale s (which contains all blocks)
                       // in the transform

   int jpos(int s) { return (TabPos(s)*Ncb);}
                       // return the starting position (column number) 
                       // of the scale s (which contains all blocks) 
                       // in the transform

   int ipos(int s, int Bi) { return NlBlockTransSize*Bi;};
                       // return the starting position (line number) 
                       // of the scale s of the block Bi

   int jpos(int s, int Bj) { return (TabPos(s)*Ncb + Bj*TabSize(s));}
                       // return the starting position (column number) 
                       // of the scale s of the block Bj

   void get_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s);
                        // extract the scale s from the ridgelet 
                        // transform ImaTrans. The return image contains
                        // the band s of all the blocks
   
   void get_scale_direction(Ifloat &ImaTrans, Ifloat &ImaScale, 
                            int s, float Angle);
                        // extract the scale s from the ridgelet 
                        // transform ImaTrans. The return image contains
                        // the band s of all the blocks in only
			// one direction by Angle
  void get_scale_direction(Ifloat &ImaTrans, Ifloat &ImaScale, int s, 
                          int NumDirect);
		        // Ditto, but take in input the angle number NumDirect
			// instead of the angle value
			// NumDirect is inside [0..2*BS-1[
			
   void put_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s);
                        // insert a scale into the ridgelet transform
                        // ImaScale must contains
                        // the band s of all the blocks

   int rid_nl () { return NlRid;}  // return the number of lines of the
                                   // ridgelet transform

   int rid_nc () { return NcRid;}  // return the number of columns 
                                   // of the ridgelet transform

   int rid_block_nl () { return Nlb;}  // return the number of blocks 
                                       // in the line direction

   int rid_block_nc () { return Ncb;}  // return the number of blocks 
                                       // in the column direction

   int rid_block_nbr () { return Nlb*Ncb;}  // return the number of blocks

   int rid_np () { return  RidPixNbr;}  // return the number of pixels
                                        // in one line of the ridgelet 
                                        // transform for one block
   
   int rid_one_block_nl () { return NlBlockTransSize;}  
                                        // return the number of lines
                                        // of one ridgelet block 
   int rid_one_block_nc () { return NcBlockTransSize;}  
                                        // return the number of columns
                                        // of one ridgelet block 
										
   int rid_pos (int s) { return TabPos(s);}
                                       // Return the starting pixel column index 
                                       // of the scale s (for a given one block)
                                       // (s=0,..,NbrScale-1)
   
   int rid_one_block_nc (int s) { return TabSize (s);}  
   int rid_size (int s) { return TabSize (s);}
                                       // Return the number of pixels
                                       // of the scale s in a given line 
                                       // (s=0,..,NbrScale-1)
                                       // (for a given one block)
 
   float rid_norm (int s) { return TabNorm (s);}
                                       // Normalization parameter
                                       // As the noise variance is not kept
                                       // rid_norm return the weight to 
                                       // apply to the noise standard deviation
                                       // in order to compare a coefficient
                                       // to its noise standard deviation

   int get_num_direction(float Angle, Bool InDegree=False);
                        // Get the direction number relative to closest angle
			// to a given angle.
			// Return value is inside [0,2*BS[.
			// if InDegree=True, the angle is given in degrees.
			
   float get_angle_direction(int NumDirect, Bool InDegree=False);
			// Return the angle relative to a given direction
			// if InDegree=True, the angle is return in degrees.
			
   Bool GetAutoNbScale; // If true the number of scale is automatically
                        // calculated by:
                        // NbrScale = fix( log( (3N/4) / log(2)));

   int NbrScale;       // Number of scales used by the 1D wavelet transform

   Bool OnlyHighFreq; // If true, the filtering is applied only on the
                      // high frequency part of the image
                      // this allows us to reduce artefacts relative to
                      // block windowing

   type_ridgelet_WTtrans RidTrans; // Ridgelet transform type

   Bool BlockOverlap;  // If True, Overlapped blocks are used.

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

   Bool WPTrans;      // Wavelet packets are applied on the first scale
                      // of the ridgelet transform
   Bool WeightBefTrans;		       
   Bool MadThreshold;  // when true, the threshold level is estimating from 
                       // the absolute median value instead of the ksigma
                       // thresholding
                       // (Default is no, and it produce poor results!)
                       // Thus parameter is used by the filtering routine
   Bool NoiseInvariantPerRotation;
                       // If True, the noise is considered as invariant
		       // per orientation.
		       
                       // A smoothness constraint can be added to the 
                       // reconstructed image. The three following parameters
                       // are used for this iterative reconstruction
   int NbrFilterIter;  // Number of iteration for the reconstruction
                       // Default is 1, no iterative reconstruction
   float RegulParam;   // Regularisation parameter. Default is 0.2
   float CvgParam;     // Convergence parameter for the reconstruction
                       // Default is 1.
   float KillFluctu;
   type_format FormatInputImag;   // data format of the input image
   Bool VarNorm;       // If set, the ridgelet coefficient are divided by sqrt(BlockSize)
                       // By default, no.
   Ifloat *TabImaSigmaStat; // if a normalization per block is required
                            // the normalized image is store in TabImaSigmaStat
                            // TabImaSigmaStat[s](i,j) corresponds to the
			    // normalization value at scale s, angle i
			    // and position j inside a given block.
			    
   void reset(SubBand1D *SB1D=NULL);
                      // Reset all internal variables

   Ridgelet(SubBand1D &SB1D)  {reset(&SB1D);}  
   Ridgelet()  {reset();}  

   void alloc(int Nlima, int Ncima, int BlockSize);
                     // initialize the variable for block transform from
                     // the input image size
                     // this routine is automatically called from transform
                     // recons, and filtering routines

   void set_block_param_from_transform(int Nlt, int Nct, int BlockSize);
                     // initialize the variable for block reconstruction from
                     // the input ridgelet transform size
                     // this routine is automatically called from  
                     // recons routine 

   void transform(Ifloat &Image, Ifloat &Transf, int BlockSize);
                      // Ridgelet transform of an image whith NbrScale
                      // aplied block per block
                      // if BlockSize=0, BlockSize=image size

   void transform(Ifloat &Image, Ifloat &Transf);
                      // Ridgelet transform of an image whith NbrScale
                      // aplied block per block
		      // The class must have been initalized before.
                      // by set_block_param or 
		      // set_block_param_from_transform.
		       
   void recons(Ifloat &Trans, Ifloat &Image, int BlockSize);
                       // Image reconstruction from the ridgelet transform
                       // applied block per block

   void recons(Ifloat &Trans, Ifloat &Image);
                       // Image reconstruction from the ridgelet transform
                       // applied block per block.
		       // The class must have been initalized before.
                       // by set_block_param or 
		       // set_block_param_from_transform.
		       		       
   void recons(Ifloat &Transf, Ifloat &Image, int BlockSize, 
                       int InputImaNl, int InputImaNc);
                       // Image reconstruction from the ridgelet transform
                       // applied block per block
                       // InputImaNl, InputImaNc are the dimensions of 
                       // the reconstructed image

   void set_tab_norm(int NbrScale);
                      // initialize the normalization table

   void set_tab_norm_with_tab(int NbrScale, fltarray & Tab);
                     // initialize the band normalization table with Tab
                     // This routine is used by the curvelet transform
                     // because the noise inside a given band is colorated

   void set_tab_norm_from_simu(Ifloat &Simu, int NbrScale);
                     // initialize the band normalization table with
		     // simulated data. 

   void filtering(Ifloat &Ima, float NoiseIma, float N_Sigma, int BlockSize=0, fltarray *TabN=NULL);
                      // Image filtering by using the ridgelet transform
                      // block per block
                      // BlockSize = 0 ==> only one block, i.e. the image
 
   void undec_wt_filtering(Ifloat &Image, float NoiseIma, float N_Sigma, int BlockSize=0);
                      // Image filtering by the ridgelet transform
                      // using an undecimated orthogonal wavelet transform
                      // if OrthoWT == True or an undecimated FFT based
                      // WT if  OrthoWT == False
                      // BlockSize = 0 ==> only one block, i.e. the image

   void info_stat(Ifloat &Transf);
                     // print some statistical information, band per band
                     // of the ridgelet transform

   void get_low_and_freq(Ifloat &Image,  Ifloat &Low, Ifloat &High, int BlockSize);
                     // decompose Image into two frequency bands: Low and High
                     //  Image = Low + Hig 

   void radon_filtering(Ifloat &Ima, float NoiseIma, float N_Sigma, int BlockSize=0);
                      // Image filtering by using the radon transform
                      // block per block 
                      // input image must be zero mean (i.e. wavelet scale)
                      // BlockSize = 0 ==> one block, i.e. the image

   void thresholding(Ifloat &Transf, float NoiseIma, float N_Sigma);
                     // threshold a ridgelet transform 
                     // by a ksigma thresholding

   void mad_threshold(Ifloat &Transf, float N_Sigma);
                      // MAD thresholding of the ridgelet transform
                      // Level = Nsigma*get_abs_median(Band) / 0.6745;
                      // Transf = one block of the ridgelet transform
		      
   void kill_fluctu(Ifloat &Ima, float NoiseIma, float N_Sigma, 
                           int BlockSize);
                    // if the input image to filter is zero mean, a post
                    // processing can be applied consisting in thresholding
                    // small values under:
                    //    T = NoiseIma*N_Sigma / sqrt((float) B) / 2.

   void init_scale(Ifloat &Transf, int s);
                    // set to zero the scale s of the transform

   void normalize_coef(Ifloat &Transf, Bool Inverse=False, float NoiseIma=1.);
                    // Normalize all coefficients so that noise variance is 
                    // equal to 1.
		    
   void get_imablocksigma(Ifloat &ImaTrans, int s, Ifloat &ImaSigmaStat, 
                          Bool UseSigma=False);
                    // Calculate the MAD for 
		    // each ridgelet pixel at a given direction and
		    // at a give position
		    // ImaSigmaStat(i,j) = OUT: MAD of the pixels values 
		    //        at position (i,j) in the scale s.
		    //        i = 0..2*BS-1
		    //        j = 0:BS-1
		    // If UseSigma == True, use the standard deviation instead
		    // of the MAD estimation.
   
   void get_blockvect(Ifloat &ImaTrans, int s, 
                      int Posi, int Posj, fltarray & VectPix);
		    // return a vector contains all pixels at a position
		    // Posi,Posj in the blocks at scale s		    
                    // VectPix(i) = OUT: pixel value 
		    //        at position (Posi,Posj) in the scale s.
		    //        in the ith block.
		    //        i = 0..rid_block_nbr()-1

   void get_anglevect(Ifloat &ImaTrans, int s, int Posi, fltarray & VectPix);
   		    // return a vector contains all pixels at a given
		    // angle (angle number Posi) in the blocks at scale s		    
                    // VectPix(i) = OUT: pixel value 
		    //        at position (Posi, Posj) in the scale s, 
		    //        where Posj = corresponds to the block
		    //        number i%NcOneBlock at pixel position 
		    //        i mod NcOneBlock 
		    //        Posi = 0..2*BS-1
   void put_anglevect(Ifloat &ImaTrans, int s, int Indi,  fltarray & VectPix);
   void mad_normalize_per_angle(Ifloat &ImaTrans);
   void init_tab_sigma_block(Ifloat &ImaTrans, Bool UseSigma=False);
   void normalize_block(Ifloat &ImaTrans, Bool InverseNormalization=False);
   void mad_normalize_per_block(Ifloat &ImaTrans, Bool UseSigma=False);
   void get_scale_without_bord(Ifloat &ImaTrans, int s, Ifloat &ImaScale);
            
    void read  (char *Name, Ifloat &Ima);
    // Read a ridgelet transform and initialize the class
    
    void write (char *Name, Ifloat &Ima);
    // write the ridgelet transform in the FITS format
    		     
   ~Ridgelet(){ delete [] TabImaSigmaStat; reset();}
};

/***********************************************************************/

#endif
