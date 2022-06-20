/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  14/03/00 
**    
**    File:  MC_Comp.h
**
*******************************************************************************
**
**    DESCRIPTION  multichannel image compression
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef __MR_COMP3D__
#define __MR_COMP3D__

// #include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_BIO.h"
#include "SB_Filter.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "IM3D_IO.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "MR_Comp.h"
#include "IM_Lut.h"

class BlockMRInfo;
                             
class MR_Comp3DData {
    short ValIdent[3];    
    Bool IntTrans;          // if set, the compression works with integer
    Bool NoiseThresholding; // if set, supress the noise
    void header_resol();    // modify the fits header taking into account
                            // the finale resolution of the decompressed image


    void encode_residu();     // encore the residual part
    void store_noise_info();  // store nosie information
    void encode_data();       // encode the data
    void variance_stab();     // apply variance stabilization
    void replace_badpix();    // replace bad pixels
    void read_data();         // read the data
    void store_header();      // store the header
    void print_info();        // print some information
    
    void set_tabblock();      // intialize table for block operations
    void get_block(int NumBlock);  // get a block 
    void put_block(int NumBlock);  // put a block
    void decode_residu();          // decode the residual part
    void read_noise_info();        // read the noise information
    void read_tinfo();             // read the general information 
    void decode_data();            // decode the data
    void read_header();            // read the header
 
    void seek_band(int NbrBand);   // seek NbrBand bands
    void seek_resol(int NbrSkipResol); // seek NbrSkipResol resolutions
    void seek_block(int Numblock);     // seek to the block number Numblock
    // void dec_data();                   // last operation before saving the  results
    void set_block_info_data();   // set information when block decoding
    void upblock(Ifloat &Low);    // update the current block
 

    public:
    // image size
    int Nl,Nc,Nf;

    // use for compression with float
    fltarray Dat;
    fltarray Residual;
    fltarray *BDat;

    // use for compression with integer
    int *IResidual;
    int *IDat;
    int *BIDat;
    
    FILE *InFile;
    FILE *OutFile;
    int ReadHd;
    fitsstruct Header;
    float MinDat;
    float MaxDat;
    float SigmaResi;
   
    // parameter for compression
    char *Cmd;
    char *File_Name_Imag;      /* input file image */
    char *File_Name_Transform; /* output file name */
    int Nbr_Plan;                  /* number of scales */
    float N_Sigma;                 /* number of sigma (for the noise) */
    float Noise_Ima;              /* noise standard deviation */
    float *TabSigmaNoise;         /* noise standard deviation (one per frame */

    type_noise Stat_Noise;        /* type of noise */
    type_transform Transform;     /* type of transform */
    Bool SupIsol;                 /* suppress isolated pixel in the support */
    int MaxIter;
    float NoiseQuantif;           /* Quantization parameter for the estimated
                                     noise map */
    float SignalQuantif;         /* Quantization parameter for the detected
                                    signal */
    int MedianWindowSize;        /* Median window size for the transform */
    Bool NoBscale;               /* If input image format = FITS
                                    and optimed method, then if NoBscale
                                    is set, the fits operation:
                                    Y = X * Bscale + Bzero is not applied */
    Bool SupNeg;                 /* Suppress negative detection */
    int KeepResi;
           // if KeepResi=1, the noise is compressed with a lost of information
           // if KeepResi=2, the noise is keeped without loosing anything
    int KeepFitsHeader;
    type_comp Comp_Method;  // Compression method
    int Resol;              // Resolution level for the decompression
    int IterRec;            // number of iterations for the decompression
    Bool KillDiag0;         // if true, supress Diagonal scale 0
    char OutputType;        // output format for the decompression
    Bool SimuNoise;         // add a simlated noise map at the decompression
    int Elstr_Size;        //  Size Structural element
    int Elstr_Shape;       //  Square=0; Cercle=1
    int N_Iter;            //  No iteration for morphomate
    Bool UseMR_for_BGR;    //  if true, multiresolution is used for
                           // background estimation
    int NpixBgr;           // background image resolution (in pixel)
                           //  (used if the bgr is calculated from
                           //   multiresolution 
    int NbrImag;           // Number of images in the file
    Bool RGBImage;         // True if it is a RGB image
    Bool RGBChromDownSample; // If true, chrominance map are downsampled by 2.
    Bool UseBlock;         // if True, image are compressed by block
    int  NbrBlock;         // Number of blocks
    int  BlockSize;        // Block size in both directions
    int  BlockOverlap;      // Overlap between two blocks
    int  L_Nl, L_Nc;       // size of the big image
    int NbrBlocNl, NbrBlocNc; // Number of blocks in both directions
    int *TabBlockNl, *TabBlockNc;
    int *TabBlockDepNl, *TabBlockDepNc;
    
    type_3d_format FormatInput; // Input data format;
    type_quant TypeQuant ;  // Type of quantification
    type_coding TypeCoder;  // Coding type
    Bool UseBudget;         // if set, a budget is fixed
    float CompressionRatio; // wanted compression ratio
    unsigned int TargetNbrByte; // wanted number of bytes
    unsigned int EncodedNbrByte; // achieved number of bytes
    Bool NoiseInData;       // Data contains noise 
    Bool ReadInPseudo;
    Bool Verbose;
    type_lift LiftingTrans; // Lifting transform parameter 
                            // (if the lifing transform is used).
   MR_Comp3DData() 
   {    
    RGBChromDownSample = True;
    ReadInPseudo=False;
    LiftingTrans = DEF_LIFT;
    Verbose = False;
    BlockOverlap=0;
    TabBlockNl=NULL;
    TabBlockNc=NULL;
    TabBlockDepNl=NULL;
    TabBlockDepNc=NULL;
    NbrBlocNl=0;
    NbrBlocNc=0;
    L_Nl = 0;
    L_Nc = 0;
    TargetNbrByte = 0;
    EncodedNbrByte = 0;
    NoiseInData = True;
    UseBudget = False;
    CompressionRatio=-1;
    IntTrans = False;
    NoiseThresholding = True;
    TypeQuant = Q_UNIFORM ;
    TypeCoder = C_HUFFMAN;
    FormatInput = F3D_UNKNOWN;
    BlockSize = -1;
    UseBlock = False;
    NbrBlock = 1;
    NbrImag = 1;
    KillDiag0=False;
    IterRec=0;
    NpixBgr = DEFAULT_NPIX_BGR;
    UseMR_for_BGR = False;
    NoBscale = False;
    SupNeg = False;
    MedianWindowSize=DEF_IPYRMED_WINDOW_SIZE;
    Cmd = NULL;
    File_Name_Imag = NULL;
    File_Name_Transform = NULL;
    Nbr_Plan=DEFAULT_CMP_NBR_SCALE;
    N_Sigma=DEFAULT_N_SIGMA;
    Noise_Ima=0.;
    Stat_Noise = DEFAULT_STAT_NOISE;
    Transform = TM_PYR_MEDIAN;
    SupIsol=False;
    MaxIter = DEFAULT_MAX_ITER_COMP;
    NoiseQuantif = DEFAULT_NOISE_QUANTIF;
    SignalQuantif = DEFAULT_SIGNAL_QUANTIF;
    KeepResi = 0;
    KeepFitsHeader = 0;
    Comp_Method = COMP_MRMEDIAN;
    MinDat=0.;
    MaxDat=0.;
    ValIdent[0]= 377;
    ValIdent[1]= 331;
    ValIdent[2]= 377;
    IResidual = NULL;
    IDat = NULL;
    BIDat = NULL;
    BDat = NULL;
    Resol = DEFAULT_RESOL; 
    SimuNoise = False;
    OutputType='u';
    Elstr_Size=DEFAULT_ELSTR_SIZE;
    Elstr_Shape=DEFAULT_ELSTR_SHAPE;
    N_Iter=0;
   }
   
   // Decompression initialization
   void init_decompress(); // open the input file and read the header

   // apply the compression
   void compress();
   
   // apply the decompression ==> the field Dat contains the decompressed image
   // the field Header contains the fits header of Dat
   void decompress();
   
   // deompress one block ==> the field Dat contains the decompressed image
   // the field Header contains the fits header of Dat
   void decompress(int NumBlock);

   // get a resolution: if the image has been compressed by block, it
   // returns the resolution of one single block
   void get_resol(int Resol, int &SizeBufResol, char * &BuffResol, int NumBlock=1);

   // increase the resolution of an image: if the image has been compressed by 
   // block, it increases the resolution of one single block
   void up_mrresol(fltarray &ImagLowResol, int NumBlock=1);
   void up_mrresol(fltarray &ImagLowResol, BlockMRInfo &BV, BlockMRInfo &BVup);

   // write the result to the disk
   void write();
   
   ~MR_Comp3DData();
};

 
#endif
