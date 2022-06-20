/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Comp.h
**
*******************************************************************************
**
**    DESCRIPTION  prototype for image compression
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

#ifndef __MR_COMP__
#define __MR_COMP__

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_IO.h"
#include "IM_BIO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"
#include "IM_Lut.h"
#include "IM_Comp.h"

#define DEFAULT_MAX_ITER_COMP 0
#define DEFAULT_EPSILON_COMP 1e-5
#define DEFAULT_SIGNAL_QUANTIF 1.5
#define DEFAULT_NOISE_QUANTIF 0.5
#define DEFAULT_SUP_ISOL False
#define DEFAULT_CMP_NBR_SCALE 6
#define DEFAULT_RESOL -1  
 
// default values for morphology mathematic compression method
#define DEFAULT_ELSTR_SIZE 5   //  Size Structural element
#define DEFAULT_ELSTR_SHAPE 1  //  Square=0; Cercle=1
#define DEFAULT_NPIX_BGR 16    //  background image resolution (in pixel)
                               //  (used if the bgr is calculated from
                               //   multiresolution  
                               
                               
 // compression decompression using  pyramidal median transform                          
void mr_compress (Ifloat &Imag, int Nbr_Plan, float &Noise, float NSigma, 
            float SignalQuant, int Budget, FILE *FileDes, 
            Ifloat &Residual, Bool SupIsol=DEFAULT_SUP_ISOL,
            Bool NoiseInData=True, int MedianWindowSize=DEF_IPYRMED_WINDOW_SIZE,
            Bool SupNeg = False, Bool Verbose = False);

void mr_decompress (Ifloat &Imag,  int Resol, FILE *FileDes, int &Nli, int &Nci, Bool Verbose = False);

void mri_compress (int *Imag, int Nl, int Nc, int Nbr_Plan, float *Noise,
                    float SignalQuantif, float  NSigma,
                    FILE *FileDes, int *Residual, 
                    Bool SupIsol=DEFAULT_SUP_ISOL,
                    Bool NoiseInData=True, Bool Rec=True, 
                    int MedianWindowSize=DEF_IPYRMED_WINDOW_SIZE,
                    Bool SupNeg = False, Bool Verbose = False);

void mri_decompress (int **Imag, int *Nl, int *Nc, int *Nbr_Plan, 
                     int Resol, FILE *FileDes, int &Nli, int &Nci, 
		     Bool Verbose = False);

// compression decompression using morphology mathematic
void morphoi_compress (int *Imag, int Nl, int Nc, FILE *FileDes, 
                       float & Noise,float N_SigmaNoise, 
                       int Elstr_Size = DEFAULT_ELSTR_SIZE, 
                       int Elstr_Shape=DEFAULT_ELSTR_SHAPE,
                       Bool UseMR_for_BGR=False, int NPixBgr=DEFAULT_NPIX_BGR,
		       Bool Verbose = False);
 
void morpho_compress (Ifloat &Imag, FILE *FileDes,   
                     float &Noise_Ima, float N_Sigma, 
                     int Elstr_Size=DEFAULT_ELSTR_SIZE,
                     int Elstr_Shape=DEFAULT_ELSTR_SHAPE,
                     Bool UseMR_for_BGR=False, int NPixBgr=DEFAULT_NPIX_BGR,
		     Bool Verbose = False);
 
void morphoi_decompress (int *Imag, int & N_Iter, int Nl, int Nc,
                         FILE *FileDes, int Elstr_Size=DEFAULT_ELSTR_SIZE,
                         int Elstr_Shape=DEFAULT_ELSTR_SHAPE, 
			 Bool Verbose = False);
 
void morpho_decompress (Ifloat &Imag, int & N_Iter, FILE *FileDes, 
                        int Elstr_Size=DEFAULT_ELSTR_SIZE, 
                        int Elstr_Shape=DEFAULT_ELSTR_SHAPE,
			Bool Verbose = False);
			
void mr_i_ortho_compress (Iint &Imag, int Nbr_Plan, FILE *FileDes, Bool Verbose= False);
void mr_i_ortho_decompress (Iint &Imag, int Resol, int & Nbr_Plan, FILE *FileDes, int &Nli, int &Nci, Bool Verbose= False);

// compression using multiresolution  and fixed budget

void mr_getquant(MultiResol & MR_Data,float *TabSigma, int budget, 
               float *TabQuantif, float *TabMean, float *TabRate=NULL,
	       Bool Verbose = False);
               
void mr_budcompress(Ifloat &Imag, type_transform Transform, 
                    int Nbr_Plan, int Budget, FILE *FileDes, Ifloat &Residual,
                    Bool KillNoise, float &Noise, float NSigma,
                    float SignalQuant, Bool SupIsol=DEFAULT_SUP_ISOL,
                    Bool SupDil = True, Bool SupNeg=False, 
                    int MedianWindowSize=DEF_IPYRMED_WINDOW_SIZE,
                    int NbrIter=0, Bool Verbose=False, Bool RecIma=True);
                                      
void mr_buddecompress (Ifloat &Imag, type_transform Transform, 
                       int Resol, FILE *FileDes, int &Nli, int &Nci, int NbrIterRec=1,
		       Bool Verbose=False);

// Quantify a value "Val" with a quantification parameter "Q"
inline int quant (float Val, float Q)
{
   int ValReturn;
   if (Val >= 0)  ValReturn = (int) (Val / Q +0.5);
   else  ValReturn = (int) (Val / Q - 0.5);
   return ValReturn;
}

// return the quantification parameter at scale s, knowing
// the quantification parameter at last and first scale by the low
// Q(s) = Q(F) * (1 + (Q(F) - Q(L)) / 2^s)
inline float get_quant_param1(int s, float QuantFirstScale, float QuantLastScale)
{
   float Diff;
   Diff =  QuantFirstScale - QuantLastScale;
   if (Diff > FLOAT_EPSILON)
             Diff =  QuantLastScale +  Diff / POW2(s);
   else Diff = QuantLastScale ;
   return Diff;
}
// return the quantification parameter at scale s, knowing
// the quantification parameter at last and first scale,
// and the number of scales N by the low
// Q(s) = Q(F) + (2^{N-1-s}-1)/(2^{N-1-1) * (Q(F) - Q(L))
inline float get_quant_param(int s, int NbrScale, 
                             float QuantFirstScale, float QuantLastScale)
{
   int s1 = NbrScale - 1 - s;
   float Diff;
   Diff =  QuantFirstScale - QuantLastScale;
   if (Diff > FLOAT_EPSILON)
         Diff = QuantLastScale + (float)(POW2(s1)-1) / (float)(POW2(NbrScale - 1)-1) * Diff;
   else Diff = QuantLastScale;
   return Diff;
}  
        
void lifting_compress(Ifloat &Imag, int Nbr_Plan, FILE *FileDes,  
                      float NSigma, float SignalQuantif,Bool UnifQuant=False,
                      Bool Verbose=False, Bool RecImag=False);
void lifting_decompress(Ifloat &Imag, int Resol, FILE *FileDes, 
                       int &Nli, int &Nci, Bool Verbose);

// return the size of an image at a lower resolution
inline int get_size_resol_ima(int N, int Resol, type_transform Transform)
{
  int s,Val=N;
  set_transform Set_Transform = SetTransform (Transform);
  switch (Set_Transform)
  {
     case TRANSF_PAVE: break;
     case TRANSF_PYR:
        for (s = 0; s < Resol; s++) Val = (Val+1) / 2;
        break;
        break;
     case TRANSF_MALLAT:
     case TRANSF_FEAUVEAU:
        for (s = 0; s < Resol; s++) Val = (Val+1) / 2;
        break;
     default: 
        cerr << "Error: unknowed transform in get_size_resol_ima " << endl;
	cerr << "StringTransform = " << StringTransform(Transform) << endl;
	cerr << "N = " << N << endl;
	cerr << "Resol = " << Resol<< endl;
	exit(-1);
	break;
  }
  return Val;
}

class BlockMRInfo;
                             
class MR_CompData {
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
    void dec_data();                   // last operation before saving the  results
    void set_block_info_data();   // set information when block decoding
    void upblock(Ifloat &Low);    // update the current block


    public:
    // image size
    int Nl,Nc;

    // use for compression with float
    Ifloat Dat;
    Ifloat Residual;
    Ifloat *BDat;
    
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
    Bool UseBlock;         // if True, image are compressed by block
    int  NbrBlock;         // Number of blocks
    int  BlockSize;        // Block size in both directions
    int  BlockOverlap;      // Overlap between two blocks
    int  L_Nl, L_Nc;       // size of the big image
    int NbrBlocNl, NbrBlocNc; // Number of blocks in both directions
    int *TabBlockNl, *TabBlockNc;
    int *TabBlockDepNl, *TabBlockDepNc;
    
    type_format FormatInput; // Input data format;
    type_quant TypeQuant ;  // Type of quantification
    type_coding TypeCoder;  // Coding type
    Bool UseBudget;         // if set, a budget is fixed
    float CompressionRatio; // wanted compression ratio
    unsigned int TargetNbrByte; // wanted number of bytes
    unsigned int EncodedNbrByte; // achieved number of bytes
    Bool NoiseInData;       // Data contains noise 
    Bool UseLiftingInsteadOfFilter; // Replace the 7/9 filter bank by the 7/9
                                    // filter decomposition
    Bool UnifQuant; // only used if UseLiftingInsteadOfFilter == True
                    // if true, use a uniform wuantization
    Bool Verbose;
    type_lift LiftingTrans; // Lifting transform parameter 
                            // (if the lifing transform is used).
   MR_CompData() 
   {
    UnifQuant = False;
    UseLiftingInsteadOfFilter = False;
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
    FormatInput = F_UNKNOWN;
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
   void up_mrresol(Ifloat &ImagLowResol, int NumBlock=1);
   void up_mrresol(Ifloat &ImagLowResol, BlockMRInfo &BV, BlockMRInfo &BVup);

   // write the result to the disk
   void write();
   
   ~MR_CompData();
};


class BlockMRInfo {
  int L_Nl, L_Nc;   // large image size
  int NbrBlock;     // number of blocks
  int BlockSize;    // block size
  int NbrBlockNl;   // number of block in lines
  int NbrBlockNc;   // number of block in columns
  
  int Resol;        // resolution 
  int Nl, Nc;       // image size at the resolution Resol
  intarray TabNl,TabNc; // size of each blocks at the resolution Resol
  intarray TabDepNl, TabDepNc; // Position of each block

  void init_param(int Nli, int Nci, int Bs);
  void set_pos();
  int FirstBlockNl, FirstBlockNc; // First and last block contained in
  int LastBlockNl, LastBlockNc;   // in the image

   public:
  Bool Verbose;
  int KeepResi;     // The noise map is also in the compressed file

  BlockMRInfo(int Nli, int Nci, int Bs);
  BlockMRInfo(BlockMRInfo &BI, int Scale);
  int check_size(int NbrScale, int MaxSizeLastScale);

  BlockMRInfo(int Indi, int Indj, BlockMRInfo &BV, int NlVisu, int NcVisu);
  BlockMRInfo(char *FileName);
  
    // block size
  int n_elem() const  { return Nl*Nc;}
  int resol() const  { return Resol;}
  
  int fb_nl() const  { return FirstBlockNl;}
  int lb_nl() const  { return LastBlockNl;}
  int fb_nc() const  { return FirstBlockNc;}
  int lb_nc() const  { return LastBlockNc;}

  // number of lines and columns of the block
  int nl() const  { return Nl;}
  int nc() const  { return Nc;}
  
  // number of blocks
  int nbr_block() {return NbrBlock;}
  int nbr_block_nl() const  { return  NbrBlockNl;}
  int nbr_block_nc() const  { return  NbrBlockNc;}
      
  // block size of a given block
  int block_nl(int Bi, int Bj) const  { return  TabNl(Bi,Bj);}
  int block_nc(int Bi, int Bj) const  { return  TabNc (Bi,Bj);}
  int pos_nl(int Bi, int Bj) const  { return  TabDepNl(Bi,Bj);}
  int pos_nc(int Bi, int Bj) const  { return  TabDepNc(Bi,Bj);}
  
  // number of lines and columns of the image
  int ima_nl() const  { return L_Nl;}
  int ima_nc() const  { return L_Nc;}
  int block_size() const  { return  BlockSize;}
  
  // image size
  int n_elem_ima() const  { return L_Nl*L_Nc;}
  
  void write(char *FileName);
  void read(char *FileName);
};

#endif
