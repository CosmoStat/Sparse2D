/*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_IO.h
**
*******************************************************************************
**
**    DESCRIPTION  FITS Include
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM_BIO_H_
#define _IM_BIO_H_

#include"Array.h"
#include"IM_Obj.h"
#include"IM_IOTools.h"
#include"IM_IO.h"


class IOInfoData {
  public:
    fitsstruct *PtrFits;
    midas_pic_des *PtrMidas;
    PICINFO *PtrGif;
    int Nl;
    int Nc;
    int Nima;
    type_format Format;
    type_data Type;
    IOInfoData()  { 
       PtrFits = NULL;
       PtrMidas = NULL;
       PtrGif = NULL;
       Nl = Nc = 0;
       Nima = 1;
       Format = F_UNKNOWN;
       Type = UNKNOWN;
    }
    void get_info(char *FileName);
    void put_info(char *FileName);
    ~IOInfoData() 
    { 
//      if (PtrFits != NULL)
// 	{
// 	   if (PtrFits->fitshead != NULL) delete PtrFits->fitshead;
// 	   delete PtrFits;
// 	}
        if (PtrMidas != NULL) delete PtrMidas; 
	if (PtrGif != NULL) delete PtrGif;
	PtrFits = NULL;
	PtrMidas = NULL;
	PtrGif = NULL;
    }
};

void io_read_block_ima(char *File_Name, Ifloat &Image, 
               int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale = False);
void io_read_block_ima(char *File_Name, Iint &Image, 
               int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale = False);
void io_write_block_ima(char *File_Name, Ifloat &Image, 
               int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale = False);
void io_write_block_ima(char *File_Name, Iint &Image, 
               int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale = False);
void io_close_block_ima(char *File_Name, IOInfoData &InfoDat);


/*********************************************************************/

class LIfloat {
  int NumBlock;    // block number in memory
  Ifloat Block;    // block data
  int Nl, Nc;      // block size;
  int Depi, Depj;  // starting points
  char *FileName;  // large image file name

  int NbrBlock;             // Number of blocks
  int BlockSize;            // Block size in both directions
  int BlockOverlap;         // Overlap between two blocks
  int L_Nl, L_Nc;           // size of the big image
  int NbrBlocNl, NbrBlocNc; // Number of blocks in both directions
  int *TabBlockNl, *TabBlockNc;
  int *TabBlockDepNl, *TabBlockDepNc;
  void set_tabblock(); 
  void reset_param()
  {
    BlockOverlap=0;
    TabBlockNl=NULL;
    TabBlockNc=NULL;
    TabBlockDepNl=NULL;
    TabBlockDepNc=NULL;
    NbrBlocNl=0;
    NbrBlocNc=0;
    Nl = 0;
    Nc = 0;
    Depi = 0;
    Depj = 0;
    L_Nl = 0;
    L_Nc = 0;
    NoBscale = False;
    Verbose = False;
    FileName = NULL;
    NumBlock = -1;
  } 
public:
  Bool NoBscale; // for fits file only. If true, Bscaling is not applied
                 // when reading and writing blocks
  Bool Verbose;  
  IOInfoData Info; // info on the file and the large image parameter
  LIfloat() 
  {
      reset_param();
  }
  type_format format_image() {return Info.Format;}
  
  void set_block(int Nb); // set to the class to block Nb
  void get_block(int Nb); // read a new block from the disk
  void next_block() { get_block(NumBlock+1);}
  void put_block();       // write the block to the disk
  
  Ifloat & block () { return Block;} // reference to block data
  float *buffer()   { return Block.buffer();} // pointer to the block buffer
  
  // get a pixel of a block
  inline float & operator() (int i)  const { return  Block(i);}
  inline float & operator() (int i, int j)  const  {return Block(i,j);}
  inline float operator() (int i, int j, type_border bord) 
                                        const {return Block(i,j,bord);}
  
  // block size
  int n_elem() const  { return Nl*Nc;}

  // number of lines and columns of the block
  int nl() const  { return Nl;}
  int nc() const  { return Nc;}
  
  // number of blocks
  int nbr_block() {return NbrBlock;}
  int nbr_block_nl() const  { return  NbrBlocNl;}
  int nbr_block_nc() const  { return  NbrBlocNc;}
    
  // block number in memory
  int which_block() {return NumBlock;}
  
  // block size of a given block
  int block_nl(int Nb) const  { return  TabBlockNl[Nb];}
  int block_nc(int Nb) const  { return  TabBlockNc[Nb];}
  
  // number of lines and columns of the image
  int ima_nl() const  { return L_Nl;}
  int ima_nc() const  { return L_Nc;}
  
  // image size
  int n_elem_ima() const  { return L_Nl*L_Nc;}
  
  // set the class from an existing image file
  void init(char *FileImageName, int BSize);
  
  // set the class and create the file
  // if Info is set before, this information is taken into account
  void create(char *FileImageName, int NbrL, int NbrC, int BlockSize);
  
  // set the class and create the file
  // make a copy from IData to Info (!!!but pointer to structure  are also
  // copied, and the structures are not duplicated)
  void create(char *FileImageName, IOInfoData &IData, int BSize);
  
  void close(); // to be done when the file has been created or updated
  ~LIfloat();
 };

/*********************************************************************/

class LIint {
  int NumBlock;    // block number in memory
  Iint Block;    // block data
  int Nl, Nc;      // block size;
  int Depi, Depj;  // starting points
  char *FileName;  // large image file name

  int NbrBlock;             // Number of blocks
  int BlockSize;            // Block size in both directions
  int BlockOverlap;         // Overlap between two blocks
  int L_Nl, L_Nc;           // size of the big image
  int NbrBlocNl, NbrBlocNc; // Number of blocks in both directions
  int *TabBlockNl, *TabBlockNc;
  int *TabBlockDepNl, *TabBlockDepNc;
  void set_tabblock(); 
  void reset_param()
  {
    BlockOverlap=0;
    TabBlockNl=NULL;
    TabBlockNc=NULL;
    TabBlockDepNl=NULL;
    TabBlockDepNc=NULL;
    NbrBlocNl=0;
    NbrBlocNc=0;
    Nl = 0;
    Nc = 0;
    Depi = 0;
    Depj = 0;
    L_Nl = 0;
    L_Nc = 0;
    NoBscale = False;
    Verbose = False;
    FileName = NULL;
    NumBlock = -1;
  } 
public:
  Bool NoBscale; // for fits file only. If true, Bscaling is not applied
                 // when reading and writing blocks
  Bool Verbose;  
  IOInfoData Info; // info on the file and the large image parameter
  LIint() 
  {
      reset_param();
  }
  type_format format_image() {return Info.Format;}
  
  void set_block(int Nb); // set to the class to block Nb
  void get_block(int Nb); // read a new block from the disk
  void next_block() { get_block(NumBlock+1);}
  void put_block();       // write the block to the disk
  
  Iint & block () { return Block;} // reference to block data
  int *buffer()   { return Block.buffer();} // pointer to the block buffer
  
  // get a pixel of a block
  inline int & operator() (int i)  const { return  Block(i);}
  inline int & operator() (int i, int j)  const  {return Block(i,j);}
  inline int operator() (int i, int j, type_border bord) 
                                        const {return Block(i,j,bord);}
  
  // block size
  int n_elem() const  { return Nl*Nc;}

  // number of lines and columns of the block
  int nl() const  { return Nl;}
  int nc() const  { return Nc;}
  
  // number of blocks
  int nbr_block() {return NbrBlock;}
  int nbr_block_nl() const  { return  NbrBlocNl;}
  int nbr_block_nc() const  { return  NbrBlocNc;}
    
  // block number in memory
  int which_block() {return NumBlock;}
  
  // block size of a given block
  int block_nl(int Nb) const  { return  TabBlockNl[Nb];}
  int block_nc(int Nb) const  { return  TabBlockNc[Nb];}
  
  // number of lines and columns of the image
  int ima_nl() const  { return L_Nl;}
  int ima_nc() const  { return L_Nc;}
  
  // image size
  int n_elem_ima() const  { return L_Nl*L_Nc;}
  
  // set the class from an existing image file
  void init(char *FileImageName, int BSize);
  
  // set the class and create the file
  // if Info is set before, this information is taken into account
  void create(char *FileImageName, int NbrL, int NbrC, int BlockSize);
  
  // set the class and create the file
  // make a copy from IData to Info (!!!but pointer to structure  are also
  // copied, and the structures are not duplicated)
  void create(char *FileImageName, IOInfoData &IData, int BSize);
  
  void close(); // to be done when the file has been created or updated
  ~LIint();
 };
 

#endif
