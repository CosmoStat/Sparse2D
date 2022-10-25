/******************************************************************************
**                   Copyright (C) 2003 by CEA 
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  28/03/2003 
**    
**    File:  IM_Block2D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition for 2D block manadgement
**    ----------- 
******************************************************************************/

#ifndef _BLOCK2D_H_
#define _BLOCK2D_H_

#include "IM_Obj.h"
  
/***********************************************************************/

class Block2D {   
   private:
   int Nc, Nl;
   int BlockImaSize;     // Ima block size  without taking account
                          // the overlapping aera
   int BlockSize;         // Ima block size 	taking account
                          // the overlapping aera		  
   int Ncb,Nlb;       // Number of blocks in the x,y directions
   int NbrBlock;   
   void reset_param() {Nc=Nl=0;BlockImaSize=BlockSize=NbrBlock; 
                       Ncb=Nlb=0;BlockOverlap=False;
       Verbose=False;WeightFirst=False;	Overlap=0.5;OverlapSize=0;}   

   public:
   Bool Verbose;
   Bool BlockOverlap;
   Bool WeightFirst;
   float Overlap;
   int OverlapSize;
   inline int nl() { return Nl;}
   inline int nc() { return Nc;}
   inline int nbr_block_nl() { return Nlb;}
   inline int nbr_block_nc() { return Ncb;}
   inline int nbr_block()    { return NbrBlock;}
   inline int block_size()   { return BlockSize;}
   inline int l_nl() {return BlockSize*Nlb;}
   inline int l_nc() {return BlockSize*Ncb;}
   Block2D() {reset_param();} 
   void alloc(int Nlc, int Ncc, int ParamBlockSize);
   void alloc(int Nlc, int Ncc, int ParamBlockSize, float Overlaping);
   void get_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock, Bool Weight=False);
   void put_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock);
   double get_weight(int PosPix, int CurrentBlock, int MaxBlockNumber); 
   void add_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock, Bool Weight=False);
   ~Block2D() {reset_param();} 
};

/***********************************************************************/

#endif
