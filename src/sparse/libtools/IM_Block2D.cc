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
**    File:  IM_Block2D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION 2D block manadgement
**    -----------  
**                 
******************************************************************************/

#include "IM_Block2D.h"

/***********************************************************************/

void Block2D::alloc(int Nlc, int Ncc,  int ParamBlockSize, float Overlapping)
{
	Overlap = Overlapping;
	if(Overlap>0) BlockOverlap=True;
	alloc(Nlc, Ncc,  ParamBlockSize);
}

/***********************************************************************/

void Block2D::alloc(int Nlc, int Ncc, int ParamBlockSize)
{
    int MaxXY;
    Nl=Nlc;
    Nc=Ncc;
    BlockSize = ParamBlockSize;
    

	if ((BlockOverlap == True) && (Overlap==0)) Overlap=0.5;
    MaxXY = MAX(Nl,Nc);
        
    if ((BlockSize==0) || (BlockSize==MaxXY))
    {
       Nlb = Ncb = 1;
       BlockSize = MaxXY;
       BlockOverlap=False;
       Overlap = 0;
       OverlapSize = 0;
    }
    else
    {
        BlockImaSize = ceil(BlockSize*(1-Overlap));
		OverlapSize = BlockSize-BlockImaSize;

       // cout << "BlockSize = " << BlockSize << endl;
 
       // cout << "BlockImaSize = " << BlockImaSize << endl;
        int ix=Nl-((Nl-OverlapSize)/BlockImaSize)*BlockImaSize;
		Nlb = ((Nl-OverlapSize)/BlockImaSize) + (ix!=OverlapSize);
		int iy=Nc-((Nc-OverlapSize)/BlockImaSize)*BlockImaSize;
		Ncb = ((Nc-OverlapSize)/BlockImaSize) + (iy!=OverlapSize);

       // Nlb =  (Nl %  BlockImaSize == 0) ? Nl / BlockImaSize: Nl / BlockImaSize + 1;
       // Ncb =  (Nc %  BlockImaSize == 0) ? Nc / BlockImaSize: Nc / BlockImaSize + 1;
    }
    
    NbrBlock = Nlb*Ncb;
    // if (BlockOverlap == False) BlockImaSize = 0;
    if (Verbose == True)
    {
       printf( "\nima size = (%2d,%2d), BlockSize = %2d, BlockImaSize = %2d, OverlapSize = %2d\n", Nl,Nc,BlockSize, BlockImaSize, OverlapSize);
       printf( "Nbr Blocks = (%2d,%2d)\n", Nlb,Ncb);
       printf( "NbrBlocks*BlockImaSize + OverlapSize = (%2d,%2d)\n", Nlb*BlockImaSize+OverlapSize, Ncb*BlockImaSize+OverlapSize);
    }
}

/***********************************************************************/

double Block2D::get_weight(int PosPix, int CurrentBlock, int MaxBlockNumber) 
{
	double ValReturn;
	double ValL=1.,ValR=1.;
    
	if((PosPix < OverlapSize) && (CurrentBlock > 0))
		ValL = (double) pow(sin((double) (PosPix+1) 
                               / (double) OverlapSize * PI/2.), 2.);
    
	if((PosPix > BlockImaSize-1) && (CurrentBlock < MaxBlockNumber-1))
		ValR = (double) pow(cos((double) (PosPix-(BlockImaSize-1)) 
                               / (double) OverlapSize * PI/2.), 2.);
	
	ValReturn = ValL*ValR;
	
	return ValReturn;
}

/*************************************************************************/


void Block2D::get_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock, Bool Weight)
// Extract the block (Bi,Bj) from Ima and put it in ImaBlock
{
   int k,l;
   int Depi = BlockImaSize*Bi;
   int Depj = BlockImaSize*Bj;
    
   if ((WeightFirst == False) || (Weight == False))
   {
      // cout << "get_block_ima NO Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++) 
      {
          int Indk = test_index_mirror(Depi+k,Nl);
          int Indl = test_index_mirror(Depj+l,Nc);
          ImaBlock(k,l) = Ima(Indk,Indl);  
      }
   }
   else
   {
      // cout << "get_block_ima  Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
      {
          int Indk = test_index_mirror(Depi+k,Nl);
          int Indl = test_index_mirror(Depj+l,Nc);
          ImaBlock(k,l) = Ima(Indk,Indl) * get_weight(k, Bi, Nlb) * get_weight(l, Bj, Ncb);  
      }
   }
   // cout << "BLOCK: " << Bi <<","<<Bj<<","<<Bk<<": Dep = "<< Depi<<","<<Depj<<","<<Depk<<endl;
}

/*************************************************************************/

void Block2D::put_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Extract the block (Bi,Bj,Bk) from ima and put it in imaBlock
{
   int k,l;
   int Depi = BlockImaSize*Bi;
   int Depj = BlockImaSize*Bj;
 
   for (k = 0; k < BlockSize; k++)
   for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                               Ima(Depi+k,Depj+l) = ImaBlock(k,l);   
   // cout << "BLOCK: " << Bi <<","<<Bj<<","<<Bk<<": Dep = "<< Depi<<","<<Depj<<","<<Depk<<endl;
}

/*************************************************************************/

void Block2D::add_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock, Bool Weight)
// Add the block (Bi,Bj) ImaBlock in Ima  with weighted values
{
   int k,l;
    int Depi = BlockImaSize*Bi;
    int Depj = BlockImaSize*Bj;
    
   if ((WeightFirst == False) || (Weight == True))
   {
     // cout << "add_block_ima Weight " << endl;
     // if (BlockOverlap == True) cout << "BlockOverlap True" << endl;
     // if (WeightFirst == False) cout << "WeightFirst False" << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                   Ima(Depi+k,Depj+l) += ImaBlock(k,l)*
	             get_weight(k, Bi, Nlb) *  get_weight(l, Bj, Ncb);  
   }
   else
   {
      // cout << "add_block_ima NO Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                   Ima(Depi+k,Depj+l) += ImaBlock(k,l);  
   }
}

/****************************************************************************/
