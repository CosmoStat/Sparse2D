/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  11/06/03 
**    
**    File:  IM1D_Dct.cc
**
**    Modification history :
**
******************************************************************************/

#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM1D_Block.h"
#include "IM1D_Dct.h"

/*********************************************************************/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
//
//                         Block1D
//
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*********************************************************************/

/*********************************************************************/
// Block1D::alloc          
/*********************************************************************/

/*
void Block1D::alloc(int NSig, int ParamBlockSize, float Overlapping, Bool UseBlockSizePow2)
{
	Overlap = Overlapping;
	if(Overlap>0) _BlockOverlap=True;
	alloc(NSig, ParamBlockSize);
}
*/
/***********************************************************************/

void Block1D::alloc (int NSig, int ParamBlockSize, Bool UseBlockSizePow2)
{
    _Nx=NSig;
    _BlockSize = ParamBlockSize;
    BLOCKSIZE_Power_of_2 = UseBlockSizePow2;

    if (_BlockSize==0) {
       _Nxb = 1;
       if (is_power_of_2(_Nx)) _BlockSize = _Nx;
       else {
          int i=1;
          while (i < _Nx) i*=2;
      _BlockSize = i;
       }
       _BlockOverlap=False;
    } else if (_BlockSize==_Nx) {
       _Nxb = 1;
       int i=1;
       if (is_power_of_2(_Nx)) _BlockSize = _Nx;
       else {
          while (i < _Nx) i*=2;
         _BlockSize = i;
       }
       _BlockOverlap=False;
    } else {
       // cout << "BlockSize = " << BlockSize << endl;
       if (_BlockOverlap == True){
               _BlockSigSize = (_BlockSize+1)/2;
               _Nxb =  (_Nx %  _BlockSigSize == 0) ?
                       _Nx / _BlockSigSize - 1 : _Nx / _BlockSigSize + 1 ;
       
       }
       else{ _BlockSigSize = _BlockSize;
             _Nxb =  (_Nx %  _BlockSigSize == 0) ?
                     _Nx / _BlockSigSize : _Nx / _BlockSigSize + 1 ;
       
       }
       if (!(is_power_of_2(_BlockSize))){
          int i=1, b=_BlockSize;
          while (i < b) i*=2;
         _BlockSize = i;
       }
    }
    
    _NbrBlock = _Nxb;
    if (_BlockOverlap == False) _BlockSigSize = 0;
    
    if (_Verbose == True) {
       cout << "Sig size = " << _Nx << ", BlockSize = " << _BlockSize << endl;
       cout << "Nbr Blocks = " << _Nxb << endl;
       if (_BlockOverlap == True)
          cout << "Block Overlap: size = " << _BlockSigSize << endl;
       else cout << "No Block Overlap" << endl;
    }
}
/* PB code introduced in 2013 for dictionary learning, make cb_mca1d having a bad behavior
 problem is not undertood.
 
void Block1D::alloc (int NSig, int ParamBlockSize, Bool UseBlockSizePow2) 
{

    _Nx=NSig;
    _BlockSize = ParamBlockSize;
    BLOCKSIZE_Power_of_2 = UseBlockSizePow2;
            
	if ((_BlockOverlap == True) && (Overlap==0)) Overlap=0.5;
    if ((_BlockSize==0) || (_BlockSize==_Nx))
    {
       _Nxb = 1;
       _BlockSize = _Nx;
       _BlockOverlap=False;
        Overlap = 0;
        OverlapSize = 0;
        
        if (BLOCKSIZE_Power_of_2 == True)
        {
           if (is_power_of_2(_Nx)) _BlockSize = _Nx;
           else 
           {
              int i=1;
              while (i < _Nx) i*=2;
	         _BlockSize = i;
           }
        }  
    }
    else 
    {
        _BlockSigSize = ceil(_BlockSize*(1-Overlap));
        OverlapSize = _BlockSize-_BlockSigSize;
        int ix=NSig-((NSig-OverlapSize)/_BlockSigSize)*_BlockSigSize;
		_Nxb = ((NSig-OverlapSize)/_BlockSigSize) + (ix!=OverlapSize);
        
        if (BLOCKSIZE_Power_of_2 == True)
        {
           if (!(is_power_of_2(_BlockSize)))
           {
              int i=1, b=_BlockSize;
              while (i < b) i*=2;
              _BlockSize = i;
              _BlockSigSize = ceil(_BlockSize*(1-Overlap));
              OverlapSize = _BlockSize-_BlockSigSize;
              int ix=NSig-((NSig-OverlapSize)/_BlockSigSize)*_BlockSigSize;
              _Nxb = ((NSig-OverlapSize)/_BlockSigSize) + (ix!=OverlapSize);
           }
        }
     }
    
    _NbrBlock = _Nxb;
 
   //  if (_Verbose == True)
    {
       cout << "InParamBlockSize = " << ParamBlockSize << ". Sig size = " << _Nx << ", BlockSize = " << _BlockSize << endl;
       cout << "Nbr Blocks = " << _Nxb << endl;
       if (_BlockOverlap == True) 
          cout << "Block Overlap: size = " << OverlapSize / 2. << ", BlockSigSize = " <<  _BlockSigSize << endl;
       else cout << "No Block Overlap" << endl;
    }
}
*/
/*********************************************************************/
// Block1D::get_weight        
/*********************************************************************/
double Block1D::get_weight (int PosPix, int CurrentBlock, int MaxBlockNumber) {

    double ValReturn=1.;
    if ( (_BlockSize % 2 == 0) ||
         ( (PosPix != 0)  && (PosPix != _BlockSize-1))) {
         
       int Center = _BlockSize / 2;
       if ((CurrentBlock != 0)  &&  (PosPix < Center))
           ValReturn =   (float) pow(sin((double) (PosPix)
                       / (double) Center*PI/2.), 2.);
       else if ((CurrentBlock != MaxBlockNumber-1)  &&  (PosPix > Center))
           ValReturn =   (float) pow(cos((double)(PosPix-Center)
                       / (double) Center *PI/2.), 2.);
    
    } else {
    
       if  ((CurrentBlock != 0) && (PosPix==0))
          ValReturn=0.;
       else if ((CurrentBlock != MaxBlockNumber-1) && (PosPix == _BlockSize-1))
      ValReturn=0.;
    }
    return ValReturn;
}

/* PB code introduced in 2013 for dictionary learning, make cb_mca1d having a bad behavior
 problem is not undertood.
double Block1D::get_weight (int PosPix, int CurrentBlock, int MaxBlockNumber)
{
	double ValReturn;
	double ValL=1.,ValR=1.;
    
	if((PosPix < OverlapSize) && (CurrentBlock > 0))
		ValL = (double) pow(sin((double) (PosPix+1) 
                               / (double) OverlapSize * PI/2.), 2.);
    
	if((PosPix > _BlockSigSize-1) && (CurrentBlock < _NbrBlock-1))
		ValR = (double) pow(cos((double) (PosPix-(_BlockSigSize-1)) 
                               / (double) OverlapSize * PI/2.), 2.);
	
	ValReturn = ValL*ValR;
	
	return ValReturn;
    
}
*/
/*********************************************************************/
// Block1D::get_block_sig        
/*********************************************************************/
/* PB code introduced in 2013 for dictionary learning, make cb_mca1d having a bad behavior
 problem is not undertood.
void Block1D::get_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock, Bool Weight) 
{
// Extract the block (Bi) from SIg and put it in SigBlock

   int Depx = _BlockSigSize*Bx;
   
    if ((_WeightFirst == False) || (Weight == False))
    {
          // cout << "get_block_ima NO Weight " << endl;
        for (int k = 0; k < _BlockSize; k++) 
        {
            int Indk = test_index_mirror(Depx+k,_Nx);
            SigBlock(k) = Sig(Indk);  
        }
   
   } 
   else 
   {
      // cout << "get_block_ima  Weight " << endl;
      for (int k = 0; k < _BlockSize; k++)      
      {
          int Indk = test_index_mirror(Depx+k,_Nx);
          SigBlock(k) = Sig(Indk) * get_weight(k, Bx, _Nxb);  
      }
   }
   //cout << "BLOCK: " << Bx <<", Dep = "<< Depx << endl;
}
*/
void Block1D::get_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock, Bool Weight) {
// Extract the block (Bi) from SIg and put it in SigBlock

   int Depx = (_BlockSize-_BlockSigSize)*Bx;
   type_border Bord = I_MIRROR;
   
   if ((_BlockOverlap == False) || (_WeightFirst == False)) {
      // cout << "get_block_ima NO Weight " << endl;
      for (int k = 0; k < _BlockSize; k++)
         SigBlock(k) = Sig(Depx+k,Bord);
   
   } else {
      // cout << "get_block_ima  Weight " << endl;
      for (int k = 0; k < _BlockSize; k++)
         SigBlock(k) = Sig(Depx+k,Bord) * get_weight(k, Bx, _Nxb);
   }
   //cout << "BLOCK: " << Bx <<", Dep = "<< Depx << endl;
}

/*********************************************************************/
// Block1D::put_block_sig          
/*********************************************************************/
/* PB code introduced in 2013 for dictionary learning, make cb_mca1d having a bad behavior
 problem is not undertood.
void Block1D::put_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock) 
{
// Extract the block (Bi) from sig and put it in SigBlock

    int Depx = _BlockSigSize*Bx;
 
   for (int k=0; k<_BlockSize; k++)
      if ((Depx+k >= 0) && (Depx+k < _Nx))   Sig(Depx+k) = SigBlock(k);

   //cout << "BLOCK: " << Bx <<", Dep = "<< Depx << endl;
}
*/

void Block1D::put_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock) {
// Extract the block (Bi) from sig and put it in SigBlock

   int Depx = (_BlockSize-_BlockSigSize)*Bx;
 
   for (int k=0; k<_BlockSize; k++)
      if ((Depx+k >= 0) && (Depx+k < _Nx))
         Sig(Depx+k) = SigBlock(k);

   //cout << "BLOCK: " << Bx <<", Dep = "<< Depx << endl;
}

/*********************************************************************/
// Block1D::add_block_sig       
/*********************************************************************/
/* PB code introduced in 2013 for dictionary learning, make cb_mca1d having a bad behavior
   problem is not undertood.
void Block1D::add_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock, Bool Weight) 
{
// Add the block (Bx) SigBlock in Sig  with weighted values

    int Depx = _BlockSigSize*Bx;
    
    if ((_WeightFirst == False) || (Weight == True))
    {   
      // cout << "add_block_sig Weight " << endl;
      // if (_BlockOverlap == True) cout << "_BlockOverlap True" << endl;
      // if (_WeightFirst == False) cout << "_WeightFirst False" << endl;
      for (int k=0; k<_BlockSize; k++)
         if ((Depx+k >= 0) && (Depx+k < _Nx)) 
            Sig(Depx+k) += SigBlock(k)*get_weight(k, Bx, _Nxb);  
   
   } 
   else 
   {
      // cout << "add_block_sig NO Weight " << endl;
      for (int k=0; k<_BlockSize; k++)
         if ((Depx+k >= 0) && (Depx+k < _Nx)) 
            Sig(Depx+k) += SigBlock(k);  
   }
}
 */

void Block1D::add_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock, Bool Weight) {
// Add the block (Bx) SigBlock in Sig  with weighted values

   int Depx = (_BlockSize-_BlockSigSize)*Bx;
    
   if ((_BlockOverlap == True) && (_WeightFirst == False)) {
   
      // cout << "add_block_sig Weight " << endl;
      // if (_BlockOverlap == True) cout << "_BlockOverlap True" << endl;
      // if (_WeightFirst == False) cout << "_WeightFirst False" << endl;
      for (int k=0; k<_BlockSize; k++)
         if ((Depx+k >= 0) && (Depx+k < _Nx))
            Sig(Depx+k) += SigBlock(k)*get_weight(k, Bx, _Nxb);
   
   } else {
   
      // cout << "add_block_sig NO Weight " << endl;
      for (int k=0; k<_BlockSize; k++)
         if ((Depx+k >= 0) && (Depx+k < _Nx))
            Sig(Depx+k) += SigBlock(k);
   }
}


/*********************************************************************/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
//
//                         LOCAL_DCT1D
//
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*********************************************************************/


/*********************************************************************/
// LOCAL_DCT1D::set_dct   
/*********************************************************************/
void LOCAL_DCT1D::set_dct (fltarray& DCTData) {

   if (DCTData.nl() != _Nxt) {
      cout << "Error: Bad signal size in LOCAL_DCT1D::set_dct ... " << endl;
      cout << "       signal size is: Nx = " << DCTData.nx() << endl;
      cout << "       Expected size is: Nx = " << _Nxt << endl;
      exit(-1);
   }
   _DCTSig = DCTData;
}

/*********************************************************************/
// LOCAL_DCT1D:: recons  
/*********************************************************************/
void LOCAL_DCT1D::recons(fltarray& Sig) {

   if (_Nx == 0) {
      cout << "Error: LOCAL_DCT1D class is not allocated ... " << endl;
      exit(-1);
   }
   
   if (Sig.nx() != _Nx)  Sig.resize(_Nx);
   sig_invblockdct (_DCTSig, Sig, _B1D, _B1DTrans);
}

/*********************************************************************/
// LOCAL_DCT1D::transform  
/*********************************************************************/
void LOCAL_DCT1D::transform(fltarray& Sig) {

   if (_Nx == 0) {
      cout << "Error: LOCAL_DC1D class is not allocated ... " << endl;
      exit(-1);
   }
   
   if (Sig.nx() != _Nx) {
      cout << "Error: Bad image size in LOCAL_DCT1D::transform ... " << endl;
      cout << "       signal size is: Nx = " << Sig.nx() << endl;
      cout << "       Expected size is: Nx = " << _Nxt << endl;
      exit(-1);
   }
   sig_blockdct (Sig, _DCTSig, _B1D, _B1DTrans);
}

/*********************************************************************/
// LOCAL_DCT1D::threshold  
/*********************************************************************/
void LOCAL_DCT1D::threshold (float Lamba_Sigma) {

   int BS = _B1DTrans.block_size();
   fltarray BlockCosSig(BS, "blocksig");
   float Coef, Level = Lamba_Sigma*SigmaNoise*COS_Sensibility;
   float Noise = SigmaNoise * norm();

   if (SuppressDctCompCont){
      for (int bx=0; bx<_B1DTrans.nbr_block_nx(); bx++) {
         _B1DTrans.get_block_sig (bx, _DCTSig, BlockCosSig);
         BlockCosSig(0) = 0;
	 _B1DTrans.put_block_sig(bx, _DCTSig, BlockCosSig);
      } 
   }

   if ((COSSupportMin < 0) || (COSSupportMin >= 0.5)) COSSupportMin = 0;
   COSSupportMin *= BS;
   if (COSSupportMin == 0) {
      for (int i=0; i < _DCTSig.nx(); i++) {
	 Coef = _DCTSig(i);
	 _DCTSig(i) = update(Coef, Level, Noise);
      }
   } else {
   
      //cout << "COSSupportMin = " << COSSupportMin << " BS = " << BS << endl;
      //INFO(LDCT.DCTSig, "BEF");
      for (int bx=0; bx<_B1DTrans.nbr_block_nx(); bx++) {
      
         _B1DTrans.get_block_sig (bx, _DCTSig, BlockCosSig);
         
         for (int i=0; i < BS; i++) {
	     int Indi = (i < BS/2) ? i+BS/2: i - BS/2;
             int Rad = Indi - BS/2;
	     if (BlockCosSig(i)<0) BlockCosSig(i)=0;
             BlockCosSig(i) =  update (BlockCosSig(i), Level, Noise);
	     
	     
	     //cout << " Rad = " << Rad  << " COSSupportMin = " << COSSupportMin << endl;
             if (Rad < COSSupportMin) BlockCosSig(i) = 0; 
 	 }
	 _B1DTrans.put_block_sig(bx, _DCTSig, BlockCosSig);
      }
      //INFO(LDCT.DCTSig, "AFTER");
   }
   
   if (Verbose) {
      int NumberDetectedCoef=0;
      for (int i=0;i<_DCTSig.nx();i++)
         if (fabs(_DCTSig(i)) > 0) NumberDetectedCoef++;
      cout << "Number detected coef : " << NumberDetectedCoef << endl;
   }

}

/*********************************************************************/
// LOCAL_DCT1D::alloc_from_trans   
/*********************************************************************/
void LOCAL_DCT1D::alloc_from_trans (int Nxi, int BS, Bool Overlapping, 
                                    Bool WeightF) {
    int Nx1 = Nxi;
    if (Overlapping == True)
       Nx1 /= 2;
    alloc(Nx1, BS, Overlapping, WeightF);
}

/*********************************************************************/
// LOCAL_DCT1D::alloc_from_trans   
/*********************************************************************/
void LOCAL_DCT1D::reset() {
   _Nx = _Nxt = 0;
   _BlockSize = 0;
   _Overlap = False;
}

/*********************************************************************/
// LOCAL_DCT1D::alloc   
/*********************************************************************/
void LOCAL_DCT1D::alloc(int Nxi, int BS, Bool Overlapping, Bool WeightF) {
   
   _Nx= Nxi;
   _BlockSize = BS;
   _Overlap = Overlapping;
   if (_Nx==0) {
      cout << "Error: LOCAL_DCT1D class is allocated for an signal of zero dimension ... " << endl;
      cout << "       Image size is: Nx = " << _Nx << endl;
      exit(-1);
   }
   
   if ((_BlockSize > _Nx) || (_BlockSize <= 0)) _BlockSize = _Nx;
   _B1D._BlockOverlap = Overlapping;
   _B1D._WeightFirst = WeightF;
   _B1D.alloc (_Nx, _BlockSize);
   _BlockSize = _B1D.block_size();
   _Nxt = _B1D.nbr_block_nx() * _BlockSize;
   _B1DTrans.alloc (_Nxt, _BlockSize);
    _DCTSig.alloc (_Nxt, "DCTSig");
}
/*********************************************************************/

float LOCAL_DCT1D::update (float CoefSol, float Threshold, 
                                  float SoftLevel) {
// L2 norm minimization with L1 norm penalt
   
   float NewCoef = hard_threshold(CoefSol, Threshold);
   if (UseNormL1 == True) {
       float SoftL =  SoftLevel*NSigmaSoft;
       float T2 = Threshold/2.;
       if (SoftL < T2)
            NewCoef = soft_threshold(NewCoef, SoftL);
       else NewCoef = soft_threshold(NewCoef, T2);
   }
   return NewCoef; 
}

/*********************************************************************/
// LOCAL_DCT1D,   sig_invblockdct    
/*********************************************************************/
void sig_invblockdct (fltarray& Trans, fltarray& Sig, int BlockSize,
                     Bool Overlap, Bool WeightFirst) {
                     
   int Nx = Trans.nx();
   if (Overlap == True) {
      Nx /=2;
   }
   sig_invblockdct (Trans, Sig, BlockSize, Nx, Overlap, WeightFirst);
}

/*********************************************************************/
// LOCAL_DCT1D,    sig_dct  
/*********************************************************************/
void sig_dct ( fltarray& Sig, fltarray& Trans, Bool Reverse) {

    int Nx = Sig.nx();
    float Norm = sqrt(2./float(Nx));
    int isign = (Reverse == False) ? 1: -1;
    float *Ptr = Trans.buffer() - 1;
    Trans = Sig;
    cosft2 (Ptr, Nx, isign);
    
    for (int i=0; i<Nx; i++) Trans(i) *= Norm;     
}

/*********************************************************************/
// LOCAL_DCT1D,     sig_invblockdct 
/*********************************************************************/
void sig_invblockdct (fltarray& Trans, fltarray& Sig,  Block1D& B1D, 
                      Block1D& B1DTrans) {
                      
   fltarray SigBlock, Block, TransBlock;
   int Nx = B1D.nx();
   int BS = B1D.block_size();
   int Nxt = B1D.nbr_block_nx() * BS;
   if (Trans.nx() != Nxt) {
       cout << "Error in sig_invblockdct: Nxt = " <<  Nxt << endl;
       cout << "   Expected values are Nxt = " <<  Trans.nx() << endl;
       exit(-1);
   }
   if (Sig.nx() != Nx) Sig.resize(Nx);
   SigBlock.alloc(BS,"SigBlock");
   TransBlock.alloc(BS,"SigTransBlock");
   Sig.init();
   for (int i=0; i<B1D.nbr_block_nx(); i++) {
      B1DTrans.get_block_sig (i, Trans, TransBlock);
      sig_dct (TransBlock, SigBlock, True);
      B1D.add_block_sig (i, Sig, SigBlock);
   }
}

/*********************************************************************/
// LOCAL_DCT1D,      sig_invblockdct
/*********************************************************************/
void sig_invblockdct (fltarray& Trans, fltarray& Sig, int BlockSize, int Nx, 
                      Bool Overlap, Bool WeightFirst) {
                      
   Block1D B1D,B1DTrans;
   fltarray SigBlock, Block, TransBlock;
   
   B1D._Verbose = False;
   B1D._BlockOverlap = Overlap;
   B1D._WeightFirst = WeightFirst;
   B1D.alloc (Nx, BlockSize);
   int BS = B1D.block_size();
   int Nxt = B1D.nbr_block_nx() * BS;
   if (Trans.nx() != Nxt) {
       cout << "Error in sig_invblockdct: Nxt = " <<  Nxt << endl;
       cout << "   Expected values are Nlt = " <<  Trans.nx() << endl;
       exit(-1);
   }
   if (Sig.nx() != Nx) Sig.resize(Nx);
   SigBlock.alloc (BS, "SigBlock");
   TransBlock.alloc (BS, "SigTransBlock");
   B1DTrans.alloc (Nxt, BlockSize);
   
   Sig.init();
   for (int i=0; i<B1D.nbr_block_nx(); i++) {
      B1DTrans.get_block_sig (i, Trans, TransBlock);
      sig_dct (TransBlock, SigBlock, True);
      B1D.add_block_sig (i, Sig, SigBlock);
   }
}

/*********************************************************************/
// LOCAL_DCT1D,  sig_blockdct    
/*********************************************************************/
void sig_blockdct (fltarray& Sig, fltarray& Trans, Block1D & B1D, 
                   Block1D & B1DTrans) {
                  
   fltarray SigBlock, Block, TransBlock;
   int BS = B1D.block_size();
   int Nxt = B1D.l_nx();
   
   if (Trans.nx() != Nxt) Trans.resize(Nxt);
   SigBlock.alloc (BS, "SigBlock");
   TransBlock.alloc (BS, "SigTransBlock");   
   for (int i=0; i<B1D.nbr_block_nx(); i++) {
      B1D.get_block_sig (i, Sig, SigBlock);
      sig_dct (SigBlock, TransBlock);
      B1DTrans.put_block_sig (i, Trans, TransBlock);
   }
}

/*********************************************************************/
// LOCAL_DCT1D,   sig_blockdct   
/*********************************************************************/
void sig_blockdct (fltarray& Sig, fltarray& Trans, int BlockSize, 
                   Bool Overlap, Bool WeightFirst) {
                   
   int Nx = Sig.nx();
   Block1D B1D,B1DTrans;
   fltarray SigBlock, Block, TransBlock;
   
   B1D._Verbose = False;
   B1D._BlockOverlap = Overlap;
   B1D._WeightFirst = WeightFirst;
   B1D.alloc(Nx, BlockSize);
   int BS = B1D.block_size();
   int Nxt = B1D.nbr_block_nx() * BS;
   
   if (Trans.nx() != Nxt) Trans.resize(Nxt);
   SigBlock.alloc(BS, "SigBlock");
   TransBlock.alloc (BS, "SigTransBlock");
   B1DTrans.alloc (Nxt, BlockSize);
   
   for (int i=0; i<B1D.nbr_block_nx(); i++) {
      if (WeightFirst == False) B1D.get_block_sig (i, Sig, SigBlock);
      else B1D.get_block_sig (i, Sig, SigBlock);
      sig_dct (SigBlock, TransBlock);
      B1DTrans.put_block_sig (i, Trans, TransBlock);
   }
}

/*********************************************************************/
// LOCAL_DCT1D  getAbsMaxTransf   
/*********************************************************************/
float LOCAL_DCT1D::getAbsMaxTransf (fltarray& Sig) {

   float MaxLevel = _DCTSig.max();
   float NormMaxLevel = MaxLevel/SigmaNoise/norm()/COS_Sensibility;
                
   if (Verbose)
      cout << "Lambda dct  : " << NormMaxLevel << endl;
      
   return NormMaxLevel;             
}

