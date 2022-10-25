/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre, Hubert Druesne
**
**    Date:  6/08/03
**    
**    File:  IM1D_ALDCT.cc
**
**    Modification history :
**
******************************************************************************/

// #include <math.h>
// #include "IM1D_ALDCT.h"
// #include "IM_Obj.h"
// #include "IM_IO.h"
#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM1D_Block.h"
#include "IM1D_Dct.h"
#include "IM1D_ALDCT.h"

/******************************************************************************/

void BASIS::alloc(int NbrBand){
    int Nb=1;
    for (int i=0; i<NbrBand; i++) Nb*=2;
    levels.alloc(Nb);
    contents = new fltarray [Nb];
    indice=0;
}

/******************************************************************************/

void ALDCT::alloc (fltarray* & Tab, int Nx, int NbrBand) {
   
    char ch[80];
    NbScale = NbrBand;
    Tab = new fltarray [NbScale];
    int N=Nx;
    int i=1;
    if (!(is_power_of_2(Nx))){
       while (i<Nx) i*=2;
       N=i;
    }
    NbPts = N;
    for (int s=0; s<NbScale; s++) {
	sprintf (ch, "ALDCT_%d", s+1);
        Tab[s].alloc (N, ch);   
    }
    base.alloc(NbScale);

}

/******************************************************************************/

void ALDCT::transform(fltarray& Signal, fltarray* TabALDCT){
   
   fltarray ReformSig(NbPts);
   for (int i=0; i<Signal.nx(); i++) ReformSig(i)=Signal(i);
   
   int x=1;
   for (int s=0; s<NbScale; s++) {
      int BlockSize=NbPts/(x);
      AL_ldct.alloc(NbPts, BlockSize, False, False);
      AL_ldct.transform(ReformSig);
      TabALDCT[s] = AL_ldct._DCTSig;
      x*=2;
   }

}
/******************************************************************************/

void ALDCT::recons(fltarray& SigRec){
   fltarray Sig(NbPts);
   int N=0;
   int indice=0;
   for (int j=0; j<base.indice; j++){
      N=base.contents[j].nx();
      Sig.resize(N);
      sig_inv(base.contents[j], Sig, N,
                             N, False, False);
      for (int i=0; i<N; i++) SigRec(indice++)=Sig(i);
   }
   base.indice=0;
}
/******************************************************************************/
void ALDCT::threshold (fltarray* TabALDCT, float NSigma, int IterNumber) {

   float Noise = SigmaNoise * norm();
   float Level = NSigma * SigmaNoise * Sensibility;
   
   for (int s=0; s<NbScale; s++){
      for (int i=0; i<NbPts; i++)
         TabALDCT[s](i) = update (TabALDCT[s](i), Level, Noise);
   } 
} 


/****************************************************************************/
inline float ALDCT::update (float CoefSol, float Threshold, 
                                  float SoftLevel) {
					 
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




/******************************************************************************/

float ALDCT::cost_L2(fltarray& Sig){
  float absu;
  float cost=0;
  for (int i=0; i<Sig.nx(); i++){
     absu=fabs(Sig(i));
     if (absu>0) cost += exp(2*log(absu));
  }
  return cost;
}

/******************************************************************************/

float ALDCT::cost_L1(fltarray& Sig){

  float cost=0;
  for (int i=0; i<Sig.nx(); i++){
     cost+=fabs(Sig(i));
  }
  return cost;
}

/******************************************************************************/

float ALDCT::ml2logl2(fltarray& Sig){
  float cost=0;
  float usq;
  for (int i=0; i<Sig.nx(); i++){
     usq=Sig(i)*Sig(i);
     if (usq > 0) cost -= usq * log(usq);
  }
  return cost;
}

/******************************************************************************/

void ALDCT::getBestBasis (fltarray* TabALDCT, fltarray& BestBasis, 
                           int Infocost) {

   float Cost = compute_base (TabALDCT, 0, NbPts/2,0, 0, Infocost);
   
   BestBasis.alloc(NbPts);
   int new_indice=0;
   for (int j=0; j<base.indice; j++) {
   
      if (Write)
         cout <<"indice: "<< j <<" : " << base.contents[j].nx()
              <<"  level:"<<base.levels(j)<< endl;
           
      for (int i=0; i<base.contents[j].nx(); i++)
	 BestBasis(new_indice++)=base.contents[j](i);
   }   
   
   if (Write)
      fits_write_fltarr ("ALDCT_BestBasis", BestBasis);
   

}      

/******************************************************************************/


float ALDCT::getAbsMaxTransf (fltarray* TabALDCT) {

   float absMaxLevel=0;

   for (int s=0; s<NbScale-1; s++) {
   
      float MaxLevel = fabs((TabALDCT[s]).maxfabs());
      float NormMaxLevel = MaxLevel/SigmaNoise/norm()/Sensibility;
                   
      if (Write)
         cout << "Lambda ALDCT (" << s+1 << ") : " << NormMaxLevel 
              << ", Max amplitude : " << MaxLevel << endl;  
                                                              
      if (absMaxLevel <= NormMaxLevel) absMaxLevel = NormMaxLevel;
   }
   if (Verbose) cout << "Lambda ALDCT : " << absMaxLevel << endl;
            
   return absMaxLevel;
}


/******************************************************************************/

float ALDCT::compute_base (fltarray* parents, int ind_left, int ind_right,
                           int level, int Dir, int infocost) {
   int SizeBlock= ind_right - ind_left;
   fltarray parent(2*SizeBlock);
   fltarray left_child(SizeBlock);
   fltarray right_child(SizeBlock);
   float left_exist=0;
   float right_exist=0;
   float Cost;
   float Bestcost;
   float parent_cost;
   float left_child_cost;
   float right_child_cost;
   float comp;
   int ind;
      
//   cout<< "" << endl;   
//   cout<< "" << endl;  
//   if (Dir==0) cout << "level left= " << level << endl;
//   else cout << "level right= " << level << endl;
    
   if (level==NbScale-2) {
         for (int i=ind_left; i<ind_right; i++){
            parent(i-ind_left)=parents[level](i);
            parent(i+SizeBlock-ind_left)=parents[level](i+SizeBlock);
	 }
         if (infocost == 0)  parent_cost=ml2logl2(parent);
         if (infocost == 1) parent_cost=cost_L1(parent);
         base.levels(base.indice)=level;
         base.contents[base.indice].alloc(SizeBlock);
         base.contents[base.indice]=parent;
         base.indice+=1;
         Bestcost=parent_cost;
//         cout <<"indice1=" <<base.indice-1<<" level:"<<level<<endl;
//         cout <<"nb point de base enregistrée ==" << parent.nx() <<endl;
   }
   
   else{

      for (int i=ind_left; i<ind_right; i++){
         parent(i-ind_left)=parents[level](i);
         parent(i+SizeBlock-ind_left)=parents[level](i+SizeBlock);
      
         left_child(i-ind_left)=parents[level+1](i);
         right_child(i-ind_left)=parents[level+1](i+SizeBlock);
      
         left_exist+=left_child(i-ind_left);
         right_exist+=right_child(i-ind_left);
      }
   
      if (infocost ==0){
        parent_cost=ml2logl2(parent);
//	cout << "parent_cost = " << parent_cost << endl;
        left_child_cost=ml2logl2(left_child);
//        cout << "left_child_cost = " << left_child_cost << endl;
        right_child_cost=ml2logl2(right_child);
//        cout << "right_child_cost = " << right_child_cost << endl;
      }
   
      if (infocost == 1){
         parent_cost=cost_L1(parent);
//         cout << "parent_cost = " << parent_cost << endl;
	 left_child_cost=cost_L1(left_child);
//         cout << "left_child_cost = " << left_child_cost << endl;
         right_child_cost=cost_L1(right_child);
//         cout << "right_child_cost = " << right_child_cost << endl;
      }
   
      if (infocost == 2){
         parent_cost=cost_L2(parent);
//	 cout << "parent_cost = " << parent_cost << endl;
	 left_child_cost=cost_L2(left_child);
//         cout << "left_child_cost = " << left_child_cost << endl;
         right_child_cost=cost_L2(right_child);
//         cout << "right_child_cost = " << right_child_cost << endl;
      }
   
      ind = base.indice;
      Cost = 0;
      if (left_child.n_elem()!=0)
         Cost+=compute_base(parents, ind_left, ind_left+SizeBlock/2,
	                level+1, 0, infocost);
      if (right_child.n_elem()!=0)
         Cost+=compute_base(parents, ind_right, ind_right+SizeBlock/2,
	                level+1, 1, infocost);
//      cout << "Cost accu =====" << Cost <<endl;
//      cout << "Cost parent =====" << parent_cost <<endl;
      
      comp=Cost-parent_cost;
      
      if (comp<0){
         Bestcost=Cost;
      }
      else{
         Bestcost=parent_cost;
	 base.indice = ind;
//         cout <<"indice2=" <<base.indice<<" level:"<<level<<endl;
         base.levels(base.indice)=level;
         base.contents[base.indice].alloc(2*SizeBlock);
         base.contents[base.indice]=parent;
//         cout <<"nb point de base enregistrée ==" << parent.nx() <<endl;
         base.indice+=1;
      }
      
   }
   return Bestcost;
}

/******************************************************************************/

void sig_inv (fltarray& Trans, fltarray& Sig, int BlockSize, int Nx, 
                      Bool Overlap, Bool WeightFirst) {
                      
   Block1D B1D,B1DTrans;
   fltarray SigBlock, Block, TransBlock;
   
   B1D._Verbose = True;
   B1D._BlockOverlap = Overlap;
   B1D._WeightFirst = WeightFirst;
   B1D.alloc (Nx, BlockSize);
   int BS = B1D.block_size();
   int Nxt = B1D.nbr_block_nx() * BS;
   if (Trans.nx() != Nxt) {
       cout << "Error in sig_inv: Nxt = " <<  Nxt << endl;
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

/******************************************************************************/
