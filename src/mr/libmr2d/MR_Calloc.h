
/*---------------------------------------------------------------------------*/
// allocator.hh
//
// Given rate/distortion curves for nSets collections of transform
// coefficients (contained in CoeffSet objects), performs a
// constrained optimization of quantizer resolutions.  An array of
// quantizer precisions, precision[i], is found so that the sum (over
// i) of weight[i]*distortion[i][precision[i]] is minimized subject to
// the constraint that the sum (over i) of cost[i][precision[i]] is
// less than or equal to the given budget.
//
// Functions:
// ----------
// optimalAllocate     Does bit allocation using an algorithm described
//                     in Y. Shoham and A. Gersho, "Efficient bit
//                     allocation for an arbitrary set of quantizers,"
//                     IEEE Transactions on Acoustics, Speech, and
//                     Signal Processing, Vol. 36, No. 9,
//                     pp. 1445-1453, Sept 1988.
//
// greedyAugment       The Shoham & Gersho algorithm doesn't yield
//                     optimal allocations for all possible budgets.
//                     The optimalAllocate routine returns the best
//                     allocation that doesn't exceed the given
//                     budget.  GreedyAugment uses marginal analysis
//                     to greedily increase individual quantizer
//                     precisions until we reach the budget.
//                     Allocations will still be a little under budget
//                     but shouldn't be by much.  Note that the header
//                     info is not included in the overall budget.
//
// print               Prints out the current allocation
//
/*---------------------------------------------------------------------------*/

#ifndef _ALLOCATOR_
#define _ALLOCATOR_

/*---------------------------------------------------------------------------*/

#include "Array.h"
 
#define NBR_QUANT 10  // number of quantizers to examine for allocation
   

class BandRateDist {
public:
int NbrBand;
int NbrQuant;
fltarray Dist;
fltarray Rate;
fltarray QuantCoeff;
fltarray TabMin,TabMax,TabMean;

BandRateDist (int NBand, int NQuant)
{
   NbrBand = NBand;
   NbrQuant = NQuant;
   Dist.alloc(NbrBand,NbrQuant);
   Rate.alloc(NbrBand,NbrQuant);
   QuantCoeff.alloc(NbrBand,NbrQuant);
   TabMin.alloc(NbrBand);
   TabMax.alloc(NbrBand);
   TabMean.alloc(NbrBand);
}
~BandRateDist ()
{
   NbrBand = NbrQuant = 0;
   Dist.free();
   Rate.free();
   TabMin.free();
   TabMax.free();
   TabMean.free();
}
};


class Allocator {
public:
  Allocator ();
  ~Allocator ();

  void optimalAllocate (BandRateDist &BRD, int nSets, 
			int budget, int augment, float *weight);

  void greedyAugment   (BandRateDist &BRD, int nSets, float bitsLeft, 
			float *weight);

  void print           (BandRateDist &BRD, int nSets, int *NdataBand);

  void allocateLambda  (BandRateDist &BRD, int nSets, 
			float lambda, float &optimalRate, 
			float &optimalDist, float *weight);
  void resetPrecision  (int nSets);

  int nSets, *precision;
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
