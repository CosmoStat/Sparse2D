

#include <cstdio>
#include <iostream>
#include <cmath>
#include "GlobalInc.h"

// static int DebugAlloc;
#include "MR_Calloc.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Allocator::Allocator ()
{
  precision = NULL;
}

/*---------------------------------------------------------------------------*/

Allocator::~Allocator ()
{
  if (precision != NULL)
    delete [] precision;
}

/*---------------------------------------------------------------------------*/

void Allocator::resetPrecision (int nSets)
{
  if (precision != NULL)
    delete [] precision;

  precision = new int [nSets];
  for (int i = 0; i < nSets; i++)
    precision[i] = 0;
}

/*---------------------------------------------------------------------------*/

void Allocator::optimalAllocate (BandRateDist &BRD, int nSets, 
				 int budget, int augment, float *weight)
{
  float bitBudget = 8*budget;

  float lambda, lambdaLow, lambdaHigh;
  float rateLow, rateHigh, currentRate;
  float distLow, distHigh, currentDist;

  resetPrecision (nSets);
  
  lambdaLow = 0.0;
  allocateLambda (BRD, nSets, lambdaLow, rateLow, distLow, weight);
  if (rateLow < bitBudget) // this uses the largest possible # of bits
  {
    //if (DebugAlloc) cout << "budget ok: " << "bitBudget = " << bitBudget << 
    //                      "rateLow = " << rateLow << endl;
    return;                //   -- if this is within the budget, do it
  }
  lambdaHigh = 1000000.0;
  float lastRateHigh = -1;
  do {
    // try to use the smallest possible # of bits
    allocateLambda (BRD, nSets, lambdaHigh, rateHigh, distHigh, weight);

    // if this is still > bitBudget, try again w/ larger lambda
    if (rateHigh > bitBudget && lastRateHigh != rateHigh) {
      lambdaLow = lambdaHigh;
      rateLow = rateHigh;
      distLow = distHigh;
      lambdaHigh *= 10.0;
    }
  } while (rateHigh > bitBudget && lastRateHigh != rateHigh);

  // give up when changing lambda has no effect on things
  if (lastRateHigh == rateHigh)
    return;

  // Note rateLow will be > rateHigh
  if (rateLow < bitBudget) 
    fprintf (stderr, "Failed to bracket bit budget = %d: rateLow = %g rateHigh = %g\n", 
	   budget, rateLow, rateHigh);
  
  while (lambdaHigh - lambdaLow > 0.01)  {
    lambda = (lambdaLow + lambdaHigh)/2.0;
    
    allocateLambda (BRD, nSets, lambda, currentRate, currentDist, weight);
    
    if (currentRate > bitBudget)
      lambdaLow = lambda;
    else
      lambdaHigh = lambda;
  }
  
  if (currentRate > bitBudget)  {
    lambda = lambdaHigh;
    allocateLambda (BRD, nSets, lambda, currentRate, currentDist, weight);
  }

  if (augment)
    greedyAugment (BRD, nSets, bitBudget-currentRate, weight);
}

/*---------------------------------------------------------------------------*/

void Allocator::allocateLambda  (BandRateDist &BRD, int nSets, 
				 float lambda, float &optimalRate, 
				 float &optimalDist, float *weight)
{
   int i, j;
   float G, minG, minRate, minDist;

   optimalRate = optimalDist = 0.0;

   // want to minimize G = distortion + lambda * rate
   
   // loop through all rate-distortion curves
   for (i = 0; i < nSets; i++)  {
     minG = minRate = minDist = Maxfloat;
     
     for (j = 0; j < BRD.NbrQuant; j++) {
       G = weight[i]*BRD.Dist(i,j) + lambda * BRD.Rate(i,j);
       if (G < minG)  {
	 minG = G;
	 minRate =  BRD.Rate(i,j);
	 minDist = weight[i]* BRD.Dist(i,j);
	 precision[i] = j;
       }
     }
     
     optimalRate += minRate;
     optimalDist += minDist;
   }
   /*
   printf ("lambda = %g  optimal rate = %g, optimal dist = %g\n",
   	    lambda, optimalRate, optimalDist); 
   for (i = 0; i < nSets; i++) printf (" Prec(%d) = %d\n", i+1, precision[i]);
   */
}

/*---------------------------------------------------------------------------*/

void Allocator::greedyAugment (BandRateDist &BRD, int nSets, 
			       float bitsLeft, float *weight)
{
  int bestSet, newPrecision = -1;
  float delta, maxDelta, bestDeltaDist, bestDeltaRate = 0;

  do {
    bestSet = -1;
    maxDelta = 0;
	    
    // Find best coeff set to augment 
    for (int i = 0; i < nSets; i++) {
      for (int j = precision[i]+1; j < BRD.NbrQuant; j++) {
	float deltaRate =  BRD.Rate(i,j) -BRD.Rate(i, precision[i]);
 	float deltaDist = -weight[i]*( BRD.Dist(i,j) - BRD.Dist(i,precision[i]));
 	
	if (deltaRate != 0 && deltaRate <= bitsLeft) {
	  delta = deltaDist / deltaRate;
	  
	  if (delta > maxDelta) {
	    maxDelta = delta;
	    bestDeltaRate = deltaRate;
	    bestDeltaDist = deltaDist;
	    bestSet = i;
	    newPrecision = j;
	  }
	}
      }
    }
    
    if (bestSet != -1) {
      precision[bestSet] = newPrecision;
      bitsLeft -= bestDeltaRate;
    }
  } while (bestSet != -1);
}

/*---------------------------------------------------------------------------*/

void Allocator::print (BandRateDist &BRD, int nSets, int *NdataBand)
{
  float totalRate = 0, totalDist = 0;
  int totalData = 0;

  printf ("Set  Precision     Rate                Distortion\n");
  for (int i = 0; i < nSets; i++) {
    printf ("%2d:     %2d       %9.2f   %5.2f   %11.2f   %7.2f\n", i, precision[i],
	    BRD.Rate(i, precision[i]), 
	    BRD.Rate(i, precision[i]) / (float) NdataBand[i], 
	    BRD.Dist(i, precision[i]),   
	    BRD.Dist(i, precision[i])/(float) NdataBand[i] );
    totalRate += BRD.Rate(i, precision[i]);
    totalDist += BRD.Dist(i, precision[i]);
    totalData +=  NdataBand[i];
  }
  float rms = sqrt(totalDist/(float)totalData);
  float psnr = 20.0 * log(255.0/rms)/log(10.0);

  printf ("\n");
  printf ("total rate = %g\n", totalRate/8.0);
  printf ("total dist = %g\n", totalDist);
  printf ("total coeffs = %d\n", totalData);
  printf ("RMS error = %g\n", rms);
  printf ("PSNR (transform domain) = %g\n", psnr);
}

/*---------------------------------------------------------------------------*/
