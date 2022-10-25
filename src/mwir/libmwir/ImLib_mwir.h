/* 
 * Filename : ImLib.h
 * 
 * Class Description
 * 1. Utils : useful functions
 */
#ifndef _IMLIB_H
#define _IMLIB_H
#include "Array.h"
#include "cdflib.h"

template <class DATATYPE>
class Utils
{
 protected:
  static double dist2D(double x, double y, double x0, double y0)
    {
      return sqrt((x - x0) * (x - x0)+(y - y0) * (y - y0));
    }

 public:
  // median value of the absolute value of a given data set
  static DATATYPE absMedian (to_array<DATATYPE, true> &data);
  // critical threshold for standard normal distribution, tailProba <= 0.5
  static double criticalThreshGauss (double tailProba);
  // cumulative of the normal distribution
  static double cumNormal (double x);
  // Anscombe transform of a data set
  static void anscombeTransform (to_array<DATATYPE, true> &data);
  static void invAnscombeTransform (to_array<DATATYPE, true> &data);

  // generate Donoho's clipped block 1D signal
  static void blockSignal (int len, double scaling, to_array<DATATYPE, true> &signal1D);
  // generate a 1D signal with a Gaussian form
  static void regularSignal (int cx, double sigma, double peakIntense, double backIntense, to_array<DATATYPE, true> &signal1D);
  // generate a 3D signal with a regional evalutionary Gaussian form
  static void regularSignal (int cx, int cy, double sigma0, double sigma1, double k, \
			     double peakIntense, double backIntense, to_array<DATATYPE, true> &signal3D, double flow[]);

  // calculate the source flow of a Gaussian form signal
  static void gaussianSourceFlow (to_array<DATATYPE, true> &signal3D, double cx, double cy, \
				  double sigma0, double sigma1, double k, double backIntense, double flow[]);

  // mean of a vector
  static double mean (DATATYPE vec[], int len);
  // standard deviation of a vector
  static double std (DATATYPE vec[], int len);
  // Normalized Integrated Square Error where the noise is Poissonian
  static double NISE (to_array<DATATYPE, true> &orig, to_array<DATATYPE, true> &estim);
  // Integrated Square Error where the noise has a spatially stable variance
  static double ISE (to_array<DATATYPE, true> &orig, to_array<DATATYPE, true> &estim);

  // return the combinatorial number C(n, k)
  static double Cnk (int n, int k);
  // discrete Binomial distribution C(n, k) p^k (1-p)^{n-k}
  static double binomial (int n, int k, double p);
  // EM 2-Mixture Binomial Parameters estimation
  // alpha is the mixture H0 portion (MUST be initialized)
  // pa    is true detection probability (MUST be initialized)
  // pi    is false detection probability (MUST be initialized)
  // NExp  is the total number of experiences
  // ndet  record for each pixel, the number of detection among NExp experiences
  // ndlen is the length of the array ndet
  // NITER is max. number of iteration
  static void EMMixBinomialParam (double &alpha, double &pa, double &pi, int NExp, to_array<DATATYPE, true> ndet[], int ndlen, int NITER);
//  // ROC area where [alpha0, alpha1] is the valid H0 portion interval
//  static double rocArea (int len, double alpha[], double pa[], double pi[], double alpha0, double alpha1);
  // ROC area
  static double rocArea (int len, double pa[], double pi[]);
};

template <class DATATYPE> 
DATATYPE Utils<DATATYPE>::absMedian (to_array<DATATYPE, true> &data)
{
	DATATYPE medv;
	int less, greater, equal;
	double min, max, guess;
	DATATYPE v, maxltguess, mingtguess;
	bool alldone = false;
	int dl = data.n_elem();
	
	min = ABS(data(0)); max = min;
	for (int i=0; i<dl; i++)
	{
		v = ABS(data(i));
		if (v<min) min = v;
		else if (v>max) max = v;
	}
	
	while (!alldone)
	{
		guess = (min+max)/2;
		less = 0;
		greater = 0;
		equal = 0;
		maxltguess = min;
		mingtguess = max;

		for (int i=0; i<dl; i++)
		{
	  	    v = ABS(data(i));
			if (v<guess)
			{
				less++;
				if (v>maxltguess) maxltguess = v;
			} 
			else if (v>guess)
			{
				greater++;
				if (v<mingtguess) mingtguess = v;
			} 
			else equal++;
		}
				
		if ((less <= (dl+1)/2) && (greater <= (dl+1)/2)) alldone = true;
		else if (less>greater) max = maxltguess;
		else min = mingtguess;
	}
	
    if (less >= (dl+1)/2) medv = maxltguess;
	else if (less+equal >= (dl+1)/2) medv = guess;
	else medv = mingtguess;
		
	return medv;	
}

template <class DATATYPE> 
double Utils<DATATYPE>::criticalThreshGauss (double tailProba)
{
      double p = MIN(1., MAX(1 - tailProba, .5));
      return MAX(0, stvaln(&p));  
}

template <class DATATYPE> 
double Utils<DATATYPE>::cumNormal (double x)
{
  double arg = x, result, ccum;
  cumnor(&arg, &result, &ccum);
  return result;
}

template <class DATATYPE> 
void Utils<DATATYPE>::anscombeTransform (to_array<DATATYPE, true> &data)
{
  int nlen = data.n_elem();
  for (int i=0; i<nlen; i++)
    data(i) = (DATATYPE) (2. * sqrt(data(i) + 3./8.));
}

template <class DATATYPE> 
void Utils<DATATYPE>::invAnscombeTransform (to_array<DATATYPE, true> &data)
{
  int nlen = data.n_elem();
  for (int i=0; i<nlen; i++)
    data(i) = (DATATYPE) ((data(i) * data(i)) / 4. - 3./8.);
}

template <class DATATYPE> 
void Utils<DATATYPE>::blockSignal (int len, double scaling, to_array<DATATYPE, true> &signal1D)
{
  signal1D.resize(len);

  for (int x=0; x<len; x++)
    {
      if      (x>=0.098*len && x<0.127*len) signal1D(x) = (DATATYPE) (15. * scaling);
      else if (x>=0.146*len && x<0.244*len) signal1D(x) = (DATATYPE) (9.  * scaling);
      else if (x>=0.264*len && x<0.390*len) signal1D(x) = (DATATYPE) (12. * scaling);
      else if (x>=0.439*len && x<0.586*len) signal1D(x) = (DATATYPE) (6.  * scaling);
      else if (x>=0.586*len && x<0.683*len) signal1D(x) = (DATATYPE) (18. * scaling);
      else if (x>=0.683*len && x<0.703*len) signal1D(x) = (DATATYPE) (8.  * scaling);
      else if (x>=0.703*len && x<0.732*len) signal1D(x) = (DATATYPE) (15. * scaling); 
      else                                  signal1D(x) = (DATATYPE) (3.  * scaling); 
    }  
}

template <class DATATYPE>
void Utils<DATATYPE>::regularSignal (int cx, double sigma, double peakIntense, double backIntense, to_array<DATATYPE, true> &signal1D)
{
  int nx = signal1D.nx();
  double sig = 2. * sigma * sigma;

  for (int x=0; x<nx; x++)
    signal1D(x) = (DATATYPE) (peakIntense * exp(-(x-cx)*(x-cx)/sig) + backIntense);
}

template <class DATATYPE> 
void Utils<DATATYPE>::regularSignal (int cx, int cy, double sigma0, double sigma1, double k, \
				     double peakIntense, double backIntense, to_array<DATATYPE, true> &signal3D, double flow[])
{
  int nx = signal3D.nx(), ny = signal3D.ny(), nz = signal3D.nz();
  double dsig = (sigma1 - sigma0) / (nz - 1), sig = sigma0;
  double vv;
  double tz = 2. * nz * nz / 16.;
  
  signal3D.resize(nx, ny, nz);
  for (int z=0; z<nz; z++)
    {
      double fl = 0;
      double txy = 2. * sig * sig;
      for (int x=0; x<nx; x++)
	for (int y=0; y<ny; y++)
	  {
	    if (dist2D(x,y,cx,cy) <= k*sig)
	      {
		vv = peakIntense * exp(-(x-cx)*(x-cx)/txy - (y-cy)*(y-cy)/txy - z*z/tz);
		fl += vv;
		vv += backIntense;
	      }
	    else vv = backIntense;
	    signal3D(x,y,z) = (DATATYPE) vv;
	  }
      flow[z] = ((DATATYPE) fl);
      sig += dsig;
    }
}

template <class DATATYPE>
void Utils<DATATYPE>::gaussianSourceFlow (to_array<DATATYPE, true> &signal3D, double cx, double cy, \
					  double sigma0, double sigma1, double k, double backIntense, double flow[])
{
  int nx = signal3D.nx(), ny = signal3D.ny(), nz = signal3D.nz();
  double dsig = (sigma1 - sigma0) / (nz - 1), sig = sigma0;

  for (int z=0; z<nz; z++)
    {
      flow[z] = 0;
      for (int x=0; x<nx; x++)
	for (int y=0; y<ny; y++)
	  if (dist2D(x,y,cx,cy) <= k*sig)
	    flow[z] += signal3D(x,y,z) - backIntense;
      sig += dsig;
    }
}

template <class DATATYPE> 
double Utils<DATATYPE>::mean (DATATYPE vec[], int len)
{
  double sum = 0.;

  for (int i=0; i<len; i++)
    sum += vec[i];

  return sum / len;
}

template <class DATATYPE> 
double Utils<DATATYPE>::std (DATATYPE vec[], int len)
{
  if (len == 1) return 0.;

  double M = Utils<DATATYPE>::mean(vec, len);
  double sum = 0.;

  for (int i=0; i<len; i++)
    sum += (vec[i] - M) * (vec[i] - M);

  return sqrt(sum / (len-1));
}

template <class DATATYPE> 
double Utils<DATATYPE>::NISE (to_array<DATATYPE, true> &orig, to_array<DATATYPE, true> &estim)
{
  int dlen = orig.n_elem();
  double sum = 0., lambda, temp;

  for (int i=0; i<dlen; i++)
    {
      lambda = orig(i);
      temp = estim(i) - lambda;
      sum += temp * temp / lambda;
    }

  return sum / dlen;
}

template <class DATATYPE> 
double Utils<DATATYPE>::ISE (to_array<DATATYPE, true> &orig, to_array<DATATYPE, true> &estim)
{
  int dlen = orig.n_elem();
  double sum = 0., temp;

  for (int i=0; i<dlen; i++)
    {
      temp = estim(i) - orig(i);
      sum += temp * temp;
    }

  return sum / dlen;
}

template <class DATATYPE> 
double Utils<DATATYPE>::Cnk (int n, int k)
{
  double cnk;
  int i;
  int mn;
  int mx;

  mn = MIN(k, n-k);

  if (mn < 0)
    cnk = 0;
  else if (mn == 0)
    cnk = 1;
  else
  {
    mx = MAX (k, n-k);
    cnk = mx + 1;
    
    for (i=2; i<=mn; i++)
      cnk = (cnk * ( mx + i )) / ((double) i);
  }

  return cnk;
}

template <class DATATYPE> 
double Utils<DATATYPE>::binomial (int n, int k, double p)
{
  return Utils<DATATYPE>::Cnk(n, k) * pow(p, (double)k) * pow(1-p, (double)(n-k)); 
}

template <class DATATYPE> 
void Utils<DATATYPE>::EMMixBinomialParam (double &alpha, double &pa, double &pi, int NExp, to_array<DATATYPE, true> ndet[], int ndlen, int NITER)
{
  long N = 0;
  int k;
  double pI, pA, t1, t2, A, B, C, D;

  for (int i=0; i<ndlen; i++)
    N += ndet[i].n_elem();

  for (int i=0; i<NITER; i++)
    {
      A = B = C = D = 0.;
      for (int l=0; l<ndlen; l++)
	{
	  int datalen = ndet[l].n_elem();
	  for (int j=0; j<datalen; j++)
	    {
	      k = (int) ndet[l](j);
	      t1 = alpha * binomial(NExp, k, pi);
	      t2 = (1-alpha) * binomial(NExp, k, pa);
	      if (t1 == 0) pI = 0.;
	      else pI = t1 / (t1 + t2);
	      pA = 1. - pI;
	      A += pI;
	      B += pA;
	      C += k * pI;
	      D += k * pA;
	    }
	}
      alpha = A / N;
      if (B == 0) pa = 0.;
      else pa = D / (NExp * B);
      if (A == 0) pi = 0.;
      else pi = C / (NExp * A);
    }
}

//template <class DATATYPE> 
//double Utils<DATATYPE>::rocArea (int len, double alpha[], double pa[], double pi[], double alpha0, double alpha1)
//{
//  int pre = -1;
//  double area = 0.;
//  
//  for (int i=0; i<len; i++)
//    if ((alpha[i] >= alpha0) && (alpha[i] <= alpha1))
//      {
//	if (pre < 0)
//	  {
//	    area = pi[0] * pa[0] / 2.;
//	    pre = i;
//	  }
//	else if (pi[i] > pi[pre])
//	  {
//	    area += (pa[pre] + pa[i]) * (pi[i] - pi[pre]) / 2.;
//	    pre = i;
//	  }
//      }
//  if (pi[pre] != 1)
//    area += (pa[pre] + 1) * (1 - pi[pre]) / 2.;
//  
//  return area;
//}

template <class DATATYPE> 
double Utils<DATATYPE>::rocArea (int len, double pa[], double pi[])
{
  int pre = -1;
  double area = 0.;
  
  for (int i=0; i<len; i++)
  {
	if (pre < 0)
	  {
	    area = pi[0] * pa[0] / 2.;
	    pre = i;
	  }
	else if (pi[i] > pi[pre])
	  {
	    area += (pa[pre] + pa[i]) * (pi[i] - pi[pre]) / 2.;
	    pre = i;
	  }
  }
  if (pi[pre] != 1)
    area += (pa[pre] + 1) * (1 - pi[pre]) / 2.;
  
  return area;
}

#endif
