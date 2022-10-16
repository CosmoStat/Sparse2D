/* 
 * Filename : ImLib.h
 * 
 * Class Description
 * 1. Utils : useful functions
 */
#ifndef _IMLIB_H
#define _IMLIB_H
#include "DefMath.h"
#include "Array.h"
#include "cdflib.h"
#include "Wavelet.h"

static double t1[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static double t2[] = {1, 0.2734375, 0.12347412109375, 0.0603594779968262, 0.0300147198140621, \
                   0.014986945054261, 0.00749092722048772, 0.00374514565088013, \
                   0.00187253308688793, 0.00093626157632392, 0.000468130167278172 \
                  };
static double t3[] = {1, 0.08447265625, 0.017338752746582, 0.00414866860955954, 0.00102616434651281, \
                   0.00025586256661736, 6.39233749606246e-5, 1.59782042631771e-5, \
                   3.99438613268932e-6, 9.9858622538761e-7, 2.49645912118706e-7 \
                  };
static double sumhh[] = {
         0, // dummy number
         0.375,
             0.15771484375,
        0.0760841369628906,
        0.0377178490161896,
        0.0188190417829901,
       0.00940455525596917,
       0.00470165753587537,
       0.00235075127553241,
       0.00117536595181247,
      0.000587681765180673,
      0.000293840731250224
};

// Antonini79 L2-normalized
static double t1Antonini79[] = {1., 1.41421356237310, 2.0, 2.82842712474619, \
                                4.0, 5.65685424949239, 8.0, 11.31370849898477, \
                                16.0, 22.62741699796953, 32.0};
static double t2Antonini79[] = {1., 0.98295365728765, 1.03060246849226, 1.05209302224409, \
                                1.05847325456386, 1.06015396256899, 1.06058070717996, \
                                1.06068786810321, 1.06071469140228, 1.06072139947288, \
                                1.06072307664033};
static double t3Antonini79[] = {1., 0.63570269427585, 0.47944255910318, 0.35015620568978, \
                                0.25001964565880, 0.17724608900296, 0.12541390863709, \
                                0.08869559519435, 0.06271983610386, 0.04435007755512, \
                                0.03136032122693};
static double sigma_g_Antonini79[] = {0., 1.02001762915889, 0.99196230954764, 1.00436911704512, \
                                      1.01088444846542, 1.01284519436391, 1.01336351749413, \
                                      1.01349529471499};
// hg                                      
static double sigma_hg_Antonini79[] ={0., 1.01128647562687, 1.01572054915842, \
                                      1.01478685097925, 1.01394498056575, 1.01364900390853, \
                                      1.01356745198182, 1.01354649513165};  
// gh
static double sigma_gh_Antonini79[] ={0., 1.01128647562687, 1.01572054915842, \
                                      1.01478685097925, 1.01394498056575, 1.01364900390853, \
                                      1.01356745198182, 1.01354649513165};   
// gg
static double sigma_gg_Antonini79[] ={0., 1.04043596379492, 0.98398922356310, \
                                      1.00875732327399, 1.02188736814924, 1.02585538774607, \
                                      1.02690561858808, 1.02717271240943};
             
// Antonini79                                      
static double Hht2[] = {0.98295365728765, \
   1.03060246849226, \
   1.05209302224409, \
   1.05847325456386, \
   1.06015396256899, \
   1.06058070717996, \
   1.06068786810321, \
   1.06071469140228, \
   1.06072139947288, \
   1.06072307664033};

static double Hgt2[] = {1.04043596379492, \
   0.96721580603298, \
   1.03962778747582, \
   1.07512056954920, \
   1.08584049097946, \
   1.08867806073051, \
   1.08939956172315, \
   1.08958081336585, \
   1.08962618736738, \
   1.08963753501908};
            
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
  // Inverse Anscombe transform of a data set
  // bias correction is recommended while reconstruct the data after estimation
  static void invAnscombeTransform (to_array<DATATYPE, true> &data, bool cbias);

  // generate Donoho's clipped block 1D signal
  static void blockSignal (int len, double scaling, to_array<DATATYPE, true> &signal1D);
  // generate a 1D signal with a Gaussian form
  static void regularSignal (int cx, double sigma, double peakIntense, double backIntense, to_array<DATATYPE, true> &signal1D);
  // generate a 3D signal with a regional evolutionary Gaussian form
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
  // ROC area where [alpha0, alpha1] is the valid H0 portion interval
  // static double rocArea (int len, double alpha[], double pa[], double pi[], double alpha0, double alpha1);
  // ROC area
  static double rocArea (int len, double pa[], double pi[]);

  // b3Spline VST coefs
  static void b3VSTCoefs (int dim, int scale, double &b, double &c);  
  // b3Spline VST coefs 2D+1D
  static void b32D1DVSTCoefs (int scalexy, int scalez, double &b, double &c);  
  // Antonini 7/9 VST coefs
  static void antonini79VSTCoefs (int dim, int scale, double &b, double &c);

  // b3Spline coef Var
  static double b3CoefVar (int dim, int scale);

  // Coupled schema
  // b = 1
  static double b3VSTCoefVar (int dim, int scale);  
  // corrected mean and var b = 1
  static double b3VSTWaveCoefCMean (int dim, int scale, double lambda);
  static double b3VSTWaveCoefCVar (int dim, int scale, double lambda);  

  // Separated schema
  // B3Spline wavelet coefficient var
  static double b3VSTSepCoefVar (int dim, int scale);
  // Antonini 7/9 wavelet coef Var
  static double antonini79SepCoefStd (int scale);
  static void antonini79SepCoefStd (int scale, double &stdhg, double &stdgh, double &stdgg);
  // Antonini 7/9 wavelet VST coef Var
  static double antonini79VSTSepCoefStd (int scale);
  static void antonini79VSTSepCoefStd (int scale, double &stdhg, double &stdgh, double &stdgg);

  // 2D+1D
  // variance of two stabilized 1D band : T(h^{scalexy}*h^{scalez-1}c0) - T(h^{scalexy}*h^{scalez}c0)
  // b = 1
  static double b32D1DVSTCoefZVar (int dim, int scalexy, int scalez);
  // variance of two stabilized 2D band
  // b = 1
  static double b32D1DVSTCoefXYVar (int dim, int scalexy, int scalez);

  // variance of direct wavelet coef
  static double b3DirWaveCoefVar (int dim, int scale, double lambda);
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
void Utils<DATATYPE>::invAnscombeTransform (to_array<DATATYPE, true> &data, bool cbias)
{
  double cb = cbias ? 1. : 0.;  
  int nlen = data.n_elem();
  for (int i=0; i<nlen; i++)
    data(i) = (DATATYPE) ((data(i) * data(i) + cb) / 4. - 3./8.);
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

template <class DATATYPE> 
void Utils<DATATYPE>::b3VSTCoefs (int dim, int scale, double &b, double &c)
{    
    if (scale > 10) throw DataException("Utils::b3VSTCoefs: VST coefs not precalculated for scale > 10");
    
    double tt1 = pow(t1[scale], dim);
    double tt2 = pow(t2[scale], dim);
    double tt3 = pow(t3[scale], dim);
    
    c = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    b = 2 * sqrt(tt1/tt2);
}

template <class DATATYPE> 
void Utils<DATATYPE>::antonini79VSTCoefs (int dim, int scale, double &b, double &c)
{    
    if (scale > 10) throw DataException("Utils::antonini79VSTCoefs: VST coefs not precalculated for scale > 10");
    
    double tt1 = pow(t1Antonini79[scale], dim);
    double tt2 = pow(t2Antonini79[scale], dim);
    double tt3 = pow(t3Antonini79[scale], dim);
    
    c = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    b = 2 * sqrt(ABS(tt1)/tt2);
}

template <class DATATYPE> 
double Utils<DATATYPE>::antonini79SepCoefStd (int scale)
{
     return sqrt(Hgt2[scale-1]);
}

template <class DATATYPE> 
void Utils<DATATYPE>::antonini79SepCoefStd (int scale, double &stdhg, double &stdgh, double &stdgg)
{
     stdhg = sqrt(Hht2[scale-1]*Hgt2[scale-1]);
     stdgh = sqrt(Hgt2[scale-1]*Hht2[scale-1]);
     stdgg = Hgt2[scale-1];
}

template <class DATATYPE> 
double Utils<DATATYPE>::antonini79VSTSepCoefStd (int scale)
{
     return sigma_g_Antonini79[scale];
}

template <class DATATYPE> 
void Utils<DATATYPE>::antonini79VSTSepCoefStd (int scale, double &stdhg, double &stdgh, double &stdgg)
{
     stdhg = sigma_hg_Antonini79[scale];
     stdgh = sigma_gh_Antonini79[scale];
     stdgg = sigma_gg_Antonini79[scale];    
}

template <class DATATYPE> 
void Utils<DATATYPE>::b32D1DVSTCoefs (int scalexy, int scalez, double &b, double &c)
{
    if ((scalexy > 10) || (scalez > 10)) throw DataException("Utils::b32D1DVSTCoefs: VST coefs not precalculated for scale > 10");
    
    double tt1 = pow(t1[scalexy], 2) * t1[scalez];
    double tt2 = pow(t2[scalexy], 2) * t2[scalez];
    double tt3 = pow(t3[scalexy], 2) * t3[scalez];
    
    c = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    b = 2 * sqrt(tt1/tt2);    
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3VSTCoefVar (int dim, int scale)
{
    if (scale > 10) throw DataException("Utils::b3VSTCoefVar: VST coefs Var not precalculated for scale > 10");

    double tt11 = pow(t1[scale-1], dim);
    double tt1  = pow(t1[scale], dim);
    double tt21 = pow(t2[scale-1], dim);
    double tt2  = pow(t2[scale], dim);
    double shh  = pow(sumhh[scale], dim);
    
    return MAX((tt21 / (4. * tt11 * tt11) + tt2 / (4. * tt1 * tt1) - shh / (2. * tt11 * tt1)), 0.); // in case of numerical error
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3VSTSepCoefVar (int dim, int scale)
{
    double sgghh = pow(t2[1],dim)+1-2*pow(3./8., dim); 
    // correlation on h is ignored : not exact since g is not separable
    
    return MAX(sgghh, 0.); // in case of numerical error
}

template <class DATATYPE> 
double Utils<DATATYPE>::b32D1DVSTCoefZVar (int dim, int scalexy, int scalez)
{   // tau_1 = 1 in B3 low-pass filter
    double tt21 = pow(t2[scalexy], 2) * t2[scalez-1];
    double tt22 = pow(t2[scalexy], 2) * t2[scalez];
    double shh = pow(t2[scalexy], 2) * sumhh[scalez];
    
    return MAX(tt21 / 4. + tt22 / 4. - shh / 2., 0.); // in case of numerical error
}

template <class DATATYPE> 
double Utils<DATATYPE>::b32D1DVSTCoefXYVar (int dim, int scalexy, int scalez)
{   // tau_1 = 1 in B3 low-pass filter
    double tt21 = pow(t2[scalexy-1], 2) * t2[scalez];
    double tt22 = pow(t2[scalexy], 2) * t2[scalez];
    double shh = pow(sumhh[scalexy], 2) * t2[scalez];
    
    return MAX(tt21 / 4. + tt22 / 4. - shh / 2., 0.); // in case of numerical error
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3VSTWaveCoefCMean (int dim, int scale, double lambda)
{
    double cpre, ccur;
    
    double tt1 = pow(t1[scale-1], dim);
    double tt2 = pow(t2[scale-1], dim);
    double tt3 = pow(t3[scale-1], dim);
    
    cpre = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    
    tt1 = pow(t1[scale], dim);
    tt2 = pow(t2[scale], dim);
    tt3 = pow(t3[scale], dim);
    
    ccur = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    
    double tt2pre = pow(t2[scale-1], dim);
    double tt2cur = pow(t2[scale], dim);

    double res = sqrt(lambda+cpre) - sqrt(lambda+ccur) - lambda*(tt2pre/(8.*pow(lambda+cpre, 1.5)) - tt2cur/(8.*pow(lambda+ccur, 1.5)));
    
    return res;
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3VSTWaveCoefCVar (int dim, int scale, double lambda)
{
    double cpre, ccur, res;
    
    double tt1 = pow(t1[scale-1], dim);
    double tt2 = pow(t2[scale-1], dim);
    double tt3 = pow(t3[scale-1], dim);
    
    cpre = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    
    tt1 = pow(t1[scale], dim);
    tt2 = pow(t2[scale], dim);
    tt3 = pow(t3[scale], dim);
    
    ccur = 7*tt2/(8*tt1) - tt3 / (2*tt2);
    
    double tt2pre = pow(t2[scale-1], dim);
    double tt2cur = pow(t2[scale], dim);
    double shh  = pow(sumhh[scale], dim);

    res = lambda*tt2pre/(4.*(lambda+cpre)) + lambda*tt2cur/(4.*(lambda+ccur));
    res -= .5*lambda*shh/(sqrt(lambda+cpre) * sqrt(lambda+ccur));
    res = MAX(res, 0.); // in case of numerical error
    
    return res;
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3CoefVar (int dim, int scale)
{
    double tt2pre = pow(t2[scale-1], dim);
    double tt2cur = pow(t2[scale], dim);
    double shh  = pow(sumhh[scale], dim);    
    double res = tt2pre + tt2cur - 2 * shh;
    
    return  MAX(res, 0);
}

template <class DATATYPE> 
double Utils<DATATYPE>::b3DirWaveCoefVar (int dim, int scale, double lambda)
{
    if (scale > 10) throw DataException("Utils::b3VSTCoefVar: VST coefs Var not precalculated for scale > 10");

    double tt21 = pow(t2[scale-1], dim);
    double tt2  = pow(t2[scale], dim);
    double shh  = pow(sumhh[scale], dim);
    
    return MAX(lambda*(tt21+tt2-2.*shh),0.); // in case of numerical error
}

#endif
