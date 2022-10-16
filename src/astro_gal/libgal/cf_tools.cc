#include "cf_tools.h"

void gaussint(float z, float *p, float *q, double *pdf)
{
  //
  //	normal distribution probabilities accurate to 1.e-15.
  //	z = no. of standard deviations from the mean.
  //	p, q = probabilities to the left & right of z.   p + q = 1.
  //       pdf = the probability density.
  //
  //       based upon algorithm 5666 for the error function, from:
  //       hart, j.f. et al, 'computer approximations', wiley 1968
  //
  //       programmer: alan miller
  //
  //	latest revision - 30 march 1986
  //
  static double 
    p0=220.2068679123761,
    p1=221.2135961699311, 
    p2=112.0792914978709,
    p3=33.91286607838300, 
    p4=6.373962203531650,
    p5=.7003830644436881, 
    p6=.3526249659989109e-01;
  
  static double
    q0=440.4137358247522,
    q1=793.8265125199484, 
    q2=637.3336333788311,
    q3=296.5642487796737, 
    q4=86.78073220294608,
    q5=16.06417757920695, 
    q6=1.755667163182642,
    q7=.8838834764831844e-1;
  static double
    cutoff=7.071, root2pi=2.506628274631001;
  
  double zabs,expntl;

  zabs = fabs(z);
  //
  //	|z| > 37.
  //
  if (zabs > 37.0)
    {
      *pdf = 0.0;
      if (z > 0.0) {
	*p = 1.0;
	*q = 0.0;
      }
      else {
	*p = 0.0;
	*q = 1.0;
      }
      return;
    }
  //
  //	|z| <= 37.
  //
  expntl = exp(-0.5*zabs*zabs);
  *pdf = expntl/root2pi;
  //
  //	|z| < cutoff = 10/sqrt(2).
  //
  if (zabs < cutoff) 
    *p = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs +
		  p2)*zabs + p1)*zabs + p0)/(((((((q7*zabs + q6)*zabs +
						  q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +
					     q0);
  
  //
  //	|z| >= cutoff.
  //
  else
    *p = *pdf/(zabs + 1.0/(zabs + 2.0/(zabs + 3.0/(zabs + 4.0/
						   (zabs + 0.65)))));
  //
  if (z < 0.0) 
    *q = 1.0 - *p;
  else
    {
      *q = *p;
      *p = 1.0 - *q;
    }
}

double invgauss (float p, int *ifault)
     //
     //	algorithm as241  appl. statist. (1988) vol. 37, no. 3
     //
     //	produces the normal deviate z corresponding to a given lower
     //	tail area of p; z is accurate to about 1 part in 10**16.
     //
     //	the hash sums below are the sums of the mantissas of the
     //	coefficients.   they are included for use in checking
     //	transcription.
     //
{
  static double zero = 0.0, one = 1.0, half = 0.5,
    split1 = 0.4250, split2 = 5.0,
    const1 = 0.180625, const2 = 1.6;
  double ppnd16,q,r;
  //
  //	coefficients for p close to 0.5
  //
  static double a0 = 3.3871328727963666080e0,
    a1 = 1.3314166789178437745e+2,
    a2 = 1.9715909503065514427e+3,
    a3 = 1.3731693765509461125e+4,
    a4 = 4.5921953931549871457e+4,
    a5 = 6.7265770927008700853e+4,
    a6 = 3.3430575583588128105e+4,
    a7 = 2.5090809287301226727e+3,
    b1 = 4.2313330701600911252e+1,
    b2 = 6.8718700749205790830e+2,
    b3 = 5.3941960214247511077e+3,
    b4 = 2.1213794301586595867e+4,
    b5 = 3.9307895800092710610e+4,
    b6 = 2.8729085735721942674e+4,
    b7 = 5.2264952788528545610e+3;
  //	hash sum ab    55.88319 28806 14901 4439
  //
  //	coefficients for p not close to 0, 0.5 or 1.
  //
  static double c0 = 1.42343711074968357734e0,
    c1 = 4.63033784615654529590e0,
    c2 = 5.76949722146069140550e0,
    c3 = 3.64784832476320460504e0,
    c4 = 1.27045825245236838258e0,
    c5 = 2.41780725177450611770e-1,
    c6 = 2.27238449892691845833e-2,
    c7 = 7.74545014278341407640e-4,
    d1 = 2.05319162663775882187e0,
    d2 = 1.67638483018380384940e0,
    d3 = 6.89767334985100004550e-1,
    d4 = 1.48103976427480074590e-1,
    d5 = 1.51986665636164571966e-2,
    d6 = 5.47593808499534494600e-4,
    d7 = 1.05075007164441684324e-9;
  //	hash sum cd    49.33206 50330 16102 89036
  //
  //	coefficients for p near 0 or 1.
  //
  static double e0 = 6.65790464350110377720e0,
    e1 = 5.46378491116411436990e0,
    e2 = 1.78482653991729133580e0,
    e3 = 2.96560571828504891230e-1,
    e4 = 2.65321895265761230930e-2,
    e5 = 1.24266094738807843860e-3,
    e6 = 2.71155556874348757815e-5,
    e7 = 2.01033439929228813265e-7,
    f1 = 5.99832206555887937690e-1,
    f2 = 1.36929880922735805310e-1,
    f3 = 1.48753612908506148525e-2,
    f4 = 7.86869131145613259100e-4,
    f5 = 1.84631831751005468180e-5,
    f6 = 1.42151175831644588870e-7,
    f7 = 2.04426310338993978564e-15;
  //	hash sum ef    47.52583 31754 92896 71629
  //
  *ifault = 0;
  q = p - half;
  if (fabs(q) <= split1) {
    r = const1 - q * q;
    ppnd16 = q * (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
		    * r + a2) * r + a1) * r + a0) /
      (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
	 * r + b2) * r + b1) * r + one);
	  return ppnd16;
  }
  else {
    if (q < zero) 
      r = p;
    else
      r = one - p;
  if (r <= zero) {
    *ifault = 1;
    ppnd16 = zero;
    return ppnd16;
  }
  r = sqrt(-log(r));
  if (r <= split2) {
    r = r - const2;
    ppnd16 = (((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
		* r + c2) * r + c1) * r + c0) /
      (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
	 * r + d2) * r + d1) * r + one);
  }
  else {
    r = r - split2;
    ppnd16 = (((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
		* r + e2) * r + e1) * r + e0) /
      (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
	 * r + f2) * r + f1) * r + one);
  }
  if (q < zero) 
    ppnd16 = - ppnd16;
  return ppnd16;
  }
}

// // Polynomial interpolation from Mumerical Recipes with some modifications
// // I.V. 26-08-99
// 
// 
// void polint(fltarray & xa, fltarray & ya, int n, float x, 
// 	    float *y, float *dy)
// {
//   int i,m,ns=0;
//   float den,dif,dift,ho,hp,w;
//   
//   fltarray c,d;
//   
//   dif=fabs(x-xa(0));
//   c.alloc(n);
//   d.alloc(n);
//   
//   for (i=0;i<n;i++) {
//     if ( (dift=fabs(x-xa(i))) < dif) {
//       ns=i;
//       dif=dift;
//     }
//     c(i)=ya(i);
//     d(i)=ya(i);
//   }
// 
//   *y=ya(ns--);
// 
//   for (m=1;m<n;m++) {
//     for (i=0;i<n-m-1;i++) {
//       ho=xa(i)-x;
//       hp=xa(i+m)-x;
//       w=c(i+1)-d(i);
//       if ( (den=ho-hp) == 0.0) cerr << "Error in routine polint\n";
//       den=w/den;
//       d(i)=hp*den;
//       c(i)=ho*den;
//     }
//     *y += (*dy=(2*ns < (n-m) ? c(ns+1) : d(ns--)));
//   }
//   d.free();
//   c.free();
// }

void cf_indexx(int n, fltarray & arr, intarray & indx)
{
  int i,indxt,ir=n-1,itemp,j,k,l=0;
  int jstack=0;
  intarray istack;
  float a;
  
  istack.alloc(NSTACK);
  
  for (j=0;j<n;j++) 
    indx(j)=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx(j);
	a=arr(indxt);
	for (i=j-1;i>=0;i--) {
	  if (arr(indx(i)) <= a) break;
	  indx(i+1)=indx(i);
	}
	indx(i+1)=indxt;
      }
      if (jstack == 0) break;
      ir=istack(jstack--);
      l=istack(jstack--);
    } else {
      k=(l+ir) >> 1;
      SWAP(indx(k),indx(l+1));
      if (arr(indx(l+1)) > arr(indx(ir))) {
	SWAP(indx(l+1),indx(ir))
	  }
      if (arr(indx(l)) > arr(indx(ir))) {
	SWAP(indx(l),indx(ir))
	  }
      if (arr(indx(l+1)) > arr(indx(l))) {
	SWAP(indx(l+1),indx(l))
	  }
      i=l+1;
      j=ir;
      indxt=indx(l);
      a=arr(indxt);
      for (;;) {
	do i++; while (arr(indx(i)) < a);
	do j--; while (arr(indx(j)) > a);
	if (j < i) break;
	SWAP(indx(i),indx(j))
	  }
      indx(l)=indx(j);
      indx(j)=indxt;
      jstack += 2;
      if (jstack > NSTACK) cerr << "NSTACK too small in indexx.\n";
      if (ir-i+1 >= j-l) {
	istack(jstack)=ir;
	istack(jstack-1)=i;
	ir=j-1;
      } else {
	istack(jstack)=j-1;
	istack(jstack-1)=l;
	l=i;
      }
    }
  }
  istack.free();
}

// Calculate mean with missing values

float xmean(fltarray &xxx, int miss)
{
  int j=0,n;
  double sum=0;
  n = xxx.n_elem();

  for (int i=0;i<n;i++) 
    {
      if (xxx(i) != miss) {
	sum += xxx(i);
	j++;
      }
    }
  if (j == 0) 
    {
      cout << "WARNING! All values are missing? " << endl;
      return xxx(0);
    }
  else
    return ((float) (sum/(double)j));
}


// Calculate st.dev. with missing values

float xstdev(fltarray &xxx, int miss)
{
  int i,j=0,n;
  double sum=0.0,std=0.0;
  n = xxx.n_elem();

  for (i=0;i<n;i++) 
    {
      if (xxx(i) != miss) {
	sum += xxx(i);
	j++;
      }
    }
  if (j == 0) 
    {
      cout << "WARNING! All values are missing? " << endl;
      return -1;
    }
  else sum /= (double) j;

  for (i=0;i<j;i++) std += (xxx(i) - sum)*(xxx(i) - sum);
  std /= (double)(j - 1.0);
  return ((float) (sqrt(std)));
}

void efron(float mmm, fltarray &xxx, int miss, float *s1, float *s2)
{
  int n,i,j,status,iii;
  fltarray yy,tmp,xx;
  intarray ix;
  float gxx,arg,tmp1,tmp2,dummy,gmean,xmin,xmax;
  double pdf;
  //  float dgxx;
  /*  float t0 = 0.16, t1 = 0.84, gt0 = -1.0, gt1 = 1.0;*/

  n = xxx.n_elem();
  yy.alloc(n);
  tmp.alloc(n);
  xx.alloc(n);
  ix.alloc(n);
  
  for (i=0,j=0;i<n;i++)
    {
      if (xxx(i) != miss)
	{
	  tmp(j) = xxx(i);
	  j++;
	}
    }
  xmin = tmp.min();  
  xmax = tmp.max();  
  
  if (mmm <= xmin || mmm >= xmax)
    {
      cerr << "Warning: the input mean is outside the distribution!\n";
      cerr << "         Using the distribution mean instead\n";
      mmm = tmp.mean();
    }
  
  // index sorting 
  
  cf_indexx(j,tmp,ix);
  
  for (i=0;i<j;i++)
    {
      xx(i) = tmp(ix(i));
      // the integral distribution G(xx)
      yy(i) = (float)i/(float)j;
    }

  // now interpolating to find G(mean)
  //  polint(xx,yy,j,mmm,&gxx,&dgxx);
  cf_locate(xx,j,mmm,&iii);
  gxx = (yy(iii) + yy(iii+1))/2.0;

  gmean = fabs(invgauss(gxx,&status));

  //left point
  arg = 2.0*gmean - 1.0;
  gaussint(arg,&tmp1,&tmp2,&pdf);
  dummy = MIN(tmp1,tmp2);

  // Now interchange the arrays for interpolation to get G^{-1}
  //polint(yy,xx,j,dummy,&gxx,&dgxx)
  cf_locate(yy,j,dummy,&iii);
  *s1 = (xx(iii) + xx(iii+1))/2.0;
  
  // right point
  arg = 2.0*gmean + 1.0;
  gaussint(arg,&tmp1,&tmp2,&pdf);
  dummy = MIN(tmp1,tmp2);
  // Now interchange the arrays for interpolation to get G^{-1}
  cf_locate(yy,j,dummy,&iii);
  *s2 = (xx(iii) + xx(iii+1))/2.0;
  //polint(yy,xx,j,dummy,&gxx,&dgxx);
  // *s2 = gxx;
}

void cf_locate(fltarray &xx, int n, float x, int *j)
{
	int ju,jm,jl;
	int ascnd;

	jl=-1;
	ju=n;
	ascnd=(xx(n-1) > xx(0));
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx(jm) == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
}
/* (C) Copr. 1986-92 Numerical Recipes Software '?t5$. */
