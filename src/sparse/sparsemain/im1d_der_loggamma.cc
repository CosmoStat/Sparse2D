/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/02/00
**    
**    File:  dct.cc
**
*******************************************************************************
**
**    DESCRIPTION  dct program
**    ----------- 
**                 
******************************************************************************/

#include "GlobalInc.h"
#include "IM_IO.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int Order=0;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
        
    fprintf(OUTMAN, "         [-o Order]\n");
    fprintf(OUTMAN, "             Order value.\n");
    fprintf(OUTMAN, "             Default is 0.\n");
    manline();
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv,"o:vzZ")) != -1) 
    {
	switch (c) 
        {
           case 'o': 
                if (sscanf(OptArg,"%d",&Order) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad order value: %s\n", OptArg);
                    exit(-1);
                }
  	        break;
          case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
       

       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}


/*******************************************************************/
/*******************************************************************/


static double s_PolyGamma(double x, int order);
 
/*********************************************************************/
/** @file polygamma.c
 * Definitions for portable math library
 */

#include<stdio.h>
#include<ctype.h>
#include<cmath>
#include<cstdio>
#include<cassert>
#include<cstdlib>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<time.h>

#define LOGDERIV_ORDER_MAX	4
#define POLYGAMMA_ORDER_MAX 	LOGDERIV_ORDER_MAX
#define DBL_EPSILON   		2.2204460492503131e-16
#define NCBIMATH_PI   		3.1415926535897932384626433832795
#define NCBIMATH_LNPI   	1.1447298858494001741434273513531
#define NCBIMATH_LN2   		0.69314718055994530941723212145818
#define DIM(A)	(sizeof(A)/sizeof((A)[0]))
#define ABS(a)	((a)>=0?(a):-(a))
// #define MAX(a,b) ((a)>=(b)?(a):(b))
// #define MIN(a,b) ((a)<=(b)?(a):(b))

/** Compute ln(abs(gamma(x))) to 10-digit accuracy
 *  @param x Point to evaluate ln(abs(gamma(x)))
 *  @return The function value
 */
static double s_LnGamma(double x)
{
   return lgamma(x);
}
/** Tabulated values of the first few factorials */
static const double kPrecomputedFactorial[] = {
       1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.,
       39916800., 479001600., 6227020800., 87178291200., 1307674368000.,
       20922789888000., 355687428096000., 6402373705728000., 
       121645100408832000., 2432902008176640000., 51090942171709440000., 
       1124000727777607680000., 25852016738884976640000., 
       620448401733239439360000., 15511210043330985984000000.,
       403291461126605635584000000., 10888869450418352160768000000., 
       304888344611713860501504000000., 8841761993739701954543616000000., 
       265252859812191058636308480000000., 8222838654177922817725562880000000.,
       263130836933693530167218012160000000., 
       8683317618811886495518194401280000000., 
       295232799039604140847618609643520000000.
       };
       
double BLAST_Factorial(int n)
{
   if (n < 0)
      return 0.0; /* Undefined! */

   if (n < DIM(kPrecomputedFactorial))
      return kPrecomputedFactorial[n];

   return exp(s_LnGamma((double)(n + 1)));
}

double BLAST_LnGammaInt(int n)
{
   if ( (n > 1) && (n < DIM(kPrecomputedFactorial) ) ) {
      return log(kPrecomputedFactorial[n-1]);
   }

   return s_LnGamma((double)n);
}

/*
   Romberg numerical integrator
   Reference:
      Francis Scheid (1968)
      Schaum's Outline Series
      Numerical Analysis, p. 115
      McGraw-Hill Book Company, New York
*/

/** Make a parametrized function appear to have only one variable */
#define F(x)  ((*f)((x), fargs))
/** Maximum number of diagonals in the Romberg array */
#define MAX_DIAGS 20

double BLAST_RombergIntegrate(double (*f) (double,void*), void* fargs, double p, double q, double eps, int epsit, int itmin)

{
   double   romb[MAX_DIAGS];   /* present list of Romberg values */
   double   h;   /* mesh-size */
   int      i, j, k, npts;
   long   n;   /* 4^(error order in romb[i]) */
   int      epsit_cnt = 0, epsck;
   double   y;
   double   x;
   double   sum;

   /* itmin = min. no. of iterations to perform */
   itmin = MAX(1, itmin);
   itmin = MIN(itmin, MAX_DIAGS-1);

   /* epsit = min. no. of consecutive iterations that must satisfy epsilon */
   epsit = MAX(epsit, 1); /* default = 1 */
   epsit = MIN(epsit, 3); /* if > 3, the problem needs more prior analysis */

   epsck = itmin - epsit;

   npts = 1;
   h = q - p;
   x = F(p);
   if (ABS(x) == HUGE_VALF)
      return x;
   y = F(q);
   if (ABS(y) == HUGE_VALF)
      return y;
   romb[0] = 0.5 * h * (x + y);   /* trapezoidal rule */
   for (i = 1; i < MAX_DIAGS; ++i, npts *= 2, h *= 0.5) {
      sum = 0.;   /* sum of ordinates for 
		     x = p+0.5*h, p+1.5*h, ..., q-0.5*h */
      for (k = 0, x = p+0.5*h; k < npts; k++, x += h) {
	 y = F(x);
	 if (ABS(y) == HUGE_VALF)
	    return y;
	 sum += y;
      }
      romb[i] = 0.5 * (romb[i-1] + h*sum); /* new trapezoidal estimate */

      /* Update Romberg array with new column */
      for (n = 4, j = i-1; j >= 0; n *= 4, --j)
	 romb[j] = (n*romb[j+1] - romb[j]) / (n-1);

      if (i > epsck) {
	 if (ABS(romb[1] - romb[0]) > eps * ABS(romb[0])) {
	    epsit_cnt = 0;
	    continue;
	 }
	 ++epsit_cnt;
	 if (i >= itmin && epsit_cnt >= epsit)
	    return romb[0];
      }
   }
   return HUGE_VALF;
}

int BLAST_Gcd(int a, int b)
{
   int   c;

   b = ABS(b);
   if (b > a)
      c=a, a=b, b=c;

   while (b != 0) {
      c = a%b;
      a = b;
      b = c;
   }
   return a;
}

int 
BLAST_Gdb3(int* a, int* b, int* c)
{
    int g;
    if (*b == 0) 
	g = BLAST_Gcd(*a, *c);
    else 
	g = BLAST_Gcd(*a, BLAST_Gcd(*b, *c));
    if (g > 1) {
	*a /= g;
	*b /= g;
	*c /= g;
    }
    return g;
}

long BLAST_Nint(double x)
{
   x += (x >= 0. ? 0.5 : -0.5);
   return (long)x;
}


double BLAST_Powi(double x, int n)
{
   double   y;

   if (n == 0)
      return 1.;

   if (x == 0.) {
      if (n < 0) {
	 return HUGE_VALF;
      }
      return 0.;
   }

   if (n < 0) {
      x = 1./x;
      n = -n;
   }

   y = 1.;
   while (n > 0) {
      if (n & 1)
	 y *= x;
      n /= 2;
      x *= x;
   }
   return y;

}

double BLAST_LnFactorial (double x) {

    if( x <= 0.0)
	return 0.0;
    else
	return s_LnGamma(x + 1.0);
	
}

double BLAST_Expm1(double   x)
{
  double    absx = ABS(x);

  if (absx > .33)
    return exp(x) - 1.;

  if (absx < 1.e-16)
    return x;

  return x * (1. + x *
	     (1./2. + x * 
	     (1./6. + x *
	     (1./24. + x * 
	     (1./120. + x *
	     (1./720. + x * 
	     (1./5040. + x *
	     (1./40320. + x * 
	     (1./362880. + x *
	     (1./3628800. + x * 
	     (1./39916800. + x *
	     (1./479001600. + 
	      x/6227020800.))))))))))));
}

/** size of the next series term that indicates convergence
    in the log and polygamma functions */
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

double BLAST_Log1p(double x)
{
   int   i;
   double   sum, y;

   if (ABS(x) >= 0.2)
      return log(x+1.);

   /* Limit the loop to 500 terms. */
   for (i=0, sum=0., y=x; i<500 ; ) {
      sum += y/++i;
      if (ABS(y) < DBL_EPSILON)
	 break;
      y *= x;
      sum -= y/++i;
      if (y < DBL_EPSILON)
	 break;
      y *= x;
   }
   return sum;
}

/** evaluate a specified-order derivative of ln(f(x))
 * @param order The order of derivative to evaluate (0...LOGDERIV_ORDER_MAX).
 *		Derivative order 0 just computes ln(f(x))
 * @param u A list of numerical values of f(x) and its derivatives, all at
 *	    the same point x, to be used within the computations 
 * @return 'order'-th derivative of ln(f(x)) or HUGE_VALF if 
 *	    order is out of range or u[0] is zero
 */
static double s_LogDerivative(int order, double* u)
{
  int    i;
  double  y[LOGDERIV_ORDER_MAX+1];
  double  value, tmp;

  if (order < 0 || order > LOGDERIV_ORDER_MAX) {
    return HUGE_VALF;
  }

  if (order > 0 && u[0] == 0.) {
    return HUGE_VALF;
  }
  for (i = 1; i <= order; i++)
    y[i] = u[i] / u[0];

  switch (order) {
  case 0:
    if (u[0] > 0.)
      value = log(u[0]);
    else {
      return HUGE_VALF;
    }
    break;
  case 1:
    value = y[1];
    break;
  case 2:
    value = y[2] - y[1] * y[1];
    break;
  case 3:
    value = y[3] - 3. * y[2] * y[1] + 2. * y[1] * y[1] * y[1];
    break;
  case 4:
    value = y[4] - 4. * y[3] * y[1] - 3. * y[2] * y[2]
      + 12. * y[2] * (tmp = y[1] * y[1]);
    value -= 6. * tmp * tmp;
    break;
  default:
    return HUGE_VALF;
  }
  return value;
}

/** auxiliary values for computation of derivative of ln(gamma(x)) */
static double _default_gamma_coef [] = {
	 4.694580336184385e+04,
	-1.560605207784446e+05,
	 2.065049568014106e+05,
	-1.388934775095388e+05,
	 5.031796415085709e+04,
	-9.601592329182778e+03,
	 8.785855930895250e+02,
	-3.155153906098611e+01,
	 2.908143421162229e-01,
	-2.319827630494973e-04,
	 1.251639670050933e-10
	 };

/** Compute a specified-order derivative of ln(gamma(x))
 *  evaluated at some point x 
 *  @param x Value at which derivative will be evaluated
 *  @param order Order of derivative (0...POLYGAMMA_ORDER_MAX)
 *  @return 'order'-th derivative of ln(gamma(x)) at specified x.
 *	    Accuracy is to 10 digits for x >= 1
 */
static double
s_GeneralLnGamma(double x, int order)
{
   int       i;
   double     xx, tx;
   double     y[POLYGAMMA_ORDER_MAX+1];
   double    tmp, value;
   double    *coef;
   const int	  xgamma_dim = DIM(_default_gamma_coef);


   xx = x - 1.;  /* normalize from gamma(x + 1) to xx! */
   tx = xx + xgamma_dim;
   for (i = 0; i <= order; ++i) {
      tmp = tx;
      /* sum the least significant terms first */
      coef = &_default_gamma_coef[xgamma_dim];
      if (i == 0) {
	 value = *--coef / tmp;
	 while (coef > _default_gamma_coef)
	    value += *--coef / --tmp;
      }
      else {
	 value = *--coef / BLAST_Powi(tmp, i + 1);
	 while (coef > _default_gamma_coef)
	    value += *--coef / BLAST_Powi(--tmp, i + 1);
	 tmp = BLAST_Factorial(i);
	 value *= (i%2 == 0 ? tmp : -tmp);
      }
      y[i] = value;
   }
   ++y[0];
   value = s_LogDerivative(order, y);
   tmp = tx + 0.5;
   switch (order) {
   case 0:
      value += ((NCBIMATH_LNPI+NCBIMATH_LN2) / 2.)
	       + (xx + 0.5) * log(tmp) - tmp;
      break;
   case 1:
      value += log(tmp) - xgamma_dim / tmp;
      break;
   case 2:
      value += (tmp + xgamma_dim) / (tmp * tmp);
      break;
   case 3:
      value -= (1. + 2.*xgamma_dim / tmp) / (tmp * tmp);
      break;
   case 4:
      value += 2. * (1. + 3.*xgamma_dim / tmp) / (tmp * tmp * tmp);
      break;
   default:
      tmp = BLAST_Factorial(order - 2) * BLAST_Powi(tmp, 1 - order)
	    * (1. + (order - 1) * xgamma_dim / tmp);
      if (order % 2 == 0)
	 value += tmp;
      else
	 value -= tmp;
      break;
   }
   return value;
}

/** Compute, to 10-digit accuracy, a specified order
 *  derivative of ln(abs(gamma(x))).
 *  @param x value at which derivative will be evaluated
 *  @param order Order of derivative (0...POLYGAMMA_ORDER_MAX)
 *		order = 0, 1, 2, corresponds to ln(gamma), 
 *		digamma, trigamma, etc. Note that the value here
 *		is one less than that suggested by the "di" and "tri" 
 *		prefixes of digamma, trigamma, etc. In other words, 
 *		it is truly the order of the derivative.
 *   @return Computed derivative value, or HUGE_VALF if order is out of range
 */
static double s_PolyGamma(double x, int order)
{
   int      k;
   double   value, tmp;
   double   y[POLYGAMMA_ORDER_MAX+1], sx;

   order+=1; // To match the notation of Abramovitz and Stegun.
   
   if (order < 0 || order > POLYGAMMA_ORDER_MAX) {
      return HUGE_VALF;
   }

   if (x >= 1.)
      return s_GeneralLnGamma(x, order);

   if (x < 0.) {
      value = s_GeneralLnGamma(1. - x, order);
      value = ((order - 1) % 2 == 0 ? value : -value);
      if (order == 0) {
	 sx = sin(NCBIMATH_PI * x);
	 sx = ABS(sx);
	 if ( (x < -0.1 && (ceil(x) == x || sx < 2.*DBL_EPSILON))
	       || sx == 0.) {
	    return HUGE_VALF;
	 }
	 value += NCBIMATH_LNPI - log(sx);
      }
      else {
	 y[0] = sin(x *= NCBIMATH_PI);
	 tmp = 1.;
	 for (k = 1; k <= order; k++) {
	    tmp *= NCBIMATH_PI;
	    y[k] = tmp * sin(x += (NCBIMATH_PI/2.));
	 }
	 value -= s_LogDerivative(order, y);
      }
   }
   else {
      value = s_GeneralLnGamma(1. + x, order);
      if (order == 0) {
	 if (x == 0.) {
	    return HUGE_VALF;
	 }
	 value -= log(x);
      }
      else {
	 tmp = BLAST_Factorial(order - 1) * BLAST_Powi(x,  -order);
	 value += (order % 2 == 0 ? tmp : - tmp);
      }
   }
   return value;
}



 
int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;
	cout << "Order = " << Order << endl;   
    }
 
    fltarray Data;    
    fits_read_fltarr(Name_Imag_In, Data, &Header);
    for (int i=0; i < Data.nx(); i++) Data(i) = s_PolyGamma(Data(i), Order);
    fits_write_fltarr(Name_Imag_Out, Data);
    exit(0);
}
