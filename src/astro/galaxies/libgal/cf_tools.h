
#ifndef	_CFUTIL_H_
#define	_CFUTIL_H_


#include <math.h>
// #include "IM_Math.h"
#include "Array.h"

// #define MIN(a,b) ((a) > (b)) ? (a) : (b) 
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void cf_indexx(int, fltarray &, intarray &);
void gaussint(float, float *, float *, double *);
double invgauss (float, int *);
// void polint(fltarray &, fltarray &, int, float, float *, float *);
float xmean(fltarray &,int);
float xstdev(fltarray &,int);
void efron(float, fltarray &, int, float *, float *);
void cf_locate(fltarray &, int, float, int *);

#endif
