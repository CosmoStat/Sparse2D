#include <vector>
#include <TempArray.h>
#include <MatrixOper.h>
#include "IM_IO.h"
#include "C_OMP.h"
#include "C_EXTRACT_PATCHES.h"
using namespace std;
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <fcntl.h>


C_EXTRACT_PATCHES::C_EXTRACT_PATCHES(dblarray &img)
{
	Image = img;
}
C_EXTRACT_PATCHES::~C_EXTRACT_PATCHES()
{
}

dblarray C_EXTRACT_PATCHES::extract(int &W, int &N)
{
	dblarray Data_out;
	int nx = Image.nx();
	int ny = Image.ny();
	intarray x_list(N),y_list(N);
	gsl_rng *rng;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	int delta = W-1;
	// Generating a random seed using  /dev/urandom
	int randomData = open("/dev/urandom", O_RDONLY);
	unsigned long int s;
	read(randomData, &s, sizeof(s));
	gsl_rng_set(rng, s);

	// Generating x and y position for patches top left pixel
	for (int i=0;i<N;i++)
	{


		x_list(i) = gsl_rng_uniform_int(rng, nx-W);
		y_list(i) = gsl_rng_uniform_int(rng, ny-W);
	}
	// Reading patches, where (i,j) is used as the top left pixel of each patch
	int ind;
	Data_out.alloc(W*W,N);
	for (int p=0;p<N;p++)
	{
		ind = 0;
		for (int j=0;j<W;j++)
			for (int i=0;i<W;i++)
			{
				Data_out(ind,p) = Image(x_list(p)+i,y_list(p)+j);
				ind++;
			}
	}
	return(Data_out);
}






