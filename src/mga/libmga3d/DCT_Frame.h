
#ifndef _DCTTFRAME_H
#define _DCTTFRAME_H

#include "Array.h"
#include "IM_IO.h"
#include "IM_Block2D.h"
#include "IM_Obj.h"
#include "writefits3d.h"
#include "MGA_Inc.h"
#include "fftw3.h"

/***********************************************************************/


class DCT_Frame
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNz;
	
	int TransNx;
	int TransNy;
	
	bool BlockOverlap;
	int BlockSize;
	int skip_order;
	Block2D B2D, B2DTrans;
	
//Methods
	
public:
	DCT_Frame();
	~DCT_Frame();
	
	void alloc(int nx, int ny, int nz, int BlockSize, bool Overlap);
	void set_skip_order(int _skip_order) {skip_order = _skip_order;}
	
// Apply the Radon transform and store the result in TabBand
	void transform(fltarray &Data, fltarray &TabBand, bool allocTB=true);
	
// Reconstruct a cube from its radon coefficients
	void recons(fltarray &TabBand, fltarray &Data);
	
// Statistic and information tools
//	void get_norm_coeff(float NSigma);

// Filtering methods
	void threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType);

};

 
/***********************************************************************/

#endif
