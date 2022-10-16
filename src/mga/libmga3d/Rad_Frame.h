
#ifndef _RADTFRAME_H
#define _RADTFRAME_H

#include "Array.h"
#include "IM_IO.h"
#include "IM_Obj.h"
#include "Usage.h"
#include "IM_Radon.h"
#include "writefits3d.h"
#include "MGA_Inc.h"

/***********************************************************************/


class Rad_Frame
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNz;
	int RadNl; // Output radon size
	int RadNc;
	Radon* TabRadon;
	
//Methods
	
public:
	Rad_Frame();
	~Rad_Frame();
	
	void alloc(int nx, int ny, int nz);
	
// Apply the Radon transform and store the result in TabBand
	void transform(fltarray &Data, fltarray &TabBand, bool allocTB);
	
// Reconstruct a cube from its radon coefficients
	void recons(fltarray &TabBand, fltarray &Data);
	
// Statistic and information tools
//	void get_norm_coeff(float NSigma);

// Filtering methods
	void threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType);

};

 
/***********************************************************************/

#endif
