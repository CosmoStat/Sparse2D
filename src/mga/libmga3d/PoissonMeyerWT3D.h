/******************************************************************************
**                   Copyright (C) 2010 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  Sept 2010
**    
**    File:  PoissonMeyerWT.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  3D POISSON MEYER WAVELET TRANSFORM 
**    ----------- 
******************************************************************************/

#ifndef _POISSON_MWT3D
#define _POISSON_MWT3D

#include "GlobalInc.h"
#include "FFTN.h"
#include "FFTN_3D.h"
#include "IM_IO.h"
#include "MGA_Inc.h"

/***********************************************************************/

class POISSON_MWT3D
{
protected:
	// Function used to compute the low pass filtering
	fltarray H;  // H filter
	inline double lowpass_window(double x)
	{
		double l,r;
		l = exp( 1-1 / (1-exp(1-1/(1-x))));
		r = exp(1-1/(1-exp(1-1/x)));
		return (l /= sqrt(l*l+r*r));
	}

	int Nx_Cube;     // Input cube x size
	int Ny_Cube;     // Input cube y size
	int Nz_Cube;     // Input cube z size
	int NbrScale;    // Nbr scales in the WT

	cfarray TF_ExtData;	// Fourier transform of the extended cube
						// If Extend = False then TF_ExtData = FFT(input_data)
	cfarray *Tabcf_WT_Band; // Wavelet transform of the cube in Fourier domain
	
	intarray TabNx;	// TabNx[s] = x size at scale s
	intarray TabNy;	// TabNy[s] = y size at scale s
	intarray TabNz;	// TabNz[s] = z size at scale s

	Bool NeedOddSize;	// If NeedOddSize==True, all scales will have odd sizes
                		// In this case, the input cube must have an odd size

	int ExtNx;      // x size of the cube to be transformed
	int ExtNy;      // = 2*floor(D_ExtNx/2)+1
	int ExtNz;

	double D_ExtNx; // x size of the cube to be transformed
	double D_ExtNy; // D_ExtNx = Nx 
	double D_ExtNz;  

	Bool Isotrop;	// If Isotrop == True, the low pass filter is isotropic. 
            	 	//    and it is not Meyer Wavelets anymore.

	void get_hfilter(fltarray &H, double DNx=0., double DNy=0., double DNz=0.);
	// Calculate the low pass filtering in the Fourier domain to be 
	// applied at each step of the wavelet decomposition.

	void get_extFourier(fltarray &Data, cfarray & TF_ExtData);
	//  Take the FFT of the input cube: Data
	//  TF_ExtData contains the FFT of the [extended] input cube.

public:
	POISSON_MWT3D();
	~POISSON_MWT3D();
	FFTN_3D FFT3D;   // public for test purposes, else protected
	
	type_pmeyer3D Type_PMW3D;

// Initialize the WT class for a given number of scales and a give cube size
	void init(int Nbr_Scale, int Nx, int Ny, int Nz, Bool IsotropWT=False, Bool WTNeedOddSize=False);

// Computes the wavelet transform  of a cube
	void transform(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc=True);
	void transform_sum(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc=True);
	void transform_coherent(fltarray &Data, fltarray * & Tab_Cube, Bool Alloc=True);
	// Data: in = input cube
	// TabWT: out = TabWT[s] is the scale s of the WT
	// If Alloc == True then Tab_Cube is allocated

// Reconstruction of a cube from its WT.
	void recons(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc=True);
	void recons_sum(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc=True);
	void recons_coherent(fltarray * & Tab_Cube, fltarray &Data, Bool Alloc=True);
	// Tab_Cube: in = input wavelet scales
	// Data: out = output reconstructed  cube
	// If Alloc == True then Data is allocated

// Size and position tools
	inline int nxcube() {return Nx_Cube;}
	inline int nycube() {return Ny_Cube;}
	inline int nzcube() {return Nz_Cube;}
	inline double nxd() {return D_ExtNx;}
	inline double nyd() {return D_ExtNy;}
	inline double nzd() {return D_ExtNz;}
	inline int nxs(int s) {return TabNx(s);}
	inline int nys(int s) {return TabNy(s);}
	inline int nzs(int s) {return TabNz(s);}
	inline Bool isotrop() {return Isotrop;}
	int nbr_scale() {return NbrScale;} // Number of scales of the WT

// Statistic and information tools
	float get_norm(int s);		// gives the normalizing coefficient at scale s
	void extract_stat(fltarray *TabBand, char* Outname, bool normalize=true);
	void normalize(fltarray *TabBand, fltarray *TabBandNorm );
	void normalize_self(fltarray *TabBand, bool inverse=false);

// Filtering methods
	void threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void wiener(fltarray *TabBand, float noise_lvl, int LocalBS);
	
// IO tools
	void write(char *Name, fltarray * & Tab_WCube, bool Normalize);
	void read(char *Name, fltarray * & Tab_WCube, bool *NormalizeInv);
	
	void test(fltarray *TabBand);
};

/***********************************************************************/
//	void write_mono(char *Name, fltarray * & Tab_WCube);
//	void read_mono(char *Name, fltarray * & Tab_WCube);

#endif
