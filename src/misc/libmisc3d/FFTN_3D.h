/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/03/01
**    
**    File:  FFTN_3D.h
**
*******************************************************************************
**
**    DESCRIPTION  3D FFT Class routine for any size of cubes
**    -----------   
**                 
******************************************************************************/

#ifndef _CFFTN_3D_H_
#define _CFFTN_3D_H_

#include "GlobalInc.h"
#include "FFTN.h"

class FFTN_3D: public FFTN
{
	private:
		inline void pix_swap(complex_f &a, complex_f &b) { complex_f temp=a;a=b;b=temp;}
		inline void pix_swap(double &a, double &b) { double temp=a;a=b;b=temp;}
		inline void pix_swap(complex_d &a, complex_d &b) { complex_d temp=a;a=b;b=temp;}


		void uncenter(fltarray& Cube);
		void uncenter(cfarray& Cube);
		void center(fltarray& Cube);
		void center(cfarray& Cube);

		void cb_shift(fltarray& Data1, fltarray& Data2, int Dx, int Dy, int Dz); 
		void cb_shift(cfarray& Data1, cfarray& Data2, int Dx, int Dy, int Dz); 


	public:
		Bool CenterZeroFreq; // if True, the zero frequency is in the middle
							// of the image, else at the left

		FFTN_3D(){CenterZeroFreq=False;} // Constructor

		// 3D FFT transform of a real signal
		void fftn3d(fltarray &Signal, complex_f *Buff, Bool Reverse=False, bool normalize=false);
		void fftn3d(fltarray &Signal, complex_d *Buff, Bool Reverse=False, bool normalize=false);
		void fftn3d(dblarray &Signal, complex_d *Buff, Bool Reverse=False, bool normalize=false);
		inline void ifftn3d(fltarray &Cube, complex_f *Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}
		inline void ifftn3d(fltarray &Cube, complex_d *Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}
		inline void ifftn3d(dblarray &Cube, complex_d *Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}

		// 3D FFT transform of a complex signal
		void fftn3d (complex_f *Buff, int Nx, int Ny, int Nz, Bool Reverse=False, bool normalize=false);
		void fftn3d (complex_d *Buff, int Nx, int Ny, int Nz, Bool Reverse=False, bool normalize=false);

		// 3D FFT transform of a real signal
		void fftn3d(fltarray &Signal, cfarray & Buff, Bool Reverse=False, bool normalize=false);
		void fftn3d(fltarray &Signal, cdarray & Buff, Bool Reverse=False, bool normalize=false);
		void fftn3d(dblarray &Signal, cdarray & Buff, Bool Reverse=False, bool normalize=false);
		inline void ifftn3d(fltarray &Cube, cfarray &Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}
		inline void ifftn3d(fltarray &Cube, cdarray &Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}
		inline void ifftn3d(dblarray &Cube, cdarray &Buff, bool normalize=false)
		{ fftn3d(Cube, Buff, True, normalize);}

		// 3D FFT transform of a complex signal
		void fftn3d (cfarray &Buff,  Bool Reverse=False, bool normalize=false);
		void fftn3d (cdarray &Buff,  Bool Reverse=False, bool normalize=false);
		inline void ifftn3d(cfarray &Buff, bool normalize=false)
		{ fftn3d(Buff, True, normalize);}
		inline void ifftn3d(cdarray &Buff, bool normalize=false)
		{ fftn3d(Buff, True, normalize);}

		void swap_buff(complex_d *Buff, int Nx, int Ny, int Nz, Bool Reverse=False);
		void swap_buff(complex_f *Buff, int Nx, int Ny, int Nz, Bool Reverse=False);
		// Swap the complex data (i.e. change the position of the zero 
		// frequency)

		void convolve(fltarray & Data, fltarray & Data2);
		// Convolve two cubes: Data = Data * Data2

		void convolve(fltarray & Data1, fltarray & Data2, fltarray & Data3)
		{
			Data3 = Data1;
			convolve(Data3, Data2);
		}
		
		~FFTN_3D() {} // deallocation
};


void convolve3d(fltarray & Data, fltarray & Gauss);

#endif

