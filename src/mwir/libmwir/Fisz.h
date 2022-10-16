/*
 * Filename : Fisz.h
 * 
 * Class Description
 * 1. FiszTransform : Fisz transform 1D, 2D and 3D; inverse transform 1D, 2D and 3D
 */
 
#ifndef _FISZ_H
#define _FISZ_H
 
#include "Wavelet.h"

using namespace std;
extern bool VERBOSE;

template <class DATATYPE>
class FiszTransform
{		
	public:
			
		// extend the data to the size of each 
		// dimension of the nearest power of two,
		// if the data size is not of power two;
		// ext is an array of 6 elements (extension information of 6 directions)
		void dataExtension (to_array<DATATYPE, true> &data, int ext[], type_border BORDERTYPE=I_MIRROR);
		// extract the original data in using the extension information
		void dataExtraction (to_array<DATATYPE, true> &data, int ext[]);
		
		// Fisz transform 1D / 2D / 3D
		// the data size MUST be power of two, 
		// otherwise, a DataSizeException will be thrown out
		void fisz1D (to_array<DATATYPE, true> &data);
		void fisz2D (to_array<DATATYPE, true> &data);
		void fisz3D (to_array<DATATYPE, true> &data);

		// Inverse Fisz transform 1D / 2D / 3D
		// data size MUST be power of two.
		// otherwise, a DataSizeException will be thrown out
		void ifisz1D (to_array<DATATYPE, true> &data);
		void ifisz2D (to_array<DATATYPE, true> &data);
		void ifisz3D (to_array<DATATYPE, true> &data);
};

template <class DATATYPE>
void FiszTransform<DATATYPE>::dataExtension (to_array<DATATYPE, true> &data, int ext[], type_border BORDERTYPE)
{
	int lext = 0, rext = 0, uext = 0, dext = 0, bext = 0, fext = 0;
	int dim = data.naxis();
	int len1 = data.nx(), len2 = data.ny(), len3 = data.nz();

	if ((dim == 1) && !is_power_of_2(len1))
	{
		int newlen = next_power_of_2(len1);
		lext = (newlen - len1) / 2; 
		rext = (newlen - len1) - lext;
		
		to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen);
		for (int x=0; x<newlen; x++)
			(*temp)(x) = data(x-lext, BORDERTYPE);
		data = *temp;
		if (temp != NULL) { delete temp; temp = NULL; }
	}
	else if ((dim == 2) && (!is_power_of_2(len1) || !is_power_of_2(len2) || len1!=len2))
	{
		int newlen1 = len1, newlen2 = len2;
		
        int maxlen = MAX(newlen1, newlen2);
		if (!is_power_of_2(maxlen))
			maxlen = next_power_of_2(maxlen);        
        newlen1 = maxlen; newlen2 = maxlen;
		
		lext = (newlen1 - len1) / 2; rext = (newlen1 - len1) - lext;
		uext = (newlen2 - len2) / 2; dext = (newlen2 - len2) - uext;
		
		to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen1, newlen2);
		for (int y=0; y<newlen2; y++)
			for (int x=0; x<newlen1; x++)
				(*temp)(x, y) = data(x-lext, y-uext, BORDERTYPE);
		data = *temp;
		if (temp != NULL) { delete temp; temp = NULL; }
	}
	else if ((dim == 3) && (!is_power_of_2(len1) || !is_power_of_2(len2) || !is_power_of_2(len3) || len1!=len2 || len2!=len3))
	{
		int newlen1 = len1, newlen2 = len2, newlen3 = len3;
		
        int maxlen = MAX(MAX(newlen1, newlen2), newlen3);
		if (!is_power_of_2(maxlen))
			maxlen = next_power_of_2(maxlen);        
        newlen1 = maxlen; newlen2 = maxlen; newlen3 = maxlen;

		lext = (newlen1 - len1) / 2; rext = (newlen1 - len1) - lext;
		uext = (newlen2 - len2) / 2; dext = (newlen2 - len2) - uext;
		bext = (newlen3 - len3) / 2; fext = (newlen3 - len3) - fext;
		
		to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(newlen1, newlen2, newlen3);
		for (int z=0; z<newlen3; z++)
			for (int y=0; y<newlen2; y++)
				for (int x=0; x<newlen1; x++)
					(*temp)(x, y, z) = data(x-lext, y-uext, z-bext, BORDERTYPE);
		data = *temp;
		if (temp != NULL) { delete temp; temp = NULL; }
	}
	
	ext[0] = lext; ext[1] = rext; ext[2] = uext; 
	ext[3] = dext; ext[4] = bext; ext[5] = fext;
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::dataExtraction (to_array<DATATYPE, true> &data, int ext[])
{
	int dim = data.naxis();
	int lext = ext[0], rext = ext[1], uext = ext[2];
	int dext = ext[3], bext = ext[4], fext = ext[5];
	
	if ((dim == 1) && ((lext != 0) || (rext != 0)))
	{
		int finallen = data.nx() - lext - rext;
		to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen);
		for (int x=0; x<finallen; x++) (*tdata)(x) = data(lext+x);
		data = (*tdata);
		if (tdata != NULL) { delete tdata; tdata = NULL; }
	}
	else if ((dim == 2) && \
		((lext != 0) || (rext != 0) || (uext != 0) || (dext != 0)))
	{
		int finallen1 = data.nx() - lext - rext;
		int finallen2 = data.ny() - uext - dext;
		to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen1, finallen2);
		for (int x=0; x<finallen1; x++) 
			for (int y=0; y<finallen2; y++)
				(*tdata)(x, y) = data(lext+x, uext+y);
		data = *tdata;
		if (tdata != NULL) { delete tdata; tdata = NULL; }
	}
	else if ((dim == 3) && ((lext != 0) || (rext != 0) || \
		(uext != 0) || (dext != 0) || \
		(bext != 0) || (fext != 0)))
	{
		int finallen1 = data.nx() - lext - rext;
		int finallen2 = data.ny() - uext - dext;
		int finallen3 = data.nz() - bext - fext;
		to_array<DATATYPE,true> *tdata = new to_array<DATATYPE, true>(finallen1, finallen2, finallen3);
		for (int x=0; x<finallen1; x++) 
		for (int y=0; y<finallen2; y++)
		for (int z=0; z<finallen3; z++)	
			(*tdata)(x, y, z) = data(lext+x, uext+y, bext+z);
		data = *tdata;
		if (tdata != NULL) { delete tdata; tdata = NULL; }
	}	
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::fisz1D (to_array<DATATYPE, true> &data)
{
    double ca, cd;
	int len = data.nx();	
	if (len <= 1) return;
	else if (!is_power_of_2(len)) 
		throw DataSizeException(len, "data size not of power of two");
	to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(len);
	
	// decomposition by Haar transform and modification of detail coefficients
	while (len > 1)
	{
		for (int x=0; x<len; x+=2)
		{
			ca = (data(x)+data(x+1)) / 2;
			(*temp)(x/2) = ca;
			if (ca == 0)
				(*temp)(x/2 + len/2) = 0;
			else
				(*temp)(x/2 + len/2) = (data(x) - data(x+1)) / (2.0 * psqrt(ca));
		}
		for (int x=0; x<len; x++) data(x) = (*temp)(x);
		len /= 2;
	} 	
	
	// reconstruction by inverse Haar transform
	int s = 1; len = data.nx();
	while (s < len)
	{
		for (int x=0; x<s; x++)
		{
			ca = data(x); cd = data(x+s);
			(*temp)(2*x) = ca + cd;
			(*temp)(2*x+1) = ca - cd;
		}
		s *= 2;
		for (int x=0; x<s; x++) data(x) = (*temp)(x);
	}
	
	if (temp != NULL) { delete temp; temp = NULL; }
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::ifisz1D (to_array<DATATYPE, true> &data)
{
    double ca, cd;
	int len = data.nx();
	if (len <= 1) return;
	else if (!is_power_of_2(len)) 
		throw DataSizeException(len, "data size not of power of two");
	to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(len);
	
	// decomposition by Haar transform
	while (len > 1)
	{
		for (int x=0; x<len; x+=2)
		{
			(*temp)(x/2) = (data(x)+data(x+1)) / 2;
			(*temp)(x/2 + len/2) = (data(x) - data(x+1)) / 2;
		}
		for (int x=0; x<len; x++) data(x) = (*temp)(x);
		len /= 2;
	} 	

	// modification of detail coefficients and reconstruction by inverse Haar transform
	int s = 1; len = data.nx();
	while (s < len)
	{
		for (int x=0; x<s; x++)
		{
			ca = data(x); cd = data(x+s);			
			(*temp)(2*x) = ca + cd * psqrt(ca);
			(*temp)(2*x+1) = ca - cd * psqrt(ca);
		}
		s *= 2;
		for (int x=0; x<s; x++) data(x) = (*temp)(x);
	}
	
	if (temp != NULL) { delete temp; temp = NULL; }
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::fisz2D (to_array<DATATYPE, true> &data)
{
	double ca, dh, dv, dd;
	int s, offsetx, offsety;
	int len1 = data.nx(), len2 = data.ny();
	if ((len1 <= 1) || (len2 <= 1)) return;
	else if (!is_power_of_2(len1))
		throw DataSizeException(len1, "data size not of power of two");
	else if (!is_power_of_2(len2))
		throw DataSizeException(len2, "data size not of power of two");	
	to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(len1, len2);
	
	// decomposition by Haar transform and modification of detail coefficients
	while ((len1 > 1) && (len2 > 1))
	{
		for (int x=0; x<len1; x+=2)
			for (int y=0; y<len2; y+=2)
			{
				ca = (data(x,y)+data(x+1,y)+data(x,y+1)+data(x+1,y+1)) / 4;
				dh = (data(x,y)+data(x+1,y)-data(x,y+1)-data(x+1,y+1)) / 4;
				dv = (data(x,y)+data(x,y+1)-data(x+1,y)-data(x+1,y+1)) / 4;
				dd = (data(x,y)+data(x+1,y+1)-data(x,y+1)-data(x+1,y)) / 4;
				
				(*temp)(x/2, y/2) = ca;
				if (ca == 0)
				{
					(*temp)(x/2 + len1/2, y/2) = 0;
					(*temp)(x/2, y/2+len2/2) = 0;
					(*temp)(x/2 + len1/2, y/2 + len2/2) = 0;
				}
				else
				{
					(*temp)(x/2 + len1/2, y/2) = dh / psqrt(ca);
					(*temp)(x/2, y/2 + len2/2) = dv / psqrt(ca);
					(*temp)(x/2 + len1/2, y/2 + len2/2) = dd / psqrt(ca);
				}
			}

		for (int x=0; x<len1; x++) 
			for (int y=0; y<len2; y++)
				data(x, y) = (*temp)(x, y);
		len1 /= 2;
		len2 /= 2;				
	}

	// reconstruction by inverse Haar transform
	s = 1; 
	len1 = data.nx(); len2 = data.ny();
	int rx = len1 / len2, ry = len2 / len1;
	int r = MIN(rx, ry);
	while ((s < len1) && (s < len2))
	{
		for (int b=0; b<r; b++)
		{
			offsetx = (rx > 1) ? s * b : 0;
			offsety = (ry > 1) ? s * b : 0;
			for (int x=0; x<s; x++)
				for (int y=0; y<s; y++)
				{
					ca = data(x+offsetx, y+offsety); 
					dh = data(x+s+offsetx, y+offsety);
					dv = data(x+offsetx, y+s+offsety);
					dd = data(x+s+offsetx, y+s+offsety);
				
					(*temp)(2*(x+offsetx), 2*(y+offsety)) = ca + dh + dv + dd;
					(*temp)(2*(x+offsetx)+1, 2*(y+offsety)) = ca + dh - dv - dd;
					(*temp)(2*(x+offsetx), 2*(y+offsety)+1) = ca - dh + dv - dd;
					(*temp)(2*(x+offsetx)+1, 2*(y+offsety)+1) = ca - dh - dv + dd;
				}
		}
		
		s *= 2;
		for (int b=0; b<r; b++)
		{
			offsetx = (rx > 1) ? s * b : 0;
			offsety = (ry > 1) ? s * b : 0;
			for (int x=0; x<s; x++) 
				for (int y=0; y<s; y++)
					data(x+offsetx, y+offsety) = (*temp)(x+offsetx, y+offsety);
		}
	}
	
	if (temp != NULL) { delete temp; temp = NULL; }
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::ifisz2D (to_array<DATATYPE, true> &data)
{
	double ca, dh, dv, dd;
	int s, offsetx, offsety;
	int len1 = data.nx(), len2 = data.ny();
	if ((len1 <= 1) || (len2 <= 1)) return;
	else if (!is_power_of_2(len1))
		throw DataSizeException(len1, "data size not of power of two");
	else if (!is_power_of_2(len2))
		throw DataSizeException(len2, "data size not of power of two");
	to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(len1, len2);
	
	// decomposition by Haar transform
	while ((len1 > 1) && (len2 > 1))
	{
		for (int x=0; x<len1; x+=2)
			for (int y=0; y<len2; y+=2)
			{
				ca = (data(x,y)+data(x+1,y)+data(x,y+1)+data(x+1,y+1)) / 4;
				dh = (data(x,y)+data(x+1,y)-data(x,y+1)-data(x+1,y+1)) / 4;
				dv = (data(x,y)+data(x,y+1)-data(x+1,y)-data(x+1,y+1)) / 4;
				dd = (data(x,y)+data(x+1,y+1)-data(x,y+1)-data(x+1,y)) / 4;
				
				(*temp)(x/2, y/2) = ca;
				(*temp)(x/2 + len1/2, y/2) = dh;
				(*temp)(x/2, y/2 + len2/2) = dv;
				(*temp)(x/2 + len1/2, y/2 + len2/2) = dd;
			}
		for (int x=0; x<len1; x++) 
			for (int y=0; y<len2; y++)
				data(x, y) = (*temp)(x, y);
		len1 /= 2;
		len2 /= 2;				
	}

	// modification of detail coefficients and reconstruction by inverse Haar transform
	s = 1; 
	len1 = data.nx(); len2 = data.ny();
	int rx = len1 / len2, ry = len2 / len1;
	int r = MIN(rx, ry);

	while ((s < len1) && (s < len2))
	{
		for (int b=0; b<r; b++)
		{
			offsetx = (rx > 1) ? s * b : 0;
			offsety = (ry > 1) ? s * b : 0;

			for (int x=0; x<s; x++)
				for (int y=0; y<s; y++)
				{
					ca = data(x+offsetx, y+offsety); 
					dh = data(x+s+offsetx, y+offsety);
					dv = data(x+offsetx, y+s+offsety);
					dd = data(x+s+offsetx, y+s+offsety);
				
					(*temp)(2*(x+offsetx), 2*(y+offsety)) = ca + (dh + dv + dd) * psqrt(ca);
					(*temp)(2*(x+offsetx)+1, 2*(y+offsety)) = ca + (dh - dv - dd) * psqrt(ca);
					(*temp)(2*(x+offsetx), 2*(y+offsety)+1) = ca + (- dh + dv - dd) * psqrt(ca);
					(*temp)(2*(x+offsetx)+1, 2*(y+offsety)+1) = ca + (- dh - dv + dd) * psqrt(ca);
				}
		}
		
		s *= 2;
		for (int b=0; b<r; b++)
		{
			offsetx = (rx > 1) ? s * b : 0;
			offsety = (ry > 1) ? s * b : 0;
			for (int x=0; x<s; x++) 
				for (int y=0; y<s; y++)
					data(x+offsetx, y+offsety) = (*temp)(x+offsetx, y+offsety);
		}
	}
	
	if (temp != NULL) { delete temp; temp = NULL; }	
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::fisz3D (to_array<DATATYPE, true> &data)
{
	int len1 = data.nx(), len2 = data.ny(), len3 = data.nz();
	int len = MIN(MIN(len1, len2), len3);
	if (len <= 1) return;
	else if (!is_power_of_2(len1))
		throw DataSizeException(len1, "data size not of power of two");
	else if (!is_power_of_2(len2))
		throw DataSizeException(len2, "data size not of power of two");
	else if (!is_power_of_2(len3))
		throw DataSizeException(len3, "data size not of power of two");
	
	int s = iilog2(len);
	to_array<DATATYPE, true> *coef = new to_array<DATATYPE, true>[7*s + 1];
	
	dblarray dfilterh(2), dfilterg(2), rfilterh(2), rfilterg(2);
	dfilterh(0) = .5; dfilterh(1) = .5; dfilterg(0) = -.5; dfilterg(1) = .5;
    rfilterh(0) = 1; rfilterh(1) = 1; rfilterg(0) = 1; rfilterg(1) = -1;
    
    // decomposition by Haar transform and modification of detail coefficients
	for (int scale=0; scale<s; scale++)
	{
		OrthogonalWaveletTransform<DATATYPE>::dwt3D(data, dfilterh, dfilterg, \
			coef[0], coef[scale+1], coef[s+scale+1], \
			coef[2*s+scale+1], coef[3*s+scale+1], coef[4*s+scale+1], \
			coef[5*s+scale+1], coef[6*s+scale+1]);
		data = coef[0];
		
		int lx = coef[0].nx(), ly = coef[0].ny(), lz = coef[0].nz();
		for (int x=0; x<lx; x++)
		for (int y=0; y<ly; y++)
		for (int z=0; z<lz; z++)
			coef[0](x, y, z) = psqrt(coef[0](x, y, z));
			
		for (int i=0; i<7; i++)
			coef[i*s+scale+1] /= coef[0];
	}
	
	// reconstruction by inverse Haar transform
	for (int scale=s-1; scale>=0; scale--)
	{
		OrthogonalWaveletTransform<DATATYPE>::idwt3D( \
			data, coef[scale+1], coef[s+scale+1], \
			coef[2*s+scale+1], coef[3*s+scale+1], coef[4*s+scale+1], \
			coef[5*s+scale+1], coef[6*s+scale+1], \
			rfilterh, rfilterg, coef[0]);
		data = coef[0];		
	}
	
	if (coef != NULL) { delete[] coef; coef = NULL; }
}

template<class DATATYPE>
void FiszTransform<DATATYPE>::ifisz3D (to_array<DATATYPE, true> &data)
{
	int len1 = data.nx(), len2 = data.ny(), len3 = data.nz();
	int len = MIN(MIN(len1, len2), len3);
	if (len <= 1) return;
	else if (!is_power_of_2(len1))
		throw DataSizeException(len1, "data size not of power of two");
	else if (!is_power_of_2(len2))
		throw DataSizeException(len2, "data size not of power of two");
	else if (!is_power_of_2(len3))
		throw DataSizeException(len3, "data size not of power of two");

	int s = iilog2(len);
	to_array<DATATYPE, true> *coef = new to_array<DATATYPE, true>[7*s + 1];
	
	dblarray dfilterh(2), dfilterg(2), rfilterh(2), rfilterg(2);
	dfilterh(0) = .5; dfilterh(1) = .5; dfilterg(0) = -.5; dfilterg(1) = .5;
    rfilterh(0) = 1; rfilterh(1) = 1; rfilterg(0) = 1; rfilterg(1) = -1;

    // decomposition by Haar transform
	for (int scale=0; scale<s; scale++)
	{
		OrthogonalWaveletTransform<DATATYPE>::dwt3D(data, dfilterh, dfilterg, \
			coef[0], coef[scale+1], coef[s+scale+1], \
			coef[2*s+scale+1], coef[3*s+scale+1], coef[4*s+scale+1], \
			coef[5*s+scale+1], coef[6*s+scale+1]);
		data = coef[0];
	}
	
	// modification of detail coefficients and reconstruction by inverse Haar transform
	for (int scale=s-1; scale>=0; scale--)
	{
		int lx = coef[0].nx(), ly = coef[0].ny(), lz = coef[0].nz();
		for (int x=0; x<lx; x++)
		for (int y=0; y<ly; y++)
		for (int z=0; z<lz; z++)
			coef[0](x, y, z) = psqrt(coef[0](x, y, z));
			
		for (int i=0; i<7; i++)
			coef[i*s+scale+1] *= coef[0];
		
		OrthogonalWaveletTransform<DATATYPE>::idwt3D( \
			data, coef[scale+1], coef[s+scale+1], \
			coef[2*s+scale+1], coef[3*s+scale+1], coef[4*s+scale+1], \
			coef[5*s+scale+1], coef[6*s+scale+1], \
			rfilterh, rfilterg, coef[0]);
		data = coef[0];		
	}
	
	if (coef != NULL) { delete[] coef; coef = NULL; }
}

#endif
