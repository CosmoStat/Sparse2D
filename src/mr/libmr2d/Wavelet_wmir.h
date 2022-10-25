/*
 * Filename : Wavelet.h
 * 
 * Class Description
 * 1. OrthogonalWaveletTransform : Wavelet filters; DWT1D, DWT2D and DWT3D; IDWT1D, IDWT2D and IDWT3D
 * 2. WaveletShrinkage : Noise estimation; Hard/Soft wavelet thresholding; wavelet hypothesis tests
 * 3. MWIRWaveletTransform : wavelet transform in using the code of Jean-Luc Starck 1D/2D/3D decimated/undecimated
 */
#ifndef _WAVELET_H
#define _WAVELET_H
#include <string>
#include <math.h>
#include "cdflib.h"
#include "Array.h"
#include "ImLib_mwir.h"
#include "MR_HaarPoisson.h"
using namespace std;
extern bool VERBOSE;

template <class DATATYPE> 
class OrthogonalWaveletTransform
{
	protected:
		// given the data (1D) length and the filter length,
		// it calculates the decomposed data length
		static int getDecompResultLength(int datal, int filterl)
		{
			if ((datal % 2 == 0) && (filterl % 2 == 0))
				return MIN(datal, (datal+filterl-2)/2);
			else 
				return MIN(datal, (datal+filterl)/2);
		}
		
		// given the data (1D) length and the filter length,
		// it calculates the reconstructed data length
		static int getReconsResultLength(int datal, int filterl)
		{
			return MAX(datal, (datal-(filterl-1)/2)*2);
		}
		
		// given a filter, it calculates its reversion
		static void revertFilter (dblarray &filter, dblarray &revfilter)
		{
			int n = filter.n_elem();
			revfilter.resize(n);
			for (int x=0; x<n; x++)
				revfilter(x) = filter(n-1-x);
		}
		
		// given a filter, it calculates the quadratic
		// mirror filter of the given one : qmfilter[n] = (-1)^n filter[1-n]
		static void quadraticMirrorFilter (dblarray &filter, dblarray &qmfilter)
		{
			int n = filter.n_elem();
			qmfilter.resize(n);
			int s = -1;
			for (int x=n-1; x>=0; x--)
			{
				qmfilter(x) = s * filter(n-1-x);
				s = -s;
			}
		}
		
		// (I)DWT 1D along a certain axis
		// only used by (I)DWT3D
		// axis = 0, 1, 2 : X, Y, Z
		static void transformXYZ (to_array<DATATYPE, true> &data, dblarray &dfilterh, dblarray &dfilterg, to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int axis, type_border BORDERTYPE);
		static void reconstructionXYZ (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, dblarray &rfilterh, dblarray &rfilterg, to_array<DATATYPE, true> &data, int axis);
		
	public:
		// orthogonal wavelet types
		static const int HAAR=0, DAUB4=4, SYML4=16;
		
		// resize the data to the same size of the reference
		static void resizeData (to_array<DATATYPE, true> &data, to_array<DATATYPE, true> &ref, type_border BORDERTYPE=I_ZERO);
		// resize the data to [nx], [nx, ny] or [nx, ny, nz] according to the dimension of the data
		static void resizeData (to_array<DATATYPE, true> &data, int nx, int ny, int nz, type_border BORDERTYPE=I_ZERO);
		
		// get the decomposition filters, i.e. \bar{h}, \bar{g}
		// where h is considered to start from the origin
		static void getWaveletDecompFilter (int filterName, dblarray &dfilterh, dblarray &dfilterg);
		// get the reconstruction filters, i.e. \tilde{h} (=h), \tilde{g} (=g)
		static void getWaveletReconsFilter (int filterName, dblarray &rfilterh, dblarray &rfilterg);
		
		// DWT 1D
		static void dwt1D (to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE=I_MIRROR);
		static void dwt1D (to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE=I_MIRROR);
		
		// IDWT 1D 
		static void idwt1D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data);
		static void idwt1D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			int filterName, \
			to_array<DATATYPE, true> &data);			

		// DWT 2D
		static void dwt2D (to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE=I_MIRROR);
		static void dwt2D (to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE=I_MIRROR);
					
		// IDWT 2D 
		static void idwt2D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data);
		static void idwt2D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			int filterName, \
			to_array<DATATYPE, true> &data);

		// DWT 3D
		static void dwt3D (to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			type_border BORDERTYPE=I_MIRROR);
		static void dwt3D (to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			type_border BORDERTYPE=I_MIRROR);
					
		// IDWT 3D 
		static void idwt3D ( \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data);
		static void idwt3D ( \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			int filterName, \
			to_array<DATATYPE, true> &data);
};

template <class DATATYPE>
const int OrthogonalWaveletTransform<DATATYPE>::HAAR;

template <class DATATYPE>
const int OrthogonalWaveletTransform<DATATYPE>::DAUB4;

template <class DATATYPE>
const int OrthogonalWaveletTransform<DATATYPE>::SYML4;

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::resizeData(to_array<DATATYPE, true> &data, to_array<DATATYPE, true> &ref, type_border BORDERTYPE)
{
	resizeData(data, ref.nx(), ref.ny(), ref.nz(), BORDERTYPE);
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::resizeData(to_array<DATATYPE, true> &data, int nx, int ny, int nz, type_border BORDERTYPE)
{
	int dnx = data.nx(), dny = data.ny(), dnz = data.nz();
	if ((nx == dnx) && (ny == dny) && (nz == dnz)) return;
	int dim = data.naxis();

	to_array<DATATYPE, true> *temp = new to_array<DATATYPE, true>(dnx, dny, dnz);
	*temp = data;
	if (dim == 1)
	{
		data.resize(nx);
		for (int x=0; x<nx; x++)
			data(x) = (*temp)(x, BORDERTYPE);
	}
	else if (dim == 2)
	{
		data.resize(nx, ny);
		for (int x=0; x<nx; x++)
			for (int y=0; y<ny; y++)
				data(x, y) = (*temp)(x, y, BORDERTYPE);
	}
	else if (dim == 3)
	{
		data.resize(nx, ny, nz);
		for (int x=0; x<nx; x++)
			for (int y=0; y<ny; y++)
				for (int z=0; z<nz; z++)
				data(x, y, z) = (*temp)(x, y, z, BORDERTYPE);		
	}
	else throw DataSizeException(dim, "unknown size of dimension");
	
	if (temp != NULL) { delete temp; temp = NULL; }
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::getWaveletDecompFilter (int filterName, dblarray &dfilterh, dblarray &dfilterg)
{
	dblarray temp;
	
	switch (filterName)
	{
		case OrthogonalWaveletTransform<DATATYPE>::HAAR:
			dfilterh.resize(2);
			dfilterh(0) = 1 / sqrt(2.);	dfilterh(1) = dfilterh(0);

			revertFilter(dfilterh, dfilterg);
			quadraticMirrorFilter(dfilterg, temp);
			revertFilter(temp, dfilterg);
			break;

		case OrthogonalWaveletTransform<DATATYPE>::DAUB4:
			dfilterh.resize(8);
			dfilterh(0) = -0.01059740178500; dfilterh(1) = 0.03288301166698;
			dfilterh(2) = 0.03084138183599; dfilterh(3) =	-0.18703481171888;
			dfilterh(4) = -0.02798376941698; dfilterh(5) = 0.63088076792959;
			dfilterh(6) = 0.71484657055254; dfilterh(7) = 0.23037781330886;

			revertFilter(dfilterh, dfilterg);
			quadraticMirrorFilter(dfilterg, temp);
			revertFilter(temp, dfilterg);			
			break;

		case OrthogonalWaveletTransform<DATATYPE>::SYML4:
			dfilterh.resize(8); 
			dfilterh(0) = -0.07576571478927; dfilterh(1) = -0.02963552764600;
			dfilterh(2) = 0.49761866763202; dfilterh(3) = 0.80373875180592;
			dfilterh(4) = 0.29785779560528; dfilterh(5) = -0.09921954357685;
			dfilterh(6) = -0.01260396726204; dfilterh(7) = 0.03222310060404;
			
			revertFilter(dfilterh, dfilterg);
			quadraticMirrorFilter(dfilterg, temp);
			revertFilter(temp, dfilterg);			
			break;
		default:
			throw DataIDException<int>(filterName, "Unknown filter name");
	}
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::getWaveletReconsFilter (int filterName, dblarray &rfilterh,  dblarray &rfilterg)
{
	dblarray temph, tempg;
	
	getWaveletDecompFilter (filterName, temph, tempg);
	revertFilter(temph, rfilterh);
	revertFilter(tempg, rfilterg);
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::transformXYZ (to_array<DATATYPE, true> &data, dblarray &dfilterh, dblarray &dfilterg, to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int axis, type_border BORDERTYPE)
{
	to_array<DATATYPE, true> *line = NULL, *linea = NULL, *lined = NULL;
	int lenx, leny, lenz, dlenx, dleny, dlenz;
	
	if (axis == 0)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
		dlenx = getDecompResultLength(lenx, dfilterh.n_elem());
		ca.resize(dlenx, leny, lenz);
		cd.resize(dlenx, leny, lenz);
		line = new to_array<DATATYPE, true>(lenx);
		linea = new to_array<DATATYPE, true>(dlenx);
		lined = new to_array<DATATYPE, true>(dlenx);
		
		for (int z=0; z<lenz; z++)
		for (int y=0; y<leny; y++)
		{
			for (int x=0; x<lenx; x++)
				(*line)(x) = data(x, y, z);
			dwt1D(*line, dfilterh, dfilterg, *linea, *lined, BORDERTYPE);
			for (int x=0; x<dlenx; x++)
			{
				ca(x, y, z) = (*linea)(x);
				cd(x, y, z) = (*lined)(x);
			}
		}
	}	
	else if (axis == 1)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
		dleny = getDecompResultLength(leny, dfilterh.n_elem());
		ca.resize(lenx, dleny, lenz);
		cd.resize(lenx, dleny, lenz);
		line = new to_array<DATATYPE, true>(leny);
		linea = new to_array<DATATYPE, true>(dleny);
		lined = new to_array<DATATYPE, true>(dleny);

		for (int z=0; z<lenz; z++)
		for (int x=0; x<lenx; x++)
		{
			for (int y=0; y<leny; y++)
				(*line)(y) = data(x, y, z);
			dwt1D(*line, dfilterh, dfilterg, *linea, *lined, BORDERTYPE);
			for (int y=0; y<dleny; y++)
			{
				ca(x, y, z) = (*linea)(y);
				cd(x, y, z) = (*lined)(y);
			}
		}
	}
	else if (axis == 2)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
		dlenz = getDecompResultLength(lenz, dfilterh.n_elem());	
		ca.resize(lenx, leny, dlenz);
		cd.resize(lenx, leny, dlenz);
		line = new to_array<DATATYPE, true>(lenz);
		linea = new to_array<DATATYPE, true>(dlenz);
		lined = new to_array<DATATYPE, true>(dlenz);

		for (int y=0; y<leny; y++)
		for (int x=0; x<lenx; x++)
		{
			for (int z=0; z<lenz; z++)
				(*line)(z) = data(x, y, z);
			dwt1D(*line, dfilterh, dfilterg, *linea, *lined, BORDERTYPE);
			for (int z=0; z<dlenz; z++)
			{
				ca(x, y, z) = (*linea)(z);
				cd(x, y, z) = (*lined)(z);
			}
		}
	}
	else throw DataException("unknown axis");
	
	if (line != NULL) { delete line; line = NULL; }
	if (linea != NULL) { delete linea; linea = NULL; }
	if (lined != NULL) { delete lined; lined = NULL; }
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::reconstructionXYZ (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, dblarray &rfilterh, dblarray &rfilterg, to_array<DATATYPE, true> &data, int axis)
{
	to_array<DATATYPE, true> *line = NULL, *linea = NULL, *lined = NULL;
	int lenx, leny, lenz, rlenx, rleny, rlenz;
	
	if (axis == 0)
	{
		lenx = ca.nx(); leny = ca.ny(); lenz = ca.nz();
		rlenx = getReconsResultLength(lenx, rfilterh.n_elem());
		line = new to_array<DATATYPE, true>(rlenx);
		linea = new to_array<DATATYPE, true>(lenx);
		lined = new to_array<DATATYPE, true>(lenx);
		data.resize(rlenx, leny, lenz);
		
		for (int z=0; z<lenz; z++)
		for (int y=0; y<leny; y++)
		{
			for (int x=0; x<lenx; x++)
			{
				(*linea)(x) = ca(x, y, z);
				(*lined)(x) = cd(x, y, z);
			}
			idwt1D(*linea, *lined, rfilterh, rfilterg, *line);
			for (int x=0; x<rlenx; x++)
				data(x, y, z) = (*line)(x);
		}
	}
	else if (axis == 1)
	{
		lenx = ca.nx(); leny = ca.ny(); lenz = ca.nz();
		rleny = getReconsResultLength(leny, rfilterh.n_elem());
		line = new to_array<DATATYPE, true>(rleny);
		linea = new to_array<DATATYPE, true>(leny);
		lined = new to_array<DATATYPE, true>(leny);
		data.resize(lenx, rleny, lenz);
		
		for (int z=0; z<lenz; z++)
		for (int x=0; x<lenx; x++)
		{
			for (int y=0; y<leny; y++)
			{
				(*linea)(y) = ca(x, y, z);
				(*lined)(y) = cd(x, y, z);
			}
			idwt1D(*linea, *lined, rfilterh, rfilterg, *line);
			for (int y=0; y<rleny; y++)
				data(x, y, z) = (*line)(y);
		}
	}
	else if (axis == 2)
	{
		lenx = ca.nx(); leny = ca.ny(); lenz = ca.nz();
		rlenz = getReconsResultLength(lenz, rfilterh.n_elem());
		line = new to_array<DATATYPE, true>(rlenz);
		linea = new to_array<DATATYPE, true>(lenz);
		lined = new to_array<DATATYPE, true>(lenz);
		data.resize(lenx, leny, rlenz);
		
		for (int y=0; y<leny; y++)
		for (int x=0; x<lenx; x++)
		{
			for (int z=0; z<lenz; z++)
			{
				(*linea)(z) = ca(x, y, z);
				(*lined)(z) = cd(x, y, z);
			}
			idwt1D(*linea, *lined, rfilterh, rfilterg, *line);
			for (int z=0; z<rlenz; z++)
				data(x, y, z) = (*line)(z);
		}
	}
	else throw DataException("unknown axis");
	
	if (line != NULL) { delete line; line = NULL; }
	if (linea != NULL) { delete linea; linea = NULL; }
	if (lined != NULL) { delete lined; lined = NULL; }
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::dwt1D (to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE)
{
	int datalen = data.nx();
	int hlen = dfilterh.n_elem();
	int glen = dfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	int reslen = getDecompResultLength(datalen, hlen);
	
	ca.resize(reslen); cd.resize(reslen);
	int offset = 0;
	int leftoffseth = (hlen % 2 == 0) ? (hlen-2) : (hlen-1);
	int leftoffsetg = glen-2;
	double coeffa, coeffd;
	for (int x=0; x<reslen; x++)
	{
		coeffa = 0, coeffd = 0;
		for (int xx=hlen-1; xx>=0; xx--)
		{
			coeffa += dfilterh(xx) * data(offset-leftoffseth+hlen-1-xx, BORDERTYPE);
			coeffd += dfilterg(xx) * data(offset-leftoffsetg+hlen-1-xx, BORDERTYPE);
		}
		ca(x) = coeffa;
		cd(x) = coeffd;	
		offset += 2;
	}
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::dwt1D (to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE)
{
	dblarray filterh, filterg;
	getWaveletDecompFilter(filterName, filterh, filterg);
	dwt1D (data, filterh, filterg, ca, cd, BORDERTYPE);
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::idwt1D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data)
{
	int calen = ca.nx(), cdlen = cd.nx();
	if (calen != cdlen)
		throw DataSizeException(calen, "the two given wavelet coefficients vectors are not of the same size");
	int hlen = rfilterh.n_elem(), glen = rfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	int reslen = getReconsResultLength(calen, hlen);
	data.resize(reslen);
	
	int hoffset = 1 - hlen % 2;
	int goffset = 1;
	double coeffa, coeffd;
	int aleftoffset = 0, dleftoffset = 0;
	for (int x=0; x<reslen; x++)
	{
		coeffa = 0; coeffd = 0;
		
		for (int xx=hlen-1-hoffset,dataindex=0; xx>=0; xx-=2,dataindex++)
			coeffa += rfilterh(xx) * ca(aleftoffset+dataindex);
		for (int xx=glen-1-goffset,dataindex=0; xx>=0; xx-=2,dataindex++)
			coeffd += rfilterg(xx) * cd(dleftoffset+dataindex);
		
		data(x) = coeffa + coeffd;
		if (hoffset % 2 == 0) aleftoffset++;
		if (goffset % 2 == 0) dleftoffset++;
		hoffset = 1 - hoffset;
		goffset = 1 - goffset;
	}
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::idwt1D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, \
			int filterName, \
			to_array<DATATYPE, true> &data)
{
	dblarray filterh, filterg;
	getWaveletReconsFilter(filterName, filterh, filterg);
	idwt1D (ca, cd, filterh, filterg, data);	
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::dwt2D (to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE)
{
	int datalenx = data.nx(), dataleny = data.ny();
	int hlen = dfilterh.n_elem(), glen = dfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	int reslenx = getDecompResultLength(datalenx, hlen);
	int resleny = getDecompResultLength(dataleny, hlen);
	ca.resize(reslenx, resleny); 	ch.resize(reslenx, resleny);
	cv.resize(reslenx, resleny); 	cd.resize(reslenx, resleny);

	to_array<DATATYPE, true> *temph = new to_array<DATATYPE, true>(reslenx, dataleny);
	to_array<DATATYPE, true> *tempg = new to_array<DATATYPE, true>(reslenx, dataleny);
	to_array<DATATYPE, true> *line = new to_array<DATATYPE, true>(datalenx);
	to_array<DATATYPE, true> *atempline = new to_array<DATATYPE, true>(reslenx);
	to_array<DATATYPE, true> *dtempline = new to_array<DATATYPE, true>(reslenx);
	for (int y=0; y<dataleny; y++)
	{
		for (int x=0; x<datalenx; x++) (*line)(x) = data(x,y);
		dwt1D(*line, dfilterh, dfilterg, *atempline, *dtempline, BORDERTYPE);
		for (int x=0; x<reslenx; x++)
		{
			(*temph)(x, y) = (*atempline)(x);
			(*tempg)(x, y) = (*dtempline)(x);
		}
	}
	if (line != NULL) { delete line; line = NULL; }
	if (atempline != NULL) { delete atempline; atempline = NULL; }
	if (dtempline != NULL) { delete dtempline; dtempline = NULL; }
	
	to_array<DATATYPE, true> *col1 = new to_array<DATATYPE, true>(dataleny);
	to_array<DATATYPE, true> *col2 = new to_array<DATATYPE, true>(dataleny);
	to_array<DATATYPE, true> *atempcol = new to_array<DATATYPE, true>(resleny);
	to_array<DATATYPE, true> *dtempcol = new to_array<DATATYPE, true>(resleny);
	for (int x=0; x<reslenx; x++)
	{
		for (int y=0; y<dataleny; y++)
		{ 
			(*col1)(y) = (*temph)(x,y);
			(*col2)(y) = (*tempg)(x,y);
		}
		dwt1D(*col1, dfilterh, dfilterg, *atempcol, *dtempcol, BORDERTYPE);		
		for (int y=0; y<resleny; y++)
		{
			ca(x,y) = (*atempcol)(y);
			ch(x,y) = (*dtempcol)(y);
		}
		dwt1D(*col2, dfilterh, dfilterg, *atempcol, *dtempcol, BORDERTYPE);		
		for (int y=0; y<resleny; y++)
		{
			cv(x,y) = (*atempcol)(y);
			cd(x,y) = (*dtempcol)(y);
		}
	}
	if (temph != NULL) { delete temph; temph = NULL; }
	if (tempg != NULL) { delete tempg; tempg = NULL; }
	if (col1 != NULL) { delete col1; col1 = NULL; }
	if (col2 != NULL) { delete col2; col2 = NULL; }
	if (atempcol != NULL) { delete atempcol; atempcol = NULL; }
	if (dtempcol != NULL) { delete dtempcol; dtempcol = NULL; }	
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::dwt2D (to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			type_border BORDERTYPE)
{
	dblarray filterh, filterg;
	getWaveletDecompFilter(filterName, filterh, filterg);
	dwt2D (data, filterh, filterg, ca, ch, cv, cd, BORDERTYPE);
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::idwt2D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data)
{
	int calenx = ca.nx(), caleny = ca.ny();
	int chlenx = ch.nx(), chleny = ch.ny();
	int cvlenx = cv.nx(), cvleny = cv.ny();
	int cdlenx = cd.nx(), cdleny = cd.ny();
	if ((calenx != chlenx) || (calenx != cvlenx) || (calenx != cdlenx))
		throw DataSizeException(calenx, "the four given wavelet coefficients matrix are not of the same size");
	if ((caleny != chleny) || (caleny != cvleny) || (caleny != cdleny))
		throw DataSizeException(caleny, "the four given wavelet coefficients matrix are not of the same size");
	int hlen = rfilterh.n_elem(), glen = rfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	int reslenx = getReconsResultLength(calenx, hlen);
	int resleny = getReconsResultLength(caleny, hlen);
	data.resize(reslenx, resleny);

	to_array<DATATYPE, true> *temph = new to_array<DATATYPE, true>(calenx, resleny);
	to_array<DATATYPE, true> *tempg = new to_array<DATATYPE, true>(calenx, resleny);
	to_array<DATATYPE, true> *cola = new to_array<DATATYPE, true>(caleny);
	to_array<DATATYPE, true> *colh = new to_array<DATATYPE, true>(chleny);
	to_array<DATATYPE, true> *colv = new to_array<DATATYPE, true>(cvleny);
	to_array<DATATYPE, true> *cold = new to_array<DATATYPE, true>(cdleny);
	to_array<DATATYPE, true> *coltemph = new to_array<DATATYPE, true>(resleny);	
	to_array<DATATYPE, true> *coltempg = new to_array<DATATYPE, true>(resleny);	
	for (int x=0; x<calenx; x++)
	{
		for (int y=0; y<caleny; y++)
		{
			(*cola)(y) = ca(x,y);
			(*colh)(y) = ch(x,y);
			(*colv)(y) = cv(x,y);
			(*cold)(y) = cd(x,y);
		}
		idwt1D(*cola, *colh, rfilterh, rfilterg, *coltemph);
		idwt1D(*colv, *cold, rfilterh, rfilterg, *coltempg);
		for (int y=0; y<resleny; y++)
		{
			(*temph)(x, y) = (*coltemph)(y);
			(*tempg)(x, y) = (*coltempg)(y);
		}
	}
	if (cola != NULL) { delete cola; cola = NULL; }
	if (colh != NULL) { delete colh; colh = NULL; }
	if (colv != NULL) { delete colv; colv = NULL; }
	if (cold != NULL) { delete cold; cold = NULL; }
	if (coltemph != NULL) { delete coltemph; coltemph = NULL; }
	if (coltempg != NULL) { delete coltempg; coltempg = NULL; }
	
	to_array<DATATYPE, true> *lineh = new to_array<DATATYPE, true>(calenx);	
	to_array<DATATYPE, true> *lineg = new to_array<DATATYPE, true>(calenx);	
	to_array<DATATYPE, true> *line = new to_array<DATATYPE, true>(reslenx);	
	for (int y=0; y<resleny; y++)
	{
		for (int x=0; x<calenx; x++)
		{
			(*lineh)(x) = (*temph)(x,y);
			(*lineg)(x) = (*tempg)(x,y);
		}
		idwt1D(*lineh, *lineg, rfilterh, rfilterg, *line);
		for (int x=0; x<reslenx; x++)
			data(x, y) = (*line)(x);
	}
	if (temph != NULL) { delete temph; temph = NULL; }
	if (tempg != NULL) { delete tempg; tempg = NULL; }
	if (lineh != NULL) { delete lineh; lineh = NULL; }
	if (lineg != NULL) { delete lineg; lineg = NULL; }
	if (line != NULL) { delete line; line = NULL; }
}

template <class DATATYPE> 
void OrthogonalWaveletTransform<DATATYPE>::idwt2D (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &ch, to_array<DATATYPE, true> &cv, to_array<DATATYPE, true> &cd, \
			int filterName, \
			to_array<DATATYPE, true> &data)
{
	dblarray filterh, filterg;
	getWaveletReconsFilter(filterName, filterh, filterg);
	idwt2D (ca, ch, cv, cd, filterh, filterg, data);
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::dwt3D ( \
			to_array<DATATYPE, true> &data, \
			dblarray &dfilterh, dblarray &dfilterg, \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			type_border BORDERTYPE)
{
	int hlen = dfilterh.n_elem(), glen = dfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	int datalenx = data.nx(), dataleny = data.ny(), datalenz = data.nz();
	int reslenx = getDecompResultLength(datalenx, hlen);
	int resleny = getDecompResultLength(dataleny, hlen);
	int reslenz = getDecompResultLength(datalenz, hlen);
	ca.resize(reslenx, resleny, reslenz); chhg.resize(reslenx, resleny, reslenz);
	chgh.resize(reslenx, resleny, reslenz); chgg.resize(reslenx, resleny, reslenz);
	cghh.resize(reslenx, resleny, reslenz); cghg.resize(reslenx, resleny, reslenz);
	cggh.resize(reslenx, resleny, reslenz); cggg.resize(reslenx, resleny, reslenz);

	// X
	to_array<DATATYPE, true> *temph = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *tempg = new to_array<DATATYPE, true>;
	transformXYZ(data, dfilterh, dfilterg, *temph, *tempg, 0, BORDERTYPE);
	
	// Y
	to_array<DATATYPE, true> *temphh = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *temphg = new to_array<DATATYPE, true>;
	transformXYZ(*temph, dfilterh, dfilterg, *temphh, *temphg, 1, BORDERTYPE);
	if (temph != NULL) { delete temph; temph = NULL; }
	
	to_array<DATATYPE, true> *tempgh = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *tempgg = new to_array<DATATYPE, true>;
	transformXYZ(*tempg, dfilterh, dfilterg, *tempgh, *tempgg, 1, BORDERTYPE);
	if (tempg != NULL) { delete tempg; tempg = NULL; }
	
	// Z
	transformXYZ(*temphh, dfilterh, dfilterg, ca, chhg, 2, BORDERTYPE);
	transformXYZ(*temphg, dfilterh, dfilterg, chgh, chgg, 2, BORDERTYPE);
	transformXYZ(*tempgh, dfilterh, dfilterg, cghh, cghg, 2, BORDERTYPE);
	transformXYZ(*tempgg, dfilterh, dfilterg, cggh, cggg, 2, BORDERTYPE);
	
	if (temphh != NULL) { delete temphh; temphh = NULL; }
	if (temphg != NULL) { delete temphg; temphg = NULL; }
	if (tempgh != NULL) { delete tempgh; tempgh = NULL; }
	if (tempgg != NULL) { delete tempgg; tempgg = NULL; }	
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::dwt3D ( \
			to_array<DATATYPE, true> &data, \
			int filterName, \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			type_border BORDERTYPE)
{
	dblarray filterh, filterg;
	getWaveletDecompFilter(filterName, filterh, filterg);
	dwt3D (data, filterh, filterg, ca, chhg, chgh, chgg, cghh, cghg, cggh, cggg, BORDERTYPE);
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::idwt3D ( \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			dblarray &rfilterh, dblarray &rfilterg, \
			to_array<DATATYPE, true> &data)
{
	int hlen = rfilterh.n_elem(), glen = rfilterg.n_elem();
	if (hlen != glen)
		throw DataSizeException(hlen, "the two given orthogonal wavelet filters are not of the same size");
	if ((ca.nx() != chhg.nx()) || (ca.nx() != chgh.nx()) || \
		(ca.nx() != chgh.nx()) || (ca.nx() != chgg.nx()) || \
		(ca.nx() != cghh.nx()) || (ca.nx() != cghg.nx()) || \
		(ca.nx() != cggh.nx()) || (ca.nx() != cggg.nx()))
		throw DataSizeException(ca.nx(), "the eight given wavelet coefficients matrix are not of the same size");
	if ((ca.ny() != chhg.ny()) || (ca.ny() != chgh.ny()) || \
		(ca.ny() != chgh.ny()) || (ca.ny() != chgg.ny()) || \
		(ca.ny() != cghh.ny()) || (ca.ny() != cghg.ny()) || \
		(ca.ny() != cggh.ny()) || (ca.ny() != cggg.ny()))
		throw DataSizeException(ca.ny(), "the eight given wavelet coefficients matrix are not of the same size");
	if ((ca.nz() != chhg.nz()) || (ca.nz() != chgh.nz()) || \
		(ca.nz() != chgh.nz()) || (ca.nz() != chgg.nz()) || \
		(ca.nz() != cghh.nz()) || (ca.nz() != cghg.nz()) || \
		(ca.nz() != cggh.nz()) || (ca.nz() != cggg.nz()))
		throw DataSizeException(ca.nz(), "the eight given wavelet coefficients matrix are not of the same size");
	int lenx = ca.nx(), leny = ca.ny(), lenz = ca.nz();
	int reslenx = getReconsResultLength(lenx, hlen);
	int resleny = getReconsResultLength(leny, hlen);
	int reslenz = getReconsResultLength(lenz, hlen);
	data.resize(reslenx, resleny, reslenz);
	
	// Z
	to_array<DATATYPE, true> *temphh = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *temphg = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *tempgh = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *tempgg = new to_array<DATATYPE, true>;
	reconstructionXYZ(ca, chhg, rfilterh, rfilterg, *temphh, 2);
	reconstructionXYZ(chgh, chgg, rfilterh, rfilterg, *temphg, 2);
	reconstructionXYZ(cghh, cghg, rfilterh, rfilterg, *tempgh, 2);
	reconstructionXYZ(cggh, cggg, rfilterh, rfilterg, *tempgg, 2);
	
	// Y
	to_array<DATATYPE, true> *temph = new to_array<DATATYPE, true>;
	to_array<DATATYPE, true> *tempg = new to_array<DATATYPE, true>;
	reconstructionXYZ(*temphh, *temphg, rfilterh, rfilterg, *temph, 1);
	if (temphh != NULL) { delete temphh; temphh = NULL; }
	if (temphg != NULL) { delete temphg; temphg = NULL; }
	
	reconstructionXYZ(*tempgh, *tempgg, rfilterh, rfilterg, *tempg, 1);
	if (tempgh != NULL) { delete tempgh; tempgh = NULL; }
	if (tempgg != NULL) { delete tempgg; tempgg = NULL; }	

	// X
	reconstructionXYZ(*temph, *tempg, rfilterh, rfilterg, data, 0);
	if (tempg != NULL) { delete tempg; tempg = NULL; }
	if (temph != NULL) { delete temph; temph = NULL; }		
}

template <class DATATYPE>
void OrthogonalWaveletTransform<DATATYPE>::idwt3D ( \
			to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &chhg, \
			to_array<DATATYPE, true> &chgh, \
			to_array<DATATYPE, true> &chgg, \
			to_array<DATATYPE, true> &cghh, \
			to_array<DATATYPE, true> &cghg, \
			to_array<DATATYPE, true> &cggh, \
			to_array<DATATYPE, true> &cggg, \
			int filterName, \
			to_array<DATATYPE, true> &data)
{
	dblarray filterh, filterg;
	getWaveletReconsFilter(filterName, filterh, filterg);
	idwt3D (ca, chhg, chgh, chgg, cghh, cghg, cggh, cggg, filterh, filterg, data);
}

//------------------------------------------------------------------------------
template <class DATATYPE, class SUPTYPE> 
class WaveletShrinkage
{
    protected:
	double fisherApproxRootSelect (double m[4], double z, double lambda)
	{
	    double z2 = z * z;
	    double c = 1./8. * (z2 + 1 - 2*lambda + sqrt(1 + 12*lambda + 4*lambda*lambda + (12*lambda+2)*z2 + z2*z2));
	    c = MAX(c, 0);

	    for (int i=0; i<4; i++)
	      if (m[i] >= c) return m[i];
	    
	    return c;
	}

    // cumulative non-central chi-square distribution
    double cumNChi (double degfree, double noncen, double lim)
    {
        int status, which1 = 1;
	    double p, q, val, df, pnonc, bound;
	    
        df = MAX(1e-100, degfree);
        val = MAX(0, lim);
        pnonc = MAX(0, noncen);
	    cdfchn(&which1, &p, &q, &val, &df, &pnonc, &status, &bound);
	    
		if ((status == 1) || (status == 2)) p = bound;
		else if (status < 0) p = 0;
		
		return p;
    }
    
    double pvalueModelHaar (double lambda1, double lambda2, double param, double coefobs)
    {
        double p = 0;
        
        if (coefobs*param >= lambda1-lambda2)
        {
            if (coefobs > 0)
                p = cumNChi (2*param*coefobs, 2*lambda2, 2*lambda1);
            else
                p = 1 - cumNChi (-2*param*coefobs+2, 2*lambda1, 2*lambda2);
        }
        else
        {
            if (coefobs >= 0)
                p = 1 - cumNChi (2*param*coefobs+2, 2*lambda2, 2*lambda1);
            else
                p = cumNChi (-2*param*coefobs, 2*lambda1, 2*lambda2);
        }
        
        return p;
    }
    
	public:
		// gaussian noise standard deviation MAD estimation 
		// sigma = median (abs(waveletdata)) / 0.6745 
		double gaussianStdDevEstim (to_array<DATATYPE, true> &waveletData);

		// Denoise gaussian noise 
		// Hard threshold = k * sigma; sup use -sigma / sigma (>=0) to indicate the insignificant/significant coef. region 
		// N is the length of the original data
		// alpha is the preset p-value of the significance level when DEFAULT is not set
		// DEFAULT indicates to use the universal threshold
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double gaussHardThreshold (to_array<DATATYPE, true> &waveletData, int scale[], \
			double alpha, double sigma, to_array<SUPTYPE, true> *sup = NULL, int N = 1, bool DEFAULT = false);

		// Denoise gaussian noise 
		// Soft threshold = k * sigma; sup use -sigma / sigma (>=0) to indicate the insignificant/significant coef. region 
		// N is the length of the original data
		// alpha is the preset p-value of the significance level when DEFAULT is not set
		// DEFAULT indicates to use the universal threshold
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double gaussSoftThreshold (to_array<DATATYPE, true> &waveletData, int scale[], \
			double alpha, double sigma, to_array<SUPTYPE, true> *sup = NULL, int N = 1, bool DEFAULT = false);

		// Denoise gaussian noise 
		// FDR threshold; sup use -sigma / sigma (>=0) to indicate the insignificant/significant coef. region
		// alpha is the preset FDR
      	// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double gaussFDRThreshold (to_array<DATATYPE, true> &waveletData, \
			double sigma, double alpha, bool indep = true, to_array<SUPTYPE, true> *sup = NULL);

		// Direct Haar hard threshold in using Kolaczyk's approximation for Poisson noise
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// scale[] is the scale along each direction
		// N is the length of the original data
		// alpha is the preset p-value of the significance level when DEFAULT is not set
		// DEFAULT indicates to use the universal threshold
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double haarKolaThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup = NULL, int N = 1, bool DEFAULT = false);

		// Direct Haar hard threshold in using modified Kolaczyk's approximation (Fisher approx.) for Poisson noise
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// scale[] is the scale along each direction
		// N is the length of the original data
		// alpha is the preset p-value of the significance level when DEFAULT is not set
		// DEFAULT indicates to use the universal threshold
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double haarMKolaThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup = NULL, int N = 1, bool DEFAULT = false);

		// Direct Haar hard threshold based on Bijaoui-Jammal's thresholding table for Poisson noise
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// scale[] is the scale along each direction
		// alpha is the preset p-value of the significance level
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double haarBJThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup = NULL); 

		// Direct Haar hard threshold based on FDR (false discovery rate) for Poisson noise
		// for a band
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// scale[] is the scale along each direction
		// alpha is the preset FDR
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double haarFDRThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, bool indep = true, to_array<SUPTYPE, true> *sup = NULL);

		// Direct Haar hard threshold based on FDR (false discovery rate) for Poisson noise
		// for the bands of a scale
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// scale[] is the scale along each direction
		// alpha is the preset FDR
		// the thresholding p-value (of the distr. of the absolute value of the coef.) is returned 
		double haarFDRThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> cd[], int len, int scale[], double alpha, bool indep = true, to_array<SUPTYPE, true> *sup = NULL);

		// Model (prior) based Haar coef. threshold for Poisson noise
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// caModel is the approximation coefs. of a certain scale, say j, of the model data
		// cdModel is the detail coefs. of the model data
		// cd is the detail coefs to be denoised
		// scale[] is the scale along each direction
		// alpha is the preset p-value of the significance level when DEFAULT is not set
		// N is the length of the original data
		// DEFAULT indicates to use the universal threshold
		// the thresholding p-value is returned 
		double haarModelThreshold (to_array<DATATYPE, true> &caModel, \
               to_array<DATATYPE, true> &cdModel, to_array<DATATYPE, true> &cd, \
               int scale[], double alpha, to_array<SUPTYPE, true> *sup = NULL, \
               int N = 1, bool DEFAULT = false);
               
		// Model (prior) based Haar coef. FDR threshold for Poisson noise
		// The transform data is assumed to use L1-normalized Haar filter,
		// i.e. \bar{h} = [1/2, 1/2] \bar{g} = [-1/2, 1/2]; sup use -1 / 1 to indicate the insignificant/significant coef. region
		// caModel is the approximation coefs. of a certain scale, say j, of the model data
		// cdModel is the detail coefs. of the model data
		// cd is the detail coefs to be denoised
		// scale[] is the scale along each direction
		// alpha is the preset FDR
		// the thresholding p-value is returned 
        double haarModelFDRThreshold (to_array<DATATYPE, true> &caModel, \
               to_array<DATATYPE, true> &cdModel, to_array<DATATYPE, true> &cd, \
               int scale[], double alpha, bool indep = true, to_array<SUPTYPE, true> *sup = NULL);               
};

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::gaussianStdDevEstim (to_array<DATATYPE, true> &waveletData)
{
	double med = Utils<DATATYPE>::absMedian(waveletData);
	
	return (med / 0.6745);
}

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::gaussHardThreshold (to_array<DATATYPE, true> &waveletData, int scale[], double alpha, double sigma, to_array<SUPTYPE, true> *sup, int N, bool DEFAULT)
{
    int dim = waveletData.naxis();
	int nx = waveletData.nx(), ny = waveletData.ny(), nz = waveletData.nz();
	double coef;
	double cthresh, thresh;
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];

	if (DEFAULT)
	  cthresh = sqrt(2 * log(N - N / POW2(ss)));
	else
	  cthresh = Utils<double>::criticalThreshGauss(alpha / 2.);
	thresh = cthresh * sigma;

	if (dim == 1)
	  {
	    if (sup != NULL) sup->resize(nx);
	    for (int x=0; x<nx; x++)
	      {
	        if (sup != NULL) (*sup)(x) = (SUPTYPE)sigma;
         	coef = waveletData(x);
          	if (ABS(coef) <= thresh)
		    { 
		      waveletData(x) = 0;
		      if (sup != NULL) (*sup)(x) = (SUPTYPE)-sigma;
            }
	      }
	  }
	else if (dim == 2)
	  {
	    if (sup != NULL) sup->resize(nx, ny);
	    for (int x=0; x<nx; x++)
	    for (int y=0; y<ny; y++)
	      {
	        if (sup != NULL) (*sup)(x, y) = (SUPTYPE)sigma;
		    coef = waveletData(x, y);
		    if (ABS(coef) <= thresh)
		    { 
		      waveletData(x, y) = 0;
		      if (sup != NULL) (*sup)(x, y) = (SUPTYPE)-sigma;
		    }
	      }
	  }
	else // dim == 3
	  {
	    if (sup != NULL) sup->resize(nx, ny, nz);
	    for (int x=0; x<nx; x++)
	    for (int y=0; y<ny; y++)
	    for (int z=0; z<nz; z++)
	    {
	        if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) sigma;
		    coef = waveletData(x, y, z);
		    if (ABS(coef) <= thresh)
		    { 
		      waveletData(x, y, z) = 0;
		      if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)-sigma;
            }
	      }
	  }
	
	return (DEFAULT ? (2.*(1-Utils<double>::cumNormal(cthresh))) : alpha);	
}

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::gaussSoftThreshold (to_array<DATATYPE, true> &waveletData, int scale[], double alpha, double sigma, to_array<SUPTYPE, true> *sup, int N, bool DEFAULT)
{
    int dim = waveletData.naxis();
	int nx = waveletData.nx(), ny = waveletData.ny(), nz = waveletData.nz();
	double coef;
	double cthresh, thresh;
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];

	if (DEFAULT)
	  cthresh = sqrt(2 * log(N - N / POW2(ss)));
	else
	  cthresh = Utils<double>::criticalThreshGauss(alpha / 2.);
	thresh = cthresh * sigma;

	if (dim == 1)
	  {
	    if (sup != NULL) sup->resize(nx);
	    
	    for (int x=0; x<nx; x++)
	      {
	        if (sup != NULL) (*sup)(x) = (SUPTYPE)sigma;
		    coef = waveletData(x);
		    if (ABS(coef) <= thresh)
		    {
		     waveletData(x) = 0;
		     if (sup != NULL) (*sup)(x) = (SUPTYPE)-sigma;
		    }
		    else if (coef > 0)
		      waveletData(x) = coef - thresh;
            else
		      waveletData(x) = coef + thresh;
	      }	
	  }
	else if (dim == 2)
	  {
	    if (sup != NULL) sup->resize(nx, ny);
	    
	    for (int x=0; x<nx; x++)
	    for (int y=0; y<ny; y++)
	      {
	        if (sup != NULL) (*sup)(x, y) = (SUPTYPE)sigma;
		    coef = waveletData(x, y);
		    if (ABS(coef) <= thresh)
		    {
		      waveletData(x, y) = 0;
		      if (sup != NULL) (*sup)(x, y) = (SUPTYPE) -sigma;
		    }
		    else if (coef > 0)
		      waveletData(x, y) = coef - thresh;
            else
		      waveletData(x, y) = coef + thresh;
	      }	
	  }
	else // dim == 3 
	  {
	    if (sup != NULL) sup->resize(nx, ny, nz);
	    
	    for (int x=0; x<nx; x++)
	    for (int y=0; y<ny; y++)
	    for (int z=0; z<nz; z++)
	      {
	        if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)sigma;
		    coef = waveletData(x, y, z);
		    if (ABS(coef) <= thresh)
		    {
		     waveletData(x, y, z) = 0;
		     if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)-sigma;
	        }
		    else if (coef > 0)
		      waveletData(x, y, z) = coef - thresh;
            else
		      waveletData(x, y, z) = coef + thresh;
	      }	
	  }
	return (DEFAULT ? (2.*(1-Utils<double>::cumNormal(cthresh))) : alpha);
}
template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::gaussFDRThreshold (to_array<DATATYPE, true> &waveletData, \
					  double sigma, double alpha, bool indep, to_array<SUPTYPE, true> *sup)
{
  int dim = waveletData.naxis();
  int nx = waveletData.nx(), ny = waveletData.ny(), nz = waveletData.nz();
  double val, result, fdrp;

  if (dim == 1)
    {
      if (sup != NULL) sup->resize(nx);
      to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx);
      for (int x=0; x<nx; x++)
	  {
	  if (sigma < 1e-20) (*pvals)(x) = 0;
	  else
	    {
	      val = ABS(waveletData(x)) / sigma;
	      result = Utils<DATATYPE>::cumNormal (val);
	      (*pvals)(x) = 1 - result;
	    }
	  }
      fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);
      
      for (int x=0; x<nx; x++)
	  {
	  if ((*pvals)(x) > fdrp)
	    {
	      if (sup != NULL) (*sup)(x) = (SUPTYPE)-sigma;
	      waveletData(x) = 0;
	    }
	  else
	    {
	      if (sup != NULL) (*sup)(x) = (SUPTYPE)sigma;
	    }
	  }
      delete pvals; pvals = NULL;
    }
  else if (dim == 2)
    {
      if (sup != NULL) sup->resize(nx, ny);
      to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny);
      for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
	  {
	  if (sigma < 1e-20) (*pvals)(x, y) = 0;
	  else
	    {	
	      val = ABS(waveletData(x, y)) / sigma;
	      result = Utils<DATATYPE>::cumNormal (val);
	      (*pvals)(x, y) = 1 - result;
	    }
	  }
      fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);
      
      for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
	  {
	  if ((*pvals)(x, y) > fdrp)
	    {
	      if (sup != NULL) (*sup)(x, y) = (SUPTYPE)-sigma;
	      waveletData(x, y) = 0;
	    }
	  else
	    {
	      if (sup != NULL) (*sup)(x, y) = (SUPTYPE)sigma;
	    }
	  }
      delete pvals; pvals = NULL;
    }
  else if (dim == 3)
    {
      if (sup != NULL) sup->resize(nx, ny, nz);
      to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny, nz);
      for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
      for (int z=0; z<nz; z++)
	  {
	  if (sigma < 1e-20) (*pvals)(x, y, z) = 0;
	  else
	    {	
	      val = ABS(waveletData(x, y, z)) / sigma;
	      result = Utils<DATATYPE>::cumNormal (val);
	      (*pvals)(x, y, z) = 1 - result;
	    }
	  }
      fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);
      
      for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
      for (int z=0; z<nz; z++)
	  {
	  if ((*pvals)(x, y, z) > fdrp)
	    {
	      if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)-sigma;
	      waveletData(x, y, z) = 0;
	    }
	  else
	    {
	      if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)sigma;
	    }
	  }
      delete pvals; pvals = NULL;
    }
  else throw DataSizeException(dim, "unknown size of dimension");    
  
  return fdrp * 2.;  // p-value of the distr. of the absolute value of the coeff  
}

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarKolaThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup, int N, bool DEFAULT)
{
	int nx, ny, nz;
	double coeff, lambda, thresh, zalpha2;
	
	int dim = ca.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double C1, C3 = POW2(ss);
	if (DEFAULT)
	{
	  C1 = 2. * MAX(log(N / C3), 0);
	  zalpha2 = sqrt(C1);
	}
	else 
	{
	  zalpha2 = Utils<double>::criticalThreshGauss(alpha / 2.);
	  C1 = zalpha2 * zalpha2;
	}
	double C0 = POW2(-1.-ss);
	double C2 = C1 * C1;

	if (dim == 1)
	{
		nx = cd.nx();
		if (sup != NULL) sup->resize(nx);

		for (int x=0; x<nx; x++)
		{
			lambda = MAX(0, ca(x));
			thresh = C0 * (C1 + sqrt(C2 + 4*C3*lambda*C1));
			thresh = MIN(INFINITY, thresh);
			coeff = cd(x);

		        if (sup != NULL) (*sup)(x) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x) = 0;
			    if (sup != NULL) (*sup)(x) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 2)
	{
		nx = cd.nx(); ny = cd.ny();
		if (sup != NULL) sup->resize(nx, ny);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
			lambda = MAX(0, ca(x, y));
			thresh = C0 * (C1 + sqrt(C2 + 4*C3*lambda*C1));
			thresh = MIN(INFINITY, thresh);
			coeff = cd(x, y);

		        if (sup != NULL) (*sup)(x, y) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y) = 0;
			    if (sup != NULL) (*sup)(x, y) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 3)
	{
		nx = cd.nx(); ny = cd.ny(); nz = cd.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
			lambda = MAX(0, ca(x, y, z));
			thresh = C0 * (C1 + sqrt(C2 + 4*C3*lambda*C1));
			thresh = MIN(INFINITY, thresh);
			coeff = cd(x, y, z);

		        if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y, z) = 0;
			    if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) -1;
			}
		}
	}
	else throw DataSizeException(dim, "unknown size of dimension");
	return (DEFAULT ? (2.*(1-Utils<double>::cumNormal(zalpha2))) : alpha);	
}

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarMKolaThreshold (to_array<DATATYPE, true> &ca, \
			to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup, int N, bool DEFAULT)
{
	int nx, ny, nz, rn;
	double coeff, pcoef[5], pm[4], pmi[4], lambda, thresh, zalpha2;
	
	int dim = ca.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double lambdajparam = POW2(ss), z2, lambdaj, lambdaj2;
	if (DEFAULT)
	{
	  z2 = 2. * MAX(log(N / lambdajparam), 0);
	  zalpha2 = sqrt(z2);
	}
	else
	{
	  zalpha2 = Utils<double>::criticalThreshGauss(alpha / 2.); 
	  z2 = zalpha2 * zalpha2;
	}

	if (dim == 1)
	{
		nx = cd.nx();
		if (sup != NULL) sup->resize(nx);

		for (int x=0; x<nx; x++)
		{
			lambda = MAX(0, ca(x));
			lambdaj = lambdajparam * lambda;
			lambdaj2 = lambdaj * lambdaj;

			pcoef[0] = (z2 + 1) * (z2 + 1) * lambdaj2 - 4 * z2 * lambdaj2 * lambdaj;
			pcoef[1] = 2 * (z2 + 1) * (z2 + 1) * lambdaj -16 * z2 * lambdaj2 - 4 * lambdaj2;
			pcoef[2] = (z2 + 1) * (z2 + 1) - 20 * z2 * lambdaj - 12 * lambdaj + 4 * lambdaj2;
			pcoef[3] = 16 * lambdaj - 8 * (z2 + 1);
			pcoef[4] = 16;

			quartic(pcoef, pm, pmi, &rn);
			thresh = fisherApproxRootSelect(pm, zalpha2, lambdaj);
			thresh /= lambdajparam;
			thresh = MIN(INFINITY, thresh);

			coeff = cd(x);

		        if (sup != NULL) (*sup)(x) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x) = 0;
			    if (sup != NULL) (*sup)(x) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 2)
	{
		nx = cd.nx(); ny = cd.ny();
		if (sup != NULL) sup->resize(nx, ny);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
			lambda = MAX(0, ca(x, y));
			lambdaj = lambdajparam * lambda;
			lambdaj2 = lambdaj * lambdaj;

			pcoef[0] = (z2 + 1) * (z2 + 1) * lambdaj2 - 4 * z2 * lambdaj2 * lambdaj;
			pcoef[1] = 2 * (z2 + 1) * (z2 + 1) * lambdaj - 16 * z2 * lambdaj2 - 4 * lambdaj2;
			pcoef[2] = (z2 + 1) * (z2 + 1) - 20 * z2 * lambdaj - 12 * lambdaj + 4 * lambdaj2;
			pcoef[3] = 16 * lambdaj - 8 * (z2 + 1);
			pcoef[4] = 16;

			quartic(pcoef, pm, pmi, &rn);
			thresh = fisherApproxRootSelect(pm, zalpha2, lambdaj);
			thresh /= lambdajparam;
			thresh = MIN(INFINITY, thresh);

			coeff = cd(x, y);

		        if (sup != NULL) (*sup)(x, y) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y) = 0;
			    if (sup != NULL) (*sup)(x, y) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 3)
	{
		nx = cd.nx(); ny = cd.ny(); nz = cd.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
			lambda = MAX(0, ca(x, y, z));
			lambdaj = lambdajparam * lambda;
			lambdaj2 = lambdaj * lambdaj;

			pcoef[0] = (z2 + 1) * (z2 + 1) * lambdaj2 - 4 * z2 * lambdaj2 * lambdaj;
			pcoef[1] = 2 * (z2 + 1) * (z2 + 1) * lambdaj - 16 * z2 * lambdaj2 - 4 * lambdaj2;
			pcoef[2] = (z2 + 1) * (z2 + 1) - 20 * z2 * lambdaj - 12 * lambdaj + 4 * lambdaj2;
			pcoef[3] = 16 * lambdaj - 8 * (z2 + 1);
			pcoef[4] = 16;

			quartic(pcoef, pm, pmi, &rn);
			thresh = fisherApproxRootSelect(pm, zalpha2, lambdaj);
			thresh /= lambdajparam;
			thresh = MIN(INFINITY, thresh);

			coeff = cd(x, y, z);

		        if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y, z) = 0;
			    if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) -1;
			}
		}
	}
	else throw DataSizeException(dim, "unknown size of dimension");  
	return (DEFAULT ? (2.*(1-Utils<double>::cumNormal(zalpha2))) : alpha);
}

template <class DATATYPE, class SUPTYPE> 
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarBJThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, to_array<SUPTYPE, true> *sup)
{
	int nx, ny, nz;
	double coeff, lambda, thresh, alpha2 = alpha / 2.;
	
	double zalpha2 = Utils<double>::criticalThreshGauss(alpha2);
	int dim = ca.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double C = POW2(ss);

	if (dim == 1)
	{
		nx = cd.nx();
		if (sup != NULL) sup->resize(nx);

		for (int x=0; x<nx; x++)
		{
			lambda = MAX(0, ca(x));
			thresh = get_harr_poisson_threshold(lambda*C, alpha2);
			thresh /= C;
			thresh = MIN(INFINITY, thresh);
			
			coeff = cd(x);

		        if (sup != NULL) (*sup)(x) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x) = 0;
			    if (sup != NULL) (*sup)(x) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 2)
	{
		nx = cd.nx(); ny = cd.ny();
		if (sup != NULL) sup->resize(nx, ny);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
			lambda = MAX(0, ca(x, y));
			thresh = get_harr_poisson_threshold(lambda*C, alpha2);
			thresh /= C;
			thresh = MIN(INFINITY, thresh);

			coeff = cd(x, y);

		        if (sup != NULL) (*sup)(x, y) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y) = 0;
			    if (sup != NULL) (*sup)(x, y) = (SUPTYPE) -1;
			}
		}
	}
	else if (dim == 3)
	{
		nx = cd.nx(); ny = cd.ny(); nz = cd.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
			lambda = MAX(0, ca(x, y, z));
			thresh = get_harr_poisson_threshold(lambda*C, alpha2);
			thresh /= C;
			thresh = MIN(INFINITY, thresh);

			coeff = cd(x, y, z);

		        if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) 1;
			if (ABS(coeff) < thresh)
			{ 
			    cd(x, y, z) = 0;
			    if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE) -1;
			}
		}
	}
	else throw DataSizeException(dim, "unknown size of dimension");
	return alpha;
}

template <class DATATYPE, class SUPTYPE>
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarFDRThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> &cd, int scale[], double alpha, bool indep, to_array<SUPTYPE, true> *sup)
{
	int nx, ny, nz;
	double coeff, lambda, fdrp;
	double lambdaj, lambdaj2, mj, p, val, df, pnonc;
	
	int dim = ca.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double lambdajparam = POW2(ss);

	if (dim == 1)
	{
		nx = cd.nx();
		if (sup != NULL) sup->resize(nx);
	        to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx);
		for (int x=0; x<nx; x++)
		{
			lambda = MAX(0, ca(x)); lambdaj = lambdajparam * lambda;
			coeff = cd(x); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);
			(*pvals)(x) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		{
		  if ((*pvals)(x) > fdrp)
		    {
		      if (sup != NULL)
		      {
			(*sup)(x) = (SUPTYPE) -1;
		      }
		      cd(x) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			(*sup)(x) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 2)
	{
	        nx = cd.nx(); ny = cd.ny();
		if (sup != NULL) sup->resize(nx, ny);
	        to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
			lambda = MAX(0, ca(x, y)); lambdaj = lambdajparam * lambda;
			coeff = cd(x, y); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);

			(*pvals)(x, y) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
		  if ((*pvals)(x, y) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			(*sup)(x, y) = (SUPTYPE) -1;
		      }
		      cd(x, y) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			(*sup)(x, y) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 3)
	{
  	        nx = cd.nx(); ny = cd.ny(); nz = cd.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);
	        to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny, nz);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
			lambda = MAX(0, ca(x, y, z)); lambdaj = lambdajparam * lambda;
			coeff = cd(x, y, z); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);

			(*pvals)(x, y, z) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
		  if ((*pvals)(x, y, z) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			(*sup)(x, y, z) = (SUPTYPE) -1;
		      }
		      cd(x, y, z) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			(*sup)(x, y, z) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else throw DataSizeException(dim, "unknown size of dimension");    

	return fdrp * 2.;  // p-value of the distr. of the absolute value of the coeff
}

template <class DATATYPE, class SUPTYPE>
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarFDRThreshold (to_array<DATATYPE, true> &ca, to_array<DATATYPE, true> cd[], int len, int scale[], double alpha, bool indep, to_array<SUPTYPE, true> *sup)
{
	int nx, ny, nz;
	double coeff, lambda, fdrp;
	double lambdaj, lambdaj2, mj, p, val, df, pnonc;
	
	int dim = ca.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double lambdajparam = POW2(ss);

	if (dim == 1)
	{
	    nx = cd[0].nx();
	    if (sup != NULL) for (int l=0; l<len; l++) sup[l].resize(nx);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(len, nx);
		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		{
			lambda = MAX(0, ca(x)); lambdaj = lambdajparam * lambda;
			coeff = cd[l](x); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);

			(*pvals)(l, x) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		{
		  if ((*pvals)(l, x) > fdrp)
		    {
		      if (sup != NULL)
		      {
			sup[l](x) = (SUPTYPE) -1;
		      }
		      cd[l](x) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			sup[l](x) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 2)
	{
	    nx = cd[0].nx(); ny = cd[0].ny();
		if (sup != NULL) for (int l=0; l<len; l++) sup[l].resize(nx, ny);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(len, nx, ny);
		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
			lambda = MAX(0, ca(x, y)); lambdaj = lambdajparam * lambda;
			coeff = cd[l](x, y); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);

			(*pvals)(l, x, y) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
		  if ((*pvals)(l, x, y) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			sup[l](x, y) = (SUPTYPE) -1;
		      }
		      cd[l](x, y) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			sup[l](x, y) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 3)
	{
  	    nx = cd[0].nx(); ny = cd[0].ny(); nz = cd[0].nz();
		if (sup != NULL) for (int l=0; l<len; l++) sup[l].resize(nx, ny, nz);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(len*nx, ny, nz);
		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
			lambda = MAX(0, ca(x, y, z)); lambdaj = lambdajparam * lambda;
			coeff = cd[l](x, y, z); mj = lambdajparam * ABS(coeff);
			df = 2. * mj; pnonc = lambdaj; val = lambdaj;
            p = cumNChi (df, pnonc, val);

			(*pvals)(l*nx+x, y, z) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int l=0; l<len; l++)
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
		  if ((*pvals)(l*nx+x, y, z) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			sup[l](x, y, z) = (SUPTYPE) -1;
		      }
		      cd[l](x, y, z) = 0;
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			sup[l](x, y, z) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else throw DataSizeException(dim, "unknown size of dimension");    

	return fdrp * 2.;  // p-value of the distr. of the absolute value of the coeff  
}

template <class DATATYPE, class SUPTYPE>
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarModelThreshold (to_array<DATATYPE, true> &caModel, \
               to_array<DATATYPE, true> &cdModel, to_array<DATATYPE, true> &cd, \
               int scale[], double alpha, to_array<SUPTYPE, true> *sup, \
               int N, bool DEFAULT)
{
	int nx, ny, nz;
    int dim = caModel.naxis();
	double ss = 0, lambda1, lambda2, coefobs, p;
	for (int si=0; si<dim; si++) ss += scale[si];
    double param = POW2(ss);

   	if (DEFAULT)
	   alpha = 2. * (1 - Utils<double>::cumNormal (sqrt(2. * MAX(log(N / param), 0))));

	if (dim == 1)
	{
        nx = caModel.nx();
		if (sup != NULL) sup->resize(nx);
		for (int x=0; x<nx; x++)
		{
            coefobs = cd(x);
            lambda1 = (caModel(x) + cdModel(x)) * param / 2;
            lambda2 = (caModel(x) - cdModel(x)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
            if (p > alpha / 2) 
            {
               cd(x) = cdModel(x);
               if (sup != NULL) (*sup)(x) = (SUPTYPE)-1;
            }
            else
            {
               if (sup != NULL) (*sup)(x) = (SUPTYPE)1;
            }
        }
    }
    else if (dim == 2)
    {
        nx = caModel.nx(); 
        ny = caModel.ny();
		if (sup != NULL) sup->resize(nx, ny);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
            coefobs = cd(x, y);
            lambda1 = (caModel(x, y) + cdModel(x, y)) * param / 2;
            lambda2 = (caModel(x, y) - cdModel(x, y)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
            if (p > alpha / 2) 
            {
               cd(x, y) = cdModel(x, y);
               if (sup != NULL) (*sup)(x, y) = (SUPTYPE)-1;
            }
            else
            {
               if (sup != NULL) (*sup)(x, y) = (SUPTYPE)1;
            }
        }        
    }
    else  if (dim == 3)
	{
        nx = caModel.nx(); 
        ny = caModel.ny();
        nz = caModel.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
            coefobs = cd(x, y, z);
            lambda1 = (caModel(x, y, z) + cdModel(x, y, z)) * param / 2;
            lambda2 = (caModel(x, y, z) - cdModel(x, y, z)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
            if (p > alpha / 2) 
            {
               cd(x, y, z) = cdModel(x, y, z);
               if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)-1;
            }
            else
            {
               if (sup != NULL) (*sup)(x, y, z) = (SUPTYPE)1;
            }
        }                
	}
	else throw DataSizeException(dim, "unknown size of dimension");    

    return alpha;
}

template <class DATATYPE, class SUPTYPE>
double WaveletShrinkage<DATATYPE, SUPTYPE>::haarModelFDRThreshold \
       (to_array<DATATYPE, true> &caModel, \
        to_array<DATATYPE, true> &cdModel, to_array<DATATYPE, true> &cd, \
        int scale[], double alpha, bool indep, to_array<SUPTYPE, true> *sup)
{
	int nx, ny, nz;
	double coefobs, lambda1, lambda2, p, fdrp;
	
	int dim = caModel.naxis();
	double ss = 0;
	for (int si=0; si<dim; si++) ss += scale[si];
	double param = POW2(ss);

	if (dim == 1)
	{
		nx = cd.nx();
		if (sup != NULL) sup->resize(nx);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx);
		for (int x=0; x<nx; x++)
		{
            coefobs = cd(x);
            lambda1 = (caModel(x) + cdModel(x)) * param / 2;
            lambda2 = (caModel(x) - cdModel(x)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
			(*pvals)(x) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		{
		  if ((*pvals)(x) > fdrp)
		    {
		      if (sup != NULL)
		      {
			    (*sup)(x) = (SUPTYPE) -1;
		      }
		      cd(x) = cdModel(x);
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			    (*sup)(x) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 2)
	{
	    nx = cd.nx(); ny = cd.ny();
		if (sup != NULL) sup->resize(nx, ny);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
            coefobs = cd(x, y);
            lambda1 = (caModel(x, y) + cdModel(x, y)) * param / 2;
            lambda2 = (caModel(x, y) - cdModel(x, y)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
			(*pvals)(x, y) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		{
		  if ((*pvals)(x, y) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			    (*sup)(x, y) = (SUPTYPE) -1;
		      }
		      cd(x, y) = cdModel(x, y);
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			    (*sup)(x, y) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else if (dim == 3)
	{
  	    nx = cd.nx(); ny = cd.ny(); nz = cd.nz();
		if (sup != NULL) sup->resize(nx, ny, nz);
	    to_array<DATATYPE, true> *pvals = new to_array<DATATYPE, true>(nx, ny, nz);
		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
            coefobs = cd(x, y, z);
            lambda1 = (caModel(x, y, z) + cdModel(x, y, z)) * param / 2;
            lambda2 = (caModel(x, y, z) - cdModel(x, y, z)) * param / 2;
            p = pvalueModelHaar (lambda1, lambda2, param, coefobs);
			(*pvals)(x, y, z) = p;
		}
		fdrp = fdr_pvalue(pvals->buffer(), pvals->n_elem(), alpha/2., indep);

		for (int x=0; x<nx; x++)
		for (int y=0; y<ny; y++)
		for (int z=0; z<nz; z++)
		{
		  if ((*pvals)(x, y, z) > fdrp)
		    { 
		      if (sup != NULL)
		      {
			    (*sup)(x, y, z) = (SUPTYPE) -1;
		      }
		      cd(x, y, z) = cdModel(x, y, z);
		    }
		  else
		    {
		      if (sup != NULL)
		      {
			(*sup)(x, y, z) = (SUPTYPE) 1;
		      }
		    }
		}
		delete pvals; pvals = NULL;
	}
	else throw DataSizeException(dim, "unknown size of dimension");    

	return fdrp * 2.;  // p-value of the distr. of the absolute value of the coeff
}

//------------------------------------------------------------------------------
 
// Given a SubBand1D configured by a wavelet filter, 
// this class do an decimated/undecimated 1D/2D/3D wavelet transform/reconstruction
// when decimated transform/reconstruction is used, the size of the signal along that direction MUST be even
class MWIRWaveletTransform
{
	protected:
		// the 1D/2D/3D transform is realized by 1D transform along different axis
		SubBand1D &sb1Dx, &sb1Dy, &sb1Dz;
		
		// do the transform/reconstruction along a certain axis
		// dist : distance between two filter coeff. ("trous" size + 1)
		// this parameter is not used when decimated transform is required
		// axis : 0,1,2 : X,Y,Z-axis;
		void transformXYZ (fltarray &data, fltarray &ca, fltarray &cd, int dist, int axis, bool dec);
		void reconstructionXYZ (fltarray &ca, fltarray &cd, fltarray &data, int dist, int axis, bool dec);

	public:
		// constructor
		MWIRWaveletTransform (SubBand1D &sb):sb1Dx(sb),sb1Dy(sb),sb1Dz(sb){}
		MWIRWaveletTransform (SubBand1D &sbx, SubBand1D &sby, SubBand1D &sbz):sb1Dx(sbx),sb1Dy(sby),sb1Dz(sbz){}
		
		// transformation/reconstruction 1D (decimated/undecimated)
		// in decimated case, the parameter of scale is not of importance
		void transform1D (fltarray &data1D, fltarray &ch, fltarray &cg, int scale, bool dec);
		void reconstruction1D (fltarray &ch, fltarray &cg, fltarray &data1D, int scale, bool dec);
		
		// transformation/reconstruction 2D (decimated/undecimated)
		// in decimated case, the parameter of scale is not of importance
		// dec[] indicates the decimate conditions along each axis
		void transform2D (fltarray &data2D, fltarray &chh, fltarray &chg, \
					fltarray &cgh, fltarray &cgg, int scale, bool dec[]);
		void reconstruction2D (fltarray &chh, fltarray &chg, fltarray &cgh, \
					fltarray &cgg, fltarray &data2D, int scale, bool dec[]);
		
		// transform/reconstruction decimated/undecimated:
		// given a 3D data volume, do a one step transform and 8 subbands are
		// outputs. scale is the scale expected after transformation or before
		// reconstruction, i.e.
		// scale <= 0,no transform
		// scale = 1, one step decomposition of data
		// scale = 2, one step decomposition of data. The data 
		// 			  is considered as in scale 1, and the "trous" size of 
		//            the filter for undecimated transform will be 2^{2-1}-1=1
		// scale = i, one step decomposition of data. The data is
		//            considered as in scale i-1, and the "trous" size of
		//            the filter for undecimated transform will be 2^{i-1}-1
		// Notice that the parameter "Step" in the function transform of 
		// SubBand1D means "trous" size + 1 = 2^{scale-1}
		// dec[] indicates the decimate conditions along each axis
		void transform3D (fltarray &data3D, fltarray &chhh, fltarray &chhg, \
		           fltarray &chgh, fltarray &chgg, fltarray &cghh, \
		           fltarray &cghg, fltarray &cggh, fltarray &cggg, int scale, bool dec[]);
		void reconstruction3D (fltarray &chhh, fltarray &chhg, \
		           fltarray &chgh, fltarray &chgg, fltarray &cghh, \
		           fltarray &cghg, fltarray &cggh, fltarray &cggg, \
		           fltarray &data3D, int scale, bool dec[]);

		// transform/reconstruction 3D along X and Y
		void transform3DXY (fltarray &data3D, fltarray &chh, fltarray &chg, \
				    fltarray &cgh, fltarray &cgg, int scale, bool dec[]);
		void reconstruction3DXY (fltarray &chh, fltarray &chg, \
					 fltarray &cgh, fltarray &cgg, fltarray &data3D, \
					 int scale, bool dec[]);
		// transform/reconstruction 3D along Z
		void transform3DZ (fltarray &data3D, fltarray &ch, fltarray &cg, \
				    int scale, bool dec);
		void reconstruction3DZ (fltarray &ch, fltarray &cg, \
					 fltarray &data3D, int scale, bool dec);
};

void MWIRWaveletTransform::transformXYZ (fltarray &data, fltarray &ca, fltarray &cd, int dist, int axis, bool dec)
{
	float *line = NULL, *linea = NULL, *lined = NULL;
	int lenx, leny, lenz;
	
	if (axis == 0)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
	        int lenxa = dec ? (lenx+1)/2 : lenx;
		int lenxd = dec ? (lenx+1)/2 : lenx;
		ca.resize(lenxa, leny, lenz);
		cd.resize(lenxd, leny, lenz);
		line = new float[lenx];
		linea = new float[lenxa];
		lined = new float[lenxd];
		
		for (int z=0; z<lenz; z++)
		for (int y=0; y<leny; y++)
		{
			for (int x=0; x<lenx; x++)
				line[x] = data(x, y, z);
			if (!dec)
			  sb1Dx.transform(lenx, line, linea, lined, dist);
			else 
			  sb1Dx.transform(lenx, line, linea, lined);
			if (lenxa == lenxd)
			  { // performance optim
			    for (int x=0; x<lenxa; x++)
			      {
				ca(x, y, z) = linea[x];
				cd(x, y, z) = lined[x];
			      }
			  }
			else
			  {
			    for (int x=0; x<lenxa; x++) ca(x, y, z) = linea[x];
			    for (int x=0; x<lenxd; x++) cd(x, y, z) = lined[x];
			  }
		}
	}	
	else if (axis == 1)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
	        int lenya = dec ? (leny+1)/2 : leny;
		int lenyd = dec ? (leny+1)/2 : leny;
		ca.resize(lenx, lenya, lenz);
		cd.resize(lenx, lenyd, lenz);
		line = new float[leny];
		linea = new float[lenya];
		lined = new float[lenyd];

		for (int z=0; z<lenz; z++)
		for (int x=0; x<lenx; x++)
		{
			for (int y=0; y<leny; y++)
				line[y] = data(x, y, z);
			if (!dec)
			  sb1Dy.transform(leny, line, linea, lined, dist);
			else 
			  sb1Dy.transform(leny, line, linea, lined);
			if (lenya == lenyd)
			  { // performance optim.
			    for (int y=0; y<lenya; y++)
			      {
				ca(x, y, z) = linea[y];
				cd(x, y, z) = lined[y];
			      }
			  }
			else
			  {
			    for (int y=0; y<lenya; y++)	ca(x, y, z) = linea[y];
			    for (int y=0; y<lenyd; y++)	cd(x, y, z) = lined[y];
			  }
		}
	}
	else if (axis == 2)
	{
		lenx = data.nx(); leny = data.ny(); lenz = data.nz();
		int lenza = dec ? (lenz+1) / 2 : lenz;
		int lenzd = dec ? (lenz+1) / 2 : lenz;
		ca.resize(lenx, leny, lenza);
		cd.resize(lenx, leny, lenzd);
		line = new float[lenz];
		linea = new float[lenza];
		lined = new float[lenzd];

		for (int y=0; y<leny; y++)
		for (int x=0; x<lenx; x++)
		{
			for (int z=0; z<lenz; z++)
				line[z] = data(x, y, z);
			if (!dec)
			  sb1Dz.transform(lenz, line, linea, lined, dist);
			else
			  sb1Dz.transform(lenz, line, linea, lined);
			if (lenza == lenzd)
			  { // performance optim.
			    for (int z=0; z<lenza; z++) 
			      {
				ca(x, y, z) = linea[z];
				cd(x, y, z) = lined[z];
			      }
			  }
			else
			  {
			    for (int z=0; z<lenza; z++) ca(x, y, z) = linea[z];
			    for (int z=0; z<lenzd; z++) cd(x, y, z) = lined[z];
			  }
		}
	}
	else throw DataException("MWIRWaveletTransform::transformXYZ: unknown axis");
	
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
}

void MWIRWaveletTransform::reconstructionXYZ (fltarray &ca, fltarray &cd, fltarray &data, int dist, int axis, bool dec)
{
	float *line = NULL, *linea = NULL, *lined = NULL;
	int lenx, leny, lenz;
	
	if (axis == 0)
	{
		leny = ca.ny(); lenz = ca.nz();
		int lenxa = ca.nx(), lenxd = cd.nx();
		lenx = dec ? lenxa+lenxd : lenxa;
		line = new float[lenx];
		linea = new float[lenxa];
		lined = new float[lenxd];
		data.resize(lenx, leny, lenz);
		
		for (int z=0; z<lenz; z++)
		for (int y=0; y<leny; y++)
		{
		  if (lenxa == lenxd)
		    { // performance optim.
		      for (int x=0; x<lenxa; x++)
			{
			  linea[x] = ca(x, y, z);
			  lined[x] = cd(x, y, z);
			}
		    }
		  else
		    {
			for (int x=0; x<lenxa; x++) linea[x] = ca(x, y, z);
			for (int x=0; x<lenxd; x++) lined[x] = cd(x, y, z);
		    }
		  if (!dec)
		    sb1Dx.recons(lenx, linea, lined, line, dist);
		  else
		    sb1Dx.recons(lenx, linea, lined, line);
		  for (int x=0; x<lenx; x++) data(x, y, z) = line[x];
		}
	}
	else if (axis == 1)
	{
		lenx = ca.nx(); lenz = ca.nz();
		int lenya = ca.ny(), lenyd = cd.ny();
		leny = dec ? lenya+lenyd : lenya;
		line = new float[leny];
		linea = new float[lenya];
		lined = new float[lenyd];
		data.resize(lenx, leny, lenz);
		
		for (int z=0; z<lenz; z++)
		for (int x=0; x<lenx; x++)
		{
		  if (lenya == lenyd)
		    { // performance optim.
		      for (int y=0; y<lenya; y++)
			{
			  linea[y] = ca(x, y, z);
			  lined[y] = cd(x, y, z);
			}
		    }
		  else
		    {
			for (int y=0; y<lenya; y++) linea[y] = ca(x, y, z);
			for (int y=0; y<lenyd; y++) lined[y] = cd(x, y, z);
		    }
		  if (!dec)
		    sb1Dy.recons(leny, linea, lined, line, dist);
		  else
		    sb1Dy.recons(leny, linea, lined, line);
		  for (int y=0; y<leny; y++) data(x, y, z) = line[y];
		}
	}
	else if (axis == 2)
	{
		lenx = ca.nx(); leny = ca.ny(); 
		int lenza = ca.nz(), lenzd = cd.nz();
		lenz = dec ? lenza+lenzd : lenza;
		line = new float[lenz];
		linea = new float[lenza];
		lined = new float[lenzd];
		data.resize(lenx, leny, lenz);
		
		for (int y=0; y<leny; y++)
		for (int x=0; x<lenx; x++)
		{
		  if (lenza == lenzd)
		    { // performance optim.
		      for (int z=0; z<lenza; z++)
			{
			  linea[z] = ca(x, y, z);
			  lined[z] = cd(x, y, z);
			}
		    }
		  else
		    {
			for (int z=0; z<lenza; z++) linea[z] = ca(x, y, z);
			for (int z=0; z<lenzd; z++) lined[z] = cd(x, y, z);
		    }
		  if (!dec)
		    sb1Dz.recons(lenz, linea, lined, line, dist);
		  else
		    sb1Dz.recons(lenz, linea, lined, line);
		  for (int z=0; z<lenz; z++) data(x, y, z) = line[z];
		}
	}
	else throw DataException("MWIRWaveletTransform::reconstructionXYZ: unknown axis");
	
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
}

void MWIRWaveletTransform::transform1D (fltarray &data1D, fltarray &ch, fltarray &cg, int scale, bool dec)
{
	if (scale <= 0) return;
	
	int lenx = data1D.nx();
	int lenxa = dec ? (lenx+1)/2 : lenx;
	int lenxd = dec ? (lenx+1)/2 : lenx;

	float *line = new float[lenx];
	float *linea = new float[lenxa];
	float *lined = new float[lenxd];

	for (int x=0; x<lenx; x++) line[x] = data1D(x);
	if (!dec)
	{
	  int dist = POW2(scale-1);
	  sb1Dx.transform(lenx, line, linea, lined, dist);
	}
	else sb1Dx.transform(lenx, line, linea, lined);

	ch.resize(lenxa); cg.resize(lenxd);
	if (lenxa == lenxd)
	  { // performance optim.
	    for (int x=0; x<lenxa; x++)
	      {
		ch(x) = linea[x];
		cg(x) = lined[x];
	      }
	  }
	else
	  {
	    for (int x=0; x<lenxa; x++) ch(x) = linea[x];
	    for (int x=0; x<lenxd; x++) cg(x) = lined[x];
	  }
	
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
}

void MWIRWaveletTransform::reconstruction1D (fltarray &ch, fltarray &cg, fltarray &data1D, int scale, bool dec)
{
	if (scale <= 0) return;
	
	int lenxa = ch.nx(), lenxd = cg.nx();
	int lenx = dec ? lenxa+lenxd : lenxa;

	float *line = new float[lenx];
	float *linea = new float[lenxa];
	float *lined = new float[lenxd];

	data1D.resize(lenx);
	
	if (lenxa == lenxd)
	  {
	    //performance optim.
	    for (int x=0; x<lenxa; x++)
	      {
		linea[x] = ch(x);
		lined[x] = cg(x);
	      }
	  }
	else
	  {
	    for (int x=0; x<lenxa; x++) linea[x] = ch(x);
	    for (int x=0; x<lenxd; x++) lined[x] = cg(x);
	  }
	if (!dec)
	{
	  int dist = POW2(scale-1);
	  sb1Dx.recons(lenx, linea, lined, line, dist);
	}
	else
	{
	  sb1Dx.recons(lenx, linea, lined, line);
	}
	for (int x=0; x<lenx; x++) data1D(x) = line[x];
	
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
}
		
void MWIRWaveletTransform::transform2D (fltarray &data2D, fltarray &chh, fltarray &chg, \
					fltarray &cgh, fltarray &cgg, int scale, bool dec[])
{
	if (scale <= 0) return;

	int lenx = data2D.nx(), leny = data2D.ny();
	int lenxa = dec[0] ? (lenx+1)/2 : lenx;
	int lenxd = dec[0] ? (lenx+1)/2 : lenx;
	int lenya = dec[1] ? (leny+1)/2 : leny;
	int lenyd = dec[1] ? (leny+1)/2 : leny;

	fltarray *temph = new fltarray(lenxa, leny);
	fltarray *tempg = new fltarray(lenxd, leny);
	float *line = new float[lenx];
	float *linea = new float[lenxa];
	float *lined = new float[lenxd];
	
	// X-axis
	int dist = POW2(scale-1);
	for (int y=0; y<leny; y++)
	{
		for (int x=0; x<lenx; x++) line[x] = data2D(x, y);
		if (!dec[0])
		  sb1Dx.transform(lenx, line, linea, lined, dist);
		else
		  sb1Dx.transform(lenx, line, linea, lined);

		if (lenxa == lenxd)
		  {
		    //performance optim.
		    for (int x=0; x<lenxa; x++)
		      {
			(*temph)(x, y) = linea[x];
			(*tempg)(x, y) = lined[x];
		      }
		  }
		else
		  {
		    for (int x=0; x<lenxa; x++) (*temph)(x, y) = linea[x];
		    for (int x=0; x<lenxd; x++) (*tempg)(x, y) = lined[x];
		  }
	}
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }

	line = new float[leny];
	linea = new float[lenya];
	lined = new float[lenyd];
	chh.resize(lenxa, lenya); chg.resize(lenxa, lenyd);
	cgh.resize(lenxd, lenya); cgg.resize(lenxd, lenyd);
	// Y-axis
	for (int x=0; x<lenxa; x++)
	{
		for (int y=0; y<leny; y++) line[y] = (*temph)(x, y);
		if (!dec[1])
		  sb1Dy.transform(leny, line, linea, lined, dist);
		else
		  sb1Dy.transform(leny, line, linea, lined);
		if (lenya == lenyd)
		  {
		    //performance optim.
		    for (int y=0; y<lenya; y++)
		      {
			chh(x, y) = linea[y];
			chg(x, y) = lined[y];
		      }
		  }
		else
		  {
		    for (int y=0; y<lenya; y++) chh(x, y) = linea[y];
		    for (int y=0; y<lenyd; y++) chg(x, y) = lined[y];
		  }
	}
	for (int x=0; x<lenxd; x++)
	{
		for (int y=0; y<leny; y++) line[y] = (*tempg)(x, y);
		if (!dec[1])
		  sb1Dy.transform(leny, line, linea, lined, dist);
		else
		  sb1Dy.transform(leny, line, linea, lined);
		if (lenya == lenyd)
		  {
		    //performance optim.
		    for (int y=0; y<lenya; y++)
		      {
			cgh(x, y) = linea[y];
			cgg(x, y) = lined[y];
		      }
		  }
		else
		  {
		    for (int y=0; y<lenya; y++) cgh(x, y) = linea[y];
		    for (int y=0; y<lenyd; y++) cgg(x, y) = lined[y];
		  }
	}

	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
	delete temph; temph = NULL;
	delete tempg; tempg = NULL;
}

void MWIRWaveletTransform::reconstruction2D (fltarray &chh, fltarray &chg, fltarray &cgh, \
					fltarray &cgg, fltarray &data2D, int scale, bool dec[])
{
	if (scale <= 0) return;
		
	int lenxa = chh.nx(), lenya = chh.ny();
	int lenxd = cgh.nx(), lenyd = chg.ny();
	int lenx = dec[0] ? lenxa+lenxd : lenxa;
	int leny = dec[1] ? lenya+lenyd : lenya;
	fltarray *temph = new fltarray(lenxa, leny);
	fltarray *tempg = new fltarray(lenxd, leny);
	float *line = new float[leny];
	float *linea = new float[lenya];
	float *lined = new float[lenyd];

	// Y-axis
	int dist = POW2(scale-1);
	for (int x=0; x<lenxa; x++)
	{
  	        if (lenya == lenyd)
		  {
		    //performance optim.
		    for (int y=0; y<lenya; y++)
		      {
			linea[y] = chh(x, y);
			lined[y] = chg(x, y);
		      }
		  }
		else
		  {
		    for (int y=0; y<lenya; y++) linea[y] = chh(x, y);
		    for (int y=0; y<lenyd; y++) lined[y] = chg(x, y);
		  }
		if (!dec[1])
		  sb1Dy.recons(leny, linea, lined, line, dist);
		else
		  sb1Dy.recons(leny, linea, lined, line);
		for (int y=0; y<leny; y++) (*temph)(x, y) = line[y];
	}
	for (int x=0; x<lenxd; x++)
	{
  	        if (lenya == lenyd)
		  {
		    //performance optim.
		    for (int y=0; y<lenya; y++)
		      {
			linea[y] = cgh(x, y);
			lined[y] = cgg(x, y);
		      }
		  }
		else
		  {
		    for (int y=0; y<lenya; y++) linea[y] = cgh(x, y);
		    for (int y=0; y<lenyd; y++) lined[y] = cgg(x, y);
		  }
		if (!dec[1])
		  sb1Dy.recons(leny, linea, lined, line, dist);
		else
		  sb1Dy.recons(leny, linea, lined, line);
		for (int y=0; y<leny; y++) (*tempg)(x, y) = line[y];
	}
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }

	line = new float[lenx];
	linea = new float[lenxa];
	lined = new float[lenxd];
	data2D.resize(lenx, leny);
	// X-axis
	for (int y=0; y<leny; y++)
	{
  	        if (lenxa == lenxd)
		  {
		    //performance optim.
		    for (int x=0; x<lenxa; x++)
		      {
			linea[x] = (*temph)(x, y);
			lined[x] = (*tempg)(x, y);
		      }
		  }
		else
		  {
		    for (int x=0; x<lenxa; x++) linea[x] = (*temph)(x, y);
		    for (int x=0; x<lenxd; x++) lined[x] = (*tempg)(x, y);
		  }
		if (!dec[0])
		  sb1Dx.recons(lenx, linea, lined, line, dist);
		else
		  sb1Dx.recons(lenx, linea, lined, line);
		for (int x=0; x<lenx; x++) data2D(x, y) = line[x];
	}
	if (line != NULL) { delete[] line; line = NULL; }
	if (linea != NULL) { delete[] linea; linea = NULL; }
	if (lined != NULL) { delete[] lined; lined = NULL; }
	delete temph; temph = NULL;
	delete tempg; tempg = NULL;
}

void MWIRWaveletTransform::transform3D (fltarray &data3D, fltarray &chhh, fltarray &chhg, \
		           fltarray &chgh, fltarray &chgg, fltarray &cghh, \
		           fltarray &cghg, fltarray &cggh, fltarray &cggg, int scale, bool dec[])
{
	if (scale <= 0) return;
	
	int dist = POW2(scale-1);
	
	// X-axis transform
	fltarray *temph = new fltarray;
	fltarray *tempg = new fltarray;
	transformXYZ(data3D, *temph, *tempg, dist, 0, dec[0]);
	
	// Y-axis transform
	fltarray *temphh = new fltarray;
	fltarray *temphg = new fltarray;
	transformXYZ(*temph, *temphh, *temphg, dist, 1, dec[1]);
	delete temph; temph = NULL;
	
	fltarray *tempgh = new fltarray;
	fltarray *tempgg = new fltarray;
	transformXYZ(*tempg, *tempgh, *tempgg, dist, 1, dec[1]);
	delete tempg; tempg = NULL;
	
	// Z-axis transform
	transformXYZ(*temphh, chhh, chhg, dist, 2, dec[2]);
	delete temphh; temphh = NULL;

	transformXYZ(*temphg, chgh, chgg, dist, 2, dec[2]);
	delete temphg; temphg = NULL;

	transformXYZ(*tempgh, cghh, cghg, dist, 2, dec[2]);
	delete tempgh; tempgh = NULL;

	transformXYZ(*tempgg, cggh, cggg, dist, 2, dec[2]);
	delete tempgg; tempgg = NULL;
}

void MWIRWaveletTransform::reconstruction3D (fltarray &chhh, fltarray &chhg, \
		           fltarray &chgh, fltarray &chgg, fltarray &cghh, \
		           fltarray &cghg, fltarray &cggh, fltarray &cggg, \
		           fltarray &data3D, int scale, bool dec[])
{
	if (scale <= 0) return;
	
	int dist = POW2(scale-1);

	// Z-axis reconstruction	
	fltarray *temphh = new fltarray;
	reconstructionXYZ(chhh, chhg, *temphh, dist, 2, dec[2]);

	fltarray *temphg = new fltarray;
	reconstructionXYZ(chgh, chgg, *temphg, dist, 2, dec[2]);

	fltarray *tempgh = new fltarray;
	reconstructionXYZ(cghh, cghg, *tempgh, dist, 2, dec[2]);

	fltarray *tempgg = new fltarray;
	reconstructionXYZ(cggh, cggg, *tempgg, dist, 2, dec[2]);
	
	// Y-axis reconstruction	
	fltarray *temph = new fltarray;
	reconstructionXYZ(*temphh, *temphg, *temph, dist, 1, dec[1]);
	delete temphh; temphh = NULL; delete temphg; temphg = NULL;
	
	fltarray *tempg = new fltarray;
	reconstructionXYZ(*tempgh, *tempgg, *tempg, dist, 1, dec[1]);
	delete tempgh; tempgh = NULL; delete tempgg; tempgg = NULL;

	// X-axis reconstruction	
	reconstructionXYZ(*temph, *tempg, data3D, dist, 0, dec[0]);
	delete temph; temph = NULL; delete tempg; tempg = NULL;
}

void MWIRWaveletTransform::transform3DXY (fltarray &data3D, fltarray &chh, fltarray &chg, \
				    fltarray &cgh, fltarray &cgg, int scale, bool dec[])
{
  if (scale <= 0) return;
  int dist = POW2(scale-1);
  
  fltarray *ch = new fltarray, *cg = new fltarray;
  
  // X-axis
  transformXYZ (data3D, *ch, *cg, dist, 0, dec[0]);

  // Y-axis
  transformXYZ (*ch, chh, chg, dist, 1, dec[1]);
  delete ch; ch = NULL;
  transformXYZ (*cg, cgh, cgg, dist, 1, dec[1]);
  delete cg; cg = NULL;
}

void MWIRWaveletTransform::reconstruction3DXY (fltarray &chh, fltarray &chg, \
				    fltarray &cgh, fltarray &cgg, fltarray &data3D, int scale, bool dec[])
{
  if (scale <= 0) return;
  int dist = POW2(scale-1);

  fltarray *ch = new fltarray, *cg = new fltarray;
  // Y-axis
  reconstructionXYZ (chh, chg, *ch, dist, 1, dec[1]);
  reconstructionXYZ (cgh, cgg, *cg, dist, 1, dec[1]);

  // X-axis
  reconstructionXYZ (*ch, *cg, data3D, dist, 0, dec[0]);
  delete ch; ch = NULL;
  delete cg; cg = NULL;
}

void MWIRWaveletTransform::transform3DZ (fltarray &data3D, fltarray &ch, fltarray &cg, \
				    int scale, bool dec)
{
  if (scale <= 0) return;
  int dist = POW2(scale-1);

  // Z-axis
  transformXYZ (data3D, ch, cg, dist, 2, dec);
}

void MWIRWaveletTransform::reconstruction3DZ (fltarray &ch, fltarray &cg, \
					 fltarray &data3D, int scale, bool dec)
{
  if (scale <= 0) return;
  int dist = POW2(scale-1);

  // Z-axis
  reconstructionXYZ (ch, cg, data3D, dist, 2, dec);
}

#endif
