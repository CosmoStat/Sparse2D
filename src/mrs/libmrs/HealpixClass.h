/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date: 19/09/08
**    
**    File:  HealpixClass.h
**
*******************************************************************************
**
**    DESCRIPTION  Class interface to the Healpix package   
**    ----------- 
**                 
******************************************************************************/

#ifndef _HealpixClass_H_
#define _HealpixClass_H_

#include <string> 
#include "Array.h"
#include "IM_IO.h"
#include "DefMath.h"
// #include "OptMedian.h"
// #include "NR.h"

// externals
#include "healpix_map.h"
#include "fitshandle.h"
#include "healpix_data_io.h"
#include "healpix_map_fitsio.h"
#include "xcomplex.h"
#include "lsconstants.h"
#include "planck_rng.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"

#define  DEF_MRS_ORDERING RING
#define  DEF_MRS_CL_LCDM "def_cl.fits"  // extracted from "wmap_comb_tt_powspec_3yr_v2.fits"

#define SIGMA2FWHM 2.3548201
#define FWHM2SIGMA (1./SIGMA2FWHM)   // Fwhm = 2 sqrt(2 log2) sigma

#define ALM_MAX_L 4200

#define REAL double

/*********************************************************************/

template<typename T> class Hmap: public Healpix_Map<T>
// healpix map class
{
   public:
      inline int nside() {  return Healpix_Map<T>::Nside(); }
      T * buffer ();
      void alloc( int NsideIn, bool nested=false );
    
      void read(char * Name, bool set_ring=true);   // read a Healpix map
      void write(char *Name);  // write a Healpix map 
      void get_one_face( dblarray &Face, int face_num); // extract one face of the healpix map
      void get_one_face( fltarray &Face, int face_num);
      void get_all_face( dblarray &AFace);  // Extract all faces: output is a cube [*,*,12]
      void get_all_face( fltarray &AFace);
      void put_one_face( dblarray &Face, int face_num); // Insert one face of the Healpix map
      void put_one_face( fltarray &Face, int face_num);
      void put_all_face( dblarray &AFace); // Insert all faces
      void put_all_face( fltarray &AFace);
      void info();  // print statistical information
      void info(char *Mes);  // print statistical information
      int xyf2ind (int ix, int iy, int f) 
       { return (    ((*this).Scheme()==NEST) ?   Healpix_Map<T>::xyf2nest(ix, iy, f): Healpix_Map<T>::xyf2ring(ix, iy, f));   }
      void ind2xyf (int ind, int & ix, int & iy, int & f) 
       { return ( ((*this).Scheme()==NEST) ? Healpix_Map<T>::nest2xyf((int64) ind, ix, iy, f): Healpix_Map<T>::ring2xyf2((int64) ind, ix, iy, f));   }

      T &operator() (int ix, int iy, int f) { return (    ((*this).Scheme()==NEST) ?  (*this) [Healpix_Map<T>::xyf2nest(ix, iy, f)]:   (*this) [Healpix_Map<T>::xyf2ring(ix, iy, f)]);   }
     const Healpix_Map<T> & operator = (const Healpix_Map<T> & Mat);
     const Healpix_Map<T> & operator += (const Healpix_Map<T> & Mat);
     const Healpix_Map<T>  & operator *= (const Healpix_Map<T> & Mat);
     const Healpix_Map<T>  & operator -= (const Healpix_Map<T> & Mat);
     const Healpix_Map<T> & operator /= (const Healpix_Map<T> & Mat);
     double var();
     double sigma();
     void init() { (*this).fill(0.);}
     int N;
     
     int import_via_alm( Healpix_Map<T> &map_in, bool fast=true );
};


#define Hfmap Hmap<float>      // float healpix map
#define Hdmap Hmap<double>  // double healpix map

// class functions

// ===========================

template<typename T> void Hmap<T>::alloc( int NsideIn, bool nested)
{
    // std::cout <<std::cout << "ALLOC IN " << std::endl;
	if( nested == false )
	{
        // std::cout << "RING: " <<NsideIn << std::endl;
		(*this).SetNside( NsideIn, DEF_MRS_ORDERING );
	}
	else
	{
		(*this).SetNside( NsideIn, NEST );
	}
    // std::cout << "END ALLOC" << std::endl;
}

// ===========================

template<typename T> T*  Hmap<T>::buffer()
{
  arr<T> A;
  T *Ret = & ((*this)[0]);
  // A = (*this).Map();
  return (Ret);	
}

// ===========================

template<typename T>  double  Hmap<T>::var() 
{
  double s=0.;
  double a = this->average();
  for (int x=0; x< this->Npix(); x++)  s+= pow(((*this)[x]-a), (double) 2.);
  return (s / (double) this->Npix() );
}
// ===========================

template<typename T>  double  Hmap<T>::sigma() 
{
  // return ( sqrt( this->var()) );
  return (  this->rms() );
}

// ===========================

template<typename T>  const  Healpix_Map<T> &  Hmap<T>::operator = (const Healpix_Map<T> & Mat) 
{
  for (int x=0; x< this->Npix(); x++) (*this)[x] = Mat[x];
  return (*this);
}


// ===========================

template<typename T>  const  Healpix_Map<T> &  Hmap<T>::operator += (const Healpix_Map<T> & Mat) 
{
  for (int x=0; x< this->Npix(); x++) (*this)[x] += Mat[x];
  return (*this);
}

// ===========================

template<typename T> const Healpix_Map<T> &  Hmap<T>::operator -= (const Healpix_Map<T> & Mat) 
{
  for (int x=0; x<this->Npix(); x++) (*this)[x] -= Mat[x];
  return (*this);
}

// ===========================

template<typename T>  const  Healpix_Map<T> &  Hmap<T>::operator *= (const Healpix_Map<T> & Mat) 
{
  for (int x=0; x< this->Npix(); x++) (*this)[x] *= Mat[x];
  return (*this);
}

// ===========================

template<typename T> const Healpix_Map<T> &  Hmap<T>::operator /= (const Healpix_Map<T> & Mat) 
{
  for (int x=0; x<this->Npix(); x++) (*this)[x] /= Mat[x];
  return (*this);
}

// ===========================


template<typename T> void Hmap<T>::read(char *Name, bool set_ring)
{
	 char *NameFits=fitsname(Name);
	 std::string infile = std::string(NameFits);
	 free(NameFits);
     read_Healpix_map_from_fits(infile,(*this) ,1,2);
     
     if( set_ring == true )
     {
     	if( this->Scheme() !=  DEF_MRS_ORDERING )
     	{
     		this->swap_scheme();
     	}
     }
}

// ===========================

template<typename T> void Hmap<T>::write(char *Name)
{
   fitshandle out;
   char *NameFits=fitsname(Name);
   remove(NameFits);
   out.create(NameFits);
   std::string infile = std::string(NameFits);
   free(NameFits);
   write_Healpix_map_to_fits (out, (*this), planckType<T>());
}
// ===========================

template<typename T> void Hmap<T>::info()
{
   T Min,Max;
   (*this).minmax(Min,Max);
   if ((*this).Scheme()==NEST) std::cout << "  NESTED: " ;
   else std::cout << "  RING: ";
   std::cout << ",   Npix = " << Healpix_Map<T>::Npix() << " , nside = " <<  Healpix_Map<T>::Nside() << std::endl;
   std::cout << "  Min = " << Min << " Max = " << Max <<  " Mean = " <<  Healpix_Map<T>::average() <<  "  sigma = " <<  sigma() << std::endl;
}

// ===========================

template<typename T> void Hmap<T>::info(char *Mes)
{
   T Min,Max;
   (*this).minmax(Min,Max);
   if ((*this).Scheme()==NEST) std::cout << Mes << ":  NESTED" << " nside = " <<  Healpix_Map<T>::Nside();
   else std::cout << Mes << ":  RING" << Healpix_Map<T>::Nside();
   std::cout << ",   Npix = " << Healpix_Map<T>::Npix() << std::endl;
   std::cout << "      Min = " << Min << ", Max = " << Max;
   std::cout << ", Mean = " << Healpix_Map<T>::average() << ", sigma = " <<  sigma() << std::endl;
 }

// ===========================

template<typename T> void Hmap<T>::get_one_face( dblarray &Face, int face_num)
{
	int Nside = nside();
	Face.resize(Nside, Nside);
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   Face(i,j) = (*this) (i, j, face_num);
	}
	
    // int xyf2nest(int ix, int iy, int face_num) const;
    // void nest2xyf(int pix, int &ix, int &iy, int &face_num) const;
}

// ===========================

template<typename T> void Hmap<T>::put_one_face( dblarray &Face, int face_num)
{
	int Nside = nside();
 	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   (*this) (i, j, face_num) = Face(i,j);
	}
}

// ===========================

template<typename T> void Hmap<T>::get_one_face( fltarray &Face, int face_num)
{
	int Nside = nside();
	Face.resize(Nside, Nside);
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   Face(i,j) = (*this) (i, j, face_num);
	}
	
    // int xyf2nest(int ix, int iy, int face_num) const;
    // void nest2xyf(int pix, int &ix, int &iy, int &face_num) const;
}

// ===========================

template<typename T> void Hmap<T>::put_one_face( fltarray &Face, int face_num)
{
	int Nside = nside();
 	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   (*this) (i, j, face_num) = Face(i,j);
	}
}

// ===========================

template<typename T> void Hmap<T>::get_all_face( dblarray &AFace)
{
	int Nside = nside();
	AFace.resize(Nside, Nside,12);
	for (int f=0; f < 12; f++)
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   AFace(i,j,f) = (*this) (i, j, f);
	}
}

// ===========================

template<typename T> void Hmap<T>::put_all_face( dblarray &AFace)
{
	int Nside = nside();
 	for (int f=0; f < 12; f++)
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   (*this) (i, j, f) = AFace(i,j,f);
	}
}

// ===========================

template<typename T> void Hmap<T>::get_all_face( fltarray &AFace)
{
	int Nside = nside();
	AFace.resize(Nside, Nside,12);
	for (int f=0; f < 12; f++)
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   AFace(i,j,f) = (*this) (i, j, f);
	}
}

// ===========================

template<typename T> void Hmap<T>::put_all_face( fltarray &AFace)
{
	int Nside = nside();
 	for (int f=0; f < 12; f++)
	for (int i=0; i < Nside; i++)
	for (int j=0; j < Nside; j++) 
	{
  	   (*this) (i, j, f) = AFace(i,j,f);
	}
}

// ===========================

template<typename T> int Hmap<T>::import_via_alm( Healpix_Map<T> &map_in, bool fast)
{
	int status;
	
	if( this->nside() >= map_in.Nside() )
	{
		bool map_in_nested = false;
		if( map_in.Scheme() == NEST )
		{
			map_in_nested = true;
			map_in.swap_scheme();
		}
		if( this->Scheme() !=  DEF_MRS_ORDERING )
     	{
     		this->swap_scheme();
     	}
		
		int Lmax_in = 3 * map_in.Nside();
		int Mmax_in = Lmax_in;
	
		arr<double> weight_T;
		weight_T.alloc( 2*map_in.Nside() );
		if( fast == false )
		{
			read_weight_ring( "/Users/Shared/software/Healpix_2.10/data/", map_in.Nside(), weight_T );
		
			for(int m=0; m <(int)  weight_T.size(); ++m)
			{
				weight_T[m]+=1;
			}
		}
		else
		{
			weight_T.fill(1);
		}
		
		Alm<xcomplex<T> > ALM_in( Lmax_in, Mmax_in );
		
		double avg=map_in.average();
    		map_in.Add(-avg);
    	
		map2alm_iter( map_in, ALM_in, 0, weight_T );
		
		map_in.Add(avg);
		
		ALM_in(0,0) += avg*sqrt(fourpi);
		
		int Lmax_this = 3*(*this).nside();
		int Mmax_this = Lmax_this;
		Alm<xcomplex<T> > ALM_this( Lmax_this, Mmax_this );
		ALM_this.SetToZero();
		
		for( int l=0; l <= Lmax_in; l++ )
		{
			for( int m=0; m <= l; m++ )
			{
				ALM_this( l, m ) = ALM_in( l, m );
			}
		}
		
		double offset = ALM_this(0,0).real()/sqrt(fourpi);
		ALM_this(0,0) = 0;
		
		alm2map( ALM_this, (*this) );
		
		(*this).Add(offset);
		
		if( map_in_nested == true )
		{
			this->swap_scheme();
			map_in.Add(-avg);
			map_in.swap_scheme();
		}
	
		status = 1;
	}
	else
	{
		(*this).Import( map_in );
		
		status = -1;
	}
	
	return status;
}

int  mrs_get_lmax (int  & Lmax, int Nside, float ZeroPadding=0.);

/*********************************************************************/

//  void mrs_write_3maps(char *Name, Hmap<REAL> & map, Hmap<REAL> & mapdth, Hmap< > & mapdph);
 
#define MAX_ZERO_PADDING 0.
#define ALM_MAX_L_TOT  (int)(ALM_MAX_L+ALM_MAX_L*MAX_ZERO_PADDING)
   
#define ALM_DEF_NITER 0
#define DEF_ALM_FAST false
//  
//   
#define AlmR Alm<xcomplex<double> >    
class CAlmR: public AlmR
{
	int AllocNside;  // Set to nside when the class is allocated 
	bool FastALM; // true if we use the mode fast Alm computation (not as accurate as the standard mode, but faster)
    // get_ring_weights (params,par,map.Nside(),weight_T);
  	arr<double> weight_T; // internal variable

	public:
		int Niter;  // Number of iterations in the ALM reconstructions
		bool Verbose; // Verbose mode
		bool Norm;  // if true, normalize the Alm coefficients such a Gaussian noise noise with variance S^2, produce Alm with variance 1/2 S^2
 	    double NormVal; // Alm Nomalization value
 	    bool UseBeamEff; // if true, multiply after transformation and before reconstruction with the effective beam BeamEff
		inline double normval() { return NormVal;}
	
		CAlmR () : AlmR() { Norm = false; UseBeamEff = false; AllocNside=0; Niter=ALM_DEF_NITER; FastALM=DEF_ALM_FAST;Verbose=false;}
 	    void alloc (int Nside, int u_lmax=0, bool Fast=DEF_ALM_FAST); 
 	    void alm_trans(Hmap<double> & Map); // Alm transformation of a Healpix map
		void alm_rec(Hmap<double> & Map, bool KeepAlm=false, int RecNside=0);   // Alm inverse transform
		void read(char *Name, int Nside, bool Fast=false); // read the Alm from a fits file
		void write(char *Name, bool Array=false);  // write the Alm to fits file
		void info(char *Name);
		void alm2powspec(PowSpec & powspec);  // compute the power spectrum from the Alm
		void wiener(PowSpec & ps_noise, PowSpec & ps_signal); // Apply a wiener filtering to the Alm, knowing the signal and noise power spectrum
		void wiener(float SigmaNoise);  // Apply a wiener filtering to the Alm, assuming white Gaussian noise
		void set_beam_eff(int Lmax0, int Lmax1);  // the the effective beam: B[0:Lamx1] = 1, B[lmax0:*] = 0 and B decrease from Lmax1 to Lmax1
		void convol(fltarray &Filter); // convol the data with a given filter
		void convol(float Fwhm);
		fltarray WienerFilter;  // Wiener filtered computed in the two wiener routines
		fltarray BeamEff;       // Beam effective array computed in set_beam_eff
 	    void  deriv(Hmap<double> & Map, Hmap<double> & mapdth,  Hmap<double> & mapdph, float fwhm_arcmin=0, int RecNside=0);
 	    double max_absalm();
 	    int hard_threshold(float T, int & MaxNonZeroL); // threshold coeff lower thant T and return the number of non zero coeff after thresholding
 	    int soft_threshold(float T, int & MaxNonZeroL);
 	    xcomplex<double> lm_soft_threshold(int l, int m, float T);
 	    void extract_median_powspec(PowSpec &powspec);
};

// =============================

void mrs_alloc_powspec(PowSpec & PS, int Lmax);
void mrs_alloc_powspec(PowSpec & PS, fltarray &Cl);
void mrs_alloc_powspec(PowSpec & PS, dblarray &Cl);

void mrs_write_powspec(char *Name, PowSpec & PS);

void mrs_read_powspec(char *Name, PowSpec & PS);

void mrs_write_powspec(char *Name, PowSpec & powspec);


#endif
