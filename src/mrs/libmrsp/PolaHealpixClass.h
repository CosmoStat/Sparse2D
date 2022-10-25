/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Olivier Fourt, with additions by F. Sureau
**
**    Date: 11/05/09 // 05/2017
**    
**    File:  PolaHealpixClass.h
**
*******************************************************************************
**
**    DESCRIPTION  Class interface to the Healpix package   
**    ----------- 
**                 
******************************************************************************/

#ifndef _PolaHealpixClass_H_
#define _PolaHealpixClass_H_

#include "Healpix_PowSpec.h"
#include "HealpixClass.h"
// #include "map_alm_qu_eb.h"
#ifdef UPDATE_HEALPIX3_60_SOLVED
#include "alm2map_pol_der1.h"
#endif

#include "healpix_map.h"
#include "healpix_map_fitsio.h"

#include "compatibility_healpix.h"

/*********************************************************************/

/* A HEALPix Polarized map (3 maps) of a given datatype */

template<typename T> class PolaHmap
{
	private:

		
		bool teb;
		bool nested;
		int Nside;
		
	public:
		Hmap<T> map_T;
		Hmap<T> map_Q;
		Hmap<T> map_U;
		
		/*! Constructs an unallocated map. */
    	PolaHmap () {}
    	    	
    	inline int get_nside() { return Nside; }
    	inline bool flag_teb() { return teb; }
    	inline bool flag_nested() { return nested; }
    	long Npix() {return (long) map_T.Npix(); }

    	void alloc( int NsideIn, bool flag_teb=false, bool flag_nested=false);
		void alloc(arr<T> &MapT, arr<T> &MapQ, arr<T> &MapU, bool flag_teb=false, bool flag_nested=false);
		void read( char * Name, bool flag_teb=false, bool set_ring=false );   // read a Healpix map
		void write( char *Name );  // write a Healpix map
		inline Hmap<T> get_map_T() { return map_T; }
		inline Hmap<T> get_map_Q() { return map_Q; }
		inline Hmap<T> get_map_U() { return map_U; }
		Hmap<T>* get_map(int index) {
			switch (index){
				case 0: return &map_T;break;
				case 1: return &map_Q;break;
				case 2: return &map_U;break;
				default:
					printf("ONLY 3 MAPS AVAILABLE (TQU)\n");
					exit(EXIT_FAILURE); 
			}
		}
		
		REAL & operator()(int PolType, int PixNum) 
		       { if (PolType == 0) { return map_T[PixNum];}
		         else if(PolType == 1) {return map_Q[PixNum];}
			 else if (PolType == 2) return map_U[PixNum];
		       } 

		// inline & Hmap<T> TT() { return map_T; }
		// inline & Hmap<T> Q() { return map_Q; }
		// inline & Hmap<T> U() { return map_U; }
		
		int set_map_T( Hmap<T> & map );
		int set_map_Q( Hmap<T> & map );
		int set_map_U( Hmap<T> & map );

		const PolaHmap<T>  & operator = (const PolaHmap<T> & Pola_map);
     	const PolaHmap<T> & operator += (const PolaHmap<T> & Pola_map);
     	const PolaHmap<T> & operator *= (const PolaHmap<T> & Pola_map);
     	const PolaHmap<T> & operator -= (const PolaHmap<T> & Pola_map);
     	const PolaHmap<T> & operator /= (const PolaHmap<T> & Pola_map);
     	
     	void swap_tqu_teb( bool fast );
     	void swap_nested_ring();
     	
     	// Deletes the old map and creates a new map  with a given \a nside and the ordering scheme \a scheme.
    	void pola_set_nside( int nside, bool flag_teb=false, bool flag_nested=false );
    	//Set to val all pixels
    	void pola_fill(double value );
    	
    	/*! Imports the map \a orig into the current map, adjusting the ordering scheme and the map resolution if necessary.
        When downgrading, \a pessimistic determines whether or not pixels are set to \a Healpix_undef when not all of the corresponding high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
		void import( PolaHmap<T> &pola_map, bool pessimistic=false );
		int import_via_alm( PolaHmap<T> &pola_map_in, bool fast=true );

     	void info();  // print statistical information
};

#define PolaHfmap PolaHmap<float>      	// float healpix map
#define PolaHdmap PolaHmap<double>  	// double healpix map

// ===========================

template<typename T> void PolaHmap<T>::alloc( int NsideIn, bool flag_teb, bool flag_nested)
{
	Nside = NsideIn;
	teb = flag_teb;
	nested = flag_nested;
	
	map_T.alloc( NsideIn, flag_nested );
	map_Q.alloc( NsideIn, flag_nested );
	map_U.alloc( NsideIn, flag_nested );
	
}

// ===========================
template<typename T> void PolaHmap<T>::alloc(arr<T> &MapT, arr<T> &MapQ, arr<T> &MapU, bool flag_teb, bool flag_nested )
{
	Healpix_Ordering_Scheme NFlag;
	if(flag_nested) NFlag=NEST; 
	else NFlag=RING;
	
	long NsideIn=sqrt(MapT.size()/12l);
	if(MapT.size() != (unsigned long) NsideIn*NsideIn*12l) printf("Input Map size not compatible with Healpix Map (%ld vs %ld for Nside = %ld)\n. No Map allocated\n",MapT.size(),NsideIn*NsideIn*12l,NsideIn);
	else if ((MapQ.size() == MapT.size())&& (MapU.size() == MapT.size())) {
		map_T.Set(MapT, NFlag);
		Nside = map_T.Nside();
		teb = flag_teb;
		nested = flag_nested;
		map_Q.Set(MapQ, NFlag);
		map_U.Set(MapU, NFlag);
	}	else printf(" Size of arrays are not compatible: %ld vs %ld vs %ld\n. No Map allocated\n", MapT.size(), MapQ.size(), MapU.size());
}


// ===========================

template<typename T> void PolaHmap<T>::read(char *Name, bool flag_teb, bool set_ring )
{
	char *	NameFits=fitsname(Name);
	 string infile = string(NameFits);
	 free(NameFits);
     // read_Healpix_map_from_fits( infile, map_T, 1, 2 );
     // read_Healpix_map_from_fits( infile, map_Q, 2, 2 );
     // read_Healpix_map_from_fits( infile, map_U, 3, 2 );
     read_Healpix_map_from_fits(infile,map_T,map_Q,map_U);

     Nside = map_T.nside();
     
     teb = flag_teb;
     
     if( map_T.Scheme() != DEF_MRS_ORDERING )
     {
     	nested = true;
     	
     	if( map_Q.Scheme() == DEF_MRS_ORDERING )
     	{
     		map_Q.swap_scheme();
     	}
     	if( map_U.Scheme() == DEF_MRS_ORDERING )
     	{
     		map_U.swap_scheme();
     	}
     }//ordering = NESTED
     else
     {
     	nested = false;
     	
     	if( map_Q.Scheme() != DEF_MRS_ORDERING )
     	{
     		map_Q.swap_scheme();
     	}
     	if( map_U.Scheme() != DEF_MRS_ORDERING )
     	{
     		map_U.swap_scheme();
     	}
     }//ordering = RING
     
     if ( (set_ring == true) && (nested == true) )
     {// Force map to ordering = RING
     	map_T.swap_scheme();
     	map_Q.swap_scheme();
     	map_U.swap_scheme();
     } 
}

// ===========================

template<typename T> void PolaHmap<T>::write(char *Name)
{
    PDT datatype = PLANCK_FLOAT32;
	fitshandle out;
	char *	NameFits=fitsname(Name);
	remove(NameFits);
	out.create(NameFits);
	
	string infile = string(NameFits);
	free(NameFits);
	write_Healpix_map_to_fits(out, map_T, map_Q, map_U, datatype);
}

// ===========================

template<typename T> int PolaHmap<T>::set_map_T( Hmap<T> & map )
{
	int status;
	
	bool flag_map_nested = false;
	
	if ( map.Scheme() != DEF_MRS_ORDERING )
    {
    	flag_map_nested = true;
    }
	
	if( (this->Nside == map.Nside()) && (this->nested == flag_map_nested) )
	{
		map_T = map;
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

template<typename T> int PolaHmap<T>::set_map_Q( Hmap<T> & map )
{
	int status;
	
	bool flag_map_nested = false;
	
	if ( map.Scheme() != DEF_MRS_ORDERING )
    {
    	flag_map_nested = true;
    }
	
	if( (this->Nside == map.Nside()) && (this->nested == flag_map_nested) )
	{
		map_Q = map;
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

template<typename T> int PolaHmap<T>::set_map_U( Hmap<T> & map )
{
	int status;
	
	bool flag_map_nested = false;
	
	if ( map.Scheme() != DEF_MRS_ORDERING )
    {
    	flag_map_nested = true;
    }
	
	if( (this->Nside == map.Nside()) && (this->nested == flag_map_nested) )
	{
		map_U = map;
	
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

template<typename T>  const  PolaHmap<T> & PolaHmap<T>::operator = (const PolaHmap<T> & Pola_map) 
{
	if( this->Nside == Pola_map.get_nside() )
	{
		this->map_T = Pola_map.get_map_T();
		this->map_Q = Pola_map.get_map_Q();
		this->map_U = Pola_map.get_map_U();
		this->teb = Pola_map.flag_teb();
		this->nested = Pola_map.flag_nested();
	}
	
	return (*this);
}

// ===========================

template<typename T>  const  PolaHmap<T> & PolaHmap<T>::operator += (const PolaHmap<T> & Pola_map) 
{
	if( (this->Nside == Pola_map.get_nside()) && (this->teb == Pola_map.flag_teb()) && (this->nested == Pola_map.flag_nested()) )
	{
		this->map_T += Pola_map.get_map_T();
		this->map_Q += Pola_map.get_map_Q();
		this->map_U += Pola_map.get_map_U();
	}
	
	return (*this);
}

// ===========================

template<typename T>  const  PolaHmap<T> & PolaHmap<T>::operator *= (const PolaHmap<T> & Pola_map) 
{
	if( (this->Nside == Pola_map.get_nside()) && (this->teb == Pola_map.flag_teb()) && (this->nested == Pola_map.flag_nested()) )
	{
		this->map_T *= Pola_map.get_map_T();
		this->map_Q *= Pola_map.get_map_Q();
		this->map_U *= Pola_map.get_map_U();
	}
	
	return (*this);
}

// ===========================

template<typename T>  const  PolaHmap<T> & PolaHmap<T>::operator -= (const PolaHmap<T> & Pola_map) 
{
	if( (this->Nside == Pola_map.get_nside()) && (this->teb == Pola_map.flag_teb()) && (this->nested == Pola_map.flag_nested()) )
	{
		this->map_T -= Pola_map.get_map_T();
		this->map_Q -= Pola_map.get_map_Q();
		this->map_U -= Pola_map.get_map_U();
	}
	
	return (*this);
}

// ===========================

template<typename T>  const  PolaHmap<T> & PolaHmap<T>::operator /= (const PolaHmap<T> & Pola_map) 
{
	if( (this->Nside == Pola_map.get_nside()) && (this->teb == Pola_map.flag_teb()) && (this->nested == Pola_map.flag_nested()) )
	{
		this->map_T /= Pola_map.get_map_T();
		this->map_Q /= Pola_map.get_map_Q();
		this->map_U /= Pola_map.get_map_U();
	}
	
	return (*this);
}

// ===========================

template<typename T> void PolaHmap<T>::swap_tqu_teb( bool fast)   // def fast=DEF_ALM_FAST
{
	int Lmax = 3 * Nside;
	int Mmax = Lmax;
	
	arr<double> weight_T;
	weight_T.alloc( 2*Nside );
	if( fast == false )
	{
       	char *HealpixFN = (char *) getenv("HEALPIX");
       	char FN[512];
       	sprintf(FN, "%s/data/", HealpixFN);
		read_weight_ring( FN, Nside, weight_T );
		
		for(int m=0; m < (int) weight_T.size(); ++m)
		{
			weight_T[m]+=1;
		}
	}
	else
	{
		weight_T.fill(1);
	}
	
	if( nested == true )
    {
    	//map_T.swap_scheme();
    	map_Q.swap_scheme();
    	map_U.swap_scheme();
    }//ALM Trans Must be done in RING scheme
	
	if( teb == false )
	{			
       	Alm<xcomplex<T> > ALM_E(Lmax,Mmax), ALM_B(Lmax,Mmax);
       				
		// map2alm_pol_iter_QU( map_Q, map_U, ALM_E, ALM_B, 0, weight_T );
        map2alm_spin(map_Q,map_U,ALM_E,ALM_B,2,weight_T,false);
    	alm2map( ALM_E, map_Q );
    	alm2map( ALM_B, map_U );
		
		teb = true;	
	}
	else
	{
    	Alm<xcomplex<T> > ALM_E(Lmax,Mmax), ALM_B(Lmax,Mmax);
    	 	
    	map2alm_iter( map_Q, ALM_E, 0, weight_T );
    	map2alm_iter( map_U, ALM_B, 0, weight_T );
    	
    	// alm2map_pol_QU( ALM_E, ALM_B, map_Q, map_U );
        alm2map_spin(ALM_E,ALM_B,map_Q,map_U,2);
		teb = false;
	}
	
	if( nested == true )
    {
    	//map_T.swap_scheme();
    	map_Q.swap_scheme();
    	map_U.swap_scheme();
    }//ALM Trans Must be done in RING scheme, move back in NESTED scheme
}

// ===========================

template<typename T> void PolaHmap<T>::swap_nested_ring()
{
	map_T.swap_scheme();
    map_Q.swap_scheme();
    map_U.swap_scheme();
    
    nested = !nested;//invert flag
}

// ===========================
template<typename T> void PolaHmap<T>::pola_fill(double value)
{
   map_T.fill(value);
   map_Q.fill(value);
   map_U.fill(value);
}

// ===========================
template<typename T> void PolaHmap<T>::pola_set_nside( int set_nside, bool set_flag_teb, bool set_flag_nested)
{
	Nside = set_nside;
	teb = set_flag_teb;
	nested = set_flag_nested;
	if( nested == false )
	{
		map_T.SetNside( set_nside, RING );
		map_Q.SetNside( set_nside, RING );
		map_U.SetNside( set_nside, RING );
	}
	else
	{
		map_T.SetNside( set_nside, NEST );
		map_Q.SetNside( set_nside, NEST );
		map_U.SetNside( set_nside, NEST );
	}
}

// ===========================

template<typename T> void PolaHmap<T>::import( PolaHmap<T> &pola_map, bool pessimistic)
{
	Hmap<T> map_T_temp, map_Q_temp, map_U_temp;
	
	map_T_temp.alloc( pola_map.get_nside(), pola_map.flag_nested() );
    map_Q_temp.alloc( pola_map.get_nside(), pola_map.flag_nested() );
   	map_U_temp.alloc( pola_map.get_nside(), pola_map.flag_nested() );
   	
   	map_T_temp = pola_map.get_map_T();
    map_Q_temp = pola_map.get_map_Q();
   	map_U_temp = pola_map.get_map_U();
   	
	if( pola_map.get_nside() == Nside ) // no up/degrading
    {
    	map_T.Import_nograde( map_T_temp );
    	map_Q.Import_nograde( map_Q_temp );
    	map_U.Import_nograde( map_U_temp );
    	
    	nested = pola_map.flag_nested();
    	teb = pola_map.flag_teb();
    }
    else if( pola_map.get_nside() < Nside ) // upgrading
    	 {
    	 	map_T.Import_upgrade( map_T_temp );
    	 	map_Q.Import_upgrade( map_Q_temp );
    	 	map_U.Import_upgrade( map_U_temp );
    	 	
    	 	nested = pola_map.flag_nested();
    		teb = pola_map.flag_teb();
    	 }
    	 else
    	 {
    	 	map_T.Import_degrade( map_T_temp, pessimistic );
    	 	map_Q.Import_degrade( map_Q_temp, pessimistic );
    	 	map_U.Import_degrade( map_U_temp, pessimistic );
    	 	
    	 	nested = pola_map.flag_nested();
    		teb = pola_map.flag_teb();
    	 }
}

// ===========================

template<typename T> int PolaHmap<T>::import_via_alm( PolaHmap<T> &pola_map_in, bool fast )
{
	int status;
	
	if( Nside >= pola_map_in.get_nside() )
	{
		Hmap<T> map_T_temp, map_Q_temp, map_U_temp;
		
		bool pola_map_in_nested = false;
		if( pola_map_in.flag_nested() == true )
		{
			pola_map_in_nested = true;
			//pola_map_in.swap_nested_ring();
		}
	
		map_T_temp.alloc( pola_map_in.get_nside(), pola_map_in.flag_nested() );
    	map_Q_temp.alloc( pola_map_in.get_nside(), pola_map_in.flag_nested() );
   		map_U_temp.alloc( pola_map_in.get_nside(), pola_map_in.flag_nested() );
   	
   		map_T_temp = pola_map_in.get_map_T();
    	map_Q_temp = pola_map_in.get_map_Q();
   		map_U_temp = pola_map_in.get_map_U();
   		
   		if( pola_map_in.flag_teb() == true )
   		{//Map in TEB
   			map_T.import_via_alm( map_T_temp, fast );
    		map_Q.import_via_alm( map_Q_temp, fast );
    		map_U.import_via_alm( map_U_temp, fast );
   		}
   		else
   		{//Map in TQU
   			if( pola_map_in_nested == true )
   			{
   				map_T_temp.swap_scheme();
   				map_Q_temp.swap_scheme();
   				map_U_temp.swap_scheme();
   			}
   				
			int Lmax_in = 3 * pola_map_in.get_nside();
			int Mmax_in = Lmax_in;
	
			arr<double> weight_T;
			weight_T.alloc( 2*pola_map_in.get_nside() );
			if( fast == false )
			{
        		char *HealpixFN = (char *) getenv("HEALPIX");
        		char FN[512];
       		 	sprintf(FN, "%s/data/", HealpixFN);
				read_weight_ring( FN, Nside, weight_T );
		
				for(int m=0; m < (int) weight_T.size(); ++m)
				{
					weight_T[m]+=1;
				}
			}
			else
			{
				weight_T.fill(1);
			}
			
			Alm<xcomplex<T> > ALM_T_in( Lmax_in, Mmax_in );
			Alm<xcomplex<T> > ALM_E_in( Lmax_in, Mmax_in );
			Alm<xcomplex<T> > ALM_B_in( Lmax_in, Mmax_in );
			
			double avg = map_T_temp.average();
			map_T_temp.Add(-avg);
    
    		map2alm_pol_iter( map_T_temp, map_Q_temp, map_U_temp, ALM_T_in, ALM_E_in, ALM_B_in, 0, weight_T );
		
			ALM_T_in(0,0) += avg*sqrt(fourpi);
			
			int Lmax_this = 3*Nside;
			int Mmax_this = Lmax_this;
			Alm<xcomplex<T> > ALM_T_this( Lmax_this, Mmax_this );
			ALM_T_this.SetToZero();
			Alm<xcomplex<T> > ALM_E_this( Lmax_this, Mmax_this );
			ALM_E_this.SetToZero();
			Alm<xcomplex<T> > ALM_B_this( Lmax_this, Mmax_this );
			ALM_B_this.SetToZero();
			
			for( int l=0; l <= Lmax_in; l++ )
			{
				for( int m=0; m <= l; m++ )
				{
					ALM_T_this( l, m ) = ALM_T_in( l, m );
					ALM_E_this( l, m ) = ALM_E_in( l, m );
					ALM_B_this( l, m ) = ALM_B_in( l, m );
				}
			}
			
			double offset = ALM_T_this(0,0).real()/sqrt(fourpi);
			ALM_T_this(0,0) = 0;
			
			alm2map_pol( ALM_T_this, ALM_E_this, ALM_B_this, map_T, map_Q, map_U );
		
   			map_T.Add(offset);
   			
   			if( pola_map_in_nested == true )
			{
				this->swap_nested_ring();
			}		
   		}
   		
    	nested = pola_map_in_nested;
    	teb = pola_map_in.flag_teb();
    	
    	status = 1;
	}
	else
	{
		(*this).import( pola_map_in );
		
		status = -1;
	}
	
	return status;
}

// ===========================

template<typename T> void PolaHmap<T>::info()
{
	if( teb == false )
	{
		if( nested == false )
		{
			cout << "Map TQU, RING Scheme, Nside = " << Nside << endl;
		}
		else
		{
			cout << "Map TQU, NESTED Scheme, Nside = " << Nside << endl;
		}
		
		cout << "Info map T:" << endl;
		map_T.info();
		cout << "Info map Q:" << endl;
		map_Q.info();
		cout << "Info map U:" << endl;
		map_U.info();
	}
	else
	{
		if( nested == false )
		{
			cout << "Map TEB, RING Scheme, Nside = " << Nside << endl;
		}
		else
		{
			cout << "Map TEB, NESTED Scheme, Nside = " << Nside << endl;
		}
		
		cout << "Info map T:" << endl;
		map_T.info();
		cout << "Info map E:" << endl;
		map_Q.info();
		cout << "Info map B:" << endl;
		map_U.info();
	}
}

// ===========================

class PolaAlmR
{
	private:
		CAlmR ALM_T;		//AlmR ALM_T; //FCS MOD
		CAlmR ALM_E;		//AlmR ALM_E; //FCS MOD
		CAlmR ALM_B;		//AlmR ALM_B; //FCS MOD
		
		int PolaNside;
		int Pola_Lmax;
		bool PolaFastALM;
		int niter;
		
		bool Norm;  // if True, normalize the Alm coefficients such a Gaussian noise noise with variance S^2, produce Alm with variance 1/2 S^2
		double NormVal; // Alm Nomalization value
				
		arr<double> weight_T; // internal variable
		
		fltarray WienerFilter_T;
		fltarray WienerFilter_E;
		fltarray WienerFilter_B;
		
	public:
		PolaAlmR(){}
		
		inline int get_nside(){ return PolaNside; }
		inline int get_lmax(){ return Pola_Lmax; }
		inline bool flag_FastALM(){ return PolaFastALM; }
		inline int get_niter(){ return niter; }
		inline double get_normval(){ return NormVal; }
		inline bool flag_NormALM(){ return Norm; }
		
		inline CAlmR* get_alm(int index) {
			switch (index){
				case 0: return &ALM_T;break;
				case 1: return &ALM_E;break;
				case 2: return &ALM_B;break;
				default:
					printf("ONLY 3 ALM COMPONENTS AVAILABLE (TEB)\n");
					exit(EXIT_FAILURE); 
			}
		}

		
		fltarray get_WienerFilter_T(){ return WienerFilter_T; }
		fltarray get_WienerFilter_E(){ return WienerFilter_E; }
		fltarray get_WienerFilter_B(){ return WienerFilter_B; }
		
		void set_niter( int nb_iter ){ niter = nb_iter; }
		void set_flag_NormALM( bool flag_norm );
		
		int set_Alm_T( AlmR & ALM );
		int set_Alm_E( AlmR & ALM );
		int set_Alm_B( AlmR & ALM );
		
		void alloc(int Nside, int u_lmax, bool Fast);
	
		void read( char *Name, int Nside, int u_lmax=0, bool Fast=DEF_ALM_FAST  ); // read the Alm from a fits file
		void read_array( char *Name );
		void write( char *Name );  // write the Alm to fits file	
		void write_array( char *Name );	// write Alm "array" format in files "ALM_T_Name.xxx" "ALM_E_Name.xxx" "ALM_B_Name.xxx"
		void write_array2( char *Name );
		
		void pola_alm_trans( PolaHmap<REAL> & PolaMap ); // Alm transformation of a Healpix map
		void pola_alm_rec( PolaHmap<REAL> & PolaMap, bool KeepAlm, int RecNside ); // Alm inverse transform
		
		void pola_first_deriv( PolaHmap<REAL> & PolaMap, PolaHmap<REAL> & PolaMap_d_theta, PolaHmap<REAL> & PolaMap_d_phi, float fwhm_arcmin, int RecNside ); // Alm inverse transform, computation of 1st order derivatives
		
		void pola_alm2powspec( Healpix_PowSpec & spec );  // compute the power spectrum from the Alm
		void pola_alm2powspec_all( Healpix_PowSpec & spec );
		
		void convol( fltarray &Filter_T, fltarray &Filter_E, fltarray &Filter_B ); // convol the Alm's with a given filter for each component
		void convol( fltarray &Filter ); // convol the Alm's with the same filter for all components
		void convol( float Fwhm );
		
		void wiener( Healpix_PowSpec & ps_noise, Healpix_PowSpec & ps_signal ); // Apply a wiener filtering to the Alm, knowing the signal and noise power spectrum
		
		double max_absalm_T();
		double max_absalm_E();
		double max_absalm_B();
		double max_absalm_TEB();

		int hard_threshold( float lambda_t, float lambda_e, float lambda_b, int & MaxNonZeroL_T, int & MaxNonZeroL_E, int & MaxNonZeroL_B );//threshold coeff lower thant T and return the number of non zero coeff after thresholding
		int soft_threshold( float lambda_t, float lambda_e, float lambda_b, int & MaxNonZeroL_T, int & MaxNonZeroL_E, int & MaxNonZeroL_B );
		
		void info();
};

// ===========================

void PolaAlmR::alloc( int Nside, int u_lmax, bool Fast )
{
	PolaNside = Nside;
	niter = ALM_DEF_NITER;
	
	Pola_Lmax=u_lmax;
    if ( Pola_Lmax <= 0 ) Pola_Lmax = 3*Nside;
    if ( Pola_Lmax > ALM_MAX_L_TOT ) Pola_Lmax = ALM_MAX_L_TOT;
    //int Max = L_Max;
    PolaFastALM = Fast;
    
    weight_T.alloc( 2*PolaNside );
    if( PolaFastALM == false )
	{
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[512];
        sprintf(FN, "%s/data/", HealpixFN);
		read_weight_ring( FN, PolaNside, weight_T );
		for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
 	}
	else weight_T.fill(1);
     
    double lm = sqrt( (double) (PolaNside)* (double) (PolaNside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    Norm = false;
        
    ALM_T.Set( Pola_Lmax, Pola_Lmax );
    ALM_E.Set( Pola_Lmax, Pola_Lmax );
    ALM_B.Set( Pola_Lmax, Pola_Lmax );
}

// ===========================

void PolaAlmR::set_flag_NormALM( bool flag_norm )
{
	Norm = flag_norm;
	/*
	double lm = sqrt( (double) (PolaNside)* (double) (Nside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    */
}

// ===========================

int PolaAlmR::set_Alm_T( AlmR & ALM )
{
	int status;
	
	if( (Pola_Lmax == ALM.Lmax()) && (Pola_Lmax == ALM.Mmax()) )
	{		
		for(int l=0; l <= ALM.Lmax(); l++)
    	{
    		for(int m=0; m <= l; m++)
    		{
	    		ALM_T(l,m) = ALM(l,m);
    		}
    	}
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;	
}

// ===========================

int PolaAlmR::set_Alm_E( AlmR & ALM )
{
	int status;
	
	if( (Pola_Lmax == ALM.Lmax()) && (Pola_Lmax == ALM.Mmax()) )
	{		
		for(int l=0; l <= ALM.Lmax(); l++)
    	{
    		for(int m=0; m <= l; m++)
    		{
	    		ALM_E(l,m) = ALM(l,m);
    		}
    	}
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;	
}

// ===========================

int PolaAlmR::set_Alm_B( AlmR & ALM )
{
	int status;
	
	if( (Pola_Lmax == ALM.Lmax()) && (Pola_Lmax == ALM.Mmax()) )
	{		
		for(int l=0; l <= ALM.Lmax(); l++)
    	{
    		for(int m=0; m <= l; m++)
    		{
	    		ALM_B(l,m) = ALM(l,m);
    		}
    	}
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;	
}

// ===========================

void PolaAlmR::read( char *Name, int Nside, int u_lmax, bool Fast)
{
	PolaNside = Nside;
	niter = ALM_DEF_NITER;
	
	Pola_Lmax=u_lmax;
    if ( Pola_Lmax <= 0 ) Pola_Lmax = 3*Nside;
    if ( Pola_Lmax > ALM_MAX_L_TOT ) Pola_Lmax = ALM_MAX_L_TOT;

    PolaFastALM = Fast;
    
    weight_T.alloc( 2*PolaNside );
    if( PolaFastALM == false )
	{
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[512];
        sprintf(FN, "%s/data/", HealpixFN);
        read_weight_ring( FN, PolaNside, weight_T );
		for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
	}
	else weight_T.fill(1);
	
	double lm = sqrt( (double) (PolaNside)* (double) (PolaNside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    Norm = false;
	
	ALM_T.Set( Pola_Lmax, Pola_Lmax );
    ALM_E.Set( Pola_Lmax, Pola_Lmax );
    ALM_B.Set( Pola_Lmax, Pola_Lmax );
    
	char *	NameFits=fitsname(Name);
	string infile = string(NameFits);
	free(NameFits);
    read_Alm_from_fits( infile, ALM_T, Pola_Lmax, Pola_Lmax, 2 );
    read_Alm_from_fits( infile, ALM_E, Pola_Lmax, Pola_Lmax, 3 );
    read_Alm_from_fits( infile, ALM_B, Pola_Lmax, Pola_Lmax, 4 );
}

// ===========================

void PolaAlmR::write( char *Name )
{
    // PDT datatype = PLANCK_FLOAT32;
	char *	NameFits=fitsname(Name);
	string infile = string(NameFits);
	remove(NameFits);
	free(NameFits);
	fitshandle out;
	out.create(infile);
	 
    write_Alm_to_fits( out, ALM_T, ALM_T.Lmax(), ALM_T.Mmax(), PLANCK_FLOAT64);
    write_Alm_to_fits( out, ALM_E, ALM_E.Lmax(), ALM_E.Mmax(), PLANCK_FLOAT64);
    write_Alm_to_fits( out, ALM_B, ALM_B.Lmax(), ALM_B.Mmax(), PLANCK_FLOAT64);
    // write_Healpix_map_to_fits (out,ALM_T,ALM_E,ALM_B,  datatype  );
}

// ===========================

void PolaAlmR::write_array( char *Name)
{
	char file_T_Name[512];
	char file_E_Name[512];
	char file_B_Name[512];
	
	strcpy( file_T_Name, "ALM_T_" );
	strcpy( file_E_Name, "ALM_E_" );
	strcpy( file_B_Name, "ALM_B_" );
	
	strcat( file_T_Name, Name );
	strcat( file_E_Name, Name );
	strcat( file_B_Name, Name );
	
	dblarray A_T;
    A_T.alloc( ALM_T.Lmax()+1, ALM_T.Lmax()+1, 2 );
    
    dblarray A_E;
    A_E.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1, 2 );
    
    dblarray A_B;
    A_B.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1, 2 );
    
    for(int l=0; l <= ALM_T.Lmax(); l++)
    {
    	for(int m=0; m <= l; m++)
    	{
    		A_T(l,m,0) = ALM_T(l,m).real();
    		A_T(l,m,1) = ALM_T(l,m).imag();
    			
    		A_E(l,m,0) = ALM_E(l,m).real();
    		A_E(l,m,1) = ALM_E(l,m).imag();
    			
    		A_B(l,m,0) = ALM_B(l,m).real();
    		A_B(l,m,1) = ALM_B(l,m).imag();
    	}
    }
	
	fits_write_dblarr( file_T_Name, A_T );
	fits_write_dblarr( file_E_Name, A_E );
	fits_write_dblarr( file_B_Name, A_B );
}

// ===========================

void PolaAlmR::write_array2( char *Name)
{
	dblarray ALM;
    ALM.alloc( 1 + 3*( Pola_Lmax+1 ), Pola_Lmax+1, 2 );
        
    for(int l=0; l <= Pola_Lmax; l++)
    {
    	for(int m=0; m <= l; m++)
    	{
    		ALM( l, m, 0 ) = ALM_T(l,m).real();
    		ALM( l, m, 1 ) = ALM_T(l,m).imag();
    			
    		ALM( l+Pola_Lmax+1, m, 0 ) = ALM_E(l,m).real();
    		ALM( l+Pola_Lmax+1, m, 1 ) = ALM_E(l,m).imag();
    			
    		ALM( l+2*Pola_Lmax+2, m, 0 ) = ALM_B(l,m).real();
    		ALM( l+2*Pola_Lmax+2, m, 1 ) = ALM_B(l,m).imag();
    	}
    }
    
    //Write option used for almtrans
    ALM( 3*(Pola_Lmax+1), 0, 0 ) = PolaNside;
    
    if( PolaFastALM == true )
    {
    	ALM( 3*(Pola_Lmax+1), 1, 0 ) = 1;
    }
    else
    {
    	ALM( 3*(Pola_Lmax+1), 1, 0 ) = 0;
    }
    
    if( Norm == true )
    {
    	ALM( 3*(Pola_Lmax+1), 2, 0 ) = 1;
    }
    else
    {
    	ALM( 3*(Pola_Lmax+1), 2, 0 ) = 0;
    }

	fits_write_dblarr( Name, ALM );	
}

// ===========================

void PolaAlmR::read_array( char *Name )
{
	dblarray ALM;
	char *	NameFits=fitsname(Name);
	fits_read_dblarr(NameFits, ALM );
	free(NameFits);
	
	Pola_Lmax = ALM.ny()-1;

	PolaNside = (int)ALM( 3*(Pola_Lmax+1), 0, 0 );
		
	niter = ALM_DEF_NITER;
	
	double file_fast_var = ALM( 3*(Pola_Lmax+1), 1, 0 );
	if( file_fast_var == 0.0 )
	{
		PolaFastALM = false;
	}
	else
	{
		PolaFastALM = true;
	}
    
    weight_T.alloc( 2*PolaNside );
    if( PolaFastALM == false )
	{
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[512];
        sprintf(FN, "%s/data/", HealpixFN);
        read_weight_ring( FN, PolaNside, weight_T );
		for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
	}
	else weight_T.fill(1);
 	
	double lm = sqrt( (double) (PolaNside)* (double) (PolaNside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    
    double file_norm_var = ALM( 3*(Pola_Lmax+1), 2, 0 );
	if( file_norm_var == 0.0 )
	{
		Norm = false;
	}
	else
	{
		Norm = true;
	}
	
	ALM_T.Set( Pola_Lmax, Pola_Lmax );
    ALM_E.Set( Pola_Lmax, Pola_Lmax );
    ALM_B.Set( Pola_Lmax, Pola_Lmax );
    
    for(int l=0; l <= Pola_Lmax; l++)
    {
    	for(int m=0; m <= l; m++)
    	{
    		ALM_T(l,m) = xcomplex<REAL> (ALM( l, m, 0 ), ALM( l, m, 1 ));
    		ALM_E(l,m) = xcomplex<REAL> ( ALM( l+Pola_Lmax+1, m, 0 ), ALM( l+Pola_Lmax+1, m, 1 ));
    			
    		ALM_B(l,m) = xcomplex<REAL> (ALM( l+2*Pola_Lmax+2, m, 0 ), ALM( l+2*Pola_Lmax+2, m, 1 ));
    	}
    }   
}

// ===========================

void PolaAlmR::pola_alm_trans( PolaHmap<REAL> & PolaMap )
{
	Hmap<REAL> map_T_temp, map_Q_temp, map_U_temp;
   	map_T_temp = PolaMap.get_map_T();
    map_Q_temp = PolaMap.get_map_Q();
   	map_U_temp = PolaMap.get_map_U();
   	
   	if( PolaMap.flag_nested() == true )
	{//Map is in NEST scheme
	    //std::cout<<"SWAPPING NESTED->RING"<<std::endl;
		map_T_temp.swap_scheme();
    	map_Q_temp.swap_scheme();
    	map_U_temp.swap_scheme();
	}

   	
   	if( PolaMap.flag_teb() == false )
	{
	   //std::cout<<"TQU TO TEB THRU ALMS"<<std::endl;
		//Map is in TQU		
		double avg = map_T_temp.average();
		map_T_temp.Add(-avg);
    	map2alm_pol_iter( map_T_temp, map_Q_temp, map_U_temp, ALM_T, ALM_E, ALM_B, niter, weight_T );
		
		ALM_T(0,0) += avg*sqrt(fourpi);
	}
	else
	{
	    //std::cout<<"TQU TO TEB THRU ALMS"<<std::endl;
		//Map is in TEB		
		double avg = map_T_temp.average();
		map_T_temp.Add(-avg);

    	map2alm_iter( map_T_temp, ALM_T, niter, weight_T);
	   	ALM_T(0,0) += avg*sqrt(fourpi);
	   	
	   	map2alm_iter( map_Q_temp, ALM_E, niter, weight_T);
	   	map2alm_iter( map_U_temp, ALM_B, niter, weight_T);
	}
	
	if( Norm == true )
	{
	    for(int l=0; l <= Pola_Lmax; l++)
	    {
	    	for(int m=0; m <= l; m++)
	    	{
	    		ALM_T(l,m) *= NormVal;
	    		ALM_E(l,m) *= NormVal;
	    		ALM_B(l,m) *= NormVal;
	    	}
	    }
	}
}

// ===========================

void PolaAlmR::pola_alm_rec( PolaHmap<REAL> & PolaMap, bool KeepAlm=false, int RecNside=0 )
{
	Hmap<REAL> map_T_temp, map_Q_temp, map_U_temp;
	
	if( RecNside != 0 )
	{
		map_T_temp.alloc( RecNside );
    	map_Q_temp.alloc( RecNside );
   		map_U_temp.alloc( RecNside );
   		
   		PolaNside = RecNside;
   		
   		PolaMap.pola_set_nside( RecNside, PolaMap.flag_teb(), PolaMap.flag_nested() );
	}
	else
	{
		map_T_temp.alloc( PolaNside );
    	map_Q_temp.alloc( PolaNside );
   		map_U_temp.alloc( PolaNside );
	}
	
	dblarray A_T, A_E, A_B;
    if( KeepAlm == true )
    {
    	A_T.alloc( Pola_Lmax+1, Pola_Lmax+1, 2 );
    	A_E.alloc( Pola_Lmax+1, Pola_Lmax+1, 2 );
    	A_B.alloc( Pola_Lmax+1, Pola_Lmax+1, 2 );
    	
    	for(int l=0; l <= Pola_Lmax; l++)
    	{
    		for(int m=0; m <= l; m++)
    		{
    			A_T(l,m,0) = ALM_T(l,m).real();
    			A_T(l,m,1) = ALM_T(l,m).imag();
    			
    			A_E(l,m,0) = ALM_E(l,m).real();
    			A_E(l,m,1) = ALM_E(l,m).imag();
    			
    			A_B(l,m,0) = ALM_B(l,m).real();
    			A_B(l,m,1) = ALM_B(l,m).imag();
    		}
    	}
    }
     
    if( Norm == true )
	{
	    for(int l=0; l <= Pola_Lmax; l++)
	    {
	    	for(int m=0; m <= l; m++)
	    	{
	    		ALM_T(l,m) *= (1./NormVal);
	    		ALM_E(l,m) *= (1./NormVal);
	    		ALM_B(l,m) *= (1./NormVal);
	    	}
	    }
	}
	
	if( PolaMap.flag_teb() == false )
	{
		//Map is in TQU
		double offset = ALM_T(0,0).real()/sqrt(fourpi);
		ALM_T(0,0) = 0.0;
		
		alm2map_pol( ALM_T, ALM_E, ALM_B, map_T_temp, map_Q_temp, map_U_temp );
		
   		map_T_temp.Add(offset);
	}
	else
	{
		//Map is in TEB
		double offset = ALM_T(0,0).real()/sqrt(fourpi);
		ALM_T(0,0) = 0.0;
		 
    	alm2map( ALM_T, map_T_temp );
    	map_T_temp.Add(offset);
    	
    	alm2map( ALM_E, map_Q_temp );
    	alm2map( ALM_B, map_U_temp );
	}
	
	if( PolaMap.flag_nested() == true )
	{
		map_T_temp.swap_scheme();
    	map_Q_temp.swap_scheme();
    	map_U_temp.swap_scheme();
	}
		
	PolaMap.set_map_T( map_T_temp );
	PolaMap.set_map_Q( map_Q_temp );
	PolaMap.set_map_U( map_U_temp );
	
	if( KeepAlm == true )
    {
     	for(int l=0; l <= Pola_Lmax; l++)
    	{
    		for(int m=0; m <= l; m++)
    		{
    			ALM_T(l,m)  = xcomplex<REAL> ( A_T(l,m,0),  A_T(l,m,1));
    			ALM_E(l,m)  = xcomplex<REAL> ( A_E(l,m,0),A_E(l,m,1));
    			ALM_B(l,m)  = xcomplex<REAL> ( A_B(l,m,0), A_B(l,m,1));
    		}
    	}
    }
}

// ===========================

void PolaAlmR::pola_alm2powspec( Healpix_PowSpec & spec )
{
	extract_Healpix_PowSpec_pola( ALM_T, ALM_E, ALM_B, spec );//Calculate TT, EE, BB and TE spectrums!
}

// ===========================

void PolaAlmR::pola_alm2powspec_all( Healpix_PowSpec & spec )
{	
	extract_Healpix_PowSpec_pola_all( ALM_T, ALM_E, ALM_B, spec );//Calculate TT, EE, BB, TE, TB and EB spectrums!
}

// ===========================

void PolaAlmR::info()
{
	fltarray A_T_re, A_T_im, A_E_re, A_E_im, A_B_re, A_B_im;
    A_T_re.alloc( ALM_T.Lmax()+1, ALM_T.Lmax()+1 );
    A_T_im.alloc( ALM_T.Lmax()+1, ALM_T.Lmax()+1 );
    A_E_re.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1 );
    A_E_im.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1 );
    A_B_re.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1 );
    A_B_im.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1 );
       
    for(int l=0; l <= Pola_Lmax; l++)
    {
    	for(int m=0; m <= l; m++)
    	{
    		A_T_re(l,m) = ALM_T(l,m).real();
    		A_T_im(l,m) = ALM_T(l,m).imag();
    			
    		A_E_re(l,m) = ALM_E(l,m).real();
    		A_E_im(l,m) = ALM_E(l,m).imag();
    			
    		A_B_re(l,m) = ALM_B(l,m).real();
    		A_B_im(l,m) = ALM_B(l,m).imag();
    	}
    }
    
    cout << "Alm T: " << " Re ==> Min = " << A_T_re.min() << ", Max = " << A_T_re.max() << ", Sig = " << A_T_re.sigma() << endl;
    cout << "Alm T: " << " Re ==> Min = " << A_T_re.min() << ", Max = " << A_T_re.max() << ", Sig = " << A_T_re.sigma() << endl;
    cout << "Alm E: " << " Re ==> Min = " << A_E_re.min() << ", Max = " << A_E_re.max() << ", Sig = " << A_E_re.sigma() << endl;
    cout << "Alm E: " << " Re ==> Min = " << A_E_re.min() << ", Max = " << A_E_re.max() << ", Sig = " << A_E_re.sigma() << endl;
    cout << "Alm B: " << " Re ==> Min = " << A_B_re.min() << ", Max = " << A_B_re.max() << ", Sig = " << A_B_re.sigma() << endl;
    cout << "Alm B: " << " Re ==> Min = " << A_B_re.min() << ", Max = " << A_B_re.max() << ", Sig = " << A_B_re.sigma() << endl;
        
}

// ===========================

void PolaAlmR::convol( fltarray &Filter_T, fltarray &Filter_E, fltarray &Filter_B )
{
	int l,m;
	
    int LMin_T = Filter_T.nx()-1;
    int LMin_E = Filter_E.nx()-1;
    int LMin_B = Filter_B.nx()-1;
    
    if( LMin_T > Pola_Lmax )
    {
    	LMin_T = Pola_Lmax;
    }
    if( LMin_E > Pola_Lmax )
    {
    	LMin_E = Pola_Lmax;
    }
    if( LMin_B > Pola_Lmax )
    {
    	LMin_B = Pola_Lmax;
    }
    	
    for( l=0; l <= LMin_T; l++)
    {
    	for( m=0; m <= l; m++)
    	{
    		ALM_T(l,m) *= Filter_T(l);
    	}
    }
	for( l=LMin_T+1; l <= Pola_Lmax; l++)
	{
		for( m=0; m <= l; m++)
		{
			ALM_T(l,m) *= 0.;
		}
	}
	
	for( l=0; l <= LMin_E; l++)
    {
    	for( m=0; m <= l; m++)
    	{
    		ALM_E(l,m) *= Filter_E(l);
    	}
    }
	for( l=LMin_E+1; l <= Pola_Lmax; l++)
	{
		for( m=0; m <= l; m++)
		{
			ALM_E(l,m) *= 0.;
		}
	}
	
	for( l=0; l <= LMin_B; l++)
    {
    	for( m=0; m <= l; m++)
    	{
    		ALM_B(l,m) *= Filter_B(l);
    	}
    }
	for( l=LMin_B+1; l <= Pola_Lmax; l++)
	{
		for( m=0; m <= l; m++)
		{
			ALM_B(l,m) *= 0.;
		}
	}
}

// ===========================

void PolaAlmR::convol( fltarray &Filter )
{
	int l,m;
	
    int LMin = Filter.nx()-1;
    
    if( LMin > Pola_Lmax )
    {
    	LMin = Pola_Lmax;
    }
    	
    for( l=0; l <= LMin; l++)
    {
    	for( m=0; m <= l; m++)
    	{
    		ALM_T(l,m) *= Filter(l);
    		ALM_E(l,m) *= Filter(l);
    		ALM_B(l,m) *= Filter(l);
    	}
    }
	for( l=LMin+1; l <= Pola_Lmax; l++)
	{
		for( m=0; m <= l; m++)
		{
			ALM_T(l,m) *= 0.;
			ALM_E(l,m) *= 0.;
			ALM_B(l,m) *= 0.;
		}
	}
}

// ===========================

void PolaAlmR::convol( float Fwhm )
{
   smoothWithGauss( ALM_T, ALM_E, ALM_B, (double) Fwhm );
}

// ===========================

void PolaAlmR::wiener( Healpix_PowSpec & ps_noise, Healpix_PowSpec & ps_signal )
{
	WienerFilter_T.resize( Pola_Lmax+1 );
    WienerFilter_T(0) = 1.;
    WienerFilter_T(1) = 1.;
    
    WienerFilter_E.resize( Pola_Lmax+1 );
    WienerFilter_E(0) = 1.;
    WienerFilter_E(1) = 1.;
    
    WienerFilter_B.resize( Pola_Lmax+1 );
    WienerFilter_B(0) = 1.;
    WienerFilter_B(1) = 1.;
	
    for( int l=2; l <= Pola_Lmax; l++ )
    {
  	    double Num = ps_signal.tt(l);
	    double Den = ps_signal.tt(l) + ps_noise.tt(l);
	    double WienFilter = (Den <= 0) ? 0: Num / Den;
	    WienerFilter_T(l) = WienFilter;
	    
	    Num = ps_signal.ee(l);
	    Den = ps_signal.ee(l) + ps_noise.ee(l);
	    WienFilter = (Den <= 0) ? 0: Num / Den;
	    WienerFilter_E(l) = WienFilter;
	    
	    Num = ps_signal.bb(l);
	    Den = ps_signal.bb(l) + ps_noise.bb(l);
	    WienFilter = (Den <= 0) ? 0: Num / Den;
	    WienerFilter_B(l) = WienFilter;
	    
     	for( int m=0; m <= l; m++ )
     	{
     		ALM_T(l,m) *= WienerFilter_T(l);
     		ALM_E(l,m) *= WienerFilter_E(l);
     		ALM_B(l,m) *= WienerFilter_B(l);
     	}
	}
}

// ===========================

double PolaAlmR::max_absalm_TEB()
{
	double Max =0.;
	Max=MAX(max_absalm_T(),max_absalm_E());
	Max=MAX(Max,max_absalm_B());
	return Max;
}

// ===========================

double PolaAlmR::max_absalm_T()
{
	double Max=0.;
	
    for( int l=2; l <= Pola_Lmax; l++ )
    {
    	for( int m=0; m <= l; m++ ) 
		{
			if( ABS( ALM_T(l,m).real()) > Max )
			{
				Max = ABS(ALM_T(l,m).real());
			}
			if( ABS( ALM_T(l,m).imag()) > Max )
			{
				Max = ABS(ALM_T(l,m).imag());
			}
		}
    }
	return Max;
}

// ===========================

double PolaAlmR::max_absalm_E()
{
	double Max=0.;
	
    for( int l=2; l <= Pola_Lmax; l++ )
    {
    	for( int m=0; m <= l; m++ ) 
		{
			if( ABS( ALM_E(l,m).real()) > Max )
			{
				Max = ABS(ALM_E(l,m).real());
			}
			if( ABS( ALM_E(l,m).imag()) > Max )
			{
				Max = ABS(ALM_E(l,m).imag());
			}
		}
    }
	return Max;
}

// ===========================

double PolaAlmR::max_absalm_B()
{
	double Max=0.;
	
    for( int l=2; l <= Pola_Lmax; l++ )
    {
    	for( int m=0; m <= l; m++ ) 
		{
			if( ABS( ALM_B(l,m).real()) > Max )
			{
				Max = ABS(ALM_B(l,m).real());
			}
			if( ABS( ALM_B(l,m).imag()) > Max )
			{
				Max = ABS(ALM_B(l,m).imag());
			}
		}
    }
	return Max;
}

// ===========================

int PolaAlmR::hard_threshold( float lambda_t, float lambda_e, float lambda_b, int & MaxNonZeroL_T, int & MaxNonZeroL_E, int & MaxNonZeroL_B )
{
	int Cpt=0;
	MaxNonZeroL_T = 0;
	MaxNonZeroL_E = 0;
	MaxNonZeroL_B = 0;
	
    for( int l=1; l <= Pola_Lmax; l++ )
    {
    	for( int m=0; m <= l; m++ )
    	{
			if( ABS(ALM_T(l,m).real()) < lambda_t )
			{
				ALM_T(l,m) = xcomplex<REAL> ( 0.,ALM_T(l,m).imag());
			}
			else 
			{
				Cpt++;
				if( MaxNonZeroL_T < l )
				{
					MaxNonZeroL_T = l;
				}
			}
			if( ABS(ALM_T(l,m).imag()) < lambda_t )
			{
				ALM_T(l,m) = xcomplex<REAL> (ALM_T(l,m).real(),0.);
			}
			else    	 
			{
				Cpt++;
				if( MaxNonZeroL_T < l )
				{
					MaxNonZeroL_T = l;
				}
			}
			
			if( ABS(ALM_E(l,m).real()) < lambda_e )
			{
				ALM_E(l,m) = xcomplex<REAL> ( 0.,ALM_E(l,m).imag());
			}
			else 
			{
				Cpt++;
				if( MaxNonZeroL_E < l )
				{
					MaxNonZeroL_E = l;
				}
			}
			if( ABS(ALM_E(l,m).imag()) < lambda_e )
			{
				ALM_E(l,m) = xcomplex<REAL> (ALM_E(l,m).real(),0.);
			}
			else    	 
			{
				Cpt++;
				if( MaxNonZeroL_E < l )
				{
					MaxNonZeroL_E = l;
				}
			}
			
			if( ABS(ALM_B(l,m).real()) < lambda_b )
			{
				ALM_B(l,m) = xcomplex<REAL> ( 0.,ALM_B(l,m).imag());
			}
			else 
			{
				Cpt++;
				if( MaxNonZeroL_B < l )
				{
					MaxNonZeroL_B = l;
				}
			}
			if( ABS(ALM_B(l,m).imag()) < lambda_b )
			{
				ALM_B(l,m) = xcomplex<REAL> (ALM_B(l,m).real(),0.);
			}
			else    	 
			{
				Cpt++;
				if( MaxNonZeroL_B < l )
				{
					MaxNonZeroL_B = l;
				}
			}
		}
    }
	return Cpt;
}

// ===========================

int PolaAlmR::soft_threshold( float lambda_t, float lambda_e, float lambda_b, int & MaxNonZeroL_T, int & MaxNonZeroL_E, int & MaxNonZeroL_B )
{
	int Cpt=0;
	MaxNonZeroL_T = 0;
	MaxNonZeroL_E = 0;
	MaxNonZeroL_B = 0;
	
	double norm_T = 0.0;
	double norm_E = 0.0;
	double norm_B = 0.0;
	
	double Coef_t, Coef_e, Coef_b;
	
    for( int l=1; l <= Pola_Lmax; l++ )
    {
    	for( int m=0; m <= l; m++ )
    	{
			norm_T = sqrt( norm( ALM_T(l,m) ) );
			norm_E = sqrt( norm( ALM_E(l,m) ) );
			norm_B = sqrt( norm( ALM_B(l,m) ) );

			if( norm_T == 0 )
			{
                ALM_T(l,m) = xcomplex<REAL> (0.,0.);
				Cpt++;
			}
			else
			{
				Coef_t = ( 1. - lambda_t / norm_T );
				if( Coef_t <= 0 ) 
				{
					Coef_t = 0.;
					ALM_T(l,m) = xcomplex<REAL> (0.,0.);
					Cpt++;
				}
				else
				{
					if( MaxNonZeroL_T < l )
					{
						MaxNonZeroL_T = l;
					}
					ALM_T(l,m) *= Coef_t;   
				} 	    	
			}
			
			if( norm_E == 0 )
			{
				ALM_E(l,m) = xcomplex<REAL> (0.,0.);
				Cpt++;
			}
			else
			{
				Coef_e = ( 1. - lambda_e / norm_E );
				if( Coef_e <= 0 ) 
				{
					Coef_e = 0.;
					ALM_E(l,m) = xcomplex<REAL> (0.,0.);
					Cpt++;
				}
				else
				{
					if( MaxNonZeroL_E < l )
					{
						MaxNonZeroL_E = l;
					}
					ALM_E(l,m) *= Coef_e;   
				} 	    	
			}
			
			if( norm_B == 0 )
			{
				ALM_B(l,m) = xcomplex<REAL> (0.,0.);
				Cpt++;
			}
			else
			{
				Coef_b = ( 1. - lambda_b / norm_B );
				if( Coef_b <= 0 ) 
				{
					Coef_b = 0.;
					ALM_B(l,m) = xcomplex<REAL> (0.,0.);
					Cpt++;
				}
				else
				{
					if( MaxNonZeroL_B < l )
					{
						MaxNonZeroL_B = l;
					}
					ALM_B(l,m) *= Coef_b;   
				} 	    	
			}
		}
    }
	return Cpt;
}

// ===========================


void PolaAlmR::pola_first_deriv( PolaHmap<REAL> & PolaMap, PolaHmap<REAL> & PolaMap_d_theta, PolaHmap<REAL> & PolaMap_d_phi, float fwhm_arcmin=0.0, int RecNside=0 )
{
	Hmap<REAL> map_T_temp, map_Q_temp, map_U_temp;
	Hmap<REAL> map_T_d_theta_temp, map_Q_d_theta_temp, map_U_d_theta_temp;
	Hmap<REAL> map_T_d_phi_temp, map_Q_d_phi_temp, map_U_d_phi_temp;
	
	if( RecNside != 0 )
	{
		map_T_temp.alloc( RecNside );
    	map_Q_temp.alloc( RecNside );
   		map_U_temp.alloc( RecNside );
   		
   		PolaNside = RecNside;
   		
   		PolaMap.pola_set_nside( RecNside, PolaMap.flag_teb(), PolaMap.flag_nested() );
   		
   		map_T_d_theta_temp.alloc( RecNside );
    	map_Q_d_theta_temp.alloc( RecNside );
   		map_U_d_theta_temp.alloc( RecNside );
   		
   		PolaMap_d_theta.pola_set_nside( RecNside, PolaMap.flag_teb(), PolaMap.flag_nested() );
   		
   		map_T_d_phi_temp.alloc( RecNside );
    	map_Q_d_phi_temp.alloc( RecNside );
   		map_U_d_phi_temp.alloc( RecNside );
   		
   		PolaMap_d_phi.pola_set_nside( RecNside, PolaMap.flag_teb(), PolaMap.flag_nested() );
	}
	else
	{
		map_T_temp.alloc( PolaNside );
    	map_Q_temp.alloc( PolaNside );
   		map_U_temp.alloc( PolaNside );
   		
   		map_T_d_theta_temp.alloc( PolaNside );
    	map_Q_d_theta_temp.alloc( PolaNside );
   		map_U_d_theta_temp.alloc( PolaNside );
   		
   		map_T_d_phi_temp.alloc( PolaNside );
    	map_Q_d_phi_temp.alloc( PolaNside );
   		map_U_d_phi_temp.alloc( PolaNside );
	}
	
	if( Norm == true )
	{
	    for(int l=0; l <= Pola_Lmax; l++)
	    {
	    	for(int m=0; m <= l; m++)
	    	{
	    		ALM_T(l,m) *= (1./NormVal);
	    		ALM_E(l,m) *= (1./NormVal);
	    		ALM_B(l,m) *= (1./NormVal);
	    	}
	    }
	}
	
	if( fwhm_arcmin > 0 )
	{
		smoothWithGauss( ALM_T, ALM_E, ALM_B, fwhm_arcmin );
	}
	
	double offset = ALM_T(0,0).real()/sqrt(fourpi);
	ALM_T(0,0) = 0.0;
	
	if( PolaMap.flag_teb() == true )
	{
		// Map is in TEB
		alm2map_der1( ALM_T, map_T_temp, map_T_d_theta_temp, map_T_d_phi_temp );
		alm2map_der1( ALM_E, map_Q_temp, map_Q_d_theta_temp, map_Q_d_phi_temp );
		alm2map_der1( ALM_B, map_U_temp, map_U_d_theta_temp, map_U_d_phi_temp );
	}
	else
	{
        // JLS May 2, 2020, the following routine cannot compile anymore with Healpix 3.60.
        // As it seems not used in other codes, it can be removed
#ifdef UPDATE_HEALPIX3_60_SOLVED
		// Map is in TQU
		alm2map_pol_der1( ALM_T, ALM_E, ALM_B, map_T_temp, map_Q_temp, map_U_temp, map_T_d_theta_temp, map_Q_d_theta_temp, map_U_d_theta_temp, map_T_d_phi_temp, map_Q_d_phi_temp, map_U_d_phi_temp );
#else
        cout << "Error in POlaHealpixCLASS: data must be TUB and TQU." << endl;
        exit(-1);
#endif
	}
	
	map_T_temp.Add(offset);
	
	if( PolaMap.flag_nested() == true )
	{
		map_T_temp.swap_scheme();
    	map_Q_temp.swap_scheme();
    	map_U_temp.swap_scheme();
    	
    	map_T_d_theta_temp.swap_scheme();
    	map_Q_d_theta_temp.swap_scheme();
    	map_U_d_theta_temp.swap_scheme();
    	
    	map_T_d_phi_temp.swap_scheme();
    	map_Q_d_phi_temp.swap_scheme();
    	map_U_d_phi_temp.swap_scheme();
	}
		
	PolaMap.set_map_T( map_T_temp );
	PolaMap.set_map_Q( map_Q_temp );
	PolaMap.set_map_U( map_U_temp );
	
	PolaMap_d_theta.set_map_T( map_T_d_theta_temp );
	PolaMap_d_theta.set_map_Q( map_Q_d_theta_temp );
	PolaMap_d_theta.set_map_U( map_U_d_theta_temp );
	
	PolaMap_d_phi.set_map_T( map_T_d_phi_temp );
	PolaMap_d_phi.set_map_Q( map_Q_d_phi_temp );
	PolaMap_d_phi.set_map_U( map_U_d_phi_temp );
}

// ===========================

void mrs_write_pola_spec( char *Name, PowSpec & PS )
{
	fltarray pola_spec( PS.Lmax()+1, PS.Num_specs() );
	
	for(int l=0; l <= PS.Lmax(); l++)
	{
		pola_spec( l, 0 ) = PS.tt(l);
		pola_spec( l, 1 ) = PS.gg(l);
		pola_spec( l, 2 ) = PS.cc(l);
		pola_spec( l, 3 ) = PS.tg(l);
	}
	char *	NameFits=fitsname(Name);
    fits_write_fltarr(NameFits, pola_spec );//File read using MRS/astron/pro/fits_read.pro	FITS_READ, 'filename.fits', spec
	free(NameFits);
}

// ===========================

#endif
