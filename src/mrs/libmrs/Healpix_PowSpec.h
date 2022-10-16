/***********************************************************************************

	Extension of HealPIX class PowSpec
	
	see files powpsec.h, alm_powspec_tools.h, powspec_fitsio.h and *.cc files
	
	enable work with the six pola spectrums
	
************************************************************************************/


#ifndef _HP_POWSPEC_
#define _HP_POWSPEC_

#include "IM_IO.h"
#include <cmath>
#include "arr.h"
#include "alm.h"
#include "planck_rng.h"
#include "xcomplex.h"
#include <string>
#include "HealpixClass.h"

class Healpix_PowSpec
{
	private:
    
    	void dealloc();
    
	public:
	arr<double> tt_, ee_, bb_, te_, tb_, eb_;
    	int num_specs;

		Healpix_PowSpec(){}
		
		/*! Constructs a \a Healpix_PowSpec with \a nspec components and a maximum
        multipole of \a lmax. \a nspec can be 1 (TT), 4 (TT,GG,CC,TG) or
        6 (TT,GG,CC,TG,TC,GC). */
	    Healpix_PowSpec( int nspec, int lmax );
	    
	    /*! Returns the number of spectral components. */
    	int Num_specs() const { return num_specs; }
    	/*! Returns the maximum \a l. */
    	int Lmax() const { return tt_.size()-1; }
    	/*! Returns the TT array (read-only). */
    	const arr<double> &tt() const { return tt_; }
    	/*! Returns the GG array (read-only). */
    	const arr<double> &ee() const { return ee_; }
   		/*! Returns the CC array (read-only). */
    	const arr<double> &bb() const { return bb_; }
    	/*! Returns the TG array (read-only). */
    	const arr<double> &te() const { return te_; }
    	/*! Returns the TC array (read-only). */
    	const arr<double> &tb() const { return tb_; }
    	/*! Returns the GC array (read-only). */
    	const arr<double> &eb() const { return eb_; }
    	
    	/*! Returns TT(l) (read-write). */
    	double &tt (int l) { return tt_[l]; }
    	/*! Returns GG(l) (read-write). */
    	double &ee (int l) { return ee_[l]; }
    	/*! Returns CC(l) (read-write). */
    	double &bb (int l) { return bb_[l]; }
    	/*! Returns TG(l) (read-write). */
    	double &te (int l) { return te_[l]; }
    	/*! Returns TC(l) (read-write). */
    	double &tb (int l) { return tb_[l]; }
    	/*! Returns GC(l) (read-write). */
    	double &eb (int l) { return eb_[l]; }
    	
    	/*! Returns TT(l) (read-only). */
    	const double &tt (int l) const { return tt_[l]; }
    	/*! Returns GG(l) (read-only). */
    	const double &ee (int l) const { return ee_[l]; }
    	/*! Returns CC(l) (read-only). */
    	const double &bb (int l) const { return bb_[l]; }
    	/*! Returns TG(l) (read-only). */
    	const double &te (int l) const { return te_[l]; }
    	/*! Returns TC(l) (read-only). */
    	const double &tb (int l) const { return tb_[l]; }
    	/*! Returns GC(l) (read-only). */
    	const double &eb (int l) const { return eb_[l]; }
      	
      	
      	// Set without deallocation	!!
		int Set_TT( arr<double> &tt_new );
		int Set_EE( arr<double> &ee_new );
		int Set_BB( arr<double> &bb_new );
		int Set_TE( arr<double> &te_new );
		int Set_TB( arr<double> &tb_new );
		int Set_EB( arr<double> &eb_new );
		
		/*! Sets the whole TT array. */
    	void Set( arr<double> &tt_new );
    	/*! Sets all components (4). */
    	void Set( arr<double> &tt_new, arr<double> &ee_new, arr<double> &bb_new, arr<double> &te_new);
		
   		/*! Sets all components (6). */
    	void Set( arr<double> &tt_new, arr<double> &ee_new, arr<double> &bb_new, arr<double> &te_new, arr<double> &tb_new, arr<double> &eb_new );
    	
    	void read( char *Name );
    	void write( char *Name );
    	
    	/* Smooths the spectrum with a Gaussian beam. \a fwhm is given in radian.
       	\note This is only implememted for 1 and 4 spectra so far. */
    	void Smooth_with_Gauss( double fwhm );
};

// ===========================

Healpix_PowSpec::Healpix_PowSpec( int nspec, int lmax )
{
	num_specs = nspec;
	
	planck_assert( (num_specs==1) || (num_specs==4) || (num_specs==6), "wrong number of spectrums" );
	
	tt_.alloc( lmax+1 );
    
    if( num_specs > 1 )
    {
    	ee_.alloc( lmax+1 );
        bb_.alloc( lmax+1 );
        te_.alloc( lmax+1 );
    }
    
    if( num_specs > 4 )
	{
        tb_.alloc( lmax+1 );
        eb_.alloc( lmax+1 );
    }
}

// ===========================

void Healpix_PowSpec::dealloc()
{
	tt_.dealloc();
	ee_.dealloc();
	bb_.dealloc();
	te_.dealloc();
	tb_.dealloc();
	eb_.dealloc();
}

// ===========================

int Healpix_PowSpec::Set_TT( arr<double> &tt_new )
{
	int status;
	
	if( tt_new.size() == tt_.size() )
	{
		tt_.transfer( tt_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

int Healpix_PowSpec::Set_EE( arr<double> &ee_new )
{
	int status;
	
	if( (ee_new.size() == ee_.size()) && (num_specs > 1) )
	{
		ee_.transfer( ee_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

int Healpix_PowSpec::Set_BB( arr<double> &bb_new )
{
	int status;
	
	if( (bb_new.size() == bb_.size()) && (num_specs > 2) )
	{
		bb_.transfer( bb_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

int Healpix_PowSpec::Set_TE( arr<double> &te_new )
{
	int status;
	
	if( (te_new.size() == te_.size()) && (num_specs > 3) )
	{
		te_.transfer( te_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

int Healpix_PowSpec::Set_TB( arr<double> &tb_new )
{
	int status;
	
	if( (tb_new.size() == tb_.size()) && (num_specs > 4) )
	{
		tb_.transfer( tb_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

int Healpix_PowSpec::Set_EB( arr<double> &eb_new )
{
	int status;
	
	if( (eb_new.size() == eb_.size()) && (num_specs > 5) )
	{
		eb_.transfer( eb_new );
		
		status = 1;
	}
	else
	{
		status = -1;
	}
	
	return status;
}

// ===========================

void Healpix_PowSpec::Set(arr<double> &tt_new)
{
	dealloc();
	num_specs = 1;
	
	tt_.transfer(tt_new);

	for(int l=0; l < (int) tt_.size(); ++l)
	{
		if( tt_[l] < 0 )
		{
			cerr << "Warning: negative values in TT spectrum" << endl;
			break;
		}
	}
}

// ===========================

void Healpix_PowSpec::Set( arr<double> &tt_new, arr<double> &ee_new, arr<double> &bb_new, arr<double> &te_new )
{
	dealloc();
	
	num_specs = 4;
	
	tt_.transfer( tt_new );
	ee_.transfer( ee_new );
	bb_.transfer( bb_new );
	te_.transfer( te_new );
	
	planck_assert( (tt_.size() == ee_.size()) && (tt_.size() == bb_.size()) && (tt_.size() == te_.size()), "PowSpec::Set: size mismatch" );
	
/*	for (int l=0; l < tt_.size(); ++l)
    {
      	planck_assert (tt_[l]>=0, "negative TT spectrum at l="+std::to_string(l));
    	planck_assert (ee_[l]>=0, "negative GG spectrum at l="+dataToString(l));
    	planck_assert (bb_[l]>=0, "negative CC spectrum at l="+dataToString(l));
    	planck_assert (abs(te_[l]<=sqrt(tt_[l]*ee_[l])), "Inconsistent T, E and TxE terms at l="+dataToString(l) );
    }*/
}

// ===========================

void Healpix_PowSpec::Set( arr<double> &tt_new, arr<double> &ee_new, arr<double> &bb_new, arr<double> &te_new, arr<double> &tb_new, arr<double> &eb_new )
{
	dealloc();
	
	num_specs = 6;
	
	tt_.transfer( tt_new );
	ee_.transfer( ee_new );
	bb_.transfer( bb_new );
	te_.transfer( te_new );
	tb_.transfer( tb_new );
	eb_.transfer( eb_new );
	
	planck_assert( (tt_.size() == ee_.size()) && (tt_.size() == bb_.size()) && (tt_.size() == te_.size()) && (tt_.size() == tb_.size()) && (tt_.size() == eb_.size()), "PowSpec::Set: size mismatch" );
	
/*	for (int l=0; l < tt_.size(); ++l)
    {
    	planck_assert (tt_[l]>=0, "negative TT spectrum at l="+to_string(l));
    	planck_assert (ee_[l]>=0, "negative GG spectrum at l="+dataToString(l));
    	planck_assert (bb_[l]>=0, "negative CC spectrum at l="+dataToString(l));
    	planck_assert (abs(te_[l]<=sqrt(tt_[l]*ee_[l])), "Inconsistent T, E and TxE terms at l="+dataToString(l) );
    }*/
}

// ===========================

void Healpix_PowSpec::write( char *Name )
{
	int size = tt_.size();
	
	fltarray spec( size, num_specs );
	
	for(int l=0; l < size; l++)
	{
		if( num_specs == 1 )
		{
			spec( l ) = tt(l);
		}
		
		if( num_specs == 4 )
		{
			spec( l, 0 ) = tt(l);
			spec( l, 1 ) = ee(l);
			spec( l, 2 ) = bb(l);
			spec( l, 3 ) = te(l);
		}
		
		if( num_specs == 6 )
		{
			spec( l, 0 ) = tt(l);
			spec( l, 1 ) = ee(l);
			spec( l, 2 ) = bb(l);
			spec( l, 3 ) = te(l);
			spec( l, 4 ) = tb(l);
			spec( l, 5 ) = eb(l);
		}
	}
	char * NameFits= fitsname(Name);//FCS added to prevent memory leak
    fits_write_fltarr(NameFits, spec );//File read using MRS/astron/pro/fits_read.pro	FITS_READ, 'filename.fits', spec
	free(NameFits);
}

// ===========================

void Healpix_PowSpec::read( char *Name )

{
	fltarray spec;
	char * NameFits= fitsname(Name);//FCS added to prevent memory leak
	fits_read_fltarr(NameFits, spec );
	free(NameFits);
	int size = spec.nx();
	num_specs = spec.ny();
	
	tt_.alloc( size );
	
	
	if( num_specs > 1 )
	{
		ee_.alloc( size );
		
        bb_.alloc( size );
        
        te_.alloc( size );
    }
    
    if( num_specs > 4 )
    {
        tb_.alloc( size );
        
        eb_.alloc( size );
    }
	
  	for(int l=0; l < size; l++)
	{
		if( num_specs == 1 )
		{
			tt(l) = spec( l );
		}
		
		if( num_specs == 4 )
		{
			tt(l) = spec( l, 0 );
			ee(l) = spec( l, 1 );
			bb(l) = spec( l, 2 );
			te(l) = spec( l, 3 );
		}
		
		if( num_specs == 6 )
		{
			tt(l) = spec( l, 0 );
			ee(l) = spec( l, 1 );
			bb(l) = spec( l, 2 );
			te(l) = spec( l, 3 );
			tb(l) = spec( l, 4 );
			eb(l) = spec( l, 5 );
		}
	}
}

// ===========================

void Healpix_PowSpec::Smooth_with_Gauss( double fwhm )
{
	//planck_assert( num_specs<=4, "not yet implemented for num_specs>4" );
	
	double sigma = fwhm*FWHM2SIGMA;
	double fact_pol = exp(2*sigma*sigma);
	
	for(int l=0; l < (int) tt_.size(); ++l)
    {
    	double f1 = exp(-.5*l*(l+1)*sigma*sigma);
    	double f2 = f1*fact_pol;
    	
    	tt_[l] *= f1*f1;
    	
    	if( num_specs > 1 )
    	{
    		ee_[l] *= f2*f2;
    		bb_[l] *= f2*f2;
    		te_[l] *= f1*f2;
    	}
    	
    	if( num_specs > 4 )
    	{
    		tb_[l] *= f1*f2;
    		eb_[l] *= f2*f2;
    	}
    }
}

// ===========================

void extract_Healpix_PowSpec( const Alm<xcomplex<double> > &alm, Healpix_PowSpec &Spec )
{
	arr<double> tt( alm.Lmax()+1 );
	
	for(int l=0; l <= alm.Lmax(); ++l)
    {
    	tt[l] = norm( alm(l,0) );
    	
    	int limit = min( l,alm.Mmax() );
    	
    	for(int m=1; m <= limit; ++m)
    	{
    		tt[l] += 2*norm(alm(l,m));
    	}
    	
    	tt[l] /= (2*l+1);
    }
  
  Spec.Set(tt);
}
  
// ===========================

void extract_cross_Healpix_PowSpec( const Alm<xcomplex<double> > &alm1, const Alm<xcomplex<double> > &alm2, Healpix_PowSpec &Spec )
{
	planck_assert( alm1.conformable(alm2), "extract_cross_Healpix_PowSpec: a_lms are not conformable" );
	
	arr<double> tt( alm1.Lmax()+1 );
	
	for(int l=0; l <= alm1.Lmax(); ++l)
    {
    	tt[l] = real(alm1(l,0))*real(alm2(l,0));
    	
    	int limit = min( l,alm1.Mmax() );
    	
    	for(int m=1; m <= limit; ++m)
    	{
    		tt[l] += 2 * (real(alm1(l,m))* real(alm2(l,m)) + imag(alm1(l,m))* imag(alm2(l,m)));
    	}
    	
    	tt[l] /= (2*l+1);
    }

	Spec.Set(tt);
}

// ===========================

void extract_Healpix_PowSpec_pola( const Alm<xcomplex<double> > &almT, const Alm<xcomplex<double> > &almE, const Alm<xcomplex<double> > &almB, Healpix_PowSpec &Spec )
{
	planck_assert( almT.conformable(almE) && almT.conformable(almB), "extract_powspec: a_lms are not conformable" );
	
	int lmax = almT.Lmax();
	
	arr<double> tt(lmax+1), ee(lmax+1), bb(lmax+1), te(lmax+1);
	
	for(int l=0; l <= lmax; ++l)
    {
    	tt[l] = norm( almT(l,0) );
    	ee[l] = norm( almE(l,0) );
    	bb[l] = norm( almB(l,0) );
    	te[l] = real( almT(l,0)*conj( almE(l,0) ) );
    	
    	int limit = min( l,almT.Mmax() );
    	
    	for(int m=1; m <= limit; ++m)
    	{
    		tt[l] += 2*norm( almT(l,m) );
    		ee[l] += 2*norm( almE(l,m) );
    		bb[l] += 2*norm( almB(l,m) );
    		te[l] += 2* real( almT(l,m)*conj( almE(l,m) ) );
    	}
    	
    	tt[l] /= (2*l+1);
    	ee[l] /= (2*l+1);
    	bb[l] /= (2*l+1);
    	te[l] /= (2*l+1);
    }
    
    Spec.Set(tt,ee,bb,te);
}
  
// ===========================

void extract_Healpix_PowSpec_pola_all( const Alm<xcomplex<double> > &almT, const Alm<xcomplex<double> > &almE, const Alm<xcomplex<double> > &almB, Healpix_PowSpec &Spec )
{
	planck_assert( almT.conformable(almE) && almT.conformable(almB), "extract_powspec: a_lms are not conformable" );
	
	int lmax = almT.Lmax();
	
	arr<double> tt(lmax+1), ee(lmax+1), bb(lmax+1), te(lmax+1), tb(lmax+1), eb(lmax+1);
	
	for(int l=0; l <= lmax; ++l)
    {
    	tt[l] = norm( almT(l,0) );
    	ee[l] = norm( almE(l,0) );
    	bb[l] = norm( almB(l,0) );
    	
    	te[l] = real( almT(l,0)*conj( almE(l,0) ) );
    	tb[l] = real( almT(l,0)*conj( almB(l,0) ) );
    	eb[l] = real( almE(l,0)*conj( almB(l,0) ) );
    	
    	int limit = min( l,almT.Mmax() );
    	
    	for(int m=1; m <= limit; ++m)
    	{
    		tt[l] += 2*norm( almT(l,m) );
    		ee[l] += 2*norm( almE(l,m) );
    		bb[l] += 2*norm( almB(l,m) );
    		
    		te[l] += 2* real( almT(l,m)*conj( almE(l,m) ) );
    		tb[l] += 2* real( almT(l,m)*conj( almB(l,m) ) );
    		eb[l] += 2* real( almE(l,m)*conj( almB(l,m) ) );
    	}
    	
    	tt[l] /= (2*l+1);
    	ee[l] /= (2*l+1);
    	bb[l] /= (2*l+1);
    	
    	te[l] /= (2*l+1);
    	tb[l] /= (2*l+1);
    	eb[l] /= (2*l+1);
    }
    
    Spec.Set( tt, ee, bb, te, tb, eb );
}

// ===========================

void create_alm_Healpix_PowSpec( const Healpix_PowSpec &Spec, Alm<xcomplex<double> > &alm, planck_rng &rng )
{
	int lmax = alm.Lmax();
	int mmax = alm.Mmax();
	const double hsqrt2 = 1/sqrt(2.);
	
	for(int l=0; l <= lmax; ++l)
    {
    	double rms_tt = sqrt( Spec.tt(l) );
    	double zeta1_r = rng.rand_gauss();
    	alm(l,0) = zeta1_r * rms_tt;
    	
    	for(int m=1; m <= min(l,mmax); ++m)
    	{
    		zeta1_r = rng.rand_gauss()*hsqrt2;
    		double zeta1_i = rng.rand_gauss()*hsqrt2;
    		
    		alm(l,m) = xcomplex<REAL>( zeta1_r*rms_tt, zeta1_i*rms_tt );
    	}
    }
}

// ===========================

void create_alm_pol_Healpix_PowSpec( const Healpix_PowSpec &Spec, Alm<xcomplex<double> > &almT, Alm<xcomplex<double> > &almE, Alm<xcomplex<double> > &almB, planck_rng &rng )
{
	int lmax = almT.Lmax();
	int mmax = almT.Mmax();
	const double hsqrt2 = 1/sqrt(2.);
	
	for(int l=0; l <= lmax; ++l)
    {
    	double rms_tt=0, rms_g1=0;
    	
    	if( Spec.tt(l) != 0 )
    	{
    		rms_tt = sqrt( Spec.tt(l) );
    		rms_g1 = Spec.te(l)/rms_tt;
    	}
    	
    	double zeta1_r = rng.rand_gauss();
    	almT(l,0) = zeta1_r * rms_tt;
    	almE(l,0) = zeta1_r * rms_g1;
    	
    	for(int m=1; m <= min(l,mmax); ++m)
    	{
    		zeta1_r = rng.rand_gauss()*hsqrt2;
    		double zeta1_i = rng.rand_gauss()*hsqrt2;
    		
    		almT(l,m) = xcomplex<REAL>(  zeta1_r*rms_tt, zeta1_i*rms_tt );
    		almE(l,m) = xcomplex<REAL>( zeta1_r*rms_g1, zeta1_i*rms_g1 );
    	}
    }
    
    for(int l=0; l <= lmax; ++l)
    {
    	double rms_g2 = 0;
    	double rms_cc = 0;
    	
    	if( Spec.tt(l) != 0 )
    	{
    		rms_g2 = Spec.ee(l) - (Spec.te(l)/Spec.tt(l))*Spec.te(l);
    		if( rms_g2 <= 0 )
    		{
    			// planck_assert( abs( rms_g2 ) <= 1e-8*abs( Spec.ee(l) ), string("Inconsistent TT, EE and TE spectra at l=")+dataToString(l) );
    			rms_g2 = 0;
    		}
    		rms_g2 = sqrt( rms_g2 );
    		rms_cc = sqrt( Spec.bb(l) );
    	}
    	
    	almE(l,0) += rng.rand_gauss()*rms_g2;
    	almB(l,0)  = rng.rand_gauss()*rms_cc;
    	
    	for(int m=1; m <= min(l,mmax); ++m)
    	{
    		double zeta2_r = rng.rand_gauss()*hsqrt2;
    		double zeta2_i = rng.rand_gauss()*hsqrt2;
    		double zeta3_r = rng.rand_gauss()*hsqrt2;
    		double zeta3_i = rng.rand_gauss()*hsqrt2;
    		
    		almE(l,m) += xcomplex<double> (zeta2_r*rms_g2,zeta2_i*rms_g2);
    		almB(l,m) = xcomplex<REAL>(  zeta3_r*rms_cc,zeta3_i*rms_cc );
    	}
    }
}

// ===========================
// ===========================

void mrsp_alloc_powspec(Healpix_PowSpec & PS, int lmax)
{
	PS.num_specs = 6;
	
	planck_assert( (PS.num_specs==1) || (PS.num_specs==4) || (PS.num_specs==6), "wrong number of spectrums" );
	
	PS.tt_.alloc( lmax+1 );
    
    if( PS.num_specs > 1 )
    {
    	PS.ee_.alloc( lmax+1 );
        PS.bb_.alloc( lmax+1 );
        PS.te_.alloc( lmax+1 );
    }
    
    if( PS.num_specs > 4 )
	{
        PS.tb_.alloc( lmax+1 );
        PS.eb_.alloc( lmax+1 );
    }
}

//===========================
/*
void mrsp_alloc_powspec(Healpix_PowSpec & PS, fltarray &Cl)
{
    int Lmax = Cl.nx()-1;
    arr<double> tt(Lmax+1);
    for (int l=0; l<= Lmax; ++l) tt[l] = Cl(l);
        PS.Set(tt);
}
*/

//===========================
/*
void mrsp_alloc_powspec(Healpix_PowSpec & PS, dblarray &Cl)
{
    int Lmax = Cl.nx()-1;
    arr<double> tt(Lmax+1);
    for (int l=0; l<=Lmax; ++l) tt[l] = Cl(l);
    PS.Set(tt);
}
*/
//===========================

#endif
