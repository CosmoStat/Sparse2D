/**************************************************************

	Compute map2alm_pol and alm2map_pol for QU and EB only
	
***************************************************************/

#ifdef UPDATE_HEALPIX3_60_SOLVED

#include <vector>
//#include "alm_healpix_tools.h"
//#include "alm_map_tools.h"

#include "alm.h"
#include "arr.h"
// #include "fftpack_support.h"
#include "ylmgen.h"
#include "healpix_map.h"
#include "xcomplex.h"

using namespace std;

/*! A class holding information about a ring of pixels in a spherical map. */

class ringinfo
{
	public:
		double theta, phi0, weight, cth, sth;
		int nph, ofs;

	ringinfo()
      : nph(0) {}
    //! Constructs a \a ringinfo object.
    //  \param theta_ colatitude of the ring (in radian)
    //  \param phi0_ longitude of the first pixel in the ring (in radian)
    //  \param weight_ weighting factor for all pixels in the ring. This is typically the surface of a pixel in sterad.
    //  \note \a weight_ is only needed for map analysis, not synthesis.
    //  \param nph_ number of pixels in the ring
    //  \param ofs_ index of the first ring pixel in the total map array (counting from zero)
	ringinfo (double theta_, double phi0_, double weight_, int nph_, int ofs_)
      : theta(theta_), phi0(phi0_), weight(weight_),
        cth(cos(theta)), sth(sin(theta)), nph(nph_), ofs(ofs_)
      {}
};

/*! A class holding information about a ring pair in a spherical map. */

class ringpair
{
	public:
		ringinfo r1, r2;

    //! Initialize the object with the ring described by \a info. The second ring is left empty.
	ringpair (const ringinfo &info)
      : r1(info) {}
    //! Initialize the object with the rings described by \a info1 and \a info2.
    // \note The colatitude of \a info2 must be \f$\pi\f$ minus the colatitude of \a info1.
	ringpair (const ringinfo &info1,const ringinfo &info2)
      : r1(info1), r2(info2)
	{
		planck_assert( approx( r1.theta, pi-r2.theta, 1e-10 ), "invalid ringpair" );
	}
};

namespace {

struct info_comparator
{
	inline bool operator()( const ringinfo &a, const ringinfo &b ){ return a.sth<b.sth; }
};

struct pair_comparator
{
	inline bool operator()(const ringpair &a, const ringpair &b)
    {
		if( a.r1.nph == b.r1.nph )
		{
			return a.r1.phi0<b.r1.phi0;
		}
		else
		{
			return a.r1.nph<b.r1.nph;
		}
    }
};

void init_lam_fact_1d( int m, arr<double> &lam_fact )
{
	for( int l=m; l < (int) lam_fact.size(); ++l)
	{
		lam_fact[l] = (l<2) ? 0. : 2*sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
	}
}

void init_lam_fact_deriv_1d( int m, arr<double> &lam_fact )
{
	lam_fact[m]=0;

	for( int l=m+1; l < (int) lam_fact.size(); ++l)
	{
		lam_fact[l] = sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
	}
}

void init_normal_l( arr<double> &normal_l )
{
	for( int l=0; l < (int) normal_l.size(); ++l)
	{
		normal_l[l] = (l<2) ? 0. : sqrt(1./((l+2.)*(l+1.)*l*(l-1.)));
	}
}

void get_chunk_info( int nrings, int &nchunks, int &chunksize )
{
	nchunks = nrings/max(100,nrings/10) + 1;
	chunksize = ( nrings + nchunks -1 )/nchunks;
}

class ringhelper
{
	private:
		double phi0_;
		arr<xcomplex<double> > shiftarr, work;
		rfft plan;
		bool norot;

		void update( int nph, int mmax, double phi0 )
		{
			norot = ( abs(phi0) < 1e-14 );
			if( !norot )
			{
				if( ( mmax != (int) shiftarr.size()-1 ) || ( !approx( phi0, phi0_, 1e-12 ) ) )
				{
					shiftarr.alloc(mmax+1);
					phi0_ = phi0;

					for( int m=0; m <= mmax; ++m)
					{
						shiftarr[m] = xcomplex<REAL> (cos(m*phi0), sin(m*phi0) );
					}
				}
			}

			if( nph != (int) plan.size() )
			{
				plan.Set(nph);
			}
			if( nph > (int) work.size() )
			{
				work.alloc(2*nph);
			}
		}

	public:
		ringhelper() : phi0_(0), norot(true) {}

		template<typename T> void phase2ring( int nph, int mmax, double phi0, const xcomplex<double> *phase, T *ring )
		{
			update( nph, mmax, phi0 );

			for( int m=1; m < nph; ++m )
			{
				work[m]=0;
			}
			work[0] = phase[0];

			if(norot)
			{
				for( int m=1; m <= mmax; ++m )
				{
					work[m%nph] += phase[m];
					work[nph-1-((m-1)%nph)] += conj(phase[m]);
				}
			}
			else
			{
				for( int m=1; m <= mmax; ++m )
				{
					xcomplex<double> tmp = phase[m]*shiftarr[m];
					work[m%nph] += tmp;
					work[nph-1-((m-1)%nph)] += conj(tmp);
				}
			}

			plan.backward_c( work );

			for( int m=0; m < nph; ++m )
			{
				ring[m] = real(work[m]);
			}
		}

		template<typename T> void phase2ring( int mmax, const xcomplex<double> *phase, const ringinfo &info, T *data )
		{
			if( info.nph > 0 )
			{
				phase2ring( info.nph, mmax, info.phi0, phase, data+info.ofs );
			}
		}

		template<typename T> void phase2pair( int mmax, const xcomplex<double> *phase1, const xcomplex<double> *phase2, const ringpair &pair, T *data )
		{
			phase2ring( mmax, phase1, pair.r1, data );
			phase2ring( mmax, phase2, pair.r2, data );
		}

		template<typename T> void ring2phase( int nph, int mmax, double phi0, double weight, const T *ring, xcomplex<double> *phase )
		{
			update( nph, mmax, -phi0 );

			for( int m=0; m < nph; ++m )
			{
				work[m] = ring[m]*weight;
			}
			
			plan.forward_c( work );

			if( norot )
			{
				for( int m=0; m <= mmax; ++m )
				{
					phase[m] = work[m%nph];
				}
			}
			else
			{
				for( int m=0; m <= mmax; ++m )
				{
					phase[m] = work[m%nph]*shiftarr[m];
				}
			}
		}

		template<typename T> void ring2phase( int mmax, const ringinfo &info, const T *data, xcomplex<double> *phase )
		{
			if( info.nph > 0 )
			{
				ring2phase( info.nph, mmax, info.phi0, info.weight, data+info.ofs, phase );
			}
		}

		template<typename T> void pair2phase( int mmax, const ringpair &pair, const T *data, xcomplex<double> *phase1, xcomplex<double> *phase2 )
		{
			ring2phase( mmax, pair.r1, data, phase1 );
			ring2phase( mmax, pair.r2, data, phase2 );
		}
};

void healpix2ringpairs( const Healpix_Base &base, const arr<double> &weight, std::vector<ringpair> &pair )
{
	pair.clear();
	int startpix, ringpix;
	double theta, wgt, phi0;
	bool shifted;
	int nside = base.Nside();

	for(int m=0; m < 2*nside-1; ++m)
    {
    	base.get_ring_info2( m+1, startpix, ringpix, theta, shifted );
		wgt = weight[m]*fourpi/base.Npix();
    	phi0 = shifted ? pi/ringpix : 0;
    	pair.push_back( ringpair( ringinfo( theta, phi0, wgt, ringpix, startpix ), ringinfo( pi-theta, phi0, wgt, ringpix, base.Npix()-startpix-ringpix ) ) );
    }

	base.get_ring_info2( 2*nside, startpix, ringpix, theta, shifted );
	wgt = weight[2*nside-1]*fourpi/base.Npix();
	phi0 = shifted ? pi/ringpix : 0;
	pair.push_back( ringpair( ringinfo( theta, phi0, wgt, ringpix, startpix ) ) );
}

void healpix2ringpairs( const Healpix_Base &base, std::vector<ringpair> &pair )
{
	arr<double> wgt( 2*base.Nside() );
	wgt.fill(1);
	healpix2ringpairs( base, wgt, pair );
}

} // namespace

/*
void info2pair( const std::vector<ringpair> &info, std::vector<ringpair> &pair )
{
	pair.clear();
	
	vector<ringinfo> info2=info;
	sort( info2.begin(), info2.end(), info_comparator() );
	
	unsigned int pos=0;
	while( pos < info2.size()-1 )
    {
		if( approx( info2[pos].cth, -info2[pos+1].cth, 1e-12 ) )
		{
			pair.push_back( ringpair( info2[pos], info2[pos+1] ) );
			pos += 2;
		}
		else
		{
			pair.push_back( ringpair( info2[pos] ) );
			++pos;
		}
	}

	if( pos < info2.size() )
	{
		pair.push_back(info2[pos]);
	}

	sort( pair.begin(), pair.end(), pair_comparator() );
}
*/

/***************************************************************************************************************************/

//void alm2map_pol_QU( /*const Alm<xcomplex<T> > &almT,*/ const Alm<xcomplex<T> > &almE, const Alm<xcomplex<T> > &almB, const vector<ringpair> &pair, /* T *mapT,*/ T *mapQ, T *mapU )
template<typename T> void alm2map_pol_QU( const Alm<xcomplex<T> > &almE, const Alm<xcomplex<T> > &almB, const std::vector<ringpair> &pair, T *mapQ, T *mapU )
{
	int lmax = almE.Lmax();
	int mmax = almE.Mmax();

	planck_assert( almE.conformable(almB), "alm2map_pol: a_lm are not conformable" );

	arr<double> normal_l (lmax+1);
	init_normal_l( normal_l );

	int nchunks, chunksize;
	get_chunk_info( pair.size(), nchunks, chunksize );

	arr2<xcomplex<double> > phas1Q(chunksize,mmax+1), phas2Q(chunksize,mmax+1), phas1U(chunksize,mmax+1), phas2U(chunksize,mmax+1);//phas1T(chunksize,mmax+1), phas2T(chunksize,mmax+1)

	for(int chunk=0; chunk < nchunks; ++chunk)
    {
		int llim = chunk*chunksize, ulim = min( llim+chunksize, int(pair.size()) );

		#pragma omp parallel
		{
			Ylmgen generator( lmax, mmax, 1e-30 );
			arr<double> Ylm;
			arr<double> lam_fact (lmax+1);
			arr<xcomplex<double>[2]> alm_tmp(lmax+1);
			int m;

			#pragma omp for schedule(dynamic,1)
			for( m=0; m <= mmax; ++m )
			{
				int m2 = m*m;
				init_lam_fact_1d( m, lam_fact );
				
				for(int l=m; l <= lmax; ++l)
				{
					//alm_tmp[l][0] = almT(l,m);
					// alm_tmp[l][0].re = almE(l,m).re *(-normal_l[l]);
                    // alm_tmp[l][0].im = almE(l,m).im *(-normal_l[l]);
                    alm_tmp[l][0] = almE(l,m) * xcomplex<REAL>(-normal_l[l],0.);
					// alm_tmp[l][1].re = almB(l,m).re*(-normal_l[l]);
                    //alm_tmp[l][1].im = almB(l,m).im*(-normal_l[l]);
                    alm_tmp[l][1] = almB(l,m) * xcomplex<REAL>(-normal_l[l],0.);
				}
				
				for(int ith=0; ith < ulim-llim; ++ith)
				{
					double cth=pair[ith+llim].r1.cth, sth=pair[ith+llim].r1.sth;
					int l;
					generator.get_Ylm( cth, sth, m, Ylm, l );

					if( l <= lmax )
					{
						double one_on_s2 = 1/(sth*sth);
						double c_on_s2 = cth * one_on_s2;
						double two_on_s2 = 2*one_on_s2;
						double m_on_s2 = m*one_on_s2;
						double twocth = 2*cth;

						if( pair[ith+llim].r2.nph > 0 )
						{
							xcomplex<double> Q1=0, Q2=0, U1=0, U2=0;//T1=0, T2=0
							double lam_lm = 0;

							if( (l-m)&1 )
							{
								//ALM2MAP_POL_MACRO_QU(Q2,Q1,U2,U1)
								double lam_lm1m = lam_lm;
								lam_lm = Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								xcomplex<double> lambda_w = xcomplex<double>(a_w*lam_lm - t1*c_on_s2,0.);
								xcomplex<double> lambda_x = xcomplex<double>(m_on_s2 * (a_x-t1),0.);
								
								// Q2.re += alm_tmp[l][0].re*lambda_w;
								// Q2.im += alm_tmp[l][0].im*lambda_w;
                                Q2 += alm_tmp[l][0] * lambda_w;
								
								// U2.re -= alm_tmp[l][1].re*lambda_w;
								// U2.im -= alm_tmp[l][1].im*lambda_w;
                                U2 -= alm_tmp[l][1] * lambda_w;
								
								// Q1.re -= alm_tmp[l][1].im*lambda_x;
								// Q1.im += alm_tmp[l][1].re*lambda_x;
                                Q1 += xcomplex<REAL>( - imag(alm_tmp[l][1])  lambda_x), real(alm_tmp[l][1] * lambda_x);
								
								// U1.re -= alm_tmp[l][0].im*lambda_x;
								// U1.im += alm_tmp[l][0].re*lambda_x;
                                U1 += xcomplex<REAL>( - imag(alm_tmp[l][0] *  lambda_x), real(alm_tmp[l][0]*  lambda_x));
								++l;
							}						
							
							for( ; l < lmax ; )
							{
								//ALM2MAP_POL_MACRO_QU(Q1,Q2,U1,U2)
								double lam_lm1m = lam_lm;
								lam_lm = Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								xcomplex<double>  lambda_w = xcomplex<double> (a_w*lam_lm - t1*c_on_s2,0.);
								xcomplex<double>  lambda_x = xcomplex<double> (m_on_s2 * (a_x-t1),0.);
								
								// Q1.re += alm_tmp[l][0].re*lambda_w;
								// Q1.im += alm_tmp[l][0].im*lambda_w;
                                Q1 += alm_tmp[l][0] * lambda_w;
								
								// U1.re -= alm_tmp[l][1].re*lambda_w;
								// U1.im -= alm_tmp[l][1].im*lambda_w;
                                U1 += alm_tmp[l][1] * lambda_w;
								
								// Q2.re -= alm_tmp[l][1].im*lambda_x;
								// Q2.im += alm_tmp[l][1].re*lambda_x;
                                Q2 += xcomplex<REAL>( - imag(alm_tmp[l][1] * lambda_x), real(alm_tmp[l][1]* lambda_x));

								// U2.re -= alm_tmp[l][0].im*lambda_x;
								// U2.im += alm_tmp[l][0].re*lambda_x;
                                U2 += xcomplex<REAL>( - imag(alm_tmp[l][0] * lambda_x), real(alm_tmp[l][0]* lambda_x));

								++l;
								
								//ALM2MAP_POL_MACRO_QU(Q2,Q1,U2,U1)
								lam_lm1m = lam_lm;
								lam_lm = Ylm[l];
								
								t1  = lam_lm1m*lam_fact[l];
								a_w = (l-m2)*two_on_s2 + l*(l-1);
								a_x = twocth*(l-1)*lam_lm;
								lambda_w = a_w*lam_lm - t1*c_on_s2;
								lambda_x = m_on_s2 * (a_x-t1);
								
								// Q2.re += alm_tmp[l][0].re*lambda_w;
								// Q2.im += alm_tmp[l][0].im*lambda_w;
                                Q2 += alm_tmp[l][0] * xcomplex<REAL>(lambda_w,0.);
								
								// U2.re -= alm_tmp[l][1].re*lambda_w;
								// U2.im -= alm_tmp[l][1].im*lambda_w;
                                U2 -= alm_tmp[l][1] * xcomplex<REAL>(lambda_w,0.);

								// Q1.re -= alm_tmp[l][1].im*lambda_x;
								// Q1.im += alm_tmp[l][1].re*lambda_x;
                                Q1 += xcomplex<REAL>( - imag(alm_tmp[l][1] * lambda_x), real(alm_tmp[l][1] * lambda_x));
								
								// U1.re -= alm_tmp[l][0].im*lambda_x;
								// U1.im += alm_tmp[l][0].re*lambda_x;
                                U1 += xcomplex<REAL>( - imag(alm_tmp[l][0] *  lambda_x), real(alm_tmp[l][0]*  lambda_x));

								++l;
							}

							if( l == lmax )
							{
								//ALM2MAP_POL_MACRO_QU(Q1,Q2,U1,U2)
								double lam_lm1m = lam_lm;
								lam_lm = Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								xcomplex<double> lambda_w = xcomplex<double>(a_w*lam_lm - t1*c_on_s2,0.);
								xcomplex<double> lambda_x = xcomplex<double>(m_on_s2 * (a_x-t1),0.);
								
								// Q1.re += alm_tmp[l][0].re*lambda_w;
								// Q1.im += alm_tmp[l][0].im*lambda_w;
                                Q1 += alm_tmp[l][0] *  lambda_w;

								// U1.re -= alm_tmp[l][1].re*lambda_w;
								// U1.im -= alm_tmp[l][1].im*lambda_w;
                                U1 -= alm_tmp[l][1] *  lambda_w;
								
								// Q2.re -= alm_tmp[l][1].im*lambda_x;
								// Q2.im += alm_tmp[l][1].re*lambda_x;
                                Q2 += xcomplex<REAL>( - imag(alm_tmp[l][1]) *  lambda_x, real(alm_tmp[l][1])*  lambda_x);

								// U2.re -= alm_tmp[l][0].im*lambda_x;
								// U2.im += alm_tmp[l][0].re*lambda_x;
                                U2 += xcomplex<REAL>( - imag(alm_tmp[l][0]) *  lambda_x, real(alm_tmp[l][0]) *  lambda_x);
								++l;
							}

							//phas1T[ith][m] = T1+T2;
							//phas2T[ith][m] = T1-T2;
							phas1Q[ith][m] =-Q1-Q2;
							phas2Q[ith][m] =-Q1+Q2;
							phas1U[ith][m] = U1+U2;
							phas2U[ith][m] = U1-U2;
						}
						else
						{
							xcomplex<double> Q1=0, U1=0;//T1=0,
							double lam_lm = 0;

							for( ; l <= lmax; )
							{
								//ALM2MAP_POL_MACRO_QU(Q1,Q1,U1,U1)
								double lam_lm1m = lam_lm;
								lam_lm = Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								xcomplex<REAL> lambda_w = xcomplex<REAL>(a_w*lam_lm - t1*c_on_s2,0.);
								xcomplex<REAL> lambda_x = xcomplex<REAL>(m_on_s2 * (a_x-t1),0.);
								
								// Q1.re += alm_tmp[l][0].re*lambda_w;
								// Q1.im += alm_tmp[l][0].im*lambda_w;
                                Q1 += alm_tmp[l][0] *  lambda_w;

								// U1.re -= alm_tmp[l][1].re*lambda_w;
								// U1.im -= alm_tmp[l][1].im*lambda_w;
                                U1 -= alm_tmp[l][1] *  lambda_w;

                                // JLS: is there a bug here: Q1 and U1 again.
								// Q1.re -= alm_tmp[l][1].im*lambda_x;
								// Q1.im += alm_tmp[l][1].re*lambda_x;
                                Q1 += xcomplex<REAL>( - imag(alm_tmp[l][1]) *  lambda_x, real(alm_tmp[l][1])*  lambda_x);

								// U1.re -= alm_tmp[l][0].im*lambda_x;
								// U1.im += alm_tmp[l][0].re*lambda_x;
                                U1 += xcomplex<REAL>( - imag(alm_tmp[l][0]) *  lambda_x, real(alm_tmp[l][0])*  lambda_x);
								++l;
							}

							//phas1T[ith][m] = T1;
							phas1Q[ith][m] =-Q1;
							phas1U[ith][m] = U1;
						}
					}
					else
					{
						//phas1T[ith][m] = phas2T[ith][m] = 0;
						phas1Q[ith][m] = phas2Q[ith][m] = 0;
						phas1U[ith][m] = phas2U[ith][m] = 0;
					}
				}
			}
		} // end of parallel region

		#pragma omp parallel
		{
			ringhelper helper;
			int ith;

			#pragma omp for schedule(dynamic,1)
			for( ith=llim; ith < ulim; ++ith )
			{
				//helper.phase2pair( mmax, phas1T[ith-llim], phas2T[ith-llim], pair[ith], mapT );
				helper.phase2pair( mmax, phas1Q[ith-llim], phas2Q[ith-llim], pair[ith], mapQ );
				helper.phase2pair( mmax, phas1U[ith-llim], phas2U[ith-llim], pair[ith], mapU );
			}
		} // end of parallel region
	}
}

template void alm2map_pol_QU( const Alm<xcomplex<float> > &almE, const Alm<xcomplex<float> > &almB, const std::vector<ringpair> &pair, float *mapQ, float *mapU );
template void alm2map_pol_QU( const Alm<xcomplex<double> > &almE, const Alm<xcomplex<double> > &almB, const std::vector<ringpair> &pair, double *mapQ, double *mapU );

//void map2alm_pol_QU( const vector<ringpair> &pair, /*const T *mapT,*/ const T *mapQ, const T *mapU, /*Alm<xcomplex<T> > &almT,*/ Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, bool add_alm )
template<typename T> void map2alm_pol_QU( const std::vector<ringpair> &pair, const T *mapQ, const T *mapU, Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, bool add_alm )
{
	planck_assert( almE.conformable(almB), "map2alm_pol: a_lm are not conformable" );

	int lmax = almE.Lmax(), mmax = almE.Mmax();

	arr<double> normal_l (lmax+1);
	init_normal_l( normal_l );

	int nchunks, chunksize;
	get_chunk_info( pair.size(), nchunks, chunksize );

	arr2<xcomplex<double> > phas1Q(chunksize,mmax+1), phas2Q(chunksize,mmax+1), phas1U(chunksize,mmax+1), phas2U(chunksize,mmax+1);//phas1T(chunksize,mmax+1), phas2T(chunksize,mmax+1)

	if( !add_alm )
    { 
    	//almT.SetToZero();
    	almE.SetToZero();
    	almB.SetToZero();
    }

	for( int chunk=0; chunk < nchunks; ++chunk )
    {
		int llim = chunk*chunksize, ulim = min( llim+chunksize, int(pair.size()) );

		#pragma omp parallel
		{
			ringhelper helper;

			int ith;
			#pragma omp for schedule(dynamic,1)
			for( ith=llim; ith < ulim; ++ith)
			{
				//helper.pair2phase( mmax, pair[ith], mapT, phas1T[ith-llim], phas2T[ith-llim] );
				helper.pair2phase( mmax, pair[ith], mapQ, phas1Q[ith-llim], phas2Q[ith-llim] );
				helper.pair2phase( mmax, pair[ith], mapU, phas1U[ith-llim], phas2U[ith-llim] );
			}
		} // end of parallel region

		#pragma omp parallel
		{
			Ylmgen generator( lmax, mmax, 1e-30 );
			arr<double> Ylm;
			arr<double> lam_fact(lmax+1);
			arr<xcomplex<double>[2] > alm_tmp(lmax+1);
			int m;

			#pragma omp for schedule(dynamic,1)
			for( m=0; m <= mmax; ++m )
			{
				init_lam_fact_1d( m, lam_fact );

				for( int l=m; l < (int) alm_tmp.size(); ++l )
				{
					alm_tmp[l][0] = alm_tmp[l][1] = 0;
				}

				for( int ith=0; ith < ulim-llim; ++ith )
				{
					int l;
					double cth=pair[ith+llim].r1.cth, sth=pair[ith+llim].r1.sth;
					generator.get_Ylm( cth, sth, m, Ylm, l );

					if( l <= lmax )
					{
						double one_on_s2 = 1/(sth*sth);
						double c_on_s2 = cth * one_on_s2;
						double two_on_s2 = 2*one_on_s2;
						double twocth = 2*cth;
						int m2 = m*m;
						double m_on_s2 = m*one_on_s2;

						if( pair[ith+llim].r2.nph > 0 )
						{
							xcomplex<double> Q1 = phas1Q[ith][m]+phas2Q[ith][m], Q2 = phas1Q[ith][m]-phas2Q[ith][m], U1 = phas1U[ith][m]+phas2U[ith][m], U2 = phas1U[ith][m]-phas2U[ith][m];
							//T1 = phas1T[ith][m]+phas2T[ith][m], T2 = phas1T[ith][m]-phas2T[ith][m]

							double lam_lm = 0;
							if( (l-m)&1 )
							{
								//MAP2ALM_POL_MACRO_QU(Q2,Q1,U2,U1)
								double lam_lm1m=lam_lm;
								lam_lm=Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								double lambda_w = a_w*lam_lm - t1*c_on_s2;
								double lambda_x = m_on_s2 * (a_x-t1);
								
								// alm_tmp[l][0].re += Q2.re*lambda_w - U1.im*lambda_x;
								// alm_tmp[l][0].im += Q2.im*lambda_w + U1.re*lambda_x;
                                alm_tmp[l][0] += xcomplex<REAL>( real(Q2) * (REAL) lambda_x - imag(U1) * (REAL) lambda_x, imag(Q2) * (REAL) lambda_x + real(U1) * (REAL) lambda_x);

                                alm_tmp[l][0] += Q2 * (REAL) lambda_w;
                                
								alm_tmp[l][1] +=  xcomplex<REAL>( U2.real()*lambda_w + Q1.imag()*lambda_x, U2.imag()*lambda_w - Q1.real()*lambda_x);
								++l;

							}
						
							for( ; l < lmax; )
							{
								//MAP2ALM_POL_MACRO_QU(Q1,Q2,U1,U2)
								double lam_lm1m=lam_lm;
								lam_lm=Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								double lambda_w = a_w*lam_lm - t1*c_on_s2;
								double lambda_x = m_on_s2 * (a_x-t1);
								
								// alm_tmp[l][0].re += Q1.re*lambda_w - U2.im*lambda_x;
								// alm_tmp[l][0].im += Q1.im*lambda_w + U2.re*lambda_x;
                                alm_tmp[l][0] += xcomplex<REAL>( Q1.real()*lambda_w - U2.imag()*lambda_x, Q1.imag()*lambda_w + U2.real()*lambda_x);
								// alm_tmp[l][1].re += U1.re*lambda_w + Q2.im*lambda_x;
								// alm_tmp[l][1].im += U1.im*lambda_w - Q2.re*lambda_x;
                                alm_tmp[l][1]+= xcomplex<REAL>(U1.real()*lambda_w + Q2.imag()*lambda_x,U1.imag()*lambda_w - Q2.real()*lambda_x);
								++l;
								
								//MAP2ALM_POL_MACRO_QU(Q2,Q1,U2,U1)
								lam_lm1m=lam_lm;
								lam_lm=Ylm[l];
								
								t1  = lam_lm1m*lam_fact[l];
								a_w = (l-m2)*two_on_s2 + l*(l-1);
								a_x = twocth*(l-1)*lam_lm;
								lambda_w = a_w*lam_lm - t1*c_on_s2;
								lambda_x = m_on_s2 * (a_x-t1);
								
								// alm_tmp[l][0].re += Q2.re*lambda_w - U1.im*lambda_x;
								// alm_tmp[l][0].im += Q2.im*lambda_w + U1.re*lambda_x;
                                alm_tmp[l][0] += xcomplex<REAL>(Q2.real()*lambda_w - U1.imag()*lambda_x, Q2.imag()*lambda_w + U1.real()*lambda_x);
								// alm_tmp[l][1].re += U2.re*lambda_w + Q1.im*lambda_x;
								// alm_tmp[l][1].im += U2.im*lambda_w - Q1.re*lambda_x;
                                alm_tmp[l][1] += xcomplex<REAL>(U2.real()*lambda_w + Q1.imag()*lambda_x,U2.imag()*lambda_w - Q1.real()*lambda_x);
								++l;
							}

							if( l == lmax )
							{
								//MAP2ALM_POL_MACRO_QU(Q1,Q2,U1,U2)
								double lam_lm1m=lam_lm;
								lam_lm=Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								double lambda_w = a_w*lam_lm - t1*c_on_s2;
								double lambda_x = m_on_s2 * (a_x-t1);
								
								// alm_tmp[l][0].re += Q1.re*lambda_w - U2.im*lambda_x;
								// alm_tmp[l][0].im += Q1.im*lambda_w + U2.re*lambda_x;
                                alm_tmp[l][0] += xcomplex<REAL>(Q1.real()*lambda_w - U2.imag()*lambda_x,Q1.imag()*lambda_w + U2.real()*lambda_x);
								// alm_tmp[l][1].re += U1.re*lambda_w + Q2.im*lambda_x;
								// alm_tmp[l][1].im += U1.im*lambda_w - Q2.re*lambda_x;
                                alm_tmp[l][1] += xcomplex<REAL>(U1.real()*lambda_w + Q2.imag()*lambda_x,U1.imag()*lambda_w - Q2.real()*lambda_x);
								++l;
							}
						}
						else
						{
							xcomplex<double> Q1 = phas1Q[ith][m], U1 = phas1U[ith][m];
							//T1 = phas1T[ith][m]
							double lam_lm = 0;

							for( ; l <= lmax; )
							{
								//MAP2ALM_POL_MACRO_QU(Q1,Q1,U1,U1)
								double lam_lm1m=lam_lm;
								lam_lm=Ylm[l];
								
								double t1  = lam_lm1m*lam_fact[l];
								double a_w = (l-m2)*two_on_s2 + l*(l-1);
								double a_x = twocth*(l-1)*lam_lm;
								double lambda_w = a_w*lam_lm - t1*c_on_s2;
								double lambda_x = m_on_s2 * (a_x-t1);
								
								// alm_tmp[l][0].re += Q1.re*lambda_w - U1.im*lambda_x;
								// alm_tmp[l][0].im += Q1.im*lambda_w + U1.re*lambda_x;
                                alm_tmp[l][0] += xcomplex<REAL>(Q1.real()*lambda_w - U1.imag()*lambda_x,Q1.imag()*lambda_w + U1.real()*lambda_x);
								// alm_tmp[l][1].re += U1.re*lambda_w + Q1.im*lambda_x;
								// alm_tmp[l][1].im += U1.im*lambda_w - Q1.re*lambda_x;
                                alm_tmp[l][1] += xcomplex<REAL>(U1.real()*lambda_w + Q1.imag()*lambda_x,U1.imag()*lambda_w - Q1.real()*lambda_x);
								++l;
							}
						}
					}
				}

				xcomplex<T> *palmE=almE.mstart(m), *palmB=almB.mstart(m);
				// *palmT=almT.mstart(m)

				for( int l=m; l <= lmax; ++l )
				{
					//palmT[l].re += alm_tmp[l][0].re;
					//palmT[l].im += alm_tmp[l][0].im;
					// palmE[l].re += alm_tmp[l][0].re*normal_l[l];
					// palmE[l].im += alm_tmp[l][0].im*normal_l[l];
                    palmE[l] += xcomplex<REAL>(alm_tmp[l][0].real() *normal_l[l], alm_tmp[l][0].imag() *normal_l[l]);
					// palmB[l].re += alm_tmp[l][1].re*normal_l[l];
					// palmB[l].im += alm_tmp[l][1].im*normal_l[l];
                    almB[l] += xcomplex<REAL>(alm_tmp[l][1].real()*normal_l[l], alm_tmp[l][1].imag()*normal_l[l]);
				}
			}
		} // end of parallel region
	}
}

template void map2alm_pol_QU( const std::vector<ringpair> &pair, const float *mapQ, const float *mapU, Alm<xcomplex<float> > &almE, Alm<xcomplex<float> > &almB, bool add_alm );
template void map2alm_pol_QU( const std::vector<ringpair> &pair, const double *mapQ, const double *mapU, Alm<xcomplex<double> > &almE, Alm<xcomplex<double> > &almB, bool add_alm );

/*********************************************************************************************************************************************************************/

//void alm2map_pol_QU( /*const Alm<xcomplex<T> > &almT,*/ const Alm<xcomplex<T> > &almE, const Alm<xcomplex<T> > &almB, /*Healpix_Map<T> &mapT,*/ Healpix_Map<T> &mapQ, Healpix_Map<T> &mapU )
template<typename T> void alm2map_pol_QU( const Alm<xcomplex<T> > &almE, const Alm<xcomplex<T> > &almB, Healpix_Map<T> &mapQ, Healpix_Map<T> &mapU )
{
	/*
	planck_assert( mapT.Scheme()==RING, "alm2map_pol: maps must be in RING scheme" );
	planck_assert( mapT.conformable( mapQ ) && mapT.conformable( mapU ), "alm2map_pol: maps are not conformable" );
	*/
	planck_assert( mapQ.Scheme()==RING, "alm2map_pol: maps must be in RING scheme" );
	planck_assert( mapQ.conformable( mapU ), "alm2map_pol: maps are not conformable" );
	
	std::vector<ringpair> pair;
	healpix2ringpairs( mapQ, pair );
	
	alm2map_pol_QU( almE, almB, pair, &mapQ[0], &mapU[0] );
}

template void alm2map_pol_QU( const Alm<xcomplex<float> > &almE, const Alm<xcomplex<float> > &almB, Healpix_Map<float> &mapQ, Healpix_Map<float> &mapU );
template void alm2map_pol_QU( const Alm<xcomplex<double> > &almE, const Alm<xcomplex<double> > &almB, Healpix_Map<double> &mapQ, Healpix_Map<double> &mapU );

//void map2alm_pol_QU( /*const Healpix_Map<T> &mapT,*/ const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, /*Alm<xcomplex<T> > &almT,*/ Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, const arr<double> &weight, bool add_alm )
template<typename T> void map2alm_pol_QU( const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, const arr<double> &weight, bool add_alm=false )
{
	/*
	planck_assert( mapT.Scheme()==RING, "map2alm_pol: maps must be in RING scheme" );
	planck_assert( mapT.conformable( mapQ ) && mapT.conformable( mapU ), "map2alm_pol: maps are not conformable" );
	planck_assert( weight.size() >= 2*mapT.Nside(), "map2alm_pol: at least one weight array has too few entries" );
	*/
	planck_assert( mapQ.Scheme()==RING, "map2alm_pol: maps must be in RING scheme" );
	planck_assert( mapQ.conformable( mapU ), "map2alm_pol: maps are not conformable" );
	planck_assert( weight.size() >= (unsigned long) 2*mapQ.Nside(), "map2alm_pol: at least one weight array has too few entries" );
	
	std::vector<ringpair> pair;
	healpix2ringpairs( mapQ, weight, pair );
	
	map2alm_pol_QU( pair, &mapQ[0], &mapU[0], almE, almB, add_alm );
}

template void map2alm_pol_QU( const Healpix_Map<float> &mapQ, const Healpix_Map<float> &mapU, Alm<xcomplex<float> > &almE, Alm<xcomplex<float> > &almB, const arr<double> &weight, bool add_alm );
template void map2alm_pol_QU( const Healpix_Map<double> &mapQ, const Healpix_Map<double> &mapU, Alm<xcomplex<double> > &almE, Alm<xcomplex<double> > &almB, const arr<double> &weight, bool add_alm );

//void map2alm_pol_iter_QU( /*const Healpix_Map<T> &mapT,*/ const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, /*Alm<xcomplex<T> > &almT,*/ Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, int num_iter, const arr<double> weight )
template<typename T> void map2alm_pol_iter_QU( const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, Alm<xcomplex<T> > &almE, Alm<xcomplex<T> > &almB, int num_iter, const arr<double> weight )
{
	map2alm_pol_QU( mapQ, mapU, almE, almB, weight );
	
	for(int iter=1; iter <= num_iter; ++iter)
    {
    	Healpix_Map<T> mapQ2( mapQ.Nside(), mapQ.Scheme(), SET_NSIDE ), mapU2( mapQ.Nside(), mapQ.Scheme(), SET_NSIDE );//mapT2( mapT.Nside(), mapT.Scheme(), SET_NSIDE )
    	
    	alm2map_pol_QU( almE, almB, mapQ2, mapU2 );
    	
    	for(int m=0; m < mapQ.Npix(); ++m)
    	{
    		//mapT2[m] = mapT[m]-mapT2[m];
    		mapQ2[m] = mapQ[m]-mapQ2[m];
    		mapU2[m] = mapU[m]-mapU2[m];
    	}
    	
    	map2alm_pol_QU( mapQ2, mapU2, almE, almB, weight, true );
    }
}

//template void map2alm_pol_iter_QU( const Healpix_Map<float> &mapQ, const Healpix_Map<float> &mapU, Alm<xcomplex<float> > &almE, Alm<xcomplex<float> > &almB, int num_iter, const arr<double> &weight );
//template void map2alm_pol_iter_QU( const Healpix_Map<double> &mapQ, const Healpix_Map<double> &mapU, Alm<xcomplex<double> > &almE, Alm<xcomplex<double> > &almB, int num_iter, const arr<double> &weight );
#endif
