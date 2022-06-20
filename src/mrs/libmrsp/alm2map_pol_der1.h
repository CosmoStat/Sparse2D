
  
#ifndef _POLDER_H_
#define _POLDER_H_

#include <vector>
#include "alm.h"
#include "arr.h"
// #include "fftpack_support.h"
// #include "ylmgen.h"
#include "healpix_map.h"
#include "xcomplex.h"

#ifndef DBL_MAX
	#define DBL_MAX 1.7976931348623158e+308
#endif
#define SMAXCHK 50 // maximum size of chunk (in number of ring pairs)
#define KvS 1.0

#define RSMAX 20
#define RSMIN -20

#define HPX_MXL0 40 
#define HPX_MXL1 1.35

#define LOG2LG 100

double rescale_tab[41]; // index from -20 to 20

template<typename T> void gen_mfac( long m_max, arr<T> &m_fact )
{
	long m;
	
	m_fact[0] = 1.0;
	
	for( m=1; m <= m_max; m++ )
	{
		m_fact[m] = m_fact[m-1]*sqrt( (double)(2*m+1)/(double)(2*m) );
	}
	
	for( m=0; m <= m_max; m++ )
	{
		m_fact[m] = log( inv_sqrt4pi * m_fact[m] ) * inv_ln2;
	}
}

template void gen_mfac( long m_max, arr<float> &m_fact );
template void gen_mfac( long m_max, arr<double> &m_fact );

template<typename T> void gen_normpol( long l_max, arr<T> &normal_l )
{
	long l;
	double fl,xx;
	
	normal_l[0] = 0.0;
	normal_l[1] = 0.0;
	
	for( l=2; l <= l_max; l++)
	{
		fl = (double)l;
		xx = (fl + 2.0) * (fl + 1.0) * fl * (fl - 1.0);
		normal_l[l] = sqrt( KvS / xx );
	}
}

void init_rescale()
{
	long smax;
	double FL_LARGE = pow( 2.0, (double)LOG2LG );
	double logOVFLOW = log( FL_LARGE );
	
	smax = (long)( log( DBL_MAX ) / logOVFLOW );
	
	if( smax > (RSMAX-1) )
	{
		cout << "Array rescale_tab too small" << '\n';
		cout << "smax: " << smax << " RSMAX: " << RSMAX << '\n';
	}
	
	for( long n=-20; n <= 20; n++ )
	{
		rescale_tab[n+20] = pow( FL_LARGE, (double)(n - RSMAX) );
	}
	rescale_tab[20] = 1.0;
}

void get_pixel_layout( long nside, long ith, double &cth, double &sth, long &nphi, long &startpix, long &kphi0 )
{
	long nrings;
	double dth1, dth2, dst1;
	
	nrings = 2*nside;
	
	//if (ith < 1 .or. ith> nrings) then
    //   prlong*,'ith out of bounds ',ith,1,nrings
    //   call fatal_error
    //endif
    
    dth1 = 1.0 / ( 3.0*nside*nside );
    dth2 = 2.0 / ( 3.0*nside );
    dst1 = 1.0 / ( sqrt(6.0) * nside );
    
    if( ith < nside )
    {
    	cth = 1.0 - (double)(ith*ith) * dth1;
		nphi = 4*ith;
		kphi0 = 1;
		sth = sin( 2.0 * asin( ith * dst1 ) );// sin(theta)
		startpix = 2*ith*( ith - 1 );

    }
    else
    {
    	cth = ( 2*nside - ith ) * dth2;
		nphi = 4*nside;
		kphi0 = ( ith + 1 - nside ) % 2;
		sth = sqrt( ( 1.0 - cth )*( 1.0 + cth ) );// sin(theta)
		startpix = 2*nside*( nside - 1 ) + ( ith - nside )*nphi;
    }
	
}

template<typename T> void gen_recfac( long l_max, long m, arr2<T> &recfac )
{
	double fm2, fl2;
	
	recfac.fill(0.0);
	fm2 = m*m;
	
	for( long l=m; l <= l_max; l++ )
	{
		fl2 = (l+1)*(l+1);
		recfac[0][l] = sqrt( ( 4.0 * fl2 - 1.0 ) / ( fl2 - fm2 ) );
		
		recfac[1][l] = 1.0 / recfac[0][l];
	}
}

template void gen_recfac( long l_max, long m, arr2<float> &recfac );
template void gen_recfac( long l_max, long m, arr2<double> &recfac );

template<typename T> void gen_lamfac( long l_max, long m, arr<T> &lam_fact )
{
	double fm2, fl, fl2;
	long l_start;
	
	lam_fact.fill(0.0);
	fm2 = m*m;
	
	l_start = MAX( 2, m+1 );
	for( long l=l_start; l <= l_max; l++ )
	{
		fl = l;
		fl2 = l*l;
		lam_fact[l] = sqrt( ( 2.0 * fl + 1.0 ) / ( 2.0 * fl - 1.0 ) * ( fl2 - fm2 ) );
	}
}

template void gen_lamfac( long l_max, long m, arr<float> &lam_fact );
template void gen_lamfac( long l_max, long m, arr<double> &lam_fact );

template<typename T> void gen_lamfac_der( long l_max, long m, arr<T> &lam_fact )
{
	double fm2, fl, fl2;
	long l_start;
	
	lam_fact.fill(0.0);
	fm2 = m*m;
	
	l_start = MAX( 1, m+1 );//different lower bound than pol. factor
	for( long l=l_start; l <= l_max; l++ )
	{
		fl = l;
		fl2 = l*l;
		lam_fact[l] = sqrt( ( 2.0 * fl + 1.0 ) / ( 2.0 * fl - 1.0 ) * ( fl2 - fm2 ) );
		//different normalization than polarization factor
	}
}

template void gen_lamfac_der( long l_max, long m, arr<float> &lam_fact );
template void gen_lamfac_der( long l_max, long m, arr<double> &lam_fact );

inline long l_min_ylm( long m, double sth )
{
	long lmin = m;
	
	if( HPX_MXL0 > 0 )
	{
		lmin = MAX( lmin, long( ( m - HPX_MXL0 )/( HPX_MXL1 * sth ) ) );
	}
	
	return lmin;
}

template<typename T> void do_lam_lm_pol( long lmax, long m, double cth, double sth, double mfac, arr2<T> recfac, arr<T> lam_fact, arr2<T> &lam_lm )
{
	long scalel, l, l_min, index_rstab;
	double log2val, dlog2lg, ovflow, unflow, corfac, lam_mm, lam_0, lam_1, lam_2, lam_lm1m, normal_m, fm2, fl, flm1, two_cth, one_on_s2, fm_on_s2, two_on_s2, c_on_s2, a_w, a_x, b_w;
	
	//define constants
    ovflow = rescale_tab[21];
    unflow = rescale_tab[19];
    l_min = l_min_ylm( m, sth );
    dlog2lg = (double)LOG2LG;
    
    fm2 = m*m;
    normal_m = ( 2.0 * m ) * ( 1 - m );
    two_cth = 2.0 * cth;
    one_on_s2 = 1.0 / ( sth * sth );
    fm_on_s2 = m * one_on_s2;
    two_on_s2 = 2.0 * one_on_s2;
    c_on_s2 = cth * one_on_s2;
    b_w = c_on_s2;
    
    //computes lamba_mm
    log2val = mfac + m*log( sth ) * inv_ln2;// log_2(lam_mm)
    scalel = long ( log2val / dlog2lg );
    
    index_rstab = 20 + MAX( scalel, RSMIN );
    corfac = rescale_tab[index_rstab];
    lam_mm = pow( 2.0, log2val - scalel * dlog2lg ); // rescaled lam_mm
    if( (m & 1) > 0 )
    {
    	lam_mm = -lam_mm; // negative for odd m
    }
    
    lam_lm.fill(0.0);
    
    // --- l = m ---
    lam_lm[0][m] = corfac * lam_mm;//Actual lam_mm 

    if( m >= l_min )
    {
    	//skip Ymm if too small
		lam_lm[1][m] = ( normal_m * lam_lm[0][m] ) * ( 0.5 - one_on_s2 );
		lam_lm[2][m] = ( normal_m * lam_lm[0][m] ) * c_on_s2;
    }

    // --- l > m ---
    lam_0 = 0.0;
    lam_1 = 1.0;
    lam_2 = cth * lam_1 * recfac[0][m];
    
    for( l = m+1; l <= lmax; l++ )
    {
    	//do recursion
		lam_lm1m = lam_lm[0][l-1] * lam_fact[l];// must be incremented even if not used
		lam_lm[0][l] = lam_2 * corfac * lam_mm;
		
		if( l >= l_min )
		{
			fl = l;
			flm1 = fl - 1.0;
			a_w =  two_on_s2 * ( fl - fm2 )  + flm1 * fl;
			a_x =  two_cth * flm1;
			lam_lm[1][l] = b_w * lam_lm1m - a_w * lam_lm[0][l];
			lam_lm[2][l] = fm_on_s2 * ( lam_lm1m - a_x * lam_lm[0][l] );
		}
		
		lam_0 = lam_1 * recfac[1][l-1];
		lam_1 = lam_2;
		lam_2 = ( cth * lam_1 - lam_0 ) * recfac[0][l];
		
		// do dynamic rescaling
		if( abs(lam_2) > ovflow )
		{
			lam_1 = lam_1*unflow;
			lam_2 = lam_2*unflow;
			scalel = scalel + 1;
			index_rstab = 20 + MAX( scalel, RSMIN );
    		corfac = rescale_tab[index_rstab];
		}
		else
		{
			if( (abs(lam_2) < unflow) && (abs(lam_2) != 0.0) )
			{
				lam_1 = lam_1*unflow;
				lam_2 = lam_2*unflow;
				scalel = scalel - 1;
				index_rstab = 20 + MAX( scalel, RSMIN );
    			corfac = rescale_tab[index_rstab];
			}
		}
    }
}

template void do_lam_lm_pol( long lmax, long m, double cth, double sth, double mfac, arr2<float> recfac, arr<float> lam_fact, arr2<float> &lam_lm );
template void do_lam_lm_pol( long lmax, long m, double cth, double sth, double mfac, arr2<double> recfac, arr<double> lam_fact, arr2<double> &lam_lm );

template<typename T> void ring_synthesis( long nside, long lmax, long mmax, arr<xcomplex<T> > datain, long nph, arr<double> &dataout, long kphi0 )
{
	long iw, ksign, m, k, kshift;
	double arg;
	xcomplex<T> dat, cplx;
	
	ksign = 1;
	kshift = pow( (double)-1, (double)kphi0 );// either 1 or -1
		
	arr<xcomplex<T> > bw;
	bw.alloc( nph );
	bw.fill( 0.0 );
	
	// all frequencies [-m,m] are wrapped in [0,nph-1]
	bw[0] = datain[0];
	for( m=1; m <= mmax; m++ )
	{
		iw = m % nph;// between 0 and nph-1  = m
		k = ( m - iw)  / nph;// number of 'turns'
		bw[iw] = bw[iw] + datain[m] * pow( (double)kshift, (double)k );// complex number
		iw = -m % nph;// between 0 and nph-1  = m
		k  = (-m - iw) / nph;// number of 'turns'
		bw[iw] = bw[iw] + conj(datain[m]) * xcomplex<REAL> (pow( (double)kshift, (double)k ), 0.);// complex number 
	}
	cout << "ring_synthesis loop for( m=1; m <= mmax; m++ ) done" << '\n';
	//     kshift**k = 1       for even turn numbers
    //               = 1 or -1 for odd  turn numbers : results from the shift in space

    //     applies the shift in position <-> phase factor in Fourier space
    dataout[0] = bw[0].real();
    for( iw=1; iw < nph/2; iw++ )
    {
    	m = ksign * iw;
    	if( kphi0 == 1 )
    	{
    		arg = m * PI / (double)(nph);
    		cplx.Set( cos(arg), sin(arg) );
    		dat = bw[iw] * cplx;
    	}
    	else
    	{
    		dat = bw[iw];
    	}
    	dataout[iw*2-1] = dat.real();
		dataout[iw*2] = dat.imag();
    }
    cout << "ring_synthesis loop for( iw=1; iw < nph/2; iw++ ) done" << '\n';
    
    iw = nph / 2;
	m = ksign * iw;
	if( kphi0 == 1 )
    {
    	arg = m * PI / (double)(nph);
    	cplx.Set( cos(arg), sin(arg) );
    	dat = bw[iw] * cplx;
    }
    else
    {
    	dat = bw[iw];
    }
    dataout[iw*2-1] = dat.real();
    
    double minv, maxv;
    dataout.minmax( minv, maxv);
    
    cout << "dataout min: " << minv << " dataout max: " << maxv << '\n';
    
    //long size = dataout.size();
    //cout << "dataout size: " << size << '\n';
    //for( long ind=0; ind < size; ind++ )
    //{
    //	cout << " index: " << ind;// << " data: " << dataout[ind];
    //}
    
    rfft plan;
    plan.Set( nph );
    plan.backward_fftpack( dataout );
    
    cout << "\nring_synthesis fft done, exiting" << '\n';
}

template void ring_synthesis( long nside, long lmax, long mmax, arr<xcomplex<float> > datain, long nph, arr<double> &dataout, long kphi0 );
template void ring_synthesis( long nside, long lmax, long mmax, arr<xcomplex<double> > datain, long nph, arr<double> &dataout, long kphi0 );

/*
template<typename T> void alm2map_pol_der1( const Alm<xcomplex<T> > Alm_T, const Alm<xcomplex<T> > Alm_E, const Alm<xcomplex<T> > Alm_B, 
											Healpix_Map<T> &mapT, Healpix_Map<T> &mapQ, Healpix_Map<T> &mapU, 
											Healpix_Map<T> &mapT_d_theta, Healpix_Map<T> &mapQ_d_theta, Healpix_Map<T> &mapU_d_theta, 
											Healpix_Map<T> &mapT_d_phi, Healpix_Map<T> &mapQ_d_phi, Healpix_Map<T> &mapU_d_phi )
{
	planck_assert( mapT.Scheme() == RING, "alm2map_pol_der1: maps must be in RING scheme" );
	planck_assert( mapT.conformable( mapQ ) && mapT.conformable( mapU ), "alm2map_pol_der1: maps are not conformable" );
	
	planck_assert( mapT.conformable( mapT_d_theta ) && mapT.conformable( mapT_d_phi ), "alm2map_pol_der1: maps are not conformable" );
	
	planck_assert( mapT_d_theta.Scheme() == RING, "alm2map_pol_der1: maps must be in RING scheme" );
	planck_assert( mapT_d_theta.conformable( mapQ_d_theta ) && mapT_d_theta.conformable( mapU_d_theta ), "alm2map_pol_der1: maps are not conformable" );
	
	planck_assert( mapT_d_phi.Scheme() == RING, "alm2map_pol_der1: maps must be in RING scheme" );
	planck_assert( mapT_d_phi.conformable( mapQ_d_phi ) && mapT_d_phi.conformable( mapU_d_phi ), "alm2map_pol_der1: maps are not conformable" );
	
	planck_assert( Alm_T.conformable( Alm_E ) && Alm_T.conformable( Alm_B ), "alm2map_pol_der1: a_lm are not conformable" );
	
	long startpix[SMAXCHK];
	long nph[SMAXCHK];
	long kphi0[SMAXCHK];
	double cth[SMAXCHK];
	double sth[SMAXCHK];
	double one_on_s[SMAXCHK];
	
	long nside = mapT.Nside();
	long npix = mapT.Npix();
	long lmax = Alm_T.Lmax();
	long mmax = Alm_T.Mmax();
	
	long nrings = 2*nside;
	long nphmx = 4*nside;
	
	long nchunks = nrings/SMAXCHK + 1;// number of chunks
	long chunksize = ( nrings + nchunks - 1 )/nchunks;// <= SMAXCHK
	
	arr<T> mfac, normal_l, lam_fact, lam_fact_der;
	arr<double> ring;
	mfac.alloc( mmax+1 );
	normal_l.alloc( lmax+1 );
	
	arr3<xcomplex<T> > b_d1;
	b_d1.alloc( 18, mmax+1, chunksize );
	
	arr2<T> recfac, dalm, lam_lm;
		
	arr<xcomplex<T> > bline;
	
	// warning compilation case OpenMP
	//
//	if( do_openmp() == false )// if (.not. do_openmp())
//	{
//		lam_fact.alloc( lmax+1 );
//		lam_fact_der.alloc( lmax+1 );
//		ring.alloc( nphmx );
//		recfac.alloc( 2, lmax+1 );
//		dalm.alloc( 6, lmax+1 );
//		lam_lm.alloc( 3, lmax+1 );
//		bline.alloc( mmax+1 );
//	}
	
	
	gen_mfac( mmax, mfac );// init mfac array
	
	init_rescale();
	
	gen_normpol( lmax, normal_l );// generate Polarization normalisation
	
	mapT.fill(0.0);
	mapQ.fill(0.0);
	mapU.fill(0.0);
	
	mapT_d_theta.fill(0.0);
	mapQ_d_theta.fill(0.0);
	mapU_d_theta.fill(0.0);
	
	mapT_d_phi.fill(0.0);
	mapQ_d_phi.fill(0.0);
	mapU_d_phi.fill(0.0);
	
	cout << "alm2map_pol_der1 init done" << '\n';
	
	long lchk, uchk, ith, ithl;
	
	for( long ichunk=0; ichunk < nchunks; nchunks++ )
	{
		lchk = ichunk * chunksize + 1;
		uchk = MIN( lchk + chunksize - 1, nrings );
		
		for( ith = lchk; ith <= uchk; ith++ )
		{
			ithl = ith - lchk;// local index
			//get pixel location information
			get_pixel_layout( nside, ith, cth[ithl], sth[ithl], nph[ithl], startpix[ithl], kphi0[ithl] );
			one_on_s[ithl]  = 1.0 / sth[ithl];
		}
		
		//for each theta, and each m, computes
		//b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)

		//lambda_mm tends to go down when m increases (risk of underflow)
		//lambda_lm tends to go up   when l increases (risk of overflow)

		b_d1.fill(0);//pad with zeros
		
		long m, ll, l_min, k, k0, k1, par_lm, l;
		double fm, f2m, fm2, lam_lm1m, cth_ring, one_on_s1, one_on_s2, cotanth, fllp1, fl, a0, xp, at, aq, derW, derX, derY, f2, f3, b0t, b0p, bx;
		double factor[2];
		
		arr<T> b_ns, b_ns_p, b_ns_t;
		//b_ns.alloc(12);// index = -3:8
		//b_ns_p.alloc(12);
		//b_ns_t.alloc(12);
		
		cout << "alm2map_pol_der1 enter first parallel loop" << '\n';

		#pragma omp parallel default(none) \
					shared( lmax, mmax, lchk, uchk, rescale_tab, normal_l, cth, sth, mfac, Alm_T, Alm_E, Alm_B, one_on_s, b_d1 ) \
					private( recfac, dalm, lam_fact, lam_fact_der, m, ll, fm, f2m, fm2, ithl, l_min, k, k0, k1, par_lm, lam_lm, \
					lam_lm1m, cth_ring, one_on_s1, one_on_s2, cotanth, factor, b_ns, b_ns_t, b_ns_p, l, fllp1, fl, a0, xp, at, aq, \
					derW, derX, derY, f2, f3, b0t, b0p, bx )
		{
			// warning compilation case OpenMP
			//if( do_openmp() == true )// if ( do_openmp())
			//{
				recfac.alloc( 2, lmax+1 );
				dalm.alloc( 6, lmax+1 );
				lam_fact.alloc( lmax+1 );
				lam_fact_der.alloc( lmax+1 );          
				lam_lm.alloc( 3, lmax+1 );
			//}
			
			b_ns.alloc(12);// index = -3:8
			b_ns_p.alloc(12);
			b_ns_t.alloc(12);
			
			//printf( "alm2map_pol_der1 1st parallel region dynamic allocation done\n" );
			
			//printf( "alm2map_pol_der1 b_ns size: %ld\n", b_ns.size() );
			
			//printf( "alm2map_pol_der1 b_d1 size1: %ld, b_d1 size2: %ld, b_d1 size3: %ld\n", b_d1.size1(), b_d1.size2(), b_d1.size3() );
			
			#pragma omp for schedule( dynamic, 1 )
			for( m=0; m <= mmax; m++ )
			{
				//generate recursion factors (recfac) for Ylm of degree m
          		gen_recfac( lmax, m, recfac );
          		//generate Ylm relation factor for degree m
          		gen_lamfac_der( lmax, m, lam_fact_der );
          		gen_lamfac( lmax, m, lam_fact );
          		f2m = 2.0 * m;
          		fm2 = m*m;
          		fm  = m;
          		
          		//printf( "alm2map_pol_der1 recfac and Ylm generated\n" );
          		
          		//extract needed alm under memory and CPU efficient form
          		for( ll=m; ll <= lmax; ll++ )
          		{
          			dalm[0][ll] = Alm_T( ll, m ).re;
          			dalm[1][ll] = Alm_T( ll, m ).im;
          			dalm[2][ll] = Alm_E( ll, m ).re*normal_l[ll];
          			dalm[3][ll] = Alm_E( ll, m ).im*normal_l[ll];
          			dalm[4][ll] = Alm_B( ll, m ).re*normal_l[ll];
          			dalm[5][ll] = Alm_B( ll, m ).im*normal_l[ll];
          		}
          		
          		//printf( "alm2map_pol_der1 alm extracted\n" );
          		
          		for( ithl=0; ithl <= (uchk-lchk); ithl ++ )
          		{
          			l_min = l_min_ylm( m, sth[ithl] );
          			
          			if( lmax >= l_min )
          			{
          				//printf( "alm2map_pol_der1 enter if( lmax >= l_min )\n" );
          				
          				//compute lam_lm(p,theta) for all l>=m
          				do_lam_lm_pol( lmax, m, cth[ithl], sth[ithl], mfac[m], recfac, lam_fact, lam_lm );

                		cth_ring = cth[ithl];
                		one_on_s1 = one_on_s[ithl];
                		one_on_s2 = one_on_s[ithl]*one_on_s[ithl];
                		cotanth = cth_ring * one_on_s[ithl];

                		b_ns.fill(0.0);
                		b_ns_t.fill(0.0);
                		b_ns_p.fill(0.0);
                		
                		for( l=l_min; l <= lmax; l++ )
                		{
                			//printf( "alm2map_pol_der1 enter for( l=l_min; l <= lmax; l++ )\n" );
                			
                			fl = l;
                			fllp1 = l*l + l;
                			par_lm = 3;//! = (-1)^(l+m)

                			if( ( (l + m) % 2 ) == 1 )
                			{
                				par_lm = -par_lm;
                			}
                			
                			//printf( "alm2map_pol_der1 par_lm value: %ld\n", par_lm );
                			
                			//printf( "alm2map_pol_der1 index l value: %ld, index m value: %ld\n", l, m );
                			
                			//--------------------------
							// f = Y_lm * a_lm;
							factor[0] = lam_lm[0][l] * dalm[0][l];
							factor[1] = lam_lm[0][l] * dalm[1][l];
							
							//printf( "alm2map_pol_der1 factor param set\n" );
							
							//printf( "alm2map_pol_der1 b_ns size: %ld\n", b_ns.size() );
							
							b_ns[par_lm+3] = b_ns[par_lm+3] + factor[0];// T even
							b_ns[par_lm+4] = b_ns[par_lm+4] + factor[1];
							
							//printf( "alm2map_pol_der1 b_ns tab T even part done\n" );
							
							b_ns[par_lm+5] = b_ns[par_lm+5] - lam_lm[1][l] * dalm[2][l];// Q, U  even
							b_ns[par_lm+6] = b_ns[par_lm+6] - lam_lm[1][l] * dalm[3][l];
							b_ns[par_lm+7] = b_ns[par_lm+7] - lam_lm[1][l] * dalm[4][l];
							b_ns[par_lm+8] = b_ns[par_lm+8] - lam_lm[1][l] * dalm[5][l];
							
							//printf( "alm2map_pol_der1 b_ns tab Q, U even part done\n" );
							
							b_ns[5-par_lm] = b_ns[5-par_lm] + lam_lm[2][l] * dalm[5][l];// Q odd
							b_ns[6-par_lm] = b_ns[6-par_lm] - lam_lm[2][l] * dalm[4][l];
							b_ns[7-par_lm] = b_ns[7-par_lm] - lam_lm[2][l] * dalm[3][l];// U odd
							b_ns[8-par_lm] = b_ns[8-par_lm] + lam_lm[2][l] * dalm[2][l];
							
							//printf( "alm2map_pol_der1 b_ns tab computed\n" );

							//-------------------------- 1st derivatives
							//printf( "alm2map_pol_der1 compute 1st derivatives\n" );
							if( l > 0 )
							{
								// df/dphi = i * m * Y_lm * a_lm
								f2 = m * lam_lm[1][l];
                      			f3 = m * lam_lm[2][l];
                      			b_ns_p[par_lm+3] = b_ns_p[par_lm+3] - m * factor[1];// warning negative index
                      			b_ns_p[par_lm+4] = b_ns_p[par_lm+4] + m * factor[0];
                      			b_ns_p[par_lm+5] = b_ns_p[par_lm+5] + f2 * dalm[3][l];
                      			b_ns_p[par_lm+6] = b_ns_p[par_lm+6] - f2 * dalm[2][l];
                      			b_ns_p[par_lm+7] = b_ns_p[par_lm+7] + f2 * dalm[5][l];
                      			b_ns_p[par_lm+8] = b_ns_p[par_lm+8] - f2 * dalm[4][l];
                      			
                     			b_ns_p[5-par_lm] = b_ns_p[5-par_lm] + f3 * dalm[4][l];// Q odd
                      			b_ns_p[6-par_lm] = b_ns_p[6-par_lm] + f3 * dalm[5][l];
                      			b_ns_p[7-par_lm] = b_ns_p[7-par_lm] - f3 * dalm[2][l];// U odd
                      			b_ns_p[8-par_lm] = b_ns_p[8-par_lm] - f3 * dalm[3][l];
                      			
                      			// dY_lm/dtheta = (l/tan(theta)*Y_lm                         -fact/sin(theta)*Y_l-1m)
                      			// dW_lm/dtheta = (l/tan(theta)*W_lm - S*m/l/sin(theta)*X_lm -fact/sin(theta)*sqrt(1-S^2/l^2)*W_l-1m
                      			// dX_lm/dtheta = (l/tan(theta)*X_lm - S*m/l/sin(theta)*W_lm -fact/sin(theta)*sqrt(1-S^2/l^2)*X_l-1m
                      			
                      			a0 = fl * cotanth;// l/tan(theta)
                      			at = lam_fact_der[l] * one_on_s1;// sqrt((2l+1)/(2l-1)*(l^2-m^2))/sin(theta)
                      			derY = a0 * lam_lm[0][l] - at * lam_lm[0][l-1];
                      			b_ns_t[3-par_lm] = b_ns_t[3-par_lm] + derY * dalm[0][l];// T odd
                      			b_ns_t[4-par_lm] = b_ns_t[4-par_lm] + derY * dalm[1][l];
							}
							
							if( l > 1 )
							{
								xp = (2*m)*one_on_s1/fl;// spin m / (l sin(theta))
								aq = at * sqrt( 1.0 - 4.0/( fl * fl ) );// at * sqrt(l^2-spin^2)/l
								derW = a0 * lam_lm[1][l] - aq * lam_lm[1][l-1] + xp * lam_lm[2][l];
                      			derX = a0 * lam_lm[2][l] - aq * lam_lm[2][l-1] + xp * lam_lm[1][l];
                      			b_ns_t[5-par_lm] = b_ns_t[5-par_lm] - derW * dalm[2][l];// Q, U  odd
                      			b_ns_t[6-par_lm] = b_ns_t[6-par_lm] - derW * dalm[3][l];
                      			b_ns_t[7-par_lm] = b_ns_t[7-par_lm] - derW * dalm[4][l];
                      			b_ns_t[8-par_lm] = b_ns_t[8-par_lm] - derW * dalm[5][l];
                      			b_ns_t[5+par_lm] = b_ns_t[5+par_lm] + derX * dalm[5][l];// Q even
                      			b_ns_t[6+par_lm] = b_ns_t[6+par_lm] - derX * dalm[4][l];
                      			b_ns_t[7+par_lm] = b_ns_t[7+par_lm] - derX * dalm[3][l];// U even
                      			b_ns_t[8+par_lm] = b_ns_t[8+par_lm] + derX * dalm[2][l];
							}
							
                		}//end loop on l
                		
                		for( k=0; k <= 2; k++ )
                		{
                			// loop on T,Q,U
                			k0 = 2*k;
                   			k1 = k0 + 1;
                   			
                   			// fields
                   			b_d1( 0+k, m, ithl).re = b_ns[k0+6] + b_ns[k0];// north=Even+Odd
                   			b_d1( 0+k, m, ithl).im = b_ns[k1+6] + b_ns[k1];
                   			b_d1( 3+k, m, ithl).re = b_ns[k0+6] - b_ns[k0];// south=Even-Odd
                   			b_d1( 3+k, m, ithl).im = b_ns[k1+6] - b_ns[k1];
                   			
                   			//printf( "alm2map_pol_der1 b_d1: T, Q, U\n" );
                   			
                   			// dfield/dtheta
                   			b_d1( 6+k, m, ithl).re = b_ns_t[k0+6] + b_ns_t[k0];// north=Even+Odd
                   			b_d1( 6+k, m, ithl).im = b_ns_t[k1+6] + b_ns_t[k1];
                   			b_d1( 9+k, m, ithl).re = b_ns_t[k0+6] - b_ns_t[k0];// south=Even-Odd
                   			b_d1( 9+k, m, ithl).im = b_ns_t[k1+6] - b_ns_t[k1];
                   			
                   			//printf( "alm2map_pol_der1 b_d1: d/dtheta T, Q, U\n" );
                   			
                   			// dfield/dphi/sin(theta)
                   			b_d1( 12+k, m, ithl).re = ( b_ns_p[k0+6] + b_ns_p[k0] ) * one_on_s1;
                   			b_d1( 12+k, m, ithl).im = ( b_ns_p[k1+6] + b_ns_p[k1] ) * one_on_s1;
                   			b_d1( 15+k, m, ithl).re = ( b_ns_p[k0+6] - b_ns_p[k0] ) * one_on_s1;
                   			b_d1( 15+k, m, ithl).im = ( b_ns_p[k1+6] - b_ns_p[k1] ) * one_on_s1;
                   			
                   			//printf( "alm2map_pol_der1 b_d1: d/dphi T, Q, U\n" );
                		}
                		                		
          			}// end if( lmax >= l_min )
          			
          		}// and loop on ithl
          			
          	}// end loop on m
          	
          	// warning compilation case OpenMP
          	//if( do_openmp() == true )// if ( do_openmp())
          	//{
          		//recfac.dealloc();
          		//dalm.dealloc();
          		//lam_fact.dealloc();
          		//lam_fact_der.dealloc();
          		//lam_lm.dealloc();
          	//}
          	//b_ns.dealloc();
			//b_ns_p.dealloc();
			//b_ns_t.dealloc();
          	
		}// end of parallel region
		
		cout << "alm2map_pol_der1 first parallel loop done" << '\n';
		
		long nphl, istart_south, istart_north;
		
		//shared mapT, mapQ, mapU, mapT_d_theta, mapQ_d_theta, mapU_d_theta, mapT_d_phi, mapQ_d_phi, mapU_d_phi
		//#pragma omp parallel default(none) \
		//			shared( nside, lmax, mmax, npix, nrings, nphmx, lchk, uchk, b_d1, nph, startpix, kphi0 ) \
		//			private( ithl, nphl, istart_north, istart_south, ith, ring, bline, k0 )
		//{
			// warning compilation case OpenMP
			//if( do_openmp() == true )// if ( do_openmp())
			//{
				ring.alloc( nphmx );
				bline.alloc( mmax );
			//}
			
			printf( "alm2map_pol_der1 2nd parallel region dynamic allocation done\n" );
			
			printf( "alm2map_pol_der1 b_d1 size1: %ld, b_d1 size2: %ld, b_d1 size3: %ld\n", b_d1.size1(), b_d1.size2(), b_d1.size3() );
			
			//#pragma omp for schedule( dynamic, 1 )
			for( ithl=0; ithl <= (uchk-lchk); ithl ++ )
			{
				printf( "alm2map_pol_der1 enter for( ithl=0; ithl <= (uchk-lchk); ithl ++ )\n" );
				nphl = nph[ithl];
				istart_north = startpix[ithl];
          		istart_south = npix - istart_north - nphl;
          		ith = ithl + lchk;
          		
          		//        ---------------------------------------------------------------
				//        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta)
				//        ---------------------------------------------------------------
				
				long ind;
				printf( "alm2map_pol_der1 enter 1st loop for( k0=0; k0 <= 2; k0++ )\n" );
				for( k0=0; k0 <= 2; k0++ )
				{
					cout << "k0: " << k0 << " mmax: " << mmax << " ithl: " << ithl << '\n';
					bline.fill( 0.0 );
					for( ind=0; ind < mmax; ind++ )
					{
						bline[ind].re = b_d1( k0, ind, ithl ).re;
						bline[ind].im = b_d1( k0, ind, ithl ).im;
					}
					printf( "\nalm2map_pol_der1 bline copy from b_d1\n" );
					cout << "nside: " << nside << " lmax: " << lmax<< " mmax: " << mmax << " nphl: " << nphl << " ithl: " << ithl << " kphi0[ithl]: " << kphi0[ithl] << '\n';
					ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // north hemisph. + equator
					printf( "alm2map_pol_der1 ring_synthesis T, Q, U north hemisph. + equator done\n" );
					for( ind=0; ind < nphl; ind++ )
					{
						switch( k0 )
						{
							case 0 :
								mapT[istart_north+ind] = ring[ind];
								break;
							case 1 :
								mapQ[istart_north+ind] = ring[ind];
								break;
							case 2:
								mapU[istart_north+ind] = ring[ind];
								break;
						}
					}
				}
				
				printf( "alm2map_pol_der1 enter 2nd loop for( k0=0; k0 <= 2; k0++ )\n" );
				for( k0=0; k0 <= 2; k0++ )
				{
					for( ind=0; ind < mmax; ind++ )
					{
						bline[ind] = b_d1( 6+k0, ind, ithl );
					}
					ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // north hemisph. + equator
					printf( "alm2map_pol_der1 ring_synthesis T, Q, U d/d_theta north hemisph. + equator done\n" );
					for( ind=0; ind < nphl; ind++ )
					{
						switch( k0 )
						{
							case 0 :
								mapT_d_theta[istart_north+ind] = ring[ind];
								break;
							case 1 :
								mapQ_d_theta[istart_north+ind] = ring[ind];
								break;
							case 2:
								mapU_d_theta[istart_north+ind] = ring[ind];
								break;
						}
					}
				}
				
				printf( "alm2map_pol_der1 enter 3rd loop for( k0=0; k0 <= 2; k0++ )\n" );				
				for( k0=0; k0 <= 2; k0++ )
				{
					for( ind=0; ind < mmax; ind++ )
					{
						bline[ind] = b_d1( 12+k0, ind, ithl );
					}
					ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // north hemisph. + equator
					printf( "alm2map_pol_der1 ring_synthesis T, Q, U d/d_phi north hemisph. + equator done\n" );
					for( ind=0; ind < nphl; ind++ )
					{
						switch( k0 )
						{
							case 0 :
								mapT_d_phi[istart_north+ind] = ring[ind];
								break;
							case 1 :
								mapQ_d_phi[istart_north+ind] = ring[ind];
								break;
							case 2:
								mapU_d_phi[istart_north+ind] = ring[ind];
								break;
						}
					}
				}
				
				if( ith < nrings )
				{
					printf( "alm2map_pol_der1 enter 4th loop for( k0=0; k0 <= 2; k0++ )\n" );
					for( k0=0; k0 <= 2; k0++ )
					{
						for( ind=0; ind < mmax; ind++ )
						{
							bline[ind] = b_d1( 3+k0, ind, ithl );
						}
						ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // south hemisph. w/o equat
						printf( "alm2map_pol_der1 ring_synthesis T, Q, U south hemisph. w/o equat done\n" );
						for( ind=0; ind < nphl; ind++ )
						{
							switch( k0 )
							{
								case 0 :
									mapT[istart_south+ind] = ring[ind];
									break;
								case 1 :
									mapQ[istart_south+ind] = ring[ind];
									break;
								case 2:
									mapU[istart_south+ind] = ring[ind];
									break;
							}
						}
					}

					printf( "alm2map_pol_der1 enter 5th loop for( k0=0; k0 <= 2; k0++ )\n" );					
					for( k0=0; k0 <= 2; k0++ )
					{
						for( ind=0; ind < mmax; ind++ )
						{
							bline[ind] = b_d1( 9+k0, ind, ithl );
						}
						ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // south hemisph. w/o equat
						printf( "alm2map_pol_der1 ring_synthesis T, Q, U d/d_theta south hemisph. w/o equat done\n" );
						for( ind=0; ind < nphl; ind++ )
						{
							switch( k0 )
							{
								case 0 :
									mapT_d_theta[istart_south+ind] = ring[ind];
									break;
								case 1 :
									mapQ_d_theta[istart_south+ind] = ring[ind];
									break;
								case 2:
									mapU_d_theta[istart_south+ind] = ring[ind];
									break;
							}
						}
					}

					printf( "alm2map_pol_der1 enter 6th loop for( k0=0; k0 <= 2; k0++ )\n" );				
					for( k0=0; k0 <= 2; k0++ )
					{
						for( ind=0; ind < mmax; ind++ )
						{
							bline[ind] = b_d1( 15+k0, ind, ithl );
						}
						ring_synthesis( nside, lmax, mmax, bline, nphl, ring, kphi0[ithl] ); // south hemisph. w/o equat
						printf( "alm2map_pol_der1 ring_synthesis T, Q, U d/d_phi south hemisph. w/o equat done\n" );
						for( ind=0; ind < nphl; ind++ )
						{
							switch( k0 )
							{
								case 0 :
									mapT_d_phi[istart_south+ind] = ring[ind];
									break;
								case 1 :
									mapQ_d_phi[istart_south+ind] = ring[ind];
									break;
								case 2:
									mapU_d_phi[istart_south+ind] = ring[ind];
									break;
							}
						}
					}					
				}// end if( ith < nrings )
								
			}// and loop on ithl
			
			// warning compilation case OpenMP
          	//if( do_openmp() == true )// if ( do_openmp())
          	//{
          		//ring.dealloc();
          		//bline.dealloc();
          	//}
						
		//}// end of parallel region
		
		cout << "alm2map_pol_der1 second parallel loop done" << '\n';
					
	}// end loop on ichunck

	//     --------------------
    //     free memory and exit
    //     --------------------
    
    // warning compilation case OpenMP
    
//	if( do_openmp() == false )// if (.not. do_openmp())
//	{
//		lam_fact.dealloc();
//		lam_fact_der.dealloc();
//		ring.dealloc();
//		recfac.dealloc();
//		dalm.dealloc();
//		lam_lm.dealloc();
//		bline.dealloc();
//	}
	
    
    //mfac.dealloc();
    //b_d1.dealloc();
    //normal_l.dealloc();
}

template void alm2map_pol_der1( const Alm<xcomplex<float> > Alm_T, const Alm<xcomplex<float> > Alm_E, const Alm<xcomplex<float> > Alm_B, 
											Healpix_Map<float> &mapT, Healpix_Map<float> &mapQ, Healpix_Map<float> &mapU, 
											Healpix_Map<float> &mapT_d_theta, Healpix_Map<float> &mapQ_d_theta, Healpix_Map<float> &mapU_d_theta, 
											Healpix_Map<float> &mapT_d_phi, Healpix_Map<float> &mapQ_d_phi, Healpix_Map<float> &mapU_d_phi );
											
template void alm2map_pol_der1( const Alm<xcomplex<double> > Alm_T, const Alm<xcomplex<double> > Alm_E, const Alm<xcomplex<double> > Alm_B, 
											Healpix_Map<double> &mapT, Healpix_Map<double> &mapQ, Healpix_Map<double> &mapU, 
											Healpix_Map<double> &mapT_d_theta, Healpix_Map<double> &mapQ_d_theta, Healpix_Map<double> &mapU_d_theta, 
											Healpix_Map<double> &mapT_d_phi, Healpix_Map<double> &mapQ_d_phi, Healpix_Map<double> &mapU_d_phi );
*/

#endif
