#include"HealpixClass.h"


double gammln(double xx)
{
       double x,y,tmp,ser;
       
       static double cof[6]={ 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
       
       int j;

       y=x=xx;
       
       tmp=x+5.5;
       
       tmp -= (x+0.5)*log(tmp);
       
       ser=1.000000000190015;
       
       for (j=0;j<=5;j++) ser += cof[j]/++y;
       
       return -tmp+log(2.5066282746310005*ser/x);
}


 arr<double> mrs_wigner3j2( double il1, double il2, arr<double> il3, arr<double> lnwpro )
{

	long l1 = (long)il1;
	long l2 = (long)il2;
	
	long n_elements = il3.size();

	arr<long> l3( n_elements );
	arr<long> L( n_elements );
	arr<long> L_2( n_elements );
	
	arr<double> c( n_elements);
	c.fill(0.0);
	
	long min = abs( l1 - l2 );
	long max = l1 + l2;
	long k;
	long index1, index2, index3, indexl;
	double lnw1, lnw2, lnw3, lnwl, lnc;
	
	for( k = 0; k < n_elements; k++ )
	{
		l3[k] = (long)il3[k];
		
		L[k] = l3[k] + l1 + l2;
		
		L_2[k] = L[k] / 2;
		
		if( ( (2*L_2[k]-L[k]) == 0 )&&( l3[k] >= min )&&( l3[k] <= max ) )
		{
			index1 = L_2[k]-l1;
			lnw1 = lnwpro[index1];
		
			index2 = L_2[k]-l2;
			lnw2 = lnwpro[index2];
			
			index3 = L_2[k]-l3[k];
			lnw3 = lnwpro[index3];
			
			indexl = L_2[k];
			lnwl = lnwpro[indexl];
			
			lnc = lnw1 + lnw2 + lnw3 - lnwl - log( L[k] + 1.0 );
			c[k] = exp( lnc );
		}
	}
	return c;	
}



//----------------------------------------------------------------------------

dblarray mrs_mmake_mll( arr<double> ell, arr<double> well, long lmax )
{
	long index;
	long n_elements_ell = ell.size();
	
	arr<long> k( 2+3*n_elements_ell);
	k.fill(0);
	long k_size = k.size();
	arr<double> lnwpro( k_size);
	lnwpro.fill(0.0);
	#pragma omp parallel default(none) shared( k, k_size, lnwpro ) private( index )
	{	
		#pragma omp for schedule( dynamic, 1 )
		for( index = 0; index < k_size; index++ )
		{
			k[index] = index;
		
			lnwpro[index] = gammln( (double)2.0*index + 1.0 ) - (double)2.0*gammln( (double)index + 1.0 );
		}
	}

	dblarray m;
	m.alloc( lmax+1, lmax+1 );
	int l1;
	int l2;
	double a;
	arr<double> mwigner;
	arr<double> c( n_elements_ell );
	dblarray b;
		
	#pragma omp parallel default(none) shared( ell, well, n_elements_ell, c ) private( index )
	{	
		#pragma omp for schedule( dynamic, 1 )
		for( index = 0; index < n_elements_ell; index++ )
		{
			c[index] = (2.0*ell[index] + 1.0)*(well[index]);
		}
	}
	
	#pragma omp parallel default(none) shared( ell, n_elements_ell, c, m, lmax, lnwpro ) private( index, l1, l2, a, b, mwigner )
	{	
		#pragma omp for schedule( dynamic, 1 )
		for( l1 = 0; l1 <= lmax; l1++ )
		{
			for( l2 = 0; l2 <= lmax; l2++ )
			{
				mwigner.alloc( n_elements_ell );
				b.alloc( n_elements_ell );

				mwigner = mrs_wigner3j2( ell[l1], ell[l2], ell, lnwpro );
			
				for( index = 0; index < n_elements_ell; index++ )
				{
					b(index) = (c[index])*(mwigner[index]);
				}
			
				a = b.total();
			
				m( l1, l2 ) = ( 2.0*ell[l2] + 1.0 )*a/( 4.0*PI );
			}
		}
	}

	return m;	
}



//----------------------------------------------------------------------------

dblarray mrs_matmask( Hmap<double> Mask, int & lmax )
{

	long Nside = Mask.Nside();
 
	CAlmR ALM;
	//ALM.Norm = false;
	ALM.alloc( Nside, lmax );
    lmax = ALM.Lmax();
	ALM.alm_trans( Mask );
	
	PowSpec SigPowSpec( 1, lmax );
	ALM.alm2powspec( SigPowSpec );
	
	// cout << "mrs_matmask" << lmax << endl;
    arr<double> ell( lmax + 1 );
	long index;
	
	//	pragma omp parallel default(none) shared( ell ) private( index )
	{	
	//	pragma omp for schedule( dynamic, 1 )
		for( index = 0; index <= lmax; index++ )
		{
			ell[index] = (double)index;
		}
	}
	
	arr<double> pm( SigPowSpec.tt() );

	dblarray MatMask = mrs_mmake_mll( ell, pm, lmax );

	return MatMask;
}

//--------------------------------------------------------------------------

void matmask_mult_cl(PowSpec & PS, dblarray & Mat, PowSpec & PSout, bool transpose=false)
{
  if (transpose == false)
  {
     
     for (int l=0; l <= PS.Lmax(); l++)
     {
        PSout.tt(l) = 0;
        for (int l1=0; l1 <= PS.Lmax(); l1++)
               PSout.tt(l) += Mat(l,l1) * PS.tt(l1);
     }
  }
  else {
      for (int l=0; l <= PS.Lmax(); l++)
      {
          PSout.tt(l) = 0;
          for (int l1=0; l1 <= PS.Lmax(); l1++)
              PSout.tt(l) += Mat(l1,l) * PS.tt(l1);
      }
  }
}

//--------------------------------------------------------------------------



void matmask_mult_tabcl(dblarray & PS, dblarray & Mat, dblarray & PSout, bool transpose=false)
{
    if (transpose == false)
    {
        for (int l=0; l < PS.nx(); l++)
        {
            PSout(l) = 0;
            for (int l1=0; l1 < PS.nx(); l1++)
                PSout(l) += Mat(l,l1) * PS(l1);
        }
    }
    else {
        for (int l=0; l < PS.nx(); l++)
        {
            PSout(l) = 0;
            for (int l1=0; l1 < PS.nx(); l1++)
                PSout(l) += Mat(l1,l) * PS(l1);
        }
    }
}

//--------------------------------------------------------------------------


/*********************************************************************/

void iter_master_decconv(PowSpec & ClData, PowSpec & ClSol, Hdmap & MaskMap, dblarray & MatMask, int Niter, int Lmax_Tot, bool Verbose)
{
    double Fsky = MaskMap.Npix();
    double Tot = 0.;
    for (int l=0; l < MaskMap.Npix(); l++) Tot += (double) MaskMap[l];
    Fsky /= Tot;
    if (Verbose == True) cout << "Fsky (%) = " << (double) 1. / Fsky * (double) 100. << endl; 
    
    
    // Compute the matrix related to the mask
    if (MatMask.naxis() == 2) 
    {
        if ((MatMask.nx() != Lmax_Tot+1) || (MatMask.ny() != Lmax_Tot+1))
        {
            cout << "Error: the input mask matrix has not the correct dimension: (" << MatMask.nx() << "," << MatMask.ny() << "), instead of " <<  Lmax_Tot << endl;
            exit(-1);
        }
    }
    else 
    {
        if (Verbose == True)  cout << " Compute the mask matrix: Lmax = "  << Lmax_Tot << endl;
        
        MatMask = mrs_matmask( MaskMap, Lmax_Tot);
        if (Verbose == True)  cout << " Computed Matrix: Nx = " << MatMask.nx() << ", Ny = " << MatMask.ny() << endl;
        // fits_write_dblarr("xx_matmask.fits", MatMask );
    }
    
    PowSpec  ClResi, ClConv;
    mrs_alloc_powspec(ClSol, Lmax_Tot);
    mrs_alloc_powspec(ClResi, Lmax_Tot);
    mrs_alloc_powspec(ClConv, Lmax_Tot);
    
    cout << Fsky << endl;
    
    for (int l=0; l <=  Lmax_Tot; l++) ClSol.tt(l) = ClData.tt(l) * Fsky;
    dblarray Z(Lmax_Tot+1);
    
    for (int i=0; i < Niter; i++)
    { 
        matmask_mult_cl(ClSol, MatMask, ClConv);
        // for (int l=0; l <=  Lmax_Tot; l++) Z(l) = ClResi.tt(l);
        // Z.info("ClConv");
        
        for (int l=0; l <=  Lmax_Tot; l++) ClResi.tt(l) = ClData.tt(l) - ClConv.tt(l);
        
        // for (int l=0; l <=  Lmax_Tot; l++) Z(l) = ClResi.tt(l);
        // Z.info("CLResi");
        
        matmask_mult_cl(ClResi, MatMask, ClConv, True);
        
        // if (Verbose == True) cout << "   ResiSpec = " << sqrt( ErrSpec / (float) ALM.Lmax() )  << endl;
        for (int l=0; l <=  Lmax_Tot; l++) ClSol.tt(l) = MAX(0, ClSol.tt(l) + ClConv.tt(l));
        
        if (Verbose == True)  
        {
            double ErrSpec=0.;
            for (int l=0; l <= Lmax_Tot; l++)  Z(l) = ClResi.tt(l);
            cout << "Iter  " << i+1 << ", Sigma_ResiSpec = " << Z.sigma() << endl;
        }        
    }
}

/*********************************************************************/

