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
**    File:  HealpixClass.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class interface to the Healpix package   
**    ----------- 
**                 
******************************************************************************/

#include "HealpixClass.h"

int  mrs_get_lmax (int  & Lmax, int Nside, float ZeroPadding)
{
    if (Lmax <= 0) Lmax = 3 * Nside;
    if (Lmax > ALM_MAX_L) Lmax = ALM_MAX_L;
    int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
    return Lmax_Tot;
}


void mrs_write_3maps(char *Name, Hmap<REAL> & map, Hmap<REAL> & mapdth, Hmap<REAL> & mapdph)
{
    fitshandle out;
    remove(fitsname(Name));
    out.create(fitsname(Name));
     // write_Healpix_map_to_fits (out,map,mapdth,mapdph, FITSUTIL<T>::DTYPE);
    write_Healpix_map_to_fits (out,map,mapdth,mapdph, planckType<REAL>());
}
 

void CAlmR::alloc(int Nside, int u_lmax, bool Fast) 
{
	AllocNside=Nside;
	FastALM=Fast;
	int L_Max=u_lmax;
	if (L_Max <=0) L_Max = 3*Nside;
	if (L_Max > ALM_MAX_L_TOT) L_Max = ALM_MAX_L_TOT; 
	int M_Max= L_Max; 
	FastALM=Fast;
	// double lm = sqrt( (double) (Nside)* (double) (Nside)*12.);
	// NormVal = sqrt( (double)(lm*(lm+1)/(4.* PI)));
	double Nelem = (double) (Nside)* (double) (Nside)*12.;
	NormVal = sqrt( Nelem /(4.* PI));

    char const* temp = getenv("HEALPIX");
    std::string DirWeight;
    if(temp != NULL)
    {
        DirWeight = std::string(temp) + std::string("/data");
        // cout < " DirWeight " << DirWeight << endl;
    }
    else
    {
        cout << "Error: HEALPIX environment variable is not defined." << endl;
        cout << "       cannot find directory $HEALPIX/data." << endl;
        exit(-1);
    }
    // cout << DirWeight << endl;
    
	Set(L_Max, M_Max);
	weight_T.alloc (2* Nside);
	BeamEff.alloc(L_Max+1);

	if (Fast == false)
	{
		// read the ring
    if (Verbose == true)  std::cout << "DirW = " << DirWeight << " NormVal = " << NormVal << std::endl;
		read_weight_ring ( DirWeight, Nside, weight_T);
		for (int m=0; m< (int) weight_T.size(); ++m) weight_T[m]+=1;
	}
	else weight_T.fill(1);
}

  // =============================
 
 void CAlmR::extract_median_powspec(PowSpec &powspec)
 {
   arr<double> tt(Lmax()+1);
   fltarray TabM(Lmax()+1);
   
   for (int l=0; l<= Lmax(); ++l)
   {
      TabM(0) = norm((*this)(l,0));
      int limit = std::min(l,(*this).Mmax());
      for (int m=1; m<=limit; ++m) TabM(m) = norm((*this)(l,m));
      
      tt[l] = get_median( TabM.buffer(), (limit+1));
   }
   powspec.Set(tt);
 }
 
 // =============================
 
 int  CAlmR::hard_threshold(float T, int & MaxNonZeroL)
 {
 	int Cpt=0;
 	MaxNonZeroL = 0;
     for (int l=1; l <= Lmax(); l++)
     for (int m=0; m <= l; m++)
     {
          REAL  Rea = real((*this)(l,m));
          REAL  Ima = imag((*this)(l,m));
 
     	 if  (ABS(Rea)  < T) Rea = 0.;
     	 else
     	 {
     	 	  Cpt++;
     	 	  if (MaxNonZeroL < l) MaxNonZeroL = l;
              // printf("%5.5f ", (*this)(l,m).re);
     	 }
     	 if  (ABS(Ima) < T) Ima = 0.;
     	 else
     	 {
     	 	  Cpt++;
     	 	  if (MaxNonZeroL < l) MaxNonZeroL = l;
             // printf("%5.5f ", (*this)(l,m).im);
     	 }
         (*this)(l,m) = xcomplex<REAL> (Rea, Ima);
     }
     // std::cout << "hard_threshold:Thres =  " << T << ", Keep " <<  Cpt <<  std::endl;
     return Cpt;
 }
 
 // =============================
 
 
 // =============================
 
 xcomplex<REAL> CAlmR::lm_soft_threshold(int l, int m, float T)
 {
 	xcomplex<REAL> ValRet = 0.;
 
    ValRet = (*this)(l,m);
    double M = sqrt(  norm( ValRet ) );
    if (M == 0)  ValRet = xcomplex<REAL> (0.,0.);
    else
    {
        if (T > M)
        {
           ValRet = xcomplex<REAL> (0.,0.);
        }
        else
        {
           double Coef = (1. - T / M);
           if (Coef <= 0)  Coef = 0.;
           ValRet  *= Coef;
        }
    }
    return ValRet;
 }
 
 // =============================
 
 
 int  CAlmR::soft_threshold(float T, int & MaxNonZeroL)
 {
  	int Cpt=0;
 	MaxNonZeroL = 0;
     for (int l=1; l <= Lmax(); l++)
     for (int m=0; m <= l; m++)
     {
         double M = norm((*this)(l,m));
         M = (M > 0) ? sqrt(M): 0.;
         if (M == 0)
     	 {
             (*this)(l,m) = xcomplex<REAL> (0., 0.);
      	 }
     	 else
     	 {
     	    if (T > M)
     	    {
                 (*this)(l,m) = xcomplex<REAL> (0., 0.);
     	    }
     	    else
     	    {
                 double Coef= (1. - T / M);
                 if (MaxNonZeroL < l) MaxNonZeroL = l;
      	        (*this)(l,m)  *= Coef;
     	    	Cpt++;
     	    }
     	}
     }
     return Cpt;
     // std::cout << " CptALM = " << Cpt << std::endl;
 }
 
 
 // =============================
 
 double CAlmR::max_absalm()
 {
     double Max=0.;
     for (int l=1; l <= Lmax(); l++)
     for (int m=0; m <= l; m++)
     {
         REAL  Rea = real((*this)(l,m));
         REAL  Ima = imag((*this)(l,m));
 
 
     	 if  (  ABS(Rea) > Max) Max = ABS(Rea);
     	 if  (  ABS(Ima) > Max) Max = ABS(Ima);
     }
 	return (Max);
 }
 // =============================
 
 void CAlmR::deriv(Hmap<REAL> & map,  Hmap<REAL> & mapdth,  Hmap<REAL> & mapdph,  float fwhm_arcmin, int RecNside)
 {
 	if (Norm == true)
 	{
 	    for (int l=0; l <= Lmax(); l++)
 	    for (int m=0; m <= l; m++)   operator()(l,m) *=  (1./NormVal);
 	}
 	
     if (fwhm_arcmin>0) smoothWithGauss (*this, fwhm_arcmin);
   
 	if  (RecNside != 0) map.SetNside(RecNside, RING);
 	else map.SetNside(AllocNside, RING);
 	if  (RecNside != 0) mapdth.SetNside(RecNside, RING);
 	else mapdth.SetNside(AllocNside, RING);
 	if  (RecNside != 0) mapdph.SetNside(RecNside, RING);
 	else mapdph.SetNside(AllocNside, RING);
  
     double offset = (*this)(0,0).real()/sqrt(fourpi);
 	(*this)(0,0) = 0;
     alm2map_der1(*this,map,mapdth,mapdph);
     for (int m=0; m<map.Npix(); ++m) map[m]+=offset;
     
  }
  
 // =============================
 
 void CAlmR::read(char *Name, int Nside, bool Fast)
 {
     int lmaxIn, dummy;
     string infile = string(fitsname(Name));
     // std::cout << "PB ALMREAD FILE " << fitsname(Name) << " " << Fast << std::endl;
 
     get_almsize (infile,lmaxIn,dummy);
     // std::cout << "ALMREAD FILE infile " << infile <<  " lmaxIn=" << lmaxIn << ",nside = " <<  Nside << std::endl;
 
     alloc (Nside, lmaxIn, Fast);
     // std::cout << " Lmax() = " << Lmax() << std::endl;
 
     read_Alm_from_fits(infile, *this, Lmax(),Mmax(),2);
     // std::cout << "ALMREAD FILE " << fitsname(Name) <<  " " << Lmax() << "  " << lmaxIn << std::endl;
      FastALM = Fast;
      AllocNside = Nside;
 }
 
 // ===========================
 
 void CAlmR::write(char *Name, bool Array)
 {
 	 if (Array == false)
 	 {
 	 string infile = string(fitsname(Name));
 	 remove(fitsname(Name));
 	 fitshandle out;
 	 out.create (infile);
      write_Alm_to_fits (out,*this, Lmax(), Mmax(), planckType<float32>());
 	 }
 	 else
 	 {
         fltarray A;
     	A.alloc(Lmax()+1,Lmax()+1,2);
     	for (int l=0; l <= Lmax(); l++)
 	    for (int m=0; m <= l; m++)
 	    {
 	    	 A(l,m,0) = real((*this)(l,m));
 	    	 A(l,m,1) = imag((*this)(l,m));
 	    }
 	    fits_write_fltarr(Name, A);
     }
 
      
 }
 // ===========================
  
 // ===========================
 
 void CAlmR::info(char *Name)
 {
         fltarray A,B;
     	A.alloc(Lmax()+1,Lmax()+1);
         B.alloc(Lmax()+1,Lmax()+1);
 
     	for (int l=0; l <= Lmax(); l++)
 	    for (int m=0; m <= l; m++)
 	    {
 	    	 A(l,m) = real((*this)(l,m));
 	    	 B(l,m) = imag((*this)(l,m));
 	    }
 	    std::cout << Name << " Re ==> Min = " << A.min() << ", Max = " << A.max() << ", Sig = " << A.sigma() << std::endl;
         std::cout << Name << " Im ==> Min = " << B.min() << ", Max = " << B.max() << ", Sig = " << B.sigma() << std::endl;
 }
 
 // ===========================

void CAlmR::alm_trans(Hmap<double> & map)
{
	//map.info();
    double avg=map.average();
    map.Add(-avg);
    if (map.Scheme()==NEST) map.swap_scheme();
    map2alm_iter(map, *this, Niter, weight_T);
    (*this)(0,0) += avg*sqrt(fourpi);
	if (Norm == true)
	{
	    if (Verbose == true ) std::cout << "alm_trans Norm = " << NormVal << " " <<  Lmax() << std::endl;
	    for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++)  (*this)(l,m) *=  NormVal;
	}
	if (UseBeamEff== true)
	{
 	    for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++)  (*this)(l,m) *=  BeamEff(l);
	}
    map.Add(avg);

	// if (Verbose == true) map.info();
}
// // ===========================
// 
void CAlmR::alm_rec(Hmap<double> & Map, bool KeepAlm, int RecNside)
{
    dblarray A;
    if (KeepAlm == true)
    {
    	A.alloc(Lmax()+1,Lmax()+1,2);
    	for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++) 
	    {
	    	 A(l,m,0) = real((*this)(l,m));
	    	 A(l,m,1) = imag((*this)(l,m));
	    }
    }
    
    if (UseBeamEff== true)
	{
 	    for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++)  (*this)(l,m) *=  BeamEff(l);
	}

 	if (Norm == true)
	{
	    for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++)   (*this)(l,m) *=  (1./NormVal);
	}
	
    double offset = (*this)(0,0).real()/sqrt(fourpi);
    
	(*this)(0,0) = 0;
	if  (RecNside != 0) Map.SetNside(RecNside, RING);
	else Map.SetNside(AllocNside, RING);

    alm2map(*this,Map);
    Map.Add(offset);
    
    if (KeepAlm == true)
    {
     	for (int l=0; l <= Lmax(); l++)
	    for (int m=0; m <= l; m++) 
	    {
	    	 (*this)(l,m) = xcomplex<double> (A(l,m,0), A(l,m,1));
	    }
    }

    // (*this)(0,0) = Z;
   // std::cout << "OUT Almrec,  "<< " Max = " <<  max_absalm() << std::endl;
}
 
 // ===========================
 
 void CAlmR::alm2powspec(PowSpec & powspec)
 {
 	  extract_powspec (*this, powspec);
 }
 // ===========================
 
 void CAlmR::convol(float Fwhm)
 {
    smoothWithGauss(*this, (double) Fwhm);
 }
 // ===========================
 
 void CAlmR::convol(fltarray &Filter)
 {
     int l,m;
     int LMin = Filter.nx()-1;
     if (LMin > Lmax()) LMin = Lmax();
     for (l=0; l <= LMin; l++)
     for (m=0; m <= l; m++)  (*this)(l,m) *= Filter(l);
 	for (l=LMin+1; l <= Lmax(); l++)
     for (m=0; m <= l; m++)  (*this)(l,m) *= 0.;
 }
 
 // ===========================
 
static void spline2(int size, int l, int lc,  fltarray & tab)
{
	tab.resize(size+1);
	for (int i=0; i <=size; i++)
	{
        float res =(double) i;
       res = 2.0 * l * res / (lc *size);
       float  d3=3;
       tab(i) = (3.0/2.0)*1.0 /12.0 * ( POW(ABS(res-2.), d3) - 4.0*  POW(ABS(res-1.),d3) + 6 *  POW(ABS(res), d3) - 4.0 * POW(ABS(res+1.),d3)+ POW(ABS(res+2.),d3));
	}
}
// // ===========================
// 
void CAlmR::set_beam_eff(int lmin, int lmax)
{
	 int l;
	 BeamEff.resize(Lmax()+1);
	 for (l=0;l <= lmin; l++)  BeamEff(l) = 1.;
	 int Np = lmax-lmin;
	 fltarray Tab;
     spline2(Np, 1,1, Tab);
	 for (l=lmin;l <= lmax; l++) BeamEff(l) = Tab(l-lmin);
	 for (l=lmax;l <= Lmax(); l++) BeamEff(l) = 0;
}
// // ===========================
 
 void CAlmR::wiener(PowSpec & ps_noise, PowSpec & ps_signal)
 {
  	WienerFilter.resize( Lmax()+1);
     WienerFilter(0) = 1.;
     WienerFilter(1) = 1.;
 	
     for (int l=2; l <= Lmax(); l++)
     {
   	    double Num = ps_signal.tt(l);
 	    double Den = ps_signal.tt(l) + ps_noise.tt(l);
 	    double WienFilter = (Den <= 0) ? 0: Num / Den;
 	    WienerFilter(l) = WienFilter;
      	for (int m=0; m <= l; m++) (*this)(l,m) *= WienFilter;
 	 }
 	// WienerFilter.info("WF: ");
 	// for (int l=0; l <= Lmax(); l++)  WienerFilter(l) = ps_signal.tt(l);
 	// WienerFilter.info("SIG: ");
 	// for (int l=0; l <= Lmax(); l++)  WienerFilter(l) = ps_noise.tt(l);
 	// WienerFilter.info("NOISE: ");
 }
 
 // ===========================
 
 void CAlmR::wiener(float SigmaNoise)
 {
    double Noise = SigmaNoise*SigmaNoise;
    PowSpec ps_data;
    alm2powspec(ps_data);
    WienerFilter.alloc( Lmax()+1);
    WienerFilter(0) = 1.;
     WienerFilter(1) = 1.;
    if (Verbose == true)  std::cout << "wiener: " << Lmax() << " " << Mmax() << std::endl;
     for (int l=2; l <= Lmax(); l++)
 	{
 	  double Num = ps_data.tt(l) - Noise;
 	  if (Num < 0) Num=0;
 	  double Den = ps_data.tt(l);
 	  double WienFilter = (Den <= 0) ? 0: Num / Den;
 	  WienerFilter(l) = WienFilter;
 	  for (int m=0; m <= l; m++)  (*this)(l,m) *= WienFilter;
     }
 }
// ===========================

 void mrs_alloc_powspec(PowSpec & PS, int Lmax)
 {
     arr<double> tt(Lmax+1);
     for (int l=0; l<= Lmax; ++l) tt[l] = 0.;
     PS.Set(tt);
 }
 
 //===========================
 
 void mrs_alloc_powspec(PowSpec & PS, fltarray &Cl)
 {
     int Lmax = Cl.nx()-1;
     arr<double> tt(Lmax+1);
     for (int l=0; l<= Lmax; ++l) tt[l] = Cl(l);
         PS.Set(tt);
 }
 
 
 //===========================
 
 void mrs_alloc_powspec(PowSpec & PS, dblarray &Cl)
 {
     int Lmax = Cl.nx()-1;
     arr<double> tt(Lmax+1);
     for (int l=0; l<=Lmax; ++l) tt[l] = Cl(l);
     PS.Set(tt);
 }
 
 
 // ===========================
 
 // ===========================
 
 void mrs_write_powspec(char *Name, PowSpec & PS)
 {
 	// fitshandle out;
     // remove(fitsname(Name));
 	// out.create (fitsname(Name));
 	// write_powspec_to_fits (out,PS,1);
 	fltarray Cl(PS.Lmax()+1);
 	for (int l=0; l <= PS.Lmax(); l++) Cl(l) = PS.tt(l);
     fits_write_fltarr(fitsname(Name), Cl);
 }
 
 // ===========================
 
 void mrs_read_powspec(char *Name, PowSpec & PS)
 {
 	fltarray Cl;
 	fits_read_fltarr(fitsname(Name), Cl);
     mrs_alloc_powspec(PS, Cl);
     
   	for (int l=0; l < Cl.nx(); l++)  if (l <= PS.Lmax()) PS.tt(l)=  Cl(l);
     for (int l=Cl.nx(); l <=PS.Lmax(); l++)   PS.tt(l) = 0;
 }
 
 // ===========================

 
 
