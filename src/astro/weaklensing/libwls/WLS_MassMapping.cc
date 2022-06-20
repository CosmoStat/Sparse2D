/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.-L.  Starck
**
**    Date: 04/07/2020
**
**    File:  WL_WLS_MassMapping.h
**
*******************************************************************************
**
**    DESCRIPTION Weak Lensing on the Sphere package
**    -----------
**
******************************************************************************/

#include "WLS_MassMapping.h"

/*********************************************************************/


// ===========================

template<typename T> void SP_WL_Field<T>::alloc( int NsideIn, bool flag_eb, bool flag_nested)
{
    Nside = NsideIn;
    FlagEB = flag_eb;
    nested = flag_nested;
    map_G1.alloc( NsideIn, flag_nested );
    map_G2.alloc( NsideIn, flag_nested );
    map_G1.fill(0.);
    map_G2.fill(0.);

}

// ===========================

template<typename T> void SP_WL_Field<T>::alloc(Hmap<T> &MapG1, Hmap<T> &MapG2, bool flag_eb, bool flag_nested )
{
    Healpix_Ordering_Scheme NFlag;
    if(flag_nested) NFlag=NEST;
    else NFlag=RING;
    
    map_G1.alloc( MapG1.Nside(), flag_nested );
    map_G2.alloc( MapG2.Nside(), flag_nested );
    map_G1.fill(0.);
    map_G2.fill(0.);
    // cout << "SP_WL_Field alloc:nside " <<MapG2.nside() <<  MapG1.nside() << endl;
    if (MapG2.nside() == MapG1.nside())
    {
        Nside = MapG1.Nside();
        FlagEB = flag_eb;
        nested = flag_nested;
        set_map_G1(MapG1);
        set_map_G2(MapG2);
    }
    else printf(" Size of arrays are not compatible: %ld vs %ld\n. No Map allocated\n", MapG1.nside(), MapG2.nside());
    
    // cout << "SP_WL_Field alloc:nside " <<map_G1.nside() <<  map_G2.nside() << endl;

}

// ===========================

template<typename T> void SP_WL_Field<T>::read(char *NameG1, char * NameG2, bool flag_eb, bool set_ring )
{
    map_G1.read(NameG1);
    map_G2.read(NameG2);
     Nside = map_G1.nside();
     FlagEB = flag_eb;
     if( map_G1.Scheme() != DEF_MRS_ORDERING )
     {
         nested = true;
         if( map_G2.Scheme() == DEF_MRS_ORDERING )
         {
             map_G2.swap_scheme();
         }
     }//ordering = NESTED
     else
     {
         nested = false;
         if( map_G2.Scheme() != DEF_MRS_ORDERING )
         {
             map_G2.swap_scheme();
         }
     }//ordering = RING
     
     if ( (set_ring == true) && (nested == true) )
     {// Force map to ordering = RING
         map_G1.swap_scheme();
         map_G2.swap_scheme();
     }
}

// ===========================

template<typename T> void SP_WL_Field<T>::write(char *NameG1, char * NameG2)
{
    map_G1.write(NameG1);
    map_G2.write(NameG2);
}


// ===========================

template<typename T> int SP_WL_Field<T>::set_map_G1( Hmap<T> & map )
{
    int status;
    
    bool flag_map_nested = false;
    
    if ( map.Scheme() != DEF_MRS_ORDERING )
    {
        flag_map_nested = true;
    }
    
    if( (this->Nside == map.Nside()) && (this->nested == flag_map_nested) )
    {
        // cout << "SET G1" << endl;
        map_G1 = map;
        status = 1;
    }
    else
    {
        status = -1;
    }
    
    return status;
}

// ===========================

template<typename T> int SP_WL_Field<T>::set_map_G2( Hmap<T> & map )
{
    int status;
    
    bool flag_map_nested = false;
    
    if ( map.Scheme() != DEF_MRS_ORDERING )
    {
        flag_map_nested = true;
    }
    
    if( (this->Nside == map.Nside()) && (this->nested == flag_map_nested) )
    {
        map_G2 = map;
        status = 1;
    }
    else
    {
        status = -1;
    }
    return status;
}

// ===========================

template<typename T>  const  SP_WL_Field<T> & SP_WL_Field<T>::operator = (SP_WL_Field<T> & S_Map)
{
    if( get_nside() == S_Map.get_nside() )
    {
        this->map_G1 = S_Map.get_map_G1();
        this->map_G2 = S_Map.get_map_G2();
        this->FlagEB = S_Map.flag_eb();
        this->nested = S_Map.flag_nested();
    }
    return (*this);
}

// ===========================

template<typename T>  const  SP_WL_Field<T> & SP_WL_Field<T>::operator += (const SP_WL_Field<T> & S_Map)
{
    if( (this->Nside == S_Map.get_nside()) && (this->FlagEB == S_Map.flag_eb()) && (this->nested == S_Map.flag_nested()) )
    {
        this->map_G1 += S_Map.get_map_G1();
        this->map_G2 += S_Map.get_map_G2();
    }
    
    return (*this);
}

// ===========================

template<typename T>  const  SP_WL_Field<T> & SP_WL_Field<T>::operator *= (const SP_WL_Field<T> & S_Map)
{
    if( (this->Nside == S_Map.get_nside()) && (this->FlagEB == S_Map.flag_eb()) && (this->nested == S_Map.flag_nested()) )
    {
        this->map_G1 *= S_Map.get_map_G1();
        this->map_G2 *= S_Map.get_map_G2();
    }
    
    return (*this);
}

// ===========================

template<typename T>  const  SP_WL_Field<T> & SP_WL_Field<T>::operator -= (const SP_WL_Field<T> & S_Map)
{
    if( (this->Nside == S_Map.get_nside()) && (this->FlagEB == S_Map.flag_eb()) && (this->nested == S_Map.flag_nested()) )
    {
        this->map_G1 -= S_Map.get_map_G1();
        this->map_G2 -= S_Map.get_map_G2();
    }
    
    return (*this);
}

// ===========================

template<typename T>  const  SP_WL_Field<T> & SP_WL_Field<T>::operator /= (const SP_WL_Field<T> & S_Map)
{
    if( (this->Nside == S_Map.get_nside()) && (this->FlagEB == S_Map.flag_eb()) && (this->nested == S_Map.flag_nested()) )
    {
        this->map_G1 /= S_Map.get_map_G1();
        this->map_G2 /= S_Map.get_map_G2();
    }
    
    return (*this);
}

// ===========================

template<typename T> void SP_WL_Field<T>::swap_gamma_eb( bool fast)   // def fast=DEF_ALM_FAST
{
    int Lmax = 3 * Nside;
    int Mmax = Lmax;
    
    cout << "swap_gamma_eb " << endl;
    if( nested == true ) cout << "swap_gamma_eb: nested" << endl;

    arr<double> weight_T;
    weight_T.alloc( 2*Nside );
    if( fast == false )
    {
           cout << "swap_gamma_eb: read weight" << endl;
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

        map_G1.swap_scheme();
        map_G2.swap_scheme();
    }//ALM Trans Must be done in RING scheme
    
   //  G1 Min = -0.0766389, Max = 0.0752903, Mean = 1.06113e-05, sigma = 0.0022289
   //  G2 Min = -0.0692511, Max = 0.0636993, Mean = -1.47991e-05, sigma = 0.00222981

    if( FlagEB == false )
    {
        Alm<xcomplex<T> > ALM_E(Lmax,Mmax), ALM_B(Lmax,Mmax);
        cout << "spin2 " << Lmax << "  " << Mmax << endl;
        map_G1.info("G1:");
        map_G2.info("G2:");
        map2alm_spin(map_G1,map_G2,ALM_E,ALM_B,2,weight_T,false);
        map_G1.fill(0);
        alm2map( ALM_E, map_G1 );
        map_G1.info("swap_gamma_eb map_G1:");
        map_G2.fill(0);
        alm2map( ALM_B, map_G2 );
        map_G2.info("swap_gamma_eb map_G2:");
        FlagEB = true;
    }
    else
    {
        Alm<xcomplex<T> > ALM_E(Lmax,Mmax), ALM_B(Lmax,Mmax);
             
        map2alm_iter( map_G1, ALM_E, 0, weight_T );
        map2alm_iter( map_G2, ALM_B, 0, weight_T );
        alm2map_spin(ALM_E,ALM_B,map_G1,map_G2,2);
        FlagEB = false;
    }
    
    if( nested == true )
    {
        //map_T.swap_scheme();
        map_G1.swap_scheme();
        map_G2.swap_scheme();
    }//ALM Trans Must be done in RING scheme, move back in NESTED scheme
}


// ===========================

template<typename T> void SP_WL_Field<T>::swap_nested_ring()
{
    map_G1.swap_scheme();
    map_G2.swap_scheme();
    nested = !nested;//invert flag
}

// ===========================

template<typename T> void SP_WL_Field<T>::shear_fill(double value)
{
   map_G1.fill(value);
   map_G2.fill(value);
}

// ===========================

template<typename T> void SP_WL_Field<T>::shear_set_nside( int set_nside, bool set_flag_eb, bool set_flag_nested)
{
    Nside = set_nside;
    FlagEB = set_flag_eb;
    nested = set_flag_nested;
    if( nested == false )
    {
        map_G1.SetNside( set_nside, RING );
        map_G2.SetNside( set_nside, RING );
    }
    else
    {
        map_G1.SetNside( set_nside, NEST );
        map_G2.SetNside( set_nside, NEST );
    }
}

// ===========================

template<typename T> void SP_WL_Field<T>::import( SP_WL_Field<T> &S_Map, bool pessimistic)
{
    Hmap<T> map_G1_temp, map_G2_temp;
    
    map_G1_temp.alloc( S_Map.get_nside(), S_Map.flag_nested() );
    map_G2_temp.alloc( S_Map.get_nside(), S_Map.flag_nested() );
       
    map_G1_temp = S_Map.get_map_G1();
    map_G2_temp = S_Map.get_map_G2();
       
    if( S_Map.get_nside() == Nside ) // no up/degrading
    {
        map_G1.Import_nograde( map_G1_temp );
        map_G2.Import_nograde( map_G2_temp );
        
        nested = S_Map.flag_nested();
        FlagEB = S_Map.flag_eb();
    }
    else if( S_Map.get_nside() < Nside ) // upgrading
         {
             map_G1.Import_upgrade( map_G1_temp );
             map_G2.Import_upgrade( map_G2_temp );
             nested = S_Map.flag_nested();
            FlagEB = S_Map.flag_eb();
         }
         else
         {
             map_G1.Import_degrade( map_G1_temp, pessimistic );
             map_G2.Import_degrade( map_G2_temp, pessimistic );
             nested = S_Map.flag_nested();
             FlagEB = S_Map.flag_eb();
         }
}

// ===========================

template<typename T> int SP_WL_Field<T>::import_via_alm( SP_WL_Field<T> &S_Map_in, bool fast )
{
    int status;
    
    if( Nside >= S_Map_in.get_nside() )
    {
        Hmap<T> map_G1_temp, map_G2_temp;
        
        bool S_Map_in_nested = false;
        if( S_Map_in.flag_nested() == true )
        {
            S_Map_in_nested = true;
        }
    
        map_G1_temp.alloc( S_Map_in.get_nside(), S_Map_in.flag_nested() );
        map_G2_temp.alloc( S_Map_in.get_nside(), S_Map_in.flag_nested() );
        map_G1_temp = S_Map_in.get_map_G1();
        map_G2_temp = S_Map_in.get_map_G2();
           
           if( S_Map_in.flag_eb() == true )
           {//Map in FlagEB
            map_G1.import_via_alm( map_G1_temp, fast );
            map_G2.import_via_alm( map_G2_temp, fast );
           }
           else
           {//Map in TQU
               if( S_Map_in_nested == true )
               {
                   map_G1_temp.swap_scheme();
                   map_G2_temp.swap_scheme();
               }
                   
            int Lmax_in = 3 * S_Map_in.get_nside();
            int Mmax_in = Lmax_in;
    
            arr<double> weight_T;
            weight_T.alloc( 2*S_Map_in.get_nside() );
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
            Alm<xcomplex<T> > ALM_E_in( Lmax_in, Mmax_in );
            Alm<xcomplex<T> > ALM_B_in( Lmax_in, Mmax_in );
            map2alm_spin(map_G1_temp,map_G2_temp,ALM_E_in,ALM_B_in,2,weight_T,false);

            int Lmax_this = 3*Nside;
            int Mmax_this = Lmax_this;
            Alm<xcomplex<T> > ALM_E_this( Lmax_this, Mmax_this );
            ALM_E_this.SetToZero();
            Alm<xcomplex<T> > ALM_B_this( Lmax_this, Mmax_this );
            ALM_B_this.SetToZero();
            
            for( int l=0; l <= Lmax_in; l++ )
            {
                for( int m=0; m <= l; m++ )
                {
                    ALM_E_this( l, m ) = ALM_E_in( l, m );
                    ALM_B_this( l, m ) = ALM_B_in( l, m );
                }
            }
            alm2map_spin(ALM_E_this,ALM_B_this,map_G1,map_G2,2);
            if( S_Map_in_nested == true )
            {
                this->swap_nested_ring();
            }
           }
        nested = S_Map_in_nested;
        FlagEB = S_Map_in.flag_eb();
        status = 1;
    }
    else
    {
        (*this).import( S_Map_in );
        status = -1;
    }
    return status;
}

// ===========================

template<typename T> void SP_WL_Field<T>::info()
{
    if( FlagEB == false )
    {
        if( nested == false )
        {
            cout << "Map TQU, RING Scheme, Nside = " << Nside << endl;
        }
        else
        {
            cout << "Map TQU, NESTED Scheme, Nside = " << Nside << endl;
        }
        cout << "Info map G1:" << endl;
        map_G1.info();
        cout << "Info map G2:" << endl;
        map_G2.info();
    }
    else
    {
        if( nested == false )
        {
            cout << "Map FlagEB, RING Scheme, Nside = " << Nside << endl;
        }
        else
        {
            cout << "Map FlagEB, NESTED Scheme, Nside = " << Nside << endl;
        }

        cout << "Info map E:" << endl;
        map_G1.info();
        cout << "Info map B:" << endl;
        map_G2.info();
    }
}

// ===========================
// ===========================
// ===========================
// ===========================
// ===========================


// ===========================

void ShearAlm::alloc( int Nside, int u_lmax, bool Fast )
{
    ShearNside = Nside;
    niter = WL_DEF_ALM_NITER;
    
    Shear_Lmax=u_lmax;
    if ( Shear_Lmax <= 0 ) Shear_Lmax = 3*Nside;
    if ( Shear_Lmax > ALM_MAX_L_TOT ) Shear_Lmax = ALM_MAX_L_TOT;
    //int Max = L_Max;
    ShearFastALM = Fast;
    
    weight_T.alloc( 2*ShearNside );
    if( ShearFastALM == false )
    {
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[512];
        sprintf(FN, "%s/data/", HealpixFN);
        read_weight_ring( FN, ShearNside, weight_T );
        for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
     }
    else weight_T.fill(1);
     
    double lm = sqrt( (double) (ShearNside)* (double) (ShearNside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    Norm = false;
    cout << " Shear_Lmax = " << Shear_Lmax << endl;

    ALM_E.Set( Shear_Lmax, Shear_Lmax );
    ALM_B.Set( Shear_Lmax, Shear_Lmax );
    
    WienerFilter_E.resize( Shear_Lmax+1 );
    WienerFilter_E(0) = 1.;
    WienerFilter_E(1) = 1.;
    
    WienerFilter_B.resize( Shear_Lmax+1 );
    WienerFilter_B(0) = 1.;
    WienerFilter_B(1) = 1.;
}

// ===========================

void ShearAlm::set_flag_NormALM( bool flag_norm )
{
    Norm = flag_norm;
    /*
    double lm = sqrt( (double) (ShearNside)* (double) (Nside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    */
}

// ===========================

int ShearAlm::set_Alm_E( AlmR & ALM )
{
    int status;
    
    if( (Shear_Lmax == ALM.Lmax()) && (Shear_Lmax == ALM.Mmax()) )
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

int ShearAlm::set_Alm_B( AlmR & ALM )
{
    int status;
    
    if( (Shear_Lmax == ALM.Lmax()) && (Shear_Lmax == ALM.Mmax()) )
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

void ShearAlm::read( char *Name, int Nside, int u_lmax, bool Fast)
{
    ShearNside = Nside;
    niter = ALM_DEF_NITER;
    
    Shear_Lmax=u_lmax;
    if ( Shear_Lmax <= 0 ) Shear_Lmax = 3*Nside;
    if ( Shear_Lmax > ALM_MAX_L_TOT ) Shear_Lmax = ALM_MAX_L_TOT;

    ShearFastALM = Fast;
    
    weight_T.alloc( 2*ShearNside );
    if( ShearFastALM == false )
    {
        char *HealpixFN = (char *) getenv("HEALPIX");
        char FN[512];
        sprintf(FN, "%s/data/", HealpixFN);
        read_weight_ring( FN, ShearNside, weight_T );
        for(int m=0; m < (int) weight_T.size(); ++m) weight_T[m]+=1;
    }
    else weight_T.fill(1);
    
    double lm = sqrt( (double) (ShearNside)* (double) (ShearNside)*12.);
    NormVal = sqrt( (double)( lm*( lm+1 )/( 4.* PI ) ) );
    Norm = false;
    
    ALM_E.Set( Shear_Lmax, Shear_Lmax );
    ALM_B.Set( Shear_Lmax, Shear_Lmax );
    
    char *    NameFits=fitsname(Name);
    string infile = string(NameFits);
    free(NameFits);
    read_Alm_from_fits( infile, ALM_E, Shear_Lmax, Shear_Lmax, 3 );
    read_Alm_from_fits( infile, ALM_B, Shear_Lmax, Shear_Lmax, 4 );
}

// ===========================

void ShearAlm::write( char *Name )
{
    // PDT datatype = PLANCK_FLOAT32;
    char *    NameFits=fitsname(Name);
    string infile = string(NameFits);
    remove(NameFits);
    free(NameFits);
    fitshandle out;
    out.create(infile);
    write_Alm_to_fits( out, ALM_E, ALM_E.Lmax(), ALM_E.Mmax(), PLANCK_FLOAT64);
    write_Alm_to_fits( out, ALM_B, ALM_B.Lmax(), ALM_B.Mmax(), PLANCK_FLOAT64);
}

// ===========================

void ShearAlm::write_array( char *Name)
{
    char file_E_Name[512];
    char file_B_Name[512];
    
    strcpy( file_E_Name, "ALM_E_" );
    strcpy( file_B_Name, "ALM_B_" );
    
    strcat( file_E_Name, Name );
    strcat( file_B_Name, Name );

    dblarray A_E;
    A_E.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1, 2 );
    
    dblarray A_B;
    A_B.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1, 2 );
    
    for(int l=0; l <= ALM_E.Lmax(); l++)
    {
        for(int m=0; m <= l; m++)
        {
            A_E(l,m,0) = ALM_E(l,m).real();
            A_E(l,m,1) = ALM_E(l,m).imag();
                
            A_B(l,m,0) = ALM_B(l,m).real();
            A_B(l,m,1) = ALM_B(l,m).imag();
        }
    }
    fits_write_dblarr( file_E_Name, A_E );
    fits_write_dblarr( file_B_Name, A_B );
}

// ===========================

void ShearAlm::shear_alm_trans( WLS_Field & SMap, bool UseSpin2)
{
    Hmap<REAL> map_G1_temp, map_G2_temp;
    map_G1_temp = SMap.get_map_G1();
    map_G2_temp = SMap.get_map_G2();
       
    if( SMap.flag_nested() == true )
    {//Map is in NEST scheme
        //std::cout<<"SWAPPING NESTED->RING"<<std::endl;
        map_G1_temp.swap_scheme();
        map_G2_temp.swap_scheme();
    }

    if( UseSpin2 == true )
    {
       //std::cout<<"TQU TO FlagEB THRU ALMS"<<std::endl;
        //Map is in G1G2
        map2alm_spin(map_G1_temp,map_G2_temp,ALM_E,ALM_B,2,weight_T,false);
    }
    else
    {
        //std::cout<<"TQU TO FlagEB THRU ALMS"<<std::endl;
        //Map is in FlagEB
        map2alm_iter( map_G1_temp, ALM_E, niter, weight_T);
        map2alm_iter( map_G2_temp, ALM_B, niter, weight_T);
    }
    
    if( Norm == true )
    {
        for(int l=0; l <= Shear_Lmax; l++)
        {
            for(int m=0; m <= l; m++)
            {
                ALM_E(l,m) *= NormVal;
                ALM_B(l,m) *= NormVal;
            }
        }
    }
    ALM_E(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_E(1,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(1,0)=xcomplex<REAL> ( 0.,0.);
    ALM_E(1,1)=xcomplex<REAL> ( 0.,0.);
    ALM_B(1,1)=xcomplex<REAL> ( 0.,0.);
}

// ===========================

void ShearAlm::shear_alm_rec( WLS_Field & SMap, bool UseSpin2, bool KeepAlm, int RecNside)
{
    Hmap<REAL> map_G1_temp, map_G2_temp;
    
    if( RecNside != 0 )
    {
        map_G1_temp.alloc( RecNside );
        map_G2_temp.alloc( RecNside );
        ShearNside = RecNside;
        SMap.shear_set_nside( RecNside, SMap.flag_eb(), SMap.flag_nested() );
    }
    else
    {
        map_G1_temp.alloc( ShearNside );
        map_G2_temp.alloc( ShearNside );
    }
    
    dblarray A_E, A_B;
    if( KeepAlm == true )
    {
        A_E.alloc( Shear_Lmax+1, Shear_Lmax+1, 2 );
        A_B.alloc( Shear_Lmax+1, Shear_Lmax+1, 2 );
        
        for(int l=0; l <= Shear_Lmax; l++)
        {
            for(int m=0; m <= l; m++)
            {
                A_E(l,m,0) = ALM_E(l,m).real();
                A_E(l,m,1) = ALM_E(l,m).imag();
                
                A_B(l,m,0) = ALM_B(l,m).real();
                A_B(l,m,1) = ALM_B(l,m).imag();
            }
        }
    }
     
    if( Norm == true )
    {
        for(int l=0; l <= Shear_Lmax; l++)
        {
            for(int m=0; m <= l; m++)
            {
                ALM_E(l,m) *= (1./NormVal);
                ALM_B(l,m) *= (1./NormVal);
            }
        }
    }
    ALM_E(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_E(1,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(1,0)=xcomplex<REAL> ( 0.,0.);
    ALM_E(1,1)=xcomplex<REAL> ( 0.,0.);
    ALM_B(1,1)=xcomplex<REAL> ( 0.,0.);
    if( UseSpin2 == true )
    {
        alm2map_spin(ALM_E,ALM_B,map_G1_temp,map_G2_temp,2);
    }
    else
    {
        alm2map( ALM_E, map_G1_temp );
        alm2map( ALM_B, map_G2_temp );
    }
    
    if( SMap.flag_nested() == true )
    {
        map_G1_temp.swap_scheme();
        map_G2_temp.swap_scheme();
    }
    SMap.set_map_G1( map_G1_temp );
    SMap.set_map_G2( map_G2_temp );
    
    if( KeepAlm == true )
    {
         for(int l=0; l <= Shear_Lmax; l++)
        {
            for(int m=0; m <= l; m++)
            {
                ALM_E(l,m)  = xcomplex<REAL> ( A_E(l,m,0),A_E(l,m,1));
                ALM_B(l,m)  = xcomplex<REAL> ( A_B(l,m,0), A_B(l,m,1));
            }
        }
    }
}

// ===========================

void ShearAlm::shear_alm2powspec( PowSpec & specE, PowSpec & specB)
{
    extract_powspec (ALM_E, specE);
    extract_powspec (ALM_B, specB);
}

// ===========================

void ShearAlm::info()
{
    fltarray A_E_re, A_E_im, A_B_re, A_B_im;
    A_E_re.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1 );
    A_E_im.alloc( ALM_E.Lmax()+1, ALM_E.Lmax()+1 );
    A_B_re.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1 );
    A_B_im.alloc( ALM_B.Lmax()+1, ALM_B.Lmax()+1 );
       
    for(int l=0; l <= Shear_Lmax; l++)
    {
        for(int m=0; m <= l; m++)
        {
            A_E_re(l,m) = ALM_E(l,m).real();
            A_E_im(l,m) = ALM_E(l,m).imag();
                
            A_B_re(l,m) = ALM_B(l,m).real();
            A_B_im(l,m) = ALM_B(l,m).imag();
        }
    }
    cout << "Alm E: " << " Re ==> Min = " << A_E_re.min() << ", Max = " << A_E_re.max() << ", Sig = " << A_E_re.sigma() << endl;
    cout << "Alm E: " << " Re ==> Min = " << A_E_re.min() << ", Max = " << A_E_re.max() << ", Sig = " << A_E_re.sigma() << endl;
    cout << "Alm B: " << " Re ==> Min = " << A_B_re.min() << ", Max = " << A_B_re.max() << ", Sig = " << A_B_re.sigma() << endl;
    cout << "Alm B: " << " Re ==> Min = " << A_B_re.min() << ", Max = " << A_B_re.max() << ", Sig = " << A_B_re.sigma() << endl;
}

// ===========================

void ShearAlm::convol(fltarray &Filter_E, fltarray &Filter_B )
{
    int l,m;
    int LMin_E = Filter_E.nx()-1;
    int LMin_B = Filter_B.nx()-1;
    
    if( LMin_E > Shear_Lmax )
    {
        LMin_E = Shear_Lmax;
    }
    if( LMin_B > Shear_Lmax )
    {
        LMin_B = Shear_Lmax;
    }
    
    for( l=0; l <= LMin_E; l++)
    {
        for( m=0; m <= l; m++)
        {
            ALM_E(l,m) *= Filter_E(l);
        }
    }
    for( l=LMin_E+1; l <= Shear_Lmax; l++)
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
    for( l=LMin_B+1; l <= Shear_Lmax; l++)
    {
        for( m=0; m <= l; m++)
        {
            ALM_B(l,m) *= 0.;
        }
    }
}

// ===========================

void ShearAlm::convol( fltarray &Filter )
{
    int l,m;
    
    int LMin = Filter.nx()-1;
    
    if( LMin > Shear_Lmax )
    {
        LMin = Shear_Lmax;
    }
        
    for( l=0; l <= LMin; l++)
    {
        for( m=0; m <= l; m++)
        {
            ALM_E(l,m) *= Filter(l);
            ALM_B(l,m) *= Filter(l);
        }
    }
    for( l=LMin+1; l <= Shear_Lmax; l++)
    {
        for( m=0; m <= l; m++)
        {
            ALM_E(l,m) *= 0.;
            ALM_B(l,m) *= 0.;
        }
    }
}

// ===========================

void ShearAlm::convol( float Fwhm )
{
   smoothWithGauss(ALM_E, (double) Fwhm );
   smoothWithGauss(ALM_B, (double) Fwhm );
}

// ===========================

void ShearAlm::set_wiener_filter( PowSpec & ps_noise, PowSpec & ps_signal )
{
    double PS,PN;
    WienerFilter_E(0) = WienerFilter_E(1) = 1;
    WienerFilter_B(0) = WienerFilter_B(1) = 1;
    int NlPS =ps_signal.Lmax()+1;
    int NlPN =ps_noise.Lmax()+1;
    for( int l=2; l <= Shear_Lmax; l++ )
    {
        PS = PN = 0;
        if (l < NlPS) PS = ps_signal.tt(l);
        if (l < NlPN) PN = ps_noise.tt(l);
        double Num = PS;
        double Den = PS + PN;
        double WienFilter = (Den <= 0) ? 0: Num / Den;
        WienerFilter_E(l) = WienFilter;
        // In theorie, we should have differents signal powspec for E and B
        // here, we reuse the E-filter for B as well, to perform the same
        // smoothing on both component.
        WienerFilter_B(l) = WienFilter;
    }
}

// ===========================

void ShearAlm::alm_mult_wiener()
{
    for( int l=0; l <= Shear_Lmax; l++ )
    for( int m=0; m <= l; m++ )
    {
        ALM_E(l,m) *= WienerFilter_E(l);
        ALM_B(l,m) *= WienerFilter_B(l);
    }
    ALM_E(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(0,0)=xcomplex<REAL> ( 0.,0.);
}

// ===========================

void ShearAlm::mult_wiener(WLS_Field & SMap, bool UseSpin2Trans, bool UseSpin2Rec)
{
    shear_alm_trans(SMap, UseSpin2Trans); // Alm transformation of a Healpix map
    for( int l=0; l <= Shear_Lmax; l++ )
    for( int m=0; m <= l; m++ )
    {
        ALM_E(l,m) *= WienerFilter_E(l);
        ALM_B(l,m) *= WienerFilter_B(l);
    }
    ALM_E(0,0)=xcomplex<REAL> ( 0.,0.);
    ALM_B(0,0)=xcomplex<REAL> ( 0.,0.);
    shear_alm_rec(SMap, UseSpin2Rec); // Alm inverse transform
}

// ===========================

void ShearAlm::wiener(WLS_Field & SMap, PowSpec & ps_noise, PowSpec & ps_signal, bool UseSpin2)
{
    set_wiener_filter(ps_noise, ps_signal);
    mult_wiener(SMap, UseSpin2, UseSpin2);
}

// ===========================

double ShearAlm::max_absalm_FlagEB()
{
    double Max =0.;
    Max=MAX(max_absalm_B(),max_absalm_E());
    return Max;
}

// ===========================

double ShearAlm::max_absalm_E()
{
    double Max=0.;
    
    for( int l=2; l <= Shear_Lmax; l++ )
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

double ShearAlm::max_absalm_B()
{
    double Max=0.;
    
    for( int l=2; l <= Shear_Lmax; l++ )
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

int ShearAlm::hard_threshold(float lambda_e, float lambda_b, int & MaxNonZeroL_E, int & MaxNonZeroL_B )
{
    int Cpt=0;
    MaxNonZeroL_E = 0;
    MaxNonZeroL_B = 0;
    
    for( int l=1; l <= Shear_Lmax; l++ )
    {
        for( int m=0; m <= l; m++ )
        {
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

int ShearAlm::soft_threshold(float lambda_e, float lambda_b, int & MaxNonZeroL_E, int & MaxNonZeroL_B )
{
    int Cpt=0;
    MaxNonZeroL_E = 0;
    MaxNonZeroL_B = 0;
    
    double norm_E = 0.0;
    double norm_B = 0.0;
    
    double Coef_t, Coef_e, Coef_b;
    
    for( int l=1; l <= Shear_Lmax; l++ )
    {
        for( int m=0; m <= l; m++ )
        {
            norm_E = sqrt( norm( ALM_E(l,m) ) );
            norm_B = sqrt( norm( ALM_B(l,m) ) );
            
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

// ===========================

void WLS_MassMapping::alloc(Hdmap &G1, Hdmap &G2, Hdmap &CovMatrix, Hdmap &MaskData, int NScale, int lmax)
{
    GammaData.alloc (G1,G2, false,false);

    WT.set_alm_iter(CAlm.get_niter());

    Nside = G1.Nside();
    Npix = G1.Npix();
    Mask.alloc(Nside, DEF_MRS_ORDERING);
    CovMat.alloc(Nside, DEF_MRS_ORDERING);
    MinCovMat = WL_INFINITE_COV_VALUE;
    FieldTemp.alloc(Nside);

    if (Verbose)
         cout << "WLS_MassMapping alloc: G1 Nside = " << Nside << ", Npix = " << Npix << endl;
 
    // CovMat is the diagonal cov matrix on each indivual shear component,
    // while CovMatrix is the diagonal cov matrix on the complex shear field
    for (int i=0; i < Npix; i++)
    {
        Mask[i] = MaskData[i];
        if ((G1[i]==0) && (G2[i]==0))
            Mask[i]=0;
        if ((CovMatrix[i] >= WL_INFINITE_COV_VALUE) || (Mask[i] == 0)) CovMat[i] = WL_INFINITE_COV_VALUE;
        else
        {
            CovMat[i] = CovMatrix[i] / 2.;
            if (CovMat[i] == 0 )
                CovMat[i] = WL_INFINITE_COV_VALUE;
            else if (CovMat[i]  < MinCovMat)
                 MinCovMat = CovMat[i];
        }
        if (CovMat[i] >= WL_INFINITE_COV_VALUE)
             Mask[i] =0;
        GammaData.map_G1[i] =G1[i];
        GammaData.map_G2[i] =G2[i];
    }
    if (Verbose)
    {
        cout << "CovMat Nside = " << Nside << ", MinCovMat= " << MinCovMat << endl;
    }
             
    WeightGamma_Wiener.alloc(Nside, DEF_MRS_ORDERING);
    for (int i=0; i < Npix; i++)
        WeightGamma_Wiener[i] = MinCovMat / CovMat[i];
    
    WeightGamma_Sparse.alloc(Nside, DEF_MRS_ORDERING);
    for (int i=0; i < Npix; i++)
        WeightGamma_Sparse[i] = sqrt(MinCovMat / CovMat[i]);

    bool fast=false;
    CAlm.alloc(Nside, lmax, fast);

    // Wavelet initialization
    int NbrScale;
    if (NScale >= 2) NbrScale=NScale;
    else  NbrScale = (int) (log((double) Nside) / log(2.)-2);
    WT.set_alm_iter(CAlm.get_niter());
    bool Tight_Frame=true;
    WT.wt_alloc(Nside, NbrScale, CAlm.get_lmax(), DEF_MRS_ORDERING, Tight_Frame);
    TabActivCoefE.alloc(Npix, NbrScale);
    TabActivCoefB.alloc(Npix, NbrScale);
    if (Verbose)
    {
        cout << "LMax = "  <<  CAlm.get_lmax() << ", Nscale = " << WT.nscale() << endl;
    }
}

// ===========================

void WLS_MassMapping::alloc(Hdmap &G1, Hdmap &G2, Hdmap &CovMatrix, int NbrScale, int lmax)
{
    Hdmap MaskData;
    MaskData.alloc(G1.nside(), DEF_MRS_ORDERING);
    MaskData.fill(1.);
    alloc(G1, G2, CovMatrix, MaskData,  NbrScale, lmax);
}

// ===========================

void WLS_MassMapping::alloc(Hdmap &G1, Hdmap &G2, float SigmaNoise, Hdmap &MaskData, int NbrScale, int lmax)
{
    double Var =SigmaNoise*SigmaNoise;
    Hdmap CovMatrix;
    CovMatrix.alloc(G1.nside(), DEF_MRS_ORDERING);
    CovMatrix.fill(Var);
    alloc(G1, G2, CovMatrix, MaskData,  NbrScale,  lmax);
}

// ===========================

void WLS_MassMapping::alloc(Hdmap &G1, Hdmap &G2, float SigmaNoise, int NbrScale, int lmax)
{
    Hdmap MaskData;
    MaskData.alloc(G1.nside(), DEF_MRS_ORDERING);
    MaskData.fill(1.);
    alloc(G1, G2, SigmaNoise, MaskData,  NbrScale, lmax);
}

// ===========================

void WLS_MassMapping::gamma2eb(WLS_Field & Gamma, WLS_Field & Kappa, bool KeepAlm)
{
    // cout << "Gamma nside = " <<Gamma.get_nside() << endl;
    if ((Kappa.get_nside() ==0) || (Kappa.get_nside() !=Gamma.get_nside()))
            Kappa.alloc( Gamma.get_nside(), false, Gamma.flag_nested());

    CAlm.shear_alm_trans(Gamma, true);
    CAlm.shear_alm_rec(Kappa, false, KeepAlm);
    // Kappa = Gamma;
    // Kappa.swap_gamma_eb(false);
}

// ===========================

void WLS_MassMapping::eb2gamma(WLS_Field & Kappa, WLS_Field & Gamma, bool KeepAlm)
{
    // cout << "Gamma nside = " <<Gamma.get_nside() << endl;
    if ((Gamma.get_nside() ==0) || (Kappa.get_nside() != Gamma.get_nside()))
            Gamma.alloc( Kappa.get_nside(), false, Kappa.flag_nested());
    else Gamma.shear_fill(0);
    
    CAlm.shear_alm_trans(Kappa, false);
    CAlm.shear_alm_rec(Gamma, true, KeepAlm);
    // Kappa = Gamma;
    // Kappa.swap_gamma_eb(false);
}

// ===========================

void WLS_MassMapping::sp_kaiser_squires(WLS_Field & Gamma, WLS_Field & Kappa, float Fwhm)
{
    // cout << "Gamma nside = " <<Gamma.get_nside() << endl;
    if ((Kappa.get_nside() ==0) || (Kappa.get_nside() !=Gamma.get_nside()))
            Kappa.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    else Kappa.shear_fill(0);

    CAlm.shear_alm_trans(Gamma, true);
    if (Fwhm > 0)
    {
        CAlm.convol(Fwhm);
    }
    CAlm.shear_alm_rec(Kappa, false);
}

// ===========================

void WLS_MassMapping::get_residual_gamma(WLS_Field & Kappa, WLS_Field & Resi_Gamma, wl_type_weight TypeWeight)
{
    eb2gamma(Kappa, FieldTemp);
    for (int i=0; i < Npix; i++)
    {
        Resi_Gamma.map_G1[i] = (GammaData.map_G1[i] - FieldTemp.map_G1[i]) * Mask[i];
        Resi_Gamma.map_G2[i] = (GammaData.map_G2[i] - FieldTemp.map_G2[i]) * Mask[i];
        if (TypeWeight == SIGMA_WEIGHT)
        {
            Resi_Gamma.map_G1[i] *= WeightGamma_Sparse[i];
            Resi_Gamma.map_G2[i] *= WeightGamma_Sparse[i];
        }
        else if (TypeWeight == VAR_WEIGH)
        {
            Resi_Gamma.map_G1[i] *=  WeightGamma_Wiener[i];
            Resi_Gamma.map_G2[i] *=  WeightGamma_Wiener[i];
        }
    }
}

// ===========================

void WLS_MassMapping::get_residual_eb(WLS_Field & Kappa, WLS_Field & Resi_EB, wl_type_weight TypeWeight)
{
    get_residual_gamma(Kappa, FieldTemp, TypeWeight);
    gamma2eb(FieldTemp, Resi_EB);
}

// ===========================

void WLS_MassMapping::wiener(WLS_Field & SMap, PowSpec & ps_noise, PowSpec & ps_signal, bool UseSpin2Trans, bool UseSpin2Rec)
{
    CAlm.set_wiener_filter(ps_noise, ps_signal);
    CAlm.mult_wiener(SMap, UseSpin2Trans, UseSpin2Rec);
}

// ===========================

void WLS_MassMapping::iter_wiener(WLS_Field & Gamma, PowSpec & ps_signal, WLS_Field & Kappa, int NiterWiener)
{
    WLS_Field Resi_EB;
    if ((Kappa.get_nside() ==0) || (Kappa.get_nside() !=Gamma.get_nside()))
            Kappa.alloc( Gamma.get_nside(), false, Gamma.flag_nested());

    Resi_EB.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    
    PowSpec ps_noise;
    mrs_alloc_powspec(ps_noise, CAlm.get_lmax());
    for (int i=0; i <= CAlm.get_lmax(); i++)
         ps_noise.tt(i) = MinCovMat;

    int lms = ps_signal.Lmax();
    
    CAlm.set_wiener_filter(ps_noise, ps_signal);
    // fits_write_dblarr( "xx_filter_wiener_E.fits", CAlm.WienerFilter_E );
    // fits_write_dblarr( "xx_filter_wiener_B.fits", CAlm.WienerFilter_E );

    // Kappa.shear_fill(0);
    for (int i=0; i < NiterWiener; i++)
    {
        // enum wl_type_weight {NO_WEIGHT, SIGMA_WEIGHT, VAR_WEIGH};
        get_residual_eb(Kappa, Resi_EB, VAR_WEIGH);
        Kappa.map_G1 += Resi_EB.map_G1;
        Kappa.map_G2 += Resi_EB.map_G2;
        // wiener(Kappa, ps_noise, ps_signal, bool UseSpin2Trans, bool UseSpin2Rec);
        CAlm.mult_wiener(Kappa, false, false);
        if (Verbose)
        {
            cout << "Iter " << i+1 << ": SigmaResi = " << Resi_EB.map_G1.sigma()  << ", " <<  Resi_EB.map_G2.sigma() << endl;
            Kappa.map_G1.info("  kappa-E");
            Kappa.map_G2.info("  kappa-B");
        }
    }
}


// ===========================

void WLS_MassMapping::sparse_reconstruction(WLS_Field & Gamma, WLS_Field & Kappa, float NSigma, int NiterSparse, int FirstDetectScale)
{
    float SigmaNoise = sqrt(MinCovMat);
    bool UseMad=false;
    bool KillLastScale=false;
    WLS_Field Resi_EB;
    if ((Kappa.get_nside() ==0) || (Kappa.get_nside() !=Gamma.get_nside()))
            Kappa.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    Resi_EB.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    
    // Get the set of Active coefficients
    
    for (int i=0; i < NiterSparse; i++)
    {
        // enum wl_type_weight {NO_WEIGHT, SIGMA_WEIGHT, VAR_WEIGH};
        get_residual_eb(Kappa, Resi_EB, SIGMA_WEIGHT);
        Kappa.map_G1 += Resi_EB.map_G1;
        Kappa.map_G2 += Resi_EB.map_G2;
        WT.hard_thresholding(Kappa.map_G1, NSigma, SigmaNoise,UseMad, KillLastScale, FirstDetectScale);
        WT.hard_thresholding(Kappa.map_G2, NSigma, SigmaNoise,UseMad, KillLastScale,FirstDetectScale);
        
        if (Verbose)
        {
            cout << "Iter " << i+1 << ": SigmaResi = " << Resi_EB.map_G1.sigma()  << ", " <<  Resi_EB.map_G2.sigma() << endl;
        }
    }
}

// ===========================

void WLS_MassMapping::get_active_coef(WLS_Field & Kappa, float NSigma, bool KillLastScale, bool OnlyPos, bool NoSparseBMode, int FirstDetectScale)
{
    float SigmaNoise = sqrt(MinCovMat);
    WT.transform(Kappa.map_G1,WT_BandLimit);
    
    cout << " FirstDetectScale = " <<FirstDetectScale << endl;
    
    for (int b=0; b < FirstDetectScale; b++)
        for (int i=0; i < WT.WTTrans.nx(); i++)
        {
                 TabActivCoefE(i,b) = 0;
                TabActivCoefB(i,b)= 0;
        }
    
    for (int b=FirstDetectScale; b < WT.nscale()-1; b++)
    {
        float Level = SigmaNoise * NSigma * WT.TabNorm(b);
        for (int i=0; i < WT.WTTrans.nx(); i++)
        {
            if (OnlyPos == true) TabActivCoefE(i,b)= (WT.WTTrans(i,b) < Level) ? 0:1;
            else TabActivCoefE(i,b)= (ABS(WT.WTTrans(i,b)) < Level) ? 0:1;
        }
    }
    
    if (NoSparseBMode == false)
    {
        WT.transform(Kappa.map_G2, WT_BandLimit);
        for (int b=FirstDetectScale; b < WT.nscale()-1; b++)
        {
            float Level = SigmaNoise * NSigma * WT.TabNorm(b);
            for (int i=0; i < WT.WTTrans.nx(); i++)
            {
                if (OnlyPos == true) TabActivCoefB(i,b)= (WT.WTTrans(i,b) < Level) ? 0:1;
                else TabActivCoefB(i,b)= (ABS(WT.WTTrans(i,b)) < Level) ? 0:1;            }
        }
    }
    
    int LastScaleVal = (KillLastScale == false) ? 1: 0;
    int b=WT.nscale()-1;
    for (int i=0; i < WT.WTTrans.nx(); i++)
    {
        TabActivCoefE(i,b) = LastScaleVal;
        TabActivCoefB(i,b) = (NoSparseBMode == false) ? LastScaleVal: 0;
    }
    
    for (int b=0; b < WT.nscale(); b++)
    {
        int Cpt = 0;
        for (int i=0; i < WT.WTTrans.nx(); i++)
            if (TabActivCoefE(i,b) > 0) Cpt++;
        cout << "BAND " << b+1 << ": Cpt = " << Cpt << ", Nx = " << WT.WTTrans.nx() << endl;
    }

}
// ===========================

void WLS_MassMapping::mult_sparse(WLS_Field & KappaSparse, bool NoSparseBMode)
{
    WT.transform(KappaSparse.map_G1, WT_BandLimit);
    for (int b=0; b < WT.nscale(); b++)
    for (int i=0; i < WT.WTTrans.nx(); i++)
        WT.WTTrans(i,b) *= TabActivCoefE(i,b);
    WT.recons(KappaSparse.map_G1, WT_BandLimit);
    if (NoSparseBMode == false)
    {
        WT.transform(KappaSparse.map_G2, WT_BandLimit);
        for (int b=0; b < WT.nscale(); b++)
        for (int i=0; i < WT.WTTrans.nx(); i++)
            WT.WTTrans(i,b) *= TabActivCoefB(i,b);
        WT.recons(KappaSparse.map_G2, WT_BandLimit);
    }
    else KappaSparse.map_G2.fill(0);
}


// ===========================

void WLS_MassMapping::mcalens(WLS_Field & Gamma, PowSpec & ps_signal, WLS_Field & Kappa, WLS_Field & KappaSparse, float NSigma, int NiterSparse, bool SparsePositivityConstraint, int FirstDetectScale)
{
    float SigmaNoise = sqrt(MinCovMat);
    bool UseMad=false;
    bool KillLastScale=false;
    bool OnlyPos = true;
    bool NoSparseBMode = true;
    
    if (Verbose)
    {
        cout << "MCALens: Niter = " << NiterSparse << ", NSigma=" <<  NSigma  << ", SparsePos=" <<  SparsePositivityConstraint << ", NoSparseBMode = " << NoSparseBMode << ", FirstDetectScale = " << FirstDetectScale+1 << endl;

    }
    WLS_Field Resi_EB, KappaWiener;
    if ((Kappa.get_nside() ==0) || (Kappa.get_nside() !=Gamma.get_nside()))
            Kappa.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    Resi_EB.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    KappaSparse.alloc( Gamma.get_nside(), false, Gamma.flag_nested());
    KappaWiener.alloc( Gamma.get_nside(), false, Gamma.flag_nested());

    PowSpec ps_noise;
    mrs_alloc_powspec(ps_noise, CAlm.get_lmax());
    for (int i=0; i <= CAlm.get_lmax(); i++)
         ps_noise.tt(i) = MinCovMat;

    int lms = ps_signal.Lmax();
    CAlm.set_wiener_filter(ps_noise, ps_signal);

    get_residual_eb(Kappa, Resi_EB, SIGMA_WEIGHT);
    get_active_coef(Resi_EB, NSigma, true, OnlyPos, NoSparseBMode, FirstDetectScale);
    
    for (int i=0; i < NiterSparse; i++)
    {
        get_residual_eb(Kappa, Resi_EB, SIGMA_WEIGHT);
        mult_sparse(Resi_EB, NoSparseBMode);
        KappaSparse.map_G1 += Resi_EB.map_G1;
        KappaSparse.map_G2 += Resi_EB.map_G2;
            
        for (int i=0; i < Npix; i++)
        {
            if (SparsePositivityConstraint)
            {
                if (KappaSparse.map_G1[i] < 0) KappaSparse.map_G1[i] = 0;
                if (KappaSparse.map_G2[i] < 0) KappaSparse.map_G2[i] = 0;
            }
            Kappa.map_G1[i] = KappaWiener.map_G1[i] + KappaSparse.map_G1[i];
            if (NoSparseBMode)
                 Kappa.map_G2[i] = KappaSparse.map_G2[i] = 0;
            else Kappa.map_G2[i] = KappaWiener.map_G2[i] + KappaSparse.map_G2[i];
        }
        get_residual_eb(Kappa, Resi_EB, VAR_WEIGH);
        KappaWiener.map_G1 += Resi_EB.map_G1;
        KappaWiener.map_G2 += Resi_EB.map_G2;
        CAlm.mult_wiener(KappaWiener, false, false);
        
        for (int i=0; i < Npix; i++)
        {
            Kappa.map_G1[i] = KappaWiener.map_G1[i] + KappaSparse.map_G1[i];
            if (NoSparseBMode)
                Kappa.map_G2[i] = KappaWiener.map_G2[i] = 0;
            else Kappa.map_G2[i] = KappaWiener.map_G2[i] + KappaSparse.map_G2[i];
        }
        
        if (Verbose)
        {
            cout << "Iter " << i+1 << ": SigmaResi = " << Resi_EB.map_G1.sigma()  << ", " <<  Resi_EB.map_G2.sigma() << endl;
        }
    }
    
    if (Verbose)
    {
         KappaSparse.map_G1.info((char*) "Sparse E mode");
         KappaWiener.map_G1.info((char*) "Wiener E mode");
         Kappa.map_G1.info((char*) "Kappa E mode");
         KappaSparse.map_G2.info((char*) "Sparse B mode");
         KappaWiener.map_G2.info((char*) "Wiener B mode");
         Kappa.map_G2.info((char*) "Kappa B mode");
    }
}

// ===========================
