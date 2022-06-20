/*##########################################################################
pySAP - Copyright (C) CEA, 2017 - 2018
Distributed under the terms of the CeCILL-B license, as published by
the CEA-CNRS-INRIA. Refer to the LICENSE file or to
http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
for details.
##########################################################################*/

#ifndef ALM_H_
#define ALM_H_
 
// Includes
#include <iostream>
#include <string>
#include <sstream>
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Prob.h"
#include "numpydata.hpp"
#include "HealpixClass.h"
 

class C_ALM {
public:
    CAlmR alm;

    // Constructor
    C_ALM(bool verbose=false);

    // Destructor
    ~C_ALM(){};
 
    // Information method
    void info();
    
    // Allocation for a given nside.
    // by default, number of scales is automatically estimated
    void alloc(int Nside, int lmax=0, bool Verbose=false)
    {
        cout << "ALM ALLOC" << Nside << " " << lmax << endl;
        bool Fast=DEF_ALM_FAST;// false
        alm.alloc(Nside,lmax,Fast);
    }

    // Wavelet transform
    void trans(py::array_t<float>& arr);
    py::array_t<float> recons();
    void hthres(float Thres);
    double maxalm() { return alm.max_absalm();}
    void use_l2norm(bool l2norm=false) {alm.Norm = l2norm;}
    py::array_t<xcomplex<REAL>>  get_tabalm();
    void put_tabalm(py::array_t<xcomplex<REAL>> arr);
    py::array_t<xcomplex<REAL>>  get_alm();
    void put_alm(py::array_t<xcomplex<REAL>> arr);
    void wiener(py::array_t<float>& psn,py::array_t<float>& pss);
    py::array_t<float> alm2spec();

private:
    std::string m_opath;
    int nb_procs;
    bool Verbose;
};

// Constructor
C_ALM::C_ALM(bool verbose)
{
    cout << "ALM INIT" << endl;
    // Define instance attributes
    this->Verbose = verbose;
    this->nb_procs=0;
    // The maximum number of threads returned by omp_get_max_threads()
    // (which is the default number of threads used by OMP in parallel
    // regions) can sometimes be far below the number of CPUs.
    // It is then required to set it in relation to the real number of CPUs
    // (the -1 is used to live one thread to the main process which ensures
    // better and more constant performances). - Fabrice Poupon 2013/03/09
    #ifdef _OPENMP
        if (nb_procs <= 0)
            this->nb_procs = omp_get_num_procs() - 1;
        else
            this->nb_procs = nb_procs;
        omp_set_num_threads(this->nb_procs);
    #endif
}

// Wavelet transform
void C_ALM::trans(py::array_t<float>& arr)
{
    Hmap<REAL> Map = array2hmap(arr);
    int Npix = arr.shape(0);
    int Nside =floor(sqrt((double) Npix / 12.));
    cout << "ALM ALLOC" << endl;
    this->alloc(Nside);
    alm.alm_trans(Map);
}
py::array_t<float> C_ALM::recons()
{
    Hmap<REAL> Map;
    alm.alm_rec(Map);
    return hmap2array(Map);
}
void C_ALM::hthres(float T)
{
    int MaxNonZeroL=0;
    int NonZeroNbr = alm.hard_threshold(T,MaxNonZeroL);
}

py::array_t<xcomplex<REAL>> C_ALM::get_alm()
{
    int Nl=alm.Lmax()+1;
    int Nalm = (alm.Lmax()+1)*(alm.Lmax()+2)/2;
    auto arr = py::array_t<xcomplex<REAL>>(Nalm);
    auto buffer = arr.request();
    xcomplex<REAL> *pointer = (xcomplex<REAL> *) buffer.ptr;
    int i=0;
    for (int l=0; l <= alm.Lmax(); l++)
    for (int m=0; m <= l; m++)
        pointer[i++] = alm(l,m);
    return arr;
}

void C_ALM::put_alm(py::array_t<xcomplex<REAL>> arr)
{
    int Nl=alm.Lmax()+1;
    int nelem= Nl*Nl;
    auto buffer = arr.request();
    xcomplex<REAL> *pointer = (xcomplex<REAL> *) buffer.ptr;
    int i=0;
    for (int l=0; l <= alm.Lmax(); l++)
    for (int m=0; m <= l; m++)
        alm(l,m) = pointer[i++];
}

py::array_t<xcomplex<REAL>> C_ALM::get_tabalm()
{
    int Nl=alm.Lmax()+1;
    int nelem= Nl*Nl;
    auto arr = py::array_t<xcomplex<REAL>>(nelem);
    auto buffer = arr.request();
    xcomplex<REAL> *pointer = (xcomplex<REAL> *) buffer.ptr;

    for (int l=0; l <= alm.Lmax(); l++)
    for (int m=0; m <= l; m++)
        pointer[l + m * Nl] = alm(l,m);
    arr.resize({Nl, Nl});
    return arr;
}

void C_ALM::put_tabalm(py::array_t<xcomplex<REAL>> arr)
{
    int Nl=alm.Lmax()+1;
    int nelem= Nl*Nl;
    auto buffer = arr.request();
    xcomplex<REAL> *pointer = (xcomplex<REAL> *) buffer.ptr;
    for (int l=0; l <= alm.Lmax(); l++)
    for (int m=0; m <= l; m++)
        alm(l,m) = pointer[l + m * Nl];
}

void C_ALM::wiener(py::array_t<float>& psn, py::array_t<float>& pss)
{
    fltarray Cls,Cln;
    int Nl = alm.Lmax()+1;
    auto buffer = pss.request();
    float *pointer = (float *) buffer.ptr;
    for (int i=0; i<Nl;i++) Cls(i) = pointer[i];
    buffer = psn.request();
    pointer = (float *) buffer.ptr;
    for (int i=0; i<Nl;i++) Cln(i) = pointer[i];
    PowSpec ps_signal,ps_noise;
    mrs_alloc_powspec(ps_signal, Cls);
    mrs_alloc_powspec(ps_noise, Cln);
    alm.wiener(ps_noise,ps_signal);
}

py::array_t<float> C_ALM::alm2spec()
{
    PowSpec ps_data;
    alm.alm2powspec(ps_data);
    int Nl = alm.Lmax()+1;
    auto arr = py::array_t<xcomplex<REAL>>(Nl);
    auto buffer = arr.request();
    xcomplex<REAL> *pointer = (xcomplex<REAL> *) buffer.ptr;
    for (int l=0; l < Nl; l++)
        pointer[l] = ps_data.tt(l);
    return arr;
}

// py::array_t<float> alm2powspec();


/*
 
py::array_t<float> C_ALM::hard_filtering(py::array_t<float> arr, float NSigma, float & SigmaNoise, bool UseMad, bool KillLastScale)
{
    if (arr.ndim() != 1)
           throw std::runtime_error("Input should be 1-D NumPy array");
    auto buffer = arr.request();
    float *pointer = (float *) buffer.ptr;
       
    int Npix = arr.shape(0);
    int Nside =floor(sqrt((double) Npix / 12.));
    if (Nside != WT.nside())
        throw std::runtime_error("Input map has not expected nside.");

    if ((!this->mr_initialized) || (Nside != WT.nside()))
    {
        int Ns=0;
        alloc(Nside, Ns);
    }
    
    // Copy the numpy array into a healpix map
    Hdmap Map;
    Map.alloc(Nside);
    for (int i=0; i < Npix; i++)  Map[i] = pointer[i];
       
    // Wavelet transform denoising
    WT.hard_thresholding(Map, NSigma,SigmaNoise,UseMad, KillLastScale);
    
    auto arr1 = py::array_t<float>(Npix);
    auto buf1 = arr1.request();
    pointer = (float *) buf1.ptr;
    for (int i=0; i<Npix; i++)
                pointer[i] = Map[i];
    return arr1;
}
*/
#endif
