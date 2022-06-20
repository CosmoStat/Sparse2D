/*##########################################################################
pySAP - Copyright (C) CEA, 2017 - 2018
Distributed under the terms of the CeCILL-B license, as published by
the CEA-CNRS-INRIA. Refer to the LICENSE file or to
http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
for details.
##########################################################################*/

#ifndef MRS_H_
#define MRS_H_
 
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
#include "MRS_Sparse.h"

class MRS {
public:
    // Constructor
    MRS(bool verbose=false);

    // Destructor
    ~MRS();
 
    // Information method
    void info();
    
    // Allocation for a given nside.
    // by default, number of scales is automatically estimated
    void alloc(int Nside, int Nscale=0, int Lmax=0, int ALM_iter=0, bool Verbose=false);

    // Wavelet transform
    py::array_t<float> uwt(py::array_t<float>& arr, int Nscale=0);

    // Reconstruction method
    py::array_t<float> iuwt(py::array_t<float>  mr_data);

    // Getter/setter functions for the input/output image path
    void set_opath(std::string path) {this->m_opath = path;}
    string get_opath() const {return m_opath;}

    py::array_t<float> get_tabnorm();

    int nside() {return WT.nside();}
    
    int nscale() {return WT.nscale();}

    void set_alm_iter(int ALM_IT) { WT.ALM_iter=ALM_IT; }
    
    void hthres_scale(int Scale, float NSigma, float & SigmaNoise, bool UseMad=false)
         { WT.hard_thresholding(Scale,NSigma,SigmaNoise, UseMad);}

    py::array_t<float>  hard_filtering(py::array_t<float> DataIn, float NSigma, float & SigmaNoise, bool UseMad=false, bool KillLastScale=false);

private:
    C_UWT WT;
    bool mr_initialized;
    std::string m_opath;
    int nb_procs;
    bool Verbose;
};

// Constructor
MRS::MRS(bool verbose)
{
    // Define instance attributes
    this->Verbose = verbose;
    this->mr_initialized=false;
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

// Destructor
MRS::~MRS(){
    mr_initialized=False;
    nb_procs=0;
    Verbose=False;
}

// int MRS::get_lmax(int  & Lmax, int Nside, float ZeroPadding)
// {
//    return mrs_get_lmax(Lmax, Nside, ZeroPadding);
// }

void MRS::info(){
    // Information message
    cout << "---------" << endl;
    cout << "Information" << endl;
    cout << "Runtime parameters:" << endl;
    cout << "  Number of procs: " << this->nb_procs << endl;
    cout << "---------" << endl;
}

// Transform method
void  MRS::alloc(int Nside, int Nscale, int LmaxIn, int ALM_iter, bool Verb)
{
    if (Verb) (*this).info();
    bool nested=false;
    float ZeroPadding=0.;
    int LM=LmaxIn;
    int Lmax = mrs_get_lmax (LM,  Nside, ZeroPadding);
    cout << "ALLOC: " << Nside << " " << Nscale << " " << LmaxIn << " " << ALM_iter << endl;
    WT.wt_alloc(Nside, Nscale, Lmax, nested);  //  DEF_MRS_ORDERING);
    WT.ALM_iter = ALM_iter;
    if (Verb) Verbose=Verb;
    this->mr_initialized = true;
}

// Transform method
py::array_t<float>  MRS::uwt(py::array_t<float>& arr, int Ns)
{
    Verbose=True;
    cout << "Input arr.ndim() = " << arr.ndim() << " " << arr.shape(0) << endl;
    if (arr.ndim() != 1)
        throw std::runtime_error("Input should be 1-D NumPy array");

    
    cout << "OK" << endl;
    exit(-1);
    
    
    auto buffer = arr.request();
    float *pointer = (float *) buffer.ptr;
 
    int Npix = arr.shape(0);
    int Nside =floor(sqrt((double) Npix / 12.));
    if ((!this->mr_initialized) || ((Ns>1) && (Ns != WT.nscale()))
            || (Nside != WT.nside()))
        alloc(Nside, Ns);
    
    if (Nside != WT.nside())
    {
        throw std::runtime_error("Input map has not expected nside.");
    }

    cout << "alloc ok. Nside = " << Nside << " " << Ns << ", Npix = " << Npix << endl;

    // Copy the numpy array into a healpix map
    Hdmap Map;
    Map.alloc(Nside);
    for (int i=0; i < Npix; i++)  Map[i] = pointer[i];
    
    if (Verbose == True)
     cout << "WTTRANS: Nside = " << Nside << ", Nscale = "<< WT.nscale() << endl;
 
    // Wavelet transform
    WT.transform(Map);
    cout << "transform ok. " << endl;

    // Recopy to numpy array
    auto arr1 = py::array_t<float>(Npix*WT.nscale());
    auto buf1 = arr1.request();
    pointer = (float *) buf1.ptr;
    for(int s=0; s < WT.nscale(); s++)
    for (int i=0; i<Npix; i++)
    {
        pointer[i + s * Npix] = WT.WTTrans(i,s);
    }
    cout << "copy ok. Npix =" << Npix << endl;

    arr1.resize({WT.nscale(), Npix});
    cout << "resize dim = " << arr1.ndim() << " " << arr1.shape(1)<< " " <<   arr1.shape(0) << endl;
    return arr1;
    
}

py::array_t<float> MRS::get_tabnorm()
{
    auto arr1 = py::array_t<float>(WT.TabNorm.nx());
    auto buf1 = arr1.request();
    float *pointer = (float *) buf1.ptr;
    for (int i=0; i<WT.TabNorm.nx(); i++)
            pointer[i] = WT.TabNorm(i);
    return arr1;
}

// Reconstruction method
py::array_t<float> MRS::iuwt(py::array_t<float>  mr_data)
{
    if (mr_data.ndim() != 2)
        throw std::runtime_error("Input should be 2-D NumPy array");
    auto buffer = mr_data.request();
    float *pointer = (float *) buffer.ptr;
    
    int Ns = mr_data.shape(0);
    int Npix = mr_data.shape(1);
    int Nside =floor(sqrt((double) Npix / 12.));
    if (this->Verbose > 0)
    {
        cout << "Starting Reconstruction" << endl;
        cout << "Runtime parameters:" << endl;
        cout << "  Number of bands: " << Ns << endl;
        cout << "  Npix per band : " << Npix << endl;
    }

    if ((!this->mr_initialized) || ((Ns>1) && (Ns != WT.nscale()))
            || (Nside != WT.nside()))
         alloc(Nside, Ns);
    
    for (int i=0; i<Npix; i++)
    {
        for(int s=0; s < WT.nscale(); s++)
            WT.WTTrans(i,s) = pointer[i + s * Npix];
    }
    
    // Wavelet reconstruction
    Hdmap Map;
    Map.alloc(Nside);
    WT.recons(Map);
    
    auto arr1 = py::array_t<float>(Npix);
    auto buf1 = arr1.request();
    pointer = (float *) buf1.ptr;
    for (int i=0; i<Npix; i++)
            pointer[i] = Map[i];

    return arr1;
}


py::array_t<float> MRS::hard_filtering(py::array_t<float> arr, float NSigma, float & SigmaNoise, bool UseMad, bool KillLastScale)
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

#endif
