/*
 * Filename : B3VSTAtrous.h
 * B3 a trous transform and reconstruction
 */
 
#ifndef _B3VSTATROUS_H
#define _B3VSTATROUS_H

#include "DefMath.h"
#include "Border.h"
#include "Wavelet.h"
#include "Atrous.h"

class B3VSTAtrous
{
      float *h;               // B3 low-pass filter
      int filterLen;          // filter length of initial scale
      type_border typeBorder; // border management of data
            
      public:

      B3VSTAtrous (type_border tb = I_MIRROR) 
      {  
          h = B3;
          filterLen = 5;
          typeBorder = tb;
      }

      // MSVST data -- VST --> vstdata
      void MSVST (fltarray &data, fltarray &vstdata, int scale, bool coupled=true);
      // Inverse MSVST vstdata -- VST^{-1} --> data (without positive projection)
      // bias correction is recommended while reconstruct the data after estimation
      void InvMSVST (fltarray &vstdata, fltarray &data, int scale, bool cbias, fltarray *intensity = NULL);

      // convolve the input with the filter h to get output
      // convolution on the specified dimension (1 = X, 2 = Y, 3 = Z)
      // step is the distance between two non zero coefs of filter
      void convol(fltarray &input, fltarray &output, int dim, int step);
      
      // one-step decomposition 1D, 2D and 3D used in coupled detection and estimation
      void transform(fltarray &data, fltarray &appr, fltarray &apprvst, fltarray &detail, int scale);
      // one-step decomposition 1D, 2D and 3D used in separated detection and estimation
      void transform(fltarray &data, fltarray &appr, fltarray &detail, int scale);
      // one-step reconstruction 1D, 2D and 3D used in coupled detection and estimation
      void recons(fltarray &apprvst, fltarray &detail, fltarray &datavst, int scale);
};

void B3VSTAtrous::MSVST (fltarray &data, fltarray &vstdata, int scale, bool coupled)
{
    double b, c;
    Utils<double>::b3VSTCoefs (data.naxis(), scale, b, c);
    vstdata.resize(data.nx(), data.ny(), data.nz());
    int n = data.n_elem();
    for (int i=0; i<n; i++)
    {
        if (coupled)
            vstdata(i) = sqrt(data(i) + c); // tau_1 = 1
        else
            vstdata(i) = b * sqrt(data(i) + c);
    }
}

void B3VSTAtrous::InvMSVST (fltarray &vstdata, fltarray &data, int scale, bool cbias, fltarray *intensity)
{
    double b, c;
    Utils<double>::b3VSTCoefs (vstdata.naxis(), scale, b, c);
    double cb = cbias ? 1. / (b * b) : 0.;
    data.resize(vstdata.nx(), vstdata.ny(), vstdata.nz());
    int n = vstdata.n_elem();
    for (int i=0; i<n; i++)
        data(i) = vstdata(i)*vstdata(i) + cb * (intensity == NULL ? 1. : MAX(0.,(*intensity)(i))/(MAX(0.,(*intensity)(i)) + c)) - c; // tau_1 = 1    
}

void B3VSTAtrous::convol(fltarray &input, fltarray &output, int dim, int step)
{    
     int nx = input.nx(), ny = input.ny(), nz = input.nz();
     output.resize(nx, ny, nz);
     
     if (dim == 1) // X axis
     {
         int leny = MAX(1, ny), lenz = MAX(1, nz);
         int index, c;
         for (int y=0; y<leny; y++)
         for (int z=0; z<lenz; z++)
         for (int x=0; x<nx; x++)
         {
             output(x, y, z) = 0;
             for (int p=0; p<filterLen; p++)
  	         {
                  c = (filterLen/2-p)*step;
                  index = x - c;
                  index = get_index(index, nx, typeBorder);
                  output(x, y ,z) += input(index, y, z) * h[filterLen-p-1];
             } 
         }     
     }
     else if (dim == 2) // Y axis
     {
         int lenz = MAX(1, nz);
         int index, c;
         for (int x=0; x<nx; x++)
         for (int z=0; z<lenz; z++)
         for (int y=0; y<ny; y++)
         {
             output(x, y, z) = 0;
             for (int p=0; p<filterLen; p++)
  	         {
                  c = (filterLen/2-p)*step;
                  index = y - c;
                  index = get_index(index, ny, typeBorder);
                  output(x, y ,z) += input(x, index, z) * h[filterLen-p-1];
             } 
         }     
     }
     else if (dim == 3) // Z axis
     {
         int index, c;
         for (int x=0; x<nx; x++)
         for (int y=0; y<ny; y++)
         for (int z=0; z<nz; z++)
         {
             output(x, y, z) = 0;
             for (int p=0; p<filterLen; p++)
  	         {
                  c = (filterLen/2-p)*step;
                  index = z - c;
                  index = get_index(index, nz, typeBorder);
                  output(x, y ,z) += input(x, y, index) * h[filterLen-p-1];
             } 
         }     
     }
     else throw DataException("B3VSTAtrous::convol: unknown axis");
}

void B3VSTAtrous::transform(fltarray &data, fltarray &appr, fltarray &apprvst, fltarray &detail, int scale)
{
    int dim = data.naxis();
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scale-1);
    appr.resize(nx, ny, nz); 
    apprvst.resize(nx, ny, nz); 
    detail.resize(nx, ny, nz);
    fltarray *vst = new fltarray(nx, ny, nz);
    
    if (dim == 1)
    {
        MSVST(data, *vst, scale-1);
        convol(data, appr, 1, step);
        MSVST(appr, apprvst, scale);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - apprvst(i);
    }
    else if (dim == 2)
    {
        MSVST(data, *vst, scale-1);
        convol(data, detail, 1, step);
        convol(detail, appr, 2, step);
        MSVST(appr, apprvst, scale);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - apprvst(i);
    }
    else if (dim == 3)
    {
        MSVST(data, *vst, scale-1);
        convol(data, appr, 1, step);
        convol(appr, detail, 2, step);
        convol(detail, appr, 3, step);
        MSVST(appr, apprvst, scale);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - apprvst(i);
    }
    else throw DataException("B3VSTAtrous::transform: unknown axis");
    
    delete vst; vst = NULL;
}

void B3VSTAtrous::transform(fltarray &data, fltarray &appr, fltarray &detail, int scale)
{
    int dim = data.naxis();
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scale-1);
    appr.resize(nx, ny, nz); 
    detail.resize(nx, ny, nz);
    fltarray *vst = new fltarray(nx, ny, nz);
    
    if (dim == 1)
    {
        MSVST(data, *vst, scale-1, false);        
        convol(*vst, appr, 1, step);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - appr(i);
        convol(data, appr, 1, step);
    }
    else if (dim == 2)
    {
        MSVST(data, *vst, scale-1, false);
        convol(*vst, detail, 1, step);
        convol(detail, appr, 2, step);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - appr(i);
        convol(data, *vst, 1, step);
        convol(*vst, appr, 2, step);
    }
    else if (dim == 3)
    {
        MSVST(data, *vst, scale-1, false);
        convol(*vst, appr, 1, step);
        convol(appr, detail, 2, step);
        convol(detail, appr, 3, step);
        for (int i=0; i<nlen; i++) detail(i) = (*vst)(i) - appr(i);
        convol(data, appr, 1, step);
        convol(appr, *vst, 2, step);
        convol(*vst, appr, 3, step);
    }
    else throw DataException("B3VSTAtrous::transform: unknown axis");
    
    delete vst; vst = NULL;
}

void B3VSTAtrous::recons(fltarray &apprvst, fltarray &detail, fltarray &datavst, int scale)
{
    datavst = apprvst;
    datavst += detail;
}

#endif
