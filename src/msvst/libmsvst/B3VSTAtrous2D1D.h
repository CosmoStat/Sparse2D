/*
 * Filename : B3VSTAtrous2D1D.h
 * B3 a trous transform and reconstruction
 */
 
#ifndef _B3VSTAtrous2D1D_H
#define _B3VSTAtrous2D1D_H

#include "DefMath.h"
#include "Border.h"
#include "PoisMSVSTException.h"
#include "ImLib.h"
#include "Atrous.h"

class B3VSTAtrous2D1D
{
      float *h;               // B3 low-pass filter
      int filterLen;          // filter length of initial scale
      type_border typeBorder; // border management of data
      
      void MSVST (fltarray &data, fltarray &vstdata, int scalexy, int scalez);
      
      public:

      B3VSTAtrous2D1D (type_border tb = I_MIRROR) 
      {  
          h = B3;
          filterLen = 5;
          typeBorder = tb;
      }
      
      // convolve the input with the filter h to get output
      // convolution on the specified dimension (1 = X, 2 = Y, 3 = Z)
      // step is the distance between two non zero coefs of filter
      void convol(fltarray &input, fltarray &output, int dim, int step);
      
      // one-step decomposition
      void transformXY(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez);
      void transformZ_Approx(fltarray &data, fltarray &appr, fltarray &vstdetail, int scalexy, int scalez);
      void transformZ_Detail(fltarray &approxpxypz, fltarray &approxpxycz, fltarray approxcxypz, fltarray approxcxycz, 
                             fltarray &vstda, fltarray &vstdd, int scalexy, int scalez);
};

void B3VSTAtrous2D1D::MSVST (fltarray &data, fltarray &vstdata, int scalexy, int scalez)
{
    double b, c;
    Utils<double>::b32D1DVSTCoefs (scalexy, scalez, b, c);
    vstdata.resize(data.nx(), data.ny(), data.nz());
    int n = data.n_elem();
    for (int i=0; i<n; i++)
        vstdata(i) = sqrt(data(i) + c);
}

void B3VSTAtrous2D1D::convol(fltarray &input, fltarray &output, int dim, int step)
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
     else throw DataException("B3VSTAtrous2D1D::convol: unknown axis");
}

void B3VSTAtrous2D1D::transformXY(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez)
{
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scalexy-1);

    appr.resize(nx, ny, nz); 
    detail.resize(nx, ny, nz);

    convol(data, detail, 1, step);
    convol(detail, appr, 2, step);
    for (int i=0; i<nlen; i++) detail(i) = data(i) - appr(i);
}

void B3VSTAtrous2D1D::transformZ_Approx(fltarray &data, fltarray &appr, fltarray &vstdetail, int scalexy, int scalez)
{
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    fltarray *temp = new fltarray(nx, ny, nz);
    int nlen = data.n_elem();
    int step = POW2(scalez-1);

    appr.resize(nx, ny, nz); 
    vstdetail.resize(nx, ny, nz);

    convol(data, appr, 3, step);
    MSVST(data, vstdetail, scalexy, scalez-1);
    MSVST(appr, *temp, scalexy, scalez);
    for (int i=0; i<nlen; i++) vstdetail(i) = vstdetail(i) - (*temp)(i);
    
    delete temp; temp = NULL;
}

void B3VSTAtrous2D1D::transformZ_Detail(fltarray &approxpxypz, fltarray &approxpxycz, 
                                        fltarray approxcxypz, fltarray approxcxycz, 
                                        fltarray &vstda, fltarray &vstdd, int scalexy, int scalez)
{
    int nx = approxpxypz.nx(), ny = approxpxypz.ny(), nz = approxpxypz.nz();
    int step = POW2(scalez-1);
    fltarray *vst1 = new fltarray, *vst2 = new fltarray;
    
    MSVST(approxpxycz, *vst1, scalexy-1, scalez);
    MSVST(approxcxycz, *vst2, scalexy, scalez);
    vstda.resize(nx, ny, nz);
    vstda = *vst1;
    vstda -= *vst2;

    MSVST(approxpxypz, *vst1, scalexy-1, scalez-1);
    MSVST(approxcxypz, *vst2, scalexy, scalez-1);
    vstdd.resize(nx, ny, nz);
    vstdd = *vst1;
    vstdd -= *vst2;
    convol(vstdd, *vst1, 3, step);
    vstdd -= *vst1;
    
    delete vst1; delete vst2;
    vst1 = NULL; vst2 = NULL;
}

#endif
