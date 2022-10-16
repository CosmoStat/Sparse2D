/*
 * Filename : Atrous2D1D.h
 * 2D+1D B3 A trous transform and reconstruction
 */
 
#ifndef _ATROUS2D1D_H
#define _ATROUS2D1D_H

#include "DefMath.h"
#include "Border.h"
#include "PoisMSVSTException.h"
#include "ImLib.h"
#include "Atrous.h"

class Atrous2D1D
{
      float *h;               // B3 low-pass filter
      int filterLen;          // filter length of initial scale
      type_border typeBorder; // border management of data
      bool ver2;              // the 2nd version of filtering
                              // h0 = h            h1 = h
                              // g0 = Id - h*h     g1 = Id
      
      public:

      Atrous2D1D (type_border tb = I_MIRROR, bool newver = false) 
      {  
          h = B3;
          filterLen = 5;
          typeBorder = tb;
          ver2 = newver;
      }
      
      // convolve the input with the filter h to get output
      // convolution on the specified dimension (1 = X, 2 = Y, 3 = Z)
      // step is the distance between two non zero coefs of filter
      void convol(fltarray &input, fltarray &output, int dim, int step);
      
      // one-step decomposition
      void transformXY(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez);
      void transformZ(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez);

      // one-step reconstruction
      void reconsXY(fltarray &appr, fltarray &detail, fltarray &data, int scalexy, int scalez);
      void reconsZ(fltarray &appr, fltarray &detail, fltarray &data, int scalexy, int scalez);
};

void Atrous2D1D::convol(fltarray &input, fltarray &output, int dim, int step)
{    
     int nx = input.nx(), ny = input.ny(), nz = input.nz();
     output.resize(nx, ny, nz);
     bool even = ((filterLen % 2) == 0);
     
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
                  c = (even && (step!=1)) ? step/2 + (filterLen/2-1-p)*step : (filterLen/2-p)*step;
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
                  c = (even && (step!=1)) ? step/2 + (filterLen/2-1-p)*step : (filterLen/2-p)*step;
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
                  c = (even && (step!=1)) ? step/2 + (filterLen/2-1-p)*step : (filterLen/2-p)*step;
                  index = z - c;
                  index = get_index(index, nz, typeBorder);
                  output(x, y ,z) += input(x, y, index) * h[filterLen-p-1];
             } 
         }     
     }
     else throw DataException("Atrous2D1D::convol: unknown axis");
}

void Atrous2D1D::transformXY(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez)
{
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scalexy-1);

    appr.resize(nx, ny, nz); 
    detail.resize(nx, ny, nz);
    
    convol(data, detail, 1, step);
    convol(detail, appr, 2, step);
    if (ver2)
    {
        fltarray *temp = new fltarray(nx, ny, nz);
        convol(appr, detail, 1, step);
        convol(detail, *temp, 2, step);        
        for (int i=0; i<nlen; i++) detail(i) = data(i) - (*temp)(i);
        delete temp; temp = NULL;
    }
    else
        for (int i=0; i<nlen; i++) detail(i) = data(i) - appr(i);
}

void Atrous2D1D::transformZ(fltarray &data, fltarray &appr, fltarray &detail, int scalexy, int scalez)
{
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scalez-1);

    appr.resize(nx, ny, nz); 
    detail.resize(nx, ny, nz);

    convol(data, appr, 3, step);
    if (ver2)
    {
        convol(appr, detail, 3, step);
        for (int i=0; i<nlen; i++) detail(i) = data(i) - detail(i);
    }
    else
        for (int i=0; i<nlen; i++) detail(i) = data(i) - appr(i);
}

void Atrous2D1D::reconsXY(fltarray &appr, fltarray &detail, fltarray &data, int scalexy, int scalez)
{
    int nx = appr.nx(), ny = appr.ny(), nz = appr.nz();
    int nlen = appr.n_elem();
    int step = POW2(scalexy-1);
    data.resize(nx, ny, nz);

    if (ver2)
    {
        fltarray *temp = new fltarray(nx, ny, nz);
        convol(appr, *temp, 1, step);            
        convol(*temp, data, 2, step);            
        data += detail;
        delete temp; temp = NULL;
    }
    else
        for (int i=0; i<nlen; i++) data(i) = appr(i) + detail(i);
}

void Atrous2D1D::reconsZ(fltarray &appr, fltarray &detail, fltarray &data, int scalexy, int scalez)
{
    int nx = appr.nx(), ny = appr.ny(), nz = appr.nz();
    int nlen = appr.n_elem();
    int step = POW2(scalez-1);
    data.resize(nx, ny, nz);

    if (ver2)
    {
        convol(appr, data, 3, step);            
        data += detail;
    }
    else
        for (int i=0; i<nlen; i++) data(i) = appr(i) + detail(i);
}

#endif
