/*
 * Filename : Atrous.h
 * A trous transform and reconstruction
 */
 
#ifndef _ATROUS_H
#define _ATROUS_H

#include "DefMath.h"
#include "Border.h"
#include "PoisMSVSTException.h"

// order = 0
static float B0[2] = { 1./2., 1./2. };
// order = 3
static float B3[5] = {1./16., 1./4., 3./8., 1./4., 1./16.};
// order = -1
static float B0BH[6] = { -1./16, 1./16., 1./2., 1./2., 1./16., -1./16. };

class SplineAtrous
{
      float* h;               // low-pass filter
      int filterLen;          // filter length of initial scale
      float corr;             // correction factor
      bool ver2;              // the 2nd version of filtering
                              // h0 = h            h1 = h
                              // g0 = Id - h*h     g1 = Id
      type_border typeBorder; // border management of data
      
      void initSplineAtrous (int splineOrder, float c, type_border tb, bool newver)
      {
          if (splineOrder == 0)
          {
              h = B0;
              filterLen = 2;
          }
          else if (splineOrder == 3)
          {
              h = B3;
              filterLen = 5;
          }
          else if (splineOrder == -1)
          {
               h = B0BH;
               filterLen = 6;
          }
          else
          {
              h = NULL;
              filterLen = 0;
          }
          corr = c;
          typeBorder = tb;
          ver2 = newver;
      }

      public:

      SplineAtrous ()
      {
          initSplineAtrous (0, 1., I_MIRROR, false);
      }
      
      SplineAtrous (int splineOrder)
      {
          initSplineAtrous (splineOrder, 1., I_MIRROR, false);
      }
      
      SplineAtrous (int splineOrder, float c, type_border tb, bool newver) 
      {  
          initSplineAtrous (splineOrder, c, tb, newver);
      }
      
      // convolve the input with the filter h to get output
      // convolution on the specified dimension (0 = X, 1 = Y, 2 = Z)
      // step is the distance between two non zero coefs of filter
      void convol(fltarray &input, fltarray &output, int dim, int step);
      
      // one-step decomposition 1D, 2D and 3D
      void transform(fltarray &data, fltarray &appr, fltarray &detail, int scale);
      // one-step reconstruction 1D, 2D and 3D
      void recons(fltarray &appr, fltarray &detail, fltarray &data, int scale);
};

void SplineAtrous::convol(fltarray &input, fltarray &output, int dim, int step)
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
     else throw DataException("SplineAtrous::convol: unknown axis");
}

void SplineAtrous::transform(fltarray &data, fltarray &appr, fltarray &detail, int scale)
{
    int dim = data.naxis();
    int nx = data.nx(), ny = data.ny(), nz = data.nz();
    int nlen = data.n_elem();
    int step = POW2(scale-1);
    appr.resize(nx, ny, nz); detail.resize(nx, ny, nz);
    
    if (dim == 1)
    {
        convol(data, appr, 1, step);
        if (ver2)
        {
            convol(appr, detail, 1, step);
            for (int i=0; i<nlen; i++) detail(i) = data(i) - detail(i);
        }
        else
            for (int i=0; i<nlen; i++) detail(i) = data(i) - appr(i);
    }
    else if (dim == 2)
    {
        float c = 1./corr;
        convol(data, detail, 1, step);
        convol(detail, appr, 2, step);
        if (ver2)
        {
            fltarray *temp = new fltarray(nx, ny, nz);
            convol(appr, detail, 1, step);
            convol(detail, *temp, 2, step);
            for (int i=0; i<nlen; i++) detail(i) = c * data(i) - (*temp)(i);
            delete temp; temp = NULL;
        }
        else
            for (int i=0; i<nlen; i++) detail(i) = c * data(i) - appr(i);
    }
    else if (dim == 3)
    {
        float c = 1./ (corr * corr);
        convol(data, appr, 1, step);
        convol(appr, detail, 2, step);
        convol(detail, appr, 3, step);
        if (ver2)
        {
            fltarray *temp = new fltarray(nx, ny, nz);
            convol(appr, detail, 1, step);
            convol(detail, *temp, 2, step);
            convol(*temp, detail, 3, step);
            for (int i=0; i<nlen; i++) detail(i) = c * data(i) - detail(i);
            delete temp; temp = NULL;
        }
        else
            for (int i=0; i<nlen; i++) detail(i) = c * data(i) - appr(i);
    }
    else throw DataException("SplineAtrous::transform: unknown axis");
}


void SplineAtrous::recons(fltarray &appr, fltarray &detail, fltarray &data, int scale)
{
    int dim = appr.naxis();
    int nx = appr.nx(), ny = appr.ny(), nz = appr.nz();
    int nlen = appr.n_elem();
    int step = POW2(scale-1);
    data.resize(nx, ny, nz);
    
    if (dim == 1)
    {
        if (ver2)
        {
            convol(appr, data, 1, step);            
            data += detail;
        }
        else
            for (int i=0; i<nlen; i++) data(i) = appr(i) + detail(i);
    }
    else if (dim == 2)
    {
        float c = corr;
        if (ver2)
        {
            fltarray *temp = new fltarray(nx, ny, nz);
            convol(appr, *temp, 1, step);            
            convol(*temp, data, 2, step);            
            for (int i=0; i<nlen; i++) data(i) = c * (data(i) + detail(i));
            delete temp; temp = NULL;
        }
        else
            for (int i=0; i<nlen; i++) data(i) = c * (appr(i) + detail(i));
    }
    else if (dim == 3)
    {
        float c = corr * corr;
        if (ver2)
        {
            fltarray *temp = new fltarray(nx, ny, nz);
            convol(appr, data, 1, step);            
            convol(data, *temp, 2, step);            
            convol(*temp, data, 3, step);            
            for (int i=0; i<nlen; i++) data(i) = c * (data(i) + detail(i));
            delete temp; temp = NULL;
        }
        else
            for (int i=0; i<nlen; i++) data(i) = c * (appr(i) + detail(i));
    }
    else throw DataException("SplineAtrous::recons: unknown axis");
}
#endif
