/******************************************************************************
**                   Copyright (C) 2005 by CEA + Valencia observatory
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Enn Saar, Jean-Luc Starck and Vicent Martinez
**
**    Date:  27/02/03
**    
**    File:  Genus3D.h
**
*******************************************************************************
**
**    DESCRIPTION  genus program
**    ----------- 
**                 
******************************************************************************/
      

#ifndef _GENUS3D_H_
#define _GENUS3D_H_

/* The 4 Minkowski functionals 'MF', calculated by the Crofton formula,
 * see Schmalzing & Buchert, ApJ 482, L1-L4, 1997 (astro-ph/9702130),
 * and Coles, Davies & Pearson, MN 281, 1375-1384, 1996.
 *  After thresholding at a given level, the input data are an integer 0-1 
 *  grid density array 'dat' (should, in principle, be Boolean, 
 *  but does it make the code faster?).
 * We march over the grid and count all filled vertices (NV), 
 * the edges from a vertex towards the the three coordinate directions (NE), 
 * the faces along the coordinate planes (NF), 
 * and the cube (NC) from the vertex (all in case if these are formed by
 * the filled vertices). Using the positive coordinate directions ensures
 * that we are counting all basics only once.
 */
 
//ES The formulas for all four functionals for a smoothed Gaussian field
//ES are in Schmalzing & Buchert, ApJ 482, L1-L4, astro-ph/9702130.
//ES I list these formulas here:
//ES MF_0=0.5-0.5*erf(nu/sqrt(2));
//ES MF_1=2./3.*lambda/sqrt(2*M_PI)*exp(-nu*nu/2);
//ES MF_2=2./3.*lambda**lambda/sqrt(2*M_PI)*nu*exp(-nu*nu/2);
//ES MF_3=lambda*lambda**lambda/sqrt(2*M_PI)*(nu*nu-1)*exp(-nu*nu/2);
//ES lambda=Lambda/sqrt(2*M_PI);
//ES \Lambda^2=-\frac{\xi''(0)}{\xi(0)}=\frac13\langle k^2\rangle,
//ES and for a Gaussian field with the power spectrum
//ES P(k)\sim k^n and a normal Gaussian filter G(r)\sim exp(-r^2/2R^2),
//ES Lambda=1./R*sqrt((n+3)/6).

void mf_curves(fltarray &dens, fltarray &ResMinFun, float MFSTEP, Bool Periodic,  
               Bool  Segmentation, float Step, Bool Verb);
//          MF(0..3) == Minkovski function V0-V3
//      MinFun(p, 0) = MF(0);   // First Minkowski functional = Area of the surface.
//      MinFun(p, 1) = MF(1);   // Second Minkowski functional = Volume enclosed by the surface.
//      MinFun(p, 2) = MF(2);   // Third Minkowski functional = Integrated mean curvature of the surface.
//      MinFun(p, 3) = MF(3);   // Fourth Minkowski functional = Euler characteristic 
//                                     Genus  = 1 - MF(3] / 2
//      MinFun(p, 4) = vf;      // Percent. of Volume larger than the level, 
//      MinFun(p, 5) = ro;      // Threshold Level
//      MinFun(p, 6) = NuXAbs;  // Nu = xerf(vf)


void genpois(fltarray &dens);
// Generate a Poisson noise

extern void convolve(fltarray & Data, fltarray & Gauss);

void read_data(fltarray & Data, char *Name_Cube_In, float BinCat=1);
//  read input data and test if it is a catalogue or a fits

void gengauss(fltarray &dens, double NuIndex, float POWA);
// Generate a randonm Gaussian field with index NuIndex and spectral amplitude POWA
#endif
