/******************************************************************************
**                   Copyright (C) 2003 by CEA + Valencia observatory
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
**    File:  genus.cc
**
*******************************************************************************
**
**    DESCRIPTION  genus program
**    ----------- 
**                 
******************************************************************************/
   
//ES are comments by ES
  
#include "Array.h"
#include "IM_IO.h"
#include "NR.h"
#include "FFTN.h"
#include "FFTN_3D.h"
#include "DefPoint.h"
#include "CUBE_Segment.h"
#include "Genus3D.h"
#include "DefFunc.h"

/***************************************/

static int genus_segmentation(intarray & dat, float NuXAbs)
{
   int N1 = dat.nx();
   int N2 = dat.ny();
   int N3 = dat.nz();
   int i,j,k;
   Bool CleanBord = False;
   int FirstLabel = 1;
   fltarray Cube(N1,N2,N3);
   intarray Segment(N1,N2,N3);
   float Threshold = FLOAT_EPSILON;
   int NbrIsol = 0;
   int NbrHole = 0;

   for(i=0;i<N1;i++)
   for(j=0;j<N2;j++)
   for(k=0;k<N3;k++) Cube(i,j,k) = dat(i,j,k);   
   cube_segment(Cube, Segment, NbrIsol,  Threshold,  CleanBord, FirstLabel);
 
   // Reverse the values.
   for(i=0;i<N1;i++) 
   for(j=0;j<N2;j++) 
   for(k=0;k<N3;k++) Cube(i,j,k) = (dat(i,j,k) == 0) ? 1: 0;
   cube_segment (Cube, Segment,  NbrHole,  Threshold,  CleanBord, FirstLabel);
   if (NuXAbs >= 0)  return (NbrHole - NbrIsol + 1);
   else return (NbrIsol - NbrHole + 1);
}

/***************************************/

static void densort(fltarray &sdens)
{
   float *Buffer = sdens.buffer() -1;
   unsigned long N = sdens.nx()*sdens.ny()*sdens.nz();
   sort(N, Buffer);
   //  qsort(sdens.buffer(), N, sizeof(float), comp);
}


/***************************************/

//ES This function is new, until the next //ES comment
/* The 4 Minkowski functionals 'MF', calculated by the Crofton formula,
 * see Schmalzing & Buchert, ApJ 482, L1-L4, 1997 (astro-ph/9702130),
 * and Coles, Davies & Pearson, MN 281, 1375-1384, 1996.
 * The input data are an integer 0-1 grid density array 'dat'
 * (should, in principle, be Boolean, but does it make the code faster?).
 * We march over the grid and count all filled vertices (NV), 
 * the edges from a vertex towards the the three coordinate directions (NE), 
 * the faces along the coordinate planes (NF), 
 * and the cube (NC) from the vertex (all in case if these are formed by
 * the filled vertices). Using the positive coordinate directions ensures
 * that we are counting all basics only once.
 */
static void minfun(intarray &dat, fltarray &MF, Bool Periodic)
{
   int N1 = dat.nx();
   int N2 = dat.ny();
   int N3 = dat.nz();
   int i,j,k,i1,j1,k1;
   int NC,NF,NE,NV; 
   int ne,nf;
   int v100,v010;
   float V;


	N1--; N2--; N3--; 		// the array was padded
	V=N1*N2*N3;
	NC=NF=NE=NV=0; // cubes, faces, edges & vertices

	for(i=0;i<N1;i++) {
		if((i1=i+1)==N1&&Periodic==True) i1=0;		
		for(j=0;j<N2;j++) {
			if((j1=j+1)==N2&&Periodic==True) j1=0;		
			for(k=0;k<N3;k++) {
				if(!dat(i,j,k)) continue;
				NV++;				// found a vertex	
				ne=nf=0;			// local edges and faces		
				if((v100=dat(i1,j,k))) ne++; // found an edge
				if((v010=dat(i,j1,k))) ne++;
				NE+=ne; 
				if((ne==2)&&dat(i1,j1,k)) nf++;	// found a face
				if((k1=k+1)==N3&&Periodic==True) k1=0;		
				if(dat(i,j,k1)) {
					NE++;
					if(v100&&dat(i1,j,k1)) nf++;
					if(v010&&dat(i,j1,k1)) nf++;
					}
				NF+=nf;
				if((nf==3)&&dat(i1,j1,k1)) NC++; // found a cube
				}
			}
		}
// The Crofton formula is NC/V, but NV/V is a better volume estimator
	MF(0)=NV/V;
// Give MF-s for the whole brick, not densities;
// these can always be normalized later, and MF-s are easier to understand
	MF(1)=(2./9.)*(-3*NC+NF);
	MF(2)=(2./9.)*(3*NC-2*NF+NE);
	MF(3)=(-NC+NF-NE+NV);		
//	MF(1)=(2./9.)*(-3*NC+NF)/V;
//	MF(2)=(2./9.)*(3*NC-2*NF+NE)/V;
//	MF(3)=(-NC+NF-NE+NV)/V;		
}
//ES end of minfun()



//ES changed the name, as said above
void mf_curves(fltarray &dens, fltarray &ResMinFun, float MFSTEP, Bool Periodic, 
                Bool  Segmentation, float Step, Bool Verb)
// ResMinFun curve 
// dens = three dimensional array
// ResMinFun = two dimensional array
//        ResMinFun(*, 0) = nu values
//ES added three Minkowski functionals more
//        ResMinFun(*, i) = Mink. functionals of order i-1, i=1..4.
{
   int Np,p,i,j,k;
//ES we need both the volume fraction vf and the normalized
//ES deviation from the mean nu.
   float ro,dd,meanro,sigmaro;
   float vf,mf=0;		// volume and mass fractions
   float prmass,m0,mtot;	//masses
   int N3,iro;
   float VFSTEP=MFSTEP;
   int Nx = dens.nx();
   int Ny = dens.ny();
   int Nz = dens.nz();
   N3=Nx*Ny*Nz;
//ES reserved space for padding
   intarray cells(Nx+1,Ny+1,Nz+1);	// Nonperiodic case needs padding by 0
   fltarray sdens(Nx,Ny,Nz);
//ES mass array
   fltarray smass(N3);
//ES Minkowski functionals
   fltarray MF(4);
//ES need to declare G here
   double G=0.;
   fltarray MinFun;
  
   sdens = dens;
   densort(sdens);
   
   Np = (int)(2/VFSTEP); 		// exact size not known
//   Np = (int)((1.-2. * VFSTEP)/(float)(VFSTEP) + 1); 
   if (Segmentation == False) MinFun.alloc(Np, 7);
   else MinFun.alloc(Np, 8);
   
   if (Verb == True)  cout << "# Np = " << Np << endl;

//ES meanro and sigmaro are needed to calculate nu
//ES We also set up a cumulative mass array
   for(i=0,meanro=sigmaro=prmass=0.;i<N3;i++) {
   	  dd=sdens(i);
	  meanro+=dd; sigmaro+=dd*dd;
	  prmass+=dd;		// sum the mass
	  smass(i)=prmass;
	  }
   meanro/=N3; sigmaro=sqrt(sigmaro/N3-meanro*meanro);
   if(Verb==True) printf("# meanro: %g  sigmaro: %g\n",meanro,sigmaro);
//ES density might be negative, total mass zero?
//ES suppose that the minimum density is 0, dens+= -sdens(0)
//ES these contortions map the mass interval to [0,1]
   m0=sdens(0);
   mtot=smass(N3-1)-N3*m0;
   for(i=0;i<N3;i++) 
		smass(i)=(smass(i)-(i+1)*m0)/mtot;
   
//ES The easiest way to pad
   for(i=0;i<=Nx;i++)
   for(j=0;j<=Ny;j++)
   for(k=0;k<=Nz;k++) cells(i,j,k) = 0;

  if (Verb == True) {
 	printf("#  mf\t  vf\t  ro\t  nu\t  MF0\t  MF1\t  MF2\t  MF3");
   if(Segmentation==True) printf("\t  G");
   printf("\n");
   fflush(stdout);
  } 
   p = 0;
//ES The loop is over the volume fraction vf, first
   for(vf=VFSTEP;vf<=1.-VFSTEP;vf+=VFSTEP) {
      iro=(int) (vf*N3);
      ro=sdens(iro);
	  mf=smass(iro);
      float NuXAbs = xerf(vf);
// NuXAbs is such that: Nu = \int_{-\infty}^NuXAb 1/\sqrt(2PI) exp{-x^2/2} dx
                  // NuXAbs corresponds to the x-axis in the 
			      // (Weinberg et all, 1997) ApJ paper.
      for(i=0;i<Nx;i++)
      for(j=0;j<Ny;j++)
      for(k=0;k<Nz;k++) cells(i,j,k) = (dens(i,j,k) >=ro) ? 1:0;
      minfun(cells,MF, Periodic); 
 //ES storing all for Minkowski functionals
      MinFun(p, 0) = MF(0);
      MinFun(p, 1) = MF(1);
      MinFun(p, 2) = MF(2);
      MinFun(p, 3) = MF(3);
      MinFun(p, 4) = vf;
      MinFun(p, 5) = ro;
      MinFun(p, 6) = NuXAbs;
    if (Segmentation == True) {
        G=genus_segmentation(cells, NuXAbs);
 		MinFun(p, 7) = G;
	}
    if (Verb == True) {
		printf("%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
		 mf,vf,ro,NuXAbs,MF(0),Step*Step*MF(1),Step*MF(2),MF(3));
    if (Segmentation == True) 
        printf("\t%.4g", G);
	 if (Verb == True) printf("\n");
	 fflush(stdout);
     }
	p++;
   }
   cout << "P = " << p << endl;
   
//   if(Verb==True) printf("# Switching over to mass fractions\n");
//ES Continue with loop over mass fraction
   for(mf=MFSTEP*(floor(mf/MFSTEP)+1);mf<=1.-MFSTEP;mf+=MFSTEP) {
	  for(i=0;i<N3;i++) if(smass(i)>=mf) break;
	  ro=sdens(i); 
	  vf=(float)i/N3;
      float NuXAbs = xerf(vf);
      for(i=0;i<Nx;i++)
      for(j=0;j<Ny;j++)
      for(k=0;k<Nz;k++) cells(i,j,k) = (dens(i,j,k) >=ro) ? 1:0;
      minfun(cells,MF, Periodic); 
 //ES storing all for Minkowski functionals
      MinFun(p, 0) = MF(0);
      MinFun(p, 1) = MF(1);
      MinFun(p, 2) = MF(2);
      MinFun(p, 3) = MF(3);
      MinFun(p, 4) = vf;
      MinFun(p, 5) = ro;
      MinFun(p, 6) = NuXAbs;
     if (Segmentation == True)
     {
        G=genus_segmentation(cells, NuXAbs);
 	    MinFun(p, 7) = G;
     }
      if (Verb == True) {
		printf("%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
		 mf,vf,ro,NuXAbs,MF(0),Step*Step*MF(1),Step*MF(2),MF(3));
     if (Segmentation == True)
        printf("\t%.4g", G);
	 printf("\n");
	 fflush(stdout);
     }
	p++;
   }
   if (Verb == True)  cout << "Number of points = " << p << endl;
   ResMinFun.alloc(p, 7);
   for (i = 0; i < p; i++)
   {
      ResMinFun(i, 0) = MinFun(i, 0);
      ResMinFun(i, 1) = MinFun(i, 1);
      ResMinFun(i, 2) = MinFun(i, 2);
      ResMinFun(i, 3) = MinFun(i, 3);
      ResMinFun(i, 4) = MinFun(i, 4);
      ResMinFun(i, 5) = MinFun(i, 5);
      ResMinFun(i, 6) = MinFun(i, 6); 
   }
}

/***************************************/

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

void make_theo_genus(fltarray &MinFun, double NuIndex, float Sigma)
{
   int i;
   double  Nu2;
   double Cst = 1./sqrt(2*PI);
   double L=1./Sigma*sqrt((NuIndex+3.)/6.);
   double Vlambda = L /sqrt(2*PI);
   double V2 = Vlambda*Vlambda;
   
   for (i = 0; i < MinFun.nx(); i++)
   {
      double nu = MinFun(i,4);
      Nu2 = nu*nu;
      MinFun(i,0) = (float) (0.5 - 0.5*erf(nu/sqrt(2.)));
      MinFun(i,1) = (float) (2./3.*Vlambda*Cst *exp(-Nu2/2));
      MinFun(i,2) = (float) (2./3.*V2*Cst* nu*exp(-Nu2/2));
      MinFun(i,3) = (float) (Vlambda*V2*Cst*(Nu2-1)*exp(-Nu2/2));
   }
}
  
 
/*********************************************************************/
 
#define ROM 10

void genpois(fltarray &dens)
// Generate a Poisson noise
{
   int i,j,k,n;
   double a;
   int Nx = dens.nx();
   int Ny = dens.ny();
   int Nz = dens.nz();    
   for(i=0;i<Nx;i++)
   for(j=0;j<Ny;j++)
   for(k=0;k<Nz;k++) 
   {
      for(n=0,a=1.;n<ROM;n++) a*=drand48();
      dens(i,j,k)=-log(a)/ROM;
   }
}

/***************************************************************************/
/*
void genvalgauss(float disp, float & x, float & y)    // Recipes 
{
    float v1,v2,r,fac;
    // double drand48();

    do {
        v1=2.0*drand48()-1; 
        v2=2.0*drand48()-1; 
        r=v1*v1+v2*v2;
        }
    while(r>=1.0);
    fac=sqrt(-2*disp*log((double) r)/r);
    x=v1*fac;
    y=v2*fac;
}
*/
/***************************************************************************/

void gengauss(fltarray &dens, double NuIndex, float POWA)
// Generate a randonm Gaussian field
{
   int i,j,k;
   float Sigma=1.;
   FFTN_3D FFT;
   FFT.CenterZeroFreq = True;

   // Create a Gaussian white noise
   simu_gaussian_noise3d (dens, Sigma);
   
   if (NuIndex > 0)
   {
      int Nx = dens.nx();
      int Ny = dens.ny();
      int Nz = dens.nz();
      int Nx2 = Nx/2;
      int Ny2 = Ny/2;
      int Nz2 = Nz/2;
   
      // Take the Fourier transform of the noise
      cfarray FD(Nx,Ny,Nz);
      FFT.fftn3d(dens, FD);
  
      // Weight each frequency following the power spectrum of
      // the randonm Gaussian field, P(nu) = nu^{-NuIndex}
      for (i=0; i < Nx; i++)
      for (j=0; j < Ny; j++)
      for (k=0; k < Nz; k++)
      {
          float u = (i -  Nx2) / (float) Nx;
          float v = (j -  Ny2) / (float) Ny;
          float w = (k -  Nz2) / (float) Nz;
          float Nu = sqrt(u*u + v*v + w*w);
	  float P = 0;
          if (Nu > 0) P = POWA*pow((double) Nu, NuIndex);
          FD(i,j,k) *= P;
      }
     
      // Take the inverse Fourier transform
      FFT.fftn3d(FD, True);
 
      // Keep the real part
      for (i=0; i < Nx; i++)
      for (j=0; j < Ny; j++)
      for (k=0; k < Nz; k++) dens(i,j,k)  = FD(i,j,k).real();
   }

}

/*********************************************************************/

//	read input data and test if it is a catalogue or a fits
void read_data(fltarray & Data, char *Name_Cube_In, float BinCat)
{
   
   if(strstr(Name_Cube_In,".fits") || strstr(Name_Cube_In,".Fits")
			|| strstr(Name_Cube_In,".FITS"))
   { 
       fits_read_fltarr(Name_Cube_In, Data);
       // io_3d_read_data(Name_Cube_In, Data);
   }
   else if (strstr(Name_Cube_In,".cat") || strstr(Name_Cube_In,".Cat")
			|| strstr(Name_Cube_In,".CAT"))
   {
      ArrayPoint Catalogue;//initialisation of catalogue
      Catalogue.read(Name_Cube_In);
      if(Catalogue.dim()!=3)
      {
         cout <<"Error: the catalogue do not have 3 dimensions ... "<<endl;
         exit(0);
      }
      int Nx=(int)((Catalogue.Pmax.x()-Catalogue.Pmin.x())/BinCat+0.5);
      int Ny=(int)((Catalogue.Pmax.y()-Catalogue.Pmin.y())/BinCat+0.5);
      int Nz=(int)((Catalogue.Pmax.z()-Catalogue.Pmin.z())/BinCat+0.5);
      Data.alloc(Nx,Ny,Nz);
      for(int i=0;i<Catalogue.np();i++)
      {
         int pos_x=(int)((Catalogue(i).x()-Catalogue.Pmin.x()) /BinCat);
         int pos_y=(int)((Catalogue(i).y()-Catalogue.Pmin.y()) /BinCat);
         int pos_z=(int)((Catalogue(i).z()-Catalogue.Pmin.z()) /BinCat);
         if(pos_x<0 || pos_x>=Nx || pos_y<0 || pos_y>=Ny || 
					pos_z<0 || pos_z>=Nz )
	 {
            cout <<"Error: Bad coordinate in the catalogue ... " << endl;
	    cout <<"       Cube size = " << Nx << "," << Ny << "," << Nz << endl;
	    cout <<"       Point number = " << i+1 << endl;
	    cout <<"        pos_x = " << pos_x << endl;
	    cout <<"        pos_y = " << pos_y << endl;
	    cout <<"        pos_z = " << pos_z << endl;
	    exit(-1);
	 }
	 Data(pos_x,pos_y,pos_z) +=1;
      }
   }
   else
   {
      cout <<"Error: unknown format file ..." << endl;
      exit(0);
   }
}

/*********************************************************************/
