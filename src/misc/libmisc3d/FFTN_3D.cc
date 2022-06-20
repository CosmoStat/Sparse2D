/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  20/08/01
**    
**    File:  FFTN_3D.h
**
*******************************************************************************
**
**    DESCRIPTION  3D FFT Class routine for any size of images
**    ----------- 
**
******************************************************************************/

#include "FFTN_3D.h"

/**********************************************************************/

void FFTN_3D::cb_shift(fltarray& Data1, fltarray& Data2, int Dx, int Dy, int Dz) {

   int Nx=Data1.nx();
   int Ny=Data1.ny();
   int Nz=Data1.nz();
    
   for (int i = 0; i < Nx; i++) 
   for (int j = 0; j < Ny; j++) 
   for (int k = 0; k < Nz; k++) 
   {
       int X = i + Dx;
       int Y = j + Dy;
       int Z = k + Dz;
       
       if (X < 0) X += Nx;
       if (X >= Nx) X -= Nx;
       if (Y < 0) Y += Ny;
       if (Y >= Ny) Y -= Ny;
       if (Z < 0) Z += Nz;
       if (Z >= Nz) Z -= Nz; 
             
       if ((X<0) || (X>=Nx)  || (Y<0) || (Y>=Ny)  || (Z<0) || (Z>=Nz)) Data2(i,j,k) = 0.0;
       else Data2(i,j,k) = Data1(X,Y,Z);
   }
}

/**********************************************************************/

void FFTN_3D::cb_shift(cfarray& Data1, cfarray& Data2, int Dx, int Dy, int Dz) {

   int Nx=Data1.nx();
   int Ny=Data1.ny();
   int Nz=Data1.nz();
    
   for (int i = 0; i < Nx; i++) 
   for (int j = 0; j < Ny; j++) 
   for (int k = 0; k < Nz; k++) 
   {
       int X = i + Dx;
       int Y = j + Dy;
       int Z = k + Dz;
       
       if (X < 0) X += Nx;
       if (X >= Nx) X -= Nx;
       if (Y < 0) Y += Ny;
       if (Y >= Ny) Y -= Ny;
       if (Z < 0) Z += Nz;
       if (Z >= Nz) Z -= Nz; 
             
       if ((X < 0) || (X >= Nx)  || (Y < 0) || (Y >= Ny) || (Z < 0) || (Z >= Nz)) Data2(i,j,k) = complex_f(0.0,0.0);
       else Data2(i,j,k) = Data1(X,Y,Z);
   }
}

/**********************************************************************/

void FFTN_3D::center(fltarray &Cube) {

   int Nx=Cube.nx();
   int Ny=Cube.ny();
   int Nz=Cube.nz();
 
   fltarray Dat(Nx,Ny,Nz,(char*)"dat");
   Dat = Cube;
   cb_shift(Dat, Cube, (Nx+1)/2, (Ny+1)/2, (Nz+1)/2);
}
/**********************************************************************/

void FFTN_3D::uncenter(fltarray &Cube) {

   int Nx=Cube.nx();
   int Ny=Cube.ny();
   int Nz=Cube.nz();
   fltarray Dat(Nx,Ny,Nz,(char*)"dat");
   Dat = Cube;
   cb_shift(Dat, Cube, -(Nx+1)/2, -(Ny+1)/2, -(Nz+1)/2);
}
/**********************************************************************/

void FFTN_3D::center(cfarray &Cube) {

   int Nx=Cube.nx();
   int Ny=Cube.ny();
   int Nz=Cube.nz();
 
   cfarray Dat(Nx,Ny,Nz,(char*)"dat");
   Dat = Cube; 
   cb_shift(Dat, Cube, (Nx+1)/2, (Ny+1)/2, (Nz+1)/2);
}
/**********************************************************************/

void FFTN_3D::uncenter(cfarray &Cube) {

   int Nx=Cube.nx();
   int Ny=Cube.ny();
   int Nz=Cube.nz();
   cfarray Dat(Nx,Ny,Nz,(char*)"dat");
   Dat = Cube;
   cb_shift(Dat, Cube, -(Nx+1)/2, -(Ny+1)/2, -(Nz+1)/2);
}
/**********************************************************************/







/******************************************************************************/
void FFTN_3D::swap_buff(complex_f *Buff, int Nx, int Ny, int Nz, Bool Reverse)
{
    int i,j,k,Indi,Indj,Indk;
    complex_f temp;

    if (Reverse == False)
    {
       if (Nx % 2 != 0) 
       {
          for (j = 0; j < Ny; j++)
	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[Nx/2+j*Nx+k*Nx*Ny];
             for (i = Nx/2; i < Nx-1; i++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+1+j*Nx+k*Nx*Ny];
             Buff[Nx-1+j*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Ny % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
 	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[i+ (Ny/2)*Nx+k*Nx*Ny];
             for (j = Ny/2; j < Ny-1; j++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+(j+1)*Nx+k*Nx*Ny];
             Buff[i+(Ny-1)*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Nz % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
	  for (j = 0; j < Ny; j++)
	  {
             temp = Buff[i+j*Nx+(Nz/2)*Nx*Ny];
             for (k = Nz/2; k < Nz-1; k++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+j*Nx+(k+1)*Nx*Ny];
             Buff[i+j*Nx+(Nz-1)*Nx*Ny] = temp;
	  }
       }
    }

 
    for (i = 0; i < Nx/2; i++)
    for (j = 0; j < Ny/2; j++)
    for (k = 0; k < Nz/2; k++)
    {
       Indi = Nx/2+i;
       Indj = Ny/2+j;
       Indk = Nz/2+k;
       pix_swap(Buff[i+j*Nx+k*Nx*Ny], Buff[Indi+Indj*Nx+Indk*Nx*Ny]);
       pix_swap(Buff[Indi+j*Nx+k*Nx*Ny], Buff[i+Indj*Nx+Indk*Nx*Ny]);
       pix_swap(Buff[i+j*Nx+Indk*Nx*Ny], Buff[Indi+Indj*Nx+k*Nx*Ny]);
       pix_swap(Buff[Indi+j*Nx+Indk*Nx*Ny], Buff[i+Indj*Nx+k*Nx*Ny]);
    }  


    if (Reverse == True)
    {
       if (Nx % 2 != 0) 
       {
          for (j = 0; j < Ny; j++)
	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[Nx-1+j*Nx+k*Nx*Ny];
             for (i = Nx-1; i > Nx/2; i--) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i-1+j*Nx+k*Nx*Ny];
             Buff[Nx/2+j*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Ny % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
 	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[i+ (Ny-1)*Nx+k*Nx*Ny];
	     for (j = Ny-1; j > Ny/2; j--) 
 	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+(j-1)*Nx+k*Nx*Ny];
             Buff[i+(Ny/2)*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Nz % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
	  for (j = 0; j < Ny; j++)
	  {
             temp = Buff[i+j*Nx+(Nz-1)*Nx*Ny];
 	     for (k = Nz-1; k > Nz/2; k--)
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+j*Nx+(k-1)*Nx*Ny];
             Buff[i+j*Nx+(Nz/2)*Nx*Ny] = temp;
	  }
       } 
       
    }
}
/******************************************************************************/

void FFTN_3D::swap_buff(complex_d *Buff, int Nx, int Ny, int Nz, Bool Reverse)
{
    int i,j,k,Indi,Indj,Indk;
    complex_d temp;

    if (Reverse == False)
    {
       if (Nx % 2 != 0) 
       {
          for (j = 0; j < Ny; j++)
	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[Nx/2+j*Nx+k*Nx*Ny];
             for (i = Nx/2; i < Nx-1; i++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+1+j*Nx+k*Nx*Ny];
             Buff[Nx-1+j*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Ny % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
 	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[i+ (Ny/2)*Nx+k*Nx*Ny];
             for (j = Ny/2; j < Ny-1; j++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+(j+1)*Nx+k*Nx*Ny];
             Buff[i+(Ny-1)*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Nz % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
	  for (j = 0; j < Ny; j++)
	  {
             temp = Buff[i+j*Nx+(Nz/2)*Nx*Ny];
             for (k = Nz/2; k < Nz-1; k++) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+j*Nx+(k+1)*Nx*Ny];
             Buff[i+j*Nx+(Nz-1)*Nx*Ny] = temp;
	  }
       }
    }

 
    for (i = 0; i < Nx/2; i++)
    for (j = 0; j < Ny/2; j++)
    for (k = 0; k < Nz/2; k++)
    {
       Indi = Nx/2+i;
       Indj = Ny/2+j;
       Indk = Nz/2+k;
       pix_swap(Buff[i+j*Nx+k*Nx*Ny], Buff[Indi+Indj*Nx+Indk*Nx*Ny]);
       pix_swap(Buff[Indi+j*Nx+k*Nx*Ny], Buff[i+Indj*Nx+Indk*Nx*Ny]);
       pix_swap(Buff[i+j*Nx+Indk*Nx*Ny], Buff[Indi+Indj*Nx+k*Nx*Ny]);
       pix_swap(Buff[Indi+j*Nx+Indk*Nx*Ny], Buff[i+Indj*Nx+k*Nx*Ny]);
    }  


    if (Reverse == True)
    {
       if (Nx % 2 != 0) 
       {
          for (j = 0; j < Ny; j++)
	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[Nx-1+j*Nx+k*Nx*Ny];
             for (i = Nx-1; i > Nx/2; i--) 
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i-1+j*Nx+k*Nx*Ny];
             Buff[Nx/2+j*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Ny % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
 	  for (k = 0; k < Nz; k++)
	  {
             temp = Buff[i+ (Ny-1)*Nx+k*Nx*Ny];
	     for (j = Ny-1; j > Ny/2; j--) 
 	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+(j-1)*Nx+k*Nx*Ny];
             Buff[i+(Ny/2)*Nx+k*Nx*Ny] = temp;
	  }
       }
       if (Nz % 2 != 0) 
       {
          for (i = 0; i < Nx; i++)
	  for (j = 0; j < Ny; j++)
	  {
             temp = Buff[i+j*Nx+(Nz-1)*Nx*Ny];
 	     for (k = Nz-1; k > Nz/2; k--)
	       Buff[i+j*Nx+k*Nx*Ny] = Buff[i+j*Nx+(k-1)*Nx*Ny];
             Buff[i+j*Nx+(Nz/2)*Nx*Ny] = temp;
	  }
       } 
       
    }
}




/******************************************************************************/


void FFTN_3D::fftn3d(fltarray &Data, complex_f *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.n_elem();
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();

    if ((Reverse == True) && (CenterZeroFreq == True))
                                   swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
    float *FFT_Dat = (float *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
       //  transform1d return the FFT whith the zero frequency at the left
    if ((Reverse == False) && (CenterZeroFreq == True)) 
                                      swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
}

 
/******************************************************************************/

void FFTN_3D::fftn3d(fltarray &Data, cfarray &Data_Out, Bool Reverse, bool normalize) {

//cerr<<"FFT3D : data.nx="<<Data.nx()<<", fft.nx="<<Data_Out.nx()<<endl;
    int Ind=0;
    int N = Data.n_elem();
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();

    float *FFT_Dat = (float *) Data_Out.buffer();
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }       
    if (CenterZeroFreq == True) uncenter(Data_Out);
    transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
    if (CenterZeroFreq == True) center(Data_Out);
}

/******************************************************************************/

void FFTN_3D::fftn3d (cfarray &Buff,  Bool Reverse, bool normalize) {
  
    int Ind=0;
    int N = Buff.n_elem();
    int Nx = Buff.nx();
    int Ny = Buff.ny();
    int Nz = Buff.nz();

    
    float *FFT_Dat = (float *) Buff.buffer();
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = (Buff(i)).real();
       FFT_Dat[Ind++] = (Buff(i)).imag();
    }       
    
    if (CenterZeroFreq == True) uncenter(Buff);
    transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
    if (CenterZeroFreq == True) center(Buff);
    
   
}
/******************************************************************************/

void FFTN_3D::fftn3d(complex_f *Dat, int Nx, int Ny, int Nz, Bool Reverse, bool normalize)
{
   float *FFT_Dat = (float *) Dat;
   if ((Reverse == True) && (CenterZeroFreq == True))
                                             swap_buff(Dat, Nx, Ny, Nz, Reverse);
   transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
   if ((Reverse == False) && (CenterZeroFreq == True))
                                             swap_buff(Dat, Nx, Ny, Nz, Reverse);
}

/******************************************************************************/
// Ditto for double
/******************************************************************************/

void FFTN_3D::fftn3d(fltarray &Data, complex_d *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.n_elem();
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();

    if ((Reverse == True) && (CenterZeroFreq == True))
                                   swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
    double *FFT_Dat = (double *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
       //  transform1d return the FFT whith the zero frequency at the left
    if ((Reverse == False) && (CenterZeroFreq == True)) 
                                      swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
}

/******************************************************************************/

void FFTN_3D::fftn3d(dblarray &Data, complex_d *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.n_elem();
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();

    if ((Reverse == True) && (CenterZeroFreq == True))
                                   swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
    double *FFT_Dat = (double *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
       //  transform1d return the FFT whith the zero frequency at the left
    if ((Reverse == False) && (CenterZeroFreq == True)) 
                                      swap_buff(Data_Out, Nx, Ny, Nz, Reverse);
}
/******************************************************************************/

void FFTN_3D::fftn3d(fltarray &Data, cdarray &Data_Out, Bool Reverse, bool normalize)
{
    fftn3d(Data,Data_Out.buffer(),Reverse, normalize);
}

/******************************************************************************/

void FFTN_3D::fftn3d(dblarray &Data, cdarray &Data_Out, Bool Reverse, bool normalize)
{
    fftn3d(Data,Data_Out.buffer(),Reverse, normalize);
}    
/******************************************************************************/

void FFTN_3D::fftn3d (cdarray &Buff,  Bool Reverse, bool normalize)
 {
    fftn3d(Buff.buffer(), Buff.nx(), Buff.ny(), Buff.nz(), Reverse, normalize);
 }
/******************************************************************************/

void FFTN_3D::fftn3d(complex_d *Dat, int Nx, int Ny, int Nz, Bool Reverse, bool normalize)
{
   double *FFT_Dat = (double *) Dat;
   if ((Reverse == True) && (CenterZeroFreq == True))
                                             swap_buff(Dat, Nx, Ny, Nz, Reverse);
   transform3d(FFT_Dat, Nx, Ny, Nz, Reverse, normalize);
   if ((Reverse == False) && (CenterZeroFreq == True))
                                             swap_buff(Dat, Nx, Ny, Nz, Reverse);
}

/******************************************************************************/
 
void FFTN_3D::convolve(fltarray & Data, fltarray & Gauss)
// Convolve two cubes: Data = Data * Gauss
{
   int i,j,k;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   cfarray FD(Nx,Ny,Nz);
   cfarray FG(Nx,Ny,Nz);
   fftn3d(Data, FD);
   fftn3d(Gauss, FG);
   FD *= FG;
   fftn3d(FD, True);
   for (i=0; i < Nx; i++)
   for (j=0; j < Ny; j++)
   for (k=0; k < Nz; k++) Data(i,j,k)  = FD(i,j,k).real();
}

/******************************************************************************/
 
void convolve3d(fltarray & Data, fltarray & Gauss)
// Convolve two cubes: Data = Data * Gauss
{
   FFTN_3D FFT;
   FFT.CenterZeroFreq = True;
   int i,j,k;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   cfarray FD(Nx,Ny,Nz);
   cfarray FG(Nx,Ny,Nz);
   FFT.fftn3d(Data, FD);
   FFT.fftn3d(Gauss, FG);
   FD *= FG;
   FFT.fftn3d(FD, True);
   for (i=0; i < Nx; i++)
   for (j=0; j < Ny; j++)
   for (k=0; k < Nz; k++) Data(i,j,k)  = FD(i,j,k).real();
}

/***************************************************************************/
