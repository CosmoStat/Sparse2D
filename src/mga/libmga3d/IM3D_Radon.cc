/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  15/10/2001 
**    
**    File:  IM_Radon3D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Radon3D transform and reconstruction
**    -----------  
**                 
******************************************************************************/

#include "IM3D_Radon.h"
#include "IM_IO.h"
#include "FFTN_1D.h"

/***************************************/

void Radon3D::alloc(int Nx_Dat, int Ny_Dat, int Nz_Dat, type_radon3d RM)
{
   Alloc = True;
   RadonMethod=RM; 
   Nx = Nx_Dat;
   Ny = Ny_Dat;
   Nz = Nz_Dat;
   CFFT3D.CenterZeroFreq = True;
   
   if (RadonMethod == RADON3D_FFT)
   {
      ProjectNbr = Nx*Nx*3;
      ProjectNbr2D = 2*Nx;
      ResolNbr = Nx;
   }
   // NormRadonFFT = (float) ResolNbr / (float) (Nl*Nc);

   if (Verbose == True)
   {
      cout << "Radon3D image size: Nl = " << ProjectNbr << " Nc = " << ResolNbr << endl;
   }
 
   make_table();
}

/***************************************/

void Radon3D::reset_param() 
{
   Alloc=False; 
   ProjectNbr = 0;
   ProjectNbr2D = 0;
   Verbose=False;
   ResolNbr=0;
   RadonMethod=DEF_RADON3D;
   NormRadonFFT=1.;
   if (SinTab != NULL) 
   {
        delete [] SinTab;
        SinTab = NULL;
   }
   if (CosTab != NULL) 
   {
      delete [] SinTab;
      CosTab = NULL;
   }
   if (TabAngle  != NULL) 
   {
      delete []  TabAngle;
      TabAngle = NULL;
   }
   if (TabAngle2  != NULL) 
   {
      delete []  TabAngle2;
      TabAngle2 = NULL;
   }
   if (TabStep != NULL) 
   {
      delete [] TabStep;
      TabStep = NULL;
   }
}

/***************************************/

// void Radon3D::num_project(int p, int & u, int & v, int & w, float &Angle1, float &Angle2)
// // number of projects: 2N^2 * (N/2) + N^2 = N^3+
// {
//    int Ind;
//    double L,D = Nx/2;
//    u = v = w = 0;
//    
//    if (p < 2*Nx*Nx)
//    { 
//       int Ind2D = p % ProjectNbr2D;  // ProjectNbr2D=2N
//       if (Ind2D <= Nx / 2)
//       {
//          u = Nx/2;
//          v = Ind2D;
// 	 Angle1 = atan2( (double)(v), u);
//       }
//       else if (Ind2D < Nx+Nx/2)
//       {
//          u = Nx - Ind2D;
//          v = Nx/2;
// 	 Angle1 = atan2( (double)(v), u);
//       }
//       else 
//       {
//           u = -Nx/2;
//           v = 2*Nx-Ind2D;
// 	  Angle1 =  atan2( (double)(v) , u);
//       }
//       w = p / ProjectNbr2D - Nx/2;
//       L = sqrt((double)(u)*u+(double)(v)*v);
//       Angle2 = atan2( (double)(w) , L);
//    }
//    else
//    {
//       w = Nx/2;
//       Ind = p - 2*Nx*Nx;
//       v = Ind / Nx - Nx/2;
//       u = Ind % Nx - Nx/2;
//       // if (u == 0) Angle1 = PI/2;
//       Angle1 =  atan2( (double)(v), (double)(u));
//       L = sqrt((double)(u)*u+(double)(v)*v);
//       // if (L == 0) Angle2 =PI/2.;
//       Angle2 = atan2( (double)(w), L);
//       
//    } 
//    if ((u == -2) && (v == -1) && (w == 4)) cout << "OK"  << endl; 
//    if ((u == 4) && (v == 0)&& (w == -4))
//         cout << "INFO L = " << L << " " << atan2( (double)(w) , L) << endl;
// }

/***************************************/

void Radon3D::make_table()
{
   int t,p;
   //float Angle1, Angle2;
   TabAngle = new double[ProjectNbr];
   TabAngle2 = new double[ProjectNbr];
   TabStep = new double[ProjectNbr];
   fltarray DD(Nx,Nx,Nx);
   double D = Nx/2;
   

   for (t=0; t <= Nx/2; t++) 
   {
       double A = atan( (double) t / D);
       TabAngle[t] = A;
       TabAngle[Nx-t] = PI/2. - A;
       TabAngle[Nx+t] = PI/2 + A;
       if (t > 0) TabAngle[ProjectNbr2D -t] = PI - A;
   }
   
   for (p=0; p < ProjectNbr; p++)
   {
      int u,v,w;
      double L;
      if (p < 2*Nx*Nx)
      { 
	 int Pos  = (p / ProjectNbr2D) * ProjectNbr2D;
         t = p % ProjectNbr2D;  // ProjectNbr2D=2N
 	 if (Pos != 0) TabAngle[Pos+t] = TabAngle[t];
         if (t <= Nx / 2)
         {
            u = Nx/2;
            v = t;
         }
         else if (t <= Nx+Nx/2)
         {
            u = Nx - t;
            v = Nx/2;
         }
         else 
         {
             u = -Nx/2;
             v = 2*Nx-t;
         }         
	 w = p / ProjectNbr2D - Nx/2;
         L = sqrt((double)(u)*u+(double)(v)*v+(double)(w)*w);
         TabStep[p] = 2.*L / (double) (ResolNbr);
      }
      else
      {
         w = Nx/2;
         t = p - 2*Nx*Nx;
         v = t / Nx - Nx/2;
         u = t % Nx - Nx/2;
         TabAngle[p] =  atan2( (double)(v), (double)(u));
	 L = sqrt((double)(u)*u+(double)(v)*v+(double)(w)*w);
   	 TabStep[p] = (2.*L) / (double) (ResolNbr-1);
      }

      L = sqrt((double)(u)*u+(double)(v)*v);
      TabAngle2[p] =  atan2( (double)(w) , L);
  
//      if (Verbose == True)
//         cout << "Proj num = "  << p+1 << " " << u  << " " << v  << " " << w  <<  "  Angles = "  << TabAngle[p]*180./PI << " " << TabAngle2[p]*180./PI  << " " <<  TabStep[p] << endl;

   } 
}

/*************************************************************************/
 
void Radon3D::transform(fltarray &Cube, Ifloat &Radon3DImage)
{
   // If the class is not allocated, 
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) || (Nx != Cube.nx()) 
        || (Ny != Cube.ny()) || (Nz != Cube.nz())) 
      alloc(Cube.nx(), Cube.ny(), Cube.nz());
 
   if ((Radon3DImage.nl() != ProjectNbr) || ( Radon3DImage.nc() != ResolNbr)) 
                                        Radon3DImage.resize(ProjectNbr,ResolNbr);
   Radon3DImage.init();
   if (Verbose == True)
   {
       cout << "RADON transform: input cube size = " << Nx << " " << Ny << " " << Nz << endl;
       cout << "RADON transform: output image size = " << Radon3DImage.nl() << " " << Radon3DImage.nc() << endl;
   }

   switch (RadonMethod)
   {
       case RADON3D_FFT:
              fft_trans(Cube, Radon3DImage);
              break;
   }
   
//   io_write_ima_float((char*)"RadT.fits", Radon3DImage);
}


/***************************************/

void Radon3D::fft_trans(fltarray &Cube, Ifloat &Radon3DImage)
{
   cfarray Data2(ResolNbr,ProjectNbr);
   fft_trans_cf(Cube, Data2);
   invfft_line(Data2,Radon3DImage);
}

/***************************************/

/*static int iiround(double point)
{
   int result;
   if (point >= (double) 0.0) result = (int) ( point + (double) 0.51);
   else result = (int) (point - 0.5);
   // printf ("iiround: %F ==> %F, %d \n", point, point+0.5, result);
   return result;
}*/

/***************************************/
//extern int numero_du_block;
void Radon3D::fft_trans_cf(fltarray &Cube, cfarray & Data2)
//  Cube: IN = input cube 
//  Data2: OUT = 1D Fourier transform of the Radon lines
//          FFT^-1(Data2(i,*)) = ith line of the Radon transform
{
	int t,r,x,y,z,x1,y1,z1;
	double u,v,w;
	if ((Alloc == False) || (Nx != Cube.nx()) || (Ny != Cube.ny()) || (Nz != Cube.nz()))
		alloc(Cube.nx(), Cube.ny(), Cube.nz(), RADON3D_FFT);
	cfarray Data1(Nx,Ny,Nz);
	int NlRadon3D = ProjectNbr;
	int NcRadon3D = ResolNbr ;
	// Cube.info("CubeT");
	if ((Data2.ny() != ProjectNbr) || ( Data2.nx() != ResolNbr)) 
		Data2.resize(ResolNbr, ProjectNbr);   
	if ((Data2.ny() != NlRadon3D) || (Data2.nx() != NcRadon3D))
	{
		cout << "Error: bad image size for the Radon3D transform ..." << endl;
		exit(-1);
	}
	//cout << "MinDat = " << Cube.min() << endl;
	//cout << "MaxDat = " << Cube.max() << endl;
   
   
FFTN_3D CFFT3D_local;

	CFFT3D_local.CenterZeroFreq = True;

//char filename[64];
//sprintf(filename,"Cube_%d.fits",numero_du_block);
//fits_write_fltarr(filename, Cube);
// identiques

	CFFT3D_local.fftn3d(Cube, Data1);

//fltarray prov(Nx,Ny,Nz);
//for (int i=0;i<Nx;i++)
//for (int j=0;j<Ny;j++)
//for (int k=0;k<Nz;k++)
//   prov(i,j,k)=Data1(i,j,k).real();
//sprintf(filename,"Data1_%d.fits",numero_du_block);
//fits_write_fltarr(filename, prov);
// data différents

	for (t=0; t < NlRadon3D; t++)
	for (r=0; r < NcRadon3D; r++) Data2(r,t) = complex_f(0.,0.);
	fltarray DD(Nx,Nx,Nx);

//for (x=0; x< Nx; x++)
//for (y=0; y< Nx; y++)
//for (z=0; z< Nx; z++) DD(x,y,z) = Data1(x,y,z).real();
//DD.info("real rec");
//fits_write_fltarr("xx_ore.fits", DD);
//for (x=0; x< Nx; x++)
//for (y=0; y< Nx; y++)
//for (z=0; z< Nx; z++) DD(x,y,z) = Data1(x,y,z).imag();
//DD.info("imag rec");
//fits_write_fltarr("xx_oim.fits", DD);
   
      
	for (t=0; t < NlRadon3D; t++)
	{
		double Angle  = TabAngle[t];
		double Angle2 = TabAngle2[t];
		double StepProj = TabStep[t];
		int rsym,DirLine = 1;
		r = 0;
		double rr = DirLine*((double) r-NcRadon3D/2.) * StepProj;
		u = rr * cos(Angle)* cos(Angle2);
		v = rr * sin(Angle)* cos(Angle2);
		w = rr * sin(Angle2);
		if (t >= 2*Nx*Nx)
		{
			u = -u;
			v = -v;
		}

		if (NcRadon3D % 2 == 0)
		{
			if ((u >= Nx/2) && (v < Nx/2) && (w < Nx/2)) DirLine = -1;
			else if ((u < Nx/2) && (v >= Nx/2) && (w < Nx/2)) DirLine = -1;
			else if ((u < Nx/2) && (v < Nx/2) && (w >= Nx/2)) DirLine = -1;
		}


		for (r=0; r <= NcRadon3D/2; r++)
		{
			rr = DirLine * ((double) r-NcRadon3D/2) * StepProj;
			w = rr * sin(Angle2);
			u = rr * cos(Angle)* cos(Angle2);
			v = rr * sin(Angle)* cos(Angle2);
			if ((NcRadon3D % 2 == 0) && (t >= 2*Nx*Nx))
			{
				u = -u;
				v = -v;

				int LastNx = (Nx % 2 == 0) ? Nx/2 - 1: Nx/2;
				if (u > LastNx) u = LastNx;
				else if (u < -Nx/2) u = -Nx/2;
				if (v > LastNx) v = LastNx;
				else if (v < -Nx/2) v = -Nx/2;
				if (w > LastNx) w = LastNx;
				else if (w < -Nx/2) w = -Nx/2;
			}
			x1 = iround(u);
			y1 = iround(v);
			z1 = iround(w);
			x = x1 + Nx/2;
			y = y1 + Ny/2;
			z = z1 + Nz/2;

			// if (r == 0)
			//  printf("Proj %d: %1d %1d %1d , A2 = %5.2f, %5.2f %5.2f %5.2f \n", t, x,y,z, Angle2*180./ PI, u,v,w);

			if ((y>=0) && (y< Ny) && (x>=0) && (x<Nx) 
				&& (z< Nz)  && (z>=0))  
			{
				Data2(r,t) =  Data1(x,y,z);
				DD(x,y,z) += 1.;      
			}
			int xs = -x1;
			int ys = -y1;
			int zs = -z1;
			x = xs + Nx/2;
			y = ys + Ny/2;
			z = zs + Nz/2;
			rsym = 2*((int)(NcRadon3D/2)) - r;
			if ((y>=0) && (y< Ny) && (x>=0) && (x<Nx) 
				&& (z>=0) && (z<Nz) && (rsym >= 0) && 
				(rsym < NcRadon3D))
			{
				// if (DataCount1(y,x) == 0)
				{
					Data2(rsym,t) = Data1(x,y,z); 
					DD(x,y,z)+=1;  
				}
			}
			else if ((rsym >= 0) && (rsym <  NcRadon3D))
			{
				cout << "Error occured during the transformation ..." << endl;
				exit(-1);
			}
		}
	}
// Data2(NcRadon3D/2, NlRadon3D/2) = Data1(Nx/2, Ny/2, Nz/2);
// DD.info("DD");
// fits_write_fltarr ("xx_dd", DD);

//cout << "DD MIN = " << DD.min() << endl;
// io_3d_write_data("xx_rad", DD);
//     Ifloat IMRAD(NlRadon3D, NcRadon3D,"RAD");
//     for (t=0; t < NlRadon3D; t++)
//     for (r=0; r < NcRadon3D; r++)
//        IMRAD(t,r) = Data2(r,t).real();
//      printf("IMRAD  MinRe = %f, MaxRe =  %f\n", min(IMRAD) , max(IMRAD));
//     io_write_ima_float("xx_real.fits", IMRAD);
//     for (t=0; t < NlRadon3D; t++)
//     for (r=0; r < NcRadon3D; r++)
//        IMRAD(t,r) = Data2(r,t).imag();
//      printf("       MinImag = %f, MaxImag =  %f\n", min(IMRAD) , max(IMRAD));
//     io_write_ima_float("xx_imag.fits", IMRAD);
}

/***************************************/

void Radon3D::invfft_line(cfarray & Data2, Ifloat &Radon3DImage)
{
   int i,j;
   int NlRadon3D = ProjectNbr;
   int NcRadon3D = ResolNbr ;
   //float Norm = NormRadonFFT;  
    Ifloat IMRAD(NlRadon3D, NcRadon3D,(char*)"RAD");

   if ((Radon3DImage.nl() != Data2.ny()) || (Radon3DImage.nc() != Data2.nx()))
      	                              Radon3DImage.resize(NlRadon3D,NcRadon3D);
					       
   FFTN_1D FFT1D;
   FFT1D.CenterZeroFreq = True;
   cdarray Line(NcRadon3D);
   for (i=0; i< NlRadon3D; i++)
   {
      for (j=0; j < NcRadon3D; j++)  Line(j) =  Data2(j,i);
      FFT1D.fftn1d (Line.buffer(), NcRadon3D, True);
      for (j=0; j < NcRadon3D; j++) Radon3DImage(i,j) = (float) (Line(j)).real();
      for (j=0; j < NcRadon3D; j++) IMRAD(i,j) = (float) (Line(j)).imag();
   } 
//   INFO(Radon3DImage, (char*)"RAD");
//   INFO(Radon3DImage, (char*)"RAD IMAG");
}


/***************************************/
//       Inverse Radon Transform
/***************************************/


void Radon3D::recons(Ifloat &Radon3DImage, fltarray &CubeRec)
{
   // INFO(Radon3DImage, "RAD");
   // If the class is not allocated,  
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) || 
        (ProjectNbr  !=  Radon3DImage.nl()) || (ResolNbr != Radon3DImage.nc()))
   {
       int Nr =  Radon3DImage.nc();
       alloc(Nr,Nr,Nr,DEF_RADON3D);
   }
   if (( CubeRec.nx() != Nx) || (CubeRec.ny() != Ny)
     || (CubeRec.nz() != Nz)) CubeRec.reform(Nx,Ny,Nz);
   CubeRec.init();
   if (Verbose == True)
   {
       cout << "RADON reconstruct: input image size = " <<  Radon3DImage.nl() << " " <<  Radon3DImage.nc() << endl;
       cout << "RADON reconstruct: output cube size = " << CubeRec.nx() << " " <<  
                   CubeRec.ny() <<  " " <<  CubeRec.nz() <<endl;
   }
   switch (RadonMethod)
   {
        case RADON3D_FFT:
              fft_rec(Radon3DImage, CubeRec);
              break;
   }
}
  
/***************************************/

void Radon3D::fft_rec(Ifloat &Radon3DImage, fltarray &Cube)
{
   int NlRadon3D = ProjectNbr;
   int NcRadon3D = ResolNbr ;
   cfarray Data1(NcRadon3D, NlRadon3D);

   fft_line(Radon3DImage, Data1);
   fft_rec_cf(Data1, Cube);
}
 
/***************************************/

void Radon3D::fft_line(Ifloat &Radon3DImage, cfarray & Data1)
{
   int i,j;
   int NlRadon3D = ProjectNbr;
   int NcRadon3D = ResolNbr ;
 
   FFTN_1D FFT1D;
   FFT1D.CenterZeroFreq = True;
   cdarray Line(NcRadon3D);
   for (i=0; i< NlRadon3D; i++)
   {
      for (j=0; j < NcRadon3D; j++) Line(j)=complex_d((double) Radon3DImage(i,j), 0);
      FFT1D.fftn1d(Line.buffer(), NcRadon3D, False);
      for (j=0; j < NcRadon3D; j++)  
          Data1(j,i) =  complex_f(Line(j).real(), Line(j).imag());
   }    
}

/***************************************/

void Radon3D::fft_rec_cf(cfarray & Data1, fltarray &Cube)
{
   int t,r,x,y,z,x1,y1,z1;
   double u,v,w;
   // If the class is not allocated,  
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) || 
        (ProjectNbr  !=  Data1.ny()) || (ResolNbr != Data1.nx()))
   {
       int Nr =  Data1.nx();
       alloc(Nr,Nr,Nr,RADON3D_FFT);
   }
   int NlRadon3D = ProjectNbr;
   int NcRadon3D = ResolNbr ;
   int SizeFT = NcRadon3D;   // sampling size in Fourier space
          
   if (( Cube.nx() != Nx) || ( Cube.ny() != Ny) || ( Cube.nz() != Nz))
       Cube.reform(Nx,Ny,Nz);
   Cube.init();
   if (Verbose == True)
   {
       cout << "RADON reconstruct: input image size = " <<  Data1.ny() << " " <<  Data1.nx() << endl;
       cout << "RADON reconstruct: output cube size = " << Cube.nx() << " " <<  Cube.ny() << " " <<  Cube.nz() <<endl;
   }
   
   cfarray Data2(SizeFT, SizeFT, SizeFT);
   fltarray DataCount(SizeFT, SizeFT, SizeFT);
   
   for (t=0; t < NlRadon3D; t++)
   {
       double Angle  = TabAngle[t];
       double Angle2 = TabAngle2[t];
       double StepProj = TabStep[t];
       int rsym,DirLine = 1;
       r = 0;
       double rr = DirLine*((double) r-NcRadon3D/2.) * StepProj;
       u = rr * cos(Angle)* cos(Angle2);
       v = rr * sin(Angle)* cos(Angle2);
       w = rr * sin(Angle2);
       if (t >= 2*Nx*Nx)
       {
	   u = -u;
           v = -v;
       }
       if (NcRadon3D % 2 == 0)
       {
          if ((u >= Nx/2) && (v < Nx/2) && (w < Nx/2)) DirLine = -1;
          else if ((u < Nx/2) && (v >= Nx/2) && (w < Nx/2)) DirLine = -1;
          else if ((u < Nx/2) && (v < Nx/2) && (w >= Nx/2)) DirLine = -1;
       }
       for (r=0; r <= NcRadon3D/2; r++)
       {
          rr = DirLine * ((double) r-NcRadon3D/2) * StepProj;
          w = rr * sin(Angle2);
 	  u = rr * cos(Angle)* cos(Angle2);
	  v = rr * sin(Angle)* cos(Angle2);
	  if ((NcRadon3D % 2 == 0) && (t >= 2*Nx*Nx))
	  {
	      u = -u;
              v = -v;
 	      int LastNx = (Nx % 2 == 0) ? Nx/2 - 1: Nx/2;
              if (u > LastNx) u = LastNx;
	      else if (u < -Nx/2) u = -Nx/2;
  	      if (v > LastNx) v = LastNx;
	      else if (v < -Nx/2) v = -Nx/2;
  	      if (w > LastNx) w = LastNx;
	      else if (w < -Nx/2) w = -Nx/2;
	 }
         x1 = iround(u);
         y1 = iround(v);
         z1 = iround(w);
         x = x1 + Nx/2;
         y = y1 + Ny/2;
         z = z1 + Nz/2;  
          
	  
         if ((y>=0) && (y< Ny) && (x>=0) && (x<Nx) && (z>=0) && (z<Nz)) 
         {
            Data2(x,y,z) +=  Data1(r,t); 
            DataCount(x,y,z) += 1.;
         }
	int xs = -x1;
	int ys = -y1;
	int zs = -z1;
	x = xs + Nx/2;
        y = ys + Ny/2;
        z = zs + Nz/2;
	rsym = 2*((int)(NcRadon3D/2)) - r;
	if ((y>=0) && (y< Ny) && (x>=0) && (x<Nx) 
	        && (z>=0) && (z<Nz) && (rsym >= 0) && 
	         (rsym < NcRadon3D))  
	{
	    Data2(x,y,z) +=  Data1(rsym,t);
	    DataCount(x,y,z) += 1.;
	}
      }
   }
//   DataCount.info("DDC");
//   fits_write_fltarr ("xx_dd1", DataCount);
    for (y=0; y< SizeFT; y++)
   for (x=0; x< SizeFT; x++)
   for (z=0; z< SizeFT; z++)
   {
       if (DataCount(x,y,z)  != 0) Data2(x,y,z) /= (double) DataCount(x,y,z);
       else Data2(x,y,z) = complex_f(0., 0.);
       // DataCount(x,y,z) = Data2(x,y,z).real();
   }
//   DataCount.info("real rec");
//   fits_write_fltarr("xx_re.fits", DataCount);
//    for (y=0; y< SizeFT; y++)
//    for (x=0; x< SizeFT; x++)
//    for (z=0; z< SizeFT; z++)
//    {
//       DataCount(x,y,z) = Data2(x,y,z).imag();
//    }
//    DataCount.info("imag rec");
//    fits_write_fltarr("xx_im.fits", DataCount);
   // Data2(SizeFT /2, SizeFT/2, SizeFT/2) =  Data1(NcRadon3D/2,NlRadon3D/2);
 
FFTN_3D CFFT3D_local;
   CFFT3D_local.CenterZeroFreq = True;
   CFFT3D_local.fftn3d(Data2,True);
   for (x=0;x<Nx;x++)
   for (y=0;y<Nx;y++)
   for (z=0;z<Nx;z++) Cube(x,y,z) = Data2(x,y,z).real();
   // Cube.info("REC");
}
 

/***************************************/
