/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000
**
**    File:  IM_Radon.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Radon transform and reconstruction
**    -----------
**
******************************************************************************/

#include "IM_Radon.h"
#include "IM_IO.h"
#include "FFTN_1D.h"
#include "FFTN_2D.h"

#define DEBUG_RAD 0

/***************************************/

void Radon::alloc(int Nl_Imag, int Nc_Imag, type_radon RM, int RadonNl, int RadonNc)
{
   Alloc = True;
   RadonMethod=RM;
   Nl = Nl_Imag;
   Nc = Nc_Imag;
   if (RadonNl <= 0) ProjectNbr = Nl_Imag+Nc_Imag;
   else ProjectNbr = RadonNl;

   if (RadonNc <= 0)
   {
      ResolNbr = Nc_Imag;
      if ((RadonMethod == RADON_FFT_2) || (RadonMethod == RADON_FSS)) ResolNbr *= 2;
   }
   else  ResolNbr = RadonNc;
   // NormRadonFFT = (float) ResolNbr / (float) (Nl*Nc);


   if (RadonMethod == RADON_FFT)
   {
      if (ProjectNbr != Nl_Imag+Nc_Imag)
      {
         cout << "Error: with this Radon method, Projections number must be equal to Nl+Nc " << endl;
         exit(-1);
      }
      if (ResolNbr != Nc)
      {
         cout << "Error: with this Radon method, Resolution number must be equal to Nc " << endl;
         exit(-1);
      }
   }
   else if (RadonMethod == RADON_FINITE )
   {
      ProjectNbr = Nl_Imag+1;
      ResolNbr = Nc_Imag;
   }
   else if ((RadonMethod == RADON_FFT_2)|| (RadonMethod == RADON_FSS))
   {
      if (Nl_Imag != Nc_Imag)
      {
	 cout << "Error: input image must be square when using this transform ... " << endl;
         exit(-1);
      }
//       if (RadonMethod == RADON_FFT_2)
//       {
// 	 if (Nl_Imag % 2 !=  1)
// 	 {
// 	    cout << "Error: the number of lines must must be an odd number when using this transform  ... " << endl;
//             exit(-1);
//          }
//       }
      if (RadonMethod == RADON_FFT_2)
      {
         FSSR.alloc(Nl_Imag);
	 ProjectNbr = FSSR.nlr();
	 ResolNbr = FSSR.ncr();
      }
      else
      {
         SSR.alloc(Nl_Imag);
         ProjectNbr = SSR.nlr();
	 ResolNbr = SSR.ncr();
      }
   }
//    if (Verbose == True)
//    {
//       cout << "ALLOC RADON: " << StringRadon(RadonMethod) << " Image size = " << Nl_Imag << " " << Nc_Imag << endl;
//       cout << "Radon image size: Nl = " << ProjectNbr << " Nc = " << ResolNbr << endl;
//    }
   make_table();
}

/***************************************/

void Radon::reset_param()
{
   Alloc=False;
   ProjectNbr = 0;
   Verbose=False;
   ResolNbr=0;
   RadonMethod=DEF_RADON;
   FilterWidth=DEF_RADON_FILTER_WIDTH;
   FilterSigma=DEF_RADON_SIGMA;
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
   if (TabStep != NULL)
   {
      delete [] TabStep;
      TabStep = NULL;
   }
}

/***************************************/

void Radon::make_table()
{
   int t,a;

   if ((RadonMethod == RADON_PROJECT_BACK) || (RadonMethod == RADON_PROJECT_FFT))
   {
      SinTab = new double [ProjectNbr];
      CosTab = new double [ProjectNbr];
      for(a=0; a < ProjectNbr; a++)
      {
         double Val = (double) a*PI/  (double) ProjectNbr;
         SinTab[a] = (float) sin(Val);
         CosTab[a] = (float) cos(Val);
      }
   }
   else if (RadonMethod != RADON_FINITE)
   {
      TabAngle = new double[ProjectNbr];
      TabStep = new double[ProjectNbr];
      double D = (int)(Nl/2);
      for (t=0; t <= Nl/2; t++)
      {
         double A = atan( (double) (t) / D);
	 // if (t == Nl/2) cout << "Angle = " << A << " " << A/PI*180. << endl;
         //double Step = sqrt(4.*t*t+Nl*Nl) / (double) (ResolNbr);
         TabAngle[t] = A;
         TabAngle[Nl-t] = PI/2. - A;
         TabAngle[Nl+t] = PI/2. + A;
         if (t > 0) TabAngle[ProjectNbr -t] = PI - A;
      }
      for (t=0; t < ProjectNbr; t++)
      {
         int u,v;
	 double L;
         if (t <= Nl / 2)
         {
            u = Nc/2;
            v = t;
         }
         else if (t <= Nc+Nl/2)
         {
            u = (Nc+1)/2 - t + Nl/2;
            v = Nl/2;
         }
         else
         {
             u = -Nc/2;
             v = (Nl+Nc)-t;
         }
 	 L = sqrt((double)(u)*u+(double)(v)*v);
         if (Nl %2 == 0) TabStep[t] = (2.*L) / (double) (ResolNbr);
	 else TabStep[t] = (2.*L) / (double) (ResolNbr);
      }
   }
}

/***************************************/

void Radon::transform(Ifloat &Data1, Ifloat &RadonImage)
{
   // If the class is not allocated,
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) || (Nl != Data1.nl()) || (Nc != Data1.nc()))
      alloc(Data1.nl(), Data1.nc());

   if (Data1.nl() != Data1.nc())
   {
      cout << "Error: this RADON transform need square images ... " << endl;
      exit(-1);
   }
   if ((RadonImage.nl() != ProjectNbr) || ( RadonImage.nc() != ResolNbr))
                                        RadonImage.resize(ProjectNbr,ResolNbr);
   RadonImage.init();
   if (Verbose == True)
   {
       cout << "RADON transform: input image size = " << Data1.nl() << " " << Data1.nc() << endl;
       cout << "                 output image size = " << RadonImage.nl() << " " << RadonImage.nc() << endl;
   }

   switch (RadonMethod)
   {
       case RADON_PROJECT_BACK:
       case RADON_PROJECT_FFT:
             project(Data1, RadonImage);
              break;
       case RADON_FSS:
               SSR.transform(Data1, RadonImage);
	       break;
       case RADON_FFT_2:
               FSSR.transform(Data1, RadonImage);
	       break;
             break;
       case RADON_FFT:
              fft_trans(Data1, RadonImage);
              break;
       case  RADON_FINITE:
              finite_radon(Data1, RadonImage);
              break;
   }
}


/***************************************/

void Radon::project(Ifloat &Data1, Ifloat &Data2)
{
   int x, y, a;
   int  Xcenter, Ycenter, Ileft, Iright;
   float val, Position, Aleft, Aright;

   /* Project each pixel in the image */
   Xcenter = Nc/2;
   Ycenter = Nl/2;
   // cout << "project " << Xcenter << " " << endl;
   // cout << " Data1 " <<  Data1.nl() << " " << Data1.nc() << endl;
   // cout << " Data2 " <<  Data2.nl() << " " << Data2.nc() << endl;

   for(y= -Ycenter; y<Ycenter; y++)
   for(x= -Xcenter; x<Xcenter; x++)
   {
      val = Data1(y+Ycenter,x+Xcenter);
      for(a=0; a< ProjectNbr; a++)
      {
          // Fast inaccurate projection
          Position = (float)((double) x * CosTab[a] + (double) y * SinTab[a]);

          // Slower more accurate projection
          // Calculate coordinates of left and right pixels
          if (Position >= 0.0)
          {
             Ileft = (int) Position;
             Iright = Ileft + 1;
          }
          else
          {
             Iright = (int) Position;
             Ileft = Iright - 1;
          }

          /* Calculate contribution to left and right pixels */
          Aleft = Iright - Position;
          Aright = Position - Ileft;
          if (Ileft<0) Ileft += ResolNbr;
          if (Iright<0) Iright += ResolNbr;
          if ((Ileft >= 0) && (Ileft  < ResolNbr))
                                          Data2(a,Ileft) += (val*Aleft);
          if ((Iright >= 0) && (Iright < ResolNbr))
                                          Data2(a,Iright) += (val*Aright);
       }
   }
   // cout << "EXIT project " << endl;
}

/***************************************/

void Radon::recons(Ifloat &RadonImage, Ifloat &ImageRec)
{
   // If the class is not allocated,
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) ||
        (ProjectNbr  !=  RadonImage.nl()) || (ResolNbr != RadonImage.nc()))
   {
       int Nr =  RadonImage.nc();
       alloc(Nr,Nr,DEF_RADON,RadonImage.nl(),RadonImage.nc());
   }
   if (( ImageRec.nl() != Nl) || ( ImageRec.nc() != Nc)) ImageRec.resize(Nl,Nc);
   ImageRec.init();
   if (Verbose == True)
   {
       cout << "RADON reconstruct: input image size = " <<  RadonImage.nl() << " " <<  RadonImage.nc() << endl;
       cout << "                   output image size = " << ImageRec.nl() << " " <<  ImageRec.nc() << endl;
   }
   switch (RadonMethod)
   {
        case RADON_PROJECT_BACK:
              back_project(RadonImage , ImageRec);
              break;
        case RADON_PROJECT_FFT:
              invfft_project(RadonImage, ImageRec);
              break;
        case RADON_FFT_2:
	      if (NbrIter < 2) FSSR.inverse(RadonImage, ImageRec);
	      else FSSR.iter_inverse(RadonImage, ImageRec, NbrIter);
	      break;
	case RADON_FSS:
	      if (NbrIter < 2) SSR.inverse(RadonImage, ImageRec); // SSR.inverse_fft(RadonImage, ImageRec);
	      else  SSR.iter_inverse(RadonImage, ImageRec, NbrIter);
	      break;
	case RADON_FFT:
              fft_rec(RadonImage, ImageRec);
              break;
        case  RADON_FINITE:
              inv_finite_radon(RadonImage, ImageRec);
              break;
   }
}

/***************************************/

void Radon::back_project(Ifloat &Data1, Ifloat &Data2)
{
   int x,y,i,j,Xcenter,Ycenter,Xdim,Ydim;
   float val, Position, Scale;

   Xdim = Ydim = Nc;
   // Project each pixel in the image
   Scale = 1.0/(Nl*Nc);
   Xcenter = Xdim/2;
   Ycenter = Ydim/2;
   for(y = -Ycenter; y<Ycenter; y++)
   for(x = -Xcenter; x<Xcenter; x++)
   {
      val = 0.0;
      for(i=0; i< ProjectNbr; i++)
      {
         Position = x * CosTab[i] + y * SinTab[i];
         j = iround(Position);
         if (j < 0) j += ResolNbr;
         if ((j>=0) && (j < ResolNbr)) val += Data1(i,j);
      }
      Data2(y+Ycenter,x+Xcenter) = val*Scale;
   }
}

/***************************************/

void Radon::invfft_project(Ifloat &RadonImage, Ifloat &Image, Bool Permut)
{
   int i,j,t,r,x,y,x1,y1;
   int NcRadon = ResolNbr ;
   int NlRadon = ProjectNbr;
   Icomplex_f Data1(NlRadon, NcRadon,(char*)"data1");
   int SizeFT  = NcRadon;   // sampling size in Fourier space
   Icomplex_f Line(1, NcRadon,(char*)"Line");
   Icomplex_f FTLine(1, NcRadon,(char*)"Line");
   for (i=0; i < NlRadon; i++)
   {
      for (j=0; j < NcRadon; j++)
      {
         int jj = (Permut == False) ? j : (j + NcRadon/2) % NcRadon;
         Line(j) = complex_f(RadonImage(i,jj), 0.);
       }
      fft1d (Line, FTLine, 1);
      for (j=0; j < NcRadon; j++)  Data1(i,j)  =  FTLine(j) ;
   }

   Icomplex_f Data2(SizeFT,  SizeFT,(char*)"data1");
   for (y=0; y < SizeFT; y++)
   for (x=0; x < SizeFT; x++)
   {
      y1 = y - SizeFT / 2;
      x1 = x - SizeFT / 2;
      if  (x1==0 && y1==0)
      {
         t =  NlRadon/2;
         r =  NcRadon/2;
      }
      else
      {
         r = (int)(sqrt((double) x1*x1 + y1*y1)+0.5);
         t = (int)(atan2((double) y1, (double) x1)*NlRadon/PI+0.5);
         if (t < 0) {t += NlRadon; r = - r;}
         r += NcRadon /2;
      }

      if ((t>=0) && (t<NlRadon) && (r>=0) && (r<NcRadon))
            Data2(y,x) = Data1(t,r);
       else Data2(y,x) = complex_f(0.,0.);
   }
   fft2d (Data2, -1);
   for (y=0; y<Nl; y++)
   for (x=0; x<Nc; x++)
   {
      i = y - Nl/2 + SizeFT/2;
      j = x - Nc/2 + SizeFT/2;
      Image(y,x) = Data2(i,j).real();
   }
}

/***************************************/

void Radon::filter(Ifloat &RadonImage)
{
   int i,j;
   float *Filter, exponent;
   int Debug = Verbose;
   int NlRadon = RadonImage.nl();
   int NcRadon = RadonImage.nc();
   int Width = FilterWidth;
   float StdDev = FilterSigma;
   if (Verbose == True)
   {
       cout << "Filtering in Fourier space of each projection of the Radon transform " << endl;
       cout << "   Filter width = " << Width << " Filter StdDev = " <<  StdDev << endl;
   }

   // Initialize the Filter
   Filter = new float [NcRadon+1];
   if (Width > NcRadon/2) Width = NcRadon/2;
   Filter[NcRadon/2] = 1;

   for (i=1; i<Width; i++) Filter[NcRadon/2+i] = Filter[NcRadon/2-i] = i;
   for (i=Width; i<=NcRadon/2; i++)
   {
      exponent = -(Width-i)*(Width-i)/(StdDev*StdDev);
      Filter[NcRadon/2+i] = Filter[NcRadon/2-i]  = Width * (float)exp((double)exponent);
      // if (Debug == True) printf("Filter[%d] = %f\n", i, Filter[i]);
   }

   if (Debug == True)
       for (i=0; i<NcRadon; i++) printf("Filter[%d] = %f\n", i, Filter[i]);

   // Filter each projection
   Icomplex_f Line(1, NcRadon,(char*)"Line");
   Icomplex_f FTLine(1, NcRadon,(char*)"Line");
   for(i = 0; i < NlRadon; i++)
   {
      for (j = 0; j < NcRadon; j++) Line(j)  = complex_f(RadonImage(i,j), 0.);
      fft1d (Line, FTLine, 1);
      for (j = 0; j < NcRadon; j++) FTLine(j) *= Filter[j];
      fft1d (FTLine, Line, -1);
      for (j=0; j < NcRadon; j++)  RadonImage(i,j) = (Line(j)).real();
   }

   delete [] Filter;
}

/***************************************/

void Radon::finite_radon(Ifloat &Image, Ifloat & Radon)
{
   int n,i,j,k,l;
   if ((Alloc == False) || (Nl != Image.nl()) || (Nc != Image.nc()))
                                   alloc(Image.nl(), Image.nc(), RADON_FINITE);

   int p = Nc;
   // float pp = p;

   for (i=0; i < p; i++)
   {
      for (l=0; l < p; l++) Radon(i,l) = 0;
      n = i;
      for (l=0; l < p; l++)
      {
         n = n - i;
         if (n < 0) n += p;
         j = n - 1;
         for (k=0; k < p; k++)
         {
            j++;
	    if (j >= p) j -= p;
 	    Radon(i,j) += Image(k,l);
         }
      }
      // for (l=0; l < p; l++) Radon(i,l) /= pp;
   }

   for (i=0; i < p; i++) Radon(p,i) = 0;
   for (l=0; l < p; l++)
   for (k=0; k < p; k++) Radon(p,l) = Radon(p,l) + Image(k,l);
   // for (i=0; i < p; i++) Radon(p,i) /= pp;
}


/***************************************/

void Radon::inv_finite_radon(Ifloat & Radon, Ifloat &Image)
{
   int n,i,j,k,l;
   if ((Alloc == False) || (Nl != Image.nl()) || (Nc != Image.nc()))
                                   alloc(Image.nl(), Image.nc(), RADON_FINITE);

   if ((Alloc == False) ||
        (ProjectNbr  != Image.nl()) || (ResolNbr !=  Image.nc()))
   {
       int Nr = Radon.nc();
       alloc(Nr,Nr, RADON_FINITE, Radon.nl(), Radon.nc());
   }

   Image.init();

   float a=0;
   int p = Nc;
   float pp = p;

   for (j=0; j < p; j++) a += Radon(0,j);
   a /= (float) (pp*pp);

   for (l=0; l < p; l++)
   {
      for (k=0; k < p; k++) Image(k,l) = Radon(p,l) / pp;
      n = -l;
      for (i=0; i < p; i++)
      {
         n += l;
         if (n >= p) n -= p;
         k = n - 1;
         for (j=0; j < p; j++)
         {
            k ++;
            if (k >= p)  k -= p;
           Image(k,l) = Image(k,l) + Radon(i,j) / pp - a;
         }
      }
   }
}


/***************************************/

void Radon::invfft_line(Icomplex_f & TF_Rad, Ifloat &RadonImage)
{
   int i,j;
   int NlRadon = ProjectNbr;
   int NcRadon = ResolNbr ;

   if ((RadonImage.nl() != TF_Rad.nl()) || (RadonImage.nc() != TF_Rad.nc()))
     	                                    RadonImage.resize(NlRadon,NcRadon);
   // io_write_ima_complex_f  ("xx_cf", TF_Rad);
   FFTN_1D FFT1D;
   FFT1D.CenterZeroFreq = True;
   cdarray Line(NcRadon);
   // Ifloat RadonImageImag(NlRadon, NcRadon, "II");
   for (i=0; i< NlRadon; i++)
   {
      for (j=0; j < NcRadon; j++) Line(j) =  TF_Rad(i,j);
      if (RadonMethod == RADON_FFT_2)  Line(NcRadon/2) *= 2.;
      // if (Line(0).imag() != 0) Line(0) = complex_d(0, 0);
      FFT1D.fftn1d (Line.buffer(), NcRadon, True);

      // for (j=0; j < NcRadon; j++) RadonImageImag(i,j) = (float) (Line(j)).imag();
      for (j=0; j < NcRadon; j++) RadonImage(i,j) = (float) (Line(j)).real();

#if DEBUG_RAD
      for (j=0; j < NcRadon; j++) Line(j) = complex_d(RadonImage(i,j), 0);
      FFT1D.fftn1d (Line.buffer(), NcRadon, False);

      if (RadonMethod == RADON_FFT_2) Line(NcRadon/2) /= 2.;
      for (j=0; j < NcRadon; j++)
           TF_Rad(i,j) = complex_f(Line(j).real(), Line(j).imag());
#endif
   }
   // INFO_X(RadonImageImag, "Imaginary Part");
}

/***************************************/

void Radon::fft_trans(Ifloat &Image, Ifloat &RadonImage)
{
   Icomplex_f TF_Rad(ProjectNbr,ResolNbr,(char*)"data1");
   fft_trans_cf(Image, TF_Rad);
   invfft_line(TF_Rad,RadonImage);

#if DEBUG_RAD
   Image.init();
   fft_rec_cf(TF_Rad, Image);
   io_write_ima_float("xx_ima_rec.fits", Image);
   printf("Rec Min = %f, max = %f \n", min(Image), max(Image));
#endif
}

/***************************************/

void Radon::fft_line(Ifloat &RadonImage, Icomplex_f & TF_Rad, Bool Permut)
{
   int i,j;
   int NlRadon = ProjectNbr;
   int NcRadon = ResolNbr ;

   FFTN_1D FFT1D;
   FFT1D.CenterZeroFreq = True;
   cdarray Line(NcRadon);
   for (i=0; i< NlRadon; i++)
   {
      for (j=0; j < NcRadon; j++) Line(j)=complex_d((double) RadonImage(i,j), 0);
      FFT1D.fftn1d(Line.buffer(), NcRadon, False);

      if (RadonMethod == RADON_FFT_2)  Line(NcRadon/2) /= 2.;
      for (j=0; j < NcRadon; j++)
           TF_Rad(i,j) = complex_f(Line(j).real(), Line(j).imag());
   }
}
/***************************************/

void Radon::fft_trans_cf(Ifloat &Image, Icomplex_f & TF_Rad)
{
   int t,r,x,y,x1,y1;
   double u,v;
   FFTN_2D FFT;
   FFT.CenterZeroFreq = True;

   if ((Alloc == False) || (Nl != Image.nl()) || (Nc != Image.nc()))
                                        alloc(Image.nl(), Image.nc(), RADON_FFT);
   Icomplex_f TF_Ima(Nl,Nc,(char*)"data1");
   int NlRadon = ProjectNbr;
   int NcRadon = ResolNbr ;
   int Nl2 = Nl >> 1;
   int Nc2 = Nc >> 1;
   int NcRadon2 = NcRadon >> 1;

   if ((TF_Rad.nl() != ProjectNbr) || ( TF_Rad.nc() != ResolNbr))
                                       TF_Rad.resize(ProjectNbr,ResolNbr);
   if ((TF_Rad.nl() != NlRadon) || (TF_Rad.nc() != NcRadon))
   {
       cout << "Error: bad image size for the Radon transform ..." << endl;
       exit(-1);
   }

#if DEBUG_RAD
     cout << "fft_trans: Nl = " << Nl << " Nc = " << Nc << endl;
     cout << "  NlRadon = " <<   NlRadon << " NcRadon  = " <<  NcRadon << endl;
#endif

   // fft2d (Image, TF_Ima, 1);
   FFT.fftn2d(Image, TF_Ima, False);
   Ifloat DataCount1(Nl,Nc,"DD");
   // io_write_ima_complex_f  ("xx_cf", TF_Ima);

#if DEBUG_RAD
   for (y=0; y< Nl; y++)
   for (x=0; x< Nc; x++)  DataCount1(y,x) = TF_Ima(y,x).real();
   io_write_ima_float("xx_fft_ima.fits", DataCount1);
#endif

   for (t=0; t < NlRadon; t++)
   for (r=0; r < NcRadon; r++) TF_Rad(t,r) = complex_f(0.,0.);

   for (t=0; t < NlRadon; t++)
   {
      double StepProj = TabStep[t];
      double Angle = TabAngle[t];
      int rsym, DirLine = 1;
      double rr = -NcRadon/2.*StepProj;
      u = rr * cos(Angle);
      v = rr * sin(Angle);
      if ((u >= Nc2) && (v < Nc2)) DirLine = -1;
      if (Angle == 0) DirLine = 1;

      // if (ABS(Angle-PI/4.) < 0.01)  cout << Angle/PI*180. << " " << x1 << " " << y1 << " " << u << " " << v << endl;

      for (r=0; r <= NcRadon2; r++)
      {
         rr = DirLine * ((double) r-(int)(NcRadon2)) * StepProj;
         u = rr * cos(Angle);
         v = rr * sin(Angle);

         x1 = iround(u);
         y1 = iround(v);
	 x = x1 + Nc2;
	 y = y1 + Nl2;
         if ((y>=0) && (y< Nl) && (x>=0) && (x<Nc))
	 {
	     // if (DataCount1(y,x) == 0)
	     {
	         TF_Rad(t,r) = TF_Ima(y,x);
	         DataCount1(y,x)+=1;
             }
         }
       // if (ABS(Angle) < 0.0001)  cout <<   " " << x1 << " " << y1 << " " << u << " " << v << endl;
	 int xs = -x1;
	 int ys = -y1;
	 x = xs + Nc2;
	 y = ys + Nl2;
         rsym = 2*((int)(NcRadon2)) - r;

 	 if ((y>=0) && (y< Nl) && (x>=0) && (x<Nc) && (rsym >= 0) &&
	      (rsym < NcRadon))
	      {
	         // if (DataCount1(y,x) == 0)
	         {
 	     	      TF_Rad(t,rsym) = TF_Ima(y, x);
	              DataCount1(y,x)+=1;
		}
              }
	 else if ((rsym >= 0) && (rsym < NcRadon))
         {
	      cout << "Error occured during the transformation ..." << endl;
	      exit(-1);
	       TF_Rad(t,rsym) = complex_f(TF_Rad(t,r).real(), -TF_Rad(t,r).imag());
	  }
// 	  if (ABS(Angle) < 0.0001)
// 	  {
// 	    cout << " r = " << r << "  " << TF_Rad(t,r) << " " << rsym << " " <<  TF_Rad(t,rsym);
// 	    cout << " x,y = (" << x1+ Nc/2 << "," << y1+ Nl/2 << ")";
// 	    cout << " xs,ys = (" << x << "," << y << ")" << endl;
// 	  }
      }
   }
   // io_write_ima_float("xx_fft_rad.fits", DataCount1);

#if DEBUG_RAD
   Ifloat DataCount(NlRadon,NcRadon,(char*)"DD");
   for (y=0; y< NlRadon; y++)
   for (x=0; x< NcRadon; x++)  DataCount(y,x) = TF_Rad(y,x).real();
   io_write_ima_float("xx_fft_rad.fits", DataCount);

   cout << "RE  TF_Rad = " << TF_Rad(NlRadon/2, NcRadon/2).real() <<
           " IM TF_Rad = " << TF_Rad(NlRadon/2, NcRadon/2).imag() << endl;
#endif
}

/***************************************/

void Radon::fft_rec_cf(Icomplex_f & TF_Rad, Ifloat &Image)
{
   int t,r,x,y,x1,y1;
   double u,v;
   FFTN_2D FFT;
   FFT.CenterZeroFreq = True;

   // If the class is not allocated,
   // or if the class was allocated for another image size
   // then we do it with default parameter
   if ((Alloc == False) ||
        (ProjectNbr  !=  TF_Rad.nl()) || (ResolNbr != TF_Rad.nc()))
   {
       int Nr =  TF_Rad.nc();
       // cout << " ALLOC ... " << endl;
       alloc(Nr,Nr,RADON_FFT,TF_Rad.nl(),TF_Rad.nc());
   }
   int NlRadon = ProjectNbr;
   int NcRadon = ResolNbr ;

   if (( Image.nl() != Nl) || ( Image.nc() != Nc)) Image.resize(Nl,Nc);
   Image.init();
   if (Verbose == True)
   {
       cout << "RADON reconstruct: input image size = " <<  TF_Rad.nl() << " " <<  TF_Rad.nc() << endl;
       cout << "RADON reconstruct: output image size = " << Image.nl() << " " <<  Image.nc() << endl;
   }

   Icomplex_f TF_Ima(Nc,Nc,"data1");
   // TF_Ima.init();
   Ifloat DataCount(Nc,Nc,"Count");

#if DEBUG_RAD
   Ifloat TT(ProjectNbr,   NcRadon,(char*)"data1");
   for (y=0; y< ProjectNbr  ; y++)
   for (x=0; x< NcRadon  ; x++)  TT(y,x) = TF_Rad(y,x).real();
   INFO_X(TT,"tt");
   io_write_ima_float("xx_fft_rad_rec.fits",  TT);
#endif

   for (t=0; t < NlRadon; t++)
   {
      double StepProj = TabStep[t];
      double Angle = TabAngle[t];
      int rsym,DirLine = 1;
      r = 0;
      double rr = -NcRadon/2.*StepProj;
      u = rr * cos(Angle);
      v = rr * sin(Angle);
      if ((u >= Nc/2) && (v < Nc/2)) DirLine = -1;
      if (Angle == 0) DirLine = 1;

      for (r=0; r <= NcRadon/2; r++)
      {
         rr = DirLine * ((double) r-(int)(NcRadon/2)) * StepProj;
         u = rr * cos(Angle);
         v = rr * sin(Angle);

         x1 = iround(u);
         y1 = iround(v);
	 x = x1 + Nc/2;
	 y =  y1 + Nc/2;
         if ((y>=0) && (y < Nl) && (x>=0) && (x<Nc))
	 {
	   // if (DataCount(y,x) == 0)
	   {
	      TF_Ima(y,x) +=  TF_Rad(t,r);
              DataCount(y,x) += 1.;
	   }
	   // TF_Ima(y,x) =  TF_Rad(t,r);
 	 }
         x1 = -x1;
	 y1 = -y1;
	 x = x1 + Nc/2;
	 y = y1 + Nc/2;
	 // rsym = NcRadon - r;
         rsym = 2*((int)(NcRadon/2)) - r;
	 // if (t == 0)
	 //  cout << "r = " << r << " X1 = " << -x1+Nc/2 << " rsym = " << rsym << " X1s = " << x1+Nc/2 << " " << y << endl;

	 if ((y>=0) && (y< Nl) && (x>=0) && (x<Nc) && (rsym >= 0)
	     && (rsym < NcRadon))
	 {
	    // if (DataCount(y,x) == 0)
	    {
	       TF_Ima(y,x) += TF_Rad(t,rsym);
               DataCount(y,x) += 1.;
	    }
	    // TF_Ima(y,x) =  TF_Rad(t,rsym);
	 }
      }
   }

   // TF_Ima(Nc /2, Nc/2) =  TF_Rad(NlRadon/2, NcRadon/2);
   // cout << "RE = " << TF_Rad(NlRadon/2, NcRadon/2).real() <<
   //         " IM = " <<TF_Rad(NlRadon/2, NcRadon/2).imag() << endl;
   // io_write_ima_float("xx_dd.fits", DataCount);

   for (y=0; y< Nc; y++)
   for (x=0; x< Nc; x++)
       if (DataCount(y,x)  > 1) TF_Ima(y,x) /= (float) DataCount(y,x);

   // cout << Nl << " " << Nc << " RE = " << TF_Ima(Nc/2, Nc/2).real() <<
   //        " IM = " << TF_Ima(Nc/2, Nc/2).imag() << endl;

#if DEBUG_RAD
  cout << "DDmin = " << min(DataCount) << " DDmax = " <<  max(DataCount) << endl;
  for (y=0; y< Nc; y++)
  for (x=0; x< Nc; x++)  DataCount(y,x) = TF_Ima(y,x).real();
  io_write_ima_float("xx_fft_ima_rec.fits", DataCount);
  cout << "FFT size = " << TF_Ima.nl() << " " <<  TF_Ima.nc() << endl;
#endif

   // fft2d (TF_Ima, -1);
   FFT.fftn2d(TF_Ima, True);

   for (y=0; y<Nl; y++)
   for (x=0; x<Nc; x++) Image(y,x) = TF_Ima(y,x).real();

#if DEBUG_RAD
   printf("Rec Min = %f, max = %f \n", min(Image), max(Image));
#endif
}

/***************************************/

void Radon::fft_rec(Ifloat &RadonImage, Ifloat &Image, Bool Permut)
{
   Icomplex_f TF_Rad(ProjectNbr, ResolNbr,(char*)"data1");
   fft_line(RadonImage, TF_Rad,Permut);
   fft_rec_cf(TF_Rad, Image);
}

/***************************************/
/*************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 4) || (File_Name_In[L-1] != 'd')
                || (File_Name_In[L-2] != 'a')
		|| (File_Name_In[L-3] != 'r')
                || (File_Name_In[L-4] != '.'))
    {
        strcat (File_Name_Out, ".rad");
    }
}

/*************************************************************************/

/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/

void Radon::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nl", (long) Nl,(char*)"Nl",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Nc",(long) Nc,(char*)"Nc",&status))
     PrintError( status );
  if ( ffpkyj(fptr, (char*)"Nlr", (long) ProjectNbr,(char*)"Number of project",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Ncr",(long) ResolNbr,(char*)"resolution",&status))
     PrintError( status );
  if ( ffpkyj(fptr, (char*)"Type_Tra", (long) RadonMethod, (char*)StringRadon(RadonMethod), &status))
	 PrintError( status );
  if ( ffpkyj(fptr, (char*)"FormatIn",(long)FormatInputImag,(char*)"Format", &status))
     PrintError( status );
}


/****************************************************************************/

void Radon::write (char *Name, Ifloat &Ima)
/* new version with fits */
{
 char filename[256];
 fitsfile *fptr;
 int status;
 //int i,j,s;
 //float *Ptr;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[3];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 //long fpixels[3];
 //long lpixels[3];

/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */

/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 naxis = 2;
 naxes[0] = radon_nc();
 naxes[1] = radon_nl();

 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // write the header of the multiresolution file
  mr_io_fill_header(fptr);

  nelements = naxes[0] * naxes[1];
  if ( ffppre(fptr, group, firstpixel, nelements, Ima.buffer(), &status) )
              PrintError( status );

 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );
// cout << " end of write fits " << endl;
}


/****************************************************************************/

void Radon::read (char *Name, Ifloat & Ima)
{
    // for fits
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT];
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];

     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;

#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif

    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) )
         PrintError( status );

    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

     nelements = naxes[0] * naxes[1];
     // cout << " begin to read " << endl;
     Ima.alloc(naxes[1], naxes[0], "read io");

    if (ffgkyj(fptr,(char*)"Nl", &mon_long, comment, &status)) PrintError( status );
    int Nli = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Nc", &mon_long, comment, &status)) PrintError( status );
    int Nci = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Nlr", &mon_long, comment, &status)) PrintError( status );
    int ProjectNbri = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Ncr", &mon_long, comment, &status)) PrintError( status );
    int ResolNbri = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
    type_radon RadTransi = (type_radon) mon_long;
    if (ffgkyj(fptr,(char*)"FormatIn", &mon_long, comment, &status)) PrintError( status );
    FormatInputImag = (type_format)mon_long;

    Ptr = Ima.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

    if ( ffclos(fptr, &status) ) PrintError( status );
    alloc(Nli, Nci, RadTransi, ProjectNbri, ResolNbri );
// cout << " end of read fits file " << endl;
}

/****************************************************************************/
