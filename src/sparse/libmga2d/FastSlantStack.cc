/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/08/02
**
**    File:  FastSlantStack.cc
**
*******************************************************************************
**
**    DESCRIPTION  Fast Slant Stack radon transform of an image
**    -----------
**
**
*******************************************************************************/


#include "FastSlantStack.h"
#include "IM_IO.h"
#include "FSS.h"

/**********************************************************************/

inline float REM(float Arg, float Denominator)
{
   int Val = (int)(Arg / Denominator);
   float Vali = (float) Val * Denominator;
   return (Arg - Vali);
}

/**********************************************************************/

void slant_stack_radon::inverse_fft(Ifloat & Radon, Ifloat &Data)
{
   Ifloat PSR(Nr, Nr, "buff");
   Icomplex_f TF_Buff(Nr,Nr, "buff");
   Icomplex_f TF_Buff1(Nr,N, "buff");
   Icomplex_f TF_Buff2(N,Nr, "buff");
   Icomplex_f TF_Rec(Nr,Nr);
   intarray TX0,TY0;
   intarray TX1,TY1;
   intarray TX2,TY2;
   Ifloat Buff1(Nr,N, "buff");
   Ifloat Buff2(N,Nr, "buff");
   TX0.alloc(2*Nr);
   TY0.alloc(2*Nr);
   TX1.alloc(2*Nr);
   TY1.alloc(2*Nr);
   TX2.alloc(2*Nr);
   TY2.alloc(2*Nr);
   TF_Buff.init();
   TF_Buff1.init();
   TF_Buff2.init();
   Buff1.init();
   Buff2.init();
   int i,j,l,N2 = (N+1) / 2;
   int Nr2 = Nr / 2;

   FFT2D.fftn2d(Data, TF_Buff);

   for (int l=0; l < N; l++)
   {
       cfarray FTLine(Nr);
       int Np,x2,y2;
       x2 = l - N2;
       y2 = Nr2;
       symgetline (x2, y2, Nr2, N2,  TX1,  TY1, Np);

       // for (i=0; i < Nr; i++) FTLine(i) = complex_f(Radon(N-1-l,i), 0.);
       for (i=0; i < Nr; i++) FTLine(i) = complex_f(Radon(N+l,Nr-1-i), 0.);
       FFT1D.fftn1d (FTLine);
       for (i=0; i < Nr; i++)
       {
          if ( (TY1(i) >= 0) && (TX1(i) > 0) && (TY1(i) < TF_Buff1.nl()) && (TX1(i) < TF_Buff1.nc()))
	  {
	     TF_Buff1(TY1(i),TX1(i)) += FTLine(i);
             Buff1(TY1(i),TX1(i)) += 1;
	  }
       }
    }
    // cout << "1" << endl;
    // io_write_ima_float("xx_buff1.fits", Buff1);

    for (l=0; l < N; l++)
    {
       cfarray FTLine(Nr);
       int Np,x2,y2;
       x2 = Nr2;
       y2 = l - N2;
       symgetline (x2, y2, N2, Nr2,  TX2,  TY2, Np);
//        if (Np != Nr)
//        {
//           cout << "PB: Np(Nr) = " << Np << endl;
// 	  exit(0);
//        }
       // for (i=0; i < Nr; i++) FTLine(i) += complex_f(Radon(N-1-l+N,i),0.);
       for (i=0; i < Nr; i++) FTLine(i) += complex_f(Radon(N-1-l,i),0.);
       FFT1D.fftn1d (FTLine);
       for (i=0; i < Nr; i++)
       {
          if ( (TY2(i) >= 0) && (TX2(i) >= 0) && (TY2(i) < TF_Buff2.nl()) && (TX2(i) < TF_Buff2.nc()))
	  {
	     TF_Buff2(TY2(i),TX2(i)) += FTLine(i);
             Buff2(TY2(i),TX2(i)) += 1;
	  }
       }
    }
    //cout << "2" << endl;
    // io_write_ima_float("xx_buff2.fits", Buff2);

    for (i = 0; i < Nr; i++)
    for (j = 0; j < N; j++)
    {
       if (Buff1(i,j) != 0) TF_Buff1(i,j) /= Buff1(i,j);
    }
    //cout << "Buff 3 " << endl;

    for (i = 0; i < N; i++)
    for (j = 0; j < Nr; j++)
    {
       if (Buff2(i,j) != 0) TF_Buff2(i,j) /= Buff2(i,j);
    }
    //cout << "Buff 4 " << endl;
    for (i = 0; i < Nr; i++)
    for (j = 0; j < Nr; j++)
    {
       int ii =  i - Nr2;
       int jj = j - Nr2;
       if (ABS(ii) >= ABS(jj)) { if (TF_Buff1(i,j/2).real() != 0) TF_Buff(i,j) = TF_Buff1(i,j/2);}
       else if (TF_Buff2(i/2,j).real() != 0) TF_Buff(i,j) = TF_Buff2(i/2,j);
    }
    //cout << "Buff 5 " << endl;
    FFT2D.ifftn2d(TF_Buff);
    // cout << "Buff 6 " << endl;

    for (i=0; i < N; i++)
    for (j=0; j < N; j++)
       Data(i,j) = TF_Buff(i+N2,j+N2).real();
    //cout << "Buff 7 " << endl;
}



/**********************************************************************/

void slant_stack_radon::iter_inverse(Ifloat & Radon, Ifloat &Data, int Niter)
{
   if (UseLeviCode == True)
   {
   cout << "iter_inverse" << endl;
       extern int NUM_OF_ITER;
       NUM_OF_ITER = Niter;
       Idouble DataI(N,N);
       Idouble ResI(Nr,Nr);
       Idouble ResR(Nr,Nr);
       Idouble DataR(N,N);
       // ResR = Radon;
       for (int i=0; i < Nr; i++)
       for (int j=0; j < Nr; j++) ResR(i,j) = Radon(i,j);

       int Verb = (Verbose == True) ?  1 : 0;
       Verb = 1;
       Inv_SS(DataR.buffer(), DataI.buffer(), ResR.buffer(), ResI.buffer(), N, Verb);
       for (int i=0; i < N; i++)
       for (int j=0; j < N; j++) Data(i,j) = DataR(i,j);
    }
   else
   {
   int i,j,Nl = Data.nl();
   int Nc = Data.nc();
   Ifloat RadonRec(Nr,Nr,"rec");
   Ifloat DataRec(Nl,Nc,"rec data");
   Ifloat  RadonAux;


   if (OptIter == True) RadonAux.alloc(Nr,Nr,"rec");
    inverse(Radon, Data);

    if (Positive == True) threshold(Data);
    for (int Iter = 0; Iter < Niter; Iter++)
    {
       transform(Data, RadonRec);
       for (i = 0; i < Nr; i++)
       for (j = 0; j < Nr; j++) RadonRec(i,j) = Radon(i,j) - RadonRec(i,j);
       // INFO_X(RadonRec, "RESI radon:");

       inverse(RadonRec, DataRec);
       if (OptIter == True)
       {
          float Alpha=0;
	  float Num=0, Den=0;
          transform(DataRec, RadonAux);
	  for (i = 0; i < Nr; i++)
          for (j = 0; j < Nr; j++)
	  {
 	     Num += RadonRec(i,j)*RadonAux(i,j);
             Den += RadonAux(i,j)*RadonAux(i,j);
  	  }
	  if (Den >  FLOAT_EPSILON) Alpha  = Num / Den;
	  else Alpha =1;
 	  if (Alpha > 10) Alpha = 10;
	  else if (Alpha < 1) Alpha = 1;
 	  if (Verbose == True)
          cout << "Iter = " << Iter + 1<< " ErrRadon = " <<  RadonRec.sigma() << " Alpha = " <<  Alpha << endl;
	  for (i = 0; i < N; i++)
          for (j = 0; j < N; j++) Data(i,j) += DataRec(i,j) * Alpha;
       }
       else
       {
           if (Verbose == True)
          cout << "Iter = " << Iter + 1<< " Err = " <<  RadonRec.sigma() << endl;
          Data += DataRec;
       }
       if (Positive == True) threshold(Data);
   }
   }
}

/**********************************************************************/

void slant_stack_radon::inverse(Ifloat & Radon, Ifloat &Data)
{
   if (UseLeviCode == True)
   {
       // cout << "RUN ASS " << endl;
       int k,l;
       float Norm = N*N*2.;
       Idouble DataI(N,N);
       Idouble ResI(Nr,Nr);
       Idouble ResR(Nr,Nr);
       Idouble DataR(N,N);
       // ResR = Radon;
       for (int i=0; i < Nr; i++)
       for (int j=0; j < Nr; j++) ResR(i,j) = Radon(i,j);
       Adj_SS(DataR.buffer(), DataI.buffer(), ResR.buffer(), ResI.buffer(), N);
       for (k=0; k < N; k++)
       for (l=0; l < N; l++) Data(k,l) = DataR(k,l) / Norm;
   }
   else
   {
   int i,j;
   Ifloat BuffY(Nr, N, "buff");
   Ifloat BuffX(N, Nr, "buff");
   TF1DDataX.init();
   TF1DDataY.init();
   Data.init();
    if ((Radon.nl() != Nr) || (Radon.nc() != Nr))
    {
       cout << "Error: the class fft_slant_stack_radon is " << endl;
       cout << "       not initialized for this radon image size.  " << endl;
       cout << "       Radon Data size = " << Radon.nl() << " " << Radon.nc() << endl;
       cout << "       Class initialization = " << Nr << " " << Nr << endl;
       exit(0);
    }

   if ((Data.nl() != N) || (Data.nc() != N)) Data.resize(N,N);
   complex_f *CFBuffY = TF1DDataY.buffer();
   complex_f *CFBuffX = TF1DDataX.buffer();
   int N2 = (N+1)/2;
   double C1 = 2.*PI / (double) N;
   for (int Dep = 0; Dep < N; Dep ++)
   {
      float Dx = - Dep +  N2;
      for (i=0; i < Nr; i++)
      {
         float Dist = N;
 	 int dd = (int) ABS(Dx);
	 int P2 = N+N2-dd;
	 int P1 = N2+dd;
         if ((i < N2-dd) || (i > N+N2+dd))
  	     for (j=0; j < N; j++)  BuffY(i,j) = BuffX(j,i) = 0;
 	 else if ((dd == 0) || ((i < P2) && (i > P1)))
	 {
	    Dist = N;
	    for (j=0; j < N; j++)
	    {
 	         BuffY(i,j) = Radon(2*N-1-Dep,Nr-1-i) / Dist;
		 BuffX(j,i) = Radon(Dep,i) / Dist;

	       // BuffY(i,j) = Radon(Dep,i) / Dist;
	       // BuffX(j,i) = Radon(Dep+N,i) / Dist;
		   // BuffY(i,j) = BuffX(j,i) = 2;
	    }
	 }
	 else
 	 {
            int ii = (i > N) ? i -P2: P1 - i;
 	    Dist = N - iround((float) ii *  (float) N / dd / 2.);
            for (j=0; j < N; j++)
            {
 	       if (((Dep < N2) && (j < Dist)  && (i < N)) ||
	           ((Dep < N2) && (j > N-Dist)  && (i >= N)) ||
	           ((Dep > N2) && (j < Dist)  && (i >= N)) ||
	           ((Dep > N2) && (j > N-Dist)  && (i < N)))
               {
	           if (Dist > 0)
		   {
		      //BuffY(i,j) = Radon(Dep,i) / Dist;
	              //BuffX(j,i) = Radon(Dep+N,i) / Dist;
		      BuffY(i,j) = Radon(2*N-1-Dep,Nr-1-i) / Dist;
		      BuffX(j,i) = Radon(Dep,i) / Dist;
		   }
		   else BuffY(i,j) = BuffX(j,i) = 0;
		   // BuffY(i,j) = BuffX(j,i) = 1;
 	       }
	       else BuffY(i,j) = BuffX(j,i) = 0;
    	    }
	 }
      }
//       if (Dep == 0) io_write_ima_float("xx_rbuffy1.fits", BuffY);
//       if (Dep == N/4) io_write_ima_float("xx_rbuffy2.fits", BuffY);
//       if (Dep == N-1) io_write_ima_float("xx_rbuffy3.fits", BuffY);
//       if (Dep == 3*N/4) io_write_ima_float("xx_rbuffy4.fits", BuffY);
//       if (Dep == 0) io_write_ima_float("xx_rbuffx1.fits", BuffX);
//       if (Dep == N/4) io_write_ima_float("xx_rbuffx2.fits", BuffX);
//       if (Dep == N-1) io_write_ima_float("xx_rbuffx3.fits", BuffX);
//       if (Dep == 3*N/4) io_write_ima_float("xx_rbuffx4.fits", BuffX);
      for (j=0; j < N; j++)
      {
         float y = Dx* ((float) j /  (float) (N-1) - 0.5);
	 for (i=0; i < Nr; i++)
	 {
	    FTLineY(i) = complex_f(BuffY(i,j),0.);
	    FTLineX(i) = complex_f(BuffX(j,i),0.);
	 }
	 FFT1D.fftn1d (FTLineX);
	 FFT1D.fftn1d (FTLineY);
         for (i=0; i< Nr; i++)
         {
            double Pos,Sre,Sim,Rre,Rim,Norm;
	    float t = REM(y, Nr);
            Pos = (double) i - (double) N + 0.5;
            Sre =  cos(C1*Pos*t);
            Sim = -sin(C1*Pos*t);
            Rre = FTLineY(i).real() * Sre - FTLineY(i).imag() * Sim;
            Rim = FTLineY(i).real() * Sim + FTLineY(i).imag() * Sre;
	    Norm = 1. / (float) N;
// 	    k = i - Nr / 2;
// 	    Norm = (k == 0) ? sqrt(1./8.) / sqrt((float) Nr):
// 	                       sqrt((float) ABS(k)/2) / (float) Nr;
	    Rre *= Norm;
	    Rim *= Norm;
            TF1DDataY(i,j) += complex_f((float) Rre, (float) Rim);
	    Rre = FTLineX(i).real() * Sre - FTLineX(i).imag() * Sim;
            Rim = FTLineX(i).real() * Sim + FTLineX(i).imag() * Sre;
	    Rre *= Norm;
	    Rim *= Norm;
	    TF1DDataX(i,j) += complex_f((float) Rre, (float) Rim);
	 }
      }
  }

  CFBuffY = TF1DDataY.buffer();
  CFBuffX = TF1DDataX.buffer();
  int N4 = Nr/4;
  for (j=0; j < N; j++)
  {
      FFT1D.fftn1d (CFBuffY, Nr, True);
      FFT1D.fftn1d (CFBuffX, Nr, True);
      for (i=0; i < N; i++)
      {
          Data(i,j) += 0.5*CFBuffY[i+N4].real();
	  Data(j,i) += 0.5*CFBuffX[i+N4].real();
      }
      CFBuffY += Nr;
      CFBuffX += Nr;
  }
  }
}

/**********************************************************************/


void slant_stack_radon::alloc(int Dim)
{
    Positive = False;
    N =  Dim;
    Nr = 2*N;
    FFT1D.CenterZeroFreq = True;
    TF1DDataY.alloc(Nr,N);
    TF1DDataX.alloc(Nr,N);
    FTLineX.alloc(Nr);
    FTLineY.alloc(Nr);
}


/**********************************************************************/

void slant_stack_radon::transform(Ifloat &Data, Ifloat & Radon)
{
   if (UseLeviCode == True)
   {
      Radon.init();
      if ((Radon.nl() != Nr) || (Radon.nc() != Nr)) Radon.resize(Nr,Nr);
      Idouble  DataI(N,N);
      Idouble  DataR(N,N);
      Idouble  ResI(Nr,Nr);
      Idouble RadonR(Nr,Nr);
      for (int j=0; j < N; j++)
      for (int i=0; i < N; i++) DataR(i,j) = Data(i,j);

      // DataR = Data;
      // cout << "RUN SS " << endl;
      SS(RadonR.buffer(), ResI.buffer(), DataR.buffer(), DataI.buffer(), N);
      for (int j=0; j < Nr; j++)
      for (int i=0; i < Nr; i++) Radon(i,j) = RadonR(i,j);
   }
   else
   {
   int i,j;

   Ifloat BuffY(Nr, N, "buff");
   Ifloat BuffX(N, Nr, "buff");
   if ((Radon.nl() != Nr) || (Radon.nc() != Nr))
                                           Radon.resize(Nr,Nr);
   Radon.init();
   // 1D FFT for each column using zero padding on the y-axis
   complex_f *CFBuffY = TF1DDataY.buffer();
   complex_f *CFBuffX = TF1DDataX.buffer();
   int N4 = Nr/4;
   for (j=0; j < N; j++)
   {
      for (i=0; i < N; i++)
      {
          // Extract one column from the data with zero padding
          CFBuffY[i+N4] = complex_f(Data(i,j), 0.);
	  // Extract one  line from the data with zero padding
	  CFBuffX[i+N4] = complex_f(Data(j,i), 0.);
      }
      FFT1D.fftn1d (CFBuffY, Nr);
      FFT1D.fftn1d (CFBuffX, Nr);
      CFBuffY += Nr;
      CFBuffX += Nr;
   }

   int N2 = (N+1)/2;
   double C1 = 2.*PI / (double) N;
   for (int Dep = 0; Dep < N; Dep ++)
   {
      float Dx = Dep -  N2;
      for (j=0; j < N; j++)
      {
         float y = Dx * ((float) j /  (float) (N-1) - 0.5);
         for (i=0; i< Nr; i++)
         {
            double Pos,Sre,Sim,Rre,Rim;
	    float t = REM(y, Nr);
            Pos = (double) i - (double) N + 0.5;
            Sre =  cos(C1*Pos*t);
            Sim = -sin(C1*Pos*t);
            Rre = TF1DDataY(i,j).real() * Sre - TF1DDataY(i,j).imag() * Sim;
            Rim = TF1DDataY(i,j).real() * Sim + TF1DDataY(i,j).imag() * Sre;
            FTLineY(i) = complex_f((float) Rre, (float) Rim);

	    Rre = TF1DDataX(i,j).real() * Sre - TF1DDataX(i,j).imag() * Sim;
            Rim = TF1DDataX(i,j).real() * Sim + TF1DDataX(i,j).imag() * Sre;
            FTLineX(i) = complex_f((float) Rre, (float) Rim);
	 }
         FFT1D.fftn1d (FTLineY, True);
	 FFT1D.fftn1d (FTLineX, True);
         for (i=0; i < Nr; i++)
	 {
	    BuffY(i,j) = (float) (FTLineY(i)).real();
	    BuffX(j,i) = (float) (FTLineX(i)).real();
	 }
      }
//       if (Dep == 0) io_write_ima_float("xx_buffy1.fits", BuffY);
//       if (Dep == N/4) io_write_ima_float("xx_buffy2.fits", BuffY);
//       if (Dep == N-1) io_write_ima_float("xx_buffy3.fits", BuffY);
//       if (Dep == 0) io_write_ima_float("xx_buffx1.fits", BuffX);
//       if (Dep == N/4) io_write_ima_float("xx_buffx2.fits", BuffX);
//       if (Dep == N-1) io_write_ima_float("xx_buffx3.fits", BuffX);
      for (i=0; i < Nr; i++)
      for (j=0; j < N; j++)
      {
        // Radon(Dep,i) += BuffY(i,j);
	// Radon(Dep+N,i) += BuffX(j,i);
	 Radon(2*N-1-Dep,Nr-1-i) += BuffY(i,j);
	 Radon(Dep,i)   += BuffX(j,i);
      }
  }
  }
}
/**********************************************************************/
/**********************************************************************/
//
// Fast Slant Stack radon transform of an image in the Fourier domain
//
/**********************************************************************/
/**********************************************************************/

void fft_slant_stack_radon::alloc(int N)
{
    Positive = False;
    Nl =  N;
    Nc =  N;
    Nl1 = Nl*2+1;
    Nc1 = (Nc % 2 == 0) ? Nc + 1: Nc;
    Nl2 = (Nl % 2 == 0) ? Nl + 1: Nl;
    Nc2 = Nc*2+1;
    Nlr = Nc1 + Nl2;
    Ncr = Nl1;
    Nl22 = Nl2 / 2;
    Nc22 = Nc2 / 2;
    Nl12 = Nl1 / 2;
    Nc12 = Nc1 / 2;
    FFT1D.CenterZeroFreq = True;
    TF_Buff1.alloc(Nl1,Nc1,"buff1");
    TF_Buff2.alloc(Nl2,Nc2,"buff2");
    TX0.alloc(2*Nl1);
    TY0.alloc(2*Nc2);
    TX1.alloc(2*Nl1);
    TY1.alloc(2*Nl1);
    TX2.alloc(2*Nc2);
    TY2.alloc(2*Nc2);
    Buff1.alloc(Nl1, Nc1, "buff");
    Buff2.alloc(Nl2, Nc2, "buff");
    TF_Rec.alloc(Nl1,Nc2, "rec buff");
}


 /******************************************************************************/

void fft_slant_stack_radon::transform(Ifloat &Data, Ifloat & Radon)
{
   int i,j;

    if ((Data.nl() != Nl) || (Data.nc() != Nc))
    {
          cout << "Error: the class fft_slant_stack_radon is " << endl;
	  cout << "       not initialized for this image size.  " << endl;
	  cout << "       Data size = " << Data.nl() << " " << Data.nc() << endl;
	  cout << "       Class initialization = " << Nl << " " << Nc << endl;
	  exit(0);
    }
    if ((Radon.nl() != Nlr) || (Radon.nc() != Ncr))
                                           Radon.resize(Nlr,Ncr);
    if ((Radon.nl() != Nlr) || (Radon.nc() != Ncr))
                                           Radon.resize(Nlr,Ncr);
    // cout << Data.nl() << " " <<  Data.nc() << endl;
    // cout << TFBuff1.nl() << " " <<  TFBuff1.nc() << endl;
    TF_Buff1.init();
    TF_Buff2.init();
    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++)
                TF_Buff1(i+Nl1/4,j) = complex_f(Data(i,j), 0.);

    for (i=0; i < Nl; i++)
    for (j=0; j < Nc; j++)
                TF_Buff2(i,j+Nc2/4) = complex_f(Data(i,j), 0.);

    fft_radon(TF_Buff1, TF_Buff2, Radon);
}

/******************************************************************************/

void fft_slant_stack_radon::fft_radon(Icomplex_f & TFBuff1,
                                      Icomplex_f & TFBuff2, Ifloat & Radon)
{
   int i;
   int Np,x2,y2;
  // Ifloat PSR(Nlr, Ncr, "buff");

   // Buff1.init();
   // Buff2.init();
    FFT2D.fftn2d(TFBuff1);
    // complex_f F1 = TFBuff1(Nl1/2,Nc1/2);
    // cout << " TF_BUFF1(0,0) = " << sqrt(F1.real()*F1.real()+ F1.imag()*F1.imag()) << endl;

    // Extract line in buffer 1
    for (int l1=0; l1 < Nc1; l1++)
    {
       // cout << l1 << endl;
       cfarray FTLine(Nl1);
       x2 = l1 - Nc12;
       y2 = Nl12;
       symgetline (x2, y2, Nl12, Nc12,  TX1,  TY1, Np);

       if (Np != Nl1)
       {
          cout << "PB: Np(Nl1) = " << Np << endl;
	  exit(0);
       }
       if ((TX1(Np/2) != Nc12) || (TY1(Np/2) != Nl12))
       {
          cout << "PB: TX1(Np/2) = " << TX1(Np/2) << endl;
	  cout << "    TY1(Np/2) = " << TY1(Np/2) << endl;
	  exit(0);
       }
       for (i=0; i < Nl1; i++) FTLine(i) = TFBuff1(TY1(i),TX1(i));
//        for (i=0; i < Nl1; i++)
//        {
//           // cout << TY1(i) << " " << TX1(i) << endl;
//           Buff1(TY1(i),TX1(i)) += 1;
// 	  // PSR(l1,i) = sqrt(FTLine(i).real()*FTLine(i).real() +
// 	  //              FTLine(i).imag()*FTLine(i).imag());
//        }
       FFT1D.fftn1d (FTLine, True);
//        if (l1 == 128)
//        {
//           cout <<  " TX1 = " <<  TX1(127) <<  " " << TX1(128) <<  " " << TX1(129)  << endl;
//           cout <<  " TY1 = " <<  TY1(127) <<  " " << TY1(128) <<  " " << TY1(129)  << endl;
//           cout <<  " 127 = " << PSR(l1,127) <<" 128 = " << PSR(l1,128) << " 129 = " << PSR(l1,129) << endl;
//        }
       for (i=0; i < Nl1; i++) Radon(l1,i) = FTLine(i).real();
    }
    // io_write_ima_float("xx_buff1.fits", Buff1);

    // Extract line in buffer 2
    FFT2D.fftn2d(TFBuff2);
    // F1 = TFBuff2(Nl2/2,Nc2/2);
    // cout << " TF_BUFF2(0,0) = " << sqrt(F1.real()*F1.real()+ F1.imag()*F1.imag()) << endl;

    for (int l2=0; l2 < Nl2; l2++)
    {
       Np = 0;
       cfarray FTLine(Nc2);
       x2 = Nc22;
       y2 = l2 - Nl22;
       symgetline (x2, y2, Nl22, Nc22,  TX2,  TY2, Np);

       if (Np != Nc2)
       {
          cout << "PB: Np(Nc2) = " << Np << endl;
	  exit(0);
       }
       if ((TX2(Np/2) != Nc22) || (TY2(Np/2) != Nl22))
       {
          cout << "PB: TX2(Np/2) = " << TX2(Np/2) << endl;
	  cout << "    TY2(Np/2) = " << TY2(Np/2) << endl;
	  cout << "   Nc22 = " <<Nc22 << "   Nl22 = " << Nl22 <<endl;
	  cout << "   x2 = " << x2 << "   y2 = " << y2 <<endl;
	  cout << "   Np = " << Np << "   Nc2 = " << Nc2 <<endl;
	  for (i=0; i < Nc2; i++) cout << " " << TX2(i);
	  cout << "Y " << endl;
	  for (i=0; i < Nc2; i++) cout << " " << TY2(i);
	  cout << "  " << endl;
	  exit(0);
       }
       for (i=0; i < Nc2; i++) FTLine(i) = TFBuff2(TY2(i),TX2(i));
//        for (i=0; i < Nc2; i++)
//        {
//           // cout << TY1(i) << " " << TX1(i) << endl;
//           Buff2(TY2(i),TX2(i)) += 1;
// 	  // PSR(l2+Nc1,i) = sqrt(FTLine(i).real()*FTLine(i).real() +
// 	  //             FTLine(i).imag()*FTLine(i).imag());
//        }
       FFT1D.fftn1d (FTLine, True);
       for (i=0; i < Nc2; i++) Radon(l2+Nc1,i) = FTLine(i).real();
    }
    // io_write_ima_float("xx_buff2.fits", Buff2);
    // io_write_ima_float("xx_psr.fits", PSR);
}

/*********************************************************************/

void fft_slant_stack_radon::inverse(Ifloat & Radon, Ifloat &Data)
{
   int i,j;
   if ((Radon.nl() != Nlr) || (Radon.nc() != Ncr))
   {
       cout << "Error: the class fft_slant_stack_radon is " << endl;
       cout << "       not initialized for this radon image size.  " << endl;
       cout << "       Radon Data size = " << Radon.nl() << " " << Radon.nc() << endl;
       cout << "       Class initialization = " << Nlr << " " << Ncr << endl;
       exit(0);
   }

   inv_fft_radon(Radon,TF_Rec);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
       Data(i,j) = TF_Rec(i+Nl1/4,j+Nc2/4).real();
}

/*********************************************************************/

void fft_slant_stack_radon::inv_fft_radon(Ifloat & Radon, Icomplex_f & TFBuff)
{
   int i,j;
   Ifloat PSR(Nlr, Ncr, "buff");
   Bool TDEBUG = False;

   TFBuff.init();
   TF_Buff1.init();
   TF_Buff2.init();
   Buff1.init();
   Buff2.init();
   for (int l1=0; l1 < Nc1; l1++)
   {
       cfarray FTLine(Nl1);
       int Np,x2,y2;
       x2 = l1 - Nc12;
       y2 = Nl12;
       symgetline (x2, y2, Nl12, Nc12,  TX1,  TY1, Np);

       if (Np != Nl1)
       {
          cout << "PB: Np(Nc2) = " << Np << endl;
	  exit(0);
       }
       for (i=0; i < Nl1; i++) FTLine(i) = complex_f(Radon(l1,i), 0.);
       FFT1D.fftn1d (FTLine);
       for (i=0; i < Nl1; i++)
       {
          TF_Buff1(TY1(i),TX1(i)) += FTLine(i);
          Buff1(TY1(i),TX1(i)) += 1;
       }
    }

    for (int l2=0; l2 < Nl2; l2++)
    {
       cfarray FTLine(Nc2);
       int Np,x2,y2;
       x2 = Nc22;
       y2 = l2 - Nl22;
       symgetline (x2, y2, Nl22, Nc22,  TX2,  TY2, Np);
       if (Np != Nc2)
       {
          cout << "PB: Np(Nc2) = " << Np << endl;
	  exit(0);
       }
       for (i=0; i < Nc2; i++) FTLine(i) = complex_f(Radon(l2+Nc1,i),0.);
       FFT1D.fftn1d (FTLine);
       for (i=0; i < Nc2; i++)
       {
          TF_Buff2(TY2(i),TX2(i)) += FTLine(i);
          Buff2(TY2(i),TX2(i)) += 1;
        }
    }

    for (i = 0; i < Nl1; i++)
    for (j = 0; j < Nc1; j++)
    {
       if (Buff1(i,j) == 0)
       {
          int v;
          int u =  2*(j - Nc12) + Nc22;
	  int ii = i - Nl12;
	  if (ii > 0) v = (ii+1) / 2 + Nl22;
	  else if (ii < 0) v = (ii-1) / 2 + Nl22;
	  else v = Nl22;
	  if (u < 0) u = 0;
	  if (v < 0) v = 0;
 	  if (u >= Nc2) u = Nc2-1;
	  if (v >= Nl2) v = Nl2-1;
	  if ((TDEBUG == True) && (Buff2(v,u) == 0))
	  {
	     cout << "PB: (u,v) = " <<u << " " << v << endl;
	     cout << "    Buff2(v,u) = " << Buff2(v,u)<< endl;
	     // exit(0);
 	  }
	  if (Buff2(v,u) != 0)  TF_Buff1(i,j) = TF_Buff2(v,u) / Buff2(v,u);
       }
       else TF_Buff1(i,j) /= Buff1(i,j);
    }
    // cout << "Buff 2 " << endl;

    for (i = 0; i < Nl2; i++)
    for (j = 0; j < Nc2; j++)
    {
       if (Buff2(i,j) == 0)
       {
          int u,v;
          int jj =  j - Nc22;
 	  if (jj > 0) u = (jj+1) / 2 + Nc12;
	  else if (jj < 0) u = (jj-1) / 2 + Nc12;
	  else u = Nc12;
	  v  = 2*(i - Nl22) + Nl12;
	  if (u < 0) u = 0;
	  if (v < 0) v = 0;
	  if (u >= Nc1) u = Nc1-1;
	  if (v >= Nl1) v = Nl1-1;
 	  if ((TDEBUG == True) && (Buff1(v,u) == 0))
	  {
	     cout << "PB: (u,v) = " <<u << " " << v << endl;
	     cout << "    Buff1(v,u) = " << Buff1(v,u)<< endl;
	     // exit(0);
 	  }
	  if (Buff1(v,u) != 0)  TF_Buff2(i,j) = TF_Buff1(v,u);
       }
       else TF_Buff2(i,j) /= Buff2(i,j);
    }

    for (i = 0; i < Nl1; i++)
    for (j = 0; j < Nc2; j++)
    {
       int ii =  i - Nl12;
       int jj = j - Nc22;
       if (ABS(ii) >= ABS(jj)) TFBuff(i,j) = TF_Buff1(i,j/2);
       else TFBuff(i,j) = TF_Buff2(i/2,j);
    }
    FFT2D.ifftn2d(TFBuff);

    // Ifloat X(Nl1,Nc2,"x");
    // for (i=0; i < Nl1; i++)
    // for (j=0; j < Nc2; j++) X(i,j) = TFBuff(i,j).real();
    // io_write_ima_float("xx_rec.fits", X);
}

/*********************************************************************/

void fft_slant_stack_radon::iter_inverse(Ifloat & Radon, Ifloat &Data, int Niter)
{
   int i,j,Iter;
   // float RegulParam = 0.0;

   if ((Radon.nl() != Nlr) || (Radon.nc() != Ncr))
   {
       cout << "Error: the class fft_slant_stack_radon is " << endl;
       cout << "       not initialized for this radon image size.  " << endl;
       cout << "       Radon Data size = " << Radon.nl() << " " << Radon.nc() << endl;
       cout << "       Class initialization = " << Nlr << " " << Ncr << endl;
       exit(0);
   }
   Ifloat RadonRec(Nlr,Ncr,"rec");
   Ifloat DataRec(Nl,Nc,"rec data");
   Ifloat Grad(Nl,Nc,"grad data");
   // RegulIma RI;

   inverse(Radon,Data);

   if (Positive == True) threshold(Data);
   for (Iter = 0; Iter < Niter; Iter++)
   {
       cout << "Iter = " << Iter + 1<< endl;
       transform(Data, RadonRec);
       if (Iter == 0) io_write_ima_float("xx_rec.fits",  RadonRec);
       if (Iter == Niter-1) io_write_ima_float("xx_rec1.fits",  RadonRec);
       // fft_radon(Data, RadonRec);
       // INFO_X(Radon, "RADON:");
       // INFO_X(RadonRec, "RADONREC:");
       for (i = 0; i < Nlr; i++)
       for (j = 0; j < Ncr; j++) RadonRec(i,j) = Radon(i,j) - RadonRec(i,j);
       INFO_X(RadonRec, "RESI radon:");
       inverse(RadonRec,  DataRec);
       Data += DataRec;
       INFO_X(Data, "Data:");
       if (Positive == True) threshold(Data);
   }
}

/**********************************************************************/
