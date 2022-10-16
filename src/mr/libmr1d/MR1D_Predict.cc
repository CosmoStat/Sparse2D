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
**    Date:  20/07/00
**    
**    File:  MR1D_Predict.cc
**    
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------    
**
**    this module contains the routines for 1D prediction
**
**************************************************************************/ 

#include "GlobalInc.h"
#include "MR1D_Obj.h"
#include "MatrixOper.h"
#include "NeuNet.h"
#include "MR1D_Predict.h"
#include "IM1D_IO.h"


/***************************************************************************/

char *StringExMethod (extrapol_type   type)
{
    switch (type)
    {
        case  EX_POLYNOM:
              return ("Polynomial prediction."); 
        case  EX_NN:
              return ("Neural Network Prediction.");
         case EX_MR_NN:
              return ("Multiresolution Neural network prediction.");
        case  EX_AR:
              return ("Autoregressive Model.");              
        case  EX_MR_AR_PRED1:
              return ("Autoregressive model per scale."); 
         case  EX_MR_AR_PRED2:
              return ("Multiresolution Autoregressive Model.");             
        default: 
              return ("Undefined method");
     }
}

/***************************************************************************/

void extrap_usage(extrapol_type Method)
{
    fprintf(OUTMAN, "        [-P predict_method]\n");
    for (int i = 0; i < NBR_EX_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                          StringExMethod((extrapol_type)i));
    fprintf(OUTMAN, "             default is %s\n", StringExMethod((extrapol_type) Method));
}


/***************************************************************************/

double old_pred_pol_deg2(fltarray &Data, int Np, int NLast, int DistPred)
// Calcule the polynome of degree 2 fitting the Nlast points
// of the data
// Data(i) = a x^2 + b x + 1    
//           and i = [Np-NLast,Np-1]
// Np: IN = Number of pixels in Data
// Data: IN = input data
// Nlast: IN = number of points in Data to be used
// DistPred: IN = Calculate the prediction at a distance DistPred+1 from
//                the last pixel (DistPred=0) for the next pixel prediction)
//                and return this value      
{
   int k;
   dblarray MAT_X(1,3);  
   dblarray W;
   int Pos = Np-1;
   int Pmin = MAX(0,Pos-NLast/2);
   int Pmax = Pos;
   int Nx = Pmax-Pmin+1;
   W.alloc(Nx);
   // cout << "Pmax = " << Pmax << " Pmin = " << Pmin << endl;
   for (k=0; k<Nx; k++) 
   {
      // double D = Nx-1;
      double x = (double) k / (double) (Nx-1);
      W(k) = x;
      // if (k == 1) cout << "k = " << k << " w = " << W(k) << " D = " << D << endl;
   }
   // cout << "pred_pol_deg2" << endl;
   // cout << "Data.nx() " << Data.nx() << " Np = " << Np << " Nx = " << Nx << endl;
   fit1d_pol_deg2(Data, MAT_X, W, Np,Nx);
   double P, x = NLast + DistPred;
   P = MAT_X(0)*x*x+ MAT_X(1)*x +MAT_X(2);
   // cout << "end pred_pol_deg2" << endl;
   return P;
}

/***************************************************************************/
 
double pred_pol_deg2(fltarray &Data, int Np, int NLast, int DistPred)
// Calcule the polynome of degree 2 fitting the Nlast points
// of the data
// Data(i) = a x^2 + b x + 1    
//           and i = [Np-NLast,Np-1]
// Np: IN = Number of pixels in Data
// Data: IN = input data
// Nlast: IN = number of points in Data to be used
// DistPred: IN = Calculate the prediction at a distance DistPred+1 from
//                the last pixel (DistPred=0) for the next pixel prediction)
//                and return this value      
{
   int k;
   dblarray MAT_X(1,3);  
   dblarray W;
   int Nx = NLast;
   int Pmin = MAX(0,Np-NLast);
   int Pmax = Np-1;
   Nx = Pmax-Pmin+1;
   W.alloc(Nx);
   for (k=0; k<Nx; k++) 
   {
      double x = (double) k / (double)(Nx-1);
      W(k) = sqrt(x*x);
   }
   fit1d_pol_deg2(Data, MAT_X, W, Np, Nx);
   double P, x = Nx+DistPred;
   P = MAT_X(0)*x*x+ MAT_X(1)*x +MAT_X(2);  
   
//    cout << "Data: " <<  NLast << " " << Data( Np-4) 
//          << " "  << Data( Np-3) << " " << " "  << Data( Np-2) 
// 	 << " "  << Data( Np-1) << endl;
//    cout << "MAT_X = " << MAT_X(0) << " " << MAT_X(1) << " " << MAT_X(2) << endl;
//    cout << "P = " << P << endl;
   return P;
}

/***************************************************************************/
 
void MR_PRED::mk_tendancy(fltarray & Data, int Np, fltarray & TendTab)
{
    tendancy_est(Data, TendTab, NbrPixTendancy, Np);
} 

/*********************************************************************/

void MR_PRED::wave_1d_trans(fltarray &Signal, MR_1D & MR_Data, int N)
/* linear wavelet transform */
{
    int i,k,indi,b;
    int Step;
    // int SizeFilter = 4;
    // float HFilter[4] = {0.125,0.375,0.375,0.125};   
    int FirstPix = -2;
    //int SizeFilter = 4; 
    //float HFilter[4] = {0.125,0.375,0.375,0.125};   
    int SizeFilter = 5;
    float HFilter[5] = {1./16.,1./4.,3./8.,1./4.,1./16.}; 
    // int SizeFilter = 2;  
    // float HFilter[2] = {0.5,0.5};
    
    // type_border Border = MR_Data.Border;
    int Nb = MR_Data.nbr_scale();
    // fltarray Scale(N);
    int NLast = NbrPixTendancy;
    dblarray MAT_Pol(1,3);
    dblarray MAT_Pol1(1,2);
    
    fltarray Data(N);
    for (i = 0; i < N; i++) Data(i) = Signal(i);
    for (b = 0; b < Nb - 1; b++)
    {
       float Coef=0.;
          // MR_Data.scale(Scale, b);
	  // cout << Scale.min() << " " << Scale.max() << endl;
	  if (TPredBord == P_POL2) fit1d_pol_deg2(Data, MAT_Pol, N, NLast);
          else if (TPredBord == P_POL1)
	               fit1d_pol_deg1(Data, MAT_Pol1, N, NLast);
		       // calc_pol_deg1(Data, N, NLast, MAT_Pol1);

          // copy the data in the cube  
          for (i = 0; i < N; i++) MR_Data(b,i) = Data(i);

          Step = iround (pow((double)2., (double) b));
	  // cout << "Scale " << b << " Step = "  << Step << endl;
	  // if (TPredBord == P_POL2) 
	  //  cout << " POL = " << MAT_Pol(0) << " " << MAT_Pol(1) << " " << MAT_Pol(2) << endl;
//           for (k=0;k<2;k++)
// 	  {
// 	     indi = k*Step;
//              double x = NLast + indi;
//              Coef = (float) (MAT_Pol(0)*x*x+ MAT_Pol(1)*x +MAT_Pol(2));		
// 	     cout << indi << ": Coef Extrapol = " << Coef << endl; 
//           }

          int BorderPos = N;
          for (i = 0; i < N; i ++)
          {
	      Data(i) = 0.;
	      // BorderPos = MAX(i+1,10);
	      for (k=FirstPix;k<SizeFilter+FirstPix;k++)
	      {
	         indi = i - k*Step;
  		 
		 if (indi < 0) Coef = MR_Data(b,-indi);
		 else if (indi < BorderPos) Coef = MR_Data(b,indi);
		 else
		 {
		    double x;
		    switch (TPredBord)
		    {
		       case P_CONT:
		         Coef = MR_Data(b,BorderPos - 1);
			 break;
		       case P_MIRROR:
 		         indi  = 2 * (BorderPos - 1) - indi;
 		         Coef = MR_Data(b,indi);
			 break;
		       case P_MIRROR_2:
 		         indi  = 2 * (BorderPos - 1) - indi;
 		         Coef = 2.*MR_Data(b,BorderPos-1) - MR_Data(b,indi);
			 break;
                       case P_POL1:
		         indi -= BorderPos;
			 x = NLast + indi;
                         Coef = (float) (MAT_Pol1(0)*x+MAT_Pol1(1));
			 break;
                       case P_POL2:
		         indi -= BorderPos;
			 x = NLast + indi;
                         Coef = (float) (MAT_Pol(0)*x*x+ MAT_Pol(1)*x +MAT_Pol(2));		
		         break;
		     }
                  }
                  Data(i) += HFilter[k-FirstPix] * Coef; 
	      }
	  }

         //  calcul the wavelet coefficients  
         for (i = 0; i < N; i++) MR_Data(b,i) -= Data (i);
    }
    
    /* copy the low resolution signal in the cube */
    for (i = 0; i < N; i++) MR_Data(Nb-1,i) = Data(i);
}


/***************************************************************************/

void MR_PRED::get_mr_ar1_model(MR_1D & MR_Data, int Np, int User_AR_Order)
   // set TabAR_OrderPerScale and TabAR_ModelPerScale
   // considering the lags.
   // if AR_AutoDetect == True calculate the AR order per scale
   //                       AR order are always calculated using lags
   //                       in the scale
   // if Stationary == True, the last smooth scale is also analyzed
   // MR_Data: IN = wavelet data
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // User_AR_Order: IN = AR order to use (only if AR_AutoDetect == False) 
{
    int i,b;
    int Nx = MR_Data.size_ima_np();
    int Step= (NoLag == False) ? 2: 1;
    int NbrScale = MR_Data.nbr_scale();
    int Nbr_AR=User_AR_Order;
    int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
    fltarray Scale(Nx);
    fltarray ARModel;

    if (AR_AutoDetect == False)
    {
        ARModel.alloc(User_AR_Order);
        for (b=0; b < NusedScale; b++) TabAR_OrderPerScale(b) = User_AR_Order;
    }
    TabAR_ModelPerScale.alloc(MaxNbrAR, NbrScale);
 
    for (b=0; b < NusedScale; b++)
    {
        int BestNbrAR=User_AR_Order; 
        MR_Data.scale(Scale, b);
	// Find the best order model
        if (AR_AutoDetect == True)
        {
           get_best_ar_model(Scale, Np, BestNbrAR, ARModel,Step,b);
           TabAR_OrderPerScale(b) = MAX(MinAR_OrderPerScale,BestNbrAR);
        }
        else get_ar_model(Scale, Np, Nbr_AR, NbrTraining, ARModel, Step);
	
// cout << " BestNbrAR = " << BestNbrAR << " " <<  TabAR_ModelPerScale.nx() << endl;
// cout << "  TabAR_ModelPerScale " << TabAR_ModelPerScale.nx()  << " " << TabAR_ModelPerScale.ny()  << endl;
        for (i=0; i < BestNbrAR; i++) TabAR_ModelPerScale(i,b) = ARModel(i);

//         if (Verbose == True)
//         {
//            cout << "SCALE " << b+1<< " Nx " << Np  << " AR order   = " <<   TabAR_OrderPerScale(b)  << endl;
//            cout << "  AR = ";
//            for (i=0; i< TabAR_OrderPerScale(b); i++) 
//                    cout <<  TabAR_ModelPerScale(i,b) << " ";
//            cout << endl;
//         }
        if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
        // Step *= 2;
    }
}

/***************************************************************************/

double MR_PRED::ar1_predict(MR_1D & MR_Data, fltarray &TabSignalPred, 
                          int  User_AR_Order)
{
  int i,b;
  int Step= (NoLag == False) ? 2: 1;
  int NbrScale = MR_Data.nbr_scale();
  int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
  int Nx= MR_Data.size_ima_np();
  double Pred=0.;
  if (LastTraining  > 0) NbrTraining  =  LastTraining-FirstTraining+1;
  else  NbrTraining = 0;
  if (EvaluationMode == True)
  {
    if (LastEvaluate > 0) NbrEvaluate = LastEvaluate-FirstEvaluate+1;
    else NbrEvaluate = 0;
  }
  else NbrEvaluate = 0;
  fltarray Scale(Nx);
  fltarray TabScalePredY(Nx);

  if (Verbose == True) 
    cout << " ar1_predict " << NbrPixTendancy << " " << NbrPixTendancy << endl;

  // Find the order of the AR to be used at each scale
  // and the AR coefficients relative the AR model at each scale 
  get_mr_ar1_model(MR_Data, NbrTraining, User_AR_Order);
  int NbrAR = (int) TabAR_OrderPerScale.max();
  Step= (NoLag == False) ? 2: 1;
  if (NoLag == False)  for (b=0; b < NbrScale-1; b++) Step *= 2;
  int FirstPos = NbrAR*Step + PredDistance;
  if (FirstTraining < FirstPos) 
  {
       FirstTraining = FirstPos;
       NbrTraining = LastTraining - FirstTraining + 1;
  }
  Step= (NoLag == False) ? 2: 1;

 
  for (b=0; b < NusedScale; b++)
  {
     // Get the AR model found by set_mr_AR_model routine
     int BestNbrAR = TabAR_OrderPerScale(b); 
     fltarray ARModel(BestNbrAR);

     for (i=0; i < BestNbrAR; i++) ARModel(i) = TabAR_ModelPerScale(i,b);
     MR_Data.scale(Scale, b);

     if (Verbose == True)
     {
        cout << "SCALE " << b+1 << " AR order   = " <<  BestNbrAR  << endl;
        cout << "  AR = ";
        for (i=0; i<  BestNbrAR; i++)  cout <<  ARModel(i) << " ";
        cout << endl;
     }
     
     // Training part prediction
     for (i=FirstTraining; i<= LastEvaluate; i++)
     {
         int Pos = i;
         TabSignalPred(Pos) += ar_prediction(Scale, BestNbrAR, Pos, ARModel, Step);
     }
     
     // Predict the wavelet coefficients
//      for (i=0; i< NbrEvaluate; i++)
//      {
//          int Pos = Nx-NbrEvaluate+i;
//          MR_Data(b,Pos) = ar_prediction(Scale, BestNbrAR, Pos, ARModel, Step);
//          TabSignalPred(Pos) += MR_Data(b,Pos);
//      }
     Pred += ar_prediction(Scale, BestNbrAR, Nx, ARModel, Step);
     if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
  }
  
  if (Stationary == False)
  {
      float PredLast=0.;  

      // Haar interpolation for the last scale
      b = NbrScale-1;
      MR_Data.scale(Scale, b);
      
      // Training part prediction
      for (i=FirstTraining; i<= LastTraining; i++)
      {
         int Pos = i;
	 int DistPred=0;
	 PredLast=0.;  // = Scale(Pos-1) for Haar prediction 
	 PredLast = extr_last_scale(Scale, Pos, DistPred);	
         TabSignalPred(Pos) += PredLast;	 
      }
           
      for (i=0; i< NbrEvaluate; i++)
      {
         int Pos  = Nx - NbrEvaluate + i;
	 int DistPred=0;
	 PredLast=0.;  // = Scale(Pos-1) for Haar prediction 
	 PredLast = extr_last_scale(Scale, Pos, DistPred);	
         TabSignalPred(Pos) += PredLast;
         MR_Data(NbrScale-1,Pos) = PredLast;
      }
      PredLast = extr_last_scale(Scale,Nx,PredDistance);
      {
          cout << "Data: T = " <<  NbrPixLastPred << ": " << Scale( Nx-4) 
            << " "  << Scale( Nx-3) << " " << " "  << Scale( Nx-2) 
	    << " "  << Scale( Nx-1) << endl;
          cout << "P = " << PredLast << endl; 
      }
      Pred += PredLast;
  }
  // cout << "End AR 1 predict " << endl;
  return Pred;
}

/***************************************************************************/

void MR_PRED::get_mr_ar2_model(fltarray & Signal, MR_1D & MR_Data, int Np,   
                               int NbrTotAR, fltarray & ARModel)
// Get the AR model to predict x(t) from the wavelet coefficients.
// MR_Data: IN = wavelet data
// Signal: IN = input signal
// Np: IN = consider only the Np first in the data set (t0 = Np)
// NbrTotAR: IN = total number of coefficients used for the prediction
//                NbrTotAR = sum_i TabAR_OrderPerScale(i)
// ARModel: OUT = AR model. ARModel(0..NbrTotAR-1)	       
{
   int NbrScale = MR_Data.nbr_scale();
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
   int b,i,j,IndData;
   int t0 = Np-1-PredDistance;
   int NbrAR = (int) TabAR_OrderPerScale.max();
   int Step= (NoLag == False) ? 2: 1;
   if (NoLag == False)  for (b=0; b < NbrScale-1; b++) Step *= 2;
   // FirstPos = first position where the AR cannot be evaluate
   //            i.e. the last pixels of the data
   int FirstPos = t0 - NbrAR*Step - PredDistance;

 // cout << " XX " << TabAR_OrderPerScale.max() << " NArParam = " << NbrTotAR  << endl;
 // cout << " FirstPos  " << FirstPos  << " t0  = " <<  t0 << endl;
 // cout << "  NbrAR  " <<   NbrAR << " Step    = " << Step   << endl;
    if (Debug == True) cout << "get_mr_ar2_model" << endl;

   // FirstPos is depending of the number of training pixels
   if ((NbrTraining > 0) && (FirstPos > NbrTraining*Step))  
                                             FirstPos = NbrTraining*Step;

   // If the number of samples is not high enough, the AR order
   // must be reduced.
   if (FirstPos < 0)
   {
       NbrAR = t0 / Step;
       FirstPos = t0 - NbrAR*Step;
       if (FirstPos < 0)
       {
          cout << "Error: Number of scales is to high ... " << endl;
          exit(-1);
       }
   }

   if (Verbose == True)
        cout << "FirstPos = " << FirstPos << " Max AR = " << NbrAR << " PredDistance = "  << PredDistance << endl;

   // cout << " Step = " << Step  << "t0 = " << t0 << " NbrAR = " << NbrAR  << " FirstPos = " << FirstPos << endl;
   // Learning values: A X = B
   NbrAR = NbrTotAR;

   dblarray MAT_X(1,NbrAR);
   dblarray MAT_A(NbrAR, FirstPos);
   dblarray InvMAT_A;
   dblarray MAT_B(1,FirstPos);
   fltarray DataTend(Signal.nx());
   fltarray TabTend(Signal.nx());
   ARModel.alloc(NbrAR);
 
   MAT_X.init();
   int IndCoef;    
    
   for (i=0; i < FirstPos; i++)
   {
      if (Debug == True) cout << "Pos " << i+1 << endl;
      tend_wt_transform(Signal, DataTend, TabTend, MR_Data, Signal.nx());
      Step= (NoLag == False) ? 2: 1;
      IndCoef = 0;
      for (b=0; b < NusedScale; b++)
      {
         for (j=0; j < TabAR_OrderPerScale(b); j++)
         {
            // IndData = t0 - i - (j+1) *Step;
            IndData = t0 - i - 1 - j*Step;
            if (IndData < 0)
            {
               cout << endl;
               cout << "ERROR: IndData = " << IndData << endl;
               exit(-1);
            }
            MAT_A (IndCoef++,i) = MR_Data(b,IndData);
         }
         if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
      }
      IndData = t0 - i + PredDistance;
      if ((IndData < 0) || (IndData >= Signal.nx()))
            {
               cout << endl;
               cout << "ERROR: IndData = " << IndData << endl;
               exit(-1);
            }
      if (Stationary == False) 
           MAT_B(i) = DataTend(IndData) - MR_Data(NbrScale-1,IndData);
      else MAT_B(i) = DataTend(IndData);
   }
   MatOper MO;
   MO.lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
   for (j=0; j < NbrAR; j++) ARModel(j) = (float) MAT_X(j);
   if (Debug == True) cout << "END get_mr_ar2_model" << endl;
}

/***************************************************************************/

void MR_PRED::tend_wt_transform(fltarray & Data, 
                              fltarray & DataTend, fltarray & TabTend,
			      MR_1D & MR_Data, int Np, Bool ForceWT)
{
   int k;
   DataTend = Data;
   if (UseTend == True)
   {
       mk_tendancy(Data,Np, TabTend);
       for (k=0; k< Np; k++) DataTend(k) -= TabTend(k);
       for (k=Np; k< Data.nx(); k++) DataTend(k) = 0;
   }
   if (UseB3SplineWT == True) 
   {
      (MR_Data.image()).init();
      wave_1d_trans(DataTend, MR_Data, Np);
   }
   else if (ForceWT == True) MR_Data.transform(DataTend);
}


/***************************************************************************/

double MR_PRED::extr_last_scale(fltarray & Scale, int Np, int DistPred)
{ 
   double PredTendVal=0.;
   int indi = Np+DistPred;
   
   switch (TPredBord)
   {
      case P_POL1:
        {
	   int x = NbrPixLastPred + DistPred;
           dblarray MAT_Pol1(1,2);	               
	   fit1d_pol_deg1(Scale, MAT_Pol1, Np, NbrPixLastPred);
           PredTendVal = (float) (MAT_Pol1(0)*x+MAT_Pol1(1));
	 }        
	  break;
      case P_POL2:
        if (NbrPixLastPred > 3)
 	    PredTendVal = pred_pol_deg2(Scale,Np, NbrPixLastPred, DistPred);
        else if (Np > 0) PredTendVal = Scale(Np-1);
	break;
      case P_CONT:
        PredTendVal = Scale(Np-1);
	break;
      case P_MIRROR:
         indi  = 2 * (Np - 1) - indi;
	 if ((indi >= 0) && (indi < Np)) PredTendVal = Scale(indi);
	 break;
      case P_MIRROR_2:
         indi  = 2 * (Np - 1) - indi;
	 if ((indi >= 0) && (indi < Np)) 
	                PredTendVal = 2.*Scale(Np-1) - Scale(indi);
         break;
   }
   return PredTendVal;
}

 
   
/***************************************************************************/

double MR_PRED::ar2_predict(fltarray & Signal, MR_1D & MR_Data, 
                          fltarray &TabSignalPred, int User_AR_Order)
{
   int i,b,k;
   int Step= (NoLag == False) ? 2: 1;
   int NbrScale = MR_Data.nbr_scale();
   int Nx= MR_Data.size_ima_np();
   // int NPredict = TabSignalPred.nx() - Nx;
   float Pred;
   int Np;
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
    
   if (Debug == True) cout << " AR2 PRED " << NbrEvaluate << endl;
   fltarray Scale(Nx);
   fltarray LastScale(Nx);
   fltarray TabTend(Nx);
   fltarray ARModel;
   fltarray DataTend(Nx);
   			         
   // Find the best AR order per scale
   if (AR_AutoDetect == True)
   {
      Step = (NoLag == False) ? 2: 1;
      for (b=0; b < NusedScale; b++)
      {
         int BestNbrAR;   
         MR_Data.scale(Scale, b);
         get_best_ar_model(Scale, NbrTraining, BestNbrAR, ARModel, Step,b);
         TabAR_OrderPerScale(b) = MAX(MinAR_OrderPerScale,BestNbrAR);
         if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
      }
   }
   else 
   {
      for (b=0; b < NusedScale; b++) TabAR_OrderPerScale(b) = User_AR_Order;
   }

    
  int NArParam = 0;
  for (b=0; b < NusedScale; b++) NArParam += TabAR_OrderPerScale(b);
   get_mr_ar2_model(Signal, MR_Data, NbrTraining, NArParam, ARModel);
   int NbrAR = (int) TabAR_OrderPerScale.max();
   Step= (NoLag == False) ? 2: 1;
   if (NoLag == False)  for (b=0; b < NbrScale-1; b++) Step *= 2;
   int FirstPos = NbrAR*Step + PredDistance;
   if (FirstTraining < FirstPos) 
   {
       FirstTraining = FirstPos;
       NbrTraining = LastTraining - FirstTraining + 1;
   }
     
   if (Verbose == True)
   {
       cout << "Total number of AR coefficients = " << NArParam << endl;
       if (Verbose == True) cout << "  Model AR = " << endl;
       int IndCoef = 0;
       for (b=0; b < NusedScale; b++)
       {
          cout << "    Scale " << b+1 << " NCoef = " <<  TabAR_OrderPerScale(b) << " Coef = ";
          for (k=0; k<  TabAR_OrderPerScale(b); k++)  cout << " " << ARModel(IndCoef++);   
          cout << endl;
       }
   }
  int IndCoef;
  // b = NbrScale-1;
  // MR_Data.scale(Scale,b);
  // if (Stationary == False) io_1d_write_data("xxl.fits", Scale);

  // if (Stationary == True) cout << "STATIONARY SIGNAL " << endl;

  // Do the evaluation 
  //for (i=FirstTraining; i<= LastTraining; i++)
  //{
  //    if ((Verbose == True) && (UseTend == True))
  //                               cout << "Training Evaluation " << i+1 << endl;
 
	 
	 
  // for (i=0; i< NbrEvaluate; i++)
  for (i=FirstTraining; i<=LastEvaluate; i++)
  {
      if ((Verbose == True) && (UseTend == True))
                                 cout << "Evaluation " << i+1 << endl;
      IndCoef = 0;
      Pred = 0;
      // Np = Nx - NbrEvaluate + i;
      Np = i;
      tend_wt_transform(Signal, DataTend, TabTend, MR_Data, Np);

      Step = (NoLag == False) ? 2: 1;
      for (b=0; b < NusedScale; b++)
      {
         for (k=1; k<= TabAR_OrderPerScale(b); k++) 
	 {
	    int IndMR =  Np-1-(Step*(k-1));
	    if (IndMR < 0)
	    {
	       cout << "Error: bad IndMR = " << IndMR << " Np = " << Np << " k = " << k << " Step = " << Step << endl;
               exit(-1);
	    }
	    Pred += ARModel(IndCoef++) * MR_Data(b, IndMR);
         }
	 if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
     }
     TabSignalPred(Np) = Pred;
     int DistPred=0;
     float PredTendVal=0.;
     if (Stationary == False) 
     {
        static int P=0;
	if (P==0) 
	{
	   if (Debug == True) cout << " No stationary ... " << endl;
	   // cout << " NbrPixTendancy = " << NbrPixTendancy << endl;
	   P=1;
	}
        b = NbrScale-1;
        MR_Data.scale(Scale,b);
        PredTendVal = extr_last_scale(Scale, Np, DistPred);
        TabSignalPred(Np) += PredTendVal;
 	LastScale(Np) = PredTendVal;
     }
     if (UseTend == True)
     {
         // PredTendVal=2.*TabTend(Np-1) - TabTend(Np-2);
	 PredTendVal=extr_last_scale(TabTend,Np,DistPred);
	 TabSignalPred(Np) += PredTendVal;
     }    
  }
  // if (Stationary == False) io_1d_write_data("xx.fits", LastScale);
  
  // Do the prediction
  Pred = 0;
  Np = Nx;
  IndCoef = 0;
  if (NbrEvaluate > 0) tend_wt_transform(Signal, DataTend, Tend, MR_Data,Nx);
      
  for (b=0; b < NusedScale; b++)
  {
     for (k=1; k<= TabAR_OrderPerScale(b); k++) 
              Pred += ARModel(IndCoef++) * MR_Data(b, Np-1-(Step*(k-1)));
     if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
  }
  b = NbrScale-1;
  MR_Data.scale(Scale,b);
  if (Stationary == False)  
  {
     float PredTendVal=0.;  // = Scale(Np-1) for Haar prediction
     PredTendVal = extr_last_scale(Scale, Np, PredDistance);	
     Pred += PredTendVal;
  }  
  if (UseTend == True)
  {
     float PredTendVal=extr_last_scale(Tend,Np,PredDistance);
     Pred += PredTendVal;
  }    
  if (Debug == True) cout << "End AR 2 predict " << Pred << endl;
  return Pred;
}

/***************************************************************************/

void MR_PRED::nn1_predict(MR_1D & MR_Data,  fltarray &TabSignalPred, int User_AR_Order)
{
  int i,b,BestNbrAR;
  int Step= (NoLag == False) ? 2: 1;
  int NbrScale = MR_Data.nbr_scale();
  int Nx= MR_Data.size_ima_np();
  double NNErr=0;

  fltarray Scale(Nx);
  fltarray TabScalePredY(Nx);

  // if (Verbose == True) cout << " nn1_predict  " << NbrTraining << " " << User_AR_Order << endl;
  TabSignalPred.init();
  for (b=0; b < NbrScale; b++)
  {
     fltarray ARModel;
     MR_Data.scale(Scale, b);

     if (AR_AutoDetect == True) 
            get_best_ar_model(Scale,  NbrTraining, BestNbrAR, ARModel, Step, b);
     NN_PREDICTION NNP;
     int N = (AR_AutoDetect == True) ? BestNbrAR: User_AR_Order;
     if (Verbose == True)
        cout << "SCALE " << b+1 << " N = " <<  N  << endl;
     NNP.TrainLWB = N;
     NNP.TrainUPB = LastTraining;
     NNP.TestUPB = LastEvaluate-1;
     NNP.TestLWB = NNP.TestUPB / 2;
     NNP.EvalLWB = FirstEvaluate;
     NNP.EvalUPB = LastEvaluate;
     NNP.alloc(N, Nx);
     NNErr = NNP.eval_predict(Scale.buffer(), TabScalePredY.buffer());
     for (i= FirstEvaluate; i <= LastEvaluate; i++)
                             TabSignalPred(i) += TabScalePredY(i);
     if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
   }
}

/***************************************************************************/

double MR_PRED::make_prediction(extrapol_type Method, fltarray &Signal,
                              fltarray &TabSignalPred, MR_1D & MR_Data, 
                              int NbrPredict, int User_AR_Order,
			      float & Inter)
{
   double Pred = 0.;
   int i,k,Nx=Signal.nx();
   TabSignalPred.alloc(NbrPredict+Nx);
   // for (i=0; i < Nx; i++) TabSignalPred(i) = Signal(i);
   TabSignalPred.init();
   int NParam=0;
   
   if (EvaluationMode == True)
   {
      // By default, separate the serie into two parts: on part for the training
      // and one part the evaluation
      if (LastTraining  <= 0) LastTraining = Nx/2 - 1;
      if (FirstEvaluate <= 0) FirstEvaluate  = Nx/2;
      if (LastEvaluate  <= 0) LastEvaluate = Nx -1;
      NbrTraining = LastTraining-FirstTraining+1;
      NbrEvaluate = LastEvaluate-FirstEvaluate+1;
      for (i=FirstEvaluate; i <= LastEvaluate; i++) TabSignalPred(i) = 0.;
   }
   else
   {
      if (LastTraining <= 0) LastTraining = Nx - 1;
      NbrTraining = LastTraining-FirstTraining+1;
      LastEvaluate = -1;
      NbrEvaluate = FirstEvaluate = 0;
   }
   TabAR_OrderPerScale.alloc(MAX(MR_Data.nbr_scale(),1));
  //  if (Method == EX_MR_AR_PRED1) NoLag = True;
   
   if (AR_AutoDetect == False)
   {
      if (User_AR_Order < 0) 
      {
           cout << "Error: User_AR_Order is not correctly set ... " << endl;
           exit(-1);
      }
      if (User_AR_Order > MaxNbrAR) MaxNbrAR = User_AR_Order;
   }
   
   if (Verbose == True) 
   {
      if (AR_AutoDetect == False) cout << "ARNbr  = " << User_AR_Order  << endl;
      cout << "Np   = " << Nx  << endl;
      cout << "NbrTraining    = " << NbrTraining   << endl;
      cout << "NbrEvaluate     = " <<  NbrEvaluate   << endl;
      if (AR_AutoDetect == True) cout << "Automatic AR order estimation " << endl;
      if (Stationary == True) cout << "Stationary signal " << endl;
      if (NoLag == True) cout << "Do not use lags " << endl;
      if (AR_AutoDetect == True)      
        switch (AR_Order_Detect)
        {
         case AR_ORDER_AIC: 
                cout << "  Use AIC method criterion " << endl; break;
         case AR_ORDER_AICC: 
                cout << "  Use AICC method criterion " << endl; break;
         case AR_ORDER_BIC: 
                 cout << "  Use BIC method criterion " << endl; break;
         default: cout << "Error: unknown criterion ... " << endl;
                  exit(-1);
         }
   }
   double NNErr=0.;
   int p,BestNbrAR=User_AR_Order;
   fltarray ARModel;
   int Nloop = MAX(1,NbrPredict);
   Bool MakeEval = EvaluationMode;
         
   switch (Method)
   {
      case EX_AR:
       {
         fltarray DataTend(Nx);
         fltarray TabTend(Nx);
	 DataTend = Signal;
	 if (UseTend == True)
	 {
	    mk_tendancy(Signal,Nx, Tend);
	    DataTend -= Tend;
	 }
	 
         for (p=0; p < Nloop; p++)
	 {
 	    PredDistance = p;
            if (AR_AutoDetect == True) 
                 get_best_ar_model(DataTend,  NbrTraining, BestNbrAR, ARModel);
            else get_ar_model(DataTend,  NbrTraining,  User_AR_Order, 0, ARModel);
	    NParam = BestNbrAR;
	    if (FirstTraining == 0)
	    {    
	        FirstTraining = BestNbrAR;
	        NbrTraining = LastTraining - FirstTraining + 1;
            }
	    // Evaluation prediction calculation
	    if (EvaluationMode == True)
	    {
              for (i=FirstTraining; i <= LastEvaluate; i++)
	      {
	         int Np = i+1;
	         if (UseTend == True) 
		 {
		    mk_tendancy(Signal,Np, TabTend);
		    for (k=0; k < Np; k++) DataTend(k) = Signal(k) - TabTend(k);
		 }
                 TabSignalPred(i) = ar_prediction(DataTend, BestNbrAR,i, ARModel);
                 if (UseTend == True)
 		    TabSignalPred(i) += 2.*TabTend(Np-1) - TabTend(Np-2-PredDistance);
               }
	    }
	    
	    // prediction
	    if (NbrPredict > 0)
	    {
	      TabSignalPred(Nx+p) = ar_prediction(DataTend, BestNbrAR, Nx, ARModel);
              if (UseTend == True)
 		 TabSignalPred(Nx+p) += 2.*Tend(Nx-1) - Tend(Nx-2-PredDistance);
	    }	    
            TabAR_OrderPerScale(0) = BestNbrAR;
            if (Verbose == True)
            {
	       // cout << " p = " << p << endl;
               cout << "AR order = " << BestNbrAR << " ARModel = ";
               for (i=0; i< BestNbrAR; i++)  cout << " " << ARModel(i); 
	       // cout << " Pred = " << TabSignalPred(Nx+p) << endl;
               cout << endl;  
            }
	    EvaluationMode = False;
	 }
	}
         break;
      case EX_MR_AR_PRED1:
         MR_Data.transform(Signal);
         for (p=0; p < Nloop; p++)
	 {
	    PredDistance = p;
            Pred = ar1_predict(MR_Data, TabSignalPred, User_AR_Order);
	    NParam = (int) TabAR_OrderPerScale.total();
	    if (NbrPredict > 0) TabSignalPred(Nx+p) = Pred;
	    EvaluationMode = False;
	 }
	 if ((Verbose == True) && (Nloop == 1))
	          cout << "Predicted value = " << Pred << endl;
         break;
      case  EX_MR_AR_PRED2:
         {
         fltarray DataTend(Signal.nx());
         tend_wt_transform(Signal, DataTend, Tend, MR_Data, Signal.nx(), True);
         for (p=0; p < Nloop; p++)
	 {
	    PredDistance = p;      
            Pred = ar2_predict(Signal, MR_Data, TabSignalPred, User_AR_Order);
	    NParam = (int) TabAR_OrderPerScale.total();
	    if (((Nx+p-1) < 0) || ((Nx+p-1) >= TabSignalPred.nx()))
	    {
	       cout << "Error: Nx = " << Nx << " p = " << p;
	       cout << " TabSignalPred.nx() = " << TabSignalPred.nx() << endl;
	       exit(-1);
	    }
	    if (NbrPredict > 0) TabSignalPred(Nx+p) = Pred;
	    // EvaluationMode = False;
	 }
	 if ((Verbose == True) && (Nloop == 1))
	          cout << "Predicted value = " << Pred << endl;
	 }	    
         break;
      case  EX_NN:
         {
            if (AR_AutoDetect == True) 
	    {
               get_best_ar_model(Signal,  NbrTraining, BestNbrAR, ARModel);
	       if (Verbose ==True)
	         cout << "Prediction from the  " << BestNbrAR << " values " << endl;
            }
	    NN_PREDICTION NNP;
            int N = (AR_AutoDetect == True) ? BestNbrAR: User_AR_Order;
            int FirstPos = N  + PredDistance;
	    NParam = N;
            if (FirstTraining < FirstPos) 
            {
               FirstTraining = FirstPos;
               NbrTraining = LastTraining - FirstTraining + 1;
            }
            NNP.TrainLWB = FirstTraining;
            NNP.TrainUPB = LastTraining;
            NNP.TestUPB = LastEvaluate-1;
            NNP.TestLWB = NNP.TestUPB / 2;
            NNP.EvalLWB = FirstTraining; // FirstEvaluate;
            NNP.EvalUPB = LastEvaluate;
            NNP.alloc(N, Nx);
            NNErr = NNP.eval_predict(Signal.buffer(),TabSignalPred.buffer(), NbrPredict);  
         }
         break;
      case EX_MR_NN:
         MR_Data.transform(Signal);
         nn1_predict(MR_Data, TabSignalPred, User_AR_Order);
         break;
      default:
         cout << "Error: not implemented extrapolation method ... " << endl;
         exit(-1);
   }
   EvaluationMode = MakeEval;
   
   // Error evaluation
   double ErrTrain=0.;
   double ErrPred=0.,Err,Pred0=0.;
   if (EvaluationMode == True)
   {
     int Cpt = 0;
     float Inter = 0.;
     for (i=FirstTraining; i <= LastTraining; i++)
     {
        Err = Signal(i) - TabSignalPred(i);
        ErrTrain += Err*Err;
      }
     ErrTrain /= (double) (NbrTraining-NParam);
     ErrTrain = sqrt(ErrTrain);
     
     for (i=FirstEvaluate; i <= LastEvaluate; i++)
     {
        Err = Signal(i) - TabSignalPred(i);
        ErrPred += Err*Err;
        Pred0 += Signal(i)*Signal(i);
     }
     ErrPred /= (double) NbrEvaluate;
     Pred0 /= (double) NbrEvaluate;
     ErrPred = sqrt(ErrPred);
     Pred0 = sqrt(Pred0);
     for (i=FirstEvaluate; i <= LastEvaluate; i++)
     {
        Err = Signal(i) - TabSignalPred(i);
        if (ABS(Err) <= 2.*ErrTrain) Cpt ++;
     }
     Inter = (float) Cpt / (float) NbrEvaluate * 100.;
     if (Verbose == True) 
     {
        cout << " FirstTraining = " << FirstTraining << " NbrTraining = " << NbrTraining << endl;
        cout << " FirstEvaluate = " << FirstEvaluate << " NbrEvaluate  = " << NbrEvaluate << endl;
        cout << " Nbr of param = " << NParam << endl;
        cout << "Pred Error for a 0 value prediction (i.e. max err.) = " << Pred0 << endl;
        cout << "Pred Error for the training part = " << ErrTrain << endl;
        cout << "Pred Error for the evaluation part = " << ErrPred << endl;
        cout << "Percentage in interval [TrainPred-2sigma,TrainPred+2sigma] = " << Inter << endl;
     }
   }
   return ErrPred;
}  

/***************************************************************************/
/***************************************************************************/

/* 
void get_ar_model_redundant(fltarray & Signal, int Np, int NbrAR, 
                  int MaxNbrTraining, 
                  fltarray & ARModel, 
                  Bool Verbose,  int Step, int NIter)
{
   int i,j,IndData;
   int t0 = Np-1;
   // int N = Signal.nx();
   int NbrTraining;
   int NFirstAR = Step - 1;
   NbrTraining = t0 - (NbrAR-NFirstAR)*Step;
   if ((MaxNbrTraining  > 0) && (NbrTraining >  MaxNbrTraining)) NbrTraining = MaxNbrTraining;

   if (NbrAR < Step) 
   {
       cout << "Error: in the redundant AR model, AR(order) must be > Step ... " << endl;
       exit(0);
   }
   if (NbrTraining < 0)
   {
       cout << "Error: Number of scales is to high ... " << endl;
       exit(-1);
   }
   // cout << "NbrTraining = " << NbrTraining << endl;
   // cout << " Step = " << Step  << "t0 = " << t0 << " NbrAR = " << NbrAR  << " NbrTraining = " << NbrTraining << endl;
 
   // Learning values: A X = B
   dblarray MAT_X(1,NbrAR);
   dblarray MAT_A(NbrAR, NbrTraining);
   dblarray InvMAT_A;
   dblarray MAT_B(1,NbrTraining);

   ARModel.alloc(NbrAR);
 
   MAT_X.init();
   for (i=0; i < NbrTraining; i++)
   {
      for (j=0; j < Step-1; j++)
      {
         IndData = t0 - i - (j+1);
         if (IndData < 0)
         {
            cout << endl;
            cout << "ERROR: IndData = " << IndData << endl;
            exit(-1);
         }
         MAT_A (j,i) = Signal(IndData);
      }
      for (j=1; j < NbrAR; j++)
      {
        IndData = t0 - i - (j+1)*Step;
        if (IndData < 0)
        {
            cout << endl;
            cout << "ERROR: IndData = " << IndData << endl;
            exit(-1);
        }
        MAT_A (j+Step,i) = Signal(IndData);
      }
      IndData = t0 - i;
      MAT_B(i) = Signal(IndData);
   } 
   lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
 
   for (j=0; j < NbrAR; j++) ARModel(j) = (float) MAT_X(j);
}
*/
