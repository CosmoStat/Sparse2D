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
**    Date:  20/07/2000
**    
**    File:  MR1D_Predict.h
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

#ifndef __MR1D_EXTRAPOL__
#define __MR1D_EXTRAPOL__

#define NBR_EX_METHOD 4 // NN abd AR not yet implemented

enum  extrapol_type {EX_AR, EX_MR_AR_PRED1, EX_MR_AR_PRED2, EX_NN, EX_MR_NN, EX_POLYNOM};
              // EX_POLYNOM ==> Polynomial prediction. Does not work correctly! 
              // EX_NN      ==>  Neural network prediction. Not implemented
              // EX_AR      ==>  Autoregressive model. 
              //                 predict x(t) from  x(t-j)     
              // EX_MR_PRED1  ==> predict the wavelet coefficient w(t) of
              //                x(t) from the wavelet coefficient w(t-j)
              //                or w(t-2^j) if lags are considered
              // EX_MR_PRED2  ==> predict x(t) from the wavelet coefficient
              //                w(t-j) or w(t-2^j) if lags are considered


char *StringExMethod (extrapol_type type);
void extrap_usage(extrapol_type Method);


#define NBR_PRED_BORD 5
enum  pred_mirror_type {P_CONT, P_MIRROR, P_MIRROR_2, P_POL1, P_POL2};
inline char * StringPredBord (pred_mirror_type type)
{
    switch (type)
    {
        case P_CONT: 
              return ("Constant border");break;
        case P_MIRROR: 
              return ("Mirror border");break;
        case P_MIRROR_2: 
              return ("Double mirror border");break;
        case P_POL1: 
              return ("Polynomial extrapolation (deg 1)");break;
        case P_POL2: 
              return ("Polynomial extrapolation (deg 2)");break;
        default:
              return ("Undefined border type");
              break;
    }
}
inline void pred_border_usage(pred_mirror_type type)
{
  int i;
    fprintf(OUTMAN, "        [-B extrapol_type]\n");
    for (i = 0; i < NBR_PRED_BORD; i++)
    {
       pred_mirror_type T = (pred_mirror_type) i;
       fprintf(OUTMAN, "              %d: %s \n",i+1, StringPredBord(T));
    }
    fprintf(OUTMAN, "             default is %s\n", StringPredBord(type));
    manline();
}
 
/***************************************************************************/
class MR_FIL_PRED;  

class MR_PRED: public AR_PREDICTION
{        
   friend class MR_FIL_PRED;

   void get_mr_ar1_model(MR_1D & MR_Data, int Np, int User_AR_Order=-1);
   // set TabAR_OrderPerScale and TabAR_ModelPerScale
   // considering the lags.
   // if AR_AutoDetect == True calculate the AR order per scale
   //                       AR order are always calculated using lags
   //                       in the scale
   // if Stationary == True, the last smooth scale is also analyzed
   // MR_Data: IN = wavelet data
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // User_AR_Order: IN = AR order to use (only if AR_AutoDetect == False)

   void get_mr_ar2_model(fltarray & Signal, MR_1D & MR_Data, int Np,   
                         int NbrTotAR, fltarray & ARModel);
   // Get the AR model to predict x(t) from the wavelet coefficients.
   // MR_Data: IN = wavelet data
   // Signal: IN = input signal
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // NbrTotAR: IN = total number of coefficients used for the prediction
   //                NbrTotAR = sum_i TabAR_OrderPerScale(i)
   // ARModel: OUT = AR model. ARModel(0..NbrTotAR-1)


   int FirstTraining; // set to 0 (the N first pixels are used for the training
   int NbrTraining;  // Number of pixels used for the AR estimation
                     // if NbrLearning = 0, all points are used   
   int NbrEvaluate;  // Number of pixels used for the evaluation of 
                     // the prediction
   dblarray TabAR_ModelPerScale; // AR model per scale

   double ar1_predict(MR_1D & MR_Data, fltarray &TabSignalPred,int User_AR_Order=-1);
                    // make the prediction of the wavelet coefficient scale
                    // per scale. The signal prediction is the reconstruction
                    // of all scales.
		    // if LastEvaluate-FirstEvaluate > 0
		    // then put in TabSignalPred the prediction values
		    // return the prediction.

   double ar2_predict(fltarray & Signal, MR_1D & MR_Data, 
                    fltarray &TabSignalPred, int User_AR_Order=-1);
                    // make the prediction of the signal from its wavelet
		    // coefficients
		    // return the prediction.
		    // if LastEvaluate-FirstEvaluate > 0
		    // then put in TabSignalPred the prediction values
		    // at position range from [FirstEvaluate,LastEvaluate]
		    
   void nn1_predict(MR_1D & MR_Data,  fltarray &TabSignalPred, int User_AR_Order);

   double extr_last_scale(fltarray & Scale, int Np, int DistPred);
                   // Extrapolation of the last scale.
		   // Scale = Input Data (0..Np-1)
		   // Np = number of pixels 
		   // DistPred = Distance to predict, i.e. prediction at
		   // at pixel position:  Np+DistPred
      
public:
   
   MR_PRED() {NoLag=False;Verbose=False;Stationary=False;MaxNbrAR=10;
            AR_AutoDetect=True;AR_Order_Detect=AR_ORDER_BIC;
            NbrTraining=0;FirstTraining=0;NbrTraining=0;NbrEvaluate=0;
            LastTraining=-1;FirstEvaluate=0;LastEvaluate=-1;
	    NbrPixLastPred=5;NbrPixTendancy=100;
	    EvaluationMode=False;Debug=False;MinAR_OrderPerScale=1;
	    TPredBord=P_MIRROR;UseB3SplineWT=False;UseTend=False;}
   Bool UseTend;        // If true, a tendancy is calculated using polynomes
                        // of degree 2.
   fltarray Tend;       // Tendancy array. Only used if UseTend == True			
   Bool UseB3SplineWT; // If true, then a B3 spline is used instead of the
                       // Haar WT    
   pred_mirror_type TPredBord;
                       // Type of border to be used with B3-spline WT
                       // Extrapolation type for the last scale
   Bool NoLag;         // if true, do not use any "hole" between pixels
   Bool Verbose;       // if true, print information
   Bool Debug;
   Bool Stationary;    // if true, make a stationary assumption
   Bool AR_AutoDetect; // if true, detect automatically the AR order
   int MaxNbrAR;       // AR maximum order value. Default is 10

   int LastTraining;   // Last pixel position for the training
   int FirstEvaluate;  // First pixel position for the evaluation
   int LastEvaluate;   // Last pixel position for the evaluation
   intarray TabAR_OrderPerScale; // Order number per scale
                                 // For a simple AR, TabAR_OrderPerScale(0)
                                 // is the used AR order
   int MinAR_OrderPerScale; // In the automatic AR estimation at each scale, we
                            // can force the AR order to be at least equal to
			    // this value. Default is 1.
   int NbrPixLastPred;  
                        // Number of pixels used for the prediction of
                        // the last scale (by a polynome of degree 2)
			// Only used in case of non stationary signal.
			// def. value is 5
   int NbrPixTendancy;  // Number of pixels used for the prediction of
                        // the tendancy by a polynome of degree 2
   Bool EvaluationMode; // If true, then a part of signal is used for
                        // prediction evaluation
			// if false, no evaluation and the whole signal
			// is used for the prediction
   double make_prediction(extrapol_type Method, fltarray &Signal,
                              fltarray &TabSignalPred, MR_1D & MR_Data, 
                              int NbrPredict, int User_AR_Order,
			      float & Inter);
			// Evaluate a prediction method
			// Method: IN = prediction method
			// Signal: IN: = input signal
			// TabSignalPred: OUT = output signal
			//    TabSignalPred(i) are the predicted values
			//      i in [FirstEvaluate .. LastEvaluate]
			//      if i not in [FirstEvaluate .. LastEvaluate]
			//      then TabSignalPred(i)  = Signal(i)
			//    pixels from 0 to LastTraining are used for
			//    estimation of prediction parameters
			//   Calculate TabAR_OrderPerScale(j) as the
			//   best AR order at scale j 
			// Inter: Out = number of pixels in the interval
			//        [Pred-2Sigma, Pred+2sigma]
   void wave_1d_trans(fltarray &Signal, MR_1D & MR_Data, int N);
   void mk_tendancy(fltarray & Data, int Np, fltarray & TendTab);
   void tend_wt_transform(fltarray & Data, 
                          fltarray & DataTend, fltarray & TabTend,
			  MR_1D & MR_Data, int Np, Bool ForceWT=False);
  ~MR_PRED(){}
};

void autocor(fltarray  & Data, fltarray &TabAuto, int Np, int NBShift);
// Calculate the autocorrelation function of a signal
// on the Np first points, with NBShift: TabAuto(0..NBShift-1): 

void tendancy_est(fltarray & Data, int Np, fltarray & Tend, int NbrPixTendancy);
// Calculate the tendancy on the Np first pixels of the signal Data
// The tendancy is estimated using a weighted polynome of degree 2 in
// a window of size NbrPixTendancy around each pixel.

double pred_pol_deg2(fltarray &Data, int Np, int NLast, int DistPred);
// Calcule the polynome of degree 2 fitting the Nlast points
// of the data
// Data(i) = a x^2 + b x + 1    
//           and i = [Np-NLast,Np-1]
// Np: IN = Number of pixels in Data
// Data: IN = input data
// Nlast: IN = number of points in Data to be used
// DistPred: IN = Calculate the prediction at a distance DistPred from
//                the last pixel (DistPred=1) for the next pixel prediction)
//                and return this value   

#endif
