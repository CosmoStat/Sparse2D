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
**    this module contains the routines for neural network prediction
**    http://www.geocities.com/CapeCanaveral/1624/bpn.html
**
**************************************************************************/ 

#ifndef __MR1D_NENET__
#define __MR1D_NENET__
 

/***************************************************************************/

typedef struct {                     /* A LAYER OF A NET:                     */
        int           Units;         /* - number of units in this layer       */
        double*         Output;        /* - output of ith unit                  */
        double *         Error;         /* - error term of ith unit              */
        double **        Weight;        /* - connection weights to ith unit      */
        double **        WeightSave;    /* - saved weights for stopped training  */
        double **        dWeight;       /* - last weight deltas for momentum     */
} LAYER;


/***************************************************************************/

class NN_PREDICTION
{
   LAYER**       Layer;         /* - layers of this net                  */
   LAYER*        InputLayer;    /* - input layer                         */
   LAYER*        OutputLayer;   /* - output layer                        */
   double          Alpha;         /* - momentum factor                     */
   double          Eta;           /* - learning rate                       */
   double          Gain;          /* - gain of sigmoid function            */
   double          Error;         /* - total net error                     */

  int NbrLayers;
  int *Units; //  = {N, 10, M};
  int N;
  int M;
  int NbrData;
  int NbrTraining;
  int NbrTest;
  int NbrEval;

  double  Mean;
  double  TrainError;
  double  TrainErrorPredictingMean;
  double  TestError;
  double  TestErrorPredictingMean;
  double *Data;
  double *PredData;
  double MinV, MaxV;

  void InitializeRandoms() {srand(4711);}
  int RandomEqualint(int Low, int High) {return rand() % (High-Low+1)+Low;}      
  double RandomEqualdouble(double Low, double High)
                   {return ((double) rand() / RAND_MAX) * (High-Low) + Low;}      

  void NormalizeData();
  void InitializeApplication();
  void GenerateNetwork();
  void RandomWeights();
  void SetInput(double* Input);
  void GetOutput(double* Output);
  void SaveWeights();
  void RestoreWeights();
  void PropagateLayer(LAYER* Lower, LAYER* Upper);
  void PropagateNet();
  void ComputeOutputError(double* Target);
  void BackpropagateLayer(LAYER* Upper, LAYER* Lower);
  void  BackpropagateNet();
  void AdjustWeights();
  void SimulateNet(double* Input, double* Output, double* Target, Bool Training);
  void TrainNet(int Epochs);
  void TestNet();
  void EvaluateNet();
  void EvaluatePred(float *Pred);
  void MakePred(float *Pred, int NbrPredict);

  public:
   Bool Verbose;       // if true, print information
   int TrainLWB;
   int TrainUPB;
   int TestLWB;
   int TestUPB;
   int EvalLWB;
   int EvalUPB;
   NN_PREDICTION(){Verbose=False;NbrLayers=3;N=30;M=1;
                     Units = new int [NbrLayers]; Units[0] = N;
                     Units[1] = 10; Units[2] = M;Data=NULL;PredData=NULL;}
   void alloc(int Nt, int Np) {NbrLayers=3;N=Nt;M=1;
                          Units = new int [NbrLayers]; Units[0] = N;
                          Units[1] = 10; Units[2] = M;NbrData= Np;
                          TrainLWB = N; 
                          if (TrainUPB == 0) TrainUPB = NbrData/3;
                          NbrTraining = (TrainUPB - TrainLWB + 1);
                          NbrTest  =  ( TestUPB -  TestLWB + 1);
                          NbrEval  =  (EvalUPB - EvalLWB + 1);
                          PredData = new double [NbrData];
                          Data = new double [NbrData];}

  double eval_predict(float *Input, float *Pred, int NbrPredict=0);
  // return the prediction error on the evaluation part
  // Input: IN = input signal
  // Pred: Out = out predicted signal
  //       Pred must be allocated with a size of Np+NbrPredict
  //       Pred[0..EvalLWB-1] == Input[0..EvalLWB-1]
  //       Pred[i] i in [EvalLWB..EvalUPB] == prediction using Input[0..i-1]
  //       Pred[i] i >= Np ==> prediction from the last date and the last
  //                    prediction

   ~NN_PREDICTION() {if (Data != NULL) delete [] Data;
                     if (PredData != NULL) delete []  PredData;
                     if (Units != NULL) delete [] Units;
                    }
};

/***************************************************************************/
 

#endif
