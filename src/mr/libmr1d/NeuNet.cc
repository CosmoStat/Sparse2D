


/******************************************************************************
                            D E C L A R A T I O N S
 ******************************************************************************/

#include "GlobalInc.h"
#include "NeuNet.h"

#define LO            0.1
#define HI            0.9
#define BIAS          1
#define sqr(x)        ((x)*(x))


/******************************************************************************
               A P P L I C A T I O N - S P E C I F I C   C O D E
 ******************************************************************************/

void NN_PREDICTION::NormalizeData()
{
  int   Pos;
        
  MinV = Data[0];
  MaxV = Data[0];
  for (Pos=1; Pos<NbrData; Pos++) 
  {
     MinV = MIN(MinV, Data[Pos]);
     MaxV = MAX(MaxV, Data[Pos]);
  }
  Mean = 0;
  for (Pos=0; Pos<NbrData; Pos++) 
  {
    PredData[Pos] = Data [Pos] = ((Data[Pos]-MinV) / (MaxV-MinV)) * (HI-LO) + LO;
    Mean += Data[Pos] / NbrData;
  }
}


void NN_PREDICTION::InitializeApplication()
{
  int   Pos, i;
  double Out, Err;

  Alpha = 0.5;
  Eta   = 0.05;
  Gain  = 1;

  NormalizeData();
  TrainErrorPredictingMean = 0;
  for (Pos=TrainLWB; Pos<=TrainUPB; Pos++) 
  {
    for (i=0; i<M; i++) 
    {
      Out = Data[Pos+i];
      Err = Mean - Out;
      TrainErrorPredictingMean += 0.5 * sqr(Err);
    }
  }
  TestErrorPredictingMean = 0;
  for (Pos= TestLWB; Pos<= TestUPB; Pos++) 
  {
    for (i=0; i<M; i++) 
    {
      Out = Data[Pos+i];
      Err = Mean - Out;
      TestErrorPredictingMean += 0.5 * sqr(Err);
    }
  }
}


/******************************************************************************
                          I N I T I A L I Z A T I O N
 ******************************************************************************/


void NN_PREDICTION::GenerateNetwork()
{
  int  l,i;

  Layer = (LAYER**) calloc(NbrLayers, sizeof(LAYER*));
   
  for (l=0; l<NbrLayers; l++) 
  {
    Layer[l] = (LAYER*) malloc(sizeof(LAYER));
      
    Layer[l]->Units      = Units[l];
    Layer[l]->Output     = (double*)  calloc(Units[l]+1, sizeof(double));
    Layer[l]->Error      = (double*)  calloc(Units[l]+1, sizeof(double));
    Layer[l]->Weight     = (double**) calloc(Units[l]+1, sizeof(double*));
    Layer[l]->WeightSave = (double**) calloc(Units[l]+1, sizeof(double*));
    Layer[l]->dWeight    = (double**) calloc(Units[l]+1, sizeof(double*));
    Layer[l]->Output[0]  = BIAS;
      
    if (l != 0) 
    {
      for (i=1; i<=Units[l]; i++) 
      {
        Layer[l]->Weight[i]     = (double*) calloc(Units[l-1]+1, sizeof(double));
        Layer[l]->WeightSave[i] = (double*) calloc(Units[l-1]+1, sizeof(double));
        Layer[l]->dWeight[i]    = (double*) calloc(Units[l-1]+1, sizeof(double));
      }
    }
  }
  InputLayer  = Layer[0];
  OutputLayer = Layer[NbrLayers - 1];
  Alpha       = 0.9;
  Eta         = 0.25;
  Gain        = 1;
}


void NN_PREDICTION::RandomWeights()
{
  int  l,i,j;
   
  for (l=1; l<NbrLayers; l++) 
  {
    for (i=1; i<=Layer[l]->Units; i++) 
    {
      for (j=0; j<=Layer[l-1]->Units; j++) 
      {
        Layer[l]->Weight[i][j] = RandomEqualdouble(-0.5, 0.5);
      }
    }
  }
}


void NN_PREDICTION::SetInput(double* Input)
{
  int  i;
   
  for (i=1; i<=InputLayer->Units; i++) 
  {
    InputLayer->Output[i] = Input[i-1];
  }
}


void NN_PREDICTION::GetOutput(double* Output)
{
  int  i;
   
  for (i=1; i<= OutputLayer->Units; i++) {
    Output[i-1] = OutputLayer->Output[i];
  }
}


/******************************************************************************
            S U P P O R T   F O R   S T O P P E D   T R A I N I N G
 ******************************************************************************/


void NN_PREDICTION::SaveWeights()
{
  int  l,i,j;

  for (l=1; l<NbrLayers; l++) 
  {
    for (i=1; i<=Layer[l]->Units; i++) 
    {
      for (j=0; j<=Layer[l-1]->Units; j++) 
      {
        Layer[l]->WeightSave[i][j] = Layer[l]->Weight[i][j];
      }
    }
  }
}


void NN_PREDICTION::RestoreWeights()
{
  int  l,i,j;

  for (l=1; l<NbrLayers; l++) 
  {
    for (i=1; i<=Layer[l]->Units; i++) 
    {
      for (j=0; j<=Layer[l-1]->Units; j++) 
      {
         Layer[l]->Weight[i][j] = Layer[l]->WeightSave[i][j];
      }
    }
  }
}


/******************************************************************************
                     P R O P A G A T I N G   S I G N A L S
 ******************************************************************************/


void NN_PREDICTION::PropagateLayer(LAYER* Lower, LAYER* Upper)
{
  int   i,j;
  double Sum;

  for (i=1; i<=Upper->Units; i++) 
  {
    Sum = 0;
    for (j=0; j<=Lower->Units; j++) 
    {
      Sum += Upper->Weight[i][j] * Lower->Output[j];
    }
    Upper->Output[i] = 1 / (1 + exp(-Gain * Sum));
  }
}


void NN_PREDICTION::PropagateNet()
{
  int  l;
   
  for (l=0; l<NbrLayers-1; l++) 
  {
    PropagateLayer(Layer[l], Layer[l+1]);
  }
}


/******************************************************************************
                  B A C K P R O P A G A T I N G   E R R O R S
 ******************************************************************************/


void NN_PREDICTION::ComputeOutputError(double* Target)
{
  int   i;
  double Out, Err;
   
  Error = 0;
  for (i=1; i<=OutputLayer->Units; i++) 
  {
    Out = OutputLayer->Output[i];
    Err = Target[i-1]-Out;
    OutputLayer->Error[i] = Gain * Out * (1-Out) * Err;
    Error += 0.5 * sqr(Err);
  }
}


void NN_PREDICTION::BackpropagateLayer(LAYER* Upper, LAYER* Lower)
{
  int   i,j;
  double Out, Err;
   
  for (i=1; i<=Lower->Units; i++) 
  {
    Out = Lower->Output[i];
    Err = 0;
    for (j=1; j<=Upper->Units; j++) 
    {
      Err += Upper->Weight[j][i] * Upper->Error[j];
    }
    Lower->Error[i] = Gain * Out * (1-Out) * Err;
  }
}


void NN_PREDICTION::BackpropagateNet()
{
  int  l;
   
  for (l=NbrLayers-1; l>1; l--) 
  {
    BackpropagateLayer(Layer[l], Layer[l-1]);
  }
}


void NN_PREDICTION::AdjustWeights()
{
  int   l,i,j;
  double Out, Err, dWeight;
   
  for (l=1; l<NbrLayers; l++) 
  {
    for (i=1; i<=Layer[l]->Units; i++)
    {
      for (j=0; j<=Layer[l-1]->Units; j++) 
      {
        Out = Layer[l-1]->Output[j];
        Err = Layer[l]->Error[i];
        dWeight = Layer[l]->dWeight[i][j];
        Layer[l]->Weight[i][j] += Eta * Err * Out + Alpha * dWeight;
        Layer[l]->dWeight[i][j] = Eta * Err * Out;
      }
    }
  }
}


/******************************************************************************
                      S I M U L A T I N G   T H E   N E T
 ******************************************************************************/


void NN_PREDICTION::SimulateNet(double* Input, double* Output, double* Target, Bool Training)
{
  SetInput(Input);
  PropagateNet();
  GetOutput(Output);
   
  ComputeOutputError(Target);
  if (Training == True) 
  {
    BackpropagateNet();
    AdjustWeights();
  }
}


void NN_PREDICTION::TrainNet(int Epochs)
{
  int   Pos, n;
  double* Output = new double[M];

  for (n=0; n<Epochs*NbrTraining; n++) 
  {
     Pos = RandomEqualint (TrainLWB, TrainUPB);
     SimulateNet(&(Data[Pos-N]), Output, &(Data[Pos]), True);
  }
}


void NN_PREDICTION::TestNet()
{
  int   Pos;
  double* Output = new double[M];

  TrainError = 0;
  for (Pos=TrainLWB; Pos<=TrainUPB; Pos++) 
  {
    SimulateNet(&(Data[Pos-N]), Output, &(Data[Pos]), False);
    TrainError += Error;
  }
  TestError = 0;
  for (Pos= TestLWB; Pos<= TestUPB; Pos++) 
  {
    SimulateNet(&(Data[Pos-N]), Output, &(Data[Pos]), False);
    TestError += Error;
  }
  if (Verbose == True)
    printf("\nNMSE is %0.3f on Training Set and %0.3f on Test Set",
             TrainError / TrainErrorPredictingMean,
             TestError / TestErrorPredictingMean);
}


void NN_PREDICTION::EvaluateNet()
{
  int   Pos;
  double* Output = new double[M];
  double* Output_ = new double[M];

  if (Verbose == True)
  {
     printf("\n\n\n");
     printf("Pos    Data    Open-Loop Prediction    Closed-Loop Prediction\n");
     printf("\n");
  }
  for (Pos=EvalLWB; Pos<=EvalUPB; Pos++) 
  {
    SimulateNet(&(Data [Pos-N]), Output,  &(Data [Pos]), False);
    SimulateNet(&(PredData[Pos-N]), Output_, &(PredData[Pos]), False);
    PredData[Pos] = Output_[0];
    if (Verbose == True)
      printf("%d       %0.3f                   %0.3f                     %0.3f\n",
               Pos+1,
               Data[Pos],
               Output [0],
               Output_[0]);
  }
}

void NN_PREDICTION::EvaluatePred(float *Pred)
{
  int   Pos;
  double* Output = new double[M];

  for (Pos=EvalLWB; Pos<=EvalUPB; Pos++) 
  {
    SimulateNet(&(Data [Pos-N]), Output,  &(Data [Pos]), False);
    Pred[Pos] = (float) Output[0];
  }
}

void NN_PREDICTION::MakePred(float *Pred, int NbrPredict)
{
  int   Pos;
  double* Output = new double[M];
  int NbrTot = NbrData+NbrPredict;
  double *P = new double [NbrTot];
  
  for (Pos=0; Pos<NbrData; Pos++) P[Pos] = Data[Pos];
 
  for (Pos=NbrData; Pos<NbrData+NbrPredict; Pos++) 
  {
     P[Pos]  = 0;
     SimulateNet(&(P[Pos-N]), Output, &(P[Pos]), False);
     P[Pos]  = (float) Output[0];
  }
  for (Pos=EvalLWB; Pos<=EvalUPB; Pos++) 
  {
    SimulateNet(&(Data [Pos-N]), Output,  &(Data [Pos]), False);
    P[Pos] = (float) Output[0];
  }
    
  double Coef = (HI-LO) / (MaxV-MinV);
  for (Pos=0; Pos<NbrData+NbrPredict; Pos++) 
       Pred[Pos]  = (float) ((P[Pos] - LO) / Coef + MinV);
       
  delete [] P;
}

/******************************************************************************
                                    M A I N
******************************************************************************/

double NN_PREDICTION::eval_predict(float *Input, float *Pred, int NbrPredict)
{
  int i;
  Bool Stop;
  double MinTestError;

  for (i=0; i < NbrData; i++) Data[i] = Input[i];

  InitializeRandoms();
  GenerateNetwork();
  RandomWeights();
  InitializeApplication();

  Stop = False;
  MinTestError = Maxfloat;
  do  
  {
    TrainNet(10);
    TestNet();
    if (TestError < MinTestError) 
    {
      // printf(" - saving Weights ...");
      MinTestError = TestError;
      SaveWeights();
    }
    else if (TestError > 1.2 * MinTestError) 
    {
      // printf(" - stopping Training and restoring Weights ...");
      Stop = True;
      RestoreWeights();
    }
  } while (Stop == False);
  TestNet();
  // EvaluateNet();
  // EvaluatePred(Pred);
  MakePred(Pred, NbrPredict);

  double Err=0.;
  for (int Pos=EvalLWB; Pos<=EvalUPB; Pos++) 
  {
     Err += (Input[Pos] - Pred[Pos])*(Input[Pos] - Pred[Pos]);
  }
  if (NbrEval > 0)
  {
     Err /= (double) NbrEval;
     Err = sqrt(Err);
  }
  return Err;
}
