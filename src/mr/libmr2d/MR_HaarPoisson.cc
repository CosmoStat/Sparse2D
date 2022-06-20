/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  1.03.99
**    
**    File:  MR_HaarPoisson.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for Poisson noise management with Haar transform
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "MR_Obj.h"

#include "MR_HaarPoisson.h"

/****************************************************************************/
double TabHaarEps[NBR_EPS]={1e-6,1e-5,1e-4,1e-3,1e-2};

double TabHaarThresholdPoisson[NBR_LAMBA][NBR_EPS] = 
{ 
{1,1,1,1,1},  // L = 0
{1,1,1,1,1},  // L = 2^-30 
{1,1,1,1,1},  // L = 2^-29 
{1,1,1,1,1},  // L = 2^-28
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},
{1,1,1,1,1},  // L = 2^-19 
{2,1,1,1,1},
{2,1,1,1,1},
{2,1,1,1,1},
{2,2,1,1,1},
{2,2,1,1,1},
{2,2,1,1,1},
{2,2,2,1,1},
{2,2,2,1,1},
{2,2,2,1,1},
{2,2,2,1,1},
{3,2,2,2,1},
{3,2,2,2,1},
{3,3,2,2,1}, 	
{3,3,3,2,2},	
{4,3,3,2,2},
{4,4,3,3,2},
{5,4,5,3,2},  // L =2^-2
{6,5,5,4,3},  // L =2^-1
{7,6,6,5,3},  // L =2^-0  = 1
{9, 8, 7, 6, 4},  // L =2^1
{12,10, 9, 7, 6},
{15,  14,  12,  10, 8},
{21,  19,  16,  14,  10},
{28,  26,  22,  19,  14},
{39,  35,  31,  26,  20},
{55,  49,  43,  36,  27},
{77,  69,  61,  50,  38},
{109,  98,  85,  71,  54},
{153, 138, 120, 100,  75},
{216, 194, 169, 141, 106},
{305, 274, 239, 199, 150},
{431, 387, 338, 281, 212},
{609, 547, 477, 397, 299},
{861, 773, 674, 560, 422},
{1217, 1092, 953, 792, 596},
{1721, 1545, 1347, 1119, 843},
{2434, 2184, 1905, 1583, 1192},
{3442, 3089, 2693, 2238, 1685},
{4868, 4368, 3809, 3165, 2383},
{6884, 6177, 5386, 4476, 3369},
{9736, 8735, 7617, 6329, 4765},
{13768, 12353, 10772, 8951, 6738},
{19471, 17469, 15234, 12658, 9529},
{27535, 24705, 21543, 17901, 13476},
{38941, 34938, 30467, 25316, 19058},
{55070, 49410, 43086, 35802, 26952},
{77881, 69876, 60933, 50631, 38115},
{110140, 98820, 86172, 71603, 53903},
{155761, 139752, 121865, 101261, 76230}// L =2^30
};

/****************************************************************************/


static double linear_interpol(double x1, double x2, double y1, double y2, double x)
{
   double Delta = ABS(x2-x1);
   if (Delta == 0) return y1;
   else
   {
       double w = ABS(x2 - x) / Delta;
       return (y1 * w + y2 * (1-w));
   }
}

/****************************************************************************/
    
double get_harr_poisson_threshold(double Lambda, double Eps)
{
   if ((Eps < TabHaarEps[0]) || (Eps > TabHaarEps[NBR_EPS-1]))
   {
      cerr << "Error: Epsilon value not pre-computed ... " << endl;
      cerr << "   Eps = " << Eps << endl;
      exit(-1);
   }
   int Eb=0;
   int Ea=0;
   double ValRet,ValEpsBef,ValEpsAfter;
   
   while (Eps > TabHaarEps[Ea]) Ea ++;
   Eb = Ea - 1;
   ValEpsBef = TabHaarEps[Eb];
   ValEpsAfter = TabHaarEps[Ea];
   // we get  ValEpsBef < Eps < ValEpsAfter
   
   if (Lambda < 0)
   {
      cerr << "Error: Lamda value must be > 0 ... " << endl;
      exit(-1);
   }
   else if (Lambda > HaarMaxExp)
   {
      cerr << "Error: Lamda value not pre-computed ... " << endl;
      exit(-1);
   }
   
   if (Lambda <= HaarInfExp) ValRet = 1;
   else
   {
      // the table is defined between 2^-30 and 2^30
      int Fi = -30;
      double ValInd = log(Lambda) / log(2.);
      double ValLambdaBef, ValLambdaAfter;
      int Lb=0;
      int La=0;
      if (ValInd < Fi) 
      {
         La = 1;
         ValLambdaBef = 0;
	 ValLambdaAfter = HaarMinExp;
      }
      else
      {
         Lb = (ValInd < 0) ? (int) ValInd - 1: (int) ValInd;
	 La = (ValInd < 0) ? (int) ValInd: (int) ValInd + 1;
	 if (La > 30) La = 30;
	 ValLambdaBef = POW2((double) Lb);
	 ValLambdaAfter = POW2((double) La);
	 Lb = Lb - Fi + 1;
	 La = La - Fi + 1;
      }
      double V1 = linear_interpol(ValLambdaBef, ValLambdaAfter, 
                     TabHaarThresholdPoisson[Lb][Eb], 
		     TabHaarThresholdPoisson[La][Eb], Lambda);
      double V2 = linear_interpol(ValLambdaBef, ValLambdaAfter, 
                     TabHaarThresholdPoisson[Lb][Ea], 
		     TabHaarThresholdPoisson[La][Ea], Lambda);		    
      ValRet = linear_interpol((double) ValEpsBef, (double) ValEpsAfter,  
	                       V1, V2, Eps);
    }
    return ValRet;
}

/****************************************************************************/
