/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  SB_Filter1D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  1D Sub-Band decomposition
**    -----------  
**                 
******************************************************************************/
 

#include "MSVST_SB_Filter1D.h"

#define FIL_DEBUG 0


static void filter_complement(float *Filter, float *ComFilter, int N, int Sign)
{
    int i, Sgn=Sign;
    for (i = 0; i < N; i++)
    {
        ComFilter[i] = Sgn * Filter[i];
        Sgn *= -1;
    }       
}

/***********************************************************************/

void SubBandFilter::reset_param()
{
   H0 = G0 = H1 = G1 = NULL;
   Size_H0 = Size_H1 = Size_G0 = Size_G1 = 0;
   Start_H0 = Start_H1 = Start_G0 = Start_G1 =0;   
   TypeFilter = SB_UNKNOWN;
   NormCoef = 0;
   Verbose = False;

}

/***********************************************************************/

SubBandFilter::SubBandFilter()
{
   reset_param();
   FilterAnaSynt FILTER(F_HAAR);
   init(FILTER, DEF_SB_NORM);
}

SubBandFilter::SubBandFilter(const SubBandFilter &sbf):SubBand1D(sbf)
{
   reset_param();
   FilterAnaSynt FILTER(sbf.TypeFilter);
   init(FILTER, sbf.TypeNorm);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(FilterAnaSynt & FILTER)
{
   reset_param();
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(FilterAnaSynt & FILTER, sb_type_norm Norm)
{
   reset_param();
   init(FILTER, Norm);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(type_sb_filter T_Filter)
{
   reset_param();
   FilterAnaSynt FILTER(T_Filter);
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(type_sb_filter T_Filter, sb_type_norm Norm)
{
   reset_param();
   FilterAnaSynt FILTER(T_Filter);
   init(FILTER, Norm);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(char *FileName)
{
   reset_param();
   FilterAnaSynt FILTER(FileName);
   init(FILTER, DEF_SB_NORM);
}

/***********************************************************************/

SubBandFilter::SubBandFilter(char *FileName, sb_type_norm Norm)
{
   reset_param();
   FilterAnaSynt FILTER(FileName);
   init(FILTER, Norm);
}


/***********************************************************************/
 
SubBandFilter::SubBandFilter(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm)
{
   sb_type_norm FilterNorm = NORM_L1; // filters h and g are supposed to be L1 normalized
   int i;
   TypeNorm = FilNorm;
   NormCoef = (TypeNorm == NORM_L1) ? 2 : 1;
      
   Size_H0 = Nh;
   Size_G0 = Ng;
   Size_H1 = Size_G0;
   Size_G1 = Size_H0;
   H0 = new float[Size_H0];
   H1 = new float[Size_H1];
   G1 = new float[Size_G1];
   G0 = new float[Size_G0];
   float NormCorrection = 1.;
   if (FilterNorm == TypeNorm) NormCorrection = 1.;
   else if (TypeNorm == NORM_L1)  NormCorrection = 1. / sqrt(2.);
   else if (TypeNorm == NORM_L2)  NormCorrection = sqrt(2.);
   for (i=0; i < Size_H0; i++) H0[i] = H[i]*NormCorrection;
   for (i=0; i < Size_G0; i++) G0[i] = G[i]*NormCorrection;
    
   Start_H0 = - Size_H0 / 2;
   Start_H1 = - (Size_H1-1) / 2;
   Start_G0 = - Size_G0 / 2;
   Start_G1 = - (Size_G1-1) / 2; 

   filter_complement(H0, G1, Size_H0, 1);
   filter_complement(G0, H1, Size_G0, -1);  

//    float T=0;
//    for (i=0; i < Size_H0; i++) T += H0[i];
//    printf("H0 = %f\n", T);
//    for (i=0; i < Size_H0; i++) printf("%f\n", H0[i]);
//
//    T=0;
//    for (i=0; i < Size_G0; i++) T+= G0[i];
//    printf("G0 = %f\n", T);
//    for (i=0; i < Size_G0; i++) printf("%f\n", G0[i]);
//    
//    T=0;
//    for (i=0; i < Size_H1; i++) T+= H1[i];
//    printf("H1 = %f\n", T);
//    for (i=0; i < Size_H1; i++) printf("%f\n", H1[i]);
//    
//    T=0;
//    for (i=0; i < Size_G1; i++) T+= G1[i];
//    printf("G1 = %f\n", T);
//    for (i=0; i < Size_G1; i++) printf("%f\n", G1[i]);
}

/***********************************************************************/

void SubBandFilter::init(FilterAnaSynt & FILTER, sb_type_norm Norm)
{
    if ((FILTER.type_filter() == U_HAAR_B3S) || (FILTER.type_filter() == U_HAAR_B3S2))
    {
        TypeNorm = Norm;
        NormCoef = (TypeNorm == NORM_L1) ? 2 : 1;
	
 	    TypeFilter = FILTER.type_filter();
	    float *F_H0 = FILTER.analysis();
	    float *F_H1 = FILTER.synthesis();
	    float *F_G0 = FILTER.analysis2();
	    float *F_G1 = FILTER.synthesis2();
	    Size_H0 = FILTER.size_analysis();
	    Size_H1 = FILTER.size_synthesis();
	    Size_G0 = FILTER.size_analysis2();
	    Size_G1 = FILTER.size_synthesis2();
     	H0 = new float[Size_H0];
      	H1 = new float[Size_H1];
       	G0 = new float[Size_G0];
        G1 = new float[Size_G1];

	    // If the wanted normalization is not this of the filter
     	// we need to correct the coefficient
      	float NormCorrection = 1;
	    if (FILTER.norm () == TypeNorm) NormCorrection = 1.;
	    else if (TypeNorm == NORM_L1)  NormCorrection = 1. / sqrt(2.);
     	else if (TypeNorm == NORM_L2)  NormCorrection = sqrt(2.);
	
	    for (int i=0; i < Size_H0; i++) H0[i] = F_H0[i]*NormCorrection;
     	for (int i=0; i < Size_H1; i++) H1[i] = F_H1[i]*NormCorrection;
      	for (int i=0; i < Size_G0; i++) G0[i] = F_G0[i]*NormCorrection;
       	for (int i=0; i < Size_G1; i++) G1[i] = F_G1[i]*NormCorrection;
	
        Start_H0 = - Size_H0 / 2;
        Start_H1 = - (Size_H1-1) / 2;
        Start_G0 = - Size_G0 / 2;
        Start_G1 = - (Size_G1-1) / 2; 
    }
    else
    {
	    int i;
     	TypeNorm = Norm;
      	NormCoef = (TypeNorm == NORM_L1) ? 2 : 1;
	
 	    TypeFilter = FILTER.type_filter();
      	float *F_H0 = FILTER.analysis();
       	float *F_H1 = FILTER.synthesis();
        Size_H0 = FILTER.size_analysis();
        Size_H1 = FILTER.size_synthesis();
        Size_G1 = Size_H0;
        Size_G0 = Size_H1;
	    H0 = new float[Size_H0];
	    H1 = new float[Size_H1];
	
	    // If the wanted normalization is not this of the filter
	    // we need to correct the coefficient
	    float NormCorrection=1;
	    if (FILTER.norm () == TypeNorm) NormCorrection = 1.;
	    else if (TypeNorm == NORM_L1)  NormCorrection = 1. / sqrt(2.);
	    else if (TypeNorm == NORM_L2)  NormCorrection = sqrt(2.);
	
	    for (i=0; i < Size_H0; i++) H0[i] = F_H0[i]*NormCorrection;
	    for (i=0; i < Size_H1; i++) H1[i] = F_H1[i]*NormCorrection;
	
	    // float SumH=0.,SumG=0.;
	    // for (i=0; i < Size_H0; i++) SumH += H0[i]*H0[i];
	    // for (i=0; i < Size_H1; i++) SumG += H1[i]*H1[i];
	    // cout << " BAND SUM = " << Size_H0 <<  " " << Size_H1 <<  " " << SumH << " " << SumG << endl;
     	// cout << " H0 " << endl;
      	// for (i=0; i < Size_H0; i++) printf("%f\n", H0[i]);
       	// cout << " H1 " << endl;
        // for (i=0; i < Size_H1; i++) printf("%f\n", H1[i]);
	
        G1 = new float[Size_G1];
        G0 = new float[Size_G0];
	
	    filter_complement(H0, G1, Size_H0, 1);
	    filter_complement(H1, G0, Size_H1, -1);    
        Start_H0 = - Size_H0 / 2;
        Start_H1 = - (Size_H1-1) / 2;
        Start_G0 = - Size_G0 / 2;
        Start_G1 = - (Size_G1-1) / 2; 
      }
#if FIL_DEBUG
    printf("Filter \n");
    for (int i = 0; i <  Size_H0; i++)
      printf("H0(%d) = %1.11f\n", i,  H0[i]);
    printf("\n");
    for (int i = 0; i <  Size_G0; i++)
      printf("G0(%d) = %1.11f\n", i, G0[i]);
    printf("\n");
    for (int i = 0; i <   Size_H1; i++)
      printf("H1(%d) = %1.11f\n", i,  H1[i]);
    printf("\n");
    for (int i = 0; i <  Size_G1; i++)
      printf("G1(%d) = %1.11f\n", i,  G1[i]);
    printf("ENDFilter \n");
#endif
};

SubBandFilter::~SubBandFilter() 
{
  if (H0 != NULL) { delete [] H0; H0 = NULL; }
  if (H1 != NULL) { delete [] H1; H1 = NULL; } 
  if (G0 != NULL) { delete [] G0; G0 = NULL; }
  if (G1 != NULL) { delete [] G1; G1 = NULL; }
  reset_param();
};

/*************************************************************/

void SubBandFilter::convol(int N, int StartF, int SizeF, float *Input, float *Filter, float *Output, int Step)
{    
    int i, p, Index;
    
    for (i = 0; i < N; i ++)
    {
        Output[i] = 0;
        for (p = 0; p < SizeF; p++)
		{
           Index = i - (StartF + SizeF - p - 1) * Step;
           Index = test_index(Index, N);
           Output[i] += Input[Index] * Filter[SizeF-p-1];
		} 
    }     
}

// convolves the data with the filter H0
void SubBandFilter::convol_h0 (int N, float *Input, float *Output)
{
    float *temp = new float[N];
    convol(N, Start_H0, Size_H0, Input, H0, temp);
    for (int i=0; i<N; i+=2)
        Output[i/2] = temp[i];
    delete [] temp;
}

// convolves the data with the filter G0
void SubBandFilter::convol_g0 (int N, float *Input, float *Output)
{
    float *temp = new float[N];
    convol(N, Start_G0, Size_G0, Input, G0, temp);
    for (int i=0; i<N; i+=2)
        Output[i/2] = temp[i];
    delete [] temp;
}

// convolves the data with the filter H1
void SubBandFilter::convol_h1 (int N, float *Input, float *Output)
{
    float *temp = new float[N];
    for (int i=0; i<N; i++)
        if (i % 2 == 0)
            temp[i] = Input[i/2];
        else 
            temp[i] = 0;
            
    convol(N, Start_H1, Size_H1, temp, H1, Output);
    delete [] temp;
}

// convolves the data with the filter G1
void SubBandFilter::convol_g1 (int N, float *Input, float *Output)
{
    float *temp = new float[N];
    for (int i=0; i<N; i++)
        if (i % 2 == 0)
            temp[i] = Input[i/2];
        else 
            temp[i] = 0;
            
    convol(N, Start_G1, Size_G1, temp, G1, Output);
    delete [] temp;
}

/****************************************************************************/

void SubBandFilter::transform (int N, float *SignalIn, float *SignalOut, float *DetailOut)
{
       convol_h0 (N, SignalIn, SignalOut);
       convol_g0 (N, SignalIn, DetailOut);
}

void SubBandFilter::recons (int N, float *SignalIn,  float *DetailIn, 
                           float *SignalOut)
{
    int i;
    float *Temp = new float [N];
    for (i=0; i < N; i++) Temp[i]=0;
    
    convol_h1 (N, SignalIn, SignalOut);
    convol_g1 (N, DetailIn, Temp);
    for (i = 0; i< N; i++) 
         SignalOut[i] = NormCoef*(SignalOut[i]+Temp[i]);
    delete [] Temp;
}

void SubBandFilter::transform (int N, float *SignalIn, float *SignalOut, 
                   float *DetailOut, int Step)
{
  convol(N, Start_H0, Size_H0, SignalIn, H0, SignalOut, Step);
  convol(N, Start_G0, Size_G0, SignalIn, G0, DetailOut, Step);
}


void SubBandFilter::recons(int N, float *SignalIn,  float *DetailIn, 
                           float *SignalOut, int Step)
{
  float *Temp = new float[N];
  convol(N, Start_H1, Size_H1, SignalIn, H1, SignalOut, Step);
  convol(N, Start_G1, Size_G1, DetailIn, G1, Temp, Step);
  for (int i = 0; i< N; i++) 
      SignalOut[i] = 0.5 * NormCoef * (SignalOut[i] + Temp[i]);
  delete [] Temp;
}
