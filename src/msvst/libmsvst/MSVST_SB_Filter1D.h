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
**    File:  SB_Filter1D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _MSVST_SB_FILTER1D_H_
#define _MSVST_SB_FILTER1D_H_

#include "Border.h"
#include "MSVST_Filter.h"

/***********************************************************************/
// calculate the size of the vector length after Scale times of decimated filtering
inline int size_resol(int Scale, int N)
{
    int Val = N;
    for (int s=0; s<Scale; s++) Val = (Val + 1) / 2;
    return Val;
}

/***********************************************************************/

// Class for sub-band decomposition with and without decimation
class SubBand1D {
    public:
       type_border Border;    // Border management type
       
       SubBand1D ()           
       {         
           Border = I_MIRROR;
       }
       
       SubBand1D (const SubBand1D &sb1d)
       {
           Border = sb1d.Border;
       }
                     
       inline int test_index(int i, int N)
       {
	      return get_index(i, N, Border);
       }           
	   
       // (int N, float *High, float *Low, float *Det)
       // sub-band decomposition with decimation
       // N = even number of pixels in the input signal High
       // High is decomposed in two sub signals Low and Det
       // Low = low resolution signal: its size is N/2
       // Det = Detail information:  its size is N/2
       virtual void transform (int , float *, float *, float *){};
      
       // (int N, float *Low, float *Det, float *High)
       // sub-band reconstruction with decimation
       // N = even number of pixels in the input signal High
       // inverse transform: the output High (size N) is reconstructed
       // from Low (size N/2) and Det, size N/2
       virtual void recons (int , float *,  float *,  float *){};
                    
       // (int N, float *High, float *Low, float *Det, int Step)
       // sub-band decomposition without decimation
       // High is decomposed in two sub signals Low and Det
       // Low = low resolution signal: its size is N
       // Det = Detail information:  its size is N
       virtual void transform (int , float *, float *, float *, int){};
       
       // (int N, float *Low, float *Det, float *High, int Step)
       // sub-band reconstruction without decimation
       // inverse transform: the output High (size N) is reconstructed
       // from Low (size N) and Det (size N)
       virtual void recons (int , float *,  float *,  float *, int){};
       
       virtual ~SubBand1D(){};
};

/***********************************************************************
* filter bank transform (one step transform)
************************************************************************/

class SubBandFilter: public SubBand1D
{
   int NormCoef;  // NormCoef = 2 for L1 normalization (multiply with the reconstruction coef.)
                  // NormCoef = 1 for L2 normalization
  
   float *H0, *G0, *H1, *G1; // filter banks
   int Size_H0, Size_H1, Size_G0, Size_G1; // filter sizes
   int Start_H0, Start_H1, Start_G0, Start_G1; // First position of the 
   type_sb_filter TypeFilter; // type of bands filter
   sb_type_norm TypeNorm; // filter normalisation
   
   Bool Verbose;

   void reset_param();  // reinitialize all the data fields to nulls or zeros
   void init(FilterAnaSynt & FILTER, sb_type_norm Norm);
   void init(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm);
   void init();
   
   public:
   SubBandFilter();
   SubBandFilter(const SubBandFilter &sbf);

   // construct with a filter type that creates an instance of FilerAnaSynt;
   // the default normalization will be required
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(type_sb_filter T_Filter);
   
   // construct with a filter type  that creates an instance of FilerAnaSynt,
   // Norm is the required normalization
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(type_sb_filter T_Filter, sb_type_norm Norm);
   
   // construct with a FilterAnaSynt;
   // the default normalization will be required
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(FilterAnaSynt &FILTER); 
   
   // construct with a FilterAnaSynt and a required normalization
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(FilterAnaSynt &FILTER, sb_type_norm Norm); 
   
   // construct with a file that creates an instance of FilterAnaSynt;
   // the default normalization will be required
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(char *FileName);
   
   // construct with a file that create an instance of FilterAnaSynt;
   // Norm is the required normalization
   // correction will be done if the required normaliztion is not the same with
   // that in the FilterAnaSynt instance
   SubBandFilter(char *FileName, sb_type_norm Norm);
   
   // construct by two given filters
   // H and G are supposed to L1 normalized;
   // correction will be done if FilNorm (the required normalization) is not L1;
   SubBandFilter(float *H, int Nh, float *G, int Ng, sb_type_norm FilNorm);
   
   // give the filter sizes 
   void getFilterSize(int &h0len, int &h1len, int &g0len, int &g1len) 
   { h0len = Size_H0; h1len = Size_H1; g0len = Size_G0; g1len = Size_G1; }
      
   // display info of fiters
   void displayFilterInfo()
   {
        cout << "Type Filter = " << TypeFilter << endl;
        cout << "H0 = " << endl;
        for (int i=0; i<Size_H0; i++) cout << H0[i] << " " << endl;
        cout << endl;
        cout << "G0 = " << endl;
        for (int i=0; i<Size_G0; i++) cout << G0[i] << " " << endl;
        cout << endl;
        cout << "H1 = " << endl;
        for (int i=0; i<Size_H1; i++) cout << H1[i] << " " << endl;
        cout << endl;
        cout << "G1 = " << endl;
        for (int i=0; i<Size_G1; i++) cout << G1[i] << " " << endl;
        cout << endl;
   }
   
   // convolve Input with Filter to give Output. All have length as N.
   void convol(int N, int StartF, int SizeF, float *Input, float *Filter, float *Output, int Step = 1);
   // convolved a signal by H0 filter with subsampling, input size N is even
   void convol_h0 (int N, float *Input, float *Output);
   // convolved a signal by G0 filter with subsampling, input size N is even
   void convol_g0 (int N, float *Input, float *Output); 
   // convolved a signal by H1 filter with oversampling, output size N is even
   void convol_h1 (int N, float *Input, float *Output);
   // convolved a signal by G1 filter with oversampling, output size N is even
   void convol_g1 (int N, float *Input, float *Output);
   
   // decimated one step transform
   // SignalIn = input signal of even size N
   // SignalOut = output smooth signal of size N/2
   // DetailOut = output wavelet coefficient of size N/2
  void transform (int N, float *SignalIn, float *SignalOut, float *DetailOut);
   
   // reconstruct a signal for the low resolution part and the 
   // associated decimated wavelet coefficient
   // N = ouput signal even size
   // SignalIn = input signal of size N/2
   // DetailIn = input wavelet coefficient of size N/2
   // SignalOut = outpout signal of size N
   void recons (int N, float *SignalIn,  float *DetailIn,  float *SignalOut);
   
   // one step undecimated subband transform
   // SignalIn = input signal of size N
   // SignalOut= output smooth signal of size N
   // DetailOut= output wavelet coefficient signal of size N
   // Step = distance between two cofficients
   void transform (int N, float *SignalIn, float *SignalOut, float *DetailOut, int Step);
   
   // reconstruct a signal for the low undecimated resolution part and the 
   // associated undecimated wavelet coefficient
   // N = input-ouput signal size
   // SignalIn = input signal of size N
   // DetailIn = input wavelet coefficient N
   // SignalOut = outpout signal of size N
   // Step = distance between two cofficients
   void recons (int N, float *SignalIn,  float *DetailIn,  float *SignalOut, int Step);
   
   ~SubBandFilter();
};

#endif
