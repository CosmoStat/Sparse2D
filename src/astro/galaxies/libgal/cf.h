/***********************************************************
**	Copyright (C) 1999 CEA
************************************************************
**
**    UNIT
**
**    Version:  1.0
**
**    Author: 	J.L. Starck
**
**    Date: 	27/08/99
**    
**    File:  	CF_Ana.h
**
************************************************************
**
**  1D,2D,3D correlation function analysis 
**  
************************************************************/


#ifndef	_CF_H_
#define	_CF_H_

// #include "IM_Math.h"
#include "Array.h"


// Pair histogram calculation between DistMin and DistMax with a given step

class CorrFunAna {
    float DistMin;   // Distance min between two-point  
    float DistMax;   // Distance max between two-point
    float SquareDistMin;   // Square of Distance min between two-point (useful for optimising time)
    float SquareDistMax;   // Square Distance max between two-point (useful for optimising time)
    float Step;      // Step value for the analysis
    int Nc;          // Number of bins = (DistMax - DistMin) / Step
    dblarray PairHisto; // Pair histogram
                       // PairHisto(0,j) = distance corresponding to the
					   //                 jth analyis point (j=0..Nc-1)
					   // PairHisto(1,j) = number of pairs in the bin
    
    int index_dist(float r);  // return the corresponding index in PairHisto
                              // to the distance r
  public:
    Bool Verbose;
    int np () { return Nc;}       // return the number of bins
    float step () { return Step;} // return the step
    float coord(int BinIndex) { return PairHisto(0,BinIndex);} 
                                  // return the distance
    float val(int BinIndex) { return PairHisto(1,BinIndex);}
                                  // return the number of bins at the distance
				  // given by PairHisto(0,BinIndex)
    // initialization either from the number of bins, either from the 
    // histogram				  			  
    CorrFunAna(float Dmin, float Dmax, int nbins);
    CorrFunAna(float Dmin, float Dmax, float StepVal);
	
    
    // reset the table PairHisto (but not the table dimension)                           
    void init() {PairHisto.init();}
    

    // find pairs
    void cf_find_pairs(ArrayPoint & Data, ArrayPoint & DataWeight, fltarray &CF_DataData);
    void cf_find_pairs(ArrayPoint & Data1, ArrayPoint & Data2, ArrayPoint & DataWeight1, ArrayPoint & DataWeight2, fltarray &CF_Data1Data2);

};


void make_histo(fltarray &dd, fltarray &rr, fltarray &dr, fltarray &Result);

void normalize_histo(ArrayPoint & DataWeight, ArrayPoint & RndWeight, fltarray &Result);


#endif


