/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: 10/30/02
**
**    Date:  02/10/30
**    
**    File:  mc_filter.h
**
*******************************************************************************/

#ifndef _MC_FILTER_H
#define _MC_FILTER_H


class mc_filter {

private:
   
public:

   mc_filter () {};
   virtual ~mc_filter () {};
   virtual mc_filter* filt_GetpFilter()=0;
   
   virtual void filt_Compute (MultiResol *TabMR, 
                              MRNoiseModel *TabMRNoise) {};
   
   virtual void filt_Transform (MultiResol  *TabMRin,  
                                MultiResol *TabMRout,
                                MRNoiseModel *TabNoiseIn=NULL, 
		                MRNoiseModel *TabNoiseOut=NULL) {
				
      if (TabNoiseIn != (MRNoiseModel*)NULL) {
         filt_TransfSignal (TabMRin, TabMRout, NULL);
         filt_ComputeNoise (TabMRin, TabNoiseIn, TabNoiseOut); 
         filt_Threshold (TabMRout, TabNoiseOut);
      }
      else filt_TransfSignal (TabMRin, TabMRout);
   };				
   
   virtual void filt_Invtransform (MultiResol *TabMRin , 
                                   MultiResol *TabMRout) {};
   
   virtual void filt_TransfSignal (MultiResol *TabMRin , 
                                   MultiResol *TabMRout,
                                   MRNoiseModel *TabNoiseIn=NULL) {};
   
   virtual void filt_Threshold  (MultiResol  *TabMRin_out, 
                                 MRNoiseModel *TabNoiseOut) {};
   
   virtual void filt_ComputeNoise  (MultiResol *TabMRin, 
                                    MRNoiseModel *TabNoiseIn, 
                                    MRNoiseModel *TabNoiseOut) {}; 
   
   virtual void filt_Print () {};
			   
};

#endif
