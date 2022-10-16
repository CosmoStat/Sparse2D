/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.1
**
**    Author: 08/23/99
**
**    Date:  01/02/09
**    
**    File:  mc1d_com.h
**
*******************************************************************************/
#ifndef _MC1D_COM_H
#define _MC1D_COM_H


#include "GlobalInc.h"
#include "IM_Obj.h"
#include "MR1D_Obj.h"

enum CorrelComputeType1d {E_LOCAL1D=0,
			  E_GLOBAL1D=1,
			  E_GLOBAL_WITHOUT_LAST1D=2,
			  E_IMPORT1D=3};
enum CorrelNoiseType1d   {E_WITHOUT1D=0,
                          E_THRESHOLD1D=1,
                          E_LINEAR1D=2};
			  //E_PROBA1D=3

const char* CorCompTransform1d (CorrelComputeType1d CorrelComp);
const char* CorNoiseTransform1d (CorrelNoiseType1d CorrelNoise);			  
void correl_compute1d (CorrelComputeType1d CorrelComp);
void correl_noise1d (CorrelNoiseType1d CorrelNoise);



class to_Param1d {
public:
   char      tc_NameIn[256];          /* input cube name */
   char      tc_NameOut[256];         /* output cube name */
   fltarray  ao_2dDataIn;             /* data in, 2D fltarray */
   fltarray* po_TabSpectre;           /* data in, 1D fltarr */
   int       i_NbPts,                 /* data in dimension */
	     i_NbSpectre,
	     i_NbIter;  
public:
   virtual void operator () (int argc, char *argv[])=0; /* read parameters */
   virtual ~to_Param1d() {};
private:
   virtual void read_spectre (char *argv[])=0; /* read images in */
   virtual void usage (char *argv[])=0;      /* trace usage */
   virtual void trace ()=0;                  /* trace parameters */
};


class to_Result1d {
public:
   fltarray    o_2dDataOut;          /* data out, 3D fltaaray */
   //to_Param*   po_Param;             /* associated param */
   MR_1D* po_TabMr1d;                /* multiresol of data IN */
public:
   virtual void res_FirstCompute ()=0;       /* first loop work */  
   virtual void res_CurrentCompute ()=0;     /* curent loop work */
   virtual void res_Recons ()=0;             /* inv reconstruct process */
   virtual void res_Write ()=0;              /* write 2D out result */
   virtual to_Param1d* res_GetpParam ()=0;   /* get to_Param pointeur */
   virtual ~to_Result1d() {};
};


class to_Iter1d {
public:
   fltarray   o_2dResidu;
   fltarray   o_2dInit;
   int        NbIter;
public:
   to_Iter1d () {NbIter=1;};
   virtual void iter_Begin ();
   virtual Bool iter_End ();
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~to_Iter1d() {};
private:
   virtual void iter_Residu ();  
   virtual to_Result1d* iter_getpResult()=0;
};



class to_SoftIter1d {
public:
   //fltarray   o_2dCurrentSol;   
   int        NbIter;
public:
   to_SoftIter1d () {NbIter=1;};
   virtual void iter_Begin ();
   virtual Bool iter_End ();
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~to_SoftIter1d() {};
private:
   virtual to_Result1d* iter_getpResult()=0;
   virtual void iter_ActionOnSol()=0;   
};


Bool control_mr1d    (int NbMultRes1d, MR_1D* TabMr1d);
Bool control_spectre (int NbSpectre, fltarray* TabSpectre);
// verify that all images have same size 
  
#endif
