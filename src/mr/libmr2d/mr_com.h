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
**    Date:  00/07/13
**    
**    File:  mc2d_com.h
**
*******************************************************************************/
#ifndef _MR_COM_H
#define _MR_COM_H


#include "IM_Obj.h"
#include "MR_Obj.h"

class to_Param {
public:
   //char     tc_NameIn[256];          /* input name */
   //char     tc_NameOut[256];         /* output name */
   //Ifloat   o_Imag;                  /* data in */   
   //int      i_NbLin,                 /* data in dimension */
   //         i_NbCol, 
//	    i_NbIter;  
public:
   virtual void operator () (int argc, char *argv[])=0; /* read parameters */
   virtual ~to_Param() {};
private:
   virtual void usage (char *argv[])=0;      /* trace usage */
   virtual void trace ()=0;                  /* trace parameters */
};


class to_Result {
public:
   //Ifloat      o_DataOut;          /* data out, 3D fltaaray */
   //to_Param*   po_Param;         /* associated param */
   //MultiResol* po_Mr2d;             /* multiresol of data IN */
public:
   virtual void res_FirstCompute ()=0;       /* first loop work */  
   virtual void res_CurrentCompute ()=0;     /* curent loop work */
   virtual void res_Recons ()=0;             /* inv reconstruct process */
   virtual void res_Write ()=0;              /* write 3D out result */
   virtual to_Param* res_GetpParam ()=0;     /* get to_Param pointeur */
   virtual ~to_Result() {};
};


class to_Iter {
public:
   //Ifloat     o_Residu;
   //Ifloat     o_Init;
   //int        NbIter;
public:
   //to_Iter () {NbIter=1;};
   virtual void iter_Begin ();
   virtual Bool iter_End ();
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~to_Iter() {};
private:
   virtual void iter_Residu ();  
   virtual to_Result* iter_getpResult()=0;
};

class to_SoftIter {
public:
   //fltarray   o_2dCurrentSol;   
   //int        NbIter;
public:
   //to_SoftIter () {NbIter=1;};
   virtual void iter_Begin ();
   virtual Bool iter_End ();
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~to_SoftIter() {};
private:
   virtual to_Result* iter_getpResult()=0;
   virtual void iter_ActionOnSol()=0;   
};
// verify that all images have same size 
  
#endif
