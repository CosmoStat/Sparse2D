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
**    File:  mc2d_com.h
**
*******************************************************************************/
#ifndef _MC2D_COM_H
#define _MC2D_COM_H


#include "IM_Obj.h"
#include "MR_Obj.h"

class to_Param {
public:
   char     tc_NameIn[256];          /* input cube name */
   char     tc_NameOut[256];         /* output cube name */
   fltarray ao_3dDataIn;             /* data in, 3D fltarray */
   Ifloat*  po_TabIma;               /* data in, 1D Ifloat */
   int      i_NbLin,                 /* data in dimension */
            i_NbCol, 
	    i_NbImage,
	    i_NbIter;  
public:
   virtual void operator () (int argc, char *argv[])=0; /* read parameters */
   virtual ~to_Param() {};
private:
   virtual void read_image (char *argv[])=0; /* read images in */
   virtual void usage (char *argv[])=0;      /* trace usage */
   virtual void trace ()=0;                  /* trace parameters */
};


class to_Result {
public:
   fltarray    o_3dDataOut;          /* data out, 3D fltaaray */
   //to_Param*   po_Param;             /* associated param */
   MultiResol* po_TabMr2d;           /* multiresol of data IN */
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
   fltarray   o_3dResidu;
   fltarray   o_3dInit;
   int        NbIter;
public:
   to_Iter () {NbIter=1;};
   virtual void iter_Begin ();
   virtual Bool iter_End ();
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~to_Iter() {};
private:
   virtual void iter_Residu ();  
   virtual to_Result* iter_getpResult()=0;
};


Bool control_mr2d (int NbMultRes2d, MultiResol* TabMr2d);
Bool control_image (int NbImage, Ifloat* TabImage);
// verify that all images have same size 
  
#endif
