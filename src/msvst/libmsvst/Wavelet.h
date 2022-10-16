/*
 * Filename : Wavelet.h
 * 
 * Class Description
 * 1. WaveletShrinkage : Wavelet hypothesis tests
 */
 
#ifndef _WAVELET1_H
#define _WAVELET1_H

double spmpar(int *i);
double exparg(int *l);
double psi(double *xx);

#include <string>
#include <math.h>
#include "cdflib.h"
#include "Array.h"
#include "DefMath.h"
#include "MSVST_SB_Filter1D.h"
#ifdef MWIR
#include "stocc.h"
#include "cdflib.h"
#include "ImLib_mwir.h"
#include "MR_HaarPoisson.h"
#include "MWIRException.h"
#include "Wavelet_wmir.h"
#include "Fisz.h"
#else
#include "ImLib_msvst.h"
#include "PoisMSVSTException.h"
#include "ImLib_msvst.h"
#include "Wavelet_msvst.h"
#endif


 #endif
