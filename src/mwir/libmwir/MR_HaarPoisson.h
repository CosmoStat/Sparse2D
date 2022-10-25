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
**    Date:  15.7.99
**    
**    File:  MR_HaarPoisson.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for Poisson noise with Haar transform  
**    ----------- 
**                 
******************************************************************************/

#ifndef _CHAARPOISSON_H_
#define _CHAARPOISSON_H_

/****************************************************************************/

//epsilon  10^-6  10^-5  10^-4  10^-3  10^-2
const int NBR_EPS = 5;
const int NBR_LAMBA = 62;
extern double TabHaarEps[NBR_EPS];  // = {1e-6,1e-5,1e-4,1e-3,1e-2};

const double HaarMaxExp = pow((double) 2., (double) 30);
const double HaarMinExp = pow((double) 2., (double) -30);
const double HaarInfExp = pow((double) 2., (double) -19);

double get_harr_poisson_threshold(double Lambda, double Eps);
// return the threshold level using Jammal-Bijaoui table
// Lambda = number of count per pixel
// Eps = probability of false detection
#endif

  

