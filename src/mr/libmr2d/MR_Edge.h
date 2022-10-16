/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  24/07/98 
**    
**    File:  MR_Edge.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#ifndef __MREDGE__
#define __MREDGE__


/*********************************************************************/ 

void mr_zero_cross_edge(MultiResol &MR_Data, Bool KillIsol = True);
/*
   supress at each scale all pixels which are not a zero crossing
   if KillIsol == True, isosalted pixels are also suppressed
*/

/*********************************************************************/ 

void mr_band_edge(MultiResol &MR_Data, Ifloat &Edge, Iint &TCont, 
                  int b,  MRNoiseModel & ModelData);
/* Calculate the edge in a band of the multiresolution transform. 
    MR_Data = in: wavelet transform of the data
    ModelData = in: noise model class
    Edge = out: detected edge at scale b		  
    TCont= out: type of contour	at scale b
*/
    
/*********************************************************************/ 
		  
void mr_get_edge(MultiResol &MR_Data, MultiResol &MR_Edge, 
                 Ifloat &ImaEdge, Iint &TCont, MRNoiseModel & ModelData);
/*
 calculate the edge of an image in multiresolution space, and
 derive an estimation of the edge map.
 The type of transform must be a PAVE transform.
 
 MR_Data = in: wavelet transform of the data
 ModelData = in: noise model class
 MR_Edge = out: detected edge per scale
 ImaEdge = out: normalized edge map
 TCont= out: type of contour
*/

/*********************************************************************/ 

Bool pix_tcont(int Tcont, int i, int j, 
                     int & i1, int & j1, int & i2, int &j2, int DistPix);
/*
  return the coordinated of the two points related to a pixel i,j
  in a contour with a angle given in Angle
  Tcont = in: Contour type
  i,j = in: pixel position
  i1,j1 = out: position of the first pixel in the contour
  i2,j2 = out: position of the second pixel in the contour
  DistPix = in: distance between two consecutive pixel
  return False if no contour is found in i,j
        1  Angle = 90
        2  Angle = 0
        3  Angle = 45
        4  Angle = 135
        5  cont =     _
	             |
        6  cont =     _
	               |
		       
        7  cont =     _|
	                                                   
        8  cont =    |_
*/

/*********************************************************************/ 

float val_contour_min(Ifloat &Data, Ifloat &ImaEdge, Iint & Angle,
                             int i, int j, float Level, int Step);

/*  return the minimum value between the 3 in the contour at a pixel i,j
    Data = in: input image
    ImaEdge = in: edge map of Data
    Angle= in: type of contour
    int i,j = in: pixel position
    float Level= in: detection level for a contour
    Step = in: 2^Step = distance between two consecutive pixels
    
*/

/*********************************************************************/ 

void im_get_edge(Ifloat &Ima, Iint &Angle, Ifloat &ImaEdge, int Step);

/* calculate the edge of a derivative image
   Ima = in: derivative image or wavelet band
   Angle = out: type of contour
   ImaEdge = out: edge map
   Step = in: 2^Step = distance between two consecutive pixels
*/


#endif




