/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12 
**    
**    File:  MR_Edge.h
**
*******************************************************************************
**
**    DESCRIPTION  Contour detection in the multiresolution space
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

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "MR_Edge.h"

const int NBR_TCONT = 8;

/***************************************/
 
// return the coordinated of the two points related to a pixel i,j
// in a contour with a angle given in Angle
// Tcont = in: Contour type
// i,j = in: pixel position
// i1,j1 = out: position of the first pixel in the contour
// i2,j2 = out: position of the second pixel in the contour
// DistPix = in: distance between two consecutive pixel
// return False if no contour is found in i,j

Bool pix_tcont(int Tcont, int i, int j, 
                     int & i1, int & j1, int & i2, int &j2, int DistPix)
{
   int Di1, Di2, Dj1, Dj2;
   Bool ValRet = True;
   Di1 =  Dj1 =  Di2 = Dj2 = 0;
   switch(Tcont)
   {
            case 1: Dj1 = -1; Dj2 = 1; break;  // Angle = 90
            case 2: Di1 = -1; Di2 = 1; break;  // Angle = 0
            case 3: Di1 = -1; Dj1 = -1; Di2= 1; Dj2=  1; break; // Angle = 45
            case 4: Di1 = -1; Dj1 =  1; Di2= 1; Dj2= -1; break; // Angle = 135
            case 5: Di1 = -1; Dj2 = -1;break; // cont =     _
	                                      //             |
            case 6: Di1 = -1; Dj2 = 1; break; // cont =       _
	                                      //             |
            case 7: Dj1 = -1; Di2 = 1; break; // cont =      _|
	                                      //             
            case 8: Dj1 =  1; Di2 = 1; break; // cont =      |_
	                                      //             
            default: ValRet = False; break; // pb
   }  
   i1 = i + Di1 * DistPix;
   i2 = i + Di2 * DistPix;
   j1 = j + Dj1 * DistPix;
   j2 = j + Dj2 * DistPix; 
   return ValRet;
}


/***************************************/
 
// return the minumum value between the 3 in the contour
float val_contour_min(Ifloat &Data, Ifloat &ImaEdge, Iint & Angle,
                             int i, int j, float Level, int Step)
{ 
   int i1,j1,i2,j2;
   float ValRet;
   int K = POW2(Step);
   
   if ((ABS(ImaEdge(i,j)) > Level) 
        && (pix_tcont(Angle(i,j), i, j, i1, j1, i2, j2, K) == True))
   {
      ValRet = Data(i,j);
      if (ABS(ValRet) > ABS(Data(i1,j1,I_MIRROR)))  ValRet = Data(i1,j1,I_MIRROR);
      if (ABS(ValRet) > ABS(Data(i2,j2,I_MIRROR)))  ValRet = Data(i2,j2,I_MIRROR);
      if (ValRet*Data(i,j) < 0) ValRet = 0.;
   }
   else ValRet = 0;
   return ValRet;
}
 
/***************************************/

void im_get_edge(Ifloat &Ima, Iint &Angle, Ifloat &ImaEdge, int Step)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int c,i,j,i1,j1,i2,j2;
   fltarray Tab(NBR_TCONT);
   int Ind_Angle, Ind_Angle_Min;
   type_border Border=I_MIRROR;
   int K = POW2(Step);
   float Min;
    

   for (i=0;i<Nl;i++)
   for (j=0;j<Nc;j++)
   {
      for (c = 1; c <= NBR_TCONT; c++)
      {
        if (pix_tcont(c, i,j,i1,j1,i2,j2,K) == True)
 	   Tab(c-1) = (Ima(i1,j1,Border) + Ima(i,j) + Ima(i2,j2,Border)) / 3.;
	else Tab(c-1) = 0.;
      }  

      //ImaEdge(i,j) = max(Tab, Ind_Angle);
      ImaEdge(i,j) = Tab.max(Ind_Angle);
      //Min = min (Tab, Ind_Angle_Min);
      Min = Tab.min (Ind_Angle_Min);
      if (ABS(Min) > ABS(ImaEdge(i,j)))
      {
         ImaEdge(i,j) = Min;
         Angle(i,j) = Ind_Angle_Min;
      }
      Angle(i,j) = Ind_Angle + 1;
   }
}

/*********************************************************************/

// calculate the edge at each scale of the wavelet transform
//  MR_Edge: output mr edge,  normalized by the noise 
 
/* Algorithm:
      for a band b in MR_Data
         call im_get_edge which do the following points:
              . average 3 values in the 8 directions for each pixel of the band  
	      . retain the max value and the related angle
	 search the edge by zero crossing in the band b 
*/	 
void mr_band_edge(MultiResol &MR_Data, Ifloat &Edge, Iint &TCont, 
                  int b,  MRNoiseModel & ModelData)
{
   float NoiseEdge;
   int i,j;
   int Nl = MR_Data.size_band_nl(b);
   int Nc = MR_Data.size_band_nc(b);  
      
   Edge.resize(Nl,Nc);
   TCont.resize(Nl,Nc);
   Ifloat ImaAux(Nl, Nc, "ImaAux");
     
   int Step = (MR_Data.band_to_scale(b) == b) ? b : 0;
       
    im_get_edge(MR_Data.band(b), TCont, ImaAux , Step);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        NoiseEdge = ModelData.sigma(b,i,j) /sqrt(3.);
        Edge(i,j) = ABS(ImaAux(i,j)) / NoiseEdge;
                 
        // test zero crossing
        if (( MR_Data(b,i,j)* MR_Data(b,i+1,j,I_MIRROR) > 0) &&
                  ( MR_Data(b,i,j)*  MR_Data(b,i+1,j+1,I_MIRROR) > 0) &&
                  ( MR_Data(b,i,j)*  MR_Data(b,i,j+1,I_MIRROR) > 0))
                                                             Edge(i,j) = 0.;
    }
            
    for (i=0;i < Nl; i++)
    for (j=0;j < Nc; j++) 
    {
        // supress isolated pixel
        if (isolated_pixel(Edge,i,j,I_MIRROR) == True) Edge(i,j) = 0.;       
    } 
}  
 



/*********************************************************************/

// calculate the edge at each scale of the wavelet transform
//  MR_Edge: output mr edge,  normalized by the noise 
//  ImaEdge: output image edge
// ImaAngle: output image abgle

/* Algorithm:
      for each band b in MR_Data
         call im_get_edge which do the following points:
              . average 3 values in the 8 directions for each pixel of the band  
	      . retain the max value and the related angle
	 search the edge by zero crossing in the band b 
*/
	 
void mr_get_edge(MultiResol &MR_Data, MultiResol &MR_Edge, 
                 Ifloat &ImaEdge, Iint &TCont, MRNoiseModel & ModelData)
{
    int i,j,b;
    int Nl = MR_Data.size_ima_nl();
    int Nc = MR_Data.size_ima_nc();
    int Nbr_Plan = MR_Data.nbr_scale();
    type_transform Transform = MR_Data.Type_Transform;
    
    if (SetTransform(Transform) != TRANSF_PAVE)
    {
       cerr << "Error: Edge detection routine works only with PAVE transform ... " << endl;
       exit(-1);
    }
    
    MR_Edge.alloc(Nl, Nc,Nbr_Plan, Transform, "MRNoiseModel");
    ImaEdge.alloc(Nl,Nc, "ImaEdge");
    TCont.alloc(Nl,Nc, "ImaAngle");
        
    Iint Angle(Nl,Nc, "Angle");
    for (b = 0; b < MR_Edge.nbr_band()-1; b++)
    {
       mr_band_edge(MR_Data, MR_Edge.band(b), Angle, b, ModelData);
 
       for (i=0;i < Nl; i++)
       for (j=0;j < Nc; j++) 
       {
          // the cumulated edge map ImaEdge is updated
          if ((MR_Edge(b,i,j) > FLOAT_EPSILON) && (ImaEdge(i,j) <  MR_Edge(b,i,j)))
          {
                 ImaEdge(i,j) =  MR_Edge(b,i,j);
                 TCont(i,j) = Angle(i,j);
          }            
       } 
    } 
}  
 
/*********************************************************************/

void mr_zero_cross_edge(MultiResol &MR_Data, Bool KillIsol)
// 
{
    int i,j,b;
    int Nl = MR_Data.size_ima_nl();
    int Nc = MR_Data.size_ima_nc();
    int Nbr_Plan = MR_Data.nbr_scale();
    type_transform Transform = MR_Data.Type_Transform;
    Ifloat ZC(Nl,Nc,"zero cross");
    
    if (SetTransform(Transform) != TRANSF_PAVE)
    {
       cerr << "Error: Edge detection routine works only with PAVE transform ... " << endl;
       exit(-1);
    }
    MR_Data.band(Nbr_Plan-1).init();
    for (b = 0; b < MR_Data.nbr_band()-1; b++)
    {
       ZC = MR_Data.band(b);
       for (i=0;i < Nl; i++)
       for (j=0;j < Nc; j++) 
       {
          // test zero crossing
          if (( MR_Data(b,i,j)* MR_Data(b,i+1,j,I_MIRROR) > 0) &&
                  ( MR_Data(b,i,j)*  MR_Data(b,i+1,j+1,I_MIRROR) > 0) &&
                  ( MR_Data(b,i,j)*  MR_Data(b,i,j+1,I_MIRROR) > 0))
                                                            ZC(i,j) = 0.; 
       }
       
       MR_Data.band(b) = ZC;
       // supress isolated pixel
       if (KillIsol == True)
       {
          for (i=0;i < Nl; i++)
          for (j=0;j < Nc; j++) 
             if (isolated_pixel(ZC,i,j,I_MIRROR) == True)  MR_Data(b,i,j) = 0.;       
       } 
     }
}  

/*********************************************************************/
