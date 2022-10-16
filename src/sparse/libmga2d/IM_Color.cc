/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  24/04/01
**    
**    File:  IM_Color.cc
**
**
*******************************************************************************
**
**    DESCRIPTION  Color conversion system
**    ----------- 
**                 
******************************************************************************/

#include "GlobalInc.h"
#include "IM_Color.h"

/***************************************************************************/
 	           
void XYZ_to_LAB (float X, float Y, float  Z, float & L, float &A, float & B)
{
   double  x0, y0, z0;
   double  x, y, z;
       
   x0 = 0.98072;
   y0 = 1.00000;
   z0 = 1.18225;
       
   if ((X / x0) > 0.008856)  x = pow ((X / x0), (1.0 / 3.0));
   else                      x = (7.787 * (X / x0)) + (16.0 / 116.0);
   if ((Y / y0) > 0.008856)  y = pow ((Y / y0), (1.0 / 3.0));
   else                      y = (7.787 * (Y / y0)) + (16.0 / 116.0);
   if ((Z / z0) > 0.008856)  z = pow ((Z / z0), (1.0 / 3.0));
   else                      z = (7.787 * (Z / z0)) + (16.0 / 116.0);
       
   L = (116.0 * y) - 16.0;
   A = 504.3 * (x - y);
   B = 201.7 * (y - z);
}
     
/***************************************************************************/
    
     
void XYZ_to_LUV (float  X, float  Y, float  Z,
                 float  & L, float  & U, float  & V)
{
    double  x0, y0, z0;
    double  u0, v0;
    double  y, u, v;
       
    x0 = 0.98072;
    y0 = 1.00000;
    z0 = 1.18225;
       
    if ((Y / y0) > 0.008856)  y = pow ((Y / y0), (1.0 / 3.0));
    else                      y = (7.787 * (Y / y0)) + (16.0 / 116.0);
       
    L = (116.0 * y) - 16.0;
       
    if ((X + (15.0 * Y) + (3.0 * Z)) == 0.0) 
    {
          U = 0;
          V = 0;
     } 
     else 
     {
         u0 =  0.2009;
         v0 =  0.4610;
         u  = (4.0 * X)  / (X  + (15.0 * Y)  + (3.0 * Z));
         v  = (9.0 * Y)  / (X  + (15.0 * Y)  + (3.0 * Z));
         U = 13.0 * L * (u - u0);
         V = 13.0 * L * (v - v0);
     }
}

/***************************************************************************/

void LUV_to_XYZ (float  L, float  U, float  V, 
                 float  & X, float  & Y, float &  Z)
{
    double  x0, y0, z0;
    double  u0, v0;
    double  u, v;
    x0 = 0.98072;
    y0 = 1.00000;
    z0 = 1.18225;
    u0 =  0.2009;
    v0 =  0.4610;
	 
    u = U / (13.*L) + u0;
    v = V / (13.*L) + v0;
    Y = y0*(pow(((L+16.) /116.), 3.));
    X = 9*u/(4*v)*Y;
    Z = Y*(12. - 3*u - 20. * v) / (4.*v);
}

/***************************************************************************/

void RGB_to_XYZ (float  R, float G, float  B,
                    float & X, float & Y, float  &Z)
{
   X = (0.412453 * R) + (0.357580 * G) + (0.180423 * B);
   Y = (0.212671 * R) + (0.715160  * G) + (0.072169  * B);
   Z = ( 0.019334 * R) + ( 0.119193 * G) + (0.950227 * B);
}

/***************************************************************************/

void XYZ_to_RGB(float  X, float Y, float  Z,
                    float & R, float & G, float  &B)
{
   R = (3.240479 * X) + (-1.537150 * Y) + (-0.498535 * Z);
   G = (-0.969256 * X) + (1.875992  * Y) + (0.041556 * Z);
   B = (0.055648 * X) + (-0.204043 * Y) + (1.057311 * Z);
}

/***************************************************************************/

void RGB_to_LUV (float  R, float  G, float  B,
                    float  & L, float  & U, float  & V)
{
   float X,  Y, Z;
   RGB_to_XYZ(R,G,B,X,Y,Z);
   XYZ_to_LUV(X,Y,Z,L,U,V);
}

/***************************************************************************/
		    
void LUV_to_RGB (float L, float  U, float  V,
                    float  & R, float  & G, float  & B)
{
   float X, Y, Z;
   LUV_to_XYZ(L,U,V,X,Y,Z);
   XYZ_to_RGB(X,Y,Z,R,G,B);
 }	

/***************************************************************************/

void rgb_to_luv(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float L,  U, V;
       RGB_to_LUV(Data(i,j,0), Data(i,j,1),Data(i,j,2), L,  U, V);
       Data(i,j,0) = L;
       Data(i,j,1) = U;
       Data(i,j,2) = V;
   }
}

/***************************************************************************/

void luv_to_rgb(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float R,G,B;
       LUV_to_RGB(Data(i,j,0), Data(i,j,1),Data(i,j,2), R,G,B);
       Data(i,j,0) = R;
       Data(i,j,1) = G;
       Data(i,j,2) = B;
   }
}

/***************************************************************************/

void RGB_to_HSV(float r, float g, float b, float & h, float & s, float & v )
{
        float min, max, delta;
        fltarray Tab(3);
	Tab(0) = r;
	Tab(1) = g;
	Tab(2) = b;
        min = Tab.min();
        max = Tab.max();
        v = max;                               // v

        delta = max - min;

        if( max != 0 )
                s = delta / max;               // s
        else {
                // r = g = b = 0                // s = 0, v is undefined
                s = 0;
                h = -1;
                return;
        }

        if( r == max )
                h = ( g - b ) / delta;         // between yellow & magenta
        else if( g == max )
                h = 2 + ( b - r ) / delta;     // between cyan & yellow
        else
                h = 4 + ( r - g ) / delta;     // between magenta & cyan

        h *= 60;                               // degrees
        if( h < 0 )
                h += 360;

}

/***************************************************************************/

void rgb_to_hsv(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float H,  S, V;
       RGB_to_HSV(Data(i,j,0), Data(i,j,1),Data(i,j,2), H,  S, V);
       Data(i,j,0) = H;
       Data(i,j,1) = S;
       Data(i,j,2) = V;
   }
}

/***************************************************************************/

void HSV_to_RGB( float h, float s, float v, float & r, float & g, float & b)
{
        int i;
        float f, p, q, t;

        if( s == 0 ) {
                // achromatic (grey)
                r = g = b = v;
                return;
        }

        h /= 60;                        // sector 0 to 5
        i = (int) floor( h );
        f = h - i;                      // factorial part of h
        p = v * ( 1 - s );
        q = v * ( 1 - s * f );
        t = v * ( 1 - s * ( 1 - f ) );

        switch( i ) {
                case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                default:                // case 5:
                        r = v;
                        g = p;
                        b = q;
                        break;
        }

}

/***************************************************************************/

void hsv_to_rgb(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float R,G,B;
       HSV_to_RGB(Data(i,j,0), Data(i,j,1),Data(i,j,2), R,G,B);
       Data(i,j,0) = R;
       Data(i,j,1) = G;
       Data(i,j,2) = B;
   }
}

/***************************************************************************/

void col_rescale(fltarray &Data, float MinVal, float MaxVal)
{
   int b,i,j;
   int Nx = Data.nx();
   int Ny = Data.ny();
   int Nz = Data.nz();
   float MinT = Data.min();
   float MaxT =  Data.max();
   float ScaleT = (MaxVal-MinVal) / (MaxT - MinT);
   for (b = 0; b < Nz; b++) 
   for (i=0; i < Ny; i++)
   for (j=0; j < Nx; j++) Data(j,i,b) = (Data(j,i ,b) - MinT)*ScaleT + MinVal;
}

/***************************************************************************/
