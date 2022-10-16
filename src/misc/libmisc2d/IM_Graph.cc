


#include "IM_Obj.h"
#include "IM_Graphics.h"


C_Trigo Trigo;

/***************************** im_draw_line *****************************/

/* Draw a line in an image */

void im_draw_line(Ifloat & Ima, double Begin_X, double Begin_Y,
                  double End_X, double End_Y, float Val)
{
   double	dx,dy, slope;
   int		ix1,iy1, ix2,iy2, ix,iy;
   int Nl = Ima.nl();
   int Nc = Ima.nc();

  dx = End_X-Begin_X;
  dy = End_Y-Begin_Y;
  if (fabs(dx) > fabs(dy))
  {
    slope = dy/dx;
    ix1 = RINT(Begin_X);
    ix2 = RINT(End_X);
    if (ix2>ix1)
    {
      for (ix=ix1+1; ix<=ix2; ix++)
        if (ix >= 0 && ix < Nc)
        {
          iy = RINT(Begin_Y+(ix-Begin_X)*slope);
          if (iy >= 0 && iy < Nl) Ima(iy,ix) = Val;
        }
    }
    else
    {
      for (ix=ix1-1; ix>=ix2; ix--)
        if (ix >= 0 && ix< Nc)
          {
          iy = RINT(Begin_Y+(ix-Begin_X)*slope);
          if (iy >=0 && iy< Nl) Ima(iy,ix) = Val;
          }
    }
  }
  else
  {
    slope = dx/dy;
    iy1 = RINT(Begin_Y);
    iy2 = RINT(End_Y);
    if (iy2>iy1)
    {
      for (iy=iy1+1; iy<=iy2; iy++)
        if (iy>=0 && iy< Nl)
        {
          ix = RINT(Begin_X+(iy-Begin_Y)*slope);
          if (ix >=0 && ix< Nc) Ima(iy,ix) = Val;
        }
    }
    else
      for (iy=iy1-1; iy>=iy2; iy--)
      {
        if (iy>=0 && iy< Nl)
        {
          ix = RINT(Begin_X+(iy-Begin_Y)*slope);
          if (ix>=0 && ix< Nc) Ima(iy,ix) = Val;
        }
      }
    }
}


/******************************** circle *********************************/

/* Draw a circle in an image */
void im_draw_circle(Ifloat &Ima, double x, double y, double r, float Val)
{
   int i;
   double Begin_X, End_X, Begin_Y, End_Y;

   Begin_X = x + r;
   Begin_Y = y;
   for (i=0; i<SIZE_TAB_TRIG; i++) 
   {
       End_X = x+r*Trigo.ctg[i];
       End_Y = y+r*Trigo.stg[i];
       im_draw_line(Ima, Begin_X, Begin_Y, End_X, End_Y, Val);
       Begin_X = End_X;
       Begin_Y = End_Y;
   }
}


/******************************** ellips *********************************/

/* Draw an ellips in an image */

void im_draw_ellips(Ifloat &Ima, double x, double y, double a,
		    double b, double theta, float Val, Bool dotflag)
{
   int		i;
   double ct, st;
   double Begin_X, End_X, Begin_Y, End_Y;

   ct = cos(PI*theta/180);
   st = sin(PI*theta/180);

  Begin_X = x+a*ct;
  Begin_Y = y+a*st;

  for (i=1; i<SIZE_TAB_TRIG; i++)
  {
     if ( (dotflag == True) && !(i&1))
     {
         Begin_X = x + a*Trigo.ctg[i]*ct - b*Trigo.stg[i]*st;
         Begin_Y = y + a*Trigo.ctg[i]*st + b*Trigo.stg[i]*ct;
     }
     else 
     {
         End_X = x + a*Trigo.ctg[i]*ct - b*Trigo.stg[i]*st;
         End_Y = y + a*Trigo.ctg[i]*st + b*Trigo.stg[i]*ct;
         im_draw_line(Ima, Begin_X, Begin_Y, End_X, End_Y, Val);
         Begin_X = End_X;
         Begin_Y = End_Y;
     }
  }
}

/***************************************************************************/

