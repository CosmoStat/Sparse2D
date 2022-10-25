/******************************************************************************
**                   Copyright (C) 1996 OCA-CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 4.1
**
**    Author: Jean-Luc Starck
**
**    Date:  92/10/09
**    
**    File:  IM_Isoc.c
**
*******************************************************************************
**
**    DESCRIPTION : creates an isophot image from an image
**    ----------    
**
**    PARAMETRES  
**    ----------  
**
**    RESULTS
**    -------
**
**    LIMITS
**    ------
**
**    CONSTRAINTS
**    -----------
**
**    BUGS
**    ----
**
******************************************************************************/

// static char sccsid[] = "@(#)Lib_Isoc_Ima.c 4.1 96/08/09 OCA-CEA @(#)";


#include "IM_Obj.h"
#include "IM_Graphics.h"

/***************************************************************************/

struct point
 {float x;
  float y;
  float I;
 };

/*************************************************************************/

static void ordonner(struct point p1, struct point p2, struct point p3,
                struct point *tri1, struct point *tri2, struct point *tri3)
{
if (p3.I>p2.I) {
               if (p3.I>p1.I) {
                              if (p2.I>p1.I) {
					     *tri1 =p1;
					     *tri2 =p2;
					     *tri3 =p3;
					     }
			      else {
				     *tri1 = p2;
				     *tri2 = p1;
				     *tri3 = p3;
				    }
			       }
		else {
		     *tri1 = p2;
		     *tri2 = p3;
		     *tri3 = p1;
		    }
		}
else {
     if (p2.I>p1.I) {
                    if (p3.I>p1.I) {
			           *tri1 = p1;
				   *tri2 = p3;
				   *tri3 = p2;
				    }
                    else {
		         *tri1 = p3;
			 *tri2 = p1;
			 *tri3 = p2;
			 }
                    }
    else {
	 *tri1 = p3;
	 *tri2 = p2;
	 *tri3 = p1;
	 }
     }


}

/*************************************************************************/

static void calcul_inter(struct point t1,struct point *t2,struct point t3)
{
    float coef;

    coef = (t2->I - t1.I) / (t3.I-t1.I);
    t2->x = coef * (t3.x-t1.x) + t1.x;
    t2->y = coef * (t3.y-t1.y) + t1.y;
}

/*************************************************************************/

static void line (int x1,int y1,int x2,int y2, float Level, Ifloat &Image)
/* 
   draw a line in the image Image between points of coordinates
   x1,y1 and x2,y2. Level specifies the level intensity of the line.
*/
{
    int j,pasx,x,y,dx,dy,cumul;
    float pasy;

    if (x1 < x2) pasx = 1;
    else pasx = -1;
    if (y1 < y2) pasy = 1;
    else pasy = -1;
    dx = abs (x2 - x1);
    dy = abs (y2 - y1);
    x= x1;
    y = y1;
    Image (x,y1) = Level;
    
    if (dx > dy)
    {
        cumul = dx / 2;
        for (j = 0; j < dx; j++)
        {
            x += pasx;
            cumul += dy;
            if (cumul > dx)
            {
                cumul -= dx;
                y += (int) pasy;
            }
            Image (x,y) = Level;
        }
    }
    else
    {
        cumul = dy / 2;
        for (j = 0; j < dy; j++)
        {
            y += (int) pasy;
            cumul += dx;
            if (cumul > dy)
            {
                cumul -= dy;
                x += pasx;
            }
            Image (x,y) = Level;
        } 
    }
}       
    

/******************************************************************************/

static void range(struct point p1,struct point  p2, Ifloat &Tab_Iso, 
            int Mod_Visu, float zoom)
{
    int x1,y1,x2,y2;
    float Level;
    /*
    x1= p1.x * zoom + 0.5 ;
    y1= p1.y * zoom + 0.5 ;
    x2= p2.x * zoom + 0.5 ;
    y2= p2.y * zoom + 0.5 ;
    */

    if (Mod_Visu == 0) Level = 1.;
    else Level = p1.I;
    x1 = (int) (p1.x * zoom + 0.5);
    y1 = (int) (p1.y * zoom + 0.5);
    x2 = (int) (p2.x * zoom + 0.5);
    y2= (int) (p2.y * zoom + 0.5);
    line (x1, y1, x2, y2, Level, Tab_Iso);
}

/*************************************************************************/

static void calcul_iso(struct point tri1,struct point tri2, struct point tri3,
                 float seuil_min, float pas,
                 Ifloat &Tab_Iso, int Mod_Visu, float zoom)
/* compute the isophots between 3 points tr1,tr2,tr3 */
{
    int nbiso1,nbiso3,nbiso,k;
    struct point inter1,inter2;

    nbiso1 = (int) ((tri1.I-seuil_min)/pas);
    nbiso3 = (int) ((tri3.I-seuil_min)/pas);
    nbiso = nbiso3 - nbiso1;
    for (k=1;k<=nbiso;k++)
    {
        inter1.I = seuil_min + pas * (nbiso1 + k);
        calcul_inter(tri1,&inter1,tri3);
        inter2.I=inter1.I;
        if (inter2.I>tri2.I) calcul_inter(tri2,&inter2,tri3);
        else calcul_inter(tri1,&inter2,tri2);
        range(inter1,inter2,Tab_Iso,Mod_Visu,zoom);
    }
}


/***************************************************************/

void im_isophot (Ifloat &Pict, Ifloat &Tab_Iso, float seuil_min, 
            float Seuil_Max, float pas, int Mod_Visu)
{
    int i,j;
//    float *ptr1,*ptr2,*ptr11,*ptr22;
    struct point cinq,p1,p2,p3,p4,tri1,tri2,tri3;
    float zoom;
    int Nl = Pict.nl();
    int Nc = Pict.nc();
    int Nl_Iso = Tab_Iso.nl();
//    int Nc_Iso = Tab_Iso.nc();

    zoom = (float)(Nl_Iso) / (float)(Nl);

    Tab_Iso.init();

/*                        DEBUT DU TRAITEMENT                         */

    for (i = 0; i < Nl - 2 ; i++)
    {
        for ( j = 0; j < Nc - 2; j++)
        {
            /* calcul des valeurs des 4 points de la fenetre */
            p1.I = Pict(i,j); 
            p1.x = i;
            p1.y = j;
            p2.I = Pict(i,j+1); 
            p2.x = i;
            p2.y = j+1;
            p3.I = Pict(i+1,j); 
            p3.x = i+1;
            p3.y = j;
            p4.I = Pict(i+1,j+1); 
            p4.x = i+1;
            p4.y = j+1;

            /* calcul des valeurs du cinquieme point */
            cinq.I = (p1.I + p2.I + p3.I + p4.I) / 4.;
            if ((cinq.I >= seuil_min) && (cinq.I < Seuil_Max))
            {
                cinq.x = i + 0.5;
                cinq.y = j + 0.5;
           
                /*   calcul des isophotes pour le premier triangle   */
                ordonner(p1,p2,cinq,&tri1,&tri2,&tri3);
                calcul_iso(tri1,tri2,tri3,seuil_min,pas,Tab_Iso,Mod_Visu,zoom);

                /*   calcul des isophotes pour le second triangle   */
                ordonner(p1,p3,cinq,&tri1,&tri2,&tri3);
                calcul_iso(tri1,tri2,tri3,seuil_min,pas,Tab_Iso,Mod_Visu,zoom);

                /*   calcul des isophotes pour le troisieme triangle   */
                ordonner(p2,p4,cinq,&tri1,&tri2,&tri3);
                calcul_iso(tri1,tri2,tri3,seuil_min,pas,Tab_Iso,Mod_Visu,zoom);

                /*   calcul des isophotes pour le quatrieme triangle   */
                ordonner(p3,p4,cinq,&tri1,&tri2,&tri3);
                calcul_iso(tri1,tri2,tri3,seuil_min,pas,Tab_Iso,Mod_Visu,zoom);
            }
        }
    }
}

/***************************************************************/

void im_quick_isophot (Ifloat &Pict, Ifloat &Tab_Iso, float seuil_min, float Seuil_Max, float pas, int Mod_Visu)
{
    int i,j,k,l;
    float  p[3],Val;
    // float zoom;
    float Level;
    int Inter, Iso;
    int Nl = Pict.nl();
    int Nc = Pict.nc();

    for (i = 0; i < Nl-1; i++)
    for (j = 0; j < Nc-1; j++)
    {
        Tab_Iso(i,j) = 0.;
        p[0] = Pict(i,j);
        p[1] = Pict(i+1,j);
        p[2] = Pict(i,j+1);
        for (k = 0; k < 2; k++)
        for (l = k+1; l < 3; l++)
            if (p[k] < p[l]) 
            {
                Val = p[l];
                p[l] = p[k];
                p[k] = Val;
            }
        Inter = (int)((p[0] - seuil_min) / pas);
        Iso = (int)((p[2] - seuil_min) / pas - Inter);
        Level = seuil_min + (Iso+Inter)*pas;
        if ((Iso != 0) && (Level > seuil_min) && (Level < Seuil_Max))
        {
            if (Mod_Visu == 0) Tab_Iso(i,j) = 1.;
            else Tab_Iso(i,j) = Level;
        }
    }
}
 
/***************************************************************/
