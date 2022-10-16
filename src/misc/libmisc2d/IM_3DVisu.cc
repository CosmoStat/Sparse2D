/******************************************************************************
**                   Copyxright (C) 1996 OCA-CEA
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
**    File:  IM3D_Visu.c
**
*******************************************************************************
**
**    DESCRIPTION : represent an image as 3 dimensional view
**    ----------    (intensity is the third dimension)
**
**    PARAMETRES  
**    ----------  
**                Angle : visualization angle 
**                        Angle = 0 .. 90 degres 
**                        90 degres corresponds an upper view
**                Increment = specify the used lines number
**                            Increment = 1 .. 10
**                            Number of displayed lines = 256 / Increment   
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

// static char sccsid[] = "@(#)Lib_Caval_Ima.c 4.1 92/10/09 OCA-CEA @(#)";


#include "IM_Obj.h"
#include "IM_Graphics.h"


/******************************************************************************/

float x,deltax; 
int xleft,xright,it,visible;
float X1,Y1;

/******************************************************************************/

static void move(float x,float y)
{
    X1 = x;
    Y1 = y;
}

/***************************************************************************/

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
    dx = ABS  (x2 - x1);
    dy = ABS (y2 - y1);
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

static void draw(float x,float y, Ifloat &Image)
{
    int x1,y1,x2,y2;
    int nx = Image.nc();
    int ny = Image.nl();

    x1 = (int)( X1 * (float) nx + 0.5);
    y1 = (int)(Y1 * (float) ny + 0.5);
    x2 = (int)(x * (float) nx + 0.5);
    y2 = (int)(y * (float) ny + 0.5);

    /* printf ("(%d,%d)->(%d,%d),  ", x1,y1,x2,y2);*/
    line (y1, x1, y2, x2, 1., Image);
    X1 = x;
    Y1 = y;
}

/******************************************************************************/

static void routine3d (float *pttab, int nx, Ifloat &Tab_Caval)
{ 
 float al,ar,em,xx,y;

 if (pttab[it-1]<0.)
  {
   if (!visible) goto END;
  }
  else
   if (visible) 
    {
     draw(x,pttab[it-1],Tab_Caval);
     goto END;
    }

 if (xleft-1<=nx-1 || !((xleft-1)%nx) || xright-1<=nx-1 || !((xright-1)%nx)) 
  {
   visible=0;
   goto END;
  }

 al=ABS(pttab[xleft-1]);
 ar=ABS(pttab[xright-1]);

 if (pttab[xleft-1] < 0.)
  {
   do
    {
     xleft-=nx+1;
    }
   while (pttab[xleft-1] < 0.);
   
   xright=xleft+1;
   if (pttab[xright-1] < 0.) xright=xleft-nx;
  }
  else
   {
    do 
     {
      xright-=nx+1;
     } 
    while (pttab[xright-1] < 0.);
    
    xleft=xright+nx;
    if (pttab[xleft-1] < 0.) xleft=xright-1;
   }

 em=ABS(pttab[xright-1])-ABS(pttab[xleft-1]);
 xx=(ABS(em-ar+al) < 0.0001) ? .5 : (al-ABS(pttab[xleft-1]))/(em-ar+al);
 y=em*xx+ABS(pttab[xleft-1]);

 if (deltax > 0.) xx=1.-xx;

 xx=x-xx*deltax;

 if (visible)
  {
   draw(xx,y,Tab_Caval);
   visible=0;    
  }
  else  
   {
    move(xx,y);
    visible=1;
    draw(x,pttab[it-1],Tab_Caval);
   }

 END:
  xx=0;
}

/******************************************************************************/

void im_3dview_1 (Ifloat &Imag, Ifloat &Tab_Caval, float angle, int Step)
{
 register int i,j,ki,kj;
 float deltay,deltav,peak,dx,ao=0.;
 float buf1,buf2;
 int incx=0;
 float **array, *tab, *tab1,*Pict;
 float fmin,fmax;
 int nx,ny;
 int Nl_Imag = Imag.nl();
 int Nc_Imag = Imag.nc();

 nx=Nc_Imag/Step; 
 ny=Nl_Imag/Step; 
 Tab_Caval.init();

 tab1 = tab = (float *)malloc (nx*ny*sizeof(float)); 
 array = (float **)malloc (ny*sizeof(float *)); 
 for (i=0;i<ny;i++)
 {
    array[i]=tab;
    tab+=nx; 
 } 

 Pict = (float *)malloc (nx*ny*sizeof(float));
 for (i=0;i<Nl_Imag;i+=Step)
 for (j=0;j<Nc_Imag;j+=Step) Pict[incx++]=Imag(i,j);
 fmax = fmin = Imag(0);
 for (i = 0; i<nx*ny; i++)
 {
    if (Pict[i] > fmax) fmax = Pict[i];
    if (Pict[i] < fmin) fmin = Pict[i];
 }

 for (i=0;i<ny;i++)
  for (j=0;j<nx;j++)
    array[i][j]=ABS(Pict[j+ny*i]);
 free ((char *) Pict);

 if  (nx<2 || ny<2) goto END;
 
 deltax=1./(nx+ny);
 deltay=deltax*sin (angle/58.);
 deltav=1.-ABS(sin (angle/58.));
 fmax=(fmax<0.0001) ? deltav : fmax-fmin;
 deltav/=fmax;

 for (j=0;j<ny;j++)
  {
   kj=ny-j-1;
   ki=0;
   peak=array[kj][ki]=deltay*(ki+kj+2) + deltav*(array[kj][ki]-fmin);
   while (++ki<=nx-1 && ++kj<=ny-1)
    {
     array[kj][ki]=deltay*(ki+kj+2) + deltav*(array[kj][ki]-fmin);
     if (array[kj][ki]>peak) peak=array[kj][ki];
      else array[kj][ki]=-ABS(array[kj][ki]);
    }
  }

 for (i=1;i<nx;i++)
  {
   ki=i;
   kj=0;
   peak=array[kj][ki]=deltay*(ki+kj+2) + deltav*(array[kj][ki]-fmin);
   while (++ki<=nx-1 && ++kj<=ny-1)
    {
     array[kj][ki]=deltay*(ki+kj+2) + deltav*(array[kj][ki]-fmin);
     if (array[kj][ki]>peak) peak=array[kj][ki];
      else array[kj][ki]=-ABS(array[kj][ki]);
    }
  }

 buf1=deltax*(nx+ny-2);
 buf2=deltay*(nx+1);
 move(buf1,buf2);

 buf1=deltax*(ny-1);
 buf2=deltay*2.;
 draw(buf1,buf2,Tab_Caval);

 buf2=deltay*(ny+1);
 draw(ao,buf2,Tab_Caval);

 for (j=0;j<ny;j+=2)
  {
   kj=ny-j-1;
   dx=deltax*(j-1);
   x=dx+deltax;
  
   buf1=deltay*(kj+2);
   move(x,buf1);
   draw(x,array[kj][0],Tab_Caval);
  
   visible=1;

   for (i=1;i<nx;i++)
    {
     it=xright=1+i+nx*kj;
     xleft=xright-1;
     x=dx+deltax*(i+1);
     routine3d (&(array[0][0]),nx, Tab_Caval);
    }

   kj--;

   if (kj<0) 
    {
     ki=nx-1;
     incx=-1;
     deltax=-deltax;
     dx=deltax*(-ki-ny);   
     x=dx+deltax;
     buf1=deltay*(ki+2);
     move(x,buf1);
     draw(x,array[0][ki],Tab_Caval);

     goto chat_huant;
    }

   visible=(array[kj][nx-1]>=0.);

   dx=deltax*(nx+j+1);
   buf1=dx-deltax;
   if (visible) move(buf1,array[kj][nx-1]);

   deltax=-deltax;

   for (i=1;i<nx;i++)
    {
     ki=nx-i;
     it=xleft=ki+nx*kj;
     xright=xleft+1;
     x=dx+deltax*(i+1);
     routine3d (&(array[0][0]),nx, Tab_Caval);
    }

   x=dx+deltax*nx;
   if (!visible) move(x,array[kj][0]);

   buf1=deltay*(kj+2);
   draw(x,buf1,Tab_Caval);

   deltax=-deltax;
  }

 incx=1;
 ki=0;

 do
  {
   dx=deltax*(ki+ny);
   deltax=-deltax;
   x=dx+deltax;
   move(x,array[0][0]);
   chat_huant :
   visible=1;

   for (j=1;j<ny;j++)
    {
     it=xleft=ki+1+nx*j;
     xright=xleft-nx;
     x=dx+deltax*(j+1);
     routine3d (&(array[0][0]),nx, Tab_Caval);
    }

   ki+=incx;
   if (ki<0 || ki>nx-1) 
    {
     free ((char *) tab1);
     free ((char *) array);
     goto END;
    }

   visible=(array[ny-1][ki]>=0.);
   
   deltax=-deltax;
   dx=deltax*(ki-1);
   x=dx+deltax;
   if (visible) move(x,array[ny-1][ki]);

   for (j=1;j<ny;j++)
    {
     kj=ny-j-1;
     it=xright=ki+1+nx*kj;
     xleft=xright+nx;
     x=dx+deltax*(j+1);
     routine3d (&(array[0][0]),nx, Tab_Caval);
    }

   x=dx+deltax*ny;
   if (!visible) move(x,array[0][ki]);

   if (!ki) 
    {
     free ((char *) tab1);
     free ((char *) array);
     goto END;
    }

   buf1=deltay*(ki+2);
   draw(x,buf1,Tab_Caval);
  
   ki+=incx;
   if (ki>nx-1) 
    {
     free ((char *) tab1);
     free ((char *) array);
     goto END;
    }
  }
 while (!ki);
 
 deltax=-deltax;
 dx=deltax*(-ki-ny);   
 x=dx+deltax;
 buf1=deltay*(ki+2);
 move(x,buf1);
 draw(x,array[0][ki],Tab_Caval);
 
 goto chat_huant;
 END:
 i=0;
}

/******************************************************************************/

void im_3dview_2 (Ifloat &Tab_Img, Ifloat &Tab_Caval, int Inc, int Mod_Visu)
{
    int j,Size,Indi,Indj;
    float i,Seuil, X_Pas,Increment,Nb_Line;
    float *Tab_Max;
    int Borne;
    int X1,Y1,X2,Y2;
    int Couleur;
    float Val,H_Facteur,L_Facteur,u,Interp;
    float Level;
    int Nl = Tab_Img.nl();
    int Nc = Tab_Img.nc();
    int Maxx = Tab_Caval.nc();
    int Maxy = Tab_Caval.nl();

    Size = Maxx * Maxy;    
    Tab_Max = new float[Maxx];
    for (j = 0; j < Maxx; j++) Tab_Max [j] = 0.;
    Tab_Caval.init();

    Increment  = (float) Inc;
    {
        float Min, Max;
 
        Min = Max = Tab_Img(0);
        for (j = 0; j < Nl*Nc; j++)
        {
            if (Min > Tab_Img(j)) Min = Tab_Img(j);
            if (Max < Tab_Img(j)) Max = Tab_Img(j);
        }
        for (j = 0; j < Nl*Nc; j++) Tab_Img(j) = (Tab_Img(j)-Min) / (Max - Min);
    }

    /* Seuil a partir duquel les lignes de l'image tab_caval1 passe
       a la seconde couleur */
    Seuil = 1. / 5.;
    
    /* parametre d'echelle compte tenu de l'angle de visualisation */
    H_Facteur = (float) Maxy / 2.;
    L_Facteur = (float) Maxx / (float) Nc;
    X_Pas = 1. / L_Facteur;
    Nb_Line = (float) Nl / ((float) Maxy / 2.);
    Borne = (int) ((Maxy /2. > Nl)? Nl : Maxy / 2.);
    
    //printf ("Borne = %d, Increment = %f, nbl = %f\n",Borne,Increment,Nb_Line);
    i = Borne + H_Facteur;
    //printf ("H_Facteur = %f, L_Facteur = %f, Limit = %f\n",H_Facteur,L_Facteur,i);
    
    for (i = 0; i < Borne; i += Increment)
    {
        /* Calcul de la position du premier point de la ligne */
        u = i * Nb_Line;
        Indi = (int) u * Nc;

        if (Indi > Nl*Nc)
        {
             printf ("(%f,%d)",Nb_Line,Indi);
        }
        Val = H_Facteur * Tab_Img(Indi) +  i;
        if (Val > Tab_Max [0]) Tab_Max [0] = Val;
        X1 = 0;
        Y1 = (int) Val;
        u = 0.;
        
        for (j = 1; j < Maxx; j ++)
        {
            /* Calcul de l'ordonnee (hauteur) */
            u += X_Pas;
            u = ((float) Nc / (float) Maxx) * (float) j;
            Indj = (int) u;
            
            if ((L_Facteur > 1.) && ((Indj + 1) < Nc))
            {
                Interp = u - Indj;
                Val  = Tab_Img (Indi+Indj);
                Val = H_Facteur * (Val + (Tab_Img (Indi+Indj+1) 
                                             - Val) * Interp) +  i;
                if ((Val >= Maxy) || (Val < 0.))
                {
                    printf ("zoom val (%d,%d)=%f ",Indi,Indj,Val);
                    exit (0);
                }
            }
            else
            {
                Val = H_Facteur * Tab_Img (Indi+Indj) + i;
                if ((Val >= Maxy) || (Val < 0.))
                {
                    printf ("calcul val %5.1f:(%d,%d)=%f ",i,Indi,Indj,Val);
                    exit (0);
                }
            }
                
            /* Tab_Max contient la position la plus haute 
               sur chaque colonne pour la gestion des faces cachees. 
               Si la nouvelle hauteur est superieure a celles deja traces
               alors on trace.
            */
            if (Val > Tab_Max [j])
            {
                Tab_Max [j] = Val;
                X2 = j;
               /* Y2 = Maxy - 1 - (Tab_Max [j] + 0.5);*/
                Y2 = (int) Val;
                if (Y2 < 0)
		{
                    printf ("Y2 < 0 : %d",Y2);
                    exit (0);
                    Y2 = 0;
		}
		if (Tab_Img (Indi+Indj) < Seuil) Couleur = VISU3D_LEVEL_1;
                else Couleur = VISU3D_LEVEL_2;
                
                if ((Couleur == VISU3D_LEVEL_1) 
                            && ((Tab_Caval (Y1*Maxx+X1) == VISU3D_LEVEL_2)
                                || (Tab_Caval (Y2*Maxx+X2) == VISU3D_LEVEL_2)))
                {
                    X1 = X1+1;
                    Y1 = Y2;
                }
                else
                {
                    Level = Tab_Img (Indi+Indj);
                    if (Y1 < 0)
		    {
                        printf ("Y1(%d,%5.1f)",j,Tab_Max [j]);
                        exit (-1);
                        Y1 = 0;
		    }
		    if (Mod_Visu == VISU_BLACK_WHITE)
		    {
		        Level = Couleur;
		    }
                    line (Y1,X1,Y2,X2,Level,Tab_Caval);
                    X1 = X2;
                    Y1 = Y2;
                }
            }
            else
            {
                X1 = j;
               /*  Y1 = Maxy - Tab_Max [j];*/
                Y1 = (int)(Tab_Max [j]);
                if (Y1 < 0)
		{
                    printf ("Y1(%d,%f)",j,Tab_Max [j]);
                    exit (-1);
                    Y1 = 0;
		}
            }
        }
    }
    
    delete [] Tab_Max;
}



