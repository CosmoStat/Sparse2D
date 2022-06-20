/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Pierre Djamdji && Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_FusionTool.cc
**
*******************************************************************************
**
**    DESCRIPTION  routines for image registration
**    -----------  
**                 
**
**
*******************************************************************************
**
** void copie_pts_fichier (char *nom_entree, float *x0, float *y0, 
**                    float *x1, float *y1, int nmax) 
** 
** char  *nom_entree  :  nom du fichier en sortie
** int    nmax   : nombre de maxima
** float *x0, *y0, *x1, *y1 : cvecteur a copier
** 
** Copie de vecteur de type float et de taille nmax dans le fichier nom_entree
** Les vecteurs sont stockes dans le fichier en :
** nmax sur la premiere ligne puis sequenciellement
** X0[i] et Y0[i] sur la premiere ligne et
** X1[i] et Y1[i] sur la seconde ligne
** et ainsi de suite pour tout les i, i=1..nmax
** 
*******************************************************************************
**
** void gauss (double **mat, double *vec1, int dimension, double *vvec2)
** 
** double **    mat          : matrice de valeur floatant
** double    *vec1, *vec2  : vecteurs en entree et resultat
** int       dimension     : dimension de la matrice
**
** Calcul d'un systeme d'equation par la methode de gauss avec recherche du
** pivot maximale en ligne et colonnes. 
** mat sera la matrice en entree, vec1 le vecteur en entree ou second membre
** du systeme d'equation et vec2 sera la solution recherchee
** 
******************************************************************************
** 
** void moindre_carre(int nbmax, 
**                    float *x_0,  float *y_0,  float *x_1,  float *y_1,
** 	              int type_equa,
**                    float *ax, float *bx, float *cx, float *dx, float *ex, 
**                    float *fx, float *gx, float *hx, float *ix, float *jx,
**                    float *ay, float *by, float *cy, float *dy, float *ey, 
**                    float *fy, float *gy, float *hy, float *iy, float *jy  )
** 
** int    nbmax  : nombre de max
** float  *x_0, *y_0, *x_1, *y_1  : vecteur de points en x et y 
**                                  respectivement pour les variables
**                                  X,Y et x',y'
** int    type_equa               : type de l'equation
** float  *ax, *bx, *cx, *dx, *ex, *fx, *gx, *hx, *ix, *jx,
**        *ay, *by, *cy, *dy, *ey, *fy, *gy, *hy, *iy, *jy   
**                                : parametre des equations polynomiales
**                                  respectivement en x et y
**
** Calcul des parametres d'equations polynomiales d'ordre 1, 2 et 3
** par moindre carre et par l'algorithme de gauss.
** l'equation pour le type 1 s'ecrit :
** x' = aX + bY + c
** y' = dX + eY + f
** 
******************************************************************************
** 
** void moindre_carre_f(int nbmax, 
**                      float *x_0,  float *y_0,  float *x_1,  float *y_1,
** 	                int type_equa,
**                      float *ax, float *bx, float *cx, float *dx, float *ex, 
**                      float *fx, float *gx, float *hx, float *ix, float *jx,
**                      float *ay, float *by, float *cy, float *dy, float *ey, 
**                      float *fy, float *gy, float *hy, float *iy, float *jy  )
**
** int    nbmax  : nombre de max
** float  *x_0, *y_0, *x_1, *y_1  : vecteur de points en x et y 
**                                  respectivement pour les variables
**                                  X,Y et x',y'
** int    type_equa               : type de l'equation
** float  *ax, *bx, *cx, *dx, *ex, *fx, *gx, *hx, *ix, *jx,
**        *ay, *by, *cy, *dy, *ey, *fy, *gy, *hy, *iy, *jy   
**                                : parametre des equations polynomiales
**                                  respectivement en x et y
**
** Calcul des parametres d'equations polynomiales d'ordre 1, 2 et 3
** par moindre carre et par l'algorithme de gauss.
** l'equation pour le type 1 s'ecrit :
** x' = aX + bY + c
** y' = dX + eY + f
** 
** NOTE : identique a moindre_carre.c
**
******************************************************************************
**
** void statis (float *tab_in, int nb_pts, char *name_mouch, char *name_tab)
**
** char   *name_mouch, *name_tab : nom du ficher mouchard et
**                                 du tableau en traitement
** float  *tab_in         : tableau de donnees en entree (vecteur)
** int    nb_pts          : nbre de points de donnees
** 
**
** Calul de statistiques elementaires sur une vecteur de donnes de type floatant
** Les donnnes en sortie sont : le min, le max, la moyenne, la variance et
** l'ecart-type. Ces donnes sont seulement ecrite sur le mouchard ainsi que
** le nom du tableau de valeur en traitement
**
******************************************************************************
**
** void transpo_coor_max(float *x_e, float *y_e, 
**                      parametre *par,
**                      int nmax, int type_equa,
**		      float **x_s, float **y_s)
**
** int        nmax, type_equa : nbre de max et type d'equation
** float      *x_e, *y_e      : vecteur en entree a transformee
** tableau_r  *x_s, *y_s      : vecteur en sortie resultant de la transformation
** parametre  *par            : valeur des parametres du modele polynomiale 
**                              utilise
** 
** Transformation des valeurs des vecteurs x_e et y_e de type float en entree
** a l'aide d'un modele polynomiale de degree 1, 2 ou 3. Les parametres du
** modele sont donnees par la structure par. Le resultat sera stocke dans
** x_s et y_s respectivement pour x_e et y_e.
** 
******************************************************************************
** 
** void trans_coor_pxl_reel_tab(int depart_x, int depart_y, 
**                             float pas_x, float pas_y,
**                             float *tab01_x, float *tab01_y , 
**                             float *tab02_x, float *tab02_y,
**                             float **tab11_x, float **tab11_y,
**                             float **tab12_x, float **tab12_y,
**                             int nmax  ) 
** 
** int        depart_x, depart_y  : coordonnees de depart en x et y
** int        nmax                : nombre de points
** float      *tab01_x, *tab01_y , 
**            *tab02_x, *tab02_y  : vecteur en coordonnees pixel a transformer
** float      pas_x, pas_y        : pas d'echantillonnage en x et y respecti-
**                                  vement
** tableau_r  *tab11_x, *tab11_y, 
**            *tab12_x, *tab12_y  : vecteur en coordonnees reelles resultant
**                                  de la transformation
** 
** 
** Transformation de vecteurs de type floatant de coordonnees pixels en
** coordonnees reelles. La transformation s'effectue sur deux systemes
** de vecteurs en X et Y
** les coordonnees de depart doivent etre en coordonnees pixel
** les coordonnees reelles de l'images d'origine seront considerees comme
** etant egale a 1. en x et y
**
******************************************************************************/

// static char sccsid[] = "@(#)MR_FusionTool.cc 3.1 96/05/02 CEA 1996 @(#)";

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Fusion.h"


/***************************************************************************/

void trans_coor_pxl_reel_tab(int depart_x, int depart_y, 
                             float pas_x, float pas_y,
                             float *tab01_x, float *tab01_y , 
                             float *tab02_x, float *tab02_y,
                             float **tab11_x, float **tab11_y,
                             float **tab12_x, float **tab12_y,
                             int nmax  ) 
{
 
 int        i ;
 float      depart_x_cr, depart_y_cr ;

 /* ----- allocation dynamique des tableaux en sortie ------ */

   *tab11_x = new float[nmax] ;
   *tab11_y = new float[nmax] ;
   *tab12_x = new float[nmax] ;
   *tab12_y = new float[nmax] ;


 /* ----- allocation dynamique des tableaux en sortie ------ */

/* transformation des coordonnees depart en coordonnees depart reelles */

   depart_x_cr = 1. + ((float)depart_x * pas_x) ;
   depart_y_cr = 1. + ((float)depart_y * pas_y) ;


   for (i = 0 ; i < nmax ; i++)
     {
      /* transformation des coordonnees pixels en coordonnees reelles */
         (*tab11_x)[i] = (float)depart_x_cr + (tab01_x[i] * pas_x) ;
         (*tab11_y)[i] = (float)depart_y_cr + (tab01_y[i] * pas_y) ;
         (*tab12_x)[i] = (float)depart_x_cr + (tab02_x[i] * pas_x) ;
         (*tab12_y)[i] = (float)depart_y_cr + (tab02_y[i] * pas_y) ;  
      }
}
/* end */

/***************************************************************************/

void transpo_coor_max(float *x_e, float *y_e, parametre *par,
                      int nmax, int type_equa, float **x_s, float **y_s)

{
  int    i;
  double   col, col2, ligne, ligne2, col3, ligne3 ;

  /* ----------   allocation dynamique   -------------- */

   *x_s = new float[nmax] ;
   *y_s = new float[nmax] ;

  for (i = 0; i < nmax; i++) 
  {
      switch(type_equa)
      {
       case T_EQUA_0:
          (*x_s)[i] = par->ax * x_e[i] - par->bx * y_e[i] + par->cx;
          (*y_s)[i] = par->bx * x_e[i] + par->ax * y_e[i] + par->cy;
          break;
       case T_EQUA_1:
          (*x_s)[i] = par->ax * x_e[i] + par->bx * y_e[i] + par->cx;
          (*y_s)[i] = par->ay * x_e[i] + par->by * y_e[i] + par->cy;
          break;
       case T_EQUA_2:
          col = x_e[i];
          col2 = x_e[i] * x_e[i];
          ligne = y_e[i];
          ligne2 = y_e[i] * y_e[i];

          (*x_s)[i] = par->ax * col2 + par->bx * ligne2 + par->cx * ligne * col + 
                  par->dx * col + par->ex * ligne + par->fx;
          (*y_s)[i] = par->ay * col2 + par->by * ligne2 + par->cy * ligne * col + 
                  par->dy * col + par->ey * ligne + par->fy;
          break;
       case T_EQUA_3:
          col  = x_e[i];
          col2 = x_e[i] * x_e[i];
          col3 = col2   * x_e[i];
          ligne = y_e[i];
          ligne2 = y_e[i] * y_e[i];
          ligne3 = ligne2 * y_e[i];

          (*x_s)[i] = par->ax * col3 + par->bx * ligne3 + par->cx * col2 + 
                  par->dx * ligne2 + par->ex * col * ligne2 + 
                  par->fx * col2 * ligne + par->gx * col * ligne + 
                  par->hx * col + par->ix * ligne + par->jx ;
          (*y_s)[i] = par->ay * col3 + par->by * ligne3 + par->cy * col2 + 
                  par->dy * ligne2 + par->ey * col * ligne2 +
                  par->fy * col2 * ligne + par->gy * col * ligne + 
                  par->hy * col + par->iy * ligne + par->jy; 
          break;
       default:
         cerr << "Error: bad type of deformation model ..." << endl;
         exit(-1);
         break;
     }
  }

}

/***************************************************************************/

void moindre_carre_f(int nbmax, 
                     float *x_0,  float *y_0,  float *x_1,  float *y_1,
	             int type_equa,
                     float *ax, float *bx, float *cx, float *dx, float *ex, 
                     float *fx, float *gx, float *hx, float *ix, float *jx,
                     float *ay, float *by, float *cy, float *dy, float *ey, 
                     float *fy, float *gy, float *hy, float *iy, float *jy  )
{
  double   zx1, zx0, zx02, zx03, zx04, zx05, zx06,
           zy1, zy0, zy02, zy03, zy04, zy05, zy06  ;
  double   msx0, msx02, msx03, msx04, msx05, msx06, 
           msy0, msy02, msy03, msy04, msy05, msy06,
           msx0_y0, msx0_y02, msx0_y03, msx0_y04, msx0_y05, 
	   msx02_y0, msx02_y02, msx02_y03, msx02_y04,
           msx03_y0, msx03_y02, msx03_y03,
           msx04_y0, msx04_y02,
           msx05_y0,
           msx1, msx1_x0, msx1_y0, msx1_x02, msx1_y02, msx1_x03, msx1_y03, 
           msx1_x0_y0, msx1_x02_y0, msx1_x0_y02,
           msy1, msy1_x0, msy1_y0, msy1_x02, msy1_y02, msy1_x03, msy1_y03, 
           msy1_x0_y0, msy1_x02_y0, msy1_x0_y02 ;



  double    **a0, **a1 ;

  double  *b0, *b1, *xs, *xs1;
  int      i, j, nb_param=0 ;


  // calcul du nombre de parametres necessaires pour l'allocation dynamique
  //  (n+1) * (n+2) /2   n etant le degre du polynome   

  if (type_equa == T_EQUA_0)
     nb_param = 3 ;
  if (type_equa == T_EQUA_1)
     nb_param = 3 ;
  if (type_equa == T_EQUA_2)
     nb_param = 6 ;
  if (type_equa == T_EQUA_3)
     nb_param = 10 ;


  /* -----------     allocation dynamique    ----------- */
   b0  = new double[nb_param] ;
   b1  = new double[nb_param] ;
   xs  = new double[nb_param] ;
   xs1 = new double[nb_param] ;

   /* allocation a deux dimensions */
  a0 = new double* [nb_param] ;
  a1 = new double* [nb_param] ;
  for (i = 0; i < nb_param; i++) 
  { 
      a0[i] = new double[nb_param] ;
      a1[i] = new double[nb_param] ;
  }

 /* -----------  fin de l'allocation dynamique  ----------- */

  msx0  = 0.0;
  msx02 = 0.0;
  msx03 = 0.0;
  msx04 = 0.0;
  msx05 = 0.0;
  msx06 = 0.0;

  msy0  = 0.0;
  msy02 = 0.0;
  msy03 = 0.0;
  msy04 = 0.0;
  msy05 = 0.0;
  msy06 = 0.0;

  msx0_y0   = 0.0;
  msx0_y02  = 0.0;
  msx0_y03  = 0.0;
  msx0_y04  = 0.0;
  msx0_y05  = 0.0;

  msx02_y0  = 0.0;
  msx02_y02 = 0.0;
  msx02_y03 = 0.0;
  msx02_y04 = 0.0;

  msx03_y0  = 0.0;
  msx03_y02 = 0.0;
  msx03_y03 = 0.0;

  msx04_y0  = 0.0;
  msx04_y02 = 0.0;

  msx05_y0  = 0.0;

  msx1       = 0.0;
  msx1_x0    = 0.0;
  msx1_y0    = 0.0;
  msx1_x02   = 0.0;
  msx1_y02   = 0.0;
  msx1_x03   = 0.0;
  msx1_y03   = 0.0;
  msx1_x0_y0  = 0.0;
  msx1_x02_y0 = 0.0;
  msx1_x0_y02 = 0.0;

  msy1       = 0.0;
  msy1_x0    = 0.0;
  msy1_y0    = 0.0;
  msy1_x02   = 0.0;
  msy1_y02   = 0.0;
  msy1_x03   = 0.0;
  msy1_y03   = 0.0;
  msy1_x0_y0 = 0.0;
  msy1_x02_y0 = 0.0;
  msy1_x0_y02 = 0.0;



  *ax = 0 ;
  *bx = 0 ;
  *cx = 0 ;
  *dx = 0 ;
  *ex = 0 ;
  *fx = 0 ;
  *gx = 0 ;
  *hx = 0 ;
  *ix = 0 ;
  *jx = 0 ;

  *ay = 0 ;
  *by = 0 ;
  *cy = 0 ;
  *dy = 0 ;
  *ey = 0 ;
  *fy = 0 ; 
  *gy = 0 ;
  *hy = 0 ;
  *iy = 0 ;
  *jy = 0 ; 

  /* ---  cas polynomiale de degre 1 de type I  ...  lineaire  --- */

  if (type_equa == T_EQUA_0) 
  {
    for (i = 0; i < nbmax; i++)
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];
 
      msx1 += zx1;
      msy1 += zy1;
      msx0 += zx0;
      msy0 += zy0;
      msx1_x0 += zx0 * zx1;
      msy1_y0 += zy0 * zy1;
      msy1_x0 += zx0 * zy1;
      msx1_y0 += zx1 * zy0;
      msx0_y0 += zx0 * zy0;
      msx02 += zx0 * zx0;
      msy02 += zy0 * zy0;
    }



    a0[0][0] = msy02  ;
    a0[0][1] = msx0_y0;
    a0[0][2] = msy0   ;

    a0[1][0] = msx0_y0;
    a0[1][1] = msx02  ;
    a0[1][2] = msx0   ;

    a0[2][0] = msy0   ;
    a0[2][1] = msx0   ;
    a0[2][2] = nbmax  ;


    b0[0]    = msy1_y0;
    b0[1]    = msy1_x0;
    b0[2]    = msy1   ;




    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cy                ----- */
    /* ------------------------------------------------------------------ */




    gauss(a0, b0, 3, xs);




    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cy = (float)xs[2];


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cx                ----- */
    /* ------------------------------------------------------------------ */


    a0[0][0] =  msx02   ;
    a0[0][1] = -msx0_y0 ;
    a0[0][2] =  msx0    ;
    b0[0]    =  msx1_x0 ;

    a0[1][0] =  msx0_y0 ;
    a0[1][1] = -msy02   ;
    a0[1][2] =  msy0    ;
    b0[1]    =  msx1_y0 ;

    a0[2][0] =  msy0    ;
    a0[2][1] =  msx0    ;
    a0[2][2] =  nbmax   ;
    b0[2]    =  msx1    ;

    gauss(a0, b0, 3, xs);

    *cx = xs[2];
    if (xs[0] >= *ax)  *ax = xs[0];
    if (xs[1] >= *bx)  *bx = xs[2];
    goto the_end ;

  }

  /* -- cas polynomiale de degre 1 de type II  ...  lineaire -- */

  if (type_equa == T_EQUA_1)
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];


      msx1 += zx1;
      msy1 += zy1;
      msx0 += zx0;
      msy0 += zy0;
      msx1_x0 += zx0 * zx1;
      msy1_y0 += zy0 * zy1;
      msy1_x0 += zx0 * zy1;
      msx1_y0 += zx1 * zy0;
      msx0_y0 += zx0 * zy0;
      msx02 += zx0 * zx0;
      msy02 += zy0 * zy0;
    }



    a0[0][0] = msx02   ;
    a0[0][1] = msx0_y0 ;
    a0[0][2] = msx0    ;

    a0[1][0] = msx0_y0 ;
    a0[1][1] = msy02   ;
    a0[1][2] = msy0    ;
 
    a0[2][0] = msx0    ;
    a0[2][1] = msy0    ;
    a0[2][2] = nbmax   ;



    b0[0]    = msx1_x0 ;
    b0[1]    = msx1_y0 ;
    b0[2]    = msx1    ;

    b1[0]    = msy1_x0 ;
    b1[1]    = msy1_y0 ;
    b1[2]    = msy1    ;


    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cx                ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ay, by, et cy                ----- */
    /* ------------------------------------------------------------------ */

    gauss(a1, b1, nb_param, xs1);



    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];

    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];


    goto the_end;

  }

  /* - cas polynomiale de degre 2   ...   non lineaire   -- */

  if (type_equa == T_EQUA_2) 
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];

      zx02 = zx0 * zx0;
      zx03 = zx02 * zx0;
      zx04 = zx03 * zx0;
      zy02 = zy0 * zy0;
      zy03 = zy02 * zy0;
      zy04 = zy03 * zy0;

      msx0 += zx0;
      msx02 += zx02;
      msx03 += zx03;
      msx04 += zx04;

      msy0 += zy0;
      msy02 += zy02;
      msy03 += zy03;
      msy04 += zy04;

      msx0_y0 += zx0 * zy0;
      msx02_y0 += zx02 * zy0;
      msx03_y0 += zx03 * zy0;
      msx0_y02 += zx0 * zy02;
      msx0_y03 += zx0 * zy03;
      msx02_y02 += zx02 * zy02;

      msx1 += zx1;
      msx1_x0 += zx1 * zx0;
      msx1_y0 += zx1 * zy0;
      msx1_x02 += zx1 * zx02;
      msx1_y02 += zx1 * zy02;
      msx1_x0_y0 += zx1 * zx0 * zy0;

      msy1 += zy1;
      msy1_x0 += zy1 * zx0;
      msy1_y0 += zy1 * zy0;
      msy1_x02 += zy1 * zx02;
      msy1_y02 += zy1 * zy02;
      msy1_x0_y0 += zy1 * zx0 * zy0;

    }


    a0[0][0] = msx04;
    a0[0][1] = msx02_y02;
    a0[0][2] = msx03_y0;
    a0[0][3] = msx03;
    a0[0][4] = msx02_y0;
    a0[0][5] = msx02;

    a0[1][0] = a0[0][1];
    a0[1][1] = msy04;
    a0[1][2] = msx0_y03;
    a0[1][3] = msx0_y02;
    a0[1][4] = msy03;
    a0[1][5] = msy02;

    a0[2][0] = a0[0][2];
    a0[2][1] = a0[1][2];
    a0[2][2] = msx02_y02;
    a0[2][3] = msx02_y0;
    a0[2][4] = msx0_y02;
    a0[2][5] = msx0_y0;

    a0[3][0] = a0[0][3];
    a0[3][1] = a0[1][3];
    a0[3][2] = a0[2][3];
    a0[3][3] = msx02;
    a0[3][4] = msx0_y0;
    a0[3][5] = msx0;

    a0[4][0] = a0[0][4];
    a0[4][1] = a0[1][4];
    a0[4][2] = a0[2][4];
    a0[4][3] = a0[3][4];
    a0[4][4] = msy02;
    a0[4][5] = msy0;

    a0[5][0] = a0[0][5];
    a0[5][1] = a0[1][5];
    a0[5][2] = a0[2][5];
    a0[5][3] = a0[3][5];
    a0[5][4] = a0[4][5];
    a0[5][5] = nbmax;



    b0[0] = msx1_x02;
    b0[1] = msx1_y02;
    b0[2] = msx1_x0_y0;
    b0[3] = msx1_x0;
    b0[4] = msx1_y0;
    b0[5] = msx1;

    b1[0] = msy1_x02;
    b1[1] = msy1_y02;
    b1[2] = msy1_x0_y0;
    b1[3] = msy1_x0;
    b1[4] = msy1_y0;
    b1[5] = msy1;


    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax,bx,cx,dx,ex,fx,gx,hx      ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);



    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ay,by,cy,dy,ey,fy            ----- */
    /* ------------------------------------------------------------------ */

    gauss(a1, b1, nb_param, xs1);



    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];
    *dx = (float)xs[3];
    *ex = (float)xs[4];
    *fx = (float)xs[5];

    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];
    *dy = (float)xs1[3];
    *ey = (float)xs1[4];
    *fy = (float)xs1[5];


    goto the_end;

  }

  /* the_end du polynome de degre 2  */


  /* -- cas polynomiale de degre 3   ...   non lineaire  -- */


  if (type_equa == T_EQUA_3) 
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];

      zx02 = zx0  * zx0;
      zx03 = zx02 * zx0;
      zx04 = zx03 * zx0;
      zx05 = zx04 * zx0;
      zx06 = zx05 * zx0;

      zy02 = zy0  * zy0;
      zy03 = zy02 * zy0;
      zy04 = zy03 * zy0;
      zy05 = zy04 * zy0;
      zy06 = zy05 * zy0;


      msx0  += zx0;
      msx02 += zx02;
      msx03 += zx03;
      msx04 += zx04;
      msx05 += zx05;
      msx06 += zx06;

      msy0  += zy0;
      msy02 += zy02;
      msy03 += zy03;
      msy04 += zy04;
      msy05 += zy05;
      msy06 += zy06;

      msx0_y0    += zx0 * zy0;
      msx0_y02   += zx0 * zy02;
      msx0_y03   += zx0 * zy03;
      msx0_y04   += zx0 * zy04;
      msx0_y05   += zx0 * zy05;

      msx02_y0   += zx02 * zy0;
      msx02_y02  += zx02 * zy02;
      msx02_y03  += zx02 * zy03;
      msx02_y04  += zx02 * zy04;

      msx03_y0   += zx03 * zy0;
      msx03_y02  += zx03 * zy02;
      msx03_y03  += zx03 * zy03;

      msx04_y0   += zx04 * zy0;
      msx04_y02  += zx04 * zy02;

      msx05_y0   += zx05 * zy0;

      msx1 += zx1;
      msx1_x0 += zx1 * zx0;
      msx1_y0 += zx1 * zy0;
      msx1_x02 += zx1 * zx02;
      msx1_y02 += zx1 * zy02;
      msx1_x03 += zx1 * zx03;
      msx1_y03 += zx1 * zy03;

      msx1_x0_y0  += zx1 * zx0  * zy0;
      msx1_x02_y0 += zx1 * zx02 * zy0;
      msx1_x0_y02 += zx1 * zx0  * zy02;

      msy1     += zy1;
      msy1_x0  += zy1 * zx0;
      msy1_y0  += zy1 * zy0;
      msy1_x02 += zy1 * zx02;
      msy1_y02 += zy1 * zy02;
      msy1_x03 += zy1 * zx03;
      msy1_y03 += zy1 * zy03;

      msy1_x0_y0 += zy1 * zx0 * zy0;
      msy1_x02_y0 += zy1 * zx02 * zy0;
      msy1_x0_y02 += zy1 * zx0  * zy02;


    }




    a0[0][0] = msx06;
    a0[0][1] = msx03_y03;
    a0[0][2] = msx05;
    a0[0][3] = msx03_y02;
    a0[0][4] = msx04_y02;
    a0[0][5] = msx05_y0;
    a0[0][6] = msx04_y0;
    a0[0][7] = msx04;
    a0[0][8] = msx03_y0;
    a0[0][9] = msx03;

    a0[1][0] = a0[0][1];
    a0[1][1] = msy06;
    a0[1][2] = msx02_y03;
    a0[1][3] = msy05;
    a0[1][4] = msx0_y05;
    a0[1][5] = msx02_y04;
    a0[1][6] = msx0_y04;
    a0[1][7] = msx0_y03;
    a0[1][8] = msy04;
    a0[1][9] = msy03;

    a0[2][0] = a0[0][2];
    a0[2][1] = a0[1][2];
    a0[2][2] = msx04;
    a0[2][3] = msx02_y02;
    a0[2][4] = msx03_y02;
    a0[2][5] = msx04_y0;
    a0[2][6] = msx03_y0;
    a0[2][7] = msx03;
    a0[2][8] = msx02_y0;
    a0[2][9] = msx02;

    a0[3][0] = a0[0][3];
    a0[3][1] = a0[1][3];
    a0[3][2] = a0[2][3];
    a0[3][3] = msy04;
    a0[3][4] = msx0_y04;
    a0[3][5] = msx02_y03;
    a0[3][6] = msx0_y03;
    a0[3][7] = msx0_y02;
    a0[3][8] = msy03;
    a0[3][9] = msy02;

    a0[4][0] = a0[0][4];
    a0[4][1] = a0[1][4];
    a0[4][2] = a0[2][4];
    a0[4][3] = a0[3][4];
    a0[4][4] = msx02_y04;
    a0[4][5] = msx03_y03;
    a0[4][6] = msx02_y03;
    a0[4][7] = msx02_y02;
    a0[4][8] = msx0_y03;
    a0[4][9] = msx0_y02;

    a0[5][0] = a0[0][5];
    a0[5][1] = a0[1][5];
    a0[5][2] = a0[2][5];
    a0[5][3] = a0[3][5];
    a0[5][4] = a0[4][5];
    a0[5][5] = msx04_y02;
    a0[5][6] = msx03_y02;
    a0[5][7] = msx03_y0;
    a0[5][8] = msx02_y02;
    a0[5][9] = msx02_y0;

    a0[6][0] = a0[0][6];
    a0[6][1] = a0[1][6];
    a0[6][2] = a0[2][6];
    a0[6][3] = a0[3][6];
    a0[6][4] = a0[4][6];
    a0[6][5] = a0[5][6];
    a0[6][6] = msx02_y02;
    a0[6][7] = msx02_y0;
    a0[6][8] = msx0_y02;
    a0[6][9] = msx0_y0;

    a0[7][0] = a0[0][7];
    a0[7][1] = a0[1][7];
    a0[7][2] = a0[2][7];
    a0[7][3] = a0[3][7];
    a0[7][4] = a0[4][7];
    a0[7][5] = a0[5][7];
    a0[7][6] = a0[6][7];
    a0[7][7] = msx02;
    a0[7][8] = msx0_y0;
    a0[7][9] = msx0;

    a0[8][0] = a0[0][8];
    a0[8][1] = a0[1][8];
    a0[8][2] = a0[2][8];
    a0[8][3] = a0[3][8];
    a0[8][4] = a0[4][8];
    a0[8][5] = a0[5][8];
    a0[8][6] = a0[6][8];
    a0[8][7] = a0[7][8];
    a0[8][8] = msy02;
    a0[8][9] = msy0;

    a0[9][0] = a0[0][9];
    a0[9][1] = a0[1][9];
    a0[9][2] = a0[2][9];
    a0[9][3] = a0[3][9];
    a0[9][4] = a0[4][9];
    a0[9][5] = a0[5][9];
    a0[9][6] = a0[6][9];
    a0[9][7] = a0[7][9];
    a0[9][8] = a0[8][9];
    a0[9][9] = nbmax;



    b0[0] = msx1_x03;
    b0[1] = msx1_y03;
    b0[2] = msx1_x02;
    b0[3] = msx1_y02;
    b0[4] = msx1_x0_y02;
    b0[5] = msx1_x02_y0;
    b0[6] = msx1_x0_y0;
    b0[7] = msx1_x0;
    b0[8] = msx1_y0;
    b0[9] = msx1;

    b1[0] = msy1_x03;
    b1[1] = msy1_y03;
    b1[2] = msy1_x02;
    b1[3] = msy1_y02;
    b1[4] = msy1_x0_y02;
    b1[5] = msy1_x02_y0;
    b1[6] = msy1_x0_y0;
    b1[7] = msy1_x0;
    b1[8] = msy1_y0;
    b1[9] = msy1;




    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }


    /* ------------------------------------------------------------------ */
    /* -----      resolution de : ax,bx,cx,dx,ex,fx,gx,hx,ix,jx     ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);


    /* ------------------------------------------------------------------ */
    /* -----    resolution de : ay,by,cy,dy,ey,fy,gy,hy,iy,jy       ----- */
    /* ------------------------------------------------------------------ */
  
    gauss(a1, b1, nb_param, xs1);


    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];
    *dx = (float)xs[3];
    *ex = (float)xs[4];
    *fx = (float)xs[5];
    *gx = (float)xs[6];
    *hx = (float)xs[7];
    *ix = (float)xs[8];
    *jx = (float)xs[9];


    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];
    *dy = (float)xs1[3];
    *ey = (float)xs1[4];
    *fy = (float)xs1[5];
    *gy = (float)xs1[6];
    *hy = (float)xs1[7];
    *iy = (float)xs1[8];
    *jy = (float)xs1[9];


    goto the_end;

  }

  /* the_end du polynome de degre 3  */



  the_end: ;




  /* ------ liberation de l'espace memoire ------ */


  delete [] xs1 ;
  delete [] xs ;
  delete [] b1 ;
  delete [] b0 ;

      /* ------ liberation de l'espace memoire 2D ------ */

 for (i = 0; i < nb_param; i++)
   { 
     delete [] a0[i] ;
   }
 delete [] a0 ;

 for (i = 0; i < nb_param; i++)
   { 
     delete [] a1[i] ;
   }
 delete [] a1 ;


  /* ------ liberation de l'espace memoire ------ */

/* End. */

}


/***********************************************************************/


void moindre_carre(int nbmax, 
                   float *x_0,  float *y_0,  float *x_1,  float *y_1,
	           int type_equa,
                   float *ax, float *bx, float *cx, float *dx, float *ex, 
                   float *fx, float *gx, float *hx, float *ix, float *jx,
                   float *ay, float *by, float *cy, float *dy, float *ey, 
                   float *fy, float *gy, float *hy, float *iy, float *jy  )

{
  double   zx1, zx0, zx02, zx03, zx04, zx05, zx06,
           zy1, zy0, zy02, zy03, zy04, zy05, zy06  ;
  double   msx0, msx02, msx03, msx04, msx05, msx06, 
           msy0, msy02, msy03, msy04, msy05, msy06,
           msx0_y0, msx0_y02, msx0_y03, msx0_y04, msx0_y05, 
	   msx02_y0, msx02_y02, msx02_y03, msx02_y04,
           msx03_y0, msx03_y02, msx03_y03,
           msx04_y0, msx04_y02,
           msx05_y0,
           msx1, msx1_x0, msx1_y0, msx1_x02, msx1_y02, msx1_x03, msx1_y03, 
           msx1_x0_y0, msx1_x02_y0, msx1_x0_y02,
           msy1, msy1_x0, msy1_y0, msy1_x02, msy1_y02, msy1_x03, msy1_y03, 
           msy1_x0_y0, msy1_x02_y0, msy1_x0_y02 ;




  double  **a1, **a0; 
  double  *b0, *b1, *xs, *xs1;
  float   *x00, *x10, *y00, *y10  ;
  float   *x11, *y11 ;
  int      i, j, nb_param=0 ;


  /* calcul du nombre de parametres necessaires pour l'allocation dynamique */
  if (type_equa == T_EQUA_0)
     nb_param = 3 ;
  if (type_equa == T_EQUA_1)
     nb_param = 3 ;
  if (type_equa == T_EQUA_2)
     nb_param = 6 ;
  if (type_equa == T_EQUA_3)
     nb_param = 10 ;


  /* -----------     allocation dynamique    ----------- */

   b0  = new double[nb_param] ;
   b1  = new double[nb_param] ;
   xs  = new double[nb_param] ;
   xs1 = new double[nb_param] ;



   x00 = new float[2000] ;
   x10 = new float[2000] ;
   y00 = new float[2000] ;
   y10 = new float[2000] ;


   x11 = new float[2000] ;
   y11 = new float[2000] ;





            /* allocation a deux dimensions */


 a0 = new double* [nb_param] ;
 for (i = 0; i < nb_param; i++)
   { 
     a0[i] = new double[nb_param] ;
   }

 a1 = new double* [nb_param] ;
 for (i = 0; i < nb_param; i++)
   { 
     a1[i] = new double[nb_param] ;
   }


 /* -----------  fin de l'allocation dynamique  ----------- */


  msx0  = 0.0;
  msx02 = 0.0;
  msx03 = 0.0;
  msx04 = 0.0;
  msx05 = 0.0;
  msx06 = 0.0;

  msy0  = 0.0;
  msy02 = 0.0;
  msy03 = 0.0;
  msy04 = 0.0;
  msy05 = 0.0;
  msy06 = 0.0;

  msx0_y0   = 0.0;
  msx0_y02  = 0.0;
  msx0_y03  = 0.0;
  msx0_y04  = 0.0;
  msx0_y05  = 0.0;

  msx02_y0  = 0.0;
  msx02_y02 = 0.0;
  msx02_y03 = 0.0;
  msx02_y04 = 0.0;

  msx03_y0  = 0.0;
  msx03_y02 = 0.0;
  msx03_y03 = 0.0;

  msx04_y0  = 0.0;
  msx04_y02 = 0.0;

  msx05_y0  = 0.0;

  msx1       = 0.0;
  msx1_x0    = 0.0;
  msx1_y0    = 0.0;
  msx1_x02   = 0.0;
  msx1_y02   = 0.0;
  msx1_x03   = 0.0;
  msx1_y03   = 0.0;
  msx1_x0_y0  = 0.0;
  msx1_x02_y0 = 0.0;
  msx1_x0_y02 = 0.0;

  msy1       = 0.0;
  msy1_x0    = 0.0;
  msy1_y0    = 0.0;
  msy1_x02   = 0.0;
  msy1_y02   = 0.0;
  msy1_x03   = 0.0;
  msy1_y03   = 0.0;
  msy1_x0_y0 = 0.0;
  msy1_x02_y0 = 0.0;
  msy1_x0_y02 = 0.0;



  *ax = 0 ;
  *bx = 0 ;
  *cx = 0 ;
  *dx = 0 ;
  *ex = 0 ;
  *fx = 0 ;
  *gx = 0 ;
  *hx = 0 ;
  *ix = 0 ;
  *jx = 0 ;

  *ay = 0 ;
  *by = 0 ;
  *cy = 0 ;
  *dy = 0 ;
  *ey = 0 ;
  *fy = 0 ; 
  *gy = 0 ;
  *hy = 0 ;
  *iy = 0 ;
  *jy = 0 ; 


  /*  cas polynomiale de degre 1 de type I  ...  lineaire  -- */


  if (type_equa == T_EQUA_0) 
  {
    for (i = 0; i < nbmax; i++)
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];
 
      msx1 += zx1;
      msy1 += zy1;
      msx0 += zx0;
      msy0 += zy0;
      msx1_x0 += zx0 * zx1;
      msy1_y0 += zy0 * zy1;
      msy1_x0 += zx0 * zy1;
      msx1_y0 += zx1 * zy0;
      msx0_y0 += zx0 * zy0;
      msx02 += zx0 * zx0;
      msy02 += zy0 * zy0;
    }



    a0[0][0] = msy02  ;
    a0[0][1] = msx0_y0;
    a0[0][2] = msy0   ;

    a0[1][0] = msx0_y0;
    a0[1][1] = msx02  ;
    a0[1][2] = msx0   ;

    a0[2][0] = msy0   ;
    a0[2][1] = msx0   ;
    a0[2][2] = nbmax  ;


    b0[0]    = msy1_y0;
    b0[1]    = msy1_x0;
    b0[2]    = msy1   ;




    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cy                ----- */
    /* ------------------------------------------------------------------ */


    gauss(a0, b0, 3, xs);


    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cy = (float)xs[2];


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cx                ----- */
    /* ------------------------------------------------------------------ */


    a0[0][0] =  msx02   ;
    a0[0][1] = -msx0_y0 ;
    a0[0][2] =  msx0    ;
    b0[0]    =  msx1_x0 ;

    a0[1][0] =  msx0_y0 ;
    a0[1][1] = -msy02   ;
    a0[1][2] =  msy0    ;
    b0[1]    =  msx1_y0 ;

    a0[2][0] =  msy0    ;
    a0[2][1] =  msx0    ;
    a0[2][2] =  nbmax   ;
    b0[2]    =  msx1    ;


    gauss(a0, b0, 3, xs);

 
    *cx = xs[2];
    if (xs[0] >= *ax)
      *ax = xs[0];

    if (xs[1] >= *bx)
      *bx = xs[2];

    goto the_end ;

  }


  /* -cas polynomiale de degre 1 de type II  ...  lineaire  - */



  if (type_equa == T_EQUA_1)
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];


      msx1 += zx1;
      msy1 += zy1;
      msx0 += zx0;
      msy0 += zy0;
      msx1_x0 += zx0 * zx1;
      msy1_y0 += zy0 * zy1;
      msy1_x0 += zx0 * zy1;
      msx1_y0 += zx1 * zy0;
      msx0_y0 += zx0 * zy0;
      msx02 += zx0 * zx0;
      msy02 += zy0 * zy0;
    }



    a0[0][0] = msx02   ;
    a0[0][1] = msx0_y0 ;
    a0[0][2] = msx0    ;

    a0[1][0] = msx0_y0 ;
    a0[1][1] = msy02   ;
    a0[1][2] = msy0    ;
 
    a0[2][0] = msx0    ;
    a0[2][1] = msy0    ;
    a0[2][2] = nbmax   ;



    b0[0]    = msx1_x0 ;
    b0[1]    = msx1_y0 ;
    b0[2]    = msx1    ;

    b1[0]    = msy1_x0 ;
    b1[1]    = msy1_y0 ;
    b1[2]    = msy1    ;


    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }



    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax, bx, et cx                ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);


    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ay, by, et cy                ----- */
    /* ------------------------------------------------------------------ */

    gauss(a1, b1, nb_param, xs1);



    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];

    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];

    goto the_end;

  }


  /* cas polynomiale de degre 2   ...   non lineaire  -- */


  if (type_equa == T_EQUA_2) 
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];

      zx02 = zx0 * zx0;
      zx03 = zx02 * zx0;
      zx04 = zx03 * zx0;
      zy02 = zy0 * zy0;
      zy03 = zy02 * zy0;
      zy04 = zy03 * zy0;

      msx0 += zx0;
      msx02 += zx02;
      msx03 += zx03;
      msx04 += zx04;

      msy0 += zy0;
      msy02 += zy02;
      msy03 += zy03;
      msy04 += zy04;

      msx0_y0 += zx0 * zy0;
      msx02_y0 += zx02 * zy0;
      msx03_y0 += zx03 * zy0;
      msx0_y02 += zx0 * zy02;
      msx0_y03 += zx0 * zy03;
      msx02_y02 += zx02 * zy02;

      msx1 += zx1;
      msx1_x0 += zx1 * zx0;
      msx1_y0 += zx1 * zy0;
      msx1_x02 += zx1 * zx02;
      msx1_y02 += zx1 * zy02;
      msx1_x0_y0 += zx1 * zx0 * zy0;

      msy1 += zy1;
      msy1_x0 += zy1 * zx0;
      msy1_y0 += zy1 * zy0;
      msy1_x02 += zy1 * zx02;
      msy1_y02 += zy1 * zy02;
      msy1_x0_y0 += zy1 * zx0 * zy0;

    }


    a0[0][0] = msx04;
    a0[0][1] = msx02_y02;
    a0[0][2] = msx03_y0;
    a0[0][3] = msx03;
    a0[0][4] = msx02_y0;
    a0[0][5] = msx02;

    a0[1][0] = a0[0][1];
    a0[1][1] = msy04;
    a0[1][2] = msx0_y03;
    a0[1][3] = msx0_y02;
    a0[1][4] = msy03;
    a0[1][5] = msy02;

    a0[2][0] = a0[0][2];
    a0[2][1] = a0[1][2];
    a0[2][2] = msx02_y02;
    a0[2][3] = msx02_y0;
    a0[2][4] = msx0_y02;
    a0[2][5] = msx0_y0;

    a0[3][0] = a0[0][3];
    a0[3][1] = a0[1][3];
    a0[3][2] = a0[2][3];
    a0[3][3] = msx02;
    a0[3][4] = msx0_y0;
    a0[3][5] = msx0;

    a0[4][0] = a0[0][4];
    a0[4][1] = a0[1][4];
    a0[4][2] = a0[2][4];
    a0[4][3] = a0[3][4];
    a0[4][4] = msy02;
    a0[4][5] = msy0;

    a0[5][0] = a0[0][5];
    a0[5][1] = a0[1][5];
    a0[5][2] = a0[2][5];
    a0[5][3] = a0[3][5];
    a0[5][4] = a0[4][5];
    a0[5][5] = nbmax;



    b0[0] = msx1_x02;
    b0[1] = msx1_y02;
    b0[2] = msx1_x0_y0;
    b0[3] = msx1_x0;
    b0[4] = msx1_y0;
    b0[5] = msx1;

    b1[0] = msy1_x02;
    b1[1] = msy1_y02;
    b1[2] = msy1_x0_y0;
    b1[3] = msy1_x0;
    b1[4] = msy1_y0;
    b1[5] = msy1;


    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }

    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ax,bx,cx,dx,ex,fx,gx,hx      ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);



    /* ------------------------------------------------------------------ */
    /* -----           resolution de : ay,by,cy,dy,ey,fy            ----- */
    /* ------------------------------------------------------------------ */

    gauss(a1, b1, nb_param, xs1);



    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];
    *dx = (float)xs[3];
    *ex = (float)xs[4];
    *fx = (float)xs[5];

    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];
    *dy = (float)xs1[3];
    *ey = (float)xs1[4];
    *fy = (float)xs1[5];


    goto the_end;

  }

  /* the_end du polynome de degre 2  */


  /* --cas polynomiale de degre 3   ...   non lineaire  ----- */


  if (type_equa == T_EQUA_3) 
  {
    for (i = 0; i < nbmax; i++) 
    {
      zx0 = x_0[i];
      zx1 = x_1[i];
      zy1 = y_1[i];
      zy0 = y_0[i];

      zx02 = zx0  * zx0;
      zx03 = zx02 * zx0;
      zx04 = zx03 * zx0;
      zx05 = zx04 * zx0;
      zx06 = zx05 * zx0;

      zy02 = zy0  * zy0;
      zy03 = zy02 * zy0;
      zy04 = zy03 * zy0;
      zy05 = zy04 * zy0;
      zy06 = zy05 * zy0;


      msx0  += zx0;
      msx02 += zx02;
      msx03 += zx03;
      msx04 += zx04;
      msx05 += zx05;
      msx06 += zx06;

      msy0  += zy0;
      msy02 += zy02;
      msy03 += zy03;
      msy04 += zy04;
      msy05 += zy05;
      msy06 += zy06;

      msx0_y0    += zx0 * zy0;
      msx0_y02   += zx0 * zy02;
      msx0_y03   += zx0 * zy03;
      msx0_y04   += zx0 * zy04;
      msx0_y05   += zx0 * zy05;

      msx02_y0   += zx02 * zy0;
      msx02_y02  += zx02 * zy02;
      msx02_y03  += zx02 * zy03;
      msx02_y04  += zx02 * zy04;

      msx03_y0   += zx03 * zy0;
      msx03_y02  += zx03 * zy02;
      msx03_y03  += zx03 * zy03;

      msx04_y0   += zx04 * zy0;
      msx04_y02  += zx04 * zy02;

      msx05_y0   += zx05 * zy0;

      msx1 += zx1;
      msx1_x0 += zx1 * zx0;
      msx1_y0 += zx1 * zy0;
      msx1_x02 += zx1 * zx02;
      msx1_y02 += zx1 * zy02;
      msx1_x03 += zx1 * zx03;
      msx1_y03 += zx1 * zy03;

      msx1_x0_y0  += zx1 * zx0  * zy0;
      msx1_x02_y0 += zx1 * zx02 * zy0;
      msx1_x0_y02 += zx1 * zx0  * zy02;

      msy1     += zy1;
      msy1_x0  += zy1 * zx0;
      msy1_y0  += zy1 * zy0;
      msy1_x02 += zy1 * zx02;
      msy1_y02 += zy1 * zy02;
      msy1_x03 += zy1 * zx03;
      msy1_y03 += zy1 * zy03;

      msy1_x0_y0 += zy1 * zx0 * zy0;
      msy1_x02_y0 += zy1 * zx02 * zy0;
      msy1_x0_y02 += zy1 * zx0  * zy02;
    }

    a0[0][0] = msx06;
    a0[0][1] = msx03_y03;
    a0[0][2] = msx05;
    a0[0][3] = msx03_y02;
    a0[0][4] = msx04_y02;
    a0[0][5] = msx05_y0;
    a0[0][6] = msx04_y0;
    a0[0][7] = msx04;
    a0[0][8] = msx03_y0;
    a0[0][9] = msx03;

    a0[1][0] = a0[0][1];
    a0[1][1] = msy06;
    a0[1][2] = msx02_y03;
    a0[1][3] = msy05;
    a0[1][4] = msx0_y05;
    a0[1][5] = msx02_y04;
    a0[1][6] = msx0_y04;
    a0[1][7] = msx0_y03;
    a0[1][8] = msy04;
    a0[1][9] = msy03;

    a0[2][0] = a0[0][2];
    a0[2][1] = a0[1][2];
    a0[2][2] = msx04;
    a0[2][3] = msx02_y02;
    a0[2][4] = msx03_y02;
    a0[2][5] = msx04_y0;
    a0[2][6] = msx03_y0;
    a0[2][7] = msx03;
    a0[2][8] = msx02_y0;
    a0[2][9] = msx02;

    a0[3][0] = a0[0][3];
    a0[3][1] = a0[1][3];
    a0[3][2] = a0[2][3];
    a0[3][3] = msy04;
    a0[3][4] = msx0_y04;
    a0[3][5] = msx02_y03;
    a0[3][6] = msx0_y03;
    a0[3][7] = msx0_y02;
    a0[3][8] = msy03;
    a0[3][9] = msy02;

    a0[4][0] = a0[0][4];
    a0[4][1] = a0[1][4];
    a0[4][2] = a0[2][4];
    a0[4][3] = a0[3][4];
    a0[4][4] = msx02_y04;
    a0[4][5] = msx03_y03;
    a0[4][6] = msx02_y03;
    a0[4][7] = msx02_y02;
    a0[4][8] = msx0_y03;
    a0[4][9] = msx0_y02;

    a0[5][0] = a0[0][5];
    a0[5][1] = a0[1][5];
    a0[5][2] = a0[2][5];
    a0[5][3] = a0[3][5];
    a0[5][4] = a0[4][5];
    a0[5][5] = msx04_y02;
    a0[5][6] = msx03_y02;
    a0[5][7] = msx03_y0;
    a0[5][8] = msx02_y02;
    a0[5][9] = msx02_y0;

    a0[6][0] = a0[0][6];
    a0[6][1] = a0[1][6];
    a0[6][2] = a0[2][6];
    a0[6][3] = a0[3][6];
    a0[6][4] = a0[4][6];
    a0[6][5] = a0[5][6];
    a0[6][6] = msx02_y02;
    a0[6][7] = msx02_y0;
    a0[6][8] = msx0_y02;
    a0[6][9] = msx0_y0;

    a0[7][0] = a0[0][7];
    a0[7][1] = a0[1][7];
    a0[7][2] = a0[2][7];
    a0[7][3] = a0[3][7];
    a0[7][4] = a0[4][7];
    a0[7][5] = a0[5][7];
    a0[7][6] = a0[6][7];
    a0[7][7] = msx02;
    a0[7][8] = msx0_y0;
    a0[7][9] = msx0;

    a0[8][0] = a0[0][8];
    a0[8][1] = a0[1][8];
    a0[8][2] = a0[2][8];
    a0[8][3] = a0[3][8];
    a0[8][4] = a0[4][8];
    a0[8][5] = a0[5][8];
    a0[8][6] = a0[6][8];
    a0[8][7] = a0[7][8];
    a0[8][8] = msy02;
    a0[8][9] = msy0;

    a0[9][0] = a0[0][9];
    a0[9][1] = a0[1][9];
    a0[9][2] = a0[2][9];
    a0[9][3] = a0[3][9];
    a0[9][4] = a0[4][9];
    a0[9][5] = a0[5][9];
    a0[9][6] = a0[6][9];
    a0[9][7] = a0[7][9];
    a0[9][8] = a0[8][9];
    a0[9][9] = nbmax;



    b0[0] = msx1_x03;
    b0[1] = msx1_y03;
    b0[2] = msx1_x02;
    b0[3] = msx1_y02;
    b0[4] = msx1_x0_y02;
    b0[5] = msx1_x02_y0;
    b0[6] = msx1_x0_y0;
    b0[7] = msx1_x0;
    b0[8] = msx1_y0;
    b0[9] = msx1;

    b1[0] = msy1_x03;
    b1[1] = msy1_y03;
    b1[2] = msy1_x02;
    b1[3] = msy1_y02;
    b1[4] = msy1_x0_y02;
    b1[5] = msy1_x02_y0;
    b1[6] = msy1_x0_y0;
    b1[7] = msy1_x0;
    b1[8] = msy1_y0;
    b1[9] = msy1;


    for (i = 0; i< nb_param;i++)
        {
         for (j = 0; j < nb_param; j++)
             {
              a1[i][j] = a0[i][j] ;
             }
         }

    /* ------------------------------------------------------------------ */
    /* -----      resolution de : ax,bx,cx,dx,ex,fx,gx,hx,ix,jx     ----- */
    /* ------------------------------------------------------------------ */

    gauss(a0, b0, nb_param, xs);


    /* ------------------------------------------------------------------ */
    /* -----    resolution de : ay,by,cy,dy,ey,fy,gy,hy,iy,jy       ----- */
    /* ------------------------------------------------------------------ */
  
    gauss(a1, b1, nb_param, xs1);

    *ax = (float)xs[0];
    *bx = (float)xs[1];
    *cx = (float)xs[2];
    *dx = (float)xs[3];
    *ex = (float)xs[4];
    *fx = (float)xs[5];
    *gx = (float)xs[6];
    *hx = (float)xs[7];
    *ix = (float)xs[8];
    *jx = (float)xs[9];


    *ay = (float)xs1[0];
    *by = (float)xs1[1];
    *cy = (float)xs1[2];
    *dy = (float)xs1[3];
    *ey = (float)xs1[4];
    *fy = (float)xs1[5];
    *gy = (float)xs1[6];
    *hy = (float)xs1[7];
    *iy = (float)xs1[8];
    *jy = (float)xs1[9];

    goto the_end;

  }
  /* the_end du polynome de degre 3  */


  the_end: ;

  /* ------ liberation de l'espace memoire ------ */

  delete [] y11 ;
  delete [] x11 ;

  delete [] y10 ;
  delete [] y00 ;
  delete [] x10 ;
  delete [] x00 ;

  delete [] xs1 ;
  delete [] xs ;
  delete [] b1 ;
  delete [] b0 ;

      /* ------ liberation de l'espace memoire 2D ------ */

 for (i = 0; i < nb_param; i++)
   { 
     delete [] a0[i] ;
   }
 delete [] a0 ;

 for (i = 0; i < nb_param; i++)
   { 
     delete [] a1[i] ;
   }
 delete [] a1 ;

  /* ------ liberation de l'espace memoire ------ */
}
/* End. */

/**************************************************************************/

void statis (float* tab_in, int nb_pts, char* name_mouch,  char*name_tab)
{
  long   j ;
  double  som, som1 ;
  double  moy, min, max, ecart, ecarq ;
  FILE    *file_mouch ;
  
  file_mouch = fopen(name_mouch,"a") ;
  if (file_mouch == NULL)
   {
    cerr << "Error: cannot open output file : " << name_mouch << endl;
    exit(-1) ;
   }


  moy   = 0.0;
  min   = 1e10;
  max   = -1e10;
  ecart = 0.0;
  ecarq = 0.0;

  som  = 0.0;
  som1 = 0.0;

    for (j = 0; j < nb_pts; j++)
    {

      if (tab_in[j] >= max) max = tab_in[j];
      
      if (tab_in[j] <= min)  min = tab_in[j];
      
      som += tab_in[j] ;
      som1 += tab_in[j] * tab_in[j];
    }
  if (nb_pts > 0)
  {
  moy   = (float)(som / nb_pts) ;
  ecarq = (float)((nb_pts * (moy) * (moy)) / (1.0 - nb_pts) + (som1 / (nb_pts - 1.0))) ;
  if (ecarq > 0) ecart = (float)sqrt((double)(ecarq));
  else ecarq = 0.;
 }

  /* ------------------------------------------------------------ */
  /* --------          sequence d'ecriture               -------- */
  /* ------------------------------------------------------------ */

  fprintf(file_mouch, "============  Statistiques ============ \n");
  fprintf(file_mouch, "   \n");
  fprintf(file_mouch, "Tableau traiter           :  %s\n", name_tab);
  fprintf(file_mouch, "Nombre de Points total    :  %d\n", nb_pts);
  fprintf(file_mouch, "  \n");
  fprintf(file_mouch, "Minimum     :  %f\n", min);
  fprintf(file_mouch, "Maximum     :  %f\n", max);
  fprintf(file_mouch, "Moyenne     :  %f\n", moy);
  fprintf(file_mouch, "Variance    :  %f\n", ecarq);
  fprintf(file_mouch, "Ecart Type  :  %f\n", ecart);
  fprintf(file_mouch, "  \n");
  fprintf(file_mouch, "  \n");

  fclose(file_mouch);


  }
/* End. */

/**************************************************************************/

void gauss(double **mat, double *vec1, int dimension, double *vec2)
{
  int     i, j, k, k1, n, sup1=0;
  double  sup, a, a1, a2, somme, epsilon;
  n = dimension;          /* dimension de la matrice  */
  epsilon = 0.000000001;


  /* ------------------------------- */
  /* ----   triangularisation   ---- */
  /* ------------------------------- */

  for (k = 1; k <= n-1; k++) 
 {
    /* --------------------------------------------------------- */
    /* ----  recherche du pivot maximal colonne par colonne ---- */
    /* --------------------------------------------------------- */

    k1 = k;
    sup = 0.0;
    for (i = k1; i <= n-1; i++)
    {
      if (fabs(mat[i-1][k-1]) >= sup) 
      {
        sup = fabs(mat[i-1][k-1]) ;
	sup1 = i;
      }
    }
    if (sup <= epsilon) 
    {
      printf("Warning: cannot invert the matrix ... \n");
      goto fin;
    }

    /* ------------------------------- */
    /* ----- fin de la recherche  ---- */
    /* ------------------------------- */




    /* ------------------------------- */
    /* -----    substitution    ------ */
    /* ------------------------------- */

    /* ------------------ */
    /* ---   lignes  ---- */
    /* ------------------ */

    for (j = 1; j <= n ; j++) 
    {
      a = mat[k1-1][j-1];
      mat[k1-1][j-1] = mat[sup1-1][j-1];
      mat[sup1-1][j-1] = a;
    }

    /* ------------------- */
    /* -- second membre -- */
    /* ------------------- */


    a = vec1[k1-1];
    vec1[k1-1] = vec1[sup1-1];
    vec1[sup1-1] = a;

    /* ------------------------------- */
    /* ---- fin de la substitution ----*/
    /* ------------------------------- */

    /* ------------------------- */
    /* --- triangularisation --- */
    /* ------------------------- */

    for (i = k1+1 ; i <= n; i++) 
    {
      a = (double)mat[i-1][k-1] / mat[k1-1][k1-1];
      mat[i-1][k1-1] = 0;
      for (j = k1+1 ; j <= n; j++)
      {
	a1 = a * mat[k1-1][j-1];
	mat[i-1][j-1] -= a1;
      }
      a2 = a * vec1[k1-1];
      vec1[i-1] -= a2;
    }


    /* ----------------------------------- */
    /* --- fin de la triangularisation --- */
    /* ----------------------------------- */
  }



  if (fabs(mat[n-1][n-1]) <= epsilon) 
  {
    printf("Resolution impossible\n");
    goto fin;
  }

  /* --------------------------------------------- */
  /* ---  resolution du systeme triangulaire  ---- */
  /* --------------------------------------------- */

  vec2[n-1] = (vec1[n-1] / mat[n-1][n-1]) ;
  if (fabs(vec2[n-1]) < epsilon)
    vec2[n-1] = 0;


  for (i = n - 1; i >= 1; i--) 
  {
    somme = 0.0;
    for (j = i + 1; j <= n; j++) 
    {
      somme += mat[i-1][j-1] / mat[i-1][i-1] * vec2[j-1];
    }
    vec2[i-1] = vec1[i-1] / mat[i-1][i-1] - somme;
    if (fabs(vec2[i-1]) < epsilon) 
    {
      vec2[i-1] = 0;
    }
  }

fin: ;

}

/* End. */


/**************************************************************************/

void copie_pts_fichier(char *nom_entree, 
                       float *x0, float *y0, float *x1, float *y1, 
                       int nmax) 
{
  FILE   *file_in ;
  int    i ;

  file_in = fopen(nom_entree,"w+") ;
  if (file_in == NULL)
   {
    printf("Error in opening file %s \n", nom_entree) ;
    exit(-1) ;
   }

  fprintf(file_in, "%d\n", nmax);
  for (i = 0 ; i < nmax ; i++)
    fprintf(file_in, "%4.2f %4.2f      %4.2f %4.2f\n", 
                                    x0[i], y0[i], x1[i], y1[i]);
 
  fclose(file_in) ;
}
/* end */

/**************************************************************************/
