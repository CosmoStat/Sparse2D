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
**    File:  MR_Registration.cc
**
*******************************************************************************
**
**    DESCRIPTION  routines for image regustration
**    -----------  
**                 
**
**
*******************************************************************************
**
** detection (Ifloat &Im_seuil, float **x_e, float **y_e,float **f_e, int *nmax)
** 
** Ifloat &Im_seuil : image seuillee
** float *x_e, *y_e : vecteur contenant les coordonnees en
**                    x et y des maximas detectes
** float *f_e : vecteur contenant la valeur pixel du maxima detecte
** int  *nmax : nbre de maxima detecte
** 
** Detection de maxima dans l'image Im_seuil. 
** Les coordonnes des maxima ainsi que leur valeur son conserve dans les
** vecteurs x_e, y_e et nmax.
** 
******************************************************************************
**
** recalage (parametre  *par, int type_equa, int type_rec,float a,
**           Ifloat &Im_in, Ifloat &Im_out)
**
** Recalage geometique de l'image en coordonnees pixel a l'aide
** d'un modele polynomiale de degree 1, 2 ou 3
**
** float  a : parametre a pour la convolution cubique
** int type_equa, type_rec : type d'equation pour le recalage et
**                           type de recalage : meme resolution ou
**                           resolution differente
** parametre  *par : parametre du modele polynomiale de correction geometrique
** Ifloat &Im_in, &Im_out: image d'entree et de sortie
**
** 
******************************************************************************
** 
** mise_en_correspondance_coor_reel(int resolution, int res_min, 
**                                    int nature_rec,
**                                    int depart_br_x, int depart_br_y, 
**                                    float pas_br_x, float pas_br_y,
**                                    float *x0_e, float *y0_e, float *x1_e, 
**                                    float *y1_e, float *xr_e, float *yr_e,
**                                    float rayon, int nmax0, int nmax1, 
**                                    int type_equa, 
**                                    float **x0_sf, float **y0_sf, 
**                                    float **x1_sf, float **y1_sf,
**                                    int *nmax,
**                                    parametre *par, parametre *par_cr ) 
** 
** 
** parametre  *par                     : parametre du polynome en sortie
**            *par_cr                  : parametre du polynome en sortie
**                                       en coordonnee reelle
** int        nmax0, nmax1     : nbre de points pour le premier systeme de 
**                               vecteur x0_e et y0_e et pour le deuxieme
**                               systeme de vecteur x1_e et y1_e
** int        nature_rec       : nature du recalage : recalage sur image de
**                               meme resolution ou de resolution differente
** int        *nmax,type_equa  : nombre de maxima apparie et type d'equation
**                               pour l'appriement
** int        resolution, res_min : resolution en cours et resolution min
** int        depart_br_x, depart_br_y : coordonnes de depart en basse 
**                                       resolution en X et Y
** float      rayon                    : rayon de recherche pour le matching
** float      pas_br_x, pas_br_y       : pas d'echantillonnage en basse res
**                                       en X et Y
** float      *x0_e, *y0_e,    : respectivement vecteur de valeurs en X et Y
**            *x1_e, *y1_e,      des points a apparier : x0_e et y0_e coord
**            *xr_e, *yr_e;      des points de l'image de reference, x1_e et
**                               y1_e coord des points de l'image de travaille,
**                               xr_e et yr_e coordonnees des points de l'image
**                               de reference prealablement transformee a l'aide
**                               du modele polynomiale par transpo_coor_max.c
** float    **x0_sf, **y0_sf,    : vecteur de point apparies en X et Y 
**          **x1_sf, **y1_sf       respectivement pour l'image de referene et
**                               celle de travail
**            
** Mise en correspondance en coordonnees reelles : Appariement de paires de 
** points entre les amers de l'image de reference et celle de travail a
** l'aide d'un modele polynomiale avec calcul du residu et du root mean square
** distance error, ainsi que des parametres du polynome de rectification
** geometrique. Le polynome est d'ordre1, 2 ou 3.
** 
** 
** PROCEDURE UTILISEE : moindre_carre
**                      gauss (utilise par la procedure moindre carre)
**                      copie_pts_fichier
** 
** 
******************************************************************************
**
** void recalage_coor_reelle(parametre *par,
**                          int type_equa, int type_rec, float a,
**                          Ifloat &Im_in, Ifloat &Im_out,
**                          int depart_x, int depart_y, 
**                          float pas_x, float pas_y)
** 
** Ifloat &Im_in, Ifloat &Im_out: image d'entree et de sortie
** float     pas_x, pas_y, a : pas d'echantillonnage en x et y et valeur
**                             du parametre a pour la convolution cubique
**           type_equa, type_rec : type d'equation pour le recalage et
**                                 type de recalage : meme resolution ou
**                                 resolution differente
**           depart_x, depart_y  : coordonnees de depart en x et y
** parametre *par                 : parametre du modele polynomiale de 
**                                  correction geometrique
** 
** Recallage geometique de l'image en coordonnees reelles ou pixel a l'aide
** d'un modele polynomiale de degree 1, 2 ou 3

******************************************************************************/

// static char sccsid[] = "@(#)MR_Registration.cc 3.1 96/05/02 CEA 1996 @(#)";

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_IO.h"
#include "MR_Fusion.h"

// #define MOUCHARD 1
/**************************************************************************/

void recalage_coor_reelle(parametre *par,
                          int type_equa, int type_rec, float a,
                          Ifloat &Im_in, Ifloat &Im_out,
                          int depart_x, int depart_y, float pas_x, float pas_y)

{

  /* ----------------------------------------------- */
  /* ------   mapping image output to input   ------ */
  /* ------   for float images only           ------ */
  /* ------   real coordinates system only    ------ */
  /* ----------------------------------------------- */


  float           x, y,
                  x_reel=0., y_reel=0.,
                  ligne_reel, col_reel,
                  ligne2_reel, col2_reel,
                  ligne3_reel, col3_reel ;
 float            depart_x_cr, depart_y_cr ;

  float	          t1, t2, t3, t4, a1, a2, a3, a4, b1, b2, b3, b4, val ;
  float	  	  delta_x, delta_x_2, delta_x_3, delta_y, delta_y_2, delta_y_3 ;
  int             nb_lig, nb_col ;

  long            compteur, nbl0;
  long            col, ligne,
                  iix, iiy, 
		  iix1, iiy1 ,
		  iix2, iiy2 ,
		  iix0, iiy0 ;

  float    val1 ;


  /*######################################################################## */
  /*#########                Traitement                             ######## */
  /*######################################################################## */

  nb_lig = Im_in.nl() ;
  nb_col = Im_in.nc() ;

  /* ----------------------------------------------------------------------- */
  /* ---    initialisation de l'image recale a zero                      --- */
  /* ---    l'image recale aura la meme dimension que                    --- */
  /* ---    l'image de reference                                         --- */
  /* ----------------------------------------------------------------------- */

  nbl0 = nb_col;


//  for (i = 0; i < nb_lig; i++)
 //   {
  //     for (j = 0; j < nb_col; j++)
   //      {
    //         Im_out[i][j] = 0 ;
    //     }
   // }




 
/* transformation des coordonnees depart en coordonnees depart reelles */

   depart_x_cr = 1. + ((float)depart_x * pas_x) ;
   depart_y_cr = 1. + ((float)depart_y * pas_y) ;




  for (ligne = 0; ligne < nb_lig; ligne++) 
  {
    compteur = (ligne) * nbl0;
    /* transformation des coordonnees pixels lignes en coordonnees reelles */
       ligne_reel = depart_y_cr + ((float)ligne * pas_y) ;
 


    for (col = 0; col < nb_col; col++) 
     {
      /* transformation des coordonnees pixels colonnes en coordonnees reelles */
         col_reel = depart_x_cr + ((float)col * pas_x) ;
      

      if (type_equa == T_EQUA_0) 
     {
	x_reel = par->ax * col_reel - par->bx * ligne_reel + par->cx ;
	y_reel = par->bx * col_reel + par->ax * ligne_reel + par->cy ;
      }

      if (type_equa == T_EQUA_1) 
      {
	x_reel = par->ax * col_reel + par->bx * ligne_reel + par->cx ;
	y_reel = par->ay * col_reel + par->by * ligne_reel + par->cy ;
      }


      if (type_equa == T_EQUA_2) 
      {
	col2_reel = col_reel * col_reel;
	ligne2_reel = ligne_reel * ligne_reel;

	x_reel = par->ax * col2_reel + 
		 par->bx * ligne2_reel +
		 par->cx * ligne_reel * col_reel +
                 par->dx * col_reel  +
		 par->ex * ligne_reel  +
		 par->fx ;
	y_reel = par->ay * col2_reel +
		 par->by * ligne2_reel +
		 par->cy * ligne_reel * col_reel + 
                 par->dy * col_reel  +
		 par->ey * ligne_reel  +
		 par->fy ;
      }

  if (type_equa == T_EQUA_3) 
  {
      col2_reel = col_reel  * col_reel;
      col3_reel = col2_reel * col_reel;

      ligne2_reel = ligne_reel  * ligne_reel ;
      ligne3_reel = ligne2_reel * ligne_reel ;

      x_reel = par->ax * col3_reel              +
	       par->bx * ligne3_reel            +
               par->cx * col2_reel              + 
	       par->dx * ligne2_reel            +
               par->ex * col_reel * ligne2_reel + 
	       par->fx * col2_reel * ligne_reel +
               par->gx * col_reel * ligne_reel  +
	       par->hx * col_reel               +
               par->ix * ligne_reel             +
    	       par->jx ;
  

      y_reel = par->ay * col3_reel              + 
	       par->by * ligne3_reel            +
               par->cy * col2_reel              + 
	       par->dy * ligne2_reel            + 
               par->ey * col_reel * ligne2_reel + 
	       par->fy * col2_reel * ligne_reel +
               par->gy * col_reel * ligne_reel  + 
	       par->hy * col_reel               +
               par->iy * ligne_reel             + 
	       par->jy ;

  }




      /* transformation des coordonnees reelles au coordonnees pixels */
         x = (x_reel - depart_x_cr) / pas_x ;        
         y = (y_reel - depart_y_cr) / pas_y ;        



           /* calcul de la valeur interpole a assigne au pixel de l'image a recale */

           iix = (int) floor(x) ;
           iiy = (int) floor(y) ;


           if (iix >= nb_col) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iix < 0) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iiy >= nb_lig) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iiy < 0) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

 
           delta_x = x - iix ;
           delta_y = y - iiy ;


           /* cas de l'interpolation du plus proche voisin */

           if (type_rec == T_REC_CLOSER_PIXEL) 
            {
	      Im_out(ligne,col) = Im_in(iiy,iix) ;
            }



          /* cas de l'interpolation bilineaire */
          if (type_rec == T_REC_BILINEAR_INTERPOL)
           {
             iiy1 = iiy + 1 ;
             iix1 = iix + 1 ;
             if (iiy1 >= nb_lig)
              {
                iiy1 = 2 * nb_lig - iiy1 - 1 ;
              }
             if (iix1 >= nb_col)
	      {
                iix1 = 2 * nb_col - iix1 - 1;
              }


	     t1 = ((1.0 - delta_x) * (double)(Im_in(iiy,iix))  ) + 
		  ((double)((Im_in(iiy,iix1))) * delta_x) ;
	     t2 = ((double)(Im_in(iiy1,iix)) * (1.0 - delta_x)) + 
		  ((double)(Im_in(iiy1,iix1)) * delta_x) ;
             val = (t1 * (1. - delta_y)) + (t2 * delta_y) ;
             /*
             if (val >= 255)
                val = 255 ;
             if (val <= 0)
                val = 0 ;
             val1 = (irint(val)) ;
             */

             val1 = val ;
 
             Im_out(ligne,col) = val1 ;
   	   }



         /* cas de la convolution cubique avec parametre a */

         if (type_rec == T_REC_CUBIC_INTERPOL)
          {
            delta_x_2 = delta_x * delta_x ;
            delta_x_3 = delta_x_2 * delta_x ;
            delta_y_2 = delta_y * delta_y ;
            delta_y_3 = delta_y_2 * delta_y ;

	    a1 = (a-2)*(delta_x_3 - 2*delta_x_2 + delta_x) ;
            a2 = a * delta_x_3 - (1+a)*delta_x_2 + 1 ;
            a3 = -a*delta_x_3 + (2*a-1)*delta_x_2 + (2-a)*delta_x ;
            a4 = (a-2) * (delta_x_2 - delta_x_3) ;

	    b1 = (a-2)*(delta_y_3 - 2*delta_y_2 + delta_y) ;
            b2 = a * delta_y_3 - (1+a)*delta_y_2 + 1 ;
            b3 = -a*delta_y_3  + (2*a - 1)*delta_y_2 + (2-a)*delta_y ;
            b4 = (a-2) * (delta_y_2 - delta_y_3) ;

	    iiy0 = iiy - 1 ;
	    iiy1 = iiy + 1 ;
	    iiy2 = iiy + 2 ;

            iix0 = iix - 1 ;
            iix1 = iix + 1 ;
            iix2 = iix + 2 ;

            if (iix0 < 0)
               iix0 = iix0 - 2 * iix0 - 1 ;
            if (iix1 >= nb_col)
               iix1 = 2 * nb_col - iix1 -1 ;
            if (iix2 >= nb_col)
               iix2 = 2 * nb_col - iix2 -1 ;

            if (iiy0 < 0)
               iiy0 = iiy0 - 2 * iiy0 - 1 ;
            if (iiy1 >= nb_lig)
               iiy1 = 2 * nb_lig - iiy1 -1 ;
            if (iiy2 >= nb_lig)
               iiy2 = 2 * nb_lig - iiy2 -1 ;

            t1 = (double)(Im_in(iiy0,iix0))*a1 + 
		 (double)(Im_in(iiy0,iix))*a2  + 
		 (double)(Im_in(iiy0,iix1))*a3 + 
		 (double)(Im_in(iiy0,iix2))*a4 ;
            t2 = (double)(Im_in(iiy,iix0))*a1  + 
		 (double)(Im_in(iiy,iix))*a2   + 
		 (double)(Im_in(iiy,iix1))*a3  + 
		 (double)(Im_in(iiy,iix2))*a4 ;
            t3 = (double)(Im_in(iiy1,iix0))*a1 + 
		 (double)(Im_in(iiy1,iix))*a2  + 
		 (double)(Im_in(iiy1,iix1))*a3 + 
		 (double)(Im_in(iiy1,iix2))*a4 ;
            t4 = (double)(Im_in(iiy2,iix0))*a1 + 
		 (double)(Im_in(iiy2,iix))*a2  + 
		 (double)(Im_in(iiy2,iix1))*a3 +
		 (double)(Im_in(iiy2,iix2))*a4 ;


	    val = (t1*b1 + t2*b2 + t3*b3 + t4*b4) ;
            /*
            if (val >= 255)
               val = 255 ;
            if (val <= 0)
               val = 0 ;
            val1 = (irint(val)) ;
            */
            val1 = val ;
	    Im_out(ligne,col) = val1 ;

           }



          suite: ;
  
      } /* end col */
   } /* end ligne */
} /* End. */

/**************************************************************************/

void mise_en_correspondance_coor_reel(int resolution, int res_min, int nature_rec,
                                      int depart_br_x, int depart_br_y, 
				      float pas_br_x, float pas_br_y,
                          	      float *x0_e, float *y0_e, float *x1_e, 
                                      float *y1_e, float *xr_e, float *yr_e,
                          	      float rayon, int nmax0, int nmax1, 
                          	      int type_equa, 
                           	      float **x0_sf, float **y0_sf, 
                                      float **x1_sf, float **y1_sf,
                          	      int *nmax,
				      parametre *par, parametre *par_cr ) 

{

  float	    *x0_s, *y0_s, *x1_s, *y1_s;
  float      *x00_s, *y00_s, *x10_s, *y10_s, *x00r_s, *y00r_s ;
  float      a00, a11, b00, b11 ;
  int        l, i, j;
  int        pt_mult , max ;
  int        nb_mauvais, nb_mauvais_redondant;
  int        iteration, dim ;
  float	     seuil_residu ;
 
  float      *x01_r , *y01_r;
  float      *x0_cr_s, *y0_cr_s, *x1_cr_s, *y1_cr_s ;
  float      a1x, b1x, c1x, d1x, e1x, f1x, g1x, h1x, i1x, j1x,
             a1y, b1y, c1y, d1y, e1y, f1y, g1y, h1y, i1y, j1y ;
  float      a2x, b2x, c2x, d2x, e2x, f2x, g2x, h2x, i2x, j2x,
             a2y, b2y, c2y, d2y, e2y, f2y, g2y, h2y, i2y, j2y ;
  float      *residu , *dif_x, *dif_y ;
  float      somme_residu, residu_total, somme_dif2, 
             somme_dif2_cr, RMSE, RMSE_cr ;
  float      res0, res1 ;

  float       col2, lig2, col3, lig3 ;
  float      a, b, c, d, e, a0, b0, z1, z2, ro_2;

#ifdef MOUCHARD
  FILE *file_residu ;
  char	     name_tab[256] ;

  char *mouch_residu = "mouchard.txt";
#endif


  if (nmax0 > nmax1) max = nmax0 ;
  else max = nmax1 ;

  dim = 5 * max ;

   // allocation dynamique des variables   

    x00_s  = new float[dim] ;
    y00_s  = new float[dim] ;
    x10_s  = new float[dim] ;
    y10_s  = new float[dim] ;
    x00r_s = new float[dim] ;
    y00r_s = new float[dim] ;

 
  l = -1 ;
  e = rayon;
  ro_2 = e * e;


  // le seuil est fixe a la demi taille de l'ondelette */
  // et a pour limite inferieure 1           

  if (rayon >= 2) seuil_residu = (rayon+1) / 2 ;  
  else seuil_residu = 1 ;


/* ----------------------------------------------------------------------- */
/* ---- (a0,b0) est le centre du cercle de recherche, ro est le rayon ---- */
/* ---- Le centre du cercle sera le point prealablement rectifie par  ---- */
/* ---- rapport a l'equation de l'etape precedente                    ---- */
/* ----------------------------------------------------------------------- */


  for (i = 0; i < nmax0; i++)
    {
      b  = y0_e[i];
      a  = x0_e[i];
      b0 = yr_e[i];
      a0 = xr_e[i];
      pt_mult = 0;

      for (j = 0; j < nmax1; j++) 
        {
          d = y1_e[j];
          c = x1_e[j];
          z1 = (c - a0) * (c - a0);
          z2 = (d - b0) * (d - b0);

          if (z1 + z2 <= ro_2) 
           {
	     pt_mult++;
	     l++;

	     x00_s[l] = a;
	     y00_s[l] = b;
	     x00r_s[l] = a0;
	     y00r_s[l] = b0;
	     x10_s[l] = x1_e[j];
	     y10_s[l] = y1_e[j] ;

             if ((l+1) == dim)
              {
                dim += (5 * max) ;
                x00_s  = (float *) realloc( x00_s,  dim * sizeof(float)) ;
                y00_s  = (float *) realloc( y00_s,  dim * sizeof(float)) ;
                x10_s  = (float *) realloc( x10_s,  dim * sizeof(float)) ;
                y10_s  = (float *) realloc( y10_s,  dim * sizeof(float)) ;
                x00r_s = (float *) realloc( x00r_s, dim * sizeof(float)) ;
                y00r_s = (float *) realloc( y00r_s, dim * sizeof(float)) ;
              }
           }
        }
    }


  *nmax = l + 1 ;
  l = 0 ;
  if (*nmax <= 2)
   {
     cerr << "WARNING: not enough detected maxima at scale " << resolution+1 << " : " << *nmax << endl;
     return;
   }

  /* -------------------------------------------------- */  
  /* ----  traitement des valeurs aux bords      ------ */
  /* ----  et elimination de ces valeurs         ------ */
  /* -------------------------------------------------- */  

  /* ----  mise a zero des points se trouvant sur les bords ---- */

  for (i = 0 ; i < *nmax; i++)
    {
      if (x00_s[i] == 0 || y00_s[i] == 0 || x10_s[i] == 0 || y10_s[i] == 0 
         ||  x00r_s[i] == 0 ||  y00r_s[i] == 0 ) 
       {  
         x00_s[i]  = 0. ;
         y00_s[i]  = 0. ;
         x10_s[i]  = 0. ;
         y10_s[i]  = 0. ;
         x00r_s[i] = 0. ;
         y00r_s[i] = 0. ;
       }
     }
  
  /* ----  elimination des points aux bords ---- */
  l = 0 ;
  for (i = 0; i < *nmax; i++) 
  {
      if (x00_s[i] != 0 && y00_s[i] !=0)
      {
         x00_s[l]  = x00_s[i]  ;
         y00_s[l]  = y00_s[i]  ;
         x10_s[l]  = x10_s[i]  ;
         y10_s[l]  = y10_s[i]  ;
         x00r_s[l] = x00r_s[i] ;
         y00r_s[l] = y00r_s[i] ;
         l++;
      }
  }
  *nmax = l ;

 
  /* -------------------------------------------------- */   
  /* ----   ecriture des valeurs en sortie apres  ----- */
  /* ----   elimination des points defectueux     ----- */
  /* ----   du pretraitement                      ----- */
  /* -------------------------------------------------- */   


   /* -----   allocation dynamique de x0_s, y0_s, x1_s, y1_s ------ */

   x0_s = new float[*nmax] ;
   y0_s = new float[*nmax] ;
   x1_s = new float[*nmax] ;
   y1_s = new float[*nmax] ;
 
   /* -----   allocation dynamique de x0_s, y0_s, x1_s, y1_s ------ */
 
  for (i = 0; i < *nmax; i++)
    {
      x0_s[i] = x00_s[i];
      y0_s[i] = y00_s[i];
      x1_s[i] = x10_s[i];
      y1_s[i] = y10_s[i];
    }

  delete [] x00_s ;
  delete [] y00_s ;
  delete [] x10_s ;
  delete [] y10_s ;
  delete [] x00r_s ;
  delete [] y00r_s ;


  if (*nmax <= 2) 
   {
    cerr << " WARNING: maximas  at border have been removed " << endl;
    cerr << "          Not enough detected maxima at scale " << resolution+1 << ": " << *nmax << endl;
    return;
   }
 

iteration = 0 ;

   /* ------------    allocation dynamique des variables  ------------ */

    x01_r  = new float[*nmax] ;
    y01_r  = new float[*nmax] ;
    residu = new float[*nmax] ;
    dif_x  = new float[*nmax] ;
    dif_y  = new float[*nmax] ;

   /* ------------    fin de l'allocation dynamique      ------------ */

recommencer:
  if (*nmax <= 2)
   {
     cerr << "WARNING: not enough detected maxima at scale " << resolution+1 << ": " << *nmax << endl;
     return;
   }
   
iteration +=  1 ;

  /* -------------------------------------------------- */   
  /* ----            calcul des residus           ----- */
  /* -------------------------------------------------- */   

  moindre_carre(*nmax, x0_s, y0_s, x1_s, y1_s,
		type_equa,
                &a1x, &b1x, &c1x, &d1x, &e1x, &f1x, &g1x, &h1x, &i1x, &j1x,
                &a1y, &b1y, &c1y, &d1y,	&e1y, &f1y, &g1y, &h1y, &i1y, &j1y);

   par->ax = a1x ;
   par->bx = b1x ;
   par->cx = c1x ;
   par->dx = d1x ;
   par->ex = e1x ;
   par->fx = f1x ;
   par->gx = g1x ;
   par->hx = h1x ;
   par->ix = i1x ;
   par->jx = j1x ;


   par->ay = a1y ;
   par->by = b1y ;
   par->cy = c1y ;
   par->dy = d1y ;
   par->ey = e1y ;
   par->fy = f1y ;
   par->gy = g1y ;
   par->hy = h1y ;
   par->iy = i1y ;
   par->jy = j1y ;

  nb_mauvais   = 0 ;  
  somme_residu = 0 ;
  somme_dif2   = 0 ;

  for (i = 0; i < *nmax; i++)
  {
    switch(type_equa)
    {
       case T_EQUA_0:
          x01_r[i] = (float)((par->ax) * x0_s[i] - 
                         (par->bx) * y0_s[i] + (par->cx)) ;
          y01_r[i] = (float)((par->bx) * x0_s[i] + 
                         (par->ax) * y0_s[i] + (par->cy)) ;
          break;
       case T_EQUA_1:
           x01_r[i] = (float)((par->ax) * x0_s[i] + 
                                (par->bx) * y0_s[i] + (par->cx)) ;
           y01_r[i] = (float)((par->ay) * x0_s[i] + 
                                (par->by) * y0_s[i] + (par->cy)) ;
          break;
       case T_EQUA_2:    
          col2 = x0_s[i] * x0_s[i];
          lig2 = y0_s[i] * y0_s[i];
          x01_r[i] = (float)((par->ax) * col2 + (par->bx) * lig2 + 
                 (par->cx) * y0_s[i] * x0_s[i] +
                 (par->dx) * x0_s[i] + (par->ex) * y0_s[i] + (par->fx)) ;
          y01_r[i] = (float)((par->ay) * col2 + (par->by) * lig2 + 
                 (par->cy) * y0_s[i] * x0_s[i] +
	         (par->dy) * x0_s[i] + (par->ey) * y0_s[i] + (par->fy)) ;
          break;
       case T_EQUA_3:     
          col2 = x0_s[i] * x0_s[i];
          lig2 = y0_s[i] * y0_s[i];
          col3 = col2 * x0_s[i];
          lig3 = lig2 * y0_s[i];
          x01_r[i] = (par->ax) * col3  +  (par->bx) * lig3  +  (par->cx) * col2  +
                 (par->dx) * lig2  +  (par->ex) * x0_s[i] * lig2  +
                 (par->fx) * y0_s[i] * col2  +  (par->gx) * y0_s[i] * x0_s[i]  +
                 (par->hx) * x0_s[i] + (par->ix) * y0_s[i] + (par->jx) ;
          y01_r[i] = (par->ay) * col3  +  (par->by) * lig3  +  (par->cy) * col2  +
                 (par->dy) * lig2  +  (par->ey) * x0_s[i] * lig2  +
                 (par->fy) * y0_s[i] * col2  +  (par->gy) * y0_s[i] * x0_s[i]  +
                 (par->hy) * x0_s[i] + (par->iy) * y0_s[i] + (par->jy) ;
          break;
       default: break;
    }
    
     dif_x[i] = (float)(x1_s[i]) - x01_r[i];
     dif_y[i] = (float)(y1_s[i]) - y01_r[i];
     residu[i] = (float)(sqrt((double)(dif_x[i] * dif_x[i] + 
                                        dif_y[i] * dif_y[i])));
     if (residu[i] > seuil_residu)
      {
	x0_s[i] = 0;
	y0_s[i] = 0;
        x1_s[i] = 0;
	y1_s[i] = 0;
        nb_mauvais++;
      }
      somme_residu += residu[i] ;
      somme_dif2   += (float)( (dif_x[i] * dif_x[i]) +  (dif_y[i] * dif_y[i]) ) ;
   } // end for
  
   if (*nmax > 0) residu_total = somme_residu / (float)*nmax ;
   else residu_total = 0;
  
  if (*nmax > 0) RMSE = sqrt(somme_dif2 / (float)(*nmax) ) ;
  else RMSE = 0;
  
  
#ifdef MOUCHARD

  file_residu = fopen(mouch_residu,"a") ;
  fprintf(file_residu,"Resolution    : %d\n",resolution) ;
  fprintf(file_residu,"----------------------------\n") ;
  fprintf(file_residu,"nb_points     =  %d\n",(*nmax)) ;
  fprintf(file_residu,"nb_mauvais    =  %d\n",nb_mauvais) ;
  fprintf(file_residu,"iteration     =  %d\n",iteration) ;
  fprintf(file_residu,"root mean square error     =  %f\n",RMSE) ;
  fclose(file_residu) ;

  sprintf(name_tab,"dif_x") ;
  statis (dif_x, *nmax, mouch_residu, name_tab) ;
  sprintf(name_tab,"dif_y") ;
  statis (dif_y, *nmax, mouch_residu,name_tab) ;
  sprintf(name_tab,"residu") ;
  statis (residu, *nmax, mouch_residu,name_tab) ;
#endif

 
/* ----------------------------------------------------------- */   
/* -----   triage des points  apres calcul des residus   ----- */
/* ----------------------------------------------------------- */   
  if (nb_mauvais != 0)
  {
     l = 0 ;
     for (i = 0; i < *nmax; i++) 
     {
       if (x1_s[i] != 0 && y1_s[i] !=0)
       {
         x0_s[l] = x0_s[i];
         y0_s[l] = y0_s[i];
         x1_s[l] = x1_s[i];
         y1_s[l] = y1_s[i];
         residu[l]  = residu[i] ;
         dif_x[l]   = dif_x[i]  ;
         dif_y[l]   = dif_y[i]  ;
         l++;
       }
     }
     *nmax = l;
     goto recommencer;
  }

 
/* ---------------------------------------------------------------------- */
/* -----   triage pour elimination des points redondants            ----- */
/* ---------------------------------------------------------------------- */

  recommencer_1 :
  nb_mauvais_redondant = 0 ;

  for (i = 0 ; i < (*nmax - 1) ; i++) 
  {
     a00 = y0_s[i] ;
     b00 = x0_s[i] ;
     res0 = residu[i] ;
     if (a00 != 0 && b00 != 0)
     {
         for (j = i+1 ; j < *nmax ; j++) 
         {
            a11 = y0_s[j] ;
            b11 = x0_s[j] ;
            res1 = residu[j] ;
            if ((a00 == a11) && (b00 == b11))
            {
                if (res0 <= res1)
                {
                    y1_s[j] = 0 ;
                    x1_s[j] = 0 ;
                    y0_s[j] = 0 ;
                    x0_s[j] = 0 ;
                    residu[j] = 0 ;
                    nb_mauvais_redondant += 1 ;
                 }
                 else
                 {
                     y1_s[i] = 0 ;
                     x1_s[i] = 0 ;
                     y0_s[i] = 0 ;
                     x0_s[i] = 0 ;
                     residu[i] = 0 ;
                     nb_mauvais_redondant += 1 ;
                  }
              }
         }
     } 
  }

if (nb_mauvais_redondant != 0)
 {
  /* ---------------------------------------------------------------- */   
  /* -----   triage des points apres reperage des redondances   ----- */
  /* ---------------------------------------------------------------- */   

  l = 0 ;

  for (i = 0; i < *nmax; i++) 
  {
    if (x1_s[i] != 0 && y1_s[i] !=0)
    {
      x0_s[l] = x0_s[i];
      y0_s[l] = y0_s[i];
      x1_s[l] = x1_s[i];
      y1_s[l] = y1_s[i];
      residu[l]  = residu[i] ;
      dif_x[l]   = dif_x[i]  ;
      dif_y[l]   = dif_y[i]  ;
      l++;
    }
  }

 *nmax = l;

 goto recommencer_1 ;

 }   

               
  /* ---------------------------------------------------------- */   
  /* -----  fin du triage apres reperage des redondances ------ */
  /* ---------------------------------------------------------- */   


  *x0_sf = new float[*nmax] ;
  *y0_sf = new float[*nmax] ;
  *x1_sf = new float[*nmax] ;
  *y1_sf = new float[*nmax] ;



  for (i=0; i< (*nmax); i++)
  {
    (*x0_sf)[i] = x0_s[i] ;
    (*y0_sf)[i] = y0_s[i] ;
    (*x1_sf)[i] = x1_s[i] ;
    (*y1_sf)[i] = y1_s[i] ;
  }

 delete [] x0_s ;
 delete [] y0_s ;
 delete [] x1_s ;
 delete [] y1_s ;




/* ---------------------------------------------------------------------- */
/* ----  calcul de l'equation finale de mise en correspondance      ----- */
/* ---------------------------------------------------------------------- */


  moindre_carre(*nmax, *x0_sf, *y0_sf, *x1_sf, *y1_sf,
		type_equa,
                &a1x, &b1x, &c1x, &d1x, &e1x, &f1x, &g1x, &h1x, &i1x, &j1x,
                &a1y, &b1y, &c1y, &d1y,	&e1y, &f1y, &g1y, &h1y, &i1y, &j1y);

  
   par->ax = a1x ;
   par->bx = b1x ;
   par->cx = c1x ;
   par->dx = d1x ;
   par->ex = e1x ;
   par->fx = f1x ;
   par->gx = g1x ;
   par->hx = h1x ;
   par->ix = i1x ;
   par->jx = j1x ;


   par->ay = a1y ;
   par->by = b1y ;
   par->cy = c1y ;
   par->dy = d1y ;
   par->ey = e1y ;
   par->fy = f1y ;
   par->gy = g1y ;
   par->hy = h1y ;
   par->iy = i1y ;
   par->jy = j1y ;


/* ---------------------------------------------------------------------- */
/* ----  fin du calcul de l'equation finale                         ----- */
/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */
/* ----  Calcul des residus issus de l'equation finale de           ----- */
/* ----  mise en correspondance geometrique                         ----- */
/* ---------------------------------------------------------------------- */


  somme_residu = 0 ;
  somme_dif2   = 0 ;



  if (type_equa == T_EQUA_0) 
  {
    for (i = 0; i < *nmax; i++)
    {
      x01_r[i] = (float)((par->ax) * (*x0_sf)[i] - 
                 (par->bx) * (*y0_sf)[i] + (par->cx)) ;
      y01_r[i] = (float)((par->bx) * (*x0_sf)[i] + 
                 (par->ax) * (*y0_sf)[i] + (par->cy)) ;
      dif_x[i] = (float)((*x1_sf)[i]) - x01_r[i];
      dif_y[i] = (float)((*y1_sf)[i]) - y01_r[i];
      residu[i] = (float)(sqrt((double)(dif_x[i] * dif_x[i] + 
                                        dif_y[i] * dif_y[i])));
      somme_residu += residu[i] ;
      somme_dif2   += (float)( (dif_x[i] * dif_x[i]) +
                               (dif_y[i] * dif_y[i]) ) ;
    }  
  }

 
 

  if (type_equa == T_EQUA_1)
  {
    for (i = 0; i < *nmax; i++) 
    {
      x01_r[i] = (float)((par->ax) * (*x0_sf)[i] + 
                 (par->bx) * (*y0_sf)[i] + (par->cx)) ;
      y01_r[i] = (float)((par->ay) * (*x0_sf)[i] + 
                 (par->by) * (*y0_sf)[i] + (par->cy)) ;
      dif_x[i] = (float)((*x1_sf)[i] - x01_r[i]);
      dif_y[i] = (float)((*y1_sf)[i] - y01_r[i]);
      residu[i] = (float)(sqrt((double)(dif_x[i] * dif_x[i]  + 
                                        dif_y[i] * dif_y[i])));
      somme_residu += residu[i] ;
      somme_dif2   += (float)( (dif_x[i] * dif_x[i]) + (dif_y[i] * dif_y[i]) ) ;
    }
  }

  
  

  if (type_equa == T_EQUA_2) 
  {
    for (i = 0; i < *nmax; i++)
    {

      col2 = (*x0_sf)[i] * (*x0_sf)[i];
      lig2 = (*y0_sf)[i] * (*y0_sf)[i];
      x01_r[i] = (float)((par->ax) * col2 + (par->bx) * lig2 + 
                 (par->cx) * (*y0_sf)[i] * (*x0_sf)[i] +
                 (par->dx) * (*x0_sf)[i] + 
                 (par->ex) * (*y0_sf)[i] + (par->fx)) ;
      y01_r[i] = (float)((par->ay) * col2 + (par->by) * lig2 + 
                 (par->cy) * (*y0_sf)[i] * (*x0_sf)[i] +
	         (par->dy) * (*x0_sf)[i] + 
                 (par->ey) * (*y0_sf)[i] + (par->fy)) ;
      dif_x[i] = (float)((*x1_sf)[i]) - x01_r[i];
      dif_y[i] = (float)((*y1_sf)[i]) - y01_r[i];
      residu[i] = (float)(sqrt((double)(dif_x[i] * dif_x[i] + 
                                        dif_y[i] * dif_y[i])));
      if (residu[i] > seuil_residu) 
      somme_residu +=  residu[i] ;
      somme_dif2   += (float)( (dif_x[i] * dif_x[i]) + 
                               (dif_y[i] * dif_y[i]) ) ;
    }
  }
  
 

  if (type_equa == T_EQUA_3) 
  {
    for (i = 0; i < *nmax; i++)
    {
      col2 = (*x0_sf)[i] * (*x0_sf)[i];
      lig2 = (*y0_sf)[i] * (*y0_sf)[i];
      col3 = col2 * (*x0_sf)[i];
      lig3 = lig2 * (*y0_sf)[i];
      x01_r[i] = (par->ax) * col3  +  (par->bx) * lig3  +  (par->cx) * col2  +
                 (par->dx) * lig2  +  (par->ex) * (*x0_sf)[i] * lig2  +
                 (par->fx) * (*y0_sf)[i] * col2  +  
                  (par->gx) * (*y0_sf)[i] * (*x0_sf)[i]  +
                 (par->hx) * (*x0_sf)[i] + (par->ix) * (*y0_sf)[i] + (par->jx) ;
      y01_r[i] = (par->ay) * col3  +  (par->by) * lig3  +  (par->cy) * col2  +
                 (par->dy) * lig2  +  (par->ey) * (*x0_sf)[i] * lig2  +
                 (par->fy) * (*y0_sf)[i] * col2  +  
                 (par->gy) * (*y0_sf)[i] * (*x0_sf)[i]  +
                 (par->hy) * (*x0_sf)[i] + (par->iy) * (*y0_sf)[i] + (par->jy) ;
      dif_x[i] = (float)((*x1_sf)[i]) - x01_r[i];
      dif_y[i] = (float)((*y1_sf)[i]) - y01_r[i];
      residu[i] = (float)(sqrt((double)(dif_x[i] * dif_x[i] + 
                                        dif_y[i] * dif_y[i])));
      somme_residu += residu[i] ;
      somme_dif2   += (float)( (dif_x[i] * dif_x[i]) + 
                               (dif_y[i] * dif_y[i]) ) ;

    }
  }

  RMSE = sqrt(somme_dif2 / (float)(*nmax) ) ;




/* ---------------------------------------------------------------------- */
/* -----                mise en correspondance finale               ----- */
/* ---------------------------------------------------------------------- */
 
   
#ifdef MOUCHARD
   file_residu = fopen(mouch_residu,"a") ;

   fprintf(file_residu,"\n") ;  
   fprintf(file_residu,"\n") ;  
   fprintf(file_residu,"\n") ;  
   fprintf(file_residu,"mise en correspondance finale \n") ;  
   fprintf(file_residu,"Resolution    : %d\n",resolution) ;
   fprintf(file_residu,"----------------------------\n") ;
   fprintf(file_residu,"nb_points     =  %d\n",(*nmax)) ;
   fprintf(file_residu,"nb_mauvais    =  %d\n",nb_mauvais) ;
   fprintf(file_residu,"iteration     =  %d\n",iteration) ;
   fprintf( file_residu,"nombre de points redondants   :%d\n",
            nb_mauvais_redondant) ;
   fprintf(file_residu,"root mean square error     =  %f\n",RMSE) ;

   fclose(file_residu) ;
   
   sprintf(name_tab,"dif_x") ;
   statis (dif_x, *nmax, mouch_residu, name_tab) ;
   sprintf(name_tab,"dif_y") ;
   statis (dif_y, *nmax, mouch_residu,name_tab) ;
   sprintf(name_tab,"residu") ;
   statis (residu, *nmax, mouch_residu,name_tab) ;
#endif


if ( resolution == res_min && nature_rec == 1)
  { 

  /* ------------------------------------------------------------------------ */
  /* ----  transformations des points d'amers des coordonnees pixels   ------ */
  /* ----  aux coordonnees reelles                                     ------ */
  /* ------------------------------------------------------------------------ */

  trans_coor_pxl_reel_tab(depart_br_x, depart_br_y, pas_br_x, pas_br_y,
                          *x0_sf, *y0_sf, *x1_sf, *y1_sf,
                          &x0_cr_s, &y0_cr_s, &x1_cr_s, &y1_cr_s,
                          *nmax ) ;



  /* ---------------------------------------------------------------------- */
  /* ----  calcul de l'equation finale de mise en correspondance      ----- */
  /* ----  en coordonnees reelles                                     ----- */
  /* ---------------------------------------------------------------------- */
   
  moindre_carre_f(*nmax, x0_cr_s, y0_cr_s, x1_cr_s, y1_cr_s,
		  type_equa,
                  &a2x, &b2x, &c2x, &d2x, &e2x, &f2x, &g2x, &h2x, &i2x, &j2x,
                  &a2y, &b2y, &c2y, &d2y, &e2y, &f2y, &g2y, &h2y, &i2y, &j2y);

  
   par_cr->ax = a2x ;
   par_cr->bx = b2x ;
   par_cr->cx = c2x ;
   par_cr->dx = d2x ;
   par_cr->ex = e2x ;
   par_cr->fx = f2x ;
   par_cr->gx = g2x ;
   par_cr->hx = h2x ;
   par_cr->ix = i2x ;
   par_cr->jx = j2x ;


   par_cr->ay = a2y ;
   par_cr->by = b2y ;
   par_cr->cy = c2y ;
   par_cr->dy = d2y ;
   par_cr->ey = e2y ;
   par_cr->fy = f2y ;
   par_cr->gy = g2y ;
   par_cr->hy = h2y ;
   par_cr->iy = i2y ;
   par_cr->jy = j2y ;


  /* ---------------------------------------------------------------------- */
  /* ----  fin du calcul de l'equation finale en cordonnees reelles   ----- */
  /* ---------------------------------------------------------------------- */


  


  /* ---------------------------------------------------------------------- */
  /* ----  calcul des residus en cordonnees reelles                   ----- */
  /* ---------------------------------------------------------------------- */

  somme_dif2_cr   = 0 ;




  if (type_equa == T_EQUA_0) 
  {
    for (i = 0; i < *nmax; i++)
    {
      x01_r[i] = (float)((par_cr->ax) * x0_cr_s[i] - 
                 (par_cr->bx) * y0_cr_s[i] + (par_cr->cx)) ;
      y01_r[i] = (float)((par_cr->bx) * x0_cr_s[i] + 
                 (par_cr->ax) * y0_cr_s[i] + (par_cr->cy)) ;
      dif_x[i] = (float)(x1_cr_s[i]) - x01_r[i];
      dif_y[i] = (float)(y1_cr_s[i]) - y01_r[i];
      somme_dif2_cr   += (float)( (dif_x[i] * dif_x[i]) + 
                         (dif_y[i] * dif_y[i]) ) ;
    }
  }

 
 

  if (type_equa == T_EQUA_1)
  {
    for (i = 0; i < *nmax; i++) 
    {
      x01_r[i] = (float)((par_cr->ax) * x0_cr_s[i] + 
                 (par_cr->bx) * y0_cr_s[i] + (par_cr->cx)) ;
      y01_r[i] = (float)((par_cr->ay) * x0_cr_s[i] + 
                 (par_cr->by) * y0_cr_s[i] + (par_cr->cy)) ;
      dif_x[i] = (float)(x1_cr_s[i] - x01_r[i]);
      dif_y[i] = (float)(y1_cr_s[i] - y01_r[i]);
      somme_dif2_cr   += (float)( (dif_x[i] * dif_x[i]) + 
                                  (dif_y[i] * dif_y[i]) ) ;
    }
  }

  
  

  if (type_equa == T_EQUA_2) 
  {
    for (i = 0; i < *nmax; i++)
    {

      col2 = x0_cr_s[i] * x0_cr_s[i];
      lig2 = y0_cr_s[i] * y0_cr_s[i];
      x01_r[i] = (float)((par_cr->ax) * col2 + (par_cr->bx) * lig2 + 
                 (par_cr->cx) * y0_cr_s[i] * x0_cr_s[i] +
                 (par_cr->dx) * x0_cr_s[i] + (par_cr->ex) * y0_cr_s[i] + 
                 (par_cr->fx)) ;
      y01_r[i] = (float)((par_cr->ay) * col2 + (par_cr->by) * lig2 + 
                 (par_cr->cy) * y0_cr_s[i] * x0_cr_s[i] +
	         (par_cr->dy) * x0_cr_s[i] + (par_cr->ey) * y0_cr_s[i] + 
                 (par_cr->fy)) ;
      dif_x[i] = (float)(x1_cr_s[i]) - x01_r[i];
      dif_y[i] = (float)(y1_cr_s[i]) - y01_r[i];
      somme_dif2_cr   += (float)( (dif_x[i] * dif_x[i]) + 
                                  (dif_y[i] * dif_y[i]) ) ;
    }
  }
  
 

  if (type_equa == T_EQUA_3) 
  {
    for (i = 0; i < *nmax; i++)
    {
      col2 = x0_cr_s[i] * x0_cr_s[i];
      lig2 = y0_cr_s[i] * y0_cr_s[i];
      col3 = col2 * x0_cr_s[i];
      lig3 = lig2 * y0_cr_s[i];
      x01_r[i] = (par_cr->ax) * col3  +  (par_cr->bx) * lig3  +  
                 (par_cr->cx) * col2  +
                 (par_cr->dx) * lig2  +  (par_cr->ex) * x0_cr_s[i] * lig2  +
                 (par_cr->fx) * y0_cr_s[i] * col2  +  
                 (par_cr->gx) * y0_cr_s[i] * x0_cr_s[i]  +
                 (par_cr->hx) * x0_cr_s[i] + (par_cr->ix) * y0_cr_s[i] + 
                 (par_cr->jx) ;
      y01_r[i] = (par_cr->ay) * col3  +  (par_cr->by) * lig3  +  
                 (par_cr->cy) * col2  +
                 (par_cr->dy) * lig2  +  (par_cr->ey) * x0_cr_s[i] * lig2  +
                 (par_cr->fy) * y0_cr_s[i] * col2  +  
                 (par_cr->gy) * y0_cr_s[i] * x0_cr_s[i]  +
                 (par_cr->hy) * x0_cr_s[i] + (par_cr->iy) * y0_cr_s[i] + 
                 (par_cr->jy) ;
      dif_x[i] = (float)(x1_cr_s[i]) - x01_r[i];
      dif_y[i] = (float)(y1_cr_s[i]) - y01_r[i];
      somme_dif2_cr   += (float)( (dif_x[i] * dif_x[i]) + 
                                  (dif_y[i] * dif_y[i]) ) ;

    }
  }


  RMSE_cr = sqrt(somme_dif2_cr / (float)(*nmax) ) ;

#ifdef MOUCHARD    
  file_residu = fopen(mouch_residu,"a") ;

  fprintf(file_residu,"\n") ;
  fprintf(file_residu,"\n") ;
  fprintf(file_residu,"\n") ;
  fprintf(file_residu,"\n") ;
  fprintf(file_residu,"Residu en coordonnees reelles\n") ;
  fprintf(file_residu,"Resolution    : %d\n",resolution) ;
  fprintf(file_residu,"----------------------------\n") ;
  fprintf(file_residu,"nb_points     =  %d\n",(*nmax)) ;
  fprintf(file_residu,"root mean square error     =  %f\n",RMSE_cr) ;

  fclose(file_residu) ;

 sprintf(name_tab,"dif_x") ;
  statis (dif_x, *nmax, mouch_residu, name_tab) ;
  sprintf(name_tab,"dif_y") ;
  statis (dif_y, *nmax, mouch_residu,name_tab) ;
#endif

  }
/* ========================== */


  delete  [] dif_x ;
  delete  [] dif_y ;
  delete  [] x01_r ;
  delete  [] y01_r ;
  delete  [] residu ;
 
}

/**************************************************************************/

void detection( Ifloat &Im_seuil, float **x_ef, float **y_ef, float **f_ef, int *nmax)
{

  float    *x_e, *y_e;
  float	   *f_e;

  int	nb_lig, nb_col ;
  int   *suivi0_x, *suivi0_y, *suivi_x, *suivi_y;
  int   lig, col, lig1, lig2, lig3;
  int   i, taille_allocation ;


 taille_allocation = 60000 ;

/* -- nombre de lignes et de colonnes de l'image seuillee -- */

 nb_lig = Im_seuil.nl() ;
 nb_col = Im_seuil.nc() ;

/* -----        allocation dynamique 1 D           ------------ */

 suivi0_x = new int[nb_col] ;
 suivi0_y = new int[nb_col] ;
 suivi_x  = new int[nb_col] ;
 suivi_y  = new int[nb_col] ;

 x_e = new float[taille_allocation] ;
 y_e = new float[taille_allocation] ;
 f_e = new float[taille_allocation] ;


/* -----      fin de l'allocation dynamique        ------------ */


  *nmax = 0;

  /* ------------------------------------------------------------------ */
  /* ----       initialisation des compteurs de suivie             ---- */
  /* ------------------------------------------------------------------ */

  for (i = 0; i < nb_col; i++) 
  {
    suivi0_x[i] = 1;
    suivi0_y[i] = 1;
  }

  /* ------------------------------------------------------------------ */
  /* ----       fin de l'initialisation                            ---- */
  /* ------------------------------------------------------------------ */



  for (lig = 0; lig < nb_lig; lig++)
  {

    lig1 = lig - 1;
    lig2 = lig;
    lig3 = lig + 1;

    if (lig1 < 0)
      lig1 = lig;

    if (lig3 >=  nb_lig)
      lig3 = lig;


  /* ------------------------------------------------------------------ */
  /* -----       debut de la procedure de suivie                   ---- */
  /* ------------------------------------------------------------------ */


    for (col = 0; col < nb_col; col++)
    {
      if (col == 0) 
      {

	suivi_x[col] = 1;
	if (lig == 1) 
        {
	  suivi_y[col] = 1;
	  goto test;
	}


 /*a*/  if (Im_seuil(lig2,col) < Im_seuil(lig1,col))
        {
	  suivi_y[col] = 0;
	  goto cont_col;
	}

	
 /*b*/	if (Im_seuil(lig2,col) > Im_seuil(lig1,col)) 
        {
	  suivi_y[col] = 1;
	  goto test;
	}

	
 /*c*/	if (Im_seuil(lig2,col) == Im_seuil(lig1,col)) 
        {
	  if (suivi0_y[col] == 0) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  } 
          else 
          {
	    suivi_y[col] = 1;
	    goto test;
	  }
	}
      }


      if (col != 0) 
      {
	
  /*1*/	 if (Im_seuil(lig2,col) < Im_seuil(lig2,col-1)) 
         {
	  suivi_x[col] = 0;
	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    goto cont_col;
	  }

   /*a*/  if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }

	  
   /*b*/  if (Im_seuil(lig2,col) > Im_seuil(lig1,col)) 
          {
	    suivi_y[col] = 1;
	    goto cont_col;
	  }

	  
   /*c*/  if (Im_seuil(lig2,col) == Im_seuil(lig1,col)) 
          {
	    if (suivi0_y[col] == 1)
	      suivi_y[col] = 1;
	    else
	      suivi_y[col] = 0;
	    goto cont_col;
	  }

	}


	
/*2*/	if (Im_seuil(lig2,col) == Im_seuil(lig2,col-1)) 
        {
	  if (suivi_x[col - 1] == 1)
	    suivi_x[col] = 1;
	  else
	    suivi_x[col] = 0;

	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    if (suivi_x[col] == 1)
	      goto test;
	    else
	      goto cont_col;
	  }

	 
   /*a*/  if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }

	  
   /*b*/  if (Im_seuil(lig2,col) > Im_seuil(lig1,col)) 
          {
	    suivi_y[col] = 1;
	    if (suivi_x[col] == 0)
	      goto cont_col;
	    else
	      goto test;
	  }

	  
  /*c*/	  if (Im_seuil(lig2,col) == Im_seuil(lig1,col)) 
          {
	    if (suivi0_y[col] == 0) 
            {
	      suivi_y[col] = 0;
	      goto cont_col;
	    } 
            else 
            {
	      suivi_y[col] = 1;
	      if (suivi_x[col] == 1)
		goto test;
	      else
		goto cont_col;
	    }
	  }

	}



	
/*3*/	if (Im_seuil(lig2,col) > Im_seuil(lig2,col-1))
        {
	  suivi_x[col] = 1;
	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    goto test;
	  }


	  
   /*a*/ if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
         {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }

	  
  /*b*/	  if (Im_seuil(lig2,col) > Im_seuil(lig1,col))
          {
	    suivi_y[col] = 1;
	    goto test;
	  }

	  
  /*c*/	  if (Im_seuil(lig2,col) == Im_seuil(lig1,col)) 
          {
	    if (suivi0_y[col] == 0) 
            {
	      suivi_y[col] = 0;
	      goto cont_col;
	    } 
            else 
            {
	      suivi_y[col] = 1;
	      goto test;
	    }
	  }

	}

      }



test:


  /* ------------------------------------------------------------------ */
  /*----         debut de la procedure de test                      --- */
  /* ------------------------------------------------------------------ */

      if (col != 0) 
      {
	if (col != (nb_col - 1)) 
        {
         /* ----------------------------------------------------------- */
	 /* ---    test sur les elements horizontaux et verticaux   --- */
         /* ----------------------------------------------------------- */

	  if (Im_seuil(lig2,col) < Im_seuil(lig2,col-1)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
            goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig2,col+1)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig3,col)) 
	    goto cont_col;

         /* ----------------------------------------------------------- */
	 /* ---            test sur les elements diagonaux          --- */
         /* ----------------------------------------------------------- */

	  if ( (Im_seuil(lig2,col) < Im_seuil(lig1,col-1)) && 
               (Im_seuil(lig2,col) < Im_seuil(lig1,col+1)) ) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) < Im_seuil(lig1,col-1))
	    goto cont_col;
	  if (Im_seuil(lig2,col) < Im_seuil(lig1,col+1))
	    goto cont_col;
	  if ( (Im_seuil(lig2,col) <= Im_seuil(lig3,col-1)) && 
               (Im_seuil(lig2,col) <= Im_seuil(lig3,col+1)) )
	    goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig3,col-1)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig3,col+1)) 
	    goto cont_col;

	} 
        else 
       {

         /* ----------------------------------------------------------- */
	 /* ---    test sur les elements horizontaux et verticaux   --- */
         /* ----------------------------------------------------------- */

	  if (Im_seuil(lig2,col) < Im_seuil(lig2,col-1)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig3,col)) 
	    goto cont_col;

         /* ----------------------------------------------------------- */
	 /* ---           test sur les elements diagonaux           --- */
         /* ----------------------------------------------------------- */
	  if (Im_seuil(lig2,col) < Im_seuil(lig1,col-1)) 
	    goto cont_col;
	  if (Im_seuil(lig2,col) <= Im_seuil(lig3,col-1)) 
	    goto cont_col;
          }
        } 
        else 
        {

        /* ----------------------------------------------------------- */
	/* ---    test sur les elements horizontaux et verticaux   --- */
        /* ----------------------------------------------------------- */
	if (Im_seuil(lig2,col) < Im_seuil(lig1,col)) 
	  goto cont_col;
	if (Im_seuil(lig2,col) <= Im_seuil(lig2,col+1))
	  goto cont_col;
	if (Im_seuil(lig2,col) <= Im_seuil(lig3,col)) 
	  goto cont_col;

        /* ----------------------------------------------------------- */
	/* ---           test sur les elements diagonaux           --- */
        /* ----------------------------------------------------------- */
	if (Im_seuil(lig2,col) < Im_seuil(lig1,col+1)) 
	  goto cont_col;
	if (Im_seuil(lig2,col) <= Im_seuil(lig3,col+1)) 
	  goto cont_col;

      }




      x_e[*nmax] = col;
      y_e[*nmax] = lig;
      f_e[*nmax] = Im_seuil(lig2,col);



      (*nmax)++;


  /* ------------------------------------------------------------------ */  
  /* ----           fin de la procedure de test                    ---- */
  /* ------------------------------------------------------------------ */

  if ((*nmax) == taille_allocation)
    {
       taille_allocation += 10000 ;
       
       x_e = (float *) realloc( x_e, taille_allocation * sizeof(float)) ;
       y_e = (float *) realloc( y_e, taille_allocation * sizeof(float)) ;

       f_e = (float *) realloc( f_e, taille_allocation * sizeof(float)) ;
    }


cont_col: ;

    }  /* end col */


    for (i = 0; i < nb_col; i++)
      suivi0_y[i] = suivi_y[i];


  }  /* end ligne */
 delete [] suivi0_x ;
 delete [] suivi0_y ;
 delete [] suivi_x ;
 delete [] suivi_y ;

 *x_ef = new float[*nmax] ;
 *y_ef = new float[*nmax] ;
 *f_ef = new float[*nmax] ;
 
 for (i=0; i< (*nmax); i++)
   {
     (*x_ef)[i] = x_e[i] ;
     (*y_ef)[i] = y_e[i] ;
     (*f_ef)[i] = f_e[i] ;
   }
delete [] x_e ;
delete [] y_e ;
delete [] f_e ;

}

/* End. */

/**************************************************************************/

void recalage(parametre  *par, int type_equa, int type_rec, float a,
              Ifloat &Im_in, Ifloat &Im_out)

{

  /* ----------------------------------------------- */
  /* ------   mapping image output to input   ------ */
  /* ------   for float images only           ------ */
  /* ----------------------------------------------- */

  double          x=0., y=0.;
  float	          t1, t2, t3, t4, a1, a2, a3, a4, b1, b2, b3, b4, val ;
  float	  	  delta_x, delta_x_2, delta_x_3, delta_y, delta_y_2, delta_y_3 ;

  int             nb_lig, nb_col ;
  int Nl = Im_out.nl();
  int Nc = Im_out.nc();
  float    val1 ;
  long            nbl0;
  long            col, ligne,
                  col2, ligne2,
                  col3, ligne3,
                  iix, iiy, 
		  iix1, iiy1 ,
		  iix2, iiy2 ,
		  iix0, iiy0 ;


  /*######################################################################## */
  /*#########                Traitement                             ######## */
  /*######################################################################## */

  nb_lig = Im_in.nl() ;
  nb_col = Im_in.nc() ;
 
  /* ----------------------------------------------------------------------- */
  /* ---    initialisation de l'image recale a zero                      --- */
  /* ---    l'image recale aura la meme dimension que                    --- */
  /* ---    l'image de reference                                         --- */
  /* ----------------------------------------------------------------------- */

  nbl0 = nb_col;



/* calcul des coordonnees du point de l'image de reference 
   correspondant au point de l'image final */

 for (ligne = 0; ligne < Nl; ligne++) 
   {

     for (col = 0; col < Nc; col++) 
       {

         if (type_equa == T_EQUA_0) 
          {
	    x = par->ax * col - par->bx * ligne + par->cx;
	    y = par->bx * col + par->ax * ligne + par->cy;
          } 

         if (type_equa == T_EQUA_1) 
          {
	    x = par->ax * col + par->bx * ligne + par->cx;
	    y = par->ay * col + par->by * ligne + par->cy;
          }


         if (type_equa == T_EQUA_2) 
          {
	    col2 = col * col;
	    ligne2 = ligne * ligne;
	    x = par->ax * col2 + 
                par->bx * ligne2 + 
                par->cx * ligne * col + 
                par->dx * col + 
                par->ex * ligne + 
                par->fx;
	    y = par->ay * col2 + 
                par->by * ligne2 + 
                par->cy * ligne * col + 
                par->dy * col + 
                par->ey * ligne + 
                par->fy;
           }


          if (type_equa == T_EQUA_3) 
           {
    
             col2 = col  * col ;
             col3 = col2 * col ;

             ligne2 = ligne * ligne ;
             ligne3 = ligne2 * ligne;

             x = par->ax * col3 		+ 
                 par->bx * ligne3 		+ 
                 par->cx * col2 		+ 
		 par->dx * ligne2 		+ 
		 par->ex * col * ligne2 	+ 
          	 par->fx * col2 * ligne 	+ 
		 par->gx * col * ligne	+ 
		 par->hx * col 		+ 
		 par->ix * ligne 		+ 
		 par->jx ;


             y = par->ay * col3 		+ 
		 par->by * ligne3 		+ 
		 par->cy * col2 		+
		 par->dy * ligne2 		+
		 par->ey * col * ligne2 	+
          	 par->fy * col2 * ligne 	+
		 par->gy * col * ligne 	+
		 par->hy * col 		+
		 par->iy * ligne 		+
		 par->jy;


           }
  



           /* calcul de la valeur interpole a assigne au pixel de l'image a recale */

           iix = (int) floor(x) ;
           iiy = (int) floor(y) ;

           
           if (iix >= nb_col) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iix < 0) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iiy >= nb_lig) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

           if (iiy < 0) 
            {
	      Im_out(ligne,col) = 0;
	      goto suite;
            }

  
           delta_x = x - iix ;
           delta_y = y - iiy ;


           /* cas de l'interpolation du plus proche voisin */

           if (type_rec == T_REC_CLOSER_PIXEL) 
            {
	      Im_out(ligne,col) = Im_in(iiy,iix) ;
            }



          /* cas de l'interpolation bilineaire */
          if (type_rec == T_REC_BILINEAR_INTERPOL)
           {
             iiy1 = iiy + 1 ;
             iix1 = iix + 1 ;
             if (iiy1 >= nb_lig)
              {
                iiy1 = nb_lig - iiy1 + nb_lig -1 ;
              }
             if (iix1 >= nb_col)
	      {
                iix1 = nb_col - iix1 + nb_col -1;
              }


	     t1 = ((1.0 - delta_x) * (double)(Im_in(iiy,iix))  ) + 
		  ((double)((Im_in(iiy,iix1))) * delta_x) ;
	     t2 = ((double)(Im_in(iiy1,iix)) * (1.0 - delta_x)) + 
		  ((double)(Im_in(iiy1,iix1)) * delta_x) ;
             val = (t1 * (1. - delta_y)) + (t2 * delta_y) ;

             /* if (val >= 255)
                val = 255 ;
             if (val <= 0)
                val = 0 ; */
            // val1 = (irint(val)) ;
             val1 = val ;
 
             Im_out(ligne,col) = val1 ;
	   }



         /* cas de la convolution cubique avec parametre a */

         if (type_rec == T_REC_CUBIC_INTERPOL)
          {
            delta_x_2 = delta_x * delta_x ;
            delta_x_3 = delta_x_2 * delta_x ;
            delta_y_2 = delta_y * delta_y ;
            delta_y_3 = delta_y_2 * delta_y ;

	    a1 = (a-2)*(delta_x_3 - 2*delta_x_2 + delta_x) ;
            a2 = a * delta_x_3 - (1+a)*delta_x_2 + 1 ;
            a3 = -a*delta_x_3 + (2*a-1)*delta_x_2 + (2-a)*delta_x ;
            a4 = (a-2) * (delta_x_2 - delta_x_3) ;

	    b1 = (a-2)*(delta_y_3 - 2*delta_y_2 + delta_y) ;
            b2 = a * delta_y_3 - (1+a)*delta_y_2 + 1 ;
            b3 = -a*delta_y_3  + (2*a - 1)*delta_y_2 + (2-a)*delta_y ;
            b4 = (a-2) * (delta_y_2 - delta_y_3) ;

	    iiy0 = iiy - 1 ;
	    iiy1 = iiy + 1 ;
	    iiy2 = iiy + 2 ;

            iix0 = iix - 1 ;
            iix1 = iix + 1 ;
            iix2 = iix + 2 ;

            if (iix0 < 0)
               iix0 = iix0 - 2 * iix0 - 1 ;
            if (iix1 >= nb_col)
               iix1 = 2 * nb_col - iix1 -1 ;
            if (iix2 >= nb_col)
               iix2 = 2 * nb_col - iix2 -1 ;

            if (iiy0 < 0)
               iiy0 = iiy0 - 2 * iiy0 - 1 ;
            if (iiy1 >= nb_lig)
               iiy1 = 2 * nb_lig - iiy1 -1 ;
            if (iiy2 >= nb_lig)
               iiy2 = 2 * nb_lig - iiy2 -1 ;

            t1 = (double)(Im_in(iiy0,iix0))*a1 + 
		 (double)(Im_in(iiy0,iix))*a2  + 
		 (double)(Im_in(iiy0,iix1))*a3 + 
		 (double)(Im_in(iiy0,iix2))*a4 ;
            t2 = (double)(Im_in(iiy,iix0))*a1  + 
		 (double)(Im_in(iiy,iix))*a2   + 
		 (double)(Im_in(iiy,iix1))*a3  + 
		 (double)(Im_in(iiy,iix2))*a4 ;
            t3 = (double)(Im_in(iiy1,iix0))*a1 + 
		 (double)(Im_in(iiy1,iix))*a2  + 
		 (double)(Im_in(iiy1,iix1))*a3 + 
		 (double)(Im_in(iiy1,iix2))*a4 ;
            t4 = (double)(Im_in(iiy2,iix0))*a1 + 
		 (double)(Im_in(iiy2,iix))*a2  + 
		 (double)(Im_in(iiy2,iix1))*a3 +
		 (double)(Im_in(iiy2,iix2))*a4 ;


	    val = (t1*b1 + t2*b2 + t3*b3 + t4*b4) ;
            /* if (val >= 255)
               val = 255 ;
            if (val <= 0)
               val = 0 ;
            val1 = (irint(val)) ; */
            val1 = val ; 
	    Im_out(ligne,col) = val1 ;

           }



          suite: ;
  
      } /* end col */
   } /* end ligne */
} /* End. */


/**************************************************************************/
