
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
**    File:  MR_Fusion.h
**
*******************************************************************************
**
**    DESCRIPTION  header image registration
**    -----------  
**
**
******************************************************************************/

#ifndef _FUSION_H_
#define _FUSION_H_

#define T_EQUA_0 0
#define T_EQUA_1 1
#define T_EQUA_2 2
#define T_EQUA_3 3

#define T_REC_CLOSER_PIXEL 0
#define T_REC_BILINEAR_INTERPOL 1
#define T_REC_CUBIC_INTERPOL 2

typedef struct {
 float ax, bx, cx, dx, ex, fx, gx, hx, ix, jx;
 float ay, by, cy, dy, ey, fy, gy, hy, iy, jy ;
 }  parametre ;

void detection( Ifloat &Im_seuil, float **x_ef, float **y_ef, float **f_ef, int *nmax);

void transpo_coor_max(float *x_e, float *y_e, 
                      parametre *par,
                      int nmax, int type_equa,
		      float **x_s, float **y_s);


void moindre_carre(int nbmax, 
                   float *x_0,  float *y_0,  float *x_1,  float *y_1,
	           int type_equa,
                   float *ax, float *bx, float *cx, float *dx, float *ex, 
                   float *fx, float *gx, float *hx, float *ix, float *jx,
                   float *ay, float *by, float *cy, float *dy, float *ey, 
                   float *fy, float *gy, float *hy, float *iy, float *jy  );

void gauss(double **mat, double *vec1, int dimension, double *vec2);

void trans_coor_pxl_reel_tab(int depart_x, int depart_y, 
                             float pas_x, float pas_y,
                             float *tab01_x, float *tab01_y , 
                             float *tab02_x, float *tab02_y,
                             float **tab11_x, float **tab11_y,
                             float **tab12_x, float **tab12_y,
                             int nmax  );

void moindre_carre_f(int nbmax, 
                     float *x_0,  float *y_0,  float *x_1,  float *y_1,
	             int type_equa,
                     float *ax, float *bx, float *cx, float *dx, float *ex, 
                     float *fx, float *gx, float *hx, float *ix, float *jx,
                     float *ay, float *by, float *cy, float *dy, float *ey, 
                     float *fy, float *gy, float *hy, float *iy, float *jy  );

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
				      parametre *par, parametre *par_cr );

void recalage_coor_reelle(parametre *par,
                          int type_equa, int type_rec, float a,
                          Ifloat &Im_in, Ifloat &Im_out,
                          int depart_x, int depart_y, float pas_x, float pas_y);


void recalage(parametre  *par,
              int type_equa, int type_rec, float a,
              Ifloat &Im_in, Ifloat &Im_out);


void copie_pts_fichier(char *nom_entree, 
                       float *x0, float *y0, float *x1, float *y1, 
                       int nmax);

void statis (float *tab_in, int nb_pts, char *name_mouch, char *name_tab);

#endif
