/******************************************************************************
**                   Copyright (C) 1995 OCA-NICE
*******************************************************************************
**
**    UNIT
**
**    Version: 
**
**    Author:	Frederic Rue 
**
**    Date:	02/07/1996
**    
**    File:	gauss.c
**
*******************************************************************************
**
**    DESCRIPTION Bibliotheque de procedures sur les gaussiennes	
**    
**
******************************************************************************/

#include	"IM_Obj.h"
#include	"IM_VisTool.h"

/***************************************************************************
** Gauss_2D: calcul valeur de la gaussienne au point de coordonnees (x,y) 
**  
**
****************************************************************************/

float Gauss_2D(float x,float y,t_Gauss_2D Gauss)

{
float	mx,my;
float	sX,sY;
float	teta;
float	amp;
float 	offset;
float	f;
float X,Y;

mx=Gauss.mx;
my=Gauss.my;
sX=Gauss.sX;
sY=Gauss.sY;
teta=Gauss.teta;
amp=Gauss.amp;
offset=Gauss.offset;


x -= mx;
y -= my;


X = x*cos(teta)+y*sin(teta);
Y = -x*sin(teta)+y*cos(teta);

f =offset+amp*exp(-(X*X/sX/sX+Y*Y/sY/sY)/2);

return(f);

}
/***********************************************
** Ajout_gauss_2D : ajoute une gaussienne a l'image
**
**
***********************************************/
void Ajout_gauss_2D(Ifloat& Im,t_Gauss_2D Gauss)
{
int i,j;
int Nl=Im.nl();
int Nc=Im.nc();

for(i=0;i<Nl;i++)
for(j=0;j<Nc;j++) Im(i,j) += Gauss_2D((float) j,(float) i,Gauss);


}
/***************************************************************************
** Calcul_param_gauss_2D: Calcul des parametres carateristiques d'une 
** 			  gaussienne a partir des moments d'ordre 1 et 2  
**		    	  de l'image
**
****************************************************************************/
void Calcul_param_gauss_2D(t_Moment_image Moment,float Max,t_Gauss_2D& Gauss)

{
double mx,my,vx,vy,vxy;
double teta;

Gauss.mx=mx=Moment.mx;
Gauss.my=my=Moment.my;
vx=Moment.vx;
vy=Moment.vy;
vxy=Moment.vxy;

if(fabs(vx-vy)<1E-7 || fabs(vxy)<1E-9)
	{
	Gauss.sX = sqrt(fabs(vx));
	Gauss.sY = sqrt(fabs(vy));
	Gauss.teta=0.;
	}
else
	{
	Gauss.teta=teta=0.5*atan(2*vxy/(vx-vy));

	Gauss.sX = sqrt(fabs( (vx + vy )/2 + vxy/sin(2*teta) ));
	Gauss.sY = sqrt(fabs( (vx + vy )/2 - vxy/sin(2*teta) ));
	
	}
Gauss.amp=Max;
Gauss.offset=0.;	
}

/***************************************************************************
** Estime_param_gauss_2D: Estimation les parametres carateristiques d'une 
** 			  gaussienne sensee representee l'image, 
**		    	  a partir des moments d'ordre 1 et 2 de l'image
**
****************************************************************************/
void Estime_param_gauss_2D(t_Gauss_2D& Gauss_est, Ifloat & Im)

{
t_Moment_image Moment;
float	Max;


/* Calcul des moments de l'image originale */

Moment= moments(Im);
Max=max(Im);

/* Calcul parametres de la gaussienne */ 

Calcul_param_gauss_2D(Moment,Max,Gauss_est);

}
	

/***************************************************************************
** Estime_param_gauss_2D_iter: Estimation les parametres carateristiques d'une 
** 			       gaussienne sensee representee l'image, 
**		    	       a partir des moments d'ordre 1 et 2 de l'image
**
****************************************************************************/
void Estime_param_gauss_2D_iter(t_Gauss_2D& Gauss_est,Ifloat & Im,int Nb_iter)

{
int i,j,n,pixel;
int Nl = Im.nl();
int Nc = Im.nc();
t_Moment_image Moment;
float	Max;
t_Gauss_2D Gauss;

double mx0,my0,sX0,sY0,teta0,amp0;
double mxn,myn,sXn,sYn,tetan,ampn;


/* Calcul des moments de l'image originale */

Moment= moments(Im);
Max= max(Im);

/* Calcul parametres de la gaussienne initiale */ 

Calcul_param_gauss_2D(Moment,Max,Gauss);

mx0=mxn=Gauss.mx;
my0=myn=Gauss.my;
sX0=sXn=Gauss.sX;
sY0=sYn=Gauss.sY;
teta0=tetan=Gauss.teta;
amp0=ampn=Gauss.amp;


/* Algorithme de Van-Cittert */

Ifloat Im_est(Nl,Nc, "gauss iter");

for(n=0;n<Nb_iter;n++)
	{
	for(i=0;i<Nc;i++)
		{
		for(j=0;j<Nl;j++)
			{
			pixel=i+j*Nc;
			Im_est(j,i) = Gauss_2D((float) i,(float) j,Gauss);
			}
		}

	Moment= moments(Im_est);
	Max= max(Im_est);
	Calcul_param_gauss_2D(Moment,Max,Gauss);

	Gauss.mx=mxn=mx0+mxn-Gauss.mx;
	Gauss.my=myn=my0+myn-Gauss.my;
	Gauss.sX=sXn=sX0+sXn-Gauss.sX;
	Gauss.sY=sYn=sY0+sYn-Gauss.sY;
	Gauss.teta=tetan=teta0+tetan-Gauss.teta;
	Gauss.amp=ampn=amp0+ampn-Gauss.amp;
	}

Gauss_est.mx=mxn;
Gauss_est.my=myn;
Gauss_est.sX=sXn;
Gauss_est.sY=sYn;
Gauss_est.teta=tetan;
Gauss_est.amp=ampn;
Gauss_est.offset=0;

}


/***************************************************************************
** Affiche_param_gauss_2D: Affiche les parametres carateristiques d'une 
** 			  gaussienne
**
****************************************************************************/
void Affiche_param_gauss_2D(const t_Gauss_2D& Gauss)
{
cout<<"Gaussian Estimation "<<endl;
cout<<"mx="<<Gauss.mx<<endl;
cout<<"my="<<Gauss.my<<endl;

cout<<"sX="<<Gauss.sX<<endl;
cout<<"sY="<<Gauss.sY<<endl;

cout<<"teta="<<Gauss.teta*RAD_DEG<<endl;

cout<<"amp="<<Gauss.amp<<endl;

cout<<"offset="<<Gauss.offset<<endl;
}
