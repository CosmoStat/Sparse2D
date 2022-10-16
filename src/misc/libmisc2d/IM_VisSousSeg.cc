//----------------------------------------------------------
//	Copyright (C) 1996 OCA-NICE
//----------------------------------------------------------
//
//    UNIT
//
//    Version: 
//
//    Author: 	Frederic RUE  
//
//    Date: 	02/07/1996
//    
//    File:	sous_seg.c
//
//----------------------------------------------------------
//
//    DESCRIPTION Fonctions et procedures permettant de	
//    		  sous-segmenter une image  
//	  
//----------------------------------------------------------
#include	"IM_Obj.h"

struct t_Pixel
	{
	int	num;
	int	xmax;
	int	ymax;
	int	label;
	int	label_prec;
	t_Pixel *Suivant;
	};

	/***************************************************
	** D: fonction inline de distance au carre  
	** 
	**
	***************************************************/
inline float D(int x1,int y1,int x2,int y2)
{
return float( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
}

	/***************************************************
	** libere  
	** 
	**
	***************************************************/
static void libere(t_Pixel* Max)
{
if(Max->Suivant != NULL)
	libere(Max->Suivant);

if(Max->Suivant == NULL)
	{
	delete Max;
	Max = NULL;
	}
}

	/***************************************************
	** Complete  
	** 
	**
	***************************************************/
static void Complete(t_Pixel* Max_prec1,t_Pixel* Max1)
{
t_Pixel* Max_prec;
t_Pixel* Max_int;
t_Pixel* Max;

Max_prec = Max_prec1;
do
	{
	Max_int = Max_prec;
	}
while((Max_prec = Max_prec->Suivant) != NULL);
Max_prec = Max_int;
	
Max = Max1;
do
	{
	Max_prec->Suivant = new t_Pixel;
	assert(Max_prec->Suivant != NULL);
	Max_prec = Max_prec->Suivant;
			
	Max_prec->num = Max->num;
	Max_prec->xmax = Max->xmax;
	Max_prec->ymax = Max->ymax;
	Max_prec->label = Max->label;
	Max_prec->label_prec = Max->label_prec;
	Max_prec->Suivant = NULL;
	}
while((Max = Max->Suivant) != NULL);

}

	/***************************************************
	** Test_pixel:	  
	**
	***************************************************/
static int Test_pixel(int x1,int y1,int x2,int y2)
{
if( (x1-1)==x2 && y1==y2 )
	return 1;

if( (x1+1)==x2 && y1==y2 )
	return 1;

if( (y1-1)==y2 && x1==x2 )
	return 1;

if( (y1+1)==y2 && x1==x2 )
	return 1;

if( (x1-1)==x2 && (y1-1)==y2 )
	return 1;

if( (x1+1)==x2 && (y1-1)==y2 )
	return 1;

if( (x1-1)==x2 && (y1+1)==y2 )
	return 1;

if( (x1+1)==x2 && (y1+1)==y2 )
	return 1;

return 0;
}

	/***************************************************
	** Test_pixel_bord:	  
	**
	***************************************************/
static int Test_pixel_bord(int x,int y,Ifloat& Im)
{
int Nl = Im.nl();
int Nc = Im.nc();

if((x-1)>=0 && Im(y,x-1)==0.)
	return 1;

if((x+1)<Nc && Im(y,x+1)==0.)
	return 1;

if((y-1)>=0 && Im(y-1,x)==0.)
	return 1;

if((y+1)<Nl && Im(y+1,x)==0.)
	return 1;

if((x-1)>=0 && (y-1)>=0 && Im(y-1,x-1)==0.)
	return 1;

if((x+1)<Nc && (y-1)>=0 && Im(y-1,x+1)==0.)
	return 1;

if((x-1)>=0 && (y+1)<Nl && Im(y+1,x-1)==0.)
	return 1;

if((x+1)<Nc && (y+1)<Nl && Im(y+1,x+1)==0.)
	return 1;

return 0;
}

	/***************************************************
	** Test_max:	test si (x,y) est un maximum local  
	** 		du domaine etiquete label 
	**
	***************************************************/
static int Test_max(int x,int y,int label2,Ifloat& Im,Ifloat& Im_seg2,int Nb_max,
	t_Pixel* Max01,t_Pixel* Max1)
{
int Nc = Im.nc();
int Nl = Im.nl();
t_Pixel* Max;
t_Pixel* Max0;

/* Test avec les max initiaux */

if(Max01 != NULL)
	{
	Max0 = Max01;
	do 
		{
		// Si max initial retour

		if( x == Max0->xmax && y == Max0->ymax)
			return 1;

		// Evite d'avoir un max voisin d'un max initial

		if(Test_pixel(x,y,Max0->xmax,Max0->ymax))
			return 0;
		}
	while((Max0= Max0->Suivant) != NULL);

	// Evite d'avoir un max autre qu'un max initial sur les bords

	if(Test_pixel_bord(x,y,Im))
		return 0;
	}
	

/* Evite d'avoir 2 max voisins */

if(Nb_max >= 1)
	{
	Max = Max1;
	do 
		{
		if(Test_pixel(x,y,Max->xmax,Max->ymax))
			return 0;
		}
	while((Max= Max->Suivant) != NULL);
	}	

/* Test sur les voisins */

if((x-1)>=0 && Im_seg2(y,x-1)==label2 && Im(y,x-1)>Im(y,x))
	return 0;

if((x+1)<Nc && Im_seg2(y,x+1)==label2 && Im(y,x+1)>Im(y,x))
	return 0;

if((y-1)>=0 && Im_seg2(y-1,x)==label2 && Im(y-1,x)>Im(y,x))
	return 0;

if((y+1)<Nl && Im_seg2(y+1,x)==label2 && Im(y+1,x)>Im(y,x))
	return 0;

if((x-1)>=0 && (y-1)>=0 && Im_seg2(y-1,x-1)==label2 && Im(y-1,x-1)>Im(y,x))
	return 0;

if((x+1)<Nc && (y-1)>=0 && Im_seg2(y-1,x+1)==label2 && Im(y-1,x+1)>Im(y,x))
	return 0;

if((x-1)>=0 && (y+1)<Nl && Im_seg2(y+1,x-1)==label2 && Im(y+1,x-1)>Im(y,x))
	return 0;

if((x+1)<Nc && (y+1)<Nl && Im_seg2(y+1,x+1)==label2 && Im(y+1,x+1)>Im(y,x))
	return 0;

return 1;
}

	/***************************************************
	** Cherche_max: 
	** 
	**
	****************************************************/
void Cherche_max(t_Pixel* Max01,t_Pixel* Max1,int& Nb_max,Ifloat& Im,Ifloat& Im_seg2)
{
int x,y;
int Nl = Im.nl();
int Nc = Im.nc();
int label2;
t_Pixel* Max;

/*** DETERMINE LES MAXIMA LOCAUX ***/

Max = Max1;
Nb_max = 0;
for(x=0;x<Nc;x++)
	{
	for(y=0;y<Nl;y++)
		{
		if(Im(y,x) > 0.)
			{
			label2 = (int) floor(Im_seg2(y,x));
			if(Test_max(x,y,label2,Im,Im_seg2,Nb_max,
						Max01,Max1) == 1)
				{
				Nb_max++;
				if(Nb_max > 1)
					{
					Max->Suivant = new t_Pixel;
					assert(Max->Suivant != NULL);
					Max = Max->Suivant;
					}
				Max->num = Nb_max;
				Max->xmax = x;
				Max->ymax = y;
				Max->label = label2;
				Max->label_prec = Max->label;
				Max->Suivant = NULL;
				}
			}
		}
	}
}	


	/*******************************************
	** Init_field: 
	** 		
	**
	********************************************/
void Init_field(Ifloat& Im,Ifloat& Im_seg2,int label2,t_Pixel *Max1)
{
int x,y;
int xmax,ymax;
int Nl = Im.nl();
int Nc = Im.nc();
int Num_max = 0;
float dmin,d;
t_Pixel* Max;

/* Init label des max */

Max = Max1;
do
	{
	Num_max++;
	Max->label = Max->label_prec = (Num_max+label2);
	}
while((Max= Max->Suivant) != NULL);

/* Init domaine */

for(x=0;x<Nc;x++)
	{
	for(y=0;y<Nl;y++)
		{
		if(Im(y,x) == 0.)
			continue;

		Max = Max1;
		dmin = Nl*Nl + Nc*Nc;
		do
			{
			xmax = Max->xmax;
			ymax = Max->ymax;

			d = D(x,y,xmax,ymax);
			if(d<dmin)
				{
				dmin = d;
				Im_seg2(y,x) = float(Max->label);
				} 
			}
		while((Max= Max->Suivant) != NULL);
		}
	}
} 

	/*************************************
	** Cherche_label: 
	** 		
	**
	**************************************/
int Cherche_label(t_Pixel* Max,Ifloat& Im,Ifloat& Im_seg2)
{
int x = Max->xmax;
int y = Max->ymax;
int label_prec = Max->label_prec;
int Nl = Im.nl();
int Nc = Im.nc();

if((x-1)>=0 && Im(y,x-1)>0. && Im_seg2(y,x-1) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y,x-1)));
	return 1;
	}

if((x+1)<Nc && Im(y,x+1)>0. && Im_seg2(y,x+1) != label_prec )
	{
	Max->label = int(floor(Im_seg2(y,x+1)));
	return 1;
	}

if((y-1)>=0 && Im(y-1,x)>0. && Im_seg2(y-1,x) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y-1,x)));
	return 1;
	}

if((y+1)<Nl && Im(y+1,x)>0. && Im_seg2(y+1,x) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y+1,x)));
	return 1;
	}

if((x-1)>=0 && (y-1)>=0 && Im(y-1,x-1)>0. && Im_seg2(y-1,x-1) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y-1,x-1)));
	return 1;
	}

if((x+1)<Nc && (y-1)>=0 && Im(y-1,x+1)>0. && Im_seg2(y-1,x+1) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y-1,x+1)));
	return 1;
	}

if((x-1)>=0 && (y+1)<Nl && Im(y+1,x-1)>0. && Im_seg2(y+1,x-1) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y+1,x-1)));
	return 1;
	}

if((x+1)<Nc && (y+1)<Nl && Im(y+1,x+1)>0. && Im_seg2(y+1,x+1) != label_prec)
	{
	Max->label = int(floor(Im_seg2(y+1,x+1)));
	return 1;
	}

return 0;
}

	/*********************************************
	** Sous_segmente: 
	** 		
	**
	*********************************************/
void Sous_segmente(Ifloat& Im,Ifloat& Im_seg2,t_Pixel* Max1,
			t_Pixel* Max_prec1,t_Pixel* Max01)
{
int x,y;
int xmax,xmax0,ymax,ymax0;
int Nl = Im.nl();
int Nc = Im.nc();
int finit;
float d,dmin,d1,d2;
t_Pixel *Max0;
t_Pixel *Max;
t_Pixel *Max_prec;
t_Pixel *Max_lie=NULL;

/*** DETERMINE LES LABELS DES MAX ***/

Max = Max1;
do
	{
	xmax = Max->xmax;
	ymax = Max->ymax;

	finit = 0;
	Max0 = Max01;
	do
		{
		xmax0 = Max0->xmax;
		ymax0 = Max0->ymax;

		if( (xmax == xmax0) && (ymax == ymax0) )
			{
			finit = 1;
			Max->label = Max0->label;
			}
		}
	while((Max0= Max0->Suivant) != NULL);

	if(finit)
		continue;

	Cherche_label(Max,Im,Im_seg2);
	}
while((Max= Max->Suivant) != NULL);

/*** DETERMINE LES DOMAINES ASSOCIES AUX MAX DONT LE LABEL A CHANGE ***/


Max = Max1;
do
	{
	if(Max->label_prec == Max->label)
			continue;

	xmax = Max->xmax;
	ymax = Max->ymax;

	/* Selectionne le maximum precedent le plus proche de label "label_prec" */
 
	Max_prec = Max_prec1;
	dmin = Nl*Nl + Nc*Nc;
	do
		{
		if(Max_prec->label != Max->label_prec || Max->label_prec == 0)
				continue;
 
		if( Im_seg2(Max_prec->ymax,Max_prec->xmax) != Max->label_prec)
			continue;
		
		d = D(xmax,ymax,Max_prec->xmax,Max_prec->ymax);

		if(d < dmin)
			{
			dmin = d;
			Max_lie = Max_prec;
			}
		}
	while((Max_prec = Max_prec->Suivant) != NULL);


	/* Ajuste les valeurs de Im_seg2 */

	for(x=0;x<Nc;x++)
		{
		for(y=0;y<Nl;y++)
			{
			if(Im(y,x) == 0.)
				continue;

			if( Im_seg2(y,x) != Max->label_prec )
				continue;

			d1 = D(xmax,ymax,x,y);
			d2 = D(Max_lie->xmax,Max_lie->ymax,x,y);

			if(d1 <= d2)
				Im_seg2(y,x) = Max->label;
			}
		}
	}
while((Max= Max->Suivant) != NULL);
}


//----------------------------------------------------------
//	sous_segmente
//
//----------------------------------------------------------
Ifloat *	sous_segmente_iter (Ifloat & Im, int Nb_iter)
{
int n;
int label = 0;
int Nb_max,Nb_max0;
t_Pixel* Max01,*Max0;
t_Pixel* Max1;
t_Pixel* Max_prec1,*Max_prec;

Ifloat* Im_seg2;
Im_seg2 = new Ifloat(Im.nl(),Im.nc(), "ss segmente iter");
Im_seg2->init (0.0);

Max_prec1 = new t_Pixel;
assert(Max_prec1 != NULL);

// Cherche les max du domaine

Max01 = NULL;
Cherche_max(Max01,Max_prec1,Nb_max,Im,*Im_seg2);
Nb_max0 = Nb_max;


// Initialise le domaine

Init_field(Im,*Im_seg2,label,Max_prec1);
		

// Domaine avec un seul max

if(Nb_max ==1)
	return Im_seg2;

// Domaine avec plusieurs max

Max01 = new t_Pixel;
assert(Max01 != NULL);

Max_prec = Max_prec1;
Max0 = Max01;
do
	{
	Max0->num = Max_prec->num;
	Max0->xmax = Max_prec->xmax;
	Max0->ymax = Max_prec->ymax;
	Max0->label = Max_prec->label;
	Max0->label_prec = Max_prec->label_prec;
	Max0->Suivant = NULL;

	Max_prec = Max_prec->Suivant;
	if(Max_prec != NULL)
		{
		Max0->Suivant = new t_Pixel;
		assert(Max0->Suivant != NULL);
		Max0 = Max0->Suivant;
		}
	}
while(Max_prec != NULL);

	// Decoupage iteratif du domaine

for(n=1;n<=Nb_iter;n++)
	{
	Max1 = new t_Pixel;
	assert(Max1 != NULL);
			
	Cherche_max(Max01,Max1,Nb_max,Im,*Im_seg2);
	Sous_segmente(Im,*Im_seg2,Max1,Max_prec1,Max01);
						

	if(Nb_max == Nb_max0)
		break;

	Complete(Max_prec1,Max1);	

	libere(Max1);
	}
		
libere(Max01);
libere(Max_prec1);

return Im_seg2;
}
	
	
	
