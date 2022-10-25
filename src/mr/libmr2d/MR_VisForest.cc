//----------------------------------------------------------
//	Copyright (C) 1996 OCA-NICE
//----------------------------------------------------------
//
//    UNIT
//
//    Version: 
//
//    Author: 	Frederic RUE & Benoit VANDAME
//
//    Date: 	03/05/96
//    
//    File:  	Foret.C
//
//----------------------------------------------------------

//----------------------------------------------------------
//	Description de la structure Foret
//
//
//	Foret est un ensemble d'Arbre. Des fonctions
//	permettent de detecter les objets les compter
//	les sauvegarder, etc ..
//
//	arbre est un pointeurs sur des Arbre,
//	Des que la fonction add_arbre est appelee
//	alors arbre est realloue avec un pointeur en plus
//	Max_arbre indique la capacite maximun du pointeur arbre
//	Nb_arbre indique le nombre effectifs d'arbres
//
//----------------------------------------------------------

#include <math.h>
#include        "IM_Math.h"
#include        "IM_Obj.h"
#include        "MR_Obj.h"
#include        "IM_VisTool.h"
#include	"MR_VisElem.h"
#include	"MR_VisTree.h"

#define MAX_PLAN_WAVELET MAX_SCALE
#define FORET_SIZE_BLOC 1000
//----------------------------------------------------------
//	Constantes pour les noms de fichier
//----------------------------------------------------------
#define		NOM_TEX		1
#define		NOM_TEX_CH	".tex"

#define		NOM_RES		2
#define		NOM_RES_CH	".mes"

#define		NOM_OBJ		3
#define		NOM_OBJ_CH	".gro.fits"

#define		NOM_ARBRE	4
#define		NOM_ARBRE_CH	".gra.fits"

#define		NOM_FITS	5
#define		NOM_FITS_CH	".fits"

#define		NOM_REC		6
#define		NOM_REC_CH	"_rec.fits"

#define		NOM_CRT_O	7
#define		NOM_CRT_O_CH	"_obj.ps"

#define		NOM_CRT_SO	8
#define		NOM_CRT_SO_CH	"_sobj.ps"

#define		NOM_TR_SEG	9
#define		NOM_TR_SEG_CH	"_seg.fits"

#define		NOM_STAR_CAT	10
#define		NOM_STAR_CAT_CH	".cat.fits"

#define		NOM_VIDE	11


//----------------------------------------------------------
//	Constructeur de Forets
//		tout est mis a zero
//----------------------------------------------------------
Foret::Foret (void)
{
	int		i;

	Type_info |= INFO_TYPE_FORET;
	F_type = FORET_ARBRE;
	F_TO = TRANSF_PAVE;
	Type_TO = TO_PAVE_BSPLINE;
	// type de la foret par defaut
	for (i=0 ; i<MAX_PLAN_WAVELET ; i++)
	{
		info[i]		= NULL;
		arbre[i]	= NULL;
		Nb_arbre[i]	= 0;
		Max_arbre[i]	= 0;
	}
}

//----------------------------------------------------------
//	Constructeur de copie
//		Seul les inforamation de taille sont copiees
//		( les arbres ne sont pas recopies depuis F )
//----------------------------------------------------------
 
Foret::Foret (Foret & F)
{
	int	i;

	get_info (F);
	Type_info |= INFO_TYPE_FORET;
	for (i=0 ; i<F_ech ; i++)
	{
		info[i] = new Info ();
		info[i]->get_info (*F.info[i]);
		arbre[i]	= NULL;
		Nb_arbre[i]	= 0;
		Max_arbre[i]	= 0;
	}
	for ( ; i<MAX_PLAN_WAVELET ; i++)
	{
		info[i]		= NULL;
		arbre[i]	= NULL;
		Nb_arbre[i]	= 0;
		Max_arbre[i]	= 0;
	}
}
//----------------------------------------------------------
//	Destructeur de Forets
//		liberation memoire de tout les arbres sauvegardes
//----------------------------------------------------------
Foret::~Foret (void)
{
	free ();
}

//----------------------------------------------------------
//	Destructeur de Forets
//		detruit tout les arbres contenus
//----------------------------------------------------------
void	Foret::free (void)
{
	// recherche des arbres dans arbre puis destruction
	// de arbre[i]
	int		i, j;

	for (j=0 ; j<MAX_PLAN_WAVELET ; j++)
	{
		for (i=0 ; i<Max_arbre[j] ; i++)
			if (arbre[j][i])
				delete arbre[j][i]->ancetre ();
		if (arbre[j])
		{
			delete [] arbre[j];
			arbre[j] = NULL;
		}
	}
	for (i=0 ; i<MAX_PLAN_WAVELET ; i++)
		if (info[i])
		{
			delete info[i];
			info[i] = NULL;
		}
	for (j=0 ; j<MAX_PLAN_WAVELET ; j++)
		Max_arbre[i] = Nb_arbre[i] = 0;

	F_type	= 0;
	F_nl	= 0;
	F_nc	= 0;
	F_TO	= 0;
	F_type	= 0;
	F_nb	= 0; 

}
//----------------------------------------------------------
//
//----------------------------------------------------------
void		Foret::operator = (Foret * F)
{
	int		i;

	free ();
	for (i=0 ; i<MAX_SCALE ; i++)
	{
		info[i] = F->info[i];
		F->info[i] = NULL;
		arbre[i] = F->arbre[i];
		F->arbre[i] = NULL;
		Nb_arbre[i] = F->Nb_arbre[i];
		F->Nb_arbre[i] = 0;
		Max_arbre[i] = F->Max_arbre[i];
		F->Max_arbre[i] = 0;
	}
	get_info (*F);
	delete F;
}

//----------------------------------------------------------
//	Foret :: retire_arbre
//----------------------------------------------------------
 
void		Foret::retire_arbre (Arbre *ar)
{
	int		j, i;

	j = ar->Num_ech;
	i = ar->Num_ind;
	assert (Nb_arbre[j] >= i);
	assert (arbre[j][i] == ar);

	arbre[j][i] = NULL;
	for (i++ ; i<Nb_arbre[j] ; i++)
		arbre[j][i-1] = arbre[j][i];
	arbre[j][Nb_arbre[j]-1] = NULL;
	Nb_arbre[j]--;
	delete ar;
}

//----------------------------------------------------------
//	Foret :: retire_arbre
//----------------------------------------------------------
void		Foret::retire_arbre (int j, int i)
{
	Arbre		*ar;

	assert (Nb_arbre[j] > i);
	ar = arbre[j][i];
	arbre[j][i] = NULL;
	for (i++ ; i<Nb_arbre[j] ; i++)
		arbre[j][i-1] = arbre[j][i];
	arbre[j][Nb_arbre[j]-1] = NULL;
	Nb_arbre[j]--;
	delete ar->ancetre ();
}

//----------------------------------------------------------
//	Foret :: add_arbre_fils
//----------------------------------------------------------
void		Foret::add_arbre_fils (Arbre *ar)
{
	int		i, x, y;
	Arbre		*pere;

	assert (F_type == FORET_ARBRE);

	// cas ou ar est en haut des echelles
	// il s'agit alors d'un ancetre
	if (ar->Num_ech+1 == F_ech)
	{
		add_arbre (ar, ar->Num_ech);
		return;
	}

	// recherche d'un pere probable
	// a l'echelle au dessus
	ar->get_pos_max (x, y);
	
	if ((F_TO == (int) TRANSF_PYR) ||
			   ((F_TO == (int) TRANSF_SEMIPYR) && (ar->Num_ech > 1)))
 	{
// 		x = pos_2_pave (x, ar->Num_ech);
// 		y = pos_2_pave (y, ar->Num_ech);
// 		x = pos_2_pyr (x, ar->Num_ech+1);
// 		y = pos_2_pyr (y, ar->Num_ech+1);
		x /= 2;
		y /= 2;
	}

	for (i=0 ; i<Nb_arbre[ar->Num_ech+1] ; i++)
	{
		pere = get_arbre (ar->Num_ech+1, i);
		if (pere && pere->is_point (x, y))
		{
			pere->add_fils (ar);
			return;
		}
	}

	// pas de pere repertorie
	// ar est un ancetre
	add_arbre (ar, ar->Num_ech);
}
 

//----------------------------------------------------------
//	Foret :: add_arbre
//		ajoute un arbre dans la liste des arbres
//		si celui-ci n'est pas deja repertoire
//----------------------------------------------------------
void		Foret::add_arbre (Arbre *ar, int ech)
{
	int	i;
	
	assert (ech < F_ech);
	// recherche si l'arbre existe deja
//	for (i=0 ; i<Nb_arbre[ech] ; i++)
//		if (arbre[ech][i] == ar)
//			return;

	// ecriture
	if (Nb_arbre[ech] >= Max_arbre[ech])
	{
		Max_arbre[ech] += FORET_SIZE_BLOC;
		if (arbre[ech] == NULL) 
		     arbre[ech] = (Arbre **) malloc (Max_arbre[ech] * sizeof (Arbre *));
		else arbre[ech] = (Arbre **)realloc (arbre[ech],
					Max_arbre[ech] * sizeof (Arbre *));
		assert (arbre[ech]);
		for (i=Nb_arbre[ech] ; i<Max_arbre[ech] ; i++)
			arbre[ech][i] = NULL;
	}
	arbre[ech][Nb_arbre[ech]] = ar;
 	Nb_arbre[ech]++;
}


//----------------------------------------------------------
//	Foret :: creat_wavelet
//		creation d'une Wavelet minimum autour de l'arbre
//		numero a
//      This routine is copy of the next one, but W0 is a parameter
//      instead of being return by the routine.
//      This has been done for compatibility with the other MVM model,
//      when ListObj class call the creat_wavelet rotoutine
//----------------------------------------------------------
void Foret::creat_wavelet (MultiResol & W, 
                           int j, int i, int &dep_x, int &dep_y, 
                           int AddBorderX, int AddBorderY)
{
 	Arbre		*ar, *pere;
	int		xm, ym, xM, yM,
			dx, dy, dxi, dyi,
			x, y, n;

	assert (j<F_ech);
	assert (i < Nb_arbre[j]);
//	assert (F_type == FORET_OBJ);


	xm = ym = MAX (F_nl, F_nc);
	xM = yM = 0;

	// option opt_sobj
	// si opt_sobj == 0 alors tous les sous-objets sont concernes
	// sinon seul l'objet sans sous-objet et reconstruit
	pere = arbre[j][i];
	for (ar=pere ; ar->Proj ; ar=ar->fils[0]);

	// recherche du cadre minimun qui autour de arbre[a]
	// et ses fils
	pere->taille (*this, xm, ym, xM, yM);

	// le cadre est choisi au double du minimum
	// et sur un multiple de 2
	dxi = xM - xm + 1;
	dyi = yM - ym + 1;
	dx = (dxi << 1) + (1<<(pere->Num_ech+3));
	dy = (dyi << 1) + (1<<(pere->Num_ech+3));

	// Code pour que l'image soit plus grande que le masque a trou 5 * 5
	dx = MAX (dx, 5*(1<<pere->Num_ech));
	dy = MAX (dy, 5*(1<<pere->Num_ech));
	dx += AddBorderX;
	dy += AddBorderY;
	// la Wavelet ne doit pas etre plus grande que l'original
	dx = MIN (1*dx, F_nc);
	dy = MIN (1*dy, F_nl);

	// Code pour que l'image soit une puissance de 2 > dx et dy
	// pour les decimations correcte
	if (F_TO == (int) TRANSF_PYR)
  	{ 
 		n = 1 << pere->Num_ech;   // pos_2_pave(1, pere->Num_ech);
		dx--;
		dx = ((dx/n)+1) * n;
		dy--;
		dy = ((dy/n)+1) * n;
	}
        else if ((F_TO == (int) TRANSF_SEMIPYR) && (pere->Num_ech > 1))
        {
                n = 1 << (pere->Num_ech-1);   // pos_2_pave(1, pere->Num_ech);
		dx--;
		dx = ((dx/n)+1) * n;
		dy--;
		dy = ((dy/n)+1) * n;
        }
        
	// recherche du vecteur de translation 
	// entre les images des Elements, et de la Wavelet
	x = xm + dxi/2 - dx/2;
	y = ym + dyi/2 - dy/2;

	// si pyramide x et y doivent etre multiple de n
	if (F_TO == (int) TRANSF_PYR)
	{
		n = 1 << pere->Num_ech; // n = pos_2_pave(1, pere->Num_ech);
		x = iround (x/(double)n);
		x *= n;
		y =  iround (y/(double)n);
		y *= n;
	}
	else if ((F_TO == (int) TRANSF_SEMIPYR) && (pere->Num_ech > 1))
        {
                n = 1 << (pere->Num_ech-1);   // pos_2_pave(1, pere->Num_ech);
		x = iround (x/(double)n);
		x *= n;
		y = iround (y/(double)n);
		y *= n;
	}
	
        W.alloc(dy, dx, pere->Num_ech+2, (type_transform) Type_TO, "creat wavelet");

	for (n=0 ; n<=pere->Num_ech ; n++)
	{
 		W.band(n).init(0.0);
  		// tout les elements sont translate de x, y
		pere->creat_wavelet (W, x, y, n);
 	}
	// mise a zero du dernier plan lisse
	n = pere->Num_ech+1;
	W.band(n).init(0.0);
	dep_x = x;
	dep_y = y;
}
//----------------------------------------------------------
//	Foret :: creat_wavelet
//		creation d'une Wavelet minimum autour de l'arbre
//		numero a
//
//----------------------------------------------------------

MultiResol *	Foret::creat_wavelet (int j, int i, int &dep_x, int &dep_y, 
                                      int AddBorderX, int AddBorderY)
{
	MultiResol		*W;
	Arbre		*ar, *pere;
	int		xm, ym, xM, yM,
			dx, dy, dxi, dyi,
			x, y, n;

	assert (j<F_ech);
	assert (i < Nb_arbre[j]);
//	assert (F_type == FORET_OBJ);


	xm = ym = MAX (F_nl, F_nc);
	xM = yM = 0;

	// option opt_sobj
	// si opt_sobj == 0 alors tous les sous-objets sont concernes
	// sinon seul l'objet sans sous-objet et reconstruit
	pere = arbre[j][i];
	for (ar=pere ; ar->Proj ; ar=ar->fils[0]);

	// recherche du cadre minimun qui autour de arbre[a]
	// et ses fils
	pere->taille (*this, xm, ym, xM, yM);

	// le cadre est choisi au double du minimum
	// et sur un multiple de 2
	dxi = xM - xm + 1;
	dyi = yM - ym + 1;
	dx = (dxi << 1) + (1<<(pere->Num_ech+3));
	dy = (dyi << 1) + (1<<(pere->Num_ech+3));

	// Code pour que l'image soit plus grande que le masque a trou 5 * 5
	dx = MAX (dx, 5*(1<<pere->Num_ech));
	dy = MAX (dy, 5*(1<<pere->Num_ech));
	dx += AddBorderX;
	dy += AddBorderY;
	// la Wavelet ne doit pas etre plus grande que l'original
	dx = MIN (1*dx, F_nc);
	dy = MIN (1*dy, F_nl);

	// Code pour que l'image soit une puissance de 2 > dx et dy
	// pour les decimations correcte
	if (F_TO == (int) TRANSF_PYR)
  	{ 
 		n = 1 << pere->Num_ech;   // pos_2_pave(1, pere->Num_ech);
		dx--;
		dx = ((dx/n)+1) * n;
		dy--;
		dy = ((dy/n)+1) * n;
	}
        else if ((F_TO == (int) TRANSF_SEMIPYR) && (pere->Num_ech > 1))
        {
                n = 1 << (pere->Num_ech-1);   // pos_2_pave(1, pere->Num_ech);
		dx--;
		dx = ((dx/n)+1) * n;
		dy--;
		dy = ((dy/n)+1) * n;
        }
        
	// recherche du vecteur de translation 
	// entre les images des Elements, et de la Wavelet
	x = xm + dxi/2 - dx/2;
	y = ym + dyi/2 - dy/2;

	// si pyramide x et y doivent etre multiple de n
	if (F_TO == (int) TRANSF_PYR)
	{
		n = 1 << pere->Num_ech; // n = pos_2_pave(1, pere->Num_ech);
		x = iround (x/(double)n);
		x *= n;
		y =  iround (y/(double)n);
		y *= n;
	}
	else if ((F_TO == (int) TRANSF_SEMIPYR) && (pere->Num_ech > 1))
        {
                n = 1 << (pere->Num_ech-1);   // pos_2_pave(1, pere->Num_ech);
		x = iround (x/(double)n);
		x *= n;
		y = iround (y/(double)n);
		y *= n;
	}
	
        W = new MultiResol (dy, dx, pere->Num_ech+2, (type_transform) Type_TO, "creat wavelet");

	for (n=0 ; n<=pere->Num_ech ; n++)
	{
 		W->band(n).init(0.0);
  		// tout les elements sont translate de x, y
		pere->creat_wavelet (*W, x, y, n);
 	}

	// mise a zero du dernier plan lisse
	n = pere->Num_ech+1;
	W->band(n).init(0.0);
	dep_x = x;
	dep_y = y;
	return W;
}

 

//----------------------------------------------------------
//	Foret :: write_wavelet
//		ecrit tous les Arbres dans W
//		en fonction de leur niveau
//----------------------------------------------------------
void	Foret::write_wavelet (MultiResol & W)
{
	int		j, i, n;
	for (j=0 ; j<F_ech ; j++)
	{
		for (i=0 ; i<Nb_arbre[j] ; i++)
			for (n=0 ; n<=j ; n++)
				arbre[j][i]->write_wavelet (W, n);
 	}
 }

//----------------------------------------------------------
//	Foret :: write_wavelet_seg
//		ecrit tous les Arbres dans W (Valeur de Num_elem)
//		en fonction de leur niveau
//----------------------------------------------------------
void	Foret::write_wavelet_seg (MultiResol & W)
{
	int		i, j, n;
	for (n=F_ech, n-- ; n>=0 ; n--)
		for (j=0 ; j<F_ech ; j++)
		{
			for (i=0 ; i<Nb_arbre[j] ; i++)
					arbre[j][i]->write_wavelet_seg (W, n);
 		}
}

//----------------------------------------------------------
//	Foret :: write_wavelet
//		ecrit l'arbre numero p (et ses fils)
//		dans la Wavelet W 
//----------------------------------------------------------
void	Foret::write_wavelet (MultiResol & W, int j, int i)
{

	assert (j < F_ech);
	assert (i < Nb_arbre[j]);

	for (int n=0 ; n<F_ech ; n++)
	{
		arbre[j][i]->write_wavelet (W, n);
 	}	
	
}

//----------------------------------------------------------
//	Foret :: read_foret
//		lecture du fichier file
//----------------------------------------------------------
void	Foret::read_foret_info (void)
{
	char		c[CHAINE], com[CHAINE];
	int		j;
	long		l;

	assert (fptr);

	for (j=0 ; j<F_ech ; j++)
	{
		info[j] = new Info ();
		sprintf (c, "P%dI_MOY", j);
		fits_read_key_dbl (fptr, c, &info[j]->Moyenne, com, &fptr_st);
		sprintf (c, "P%dI_SIG", j);
		fits_read_key_dbl (fptr, c, &info[j]->Sigma, com, &fptr_st);
		sprintf (c, "P%dI_FOND", j);
		fits_read_key_dbl (fptr, c, &info[j]->Fond, com, &fptr_st);
		sprintf (c, "P%dI_BRT", j);
		fits_read_key_dbl (fptr, c, &info[j]->Bruit, com, &fptr_st);
		sprintf (c, "P%dI_SE", j);
		fits_read_key_dbl (fptr, c, &info[j]->Seuil, com, &fptr_st);
		sprintf (c, "P%dI_NBST", j);
		fits_read_key_dbl (fptr, c, &info[j]->Nbpoint, com, &fptr_st);

		sprintf (c, "P%dE_NB", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_elem = (int)l;
		sprintf (c, "P%dE_DEP", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_elem_dep = (int)l;
		sprintf (c, "P%dE_SS", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_elem_ss = (int)l;
		sprintf (c, "P%dE_NET", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_elem_net = (int)l;

		sprintf (c, "P%dO_NB", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_obj = (int)l;
		sprintf (c, "P%dO_S_NB", j);
		fits_read_key_lng (fptr, c, &l, com, &fptr_st);
		info[j]->Nb_s_obj = (int)l;
	}
	erreur_fits ();
}


//----------------------------------------------------------
//	Foret :: write_foret
//		ecrit la foret dans le fichier file
//----------------------------------------------------------
void	Foret::write_foret_info (void)
{
	char		c[CHAINE];
	int		j;

	assert (fptr);

	for (j=0 ; j<F_ech ; j++)
	{
		sprintf (c, "P%dI_MOY", j);
		fits_write_key_dbl (fptr, c, info[j]->Moyenne, 14, (char*)"", &fptr_st);
		sprintf (c, "P%dI_SIG", j);
		fits_write_key_dbl (fptr, c, info[j]->Sigma, 14, (char*)"", &fptr_st);
		sprintf (c, "P%dI_FOND", j);
		fits_write_key_dbl (fptr, c, info[j]->Fond, 14, (char*)"", &fptr_st);
		sprintf (c, "P%dI_BRT", j);
		fits_write_key_dbl (fptr, c, info[j]->Bruit, 14, (char*)"", &fptr_st);
		sprintf (c, "P%dI_SE", j);
		fits_write_key_dbl (fptr, c, info[j]->Seuil, 14, (char*)"", &fptr_st);
		sprintf (c, "P%dI_NBST", j);
		fits_write_key_dbl (fptr, c, info[j]->Nbpoint, 14, (char*)"", &fptr_st);

		sprintf (c, "P%dE_NB", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_elem, (char*)"", &fptr_st);
		sprintf (c, "P%dE_DEP", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_elem_dep, (char*)"", &fptr_st);
		sprintf (c, "P%dE_SS", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_elem_ss, (char*)"", &fptr_st);
		sprintf (c, "P%dE_NET", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_elem_net, (char*)"", &fptr_st);

		sprintf (c, "P%dO_NB", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_obj, (char*)"", &fptr_st);
		sprintf (c, "P%dO_S_NB", j);
		fits_write_key_lng (fptr, c, (long)info[j]->Nb_s_obj, (char*)"", &fptr_st);
	}
	erreur_fits (); 
}



//----------------------------------------------------------
//	Foret :: read_foret
//		lecture du fichier file
//----------------------------------------------------------
//----------------------------------------------------------
//		get_nom
//----------------------------------------------------------
void		get_nom (char * S, char * dest, int type)
{
	char		*pc;
	char		src[CHAINE];


	strcpy (src, S);
	pc = strstr (src, ".fits\0");
	if (pc)		*pc = 0;
	pc = strstr (src, ".fts\0");
	if (pc)		*pc = 0;

	switch (type)
	{
		case NOM_TEX	:	sprintf (dest, "%s%s", src, NOM_TEX_CH);
					break;
		case NOM_OBJ	:	sprintf (dest, "%s%s", src, NOM_OBJ_CH);
					break;
		case NOM_ARBRE	:	sprintf (dest, "%s%s", src, NOM_ARBRE_CH);
					break;
		case NOM_RES	:	sprintf (dest, "%s%s", src, NOM_RES_CH);
					break;
		case NOM_FITS	:	sprintf (dest, "%s%s", src, NOM_FITS_CH);
					break;
		case NOM_REC	:	sprintf (dest, "%s%s", src, NOM_REC_CH);
					break;
		case NOM_CRT_O	:	sprintf (dest, "%s%s", src, NOM_CRT_O_CH);
					break;
		case NOM_CRT_SO	:	sprintf (dest, "%s%s", src, NOM_CRT_SO_CH);
					break;
		case NOM_STAR_CAT:	sprintf (dest, "%s%s", src, NOM_STAR_CAT_CH);
					break;
		case NOM_TR_SEG	:	sprintf (dest, "%s%s", src, NOM_TR_SEG_CH);
					break;
		default		:	strcpy (dest, src);
					break;
	}
}

void	Foret::read_foret (char * name, int type)
{
	Arbre		*a, *t;
	int		n, i, p;
	char		file[CHAINE];

	// ouverture du fichier fits
	switch (type)
	{
		case FORET_ARBRE	:	type = NOM_ARBRE;
						break;
		default			:	type = NOM_OBJ;
						break;
	}
	get_nom (name, file, type);
	open_image_read (file);
	assert (Type_info & INFO_TYPE_FORET);
	Element::Stop = F_stop;
	// Wavelet::init_type_TO (F_TO);

	// lecture des informations sur l'image
	read_foret_info ();

	//lecture des elements
	for (i=0, n=1 ; i<F_nb ; i++)
	{
		// creation d'un arbre et lecture de celui-ci
		a = new Arbre ();
		a->read_arbre (*this, n);

		// l'arbre lu est ajoute
		t = a;
		for (p=0 ; t->Proj ; p++)
			t = t->fils[0];
		add_arbre (a, t->Num_ech);
	}
	close_image ();
}

//----------------------------------------------------------
//	Foret :: write_foret
//		ecrit la foret dans le fichier file
//----------------------------------------------------------
void	Foret::write_foret (char * name)
{
	char		file[CHAINE];
	int		i, n, j;
	long		naxes[1];

	// nom du fichier fits
	if (F_type == FORET_ARBRE)
		get_nom (name, file, NOM_ARBRE);
	else
		get_nom (name, file, NOM_OBJ);

	// creation du fichier fits
	open_image_write (file);
	// taille du fichier inconnu
	naxes[0] = 1;
	fits_create_img (fptr, -32, 1, naxes, &fptr_st);

	// calcul le nombre d'elements
	for (i=0, F_nb=0 ; i<F_ech ; i++)
		F_nb += Nb_arbre[i];

	// point d'arret de element
	F_stop = Element::Stop;
	write_image_info ();
	erreur_fits ();

	// ecriture des infos
	write_foret_info ();

	// ecriture des elements (bitmap)
	for (j=0, n=1 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			arbre[j][i]->write_arbre (*this, n);

	// actualisation de la taille du fichier fits
	fits_modify_key_lng (fptr, "NAXIS1", n, "&", &fptr_st);

	// fermeture du fichier
	close_image ();
}

//----------------------------------------------------------
//	Foret :: retire_isole (void)
//		retire les elements isoles (ni pere ni fils)
//		dont le plans est >= dep
//		defaut : dep = 0
//----------------------------------------------------------
void	Foret::retire_isole (int dep)
{
	int		i, j;
	assert (F_type == FORET_ARBRE);

	for (j=dep ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			if (!arbre[j][i]->pere && !arbre[j][i]->Nb_fils)
			{
				delete arbre[j][i];
				arbre[j][i] = NULL;
				info[j]->Nb_elem --;
				info[j]->Nb_elem_net ++;
			}

	// reorganise le tableau des arbres
	// il peut y avoir des arbre NULL a retirer
	for (j=dep ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			if (!arbre[j][i])
			{
				// si arbre[i] est NULL
				// il est remplace par le dernier de la liste
				// puis la boucle reprend un increment avant
				arbre[j][i] = arbre[j][Nb_arbre[j]-1];
				arbre[j][Nb_arbre[j]-1] = NULL;
				Nb_arbre[j] --;
				i--;
			}

	// numerote les arbres car il ne sont plus dans l'ordre
	numerote_arbre ();
}

//----------------------------------------------------------
//	Foret :: retire_isole (void)
//		retire les elements isoles (ni pere ni fils)
//		dont le plans est > dep
//		defaut : dep = 0
//----------------------------------------------------------
void	Foret::retire_petit (int taille)
{
	int		i, j;
	assert (F_type == FORET_ARBRE);

	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			if (arbre[j][i]->Nb_pixel <= taille)
			{
				delete arbre[j][i];
				arbre[j][i] = NULL;
				info[j]->Nb_elem --;
				info[j]->Nb_elem_net ++;
			}

	// reorganise le tableau des arbres
	// il peut y avoir des arbre NULL a retirer
	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			if (!arbre[j][i])
			{
				// si arbre[i] est NULL
				// il est remplace par le dernier de la liste
				// puis la boucle reprend un increment avant
				arbre[j][i] = arbre[j][Nb_arbre[j]-1];
				arbre[j][Nb_arbre[j]-1] = NULL;
				Nb_arbre[j] --;
				i--;
			}

	// numerote les arbres car les numero ne sont plus consecutifs
	numerote_arbre ();
}


//----------------------------------------------------------
//	Foret :: identif_obj (void)
//		recherche les objets dans tous les arbres
//		les objets detectes sont alors numerotes
//----------------------------------------------------------

void	Foret::identif_obj (float e, float dmax)
{
	int		i, j;

	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			arbre[j][i]->identif_obj(0, (int)e, *this, dmax);

	numerote_obj ();
}


//----------------------------------------------------------
//	Foret::is_s_obj
//----------------------------------------------------------
int	Foret::is_s_obj (int j, int i)
{
	assert (F_type == FORET_OBJ);
	return arbre[j][i]->S_obj;
}

//----------------------------------------------------------
//	Foret::is_s_obj
//----------------------------------------------------------
int	Foret::is_s_obj (int j, int i, int max)
{
	Arbre		*sobj;

	assert (F_type == FORET_OBJ);

	for (sobj=arbre[j][i] ; sobj->Proj ; sobj=sobj->fils[0])
	if (!sobj->S_obj)
		return 0;

	Arbre		*ar_p;
	ar_p = get_pere_obj (arbre[j][i]);
	if (ar_p->Num_ech < max)
		return arbre[j][i]->S_obj;
	return 0;
}


//----------------------------------------------------------
//	Foret :: numerote_obj (void)
//		Nemerote les objets en commencant par 1
//		a chaque niveau
//		seul le champ Num_obj des arbres est modifie
//----------------------------------------------------------
void	Foret::numerote_obj (void)
{
	int		i, j, numero[MAX_SCALE];
	Arbre		*a;

	// numerote les objets detectes dans une
	// foret de type arbre uniquement
	assert (F_type == FORET_ARBRE);

	// mise a zero du tableau numero
	for (i=0 ; i<MAX_PLAN_WAVELET ; i++)
		numero[i] = 0;

	// parcour de toutes les echelles
	for (j=F_ech, j-- ; j>=0 ; j--)
		for (i=0 ; i<Nb_arbre[j] ; i++)
		{
			a = arbre[j][i]->ancetre ();
			a->numerote_obj (numero);
		}

	// mise a jour de Nb_elem
	for (i=0 ; i<F_ech ; i++)
		info[i]->Nb_obj = numero[i]; 
}

//----------------------------------------------------------
//	Foret :: numerote_arbre (void)
//		Nemerote les arbres en commencant par 1
//		a chaque niveau
//		seul le champ Num_ind dess arbres est modifie
//----------------------------------------------------------
void	Foret::numerote_arbre (void)
{
	int		i, j, numero[MAX_SCALE];
	Arbre		*a;

	// mise a zero du tableau numero
	for (i=0 ; i<MAX_PLAN_WAVELET ; i++)
		numero[i] = 0;

	// parcour de toutes les echelles
	for (j=F_ech, j-- ; j>=0 ; j--)
		for (i=0 ; i<Nb_arbre[j] ; i++)
		{
			a = arbre[j][i]->ancetre ();
			a->numerote_arbre (numero);
		}

	// mise a jour de Nb_elem
	for (i=0 ; i<F_ech ; i++)
		info[i]->Nb_elem = numero[i];
}


 //----------------------------------------------------------
//	Foret :: get_nb_obj
//		retourne le nombre d'objets
//----------------------------------------------------------
int	Foret::get_nb_obj (void)
{
	int	rep = 0;
	int	i;
	for (i=0 ; i<F_ech ; i++)
		rep += info[i]->Nb_obj;
	return rep;
}

//----------------------------------------------------------
//	Foret :: get_nb_obj
//		retourne le nombre d'objets
//----------------------------------------------------------
int	Foret::get_nb_obj (int p)
{
	return Nb_arbre[p];
}

//----------------------------------------------------------
//	Foret :: get_nb_s_obj
//		retourne le nombre d'objets dans le plan ech,
//		dans tout les arbres et les fils
//----------------------------------------------------------
int	Foret::get_nb_s_obj (int ech)
{
	return info[ech]->Nb_s_obj;
}

//----------------------------------------------------------
//	Foret :: get_nb_s_obj
//		retourne le nombre d'objets
//----------------------------------------------------------
int	Foret::get_nb_s_obj (void)
{
	int	rep = 0, i;
	for (i=0 ; i<F_ech ; i++)
		rep += info[i]->Nb_s_obj;
	return rep;
}



//----------------------------------------------------------
//	Foret :: creat_obj (void)
//		retourne une Foret, ou chaque arbre est un
//		objet avec, s'il existe : son pere.
//	retour : pointeur sur une Foret
//----------------------------------------------------------
Foret *	Foret::creat_obj (void)
{
 	int		i, j;
	Foret		*rep;

	// creation d'une nouvelle foret
	rep = new Foret (*this);
	rep->F_type = FORET_OBJ;

	// ajout de tous les objets
	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			arbre[j][i]->creat_obj (*rep);

	// numerotation des arbres
	rep->numerote_arbre ();

	return rep;
}

//----------------------------------------------------------
//	Foret :: get_obj
//----------------------------------------------------------
Arbre *		Foret::get_obj (int p, int n, Arbre * a)
{
	Arbre		*rep;
	int		i, j;

//	assert (F_type == FORET_OBJ);
//	assert (!a);

	if (F_type == FORET_OBJ && !a)
	{
		// pas d'bojet au niveau p
		if (n >= Nb_arbre[p])
			return NULL;
		for (rep=arbre[p][n] ; rep->Proj ; rep=rep->fils[0]);
			return rep;
	}
	else
	{
		for (j=0 ; j<F_ech ; j++)
			for (i=0 ; i<Nb_arbre[j] ; i++)
			{
				rep = arbre[j][i]->get_obj (p, n);
				if (rep && rep!=a)
					return rep;
			}
	}
	return NULL;
}

 
//----------------------------------------------------------
//	Foret :: get_arbre
//		retorune l'arbre du plan p ayant
//		le numero d'element n
//----------------------------------------------------------
Arbre *		Foret::get_arbre (int p, int n, Arbre * a)
{
	Arbre		*rep;
	int		i, j;
	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
		{
			rep = arbre[j][i]->get_arbre (p, n);
			if (rep && rep != a)
				return rep;
		}
	return NULL;
}

//----------------------------------------------------------
//	Foret :: get_pere_obj
//		retorune l'arbre du plan p ayant
//		le numero d'element n
//----------------------------------------------------------
Arbre *		Foret::get_pere_obj (Arbre * a)
{
	if (!a->S_obj)
		return NULL;

	Arbre	*rep;
	rep = get_obj (a->Num_ech, a->Num_obj, a);
	if (rep && rep->pere)
		return rep->pere->get_pere_obj ();
	else
		return NULL;
}

//----------------------------------------------------------
//	Foret :: sous_segmnete
//----------------------------------------------------------

void		Foret::sous_segmente (int dep, int min_pixel)
{
	assert (F_type == FORET_ARBRE);
	int		i, j,x,y;
	Ifloat		*I, *Iss;

	for (j=dep ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
		{
			if (arbre[j][i]->Nb_pixel <= min_pixel)
				continue;
			I = arbre[j][i]->creat_image (x,y);
			Iss = sous_segmente_iter (*I, 20);
			arbre[j][i]->Num_ech = j;
			arbre[j][i]->Num_ind = i;
			arbre[j][i]->change_arbre (*I, *Iss, *this);
			delete Iss;
			delete I;
		}
	numerote_arbre ();
}

 

//----------------------------------------------------------
//	Foret :: creation_foret
//		creation d'une foret a partir de Wavelet
//		Wi est la Wavelet initiale
//		Ws est la Wavelet seuillee
//----------------------------------------------------------
void	Foret::creation_foret ( MultiResol & Wi)
{
	int		i, j;
	Arbre***	init = NULL;

	free ();
	F_TO = Wi.Set_Transform;
	F_nl = Wi.size_ima_nl();
	F_nc = Wi.size_ima_nc();
	F_ech =  Wi.nbr_band() - 1;
	F_type = FORET_ARBRE;
 	Type_TO = Wi.Type_Transform;
 	
         if ((Wi.Set_Transform == TRANSF_PYR)
             ||  (Wi.Set_Transform == TRANSF_SEMIPYR)) F_wave = 1;
        else F_wave = 0;

// cout << "Nb_ech = " <<  F_ech << endl;
        for (i=0 ; i< F_ech ; i++) info[i] = new Info();
                
 	// construction des arbres avec les relations inter-echelles
	// init est un poiteur sur les differents plans
	init = construit_arbre (Wi, *this);

	for (i=0 ; i< F_ech; i++)
		info[i]->Nb_elem = info[i]->Nb_elem_dep;
 		
	// compte le nombre d'elements dans le plan
// 	for (i=0 ; i<F_ech ; i++)
// 	{
// 		info[i] = new Info ();
// 		info[i]->get_info (*Wi.plan[i]);
// 		info[i]->Nb_elem = info[i]->Nb_elem_dep;
// 	}

	// dans la foret seul les ancetres sont sauvegardes
	for (i=0 ; i<F_ech ; i++)
		for (j=0 ; j<info[i]->Nb_elem ; j++)
		{
			if (!init[i][j]->pere)
				add_arbre (init[i][j], init[i][j]->Num_ech);
		}

	// liberation memoire des pointeurs de pointeurs
	for (i=0 ; i<F_ech ; i++)
		if (init[i])
			delete [] init[i];
	delete [] init;

	// numerotation des arbres en fonction des niveaux
	numerote_arbre ();
}

   	
 // 
//  
// 	int		i, j;
// 	Type_TO = Wi.Type_Transform;
//         Set_Trans = Wi.Set_Transform;
//         if (Wi.Set_Transform == TRANSF_PYR) Type_Wave = 1;
//         else Type_Wave = 0;
// 
//  	Nl_init = Wi.size_ima_nl();
// 	Nc_init = Wi.size_ima_nc();
// 	Nb_ech = Wi.nbr_scale() - 1;
// 	Arbre***	init;
// 	for (i=0 ; i<Nb_ech ; i++)
// 	{
// 		info[i] = new Info ();
// 	}
// 
//  	init = construit_arbre (Wi, *this);
//  
//  	for (i=0 ; i<Nb_ech ; i++)
// 		info[i]->Nb_elem = info[i]->Nb_elem_dep;
//  
//  	for (i=0 ; i<Nb_ech ; i++)
// 		for (j=0 ; j<info[i]->Nb_elem ; j++)
// 			add_arbre (init[i][j]->ancetre ());
//  
//  	for (i=0 ; i<Nb_ech ; i++)
// 		if (init[i])
// 			delete [] init[i];
// 	delete [] init;
 
 
/**************************************************************************/

void mr_make_graph(MultiResol & W, Foret & F, Bool SousSegment, Bool KillIsolObj)
{
   int FirstScale_to_Subsegment = 3;
   int MinSize_for_Subsegment = 9;
   
// cout << "MR graph 1 " << endl;
    F.creation_foret (W);
// cout << "MR graph 2 " << endl;

    // Creation graphe de connection inter-echelle
    if (KillIsolObj == True) F.retire_isole(-1); // elimine les elements isoles
// cout << "MR graph 3 " << endl;

 
// cout << "MR graph 4 " << endl;

    // F.numerote_arbre ();
    // sous-segmentation des domaines
    if (SousSegment == True) F.sous_segmente(FirstScale_to_Subsegment, 
                                              MinSize_for_Subsegment); 
}

/**************************************************************************/

void mr_make_obj_tree(Foret & F, Foret * &F_obj, Bool UseEnerg, Bool Align, float Dmax)
{
     
   // Identification des objets
   if (UseEnerg == True) F.identif_obj (1., Dmax);
   else F.identif_obj (0., Dmax);

   //creation de la foret des objets
   if (Align == False) F_obj = F.creat_obj ();
   else F_obj = F.creat_obj_aligne (Dmax);
}

/**************************************************************************/

void mr_make_obj_tree(MultiResol & W, Foret * &F_obj, Bool SousSegment, 
                  Bool KillIsolObj, Bool UseEnerg, Bool Align, float Dmax)
{
    Foret F;

    mr_make_graph(W, F, SousSegment, KillIsolObj);
    mr_make_obj_tree(F, F_obj, UseEnerg,Align,Dmax);
}

/**************************************************************************/


//----------------------------------------------------------
//	Foret::creat_obj_aligne
//		creation d'un graphe de type OBJET
//		avec recherche des alignements
//----------------------------------------------------------
Foret *		Foret::creat_obj_aligne (float dmax)
{
	int		i, j;
	Foret		*rep;

	// creation d'une nouvelle foret
	rep = new Foret (*this);
	rep->F_type = FORET_OBJ;

	// ajout de tous les objets
	for (j=0 ; j<F_ech ; j++)
		for (i=0 ; i<Nb_arbre[j] ; i++)
			arbre[j][i]->creat_obj_aligne (*rep, dmax);

	// numerotation des arbres
	rep->numerote_arbre ();

	return rep;
}

/**************************************************************************/
