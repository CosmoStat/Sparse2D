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
//    Date: 	05/04/96
//    
//    File:  	Arbre.C
//
//----------------------------------------------------------
//
//
//
//	Fonctions definies dans ce fichier:
//
//
//	Arbre::Arbre (void)
//	Arbre::Arbre (Element & A)
//	Arbre::Arbre (Arbre & A)
//	Arbre::~Arbre (void)
//	void	Arbre::free (void)
//	void	Arbre::add_fils (Arbre *f)
//	void	Arbre::retire_fils (Arbre *f)
//	Arbre *		Arbre::ancetre (void)
//	static void	init_arbre_plan (Wavelet & Wi, Wavelet & Ws,
//						Arbre *** rep, int i)
//	Arbre***	construit_arbre (Wavelet & W)
//	void	Arbre::write_arbre (FILE * f)
//	void	Arbre::read_arbre (FILE * f)
//	void	Arbre::creat_wavelet (Wavelet & W, int px, int py, int p)
//	void	Arbre::write_wavelet (Wavelet & W, int p)
//	Arbre *		Arbre::retire_p_e (void)
//	void		Arbre::identif_obj (void)
//	void	Arbre::numerote_obj (int ech, int & val)
//	void	Arbre::get_obj (Foret & F)
//	void	Arbre::taille (int & xmin, int & ymin, int & xmax, int & ymax)
//	void	Arbre::write_arb (FILE * f, Foret & g)
//	void	Arbre::write_obj (FILE * f, Foret & g)
//
//
//----------------------------------------------------------

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_VisTool.h"
#include "MR_VisElem.h"
#include "MR_VisTree.h"
#include "IM_Graphics.h"

#define PRINT 0
#define WRITE_SEGMENTE 0
#define DEBUG_VIS 0

Bool VisionCleanBord = False;
Bool KeepAll = False;


//----------------------------------------------------------
//	Variables globales
//----------------------------------------------------------

//----------------------------------------------------------
//	Constructeur de Arbre
//		constructeur par defaut
//		toutes les donnes sont mise a nul
//----------------------------------------------------------



Arbre::Arbre (void)
// appel du constructeur de Element
	: Element ()
{
	Nb_fils = Max_fils = 0;
	pere = NULL;
	fils = NULL;
}

//----------------------------------------------------------
//	Constructeur de Arbre
//		cree un nouvel arbre a partir de A
//		les fils existant seront aussi recres selon opt
//	Entree: A arbre a dupliquer
//		opt :	opt = 1 => copie des fils et du bitmap
//			opt = 0 => pas de copie de bitmap ni des fils
//			opt = 2 => copie du bitmap sans les fils
//----------------------------------------------------------
Arbre::Arbre (Arbre & A, int opt)
	: Element (A, opt)
{
	fils = NULL;
	pere = NULL;
	Nb_fils = Max_fils = 0;

	if (opt == 1)
		// construction des fils selon les fils de A
		for (int i=0 ; i<A.Nb_fils ; i++)
			// creation des fils
			add_fils (new Arbre (*A.fils[i], opt));	
}

//----------------------------------------------------------
//	Destructeur de la classe Arbre
//	l'arbre est detruit ainsi que ces fils
//	( le destructeur de Element est appele automatiquement )
//----------------------------------------------------------
Arbre::~Arbre (void)
{
   // liberation de la memoire : pointeur fils
   free ();
}

//----------------------------------------------------------
//	Arbre :: free
//		efface toute les donnees de 'arbre
//		les fils sontdetruit;
//----------------------------------------------------------
void	Arbre::free (void)
{
	// si l'arbre a un pere, alors il faut le retirer de la liste
	if (pere)
		pere->retire_fils (this);
	// detruit tous les fils existant
	int i;
	for (i=0 ; i<Nb_fils ; i++)
		if (fils[i])
			delete (fils[i]);
	// desalocation du pointeur de fils
	// ::free indique la fonction free de malloc.h
	if (fils)
		::free (fils);
	Max_fils = 0;
	// liberation memoire de Element
	Element::free ();
	fils = NULL;
	Max = Nb_fils = 0;
}

//----------------------------------------------------------
//	Arbre :: add_fils
//		ajoute un fils a l'arbre
//	entree: f pointeur sur  l'arbre fils a ajoute
//----------------------------------------------------------
void	Arbre::add_fils (Arbre *f)
{
	// ajoute un fils, le tableu de pointeurs: fils est realoue
	// si celui-ci est plein.
	// Max_fils = Taille maximun du tableau fils
	// Nb_fils = Nombre de fils effectifs dans le tableau fils
	if (Nb_fils >= Max_fils)
	{
		Max_fils ++;
                if (fils == NULL) fils = (Arbre **)malloc (Max_fils * sizeof (Arbre *));
		else fils = (Arbre **)realloc (fils, Max_fils * sizeof (Arbre *));
		assert (fils);
	}
	fils[Nb_fils] = f;

	// mise a jour du pere de f
	fils[Nb_fils]->pere = this;
	Nb_fils ++;
}

//----------------------------------------------------------
//	Arbre :: retire_fils
//		retire le fils f si celui-ci est connu
//----------------------------------------------------------
void	Arbre::retire_fils (Arbre *f)
{
	int	i;
	// cherche le fils a retire
	for (i=0 ; i<Nb_fils ; i++)
		if (fils[i] == f)
			break;

	// le fils f n'existe pas
	if (i == Nb_fils)
		return;

	// decale le tableau de fils vers le bas
	fils[i] = NULL;
	for (i++ ; i<Nb_fils ; i++)
		fils[i-1] = fils[i];
	fils[Nb_fils-1] = NULL;
	Nb_fils --;
	f->pere = NULL;
}

//----------------------------------------------------------
//	Arbre :: ancetre
//		retourne le plus haut des peres
//----------------------------------------------------------
Arbre *		Arbre::ancetre (void)
{
	// si le pere existe alors appel de pere->ancetre
	if (pere)
		return pere->ancetre ();
	// si pere = NULL; c'est l'ancetre
	else
		return this;
}

//----------------------------------------------------------
//	init_arbre_plan
//		construction des arbres
//		rep est un pointeur de Ws.Nb_plan de 
//		Ws.plan->Num_obj Arbre
//	sauvegarde de la segmentation Ws avec les valeurs
//	des coeffs d'ondelettes associes au plan i
//----------------------------------------------------------
static int init_arbre_plan (MultiResol & Wi, Arbre *** rep, Ifloat & Ima, int i)
{
	int		x, y, j, index, Nb_elem_dep;
	// Wi et Ws doivent ettre de meme type

	//deswap le plan i de Wi et Ws
        im_segment(Wi.band(i), Ima, Nb_elem_dep, (float) 0., VisionCleanBord);
        
#if WRITE_SEGMENTE
        char NameIma[80];
cout << "Segmentation: " << "scale " << i+1 << ", Max = " << Nb_elem_dep << endl;
        sprintf(NameIma, "scale_%d.fits", i+1);
        io_write_ima_float(NameIma, Ima);
#endif
 
	if (Nb_elem_dep)
	{
		// allocation des pointeurs sur les arbres au niveau i
		rep[i] = new Arbre * [Nb_elem_dep];
		// allocation individuelle des Arbres
		for (j=0 ; j<Nb_elem_dep ; j++)
			rep[i][j] = new Arbre ();
	}
	// cas ou il n'y a pas d'elements
	else
		rep[i] = NULL;

	for (y=0; y<Wi.size_band_nl(i) ; y++)
		for (x=0 ; x<Wi.size_band_nc(i) ; x++)
		{
			index = (int) Ima(y,x);
			// si il y a un element segmente
			if (index && index <= Nb_elem_dep)
			{
 				// index = numero de l'element a partir de 1
				index--;
                                if (rep[i][index]->Num_ind == -1)
				{
					// positionne l'element
					rep[i][index]->set_pos (x, y);
					rep[i][index]->Num_ech = i;
					rep[i][index]->Num_ind = index;
				}				
				// positionne l'element
 				rep[i][index]->add_pixel (x, y, Wi(i,y,x));
 			}
		}
// cout << "Init arbre END" << endl;

  return Nb_elem_dep;

}

//----------------------------------------------------------
//	Construit_arbre
//		construction des arbres a partir d'une ondelette W
//		celle ci sera seuille selon W.plan[x]->Seuil
//	retour
//		pointeur sur un tableau 2D d'Arbre
//		par exemple rep[2][5] donne l'adresse du
//		5 iem element du niveau 2 (en partant de 0)
//----------------------------------------------------------
Arbre***	construit_arbre (MultiResol & W, Foret & F)
{
	int    i, j, index, x, y;
	Arbre	***rep;
	rep = new Arbre ** [W.nbr_scale()-1];
        int Nl = W.size_ima_nl();
	int Nc = W.size_ima_nc();
        Ifloat Ima(Nl, Nc, "construit arbre");

	// initialisdation des Arbre du premier plan
// cout << "construit arbre 1 " << endl;
 	(F.info[0])->Nb_elem_dep =  init_arbre_plan (W, rep, Ima, 0);
// cout << "construit arbre 2 " << endl;
	// creation des liens inter plans
	for (i=1 ; i< W.nbr_scale()-1 ; i++)
	{
           int Nls = W.size_scale_nl(i);
           int Ncs = W.size_scale_nc(i);

		// initialise le plan i;
                Ima.resize(Nls,Ncs);
		F.info[i]->Nb_elem_dep = init_arbre_plan (W, rep, Ima, i);
// cout << "construit arbre 3 : " << i <<  endl;

		// recherche des liens inter plan par
		// la recherche de maxima
		for (j=0 ; j< F.info[i-1]->Nb_elem_dep ; j++)
		{
			rep[i-1][j]->get_pos_max (x, y);
			if ((W.Set_Transform == TRANSF_PYR) ||
			   ((W.Set_Transform == TRANSF_SEMIPYR) && (i > 1)))
			{
				x /= 2;
				y /= 2;
				index = (int)Ima(y,x);
			}
			else
				index = (int)Ima(y,x);
			if (index && index <= F.info[i]->Nb_elem_dep)
			{
				if (rep[i][index-1] == NULL)
					assert (0);
				rep[i][index-1]->add_fils (rep[i-1][j]);
			}
		}
	}

	// retourne la structure
	return rep;
}

//----------------------------------------------------------
//	Arbre :: write_arbre
//		ecrit les donne d'un arbre dans le flux f
//		ainsi que ses fils
//----------------------------------------------------------
 
void	Arbre::write_arbre (Info & I, int & pos)
{
	float		f;
	int		i;

	assert (I.fptr);
	// ecriture de l'element
	write_element (I, pos);

	// ecriture des fils
	f = (float)Nb_fils;
	fits_write_img (I.fptr, TFLOAT, pos, 1, &f, &I.fptr_st);
	pos ++;
	for (i=0 ; i<Nb_fils ; i++)
		fils[i]->write_arbre (I, pos);
}

//----------------------------------------------------------
//	Arbre :: read_arbre
//		lecture d'un arbre et de ses fils depuis
//		un flux de fichier
//----------------------------------------------------------

void	Arbre::read_arbre (Info & I, int & pos)
{
	int		i, n, isnul;
	float		f, nul;
	Arbre		*a;

	assert (I.fptr);
	free ();
	read_element (I, pos);
	nul = 0.0;
	fits_read_img (I.fptr, TFLOAT, pos, 1, &nul, &f, &isnul ,&I.fptr_st);
	pos ++;
	n = (int)f;
	for (i=0 ; i<n ; i++)
	{
		a = new Arbre ();
		a->read_arbre (I, pos);
		add_fils (a);
	}
} 

//----------------------------------------------------------
//	Arbre :: creat_wavelet
//		creation de la Wavelet minimum autour
//		de l'arbre et ses fils
//	parametre
//		px = translation en x de la wavelet
//		py = translation en y de la wavelet
//		p = numero du plan a ecrire
//----------------------------------------------------------
void	Arbre::creat_wavelet (MultiResol & W, int px, int py, int p)
{
	if (Num_ech > p)
		for (int i=0 ; i<Nb_fils ; i++)
			fils[i]->creat_wavelet (W, px, py, p);
	else if (p == Num_ech)
	{
		if (Num_ech > W.nbr_scale())
			assert (0);
		if (W.Set_Transform == TRANSF_PYR) 
		{
		    px >>= p; 
		    py >>= p;
		}
		else if  ((W.Set_Transform == TRANSF_SEMIPYR) && (p > 1))
		{
		    px >>= p-1; 
		    py >>= p-1;
		}
 		write_image (W.band(Num_ech), px, py);
	}
}

//----------------------------------------------------------
//	Arbre :: write_wavelet_ech
//		ecrit l'element (la valeur de Num_elem) dans la Wavelet W si
//		cette element et au niveau p
//		sinon : appel recurssif
//----------------------------------------------------------

void	Arbre::write_wavelet_ech (MultiResol & W, int p)
{
 	if (Num_ech > W.nbr_scale() - 1)
		assert (0);
	if (Num_ech > p)
		for (int i=0 ; i<Nb_fils ; i++)
			fils[i]->write_wavelet_ech (W, p);
	else if (Num_ech == p)
		write_image_ech (W.band(Num_ech));
}

//----------------------------------------------------------
//	Arbre :: write_wavelet_seg
//		ecrit l'element (la valeur de Num_elem) dans la Wavelet W si
//		cette element et au niveau p
//		sinon : appel recurssif
//----------------------------------------------------------
void	Arbre::write_wavelet_seg (MultiResol & W, int p)
{
	if (Num_ech > W.nbr_scale() - 1)
		assert (0);
	if (Num_ech > p)
		for (int i=0 ; i<Nb_fils ; i++)
			fils[i]->write_wavelet_seg (W, p);
	else if (Num_ech == p)
		write_image_seg (W.band(Num_ech));
}

//----------------------------------------------------------
//	Arbre :: write_wavelet
//		ecrit l'element dans la Wavelet W si
//		cette element et au niveau p
//		sinon : appel recurssif
//----------------------------------------------------------
void	Arbre::write_wavelet (MultiResol & W, int p)
{
	if (Num_ech > W.nbr_scale() - 1)
		assert (0);
	if (Num_ech > p)
		for (int i=0 ; i<Nb_fils ; i++)
			fils[i]->write_wavelet (W, p);
	else if (Num_ech == p)
		write_image (W.scale(Num_ech));
}


//----------------------------------------------------------
//	Arbre :: indentif_obj
//		identification des objets par maximun locales
//		les objets ne sont pas numerotes dans cette fonction
//		les sous objets sont identifies
//----------------------------------------------------------
void		Arbre::identif_obj(int so, int erg, Foret & F, float dmax)  
{
	float	val = 0.0, d, dmin = 1e10,ft;
	int	i, dx, dy, xp, yp, xf, yf, fm = -1;
	Ifloat	*Ip, *I, *temp;

	// BUG: small structures with a father (certainly real objects)
        //      disapear with the followin line. Replaced by the next one.
	// if (Num_ech && Nb_pixel <= (1<<Num_ech)) val = 1e10;

	// Test the father before deleting an object.
	// if a father exists, the object is not eliminated,
	// even if the size of the structure is small.
	if (Num_ech && Nb_pixel <= (1<<Num_ech) && !pere)
	{
	   val = 1e10;
	   
	   // cout << "Xpos = " << Xpos << "Ypos = " << Ypos <<  endl;
	   // cout << "Num_ech = " << Num_ech << endl;
	   // cout << "Nb_pixel = " << Nb_pixel << endl;
	   // cout << "(1<<Num_ech) = " << (1<<Num_ech) << endl;
	}

	// nouvelle regle a implanter:
	// si so != 0 (donc il y a un objet aux superieure)
	// tenir compte de l'alignement du Max en cours avec celui du pere
	// ou ses ancetres
        if (so)
        {
           int	x, y, xp, yp;
	   float d;

           //calcul la distance entre les max
	   get_pos_max (x, y);
	   pere->get_pos_max (xp, yp);
	   if (test_mvm_trans((set_transform) F.F_TO, Num_ech))
 	   {
		xp *= 2;
		yp *= 2;
	   }
	   d = sqrt((float)((x-xp)*(x-xp) + (y-yp)*(y-yp)));

	   // si les max sont alignes alors
	   // pas de sous objet
	   if (d <= dmax)
	   {
		 val = 1e10;
	   }
         }
        
	// recherche par rapport au pere
	if (pere &&  !val)
	{
	    int tx,ty;

		// dans le cas version ou les fils ne sont pas projetes
		// sur le pere
		// il faut faire: if (pere->Val_max > Val_max)
		//  then *this n'est pas une objet
		if (pere->Nb_fils > 1)
		{
 		Ip = pere->creat_image (tx,ty);
		first_pixel (*Ip, dx, dy);
		xp = pere->Xpos - dx;
		yp = pere->Ypos - dy;
 	        if (test_mvm_trans((set_transform) F.F_TO, Num_ech))
  		{
				temp = zoom_2 (*Ip);
				delete Ip;
				Ip = temp;
				xp *= 2;
				yp *= 2;
		}
		I = creat_image (tx,ty);
		first_pixel (*I, dx, dy);
		xf = Xpos - dx;
		yf = Ypos - dy;
		masque (*Ip, *I, xf - xp, yf - yp);
		val = max (*Ip);
		delete Ip;
		delete I;
		}
	}
	// recherche par rapport aux fils
	if (erg)
		ft = Val_max*(1<<(Num_ech+1));
	else
		ft = Val_max;
	
		
	if (val <  ft)
	{
	        // recherche le fils le plus proche spatialement
		for (i=0 ; i<Nb_fils ; i++)
		{
			get_pos_max (xp, yp);
			fils[i]->get_pos_max (xf, yf);
	                if (test_mvm_trans((set_transform) F.F_TO, Num_ech-1))
		    	{
				xp *= 2;
				yp *= 2;
				dx = xp - xf;
				dy = yp - yf;
			}
			else
			{
				dx = xp - xf;
				dy = yp - yf;
			}
			d = (float)dx * dx + (float)dy * dy;
			if (dmin > d)
			{
				dmin = d;
				fm = i;
			}
		}
	}

	// test sur les fils selon les criteres suivants
	// le fils doit etre assez grand si l element en cours
	// est au moins a l echelle 2
	if (fm != -1 && (Num_ech>=2) && fils[fm]->Nb_pixel <= (1<<Num_ech-1))
		fm = -1;

	if (fm != -1)
	{
		val = fils[fm]->Val_max;
		if (erg)
			val *= 1<<Num_ech;
	}
        if (erg)
		ft = Val_max*(1<<(Num_ech+1));
	else
		ft = Val_max;
	
	// nouvelle objet detecte			  
	if (val < ft)
	{
		if (so)
		{
			S_obj = so;
			F.info[Num_ech]->Nb_s_obj++;
		}
 		so ++;
		F.info[Num_ech]->Nb_obj ++;
		// numero temperaire: Foret.Identif_obj refera le travail
		Num_obj = F.info[Num_ech]->Nb_obj;
 	}

	// appel recurssif sur les fils
	for (i=0 ; i<Nb_fils ; i++)
		fils[i]->identif_obj (so, erg, F, dmax);
}

//----------------------------------------------------------
//	Arbre :: numerote_obj
//----------------------------------------------------------
void	Arbre::numerote_obj (int * num)
{
	int		i;
	if (Num_obj != -1)
		Num_obj = num[Num_ech]++;
	for (i=0 ; i<Nb_fils ; i++)
		fils[i]->numerote_obj (num);
}

 
//----------------------------------------------------------
//	Arbre :: numerote_arbre
//----------------------------------------------------------
void	Arbre::numerote_arbre (int * num)
{
	int		i;
	Num_ind = num[Num_ech]++;
	for (i=0 ; i<Nb_fils ; i++)
		fils[i]->numerote_arbre (num);
}


//----------------------------------------------------------
//	Arbre :: creat_obj
//		creation de la Foret des objets a partir
//		des arbres s'il sont des objets
//----------------------------------------------------------
void	Arbre::creat_obj (Foret & F)
{
	Arbre	*A, *P;
	Ifloat	*Ip, *I, *temp;
	int	xp, yp, xf, yf, dx, dy;
	if (Num_obj != -1)
	{
		A = new Arbre (*this, 1);
		if (pere)  
		{
			P = new Arbre (*pere, 2);
			P->Proj = 1;
			if (pere->Nb_fils == 1 && !S_seg)
			{
				P->add_fils (A);
				F.add_arbre (P, Num_ech);
			}
			else
			{
			   int tx,ty;
				Ip = P->creat_image (tx,ty);
				first_pixel (*Ip, dx, dy);
				xp = P->Xpos - dx;
				yp = P->Ypos - dy;
		                if (test_mvm_trans((set_transform) F.F_TO, Num_ech))  
   					{
						temp = zoom_2 (*Ip);
						delete Ip;
						Ip = temp;
						xp *= 2;
						yp *= 2;
					}
				I = creat_image (tx,ty);
				first_pixel (*I, dx, dy);
				xf = Xpos - dx;
				yf = Ypos - dy;
				masque (*Ip, *I, xf - xp, yf - yp);
				delete I;
				if (test_mvm_trans((set_transform) F.F_TO, Num_ech))
 				{
						temp = decime_2 (*Ip);
						delete Ip;
						Ip = temp;
						xp /= 2;
						yp /= 2;
				}
				first_pixel (*Ip, dx, dy);
				P->read_image (*Ip, dx, dy);
				delete Ip;
				P->Xpos = xp + dx;
				P->Ypos = yp + dy;
				// test si le projet sur le pere n'est pas nul
				if (P->Nb_pixel)
				{
					P->add_fils (A);
					F.add_arbre (P, Num_ech);
				}
				else
				{
					delete P;
					F.add_arbre (A, Num_ech);
				}
			}
		}
		else
			F.add_arbre (A, Num_ech);
	}
	for (int i=0 ; i<Nb_fils ; i++)
		fils[i]->creat_obj (F);
}

//----------------------------------------------------------
//	Arbre :: taille_obj
//		retorune la taille d'un objet et de ses fils
//----------------------------------------------------------
void	Arbre::taille (Foret & F, int & xmin, int & ymin, int & xmax, int & ymax)
{
	int	xm, ym, xM, yM;
	taille_x (xm, xM);
	taille_y (ym, yM);
	xm += Xpos;
	xM += Xpos;
	ym += Ypos;
	yM += Ypos;
 	if (F.F_TO == (int) TRANSF_PYR)
 	{
 		xm <<= Num_ech;
		ym <<= Num_ech;
		xM <<= Num_ech;
		yM <<= Num_ech;
	}
	else if ((F.F_TO == (int) TRANSF_SEMIPYR) && (Num_ech>1))
	{
		xm <<= Num_ech-1;
		ym <<= Num_ech-1;
		xM <<= Num_ech-1;
		yM <<= Num_ech-1;
	}
	xmin = MIN (xmin, xm);
	ymin = MIN (ymin, ym);
	xmax = MAX (xmax, xM);
	ymax = MAX (ymax, yM);
	for (int i=0 ; i<Nb_fils ; i++)
		fils[i]->taille (F, xmin, ymin, xmax, ymax);
}

//----------------------------------------------------------
//	Arbre :: get_arbre
//		retourne l'arbre de niveau p et de Num_elem p
//----------------------------------------------------------
Arbre *		Arbre::get_arbre (int p, int n)
{
	Arbre	*rep;
	if (p < Num_ech)
		for (int i=0 ; i<Nb_fils ; i++)
		{
			rep = fils[i]->get_arbre (p, n);
			if (rep)
				return rep;
		}
	else if (p == Num_ech && Num_ind == n)
		return this;
	else
		return NULL;
	return NULL;
}

//----------------------------------------------------------
//	Arbre :: get_obj
//		retourne l'arbre de niveau p et de Num_obj n
//----------------------------------------------------------
Arbre *		Arbre::get_obj (int p, int n)
{
	Arbre	*rep;
	if (p < Num_ech)
		for (int i=0 ; i<Nb_fils ; i++)
		{
			rep = fils[i]->get_obj (p, n);
			if (rep)
				return rep;
		}
	else if (p == Num_ech && Num_obj == n)
		return this;
	else
		return NULL;
	return NULL;
}

//----------------------------------------------------------
//	Arbre :: get_pere_obj
//		retourne le pere de l'objet demande
//----------------------------------------------------------
Arbre *		Arbre::get_pere_obj (void)
{
	if (Num_obj != -1)
		return this;
	else if (pere)
		return pere->get_pere_obj ();
	else
		return NULL;
}

//----------------------------------------------------------
//	Arbre :: change_arbre
//----------------------------------------------------------
int		Arbre::change_arbre (Ifloat & I, Ifloat & Iss, Foret & F)
{
	Arbre	**temp;
	int	n_ss, tx, ty, x, y, i, index=0;
	n_ss = (int) max(Iss);
	if (n_ss == 1)
		return n_ss;

	// calcul du vecteur de translation du bord de l'image
	first_pixel (I, x, y);
	tx = Xpos - x;
	ty = Ypos - y;

	// creation des nouveau fils
	temp = new Arbre *[n_ss];
	for (i=0 ; i<n_ss ; i++)
	{
		temp[i] = new Arbre ();
		temp[i]->Num_ech = Num_ech;
		temp[i]->read_image (I, Iss, i+1);
		temp[i]->Xpos += tx;
		temp[i]->Ypos += ty;
		temp[i]->S_seg = 1;
		if (!temp[i]->Nb_pixel)
		{
			delete temp[i];
			temp[i] = NULL;
			continue;
		}
 		if (i)
		{
			F.info[Num_ech]->Nb_elem ++;
			F.info[Num_ech]->Nb_elem_ss ++;
		}
 		F.add_arbre_fils (temp[i]);
	}

	// cas des fils : creation des nouveaux liens
	for (i=0 ; i<Nb_fils ; i++)
	{
		fils[i]->get_pos_max (x, y);
		if ((F.F_TO == (int) TRANSF_PYR) 
		               || ((F.F_TO == (int) TRANSF_SEMIPYR) && (Num_ech)))
		{
		  x/=2;
		  y/=2;
 		}
		x -= tx;
		y -= ty;
		if (y>=0 && x>=0)
			index = (int)Iss(y,x);
		else
			assert (0);
		index--;
		temp[index]->add_fils (fils[i]);
		fils[i] = NULL;
	}
	Nb_fils = 0;
	if (pere)
		pere->retire_fils (this);
	else
		F.retire_arbre (this);
	delete temp;
	return n_ss;
}

//----------------------------------------------------------
//	Arbre::creat_arbre_aligne (float)
//		creation d'un nouvel arbre ascendant
//		ne contenant que les peres et ascendants
//		alligne entre eux.
//----------------------------------------------------------
Arbre	*	Arbre::creat_arbre_aligne (float dmax, Foret &F)
{
	Arbre		*rep, *temp;
	int		x, y, xp, yp;
	float		d;

	rep = new Arbre (*this, 2);

	// fin de recursivite
	if (!pere)
		return rep;

	//calcul la distance entre les max
	get_pos_max (x, y);
	pere->get_pos_max (xp, yp);
  	if ((F.F_TO == (int) TRANSF_PYR) 
		               || ((F.F_TO == (int) TRANSF_SEMIPYR) && (Num_ech)))
	{
		xp *= 2;
		yp *= 2;
	}
	d = (x-xp)*(x-xp) + (y-yp)*(y-yp);

	// si les max sont alignes alors
	// appel recursif
	if (d <= dmax)
	{
		temp = pere->creat_arbre_aligne (dmax, F);
		temp->add_fils (rep);
	}

	return rep;
}

//----------------------------------------------------------
//	Arbre::creat_obj_aligne
//		creation des objets avec alignement
//----------------------------------------------------------
void		Arbre::creat_obj_aligne (Foret & F, float dmax)
{
	Arbre		*temp, *temp2, *a;
	int		i;

	if (Num_obj != -1)
	{
		if (pere)
		{
			temp = pere->creat_arbre_aligne (dmax, F);
			for (a=temp->ancetre () ; a->Nb_fils==1 ; a=a->fils[0])
				a->Proj = 1;
			a->Proj = 1;
			temp2 = new Arbre (*this, 1);
			temp->add_fils (temp2);
			temp = temp->ancetre ();
		}
		else
			temp = new Arbre (*this, 1);
		F.add_arbre (temp, Num_ech);
	}

	for (i=0 ; i<Nb_fils ; i++)
		fils[i]->creat_obj_aligne (F, dmax);
}

