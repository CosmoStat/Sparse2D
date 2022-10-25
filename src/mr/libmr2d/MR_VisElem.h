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
//    Date: 	28/02/96
//    
//    File:  	Element.h
//
//----------------------------------------------------------

#ifndef		_ELEMENT_H_
#define	_ELEMENT_H_

#define		FIN_LIGNE	0.0
#define		ELEMENT_SIZE_BLOC	20

class Info;
class Arbre;

//----------------------------------------------------------
//	classe Element
//----------------------------------------------------------

class	Element {
	public:
		static float	Stop;	// Valeur du pixel de stop
	public:
		int	Num_ech;	// Numero de l'echelle de l'element
		int	Num_ind;	// Numero normalise de l'element
		int	Num_obj;	// Numero de l'objet
		int	S_obj;		// flag de sous objet
		int	Proj;		// flag de projection sur le pere
		int	S_seg;		// flag de sous segmentation

		float	*image;		// pointeur sur une image en ligne float
		int	Nb_pixel;	// Nombre de pixel
		int	Taille;		// Taille en utile en float de image
		int	Max;		// Taille maximum en float de image

		int	Xpos, Ypos;	// Position du pixel image[0]
		int	Xset, Yset;	// Position courante relative a Pos
		int	Xmax, Ymax;	// Position relative du maximun
		float	Val_max;	// Valeur du maximum en float

	public:
		Element (void);		// Conmstructeur m'est a zero les params
		Element (Element &, int);	// Constructeur de copie
		~Element (void);	// Destructeur

		// libere la place memoire de Elememt
		void		free (void);

		// retourne la taille en x de l'element
		int		taille_x (int &, int &);

		// retoruen la taille en y de l'element
		int		taille_y (int &, int &);

		// ajpoute un pixel a la position courante Xset et Yset
		void		add_pixel (int, int, float);

		// ajoute unbe borne de controle dans image
		void		add_borne (float, float);

		// ajoute un float dans le buffer image
		void		add (float);

		// retorune la valeur du pixel en absolue
		float		get_pixel (int, int);

		// retourne la position en absolue du maximum
		void		get_pos_max (int &, int &);

		// modifie la positon Xpos, et Ypos
		void		set_pos (int, int);

		// ecrit les elements dans une image
		void		write_image (Ifloat &, int = 0, int = 0);

		// ecrit les elements seuilles dans une image
		void		write_image_seg (Ifloat &, int = 0, int = 0);

		// ecrit les elements seuilles en fonction de l'echelle
		void		write_image_ech (Ifloat &, int = 0, int = 0);

		// lit une image pour la stoker dans element
		void		read_image (Ifloat &, int &, int &);

		// lit une image pour la stoker dans element
		void		read_image (Ifloat &, Ifloat &, int);

		// ecrit l'element sur un flux texte - binaire
 		void		write_element (Info &, int &);

		// lit l'element sur un flux texte - binaire
 		void		read_element (Info &, int &);

		// forme une image autour de l'element
		Ifloat *		creat_image (int &, int &);

		// copie l'Element
		void		put_element (Element &);

		// test si x y apartient a l'element
		int		is_point (int, int);
	};


//----------------------------------------------------------
//	Definitions des constantes
//----------------------------------------------------------
#define		INFO_TYPE_IM		1
#define		INFO_TYPE_HIST		2
#define		INFO_TYPE_ELEM		4
#define		INFO_TYPE_OBJ		8
#define		INFO_TYPE_TRANS		16
#define		INFO_TYPE_WCS		32
#define		INFO_TYPE_TO		64
#define		INFO_TYPE_EXP		128
#define		INFO_TYPE_FORET		256
#define		INFO_TYPE_SITEOBS	512
#define		INFO_TYPE_REC		1024


//----------------------------------------------------------
//	Classe Info
//		Information relatives aux :
//			- seuillages
//			- nombres d'elements
//			- nombres d'objet
//		Classe de base de Image
//---------------------------------------------------------- 
class Info {
	public:
	// Type des informations a sauvegarder dans les flux
	int		Type_info;

	// WCS (coordonnees celestes)
	// WorldCoor	Wcs;

	// information sur le type de donnes en ondellettes
	int		Nb_ech;		// nombre d'echelle
	int		Num_ech;	// numero de l'echelle
	int		Type_TO;	// Type de la TO
                                        // nom de la transformee ondelett
                                        // ex TO_PAVE_BSPLINE
	// information sur le seuillage
	double		Moyenne;	// Moyenne de l'image
	double		Sigma;		// eccart type de l'image
	double		Fond;		// Moyenne convergeante
	double		Bruit;		// Ecart type convergeant
	double		Seuil;		// Seuillage de l'image
	double		Nbpoint;	// nombre de poinr pour les stats

	// information pour les histogrammes
	double		Max;
	double		Min;		// Valeur min contenu dans les stats
	double		Dl_classe;	// Largeur d'une classe
	int		Nb_classe;	// nombre de classe

	// nombre d'elements a chaque plan
	int		Nb_elem;	// nombre d'elements au total
	int		Nb_elem_dep;	// Nombre d'elements au depart
	int		Nb_elem_ss;	// nombre d'elements en plus apres sous seg
	int		Nb_elem_net;	// nombre d'elements apres netoyage

	// nombre d'objets
	int		Nb_obj;		// nombre d'objets
	int		Nb_s_obj;	// nombre de sous objets

	//vecteur de translation
	int		dx;		// translation en x
	int		dy;		// translation en y

	// information pour le temps de pose et gain de la camera
	double		Expause;	// temps d'exposition en seconde
	double		ADU;		// gain CCD en e-
	double		RON;		// bruit de lecture CCD

	// information pour la foret
	int		F_type;		// type de la foret
	int		F_nl;		// taille de l'image
	int		F_nc;		// taille de l'image
	int		F_wave;		// nature pave ou non pave 
	int		F_TO;		// type de l'ondelettes
	                                // TRANSF_PYR TRANSF_SEMIPYR 
	                                // ou  TRANSF_PAVE
	int		F_ech;		// nombre de plans
	int		F_nb;		// nombre d'abre
	float		F_stop;		// point d'arret de Element

	// information sur les coordonnee du lieux d'observations
	double		Obs_lat;	// latitude du lieux
	double		Obs_long;	// longitude
	double		Obs_alt;	// altitude

	// gestion des fichier fits
	int		fptr_st;	// flag : erreur FITS
	fitsfile	*fptr;		// pointeur sur le fichier

	public:

		Info (void);
		Info (char *);
		~Info (void);

		// mise a 0 des champs
		void		init (void);

		// I/O des infos dans un flux texte
		void		write_info (FILE *);
		void		read_info (FILE *);
		void		get_info (Info &);
		void		get_info (Info &, int);

		// init
		void		set_siteobs (double, double, double=0.0);

		// test si les variables de histos sont valides
		int		test_def_histo (void);

		// gestion des I/O fits
		// void		read_image_wcs (void);
		void		read_image_stat (void);

		// void		write_image_wcs (void);
		void		write_image_stat (void);

		void		read_image_info (void);
		void		write_image_info (void);

		void		erreur_fits (void);

		void		open_image_read (char *);
		void		open_image_write (char *);
		void		close_image (void);
	};



#endif
