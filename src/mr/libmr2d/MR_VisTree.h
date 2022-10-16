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
//    Date: 	12/03/96
//    
//    File:  	Arbre.h
//
//----------------------------------------------------------

#ifndef		_ARBRE_H_
	#define	_ARBRE_H_

#include	"MR_VisElem.h"
#include    "IM_Graphics.h"

#define NBR_RECOBJ_METHOD 4
#define XMM_PSF_NL 64
#define XMM_PSF_NC 64

enum type_objrec_method {GRAD_FIX_STEP, GRAD_OPTI_STEP, GRAD_CONJUG, GRAD_PSF, GRAD_PSF_XMM, RECM_UNDEFINED=-1};

#define DEF_RECOBJ_METHOD GRAD_CONJUG
#define DEF_MAX_ITER_RECOBJ 10
#define DEF_ERROR_RECOBJ 1e-5

inline const char * StringRecObj (type_objrec_method type)
{
    switch (type)
    {
        case GRAD_FIX_STEP: 
           return ("reconstruction from the fixed step gradient method");break;
        case GRAD_OPTI_STEP: 
           return ("reconstruction from the optimum step gradient method");break;
        case GRAD_CONJUG: 
           return ("reconstruction from the conjugate gradient method");break;
        case GRAD_PSF: 
           return ("reconstruction using the PSF");break;
	case GRAD_PSF_XMM: 
           return ("reconstruction using the space-variant XMM PSF");break;
        case RECM_UNDEFINED: 
              return ("Unknown reconstruction method ");break;
    }
    return ("Unknown reconstruction method ");
}

class Foret;


inline Bool test_mvm_trans(set_transform Set_Transform, int p)
{
   Bool ValRet = False;
   if ((Set_Transform == TRANSF_PYR) ||
       ((Set_Transform == TRANSF_SEMIPYR) && (p > 0))) ValRet = True;
   return ValRet;
}


//----------------------------------------------------------
//	Classe Arbre
//----------------------------------------------------------
class	Arbre : public Element {
	public:

		int	Nb_fils;
		int	Max_fils;
		Arbre	*pere;
		Arbre	**fils;
	public:
		Arbre (void);
		Arbre (Arbre &, int);
		Arbre (Element &);
		~Arbre (void);
		// liberation memoire des donnes et des fils
		void		free (void);

		// recherche le plus haut des peres

		Arbre *		ancetre (void);

		// ajoute un fils a un arbre
		void		add_fils (Arbre *);

		// retire un fils a un arbre
		void		retire_fils (Arbre *);

		// ecrit un arbre dans un flux
 		void		write_arbre (Info &, int &);

		// lit un arbre a partir d'un flux
 		void		read_arbre (Info &, int &);

 
		// ecrit le bitmap dans un Wavelet au niveau demande
		void		write_wavelet (MultiResol &, int);

		// ecrit la wavelet segmantee au niveau demande
		void		write_wavelet_seg (MultiResol &, int);

		// ecrit la wavelet segmantee au niveau demande
		void		write_wavelet_ech (MultiResol  &, int);

		// identification des objets
 		void		identif_obj (int, int, Foret &, float dmax=0);

		// identification des sous objets
		void		identif_s_obj (int);
 
		// est ce un sous objet
		int		is_s_obj (void);

		// nemerotaion des objets
		void		numerote_obj (int *);

		// numerotaion des arbres en fonction de Num_elem
		void		numerote_arbre (int *);

 	
		// creation d'une foret contenant que les objets
		void		creat_obj (Foret &);

		// retourne la taille minimun autour d'un arbre et ses fils
		void		taille (Foret & F, int &, int &, int &, int &);

		// creation d'une wavelet minimun au niveau demande
		void		creat_wavelet (MultiResol &, int, int, int);

		// retourne le graphe des objets, compatibilite avec Fred
		void		get_arbre_obj (Foret &, Arbre *);

		// retourne l'arbre au niveau et du numero demande
		Arbre *		get_arbre (int, int);

		// retourne l'objet au niveau et du numero demande
		Arbre *		get_obj (int, int);

		// retorune le pere objet
		Arbre *		get_pere_obj (void);

		// change l'arbre en relisant les images TO et seg
		int		change_arbre (Ifloat &, Ifloat &, Foret &);
		
                Arbre* creat_arbre_aligne (float dmax, Foret &F);
                void creat_obj_aligne (Foret & F, float dmax);
	};

#define		FORET_ARBRE	1
#define		FORET_OBJ	2

//----------------------------------------------------------
//	Classe Foret
//----------------------------------------------------------
class	Foret: public Info{
	public:
  
		// information sur le suillage
		Info	*info[MAX_SCALE];

		// information sur le nombre d'arbres
		int	Max_arbre[MAX_SCALE];
		int	Nb_arbre[MAX_SCALE];
		Arbre	**arbre[MAX_SCALE];
	public:
		Foret (void);
		Foret (Foret &);
		~Foret (void);
		
		void		operator = (Foret *);
	
	        // liberation memoire des arbres
		void		free (void);

		// ajoute un arbre a la foret
		void		add_arbre (Arbre *, int);

		// ajoute un arbre a la foret en recherchant son pere
		void		add_arbre_fils (Arbre *);

		// retire un arbre a la foret
		void		retire_arbre (int, int);
		
		// retire un arbre a la foret
		void		retire_arbre (Arbre *);

		// ecrit la foret dans un fichier
		void		write_foret (char *);

		// lit une foret depuis un fichier
		void		read_foret (char *, int = FORET_OBJ);

		// creation de la foret a partir des plans segmentes et initiaux
		void		creation_foret (MultiResol &);

		// ecrit les arbres dans une Wavelet deja allouee
		void		write_wavelet (MultiResol &);

		// ecrit une Wavelet segmentee deja allouee
		void		write_wavelet_seg (MultiResol &);

		// ecrit la WAvelet au niveau demande
		void		write_wavelet (MultiResol &, int, int);

		// creation d'une Wavelet minimum autour d'un arbre demande
		MultiResol *	creat_wavelet (int , int, 
		                               int &dep_x, int &dep_y,
					       int AddBorderX=0,
					       int AddBorderY=0);
		void creat_wavelet (MultiResol & W, int , int, 
		                               int &dep_x, int &dep_y,
					       int AddBorderX=0,
					       int AddBorderY=0);

		// retire les elements isole dans la foret
		void		retire_isole (int);

	        // retire les elements petit dans la foret
		void		retire_petit (int);

	        // identification des objets et numeroation
		void		identif_obj (float, float dmax=0.);
 
		// retourne si l'objet et un sous objet
		int		is_s_obj (int, int);

	        // retourne si l'objet et un sous objet
		int		is_s_obj (int, int, int);

	        // numerote les objets
		void		numerote_obj (void);

		// numerote les arbres
		void		numerote_arbre (void);

		// retourne le nombre d'objets au niveau demande
		int		get_nb_obj (int);

		// retourne le nombre d'objets a tous les niveaux
		int		get_nb_obj (void);

		// retourne le nombre de sous objets au niveau demande
		int		get_nb_s_obj (int);

		// retourne le nombre de sous objets a tous les niveaux
		int		get_nb_s_obj (void);

 
		// creation d'une foret d'objets
		Foret		*creat_obj (void);

		// retourne l'arbre au niveau et du numero demande
		Arbre *		get_arbre (int, int, Arbre * = NULL);

		// retourne l'objet au niveau et du numero demande
		Arbre *		get_obj (int, int, Arbre * = NULL);

 
		// retourne le numero de l'objet dans le foret
		int		get_nf_obj (int, int);

		// retorune l'objet pere du sous objet demande
		Arbre *		get_pere_obj (Arbre *);

		// sous segmentation de la foret
		void		sous_segmente (int, int);

 		// lecture des informations sur fichier fits
		void		read_foret_info (void);

		// lecture des informations sur fichier fits
		void		write_foret_info (void);
		
                Foret *	creat_obj_aligne (float dmax);

 	};

//----------------------------------------------------------
//	autres fonctions definies dans Foret.C et Arbre.C
//----------------------------------------------------------

Arbre***	construit_arbre ( MultiResol &, Foret & F);
void		write_obj (Foret &, Foret &, char *);
Foret *		init_arbre_image (Ifloat &, Ifloat &);


void mr_make_graph(MultiResol & W, Foret & F, Bool SousSegment, Bool KillIsolObj=True);
void mr_make_obj_tree(Foret & F, Foret * & F_obj, Bool UseEnerg=False,
                     Bool Align=False, float Dmax=1);
void mr_make_obj_tree(MultiResol & W, Foret * & F_obj, Bool SousSegment, 
                      Bool KillIsolObj=True, Bool UseEnerg=False,
                      Bool Align=False, float Dmax=1);
void mr_sm_recons_obj (MultiResol& MR_Data, Ifloat &Imag, int Nb_iter,
		       double& Error, float eps);

void mr_psf_recons_obj (MultiResol& MR_Data, Ifloat &Imag, Ifloat &Psf,
                        int Nb_iter, double& Error, float eps, Bool Deconv,
			float Alpha, float Zoom);

void mr_recons_obj (MultiResol& W, Ifloat &Im_rec, type_objrec_method Meth, 
                     int Nb_iter, double& Erreur, float eps, 
		     Ifloat & Psf, float Fwhm, Ifloat & IMGauss, 
		     Bool Deconv=False, float Zoom=1.);

#endif
