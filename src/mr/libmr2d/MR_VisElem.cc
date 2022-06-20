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
//    File:  	Elemnt.C
//
//----------------------------------------------------------

#include        <unistd.h>
#include        "IM_Obj.h"
#include	"MR_Obj.h"
#include        "IM_VisTool.h"
#include	"MR_VisElem.h"
#include	"MR_VisTree.h"

#define PRINT 0
#define ELEMENT_STOP -1.e20

//----------------------------------------------------------
//	Variable globale
//----------------------------------------------------------
float		Element::Stop = ELEMENT_STOP;


//----------------------------------------------------------
//	Constructeur de Element
//----------------------------------------------------------

Element::Element (void)
{
	image = NULL;
	Num_ech = Num_ind = Num_obj = -1;
	S_obj = Proj = 0;
	Taille = Nb_pixel = Max = S_seg = 0;
	Xset = Yset = Xpos = Ypos = Xmax = Ymax = 0;
	Val_max = 0.0;
}
//----------------------------------------------------------
//	Constructeur de copie de Element
//----------------------------------------------------------

Element::Element (Element & E, int opt)
{
	if (opt)
	{
		image = (float *)malloc (E.Taille * sizeof (float));
		Max = Taille = E.Taille;
		Nb_pixel = E.Nb_pixel;
		for (int i=0 ; i<Taille ; i++)
			image[i] = E.image[i];
	}
	else
	{
		image = NULL;
		Max = Taille = Nb_pixel = 0;
	}
	Num_ech = E.Num_ech;
	Num_ind = E.Num_ind;
	Num_obj = E.Num_obj;
	S_obj = E.S_obj;
	Proj = E.Proj;
	Xset = E.Xset;
	Yset = E.Yset;
	Xpos = E.Xpos;
	Ypos = E.Ypos;
	Xmax = E.Xmax;
	Ymax = E.Ymax;
	Val_max = E.Val_max;
	S_seg = E.S_seg;
}
 

//----------------------------------------------------------
//	Destructeur de Element
//----------------------------------------------------------
Element::~Element (void)
{
	free ();
}

//----------------------------------------------------------
//	Element :: free
//----------------------------------------------------------
 
void	Element::free (void)
{
	if (image)
		::free (image);
	else
		return;
	image = NULL;
//	Num_ech = Num_ind = Num_obj = S_obj = 0;
	Taille = Nb_pixel = Max = 0;
	Xset = Yset = Xpos = Ypos = Xmax = Ymax = 0;
	Val_max = 0.0;
}

//----------------------------------------------------------
//	Element :: add_borne
//----------------------------------------------------------

void	Element::add_borne (float x, float y)
{
	if (Taille + 3 > Max)
	{
		Max += ELEMENT_SIZE_BLOC;
		if (image == NULL) image=(float*) malloc (Max * sizeof (float));
		else image = (float*) realloc (image, Max * sizeof (float));
		assert (image);
	}
	image[Taille++] = Stop;
	image[Taille++] = x;
	image[Taille++] = y;
	Xset = (int)x;
	Yset = (int)y;
}


//----------------------------------------------------------
//	Element :: add (float f)
//----------------------------------------------------------
void	Element::add (float f)
{
	if (Taille >= Max)
	{
		Max += ELEMENT_SIZE_BLOC;
		if (image == NULL) image=(float*) malloc (Max * sizeof (float));
                else image = (float*) realloc (image, Max * sizeof (float));
 		assert (image);
	}
	image[Taille++] = f;
	if (f > Val_max)
	{
		Xmax = Xset;
		Ymax = Yset;
		Val_max = f;
	}
	Xset++;
	Nb_pixel++;
}


//----------------------------------------------------------
//	Element :: add_pixel
//----------------------------------------------------------
 
void	Element::add_pixel (int x, int y, float v)
{
	if (!v)
		return;
	x -= Xpos;
	y -= Ypos;
	if (y == Yset)
	{
		if (x == Xset)
			add (v);
		else if (x > Xset)
		{
			add_borne (x, y);
			add (v);
		}
	}
	else
	{
		add_borne (x, y);
		add (v);
	}
}

//----------------------------------------------------------
//	Element :: taille_x
//		retourne la taille en X de l'image
//----------------------------------------------------------

int	Element::taille_x (int & min, int & max)
{
	int	x, i;
	if (!Nb_pixel)
		return 0;
	min = 0;
	max = 0;
	x = -1;
	for (i=0 ; i<Taille ; i++)
	{
		if (image[i] != Stop)
		{
			x++;
			max = MAX (x, max);
		}
		else
		{
			x = (int)image[++i];
			min = MIN (x, min);
			max = MAX (x, max);
			x--;
			i++;
		}
	}
	return max - min + 1;
} 

//----------------------------------------------------------
//	Element :: taille_y
//		retourne la taille en Y de l'image
//----------------------------------------------------------

int	Element::taille_y (int & min, int & max)
{
	int	y, i;
	if (!Nb_pixel)
		return 0;
	max = y = min = 0;
	for (i=0 ; i<Taille ; i++)
		if (image[i] == Stop)
		{
			i += 2;
			y = (int)image[i];
			min = MIN (y, min);
			max = MAX (y, max);
		}
	return max - min + 1;
} 

//----------------------------------------------------------
//	Element :: creat_image
//----------------------------------------------------------
Ifloat *	Element::creat_image (int & ddx, int & ddy)
{
	Ifloat *Im;
	int	minx, maxx, miny, maxy;
	int	x, y, dx, dy, i;
	dx = taille_x (minx, maxx);
	dy = taille_y (miny, maxy);
	Im = new Ifloat (dy + 2, dx + 2, "creat_image");
	Im->init (0.0);
	x = -minx + 1;
	y = -miny + 1;
 	ddx = Xpos + minx - 1;
	ddy = Ypos + miny - 1;	
	
	for (i=0 ; i<Taille ; i++)
	{
		if (image [i] != Stop)  
			(*Im)(y,x++) = image[i];
		else
		{
			i++;
			x = (int)image[i++] - minx + 1;
			y = (int)image[i] - miny + 1;
		}
	}
	return Im;
}

//----------------------------------------------------------
//	Element :: write_image (Image & Im)
//----------------------------------------------------------
void	Element::write_image (Ifloat & Im, int xp, int yp)
{
	int	minx, maxx, miny, maxy;
	int	x, y, dx, dy, i;
	dx = taille_x (minx, maxx);
	dy = taille_y (miny, maxy);
	xp = Xpos - xp;
	yp = Ypos - yp;
	if (xp + maxx > Im.nc())
		assert (0);
	if (yp + maxy > Im.nl())
		assert (0);
	x = xp;
	y = yp;
	for (i=0 ; i<Taille ; i++)
	{
		if (image[i] != Stop)  
			Im(y,x++) = image[i];
		else
		{
			i++;
			x = xp + (int)image[i++];
			y = yp + (int)image[i];
		}
	}
}

//----------------------------------------------------------
//	Element :: write_image_seg (Image & Im)
//----------------------------------------------------------
void	Element::write_image_seg (Ifloat & Im, int xp, int yp)
{
	int	minx, maxx, miny, maxy;
	int	x, y, dx, dy, i;
	dx = taille_x (minx, maxx);
	dy = taille_y (miny, maxy);
	xp = Xpos - xp;
	yp = Ypos - yp;
	if (xp + maxx >= Im.nc())
		assert (0);
	if (yp + maxy >= Im.nl())
		assert (0);
	x = xp;
	y = yp;
	for (i=0 ; i<Taille ; i++)
	{
		if (image[i] != Stop)  
			Im(y,x++) = Num_ind;
		else
		{
			i++;
			x = xp + (int)image[i++];
			y = yp + (int)image[i];
		}
	}
}

//----------------------------------------------------------
//	Element :: write_image_seg (Image & Im)
//----------------------------------------------------------
void	Element::write_image_ech (Ifloat  & Im, int xp, int yp)
{
 	int	minx, maxx, miny, maxy;
	int	x, y, dx, dy, i;
	dx = taille_x (minx, maxx);
	dy = taille_y (miny, maxy);
	xp = Xpos - xp;
	yp = Ypos - yp;
	if (xp + maxx >= Im.nc())
		assert (0);
	if (yp + maxy >= Im.nl())
		assert (0);
	x = xp;
	y = yp;
	for (i=0 ; i<Taille ; i++)
	{
		if (image[i] != Stop)
		{
			Im(y,x) = Num_ech + 1.0;
			x++;
		}
		else
		{
			i++;
			x = xp + (int)image[i++];
			y = yp + (int)image[i];
		}
	}
}

//----------------------------------------------------------
//	Element :: get_pos_max
//----------------------------------------------------------
void	Element::get_pos_max (int & x, int & y)
{
	x = Xmax + Xpos;
	y = Ymax + Ypos;
}

//----------------------------------------------------------
//	Element :: set_pos
//----------------------------------------------------------
void	Element::set_pos (int x, int y)
{
	if (!Xpos && !Ypos)
	{
		Xpos = x;
		Ypos = y;
	}
}

//----------------------------------------------------------
//	Element :: is_point
//----------------------------------------------------------
int	Element::is_point (int x, int y)
{
	int	xx, yy, i;

	for (xx=yy=i=0 ; i<Taille ; i++)
	{
		if (image[i] == Stop )
		{
			xx = (int)image[++i];
			yy = (int)image[++i];
			continue;
		}
		if (xx+Xpos == x && yy+Ypos == y)
			return 1;
		xx++;
	}
	return 0;
}

//----------------------------------------------------------
//	Element :: read_image
//----------------------------------------------------------
void	Element::read_image (Ifloat & Im, int &ddx, int &ddy)
{
	free ();
	int	x, y;
	for (y=0 ; y<Im.nl() ; y++)
		for (x=0 ; x<Im.nc() ; x++)
			if (Im(y,x) != 0.0)
			{
				set_pos (x, y);
				add_pixel (x, y, Im(y,x));
			}
	// translation de l'image
	Xpos += ddx;
	Ypos += ddy;
}

//----------------------------------------------------------
//	Element :: read_image
//----------------------------------------------------------
void	Element::read_image (Ifloat & I, Ifloat & Iss, int n)
{
	free ();
	int	x, y;
	for (y=0 ; y<I.nl() ; y++)
		for (x=0 ; x<I.nc() ; x++)
			if (Iss(y,x) == (float)n
				&& I(y,x) != 0.0)
			{
				set_pos (x, y);
				add_pixel (x, y, I(y,x));
			}
}

//----------------------------------------------------------
//	Element :: write_element
//----------------------------------------------------------

void	Element::write_element (Info & I, int & pos)
{
	float		temp[13];
  
	temp[0] = (float)Num_ech;
	temp[1] = (float)Num_obj;
	temp[2] = (float)S_obj;
	temp[3] = (float)Proj;
	temp[4] = (float)S_seg;
	temp[5] = (float)Num_ind;
	temp[6] = (float)Nb_pixel;
	temp[7] = (float)Taille;
	temp[8] = (float)Xpos;
	temp[9] = (float)Ypos;
	temp[10] = (float)Xmax;
	temp[11] = (float)Ymax;
	temp[12] = Val_max;

	fits_write_img (I.fptr, TFLOAT, pos, 13, temp, &I.fptr_st);
	pos += 13;
	if (image && Taille)
	{
		fits_write_img (I.fptr, TFLOAT, pos, Taille, image, &I.fptr_st);
		pos += Taille;
	}
	I.erreur_fits();
        
} 

//----------------------------------------------------------
//	Element :: read_element
//----------------------------------------------------------
void	Element::read_element (Info & I, int & pos)
{
	float		temp[13];
	float		nul;
	int		isnul;
 
 	free ();
	nul = 0;
	fits_read_img (I.fptr, TFLOAT, pos, 13, &nul, &temp, &isnul, &I.fptr_st);
	pos += 13;
	 
	Num_ech =	(int)temp[0];
	Num_obj =	(int)temp[1];
	S_obj =		(int)temp[2];
	Proj =		(int)temp[3];
	S_seg =		(int)temp[4];
	Num_ind =	(int)temp[5];
	Nb_pixel =	(int)temp[6];
	Taille =	(int)temp[7];
	Xpos =		(int)temp[8];
	Ypos =		(int)temp[9];
	Xmax =		(int)temp[10];
	Ymax =		(int)temp[11];
	Val_max =	temp[12];
	Max =		Taille;
	if (Taille)
	{
 		image =  (float *) malloc (Taille * sizeof (float));
		assert (image);
		fits_read_img (I.fptr, TFLOAT, pos, Taille, &nul, image,
				&isnul, &I.fptr_st);
		pos +=  Taille;
	}
        I.erreur_fits();
} 

//----------------------------------------------------------
//	Element :: put_element
//----------------------------------------------------------
void	Element::put_element (Element & E)
{
	free ();
	image = E.image;	E.image = NULL;
	Taille = E.Taille;	E.Taille = 0;
	Max = E.Max;		E.Max = 0;
	Nb_pixel = E.Nb_pixel;	E.Nb_pixel = 0;
	Num_ech = E.Num_ech;
	Num_ind = E.Num_ind;
	Num_obj = E.Num_obj;
	S_obj = E.S_obj;
	S_seg = E.S_seg;
	Xset = E.Xset;
	Yset = E.Yset;
	Xpos = E.Xpos;
	Ypos = E.Ypos;
	Xmax = E.Xmax;
	Ymax = E.Ymax;
	Val_max = E.Val_max;
	Proj = E.Proj;
	E.free ();
}

//----------------------------------------------------------
//	constructeur de Info
//----------------------------------------------------------
Info::Info (void)
{
	init ();
}

//----------------------------------------------------------
//	constructeur de info avec un non de fichier
//----------------------------------------------------------
Info::Info (char * name)
{
	init ();
	open_image_read (name);
}

//----------------------------------------------------------
//	Info :: init
//		mise a 0 des champs de Info
//----------------------------------------------------------
void	Info::init (void)
{
	Type_info	= 0;

	// memset (&Wcs, 0, sizeof Wcs);

	Nb_ech		= 0;
	Num_ech		= 0;
	Type_TO		= 0;

	Moyenne		= 0.0;
	Sigma		= 0.0;
	Fond		= 0.0;
	Bruit		= 0.0;
	Seuil		= 0.0;
	Nbpoint		= 0.0;

	Nb_elem		= 0;
	Nb_elem_dep	= 0;
	Nb_elem_net	= 0;
	Nb_elem_ss	= 0;

	Nb_obj		= 0;
	Nb_s_obj	= 0;

	dx		= 0;
	dy		= 0;

	Max		= 0.0;
	Min		= 0.0;
	Nb_classe	= 0;
	Dl_classe	= 0.0;

	Expause		= 0.0;
	ADU		= 0.0;
	RON		= 0.0;

	F_type		= 0;
	F_nl		= 0;
	F_nc		= 0;
	F_wave		= 0;
	F_TO		= 0;
	F_ech		= 0;
	F_nb		= 0;
	F_stop		= 0;

	Obs_lat		= 0.0;
	Obs_long	= 0.0;
	Obs_alt		= 0.0;

	fptr		= NULL;
	fptr_st		= 0;
}

 

//----------------------------------------------------------
//	destructeur de Info
//----------------------------------------------------------
Info::~Info (void)
{
}

//----------------------------------------------------------
//	get_info
//		recopie les infos de Im sur *this
//----------------------------------------------------------
void	Info::get_info (Info & Im)
{
	Type_info	= Im.Type_info;

	Nb_ech		= Im.Nb_ech;
	Num_ech		= Im.Num_ech;
	Type_TO		= Im.Type_TO;

	Moyenne		= Im.Moyenne;
	Sigma		= Im.Sigma;
	Fond		= Im.Fond;
	Bruit		= Im.Bruit;
	Seuil		= Im.Seuil;
	Nbpoint		= Im.Nbpoint;

	Nb_elem		= Im.Nb_elem;
	Nb_elem_dep	= Im.Nb_elem_dep;
	Nb_elem_net	= Im.Nb_elem_net;
	Nb_elem_ss	= Im.Nb_elem_ss;
	Nb_obj		= Im.Nb_obj;
	Nb_s_obj	= Im.Nb_s_obj;

	Nb_obj		= Im.Nb_obj;
	Nb_s_obj	= Im.Nb_s_obj;

	dx		= Im.dx;
	dy		= Im.dy;

	Max		= Im.Max;
	Min		= Im.Min;
	Nb_classe	= Im.Nb_classe;
	Dl_classe	= Im.Dl_classe;

	Expause		= Im.Expause;
	ADU		= Im.ADU;
	RON		= Im.RON;

	F_type		= Im.F_type;
	F_nl		= Im.F_nl;
	F_nc		= Im.F_nc;
	F_wave		= Im.F_wave;
	F_TO		= Im.F_TO;
	F_ech		= Im.F_ech;
	F_nb		= Im.F_nb;
	F_stop		= Im.F_stop;

	Obs_lat		= Im.Obs_lat;
	Obs_long	= Im.Obs_long;
	Obs_alt		= Im.Obs_alt;

	// memcpy (&Wcs, &Im.Wcs, sizeof Wcs);
}

//----------------------------------------------------------
//	get_info
//		recopie les infos de Im sur *this selon le type
//----------------------------------------------------------
void	Info::get_info (Info & Im, int type)
{
	if (type & INFO_TYPE_TO && Im.Type_info & INFO_TYPE_TO)
	{
		Nb_ech		= Im.Nb_ech;
		Num_ech		= Im.Num_ech;
		Type_TO		= Im.Type_TO;
		Type_info	|= INFO_TYPE_TO;
	}

	if (type & INFO_TYPE_IM && Im.Type_info & INFO_TYPE_IM)
	{
		Moyenne		= Im.Moyenne;
		Sigma		= Im.Sigma;
		Fond		= Im.Fond;
		Bruit		= Im.Bruit;
		Seuil		= Im.Seuil;
		Nbpoint		= Im.Nbpoint;
		Type_info	|= INFO_TYPE_IM;
	}

	if (type & INFO_TYPE_ELEM && Im.Type_info & INFO_TYPE_ELEM)
	{
		Nb_elem		= Im.Nb_elem;
		Nb_elem_dep	= Im.Nb_elem_dep;
		Nb_elem_net	= Im.Nb_elem_net;
		Nb_elem_ss	= Im.Nb_elem_ss;
		Nb_obj		= Im.Nb_obj;
		Nb_s_obj	= Im.Nb_s_obj;
		Type_info	|= INFO_TYPE_ELEM;
	}

	if (type & INFO_TYPE_OBJ && Im.Type_info & INFO_TYPE_OBJ)
	{
		Nb_obj		= Im.Nb_obj;
		Nb_s_obj	= Im.Nb_s_obj;
		Type_info	|= INFO_TYPE_OBJ;
	}

	if (type & INFO_TYPE_TRANS && Im.Type_info & INFO_TYPE_TRANS)
	{
		dx		= Im.dx;
		dy		= Im.dy;
		Type_info	|= INFO_TYPE_TRANS;
	}

	if (type & INFO_TYPE_HIST && Im.Type_info & INFO_TYPE_HIST)
	{
		Max		= Im.Max;
		Min		= Im.Min;
		Nb_classe	= Im.Nb_classe;
		Dl_classe	= Im.Dl_classe;
		Type_info	|= INFO_TYPE_HIST;
	}

	if (type & INFO_TYPE_FORET && Im.Type_info & INFO_TYPE_FORET)
	{
		F_type		= Im.F_type;
		F_nl		= Im.F_nl;
		F_nc		= Im.F_nc;
		F_wave		= Im.F_wave;
		F_TO		= Im.F_TO;
		F_ech		= Im.F_ech;
		F_nb		= Im.F_nb;
		F_stop		= Im.F_stop;
		Type_info	|= INFO_TYPE_FORET;
	}

	if (type & INFO_TYPE_EXP && Im.Type_info & INFO_TYPE_EXP)
	{
		Expause		= Im.Expause;
		ADU		= Im.ADU;
		RON		= Im.RON;
		Type_info	|= INFO_TYPE_EXP;
	}

// 	if (type & INFO_TYPE_WCS && Im.Type_info & INFO_TYPE_WCS)
// 	{
// 		memcpy (&Wcs, &Im.Wcs, sizeof Wcs);
// 		Type_info	|= INFO_TYPE_WCS;
// 	}

	if (type & INFO_TYPE_SITEOBS && Im.Type_info & INFO_TYPE_SITEOBS)
	{
		Obs_lat		= Im.Obs_lat;
		Obs_long	= Im.Obs_long;
		Obs_alt		= Im.Obs_alt;
		Type_info	|= INFO_TYPE_SITEOBS;
	}
}


//----------------------------------------------------------
//	Info::open_image_read
//		ouverture d'un fichier fits en lecture
//----------------------------------------------------------
void	Info::open_image_read (char * name)
{
	// test si l'image n'est pas deja en lecture
	//if (fptr && fptr->writemode == READONLY)
	//	return;
	// cas ou l'image est en ecriture
	// if (fptr && fptr->writemode != READONLY)
	//	close_image ();

	// ouverture en lecture seul
	// get_nom (name, file , NOM_FITS);
	fits_open_file (&fptr,  fitsname(name), READONLY, &fptr_st);
	erreur_fits ();
	read_image_info ();
}
//----------------------------------------------------------
//	Info::read_image_info
//		lis tous champs relatifs a la classe, contenus
//		dans le fichier fits
//----------------------------------------------------------
void	Info::read_image_info (void)
{
	// read_image_wcs ();
	read_image_stat ();
}
//----------------------------------------------------------
//	Info::write_image_info
//		ecrit tous les champs valides de la classe
//		dans un fichier fits
//----------------------------------------------------------
void	Info::write_image_info (void)
{
	// write_image_wcs ();
	write_image_stat ();
}

//----------------------------------------------------------
//	Info::open_image_write
//		ouverture d'un fichier fits en lecture
//----------------------------------------------------------
void	Info::open_image_write (char * name)
{
	char		file[CHAINE];
        FILE *F;
        
	// test si l'image n'est pas deja en lecture
	// PB: writemode does not exist anymore in cfitsio
	// if (fptr && fptr->writemode == READWRITE) return;
	// cas ou l'image est en ecriture
	// if (fptr && fptr->writemode != READWRITE) close_image ();

	// ouverture en lecture seul
	// get_nom (name, file , NOM_FITS);
	F = fopen (file,"r");
	if (F)  
	{
	  fclose(F);
 	  unlink (file);
	}
	fits_create_file (&fptr, fitsname(name), &fptr_st);
//	action manuelle
//	write_image_info ();
	erreur_fits ();
}

//----------------------------------------------------------
//	Info::close_image
//		fermeture de l'image
//----------------------------------------------------------
void	Info::close_image (void)
{
	if (!fptr)
		return;
	fits_close_file (fptr, &fptr_st);
	fptr = NULL;
	erreur_fits ();
}

//----------------------------------------------------------
//	Info::erreur_fits
//		gestion des erreurs fits
//----------------------------------------------------------
void	Info::erreur_fits (void)
{
	char		err[CHAINE];

	if (!fptr_st)
		return;
	fits_get_errstatus (fptr_st, err);
	fprintf (stderr, "FITS ERROR : %s.\n", err);
	for (fits_read_errmsg (err) ; err[0] ; fits_read_errmsg (err))
		fprintf (stderr, "           : %s.\n", err);
		
	exit (-1);
}

 

//----------------------------------------------------------
//	Info::read_image_stat
//		lecture des statistiques dans les mots clefs
//		de l'image
//------------------------;----------------------------------
void	Info::read_image_stat (void)
{
	char		com[CHAINE];
	long		lng;

	//----------------------
	fits_read_key_dbl (fptr, (char*)"IM_MOY", &Moyenne, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"IM_SIGMA", &Sigma, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"IM_FOND", &Fond, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"IM_BRUIT", &Bruit, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"IM_SEUIL", &Seuil, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"IM_NBST", &Nbpoint, com, &fptr_st);
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_IM;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_dbl (fptr, (char*)"HIS_MAX", &Max, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"HIS_MIN", &Min, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"HIS_DELT", &Dl_classe, com, &fptr_st);
	fits_read_key_lng (fptr, (char*)"HIS_NBCL", &lng, com, &fptr_st);
	Nb_classe = !fptr_st ? (int)lng : 0;
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_HIST;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_lng (fptr, (char*)"ELE_NB", &lng, com, &fptr_st);
	Nb_elem = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"ELE_DEP", &lng, com, &fptr_st);
	Nb_elem_dep = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"ELE_SS", &lng, com, &fptr_st);
	Nb_elem_ss = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"ELE_NET", &lng, com, &fptr_st);
	Nb_elem_net = !fptr_st ? (int)lng : 0;
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_ELEM;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_lng (fptr, (char*)"OBJ_NB", &lng, com, &fptr_st);
	Nb_obj = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"OBJ_NBSO", &lng, com, &fptr_st);
	Nb_s_obj = !fptr_st ? (int)lng : 0;
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_OBJ;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_lng (fptr, "TRANS_DX", &lng, com, &fptr_st);
	dx = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, "TRANS_DY", &lng, com, &fptr_st);
	dy = !fptr_st ? (int)lng : 0;
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_TRANS;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_lng (fptr, (char*)"TO_ECHM", &lng, com, &fptr_st);
	Nb_ech = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"TO_ECH", &lng, com, &fptr_st);
	Num_ech = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"TO_TYPE", &lng, com, &fptr_st);
	Type_TO = !fptr_st ? (int)lng : 0;
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_TO;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_dbl (fptr, (char*)"EXPTIME", &Expause, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"CCD-GAIN", &ADU, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"CCD-REAN", &RON, com, &fptr_st);
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_EXP;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_lng (fptr, (char*)"F_TYPE", &lng, com, &fptr_st);
	F_type = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"F_NL", &lng, com, &fptr_st);
	F_nl = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr,(char*) "F_NC", &lng, com, &fptr_st);
	F_nc = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"F_TO", &lng, com, &fptr_st);
	F_TO = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"F_wave", &lng, com, &fptr_st);
	F_wave = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"F_ECH", &lng, com, &fptr_st);
	F_ech = !fptr_st ? (int)lng : 0;
	fits_read_key_lng (fptr, (char*)"F_NB", &lng, com, &fptr_st);
	F_nb = !fptr_st ? (int)lng : 0;
	fits_read_key_flt (fptr, (char*)"F_STOP", &F_stop, com, &fptr_st);
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_ELEM;
	else
		erreur_fits ();

	//----------------------
	fits_read_key_dbl (fptr, (char*)"OBS-LAT", &Obs_lat, com, &fptr_st);
	fits_read_key_dbl (fptr, (char*)"OBS-LONG", &Obs_long, com, &fptr_st);
	fits_read_key_dbl (fptr,(char*) "OBS-ALT", &Obs_alt, com, &fptr_st);
	if (fptr_st == KEY_NO_EXIST)
		fptr_st = 0;
	else if (!fptr_st)
		Type_info |= INFO_TYPE_SITEOBS;
	else
		erreur_fits ();
}

//----------------------------------------------------------
//	Info::set_siteobs
//----------------------------------------------------------
void	Info::set_siteobs (double lat, double lng, double alt)
{
	Obs_lat		= lat;
	Obs_long	= lng;
	Obs_alt		= alt;
	Type_info |= INFO_TYPE_SITEOBS;
}

 

//----------------------------------------------------------
//	Info::write_image_stat
//		ecriture des informations statistiques
//		dans un fichier FITS
//----------------------------------------------------------
void	Info::write_image_stat (void)
{
	if (Type_info & INFO_TYPE_IM)
	{
		fits_write_key_fixdbl (fptr, (char*)"IM_MOY", Moyenne, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"IM_SIGMA", Sigma, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"IM_FOND", Fond, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"IM_BRUIT", Bruit, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"IM_SEUIL", Seuil, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"IM_NBST", Nbpoint, 10, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_HIST)
	{
		fits_write_key_fixdbl (fptr, (char*)"HIS_MAX", Max, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"HIS_MIN", Min, 10, (char*)"", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"HIS_DELT", Dl_classe, 10, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"HIS_NBCL", (long)Nb_classe, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_ELEM)
	{
		fits_write_key_lng (fptr, (char*)"ELE_NB", (long)Nb_elem, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"ELE_DEP", (long)Nb_elem_dep, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"ELE_SS", (long)Nb_elem_ss, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"ELE_NET", (long)Nb_elem_net, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_OBJ)
	{
		fits_write_key_lng (fptr, (char*)"OBJ_NB", (long)Nb_obj, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"OBJ_NBSO", (long)Nb_s_obj, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_TRANS)
	{
		fits_write_key_lng (fptr, (char*)"TRANS_DX", (long)dx, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"TRANS_DY", (long)dy, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_TO)
	{
		fits_write_key_lng (fptr, (char*)"TO_ECHM", (long)Nb_ech, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"TO_ECH", (long)Num_ech, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"TO_TYPE", (long)Type_TO, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_EXP)
	{
		fits_write_key_fixdbl (fptr, (char*)"EXPTIME", Expause, 1, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"EXPTIME", "s", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"CCD-GAIN", ADU, 3, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"CCD-GAIN", "ADU/e-", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"CCD-REAN", RON, 3, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"CCD-REAN", "e-", &fptr_st);
	}
	if (Type_info & INFO_TYPE_FORET)
	{
		fits_write_key_lng (fptr, (char*)"F_TYPE", (long)F_type, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_NL", (long)F_nl, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_NC", (long)F_nc, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_TO", (long)F_TO, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_WAVE", (long)F_wave, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_ECH", (long)F_ech, (char*)"", &fptr_st);
		fits_write_key_lng (fptr, (char*)"F_NB", (long)F_nb, (char*)"", &fptr_st);
		fits_write_key_flt (fptr, (char*)"F_STOP", F_stop, 10, (char*)"", &fptr_st);
	}
	if (Type_info & INFO_TYPE_SITEOBS)
	{
		fits_write_key_fixdbl (fptr, (char*)"OBS-LAT", Obs_lat, 10, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"OBS-LAT", "deg", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"OBS-LONG", Obs_long, 10, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"OBS-LONG", "deg", &fptr_st);
		fits_write_key_fixdbl (fptr, (char*)"OBS-ALT", Obs_alt, 10, (char*)"", &fptr_st);
		fits_write_key_unit (fptr, (char*)"OBS-ALT", "m", &fptr_st);
	}
	erreur_fits ();
}
