

#ifndef		_VISTOOL_H_
#define	_VISTOOL_H_

#define CHAINE 80
const float RAD_DEG=180/PI;

/* size and offset in millimiter of the poscript
   image created by the routine mr_write_obj_ps
*/
#define SIZE_IMA_PS_X 160
#define SIZE_IMA_PS_Y 160
#define OFFSET_IMA_PS_X 25
#define OFFSET_IMA_PS_Y 70

struct t_Gauss_2D
	{
	double mx;
	double my;
	double sX;
	double sY;
	double teta;
	double amp;
	double offset;
	};

struct t_Moment_image
	{
	double mx;
	double my;
	double vx;
	double vy;
	double vxy;
	};

Ifloat *		decime_2 (Ifloat & Im);
Ifloat *		zoom_2 (Ifloat & Im, int dx=0, int dy=0);
void complete_bord (Ifloat & Ima, int dx=0, int dy=0);
void masque (Ifloat & Ima, Ifloat& Im);
void masque (Ifloat &Ima, Ifloat& Im, int dx, int dy);
int first_pixel (Ifloat &Ima, int & x, int & y);
Ifloat *	sous_segmente_iter (Ifloat & Im, int Nb_iter);
void add_image (Ifloat & Ima1, Ifloat& Im, int dx, int dy);
t_Moment_image moments(Ifloat &Ima);
float Gauss_2D(float x,float y,t_Gauss_2D Gauss);
void Calcul_param_gauss_2D(t_Moment_image Moment,float Max,t_Gauss_2D& Gauss);
void Estime_param_gauss_2D(t_Gauss_2D& Gauss_est, Ifloat & Im);
void Affiche_param_gauss_2D(const t_Gauss_2D& Gauss);
void convol_5_trou (Ifloat & Im, Ifloat & Im_out, int plan);
void segmente (Ifloat & Im, Ifloat &rep, int & Nb_elem_dep, float S);


#endif
