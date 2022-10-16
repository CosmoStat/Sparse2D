

#ifndef	_GRAPHICS_H_
#define	_GRAPHICS_H_

#define VISU3D_LEVEL_1 100
#define VISU3D_LEVEL_2 100

#define VISU_BLACK_WHITE 0
#define VISU_MULTILEVEL 1

#define SIZE_TAB_TRIG 37
#define SIZE_TAB_TRIG_2 18

#define	RINT(x)	(int)(floor(x+0.5))

class C_Trigo {
     public:
      C_Trigo()
      {
         for (int i=0; i<SIZE_TAB_TRIG; i++)
         {
            ctg[i] = cos(i*PI/SIZE_TAB_TRIG_2);
            stg[i] = sin(i*PI/SIZE_TAB_TRIG_2);
         }
      };
      double ctg[SIZE_TAB_TRIG];
      double stg[SIZE_TAB_TRIG];
      ~C_Trigo(){};
};

void im_draw_line(Ifloat & Ima, double Begin_X, double Begin_Y,
                  double End_X, double End_Y, float Val);
void im_draw_circle(Ifloat &Ima, double x, double y, double r, float Val);
void im_draw_ellips(Ifloat &Ima, double x, double y, double a,
		    double b, double theta, float Val, Bool dotflag=False);

void im_isophot (Ifloat &Pict, Ifloat &Tab_Iso, float seuil_min, 
            float Seuil_Max, float pas, int Mod_Visu=VISU_BLACK_WHITE);
void im_quick_isophot (Ifloat &Pict, Ifloat &Tab_Iso, float seuil_min, float Seuil_Max, float pas, int Mod_Visu=VISU_BLACK_WHITE);

void im_3dview_1 (Ifloat &Imag, Ifloat &Tab_Caval, float angle, int Step);
void im_3dview_2 (Ifloat &Tab_Img, Ifloat &Tab_Caval, 
                 int Inc=1, int Mod_Visu=VISU_BLACK_WHITE);
int im_create_ps_file (Ifloat &Pict, char *File_Name);


#ifdef VISUX11
void  im_x11_disp_ima (Ifloat & Pict, int H, int  L, char *Mes1, char *Mes2);
#endif


#endif
