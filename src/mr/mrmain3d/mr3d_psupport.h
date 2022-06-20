#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"

#include "Atrou3DFil.h"
#include "Mr3d_FewEvent.h"
#include "Atrou3D.h"
#include "DefPoint.h"
#include "CPca.h"


#define  CUBE_FABS(x)    (fabs(x) * fabs(x) * fabs(x))
#define  MIN3(x,y,z)	(MIN (x , MIN (y , z)))




extern int  GetOpt(int argc, char *const*argv, char *opts);
void Projection_V0(fltarray & Event_of_V0, intarray & Catalogue_Event, 
		FewEventPoisson & FEP);
void Projection_V0(ArrayPoint & Catalogue, fltarray & Event_of_V0, 
		intarray & Catalogue_Event);
void read_data(fltarray & Event_of_V0, intarray & Catalogue_Event ,
		fltarray & Filtred_Input, Bool Filter_Def,
		FewEventPoisson & FEP,int & Nx,int & Ny,int & Nz);
void init_Nb_Event(intarray & Catalogue_Event,intarray *& Nb_Event,
		FewEventPoisson & FEP,ATROUS_3D_FIL & WT_tools);
void init_support(intarray *& Nb_Event,intarray *& Support, 
		ATROUS_3D_FIL & WT_tools, FewEventPoisson & FEP);

extern int  OptInd;
extern char *OptArg;
Bool WriteAll = False;
Bool KeepPositivSup = False;
char Name_Cube_In[80]; // input file image 
char Name_Support[80];
char AnaFileName[80]="AnaFile";
char TabAnaFileName[80]="TabAnaFile";
char Name_Abaque_Default[80]="Abaque.fits";
char Name_Abaque[80];
char Filtred_FileName[80];
int Abaque_N_Scale = 0;
int Nbr_Plan = 10;
int MinEvent = 5;
int FirstScale = 0;
Bool ObjAna = False;
Bool TabAsciiRes = False;
Bool AsciiRes = False;
Bool Verbose = False;
Bool Approx= True;
Bool Filter_Def=False;
type_border Border=I_ZERO;
Bool WriteSupport=False;
Bool ReadSupport=False;
float epsilon=1e-3;
Bool FitsMod;
float bin_cat=1;
Bool SinglePoint=True;
Bool RemoveBorder=False;
