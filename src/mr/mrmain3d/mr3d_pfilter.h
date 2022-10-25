 
#include "IM_Obj.h"
#include "IM3D_IO.h"
#include "Atrou3D.h"
#include "Atrou3DFil.h"
#include "IM_IO.h"			
#include "Mr3d_FewEvent.h"
// #include "MR_Obj.h"
#include "DefPoint.h"

#define  CUBE_FABS(x)    (fabs(x) * fabs(x) * fabs(x))
#define  MIN3(x,y,z)	(MIN (x , MIN (y , z)))


extern int  OptInd;
extern char *OptArg;

		
class User_Param {
private :
	void put_default_val();
public :
	User_Param(int argc, char *argv[]);
	void usage(char *argv[]);
	void WriteOptProc();
	
	
	char Name_Cube_In[256]; /* input file image */
	char Name_Cube_Out[256]; /* output file name */
	int Nbr_Plan;  
	int FirstScaleDetect;
	int NEGFirstScaleDetect;
	int MinEvent;
	int Rythm;
	int NAutoConv;
	int NbIter;
	type_border Border;
	float epsilon;
	float bin_cat;
	float N_Sigma;
	Bool Verbose;
	Bool ProjV0;
	Bool WriteRes;
	Bool WriteDecomp;
	Bool InfoBand; 
	Bool Read_Histo;
	Bool WriteHisto;
	Bool SinglePoint;
	Bool OnlyPosVal;
	Bool PositiveImage;
	Bool WriteOpt;
	Bool AddLastScale;
	Bool ReadSupport;
	Bool WriteSupport;
	Bool Approx;
	Bool Adjoint;
	Bool RemoveBorder;
	Bool FitsMod;
	Bool CFAIter;
};	


	
class Iteration {
public :
	
	Iteration(fltarray & Event_of_V0,fltarray & Output_Cube,
		intarray *& Support,intarray *& Nb_Event_In_Boxes, 
		fltarray & Mean_Nb_Ev, ATROUS_3D_FIL & WT_tools,
		FewEventPoisson & FEP, User_Param & Param);
	void operator++();
	void operator++(int x);
	void init();
	bool test_end(){return(_CurrentIter < _NbIter);}
	~Iteration();
private :
	int _Nx,_Ny,_Nz;
	int _NbIter;
	int _CurrentIter;
	User_Param * _Param;
	fltarray * _Output_Cube;
	intarray * _Support;
	fltarray * _Event_of_V0;
	FewEventPoisson * _FEP;
	ATROUS_3D_FIL * _WT_tools;
	fltarray _Residue;
	fltarray _Result;
	intarray * _Nb_Event_In_Boxes;
	fltarray * _Mean_Nb_Ev;
};

class Iteration_CFA {
public :
	
	Iteration_CFA(fltarray & Event_of_V0,fltarray & Output_Cube,
		intarray *& Support,intarray *& Nb_Event_In_Boxes, 
		fltarray & Mean_Nb_Ev, ATROUS_3D_FIL & WT_tools,FewEventPoisson & FEP, 
		User_Param & Param);
	void operator++();
	void operator++(int x);
	void init();
	bool test_end(){return(_lambda >= 0.);}
	~Iteration_CFA();
private :
	int _Nx,_Ny,_Nz;
	int _CurrentIter;
	User_Param * _Param;
	fltarray * _Output_Cube;
	intarray * _Support;
	fltarray * _Event_of_V0;
	FewEventPoisson * _FEP;
	fltarray * _Mean_Nb_Ev;
	ATROUS_3D_FIL * _WT_tools;
	fltarray *_Wavelet_Coef_Init;
	intarray * _Nb_Event_In_Boxes;
	float _Lmax;
	int _NbIter;
	float _delta;
	float _local_threshold;
	float _lambda;
	fltarray _tab_sigma;
};


extern int  GetOpt(int argc, char *const*argv, char *opts);


void read_data(fltarray & Event_of_V0, intarray & Catalogue_Event ,
		FewEventPoisson & FEP,int & Nx,int & Ny,int & Nz,
		User_Param & Param);
void Projection_V0(fltarray & Event_of_V0, intarray & Catalogue_Event,
		FewEventPoisson & FEP, User_Param & Param);
void Projection_V0(ArrayPoint & Catalogue, fltarray & Event_of_V0, 
		intarray & Catalogue_Event, User_Param & Param);
void init_Nb_Event(intarray *& Catalogue_Event,intarray *& Nb_Event,
		fltarray & Mean_Nb_Ev, FewEventPoisson & FEP, 
		ATROUS_3D_FIL & WT_tools, User_Param & Param);
void init_support(intarray *& Nb_Event,intarray *& Support,fltarray &
		Mean_Nb_Ev, ATROUS_3D_FIL & WT_tools,FewEventPoisson & FEP, 
		int Nx,int Ny, int Nz,
		User_Param & Param);
void write_result(fltarray & Output_Cube,User_Param & Param);
void write_result_inter(fltarray & Output_Cube,int s,User_Param & Param,char text[]);
