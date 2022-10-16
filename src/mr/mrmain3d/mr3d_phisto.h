 
//#include "IM_Obj.h"
//#include "IM3D_IO.h"
//#include "Atrou3D.h"
//#include "IM_IO.h"			
#include "Mr3d_FewEvent.h"
//#include "MR_Obj.h"
//#include "DefPoint.h"
//#include "3D_tools.h"


//#define  CUBE_FABS(x)    (fabs(x) * fabs(x) * fabs(x))
//#define  MIN3(x,y,z)	(MIN (x , MIN (y , z)))


extern int  OptInd;
extern char *OptArg;

		
class User_Param {
private :
	void put_default_val();
public :
	User_Param(int argc, char *argv[]);
	void usage(char *argv[]);
	void WriteOptProc();
	
	int Rythm;
	int NAutoConv;
	Bool Verbose;
	Bool WriteOpt;
	Bool WriteRes;
	Bool DefaultOpt;
};	



extern int  GetOpt(int argc, char *const*argv, char *opts);
