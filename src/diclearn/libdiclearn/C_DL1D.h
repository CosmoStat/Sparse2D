#ifndef C_DL1D_H
#define C_DL1D_H

#include <TempArray.h>

class C_DL1D
{
public:
	C_DL1D(dblarray &training_set);
	~C_DL1D();
	dblarray dl1d(dblarray &training_set, dblarray &initD, int IterationNumber,int SparsityTarget,double ErrorTarget,bool Verb);
	dblarray dl1d(dblarray &training_set, dblarray &initD, int IterationNumber,int SparsityTarget,double ErrorTarget){return dl1d(training_set,initD,IterationNumber,SparsityTarget,ErrorTarget,False);}
	dblarray dl1d(dblarray &training_set, dblarray &initD, int IterationNumber,double ErrorTarget){return dl1d(training_set,initD,IterationNumber,initD.ny(),ErrorTarget,False);}
	dblarray dl1d(dblarray &training_set, dblarray &initD, int IterationNumber,int SparsityTarget){return dl1d(training_set,initD,IterationNumber,SparsityTarget,-1.,false);}

	/*dblarray dl1d(dblarray &training_set,int Na,);*/

	double eps;
protected:
	int Npix; // number of pixels per atom
	int Na; // number of atoms
private:
};

#endif // C_DL1D_H
