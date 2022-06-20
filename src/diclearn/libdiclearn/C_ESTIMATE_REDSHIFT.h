#ifndef C_ESTIMATE_REDSHIFT_H
#define C_ESTIMATE_REDSHIFT_H

#include <TempArray.h>

class C_ESTIMATE_REDSHIFT
{
public:
	C_ESTIMATE_REDSHIFT(dblarray &Dictionary);
	~C_ESTIMATE_REDSHIFT();
	dblarray estimate_redshift(dblarray &Data,intarray &Shift, int SparsityTarget, double ErrorTarget, dblarray &rms,bool Verb);
	dblarray estimate_redshift(dblarray &Data,intarray &Shift, int SparsityTarget , dblarray &rms,bool Verb){return estimate_redshift(Data,Shift,SparsityTarget,-1,rms,Verb);}
	dblarray estimate_redshift(dblarray &Data,intarray &Shift, double ErrorTarget , dblarray &rms,bool Verb){return estimate_redshift(Data,Shift,Npix,ErrorTarget,rms,Verb);}
	int Npix,Na;
	//	dblarray denoise_image(dblarray &NoisyImage, dblarray &Dico, int OverlapNumber,int SparsityTarget,intarray full_image_dim){return denoise_image(NoisyImage,Dico,OverlapNumber,SparsityTarget,-1.,full_image_dim,False);}

protected:
	dblarray Dico;
private:
};

#endif // C_ESTIMATE_REDSHIFT.h

