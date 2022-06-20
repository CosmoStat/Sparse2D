#ifndef C_DENOISE_IMAGE_H
#define C_DENOISE_IMAGE_H

#include <TempArray.h>

class C_DENOISE_IMAGE
{
public:
	C_DENOISE_IMAGE(dblarray &image,dblarray &Dico);
	~C_DENOISE_IMAGE();
	dblarray denoise_image(dblarray &NoisyImage, dblarray &Dico, int OverlapNumber,int SparsityTarget,double ErrorTarget,intarray full_image_dim,bool Verb);
	dblarray denoise_image(dblarray &NoisyImage, dblarray &Dico, int OverlapNumber,int SparsityTarget,double ErrorTarget,intarray full_image_dim){return denoise_image(NoisyImage,Dico,OverlapNumber,SparsityTarget,ErrorTarget,full_image_dim,False);}
	dblarray denoise_image(dblarray &NoisyImage, dblarray &Dico, int OverlapNumber,double ErrorTarget,intarray full_image_dim){return denoise_image(NoisyImage,Dico,OverlapNumber,Dico.ny(),ErrorTarget,full_image_dim,False);}
	dblarray denoise_image(dblarray &NoisyImage, dblarray &Dico, int OverlapNumber,int SparsityTarget,intarray full_image_dim){return denoise_image(NoisyImage,Dico,OverlapNumber,SparsityTarget,-1.,full_image_dim,False);}
	dblarray extract_patches(dblarray &NoisyImage,int &OverlapNumber);
	intarray compute_full_image_size(dblarray &Image, int &OverlapNumber);
protected:
	int Npix;
private:
};

#endif // C_DENOISE_IMAGE_H

