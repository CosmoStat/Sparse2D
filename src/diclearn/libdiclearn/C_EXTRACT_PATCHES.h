#ifndef C_EXTRACT_PATCHES_H
#define C_EXTRACT_PATCHES_H

#include <TempArray.h>

class C_EXTRACT_PATCHES
{
public:
	C_EXTRACT_PATCHES(dblarray &Image);
	~C_EXTRACT_PATCHES();
	dblarray extract(int &W, int &P);
protected:
	dblarray Image;
private:
};

#endif // C_DENOISE_IMAGE_H

