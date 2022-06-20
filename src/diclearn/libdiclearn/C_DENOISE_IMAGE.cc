#include "C_DL1D.h"
#include <vector>
#include <MatrixOper.h>
#include "IM_IO.h"
#include "C_OMP.h"
#include "C_DENOISE_IMAGE.h"
using namespace std;
#include <math.h>


C_DENOISE_IMAGE::C_DENOISE_IMAGE(dblarray &NoisyImage,dblarray &Dico)
{
	Npix = Dico.nx(); // atom size
}
C_DENOISE_IMAGE::~C_DENOISE_IMAGE()
{
}

intarray C_DENOISE_IMAGE::compute_full_image_size(dblarray &Image, int &OverlapNumber)
{
	dblarray full_image; // image extended using symmetric boundary conditions
	double patch_size = sqrt(Npix); // patch size, square root of atom size
	int delta = patch_size - 1; // size of the band to around Image to build full_image
	intarray full_image_dim(2);
	full_image.alloc(Image.nx() + 2*delta, Image.ny() + 2*delta);
	full_image_dim(0) = full_image.nx();
	full_image_dim(1) = full_image.ny();
	return full_image_dim;



}
dblarray C_DENOISE_IMAGE::extract_patches(dblarray &Image,int &OverlapNumber)
{
	dblarray stacked_patches; // patches stacked as columns
	dblarray full_image; // image extended using symmetric boundary conditions
	double patch_size = sqrt(Npix); // patch size, square root of atom size
	int delta = patch_size - 1; // size of the band to around Image to build full_image
	int stackInd,patchInd; // indices used when stacking patches
	int Npatch = 0; // Number of patches extracted
	if (ceilf(patch_size) != patch_size)
		cout << "Error : atom size is not a square" << endl;
	else
	{
		// building full_image using symmetric boundary conditions
		full_image.alloc(Image.nx() + 2*delta, Image.ny() + 2*delta);
		// filling the center of full_image with the pixels of Image
		for (int i=0;i<Image.nx();i++)
			for (int j=0;j<Image.ny();j++)
				full_image(i+delta,j+delta) = Image(i,j);
		// filling top left corner
		for (int i=0;i<delta;i++)
			for (int j=0;j<delta;j++)
				full_image(i,j) = Image(delta-i-1,delta-j-1);
		// filling top right corner. BEWARE, start at the top right pixel of full_image
		for (int i=0;i<delta;i++)
			for (int j=0;j<delta;j++)
				full_image(i,2*delta+Image.ny()-j-1) = Image(delta-i-1,Image.ny()-delta+j);
		// filling bottom left corner. Once again, start at the bottom left pixel
		for (int i=0;i<delta;i++)
			for (int j=0;j<delta;j++)
				full_image(Image.nx()+2*delta-i-1,j) = Image(Image.nx()-delta+i,delta-j-1);
		// filling bottom right corner. Once again, start at the bottom right pixel
		for (int i=0;i<delta;i++)
			for (int j=0;j<delta;j++)
				full_image(Image.nx()+2*delta-i-1,Image.ny()+2*delta-j-1) = Image(Image.nx()-delta+i,Image.ny()-delta+j);

		// filling top side
		for (int i=0;i<delta;i++)
			for (int j=0;j<Image.ny();j++)
				full_image(i,delta+j) = Image(delta -i -1,j);
		// filling bottom side
		for (int i=0;i<delta;i++)
			for (int j=0;j<Image.ny();j++)
				full_image(Image.nx()+2*delta-i-1,delta+j) = Image(Image.nx()-delta+i,j);
		// filling left side
		for (int i=0;i<Image.nx();i++)
			for (int j=0;j<delta;j++)
				full_image(delta+i,j) = Image(i,delta-j-1);
		// filling bottom side
		for (int i=0;i<Image.nx();i++)
			for (int j=0;j<delta;j++)
				full_image(delta+i,Image.ny()+2*delta-j-1) = Image(i,Image.ny()-delta+j);
		// Computing number of patches to be extracted
		for (int i=0;i<full_image.nx()-delta;i+=OverlapNumber)
			for (int j=0;j<full_image.ny()-delta;j+=OverlapNumber)
				Npatch++;
		// Reading patches, where (i,j) is used as the top left pixel of each patch
		stacked_patches.alloc(patch_size*patch_size,Npatch);
		patchInd = 0;
		for (int i=0;i<full_image.nx()-delta;i+=OverlapNumber)
		{
			for (int j=0;j<full_image.ny()-delta;j+=OverlapNumber)
			{
				// stacking patch (i,j)
				stackInd = 0;
				for (int q=0;q<patch_size;q++)
					for (int p=0;p<patch_size;p++)
					{
						stacked_patches(stackInd,patchInd) = full_image(i+p,j+q);

						stackInd++;
					}
				patchInd++;
			}
		}
	}
	return stacked_patches;
}
dblarray C_DENOISE_IMAGE::denoise_image(dblarray &stacked_patches, dblarray &Dico, int OverlapNumber,int SparsityTarget,double ErrorTarget,intarray full_image_size,bool Verb)
{
	dblarray X; // sparse coding coefficients matrix
	dblarray weights; // weights use to average patches. The weight of a pixel in the image is the number of patches where it is present.
	dblarray denoised_patches(stacked_patches.nx(),stacked_patches.ny());
	int Npatch = stacked_patches.ny(); // number of stacked patches
	dblarray patch_mean(Npatch); // mean of all patches
	dblarray temp_patch(Npix); // temporary patch, used for mean removal
	int patch_size = sqrt(Npix); // length of each patch
	int patchInd, stackInd; // indices used for reading stacked patches
	int delta = patch_size - 1; // size of the band to around Image to build full_image
	int nx = full_image_size(0);
	int ny = full_image_size(1);
	dblarray full_denoised_image(nx,ny);
	dblarray denoised_image(nx - 2*delta,ny - 2*delta);

	//dblarray Image_out
	// Builder coder object
	C_OMP coder(Dico);
	// Removing mean from each patch
	for (int i=0;i<Npatch;i++)
	{
		// loading current patch
		for (int j=0;j<Npix;j++)
			temp_patch(j) = stacked_patches(j,i);

		patch_mean(i) = temp_patch.mean();
		// Removing mean
		for (int j=0;j<Npix;j++)
			stacked_patches(j,i) = stacked_patches(j,i) - patch_mean(i);
	}
	// Sparse coding stacked patches in the dictionary
	// X.alloc(Dico.ny(),Npatch);


	X = coder.omp(stacked_patches,SparsityTarget,ErrorTarget,Verb);
	// Denoising stacked patches by computing DX
	for (int p=0;p<Npix;p++)
		for (int q=0;q<Npatch;q++)
		{
			denoised_patches(p,q) = 0;
			for (int k=0;k<Dico.ny();k++)
				denoised_patches(p,q) += Dico(p,k)*X(k,q);
		}

	// Initializing weights to 0
	weights.alloc(nx,ny);
	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			weights(i,j) = 0;
	// Averaging patches back together to build an image
	patchInd = 0;
	// Initializing denoised image with 0
	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			full_denoised_image(i,j) = 0;
	for (int i=0;i<nx-delta;i+=OverlapNumber)
	{
		for (int j=0;j<ny-delta;j+=OverlapNumber)
		{
			// Reading patch (i,j)
			stackInd = 0;
			for (int q=0;q<patch_size;q++)
				for (int p=0;p<patch_size;p++)
				{
					// full_denoised_image(i+p,j+q) += denoised_patches(stackInd,patchInd);
					full_denoised_image(i+p,j+q) += denoised_patches(stackInd,patchInd) + patch_mean(patchInd);
					weights(i+p,j+q) = weights(i+p,j+q) + 1;
					stackInd++;
				}
			patchInd++;
		}
	}
	// Extracting the center of full_image to remove the symetric boundaries condition pixels
	for (int i=0;i<denoised_image.nx();i++)
		for (int j=0;j<denoised_image.nx();j++)
			denoised_image(i,j) = full_denoised_image(delta+i,delta+j)/weights(delta+i,delta+j);
	return denoised_image;

}












