Dictionary Learning Package
Simon Beckouche
simon@beckouche.fr or simon.beckouche@cea.fr
Service d'Astrophysique, CEA Saclay, April 2014



The packages contains the following executables: 
+ dl_omp            (compute a sparse representation of a vector in an overcomplete dictionary using the Orthogonal Matching Pursuit algorithm)
+ dl_dl1d           (learn a dictionary from 1d or 2d training set using the Method of Optimal Direction alrogithm)
+ dl_denoise_image  (extract, denoise and average patches to denoise an image)
+ dl_extract_patches (extract patches at random position in an image)


0) Quick use
If you already know or do not care about the details, here is a list of command to use to run the functions of this package. Place a noicy image Image.fits in the same directory as the executables. You get a file denoisedImage.fits as output. 
./dl_extract_patches -v -W 10 -N 200 Image.fits inDico.fits  [extracting from the noisy image a set of 200 patches of 10 by 10 pixels to be used as initial dictionary]
./dl_extract_patches -v -W 10 -N 2000 Image.fits training.fits [extracting from the noisy image a set of 5000 patches of 10 by 10 pixels used as training set]
./dl_dl1d -E 330625 -v -I 20 inDico.fits training.fits learnedDico.fits [performing 20 iteration of dictionary learning. The value for the -E parameter should be (1.15*sigma*patch_width)^2]
./dl_denoise_image -v -E 330625 -Q 1 learnedDico.fits Image.fits denoisedImage.fits [denoising the noisy image by denoising and averaging every overlaping patches in the image]


1) dl_omp usage
dl_omp computes a sparse representation of a vector in an overcomplete dictionary using the Orthogonal Matching Pursuit algorithm.
You can run dl_omp using the syntax
'dl_omp [-v] [-S SparsityTarget or -E ErrorTarget] dictionary sample coefficient'

+ dictionary : fits file containing the dictionary used to compute the sparse decomposition with atoms as columns. It should at least contain as many atoms as the dimension of each atom to ensure stable results.
+ sample : fits file containing the sample which you want to sparse code in the dictionary, placed as columns. Each sample should be of the same length than dictionary atoms.
+ coefficient : fits file containing the computed sparse coefficient of each sample placed as column.

You can use the following options :
+ "-v" enables verbose mode.
+ "-S SparsityTarget" specifies the maximum number of nonzero coefficients to estimate (integer between 0 and the dimension of a sample).
+ "-E ErrorTarget" specifies the accuracy below which the sparse coding stops (real positive number). When denoising a sample of size n contaminated with a white gaussian noise of std sigma, you should use ErrorTarget = (C*sigma)^2 * n where we usually take C = 1.15.
Default values are without verbose mode, and using SparsityTarget = 5.

A few examples :
dl_omp -v -E 0.1 Dico.fits Signal.fits Coeff.fits
dl_omp -v -S 3 Dico.fits Signal.fits Coeff.fits

2) dl_1d
dl_1d learns a sparsifying overcomplete dictionary from 1d data using the MOD algorithm. It can be used on 2d data by stacking 2d samples into 1d samples. It ensures all atoms are used by replacing unused atoms with random traning sample during the learning process. 
You can run dl_dl1d using the syntax
'dl_dl1d [-v] [-S SparsityTarget or -E ErrorTarget] [-I IterationNumber] initial_dictionary training_set learned_dictionary'

+ initial_dictionary : fits file containing the dictionary used as initialization of the learning process. It should at least contain as many atoms as the dimension of each atom to ensure stable results.
If you input a dictionary with more rows than columns, it will considered that atoms are placed as columns.
+ training_set : fits file containing training samples placed as columns. Each sample should be of the same length than initial_dictionary atoms. Note that the mean of each training sample is removed before learning. 
+ learned_dictionary : fits file containing the learned dictionary. The learned dictionary is of the same dimensions than initial_dictionary, and has atoms which norm doesn't exceed 1. 

You can use the following options :
+ "-v" enables verbose mode.
+ "-S SparsityTarget" specifies the maximum number of nonzero coefficients to estimate during the sparse coding step (similar to the case dl_omp).
+ "-E ErrorTarget" specifies the accuracy below which the sparse coding stops (similar to the case dl_omp). 
+ "-I IterationNumber" specifies the number of iterations for the MOD algorithm (positive integer).
Default values are without verbose mode, using SparsityTarget = 5 and IterationNumber = 20.

A few examples :
dl_dl1d -v -S 3 -I 20 Dictionary.fits TrainingSet.fits LearnedDictionary.fits
dl_dl1d -v -E 0.42 -I 20 Dictionary.fits TrainingSet.fits LearnedDictionary.fits

3) dl_denoise_image
dl_denoise_image denoise and image by extracting a regular grid of patches from the noisy image, denoising them using a sparsifying dictionary, and averaging them back together. Patches don't have to be square.
You can run dl_denoise_image using the syntax
'dl_denoise_image [-v] [-S SparsityTarget or -E ErrorTarget] [-O OverlapNumber] dictionary NoisyImage DenoisedImage'

+ dictionary : fits file containing the dictionary used for denoising with atoms as columns. It should at least contain as many atoms as the dimension of each atom to ensure stable results.
+ NoisyImage : fits file containing the image to denoise
+ DenoisedImage : fits file containing the denoised image

You can use the following options :

+ "-v" enables verbose mode.
+ "-S SparsityTarget" specifies the maximum number of nonzero coefficients to estimate for denoising (similar to the case dl_omp). Using this option is not recommended when processing noisy data. 
+ "-E ErrorTarget" specifies the accuracy below which the sparse coding stops (similar to the case dl_omp). To denoise an image contaminated by a white gaussian noise of std sigma, you should use use ErrorTarget = (C*sigma)^2 * n where we usually take C = 1.15.
+ "-O OverlapNumber" specifies the number of pixels between each extracted patch (integer between 1 and the minimum of both dimensions of patch size). OverlapNumber = 1 means that all patches are extracted and denoised, OverlapNumber = patch_size/2 means half of the patch are extracted and denoised. Larger valuer of this parameter means faster processing, more block artifacts and less blurring.
Default values are without verbose mode, using SparsityTarget = 5 and OverlapNumber = 1. 

A few examples :
dl_denoise_image -v -E 0.42 -Q 1 Dictionary.fits NoisyImage.fits DenoisedImage.fits
dl_denoise_image -v -S 3 -Q 1 Dictionary.fits NoisyImage.fits DenoisedImage.fits

4) dl_extract_patches
dl_extract_patches extracts a given number of square patches of a given width from an image at random position. They can be used to train a dictionary afterwards. Patches are extracted at position picked uniformly at random, with a new random seed at each execution. 
You can run dl_extract_patches using the syntax
'dl_extract_patches [-v] [-W PatchWidth] [-N NumberOfPatches] Image Patches'

+ Image : fits file containing the image. The only restriction is it being larger than the patch width. 
+ Patches : fits file containing an image with as column the extracted patches stacked column by column in 1d vectors to be used with dl1d. 


You can use the following options :

+ "-v" enables verbose mode.
+ "-W PatchWidth" specifies the width (in pixels) of each patch. 
+ "-N NumberOfPatches" specifies how many patches should be extracted.

Default values are without verbose mode.

An example : 
dl_extract_patches -v -W 10 -N 5000 Image.fits Patches.fits

