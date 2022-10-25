#include "C_DL1D.h"
#include <vector>
#include <MatrixOper.h>
#include "IM_IO.h"
#include "C_OMP.h"
using namespace std;
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <fcntl.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_linalg.h>
#include <armadillo>
using namespace arma;
dblarray TabTrain;


C_DL1D::C_DL1D(dblarray &training_set)
{
	TabTrain = training_set; // training set, with samples as column
	Npix = TabTrain.nx(); // atom/training sample length
	eps = 2.220446049250313e-16;
}
C_DL1D::~C_DL1D()
{
}



dblarray C_DL1D::dl1d(dblarray &training_set, dblarray &initD, int IterationNumber,int SparsityTarget,double ErrorTarget,bool Verb)
{
	Na = initD.ny(); // number of atoms in dictionary
	int Ntrain = TabTrain.ny(); // number of training samples
	dblarray sample(Npix); // training sample
	double Amean,sample_mean; // atom and training sample mean
	dblarray D = initD; // learned dictionary initialized with initD
	dblarray X(Na,Ntrain);	 // sparse coding coefficients of TabTrain
	dblarray Xt(Ntrain,Na); // transpose of X
	dblarray Xtpinv(Na,Ntrain); // pseudoinverse of Xt
	dblarray Xpinv(Ntrain,Na); // X pseudo-inverse
	dblarray Scurrent(Npix); // Current sample
	dblarray Acurrent(Npix); // Current atom, during normalization
	double Anorm; // Current atom norm
	double minT; // eigenvalue threshold used to compute pseudo inverse
	int atom_usage; // number of sample using a given atom
	int new_sample_ind; // index of the new sample used in place of an useless atom
	int atoms_replaced; // number of unused atoms replaced with training samples
	MatOper terminator;
	double average_sparsity = 0;
	double average_error = 0;
	dblarray sparse_approx(Npix,Ntrain); // sparse approximation of TabTrain in the current dictionary
	dblarray sparse_sample(Npix); // sparse approximation of current training sample
	mat armaXt(Ntrain,Na);
	dblarray At(Na,Ntrain);
	dblarray Vt(Na,Na);
	dblarray V(Na,Na);
	dblarray S(Na,Na);
	mat armaInv;
	gsl_rng *rng;
	const gsl_rng_type * T;
	T = gsl_rng_default;
	rng = gsl_rng_alloc (T);
	/*gsl_matrix * gslA = gsl_matrix_alloc(Ntrain,Na);
		gsl_matrix * gslV = gsl_matrix_alloc(Na,Na);
		gsl_vector * gslS = gsl_vector_alloc(Na);
		gsl_vector * work = gsl_vector_alloc(Na);*/

	//Removing atoms mean and normalizing atoms with norm > 1
	for (int m=0;m<Na;m++)
	{
		for (int p=0;p<Npix;p++)
			Acurrent(p) = D(p,m); // reading atom p
		Amean = Acurrent.mean();
		for (int p=0;p<Npix;p++) // computing atom mean
			Acurrent(p) = Acurrent(p) - Amean; // removing atom mean
		Anorm = sqrt(Acurrent.energy()); // computing 0-mean atom norm
		if (Anorm > 1)
		{
			for (int p=0;p<Npix;p++)
				D(p,m) = (D(p,m) - Amean) / Anorm; // removing mean and normalizing
		}
		else 	for (int p=0;p<Npix;p++) D(p,m) = D(p,m) - Amean; // only removing mean
	}
	C_OMP coder(D); // building OMP sparse coder
	// Removing training sample mean
	for (int k=0; k<Ntrain; k++)
	{
		for (int i=0; i<Npix; i++) // reading sample k
			sample(i) = TabTrain(i,k);
		sample_mean = sample.mean(); // computing sample mean
		for (int i=0; i<Npix; i++)
			TabTrain(i,k) = sample(i) - sample_mean;  // subtracting mean
	}
	// iterating sparse coding and dictionary update steps
	if (Verb == true)
	{
		cout << "Learning dictionary of " << Na << " atoms of " << Npix << " pixels" << endl;
		cout << "Starting Dictionary Learning for " << IterationNumber << " iterations" << endl;
		cout << "Using OMP with SparsityTarget = "<< SparsityTarget << " and ErrorTarget = " << ErrorTarget << endl;
	}
	for (int i=0;i<IterationNumber;i++)
	{
		if (Verb == true)
		{
			if (IterationNumber > 1000)
			{
				if ((i+1)%(IterationNumber/10) == 1)
					cout << "DL iteration " << i+1 << " / " << IterationNumber << ", average sparsity " << average_sparsity << ", average error " << average_error << endl;
			}
			else cout << "DL iteration " << i+1 << " / " << IterationNumber << ", average sparsity " << average_sparsity << ", average error " << average_error << endl;
		}
		// sparse coding training sample in dictionary
		X = coder.omp(TabTrain,SparsityTarget,ErrorTarget,False);
		if (Verb == true)
		{
			// Computing average sparsity given sparse encoding coefficients
			average_sparsity = 0;
			for (int i=0;i<X.nx();i++)
				for (int j=0;j<X.ny();j++)
					if (X(i,j) !=0)
						average_sparsity++;
			average_sparsity = average_sparsity / (Ntrain);
		}

		// Computing sparse coefficients matrix pseudo inverse for dictionary update
		if (Verb == True)
		{
			cout << "Sparse coding complete, average sparsity " << average_sparsity << endl;
			cout << "Updating dictionary ..." << endl;
		}
		//		terminator.inv_mat_svd(X,Xpinv,minT);  previous method, too slow for large data
		// Transposing matrix X before computing its pseudoinverse
		for (int p=0;p<X.nx();p++)
			for (int q=0;q<X.ny();q++)
				Xt(q,p) = X(p,q);
		// Filling matrix armaXt with coefficients from Xt
		for (int p=0;p<Ntrain;p++)
			for (int q=0;q<Na;q++)
				armaXt(p,q) = Xt(p,q);
		//Chosing eigenvalue threhsold value
		minT = Ntrain*eps * sqrt(Xt.energy());
		// Computing pseudoinverse
		armaInv = pinv(armaXt,minT);
		for (int p=0;p<Na;p++)
			for (int q=0;q<Ntrain;q++)
				Xtpinv(p,q) = armaInv(p,q);

		// gsl_matrix_set(gslA,p,q,Xt(p,q));
		// gsl_matrix_free(gslA);
		// gslA = gsl_matrix_alloc(Ntrain,Na);
		// for (int p=0;p<Ntrain;p++)
		// 	for (int q=0;q<Na;q++)
		// 		gsl_matrix_set(gslA,p,q,Xt(p,q));
		// Computing SVD from A
		// gsl_linalg_SV_decomp (gslA,gslV,gslS,work);
		// gsl_linalg_SV_decomp_jacobi (gslA,gslV,gslS);
		// Thresholding/inverting eigenvalues
		// S.init(0);
		/*for (int q=0;q<Na;q++)
			if (gsl_vector_get(gslS,q)<minT)
				S(q,q) = 0;
			else
				S(q,q) = 1/(gsl_vector_get(gslS,q));
		// Reading remaining gsl matrices
		for (int p=0;p<Ntrain;p++)
			for (int q=0;q<Na;q++)
				At(q,p) = gsl_matrix_get(gslA,p,q);
		for (int p=0;p<Na;p++)
			for (int q=0;q<Na;q++)
				V(p,q) = gsl_matrix_get(gslV,p,q);
		// Computing pseudo inverse by multiplying matrices
		Xtpinv = mult(V,mult(S,At));*/

		// Applying MOD dictionary update
		for (int p=0;p<D.nx();p++)
			for (int q=0;q<D.ny();q++)
			{
				D(p,q) = 0;
				for (int k=0;k<Ntrain;k++)
					D(p,q) += TabTrain(p,k)*Xtpinv(q,k);
			}
		// Computing sparse approximation and average quadratic error
		average_error = 0;
		sparse_approx = mult(D,X);
		for (int k=0;k<Ntrain;k++)
		{
			for (int np=0;np<Npix;np++)
				sparse_sample(np) = sparse_approx(np,k) - TabTrain(np,k);
			average_error =+ sqrt(sparse_sample.energy())/Ntrain;
		}
		// Throwing away unused atoms and replacing them by random training samples
		atoms_replaced = 0;
		for (int m=0;m<Na;m++)
		{
			atom_usage = 0;
			for (int p=0;p<Ntrain;p++)
				if (X(m,p)!=0) atom_usage++;
			if (atom_usage == 0)
			{
				new_sample_ind	= gsl_rng_uniform_int(rng, Ntrain-1);
				for (int k=0;k<Npix;k++)
					Acurrent(k) = TabTrain(k,new_sample_ind);
				Amean = Acurrent.mean();
				for (int k=0;k<Npix;k++)
					Acurrent(k) = Acurrent(k) - Amean;
				Anorm = sqrt(Acurrent.energy()); // computing 0-mean atom norm
				if (Anorm > 1)
				{
					for (int k=0;k<Npix;k++)
						D(k,m) = Acurrent(k) / Anorm; // removing mean and normalizing
				}
				else
					for (int k=0;k<Npix;k++)
						D(k,m) = Acurrent(k); // only removing mean
				atoms_replaced++;
			}
		}
		if (Verb == True && atoms_replaced !=0)
			if (IterationNumber > 0)
				if ((i+1)%(IterationNumber/10) == 1) cout << "Replaced " << atoms_replaced << " unused atoms with random training samples" << endl;
				else
					cout << "Replaced " << atoms_replaced << " unused atoms with random training samples" << endl;

		// Normalizing atoms with norm > 1
		for (int m=0;m<Na;m++)
		{
			for (int p=0;p<Npix;p++)
				Acurrent(p) = D(p,m);
			Anorm = sqrt(Acurrent.energy());
			if (Anorm > 1)
				for (int p=0;p<Npix;p++)
					D(p,m) = D(p,m) / Anorm;
		}



		// Updating sparse coder with new version of dictionary
		coder.update_dictionary(D);
	}
	return D;
}











