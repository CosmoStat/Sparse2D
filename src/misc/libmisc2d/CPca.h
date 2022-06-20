/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  97/03/18
**    
**    File:  PCA.cc
**
*******************************************************************************
**
**    DESCRIPTION  do a principal component analysis  
**    ----------- 
**
***************************************************************************/

#ifndef _CPCA_H_
#define _CPCA_H_

// principal component analysis class

class CPCA {
   int N_Elem;           // Number of element in both directions
public:
   fltarray CorrelMat;   // correlation matrix (2D)
                         // CorrelMat(j,i) = index at column j and line i
   fltarray DiagCorMat;  // diagonalized correlation matrix (2D)
   fltarray EigenValue;  // Eigen Values of the matrix  (1D)
   fltarray EigenVector; // Eigen vectors (2D)
                         // EigenVector(*,p) = eigen vector number p
                         //     p >=0 and p < Nbr_Image
                         //    eigen vectors are sorted in reverse order
                         //    following eigen values
   Bool Verbose;         // print some infomation
   Bool NormAna;         // if true, coefficient are divided by the root square
                         // of the eigen values (reduced coef ...)
   Bool SortEigen;       // if true, eigen vector are sorted
   CPCA() { N_Elem = 0; Verbose = True; NormAna=False;}
   CPCA(int Dim); // create the objet with a given dimension
   void init() {N_Elem = 0; Verbose = True; NormAna=False;}
   void alloc (int Dim); // allocate the objet with a given dimension
    
   void print ();        // print the correlation matrix, the eigen values 
                         // and the eigen vectors
  
   void compute_eigen(fltarray & CorrelMatrix); 
                        // diagonalize the correlation matrix and compute
                        // the eigen values and the eigen vector
                        
   void transform (fltarray &Vin, fltarray &Vout);
                       // transform a vector using the matrix defined
                       // by the eigen vector
   
   void transform_noise (fltarray &Vin, fltarray &Vout);
                       // transform a noise vector using the matrix defined
                       // by the eigen vector 
		                                                 
   void invtransform (fltarray &Vin, fltarray &Vout, int NbrEigenVect=-1);
                       // apply the inversed transform
                       // by default all eigen vector are taken
                       // if NbrEigenVect > 0, only the first eigen vectors
                       // are used for the reconstruction
                       
   ~CPCA(){N_Elem = 0; Verbose = False;};
};


void diag_matrix(fltarray & CorrelMat, fltarray & DiagCorMat,
                 fltarray & EigenValue, fltarray & EigenVector,
                 int & nrot, Bool SortEigen=True);
void diag_verif(fltarray &CorrelMat, int N_Elem,  fltarray & EigenValue,  
            fltarray & EigenVector);                 
#endif
