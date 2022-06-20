#ifndef C_OMP_H
#define C_OMP_H

#include <TempArray.h>

class C_OMP
{
    public:
        C_OMP(dblarray &dictionary);
        ~C_OMP();
        dblarray omp(dblarray &input_vector,int SparsityTarget, double ErrorTarget, bool Verb);
        dblarray omp(dblarray &input_vector,int SparsityTarget, double ErrorTarget){return omp(input_vector,SparsityTarget,ErrorTarget,False);}
        dblarray omp(dblarray &input_vector, double ErrorTarget){return omp(input_vector,Npix,ErrorTarget,False);}
        dblarray omp(dblarray &input_vector, int SparsityTarget){return omp(input_vector,SparsityTarget,-1,False);}
        void update_dictionary(dblarray &new_dictionary);
        void update_scaling(dblarray &rescale);

      int Na;  // number of atoms
      int Npix; // number of pixels per atom
      double eps; // constant used to determine minT to control SVD eigenvalues threshold. Based on Matlab eps.

    protected:
              dblarray TabDico;
              dblarray rescale_factor;
    private:
};

#endif // C_OMP_H
