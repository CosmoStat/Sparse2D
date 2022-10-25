/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/12
**    
**    File:  CPca.cc
**
*******************************************************************************
**
**    DESCRIPTION  class for principal component analysis  
**    ----------- 
**                 
******************************************************************************/
 
#include "GlobalInc.h"
#include "CPca.h"
#include "NR.h"

/************************************************************/
float *nr_vector(int nl,int nh)
{
    float *v;

    v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}
/************************************************************/

void diag_matrix(fltarray & CorrelMat, fltarray & DiagCorMat,
                 fltarray & EigenValue, fltarray & EigenVector,
                 int & nrot, Bool SortEigen)
{
       int i,j;
       float **c; //intercorrelation matrix
       float **v; //intercorrelation eigenvectors matrix 
       float *d;  //intercorrelation eigenvalues vector 
       int N_Elem = CorrelMat.nx();
       
       // allocate intercorrelation matrix 
       c=matrix(1,N_Elem,1,N_Elem);   
       for(i=1; i <= N_Elem; i++)
       for(j=i; j <= N_Elem; j++)  {
          c[i][j]=c[j][i]= CorrelMat(j-1,i-1);
          //cout << i << "," << j << "," << c[j][i] << endl;
       }
        
       // allocate output eigenvectors 
       v=matrix(1, N_Elem, 1, N_Elem);
 
      // allocate output eigenvalues 
      d= (float *) nr_vector(1, N_Elem);

      // compute eigenvectors and eigenvalues of c
      jacobi(c, N_Elem, d, v, &nrot);
      
       // sort eigenvalues 
      if (SortEigen == True) eigsrt(d,v,N_Elem);
           
      for(i=0; i < N_Elem; i++) EigenValue(i) = d[i+1];
      for(i=0; i < N_Elem; i++)
      for(j=0; j < N_Elem; j++)
      {
         EigenVector(i,j) = v[i+1][j+1];
         DiagCorMat(j,i) = c[i+1][j+1];
      }
}


/*************************************************/                 

void diag_verif(fltarray &CorrelMat, int N_Elem,  fltarray & EigenValue,  
            fltarray & EigenVector)
{ 
   fltarray  y(N_Elem);
   int i,j,Vect;

   // apply the matrix on each eigen vector
   for (Vect=0; Vect < N_Elem; Vect++)
   {      
      // matrix operation
      for (i=0; i < N_Elem; i++)
      {
	  y(i)=0.;
	  for (j=0; j < N_Elem; j++) y(i)+= CorrelMat(j,i)*EigenVector(j,Vect);
 	  if (EigenValue(Vect) !=0.) y(i) /=  EigenValue(Vect);
	  else y(i)=0;
      }
      
      cout << "-------------------------------------" << endl; ;
      cout << "EigenValue[" << Vect+1 << "]= " <<  EigenValue(Vect) << endl; ;
      cout << "EigenVector[" << Vect+1 << "]= " ;
      for (i=0;i<N_Elem;i++) cout << EigenVector(i,Vect) << " ";
      cout << endl;
      
      cout << "EigenVerif[" << Vect+1 << "]= " ;
      for (i=0;i<N_Elem;i++)  cout << y(i) << " ";
      cout << endl;
    }

    for (i=0; i < N_Elem; i++)
    {
       y(i)=0.;
       for (j=0; j < N_Elem; j++)
           y(i)+= EigenVector(i,j)* EigenVector(i,j)*EigenValue(j);
       cout << " j = " << i+1 << ": sum_i V^2(i,j) Eigen(j) = " << y(i) << endl;
    }    
}


/*************************************************/                 


void CPCA::alloc(int Dim)
{
      N_Elem = Dim;
      CorrelMat.alloc(N_Elem, N_Elem);
      DiagCorMat.alloc(N_Elem, N_Elem);
      EigenVector.alloc(N_Elem, N_Elem);
      EigenValue.alloc(N_Elem);
      Verbose = False;
      NormAna = False;
      SortEigen =True;
}

/*************************************************/                 

CPCA::CPCA(int Dim)
{
   alloc(Dim);
}

/************************************************************/

void CPCA::print ()
{
   int i,j;
      
   // print intercorrelation matrix
   cout << endl << "correlation matrix" << endl;
   for (i=0; i< N_Elem; i++)
   {
 	 for (j=0; j < N_Elem; j++) 
	 {
	    cout.width(8);  // see p343 c++ stroustrup
	    // cout.fill(' ');
	    cout.setf(ios::right,ios::adjustfield);
	    cout.setf(ios::fixed,ios::floatfield);
	    cout.precision(2);
	    cout << CorrelMat(j,i) << " " ;
	 }
	 cout << endl;
   }
      
 
   for (int Vect=0; Vect < N_Elem; Vect++)
   {
      cout << "-------------------------------------" << endl; ;
      cout << "EigenValue[" << Vect+1 << "]= " <<  EigenValue(Vect) << endl; ;
      cout << "EigenVector[" << Vect+1 << "]= " ;
      for (i=0;i<N_Elem;i++) cout << EigenVector(i,Vect) << " ";
      cout << endl;
   }
   cout.precision(6);
// #ifndef WINDOWS
//   cout.setf(0,ios::floatfield); // reset to default
// #endif
}

/************************************************************/

void CPCA::compute_eigen(fltarray & CorrelMatrix)
{
    int nrot;
       
    CorrelMat = CorrelMatrix;
 
    diag_matrix(CorrelMat, DiagCorMat, EigenValue, EigenVector, nrot,SortEigen);
    if (Verbose == True)
    {
        // cout << "Jacobi nrot: " << nrot << endl;
    //    print();
        // cout << "verification" << endl;
        // diag_verif(CorrelMat , N_Elem,  EigenValue,  EigenVector);
    }
}       

/************************************************************/

void CPCA::transform (fltarray &Vin, fltarray &Vout)
{
     int i,j;
     for (i = 0; i < N_Elem; i++)
     {
 	 Vout (i)=0.;
 	 if (NormAna == True)
 	 {
	    for (j=0; j <  N_Elem; j++) 
	    {
	       if (EigenValue(j) > 0)
	          Vout(i) += EigenVector(j,i)* Vin(j)  / sqrt(EigenValue(j));
	    }
         }
         else
 	 {
	    for (j=0; j <  N_Elem; j++)  {
            Vout(i) += EigenVector(j,i)* Vin(j); 
            } 
         }
     }        
}


void CPCA::transform_noise (fltarray &Vin, fltarray &Vout) {
     
   int i,j;
   for (i = 0; i < N_Elem; i++) {
      Vout (i)=0.;
      if (NormAna == True) {
         for (j=0; j <  N_Elem; j++) {
	    if (EigenValue(j) > 0)
	       Vout(i) += pow((double)  EigenVector(j,i),2.) * Vin(j)  / (EigenValue(j));
	 }
      } else {
	 for (j=0; j <  N_Elem; j++)  {
            Vout(i) += pow((double) EigenVector(j,i),2.) * Vin(j); 
         } 
      }
   }        
} 
/*************************************************/
  
void CPCA::invtransform (fltarray &Vin, fltarray &Vout, int NbrEigenVect)
{
    int i,j;
    int N = (NbrEigenVect > 0) ? NbrEigenVect : N_Elem;
    if (N > N_Elem) N = N_Elem;
    
     for (i = 0; i < N_Elem; i++)
     {
 	Vout (i)=0.;
 	for (j=0; j < N; j++)  {
	   Vout(i) += EigenVector(i, j)* Vin(j);
	}
  	if (NormAna == True)   Vout(i) *= sqrt(EigenValue(i));
     }        
}

/*************************************************/


