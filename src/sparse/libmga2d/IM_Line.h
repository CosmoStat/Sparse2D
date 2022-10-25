/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/09/02 
**    
**    File:  IM_Line.h
**
*******************************************************************************
**
**    DESCRIPTION:  Line detection in an image
**    ----------- 
**
******************************************************************************/

#ifndef __IM_LINE__
#define __IM_LINE__

#define MAXLINE 512 // Maximum line size

void getline (int x1,int y1,int x2,int y2, intarray & TX, intarray & TY, int & N);
// get the pixel coordinates of a digital line segment from (x1,y1) to (x2,y2)
// x1: IN = x-coordinate of the first pixel
// y1: IN = y-coordinate of the first pixel
// x2: IN = x-coordinate of the last pixel
// y2: IN = y-coordinate of the last pixel
// N:  OUT = number of pixel in the line segment
// TX: OUT = x-coordinate array => TX(0..N-1)   
// TY: OUT = x-coordinate array => TY(0..N-1)
// WARNING:  TX,TY must be allocated before the call with a size
//           larger than N

void symgetline (int x2, int y2, int Nl2, int Nc2, intarray & TX, intarray & TY, int & Np);
// get the pixel coordinates of a digital line segment from (-x2,-y2) to (x2,y2) in an
// image of size Nl2, Nc2. The line is symmetric around the pixel (0,0)
// x2: IN = x-coordinate of the first pixel
// y2: IN = y-coordinate of the first pixel
// Np:  OUT = number of pixel in the line segment
// TX: OUT = x-coordinate array => TX(0..Np-1)   
// TY: OUT = x-coordinate array => TY(0..Np-1)
// WARNING:  TX,TY must be allocated before the call with a size
//           larger than N
//    TX(0) = x2 and TX(Np-1) = -x2
//    TY(0) = y2 and TY(Np-1) = -y2


/******************************************************************************/
// Atom definition
// An atom is defined by a list of pixel positions
class Atom 
{
     int Np;   // Number of pixels
     intarray TabX; // TabX(0..Np-1) = X position of the pixels
     intarray TabY; // TabY(0..Np-1) = Y position of the pixels
     
   public:
     void setpos(int Num, int X, int Y) { TabX(Num) = X; TabY(Num) = Y;}
              // Set the pixel number Num at the position
	      // defined by X and Y
     int np () { return Np;}
             // return the number of pixel positions
     Atom(){}
     void alloc(int N) {Np = N; TabX.alloc(N); TabY.alloc(N);}
            // Allocate an atom with N position
     int  xpos(int Num) { return TabX(Num);}
            // the x position of the pixel number Num
     int  ypos(int Num) { return TabY(Num);}
            // the y position of the pixel number Num
      ~Atom(){}
};

/****************************************************************************/

// This class contains a list of atoms and the operators
// to be applied to an image.
// Atoms are hard coded, they are defined by a set of lines
// of different sizes and different orientations
class ListAtom
{
   int Na;         // Number of atoms
   Atom *TabAtom;  // Array of atoms
   public: 
     Bool Verbose;
     Bool Norm;   // If Norm == True, the correlation values are normalized
                  // by sqrt(Np), where Np is the number of pixel in the atom.
                  // Default is true.
     int Step;  // distance between two pixels (default is 1).
     int na () { return Na; } // return the number of atoms
     Atom & at(int NumAtom) {return TabAtom[NumAtom];}
                              // the atom defined by the number NumAtom
     void get_best_correl(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma);
          // Correllation of all atoms at all positions of the image Ima,
	  // best correlation values are stored in Conv, and 
	  // best atoms in AtomIma.
	  // Ima -- Input image
	  // Conv -- Output: Conv(i,j) = best correlation value at pixel
	  //         position (i,j).
	  //          Conv(i,j) is normalized by sqrt(Np) where Np is
	  //          the number of pixels in the best atom.
	  // AtomIma -- Output: AtomIma(i,j) is the best atom at pixel
	  //         position (i,j).
     float atom_pos_correl(Ifloat &Ima, int AtNumber, int i, int j);
          // return the correlation value of the atom AtNumber
	  // with the image at position (i,j).
     void order_coef(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                          unsigned long * & TabInd);
          // ordered the coefficients Conv(i,j)
	  // Ima: InOut -- Selected atoms are subtracted from the images.
	  // Conv: InOut -- Conv is reevaluated by substraction from
	  //                the image the selected atoms, starting from
	  //                the highest fit.
	  // AtomIma: In -- AtomIma(i,j) is the best atom at pixel
	  //         position (i,j).               
          // TabInd: Out -- Ordering of the atoms following their fit.
	  //         TabInd[N-i-1] = ith best correlation at position
	  //                     IM = TabInd[N] / Nl
	  //                     JM = TabInd[N] % Nl
      
      inline void get_index_pos(int Ind, unsigned long * & TabInd, int Nl, 
                           int & IM, int &JM)
      {
         int MaxInd = TabInd[Ind] - 1;
         IM = MaxInd/ Nl;
	 JM = MaxInd - IM*Nl;
      }
      
      void add_atom(Ifloat &Ima, int i, int j, float ConvVal, int AtNumber );
      // add the image at the position (i,j), the atom AtNumber with an 
      // amplitude ConvVal.

  
      void add_nfirst_atoms(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                        unsigned long * & TabInd, int NFirst );   
      // add the image the Nfirst atoms.
      // If Norm == True, the values in Conv are  normalized by sqrt(Np),
      // where Np is the number of pixel in the atom.  

      void add_nsigma_atoms(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                        unsigned long * & TabInd, float NSigma, 
			float SigmaNoise);          
     
     void print()
     {
        cout << "Number of atoms = " << na() << endl;
        for (int a=0; a < na(); a++)
        {
           cout << "ATOM " << a+1 << " Np = " << TabAtom[a].np() << endl;
	   for (int i=0; i < TabAtom[a].np(); i++)
           cout << "X = " << TabAtom[a].xpos(i) << " Y = " << TabAtom[a].ypos(i) << endl;
        }
    }
         
    ListAtom(); 
    ~ListAtom() {};
};

/******************************************************************************/
//
//              LINE DETECTION IN AN IMAGE
// Example:
//   ListAtom LA;
//      LA.get_best_correl(Band,Conv,AtomIma);
//      LA.order_coef(Band, Conv, AtomIma, TabInd);
//      LA.add_nsigma_atoms(BandRec, Conv,  AtomIma, TabInd, 10., SigmaNoise);
//      
// 
/******************************************************************************/
    
#endif
