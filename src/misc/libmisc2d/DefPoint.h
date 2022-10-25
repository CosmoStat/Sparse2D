/***********************************************************
**	Copyright (C) 1999 CEA
************************************************************
**
**    UNIT
**
**    Version:  1.0
**
**    Author: 	J.L. Starck
**
**    Date: 	27/08/99
**    
**    File:  	DefPoint.h
**
************************************************************
**
**  Point 1D,2D,3D definition 
**  Array of point Definition
**  
************************************************************/


#ifndef	_DEFPOINT_H_
#define	_DEFPOINT_H_

// #include "IM_Math.h"
#include "Array.h"
#define D2R (M_PI/180.0)
#define NBR_TYPE_COORD 2
#define TCOORD_XYZ 1
#define TCOORD_LON_LAT 2
 
#define NBR_RND_CAT 5 
//enum type_random_cat {RND_CAT_USER, RND_CAT_LAMBDA_CDM,
//                      RND_CAT_IRAS, RND_CAT_IRAS_NORTH, RND_CAT_IRAS_SOUTH, 
//		      RND_CAT_UNDEFINED=-1};

inline char * StringCoord (int type)
{
    switch (type)
    {
        case TCOORD_XYZ: 
			return ((char*) "XYZ coordinate");break;
        case TCOORD_LON_LAT: 
			return ((char*) "longitude-latitude coordinate");break;
		default:
			return ((char*) "Undefined coordinate type");
			break;
    }
}

/* inline char * StringRNDCat (type_random_cat type)
{
    switch (type)
    {
        case RND_CAT_USER: 
              return ("User defined in the header of the catalog");break;
        case RND_CAT_LAMBDA_CDM: 
              return ("Lambda CDM simulation");break;
        case RND_CAT_IRAS_NORTH: 
              return ("IRAS North 1.2 Jy");break;
        case RND_CAT_IRAS_SOUTH: 
              return ("IRAS South 1.2 Jy");break;
        case RND_CAT_IRAS: 
              return ("IRAS North and South 1.2 Jy");break;
	default:
              return ("Undefined random catalogue type");
              break;
    }
}
*/

// Point definition
class Point {
     int Dim;        // Point dimension
     fltarray Coord; // Coordinate array
    public:
    Point() {Dim=0;}
    Point(int Dimension) {alloc(Dimension);}
    void alloc(int Dimension) {Dim=Dimension; Coord.alloc(Dim);}
    int dim () const  {return Dim;}       // return the dimension
    float & x() const { return Coord(0);} // return first coordinate
    float & y() const { return Coord(1);} // return second coordinate
    float & z() const { return Coord(2);} // return third coordinate
    float & axis(int i) const { return Coord(i);} 
                                   // return the ith coordinate
				   // i = 0 .. Dim-1
    //  definition of the "=" operator
    const Point & operator = (const Point & P)
                  { for (int i=0; i < Dim; i++) Coord(i) = P.axis(i);
		    return *this;}
    void random (float Min=0., float Max=1.); // create a random point
                                              // with all coordinate between
					      // Min and Max
    void random(Point & PMin, Point & PMax);  // idem, coordinate between
                                              // Pmin abd Pmax
    void read (FILE *File);  // read a point from a file
                             // number of float values read = Dim
    void write (FILE *File); 			     
    void print () const;     // print to std the point	   
    void trace () const;     // print to std the point	   
};

// return the distance between two points
float dist(const Point &P1, const Point &P2);
float sphdist(const Point &P1, const Point &P2);
// return the (square of) distance between two points
float squaredist(const Point &P1, const Point &P2);
float squaredist(const Point &P1, const Point &P2,float SquareDistMax);
float squaresphdist(const Point &P1, const Point &P2);
// Array of point definition

class ArrayPoint {
   Point *TabPoint; // Array of points
   int Np;          // Number of points
   int Dim;         // Dimension space 1,2 or 3
   public:
    int TCoord; // coordinate system
    int BootCoord[3]; // BootCoord[i] equal 1 if the coodinate must be 
                     // bootstraped, and 0 otherwise
    Point Pmin;    // minimum of the array point
    Point Pmax;    // maximum of the array point
    
    ArrayPoint(){Np=0;Dim=0;TabPoint=NULL;TCoord=0;BootCoord[0]=0;
                 BootCoord[1]=0;BootCoord[2]=0;}
    ArrayPoint(int Dimension, int N);
    void alloc(int Dimension, int N);
    
    //  definition of the "=" operator
    const ArrayPoint & operator = (const  ArrayPoint &Tab)
                  { for (int i=0; i < Np; i++) TabPoint[i] = Tab(i);
		    Pmin = Tab.Pmin; Pmax = Tab.Pmax;
		    return *this;}	
    // return a point i=0..N-1   	    
    inline Point & operator () (int i)  const { return  TabPoint[i];}	
    int dim () const  {return Dim;}  // return the dimension
    int coord() const  {return TCoord;}  // return the coordinate type
    int np() const { return Np;}     // return the number of points
    void print(char *Mes =NULL);     // print the full array to stdout
    
    void random(ArrayPoint & Data);
    // creates a random catalogue: The array must be first allocated.
        
    void write(char *FileName, Bool Verbose=False);
    void read(char *FileName, Bool Verbose=False, Bool ReadSimu=False);
    // read an array point from a file
    // Data format = Dim NumberofPoints CoordinateType
    //               MinAxis1 MaxAxis1 BootAxisi
    //               ...
    //               MinAxisi MaxAxisi BootAxisi
    //               coordinate point 1
    //               ...
    //               coordinate point N
    

    void toxyz();
    // convert to rectangular coordinates

    void minmax(Point & PMi, Point & PMa);
                                            // return the min and max of the 
					    // array
    ~ArrayPoint() { if (TabPoint != NULL) delete [] TabPoint;Dim=Np=0;}
    void trace () const;     // print to std the point	   

};



void bootstrap(ArrayPoint & Data, ArrayPoint & BootStrapData);
// make a boot strap on data and store the result in BootStrapData
void bootstrapxyz(ArrayPoint & Data, ArrayPoint & BootStrapData);
// idem but the bootstrap is done separately on each axis.

#endif


