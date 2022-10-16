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
**    File:  	DefPoint.cc
**
************************************************************
**
**  Point 1D,2D,3D definition 
**  Array of point Definition
**  
************************************************************/


#include "DefPoint.h"
    
/*****************************************************************/

void Point::random(Point & PMin, Point & PMax) 
{ 
   for (int i=0; i < Dim; i++)  
      Coord(i) = get_random(PMin.axis(i), PMax.axis(i));
        // Coord(i)= drand48()*(PMax.axis(i)-PMin.axis(i)) + PMin.axis(i);
        // get_random(PMin.axis(i), PMax.axis(i));
}  

/*****************************************************************/

void Point::random(float Min, float Max) 
{ 
   for (int i=0; i < Dim; i++)  Coord(i)= get_random(Min,Max);
    //      Coord(i)= drand48()*(Max-Min) + Min;
   // Coord(i)= get_random(Min,Max);
}   
    
/*****************************************************************/

void Point::print () const 
{
   for (int i=0; i < Dim; i++) cout << Coord(i) << " ";
   cout << endl;
}

/*****************************************************************/

void Point::read (FILE *File)
{
   int i,Nb=0;
   float x,y,z;
   switch (Dim) 
   {
      case 1 : 
         Nb = fscanf(File, "%f\n", &x);
	 Coord(0) = x;
         break;
      case 2 : 
         Nb = fscanf(File, "%f\t%f\n", &x,&y);
         Coord(0) = x;Coord(1) = y;
         break;
      case 3 : 
         Nb = fscanf(File, "%f\t%f\t%f\n", &x,&y,&z);
         Coord(0) = x;Coord(1) = y; Coord(2) = z;
         break;
      default:
         for (i=0; i < Dim; i++) 
	 {
            Nb += fscanf(File, "%f\t",&x);
            Coord(i) = x;
         }
         break;
    }
    if (Nb != Dim)
    {
       cerr << "Error: cannot read data value ... " << endl;
       exit(-1);
    }
}

/*****************************************************************/

void Point::write (FILE *File)
{
   int i,Nb=0;
   switch (Dim) 
   {
      case 1 : 
         Nb = fprintf(File, "%f\n", Coord(0));
         break;
      case 2 : 
         Nb = fprintf(File, "%f\t%f\n", Coord(0), Coord(1));
         break;
      case 3 : 
         Nb = fprintf(File, "%f\t%f\t%f\n", Coord(0), Coord(1), Coord(2));
         break;
      default:
         for (i=0; i < Dim; i++)  Nb += fprintf(File, "%f\t",Coord(i));
         fprintf(File, "\n");
         break;
    }
    if (Nb < 0)
    {
       cerr << "Error: cannot write data value ... " << endl;
       exit(-1);
    }
}
/*****************************************************************/

void Point::trace() const {

   std::cout << "p(" << Dim << ") := ["; 
   for( int i=0; i<Dim-1; i++ )
      std::cout << Coord(i) << ",";
   std::cout << Coord(Dim-1) << "]" << std::endl;
   
}
/*****************************************************************/

float  dist(const Point &P1, const Point &P2)
{
    double Sum=0.;
 
    for (int i=0; i < P1.dim(); i++) 
      Sum += ((double) P1.axis(i)-P2.axis(i)) * ((double) P1.axis(i)-P2.axis(i));
    return ((float) (sqrt(Sum)));
} 

/*****************************************************************
 * Spherical distance. For the special case of angular coordinates and 2D
 *****************************************************************/

float  sphdist(const Point &P1, const Point &P2)
{
    double tmp=0.;
 
    tmp = sin(P1.y()*D2R)*sin(P2.y()*D2R) + cos(P1.y()*D2R)*cos(P2.y()*D2R)*
      cos((P1.x()-P2.x())*D2R);
    return ((float) (acos(tmp)/D2R));

} 


/*****************************************************************/

float  squaredist(const Point &P1, const Point &P2)
{
    float Sum=0.0;
	int Dim=P1.dim();
	float temp;

    for (int i=0; i < Dim; i++)
	{
		temp=P1.axis(i)-P2.axis(i);
        Sum += temp*temp;
	}
    return Sum;
}


/*****************************************************************/

float  squaredist(const Point &P1, const Point &P2, float SquareDistMax)
{
    float Sum=0.0;
	int Dim=P1.dim();
	float temp;
	
    for (int i=0; i < Dim; i++)
	{
		temp=P1.axis(i)-P2.axis(i);
        Sum += temp*temp;
		if(Sum>SquareDistMax) break;
	}
    return Sum;
}

/*****************************************************************
 * Square of Spherical distance. For the special case of angular coordinates and 2D
 *****************************************************************/

float  squaresphdist(const Point &P1, const Point &P2)
{
    float tmp=0.;

    tmp = sin(P1.y()*D2R)*sin(P2.y()*D2R) + cos(P1.y()*D2R)*cos(P2.y()*D2R)*
      cos((P1.x()-P2.x())*D2R);
    return (acos(tmp)/D2R)*(acos(tmp)/D2R);

}




/*****************************************************************/

ArrayPoint::ArrayPoint(int Dimension, int N)
{  
   Np=0;
   Dim=0;
   TabPoint=NULL;
   alloc(Dimension, N);
   TCoord=1;BootCoord[0]=0;
   BootCoord[1]=0;BootCoord[2]=0;
}

/*****************************************************************/

void ArrayPoint::alloc(int Dimension, int N)
{    
    Dim = Dimension;
    Np = N;
    if (TabPoint != NULL) delete [] TabPoint;
    TabPoint = new Point[Np];
    Pmin.alloc(Dim);
    Pmax.alloc(Dim);
    for (int i=0; i < Np; i++) TabPoint[i].alloc(Dim);
}

/*****************************************************************/

void ArrayPoint::print(char *Mes)
{   
   if (Mes != NULL) cout << "PRINT " << Mes << endl;
   for (int i=0; i < Np; i++)  TabPoint[i].print();
}

/*****************************************************************/



/*****************************************************************

void ArrayPoint::random(ArrayPoint &TabData, type_random_cat RndCat)
{
   switch(RndCat)
   {
      case RND_CAT_USER: random(TabData);break;
      case RND_CAT_LAMBDA_CDM: random_lcdm(TabData);break;
      case RND_CAT_IRAS_NORTH: random_sub_iras(TabData, False);break;
      case RND_CAT_IRAS_SOUTH: random_sub_iras(TabData, True);break;
      case RND_CAT_IRAS: random_iras(TabData);break;
      default:
	  cerr << "Error: bad random catalogue type ... " << endl;
	  exit(-1);
	  break;
   }
}

*****************************************************************/
/*****************************************************************/

void ArrayPoint::random(ArrayPoint &TabData)
{
  ArrayPoint TmpData;
  TmpData.alloc(Dim,Np);

  int j;

  for (int d=0; d < TabData.dim(); d++)
    {
      switch (TabData.BootCoord[d])
	{
	case 0:
	  // Uniform between min,max
	  for (int i=0; i < Np; i++){
	    TabPoint[i].axis(d) = get_random(Pmin.axis(d),Pmax.axis(d));
	  }
	  break;
	case 1:
	  // Bootstrapped
	  for (int i=0; i < Np; i++)
	    {
	      j = (int) (get_random() * (TabData.np()-1));
	      TabPoint[i].axis(d) = TabData(j).axis(d);
	    }
	  break;
	case 2:
	  // Uniform on sphere -latitude
	  for (int i=0; i < Np; i++)
	    {
	      /*TabPoint[i].axis(d) = asin(cos(Pmin.axis(d)*D2R) +
					(cos(Pmax.axis(d)*D2R)-
					 cos(Pmin.axis(d)*D2R)*
					 get_random(0.0,1.0)))/D2R;*/
            TabPoint[i].axis(d)=asin(sin(Pmin.axis(d)*D2R) +
					(sin(Pmax.axis(d)*D2R)-
					 sin(Pmin.axis(d)*D2R))*
					 get_random(0.0,1.0))/D2R;
	    }
      break;
    case 3:
      // Uniform on sphere -distance
      for (int i=0; i < Np; i++)
	    {
            TabPoint[i].axis(d)=POW((POW(Pmax.axis(d),float(3.0))-POW(Pmin.axis(d),float(3.0)))*get_random(0.0,1.0)
                                +POW(Pmin.axis(d),float(3.0)),float(1/3.0));
	    }
	}
    }
}

/*****************************************************************

void ArrayPoint::random_uni(Point & PMi, Point & PMa)
{   
   Pmin = PMi;
   Pmax = PMa;
   for (int i=0; i < Np; i++)  
      TabPoint[i].random(Pmin,Pmax);
}

*****************************************************************/


void ArrayPoint::minmax(Point & PMi, Point & PMa)
{   
   int d,i;
   
   for (d=0; d < PMi.dim(); d++)
   {
      PMi.axis(d) = TabPoint[0].axis(d);
      PMa.axis(d) = TabPoint[0].axis(d);
   }
   for (i=0; i < Np; i++)
   for (d=0; d < PMi.dim(); d++)
   {
      if (PMi.axis(d) > TabPoint[i].axis(d)) PMi.axis(d) = TabPoint[i].axis(d);
      if (PMa.axis(d) < TabPoint[i].axis(d)) PMa.axis(d) = TabPoint[i].axis(d);
   }	
   
    PMi.print();
    PMa.print();
}

/*****************************************************************/

void ArrayPoint::read(char *FileName, Bool Verbose, Bool ReadSimu)
{
    int Naxis,N,i=0,Boot;
    
    FILE *FileDes = fopen(FileName,"r");
    if (FileDes == NULL) 
    {
       cerr << "Error: cannot open file " <<  FileName <<   endl;
       exit(-1);
    }
    if (fscanf(FileDes,"%d %d %d\n",&N,&Naxis,&TCoord) != 3)
    {
      cerr << "Error: cannot get the dimension and the size of the catalogue " <<  endl;
      exit(-1);
    }
    
    if ((Naxis < 1) || (Naxis > 3))
    {
       cerr << "Error: input data must be of dimension 1,2, or 3 ... " << endl;
       cerr << "Now it is " << Naxis << endl;
       exit(-1);
    }
    if ((TCoord < 1) || (TCoord > NBR_TYPE_COORD))
    {
       cerr << "Error: bad type of coordinate: " << TCoord << endl;
       cerr << "       0 < Coord_Type < " << NBR_TYPE_COORD+1 << endl;
       exit(-1);
    }
    
    if (Verbose == True) 
    {
       cout << "Array Dimension   = " <<  Naxis << endl;
       cout << "Number of data points = " <<  N << endl;
       cout << "Coordinate type = " << StringCoord(TCoord) << endl;
    }
    if ((Dim != Naxis) || (N != Np)) alloc(Naxis,N);
    if (ReadSimu == True)
    {
       if (Verbose == True) 
                 cout << "User supplied random catalogs " << endl;
    }
    else
    {
      for (i=0; i < Dim; i++)
      {
	float Min,Max;
	if (fscanf(FileDes,"%f\t%f\t%d\n",&Min,&Max,&Boot) != 3)
	  {
	    cerr << "Error: bad min and max range values ..." <<  endl;
	    exit(-1);
	  } 
	
	// Case of uniform on the sphere for the latitude/declination
	if (Boot == 2 && i==1) {
	  if (Max > 90.0) Max = 90.0;
	  if (Min < -90.0) Min = -90.0;
	}
	
	Pmin.axis(i) = Min;
	Pmax.axis(i) = Max;
	BootCoord[i] = Boot;
	if (Verbose == True) 
	  {
	    cout << "Axis " << i+1 << " Min = " << Min << " Max = " << Max;
	    
			switch(BootCoord[i])
			{
				case 0:
					cout << "  Uniform axis " << endl;
					break;
				case 1:
					cout << "  Bootstrapped axis " << endl;
					break;
				case 2:
					cout << "  Uniform on the sphere - latitude " << endl;
					break;
				case 3:
					cout << "  Uniform on the sphere - distance " << endl;
					break;
				default:
					cout << " Unknown generation method for the axis. Exitting!" << endl;
					exit(-1);
			}
	  }
      }
    }
    i=0;
    
    while (i < N)	
      {
	TabPoint[i].read(FileDes);
	// cout << i+1 << " D = ";
	// TabPoint[i].print();
	i++;
      }
    if (i != N)
      {
	cerr << "Error: bad number of read points from the file ... " << endl;
	cerr << "       Number of read points = " << i << endl;
	exit(-1);	
      }
    
    fclose(FileDes);

    
}
/*****************************************************************/
 
void ArrayPoint::toxyz()
{
  //conversion longitude, latitude, distance to cartesian coordinates
    ArrayPoint TmpData;
    TmpData.alloc(Dim,Np);

    for (int i=0; i < Np; i++)
      TmpData(i) = TabPoint[i];

    for (int i=0; i < Np; i++) {
      TabPoint[i].x() = TmpData(i).z()*cos(TmpData(i).y()*D2R)*
	cos(TmpData(i).x()*D2R);
      TabPoint[i].y() = TmpData(i).z()*cos(TmpData(i).y()*D2R)*
	sin(TmpData(i).x()*D2R);
      TabPoint[i].z() = TmpData(i).z()*sin(TmpData(i).y()*D2R);
	}
}
/*****************************************************************/
 
void ArrayPoint::write(char *FileName, Bool Verbose)
{
    int i=0;
    
    FILE *FileDes = fopen(FileName,"w");
    if (FileDes == NULL) 
    {
       cerr << "Error: cannot open file " <<  FileName <<   endl;
       exit(-1);
    }
    // Labatie line: if (fprintf(FileDes,"%d\t%d\t%d\n",Dim,Np,TCoord) < 0)
    if (fprintf(FileDes,"%d\t%d\t%d\n",Np,Dim,TCoord) < 0)
    {
       cerr << "Error: cannot write on file " <<  FileName << endl;
       exit(-1);
    }
    for (i=0; i < Dim; i++)
    {
       if (fprintf(FileDes,"%f\t%f\t%d\n",Pmin.axis(i),Pmax.axis(i),BootCoord[i]) < 0)
       {
          cerr << "Error: cannot write on file: " << FileName << endl;
           exit(-1);
       } 
    }
    i=0;
    while (i < Np)	
    {
      TabPoint[i].write(FileDes);
        i++;
    }
    fclose(FileDes);
}

/*****************************************************************/

void bootstrap(ArrayPoint & Data, ArrayPoint & BootStrapData)
{
   int i,j;
   int N = Data.np();
   for (i=0; i < N; i++)
   {
      j = (int) (get_random() * (N-1));
      if ((j<0) || (j >= N)) 
      {
         cerr << "Error: bad random value = " << j << endl;
	 exit(-1);
      }
      BootStrapData(i) = Data(j);
   }
   BootStrapData.Pmin = Data.Pmin;
   BootStrapData.Pmax = Data.Pmax;
}

/*****************************************************************/

void bootstrapxyz(ArrayPoint & Data, ArrayPoint & BootStrapData)
{
   int i,j;
   int Nb = BootStrapData.np();
   int Nd = Data.np();
   
   for (i=0; i < Nb; i++)
   {
      for (int d=0; d < Data.dim(); d++)
      {
         j = (int) (get_random() * (Nd-1));
	 if ((j<0) || (j >= Nd)) 
         {
            cerr << "Error: bad random value = " << j << endl;
	    exit(-1);
         }
	 BootStrapData(i).axis(d) = Data(j).axis(d);
      }
   }
   BootStrapData.Pmin = Data.Pmin;
   BootStrapData.Pmax = Data.Pmax;
}

/*****************************************************************/
void ArrayPoint::trace() const {

   std::cout << "ArrayPoint infos" << std::endl;
   std::cout << "  " << "- nb pts : " << Np << std::endl;
   std::cout << "  " << "- dim    : " << Dim << std::endl;
   std::cout << "  " << "- coord min : "; Pmin.trace();
   std::cout << "  " << "- coord max : "; Pmin.trace();
   


}
