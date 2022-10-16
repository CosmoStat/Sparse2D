/******************************************************************************
**
**    DESCRIPTION  cube segmentation  
**    ----------- 
**                 
**  void cube_segment (fltarray &cube, fltarray &Segment, int &NLabel, float Level)
**
**
**  cube : input cube we want to segment
**  Segment: output segmented cube
**  NLabel: Number of labels or regions
**  Level: segmentation level
**
*******************************************************************************/
#include "CUBE_Segment.h"
#include "IM_IO.h"
#include "IM3D_IO.h"


//----------------------------------------------------------
//	Pile
//----------------------------------------------------------
Pile::Pile(void)
{
	data = NULL;
	num = 0;
	bloc = 0;
}

//----------------------------------------------------------
//	Pile
//----------------------------------------------------------
Pile::Pile (int n)
{
	data = NULL;
	bloc = n;
	if (data != NULL) data = (int *)realloc (data, bloc * sizeof (int));
	else data = (int *) malloc (bloc * sizeof (int));

	num = 0;
}

//----------------------------------------------------------
//	Pile :: destructeur
//----------------------------------------------------------
Pile::~Pile (void)
{
	if (data) free (data);
}

//----------------------------------------------------------
//	Pile :: max retourne le plus grand ellement
//----------------------------------------------------------
int	Pile::max (void)
{
	int	max = -100000;
	for (int i=0 ; i<num ; i++)
		max = (max > data[i]) ? max : data[i];
	return max;
}

//----------------------------------------------------------
//	Pile :: allocation de la memoire de data
//----------------------------------------------------------
void	Pile::alloue (void)
{
	bloc += PILE_SIZE_BLOC;
	if (data != NULL) data = (int *)realloc (data, bloc * sizeof (int));
	else data = (int *) malloc (bloc * sizeof (int));
}

//----------------------------------------------------------
//	Pile :: ajoute une donnee
//----------------------------------------------------------
void	Pile::add_data (int d)
{
	if (num >= bloc)
		alloue ();
	data[num++] = d;
}

//----------------------------------------------------------
//	Pile :: trie les donnes
//----------------------------------------------------------

void	Pile::trie (int max)
{
	int	i, c;
	Pile	histo (max);
	for (i=0 ; i<num ; i++) histo.data[i] = 0;
	for (i=0 ; i<num ; i++) histo.data[data[i]]++;
	c = 1;
	for (i=1 ; i<num ; i++)
		if ((histo.data[i] != 0) && (i != 0))
			remplace (i, c++);
}

//----------------------------------------------------------
//	Pile :: remplace d par f
//----------------------------------------------------------
void	Pile::remplace (int d, int f)
{
	for (int i=0 ; i<num ; i++)
		if (data[i] == d)
			data[i] = f;
}



////----------------------------------------------------------
////----------------------------------------------------------
////----------------------------------------------------------
////----------------------------------------------------------
void cube_segment (fltarray & cube, intarray & Segment, 
	int &nb_seg, float S, Bool CleanBord, int FirstLabel)
{
    int Nx=cube.nx();
    int Ny=cube.ny();
    int Nz=cube.nz();
    
    int Nvoisin = 13;
    int TabNeighbour[Nvoisin];
    Pile P;
    int n = 1;
    P.add_data (0);

    if ((Segment.nx() != Nx) || (Segment.ny() != Ny)|| (Segment.nz() != Nz)) 
    	Segment.resize(Nx,Ny,Nz);
    Segment.init();

    // image segmentation
    for (int k=0; k<Nz; k++)
    for (int j=0; j<Ny; j++)
    for (int i=0; i<Nx; i++)
    {
       if (ABS(cube(i,j,k)) > S)  
       {
          Bool NewLabel= True;
          int UseLabel = 0;
          TabNeighbour[0] = Segment(i-1, j-1, k-1, I_ZERO);
	  TabNeighbour[1] = Segment(i  , j-1, k-1, I_ZERO);
	  TabNeighbour[2] = Segment(i+1, j-1, k-1, I_ZERO);
	  TabNeighbour[3] = Segment(i-1, j, k-1, I_ZERO);
	  TabNeighbour[4] = Segment(i,   j, k-1, I_ZERO);
	  TabNeighbour[5] = Segment(i+1, j, k-1, I_ZERO);
	  TabNeighbour[6] = Segment(i-1, j+1, k-1, I_ZERO);
	  TabNeighbour[7] = Segment(i,   j+1, k-1, I_ZERO);
	  TabNeighbour[8] = Segment(i+1, j+1, k-1, I_ZERO);
	  
	  TabNeighbour[9]  = Segment(i-1, j-1, k, I_ZERO);
 	  TabNeighbour[10] = Segment(i, j-1,   k, I_ZERO);
	  TabNeighbour[11] = Segment(i+1, j-1, k, I_ZERO);
	  TabNeighbour[12] = Segment(i-1, j, k, I_ZERO);
 
          // Find the label to apply to the pixel
          for (int l=0; l < Nvoisin; l++) 
          { 
              if (TabNeighbour[l] != 0) 
              {
                 int L = TabNeighbour[l];
		 int ColNeighbour = P.data[L];
                 NewLabel = False;
                 if (L > n)
                 {
                      cout << "Error: Segmentation problem ... " << endl;
                      cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                      exit(-1);
                 }
                 // TabNeighbour[l] = P.data[L];
                 if (UseLabel == 0) UseLabel = ColNeighbour;
                 else if (UseLabel > ColNeighbour) UseLabel = ColNeighbour;
              }
          }
          // New label
          if (NewLabel == True)
	  {
             Segment(i,j,k) = n;
	     P.add_data (n++);
 	  }
	  else 
          {
 	      Segment(i,j,k) = UseLabel;
              // neighbourhood label must have the same label
              // and therefore may need to be changed
              for (int l=0; l < Nvoisin; l++)
              {
                  if (TabNeighbour[l] != 0)
		  {
		     int L = TabNeighbour[l];
		     int ColNeighbour = P.data[L]; 
		     if (ColNeighbour != UseLabel)
                     {
                        if (L > n)
                        {
                           cout << "Error: Segmentation problem .... " << endl;
                           cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                           exit(-1);
                        }
                        P.remplace (P.data[L], UseLabel);
		     }
                  }
              }
           }
       }
   }

   if (CleanBord == True)
   {
      // cout << "CLEAN BORD" << endl;
      // plan (*,*,0) et plan (*,*,Nz-1)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Ny ; j++)
      {
         if (Segment(i,j,0) != 0)  P.remplace((int) Segment(i,j,0) ,0);
         if (Segment(i,j,Nz-1) != 0)  P.remplace((int) Segment(i,j,Nz-1) ,0);
      } 
      // plan (*,0,*) et plan (*,Ny-1,*)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Nz ; j++)
      {
         if (Segment(i,0,j) != 0)  P.remplace((int) Segment(i,0,j) ,0);
         if (Segment(i,Ny-1,j) != 0) P.remplace((int) Segment(i,Ny-1,j),0);
      }
      // plan (0,*,*) et plan (Nx-1,*,*)
      for (int i=0 ; i< Ny ; i++)
      for (int j=0 ; j< Nz ; j++)
      {
         if (Segment(0,i,j) != 0)  P.remplace((int) Segment(0,i,j), 0);
         if (Segment(Nx-1,i,j) != 0) P.remplace((int) Segment(Nx-1,i,j), 0);
      }
    } // END CLEAN bord

    P.trie (n);
    int FL = FirstLabel-1;
    for (int i=0 ; i < Nx; i++)
    for (int j=0 ; j < Ny; j++)
    for (int k=0 ; k < Nz; k++)
    {
        int L = Segment(i,j,k);
        if (L > n)
        {
            cout << "Error: Segmentation problem ..... " << endl;
            cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
            exit(-1);
        }
	Segment(i,j,k) = P.data[L];
	if (Segment(i,j,k) != 0) Segment(i,j,k) += FL;
    }
    nb_seg = P.max();
}



////----------------------------------------------------------
////----------------------------------------------------------
////----------------------------------------------------------
////----------------------------------------------------------
void cube_segment (intarray & cube, intarray & Segment, 
	int &nb_seg, float S, Bool CleanBord, int FirstLabel)
{
    int Nx=cube.nx();
    int Ny=cube.ny();
    int Nz=cube.nz();
    

    int TabNeighbour[10];
    float old_colour;
    Pile P;
    int n = 1;
    P.add_data (0);

    if ((Segment.nx() != Nx) || (Segment.ny() != Ny)|| (Segment.nz() != Nz)) 
    	Segment.resize(Nx,Ny,Nz);
    Segment.init();

    // First pixel cube(0, 0, 0)
    if (ABS(cube(0,0,0)) > S){
	Segment(0,0,0) = 1;
	P.add_data (n++);
    }
 	
    // X axe cube(1:Nx , 0, 0)
    for (int i=1 ; i< Nx ; i++) 
    	if (abs(cube(i,0,0)) > S){ 
        	if (Segment(i-1,0,0) != 0) 
			Segment(i,0,0) = Segment(i-1,0,0);
        	else{
           		Segment(i,0,0) = n;
           		P.add_data (n++);
		} 
    	}

    // Y axe cube(0, 1:Ny, 0)
    for (int i=1 ; i< Ny ; i++) 
    	if (abs(cube(0,i,0)) > S){
        	if (ABS(Segment(0,i-1,0)) != 0 ) 
			Segment(0,i,0) = Segment(0,i-1,0);
        	else{
	    		Segment(0,i,0) = n;
            		P.add_data (n++);
		}
    	}
	
    // Z axe cube(0, 0, 1:Nz)
    for (int i=1 ; i< Nz ; i++) 
    	if (abs(cube(0,0,i)) > S){
        	if (ABS(Segment(0,0,i-1)) != 0 ) 
			Segment(0,0,i) = Segment(0,0,i-1);
        	else{
	    		Segment(0,0,i) = n;
            		P.add_data (n++);
		}
    	}
    
    // image segmentation
    for (int i=1; i<Nx; i++)
    for (int j=1; j<Ny; j++)
    for (int k=1; k<Nz; k++)
    {
       if (abs(cube(i,j,k)) > S)  
       {
          Bool NewLabel= True;
          int UseLabel = 0;

          TabNeighbour[0] = Segment(i-1,j-1,k-1);
	  TabNeighbour[1] = Segment(i,j-1,k-1);
	  TabNeighbour[2] = (i+1<Nx) ? Segment(i+1,j-1,k-1): 0;
	  TabNeighbour[3] = Segment(i-1,j,k-1);
	  TabNeighbour[4] = Segment(i,j,k-1);
	  TabNeighbour[5] = (j+1<Ny) ? Segment(i-1,j+1,k-1): 0;
	  TabNeighbour[6] = Segment(i-1,j-1,k);
	  TabNeighbour[7] = Segment(i,j-1,k);
	  TabNeighbour[8] = Segment(i-1,j,k);
	  TabNeighbour[9] = (k+1<Nz) ? Segment(i-1,j-1,k+1): 0;


          // Find the label to apply to the pixel
          for (int l=0; l < 10; l++) 
          { 
              if (TabNeighbour[l] != 0) 
              {
                 int L = TabNeighbour[l];
                 NewLabel = False;
                 if (L > n)
                 {
                      cout << "Error: Segmentation problem ... " << endl;
                      cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                      exit(-1);
                 }
                 TabNeighbour[l] = P.data[L];
                 if (UseLabel == 0) UseLabel = TabNeighbour[l];
                 else if (UseLabel > TabNeighbour[l]) UseLabel = TabNeighbour[l];
              }
          }

          // New label
          if (NewLabel == True)
	  {
             Segment(i,j,k) = n;
	     P.add_data (n++);
	  }
	  else 
          {
	      Segment(i,j,k) = UseLabel;
              // neighbourhood label must have the same label
              // and therefore may need to be changed
              for (int l=0; l < 10; l++)
              {
                  if ((TabNeighbour[l] != 0) && (TabNeighbour[l] != UseLabel))
                  {
                     int L = TabNeighbour[l];
                     if (L > n)
                     {
                        cout << "Error: Segmentation problem .... " << endl;
                        cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                        exit(-1);
                     }
                     P.remplace (P.data[L], UseLabel);
                  }
              }
           }
       }
   }

   if (CleanBord == True)
   {
      // plan (*,*,0)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Ny ; j++)
       if (Segment(i,j,0) != 0) 
      {
          old_colour = Segment(i,j,0);        
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
       
      }
       // plan (*,*,Nz-1)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Ny ; j++)
       if (Segment(i,j,Nz-1) != 0) 
      {
          old_colour = Segment(i,j,Nz-1);         
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
       
      }
 
       // plan (*,0,*)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Nz ; j++)
       if (Segment(i,0,j) != 0) 
      {
          old_colour = Segment(i,0,j);        
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
       
      }
       
       // plan (*,Ny-1,*)
      for (int i=0 ; i< Nx ; i++)
      for (int j=0 ; j< Nz ; j++)
       if (Segment(i,Ny-1,j) != 0) 
      {
          old_colour = Segment(i,Ny-1,j);       
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
       
      }
      
       // plan (0,*,*)
      for (int i=0 ; i< Ny ; i++)
      for (int j=0 ; j< Nz ; j++)
       if (Segment(0,i,j) != 0) 
      {
          old_colour = Segment(0,i,j);        
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
      }
      
      // plan (Nx-1,*,*)
      for (int i=0 ; i< Ny ; i++)
      for (int j=0 ; j< Nz ; j++)
       if (Segment(Nx-1,i,j) != 0) 
      {
          old_colour = Segment(Nx-1,i,j);        
          P.remplace((int)old_colour,0);
          for(int a=0; a< Nx; a++)
	  for(int b=0; b< Ny; b++)
	  for(int c=0; c< Nz; c++)
          	if (Segment(a,b,c) == old_colour) Segment(a,b,c) = 0;
      }
      
    } // END CLEAN bord

    P.trie (n);
    int FL = FirstLabel-1;
    for (int i=0 ; i < Nx; i++)
    for (int j=0 ; j < Ny; j++)
    for (int k=0 ; k < Nz; k++)
    {
        int L = Segment(i,j,k);
        if (L > n)
        {
            cout << "Error: Segmentation problem ..... " << endl;
            cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
            exit(-1);
        }
	Segment(i,j,k) = P.data[Segment(i,j,k)];
	if (Segment(i,j,k) != 0) Segment(i,j,k) += FL;
    }
    nb_seg = P.max();
 }
