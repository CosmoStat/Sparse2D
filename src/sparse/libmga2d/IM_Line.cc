/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/08/02 
**    
**    File:  IM_Line.cc
**
*******************************************************************************
**
**    DESCRIPTION  Line detection in an image
**    -----------  
**                 
**
*******************************************************************************/

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_Line.h"
#include "NR.h"

/******************************************************************************/

void symgetline (int x2, int y2, int Nl2, int Nc2, intarray & TX, intarray & TY, int & Np)
{
   int x1 = 0;
   int y1 = 0;
   int N = MAX(TX.nx(),TY.nx());
   intarray TX0(N);
   intarray TY0(N);
   int N2 = MAX(Nl2,Nc2);
   getline (x1,y1,x2,y2,TX0, TY0,Np);
   for (int i = 0; i < Np; i++)
   {
       if ((i!= 0) && (TX0(i) == 0) && (TY0(i) == 0))
 	    TX0(i) = TY0(i) = 1;
       TY(N2-i) =  -TY0(i) + Nl2;
       TX(N2-i) =  -TX0(i) + Nc2;
       TY(N2+i) =  TY0(i) + Nl2;
       TX(N2+i) =  TX0(i) + Nc2;
   }
   Np = 2*Np-1;   
}

/******************************************************************************/
 
void getline (int x1,int y1,int x2,int y2, intarray & TX, intarray & TY, int & N)
/* 
   draw a line in the image Image between points of coordinates
   x1,y1 and x2,y2. 
*/
{
    int j,pasx,x,y,dx,dy,cumul;
    float pasy;
    N = 0;
    if (x1 < x2) pasx = 1;
    else pasx = -1;
    if (y1 < y2) pasy = 1;
    else pasy = -1;
    dx = abs (x2 - x1);
    dy = abs (y2 - y1);
    x= x1;
    y = y1;
    TX(N) = x;
    TY(N) = y1;
    N++;
     
    if (dx > dy)
    {
        cumul = dx / 2;
        for (j = 0; j < dx; j++)
        {
            x += pasx;
            cumul += dy;
            if (cumul > dx)
            {
                cumul -= dx;
                y += (int) pasy;
            }
 	    TX(N) = x; TY(N) = y; N++;
        }
    }
    else
    {
        cumul = dy / 2;
        for (j = 0; j < dy; j++)
        {
            y += (int) pasy;
            cumul += dx;
            if (cumul > dy)
            {
                cumul -= dy;
                x += pasx;
            }
 	    TX(N) = x; TY(N) = y; N++;
        } 
    }
}   
    
/******************************************************************************/

static void getline (int x1,int y1,int x2,int y2, Atom & AT)
/* 
   draw a line in the image Image between points of coordinates
   x1,y1 and x2,y2. 
*/
{
    int j,pasx,x,y,dx,dy,cumul;
    float pasy;
    intarray TX(MAXLINE);
    intarray TY(MAXLINE);
    int N = 0;
    if (x1 < x2) pasx = 1;
    else pasx = -1;
    if (y1 < y2) pasy = 1;
    else pasy = -1;
    dx = abs (x2 - x1);
    dy = abs (y2 - y1);
    x= x1;
    y = y1;
    TX(N) = x;
    TY(N) = y1;
    N++;
     
    if (dx > dy)
    {
        cumul = dx / 2;
        for (j = 0; j < dx; j++)
        {
            x += pasx;
            cumul += dy;
            if (cumul > dx)
            {
                cumul -= dx;
                y += (int) pasy;
            }
 	    TX(N) = x; TY(N) = y; N++;
        }
    }
    else
    {
        cumul = dy / 2;
        for (j = 0; j < dy; j++)
        {
            y += (int) pasy;
            cumul += dx;
            if (cumul > dy)
            {
                cumul -= dy;
                x += pasx;
            }
 	    TX(N) = x; TY(N) = y; N++;
        } 
    }
    AT.alloc(N);
    for (int i=0; i < N; i++)
    {
       AT.setpos(i, TX(i), TY(i));
    }
}       
    
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
    
ListAtom::ListAtom()  
{
    Verbose = False;
    Norm = True;
    Step = 1;
    int i = 0,j;
    int X, Y;
    Na = 45;
    TabAtom = new Atom [Na];
       
    // point atom
    TabAtom[0].alloc(1);         
    TabAtom[0].setpos(0, 0, 0);  

    // Line in block of size 3: 4 directions
    int SizeBlock = 3;
    int S2 = SizeBlock / 2;
    i++; X = 1; Y = 0;
    getline (X, Y, -X, -Y, TabAtom[i]);
    i++; X = 1; Y = 1;
    getline (X, Y, -X, -Y, TabAtom[i]);
    i++; X = 0; Y = 1;
    getline (X, Y, -X, -Y, TabAtom[i]);
    i++; X = -1; Y = 1;
    getline (X, Y, -X, -Y, TabAtom[i]);
    if (Verbose == True) cout << " END BLOCK 3: Ni = " << i << endl;
       
    // Line in block of size 5:  8 directions
    SizeBlock = 5;
    S2 = SizeBlock / 2;
    for (j=0; j <= S2; j++)
    {
       i++; X = S2; Y = j;
       getline (X, Y, -X, -Y, TabAtom[i]);
    }
    for (j=S2-1; j >= -S2; j--)
    {
        i++; X = j; Y = S2;
        getline (X, Y, -X, -Y, TabAtom[i]);
    }
    for (j=S2-1; j > 0; j--)
    {
        i++; X = -S2; Y = j;
        getline (X, Y, -X, -Y, TabAtom[i]);
    }
    if (Verbose == True)    cout << " END BLOCK 5: Ni = " << i << endl;

    // Line in block of size 17:  32 directions
    SizeBlock = 17;
    S2 = SizeBlock / 2;
    for (j=0; j <= S2; j++)
    {
        i++; X = S2; Y = j;
        getline (X, Y, -X, -Y, TabAtom[i]);
    }
    for (j=S2-1; j >= -S2; j--)
    {
        i++; X = j; Y = S2;
        getline (X, Y, -X, -Y, TabAtom[i]);
    }
    for (j=S2-1; j > 0; j--)
    {
        i++; X = -S2; Y = j;
        getline (X, Y, -X, -Y, TabAtom[i]);
    }
    if (Verbose == True) cout << " END BLOCK 5: Ni = " << i << endl;       
       
}
/******************************************************************************/
    
float ListAtom::atom_pos_correl(Ifloat &Ima, int AtNumber, int i, int j)
{
    int Nl = Ima.nl();
    int Nc = Ima.nc();
    float C = 0.;
    int Np = TabAtom[AtNumber].np();
    for (int p=0; p < Np; p++)           
    {
       int Posi = i+TabAtom[AtNumber].ypos(p)*Step;
       int Posj = j+TabAtom[AtNumber].xpos(p)*Step;
       if ((Posi >= 0) && (Posi < Nl) && (Posj >= 0) && (Posj < Nc))
 	      C +=  Ima(Posi,Posj);           
       else   C +=  Ima(Posi,Posj,I_MIRROR);
    }
    if (Norm == True) C /=  sqrt((float) Np);
    return C;   
}

/******************************************************************************/

void ListAtom::get_best_correl(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma)
{
   int i,j,a,Ind;
   float Val, Max, AMax, AVal;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   if ((Conv.nl() != Nl) || (Conv.nc() != Nc)) Conv.resize(Nl,Nc);
   if ((AtomIma.nl() != Nl) || (AtomIma.nc() != Nc)) AtomIma.resize(Nl,Nc);  
   
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
      Max = atom_pos_correl(Ima,0,i,j);
      AMax = ABS(Max);
      Ind = 0;
      for (a=1; a < Na; a++)
      {
  	  Val = atom_pos_correl(Ima,a,i,j);
	  AVal = ABS(Val);
	  if (AVal > AMax)
	  {
	     Max = Val;
	     AMax = AVal;
	     Ind = a;
	  }
      }
      Conv(i,j) = Max;
      AtomIma(i,j) = Ind;
   }
}

/******************************************************************************/

void ListAtom::add_atom(Ifloat &Ima, int i, int j, float ConvVal, int AtNumber)
{
    int Nl = Ima.nl(), Nc = Ima.nc();
    int Np = TabAtom[AtNumber].np();
    float Coef = (Norm == True) ? ConvVal / sqrt((float) Np): ConvVal / (float) Np;
    
    if (Verbose == True)
      cout << "(" << i << "," << j << ") = " <<  ConvVal << "  " << Coef << endl;
    for (int p=0; p < TabAtom[AtNumber].np(); p++) 
    {
       int Posi = i+TabAtom[AtNumber].ypos(p)*Step;
       int Posj = j+TabAtom[AtNumber].xpos(p)*Step;
       if ((Posi >= 0) && (Posi < Nl) && (Posj >= 0) && (Posj < Nc))
 	      Ima(Posi,Posj) += Coef;           
    }
}

/******************************************************************************/

void ListAtom::add_nfirst_atoms(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                        unsigned long * & TabInd, int NFirst)
{
   int i,Ind,IM,JM;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int  N = Nl*Nc;
   for (i = 0; i < NFirst; i++)
   {
      Ind = N - i;
      get_index_pos(Ind,TabInd,Nl,IM,JM);
      add_atom(Ima, IM,JM, Conv(IM,JM), AtomIma(IM,JM));
   }
}			

/******************************************************************************/

void ListAtom::add_nsigma_atoms(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                        unsigned long * & TabInd, float NSigma, 
			float SigmaNoise)
{
   int IM,JM;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int  N = Nl*Nc;
   int Ind = N;
   int Cpt = 0;
   float Level = NSigma*SigmaNoise;
   Bool Stop = False;
   
   do
   {
      get_index_pos(Ind,TabInd,Nl,IM,JM);
      if (ABS(Conv(IM,JM)) > Level)
          add_atom(Ima, IM,JM, Conv(IM,JM), AtomIma(IM,JM));
      else Stop = True;
      Ind --;
      Cpt ++;
   } while ((Ind > 0) && (Stop == False));
   cout << "Level = " << Level << " Number of significant coefficients = " << Cpt << " ==> " << ((float) Cpt / (float) N) * 100. << endl;	
}	
	
/******************************************************************************/

void ListAtom::order_coef(Ifloat &Ima, Ifloat &Conv, Iint &AtomIma,
                          unsigned long * & TabInd)
{
   int i;
   int IM, JM;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   float *Buff = Conv.buffer()-1;
   unsigned long N = Nl*Nc;
   int Ind = N;
   if (TabInd != NULL) delete [] TabInd;
   TabInd = new  unsigned long [N+1];
   Buff = new  float [N+1];
   for (i=0; i < (int) N; i++) Buff[i+1] = ABS(Conv(i));
   indexx(N, Buff, TabInd);
   get_index_pos(Ind, TabInd, Nl, IM, JM);
   do
   {
       // subtract the atom to the image
       add_atom(Ima, IM, JM,  - Conv(IM,JM), AtomIma(IM,JM)); 
       Ind --;
       if (Ind >= 1)
       {
          get_index_pos(Ind, TabInd, Nl, IM, JM);
 	  Conv(IM,JM) = atom_pos_correl(Ima,AtomIma(IM,JM),IM,JM);
	  Buff[TabInd[Ind]] = ABS(Conv(IM,JM));
       }
   } while (Ind > 0);
   indexx(N, Buff, TabInd);
   delete [] Buff;
}

/******************************************************************************/
