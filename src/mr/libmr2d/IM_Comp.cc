/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  IM_Comp.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_Comp.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_Comp.h"

/*********************************************************************/
char * StringComp (type_comp type)
{
    switch (type)
    {
       case COMP_INT_LIFTING: 
              return ("Compression by lifting scheme");break;
       case COMP_MRMEDIAN: 
            return ("Compression by the pyramidal median transform"); 
            break;
       case COMP_O_MRMEDIAN: 
            return ("Compression by the pyramidal median transform using integer"); 
            break;       
       case COMP_MORPHO: 
            return ("Compression by Mathematical morphology");
            break;
       case COMP_O_MORPHO: 
            return ("Compression by Mathematical morphology using integer");
            break;
        case COMP_HAAR: 
            return ("Compression by Haar transform");
            break;
       case COMP_MALLAT: 
            return ("Compression by Mallat-Daubechies transform");
            break;
       case COMP_FEAUVEAU: 
            return ("Compression by the Feauveau transform");
            break;
       case COMP_MIN_MAX: 
            return ("Compression by the min-max transform (or G transform)");
            break;  
       case COMP_O_MIN_MAX: 
            return ("Compression by the min-max transform (or G transform) using integer");
            break;              
       case COMP_WT_PMT:
            return ("Compression by mixed WT and PMT method");
            break;
       default:
            cerr << "Unknown compression method ... " << endl;
            exit(-1);
    }
}


/*********************************************************************/

float im_comp_entrop (Ifloat &Pict)
{
    int i,j,Nbr_Val,ind;
    int *Tab_Histo=NULL;
    float Prob_Ev;
    float Min,Max,Entr;
    int Nl = Pict.nl();
    int Nc = Pict.nc();
    int Size = Nl*Nc;

    Min = min(Pict);
    Max = max (Pict);

    /* Calcul de l'entropie */
    Nbr_Val = (int) (Max - Min) + 1;

    Tab_Histo = new int [Nbr_Val];    
    for (i = 0; i < Nbr_Val; i++) Tab_Histo [i] = 0;
    
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        ind = iround(Pict(i,j) - Min);
        Tab_Histo [ind] ++;
    }
    Entr = 0.;

    for (i = 0; i < Nbr_Val; i++)
    {
        if (Tab_Histo [i] > 0)
        {
            Prob_Ev = (float) Tab_Histo [i] / (float) Size;
            Entr += - Prob_Ev * log (Prob_Ev) / log(2.);
        }
    }

    delete [] Tab_Histo;
    return (Entr);
}

/****************************************************************************/

int im_nbr_level (Ifloat &Pict)
{
    int Min,Max,Nbr_Val;
    int Puis = 0, Max_Val = 1;

    Min = iround (min(Pict));
    Max = iround (max (Pict));

    Nbr_Val = Max - Min + 1;
    while (Max_Val < Nbr_Val)
    {
        Max_Val *= 2;
        Puis ++;
    }
    return (Puis);
}

/****************************************************************************/

float im_comp_entrop (int *Pict, int Size)
{
    int i,Nbr_Val,ind;
    int *Tab_Histo=NULL;
    float Prob_Ev;
    int Min,Max;
    float Entr;

    Min = Max = Pict [0];
    for (i = 1; i < Size; i++)
    {
        if (Pict [i] > Max) Max = Pict [i];
        if (Pict [i] < Min) Min = Pict [i];
    }

    /* Calcul de l'entropie */
    Nbr_Val = Max - Min + 1;
    Tab_Histo = new int [Nbr_Val];
    
    for (i = 0; i < Nbr_Val; i++) Tab_Histo [i] = 0;
    
    for (i = 0; i < Size; i++)
    {
        ind = Pict [i] - Min;
        Tab_Histo [ind] ++;
    }
    Entr = 0.;
    for (i = 0; i < Nbr_Val; i++)
    {
        if (Tab_Histo [i] > 0)
        {
            Prob_Ev = (float) Tab_Histo [i] / (float) Size;
            Entr += - Prob_Ev * log (Prob_Ev) / log(2.);
        }
    }
    if (Tab_Histo != NULL) delete [] Tab_Histo;
    return (Entr);
}

/****************************************************************************/
