
/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: 10/30/02
**
**    Date:  02/11/27 
**    
**    File:  IM_JoinUpMax.h
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
#include "IM_Obj.h"
#include "MR_Obj.h"

#define VAL_INTER 1
#define FIRST_LEVEL 20

typedef enum { E_ZERO                        =0,
                   E_ONE_TREATED                 =1,
	               E_ONE_NOT_TREATED             =2,
                   E_TWO_NOT_TREATED             =3,
                   E_TWO_TREATED                 =4,
                   E_TWO_TREATED_AND_NOT_TREATED =5,
                   E_MORE                        =6} 
	te_TypePix;
	
	
class JoinUpMax {

private:


	
	Ifloat* po_Data;
 	Iint*   po_TreatedPixel;
    int     i_ValInter;
    int     i_FirstLevel;
	int     i_CurrentLevel;
    int     Nl, Nc;
    Bool    InterNotTreated;

public:

	JoinUpMax (Ifloat& pro_DataIn, Bool pe_InterNotTreated=True, 
               int pi_ValInter=VAL_INTER, int pi_FirstLevel=FIRST_LEVEL);
    void Compute ();

private:
    void Test_Pixel (int i, int j);
    te_TypePix How_Much_Max_Around (int i, int j, int& pi_NbTreatedPix, 
                                    Iint& po_CoordTreatedPix, 
                                    int& pi_NbNotTreatedPix, 
                                    Iint& po_CoordNotTreatedPix);
    void IsPixelMaxAndTreated (int i, int j, int& pi_NbTreatedPix, 
                                    Iint& po_CoordTreatedPix, 
                                    int& pi_NbNotTreatedPix, 
                                    Iint& po_CoordNotTreatedPix);

	//JoinUpMax (JoinUpMax& po_JoinUpMax) {};
};
           

void mr2d_JoinUpMax (MultiResol& po_Mr2dDataIn, MultiResol& po_Mr2dDataOut,
                     Bool pe_InterNotTreated=True, int pi_ValInter=VAL_INTER, 
                     int pi_FirstLevel=FIRST_LEVEL);
