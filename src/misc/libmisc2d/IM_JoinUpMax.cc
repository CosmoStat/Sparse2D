/*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: 10/30/02
**
**    Date:  02/10/30 
**    
**    File:  IM_JoinUpMax.cc
**
*******************************************************************************
**
**    DESCRIPTION routines for noise treatment 
**    -----------  
**                 
*******************************************************************************/
#define NB_DIRECTION 8
#define NB_COORD 2
#define X_COORD 0
#define Y_COORD 1


#include "IM_JoinUpMax.h"
extern Bool Verbose;



JoinUpMax::JoinUpMax (Ifloat& pro_DataIn, Bool pe_InterNotTreated, 
                      int pi_ValInter, int pi_FirstLevel){

	// init attributs
	i_ValInter = pi_ValInter;
	i_FirstLevel = pi_FirstLevel;
	i_CurrentLevel = pi_FirstLevel;
    po_Data = &pro_DataIn;
    Nl = po_Data->nl(); 
    Nc = po_Data->nc();
    po_TreatedPixel = new Iint (Nl, Nc, "");
    po_TreatedPixel->init (0);
    InterNotTreated = pe_InterNotTreated;
}

void JoinUpMax::Compute () {

    for (int i=0; i<Nl; i++) {
		for (int j=0; j<Nc; j++) {

			// if pixel level is max (>-NOT_USED) and pixel not treated
            if (    (*po_Data)(i,j) > - NOT_USED
                 && (*po_TreatedPixel)(i,j) == 0) {

				this->Test_Pixel (i, j);	
			}
		}
	}
}
 

void JoinUpMax::Test_Pixel (int i, int j) {
	
	// local var
	int ai_NbTreatedPix=0, ai_NbNotTreatedPix=0;
    Iint ao_CoordTreatedPix(NB_DIRECTION, NB_COORD, "");
    Iint ao_CoordNotTreatedPix(NB_DIRECTION, NB_COORD, "");

	// pixel around point (i,j)
	te_TypePix ae_NbMaxPixelAround = this->How_Much_Max_Around (i, j, 
                                     ai_NbTreatedPix, ao_CoordTreatedPix,
                                     ai_NbNotTreatedPix, ao_CoordNotTreatedPix);

	// pixel (i,j) is treated
	(*po_TreatedPixel)(i,j) = 1;

	// for all max pixel around (i,j)
	switch (ae_NbMaxPixelAround) {

	case E_ZERO :                           // a max line of one point (i,j)
        i_CurrentLevel++;
        (*po_Data)(i,j) = i_CurrentLevel;
        break;

	case E_ONE_TREATED :                    // end a max line (of one point if 
                                            // the neighbour is an inter
        if (    !InterNotTreated 
             && (*po_Data)(ao_CoordTreatedPix(0,X_COORD),
                           ao_CoordTreatedPix(0,Y_COORD)) == 1) {
			i_CurrentLevel++; // if neighbour is inter => begin a line
		}
		(*po_Data)(i,j) = i_CurrentLevel;
        break;

	case E_ONE_NOT_TREATED :                // begin a max line
  		i_CurrentLevel++;
        (*po_Data)(i,j) = i_CurrentLevel;             
		this->Test_Pixel (ao_CoordNotTreatedPix(0,X_COORD), 
                          ao_CoordNotTreatedPix(0,Y_COORD));
        break;

	case E_TWO_NOT_TREATED :                // begin a max line with two branch
		                                    // if pix are deconnected
		// si les deux pixels sont connexes on a une inetrsection, sinon on
        // a une chaine vu par un maillon interieur
        // 2 pixels ne sont pas connexes si la dif d'ordonnees est = a 2 
        // en x ou en y
		if (InterNotTreated) {
  			i_CurrentLevel++;
			(*po_Data)(i,j) = i_CurrentLevel;
			for (int i=0;i<ai_NbNotTreatedPix;i++)
				this->Test_Pixel (ao_CoordNotTreatedPix(i,X_COORD), 
                                  ao_CoordNotTreatedPix(i,Y_COORD));
		} else {

        	if (   ABS(ao_CoordNotTreatedPix(0,X_COORD)-
                       ao_CoordNotTreatedPix(1,X_COORD)) > 1	
                || ABS(ao_CoordNotTreatedPix(0,Y_COORD)-
                       ao_CoordNotTreatedPix(1,Y_COORD)) > 1) {
	
        		(*po_Data)(i,j) = i_CurrentLevel;
				this->Test_Pixel (ao_CoordNotTreatedPix(0,X_COORD), 
                                  ao_CoordNotTreatedPix(0,Y_COORD));
				this->Test_Pixel (ao_CoordNotTreatedPix(1,X_COORD), 
                                  ao_CoordNotTreatedPix(1,Y_COORD));
        	} else {
				// intersection
        		(*po_Data)(i,j) = i_ValInter;
			}
		}
        break;

	case E_TWO_TREATED_AND_NOT_TREATED :    // in a max line
		// si les deux pixels sont connexes on a une inetrsection, sinon on
        // a une chaine 
        // 2 pixels ne sont pas connexes si la dif d'ordonnees est = a 2 
        // en x ou en y
		if (InterNotTreated) {
			(*po_Data)(i,j) = i_CurrentLevel;
			for (int i=0;i<ai_NbNotTreatedPix;i++)
				this->Test_Pixel (ao_CoordNotTreatedPix(i,X_COORD), 
                                  ao_CoordNotTreatedPix(i,Y_COORD));
		} else {

        	if (    ABS(ao_CoordNotTreatedPix(0,0)-
                        ao_CoordTreatedPix(0,X_COORD)) > 1
         	    ||  ABS(ao_CoordNotTreatedPix(0,1)-
                        ao_CoordTreatedPix(0,Y_COORD)) > 1) {

				// if treated point is an inter => begin a new line
        		if ( (*po_Data)(ao_CoordTreatedPix(0,X_COORD),
                                ao_CoordTreatedPix(0,Y_COORD)) == 1) {
					i_CurrentLevel++; // if neighbour is inter => begin a line
				}
        		(*po_Data)(i,j) = i_CurrentLevel;
				this->Test_Pixel (ao_CoordNotTreatedPix(0,X_COORD), 
                                  ao_CoordNotTreatedPix(0,Y_COORD));
        	} else {
				// intersection
        		(*po_Data)(i,j) = i_ValInter;
			}
		}
        break;

	case E_MORE : 
    case E_TWO_TREATED :                          // intersection
		if (InterNotTreated) {
			if (ai_NbTreatedPix == 0)  i_CurrentLevel++; // begin max line
			if (ai_NbNotTreatedPix == 0) {         // end of max line
        		(*po_Data)(i,j) = i_CurrentLevel;
            } else {                           // at least one  branch is not 
                                               // traited
				(*po_Data)(i,j) = i_CurrentLevel;
				for (int i=0;i<ai_NbNotTreatedPix;i++)
					this->Test_Pixel (ao_CoordNotTreatedPix(i,X_COORD), 
                                      ao_CoordNotTreatedPix(i,Y_COORD));
			}
		} else {
        	(*po_Data)(i,j) = i_ValInter;
		}
		break;

	default : 
		cout << "Pb in case of marquer_pixel" << endl;
		cout << (int)ae_NbMaxPixelAround << endl;
		exit(-1); break;
 
	}
}


te_TypePix JoinUpMax::How_Much_Max_Around (int i, int j, 
                         int& pi_NbTreatedPix, Iint& po_CoordTreatedPix, 
                         int& pi_NbNotTreatedPix, Iint& po_CoordNotTreatedPix) {


	// test SO pixel
    if (i!=0 && j!=0) IsPixelMaxAndTreated (i-1, j-1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test O pixel
    if (i!=0) IsPixelMaxAndTreated (i-1, j, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test NO pixel
    if (i!=0 && j<Nc-1) IsPixelMaxAndTreated (i-1, j+1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test N pixel 
    if (j<Nc-1) IsPixelMaxAndTreated (i, j+1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test NE pixel 
    if (i<Nl-1 && j<Nc-1) IsPixelMaxAndTreated (i+1, j+1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test E pixel 
    if (i<Nl-1) IsPixelMaxAndTreated (i+1, j, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test SE pixel 
    if (i<Nl-1 && j!=0) IsPixelMaxAndTreated (i+1, j-1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// test S pixel 
    if (j!=0) IsPixelMaxAndTreated (i, j-1, pi_NbTreatedPix, 
                 po_CoordTreatedPix, pi_NbNotTreatedPix, po_CoordNotTreatedPix);

	// return code
	if (pi_NbTreatedPix==0 && pi_NbNotTreatedPix==0) return E_ZERO;
    if (pi_NbTreatedPix==1 && pi_NbNotTreatedPix==0) return E_ONE_TREATED;
    if (pi_NbTreatedPix==0 && pi_NbNotTreatedPix==1) return E_ONE_NOT_TREATED;
    if (pi_NbTreatedPix==2 && pi_NbNotTreatedPix==0) return E_TWO_TREATED;
    if (pi_NbTreatedPix==0 && pi_NbNotTreatedPix==2) return E_TWO_NOT_TREATED;
    if (pi_NbTreatedPix==1 && pi_NbNotTreatedPix==1) 
 		return E_TWO_TREATED_AND_NOT_TREATED;
	return E_MORE;

}


void JoinUpMax::IsPixelMaxAndTreated (int i, int j, 
                         int& pi_NbTreatedPix, Iint& po_CoordTreatedPix, 
                         int& pi_NbNotTreatedPix, Iint& po_CoordNotTreatedPix) {

	// if pixel is  max
	if ( (*po_Data)(i,j) > -NOT_USED) {

		if ((*po_TreatedPixel)(i,j) > 0) {
			po_CoordTreatedPix(pi_NbTreatedPix,0)=i;
			po_CoordTreatedPix(pi_NbTreatedPix,1)=j;
			pi_NbTreatedPix++;
		} else {
			po_CoordNotTreatedPix(pi_NbNotTreatedPix,0)=i;
			po_CoordNotTreatedPix(pi_NbNotTreatedPix,1)=j;
			pi_NbNotTreatedPix++;
		}
	}
}


void mr2d_JoinUpMax (MultiResol& po_Mr2dDataIn, MultiResol& po_Mr2dDataOut,
                     Bool pe_InterNotTreated, int pi_ValInter, int pi_FirstLevel) {

	// inter reconstruction file
	char atc_FileName[20];
	int s=0;

    switch (SetTransform(po_Mr2dDataIn.Type_Transform)) {

	case TRANSF_PAVE: 
	case TRANSF_SEMIPYR:
	case TRANSF_PYR: 

		// for all scale (nbr_scale-1 values)
  		for (s=0; s<po_Mr2dDataIn.nbr_scale()-1; s++) {
			Ifloat ao_MaxCurrentScale = po_Mr2dDataIn.extract_band(s);

			if (Verbose) {
				sprintf (atc_FileName, "mr2d_MaxCurrentScale_%d", s);
				io_write_ima_float (atc_FileName, ao_MaxCurrentScale);
			}

			JoinUpMax ao_JoinUpMax (ao_MaxCurrentScale, pe_InterNotTreated, 
                                    pi_ValInter, pi_FirstLevel);
    		ao_JoinUpMax.Compute();
			po_Mr2dDataOut.insert_band (ao_MaxCurrentScale, s);

			if (Verbose) {
				sprintf (atc_FileName, "mr2d_JoinUpMax_%d", s);
				io_write_ima_float (atc_FileName, ao_MaxCurrentScale);
			}
		}
		break;

   
	case TRANSF_DIADIC_MALLAT:

		// for all modulus (nbr_scale-1 values)
  		for (s=0; s<po_Mr2dDataIn.nbr_scale()-1; s++) { 
			Ifloat ao_MaxModulus = po_Mr2dDataIn.extract_band(2*s);

			if (Verbose) {
				sprintf (atc_FileName, "mr2d_MaxModulus_%d", s);
				io_write_ima_float (atc_FileName, ao_MaxModulus);
			}

			JoinUpMax ao_JoinUpMax (ao_MaxModulus, pe_InterNotTreated, 
                                    pi_ValInter, pi_FirstLevel);
    		ao_JoinUpMax.Compute();
			po_Mr2dDataOut.insert_band (ao_MaxModulus, 2*s);

			if (Verbose) {
				sprintf (atc_FileName, "mr2d_JoinUpMax_%d", s);
				io_write_ima_float (atc_FileName, ao_MaxModulus);
			}
		}
		break;

	case TRANSF_FEAUVEAU:
	case TRANSF_MALLAT: 
	default: 
		cout << "Not yet implemented" << endl; 
		exit(-1); break;
    }
}
