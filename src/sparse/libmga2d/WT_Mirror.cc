 /*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/08/98
**
**    File:  SB_Filter.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Sub-Band decomposition
**    -----------
**
******************************************************************************/


#include "IM_Obj.h"
#include "WT_Mirror.h"

/****************************************************************************/

float MIRROR_2D_WT::get_val(Ifloat *TabTrans, int b, int i, int j)
{
   // int Nlb = size_band_nl(b);
   // int Ncb = size_band_nc(b);
   int IndBand = TabIndBand(b);
   int Depi = TabDepi(b);
   int Depj = TabDepj(b);
   return (TabTrans[IndBand])(Depi+i,Depj+j);
}

/****************************************************************************/

void MIRROR_2D_WT::put_val(Ifloat *TabTrans, int b, int i, int j, float Val)
{
   // int Nlb = size_band_nl(b);
   // int Ncb = size_band_nc(b);
   int IndBand = TabIndBand(b);
   int Depi = TabDepi(b);
   int Depj = TabDepj(b);
   (TabTrans[IndBand])(Depi+i,Depj+j) = Val;
}

/****************************************************************************/

float & MIRROR_2D_WT::operator() (int b, int i, int j)
{
//   int Nlb = size_band_nl(b);
//   int Ncb = size_band_nc(b);
   int IndBand = TabIndBand(b);
   int Depi = TabDepi(b);
   int Depj = TabDepj(b);
//   cout << "Depi = " << Depi << " i = " << i <<  " Nlb = " << Nlb << endl;
//   cout << "Depj = " << Depj << " j = " << j <<  " Ncb = " << Ncb << endl;
//   cout << "TabTransWT[IndBand]).nl() " << (TabTransWT[IndBand]).nl() << endl;
//   cout << "TabTransWT[IndBand]).nc() " << (TabTransWT[IndBand]).nc() << endl;
   return (TabTransWT[IndBand])(Depi+i,Depj+j);
}

/****************************************************************************/

void MIRROR_2D_WT::get_band(Ifloat &Band, int b)
{
    get_band(TabTransWT, Band, b);
}

/****************************************************************************/

void MIRROR_2D_WT::put_band(Ifloat &Band, int b)
{
    put_band(TabTransWT, Band, b);
}
/****************************************************************************/

void MIRROR_2D_WT::get_band(Ifloat *TabTrans, Ifloat &Band, int b)
{
   int Nlb = size_band_nl(b);
   int Ncb = size_band_nc(b);
   int IndBand = TabIndBand(b);
   int Depi = TabDepi(b);
   int Depj = TabDepj(b);
   Band.resize(Nlb,Ncb);
//   cout << "IndBand = " << IndBand << endl;
//   cout << "TabTrans[IndBand]).nl() " << (TabTrans[IndBand]).nl() << endl;
//   cout << "TabTrans[IndBand]).nc() " << (TabTrans[IndBand]).nc() << endl;

   for (int i=0; i < Nlb; i++)
   for (int j=0; j < Ncb; j++) Band(i,j) = (TabTrans[IndBand])(Depi+i,Depj+j);
}

/****************************************************************************/

void MIRROR_2D_WT::put_band(Ifloat *TabTrans, Ifloat &Band, int b)
{
   int Nlb = size_band_nl(b);
   int Ncb = size_band_nc(b);
   int IndBand = TabIndBand(b);
   int Depi = TabDepi(b);
   int Depj = TabDepj(b);
   for (int i=0; i < Nlb; i++)
   for (int j=0; j < Ncb; j++)  (TabTrans[IndBand])(Depi+i,Depj+j) = Band(i,j);
}

/****************************************************************************/

void MIRROR_2D_WT::free()
{
    if (Nbr_Plan != 0) delete [] TabTransWT;
}

/****************************************************************************/

void MIRROR_2D_WT::set_tabdec(int NumUndec, int Nbr_Plan)
{
    int i;
    int NU = (NumUndec < 0) ? Nbr_Plan: NumUndec;
    TabDec = new Bool [Nbr_Plan];
    for (i=0; i < MIN(NU,Nbr_Plan); i++) TabDec[i] = False;
    for (i=MIN(NU,Nbr_Plan); i < Nbr_Plan; i++) TabDec[i] = True;
}

/****************************************************************************/

void MIRROR_2D_WT::alloc (int Nl, int Nc, int Nbr_Plan,int NumUndec)
{
    set_tabdec(NumUndec,Nbr_Plan);
    make_alloc (Nl,Nc,Nbr_Plan);
}
/****************************************************************************/

void MIRROR_2D_WT::alloc (int Nl, int Nc, int Nbr_Plan, Bool *TabDecBand)
{
   TabDec = new Bool [Nbr_Plan];
   for (int s=0; s < Nbr_Plan; s++) TabDec[s] = TabDecBand[s];
   make_alloc (Nl, Nc, Nbr_Plan);
}
/****************************************************************************/

void MIRROR_2D_WT::make_alloc (int Nl, int Nc, int NbrPlan)
{
    char ch[80];
    int Nl_s = Nl;
    int Nc_s = Nc;
    int NbrBand_per_Resol = 3;
    int Nlb, Ncb, s, s1, s2;
    int NStep = NbrPlan-2;

    NlIma = Nl;
    NcIma = Nc;
    Nbr_Plan = NbrPlan;
    NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;
    TabTransWT = new Ifloat [NbrBand];
    // cout << "Number of WT bands = " <<  NbrBand << " " << Nbr_Plan << endl;

    // Calculate the total number of bands
     if (Nbr_Plan > 2)
     {
        int Nb0 = 0;
        int Nb1 = 0;
        int Nb2 = 0;
        for (s = NStep; s > 0; s--)
	{
 	   Nb0 += s;
	   Nb1 += s;
	   Nb2 += 3;
        }
	NbrTotBand = NbrBand + Nb0 + Nb1 + Nb2;
     }
     else NbrTotBand = NbrBand;
     TabNl.alloc(NbrTotBand);
     TabNc.alloc(NbrTotBand);
     TabDepi.alloc(NbrTotBand);
     TabDepj.alloc(NbrTotBand);
     TabIndBand.alloc(NbrTotBand);
     TabResolBand.alloc(NbrTotBand);

     //int NbrBand2Resol2 = NbrBand_per_Resol*(Nbr_Plan-2)+1;
     //int FirstIndResol2 = NbrTotBand - NbrBand2Resol2;
     // cout << "Number of bands = " <<  NbrTotBand << " " << FirstIndResol2 << endl;
    int IndBand=0;
    for (s = 0; s < NbrBand-1; s+=3)
    {
        Bool Dec = TabDec[s/3];
	Nlb = Nl_s;
	Ncb = Nc_s;

	// *********************************
	// Horizontal band
	// *********************************
        if (Dec == True)
        {
           Nlb = (Nl_s+1)/2;
           Ncb =  Nc_s/2;
        }
	sprintf (ch, "band_%d", s+1);
        TabTransWT[s].alloc (Nlb, Ncb, ch);
	if ((s == 0) && (Nbr_Plan > 2))
	{
           int NStep = Nbr_Plan-2;
	   int Depj = 0;
           for (s1 = NStep; s1 > 0; s1--)
	   {
	      int Ncb2 = (Ncb+1) / 2;
 	      TabNl(IndBand) = Nlb;
	      TabNc(IndBand) = Ncb2;
	      TabDepi(IndBand) = 0;
	      TabDepj(IndBand) = Depj;
	      TabIndBand(IndBand) = 0;
	      TabResolBand(IndBand) = s1-NStep;
	      Ncb = Ncb2;
	      Depj += Ncb2;
	      IndBand++;
           }
	   TabNl(IndBand) = Nlb;
	   TabNc(IndBand) = Ncb;
	   TabDepi(IndBand) = 0;
	   TabDepj(IndBand) = Depj;
	   TabIndBand(IndBand) = 0;
	   TabResolBand(IndBand) = -NStep;
	   IndBand++;
	}
	else
	{
	   TabNl(IndBand) = Nlb;
	   TabNc(IndBand) = Ncb;
	   TabDepi(IndBand) = 0;
	   TabDepj(IndBand) = 0;
	   TabIndBand(IndBand) = s;
	   TabResolBand(IndBand) = s/3;
	   IndBand++;
 	}

        // *********************************
	// Vertical band
	// *********************************

	Nlb = Nl_s;
	Ncb = Nc_s;
	if (Dec == True)
        {
            Nlb = Nl_s/2;
            Ncb = (Nc_s+1)/2;
	}
        sprintf (ch, "band_%d", s+2);
        TabTransWT[s+1].alloc (Nlb, Ncb, ch);
	if ((s == 0) && (Nbr_Plan > 2))
	{
           int NStep = Nbr_Plan-2;
	   int Depi = 0;
           for (s1 = NStep; s1 > 0; s1--)
	   {
	      int Nlb2 = (Nlb+1) / 2;
 	      TabNl(IndBand) = Nlb2;
	      TabNc(IndBand) = Ncb;
	      TabDepi(IndBand) = Depi;
	      TabDepj(IndBand) = 0;
	      TabIndBand(IndBand) = 1;
 	      TabResolBand(IndBand) = s1-NStep;
	      Nlb = Nlb2;
	      Depi += Nlb2;
	      IndBand++;
           }
	   TabNl(IndBand) = Nlb;
	   TabNc(IndBand) = Ncb;
	   TabDepi(IndBand) = Depi;
	   TabDepj(IndBand) = 0;
	   TabIndBand(IndBand) = 1;
	   TabResolBand(IndBand) = -NStep;
	   IndBand++;
	}
	else
	{
	   TabNl(IndBand) = Nlb;
	   TabNc(IndBand) = Ncb;
	   TabDepi(IndBand) = 0;
	   TabDepj(IndBand) = 0;
	   TabIndBand(IndBand) = s+1;
	   TabResolBand(IndBand) = s/3;
	   IndBand++;
 	}



        // *********************************
	// Diagonal band
	// *********************************
	Nlb = Nl_s;
	Ncb = Nc_s;
	if (Dec == True)
        {
            Nlb = Nl_s/2;
            Ncb = Nc_s/2;
	}
        sprintf (ch, "band_%d", s+3);
        TabTransWT[s+2].alloc(Nlb, Ncb,ch);
	if ((s == 0) && (Nbr_Plan > 2))
	{
           int NStep = Nbr_Plan-2;
	   int Depi = 0;
	   int Depj = 0;
           for (s1 = NStep; s1 > 0; s1--)
	   {
	      int Nlb2 = (Nlb+1) / 2;
	      int Ncb2 = (Ncb+1) / 2;
 	      TabNl(IndBand) = Nlb2;  // Smooth
	      TabNc(IndBand) = Ncb2;
	      TabDepi(IndBand) = Depi;
	      TabDepj(IndBand) = Depj;
	      TabIndBand(IndBand) = 2;
 	      TabResolBand(IndBand) = s1-NStep-1;
	      IndBand++;
	      if (s1 < 2)
  	      {
  	         TabNl(IndBand) = Nlb2;    // Horizontal
	         TabNc(IndBand) = Ncb/2;
	         TabDepi(IndBand) = Depi;
	         TabDepj(IndBand) = Depj+Ncb2;
	         TabIndBand(IndBand) = 2;
		 TabResolBand(IndBand) = s1-NStep-1;
 	         IndBand++;
		 TabNl(IndBand) = Nlb/2;  // Vertical
	         TabNc(IndBand) = Ncb2;
	         TabDepi(IndBand) = Depi+Nlb2;
	         TabDepj(IndBand) = Depj;
	         TabIndBand(IndBand) = 2;
		 TabResolBand(IndBand) = s1-NStep-1;
	         IndBand++;
		 TabNl(IndBand) = Nlb/2;  // Diagonal
	         TabNc(IndBand) = Ncb/2;
	         TabDepi(IndBand) = Depi+Nlb2;
	         TabDepj(IndBand) = Depj+Ncb2;
	         TabIndBand(IndBand) = 2;
		 TabResolBand(IndBand) = s1-NStep-1;
	         IndBand++;
 	      }
              else
	      {
	         // Horizontal
		 int H_Nlb = Nlb2;
		 int H_Ncb = Ncb/2;
 		 int H_Depi = Depi;
 		 int H_Depj = Depj+Ncb2;
 	         for (s2 = s1-1; s2 > 0; s2--)
		 {
 		     int H_Ncb2 = H_Ncb/2;
  	             TabNl(IndBand) = H_Nlb;
	             TabNc(IndBand) = H_Ncb2;
	             TabDepi(IndBand) = H_Depi;
	             TabDepj(IndBand) = H_Depj;
	             TabIndBand(IndBand) = 2;
		     TabResolBand(IndBand) = -NStep+s2-1;
 	             H_Depj += H_Ncb2;
	             IndBand++;
		     H_Ncb /= 2;
                 }
	         TabNl(IndBand) = H_Nlb;
	         TabNc(IndBand) = H_Ncb;
	         TabDepi(IndBand) = H_Depi;
	         TabDepj(IndBand) = H_Depj;
	         TabIndBand(IndBand) = 2;
		 TabResolBand(IndBand) = -NStep;
	         IndBand++;

 		 // Vertical
		 int V_Nlb = Nlb/2;
		 int V_Ncb = Ncb2;
 		 int V_Depi = Depi+Nlb2;
 		 int V_Depj = Depj;
	         for (s2 = s1-1; s2 > 0; s2--)
		 {
		     int V_Nlb2 = V_Nlb/2;
  	             TabNl(IndBand) = V_Nlb2;
	             TabNc(IndBand) = V_Ncb;
	             TabDepi(IndBand) = V_Depi;
	             TabDepj(IndBand) = V_Depj;
	             TabIndBand(IndBand) = 2;
 		     TabResolBand(IndBand) = -NStep+s2-1;
	             V_Depi += V_Nlb2;
	             IndBand++;
		     V_Nlb /= 2;
		 }
		 TabNl(IndBand) = V_Nlb;
	         TabNc(IndBand) = V_Ncb;
	         TabDepi(IndBand) = V_Depi;
	         TabDepj(IndBand) = V_Depj;
	         TabIndBand(IndBand) = 2;
		 TabResolBand(IndBand) = -NStep;
	         IndBand++;
	      }
  	      Depi += Nlb2;
	      Depj += Ncb2;
	      Nlb = Nlb/2;
	      Ncb = Ncb/2;
            }
	}
	else
	{
	   TabNl(IndBand) = Nlb;
	   TabNc(IndBand) = Ncb;
	   TabDepi(IndBand) = 0;
	   TabDepj(IndBand) = 0;
	   TabIndBand(IndBand) = s+2;
           TabResolBand(IndBand) = s/3;
	   IndBand++;
 	}

	if (Dec == True)
        {
            Nl_s = (Nl_s+1)/2;
            Nc_s = (Nc_s+1)/2;
        }
     }

     // Last scale
     s = NbrBand-1;
     IndBand = NbrTotBand-1;
     Nlb = Nl_s;
     Ncb = Nc_s;
     sprintf (ch, "band_%d", s+1);
     TabTransWT[s].alloc(Nlb, Ncb, ch);
     TabNl(IndBand) = Nlb;
     TabNc(IndBand) = Ncb;
     TabDepi(IndBand) = 0;
     TabDepj(IndBand) = 0;
     TabIndBand(IndBand) = s;
     TabResolBand(IndBand) = s/3;

//      for (s=0; s < NbrTotBand; s++)
//      {
//         cout << "Band " << s+1 << " Nl = " << TabNl(s) << " Nc = " << TabNc(s);
// 	cout << " Depi = " << TabDepi(s) << " Depj = " << TabDepj(s);
// 	cout << " IndBand = " << TabIndBand(s) << endl;
//      }
}

/****************************************************************************/

void MIRROR_2D_WT::get_tab_mirrorband(Ifloat * & TabMirrorBand)
{
   get_tab_mirrorband(TabTransWT, TabMirrorBand);
}

/****************************************************************************/

void MIRROR_2D_WT::get_tab_mirrorband(Ifloat *TabBand,
                                       Ifloat * & TabMirrorBand)
{
   int s;
   TabMirrorBand = new Ifloat [NbrTotBand];
   for (s = 0; s < NbrTotBand; s++)
        get_band(TabBand,TabMirrorBand[s], s);
}

/****************************************************************************/

void MIRROR_2D_WT::put_tab_mirrorband(Ifloat * & TabMirrorBand)
{
   put_tab_mirrorband(TabTransWT, TabMirrorBand);
}

/****************************************************************************/

void MIRROR_2D_WT::put_tab_mirrorband(Ifloat *TabBand,
                                       Ifloat * & TabMirrorBand)
{
   int s;
   for (s = 0; s < NbrTotBand; s++) put_band(TabBand,TabMirrorBand[s], s);
}

/****************************************************************************/

void MIRROR_2D_WT::transform(Ifloat &Imag)
{
   int i,j;
   // Step = 1;
   // int NbrBand = 3*(Nbr_Plan-1)+1;
   int Nl = Imag.nl();
   int Nc = Imag.nc();

   if ((Nl != NlIma) || (Nc != NcIma))
   {
      cout << "Error: the class in not initializaed for this image size ... " << endl;
      cout << "       Class: NlIma = " << NlIma << " NcIma = " << NcIma << endl;
      cout << "       Image: Nl = " << Nl << " Nc = " << Nc << endl;
   }

   int Nls = Nl;
   int Ncs = Nc;
   Ifloat Ima_Aux(Nl,Nc,"aux");
   Ifloat ImaH,ImaV,ImaD,ImaS;
   int NbStep = Nbr_Plan-1;
   // TabTrans = (Ptr_HD->alloc)(

   (*Ptr_HD).transform(Imag,TabTransWT,Nbr_Plan,TabDec);
   Ptr_SB1D->DistPix = (TabDec[0] == True) ? 1: 2;
   // INFO_X(TabTrans[0], "BAND HORIZONTAL");
   // INFO_X(TabTrans[1], "BAND VERTICAL");
   // INFO_X(TabTrans[2], "BAND DIAGONAL");
   if (NbStep > 1)
   {
      //io_write_ima_float("t1.fits", TabTrans[0]);
      (*Ptr_LC).transform(TabTransWT[0], NbStep, 0);
      //io_write_ima_float("t2.fits", TabTrans[0]);
      (*Ptr_LC).transform(TabTransWT[1], 0, NbStep);
   }
   Ima_Aux = TabTransWT[2];
   Nls = Ima_Aux.nl();
   Ncs = Ima_Aux.nc();
   int Depi=0;
   int Depj=0;

   for (int s = 0; s < NbStep-1; s++)
   {
       int Nl2 = (Nls+1)/2;
       int Nc2 = (Ncs+1)/2;
       // cout << "Scale " << s+1;
       // cout << "  Imag: Nl = " << Ima_Aux.nl() << " Nc = " << Ima_Aux.nc() <<  " Depi = " <<  Depi <<  " Depj = " <<  Depj <<  endl;
       ImaH.alloc(Nl2, Ncs/2);
       ImaV.resize(Nls/2, Nc2);
       ImaD.resize(Nls/2, Ncs/2);
       ImaS.resize(Nl2, Nc2);
       // cout << "   STEP OWT " << endl;
       (*Ptr_SBT).transform2d(Ima_Aux, False, &ImaH, &ImaV, &ImaD, &ImaS);
       if (NbStep-s-1 > 1)
       {
           // cout << "      --> Horizontal "  << " NStep = " << NbStep-s-1 << ImaH.nl() << " " << ImaH.nc() << endl;
          (*Ptr_LC).transform(ImaH, NbStep-s-1, 0);
           // cout << "      --> Vertical " << ImaV.nl() << " " << ImaV.nc() << endl;
 	  (*Ptr_LC).transform(ImaV, 0, NbStep-s-1);
       }
       // cout << "   COPY SMOOTH Nl = " << Nl2 << " Nc = " << Nc2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nl2; i++)
       for (j = 0; j < Nc2; j++)
                 (TabTransWT[2])(i+Depi,j+Depj) = ImaS(i,j);

       // cout << "   COPY Horizontal " << Nl2 << " " << Ncs/2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nl2; i++)
       for (j = 0; j < Ncs/2; j++)
                      (TabTransWT[2])(i+Depi,j+Depj+Nc2) = ImaH(i,j);

       // cout << "   COPY Vertical " << Nls/2 << " " << Nc2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nls/2; i++)
       for (j = 0; j < Nc2; j++)
                      (TabTransWT[2])(i+Depi+Nl2,j+Depj) = ImaV(i,j);

       // cout << "   COPY Diag " << Nls/2 << " " << Ncs/2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nls/2; i++)
       for (j = 0; j < Ncs/2; j++)
                      (TabTransWT[2])(i+Depi+Nl2,j+Depj+Nc2) = ImaD(i,j);
       Nls = Nls/2;
       Ncs = Ncs/2;
       Depi += Nl2;
       Depj += Nc2;
       Ima_Aux.resize(Nls,Ncs);
       Ima_Aux = ImaD;
       // INFO_X(ImaD, "New Diag");
    }
    Ptr_SB1D->DistPix = 1;

}

/****************************************************************************/

void MIRROR_2D_WT::recons (Ifloat &Imag)
{
   int i,j;
   // s,Step = 1;
   // int NbrBand = 3*(Nbr_Plan-1)+1;
   int Nl = NlIma;
   int Nc = NcIma;
   Imag.resize(Nl,Nc);
   int Nls = Nl;
   int Ncs = Nc;
   Ifloat Ima_Aux(Nl,Nc,"aux");
   Ifloat ImaH,ImaV,ImaD,ImaS;
   int NbStep = Nbr_Plan-1;

   Nls = TabTransWT[2].nl();
   Ncs = TabTransWT[2].nc();
   int Depi=0;
   int Depj=0;
   Ptr_SB1D->DistPix = (TabDec[0] == True) ? 1: 2;

   // for (int s = 0; s < NbStep-1; s++)
   for (int s = NbStep-2; s >= 0; s--)
   {
      int Nl2,Nc2;
      Nls = TabTransWT[2].nl();
      Ncs = TabTransWT[2].nc();
      Depi = Depj = 0;
      for (i = 0; i < s; i++)
      {
         Nl2 = (Nls+1)/2;
         Nc2 = (Ncs+1)/2;
         Nls = Nls/2;
         Ncs = Ncs/2;
         Depi += Nl2;
         Depj += Nc2;
      }
      Nl2 = (Nls+1)/2;
      Nc2 = (Ncs+1)/2;
      Ima_Aux.resize(Nls,Ncs);
      // cout << "Scale " << s+1;
      // cout << "  Imag: Nl = " << Ima_Aux.nl() << " Nc = " << Ima_Aux.nc() <<  " Depi = " <<  Depi <<  " Depj = " <<  Depj <<  endl;
      ImaH.resize(Nl2, Ncs/2);
      ImaV.resize(Nls/2, Nc2);
      ImaD.resize(Nls/2, Ncs/2);
      ImaS.resize(Nl2, Nc2);


       // cout << "   COPY SMOOTH Size = " << Nl2 << " " << Nc2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nl2; i++)
       for (j = 0; j < Nc2; j++)
                     ImaS(i,j) =(TabTransWT[2])(i+Depi,j+Depj);

       // cout << "   COPY Horizontal Size = " << Nl2 << " " << Ncs/2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nl2; i++)
       for (j = 0; j < Ncs/2; j++)
                    ImaH(i,j) = (TabTransWT[2])(i+Depi,j+Depj+Nc2);

       // cout << "   COPY Vertical Size = " << Nls/2 << " " << Nc2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       for (i = 0; i < Nls/2; i++)
       for (j = 0; j < Nc2; j++)
                     ImaV(i,j) = (TabTransWT[2])(i+Depi+Nl2,j+Depj);

       // cout << "   COPY Diag Size = " << Nls/2 << " " << Ncs/2 << " Depi = " << Depi << " Depj = " << Depj << endl;
       if (s == NbStep-2)
       {
          for (i = 0; i < Nls/2; i++)
          for (j = 0; j < Ncs/2; j++)
                     ImaD(i,j) = (TabTransWT[2])(i+Depi+Nl2,j+Depj+Nc2);
       }
       // INFO_X(ImaD, "New Diag");

       if (NbStep-s-1 > 1)
       {
           // cout << "      --> Horizontal " << " NStep = " << NbStep-s-1  << ImaH.nl() << " " << ImaH.nc() << endl;
          (*Ptr_LC).recons(ImaH, NbStep-s-1, 0);
           // cout << "      --> Vertical " << ImaV.nl() << " " << ImaV.nc() << endl;
 	  (*Ptr_LC).recons(ImaV, 0, NbStep-s-1);
       }

       // cout << "   STEP OWT " << endl;
       (*Ptr_SBT).recons2d(Ima_Aux, False, &ImaH, &ImaV, &ImaD, &ImaS);

       ImaD.resize(Ima_Aux.nl(),Ima_Aux.nc());
       if (s != 0) ImaD = Ima_Aux;
       else TabTransWT[2] = Ima_Aux;
    }

   // Ptr_SB1D->DistPix = 2;
   if (NbStep > 1)
   {
      // io_write_ima_float("t1.fits", TabTransWT[0]);
      (*Ptr_LC).recons(TabTransWT[0], NbStep, 0);
      // io_write_ima_float("t2.fits", TabTransWT[0]);
      (*Ptr_LC).recons(TabTransWT[1], 0, NbStep);
   }
   Ptr_SB1D->DistPix = 1;
   // INFO_X(TabTransWT[0], "BAND HORIZONTAL");
   // INFO_X(TabTransWT[1], "BAND VERTICAL");
   // INFO_X(TabTransWT[2], "BAND DIAGONAL");
   (*Ptr_HD).recons(TabTransWT, Imag ,Nbr_Plan,TabDec);
}

/****************************************************************************/
