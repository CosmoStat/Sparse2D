/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  10/07/04 
**    
**    File:  
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:   
**    ----------- 
**
******************************************************************************/

#ifndef _WPACKETS_H_
#define _WPACKETS_H_

#include "SB_Filter.h"
#include "map"
#include "vector"

struct Pt {int x; int y;};
struct InfoBlock {int IndImage; int Num; Pt BeginCoord; Pt Size;};
struct CmpInfoBlock :
   public binary_function <const InfoBlock, const InfoBlock, bool> {
   bool operator () (const InfoBlock& Key1, const InfoBlock& Key2) const {
      return (Key1.Num < Key2.Num ? true : false);
   }
};
struct InfoTrans {InfoBlock Parent; vector<InfoBlock> Child;};
typedef  map <InfoBlock,InfoTrans,CmpInfoBlock> MapInfo;
typedef  map <int,Ifloat*> MapTransf;
 

class WPTransf_2D {

   friend class WPACKETS_2D;

public:
	 WPTransf_2D();
	 ~WPTransf_2D() {}
         int nbr_band() const;
	 int size_band_nl(int b) const;
	 int size_band_nc(int b) const;
         float & operator() (int NumBand, int i, int j);
 	 void get_band (int NumBand, Ifloat & Band);          
	 void put_band (int NumBand, Ifloat & Band);
         void info_band (int NumBand);
         void info(std::string name=std::string(""), Bool write=False);

private:

        // attributes
        //        
        MapTransf  _Transf;
	
        intarray _TabNl;
	intarray _TabNc;
	intarray _TabDepi;
	intarray _TabDepj;
	intarray _TabIndIma;
	int _WPNbrBand;
       
        //internal functions
        //
        WPTransf_2D (const WPTransf_2D& dummy);
};


class WPACKETS_2D {

public:

         // ctor, dtor
         //
         WPACKETS_2D (SubBand1D &SB1D);
	 WPACKETS_2D (SubBand1D &SB1D_Line, SubBand1D &SB1D_Col);
         ~WPACKETS_2D ();
         
         void alloc (int Nl, int Nc, int NbrPlan, int NbrUndecimatedScale);
         void free();
         
         void transform (Ifloat& Imag, WPTransf_2D*& TabTrans);
	 void recons (WPTransf_2D*& TabTrans, Ifloat& Imag);
	 void threshold (WPTransf_2D* TabTrans, float LambdaSigma);
	 void threshold (WPTransf_2D* TabTrans, float SigmaNoise, 
                         float LambdaSigma);
         
         void trace() const;
         void write_transf(string& name) const;   
         void trace_all_blocks() const;
	 void set_verbose (const Bool Flag);
	 void set_WP (const Bool Flag);   
        
         
private:

        // prov
        WPTransf_2D* _pWP2D;

        // attributes
        //
        SubBand1D *_Ptr_SB1D_Line; // Pointer to the 1D subband decomposition along lines
	SubBand1D *_Ptr_SB1D_Col; // Pointer to the 1D subband decomposition along columns
        int _NlIma;
        int _NcIma;
        
        int _NbPlan;
        int _NbUndec;
        int _CurrentScale;
        int _SizeTransf;
        vector < MapInfo > _Infos;
	Bool _WP;
	intarray _Num;
	Bool _Verbose;
        float _DataSigmaNoise;
        
        Bool _SetInfoBand;
        
        //internal functions
        //
        void wpRecursifAnalysis (Ifloat& Plane, const Pt& LeftDownCoord, 
                                 int IndImage, int Num);
        void wpRecursifDecimated (Ifloat& Plane, const Pt& LeftDownCoord, 
                                 int IndImage, int Num);
        void wpRecursifUndecimated (Ifloat& Plane, int IndImage, int Num);
        
        void wp_Recons (WPTransf_2D*& TabTrans, int _CurrentScale, 
                        const InfoTrans& InfoBlock, int Step);  
        
        void save_res_dec (Ifloat const * Res, const Pt& LeftDownCoord, 
                           int IndImage);
        void save_res_undec (Ifloat const * Res, int IndImage);
                                                                        
        Bool is_scale_undecimated (int scale);
        InfoTrans get_info_block (int Scale, InfoBlock Parent) const ;
        void set_info (InfoBlock parentBlock, 
                       vector<InfoBlock>& childBlock);
        unsigned int get_block_number (int Scale) const;
 
        void trace_info_block(const InfoTrans& Inf, int Scale) const;
     	 
        void ortho_threshold (WPTransf_2D* TabTrans, float LambdaSigma);
        void wp_threshold (WPTransf_2D* TabTrans, float LambdaSigma);        
        
        void init_struct ();
        void wpRecursifAlloc(Ifloat& Plane, int IndImage, int Num);

        void set_info_band();
        void init_infos ();
};

/****************************************************************************/
 
#endif
