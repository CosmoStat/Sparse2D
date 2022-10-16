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

#ifndef _MRSAT_H_
#define _MRSAT_H_

#include "map"
#include "vector"
#include "MR_Obj.h"
#include "ostream"


struct InfoCoef {
   float Coef;
   float Sigma;
   float SNR;
};

struct WaveletPos {
   int xPos;
   int yPos;
   int scale;
};



class TreeCoef {
public:
   TreeCoef (WaveletPos Pos, vector<float>& w,  float SigmaImag, 
             Bool Overlap, Bool CteFalseDetect, ofstream* os, 
	     Bool ComputeSNR = True);
   void get_all_coef (std::vector <float>& Coef);
   void get_all_normalized_coef (std::vector <float>& NormCoef);
   void get_last_level (std::vector <float>& Coef);
   int get_info_threshold() const;
   void get_info_dir_threshold(std::vector <int>& dSup) const;
   void threshold (float LambdaSigma);
   void set_level_max_snr (int level);
   void set_verbose (const Bool Flag);
   void set_info (const Bool Flag);
   void set_detect_only_positive (const Bool Flag);
private:
   void comp_level( Bool ComputeSNR );
private:
   WaveletPos _Pos;
   int _DirNumber;
   int _LevelNb;
   int _levelMaxSNR;
   int _InfoThreshold;
   vector<int> _InfoDirThreshold;
   Bool _Overlap;
   Bool _CteFalseDetect;
   Bool _ImposedMaxSnrLevel;
   Bool _DetectOnlyPositive;
   float _SigmaImag;
   ofstream* _outFile;
   Bool _Verbose;
   Bool _Info;
   std::vector < std::vector <InfoCoef> > _Tree;
   std::vector < std::vector <float> > _LocalSNR;
   std::vector <float> _MaxSNR;
   void trace_tree (ofstream& of);
   void info ();
   int get_num_dir (int level, int direction) const;
   int get_level_max_snr () const;
   float get_sigma (int level, int scale, int direction) const;
   float get_current_snr (float absCoef, float normSigma, float imagSigma, 
                          int dirNumber) const;
   float nb_direction_at_direction_number( int dirNum ) const;
   Bool is_only_positive( int BlockId, int CoefNumber ) const;
   void file_info1 (int nbBloc, int nbCoef, float LambdaSigma) const;
   void file_info2 (int currentBloc, float LambdaSigma) const;
   void file_info3 (string txt, int currentBloc, int nbCoef) const;
   void file_info4 ();
};




class Direction {

public:
   Direction ();
   virtual void init (int Id, float DeltaPhi, float AngleMin, float AngleMax)=0;
   virtual void comp_support (int Nl, int Nc)=0; 
   virtual Bool is_point_in_support (float x, float y, float& coef)=0;
   virtual Bool is_point_in_overlap (float x, float y, float& coef)=0;//to be removed
   virtual float func_overlap (double phi, float phiMin, float phiMax)=0;//to be removed
   virtual float func (double phi, float phiMin, float phiMax)=0;
   virtual void info()=0;
   void write_support (string fileName, int numDir);
   void write_coef (string fileName, int numDir);
   void set_overlap (Bool Flag);
   Iint* get_support ();
   Ifloat* get_coef ();
   
   void set_verbose (const Bool Flag);
   
protected:
   int   _NlIma;
   int   _NcIma;
   int   _IdDirection;
   double _PhiMin1, _PhiMin2;
   double _PhiMax1, _PhiMax2;
   double _PhiMin1Overlap, _PhiMin2Overlap;
   double _PhiMax1Overlap, _PhiMax2Overlap;
   Bool  _Overlap;
   float _DeltaPhi;
   Iint* _Support;
   Ifloat* _Coef;
   Bool   _Verbose;
};


class Dir2 : public Direction {

public:
   Dir2 ();
   void init (int Id, float DeltaPhi, float AngleMin, float AngleMax);
   void comp_support (int Nl, int Nc); 
   Bool is_point_in_support (float x, float y, float& coef);
   Bool is_point_in_overlap (float x, float y, float& coef);//to be removed
   float func_overlap (double phi, float phiMin, float phiMax);//to be removed
   float func (double phi, float phiMin, float phiMax);
   void info();
};

class Dir1 : public Direction {

public:
   Dir1 ();
   void init (int Id, float DeltaPhi, float AngleMin, float AngleMax);
   void comp_support (int Nl, int Nc); 
   Bool is_point_in_support (float x, float y, float& coef);
   Bool is_point_in_overlap (float x, float y, float& coef);//to be removed
   float func_overlap (double phi, float phiMin, float phiMax);//to be removed
   float func (double phi, float phiMin, float phiMax);
   float support_normalisation();
   void info();
};


class DirManager {

public:

   enum DirectionType { STANDARD, NEW };
    
   DirManager();
   void init_direction( int DirectionNumber, Bool Overlap, DirectionType DirectionType,
                        Bool IncreaseDirNumber, int NlIma, int NcIma, int NbPlan);   
   
   int get_direction_number (int scale) const;
   int get_direction_number () const;
   int get_num_dir (int scale, int direction) const;
   int get_num_plane (int scale, int direction) const;
   
   Ifloat* get_coef( int NumDir );
   Iint* get_support( int NumDir );   
   float support_normalisation();
   DirectionType get_dir_type() const;
   
   void set_verbose (const Bool Flag);
   void set_debug (const Bool Flag);
   Bool get_overlap () const;
      
private:
   Direction*  _Direction;
   int         _DirectionNumber;
   int         _NlIma;
   int         _NcIma;
   int         _NbPlan;
   Bool        _Overlap;
   DirectionType _DirectionType;
   Bool        _IncreaseDirNumber;
   Bool        _TestDirection;
   Bool        _Verbose;
   Bool        _Debug;
};


class MR_Sat {

public:

         // ctor, dtor
         //
         MR_Sat ();
         ~MR_Sat ();

         void alloc (int Nl, int Nc, int NbrPlan, int DirectionNumber, 
	             Bool Overlap, DirManager::DirectionType DirectionType );
         void free();

         void transform (Ifloat& Imag);
	 void recons (Ifloat& Imag);
	 void threshold (float LambdaSigma, Bool Filter = True);
         void set_simulated_sigma() const;

         void trace() const;
	 void set_verbose (const Bool Flag);
	 void set_debug (const Bool Flag);
	 void set_border (const type_border Bord);
	 void set_norm (const sb_type_norm Norm);
	 void set_comp_sigma (const Bool Flag);
	 void set_suppress_pos (const Bool Flag);
	 void set_suppress_last_scale (const Bool Flag);
	 void set_suppress_isol_pixel (const Bool Flag);
         void set_imposed_level_of_max_snr (const int level);
         void set_sigma_image (const float SigmaImag);
         void set_write_support (const Bool Flag);
         void set_write_numsim (const Bool Flag);
         void set_recons_in_fourier (const Bool Flag);
         void set_increase_direction_number (const Bool Flag);
	 void set_cte_false_detect (const Bool Flag);
	 void set_write_norm_transf (const Bool Flag);
	 void set_test_direction (const Bool Flag);
         void write_transf (string fileName);
         void set_detect_only_positive (const Bool Flag);
         

private :

         void fft_transform (Ifloat& Imag);
         void fft_recons (Ifloat& Imag);
         void write_fft_transf (string fileName);
         void complex_info (Icomplex_f& cmpVect, int scale=-1, int dir=-1);
	 
	 void direction_transform ();
         //void init_direction (Bool Overlap);   
	 void direction_recons();
	 void tree_transform(float LambdaSigma=0., Bool filter=True);
	 void kill_isol_pixel(intarray& support);
         
         //int get_direction_number (int scale) const;
         //int get_direction_number () const;
         //int get_num_dir (int scale, int direction) const;
         //int get_num_plane (int scale, int direction) const;
         
         float simulated_sigma(fltarray& data, int d) const;
	 
private:


        // attributes
        //
        int _NlIma;
        int _NcIma;

        int _NbPlan;
        int _DirectionNumber;
	int _SizeTransf;
	Bool _Verbose;
	Bool _Debug;
	Bool _CompSigma;
        Bool _SuppressPos;
        Bool _SuppressLastScale;
        Bool _SuppressIsolPix;
        Bool _IncreaseDirNumber;
	Bool _Overlap;
	Bool _CteFalseDetect;
	Bool _WriteNormTransf;
	Bool _TestDirection;
	Bool _DetectOnlyPositive;
        float _SigmaImag;
        int _ImposedSNRMaxLevel;
        Bool _WriteSupport;
        Bool _WriteNumSim;
        Bool _ReconsInFourier;
	type_border _Border;
	sb_type_norm _Norm;
	type_transform _Transform;
	
	DirManager* _DirMng;
	
	// transf
	MultiResol _RealFftTrans;
	Icomplex_f* _FftTransf;
	//Direction* _Direction;
	Ifloat* _DirTransf;
	Ifloat* _DirRecons;
        Ifloat _LastPlane;

};

/****************************************************************************/

#endif
