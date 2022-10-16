/*******************************************************************************
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
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  
**    -----------  
**                 
******************************************************************************/


#include "IM_Obj.h" 
#include "IM_IO.h"
#include "FFTN_2D.h"

#include "MR_Sat.h"
#include "MR_Obj.h"
#include "fstream"

#include "x_sigma.cc"

extern void wave_cf_mult_strap (Icomplex_f &Ima_b, Icomplex_f &Ima_h, int s, 
              int N, float Fc, type_transform Type_Transform, Bool Down);
              
static double sym_rep( double u );
static double get_lambda_sig( double prob );
    
                       
//-------------------------------------------------------------------------

MR_Sat::MR_Sat () {
//-------------------------------------------------------------------------
   _Verbose = False;
   _Debug = False;
   _SuppressPos = True;
   _SuppressLastScale = False;
   _SuppressIsolPix = False;
   _ImposedSNRMaxLevel = -1;
   _SigmaImag = 0;
   _CompSigma = False;
   _WriteSupport = False;
   _WriteNumSim = False;
   _ReconsInFourier = False;
   _IncreaseDirNumber = False;
   _Overlap = False;
   _CteFalseDetect = False;
   _WriteNormTransf = False;
   _TestDirection = False;
   _DetectOnlyPositive = False;
}

//-------------------------------------------------------------------------

MR_Sat::~MR_Sat () {
//-------------------------------------------------------------------------
   free();
}


//-------------------------------------------------------------------------
void
MR_Sat::free () {
//-------------------------------------------------------------------------
}

//-------------------------------------------------------------------------
void
MR_Sat::alloc (int Nl, int Nc, int NbrPlan, int DirectionNumber,
               Bool Overlap, DirManager::DirectionType DirectionType) {
//-------------------------------------------------------------------------
   _NlIma = Nl;
   _NcIma = Nc;
   _NbPlan = NbrPlan;
   _DirectionNumber = DirectionNumber;
   //_Transform = TO_PAVE_FFT;
   _Transform = TO_PAVE_BSPLINE;
   //init_direction (Overlap);
   
   // direction manager
   _DirMng = new DirManager();
   _DirMng->set_debug( _Debug );
   _DirMng->init_direction( DirectionNumber, Overlap, DirectionType, 
                            _IncreaseDirNumber, Nl, Nc, NbrPlan );
   if (!_IncreaseDirNumber) _SizeTransf = (_NbPlan-1)*_DirMng->get_direction_number();
   else                     _SizeTransf = _DirMng->get_direction_number(); 
   _DirTransf = new Ifloat[_SizeTransf];
   for (int i=0;i<_SizeTransf;i++) _DirTransf[i].alloc(_NlIma,_NcIma);
   if (_Debug) cout << "  direction number = " << _SizeTransf << endl;
   
   // initialisations
   _FftTransf = new Icomplex_f[_NbPlan];
   for (int i=0;i<_NbPlan;i++) _FftTransf[i].alloc(_NlIma,_NcIma);
   _DirRecons = new Ifloat[_NbPlan];
   for (int i=0;i<_NbPlan;i++) _DirRecons[i].alloc(_NlIma,_NcIma);
   _LastPlane.alloc(_NlIma,_NcIma);
}

//-------------------------------------------------------------------------
void
MR_Sat:: transform (Ifloat& Imag) {
//-------------------------------------------------------------------------
   
   // fourrier transform
   if (_Verbose) cout << "Fourier transform..." << endl;
   if (_Debug) {
      _RealFftTrans.alloc (_NlIma, _NcIma, _NbPlan, TO_PAVE_FFT, "FftTrans");
      _RealFftTrans.Verbose = _Verbose;
      _RealFftTrans.Border = _Border;
   }
   fft_transform (Imag);   
   //if (_Debug) _RealFftTrans.write ("fft_transf");
   if (_Debug) write_fft_transf ("mr_fft");
   
   // compute direction plane
   if (_Verbose) cout << "Directional transform..." << endl;
   direction_transform();
   
}

//-------------------------------------------------------------------------
void
MR_Sat:: recons (Ifloat& Imag) {
//-------------------------------------------------------------------------

   // direction recons
   if (_Verbose) cout << "Directional recons..." << endl;
   for (int s=_NbPlan-2;s<=0;s--)
   for (int i=0;i<_NlIma;i++)
   for (int j=0;j<_NcIma;j++) 
      (_FftTransf[s])(i,j) = complex_f(0.,0.); // provisoire
   
   direction_recons();

   // fft recons
   if (_Verbose) cout << "Fourier recons..." << endl;
   fft_recons (Imag);
   
   if (!_SuppressPos) 
      for (int i=0;i<_NlIma;i++)
      for (int j=0;j<_NcIma;j++) 
         if (Imag(i,j) < 0)  Imag(i,j)=0;     
}

//-------------------------------------------------------------------------
void
MR_Sat:: threshold (float LambdaSigma, Bool Filter) {
//-------------------------------------------------------------------------
   if (_Verbose) cout << "Threshold..." << endl;
   tree_transform( LambdaSigma, Filter );
}



//-------------------------------------------------------------------------
void
MR_Sat:: trace () const {
//-------------------------------------------------------------------------
   cout << "MR Spatial adaptative transform classe attributes" << endl;
   cout << "=================================================" << endl;
   cout << " - plan number : " << _NbPlan << endl;
   cout << " - size image : (" << _NlIma << "*" << _NcIma << ")" << endl;
}
 



//-------------------------------------------------------------------------
void
MR_Sat:: set_verbose (const Bool Flag) {
//-------------------------------------------------------------------------
   _Verbose = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_debug (const Bool Flag) {
//-------------------------------------------------------------------------
   _Debug = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_comp_sigma (const Bool Flag) {
//-------------------------------------------------------------------------
   _CompSigma = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_norm (const sb_type_norm Norm) {
//-------------------------------------------------------------------------
   _Norm = Norm;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_border (const type_border Bord) {
//-------------------------------------------------------------------------
   _Border = Bord;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_suppress_pos (const Bool Flag) {
//-------------------------------------------------------------------------
   _SuppressPos = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_suppress_last_scale (const Bool Flag) {
//-------------------------------------------------------------------------
   _SuppressLastScale = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_suppress_isol_pixel (const Bool Flag) {
//-------------------------------------------------------------------------
   _SuppressIsolPix = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_imposed_level_of_max_snr (const int level) {
//-------------------------------------------------------------------------
   _ImposedSNRMaxLevel = level;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_sigma_image (const float SigmaImage) {
//-------------------------------------------------------------------------
   _SigmaImag = SigmaImage;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_write_support (const Bool Flag) {
//-------------------------------------------------------------------------
  _WriteSupport = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_write_numsim (const Bool Flag) {
//-------------------------------------------------------------------------
  _WriteNumSim = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_recons_in_fourier (const Bool Flag) {
//-------------------------------------------------------------------------
  _ReconsInFourier = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_increase_direction_number (const Bool Flag) {
//-------------------------------------------------------------------------
  _IncreaseDirNumber = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_cte_false_detect (const Bool Flag) {
//-------------------------------------------------------------------------
  _CteFalseDetect = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_write_norm_transf (const Bool Flag) {
//-------------------------------------------------------------------------
  _WriteNormTransf = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_test_direction (const Bool Flag) {
//-------------------------------------------------------------------------
  _TestDirection = Flag;
}
//-------------------------------------------------------------------------
void
MR_Sat:: set_detect_only_positive (const Bool Flag) {
//-------------------------------------------------------------------------
   _DetectOnlyPositive = Flag;
}



//-------------------------------------------------------------------------
void
MR_Sat:: fft_transform (Ifloat& Imag) {
//-------------------------------------------------------------------------
    
   Icomplex_f Ima_h_cf(_NlIma, _NcIma);
   Icomplex_f Ima_b_cf(_NlIma, _NcIma);
   //float Fc = 0.5;
   
   switch (_Transform) {
   
      /*case TO_PAVE_FFT:         
         fft2d (Imag, Ima_b_cf, 1);
         for (int s=0; s<_NbPlan-1; s ++) {
            Ima_h_cf = Ima_b_cf;
            wave_cf_mult_strap (Ima_b_cf, Ima_h_cf, s, _NlIma, Fc, 
                                _Transform, True);
            _FftTransf[s] += Ima_h_cf;
            if (_Debug) {
               fft2d (Ima_h_cf, -1);
               real( _RealFftTrans.band(s), Ima_h_cf);
	    }
         }
	 _FftTransf[_NbPlan-1] += Ima_b_cf;
	 if (_Debug) {
	    fft2d (Ima_b_cf, -1);
	    real(_RealFftTrans.scale(_NbPlan-1), Ima_b_cf);
	 }
         break;
      */
         
      case TO_PAVE_BSPLINE: {
         MultiResol mr; 
         FFTN_2D cf; cf.CenterZeroFreq=True;
         mr.alloc (_NlIma, _NcIma, _NbPlan, TO_PAVE_BSPLINE, "");
         mr.transform (Imag);
         for (int s=0; s<_NbPlan-1; s ++) {
            if (_Debug) {
               cout << "scale " << s << endl;
               cout << "=======" << endl;
            }
            Ifloat band = mr.extract_band(s);
            if (_Debug) {
               cout << "  wavelet band" << endl;
               cout << "    "; band.info();
            }
            cf.fftn2d (band, _FftTransf[s], False);
            if (_Debug) {
               cout << "  fourier tarnsform " << endl;
               complex_info (_FftTransf[s]);
            }
            
         }
         _LastPlane = mr.extract_band(_NbPlan-1);
         break;
      }
             
      default:
         fprintf (stderr,"Error: proc. fft_transform.\n");
         fprintf (stderr,"This transform is not computed this procedure\n");
         break;
   }  
}



//-------------------------------------------------------------------------
void
MR_Sat:: fft_recons (Ifloat& Imag) {
//-------------------------------------------------------------------------
    
   switch (_Transform) {
        
      case TO_PAVE_FFT: {
         int NbPlan = _SuppressLastScale == True ? _NbPlan-2: _NbPlan-1;
         for (int s=NbPlan;s>=0;s--) {
	    Ifloat inter(_NlIma, _NcIma);
	    fft2d (_FftTransf[s], -1);
	    real (inter, _FftTransf[s]);
	    Imag += inter;
	 }
         if (_Debug) {
	    Ifloat debugImag;
            debugImag = _RealFftTrans.scale(_NbPlan-1);
            for (int s=_NbPlan-2; s>=0; s--) 
               debugImag += _RealFftTrans.scale(s);
	    io_write_ima_float("debugRecons", debugImag);
	 }
         break;
      }
      
      case TO_PAVE_BSPLINE: {
         MultiResol mr; 
         FFTN_2D cf; cf.CenterZeroFreq=True;
         mr.alloc (_NlIma, _NcIma, _NbPlan, TO_PAVE_BSPLINE, "");
         for (int s=0; s<_NbPlan-1; s ++) {
            if (_Debug) {
               cout << "scale " << s << endl;
               cout << "=======" << endl;
            }
            if (_ReconsInFourier) {
 	       Ifloat inter(_NlIma, _NcIma);
               if (_Debug) {
                  cout << "  fourier tarnsform " << endl;
                  complex_info (_FftTransf[s]);
               }
               cf.fftn2d (_FftTransf[s], True);
	       real (inter, _FftTransf[s]);
               if (_Debug) {
                  cout << "  wavelet band" << endl;
                  cout << "    "; inter.info();
               }
               mr.insert_band(inter, s);
           } else {
               if (_Debug) {
                  cout << "  wavelet band" << endl;
                  cout << "    "; (_DirRecons[s]).info();
               }
               mr.insert_band(_DirRecons[s], s);
            }
         }
         mr.insert_band(_LastPlane, _NbPlan-1);
         mr.recons (Imag);
         break;
      }
	 
      default:
         fprintf (stderr,"Error: proc. fft_transform.\n");
         fprintf (stderr,"This transform is not computed this procedure\n");
         break;      
   }
}         


//-------------------------------------------------------------------------
void 
MR_Sat:: write_fft_transf (string fileName) {
//-------------------------------------------------------------------------

   for (int i=0;i<_NbPlan;i++) {
      char name[80];
      sprintf (name, "%s%d", fileName.c_str(),i); 
      Ifloat realPart(_NlIma, _NcIma); real (realPart, _FftTransf[i]);
      Ifloat imagPart(_NlIma, _NcIma); imag (imagPart, _FftTransf[i]);
      io_write_ima_complex_f ((char*)name, _FftTransf[i]);
   }
}


//-------------------------------------------------------------------------
void 
MR_Sat:: write_transf (string fileName) {
//-------------------------------------------------------------------------

   for (int s=0;s<_NbPlan-1;s++) 
   for (int d=0;d<_DirectionNumber;d++) {
      char name[80];
      sprintf (name, "%s_s%dd%d", fileName.c_str(),s,d); 
      io_write_ima_float ((char*)name, _DirTransf[_DirMng->get_num_plane(s,d)]);
   }
}
  
         

/*
//-------------------------------------------------------------------------
void 
MR_Sat:: init_direction (Bool Overlap) {
//-------------------------------------------------------------------------
   if (_Debug) std::cout << "MR_Sat:: init_direction" << std::endl;
   _Overlap = Overlap;
   int dirNumber = get_direction_number();
   //std::cout << "total number of direction : " << dirNumber << endl;
   fltarray testDir;
   if (_Debug) cout << "total number of direction : " << dirNumber << endl;
//select direction type
   _Direction = new Dir2[dirNumber];
   if (!_IncreaseDirNumber) {
      if (_TestDirection) testDir.alloc (_NlIma, _NcIma);
      //_Direction = new Direction[_DirectionNumber];
      float deltaPhi = PI / _DirectionNumber;
      float initPhi = -PI / (2.*dirNumber);
      for (int d=0;d<dirNumber;d++) {
         _Direction[d].set_overlap (Overlap);
         double decal    = (_Overlap == false ? deltaPhi : deltaPhi/2);
         double begAngle = d*decal-initPhi;
         double endAngle = d*decal-initPhi+deltaPhi;
         //std::cout << "[" << begAngle << "," << endAngle << "]" << std::endl;
         _Direction[d].init (d, deltaPhi, begAngle, endAngle); 
         _Direction[d].comp_support (_NlIma, _NcIma);
         if (_Verbose) {std::cout << "  "; _Direction[d].info();}
         if (_Debug) {
            _Direction[d].write_support("dir_support",d);
            _Direction[d].write_coef("dir",d);
         }
         if (_TestDirection) {
            Ifloat* sup = _Direction[d].get_coef();
            for (int i=0;i<_NlIma;i++)
            for (int j=0;j<_NcIma;j++) {
	       testDir(i,j) += (*sup)(i,j);
	    }
         }
      } 
   } else {
      if (_TestDirection) testDir.alloc(_NlIma, _NcIma, _NbPlan-1);
      for (int s=0;s<_NbPlan-1;s++) {
         int dirNumberAtScale = get_direction_number(s);
         float deltaPhi = PI / 
               (_Overlap == false ? dirNumberAtScale: dirNumberAtScale/2);
         float initPhi = -PI / (2.*dirNumberAtScale);
         for (int d=0;d<dirNumberAtScale;d++) {         
            int numDir = get_num_dir (s,d);
            //float deltaPhi = PI / dirNumberAtScale;
            //float initPhi = -PI / (2.*dirNumberAtScale);
            _Direction[numDir].set_overlap (Overlap);
            double decal    = (_Overlap == false ? deltaPhi : deltaPhi/2);
            double begAngle = d*decal-initPhi;
            double endAngle = d*decal-initPhi+deltaPhi;
            //std::cout << "[" << begAngle << "," << endAngle << "]" << std::endl;
            _Direction[numDir].init (d, deltaPhi, begAngle, endAngle); 
            //_Direction[numDir].init (d, deltaPhi, d*deltaPhi-initPhi, 
            //                   (d+1)*deltaPhi-initPhi); 
            _Direction[numDir].comp_support (_NlIma, _NcIma);
            if (_Verbose) {std::cout << "  "; _Direction[numDir].info();}
            if (_Debug) {
               _Direction[numDir].write_support("dir_support",numDir);
               _Direction[numDir].write_coef("dir",numDir);
            }
            if (_TestDirection) {
               Ifloat* sup = _Direction[d].get_coef();
               for (int i=0;i<_NlIma;i++)
               for (int j=0;j<_NcIma;j++) {
	          testDir(i,j,s) += (*sup)(i,j);
	       }
	    }
	 }
      }
   }
   if (_TestDirection) fits_write_fltarr ("dir_support.fits", testDir);
   if (_Verbose) std::cout << std::endl;
}
*/
//-------------------------------------------------------------------------
void 
MR_Sat:: direction_transform () {
//-------------------------------------------------------------------------
   if (_Debug) std::cout << "MR_Sat:: direction_transform" << std::endl;
   FFTN_2D cf; cf.CenterZeroFreq=True;
   for (int s=0;s<_NbPlan-1;s++) {
      if (_Debug) {
         cout << "scale " << s << endl;
         cout << "=======" << endl;
         cout << "  number of direction : " << _DirMng->get_direction_number(s) <<endl;
         cout << "  Fourier transform" << endl;
         complex_info (_FftTransf[s]);
      }
      for (int d=0;d<_DirMng->get_direction_number(s);d++) {
         if (_Debug) {
            cout << "  direction : " << d << endl;
            cout << "  --------------" << endl;
         }
         int numPlane = _DirMng->get_num_plane(s,d);
         int numDir = _DirMng->get_num_dir(s,d);
         //Iint* support = _DirMng->get_support( numDir ); //Iint* support = _Direction[numDir].get_support();
         Ifloat* coef = _DirMng->get_coef( numDir );//Ifloat* coef = _Direction[numDir].get_coef();
         Icomplex_f prov(_NlIma, _NcIma);
         for (int i=0;i<_NlIma;i++)
         for (int j=0;j<_NcIma;j++) {
	    //if ((*support)(i,j) != 0)
            //   prov(i,j) = (_FftTransf[s])(i,j);
	    //else 
	    //   prov(i,j) = complex_f(0.,0.);
//*****************************************************************************/
//MOVE to recons
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //prov(i,j) = (_FftTransf[s])(i,j) * (*coef)(i,j);
	    if ((*coef)(i,j) == 0) 
               prov(i,j) = complex_f(0.,0.);
	    else
	       prov(i,j) = (_FftTransf[s])(i,j);
         }
         if (_Debug) {
            cout << "    Fourier transform" << endl;
            complex_info (prov);
         }
         cf.fftn2d (prov, True);
         if (_Debug) {
            cout << "    Inv Fourier transform" << endl;
            complex_info (prov);
         }
         real (_DirTransf[numPlane], prov);
         if (_Debug) {
            std::cout << "    direction " << endl;
            std::cout << "    "; _DirTransf[numPlane].info();
         }
      }
   }  
   if (_CompSigma) tree_transform(0., False); 
}

 
//-------------------------------------------------------------------------
void 
MR_Sat:: direction_recons () {
//-------------------------------------------------------------------------
   if (_Debug) std::cout << "MR_Sat:: direction_recons" << std::endl;
  
   FFTN_2D cf; cf.CenterZeroFreq=True;
   for (int s=_NbPlan-2;s>=0;s--) {
      (_FftTransf[s]).init(complex_f(0.,0.));
      if (_Debug) {
         cout << "scale " << s << endl;
         cout << "=======" << endl;
         cout << "  number of direction : " << _DirMng->get_direction_number(s) <<endl;
      }
      for (int d=0;d<_DirMng->get_direction_number(s);d++) {
         if (_Debug) {
            cout << "  direction : " << d << endl;
            cout << "  --------------" << endl;
         }
         int numPlane = _DirMng->get_num_plane(s,d);
         int numDir = _DirMng->get_num_dir(s,d);
         if (_ReconsInFourier) {
            Icomplex_f prov(_NlIma, _NcIma);
            if (_Debug) {
              std::cout << "    direction " << endl;
              std::cout << "    "; _DirTransf[numPlane].info();
            }
            for (int i=0;i<_NlIma;i++)
            for (int j=0;j<_NcIma;j++) {
               prov(i,j) = complex_f((_DirTransf[numPlane])(i,j), 0.);
            }
            cf.fftn2d (prov, False);
            if (_Debug) {
               cout << "    Fourier transform (before dir threshold)" << endl;
               complex_info (prov);
            }
            Iint* support = _DirMng->get_support( numDir ); //Iint* support = _Direction[numDir].get_support();
            Ifloat* coef = _DirMng->get_coef( numDir ); //Ifloat* coef = _Direction[numDir].get_coef();
            if ( ! _DirMng->get_overlap() ) {
	       for (int i=0;i<_NlIma;i++)
               for (int j=0;j<_NcIma;j++) {
	       
	          if( _DirMng->get_dir_type() == DirManager::STANDARD ) {
	             if ((*support)(i,j) == 0) 
                        prov(i,j) = complex_f(0.,0.);
                  }
		  if( _DirMng->get_dir_type() == DirManager::NEW ) {
	             if ((*support)(i,j) == 0) 
                        prov(i,j) = complex_f(0.,0.);
	 	  }
               }
	    } else {
	       for (int i=0;i<_NlIma;i++)
               for (int j=0;j<_NcIma;j++) {
	          prov(i,j) = prov(i,j) * (*coef)(i,j);
	          if( _DirMng->get_dir_type() == DirManager::STANDARD ) {
		  }
		  if( _DirMng->get_dir_type() == DirManager::NEW ) {
		     if ((*support)(i,j) == 0)  
		        prov(i,j) = complex_f(0.,0.);   
		  }
               }	    
	    }
            if (_Debug) {
               cout << "    Fourier transform (after dir threshold)" << endl;
               complex_info (prov);
            }
            _FftTransf[s] += prov;
         
         } else {
            _DirRecons[s] += _DirTransf[numPlane];
            if (_Debug) {
               //std::cout << "  (s=" << s << ",d=" << d << ") => "; 
               std::cout << "    direction  " << endl;;
               _DirTransf[numPlane].info();
            }  
         }
      } 
      float coefNorm = _DirMng->support_normalisation();
      if (_ReconsInFourier) {
          if ( _DirMng->get_overlap() ) {
             for (int i=0;i<_NlIma;i++)
             for (int j=0;j<_NcIma;j++) 
               (_FftTransf[s])(i,j) = (_FftTransf[s])(i,j) / complex_f(coefNorm, 0.);
	 }
         if (_Debug) {
            cout << "  Fourier transform" << endl;
            complex_info (_FftTransf[s]);
         }
      } else {
         if ( _DirMng->get_overlap() ) {
             for (int i=0;i<_NlIma;i++)
             for (int j=0;j<_NcIma;j++) 
                (_DirRecons[s])(i,j) = (_DirRecons[s])(i,j) / coefNorm;
	 }
         if (_Debug) {
            cout << "  Global direction at scale " << s << endl;
	    cout << "  ---------------------------" << endl;
            _DirRecons[s].info();
         }
      }
                   
      if (_Transform == TO_PAVE_FFT) {      
         for (int i=0;i<_NlIma;i++)
         for (int j=0;j<_NcIma;j++) {
            (_FftTransf[s])(i,j) = complex_f((_DirRecons[s])(i,j), 0.);
         }
	 fft2d (_FftTransf[s], 1);
      } 
   }    
}
         
//-------------------------------------------------------------------------
void 
MR_Sat:: complex_info (Icomplex_f& cmpVect, int scale, int dir) {
//-------------------------------------------------------------------------
   if (scale != -1) cout << "scale:" << scale;
   if (dir != -1) cout << ", dir:" << dir;
   Ifloat re(_NlIma, _NcIma);
   Ifloat im(_NlIma, _NcIma);
   cout << "    "; real (re, cmpVect); re.info("real part");
   cout << "    "; imag (im, cmpVect); im.info("imag part");

}


//-------------------------------------------------------------------------
void 
MR_Sat:: tree_transform (float LambdaSigma, Bool filter) {
//-------------------------------------------------------------------------
   if (_Debug) std::cout << "MR_Sat:: tree_transform" << std::endl;

   ofstream* outFile = (ofstream*)NULL;
   if (_Debug) outFile = new ofstream ("info.dat");
   intarray support (_NlIma,_NcIma,_NbPlan);
   for (int s=0;s<_NbPlan-1;s++) {
      fltarray coef, normCoef, sup;
      int totalDir = 2*_DirMng->get_direction_number(s)-1;
      if (_CompSigma)       coef.alloc(_NlIma, _NcIma, totalDir);
      if (_WriteNormTransf) normCoef.alloc(_NlIma, _NcIma, totalDir);
      if (_WriteSupport)    sup.alloc(_NlIma, _NcIma, totalDir);
   
      float thresholdLevel = (s == 0 ? LambdaSigma+1 : LambdaSigma);
      
      for (int i=0;i<_NlIma;i++)
      for (int j=0;j<_NcIma;j++) {   
   
         vector<float> w;
         for (int d=0;d<_DirMng->get_direction_number(s);d++) {
            w.push_back((_DirTransf[_DirMng->get_num_plane(s,d)])(i,j));
         }
         WaveletPos Pos={i,j,s};
	 Bool computeSNR = True;
	 if( !filter ) computeSNR = False; //do not compute SNR if no filter
         TreeCoef tree( Pos, w, _SigmaImag, _DirMng->get_overlap(), _CteFalseDetect, 
	                outFile, computeSNR );
         tree.set_info(_Debug);
         tree.set_detect_only_positive( _DetectOnlyPositive );
         if (_CompSigma) {
            std::vector <float> w; tree.get_all_coef(w);
            for (unsigned int d=0;d<w.size();d++)
               coef(j,i,d) = w[d];
         }
         if (_WriteNormTransf) {
           std::vector <float> w; tree.get_all_normalized_coef(w);
           for (unsigned int d=0;d<w.size();d++)
               normCoef(j,i,d) = w[d];
         }
	        
         if (filter && !_CompSigma) {
         
            if (_ImposedSNRMaxLevel >= 0) 
               tree.set_level_max_snr (_ImposedSNRMaxLevel);
            tree.threshold(thresholdLevel);
         
	    std::vector<float> w; tree.get_last_level(w);
            for (int d=0;d<_DirMng->get_direction_number(s);d++) {
               (_DirTransf[_DirMng->get_num_plane(s,d)])(i,j) = w[d];
            }
            support(j,i,s) = tree.get_info_threshold();
            
            if (_WriteSupport) {
               std::vector<int> dSup(totalDir); 
               tree.get_info_dir_threshold(dSup);
               for (int d=0;d<totalDir;d++) sup(j,i,d) = dSup[d];
            }
         }
      }
      if (_CompSigma) {
         char name[80];sprintf(name, "transf_s%d.fits", s); 
         fits_write_fltarr(name,coef); 
      }
      if (_WriteNormTransf) {
         char name[80];sprintf(name, "norm_transf_s%d.fits", s); 
         fits_write_fltarr(name,normCoef); 
      }
      if (_WriteSupport) {
         char name[80];sprintf(name, "support_s%d.fits", s); 
         fits_write_fltarr(name,sup); 
      }
   }
   
   if (_SuppressIsolPix) kill_isol_pixel (support);
   
   if (_Debug) outFile->close();
   if (_WriteSupport) fits_write_intarr("support.fits", support);
}

//-------------------------------------------------------------------------
void 
MR_Sat:: kill_isol_pixel(intarray& support) {
//-------------------------------------------------------------------------

   for (int s=0;s<_NbPlan-1;s++) {
      
      for (int i=0;i<_NcIma;i++)
      for (int j=0;j<_NlIma;j++) { 
      
         if (support(s,i,j) != 0) {
         
            if (      support(s,i+1,j) == 0  and  support(s,i-1,j) == 0
                 and  support(s,i,j+1) == 0  and  support(s,i,j-1) == 0) {
                 
               support(s,i,j) = 0;
               for (int d=0;d<_DirectionNumber;d++)
                  (_DirTransf[_DirMng->get_num_plane(s,d)])(j,i) = 0; 
            }
         }
      }
   }
}


/*
//-------------------------------------------------------------------------
int 
MR_Sat:: get_direction_number(int scale) const {
//-------------------------------------------------------------------------
   int dirNumberAtCurrentScale=0;
   if (_IncreaseDirNumber)
      dirNumberAtCurrentScale = _DirectionNumber * (int)pow (2., (scale+1)/2);
   else
      dirNumberAtCurrentScale = _DirectionNumber;
   if (_Overlap) dirNumberAtCurrentScale *= 2;
   return dirNumberAtCurrentScale;
}

//-------------------------------------------------------------------------
int 
MR_Sat:: get_direction_number() const {
//-------------------------------------------------------------------------
   int totalDifferentDirNumber=0;
   if (_IncreaseDirNumber) {
      for (int s=0;s<_NbPlan-1;s++)
         totalDifferentDirNumber += get_direction_number(s);
   } else totalDifferentDirNumber = get_direction_number(0); //same number at all scale
   return totalDifferentDirNumber;
}


//-------------------------------------------------------------------------
int 
MR_Sat:: get_num_dir (int scale, int direction) const {
//-------------------------------------------------------------------------
   if (_IncreaseDirNumber) {
      int dirNumber=0;
      if (scale == 0) return direction;
      for (int s=0;s<scale;s++)
         dirNumber += get_direction_number(s);
      dirNumber += direction;
      return dirNumber;
   }
   return direction;
}

//-------------------------------------------------------------------------
int 
MR_Sat:: get_num_plane (int scale, int direction) const {
//-------------------------------------------------------------------------
   if (_IncreaseDirNumber) {
      int dirNumber=0;
      if (scale == 0) return direction;
      for (int s=0;s<scale;s++)
         dirNumber += get_direction_number(s);
      dirNumber += direction;
      return dirNumber;
   }
   return scale*get_direction_number(0) + direction;
}
*/

//-------------------------------------------------------------------------
void 
MR_Sat:: set_simulated_sigma () const {
//-------------------------------------------------------------------------

   if (_Debug )
      std::cout << "set simulated sigma" << std::endl;
      
   for (int s=0; s<_NbPlan-1;s++) {
      fltarray data;
      char name[80];sprintf(name, "transf_s%d.fits", s); 
      fits_read_fltarr(name, data);
      int d;
      int ndDir = data.nz();
      fltarray sigmaMap( ndDir );
      if ( ! _DirMng->get_overlap() ) {
         switch (ndDir) {
            case 7 : 
               for (d=0;d<ndDir;d++) x_sigma_4d[s][d] = simulated_sigma(data, d);
               break;
            case 15 : 
               for (d=0;d<ndDir;d++) x_sigma_8d[s][d] = simulated_sigma(data, d);
               break;
            case 31 : 
               for (d=0;d<ndDir;d++) x_sigma_16d[s][d] = simulated_sigma(data, d);
               break;
            case 63:
               for (d=0;d<ndDir;d++) x_sigma_32d[s][d] = simulated_sigma(data, d);
               break;
            default : 
               std::cout << "simualtion for " << ndDir << " directions ";
               std::cout << "not yet implementd \n" << std::endl;
               exit(-1);
         }
      } else {
         switch (ndDir) {
            case 7 :
               for (d=0;d<ndDir;d++) x_sigma_overlap_4d[s][d] = simulated_sigma(data, d);
               break;
            case 15 :
               for (d=0;d<ndDir;d++) x_sigma_overlap_8d[s][d] = simulated_sigma(data, d);
               break;
            case 31 :
               for (d=0;d<ndDir;d++) x_sigma_overlap_16d[s][d] = simulated_sigma(data, d);
               break;
            case 63:
               for (d=0;d<ndDir;d++) x_sigma_overlap_32d[s][d] = simulated_sigma(data, d);
               break;
           default : 
               std::cout << "simualtion for " << ndDir << " directions ";
               std::cout << "not yet implementd \n" << std::endl;
               exit(-1);
          }
      } 
      if (_WriteNumSim) {
         for( d=0 ; d<ndDir ; d++ ) sigmaMap(d) = simulated_sigma( data, d );
         char name[80];sprintf( name, "numsim_sigma_s%d.fits", s ); 
         fits_write_fltarr( name, sigmaMap ); 
         sprintf( name, "numsim_transf_s%d.fits", s ); 
         fits_write_fltarr( name, data );
      }
   }
}

//-------------------------------------------------------------------------
float 
MR_Sat:: simulated_sigma(fltarray& data, int d) const {
//-------------------------------------------------------------------------

   int nx=data.nx(), ny=data.ny();
   fltarray prov (nx, ny);
   for(int i=0; i<nx; i++)
   for(int j=0; j<ny; j++)  prov(i,j) = data(i,j,d);
   return prov.sigma();
}






//=========================================================================
//*************************************************************************
//
// ***************      Direction Manager  class     **********************
//
//*************************************************************************
//=========================================================================

//-------------------------------------------------------------------------

DirManager:: DirManager () {
//-------------------------------------------------------------------------
   _DirectionNumber = 0;
   _NlIma = 0;
   _NcIma = 0;
   _NbPlan = 0;
   _Overlap = False;
   _IncreaseDirNumber = False;
   _TestDirection = True;
   _Verbose = False;
   _Debug = False;
}


//-------------------------------------------------------------------------
void
DirManager:: init_direction( int DirectionNumber, Bool Overlap, 
                             DirectionType DirectionType, 
			     Bool IncreaseDirNumber, int NlIma, 
			     int NcIma, int NbPlan ) {
//-------------------------------------------------------------------------
   
   if( _Debug ) std::cout << "DirManager:: init_direction" << std::endl;
   _DirectionNumber = DirectionNumber;
   _NlIma = NlIma;
   _NcIma = NcIma;
   _NbPlan = NbPlan;
   _Overlap = Overlap;
   _DirectionType = DirectionType;
   _IncreaseDirNumber = IncreaseDirNumber;
   int dirNumber = get_direction_number();
   //std::cout << "total number of direction : " << dirNumber << endl;
   fltarray testDir;
   if (_Debug) cout << "total number of direction : " << dirNumber << endl;
//select direction type
   if( _DirectionType == STANDARD ) {
      _Direction = new Dir2[dirNumber];
      if (_Debug) cout << "direction type : STANDARD" << endl;
   }
   if( _DirectionType == NEW ) {
      _Direction = new Dir1[dirNumber];
      if (_Debug) cout << "direction type : NEW" << endl;
   }
   if (!_IncreaseDirNumber) {
      if (_TestDirection) testDir.alloc (_NlIma, _NcIma);
      //_Direction = new Direction[_DirectionNumber];
      float deltaPhi = PI / _DirectionNumber;
      float initPhi = -PI / (2.*dirNumber);
      for (int d=0;d<dirNumber;d++) {
         _Direction[d].set_overlap (Overlap);
	 double decal = 0;
	 if( _DirectionType == STANDARD )
            decal = (_Overlap == false ? deltaPhi : deltaPhi/2);
	 if( _DirectionType == NEW )  
	    decal = deltaPhi;
         double begAngle = d*decal-initPhi;
	 double endAngle = 0.;
	 if( _DirectionType == STANDARD )
            endAngle = d*decal-initPhi+deltaPhi;
	 if( _DirectionType == NEW )  
	    endAngle = (d+1)*deltaPhi-initPhi;   
         //std::cout << "[" << begAngle << "," << endAngle << "]" << std::endl;
         _Direction[d].init (d, deltaPhi, begAngle, endAngle); 
         _Direction[d].comp_support (_NlIma, _NcIma);
         if (_Verbose) {std::cout << "  "; _Direction[d].info();}
         if (_Debug) {
            _Direction[d].write_support("dir_support",d);
            _Direction[d].write_coef("dir",d);
         }
         if (_TestDirection) {
            Ifloat* sup = _Direction[d].get_coef();
            for (int i=0;i<_NlIma;i++)
            for (int j=0;j<_NcIma;j++) {
	       testDir(i,j) += (*sup)(i,j);
	    }
         }
      } 
   } else {
      if (_TestDirection) testDir.alloc(_NlIma, _NcIma, _NbPlan-1);
      for (int s=0;s<_NbPlan-1;s++) {
         int dirNumberAtScale = get_direction_number(s);
         float deltaPhi = 0.;
	 if( _DirectionType == STANDARD )
	    deltaPhi = PI / 
	       (_Overlap == false ? dirNumberAtScale: dirNumberAtScale/2);
	 if( _DirectionType == NEW )  
	    deltaPhi = PI / dirNumberAtScale;
         float initPhi = -PI / (2.*dirNumberAtScale);
         for (int d=0;d<dirNumberAtScale;d++) {         
            int numDir = get_num_dir (s,d);
            //float deltaPhi = PI / dirNumberAtScale;
            //float initPhi = -PI / (2.*dirNumberAtScale);
            _Direction[numDir].set_overlap (Overlap);
	    double decal = 0;
	    if( _DirectionType == STANDARD )
               decal = (_Overlap == false ? deltaPhi : deltaPhi/2);
	    if( _DirectionType == NEW )  
	       decal = deltaPhi;
            double begAngle = d*decal-initPhi;
	    double endAngle = 0.;
	    if( _DirectionType == STANDARD )
               endAngle = d*decal-initPhi+deltaPhi;
	    if( _DirectionType == NEW )  
	       endAngle = (d+1)*deltaPhi-initPhi;   
            //std::cout << "[" << begAngle << "," << endAngle << "]" << std::endl;
            _Direction[numDir].init (d, deltaPhi, begAngle, endAngle); 
            //_Direction[numDir].init (d, deltaPhi, d*deltaPhi-initPhi, 
            //                   (d+1)*deltaPhi-initPhi); 
            _Direction[numDir].comp_support (_NlIma, _NcIma);
            if (_Verbose) {std::cout << "  "; _Direction[numDir].info();}
            if (_Debug) {
               _Direction[numDir].write_support("dir_support",numDir);
               _Direction[numDir].write_coef("dir",numDir);
            }
            if (_TestDirection) {
               Ifloat* sup = _Direction[d].get_coef();
               for (int i=0;i<_NlIma;i++)
               for (int j=0;j<_NcIma;j++) {
	          testDir(i,j,s) += (*sup)(i,j);
	       }
	    }
	 }
      }
   }
   if (_TestDirection) fits_write_fltarr ("dir_support.fits", testDir);
   if (_Verbose) std::cout << std::endl;
}


//-------------------------------------------------------------------------
Ifloat*
DirManager:: get_coef( int NumDir ) {
//-------------------------------------------------------------------------
   return _Direction[NumDir].get_coef();
}

//-------------------------------------------------------------------------
float 
DirManager:: support_normalisation () {
//-------------------------------------------------------------------------
   if( _DirectionType == STANDARD ) return 1.0;
   if( _DirectionType == NEW ) return 2.0;
   cout << "undefined overlap type in support_normalisation" << endl;
   exit(-1);
   return 0;
}

//-------------------------------------------------------------------------
Iint*
DirManager:: get_support( int NumDir ) {
//-------------------------------------------------------------------------
   return _Direction[NumDir].get_support();
}

//-------------------------------------------------------------------------
DirManager::DirectionType 
DirManager:: get_dir_type() const {
//-------------------------------------------------------------------------
   return _DirectionType;
}

//-------------------------------------------------------------------------
int 
DirManager:: get_direction_number( int scale ) const {
//-------------------------------------------------------------------------
   int dirNumberAtCurrentScale=0;
   if (_IncreaseDirNumber)
      dirNumberAtCurrentScale = _DirectionNumber * (int)pow (2., (scale+1)/2);
   else
      dirNumberAtCurrentScale = _DirectionNumber;
   if (_Overlap) {
      if( _DirectionType == STANDARD) dirNumberAtCurrentScale *= 2;
      if( _DirectionType == NEW) dirNumberAtCurrentScale *= 1;
   }
   return dirNumberAtCurrentScale;
}

//-------------------------------------------------------------------------
int 
DirManager:: get_direction_number() const {
//-------------------------------------------------------------------------
   int totalDifferentDirNumber=0;
   if (_IncreaseDirNumber) {
      for (int s=0;s<_NbPlan-1;s++)
         totalDifferentDirNumber += get_direction_number(s);
   } else 
      totalDifferentDirNumber = get_direction_number(0); //same number at all scale
   return totalDifferentDirNumber;
}


//-------------------------------------------------------------------------
int 
DirManager:: get_num_dir (int scale, int direction) const {
//-------------------------------------------------------------------------
   if (_IncreaseDirNumber) {
      int dirNumber=0;
      if (scale == 0) return direction;
      for (int s=0;s<scale;s++)
         dirNumber += get_direction_number(s);
      dirNumber += direction;
      return dirNumber;
   }
   return direction;
}

//-------------------------------------------------------------------------
int 
DirManager:: get_num_plane (int scale, int direction) const {
//-------------------------------------------------------------------------
   if (_IncreaseDirNumber) {
      int dirNumber=0;
      if (scale == 0) return direction;
      for (int s=0;s<scale;s++)
         dirNumber += get_direction_number(s);
      dirNumber += direction;
      return dirNumber;
   }
   return scale*get_direction_number(0) + direction;
}

//-------------------------------------------------------------------------
void
DirManager:: set_verbose (const Bool Flag) {
//-------------------------------------------------------------------------
   _Verbose = Flag;
}
//-------------------------------------------------------------------------
void
DirManager:: set_debug (const Bool Flag) {
//-------------------------------------------------------------------------
   _Debug = Flag;
}

//-------------------------------------------------------------------------
Bool
DirManager:: get_overlap() const {
//-------------------------------------------------------------------------
   return _Overlap;
}


//=========================================================================
//*************************************************************************
//
// ***************           Direction class         **********************
//
//*************************************************************************
//=========================================================================


//-------------------------------------------------------------------------

Direction:: Direction () : _IdDirection(-1), _PhiMin1(0.), _PhiMin2(0.),
                           _PhiMax1(0.), _PhiMax2(0.), _Overlap(False) {
//-------------------------------------------------------------------------
   _Verbose = False;
}
    

//-------------------------------------------------------------------------
void 
Direction:: write_support (string fileName, int numDir) {
//-------------------------------------------------------------------------
   char name[80];
   sprintf (name, "%s%d_%d", fileName.c_str(), numDir, _IdDirection); 
   io_write_ima_int ((char*)name, *_Support);  
}

//-------------------------------------------------------------------------
void 
Direction:: write_coef (string fileName, int numDir) {
//-------------------------------------------------------------------------
   char name[80];
   sprintf (name, "%s%d_%d", fileName.c_str(), numDir, _IdDirection); 
   io_write_ima_float ((char*)name, *_Coef);  
}

//-------------------------------------------------------------------------
void 
Direction:: set_overlap (Bool Flag) {
//-------------------------------------------------------------------------
   _Overlap = Flag;
   if (_Overlap && _Verbose) std::cout << "Overlap " << endl;
}

//-------------------------------------------------------------------------
Iint*
Direction:: get_support () {
//-------------------------------------------------------------------------
   return _Support;
}

//-------------------------------------------------------------------------
Ifloat*
Direction:: get_coef () {
//-------------------------------------------------------------------------
   return _Coef;
}

//-------------------------------------------------------------------------
void
Direction:: set_verbose (const Bool Flag) {
//-------------------------------------------------------------------------
   _Verbose = Flag;
}



//=========================================================================
//*************************************************************************
//
// ***************           Direction class         **********************
//
//*************************************************************************
//=========================================================================


//-------------------------------------------------------------------------

Dir1:: Dir1 () : Direction() {
//-------------------------------------------------------------------------
}


//-------------------------------------------------------------------------
void 
Dir1:: init (int Id, float DeltaPhi, float AngleMin, float AngleMax) {
//-------------------------------------------------------------------------
   _IdDirection = Id;
   _DeltaPhi = DeltaPhi;
   _PhiMin1 = AngleMin;
   _PhiMax1 = AngleMax;
   if (_PhiMax1 >= PI) _PhiMax1 -= (2*PI);
   _PhiMin2 = _PhiMin1-PI;
   _PhiMax2 = _PhiMax1-PI;
   if (_PhiMax2 <= -PI) _PhiMax2 += (2*PI);
   if (_Overlap) {
      _PhiMin1Overlap = _PhiMin1 - _DeltaPhi;
      _PhiMax1Overlap = _PhiMax1 + _DeltaPhi;
      if (_PhiMax1Overlap >= PI) _PhiMax1Overlap -= (2*PI);
      _PhiMin2Overlap = _PhiMin2 - _DeltaPhi;
      _PhiMax2Overlap = _PhiMax2 + _DeltaPhi;
      if (_PhiMin2Overlap <= -PI) _PhiMin2Overlap += (2*PI);
      
   }
} 

//-------------------------------------------------------------------------
Bool 
Dir1:: is_point_in_support (float x, float y, float& coef) {
//-------------------------------------------------------------------------

   double phi = atan2 ((double)y, (double)x);
   coef = 0;
   //if (   (phi >=0 &&  (_PhiMin <= phi    && phi < _PhiMax)) 
   //    || (phi < 0 &&  (_PhiMin <= phi+PI && phi+PI < _PhiMax)))
   if (    (_PhiMin1 < _PhiMax1)
        && (   (_PhiMin1 <= phi    && phi < _PhiMax1)
            || (_PhiMin2 <= phi    && phi < _PhiMax2))) {
      coef = func( phi, _PhiMin1, _PhiMax1 );
      return True;
   }
   if (    (_PhiMin1 > _PhiMax1) 
        && (   (_PhiMin2 <= phi    && phi < _PhiMax2)
            || (_PhiMin1 <= phi    && phi <= PI)
            || (-PI < phi    && phi <_PhiMax1))) {
      coef = func( phi, _PhiMin1, _PhiMax1 );
      return True;
   }
   return False;
}

//-------------------------------------------------------------------------
Bool 
Dir1:: is_point_in_overlap (float x, float y, float& coef) {
//-------------------------------------------------------------------------

   double phi = atan2 ((double)y, (double)x);
   if (_PhiMin1Overlap < _PhiMin1) {
      if (_PhiMin1Overlap <= phi  && phi < _PhiMin1) {
         coef = func_overlap (phi, _PhiMin1Overlap, _PhiMin1);
         return True;
      }
      if (   (_PhiMin2Overlap < _PhiMin2) 
          && (_PhiMin2Overlap <= phi  && phi < _PhiMin2)) {
         coef = func_overlap (phi, _PhiMin2Overlap, _PhiMin2);
         return True;
      }
      if (   (_PhiMin2Overlap > _PhiMin2) 
          && (   (_PhiMin2Overlap <= phi && phi <= PI)
              || (-PI < phi    && phi < _PhiMin2))) {
         coef = func_overlap (phi, _PhiMin2Overlap-2*PI, _PhiMin2);
         return True;
      }
   }
   if (_PhiMin1Overlap > _PhiMin1) {
      if (_PhiMin2Overlap <= phi    && phi < _PhiMin2) {
         coef = func_overlap (phi, _PhiMin1Overlap, _PhiMin1);
         return True;
      }
      if (   (_PhiMin1Overlap <= phi && phi <= PI)
          || (-PI < phi    && phi < _PhiMin1)) {
         coef = func_overlap (phi, _PhiMin1Overlap-2*PI, _PhiMin1);
         return True;
      }
   }
   
   
   if (_PhiMax1Overlap > _PhiMax1) {     
      if (_PhiMax1 <= phi  && phi < _PhiMax1Overlap) {
         coef = func_overlap (phi, _PhiMax1Overlap, _PhiMax1);
         return True;
      }
      if (   (_PhiMax2Overlap > _PhiMax2) 
          && (_PhiMax2 <= phi  && phi < _PhiMax2Overlap)) {
         coef = func_overlap (phi, _PhiMax2Overlap, _PhiMax2);
         return True;
      }
      if (   (_PhiMax2Overlap < _PhiMax2) 
          && (   (_PhiMax2 <= phi && phi <= PI)
              || (-PI < phi    && phi < _PhiMax2Overlap))) {
         coef = func_overlap (phi, _PhiMax2Overlap+2*PI, _PhiMax2);
         return True;
      }
   }
   if (_PhiMax1Overlap < _PhiMax1) {
      if (_PhiMax2 <= phi && phi < _PhiMax2Overlap) {
         coef = func_overlap (phi, _PhiMax2Overlap, _PhiMax2);
         return True; 
      }
      if (   (_PhiMax1 <= phi    && phi <= PI)
          || (-PI < phi    && phi < _PhiMax1Overlap)) {
         coef = func_overlap (phi, _PhiMax1Overlap+2*PI, _PhiMax1);
         return True;
      } 
   }        
   return False;
}



//-------------------------------------------------------------------------
void 
Dir1:: info () {
//-------------------------------------------------------------------------
   if (!_Overlap)
   cout << "Direction number : " << _IdDirection << " , angle : "
        << "[" << _PhiMin1*180/PI << "," << _PhiMax1*180/PI 
        << "[ U [ " << _PhiMin2*180/PI
	<< "," << _PhiMax2*180/PI<< "[" << endl;
   if (_Overlap) {
      cout << "Direction number : " << _IdDirection << endl;;
      cout << "  first part   : "
           << "[" << _PhiMin1Overlap*180/PI << "," << _PhiMin1*180/PI << "[ U [" 
           << _PhiMin1*180/PI << "," << _PhiMax1*180/PI << "[ U ["
           << _PhiMax1*180/PI << "," << _PhiMax1Overlap*180/PI<< "[ " << endl;
      cout << "  second part : "    
           << "[" << _PhiMin2Overlap*180/PI << "," << _PhiMin2*180/PI << "[ U [" 
           << _PhiMin2*180/PI << "," << _PhiMax2*180/PI << "[ U ["
           << _PhiMax2*180/PI << "," << _PhiMax2Overlap*180/PI << "[" << endl;
   
   }
}

//-------------------------------------------------------------------------
float 
Dir1:: func_overlap (double phi, float phiMin, float phiMax) {
//-------------------------------------------------------------------------
   float a = PI / 2. / (phiMin - phiMax);
   float b = -PI * phiMax / 2. / (phiMin - phiMax); 
   return pow (cos (a*phi+b), 2.);
}

//-------------------------------------------------------------------------
float 
Dir1:: func (double phi, float phiMin, float phiMax) {
//-------------------------------------------------------------------------
   return 1.;
}

//-------------------------------------------------------------------------
void 
Dir1:: comp_support (int Nl, int Nc) {
//-------------------------------------------------------------------------
   _NlIma = Nl; _NcIma = Nc;
   _Support = new Iint(_NlIma,_NcIma); 
   _Coef = new Ifloat(_NlIma,_NcIma); 
   for (int i=0;i<_NlIma;i++)  
   for (int j=0;j<_NcIma;j++) {
      if (is_point_in_support (j-_NcIma/2, i-_NlIma/2, (*_Coef)(i,j))) {// inverse lines(i) and col(j)
         (*_Support)(i,j)=1;
      }
      else if (   _Overlap
               && is_point_in_overlap (j-_NcIma/2, i-_NlIma/2, (*_Coef)(i,j))) {
         (*_Support)(i,j)=2;
      }
   } 
   if (_Verbose) {
      int nbPts=0;
      for (int i=0;i<_NlIma;i++)  
      for (int j=0;j<_NcIma;j++) 
      if ((*_Support)(i,j) == 1) nbPts++;
      std::cout << "number of point in direction " << _IdDirection  
                <<" : " << nbPts << std::endl;
   }
}





 
//=========================================================================
//*************************************************************************
//
// ***************           Direction class         **********************
//
//*************************************************************************
//=========================================================================


//-------------------------------------------------------------------------

Dir2:: Dir2 () : Direction() {
//-------------------------------------------------------------------------
}
    
//-------------------------------------------------------------------------
void 
Dir2:: init (int Id, float DeltaPhi, float AngleMin, float AngleMax) {
//-------------------------------------------------------------------------
   _IdDirection = Id;
   _DeltaPhi = DeltaPhi;
   _PhiMin1 = AngleMin;
   _PhiMax1 = AngleMax;
   if (_PhiMax1 >= PI) _PhiMax1 -= (2*PI);
   _PhiMin2 = _PhiMin1-PI;
   _PhiMax2 = _PhiMax1-PI;
   if (_PhiMax2 <= -PI) _PhiMax2 += (2*PI);
   if (_Overlap) {
   /*   _PhiMin1Overlap = _PhiMin1 - _DeltaPhi;
      _PhiMax1Overlap = _PhiMax1 + _DeltaPhi;
      if (_PhiMax1Overlap >= PI) _PhiMax1Overlap -= (2*PI);
      _PhiMin2Overlap = _PhiMin2 - _DeltaPhi;
      _PhiMax2Overlap = _PhiMax2 + _DeltaPhi;
      if (_PhiMin2Overlap <= -PI) _PhiMin2Overlap += (2*PI);
   */   
   }
} 

//-------------------------------------------------------------------------
Bool 
Dir2:: is_point_in_support (float x, float y, float& coef) {
//-------------------------------------------------------------------------

   double phi = atan2 ((double)y, (double)x);
   //if (   (phi >=0 &&  (_PhiMin <= phi    && phi < _PhiMax)) 
   //    || (phi < 0 &&  (_PhiMin <= phi+PI && phi+PI < _PhiMax)))
   /*if (    (_PhiMin1 < _PhiMax1)
        && (   (_PhiMin1 <= phi    && phi < _PhiMax1)
            || (_PhiMin2 <= phi    && phi < _PhiMax2))) 
      return True;
   if (    (_PhiMin1 > _PhiMax1) 
        && (   (_PhiMin2 <= phi    && phi < _PhiMax2)
            || (_PhiMin1 <= phi    && phi <= PI)
            || (-PI < phi    && phi <_PhiMax1))) 
      return True;
   return False;*/
   
   
   if (_PhiMin1 < _PhiMax1) {
      if (_PhiMin1 <= phi    && phi < _PhiMax1) {
         if (_Overlap) coef = func (phi, _PhiMin1, _PhiMax1);
         else coef=1.;
         return True;
      }
      if (_PhiMin2 <= phi    && phi < _PhiMax2) {
         if (_Overlap) coef = func (phi, _PhiMin2, _PhiMax2);
         else coef=1.;
         return True;
      }
   }
   if (_PhiMin1 > _PhiMax1) {
      if (_PhiMin1 <= phi    && phi <= PI) {
         if (_Overlap) coef = func (phi, _PhiMin1, _PhiMax1+2*PI);
         else coef=1.;
         return True;
      }
      if (-PI < phi    && phi <_PhiMax1) {
         if (_Overlap) coef = func (phi+2*PI, _PhiMin1, _PhiMax1+2*PI);
         else coef=1.;
         return True;
      }
      if (_PhiMin2 <= phi    && phi < _PhiMax2) {
         if (_Overlap) coef = func (phi, _PhiMin2, _PhiMax2);
         else coef=1.;
         return True;
      }
   }
   return False;
}

//-------------------------------------------------------------------------
Bool 
Dir2:: is_point_in_overlap (float x, float y, float& coef) {
//-------------------------------------------------------------------------

   double phi = atan2 ((double)y, (double)x);
   if (_PhiMin1Overlap < _PhiMin1) {
      if (_PhiMin1Overlap <= phi  && phi < _PhiMin1) {
         coef = func_overlap (phi, _PhiMin1Overlap, _PhiMin1);
         return True;
      }
      if (   (_PhiMin2Overlap < _PhiMin2) 
          && (_PhiMin2Overlap <= phi  && phi < _PhiMin2)) {
         coef = func_overlap (phi, _PhiMin2Overlap, _PhiMin2);
         return True;
      }
      if (   (_PhiMin2Overlap > _PhiMin2) 
          && (   (_PhiMin2Overlap <= phi && phi <= PI)
              || (-PI < phi    && phi < _PhiMin2))) {
         coef = func_overlap (phi, _PhiMin2Overlap-2*PI, _PhiMin2);
         return True;
      }
   }
   if (_PhiMin1Overlap > _PhiMin1) {
      if (_PhiMin2Overlap <= phi    && phi < _PhiMin2) {
         coef = func_overlap (phi, _PhiMin1Overlap, _PhiMin1);
         return True;
      }
      if (   (_PhiMin1Overlap <= phi && phi <= PI)
          || (-PI < phi    && phi < _PhiMin1)) {
         coef = func_overlap (phi, _PhiMin1Overlap-2*PI, _PhiMin1);
         return True;
      }
   }
   
   
   if (_PhiMax1Overlap > _PhiMax1) {     
      if (_PhiMax1 <= phi  && phi < _PhiMax1Overlap) {
         coef = func_overlap (phi, _PhiMax1Overlap, _PhiMax1);
         return True;
      }
      if (   (_PhiMax2Overlap > _PhiMax2) 
          && (_PhiMax2 <= phi  && phi < _PhiMax2Overlap)) {
         coef = func_overlap (phi, _PhiMax2Overlap, _PhiMax2);
         return True;
      }
      if (   (_PhiMax2Overlap < _PhiMax2) 
          && (   (_PhiMax2 <= phi && phi <= PI)
              || (-PI < phi    && phi < _PhiMax2Overlap))) {
         coef = func_overlap (phi, _PhiMax2Overlap+2*PI, _PhiMax2);
         return True;
      }
   }
   if (_PhiMax1Overlap < _PhiMax1) {
      if (_PhiMax2 <= phi && phi < _PhiMax2Overlap) {
         coef = func_overlap (phi, _PhiMax2Overlap, _PhiMax2);
         return True; 
      }
      if (   (_PhiMax1 <= phi    && phi <= PI)
          || (-PI < phi    && phi < _PhiMax1Overlap)) {
         coef = func_overlap (phi, _PhiMax1Overlap+2*PI, _PhiMax1);
         return True;
      } 
   }        
   return False;
}



//-------------------------------------------------------------------------
void 
Dir2:: info () {
//-------------------------------------------------------------------------
   //if (!_Overlap)
   cout << "Direction number : " << _IdDirection << " , angle : "
        << "[" << _PhiMin1*180/PI << "," << _PhiMax1*180/PI 
        << "[ U [ " << _PhiMin2*180/PI
	<< "," << _PhiMax2*180/PI<< "[" << endl;
   /*if (_Overlap) {
      cout << "Direction number : " << _IdDirection << endl;;
      cout << "  first part   : "
           << "[" << _PhiMin1Overlap*180/PI << "," << _PhiMin1*180/PI << "[ U [" 
           << _PhiMin1*180/PI << "," << _PhiMax1*180/PI << "[ U ["
           << _PhiMax1*180/PI << "," << _PhiMax1Overlap*180/PI<< "[ " << endl;
      cout << "  second part : "    
           << "[" << _PhiMin2Overlap*180/PI << "," << _PhiMin2*180/PI << "[ U [" 
           << _PhiMin2*180/PI << "," << _PhiMax2*180/PI << "[ U ["
           << _PhiMax2*180/PI << "," << _PhiMax2Overlap*180/PI << "[" << endl;
   }*/
}

//-------------------------------------------------------------------------
float 
Dir2:: func_overlap (double phi, float phiMin, float phiMax) {
//-------------------------------------------------------------------------
   float a = PI / 2. / (phiMin - phiMax);
   float b = -PI * phiMax / 2. / (phiMin - phiMax); 
   return pow (cos (a*phi+b), 2.);
}

//-------------------------------------------------------------------------
float 
Dir2:: func (double phi, float phiMin, float phiMax) {
//-------------------------------------------------------------------------
   if (!_Overlap) return 1.;
   else {
      double alpha = PI / (phiMax - phiMin) * (phiMax - phi);
      return pow (sin (alpha), 2.);
   }
}

//-------------------------------------------------------------------------
void 
Dir2:: comp_support (int Nl, int Nc) {
//-------------------------------------------------------------------------
   _NlIma = Nl; _NcIma = Nc;
   _Support = new Iint(_NlIma,_NcIma); 
   _Coef = new Ifloat(_NlIma,_NcIma); 
   for (int i=0;i<_NlIma;i++)  
   for (int j=0;j<_NcIma;j++) {
      if (is_point_in_support (j-_NcIma/2, i-_NlIma/2, (*_Coef)(i,j))) {// inverse lines(i) and col(j)
         (*_Support)(i,j)=1;
         //(*_Coef)(i,j)=1.;
      }
      /*else if (   _Overlap
               && is_point_in_overlap (j-_NcIma/2, i-_NlIma/2, (*_Coef)(i,j))) {
         (*_Support)(i,j)=2;
      }*/
   } 
   if (_Verbose) {
      int nbPts=0;
      for (int i=0;i<_NlIma;i++)  
      for (int j=0;j<_NcIma;j++) 
      if ((*_Support)(i,j) == 1) nbPts++;
      std::cout << "number of point in direction " << _IdDirection  
                <<" : " << nbPts << std::endl;
   }
}
 


















//=========================================================================
//*************************************************************************
//
// ***************              TreeCoef             **********************
//
//*************************************************************************
//=========================================================================


//-------------------------------------------------------------------------

TreeCoef:: TreeCoef (WaveletPos Pos, vector<float>& w, float SigmaImag,  
                     Bool Overlap, Bool CteFalseDetect, ofstream* os, 
	             Bool ComputeSNR)  
		     : _Pos(Pos), _Overlap(Overlap), 
		       _CteFalseDetect(CteFalseDetect) {
//-------------------------------------------------------------------------
   _Overlap=False;   
   _Verbose = False;
   _Info = False;
   _ImposedMaxSnrLevel = False;
   _DetectOnlyPositive = False;
   _DirNumber = w.size();
   _levelMaxSNR = -1;
   _InfoThreshold = 0;
   _SigmaImag = SigmaImag;
   _LevelNb = iround(log((double)_DirNumber)/log(2.) + 1);
   _Tree.resize(_LevelNb);
   _LocalSNR.resize(_LevelNb);
   _outFile = os;
   for (int lv=0;lv<_LevelNb;lv++) {
      (_Tree[lv]).resize (int(pow(2.0,lv)));
      (_LocalSNR[lv]).resize (int(pow(2.0,lv)));
   }
   for (int i=0;i<_DirNumber;i++) {
      (_Tree[_LevelNb-1])[i].Coef = w[i];
   }
   _MaxSNR.reserve(_LevelNb);
   for (int lv=0;lv<_LevelNb;lv++)
      _MaxSNR[lv]=0.;
 
   
   comp_level(ComputeSNR);
   if (_Verbose) {info();}
   
}


//-------------------------------------------------------------------------
void
TreeCoef:: comp_level( Bool ComputeSNR ) {
//-------------------------------------------------------------------------

   for (int lv=_LevelNb-2;lv>=0;lv--) {
      for (int i=0;i<int(pow(2.0,lv));i++) {
         (_Tree[lv][i]).Coef =   (_Tree[lv+1][2*i]).Coef 
                               + (_Tree[lv+1][2*i+1]).Coef;
			       
         if( ComputeSNR ) {
            int indDir = get_num_dir (lv, i);
//cout << "i:" << _Pos.xPos << ", j:" << _Pos.yPos << ", s:" << _Pos.scale << endl;
//cout << "lev:" << lv << ", pt:" << i << ", ind dir:" << indDir << endl;
//cout << "val:" << fabs((_Tree[lv][i]).Coef) 
//     << ", thr:" <<get_sigma(_LevelNb, _Pos.scale, indDir)  
//     << ", sigimag:" <<_SigmaImag; 
            _LocalSNR[lv][i] = get_current_snr( fabs((_Tree[lv][i]).Coef), 
                                                get_sigma(_LevelNb, _Pos.scale, indDir),
                                                _SigmaImag, lv);
//cout << "curSNR:" << currentSNR << endl;
            _MaxSNR[lv] = max(_MaxSNR[lv], _LocalSNR[lv][i]);
	 }
      }
//cout << "maxSNR:" << _MaxSNR[lv] << endl;
//cout << endl;
   }
   
   if( ComputeSNR ) {
      int lastLevel = _LevelNb-1;
      for (int i=0;i<int(pow(2.0,lastLevel));i++) {
         int indDir = get_num_dir (lastLevel, i);
//cout << "i:" << _Pos.xPos << ", j:" << _Pos.yPos << ", s:" << _Pos.scale << endl;
//cout << "lev:" << lastLevel << ", pt:" << i << ", ind dir:" << indDir << endl;
//cout << "val:" << fabs((_Tree[lastLevel][i]).Coef)
//     << ", thr:" << get_sigma(lastLevel, _Pos.scale, indDir)  
//     << ", sigimag:" <<_SigmaImag; 

         _LocalSNR[lastLevel][i] = get_current_snr( fabs((_Tree[lastLevel][i]).Coef), 
                                             get_sigma(_LevelNb, _Pos.scale, indDir),
                                             _SigmaImag, lastLevel);
//cout << "curSNR:" << currentSNR << endl;
         _MaxSNR[lastLevel] = max(_MaxSNR[lastLevel], _LocalSNR[lastLevel][i]);
      }
//cout << "maxSNR:" << _MaxSNR[lastLevel] << endl;
//cout << endl;
       
      _levelMaxSNR = get_level_max_snr();
//cout << "maxSNR:" << _levelMaxSNR << endl;  
   } 
}


//-------------------------------------------------------------------------
float TreeCoef:: 
get_current_snr (float absCoef, float normSigma, float imagSigma, 
                 int dirNumber) const {
//-------------------------------------------------------------------------
   return absCoef / normSigma / imagSigma;
}

//-------------------------------------------------------------------------
void
TreeCoef:: trace_tree (ofstream& of) {
//-------------------------------------------------------------------------

   of << "TreeCoef [i:" << _Pos.xPos << ",j:" << _Pos.yPos
            << ",s:" << _Pos.scale << "] " << std::endl;
   for (int lv=_LevelNb-1;lv>=0;lv--) { 
      of << "  level " << lv << " (nb_direction=" << pow(2.,lv) << ") : ";
      for (int i=0;i<int(pow(2.0,lv));i++) { 
         of << "[d:" << i << "]=" << (_Tree[lv][i]).Coef << " ";
      }
      of << ",  max SNR : " << _MaxSNR[lv];
      of << std::endl;
   }
   if (_levelMaxSNR != -1) {
     of << "  level of max snr is "; 
     if (!_ImposedMaxSnrLevel) of << ": " << _levelMaxSNR << std::endl;
     else  of << "imposed to level " << _levelMaxSNR << std::endl;
  }
}


//-------------------------------------------------------------------------
void
TreeCoef:: info () {
//-------------------------------------------------------------------------

   std::cout << "Pos : i:" << _Pos.xPos << ", j:" << _Pos.yPos 
             << ", s:" << _Pos.scale;
   std::cout << ", dir number : " << _DirNumber << ", level number : "
             << _LevelNb << std::endl;
}


//-------------------------------------------------------------------------
void
TreeCoef:: get_all_coef (std::vector <float>& Coef) {
//-------------------------------------------------------------------------
   Coef.clear();
   for (int lv=0;lv<_LevelNb;lv++)    
   for (int i=0;i<int(pow(2.0,lv));i++) 
      Coef.push_back ((_Tree[lv][i]).Coef);
}

//-------------------------------------------------------------------------
void
TreeCoef:: get_all_normalized_coef (std::vector <float>& NormCoef) {
//-------------------------------------------------------------------------
   NormCoef.clear();
   for (int lv=0;lv<_LevelNb;lv++)    
   for (int i=0;i<int(pow(2.0,lv));i++) {
      if ((_Tree[lv][i]).Coef < 0)
         NormCoef.push_back (-(_LocalSNR[lv][i]));
      else
         NormCoef.push_back ((_LocalSNR[lv][i]));
   }
}

//-------------------------------------------------------------------------
void
TreeCoef:: get_last_level (std::vector <float>& Coef) {
//-------------------------------------------------------------------------
   int lastLevel =  _LevelNb-1;   
   for (int i=0;i<int(pow(2.0,lastLevel));i++) 
      Coef.push_back ((_Tree[lastLevel][i]).Coef);
}

//-------------------------------------------------------------------------
int  
TreeCoef:: get_num_dir (int level, int direction) const  {
//-------------------------------------------------------------------------
   if (level == 0) return 0;
   return (int)(pow(2.,level))-1 + direction;
}

//-------------------------------------------------------------------------
int  
TreeCoef:: get_info_threshold () const  {
//-------------------------------------------------------------------------
   return _InfoThreshold;
}


//-------------------------------------------------------------------------  
void
TreeCoef::  get_info_dir_threshold(std::vector <int>& dSup) const {
//-------------------------------------------------------------------------

  for (unsigned int i=0; i<_InfoDirThreshold.size(); i++)
     dSup[_InfoDirThreshold[i]]=1;
}

//-------------------------------------------------------------------------   
float 
TreeCoef:: nb_direction_at_direction_number( int dirNum ) const {
//-------------------------------------------------------------------------
  float ret=0;
  if( dirNum == 0) ret=1;
  else
  if( dirNum == 1 || dirNum == 2 ) ret=2;
  else
  if( dirNum > 2  && dirNum < 7 ) ret=4;
  else
  if( dirNum > 6  && dirNum < 15 ) ret=8;
  else
  if( dirNum > 14  && dirNum < 31 ) ret=16;
  else
  if( dirNum > 30  && dirNum < 63 ) ret=32;
  return ret;
}

//-------------------------------------------------------------------------
void
TreeCoef:: threshold (float LambdaSigma) {
//-------------------------------------------------------------------------
   
   if (_Info) trace_tree(*_outFile);
  
   int nbBloc = (int)pow(2.,_levelMaxSNR);
   int nbCoef = (int)pow(2.,_LevelNb-(_levelMaxSNR+1));
   int lastLevel = _LevelNb-1;
   if (_Info) file_info1 (nbBloc, nbCoef, LambdaSigma);
                       
   Bool sup = False; _InfoDirThreshold.clear();
   for (int nb=0;nb<nbBloc;nb++) {
      //int indDir = get_num_dir (_levelMaxSNR, nb);
      //float currentSNR = get_current_snr( fabs(   (_Tree[_levelMaxSNR][nb]).Coef), 
      //                                    get_sigma(_LevelNb,_Pos.scale,indDir),
      //                                    _SigmaImag, _levelMaxSNR); 
//cout << "nb:" << nb << ", indDir:" << indDir << ", currentSNR:" 
//     << currentSNR << ", _LevelNb:" << _LevelNb << endl;
//cout << "val:" << fabs(   (_Tree[_levelMaxSNR][nb]).Coef) 
//     << ", sig:" << get_sigma(_LevelNb,_Pos.scale,indDir)
//     << ", sigimag:" << _SigmaImag << endl;

      float localLambdaSigma;
      if( _CteFalseDetect ) {
      
	 double prob = sym_rep( LambdaSigma );
	 float nbDir = nb_direction_at_direction_number( nb );
	 localLambdaSigma = get_lambda_sig( 1.- (1.-prob)/nbDir );
	 
      } else {
         localLambdaSigma = LambdaSigma;
      }

      float currentSNR = _LocalSNR[_levelMaxSNR][nb];
      if (_Info) {
         file_info2 (nb, localLambdaSigma); 
         file_info3 ("      in  coef : ", nb, nbCoef);
      }
     
      if (currentSNR < localLambdaSigma) {
         for (int nc=0;nc<nbCoef;nc++) {
            (_Tree[lastLevel][nb*nbCoef + nc]).Coef = 0.;   
         }
      } else {
         //if( _DetectOnlyPositive && _Tree[_levelMaxSNR][nb].Coef < 0) {
         if( _DetectOnlyPositive && !is_only_positive( nb, nbCoef ) ) {
            if (_Info) *_outFile << "      ===> DONT KEEP COEF (ONLY POSITIVE) !!!!!!!!"  
	                         << std::endl;
	    for (int nc=0;nc<nbCoef;nc++) {
               (_Tree[lastLevel][nb*nbCoef + nc]).Coef = 0.; 
	    }  
	 } else {
            _InfoThreshold = (int)pow(2.,_levelMaxSNR);
	    if (_Info) *_outFile << "      ===> KEEP COEF !!!!!!!!"  << std::endl;
            sup = True;
            _InfoDirThreshold.push_back ( iround(pow(2., _levelMaxSNR)-1+nb));
	 }
         /*if (_levelMaxSNR == _LevelNb-1)
            _InfoDirThreshold.push_back (nb);
         else {
            int iPos=0;
            for (int il=_LevelNb-1; il>=_levelMaxSNR+1; il--)
               iPos += pow(2.,il);
            _InfoDirThreshold.push_back (iPos+nb);
            std::cout << "lev:" << _levelMaxSNR << ", num:"
                      << nb << ", iPos:" << iPos+nb << std::endl;
         }*/ 
         
      }
      
      if (_Info) file_info3 ("      out coef : ", nb, nbCoef);
   } 
   if (_Info && sup) *_outFile << "    ==> in SUP" << std::endl;
   if (_Info) file_info4();
}

//-------------------------------------------------------------------------
int  
TreeCoef:: get_level_max_snr () const  {
//-------------------------------------------------------------------------
   int maxLevel=0;
   for (int lv=0;lv<_LevelNb;lv++) 
      if (_MaxSNR[lv] > _MaxSNR[maxLevel])
         maxLevel =  lv;
   return maxLevel;
}

//-------------------------------------------------------------------------
void  
TreeCoef:: set_level_max_snr (int level) {
//-------------------------------------------------------------------------
   if (level < 0 || level > _LevelNb-1) level=0;
   _levelMaxSNR = level;
   _ImposedMaxSnrLevel=True;
}

//-------------------------------------------------------------------------   
float 
TreeCoef:: get_sigma (int level, int scale, int direction) const {
//-------------------------------------------------------------------------
   
   switch (level) {

      case  2: if (_Overlap) return  x_sigma_overlap_2d[scale][direction]; 
               else return  x_sigma_2d[scale][direction];
               break;
      case  3: if (_Overlap) return  x_sigma_overlap_4d[scale][direction]; 
               else return  x_sigma_4d[scale][direction]; 
               break;
      case  4: if (_Overlap) return  x_sigma_overlap_8d[scale][direction]; 
               else return  x_sigma_8d[scale][direction]; 
               break;
      case  5: if (_Overlap) return  x_sigma_overlap_16d[scale][direction]; 
               else return x_sigma_16d[scale][direction]; 
               break;
      case  6: if (_Overlap) return  x_sigma_overlap_32d[scale][direction]; 
               else return x_sigma_32d[scale][direction]; 
               break;
      case  7: if (_Overlap) return  x_sigma_overlap_64d[scale][direction]; 
               else return x_sigma_64d[scale][direction]; 
               break;
      default :
         cout << level << " is too big!!!" << endl;
         cout << "Not yet implmented (in TreeCoef::get_sigma)" << endl;
         exit(-1);
   }
}

 

//-------------------------------------------------------------------------
void  
TreeCoef:: file_info1 (int nbBloc, int nbCoef, float LambdaSigma) const {
//-------------------------------------------------------------------------
   *_outFile << "  number of blocs at SNR max level = " << nbBloc;
   *_outFile << ", number of coef for each block = " << nbCoef << std::endl;
   *_outFile << "  threshold procedure (k=" << LambdaSigma; 
   *_outFile << "), sigma image=" << _SigmaImag << std::endl;             
}

//-------------------------------------------------------------------------
void  
TreeCoef:: file_info2 (int currentBloc, float LambdaSigma) const {
//-------------------------------------------------------------------------
   int indDir = get_num_dir (_levelMaxSNR, currentBloc);
   float coef = fabs((_Tree[_levelMaxSNR][currentBloc]).Coef);
   float sigma = get_sigma(_LevelNb,_Pos.scale,indDir);
   
   *_outFile << "    block number : " << currentBloc << ",   local SNR : " 
             << _LocalSNR[_levelMaxSNR][currentBloc] << " (compare to k=" 
	     << LambdaSigma << ") [local SNR = " << coef << "/("
	     << sigma << "*" << _SigmaImag;
   *_outFile << ")]" << std::endl; 
}

//-------------------------------------------------------------------------
void  
TreeCoef:: file_info3 (string txt, int currentBloc, int nbCoef) const {
//-------------------------------------------------------------------------
   int lastLevel = _LevelNb-1;
   *_outFile << txt << "{ ";
   for (int nc=0;nc<nbCoef;nc++) 
      *_outFile << "[" << currentBloc*nbCoef + nc << "]" 
                       << (_Tree[lastLevel][currentBloc*nbCoef + nc]).Coef 
                       << " ";   
   *_outFile << "}" << std::endl;
}

//-------------------------------------------------------------------------
void  
TreeCoef:: file_info4 () {
//-------------------------------------------------------------------------
   std::vector<float> w; get_last_level(w);
    *_outFile << "    out coef : { ";
   for (unsigned int i=0;i<w.size();i++) 
      *_outFile << "[" << i << "]" << w[i] << " ";   
   *_outFile << "}" << std::endl; 
   *_outFile << " " << std::endl; 
}
//-------------------------------------------------------------------------
void
TreeCoef:: set_verbose (const Bool Flag) {
//-------------------------------------------------------------------------
   _Verbose = Flag;
}
//-------------------------------------------------------------------------
void
TreeCoef:: set_info (const Bool Flag) {
//-------------------------------------------------------------------------
   _Info = Flag;
}

//-------------------------------------------------------------------------
void
TreeCoef:: set_detect_only_positive (const Bool Flag) {
//-------------------------------------------------------------------------
   _DetectOnlyPositive = Flag;
}

//-------------------------------------------------------------------------
Bool 
TreeCoef:: is_only_positive( int BlockId, int CoefNumber ) const {
//-------------------------------------------------------------------------
   for( int nc=0;nc<CoefNumber;nc++ ) {
      if( (_Tree[_LevelNb-1][BlockId*CoefNumber + nc]).Coef < 0. )
         return False;
   }  
   return True;
}





//=========================================================================
//*************************************************************************
//
// ***************            Static func            **********************
//
//*************************************************************************
//=========================================================================

static double repartition_func( double u ) {

   double a  =  0.2316419;
   double b1 =  0.319381530;
   double b2 = -0.356563782;
   double b3 =  1.781477937;
   double b4 = -1.821255978;
   double b5 =  1.330274429;
   
   double t = 1. / ( 1. + a*u );

   double f = b1*t+b2*pow(t,2.)+b3*pow(t,3.)+b4*pow(t,4.)+b5*pow(t,5.);
   f = 1. - 1./sqrt(2*PI)*exp(-pow(u,2.)/2.) * f;

   return f;
}


static double sym_rep( double u ) {
   return 1. - (1 - repartition_func( u ))*2.;
}

static double get_lambda_sig( double prob ) {

   double t1 = 1.;
   double prob_t1 = sym_rep( t1 ); 
   double f1 = prob_t1 - prob;

   double t3 = 10.;
   double prob_t3 = sym_rep( t3 ); 
   double f3 = prob_t3 - prob;
   
   unsigned int nbIter = 20;   
   
   for( unsigned int i=0; i<nbIter; i++ ) {
   
      double t2 = (t1+t3)/2.;
      double prob_t2 = sym_rep( t2 ); 
      double f2 = prob_t2 - prob;
   
      if( f2 >= 0 ) {
         t3 = t2;
	 f3 = f2;
      } else {
         t1 = t2;
	 f1 = f2;
      }
   }
   
   return (t1+t3)/2.;
}












