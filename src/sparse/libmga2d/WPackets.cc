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
#include "IM_Sigma.h" // detect_noise_sigma

#include "WPackets.h"
 
 
//-------------------------------------------------------------------------

WPACKETS_2D::WPACKETS_2D (SubBand1D &SB1D) {
//-------------------------------------------------------------------------
   _Ptr_SB1D_Line = &SB1D; 
   _Ptr_SB1D_Col = &SB1D; 
   _pWP2D = new WPTransf_2D();
}


//-------------------------------------------------------------------------

WPACKETS_2D::WPACKETS_2D (SubBand1D &SB1D_Line, SubBand1D &SB1D_Col) {
//-------------------------------------------------------------------------
   _Ptr_SB1D_Line = &SB1D_Line; 
   _Ptr_SB1D_Col = &SB1D_Col; 
   _pWP2D = new WPTransf_2D();
}

//-------------------------------------------------------------------------

WPACKETS_2D::~WPACKETS_2D () {
//-------------------------------------------------------------------------
   free();
}

 
//-------------------------------------------------------------------------
void
WPACKETS_2D::free () {
//-------------------------------------------------------------------------
   for (MapTransf::iterator it=_pWP2D->_Transf.begin();
        it!=_pWP2D->_Transf.end();++it) 
      delete (_pWP2D->_Transf[it->first]);
}

//-------------------------------------------------------------------------
void
WPACKETS_2D::alloc (int Nl, int Nc, int NbrPlan, int NbrUndecimatedScale) {
//-------------------------------------------------------------------------
   _NlIma = Nl;
   _NcIma = Nc;
   _NbPlan = NbrPlan;
   _SetInfoBand = False;
   if (NbrUndecimatedScale >= NbrPlan) NbrUndecimatedScale = NbrPlan-1;
   else if (NbrUndecimatedScale < 0) NbrUndecimatedScale = 0;
   _NbUndec = NbrUndecimatedScale;
   _CurrentScale = -1;
   _Num.alloc (_NbPlan);
   _Verbose = False;
   _WP = True;
   _SizeTransf = (int) pow (4.,_NbUndec);
   _Infos.resize(_NbPlan);
   init_struct();
   if (_Verbose) trace();
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: init_infos () {
//-------------------------------------------------------------------------
   for (int p=0;p<_NbPlan;p++)
      _Infos[p].clear();
      
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: init_struct () {
//-------------------------------------------------------------------------
   Pt leftDownCoord; leftDownCoord.x=0; leftDownCoord.y=0;
   Ifloat* Imag = new Ifloat(_NlIma,_NcIma);
   if (_NbUndec==0) _pWP2D->_Transf.insert(make_pair(0,Imag));
   else wpRecursifAlloc(*Imag, 0, 0);
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: wpRecursifAlloc (Ifloat& Plane, int IndImag, int Num) {
//-------------------------------------------------------------------------

   // decrease current plan
   _CurrentScale++;
   
   if ((_CurrentScale != _NbPlan-1) && (is_scale_undecimated (_CurrentScale))) {
    
      int aNl = Plane.nl(), aNc = Plane.nc();
      int Delta = _SizeTransf / (int) pow (4.,_CurrentScale+1);
      
      for (int k=0;k<4;k++) 
         _pWP2D->_Transf[IndImag+k*Delta] = new Ifloat(aNl, aNc);
                     
      int localNum=0;
      for (int i=0;i<3;i++) {
         localNum = ++_Num(_CurrentScale);
         if (_WP) wpRecursifAlloc (Plane, IndImag+i*Delta, localNum);
      }        
   
      localNum = ++_Num(_CurrentScale);
      wpRecursifAlloc (Plane, IndImag+3*Delta, localNum);
   }
   
   // update current plan 
   _CurrentScale--;
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: transform (Ifloat& Imag, WPTransf_2D*& TabTrans) {
//-------------------------------------------------------------------------
   if (Imag.nl() != _NlIma || Imag.nc() != _NcIma) exit(-1);
   Pt leftDownCoord; leftDownCoord.x=0; leftDownCoord.y=0; 
   int IndImage=0, Num=0; _CurrentScale = -1;   
   wpRecursifAnalysis(Imag, leftDownCoord, IndImage, Num);
   if (_SetInfoBand == False) set_info_band();
   TabTrans = _pWP2D; 
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: recons (WPTransf_2D*& TabTrans, Ifloat& Imag) {
//-------------------------------------------------------------------------
      
   if (Imag.nl() != _NlIma || Imag.nc() != _NcIma) exit(-1);
   if (_Verbose) 
      cout << "decimated scale number :  " << (_NbPlan-1)-_NbUndec 
           << ", undecimates scale number : " << _NbUndec << endl;

   _CurrentScale = _NbPlan;
   int Step = (int) pow (2., _NbUndec);
            
   for (int d=0;d<(_NbPlan-1)-_NbUndec;++d) {
   
      _CurrentScale--;

      if (_Verbose)
         cout << "recons decimated scale " << _CurrentScale 
              << " to scale " << _CurrentScale-1 << endl;
              
      MapInfo::iterator it;
      for (it = (_Infos[_CurrentScale-1]).begin(); 
           it != (_Infos[_CurrentScale-1]).end();++it)      
         wp_Recons (TabTrans, _CurrentScale, it->second, Step);
   }
   
   for (int u=0;u<_NbUndec;++u) {
    
      _CurrentScale--; 
      Step /= 2;

      if (_Verbose) 
         cout << "recons undecimated from scale " << _CurrentScale 
              << " to scale " << _CurrentScale-1 << endl;
     
      MapInfo::iterator it;
      for (it = (_Infos[_CurrentScale-1]).begin(); 
           it != (_Infos[_CurrentScale-1]).end();++it)      
         wp_Recons (TabTrans, _CurrentScale, it->second, Step);
   }
   
   for (int i=0;i<_NlIma;++i)
   for (int j=0;j<_NcIma;++j)
      Imag(i,j) = (*TabTrans->_Transf[0])(i,j);
}
 

//-------------------------------------------------------------------------
void
WPACKETS_2D:: wpRecursifAnalysis (Ifloat& Plane, const Pt& LeftDownCoord,
                                  int IndImag, int Num) {
//-------------------------------------------------------------------------

   // decrease current plan
   _CurrentScale++;
   
   if (_CurrentScale == _NbPlan-1) {
   
      if (_Verbose) {
         cout << "working on last scale (" << _CurrentScale
              << ") : plane size : (" << Plane.nl() << "*" << Plane.nc() 
	      << ")" << " plane pos : (" << LeftDownCoord.x << ","
	      << LeftDownCoord.y << ")" << endl; 
      }
      
   } else {
   
      if (_Verbose) {
         cout << "working on scale " << _CurrentScale << endl;
         cout << "current plane size : (" << Plane.nl() << "*" 
	      << Plane.nc() << ")" << endl;
      }
           
      if (is_scale_undecimated (_CurrentScale)) {
         if (_Verbose) 
	    cout << "go undecimated at scale " << _CurrentScale << endl;
         wpRecursifUndecimated (Plane, IndImag, Num);
      } else {        
         if (_Verbose) 
	    cout << "go decimated at scale " << _CurrentScale << endl;
         wpRecursifDecimated (Plane, LeftDownCoord, IndImag, Num);
      }
   }
   
   // update current plan 
   _CurrentScale--;
}

//-------------------------------------------------------------------------
void 
WPACKETS_2D:: wpRecursifDecimated (Ifloat& Plane, const Pt& LeftDownCoord, 
                                   int IndImage, int Num) {
//-------------------------------------------------------------------------

      int aNl = Plane.nl(), aNc = Plane.nc();
      SubBand2D LocalSBT(*_Ptr_SB1D_Line, *_Ptr_SB1D_Col);
      
      int Step = (int) pow (2., _NbUndec);
      _Ptr_SB1D_Line->DistPix = Step;
      
      Ifloat* provPlane = new Ifloat[4];
      provPlane[0].alloc((aNl+1)/2,aNc/2);
      provPlane[1].alloc(aNl/2,(aNc+1)/2);
      provPlane[2].alloc(aNl/2,aNc/2);
      provPlane[3].alloc((aNl+1)/2,(aNc+1)/2);
      
      LocalSBT.transform2d (Plane, False, &(provPlane[0]), &(provPlane[1]),  
	                             &(provPlane[2]), &(provPlane[3]));
                       
      save_res_dec (provPlane, LeftDownCoord, IndImage);
      
      Pt sizeBlock = {aNl, aNc};
      InfoBlock parentBlock = {IndImage, Num, LeftDownCoord, sizeBlock};
      vector<InfoBlock> childBlock;
   
      Pt localLeftDown; sizeBlock.x = (aNl+1)/2; sizeBlock.y =  aNc/2;
      localLeftDown.x = LeftDownCoord.x+0;
      localLeftDown.y = LeftDownCoord.y+provPlane[3].nc(); 
      int localNum = ++_Num(_CurrentScale);
      if (_WP) wpRecursifAnalysis (provPlane[0], localLeftDown, IndImage, localNum); 
      InfoBlock  child1 =  {IndImage, localNum, localLeftDown, sizeBlock};
      childBlock.push_back(child1);
      
      sizeBlock.x = aNl/2;  sizeBlock.y = (aNc+1)/2;
      localLeftDown.x = LeftDownCoord.x+provPlane[3].nl();
      localLeftDown.y = LeftDownCoord.y+0;     
      localNum = ++_Num(_CurrentScale);
      if (_WP) wpRecursifAnalysis (provPlane[1], localLeftDown, IndImage, localNum); 
      InfoBlock  child2 =  {IndImage, localNum, localLeftDown, sizeBlock};
      childBlock.push_back(child2);
         
      sizeBlock.x = aNl/2; sizeBlock.y = aNc/2;
      localLeftDown.x = LeftDownCoord.x+provPlane[3].nl();
      localLeftDown.y = LeftDownCoord.y+provPlane[3].nc();     
      localNum = ++_Num(_CurrentScale);
      if (_WP) wpRecursifAnalysis (provPlane[2], localLeftDown, IndImage, localNum);
      InfoBlock  child3 =  {IndImage, localNum, localLeftDown, sizeBlock};
      childBlock.push_back(child3);
          
      sizeBlock.x = (aNl+1)/2; sizeBlock.y = (aNc+1)/2;
      localLeftDown.x = LeftDownCoord.x+0;
      localLeftDown.y = LeftDownCoord.y+0;     
      localNum = ++_Num(_CurrentScale);
      InfoBlock  child4 =  {IndImage, localNum, localLeftDown, sizeBlock};
      childBlock.push_back(child4);
      wpRecursifAnalysis (provPlane[3], localLeftDown, IndImage, localNum);
      
      if (_SetInfoBand == False) set_info (parentBlock, childBlock);
      if (_Verbose) 
         trace_info_block (get_info_block (_CurrentScale,parentBlock), 
	                   _CurrentScale);
      delete [] (provPlane);
}   


//-------------------------------------------------------------------------
void 
WPACKETS_2D:: wpRecursifUndecimated (Ifloat& Plane, int IndImage, int Num) {
//-------------------------------------------------------------------------
    
     int aNl = Plane.nl(), aNc = Plane.nc();
    // SubBand2D LocalSBT(*_Ptr_SB1D);
     SubBand2D LocalSBT(*_Ptr_SB1D_Line, *_Ptr_SB1D_Col);

     Ifloat* provPlane = new Ifloat[4];
     for (int i=0;i<4;++i) provPlane[i].alloc(aNl,aNc);
     
     int Step = (int) pow (2.,_CurrentScale);
     
     LocalSBT.transform2d (Plane, (provPlane[0]), (provPlane[1]),  
	                          (provPlane[2]), (provPlane[3]), Step);
                                     
     save_res_undec (provPlane, IndImage);
     
     int Delta = _SizeTransf / (int) pow (4.,_CurrentScale+1);
     if (_Verbose) {
        cout << "Delta = " << Delta << " at scale " << _CurrentScale << endl;
        cout << "indice image : ";
        for (int i=0;i<4;++i) cout << IndImage+i*Delta << " ";
        cout << endl; 
     }
           
     Pt leftDownCoord = {0, 0}; Pt sizeBlock = {aNl, aNc};
     InfoBlock parentBlock = {IndImage, Num, leftDownCoord, sizeBlock};
     vector<InfoBlock> childBlock;
     
     int localNum=0;
     for (int i=0;i<3;i++) {
        localNum = ++_Num(_CurrentScale);
        if (_WP) wpRecursifAnalysis (provPlane[i], leftDownCoord, 
	                             IndImage+i*Delta, localNum);
        InfoBlock  child =  {IndImage+i*Delta, localNum, leftDownCoord, sizeBlock};
        childBlock.push_back(child);
     }        
   
     localNum = ++_Num(_CurrentScale);
     wpRecursifAnalysis (provPlane[3], leftDownCoord, IndImage+3*Delta, localNum);
     InfoBlock  child3 =  {IndImage+3*Delta, localNum, leftDownCoord, sizeBlock};
     childBlock.push_back(child3);        
     
     if (_SetInfoBand == False) set_info (parentBlock, childBlock);
     if (_Verbose) 
        trace_info_block (get_info_block (_CurrentScale,parentBlock), 
	                  _CurrentScale);
     delete [] (provPlane);
}  

//-------------------------------------------------------------------------
void
WPACKETS_2D:: trace () const {
//-------------------------------------------------------------------------
   cout << "Wavelet Packets classe attributes" << endl;
   cout << "=================================" << endl;
   cout << " - plan number : " << _NbPlan << endl;
   cout << " - size image : (" << _NlIma << "*" << _NcIma << ")" << endl;
   cout << " - undec plane number : " << _NbUndec << endl;
   cout << " - size transf : " << _SizeTransf << endl;
   
}


//-------------------------------------------------------------------------
void
WPACKETS_2D:: save_res_dec (Ifloat const * Res, const Pt &LeftDown, int IndImage) {
//-------------------------------------------------------------------------
   
   for (int i=0;i<Res[0].nl();++i)
   for (int j=0;j<Res[0].nc();++j) 
      (*_pWP2D->_Transf[IndImage])(LeftDown.x+i,LeftDown.y+j+Res[3].nc())=Res[0](i,j);
      
   for (int i=0;i<Res[1].nl();++i)
   for (int j=0;j<Res[1].nc();++j) 
      (*_pWP2D->_Transf[IndImage])(LeftDown.x+Res[3].nl()+i,LeftDown.y+j)=Res[1](i,j);     
   
   for (int i=0;i<Res[2].nl();++i)
   for (int j=0;j<Res[2].nc();++j) 
      (*_pWP2D->_Transf[IndImage])(LeftDown.x+Res[3].nl()+i,LeftDown.y+j+Res[3].nc())=Res[2](i,j); 
   
   for (int i=0;i<Res[3].nl();++i)
   for (int j=0;j<Res[3].nc();++j) 
      (*_pWP2D->_Transf[IndImage])(LeftDown.x+i,LeftDown.y+j)=Res[3](i,j);
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: save_res_undec (Ifloat const * Res, int IndImage) {
//-------------------------------------------------------------------------
        
   int Delta = _SizeTransf / (int) pow (4.,_CurrentScale+1);
      
   for (int k=0;k<4;k++) {    
      for (int i=0;i<Res[k].nl();++i)
      for (int j=0;j<Res[k].nc();++j) 
         (*_pWP2D->_Transf[IndImage+k*Delta])(i,j)=Res[k](i,j);
   }  
}


//-------------------------------------------------------------------------
void
WPACKETS_2D:: write_transf(string& name) const {
//-------------------------------------------------------------------------
   for (MapTransf::const_iterator it=_pWP2D->_Transf.begin();
        it!=_pWP2D->_Transf.end();++it) {
      char apc[80]; sprintf(apc,"%s_%d.fits",name.c_str(), it->first);
      io_write_ima_float(apc, *it->second);
   }
}
 
//-------------------------------------------------------------------------
Bool
WPACKETS_2D:: is_scale_undecimated(int Scale) {
//-------------------------------------------------------------------------
   return (Scale < _NbUndec ? True : False);
}

//-------------------------------------------------------------------------
void 
WPACKETS_2D:: set_info (InfoBlock parentBlock, 
                        vector<InfoBlock>& childBlock){
//-------------------------------------------------------------------------                       

   if (_Verbose) 
      cout << "add on pair at (scale:" << _CurrentScale << ",IdParent:"
           << parentBlock.IndImage << "), coord parent : (" 
	   << parentBlock.BeginCoord.x << "," 
	   << parentBlock.BeginCoord.y << ")" << endl;
   InfoTrans inf;
   inf.Parent = parentBlock; inf.Child = childBlock;
   _Infos[_CurrentScale].insert(make_pair(parentBlock,inf));
}

//-------------------------------------------------------------------------
InfoTrans 
WPACKETS_2D:: get_info_block (int Scale, InfoBlock Parent) const {
//-------------------------------------------------------------------------                       
   return (_Infos[Scale]).find(Parent)->second;
}

//-------------------------------------------------------------------------
unsigned int 
WPACKETS_2D:: get_block_number (int Scale) const {
//-------------------------------------------------------------------------                       
   return (_Infos[Scale]).size();
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: trace_info_block (const InfoTrans& Inf, int Scale) const {
//-------------------------------------------------------------------------
   
   cout << " --> set_info on sacle " << Scale << endl;
   cout << "    ----  parent  [ id=" << Inf.Parent.IndImage  << ", num="
        << Inf.Parent.Num << ", size=("
        << "size=(" << Inf.Parent.Size.x << "*" 
        << Inf.Parent.Size.y << ")" << ", coord=(" 
        << Inf.Parent.BeginCoord.x << "," 
        << Inf.Parent.BeginCoord.y << ") ]" << endl;
   
        
   for (unsigned int i=0;i<Inf.Child.size();++i) {
      cout << "    ----  child " << i << " [ id " << Inf.Child[i].IndImage 
           << ", num:" << Inf.Child[i].Num << ", size=(" << Inf.Child[i].Size.x 
	   << "*" << Inf.Child[i].Size.y << ")" << ", coord=(" 
           << Inf.Child[i].BeginCoord.x << ","
           << Inf.Child[i].BeginCoord.y << ") ]";  
	   
      Ifloat prov(Inf.Child[i].Size.x,Inf.Child[i].Size.y);
      Ifloat* plane = _pWP2D->_Transf.find(Inf.Child[i].IndImage)->second;
      for (int x=0;x<Inf.Child[i].Size.x;++x)
      for (int y=0;y<Inf.Child[i].Size.y;++y)
         prov(x,y) = (*plane) (Inf.Child[i].BeginCoord.x+x,
	                       Inf.Child[i].BeginCoord.y+y);
       
      cout << ", mean=" << prov.mean() << ", sigma:" << prov.sigma()
           << ", min=" << prov.min() << ", max:" << prov.max() << endl;
   }
}


//-------------------------------------------------------------------------
void
WPACKETS_2D:: trace_all_blocks() const {
//-------------------------------------------------------------------------
   
   MapInfo::const_iterator it;
   cout << "information on all blocks" << endl;
   cout << "=========================" << endl;
   for (int i=0;i<_NbPlan-1;++i) {
      cout << "==> Scale " << i << endl;
      for (it = (_Infos[i]).begin(); it != (_Infos[i]).end();++it) 
         trace_info_block (get_info_block (i,it->first ), i);
   }
}



//-------------------------------------------------------------------------
void
WPACKETS_2D:: wp_Recons (WPTransf_2D*& TabTrans, int _CurrentScale, 
                         const InfoTrans& InfoBlock, int Step) {
//-------------------------------------------------------------------------
                
   if (_Verbose) trace_info_block (InfoBlock, _CurrentScale-1);
   // SubBand2D LocalSBT(*_Ptr_SB1D); 
   SubBand2D LocalSBT(*_Ptr_SB1D_Line, *_Ptr_SB1D_Col);

   if (_Verbose) cout << "Step = " << Step << endl;
   
   Ifloat* localTab =  new Ifloat[4];
   for (int i=0;i<4;++i) {
      int xSize = InfoBlock.Child[i].Size.x;
      int ySize = InfoBlock.Child[i].Size.y;
      localTab[i].alloc(xSize, ySize);
      for (int x=0;x<xSize;++x)
      for (int y=0;y<ySize;++y) 
         localTab[i](x,y) = (*TabTrans->_Transf[InfoBlock.Child[i].IndImage])
         (InfoBlock.Child[i].BeginCoord.x+x,InfoBlock.Child[i].BeginCoord.y+y);
   }
   Ifloat UpPlane (InfoBlock.Parent.Size.x, InfoBlock.Parent.Size.y);
   
   if (is_scale_undecimated(_CurrentScale-1)) {
      LocalSBT.recons2d ((localTab[0]), (localTab[1]), 
                    (localTab[2]), (localTab[3]), (UpPlane), Step);
   } else {
      _Ptr_SB1D_Line->DistPix = Step;
      LocalSBT.recons2d (UpPlane, False, 
                    &(localTab[0]), &(localTab[1]), 
                    &(localTab[2]), &(localTab[3]));                    
   }  
   
   delete [] localTab;             
                       
   for (int x=0;x<InfoBlock.Parent.Size.x;++x)       
   for (int y=0;y<InfoBlock.Parent.Size.y;++y)       
      (*TabTrans->_Transf[InfoBlock.Parent.IndImage])
             (InfoBlock.Parent.BeginCoord.x+x,
              InfoBlock.Parent.BeginCoord.y+y)=UpPlane(x,y); 
}

	
	
//-------------------------------------------------------------------------
void
WPACKETS_2D:: set_WP (const Bool Flag) {
//-------------------------------------------------------------------------
   _WP = Flag;
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: set_verbose (const Bool Flag) {
//-------------------------------------------------------------------------
   _Verbose = Flag;
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: threshold (WPTransf_2D* TabTrans, float LambdaSigma) 
{
//-------------------------------------------------------------------------
   if (_WP) wp_threshold (TabTrans, LambdaSigma);
   else ortho_threshold (TabTrans, LambdaSigma);
}

//-------------------------------------------------------------------------
void
WPACKETS_2D:: threshold (WPTransf_2D* TabTrans, float SigmaNoise, float LambdaSigma) 
{
//-------------------------------------------------------------------------
   _DataSigmaNoise = SigmaNoise;
   if (_WP) wp_threshold (TabTrans, LambdaSigma);
   else ortho_threshold (TabTrans, LambdaSigma);
}


//-------------------------------------------------------------------------
void
WPACKETS_2D:: ortho_threshold (WPTransf_2D* TabTrans, float LambdaSigma) {
//-------------------------------------------------------------------------

   for (int s=0; s<_NbPlan-1; s++) {
   
      float NSig = (s == 0) ? LambdaSigma + 1: LambdaSigma;
      if (LambdaSigma == 0) NSig=0;
      float thresholdLevel = NSig*_DataSigmaNoise;
      
      int numBlock=0;
      for (MapInfo::iterator it = (_Infos[s]).begin(); 
           it != (_Infos[s]).end();++it,++numBlock) {
                   
         InfoTrans InfoBlock = it->second;      
 
         for (int i=0;i<3;++i) {
         
            int xSize = InfoBlock.Child[i].Size.x;
            int ySize = InfoBlock.Child[i].Size.y; 
            int xBeg = InfoBlock.Child[i].BeginCoord.x;     
            int yBeg = InfoBlock.Child[i].BeginCoord.y; 
            int ind =  InfoBlock.Child[i].IndImage;
        
            if (_Verbose)
               cout << "scale:" << s << ", block;" << ind 
                    << ", pos:(" << xBeg << "," << yBeg << "), size:(" 
                    << xSize << "," << ySize << "), Threshold level:" 
                    << thresholdLevel << endl;       
                      
            for (int x=0;x<xSize;++x)
            for (int y=0;y<ySize;++y) {
         
               float coef = (*TabTrans->_Transf[ind])(xBeg+x,yBeg+y);
               (*TabTrans->_Transf[ind])(xBeg+x,yBeg+y) = 
                  (fabs(coef)<thresholdLevel ? 0.:coef);
                  
            }
         }      
      }        
   }
   string s = string("afterfilt");
   write_transf(s);
}


//-------------------------------------------------------------------------
void
WPACKETS_2D:: wp_threshold (WPTransf_2D* TabTrans, float LambdaSigma) {
//-------------------------------------------------------------------------
   
   
   float NSig = LambdaSigma; //unnecessary
   float thresholdLevel = NSig*_DataSigmaNoise;
   
   for (int s=0; s < TabTrans->nbr_band()-1; s++) {
      for (int i=0; i <  TabTrans->size_band_nl(s); i++)
      for (int j=0; j <  TabTrans->size_band_nc(s); j++) {
       
          float coef = (*TabTrans)(s,i,j);
          (*TabTrans)(s,i,j) = (fabs(coef)<thresholdLevel ? 0.:coef);          
      }
   }      
   //string s = string("afterfilt");
   //write_transf(s);
}

//-------------------------------------------------------------------------
void 
WPACKETS_2D::set_info_band() {
//-------------------------------------------------------------------------

   int scale=_NbPlan-2;
   int numBlock=0; 
    
   _SetInfoBand = True;
   int Num = 0;
   for (MapInfo::iterator it = (_Infos[scale]).begin(); 
        it != (_Infos[scale]).end();++it,++numBlock) {
       InfoTrans InfoBlock = it->second;   
       Num+=4;
   }
   _pWP2D->_WPNbrBand = Num;
   
   _pWP2D->_TabNl.alloc(_pWP2D->_WPNbrBand);
   _pWP2D->_TabNc.alloc(_pWP2D->_WPNbrBand);
   _pWP2D->_TabDepi.alloc(_pWP2D->_WPNbrBand);
   _pWP2D->_TabDepj.alloc(_pWP2D->_WPNbrBand);
   _pWP2D->_TabIndIma.alloc(_pWP2D->_WPNbrBand);
   
   scale=_NbPlan-2; Num = 0;
   for (MapInfo::iterator it = (_Infos[scale]).begin(); 
        it != (_Infos[scale]).end();++it,++numBlock) {
        
       InfoTrans InfoBlock = it->second;   
       for (int i=0;i<4;++i) {
       
          _pWP2D->_TabNl(Num) = InfoBlock.Child[i].Size.x;
          _pWP2D-> _TabNc(Num) = InfoBlock.Child[i].Size.y; 
          _pWP2D-> _TabDepi(Num) = InfoBlock.Child[i].BeginCoord.x;     
          _pWP2D->_TabDepj(Num) = InfoBlock.Child[i].BeginCoord.y; 
          _pWP2D->_TabIndIma(Num) = InfoBlock.Child[i].IndImage;
          if (_Verbose) {
             cout << "Band : " << Num << " Nl " << _pWP2D->_TabNl(Num) 
                  << " Nc = " << _pWP2D->_TabNc(Num) <<  " IndTabIma = " 
                  << _pWP2D->_TabIndIma(Num);
	     cout <<  " Depi " << _pWP2D->_TabDepi(Num) <<  " Depj " 
                  << _pWP2D->_TabDepj(Num)  << endl;
          }
	  Num ++;
       }
   }
}



//-------------------------------------------------------------------------
/************************************************************************/

//-------------------------------------------------------------------------

WPTransf_2D:: WPTransf_2D() : _WPNbrBand(-1) {
//-------------------------------------------------------------------------
}

//-------------------------------------------------------------------------
int
WPTransf_2D:: nbr_band() const {
//-------------------------------------------------------------------------
   return _WPNbrBand;
}
 
//-------------------------------------------------------------------------
int
WPTransf_2D:: size_band_nl(int b) const {
//-------------------------------------------------------------------------
   return _TabNl(b);
}

//-------------------------------------------------------------------------
int
WPTransf_2D:: size_band_nc(int b) const {
//-------------------------------------------------------------------------
   return _TabNc(b);
}

//-------------------------------------------------------------------------
float&
WPTransf_2D:: operator() (int b, int i, int j) {
//-------------------------------------------------------------------------
   
   if (    (b < 0) || (b >= _WPNbrBand) 
        || (i < 0) || (i >= size_band_nl(b)) 
        || (j < 0) || (j >= size_band_nc(b))) {
        
       cout << "Error: request coef ( " << b << "," << i << "," 
            <<  j << ")" << endl;
       cout << "       NbrBand = " << _WPNbrBand << endl;
       if ((b >= 0) && (b < _WPNbrBand))
            cout << "       Nl = " << size_band_nl(b) << endl;
       if ((b >= 0) && (b < _WPNbrBand))
            cout << "       Nc = " << size_band_nc(b) << endl;
       exit(-1);
   }
   return  (*_Transf[_TabIndIma(b)])(_TabDepi(b)+i,_TabDepj(b)+j);
}

//-------------------------------------------------------------------------
void 
WPTransf_2D:: get_band (int b, Ifloat & Band) {
//-------------------------------------------------------------------------
   
   if ( (b < 0) || (b >= _WPNbrBand)) {
       cout << "Error: bad band number " << b << endl;
       cout << "       NbrBand = " << _WPNbrBand << endl;
       exit(-1);
   }
   if  ((Band.nl() != size_band_nl(b)) || (Band.nc() != size_band_nc(b)))
      Band.alloc(size_band_nl(b), size_band_nc(b));
      
   for (int i=0; i < size_band_nl(b); i++)
   for (int j=0; j < size_band_nc(b); j++) 
      Band(i,j) =  (*this)(b,i,j);
}

//-------------------------------------------------------------------------
void 
WPTransf_2D:: put_band (int b, Ifloat & Band) {
//-------------------------------------------------------------------------
   
   if ( (b < 0) || (b >= _WPNbrBand)) {
       cout << "Error: bad band number " << b << endl;
       cout << "       NbrBand = " << _WPNbrBand << endl;
       exit(-1);
   }
   if  ((Band.nl() != size_band_nl(b)) || (Band.nc() != size_band_nc(b))) {
       cout << "Error: bad band size " << Band.nl() << " " << Band.nc() << endl;
       cout << "      namd size = " << size_band_nl(b) << " " << size_band_nc(b) << endl;
       exit(-1);
   }
   for (int i=0; i < size_band_nl(b); i++)
   for (int j=0; j < size_band_nc(b); j++) 
      (*this)(b,i,j) = Band(i,j);
      
}

//-------------------------------------------------------------------------
void 
WPTransf_2D:: info_band (int b) {
//-------------------------------------------------------------------------
   
   if ( (b < 0) || (b >= _WPNbrBand)) {
       cout << "Error: bad band number " << b << endl;
       cout << "       NbrBand = " << _WPNbrBand << endl;
       exit(-1);
   }
   Ifloat tab; get_band (b, tab);
   cout << "  band number " << b << " : size=(" << _TabNl(b) << ","
        << _TabNc(b) << "), pos=(" << _TabDepi(b) << "," << _TabDepj(b)
        << "), mean=" << tab.mean() << ", sigma=" << tab.sigma() 
        << ", min=" << tab.min() << ", max=" << tab.max() << endl;
}

//-------------------------------------------------------------------------
void 
WPTransf_2D:: info(std::string name, Bool write) {
//-------------------------------------------------------------------------
      
   cout << "info WP result  [nb band : " << _WPNbrBand << "]" << endl;
   for (int b=0;b<_WPNbrBand;b++) 
      info_band (b);
}

