/***********************************************************
**	Copyright (C) 1999 CEA
************************************************************
**
**    UNIT
**
**    Version:  1.0
**
**    Author: 	J.L. Starck
**
**    Date: 	27/05/99
**    
**    File:  	MVM.h
**
************************************************************
**
** Multiscale Vision Model Definition File
**
**  3 classes are need
**   a structure definition: struct2d
**        set of detected wavelet coefficients at the same scale, and
**        connected.
**
**   an object definition: info_obj2d
**        set of structure connected by some relations
**
**   MR2D_Seg: detection class
**             MR2D_Seg.segment find all objects contained in the original image
**             MR2D_Seg.get_mrobj extract the wavelet coefficinet belonging 
**                                to a given object. 
**  
************************************************************/


#ifndef	_JMVM_H_
#define	_JMVM_H_

#include "IM_Math.h"
#include "NR.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "IM_VisTool.h"
#include "IM_Graphics.h"
 
extern void rec_iter_grad_conj (MultiResol& W0, Ifloat &Im_rec, int Nb_iter,
				double& Erreur, float eps);


#define NBR_MVM 3
enum t_mvm {MVM_NO,MVM_OVERLAP,MVM_RUEBIJ,MVM_SUBOBJ,MVM_OVER_AND_ALL};
#define DEF_MVM MVM_RUEBIJ

inline const char * StringMVM (t_mvm type)
{
    switch (type)
    {
    case  MVM_NO: return ("No vision model"); break;
    case  MVM_OVERLAP: return ("Blinded Objects"); break;
    case  MVM_SUBOBJ: return ("Embedded Objects"); break;
    case  MVM_OVER_AND_ALL:return ("Blinded + Embedded Objects"); break;
    case  MVM_RUEBIJ:return ("Rue-Bijaoui Vision Model for blinded + embedded Objects"); break;
    default: return ("Error: bad vision model");
    }
}

inline int pos_in_ima(int x, int Scale, set_transform Set)
{
   int Ret=x;
   if (Set == TRANSF_PYR)  Ret <<= Scale; 
   else if  ((Set == TRANSF_SEMIPYR) && (Scale > 1)) Ret <<= Scale-1; 
   return Ret;
}

inline int pos_in_transform(int x, int Scale, set_transform Set)
{
   int Ret=x;
   if (Set == TRANSF_PYR)  Ret = size_ima_resol(Scale, x);
   else if  ((Set == TRANSF_SEMIPYR) && (Scale > 1)) Ret = size_ima_resol(Scale-1, x);
   return Ret;
}

// structure (a structure is defined for each connected group 
// of pixels detected at a given scale)
// A object is contained several structures
class struct2d
{
  void reset()
  {
      NumStruct=0;
      NumConnectObj=0;
      Scale=0;
      NumScale=0;
      Depx = Depy= -1;
      Maxx= Maxy=0;
      ValMax=0;
      Sizex= Sizey=0;
      Npix=0;
      Noise=0;
      SubStruct=False;
      NSubStruct=0;
      Flux=0;
      BadStruct=False;
      Father=0;
      NbrTotalSub=0;
  }
  public:
   struct2d() {reset();}
   int Father;       // Father structure 
   Bool BadStruct;   // True if the structure is killed
   int NumStruct;        // structure number (between 1.. nbr of structures)
   int NumConnectObj; // object number (between 1.. nbr of object)
   int Scale;        // scale of the structure
   int NumScale;     // scale where the maximum of the object is found 
   int Depx;         // starting pixel of the structure
   int Depy;         // starting pixel of the structure
   int Maxx;         // position of the maximum of the structure
   int Maxy;         // position of the maximum of the structure
   float ValMax;     // maximum value of the structure
   int Sizex;        // size of the structure
   int Sizey;        // size of the structure
   double Flux;      // integrate flux
   int Npix;         // Number of pixels in the structure
   float Noise;      // noise at this position
   int NSubStruct;   // Number of sub-structures at the next scale   
   int NbrTotalSub;  // Total number of sub-structure in all scales
   Bool SubStruct;      // equal to True, if the structure is a sub-structure
   ~struct2d() {reset();}
};

class info_obj2d {
  void reset()
  {   
    DepImaY =0;
    DepImaX=0;
    EndImaX=0;   
    EndImaY=0;              
    SizeImaX=0;
    SizeImaY=0;
    RootScale=0;
    RootIndex=0;
    Max=0;
    ScaleMax=0;
    PosMaxImaX=0;
    PosMaxImaY=0;
    if (TabIndStruct.n_elem()  != 0) TabIndStruct.free();
    NbrStruct=0;
    ObjNumberInScale=0;
    SubObj=False;
    NbrSubObj=0;
    NumFatherObj=0;
   }
  public:
   info_obj2d() {reset();}
   int NbrStruct;     // number of structure in the object. 
   int DepImaX;       // starting x position in the image
   int DepImaY;       // starting y position in the image
   int EndImaX;         // ending pixel of the structure
   int EndImaY;         // ending pixel of the structure
   int SizeImaX;      // Image size in x direction
   int SizeImaY;      // Image size in y direction
   int RootScale;     // Root structure scale 
   int RootIndex;     // Root structure index
   float Max;          // maximum value of the object
   int ScaleMax;       // ScaleMax  gives the scale where the maximum is found 
                       // scales are ranged between 0 and Nbr_Plan-2
   int PosMaxImaX;     //  position X where the maximum is found
   int PosMaxImaY;     //  position Y where the maximum is found
   int ObjNumberInScale; // object number at the object scale
   Bool SubObj;      // equal to True, if the structure is the 
                     // starting point of a sub-object
   int NumFatherObj; // father object if it is a sub-object, 		     
   int NbrSubObj;    // Number of subobj included in this object
   intarray TabIndStruct;   // TabObjIndStruct[1..NStruct]
                            // gives the index in table TabStruct.
   intarray TabIndSubObj;   // SubObj index table
  ~info_obj2d() 
  {	
     reset();
  }    
};


class MR2D_Seg {
   
   MultiResol Segment;     // Multiresolution object
                           // Segment(s,i,j) gives the structure number
                           // of the wavelet coefficient (s,i,j)
                           // if Segment(s,i,j) = 0, nothing is detected
                           // at this scale
   struct2d *TabStruct; // Array of structure  
                            // TabStruct[0] is not used
                            // TabStruct[1..NbrStruct] gives information
                            // about each detected structure
   int TabNbrStruct[MAX_SCALE]; // number of structure per scale
   int NbrStruct;           // number of structures
   int NbrSubStruct;        // Number of sub-structures

   info_obj2d *TabObj; 
   info_obj2d *TabSubObj;
   int TabNbrObjPerScale[MAX_SCALE]; // number of object per scale
   intarray TabScaleObj;   //  object list per scale
                           // TabScaleObj(s,i)
			   // object number i at scale s
			   // TabObj[ TabScaleObj(s,i) ] is the 
			   // corresponding object
			   // s = 0 .. Nbr_Plan-2
			   // i = 1 .. TabNbrObjPerScale[s]
   int Nbr_Plan;           // number of scale
   int Nl;                 // number of lines 
   int Nc;                 // number of columns
   int NbrObj;              // number of objects: one object can contain
                            // several structure: NbrObj <= NbrStruct
   int NbrSubObj;           // Number of sub-obj	
   int NbrTotObj;           // Number of object + sub object		    
   intarray TabRel;         //  TabRel(i,0) = number of sub-structures
                            //                of the structure i
 			    //  TabRel(i,j) = structure number of the jth 
			    //                sub-structure
 
   void find_structures(MultiResol & MR_Data, MRNoiseModel & NoiseModel);
                            // Find the structures in the multiresolution
			    // data set
			    
   void connect_struct(); // set the relation table TabRel by
                           // connecting the structure	
			    // two structures S(j),S(j+1) are connected
			    // if the position of the max of S(j)
			    // belongs to S(j+1) 		  
			    			    
   void change_color_tree(int StructFather, float *TabColor, int NewColor);
                             // Change recursively, starting from StructFather, 
			     // the object index number of a tree.
			     // The new  object index is NewColor.
			     			    
   void create_graph(float *TabColor, MRNoiseModel &NM);
                            // Make the graph:
                            //  Each structure are grouped into a single set. 
			    //  Each set receive a color (index).
                            //  TabColor[i] = color (or index) of the structure i
                            //  if KillIsolatedObj = True then goup with a 
			    //     single structure are killed
                            //  if VisionModel = MVM_OVERLAP then kill all 
			    //     structures (and their fathers)  which have 
			    //     several sub-structures. We consider that
                            //     sub-structures merge to make a bigger one.

   void assign_struct_to_obj(int NumStruct, int NumObj, info_obj2d *TO);
                           // assign a structure to an object
			   // out: TO[NumObj].NbrStruct
			   //       ...
			   
   void assign_tree_to_obj(int RootTree, int NumObj, info_obj2d *TO);
                           // assign recursively a tree to an object

   void set_size_obj(int NumObj, info_obj2d *TO);
			  // calculate the size of an object
     		          // Image size must be larger than the structure sizes     
    		          // we add a border around the image. The size of the border depends
    		          // also on the last scale of the object
     		          // out: TO[NumObj].DepImaX, TO[NumObj].EndImaX
     		          //      TO[NumObj].DepImaY, TO[NumObj].EndImaY
    		          //      TO[NumObj].SizeImaX, TO[NumObj].SizeImaY
     
   
   void get_tree_struct_number(int RootTree, intarray &TS, int &IndS);
    		          // out: for the N structures of the tree starting 
			  //      at RootTree, 
    		          //      set TS(i) = index of the ith structure of 
			  //      the tree, i = 1..N  
			  //      IndS = N 
	  
   void count_tree(float * & TabColor, float * & TabNewColor, int & ncol); 
                            // Count the number of objects.
                            // Change the index object table in order to a have object index 
                            // between 1..NbrObject
                            // in:  TabColor = object index table
			    //         TabColor[StructNumber] = Object number
                            // out: ncol = total number of object
                            //      TabNewColor = new object index of the 
			    //                    ith structure

   void get_info_obj(float * & TabColor);
                            // 
			    // out: TabStruct[i].NumConnectObj
			    //           object number of the ith structure
			    //       TabObj = object table
			    //       TabSubObj  = sub-object table


   Bool is_sub_struct(int StructFather, int CurrentStruct, int StructSub);
                           // return True if the CurrentStruct is a 
			   // sub-structure
			   // A sub-structure is defined by
			       // 1) S(j) > S(j+1) 
			       // 2) S(j) > S(j-1) 
	                       // 3) DIST[PosMax(j+1) - PosMax(j)] > DistMin
	                       // 4) Flux[ S(j)] > PercentFlux*Flux[ S(j+1)]
   void find_sub_struct();
  /*
   1) segmente each band of the multiscale transform 
     if (VisionCleanBord == True) then structures on the border are killed.
     out: PosFirst[s] = first label of the first structure of band s
     out: TabNbrStruct[s] = number of structures in band s
     out: NbrStruct = total number of structures
     out: Segment.band(s) contains the segmented image with label 
          between PosFirst[s] and PosFirst[s]+TabNbrStruct[s]-1
   2) find information about each structure	  
     out:TabStruct[1..NbrStruct]
              TabStruc[k].Scale = scale of the structure
	      TabStruct[Ind].Npix = number of pixels of the structures
	      ...
  
  
  */ 
   int pos_ima(int x, int Scale)
   {
      return pos_in_ima(x, Scale, Segment.Set_Transform);
   }
   
   int pos_transform(int x, int Scale)
   {
      return pos_in_transform(x, Scale, Segment.Set_Transform);
   }
   

   int size_min_struct(int Band)
   {
      int SizeMin=2;
      if  ((Segment.Set_Transform == TRANSF_SEMIPYR) && (Band > 0)) SizeMin = 4; 
      else if  (Segment.Set_Transform == TRANSF_PAVE) SizeMin = 1 << (Band+1);
      return SizeMin;
   }
   
   int size_max_struct(int Band)
   {
      int Nls = Segment.size_band_nl(Band);
      int Ncs = Segment.size_band_nc(Band);
      int SizeMax = (int) (Nls*Ncs*0.5);
      return SizeMax;
   }
   
   Bool pyr_nextscale(int s)
   {
      Bool Val = False;
      if ((Segment.Set_Transform ==  TRANSF_PYR ) ||
	   ((Segment.Set_Transform == TRANSF_SEMIPYR) && (s > 0))) Val = True;
      return Val;
   }
   Bool pyr_upscale(int s)
   {
      Bool Val = False;
      if ((Segment.Set_Transform ==  TRANSF_PYR ) ||
	   ((Segment.Set_Transform == TRANSF_SEMIPYR) && (s > 1))) Val = True;
      return Val;
   }
   int pos_next(int X, int s)
   {
       int Xn = X;
       if (pyr_nextscale(s) == True) Xn /= 2;
       return Xn; 
   }
   int pos_up(int X, int s)
   {
       int Xn = X;
       if (pyr_upscale(s) == True) Xn *= 2;
       return Xn; 
   }
   // Bool obj_nb_struct_ok(int NumObj, int MinStructNumber=1);
                                 // return true if the object contains
				 // more than MinStructNumber structures
   // Bool obj_size_ok(int NumObj); // return false if the object contains
                                 // only one structure and if the 
				 // structure size (in pixel) > 2^(Scale+1)
				 //   with Scale = 0..NbrScale-1
   public:
      Bool MVMDebug;            // Debug information
      float DistMin;         // Minimum Distance between the maxima
                            // to define sub-object
      Bool KillSmallIsolatedObj; // Delete object which are small and isolated
      Bool KillIsolatedObj; // Delete object which are  isolated
      float PercentFlux;    // minimum flux percentage to define a 
                            // a subobject                      
      Bool VisionCleanBord; // Object touching the border are killed
      Bool Verbose;        // Print processing results
      t_mvm VisionModel;   // Type of Multiscale Vision Model
      int FirstScale;      // First scale to be used for the reconstruction
      int LastScale;       // Last scale to be used for the reconstruction
      Bool PosDetect;      // Reconstruct only positive structrure
      Bool NegDetect;      // Reconstruct only negative structrure
      int AddBorderX;      // Add a border around the reconstructed 
      int AddBorderY;      // object.
      int MinSizeObjX;     // Minimum size for the objects to be reconstructed
      int MinSizeObjY;     
      int FirstObjectScale; // Minimum scale for an object to be defined
                            // i.e. if the PSF is large, structure at the
                            // first scale cannot defined an object.
                            // It means that if FirstObjectScale > 0
                            // overlapped is not treated for scale lower that
                            // FirstObjectScale (def FirstObjectScale value is 1)
                            // (object at the first scale are not considered)
                            // This has the advantage to be more robust to the
                            // noise.
      MR2D_Seg();
      MR2D_Seg(MultiResol  & MR_Data);
      void init(MultiResol & MR_Data);
      void print_info_struct (int NumStruct);
      void print_info_obj(Bool InfoStruct=False);
      int nbr_obj() { return NbrTotObj; };
      int nbr_sub_obj() { return NbrSubObj; };
      int nbr_obj_in_scale(int s) { return  TabNbrObjPerScale[s]; };
      int nbr_struct() { return NbrStruct; };
            
      struct2d  *  struct_in_obj(int NumObj, int NumStructInObj);
      struct2d  *  struct_in_obj(int Scale, int NumObjInScale, int NumStructInObj);
      struct2d  *  get_struct_info(int NumStruct);
      info_obj2d * get_obj_info(int NumObj); 
                  // return information about the object NumObj
		  // NumObj = 1 .. NbrObj
		  
      info_obj2d * get_obj_info(int Scale, int NumObjInScale);
                  // return information about the object NumObjInScale
		  // of the scale Scale
		  // NumObjInScale = 1 .. TabNbrObjPerScale[Scale]
		  
      void segment(MultiResol & MR_Data, MRNoiseModel & NoiseModel);
                  // MR_Data: input = multiresolution transform of the data
                  // MRNoiseModel: input = noise model
                  // set the field of the structure
		  
		  /* 
		   segment routine create the object graph of a multiscale
		   transform by calling the following sub-routines:
		      1) call find_structures: segmente the multiscale 
		         transform, and store in TabStruct information
			 about each structure                    
		      2) call connect_struct
		         link structure and 
		      3) create_graph
		      4) count_tree
		      5) get_info_obj
		      
		  */
      void get_mrobj(int NumObj, MultiResol & MR_Data, MultiResol & MR_Obj);
                  // extract the Wavelet transform of the 
		  // object number  NumObj
		  // MR_Data: input = multiresolution transform of the data
		  // NumObj: input = object number to extract
		  // MR_Obj: output = wavelet coefficients of the object
      void get_mrobj(int Scale, int NumObjInScale, 
                     MultiResol & MR_Data, MultiResol & MR_Obj);
      	           // extract the Wavelet transform of the 
		   // object number NumObjInScale of the scale Scale
		  
      void make_ima(MultiResol & MR_Data, Ifloat &Im_rec);
      void make_graph(Bool KillSize, Bool KillIsol, int MinStr);
      void free();
      void clean(); // reset detection structures without deallocation
      ~MR2D_Seg();
};

#endif
