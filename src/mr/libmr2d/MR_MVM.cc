//----------------------------------------------------------
//	Copyright (C) 1999 CEA
//----------------------------------------------------------
//
//    UNIT
//
//    Version:  1.0
//
//    Author: 	J.L. Starck
//
//    Date: 	27/05/99
//    
//    File:  	MVM.cc
//
//----------------------------------------------------------

#include"MR_MVM.h"

#define DEBUG_SEG 0 
#define DEBUG_STRUCT 0
#define DEBUG_GET_INFO 0

/***********************************************************************/

void MR2D_Seg::clean()
{
   FirstScale = 0;
   LastScale = 0;
   PosDetect = False;
   NegDetect = False;
   Verbose = False;
   Nl = 0;
   Nc = 0;
   Nbr_Plan = 0;
   LastScale = 0;
   
   NbrStruct = 0;
   NbrSubStruct = 0;
   NbrObj = 0;
   NbrSubObj = 0;
   NbrTotObj = 0;
   VisionCleanBord = True;
   TabStruct = NULL;
   TabObj = NULL;
   DistMin=1;
   PercentFlux=0.1;
   for (int i=0; i < MAX_SCALE; i++) TabNbrObjPerScale[i]=0;
   KillSmallIsolatedObj = False;
   KillIsolatedObj = False;
   VisionModel=MVM_NO;
   AddBorderX = AddBorderY = 0;
   MinSizeObjX = MinSizeObjY  = 0;
   FirstObjectScale = 1;
   MVMDebug = False;
}

/***********************************************************************/

void MR2D_Seg::init(MultiResol  & MR_Data)
{
   clean();
   Nl = MR_Data.size_ima_nl();
   Nc = MR_Data.size_ima_nc();
   Nbr_Plan = MR_Data.nbr_band();
   LastScale = Nbr_Plan-1;
   Segment.alloc(Nl, Nc, MR_Data.nbr_scale(), MR_Data.Type_Transform, "Segmentation");
}

/***********************************************************************/

MR2D_Seg::MR2D_Seg()
{
   clean();
}

/***********************************************************************/

MR2D_Seg::MR2D_Seg(MultiResol & MR_Data)
{
   this->init(MR_Data);
}

/***********************************************************************/

void MR2D_Seg::free()
{
   if (Nbr_Plan > 0)
   {
      if (TabObj != NULL) delete [] TabObj;
      if (TabStruct != NULL) delete [] TabStruct;
      Segment.free();
      clean();
   }
}
 
 
/***********************************************************************/

MR2D_Seg::~MR2D_Seg()
{
   this->free();
}

/***********************************************************************/

void MR2D_Seg::find_structures(MultiResol & MR_Data, MRNoiseModel & NoiseModel)
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

{
  int i,j,k,s;
  int NbrScale = Nbr_Plan;
  // int Nstruc = 0;

  // Segmentation scale by scale 
  int *PosFirst = new int [NbrScale];
  for (s=0; s < NbrScale; s++) Segment.band(s).init();
  NbrStruct = 0;
  PosFirst[0]=1;
  
  
  for (s=0; s < NbrScale-1; s++)
  { 
     // segment band s
     //   PosFirst[s] = first label used par im_segment
     //                 structure index in band s will be in 
     //                 the interval 
     //                      [ PosFirst[s], PosFirst[s]+TabNbrStruct[s]-1 ]
     im_segment(MR_Data.band(s), Segment.band(s), TabNbrStruct[s], 
                (float) 0., VisionCleanBord, PosFirst[s]);
     if (Verbose == True)
           cout << " Scale " << s+1 << " Nb = " << TabNbrStruct[s] << endl;
  
     // If we wabt to sub-segment by the Rue-Bijaoui algorithm
     //Ifloat *IS;
     //IS = sous_segmente_iter (MR_Data.band(s), 10);
     //Segment.band(s) = *IS;
     //delete IS;
     PosFirst[s+1] = PosFirst[s] + TabNbrStruct[s];    
     NbrStruct += TabNbrStruct[s];
  }
  
  
   if (MVMDebug == True)
       cout << " step1 number of detected structures = " <<  NbrStruct << endl;
 
   // 1 dim object table allocation
   if (NbrStruct != 0) 
   {
      TabStruct = new struct2d [NbrStruct+1];
      struct2d *TS;
      int ind_obj = 1;
 
      // initialize the structures
      for (s=0; s < NbrScale-1; s++)
      {
#if DEBUG_SEG
cout << " Scale = " <<  s+1 << " Ns = " << TabNbrStruct[s] << endl;
#endif
          // initialize the structures
          //  TabStruct[i].Scale = scale of the structure
          //  TabStruct[i].NumStruct = Structure index
          //  TabStruct[i].NumConnectObj = corresponding color
	  //                                 in Segment.band(s)
	  TS = TabStruct + PosFirst[s];
          for (k=0; k < TabNbrStruct[s]; k++)
	  {
	     TS[k].Scale = s;
 	     TS[k].NumStruct = ind_obj;
	     TS[k].NumConnectObj = k+1;
             ind_obj++;
	  }
// #if DEBUG_SEG
// cout << " fill the structure array: FP = " <<  PosFirst[s]<<  endl;
// #endif      
          // fill the object table from the segmented data
	  if (TabNbrStruct[s] != 0)
	  {
            for (i=0; i < MR_Data.size_band_nl(s); i++)
	    for (j=0; j < MR_Data.size_band_nc(s); j++)
	    {
	      int Ind = (int) Segment(s,i,j);
	      if (Ind > NbrStruct)
	      {
	          cerr << "Error: bad index value = " << Ind << endl;
		  exit(-1);
	      }
              if (Ind > 0)
	      {
		 TabStruct[Ind].Npix ++;
		 if (TabStruct[Ind].Npix == 1)
		 {
		    TabStruct[Ind].Sizex = 1;
		    TabStruct[Ind].Sizey = 1;
		    TabStruct[Ind].ValMax = MR_Data(s,i,j);
		    TabStruct[Ind].Noise = NoiseModel.sigma(s,i,j);
		    TabStruct[Ind].Maxx = j;
		    TabStruct[Ind].Maxy = i;
		    TabStruct[Ind].Depx = j;
		    TabStruct[Ind].Depy = i;
		    TabStruct[Ind].Flux = MR_Data(s,i,j);
		}
		else
		{
		   TabStruct[Ind].Flux += MR_Data(s,i,j);
 		   if (MR_Data(s,i,j) > TabStruct[Ind].ValMax)
		   {
		        TabStruct[Ind].ValMax = MR_Data(s,i,j); 
			TabStruct[Ind].Noise = NoiseModel.sigma(s,i,j); 
			TabStruct[Ind].Maxx = j;
		        TabStruct[Ind].Maxy = i; 
		   }  
		   if (i < TabStruct[Ind].Depy)
		           TabStruct[Ind].Depy = i;
		   if (j < TabStruct[Ind].Depx)
		           TabStruct[Ind].Depx = j;
		   if (j-TabStruct[Ind].Depx + 1 > TabStruct[Ind].Sizex)
		      TabStruct[Ind].Sizex = j-TabStruct[Ind].Depx + 1;
		   if (i-TabStruct[Ind].Depy + 1 > TabStruct[Ind].Sizey)
		      TabStruct[Ind].Sizey = i-TabStruct[Ind].Depy + 1;
	       }
            }
	 }
        }
      }
  }
  delete [] PosFirst;	      	   
}


/***********************************************************************/

void MR2D_Seg::connect_struct()
/* 
   1) set the TabRel relation using maxima position
    If the structure at scale s has N sub-structures at scale s-1,
    then TabRel(i,k) = index of the kth substructure
  
   2) structures with a size larger than size_max_struct(s) are killed
   
   3) set few fields of TabStruct:
          TabStruct[i].Father      = index of the father structure (at s+1)
          TabStruct[i].NSubStruct  = number of sub-structures  (at s-1)
	  TabStruct[i].NbrTotalSub = number of sub-structures from s=0 to s-1 
*/
{
  int i,s;
  int NbrScale = Nbr_Plan;
  int posmaxx,posmaxy,indrel;
  TabRel.init();
  
  if (Verbose == True) cout << "connect_struct ..." << endl;
        
  // Kill the bad structure (too big or too small)
  for (i = 1; i <= NbrStruct; i++)
  { 
     int Scale = TabStruct[i].Scale;
     if (TabStruct[i].Npix > size_max_struct(Scale))  
     {
        TabStruct[i].BadStruct = True;
        if (MVMDebug == True)
        {
           cout << "Kill too big struct... " << i << " Scale = " << Scale << endl;
           cout << "    Npix = " << TabStruct[i].Npix << " size max = " << size_max_struct(Scale) << endl;	
        }
     }
  }
  
  NbrSubStruct=0;
  for (i = 1; i <= NbrStruct; i++)
  {
     if (TabStruct[i].BadStruct == False)
     {
       posmaxx = TabStruct[i].Maxx;
       posmaxy = TabStruct[i].Maxy;
       s = TabStruct[i].Scale;
       if (s < NbrScale-2)
       {
         // if the position maximum of the maximum of a structure
	 // belongs to a structure at the next scale, the structures
	 // are connected
         int NextPosMaxX=pos_next(posmaxx,s);
	 int NextPosMaxY=pos_next(posmaxy,s);
 	 int IndNext = (int) Segment(s+1, NextPosMaxY, NextPosMaxX);
 	 // IndNext contains either the structure number or 0
 	 if (IndNext > NbrStruct)
	 {
	    cerr << "Error: band NbrStruct = " << NbrStruct << endl;
	    exit(-1);
	 }
	 
	 // test if a father exists
         if ((IndNext > 0) && (TabStruct[IndNext].BadStruct == False))
         {
  	     TabStruct[i].Father = IndNext;
 	     TabRel(IndNext,0) ++;
	     indrel = TabRel(IndNext,0);
	     TabRel(IndNext,indrel) = i;
	     TabStruct[IndNext].NSubStruct ++;
	     TabStruct[IndNext].NbrTotalSub += TabStruct[i].NbrTotalSub+1;
         }
       }
    }
  }

  if (MVMDebug == True)
  {
     cout << "print interscale relation table " << endl;
     for (i = 1; i <= NbrStruct; i++)
     {
       int colrel = TabRel(i, 0);
       if (colrel != 0)
       {
         cout << "Connect struct " << i << " nc = " << colrel << " to : ";  
         for (s=1; s <= colrel; s++) cout << TabRel(i, s) <<  " ," ;
	 cout << endl;
       }
     }
  }
}  



/***********************************************************************/

void MR2D_Seg::change_color_tree(int StructFather, float *TabColor, int NewColor)
// Change recursively, starting from StructFather, 
// the object index number of a tree.
// The new  object index is NewColor
{
   int i,Ind;
   TabColor[StructFather] = NewColor;
   if (TabStruct[StructFather].NSubStruct > 0)
   {
      // Loop on each sub-structure
      for (i = 1; i <= TabStruct[StructFather].NSubStruct; i++)
      {
         Ind = TabRel(StructFather,i);
	 change_color_tree(Ind, TabColor, NewColor);
      }
   }
}

/***********************************************************************/

void MR2D_Seg::create_graph(float *TabColor, MRNoiseModel &NM)
// Make the graph:
//  Each structure are grouped into a single set. Each set receive a color
//  (index).
//  TabColor[i] = color (or index) of the structure i
//   if KillIsolatedObj = True then goup with a single structure are killed
//   if VisionModel = MVM_OVERLAP then kill all structures (and their fathers)
//                    which have several sub-structures. We consider that
//                    sub-structures merge to make a bigger one.
//        
{
  int i;
  
  if (Verbose == True)
       cout << "create the object graph ..." << endl;
 
  TabColor[0]=0;     
  for (i = 1; i <= NbrStruct; i++) 
      TabColor[i] = (TabStruct[i].BadStruct == False) ? i: -1;
 
  NbrSubStruct=0;
  for (i = 1; i <= NbrStruct; i++)
  {
    if (TabStruct[i].BadStruct == False)
    {
        int IndNext = TabStruct[i].Father;
	
        // Root Structure ==> IndNext == 0
	// If it is an isolated object (no father and no sub-struct)
	if ((IndNext == 0) && (TabStruct[i].NSubStruct == 0))
	{
 	   // test if it is an isolated object
           // At the first scale, isolated object are always killed
           // because of the noise which dominates.
           // At other scales, we kill isolated object only if it is specified
  	   if ((KillIsolatedObj == True) ||
               (TabStruct[i].Scale == 0) ||
 	    ((KillSmallIsolatedObj == True) && 
	           (TabStruct[i].Npix < size_min_struct(TabStruct[i].Scale))))
           {
 	       TabStruct[i].BadStruct = True;
	       TabColor[i] = -1;
              if (MVMDebug == True)
                cout <<  "KILL ISO OBJ: Struct Number = " << i << "  Npix = " << TabStruct[i].Npix << " size min = " << size_min_struct(TabStruct[i].Scale) << endl;
 	   }
 	}
 	else
	{
     	   switch(VisionModel)
 	   {
	     case MVM_NO:
	        change_color_tree(i, TabColor, i);
	        if ((IndNext > 0) &&  (TabStruct[IndNext].NSubStruct > 1))
                {
	           TabStruct[i].SubStruct = True;
		   NbrSubObj ++;
                }
	        break;
	     case MVM_OVERLAP:
                {
                int ss, IndSub, NbrSignifSub = 0;
                int SNpix = TabStruct[i].Sizex*TabStruct[i].Sizey;
 	        float FluxSub=0;
	        int si, SubNpix, MNpix=0;

                // Calculate first some informations to be used in order
                // to determine if we are in case of overlapping objects
	        for (ss = 1; ss <= TabStruct[i].NSubStruct; ss++)
	        {
                   // Count the number of substructures which are significant
                   // at a level 5sigma (substructures at a lower level
                   // are not considered for overlapping. They are included
                   // in the father structure
	           IndSub = TabRel(i,ss);
                   int ScaleSub = TabStruct[IndSub].Scale;
                   float MaxSub = TabStruct[IndSub].ValMax;
                   int PosMaxXSub = TabStruct[IndSub].Maxx;
                   int PosMaxYSub = TabStruct[IndSub].Maxy;
                   float NoiseLevel= NM.sigma(ScaleSub,PosMaxYSub,PosMaxXSub);
                   if (MaxSub > 5.*NoiseLevel) NbrSignifSub++;

 	           // Calculate the integrate flux of all sub-structures
                   // if the integrated flux in the substructed is 
                   // larger than half the flux of the structure,
                   // then most of the flux exists already at a finer 
                   // scale, and we consider it is an overlap.
                   // If the structure size is larger than 16 times the
                   // largest substructure, we also consider it is an
                   // overlap.
                   SubNpix=TabStruct[IndSub].Sizex*TabStruct[IndSub].Sizey;
 	           float Norm = 1.;
	           if (pyr_nextscale(TabStruct[IndSub].Scale) == True) Norm = 4.;
	           FluxSub += TabStruct[IndSub].Flux*Norm;
                   if (MNpix < SubNpix) MNpix = SubNpix;
                }

	        if ((NbrSignifSub  > 1) 
                    && (TabStruct[i].Scale > FirstObjectScale)
                    && ((FluxSub > 0.5*TabStruct[i].Flux) || 
                                                     (SNpix > 16 * MNpix)))
                {
 	               if (MVMDebug == True) 
                        cout << " FluxSub = " << FluxSub  << " F = " << TabStruct[IndSub].Flux << endl;
		       ss=TabStruct[i].Father;
		       si=i;
		       TabStruct[i].BadStruct = True;
	               TabColor[i] = -1;
                       if (MVMDebug == True)
                           cout <<  "KILL OVER OBJ: Struct Number = " << i << " nbr sub = " << TabStruct[i].NSubStruct << endl;
 	   
	              while ((ss > 0) && 
	                  (TabStruct[ss].ValMax < TabStruct[si].ValMax))
	              {
                       if (MVMDebug == True)
                         cout <<  "KILL OVER FATHER OBJ: Struct Number = " << ss << endl;
 
	               TabStruct[ss].BadStruct = True;
		       TabColor[ss] = -1;
		       si = ss;
		       ss = TabStruct[si].Father;
	              }
                      if (MVMDebug == True) cout <<  "END KILL OVERLAP" << endl;
		}
		else change_color_tree(i, TabColor, i);
                }
	        break;
	     default: change_color_tree(i, TabColor, i);
	   }
        }
      }
   }
}  

/*
	    else
	    {
	       cout << "NOSUB Structure " << i <<  " FluxSub = " << FluxSub << " Flux = " << TabStruct[i].Flux  << endl;
	    }

       }*/

	       /*
	 if ((IndNext > 0) &&  (TabStruct[IndNext].NSubStruct > 1))
         {
	    // search if it is a substructure
 	    if (StrucInf.NSubStruct == 0)
	        StrucInf.SubStruct = is_sub_struct(IndNext,i,0);
	    else 
	    {
	        ss=1;
 		while ((ss <= StrucInf.NSubStruct) && (StrucInf.SubStruct == False))
		{
		   IndSub = TabRel(i,ss);
		   StrucInf.SubStruct = is_sub_struct(IndNext,i,IndSub);
		   ss++;
		}
  	     }
	     if (StrucInf.SubStruct == True)
	     {
	         NbrSubStruct++;
		 cout << " TT = " << StrucInf.BadStruct << endl;
		 cout << " TTD = " << TabStruct[i].BadStruct << endl;
	     }
	     
	     // Gives the color of the father to the structures
	     // and its sub-structrure
 	     if (StrucInf.SubStruct == False) 
	                         change_color_tree(i, TabColor, IndNext);
         } */
/***********************************************************************/


void MR2D_Seg::count_tree(float * & TabColor, float * & TabNewColor, int & ncol)
// Count the number of objects.
// Change the index object table in order to a have object index 
// between 1..NbrObject
// in:  TabColor = object index table
//                  TabColor[i] = object index of the ith structure
// out: ncol = total number of object
//      TabNewColor = new object index of the ith structure
{
  // count the number of objects
  int i,j;
  unsigned long *TabIndex = new unsigned long [NbrStruct+1];
  float *TabOrderColor = new float [NbrStruct+1];
  for (i = 0; i <= NbrStruct; i++) 
  {
     TabNewColor[i] = i;
     TabOrderColor[i] = TabColor[i];
  }
  unsigned long No = NbrStruct;
  indexx(No, TabColor, TabIndex);
  sort(No, TabOrderColor);
#if DEBUG_SEG
cout << "delete the hole in TabColor  ..." << endl;
  for (i = 1; i <= NbrStruct; i++) 
  {
     cout << TabOrderColor[i] << " ";
  }
  cout << endl;
#endif
  i=1;
  ncol = 1;
  int Col_ind;
  int Col_obj;
  do
  {
      // delete the "hole" in TabColor by creation of a new TabColor of
      // name TabNewColor
      Col_ind = TabIndex[i];
      if ((Col_ind <= 0) || (Col_ind > NbrStruct))
      {
         cerr << "Error: PB index = " << Col_ind << endl;
	 exit(-1);
      }
      Col_obj = (int) TabColor[ Col_ind ];
      if (Col_obj > 0)
      {
          for (j = 1; j <= NbrStruct; j++) 
                      if (TabColor[j] == Col_obj) TabNewColor[j] = ncol;

         // go to the next color
         Bool Stoploop=False;
         do
         {
           i++;
           if (i > NbrStruct) Stoploop = True;
           else
           {
             if (TabOrderColor[i] != TabOrderColor[i-1])
             {
                Stoploop = True;
                ncol++;
             }
             else Stoploop = False;
           }
        } while ((i <= NbrStruct) && (Stoploop == False));
     }
     else i++;
  } while (i <= NbrStruct);
 
 delete [] TabIndex; 
 delete [] TabOrderColor;
}

/***********************************************************************/

void MR2D_Seg::assign_struct_to_obj(int NumStruct, int NumObj, info_obj2d *TO)
// assign a structure to an object
// out: TO[NumObj].NbrStruct
{
    int i = NumStruct;
    int Scale = TabStruct[i].Scale;

//  if (NumStruct == 46)
//      cout << "  46: " << NumObj << " " << TabStruct[i].ValMax << endl;
    
    // Add a new sub-structure to the object
    TO[NumObj].NbrStruct += 1;
    
    // Update the RootScale field
    if ((NbrStruct == 1) || (Scale >= TO[NumObj].RootScale))
    {
	TO[NumObj].RootIndex = i;
	TO[NumObj].RootScale = Scale;
    }
    
    // if the structure maximum is larger than object maximum
    // then store information about the maximum
    if (ABS(TabStruct[i].ValMax) > ABS(TO[NumObj].Max ))
    {
        TO[NumObj].Max =  TabStruct[i].ValMax;
        TO[NumObj].PosMaxImaX = pos_ima(TabStruct[i].Maxx, Scale);
        TO[NumObj].PosMaxImaY = pos_ima(TabStruct[i].Maxy, Scale);
        TO[NumObj].ScaleMax = Scale;
        TabStruct[i].NumScale = Scale;
    }
    
    // Update DepImaX, DepImaY, EndImaX, EndImaY
    int Depx = pos_ima(TabStruct[i].Depx, Scale);
    int Depy = pos_ima(TabStruct[i].Depy, Scale);
    int Sx = pos_ima(TabStruct[i].Sizex, Scale);
    int Sy = pos_ima(TabStruct[i].Sizey, Scale);
    int Lastx = Depx + Sx - 1;
    int Lasty = Depy + Sy - 1;
    if (TO[NumObj].NbrStruct == 1)
    {
       TO[NumObj].DepImaX = Depx;
       TO[NumObj].DepImaY = Depy;
       TO[NumObj].EndImaX = Lastx;
       TO[NumObj].EndImaY = Lasty;
    }
    else
    {
       if (TO[NumObj].DepImaX > Depx) TO[NumObj].DepImaX = Depx;
       if (TO[NumObj].DepImaY > Depy) TO[NumObj].DepImaY = Depy;
       if (TO[NumObj].EndImaX < Lastx) TO[NumObj].EndImaX=  Lastx;
       if (TO[NumObj].EndImaY < Lasty) TO[NumObj].EndImaY=  Lasty;
   }
}

/***********************************************************************/

void MR2D_Seg::set_size_obj(int NumObj, info_obj2d *TO)
{
     // Image size must be larger than the structure sizes     
     // we add a border around the image. The size of the border depends
     // also on the last scale of the object
     // out: TO[NumObj].DepImaX, TO[NumObj].EndImaX
     //      TO[NumObj].DepImaY, TO[NumObj].EndImaY
     //      TO[NumObj].SizeImaX, TO[NumObj].SizeImaY
     
     int BordSize = 2+POW2(TO[NumObj].RootScale);
     int BordSizeX = BordSize + AddBorderX;
     int BordSizeY = BordSize + AddBorderY;
     int FullBordX = (BordSize+BordSizeX)*2;
     int FullBordY = (BordSize+BordSizeY)*2;
     int SOx,SOy;

     // We have to put a border around the object in order to have
     // a good reconstruction quality.
     // The border is dependant on the scale where the object appears
     // We fix:   BordSize = 2 + POW2(TO[NumObj].RootScale)
     // AddBorderX,AddBorderY are user parameter which allows us
     // to increase the border size.
     // MinSizeObjX,MinSizeObjY specify a minimum size for the
     // object to be reconstructed (user parameter)
     SOx = TO[NumObj].EndImaX - TO[NumObj].DepImaX + 1;
     TO[NumObj].SizeImaX = MAX(MinSizeObjX, SOx+FullBordX);
     FullBordX = TO[NumObj].SizeImaX - SOx;
     TO[NumObj].DepImaX -= FullBordX/2;
     TO[NumObj].EndImaX = TO[NumObj].DepImaX+TO[NumObj].SizeImaX-1;

     SOy = TO[NumObj].EndImaY - TO[NumObj].DepImaY + 1;
     TO[NumObj].SizeImaY = MAX(MinSizeObjY, SOy+FullBordY);
     FullBordY = TO[NumObj].SizeImaY - SOy;
     TO[NumObj].DepImaY -= FullBordY/2;
     TO[NumObj].EndImaY = TO[NumObj].DepImaY+TO[NumObj].SizeImaY-1;

//      int Dep;
//      Dep = TO[NumObj].DepImaX-BordSizeX;
//      if (Dep >= 0) TO[NumObj].DepImaX = Dep;
//      else
//      {
//         TO[NumObj].DepImaX = 0;
// 	TO[NumObj].EndImaX -= Dep;
//      }
// 
//      Dep = TO[NumObj].DepImaY-BordSizeY;
//      if (Dep >= 0) TO[NumObj].DepImaY = Dep;
//      else
//      {
//         TO[NumObj].DepImaY = 0;
// 	TO[NumObj].EndImaY -= Dep;
//      }
//      Dep = TO[NumObj].EndImaX+BordSizeX;
//      if (Dep < Nc) TO[NumObj].EndImaX = Dep;
//      else
//      {
//         TO[NumObj].EndImaX = Nc-1;
// 	TO[NumObj].DepImaX -= Dep - Nc + 1;
// 	if (TO[NumObj].DepImaX < 0) TO[NumObj].DepImaX = 0;
//      }
//      Dep = TO[NumObj].EndImaY+BordSizeY;
//      if (Dep < Nl) TO[NumObj].EndImaY = Dep;
//      else
//      {
//         TO[NumObj].EndImaY = Nl-1;
// 	TO[NumObj].DepImaY -= Dep - Nl + 1;
// 	if (TO[NumObj].DepImaY < 0) TO[NumObj].DepImaY = 0;
//      }
//
//     TO[NumObj].SizeImaX = TO[NumObj].EndImaX - TO[NumObj].DepImaX + 1;
//    TO[NumObj].SizeImaY = TO[NumObj].EndImaY - TO[NumObj].DepImaY + 1;
}

/***********************************************************************/

void MR2D_Seg::get_tree_struct_number(int RootTree, intarray &TS, int &IndS)
// out: for the N structures of the tree starting at RootTree, 
//      set TS(i) = index of the ith structure of the tree, i = 1..N
{
   int s,NbrRel;
   
   IndS++;
   TS(IndS) = RootTree;
   NbrRel = TabRel(RootTree,0);
   for (s = 1; s <= NbrRel; s++)
   {
      int indrel = TabRel(RootTree,s);
      if (indrel > 0) get_tree_struct_number(indrel, TS, IndS);
   }
}

/***********************************************************************/

void MR2D_Seg::assign_tree_to_obj(int RootTree, int NumObj, info_obj2d *TO)
// Assign recursively a tree to an object
{
   int s,NbrRel;
   
   assign_struct_to_obj(RootTree,NumObj,TO);
   NbrRel = TabRel(RootTree,0);
   for (s = 1; s <= NbrRel; s++)
   {
      int indrel = TabRel(RootTree,s);
      if (indrel > 0) assign_tree_to_obj(indrel, NumObj, TO);
   }
}

/***********************************************************************/

void MR2D_Seg::get_info_obj(float * & TabColor)
{
  int i,NumObj,NumSubObj=0;
   
  // search for each object its maximum value, its position, and its scale
  // count the number of structures per object
  // search the smallest image size which contain all the structures
  NbrTotObj = NbrObj+NbrSubObj;
  TabObj = new info_obj2d[NbrObj+NbrSubObj+1];
  TabScaleObj.alloc(Nbr_Plan-1,NbrTotObj+1);
  TabSubObj = TabObj+NbrObj;
  
  for (i = 1; i <= NbrStruct; i++)
  {
     if (TabStruct[i].BadStruct == False)
     {
         NumObj = (int) TabColor[i];
#if DEBUG_GET_INFO
cout << "Analyse struct " << i << " obj = " << NumObj << endl;
#endif         
         if ((NumObj < 1) || (NumObj > NbrObj))
         {
           cerr << "Error: bad object number = " << NumObj << endl;
	   exit(-1);
         }         
	 assign_struct_to_obj(i,NumObj,TabObj);
         TabStruct[i].NumConnectObj = NumObj;
         if (TabStruct[i].SubStruct == True)
	 {
	     // if there is no multiscale vision, sub-structure are
	     // automatically assign to the object of the father structure
	     if (VisionModel == MVM_NO)
	     {
 		TabObj[NumObj].NbrSubObj ++;
		NumSubObj ++;
		assign_tree_to_obj(i, NumSubObj, TabSubObj);
		TabSubObj[NumSubObj].NumFatherObj = NumObj;
		TabSubObj[NumSubObj].SubObj = True;
             }
	 }
      }
  }
  
#if DEBUG_GET_INFO
cout << "set_obj_info ... NbrObj = " << NbrObj << " SubObj = " << NbrSubObj << endl;
#endif 
  intarray TS(NbrTotObj+1);

  for (NumObj = 1; NumObj <= NbrTotObj; NumObj++)
  {
     int Ns = TabObj[NumObj].NbrStruct;
     int ScaleObj = TabObj[NumObj].ScaleMax;
     int NumObjInScale = TabNbrObjPerScale[ScaleObj]+1;
#if DEBUG_GET_INFO
cout << " Obj " << NumObj << "Nb struc = " << Ns << "Scale max = " <<  ScaleObj << endl;
cout <<  "   NbSub = " << TabObj[NumObj].NbrSubObj << endl;
cout << "    TabNbrObjPerScale[ScaleObj] = " << TabNbrObjPerScale[ScaleObj] << endl;
#endif 
     // allocate the table for structure index
     // PB a bug appears when allocating with Ns+1, with Ns = 1
     // that the reason why the fix a minimum a 32 fot the allocation
     // NOT RESOLVED
     TabObj[NumObj].TabIndStruct.alloc(MAX(32,Ns+1));
     TabObj[NumObj].TabIndSubObj.alloc(MAX(32,TabObj[NumObj].NbrSubObj+1));
     
#if DEBUG_GET_INFO
cout << " allocOK " << endl;
#endif 
    if (TabObj[NumObj].SubObj == True)
    {
        int FatherObj = TabObj[NumObj].NumFatherObj;
        TS(FatherObj) ++;
        TabObj[FatherObj].TabIndSubObj(TS(FatherObj)) = NumObj;
    }
          
     // fill TabNbrObjPerScale table
     // TabNbrObjPerScale[ScaleObj] is the current object number
     // at the scale ScaleObj
     TabObj[NumObj].ObjNumberInScale = NumObjInScale;
     TabNbrObjPerScale[ScaleObj] = NumObjInScale;
     
     // set the index table of the object per scale
#if DEBUG_GET_INFO
    cout << "TabScaleObj: " << TabScaleObj.nx() << " " << TabScaleObj.ny() << endl;
    cout << "      ScaleObj " << ScaleObj  << "  NumObjInScale " <<    NumObjInScale << endl;  
    if (ScaleObj >= TabScaleObj.nx()) exit(-1);
    if (NumObjInScale >= TabScaleObj.ny()) exit(-1);
#endif 

     TabScaleObj(ScaleObj, NumObjInScale) = NumObj;

    // set the size of the object
    set_size_obj(NumObj,TabObj);
    
    int Si=0;
    int Ri = TabObj[NumObj].RootIndex;
    // search the structure index of the object. Put them in TabIndStruct
    get_tree_struct_number(Ri, TabObj[NumObj].TabIndStruct, Si);
  } 
  
#if DEBUG_GET_INFO
cout << " SubAna " << endl;
#endif 
//   for (NumSubObj = 1; NumSubObj <= NbrSubObj; NumSubObj++)
//   {
//       NumObj = TabSubObj[NumSubObj].NumFatherObj;
//       TS(NumObj) ++;
//       TabObj[NumObj].TabIndSubObj(TS(NumObj)) = NumSubObj;
//       set_size_obj(NumSubObj,TabSubObj);
//       int Si=0;
//       int Ri = TabSubObj[NumSubObj].RootIndex;
//       TabSubObj[NumSubObj].TabIndStruct.alloc(MAX(32,TabSubObj[NumSubObj].NbrStruct));
//       get_tree_struct_number(Ri, TabSubObj[NumSubObj].TabIndStruct, Si);
//   }
  
#if DEBUG_GET_INFO
cout << "set TabObjIndStruct ..." << endl;
#endif

  // put in TabObj[NumObj].TabIndStruct all connected structures index 
  // Tabind[o] = index in the table TabIndStruct of the next 
  //             structure index to insert for the object o
//   int  *Tabind = new int [NbrObj+1];
//   for (i = 0; i <= NbrObj; i ++) Tabind[i] = 1;
//      
//   for (i = 1; i <= NbrStruct; i ++)
//   if (TabStruct[i].BadStruct == False)
//   {
//        
//      NumObj = TabStruct[i].NumConnectObj;
//      if (NumObj > NbrObj)
//      {
//         cout << "PB set_obj_info: NbrObj = " << NbrObj << " NumObj = " << NumObj << endl;
// 	exit(-1);
//      }
// #if DEBUG_GET_INFO
// cout << "NumObj=" << NumObj << " ind = " << Tabind[NumObj] << " ns = " << TabObj[NumObj].NbrStruct << endl;
// #endif     
//      TabObj[NumObj].TabIndStruct(Tabind[NumObj])=i;
//      Tabind[NumObj]++;
//   }
#if DEBUG_GET_INFO
cout << "END set_obj_info..." << endl;
#endif

   // delete Tabind;
}

/***********************************************************************/

void MR2D_Seg::segment(MultiResol & MR_Data, MRNoiseModel & NoiseModel)
{
  float *TabColor;
  
  if (Nl != MR_Data.size_ima_nl())
  {
   cerr << "Error: MR2D_Seg::segment  incompatible Nl in MR_Data ... " << endl;
  }
  if (Nc != MR_Data.size_ima_nc())
  {
   cerr << "Error: MR2D_Seg::segment  incompatible Nc in MR_Data ... " << endl;
  }
  if (Nbr_Plan != MR_Data.nbr_scale())
  {
    cerr << "Error: MR2D_Seg::segment  incompatible Nbr_Plan in MR_Data ... " << endl;
  }
  
  if (TabStruct != NULL) delete [] TabStruct;
  if (TabObj  != NULL) delete [] TabObj;

  // Find the structures
  find_structures(MR_Data, NoiseModel);

  if (MVMDebug == True) Segment.write("xx_seg.mr");

  // interscale relation taking into account subobject relations
  // TabColor[i] = object number of the structure i
  // Ncol = number of different objects
   if (NbrStruct != 0) 
   {
      // Set the connection table TabRel
      TabRel.alloc(NbrStruct+1,NbrStruct+1);
      connect_struct ();
      
      // for (int i=1; i <= NbrStruct; i++) print_info_struct(i);
      
      float *TabColorStruct = new float [NbrStruct+1];
      create_graph(TabColorStruct, NoiseModel);
      
#if DEBUG_SEG
 cout << "search sub-structure ..." << " NbrStruct = " << NbrStruct << endl;
#endif
      // find_sub_struct();

      TabColor = new float [NbrStruct+1];
      NbrObj=0;
      // for (int i = 0; i <= NbrStruct; i++) TabColor[i] = i;
      // NbrObj=NbrStruct;
      // Set the connection
      // TabColor is filled in count_tree
      // NbrObj contains the number of trees
      // TabColor[i] = color of the structure number. i = 1..Nstruc
      count_tree (TabColorStruct, TabColor, NbrObj);
      delete [] TabColorStruct;

      if (Verbose == True) cout << "Number of object = " << NbrObj << endl;

      get_info_obj(TabColor);
 
#if DEBUG_SEG
     print_info_obj(True);
     cout << "END SEGMENTATION..." << endl;
#endif

   if (MVMDebug == True)
   {
     for (int b=0; b < MR_Data.nbr_band()-1; b++)
     {
      for (int i=0; i < MR_Data.size_band_nl(b); i++)
      for (int j=0; j < MR_Data.size_band_nc(b); j++)
      {
	  int Ind = (int) Segment(b,i,j);
          if (Ind != 0)
          {
             if (TabStruct[Ind].BadStruct == True) Ind = -1;
             else
             {
                Ind = TabStruct[Ind].NumConnectObj;
             }
          }
          // Segment(b,i,j) = Ind;
          Segment(MR_Data.nbr_band()-1,i,j) += Ind;
      }
     }
     Segment.write("xx_obj.mr");
   }
   delete [] TabColor;
  }
}

/***********************************************************************/

void MR2D_Seg::make_ima(MultiResol & MR_Data, Ifloat &Im_rec)
{
    Ifloat Ima;
    //double Err;
    MultiResol MR_Obj;
    for (int o=1; o <= NbrObj; o++)
    {
       // cout << "Rec Obj " << o << endl;
       get_mrobj(o, MR_Data, MR_Obj);
       // cout << "resize " <<  TabObj[o].SizeImaY << " " << TabObj[o].SizeImaX << endl;
       Ima.resize(TabObj[o].SizeImaY, TabObj[o].SizeImaX);
       MR_Obj.recons(Ima);
       // rec_iter_grad_conj (MR_Obj, Ima, 10, Err, 1e-5);
       // INFO(Ima,"obj");
       // cout << "add ima" << TabObj[o].DepImaX << " " << TabObj[o].DepImaY << endl;
       add_image (Im_rec, Ima, TabObj[o].DepImaX, TabObj[o].DepImaY);
       MR_Obj.free();
    }
}

/***********************************************************************/

void MR2D_Seg::print_info_struct(int i)
{
  cout << "   Structure number " << TabStruct[i].NumStruct;
  if (TabStruct[i].SubStruct) cout << "  Sub-Struct: True " << endl;
  else 
  cout << "       Sub-Struct: False " << endl;
  cout << "       Scale " << TabStruct[i].Scale;
  cout << " NumScale " << TabStruct[i].NumScale;
  cout << " Depx " << TabStruct[i].Depx << endl;
  cout << "       Maxx " << TabStruct[i].Maxx;
  cout << " Maxy " << TabStruct[i].Maxy << endl;
  cout << " ValMax " << TabStruct[i].ValMax;
  cout << " Sizex " << TabStruct[i].Sizex; 
  cout << " Flux " << TabStruct[i].Flux; 
  cout << " NSub " << TabStruct[i].NSubStruct;
  cout << " Noise " << TabStruct[i].Noise << endl << endl;
}


/***********************************************************************/

struct2d  * MR2D_Seg::struct_in_obj(int NumObj, int NumStruct)
{
   info_obj2d *InfObj;
   if ((NumObj < 1) || (NumObj > NbrTotObj))
   {
      cerr << "Bad object number ... " << endl;
      exit(-1);
   }
   InfObj = get_obj_info(NumObj);
   if ((NumStruct < 1) || (NumStruct > InfObj->NbrStruct))
   {
      cerr << "Error: Bad structure number ... " << endl;
      exit(-1);
   }
   int Ind = (InfObj->TabIndStruct)(NumStruct);
   return (TabStruct+Ind);
}

/***********************************************************************/

struct2d  *  MR2D_Seg::struct_in_obj(int Scale, int NumObjInScale, int NumStruct)
{
   if ((Scale < 0) || (Scale >= Nbr_Plan-1))
   {
      cerr << "Error: Bad scale number in struct_in_obj ... " << endl;
      exit(-1);
   }
   if ((NumObjInScale < 1) || (NumObjInScale >  nbr_obj_in_scale(Scale)))
   {
      cerr << "Error: Bad scale object number in struct_in_obj ... " << endl;
      exit(-1);
   }
   int NumObj = TabScaleObj(Scale, NumObjInScale); 
   if ((NumStruct < 1) || (NumStruct > TabObj[NumObj].NbrStruct))
   {
      cerr << "Error: Bad structure number ... " << endl;
      exit(-1);
   }
   int Ind = TabObj[NumObj].TabIndStruct(NumStruct);
   return (TabStruct+Ind);
}
/***********************************************************************/

info_obj2d * MR2D_Seg::get_obj_info(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrTotObj))
   {
      cerr << "Error: Bad object number ... " << endl;
      exit(-1);
   }
      return (TabObj+NumObj);
}

/***********************************************************************/

struct2d  * MR2D_Seg::get_struct_info(int NumStruct)
{
   if ((NumStruct < 1) || (NumStruct > nbr_struct() ))
   {
      cerr << "Error: Bad structure number ... " << NumStruct << endl;
      exit(-1);
   }
   return (TabStruct + NumStruct);
}

/***********************************************************************/

info_obj2d * MR2D_Seg::get_obj_info(int Scale, int NumObjInScale)
{
   if ((Scale < 0) || (Scale >= Nbr_Plan-1))
   {
      cerr << "Error: Bad scale number in get_obj_info ... " << endl;
      exit(-1);
   }
   if ((NumObjInScale < 1) || (NumObjInScale >  nbr_obj_in_scale(Scale)))
   {
      cerr << "Error: Bad scale object number in get_obj_info ... " << endl;
      exit(-1);
   }
   int Ind = TabScaleObj(Scale, NumObjInScale);
   
   if ((Ind < 1) || (Ind > NbrTotObj))
   {
      cerr << "Error: Bad index number in get_obj_info ... " << endl;
      cerr << "        Scale = " << Scale  << endl;
      cerr << "        NumObjInScale = " <<  NumObjInScale << endl;
      cerr << "         Ind = " <<   Ind << endl;
      exit(-1);
   }
   return (TabObj + Ind);
}

/***********************************************************************/

void MR2D_Seg::print_info_obj(Bool InfoStruct)
{  
  int i,s,sc;
  info_obj2d *InfObj;
  cout << "Number of objects = " << NbrObj << endl;
  cout << "number of structures = " <<  NbrStruct << endl;
  for (s=0; s < Nbr_Plan-1; s++) 
      cout << "  Scale " << s+1 << " Nstruct = " << TabNbrStruct[s] << " Nobj = " << nbr_obj_in_scale(s) << endl;
     
  cout << endl;
  for (sc = Nbr_Plan-2; sc >= 0; sc --) 
  {
     int No = nbr_obj_in_scale(sc);
     for (s = 1; s <= No; s ++) 
     {
         InfObj = get_obj_info(sc, s);
	 cout << "Object number: " << TabScaleObj(sc, s) << endl;
         cout << "Object number in scale: " << sc+1 << " " << s << endl;
	 if (InfObj->SubObj == True)
	 cout << "  Sub-Obj: True " << endl;
         else cout << "  Sub-Obj: False " << endl;
         cout << "  Nbr structure: " << InfObj->NbrStruct   << endl;
         cout << "  From x " << InfObj->DepImaX << " to " << InfObj->EndImaX;
         cout << ", from y " << InfObj->DepImaY << " to " << InfObj->EndImaY;
         cout << ", Scale = " << InfObj->ScaleMax+1 << endl;
         cout << "  ValMax = " << InfObj->Max  << endl;
         cout << "  Nl = " << InfObj->SizeImaY << " Nc = " <<  InfObj->SizeImaX << endl;

         if (InfoStruct == True) 
              for (i = 1; i <= InfObj->NbrStruct; i++) 
  		print_info_struct(InfObj->TabIndStruct(i));
          cout << endl<< endl;
     }
  }
}

/***********************************************************************/

void MR2D_Seg::get_mrobj(int Scale, int NumObjInScale, 
                     MultiResol & MR_Data, MultiResol & MR_Obj)
{
   if ((Scale < 0) || (Scale >= Nbr_Plan-1))
   {
      cerr << "Bad scale number in get_mrobj ... " << endl;
      exit(-1);
   }
   if ((NumObjInScale < 1) || (NumObjInScale >  nbr_obj_in_scale(Scale)))
   {
      cerr << "Bad scale object number in get_mrobj ... " << endl;
      exit(-1);
   }
   int Ind = TabScaleObj(Scale, NumObjInScale);
   get_mrobj(Ind,  MR_Data, MR_Obj);
}

/***********************************************************************/

void MR2D_Seg::get_mrobj(int NumObj, MultiResol & MR_Data, MultiResol & MR_Obj)
{
   struct2d *StrucInf;
   info_obj2d *ObjInf;
   int Nco, Nlo, i,j,ii,jj,st,Np;
   
   if (NbrObj < 1)
   {
      cout << "Error: no detected object ... " << endl;
      exit(-1);
   }
   if ((NumObj < 1) || (NumObj > NbrTotObj))
      {
      cout << "Error: bad object number ... " << NumObj << endl;
      exit(-1);
   }
   
   ObjInf = get_obj_info(NumObj);
   Nco = ObjInf->SizeImaX;
   Nlo = ObjInf->SizeImaY;
   Np = ObjInf->RootScale+2;   
#if DEBUG_SEG
 cout << "Get_mrobj : Obj = " << NumObj << "  Scale = " << Np << " Nlo = " << Nlo << " Nco = " << Nco;
 cout << " NbrStruct = " << ObjInf->NbrStruct << endl;
#endif
   MR_Obj.alloc(Nlo, Nco, MAX(3,Np), Segment.Type_Transform, "MR_Obj");
 
   for (st=1; st <= ObjInf->NbrStruct; st++)
   {
       StrucInf = struct_in_obj(NumObj,st);
       int Ex=StrucInf->Depx+StrucInf->Sizex;
       int Ey=StrucInf->Depy+StrucInf->Sizey;
       int s = StrucInf->Scale;
       int DepXs = pos_transform(ObjInf->DepImaX, s);
       int DepYs = pos_transform(ObjInf->DepImaY, s);
       int IndS = ObjInf->TabIndStruct(st);
       
#if DEBUG_SEG
 cout << "  Struc = " << IndS << "  Scale = " << s+1 << " Depy = " << StrucInf->Depy << " Depx = " << StrucInf->Depx;
 cout << "    Ex = " << Ex << " Ey = " << Ey << endl;
#endif       
       for (i=StrucInf->Depy; i < Ey; i++)
       for (j=StrucInf->Depx; j < Ex; j++)
       {
          if (Segment(s,i,j) == IndS)
	  { 
	     ii = i - DepYs;
	     jj = j - DepXs;
	     if ( (ii < 0) || (ii >= MR_Obj.size_band_nl(s)) || 
	          (jj < 0) || (jj >= MR_Obj.size_band_nc(s)))
	     {
	         cout << "Errror: bad index in get_mrobj ... " << endl;
		 cout << "i,j = " << i << " " << j << " NlxNc = " << MR_Obj.size_band_nl(s) << " " << MR_Obj.size_band_nc(s) << endl;
		 exit(-1);
             }
 	     MR_Obj(s,ii,jj) = MR_Data(s,i,j);
	  }
       }
   }
#if DEBUG_SEG
 cout << "End get_mrobj " << endl;
#endif
}
 
/***********************************************************************/
/*
Bool MR2D_Seg::is_sub_struct(int StructFather, int CurrentStruct, int StructSub)
{ 
    Bool ValRet = False;
    float ValBef = (StructFather == 0) ? 0: TabStruct[StructFather].ValMax;
    float ValAfter = (StructSub == 0) ? 0:TabStruct[CurrentStruct].ValMax;
 
    float ValCurr = TabStruct[CurrentStruct].ValMax;

    // FIRST and SECOND LAWS:  
    
    if (ABS(ValCurr) > ABS(ValBef) && (ABS(ValCurr) > ABS(ValAfter)))
    {
       if (StructFather == 0) ValRet = False;
       else
       {
           int SF = StructFather;
	   int CS = CurrentStruct;
           int Px = TabStruct[SF].Maxx;
           int Py = TabStruct[SF].Maxy;  
           int PPx = pos_next(TabStruct[CS].Maxx, TabStruct[CS].Scale);
           int PPy = pos_next(TabStruct[CS].Maxy, TabStruct[CS].Scale);        	   
 cout << DistMin << " Px = " << Px << " PPx = " << PPx << " Py = " << Py << " PPy = " << PPy << endl;
 cout << PercentFlux  << " " << TabStruct[CS].Flux << " " << TabStruct[SF].Flux << endl;
 	   // THIRD LAW: test if the distance is larger than DistMin
 	   if ((ABS(PPx-Px) > DistMin) || (ABS(PPy-Py) > DistMin)) 
	   {
	      ValRet = True;   
	      //cout << "SUBOK" << endl; 
           }
	   //else cout << "KOKO" << endl;
	   //  Fourth LAW: test the flux
	   if ((PercentFlux > 0) && (ValRet == True)) 
	   {
	      float Norm = 1;
	      if (pyr_nextscale(TabStruct[CS].Scale) == True) Norm = 0.25;
              if (TabStruct[CS].Flux < TabStruct[SF].Flux*Norm*PercentFlux) ValRet = False;
 	   }
 	}
    }
    if (ValRet  == True) cout << "SUB" << endl; 
    else cout << "NOSUB" << endl;
    return ValRet;
}
*/
/***********************************************************************/

// Bool MR2D_Seg::obj_size_ok(int NumObj)
// {
//     Bool Val=True;
//     if (TabRel(NumObj,0) == 1)
//     {
//        int i = TabRel(i,1);
//        int SizeMin = 1<< (TabStruct[i].Scale+1);
//        if (TabStruct[i].Npix < SizeMin) Val=False;
//     }
//     return Val;
// }


/***********************************************************************/

// Bool MR2D_Seg::obj_nb_struct_ok(int NumObj, int MinStructNumber)
// {
//     Bool Val=True;
//     if (TabRel(NumObj,0) < MinStructNumber) Val=False;
//     return Val;
// }


/***********************************************************************/
/*
void MR2D_Seg::make_graph(Bool KillSize, Bool KillIsol, int MinStr) 
{ 
   int No;
   // obj2d *TabNewObj = new obj2d;
   int NewNbrObj=0;
     
   for (No=1; No <= NbrObj; No++)
   {
       int NbrStartNewObj = 0;
       int *TabStartNewObj = new int [TabNbrStructPerObj[No]];
       Bool BadObj = False;
         
       if ((KillIsol == True) && (obj_nb_struct_ok(No,MinStr) == False))
            BadObj = True;
       else if ((KillSize == True) && (obj_size_ok(No) == False))
            BadObj = True;

      // Count the number of subobj at each step
      // Need to known in how many objects this object has to be split
      // This can be done recursevely
      if (BadObj == False)
      {
      
      }
         
   } 
}

*/
/***********************************************************************/


/*
void MR2D_Seg::tree_sub_object(int StructFather, int CurrentStruct, Bool *TabPass)
{
     float Val = TabStruct[CurrentStruct].ValMax;
     int NbrRel,s,i,indrel;      
    
     // in order to not do it again
     TabPass[CurrentStruct] = True;
     
     NbrRel = TabRel(CurrentStruct,0); // Number of substructure 

     // No Father on this structure
     if (StructFather == 0)
     {
         
     }
     else
     {
         indrel = TabRel(i,s);
 	 float Max_Before = TabStruct[StructFather].ValMax;
	 
     }
 	     
     
     // look if any of the substructure contain a sub-object       
     for (s = 1; s <= NbrRel; s++)
     {
         indrel = TabRel(i,s);
	 tree_sub_object(CurrentStruct, indrel, TabPass);
     }     
}

*/  
/***********************************************************************/

/*
void MR2D_Seg::find_sub_struct()
{
  int i,s;
   
  // define if a structure S(j) is a starting point of a sub-object
  // to define a subobject, the laws are:
  // 1)  S(j) > S(j+1)
  // 2)  S(j) > S(j-1)
  // 3)  DIST[PosMax(j) - PosMax(j+1)] > DistMin
  // 4)  Flux(S(j)) > Percent*Flux(S(j+1))
  for (i = 1; i <= NbrStruct; i++)
  {
     float Val = TabStruct[i].ValMax;
     int Px = TabStruct[i].Maxx;
     int Py = TabStruct[i].Maxy;  
     int NbrRel = TabRel(i,0);
     
     // For each stucture at a better resolution connect to 
     // to the structure i
     for (s = 1; s <= NbrRel; s++)
     {
        int indrel = TabRel(i,s);
 	float Max_Before = TabStruct[indrel].ValMax;

        // FIRST LAW:  
        if (ABS(Val) < ABS(Max_Before))
	{
           int PPx = pos_next(TabStruct[indrel].Maxx, TabStruct[indrel].Scale);
           int PPy = pos_next(TabStruct[indrel].Maxy, TabStruct[indrel].Scale);        	   
	   Bool Sub=True;
	   
 	   // // THIRD LAW: test if the distance is larger than DistMin
 	   if ((ABS(PPx-Px) <= DistMin) && (ABS(PPy-Py) <= DistMin))  
	                   Sub = False;	
           //  test the flux
	   else  
	   {
	      float Norm = 1;
	      if (pyr_nextscale(TabStruct[indrel].Scale) == True) Norm = 4;
              if (TabStruct[i].Flux*Norm*PercentFlux < TabStruct[indrel].Flux)
	          Sub = False;
	   } 
	   
	   if (Sub == True)
	   TabStruct[s].SubStruct = Sub;
          // if (Sub ==True) TabStruct[indrel].SubObj = False;
        }
     } 
  }  
}

*/

/*
void MR2D_Seg::split_tree(float *TabColor)
{
  int i,j,s;
  int  NbrRel,NbrScale = Nbr_Plan;
  int StructFather,StructSub;
  Bool LocalMax;
  
  for (i = 1; i <= NbrStruct; i++)
   if (TabStruct[i].BadStruct == False)
   {
     NbrRel = TabRel(i,0); // Number of substructure 
     s = 0;
     LocalMax=False;
     while ((s < NbrRel) && (LocalMax == False))
     {
        s ++;
        StructFather = TabStruct[i].Father;
	StructSub = TabRel(i,s);
        LocalMax = local_max_in_obj(StructFather, i, StructSub);
     }
     
     if (LocalMax == True)  
     {
        if (TabStruct[i].Father == 0) RootObj = True;
	else
	
     
   }

}
*/

/***********************************************************************/
