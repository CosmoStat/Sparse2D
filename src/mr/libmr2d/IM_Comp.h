/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: J.L. Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_Comp.h
**
*******************************************************************************
**
**    DECRIPTION    Header for image compression
**    ---------- 
**
*****************************************************************************/

#ifndef __IM_COMP__
#define __IM_COMP__

#define NBR_COMP_METHOD 7
enum type_comp {COMP_MRMEDIAN, COMP_MORPHO, COMP_HAAR, COMP_MIN_MAX,
                COMP_MALLAT, COMP_FEAUVEAU, COMP_WT_PMT, 
	        COMP_O_MRMEDIAN, COMP_O_MORPHO,COMP_O_MIN_MAX, COMP_INT_LIFTING};

enum type_quant {Q_UNIFORM, Q_NON_UNIFORM};
enum type_coding {C_HUFFMAN, C_ARITHM};

char * StringComp (type_comp type);


inline Bool work_with_int(type_comp  Comp_Method)
{
   Bool IsInt = False;
   if (  (Comp_Method == COMP_INT_LIFTING ) || 
         (Comp_Method == COMP_O_MORPHO) || 
        (Comp_Method == COMP_O_MIN_MAX ) || 
        (Comp_Method == COMP_O_MRMEDIAN )) IsInt = True;
       //  (Comp_Method == COMP_MORPHO))  IsInt = True;
   return IsInt;
}

inline Bool noise_threshold(type_comp Comp_Method)
{
   Bool Th = True;
   if (  (Comp_Method == COMP_INT_LIFTING ) || 
 	 (Comp_Method == COMP_O_MORPHO) ||
         (Comp_Method == COMP_MORPHO) ||
         (Comp_Method == COMP_O_MIN_MAX))  Th = False;  
   return Th;
}

void encode(FILE *outfile, int *a, int nx, int ny, int scale);
/* outfile=	Output file	
   a[]=  input  array (nx,ny)	
   nx,ny= size of H-transform array		
   scale= scale factor for digitization	*/


void decode(FILE *infile, int  **a, int  *nx, int  *ny, int  *scale);
/*FILE *infile;  input file 
int  **a;  address of output array [nx][ny]	
int  *nx,*ny;	 size of output array	
int  *scale;  scale factor for digitization */


#endif
