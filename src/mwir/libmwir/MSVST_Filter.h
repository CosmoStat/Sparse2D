/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  Filter.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION : This class defines a filter bank
**    -----------   It can be created by different way:
**                 
**                  1) By using a predefined filter bank 
**                     ex: F_MALLAT_7_9,F_DAUBE_4, ...
**                  2) By using the file name if the file which
**                     contains the filter bank values to use.
**                     The format of the file is wvf format.
**                     see the Bath Wavelet Warehouse for more information:
**                          http://dmsun4.bath.ac.uk/wavelets/warehouse.html
**                    File Format (.wvf)
**                    General format to describe wavelet filter coefficients.
**
**                    [Range low] [Range high]
**                    [Analysis LP filter coefficients]
**                    .
**                    .
**                    .
**                    [Range low] [Range high]
**                    [Synthesis LP filter coefficients]
**                    .
**                    .
**                    .
**           
**                  3) By using the default user file name: mr1.wvf
**                     In this case, the file name must be in the 
**                     current directory
**                  4) By using the environment variable: CEA_FILTER
**                     
******************************************************************************/

#ifndef _DESIGN_MSVST_FILTER_H_
#define _DESIGN_MSVST_FILTER_H_

#include "GlobalInc.h"

// Type of filters

#define NBR_SB_FILTER 13 // total number of known subband (SB) filters 
enum type_sb_filter 
{
  SB_UNKNOWN, 
  F_HAAR,            // orthogonal Haar 2 points
  F_BI2HAAR,         // biorthogonal Haar 2/6
  F_BI4HAAR,         // biorthogonal Haar 2/10
  U_HAAR_B3S,
  U_HAAR_B3S2,
  F_DAUBE_4,         // orthogonal daubechies 4 points (order 2)
  F_DAUBE_8,         // orthogonal daubechies 8 points (order 4)
  F_3_5,
  F_5_3,
  F_MALLAT_7_9, 
  F_MALLAT_9_7,
  F_ODEGARD_9_7,
  F_USER 
};

#define DEF_SB_FILTER F_MALLAT_7_9            // Default user filter 
#define DEF_USER_FILTER_FILE_NAME  "mr1.wvf"  // Default user filter file name
#define USER_FILER_FILE_NAME    "CEA_FILTER"  // Environment variable name
                                              // indicating the filter file name

extern char* UserFilterFileName;

// given a subband filter type, return its description in a string
char* StringSBFilter (type_sb_filter  type);
// some basic usage info. of subband filters
void sb_usage(type_sb_filter Filter);
// user gives a string of the format "%d,[%s]". the first number is the type of the subband filter;
// the second is the UserFilterFileName
type_sb_filter get_filter_bank(char *UserArg);


enum sb_type_norm {NORM_L1, NORM_L2};
// default normalization.
#define DEF_SB_NORM NORM_L1  

class FilterAnaSynt
{
  sb_type_norm TypeNorm;       // normalization type of the filter
  // usually the normalization is L2;
  // in the case of F_3_5 and F_5_3 we use a L1 normalization.
  type_sb_filter TypeFilter;   // subband filter type
  float *Analysis, *Synthesis, \
    *Analysis2, *Synthesis2; // pointer which point to a float array of filter coefficients

  int Size_Ana, Size_Synt, \
    Size_Ana2, Size_Synt2;     // length of the analyse filter and that of the synthetic filter

  int Start_Ana, Start_Synt, \
    Start_Ana2, Start_Synt2;   // the discrete time point where the filter begins, usually zero but can be -length_of_filter/2 in the customized case
  void read_from_file(char *FileName);
  // read the filter H0 and H1 from a file
  // File format is:
  //  [Range low] [Range high]
  //  [Analysis LP filter coefficients]
  //  .
  //  .
  //  .
  //  [Range low] [Range high]
  //  [Synthesis LP filter coefficients]
  void reset_param(); // initialize all the parameters to zeros or nulls
 public:
  FilterAnaSynt() {reset_param();} // construct with SB_UNKNOWN filter type
  FilterAnaSynt(type_sb_filter T_Filter); // construct with a given filter type
  FilterAnaSynt(char *FileName); // construct with a filter file
  // filter type will be F_USER
  // check order : this->FilterFileName, UserFilterFileName
  // DEF_USER_FILTER_FILE_NAME, USER_FILER_FILE_NAME
  ~FilterAnaSynt(); // destruct with all the data fields of this class reinitialized
  
  type_sb_filter type_filter()  {return  TypeFilter;}
  sb_type_norm norm()   {return TypeNorm;}
  float *analysis()     {return Analysis;}
  float *synthesis()    {return Synthesis;}
  float *analysis2()     {return Analysis2;}
  float *synthesis2()    {return Synthesis2;}
  int size_analysis()   {return Size_Ana;}
  int size_synthesis()  {return Size_Synt;}
  int size_analysis2()   {return Size_Ana2;}
  int size_synthesis2()  {return Size_Synt2;}
  int start_analysis()  {return Start_Ana;}
  int start_synthesis() {return Start_Synt;}
  int start_analysis2()  {return Start_Ana2;}
  int start_synthesis2() {return Start_Synt2;}
    
  // gives all the data fields of this class 
  // the corresponding values, given the subband filter type
  void alloc(type_sb_filter T_Filter); 

  Bool Verbose;
  char *FilterFileName; // User filter bank file name
  // This field is used when T_Filter=F_USER
};


#endif
