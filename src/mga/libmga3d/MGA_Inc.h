
#ifndef _MGA_INC_
#define _MGA_INC_

#include "writefits3d.h"
#include <fstream>


// type enumerations and strings associated
	enum type_3dcurvelet_block {BL3D_CONST,BL3D_UP,BL3D_UP2,BL3D_DOWN, BL3D_DOWN2};
	#define NBR_3DCUR_TYPE_BLOCK 5
	#define DEF_3DCUR_NBR_SCALE 4 // Default number of scale in the a trous algo
	#define DEF_3DCUR_BLOCK_SIZE 17  // Default block size at the first scale
	#define DEF_3DCUR_TYPE_BLOCK BL3D_UP2
	const char * StringBlockType (type_3dcurvelet_block type);

	// enum filter_type {FT_UNKNOWN, FT_HARD, FT_SOFT, FT_WIENER, FT_FDR, FT_SBT, FT_CONTRAST};// FilterType
	// char* string_filter_type(filter_type t);

	enum type_wavelet3D {W3D_UNKNOWN, W3D_ATROU, W3D_MEYER, W3D_PMEYER};
	#define NBR_TYPE_W3D 3
	#define DEF_TYPE_W3D W3D_MEYER
	enum type_pmeyer3D {PMW3D_NONE, PMW3D_COHERENT, PMW3D_SUM};
	#define DEF_TYPE_PMW3D PMW3D_SUM

	enum type_extraction {EXTRACT_B, EXTRACT_F, EXTRACT_FB};
	#define DEF_TYPE_EXTRACTION EXTRACT_FB

// String manipulation
	char *add_fits(char * NameStep);
	char *add_mr(char * NameStep);
	
#endif
