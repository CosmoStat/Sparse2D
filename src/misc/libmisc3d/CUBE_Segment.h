#ifndef _SEG3D_H_
#define _SEG3D_H_

#include "IM_Obj.h"
#include "IM3D_IO.h"
#include "IM_IO.h"			
#define		PILE_SIZE_BLOC		100000


void cube_segment (fltarray & cube, intarray & Segment, 
	int &nb_seg, float S, Bool CleanBord, int FirstLabel);	
void cube_segment (intarray & cube, intarray & Segment, 
	int &nb_seg, float S, Bool CleanBord, int FirstLabel);	
	
class Pile {
	public:
		int * data;
		int num; 
		int bloc;
		Pile(void);
		Pile (int n);
		~Pile (void);
		int	max (void);
		void	alloue (void);
		void	add_data (int d);
		void	trie (int max);
		void	remplace (int d, int f);
};

#endif
