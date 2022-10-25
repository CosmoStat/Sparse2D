
#include <omp.h>
#include <iostream>
// #include "Array.h"
#include "Array.h"

#ifdef _OPENMP
    int inner_loop_threads;
    int outer_loop_threads;
#endif

/*********************************************************************/
 
inline static void openmp_params() {
	
   #ifdef _OPENMP    
 	extern int inner_loop_threads;
  	extern int outer_loop_threads;
	    cout << "Application was compiled with OpenMP support," << endl;
    	if (omp_get_max_threads() == 1)
      	cout << "but running with one process only." << endl;
    	else
    	  cout << "running with up to " << omp_get_max_threads()
           << " processes." << endl;
           if (omp_get_nested()){
           		cout << "Nested Parallelism allowed" <<endl;
           		cout << "running outer loops with " << outer_loop_threads
           << " processes." << endl;
           	cout << "running inner loops with " << inner_loop_threads
           << " processes." << endl;

           } else cout << "running with" << outer_loop_threads
           << " processes." << endl;
    #else
     cout << "Application was compiled without OpenMP support;" << endl
         << "running in scalar mode." << endl;
    #endif
  }; 
  
  
 inline  static void set_nb_threads(int outer, int inner) {
     #ifdef _OPENMP    

  	extern int inner_loop_threads;
  	extern int outer_loop_threads;
	
  	int nthr=omp_get_max_threads();

	if(outer * inner > nthr) {
		printf("Too many threads required (%d * %d vs %d procs)\n",outer,inner,nthr);
		exit(EXIT_FAILURE);
	} else {
		if(outer <= 0) outer_loop_threads=1; else outer_loop_threads=outer;
		if(inner <=0)  inner_loop_threads=1; else inner_loop_threads=inner;
	}
	omp_set_num_threads(outer_loop_threads*inner_loop_threads);	
    #endif

  };
  
  static void auto_set_inner_outer_threads(int outer_max, int inner_max,int nb_procs=0) {
	//Simple rules for setting the number of outer/inner threads:
	//First set the maximal number of procs available if not set explicitely:
	//1) if nprocs  < 4, then we use all of them 
	//2) else use only half of them (we do not need so much computing power)
	//Then dispatch the procs in outer loops and inner loops knowing that:
	//outer_max (inner_max) is the maximal number of threads useful for outer (inner) loops, set to the nb of procs available  if equal to 0.
	//The criterion is the following:
	//Maximize the number of procs used, and if several combinations, maximize the number of proc for outer loops threads (where we expect more computing time gain)
  	int nbprocs,nthr;	
  	omp_set_nested(1);
	
	nbprocs=omp_get_num_procs();	
	
	if(nb_procs==0) {
		if(nbprocs<4) omp_set_num_threads(nbprocs);	 
		else omp_set_num_threads((int) (nbprocs/2)); //Use only 50% of procs
	} else {
		if(nb_procs > nbprocs) 	{
			printf("Too many procs required (%d required vs %d avail)\n",nb_procs,nbprocs);
			exit(EXIT_FAILURE);
		} else {
			nbprocs=nb_procs;
			omp_set_num_threads((int) (nbprocs)); 
		}
	}
	nthr=omp_get_max_threads();
	if(!omp_get_nested()) set_nb_threads(nthr,1);
	else {
		if((outer_max > nthr) || (outer_max <=1)) outer_max=nthr;
		if((inner_max > nthr) || (inner_max <=1)) inner_max=nthr;
	
		int k,divisor,inner,outer;
		intarray TotRed(outer_max);
		int MaxPotRed=0;
		for(k=outer_max-1;k>=0;k--) {
			divisor=(int) floor(nthr /(k+1));
			if((divisor <= inner_max) && (divisor*(k+1) > MaxPotRed)) {
				MaxPotRed=divisor*(k+1);
				outer=k+1;
				inner=divisor;
			}
		} 
		set_nb_threads(outer,inner);
	}
	printf("Use %d outer threads and %d inner threads for nb procs=%d\n",outer_loop_threads,inner_loop_threads,nbprocs);
};
