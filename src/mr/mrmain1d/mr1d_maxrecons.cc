/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.1
**
**    Author: 
**
**    Date:  98/07/30
**    
**    File:  mr1d_maxrecons
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution recons from max of signal transform 
**    ----------- 
**                 
**   Usage: mr1d_maxrecons options image output
**
**   where options =
**
**         [-i number_of_iterations]
**              Maximum number of iterations
**
**         [-z type_of_constraint]
**              1: E_STEPH_MALLAT
**              2: E_PHIL
**              default : E_STEPH_MALLAT
**
**         [-v]
**				Verbose 
**
******************************************************************************/

 
#include "GlobalInc.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "MR1D_Obj.h"

#define DEF_ITER 10
#define OK 0
#define NOT_OK 1
#define EPS 1e-07

char stc_NameImagIn[256];      /* input file image */
char stc_NameImagOut[256];     /* output file image */

//int si_NbrPlan=0;              /* nb plan of decomp */
//Bool se_SetPlan = False;  
Bool se_SetIter = False;  
int si_MaxIter=DEF_ITER;  
Bool Verbose = False;     

type_trans_1d se_Transform=TO1_PAVE_B3SPLINE;

typedef enum { E_SANS         =0,
               E_STEPH_MALLAT =1,
	       E_PHIL         =2}
te_constraint;
static te_constraint se_Constraint=E_STEPH_MALLAT;

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

/*********************************************************************/

//static int max_scale_number (int pi_Nc) {
//
//    int ai_Nmin = pi_Nc;
//    int ai_ScaleMax;
//
//	ai_ScaleMax = iround(log((float)ai_Nmin / (float)MAX_SIZE_LAST_SCALE) 
//                                            / log(2.)+ 1.);
//	return (ai_ScaleMax);
//}

/*********************************************************************/

static void usage(char *argv[])
{

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-i number_of_iterations]\n");
    fprintf(OUTMAN, "              number of iteration\n");
    fprintf(OUTMAN, "              default is %d\n", si_MaxIter);
    manline();

//    fprintf(OUTMAN, "         [-z type_of_constraint] \n");
//    fprintf(OUTMAN, "              1: E_STEPH_MALLAT  \n");
//    fprintf(OUTMAN, "              2: E_PHIL          \n");
//    fprintf(OUTMAN, "              default : E_STEPH_MALLAT \n");
//    manline();

    verbose_usage();
    manline();

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[]) {
    int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"i:z:v")) != -1) {
	
		switch (c) {
//		case 't': /* -d <type> type of transform */
//			if (sscanf(OptArg,"%d",&c ) != 1) {
//		    		fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
//	            		exit(-1);  
//			}
//                	if ((c > 0) && (c <= NBR_TRANS_1D)) se_Transform = (type_trans_1d) (c-1);
//                	else {
//		    		fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
//	            		exit(-1);
// 			}
//			break;
		case 'i': /* -i < Number of iterations> */
			if (sscanf(OptArg,"%d",&si_MaxIter) != 1) {
				fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
				usage(argv);
			}
			if (si_MaxIter <= 0)  si_MaxIter = DEF_ITER;
			se_SetIter=True;
			break;
		case 'z': /* -z <type_of_constraint> */
			if (sscanf(OptArg,"%d",(int*)&se_Constraint) != 1) {
				fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
				usage(argv);
			}
			if ((se_Constraint < 0) ||  (se_Constraint > 2)) 
				se_Constraint = E_STEPH_MALLAT;
			break;
		case 'v':
			Verbose = True;
			break;
		case '?':
			usage(argv);
		}
	}

	/* get optional input file names from trailing parameters and open files */
	if (OptInd < argc) strcpy(stc_NameImagIn, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(stc_NameImagOut, argv[OptInd++]);
	else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}





/******************************************************************************/
//extern void wave_1d_B3deriv_atrou (fltarray &Signal, fltarray & W_1D,   
//                                   int Nbr_Plan, type_border Border);
//extern void wave_1d_rec_B3deriv_atrou (fltarray &W_1D, fltarray & Signal, 
//                                        int Nbr_Plan, type_border Border);
/******************************************************************************/



void constraint_one (int pi_Scale,            // current scale
                     int pi_Begining,         // ind of begining
		     int pi_End,              // ind of end
                     MR_1D& po_Mr1dData,      // Data IN
                     MR_1D& po_Mr1dRecData){  // Rec Data OUT

  // sign of two consecutive maxima
  int ai_SignFisrtMax = (po_Mr1dData(pi_Scale,pi_Begining) >= 0) ? 1 : -1;
  int ai_SignSecondMax = (po_Mr1dData(pi_Scale,pi_End) >= 0) ? 1 : -1;

  if (ai_SignFisrtMax == ai_SignSecondMax) {
    // same sign for po_Mr1dRecData than max
    for (int i=pi_Begining+1; i<pi_End; i++)
      po_Mr1dRecData(pi_Scale,i) = ABS(po_Mr1dRecData(pi_Scale,i))*ai_SignFisrtMax;
  } else {
    for (int i=pi_Begining+2; i<pi_End; i++) {
      if (    (ai_SignFisrtMax > 0)             // first max > 0 => g() decrease
	   && (po_Mr1dRecData(pi_Scale,i) > po_Mr1dRecData(pi_Scale,i-1))) {
	    //cout << "Pb de sens de variation" << endl;
	    	    
	    po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1);
	    
	    /*int k=i+1;
	    float af_DemiSomme = (po_Mr1dRecData(pi_Scale,i-1) + po_Mr1dRecData(pi_Scale,k))/2;
	    while (af_DemiSomme > po_Mr1dRecData(pi_Scale,i) && k<pi_End) {
	       k++;
	       af_DemiSomme = (po_Mr1dRecData(pi_Scale,i-1) + po_Mr1dRecData(pi_Scale,k))/2;
	    } 
	    po_Mr1dRecData(pi_Scale,i) = af_DemiSomme;*/
	    
      } 
      if (    (ai_SignFisrtMax < 0)             // first max < 0 => g() increase
	   && (po_Mr1dRecData(pi_Scale,i) < po_Mr1dRecData(pi_Scale,i-1))) {
	    //cout << "Pb de sens de variation" << endl;
	    
	    po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1);
	    
	    /*int k=i+1;
	    float af_DemiSomme = (po_Mr1dRecData(pi_Scale,i-1) + po_Mr1dRecData(pi_Scale,k))/2;
	    while (af_DemiSomme < po_Mr1dRecData(pi_Scale,i) && k<pi_End) {
	       k++;
	       af_DemiSomme = (po_Mr1dRecData(pi_Scale,i-1) + po_Mr1dRecData(pi_Scale,k))/2;
	    } 
	    po_Mr1dRecData(pi_Scale,i) = af_DemiSomme;*/	    
	    
	    
      } 
    }
  }
  
  if (po_Mr1dData(pi_Scale,pi_End) > po_Mr1dData(pi_Scale,pi_Begining)) {
    for (int i=pi_Begining+1; i<pi_End; i++) {
      if (po_Mr1dRecData(pi_Scale,i) < po_Mr1dData(pi_Scale,pi_Begining))
        po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_Begining);
      if (po_Mr1dRecData(pi_Scale,i) > po_Mr1dData(pi_Scale,pi_End))      
         po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_End);
    } 
  } 
  else if (po_Mr1dData(pi_Scale,pi_End) < po_Mr1dData(pi_Scale,pi_Begining)) {
    for (int i=pi_Begining+1; i<pi_End; i++) {
      if (po_Mr1dRecData(pi_Scale,i) > po_Mr1dData(pi_Scale,pi_Begining))
        po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_Begining);
      if (po_Mr1dRecData(pi_Scale,i) < po_Mr1dData(pi_Scale,pi_End))      
         po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_End);
    } 
  }       
}

void constraint_two (int pi_Scale,            // current scale
                     int pi_Begining,         // ind of begining
		             int pi_End,              // ind of end
                     MR_1D& po_Mr1dData,      // Data IN
                     MR_1D& po_Mr1dRecData){  // Rec Data OUT
 
  // sign of first max - second max
  int ai_sign = 
    (po_Mr1dData(pi_Scale,pi_Begining)-po_Mr1dData(pi_Scale,pi_End) >= 0) ? 1 : -1;

  // for each point
  for (int i=pi_Begining+2; i<pi_End; i++) {

    if ( (ai_sign > 0) &&
         (po_Mr1dRecData(pi_Scale,i) > po_Mr1dRecData(pi_Scale,i-1))) {
      //cout << "Pb de sens de variation" << endl;
      // if level at i+1 < at level at i-1=> interpolate level at i = (lev(i-1)+lev(i+1))/2
      if (po_Mr1dRecData(pi_Scale,i+1) < po_Mr1dRecData(pi_Scale,i-1))
         po_Mr1dRecData(pi_Scale,i) = (po_Mr1dRecData(pi_Scale,i-1)+po_Mr1dRecData(pi_Scale,i+1))/2.;
      else
         po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1);
    } 

    if ( (ai_sign < 0) &&
         (po_Mr1dRecData(pi_Scale,i) < po_Mr1dRecData(pi_Scale,i-1))) {
      //cout << "Pb de sens de variation" << endl;
      // if level at i+1 > at level at i-1=> interpolate level at i = (lev(i-1)+lev(i+1))/2
      if (po_Mr1dRecData(pi_Scale,i+1) > po_Mr1dRecData(pi_Scale,i-1))
         po_Mr1dRecData(pi_Scale,i) = (po_Mr1dRecData(pi_Scale,i-1)+po_Mr1dRecData(pi_Scale,i+1))/2.;
      po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1);
    } 
  }
}




// interpolate between two consecutive maxima modulus
void Ortho_Proj_Operator (int pi_Scale,            // current scale
                          int pi_Begining,         // ind of begining
	                  int pi_End,              // ind of end
                          MR_1D& po_Mr1dData,      // Data IN
                          MR_1D& po_Mr1dRecData){  // Rec Data OUT
			  
	//cout << "sca:"<<pi_Scale<<",beg:"<<pi_Begining<<",end:"<<pi_End<<endl;
			  
  int pi_NumScale = pi_Scale+1;   // pi_NumScale goes from 1 to Max, 
                                  // and pi_Scale from 0 to max-1
  pi_End = pi_End - pi_Begining;
  int ai_Dec = pi_Begining;
  pi_Begining = 0;

  // coef of interpolate function
  double af_Step1 = pow((double)2., (double)(-pi_NumScale))*pi_Begining;
  double af_Step2 = pow((double)2., (double)(-pi_NumScale))*pi_End;
  double ad_Delta = exp(af_Step1) * exp(-af_Step2) 
                  - exp(af_Step2) * exp(-af_Step1);
		  
  int ai_FirstInd = pi_Begining+ai_Dec;
  int ai_LastInd  = pi_End+ai_Dec;
		  
  double ad_DeltAlpha = 
   (po_Mr1dData(pi_Scale, ai_FirstInd) - po_Mr1dRecData(pi_Scale, ai_FirstInd))
    * exp(-af_Step2) -  exp(-af_Step1) *
   (po_Mr1dData(pi_Scale, ai_LastInd) - po_Mr1dRecData(pi_Scale, ai_LastInd));
  
  double ad_DeltBeta  = 
    (po_Mr1dData(pi_Scale, ai_LastInd) - po_Mr1dRecData(pi_Scale, ai_LastInd))
     * exp(af_Step1) - exp(af_Step2) *
    (po_Mr1dData(pi_Scale, ai_FirstInd) - po_Mr1dRecData(pi_Scale, ai_FirstInd));
  
  double ad_Alpha = ad_DeltAlpha / ad_Delta;
  double ad_Beta  = ad_DeltBeta  / ad_Delta;

  // interpolate between pi_Begining and pi_End
  for (int i=pi_Begining+1; i<pi_End; i++) {
      double ad_inter             = exp (pow((double)2., (double)(-pi_NumScale))*i)
                                  * ad_Alpha + ad_Beta *
                                  exp (- pow((double)2., (double)(-pi_NumScale))*i);
      //double ad_inter1 = po_Mr1dRecData(pi_Scale,i+ai_Dec);
      po_Mr1dRecData(pi_Scale,i+ai_Dec) += ad_inter;
      //cout <<"Data IN  scale:"<<pi_Scale <<", niv:"<<po_Mr1dData(pi_Scale,i+ai_Dec)<<endl;
      //cout <<"Data OUT scale:" << pi_Scale <<", niv:"<<ad_inter1+ad_inter;
  }


  if (se_Constraint == E_STEPH_MALLAT) 
    constraint_one (pi_Scale, pi_Begining+ai_Dec, pi_End+ai_Dec,
                    po_Mr1dData, po_Mr1dRecData);

  else if (se_Constraint == E_PHIL) 
    constraint_two (pi_Scale, pi_Begining+ai_Dec, pi_End+ai_Dec,
                    po_Mr1dData, po_Mr1dRecData);
		    
  else 
    cout << "No constraint!" << endl;		 

}



void LinearOrthoProj (int pi_Scale,            // current scale
                      int pi_Begining,         // ind of begining
	              int pi_End,              // ind of end
                      MR_1D& po_Mr1dData,      // Data IN
                      MR_1D& po_Mr1dRecData){  // Rec Data OUT



  float a = (po_Mr1dData(pi_Scale, pi_End)-po_Mr1dData(pi_Scale, pi_Begining)) / (pi_End-pi_Begining);
  float b = po_Mr1dData(pi_Scale, pi_Begining) - a * pi_Begining;
   
   
  // interpolate between pi_Begining and pi_End
  for (int i=pi_Begining+1; i<pi_End; i++) {

      po_Mr1dRecData(pi_Scale,i) = a*i + b;

  }
}


/*void interpolate  (int pi_NbrPlan,           // number of plan
                   int pi_Nl,                // number of elt
                   MR_1D& po_Mr1dData,       // Data IN
                   MR_1D& po_Mr1dRecData){   // Rec Data OUT

  // cal index of maxima of Mr1dData
  for (int s=0; s<pi_NbrPlan; s++) {         // for each scale
    int ai_NbrMaxCurrentScale = 0;           // number of max at scale s
    for (int k=0; k<pi_Nl; k++) {            // for each point, count number of max
      if (po_Mr1dData(s,k) < -EPS || po_Mr1dData(s,k) > EPS) ai_NbrMaxCurrentScale++;
    }
    if (ai_NbrMaxCurrentScale > 0) {
      intarray ao_VectIndex(ai_NbrMaxCurrentScale);
      int j=0;                               // local ind of ao_VectIndex
      for (int i=0; i<pi_Nl; i++) {          // for each point
	    if (po_Mr1dData(s,i) < -EPS || po_Mr1dData(s,i) > EPS) ao_VectIndex(j++)=i;
      }
      //cout << "number of max :" << ai_NbrMaxCurrentScale << " , scale :" << s << endl;
      // for each index of local maxima
      for (int l=0; l<ai_NbrMaxCurrentScale-1;l++){
	    Ortho_Proj_Operator (s, ao_VectIndex(l), ao_VectIndex(l+1),
                             po_Mr1dData, po_Mr1dRecData);
      }
    }
  }
}*/
 
void init_last_scale (int pi_NbrPlan,           // number of plan
                      int pi_Nl,                // number of elt
                      MR_1D& po_Mr1dData,       // Data IN
                      MR_1D& po_Mr1dRecData){   // Rec Data OUT
   for (int i=0; i<pi_Nl; i++) 
     po_Mr1dRecData (pi_NbrPlan-1,i) = po_Mr1dData(pi_NbrPlan-1,i);
}

/*void init_max (int pi_NbrPlan,           // number of plan
               int pi_Nl,                // number of elt
               MR_1D& po_Mr1dData,       // Data IN
               MR_1D& po_Mr1dRecData){   // Rec Data OUT
  for (int s=0; s<pi_NbrPlan-1; s++)
    for (int i=0; i<pi_Nl; i++) 
      if (po_Mr1dData(s,i) < -EPS || po_Mr1dData(s,i) > EPS) {
	//cout << "IN:" <<po_Mr1dData(s,i)<<", OUT:"<<po_Mr1dRecData(s,i);
	po_Mr1dRecData(s,i) = po_Mr1dData(s,i) ;
      }
}*/

/*void cut_local_max2 (int pi_Scale,              
                    int pi_Begining, 
		    int pi_End,           
                    MR_1D& pro_Mr1dRecData) { 
		    
	       
  float af_LevelBeg = pro_Mr1dRecData(pi_Scale, pi_Begining);
  float af_LevelEnd = pro_Mr1dRecData(pi_Scale, pi_End);
  
  int ai_sign = (af_LevelBeg - af_LevelEnd >= 0) ? -1 : +1;
  
  Bool ae_EndLoop=false;
  while (!ae_EndLoop) {

    ae_EndLoop = true;
    for (int i=pi_Begining+1; i<=pi_End; i++ ) {
     
      if (   (pro_Mr1dRecData(pi_Scale,i) > pro_Mr1dRecData(pi_Scale,i-1) && ai_sign == -1)
           || (pro_Mr1dRecData(pi_Scale,i) < pro_Mr1dRecData(pi_Scale,i-1) && ai_sign == +1)) {
      
        // onr modif	       
        ae_EndLoop = true;      
	
        // init val at i
        pro_Mr1dRecData(pi_Scale,i) = (pro_Mr1dRecData(pi_Scale,i-1)+
                                       pro_Mr1dRecData(pi_Scale,i-1)) / 2.;
      }			     
    }
  }      	       
}

void cut_local_max (int pi_Scale,              
                    int pi_Begining, 
		    int pi_End,           
                    MR_1D& pro_Mr1dRecData) { 
		    
  int ai_IndFirstLocalMax = pi_Begining;
  int ai_IndLastLocalMax = pi_End;
		    
  float af_LevelBeg = pro_Mr1dRecData(pi_Scale, pi_Begining);
  float af_LevelEnd = pro_Mr1dRecData(pi_Scale, pi_End);
  
  int ai_sign = (af_LevelBeg - af_LevelEnd >= 0) ? -1 : +1;
  
  for (int i=pi_Begining+1; i<=pi_End; i++ ) {
 
    if (   (pro_Mr1dRecData(pi_Scale,i) > pro_Mr1dRecData(pi_Scale,i-1) && ai_sign == -1)
        || (pro_Mr1dRecData(pi_Scale,i) < pro_Mr1dRecData(pi_Scale,i-1) && ai_sign == +1)) {
      
      // ind first local max
      ai_IndFirstLocalMax = i-1;
      
      // search for Last local max
      int j = pi_End;
      while (   (pro_Mr1dRecData(pi_Scale,j-1) <= pro_Mr1dRecData(pi_Scale,j) && ai_sign == -1)
              ||(pro_Mr1dRecData(pi_Scale,j-1) >= pro_Mr1dRecData(pi_Scale,j) && ai_sign == +1)) {
        j--;
      }
      ai_IndLastLocalMax = j;
      
      // compute mean
      float af_Som=0;
      for (int k=ai_IndFirstLocalMax; k<ai_IndLastLocalMax+1; k++)
        af_Som += pro_Mr1dRecData(pi_Scale,k);
      float af_Mean = af_Som / (ai_IndLastLocalMax-ai_IndFirstLocalMax+1);
      
      // init with the min
      for (int k=ai_IndFirstLocalMax; k<ai_IndLastLocalMax+1; k++)
        pro_Mr1dRecData(pi_Scale,k) = af_Mean;     

      // left border
      j=ai_IndFirstLocalMax;
      while (   (pro_Mr1dRecData(pi_Scale,j) <= af_Mean && ai_sign == -1)
              ||(pro_Mr1dRecData(pi_Scale,j) >= af_Mean && ai_sign == +1)) {
	 pro_Mr1dRecData(pi_Scale,j) = af_Mean;
	 j--;
      }	      
	      
      // rigth border
      j=ai_IndLastLocalMax;
      while (   (pro_Mr1dRecData(pi_Scale,j) >= af_Mean && ai_sign == -1)
              ||(pro_Mr1dRecData(pi_Scale,j) <= af_Mean && ai_sign == +1)) {
	 pro_Mr1dRecData(pi_Scale,j) = af_Mean;
	 j++;
      }	      	    
	    
      break;
    }
  }	  
}*/	  


void search_max (int pi_Ny,               // number of plan
                 int pi_Nx,               // number of elt
                 MR_1D& pro_Mr1dData, 
                 intarray& pro_NumberMax) {
		 
		 
  for (int s=0; s<pi_Ny-1; s++) {
    pro_NumberMax(s)++; //first point
    for (int i=0; i<pi_Nx; i++) 
      if (pro_Mr1dData(s,i) < -EPS || pro_Mr1dData(s,i) > EPS) 
        pro_NumberMax(s)++;
	  
    pro_NumberMax(s)++; //last point
    if (Verbose)
       cout << "Number max at scale " << s << " : " << pro_NumberMax(s) << endl;  
  }
}

void montone_function (int pi_Scale,    // Scale number 
               int pi_Nx,               // number of elt
	       intarray& pro_NumberMax, // Number of max 
               intarray** ppo_MaxInd,   // loc of max 	       
               MR_1D& po_Mr1dData,      // Data IN
               MR_1D& po_Mr1dRecData) { // Rec Data OUT 
	       
  // interpolate between max at scale s
  for (int l=0; l<pro_NumberMax(pi_Scale)-1; l++) {
  
    int pi_Begining = (**(ppo_MaxInd+pi_Scale))(l);
    int pi_End      = (**(ppo_MaxInd+pi_Scale))(l+1);
  
    if (po_Mr1dData(pi_Scale,pi_End) > po_Mr1dData(pi_Scale,pi_Begining)) {
      for (int i=pi_Begining+1; i<pi_End; i++) {
      
        if (po_Mr1dRecData(pi_Scale,i) > po_Mr1dData(pi_Scale,pi_End))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_End);
        if (po_Mr1dRecData(pi_Scale,i) < po_Mr1dData(pi_Scale,pi_Begining))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_Begining);        
	
	/*if (po_Mr1dRecData(pi_Scale,i) < po_Mr1dRecData(pi_Scale,i-1))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1) ;*/
      } 
      //cut_local_max2 (pi_Scale, pi_Begining, pi_End, po_Mr1dRecData);
        
    } else if (po_Mr1dData(pi_Scale,pi_End) < po_Mr1dData(pi_Scale,pi_Begining)) {
      for (int i=pi_Begining+1; i<pi_End; i++) {
      
        if (po_Mr1dRecData(pi_Scale,i) < po_Mr1dData(pi_Scale,pi_End))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_End);
        if (po_Mr1dRecData(pi_Scale,i) > po_Mr1dData(pi_Scale,pi_Begining))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dData(pi_Scale,pi_Begining);	  
	   
        /*if (po_Mr1dRecData(pi_Scale,i) > po_Mr1dRecData(pi_Scale,i-1))
	  po_Mr1dRecData(pi_Scale,i) = po_Mr1dRecData(pi_Scale,i-1) ;*/	         
      } 
      //cut_local_max2 (pi_Scale, pi_Begining, pi_End, po_Mr1dRecData);
    }
  }           	       
}




void init_max (int pi_Ny,               // number of plan
               int pi_Nx,               // number of elt
	       intarray& pro_NumberMax, // Number of max 
               intarray** ppo_MaxInd,    // loc of max 	       
               MR_1D& po_Mr1dData,      // Data IN
               MR_1D& po_Mr1dRecData){  // Rec Data OUT 
	       	       
  for (int s=0; s<pi_Ny-1; s++) {	       
      
    for (int l=0; l<pro_NumberMax(s); l++)
      if (pro_NumberMax(s) > 0)
        po_Mr1dRecData (s, (**(ppo_MaxInd+s))(l)) = po_Mr1dData (s, (**(ppo_MaxInd+s))(l));
  

    montone_function (s, pi_Nx, pro_NumberMax, ppo_MaxInd, po_Mr1dData, po_Mr1dRecData);
  
    //for (int l=0; l<pro_NumberMax(s); l++)
    //  if (pro_NumberMax(s) > 0)
    //    po_Mr1dRecData (s, (**(ppo_MaxInd+s))(l)) = po_Mr1dData (s, (**(ppo_MaxInd+s))(l));
  } 
}	       
	       
void interpolate  (int pi_Ny,               // number of plan
                   int pi_Nx,               // number of elt
	           intarray& pro_NumberMax, // Number of max 
                   intarray** ppo_MaxInd,    // loc of max 	     		   
                   MR_1D& po_Mr1dData,      // Data IN
                   MR_1D& po_Mr1dRecData){  // Rec Data OUT

  // cal index of maxima of Mr1dData
  for (int s=0; s<pi_Ny-1; s++) {         // for each scale
  
    // interpolate between max at scale s
    for (int l=0; l<pro_NumberMax(s)-1; l++) {
    
      int ai_begin = (**(ppo_MaxInd+s))(l);
      int ai_end   = (**(ppo_MaxInd+s))(l+1);
      
      //if (ai_end-ai_begin == 2) {
      //   po_Mr1dRecData(s,ai_begin+1) = (po_Mr1dRecData(s,ai_begin)+po_Mr1dRecData(s,ai_end))/2.;
      //} else 
      if (ai_end-ai_begin >= 2) {   
    
         Ortho_Proj_Operator (s, 
                              (**(ppo_MaxInd+s))(l), (**(ppo_MaxInd+s))(l+1),
      		              po_Mr1dData, po_Mr1dRecData);
      
      //LinearOrthoProj (s, 
      //                 (**(ppo_MaxInd+s))(l), (**(ppo_MaxInd+s))(l+1),
      //		       po_Mr1dData, po_Mr1dRecData);
      }
    }
  }
}	       
	            

/******************************************************************************/



int main(int argc, char *argv[])
{

  // local var
  fltarray  ao_2dData;
  fitsstruct FitsHeader;
  int s;
  
  //char atc_FileName[20];

  /* Get command line arguments, open input file(s) if necessary */
  lm_check(LIC_M1D);
  transinit(argc, argv);
 
  // read the data
  // -------------
  // fits_read_fltarr (stc_NameImagIn, ao_2dData, &FitsHeader);
  type_1d_format  DatFormat = io_detect_1dformat(stc_NameImagIn);
  if (DatFormat == F1D_FITS) 
      fits_read_fltarr (stc_NameImagIn, ao_2dData, &FitsHeader);
  else io_read2d_ascii(stc_NameImagIn, ao_2dData);

  int ai_Nx = ao_2dData.nx();
  int ai_Ny = ao_2dData.ny();
  
  // change format of data, and init
  MR_1D ao_Mr1dData (ai_Nx, se_Transform, "Mr1d Transform", ai_Ny);  
  ao_Mr1dData.Border = I_MIRROR;  
  for (s=0; s<ai_Ny; s++)
    for (int i=0; i<ai_Nx; i++)
      ao_Mr1dData(s,i) = ao_2dData (s*ai_Nx+i);  
      
      
  // search ind of max
  // ----------------- 
  intarray ao_NumberMax (ai_Ny);
  search_max (ai_Ny, ai_Nx, ao_Mr1dData, ao_NumberMax);
  //intarray* atpo_MaxInd [ai_Ny];
  intarray** atpo_MaxInd = new intarray* [ai_Ny];
  for (s=0; s<ai_Ny-1; s++) {
    atpo_MaxInd [s] = new intarray (ao_NumberMax(s));
    (*atpo_MaxInd [s])(0) = 0; // first max
    int l=1;
    for (int i=0; i<ai_Nx; i++) {
       if (ao_Mr1dData(s,i) < -EPS || ao_Mr1dData(s,i) > EPS)
	  (*atpo_MaxInd [s])(l++) = i;
    }
    (*atpo_MaxInd [s])(l++) = ai_Nx-1; // last max
  } 


  // fltarray of recons signal
  fltarray  ao_RecData (ai_Nx);

  // local mr1d transform of signal and recons signal
  MR_1D ao_Mr1dRecData (ai_Nx, se_Transform, "", ai_Ny);
  ao_Mr1dRecData.Border = I_MIRROR;  

  // trace
  if (Verbose) {
    cout << "Number of scales = " << ao_Mr1dData.nbr_scale() << endl;
    //cout << "CONSTRAINT TYPE = " << (int)se_Constraint << endl;
  }
  
  
  // raz of all coef of recons signal (except last scale)
  for (s=0; s<ai_Ny; s++)
    for (int i=0; i<ai_Nx; i++)
      ao_Mr1dRecData(s,i) = 0.0;
	
  // int last scale
  init_last_scale (ai_Ny, ai_Nx, ao_Mr1dData, ao_Mr1dRecData);

  //fits_write_fltarr ("InitLastScale", ao_Mr1dRecData.image());
  //fits_write_fltarr ("LastScaleAndMax", ao_Mr1dData.image());      

  // --------------------------------------------------- 
  // ---------------------------------------------------
  // iteration for recons
  for (int ai_Iter=0; ai_Iter<si_MaxIter; ai_Iter++){

    // --------------------------------------------------- 
    // ---------------------------------------------------
    interpolate (ai_Ny, ai_Nx, ao_NumberMax, (intarray**)atpo_MaxInd, ao_Mr1dData, ao_Mr1dRecData);

    // modify the modulus maximum of Mr1dRecDat : 
    // max (Mr1dRecData(s,i)) =  max (Mr1dData(s,i))
    init_max (ai_Ny, ai_Nx, ao_NumberMax, (intarray**)atpo_MaxInd, ao_Mr1dData, ao_Mr1dRecData);

    // int last scale
    init_last_scale (ai_Ny, ai_Nx, ao_Mr1dData, ao_Mr1dRecData);

    // --------------------------------------------------- 
    // ---------------------------------------------------
    // invers and foward transform of Mr1dRecDat
    ao_Mr1dRecData.recons (ao_RecData);
    if (Verbose) {
      //sprintf (atc_FileName, "fic_%d", ai_Iter);
      //fits_write_fltarr(atc_FileName, ao_RecData);
    }
    ao_Mr1dRecData.transform (ao_RecData);  
    if (Verbose) {
      //sprintf (atc_FileName, "mr1d_%d", ai_Iter);
      //fits_write_fltarr (atc_FileName, ao_Mr1dRecData.image());
    }

    if (Verbose) {
      cout << "number iteration :" << ai_Iter << endl;;
    }
  }
      
  // --------------------------------------------------- 
  // ---------------------------------------------------
  // modify the modulus maximum of Mr1dRecDat : 
  // max (Mr1dRecData(s,i)) =  max (Mr1dData(s,i))
  init_max (ai_Ny, ai_Nx, ao_NumberMax, (intarray**)atpo_MaxInd, ao_Mr1dData, ao_Mr1dRecData);

  // int last scale
  init_last_scale (ai_Ny, ai_Nx, ao_Mr1dData, ao_Mr1dRecData);

  // recons of data
  ao_Mr1dRecData.recons (ao_RecData);
  io_1d_write_data(stc_NameImagOut, ao_RecData);
  //fits_write_fltarr ("EndMr1d", ao_Mr1dRecData.image());
  
  exit(0);
} 

