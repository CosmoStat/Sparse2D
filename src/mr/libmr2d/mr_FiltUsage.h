#ifndef _MR_FILTUSAGE_H
#define _MR_FILTUSAGE_H


extern char* OptArg;

//******************************/
// Option in filter prog
//******************************/


inline void poisson_usage()
{
    fprintf(OUTMAN, "         [-a ascii_file] or [-I image_file] \n");
    fprintf(OUTMAN, "             a & I options can't be used together, \n");
    fprintf(OUTMAN, "             and one must be set. \n");
}

inline void verbose () {
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose. \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}
  
inline void use_adjoint(){
    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "              Do not use adjoint operator for reconstruction. \n");
    fprintf(OUTMAN, "              Default is yes\n"); 
}

inline void corected_imag () {

    fprintf(OUTMAN, "         [-c image file name] \n");
    fprintf(OUTMAN, "              image used in filtering algo\n");
 
}
  
inline void border_usage (type_border e_Border) {
    fprintf(OUTMAN, "         [-b type_border = 1 to 4]\n");
    fprintf(OUTMAN, "              give the type of border {1:I_CONT|2:I_MIRROR|3:I_ZERO|4:I_PERIOD}\n");
    fprintf(OUTMAN, "              ( default=2:I_MIRROR ).\n");
}

inline void readint (int& pi_IntRead, char* ppc_msg) {
  if (sscanf(OptArg,"%d",&pi_IntRead) != 1) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit(-1);
  }
}

inline void readfloat (float& pf_FloatRead, char* ppc_msg) {
  if (sscanf(OptArg,"%f",&pf_FloatRead) != 1) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit(-1);
  }
}

//inline void readchar (char*& pt_CharRead, char* ppc_msg) {
//  if (sscanf(OptArg,"%s",&pt_CharRead) != 1) {
//    fprintf(OUTMAN, ppc_msg, OptArg);
//    exit(-1);
//  }
//}

inline Bool testborn (int pi_Min, int pi_Max, int pi_Val, char* ppc_msg) {
  if ((pi_Val<pi_Min) && (pi_Val>pi_Max)) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit (-1);
  }
  return (True);
}

inline Bool testborn (float pf_Min, float pf_Max, float pf_Val, char* ppc_msg) {
  if ((pf_Val<pf_Min) && (pf_Val>pf_Max)) {
    fprintf(OUTMAN, ppc_msg, OptArg);
    exit (-1);
  }
  return (True);
}









#endif
