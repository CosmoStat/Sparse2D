#ifndef _PCA_USAGE_H
#define _PCA_USAGE_H

#include "MR1D_Obj.h"
#include "MR_PCA.h"

extern char* OptArg;

//******************************/
// Option in wk_...
//******************************/

inline void wk_noise_usage () {
   fprintf(OUTMAN, "         [-m type_of_noise]\n");
   for (int i = 0; i < NBR_NOISE-1; i++)              // !!! no posisson noise
       fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
       fprintf(OUTMAN, "             Default is Gaussian noise\n");
}

inline void noise_select () {
    fprintf(OUTMAN ,"         [-M]\n");
    fprintf(OUTMAN, "             Noise Model is used.\n");
   	fprintf(OUTMAN, "               Default is no.\n");
}

inline void write_eigen_on_disk () {
    fprintf(OUTMAN ,"         [-w]\n");
    fprintf(OUTMAN, "             Write to the disk the eigen vect images.\n");
    fprintf(OUTMAN, "               file names are: wk_ev_x \n");
    fprintf(OUTMAN, "               where x is the eigen vector number\n"); 
	fprintf(OUTMAN, "               Default is no.\n");
}

inline void nbr_iteration (int NbrIter) {
    fprintf(OUTMAN ,"         [-i NbrIteration]\n");
    fprintf(OUTMAN, "             Number of iterations.\n");
    fprintf(OUTMAN, "             Default is %d.\n", NbrIter);
}

inline void subtract_mean () {
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Do not subtract the mean. \n");
	fprintf(OUTMAN, "               Default is no.\n");
}

inline void nb_eigen_used () {
    fprintf(OUTMAN, "         [-F NbrEigenVect]\n");
    fprintf(OUTMAN, "             Number of eigen vector used for the reconstruction. \n");
    fprintf(OUTMAN, "             Default is set to the number of images\n");   
}

inline void nb_eigen_1_used () {
    fprintf(OUTMAN, "         [-F NbrEigenVect]\n");
    fprintf(OUTMAN, "             Number of eigen vector used for the reconstruction. \n");
    fprintf(OUTMAN, "             Default is 1.\n");   
}

inline void vector_not_used () {
    fprintf(OUTMAN, "         [-K EigenVect_Number]\n");
    fprintf(OUTMAN, "             Eigen vector number which will not be used for the reconstruction. \n"); 
}

inline void verbose () {
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose. \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}
  
inline void write_support_mr2d () {
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "              Write support Mr2d. \n");
    fprintf(OUTMAN, "              file neame is mr2d_Support_<n Imag>.mr \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void write_support_mr1d () {
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "              Write support Mr1d. \n");
    fprintf(OUTMAN, "              file neame is mr1d_Support_<n Spectre>.mr \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void write_support_pca () {
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "              Write support PCA. \n");
    fprintf(OUTMAN, "              file neame is pca_Support_<n Imag>.mr \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void write_support_wth () {
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "              Write support WTH. \n");
    fprintf(OUTMAN, "              file neame is wth_Support_<n Imag>.mr \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void dilate_support_pca () {
    fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "              Dilate support PCA. \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void destroy_rings_pca () {
    fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "              Destroy rings in PCA support. \n");
    fprintf(OUTMAN, "              Default is no\n"); 
}

inline void destroy_isol_pix () {
    fprintf(OUTMAN, "         [-u]\n");
    fprintf(OUTMAN, "              Destroy isol pixel in Mr2d support. \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void nsigma_pca  (float def_sigma_pca) {
    fprintf(OUTMAN, "         [-S nsigma]\n");
    fprintf(OUTMAN, "             Thresholding WT-KLT coeff at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             Default is %2.0f.\n", def_sigma_pca);
}

inline void nsigma_wth  (float def_sigma_wth) {
    fprintf(OUTMAN, "         [-S nsigma]\n");
    fprintf(OUTMAN, "             Thresholding wth at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             Default is %2.0f.\n", def_sigma_wth);
}

inline void nsigma_mr2d (float def_sigma_mr2d) {
    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "             Significant wavelet coeff larger than nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             Default is %2.0f.\n", def_sigma_mr2d);
}

inline void nsigma_mr1d (float def_sigma_mr1d) {
    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "             Thresholding mr1d at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             Default is %2.0f.\n", def_sigma_mr1d);
}

inline void positiv_imag_constraint () {
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "              Positivity constraint. \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void max_imag_constraint () {
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "              Maximun image constraint (255). \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void positiv_spectre_constraint () {
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "              Positivity constraint. \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void transform_usage(type_trans_1d Transform) {
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (int i=0; i<NBR_TRANS_1D; i++)
       fprintf(OUTMAN, "              %d: %s \n",i+1,
                       StringTransf1D((type_trans_1d)i));
    fprintf(OUTMAN, "             Default is %s\n",  
                    StringTransf1D((type_trans_1d)Transform));
}

inline void wk_conv_param (float pf_CvgParam) {
    fprintf(OUTMAN, "         [-C CvgParam]\n");
    fprintf(OUTMAN, "              Convergence parameter \n");
    fprintf(OUTMAN, "              Default is %f\n", pf_CvgParam);
}

inline void wk_regul_param (float pf_RegulVal) {
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              Default is %f\n", pf_RegulVal);
}

inline void wk_type_opt () {
    fprintf(OUTMAN, "         [-U Type_of_Regularization]\n");
    fprintf(OUTMAN, "              1: Use a fixed user Alpha value.\n");
    //fprintf(OUTMAN, "              2: Estimate the optimal Alpha.\n");
    fprintf(OUTMAN, "              2: Estimate one  Alpha value per band.\n");
    fprintf(OUTMAN, "              Default is 1.\n");
}

inline void wk_model_edge () {
    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "              Use the edge as a-priori model information.\n");
    fprintf(OUTMAN, "              Default is no.\n");
}

inline void wk_data_snr () {
    fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "              Alpha is modified using the data SNR.\n");
    fprintf(OUTMAN, "              Default is no.\n");
}

inline void wk_rms_map () {
    fprintf(OUTMAN, "         [-M RMS_Map_File_Name]\n");
    fprintf(OUTMAN, "              RMS Map.  If this Option is set, \n");
    fprintf(OUTMAN, "              the noise model is automatically fixed to:\n");
    fprintf(OUTMAN, "                 %s\n", StringNoise(NOISE_NON_UNI_ADD));
}

inline void wk_size_block_sig_clip (int SizeBlock) {
    fprintf(OUTMAN, "         [-S SizeBlock]\n");
    fprintf(OUTMAN, "             Size of the  blocks used for local variance estimation.\n");
    fprintf(OUTMAN, "             Default is %d \n", SizeBlock);
}

inline void wk_nb_iter_sig_clip (int NiterClip) {
    fprintf(OUTMAN, "         [-N NiterSigmaClip]\n");
    fprintf(OUTMAN, "             iteration number used for local variance estimation.\n");
    fprintf(OUTMAN, "             Default is %d \n", NiterClip);
}

inline void wk_nb_iter_entrop (int NbrIter) {
    fprintf(OUTMAN, "         [-i MaxIter]\n");
    fprintf(OUTMAN, "              Maximum number of iterations.\n");
    fprintf(OUTMAN, "              Default is %d.\n", NbrIter);  
}
  
inline void correl_compute (CorrelComputeType CorrelComp) {
    fprintf(OUTMAN, "         [-x CorrelMat_Method]\n");
    fprintf(OUTMAN, "              0 : %s \n", CorCompTransform(E_LOCAL));
    fprintf(OUTMAN, "              1 : %s \n", CorCompTransform(E_GLOBAL));
    fprintf(OUTMAN, "              2 : %s \n", CorCompTransform(E_GLOBAL_WITHOUT_LAST));   
    fprintf(OUTMAN, "              3 : %s \n", CorCompTransform(E_IMPORT));    
    fprintf(OUTMAN, "              Default is %s.\n",CorCompTransform(CorrelComp));  
}
  
inline void correl_noise (CorrelNoiseType CorrelNoise) {
    fprintf(OUTMAN, "         [-y Noise correlation type]\n");
    fprintf(OUTMAN, "              0 : %s \n", CorNoiseTransform(E_WITHOUT));
    fprintf(OUTMAN, "              1 : %s \n", CorNoiseTransform(E_THRESHOLD));
    fprintf(OUTMAN, "              2 : %s \n", CorNoiseTransform(E_LINEAR));
    //fprintf(OUTMAN, "              3 : %s \n", CorNoiseTransform(E_PROBA));       
    fprintf(OUTMAN, "              Default is %s.\n",CorNoiseTransform(CorrelNoise));  
}

inline void import_correl_file () {
    fprintf(OUTMAN, "         [-O Input_CorrelMat_FileName]\n");
    fprintf(OUTMAN, "              Used only with -x3 option \n");
    fprintf(OUTMAN, "              Default is no. \n");
}

inline void write_correl_matrix () {
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "              Write the correlation matrix on the disk. \n");
    fprintf(OUTMAN, "              File name is Correl_Matrix. \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void normalize_correl_matrix () {
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "              Use the covariance matrix instead of the. \n");
    fprintf(OUTMAN, "              correlation matrix. \n");
    fprintf(OUTMAN, "              Default is no.\n"); 
}

inline void wk_filter (type_sb_filter SB_Filter) {
    fprintf(OUTMAN, "         [-T Type of filter]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
       fprintf(OUTMAN, "              %d: %s \n",i, StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "              Default is %s\n\n", StringSBFilter((type_sb_filter) SB_Filter));	       
    fprintf(OUTMAN, " \n");
    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "              Use a L2 normalization. Default is L1.\n");
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











#endif
