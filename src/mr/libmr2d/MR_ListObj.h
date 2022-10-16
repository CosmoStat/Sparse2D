

#ifndef		_LISTOBJ2D_H_
#define	_LISTOBJ2D_H_

#include "IM_Math.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "IM_VisTool.h"
#include "MR_VisTree.h"
#include "IM_Graphics.h"
#include "MR_MVM.h"
 
class Object_2D {
     public:
        Object_2D()
        {
           SigmaX = 0.;    
           SigmaY = 0.;    
           PosX = 0.;      
           PosY = 0.;      
           Ra = 0.;
           Dec = 0.;
           Angle = 0.;
           Flux = 0.;
           Magnitude = 0.;
           ErrorRec = 0.;
           ValPixMax = 0.;
           DepX = DepY = 0;
           Image = NULL;
           SubObj = False;
           Perimeter = 0;
           Surface = 0;
           Morpho = 0.;
           Xmax = -1;
           Ymax=-1;
           Moment2XX = 0.;
           Moment2YY = 0.;
           Moment2XY = 0.;
           Mean = 0.;
           PosMaxCoef_X = 0;
           PosMaxCoef_Y = 0;
           ValMaxCoef = 0;
           SNR_ValMaxCoef = 0.;
           ErrorFlux = 0.;
	   FluxObjRec = 0.;
	   Extended = False;
	   Isotrop = True;
        };
        Bool Isotrop;           // True if the object is isotropic
        int NumObj;             // object number (starting from 1)
        int NumScale;           // detection scale number
        int NumObjScale;        // object number at the detection scale
        int PosMaxCoef_X;       // Position (X) of the wavelet coef maximum
        int PosMaxCoef_Y;       // Position (Y) of the wavelet coef maximum
        float ValMaxCoef;       // Value of the wavelet coef maximum
        float SNR_ValMaxCoef;   // SNR of the maximum wavelet coef
        float SNR_Obj;          // SNR of source
        Bool SubObj;            // True if the object is a sub-object
	int FatherScale;        // Father scale
	int FatherNumberInScale;// Father obect number inside the scale
        int Nlo, Nco;           // Image size which contains the object
        int Nli, Nci;           // Image size of the original image (input 
                                // image)
        double SigmaX, SigmaY;  // standard deviation in both directions
        double PosX, PosY;      // position of the center (barycenter)
                                // of the object
        double Angle;           // angle between the main axis of the object
                                // and the x-axis (degrees)
        double Flux;            // total flux
	double FluxObjRec;      // total flux obtained from the recons. object
        double ErrorFlux;       // Error on the total flux
        double Magnitude;       // magnitude
        double Ra, Dec;         // right ascension and declinaison
        double ValPixMax;       // maximum value of the object
        double ErrorRec;        // reconstruction error
	Bool Extended;          // True if the object is extended
        Ifloat *Image;          // pointer to the object image
        int DepX, DepY;         // Position of the object image in the
                                // input image
        int Perimeter;          // Object Perimeter
        int Surface;            // Object surface
        float Morpho;           // Morphological paramter
                                //   = Mean deviation of shape from sphericity
        int Xmax, Ymax;         // Position of the maximum;
        double Moment2XX, Moment2YY, Moment2XY; // second orders moments
        double Mean;            // Mean object value (flux / surface)
        void put_in_image (Ifloat &Ima_Big) { 
              if (Image != NULL) add_image (Ima_Big, *Image, DepX,DepY);};
        void get_box_obj (Ifloat &Ima_Big, Ifloat &ImaObj);
        void draw (Ifloat &Ima_Big, Bool DotFlag=False, float Val=0.);
        void simu (Ifloat &Ima_Big);
        ~Object_2D (){if (Image != NULL) delete Image;};
   };

double get_aperture_flux_obj(Ifloat &Data, Ifloat &Bgr, 
                               double  PosX,  double PosY, 
			       double SigmaX, double SigmaY, float Ksigma=3);
// Return the flux of an object at position PosX,PosY and of size
// SigmaX,SigmaY using the aperture photometry method
// Ksigma fixes the size of the box in which the flux is integrated
//        Box size =  Ksigma*MAX(SigmaX,SigmaY)
// Data = in: original data
// Bgr = in: estimated background


#define NBR_BGR  3
enum t_bgr {BGR_IMA, BGR_MR, BGR_VALUE};
inline const char * StringBgr (t_bgr type)
{
    switch (type)
    {
    case BGR_IMA: return ("Aperture photometry using a background image model"); break;
    case BGR_MR: return ("Aperture photometry using an estimated background image"); break;
    case BGR_VALUE:return ("Aperture photometry using a constant background"); break;
    default: return ("Error: bad background estimation method");
    }
}
	  
class ListObj2D {
   public:
   t_mvm VisionModel; // Type of Vision model
   Foret *F_obj;
   MR2D_Seg *MVM;
   ListObj2D();
   Bool Verbose;      // Verbose parameter
   int FirstUseScale; // First used scale for the detection
   int LastUseScale;  // Last used scale for the detection
   int Nb_TotalObj;  // Total number of objects
   int Nb_TotalSubObj;  // Total number of objects
   float FluxMult;   // Value to be multiplied to the flux
   char *NameImagIn; // input image name
   int Nl, Nc;       // input image size
   int NbSubObj;     // subobject numbers (for MVM detetion)
   Bool KeepIma;     // if KeepIma == True the small image of each object
                     // are kept
   type_objrec_method ReconsMethod;  // Object reconstruction method
   int Nb_iter_rec;    // Number of iterations per object reconstruction 
                       // for MVM reconstruction
   double ErrorRec;   // Error parameter for MVM reconstruction
   float MeanMorpho;  // mean morphological parameter
   float Signif;      // Percentage of pixels belonging to an object
   Object_2D *TabObj; // Object parameter table (from 0 to Nb_TotalObj-1)
   Ifloat Ima_SumObj; 
   Bool IsotropSeparation; // True if we want to separe isotropic objects from the others
   float AnitropicSeparationParam; // Ratio (def 1.5) SigmaX/SigmaY which we use for the 
                                   // discrimination between anisotropic object and iso.
   Ifloat Ima_SumObjAnisotopic; // Contains  the anisotropic objects.
   Ifloat Psf;        // Psf used for the reconstruction
   Ifloat IMGauss;    // Gaussian image 
   float ZoomPSF;     // if (ZoomPSF > 1) then the PSF is oversampled
   Bool  Deconv;      // if True, restore deconvolved objects
   float Fwhm;        // If Deconv == True and Fwhm > 0, Fwhm
                      // limits the resolution to achieve

   Bool AperPhot;    // if True, use an aperture photometry
   t_bgr BgrMethod;  // Background image estimation method
                     // used by the aperture photometry
   float BackgroundValue; // Background value in case of aperture photmetry
                          // with constant background
   Ifloat BgrModel;  // Background image model 
                     // used by the aperture photmetry in case of
		     // background estimation method BGR_IMA and BGR_MR	
   float KSigmaAper; // K-sigma = size box for aperture photometry
                     // Ksigma fixes the size of the box in which the flux 
                     // is integrated:
                     //          Box size =  Ksigma*Sigma
                     //  where Sigma is the second moment of the object
   int BGRSize;      // low resolution image for background calculation	  
		          
   void create_list(Ifloat & Data, MRNoiseModel & ModelData, 
                    MultiResol & MR_Data, 
                    Bool WriteAllObj = False, Bool InfoSubOb=False,
                    Bool WriteObjFullSize = False);

//    void create_list(MR2D_Seg & MVM, Ifloat & Data, MultiResol & MR_Data,
//                     MRNoiseModel & ModelData, 
//                     Bool WriteAllObj = False, Bool InfoSubOb=False);
// 
// 
// 
//    void create_list(Ifloat &Ima, Ifloat &ImaSegment, int Nb_TotalObj, Bool WriteAllObj = False);              

   void recons_obj(MultiResol& W, Ifloat &Im_rec, double & Erreur,
                   double & Flux, t_Gauss_2D & Gauss, 
		   float PosX=-1, float PosY=-1);
		   
   void init_background(Ifloat &Data);
   // background initialization from the data
   
   double get_aperture_flux_obj(Ifloat &Data,  
                   double  PosX,  double PosY, double SigmaX, double SigmaY);
   // Return the aperture photometry of an object a position PosX,PosY and
   // of size SigmaX,SigmaY
   
   ~ListObj2D ();
};

/*
void mr_write_obj_ascii(char *name_mes, Object_2D *TabObj, int NbObj);
void mr_write_obj_tex(char *name_tex, char *Name_Imag_In, char *name_mes);
*/
#define NBR_ORDER 3
enum t_order {O_Obj, O_Ra, O_SNR};
inline const char * StringOrder (t_order type)
{
    switch (type)
    {
    case  O_Obj: return ("tex table ordered by object number"); break;
    case  O_Ra: return ("tex table ordered by the Right Ascension"); break;
    case  O_SNR:return ("tex table ordered by object SNR"); break;
    default: return ("Error: bad order method");
    }
}


inline void obj_rec_methode_usage(type_objrec_method ReconsMethod)
{
    fprintf(OUTMAN, "          [-M object_reconstruction_method]\n");
    for (int i = 0; i < NBR_RECOBJ_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                         StringRecObj((type_objrec_method)i));
    fprintf(OUTMAN, "              default is: %s.\n",StringRecObj(ReconsMethod));
   
}
inline void tab_order_usage(t_order Order)
{
    fprintf(OUTMAN, "\n        [-C RADEC_Table_Order]\n");
    for (int i = 0; i < NBR_ORDER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringOrder((t_order )i));
    fprintf(OUTMAN, "              default is: %s.\n", StringOrder(Order));
}

void mr_write_obj_ascii(ListObj2D & Objects, char *Name_AsciiFile, 
                        Bool OptNumScale = False,  Bool WriteSubObj = False, Bool OnlyIsotrop=True);
void mr_write_obj_tex(ListObj2D & Objects, char *Name_TexFile, 
                      Bool OptNumScale = False, Bool WriteSubObj = False, Bool OnlyIsotrop=True);
void mr_write_obj_ps(ListObj2D & Objects, char *Name_PsFile, 
                      Bool OptNumScale= False, Bool WriteSubObj = False, Bool OnlyIsotrop=True);
void mr_write_radec_tex(ListObj2D & Objects, char *nom_tex, 
                         t_order Order=O_Obj, Bool OptNumScale= False,
			 Bool WriteSubObj = False, Bool OnlyIsotrop=True);
#endif
