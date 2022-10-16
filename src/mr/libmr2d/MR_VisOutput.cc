

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_VisTool.h"
#include "MR_VisTree.h"
#include "MR_ListObj.h"
#include "NR.h"

#include <fstream>

//----------------------------------------------------------
//	Init_ps
//----------------------------------------------------------

void Init_ps(ofstream& Fic,int ox,int oy,int lx,int ly)
{
	/* entete et variables */

Fic<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;     
Fic<<"/mm{72 mul 25.4 div}def"<<endl;     
Fic<<"/ox "<<ox<<" mm "<<"def"<<endl;
Fic<<"/oy "<<oy<<" mm "<<"def"<<endl;
Fic<<"/lx "<<lx<<" mm "<<"def"<<endl;
Fic<<"/ly "<<ly<<" mm "<<"def"<<endl;     
Fic<<"/d 1 mm def"<<endl;
Fic<<"0.5 setlinewidth"<<endl;    
Fic<<"/Times-Roman findfont"<<endl;     
Fic<<"8 scalefont"<<endl;
Fic<<"setfont"<<endl;

	/* procedure plus */     
Fic<<endl;
Fic<<"/plus"<<endl;
Fic<<"{"<<endl;    
Fic<<"d neg 0 rmoveto"<<endl;
Fic<<"d 2 mul 0 rlineto"<<endl;
Fic<<"d neg d neg rmoveto"<<endl;
Fic<<"0 d 2 mul rlineto"<<endl;
Fic<<"stroke"<<endl;
Fic<<"}def"<<endl;     

	/* procedure point */

Fic<<endl;
Fic<<"/point"<<endl;
Fic<<"{"<<endl;
Fic<<"moveto"<<endl;
Fic<<"gsave plus grestore"<<endl;
Fic<<"d 1 mul dup rmoveto"<<endl;
Fic<<"show"<<endl;
Fic<<"}def"<<endl;     

	/* cadre */

Fic<<endl;
Fic<<"ox oy moveto"<<endl;
Fic<<"0 ly rlineto"<<endl;
Fic<<"lx 0 rlineto"<<endl;
Fic<<"0 ly neg rlineto"<<endl;
Fic<<"lx neg 0 rlineto"<<endl;
Fic<<"closepath"<<endl;
Fic<<"stroke"<<endl;     
}

//----------------------------------------------------------
//	init_tex
//----------------------------------------------------------
void	init_tex (FILE *f, char * name)
{
	assert (f);
	fprintf (f, "\\documentstyle[11pt]{article}\n");
	fprintf (f, "\\setlength{\\hoffset}{-15mm}\n");
	fprintf (f, "\\setlength{\\voffset}{-20mm}\n");
	fprintf (f, "\\setlength{\\textwidth}{160mm}\n");
	fprintf (f, "\\setlength{\\textheight}{220mm}\n");
	fprintf (f, "\\title{Detections list}\n");
	if (name)
		fprintf (f, "\\author{Image : %s}\n", name);
	fprintf (f, "\n");
	fprintf (f, "\\begin{document}\n");
	fprintf (f, "\\maketitle\n");
}

//----------------------------------------------------------
//	term_tab_tex
//----------------------------------------------------------

void		term_tab_tex (FILE * f)
{
	fprintf (f, "\\hline\n");
	fprintf (f, "\\end{tabular}\n");
	fprintf (f, "\\end{table}\n");
	fprintf (f, "\\clearpage\n");
	fprintf (f, "\\newpage\n");
}

//----------------------------------------------------------
//	init_tab_rec_tex
//----------------------------------------------------------

void		init_tab_rec_tex (FILE * f)
{
	fprintf (f, "\\begin{table}[h]\n");
	fprintf (f, "\\begin{tabular}{||c||c|c|c|c|c|c|c||}\n");
	fprintf (f, "\\hline\n");
	fprintf (f, "Object & x & y & $\\sigma_x$ & $\\sigma_y$ & $\\theta$ & $I_{max} $ & Flux \\\\\n");
	fprintf (f, "\\hline\n");
	
}

//----------------------------------------------------------
//		rec_tex
//----------------------------------------------------------

void	rec_tex (FILE * f, char *name)
{
	FILE		* fstat;
	int		Num[2], Nb_ligne;
	float		Data[8], x, y;
	char		ligne[10*CHAINE];

 	fstat = fopen (name, "r");
 
	Nb_ligne = 0;
	init_tab_rec_tex (f);

	while (fgets (ligne, 10*CHAINE, fstat))
	{
//		fgets (ligne, 10*CHAINE, fstat);
		sscanf (ligne, "%d %d", &Num[0], &Num[1]);

		fgets (ligne, 10*CHAINE, fstat);
		sscanf (ligne, "%f %f %f %f ",
					&Data[0], &Data[1],
					&Data[2], &Data[3]);
		fgets (ligne, 10*CHAINE, fstat);
		sscanf (ligne, "%f %f %f %f",
					&Data[4], &Data[5],
					&Data[6], &Data[7]);
		Nb_ligne ++;
		if (Nb_ligne == 34)
		{
			term_tab_tex (f);
			init_tab_rec_tex (f);
			Nb_ligne=0;
		}

		x = Data[0];
		y = Data[1];
		fprintf (f, "\\hline\n");
		fprintf (f, "%d-%d & %.2f & %.2f & ", Num[0], Num[1], x, y);
		fprintf (f, "%.2f & %.2f & %.1f & %.2f & %.2f  \\\\\n", 
			Data[2], Data[3], Data[4], Data[5], Data[6]);
	}
	term_tab_tex (f);
	fprintf (f, "\\end{document}\n");
	fclose(fstat);
}


/**********************************************************************/

void mr_write_obj_tex(char *nom_tex, char *Name_Imag_In, char *nom_mes)
{
   FILE *Fic_tex;

   Fic_tex = fopen(nom_tex, "w");
   init_tex (Fic_tex, Name_Imag_In);
   rec_tex (Fic_tex, nom_mes);
   fclose(Fic_tex);
}

/**********************************************************************/

void mr_write_obj_ascii(ListObj2D & Objects, char *nom_mes, Bool OptNumScale, Bool WriteSubObj, Bool OnlyIsotrop)
{
   FILE *Fic_mes;
   int i;
   int NbObj = Objects.Nb_TotalObj;
   Object_2D *TabObj = Objects.TabObj;
  
   Fic_mes = fopen(nom_mes, "w");
   if(!Fic_mes)
   {
      cerr << "Error: cannot open file "<< nom_mes <<endl;
      exit(1);
   }
   for (i=0 ; i < NbObj; i++)
   {
      if ((Objects.IsotropSeparation == False) ||
           ((TabObj[i].Isotrop == True) && (OnlyIsotrop == True)) ||
	   ((TabObj[i].Isotrop == False) && (OnlyIsotrop == False)))
      {   
        if (((TabObj[i].SubObj == False) && (WriteSubObj == False)) ||
          ((TabObj[i].SubObj == True) && (WriteSubObj == True)))
        {
          // write to file
          if (OptNumScale == True) 
           fprintf(Fic_mes, "%d %d\n", TabObj[i].NumScale, TabObj[i].NumObjScale);
          else 
           fprintf(Fic_mes, "%d\n", TabObj[i].NumObj);

          fprintf(Fic_mes, "%f %f %f %f\n",TabObj[i].PosX, TabObj[i].PosY, 
                                       TabObj[i].SigmaX, TabObj[i].SigmaY);
          fprintf(Fic_mes, " %f %f", TabObj[i].Angle, TabObj[i].ValPixMax);
          fprintf(Fic_mes, " %f %f\n",  TabObj[i].Flux, TabObj[i].Magnitude); 
          fprintf(Fic_mes, " %f %f %f\n",  TabObj[i].ErrorFlux,
                                  TabObj[i].SNR_ValMaxCoef, TabObj[i].SNR_Obj);
          fprintf(Fic_mes, " %d %d %d\n", TabObj[i].PosMaxCoef_X, TabObj[i].PosMaxCoef_Y, TabObj[i].Surface);
       }
     }
   }
   fclose(Fic_mes);
}

/**********************************************************************/

void mr_write_obj_tex(ListObj2D & Objects, char *nom_tex, Bool OptNumScale, Bool WriteSubObj, Bool OnlyIsotrop)
{
   FILE *Fic_tex;
   char *Name_Imag_In = Objects.NameImagIn;
   int i,Nb_ligne=0;
   Object_2D *TabObj = Objects.TabObj;
   int NbObj = Objects.Nb_TotalObj;
   float FluxMult = Objects.FluxMult;
   double Flux;
   
   Fic_tex = fopen(nom_tex, "w");
   init_tex (Fic_tex, Name_Imag_In);
   fprintf (Fic_tex, "\\begin{table}[h]\n");
   fprintf (Fic_tex, "\\begin{tabular}{||c||c|c|c|c|c|c|c|c||}\n");
   fprintf (Fic_tex, "\\hline\n");
   fprintf (Fic_tex,
"Object & x & y & $\\sigma_x$ & $\\sigma_y$ & $Angle$ & $I_{max} $ & Flux & SNR\\_WaveCoef \\\\\n");
   fprintf (Fic_tex, "\\hline\n");

   for (i=0; i < NbObj; i++)
   {
      if ((Objects.IsotropSeparation == False) ||
           ((TabObj[i].Isotrop == True) && (OnlyIsotrop == True)) ||
	   ((TabObj[i].Isotrop == False) && (OnlyIsotrop == False)))
     {
        if (((TabObj[i].SubObj == False) && (WriteSubObj == False)) ||
          ((TabObj[i].SubObj == True) && (WriteSubObj == True)))
       {
         fprintf (Fic_tex, "\\hline\n");
         if (OptNumScale == True)
               fprintf (Fic_tex, "%d-%d & %.2f & %.2f & ", 
                                  TabObj[i].NumScale, TabObj[i].NumObjScale,
                                  TabObj[i].PosX, TabObj[i].PosY);
         else fprintf (Fic_tex, "%d & %.2f & %.2f & ", 
                  TabObj[i].NumObj, TabObj[i].PosX, TabObj[i].PosY);
         Flux = TabObj[i].Flux * FluxMult;
         fprintf (Fic_tex, "%.2f & %.2f & %.1f & %.2f & %.2f & %.2f \\\\\n", 
			TabObj[i].SigmaX, TabObj[i].SigmaY, 
			TabObj[i].Angle, TabObj[i].ValPixMax, Flux,
			TabObj[i].SNR_ValMaxCoef);
         Nb_ligne ++;
         if (Nb_ligne == 34)
         {
            term_tab_tex (Fic_tex);
            fprintf (Fic_tex, "\\begin{table}[h]\n");
            fprintf (Fic_tex, "\\begin{tabular}{||c||c|c|c|c|c|c|c|c||}\n");
            fprintf (Fic_tex, "\\hline\n");
            fprintf (Fic_tex,
    "Object & x & y & $\\sigma_x$ & $\\sigma_y$ & $Angle$ & $I_{max} $ & Flux & SNR\\_WaveCoef \\\\\n");
   	    fprintf (Fic_tex, "\\hline\n");
            Nb_ligne=0;
         }
       }
     }
   }
   term_tab_tex (Fic_tex);
   fprintf (Fic_tex, "\\end{document}\n");
   fclose(Fic_tex);
}

/**********************************************************************/

void mr_write_radec_tex(ListObj2D & Objects, char *nom_tex, 
                        t_order Order, Bool OptNumScale, Bool WriteSubObj, Bool OnlyIsotrop)
{
   FILE *Fic_tex;
   char *Name_Imag_In = Objects.NameImagIn;
   int i,j,Nb_ligne=0;
   Object_2D *TabObj = Objects.TabObj;
   int NbObj = Objects.Nb_TotalObj;
   float *Tab = new float [NbObj+1];
   unsigned long *TabInd = new unsigned long[NbObj+1];
   float FluxMult = Objects.FluxMult;
   double Flux;
   
   Fic_tex = fopen(nom_tex, "w");
   init_tex (Fic_tex, Name_Imag_In);
   fprintf (Fic_tex, "\\begin{table}[h]\n");
   fprintf (Fic_tex, "\\begin{tabular}{||c||c|c|c|c|c|c||}\n");
   fprintf (Fic_tex, "\\hline\n");
   if (WriteSubObj == False)
        fprintf (Fic_tex, "Object & Ra & Dec  & Ra (H M S) & Dec (D M S) & Flux & SNR\\_WaveCoef \\\\\n");
   else fprintf (Fic_tex, "Sub-Object & Ra & Dec  & Ra (H M S) & Dec (D M S) & Flux & SNR\\_WaveCoef \\\\\n");
 
   fprintf (Fic_tex, "\\hline\n");

   Tab[0] = 0.;
   
   
   if (Order == O_Ra )
     for (i=0; i < NbObj; i++) Tab[i+1] = TabObj[i].Ra;
   else if (Order == O_SNR )
     for (i=0; i < NbObj; i++) Tab[i+1] = (float) TabObj[i].SNR_ValMaxCoef;
     
   if (Order != O_Obj) indexx ((unsigned long) NbObj, Tab, TabInd);
   
   for (j=1; j <= NbObj; j++)
   {
      double Ra, Dec, xsec, xsc;
      int ihr, imin, ideg, imn;
      
      if (Order == O_Obj) i=j;
      else if (Order == O_SNR) i = TabInd[NbObj-j+1]-1;
      else i = TabInd[j]-1;
      
      if ((Objects.IsotropSeparation == False) ||
           ((TabObj[i].Isotrop == True) && (OnlyIsotrop == True)) ||
	   ((TabObj[i].Isotrop == False) && (OnlyIsotrop == False)))
     {
     
       if (((TabObj[i].SubObj == False) && (WriteSubObj == False)) ||
          ((TabObj[i].SubObj == True) && (WriteSubObj == True)))     
       {       
        Ra= TabObj[i].Ra;
        Dec= TabObj[i].Dec;
        radec(Ra, Dec, ihr, imin, xsec, ideg, imn, xsc);
        Flux = TabObj[i].Flux * FluxMult;
        fprintf (Fic_tex, "\\hline\n");
        if (OptNumScale == True)
               fprintf (Fic_tex, "%d-%d & ", 
                                  TabObj[i].NumScale, TabObj[i].NumObjScale);
        else fprintf (Fic_tex, "%d &  ", TabObj[i].NumObj);
        fprintf (Fic_tex, "%.5f & %.5f & %2d  %2d  %2.2f & %2d  %2d  %2.2f &  %.2f & %g \\\\\n", 
			TabObj[i].Ra, TabObj[i].Dec, 
			ihr, imin, xsec, ideg, imn, xsc,
			Flux, TabObj[i].SNR_ValMaxCoef);
			
        Nb_ligne ++;
        if (Nb_ligne == 34)
        {
         term_tab_tex (Fic_tex);
         fprintf (Fic_tex, "\\begin{table}[h]\n");
         fprintf (Fic_tex, "\\begin{tabular}{||c||c|c|c|c|c|c||}\n");
         fprintf (Fic_tex, "\\hline\n");
         fprintf (Fic_tex, "Object & Ra & Dec & Ra (H M S) & Dec (D M S) & Flux &  SNR\\_WaveCoef \\\\\n");
	 fprintf (Fic_tex, "\\hline\n");
         Nb_ligne=0;
        }
      }
    }
   }
   term_tab_tex (Fic_tex);
   fprintf (Fic_tex, "\\end{document}\n");
   fclose(Fic_tex);
   delete Tab;
   delete TabInd;
}

/**********************************************************************/

void mr_write_obj_ps(ListObj2D & Objects, char *nom_psobj, Bool OptNumScale, Bool WriteSubObj, Bool OnlyIsotrop)
{
  // poscript output
  Object_2D *TabObj = Objects.TabObj;
  char *Name_Imag_In = Objects.NameImagIn;
  int Nl = Objects.Nl;;
  int Nc = Objects.Nc;
  int Nb_obj = Objects.Nb_TotalObj;
  int x,y,xmax,ymax;
  int Num_obj;
  ofstream Fic_obj;
  int Depx = OFFSET_IMA_PS_X;  // in millimeter
  int Depy = OFFSET_IMA_PS_Y;
  int Sizex = SIZE_IMA_PS_X;
  int Sizey = SIZE_IMA_PS_Y;
  float Ratio = (float) Nl / (float) Nc;

  if (Nl > Nc) Sizex = (int) ( (float) Sizex / Ratio);
  else if (Nc > Nl) Sizey = (int) ( (float) Sizey * Ratio);
  
  Fic_obj.open(nom_psobj);
  if(!Fic_obj)
  {
      cerr<<"Error: Unable to open "<< nom_psobj << endl;
      exit(1);
  }

   // INIT FICHIER PS
   Init_ps(Fic_obj,  Depx, Depy, Sizex, Sizey);

   // ECRITURE COORDONNES OBJETS
   for (int i = 0; i < Nb_obj;i++)
   {
      if ((Objects.IsotropSeparation == False) ||
           ((TabObj[i].Isotrop == True) && (OnlyIsotrop == True)) ||
	   ((TabObj[i].Isotrop == False) && (OnlyIsotrop == False)))
     {
     
      if (((TabObj[i].SubObj == False) && (WriteSubObj == False)) ||
          ((TabObj[i].SubObj == True) && (WriteSubObj == True)))     
         {
	   Num_obj = TabObj[i].NumObj;
  	   xmax = (int) TabObj[i].PosX;
	   ymax = (int) TabObj[i].PosY;

	   x = xmax * Sizex / (Nc-1) + Depx;
	   y = ymax * Sizey / (Nl-1) + Depy;

           if (OptNumScale == False)
	       Fic_obj<<"("<<Num_obj<<") "<<x<<" mm ";
	   else
	      Fic_obj<<"("<< TabObj[i].NumScale<<"-"<< TabObj[i].NumObjScale<<") "<<x<<" mm ";
	   Fic_obj<<y<<" mm point"<<endl;
	}
     }
   }

   // ECRITURE TITRE 

   Fic_obj<<"/Times-Roman findfont"<<endl;     
   Fic_obj<<"18 scalefont"<<endl;
   Fic_obj<<"setfont"<<endl;

   Fic_obj<<"10 mm 270 mm moveto"<<endl;
   Fic_obj<<"("<< Name_Imag_In <<" (objects)   "<<Nb_obj<<")"<<" show"<<endl;

   // FERMETURE FICHIERS PS

   Fic_obj<<"showpage"<<endl;
   Fic_obj.close();
}
