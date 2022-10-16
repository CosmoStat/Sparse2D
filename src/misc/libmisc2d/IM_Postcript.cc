

#include "DefMath.h"
#include "IM_Graphics.h"
#include "SoftInfo.h"

#define PS_ADD_DATE 0
#if PS_ADD_DATE
#include <sys/time.h>
#include <sys/param.h>
extern "C" {
     int gettimeofday(struct timeval *tp, struct timezone *tzp);
}
#endif

int im_create_ps_file (Ifloat &Pict, char *File_Name)
{   
    extern softinfo Soft;

#if PS_ADD_DATE
    struct timeval tm_val;
    struct timezone tm_zon;
#endif

    int i,j,h1,h2,xphys,yphys;
    float x1,x2,new_max,new_min,new_max_min;
    FILE *fp_print;
    int print_width=256,print_height=256;
    int upp_threshold = 100, low_threshold=0;
    float maxi,mini;
    int Nl= Pict.nl();
    int Nc= Pict.nc();

    maxi = mini = Pict(0);
    for (i = 0; i<Nl*Nc; i++)
    {
        if (Pict(i) > maxi) maxi = Pict(i);
        if (Pict(i) < mini) mini = Pict(i);
    }
    fp_print = fopen(File_Name,"wt");
    if (fp_print == NULL) return(-1);
    xphys = 300-print_width/2;
    yphys = 420-print_height/2;

    new_max = (float) upp_threshold*(maxi-mini)/100. + mini;
    new_min = (float) low_threshold*(maxi-mini)/100. + mini;
    new_max_min = new_max-new_min;
#if PS_ADD_DATE
    gettimeofday(&tm_val,&tm_zon);
#endif

    fprintf(fp_print,"%%!PS-Adobe-1.0\n");
    fprintf(fp_print,"%%%%Creator: %s\n",Soft.banner());
    fprintf(fp_print,"%%%%Title: %s\n",File_Name);
#if PS_ADD_DATE
    fprintf(fp_print,"%%%%CreationDate: %s",ctime(&(tm_val.tv_sec)));
#endif
    fprintf(fp_print,"%%%%Pages: 1\n");
    fprintf(fp_print,"%%%%DocumentFonts: (atend)\n");
    fprintf(fp_print,"%%%%BoundingBox: %d %d %d %d\n",
      xphys,yphys,xphys+print_width,yphys+print_height);
    fprintf(fp_print,"%%%%EndComments\n");
    fprintf(fp_print,"%%begin(plot) --> For TeX include only\n");
    fprintf(fp_print,"/m {moveto} def  /l {lineto} def  /s {show} def\n");
    fprintf(fp_print,"/rasterimage <");

    for (i=0; i<Nl; i+=2)
    {
        for (j=Nc; j>0; j-=4)
        {
            x1 = Pict(i,j-1);
            if (x1 > new_max)  x1 = new_max;
            if (x1 < new_min)  x1 = new_min;
            x1 = 15.*(x1-new_min)/new_max_min;
            h1 = 15-(int) x1;
            x2 = Pict(i,j-3);
            if (x2 > new_max)  x2 = new_max;
            if (x2 < new_min)  x2 = new_min;
            x2 = 15.*(x2-new_min)/new_max_min;
            h2 = 15- (int) x2;
            fprintf(fp_print,"%x%x ",h1,h2);
        }
        if (i != (Nl-1))  fprintf(fp_print,"\n");
    }
    fprintf(fp_print,"> def\n");
    fprintf(fp_print,"%%%%EndProlog\n");
    fprintf(fp_print,"%%%%Page: 1\n");
    fprintf(fp_print,"newpath\n");
    fprintf(fp_print,"%d %d translate\n",xphys,yphys);
    fprintf(fp_print,"0 0 m  0 %d l  %d %d l  %d 0 l  0 0 l\n",
                   print_height,print_width,print_height,print_width);
    fprintf(fp_print,"1 setlinewidth\nstroke\n");
    fprintf(fp_print,"newpath\n");
    fprintf(fp_print,"%d  %d  scale\n",print_width,print_height);
    fprintf(fp_print,"%d  %d  4 [%d 0 0 %d 0 0]{rasterimage} image\n",
                    Nl/2,Nc/2,Nl/2,Nc/2);
    fprintf(fp_print,"stroke\n");
    fprintf(fp_print,"%%end(plot) --> For TeX include only\n");
    fprintf(fp_print,"showpage\n");
    fclose(fp_print);
    return(0);
}

/***************************************************************/
