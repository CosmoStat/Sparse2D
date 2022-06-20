/******************************************************************************
**                   Copyright (C) 2003 by CEA + Valiencoa observatory
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck and Vicent matinex
**
**    Date:  27/02/03
**    
**    File:  coxline.cc
**
*******************************************************************************
**
**    DESCRIPTION simulation of Cox line segment point process
**    ----------- 
**                 
******************************************************************************/
  
#include "Array.h"
#include "IM_IO.h"
#include "NR.h"
#include "FFTN.h"
#include "FFTN_3D.h"
#include "DefPoint.h"

#include <ctime>
extern int  OptInd;
extern char *OptArg;
char Name_Imag_Out[256]; /* output file name */

int WindowSize = 100;
float SegmentLength = 10;
float IntSegmentCenter = 0.001;
float IntPointSegment = 0.6;
int MaxSimu = 10000;
Bool Verbose=False;

/***************************************************************************/
  
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options out_catalogue \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-D Dimension]\n");
    fprintf(OUTMAN, "            Cube size dimension in the simulation.\n");
    fprintf(OUTMAN, "            Default is %d.\n",WindowSize);
    fprintf(OUTMAN, "         [-l SegmentLength]\n");
    fprintf(OUTMAN, "             Segment length. Default is %5.2f.\n", SegmentLength);
    fprintf(OUTMAN, "         [-c IntSegmentCenter]\n");
    fprintf(OUTMAN, "             Intensity of segment centers. Default is %f.\n", IntSegmentCenter);
    fprintf(OUTMAN, "         [-p IntPointSegment]\n");
    fprintf(OUTMAN, "             Intensities of points on the segments. Default is %f.\n", IntPointSegment);
    fprintf(OUTMAN, "         [-M MaxNpSimu]\n");
    fprintf(OUTMAN, "             Maximum number of points in the simulation. Default is %d.\n", MaxSimu);
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c; 
     /* get options */
    while ((c = GetOpt(argc,argv,"M:D:l:c:p:v")) != -1) 
    {
	switch (c) 
        { 
           case 'M': 
	        if (sscanf(OptArg,"%d",&MaxSimu) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of points parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'D': 
	        if (sscanf(OptArg,"%d",&WindowSize) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad window size parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'l': 
	        if (sscanf(OptArg,"%f",&SegmentLength) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad segment length: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'c':
	        if (sscanf(OptArg,"%f",&IntSegmentCenter) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad intensity parameters: %s\n", OptArg);
                    exit(-1);
                }
                break; 
           case 'p':  
	        if (sscanf(OptArg,"%f",&IntPointSegment) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad intensity parameters: %s\n", OptArg);
                    exit(-1);
                }
                break; 
	    case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    int k;
    float a,b,c,d,pa,pp,p,s,u,v,w,w1,w2,r1,r2,r3,xx,yy,zz;
    int i,n,na;
    float *x,*y,*z;
    // FILE *puntos;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);  
    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;   
        cout << "Cube size for the simulated data = " << WindowSize << endl;
	cout << "SegmentLength = " << SegmentLength    << endl; 
	cout << "IntSegmentCenter = " << IntSegmentCenter    << endl;  
	cout << "IntPointSegment = " << IntPointSegment    << endl; 
    }
        
    x = new float [MaxSimu];
    y = new float [MaxSimu];
    z = new float [MaxSimu];
    
//     if((puntos=fopen(Name_Imag_Out,"wb"))==NULL) 
//     {
//        fprintf(stderr,"Cannot open file %s:\n", Name_Imag_Out); 
//        exit(1); 
//     }
    
    srand48(time(NULL));
    a=b=c=WindowSize;	// window sides  
    d=SegmentLength;        // segment length
    pa=IntSegmentCenter;              // intensity of segment centres,
    pp=IntPointSegment;                // intensities of points on the segments  
	
	
    p=pa*(a+2*d)*(b+2*d)*(c+2*d);
    s=0.; na=0; 
    while((s-=log(drand48()))<=p) na++;
    for(n=0,i=1;i<=na;i++) 
    {
       u=drand48()*(a+2*d)-d; v=drand48()*(b+2*d)-d;
       w=drand48()*(c+2*d)-d;
       w1=drand48()*M_PI*2; w2=drand48(); w2=atan(sqrt(1-w2*w2)/w2);
       r1=sin(w2)*cos(w1); r2=sin(w2)*sin(w1); r3=cos(w2);
       if(drand48()<0.5) r3=-r3;
       s=0.;
       while((s-=log(drand48())/pp)<=d) 
       {
          xx=u+s*r1; 
	  if(xx<0 || xx>a) continue;
          yy=v+s*r2; 
	  if(yy<0 || yy>b) continue;
          zz=w+s*r3; 
	  if(zz<0 || zz>c) continue;
          x[n]=xx; 
          y[n]=yy; 
          z[n]=zz;
	  n++;
                // fprintf(puntos,"%8.4f %8.4f %8.4f\n",xx,yy,zz);
	}
     }
     
     // fclose(puntos);
     int Np = n;
     ArrayPoint AP(3, Np);
     AP.Pmin.x() = 0;
     AP.Pmin.y() = 0;
     AP.Pmin.z() = 0;
     AP.Pmax.x() = WindowSize;
     AP.Pmax.y() = WindowSize;
     AP.Pmax.z() = WindowSize;
     for(i=0; i< Np;i++) 
     {
        AP(i).x() = x[i];
	AP(i).y() = y[i];
	AP(i).z() = z[i];
     }
     AP.write(Name_Imag_Out);
    // printf("%5d\n",n);
    exit(0);
}
