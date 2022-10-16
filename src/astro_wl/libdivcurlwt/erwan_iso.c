#include <stdlib.h>
#include <stdio.h>
#include <math.h>


// void tfodi2d(float c[2][N][N],float ddn[2][N][N])
void div_isotrop_trans(fltarray & Tabc, fltarray & Tabd, fltarray & Tabddn)
{
  int N = Tabc.axis(1);
  float ***c;
  float ***ddn;
  float ***d;
  // get_ptr3d(Tabc, c);
  init_ptr3d(Tabc, c);
  put2_ptr3d(Tabc, c);
   
 float *k1 = new float[N];
 float *k2 = new float [N];
 Tabd.resize(N,N,2);
 // get_ptr3d(Tabd, d);
 init_ptr3d(Tabd, d);
   
 Tabddn.resize(N,N,2);
 init_ptr3d(Tabddn, ddn);
  // get_ptr3d(Tabddn, ddn);
 
 int j,j1,j2,i1,i2;
 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=c[0][i1][i2];
   d[1][i1][i2]=c[1][i1][i2];
  }
 } 
 
  int Nscale=1;
 for (j=N/2;j>0;j=j/2) {
 Nscale++;
  for (i2=0;i2<2*j;i2++) {
   for (i1=0;i1<2*j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=h[1][1]*k1[(2*i1-1+2*j)%(2*j)]+h[1][2]*k1[2*i1]+h[1][3]*k1[2*i1+1]+h[1][4]*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=g[1][1]*k1[(2*i1-1+2*j)%(2*j)]+g[1][2]*k1[2*i1]+g[1][3]*k1[2*i1+1]+g[1][4]*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=h[0][0]*k2[(2*i1-2+2*j)%(2*j)]+h[0][1]*k2[(2*i1-1+2*j)%(2*j)]+h[0][2]*k2[2*i1]+h[0][3]*k2[2*i1+1]+h[0][4]*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=g[0][2]*k2[2*i1]+g[0][3]*k2[2*i1+1]+g[0][4]*k2[(2*i1+2)%(2*j)];
   }
  }
  for (i1=0;i1<2*j;i1++) {
   for (i2=0;i2<2*j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=h[0][0]*k1[(2*i2-2+2*j)%(2*j)]+h[0][1]*k1[(2*i2-1+2*j)%(2*j)]+h[0][2]*k1[2*i2]+h[0][3]*k1[2*i2+1]+h[0][4]*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=g[0][2]*k1[2*i2]+g[0][3]*k1[2*i2+1]+g[0][4]*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=h[1][1]*k2[(2*i2-1+2*j)%(2*j)]+h[1][2]*k2[2*i2]+h[1][3]*k2[2*i2+1]+h[1][4]*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=g[1][1]*k2[(2*i2-1+2*j)%(2*j)]+g[1][2]*k2[2*i2]+g[1][3]*k2[2*i2+1]+g[1][4]*k2[(2*i2+2)%(2*j)];
   }
  }
 }  
  
	   
 // cout << "Nscale = " << Nscale << endl;
 for (j=N/2;j>0;j=j/2) {
  for (j1=j;j1>j/2;j1=j1/2) {
   for (j2=j;j2>j/2;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      // if (( i1 == 1) && ( i2 == 1)) cout << " j = " << j << " j1 = " << j1 << " j2 = " << j2 << " " << float(j1) / (j1*j1+j2*j2) << " " << float(j2) / (j1*j1+j2*j2) << endl;
      ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2] -  j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
      ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2] +  j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j;i2++) {
    ddn[0][i1][i2]=-d[1][i1][i2]/4.;
    ddn[1][i1][i2]=d[0][i1][i2]+d[1][i1][i2]/4.-d[1][i1][(i2-1+j)%(j)]/4.;
    ddn[0][i2][i1]=-d[0][i2][i1]/4.;
    ddn[1][i2][i1]=d[1][i2][i1]+d[0][i2][i1]/4.-d[0][(i2-1+j)%(j)][i1]/4.;
   }
  }
 }
 
 ddn[0][1][0]=d[1][1][0];
 ddn[1][1][0]=d[0][1][0];
 ddn[0][0][1]=d[0][0][1];
 ddn[1][0][1]=d[1][0][1];
 ddn[0][0][0]=d[0][0][0];
 ddn[1][0][0]=d[1][0][0];
 
 get2_ptr3d(Tabd,  d); 
 get2_ptr3d(Tabddn,  ddn);
}

/**************************************************************/

void div_isotrop_rec(fltarray & Tabddn, fltarray & Tabc)
/* void tfodiinv2d(float c[2][N][N],float ddn[2][N][N]) */
{
  int N = Tabddn.axis(1);
  float ***c;
  float ***ddn;
  float ***d;
  
  init_ptr3d(Tabddn, ddn);
  init_ptr3d(Tabddn, c);   
  init_ptr3d(Tabddn, d);
  
  put2_ptr3d(Tabddn, ddn);

 float *k1 = new float[N];
 float *k2 = new float [N];
 
 int j,j1,j2,i1,i2,l;
 float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
 float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 c[0][0][0]=ddn[0][0][0];
 c[1][0][0]=ddn[1][0][0];
 c[1][1][0]=ddn[0][1][0];
 c[0][1][0]=ddn[1][1][0];
 c[0][0][1]=ddn[0][0][1];
 c[1][0][1]=ddn[1][0][1];
 for (j=N/2;j>0;j=j/2) {
  for (j1=j;j1>j/2;j1=j1/2) {
   for (j2=j;j2>j/2;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
      c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j;i2++) {
    c[1][i1][i2]=-4.*ddn[0][i1][i2];
    c[0][i1][i2]=ddn[1][i1][i2]+ddn[0][i1][i2]-ddn[0][i1][(i2-1+j)%(j)];
    c[0][i2][i1]=-4.*ddn[0][i2][i1];
    c[1][i2][i1]=ddn[1][i2][i1]+ddn[0][i2][i1]-ddn[0][(i2-1+j)%(j)][i1];
   }
  }
 }
 
 for (i1=0;i1<N;i1++) {
  k1[i1]=0.;
  k2[i1]=0.;
 }
 
 for (j=1;j<N;j=2*j) {
  for (i2=0;i2<2*j;i2++) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     k1[(2*i1+l+2*j)%(2*j)]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][j+i1][i2]+k1[(2*i1+l+2*j)%(2*j)];
     k2[(2*i1+l+2*j)%(2*j)]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][j+i1][i2]+k2[(2*i1+l+2*j)%(2*j)];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=k1[i1];
    k1[i1]=0.;
    c[1][i1][i2]=k2[i1];
    k2[i1]=0.;
   }
  }
  
  for (i1=0;i1<2*j;i1++) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     k1[(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][i1][j+i2]+k1[(2*i2+l+2*j)%(2*j)];
     k2[(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][i1][j+i2]+k2[(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=k1[i2];
    k1[i2]=0.;
    c[1][i1][i2]=k2[i2];
    k2[i2]=0.;
   }
  }
 }
 get2_ptr3d(Tabc,  c);
}


/**************************************************************/
/* void tfogi2d(float c[2][N][N],float ddn[2][N][N]) */
void rot_isotrop_trans(fltarray & Tabc, fltarray & Tabddn)
{
  fltarray Tabd;
  int N = Tabc.axis(1);
  float ***c;
  float ***ddn;
  float ***d;
  // get_ptr3d(Tabc, c);
  init_ptr3d(Tabc, c);
  put2_ptr3d(Tabc, c);
   
 float *k1 = new float[N];
 float *k2 = new float [N];
 Tabd.resize(N,N,2);
 // get_ptr3d(Tabd, d);
 init_ptr3d(Tabd, d);
   
 Tabddn.resize(N,N,2);
 init_ptr3d(Tabddn, ddn);
  // get_ptr3d(Tabddn, ddn);
  
 int j,j1,j2,i1,i2;
 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=c[0][i1][i2];
   d[1][i1][i2]=c[1][i1][i2];
  }
 }
  
 for (j=N/2;j>0;j=j/2) {
  for (i2=0;((i2<2*j)&&(i2<N));i2++) {
   for (i1=0;i1<2*j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=h[0][0]*k1[(2*i1-2+2*j)%(2*j)]+h[0][1]*k1[(2*i1-1+2*j)%(2*j)]+h[0][2]*k1[2*i1]+h[0][3]*k1[2*i1+1]+h[0][4]*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=g[0][1]*k1[(2*i1-1+2*j)%(2*j)]+g[0][2]*k1[2*i1]+g[0][3]*k1[2*i1+1]+g[0][4]*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=h[1][0]*k2[(2*i1-2+2*j)%(2*j)]+h[1][1]*k2[(2*i1-1+2*j)%(2*j)]+h[1][2]*k2[2*i1]+h[1][3]*k2[2*i1+1]+h[1][4]*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=g[1][1]*k2[(2*i1-1+2*j)%(2*j)]+g[1][2]*k2[2*i1]+g[1][3]*k2[2*i1+1]+g[1][4]*k2[(2*i1+2)%(2*j)];
   }
  }
  for (i1=0;((i1<2*j)&&(i1<N));i1++) {
   for (i2=0;i2<2*j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=h[1][0]*k1[(2*i2-2+2*j)%(2*j)]+h[1][1]*k1[(2*i2-1+2*j)%(2*j)]+h[1][2]*k1[2*i2]+h[1][3]*k1[2*i2+1]+h[1][4]*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=g[1][1]*k1[(2*i2-1+2*j)%(2*j)]+g[1][2]*k1[2*i2]+g[1][3]*k1[2*i2+1]+g[1][4]*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=h[0][0]*k2[(2*i2-2+2*j)%(2*j)]+h[0][1]*k2[(2*i2-1+2*j)%(2*j)]+h[0][2]*k2[2*i2]+h[0][3]*k2[2*i2+1]+h[0][4]*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=g[0][1]*k2[(2*i2-1+2*j)%(2*j)]+g[0][2]*k2[2*i2]+g[0][3]*k2[2*i2+1]+g[0][4]*k2[(2*i2+2)%(2*j)];
   }
  }
 }
 
 for (j=N/2;j>0;j=j/2) {
  for (j1=j;j1>j/2;j1=j1/2) {
   for (j2=j;j2>j/2;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
      ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j;i2++) {
    ddn[0][i1][i2]=d[1][i1][i2]-d[0][i1][i2]/4.+d[0][i1][(i2-1+j)%(j)]/4.;
    ddn[1][i1][i2]=d[0][i1][i2]/4.;
    ddn[0][i2][i1]=d[0][i2][i1]-d[1][i2][i1]/4.+d[1][(i2-1+j)%(j)][i1]/4.;
    ddn[1][i2][i1]=d[1][i2][i1]/4.;
   }
  }
 }
 
 ddn[0][1][0]=d[1][1][0];
 ddn[1][1][0]=d[0][1][0];
 ddn[0][0][1]=d[0][0][1];
 ddn[1][0][1]=d[1][0][1];
 ddn[0][0][0]=d[0][0][0];
 ddn[1][0][0]=d[1][0][0];
 
  get2_ptr3d(Tabddn,  ddn);

}


/*************************************************************/

/* void tfogiinv2d(float c[2][N][N],float ddn[2][N][N]) */
void rot_isotrop_rec(fltarray & Tabddn, fltarray & Tabc)
/* void tfodiinv2d(float c[2][N][N],float ddn[2][N][N]) */
{
  int N = Tabddn.axis(1);
  float ***c;
  float ***ddn;
  float ***d;
  
  init_ptr3d(Tabddn, ddn);
  init_ptr3d(Tabddn, c);   
  init_ptr3d(Tabddn, d);
  
  put2_ptr3d(Tabddn, ddn);

 float *k1 = new float[N];
 float *k2 = new float [N];
 int j,j1,j2,i1,i2,l;
 float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
 float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 c[0][0][0]=ddn[0][0][0];
 c[1][0][0]=ddn[1][0][0];
 c[1][1][0]=ddn[0][1][0];
 c[0][1][0]=ddn[1][1][0];
 c[0][0][1]=ddn[0][0][1];
 c[1][0][1]=ddn[1][0][1];
 for (j=N/2;j>0;j=j/2) {
  for (j1=j;j1>j/2;j1=j1/2) {
   for (j2=j;j2>j/2;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
      c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j;i2++) {
    c[0][i1][i2]=4.*ddn[1][i1][i2];
    c[1][i1][i2]=ddn[0][i1][i2]+ddn[1][i1][i2]-ddn[1][i1][(i2-1+j)%(j)];
    c[0][i2][i1]=ddn[0][i2][i1]+ddn[1][i2][i1]-ddn[1][(i2-1+j)%(j)][i1];
    c[1][i2][i1]=4.*ddn[1][i2][i1];
   }
  }
 }
 
 for (i1=0;i1<N;i1++) {
  k1[i1]=0.;
  k2[i1]=0.;
 }
 
 for (j=1;j<N;j=2*j) {
  for (i2=0;((i2<2*j)&&(i2<N));i2++) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     k1[(2*i1+l+2*j)%(2*j)]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][j+i1][i2]+k1[(2*i1+l+2*j)%(2*j)];
     k2[(2*i1+l+2*j)%(2*j)]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][j+i1][i2]+k2[(2*i1+l+2*j)%(2*j)];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=k1[i1];
    k1[i1]=0.;
    c[1][i1][i2]=k2[i1];
    k2[i1]=0.;
   }
  }
  
  for (i1=0;((i1<2*j)&&(i1<N));i1++) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     k1[(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][i1][j+i2]+k1[(2*i2+l+2*j)%(2*j)];
     k2[(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][i1][j+i2]+k2[(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=k1[i2];
    k1[i2]=0.;
    c[1][i1][i2]=k2[i2];
    k2[i2]=0.;
   }
  }
 }
 
  get2_ptr3d(Tabc,  c);

}





