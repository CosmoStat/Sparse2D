// tfodanp2d
void div_ani_trans(fltarray & Tabc, fltarray & Tabd, fltarray & Tabddn)
{
  int N = Tabc.axis(1);
  float ***c;
  float ***ddn;
  float ***d;
  set_ptr3d(Tabc, c);
   
   cout << " div_ani_trans .... " << endl;
   Tabc.info("in");
   
 float *k1 = new float[N];
 float *k2 = new float [N];
 Tabd.resize(N,N,2);
 set_ptr3d(Tabd, d);
   
 Tabddn.resize(N,N,2);
 set_ptr3d(Tabddn, ddn);
 
 int j,j1,j2,i1,i2;
 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   k1[i1]=c[0][i1][i2];
   k2[i1]=c[1][i1][i2];
  }
  for (j=N/2;j>0;j=j/2) {
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=-1./4.*k1[(2*i1-1+2*j)%(2*j)]+3./4.*k1[2*i1]+3./4.*k1[2*i1+1]-1./4.*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=1./8.*k1[(2*i1-1+2*j)%(2*j)]-3./8.*k1[2*i1]+3./8.*k1[2*i1+1]-1./8.*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=-1./8.*k2[(2*i1-2+2*j)%(2*j)]+1./4.*k2[(2*i1-1+2*j)%(2*j)]+3./4.*k2[2*i1]+1./4.*k2[2*i1+1]-1./8.*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=-1./4.*k2[2*i1]+1./2.*k2[2*i1+1]-1./4.*k2[(2*i1+2)%(2*j)];
   }
   for (i1=0;i1<j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
  }
 }
 
 
 for (i1=0;i1<N;i1++) {
  for (i2=0;i2<N;i2++) {
   k1[i2]=d[0][i1][i2];
   k2[i2]=d[1][i1][i2];
  }
  for (j=N/2;j>0;j=j/2) {
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=-1./8.*k1[(2*i2-2+2*j)%(2*j)]+1./4.*k1[(2*i2-1+2*j)%(2*j)]+3./4.*k1[2*i2]+1./4.*k1[2*i2+1]-1./8.*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=-1./4.*k1[2*i2]+1./2.*k1[2*i2+1]-1./4.*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=-1./4.*k2[(2*i2-1+2*j)%(2*j)]+3./4.*k2[2*i2]+3./4.*k2[2*i2+1]-1./4.*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=1./8.*k2[(2*i2-1+2*j)%(2*j)]-3./8.*k2[2*i2]+3./8.*k2[2*i2+1]-1./8.*k2[(2*i2+2)%(2*j)];
   }
   for (i2=0;i2<j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
  }
 } 
 fits_write_fltarr("adiv.fits", Tabd);
 
 for (j1=N/2;j1>0;j1=j1/2) {
  for (j2=N/2;j2>0;j2=j2/2) {
   for (i1=0;i1<j1;i1++) {
    for (i2=0;i2<j2;i2++) {
     ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
    }
   }
  }
 }
 
  fits_write_fltarr("adivbe.fits", Tabddn);

 for (i1=1;i1<N;i1++) {
  ddn[1][i1][0]=d[1][i1][0];
  ddn[0][i1][0]=d[0][i1][0];
  ddn[1][0][i1]=d[0][0][i1];
  ddn[0][0][i1]=d[1][0][i1];
 }
 ddn[0][0][0]=d[0][0][0];
 ddn[1][0][0]=d[1][0][0];
 
 fits_write_fltarr("adive.fits", Tabddn);
}

/*********************************************************************/

void div_ani_trans(fltarray & Tabc, fltarray & Tabddn)
{
   fltarray Tabd;
   div_ani_trans(Tabc, Tabd, Tabddn);
}

/*********************************************************************/

void div_ani_rec(fltarray & Tabddn, fltarray & Tabc)
{
 
  int N = Tabddn.axis(1);
  float ***c;
  float ***cw;
  float ***ddn;
    
  Tabc.resize(N,N,2);
  Tabddn.resize(N,N,2);
  set_ptr3d(Tabddn, ddn);
  
  init_ptr3d(Tabddn, cw);   
  set_ptr3d(Tabc, c);
   
   int j,j1,j2,i1,i2,l;
   float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
   float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 
 c[0][0][0]=ddn[0][0][0];
 c[1][0][0]=ddn[1][0][0];
 for (i1=1;i1<N;i1++) {
  c[0][i1][0]=ddn[0][i1][0];
  c[1][i1][0]=ddn[1][i1][0];
  c[1][0][i1]=ddn[0][0][i1];
  c[0][0][i1]=ddn[1][0][i1];
 }
 for (j1=N/2;j1>0;j1=j1/2) {
  for (j2=N/2;j2>0;j2=j2/2) {
   for (i1=0;i1<j1;i1++) {
    for (i2=0;i2<j2;i2++) {
     c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
     c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
    }
   }
  }
 }
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   cw[0][i1][i2]=0.;
   cw[1][i1][i2]=0.;
  }
 }
 for (i1=0;i1<N;i1++) {
  for (j=1;j<N;j=2*j) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     cw[0][i1][(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][i1][j+i2]+cw[0][i1][(2*i2+l+2*j)%(2*j)];
     cw[1][i1][(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][i1][j+i2]+cw[1][i1][(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=cw[0][i1][i2];
    cw[0][i1][i2]=0.;
    c[1][i1][i2]=cw[1][i1][i2];
    cw[1][i1][i2]=0.;
   }
  }
 }
 for (i2=0;i2<N;i2++) {
  for (j=1;j<N;j=2*j) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     cw[0][(2*i1+l+2*j)%(2*j)][i2]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][j+i1][i2]+cw[0][(2*i1+l+2*j)%(2*j)][i2];
     cw[1][(2*i1+l+2*j)%(2*j)][i2]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][j+i1][i2]+cw[1][(2*i1+l+2*j)%(2*j)][i2];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=cw[0][i1][i2];
    cw[0][i1][i2]=0.;
    c[1][i1][i2]=cw[1][i1][i2];
    cw[1][i1][i2]=0.;
   }
  }
 } 
}

/*********************************************************************/

/* Transformee en ondelettes gradients (projection oblique sur Hrot) */

void rot_ani_trans(fltarray & Tabc, fltarray & Tabd, fltarray & Tabddn)
{
  int N = Tabc.axis(1);
  float ***c;
  float ***d;
  float ***dg;
  float *k1 = new float[N];
  float *k2 = new float [N];
 
  set_ptr3d(Tabc, c);
  Tabd.resize(N,N,2);
  Tabddn.resize(N,N,2);
  set_ptr3d(Tabd, d);
  set_ptr3d(Tabddn, dg);
 
 int j,j1,j2,i1,i2;

 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   k1[i1]=c[0][i1][i2];
   k2[i1]=c[1][i1][i2];
  }
  for (j=N/2;j>0;j=j/2) {
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=-1./8.*k1[(2*i1-2+2*j)%(2*j)]+1./4.*k1[(2*i1-1+2*j)%(2*j)]+3./4.*k1[2*i1]+1./4.*k1[2*i1+1]-1./8.*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=-1./4.*k1[2*i1]+1./2.*k1[2*i1+1]-1./4.*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=-1./4.*k2[(2*i1-1+2*j)%(2*j)]+3./4.*k2[2*i1]+3./4.*k2[2*i1+1]-1./4.*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=1./8.*k2[(2*i1-1+2*j)%(2*j)]-3./8.*k2[2*i1]+3./8.*k2[2*i1+1]-1./8.*k2[(2*i1+2)%(2*j)];
   }
   for (i1=0;i1<j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
  }
 }
 for (i1=0;i1<N;i1++) {
  for (i2=0;i2<N;i2++) {
   k1[i2]=d[0][i1][i2];
   k2[i2]=d[1][i1][i2];
  }
  for (j=N/2;j>0;j=j/2) {
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=-1./4.*k1[(2*i2-1+2*j)%(2*j)]+3./4.*k1[2*i2]+3./4.*k1[2*i2+1]-1./4.*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=1./8.*k1[(2*i2-1+2*j)%(2*j)]-3./8.*k1[2*i2]+3./8.*k1[2*i2+1]-1./8.*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=-1./8.*k2[(2*i2-2+2*j)%(2*j)]+1./4.*k2[(2*i2-1+2*j)%(2*j)]+3./4.*k2[2*i2]+1./4.*k2[2*i2+1]-1./8.*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=-1./4.*k2[2*i2]+1./2.*k2[2*i2+1]-1./4.*k2[(2*i2+2)%(2*j)];
   }
   for (i2=0;i2<j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
  }
 }
  fits_write_fltarr("acur.fits", Tabd);

 for (j1=N/2;j1>0;j1=j1/2) {
  for (j2=N/2;j2>0;j2=j2/2) {
   for (i1=0;i1<j1;i1++) {
    for (i2=0;i2<j2;i2++) {
     dg[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     dg[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
    }
   }
  }
 }
 for (i1=1;i1<N;i1++) {
  dg[0][i1][0]=d[0][i1][0];
  dg[1][i1][0]=d[1][i1][0];
  dg[0][0][i1]=d[1][0][i1];
  dg[1][0][i1]=d[0][0][i1];
 }
 dg[0][0][0]=d[0][0][0];
 dg[1][0][0]=d[1][0][0];
 
}

/*********************************************************************/
 
void rot_ani_trans(fltarray & Tabc, fltarray & Tabddn)
{
   fltarray Tabd;
   rot_ani_trans(Tabc, Tabd, Tabddn);
}


/*********************************************************************/

/* Transformee inverse partant des ondelettes gradients */

void rot_ani_rec(fltarray & Tabddn, fltarray & Tabc)
{
  int N = Tabddn.axis(1);
  float ***c;
  float ***dg;
  float ***cw;
  
  Tabc.resize(N,N,2);
  Tabddn.resize(N,N,2);
  set_ptr3d(Tabc, c);   
  init_ptr3d(Tabddn, cw);
  set_ptr3d(Tabddn, dg);

 float *k1 = new float[N];
 float *k2 = new float [N];
   
 int j,j1,j2,i1,i2,l;
 float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
 float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 
 c[0][0][0]=dg[0][0][0];
 c[1][0][0]=dg[1][0][0];
 for (i1=1;i1<N;i1++) {
  c[0][i1][0]=dg[0][i1][0];
  c[1][i1][0]=dg[1][i1][0];
  c[0][0][i1]=dg[1][0][i1];
  c[1][0][i1]=dg[0][0][i1];
 }
 for (j1=N/2;j1>0;j1=j1/2) {
  for (j2=N/2;j2>0;j2=j2/2) {
   for (i1=0;i1<j1;i1++) {
    for (i2=0;i2<j2;i2++) {
     c[0][j1+i1][j2+i2]=j1*dg[1][j1+i1][j2+i2]+j2*dg[0][j1+i1][j2+i2];
     c[1][j1+i1][j2+i2]=j2*dg[1][j1+i1][j2+i2]-j1*dg[0][j1+i1][j2+i2];
    }
   }
  }
 }
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   cw[0][i1][i2]=0.;
   cw[1][i1][i2]=0.;
  }
 }
 for (i1=0;i1<N;i1++) {
  for (j=1;j<N;j=2*j) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     cw[0][i1][(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][i1][j+i2]+cw[0][i1][(2*i2+l+2*j)%(2*j)];
     cw[1][i1][(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][i1][j+i2]+cw[1][i1][(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=cw[0][i1][i2];
    cw[0][i1][i2]=0.;
    c[1][i1][i2]=cw[1][i1][i2];
    cw[1][i1][i2]=0.;
   }
  }
 }
 for (i2=0;i2<N;i2++) {
  for (j=1;j<N;j=2*j) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     cw[0][(2*i1+l+2*j)%(2*j)][i2]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][j+i1][i2]+cw[0][(2*i1+l+2*j)%(2*j)][i2];
     cw[1][(2*i1+l+2*j)%(2*j)][i2]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][j+i1][i2]+cw[1][(2*i1+l+2*j)%(2*j)][i2];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=cw[0][i1][i2];
    cw[0][i1][i2]=0.;
    c[1][i1][i2]=cw[1][i1][i2];
    cw[1][i1][i2]=0.;
   }
  }
 }
}

/*********************************************************************/

/* Decomposition de Helmholtz par ondelettes splines a divergence nulle et gradient 
 void div_ani_trans(fltarray & Tabc, fltarray & Tabd, fltarray & Tabddn)

  void hodge(float vn[2][N][N]  *donnee*,float d[2][N][N] *resultat*,float vn0[2][N][N] *first guess*,float d0[2][N][N]*first guess*,
   int  aIt nb d'it*)
*/

void hodge(fltarray & Tabvn , fltarray & Tabd  ,   fltarray & Tabvn0  , fltarray & Tabd0 , int It )
{
  int N = Tabvn.axis(1);
  float ***vn;
  float ***vn0;
  float ***d;
  float ***d0;
 // init_ptr3d(Tabvn, vn);
//  put2_ptr3d(Tabvn, vn);
//  init_ptr3d(Tabvn0, vn0);
//  put2_ptr3d(Tabvn0, vn0);
//  init_ptr3d(Tabd0, d0);
//  put2_ptr3d(Tabvn0, vn0);  
  set_ptr3d(Tabvn, vn);
  set_ptr3d(Tabvn0, vn0);
  set_ptr3d(Tabd0, d0);
  set_ptr3d(Tabd, d);
  fltarray Tabc,Tabcn,Tabdn;
  float ***c;
  float ***cn;
  float ***dn;
  Tabc.resize(N,N,2);
  set_ptr3d(Tabc, c);
  Tabcn.resize(N,N,2);
  set_ptr3d(Tabcn, cn);
  Tabdn.resize(N,N,2);
  set_ptr3d(Tabdn, dn);
  int j,i,i1,i2;
  /* temp */ float mx;
 
  // float c[2][N][N],cn[2][N][N],dn[2][N][N];
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=0.;
   d[1][i1][i2]=0.;
   mx=vn[0][i1][i2];vn[0][i1][i2]=vn[1][i1][i2];vn[1][i1][i2]=mx;
  }
 }
  
 for (i=0;i<It;i++) 
 {
  cout << " Iter " << i+1 << endl;
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
   c[0][i1][i2]=5./8.*(vn[0][i1][i2]+vn[0][(i1+1)%N][i2])-1./8.*(vn[0][(i1-1+N)%N][i2]+vn[0][(i1+2)%N][i2]);
   c[1][i1][i2]=5./8.*(vn[1][i1][i2]+vn[1][i1][(i2+1)%N])-1./8.*(vn[1][i1][(i2-1+N)%N]+vn[1][i1][(i2+2)%N]);
   }
  }
  // tfodanp2d(c,dn);
  div_ani_trans(Tabc, Tabdn);
   for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    d[0][i1][i2]=d[0][i1][i2]+dn[0][i1][i2];
    d[1][i1][i2]=d[1][i1][i2]+dn[1][i1][i2];
    if ((i1+i2)!=0) {
     d[1][i1][i2]=d[1][i1][i2]-dn[1][i1][i2];
     dn[1][i1][i2]=0.;
    }
   }
  }
   // tfodanpinv2d(c,dn);
  div_ani_rec(Tabdn, Tabc);

  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    vn[0][i1][i2]=vn[0][i1][i2]-1./2.*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2]);
    vn[1][i1][i2]=vn[1][i1][i2]-1./2.*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2]);
   }
  }
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
   c[0][i1][i2]=5./8.*(vn[0][i1][i2]+vn[0][i1][(i2+1)%N])-1./8.*(vn[0][i1][(i2-1+N)%N]+vn[0][i1][(i2+2)%N]);
   c[1][i1][i2]=5./8.*(vn[1][i1][i2]+vn[1][(i1+1)%N][i2])-1./8.*(vn[1][(i1-1+N)%N][i2]+vn[1][(i1+2)%N][i2]);
   }
  }
  // tfoganp2d(c,dn);
  rot_ani_trans(Tabc, Tabdn);

  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    if ((i1+i2)!=0) {
     d[1][i1][i2]=d[1][i1][i2]+dn[1][i1][i2];
    }
    dn[0][i1][i2]=0.;
    if ((i1+i2)==0) dn[1][i1][i2]=0.;
   }
  }
  // tfoganpinv2d(c,dn);
  rot_ani_rec( Tabdn,  Tabc);
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    vn[0][i1][i2]=vn[0][i1][i2]-1./2.*(c[0][i1][(i2-1+N)%N]+c[0][i1][i2]);
    vn[1][i1][i2]=vn[1][i1][i2]-1./2.*(c[1][(i1-1+N)%N][i2]+c[1][i1][i2]);
   }
  }
 } // end loop iter
 
 
 /*actualisation de la donnée initiale*/
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d0[0][i1][i2]=d[0][i1][i2];
   d0[1][i1][i2]=0.;
   dn[1][i1][i2]=d[1][i1][i2];
   dn[0][i1][i2]=0.;
  }
 }
 d0[1][0][0]=d[1][0][0];
 dn[1][0][0]=0.;
 // tfodanpinv2d(c,d0);
 div_ani_rec(Tabd0, Tabc);

 // tfoganpinv2d(cn,dn);
 rot_ani_rec( Tabdn,  Tabcn);
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   vn0[0][i1][i2]=0.5*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2])+0.5*(cn[0][i1][(i2-1+N)%N]+cn[0][i1][i2]);
   vn0[1][i1][i2]=0.5*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2])+0.5*(cn[1][(i1-1+N)%N][i2]+cn[1][i1][i2]);
   d0[0][i1][i2]=d[0][i1][i2];
   d0[1][i1][i2]=d[1][i1][i2];
  }
 }
                fltarray  E, B;
				E.resize(N, N);
				B.resize(N, N);
			    for (int i1=0; i1 < N; i1++)
		        for (int i2=0; i2 < N; i2++)
		        {
		           E(i2,i1) =  -0.5*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2])+0.5*(cn[0][i1][(i2-1+N)%N]+cn[0][i1][i2]);
			       B(i2,i1) =  -0.5*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2])+0.5*(cn[1][(i1-1+N)%N][i2]+cn[1][i1][i2]);
		        }
                fits_write_fltarr("Emode", B);
	            fits_write_fltarr("Bmode", E);
}

/*********************************************************************/

void hodge(fltarray & Tabvn , fltarray & Tabd  , int It )
{
   int N = Tabvn.axis(1);
   fltarray Tabvn0;
   fltarray Tabd0;
    Tabvn0.resize(N,N,2);
    Tabd0.resize(N,N,2);
    Tabvn0.init();
    Tabd0.init();
    hodge( Tabvn ,  Tabd  ,  Tabvn0  ,   Tabd0 , It );

}


/*********************************************************************/

void quasi_interpol_div(fltarray & Tabp, fltarray & Tabc)
{
 
  int N = Tabp.axis(1);
  float ***p;
  float ***c;
   
  init_ptr3d(Tabp, p);   
  put2_ptr3d(Tabp, p);
   
  Tabc.resize(N,N,2);
  init_ptr3d(Tabc,c);
   
   int i1,i2;
   float f[4]={-1./8.,5./8.,5./8.,-1./8.};
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++) 
  c[0][i1][i2]=f[0]*p[0][(i1+N-1)%N][i2]+f[1]*p[0][i1][i2]+f[2]*p[0][(i1+1)%N][i2]+f[3]*p[0][(i1+2)%N][i2];
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++)
  c[1][i1][i2]=f[0]*p[1][i1][(i2+N-1)%N]+f[1]*p[1][i1][i2]+f[2]*p[1][i1][(i2+1)%N]+f[3]*p[1][i1][(i2+2)%N];
  
 get2_ptr3d(Tabc,  c);
 
}

/*********************************************************************/

void quasi_interpol_curl(fltarray & Tabp, fltarray & Tabc)
{
 
  int N = Tabp.axis(1);
  float ***p;
  float ***c;
   
  init_ptr3d(Tabp, p);   
  put2_ptr3d(Tabp, p);
   
  Tabc.resize(N,N,2);
  init_ptr3d(Tabc,c);
   
   int i1,i2;
   float f[4]={-1./8.,5./8.,5./8.,-1./8.};
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++) 
  c[0][i1][i2]=f[0]*p[0][i1][(i2+N-1)%N]+f[1]*p[0][i1][i2]+f[2]*p[0][i1][(i2+1)%N]+f[3]*p[0][i1][(i2+2)%N];
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++)
  c[1][i1][i2]=f[0]*p[1][(i1+N-1)%N][i2]+f[1]*p[1][i1][i2]+f[2]*p[1][(i1+1)%N][i2]+f[3]*p[1][(i1+2)%N][i2];
  
 get2_ptr3d(Tabc,  c);
 
}

/*********************************************************************/

void point_value_div(fltarray & Tabc, fltarray & Tabp)
{
 
  int N = Tabc.axis(1);
  float ***p;
  float ***c;
   
  init_ptr3d(Tabc, c);   
  put2_ptr3d(Tabc, c);
   
  Tabp.resize(N,N,2);
  init_ptr3d(Tabp,p);
   
   int i1,i2;
   float f[2]={1./2.,1./2.};
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++) 
  p[0][i1][i2]=f[0]*c[0][(i1+N-1)%N][i2]+f[1]*c[0][i1][i2];
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++)
  p[1][i1][i2]=f[0]*c[1][i1][(i2+N-1)%N]+f[1]*c[1][i1][i2];
  
 get2_ptr3d(Tabp,  p);
 
}

/*********************************************************************/

void point_value_curl(fltarray & Tabc, fltarray & Tabp)
{
 
  int N = Tabc.axis(1);
  float ***p;
  float ***c;
   
  init_ptr3d(Tabc, c);   
  put2_ptr3d(Tabc, c);
   
  Tabp.resize(N,N,2);
  init_ptr3d(Tabp,p);
   
   int i1,i2;
   float f[2]={1./2.,1./2.};
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++) 
  p[0][i1][i2]=f[0]*c[0][i1][(i2+N-1)%N]+f[1]*c[0][i1][i2];
  
 for (i2=0;i2<N;i2++) 
 for (i1=0;i1<N;i1++)
  p[1][i1][i2]=f[0]*c[1][(i1+N-1)%N][i2]+f[1]*c[1][i1][i2];
  
 get2_ptr3d(Tabp,  p);
 
}

/*********************************************************************/
