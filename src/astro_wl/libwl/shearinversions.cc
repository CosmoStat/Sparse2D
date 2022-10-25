/*
 * cln.h - This file is part of MRS3D
 * Created on 16/05/11
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Copyright 2012 CEA
 *
 * This software is a computer program whose purpose is to apply mutli-
 * resolution signal processing algorithms on spherical 3D data.
 *
 * This software is governed by the CeCILL  license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 */

#include "shearinversions.h"

ShearInversions::ShearInversions(int NZeroPadded)
{
    NpixFFT=NZeroPadded;
    // Normalization factor for the fft
    fftFactor =1.0/(((double)NpixFFT)*NpixFFT);

    fftFrame1  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NpixFFT*NpixFFT);
    fftFrame2  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NpixFFT*NpixFFT);
    fftFrame3  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NpixFFT*NpixFFT);
    fftFrame4  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NpixFFT*NpixFFT);

    // Create the FFTW3 plans
    plan_forward1  = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame1, fftFrame1, FFTW_FORWARD,  FFTW_MEASURE);
    plan_backward1 = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame1, fftFrame1, FFTW_BACKWARD, FFTW_MEASURE);
    plan_forward2  = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame2, fftFrame2, FFTW_FORWARD,  FFTW_MEASURE);
    plan_backward2 = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame2, fftFrame2, FFTW_BACKWARD, FFTW_MEASURE);
    plan_forward3  = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame3, fftFrame3, FFTW_FORWARD,  FFTW_MEASURE);
    plan_forward4  = fftw_plan_dft_2d(NpixFFT, NpixFFT, fftFrame4, fftFrame4, FFTW_FORWARD,  FFTW_MEASURE);
    
}


ShearInversions::~ShearInversions()
{
    fftw_free(fftFrame1);
    fftw_free(fftFrame2);
    fftw_free(fftFrame3);
    fftw_free(fftFrame4);
    fftw_destroy_plan(plan_forward1);
    fftw_destroy_plan(plan_forward2);
    fftw_destroy_plan(plan_backward1);
    fftw_destroy_plan(plan_backward2);
    fftw_destroy_plan(plan_forward3);
    fftw_destroy_plan(plan_forward4);
}


void ShearInversions::gamma2kappa(double* pt_gamma1, double* pt_gamma2, double* pt_kappa, int Nx, int Ny, int Nz)
{
    for(int z =0; z <Nz ; z++) {
        double k1,k2,k1k1,k2k2,k1k2;
        double denom;

        // Position the data inside the zero Padded Frame
        for(int x=0; x < NpixFFT ; x++)
            for(int y=0; y < NpixFFT ; y++) {
                if(x < Nx && y< Ny) {
                    fftFrame1[x*NpixFFT +y][0] = pt_gamma1[z *Nx*Ny + x*Ny +y];
                    fftFrame2[x*NpixFFT +y][0] = pt_gamma2[z *Nx*Ny + x*Ny +y];
                } else {
                    fftFrame1[x*NpixFFT +y][0] = 0.0;
                    fftFrame2[x*NpixFFT +y][0] = 0.0;
                }
                fftFrame1[x*NpixFFT +y][1] = 0.0;
                fftFrame2[x*NpixFFT +y][1] = 0.0;
            }

        fftw_execute(plan_forward1);
        fftw_execute(plan_forward2);
        for(int x=0; x < NpixFFT; x++) {
            k1 = (x <= NpixFFT/2 ? x : x - NpixFFT);
            for(int y=0; y < NpixFFT ; y++) {
                if(x+y == 0) continue;
                long pos = x*(NpixFFT) + y;

                k2 = (y <= NpixFFT/2 ? y : y - NpixFFT);
                k1k1=k1*k1;
                k2k2=k2*k2;
                k1k2=k1*k2;
                denom = 1.0/(k1*k1 + k2*k2);
                fftFrame1[pos][0] = fftFrame1[pos][0]*((k2k2 - k1k1)*denom) + fftFrame2[pos][0]*(2.0*k1k2*denom);
                fftFrame1[pos][1] = fftFrame1[pos][1]*((k2k2 - k1k1)*denom) + fftFrame2[pos][1]*(2.0*k1k2*denom);
            }
        }
        fftFrame1[0][0] = 0.0;
        fftFrame1[0][1] = 0.0;

        fftw_execute(plan_backward1);

        for(int x=0; x < Nx ; x++)
            for(int y=0; y < Ny ; y++)
                pt_kappa[z*Nx*Ny + x*Ny +y] = fftFrame1[x*NpixFFT +y][0]*fftFactor;
    }
}

void ShearInversions::gamma2kappa(dblarray& gamma1, dblarray& gamma2, dblarray& kappa)
{
    double *pt_gamma1 = gamma1.buffer();
    double *pt_gamma2 = gamma2.buffer();
    double *pt_kappa  = kappa.buffer();

    int Nx = gamma1.nx();
    int Ny = gamma1.ny();
    int Nz = gamma1.nz();

    gamma2kappa(pt_gamma1,pt_gamma2,pt_kappa,Nx,Ny,Nz+1);
}

void ShearInversions::kappa2gamma(double* pt_kappa, double* pt_gamma1, double* pt_gamma2, int Nx, int Ny, int Nz)
{
    for(int z =0; z <Nz ; z++) {
        double k1,k2,k1k1,k2k2,k1k2;
        double denom;

        // Position the data inside the zero Padded Frame
        for(int x=0; x < NpixFFT ; x++)
            for(int y=0; y < NpixFFT ; y++) {
                if(x < Nx && y< Ny) {
                    fftFrame1[x*NpixFFT +y][0] = pt_kappa[z *Nx*Ny + x*Ny +y];
                } else {
                    fftFrame1[x*NpixFFT +y][0] = 0.0;
                }
                fftFrame1[x*NpixFFT +y][1] = 0.0;
            }


        // Compute the Fourier Transform
        fftw_execute(plan_forward1);

        // Performing the inversion in Fourier Space
        for(int x=0; x < NpixFFT ; x++) {
            k1 = (x <= NpixFFT/2 ? x : x - NpixFFT);
            for(int y=0; y < NpixFFT ; y++) {
                if(x+y == 0) continue;
                long pos = x*NpixFFT + y;


                k2 = (y <= NpixFFT/2 ? y : y - NpixFFT);
                k1k1=k1*k1;
                k2k2=k2*k2;
                k1k2=k1*k2;
                denom = 1.0/(k1*k1 + k2*k2);
                fftFrame2[pos][0] = fftFrame1[pos][0] * (2.0*k1k2*denom);
                fftFrame2[pos][1] = fftFrame1[pos][1] * (2.0*k1k2*denom);
                fftFrame1[pos][0] *= ((k2k2 - k1k1)*denom);
                fftFrame1[pos][1] *= ((k2k2 - k1k1)*denom);
            }
        }


        fftFrame1[0][0] = 0.0;
        fftFrame1[0][1] = 0.0;
        fftFrame2[0][0] = 0.0;
        fftFrame2[0][1] = 0.0;

        //Performing the inverse Fourier transforms
        fftw_execute(plan_backward1);
        fftw_execute(plan_backward2);

        // Extract gamma1 and gamma2
        for(int x=0; x < Nx ; x++)
            for(int y=0; y < Ny ; y++) {
                pt_gamma1[z *Nx*Ny + x*Ny +y] = fftFrame1[x*NpixFFT +y][0]*fftFactor;
                pt_gamma2[z *Nx*Ny + x*Ny +y] = fftFrame2[x*NpixFFT +y][0]*fftFactor;
            }
    }

}

void ShearInversions::kappa2gamma(to_array< double, true >& kappa, to_array< double, true >& gamma1, to_array< double, true >& gamma2)
{
    double *pt_gamma1 = gamma1.buffer();
    double *pt_gamma2 = gamma2.buffer();
    double *pt_kappa  = kappa.buffer();

    int Nx = kappa.nx();
    int Ny = kappa.ny();
    int Nz = kappa.nz();

    kappa2gamma(pt_kappa,pt_gamma1,pt_gamma2,Nx,Ny,Nz+1);
}

void ShearInversions::kappa2flexion(double pixel_size, to_array< double, true >& kappa, to_array< double, true >& F1, to_array< double, true >& F2)
{
    double *pt_F1 = F1.buffer();
    double *pt_F2 = F2.buffer();
    double *pt_kappa  = kappa.buffer();

    int Nx = kappa.nx();
    int Ny = kappa.ny();
    int Nz = kappa.nz();

    kappa2flexion(pixel_size, pt_kappa,pt_F1,pt_F2,Nx,Ny,Nz+1);
}

void ShearInversions::kappa2flexion(double pixel_size, double* pt_kappa, double* pt_F1, double* pt_F2, int Nx, int Ny, int Nz)
{
    for(int z =0; z <Nz ; z++) {
        double k1,k2,k1k1,k2k2,k1k2;
        double denom;
        double temp;
        // Position the data inside the zero Padded Frame
        for(int x=0; x < NpixFFT ; x++)
            for(int y=0; y < NpixFFT ; y++) {
                if(x < Nx && y< Ny) {
                    fftFrame1[x*NpixFFT +y][0] = pt_kappa[z *Nx*Ny + x*Ny +y];
                } else {
                    fftFrame1[x*NpixFFT +y][0] = 0.0;
                }
                fftFrame1[x*NpixFFT +y][1] = 0.0;
            }


        // Compute the Fourier Transform
        fftw_execute(plan_forward1);

        // Performing the inversion in Fourier Space
        for(int x=0; x < NpixFFT ; x++) {
            k1 = (x <= NpixFFT/2 ? x : x - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT );
            for(int y=0; y < NpixFFT ; y++) {
                //if(x+y == 0) continue;
                long pos = x*NpixFFT + y;

                k2 = (y <= NpixFFT/2 ? y : y - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT );
                fftFrame2[pos][0] = -fftFrame1[pos][1] * k2;
                fftFrame2[pos][1] = fftFrame1[pos][0] * k2;
                temp = fftFrame1[pos][0];
                fftFrame1[pos][0] = - fftFrame1[pos][1]*k1;
                fftFrame1[pos][1] = temp *k1;
            }
        }

/*
        fftFrame1[0][0] = 0.0;
        fftFrame1[0][1] = 0.0;
        fftFrame2[0][0] = 0.0;
        fftFrame2[0][1] = 0.0;*/

        //Performing the inverse Fourier transforms
        fftw_execute(plan_backward1);
        fftw_execute(plan_backward2);

        // Extract gamma1 and gamma2
        for(int x=0; x < Nx ; x++)
            for(int y=0; y < Ny ; y++) {
                pt_F1[z *Nx*Ny + x*Ny +y] = fftFrame1[x*NpixFFT +y][0]*fftFactor;
                pt_F2[z *Nx*Ny + x*Ny +y] = fftFrame2[x*NpixFFT +y][0]*fftFactor;
            }
    }
}


void ShearInversions::flexion2kappa(double pixel_size,to_array< double, true >& F1, to_array< double, true >& F2, to_array< double, true >& kappa)
{
    double *pt_F1 = F1.buffer();
    double *pt_F2 = F2.buffer();
    double *pt_kappa  = kappa.buffer();

    int Nx = F1.nx();
    int Ny = F1.ny();
    int Nz = F1.nz();

    flexion2kappa(pixel_size, pt_F1,pt_F2,pt_kappa,Nx,Ny,Nz+1);
}

void ShearInversions::flexion2kappa(double pixel_size, double* pt_F1, double* pt_F2, double* pt_kappa, int Nx, int Ny, int Nz)
{

    for(int z =0; z <Nz ; z++) {
        double k1,k2,k1k1,k2k2,k1k2;
        double denom;

        // Position the data inside the zero Padded Frame
        for(int x=0; x < NpixFFT ; x++)
            for(int y=0; y < NpixFFT ; y++) {
                if(x < Nx && y< Ny) {
                    fftFrame1[x*NpixFFT +y][0] = pt_F1[z *Nx*Ny + x*Ny +y];
                    fftFrame2[x*NpixFFT +y][0] = pt_F2[z *Nx*Ny + x*Ny +y];
                } else {
                    fftFrame1[x*NpixFFT +y][0] = 0.0;
                    fftFrame2[x*NpixFFT +y][0] = 0.0;
                }
                fftFrame1[x*NpixFFT +y][1] = 0.0;
                fftFrame2[x*NpixFFT +y][1] = 0.0;
            }

        fftw_execute(plan_forward1);
        fftw_execute(plan_forward2);
        for(int x=0; x < NpixFFT; x++) {
            k1 = (x <= NpixFFT/2 ? x : x - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT );
            for(int y=0; y < NpixFFT ; y++) {
                if(x+y == 0) continue;
                long pos = x*(NpixFFT) + y;
                k2 = (y <= NpixFFT/2 ? y : y - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT );

                denom = 1.0/(k1*k1 + k2*k2);
                double dum = denom*(fftFrame1[pos][1]*k1 + fftFrame2[pos][1]*k2);

                fftFrame1[pos][1] = -denom*(fftFrame1[pos][0]*k1 + fftFrame2[pos][0]*k2);
                fftFrame1[pos][0] = dum;
            }
        }
        fftFrame1[0][0] = 0.0;
        fftFrame1[0][1] = 0.0;

        fftw_execute(plan_backward1);

        for(int x=0; x < Nx ; x++)
            for(int y=0; y < Ny ; y++)
                pt_kappa[z*Nx*Ny + x*Ny +y] = fftFrame1[x*NpixFFT +y][0]*fftFactor;
    }

}

void ShearInversions::gamma_flexion2kappa(double pixel_size, to_array< double, true >& gamma1, to_array< double, true >& gamma2, to_array< double, true >& F1, to_array< double, true >& F2, to_array< double, true >& kappa)
{
    double *pt_gamma1 = gamma1.buffer();
    double *pt_gamma2 = gamma2.buffer();
    double *pt_F1 = F1.buffer();
    double *pt_F2 = F2.buffer();
    double *pt_kappa  = kappa.buffer();

    int Nx = F1.nx();
    int Ny = F1.ny();
    int Nz = F1.nz();
    
    double sigma_g = gamma1.sigma();
    double sigma_f = F1.sigma();

    gamma_flexion2kappa(pixel_size, sigma_g, sigma_f, pt_gamma1, pt_gamma2, pt_F1,pt_F2,pt_kappa,Nx,Ny,Nz+1);
}

void ShearInversions::gamma_flexion2kappa(double pixel_size, double sigma_g, double sigma_f, double* gamma1, double* gamma2, double* F1, double* F2, double* kappa, int Nx, int Ny, int Nz)
{
    double sig_frac= pow(sigma_f,2)/pow(sigma_g,2);
    for(int z =0; z <Nz ; z++) {
        double k1,k2,k1k1,k2k2,k1k2, ksqr;
        double denom;

        // Position the data inside the zero Padded Frame
        for(int x=0; x < NpixFFT ; x++)
            for(int y=0; y < NpixFFT ; y++) {
                if(x < Nx && y< Ny) {
                    fftFrame1[x*NpixFFT +y][0] = F1[z *Nx*Ny + x*Ny +y];
                    fftFrame2[x*NpixFFT +y][0] = F2[z *Nx*Ny + x*Ny +y];
                    fftFrame3[x*NpixFFT +y][0] = gamma1[z *Nx*Ny + x*Ny +y];
                    fftFrame4[x*NpixFFT +y][0] = gamma2[z *Nx*Ny + x*Ny +y];
                } else {
                    fftFrame1[x*NpixFFT +y][0] = 0.0;
                    fftFrame2[x*NpixFFT +y][0] = 0.0;
                    fftFrame3[x*NpixFFT +y][0] = 0.0;
                    fftFrame4[x*NpixFFT +y][0] = 0.0;
                }
                fftFrame1[x*NpixFFT +y][1] = 0.0;
                fftFrame2[x*NpixFFT +y][1] = 0.0;
                fftFrame3[x*NpixFFT +y][1] = 0.0;
                fftFrame4[x*NpixFFT +y][1] = 0.0;
            }

        fftw_execute(plan_forward1);
        fftw_execute(plan_forward2);
        fftw_execute(plan_forward3);
        fftw_execute(plan_forward4);
        
        for(int x=0; x < NpixFFT; x++) {
            k1 = (x <= NpixFFT/2 ? x : x - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT);
            for(int y=0; y < NpixFFT ; y++) {
                if(x == 0 && y == 0) continue;
                long pos = x*(NpixFFT) + y;
                k2 = (y <= NpixFFT/2 ? y : y - NpixFFT)*2.0*M_PI/pixel_size/((double) NpixFFT);

                k1k1 = k1*k1; k2k2 = k2*k2; k1k2 = k1*k2; ksqr = k1k1 + k2k2;
                denom = 1.0/(ksqr + sig_frac);
                double frame1Real = fftFrame1[pos][0];
                double frame1Imag = fftFrame1[pos][1];
                
                fftFrame1[pos][0] = denom * (  k1 * frame1Imag + k2 * fftFrame2[pos][1] + sig_frac*( (k2k2 - k1k1)/ksqr * fftFrame3[pos][0] + 2.0*k1k2/ksqr * fftFrame4[pos][0] ) ); 
                fftFrame1[pos][1] = denom * ( -k1 *frame1Real  - k2 * fftFrame2[pos][0] + sig_frac*( (k2k2 - k1k1)/ksqr * fftFrame3[pos][1] + 2.0*k1k2/ksqr * fftFrame4[pos][1] )  );
            }
        }
        fftFrame1[0][0] = 0.0;
        fftFrame1[0][1] = 0.0;

        fftw_execute(plan_backward1);

        for(int x=0; x < Nx ; x++)
            for(int y=0; y < Ny ; y++)
                kappa[z*Nx*Ny + x*Ny +y] = fftFrame1[x*NpixFFT +y][0]*fftFactor;
    }

    
    
    
}


