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

#ifndef SHEARINVERSIONS_H
#define SHEARINVERSIONS_H

#include <TempArray.h>
#include <fftw3.h>

class ShearInversions
{

public:
  ShearInversions(int NZeroPadded);
  ~ShearInversions();

  void gamma2kappa(dblarray &gamma1, dblarray &gamma2, dblarray &kappa);
  void gamma2kappa(double *gamma1, double *gamma2, double *kappa, int Nx, int Ny, int Nz);

  void flexion2kappa(double pixel_size, dblarray &F1, dblarray &F2, dblarray &kappa);
  void flexion2kappa(double pixel_size, double *F1, double *F2, double *kappa, int Nx, int Ny, int Nz);

  void kappa2gamma(dblarray &kappa, dblarray &gamma1, dblarray &gamma2);
  void kappa2gamma(double *kappa,double *gamma1, double *gamma2,  int Nx, int Ny, int Nz);

  void kappa2flexion(double pixel_size, dblarray &kappa, dblarray &F1, dblarray &F2);
  void kappa2flexion(double pixel_size, double *kappa,double *F1, double *F2,  int Nx, int Ny, int Nz);

  void gamma_flexion2kappa(double pixel_size, dblarray &gamma1, dblarray &gamma2, dblarray &F1, dblarray &F2, dblarray &kappa);
  void gamma_flexion2kappa(double pixel_size, double sigma_g, double sigma_f, double* gamma1, double* gamma2, double* F1, double* F2, double* kappa, int Nx, int Ny, int Nz);

private:

  int NpixFFT;
  double fftFactor;
  fftw_complex* fftFrame1;
  fftw_complex* fftFrame2;
  fftw_complex* fftFrame3;
  fftw_complex* fftFrame4;
  fftw_plan plan_forward1, plan_backward1,plan_forward2, plan_backward2,plan_forward3,plan_forward4;

};

#endif // SHEARINVERSIONS_H
