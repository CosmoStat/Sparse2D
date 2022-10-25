/*##########################################################################
pySAP - Copyright (C) CEA, 2017 - 2018
Distributed under the terms of the CeCILL-B license, as published by
the CEA-CNRS-INRIA. Refer to the LICENSE file or to
http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
for details.
##########################################################################*/

#ifndef NUMPYHEALPIX_H_
#define NUMPYHEALPIX_H_

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include "IM_Obj.h"
#include "IM_IO.h"
#include <iostream>
#include "HealpixClass.h"

namespace py = pybind11;

// Helper function for fast image to array conversion
py::array_t<float> hmap2array(Hmap<REAL> & Map){
    int  Npix= Map.Npix();
    auto arr1 = py::array_t<float>(Npix);
    auto buf1 = arr1.request();
    float *pointer = (float *) buf1.ptr;
    for (int i=0; i<Npix; i++)
            pointer[i] = Map[i];
    return arr1;
}

// Helper function for fast arrat to image conversion
Hmap<REAL> array2hmap(py::array_t<float> &array){
  if (array.ndim() != 1)
    throw std::runtime_error("Input should be 1-D NumPy array");

  auto buffer = array.request();
  float *pointer = (float *) buffer.ptr;
  int Npix = array.shape(0);
  int Nside =floor(sqrt((double) Npix / 12.));
  Hdmap Map;
  Map.alloc(Nside);
  for (int i=0; i<Npix; i++)
            Map[i] = pointer[i];
  return Map;
}

#endif
