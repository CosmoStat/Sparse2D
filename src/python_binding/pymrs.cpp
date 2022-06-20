/*##########################################################################
pySAP - Copyright (C) CEA, 2017 - 2018
Distributed under the terms of the CeCILL-B license, as published by
the CEA-CNRS-INRIA. Refer to the LICENSE file or to
http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
for details.
##########################################################################*/

// Includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "mrs.hpp"
#include "alm.hpp"

// Defines a python module which will be named "pysparse"
PYBIND11_MODULE(pymrs, module)
{
  module.doc() = "Python bindings for MRS";
  module.attr("__version__") = "0.1.0";

      py::class_<MRS>(module, "MRS")
    .def(py::init<bool>(),
         py::arg("verbose")=(int)(0)
     )
    .def("alloc", &MRS::alloc,py::arg("nside")=(int)(0),
           py::arg("nscale")=(int)(0),
           py::arg("Lmax")=(int)(0),
           py::arg("ALM_iter")=(int)(0),
           py::arg("Verbose")=(bool)(false))
    .def("info", &MRS::info)
    .def("uwt", &MRS::uwt, py::arg("arr"), py::arg("nscale")=(int)(0))
    .def("iuwt", &MRS::iuwt, py::arg("mr_data"))
    .def("get_tabnorm", &MRS::get_tabnorm)
    .def("hard_filtering", &MRS::hard_filtering, py::arg("arr"),
         py::arg("NSigma")=(float)(0),
         py::arg("SigmaNoise")=(float)(0),
         py::arg("UseMad")=(bool)(false),
         py::arg("KillLastScale")=(bool)(false))

    .def("nside", &MRS::nside)
    
    .def("nscale", &MRS::nscale)
    
    .def("set_alm_iter", &MRS::set_alm_iter, py::arg("ALM_IT")=(int)(0))

    .def("hthres_scale", &MRS::hthres_scale, py::arg("scale")=(int)(0),
           py::arg("NSigma")=(float)(0),
           py::arg("SigmaNoise")=(float)(0),
         py::arg("UseMad")=(bool)(false));
    

    py::class_<C_ALM>(module, "C_ALM")
    .def(py::init<bool>(),
         py::arg("verbose")=(bool)(false)
     )
    .def("alloc", &C_ALM::alloc, py::arg("nside")=(int)(0),
         py::arg("Lmax")=(int)(0),
         py::arg("Verbose")=(bool)(false))
    .def("trans", &C_ALM::trans, py::arg("arr"))
    .def("hthres", &C_ALM::hthres, py::arg("Thres")=(float)(0))
    .def("maxalm", &C_ALM::maxalm)
    .def("use_l2norm", &C_ALM::use_l2norm, py::arg("l2norm")=(bool)(false))
    .def("get_tabalm", &C_ALM::get_tabalm)
    .def("put_tabalm", &C_ALM::put_tabalm, py::arg("arr"))
    .def("get_alm", &C_ALM::get_tabalm)
    .def("put_alm", &C_ALM::put_tabalm, py::arg("arr"))
    .def("wiener", &C_ALM::wiener, py::arg("psn"), py::arg("pss"))
    .def("alm2spec", &C_ALM::alm2spec);
//     .def("recons", &C_ALM::recons);
}
