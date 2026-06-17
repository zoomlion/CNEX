// MIT License
//
// Copyright (c) 2024 JiangminZheng
//
// Python module binding for validator_core.hpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "validator_core.hpp"

PYBIND11_MODULE(validator, m) {
    m.doc() = "A Python module for validating DNA reads using compressed mer query";

    py::class_<MerQueryManager>(m, "MerQueryManager")
        .def(py::init<>())
        .def("add_mer", &MerQueryManager::add_mer, "Add a mer to the query map",
             py::arg("mer"), py::arg("id"), py::arg("loci"))
        .def("load_from_file", &MerQueryManager::load_from_file,
             "Load mers from TSV file directly (C++ fast path)",
             py::arg("path"))
        .def("get_mer_size", &MerQueryManager::get_mer_size, "Get k-mer size")
        .def("size", &MerQueryManager::size, "Number of mers in the table");

    m.def("validate_read", &validate_read, "Validate a read sequence using compressed mer query",
          py::arg("seq"), py::arg("compressed_mer_query"), py::arg("mer_size"), py::arg("min_c"));
    m.def("dna_encoder", &dna_encoder, "Encode a DNA string to an integer",
          py::arg("dna"));
    m.def("encode_mer_at", &encode_mer_at, "Encode a k-mer at a position in a string (zero-alloc)",
          py::arg("seq"), py::arg("pos"), py::arg("k"), py::arg("result"));
}
