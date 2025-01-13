import os
import pybind11
from setuptools import setup, Extension, find_packages

pybind11_include = pybind11.get_include()
# run with python3 setup.py build_ext --inplace

robin_map_dir = os.path.join(os.getcwd(), 'robin-map')

ext_modules = [
    Extension(
        'hip.debruijn', 
        sources=['hip/debruijn.cpp'], 
        include_dirs=[pybind11_include], 
        language='c++', 
        extra_compile_args=['-std=c++17', '-O3'], 
    ), 
    Extension(
        'hip.validator', 
        sources=['hip/validator.cpp'], 
        include_dirs=[pybind11_include, f"{robin_map_dir}/include"], 
        language='c++', 
        extra_compile_args=['-std=c++17', '-O3'], 
    ), 
    Extension(
        'hip.tablemer', 
        sources=['hip/tablemer.cpp'], 
        include_dirs=[pybind11_include, f"{robin_map_dir}/include"], 
        language='c++', 
        extra_compile_args=['-std=c++17', '-O3'], 
    )
]

setup(
    name='hip',
    version='0.1',
    ext_modules=ext_modules, 
    install_requires=['pybind11'], 
)
