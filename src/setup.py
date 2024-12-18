from setuptools import setup, Extension, find_packages
import pybind11

pybind11_include = pybind11.get_include()
# run with python3 setup.py build_ext --inplace

ext_modules = [
    Extension(
        'hip.debruijn', 
        sources=['hip/debruijn.cpp'], 
        include_dirs=[pybind11_include], 
        language='c++', 
        extra_compile_args=['-std=c++11', '-O3'], 
    ), 
    Extension(
        'hip.hpMap', 
        sources=['hip/hpMap.cpp'], 
        include_dirs=[pybind11_include], 
        language='c++', 
        extra_compile_args=['-std=c++11', '-O3'], 
    )
]

setup(
    name='hip',
    version='0.1',
    ext_modules=ext_modules, 
    install_requires=['pybind11'], 
)
