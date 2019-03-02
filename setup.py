
import sys, os, glob
import setuptools
from setuptools import setup, Extension
import subprocess

os.system('ln -h -s ../include bquad/include')

sources = glob.glob(os.path.join('src','*.c'))
headers = glob.glob(os.path.join('include','*.h'))
try:
    cflags = subprocess.check_output(['gsl-config', '--cflags'], universal_newlines=True).split()
    lflags = subprocess.check_output(['gsl-config', '--libs'], universal_newlines=True).split()
except OSError:
    raise Exception("Error: must have GSL installed and gsl-config working")

ext=Extension("bquad._bquad",
              sources,
              depends=headers,
              include_dirs=['include'],
              extra_compile_args=[os.path.expandvars(flag) for flag in cflags],
              extra_link_args=[os.path.expandvars(flag) for flag in lflags])

dist = setup(name="bquad",
             author="Tom McClintock",
             author_email="mcclintock@bnl.gov",
             description="Modules for modeling galaxy clusters and related systematics.",
             license="MIT License",
             url="https://github.com/tmcclintock/Bessel_Quadrature_Integrator",
             packages=['bquad'],
             package_data={'bquad' : headers },
             ext_modules=[ext],
             install_requires=['cffi','numpy'],
             setup_requires=['pytest_runner'],
             tests_require=['pytest']
)

#setup.py doesn't put the .so file in the bquad directory, 
#so this bit makes it possible to
#import bquad from the root directory.  
#Not really advisable, but everyone does it at some
#point, so might as well facilitate it.
build_lib = glob.glob(os.path.join('build','*','bquad','_bquad*.so'))
if len(build_lib) >= 1:
    lib = os.path.join('bquad','_bquad.so')
    if os.path.lexists(lib): os.unlink(lib)
    os.link(build_lib[0], lib)
