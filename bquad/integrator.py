import numpy as np
import cffi, glob, os

bquad_dir = os.path.dirname(__file__)
include_dir = os.path.join(bquad_dir,'include')
lib_file = os.path.join(bquad_dir,'_bquad.so')
# Some installation (e.g. Travis with python 3.x)
# name this e.g. _bquad.cpython-34m.so,
# so if the normal name doesn't exist, look for something else.
if not os.path.exists(lib_file):
    alt_files = glob.glob(os.path.join(os.path.dirname(__file__),'_bquad*.so'))
    if len(alt_files) == 0:
        raise IOError("No file '_bquad.so' found in %s"%bquad_dir)
    if len(alt_files) > 1:
        raise IOError("Multiple files '_bquad*.so' "+\
                      "found in %s: %s"%(bquad_dir,alt_files))
    lib_file = alt_files[0]

_ffi = cffi.FFI()
for file_name in glob.glob(os.path.join(include_dir,'*.h')):
    _ffi.cdef(open(file_name).read())
_lib = _ffi.dlopen(lib_file)

class Integrator(object):
    """
    Performs integrals using the Bessel-Quadrature rule.
    """
    def __init__(self):
        pass
