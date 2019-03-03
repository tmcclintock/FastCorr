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

def _dc(x):
    """Cast an array correctly with cffi."""
    if isinstance(x, list): x = np.asarray(x, dtype=np.float64, order='C')
    return _ffi.cast('double*', x.ctypes.data)

def general_transform(x, t, F, nu, N = 1000, h = 1e-3):
    """Perform a transform of the type
    f(x) = \int_0^\inf dt/2pi^2 F(t/x) J_nu(t).
    """
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None] #makes x 1D
        scalar_input = True
    if x.ndim > 1:
        raise Exception("x cannot be a >1D array.")

    f = np.zeros_like(x)
    _lib.bquad_transform(_dc(x), len(x), _dc(t), _dc(F), len(t),
                         _dc(f), nu, N, h)
    if scalar_input:
        return np.squeeze(f)
    return f

