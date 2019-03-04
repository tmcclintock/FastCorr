import numpy as np
import scipy.interpolate as interp
import cffi, glob, os

fastcorr_dir = os.path.dirname(__file__)
include_dir = os.path.join(fastcorr_dir,'include')
lib_file = os.path.join(fastcorr_dir,'_fastcorr.so')
# Some installation (e.g. Travis with python 3.x)
# name this e.g. _fastcorr.cpython-34m.so,
# so if the normal name doesn't exist, look for something else.
if not os.path.exists(lib_file):
    alt_files = glob.glob(os.path.join(os.path.dirname(__file__),'_fastcorr*.so'))
    if len(alt_files) == 0:
        raise IOError("No file '_fastcorr.so' found in %s"%fastcorr_dir)
    if len(alt_files) > 1:
        raise IOError("Multiple files '_fastcorr*.so' "+\
                      "found in %s: %s"%(fastcorr_dir,alt_files))
    lib_file = alt_files[0]

_ffi = cffi.FFI()
for file_name in glob.glob(os.path.join(include_dir,'*.h')):
    _ffi.cdef(open(file_name).read())
_lib = _ffi.dlopen(lib_file)

def _dc(x):
    """Cast an array correctly with cffi."""
    if isinstance(x, list): x = np.asarray(x, dtype=np.float64, order='C')
    return _ffi.cast('double*', x.ctypes.data)

def general_transform(x, t, F, nu, N = 500, h = 0.005):
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

def general_spherical_transform(x, t, F, nu, N = 500, h = 0.005, reverse=False):
    """Perform a transform of the type
    f(x) = \int_0^\inf dk/2pi^2 k^2 F(k) j_nu(kx)
    
    or, substituting t = kx
    
    f(x) = \int_0^\inf dt/2pi^2 (t/x)^2 F(t/x) j_nu(t).

    If reverse is true, we include a factor of 8pi^3.
    """
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None] #makes x 1D
        scalar_input = True
    if x.ndim > 1:
        raise Exception("x cannot be a >1D array.")

    f = np.zeros_like(x)
    t2F = t*t*F
    _lib.bquad_spherical_transform(_dc(x), len(x), _dc(t), _dc(t2F), len(t),
                                   _dc(f), nu, N, h)
    if reverse:
        f*=8*np.pi**3
    
    if scalar_input:
        return np.squeeze(f/x)
    return f/x

def fftlog_spherical_transform(x, t, F, reverse=False):
    """Uses FFTlog to compute the discrete fourier transform
    of the type
    f(x) = \int_0^\inf dk/2pi^2 k^2 F(k) j_0(kx)
    """
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None] #makes x 1D
        scalar_input = True
    if x.ndim > 1:
        raise Exception("x cannot be a >1D array.")

    #Make a new array, since this routine does a DFT
    xtemp = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), len(t))
    ftemp = np.zeros_like(xtemp)

    if reverse:
        _lib.xi2pk(len(xtemp), _dc(t), _dc(F), _dc(xtemp), _dc(ftemp))
    else:
        _lib.pk2xi(len(xtemp), _dc(t), _dc(F), _dc(xtemp), _dc(ftemp))

    #Re-interpolate at the x values we want
    spl = interp.interp1d(np.log(xtemp), ftemp, fill_value="extrapolate")
    f = spl(np.log(x))
    
    #if reverse:
    #    f*=8*np.pi**3
    
    if scalar_input:
        return np.squeeze(f)
    return f
