/****************************************************************

This is the famous FFTLog. 

First imlplemented by the living legend Andrew Hamilton:

http://casa.colorado.edu/~ajsh/FFTLog/

This version is a C version that was adapted from the C++ version found
in Copter JWG Carlson, another big loss for the cosmology community.

https://github.com/jwgcarlson/Copter

I've transformed this from C++ to C99 as the lowest common denominator
and provided bindings for C++ and python.

These are the C++ bindings

*****************************************************************/

/* Compute the correlation function xi(r) from a power spectrum P(k), sampled
 * at logarithmically spaced points k[j]. */
void pk2xi(int N,  const double k[],  const double pk[], double r[], double xi[]);

/* Compute the power spectrum P(k) from a correlation function xi(r), sampled
 * at logarithmically spaced points r[i]. */
void xi2pk(int N,  const double r[],  const double xi[], double k[], double pk[]);

/* Compute the function
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
 * in this notation.  The input k-values must be logarithmically spaced.  The
 * resulting xi_l^m(r) will be evaluated at the dual r-values
 *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
void fftlog_ComputeXiLM(double l, double m, int N, const double k[],  const double pk[],
			double r[], double xi[]);

/* Compute the function
 *   \xi_\alpha(\theta) = \int_0^\infty \frac{d\ell}{2\pi} \ell J_\alpha(\ell\theta) C_\ell
 * The input l-values must be logarithmically spaced.  The
 * resulting xi_alpha(th) will be evaluated at the dual th-values
 *   th[0] = 1/l[N-1], ..., th[N-1] = 1/l[0]. */
void fftlog_ComputeXi2D(double bessel_order,int N,const double l[],const double cl[],
			double th[], double xi[]);
