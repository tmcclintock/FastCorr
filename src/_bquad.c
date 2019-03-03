#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_bessel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.141592653589793
#define PI_2 1.5707963267948966 //pi/2
#define PI_times_2 6.28318530718 //2*pi
#define inv_PI_times_2 0.15915494309 //1/(2*pi)

/*Compute pi times the zeros of J_nu where nu is an integer.
 */
int compute_bessel_zeros(double nu, int N, double*xi){
  int i;
  for(i = 0; i < N; i++){
    if (nu == 0){
      xi[i] = gsl_sf_bessel_zero_J0(i+1);
    }
    else if (nu == 1){
      xi[i] = gsl_sf_bessel_zero_J1(i+1);
    }
    else{//nu is not 0 or 1
      xi[i] = gsl_sf_bessel_zero_Jnu(nu, i+1);
    }
  }
  return 0;
}

/*Compute weights of the quadrature routine.
 */
int compute_weights(double nu, int N, double*xi, double*w){
  int i;
  double two_over_pi2 = 2./(PI*PI);
  for(i = 0; i < N; i++){
    if (nu == 0){
      w[i] = two_over_pi2 / (xi[i] * gsl_sf_bessel_J1(PI*xi[i]));
    }else{ //nu > 0
      w[i] = two_over_pi2 / (xi[i] * gsl_sf_bessel_Jnu(nu+1, PI*xi[i]));
    }
  }
  return 0;
}

/* Compute psi and its derivative (dpsi) given N and h.
 */
int compute_Jnu_psi_dpsi(double nu, int N, double h,
			 double*Jnu, double*xi, double*psi, double*dpsi){
  int i;
  double hxi;
  double pi_over_h = PI/h;
  int nu_int = (int)nu; //for speed
  for(i = 0; i < N; i++){
    hxi = h*xi[i];
    psi[i] = hxi * tanh(PI_2 * sinh(hxi));
    dpsi[i] = (PI*hxi*cosh(hxi) + sinh(PI*sinh(hxi)))/(1+cosh(PI*sinh(hxi)));
    if (dpsi[i]!=dpsi[i]) dpsi[i] = 1.0;//avoids NaNs
    if (nu == 0){
      Jnu[i] = gsl_sf_bessel_J0(pi_over_h * psi[i]);
    }else if(nu == 1){
      Jnu[i] = gsl_sf_bessel_J1(pi_over_h * psi[i]);
    }else if(ceil(nu) == nu){//nu is not 0 or 1 but is an integer
      Jnu[i] = gsl_sf_bessel_Jn(nu_int, pi_over_h * psi[i]);
    }else{//nu is not 0 or 1 and is not an integer
      Jnu[i] = gsl_sf_bessel_Jnu(nu, pi_over_h * psi[i]);
    }
  }
  return 0;
}

/* Perform the transform.
 * This specific function transforms functions of the form
 * \int_0^\inf dk/2PI P(k/r)*J_nu(k)
 * If used for hankel transforms or 3D Fourier transforms
 * then the extra factors of k must be handled externally
 * as well as the change of variables from k to k*r.
 */
int bquad_transform(double*r, int Nr, double*k, double*P, int Nk,
		    double*output, double nu, int N, double h){
  int i, j;
  double ki, P_temp, r_inv, sum;
  double pi_over_h = PI/h;

  //WE ASSUME F IS A POWER LAW PAST ITS ENDPOINTS
  double alpha_low = log(P[1]/P[0])/log(k[1]/k[0]);
  double A_low = P[0]/pow(k[0], alpha_low);
  double alpha_hi = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
  double A_hi = P[Nk-1]/pow(k[Nk-1], alpha_hi);

  static int init_flag = 0;
  static gsl_spline*Pspl = NULL;
  static gsl_interp_accel*acc = NULL;
  static double h_old = -1;
  static int N_max = -1;
  static double nu_old = -9999;
  static double*xi = NULL;
  static double*w = NULL;
  static double*Jnu = NULL;
  static double*psi = NULL;
  static double*dpsi = NULL;

  //Create the spline and accelerator
  if (init_flag == 0){
    Pspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
    acc = gsl_interp_accel_alloc();
  }
  gsl_spline_init(Pspl, k, P, Nk);

  //Compute things - also do this whenever h or N changes
  if ((init_flag == 0) || (h_old != h) || (N > N_max) || (nu_old != nu)){

    //These terms change only when N or nu changes
    if ((init_flag == 0) || (N > N_max) || (nu_old != nu)){
      if (xi!=NULL) free(xi);
      if (w!=NULL) free(w);
      xi = malloc(N*sizeof(double));
      w = malloc(N*sizeof(double));
      compute_bessel_zeros(nu, N, xi);
      compute_weights(nu, N, xi, w);
    }
    //These terms also change when h changes
    if (Jnu!=NULL) free(Jnu);
    if (psi!=NULL) free(psi);
    if (dpsi!=NULL) free(dpsi);
    Jnu = malloc(N*sizeof(double));
    psi = malloc(N*sizeof(double));
    dpsi = malloc(N*sizeof(double));

    compute_Jnu_psi_dpsi(nu, N, h, Jnu, xi, psi, dpsi);
    
    h_old = h;
    N_max = N;
    nu_old = nu;
    init_flag = 1;
  }

  //Sum terms appropriately into the output array
  for(j = 0; j < Nr; j++){
    r_inv = 1./r[j];
    sum = 0;
    for(i = 0; i < N; i++){
      ki = pi_over_h * psi[i] * r_inv;
      if (ki < k[0]){
	P_temp = A_low * pow(ki, alpha_low);
      }else if(ki > k[Nk-1]){
	P_temp = A_hi * pow(ki, alpha_hi);
      }else{
	P_temp = gsl_spline_eval(Pspl, ki, acc);
      }
      sum += PI*w[i]*P_temp*Jnu[i]*dpsi[i];
    }
    output[j] = inv_PI_times_2 * sum;
  }

  return 0;
}
