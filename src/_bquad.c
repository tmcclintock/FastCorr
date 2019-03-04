#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_bessel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Predefined constants for optimization
#define PI 3.14159265358979323846
#define PI_2 1.5707963267948966 /// pi/2
#define PI_times_2 6.28318530718 // 2*pi
#define inv_PI 0.31830988618 // 1/pi
#define inv_PI_times_2 0.15915494309 // 1/(2*pi)
#define sqrt2 1.41421356237 // sqrt(2)
#define sqrtPI 1.77245385091 // sqrt(pi)
#define inv_sqrtPI 0.56418958354 // 1/sqrt(pi)
#define two_over_PI2 0.20264236728 // 2/pi^2

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

/* Compute the bessel function J_nu, 
 * psi, and its derivative (dpsi) given N and h.
 */
int compute_Jnu_psi_dpsi(double nu, int N, double h,
			 double*Jnu, double*xi, double*psi, double*dpsi){
  int i;
  double hxi;
  double pi_over_h = PI/h;
  int nu_int = (int)nu;
  int nu_is_int = (ceil(nu) == nu); //nu is an integer
  for(i = 0; i < N; i++){
    hxi = h*xi[i];
    psi[i] = hxi * tanh(PI_2 * sinh(hxi));
    dpsi[i] = (PI*hxi*cosh(hxi) + sinh(PI*sinh(hxi)))/(1+cosh(PI*sinh(hxi)));
    if (dpsi[i]!=dpsi[i]) dpsi[i] = 1.0;//avoids NaNs
    if (nu == 0){
      Jnu[i] = gsl_sf_bessel_J0(pi_over_h * psi[i]);
    }else if(nu == 1){
      Jnu[i] = gsl_sf_bessel_J1(pi_over_h * psi[i]);
    }else if(nu_is_int){//nu is not 0 or 1 but is an integer
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
 * Extra factors of k must be handled externally
 * as well as the change of variables from k to k*r.
 */
int bquad_transform(double*r, int Nr, double*k, double*P, int Nk,
		    double*output, double nu, int N, double h){
  int i, j;
  double ki, P_temp, r_inv, sum;
  double pi_over_h = PI/h;

  //WE ASSUME P IS A POWER LAW PAST ITS ENDPOINTS
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

///////
//spherical transform
//////

/*Compute pi times the zeros of J_nu where nu is an integer.
 */
int compute_spherical_bessel_zeros(int nu, int N, double*xi){
  int i;
  for(i = 0; i < N; i++){
    //Note: GSL returns PI the root we want
    if (nu == 0){
      xi[i] = i+1.;
    }else{//nu > 0
      xi[i] = inv_PI * gsl_sf_bessel_zero_Jnu((nu+1)/2, i+1);
    }
  }
  return 0;
}

/*Compute weights of the quadrature routine.
 */
int compute_spherical_weights(int nu, int N, double*xi, double*w){
  int i;
  double x;
  for(i = 0; i < N; i++){
    if (nu == 0){
      w[i] = 1.;//gsl_sf_bessel_y0(PI*xi[i])/gsl_sf_bessel_j1(PI*xi[i]);
    }else if (nu == 1){
      x = PI*xi[i];
      w[i] = -(x*cos(x) + x*x*sin(x))/((3-x*x)*sin(x) - 3*x*cos(x));
    }else{ //nu > 1
      w[i] = gsl_sf_bessel_yl(nu, PI*xi[i])/gsl_sf_bessel_jl(nu+1, PI*xi[i]);
    }
  }
  return 0;
}

/* Compute the spherical bessel function j_nu, 
 * psi, and its derivative (dpsi) given N and h.
 */
int compute_jnu_psi_dpsi(int nu, int N, double h,
			 double*jnu, double*xi, double*psi, double*dpsi){
  int i;
  double t, PIsinht;
  double pi_over_h = PI/h;
  for(i = 0; i < N; i++){
    t = h*xi[i];
    PIsinht = PI*sinh(t);
    
    psi[i] = t * tanh(0.5*PIsinht);
    dpsi[i] = (PI*t*cosh(t) + sinh(PIsinht))/(1+cosh(PIsinht));
    if (dpsi[i]!=dpsi[i]) dpsi[i] = 1.0;//avoids NaNs

    if (nu == 0){
      jnu[i] = gsl_sf_bessel_j0(pi_over_h * psi[i]);
    }else if(nu == 1){
      jnu[i] = gsl_sf_bessel_j1(pi_over_h * psi[i]);
    }else if(nu == 2){
      jnu[i] = gsl_sf_bessel_j2(pi_over_h * psi[i]);
    }else{//nu is not 0, 1, or 2
      jnu[i] = gsl_sf_bessel_jl(nu, pi_over_h * psi[i]);
    }
  }
  return 0;
}

/* Perform the transform.
 * This specific function transforms functions of the form
 * \int_0^\inf dk/2PI P(k/r)*j_nu(k)
 * Extra factors of k must be handled externally
 * as well as the change of variables from k to k*r.
 */
int bquad_spherical_transform(double*r, int Nr, double*k, double*P, int Nk,
			      double*output, int nu, int N, double h){
  //printf("nu = %d\nN = %d\nh = %f\n", nu, N, h);
  int i, j;
  double ki, P_temp, r_inv, sum;
  double pi_over_h = PI/h;

  //WE ASSUME P IS A POWER LAW PAST ITS ENDPOINTS
  double alpha_low = log(P[1]/P[0])/log(k[1]/k[0]);
  double A_low = P[0]/pow(k[0], alpha_low);
  double alpha_hi = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
  double A_hi = P[Nk-1]/pow(k[Nk-1], alpha_hi);

  static int init_flag = 0;
  static int Nk_old = -999;
  static gsl_spline*Pspl = NULL;
  static gsl_interp_accel*acc = NULL;
  static double h_old = -1;
  static int N_max = -1;
  static double nu_old = -9999;
  static double*xi = NULL;
  static double*w = NULL;
  static double*jnu = NULL;
  static double*psi = NULL;
  static double*dpsi = NULL;

  //Create the spline and accelerator
  if ((init_flag == 0) || (Nk != Nk_old)){
    if (Pspl!=NULL) gsl_spline_free(Pspl);
    if (acc!=NULL) gsl_interp_accel_free(acc);
    Pspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
    acc = gsl_interp_accel_alloc();
    Nk_old = Nk;
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
      compute_spherical_bessel_zeros(nu, N, xi);
      compute_spherical_weights(nu, N, xi, w);
    }
    //These terms also change when h changes
    if (jnu!=NULL) free(jnu);
    if (psi!=NULL) free(psi);
    if (dpsi!=NULL) free(dpsi);
    jnu = malloc(N*sizeof(double));
    psi = malloc(N*sizeof(double));
    dpsi = malloc(N*sizeof(double));

    compute_jnu_psi_dpsi(nu, N, h, jnu, xi, psi, dpsi);
    
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
      sum += w[i]*P_temp*jnu[i]*dpsi[i];
    }
    output[j] = inv_PI_times_2 * sum;
  }

  return 0;
}
