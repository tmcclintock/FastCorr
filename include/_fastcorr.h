//int transform_at_nu(double*t, int Nt, double*x, double*Fx,
//		    int Nx, double nu, double*f);

int compute_bessel_zeros(int nu, int N, double*pi_xi);
int compute_weights(double nu, int N, double*pi_xi, double*w);
int compute_Jnu_psi_dpsi(double nu, int N, double h,
			 double*Jnu, double*xi, double*psi, double*dpsi);
int bquad_transform(double*r, int Nr, double*k, double*P, int Nk,
		    double*output, double nu, int N, double h);
int bquad_spherical_transform(double*r, int Nr, double*k, double*P, int Nk,
			      double*output, int nu, int N, double h);
