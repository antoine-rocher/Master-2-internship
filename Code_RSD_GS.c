#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
/************************************************************************************************************************************************************************/
/************************************************************************************************************************************************************************/

                        //VARIABLE GLOBALE
int size = 100;
int limit = 100;   
gsl_interp_accel *acc[25];
gsl_spline *spline[25];
double	fg,										// growth factor
			 	b,											//Bias factor
			 	smin, 									// Minimum of s range
			 	smax, 									// Maximum of s range
			 	sigv, 									// Value of Sigma_v (define in sigma_v function)
			 	sig8, 									// Value of sigma 8 (define in main function)
				kmin,   								// minimum k in the power spectrum file (define in main function)
			 	kmax, 									// minimum k in the power spectrum file (define in main function)
			 	lk, 									    // lk = kmax-kmin lenth for the integral calculation (define in main function)
			 	f1,										//Bias factor CLPT
			 	f2,										//Bias factor CLPT
				sig_shift;							
	const int val=160;												// Nb of lines in the files Xi_r_CLPT, Sigma_CLPT and V12_CLPT
struct my_f_params { double a; double b; };

/*********************************************************************************************************************************************************************/
/*******************\\ Gauss-Legendre integral quadrature\\***********************************************************************************************************/
																															
static const double x[] = {
    1.56289844215430828714e-02,    4.68716824215916316162e-02,
    7.80685828134366366918e-02,    1.09189203580061115002e-01,
    1.40203137236113973212e-01,    1.71080080538603274883e-01,
    2.01789864095735997236e-01,    2.32302481844973969643e-01,
    2.62588120371503479163e-01,    2.92617188038471964730e-01,
    3.22360343900529151720e-01,    3.51788526372421720979e-01,
    3.80872981624629956772e-01,    4.09585291678301542532e-01,
    4.37897402172031513100e-01,    4.65781649773358042251e-01,
    4.93210789208190933576e-01,    5.20158019881763056670e-01,
    5.46597012065094167460e-01,    5.72501932621381191292e-01,
    5.97847470247178721259e-01,    6.22608860203707771585e-01,
    6.46761908514129279840e-01,    6.70283015603141015784e-01,
    6.93149199355801965946e-01,    7.15338117573056446485e-01,
    7.36828089802020705530e-01,    7.57598118519707176062e-01,
    7.77627909649495475605e-01,    7.96897892390314476375e-01,
    8.15389238339176254384e-01,    8.33083879888400823522e-01,
    8.49964527879591284320e-01,    8.66014688497164623416e-01,
    8.81218679385018415547e-01,    8.95561644970726986709e-01,
    9.09029570982529690453e-01,    9.21609298145333952679e-01,
    9.33288535043079545942e-01,    9.44055870136255977955e-01,
    9.53900782925491742847e-01,    9.62813654255815527284e-01,
    9.70785775763706331929e-01,    9.77809358486918288561e-01,
    9.83877540706057015509e-01,    9.88984395242991747997e-01,
    9.93124937037443459632e-01,    9.96295134733125149166e-01,
    9.98491950639595818382e-01,    9.99713726773441233703e-01
};

static const double A[] = {
    3.12554234538633569472e-02,    3.12248842548493577326e-02,
    3.11638356962099067834e-02,    3.10723374275665165874e-02,
    3.09504788504909882337e-02,    3.07983790311525904274e-02,
    3.06161865839804484966e-02,    3.04040795264548200160e-02,
    3.01622651051691449196e-02,    2.98909795933328309169e-02,
    2.95904880599126425122e-02,    2.92610841106382766198e-02,
    2.89030896011252031353e-02,    2.85168543223950979908e-02,
    2.81027556591011733175e-02,    2.76611982207923882944e-02,
    2.71926134465768801373e-02,    2.66974591835709626611e-02,
    2.61762192395456763420e-02,    2.56294029102081160751e-02,
    2.50575444815795897034e-02,    2.44612027079570527207e-02,
    2.38409602659682059633e-02,    2.31974231852541216230e-02,
    2.25312202563362727021e-02,    2.18430024162473863146e-02,
    2.11334421125276415432e-02,    2.04032326462094327666e-02,
    1.96530874944353058650e-02,    1.88837396133749045537e-02,
    1.80959407221281166640e-02,    1.72904605683235824399e-02,
    1.64680861761452126430e-02,    1.56296210775460027242e-02,
    1.47758845274413017686e-02,    1.39077107037187726882e-02,
    1.30259478929715422855e-02,    1.21314576629794974079e-02,
    1.12251140231859771176e-02,    1.03078025748689695861e-02,
    9.38041965369445795116e-03,    8.44387146966897140266e-03,
    7.49907325546471157895e-03,    6.54694845084532276405e-03,
    5.58842800386551515727e-03,    4.62445006342211935096e-03,
    3.65596120132637518238e-03,    2.68392537155348241939e-03,
    1.70939265351810523958e-03,    7.34634490505671730396e-04
};

#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS

double 
  Gauss_Legendre_Integration2_100pts(double a, double b, double (*f)(double, void *), void *prms)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   const double *pA = &A[NUM_OF_POSITIVE_ZEROS - 1];

   for (; px >= x; pA--, px--) {
      dum = c * *px;
      integral += *pA * ( (*f)(d - dum,prms) + (*f)(d + dum,prms) );
   }

   return c * integral;
}

void Gauss_Legendre_Zeros_100pts( double zeros[] ) {
   
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   double *pz = &zeros[NUM_OF_ZEROS - 1];

   for (; px >= x; px--)  {
      *(zeros++) = - *px;
      *(pz--) = *px;
   }   
}

void Gauss_Legendre_Coefs_100pts( double coefs[] ) {

   const double *pA = &A[NUM_OF_POSITIVE_ZEROS - 1];
   double *pc = &coefs[NUM_OF_ZEROS - 1];

   for (; pA >= A; pA--)  {
      *(coefs++) =  *pA;
      *(pc--) = *pA;
   }   
}

/***********************************************************************************************************************************************************************/
/********************************\\ Interpolation Function \\***********************************************************************************************************/

double Pm(double k)
{
    if (k>=kmin && k<=kmax)  return gsl_spline_eval (spline[1], k, acc[1]);
    else return 0;
}

double Xim_interp(double s)
{
	if (s>=smin && s<=smax-1)	return gsl_spline_eval (spline[2], s, acc[2]);
    else	return 0;
}

double V12_interp(double s)
{	 
	if (s>=smin && s<=smax-1)	return gsl_spline_eval (spline[3], s, acc[3]);
    else return 0;
}

double Psiper_interp(double s)
{
	if (s>=smin && s<=smax-1) 	return gsl_spline_eval (spline[4], s, acc[4]);
    else return 0;

}

double Psipar_interp(double s)
{
	if (s>=smin && s<=smax-1)	return gsl_spline_eval (spline[5], s, acc[5]);
    else    return 0;
}

double Xi_R_CLPT(double s)
{
	if (s>=smin && s<=smax)	return	 gsl_spline_eval (spline[6], s, acc[6])
														+ f1*gsl_spline_eval (spline[7], s, acc[7])
														+ f2*gsl_spline_eval (spline[8], s, acc[8])
														+ f1*f1*gsl_spline_eval (spline[9], s, acc[9])
														+ f1*f2*gsl_spline_eval (spline[10], s, acc[10])
														+ f2*f2*gsl_spline_eval (spline[11], s, acc[11]);
    else	return 0;
}

double V12_CLPT(double s)
{
	if (s>=smin && s<=smax)	return 	gsl_spline_eval (spline[12], s, acc[12])
														+ f1*gsl_spline_eval (spline[13], s, acc[13])
														+ f2*gsl_spline_eval (spline[14], s, acc[14])
														+ f1*f1*gsl_spline_eval (spline[15], s, acc[15])
														+ f1*f2*gsl_spline_eval (spline[16], s, acc[16]);
    else return 0;
}

double Sig_par_CLPT(double s)
{
	if (s>=smin && s<= smax) 	return 	gsl_spline_eval (spline[17], s, acc[17])
														+ f1*gsl_spline_eval (spline[18], s, acc[18])
														+ f2*gsl_spline_eval (spline[19], s, acc[19])
														+ f1*f1*gsl_spline_eval (spline[20], s, acc[20]);
    else return 0;
}

double Sig_per_CLPT(double s)
{
	if (s>=smin && s<=smax)	return 	gsl_spline_eval (spline[21], s, acc[21])
														+ f1*gsl_spline_eval (spline[22], s, acc[22])
														+ f2*gsl_spline_eval (spline[23], s, acc[23])
														+ f1*f1*gsl_spline_eval (spline[24], s, acc[24]);
    else    return 0;
}

/**********************************************************************************************************************************************************************/
/*************\\ Gaussian Function\\***********************************************************************************************************************************/
     
double gaus (double x, double moy, double var2)
{
    return 1./(sqrt(2.*M_PI*var2))*exp(-pow((x-moy),2)/(2.*var2));
}

/**********************************************************************************************************************************************************************/
/*************\\Decomposition of s and r\\*****************************************************************************************************************************/
                                                                                                                              
double spar (double s, double mu_s)
{
    return s*mu_s;
}

double rperp (double s, double spar) // rp = sp (perp)
{
    return sqrt(s*s-spar*spar);
}

double r_real(double rp, double rpar)
{
    return sqrt(rp*rp+rpar*rpar);   
}

double mu_r(double rpar, double r)
{
    return rpar/r;
}

/***********************************************************************************************************************************************************************/
/***************\\ Sigma_8\\********************************************************************************************************************************************/

double fS8 (double k, void * params) 
{
  double x = k*8;                   
  return Pm(k)*k*k*(sin(x)-x*cos(x))*(sin(x)-x*cos(x))/pow(x,6);
}

double S8 (void)   
{
    double result, error;   
                   
    gsl_function F;
    F.function = &fS8; 
    
    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(size);
    gsl_integration_cquad(&F, smin, smax, 0, 1e-12, w, &result, &error, NULL);
    gsl_integration_cquad_workspace_free(w);
    
    return 9./(2.*M_PI*M_PI)*result;
}

/*********************************************************************************************************************************************************************/
/****************\\Normalisation of P_m\\*****************************************************************************************************************************/

double Pm_norm (double s)
{
	return Pm(s)/sig8;
}

/***********************************************************************************************************************************************************************/
/*****************\\ V_12(r) Function\\*********************************************************************************************************************************/

 //Decomposition of bessel function j1
double vsin (double k, void * params) 
{
  double alpha = *(double *) params;               
  return Pm_norm(k)*alpha/k; 
}

double vcos (double k, void * params) 
{
  double alpha = *(double *) params;
  return Pm_norm(k)*alpha;
}
 
double V12(double r) 
{
	if(r>=.1){
		double result1, error1, result2, error2, alpha, alpha1, err;
		if(r<1) err = 0.1;
			else err= 0.1;
			
		// Vsin calculation 
		alpha = 1./(r*r);
		gsl_function F1;
		F1.function = &vsin;
		F1.params = &alpha;   
		
		gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(size);
		gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(size);    
		gsl_integration_qawo_table *t1 = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_SINE, size);
		gsl_integration_qawf(&F1, kmin, err, limit, w1, w2, t1, &result1, &error1);
		
		// Vcos calculation    
		alpha1 = 1./r;
		gsl_function F2;
		F2.function = &vcos;
		F2.params = &alpha1;   
		 
		gsl_integration_workspace *w3 = gsl_integration_workspace_alloc(size);
		gsl_integration_workspace *w4 = gsl_integration_workspace_alloc(size);    
		gsl_integration_qawo_table *t2 = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_COSINE, size);
		gsl_integration_qawf(&F2, kmin, err, limit, w3, w4, t2, &result2, &error2);
		
		gsl_integration_workspace_free(w1);
		gsl_integration_workspace_free(w2);    
		gsl_integration_qawo_table_free(t1);  
		gsl_integration_workspace_free(w3);
		gsl_integration_workspace_free(w4);    
		gsl_integration_qawo_table_free(t2);    
		return result2-result1;
	}
	else return 0;
}

/************************************************************************************************************************************************************************/
/*************\\Sigma_12(mu,r) functions\\*******************************************************************************************************************************/

// Sigma_v
double fsigmav(double k, void * params)
{
    double alpha = *(double *) params;
    return alpha*Pm_norm(k);
}

double Sigmav(void)
{
    double result,error;
    double alpha = 1./3.;   
    
    gsl_function F;
    F.function = &fsigmav;
    F.params = &alpha;  
    
    gsl_integration_workspace *w =  gsl_integration_workspace_alloc(size);
    gsl_integration_qag(&F, kmin, kmax, 0, 1e-5, limit, 6, w, &result, &error); 
    gsl_integration_workspace_free(w);
    
    return result;
}

//	Psiper(r)  

 //Decomposition of bessel function j1
double Psipersin(double k, void * params)
{
    double alpha = *(double *) params;
    return alpha*Pm_norm(k)/(pow(k,3));
}

double Psipercos(double k, void * params)
{
    double alpha = *(double *) params;
    return alpha*Pm_norm(k)/(k*k);
}

double Psiper(double r)			// We adapt the error at lower scale to compute the integral
{
	if (r>=0.8){
		double err; 
//		for sin part
		double result1,error1; 
		double alpha = 1./pow(r,3);
//		for cos part
		double result2,error2; 
		double alpha1 = 1./(r*r);
		
		if(r<10) {									
			if (r<2.1) err = 10;
			else err= 1;}   
		else err=0.5;  
			
// 	Psipersin calculation    
		gsl_function F1;
		F1.function = &Psipersin;
		F1.params = &alpha;    
			
		gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(size);
		gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(size);    
		gsl_integration_qawo_table *t1 = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_SINE, size);
		gsl_integration_qawf(&F1, kmin, err, limit, w1, w2, t1, &result1, &error1);    
		
		// Psipercos calculation 
		gsl_function F2;
		F2.function = &Psipercos;
		F2.params = &alpha1;         
			   
		gsl_integration_workspace *w3 = gsl_integration_workspace_alloc(size);
		gsl_integration_workspace *w4 = gsl_integration_workspace_alloc(size);    
		gsl_integration_qawo_table *t2 = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_COSINE, size);
		gsl_integration_qawf(&F2, kmin, err, limit, w3, w4, t2, &result2, &error2);    
		
		gsl_integration_workspace_free(w1);
		gsl_integration_workspace_free(w2);    
		gsl_integration_qawo_table_free(t1);  
		gsl_integration_workspace_free(w3);
		gsl_integration_workspace_free(w4);    
		gsl_integration_qawo_table_free(t2);
		
		return (result1-result2);}
    else return 0;
}

//	 Psipar(r)
double psiparf(double k, void * params)
{
    double alpha = *(double *) params; 
    return alpha*Pm_norm(k)/k;
}

double Psipar(double r)
{
    double result,error;
    double alpha = 1./r;
    
    gsl_function F;
    F.function = &psiparf;
    F.params = &alpha;
    
    gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(size);
    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(size);    
    gsl_integration_qawo_table *t1 = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_SINE, size);
    gsl_integration_qawf(&F, kmin, 1e-1, limit, w1, w2, t1, &result, &error);
    
    gsl_integration_workspace_free(w1);
    gsl_integration_workspace_free(w2);
    gsl_integration_qawo_table_free(t1);
    
    return (result-2.*Psiper(r));
}

/********Sigma_12(mu,r)***************************************************************************************************************************************************/

double Sigma12 (double mu, double r, double sigv) 				// from linear (Gaussian streaming) model
{
    return 2.*(sigv-mu*mu*Psipar_interp(r)-(1.-(mu*mu))*Psiper_interp(r));
}

double Sigma12_CLPT (double mu, double r)							// from CLPT prediction
{
    return (mu*mu*Sig_par_CLPT(r)+(1.-mu*mu)*Sig_per_CLPT(r))/(1.+Xi_R_CLPT(r))-27.;
}
/***********************************************************************************************************************************************************************/
/************\\Correlation function in real space\\*********************************************************************************************************************/
              
double fXim(double k, void * params)
{
    double alpha = *(double *) params;
    double f = alpha*Pm_norm(k)*k;
    return f;
}

double Xim (double r)
{
    double result, error;
    double alpha = 1./(2.*M_PI*M_PI*r);
    
    gsl_function F;
    F.function = &fXim;
    F.params = &alpha;       
     
    gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(200);
    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(200);    
    gsl_integration_qawo_table *t = gsl_integration_qawo_table_alloc(r, lk, GSL_INTEG_SINE, size);
    gsl_integration_qawf(&F, kmin, 1e-3, 200, w1, w2, t, &result, &error);
    
    gsl_integration_workspace_free(w1);
    gsl_integration_workspace_free(w2);
    gsl_integration_qawo_table_free(t);
    return result;    
}

/***********************************************************************************************************************************************************************/
/*************\\Correlation function in z-space\\***********************************************************************************************************************/
																													
double fXis (double y, void *p)
{
    struct my_f_params * params = (struct my_f_params *)p;
    double spp = (params->a);																			// parallel component of s
    double rp = (params->b);																				// perpendicular componant of s
    double r = r_real(rp, y);                        															// norm of r
	double v = fg*b/(M_PI*M_PI)*V12_interp(r); 											// V_12(r)
	double mu_r = y/r;															
	double x = spp-y;																
	double moy = mu_r*v;																					// average of the gaussian distribution		
	double var = pow(fg,2)/(2.*M_PI*M_PI)*Sigma12(mu_r,r,sigv); 			// Variance Sigma_12
	if (var>0)	return (1.+Xim_interp(r))*gaus(x, moy, var); 					
	else				return 0;
}

double Xis (double sp, double spi)   //Compute in sperp, spar
{
	struct my_f_params params = {spi, sp};  
    double result = Gauss_Legendre_Integration2_100pts(-200, 200, &fXis, &params);
    return result -1.;
}

/*********************************************************************************************************************************************************************/
/****************\\Correlation function in z-space for CLPT prediction\\**********************************************************************************************/

double fXis_CLPT (double y, void *p)
{
    struct my_f_params * params = (struct my_f_params *)p;
    double spp = (params->a);
    double rp = (params->b);
    double r = r_real(rp, y);                        					
	double v = fg*V12_CLPT(r)/(1.+Xi_R_CLPT(r));
	double mu_r = y/r;
	double x = spp-y;
	double moy = mu_r*v;
	double var = pow(fg,2)*Sigma12_CLPT(mu_r,r);
	if (var>0)	return (1.+Xi_R_CLPT(r))*gaus(x, moy, var); 
	else				return 0;
}

double Xis_CLPT (double sp, double spi)   
{
	struct my_f_params params = {spi, sp};  
    double result = Gauss_Legendre_Integration2_100pts(-200, 200, &fXis_CLPT, &params);
    return result -1.;	
}

/**********************************************************************************************************************************************************************/
/*********\\Legendre Multipole\\***************************************************************************************************************************************/
													
double fmultipole(double mu, void * p) 
{
    struct my_f_params * params = (struct my_f_params *)p;
    double s = (params->a);
    double l = (params->b);
    double rp = s*sqrt(1-mu*mu);			// Decomposition of S in spar (=spi) and sp=rp
    double spi = s*mu;
    double Xi_s= Xis(rp,spi);
	return Xi_s*gsl_sf_legendre_Pl(l,mu);
}

double multipole(double s, double l)   
{
    double result;
    struct my_f_params params = {s, l}; 
	result = Gauss_Legendre_Integration2_100pts(-1, 1, &fmultipole, &params);
	return (2.*l+1.)/2.*result;
}

/********************************************************************************************************************************************************************/
/****************\\Legendre Multipole with CLPT prediction\\*********************************************************************************************************/
														
double fmultipole_CLPT(double mu, void * p) 
{
    struct my_f_params * params = (struct my_f_params *)p;
    double s = (params->a);
    double l = (params->b);
    double rp = s*sqrt(1.-mu*mu);
    double spi = s*mu;												// Decomposition of S in spar (=spi) and sp=rp
    double Xi_s= Xis_CLPT(rp,spi);
	return Xi_s*gsl_sf_legendre_Pl(l,mu);
}

double multipole_CLPT(double s, double l)   
{
    double result;
    struct my_f_params params = {s, l}; 
	result = Gauss_Legendre_Integration2_100pts(-1, 1, &fmultipole_CLPT, &params);
	return (2.*l+1.)/2.*result;
}
/************************************************************************************************************************************************************************/
/***********\\ Interpolation \\*****************************************************************************************************************************/

void interpole(int n, char ficher[100],int vmax)
{
	double T_x[vmax];
	double T_y[vmax];
	FILE* f;
	int i = 0;
	f =fopen(ficher, "r");
	for(i=0; i < vmax; i++) fscanf(f, "%lf %lf\n", &T_x[i], &T_y[i]);
	acc[n] = gsl_interp_accel_alloc ();
    spline[n] = gsl_spline_alloc(gsl_interp_cspline, vmax);
    gsl_spline_init (spline[n], T_x, T_y, vmax);
    fclose(f);
}

void interpole_Xi()
{
	FILE *fxi;
	int i;
	double Xi_x[val], Xi_f0[val], Xi_f1[val], Xi_f2[val], Xi_f1_2[val], Xi_f1_f2[val], Xi_f2_2[val];
	
	fxi = fopen("data/Xi_r_CLPT.dat", "r");
	for(i=0; i < val; i++) fscanf(fxi, "%lf %lf %lf %lf %lf %lf %lf\n",&Xi_x[i], &Xi_f0[i], &Xi_f1[i], &Xi_f2[i], &Xi_f1_2[i], &Xi_f1_f2[i], &Xi_f2_2[i]); 
    acc[6] = gsl_interp_accel_alloc ();
    spline[6] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[6], Xi_x, Xi_f0, val);
   	acc[7] = gsl_interp_accel_alloc ();
    spline[7] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[7], Xi_x, Xi_f1, val);
    acc[8] = gsl_interp_accel_alloc ();
    spline[8] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[8], Xi_x, Xi_f2, val);
    acc[9] = gsl_interp_accel_alloc ();
    spline[9] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[9], Xi_x, Xi_f1_2, val);
    acc[10] = gsl_interp_accel_alloc ();
    spline[10] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[10], Xi_x, Xi_f1_f2, val);
    acc[11] = gsl_interp_accel_alloc ();
    spline[11] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[11], Xi_x, Xi_f2_2, val);
    fclose(fxi);
}

void interpole_V12()
{
	FILE *fv12;
	int i;
	double V12_x[val], V12_f0[val], V12_f1[val], V12_f2[val], V12_f1_2[val], V12_f1_f2[val];
	
	fv12 = fopen("data/V_12_CLPT.dat", "r");
	for(i=0; i < val; i++) fscanf(fv12, "%lf %lf %lf %lf %lf %lf\n",&V12_x[i], &V12_f0[i], &V12_f1[i], &V12_f2[i], &V12_f1_2[i], &V12_f1_f2[i]);
    acc[12] = gsl_interp_accel_alloc ();
    spline[12] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[12], V12_x, V12_f0, val);
   	acc[13] = gsl_interp_accel_alloc ();
    spline[13] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[13], V12_x, V12_f1, val);
    acc[14] = gsl_interp_accel_alloc ();
    spline[14] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[14], V12_x, V12_f2, val);
    acc[15] = gsl_interp_accel_alloc ();
    spline[15] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[15], V12_x, V12_f1_2, val);
    acc[16] = gsl_interp_accel_alloc ();
    spline[16] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[16], V12_x, V12_f1_f2, val);
    fclose(fv12);
}
/*
void interpole_sigma()
{	
	int i, val =150;
	double S_x[val], S_par[val], S_per[val];
	FILE* fx;
   	fx = fopen("../../CLPT/Sigma.dat", "r");
	for(i=0; i < val; i++) fscanf(fx, "%lf %lf %lf\n", &S_x[i], &S_par[i], &S_per[i]);
	acc[8] = gsl_interp_accel_alloc ();
    spline[8] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[8], S_x, S_par, val);
    acc[9] = gsl_interp_accel_alloc ();
    spline[9] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[9], S_x, S_per, val);
    fclose(fx);
}*/

void interpole_sigma()
{	
	int i;
	double S_x[val], S_par_f0[val], S_par_f1[val], S_par_f2[val], S_par_f1_2[val], S_per_f0[val], S_per_f1[val], S_per_f2[val], S_per_f1_2[val];
	FILE* fsig;
	
   	fsig = fopen("data/Sigma_12_CLPT.dat", "r");
	for(i=0; i<val; i++) fscanf(fsig, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &S_x[i], &S_par_f0[i], &S_par_f1[i], &S_par_f2[i], &S_par_f1_2[i], &S_per_f0[i], &S_per_f1[i], &S_per_f2[i], &S_per_f1_2[i]);
   acc[17] = gsl_interp_accel_alloc ();
    spline[17] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[17], S_x, S_par_f0, val);
   	acc[18] = gsl_interp_accel_alloc ();
    spline[18] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[18], S_x, S_par_f1, val);
    acc[19] = gsl_interp_accel_alloc ();
    spline[19] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[19], S_x, S_par_f2, val);
    acc[20] = gsl_interp_accel_alloc ();
    spline[20] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[20], S_x, S_par_f1_2, val);
       acc[21] = gsl_interp_accel_alloc ();
    spline[21] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[21], S_x, S_per_f0, val);
   	acc[22] = gsl_interp_accel_alloc ();
    spline[22] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[22], S_x, S_per_f1, val);
    acc[23] = gsl_interp_accel_alloc ();
    spline[23] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[23], S_x, S_per_f2, val);
    acc[24] = gsl_interp_accel_alloc ();
    spline[24] = gsl_spline_alloc(gsl_interp_cspline, val);
    gsl_spline_init (spline[24], S_x, S_per_f1_2, val);
    fclose(fsig);
}

/***********************************************************************************************************************************************************************/
/************\\Count the number of ligne in a file\\********************************************************************************************************************/

int compte(FILE *fichier)
{
	int c;
	int nLignes = 0;
	int c2 = '\0';
	while((c=fgetc(fichier)) != EOF){	if(c=='\n')	nLignes++;	c2 = c;}
	if(c2 != '\n') nLignes++;  
	return nLignes;
}

/************************************************************************************************************************************************************************/
/******\\Write in a file\\***********************************************************************************************************************************************/

void write(char f[], double func(double))
{    
	double spi;
	FILE* fi;
    fi = fopen(f,"w+"); 			  
    for(spi=smin; spi<smax; spi+=0.5)	fprintf(fi,"%le %le\n", spi, func(spi));
    fclose(fi);
}

/******************************************************************************************************************************************************************/
/***************\\MAIN FUNCTION\\**********************************************************************************************************************************/
int main (int argc, char *argv[])
{	
	double spi,sp,li;
    int i, nLignes;
	FILE *fi, *par, *ps, *f;
	char file[BUFSIZ];
/*	
	sig_shift = 0;
	b = 1;
	fg = 0.820165;
	f1 = atof(argv[1]);
	f2 = atof(argv[2]);
	smin=2.5;
	smax=158;
	val=160;
//	printf("f1=%lf f2=%lf\n",f1,f2);
//	char* file = "power_spec_EBOSS.txt";*/
	par = fopen(argv[1],"r");  
	
// Check parameter file
	if (argc ==2 && par!=0) EXIT_SUCCESS;											
	else {
		printf("Wrong parameter file\n"); 
		return EXIT_FAILURE;
	}
			
// Read parameter file
	fscanf(par,"%*s %lf\n", &smin);
	fscanf(par,"%*s %lf\n", &smax);
	fscanf(par,"%*s %lf\n", &b);
	fscanf(par,"%*s %lf\n", &fg);
	fscanf(par,"%*s %lf\n", &f1);
	fscanf(par,"%*s %lf\n", &f2);
	fscanf(par,"%*s %lf\n", &sig_shift);
	fscanf(par,"%*s %s\n",file);
	
//Check smin and smax
	li = (smax-smin)/0.5;
	if (li<0){ 
		printf("wrong parameters smin/smax");
		return EXIT_FAILURE;
		}
		
// Count number of lines
	ps = fopen(file, "r"); 
	if (ps!=0){
		nLignes = compte(ps);
		} 
	else {
		printf("Error openning power spectrum file.\n");
		return 1;
	}

// Read kmin and kmax	
	rewind(ps);
	kmin=1e30;
	kmax=0;		
	for(i=0; i < nLignes; i++) {
		double va;
		fscanf(ps, "%lf %*f\n", &va);
		if (va<kmin) kmin=va;
		if (va>kmax) kmax=va;
	}
	fclose(ps);
	lk = kmax-kmin;	
// Interpolation of power spectrum file
	interpole(1,file, nLignes);    

	sig8=1.;
	sigv = Sigmav();
//	printf("s8=%lf\n",sig8);

// Compute and save functions in text files	
	write("data/Xim.dat",Xim);
    write("data/V12.dat",V12);
    write("data/Psiper.dat",Psiper);
    write("data/Psipar.dat",Psipar);

    interpole(2,"data/Xim.dat", li); 
    interpole(3,"data/V12.dat", li); 		
    interpole(4,"data/Psiper.dat", li); 
    interpole(5,"data/Psipar.dat", li); 

	interpole_Xi();
	interpole_V12();
	interpole_sigma();
	


//	Multipoles calculation 
	fi = fopen("data/multipole.dat","w+"); 	;     
	printf("Calcul des multipoles\n");
	for(spi=smin; spi<smax; spi+=1){     		                
		double xi0=sig8*multipole(spi,0); 
		double xi2=sig8*multipole(spi,2);
		double xi4=sig8*multipole(spi,4);          
		fprintf(fi,"%le %le %le %le\n", spi, xi0, xi2, xi4);
	}
	fclose(fi);		

	fi = fopen("data/multipole_CLPT.dat","w+");     
	printf("Calcul des  multipoles CLPT\n");
	for(spi=smin; spi<smax; spi+=1){     		                
		double xi0=multipole_CLPT(spi,0); 
		double xi2=multipole_CLPT(spi,2);
		double xi4=multipole_CLPT(spi,4);          
		//printf("%lf %lf %lf ", spi, xi0, xi2);
		fprintf(fi,"%lf %lf %lf %lf\n", spi, xi0, xi2,xi4);
	}
	fclose(fi);

	//Correlation function Xi_s(spar,sper)
	f = fopen("data/Xi_s.dat","w+");       
    printf("Calcul de Xis\n");	  
    for(spi=smin; spi<smax; spi+=0.1){   
   		for(sp=smin; sp<smax; sp+=0.1){                         
		 	double V=Xis(sp,spi);           
		 	fprintf(f,"%le %le %le\n",sp, spi, V);
     	}
     }
    fclose(f);

	f = fopen("data/Xi_s_CLPT.dat","w+");       
    printf("Calcul de Xis_CLPT\n");       			  
    for(spi=smin; spi<smax; spi+=0.1){   
   		for(sp=smin; sp<smax; sp+=0.1){                         
		 	double V=Xis_CLPT(sp,spi);           
		 	fprintf(f,"%le %le %le\n",sp, spi, V);
     	}
     }
    fclose(f);

    for (i=0; i<25; i++){
    	gsl_spline_free (spline[i]);
   		gsl_interp_accel_free (acc[i]);
   	}
   		
    return 0;
}



