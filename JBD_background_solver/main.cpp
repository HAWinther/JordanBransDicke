#include <iostream>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <iomanip>

//=============================================================================
//                                                                           //
// Quick code to compute the background evolution for the Jordan-Brans-Dicke //
// model with a cosmological constant.                                       //
//                                                                           //
// Requires the GSL library                                                  //
//                                                                           //
// Compile as g++ main.cpp -I/path/include/gsl -L/path/lib  -lgsl -o JBD     //
// Run as ./JBD outputfilename omega_BD Omega_m0                             //
//                                                                           //
// The variables used are:                                                   // 
//   x  = log(a)                                                             //
//   y  = log(phi / phi0)                                                    //
//   dy = dy/dx                                                              // 
//   z  = H a^3 e^y dy/dx                                                    //
//                                                                           //
// The equations solved are:                                                 //
//                                                                           //
// Hubble-equation (E(a) = H / H0):                                          //
//   E(a) = sqrt( exp(-y) * (Omegar0/a4 + Omegam0/a3 + Omegav0)              //
//                / (1.0 + dy - w/6 * dy^2)                                  //
//              )                                                            //
//        = sqrt( exp(-y) * (Omegar0/a4 + Omegam0/a3 + Omegav0)              //
//                + z^2/12.0 * (2.0w + 3.0) * exp(-2y)/a6                    //
//              ) - z/2.0 * exp(-y)/a3                                       //
//                                                                           //
//  Klein-Gordon equation:                                                   //
//    phi'' + 3Hphi' = 8piG/(2w+3) (rho_tot - 3P_tot)                        //
//  where rho_tot = rho_m + rho_v and P_tot = -rho_v  which can be written   //
//    dy/dx  = z * exp(-y)/[a3 E(a)]                                         //
//    dz/dx  = 3a3/[E(a) (2w+3)] * [Omegam0/a3 + 4*Omegav0]                  //
//                                                                           //
//=============================================================================

//=============================================================================
// List of functions defined below
//=============================================================================

double HubbleFunction(double x, double y, double dy);
double HubbleFunction_of_z(double x, double y, double z);
double HLCDM(double x);
double OmegaM(double x, double y, double dy);
double OmegaR(double x, double y, double dy);
double OmegaPhi(double x, double y, double dy);
double GeffOverG(double y);
int    ode_JBD(double x, const double y[], double dydx[], void *params);
void   solve_ode_JBD(double yi, double dyi, double *x_arr, double *y_arr, double *dy_arr);
void   find_correct_IC_using_bisection();

//=============================================================================
// Cosmological and model parameters
//=============================================================================
struct Parameters {
  double Omegam0; // kappa^2 rhom0/[3 H^2 phi0]
  double Omegar0; // kappa^2 rhor0/[3 H^2 phi0]
  double Omegav0; // kappa^2 Lambda/[3 H^2 phi0]
  double w;       // The JBD parameters
  
  double zini;    // Starting redshift for integration
  double epsilon; // Convergence criterion
  double dy0;     // Value of dy/dx(x=0)
  double yi;      // Initial condition yi(zini)

  int npts;       // Number of points to store
  double *x_arr;  // x  = log(a) array
  double *y_arr;  // y  = log(phi/phi0) array
  double *dy_arr; // dy = dlog(phi)/dloga array
} gg;

//=============================================================================
// The Hubble function in terms of x = log(a), y = log(phi/phi0) and dy = dy/dx
//=============================================================================
double HubbleFunction(double x, double y, double dy){
  return sqrt( exp(-y) * (gg.Omegar0 * exp(-4*x) + gg.Omegam0 * exp(-3*x) + gg.Omegav0) / (1.0 + dy - gg.w/6 * dy * dy) );
}

//=============================================================================
// The Hubble function in terms of x = log(a), y = log(phi/phi0) and z = H a^3 e^y dy/dx
//=============================================================================
double HubbleFunction_of_z(double x, double y, double z){
  return sqrt( exp(-y) * (gg.Omegar0 * exp(-4*x) + gg.Omegam0 * exp(-3*x) + gg.Omegav0)  + z*z/12.0 * (2.0*gg.w + 3.0) * exp(-2*y-6*x) ) - z/2.0 * exp(-y-3*x);
}

//=============================================================================
// Hubble function for LCDM with the same value of Omega_m0 and Omega_r0
//=============================================================================
double HLCDM(double x){
  return sqrt(gg.Omegar0 * exp(-4*x) + gg.Omegam0 * exp(-3*x) + 1.0 - gg.Omegar0 - gg.Omegam0);
}

//=============================================================================
// Density parameters
//=============================================================================
double OmegaM(double x, double y, double dy){
  return gg.Omegam0 * exp(-3*x) / pow( HubbleFunction(x, y, dy) , 2);
}
double OmegaR(double x, double y, double dy){
  return gg.Omegar0 * exp(-4*x) / pow( HubbleFunction(x, y, dy) , 2);
}
double OmegaPhi(double x, double y, double dy){
  return 1.0 - OmegaM(x,y,dy) - OmegaR(x,y,dy);
}

//=============================================================================
// GeffOverG = phi0/phi = e^{-y}
//=============================================================================
double GeffOverG(double y){
  return exp(-y);
}

//=============================================================================
// The ODE for the JBD scalar field
//=============================================================================
int ode_JBD(double x, const double y[], double dydx[], void *params){
  double Hubble = HubbleFunction_of_z(x, y[0], y[1]);

  dydx[0] = y[1] * exp(- y[0] - 3*x) / Hubble;
  dydx[1] = 3.0*exp(3.0*x) / (2*gg.w+3) / Hubble * (gg.Omegam0 * exp(-3*x) + 4.0 * gg.Omegav0);

  return GSL_SUCCESS;
}

//=============================================================================
// Solve the ODE for the initial condition (y, dy/dx) = (yi, dyi) 
// Stores the values found in the provided arrays
//=============================================================================
void solve_ode_JBD(double yi, double dyi, double *x_arr, double *y_arr, double *dy_arr){
  const double xini = log(1.0 / (1.0 + gg.zini));
  const double xend = log(1.0);
  const double deltax = (xend-xini)/double(gg.npts-1);

  // Set up ODE system
  gsl_odeiv2_system sys_JBD = {ode_JBD, NULL, 4, NULL};
  gsl_odeiv2_driver * ode_JBD_driver = gsl_odeiv2_driver_alloc_y_new (&sys_JBD,  gsl_odeiv2_step_rk2, 1e-8, 1e-8, 0.0);

  // Set IC
  double y[2] = {yi, HubbleFunction(xini, yi, dyi) * exp(3.0*xini + yi) * dyi }; 

  double ode_x = xini;

  // Set IC in array
  x_arr[0]  = xini;
  y_arr[0]  = y[0];
  dy_arr[0] = y[1];

  // Solve the ODE
  for(int i = 1; i < gg.npts; i++){
    double xnow = xini + i * deltax;

    int status = gsl_odeiv2_driver_apply(ode_JBD_driver, &ode_x, xnow, y);
    if(status != GSL_SUCCESS){
      printf("Error in integrating at x = %f  y = %f\n", xnow, y[0]);
      exit(1);
    }

    x_arr[i]  = xnow;
    y_arr[i]  = y[0];
    dy_arr[i] = y[1] / HubbleFunction_of_z(xnow, y[0], y[1]) * exp(-3.0*xnow - y[0]);
  }
}

//=============================================================================
// Find correct initial conditions using bisection
//=============================================================================
void find_correct_IC_using_bisection(){
  const int npts = gg.npts;
  double *x_arr  = gg.x_arr;
  double *y_arr  = gg.y_arr;
  double *dy_arr = gg.dy_arr;

  // Find the correct initial condition
  double philow  = 0.0;
  double phihigh = 1.0;
  double phinow, yi, dyi;

  int istep = 0;
  while(1){
    ++istep;

    // Current value for phii
    phinow = (philow+phihigh)/2.0;

    // Solve ODE
    yi  = gg.yi = log(phinow);
    dyi = 0.0;
    solve_ode_JBD(yi, dyi, x_arr, y_arr, dy_arr);

    // Check for convergence
    double phiphi0 = exp(y_arr[npts-1]);
    if( fabs(phiphi0 - 1.0) < gg.epsilon && 
        fabs(gg.dy0 / dy_arr[npts-1] - 1.0) < gg.epsilon && 
        fabs(gg.Omegav0 + gg.Omegam0 + gg.Omegar0 + gg.w/6.0 * gg.dy0 * gg.dy0 - gg.dy0 - 1.0) < gg.epsilon ) {
      std::cout << "Convergence found after " << istep << " iterations" << std::endl;
      break;
    }

    // Set new values of dy0 and Omegav0
    gg.dy0     = dy_arr[npts-1];
    gg.Omegav0 = 1.0 - gg.Omegam0 - gg.Omegar0 - gg.w/6.0 * gg.dy0 * gg.dy0 + gg.dy0;
   
    // Bisection step
    if(phiphi0 < 1.0){
      philow = phinow;
    } else {
      phihigh = phinow;
    }
  }
}

int main(int argc, char **argv){

  // Read parameters from file
  if(argc < 7){
    std::cout << "Run as ./JBD outputname w_BD Omegam0 Omegar0 zini npts" << std::endl;
    exit(1);
  }

  // Set parameters and allocate arrays
  std::string filename = std::string(argv[1]);
  gg.w       = atof(argv[2]);
  gg.Omegam0 = atof(argv[3]);
  gg.Omegar0 = atof(argv[4]);
  gg.zini    = atof(argv[5]);
  gg.npts    = atoi(argv[6]);;
  gg.epsilon = 1e-6;
  gg.x_arr   = new double[gg.npts];
  gg.y_arr   = new double[gg.npts];
  gg.dy_arr  = new double[gg.npts];

  std::cout << "===================" << std::endl;
  std::cout << "Parameters:        " << std::endl;
  std::cout << "===================" << std::endl;
  std::cout << "Omegam0:     " << gg.Omegam0 << std::endl;
  std::cout << "Omegar0:     " << gg.Omegar0 << std::endl;
  std::cout << "w:           " << gg.w       << std::endl;
  std::cout << "zini:        " << gg.zini    << std::endl;
  std::cout << "npts:        " << gg.npts    << std::endl;
  std::cout << "Filename:    " << filename   << std::endl;
  std::cout << std::endl;

  // Find the correct IC
  // After this is run we have the correct x,y,dy in gg.*_arr
  find_correct_IC_using_bisection();

  // Output the results to file
  std::ofstream fp(filename.c_str());
  fp << "#      log(a)             GeffOverG(a)       H(a)    [Model: JBD with w_BD = " << gg.w << "]" << std::endl;
  fp << std::setw(14) << gg.npts << std::endl;
  for(int i = 0; i < gg.npts; i++){

    fp << std::setw(14) << std::setprecision(10) << gg.x_arr[i] << "   " << 
          std::setw(14) << std::setprecision(10) << GeffOverG(gg.y_arr[i]) << "   " << 
          std::setw(14) << std::setprecision(10) << HubbleFunction(gg.x_arr[i], gg.y_arr[i], gg.dy_arr[i]) << "   " <<
          std::endl;
  }

  // Free up memory
  delete[] gg.x_arr;
  delete[] gg.y_arr;
  delete[] gg.dy_arr;
}
