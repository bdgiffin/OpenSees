#ifndef WIND_FIELD_H
#define WIND_FIELD_H

#include "Parameters.h"
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <math.h>


// ======================================================================== //


// Abstract base class for a generic wind field model
class WindField {
public:

  // ------------------- Declare public member functions ------------------ //

  void WindField(void) {} // Empty default constructor method

  // Pure virtual method to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) = 0;
  
  // ---------------------------------------------------------------------- //

}; // WindField


// ======================================================================== //


// Example derived class for a vortex wind field (Baker & Sterling, 2017)
// https://www.sciencedirect.com/science/article/pii/S0167610517301174
class BakerSterlingVortex : public WindField {

  // ------------------- Declare public member functions ------------------ //

  // Parameterized constructor method
  void BakerSterlingVortex(const Parameters& parameters) : WindField() {
    if (parameters.count("initial_center") > 0) {
      std::vector<double> initial_center = parameters["initial_center"];
      xc = initial_center[0];
      yc = initial_center[1];
      zc = initial_center[2];
    }
  } // BakerSterlingVortex()

  // Virtual method implementation to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) {

    // Define the shifted center of the vortex at the current evaluation time
    double xct = xc + vxc*time; // [m]
    double yct = yc + vyc*time; // [m]
    double zct = zc + vzc*time; // [m]

    // Loop over all evaluation points
    for (int i=0; i<num_points; i++) {

      // Compute normalized radial and height coordinates
      double rp = std::sqrt((x[i]-xct)*(x[i]-xct) + (y[i]-yct)*(y[i]-yct));
      double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
      double cosp = (x[i]-xct)*inv_rp;
      double sinp = (y[i]-yct)*inv_rp;
      double rbar = rp/rm;
      double zbar = (z[i]-zct)/zm;

      // Compute normalized radial, tangential, and vertical velocity of the vortex
      double one_rbar2 = 1.0 + rbar*rbar;
      double one_zbar2 = 1.0 + zbar*zbar;
      double log_one_zbar2 = std::log(one_zbar2);
      double Ubar = -4.0*rbar*zbar/(one_rbar2*one_zbar2);
      double Vbar = K*std::pow(rbar,gamma-1.0)*std::pow(log_one_zbar2,0.5*gamma)/std::pow(one_rbar2,0.5*gamma);
      double Wbar = 4.0*delta*log_one_zbar2/(one_rbar2*one_rbar2);

      // Compute the x,y,z components of the fluid velocity
      vx[i] = vxc + Um*(Ubar*cosp-Vbar*sinp);
      vy[i] = vyc + Um*(Ubar*sinp+Vbar*cosp);
      vz[i] = vzc + Um*Wbar;

      // Assume the density is constant
      rhof[i] = rho0;

    } // for i=1,...,num_points

  } // get_fluid_velocity_and_density()
  
  // -------------------- Declare private data members -------------------- //

  // Define reference values and constants for use in dimensionless evaluations
  double Um = 100.0; // [m/s] reference radial velocity
  double rm = 0.1;   // [m]   reference radius
  double zm = 10.0;  // [m]   reference height
  double S = 2.0;    // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
  double K = S*(2.0/std::log(2.0));
  double gamma = 2.0;
  double delta = zm/rm;
  double rho0 = 1.293; // [kg/m^3] reference density of air at STP

  // Define the center of the vortex, and its translational velocity
  double xc  = 10.0; // [m]
  double yc  = 0.0;  // [m]
  double zc  = 0.0;  // [m]
  double vxc = 0.0;  // [m/s]
  double vyc = 0.0;  // [m/s]
  double vzc = 0.0;  // [m/s]

}; // BakerSterlingVortex


// ======================================================================== //


// Factory method to create a new WindField model from (generic) parameterized inputs
WindField* new_WindField(std::string type, const Parameters& parameters) {

  // Attempt to create a new wind field model
  if (type == "BakerSterlingVortex") {
    return new BakerSterlingVortex(parameters);
  } else {
    std::cerr << "ERROR in `new_wind_model()`; unrecognized WindField type: " << type << std::endl;
    return nullptr;
  }
  
} // new_wind_model()


// ======================================================================== //


#endif /* WIND_FIELD_H */
