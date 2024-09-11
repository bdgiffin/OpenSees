#ifndef WIND_FIELD_H
#define WIND_FIELD_H

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

  void BakerSterlingVortex(void) : WindField() {} // Parameterized constructor method

  // Virtual method implementation to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) {
    
    // Define reference values and constants for use in dimensionless evaluations
    double Um = 100.0; // [m/s] reference radial velocity
    double rm = 0.1;   // [m] reference radius
    double zm = 10.0;  // [m] reference height
    double S = 2.0;    // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
    double K = S*(2.0/std::log(2.0));
    double gamma = 2.0;
    double delta = zm/rm;

    // Define the center of the vortex
    double xc = 10.0; // [m]
    double yc = 0.0;  // [m]
    double zc = 0.0;  // [m]

    // Loop over all evaluation points
    for (int i=0; i<num_points; i++) {

      // Compute normalized radial and height coordinates
      double rp = std::sqrt((x[i]-xc)*(x[i]-xc) + (y[i]-yc)*(y[i]-yc));
      double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
      double cosp = (x[i]-xc)*inv_rp;
      double sinp = (y[i]-yc)*inv_rp;
      double rbar = rp/rm;
      double zbar = (z[i]-zc)/zm;

      // Compute normalized radial, tangential, and vertical velocity of the vortex
      double one_rbar2 = 1.0 + rbar*rbar;
      double one_zbar2 = 1.0 + zbar*zbar;
      double log_one_zbar2 = std::log(one_zbar2);
      double Ubar = -4.0*rbar*zbar/(one_rbar2*one_zbar2);
      double Vbar = K*std::pow(rbar,gamma-1.0)*std::pow(log_one_zbar2,0.5*gamma)/std::pow(one_rbar2,0.5*gamma);
      double Wbar = 4.0*delta*log_one_zbar2/(one_rbar2*one_rbar2);

      // Compute the x,y,z components of the fluid velocity
      vx[i] = Um*(Ubar*cosp-Vbar*sinp);
      vy[i] = Um*(Ubar*sinp+Vbar*cosp);
      vz[i] = Um*Wbar;

      // Assume the density is constant
      rhof[i] = 1.293; // [kg/m^3] fluid density (of air at STP)
      
    } // for i=1,...,num_points
    
  } // get_fluid_velocity_and_density()
  
  // ---------------------------------------------------------------------- //

}; // BakerSterlingVortex


// ======================================================================== //


#endif /* WIND_FIELD_H */
