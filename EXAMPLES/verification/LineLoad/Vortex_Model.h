#ifndef VORTEX_MODEL_H
#define VORTEX_MODEL_H

#include <vector>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <iostream>

class VortedModel {
public:

VortedModel(){ }

void define_tor(double Um_in, double Sm_in, double rm_in, double z_in){
  Um = Um_in;
  S = Sm_in;
  rm = rm_in;
  zm = z_in;
  
};

void position(double xw, double yw, double zw, double &vxw, double &vyw, double &vzw ,std::string type){
 double vxf, vyf, vzf, rhof;
 if (type =="Baker"){    
   get_fluid_velocity_and_density(vxf,vyf,vzf,rhof,xw,yw,zw);
 };
 vxw =vxf;
 vyw =vxf;
 vzw = vzf;
}



void get_fluid_velocity_and_density(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
    // Define reference values and constants for use in dimensionless evaluations
    double K = S*(2.0/std::log(2.0));
    double gamma = 2.0;
    double delta = zm/rm;

    // Define the center of the vortex
    double xc = 10.0; // [m]
    double yc = 0.0; // [m]
    double zc = 0.0; // [m]

    // Compute normalized radial and height coordinates
    double rp = std::sqrt((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc));
    double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
    double cosp = (xp-xc)*inv_rp;
    double sinp = (yp-yc)*inv_rp;
    double rbar = rp/rm;
    double zbar = (zp-zc)/zm;

    // Compute normalized radial, tangential, and vertical velocity of the vortex
    double one_rbar2 = 1.0 + rbar*rbar;
    double one_zbar2 = 1.0 + zbar*zbar;
    double log_one_zbar2 = std::log(one_zbar2);
    double Ubar = -4.0*rbar*zbar/(one_rbar2*one_zbar2);
    double Vbar = K*std::pow(rbar,gamma-1.0)*std::pow(log_one_zbar2,0.5*gamma)/std::pow(one_rbar2,0.5*gamma);
    double Wbar = 4.0*delta*log_one_zbar2/(one_rbar2*one_rbar2);
    double Pbar = 8.0*std::pow(rbar,2)*zbar/std::pow(one_rbar2*one_zbar2,2)-4.15*S*S*std::pow(log_one_zbar2,2)/one_rbar2 -
                  4*log_one_zbar2*one_zbar2/std::pow(one_rbar2*one_zbar2,2);
                  

    // Compute the x,y,z components of the fluid velocity
    vxf = Um*(Ubar*cosp-Vbar*sinp);
    vyf = Um*(Ubar*sinp+Vbar*cosp);
    vzf = Um*Wbar;

    // Assume the density is constant
    rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
      
  } // get_fluid_velocity_and_density()


private:
 double Um, S, rm, zm;


};

#endif