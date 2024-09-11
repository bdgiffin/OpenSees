#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <math.h>


// ======================================================================== //


// A collection of particles
struct Particles {

  // ------------------- Declare public member functions ------------------ //

  // Take a single time step to update the state of all particles to the specified analysis time
  void time_step(double new_time) {
    
    // determine the current time step size, and update the current time
    double dt = new_time - time;
    time = new_time;
    
    // include other forces acting on the current particle:

    // include fluid drag force acting on the current particle:
    apply_drag_forces(time);

    // loop over all particles
    for (int i=0; i < num_particles; i++) {

      // include contact interaction forces with nearby members:
      if (num_members > 0) {
	find_and_apply_contact_forces(fx[i],fy[i],fz[i],0.5*d[i],xt,yt,zt);
      }

      // include contact force with the ground acting on the current particle:
      double zg = 0.0;
      double dz = zt - zg;

      // compute the contact force acting on the particle
      fz[i] -= stiffness*std::min(0.0,dz);
      
      // update the current particle's velocity and displacement
      double dt_inv_m = dt/m[i];
      vx[i] += dt_inv_m*fx[i];
      vy[i] += dt_inv_m*fy[i];
      vz[i] += dt_inv_m*fz[i];
      ux[i] += dt*vx[i];
      uy[i] += dt*vy[i];
      uz[i] += dt*vz[i];
      
    } // for(i=0...num_nodes)
    
  } // time_step()
  
  // ---------------------------------------------------------------------- //

  // Zero forces acting on all particles
  void zero_forces(void) {

    // Loop over all particles and zero the applied forces
    for (int i=0; i<num_particles; i++) {
      fx[i] = 0.0;
      fy[i] = 0.0;
      fz[i] = 0.0;

      // include gravitational force acting on the current particle
      fz[i] += -m[i]*g;
      
    } // for i=1,...,num_particles
    
  } // zero_forces()
  
  // ---------------------------------------------------------------------- //

  // apply drag forces
  void apply_drag_forces(WindField* wind_model, double time) {

    // declare persistent static data arrays
    static std::vector<double> vxf(num_particles);
    static std::vector<double> vyf(num_particles);
    static std::vector<double> vzf(num_particles);
    static std::vector<double> rhof(num_particles);
    
    // determine the fluid velocity and density at the current position of the particle
    wind_model->get_fluid_velocity_and_density(num_particles,time,x.data(),y.data(),z.data(),
					       vxf.data(),vyf.data(),vzf.data(),rhof.data());

    // Loop over all particles
    for (int i=0; i<num_particles; i++) {

      // determine the velocity of the particle relative to the fluid
      double vx_rel = vx[i] - vxf[i];
      double vy_rel = vy[i] - vyf[i];
      double vz_rel = vz[i] - vzf[i];
      double vmag_rel = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel);
    
      // compute and apply the drag force
      const double pi = 2.0*std::acos(0.0);
      double drag_coefficient = 0.47; // for an assumed spherical particle
      double area = 0.5*pi*rp*rp;
      double fdrag = -0.5*drag_coefficient*area*rhof[i]*vmag_rel;
      fx[i] += fdrag*vx_rel;
      fy[i] += fdrag*vy_rel;
      fz[i] += fdrag*vz_rel;

    } // for i=1,...,num_particles
    
  } // apply_drag_forces()
  
  // ---------------------------------------------------------------------- //

  

  // --------------------- Declare public data members -------------------- //

  // Common constants defined for all particles
  int num_particles; // The total number of particles
  double stiffness;  // Contact spring stiffness for all particles

  // Data defined separately for each particle
  std::vector<double> mass;       // The masses defined for all particles
  std::vector<double> diameter;   // The diameters of all particles
  std::vector<double>  x,  y,  z; // The current spatial coordinates of all particles
  std::vector<double> vx, vy, vz; // The current velocity of all particles
  std::vector<double> fx, fy, fz; // The current forces applied to all particles
  
  // ---------------------------------------------------------------------- //

}; // Particles


// ======================================================================== //


#endif /* PARTICLES_H */
