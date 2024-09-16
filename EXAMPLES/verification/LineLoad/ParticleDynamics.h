#ifndef PARTICLE_DYNAMICS_H
#define PARTICLE_DYNAMICS_H

#include "WindField.h"
#include "Structure.h"
#include "Particles.h"
#include "Parameters.h"
#include "Logger.h"
#include <string>
#include <vector>
#include <set>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <iostream>


// ======================================================================== //


// Main simulation driver object
struct ParticleDynamics {

  // ------------------- Declare public member functions ------------------ //
  
  // Default constructor
  ParticleDynamics(void) { }

  bool initialized(void) { return is_initialized; }

  // Initialize the particle dynamics simulation
  void initialize() {

    DEBUG(std::cout << "Initializing ParticleDynamics object" << std::endl;)

    // Initialize the global gravitational acceleration constant
    gz = -9.8; // [m/s^2]
    
    // Set the starting analysis time to zero
    time = 0.0; // [s]

    // Initialize all model sub-components:

    // Set WindField model parameters
    //std::string wind_model_type = "BakerSterlingVortex";
    //Parameters wind_model_parameters;
    //wind_model_parameters["initial_center"]  = std::vector<double>({ 10.0, 0.0, 0.0 }); // [m]
    //wind_model_parameters["center_velocity"] = std::vector<double>({  0.0, 0.0, 0.0 });  // [m/s]
    
    // Initialize the WindField model
    //wind_model = new_WindField(wind_model_type,wind_model_parameters);
    wind_model->initialize();

    // Initialize the debris Particles
    debris.initialize(wind_model);

    // Initialize the members in the Structure
    members.initialize();

    // Set the "is_initialized" flag
    is_initialized = true;
    
  } // initialize()
  
  // ---------------------------------------------------------------------- //
  
  // Update the simulation state to the indicated analysis time
  void update_state(double time_in) {

    // conditionally initialize the simulation state
    if (!is_initialized) initialize();

    // conditionally update the simulation state to the new analysis time
    if (time_in > time) {
      DEBUG(std::cout << "  Updating ParticleDynamics to new requested time state: " << time_in << " (old time: " << time << ")" << std::endl;)
    
      // determine the current time increment (the difference between the old and new times)
      double t0 = time;
      double time_increment = time_in - time;

      // determine the maximum allowable (stable) time step size
      double dt_max = debris.stable_time_step();

      // subdivide the increment into sub-steps based upon the stable time step size
      int Nsub_steps = std::ceil(time_increment/dt_max);
      double dt_sub = time_increment/Nsub_steps;

      // loop over sub-steps and take time steps
      DEBUG(std::cout << "  ParticleDynamics - number of sub-steps for the current time increment: " << Nsub_steps << std::endl;)
      for (int i=0; i < Nsub_steps; i++) {
	time_step(t0 + (i+1)*dt_sub);
      }

      // store the new analysis time
      time = time_in;

    }
    
  } // update_state()

  // ---------------------------------------------------------------------- //

private:
  
  // ------------------- Declare private member functions ----------------- //

  // Take a single time step to update the simulation state to the specified analysis time
  void time_step(double new_time) {
    
    // determine the current time step size, and update the current time
    double dt = new_time - time;
    time = new_time;

    // zero-initialize the forces acting on the Particles and the Structure for the current time step
    debris.zero_forces();
    members.zero_forces();

    // apply gravitational forces to the Particles
    debris.apply_gravitational_forces(0.0,0.0,gz);

    // apply drag forces to the Particles and the Structure
    debris.apply_drag_forces(wind_model,time);
    members.apply_drag_forces(wind_model,time);

    // apply contact forces between the Particles and the Structure
    for (int i=0; i<debris.num_particles; i++) {
      members.find_and_apply_contact_forces(debris.contact_stiff,debris.radius[i],
    					    debris.x[i],debris.y[i],debris.z[i],
    				            debris.fx[i],debris.fy[i],debris.fz[i]);
    }

    // apply boundary forces to the Particles
    debris.apply_boundary_forces(0.0,0.0,0.0,0.0,0.0,1.0);

    // update the positions and velocities of the Particles
    debris.integrate_equations_of_motion(dt);
    
  } // time_step()

  // ---------------------------------------------------------------------- //

public:
    
  // --------------------- Declare public data members -------------------- //

  double time; // Current analysis time
  double gz;   // Gravitational acceleration constant
  bool   is_initialized = false; // Initialization flag

  WindField* wind_model; // Wind field model interface
  Particles  debris;     // Compact debris represented as spherical particles
  Structure  members;    // Structural members represented as cylindrical rods
  
}; // ParticleDynamics


// ======================================================================== //


#endif /* PARTICLE_DYNAMICS_H */
