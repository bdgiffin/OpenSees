#include "ParticleDynamics.h"
#include <stdio.h>

// global instance of the particle dynamics object
ParticleDynamics particle_dynamics;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {

  // The main API function called by OpenSees to initialize the 
  void OPS_InitializeLineLoad(void) {

    // check to see if this is the first time this function is being called
    if (!particle_dynamics.initialized()) {
      // initialize the particle dynamics object and define randomized particle positions
      particle_dynamics.initialize();
    }
    
  } // OPS_InitializeLineLoad
  
  // ------------------------------------------------------------------------ //

  // The API function called by OpenSees to initialize a new LineLoad element
  void OPS_DefineLineLoadSegment(const double* coordinates, int element_tag, double radius) {
    // (input)  coordinates[2*3]: the nodal coordinates of the current element
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  radius:           the effective radius of the current element

    // check to see if this element was defined previously
    if (!particle_dynamics.member_exists(element_tag)) {
      particle_dynamics.define_cylindrical_member(coordinates, element_tag, radius);
    }
    
  } // OPS_DefineLineLoadSegment
  
  // ------------------------------------------------------------------------ //

  // The main API function called by OpenSees to apply loads to the current LineLoad element at a requested analysis time
  void OPS_ApplyLineLoad(double* forces, const double* coordinates, int element_tag, double radius, double time) {
    // (output) forces[2*3]:      the forces applied to the nodes of the current element
    // (input)  coordinates[2*3]: the nodal coordinates of the current element
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  radius:           the effective radius of the current element
    // (input)  time:             the current analysis time

    // conditionally update the simulation state to the indicated analysis time
    particle_dynamics.update_state(time);

    // apply drag load to the current member
    particle_dynamics.apply_drag_load(forces,coordinates,element_tag,radius,time);
    
  } // OPS_ApplyLineLoad
  
} // extern "C"

// ======================================================================== //
