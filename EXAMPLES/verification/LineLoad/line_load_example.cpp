#include "ParticleDynamics.h"
#include <stdio.h>

// global instance of the particle dynamics object
ParticleDynamics particle_dynamics;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {

  // The main API function called by OpenSees to initialize the external module
  void OPS_InitializeLineLoad(void) {

    // check to see if this is the first time this function is being called
    if (!particle_dynamics.initialized()) {
      // initialize the particle dynamics object and define randomized particle positions
      particle_dynamics.initialize();
    }
    
  } // OPS_InitializeLineLoad
  
  // ------------------------------------------------------------------------ //

  // The API function called by OpenSees to initialize a new LineLoad element
  void OPS_DefineLineLoadSegment(int element_tag, double radius, const double* coordinates) {
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  radius:           the effective radius of the current element
    // (input)  coordinates[2*3]: the nodal coordinates of the current element

    // define a new member in the Structure (if it doesn't already exist)
    particle_dynamics.members.define_member(coordinates, element_tag, radius);
    
  } // OPS_DefineLineLoadSegment
  
  // ------------------------------------------------------------------------ //

  // The main API function called by OpenSees to apply loads to the current LineLoad element at a requested analysis time
  void OPS_ApplyLineLoad(double time, int element_tag, const double* coordinates, double* forces) {
    // (input)  time:             the current analysis time
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  coordinates[2*3]: the updated nodal coordinates of the current element
    // (output) forces[2*3]:      the forces applied to the nodes of the current element

    // conditionally update the simulation state to the indicated analysis time
    particle_dynamics.update_state(time);

    // get loads applied to the requested element whose tag is specified
    particle_dynamics.members.get_applied_forces(element_tag,coordinates,forces);
    
  } // OPS_ApplyLineLoad
  
} // extern "C"

// ======================================================================== //
