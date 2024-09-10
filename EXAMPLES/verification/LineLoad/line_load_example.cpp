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

// Apply a simple drag load to a given structural member
void apply_drag_load(double* forces, const double* coordinates, int element_tag, double radius, double time) {

    // get the mid-point position along the length of the element
    double icrd[3] = { 0.5*(coordinates[0]+coordinates[3]),
		       0.5*(coordinates[1]+coordinates[4]),
		       0.5*(coordinates[2]+coordinates[5]) };
    
    // get the current element length and orientation vector
    double lambda[3] = { coordinates[3]-coordinates[0],
                         coordinates[4]-coordinates[1],
		         coordinates[5]-coordinates[2] };
    double length = std::sqrt(lambda[0]*lambda[0] + lambda[1]*lambda[1] + lambda[2]*lambda[2]);
    lambda[0] /= length;
    lambda[1] /= length;
    lambda[2] /= length;

    // -------------------- EXAMPLE DRAG FORCE CALCULATION -------------------- //

    // set the (hard-coded for this example) ambient velocity and density at the current location
    double drag_coeff = 1.0;
    double area = 2.0*radius;
    double density = 0.001;
    double v_ref = 1.0;
    double time_scaling = 1.0 - std::exp(-time);
    double velocity[3];
    velocity[0] = v_ref*time_scaling*(-icrd[1]);
    velocity[1] = v_ref*time_scaling*(+icrd[0]);
    velocity[2] = v_ref*time_scaling*(1.0 - std::exp(-icrd[2]));

    // get the projected wind velocity, removing the axial component
    double axial_velocity = lambda[0]*velocity[0] + lambda[1]*velocity[1] + lambda[2]*velocity[2];
    velocity[0] -= lambda[0]*axial_velocity;
    velocity[1] -= lambda[1]*axial_velocity;
    velocity[2] -= lambda[2]*axial_velocity;
    double norm_velocity = std::sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

    // compute the total drag load
    double norm_force = 0.5*density*drag_coeff*area*length*area*norm_velocity;
    double drag_force[3] = { norm_force*velocity[0], norm_force*velocity[1], norm_force*velocity[2] };

    // apply the drag force evenly between the two end-points of the member
    forces[0] = 0.5*drag_force[0]; // x-force at node 1
    forces[1] = 0.5*drag_force[1]; // y-force at node 1
    forces[2] = 0.5*drag_force[2]; // z-force at node 1
    forces[3] = 0.5*drag_force[0]; // x-force at node 2
    forces[4] = 0.5*drag_force[1]; // y-force at node 2
    forces[5] = 0.5*drag_force[2]; // z-force at node 2
		
    // ------------------------------------------------------------------------ //
  
} // applyDragLoad()

// ======================================================================== //
