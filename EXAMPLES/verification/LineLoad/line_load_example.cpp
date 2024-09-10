#include <stdio.h>
#include <math.h>
#include <iostream>
#include <set>
#include "Vortex_Model.h"
VortedModel Vm;
// global static set of all unique externally defined element tags
static std::set<int> unique_element_tags;

// Define all C API functions within the following block:
extern "C" {

  // The main API function called by the OpenSees LineLoad element type
  void OPS_ApplyLineLoad(double* forces, const double* coordinates, int element_tag, double radius, double time) {
    // (output) forces[2*3]:      the forces applied to the nodes of the current element
    // (input)  coordinates[2*3]: the nodal coordinates of the current element
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  radius:           the effective radius of the current element
    // (input)  time:             the current analysis time

    // check to see if this element is having loads applied to it for the first time,
    // and keep track of 
    if (!unique_element_tags.count(element_tag)) {
      unique_element_tags.insert(element_tag);
      std::cout << "Applying line load to new member " << element_tag << std::endl;
    } else {
      std::cout << "Applying line load to existing member " << element_tag << std::endl;
    }

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

    double Um = 74; //m/s
    double Sm = 2;
    double rm = 200; //m
    double z = 1.9; //m


    double velocity[3];
    Vm.define_tor(Um, Sm, rm, z);
    Vm.position(icrd[0], icrd[1], icrd[2], velocity[0], velocity[1], velocity[2] ,"Baker");
  

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
    
  } // OPS_ApplyLineLoad
  
} // extern "C"
