#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "SpatialHash.h"
#include <vector>


// ======================================================================== //


// A space truss structure comprised of cylindrical members
struct Structure {
  
  // ------------------- Declare public member functions ------------------ //
  
  void apply_drag_forces(WindField* wind_model, double time) {
    
  } // apply_drag_forces()
    
  // --------------------- Declare public data members -------------------- //

  // Common constants defined for all members
  int         num_members;  // The total number of members
  SpatialHash members_hash; // Spatial hash to help search for the candidate members that may be in contact with a given particle

  // Data defined separately for each member
  std::vector<int>    tag_to_index;  // Mapping from (global) element tag to (local) member index
  std::vector<int>    element_tag;   // The element tags associated with all members
  std::vector<double> radius;        // The radii of all members
  std::vector<double>  x1,  y1,  z1; // The coordinates of the first  joint for all members
  std::vector<double>  x2,  y2,  z2; // The coordinates of the second joint for all members
  std::vector<double> fx1, fy1, fz1; // The forces applied to the first  joint for all members
  std::vector<double> fx2, fy2, fz2; // The forces applied to the second joint for all members
  
  // ---------------------------------------------------------------------- //
  
}; // Structure


// ======================================================================== //


#endif /* STRUCTURE_H */
