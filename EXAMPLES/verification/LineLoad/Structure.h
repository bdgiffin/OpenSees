#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "SpatialHash.h"
#include "WindField.h"
#include <algorithm>
#include <vector>
#include <math.h>


// ======================================================================== //


// Object to map integer keys to indices
struct IndexMap {
  
  void initialize(std::vector<int>& indices_to_keys) {
    starting_key = *std::min_element(indices_to_keys.begin(),indices_to_keys.end());
    int  end_key = *std::max_element(indices_to_keys.begin(),indices_to_keys.end());
    shifted_keys_to_indices.resize(end_key-starting_key+1);
    for (int index=0; index < indices_to_keys.size(); index++) {
      shifted_keys_to_indices[indices_to_keys[index]-starting_key] = index;
    }
  } // initialize()

  int operator[](int key) { return shifted_keys_to_indices[key-starting_key]; }
  
  int starting_key;
  std::vector<int> shifted_keys_to_indices;
  
}; // IndexMap


// ======================================================================== //


// A space truss structure comprised of cylindrical members
struct Structure {
  
  // ------------------- Declare public member functions ------------------ //

  // Default empty constructor
  void Structure(void) { }
  
  // ---------------------------------------------------------------------- //

  // Legacy method to define all structural members at once
  void define_members(size_t n_members, size_t n_nodes_per_member, int *connectivity,
		      size_t n_joints, double *x_in, double *y_in, double *z_in) {

    // append new data for the current member to the stored data arrays
    for (int i=0; i < n_members; i++) {
      // shift joint indicies by 1, assuming data is coming from an Exodus file using 1-based indexing
      int id1 = connectivity[n_nodes_per_member*i+0]-1;
      int id2 = connectivity[n_nodes_per_member*i+1]-1;
      double new_coordinates[6] = { x_in[id1], y_in[id1], z_in[id1], x_in[id2], y_in[id2], z_in[id2] };
      define_member(&new_coordinates[0],i);
    } // for i=1,...,num_members
    
  } // define_member

  // ---------------------------------------------------------------------- //

  // Define a single new structural member
  void define_member(const double* new_coordinates, int new_element_tag, double new_radius = 0.0, double new_drag_coeff = 1.0) {

    // append new data for the current member to the stored data arrays
    num_members++;
    element_tag.push_back(new_element_tag);
    radius.push_back(new_radius);
    drag_coeff.push_back(new_drag_coeff);
    x1.push_back(new_coordinates[0]);
    y1.push_back(new_coordinates[1]);
    z1.push_back(new_coordinates[2]);
    x2.push_back(new_coordinates[3]);
    y2.push_back(new_coordinates[4]);
    z2.push_back(new_coordinates[5]);
    
  } // define_member

  // ---------------------------------------------------------------------- //

  // Initialize the Structure object (assuming all members have been defined)
  void initialize(void) {
    
    // initialize stored data for all members
    fx1.resize(num_members);
    fy1.resize(num_members);
    fz1.resize(num_members);
    fx2.resize(num_members);
    fy2.resize(num_members);
    fz2.resize(num_members);

    // create inverse mapping from element tags to member indices
    tag_to_index.initialize(element_tag);

    // .................................................................... //

    // get grid dimensions for spatial hashing
    double xmin = std::min(*std::min_element(x1.begin(),x1.end()),
			   *std::min_element(x2.begin(),x2.end()));
    double ymin = std::min(*std::min_element(y1.begin(),y1.end()),
			   *std::min_element(y2.begin(),y2.end()));
    double zmin = std::min(*std::min_element(z1.begin(),z1.end()),
			   *std::min_element(z2.begin(),z2.end()));
    double xmax = std::max(*std::max_element(x1.begin(),x1.end()),
			   *std::max_element(x2.begin(),x2.end()));
    double ymax = std::max(*std::max_element(y1.begin(),y1.end()),
			   *std::max_element(y2.begin(),y2.end()));
    double zmax = std::max(*std::max_element(z1.begin(),z1.end()),
			   *std::max_element(z2.begin(),z2.end()));

    // set a fixed number of grid cells on all sides
    int    Nxyz  = 100; // (this should be a user-defined parameter)
    double dedge = ((xmax - xmin) + (ymax - ymin) + (zmax - zmin))/Nxyz;

    // expand the grid dimensions by one extra layer of grid cells in all directions
    xmin -= dedge;
    ymin -= dedge;
    zmin -= dedge;
    xmax += dedge;
    ymax += dedge;
    zmax += dedge;

    // determine the number of grid cells per side
    int Nx = (xmax - xmin)/dedge;
    int Ny = (ymax - ymin)/dedge;
    int Nz = (zmax - zmin)/dedge;

    // create the spatial hash
    members_hash.set_dimensions(xmin,ymin,zmin,xmax,ymax,zmax,Nx,Ny,Nz);
    for (int i=0; i < num_members; i++) {
      members_hash.insert_segment(i,x1[i],y1[i],z1[i],x2[i],y2[i],z2[i]);
    } // for i=1,...,num_members

    // .................................................................... //
    
  } // initialize()
  
  // ---------------------------------------------------------------------- //

  // Zero forces acting on all members
  void zero_forces(void) {

    // Loop over all members and zero the applied forces
    for (int i=0; i<num_members; i++) {
      fx1[i] = 0.0; // x-force at node 1
      fy1[i] = 0.0; // y-force at node 1
      fz1[i] = 0.0; // z-force at node 1
      fx2[i] = 0.0; // x-force at node 2
      fy2[i] = 0.0; // y-force at node 2
      fz2[i] = 0.0; // z-force at node 2
    } // for i=1,...,num_members
    
  } // zero_forces()
  
  // ---------------------------------------------------------------------- //
  
  void apply_drag_forces(WindField* wind_model, double time) {

    // declare persistent static data arrays
    static std::vector<double> xmid(num_members);
    static std::vector<double> ymid(num_members);
    static std::vector<double> zmid(num_members);
    static std::vector<double> vxf(num_members);
    static std::vector<double> vyf(num_members);
    static std::vector<double> vzf(num_members);
    static std::vector<double> rhof(num_members);

    // loop over all members
    for (int i=0; i < num_members; i++) {
      
      // determine the mid-point position along the length of the current member
      xmid[i] = 0.5*(x1[i]+x2[i]);
      ymid[i] = 0.5*(y1[i]+y2[i]);
      zmid[i] = 0.5*(z1[i]+z2[i]);
    
    } // for i=0,...,num_members
    
    // determine the fluid velocity and density at the current positions of all members
    wind_model->get_fluid_velocity_and_density(num_members,time,xmid.data(),ymid.data(),zmid.data(),
					       vxf.data(),vyf.data(),vzf.data(),rhof.data());

    // loop over all members
    for (int i=0; i < num_members; i++) {
    
      // get the current element length and orientation vector
      double lambda[3] = { (x2[i]-x1[i]), (y2[i]-y1[i]), (z2[i]-z1[i]) };
      double length = std::sqrt(lambda[0]*lambda[0] + lambda[1]*lambda[1] + lambda[2]*lambda[2]);
      lambda[0] /= length;
      lambda[1] /= length;
      lambda[2] /= length;

      // get the (relative) wind velocity at the current member's mid-point location (assuming the member is stationary)
      double velocity[3] = { vxf[i], vyf[i], vzf[i] };

      // get the projected wind velocity, removing the axial component
      double axial_velocity = lambda[0]*velocity[0] + lambda[1]*velocity[1] + lambda[2]*velocity[2];
      velocity[0] -= lambda[0]*axial_velocity;
      velocity[1] -= lambda[1]*axial_velocity;
      velocity[2] -= lambda[2]*axial_velocity;
      double norm_velocity = std::sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

      // compute the total drag load
      double area = 2.0*radius[i]*length;
      double norm_force = 0.5*rhof[i]*drag_coeff[i]*area*norm_velocity;
      double drag_force[3] = { norm_force*velocity[0], norm_force*velocity[1], norm_force*velocity[2] };

      // apply the drag force evenly between the two end-points of the member
      fx1[i] += 0.5*drag_force[0]; // x-force at node 1
      fy1[i] += 0.5*drag_force[1]; // y-force at node 1
      fz1[i] += 0.5*drag_force[2]; // z-force at node 1
      fx2[i] += 0.5*drag_force[0]; // x-force at node 2
      fy2[i] += 0.5*drag_force[1]; // y-force at node 2
      fz2[i] += 0.5*drag_force[2]; // z-force at node 2
    
    } // for i=0,...,num_members
    
  } // apply_drag_forces()
  
  // ---------------------------------------------------------------------- //
  
  // find and apply contact interaction forces between all members and a single particle "p"
  void find_and_apply_contact_forces(double contact_stiff, double rp,
				     double xp, double yp, double zp,
			             double& fxp, double& fyp, double& fzp) {
    
    // find the range of grid cells that overlap with the bounding box surrounding the current particle
    int i1,j1,k1,i2,j2,k2;
    if (members_hash.find_grid_cells_overlapping_bounding_box(xp-rp,yp-rp,zp-rp,
							      xp+rp,yp+rp,zp+rp,
							      i1,j1,k1,i2,j2,k2)) {
      // loop over the full range of grid cells that overlap with the particle's bounding box
      std::set<int> segment_ids;
      for (int i=i1; i<=i2; i++) {
	for (int j=j1; j<=j2; j++) {
	  for (int k=k1; k<=k2; k++) {
	    // find the list of all segments belonging to the indicated grid index
	    int* nearby_segment_ids = nullptr;
	    int Nsegments = members_hash.find_segments_in_grid_cell(i,j,k,nearby_segment_ids);

	    // include contact interaction forces with nearby (unique) members:
	    for (int s=0; s < Nsegments; s++) segment_ids.insert(nearby_segment_ids[s]);
	  }
	}
      }
      // apply contact forces between the current particle and the found (unique) segments
      for (int s : segment_ids) {
	apply_contact_force(s,contact_stiff,rp,xp,yp,zp,fxp,fyp,fzp);
      }
    }
    
  } // find_and_apply_contact_forces()

  // ---------------------------------------------------------------------- //
  
  // compute contact interaction force between a single member "s" and a particle "p"
  void apply_contact_force(int s, double contact_stiff, double rp, double xp, double yp, double zp,
			   double& fxp, double& fyp, double& fzp) {
    
    // get the shifted coordinates of the joints of the current member
    // measured relative to the current particle's position
    double jx1 = x1[s] - xp;
    double jy1 = y1[s] - yp;
    double jz1 = z1[s] - zp;
    double jx2 = x2[s] - xp;
    double jy2 = y2[s] - yp;
    double jz2 = z2[s] - zp;

    // compute the unit tangent vector relative to the current member
    double tx = jx2 - jx1;
    double ty = jy2 - jy1;
    double tz = jz2 - jz1;
    double inv_jl2 = 1.0/(tx*tx+ty*ty+tz*tz);
    tx *= inv_jl2;
    ty *= inv_jl2;
    tz *= inv_jl2;

    // compute projected normalized coordinates of each joint on the tangent line of the member
    double xi1 = tx*jx1 + ty*jy1 + tz*jz1;
    double xi2 = tx*jx2 + ty*jy2 + tz*jz2;

    // determine if the particle's projected position on the member lies along its length
    if (xi1*xi2 >= 0.0) {
      if (abs(xi1) < abs(xi2)) {
	xi1 = 1.0;
	xi2 = 0.0;
      } else {
	xi1 = 0.0;
	xi2 = 1.0;
      }
    } else {
      xi1 = 1.0-abs(xi1);
      xi2 = 1.0-abs(xi2);
    }

    // determine the directed shortest distance from the particle to the member
    double dx = xi1*jx1 + xi2*jx2;
    double dy = xi1*jy1 + xi2*jy2;
    double dz = xi1*jz1 + xi2*jz2;
    double dl = std::sqrt(dx*dx + dy*dy + dz*dz);

    // compute the contact force acting on the particle
    double fc = contact_stiff*std::min(0.0,dl-rp)/dl;
    fxp += fc*dx;
    fyp += fc*dy;
    fzp += fc*dz;

    // compute the contact force acting on the member (distributed to the joints)
    fx1[s] -= xi1*fc*dx;
    fy1[s] -= xi1*fc*dy;
    fz1[s] -= xi1*fc*dz;
    fx2[s] -= xi2*fc*dx;
    fy2[s] -= xi2*fc*dy;
    fz2[s] -= xi2*fc*dz;
    
  } // apply_contact_force()
    
  // --------------------- Declare public data members -------------------- //

  // Common constants defined for all members
  int         num_members;  // The total number of members
  SpatialHash members_hash; // Spatial hash to help search for the candidate members that may be in contact with a given particle
  IndexMap    tag_to_index; // Mapping from (global) element tag to (local) member index

  // Data defined separately for each member
  std::vector<int>    element_tag;   // The element tags associated with all members
  std::vector<double> radius;        // The radii of all members
  std::vector<double> drag_coeff;    // The drag coefficients of all members
  std::vector<double>  x1,  y1,  z1; // The coordinates of the first  joint for all members
  std::vector<double>  x2,  y2,  z2; // The coordinates of the second joint for all members
  std::vector<double> fx1, fy1, fz1; // The forces applied to the first  joint for all members
  std::vector<double> fx2, fy2, fz2; // The forces applied to the second joint for all members
  
  // ---------------------------------------------------------------------- //
  
}; // Structure


// ======================================================================== //


#endif /* STRUCTURE_H */
