#ifndef PARTICLE_DYNAMICS_H
#define PARTICLE_DYNAMICS_H

#include "SpatialHash.h"
#include <vector>
#include <set>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <iostream>

// ======================================================================== //

// Abstract base class for a generic wind field model
class WindField {
public:

  // ------------------- Declare public member functions ------------------ //

  
  void WindField(void) {} // Empty default constructor method
  
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) = 0;
  
  // ---------------------------------------------------------------------- //

}; // WindField

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

// Example derived class for a straight-line wind field
class StraightLineWindField : public WindField {

  // ------------------- Declare public member functions ------------------ //

  void StraightLineWindField(void) : WindField() {} // Parameterized constructor method
  
  // ---------------------------------------------------------------------- //

}; // StraightLineWindField

// ======================================================================== //

// A space truss structure comprised of cylindrical members
struct Elements {
  
  // ------------------- Declare public member functions ------------------ //
  
  void apply_drag_loads(WindField* wind_model) {
    
  } // apply_drag_loads()
    
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
  
}; // Elements

// ======================================================================== //

// A collection of particles
struct Particles {

  // ------------------- Declare public member functions ------------------ //

  void apply_drag_loads(WindField* wind_model) {
    
  } // apply_drag_loads()

  // --------------------- Declare public data members -------------------- //

  // Common constants defined for all particles
  int num_particles; // The total number of particles
  double stiffness;  // Contact spring stiffness for all particles

  // Data defined separately for each particle
  std::vector<double> mass;       // The masses defined for all particles
  std::vector<double> diameter;   // The diameters of all particles
  std::vector<double>  x,  y,  z; // The current spatial coordinates of all particles
  std::vector<double> ux, uy, uz; // The total displacement of all particles
  std::vector<double> vx, vy, vz; // The current velocity of all particles
  std::vector<double> fx, fy, fz; // The current forces applied to all particles
  
  // ---------------------------------------------------------------------- //

}; // Particles

// ======================================================================== //

// Main simulation driver object
class ParticleDynamics {
public:

  // ------------------- Declare public member functions ------------------ //
  
  // Default constructor
  ParticleDynamics(void) { }

  // Initialize the particle dynamics simulation
  void initialize() {
    // Set the starting analysis time to zero
    time = 0.0; // [s]
    dt_max = std::numeric_limits<double>::max();

    // Initialize the vortex model parameters
    // (need to put something here)

    // Initialize the gravitational constant and the particle's contact spring stiffness parameter
    g         = 9.8;    // [m/s^2]
    stiffness = 1.0e+2; // [N/m] = [kg/s^2]

    // Initialize the total number of nodes (debris particles), structural members (elements), and joint (element nodes)
    num_nodes   = 0;
    num_members = 0;
    num_joints  = 0;

    // Initialize the number and randomized location of all compact particles
    define_particles();
    
  } // initialize()
  
  // ---------------------------------------------------------------------- //

  // Define and initialize the list of all compact particles
  void define_particles(size_t n_particles, double *m_in, double *d_in, double *x_in, double *y_in, double *z_in) {
    // initialize stored data for all particles
    num_nodes = n_particles;
    x.resize(num_nodes);
    y.resize(num_nodes);
    z.resize(num_nodes);
    
    m.resize(num_nodes);
    d.resize(num_nodes);
    vx.resize(num_nodes);
    vy.resize(num_nodes);
    vz.resize(num_nodes);
    fx.resize(num_nodes);
    fy.resize(num_nodes);
    fz.resize(num_nodes);
    
    ux.resize(num_nodes);
    uy.resize(num_nodes);
    uz.resize(num_nodes);

    for (int i=0; i<num_nodes; i++) {
      x[i] = x_in[i];
      y[i] = y_in[i];
      z[i] = z_in[i];

      m[i] = m_in[i];
      d[i] = d_in[i];
      vx[i] = 0.0;
      vy[i] = 0.0;
      vz[i] = 0.0;
      fx[i] = 0.0;
      fy[i] = 0.0;
      fz[i] = 0.0;
      
      ux[i] = 0.0;
      uy[i] = 0.0;
      uz[i] = 0.0;

      // initialize the particle's velocity as a specified fraction of the surrounding fluid velocity
      double fraction = 0.8;
      double rhof;
      get_fluid_velocity_and_density(vx[i],vy[i],vz[i],rhof,x[i],y[i],z[i]);
      vx[i] *= fraction;
      vy[i] *= fraction;
      vz[i] *= fraction;

      // restrict the maximum stable time step
      double frequency = std::sqrt(stiffness/m[i]); // [1/s]
      dt_max = std::min(dt_max,0.5/frequency); // [s]
    }
  } // define_particles()

  // ---------------------------------------------------------------------- //

  // Define and initialize a single new structural member
  void define_cylindrical_member(const double* coordinates, int element_tag, double radius) {

  } // 

  // ---------------------------------------------------------------------- //
  
  // Define and initialize the list of all structural members
  void define_members(size_t n_members, size_t n_nodes_per_member, int *connectivity,
		      size_t n_joints, double *jx_in, double *jy_in, double *jz_in) {
    // initialize stored data for all members
    num_members = n_members;
    id1.resize(num_members);
    id2.resize(num_members);
    for (int i=0; i < num_members; i++) {
      // shift joint indicies by 1, assuming data is coming from an Exodus file using 1-based indexing
      id1[i] = connectivity[n_nodes_per_member*i+0]-1;
      id2[i] = connectivity[n_nodes_per_member*i+1]-1;
    }
    
    // initialize stored data for all joints
    num_joints  = n_joints;
    jx.resize(num_joints);
    jy.resize(num_joints);
    jz.resize(num_joints);
    for (int i=0; i < num_joints; i++) {
      jx[i] = jx_in[i];
      jy[i] = jy_in[i];
      jz[i] = jz_in[i];
    }

    // get grid dimensions for spatial hashing
    double xmin = +std::numeric_limits<double>::max();
    double ymin = +std::numeric_limits<double>::max();
    double zmin = +std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    double zmax = -std::numeric_limits<double>::max();
    for (int i=0; i < num_joints; i++) {
      xmin = std::min(xmin,jx[i]);
      ymin = std::min(ymin,jy[i]);
      zmin = std::min(zmin,jz[i]);
      xmax = std::max(xmax,jx[i]);
      ymax = std::max(ymax,jy[i]);
      zmax = std::max(zmax,jz[i]);
    }

    // set a fixed number of grid cells on all sides
    int Nxyz = 100;
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

    // create spatial hash
    members_hash.set_dimensions(xmin,ymin,zmin,xmax,ymax,zmax,Nx,Ny,Nz);
    for (int j=0; j < num_members; j++) {
      members_hash.insert_segment(j,jx[id1[j]],jy[id1[j]],jz[id1[j]],jx[id2[j]],jy[id2[j]],jz[id2[j]]);
    }
    
  } // define_members()

  // ---------------------------------------------------------------------- //
  
  // Update the simulation time state
  void update_state(double time_in, double *ux_in, double *uy_in, double *uz_in,
		                    double *vx_in, double *vy_in, double *vz_in,
		                    double *fx_in, double *fy_in, double *fz_in) {
    
    // determine the current time increment, and subdivide the increment into sub-steps
    double t0 = time;
    double time_increment = time_in - time;
    int Nsub_steps = std::ceil(time_increment/dt_max);
    double dt_sub = time_increment/Nsub_steps;

    // loop over sub-steps and take time steps
    std::cout << "  Number of sub-steps for the current time increment: " << Nsub_steps << std::endl;
    for (int i=0; i < Nsub_steps; i++) {
      time_step(t0 + (i+1)*dt_sub);
    }
    
    // loop over all particles
    for (int i=0; i < num_nodes; i++) {
    
      // output the current updated displacement, velocity, and force acting on the particle
      ux_in[i] = ux[i];
      uy_in[i] = uy[i];
      uz_in[i] = uz[i];
      vx_in[i] = vx[i];
      vy_in[i] = vy[i];
      vz_in[i] = vz[i];
      fx_in[i] = fx[i];
      fy_in[i] = fy[i];
      fz_in[i] = fz[i];
      
    } // for(i=0...num_nodes)
    
  } // update_state()

  // ---------------------------------------------------------------------- //
  
  // Get the current state of the fluid
  void get_fluid_state(int Npoints, double *x_in, double *y_in, double *z_in,
		       double *vx_in, double *vy_in, double *vz_in, double *rhof_in) {
    
    // loop over all evaluation points
    for (int i=0; i < Npoints; i++) {
      // compute the current fluid velocity and density at the current point
      get_fluid_velocity_and_density(vx_in[i],vy_in[i],vz_in[i],rhof_in[i],x_in[i],y_in[i],z_in[i]);
    } // for(i=0...Npoints)
    
  } // get_fluid_state()

private:
  
  // ------------------- Declare private member functions ----------------- //

  // Take a single time step to update the simulation state to the specified analysis time
  void time_step(double new_time) {
    
    // determine the current time step size, and update the current time
    double dt = new_time - time;
    time = new_time;
    
    // loop over all particles
    for (int i=0; i < num_nodes; i++) {
      // determine the particle's position at the current time
      double xt = x[i] + ux[i];
      double yt = y[i] + uy[i];
      double zt = z[i] + uz[i];
      
      // zero the total forces acting on the current particle
      fx[i] = 0.0;
      fy[i] = 0.0;
      fz[i] = 0.0;

      // include gravitational force acting on the current particle
      fz[i] += -m[i]*g;

      // include other forces acting on the current particle:

      // include fluid drag force acting on the current particle:
      apply_particle_drag_force(fx[i],fy[i],fz[i],0.5*d[i],xt,yt,zt,vx[i],vy[i],vz[i]);

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
  
  void find_and_apply_contact_forces(double& fxp, double& fyp, double& fzp, double rp, double xp, double yp, double zp) {
    // find the range of grid cells that overlap with the bounding box surrounding the current particle
    int i1,j1,k1,i2,j2,k2;
    if (members_hash.find_grid_cells_overlapping_bounding_box(xp-rp,yp-rp,zp-rp,xp+rp,yp+rp,zp+rp,i1,j1,k1,i2,j2,k2)) {
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
      for (int segment_id : segment_ids) {
	apply_contact_force(fxp,fyp,fzp,rp,xp,yp,zp,
			    jx[id1[segment_id]],jy[id1[segment_id]],jz[id1[segment_id]],  // coordinates of joint 1
			    jx[id2[segment_id]],jy[id2[segment_id]],jz[id2[segment_id]]); // coordinates of joint 2
      }
    }
    
  } // find_and_apply_contact_forces()

  // ---------------------------------------------------------------------- //
  
  // compute contact interaction force between a compact spherical particle and a frame member
  void apply_contact_force(double& fxp, double& fyp, double& fzp, double rp, double xp, double yp, double zp,
			   double jx1, double jy1, double jz1, double jx2, double jy2, double jz2) {
    // get the shifted coordinates of the joints of the current member
    // measured relative to the current particle's position
    jx1 -= xp;
    jy1 -= yp;
    jz1 -= zp;
    jx2 -= xp;
    jy2 -= yp;
    jz2 -= zp;

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
    double fc = stiffness*std::min(0.0,dl-rp)/dl;
    fxp += fc*dx;
    fyp += fc*dy;
    fzp += fc*dz;
    
  } // apply_contact_force()

  // ---------------------------------------------------------------------- //
  
  void apply_particle_drag_force(double& fxp, double& fyp, double& fzp, double rp,
				 double xp, double yp, double zp, double vxp, double vyp, double vzp) {
    // determine the fluid velocity and density at the current position of the particle
    double vxf, vyf, vzf, rhof;
    get_fluid_velocity_and_density(vxf,vyf,vzf,rhof,xp,yp,zp);

    // determine the velocity of the particle relative to the fluid
    double vx_rel = vxp - vxf;
    double vy_rel = vyp - vyf;
    double vz_rel = vzp - vzf;
    double vmag_rel = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel);
    
    // compute and apply the drag force
    const double pi = 2.0*std::acos(0.0);
    double drag_coefficient = 0.47; // for an assumed spherical particle
    double area = 0.5*pi*rp*rp;
    double fdrag = -0.5*drag_coefficient*area*rhof*vmag_rel;
    fxp += fdrag*vx_rel;
    fyp += fdrag*vy_rel;
    fzp += fdrag*vz_rel;
    
  } // apply_particle_drag_force()

  // ---------------------------------------------------------------------- //
  
  // Implementation of vortex model of Baker and Sterling (2017)
  // https://www.sciencedirect.com/science/article/pii/S0167610517301174
  void get_fluid_velocity_and_density(double& vxf, double& vyf, double& vzf, double& rhof,
			  	      double xp, double yp, double zp) {
    // Define reference values and constants for use in dimensionless evaluations
    double Um = 100.0; // [m/s] reference radial velocity
    double rm = 0.1; // [m] reference radius
    double zm = 10.0; // [m] reference height
    double S = 2.0; // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
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

    // Compute the x,y,z components of the fluid velocity
    vxf = Um*(Ubar*cosp-Vbar*sinp);
    vyf = Um*(Ubar*sinp+Vbar*cosp);
    vzf = Um*Wbar;

    // Assume the density is constant
    rhof = 1.293; // [kg/m^3] fluid density (of air at STP)
      
  } // get_fluid_velocity_and_density()

  // -------------------- Declare private data members -------------------- //

  double time;   // Current analysis time
  double dt_max; // Maximum allowable time step size to maintain numerical stability
  double g;      // Gravitational acceleration constant

  Particles compact_debris; // Compact debris represented as spherical particles
  Elements  members;        // Structural members represented as cylindrical rods
  
}; // ParticleDynamics

// ======================================================================== //

#endif /* PARTICLE_DYNAMICS_H */
