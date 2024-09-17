# Module for OpenSees
from openseespy.opensees import *

# Other needed Python packages
import numpy as np
import ParticleDynamics

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Define randomized spherical particle parameters
n_particles = 1000
particle_density         =  0.5 # [kg/m^3] (roughly the density of wood)
particle_min_diameter    = 0.01 # [m]
particle_diameter_range  =  1.0 # [m]
particle_cylinder_radius = 10.0 # [m]
particle_cylinder_height = 40.0 # [m]
particle_cylinder_center = [10.0,0.0,0.0] # [m,m,m]
random_seed = 1
ParticleDynamics.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)

# Create the parameterized wind field model (Baker Sterling Vortex)
wind_field_params = np.zeros(12)
wind_field_params[0]  = 100.0 # [m/s]      Um: reference radial velocity
wind_field_params[1]  = 1.0   # [m]        rm: reference radius
wind_field_params[2]  = 4.0   # [m]        zm: reference height
wind_field_params[3]  = 2.0   #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0   #         gamma: 
wind_field_params[5]  = 1.293 # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = 10.0  # [m]       xc0: x-position of the vortex center
wind_field_params[7]  = 0.0   # [m]       yc0: y-position of the vortex center
wind_field_params[8]  = 0.0   # [m]       zc0: z-position of the vortex center
wind_field_params[9]  = 0.0   # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0   # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0   # [m/s]     vzc: z-velocity of the vortex center
ParticleDynamics.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)

# Create a Rankine vortex
#wind_field_params = np.zeros(12)
#wind_field_params[0]  = 100.0 # [m/s]      Um: reference tangential velocity
#wind_field_params[1]  = 10.0  # [m]        rm: reference outer radius
#wind_field_params[2]  = 1.0   # [m]        rc: reference core radius
#wind_field_params[3]  = 1.0   #             E: decay index
#wind_field_params[4]  = 0.0   #       (unused) 
#wind_field_params[5]  = 1.293 # [kg/m^3] rho0: reference density of air at STP
#wind_field_params[6]  = 10.0  # [m]       xc0: x-position of the vortex center
#wind_field_params[7]  = 0.0   # [m]       yc0: y-position of the vortex center
#wind_field_params[8]  = 0.0   # [m]       zc0: z-position of the vortex center
#wind_field_params[9]  = 0.0   # [m/s]     vxc: x-velocity of the vortex center
#wind_field_params[10] = 0.0   # [m/s]     vyc: y-velocity of the vortex center
#wind_field_params[11] = 0.0   # [m/s]     vzc: z-velocity of the vortex center
#ParticleDynamics.API.define_wind_field(b"RankineVortex",wind_field_params)

# ----------------------------------
# Start of OpenSees model generation
# ----------------------------------

# --------------------------------------------------------------------------------------------------
# Example 3D linear dynamic truss structure subjected to time-varying vortex wind and debris loading
# all units are in Newtons, meters, seconds
#

# SET UP ----------------------------------------------------------------------------

wipe()				               # clear opensees model
model('basic', '-ndm', 3, '-ndf', 3)	       # 3 dimensions, 3 dof per node
# file mkdir data 			       # create data directory

# define GEOMETRY -------------------------------------------------------------

# nodal coordinates:
node(1, 0.0, 0.0, 0.0) # [m] (node, X, Y, Z)
node(2, 3.0, 0.0, 0.0) # [m]
node(3, 0.0, 3.0, 0.0) # [m]
node(4, 0.0, 0.0, 3.0) # [m]

# Single point constraints -- Boundary Conditions
fix(1, 1, 1, 1) # node DX DY DZ
fix(2, 0, 1, 1)
fix(3, 0, 0, 1)

# define MATERIAL -------------------------------------------------------------

# nodal masses:
mass(1, 20.0, 20.0, 20.0) # [kg] node#, Mx My Mz, Mass=Weight/g.
mass(2, 20.0, 20.0, 20.0)
mass(3, 20.0, 20.0, 20.0)
mass(4, 20.0, 20.0, 20.0)

# define materials
uniaxialMaterial("Elastic", 1, 200.0e+9) # [kg*m/s^2] modulus of elasticity of steel (200 GPa)

# Define ELEMENTS -------------------------------------------------------------

# define truss element connectivity
# (cross-sectional area 1in^2 = 0.00065m^2)
element("Truss", 1, 1, 2, 0.00065, 1) # [m^2] (Truss, TrussID, node1, node2, area, material)
element("Truss", 2, 1, 3, 0.00065, 1) # [m^2]
element("Truss", 3, 1, 4, 0.00065, 1) # [m^2]
element("Truss", 4, 2, 3, 0.00065, 1) # [m^2]
element("Truss", 5, 3, 4, 0.00065, 1) # [m^2]
element("Truss", 6, 4, 2, 0.00065, 1) # [m^2]

# define LineLoad elements
# (effective radius 2in = 0.05m)
element("LineLoad",  7, 1, 2, 0.05, ParticleDynamics.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
element("LineLoad",  8, 1, 3, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad",  9, 1, 4, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 10, 2, 3, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 11, 3, 4, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 12, 2, 4, 0.05, ParticleDynamics.library_name) # [m]

# RECORDER -------------------------------------------------------------

# output position-velocity-displacement (PVD) data
recorder('PVD', 'LineLoadTest_PVD', 'disp', 'reaction' ,'unbalancedLoad')

# DYNAMIC analysis -------------------------------------------------------------

# create TimeSeries
timeSeries("Linear", 1)

# create a plain load pattern
pattern("Plain", 1, 1)

# set damping based on first eigen mode
#freq = eigen('-fullGenLapack', 1)[0]**0.5
#dampRatio = 0.02
#rayleigh(0., 0., 0., 2*dampRatio/freq)

# create the analysis
#wipeAnalysis()			 # clear previously-define analysis parameters
constraints('Plain')    	 # how it handles boundary conditions
numberer("RCM")                  # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral')            # how to store and solve the system of equations in the analysis
algorithm('Linear')	         # use Linear algorithm for linear analysis
integrator('Newmark', 0.5, 0.25) # determine the next time step for an analysis
analysis('Transient')            # define type of analysis: time-dependent

# RUN analysis -------------------------------------------------------------

# perform the analysis
time = 0.0 # [s] starting time
dt   = 0.01 # [s] time increment
ParticleDynamics.output_state(time)
for step_id in range(1,100):
    time = time + dt
    analyze(1,dt) # apply 1 time step of size dt in the opensees analysis
    ParticleDynamics.output_state(time)

# finalize the ParticleDynamics module (close the Exodus files)
ParticleDynamics.finalize()
