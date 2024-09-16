# Module for OpenSees
from openseespy.opensees import *

# Other needed Python packages
import numpy as np
import ParticleDynamics

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Define randomized spherical particle parameters
n_particles = 100
particle_density         =  0.5 # [kg/m^3] (roughly the density of wood)
particle_min_diameter    = 0.01 # [m]
particle_diameter_range  =  1.0 # [m]
particle_cylinder_radius = 10.0 # [m]
particle_cylinder_height = 40.0 # [m]
particle_cylinder_center = [10,0,0] # [m,m,m]
random_seed = 2
ParticleDynamics.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)

# Create the parameterized wind field model
wind_field_params = np.zeros(12)
wind_field_params[0]  = 100.0 # [m/s]      Um: reference radial velocity
wind_field_params[1]  = 0.1   # [m]        rm: reference radius
wind_field_params[2]  = 10.0  # [m]        zm: reference height
wind_field_params[3]  = 2.0   #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0   #         gamma: 
wind_field_params[5]  = 1.293 # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = 10.0  # [m]        xc: x-position of the vortex center
wind_field_params[7]  = 0.0   # [m]        yc: y-position of the vortex center
wind_field_params[8]  = 0.0   # [m]        zc: z-position of the vortex center
wind_field_params[9]  = 0.0   # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0   # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0   # [m/s]     vzc: z-velocity of the vortex center
ParticleDynamics.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)

# ----------------------------------
# Start of OpenSees model generation
# ----------------------------------

# remove existing model
wipe()

# set modelbuilder
model('basic', '-ndm', 3, '-ndf', 3)

# create nodes
node(1, 0.0, 0.0, 0.0) # [m]
node(2, 3.0, 0.0, 0.0) # [m]
node(3, 0.0, 3.0, 0.0) # [m]
node(4, 0.0, 0.0, 3.0) # [m]

# set boundary condition
fix(1, 1, 1, 1)
fix(2, 0, 1, 1)
fix(3, 0, 0, 1)

# define materials
uniaxialMaterial("Elastic", 1, 200.0e+9) # [kg*m/s^2] modulus of elasticity of steel (200 GPa)

# define elements
element("Truss", 1, 1, 2, 0.00065, 1) # [m^2] (roughly 1in^2 = 0.00065m^2)
element("Truss", 2, 1, 3, 0.00065, 1) # [m^2]
element("Truss", 3, 1, 4, 0.00065, 1) # [m^2]
element("Truss", 4, 2, 3, 0.00065, 1) # [m^2]
element("Truss", 5, 3, 4, 0.00065, 1) # [m^2]
element("Truss", 6, 4, 2, 0.00065, 1) # [m^2]

# add LineLoad elements - command: LineLoad LineLoadID node1 node2 radius lib
element("LineLoad",  7, 1, 2, 0.05, ParticleDynamics.library_name) # [m] (roughly 2in = 0.05m)
element("LineLoad",  8, 1, 3, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad",  9, 1, 4, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 10, 2, 3, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 11, 3, 4, 0.05, ParticleDynamics.library_name) # [m]
element("LineLoad", 12, 2, 4, 0.05, ParticleDynamics.library_name) # [m]

# create TimeSeries
timeSeries("Linear", 1)

# create a plain load pattern
pattern("Plain", 1, 1)

# (loads applied to the structure should come from ParticleDynamics module)

# ------------------------------
# Start of analysis generation
# ------------------------------

# create SOE
system("BandSPD")

# create DOF number
numberer("RCM")

# create constraint handler
constraints("Plain")

# create integrator
integrator("LoadControl", 1.0)

# create algorithm
algorithm("Linear")

# create analysis object
analysis("Static")

# output position-velocity-displacement (PVD) data
recorder('PVD', 'LineLoadTest_PVD', 'disp', 'reaction' ,'unbalancedLoad')

# perform the analysis
dt = 0.1 # [s]
for time in range(1,100):
    analyze(1,dt)
    ParticleDynamics.output_state(time)

# finalize the ParticleDynamics module (close the Exodus files)
ParticleDynamics.finalize()
