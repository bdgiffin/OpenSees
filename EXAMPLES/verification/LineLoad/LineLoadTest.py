# =============================================================================
# Module for OpenSees
# =============================================================================

# NOTE: THE LOCALLY MODIFIED VERSION OF OPENSEES WITH THE LINELOAD ELEMENT MUST BE USED
from openseespy.opensees import *
import vfo.vfo as vfo

# Module for ParticleDynamics
import ParticleDynamics

# Python package for reading/writing data in the Exodus mesh database format
# NOTE: PYEXODUS V0.1.5 NEEDS TO BE MODIFIED TO WORK CORRECTLY WITH PYTHON 3.12
# https://pypi.org/project/pyexodus/
import pyexodus

# Other needed Python packages
import sys
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser

# =============================================================================
# Units and constants
# =============================================================================

#os.system('SI')

# Independent units
m = 1  # meters
kg = 1  # kilograms
N = 1  # Newton
KN = 1000 * N  # kilonewtons
sec = 1  # seconds
  

# Dependent units
inch = 0.0254*m  # meters
kip = 4448.22*N  # Newtons (1 kip = 4448.22 N)
kips2in = 175.126836*kg/m  # kg/m (since 1 kip/in = 175.126836 kg/m)
sq_in = inch * inch  # square meters (1 in^2 = 0.00064516 m^2)
ksi = kip / sq_in  # Pascals (force per unit area) Pascals (since 1 ksi = 6.89476 MPa = 6.89476 * 10^6 Pa)
ft = inch/12  # meters (since 1 ft = 0.3048 meters)
mm = 0.001 * m  # meters (since 1 mm = 0.001 meters)

# Constants
g = 9.81 * m / (sec * sec)  # acceleration due to gravity in m/s^2
pi = math.acos(-1)  # value of pi


# =============================================================================
# Command-Line for providing exo file
# =============================================================================

# # read any input arguments the user may have provided
# parser = ArgumentParser()
# parser.add_argument("-f", "--file", dest="filename",
#                     help="input Exodus model file for the frame structure", metavar="FILE")
# parser.add_argument("-q", "--quiet",
#                     action="store_false", dest="verbose", default=True,
#                     help="don't print status messages to stdout")
# args = parser.parse_args()

# # ---------------------------------------------------------------------------- #

# # check to make sure that the user has specified an Exodus file to define the geometry of the structure
# if args.filename is None:
#     print("ERROR: a valid Exodus model file must be specified to define the frame structure's geometry.")
#     quit()

# # read data from the Exodus model file for the frame structure
# exoin = pyexodus.exodus(file=args.filename, mode='r', array_type='numpy', title=None, numDims=None, numNodes=None, numElems=None, numBlocks=None, numNodeSets=None, numSideSets=None, io_size=0, compression=None)
# x_in,y_in,z_in = exoin.get_coords() # [m] (assumed units/dimensions of structure expressed in meters)
# n_joints = len(x_in)
# connect_in,n_members,n_nodes_per_member = exoin.get_elem_connectivity(id=1)
# exoin.close()

# # do some model pre-processing: identify the joints with supports
# supports = []
# for i in range(0,n_joints):
#     if (abs(z_in[i]) < 1.0e-6):
#         supports.append(i)

# #print(connect_in[:5])


# =============================================================================
# Command-Line for reading CSV file for layer use
# =============================================================================
#read any input arguments the user may have provided
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="input csv file with points and connectivity frame structure", metavar="FILE")
parser.add_argument("-q", "--quiet",
                    action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
args = parser.parse_args()
# check to make sure that the user has specified an Exodus file to define the geometry of the structure
if args.filename is None:
    print("ERROR: a valid CSV file must be specified to define the frame structure's geometry.")
    quit()

def read_points_and_connectivity_from_txt(filename):
    points = []
    connectivity = []
    layer = []  # Will store connectivity for each layer as a list of lists
    current_layer = None  # Track the current layer
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            reading_points = False
            reading_connectivity = False

            for line in lines:
                line = line.strip()
                if line == "Points:":
                    reading_points = True
                    reading_connectivity = False
                elif line == "Connectivity:":
                    reading_points = False
                    reading_connectivity = True
                elif reading_points:
                    coords = line.split(',')
                    if len(coords) == 3:
                        x, y, z = map(float, coords)
                        points.append((x, y, z))
                elif reading_connectivity:
                    if line.startswith("Layer name:"):
                        layer_name = line.split(":")[1].strip()
                        layer.append([])
                        current_layer = len(layer) - 1  # Index of the current layer
                    else:
                        element_nodes = list(map(int, line.split()))
                        element_nodes = [node_id for node_id in element_nodes]  #Change this for TT because it sarts with 1 
                        connectivity.append(element_nodes)
                        # Append element nodes to the current layer in connectivity2
                        layer[current_layer].append(element_nodes)
    except Exception as e:
        print(f"Error reading from {filename}: {e}")
    return points, connectivity, layer

txt_filename = args.filename  #For Transmission tower change line 172 element_nodes = [node_id for node_id in element_nodes]]
points, connect_in, layer_in = read_points_and_connectivity_from_txt(txt_filename)
n_members = len(connect_in)
n_nodes_per_member = 2
x_in =[row[0] for row in points]
y_in = [row[1] for row in points]
z_in = [row[2] for row in points]
n_joints = len(x_in)

supports = []
for i in range(0,n_joints):
    if (abs(z_in[i]) < 1.0e-6):
        supports.append(i)


# =============================================================================
# Command-Line for Defining properties
# =============================================================================

# define the cross-sectional area, effective member radius, and mass density assigned to all members
modulus_of_elasticity = 200.0e+9 # [kg*m/s^2] modulus of elasticity of steel (200 GPa)
poissons_ratio        = 0.3      # dimensionless
steel_mass_density    = 7850.0   # [kg/m^3] (mass density of steel)
cross_sectional_area  = 0.00065  # (cross-sectional area 1in^2 = 0.00065m^2)
radius_of_gyration    = 0.05     # (effective radius 2in = 0.05m)

shear_modulus        = 0.5*modulus_of_elasticity/(1.0+poissons_ratio)
mass_per_unit_length = steel_mass_density*cross_sectional_area
polar_moment_of_area = cross_sectional_area*radius_of_gyration*radius_of_gyration
moment_of_area_x     = 0.5*polar_moment_of_area
moment_of_area_y     = 0.5*polar_moment_of_area
        
#lump the nodal masses to the joints of the structure
lumped_mass = np.zeros(n_joints)
for i in range(0,n_members):
   i1 = connect_in[i][0]-1
   i2 = connect_in[i][1]-1
   dx = x_in[i2] - x_in[i1]
   dy = y_in[i2] - y_in[i1]
   dz = z_in[i2] - z_in[i1]
   member_length = math.sqrt(dx*dx + dy*dy + dz*dz)
   member_mass   = steel_mass_density*cross_sectional_area*member_length
   lumped_mass[i1] = lumped_mass[i1] + 0.5*member_mass
   lumped_mass[i2] = lumped_mass[i2] + 0.5*member_mass

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Define randomized spherical particle parameters
n_particles = 100
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

# =============================================================================
# SET UP ----------------------------------------------------------------------------
# =============================================================================

wipe()				               # clear opensees model
model('basic', '-ndm', 3, '-ndf', 6)	       # 3 dimensions, 3 dof per node
# file mkdir data 			       # create data directory

# define GEOMETRY -------------------------------------------------------------

# nodal coordinates:
for i in range(0,n_joints):
    node(i+1, x_in[i], y_in[i], z_in[i]) # [m] (node, X, Y, Z)

# Single point constraints -- Boundary Conditions
for isupport in supports:
    fix(isupport+1, 1, 1, 1, 0, 0, 0) # node DX DY DZ RX RY RZ

# define MATERIAL -------------------------------------------------------------

# nodal masses: (only needed if masses are not already computed by the element)
for i in range(0,n_joints):
    mass(i+1, lumped_mass[i], lumped_mass[i], lumped_mass[i]) # [kg] node#, Mx My Mz, Mass=Weight/g.

# define materials
# =============================================================================
# Rough Elements assign
# =============================================================================
Es = 29000 * ksi  # Steel Young's Modulus
nu = 0.3  # Poisson's ratio
Gs = Es / (2 * (1 + nu))  # Torsional stiffness Modulus
J = 10  # Large torsional stiffness

# =============================================================================
#Transformation absed on the orientation
# =============================================================================
Transf = [1,2,3]
geomTransf('Linear', Transf[0], 0, 0, 1)
geomTransf('Linear', Transf[1], 0, 0, 1)
geomTransf('Linear', Transf[2], 0, 0, 1)
#uniaxialMaterial("Elastic", 1, 200.0e+9) # [kg*m/s^2] modulus of elasticity of steel (200 GPa)

# Define ELEMENTS -------------------------------------------------------------

# =============================================================================
# Defining Fiber Section
# =============================================================================

Fy = 60.0 * ksi
Es = 29000 * ksi  # Steel Young's Modulus
nu = 0.3
Bs = 0.01
R0 = 18
cR1 = 0.925
cR2 = 0.15
matIDhard = 2
matType = 'Steel02'
# Function to define uniaxial material in Python
uniaxialMaterial(matType, matIDhard, Fy, Es, Bs, R0, cR1, cR2)

#Properties of L-Section
#Main Lega L150*150*14
#Unit Wt = 310*N/m^3
#Area = 4004 mm^2
# Radius of Gyration = 46.308*mm
# Ix = Iy = 845.4*cm^4
 # Wx=Wy = 78.33*cm^3
secTag = 2
BreID = 2
Lfiber = 20
Sfiber = 3
Thick = 14*mm
Length = 150*mm
Ly1= -Thick/2
Hy1= -Thick/2
Ly2= Length-Thick/2
Hy2= Thick/2    

def FiberCreation(secTag,matIDhard,Sfiber,Lfiber,Ly1,Hy1,Ly2,Hy2):
    section('Fiber',secTag,'-GJ', 1.0e10)
    patch('rect', matIDhard, Sfiber, Lfiber,Hy1,Ly1,Hy2,Ly2)
    patch('rect', matIDhard, Lfiber, Sfiber,-Ly1,Hy1,Ly2,Hy2)


    # SecTagTorsion = 4
    # uniaxialMaterial('Elastic', SecTagTorsion, 1.0e12 )

    # fib_sec_1 = [['section', 'Fiber', secTag, '-torsion', SecTagTorsion],
    #         ['patch', 'rect', matIDhard, Sfiber, Lfiber,Hy1,Ly1,Hy2,Ly2],
    #         ['patch', 'rect', matIDhard, Lfiber, Sfiber,-Ly1,Hy1,Ly2,Hy2],
    #         ]
    # opsv.fib_sec_list_to_cmds(fib_sec_1)   
    # matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
    # opsv.plot_fiber_section(fib_sec_1 , matcolor=matcolor)
    # plt.axis('equal')
    # plt.show()  
# Function to create nodes and elements in OpenSeesPy
     
FiberCreation(secTag,matIDhard,Sfiber,Lfiber,Ly1,Hy1,Ly2,Hy2)
QLsection = 310*N/pow(m,3)
QDLsection = QLsection*(Length*Thick+(Length-Thick)*Thick)
#print(len(connectivity))

Radius = 0.5 # This is supposed radius for the member to calculate wind and debris forces  
#Mass_Den = 
k = 0
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
         #print(f"Index {i+1}: {connectivity[i]}")  
         element('elasticBeamColumn', k+1, layer_in[i][j][0], layer_in[i][j][1], secTag, Transf[i])
         element("LineLoad", k+1+n_members, layer_in[i][j][0], layer_in[i][j][1], radius_of_gyration, ParticleDynamics.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
         k+=1

   
# for i in range(0,n_members):
#     #print(f"Index {i+1}: {connectivity[i]}")  
#     element('elasticBeamColumn', i+1,  int(connect_in[i][0]), int(connect_in[i][1]), secTag, Transf)
        
# # define LineLoad elements
# for i in range(0,n_members):
#     element("LineLoad", i+1+n_members, int(connect_in[i][0]), int(connect_in[i][1]), radius_of_gyration, ParticleDynamics.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)

# =============================================================================
# RECORDER -------------------------------------------------------------
# =============================================================================

x = len(layer_in[0])
y = len(layer_in[0]) + len(layer_in[1])
z = len(layer_in[0]) + len(layer_in[1]) + len(layer_in[2])

vfo.plot_model(
    elementgroups=[
        [list(range(1, x)), list(range(x, y)), list(range(y, z))],
        ["red", "blue", "green"]
    ],
    show_nodes='yes',
    show_nodetags='no',
    show_eletags='no',
    font_size=15,
    setview='3D',
    line_width=3
)

# vfo.plot_model(show_nodes='yes', show_nodetags='no', show_eletags='no', font_size=15, setview='3D', elementgroups=None, line_width=3)  

output_directory = 'LineLoadTest_PVD'
# output position-velocity-displacement (PVD) data
recorder('PVD', output_directory, 'disp', 'reaction' ,'unbalancedLoad')

# create a matching directory in which to dump the output data
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# =============================================================================
# DYNAMIC analysis -------------------------------------------------------------
# =============================================================================

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
