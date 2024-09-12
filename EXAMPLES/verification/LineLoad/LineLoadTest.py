# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_int

# Module for OpenSees
from openseespy.opensees import *

# Note: look in the following location when replacing existing OpenSeesPy
# /opt/homebrew/lib/python3.12/site-packages/openseespy/__init__.py

# Other needed Python packages
import sys
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser

# ---------------------------------------------------------------------------- #

# Load the pre-compiled external C/C++ "shared object" libraries
libpd = CDLL("./line_load_example.so")

# Define types to convert Numpy arrays into C arrays:

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
libpd.define_particles.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
libpd.define_particles.restype  = None

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Spherical particle parameters
n_particles = 100
particle_density = 0.5 # [kg/m^3] (roughly the density of wood)
particle_min_diameter = 0.01 # [m]
particle_diameter_range = 1.0 # [m]
particle_cylinder_radius = 10.0 # [m] (diameter of cylinder in which particles will be randomly distributed during initialization)
particle_cylinder_height = 40.0 # [m] (height of cylinder in which particles will be randomly distributed during initialization)
particle_cylinder_center = [10,0,0] # [m,m,m] (x,y,z coordinate center of cylinder)

# Generate a random collection of particles
diameters = particle_min_diameter + particle_diameter_range*np.random.rand(n_particles) # [m] (ranging from 10-20 cm in diameter)
masses = np.zeros(n_particles)
for i in range(0,n_particles):
    iradius = 0.5*diameters[i]
    masses[i] = particle_density*(4.0/3.0)*math.pi*iradius*iradius*iradius # [kg] (assuming roughly spherical shape)

position_x = np.zeros(n_particles)
position_y = np.zeros(n_particles)
for i in range(0,n_particles):
    iradial_position = particle_cylinder_radius*np.random.rand(1)
    icircum_position = 2.0*math.pi*np.random.rand(1)
    position_x[i] = particle_cylinder_center[0] + iradial_position[0]*math.cos(icircum_position[0])
    position_y[i] = particle_cylinder_center[1] + iradial_position[0]*math.sin(icircum_position[0])
position_z = particle_cylinder_center[2] + particle_cylinder_height*np.random.rand(n_particles)

# Call particle_dynamics initialization API function
libpd.define_particles(n_particles,masses,diameters,position_x,position_y,position_z)

# ----------------------------------
# Start of OpenSees model generation
# ----------------------------------

# remove existing model
wipe()

# set modelbuilder
model('basic', '-ndm', 3, '-ndf', 3)

# create nodes
node(1,   0.0,   0.0,   0.0)
node(2, 100.0,   0.0,   0.0)
node(3,   0.0, 100.0,   0.0)
node(4,   0.0,   0.0, 100.0)

# set boundary condition
fix(1, 1, 1, 1)
fix(2, 0, 1, 1)
fix(3, 0, 0, 1)

# define materials
uniaxialMaterial("Elastic", 1, 3000.0)

# define elements
element("Truss", 1, 1, 2, 1.0, 1)
element("Truss", 2, 1, 3, 1.0, 1)
element("Truss", 3, 1, 4, 1.0, 1)
element("Truss", 4, 2, 3, 1.0, 1)
element("Truss", 5, 3, 4, 1.0, 1)
element("Truss", 6, 4, 2, 1.0, 1)

# add LineLoad elements - command: LineLoad LineLoadID node1 node2 radius lib
element("LineLoad",  7, 1, 2, 0.5, "line_load_example.so")
element("LineLoad",  8, 1, 3, 0.5, "line_load_example.so")
element("LineLoad",  9, 1, 4, 0.5, "line_load_example.so")
element("LineLoad", 10, 2, 3, 0.5, "line_load_example.so")
element("LineLoad", 11, 3, 4, 0.5, "line_load_example.so")
element("LineLoad", 12, 2, 4, 0.5, "line_load_example.so")

# create TimeSeries
timeSeries("Linear", 1)

# create a plain load pattern
pattern("Plain", 1, 1)

# Create the nodal load - command: load nodeID xForce yForce
#load(4, 100.0, -50.0)

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

# perform the analysis
analyze(10)

ux = nodeDisp(4,1)
uy = nodeDisp(4,2)
