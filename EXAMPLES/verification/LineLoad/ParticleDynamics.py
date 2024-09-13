# Module for calling the C/C++ ParticleDynamicsAPI functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_int, c_char_p

# Other needed Python packages
import sys
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser

# ---------------------------------------------------------------------------- #

# Module initialization:

# Load the pre-compiled external C/C++ "shared object" libraries
library_name = "./ParticleDynamicsAPI.so"
API = CDLL(library_name)

# Define types to convert Numpy arrays into C arrays:

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
API.define_wind_field.argtypes = [c_char_p, ND_POINTER_1]
API.define_wind_field.restype  = None
API.define_particles.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_particles.restype  = None
API.define_members.argtypes = [c_size_t, c_size_t, NI_POINTER_2, c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_members.restype  = None
API.get_particle_field_data.argtypes = [ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_particle_field_data.restype  = None
API.get_wind_field_data.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_wind_field_data.restype  = None

# ---------------------------------------------------------------------------- #

# Generate a randomized collection of particles with variable diameters and positions,
# and initialize the ParticleDynamics object with these particles prior to initialization
def create_random_particles(n_particles,density,min_diameter,diameter_range,cylinder_radius,cylinder_height,cylinder_center,random_seed):
    # n_particles     [int]    The total number of spherical debris particles to be defined
    # density         [kg/m^3] The constant mass density assigned to all particles
    # min_diameter    [m]      The minimum particle diameter
    # diameter_range  [m]      The range of random particle diameters
    # cylinder_radius [m]      The diameter of the cylinder in which particles will be randomly distributed
    # cylinder_height [m]      The height of the cylinder in which particles will be randomly distributed
    # cylinder_center [m,m,m]  The x,y,z coordinate center of cylinder at the time of initialization
    # random_seed     [int]    The random number generator seed value, to ensure reproducibility

    # Instantiate a random number generator (rng) initialized with the provided random_seed value
    rng = np.random.default_rng(random_seed)

    # Generate a random collection of spherical particles with variable diameters but constant density
    diameters = min_diameter + diameter_range*rng.random(n_particles) # [m]
    masses = np.zeros(n_particles)
    for i in range(0,n_particles):
        iradius = 0.5*diameters[i]
        masses[i] = density*(4.0/3.0)*math.pi*iradius*iradius*iradius # [kg] (assuming roughly spherical shape)

    # Randomize the initial positions of all particles
    position_x = np.zeros(n_particles)
    position_y = np.zeros(n_particles)
    for i in range(0,n_particles):
        iradial_position = cylinder_radius*rng.random(1)
        icircum_position = 2.0*math.pi*rng.random(1)
        position_x[i] = cylinder_center[0] + iradial_position[0]*math.cos(icircum_position[0])
        position_y[i] = cylinder_center[1] + iradial_position[0]*math.sin(icircum_position[0])
    position_z = cylinder_center[2] + cylinder_height*rng.random(n_particles)

    # Call particle_dynamics initialization API function
    API.define_particles(n_particles,masses,diameters,position_x,position_y,position_z)

# ---------------------------------------------------------------------------- #
