import os
import sys
import math
import numpy as np
import time as timer
from datetime import timedelta
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import pyexodus  # https://pypi.org/project/pyexodus/ PYEXODUS V0.1.5 WITH PYTHON 3.12

from Modules.SiUnits import * #Importing SI units
# Append the location of the locally installed SWIRL package to sys.path
sys.path.append("/root/Research/SWIRL/install/package/")
import SWIRL

n_particles = 0
particle_density         =  0.5*kg/pow(m,3) # [kg/m^3] (roughly the density of wood)
particle_min_diameter    = 0.01*m # [m]
particle_diameter_range  =  1.0*m # [m]
particle_cylinder_radius = 50.0*m # [m]
particle_cylinder_height = 25.25*m # [m]
particle_cylinder_center = [0.0*m,0.0*m,0.0*m] # [m,m,m]
random_seed = 1
SWIRL.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)

# Create the parameterized wind field model (Baker Sterling Vortex)
wind_field_params = np.zeros(12)
wind_field_params[0]  = 40*m/sec # [m/s]      Um: reference radial velocity
Vm = 10*m/sec                                 #Vm : maximum Cirumfrential velocity 
wind_field_params[1]  = 100.0*m   # [m]        rm: reference radius
wind_field_params[2]  = 30.0*m  # [m]        zm: reference height
wind_field_params[3]  =  2  #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0   #         gamma: 
wind_field_params[5]  = 1.293*kg/pow(m,3) # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = 0.0*m  # [m]       xc0: x-position of the vortex center
wind_field_params[7]  = 0.0*m   # [m]       yc0: y-position of the vortex center
wind_field_params[8]  = 0.0 *m  # [m]       zc0: z-position of the vortex center
wind_field_params[9]  = 0.0*m/sec   # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0*m/sec   # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0*m/sec   # [m/s]     vzc: z-velocity of the vortex center
SWIRL.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)
points = 10
xx =np.zeros(points,dtype='float64')
xy = np.zeros(points,dtype='float64')
xz = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype='float64')
vx=np.zeros(points,dtype='float64')
vy =np.zeros(points,dtype='float64')
vz=np.zeros(points,dtype='float64')

rho = np.zeros(points,dtype='float64')
rho = rho+1.293*kg/pow(m,3)
SWIRL.API.get_wind_field_data(points,xx,xy,xz,vx,vy,vz,rho)

print(vx)
