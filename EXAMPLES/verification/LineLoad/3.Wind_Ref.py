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
particle_cylinder_radius = 200.0*m # [m]
particle_cylinder_height = 25.25*m # [m]
particle_cylinder_center = [70.0*m,0.0*m,0.0*m] # [m,m,m]
random_seed = 1
SWIRL.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)

# Create the parameterized wind field model (Baker Sterling Vortex)
wind_field_params = np.zeros(12)
wind_field_params[0]  = 40*m/sec # [m/s]      Um: reference radial velocity
Vm = 10*m/sec                                 #Vm : maximum Cirumfrential velocity 
wind_field_params[1]  = 100.0*m   # [m]        rm: reference radius
wind_field_params[2]  = 25.5*m  # [m]        zm: reference height
wind_field_params[3]  =  2  #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0   #         gamma: 
wind_field_params[5]  = 1.293*kg/pow(m,3) # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = particle_cylinder_center[0]  # [m]       xc0: x-position of the vortex center
wind_field_params[7]  = particle_cylinder_center[1]  # [m]       yc0: y-position of the vortex center
wind_field_params[8]  = particle_cylinder_center[2] # [m]       zc0: z-position of the vortex center
wind_field_params[9]  = 0.0*m/sec   # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0*m/sec   # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0*m/sec   # [m/s]     vzc: z-velocity of the vortex center
SWIRL.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)
xx = np.array([0.0, 0.0, 0.0, 0.0, *[3.35] * 4, *[1.52/2] * 4, *[-3.35] * 4, *[-1.52/2] * 4, *[3.35] * 4, *[1.52/2] * 4, *[-3.35] * 4, *[-1.52/2] * 4])
xy = np.array([0.0, 0.0, 0.0, 0.0, *[3.35] * 4, *[1.52/2] * 4, *[-3.35] * 4, *[-1.52/2] * 4, *[-3.35] * 4, *[-1.52/2] * 4, *[3.35] * 4, *[1.52/2] * 4])
xz = np.array([*[15.24, 18.44, 21.64, 25.25] * 9])
points1 = len(xx)
vx=np.zeros(points1,dtype='float64')
vy =np.zeros(points1,dtype='float64')
vz=np.zeros(points1,dtype='float64')
rho = np.zeros(points1,dtype='float64')
rho = rho+ wind_field_params[5]
time = 0.0 # [s] starting time
dt = 0.001
SWIRL.API.update_state(time) # required for initialization
SWIRL.output_state(time)
# for step_id in range(1,1000):
#     time = time + dt
#     SWIRL.API.update_state(time)
#     SWIRL.output_state(time)

SWIRL.API.get_wind_field_data(points1,xx,xy,xz,vx,vy,vz,rho)
radial_wind = np.sqrt(np.mean(vx)**2 + np.mean(vy)**2)
print(vx)
print(vy)
print("The intensity measure of the wind is",radial_wind)