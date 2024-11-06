import openseespy.opensees as op
import math
import vfo.vfo as vfo
import opsvis as opsv
import matplotlib.pyplot as plt
import numpy as np

import pyexodus
import sys
import os
import math
import time as timer
from datetime import timedelta
from argparse import ArgumentParser

from Modules.SiUnits import * #Importing SI units

# read any input arguments the user may have provided
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="input Exodus model file for the frame structure", metavar="FILE")
parser.add_argument("-q", "--quiet",
                    action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
args = parser.parse_args()

# =============================================================================
# Command-Line for reading CSV file for layer
# =============================================================================
 # read any input arguments the user may have provided
if args.filename.lower().endswith('.csv'):
    file_format = "csv"
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
    x_in =[round(row[0]*m,2)for row in points]
    y_in = [round(row[1]*m,2) for row in points]
    z_in = [round(abs(row[2])*m,2)for row in points]
    n_joints = len(x_in)
    # p_array = np.array([[x_in[i], y_in[i], z_in[i]] for i in range(n_joints)])
    # p_array = p_array[p_array[:, 2].argsort()]
    
else:
    print("Provide a vaid file format csv")
    exit()


supports = []
for i in range(0,n_joints):
    if (abs(z_in[i]) < 1.0e-6):
        supports.append(i)
# =============================================================================
# Define nodes and fix supports
# =============================================================================
op.wipe()				               
op.model('basic', '-ndm', 3, '-ndf', 6)	      

for i in range(0,n_joints):
    op.node(i+1, x_in[i], y_in[i], z_in[i]) 

for isupport in supports:
    op.fix(isupport+1, 1, 1, 1, 1, 1, 1)
    
    
Thick = np.array([12*mm,12*mm,10*mm,6*mm,5*mm,5*mm,5*mm]) #should be based on number of layers
Length = np.array([130*mm,110*mm,100*mm,60*mm,45*mm,50*mm,50*mm]) #should be based on number of layers
Ly1= -Thick/2
Hy1= -Thick/2
Ly2= Length-Thick/2
Hy2= Thick/2    
Es = 200* Gpa 
nu = 0.3

QLsection = 310*N/pow(m,3)
steel_mass_density    = 7850.0*kg/(m**3)  # [kg/m^3] (mass density of steel)
shear_modulus        = 0.5*Es /(1.0+nu)
radius_of_gyration    = 0.05*m     # (effective radius 2in = 0.05m)
cross_sectional_area = np.zeros(len(layer_in))
mass_per_unit_length =np.zeros(len(layer_in))
polar_moment_of_area = np.zeros(len(layer_in))
moment_of_area_x = np.zeros(len(layer_in))
moment_of_area_y =np.zeros(len(layer_in))
lumped_mass = np.zeros(n_joints)


for i in range(len(layer_in)):
   cross_sectional_area[i]  = Length[i]*Thick[i]+(Length[i]-Thick[i])*Thick[i]  # (cross-sectional area 1in^2 = 0.00065m^2)
   mass_per_unit_length[i] = steel_mass_density*cross_sectional_area[i]
   polar_moment_of_area[i] = cross_sectional_area[i]*radius_of_gyration*radius_of_gyration
   moment_of_area_x[i]     = 0.5*polar_moment_of_area[i]
   moment_of_area_y[i]     = 0.5*polar_moment_of_area[i]


# Elastic Section Use only when not using Fiber Section
secTag = [1,2,3,4,5,6,7]
for i in range(len(layer_in)):
    op.section('Elastic', secTag[i], Es, cross_sectional_area[i], moment_of_area_x[i], moment_of_area_y[i], shear_modulus, polar_moment_of_area[i])

k = 0
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
         op.element('Truss', k+1, layer_in[i][j][0], layer_in[i][j][1], secTag[i],'-rho', mass_per_unit_length[i])
         k+=1
         
# array to store node ID and z_in value
push = np.array([[i + 1, z_in[i]] for i in range(n_joints)])
# Sort the array by (z_in values)
push = push[push[:, 1].argsort()]
Height = push[-1][1]-push[0][1]
print("Height of the tower is ",Height)

# =============================================================================
#For Hindgeing the coss memebers 
# =============================================================================
from collections import defaultdict

def find_midpoint(element):
    # Get the nodes connected to the element
    nodes = op.eleNodes(element)
    
    # Initialize sums for coordinates
    x_sum, y_sum, z_sum = 0.0, 0.0, 0.0

    # Loop through each node and sum up the coordinates
    for node in nodes:
        coord = op.nodeCoord(node)
        x_sum += coord[0]
        y_sum += coord[1]
        z_sum += coord[2]

    # Calculate the number of nodes
    num_nodes = len(nodes)

    # Calculate the midpoint (centroid)
    midpoint = (x_sum / 2, y_sum / 2, z_sum /2)
    
    return midpoint

def find_element_pairs_with_matching_z():
    element_pairs = []
    elements = op.getEleTags()

    # Loop through each pair of elements
    for i in range(len(elements)):
        element_1 = elements[i]
        midpoint1 = find_midpoint(element_1)
        
        for j in range(i + 1, len(elements)):  # Only compare elements ahead in the list
            element_2 = elements[j]
            midpoint2 = find_midpoint(element_2)
            if midpoint1[0] == midpoint2[0] and midpoint1[1] == midpoint2[1] and midpoint1[2] == midpoint2[2]:
                element_pairs.append((element_1, element_2 , midpoint1))

    return element_pairs

# Find and print pairs of elements
matching_pairs = find_element_pairs_with_matching_z()
print("Pairs of elements with matching z-coordinates at at least 2 nodes but differing x/y coordinates:")
for pair in matching_pairs:
    print(pair)


 # =============================================================================
#For Hindgeing the coss memebers 
# =============================================================================        
         
         
x = [0] * (len(layer_in) + 1)
x[0] = 1

for i in range(len(layer_in)):
    x[i+1] = x[i] + len(layer_in[i])

element_ranges = [list(range(x[i], x[i+1])) for i in range(len(layer_in))]
colors = ["red", "blue", "green", "yellow", "cyan", "magenta", "orange"]

vfo.plot_model(
    elementgroups=[element_ranges, colors[:len(element_ranges)]],
    show_nodes='yes',
    show_nodetags='yes',
    show_eletags='no',
    font_size=15,
    setview='3D',
    line_width=3
)

        
         
IDctrlNode = int(push[-1][0]) ;					# node where disp is read for disp control
IDctrlDOF = 1;					# degree of freedom read for disp control (1 = x displacement)
Dmax = 0.1*Height;		# maximum displacement of pushover: 10% roof drift
Dincr = 0.01;				# displacement increment

#New analysis 
Constrant_Type = "Plain" #Mainly used for rigid diaphragm
numberer_Type = "RCM"
system_type = "FullGeneral"
algorith_type = "ModifiedNewton"
Integrator_type = "Newmark"
N_Gamma = 0.5
N_Beta = 0.25
analysis_type = "Transient"
Tol = 1.0e-3
maxNumIter = 5
testTypeDynamic = "NormUnbalance"
