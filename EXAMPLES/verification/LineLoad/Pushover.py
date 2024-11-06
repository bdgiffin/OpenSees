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
       
#=============================================================================
# define MATERIAL -------------------------------------------------------------
# =============================================================================
Fy = 345 * Mpa
Es = 200* Gpa  # Steel Young's Modulus
nu = 0.3
Gs = Es / (2 * (1 + nu))  # Torsional stiffness Modulus
J = 10  # Large torsional stiffness
Bs = 0.01
R0 = 20
cR1 = 0.925
cR2 = 0.15
a1 = 0.39
a2 = 1.0
a3 = 0.029
a4 = 1.0
matIDhard = 2
matType = 'Steel02'
# Function to define uniaxial material in Python
op.uniaxialMaterial(matType, matIDhard, Fy, Es, Bs, R0, cR1, cR2)

# =============================================================================
#Transformation based on the orientation
# =============================================================================
Trans_Type = "Corotational" #'Linear 
Transf = [1,2,3,4,5,6,7]

op.geomTransf(Trans_Type, Transf[0], 0, 0, 1)     #main vertical members
op.geomTransf(Trans_Type, Transf[1], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[2], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[3], 0, 1, 1)  
op.geomTransf(Trans_Type, Transf[4], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[5], 0, 0, 1)  #Cross Members
op.geomTransf(Trans_Type, Transf[6], 0, 0, 1)  #Cross Members

# =============================================================================
# Define Section Sizes
# =============================================================================
secTag = [1,2,3,4,5,6,7]
BreID = 2
Lfiber = 20
Sfiber = 3
Thick = np.array([12*mm,12*mm,10*mm,6*mm,5*mm,5*mm,5*mm]) #should be based on number of layers
Length = np.array([130*mm,110*mm,100*mm,60*mm,45*mm,50*mm,50*mm]) #should be based on number of layers
Ly1= -Thick/2
Hy1= -Thick/2
Ly2= Length-Thick/2
Hy2= Thick/2    

# =============================================================================
# Define Weights of the members
# =============================================================================

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
# for i in range(len(layer_in)):
#     op.section('Elastic', secTag[i], Es, cross_sectional_area[i], moment_of_area_x[i], moment_of_area_y[i], shear_modulus, polar_moment_of_area[i])

# =============================================================================
# Defining Hysteric Model 
# =============================================================================
# Paper for hysteric Modeling https://www.tandfonline.com/doi/pdf/10.1080/15732479.2019.1673783
#From paper https://www.researchgate.net/profile/Davoud-Nezamolmolki/publication/349346125_The_Effect_of_Nonlinear_Behavior_of_Bolted_Connections_on_Dynamic_Analysis_of_Steel_Transmission_Towers/links/60b8672aa6fdcc22eacf4aca/The-Effect-of-Nonlinear-Behavior-of-Bolted-Connections-on-Dynamic-Analysis-of-Steel-Transmission-Towers.pdf


# =============================================================================
# Plot Hysteric Graph for joint slippage just for refence not used in opensees model 
# =============================================================================
# x =np.array([0*mm, 46.95/139.0*mm, 0.85*mm,1.85*mm,1.16*mm]) #Single Leg bolted connection B
# y = np.array([0*KN,46.95*KN,46.95*KN,168.2*KN,201.6*KN]) #Single Leg bolted connection B

# # Initialize x1 with the same length as x and accumulate the values
# x1 = np.zeros(len(x))
# x1[0] = x[0]  # First element
# for i in range(1, len(x)):
#     x1[i] = x[i] + x1[i-1]  # Accumulate the values

# # Combined Hyseric Plot
# plt.plot(x1, y, linewidth=3, label='Main Plot')

# # Additional plots for 3 hysteric Model 
# plt.plot(x1[0:4], [y[0], y[1], y[2], 0], label='Plot 1', linewidth=2) #1st hysteric Model
# plt.plot(x1[0:5], [0, 0, 0, y[3], 0], label='Plot 2', linewidth=2) #2nd hysteric Model
# plt.plot(x1[0:5], [0, 0, 0, 0, y[4]], label='Plot 3', linewidth=2) #3rd hysteric Model

# # Plot dotted vertical lines from each y value to the x-axis
# for i in range(len(x1)):
#     plt.vlines(x=x1[i], ymin=0, ymax=y[i], colors='gray', linestyles='dotted')

# # Add labels to each section
# for i in range(len(x1)-1):
#     plt.text((x1[i]+x1[i+1])/2, -10 , f' Phase{i+1}', ha='center', fontsize=9)


# # Labels and title
# plt.xlabel('Displacement (mm)', fontsize=14)
# plt.ylabel('Force (KN)', fontsize=14)
# plt.title('Plot of Hysteric Models',fontsize=20)

# # Add a legend
# plt.legend(["Hysteric Model combined","1st Hysteric Model","2nd Hysteric Model","3rd Hysteric Model"],fontsize = 12 )

# # Show the plot
# plt.show()

# =============================================================================
# Defining Hysteric Model used in opensees model
# =============================================================================
# hyst_Tag = [5,6,7]
# parallel = 8
# p1 = [x1[1],y[1]]
# p2 = [x1[2],y[2]]
# p3 = [x1[3],0]
# print(  y[1],x1[1],
#         y[2],x1[2],
#         0,x1[3],
#         -y[1],-x1[1],
#         -y[2],-x1[2],
#         0,-x1[3],
#                     )
# def hysteric(y1,x1,y2,x2,y3,x3,tag):
#     op.uniaxialMaterial('Hysteretic',tag,
#                     y1,x1,
#                     y2,x2,
#                     y3,x3,
#                     -y1,-x1,
#                     -y2,-x2,
#                     -y3,-x3,
#                     0,0,
#                     0,0
#                     )
    
# hysteric(y[1],x1[1],y[2],x1[2],0,x1[3],hyst_Tag[0])
# hysteric(0,x1[2],y[3],x1[3],0,x1[4],hyst_Tag[1])
# #hysteric( 0,x1[3],y[4],x1[4],alpha*y[1],x1[4],hyst_Tag[2])
# op.uniaxialMaterial('Parallel',parallel,hyst_Tag[0],hyst_Tag[1],hyst_Tag[2])
      
      
      
# hinge_nodes = list(range(135, 159))
# hinge_nodes.extend([54,73,87,93])#64,76,90,95,130,78,129,96])
# # print(hinge_nodes)
# for hinge in hinge_nodes:
#     op.fix(hinge, 1, 1, 1, 0, 0, 0) 
# =============================================================================
# Defining Nodal Mass 
# =============================================================================
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
            i1 = layer_in[i][j][0]-1
            i2 = layer_in[i][j][1]-1
            dx = x_in[i2] - x_in[i1]
            dy = y_in[i2] - y_in[i1]
            dz = z_in[i2] - z_in[i1]
            member_length = math.sqrt(dx*dx + dy*dy + dz*dz)
            member_mass   = steel_mass_density*cross_sectional_area[i]*member_length
            lumped_mass[i1] = lumped_mass[i1] + 0.5*member_mass
            lumped_mass[i2] = lumped_mass[i2] + 0.5*member_mass

# nodal masses: (only needed if masses are not already computed by the element)
op.timeSeries("Linear", 1)
op.pattern("Plain", 1, 1)
for i in range(0,n_joints):
     op.mass(i+1, lumped_mass[i], lumped_mass[i], 0.0) # [kg] node#, Mx My Mz, Mass=Weight/g.
     op.load(i+1, 0.0, 0.0, -lumped_mass[i]*g)
# =============================================================================
# Defining Fiber Section
# =============================================================================
from Modules.Fiber import *  
for i in range(len(layer_in)):
    FiberCreation(secTag[i],matIDhard,Sfiber,Lfiber,Ly1[i],Hy1[i],Ly2[i],Hy2[i])
# =============================================================================
# Defining Element Section
# =============================================================================
k = 0
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
         op.element('nonlinearBeamColumn', k+1, layer_in[i][j][0], layer_in[i][j][1], 10, secTag[i], Transf[i])
         k+=1
# =============================================================================
# RECORDER -------------------------------------------------------------
# =============================================================================
# Ensure dataDir is a valid directory path
dataDir = "output"
if not os.path.exists(dataDir):
    os.makedirs(dataDir)
    print(f"Directory '{dataDir}' created.")
    
# array to store node ID and z_in value used for finding the topmost and botom nodes
push = np.array([[i + 1, z_in[i]] for i in range(n_joints)])
# Sort the array by (z_in values)
push = push[push[:, 1].argsort()]
Height = push[-1][1]-push[0][1]
print("Height of the tower is ",Height)

# =============================================================================
# For 3-d Visualization 
# =============================================================================
# x = [0] * (len(layer_in) + 1)
# x[0] = 1
# for i in range(len(layer_in)):
#     x[i+1] = x[i] + len(layer_in[i])

# element_ranges = [list(range(x[i], x[i+1])) for i in range(len(layer_in))]
# colors = ["red", "blue", "green", "yellow", "cyan", "magenta", "orange"]

# vfo.plot_model(
#     elementgroups=[element_ranges, colors[:len(element_ranges)]],
#     show_nodes='yes',
#     show_nodetags='yes',
#     show_eletags='no',
#     font_size=15,
#     setview='3D',
#     line_width=3
# )
# exit()
# =============================================================================
# Recorder for max displacement at the top and base reaction
# =============================================================================
free_file = os.path.join(dataDir, "DFree.out")
fixed_file = os.path.join(dataDir, "DFixed.out")
react = os.path.join(dataDir, "RXN.out")
# Define a recorder for the drift
op.recorder("Node", '-file', free_file, 'time', '-node', int(push[-1][0]),'-precision',3, '-time' ,'-dof', 1, 'disp')
op.recorder("Node", '-file', fixed_file, 'time', '-node', int(push[0][0]),'-precision',3, '-time' ,'-dof', 1, 'disp')
op.recorder("Node", '-file', react, 'time', '-node', int(push[0][0]), '-precision',3,'-time', '-dof', 1, 'reaction')
# =============================================================================
# Applying Lateral load pattern  ----------------------------------------------
# =============================================================================
#file to store pushover data
height_ranges = [2, 5, 15, 25, 36]
load_push = [10,20,30,40,50]

unique_z_in = np.unique(push[:, 1])
op.timeSeries('Linear',2)
op.pattern('Plain',2, 2)
for idx, z_val in enumerate(unique_z_in):
    nodes_at_same_z = push[push[:, 1] == z_val][:, 0]  # Extract node IDs for the current z_in value
    
    for i in range(len(height_ranges)):
        if z_val<=height_ranges[i]:
            check = i
            break 
    load_per_node = load_push[check]/len(nodes_at_same_z)  
    # print(z_val,check,load_push[check],len(nodes_at_same_z),load_per_node)
    # print(nodes_at_same_z)
    for node in nodes_at_same_z:
        op.load(int(node), load_per_node, 0, 0 , 0, 0, 0)  # Apply load in the z direction

    # Optionally increase the base load for the next iteration
      # Increase the base load for higher z_in values


IDctrlNode = int(push[-1][0]) ;# node where disp is read for disp control
IDctrlDOF = 1;# degree of freedom read for disp control (1 = x displacement)
Dmax = 0.8;	# maximum displacement of pushover:
Dincr = 0.01;# displacement increment

# =============================================================================
# Pushover Analysis
# =============================================================================
Constrant_Type = "Plain" 
numberer_Type = "RCM"
system_type = "BandGeneral"
algorith_type = "ModifiedNewton"
Integrator_type = "Newmark"
N_Gamma = 0.5
N_Beta = 0.25
analysis_type = "Transient"
Tol = 1.0e-3
maxNumIter = 5
testTypeDynamic = "NormUnbalance"

op.constraints(Constrant_Type)    # how it handles boundary conditions
op.numberer(numberer_Type)        # renumber dof's to minimize band-width (optimization), if you want to
op.system(system_type)            # how to store and solve the system of equations in the analysis
op.test(testTypeDynamic, Tol, 500) #determine if convergence has been achieved at the end of an iteration step
op.algorithm(algorith_type)	   # use Linear algorithm for linear analysis
op.integrator("DisplacementControl", IDctrlNode, IDctrlDOF, Dincr)
op.analysis('Static')
Nsteps = int(Dmax/Dincr)
op.analyze(Nsteps)
print('Nsteps')
print("Pushover Complete")
# =============================================================================
# Plot the Pushover Curve -----------------------------------------------------
# =============================================================================

time_series1 = np.loadtxt(free_file)
time_series2 = np.loadtxt(fixed_file)
Base = np.loadtxt(react)
drift = (time_series1[:, 1]) # 3.0 is the height between node 1 and node 2
BaseRxn = abs(Base[:,1])
plt.plot(drift, BaseRxn/1000, color='blue', linewidth=2, label='Base Rxn vs Displacement')
plt.ylabel("Base Reaction (kN)", fontsize=14)
plt.xlabel("Displacement at the Top (m)", fontsize=14)
plt.title("Pushover Curve", fontsize=16)
plt.legend(loc='best', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)  
plt.tight_layout()
plt.show()