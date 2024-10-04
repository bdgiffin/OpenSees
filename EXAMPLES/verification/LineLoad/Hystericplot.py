import openseespy.opensees as op
import matplotlib.pyplot as plt
import math
# Units and constants (unchanged)
# ...
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

# Initialize model
op.wipe()
op.model('basic', '-ndm', 1, '-ndf', 1)
op.node(1, 0.0)
op.node(2, 0.0)

# Fix node 1 and leave node 2 free
op.fix(1, 1)

# Define Steel02 material
Fy = 60.0 * ksi
Es = 29000 * ksi  # Steel Young's Modulus
nu = 0.3
Bs = 0.01
R0 = 18
cR1 = 0.925
cR2 = 0.15
matIDhard = 1
matType = 'Steel02'
op.uniaxialMaterial(matType, matIDhard, Fy, Es, Bs, R0, cR1, cR2)

# Define zeroLength element
op.element('zeroLength', 1, 1, 2, '-mat', matIDhard, '-dir', 1)

# Define time series for cyclic loading
op.timeSeries('Linear', 1)
# Associate time series with a Plain pattern for displacement control
op.pattern('Plain', 1)
op.load(2,20)

# Define analysis parameters
op.system('BandGeneral')
op.numberer('Plain')
op.constraints('Plain')
op.integrator('LoadControl', 0.1)  # Changed to LoadControl
op.algorithm('Newton')
op.test('NormDispIncr', 1.0e-8, 10)
op.analysis('Static')

# Perform analysis and collect stress-strain data
strain_list = []
stress_list = []

num_steps = 50  # Increased number of steps for smoother curve
for step in range(num_steps):
    # Analyze the next step
    op.analyze(1)

    # Get strain (displacement at node 2)
    strain = op.nodeDisp(2, 1)

    # Get stress (material stress in zeroLength element)
    stress = op.eleResponse(1, 'material', 'stress')

    strain_list.append(strain)
    stress_list.append(stress)

# Plotting the hysteresis curve using matplotlib
# plt.figure(figsize=(10, 6))
# plt.plot(strain_list, stress_list, '-', label='Steel02 Hysteresis')
# plt.xlabel('Strain')
# plt.ylabel('Stress (Pa)')
# plt.title('Steel02 Hysteresis Curve')
# plt.grid(True)
# plt.legend()
# plt.show()


print(stress_list[:5])