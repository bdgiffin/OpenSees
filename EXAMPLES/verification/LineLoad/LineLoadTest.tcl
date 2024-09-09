# written: bdg
#
# purpose: ensure that the LineLoad element and load pattern are functioning

#create the ModelBuilder object
model basic -ndm 3 -ndf 3

# build the model

# add nodes - command: node nodeId xCrd yCrd zCrd
node 1   0.0   0.0   0.0
node 2 100.0   0.0   0.0
node 3   0.0 100.0   0.0
node 4   0.0   0.0 100.0

# add material - command: material <matType> matID <matArgs>
uniaxialMaterial Elastic 1 3000

# add truss elements - command: truss trussID node1 node2 A matID
element truss 1 1 2 1.0 1
element truss 2 1 3 1.0 1
element truss 3 1 4 1.0 1
element truss 4 2 3 1.0 1
element truss 5 3 4 1.0 1
element truss 6 4 2 1.0 1

# add LineLoad elements - command: LineLoad LineLoadID node1 node2 radius lib
element LineLoad  7 1 2 0.5 line_load_example
element LineLoad  8 1 3 0.5 line_load_example
element LineLoad  9 1 4 0.5 line_load_example
element LineLoad 10 2 3 0.5 line_load_example
element LineLoad 11 3 4 0.5 line_load_example
element LineLoad 12 2 4 0.5 line_load_example

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt? zRestrnt?
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1

# create the load pattern
pattern Plain 1 Linear {
    eleLoad -range 7 12 -type LineLoader
}

# build the components for the analysis object
system BandSPD
constraints Plain
integrator LoadControl 1.0 
algorithm Linear
numberer RCM

# create the analysis object 
analysis Static 

# create a Recorder object for the nodal displacements at node 4
recorder Node -file Example.out -load -nodes 4 -dof 1 2 3 disp

# perform the analysis
analyze 10

# print the results at node 4 and at all elements
print node 4
print ele

exit
