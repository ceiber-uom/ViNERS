
# \source\
Optional folder that contains raw data prior to conversion to the format contained in the primary data folder.
This folder contains the following

## \source\fascicles\
splinedata from microscopy sections or generated from the literature suitable for import into GMSH using mesh.insert_gmsh_fascicles
- Bertrand.splines.dat
- Havton.splines.dat
- Hulesbosch.splines.dat

## \source\axons\
axon population data from EM or from the literature, needed for +models.axon_population

## \source\mesh\
default generic GMSH .geo model framework, needed for +models.nerve_model.These files get translated to a gmsh .geo file using tools.makeFromTemplate(). geometric discriptions of the model features are implemented in +mesh
Files:
- pelvic_nerve.geo.template
- pelvic_nerve.geo.variable

