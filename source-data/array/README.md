
# \~/source/array README

This folder contains .template and .variable files for generating .geo files compatabile with GMSH (https://gmsh.info) for different electrode array designs. This folder also contains some example files .geo files (.geo.example) of expected outputs, which can be converted to meshes by GMSH. 

ViNERS uses its inbuilt templating engine (tools.makeFromTemplate) to parse these template files into .geo files (which usually get thrown away after mesh generation)

For help developing your own .geo files, please check out https://gmsh.info/doc/texinfo/gmsh.html#Volumes for geometry primatives, as well as `mesh.insert_gmsh_electrodes` and `mesh.insert_gmsh_fascicles`. There is room for improvement with both of these functions in terms of interacting more robustly with experimenter-specified electrode array geometry

You are not obligated to use mesh.insert_gmsh_electrodes (and indeed this function is unlikely to work with your novel new electrode array design) but I strongly recommend using mesh.insert_gmsh_fascicles. 