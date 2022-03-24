# FEBioVTK
Plugin for FEBio that creates a VTK plot file e.g. for visualization in ParaView.

Requires that VTK version 8 is installed.

# TODO
* [ ] Select which fields to output (Currently list of names is hardcoded. It is difficult to access list of names from febio-plotfile in a plugin)
* [ ] When to plot as specified in control section. (Currently plots when simulation has been solved)