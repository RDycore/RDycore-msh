## `Meshing for RDycore Solvers`

A set of python scripts to generate unstructured meshes for the `RDycore` project.

### `Required Software`

To run the meshing workflows you will need:

    * A working `python >=3.6` installation.
    * The `numpy`, `scipy`, `scikit-image`, `inpoly` and `jigsawpy` packages (see `requirements.txt`).

### `Defining a Configuration`

New meshes/domains can be added by defining a new top-level script. The `mid-atlantic.py` config. can be used as an example, defining the following routines to build inputs for the `jigsaw` meshing library:

    * build_proj: setup a stereographic projection to map between local and spherical coordinate systems.
    * build_geom: load a set of watershed and river geometry datasets to define the domain geometry.
    * build_init: extract any "fixed-points" from the geometry to be included as mesh initial conditions.
    * build_spac: define a nonuniform mesh-spacing function h(x) to control mesh resolution.
    * build_mesh: call the jigsaw library to build the corresponding unstructured mesh.

The `mid_atlantic.py` example meshes a number of watersheds adjacent to the US mid-Atlantic coastline, using river geometry derived from the HydroSHEDS product, and watershed geometries provided by the USGS National Hydrography Dataset. Stream geometry is "burned-into" the mesh through the alignment of cell edges and mesh resolution is adapted to resolve stream-adjacent regions.


### `Required Datasets`

Underlying geographical datasets (watersheds, stream networks, DEMs, etc) are attached to the current `release` of this repository and must be downloaded before workflows can be run.

### `TODO`

Development tasks (in no particular order):

    * Are we using the best datasets? Probably not: investigate VoTe, etc.
    * River geometry is typically pixel-aligned: is it useful to smooth raw polylines via splines, etc?
    * What is an optimal set of resolution heuristics to resolve rivers, floodplains, steep topography, etc?
    * Quad-based (anisotropic) meshing for river channels.
    * Improved optimisation of cells/vertices attached to interior constraints.



