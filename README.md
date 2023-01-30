## `Meshing for RDycore`

A set of python scripts to generate unstructured meshes for the `RDycore` project.

### `Required Software`

To run the meshing workflows you will need:

- A working python >=3.6 installation.
- A number of python packages (see requirements.txt).
- Paraview (or equiv.), for visualisation.

### `Defining a Configuration`

New meshes/domains can be added by defining a new top-level script. The `mid-atlantic.py` config. can be used as an example, defining the following routines to build inputs for the [`jigsaw` meshing library](https://github.com/dengwirda/jigsaw):

- build_proj: make a stereographic projection to map local <=> spherical coord. systems.
- build_geom: load watershed and river geometry datasets to define the domain geometry.
- build_init: define "fixed-points" to be included as mesh initial conditions.
- build_spac: define a nonuniform mesh-spacing function h(x) to control mesh resolution.
- build_mesh: call jigsaw to build the corresponding unstructured mesh.

The `mid_atlantic.py` example meshes a number of watersheds adjacent to the US mid-Atlantic coastline, using river geometry derived from the [HydroSHEDS](https://www.hydrosheds.org/) product, and watershed geometries provided by the [USGS National Hydrography Dataset](https://www.usgs.gov/national-hydrography/national-hydrography-dataset). Stream geometry is "burned-into" the mesh through the alignment of cell edges and mesh resolution is adapted to resolve stream-adjacent regions. Output mesh files are stored in the `msh` and `vtk` directories.

### `Required Datasets`

Underlying geographical datasets (watersheds, stream networks, DEMs, etc) are attached to the current `release` of this repository and must be downloaded before workflows can be run. The `get_datasets.py` script can be used to download and unpack all datasets into an `RDycore-msh` installation.

### `Example Installation`

A conda environment can be used to wrap dependencies for the `RDycore-msh` scripts, allowing the workflows to be run on laptops and HPC machines alike:
````
conda create --name rdycore-msh python=3.10
conda activate rdycore-msh
````
Install python dependencies into the `rdycore-msh env`:
````
pip install netCDF4 fiona geojson
pip install numpy scipy scikit-image
pip install inpoly
pip install githubrelease
````
The `jigsaw` library can be installed into the `rdycore-msh env` from src., ensuring the latest version is linked:
````
git clone https://github.com/dengwirda/jigsaw-python.git
cd jigsaw-python
python setup.py build_external install
cd ..
conda deactivate
````
The `rdycore-msh env` is now built. If you need to delete it, use `conda env remove --name rdycore-msh`.

The `rdycore-msh env` should be activated whenever you run `RDycore-msh` workflows, for example: 
````
git clone https://github.com/rdycore/rdycore-msh.git
conda activate rdycore-msh
python rdycore-msh/get_datasets.py
python rdycore-msh/mid_atlantic.py
````

### `TODO`

Development tasks (in no particular order):

- Are we using the best datasets? Probably not: investigate VoTe, etc.
- River geometry is pixel-aligned: smooth raw polylines via splines, etc?
- Optimal h(x) heuristics to resolve rivers, floodplains, steep topography, etc?
- Quad-based (anisotropic) meshing for river channels.
- Improved optimisation of cells/vertices attached to interior constraints.



