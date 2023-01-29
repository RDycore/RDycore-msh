
import os
import copy
import numpy as np
import netCDF4 as nc
from scipy import spatial
from skimage.filters import gaussian

from util.loadshp import loadshp
from util.utility import addpoly, addline, innerto, fixgeom, dist_to
import jigsawpy

NAME = "mid_atlantic"

RSPH = +6371220.0  # earth radius [m]

XMID = -75.2316  # projection centre: lon. [deg]
YMID = +39.1269  # projection centre: lat. [deg]

HMAX = +2500.0  # max. mesh spacing [m]
HMIN = +250.0  # min. mesh spacing [m]
DHDX = +0.1  # limit on spacing gradation |dh/dx|
XFLD = +1000.0  # floodplain width [m]

RTAG = +100  # ID flag for river polylines

opts = jigsawpy.jigsaw_jig_t()  # user-opts
proj = jigsawpy.jigsaw_prj_t()  # coord. projection

temp = jigsawpy.jigsaw_msh_t()

# mesh objects on sphere
gsph = jigsawpy.jigsaw_msh_t()  # geometry
isph = jigsawpy.jigsaw_msh_t()  # initial points
ssph = jigsawpy.jigsaw_msh_t()  # spacing h(x)
msph = jigsawpy.jigsaw_msh_t()  # generated mesh

# mesh objects in projection
gprj = jigsawpy.jigsaw_msh_t()
iprj = jigsawpy.jigsaw_msh_t()
sprj = jigsawpy.jigsaw_msh_t()
mprj = jigsawpy.jigsaw_msh_t()

HERE = os.path.abspath(os.path.dirname(__file__))
VDIR = os.path.join(HERE, "vtk")
MDIR = os.path.join(HERE, "msh")


def narvis_tol(feat):               # filer enclosed streams

    return feat["properties"]["UP_CELLS"] > 1000


def build_proj():
    """
    Setup a local stereographic projection to rotate between
    spherical and local coordinate systems.

    """
    proj.prjID = "stereographic"
    proj.radii = RSPH
    proj.xbase = XMID * np.pi / +180.  # lon,lat
    proj.ybase = YMID * np.pi / +180.

    return


def build_geom():
    """
    Extract the river + watershed geometry to be represented

    """
    print("*building geom.")

#------------------------------------ set watershed boundary
    filename = [
    "NHD_H_0204_HU4_Shape/Shape/WBDHU4.shp",
    "NHD_H_0205_HU4_Shape/Shape/WBDHU4.shp",
    "NHD_H_0206_HU4_Shape/Shape/WBDHU4.shp",
    "NHD_H_0207_HU4_Shape/Shape/WBDHU4.shp",
    "NHD_H_0208_HU4_Shape/Shape/WBDHU4.shp"
        ]

#   add watershed geometry as "exterior" polygons, each sets
#   a polygon to be meshed
#   each watershed polygon is assigned an ID, such that
#   corresponding mesh cells can later be identified
#       mesh.tria3["IDtag"] = ID
#   IDs are currently just an index into the filenames above

#   multipolygons handled using jigsaw's msh_t.bound

    gsph.mshID = "ellipsoid-mesh"
    gsph.radii = RSPH * np.ones(3)
    for fpos in range(len(filename)):
        fdir = os.path.join(
            HERE, "data", filename[fpos])
        loadshp(fdir, temp)
        temp.point["coord"] *= np.pi / +180.
        addpoly(gsph, temp, fpos + 1)

#------------------------------------ approx. stream network
    filename = os.path.join(
        HERE, "data", "namerica_rivers", "narivs.shp")

    loadshp(filename, temp, narvis_tol)
    temp.point["coord"] *= np.pi / +180.

#-- keep rivers inside watersheds
    itag = innerto(temp.vert2["coord"], gsph)

    keep = np.logical_and.reduce((
        itag[temp.edge2["index"][:, 0]] > +0,
        itag[temp.edge2["index"][:, 1]] > +0
    ))
    temp.edge2 = temp.edge2[keep]

    fixgeom(temp)  # remove duplicate nodes and reindex

#   add river geometry as "interior" polylines - constraints
#   that mesh edges will follow (stream burning), but not
#   part of any bounding polygons

    addline(gsph, temp, RTAG)

#------------------------------------ rotate to local coord.
    gprj.mshID = copy.deepcopy(gsph.mshID)
    gprj.point = copy.deepcopy(gsph.point)
    gprj.edge2 = copy.deepcopy(gsph.edge2)
    gprj.bound = copy.deepcopy(gsph.bound)
    jigsawpy.project(gprj, proj, "fwd")

    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_geom_prj.vtk"), gprj)
    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_geom_sph.vtk"), gsph)

    return


def build_init():
    """
    Extract any nonmanifold nodes in the river network to be
    used as fixed points in the mesh-gen.

    """
    print("*building init.")

#------------------------------------ get topological degree
    flag = gprj.edge2["IDtag"] == RTAG  # rivers

    erof = gprj.edge2["index"][flag]
    ndeg = np.bincount(np.reshape(erof, erof.size))

#   add any river polyline vertices as initial points if the
#   topo.-degree!=2, e.g. at confluences or endpoints

    take = np.full(gprj.point.size, False, dtype=bool)
    take[ndeg != 2] = True
    take[ndeg == 0] = False

    iprj.mshID = gprj.mshID
    iprj.point = gprj.point[take]

#------------------------------------ rotate to major coord.
    isph.mshID = copy.deepcopy(iprj.mshID)
    isph.point = copy.deepcopy(iprj.point)
    isph.edge2 = copy.deepcopy(iprj.edge2)
    jigsawpy.project(isph, proj, "inv")

    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_init_prj.vtk"), iprj)
    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_init_sph.vtk"), isph)

    return


def build_spac():
    """
    Construct a mesh-spacing function h(x) for the domain.

    """
    print("*building spac.")

#   a DEM dataset, e.g.: 
#       gebco.net/data_and_products/gridded_bathymetry_data/
    data = nc.Dataset(os.path.join(
        HERE, "data", "GEBCO",
    "gebco_2022_n44_s34_w-81_e-72.nc"), "r")

    ssph.mshID = "ellipsoid-grid"
    ssph.radii = RSPH * np.ones(3)
    ssph.xgrid = np.asarray(
        data["lon"][:], dtype=ssph.REALS_t)
    ssph.xgrid*= np.pi / 180.
    ssph.ygrid = np.asarray(
        data["lat"][:], dtype=ssph.REALS_t)
    ssph.ygrid*= np.pi / 180.    

#------------------------------------ dist.-to-rof heuristic
    dist = dist_to(ssph, gsph, RTAG)

    ssph.value = dist[:]
    ssph.value[dist <= XFLD * 1.] = HMIN
    ssph.value = np.minimum(ssph.value, HMAX)
    ssph.value = np.maximum(ssph.value, HMIN)

#------------------------------------ set limits on |dh/dx|.
    ssph.slope = DHDX * np.ones(
        (ssph.ygrid.size, ssph.xgrid.size))

    ssph.slope[dist <= XFLD * 2.]*= .50

    jigsawpy.lib.marche(opts, ssph)

#------------------------------------ rotate to local coord.
    sprj.mshID = copy.deepcopy(ssph.mshID)
    sprj.xgrid = copy.deepcopy(ssph.xgrid)
    sprj.ygrid = copy.deepcopy(ssph.ygrid)
    sprj.value = copy.deepcopy(ssph.value)
    jigsawpy.project(sprj, proj, "fwd")

    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_spac_prj.vtk"), sprj)
    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_spac_sph.vtk"), ssph)

    return


def build_mesh():
    """
    Given the objects defined above, mesh the full domain.

    """
    print("*building mesh.")

#------------------------------------ set jigsaw's file obj.
    opts.geom_file = \
        os.path.join(MDIR, NAME + "_geom_prj.msh")
    opts.init_file = \
        os.path.join(MDIR, NAME + "_init_prj.msh")
    opts.hfun_file = \
        os.path.join(MDIR, NAME + "_spac_prj.msh")

    opts.mesh_file = \
        os.path.join(MDIR, NAME + "_mesh_prj.msh")

    opts.jcfg_file = \
        os.path.join(MDIR, NAME + "_opts.jig")

    jigsawpy.savemsh(opts.geom_file, gprj)
    jigsawpy.savemsh(opts.init_file, iprj)
    jigsawpy.savemsh(opts.hfun_file, sprj)

#------------------------------------ set jigsaw's useropts.
#   e.g.: github.com/dengwirda/jigsaw/wiki/jig-file-format
    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")  # null spacing lim
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2  # topologically 2-dim.
    opts.mesh_rad2 = +1.20  # relax radius-edge constraint

   #opts.optm_iter = +0  # to disable mesh optim.

   #opts.verbosity = +1  # to show more output

#------------------------------------ mesh-gen. using jigsaw
    jigsawpy.cmd.jigsaw(opts, mprj)

#------------------------------------ rotate to major coord.
    msph.mshID = copy.deepcopy(mprj.mshID)
    msph.point = copy.deepcopy(mprj.point[:])
    msph.edge2 = copy.deepcopy(mprj.edge2[:])
    msph.tria3 = copy.deepcopy(mprj.tria3[:])
    jigsawpy.project(msph, proj, "inv")

    jigsawpy.savevtk(os.path.join(
        MDIR, NAME + "_mesh_sph.msh"), msph)

    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_mesh_prj.vtk"), mprj)
    jigsawpy.savevtk(os.path.join(
        VDIR, NAME + "_mesh_sph.vtk"), msph)

    return


if (__name__ == "__main__"):
#----------------- build input objects and mesh using jigsaw
    build_proj()
    build_geom()
    build_init()
    build_spac()
    build_mesh()



