
import copy
import numpy as np
from scipy import spatial
from scipy import interpolate
import jigsawpy

from inpoly import inpoly2


def addpoly(geom, poly, itag):
    """
    ADDPOLY: add new closed polygon POLY to mst_t obj. GEOM.

    The POLY.POINT + POLY.EDGE2 arrays are appended to GEOM,
    and a new "loop" added to GEOM.BOUND. This new loop is
    assigned ID = ITAG.

    """
    temp = jigsawpy.jigsaw_msh_t()

    temp.point = poly.point
    temp.edge2 = poly.edge2

    zipmesh(temp)                       # ensure compressed

    temp.bound = np.empty(
        (poly.edge2.size), dtype=geom.BOUND_t)
    temp.bound["index"] = \
        np.arange(0, poly.edge2.size)
    temp.bound["cells"] = \
        jigsawpy.jigsaw_def_t.JIGSAW_EDGE2_TAG
    temp.bound["IDtag"] = itag

    temp.edge2["index"] += geom.vert2.size
    temp.bound["index"] += geom.edge2.size

    geom.point = np.concatenate(
        (geom.vert2, temp.vert2), axis=0)
    geom.edge2 = np.concatenate(
        (geom.edge2, temp.edge2), axis=0)
    geom.bound = np.concatenate(
        (geom.bound, temp.bound), axis=0)

    return


def addline(geom, line, itag):
    """
    ADDLINE: push new open polyline LINE to mst_t obj. GEOM.

    The LINE.POINT + LINE.EDGE2 arrays are appended to GEOM.
    The new edges are assigned ID = ITAG.

    """
    temp = jigsawpy.jigsaw_msh_t()

    temp.point = line.point
    temp.edge2 = line.edge2

    temp.edge2["IDtag"] = itag

    zipmesh(temp)                       # ensure compressed

    temp.edge2["index"] += geom.vert2.size

    geom.point = np.concatenate(
        (geom.vert2, temp.vert2), axis=0)
    geom.edge2 = np.concatenate(
        (geom.edge2, temp.edge2), axis=0)

    return


def zipmesh(mesh):
    """
    ZIPMESH: "zip" a mst_t obj., pruning any unused points /
    cells, and compressing cell indexing.

    """
    used = np.full(
        mesh.point.size, False, dtype=bool)

#---------------------------------- flag nodes used in cells
    if (mesh.edge2 is not None and
            mesh.edge2.size > +0):

        used[mesh.edge2[
            "index"].reshape(-1)] = True

    if (mesh.tria3 is not None and
            mesh.tria3.size > +0):

        used[mesh.tria3[
            "index"].reshape(-1)] = True

    if (mesh.quad4 is not None and
            mesh.quad4.size > +0):

        used[mesh.quad4[
            "index"].reshape(-1)] = True

#---------------------------------- re-index cells, compress
    redo = np.full(
        mesh.point.size, 0, dtype=np.int32)

    redo[used] = np.arange(
        0, np.count_nonzero(used))

    if (mesh.edge2 is not None and
            mesh.edge2.size > +0):

        mesh.edge2["index"] = \
            redo[mesh.edge2["index"]]

    if (mesh.tria3 is not None and
            mesh.tria3.size > +0):

        mesh.tria3["index"] = \
            redo[mesh.tria3["index"]]

    if (mesh.quad4 is not None and
            mesh.quad4.size > +0):

        mesh.quad4["index"] = \
            redo[mesh.quad4["index"]]

#---------------------------------- prune any un-used points
    mesh.point = mesh.point[used]

    return


def fixgeom(geom):
    """
    FIXGEOM: remove any duplicate vertices and reindex cells

    """
    xpos = np.round(geom.point["coord"], decimals=9)

    __, ifwd, jinv = np.unique(
        xpos, axis=0, 
        return_index=True, return_inverse=True)

    geom.point = geom.point[ifwd]
    geom.edge2["index"] = jinv[geom.edge2["index"]]


    indx = geom.edge2["index"]

    __, ifwd, jinv = np.unique(
        indx, axis=0, 
        return_index=True, return_inverse=True)

    geom.edge2 = geom.edge2[ifwd]

    return


def innerto(vert, geom):
    """
    INNERTO: return ID-tags of polygons enclosing each point
    in VERT. Return -1 for exterior points.

    """
    itag = np.full(
        vert.shape[+0], -1, dtype=np.int32)

    imin = np.min(geom.bound["IDtag"])
    imax = np.max(geom.bound["IDtag"])

    for ipos in range(imin, imax + 1):

#---------------------------------- inpolygon for POLY=ITAG
        indx = geom.bound["index"][
            geom.bound["IDtag"] == ipos]

        xpts = geom.vert2["coord"]
        edge = geom.edge2["index"][indx]

        mask, _ = inpoly2(vert, xpts, edge)

        itag[mask] = ipos

    return itag


def sphdist(rsph, xone, yone, xtwo, ytwo):
    """
    SPHDIST: return the distance from the points [XONE,YONE]
    to the point [XONE,YONE] on a sphere of radii RSPH.

    """
    dlon = .5 * (xone - xtwo)
    dlat = .5 * (yone - ytwo)

    dist = 2. * rsph * np.arcsin(np.sqrt(
        np.sin(dlat) ** 2 +
        np.sin(dlon) ** 2 * np.cos(yone) * np.cos(ytwo)
    ))

    return dist


def rd_dist(xone, xtwo):
    """
    RD_DIST: return the distance between the points XONE and
    XTWO in R^D. Points are N-by-D arrays.

    """
    return np.sqrt(np.sum((xtwo - xone) ** 2, axis=1,
                          keepdims=True))


def dist_to(grid, poly, ptag=None):
    """
    DIST-TO: return the (geodesic) dist. from the mst_t obj.
    GRID to the polyline obj. POLY.

    Set PTAG to filter edges in POLY via IDs.

    """
    opts = jigsawpy.jigsaw_jig_t()
    near = jigsawpy.jigsaw_msh_t()

    dist = copy.deepcopy(grid)
    geom = copy.deepcopy(poly)

    if (grid.mshID.lower() == "ellipsoid-grid"):
#------------------------------------ find dist. to geo: S^2
        xdel = min(np.diff(grid.xgrid))
        ydel = min(np.diff(grid.ygrid))
        dlen = min(xdel, ydel) * np.mean(grid.radii)

        divgeom(geom, dlen)

        if (ptag is None):
            mark = np.full(
                geom.point.size, True, dtype=bool)
        else:
            mark = geom.edge2["IDtag"] == ptag
            mark = np.unique(
                geom.edge2["index"][mark])

        tree = spatial.cKDTree(
            jigsawpy.S2toR3(
            geom.radii, geom.point["coord"][mark]))

        xgrd, ygrd = \
            np.meshgrid(dist.xgrid, dist.ygrid)

        xpos = np.vstack(
            (xgrd.flatten(), ygrd.flatten())).T

        xpos = jigsawpy.S2toR3(grid.radii, xpos)

        dval, _ = tree.query(
            xpos, distance_upper_bound=dlen * 2.0)

        dist.slope = np.ones((
            dist.ygrid.size, dist.xgrid.size))

        dist.value = \
            np.reshape(dval, dist.slope.shape)

        dist.value = np.minimum(
            dist.value, np.mean(grid.radii) * 8.0)

    else:
        raise ValueError("Unsupported grid type.")

#-- solve |grad(d)| = 1., via jigsaw's fast-marching solver.
#-- returns dist.value as (a PDE-based approximation) to the 
#-- geodesic distance

    jigsawpy.lib.marche(opts, dist)

    return dist.value


def divgeom(geom, spac):
    """
    DIVGEOM: subdivide the edges in the msh_t object GEOM to
    satisfy the spacing threshold SPAC.

    The SPAC. param. may either be a scalar length or a full
    mesh spacing definition (a msh_t obj).

    The GEOM. object is modified in-place.

    """
    kind = geom.mshID.lower()

    if (not isinstance(spac, jigsawpy.jigsaw_msh_t)):
        spac = float(spac)
    else:
        if ("-grid" not in spac.mshID.lower()):
            raise Exception(
                "Unsupported SPAC.MSHID type")

        hfun = interpolate.RectBivariateSpline(
            spac.ygrid, spac.xgrid, spac.value)
    
    while True:                             # while too long

        vert = geom.point["coord"]
        cell = geom.edge2["index"]

    #-------------------------------- eval. spacing at nodes
        if (isinstance(spac, float)):
            hval = np.full(
                (vert.shape[0]), spac, dtype=float)
        else:
            hval = hfun(
                vert[:, 1], vert[:, 0], grid=False)

    #-------------------------------- eval. edge length, x^d
        if (kind == "ellipsoid-mesh"):

            xpos = jigsawpy.S2toR3(
                geom.radii, geom.point["coord"])

            elen = sphdist(
                np.mean(geom.radii),
                vert[cell[:, 0], 0],
                vert[cell[:, 0], 1],
                vert[cell[:, 1], 0],
                vert[cell[:, 1], 1])

        if (kind == "euclidean-mesh"):

            xpos = np.array(vert[:], copy=True)

            elen = rd_dist(
                vert[cell[:, 0], :],
                vert[cell[:, 1], :])

        hmid = 0.5 * (
            hval[cell[:, 0]] + hval[cell[:, 1]]
        )

        mask = elen >= hmid * 4. / 3.       # TRUE to subdiv

        ndiv = np.count_nonzero(mask)

        if (not np.any(mask)): break

        xmid = 0.5 * (
            xpos[cell[:, 0]] + xpos[cell[:, 1]]
        )

    #-------------------------------- subdiv. edge at middle
        if (kind == "ellipsoid-mesh"):

            newx = np.zeros(
                (ndiv), dtype=geom.point.dtype)

            newx["coord"] = R3toS2(
                geom.radii, xmid[mask])
            
        if (kind == "euclidean-mesh"):

            newx = np.zeros(
                (ndiv), dtype=geom.point.dtype)

            newx["coord"] = xmid[mask]

    #-------------------------------- re-index for new edges
        inew = np.arange(
            +0, ndiv) + geom.point.size

        new1 = np.empty(
            (ndiv), dtype=geom.EDGE2_t)
        new1["index"][:, 0] = \
            geom.edge2["index"][mask, 0]
        new1["index"][:, 1] = inew
        new1["IDtag"] = \
            geom.edge2["IDtag"][mask]
        
        new2 = np.empty(
            (ndiv), dtype=geom.EDGE2_t)
        new2["index"][:, 0] = inew
        new2["index"][:, 1] = \
            geom.edge2["index"][mask, 1]
        new2["IDtag"] = \
            geom.edge2["IDtag"][mask]

        geom.edge2[mask] = new1

    #-------------------------------- add new cells to GEOM.
        geom.point = np.concatenate(
            (geom.point, newx), axis=0)

        geom.edge2 = np.concatenate(
            (geom.edge2, new2), axis=0)

    return


def R3toS2(radii, E3):
    """
    R3TOS2: return the LON-LAT coord's associated with a set
    of points in R^3. A (geocentric) projection is first
    done to ensure the points lie on the ellipsoidal surface,
    with the projected points then transformed to [LON,LAT]
    pairs. 

    """   
    PP = .5 * E3

    ax = PP[:, 0] ** 1 / radii[0] ** 1
    ay = PP[:, 1] ** 1 / radii[1] ** 1
    az = PP[:, 2] ** 1 / radii[2] ** 1

    aa = ax ** 2 + ay ** 2 + az ** 2

    bx = PP[:, 0] ** 2 / radii[0] ** 2
    by = PP[:, 1] ** 2 / radii[1] ** 2
    bz = PP[:, 2] ** 2 / radii[2] ** 2

    bb = bx * 2. + by * 2. + bz * 2.

    cx = PP[:, 0] ** 1 / radii[0] ** 1
    cy = PP[:, 1] ** 1 / radii[1] ** 1
    cz = PP[:, 2] ** 1 / radii[2] ** 1

    cc = cx ** 2 + cy ** 2 + cz ** 2
    cc = cc - 1.0

    ts = bb * bb - 4. * aa * cc

    ok = ts >= .0

    AA = aa[ok]; BB = bb[ok]; CC = cc[ok]; TS = ts[ok]

    t1 = (-BB + np.sqrt(TS)) / AA / 2.0
    t2 = (-BB - np.sqrt(TS)) / AA / 2.0

    tt = np.maximum(t1, t2)
    
    P3 = np.zeros(E3.shape, dtype=float)
    P3[ok, 0] = (1. + tt) * PP[ok, 0]
    P3[ok, 1] = (1. + tt) * PP[ok, 1]
    P3[ok, 2] = (1. + tt) * PP[ok, 2]

    return jigsawpy.R3toS2(radii, P3)



