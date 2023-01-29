
import numpy as np
import fiona
import argparse

from jigsawpy import jigsaw_msh_t, savemsh


def nlevels(lists):
#---------------------------------- find list nesting lvl.
    if isinstance(lists, list):
        if(len(lists) == 0):
            nlev = +1

        else:
            nlev = +1 + max([
                nlevels(l) for l in lists])
    else:
        nlev = +0

    return nlev


def coreshp(data, nset, eset, nobj, last):
    """
    CORESHP: load a shape file into a jigsaw msh_t object.

    """

    if (nlevels(data) > +1):

#---------------------------------- iterate into list lvl.
        for next in data:
            nobj, last = coreshp(
                next, nset, eset, nobj, last)

    elif (nlevels(data) == +1):

#---------------------------------- read last level coord.
        npts = len(data)

        temp = jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (npts + 0), dtype=temp.VERT2_t)
        temp.edge2 = np.zeros(
            (npts - 1), dtype=temp.EDGE2_t)

        temp.edge2["IDtag"][:] = nobj

        nobj = nobj + 1

        indx = np.arange(0, npts - 1) + last

        last = last + npts

        temp.vert2["coord"][:, 0] = \
            [ipos[0] for ipos in data]

        temp.vert2["coord"][:, 1] = \
            [ipos[1] for ipos in data]

        temp.edge2["index"][:, 0] = indx + 0
        temp.edge2["index"][:, 1] = indx + 1

        nset.append(temp.vert2)
        eset.append(temp.edge2)

    else:

        npts = 1

        temp = jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (npts + 0), dtype=temp.VERT2_t)

        temp.vert2["coord"][:, 0] = data[0]
        temp.vert2["coord"][:, 1] = data[1]

        nset.append(temp.vert2)

    return nobj, last


def loadshp(name, mesh, filt=None):
    """
    LOADSHP: load a shape file into a jigsaw msh_t object.

    """

    if (not isinstance(name, str)):
        raise Exception("Incorrect type: NAME.")

    if (not isinstance(mesh, jigsaw_msh_t)):
        raise Exception("Incorrect type: MESH.")

    nobj = +0; last = +0; nset = []; eset = []

    for feat in fiona.open(name):

        if (filt is not None):
            if (not filt(feat)): continue

        data = feat["geometry"]["coordinates"]

        nobj, last = \
            coreshp(data, nset, eset, nobj, last)

    if (len(nset) > 0):
        mesh.vert2 = np.concatenate(nset, axis=0)
    
    if (len(eset) > 0):
        mesh.edge2 = np.concatenate(eset, axis=0)

    return


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--src_file", dest="src_file", type=str,
        required=True, help="Input shape file to load")

    parser.add_argument(
        "--dst_file", dest="dst_file", type=str,
        required=True, help="Output .msh file to save")

    args = parser.parse_args()

    mesh = jigsaw_msh_t()

    loadshp(name=args.src_file, mesh=mesh)

    savemsh(name=args.dst_file, mesh=mesh)
