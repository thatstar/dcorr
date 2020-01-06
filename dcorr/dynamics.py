import math
import numpy as np
from numba import njit
from dcorr.dump import read_dump, window_iter
from dcorr.lebedev import lebedev_grid


@njit("f8(f8[:,:], f8[:,:])")
def sisf_inner(dpos, qgrid):
    s = 0.0
    for i in range(qgrid.shape[0]):
        for j in range(dpos.shape[0]):
            s += qgrid[i, 3]*math.cos(np.sum(qgrid[i, 0:3]*dpos[j, 0:3]))

    return s/dpos.shape[0]


def sisf_one(window, qgrid, itype=0):
    res = []
    for i, w in enumerate(window):
        if i == 0:
            tref = w['time']
            posref = w['positions']
        t = w['time']
        pos = w['positions']
        dpos = pos - posref
        if itype > 0:
            assert itype <= w['types'].max()
            mask = w['types'] == itype
            dpos = dpos[mask]
        res.append([t - tref, sisf_inner(dpos, qgrid)])

    return np.array(res)


def sisf(dumpfile, nt, ndt, itype=0, qmax=1.0, nq=6, dt=0.002, maxframes=0):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    qgrid = lebedev_grid(nq)
    qgrid[:, 0:3] = qmax*qgrid[:, 0:3]
    ic = 0
    res = []
    for window in window_iter(dump, width=nt, stride=ndt):
        istart = window[0]['index']
        iend = window[-1]['index']
        print("SISF: {:3d}-th average, [{:5d} to {:5d}].".format(ic, istart, iend))
        s = sisf_one(window, qgrid, itype)
        res.append(s)
        ic += 1
    t = res[0][:, 0]
    s = 0
    for i in res:
        s += i[:,1]
    s /= len(res)

    return np.c_[t, s]
