import math
import numpy as np
from numba import njit
from dcorr.dump import read_dump, window_iter
from dcorr.lebedev import lebedev_grid


def type_masker(itype=0):
    def masker(frame):
        if itype > 0:
            assert itype <= frame['types'].max()
            mask = frame['types'] == itype
        else:
            mask = np.array([True for i in frame['types']])
        return mask
    return masker


@njit("f8(f8[:,:], f8[:,:])")
def sisf_inner(dpos, qgrid):
    s = 0.0
    for i in range(qgrid.shape[0]):
        for j in range(dpos.shape[0]):
            s += qgrid[i, 3]*math.cos(np.sum(qgrid[i, 0:3]*dpos[j, 0:3]))

    return s


def dynamics_one(window, qgrid, rtol, masker):
    res = []
    rsqtol = rtol*rtol
    for i, w in enumerate(window):
        if i == 0:
            tref = w['time']
            posref = w['positions']
            mask = masker(w)
            nmasked = mask.sum()
            assert nmasked > 0
        tcurr = w['time']
        pos = w['positions']
        dpos = pos - posref
        dpos = dpos[mask]
        t = tcurr - tref
        sisf = sisf_inner(dpos, qgrid)
        dsq = np.square(dpos).sum(axis=1)
        msd = dsq.mean()
        if np.all(dsq == 0):
            alpha2 = 0
        else:
            alpha2 = (3/5)*(np.square(dsq).mean()/msd**2) - 1
        qt = (dsq <= rsqtol).sum()
        res.append([t, sisf, msd, alpha2, qt])
    return nmasked, np.array(res)


def dynamics(dumpfile, nt, ndt, masker=type_masker(), qmax=1.0, nq=6, rtol=1.0, dt=0.002, maxframes=0):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    qgrid = lebedev_grid(nq)
    qgrid[:, 0:3] = qmax*qgrid[:, 0:3]
    ic = 0
    res = []
    for window in window_iter(dump, width=nt, stride=ndt):
        istart = window[0]['index']
        iend = window[-1]['index']
        print("DYNAMICS: {:3d}-th average, [{:5d} to {:5d}].".format(ic, istart, iend))
        nmasked, s = dynamics_one(window, qgrid, rtol, masker=masker)
        res.append(s)
        ic += 1
    t = res[0][:, 0]
    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    for i in res:
        s1 += i[:,1]
        s2 += i[:,2]
        s3 += i[:,3]
        s4 += i[:,4]
    s1 /= len(res)
    s2 /= len(res)
    s3 /= len(res)
    s4 /= len(res)
    s1x4 = 0
    s4x4 = 0
    for i in res:
        s1x4 += np.square(i[:,1]) - np.square(s1)
        s4x4 += np.square(i[:,4]) - np.square(s4)
    s1 /= nmasked
    s4 /= nmasked
    s1x4 /= len(res)*nmasked
    s4x4 /= len(res)*nmasked

    return np.c_[t, s1, s1x4, s2, s3, s4, s4x4]
