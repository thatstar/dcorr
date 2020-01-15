import sys
import numpy as np
import numexpr as ne
from dcorr.dump import read_dump, window_iter
from dcorr.voronoi import voronoi_indices


def type_masker(itype=0):
    def masker(frame):
        if itype > 0:
            assert itype <= frame['types'].max()
            mask = frame['types'] == itype
        else:
            mask = np.array([True for i in frame['types']])
        return mask
    return masker


def dynamics_one(window, mask, ncorr, qmax, rsqtol):
    res = np.empty((ncorr - 1, 5))
    res.fill(np.nan)
    for i in range(1, len(window)):
        dt = window[i]['time'] - window[0]['time']
        dpos = window[i]['positions'] - window[0]['positions']
        dpos = dpos[mask]
        sisf = ne.evaluate("cos(qmax*dpos)").mean(axis=1).sum()
        dsq = np.square(dpos).sum(axis=1)
        qt = ne.evaluate("dsq <= rsqtol").sum()
        msd = dsq.mean()
        alpha2 = (3/5)*(np.square(dsq).mean()/msd**2) - 1
        res[i - 1, 0] = dt
        res[i - 1, 1] = sisf
        res[i - 1, 2] = qt
        res[i - 1, 3] = msd
        res[i - 1, 4] = alpha2

    return res


def dynamics(dumpfile, ncorr, nshift, masker=type_masker(), qmax=1.0, rtol=1.0, dt=0.002, maxframes=0):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    rsqtol = rtol*rtol
    ic = 0
    s = []
    for window in window_iter(dump, width=ncorr, stride=nshift):
        istart = window[0]['index']
        iend = window[-1]['index']
        print("DYNAMICS: {:3d}-th average, [{:5d} to {:5d}].".format(ic, istart, iend))
        sys.stdout.flush()
        mask = masker(window[0])
        nmasked = mask.sum()
        assert nmasked > 0
        s.append(dynamics_one(window, mask, ncorr, qmax, rsqtol))
        ic += 1
    s = np.stack(s, axis=2)
    res = np.zeros((ncorr - 1, 7))
    res[:, 0] = s[:, 0, 0]
    res[:, 1] = np.nanmean(s[:, 1, :], axis=1)/nmasked
    res[:, 2] = (np.nanmean(s[:, 1, :]**2, axis=1) - np.nanmean(s[:, 1, :], axis=1)**2)/nmasked
    res[:, 3] = np.nanmean(s[:, 2, :], axis=1)/nmasked
    res[:, 4] = (np.nanmean(s[:, 2, :]**2, axis=1) - np.nanmean(s[:, 2, :], axis=1)**2)/nmasked
    res[:, 5] = np.nanmean(s[:, 3, :], axis=1)
    res[:, 6] = np.nanmean(s[:, 4, :], axis=1)

    return res


def find_window_width(dumpfile, X4time, dt, maxframes):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    t = np.array([abs(d['time'] - X4time) for d in dump])

    return int(t.argmin() + 1)


def mobility(dumpfile, width, nshift=1, masker=type_masker(), nbins=100, dt=0.002, maxframes=0):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    dsq = []
    for window in window_iter(dump, width=width, stride=nshift):
        if len(window) == width:
            dpos = window[-1]['positions'] - window[0]['positions']
            dsq.append(np.square(dpos[masker(window[0]), :]).sum(axis=1))
    dsq = np.concatenate(dsq)
    bins = np.logspace(np.log10(1e-2), np.log10(dsq.max()), nbins + 1)
    h, bins = np.histogram(dsq, bins=bins)
    p = 0.5*(bins[1::] + bins[0:-1])
    h = 100*h/len(dsq)
    c = np.cumsum(h)
    return np.c_[p, h, c]


def mobility_analysis(dumpfile, width, slowcut, fastcut, nshift=1, masker=type_masker(), dt=0.002, maxframes=0):
    dump = read_dump(dumpfile, maxframes=maxframes, dt=dt)
    voroinds = []
    mask0 = []
    mask1 = []
    ic = 0
    for window in window_iter(dump, width=width, stride=nshift):
        if len(window) == width:
            istart = window[0]['index']
            iend = window[-1]['index']
            print("Voronoi: {:3d}-th window, [{:5d} to {:5d}].".format(ic, istart, iend))
            dpos = window[-1]['positions'] - window[0]['positions']
            ma = masker(window[0])
            dsq = np.square(dpos[ma, :]).sum(axis=1)
            mask0.append(dsq <= slowcut)
            mask1.append(dsq >= fastcut)
            voroinds.append(voronoi_indices(window[0])[ma, :])
            ic += 1
    voroinds = np.vstack(voroinds)
    mask0 = np.hstack(mask0)
    mask1 = np.hstack(mask1)
    va, ca = np.unique(voroinds, axis=0, return_counts=True)
    res = {}
    for i in range(va.shape[0]):
        key = "{}-{}-{}-{}".format(int(va[i,0]), int(va[i,1]), int(va[i,2]), int(va[i,3]))
        res[key] = [ca[i]/ca.sum(), 0, 0]
    v0, c0 = np.unique(voroinds[mask0, :], axis=0, return_counts=True)
    for i in range(v0.shape[0]):
        key = "{}-{}-{}-{}".format(int(v0[i,0]), int(v0[i,1]), int(v0[i,2]), int(v0[i,3]))
        res[key][1] = c0[i]/mask0.sum()
    v1, c1 = np.unique(voroinds[mask1, :], axis=0, return_counts=True)
    for i in range(v1.shape[0]):
        key = "{}-{}-{}-{}".format(int(v1[i,0]), int(v1[i,1]), int(v1[i,2]), int(v1[i,3]))
        res[key][2] = c1[i]/mask1.sum()
    print("Numbers: {} {} {}".format(ca.sum(), mask0.sum(), mask1.sum()))

    return res
