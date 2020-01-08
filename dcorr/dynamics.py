import sys
import numpy as np
import numexpr as ne
from dcorr.dump import read_dump, window_iter


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

