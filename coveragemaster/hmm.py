"""Hidden Markov Model + Haar wavelet engine for CNV calling."""
import sys

import pywt
from numpy import (
    argmax, array, ceil, median, ones, repeat, where, zeros, around,
)
from scipy.stats import norm

from coveragemaster.log import logreport


class _Observations:
    def __init__(self, o_states):
        self.o_states = o_states


def _emission(obs, O, nstd=7):
    """Emission probability vector for one observation value."""
    if 0.9 < obs < 1.1:
        return array([0.5, 0.0, 0.0])
    os1 = obs if obs < 2 else 2
    res = []
    for state in O.o_states:
        if state == max(O.o_states):
            p = 0.5 - abs(norm(state, nstd).cdf(os1) - 0.5) if state < 2 else 0.5
        else:
            p = 0.5 - abs(norm(state, nstd).cdf(os1) - 0.5)
        res.append(p)
    return array(res)


def _downsample(signal, lev):
    if lev == 0:
        return signal
    ds = pywt.downcoef("a", signal + 1e-5, "haar", level=lev)
    ds = ds / median(ds) * median(signal)
    return ds


def _upsample_trigger(signal, ln):
    signal = around(signal, 1)
    signal_pos = where(signal != 1)[0]
    try:
        signal[signal_pos - 1] = 0.5
        signal[signal_pos + 1] = 0.5
    except IndexError:
        pass
    return repeat(signal.tolist(), ceil(ln / len(signal)))


def _upsample(signal, ln):
    try:
        signal = signal.tolist()
    except AttributeError:
        pass
    return repeat([1] + signal, ceil(ln / len(signal)))


def HMM_long(signal, std_o, gene, booster=1, lev=5, LOGFILE="", mask=None, minlev=0, orgmedian=1):
    """Multi-resolution Viterbi HMM with wavelet zoom.

    States: 1.0 = diploid normal, 1.5 = duplication, 0.5 = deletion.
    """
    if len(signal) > 1e5:
        lev = 12

    states = [1, 1.5, 0.5]
    Obs = _Observations(states)
    n = len(states)
    trans = zeros((n, n))
    a = 4.5e-6
    for i in range(n):
        for j in range(n):
            trans[i, j] = a * booster if i != j else 1 - a * booster

    v = (1 / 3, 1 / 3, 1 / 3)
    ZOOM = True
    BEGIN = True
    trigger = None

    while ZOOM:
        pos_states = _downsample(signal, lev)
        std = _downsample(std_o, lev)
        path_max = []
        _em_v = [_emission(pos_states[t], Obs, std[t]) for t in range(len(pos_states))]

        if BEGIN:
            trigger = array([
                1.5 * (em[0] < em[1]) or 0.5 * (em[0] < em[2]) or 1
                for em in _em_v
            ], dtype=float)
            BEGIN = False

        if not sum(abs(trigger - 1)):
            trigger = _upsample(trigger, len(signal))[: len(signal)]
            logreport(f"{gene}: nothing there", logfile=LOGFILE)
            return trigger, states

        for t in range(1, len(pos_states)):
            _v = v * trans
            _m = argmax(_v, 1)
            path_max.append(_m)
            _em = _em_v[t]
            v = array([_v[i, _m[i]] * _em[i] for i in range(len(_m))])
            v = v / max(v)

        M = argmax(v)
        i = len(pos_states) - 2
        final_path = [states[M]]
        while i > -1:
            M = path_max[i][M]
            final_path.append(states[M])
            i -= 1

        approx = _upsample(final_path[::-1], len(signal))[: len(signal)]

        if mask is not None:
            trigger = _upsample_trigger(trigger, len(signal))[: len(signal)]
            if sum(trigger * mask) and sum((approx - 1) * trigger * mask) == 0:
                signal = signal * (trigger != 1) + 1 * (trigger == 1)
                lev -= 1
                if lev < minlev:
                    break
                logreport(f"{gene}: Zooming level {lev}", logfile=LOGFILE)
            else:
                ZOOM = False
        else:
            break

        if lev > 6:
            break

    return approx, states
