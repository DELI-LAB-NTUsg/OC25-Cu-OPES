"""Shared bits for the constrained-z eSEN trial: model, the target and its collective variables.

height h  = z(target) - mean z(other top-layer Cu)            (as PLUMED height.z, target above the top layer)
cno       = sum over water O of (1 - (r/2.6)^12) / (1 - (r/2.6)^24)   (PLUMED COORDINATION R_0=2.6 NN=12 MM=24)
cn        = sum over other Cu of (1 - (r/3.0)^6) / (1 - (r/3.0)^12)   (PLUMED cn, R_0=3.0 NN=6 MM=12)
"""
import numpy as np


def calculator():
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    return FAIRChemCalculator(pretrained_mlip.get_predict_unit("esen-sm-conserving-all-oc25", device="cuda"))


def read_target(path):
    """target_<name>.txt (written at equilibration) -> target, C (or None), other top-layer Cu."""
    d = {}
    for line in open(path):
        k, *v = line.split()
        d[k] = [int(x) for x in v]
    target = d["target_ase"][0]
    top = [i for i in d["top_layer_ase"] if i != target]
    return target, d.get("C_ase", [None])[0], np.array(top)


class CVs:
    def __init__(self, atoms, target, carbon, top):
        sym = atoms.numbers
        self.t, self.c, self.top = target, carbon, top
        self.cu = np.array([i for i in np.where(sym == 29)[0] if i != target])
        o = np.where(sym == 8)[0]
        if carbon is not None:                                  # drop the CO oxygen
            o = o[atoms.get_distances(carbon, o, mic=True) > 1.4]
        self.o = o

    @staticmethod
    def switch(r, r0, n, m):
        x = r / r0
        return np.where(np.isclose(x, 1.0), n / m, (1 - x ** n) / (1 - x ** m)).sum()

    def __call__(self, atoms):
        p = atoms.positions
        h = p[self.t, 2] - p[self.top, 2].mean()
        cno = self.switch(atoms.get_distances(self.t, self.o, mic=True), 2.6, 12, 24)
        cn = self.switch(atoms.get_distances(self.t, self.cu, mic=True), 3.0, 6, 12)
        dc = atoms.get_distance(self.t, self.c, mic=True) if self.c is not None else np.nan
        return h, cno, cn, dc
