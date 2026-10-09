from ase import units
from ase.constraints import FixAtoms
from ase.io import read
from ase.md.langevin import Langevin
from ase.calculators.plumed import Plumed
from fairchem.core import FAIRChemCalculator, pretrained_mlip
import time

# ---- settings ----
structure_file = "eq_cu100_solv_00_w4Cs_nCO.xyz"
plumed_file = "OPES_calc_nCO_Cs.dat"
timestep = 0.5 * units.fs
temperature_K = 300
friction = 0.002 / units.fs
nsteps = 20_000       # ~10 ps at 0.5 fs
loginterval = 100

start_time = time.time()

atoms = read(structure_file) # Frozen atoms carried on
atoms.pbc = True

print("Fixed atoms:", sum(len(c.index) for c in atoms.constraints))

base_calc = FAIRChemCalculator(
    pretrained_mlip.get_predict_unit("esen-sm-conserving-all-oc25", device="cuda")
)

setup = open(plumed_file).read().splitlines()
atoms.calc = Plumed(
    calc=base_calc,
    input=setup,
    timestep=timestep,
    atoms=atoms,
    kT=units.kB * temperature_K,
)

print(f"Starting production OPES run... ({nsteps * timestep / units.fs / 1000:.1f} ps)")
dyn = Langevin(
    atoms,
    timestep=timestep,
    temperature_K=temperature_K,
    friction=friction,
    trajectory="prod_Cu-Cu.traj",
    logfile="prod_Cu-Cu.log",
    loginterval=loginterval,
)
dyn.run(nsteps)

atoms.write("final_Cu-Cu.xyz")

elapsed = time.time() - start_time
print(f"Production run completed in {elapsed:.2f} seconds")
print("T =", round(atoms.get_temperature(), 1), "K")