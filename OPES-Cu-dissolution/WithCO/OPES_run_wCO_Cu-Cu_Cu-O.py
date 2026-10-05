from ase import units
from ase.constraints import FixAtoms
from ase.io import read
from ase.md.langevin import Langevin
from ase.calculators.plumed import Plumed
from fairchem.core import FAIRChemCalculator, pretrained_mlip
import time

# ---- settings ----
structure_file = "EQ_wCO.xyz"
plumed_file = "OPES_calc_wCO_Cu-Cu_Cu-O.dat"
z_fixed = 15.0          # bottom two Cu layers fixed
timestep = 0.5 * units.fs
temperature_K = 300
friction = 0.002 / units.fs
nsteps = 1_000_000       # ~500 ps at 0.5 fs
loginterval = 100

start_time = time.time()

atoms = read(structure_file)
atoms.pbc = True

z = atoms.positions[:, 2]
atoms.set_constraint(FixAtoms(mask=z < z_fixed))

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
    trajectory="prod_wCO_Cu-Cu_Cu-O.traj",
    logfile="prod_wCO_Cu-Cu_Cu-O.log",
    loginterval=loginterval,
)
dyn.run(nsteps)

atoms.write("final_wCO_Cu-Cu_Cu-O.xyz")

elapsed = time.time() - start_time
print(f"Production run completed in {elapsed:.2f} seconds")
print("T =", round(atoms.get_temperature(), 1), "K")