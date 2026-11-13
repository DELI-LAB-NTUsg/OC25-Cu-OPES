# Cu dissolution with OC25: package of 2026-11-13

Three parts:

| Folder | What | Status |
|---|---|---|
| `A_constrained_MD/` | constrained-height MD: eSEN-OC25 equilibration (100 ps per window), then VASP (20 ps per window); Cu(100) 4×4, ± CO ± Cs | ready |
| `B_OPES/` | final eSEN-OC25 OPES_METAD_EXPLORE runs on Cu(100), neutral, ± CO, 10 walkers × 12 ns | ready |
| `C_single_point_DFT/` | DFT single points on 33 frames from your previous campaign (Cu(100) 8×8 and Cu(310), ± Cs ± CO) for parity | ready |

---

## A. Constrained MD (`A_constrained_MD/`)

Free energy of a Cu adatom going into solution by thermodynamic integration: the target Cu's height h is held fixed
(its z frozen, x and y free) in a series of windows, and F(h) = −∫⟨F_z⟩ dh, with F_z the force on the frozen z.
For a single Cartesian coordinate no blue-moon correction is needed, so F(h) compares directly with the OPES F(h).

### Systems (`esen/inputs/`)
Cu(100) 4×4, 4 layers (64 Cu, bottom 2 layers = 32 Cu fixed), 10.11 × 10.11 Å, 12 Å of water. The target is an
adatom (h = 1.75 Å) on the hollow farthest from its own vacancy (an adatom–vacancy pair). Each cell was equilibrated
100 ps with eSEN-OC25 at h = 1.75 Å (target z fixed, CO tethered in wCO):

| name | CO | Cs | waters | atoms | target (ASE, 0-based) | C |
|---|---|---|---|---|---|---|
| `4x4ad_nCO_00` | – | – | 41 | 187 | 45 | – |
| `4x4ad_wCO_00` | on the adatom | – | 37 | 177 | 45 | 64 |
| `4x4adcs_nCO_00` | – | 1 (σ = −15.7 µC/cm²) | 37 | 176 | 45 | – |
| `4x4adcs_wCO_00` | on the adatom | 1 | 33 | 166 | 45 | 64 |

The water count differs because waters overlapping the CO and the Cs were removed when the cells were built; the
water surface is 12.3–13.4 Å above the Cu surface in all four, well above the highest window (6.25 Å).

(`EQ_<name>.xyz` = structure with velocities and fixed atoms; `target_<name>.txt` = target, carbon, top-layer Cu.)
h is the target's height above the mean z of the other 15 top-layer Cu.

### Windows
**7 heights, 0.75 Å apart: h = 1.75, 2.50, 3.25, 4.00, 4.75, 5.50, 6.25 Å**, sampled on the way up and on the way down:
13 windows per system (k00–k06 up 1.75 → 6.25 Å, k07–k12 down 5.50 → 1.75 Å), **52 windows in total**.
Why both directions: in our eSEN tests the water shell around the Cu lags behind at h ≈ 2.8–3.3 Å; going up alone
overestimates F, going down underestimates it, and their average matched OPES within 0.02–0.04 eV.

### Constraint and biases (no PLUMED files: set inside `scan.py` for eSEN, written by `make_vasp_inputs.py` for VASP)
| | Applies to | eSEN (step 1) | VASP (step 2) |
|---|---|---|---|
| Target height | all | constraint, not a bias: target z frozen (`FixCartesian`), x and y free | selective dynamics, target `T T F` |
| Bottom 2 Cu layers | all | fixed (`FixAtoms`) | `F F F` |
| CO on the target | wCO only | PLUMED `c_cu: DISTANCE ATOMS=<C>,<target>`, `UPPER_WALLS ARG=c_cu AT=2.3 KAPPA=2000.0` | `ICONST`: `R <C> <target> 4`; `FBIAS_R0 = 2.3`, `FBIAS_A = 3.0`, `FBIAS_D = 100`, i.e. V = A / (1 + exp[−D(ξ/ξ₀ − 1)]) |
| Other biases (OPES, metadynamics) | – | none | none |

The two CO walls differ in shape (PLUMED: harmonic beyond 2.3 Å; VASP: a smooth 3 eV step at 2.3 Å); both are ~0 while
the CO is bound (C–target 1.85–2.15 Å).

### Step 1: eSEN equilibration, one chain per system (`esen/`)
```bash
cd esen
python -u scan.py 4x4ad_nCO_00   100 1.75 updown 0.75 6.25
python -u scan.py 4x4ad_wCO_00   100 1.75 updown 0.75 6.25
python -u scan.py 4x4adcs_nCO_00 100 1.75 updown 0.75 6.25
python -u scan.py 4x4adcs_wCO_00 100 1.75 updown 0.75 6.25
```
- Each chain: 10 ps run-in at 1.75 Å, then 13 windows of 100 ps. Every window starts from the end of the previous one
  with the target moved by 0.75 Å, so the water follows the Cu. The move is gradual: 0.05 Å at a time with 20 fs of
  MD after each (0.3 ps per window change, not recorded), and in wCO the CO (C and O) moves together with the target,
  so the Cu–C bond is never compressed or stretched against the wall.
- eSEN-sm-conserving-all-oc25, Langevin 300 K, friction 5 ps⁻¹, 0.5 fs, normal H masses; target z fixed
  (`FixCartesian`). wCO: PLUMED upper wall on C–target (2.3 Å, κ 2000 eV/Å²; needs py-plumed). Water watchdog attached.
- Outputs: `scan_<name>.csv` (every 2 fs: window, direction, h, F_z, cno, cn, C–Cu, E_pot) and
  `win_<name>_k00…k12.xyz` (last frame of each window, positions + velocities: the VASP starting points).
- **If the water watchdog stops a chain**, restart it as last time: roll back at least 50 ps, new thermostat seed.
  The watchdog prints the MD step at which the water was first broken; with k the last window in `scan_<name>.csv`,
  that is t = step / 2000 ps into window k. Roll back to the end of a window:
  - broken at t ≥ 50 ps: rerun from window k, i.e. `k0 = k`
  - broken at t < 50 ps: rerun from window k − 1, i.e. `k0 = k − 1` (if that is −1, the chain starts again from the
    beginning, run-in included)

  ```bash
  python -u scan.py 4x4ad_wCO_00 100 1.75 updown 0.75 6.25 <k0> <n>     # n = restart number: 1, 2, 3
  ```
  This continues from `win_<name>_k<k0−1>.xyz` with new velocities and thermostat seeds, rewrites the window files
  from k0 on, and logs to `scan_<name>_r<n>.csv` (the first run's csv is kept). Please log each restart (system,
  window, time, printed atom indices) and keep the log of the stopped run. If the same chain breaks three times, stop
  it and let us know.

### Step 2: VASP, 20 ps per window (`vasp/`)
```bash
cd vasp
for n in 4x4ad_nCO_00 4x4ad_wCO_00 4x4adcs_nCO_00 4x4adcs_wCO_00; do python3 make_vasp_inputs.py $n ../esen; done
bash make_potcar.sh          # set PP inside first; Cs -> Cs_sv
```
This makes `vasp/<name>/k<kk>_<up|down>_h<h>/` with POSCAR, INCAR, KPOINTS, (ICONST), `index_map.txt`; each window
starts from its eSEN frame (positions and velocities).
- **20 ps per window**: `NSW = 40000` at `POTIM = 0.5` fs.
- OC25 level of theory (RPBE `GGA = RP`, D3 `IVDW = 11`, 400 eV, σ 0.1, dipole correction, `EDIFF = 1E-5`,
  KPOINTS 4×4×1) + `IBRION = 0`, `MDALGO = 3` (Langevin, `LANGEVIN_GAMMA = 5` per species), `TEBEG = 300`,
  `ISIF = 2`, `ISYM = 0`. Parallelisation tags are left to you.
- Selective dynamics: target `T T F`, the 32 bottom Cu `F F F`, all else free.
- wCO: `ICONST` `R <C> <target> 4` + `FBIAS_R0 = 2.3`, `FBIAS_A = 3.0`, `FBIAS_D = 100` (Fermi-type upper wall on
  C–target; needs a VASP build with `-Dtbdyn`).
- To continue a window: use its CONTCAR as it is (the predictor-corrector block is fine within the same window).

What we checked in short local VASP tests (6.5.0 and 6.6.1) on these cells: selective dynamics with
`MDALGO = 3` keeps the target z fixed; the force on the frozen z is printed in OUTCAR and matched eSEN on the same
frame (−0.067 vs −0.050 and −0.328 vs −0.338 eV/Å); ICONST/FBIAS is accepted with `MDALGO = 3` and acts on the
dynamics, but its force is **not** included in OUTCAR's forces (we correct for it in the analysis from REPORT); the
velocities in the POSCAR are read.

### Please send back
- eSEN: `scan_<name>.csv` (and any `scan_<name>_r*.csv`), `win_<name>_k*.xyz`, logs, restart log.
- VASP, per window: OUTCAR (or vasprun.xml), XDATCAR, OSZICAR, REPORT (wCO: C–target distance every step), CONTCAR,
  INCAR, and the job log.

---

## B. OPES (`B_OPES/`)

### Systems
New cells, built with our structure generator and equilibrated with eSEN-OC25 (FIRE to 0.05 eV/Å, then 100 ps
Langevin NVT at 300 K, 0.5 fs; water checked intact in every cell):

- Cu(100) 8×8, 6 layers (384 Cu; the bottom 2 layers, 128 Cu, fixed via `move_mask` in the xyz), a = 3.5738 Å,
  cell 20.22 × 20.22 Å, 12 Å of water, 12 Å vacuum above it, **neutral (no Cs)**.
- 5 independent starts (different water packings) per condition:
  - nCO: adsorbed H₂O + 163 waters = 164 H₂O, 876 atoms
  - wCO: CO + 163 waters, 875 atoms
- **10 independent walkers**, one per structure: `structures/EQ_cu100_8x8_<nCO|wCO>_<00–04>.xyz`, each with its
  own PLUMED file `plumed/plumed_cu100_8x8_<nCO|wCO>_<00–04>.dat`.

### Target Cu (`targets.csv`)
- wCO: the top-layer Cu nearest the CO carbon after equilibration (the CO equilibrated freely and sits on bridge or
  hollow sites; C–target 1.96–2.15 Å).
- nCO start s uses **the same Cu index as wCO start s** (same lattice site; the slab is written in the same order).
- Targets (PLUMED, 1-based): starts 00–04 → 56, 32, 86, 212, 146. Carbon (wCO): 385.

### PLUMED (one file per walker; no comment lines on purpose)
| | |
|---|---|
| Method | `OPES_METAD_EXPLORE`, BARRIER 3.0 eV, TEMP 300, PACE 200 steps (= 100 fs at 0.5 fs), SIGMA_MIN 0.05, 0.05 |
| Biased CVs | `height.z` (target above the COM of the other 63 top-layer Cu) and `cno` (Cu–water O coordination, R_0 2.6, NN 12, MM 24) |
| Recorded | `cn` (Cu–Cu coordination, R_0 3.0, NN 6, MM 12), `c_cu` (C–target distance, wCO) |
| Lower wall | `height.z` ≥ −0.75 Å, κ 500 eV/Å² (stiff floor, tested: keeps the target out of the slab) |
| Upper walls | `height.z` ≤ 10 Å, κ 10 eV/Å² (water top is 13.4–14.6 Å above the surface); wCO: `c_cu` ≤ 2.3 Å, κ 2000 eV/Å² (CO tethered to the target) |
| States | `State.data`, written every 200 000 steps (= 100 ps) with `STORE_STATES` |
| COLVAR | every 100 steps (50 fs): cn, (c_cu), cno, height.z, opes.*, uwall.bias, lwall.bias |

The walls are separate bias actions (EXPLORE does not take `EXTRA_BIAS`), as in the earlier EXPLORE runs.

### MD
- Model: eSEN-sm-conserving-all-oc25. Langevin NVT, 300 K, **timestep 0.5 fs** (PACE and STATE_WSTRIDE above assume
  0.5 fs), same thermostat and water watchdog as your previous campaign.
- **Length: 12 ns per walker** (24 000 000 steps).
- **Please keep these runs extendable.** Keep the checkpoint (last frame with velocities) and the full OPES state
  (`State.data`, `Kernels.data`) of every walker at 12 ns, so we can continue walkers that are not converged instead
  of starting again. 12 ns is a multiple of `STATE_WSTRIDE` (100 ps), so a state is stored exactly at the end.
- Please keep: COLVAR, Kernels.data, State.data (all stored states), all-atom frames **every 10 ps** (1 ps is not
  needed), MD logs, water checks.

### Water check (`B_OPES/water_check.py`, the same file as in the last package; needs ASE and SciPy)
eSEN can split a water at the leaving Cu (Cu–H in solution + OH⁻ on the surface; 5 events in the last campaign, all
Cu(310) without Cs). We want these runs restricted to intact water:
1. Before launching (all 10 structures print `clean`, checked):
   `python water_check.py structures/EQ_cu100_8x8_nCO_00.xyz`
2. During the run, attached to the dynamics:
   ```python
   from water_check import WaterWatchdog
   dyn.attach(WaterWatchdog(atoms, dyn), interval=2000)   # every 1 ps at 0.5 fs
   ```
   It stops a walker on any H₂, or when water stays broken for 10 consecutive checks (10 ps), and reports the step at
   which the water was first seen broken.
3. Restart rule as last time: roll back at least 50 ps before that step, new thermostat seed, log every break and
   restart; stop a walker after three breaks and let us know.

---

## C. Single-point DFT (`C_single_point_DFT/`)

33 frames picked from your previous campaign (Cu(100) 8×8 4-layer and Cu(310), nCs/4Cs × nCO/wCO, standard OPES),
816–877 atoms each. `frames.csv` gives, per frame, the source walker and time, the target height and why it was chosen:
- per condition (8 conditions): a frame leaving the surface (h 3–4 Å), a dissolved Cu (h 5–7.5 Å), and a surface
  frame (Cu(100) adatom / Cu(310) raised step site) — 24 frames
- trap states (target in a raised site next to a popped-up adatom; fully hydrated Cu at h ≈ 4 Å; Cu perched on a
  raised Cu; CO "pillar") — 4 frames
- lattice references, Cu(100) without Cs, without and with CO — 2 frames
- the eSEN water split on Cu(310) without Cs (2 ps before the split, and two frames at the split): **these frames
  intentionally contain a Cu–H in solution and an OH⁻** — 3 frames

Files:
- `frames/<id>_<walker>_<time>ps/POSCAR` (species sorted, line 6 = element list) and `KPOINTS` (OC25 rule
  round(40/a), round(40/b), 1 → 2×2×1 for every frame)
- `INCAR` (one for all): OC25 level of theory — RPBE (`GGA = RP`), D3 zero damping (`IVDW = 11`), 400 eV, `PREC =
  Normal`, Gaussian smearing σ 0.1, dipole correction along z, `EDIFF = 1E-6`, `ISYM = 0`, `LAECHG = .TRUE.` (for Bader
  charges), `AMIN = 0.01` (as OC25). `ALGO = Normal` (OC25 used All; same converged result).
- `make_potcar.sh`: set `PP` to your PBE PAW folder, then run once in `C_single_point_DFT/`. It copies the INCAR into
  every frame folder and builds the POTCAR (plain PBE for Cu, C, O, H; `Cs_sv` for Cs).
- `frames.txt`: the 33 folders, one per line.

Please send back per frame: OUTCAR (energy, forces), and if possible AECCAR0, AECCAR2 and CHGCAR (Bader).
