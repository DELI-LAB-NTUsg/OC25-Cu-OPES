# Cu dissolution OPES: fixed structures and a short OPES method test

Thank you for the 10.25 ns six-layer data. Two things came out of it, and this package addresses both.

## 1. The previous starting structures had broken water

All four `EQ_*.xyz` structures we sent on 2026-10-09 already contained H2 and OH- at t = 0
(1 H2 + 2 OH- per cell; 2 H2 + 4 OH- in 4Cs_nCO). The water had been packed at far too high a density
in our structure generation (nearest O–O about 1.65 Å instead of about 2.8 Å), and the geometry
optimisation turned the overlapping molecules into H2 + OH-. No new H2 formed during your OPES runs.
The OH- bound to the target Cu in many frames, so those runs, especially the Cs ones, are not usable for barriers.

The generator is fixed, and the new structures in this package were checked after packing,
relaxation and equilibration: no H2, OH- or H3O+, and no overlapping atoms.

| Condition | File | Atoms | Water | Cs | CO |
|---|---|---:|---:|---:|---|
| `nCs_nCO` | `nCs_nCO_EQ.xyz` | 723 | 145 | 0 | no |
| `nCs_wCO` | `nCs_wCO_EQ.xyz` | 722 | 144 | 0 | yes |
| `4Cs_nCO` | `4Cs_nCO_EQ.xyz` | 727 | 145 | 4 | no |
| `4Cs_wCO` | `4Cs_wCO_EQ.xyz` | 726 | 144 | 4 | yes |

Six Cu(100) layers, 8×6, 288 Cu; the bottom two layers (96 Cu) are fixed. Equilibration: FIRE relaxation
to fmax 0.05 eV/Å, then 100 ps Langevin NVT at 300 K, 0.5 fs, friction 0.01/fs, with esen-sm-conserving-all-oc25.

## 2. Request: 5 ns test of standard OPES vs OPES_EXPLORE

The second campaign sampled about 8× fewer basin transitions than the first (about 3 vs about 25 per ns per
walker, counted the same way in both). The most likely main reason is the OPES variant: the first campaign
used `OPES_METAD_EXPLORE` and the second `OPES_METAD`. Standard OPES converges its bias quickly and then
stops pushing, so walkers stayed in one region for ns at a time and the walkers of one condition disagree on the barrier by
up to about 1 eV. Before choosing a method for the production runs, we would like to compare both directly:

| | `OPES_METAD` | `OPES_METAD_EXPLORE` |
|---|---|---|
| `nCs_nCO` | 5 ns per walker | 5 ns per walker |
| `nCs_wCO` | 5 ns per walker | 5 ns per walker |
| `4Cs_nCO` | 5 ns per walker | 5 ns per walker |
| `4Cs_wCO` | 5 ns per walker | 5 ns per walker |

Please run all four conditions with both methods. Use as many walkers as your compute allows, with the same
number in every arm, and more than one so we can check whether independent walkers agree.

**Please keep these runs extendable.** Keep the checkpoint and the full OPES state (`State.data`, `Kernels.data`)
of every walker at 5 ns. We will then continue the better method from 5 ns to the production length instead of
starting again, so the test runs become the first part of the production runs.

Apart from the OPES action, the settings are those of your second campaign, with one change:
the upper height wall moves from 12 Å to 10 Å. In your data the free energy is flat above about 6 Å, and about
25% of the time was spent above 7.5 Å, so this saves sampling without losing information.

- 300 K Langevin, 0.5 fs, friction 0.005/fs, esen-sm-conserving-all-oc25, `charge=0`
- biased CVs `height.z` and `cno`; BARRIER 3 eV, PACE 200 steps (100 fs), SIGMA_MIN 0.05,0.05
- walls: height −1 Å (lower) and **10 Å** (upper), κ = 10; CO systems also the C–target Cu upper wall at 2.3 Å, κ = 2000
- separate random seeds per walker

### PLUMED input

- `plumed_nCO.dat`: systems without CO
- `plumed_wCO.dat`: systems with CO (adds the C–target Cu wall and prints `c_cu`)

Both are written for the standard arm. For the explore arm, change `OPES_METAD` to `OPES_METAD_EXPLORE`
on the `opes:` line; nothing else changes. Fill in the `<...>` placeholders from the table below (1-based).

### Atom indices (1-based, PLUMED)

| Group | `nCs_nCO` | `nCs_wCO` | `4Cs_nCO` | `4Cs_wCO` |
|---|---|---|---|---|
| Target Cu `<TARGET_CU>` | 275 | 275 | 5 | 5 |
| Other Cu `<OTHER_CU>` | 1-274,276-288 | 1-274,276-288 | 1-4,6-288 | 1-4,6-288 |
| Top-layer Cu, excl. target `<TOP_CU>` | 5-269:6,281-287:6 | 5-269:6,281-287:6 | 11-287:6 | 11-287:6 |
| Water O `<WATER_O>` | 291,292-721:3 | 291-720:3 | 291,292-721:3 | 291-720:3 |
| CO carbon `<CARBON>` | – | 289 | – | 289 |
| CO oxygen | – | 290 | – | 290 |
| Cs | – | – | 724-727 | 723-726 |

The atom order differs from the 2026-10-09 structures, so please use these indices rather than the previous ones.
In the CO systems the target is the Cu nearest the CO carbon after equilibration. In the systems without CO the
target is the same lattice site as in the CO partner, so each nCO/wCO pair starts from an equivalent site.
In the systems without CO, water O 291 is the adsorbed water (atoms 289–291, H H O); it is included in `waterO`.

## Requested output

Same as the last bundle, which worked very well:

- full COLVAR (every 50 fs), `Kernels.data`
- complete `State` snapshots every 0.5 ns, plus the final state
- all-atom trajectory every 20 ps
- the exact PLUMED input and run script per walker

## Water check

`water_check.py` (included, needs ASE and SciPy) checks that the water stays intact.

1. Before launching, on each starting structure. It should print `clean`:
   ```
   python water_check.py EQ_nCs_wCO.xyz
   ```
2. During the run, please attach the watchdog to the dynamics in your driver:
   ```python
   from water_check import WaterWatchdog
   dyn.attach(WaterWatchdog(atoms, dyn), interval=2000)   # every 1 ps at 0.5 fs
   ```
   It stops the walker (RuntimeError) as soon as an H2 forms, or when water has been broken (free H, OH-, H3O+)
   in 10 consecutive checks, i.e. 10 ps. Single-frame proton hops between waters do not trigger it.
   It prints the MD step at which the water was first seen broken.

   If a walker stops, please restart it from just before that step: roll the walker back to an earlier
   point where the water is intact (truncating COLVAR, Kernels, State and trajectory to that time, as you did for
   the 4Cs_wCO/w2 restart), and continue with a new thermostat seed. Please log each restart (walker, time,
   printed atom indices) and include the stopped frame in the bundle. If the same walker breaks three times,
   stop it and let us know: that would point to systematic chemistry at the dissolving Cu rather than a one-off.
3. Optionally, on the 20 ps trajectory afterwards: `python water_check.py --every 5 <walker>_every20ps.extxyz.gz`

## What we will compare, and next steps

For each arm: barrier vs simulation time from the 0.5 ns states, agreement between walkers,
transitions per ns, and surface roughness / adatoms. We hope to have this done by the end of this week. The walkers
of the chosen method would then be continued from 5 ns to the production length, and we may add runs with
different numbers of Cs.
