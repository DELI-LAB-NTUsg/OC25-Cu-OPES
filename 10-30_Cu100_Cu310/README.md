# Cu dissolution OPES: Cu(100) and Cu(310), with corrected walls

Thank you for the 5 ns OPES vs EXPLORE test. This package contains the structures and PLUMED inputs for the next
runs: two surfaces, Cu(100) and the stepped Cu(310), each without and with CO and without and with 4 Cs, and new
walls that keep the dissolving Cu out of the slab.

## 1. Why the walls change

In the 5 ns test the lower height wall (−1 Å, κ = 10 eV/Å²) did not hold in the EXPLORE arm. `OPES_METAD_EXPLORE`
cannot take the walls into its bias, so the bias kept pushing against the wall:

- the target Cu was below −1 Å for 0.8–3.7 % of the time, down to −1.8 Å, with up to 6.4 eV in the wall and up to
  12–13 Cu neighbours (buried in the slab);
- all 10 water-watchdog events (including the walker that stopped) happened right after the target had been pushed
  under the wall; the water split on the strained target;
- 24 of the 29 lasting adatoms formed while the target was forced under the wall, by a neighbour being pushed up.

Standard OPES has the walls in `EXTRA_BIAS` and never went below −1.2 Å.

We tested a stiff floor (−0.75 Å, κ = 500 eV/Å²) in two short local EXPLORE runs on your 6-layer Cu(100) structures
(nCs_nCO, 4Cs_nCO, 2.7–2.9 ns): the target never went more than 0.03 Å past it (at most 0.34 eV in the wall), the
walkers kept crossing between lattice and solution (more than 3 round trips per ns), and no adatom lasted more than
85 ps. −0.75 Å removes nothing physical: from the reweighted standard-OPES data, a lattice Cu is 0.9 eV up at −0.75 Å.

**Cu(310) needs different walls.** Its (310) atomic planes are only 0.57 Å apart, and three rows are exposed to the
water: the step edge (where the target sits), the terrace middle (~0.6 Å lower) and the step foot (~1.1 Å lower).
A −0.75 Å floor would block real surface sites. From 4 × 90 ps of unbiased MD at 300 K, per row
(18 Cu each; height below the step-edge row; cn as defined in the PLUMED files):

| Site | cn (0.1–99.9 %) | height (0.1–99.9 %) |
|---|---|---|
| step edge (target's row) | 6.6–8.5 | −0.32 to +0.35 Å |
| terrace middle (exposed) | 8.4–10.2 | −0.93 to −0.20 Å |
| step foot (exposed) | 9.4–11.2, max 11.31 | −1.43 to −0.69 Å, min −1.52 |
| first subsurface | 11.1–12.6 | −1.96 to −1.27 Å |

So on Cu(310) two walls are set just past every exposed surface configuration: height ≥ −1.6 Å (κ = 500) and
cn ≤ 11.5 (κ = 50 per cn², so 0.25 past it costs ~3 eV). The height reference is the step-edge row.

| | Cu(100) | Cu(310) |
|---|---|---|
| Height reference `topCu` | the other 63 top-layer Cu | the other 17 step-edge Cu |
| Lower height wall | −0.75 Å, κ = 500 eV/Å² | −1.6 Å, κ = 500 eV/Å² |
| Upper height wall | 10 Å, κ = 10 eV/Å² | 10 Å, κ = 10 eV/Å² |
| Upper cn wall | – | 11.5, κ = 50 |
| C–target Cu wall (wCO) | 2.3 Å, κ = 2000 eV/Å² | 2.3 Å, κ = 2000 eV/Å² |

## 2. Structures

| Surface | Condition | File | Atoms | Water | Cs | CO |
|---|---|---|---:|---:|---:|---|
| Cu(100) | `nCs_nCO` | `EQ_cu100_nCs_nCO.xyz` | 871 | 205 | 0 | no |
| Cu(100) | `nCs_wCO` | `EQ_cu100_nCs_wCO.xyz` | 873 | 205 | 0 | yes |
| Cu(100) | `4Cs_nCO` | `EQ_cu100_4Cs_nCO.xyz` | 875 | 205 | 4 | no |
| Cu(100) | `4Cs_wCO` | `EQ_cu100_4Cs_wCO.xyz` | 877 | 205 | 4 | yes |
| Cu(310) | `nCs_nCO` | `EQ_cu310_nCs_nCO.xyz` | 855 | 177 | 0 | no |
| Cu(310) | `nCs_wCO` | `EQ_cu310_nCs_wCO.xyz` | 845 | 173 | 0 | yes |
| Cu(310) | `4Cs_nCO` | `EQ_cu310_4Cs_nCO.xyz` | 826 | 166 | 4 | no |
| Cu(310) | `4Cs_wCO` | `EQ_cu310_4Cs_wCO.xyz` | 816 | 162 | 4 | yes |

- **Cu(100):** 8 × 8, 4 layers, 256 Cu, bottom 2 layers (128 Cu) fixed; cell 20.23 × 20.23 × 42.36 Å.
  4 Cs: σ ≈ −15.7 µC/cm².
- **Cu(310):** 18 (310) planes × 18 Cu, 324 Cu, bottom 8 planes (144 Cu) fixed; cell 21.46 × 16.96 × 46.94 Å.
  The target is a step-edge Cu. 4 Cs: σ ≈ −17.6 µC/cm².
- Fixed atoms are in the `move_mask` column of each file.
- Already Equilibrated: FIRE to fmax 0.05 eV/Å, then 100 ps Langevin NVT at 300 K, 0.5 fs, friction 0.01/fs, with
  the water watchdog. All eight are water-clean (`water_check.py` prints `clean`).
- Each nCO/wCO pair uses the same target Cu: 74 for Cu(100) nCs, 178 for Cu(100) 4Cs, 311 for all four Cu(310).
  In the Cu(100) CO systems the CO moved during equilibration (no tether there) and now bridges two top-layer Cu;
  the target is the Cu nearest the carbon, and the CO-free partner uses the same Cu. In Cu(310) nCs_wCO the CO
  bridges the target and its neighbour Cu 38 (C–Cu 2.03 and 1.93 Å); in Cu(310) 4Cs_wCO it sits atop the target
  (1.82 Å). The Cs moved during equilibration too, so the nearest Cs to the target differs within the 4Cs pairs
  (Cu(100): 4.2 Å without CO, 7.4 Å with CO).

## 3. PLUMED input

| File | Use |
|---|---|
| `plumed_cu100_nCO.dat` | Cu(100) without CO |
| `plumed_cu100_wCO.dat` | Cu(100) with CO (adds the C–target Cu wall and prints `c_cu`) |
| `plumed_cu310_nCO.dat` | Cu(310) without CO (adds the cn wall) |
| `plumed_cu310_wCO.dat` | Cu(310) with CO |

The files are written for standard `OPES_METAD`, with every wall in `EXTRA_BIAS`. For an EXPLORE arm, change
`OPES_METAD` to `OPES_METAD_EXPLORE` on the `opes:` line and delete the `EXTRA_BIAS=...` keyword (EXPLORE does not
accept it); the walls stay active. Everything else is as in the 5 ns test: biased CVs `height.z` and `cno`,
BARRIER 3 eV, PACE 200 steps (100 fs), SIGMA_MIN 0.05,0.05. Full OPES states are stored every 100 ps
(`STATE_WSTRIDE=200000` with `STORE_STATES`). Fill in the `<...>` placeholders from the table below.

### Atom indices (1-based, PLUMED)

| Structure | Target Cu `<TARGET_CU>` | Other Cu `<OTHER_CU>` | Surface Cu, excl. target `<TOP_CU>` | Water O `<WATER_O>` | CO C `<CARBON>` | CO O | Cs |
|---|---|---|---|---|---|---|---|
| `EQ_cu100_nCs_nCO.xyz` | 74 | 1-73,75-256 | 2-70:4,78-254:4 | 257-869:3 | – | – | – |
| `EQ_cu100_nCs_wCO.xyz` | 74 | 1-73,75-256 | 2-70:4,78-254:4 | 259-871:3 | 257 | 258 | – |
| `EQ_cu100_4Cs_nCO.xyz` | 178 | 1-177,179-256 | 2-174:4,182-254:4 | 257-869:3 | – | – | 872-875 |
| `EQ_cu100_4Cs_wCO.xyz` | 178 | 1-177,179-256 | 2-174:4,182-254:4 | 259-871:3 | 257 | 258 | 874-877 |
| `EQ_cu310_nCs_nCO.xyz` | 311 | 1-310,312-324 | 5-293:18 | 325-853:3 | – | – | – |
| `EQ_cu310_nCs_wCO.xyz` | 311 | 1-310,312-324 | 5-293:18 | 327-843:3 | 325 | 326 | – |
| `EQ_cu310_4Cs_nCO.xyz` | 311 | 1-310,312-324 | 5-293:18 | 325-820:3 | – | – | 823-826 |
| `EQ_cu310_4Cs_wCO.xyz` | 311 | 1-310,312-324 | 5-293:18 | 327-810:3 | 325 | 326 | 813-816 |

Each group was checked against its structure: `otherCu` is every Cu except the target; `topCu` is the full surface
row without the target (Cu(100): the 63 other top-layer Cu; Cu(310): the 17 other step-edge Cu); `waterO` is every
water oxygen and nothing else. At t = 0 every wall is at zero energy: the target starts at height −0.21 to
+0.26 Å, cn 7.5–7.8 on Cu(310) and 8.4–9.5 on Cu(100), and C–target Cu 1.82–2.03 Å in the CO systems.

## 4. Requested runs

| | `nCs_nCO` | `nCs_wCO` | `4Cs_nCO` | `4Cs_wCO` |
|---|---|---|---|---|
| Cu(100) | ✓ | ✓ | ✓ | ✓ |
| Cu(310) | ✓ | ✓ | ✓ | ✓ |

- as many walkers as your compute allows, with the same number for every system, more than one so we can check
  whether independent walkers agree, and separate seeds
- 5 ns per walker, kept extendable (checkpoint, `State.data`, `Kernels.data`)
- 300 K Langevin, 0.5 fs, friction 0.005/fs, esen-sm-conserving-all-oc25, `charge=0`, as in the 5 ns test.
  We are still investigating the performance of UMA-1.2.1 before changing models, so please stay on eSEN for now.

## 5. Requested output

- full COLVAR (every 50 fs) and `Kernels.data`
- the `State.data` with the states every 100 ps, plus the final state
- all-atom trajectory every 1 ps (positions only is enough; this lets us time adatom formation against the COLVAR)
- the exact PLUMED input and run script per walker, and the water-watchdog log

## 6. Water check

`water_check.py` (needs ASE and SciPy), as in the last package.

1. Before launching, on each structure (should print `clean`): `python water_check.py EQ_cu310_nCs_wCO.xyz`
2. During the run, attached to the dynamics:
   ```python
   from water_check import WaterWatchdog
   dyn.attach(WaterWatchdog(atoms, dyn), interval=2000)   # every 1 ps at 0.5 fs
   ```
   It stops a walker on any H2, or when water stays broken for 10 consecutive checks (10 ps).
3. **One change to the restart rule:** please roll back at least 50 ps before the step at which the water was first
   seen broken (in the 5 ns test the rollbacks were 0.5–2.5 ps, and the stopped walker restarted three times from the
   same strained configuration). New thermostat seed, log each restart; stop a walker after three breaks and let us know.
