# Cu dissolution enhanced-sampling campaign: collaborator data snapshot response

Based on the data provided in the data snapshot, a few concerns have raised about the current process, due to which we suggest testing a few changes to the simulation protocol. Overall, barriers can be calculated from the data, which we recommend doing using the `FES_from_State.py`, which was taken directly from the PLUMED github. For plotting as 1D vs. the height we have used the provided `barriers.py`. 

### Surface Roughness
The nature of re-entry onto the lattice is worrying. Using `adatoms.py`, we analysed the number of adatoms and the surface roughenss, this seems to increase very fast and plateaus later. We wonder whether this roughness is part of physically reasonable behaviour, or whether OPES biasing is so large, that it's breaking the lattice to a level that's not physically reasonable. This question is relevant as the nature of the surface dictates the energetics in part. This is also relevant when, in future, we test different facets. Another part of the problem is that OPES was biasing to also sample the highly negative heights of the target Cu, effectively drilling into the surface and exposing the fixed Cu's. To combat this we offer a few suggested changes. Please note that surface roughess is expected, we hope to find out if this converges to a set value, and want to ensure that the roughness happens to a physical level.
#### Suggested changes for adatoms:
- Enforce a lower wall on the height (suggested -1A) to ensure the atom does not go into the second layer.
- New files provided (in `nCs_nCO` and `nCs_wCO`), which have 6 Cu layers, in stead of 4, while keeping the bottom 2 fixed again, in order to mitigate the fixed atoms becoming uncovered.

### New Frontier
Alongside testing how the system can be optimised to ensure the surface is not unphysically rough, we suggest a new frontier to run, which includes four caesium atoms to induce surface charging. This is to see how the barrier changes in charged conditions, matching Vavra, et al.
- Does a negative surface charge impact dissolution?

## Things to note
- New configurations (adding Caesium, moving to 6 layers) change the identies of the atoms that need to be biased, please follow the table below (all 1-indexed for PLUMED use):

| Object | `EQ_4Cs_nCO.xyz` | `EQ_4Cs_wCO.xyz` | `EQ_nCs_nCO.xyz` | `EQ_nCs_wCO.xyz` |
| ------ | ---------- | ---------- | ---------- | ---------- |
| **C**-O | - | 289 | - | 289 |
| C-**O** | - | 290 | - | 290 |
| Target Cu | 185 | 185 | 77 | 77 |
| Other Cu | 1-184, 186-288 | 1-184, 186-288 | 1-76, 78-288 | 1-76, 78-288 |
| topCu | 5-179:6, 191-287:6 | 5-179:6, 191-287:6 | 5-71:6, 83-287:6 | 5-71:6, 83-287:6 |
| WaterO | 291, 292-718:3 | 291-717:3 | 291, 292-718:3 | 291-717:3 |
| Cs | 721-724 | 720-723 | - | - |

- The 8A upper wall seems low, as we do want to try and simulate real solvation dynamics so that barrier can be increased to 12A.
- All .xyz files have been equilibrated. 
- The lower wall at -1A is of an untested magnitude, so please feel free to change it if it hinders effective, natural and physical dissolution.

## Format Requests
Most of what was sent was very good, including the complete State and COLVAR files, for calculating the barriers. However, the following additions are wished:
- Direct copies of the OPES calculator files. (`OPES_[...].dat`)
- Some trajectories are necessary with a sampling of more than 50 ps, with multiple events, so that we can study how the other surface atoms move during events. Not all runs have to be this way, but a variety of runs that, after equilibration, show multiple events with a higher resolution are wished.
- Some runs need to be longer aswell, in the previous runs taking 7 ns we see some apparent plateauing of the roughness and the number of adatoms, we want to see if longer runs (e.g. towards 12 ns if possible), this plateauing still takes place. If it is a plateau, this would likely mean the energetics only converge after a level of disorder is applied to the lattice, likely leading to long equilibration time in future.