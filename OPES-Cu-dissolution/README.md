# OC25 Application on Cu dissolution

### Research Summary

Testing the mechanism proposed by Vavra, et al. in [Nature Catalysis](https://www.nature.com/articles/s41929-023-01070-8), where CO adsorption to a Cu electrode weakens Cu-Cu bonding and drives Cu dissolution as a [CuCO]+ complex during the CO2 reduction reaction. Here, we use OC25 with OPES-METAD enhanced sampling to test the mechanism and quantify the barrier. A single Cu(100) is biased towards dissolution, either with or without a CO adsorbed. The possible biases are:
- Only the Cu-Cu CN
- Both the Cu-Cu and the Cu-O CN


Here, input files are provided for sampling the same Cu both with and without CO adsorbed and for biasing either only the Cu-Cu CN or both the Cu-Cu and the Cu-O CN as CVs.

Currently, we are hoping to see whether variation of the barrier and the pace can lead to more dissolution and readsorption events. We also want to test whether biasing only the Cu-Cu CN or both Cu-Cu and Cu-O CN works better. The first effects we've seen on our tests show that biasing only Cu-Cu is slightly better, but it is interesting to see how this image changes over larger timescales.

Untill now we've currently always set the pass at 200 (100fs). We've found a barrier of 0.5 eV as too low and have been running at 2 eV since. 

### Goals:
- Test a range of barriers (e.g. 0.5-5eV)
- Test a range of pace (e.g. 50fs-500fs)
- Obtain a large number where the target Cu is both in solution (Cu-Cu CN $\approx0$) and in the lattice (Cu-Cu CN $\approx8.5$), every event going from one configuration to the other helps quantify the barrier. The target Cu may spend a long time on the surface (Cu-Cu CN $\approx 4$).

#### Practical notes:
- Barrier and PACE live in the .dat calculator files (BARRIER=, PACE=), not as run-script arguments.
- Filenames are currently fixed per CV choice (prod_Cu-Cu.traj, COLVAR_Cu-Cu, Kernels_Cu-Cu.data, ...), filenames need to be changed or separate directories made to separate runs.



### Requirments
Requires Python, ASE, fairchem OC25 and PLUMED METAD.
GPU required using device="cuda".

### File description
- `COLVAR_Check.py` file in which the COLVAR file and a title can be specified to plot Cu-Cu CN histogram, Cu-Cu CN vs time, target Cu height vs time and OPES Bias vs time.
#### /WithCO
In the runs with CO, the distance between the target Cu and the Carbon is limited to 2.3 A. 
- `EQ_wCO.xyz` equilibriated starting structure with CO.
- `OPES_calc_wCO_Cu-Cu.dat` OPES calculator file with only Cu-Cu CN as CV, writes kernels in `Kernels_Cu-Cu.data`, Sigma file in `State_Cu-Cu.data` and collective variables log (incl. Cu-Cu CN, Cu-O CN, Cu-C distance and target Cu height) in `COLVAR_Cu-Cu`
- `OPES_run_wCO_Cu-Cu.py` Python file running on above structure and OPES calculator in NVT for 500ps. Writes trajectory file in `prod_Cu-Cu.traj`, log in `prod_Cu-Cu.log` and the final configuration in `final_Cu-Cu.xyz`.
- `OPES_calc_wCO_Cu-Cu_Cu-O.dat` OPES calculator file with both Cu-Cu CN and Cu-O as CVs, writes kernels in `Kernels_Cu-Cu_Cu-O.data`, Sigma file in `State_Cu-Cu_Cu-O.data` and collective variables log (incl. Cu-Cu CN, Cu-O CN, Cu-C distance and target Cu height) in `COLVAR_Cu-Cu_Cu-O`
- `OPES_run_wCO_Cu-Cu_Cu-O.py` Python file running on above structure and OPES calculator in NVT for 500ps. Writes trajectory file in `prod_wCO_Cu-Cu_Cu-O.traj`, log in `prod_wCO_Cu-Cu_Cu-O.log` and the final configuration in `final_wCO_Cu-Cu_Cu-O.xyz`.

#### /WithoutCO
- `EQ_nCO.xyz` equilibriated starting structure without CO.
- `OPES_calc_nCO_Cu-Cu.dat` OPES calculator file with only Cu-Cu CN as CV, writes kernels in `Kernels_Cu-Cu.data`, Sigma file in `State_Cu-Cu.data` and collective variables log (incl. Cu-Cu CN, Cu-O CN and target Cu height) in `COLVAR_Cu-Cu`
- `OPES_run_nCO_Cu-Cu.py` Python file running on above structure and OPES calculator in NVT for 500ps. Writes trajectory file in `prod_nCO_Cu-Cu.traj`, log in `prod_nCO_Cu-Cu.log` and the final configuration in `final_nCO_Cu-Cu.xyz`.
- `OPES_calc_nCO_Cu-Cu_Cu-O.dat` OPES calculator file with both Cu-Cu CN and Cu-O as CVs, writes kernels in `Kernels_Cu-Cu_Cu-O.data`, Sigma file in `State_Cu-Cu_Cu-O.data` and collective variables log (incl. Cu-Cu CN, Cu-O CN and target Cu height) in `COLVAR_Cu-Cu_Cu-O`
- `OPES_run_nCO_Cu-Cu_Cu-O.py` Python file running on above structure and OPES calculator in NVT for 500ps. Writes trajectory file in `prod_nCO_Cu-Cu_Cu-O.traj`, log in `prod_nCO_Cu-Cu_Cu-O.log` and the final configuration in `final_nCO_Cu-Cu_Cu-O.xyz`.