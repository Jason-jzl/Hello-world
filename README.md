# Logic-based MATLAB model of endothelial signaling in PAH

This repository contains the MATLAB code for class project:
**"A logic-based MATLAB model of endothelial signaling in pulmonary arterial hypertension"**.

## Files

- `simulate_PAH_model.m`  
  Defines the continuous logic ODE system for the endothelial signaling network.
  Nodes include LaminarFlow, OscFlow, KLF2, BMP4, SMAD1_5, BMPR2, PKA, eNOS_Phos, NO,
  and PASMC_Proliferation.

- `run_PAH_project.m`  
  Wrapper script that:
  1. sets input parameters (BMPR2_level, shear_level),
  2. simulates three baseline scenarios (normal, low BMPR2, low shear),
  3. performs a 10×10 parameter scan over BMPR2_level and shear_level,
  4. generates the figures used in the report.


## How to run

1. Open MATLAB and add this folder to the MATLAB path.
2. Open `run_PAH_project.m`.
3. Run the script.  
   It will:
   - call `simulate_PAH_model.m`,
   - simulate the three baseline scenarios,
   - compute steady states,
   - run the parameter scan,
   - and generate all plots.


