------------------------------------------------------------
README - pmm-qmmm_input_generator.py
------------------------------------------------------------

Program requirements
------------------------------------------------------------
This script generates input files for PMM-QMMM simulations
using QM/MM data from ORCA and MD trajectories.

In the working directory, you must have the following
folders for each electronic state:

  S0  -> ground state
  S1  -> first excited state
  S2  -> second excited state
  ...

Each folder must contain the ORCA output file:

  ex_state.out

Python version: 3.6 or higher

Required Python modules: os, sys, numpy

------------------------------------------------------------
Running the program
------------------------------------------------------------

Important warning  
The PMM-QMMM method may produce unreliable results for:
- Triatomic molecules
- Perfectly linear molecules
Please ensure that your system is suitable before proceeding.


Run the program with:

  python3 pmm-qmmm_input_generator.py

You will be prompted to enter:

- The number of electronic excited states (0 = ground state only)
- File paths for PDB, ORCA forcefield, point charges, and MD trajectory
- Initial and final frame indices
- Charge of the Quantum Center (QC)
- Atom indices for the QC and optional QC subregion
- Whether covalent bonds exist between QC and the environment

The program validates all indices and file paths.

------------------------------------------------------------
Required inputs
------------------------------------------------------------
Prompt shown on screen:                       Description
---------------------------------------------------------------------
Enter the PDB file name                       -> PDB structure file
Enter the ORCA forcefield file (*prms)        -> ORCA forcefield parameters
Enter the charge of the Quantum Center (QC)   -> Total QC charge
Enter the ORCA point-charges file             -> ORCA point-charges (.pc.tmp)
Enter the MD trajectory file in XTC format    -> MD trajectory (.xtc)
Enter the initial frame (0 for the first)     -> Starting frame index
Enter the last frame                          -> Last frame index
Enter the indices of atoms in the QC          -> QC atom indices (starting from 0)
Enter indices of atoms for QC subregion       -> Optional subregion for fitting
Enter indices of atoms bonded to QC           -> If covalent bonds exist
Enter indices of QC fragment outside the QC   -> If covalent bonds exist
Enter the number of electronic excited states -> Total number of excited states

------------------------------------------------------------
Output
------------------------------------------------------------
The program generates the following files:

  MD_charges.txt       -> Simulation charges for MD
  QC_indexes.txt       -> QC atom indices with masses
  env_indexes.txt      -> Environment atom indices
  file_QC.txt          -> QM geometry of the QC
  input-pmm-QMMM       -> Complete input file for PMM-QMMM

At the end of execution, you will see:

  The input files have been successfully generated!

------------------------------------------------------------
Notes
------------------------------------------------------------
- All required input files must exist in the current working directory,
  or the full path must be provided.
- Atom indices start from 0.
- For any element not listed in the atomic mass dictionary,
  the user must manually add its atomic mass at the beginning of the script.
- Input indices are validated: duplicates, out-of-range values,
  or non-numeric entries will trigger an error message and
  prompt for re-entry.
- The script requires Python 3.6 or higher, due to the use of f-strings.
- The program is designed for QM/MM workflows using ORCA
  and molecular dynamics trajectory files in XTC format.
- The code is designed to post-process output files from ORCA TD-DFT calculations.
------------------------------------------------------------

