------------------------------------------------------------
README – PMM-QM/MM
------------------------------------------------------------

Program requirements
------------------------------------------------------------
PMM-QM/MM is a Fortran program implementing the Perturbed Matrix Method (PMM)
for QM/MM calculations on molecular systems along MD trajectories.

The program allows both:

- Standard PMM calculations: unperturbed states are computed in a gas-phase QM calculation
  (This is achieved by setting the point-charges file to contain only a single value "0")

- PMM-QM/MM calculations: unperturbed states are computed in the presence of environment point charges,
  i.e., a QM/MM calculation in the additive scheme using electrostatic embedding

Based on the selection made in pmm-qmmm_input_generator.py, the program computes:

- Perturbed electric dipole moment and perturbed electronic energy
  (if "Perturbed Electronic State Properties" is selected)

- Perturbed transition electric dipole moments and perturbed arrival-state energies
  (if "Perturbed Properties for Spectral Absorption Bands")

For compilation, the following are required:

- Fortran 90 compiler (tested with gfortran)
- LAPACK library (liblapack.a)
- BLAS library (librefblas.a)
- XDR file reader for XTC trajectories (libxdrf.a)

Compilation example:

  gfortran PMM-QMMM.f90 liblapack.a librefblas.a libxdrf.a -o PMM-QMMM.x

------------------------------------------------------------
Running the program
------------------------------------------------------------

Run the program with:

  ./PMM-QMMM.x

The program reads the file input-pmm-QMMM, typically generated automatically
by pmm-qmmm_input_generator.py, and computes frame-by-frame perturbed QM/MM
properties along the MD trajectory.

------------------------------------------------------------
Required inputs
------------------------------------------------------------
The input-pmm-QMMM file must contain:

- ORCA point-charges file (.pc.tmp)
  (If the file contains only the value "0", a standard PMM calculation is performed;
   otherwise, a PMM-QM/MM calculation is carried out.)
- QM geometry file (QC atoms only)
- MD trajectory file in XTC format
- Initial frame index (0 = first frame)
- Last frame index
- Number of atoms in the system
- MD simulation charges file
- QC index file with masses
- Environment index file (all atoms outside QC)
- Charge of the Quantum Center (QC)
- Number of electronic states (including ground state)
- Number of atoms in the QC
- Electric dipole moment matrices (x, y, z components)
- Ground-state and excited-state energies
- Indices for the calculation type:
  (0 = Perturbed Electronic State Properties,
   1 = Perturbed Properties for Spectral Absorption Bands)

------------------------------------------------------------
Output
------------------------------------------------------------

The program produces the following output files for each electronic state:

1) output_Si.xvg

Contains, for each MD frame:

frame   Energy(a.u.)   mux(a.u.)   muy(a.u.)   muz(a.u.)   |mu|(a.u.)

Interpretation:

If "Perturbed Electronic State Properties" is selected:

Energy (a.u.) = perturbed energy of the i-th electronic state

mux, muy, muz, |mu| = components and magnitude of the perturbed electric dipole moment of the i-th electronic
                      state

If "Perturbed Properties for Spectral Absorption Bands" is selected:

Energy (a.u.) = perturbed energy of the i-th arrival (excited) electronic state

mux, muy, muz, |mu| = components and magnitude of the perturbed transition dipole moment from the ground state
                      to the i-th electronic state

All energy and dipole values are expressed in atomic units.

------------------------------------------------------------

2) dipole_fluctuations_Si.xvg

Contains dimensionless ratios describing how the perturbed dipole of the i-th elecetronic state deviates from
the reference (unperturbed) QM dipole:

frame   mux/mux0   muy/muy0   muz/muz0

Where:

mux/mux0 = perturbed dipole component / reference dipole component

The reference dipole (mux0, muy0, muz0) is taken from the QM or QM/MM calculation
used as the basis to construct the unperturbed Hamiltonian.

The definition adapts automatically:

If "Perturbed Electronic State Properties" is selected, the ratios refer to transition dipole moments from the
ground state to the i-th electronic state.

If "Perturbed Properties for Spectral Absorption Bands" is selected, the ratios refer to the i-th electronic
state dipole moments.

------------------------------------------------------------

3) unpert_state_coefficients_to_Si_pert.xvg

Contains the coefficients describing the expansion of the i-th perturbed electronic
state over the unperturbed (gas-phase) electronic states:

|Ψ_perturbed⟩_i = Σ_l c_{i,l} |Ψ_l^0⟩

For each MD frame, the file contains:

frame   c_{i,0}   c_{i,1}   c_{i,2}   ...

where c_{i,l} is the coefficient of the l-th unperturbed electronic state
contributing to the i-th perturbed electronic state. The number of coefficients
depends on the number of unperturbed states included in the calculation.

These coefficients quantify the contribution of each unperturbed electronic state to the selected i-th perturbed
electronic state and can be used to analyze state mixing induced by the environment.

------------------------------------------------------------
Workflow
------------------------------------------------------------
1. Initialization
   - Reads input files and QC/environment coordinates.
   - Computes centers of mass (COM) of QC and optional QC subregion.
   - Shifts QC coordinates relative to subregion COM.

2. Electric Field and Potential
   - Computes electric field and potential at the QC COM due to environment charges.
   - Corrects fields by subtracting reference QC values.

3. Trajectory Processing
   - Reads XTC frames and computes QC COM for each frame.
   - Aligns QC structure to reference geometry via Euler angle fitting.
   - Rotates QM dipoles into the MD frame.

4. Hamiltonian Construction
   - Builds Hamiltonian including diagonal energies and QC/environment interactions.
   - Diagonalizes Hamiltonian using LAPACK DSYEV to obtain perturbed eigenvalues
     and eigenvectors.

5. Dipole Evaluation
   - Projects the dipole operator in the basis of the perturbed eigenvectors.
   - Evaluates perturbed transition dipoles or perturbed state dipoles.

6. Output
   - Writes perturbed energies and dipole moments to output_Si.xvg:
     state dipoles for "Perturbed Electronic State Properties" calculations,
     transition dipoles for "Perturbed Properties for Spectral Absorption Bands".
   - Writes dipole fluctuation ratios to dipole_fluctuations_Si.xvg:
     state dipole ratios for "Perturbed Electronic State Properties" calculations,
     transition dipole ratios for "Perturbed Properties for Spectral Absorption Bands".
   - Writes the expansion coefficients of the unperturbed electronic states
     contributing to the i-th perturbed state to
     unpert_state_coefficients_to_Si_pert.xvg.

------------------------------------------------------------
Notes
------------------------------------------------------------
- Units: distances in Å and converted in Bohr, energies in eV or Hartree.
- QC = Quantum Center, environment = all atoms outside QC.
- Euler angle fitting ensures optimal alignment of QC reference and MD frame QC atoms.
- Ensure all input files exist in the working directory or provide full paths.
- The program requires valid input files generated by
  pmm-qmmm_input_generator.py (highly recommended) or an equivalent workflow.
- The number of electronic states includes the ground state.

------------------------------------------------------------
Authorship and Credits
------------------------------------------------------------
This code is based on a previous PMM implementation originally developed by
Prof. Massimiliano Aschi.

The present version has been modified, extended, and maintained by
Gianluca Dell'Orletta, under the supervision of Prof. Massimiliano Aschi and Prof. Isabella Daidone.

------------------------------------------------------------
Citation
------------------------------------------------------------
If you use PMM-QM/MM in published work, please cite:

Please also cite the original PMM methodology paper:

Aschi, M.; Spezia, R.; Di Nola, A.; Amadei, A. 
A first-principles method to model perturbed electronic wavefunctions: the effect of an external homogeneous electric field. 
Chem. Phys. Lett. 2001, 344, 374–380
doi.org/10.1016/S0009-2614(01)00638-8
