# Script to calculate diffusion coefficient from a LAMMPS log file using MDPropTools.

from mdproptools.dynamical.diffusion import Diffusion

log_pattern = "log.lammps" # Name of the LAMMPS log file containing MSD data

DIFF = Diffusion()
Msd = DIFF.get_msd_from_log(log_pattern) # Extract MSD from the LAMMPS log file

Diff = DIFF.calc_diff(Msd, save=True, plot=True) # Calculate diffusion coefficient, Results are saved to file and plots are generated