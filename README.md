# CO₂ER Database

## Contents
* [Overview](#overview)
* [COSMO-RS](#cosmo-rs)
  * [Pure_solvent_properties_data.csv](#pure_solvent_properties_data.csv)
  * [Solubility_co2.csv](#solubility_co2csv)
  * [Cosmo_files](#cosmo_files)
* [Classical Molecular Dynamics](#classical-molecular-dynamics)
  * [Coordination_Numbers (CN)](#coordination_numbers-cn)
  * [Diffusion](#diffusion)
  * [Force_Fields](#force_fields)
  * [MD_Data_files](#md_data_files)
  * [MD_templates](#md_templates)
  * [Radial_Distribution_Functions (RDF)](#radial_distribution_functions-rdf)
  * [Solvent_data.csv](#solvent_datacsv)
* [Machine Learning](#machine-learning)
* [DFT Computed Properties](#dft-computed-properties)
  * [Binding_Energy](#binding_energy)
  * [Electron_Affinity](#electron_affinity)
  * [Ionization_Potentials](#ionization_potentials)
  * [DFT_Computed_Properties.csv](#dft_computed_propertiescsv)
* [Scripts](#scripts)
* [XYZ Files](#xyz-files)
* [Updates](#updates)
* [Licensing](#licensing)
* [Contact](#contact)

---

## Overview

The CO₂ER database is a publicly available repository of computed molecular properties and simulation-ready inputs for the study of organic solvents in electrochemical CO₂ reduction (CO₂ER). The data spans quantum chemical (DFT), thermodynamic (COSMO-RS), and classical molecular dynamics (MD) calculations and covers solvents from 11 chemical families.

This resource includes DFT-optimized structures, binding energies, electronic descriptors (HOMO, LUMO, IP, EA), MD-derived solvation structure and diffusion data, and COSMO-RS-predicted properties such as CO₂ solubility and viscosity. It is intended to support high-throughput screening, simulation, and machine learning model development for CO₂ER electrolyte design.

---

## COSMO-RS

### Pure_solvent_properties_data.csv
Contains 24 solvent properties predicted using COSMO-RS, including melting point, flash point, viscosity, van der Waals volume, polarity, and more. These properties are useful for assessing solvent stability and CO₂ compatibility.

### Solubility_co2.csv
Estimated solubility of CO₂ in pure organic solvents, as predicted by COSMO-RS. Data is critical for evaluating the absorption potential of solvents in CO₂ER.

### Cosmo_files
DFT-derived COSMO files used as input for COSMO-RS simulations. Prepared using Gaussian via MISPR workflows.

---

## Classical Molecular Dynamics

### Coordination_Numbers (CN)
Raw CSV files with coordination numbers between CO₂ (carbon atom) and surrounding electrolyte components including solvent, salt cation (TBA⁺), and anion (BF₄⁻).

### Diffusion
CSV files reporting the diffusion coefficients of each species (TBA⁺, BF₄⁻, solvent, CO₂) in electrolyte systems with 0.1 M [TBA⁺][BF₄⁻] and 0.1 M CO₂. Includes standard deviations and R² values.

### Force_Fields
OPLS-AA-compatible JSON force field files for each solvent. These are required for running classical MD simulations.

### MD_Data_files
LAMMPS-compatible `.data` files describing initial atomic positions, atom types, and topology for electrolyte systems.

### MD_templates
Reusable LAMMPS input scripts for various simulation stages: energy minimization, NPT/NVT equilibration, melting, quenching, and production.

### Radial_Distribution_Functions (RDF)
CSV files containing RDF data between CO₂ atoms (C and O) and the center of mass of TBA⁺, BF₄⁻, and solvent molecules. Useful for analyzing solvation shell structure.

### Solvent_data.csv
Metadata file summarizing MD setup details: number of molecules, atom types, initial seeds, molecular weights, and other relevant system parameters.

---

## Machine Learning

This section (coming soon) will include pre-trained ML models for predicting solvent properties relevant to CO₂ER: CO₂ solubility, viscosity, ionization potential, and electron affinity. 

Each model will be trained on computed datasets using four classes of descriptors:
- RDKit molecular descriptors
- Mordred chemical descriptors
- DFT-based descriptors (e.g., HOMO, LUMO, charge)
- COSMO-RS-derived thermodynamic features

---

## DFT Computed Properties

### Binding_Energy
DFT-calculated binding energies for solvent–CO₂ complexes. Stored in raw JSON format, including all metadata and convergence details.

### Electron_Affinity
DFT-computed electron affinity values for each solvent. Includes frontier orbital data and solvent geometries.

### Ionization_Potentials
DFT-computed ionization potentials for each solvent. Includes both vertical and adiabatic IP values, where applicable.

### DFT_Computed_Properties.csv
Consolidated table of 13 molecular properties from DFT, including HOMO/LUMO energies, SCF energy, dipole moment, polarizability, and Mulliken charges.

---

## Scripts

- `run_md.py`: Automates LAMMPS-based MD simulations using MISPR for systems with [TBA⁺][BF₄⁻] and CO₂ in various solvents.
- `run_ip_ea.py`: Automates IP and EA calculations using Gaussian via MISPR.
- `run_be.py`: Automates binding energy calculations between solvents and CO₂.

---

## XYZ Files

XYZ format geometry files of DFT-optimized solvent molecules. Useful for visualization or preparing inputs for other simulation workflows.

---

## Updates
See [`updates.md`](./updates.md) for changelog and dataset version history.

---

## Licensing
This repository is released under the MIT License and made publicly available under the Creative Commons Attribution 4.0 (CC BY 4.0) license. You may copy, distribute, and adapt the content with appropriate credit.

---

## Contact
For questions or suggestions, please contact:

**Kuldeepsinh Raj**  
Ph.D. Student, Materials Science and Engineering  
Stony Brook University  
[kuldeepsinh.raj@stonybrook.edu](mailto:kuldeepsinh.raj@stonybrook.edu)
