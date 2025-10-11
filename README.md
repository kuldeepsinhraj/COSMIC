# COSMIC(CO₂ Solvent Materials Informatics Collection) Database
<img width="875" height="473" alt="cosmic_logo" src="https://github.com/user-attachments/assets/b8161d00-06da-4626-9f0d-c6aee2e23009" />


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
  * [Data](#data)
    * [Chemical Class](#chemical-class)
    * [Data_computed_properties](#data_computed_properties)
    * [Descriptors](#descriptors)
    * [df_filtered.csv](#df_filteredcsv)
    * [df_filtered_ip.csv](#df_filtered_ipcsv)
    * [filtered_ea_clean.csv](#filtered_ea_cleancsv)
  * [src](#src)
    * [catboost_model_ea.cbm](#catboost_model_eacbm)
    * [catboost_model_ea_dft.cbm](#catboost_model_ea_dftcbm)
    * [cgb_ea_model.json](#cgb_ea_modeljson)
    * [data_merge.ipynb](#data_mergeipynb)
    * [model_catboost_ea.ipynb](#model_catboost_eaipynb)
    * [model_rf_ip.ipynb](#model_rf_ipipynb)
    * [model_viscosity_xgb.ipynb](#model_viscosity_xgbipynb)
    * [model_xgboost_solubility.ipynb](#model_xgboost_solubilityipynb)
    * [preprocess.ipynb](#preprocessipynb)
    * [preprocess_ea.ipynb](#preprocess_eaipynb)
    * [preprocess_ip.ipynb](#preprocess_ipipynb)
    * [xgb_ea_model.json](#xgb_ea_modeljson)
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

The CO₂ER database is a publicly available repository of computed molecular properties and simulation-ready inputs for the study of organic solvents in electrochemical CO₂ reduction (CO₂ER). The data spans quantum chemical (DFT), thermodynamic (COSMO-RS), and classical molecular dynamics (MD) calculations and covers solvents from 11 chemical classes.

The dataset includes DFT-optimized geometries (in XYZ format) and a range of DFT-computed properties such as binding energies, partial atomic charges, polarizabilities, ionization potentials, electron affinities, and HOMO–LUMO energies. It also contains MD-derived properties like radial distribution functions, coordination numbers, and diffusion coefficients, as well as COSMO-RS-predicted properties such as CO₂ solubility, viscosity, and more. This resource is designed to support high-throughput screening, molecular simulation, and machine learning model development for CO₂ER electrolyte discovery.

---

## COSMO-RS

### Pure_solvent_properties_data.csv
Contains 24 solvent properties predicted using COSMO-RS, including melting point, flash point, viscosity, van der Waals volume, polarity, and more.

### Solubility_co2.csv
Estimated solubility of CO₂ in pure organic solvents, as predicted by COSMO-RS.

### Cosmo_files
DFT-derived COSMO files used as input for COSMO-RS simulations. Prepared using DFT Gaussian calculation via MISPR workflows.

---

## Classical Molecular Dynamics

### Coordination_Numbers (CN)
Raw CSV files with coordination numbers between CO₂ (carbon atom) and surrounding electrolyte solvent.

### Diffusion
CSV files reporting the diffusion coefficients of each species (TBA⁺, BF₄⁻, solvent, CO₂) in electrolyte systems with 0.1 M [TBA⁺][BF₄⁻] and 0.1 M CO₂. Includes standard deviations and R² values.

### Force_Fields
json files containing the solvent force field parameters (OPLS/AA) used to run the MD simulations. 

### MD_Data_files
LAMMPS-compatible data files describing initial atomic positions, atom types, and topology for electrolyte systems.

### MD_templates
LAMMPS template files used to run the following MD steps: energy minimization, NPT equilibration, melting and quenching, and NVT production.

### Radial_Distribution_Functions (RDF)
CSV files containing RDF data between CO₂ atoms (C and O) and the center of mass of TBA⁺, BF₄⁻, and solvent molecules. Useful for analyzing solvation shell structure.

### Solvent_data.csv
Metadata file summarizing MD setup details: number of molecules, atom types, initial seeds, molecular weights, and other relevant system parameters.

---

## Machine Learning

This section includes pre-trained ML models for predicting solvent properties relevant to CO₂ER: CO₂ solubility, viscosity, ionization potential, and electron affinity. 

### Data

#### Chemical Class
Folder containing `chemical_class.csv`, which has a detailed description of the chemical class for each molecule.

#### Data_computed_properties
Folder containing `data_computed_properties.csv` with computed properties used in ML workflows.

#### Descriptors
- RDKit molecular descriptors  
- Mordred chemical descriptors  
- `dft_descriptors.csv`: DFT-based descriptors  
- `cosmors_descriptors.csv`: COSMO-RS-derived descriptors

#### df_filtered.csv
Filtered data containing ~400–500 columns from the main `df_merged.csv`.

#### df_filtered_ip.csv
Filtered data for ionization potential containing ~100 columns from `df_filtered.csv`.

#### filtered_ea_clean.csv
Filtered data for electron affinity containing ~100 columns from `df_filtered.csv`.

### src

#### catboost_model_ea.cbm
Saved EA model.

#### catboost_model_ea_dft.cbm
Saved EA model with DFT descriptors.

#### cgb_ea_model.json
Saved EA model in JSON format.

#### data_merge.ipynb
Script to merge data of different molecules and their respective descriptors.

#### model_catboost_ea.ipynb
Script for EA model using the CatBoost algorithm.

#### model_rf_ip.ipynb
Script for IP model using the Random Forest algorithm.

#### model_viscosity_xgb.ipynb
Script for Viscosity model using the XGBoost algorithm.

#### model_xgboost_solubility.ipynb
Script for Solubility model using the XGBoost algorithm.

#### preprocess.ipynb
Script for initial preprocessing to create `df_filtered.csv`.

#### preprocess_ea.ipynb
Script for further preprocessing of EA data to create `filtered_ea_clean.csv`.

#### preprocess_ip.ipynb
Script for further preprocessing of IP data to create `df_filtered_ip.csv`.

#### xgb_ea_model.json
EA model serialized for XGBoost.



---

## DFT Computed Properties

### Binding_Energy
DFT-calculated binding energies for solvent–CO₂ complexes. Stored in raw JSON format, including all metadata.

### Electron_Affinity
Raw JSON files containing DFT calculations of the solvent electron affinity, including all associated metadata.

### Ionization_Potentials
Raw JSON files containing DFT calculations of the solvent ionization potentials, including all associated metadata.

### DFT_Computed_Properties.csv
Consolidated table of 13 molecular properties from DFT, including HOMO/LUMO energies, SCF energy, dipole moment, polarizability, charges and more.

---

## Scripts

- `run_md.py`: script for running the automated MD simulations of the electrolyte systems composed of 0.1 M [TBA⁺][BF₄⁻] salt and 0.1 M CO₂ in various solvents system using MISPR
- `run_ip_ea.py`: script for running the automated IP and EA simulations of the solvents system using MISPR.
- `run_be.py`: script for running the automated binding energy simulations solvents system using MISPR.
---

## XYZ Files

XYZ format geometry files of DFT-optimized solvent molecules. 

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
