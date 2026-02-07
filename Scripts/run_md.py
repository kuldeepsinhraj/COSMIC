"""
Script for running bulk electrolyte MD simulations using the MISPR LAMMPS workflow. 

Note: Paths are relative to the working directory where this script is launched.
"""

import os
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from fireworks import LaunchPad, Workflow
from mispr.lammps.workflows.base import lammps_workflow
from pymatgen.core.structure import Molecule

lpad = LaunchPad.auto_load()  # FireWorks LaunchPad

working_dir = os.getcwd()

solv_pdb_path = os.path.abspath(os.path.join("..", "pdb_solvents"))
solute_pdb_path = os.path.abspath(os.path.join("..", "pdb_solutes"))

solv_pdb_files = os.listdir(solv_pdb_path)

solv_df = pd.read_csv(os.path.join(working_dir, "..", "solvent_data.csv")) # Load solvent metadata

 # Gather solvent PDB files
for file in solv_pdb_files:
    solv_name = file.split(".")[0].strip()
    
    # Pull system-specific parameters
    seed = int(solv_df.loc[solv_df["name"] == solv_name, ["position_seed"]].to_numpy()[0][0])
    solv_num = int(solv_df.loc[solv_df["name"] == solv_name, ["num_mols"]].to_numpy()[0][0])
    charge_scaling = 1 / solv_df.loc[solv_df["name"] == solv_name, ["n"]].to_numpy()[0][0]
    solv_atom_types = int(solv_df.loc[solv_df["name"] == solv_name, ["num_types"]].to_numpy()[0][0])
    
    #Build LAMMPS group definitions based on atom type layout
    GROUPS = "group TBA type 1 2 3\ngroup SOL type " + " ".join([str(4 + i) for i in range(solv_atom_types)]) + "\ngroup BF type " + " ".join([str(4 + solv_atom_types + i) for i in range(2)]) + "\ngroup CO type " + " ".join([str(6 + solv_atom_types + i) for i in range(2)])

    output_dir = os.path.join(working_dir, "..", "outputs", solv_name)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # STATIC
    num_mols = {"tba": 20,
                "bf4": 20,
                "solv": solv_num,
                "co2": 20}

    box_side = 69.3
    origin = [0.25, 0.25, 0.25]

    SYS_TEMP = 298.15 # K
    SYS_PRESS = 1 # atm

    # define cation
    tba = Molecule.from_file(os.path.join(solute_pdb_path, "tba.pdb"))
    tba.set_charge_and_spin(1, 1)

    # define anion
    bf4 = Molecule.from_file(os.path.join(solute_pdb_path, "bf4.pdb"))
    bf4.set_charge_and_spin(-1, 1)

    # define solvent
    solv = Molecule.from_file(os.path.join(solv_pdb_path, file))
    solv.set_charge_and_spin(0, 1)

    # define CO2
    co2 = Molecule.from_file(os.path.join(solute_pdb_path, "co2.pdb"))
    co2.set_charge_and_spin(0, 1)

    system_species_data = {
        "solv": {
            "molecule": solv,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_opls",
            "ff_param_data": os.path.join(solv_pdb_path, file), #file_name_or_path
            "mol_mixture_type": "Solvents",
            "mixture_data": num_mols["solv"],
            "save_ff_to_file": True,
        },
        "bf4": {
            "molecule": bf4,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_opls",
            "ff_param_data": os.path.join(solute_pdb_path, "bf4.pdb"),
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["bf4"],
            "save_ff_to_file": True,
        },
        "tba": {
            "molecule": tba,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_opls",
            "ff_param_data": os.path.join(solute_pdb_path, "tba.pdb"),
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["tba"],
            "save_ff_to_file": True,
        },
        "co2": {
            "molecule": co2,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_opls",
            "ff_param_data": os.path.join(solute_pdb_path, "co2.pdb"),
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["co2"],
            "save_ff_to_file": True,
        },
    }

    # LAMMPS recipe
    recipe = [["emin", ["template_filename", "emin_gaff"]],
            ["npt", ["template_filename", "npt"]],
            ["melt", ["template_filename", "nvt"]],
            ["quench", ["template_filename", "nvt"]],
            ["nvt", ["template_filename", "nvt"]],]

    # LAMMPS settings
    SPEC_BOND_STYLE = "lj/coul 0.0 0.0 0.5"
    DIHED_STYLE = "fourier"
    PAIR_MOD = "geometric"
    COMPUTES = "compute TBA_msd TBA msd com yes \ncompute SOL_msd SOL msd com yes \ncompute BF_msd BF msd com yes \ncompute CO_msd CO msd com yes"
    THERMO_COMP = "c_TBA_msd[4] c_SOL_msd[4] c_BF_msd[4] c_CO_msd[4]"

    recipe_settings = [{"special_bonds_style": SPEC_BOND_STYLE,
                        "dihedral_style": DIHED_STYLE,
                        "pair_modify_value": PAIR_MOD,
                        "data_filename": "../data.mixture",
                        "restart_finalname": "restart.emin.restart"},
                        {"special_bonds_style": SPEC_BOND_STYLE,
                        "dihedral_style": DIHED_STYLE,
                        "pair_modify_value": PAIR_MOD,
                        "temperature_initial": SYS_TEMP,
                        "temperature_final": SYS_TEMP,
                        "pressure_initial": SYS_PRESS,
                        "pressure_final": SYS_PRESS,
                        "restart_filename": "../emin/restart.emin.restart",
                        "restart_final_filename": "restart.npt.restart",
                        "dump_modify_logic": "",
                        "run": 2000000},
                        {"special_bonds_style": SPEC_BOND_STYLE,
                        "dihedral_style": DIHED_STYLE,
                        "pair_modify_value": PAIR_MOD,
                        "restart_filename": "../npt/restart.npt.restart",
                        "temperature_initial": 500.0,
                        "temperature_final": 500.0,
                        "dump_modify_logic": "",
                        "run": 1000000,
                        "restart_final_filename": "restart.melt_final.restart"},
                        {"special_bonds_style": SPEC_BOND_STYLE,
                        "dihedral_style": DIHED_STYLE,
                        "pair_modify_value": PAIR_MOD,
                        "restart_filename": "../melt/restart.melt_final.restart",
                        "temperature_initial": 500.0,
                        "temperature_final": SYS_TEMP,
                        "dump_modify_logic": "",
                        "run": 2000000,
                        "restart_final_filename": "restart.quench_final.restart"},
                        {"special_bonds_style": SPEC_BOND_STYLE,
                        "dihedral_style": DIHED_STYLE,
                        "pair_modify_value": PAIR_MOD,
                        "restart_filename": "../quench/restart.quench_final.restart",
                        "group_definitions": GROUPS,
                        "compute_definitions": COMPUTES,
                        "thermo_style_compute": THERMO_COMP,
                        "dump_modify_logic": "",
                        "reset_timestep_logic": "",
                        "temperature_initial": SYS_TEMP,
                        "temperature_final": SYS_TEMP,
                        "run": 10000000,
                        "restart_final_filename": "restart.nvt_10-ns.restart"},]

    # SLURM / queue adapter settings
    nnodes = 4 # nodes
    ntasks_node = 28 # ntasks_per_node
    ntasks = int(nnodes * ntasks_node) # ntasks
    mail_user = "your_email_id" #provide your email id
    mail_type = "ALL"
    partition = "rajput-28core"
    walltime = "48:00:00"
    pre_rocket = "module load /gpfs/projects/RajputGroup/modulefiles/lammps/gcc12.1/mvapich2/23Jun2022"

    qadapter = [{"queue": partition,
                    "walltime": walltime,
                    "job_name": "e",
                    "nodes": nnodes,
                    "ntasks": ntasks,
                    "ntasks_per_node": ntasks_node,
                    "mail_user": mail_user,
                    "mail_type": mail_type,
                    "pre_rocket": pre_rocket},
                {"queue": partition,
                    "walltime": walltime,
                    "job_name": "p",
                    "nodes": nnodes,
                    "ntasks": ntasks,
                    "ntasks_per_node": ntasks_node,
                    "mail_user": mail_user,
                    "mail_type": mail_type,
                    "pre_rocket": pre_rocket},
                {"queue": partition,
                    "walltime": walltime,
                    "job_name": "m",
                    "nodes": nnodes,
                    "ntasks": ntasks,
                    "ntasks_per_node": ntasks_node,
                    "mail_user": mail_user,
                    "mail_type": mail_type,
                    "pre_rocket": pre_rocket},
                {"queue": partition,
                    "walltime": walltime,
                    "job_name": "q",
                    "nodes": nnodes,
                    "ntasks": ntasks,
                    "ntasks_per_node": ntasks_node,
                    "mail_user": mail_user,
                    "mail_type": mail_type,
                    "pre_rocket": pre_rocket},
                {"queue": partition,
                    "walltime": walltime,
                    "job_name": "v10",
                    "nodes": nnodes,
                    "ntasks": ntasks,
                    "ntasks_per_node": ntasks_node,
                    "mail_user": mail_user,
                    "mail_type": mail_type,
                    "pre_rocket": pre_rocket},]

    # Create workflow
    wf = lammps_workflow(
        system_species_data=system_species_data,
        system_mixture_type="number of molecules",
        box_data=box_side,
        origin=origin,
        position_seed=seed,
        box_data_type="cubic",
        working_dir=output_dir,
        recipe=recipe,
        recipe_settings=recipe_settings,
        recipe_qadapter=qadapter,
        scale_charges=True,
        charge_scaling_factor=charge_scaling
    )

    lpad.add_wf(wf) #Submit to LaunchPad
