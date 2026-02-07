# Coordination Number (CN) calculation using MISPR analysis workflow

import os
import json
import pandas as pd
from fireworks import LaunchPad, Workflow
from mispr.lammps.workflows.base import lammps_analysis_fws

# Load FireWorks launchpad
lpad = LaunchPad.auto_load()

# Load solvent list and pair-specific CN cutoffs (r_cut) from CSV
working_dir = os.getcwd()
solv_df = pd.read_csv(os.path.join(working_dir, "..", "solvent_data.csv"))
r_cut_df = pd.read_csv(os.path.join(working_dir, "..", "r_cut.csv"))
r_cut_dict = r_cut_df.set_index("solvent").to_dict(orient="index")
solv_output_files = solv_df["name"].tolist()

incomplete_systems = []

for solv_name in solv_output_files:
    print(f"Processing: {solv_name}")

    # Get pair-specific r_cut values for this solvent system
    try:
        rcut_row = r_cut_dict[solv_name]
        r_cut = [
            float(rcut_row["C_CO2 - BF4"]),
            float(rcut_row["C_CO2 - Solvent"]),
            float(rcut_row["C_CO2 - TBA+"]),
            float(rcut_row["O_CO2 - BF4"]),
            float(rcut_row["O_CO2 - Solvent"]),
            float(rcut_row["O_CO2 - TBA+"]),
        ]
    except KeyError:
        print(f"[Warning] r_cut not found for {solv_name} — skipping")
        incomplete_systems.append(solv_name)
        continue

    # Paths
    nvt_dir = os.path.join(working_dir, "..", "outputs", solv_name, "nvt")
    output_dir = os.path.join(working_dir, "..", "outputs", solv_name)

    if "restart.nvt_10-ns.restart" not in os.listdir(nvt_dir):
        incomplete_systems.append(solv_name)
        continue
    
    # Read FW.json
    with open(os.path.join(nvt_dir, "FW.json"), "r") as f:
        fw_json = json.load(f)

    num_mols_list = fw_json["spec"]["num_mols_list"]
    num_atoms_per_mol = fw_json["spec"]["num_atoms_per_mol"]
    mass_list = fw_json["spec"]["default_masses"]
    num_types = len(mass_list)

    # Custom molecule ordering: TBA → Solvent → BF₄⁻ → CO₂
    solv_atom_types = int(solv_df.loc[solv_df["name"] == solv_name, "num_types"].values[0])

    co2_o_atom_type = 6 + solv_atom_types
    co2_c_atom_type = 6 + solv_atom_types + 1
    cation_mol_type = 1
    solvent_mol_type = 2
    anion_mol_type = 3

    partial_relations_mol = [
        [co2_c_atom_type, co2_c_atom_type, co2_c_atom_type,
         co2_o_atom_type, co2_o_atom_type, co2_o_atom_type],
        [anion_mol_type, solvent_mol_type, cation_mol_type,
         anion_mol_type, solvent_mol_type, cation_mol_type]
    ]

    nnodes = 1
    ntasks_node = 1
    ntasks = int(nnodes * ntasks_node)

    qadapter = {
        "queue": "rajput-28core",
        "walltime": "48:00:00",
        "job_name": f"cn_{solv_name}",
        "nodes": nnodes,
        "ntasks": ntasks,
        "ntasks_per_node": ntasks_node,
        "mail_user": "your_email_id",  # replace with your email
        "mail_type": "ALL",
        "pre_rocket": "module load /gpfs/projects/RajputGroup/modulefiles/lammps/gcc12.1/mvapich2/23Jun2022"
    }

    analysis_list = ["cn"]
    analysis_settings = [{
        "cn_type": "molecular",
        "r_cut": r_cut,  
        "bin_size": 0.5, 
        "mass": mass_list,
        "partial_relations": partial_relations_mol,
        "filename": os.path.join(nvt_dir, "dump.nvt.*.dump"),
        "num_mols": num_mols_list,
        "num_atoms_per_mol": num_atoms_per_mol,
        "path_or_buff": "cn_co2_com_molecular.csv"
    }]

    fws, links_dict = lammps_analysis_fws(analysis_list, analysis_settings, output_dir, qadapter=qadapter)
    wf = Workflow(fws, links_dict, name=f"cn_{solv_name}")
    lpad.add_wf(wf)

print("Incomplete systems:", incomplete_systems)

