# Script to submit IP and EA Gaussian workflows using MISPR

import os
from os.path import isfile, join
from pathlib import Path
from fireworks import LaunchPad
from mispr.gaussian.workflows.base.ip_ea import get_ip_ea

working_dir = os.getcwd()

lpad = LaunchPad.auto_load()

func = "wB97XD"
basis = "Def2SVP"


ip_dir = os.path.abspath(os.path.join(working_dir, "..", "ip"))
Path(ip_dir).mkdir(parents=True, exist_ok=True)

ea_dir = os.path.abspath(os.path.join(working_dir, "..", "ea"))
Path(ea_dir).mkdir(parents=True, exist_ok=True)

xyz_dir = os.path.abspath(os.path.join(working_dir, "..", "xyz")) # Directory containing input XYZ files
molecule_files = [f for f in os.listdir(xyz_dir) if isfile(join(xyz_dir, f))]

for molecule in molecule_files:
    mol_name = molecule.split(".")[0]

    ip_mol_dir = os.path.join(ip_dir, mol_name)
    Path(ip_mol_dir).mkdir(parents=True, exist_ok=True)

    ea_mol_dir = os.path.join(ea_dir, mol_name)
    Path(ea_mol_dir).mkdir(parents=True, exist_ok=True)

    charge = 0
    spin = 1

    # Gaussian DFT settings 
    opt_gaussian_inputs = {
        "functional": func,
        "basis_set": basis,
        "route_parameters": {
            "Opt": "(calcfc, tight)",
            "SCF": "(tight, xqc)",
            "int": "ultrafine",
            "NoSymmetry": None,
            "test": None,
        },
        "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "28"},
    }

    freq_gaussian_inputs = {
        "functional": func,
        "basis_set": basis,
        "route_parameters": {
            "Freq": None,
            "iop(7/33=1)": None,
            "NoSymmetry": None,
            "test": None,
        },
        "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "28"},
    }
   
   # IP calculation
    ip_wf = get_ip_ea(
        mol_operation_type="get_from_file",
        mol=os.path.join(xyz_dir, molecule),
        ref_charge=charge,
        single_step=True,
        num_electrons=1,
        solvent_gaussian_inputs=f"(SMD, Solvent=Acetone)",
        states=["cation"],
	    phases=["solution"],
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        name="ip_calculation",
        working_dir=ip_mol_dir,
        save_to_file=True,
    )
    lpad.add_wf(ip_wf)
    
    # EA calculation
    ea_wf = get_ip_ea(
        mol_operation_type="get_from_file",
        mol=os.path.join(xyz_dir, molecule),
        ref_charge=charge,
        single_step=True,
        num_electrons=1,
        solvent_gaussian_inputs=f"(SMD, Solvent=Acetone)",
        states=["anion"],
	    phases=["solution"],
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        name="ea_calculation",
        working_dir=ea_mol_dir,
        save_to_file=True,
    )
    lpad.add_wf(ea_wf)
