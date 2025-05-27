import os
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.optimize import fire2 # Import ASE's fire2 optimizer
from ase.filters import ExpCellFilter # CORRECTED: Import from ase.constraints
import numpy as np # Still useful for np.linspace, if used elsewhere

# --- System and Pseudopotential Definitions ---
crystals = ['TiN', 'ZrN', 'NbN', 'ScN', 'VN']
pseudopotentials = {'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
                    'Zr': 'Zr_pbe_v1.uspp.F.UPF', 
                    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
                    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                    'V': 'v_pbe_v1.4.uspp.F.UPF',
                    'N': 'N.pbe-n-radius_5.UPF'}

# --- DFT Parameters for Optimization ---
fixed_ecutwfc_ev = 650  # eV
fixed_kpts = (7, 7, 7)
eV_to_Ry = 1.0 / 13.6057
fixed_ecutwfc_ry = fixed_ecutwfc_ev * eV_to_Ry

# --- ASE Optimizer Settings ---
fmax_threshold = 0.01 # eV/Angstrom (as specified)

optimization_results = {} # Dictionary to store final results

for crystal in crystals:
    # Create a dedicated output directory for each crystal's QE calculations
    output_dir = f'{crystal}_ase_calculations'
    os.makedirs(output_dir, exist_ok=True) 

    # Read the initial crystal structure
    initial_atoms = read(f'{crystal}.cif')
    initial_a = initial_atoms.get_cell()[0, 0] # Initial lattice constant

    # Prepare the pseudopotentials dictionary for the current crystal
    calc_pseudopotentials = {}
    for symbol in initial_atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            raise ValueError(f"Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary.")
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    # --- Quantum ESPRESSO Input Settings (now just for SCF calculations) ---
    input_settings = {
        'control': {
            'calculation': 'scf', # Changed to scf; ASE BFGS controls relaxation loop
            'prefix': crystal,
            'outdir': output_dir,
            'tprnfor': True,           # Print forces in output
            'tstress': True,           # Print stress in output
        },
        'system': {
            'ecutwfc': fixed_ecutwfc_ry,
            'ecutrho': fixed_ecutwfc_ry * 8,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.01,
            'input_dft': 'pbe',
            'vdw_corr': 'dft-d3',
            'dftd3_version': 4,
        },
        'electrons': {
            'conv_thr': 1.0e-7,
            'mixing_beta': 0.7,
        },
        # Removed: 'ions' and 'cell' sections (as per request)
    }

    # Initialize the Espresso calculator
    calc = Espresso(pseudopotentials=calc_pseudopotentials,
                    input_data=input_settings,
                    kpts=fixed_kpts)
    
    # Assign the calculator to the atoms object
    initial_atoms.calc= calc

    final_energy = "Failed"
    relaxed_lattice_constant = "Failed"
    status = "Failed"

    try:
        # --- Set up ASE's  Optimizer ---
        # ExpCellFilter wraps the atoms object to allow cell optimization
        # The optimizer will operate on 'f' which represents the atoms object and its cell
        f = ExpCellFilter(initial_atoms) 
        
        # Define trajectory and log file names
        log_filename = f'{crystal}_ase_opt.log'
        traj_filename = f'{crystal}_ase_opt.traj'

        # Initialize the  optimizer
        # The 'f' argument is the atoms object (or filter) to be optimized

        opt = fire2(f, trajectory=traj_filename, logfile=log_filename)

        # It will call the calculator (Espresso) at each step.
        opt.run(fmax=fmax_threshold, steps=2000)

        # After optimization, the initial_atoms object is updated with the relaxed structure
        final_energy = initial_atoms.get_potential_energy() # Get final energy from relaxed structure
        relaxed_cell = initial_atoms.get_cell()
        relaxed_lattice_constant = relaxed_cell[0, 0] # Assuming cubic cell

        # Save the final optimized structure to a new CIF file
        relaxed_filename = f'{crystal}_relaxed_ase_bfgs.cif'
        write(relaxed_filename, initial_atoms)

        status = "Converged"
        print(f"  {crystal} ASE BFGS Optimization Completed Successfully.")
        print(f"    Final Energy: {final_energy:.4f} eV")
        print(f"    Relaxed Lattice Constant: {relaxed_lattice_constant:.4f} Å")
        print(f"    Relaxed structure saved to {relaxed_filename}")
        print(f"    Optimization trajectory saved to {traj_filename}")
        print(f"    Optimization log saved to {log_filename}")

    except Exception as e:
        status = f"Failed ({e})"
        print(f"  {crystal} ASE BFGS Optimization Failed. Error: {e}")
