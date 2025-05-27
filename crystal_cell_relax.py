import os
import sys 
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.optimize import fire2 # Changed from BFGS to Fire2
from ase.constraints import ExpCellFilter 
import numpy as np 

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
fmax_threshold = 0.01 # eV/Angstrom 
max_opt_steps = 2000  # Maximum number of optimization steps (as requested)

# --- Check for Quantum ESPRESSO Executable ---
espresso_command = os.environ.get('ASE_ESPRESSO_COMMAND')
if not espresso_command:
    print("CRITICAL ERROR: 'ASE_ESPRESSO_COMMAND' environment variable is not set.")
    print("Please set it to the command that executes Quantum ESPRESSO's pw.x,")
    print("e.g., 'export ASE_ESPRESSO_COMMAND=\"pw.x\"' or 'export ASE_ESPRESSO_COMMAND=\"mpirun -np 4 pw.x\"'.")
    sys.exit(1) 

if 'pw.x' not in espresso_command:
    print("WARNING: 'ASE_ESPRESSO_COMMAND' does not contain 'pw.x'. Make sure it's correctly configured.")

# --- SET OMP_NUM_THREADS to 1 ---
os.environ['OMP_NUM_THREADS'] = '1' 
print(f"Setting OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']} for Quantum ESPRESSO calculations.")


print(f"Starting Optimization using ASE's Fire2 algorithm.")
print(f"  Quantum ESPRESSO Command: {espresso_command}")
print(f"  Fixed DFT Parameters: ecutwfc = {fixed_ecutwfc_ev} eV, kpts = {fixed_kpts}")
print(f"  ASE Fire2 Optimization Criteria: fmax = {fmax_threshold} eV/Angstrom, max_steps = {max_opt_steps}")
print("="*80)

optimization_results = {} 

for crystal in crystals:
    output_dir = f'{crystal}_ase_fire2_calculations' # Changed output directory name
    os.makedirs(output_dir, exist_ok=True) 

    try:
        initial_atoms = read(f'{crystal}.cif')
        initial_a = initial_atoms.get_cell()[0, 0] 
        print(f"\nProcessing {crystal}. Initial structure loaded from {crystal}.cif")
        print(f"  Initial Lattice Constant: {initial_a:.4f} Å")
    except FileNotFoundError:
        print(f"ERROR: {crystal}.cif not found. Skipping {crystal}.")
        optimization_results[crystal] = {'status': 'CIF not found'}
        continue 
    except Exception as e:
        print(f"ERROR: Could not read {crystal}.cif - {e}. Skipping {crystal}.")
        optimization_results[crystal] = {'status': f'Error reading CIF: {e}'}
        continue 

    calc_pseudopotentials = {}
    for symbol in initial_atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            print(f"ERROR: Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary for {crystal}. Skipping.")
            optimization_results[crystal] = {'status': f'Missing pseudopotential for {symbol}'}
            continue 
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    input_settings = {
        'control': {
            'calculation': 'scf', 
            'prefix': crystal,
            'outdir': output_dir,
            # Removed tprnfor and tstress as per request
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
            'conv_thr': 1.0e-7, # Tightened SCF convergence slightly, as requested
            'mixing_beta': 0.7,
        },
    }

    calc = Espresso(command=espresso_command, 
                    pseudopotentials=calc_pseudopotentials,
                    input_data=input_settings,
                    kpts=fixed_kpts)
    
    initial_atoms.set_calculator(calc) 

    final_energy = "Failed"
    relaxed_lattice_constant = "Failed"
    status = "Failed"

    try:
        f = ExpCellFilter(initial_atoms) 
        
        log_filename = f'{crystal}_ase_opt.log' # Log file name
        traj_filename = f'{crystal}_ase_opt.traj' # Trajectory file name

        opt = Fire2(f, trajectory=traj_filename, logfile=log_filename) # Changed to Fire2

        opt.run(fmax=fmax_threshold, steps=max_opt_steps) # Added steps limit

        final_energy = initial_atoms.get_potential_energy() 
        relaxed_cell = initial_atoms.get_cell()
        relaxed_lattice_constant = relaxed_cell[0, 0] 

        relaxed_filename = f'{crystal}_relaxed_ase_fire2.cif' # Changed output CIF name
        write(relaxed_filename, initial_atoms)

        status = "Converged"
        print(f"  {crystal} ASE Fire2 Optimization Completed Successfully.") # Updated message
        print(f"    Final Energy: {final_energy:.4f} eV")
        print(f"    Relaxed Lattice Constant: {relaxed_lattice_constant:.4f} Å")
        print(f"    Relaxed structure saved to {relaxed_filename}")
        print(f"    Optimization trajectory saved to {traj_filename}")
        print(f"    Optimization log saved to {log_filename}")

    except Exception as e:
        status = f"Failed ({e})"
        print(f"  {crystal} ASE Fire2 Optimization Failed. Error: {e}") # Updated message
    
    optimization_results[crystal] = {
        'initial_a': initial_a,
        'final_energy': final_energy,
        'relaxed_lattice_constant': relaxed_lattice_constant,
        'status': status
    }

print("\n" + "="*80)
print("All ASE Fire2 Optimizations Attempted.") # Updated message
print("Final Summary of Results:")
for crystal, results in optimization_results.items():
    if results['status'] == "Converged":
        print(f"  {crystal}: Initial a = {results['initial_a']:.4f} Å, Final a = {results['relaxed_lattice_constant']:.4f} Å, Energy = {results['final_energy']:.4f} eV")
    else:
        print(f"  {crystal}: Calculation {results['status']}")
print("="*80)