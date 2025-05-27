import matplotlib.pyplot as plt
from ase.io import read
from ase.calculators.espresso import Espresso
import numpy as np
import os
import csv

# --- System and Pseudopotential Definitions ---
crystals = ['TiN', 'ZrN', 'NbN', 'ScN', 'VN']
pseudopotentials = {'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
                    'Zr': 'Zr_pbe_v1.uspp.F.UPF', # Corrected capitalization for Zr UPF
                    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
                    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                    'V': 'v_pbe_v1.4.uspp.F.UPF',
                    'N': 'N.pbe-n-radius_5.UPF'}

# --- Fixed DFT Parameters (from prior convergence tests) ---
fixed_ecutwfc_ev = 600
fixed_kpts = (7, 7, 7) # Using (7,7,7) as the highest tested k-point mesh
eV_to_Ry = 1.0 / 13.6057
fixed_ecutwfc_ry = fixed_ecutwfc_ev * eV_to_Ry

# --- Lattice Constant Scan Parameters ---
# Generate scaling factors around the initial lattice constant
# E.g., from 96% to 104% of initial, with 7 steps
lattice_constant_scaling_factors = np.linspace(0.96, 1.04, 7).round(3) 

# Dictionary to store energies and actual lattice constants for CSV output
# Structure: {crystal_name: {actual_lattice_constant: energy, ...}}
# We need to map actual_lattice_constant to energy for each crystal
# It's better to store a list of (actual_lattice_constant, energy) pairs per crystal
crystal_scan_results = {crystal: [] for crystal in crystals}


print(f"Starting Lattice Constant Scan with fixed ecutwfc = {fixed_ecutwfc_ev} eV and kpts = {fixed_kpts}.")
print("Scanning lattice constant scaling factors: " + ", ".join(map(str, lattice_constant_scaling_factors)))
print("="*80)

for crystal in crystals:
    # Create an output directory for each crystal's Quantum ESPRESSO calculations
    output_dir = f'{crystal}_lattice_scan_calculations'
    os.makedirs(output_dir, exist_ok=True) 

    # Read the initial crystal structure
    initial_atoms = read(f'{crystal}.cif')
    
    # Get the initial lattice constant (assuming cubic or effectively cubic conventional cell)
    initial_a = initial_atoms.get_cell()[0, 0]
    print(f"\nProcessing {crystal}. Initial lattice constant (from CIF): {initial_a:.4f} Å")

    # Prepare the pseudopotentials dictionary for the current crystal
    calc_pseudopotentials = {}
    for symbol in initial_atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            raise ValueError(f"Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary.")
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    for scaling_factor in lattice_constant_scaling_factors:
        atoms_to_calc = initial_atoms.copy()
        
        scaled_a = initial_a * scaling_factor
        atoms_to_calc.set_cell(atoms_to_calc.get_cell() * scaling_factor, scale_atoms=True)
        
        print(f"  Calculating {crystal} at scaled lattice constant: {scaled_a:.4f} Å (Factor: {scaling_factor})")

        # Define the input settings for the Quantum ESPRESSO calculation
        input_settings = {
            'control': {
                'calculation': 'scf',
                'prefix': crystal,
                'outdir': output_dir,
            },
            'system': {
                'ecutwfc': fixed_ecutwfc_ry,
                'ecutrho': fixed_ecutwfc_ry * 8,
                'occupations': 'smearing',
                'smearing': 'gaussian',
                'degauss': 0.01,
                'input_dft': 'pbe',     # Specify PBE functional
                'vdw_corr': 'dft-d3',   # Apply Grimme D3 empirical dispersion correction
                'dftd3_version' : 4
            },
            'electrons': {
                'conv_thr': 1.0e-8,
            },
        }

        # Initialize the Espresso calculator
        calc = Espresso(pseudopotentials=calc_pseudopotentials,
                        input_data=input_settings,
                        kpts=fixed_kpts)
        
        atoms_to_calc.set_calculator(calc)

        try:
            energy = atoms_to_calc.get_potential_energy()
            # Store (actual_lattice_constant, energy) for the current crystal
            crystal_scan_results[crystal].append((scaled_a, energy))
            print(f"    Energy: {energy:.4f} eV")
        except Exception as e:
            # Store (actual_lattice_constant, "Failed")
            crystal_scan_results[crystal].append((scaled_a, "Failed"))
            print(f"    Calculation failed for {crystal} at scaling factor {scaling_factor}. Error: {e}")

print("\n" + "="*80)
print("All lattice constant scan calculations complete. Writing results to CSV.")
print("="*80)

# --- Output results to a CSV file ---
csv_filename = "lattice_constant_scan_results.csv" # Renamed to avoid confusion

with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Prepare header: 'Lattice Constant (Å)', 'TiN_Energy (eV)', 'ZrN_Energy (eV)', ...
    header = ['Lattice Constant (Å)']
    for crystal in crystals:
        header.append(f'{crystal}_Energy (eV)')
    writer.writerow(header)

    # Prepare data rows. We need to gather all unique lattice constants across all crystals
    # and then populate the energy for each crystal at those points.
    # Since we are using common scaling factors, the order of lattice constants will be consistent.
    
    # Iterate through the scaling factors to get the actual lattice constants
    # and then corresponding energies for each crystal
    for i, scaling_factor in enumerate(lattice_constant_scaling_factors):
        # We'll use the TiN lattice constant as the primary 'Lattice Constant (Å)' for the row,
        # assuming the scaling factors ensure comparable 'steps' across crystals.
        # It's more accurate to say each cell's actual 'a' is on that row.
        
        row_data = []
        
        # Get the actual lattice constant for the first crystal (e.g., TiN) for this row's primary column
        # and then append energies for all crystals.
        # This assumes that the `i`-th element in each crystal's `crystal_scan_results` corresponds
        # to the same scaled lattice constant point.
        
        # Add the primary lattice constant for the row (e.g., from TiN's actual scaled_a)
        if len(crystal_scan_results[crystals[0]]) > i:
             primary_latt_const = crystal_scan_results[crystals[0]][i][0] # (scaled_a, energy)
             row_data.append(f"{primary_latt_const:.4f}")
        else:
             row_data.append("N/A")

        for crystal in crystals:
            if len(crystal_scan_results[crystal]) > i:
                energy_val = crystal_scan_results[crystal][i][1] # (scaled_a, energy)
                if isinstance(energy_val, float):
                    row_data.append(f"{energy_val:.4f}")
                else:
                    row_data.append(energy_val) # Failed/N/A
            else:
                row_data.append("N/A") # Missing data for this point

        writer.writerow(row_data)

print(f"\nResults saved to {csv_filename}")
print("="*80)