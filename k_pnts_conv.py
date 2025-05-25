from ase.io import read
from ase.calculators.espresso import Espresso
import os
import csv # For writing CSV files

# --- System and Pseudopotential Definitions ---
crystals = ['TiN', 'ZrN', 'NbN', 'ScN', 'VN']
pseudopotentials = {'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
                    'Zr': 'Zr_pbe_v1.uspp.F.UPF',
                    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
                    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                    'V': 'v_pbe_v1.4.uspp.F.UPF',
                    'N': 'N.pbe-n-radius_5.UPF'}

# --- Fixed Parameters for K-point Convergence ---
fixed_ecutwfc_ev = 600 # Fixed energy cutoff for k-point convergence
eV_to_Ry = 1.0 / 13.6057  # Conversion factor from eV to Ry
fixed_ecutwfc_ry = fixed_ecutwfc_ev * eV_to_Ry

# Define the k-point mesh values to test
kpts_values = [(i, i, i) for i in range(3, 8)] # (3,3,3) to (7,7,7)

# Dictionary to store energies for CSV output
# Structure: {kpt_tuple: {crystal_name: energy, ...}}
energies_data_for_csv = {kpt: {} for kpt in kpts_values}

print(f"Starting K-point convergence test with fixed ecutwfc = {fixed_ecutwfc_ev} eV.")
print(f"Testing k-point meshes: {kpts_values}")
print("="*80)

for crystal in crystals:
    # Create an output directory for each crystal's Quantum ESPRESSO calculations
    # This keeps output files organized by crystal
    output_dir = f'{crystal}_k_point_calculations'
    os.makedirs(output_dir, exist_ok=True) # Creates the directory if it does not already exist

    # Read the crystal structure from a .cif file
    # Ensure your .cif files (e.g., TiN.cif, ZrN.cif) are in the script's directory
    atoms = read(f'{crystal}.cif')

    # Prepare the pseudopotentials dictionary for the current crystal
    calc_pseudopotentials = {}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            raise ValueError(f"Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary.")
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    for kpt_tuple in kpts_values:
        # Define the input settings for the Quantum ESPRESSO calculation
        input_settings = {
            'control': {
                'calculation': 'scf',  # Self-consistent field calculation
                'prefix': crystal,     # Prefix for output files (e.g., TiN.scf.out)
                'outdir': output_dir,  # Specify the output directory for QE files
            },
            'system': {
                'ecutwfc': fixed_ecutwfc_ry, # Use the fixed energy cutoff
                'ecutrho': fixed_ecutwfc_ry * 8, # Charge density cutoff
                'occupations': 'smearing',
                'smearing': 'gaussian',
                'degauss': 0.01,
                'input_dft': 'pbe',         # Specify PBE functional
                'vdw_corr': 'dft-d3',
                'dftd3_version' : 4 # Apply Grimme D3 empirical dispersion correction
            },
            'electrons': {
                'conv_thr': 1.0e-8,    # Electronic convergence threshold
            },
        }

        # Initialize the Espresso calculator
        calc = Espresso(pseudopotentials=calc_pseudopotentials,
                        input_data=input_settings,
                        kpts=kpt_tuple) # Use the current k-point tuple
        
        atoms.set_calculator(calc)

        try:
            # Run the calculation and get the potential energy
            energy = atoms.get_potential_energy()
            energies_data_for_csv[kpt_tuple][crystal] = energy
            print(f"Calculated {crystal} with kpts: {kpt_tuple}, Energy: {energy:.4f} eV")
        except Exception as e:
            # Handle cases where the calculation might fail
            energies_data_for_csv[kpt_tuple][crystal] = "Failed"
            print(f"Calculation failed for {crystal} with kpts: {kpt_tuple}. Error: {e}")

print("\n" + "="*80)
print("All calculations complete. Writing results to CSV.")
print("="*80)

# --- Output results to a CSV file ---
csv_filename = "k_pnts.csv"

with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write the header row
    header = ['k_points'] + crystals
    writer.writerow(header)

    # Write data rows
    # Sort kpts_values for consistent output order
    sorted_kpts_values = sorted(kpts_values) 
    for kpt_tuple in sorted_kpts_values:
        # Represent k-point tuple as a string, e.g., "(3, 3, 3)"
        row_data = [str(kpt_tuple)] 
        for crystal in crystals:
            energy = energies_data_for_csv[kpt_tuple].get(crystal, "N/A")
            if isinstance(energy, float):
                row_data.append(f"{energy:.4f}") # Format float for CSV
            else:
                row_data.append(energy) # Keep "Failed" or "N/A" as string
        writer.writerow(row_data)

print(f"\nResults saved to {csv_filename}")
print("="*80)