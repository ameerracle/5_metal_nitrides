from ase.io import read, write
from ase.visualize import view
from ase.calculators.espresso import Espresso
from ase.optimize import BFGS # For optimization
import os
import csv # For writing CSV files

crystals = ['TiN', 'ZrN', 'NbN', 'ScN', 'VN']
pseudopotentials = {'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
                    'Zr': 'Zr_pbe_v1.uspp.F.UPF',
                    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
                    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                    'V': 'v_pbe_v1.4.uspp.F.UPF',
                    'N': 'N.pbe-n-radius_5.UPF'}

ecutwfc_values_ev = range(500, 651, 25) # From 500 to 600 eV, in steps of 25 eV
eV_to_Ry = 1.0 / 13.6057  # Conversion factor from eV to Ry

# Dictionary to store energies for CSV output
# The structure will be: {ecutwfc: {crystal_name: energy, ...}}
energies_data_for_csv = {ecut: {} for ecut in ecutwfc_values_ev}

for crystal in crystals:
    # Create an output directory for each crystal's Quantum ESPRESSO calculations
    output_dir = f'{crystal}_calculations'
    os.makedirs(output_dir, exist_ok=True) # Creates the directory if it does not already exist

    # Read the crystal structure from a .cif file
    # Note: If your .cif files are in a different directory, adjust the path here.
    atoms = read(f'{crystal}.cif')

    # Prepare the pseudopotentials dictionary for the current crystal
    calc_pseudopotentials = {}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            raise ValueError(f"Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary.")
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    for ecutwfc_ev in ecutwfc_values_ev:
        # Convert the wavefunction cutoff energy from eV to Ry for Quantum ESPRESSO
        ecutwfc_ry = ecutwfc_ev * eV_to_Ry

        # Define the input settings for the Quantum ESPRESSO calculation
        # Adding PBE (input_dft) and Grimme D3 correction (vdw_corr='D3')
        input_settings = {
            'control': {
                'calculation': 'scf',  # Self-consistent field calculation
                'prefix': crystal,     # Prefix for output files (e.g., TiN.scf.out)
                'outdir': output_dir,  # Specify the output directory for QE files
            },
            'system': {
                'ecutwfc': ecutwfc_ry,      # Wavefunction cutoff energy in Ry
                'ecutrho': ecutwfc_ry * 8,  # Charge density cutoff (typically 4-8 times ecutwfc)
                'occupations': 'smearing',  # Method for smearing electron occupations
                'smearing': 'gaussian',     # Type of smearing (Gaussian distribution)
                'degauss': 0.01,            # Smearing width in Ry
                'input_dft': 'pbe',         # Specify PBE functional
                'vdw_corr': 'dft-d3',
                'dftd3_version' : 4 # Apply Grimme D3 empirical dispersion correction
            },
            'electrons': {
                'conv_thr': 1.0e-8,         # Electronic convergence threshold
            },
        }

        # Define the k-point mesh for the Brillouin zone sampling
        kpts = (5, 5, 5)

        # Initialize the Espresso calculator with the pseudopotentials, input settings, and k-points
        calc = Espresso(pseudopotentials=calc_pseudopotentials,
                        input_data=input_settings,
                        kpts=kpts)
        
        # Attach the calculator to the atoms object
        atoms.set_calculator(calc)

        try:
            # Run the calculation and get the potential energy
            # Since we're doing SCF for convergence, we just get energy.
            # If you wanted to fully optimize geometry, you'd use BFGS here.
            # For this cutoff convergence, we assume fixed geometry.
            energy = atoms.get_potential_energy()
            energies_data_for_csv[ecutwfc_ev][crystal] = energy
            print(f"Calculated {crystal} with ecutwfc: {ecutwfc_ev} eV, Energy: {energy:.4f} eV")
        except Exception as e:
            # Handle cases where the calculation might fail
            energies_data_for_csv[ecutwfc_ev][crystal] = "Failed"
            print(f"Calculation failed for {crystal} with ecutwfc: {ecutwfc_ev} eV. Error: {e}")

# --- Output results to a CSV file ---
csv_filename = "crystal_e_convergence.csv"

with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write the header row
    header = ['ecutwfc (eV)'] + crystals
    writer.writerow(header)

    # Write data rows
    for ecutwfc_ev in ecutwfc_values_ev:
        row_data = [ecutwfc_ev]
        for crystal in crystals:
            energy = energies_data_for_csv[ecutwfc_ev].get(crystal, "N/A")
            if isinstance(energy, float):
                row_data.append(f"{energy:.4f}") # Format float for CSV
            else:
                row_data.append(energy) # Keep "Failed" or "N/A" as string
        writer.writerow(row_data)

print(f"\nResults saved to {csv_filename}")
print("="*80)