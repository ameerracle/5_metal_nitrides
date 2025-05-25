from ase.io import read, write
from ase.visualize import view
from ase.calculators.espresso import Espresso
import os

crystals = ['TiN', 'ZrN', 'NbN', 'ScN', 'VN']
pseudopotentials = {'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
                    'Zr': 'zr_pbe_v1.uspp.F.UPF',
                    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
                    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                    'V': 'v_pbe_v1.4.uspp.F.UPF',
                    'N': 'N.pbe-n-radius_5.UPF'}

ecutwfc_values_ev = range(350, 551, 25)  # 350, 375, 400, ..., 550
eV_to_Ry = 1.0 / 13.6057  # Conversion factor from eV to Ry

# Dictionary to store energies for table printing
energies_table = {ecut: {} for ecut in ecutwfc_values_ev}

for crystal in crystals:
    # Create an output directory for each crystal
    output_dir = f'{crystal}_calculations'
    os.makedirs(output_dir, exist_ok=True) # Creates directory if it doesn't exist

    atoms = read(f'{crystal}.cif')

    # Ensure all necessary pseudopotentials are in the dictionary for the calculator
    calc_pseudopotentials = {}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            raise ValueError(f"Pseudopotential for {symbol} not found in the 'pseudopotentials' dictionary.")
        calc_pseudopotentials[symbol] = pseudopotentials[symbol]

    for ecutwfc_ev in ecutwfc_values_ev:
        # Convert to Ry for each cutoff value
        ecutwfc_ry = ecutwfc_ev * eV_to_Ry

        input_settings = {
            'control': {
                'calculation': 'scf',
                'prefix': crystal,
                'outdir': output_dir, # Set the output directory for QE
            },
            'system': {
                'ecutwfc': ecutwfc_ry,  # Wavefunction cutoff in Ry
                'ecutrho': ecutwfc_ry * 8,  # Charge density cutoff
                'occupations': 'smearing',
                'smearing': 'gaussian',
                'degauss': 0.01,
            },
            'electrons': {
                'conv_thr': 1.0e-8, # Example: add an electronic convergence threshold
            },
        }

        # Define k-points
        kpts = (5, 5, 5)

        # Initialize and run the calculator
        calc = Espresso(pseudopotentials=calc_pseudopotentials,
                        input_data=input_settings,
                        kpts=kpts)
        atoms.set_calculator(calc)

        # THIS IS THE CRUCIAL INDENTATION FIX:
        # The try-except block needs to be inside the inner loop
        # so that a calculation is attempted for each crystal and each ecutwfc.
        try:
            energy = atoms.get_potential_energy()
            energies_table[ecutwfc_ev][crystal] = energy
            print(f"Calculated {crystal} with ecutwfc: {ecutwfc_ev} eV, Energy: {energy:.4f} eV")
        except Exception as e:
            energies_table[ecutwfc_ev][crystal] = "Failed"
            print(f"Calculation failed for {crystal} with ecutwfc: {ecutwfc_ev} eV. Error: {e}")

# --- Output results to a file ---
output_filename = "crystal_e.txt"

with open(output_filename, 'w') as f: # Open the file in write mode ('w')
    f.write("\n" + "="*80 + "\n")
    f.write("Final Results:\n")
    f.write("="*80 + "\n")

    # Header
    header = f"{'ecutwfc (eV)':<15}" + "".join([f"{crystal:>15}" for crystal in crystals])
    f.write(header + "\n") # Write header to file
    f.write("-" * len(header) + "\n") # Write separator to file

    # Data rows
    for ecutwfc_ev in ecutwfc_values_ev:
        row = f"{ecutwfc_ev:<15}"
        for crystal in crystals:
            energy = energies_table[ecutwfc_ev].get(crystal, "N/A")
            if isinstance(energy, float):
                row += f"{energy:>15.4f}"
            else:
                row += f"{energy:>15}"
        f.write(row + "\n") # Write each data row to file

    f.write("="*80 + "\n") # Write final separator to file

print(f"\nResults saved to {output_filename}") # Inform the user that the file has been created
print("="*80)