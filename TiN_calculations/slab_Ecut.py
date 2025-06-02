# ecutwfc_convergence.py
from ase.io import read
from ase.calculators.espresso import Espresso
import os
import csv 

# Define Pseudopotential Filenames based on PP_list_SSP.txt 
# Make sure these filenames match exactly what's in your PP_list_SSP.txt
pseudopotential_files = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
    'Zr': 'Zr_pbe_v1.uspp.F.UPF',
    'N':  'n_pbe_v1.2.uspp.F.UPF'
}

# Define the list of materials to process
materials_to_process = ['TiN']

# Quantum ESPRESSO Calculator Settings (common for all tests)
# Basic parameters for a self-consistent field (SCF) calculation
input_data = {
    'control': {
        'calculation': 'scf',
        'outdir': './QE_output',
    },
    'system': {

        'input_dft': 'pbe', # Specify PBE functional
        'vdw_corr': 'dft-d3',
        'smearing': 'gaussian', # Methfessel-Paxton smearing
        'degauss': 0.02, # Smearing width in Ry
    },
    'electrons': {
        'conv_thr': 1.0e-6, # SCF convergence threshold in Ry
        'mixing_mode': 'local-TF', # Use Thomas-Fermi mixing for metals
        'mixing_beta': 0.3,  # Adjust mixing beta
        'electron_maxstep': 100, # Directory for QE temporary files and output
    }
}

for test_material in materials_to_process:
    print(f"\n--- Processing Material: {test_material} ---")

    # Load the already-made slab structure from the current directory
    slab_filename = 'TiN_slab_D3.in'
    print(f"Loading slab from: {slab_filename}")
    
    slab = read(slab_filename)

    # Update prefix for the current material
    input_data['control']['prefix'] = test_material

    # Ecutwfc Convergence Test
    print(f"\nStarting ecutwfc convergence test for {test_material}...")
    # Ecutwfc range limited from 35 to 60 as requested
    ecutwfc_range = [40] 

    # Use a reasonably dense k-point mesh (e.g., one from the converged k-point test, or start with 7x7x1)
    converged_kpts = (3, 3, 1) # Retaining this value as per previous context

    ecut_energies = []
    for ecut in ecutwfc_range:
        print(f"    Running with ecutwfc: {ecut} Ry and k-points: {converged_kpts}")
        
        # Adjust ecutrho based on ecutwfc
        input_data['system']['ecutrho'] = 300 # Corrected to 8x as per initial input_data setup
        
        calc = Espresso(
            input_data=input_data,
            kpts=converged_kpts,
            ecutwfc=ecut,
            pseudopotentials=pseudopotential_files # Using the specific PP files now
        )
        slab.set_calculator(calc)
        
        energy = slab.get_potential_energy() # Get total energy in eV
        ecut_energies.append(energy)
        print(f"    Total energy: {energy:.4f} eV")

    # Save ecutwfc convergence data to CSV in the current directory
    csv_ecutwfc_filename = f'{test_material}_ecutwfc_convergence.csv'
    with open(csv_ecutwfc_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Ecutwfc (Ry)', 'Total Energy (eV)'])
        for i, ecut in enumerate(ecutwfc_range):
            writer.writerow([ecut, ecut_energies[i]])
    print(f"Ecutwfc convergence data saved to {csv_ecutwfc_filename}")

print("\nEcutwfc convergence tests completed for specified materials. Data saved to CSV files for plotting.")