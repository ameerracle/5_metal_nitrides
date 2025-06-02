# kpoint_convergence.py
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
    'N':  'N.pbe-n-radius_5.UPF'
}

# Define the list of materials to process
materials_to_process = ['TiN']

# Fixed Ecutwfc for k-point convergence test
fixed_ecutwfc = 40 # Ry, a reasonably high value for initial k-point test

# Quantum ESPRESSO Calculator Settings (common for all tests)
# Basic parameters for a self-consistent field (SCF) calculation
input_data = {
    'control': {
        'calculation': 'scf',
        'outdir': './QE_output',
    },
    'system': {
        'smearing': 'mv', # Methfessel-Paxton smearing
        'degauss': 0.01, # Smearing width in Ry
        'input_dft': 'pbe', # Specify PBE functional
        'vdw_corr': 'dft-d3',
        'dftd3_version' : 4 # Apply Grimme D3 empirical dispersion correction
    },
    'electrons': {
        'conv_thr': 1.0e-6, # SCF convergence threshold in Ry
        'mixing_mode': 'TF', # Use Thomas-Fermi mixing for metals
        'mixing_beta': 0.4,  # Adjust mixing beta
        'electron_maxstep': 200, # Directory for QE temporary files and output
    }
}

for test_material in materials_to_process:
    print(f"\n--- Processing Material: {test_material} ---")

    # Load the already-made slab structure from the current directory
    slab_filename = f'{test_material}_slab_100.traj'
    print(f"Loading slab from: {slab_filename}")
    
    slab = read(slab_filename)

    # Update prefix for the current material
    input_data['control']['prefix'] = test_material
    input_data['system']['ecutrho'] = 9 * fixed_ecutwfc # Set initial ecutrho for k-point test

    # K-point Convergence Test
    print(f"\nStarting k-point convergence test for {test_material}...")
    # K-point range limited from (3,3,1) to (7,7,1) as requested
    k_points_range = [(3, 3, 1), (5, 5, 1), (4, 4, 1), (6, 6, 1)] 

    k_point_energies = []
    for kpt in k_points_range:
        print(f"    Running with k-points: {kpt} and ecutwfc: {fixed_ecutwfc} Ry")
        calc = Espresso(
            input_data=input_data,
            kpts=kpt,
            ecutwfc=fixed_ecutwfc,
            pseudopotentials=pseudopotential_files # Using the specific PP files now
        )
        slab.set_calculator(calc)
        
        energy = slab.get_potential_energy() # Get total energy in eV
        k_point_energies.append(energy)
        print(f"    Total energy: {energy:.4f} eV")

    # Save k-point convergence data to CSV in the current directory
    csv_kpoint_filename = f'{test_material}_kpoint_convergence.csv'
    with open(csv_kpoint_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['K-point Mesh (nx)', 'K-point Mesh (ny)', 'K-point Mesh (nz)', 'Total Energy (eV)'])
        for i, kpt in enumerate(k_points_range):
            writer.writerow([kpt[0], kpt[1], kpt[2], k_point_energies[i]])
    print(f"K-point convergence data saved to {csv_kpoint_filename}")

print("\nK-point convergence tests completed for specified materials. Data saved to CSV files for plotting.")