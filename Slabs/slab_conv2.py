from ase.io import read, write
from ase.atoms import Atoms 
from ase.build import surface 
from ase.visualize import view 
from ase.constraints import FixAtoms 
from ase.calculators.espresso import Espresso
import os
import csv 

# --- Define Pseudopotential Filenames based on PP_list_SSP.txt ---
# Make sure these filenames match exactly what's in your PP_list_SSP.txt
pseudopotential_files = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
    'Zr': 'Zr_pbe_v1.uspp.F.UPF',
    'N':  'N.pbe-n-radius_5.UPF'
}

# Define the list of materials to process (all except TiN, as requested)
materials_to_process = ['TiN']

# Path to the directory where the slab files are located
slab_directory = 'Slabs' 

# Fixed Ecutwfc for k-point convergence test
fixed_ecutwfc = 45 # Ry, a reasonably high value for initial k-point test

# Quantum ESPRESSO Calculator Settings (common for all tests, adjust only ecutrho in loop)
# Basic parameters for a self-consistent field (SCF) calculation
input_data = {
    'control': {
        'calculation': 'scf',
        'outdir': './QE_output',


    },
    'system': {
        # 'ecutrho' will be adjusted within each test loop
        'smearing': 'mv', # Methfessel-Paxton smearing
        'degauss': 0.01, # Smearing width in Ry
        'input_dft': 'pbe',         # Specify PBE functional
        'vdw_corr': 'dft-d3',
        'dftd3_version' : 4 # Apply Grimme D3 empirical dispersion correction
        # 'prefix' will be set inside the loop for each material
    },
    'electrons': {
        'conv_thr': 1.0e-6, # SCF convergence threshold in Ry
        'mixing_mode': 'TF', # <-- Add this line: Use Thomas-Fermi mixing for metals
        'mixing_beta': 0.4,  # <-- Add this line: Adjust mixing beta (e.g., 
        'electron_maxstep': 200, # Directory for QE temporary files and output
    }
}

for test_material in materials_to_process:
    print(f"\n--- Processing Material: {test_material} ---")

    # Load the already-made slab structure using os.path.join for OS-agnosticism
    slab_filename = os.path.join(slab_directory, f'{test_material}_slab_100.traj')
    print(f"Loading slab from: {slab_filename}")
    
    slab = read(slab_filename)

    # Update prefix for the current material
    input_data['control']['prefix'] = test_material
    input_data['system']['ecutrho'] = 8 * fixed_ecutwfc # Set initial ecutrho for k-point test

    # --- 2. K-point Convergence Test ---
    print(f"\nStarting k-point convergence test for {test_material}...")
    # K-point range limited from (3,3,1) to (7,7,1) as requested
    k_points_range = [(3, 3, 1), (4, 4, 1), (5, 5, 1), (6, 6, 1)] 

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

    # --- 3. Ecutwfc Convergence Test ---
    print(f"\nStarting ecutwfc convergence test for {test_material}...")
    # Ecutwfc range limited from 35 to 60 as requested
    ecutwfc_range = [40, 45, 50, 55, 60] 

    # Use a reasonably dense k-point mesh (e.g., one from the converged k-point test, or start with 7x7x1)
    converged_kpts = (5, 5, 1) # Retaining this value as per previous context

    ecut_energies = []
    for ecut in ecutwfc_range:
        print(f"    Running with ecutwfc: {ecut} Ry and k-points: {converged_kpts}")
        
        # Adjust ecutrho based on ecutwfc
        input_data['system']['ecutrho'] = 8 * ecut # Corrected to 8x as per initial input_data setup
        
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

print("\nAll convergence tests completed for specified materials. Data saved to CSV files for plotting.")