from ase.io import read, write
from ase.build import surface
from ase.constraints import FixAtoms
from ase.calculators.espresso import Espresso
import os
import csv # Import the csv module for saving data

# --- Define Pseudopotential Filenames based on PP_list_SSP.txt ---
# Make sure these filenames match exactly what's in your PP_list_SSP.txt
pseudopotential_files = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V': 'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
    'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
    'Zr': 'Zr_pbe_v1.uspp.F.UPF',
    'N': 'N.pbe-n-radius_5.UPF'
}

# --- Parameters you'll need to define ---
# Path to your pseudopotentials directory on your server
# (This should be the directory where the above .UPF files are located)
pseudo_dir = '/path/to/your/pseudopotentials/directory' # <--- IMPORTANT: UPDATE THIS PATH

# Choose the slab to test (e.g., TiN)
test_material = 'TiN' # You can change this to 'VN', 'ScN', 'NbN', or 'ZrN'
fixed_index_dir = './Fixed_Index_Crystals'
slab_dir = './Slabs' # Directory where slabs will be generated and saved for convergence tests
os.makedirs(slab_dir, exist_ok=True)

# Load the reordered crystal
reordered_crystal_path = os.path.join(fixed_index_dir, f'{test_material}_reordered.xyz')

if not os.path.exists(reordered_crystal_path):
    print(f"Error: Reordered crystal for {test_material} not found at {reordered_crystal_path}.")
    print("Please ensure you ran the previous script to fix atom indexing and save .xyz files.")
    exit() # Exit if the required file is not found

crystal = read(reordered_crystal_path)

# --- Re-generate the slab with 2 layers (as per your last request) ---
# Define Miller index (e.g., (100))
miller_index = (1, 0, 0)

# Create the slab structure: 2 layers, 10.0 Angstrom vacuum
slab = surface(crystal, miller_index, layers=2, vacuum=10.0)

# Create a 2x2 supercell
slab = slab * (2, 2, 1)

# Apply constraint: Fix the bottom 1 layer (for 2 total layers, 1 fixed, 1 unfixed)
z_coords = slab.get_positions()[:, 2]
unique_z = sorted(list(set(z_coords)))
fixed_z_threshold = unique_z[0] # The bottom-most layer

fixed_indices = [atom.index for atom in slab if atom.position[2] <= fixed_z_threshold]
slab.set_constraint(FixAtoms(indices=fixed_indices))

# Save this specific slab for consistency (optional, but good for tracking)
write(os.path.join(slab_dir, f'{test_material}_slab_test.traj'), slab)
print(f"Slab for {test_material} generated with {slab.get_global_number_of_atoms()} atoms for convergence testing.")
print(f"Number of fixed atoms: {len(fixed_indices)}")

# --- Quantum ESPRESSO Calculator Settings (common for both tests) ---
# Basic parameters for a self-consistent field (SCF) calculation
input_data = {
    'control': {
        'calculation': 'scf',
        'pseudo_dir': pseudo_dir,
        'outdir': './QE_output', # Directory for QE temporary files and output
        'prefix': test_material,
        'tprnfor': True, # Print forces
        'tstress': True, # Print stress
    },
    'system': {
        'ecutrho': 4 * 40, # Common for norm-conserving, adjust for USPP/PAW (e.g., 8-12 * ecutwfc)
        'smearing': 'mv', # Methfessel-Paxton smearing
        'degauss': 0.01, # Smearing width in Ry
    },
    'electrons': {
        'conv_thr': 1.0e-8, # SCF convergence threshold in Ry
        'mixing_beta': 0.3, # Mixing factor for charge density
    }
}

# Ensure the output directory exists
os.makedirs(input_data['control']['outdir'], exist_ok=True)

# --- 2. K-point Convergence Test ---
print("\nStarting k-point convergence test...")
k_points_range = [(3, 3, 1), (4, 4, 1), (5, 5, 1), (6, 6, 1), (7, 7, 1), (8, 8, 1), (9, 9, 1)] # Example range
fixed_ecutwfc = 40 # Ry, a reasonably high value for initial k-point test

k_point_energies = []
for kpt in k_points_range:
    print(f"  Running with k-points: {kpt} and ecutwfc: {fixed_ecutwfc} Ry")
    
    calc = Espresso(
        input_data=input_data,
        kpts=kpt,
        ecutwfc=fixed_ecutwfc,
        pseudopotentials=pseudopotential_files # Using the specific PP files now
    )
    slab.set_calculator(calc)
    
    try:
        energy = slab.get_potential_energy() # Get total energy in eV
        energy_per_atom = energy / slab.get_global_number_of_atoms() # Energy per atom
        k_point_energies.append(energy_per_atom)
        print(f"  Total energy per atom: {energy_per_atom:.4f} eV")
    except Exception as e:
        print(f"  Calculation failed for k-points {kpt}: {e}")
        k_point_energies.append(None) # Mark as failed

# Save k-point convergence data to CSV
csv_kpoint_filename = os.path.join(slab_dir, f'{test_material}_kpoint_convergence.csv')
with open(csv_kpoint_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['K-point Mesh (nx)', 'K-point Mesh (ny)', 'K-point Mesh (nz)', 'Energy per Atom (eV)'])
    for i, kpt in enumerate(k_points_range):
        if k_point_energies[i] is not None:
            writer.writerow([kpt[0], kpt[1], kpt[2], k_point_energies[i]])
        else:
            writer.writerow([kpt[0], kpt[1], kpt[2], 'Failed'])
print(f"K-point convergence data saved to {csv_kpoint_filename}")

# --- 3. Ecutwfc Convergence Test ---
print("\nStarting ecutwfc convergence test...")
ecutwfc_range = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70] # Example range in Ry

# Use a reasonably dense k-point mesh (e.g., one from the converged k-point test, or start with 6x6x1)
# IMPORTANT: After your k-point test, replace this with the actual converged k-point mesh!
converged_kpts = (6, 6, 1) # <--- UPDATE THIS after running K-point test

ecut_energies = []
for ecut in ecutwfc_range:
    print(f"  Running with ecutwfc: {ecut} Ry and k-points: {converged_kpts}")
    
    # Adjust ecutrho based on ecutwfc, assuming a 4x ratio for now (typical for norm-conserving PPs)
    # If using USPP/PAW, you might need 8-12x. Check your pseudopotential type from PP_list_SSP.txt.
    input_data['system']['ecutrho'] = 4 * ecut
    
    calc = Espresso(
        input_data=input_data,
        kpts=converged_kpts,
        ecutwfc=ecut,
        pseudopotentials=pseudopotential_files # Using the specific PP files now
    )
    slab.set_calculator(calc)
    
    try:
        energy = slab.get_potential_energy() # Get total energy in eV
        energy_per_atom = energy / slab.get_global_number_of_atoms() # Energy per atom
        ecut_energies.append(energy_per_atom)
        print(f"  Total energy per atom: {energy_per_atom:.4f} eV")
    except Exception as e:
        print(f"  Calculation failed for ecutwfc {ecut}: {e}")
        ecut_energies.append(None) # Mark as failed

# Save ecutwfc convergence data to CSV
csv_ecutwfc_filename = os.path.join(slab_dir, f'{test_material}_ecutwfc_convergence.csv')
with open(csv_ecutwfc_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Ecutwfc (Ry)', 'Energy per Atom (eV)'])
    for i, ecut in enumerate(ecutwfc_range):
        if ecut_energies[i] is not None:
            writer.writerow([ecut, ecut_energies[i]])
        else:
            writer.writerow([ecut, 'Failed'])
print(f"Ecutwfc convergence data saved to {csv_ecutwfc_filename}")

print("\nConvergence tests completed. Data saved to CSV files for plotting.")