from ase.io import read, write
from ase.atoms import Atoms # Keep if potentially used for other parts, though not for direct slab making here
from ase.build import surface # Keep if needed for *future* slab generation, but not used in this specific run
from ase.visualize import view # Keep if user wants to view the loaded slab
from ase.constraints import FixAtoms # Keep if constraints are loaded from .traj or needed for future ops
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


# Choose the slab to test (e.g., TiN)
test_material = 'TiN' # You can change this to 'VN', 'ScN', 'NbN', or 'ZrN'

slab = read('TiN_slab_100.traj')

# Optional: View the loaded slab in ASE GUI (uncomment if you want to see it)
# view(slab) 

fixed_ecutwfc = 40 # Ry, a reasonably high value for initial k-point test
# --- Quantum ESPRESSO Calculator Settings (common for both tests) ---
# Basic parameters for a self-consistent field (SCF) calculation
input_data = {
    'control': {
        'calculation': 'scf',
        'outdir': './QE_output', # Directory for QE temporary files and output
        'prefix': test_material,
    },
    'system': {
        'ecutrho': 8 * fixed_ecutwfc, # Common for norm-conserving, adjust for USPP/PAW (e.g., 8-12 * ecutwfc)
        'smearing': 'mv', # Methfessel-Paxton smearing
        'degauss': 0.01, # Smearing width in Ry
    },
    'electrons': {
        'conv_thr': 1.0e-7, # SCF convergence threshold in Ry
    }
}

# Ensure the output directory for QE exists
os.makedirs(input_data['control']['outdir'], exist_ok=True)

# --- 2. K-point Convergence Test ---
print("\nStarting k-point convergence test...")
k_points_range = [(3, 3, 1), (4, 4, 1), (5, 5, 1), (6, 6, 1), (7, 7, 1), (8, 8, 1), (9, 9, 1)] # Example range


k_point_energies = []
for kpt in k_points_range:
    print(f"  Running with k-points: {kpt} and ecutwfc: {fixed_ecutwfc} Ry")
    
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
        print(f"  Total energy per atom: {energy_per_atom:.4f} eV")
    except Exception as e:
        print(f"  Calculation failed for k-points {kpt}: {e}")
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
    print(f"  Running with ecutwfc: {ecut} Ry and k-points: {converged_kpts}")
    
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
        print(f"  Total energy per atom: {energy_per_atom:.4f} eV")
    except Exception as e:
        print(f"  Calculation failed for ecutwfc {ecut}: {e}")
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