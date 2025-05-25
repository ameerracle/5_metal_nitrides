from ase.io import read
from ase.visualize import view
from ase.calculators.espresso import Espresso

PSFs = ['Li2S', 'Li2S2', 'Li2S4', 'Li2S6', 'Li2S8', 'S8']
pseudopotentials = {'Li': 'li_pbe_v1.4.uspp.F.UPF',
                  'S': 's_pbe_v1.4.uspp.F.UPF'}

# Define ecutwfc values in eV (350 to 550 eV in steps of 25 eV)
ecutwfc_values_ev = range(350, 551, 25)  # 350, 375, 400, ..., 550
eV_to_Ry = 1.0 / 13.6057  # Conversion factor from eV to Ry

# K-points setup
kpts = (1, 1, 1)

# Dictionary to store results
results = {psf: {} for psf in PSFs}

# Print header
print("\n" + "="*80)
print(f"{'ecutwfc (eV)':<15}" + "".join([f"{psf:>15}" for psf in PSFs]))
print("-"*80)

# Iterate over ecutwfc values
for ecutwfc_ev in ecutwfc_values_ev:
    # Convert to Ry
    ecutwfc_ry = ecutwfc_ev * eV_to_Ry
    
    # Update input data with current ecutwfc
    input_data = {
        'control': {
            'calculation': 'scf',
        },
        'system': {
            'ecutwfc': ecutwfc_ry,  # Wavefunction cutoff in Ry
            'ecutrho': ecutwfc_ry * 8,  # Charge density cutoff
            'occupations': 'smearing',
            'smearing': 'marzari-vanderbilt',
            'degauss': 0.005,  # Smearing width
        },
    }
    
    # Print current ecutwfc value
    print(f"{ecutwfc_ev:<15}", end="")
    
    # Iterate over all PSFs
    for psf in PSFs:
        try:
            filelist = f'{psf}_D3.out'
            
            # Read the atoms from the file
            atoms = read(filelist) 
            atoms.center(vacuum=7)
            
            # Initialize the Espresso calculator for the current structure
            calc = Espresso(
                pseudopotentials=pseudopotentials,
                input_data=input_data,
                kpts=kpts,
            )
            
            # Attach the calculator to the atoms object
            atoms.set_calculator(calc)
            
            # Run the calculation
            energy = atoms.get_potential_energy()
            results[psf][ecutwfc_ev] = energy
            print(f"{energy:>15.4f}", end="")
            
        except Exception as e:
            print(f"{'Error':>15}", end="")
            print(f"\nError processing {psf} with ecutwfc={ecutwfc_ev}eV: {str(e)}")
    
    print()  # New line after each ecutwfc value

# Print final results in a nice table
print("\n" + "="*80)
print("Final Results:")
print("="*80)
print(f"{'ecutwfc (eV)':<15}" + "".join([f"{psf:>15}" for psf in PSFs]))
print("-"*80)

for ecutwfc_ev in ecutwfc_values_ev:
    print(f"{ecutwfc_ev:<15}", end="")
    for psf in PSFs:
        if ecutwfc_ev in results[psf]:
            print(f"{results[psf][ecutwfc_ev]:>15.4f}", end="")
        else:
            print(f"{'N/A':>15}", end="")
    print()

print("-"*80)
print("All calculations completed.")
    
