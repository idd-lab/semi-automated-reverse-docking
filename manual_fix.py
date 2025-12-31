import os
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_B"
INPUT_DIR = os.path.join(BASE_PATH, "cleaned_pdbs")
OUTPUT_DIR = os.path.join(BASE_PATH, "fixed_pdbs")

def fix_specific_pdb(user_input):
    # Handle inputs like "Q96DI7" or "Q96DI7.pdb"
    clean_id = user_input.strip().replace(".pdb", "")
    filename = f"{clean_id}.pdb"
    
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)

    if not os.path.exists(input_path):
        print(f"❌ ERROR: File not found: {input_path}")
        return

    print(f"\n--- Processing {clean_id} ---")
    try:
        print("1. Loading PDB...")
        fixer = PDBFixer(filename=input_path)

        # === THE CRITICAL FIX ===
        # We manually empty this dictionary. 
        # This stops PDBFixer from trying to rebuild huge missing loops (the cause of the freeze)
        # AND prevents the 'AttributeError' crash.
        print("2. Applying Safety Override (Skipping Loop Building)...")
        fixer.missingResidues = {} 

        print("3. Finding missing atoms...")
        fixer.findMissingAtoms()

        print("4. Adding missing atoms...")
        fixer.addMissingAtoms()

        print("5. Adding hydrogens (pH 7.4)...")
        fixer.addMissingHydrogens(pH=7.4)

        print("6. Saving...")
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"✅ SUCCESS! Saved to: {output_path}")

    except Exception as e:
        print(f"❌ FAILED: {str(e)}")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"MANUAL PDB REPAIR TOOL")
    print(f"Input Directory: {INPUT_DIR}")
    print(f"Output Directory: {OUTPUT_DIR}")
    print("---------------------------------------------------")
    print("Type the PDB ID (e.g., Q96DI7) and press Enter.")
    print("Type 'exit' or 'q' to quit.")
    print("---------------------------------------------------")

    while True:
        try:
            user_input = input("\n> Enter Filename/ID: ")
            
            if user_input.lower() in ['exit', 'q', 'quit']:
                print("Exiting...")
                break
            
            if not user_input.strip():
                continue

            fix_specific_pdb(user_input)
            
        except KeyboardInterrupt:
            print("\nExiting...")
            break

if __name__ == "__main__":
    main()