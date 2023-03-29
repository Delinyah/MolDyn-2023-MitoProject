import os
import subprocess

root_folder = "~/Project/OMM/structures-outer-new"

for folder in os.listdir(os.path.expanduser(root_folder)):
    folder_path = os.path.join(root_folder, folder)
    
    if os.path.isdir(folder_path):
        gro_file = os.path.join(folder_path, "md.gro")
        xtc_file = os.path.join(folder_path, "md.xtc")
        
        if os.path.isfile(gro_file) and os.path.isfile(xtc_file):
            print(f"Running fingerprinting for: {folder}")
            subprocess.run(["python", "general_fingerprinting.py", gro_file, xtc_file])
