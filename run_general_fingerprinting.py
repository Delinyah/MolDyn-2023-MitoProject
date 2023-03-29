#!/usr/bin/env python2

import os
import subprocess

root_folder = "~/Project/OMM/structures-outer-new"

for folder in os.listdir(os.path.expanduser(root_folder)):
    folder_path = os.path.join(os.path.expanduser(root_folder), folder)

    if os.path.isdir(folder_path):
        gro_file = None
        xtc_file = None

        for file in os.listdir(folder_path):
            if file == "md.gro":
                gro_file = os.path.join(folder_path, file)
            elif file == "md.xtc":
                xtc_file = os.path.join(folder_path, file)

        if gro_file and xtc_file:
            print("Running fingerprinting for: {}".format(folder))
            subprocess.call(["python2", "general_fingerprinting.py", gro_file, xtc_file])
        else:
            print("No matching files found in: {}".format(folder))
    else:
        print("Not a directory: {}".format(folder))
