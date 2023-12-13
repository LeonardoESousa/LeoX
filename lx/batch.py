#!/usr/bin/env python3
import os
import sys
import shutil
import lx.tools

BATCH_FILE = sys.argv[1]
GAUSSIAN = sys.argv[2]
NUM = int(sys.argv[3])

def run_batch():
    scripts = [i for i in os.listdir(".") if ".sh" in i]
    for file in scripts:
        shutil.copy(file, "Geometries")
    os.chdir("Geometries")
    files = [i for i in os.listdir(".") if "Geometr" in i and ".com" in i]
    factor = lx.tools.set_factor(files[0])
    the_watcher = lx.tools.Watcher('.', counter=factor)
    the_watcher.run(BATCH_FILE, GAUSSIAN, NUM)
