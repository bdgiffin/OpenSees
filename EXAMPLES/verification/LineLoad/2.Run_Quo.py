import os
import pathlib
from Params import *
import sys
args = " ".join(sys.argv[1:])  # Capture all arguments passed to the script
command = f'powershell.exe wsl python3.12 2.Model_QuoFEM.py {args}'
os.system(command)