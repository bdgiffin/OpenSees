import os
import pathlib
from Params1 import *
import sys
args = " ".join(sys.argv[1:])  # Capture all arguments passed to the script
command = f'powershell.exe wsl "source ~/.myvenv/bin/activate ; python 2.Model_QuoFEM.py {args}"'
os.system(command)