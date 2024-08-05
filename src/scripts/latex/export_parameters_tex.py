import src.constants
import subprocess

def copy2clip(txt):
    cmd='echo '+txt.strip()+'|clip'
    return subprocess.check_call(cmd, shell=True)

parameters = src.constants.BISTABLE_STRONG_PDMP_PARAMS

string = ""
for key, val in parameters.items():
    key_by_part = key.split("_")
    if len(key_by_part) == 1:
        key = key_by_part[0]
    elif len(key_by_part) == 2:
        key = f"{key_by_part[0]}_{{{key_by_part[1]}}}"
    else: 
        key = f"{key_by_part[0]}_{{{key_by_part[1]}}}^{{{key_by_part[2]}}}"
    print(key)
    string = f"{string} ${key} = {val},$"
copy2clip(string[:-1])
print(string)
