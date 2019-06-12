import multiprocessing
import subprocess
import sys

def __os():
    if sys.platform.startswith("win"):
        return "WINDOWS"
    elif sys.platform.startswith("cygwin"):
        return "CYGWIN"
    elif sys.platform.startswith("linux"):
        return "LINUX"
    elif sys.platform.startswith("darwin"):
        return "OSX"
    
def __which(exe):
    platform = __os()
    if platform == "WINDOWS":
        cmd = "where.exe"
    else:
        cmd = "which"
        
    try:
        path = subprocess.check_output([cmd, exe], universal_newlines=True).strip()
    except subprocess.CalledProcessError:
        path = ""
        
    return path

def proc_err_output(output):
    lines = output.split("\n")
    idx = None
    for l in lines:
        if l.startswith("ERROR"):
            idx = lines.index(l)
    if idx is None:
        return output
    else:
        errline = lines[idx]
        # strip "ERROR: "
        lines[idx] = errline[errline.index(" ")+1:]
        # strip "Bertini will now exit due to this error"
        return "\n".join(lines[idx:-1])

BERTINI = __which("bertini")
MPIRUN  = __which("mpirun")
PCOUNT  = multiprocessing.cpu_count()
