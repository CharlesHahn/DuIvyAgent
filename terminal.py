""" 

"""


import subprocess


def run_terminal(command: str): 
    """run the terminal command and output status, output, error"""
    res = {"returncode": -1, "output": "", "error": ""}
    print(f"\nrunning in terminal: \n {command}")
    p = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = p.communicate()
    p.wait()
    output, error = output.decode(), error.decode()
    res["returncode"] = p.returncode
    res["output"] = output
    res["error"] = error

    return res 
