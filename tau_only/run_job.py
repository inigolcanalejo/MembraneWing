import subprocess

tau_path = "/home/inigo/software/TAU/taudir_repos.2019.13.08/bin/py_turb1eq/tau.py"
tau_script = "run_tau.py"
tau_input_file = "airfoil_Structured.cntl"
tau_log_file = "log_TAU.out"
p = subprocess.Popen(
    ["python", tau_path, tau_script, tau_input_file, tau_log_file])
p.communicate()