import subprocess

tau_path = "/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq/tau.py"
tau_script = "run_tau.py"
tau_input_file = "airfoil_Structured.cntl"
tau_log_file = "log_TAU.out"
p = subprocess.Popen(
    ["python", tau_path, tau_script, tau_input_file, tau_log_file])
p.communicate()