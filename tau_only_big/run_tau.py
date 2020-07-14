# import modules for TAU
import PyPara
import PyPrep
import PySolv
import shutil, sys
sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
from tau_python import tau_mpi_rank
from tau_python import tau_mpi_nranks
from tau_python import tau_solver_unsteady_get_physical_time
rank = tau_mpi_rank()

# Definition of the parameter file
para_path='/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only_big/airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Init Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)

Prep.run(write_dualgrid=0,free_primgrid=False)
Solver.init()
Solver.outer_loop()
Solver.output()
Solver.finalize()
Para.free_parameters()
tau("exit")