# -*- coding: utf-8 -*-
#!/usr/bin/python

import time
import shutil
import sys

sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
sys.path.append('/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only')

# import modules for TAU
import PyPara
import PyPrep
import PySolv
from tau_python import *

# Definition of the parameter file
para_path='/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only/airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Init Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)

grid = Para.get_para_value("Primary grid filename") # Primary grid filename
#n_outer = int(Para.get_para_value('Unsteady physical time steps'))

# Prep + DataSet
# Read Parameter file with Para already done
# preprocessing to create dual grid structure
Prep.run(write_dualgrid=0,free_primgrid=False)

# solve
Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1) # flow solver  print "time step check out = %d"% this_step_out
Solver.outer_loop()
# Solver.output()
tau_plt_init_tecplot_params(para_path_mod)
tau_solver_write_output_conditional()
print 'Solve ok'

Solver.finalize()
tau_free_dualgrid()
tau_free_prims()
Para.free_parameters()