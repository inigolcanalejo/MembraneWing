# -*- coding: utf-8 -*-
#!/usr/bin/python

import time
import shutil
import sys
import array
import string
import scipy
import fileinput
import sys
import subprocess
import time
import numpy as np
import matplotlib.pylab as plt
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os

sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
sys.path.append('/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only')

# import modules for TAU
import PyPara
import PyDeform
import PyPrep
import PySurfDeflect
import PySolv
import PyDataSet
from tau_python import tau_mpi_rank
from tau_python import tau_mpi_nranks
from tau_python import *
from tau_python import tau_msg
rank = tau_mpi_rank()

# Definition of the parameter file
para_path='/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only/airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

#-------------------------------------------------------------------------------
# Init Tau python classes
#-------------------------------------------------------------------------------
Para = PyPara.Parafile(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)
DataSetList = PyDataSet.DataSetList()

grid = Para.get_para_value("Primary grid filename") # Primary grid filename
#n_outer = int(Para.get_para_value('Unsteady physical time steps'))
surfaces = ["MEMBRANE"]

#------------------------------------------------------------------
# Prep + DataSet
#-------------------------------------------------------------------
Prep.run(write_dualgrid=1,free_primgrid=False) # Read Parameter file with Para already done
                                                    # preprocessing to create dual grid structure

tau_parallel_sync()
#-----------------------------------------------------------------
# solve
#-----------------------------------------------------------------
this_step_out=0
n_outer_in=1 # Number of iterations for one coupling iteration
n_outer_out=1
fluidIter=50
for this_step_out in range(n_outer_out):
  this_step=0

  Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1) # flow solver  print "time step check out = %d"% this_step_out
  Solver.outer_loop()
  Solver.output()
  tau_plt_init_tecplot_params(para_path_mod)
  tau_solver_write_output_conditional()
  tau_parallel_sync()
  if tau_mpi_rank() == 0:
    print 'Solve ok'
  tau_parallel_sync()
  time.sleep(20)

  DataSetList.write_output()

  tau_parallel_sync()
  Solver.finalize()
  tau_free_dualgrid()
  tau_free_prims()
  Para.free_parameters()

DataSetList.free_data()

'''
      Functions.meshDeformation(NodesNr,nodes,dispTau,dispTauOld)
      print "deformationstart"

      for i in xrange(0,3*NodesNr):
          dispTauOld[i]=dispTau[i]

      print "afterDeformation"


      Solver.output()
      Solver.finalize()
      tau_free_dualgrid()
      tau_free_prims()
      Para.free_parameters()

DataSetList.free_data()
'''


