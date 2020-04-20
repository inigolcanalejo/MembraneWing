# -*- coding: utf-8 -*-
#!/work/piquee/Softwares/Python.2.7/anaconda2/bin/python
#*******************************************************************************
# $Revision: $
# $Date: $
# $Author: $
#******************************************************************************/
#
#******************************************************************************/
#
#
#*******************************************************************************
# import modules for Python
#*******************************************************************************

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
# import matplotlib.pylab as plt
# from mpi4py import MPI
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os

#*******************************************************************************
# path
#*******************************************************************************
sys.path.append('/work/piquee/Softwares/TAU/TAU_OPENMPI/taubin_release.2016.1.0_OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
sys.path.append('/work/piquee/coupling_TAU_Carat/SuperMUC_Script')
sys.path.append('/work/piquee/Softwares/pyEmpire/EMPIRE-Core/Empire')

#*******************************************************************************
# import modules for TAU
#*******************************************************************************

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

#*******************************************************************************
# import modules for Empire and the functions
#*******************************************************************************
from EMPIRE import *

import Functions_mod4
Functions = Functions_mod4
initEMPIRE('TAUclient')

#*******************************************************************************
# Options for CARAT 
#*******************************************************************************
#isDual=0
#oppSurfNormal=1
#enforceConsistency=1

#ping.recvMesh('myMeshCarat','ping'); # 
#caratMeshObj=ping.mesh('myMeshCarat') # 
#numElemsCaratMesh = caratMeshObj.numElems.value
#numNodesCaratMesh = caratMeshObj.numNodes.value
#nodesCaratMesh = caratMeshObj.nodes
#-------------------------------------------------------------------------------
# Definition of the parameter file
#-------------------------------------------------------------------------------

path = '/work/piquee/coupling_TAU_Carat/SuperMUC_Script'
TAU_path = '/work/piquee/Softwares/TAU/TAU_INTEL_2016_1_0/bin'
para_path=path + '/airfoil_Structured.cntl'

#-------------------------------------------------------------------------------
# Don't modify original parafile. That is optional.
#-------------------------------------------------------------------------------
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

#-------------------------------------------------------------------------------
# Init Tau python classes + get the informations necessary for the preprocessing 
#-------------------------------------------------------------------------------
Para = PyPara.Parafile(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)

# Primary grid filename
grid = Para.get_para_value("Primary grid filename") # 
# Surfaces where Cp has to be written 
surfaces = ["MEMBRANE"] # 
# get the name of the file where the deformation of the TAU Mesh are saved
deformfile = Para.get_para_value('RBF basis coordinates and deflections filename') # 

#------------------------------------------------------------------
# Preprocessing 
#-------------------------------------------------------------------
Prep.run(write_dualgrid=1,free_primgrid=False) # Read Parameter file with Para already done
                                               # preprocessing to create dual grid structure

