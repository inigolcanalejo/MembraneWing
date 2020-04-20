# -*- coding: utf-8 -*-
#!/usr/bin/python
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
sys.path.append('/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin/py_turb1eq') 
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
isDual=0
oppSurfNormal=1
enforceConsistency=1

ping.recvMesh('myMeshCarat','ping'); # 
caratMeshObj=ping.meshMap['myMeshCarat'] # 

numElemsCaratMesh = caratMeshObj.numElems.value
numNodesCaratMesh = caratMeshObj.numNodes.value
nodesCaratMesh = caratMeshObj.nodes
#-------------------------------------------------------------------------------
# Definition of the parameter file
#-------------------------------------------------------------------------------

path = '/work/piquee/coupling_TAU_Carat/SuperMUC_Script'
TAU_path = '/work/piquee/Softwares/TAU/TAU_2016.2/2016.2.0/bin' 
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

#------------------------------------------------------------------
# Convert the initial TAU Mesh in '.dat' to find the information necessary  
#-------------------------------------------------------------------
# modify 'Tautoplt.cntl' with the initial TAU mesh name,
fname_tau, dummy_sol_path, grid_path = Functions.readTautoplt(path + '/Tautoplt.cntl', path + '/Tautoplt_initial.cntl',path + '/Mesh/airfoil_Structured_scaliert.grid', '(none)' ,path) #
print 'command = ' + TAU_path + '/tau2plt ' + path + '/Tautoplt.cntl'

# execute the command tau2plt - creation of the TAU mesh in a '.dat' file
subprocess.call(TAU_path + '/tau2plt ' + path + '/Tautoplt.cntl' ,shell=True) #

#definition of the name of the TAU mesh + '.dat'
fname_mesh = path + '/Mesh/airfoil_Structured_scaliert.grid.dat'# 

#find how many nodes and elements are in the part 'membrane' in the TAU mesh
NodesNr_Mesh, ElemsNr_Mesh = Functions.readElement_Mesh(fname_mesh,'"MEMBRANE"') #  

# Find and write in a table the connectivity between the elements of the TAU Mesh 
elemTable_Mesh = Functions.readElement(fname_mesh,ElemsNr_Mesh)# 

# Write the Information in a file 'ElemTable_mesh'
with open('ElemTable_mesh', 'w')  as f:											
  for i in xrange(0,ElemsNr_Mesh):											
    f.write('%d\t%f\t%f\t%f\t%f\n'%(i,elemTable_Mesh[i,0],elemTable_Mesh[i,1],elemTable_Mesh[i,2],elemTable_Mesh[i,3]))

#-----------------------------------------------------------------
# solve
#-----------------------------------------------------------------
n_outer_in=1 # 
n_outer_out=150 # Number of external coupling iterations
fluidIter=200 # Number of TAU iterations 
for this_step_out in range(n_outer_out):
  
  this_step=0
  Para = PyPara.Parafile(para_path_mod)
  Prep = PyPrep.Preprocessing(para_path_mod)     
  grid = Para.get_para_value("Primary grid filename") 
  Prep.run(write_dualgrid=False, free_primgrid=False) 
  print this_step_out
  
  # flow solver 
  print "time step check out = %d"% this_step_out
  Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1) 
  Solver.outer_loop() 
  Solver.output()
  tau_plt_init_tecplot_params(para_path_mod)  
  tau_solver_write_output_conditional()
  tau_parallel_sync()

 # find the solution file name in '/Outputs'
  print 'FindFile fname '
  List_file = Functions.findFile( path +'/Outputs/*')
  print 'List_file=', List_file
  fname = Functions.findFname(List_file,path,this_step_out, '/Outputs/airfoilSol.pval.unsteady_i=')
  fname_s = fname[0:fname.find('pval')]+ 'surface.' + fname[fname.find('pval'):len(fname)]
  print 'fname = ', fname

 # find the mesh file name in '/Mesh'
  if this_step_out == 0:
    print 'FindFile fname mesh '
    List_file_mesh = Functions.findFile(path + '/Mesh/*')
    print 'List_file_mesh=', List_file_mesh
    fname_mesh = Functions.findFname_0(List_file_mesh,path + '/Mesh/','airfoil_Structured_scaliert.grid')
    print 'fname_mesh = ', fname_mesh
  else:
    print 'FindFile fname mesh '
    List_file_mesh = Functions.findFile(path + '/Mesh/*')
    print 'List_file_mesh=', List_file_mesh
    fname_mesh = Functions.findFname(List_file_mesh,path + '/Mesh/',this_step_out-1 ,'airfoil_Structured_scaliert.grid.def.')
    print 'fname_mesh = ', fname_mesh

  # write the solution data of the time step this_step_out in .dat
  print 'Write Data with TAU2plt'
  fnameTAU_mod, sol_path, grid_path = Functions.readTautoplt(path + '/Tautoplt.cntl', path + '/Tautoplt_initial.cntl',fname_mesh,fname,path)
  subprocess.call(TAU_path + '/tau2plt ' + path + '/Tautoplt.cntl' ,shell=True)
  print TAU_path + '/tau2plt ' + path + '/Tautoplt.cntl'
  print 'Data written'
  print 'Number of line in fname = ', fname_s + '.dat'
  it = Functions.numberLines(fname_s + '.dat')
  print 'it =', it


  # Read pressure distribution from file
  NodesNr_Sol,ElemsNr_Sol,X,Y,Z,CP,P,elemTable_Sol,liste_number=Functions.readPressure(fname_s + '.dat',it,0)

  elemTable = elemTable_Sol.astype(int)
  NodesNr = NodesNr_Sol
  ElemsNr = ElemsNr_Sol

  with open('xpNode','w') as f:
    for i in xrange(0,NodesNr):											
      f.write('%d\t%f\t%f\t%f\t%f\n'%(i,X[i],Y[i],Z[i],P[i]))	

  nodes,nodesID,elems,numNodesPerElem=Functions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)
  print 'aaaaaaaaaaaaaaaaaaaaa'
  if (this_step_out==0):  
      TAUclient.setMesh('myMeshTau', NodesNr, ElemsNr, nodes, nodesID, numNodesPerElem, elems, 'TAUclient')
      print 'bbbbbbbbbbbbbbbbbbbbbbb'

  # calculating cp at the center of each interface element    
  pCell=Functions.calcpCell(ElemsNr,P,X,elemTable)
  print 'calcCpCell'
  # calculating element area and normal vector
  area,normal = Functions.calcAreaNormal(ElemsNr,elemTable,X,Y,Z,fluidIter*(this_step_out+1))
  print 'calcAreaNormal'
  # calculating the force vector
  forcesTauNP = Functions.calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fluidIter*(this_step_out+1))
  print 'calcFluidForceVector'
  print 'ccccccccccccccccccc'

  #####################################
  ##       force Mapping		 ###
  #####################################

  ping.addDataFieldToMesh('myMeshCarat','forcesCarat');
  forcesCaratNP=np.zeros(3*numNodesCaratMesh)
  ping.setDataField('forcesCarat',forcesCaratNP);
  print 'ping.setDataField'
  import sys
  sys.stdout.flush()
  fC = ping.dataField('forcesCarat')  
  TAUclient.setDataField('forcesTau',forcesTauNP)  
  fT=TAUclient.dataField('forcesTau')
  mapperName1 = 'mapperForce'+str(this_step_out+1)  
  mapper = FEMapper.FEMapper(mapperName1,'Mortar',TAUclient,'myMeshTau',ping,'myMeshCarat', isDual, oppSurfNormal, enforceConsistency)
  print 'FEMapper.FEMapper'
  mapper.doConsistentMapping(TAUclient,fT, ping, fC)
  print 'mapper.doConsistentMapping'
  ping.sendDataField('forcesCarat','ping')
  
  #####################################
  ##       displacement Mapping	 ###
  #####################################
  ping.recvDataField('displacementsCarat', 'ping')
  dispTemp = ping.dataField('displacementsCarat').array
  dC = ping.dataField('displacementsCarat')
  displacementTauNP=np.zeros(3*NodesNr)
  TAUclient.setDataField('displacementsTau',displacementTauNP);
  dT = TAUclient.dataField('displacementsTau')
  mapperName2 = 'mapperDisp'+str(this_step_out+1)
  mapper = FEMapper.FEMapper(mapperName2,'Mortar',ping,'myMeshCarat',TAUclient,'myMeshTau', isDual, oppSurfNormal, enforceConsistency)
  mapper.doConsistentMapping(ping,dC, TAUclient, dT)
  print "displacementMappingDone"
  dispTau = TAUclient.dataField('displacementsTau').array
  print 'this_step_out= ', this_step_out
  if(this_step_out==0):
    dispTauOld=np.zeros(3*NodesNr)
    dispTau_transpose = np.transpose(dispTau)
    print 'dispTau =', dispTau_transpose

  ##-------------------------------------------------------------------------------
  ##  Deformation
  ##-------------------------------------------------------------------------------	
  if tau_mpi_rank() == 0:
    print "deformationstart"
    [ids,coordinates,globalID,coords]=Functions.meshDeformation(NodesNr,nodes,dispTau,dispTauOld,0)
    PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])
    print "afterPySurfDeflect"


  Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1) 


  for i in xrange(0,3*NodesNr):
    dispTauOld[i]=dispTau[i]
  print "afterDeformation"

  tau_parallel_sync()
  Solver.finalize()
  tau_free_dualgrid()
  tau_free_prims()
  Para.free_parameters()





