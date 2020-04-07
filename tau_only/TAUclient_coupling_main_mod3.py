# -*- coding: utf-8 -*-
#!/usr/bin/python
#*******************************************************************************
# $Revision: $
# $Date: $
# $Author: $
#******************************************************************************/
#
#******************************************************************************/
#*******************************************************************************
# import modules & load objects used below
#*******************************************************************************

import time
import shutil
import sys
import array
import string
import scipy
import fileinput
import sys

sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
sys.path.append('TBD')
sys.path.append('/home/inigo/software/pyEmpire/EMPIRE-Core/Empire')

import Functions_mod
Functions = Functions_mod
#from Functions import *
from EMPIRE import *
import time
import numpy as np
import matplotlib.pylab as plt

print "1 test"
initEMPIRE('TAUclient')

isDual=0
oppSurfNormal=0
enforceConsistency=1

print "2 test"

ping.recvMesh('myMeshCarat','ping');
caratMeshObj=ping.mesh('myMeshCarat')

print "3 test"

numElemsCaratMesh = caratMeshObj.numElems.value
numNodesCaratMesh = caratMeshObj.numNodes.value
nodesCaratMesh = caratMeshObj.nodes

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

from distutils.version import StrictVersion
import scipy

from scipy.io import netcdf
import math as m
import os
from tau_python import tau_msg

#-------------------------------------------------------------------------------
# Definition of the parameter file
para_path='TBD'
ioname = 'test2'

#-------------------------------------------------------------------------------
# Don't modify original parafile. That is optional.
#-------------------------------------------------------------------------------
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

# for now I need to read in certain data - sigma (surface deflection)
deformfile = Para.get_para_value('RBF basis coordinates and deflections filename')


#if StrictVersion(scipy.__version__).__cmp__('0.8.0') < 0:
 #   raise ImportError('At least SciPy version \'0.8.0\' is needed!')

#------------------------------------------------------------------
# Prep + DataSet
#-------------------------------------------------------------------
Prep.run(write_dualgrid=False, free_primgrid=False) # Read Parameter file with Para already done
                                                    # preprocessing to create dual grid structure

DS = PyDataSet.DataSet(dataset_identifier = "name", output_functions = "surface",dataset_type ="surface",surf_def  = "name",surf_zone_list  = surfaces)
for i in range(0, len(surfaces)):
    DS.define_output(output_name       = surfaces[i],
			output_period     = 1,
                        output_variables  = ["cp"],
                        output_gather     = 0)

DataSetList.store_dataset(DS)

#-----------------------------------------------------------------
# solve
#-----------------------------------------------------------------
this_step_out=0
n_outer_in=1 # Number of iterations for one coupling iteration
n_outer_out=2
#outputInterval = 5
fluidIter=TBD
for this_step_out in range(n_outer_out):
  
  if this_step_out==1:
    search_text = 'Primary grid filename: TBD'
    replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)
    for line in fileinput.input('TBD',inplace=1):
      line = line.replace(search_text, replace_text)
      sys.stdout.write(line)
    Para = PyPara.Parafile(para_path_mod)
    Prep = PyPrep.Preprocessing(para_path_mod)     
    grid = Para.get_para_value("Primary grid filename") # Primary grid filename
    Prep.run(write_dualgrid=False, free_primgrid=False) # Read Parameter file with Para already done
    print this_step_out
    print "step1"

  if this_step_out>1:
    search_text = 'Primary grid filename: TBD.def.' + str(this_step_out-1)
    replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)
    for line in fileinput.input('TBD',inplace=1):
      line = line.replace(search_text, replace_text)
      sys.stdout.write(line)
    Para = PyPara.Parafile(para_path_mod)
    Prep = PyPrep.Preprocessing(para_path_mod)
    grid = Para.get_para_value("Primary grid filename") # Primary grid filename
    Prep.run(write_dualgrid=False, free_primgrid=False) # Read Parameter file with Para already done
    print grid

  Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1) # flow solver  print "time step check out = %d"% this_step_out
  print "time step check out = %d"% this_step_out
  
  Solver.outer_loop() # flow solver
  tau_plt_init_tecplot_params(para_path)  
  tau_solver_write_output_conditional()
  DataSetList.write_output()
   
  List_file = Functions.findFile('/home/inigo/simulations/membranProjekt/Ploetz_FSI_instationaer_U20_AOA6/Outputs/*.plt')
  print List_file
  print this_step_out
  fname = Functions.findFname(List_file,this_step_out,'/home/inigo/simulations/membranProjekt/Ploetz_FSI_instationaer_U20_AOA6/Outputs/')
  print fname
  # Read pressure distribution from file
  NodesNr,ElemsNr,X,Y,Z,P,elemTable=Functions.readPressure(this_step_out,fname)
  elemTable=elemTable.astype(int)
    #temporary write to file
  f = open('xpNode', 'w')
  for i in xrange(0,NodesNr):
    f.write('%d\t%f\t%f\t%f\t%f\n'%(i,X[i],Y[i],Z[i],P[i]))
  f.close()  
  f = open('ElemTable', 'w')
  for i in xrange(0,ElemsNr):
    f.write('%d\t%f\t%f\t%f\t%f\n'%(i,elemTable[i,0],elemTable[i,1],elemTable[i,2],elemTable[i,3]))
  f.close()
 ##P[i]=0
    
      
    
  # Interface Mesh_Fluid side
  nodes,nodesID,elems,numNodesPerElem=Functions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)
 
  if (this_step_out==0):  
    TAUclient.setMesh('myMeshTau', NodesNr, ElemsNr, nodes, nodesID, numNodesPerElem, elems,'TAUclient')
 
  # calculating cp at the center of each interface element
  pCell=Functions.calcpCell(ElemsNr,P,X,elemTable)

  # calculating element area and normal vector
  area,normal = Functions.calcAreaNormal(ElemsNr,elemTable,X,Y,Z,fluidIter*(this_step_out+1))

  # calculating the force vector
  forcesTauNP = Functions.calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fluidIter*(this_step_out+1))
  ####################################
  #       forece Mapping		 ###
  ####################################
  
  forcesCaratNP=np.zeros(3*numNodesCaratMesh)
  ping.setDataField('forcesCarat',forcesCaratNP);
  fC = ping.dataField('forcesCarat')  
  TAUclient.setDataField('forcesTau',forcesTauNP)  
  fT=TAUclient.dataField('forcesTau')
  mapperName1 = 'mapperForce'+str(this_step_out+1)
  
  mapper = FEMapper.FEMapper(mapperName1,'Mortar',TAUclient,'myMeshTau',ping,'myMeshCarat', isDual, oppSurfNormal, enforceConsistency)
    
  mapper.doConsistentMapping(TAUclient,fT, ping, fC)
  
  ping.sendDataField('forcesCarat','ping')
  print "forceMappingDone"
  ####################################
  #       displacement Mapping	 ###
  ####################################
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
  if(this_step_out==0):
    dispTauOld=np.zeros(3*NodesNr)
  
 
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



