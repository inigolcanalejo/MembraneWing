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
import subprocess
import time
import numpy as np
import matplotlib.pylab as plt
#from mpi4py import MPI
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os

sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')
sys.path.append('TBD')
sys.path.append('/home/inigo/software/pyEmpire/EMPIRE-Core/Empire')

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

from EMPIRE import *

import Functions_mod3
Functions = Functions_mod3

initEMPIRE('TAUclient')

isDual=False
oppSurfNormal1=True
oppSurfNormal2=True
enforceConsistency1=0
enforceConsistency2=0


if (rank == 0):
 ping.recvMesh('myMeshCarat','ping');
 caratMeshObj=ping.mesh('myMeshCarat')

 numElemsCaratMesh = caratMeshObj.numElems.value
 numNodesCaratMesh = caratMeshObj.numNodes.value
 nodesCaratMesh = caratMeshObj.nodes

#-------------------------------------------------------------------------------
# Definition of the parameter file
para_path='TBD'
ioname = 'test2'
TAU_path = '/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/'
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
deformfile = Para.get_para_value('RBF basis coordinates and deflections filename')


#------------------------------------------------------------------
# Prep + DataSet
#-------------------------------------------------------------------
Prep.run(write_dualgrid=1,free_primgrid=False) # Read Parameter file with Para already done
                                                    # preprocessing to create dual grid structure

DS = PyDataSet.DataSet(dataset_identifier = "name", output_functions = "surface",dataset_type ="surface",surf_def  = "name",surf_zone_list  = surfaces)
for i in range(0, len(surfaces)):
    DS.define_output(output_name       = surfaces[i],
			output_period     = 1,
                        output_variables  = ["cp"],
                        output_gather     = 1)

DataSetList.store_dataset(DS)
tau_parallel_sync()
#-----------------------------------------------------------------
# solve
#-----------------------------------------------------------------
this_step_out=0
n_outer_in=1 # Number of iterations for one coupling iteration
n_outer_out=TBD
fluidIter=TBD
for this_step_out in range(n_outer_out):
  this_step=0

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
###

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
  if tau_mpi_rank() == 0:
    print 'Data wrote'
  tau_parallel_sync()

  if tau_mpi_rank() == 0:
    List_file = Functions.findFile('/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only/Outputs/*.plt')
    print List_file
    fname = Functions.findFname(List_file,'/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only',this_step_out)
    print fname
    it = Functions.numberLines(fname)

    # Read pressure distribution from file
    NodesNr,ElemsNr,X,Y,Z,P,elemTable,liste_number=Functions.readPressure(fname,it,0)
    elemTable = elemTable.astype(int)

    with open('xpNode','w') as f:
      for i in xrange(0,NodesNr):
	f.write('%d\t%f\t%f\t%f\t%f\n'%(i,X[i],Y[i],Z[i],P[i]))

    with open('ElemTable', 'w')  as f:
      for i in xrange(0,ElemsNr):
	f.write('%d\t%f\t%f\t%f\t%f\n'%(i,elemTable[i,0],elemTable[i,1],elemTable[i,2],elemTable[i,3]))

    nodes,nodesID,elems,numNodesPerElem=Functions.interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z)


  if (this_step_out==0 and rank == 0):
      TAUclient.setMesh('myMeshTau', NodesNr, ElemsNr, nodes, nodesID, numNodesPerElem, elems, 'TAUclient')
      # calculating cp at the center of each interface element
      pCell=Functions.calcpCell(ElemsNr,P,X,elemTable)
      # calculating element area and normal vector
      area,normal = Functions.calcAreaNormal(ElemsNr,elemTable,X,Y,Z,fluidIter*(this_step_out+1))

    # calculating the force vector
      forcesTauNP = Functions.calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fluidIter*(this_step_out+1))
      print 'calcFluidForceVector'
    ####################################
    #       force Mapping		 ###
    ####################################
  if(rank == 0):
      forcesCaratNP=np.zeros(3*numNodesCaratMesh)
      ping.setDataField('forcesCarat',forcesCaratNP);
      print 'ping.setDataField'
      fC = ping.dataField('forcesCarat')
      TAUclient.setDataField('forcesTau',forcesTauNP)
      fT=TAUclient.dataField('forcesTau')
      mapperName1 = 'mapperForce'+str(this_step_out+1)
      #mapper = FEMapper.FEMapper(mapperName1,'Mortar',TAUclient,'myMeshTau',ping,'myMeshCarat')#, isDual, oppSurfNormal, enforceConsistency)
      mapper = FEMapper.FEMapper(mapperName1,'Mortar',TAUclient,'myMeshTau',ping,'myMeshCarat', isDual, oppSurfNormal1, enforceConsistency1)
      print 'FEMapper.FEMapper'
      mapper.doConsistentMapping(TAUclient,fT, ping, fC)
      print 'mapper.doConsistentMapping'
      ping.sendDataField('forcesCarat','ping')
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
      mapper = FEMapper.FEMapper(mapperName2,'Mortar',ping,'myMeshCarat',TAUclient,'myMeshTau', isDual, oppSurfNormal2, enforceConsistency2)

      mapper.doConsistentMapping(ping,dC, TAUclient, dT)
      print "displacementMappingDone"
      dispTau = TAUclient.dataField('displacementsTau').array
      print 'TAU line 240'
      print 'this_step_out= ', this_step_out
      if(this_step_out==0):
	  print 'in if schleife'#
	  dispTauOld=np.zeros(3*NodesNr)
	  print 'noch in if...'#
      print 'schon draußen'#
      dispTau_transpose = np.transpose(dispTau)
      print 'dispTau =', dispTau_transpose

  #-------------------------------------------------------------------------------
  # Gather Process + Deformation
  #-------------------------------------------------------------------------------
  if tau_mpi_rank() == 0:
    print 'gather'
    subprocess.call(TAU_path + 'gather ' + para_path_mod ,shell=True)
    print 'finish gather'
    print "deformationstart"
    [ids,coordinates,globalID,coords]=Functions.meshDeformation(NodesNr,nodes,dispTau,dispTauOld,0)
    PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])
    print "afterPySurfDeflect"

  tau_parallel_sync()
  Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1)

  if tau_mpi_rank() == 0:
    for i in xrange(0,3*NodesNr):
      dispTauOld[i]=dispTau[i]
    print "afterDeformation"

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


