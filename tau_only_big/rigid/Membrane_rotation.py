# -*- coding: utf-8 -*-
# ******************************************************************************
# Import modules & load objects used below
# ******************************************************************************
import sys

#sys.path.append('/hppfs/work/pn69ni/di73jef3/Softwares/TAU/taudir_release.2018.1.0_TMC_impi_1.6_netcdf_4_SPMUC_python_likeHW/bin/py_turb1eq')
sys.path.append('/home/inigo/software/TAU/TAU_2016.2/2016.2.0/bin/py_turb1eq')

import PyPara
import PySolv
import PyPrep
import PyCopyCluster
import tau
import tau_python
import numpy
import PyMotionExternalDelegate
from functools import reduce
import scipy.io as sio
import subprocess
import time


#from mpi4py import MPI
import scipy
from datetime import datetime

from tau_python import *
from tau_python import tau_mpi_rank
from tau_python import tau_msg
from tau_python import tau_mpi_nranks
from tau_python import tau_solver_unsteady_get_physical_time
rank = tau_mpi_rank()

#-------------------------------------------------------------------------------
# Create an instance of PyCopyCluster Class and initialize it:
# -> Handles file-transfer of grid and restart files for standard directory
#    structure
#-------------------------------------------------------------------------------
path = '/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only_big/rigid'
para_path=path + '/Membrane_rotation.cntl'
NCores=2
TAU_path='/home/inigo/software/TAU/TAU_2016.2/2016.2.0/bin'
# mpirun_path='/work/piquee/Softwares/openmpi_1.6.4/bin'
# --- Prepare parameters for unsteady simulation ---
nodeNames        = ["MEMBRANE"]

RestartTime = 1.5
nStepsMAX = 1000 # maximal timesteps
deltaT           = 0.005 # unsteady physical time step size
nTimestepsperRUN       = 7 # number of unsteady time steps per supermuc run
nTimesteps = 207 # updating number of maximum stepnumber for run -> nTimestepsperRUN + nStepsStart
nInnerIters      = 300 # inner iterations for dual time stepping
outputPeriod     = 1 # external output period (needed if output period in para_path not working properly)
pitchDeg = 0 # mean pitch angle


# --- Load external excitation files ---
# --- thetaDeg -> pitch angle per timestep in deg
# --- thetaRate -> pitch rate per timestep in deg/s
thetaDeg=numpy.loadtxt(path + '/signal/APRBSDeg_membrane.dat')
thetaRate=numpy.loadtxt(path + '/signal/APRBSRate_membrane.dat')


class MotionStringGenerator(object):
	"""Auxiliary class to generate TAU motion strings for a translatory x-motion
	and/or a pitching oscillation.
	"""

	def __init__(self, deltaT=deltaT, pitchDeg=pitchDeg, thetaDeg=thetaDeg, thetaRate=thetaRate):
		self.deltaT       = deltaT
		self.pitchDeg = pitchDeg
		self.thetaDeg = thetaDeg
		self.thetaRate = thetaRate


	def GetMotionString(self,step):
		self.time = step*self.deltaT
		self.thetaInstant     = self.thetaDeg[step]
		self.pitchFreq    = self.thetaRate[step]

		p     = 0.
		q     = numpy.deg2rad(self.pitchFreq)
		r     = 0.
		phi   = 0.
		theta = numpy.deg2rad(self.pitchDeg) + numpy.deg2rad(self.thetaInstant)
		psi   = 0.
		u     = 0
		v     = 0.
		w     = 0.
		dx    = 0.
		dy    = 0.
		dz    = 0.
		motionString=" ".join(map(str, [p,q,r,phi,theta,psi,u,v,w,dx,dy,dz]))

		return motionString

	def __call__(self, step):
		return self.GetMotionString(step)

def PrintBlockHeader(header):
 	tau_python.tau_msg("\n" + 50 * "*" + \
 					   "\n" + "* %s\n" %header + \
 					   50*"*" + "\n")

#--- Instanciate required modules ---
MyTauMotionDelegate = PyMotionExternalDelegate.MotionExternalDelegate(nodeNames)

MyMotionStringGenerator = \
			MotionStringGenerator(deltaT, pitchDeg, thetaDeg, thetaRate)

Para   = PyPara.Parafile(para_path)
Prep   = PyPrep.Preprocessing(para_path)
Solver = PySolv.Solver(para_path, delegate=MyTauMotionDelegate)

#------------------------------------------------------------------
# Partitioning
#-------------------------------------------------------------------
#if (rank == 0):
#  print 'Partitioning start', str(datetime.now())
#  print mpirun_path + '/mpirun -np ' + '%s '%NCores + TAU_path + '/ptau3d.subgrids ' + para_path + ' out.partitioning use_mpi'
#  subprocess.call(mpirun_path + '/mpirun -np ' + '%s '%NCores + TAU_path + '/ptau3d.subgrids ' + para_path + ' out.partitioning use_mpi',shell=True)
#  subprocess.call('sleep 5' ,shell=True)
#  print 'Partitioning end'
#tau_parallel_sync()


# --- Prepare parameters for unsteady simulation ---
Para.update({"Unsteady physical time step size" : deltaT,
	     "Unsteady physical time steps": 1,
	     "Unsteady inner iterations per time step": nInnerIters,
	     "Unsteady allow external control over progress in time (0/1)": 1,
             "Unsteady show pseudo time steps (0/1)": 1,
	     "Unsteady enable external control over progress in time (0/1)": 1})



deltaT = float(Para.get_para_value("Unsteady physical time step size"))
nStepsStart = int(float(Para.get_para_value('Finished time steps','single_key')))

if (rank == 0):
  PrintBlockHeader("Restart from stepnumber: %d" %(nStepsStart))
tau_parallel_sync()

# ******************************************************************************
# Start computation
# ******************************************************************************
if nStepsStart < nStepsMAX:
	nTimesteps = nStepsStart+nTimestepsperRUN

# ----------------------------------------------------------------------
	# --- Init motion stack ---
	motionString = MyMotionStringGenerator(0)
	#MyTauMotionDelegate.UpdateMotion("MEMBRANE", motionString)
	MyTauMotionDelegate.PushMotion()
	if (rank == 0):
		PrintBlockHeader("Inital Motionstring: %s" %(motionString))
	tau_parallel_sync()
# ----------------------------------------------------------------------
	# --- Perform preprocessing ---
	if (rank == 0):
		PrintBlockHeader("PREPROCESSOR")
	tau_parallel_sync()
	Prep.run(write_dualgrid=0, free_primgrid=1)

# ----------------------------------------------------------------------
	# --- Execute solver for a couple of unsteady time steps ---
	if (rank == 0):
		PrintBlockHeader("SOLVER")
	tau_parallel_sync()
	Solver.init()
# ----------------------------------------------------------------------
	for iTimestep in range(nStepsStart+1,nTimesteps):
		time = iTimestep * deltaT
		if (rank == 0):
			PrintBlockHeader("time step %s" %(str(iTimestep)))
			PrintBlockHeader("time %s" %(str(time)))
		tau_parallel_sync()
		motionString = MyMotionStringGenerator(iTimestep - nStepsStart)
		MyTauMotionDelegate.UpdateMotion("MEMBRANE", motionString)
		MyTauMotionDelegate.InitExchange()
		if (rank == 0):
			PrintBlockHeader("time step %s" %(str(iTimestep)))
		tau_parallel_sync()

		Solver.outer_loop()

		if iTimestep%outputPeriod == 0:
			tau_python.tau_solver_output_field()
			tau_python.tau_solver_output_surface()
# ----------------------------------------------------------------------
	# --- Shut down everything ---
	Solver.output()
	ActualTime = (tau_solver_unsteady_get_physical_time()-RestartTime)
	stepsDONE = ActualTime/deltaT
	Para.update({"Finished time steps":stepsDONE})
	Solver.finalize()

# ------------------------------------------------------------------------------
# Exit python
# ------------------------------------------------------------------------------
