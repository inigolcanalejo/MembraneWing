# -*- coding: utf-8 -*-
#*******************************************************************************
# $Revision: $
# $Date: $
# $Author: $
#******************************************************************************/
#
#******************************************************************************/
#*******************************************************************************
# import modules Python, TAU, and initialize the python classes 
#*******************************************************************************
import time
import shutil
import sys
import array
import string
import scipy
import numpy as np 
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os
import glob 
import re 

sys.path.append('/work/piquee/Softwares/TAU/TAU_OPENMPI/taubin_release.2016.1.0_OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')

import PyPara
import PyDeform
import PyPrep
import PySurfDeflect

from tau_python import tau_msg
from tau_python import *

# Need the parameter file to initialize the python classes
para_path= '/work/piquee/coupling_TAU_Carat/SuperMUC_Script/airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
# shutil.copy(para_path, para_path_mod) # delete if working

#*******************************************************************************
# definition of functions
#*******************************************************************************
# Find all the data in the file 'path'
def findFile(path):
    List_file = glob.glob(path)
    return List_file

# Find the data with the name 'word' in the file 'path'
def findFname_0(List_file,path, word):
    for file in List_file:
        if file == '%s'%path + '%s'%word:
            print file
            return file

# Find data with the name 'word' and 'this_step_out' in the file 'path'
def findFname(List_file,path,this_step_out, word):
    this_step_out += 1
    for file in List_file:
        if file.startswith('%s'%path + '%s'%word + '%s'%this_step_out):
            print file
            return file

# Read 'fname'-solution data- and give the number of lines inside - The output parameters are saved in block
# format. Thefore it is not equal to the number of nodes.
def numberLines(fname):
    with open(fname,'r') as f:
      header1 = f.readline()  #liest document linie für linie durch 
      header2 = f.readline()  #5 readline- Annahme fünf spalten 
      header3 = f.readline()
      header4 = f.readline()
      header5 = f.readline()
      it = 0
      line = f.readline()
      while line:
          if 'E' in line:
              it = it+1
              line = f.readline()
          else:
	      break
    return it


# Modify the 'Tautoplt.cntl' from 'Tautoplt_initial.cntl' file to write the solutions, or the mesh in '.dat'
def readTautoplt(fname_mod,fname_o,grid_path,sol_path,path):
    fs = open(fname_o,'r+')
    fd = open(fname_mod,'w')
    line = fs.readline()
    while line:
      if 'Primary grid filename:' in line:
	  line = 'Primary grid filename:' + grid_path + ' \n'
	  fd.write(line) 
          print line
	  line = 'Boundary mapping filename:' + path + '/airfoil_Structured.cntl' + ' \n'
	  fd.write(line)  
          print line
          line = 'Restart-data prefix:' + sol_path + ' \n'
          fd.write(line)
          print line
          line = fs.readline()
	  line = fs.readline()
	  line = fs.readline()
	  print line
      else:
          fd.write(line)
          line = fs.readline()
    fd.close()
    fs.close()
    return fname_mod, sol_path, grid_path


# Find how many nodes and elements are in the part 'membrane' in the TAU mesh ONLY for the MEMBRANE part
def readElement_Mesh(fname, word_Membrane):
    with open(fname,'r') as f:
      line = f.readline()
      while line:
	if 'ZONE T=' + word_Membrane in line:
	  header2 = f.readline()  #5 readline- Annahme fünf spalten 
	  d=[int(s) for s in re.findall(r'\b\d+\b', header2)]  #find all sucht muster im document 
	  NodesNr_Mesh = d[0]
	  print 'NodesNr_Mesh', NodesNr_Mesh
	  ElemsNr_Mesh = d[1]
	  print 'ElemsNr_Mesh', ElemsNr_Mesh
          break
	else:
	  line = f.readline()
    return NodesNr_Mesh, ElemsNr_Mesh


# Write the elements in the TAU mesh in a Table
def readElement(fname,elemsNr_Mesh):
    with open(fname,'r') as f:
      header1 = f.readline()  #liest document linie für linie durch 
      header2 = f.readline()  #5 readline- Annahme fünf spalten 
      header3 = f.readline()
      header4 = f.readline()
      header5 = f.readline()
      line = f.readline()  #liest document linie für linie durch 
      elemTable_Mesh = np.zeros([elemsNr_Mesh,4])
      k = 0
      while line:
	if 'E' in line:
	    line = f.readline()
        else:
	    for i in range(elemsNr_Mesh):
              for elem in line.split():
	        elemTable_Mesh[k,0]=float(line.split()[0]) 
                elemTable_Mesh[k,1]=float(line.split()[1])
                elemTable_Mesh[k,2]=float(line.split()[2])
                elemTable_Mesh[k,3]=float(line.split()[3])
              k=k+1
	      line = f.readline()
            break 
    return elemTable_Mesh

# Find how many nodes and elements are in the TAU mesh
def readElement_Sol(fname):
    with open(fname,'r') as f:
      header1 = f.readline()  #liest document linie für linie durch 
      header2 = f.readline()  #5 readline- Annahme fünf spalten 
      header3 = f.readline()
      header4 = f.readline()
      header5 = f.readline()
      d=[int(s) for s in re.findall(r'\b\d+\b', header4)]  #find all sucht muster im document 
      NodesNr_Sol = d[0]
      print 'NodesNr_Sol', NodesNr_Sol
      ElemsNr_Sol = d[1]
      print 'ElemsNr_Sol', ElemsNr_Sol
    return NodesNr_Sol, ElemsNr_Sol


# Read Cp from the solution file and calculate 'Pressure' on the nodes of TAU Mesh
def readPressure(fname,it,error):
    with open(fname,'r') as f:
      header1 = f.readline()  #liest document linie für linie durch 
      header2 = f.readline()  #5 readline- Annahme fünf spalten 
      header3 = f.readline()
      header4 = f.readline()
      header5 = f.readline()

      d=[int(s) for s in re.findall(r'\b\d+\b', header4)]  #find all sucht muster im document 
      NodesNr = d[0]
      ElemsNr = d[1]

      # write X,Y,Z,CP of the document in a vector = liste_number
      liste_number = []
      for i in xrange(it):
        line = f.readline()
        for elem in line.split():
          liste_number.append(float(elem))
	
      # write ElemTable of the document 
      elemTable_Sol = np.zeros([ElemsNr,4])
      k = 0
      line = f.readline()
      while line:
        elemTable_Sol[k,0]=float(line.split()[0]) 
        elemTable_Sol[k,1]=float(line.split()[1])
        elemTable_Sol[k,2]=float(line.split()[2])
        elemTable_Sol[k,3]=float(line.split()[3])
	k=k+1
        line = f.readline()
	
      # reshape content in X, Y, Z, Cp
      X=liste_number[0:NodesNr]
      Y=liste_number[NodesNr:2*NodesNr]
      Z=liste_number[2*NodesNr:3*NodesNr]
      CP=liste_number[5*NodesNr:6*NodesNr]
      X=X[0:NodesNr-error]
      Y=Y[0:NodesNr-error]
      Z=Z[0:NodesNr-error]
      CP=CP[0:NodesNr-error]
      P=[x*20*20*0.5*1.2 for x in CP]
      P=np.array(P)

    return NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable_Sol,liste_number


# 
def interfaceMeshFluid(NodesNr,ElemsNr,elemTable,X,Y,Z):
  nodes=np.zeros(NodesNr*3) # array to store the coordinates of the nodes in the fluid mesh: x1,y1,z1,x2,y2,z2,...
  nodesID=np.zeros(NodesNr) # array to store the IDs of the nodes in the fluid mesh: IDnode1, IDnode2,...
  elems =np.zeros(4*ElemsNr)# array to store the element table
  numNodesPerElem=np.zeros(ElemsNr)

  for i in xrange(0,NodesNr):
    nodes[3*i+0] = X[i]
    nodes[3*i+1] = Y[i]
    nodes[3*i+2] = Z[i]
    #print "nodes1 %f %f %f"%(nodes[3*i+0],nodes[3*i+1],nodes[3*i+2])

  for i in xrange(0,NodesNr):
    nodesID[i] = i+1

  for i in xrange(0,ElemsNr): 
    elems[i*4+0]=elemTable[i,0]     
    elems[i*4+1]=elemTable[i,1]
    elems[i*4+2]=elemTable[i,2]
    elems[i*4+3]=elemTable[i,3]
    numNodesPerElem[i]=4;  

  return nodes,nodesID,elems,numNodesPerElem


# Calculate the Pressure on the elements from the pressure on the nodes
def calcpCell(ElemsNr,P,X,elemTable):
  pCell = np.zeros(ElemsNr); # cp for interface elements
  print 'len(elemTable) = ', len(elemTable)
  with open('xp','w') as f:
    for i in xrange(0,ElemsNr): 
      pCell[i] = 0.25* (P[elemTable[i,0]-1] + P[elemTable[i,1]-1] + P[elemTable[i,2]-1] + P[elemTable[i,3]-1]);
      x= 0.25* (X[elemTable[i,0]-1] + X[elemTable[i,1]-1] + X[elemTable[i,2]-1] + X[elemTable[i,3]-1]);
      f.write('%d\t%f\t%f\n'%(i,x,pCell[i]))
    
  return pCell


# Calculate the area and the normal of the area for each cell - element of TAU Mesh
def calcAreaNormal(ElemsNr,elemTable,X,Y,Z,fNumber):
  area = np.zeros(ElemsNr); # area of each interface element
  normal = np.zeros([ElemsNr,3]) # element normal vector
  A = np.zeros(3);
  B = np.zeros(3);
  for i in xrange(0,ElemsNr):
    A[0] = X[elemTable[i,2]-1] - X[elemTable[i,0]-1]
    A[1] = Y[elemTable[i,2]-1] - Y[elemTable[i,0]-1]
    A[2] = Z[elemTable[i,2]-1] - Z[elemTable[i,0]-1]
    B[0] = X[elemTable[i,3]-1] - X[elemTable[i,1]-1]
    B[1] = Y[elemTable[i,3]-1] - Y[elemTable[i,1]-1]
    B[2] = Z[elemTable[i,3]-1] - Z[elemTable[i,1]-1]
    a=np.cross(A,B)   #np quasi ein matematischer provider
    norm=np.linalg.norm(a)
    a=a/norm
    normal[i,0]=a[0]
    normal[i,1]=a[1]
    normal[i,2]=a[2]
    A[0] = X[elemTable[i,1]-1] - X[elemTable[i,0]-1]
    A[1] = Y[elemTable[i,1]-1] - Y[elemTable[i,0]-1]
    A[2] = Z[elemTable[i,1]-1] - Z[elemTable[i,0]-1]
    B[0] = X[elemTable[i,3]-1] - X[elemTable[i,0]-1]
    B[1] = Y[elemTable[i,3]-1] - Y[elemTable[i,0]-1]
    B[2] = Z[elemTable[i,3]-1] - Z[elemTable[i,0]-1]
    a=np.cross(A,B)
    norm=np.linalg.norm(a)
    area[i] =  norm # hold only for 2d case, MUST be modified for 3d interfaces
  f_name = 'Outputs/Area_'+str(fNumber)+'.dat'
  with open(f_name,'w') as fwrite:
    for i in xrange(0,len(area[:])):
      fwrite.write("%f\n" % (area[i]))
  return area, normal

		    
# Calculate the Vector Force 
def calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fNumber):
  forcesTauNP = np.zeros(NodesNr*3)
  for i in xrange(0,ElemsNr):
    #p= cpCell[i] * q
    p=pCell[i]
    Fx = p * area[i] * normal[i,0]
    Fy = p * area[i] * normal[i,1]
    Fz = p * area[i] * normal[i,2]
    #print 'test Fx, Fy, Fz', Fx, Fy, Fz
    forcesTauNP[3*(elemTable[i,0]-1)+0] += 0.25 * Fx
    forcesTauNP[3*(elemTable[i,0]-1)+1] += 0.25 * Fy
    forcesTauNP[3*(elemTable[i,0]-1)+2] += 0.25 * Fz
    forcesTauNP[3*(elemTable[i,1]-1)+0] += 0.25 * Fx
    forcesTauNP[3*(elemTable[i,1]-1)+1] += 0.25 * Fy
    forcesTauNP[3*(elemTable[i,1]-1)+2] += 0.25 * Fz
    forcesTauNP[3*(elemTable[i,2]-1)+0] += 0.25 * Fx
    forcesTauNP[3*(elemTable[i,2]-1)+1] += 0.25 * Fy
    forcesTauNP[3*(elemTable[i,2]-1)+2] += 0.25 * Fz
    forcesTauNP[3*(elemTable[i,3]-1)+0] += 0.25 * Fx
    forcesTauNP[3*(elemTable[i,3]-1)+1] += 0.25 * Fy
    forcesTauNP[3*(elemTable[i,3]-1)+2] += 0.25 * Fz
  f_name = 'Outputs/ForcesTauNP_'+str(fNumber)+'.dat'
  with open(f_name,'w') as fwrite:
    for i in xrange(0,len(forcesTauNP[:])):
      fwrite.write("%f\n" % (forcesTauNP[i]))
  return forcesTauNP

		    
# Execute the Mesh deformation of TAU  
def meshDeformation(NodesNr,nodes,dispTau,dispTauOld,error):

  disp=np.zeros([NodesNr,3])#NodesNr
  for i in xrange(0,NodesNr):#NodesNr
    disp[i,0]=1*(dispTau[3*i+0]-dispTauOld[3*i+0])
    disp[i,1]=1*(dispTau[3*i+1]-dispTauOld[3*i+1])
    disp[i,2]=1*(dispTau[3*i+2]-dispTauOld[3*i+2])
  Para = PyPara.Parafile(para_path_mod)
  ids, coordinates = PySurfDeflect.read_tau_grid(Para)
  coords=np.zeros([NodesNr-error,6])#NodesNr

  for i in xrange(0,NodesNr-error):
    coords[i,0]=coordinates[0,i]
    coords[i,1]=coordinates[1,i]
    coords[i,2]=coordinates[2,i]  

  globalID = np.zeros(NodesNr-error)  
  for i in xrange(0,NodesNr-error):
    xi = coords[i,0]
    yi = coords[i,1]
    zi = coords[i,2]
    for k in xrange(0,NodesNr-error):
      xk = nodes[3*k+0]
      yk = nodes[3*k+1]
      zk = nodes[3*k+2]
      dist2 = (xi-xk)*(xi-xk) + (yi-yk)*(yi-yk) + (zi-zk)*(zi-zk)

      if dist2<0.00001:
	#K=k    
	globalID[i]=k
	globalID = globalID.astype(int)
    #print "%d found %d" % (i,K)
    coords[i,3] = disp[globalID[i],0]
    coords[i,4] = disp[globalID[i],1]
    coords[i,5] = disp[globalID[i],2]

  fname_new = 'interface_deformfile.nc'
  ncf = netcdf.netcdf_file(fname_new, 'w')
  # define dimensions
  nops = 'no_of_points'
  number_of_points = len(ids[:])
  ncf.createDimension(nops, number_of_points)
  # define variables
  gid = ncf.createVariable('global_id', 'i', (nops,))
  ncx = ncf.createVariable('x', 'd', (nops,))
  ncy = ncf.createVariable('y', 'd', (nops,))
  ncz = ncf.createVariable('z', 'd', (nops,))
  ncdx = ncf.createVariable('dx', 'd', (nops,))
  ncdy = ncf.createVariable('dy', 'd', (nops,))
  ncdz = ncf.createVariable('dz', 'd', (nops,))
  # write data
  gid[:] = ids
  ncx[:] = coords[:,0]
  ncy[:] = coords[:,1]
  ncz[:] = coords[:,2]
  ncdx[:] = coords[:,3]
  ncdy[:] = coords[:,4]
  ncdz[:] = coords[:,5]
  ncf.close()

  return ids, coordinates, globalID, coords