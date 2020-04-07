# -*- coding: utf-8 -*-
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
import numpy as np #numpy ist ein Grundpaket enthält bspw. eine leistungsfähiges Arrayobjekt
from distutils.version import StrictVersion
import scipy
from scipy.io import netcdf
import math as m
import os
import glob 
import re 

sys.path.append('/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq')

import PyPara
import PyDeform
import PyPrep
import PySurfDeflect
import PySolv
import PyDataSet 

from tau_python import tau_msg
from tau_python import *

para_path= 'TBD'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

Para = PyPara.Parafile(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)

def findFile(path):
    List_file = glob.glob(path)
    return List_file

def findFname(List_file,path,this_step_out):    			
    this_step_out +=1
    for file in List_file:
        if file.startswith('%s'%path + '/Outputs/airfoilSol.MEMBRANE_i=' + '%s'%this_step_out):   #'%e'+'%d',%path %this_step_out
            print file
            return file

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
              a = line
          line = f.readline()
    return it

def readPressure(fname,it,error):
    with open(fname,'r') as f:
      header1 = f.readline()  #liest document linie für linie durch 
      header2 = f.readline()  #5 readline- Annahme fünf spalten 
      header3 = f.readline()
      header4 = f.readline()
      header5 = f.readline()

      d=[int(s) for s in re.findall(r'\b\d+\b', header4)]  #find all sucht muster im document 
      NodesNr = d[0]
      print 'NodesNr', NodesNr
      ElemsNr = d[1]
      print 'ElemsNr', ElemsNr

# write X,Y,Z,CP of the document in a vector = liste_number
      liste_number = []
      for i in xrange(it):
        line = f.readline()
        for elem in line.split():
          liste_number.append(float(elem))
      
# write ElemTable of the document 
      elemTable = np.zeros([ElemsNr,4])
      k = 0
      line = f.readline()
      while line:
        elemTable[k,0]=float(line.split()[0]) 
        elemTable[k,1]=float(line.split()[1])
        elemTable[k,2]=float(line.split()[2])
        elemTable[k,3]=float(line.split()[3])
	k=k+1
        line = f.readline()

# reshape content in X, Y, Z, Cp      
      X=liste_number[0:NodesNr]
      Y=liste_number[NodesNr:2*NodesNr]
      Z=liste_number[2*NodesNr:3*NodesNr]
      CP=liste_number[3*NodesNr:4*NodesNr]
      X=X[0:NodesNr-error]
      Y=Y[0:NodesNr-error]
      Z=Z[0:NodesNr-error]
      CP=CP[0:NodesNr-error]
      P=[x*20*20*0.5*1.2 for x in CP]
      P=np.array(P)

    return NodesNr,ElemsNr,X,Y,Z,P,elemTable,liste_number

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


def calcpCell(ElemsNr,P,X,elemTable):
  pCell = np.zeros(ElemsNr); # cp for interface elements
  with open('xp','w') as f:
 # f = open('xp', 'w')
    for i in xrange(0,ElemsNr): 
      pCell[i] = 0.25* (P[elemTable[i,0]-1] + P[elemTable[i,1]-1] + P[elemTable[i,2]-1] + P[elemTable[i,3]-1]);
      x= 0.25* (X[elemTable[i,0]-1] + X[elemTable[i,1]-1] + X[elemTable[i,2]-1] + X[elemTable[i,3]-1]);
      f.write('%d\t%f\t%f\n'%(i,x,pCell[i]))
    
  return pCell
 
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

def calcFluidForceVector(ElemsNr,elemTable,NodesNr,pCell,area,normal,fNumber):
  forcesTauNP = np.zeros(NodesNr*3)
  for i in xrange(0,ElemsNr):
    #p= cpCell[i] * q
    p=pCell[i]
    #p=200
    Fx = p * area[i] * normal[i,0]
    Fy = p * area[i] * normal[i,1]
    Fz = p * area[i] * normal[i,2]
    ####Fehler für alle i
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

  
def meshDeformation(NodesNr,nodes,dispTau,dispTauOld,error):

  disp=np.zeros([NodesNr,3])#NodesNr
  for i in xrange(0,NodesNr):#NodesNr
    disp[i,0]=1*(dispTau[3*i+0]-dispTauOld[3*i+0])
    disp[i,1]=1*(dispTau[3*i+1]-dispTauOld[3*i+1])
    disp[i,2]=1*(dispTau[3*i+2]-dispTauOld[3*i+2])
     
  ids, coordinates = PySurfDeflect.read_tau_grid(Para)
  print 'richtig geil'
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
      #print"%d %d %f %f %f"%(i,k,xi,xk,dist2)
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

'''
def meshDeformation(NodesNr,nodes,dispTau,dispTauOld):

  disp=np.zeros([NodesNr,3])
  for i in xrange(0,NodesNr):
    disp[i,0]=1*(dispTau[3*i+0]-dispTauOld[3*i+0])
    disp[i,1]=1*(dispTau[3*i+1]-dispTauOld[3*i+1])
    disp[i,2]=1*(dispTau[3*i+2]-dispTauOld[3*i+2])
    #disp[i,0]=0
    #disp[i,1]=0
    #disp[i,2]=0.5*0.005*(-(nodes[3*i+0]-50)*(nodes[3*i+0]-50)+2500)  
    
  ids, coordinates = PySurfDeflect.read_tau_grid(Para)

  coords=np.zeros([NodesNr,6])
  for i in xrange(0,NodesNr):
    coords[i,0]=coordinates[0,i]
    coords[i,1]=coordinates[1,i]
    coords[i,2]=coordinates[2,i]  

  globalID = np.zeros(NodesNr)  
  for i in xrange(0,NodesNr):
    xi = coords[i,0]
    yi = coords[i,1]
    zi = coords[i,2]
    for k in xrange(0,NodesNr):
      xk = nodes[3*k+0]
      yk = nodes[3*k+1]
      zk = nodes[3*k+2]
      dist2 = (xi-xk)*(xi-xk) + (yi-yk)*(yi-yk) + (zi-zk)*(zi-zk)
      #print"%d %d %f %f %f"%(i,k,xi,xk,dist2)elemTable.astype(int)
      if dist2<0.00001:
	#K=k    
	globalID[i]=k
	globalID = globalID.astype(int)
    #print "%d found %d" % (i,K)
    coords[i,3] = disp[globalID[i],0]
    coords[i,4] = disp[globalID[i],1]
    coords[i,5] = disp[globalID[i],2]
  # print"check2 %d %f %f"% (i, coords[i,0], coords[i,5])
    
  #ids = numpy.matrix(ids)
  #ids = ids.T
  #coords = numpy.matrix(coords)

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

  PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])
  Deform.run(read_primgrid=1, write_primgrid=1, read_deformation=0, field_io=1)
'''