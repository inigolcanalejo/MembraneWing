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
import matplotlib.pylab as plt
sys.path.append('/home/hpc/pr86fi/di73jef/TAU/TAU_INTEL_2016_1_0/bin/py_turb1eq')
import PyPara
import PyDeform
import PyPrep
import PySurfDeflect
import PySolv
import PyDataSet 

para_path= 'TBD'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

Para = PyPara.Parafile(para_path_mod)
Deform = PyDeform.Deformation(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
from distutils.version import StrictVersion
import scipy

from scipy.io import netcdf
import math as m
import os
import glob 
from tau_python import tau_msg
import re 

def findFile(path):
    List_file = glob.glob(path)
    return List_file

def findFname(List_file,this_step_out,path):
    this_step_out +=1
    for file in List_file:
        if file.startswith(path + 'airfoilSol.MEMBRANE_i=%d'%this_step_out):
            print file
            return file

def readPressure(this_step_out,fname):
    f = open(fname,'r')
    header1 = f.readline()  #liest document linie für linie durch 
    header2 = f.readline()  #5 readline- Annahme fünf spalten 
    header3 = f.readline()
    header4 = f.readline()
    header5 = f.readline()

    d=[int(s) for s in re.findall(r'\b\d+\b', header4)]  #findall sucht muster im document 
    NodesNr = d[0]
    ElemsNr = d[1]
    DataLines = int(np.ceil(NodesNr/5.0)) #ceil sucht den kleinsten integer das noch grösser als (x) ist
    X=np.zeros(NodesNr) #gibt einen array aus nullen zurrück mit der grösse NodesNr
    Y=np.zeros(NodesNr) 
    Z=np.zeros(NodesNr)
    CP=np.zeros(NodesNr)
  #reading x coordinates
    for i in xrange(0,DataLines): #xrange ist das selbe wie range mit dem unterschied bei xrange ist die speicherbelegung unabhängig von der grösse 
        header = f.readline()
        final_list = []
        for elem in header.split():
            final_list.append(float(elem))
            for j in xrange(0,len(final_list)):
	        X[5*i+j]=final_list[j]
	
  #readeing the empty line
  #header = f.readline()
  #reading y coordinates
    for i in xrange(0,DataLines):
        header = f.readline()
        final_list = []
        for elem in header.split():
            final_list.append(float(elem))
            for j in xrange(0,len(final_list)):
	        Y[5*i+j]=final_list[j]

  #readeing the empty line
  #header = f.readline()
  #reading z coordinates
    for i in xrange(0,DataLines):
        header = f.readline()
        final_list = []
        for elem in header.split():
            final_list.append(float(elem))
            for j in xrange(0,len(final_list)):
	        Z[5*i+j]=final_list[j]

  #readeing the empty line
  #header = f.readline()
  #reading cp
    for i in xrange(0,DataLines):
        header = f.readline()
        final_list = []
        for elem in header.split():
            final_list.append(float(elem))
            for j in xrange(0,len(final_list)):
	        CP[5*i+j]=final_list[j]
	        P=[x*20*20*0.5*1.2 for x in CP]
                P=np.array(P)
    elemTable = np.zeros([ElemsNr,4])
  #readeing the empty line
  #header = f.readline()
    for i in xrange(0,ElemsNr):
        header = f.readline()
        final_list = []
        for elem in header.split():
            final_list.append(float(elem))
        elemTable[i,0]=final_list[0]
        elemTable[i,1]=final_list[1]
        elemTable[i,2]=final_list[2]
        elemTable[i,3]=final_list[3]  
  
    return NodesNr,ElemsNr,X,Y,Z,P,elemTable




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
  pCell = np.zeros(ElemsNr) # cp for interface elements
  f = open('xp', 'w')
  for i in xrange(0,ElemsNr): 
    pCell[i] = 0.25* (P[int(elemTable[i,0]-1)] + P[int(elemTable[i,1]-1)] + P[int(elemTable[i,2]-1)] + P[int(elemTable[i,3]-1)]);
    x= 0.25* (X[elemTable[i,0]-1] + X[elemTable[i,1]-1] + X[elemTable[i,2]-1] + X[elemTable[i,3]-1]);
    f.write('%d\t%f\t%f\n'%(i,x,pCell[i]))
    #print "pCell\t%d\t%f"%(i,pCell[i])
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
  fwrite = open(f_name,'w')
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
  fwrite = open(f_name,'w')
  for i in xrange(0,len(forcesTauNP[:])):
    fwrite.write("%f\n" % (forcesTauNP[i]))
  return forcesTauNP


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




