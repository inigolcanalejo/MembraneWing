# -*- coding: utf-8 -*-
import itertools
def findFileName0(list_of_interface_file_paths,working_path,word): 
    for file in list_of_interface_file_paths:
        if file.startswith('%s'%working_path + '%s'%word): ####### i would like to make it general #######
            print file
            return file

def findFileName(list_of_interface_file_paths,working_path,word,this_step_out): 
    this_step_out +=1
    for file in list_of_interface_file_paths:
        if file.startswith('%s'%working_path + '%s'%word + '%s'%this_step_out): ####### i would like to make it general #######
            print file
            return file

def findInterfaceFileNumberOfLines(fname):
    with open(fname,'r') as f:
        it = 0
        pattern = {'E+', 'E-'}
        for line in f:
            if 'E+' in line or 'E-' in line:
                it = it+1
    return it

def PrintBlockHeader(header):
 	tau_python.tau_msg("\n" + 50 * "*" + "\n" + "* %s\n" %header + 50*"*" + "\n")

def readTautoplt(fname_mod, fname_o, mesh_iteration, interface_file_name, para_path_mod):
    fs = open(fname_o,'r+')
    fd = open(fname_mod,'w')
    line = fs.readline()
    while line:
        if 'Primary grid filename:' in line:
	    line = 'Primary grid filename:' + mesh_iteration + ' \n'
	    fd.write(line) 
            print line
            line = fs.readline()
        if 'Boundary mapping filename:' in line:
	    line = 'Boundary mapping filename:' + para_path_mod + ' \n'
	    fd.write(line)  
            print line
            line = fs.readline()
        if 'Restart-data prefix:' in line:
            line = 'Restart-data prefix:' + interface_file_name + ' \n'
            fd.write(line)
            print line
            line = fs.readline()
        else:
            line = fs.readline()
            fd.write(line)
    fd.close()
    fs.close()


# Read Cp from the solution file and calculate 'Pressure' on the nodes of TAU Mesh
def readPressure(interface_file_name,interface_file_number_of_lines,error,velocity):
    with open(fname,'r') as f:
        header1 = f.readline()  #liest document linie für linie durch 
        header2 = f.readline() 
        print "Careful ---- headers 2 = ", header2 #5 readline- Annahme fünf spalten 
        header2_split = header2.split()
        pos_X = header2_split.index('"x"')
        pos_Y = header2_split.index('"y"')
        pos_Z = header2_split.index('"z"')
        pos_Cp = header2_split.index('"cp"')
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
        X=liste_number[(pos_X-2)*NodesNr:(pos_X-2+1)*NodesNr]
        Y=liste_number[(pos_Y-2)*NodesNr:(pos_Y-2+1)*NodesNr]
        Z=liste_number[(pos_Z-2)*NodesNr:(pos_Z-2+1)*NodesNr]
        CP=liste_number[(pos_Cp-2)*NodesNr:(pos_Cp-2+1)*NodesNr]
        X=X[0:NodesNr-error]
        Y=Y[0:NodesNr-error]
        Z=Z[0:NodesNr-error]
        CP=CP[0:NodesNr-error]
        P=[x*velocity*velocity*0.5*1.2 for x in CP]
        P=np.array(P)

    return NodesNr,ElemsNr,X,Y,Z,CP,P,elemTable_Sol,liste_number
