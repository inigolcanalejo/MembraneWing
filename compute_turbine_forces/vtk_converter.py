output_filename = 'airfoilSol.surface.pval.5010.MEMBRANE_UP.dat'
vtk_output_filename = 'airfoilSol.surface.pval.5010.MEMBRANE_UP.vtk'

import re, sys, json, os
import numpy as np

working_path = os.getcwd() + '/'
with open(working_path + 'input/tau_settings.json') as json_file:
    tau_settings = json.load(json_file)

start_step = tau_settings["start_step"]
tau_path = tau_settings["tau_path"]
sys.path.append(tau_settings["kratos_path"])
sys.path.append(tau_path + "py_turb1eq/")

# tau_functions can only be imported after appending kratos' path
import tau_functions_Steady as TauFunctions

# Read output file
position_info, mesh_info, nodal_data, elem_connectivities = TauFunctions.ReadInterfaceFile(output_filename)

# Read mesh info
NodesNr = mesh_info[0]
print("NodesNr = ", NodesNr)

X, Y, Z = TauFunctions.SaveCoordinatesList(nodal_data, position_info, NodesNr)

velocity = 41.09
P = TauFunctions.SavePressure(nodal_data, position_info, NodesNr, velocity)
cp_position = position_info.index('"cp"') - 2
CP = nodal_data[(cp_position)*NodesNr:(cp_position+1)*NodesNr]

# calculating the force vector
fluid_forces = TauFunctions.CalculateNodalFluidForces(X, Y, Z, P, elem_connectivities)

cell_normals = np.zeros(int(3*len(elem_connectivities)/4))
cell_areas = np.zeros(int(len(elem_connectivities)/4))

point_normals = np.zeros(3*len(X))

for cell in range(int(len(elem_connectivities)/4)):
    # Get the node ids of the cell
    node_ids = TauFunctions.GetCellNodeIds(elem_connectivities, cell)
    normal = TauFunctions.CalculateCellNormal(X, Y, Z, node_ids)
    cell_areas[cell] = TauFunctions.CalculateCellArea(X, Y, Z, node_ids)
    for component in range(3):
        cell_normals[3*cell+component] = normal[component]
    for node in range(4):
        for component in range(3):
            point_normals[3*node_ids[node]+component] += 0.25 * normal[component]

# Write vtk file
with open(vtk_output_filename, 'w') as vtk_file:
    # write header
    vtk_file.write('# vtk DataFile Version 4.0\n')
    vtk_file.write('vtk output\n')
    vtk_file.write('ASCII\n')
    vtk_file.write('DATASET UNSTRUCTURED_GRID\n')

    # write points
    header5 = 'POINTS ' + str(NodesNr) + ' float\n'
    vtk_file.write(header5)
    for node in range(len(X)):
        coordinates = str(X[node]) + ' '
        coordinates += str(Y[node]) + ' '
        coordinates += str(Z[node]) + '\n'
        vtk_file.write(coordinates)

    # write cells
    header6 = '\nCELLS ' + str(mesh_info[1]) + ' ' + str(mesh_info[1]*5) + '\n'
    vtk_file.write(header6)
    for element in range(mesh_info[1]):
        connectivity = '4 '
        connectivity += str(elem_connectivities[element*4+0]-1) + ' '
        connectivity += str(elem_connectivities[element*4+1]-1) + ' '
        connectivity += str(elem_connectivities[element*4+2]-1) + ' '
        connectivity += str(elem_connectivities[element*4+3]-1) + '\n'
        vtk_file.write(connectivity)

    # write cell types
    header7 = '\nCELL_TYPES ' + str(mesh_info[1]) + '\n'
    vtk_file.write(header7)
    for _ in range(mesh_info[1]):
        vtk_file.write('9\n')

    # write point data
    header8 = '\nPOINT_DATA ' + str(mesh_info[0]) + '\n'
    vtk_file.write(header8)
    vtk_file.write('SCALARS CP float\n')
    vtk_file.write('LOOKUP_TABLE my_table\n')
    # for i in range(len(CP)):
    #     cp = str(CP[i]) + '\n'
    #     vtk_file.write(cp)

    for i in range(len(CP)):
        if CP[i] < 0.0:
            vtk_file.write('-5.0\n')
        else:
            vtk_file.write('5.0\n')

    vtk_file.write('\nVECTORS POINT_NORMALS float\n')
    for i in range(len(X)):
        point_normal = str(point_normals[3*i]) + ' '
        point_normal += str(point_normals[3*i]) + ' '
        point_normal += str(point_normals[3*i]) + '\n'
        vtk_file.write(point_normal)

    vtk_file.write('\nFIELD FieldData 1\n')
    header9 = 'FORCES 3 ' + str(mesh_info[0]) + ' float\n'
    vtk_file.write(header9)
    for i in range(len(X)):
        force = str(fluid_forces[3*i]) + ' '
        force += str(fluid_forces[3*i+1]) + ' '
        force += str(fluid_forces[3*i+2]) + '\n'
        vtk_file.write(force)

    header10 = '\nCELL_DATA ' + str(mesh_info[1]) + '\n'
    vtk_file.write(header10)
    # vtk_file.write('NORMALS cell_normals float\n')
    vtk_file.write('VECTORS NORMALS float\n')
    # header11 = 'NORMALS 3 ' + str(mesh_info[1]) + ' float\n'
    # vtk_file.write(header11)
    for cell in range(int(len(elem_connectivities)/4)):
        cell_normal = str(cell_normals[3*cell]) + ' '
        cell_normal += str(cell_normals[3*cell]) + ' '
        cell_normal += str(cell_normals[3*cell]) + '\n'
        vtk_file.write(cell_normal)

    # vtk_file.write('\nSCALARS areas float\n')
    # vtk_file.write('LOOKUP_TABLE areas_table\n')
    # for cell in range(int(len(elem_connectivities)/4)):
    #     area = str(cell_areas[cell]) + '\n'
    #     vtk_file.write(area)







