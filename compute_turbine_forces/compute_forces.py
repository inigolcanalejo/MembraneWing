import re, sys, json, os
import numpy as np
import matplotlib.pyplot as plt

working_path = os.getcwd() + '/'
with open(working_path + 'input/tau_settings.json') as json_file:
    tau_settings = json.load(json_file)

start_step = tau_settings["start_step"]
tau_path = tau_settings["tau_path"]
sys.path.append(tau_settings["kratos_path"])
sys.path.append(tau_path + "py_turb1eq/")

# tau_functions can only be imported after appending kratos' path
import tau_functions_Steady as TauFunctionsSteady
import tau_functions as TauFunctions

def ComputeThrustAndTorque(output_filename):
    print('Computing: ', output_filename)
    # Read output file
    position_info, mesh_info, nodal_data, elem_connectivities = TauFunctionsSteady.ReadInterfaceFile(output_filename)

    # Read mesh info
    NodesNr = mesh_info[0]
    print("NodesNr = ", NodesNr)

    X, Y, Z = TauFunctionsSteady.SaveCoordinatesList(nodal_data, position_info, NodesNr)

    velocity = 41.09
    P = TauFunctionsSteady.SavePressure(nodal_data, position_info, NodesNr, velocity)

    # calculating the force vector
    fluid_forces = TauFunctionsSteady.CalculateNodalFluidForces(X, Y, Z, P, elem_connectivities)

    # list containing "y" coordinate of each section (not ordered)
    section_y = []
    section_y.append(Y[0])

    # list containing thrust of each section
    section_thrust = [0]
    section_torque = [0]

    # loop over y coordinates of the nodes
    for i in range(len(Y)):
        # initialize flag
        section_found = False

        # loop over already found sections
        for j in range(len(section_y)):
            # if distance is below tolerance, the section already exists
            if abs(Y[i]-section_y[j]) < 1.99e-3:
                section_found = True

                # add force_x contribution to section thrust
                section_thrust[j] += fluid_forces[3*i]
                section_torque[j] += fluid_forces[3*i+2]*abs(Y[i])

        # if section does not exist, add y and force_x to the lists of sections
        if not section_found:
            # add y
            section_y.append(Y[i])
            # add force_x
            section_thrust.append(fluid_forces[3*i])
            section_torque.append(fluid_forces[3*i+2]*abs(Y[i]))

    # list containing relative value of y coordinate of each section
    section_y_abs = []
    for i in range(len(section_y)):
        section_y_abs.append(abs(section_y[i]))

    # sort both lists
    sorted_section_thrust = [i for _,i in sorted(zip(section_y_abs,section_thrust))]
    sorted_section_torque = [i for _,i in sorted(zip(section_y_abs,section_torque))]
    section_y_abs.sort()

    '''
    # compute distance between sections
    delta_y = [0.5*(section_y_abs[1]-section_y_abs[0])]

    for i in range(len(section_y_abs)-2):
        delta_y.append(0.5*(section_y_abs[i+2]-section_y_abs[i]))

    delta_y.append(0.5*(section_y_abs[len(section_y_abs)-1]-section_y_abs[len(section_y_abs)-2]))

    # compute thrust distribution
    section_thrust_distribution = []
    for i in range(len(sorted_section_thrust)):
        section_thrust_distribution.append(sorted_section_thrust[i]/delta_y[i])
    '''

    return sorted_section_thrust, sorted_section_torque, section_y_abs

# Main script
step = 1
output_file_pattern = 'airfoilSol.pval.deform_i='
sub_step = 0
TauFunctions.ChangeFormat(working_path, step, "MEMBRANE_UP", "MEMBRANE_DOWN", output_file_pattern, sub_step)
TauFunctions.ChangeFormat(working_path, step, "LE_FF", "TE_FF", output_file_pattern, sub_step)

input_filename = 'Outputs/airfoilSol.surface.pval.deform_i=1_t=2.315e-03_subiter=0'
output_filename_up = input_filename + '.MEMBRANE_UP.dat'
output_filename_down = input_filename + '.MEMBRANE_DOWN.dat'
output_filename_le = input_filename + '.LE_FF.dat'
output_filename_te = input_filename + '.TE_FF.dat'

thrust_up, torque_up, section_y_abs_up = ComputeThrustAndTorque(output_filename_up)
thrust_down, torque_down, section_y_abs_down = ComputeThrustAndTorque(output_filename_down)
thrust_le, torque_le, section_y_abs_le = ComputeThrustAndTorque(output_filename_le)
thrust_te, torque_te, section_y_abs_te = ComputeThrustAndTorque(output_filename_te)

total_thrust = []
total_torque = []

for i in range(len(thrust_up)):
    total_thrust.append(thrust_up[i] + thrust_down[i] + thrust_le[i] + thrust_te[i])
    total_torque.append(torque_up[i] + torque_down[i] + torque_le[i] + torque_te[i])

# plotting thrust and torque distributions
plt.title('Torque and Thrust radial distributions')
plt.xlabel('Radius [m]')

# Thrust
# plt.plot(section_y_abs_up, thrust_up, label='Thrust_up [N]')
# plt.plot(section_y_abs_up, thrust_down, label='Thrust_down [N]')
# plt.plot(section_y_abs_up, thrust_le, label='Thrust_le [N]')
# plt.plot(section_y_abs_up, thrust_te, label='Thrust_te [N]')
# plt.plot(section_y_abs_up, total_thrust, label='Thrust_total [N]')

# Torque
plt.plot(section_y_abs_up, torque_up, label='Torque_up [Nm]')
plt.plot(section_y_abs_up, torque_down, label='Torque_down [Nm]')
plt.plot(section_y_abs_up, torque_le, label='Torque_le [Nm]')
plt.plot(section_y_abs_up, torque_te, label='Torque_te [Nm]')
plt.plot(section_y_abs_up, total_torque, label='Torque_total [Nm]')
plt.legend(loc="upper left")
plt.show()