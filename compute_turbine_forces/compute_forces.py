import re, sys, json, os
import numpy as np
import matplotlib.pyplot as plt
import tau_functions_Steady as TauFunctionsSteady

working_path = os.getcwd() + '/'

# This function computes thrust and torque
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
    min_y = min(section_y)
    section_y_abs = []
    for i in range(len(section_y)):
        section_y_abs.append(abs(section_y[i]/min_y))
        # print(abs(section_y[i])/min_y)

    # sort both lists
    sorted_section_thrust = [i for _,i in sorted(zip(section_y_abs,section_thrust))]
    sorted_section_torque = [i for _,i in sorted(zip(section_y_abs,section_torque))]
    section_y_abs.sort()

    #'''
    # compute distance between sections
    delta_y = [0.5*(section_y_abs[1]-section_y_abs[0])]

    for i in range(len(section_y_abs)-2):
        delta_y.append(0.5*(section_y_abs[i+2]-section_y_abs[i]))

    delta_y.append(0.5*(section_y_abs[len(section_y_abs)-1]-section_y_abs[len(section_y_abs)-2]))

    # compute thrust distribution
    section_thrust_distribution = []
    section_torque_distribution = []
    for i in range(len(sorted_section_thrust)):
        section_thrust_distribution.append(sorted_section_thrust[i]/delta_y[i])
        section_torque_distribution.append(sorted_section_torque[i]/delta_y[i])
    #'''

    # compute total thrust and torque
    total_thrust = sum(sorted_section_thrust)
    total_torque = sum(sorted_section_torque)

    return section_thrust_distribution, section_torque_distribution, section_y_abs, total_thrust, total_torque

# This function writes a results file with the distributed thrust or torque
def WriteRadialDistributionResultsFile(output_filename, force, y):
    with open(output_filename, 'w') as output_file:
        for i in range(len(force)):
            output_file.write('{0:15f} {1:15f}\n'.format(y[i], force[i]))

# This function plots thrust and torque
def PlotThrustAndTorque(step, le_name, te_name, model_name, axs):
    output_file_pattern = 'airfoilSol.pval.'
    TauFunctionsSteady.ChangeFormat(working_path, step, "MEMBRANE_UP", "MEMBRANE_DOWN", output_file_pattern)
    TauFunctionsSteady.ChangeFormat(working_path, step, le_name, te_name, output_file_pattern)

    input_filename = 'Outputs/airfoilSol.surface.pval.' + str(step)
    output_filename_up = input_filename + '.MEMBRANE_UP.dat'
    output_filename_down = input_filename + '.MEMBRANE_DOWN.dat'
    output_filename_le = input_filename + '.' + le_name + '.dat'
    output_filename_te = input_filename + '.' + te_name + '.dat'

    thrust_up, torque_up, section_y_abs_up, total_thrust_up, total_torque_up = ComputeThrustAndTorque(output_filename_up)
    thrust_down, torque_down, section_y_abs_down, total_thrust_down, total_torque_down = ComputeThrustAndTorque(output_filename_down)
    thrust_le, torque_le, section_y_abs_le, total_thrust_le, total_torque_le = ComputeThrustAndTorque(output_filename_le)
    thrust_te, torque_te, section_y_abs_te, total_thrust_te, total_torque_te = ComputeThrustAndTorque(output_filename_te)

    thrust_distribution = []
    torque_distribution = []
    plot_y = []

    number_of_cut_sections = 19
    for i in range(len(thrust_up)-number_of_cut_sections):
        plot_y.append(section_y_abs_up[i])
        thrust_distribution.append(thrust_up[i] + thrust_down[i] + thrust_le[i] + thrust_te[i])
        torque_distribution.append(torque_up[i] + torque_down[i] + torque_le[i] + torque_te[i])

    # Write radial distribution files to import in latex
    thrust_output_filename = "thrust_radial_distribution_" + model_name + ".dat"
    torque_output_filename = "torque_radial_distribution_" + model_name + ".dat"
    WriteRadialDistributionResultsFile(thrust_output_filename, thrust_distribution, plot_y)
    WriteRadialDistributionResultsFile(torque_output_filename, torque_distribution, plot_y)

    total_thrust = total_thrust_up + total_thrust_down + total_thrust_le + total_thrust_te
    total_torque = total_torque_up + total_torque_down + total_torque_le + total_torque_te

    print('total_thrust_' + model_name + ' = ', round(total_thrust,2)*2.0, ' [N]')
    print('total_torque_' + model_name + ' = ', round(total_torque,2)*2.0, ' [Nm]')

    axs[0].plot(plot_y, thrust_distribution, label= model_name)
    axs[1].plot(plot_y, torque_distribution, label= model_name)
    # axs[0].plot(section_y_abs_up, thrust_up, label= 'thrust_up_' + model_name)
    # axs[0].plot(section_y_abs_up, thrust_down, label= 'thrust_down' + model_name)
    # axs[0].plot(section_y_abs_up, thrust_le, label= 'thrust_le' + model_name)
    # axs[0].plot(section_y_abs_up, thrust_te, label= 'thrust_te' + model_name)
    # plt.plot(section_y_abs_up, thrust_up, label= 'thrust_up_')
    # plt.plot(section_y_abs_up, thrust_down, label= 'thrust_down_')
    # plt.plot(section_y_abs_up, thrust_le, label= 'thrust_le_')
    # plt.plot(section_y_abs_up, thrust_te, label= 'thrust_te_')
    # plt.plot(section_y_abs_up, torque_up, label= 'torque_up_')
    # plt.plot(section_y_abs_up, torque_down, label= 'torque_down_')
    # plt.plot(section_y_abs_up, torque_le, label= 'torque_le_')
    # plt.plot(section_y_abs_up, torque_te, label= 'torque_te_')

###########################################################################
# Main script
###########################################################################

# preparing thrust and torque plots
fig, axs = plt.subplots(2)
# axs = 0

# Original case
step_original = 14001
le_name_original = "LEADING_EDGE"
te_name_original = "TRAILING_EDGE"
model_name_original = "original_cfd"
PlotThrustAndTorque(step_original, le_name_original, te_name_original, model_name_original, axs)

# Formfound case
step_formfound = 19000
le_name_formfound = "LE_FF"
te_name_formfound = "TE_FF"
model_name_formfound = "formfound_cfd"
PlotThrustAndTorque(step_formfound, le_name_formfound, te_name_formfound, model_name_formfound, axs)

# FSI case
step_fsi = 35000
le_name_fsi = "LE_FF"
te_name_fsi = "TE_FF"
model_name_fsi = "fsi"
PlotThrustAndTorque(step_fsi, le_name_fsi, te_name_fsi, model_name_fsi, axs)

# plotting thrust and torque distributions
fig.suptitle('Thrust and torque radial distributions')
axs[0].set(ylabel='Thrust [N/m]')
axs[1].set(ylabel='Torque [N]')
axs[0].legend(loc="upper left")
axs[0].grid()
axs[1].grid()
for ax in axs.flat:
    ax.set(xlabel='r/R [-]')
plt.show()

# plotting option 2:
# plt.title('Thrust and torque radial distributions')
# plt.legend(loc="upper left")
# plt.grid()
# plt.show()