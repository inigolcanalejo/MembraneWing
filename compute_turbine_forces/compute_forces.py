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

def WriteRadialDistributionResultsFile(output_filename, force, y):
    with open(output_filename, 'w') as output_file:
        for i in range(len(force)):
            output_file.write('{0:15f} {1:15f}\n'.format(y[i], force[i]))

# Main script CFD Formfound
step = 19000
output_file_pattern = 'airfoilSol.pval.'
sub_step = 0
TauFunctionsSteady.ChangeFormat(working_path, step, "MEMBRANE_UP", "MEMBRANE_DOWN", output_file_pattern)
TauFunctionsSteady.ChangeFormat(working_path, step, "LE_FF", "TE_FF", output_file_pattern)

input_filename = 'Outputs/airfoilSol.surface.pval.19000'
output_filename_up = input_filename + '.MEMBRANE_UP.dat'
output_filename_down = input_filename + '.MEMBRANE_DOWN.dat'
output_filename_le = input_filename + '.LE_FF.dat'
output_filename_te = input_filename + '.TE_FF.dat'

thrust_up, torque_up, section_y_abs_up, total_thrust_up, total_torque_up = ComputeThrustAndTorque(output_filename_up)
thrust_down, torque_down, section_y_abs_down, total_thrust_down, total_torque_down = ComputeThrustAndTorque(output_filename_down)
thrust_le, torque_le, section_y_abs_le, total_thrust_le, total_torque_le = ComputeThrustAndTorque(output_filename_le)
thrust_te, torque_te, section_y_abs_te, total_thrust_te, total_torque_te = ComputeThrustAndTorque(output_filename_te)

thrust_distribution_formfinding = []
torque_distribution_formfinding = []
plot_y = []

number_of_cut_sections = 19
for i in range(len(thrust_up)-number_of_cut_sections):
    plot_y.append(section_y_abs_up[i])
    thrust_distribution_formfinding.append(thrust_up[i] + thrust_down[i] + thrust_le[i] + thrust_te[i])
    torque_distribution_formfinding.append(torque_up[i] + torque_down[i] + torque_le[i] + torque_te[i])

# Write radial distribution files to import in latex
thrust_output_filename_formfound = "thrust_radial_distribution_formfound_cfd.dat"
torque_output_filename_formfound = "torque_radial_distribution_formfound_cfd.dat"
WriteRadialDistributionResultsFile(thrust_output_filename_formfound, thrust_distribution_formfinding, plot_y)
WriteRadialDistributionResultsFile(torque_output_filename_formfound, torque_distribution_formfinding, plot_y)

total_thrust_formfinding = total_thrust_up + total_thrust_down + total_thrust_le + total_thrust_te
total_torque_formfinding = total_torque_up + total_torque_down + total_torque_le + total_torque_te

print('total_thrust_formfinding = ', round(total_thrust_formfinding,2)*2.0, ' [N]')
print('total_torque_formfinding = ', round(total_torque_formfinding,2)*2.0, ' [Nm]')

''' test
print('total_thrust_formfinding = ', round(sum(thrust_distribution_formfinding),2), ' [N]')
print('total_torque_formfinding = ', round(sum(torque_distribution_formfinding),2), ' [Nm]')
#'''


# Main script CFD ORIGINAL
step = 14001
output_file_pattern = 'airfoilSol.pval.'
sub_step = 0
TauFunctionsSteady.ChangeFormat(working_path, step, "MEMBRANE_UP", "MEMBRANE_DOWN", output_file_pattern)
TauFunctionsSteady.ChangeFormat(working_path, step, "LEADING_EDGE", "TRAILING_EDGE", output_file_pattern)

input_filename = 'Outputs/airfoilSol.surface.pval.14001'
output_filename_up = input_filename + '.MEMBRANE_UP.dat'
output_filename_down = input_filename + '.MEMBRANE_DOWN.dat'
output_filename_le = input_filename + '.LEADING_EDGE.dat'
output_filename_te = input_filename + '.TRAILING_EDGE.dat'

thrust_up_original,   torque_up_original,   section_y_abs_up_original,   total_thrust_up_original,    total_torque_up_original   = ComputeThrustAndTorque(output_filename_up)
thrust_down_original, torque_down_original, section_y_abs_down_original, total_thrust_down_original,  total_torque_down_original = ComputeThrustAndTorque(output_filename_down)
thrust_le_original,   torque_le_original,   section_y_abs_le_original,   total_thrust_le_original,    total_torque_le_original   = ComputeThrustAndTorque(output_filename_le)
thrust_te_original,   torque_te_original,   section_y_abs_te_original,   total_thrust_te_original,    total_torque_te_original   = ComputeThrustAndTorque(output_filename_te)

thrust_distribution_original = []
torque_distribution_original = []
plot_y_original = []

number_of_cut_sections = 12
for i in range(len(thrust_up_original)-number_of_cut_sections):
    plot_y_original.append(section_y_abs_down_original[i])
    thrust_distribution_original.append(thrust_up_original[i] + thrust_down_original[i] + thrust_le_original[i] + thrust_te_original[i])
    torque_distribution_original.append(torque_up_original[i] + torque_down_original[i] + torque_le_original[i] + torque_te_original[i])

# Write radial distribution files to import in latex
thrust_output_filename_original = "thrust_radial_distribution_original_cfd.dat"
torque_output_filename_original = "torque_radial_distribution_original_cfd.dat"
WriteRadialDistributionResultsFile(thrust_output_filename_original, thrust_distribution_original, plot_y_original)
WriteRadialDistributionResultsFile(torque_output_filename_original, torque_distribution_original, plot_y_original)

total_thrust_original = total_thrust_up_original + total_thrust_down_original + total_thrust_le_original + total_thrust_te_original
total_torque_original = total_torque_up_original + total_torque_down_original + total_torque_le_original + total_torque_te_original

print('total_thrust_original = ', round(total_thrust_original,2)*2.0, ' [N]')
print('total_torque_original = ', round(total_torque_original,2)*2.0, ' [Nm]')

''' test
print('total_thrust_original = ', round(sum(thrust_distribution_original),2), ' [N]')
print('total_torque_original = ', round(sum(torque_distribution_original),2), ' [Nm]')
#'''

# Main script FSI
step = 27000
output_file_pattern = 'airfoilSol.pval.'
sub_step = 0
TauFunctionsSteady.ChangeFormat(working_path, step, "MEMBRANE_UP", "MEMBRANE_DOWN", output_file_pattern)
TauFunctionsSteady.ChangeFormat(working_path, step, "LE_FF", "TE_FF", output_file_pattern)

input_filename = 'Outputs/airfoilSol.surface.pval.27000'
output_filename_up = input_filename + '.MEMBRANE_UP.dat'
output_filename_down = input_filename + '.MEMBRANE_DOWN.dat'
output_filename_le = input_filename + '.LE_FF.dat'
output_filename_te = input_filename + '.TE_FF.dat'

thrust_up_fsi,   torque_up_fsi,   section_y_abs_up_fsi,   total_thrust_up_fsi,    total_torque_up_fsi   = ComputeThrustAndTorque(output_filename_up)
thrust_down_fsi, torque_down_fsi, section_y_abs_down_fsi, total_thrust_down_fsi,  total_torque_down_fsi = ComputeThrustAndTorque(output_filename_down)
thrust_le_fsi,   torque_le_fsi,   section_y_abs_le_fsi,   total_thrust_le_fsi,    total_torque_le_fsi   = ComputeThrustAndTorque(output_filename_le)
thrust_te_fsi,   torque_te_fsi,   section_y_abs_te_fsi,   total_thrust_te_fsi,    total_torque_te_fsi   = ComputeThrustAndTorque(output_filename_te)

thrust_distribution_fsi = []
torque_distribution_fsi = []
plot_y_fsi = []

number_of_cut_sections = 19
for i in range(len(thrust_up_fsi)-number_of_cut_sections):
    plot_y_fsi.append(section_y_abs_down_fsi[i])
    thrust_distribution_fsi.append(thrust_up_fsi[i] + thrust_down_fsi[i] + thrust_le_fsi[i] + thrust_te_fsi[i])
    torque_distribution_fsi.append(torque_up_fsi[i] + torque_down_fsi[i] + torque_le_fsi[i] + torque_te_fsi[i])

# Write radial distribution files to import in latex
thrust_output_filename_fsi = "thrust_radial_distribution_fsi_cfd.dat"
torque_output_filename_fsi = "torque_radial_distribution_fsi_cfd.dat"
WriteRadialDistributionResultsFile(thrust_output_filename_fsi, thrust_distribution_fsi, plot_y_fsi)
WriteRadialDistributionResultsFile(torque_output_filename_fsi, torque_distribution_fsi, plot_y_fsi)

TOTAL_THRUST_FSI = total_thrust_up_fsi + total_thrust_down_fsi + total_thrust_le_fsi + total_thrust_te_fsi
TOTAL_TORQUE_FSI = total_torque_up_fsi + total_torque_down_fsi + total_torque_le_fsi + total_torque_te_fsi

print('TOTAL_THRUST_FSI = ', round(TOTAL_THRUST_FSI,2)*2.0, ' [N]')
print('TOTAL_TORQUE_FSI = ', round(TOTAL_TORQUE_FSI,2)*2.0, ' [Nm]')

''' test
print('TOTAL_THRUST_FSI = ', round(sum(total_thrust_FSI),2), ' [N]')
print('TOTAL_TORQUE_FSI = ', round(sum(total_torque_FSI),2), ' [Nm]')
#'''


# plotting thrust and torque distributions
plt.title('Torque and Thrust radial distributions')
plt.xlabel('r/R [-]')
#plt.xlabel('Radius [m]')

# Thrust
plt.plot(plot_y, thrust_distribution_formfinding, label='Thrustl Formfinding [N/m]')
plt.plot(plot_y_original, thrust_distribution_original, label='Thrust NASA [N/m]')
plt.plot(plot_y_fsi, thrust_distribution_fsi, label='Thrust Membrane Blade [N/m]')

# Torque
# plt.plot(plot_y, torque_distribution_formfinding, label='Torque_total Formfinding [N]')
# plt.plot(plot_y_original, torque_distribution_original, label='Torque_total NASA [N]')
# plt.plot(plot_y_fsi, torque_distribution_fsi, label='Torque Membrane Blade [N]')

plt.legend(loc="upper left")
plt.show()

# FORMFINDING CFD
# Thrust
# plt.plot(section_y_abs_up, thrust_up, label='Thrust_up [N]')
# plt.plot(section_y_abs_up, thrust_down, label='Thrust_down [N]')
# plt.plot(section_y_abs_up, thrust_le, label='Thrust_le [N]')
# plt.plot(section_y_abs_up, thrust_te, label='Thrust_te [N]')
# plt.plot(plot_y, thrust_distribution_formfinding, label='Thrust_total Formfinding [N/m]')

# Torque
# plt.plot(section_y_abs_up, torque_up, label='Torque_up [N]')
# plt.plot(section_y_abs_up, torque_down, label='Torque_down [N]')
# plt.plot(section_y_abs_up, torque_le, label='Torque_le [N]')
# plt.plot(section_y_abs_up, torque_te, label='Torque_te [N]')
# plt.plot(plot_y, torque_distribution_formfinding, label='Torque_total Formfinding [N]')
# plt.legend(loc="upper left")
# plt.show()

# NASA-Ames Phase VI CFD
# Thrust
# plt.plot(section_y_abs_up, thrust_up, label='Thrust_up [N]')
# plt.plot(section_y_abs_up, thrust_down, label='Thrust_down [N]')
# plt.plot(section_y_abs_up, thrust_le, label='Thrust_le [N]')
# plt.plot(section_y_abs_up, thrust_te, label='Thrust_te [N]')
# plt.plot(plot_y_original, thrust_distribution_formfinding, label='Thrust_total NASA [N/m]')

# Torque
# plt.plot(section_y_abs_up, torque_up, label='Torque_up [N]')
# plt.plot(section_y_abs_up, torque_down, label='Torque_down [N]')
# plt.plot(section_y_abs_up, torque_le, label='Torque_le [N]')
# plt.plot(section_y_abs_up, torque_te, label='Torque_te [N]')
# plt.plot(plot_y_original, torque_distribution_formfinding, label='Torque_total NASA [N]')
# plt.legend(loc="upper left")
# plt.show()