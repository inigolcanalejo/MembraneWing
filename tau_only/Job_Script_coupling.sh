#
# . /etc/profile
# . /etc/profile.d/modules.sh
#
# module unload mpi.intel
# module unload intel
# module load intel/15.0
# module load mpi.ompi/1.8/intel
# module load mpi4py/1.3.1
# module load tempdir/1.0
# module unload gcc
# module load cmake/3.4
# module load mkl/11.3
# module load geos/3.3.3
# module load  obspy/0.9.2
# module unload python
# module load python/2.7.6

#Definition

Work_Dir=/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only
Casename=airfoil_Structured

#*********************************
Unsteady_physical_time_steps=3
Unsteady_inner_iterations_per_time_step=2
n_outer_out=1
fluidIter=50
Maximal_step_number=200
Minimum_residual=1e-16
Step_size=0.005
AOA=6
#*********************************

Para_Path=$Work_Dir/$Casename.cntl
Grid_File=$Work_Dir/Mesh/airfoil_Structured_scaliert.grid
Output=$Work_Dir/Outputs
Kopl=$Work_Dir/TAUclient_coupling_main_mod_mitCarat_test3.py
Funct=$Work_Dir/Functions_mod3.py
Tautoplt=$Work_Dir/Tautoplt.cntl
#*********************************

rm -r -f Output_Files
rm -r -f Ergebnisse


#Parameterdatei

sed 's|'"Primary grid filename: TBD"'|'"Primary grid filename: $Grid_File"'|g' -i /$Para_Path
sed 's|'"Output files prefix: TBD"'|'"Output files prefix: $Output/airfoilSol"'|g' -i /$Para_Path
sed 's|'"Maximal time step number: TBD"'|'"Maximal time step number: $Maximal_step_number"'|g' -i /$Para_Path
sed 's|'"Minimum residual: TBD"'|'"Minimum residual: $Minimum_residual"'|g' -i /$Para_Path
sed 's|'"Unsteady physical time steps: TBD"'|'"Unsteady physical time steps: $Unsteady_physical_time_steps"'|g' -i /$Para_Path
sed 's|'"Unsteady inner iterations per time step: TBD"'|'"Unsteady inner iterations per time step: $Unsteady_inner_iterations_per_time_step"'|g' -i /$Para_Path
sed 's|'"Angle alpha (degree): TBD"'|'"Angle alpha (degree): $AOA"'|g' -i /$Para_Path
sed 's|'"Unsteady physical time step size: TBD"'|'"Unsteady physical time step size: $Step_size"'|g' -i /$Para_Path

#Kopplungsdatei

sed 's|'"para_path='TBD'"'|'"para_path='$Para_Path'"'|g' -i /$Kopl
sed 's|'"n_outer_out=TBD"'|'"n_outer_out=$n_outer_out"'|g' -i /$Kopl
sed 's|'"fluidIter=TBD"'|'"fluidIter=$fluidIter"'|g' -i /$Kopl
sed 's|'"sys.path.append('TBD')"'|'"sys.path.append('$Work_Dir')"'|g' -i /$Kopl
sed 's|'"search_text = 'Primary grid filename: TBD'"'|'"search_text = 'Primary grid filename: $Grid_File'"'|g' -i /$Kopl
sed 's|'"replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)"'|'"replace_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out)"'|g' -i /$Kopl
sed 's|'"for line in fileinput.input('TBD',inplace=1):"'|'"for line in fileinput.input('$Para_Path',inplace=1):"'|g' -i /$Kopl
sed 's|'"search_text = 'Primary grid filename: TBD.def.' + str(this_step_out-1)"'|'"search_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out-1)"'|g' -i /$Kopl
sed 's|'"replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)"'|'"replace_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out)"'|g' -i /$Kopl
sed 's|'"for line in fileinput.input('TBD',inplace=1):"'|'"for line in fileinput.input('$Para_Path',inplace=1):"'|g' -i /$Kopl
sed 's|'"List_file = Functions.findFile('TBD')"'|'"List_file = Functions.findFile('$Output/*.plt')"'|g' -i /$Kopl

#Functions_mod

sed 's|'"para_path= 'TBD'"'|'"para_path= '$Para_Path'"'|g' -i /$Funct
sed 's|'"        if file.startswith('TBD/airfoilSol.MEMBRANE_i=%d'%this_step_out):"'|'"        if file.startswith('$Output/airfoilSol.MEMBRANE_i=%d'%this_step_out):"'|g' -i /$Funct

#Tautoplt
sed 's|'"Primary grid filename: TBD"'|'"Primary grid filename: $Grid_File.def.1"'|g' -i /$Tautoplt
sed 's|'"Boundary mapping filename: TBD"'|'"Boundary mapping filename: $Para_Path"'|g' -i /$Tautoplt
sed 's|'"Restart-data prefix: TBD"'|'"Restart-data prefix: $Work_Dir/Ergebnisse/Outputs/airfoilSol.pval.unsteady_i=1_t=5.000e-03"'|g' -i /$Tautoplt
#**********************************

mkdir Outputs

#FSI beginnen


#source /home/inigo/software/pyEmpire/EMPIRE-Core/etc/bashrc.sh ICC

#cd /home/inigo/software/TAU/tau_only/

mpirun -n 1  /home/inigo/software/carat/src/carat CARAT/caratInput_FSI > log_Carat_neu_1 &

sleep 10s

mpirun -n 1  python /home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/py_turb1eq/tau.py TAUclient_coupling_main_mod_mitCarat_test3.py airfoil_Structured.cntl log_TAU_neu_1.out

#**************************************

#Clean files

#***************************************

mkdir Ergebnisse
mkdir Output_Files

mv $Work_Dir/walldistances_matrix $Work_Dir/Output_Files
mv $Work_Dir/xp $Work_Dir/Output_Files
mv $Work_Dir/xpNode $Work_Dir/Output_Files
mv $Work_Dir/deformation_file.tec  $Work_Dir/Output_Files
mv $Work_Dir/ElemTable  $Work_Dir/Output_Files
mv $Work_Dir/interface_deformfile.nc $Work_Dir/Output_Files
mv $Work_Dir/ping.port $Work_Dir/Output_Files
mv $Work_Dir/RBF_matrix $Work_Dir/Output_Files
mv $Work_Dir/sta_geo_nonlin_load_disp_curve.dat $Work_Dir/Output_Files
mv $Work_Dir/log_Carat_neu_1 $Work_Dir/Ergebnisse
mv $Work_Dir/log_TAU_neu_1.out.Python-Interface.stderr $Work_Dir/Ergebnisse
mv $Work_Dir/log_TAU_neu_1.out.Python-Interface.stdout $Work_Dir/Ergebnisse
mv $Work_Dir/out.err $Work_Dir/Ergebnisse
mv $Work_Dir/out.log $Work_Dir/Ergebnisse
mv $Work_Dir/Outputs $Work_Dir/Ergebnisse

/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/tau2plt Tautoplt.cntl
#Tautoplt.cntl muss bei änderung von n_outer_out angepasst werden!!!!!!!

#Parameterdatei

sed 's|'"Primary grid filename: $Grid_File.*"'|'"Primary grid filename: TBD"'|g' -i /$Para_Path
sed 's|'"Output files prefix: $Output/airfoilSol"'|'"Output files prefix: TBD"'|g' -i /$Para_Path
sed 's|'"Maximal time step number: $Maximal_step_number"'|'"Maximal time step number: TBD"'|g' -i /$Para_Path
sed 's|'"Minimum residual: $Minimum_residual"'|'"Minimum residual: TBD"'|g' -i /$Para_Path
sed 's|'"Unsteady physical time steps: $Unsteady_physical_time_steps"'|'"Unsteady physical time steps: TBD"'|g' -i /$Para_Path
sed 's|'"Unsteady inner iterations per time step: $Unsteady_inner_iterations_per_time_step"'|'"Unsteady inner iterations per time step: TBD"'|g' -i /$Para_Path
sed 's|'"Angle alpha (degree): $AOA"'|'"Angle alpha (degree): TBD"'|g' -i /$Para_Path
sed 's|'"Unsteady physical time step size: $Step_size"'|'"Unsteady physical time step size: TBD"'|g' -i /$Para_Path

#Kopplungsdatei

sed 's|'"para_path='$Para_Path'"'|'"para_path='TBD'"'|g' -i /$Kopl
sed 's|'"n_outer_out=$n_outer_out"'|'"n_outer_out=TBD"'|g' -i /$Kopl
sed 's|'"fluidIter=$fluidIter"'|'"fluidIter=TBD"'|g' -i /$Kopl
sed 's|'"sys.path.append('$Work_Dir')"'|'"sys.path.append('TBD')"'|g' -i /$Kopl
sed 's|'"search_text = 'Primary grid filename: $Grid_File'"'|'"search_text = 'Primary grid filename: TBD'"'|g' -i /$Kopl
sed 's|'"replace_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out)"'|'"replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)"'|g' -i /$Kopl
sed 's|'"for line in fileinput.input('$Para_Path',inplace=1):"'|'"for line in fileinput.input('TBD',inplace=1):"'|g' -i /$Kopl
sed 's|'"search_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out-1)"'|'"search_text = 'Primary grid filename: TBD.def.' + str(this_step_out-1)"'|g' -i /$Kopl
sed 's|'"replace_text = 'Primary grid filename: $Grid_File.def.' + str(this_step_out)"'|'"replace_text = 'Primary grid filename: TBD.def.' + str(this_step_out)"'|g' -i /$Kopl
sed 's|'"for line in fileinput.input('$Para_Path',inplace=1):"'|'"for line in fileinput.input('TBD',inplace=1):"'|g' -i /$Kopl
sed 's|'"List_file = Functions.findFile('$Output/*.plt')"'|'"List_file = Functions.findFile('TBD')"'|g' -i /$Kopl

#Functions_mod

sed 's|'"para_path= '$Para_Path'"'|'"para_path= 'TBD'"'|g' -i /$Funct
sed 's|'"        if file.startswith('$Output/airfoilSol.MEMBRANE_i=%d'%this_step_out):"'|'"        if file.startswith('TBD/airfoilSol.MEMBRANE_i=%d'%this_step_out):"'|g' -i /$Funct

#TAUtoplt

sed 's|'"Primary grid filename: $Grid_File.def.1"'|'"Primary grid filename: TBD"'|g' -i /$Tautoplt
sed 's|'"Boundary mapping filename: $Para_Path"'|'"Boundary mapping filename: TBD"'|g' -i /$Tautoplt
sed 's|'"Restart-data prefix: $Work_Dir/Ergebnisse/Outputs/airfoilSol.pval.unsteady_i=1_t=5.000e-03"'|'"Restart-data prefix: TBD"'|g' -i /$Tautoplt

