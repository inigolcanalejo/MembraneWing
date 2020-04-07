#Definition
Work_Dir=/home/inigo/simulations/membranProjekt/works/works180718/Ploetz_FSI_instationaer_U20_AOA6
Casename=airfoil_Structured

#*********************************
Unsteady_physical_time_steps=200
Unsteady_inner_iterations_per_time_step=1
n_outer_out=2
fluidIter=50
Maximal_step_number=200
Minimum_residual=1e-16
Step_size=0.005
AOA=6
#*********************************

Para_Path=$Work_Dir/$Casename.cntl
Grid_File=$Work_Dir/Mesh/airfoil_Structured_scaliert.grid

Tautoplt=$Work_Dir/TautopltTest.cntl
#*********************************

#Tautoplt
sed 's|'"Primary grid filename: TBD"'|'"Primary grid filename: $Grid_File.def.100"'|g' -i /$Tautoplt
sed 's|'"Boundary mapping filename: TBD"'|'"Boundary mapping filename: $Para_Path"'|g' -i /$Tautoplt
sed 's|'"Restart-data prefix: TBD"'|'"Restart-data prefix: $Work_Dir/Ergebnisse/Outputs/airfoilSol.pval.unsteady_i=100_t=5.000e-01"'|g' -i /$Tautoplt
#**********************************

/home/inigo/software/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/taubin_svn19618.OPENMPI1.6.4_Python2.7.5/bin/tau2plt TautopltTest.cntl
#TautopltTest.cntl muss bei änderung von n_outer_out angepasst werden!!!!!!!

#TAUtoplt
sed 's|'"Primary grid filename: $Grid_File.def.100"'|'"Primary grid filename: TBD"'|g' -i /$Tautoplt
sed 's|'"Boundary mapping filename: $Para_Path"'|'"Boundary mapping filename: TBD"'|g' -i /$Tautoplt
sed 's|'"Restart-data prefix: $Work_Dir/Ergebnisse/Outputs/airfoilSol.pval.unsteady_i=100_t=5.000e-01"'|'"Restart-data prefix: TBD"'|g' -i /$Tautoplt

