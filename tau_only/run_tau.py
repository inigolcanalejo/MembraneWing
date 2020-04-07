# import modules for TAU
import PyPara
import PySolv

# Definition of the parameter file
para_path='/media/inigo/10740FB2740F9A1C/simulations/MembraneWing/tau_only/airfoil_Structured.cntl'

# Init Tau python classes
Para = PyPara.Parafile(para_path)
Solver = PySolv.Solver(para_path)

Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
Solver.outer_loop()
Solver.output()
Solver.finalize()
Para.free_parameters()