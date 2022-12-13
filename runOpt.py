"""
AWAKE Run 2 Electron Line Model
R. Ramjiawan
Jun 2022
Track beam through line and extract beam parameters at merge-point
"""

import OptEnv as opt_env
import errorOptEnv as errorEnv
import plot_save_output as plot
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec
import pickle


# Initial values for quadrupoles (q), sextupoles (s), octupoles (o) and distances (a)
# nominal

a0 = 0
a1 = 0
a2 = 0
a3 = 0
a4 = 0
a5 = 0
a6 = 0
a7 = 0
a8 = 0
a9 = 0

# #ptc 5 um
q0 = 1.3396464743085914
q1 = 4.613669120362551
q2 = 7.063555777716582
q3 = 5.05436655065204
q4 = -4.89808966183322
q5 = -4.439862101295262
s0 = 52.646836071756844
s1 = -272.7838088808968
s2 = 77.04320547301776
s3 = -31.279857090658226
s4 = -45.210626979586145
s5 = -274.5906112848346
o0 = -4138.673920773547
o1 = 8138.8390297326505
o2 = 1292.4421785900909
o3 = 1403.4746408686688

# # enter the values here over which to optimise, otherwise hard-code them into MAD-X file
x = {
    0: {'name': 'quad0', 'strength': q0, 'type': 'quadrupole', 'norm': 7},
    1: {'name': 'quad1', 'strength': q1, 'type': 'quadrupole', 'norm': 7},
    2: {'name': 'quad2', 'strength': q2, 'type': 'quadrupole', 'norm': 10},
    3: {'name': 'quad3', 'strength': q3, 'type': 'quadrupole', 'norm': 7},
    4: {'name': 'quad4', 'strength': q4, 'type': 'quadrupole', 'norm': 7},
    5: {'name': 'quad5', 'strength': q5, 'type': 'quadrupole', 'norm': 7},
    6: {'name': 'sext0', 'strength': s0, 'type': 'sextupole', 'norm': 1000},
    7: {'name': 'sext1', 'strength': s1, 'type': 'sextupole', 'norm': 1000},
    8: {'name': 'sext2', 'strength': s2, 'type': 'sextupole', 'norm': 1000},
    9: {'name': 'sext3', 'strength': s3, 'type': 'sextupole', 'norm': 1000},
    10: {'name': 'sext4', 'strength': s4, 'type': 'sextupole', 'norm': 1000},
    11: {'name': 'sext5', 'strength': s5, 'type': 'sextupole', 'norm': 1000},
    12: {'name': 'oct0', 'strength': o0, 'type': 'octupole', 'norm': 5000},
    13: {'name': 'oct1', 'strength': o1, 'type': 'octupole', 'norm': 10000},
    14: {'name': 'oct2', 'strength': o2, 'type': 'octupole', 'norm': 3000},
    15: {'name': 'oct3', 'strength': o3, 'type': 'octupole', 'norm': 3000},
    # 15: {'name': 'dist0', 'strength': a0, 'type': 'distance', 'norm': 0.2},
    # 16: {'name': 'dist1', 'strength': a1, 'type': 'distance', 'norm': 0.2},
    # 17: {'name': 'dist2', 'strength': a2, 'type': 'distance', 'norm': 0.2},
    # 18: {'name': 'dist3', 'strength': a3, 'type': 'distance', 'norm': 0.1},
    # 19: {'name': 'dist4', 'strength': a4, 'type': 'distance', 'norm': 0.2},
    # 20: {'name': 'dist5', 'strength': a5, 'type': 'distance', 'norm': 0.3},
    # 21: {'name': 'dist6', 'strength': a6, 'type': 'distance', 'norm': 0.3},
    # 22: {'name': 'dist7', 'strength': a7, 'type': 'distance', 'norm': 0.3},
    # 23: {'name': 'dist8', 'strength': a8, 'type': 'distance', 'norm': 0.3}
}

# Specify parameters for optimisation
solver = 'Powell'
n_iter = 20
n_particles = 10000  # Used to generate distribution to track
foil_w = 100e-6      # Thickness of scattering foil
init_dist = []
thin = False         # Use thin track vs. PTC track
file = 'distr/Ellipse_150MeV_nominal.tfs'

# Initialise environment
env = opt_env.kOptEnv(solver, n_particles, n_iter, init_dist, foil_w, x, thin=thin)

# Initialise input distribution
var = []
f = open(file, 'r')  # initialize empty array
for line in f:
    var.append(
        line.strip().split())
f.close()
init_dist = np.array(var)[0:n_particles, 0:6].astype(np.float)
env.init_dist = init_dist
del var

# Either use optimiser (solution) or just output as is (step)
# If don't use step, script will run with values from general_tt43_python
if solver == "pyMOO":
    env.step_MO(env.norm_data([y['strength'] for y in x.values()]))
    plot = plot.Plot(env.madx, env.x_best, x, init_dist, foil_w, env.output_all, env.x_all)
    plot.twiss()
else:
    env.step(env.norm_data([y['strength'] for y in x.values()]))

# Optimise
if solver != "pyMOO":
    # solution = minimize(env.step, env.norm_data([y['strength'] for y in x.values()]), method=solver, options={'maxfev':n_iter})
    plot = plot.Plot(env.madx, env.x_best, x, init_dist, foil_w, env.output_all, env.x_all)
    plot.twiss()
    plot.plotmat_twiss()
    plot.plot1()
    env.step(env.norm_data([y['strength'] for y in x.values()]))
    plot.error()
    plot.plotmat()
else:
    from pymoo.model.problem import Problem
    from pymoo.algorithms.nsga2 import NSGA2
    from pymoo.algorithms.so_genetic_algorithm import GA
    from pymoo.factory import get_sampling, get_crossover, get_mutation
    from pymoo.optimize import minimize
    from pymoo.factory import get_termination
    from pymoo.visualization.scatter import Scatter

    x_0 = env.norm_data([y['strength'] for y in x.values()])
    norm_vect = env.norm_data([y['norm'] for y in x.values()])  # normalise
    n_obj = 1     # number of objective functions

    class MatchingProblem(opt_env.kOptEnv, Problem):
        def __init__(self,
                     norm_vect,
                     x_0,
                     n_var=len(x_0),
                     n_obj=n_obj,
                     n_constr=0,
                     xl=None,
                     xu=None):
            opt_env.kOptEnv.__init__(self, solver, n_particles, n_iter, init_dist, foil_w, x, thin=thin)
            Problem.__init__(self,
                             n_var=len(x_0),
                             n_obj=n_obj,
                             n_constr=n_constr,
                             xl=-np.ones(np.shape(norm_vect)),
                             xu=np.ones(np.shape(norm_vect)))

        def _evaluate(self, x_n, out, *args, **kwargs):
            f = []
            for j in range(x_n.shape[0]):
                y_raw_all, y_raw_single = self.step_MO(x_n[j, :])

                if self.n_obj == 1:
                    f.append(y_raw_single)
                else:
                    f.append(y_raw_all)
            out["F"] = np.vstack(f)

    problem = MatchingProblem(
        norm_vect=norm_vect,
        n_var=len(x_0),
        n_obj=n_obj,
        n_constr=0,      # define constraints if required
        x_0=x_0,
        xl=-np.ones(np.shape(norm_vect)),
        xu=np.ones(np.shape(norm_vect)))

    # problem.evaluate(np.vstack([problem.x_0, problem.x_0, -np.ones_like(problem.x_0)]))

    algorithm = GA(
        pop_size=200,
        n_offsprings=200,
        sampling=get_sampling("real_lhs"),
        crossover=get_crossover("real_sbx", prob=0.9, eta=15),
        mutation=get_mutation("real_pm", eta=30),
        eliminate_duplicates=True
    )

    # termination = MultiObjectiveDefaultTermination(
    #     x_tol=1e-8,
    #     cv_tol=1e-6,
    #     f_tol=1e-7,
    #     nth_gen=5,
    #     n_last=30,
    #     n_max_gen=50000,
    #     n_max_evals=200000
    # )
    termination = get_termination("n_eval", n_iter)

    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=1,
                   copy_algorithm=False,
                   verbose=True)

    ps = problem.pareto_set(use_cache=False, flatten=False)
    pf = problem.pareto_front(use_cache=False, flatten=False)

